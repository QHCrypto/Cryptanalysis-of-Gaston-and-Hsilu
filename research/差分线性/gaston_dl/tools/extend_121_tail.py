#!/usr/bin/env python3
"""Append one linear round to the saved 1+2+1 Gaston candidate.

Exhaustive single-characteristic optimization for the fixed tail entrance;
complete endpoint DATF/HNW checks cover one deterministic representative per
selected extension, not all tied final masks or a complete linear hull.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import math
import subprocess
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

from fixed_tail_search import RC, S, checked_round, column, dot, embed, lat, linear, transpose
from linear_seed_scan import active, encode
from datf_paper_reassess import evaluate_endpoint
from hnw26_compare import run_job, save, fraction

ROOT = Path(__file__).resolve().parents[1]
MASK64 = (1 << 64) - 1


def pack(state):
    return sum(w << (64*i) for i, w in enumerate(state))


def unpack(value):
    return tuple((value >> (64*i)) & MASK64 for i in range(5))


def decode(words):
    return tuple(int(w, 16) for w in words)


def words(state):
    return [f'{w:016x}' for w in state]


def value(q):
    q = Fraction(q)
    return {'numerator': str(q.numerator), 'denominator': str(q.denominator),
            'log2_abs': math.log2(abs(q.numerator))-math.log2(q.denominator) if q else None}


class InverseTranspose:
    def __init__(self):
        self.basis = {}
        for i in range(320):
            vector, witness = pack(transpose(unpack(1 << i))), 1 << i
            while vector:
                pivot = vector.bit_length()-1
                if pivot not in self.basis:
                    self.basis[pivot] = vector, witness
                    break
                row, other = self.basis[pivot]
                vector ^= row
                witness ^= other
        assert len(self.basis) == 320, 'L must be invertible'

    def __call__(self, target):
        vector, witness = pack(target), 0
        while vector:
            row, other = self.basis[vector.bit_length()-1]
            vector ^= row
            witness ^= other
        result = unpack(witness)
        assert transpose(result) == target
        return result


def assemble(positions, patterns):
    return tuple(sum(((v >> i) & 1) << j for j, v in zip(positions, patterns)) for i in range(5))


def solve(first_input, original_output):
    inverse = InverseTranspose()
    table = [[lat(a, b) for b in range(32)] for a in range(32)]
    assert sorted(S) == list(range(32))
    assert all(table[0][b] == (32 if b == 0 else 0) for b in range(32))
    assert all(table[a][0] == (32 if a == 0 else 0) for a in range(32))
    positions = [j for j in range(64) if column(first_input, j)]
    options = [[b for b in range(1, 32) if table[column(first_input, j)][b]] for j in positions]
    records = []
    basis_checks = 0
    for patterns in itertools.product(*options):
        b1 = assemble(positions, patterns)
        a2 = inverse(b1)
        # Independent check against forward L, including both rho layers.
        for i in range(320):
            assert dot(linear(unpack(1 << i)), a2) == ((pack(b1) >> i) & 1)
            basis_checks += 1
        chosen, local_choices, maximum_numerators = [], [], []
        for j in range(64):
            a = column(a2, j)
            maximum = max(abs(v) for v in table[a])
            candidates = [b for b in range(32) if abs(table[a][b]) == maximum]
            # Reproducible representative among equally good final masks.
            chosen.append(min(candidates, key=lambda b: (b.bit_count(), b)))
            maximum_numerators.append(maximum)
            local_choices.append({'column': j, 'input_mask': a, 'optimal_output_masks': candidates,
                                  'maximum_abs_LAT_numerator': maximum})
        b2 = assemble(range(64), chosen)
        r1, r2 = checked_round(first_input, b1), checked_round(a2, b2)
        total_w2 = r1['w_squared'] + r2['w_squared']
        bare_sign = r1['sign'] * r2['sign']
        phase = ((first_input[0] & RC[10]).bit_count() + (a2[0] & RC[11]).bit_count()) & 1
        correlation = Fraction(bare_sign * (-1 if phase else 1), 2**(total_w2//2))
        records.append({'preserves_original_last_mask': b1 == original_output,
                        'first_input': encode(first_input), 'first_output': encode(b1),
                        'second_input': encode(a2), 'final_output': encode(b2),
                        'first_round_squared_weight': r1['w_squared'],
                        'added_round_squared_weight': r2['w_squared'],
                        'tail_squared_weight': total_w2,
                        'added_round_active_columns': active(a2),
                        'added_round_active_bits': sum(w.bit_count() for w in a2),
                        'bare_LAT_sign': bare_sign, 'constant_phase': phase,
                        'tail_correlation': value(correlation),
                        'tail_squared_correlation': value(correlation**2),
                        'optimal_final_mask_count': math.prod(len(c['optimal_output_masks']) for c in local_choices),
                        'optimal_final_masks_factorization': local_choices,
                        'selected_output_rule': 'minimum bit weight per column, then smallest five-bit integer'})
    records.sort(key=lambda r: (r['tail_squared_weight'], r['first_output']['words']))
    for i, r in enumerate(records, 1):
        r['rank'] = i
    audit = {'status': 'PASS', 'inverse_matrix_rank': len(inverse.basis),
             'forward_L_basis_checks': basis_checks, 'first_active_columns': positions,
             'first_chi_output_choices_per_active_column': list(map(len, options)),
             'first_chi_output_masks_exhausted': len(records),
             'tail_squared_weight_histogram': dict(sorted(Counter(r['tail_squared_weight'] for r in records).items())),
             'added_active_columns_histogram': dict(sorted(Counter(r['added_round_active_columns'] for r in records).items())),
             'coverage': 'Every nonzero first-chi output for the fixed first-chi input; every last-chi mask optimized by independent local LAT maxima.',
             'optimality_scope': 'Single two-round linear characteristic with this fixed entrance; not the full linear hull or full differential-linear correlation.'}
    return records, audit


def evaluate_datf(job, folder):
    cached = folder / (job['id'] + '.summary.json')
    if cached.exists():
        result = json.loads(cached.read_text())
        with gzip.open(result['artifact'], 'rt') as f:
            old = json.load(f)['result']
        assert all(old[k] == v for k, v in job.items())
        return result
    return evaluate_endpoint(job, str(folder), 1_000_000, 32)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=ROOT / 'results/extend_121_tail')
    args = p.parse_args()
    out = args.output.resolve()
    (out / 'datf_jobs').mkdir(parents=True, exist_ok=True)
    source_path = ROOT / 'results/datf_pools/merged/merged_results.json'
    source = json.loads(source_path.read_text())
    row = next(r for r in source['rows'] if r['id'] == 'D1_M2_L1_pool')
    candidate = row['selected'] or row['best']
    with gzip.open(ROOT / 'results/datf_pools/pools/linear_one.jsonl.gz', 'rt') as f:
        tail = next(r for r in map(json.loads, f) if r['id'] == candidate['linear_id'])
    first_input, original_output = decode(tail['input']), decode(tail['output'])
    assert original_output == decode(candidate['output_mask'])
    assert transpose(first_input) == decode(candidate['middle_output_mask'])
    config = {'source_candidate': candidate, 'original_linear_tail': tail,
              'first_L_omitted': True, 'target_split': '1+2+2', 'rounds': 5, 'start_round_index': 7,
              'constant_convention': 'Use the last five reference rounds. Original four-round trail geometry is shifted one reference round earlier; no thirteenth round constant is invented.',
              'representative_selection': 'For each first output, minimize added absolute-correlation weight; break local output ties by bit weight then integer value.',
              'source_hashes': {str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in
                               [source_path, ROOT / 'results/datf_pools/pools/linear_one.jsonl.gz', Path(__file__),
                                ROOT / 'tools/fixed_tail_search.py', ROOT / 'tools/linear_seed_scan.py',
                                ROOT / 'build/hnw_endpoint', ROOT / 'build/datf_middle_cli']}}
    cp = out / 'configuration.json'
    if cp.exists():
        assert json.loads(cp.read_text()) == config, 'Configuration changed; use a new output directory.'
    else:
        save(cp, config)
    began = time.perf_counter()
    records, audit = solve(first_input, original_output)
    strict = next(r for r in records if r['preserves_original_last_mask'])
    best = records[0]
    assert len(records) == 100
    assert strict['tail_squared_weight'] == 132 and best['tail_squared_weight'] == 124
    save(out / 'all_extensions.json', {'records': records, 'audit': audit})
    print('SCAN', len(records), 'strict', strict['tail_squared_weight'], 'best', best['tail_squared_weight'], flush=True)
    # The shifted middle has internal RC[9]=0x69 instead of old RC[10]=0x5a.
    middle_input = decode(candidate['middle_input_difference'])
    middle_mask = decode(candidate['middle_output_mask'])
    rows = [' '.join(['2', *words(middle_input), *words(middle_mask), f'{rc:x}']) for rc in (RC[10], RC[9])]
    reply = subprocess.check_output([str(ROOT / 'build/datf_middle_cli')], input='\n'.join(rows)+'\n', text=True)
    middle_controls = []
    for text, rc in zip(reply.splitlines(), (RC[10], RC[9])):
        sign, weight, phase = text.split()
        middle_controls.append({'second_constant': f'{rc:02x}', 'correlation': value(Fraction(int(sign), 2**int(weight))),
                                'constant_phase_mask': phase})
    middle_job = {'id': 'middle_shifted', 'kind': 'middle_control', 'rounds': 2, 'start_round_index': 8,
                  'input_difference': {'words': words(middle_input)}, 'output_mask': {'words': words(middle_mask)}}
    middle_datf = evaluate_datf(middle_job, out / 'datf_jobs')
    assert middle_datf['status'] == 'completed'
    assert fraction(middle_datf['correlation']) == fraction(middle_controls[1]['correlation']) == Fraction(1, 4)
    assert fraction(middle_controls[0]['correlation']) == fraction(candidate['middle_DATF_correlation'])
    results = []
    for name, extension in [('strict_append', strict), ('reoptimized_tail', best)]:
        mask = words(decode(extension['final_output']['words']))
        job = {'id': name, 'kind': 'full_endpoint', 'rounds': 5, 'start_round_index': 7,
               'input_difference': {'words': candidate['input_difference']}, 'output_mask': {'words': mask}}
        datf = evaluate_datf(job, out / 'datf_jobs')
        hj = {'id': name, 'candidate': {'R': 5, 'start_round_index': 7,
                                      'input_difference': candidate['input_difference'], 'output_mask': mask}}
        hnw = {mode: run_job(hj, mode, out) for mode in ('full', 'zero')}
        splice = (fraction(candidate['differential_probability']) * fraction(middle_controls[1]['correlation']) *
                  fraction(extension['tail_squared_correlation']))
        results.append({'id': name, 'extension': extension, 'splice_score': value(splice), 'DATF': datf, 'HNW': hnw})
        print('EVALUATED', name, 'score', str(splice), 'DATF', datf['correlation'], flush=True)
    # Recheck the old endpoint after shifting its reference interval one round
    # earlier, so the five-round construction has a clear constant boundary.
    rebased_job = {'id': 'old_endpoint_shifted', 'kind': 'full_endpoint_control', 'rounds': 4, 'start_round_index': 7,
                   'input_difference': {'words': candidate['input_difference']},
                   'output_mask': {'words': candidate['output_mask']}}
    rebased = evaluate_datf(rebased_job, out / 'datf_jobs')
    document = {'status': 'completed', 'source_row_id': row['id'], 'source_candidate': candidate,
                'source_full_DATF': row['selected_correlation'], 'rounds': 5, 'start_round_index': 7,
                'first_L_omitted': True, 'middle_controls': middle_controls,
                'shifted_middle_DATF': middle_datf, 'shifted_old_endpoint_DATF': rebased,
                'results': results, 'search_audit': audit, 'elapsed_seconds': time.perf_counter()-began,
                'full_endpoint_coverage': 'One deterministic optimal last-mask representative for strict append and one for the globally best tail with fixed entrance. Other equally weighted output masks were not endpoint-evaluated.',
                'interpretation': 'Optimality concerns single-characteristic LAT products. DATF and HNW are full-endpoint model estimates, not exact full-permutation correlations.'}
    save(out / 'results.json', document)
    save(out / 'validation.json', {'status': 'PASS', 'linear_scan': audit,
          'two_independent_middle_implementations_agree': True,
          'full_DATF_factor_audits': {r['id']: r['DATF']['factor_audit_status'] for r in results},
          'full_HNW_first_DDT_branches_per_endpoint_and_mode': 4,
          'full_HNW_branch_probability_and_contribution_sums_checked': True,
          'dag_direct_permutation_samples_per_DATF_job': 32})
    print('RESULT', out / 'results.json', flush=True)


if __name__ == '__main__':
    main()
