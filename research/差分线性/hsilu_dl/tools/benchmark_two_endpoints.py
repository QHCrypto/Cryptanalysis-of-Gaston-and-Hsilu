#!/usr/bin/env python3
"""Time fresh DATF engines and HNW processes for H4d and H5.

No DATF memoized state is shared between repetitions or endpoints. HNW uses
the existing rd=1, rl=0 local-column estimator, including multi-column product
readout. These timings do not concern exact full-permutation correlations.
"""
import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import time
from fractions import Fraction
from pathlib import Path

from hsilu import add, construct
from datf import support
from datf_k06 import value_record

ROOT = Path(__file__).resolve().parents[1]


def fraction(record):
    return Fraction(int(record['numerator']), int(record['denominator']))


def timing(values):
    return {'min': min(values), 'median': statistics.median(values), 'max': max(values)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repetitions', type=int, default=5)
    args = parser.parse_args()
    assert args.repetitions >= 1
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = json.loads((ROOT/'results/research_summary.json').read_text())
    jobs = [next(r for r in source[key] if r['id'] == identifier)
            for key, identifier in [('four_DATF', 'R4_pool_000001'),
                                    ('five_DATF', 'R5_pool_000000')]]
    report = {'method_scope': 'ordinary DATF and HNW local-column approximations; not exact permutation correlations',
              'boundary': 'same Hsilu core as archived results: chi first, final L omitted; first RC absorbed into uniform state',
              'datf_timing_scope': 'fresh symbolic DAG construction, target construction, and estimate; Python startup/import excluded',
              'hnw_timing_scope': 'fresh process wall time includes startup, tables, propagation and JSON output; internal propagation time also retained',
              'repetitions': args.repetitions, 'parallel_workers': 1,
              'load_average_at_start': os.getloadavg(), 'platform': platform.platform(),
              'results': []}
    cpuinfo = Path('/proc/cpuinfo').read_text()
    report['cpu_model'] = next(line.split(':', 1)[1].strip() for line in cpuinfo.splitlines()
                               if line.startswith('model name'))
    paths = [Path(__file__), ROOT/'tools/hsilu.py', ROOT/'build/hnw_endpoint',
             ROOT.parent/'gaston_dl/tools/datf.py']
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for job in jobs:
        identifier, rounds = job['id'], job['rounds']
        delta = tuple(int(w, 16) for w in job['input_difference']['words'])
        mask = tuple(int(w, 16) for w in job['core_output_mask']['words'])
        row = {'id': identifier, 'rounds': rounds, 'start_round_index': 12-rounds,
               'input_difference': job['input_difference'], 'output_mask': job['core_output_mask'],
               'archived_DATF_seconds': job['elapsed_seconds'], 'DATF_runs': [], 'HNW': {}}
        for repeat in range(args.repetitions):
            began = time.perf_counter()
            engine, _, difference = construct(delta, rounds, max_nodes=500000)
            target = add(*(difference[i][j] for i in range(5) for j in range(64) if mask[i] >> j & 1))
            constructed = time.perf_counter()
            trace = []
            rho = engine.estimate(target, label=identifier, trace=trace)
            finished = time.perf_counter()
            assert rho == fraction(job['correlation'])
            assert trace == job['trace']
            row['DATF_runs'].append({'construct_seconds': constructed-began,
                'estimate_seconds': finished-constructed, 'total_seconds': finished-began,
                'statistics': dict(engine.stats)})
        row['DATF_correlation'] = value_record(rho)
        row['DATF_trace'] = trace
        row['target_terms'], row['target_variables'] = len(target), support(target).bit_count()
        row['DATF_seconds'] = {key: timing([r[key] for r in row['DATF_runs']])
                               for key in ['construct_seconds', 'estimate_seconds', 'total_seconds']}
        for mode in ['full', 'zero']:
            runs = []
            for repeat in range(args.repetitions):
                artifact = out/f'{identifier}_{mode}_{repeat+1:02d}.json'
                command = [str(ROOT/'build/hnw_endpoint'), mode, str(rounds), str(12-rounds),
                           *job['input_difference']['words'], *job['core_output_mask']['words'], str(artifact)]
                began = time.perf_counter()
                subprocess.run(command, check=True, capture_output=True, text=True)
                elapsed = time.perf_counter()-began
                result = json.loads(artifact.read_text())
                assert result['status'] == 'completed'
                assert result['completed_branches'] == result['expected_first_DDT_branches']
                assert fraction(result['probability_mass']) == 1
                assert sum((fraction(b['contribution']) for b in result['branches']), Fraction(0)) == fraction(result['correlation'])
                if runs:
                    assert result['correlation'] == runs[0]['correlation']
                runs.append({'wall_seconds': elapsed, 'internal_seconds': result['elapsed_seconds'],
                             'correlation': result['correlation'], 'branches': result['completed_branches'],
                             'absorbed_histogram': result['absorbed_histogram'], 'artifact': str(artifact)})
            row['HNW'][mode] = {'runs': runs, 'correlation': runs[0]['correlation'],
                               'branches': runs[0]['branches'],
                               'wall_seconds': timing([r['wall_seconds'] for r in runs]),
                               'internal_seconds': timing([r['internal_seconds'] for r in runs])}
        report['results'].append(row)
        (out/'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
        print(json.dumps({'id': identifier, 'DATF_seconds': row['DATF_seconds'],
            'HNW': {mode: {key: data[key] for key in ['correlation', 'branches', 'wall_seconds', 'internal_seconds']}
                    for mode, data in row['HNW'].items()}}, ensure_ascii=False), flush=True)
    print('SAVED', out/'summary.json', flush=True)


if __name__ == '__main__':
    main()
