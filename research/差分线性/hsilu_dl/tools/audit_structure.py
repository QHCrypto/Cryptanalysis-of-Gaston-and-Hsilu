#!/usr/bin/env python3
"""Finite structural certificates, independent of DATF bias estimates."""
import argparse
import itertools
import json
import random
from pathlib import Path
from hsilu import *
from sparse_prefixes import rank, column, rank_certificate, DDT, middle

ROOT = Path(__file__).resolve().parents[1]


def decode(record):
    return tuple(int(x, 16) for x in record['words'])


def cone_linear(a):
    # OR over-approximates dependency through XOR. Cancellation can only
    # remove dependencies, so a missing bit proves absence of dependency.
    x = list(a)
    for u, v, n in zip(N1, N2, SH):
        x.append(x[u] | rol(x[v], n))
    return tuple(rol(x[13+i], ROT[i]) for i in range(5))


def cone_chi(a):
    return tuple(a[i] | a[(i+1) % 5] | a[(i+2) % 5] for i in range(5))


def dependency(delta, rounds):
    stages = []
    for r in range(rounds):
        if r:
            delta = cone_linear(delta)
        delta = cone_chi(delta)
        stages.append(encode(delta))
    return stages


def canonical(alpha, beta):
    return min((tuple(rol(v, s) for v in alpha), tuple(rol(v, s) for v in beta))
               for s in range(64))


def bridges_audit():
    invt = inverse_function(transpose)
    found, tested, hits = set(), 0, {}
    # If two positive activities sum to <=5, at least one is <=2.
    # Fix an active column at zero by the exact rotation symmetry of L.
    for side, f in [('beta', transpose), ('alpha', invt)]:
        table = {(j, p): f(embed(j, p)) for j in range(64) for p in range(1, 32)}
        hits[side] = 0
        for p in range(1, 32):
            fixed, other = embed(0, p), table[0, p]
            tested += 1
            if active(other) <= 4:
                alpha, beta = (other, fixed) if side == 'beta' else (fixed, other)
                found.add(canonical(alpha, beta)); hits[side] += 1
        for j, p, q in itertools.product(range(1, 64), range(1, 32), range(1, 32)):
            fixed = tuple(x ^ y for x, y in zip(embed(0, p), embed(j, q)))
            other = tuple(x ^ y for x, y in zip(table[0, p], table[j, q]))
            tested += 1
            if active(other) <= 3:
                alpha, beta = (other, fixed) if side == 'beta' else (fixed, other)
                found.add(canonical(alpha, beta)); hits[side] += 1
    old = json.loads((ROOT/'results/baseline/linear_bridges.json').read_text())
    assert found == {canonical(decode(r['alpha']), decode(r['beta'])) for r in old['bridges']}
    records = []
    # Descending canonical order preserves the identifiers in the saved pool.
    for alpha, beta in sorted(found, reverse=True):
        assert transpose(beta) == alpha
        # Independent inner products using the FORWARD linear map.
        for i in range(320):
            b = unpack(1 << i)
            assert (pack(linear(b)) & pack(beta)).bit_count() % 2 == (pack(alpha) >> i) & 1
        records.append({'alpha': encode(alpha), 'beta': encode(beta)})
    assert tested == 121148 and len(records) == 3
    return {'tested': tested, 'raw_hits': hits, 'rotation_classes': len(found),
            'scope': 'All nonzero alpha=L^T(beta) with total active columns <=5, modulo common rotation',
            'completeness_argument': 'At least one side has <=2 columns; enumerate both directions. Invertibility excludes zero sides.',
            'minimum_total_activity': min(active(a)+active(b) for a, b in found),
            'bridges': records, 'forward_basis_checks': len(records)*320,
            'constant_limitation': 'Rotation is exact for L, not asserted for full fixed-RC correlations.'}


def sparse_rank_audit():
    source = json.loads((ROOT/'results/sparse_diff/candidates.json').read_text())
    forward_columns = [linear(unpack(1 << i)) for i in range(320)]
    records = []
    pool = json.loads((ROOT/'results/baseline/linear_pool.json').read_text())
    checked_sources, tested_middle, matches = 0, 0, 0
    for p in source['prefixes']:
        x0, x1 = [decode(p['states'][k]) for k in ('x0', 'x1')]
        states = tuple(decode(p['states'][k]) for k in ('x0','y0','x1','y1'))
        assert linear(states[1]) == states[2]
        for r in (0,1):
            counts = [DDT[column(states[2*r],j)][column(states[2*r+1],j)] for j in range(64)]
            assert all(counts)
            assert sum(5-(n.bit_length()-1) for n in counts) == p['weights'][r]
        for path in p['sources']:
            lines = Path(path).read_text().splitlines(); original = {}
            for k,line in enumerate(lines):
                if line in ('x0','y0','x1','y1'):
                    original[line] = tuple(int(w.replace('.','0'),2) for w in lines[k+1:k+6])
            canon = min(tuple(tuple(rol(w,s) for w in original[k]) for k in ('x0','y0','x1','y1')) for s in range(64))
            assert canon == states
            checked_sources += 1
        before = linear(states[3])
        for s in range(64):
            delta = tuple(rol(w,s) for w in before)
            for head in pool['heads']:
                matches += middle(delta, decode(head['join_mask'])) != 0
                tested_middle += 1
        # Build affine derivative coefficients by finite differences of the
        # actual bitwise chi, not by the transposed symbolic formula.
        base = tuple(a ^ b for a, b in zip(chi((0,)*5), chi(x1)))
        images = []
        for y in forward_columns:
            derivative = tuple(a ^ b ^ c for a, b, c in zip(
                chi(y), chi(tuple(a ^ b for a, b in zip(y, x1))), base))
            images.append(pack(derivative))
        equations = [sum(((v >> j) & 1) << i for i, v in enumerate(images)) for j in range(320)]
        rows = [r for r in equations if r]
        activity = x0[0] | x0[1] | x0[2] | x0[3] | x0[4]
        fixed = pack((activity,)*5)
        full, outside = rank(rows), rank([r & ~fixed for r in rows])
        expected = rank_certificate(x0, x1)
        assert full == expected['full_constraint_rank'] == p['weights'][1]
        assert outside == expected['outside_first_active_columns_rank'] == full
        records.append({'id': p['id'], 'rank': full, 'outside_rank': outside,
                        'probability': '1/' + str(1 << sum(p['weights']))})
    assert len(records) == 10
    assert checked_sources == 18 and tested_middle == 245760 and matches == 0
    return {'method': 'Forward L basis and finite differences of chi; no L^T used to construct matrices',
            'prefix_certificates': records, 'all_exact_product_probabilities': True,
            'current_source_files_checked': checked_sources,
            'middle_matches_recomputed': matches, 'middle_connections_recomputed':tested_middle}


def structural_audit():
    invt = inverse_function(transpose)
    families = []
    for rounds, lane, col, missing in [(2, 0, 0, [(0, 1)]), (3, 3, 0, [(0, 29), (1, 29)])]:
        delta = embed(col, 1 << lane)
        stages = dependency(delta, rounds)
        final = decode(stages[-1])
        assert all(not (final[i] >> j & 1) for i, j in missing)
        e, _, d = construct(delta, rounds)
        assert all(d[i][j] == ZERO for i, j in missing)
        masks = []
        for subset in range(1, 1 << len(missing)):
            mask = [0]*5
            for k, (i, j) in enumerate(missing):
                if subset >> k & 1:
                    mask[i] ^= 1 << j
            masks.append({'core': encode(mask), 'full': encode(invt(mask))})
        # The conservative dependency proof holds for every input position;
        # constants never add a dependency on the varied input bit.
        for shift in range(64):
            rotated = tuple(rol(v, shift) for v in delta)
            cone = decode(dependency(rotated, rounds)[-1])
            assert all(not (cone[i] >> ((j+shift) % 64) & 1) for i, j in missing)
        families.append({'rounds': rounds, 'input': encode(delta), 'stages': stages,
                         'unreachable_target_bits': missing, 'masks': masks,
                         'rho': '1', 'all_64_rotations_certified': True,
                         'scope': 'Pointwise zero derivative for every base state; arbitrary RC values'})
    assert [r['bit_weight'] for r in families[1]['stages']] == [3, 66, 318]
    return families


def boundary_checks(families):
    selections = json.loads((ROOT/'results/baseline/pool_selection.json').read_text())
    invt = inverse_function(transpose)
    jobs = [(x['rounds'], decode(x['input_difference']), decode(x['core_output_mask'])) for x in selections]
    jobs += [(x['rounds'], decode(x['input']), decode(m['core'])) for x in families for m in x['masks']]
    rng = random.Random(2026101820)
    count = 0
    for rounds, delta, mask in jobs:
        fullmask = invt(mask)
        assert transpose(fullmask) == mask
        for _ in range(16):
            s = tuple(rng.getrandbits(64) for _ in range(5))
            t = tuple(a ^ b for a, b in zip(s, delta))
            def full(a):
                for rc in RC[12-rounds:]:
                    a = linear(chi((a[0] ^ rc, *a[1:])))
                return a
            out = tuple(a ^ b for a, b in zip(full(s), full(t)))
            shifted = (s[0] ^ RC[12-rounds], *s[1:])
            dcore = tuple(a ^ b for a, b in zip(core(shifted, rounds),
                         core(tuple(a ^ b for a, b in zip(shifted, delta)), rounds)))
            assert (pack(out) & pack(fullmask)).bit_count() % 2 == (pack(dcore) & pack(mask)).bit_count() % 2
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', type=Path)
    parser.add_argument('--export-bridges', type=Path,
                        help='Independently re-enumerate and export seeds for a fresh pool directory')
    args = parser.parse_args()
    if args.export_bridges:
        data = bridges_audit()
        args.export_bridges.parent.mkdir(parents=True, exist_ok=True)
        with args.export_bridges.open('x') as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('EXPORTED', args.export_bridges)
    if not args.output:
        if not args.export_bridges:
            parser.error('Provide --output or --export-bridges')
        return
    families = structural_audit()
    data = {'status': 'PASS', 'structural_families': families,
            'linear_bridges': bridges_audit(), 'sparse_prefix_ranks': sparse_rank_audit(),
            'full_core_endpoint_parity_checks': boundary_checks(families)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
    print('PASS structural 2/3 rounds, 121148 linear cases, 10 forward rank certificates, full/core masks')


if __name__ == '__main__':
    main()
