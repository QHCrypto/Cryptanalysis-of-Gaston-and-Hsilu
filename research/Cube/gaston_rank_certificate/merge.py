"""Combine non-overlapping interval sums; only full coverage certifies rank."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from gaston_6r_ppm_transfer_audit import rotate_output
from structure import rank


def merge(paths):
    records = [json.loads(Path(p).read_text()) for p in paths]
    records.sort(key=lambda r: r['begin'])
    metadata = ['implementation', 'word_bits', 'rounds', 'period', 'cube_dimension', 'iv_hex', 'total_points']
    first = records[0]
    sums = [[0, 0] for _ in range(3)]
    covered = 0
    previous_end = 0
    gap = False
    for r in records:
        if any(r[k] != first[k] for k in metadata):
            raise ValueError('incompatible interval metadata')
        if r['begin'] < previous_end:
            raise ValueError('overlapping intervals')
        if r['begin'] != previous_end:
            gap = True
        if r['count'] <= 0 or r['begin'] + r['count'] > first['total_points']:
            raise ValueError('invalid interval')
        previous_end = r['begin'] + r['count']
        covered += r['count']
        for k in range(3):
            for w in range(2):
                sums[k][w] ^= int(r['sums_hex'][k][w], 16)
    complete = not gap and covered == first['total_points']
    result = {k: first[k] for k in metadata}
    result.update({'files': list(map(str, paths)), 'covered_points': covered,
                   'coverage_complete': complete,
                   'sums_hex': [[f'{x:016x}' for x in row] for row in sums]})
    if complete:
        v = [(sums[a+1][0] ^ sums[0][0]) | ((sums[a+1][1] ^ sums[0][1]) << 64)
             for a in range(2)]
        columns = [rotate_output(x, j, 64) for x in v for j in range(64)]
        residues = [sum((((x >> (64*w)) & ((1 << 64)-1)).bit_count() & 1) << w
                        for w in range(2)) for x in v]
        actual_rank = rank(columns)
        assert actual_rank >= 64 * rank(residues)
        assert (actual_rank == 128) == (rank(residues) == 2)
        result.update({'representative_columns_hex': [f'{x:032x}' for x in v],
                       'residue_columns': residues, 'residue_rank': rank(residues),
                       'matrix_rank': actual_rank,
                       'rank_lower_bound_from_residue_only': 64 * rank(residues)})
    else:
        result['rank_certificate'] = None
        result['reason'] = 'Partial sums do not determine any full-cube matrix entry.'
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('intervals', nargs='+')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = merge(args.intervals)
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
