"""Direct three-cube interpolation, including constants; for small validation runs."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'reference'))
from gaston import IV
from parity_cube import cube_sum, gf2_rank, matrix_from_columns


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rounds', type=int, default=3, choices=range(1,7))
    p.add_argument('--width', type=int, default=64)
    p.add_argument('--result', type=Path, required=True)
    args = p.parse_args()
    period = 1 << args.rounds
    if not period <= args.width <= 64 or args.width % period:
        p.error('2^rounds must divide width, with width <= 64')
    points = 1 << (period-1)
    print(f'Three complete cube sums: {points} points each, {3*points} permutation evaluations.', flush=True)
    start = time.monotonic()
    sums = [cube_sum(key, args.rounds, args.width) for key in (0, 1, 1 << args.width)]
    columns = [sums[1] ^ sums[0], sums[2] ^ sums[0]]
    matrix = matrix_from_columns(columns, args.width)
    h = lambda value: f'{value:0{args.width//2}x}'
    result = dict(cube_kind='even_parity_periodic',rounds=args.rounds,word_bits=args.width,
                  cube_period=period,cube_dimension=period-1,iv_hex=f'{IV:016x}',
                  evaluations=3*points,complete=True,sums_hex=[h(x) for x in sums],
                  constant_term_hex=h(sums[0]),representative_columns_hex=[h(x) for x in columns],
                  matrix=matrix,matrix_rank=gf2_rank(matrix),seconds=time.monotonic()-start,
                  scope='Complete Boolean cube sums with IV, round constants and all chi terms. Values apply only to the recorded round count and width.')
    args.result.parent.mkdir(parents=True,exist_ok=True)
    temp = args.result.with_name(args.result.name+'.tmp')
    temp.write_text(json.dumps(result,indent=2)+'\n')
    temp.replace(args.result)
    print(json.dumps({k:v for k,v in result.items() if k!='matrix'},indent=2))


if __name__ == '__main__':
    main()
