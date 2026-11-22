"""Exact MP/MILP for two representative output rows of a parity cube.

Completed row coefficients generate two columns and the full block-circulant
matrix. This is not a claim that two solver calls finish path enumeration.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'reference'))
from parity_cube import gf2_rank, matrix_from_rows


def summarize(directory, rounds, width, key_split, individual=False):
    rows = [[None]*(2*width) for _ in range(2)]
    certified = []
    for b in range(2):
        if individual:
            for k in range(2*width):
                path = directory/f'key_{k}_row_{b}.json'
                if not path.exists():
                    continue
                result = json.loads(path.read_text())
                if (result['cube_kind'] != 'even_parity_periodic' or result['rounds'] != rounds
                        or result['word_bits'] != width or result['output_bit'] != b*width
                        or result['key_index'] != k or result['key_degree_split'] != bool(key_split)):
                    raise ValueError(f'Incompatible result: {path}')
                if result['complete']:
                    value = result['coefficient']
                    if value not in (0, 1):
                        raise ValueError('invalid certified coefficient')
                    rows[b][k] = value
            if all(x is not None for x in rows[b]):
                certified.append(b)
            continue
        path = directory/f'row_{b}.json'
        if not path.exists():
            continue
        result = json.loads(path.read_text())
        if (result['cube_kind'] != 'even_parity_periodic' or result['rounds'] != rounds
                or result['word_bits'] != width or result['output_bit'] != b*width
                or result['key_index'] != -1 or result['key_degree_split'] != bool(key_split)):
            raise ValueError(f'Incompatible result: {path}')
        if result['complete']:
            rows[b] = [result['coefficients'][str(k)] for k in range(2*width)]
            if any(x not in (0, 1) for x in rows[b]):
                raise ValueError('invalid certified coefficients')
            certified.append(b)
    matrix = matrix_from_rows(rows, width)
    complete = len(certified) == 2
    result = dict(cube_kind='even_parity_periodic', rounds=rounds, word_bits=width,
                  cube_period=1 << rounds, cube_dimension=(1 << rounds)-1,
                  parameterization='N1=0; N0 repeats a block with bit0=sum(u_j), bitj=u_j',
                  key_degree_split=bool(key_split), representative_output_rows=rows,
                  certified_representative_rows=certified,
                  representative_columns=[[matrix[i][a*width] for i in range(2*width)] for a in range(2)],
                  matrix=matrix, complete=complete,
                  unknown_entries=sum(x is None for row in matrix for x in row),
                  matrix_rank=gf2_rank(matrix) if complete else None,
                  certified_one_terms=[[i,k] for i,row in enumerate(matrix) for k,x in enumerate(row) if x == 1],
                  constant_term=None, constant_term_computed=False,
                  scope='Only exhaustive path parity or proven zero supplies coefficients. The constant term needs a full zero-key cube sum or a separate exact computation.')
    temp = directory/'matrix.json.tmp'
    temp.write_text(json.dumps(result, indent=2)+'\n')
    temp.replace(directory/'matrix.json')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rounds', type=int, default=6, choices=range(1,7))
    p.add_argument('--width', type=int, default=64)
    p.add_argument('--mode', choices=('build','screen','count'), default='count')
    p.add_argument('--output-words', default='0,1')
    p.add_argument('--keys', default='all', help='all: sequential k0,...; comma-separated indices; joint: legacy joint search')
    p.add_argument('--key-split', type=int, choices=(0,1), default=1)
    p.add_argument('--pool', type=int, default=256)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--directory', type=Path)
    p.add_argument('--summarize-only', action='store_true')
    args = p.parse_args()
    period = 1 << args.rounds
    if not period <= args.width <= 64 or args.width % period:
        p.error('2^rounds must divide width, with width <= 64')
    words = [int(x) for x in args.output_words.split(',')]
    if not words or len(set(words)) != len(words) or any(x not in (0,1) for x in words):
        p.error('output-words must be 0, 1, or 0,1')
    individual = args.keys != 'joint'
    try:
        keys = list(range(2*args.width)) if args.keys == 'all' else ([-1] if not individual else [int(x) for x in args.keys.split(',')])
    except ValueError:
        p.error('keys must be all, joint, or comma-separated indices')
    if individual and (not keys or len(set(keys)) != len(keys) or any(k < 0 or k >= 2*args.width for k in keys)):
        p.error('invalid key indices')
    suffix = '_individual' if individual else ''
    directory = args.directory or ROOT/f'results/parity/r{args.rounds}_w{args.width}_split{args.key_split}_{args.mode}{suffix}'
    directory.mkdir(parents=True, exist_ok=True)
    def refresh():
        result = summarize(directory, args.rounds, args.width, args.key_split, individual)
        print(json.dumps({k:result[k] for k in ('complete','certified_representative_rows',
                                              'unknown_entries','matrix_rank','constant_term_computed')}), flush=True)
    refresh()
    if args.summarize_only:
        return 0
    for k, b in ((k, b) for k in keys for b in words):
        stem = f'key_{k}_row_{b}' if individual else f'row_{b}'
        result = directory/f'{stem}.json'
        checkpoint = directory/f'{stem}.trails'
        if result.exists() and json.loads(result.read_text())['complete']:
            continue
        cmd = [str(ROOT/'gaston_parity_cube'),'--rounds',str(args.rounds),'--width',str(args.width),
               '--key',str(k) if individual else 'all','--output-bit',str(b*args.width),'--key-split',str(args.key_split),
               '--mode',args.mode,'--pool',str(args.pool),'--threads',str(args.threads),
               '--result',str(result)]
        if args.mode == 'count':
            cmd += ['--checkpoint',str(checkpoint)]
            if checkpoint.exists():
                cmd += ['--resume']
        child = subprocess.Popen(cmd)
        try:
            code = child.wait()
        except KeyboardInterrupt:
            code = child.wait()  # The foreground process group already received SIGINT.
            refresh()
            return 130
        refresh()
        if code:
            return code
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
