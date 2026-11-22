"""Sequential exact computation of the 128 x 128 linear-superpoly matrix.

No solve timeout. Checkpoints permit restart; unknown entries are JSON null.
Gurobi logs inherit the terminal. Ctrl+C asks the current solver to checkpoint.
"""
import argparse
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outputs', default='0-127')
    parser.add_argument('--pool', type=int, default=256)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--key-split', type=int, choices=(0,1), default=1)
    args = parser.parse_args()
    outputs = []
    for part in args.outputs.split(','):
        if '-' in part:
            a,b = map(int, part.split('-')); outputs += list(range(a,b+1))
        else:
            outputs.append(int(part))
    if not outputs or any(i<0 or i>127 for i in outputs) or len(outputs)!=len(set(outputs)):
        parser.error('outputs must be distinct indices in 0..127')
    config = json.loads((ROOT/'cube_config.json').read_text())
    # Split and legacy circuits enumerate different paths; keep checkpoints apart.
    directory = ROOT / ('results/matrix_split1' if args.key_split else 'results/matrix')
    directory.mkdir(exist_ok=True)
    nonce = ','.join(map(str, config['nonce_indices']))

    def summarize():
        rows = [[None]*128 for _ in range(128)]
        complete = []
        for p in directory.glob('output_*.json'):
            r = json.loads(p.read_text())
            if (r['rounds']!=6 or r['word_bits']!=64 or r['nonce_indices']!=config['nonce_indices']
                    or r['key_index']!=-1 or r.get('key_degree_split',False)!=bool(args.key_split)):
                raise ValueError('matrix output does not match cube_config.json')
            if r['complete']:
                j = r['output_bit']
                rows[j] = [r['coefficients'][str(k)] for k in range(128)]
                complete.append(j)
        result = dict(nonce_indices=config['nonce_indices'], matrix=rows,
                      certified_output_rows=sorted(complete),
                      certified_one_entries=[[j,k] for j,row in enumerate(rows) for k,v in enumerate(row) if v==1],
                      complete=len(complete)==128, key_degree_split=bool(args.key_split),
                      unknown_entries=sum(v is None for row in rows for v in row))
        temp = directory/'matrix.json.tmp'
        temp.write_text(json.dumps(result, indent=2)+'\n')
        temp.replace(directory/'matrix.json')

    summarize()
    for output in outputs:
        result = directory/f'output_{output:03d}.json'
        checkpoint = directory/f'output_{output:03d}.trails'
        if result.exists() and json.loads(result.read_text())['complete']:
            continue
        cmd = [str(ROOT/'gaston_cube'), '--rounds','6','--key','all',
               '--nonce',nonce,'--output-bit',str(output),'--mode','count',
               '--key-split',str(args.key_split),
               '--threads',str(args.threads),'--pool',str(args.pool),
               '--result',str(result),'--checkpoint',str(checkpoint)]
        if checkpoint.exists():
            cmd.append('--resume')
        child = subprocess.Popen(cmd)
        try:
            status = child.wait()
        except KeyboardInterrupt:
            # Same terminal process group already receives SIGINT.
            status = child.wait()
            summarize()
            return 130
        summarize()
        if status:
            return status
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
