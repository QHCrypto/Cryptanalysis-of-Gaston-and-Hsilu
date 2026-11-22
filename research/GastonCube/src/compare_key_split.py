"""Compare rooted MILP sizes for the fixed research cube, without claiming speed."""
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    config = json.loads((ROOT/'cube_config.json').read_text())
    directory = ROOT/'results/key_split_comparison'
    directory.mkdir(exist_ok=True)
    nonce = ','.join(map(str, config['nonce_indices']))
    rows = []
    for key in ('0', 'all'):
        for split in (0, 1):
            result = directory/f'k{key}_split{split}.json'
            subprocess.run([str(ROOT/'gaston_cube'), '--rounds','6','--key',key,
                            '--nonce',nonce,'--output-bit','0','--key-split',str(split),
                            '--mode','build','--result',str(result)], check=True)
            r = json.loads(result.read_text())
            rows.append({k:r[k] for k in ('key_index','key_degree_split','circuit_nodes',
                                         'binary_variables','constraints')})
    (directory/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__ == '__main__':
    main()
