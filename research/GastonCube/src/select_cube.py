"""Choose among explicit coordinate cubes by a declared MILP-size proxy."""
import argparse
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = {
    'n0_consecutive': list(range(63)),
    'n1_consecutive': list(range(64, 127)),
    'two_contiguous_blocks': list(range(32)) + list(range(64, 95)),
    'even_coordinates': list(range(0, 126, 2)),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--key-split', type=int, choices=(0,1), default=1)
    args = parser.parse_args()
    rows = []
    for name, nonce in CANDIDATES.items():
        result = ROOT / 'results' / (name + '_size.json')
        cmd = [str(ROOT / 'gaston_cube'), '--rounds', '6', '--key', '0',
               '--nonce', ','.join(map(str, nonce)), '--output-bit', '0',
               '--mode', 'build', '--key-split', str(args.key_split), '--result', str(result)]
        subprocess.run(cmd, check=True)
        data = json.loads(result.read_text())
        rows.append(dict(name=name, nonce_indices=nonce,
                         binary_variables=data['binary_variables'],
                         constraints=data['constraints'], circuit_nodes=data['circuit_nodes']))
    rows.sort(key=lambda r: (r['binary_variables'], r['constraints']))
    config = dict(research_object='Gaston-AEAD128 reduced initialization',
                  rounds=6, word_bits=64, iv_hex='00001000808c0001',
                  input_words=['IV','K0','K1','N0','N1'],
                  observable='first ciphertext block; empty AD and zero first plaintext block',
                  output_bits=list(range(128)), key_bits=list(range(128)),
                  nonce_indices=rows[0]['nonce_indices'], noncube_assignment=0,
                  bit_order='LSB first within each little-endian 64-bit word',
                  criterion='Minimum binary variables among four fixed cubes for key k0 and output S0[0]. This is a size proxy, not a proven minimum runtime.',
                  model_key_degree_split=bool(args.key_split),
                  candidates=rows, full_coefficient_matrix_complete=False)
    (ROOT/'cube_config.json').write_text(json.dumps(config, indent=2)+'\n')
    print(json.dumps({k:v for k,v in config.items() if k not in ('output_bits','key_bits','candidates')}, indent=2))


if __name__ == '__main__':
    main()
