#!/usr/bin/env python3
"""Regression checks of certified round-two initialization versus the legacy updater."""
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    data = json.loads((ROOT / 'results/fixed_tail_extension/results.json').read_text())
    row = next(r for r in data['rows'] if r['id'] == 'D1_M1_D3_chi1_to_3')
    c = row['candidates'][0]
    prefix = c['join_mask']['words']
    delta = c['input_difference']['words']
    checks = 0
    indices = [0, 1, 3, 15, 255, 1023, 12345, 32768, 49151, 65535]
    with tempfile.TemporaryDirectory(prefix='gaston-hnw-') as tmp:
        for mode in ('full', 'zero'):
            for mask in (prefix, ['1', '0', '0', '0', '0']):
                for index in indices:
                    pair = []
                    for fast in (False, True):
                        path = Path(tmp) / 'job.json'
                        command = [str(ROOT / 'build/hnw_endpoint'), mode, '2', '7', *delta, *mask,
                                   str(path), str(index), str(index + 1)]
                        if fast:
                            command += prefix
                        subprocess.run(command, check=True, capture_output=True)
                        pair.append(json.loads(path.read_text()))
                    a, b = pair
                    assert a['correlation'] == b['correlation']
                    assert a['probability_mass'] == b['probability_mass']
                    assert all(b['branches'][0][k] == v for k, v in a['branches'][0].items())
                    if mask == prefix:
                        assert b['correlation'] == b['prefix_correlation']
                    checks += 1
    result = {'status': 'PASS', 'legacy_vs_certified_initialization_cases': checks,
              'sampled_global_branch_indices': indices, 'rounds': 2,
              'readouts': ['dense prefix mask', 'single output bit'], 'modes': ['full', 'zero']}
    print(json.dumps(result))


if __name__ == '__main__':
    main()
