#!/usr/bin/env python3
"""Reproduce the complete support proof, or check the saved computation.

All enumeration is sequential and single-threaded, without time limits.
--recompute replaces the named 20261024 generated reports, preserving logs in
a fresh directory; it does not start, stop, or change a MILP solver.
"""
import argparse
from datetime import datetime
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recompute', action='store_true',
                        help='rebuild and rerun both complete enumerations; otherwise check saved outputs')
    args = parser.parse_args()
    if args.recompute:
        logs = ROOT / 'results' / ('support_reproduction_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
        logs.mkdir()
        commands = []

        def run(command, name):
            print('Running:', ' '.join(map(str, command)), flush=True)
            start = time.monotonic()
            result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
            elapsed = time.monotonic() - start
            (logs / (name + '.stdout')).write_text(result.stdout)
            (logs / (name + '.stderr')).write_text(result.stderr)
            commands.append({'command': list(map(str, command)), 'seconds': elapsed,
                             'exit_code': result.returncode})
            (logs / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
            result.check_returncode()
            return result.stdout, elapsed

        run(['make', 'support-proof-binaries'], 'build')
        run([sys.executable, 'src/build_support_reference.py'], 'reference')
        for k, name in ((6, 'six'), (7, 'seven'), (8, 'eight')):
            raw, _ = run(['src/support_rank_probe', '--full', str(k)], f'primary_{k}')
            data = json.loads(raw)
            (ROOT / f'results/support_rank_{name}_20261024.json').write_text(json.dumps(data, indent=2) + '\n')
            raw, elapsed = run(['src/verify_support_ranks', 'results/support_reference_matrices_20261024.txt', str(k)],
                               f'independent_{k}')
            data = json.loads(raw)
            data.update(seconds=elapsed, middle_active_columns=k, complete_anchored_enumeration=True)
            (ROOT / f'results/support_rank_{name}_independent_20261024.json').write_text(json.dumps(data, indent=2) + '\n')
        run([sys.executable, 'src/analyze_support_spaces.py'], 'nullspaces')
        raw, _ = run([sys.executable, 'src/verify_support_proof.py'], 'proof')
        print(raw)
        print('Reproduction logs:', logs)
    else:
        subprocess.run([sys.executable, 'src/verify_support_proof.py'], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
