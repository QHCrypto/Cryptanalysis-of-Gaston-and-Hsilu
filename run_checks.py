#!/usr/bin/env python3
"""Run documented checks in a temporary copy; preserve archived research inputs."""
import argparse
import datetime
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time

PACKAGE = Path(__file__).resolve().parent
PY = sys.executable

def jobs(profile):
    # Each entry is (label, [(working directory below research, argv), ...]).
    paper = 'Gaston_Hsilu_CTP/paper_rewriting_output'
    h3 = '代数攻击/Hsilu3r_new/validation_current_prefix'
    if profile == 'quick':
        return [
            ('parameters', [('nist_compatibility_20261109',['make','native.so']),
                            (paper,[PY,'verify_section2.py'])]),
            ('differential', [(paper+'/ctp_reviewer_20261117',[PY,'independent_checks.py']),
                              ('Hsilu_trail/search_below82_20261117',[PY,'verify.py']),
                              ('gaston_active_pattern/two_round_outside_proof',[PY,'verify.py'])]),
            ('gaston_support', [('Gaston_差分路径重量下界',['make','support-proof-binaries']),
                                ('Gaston_差分路径重量下界',[PY,'src/verify_support_proof.py'])]),
            ('differential_linear', [('差分线性/gaston_dl',['make','test','build/hnw_endpoint']),
                                     ('差分线性/gaston_dl',['./build/hnw_endpoint','--selftest']),
                                     ('差分线性/hsilu_dl',['make','test','build/hnw_endpoint']),
                                     ('差分线性/hsilu_dl',['./build/hnw_endpoint','--selftest']),
                                     (paper+'/ctp_reviewer_20261117',[PY,'dl_check.py'])]),
            ('preimage_structure', [(h3,[PY,'verify_prefix.py']),
                                    (h3,['g++','-O3','-std=c++17','enumerate_branches.cpp','-o','enumerate_branches']),
                                    ('..',[PY,'scripts/check_saved_results.py','branches'])]),
            ('cubes', [('Cube/gaston_rank_certificate',['make','certificate']),
                       ('..',[PY,'scripts/check_saved_results.py','cubes']),
                       ('Cube/hsilu_quartic_screen',[PY,'verify.py']),
                       ('Cube',[PY,'src/hsilu_7r_batch_audit.py'])]),
            ('preimage_accounting', [('..',[PY,'scripts/preimage_accounting.py'])]),
        ]
    if profile == 'datf':
        return [('datf_algorithm5', [(paper,[PY,'evaluate_datf_algorithm5.py','--target','G4a']),
                                    (paper,[PY,'evaluate_datf_algorithm5.py','--target','H4a']),
                                    ('..',[PY,'scripts/check_saved_results.py','datf'])])]
    raise ValueError(profile)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=['quick','datf'], default='quick')
    parser.add_argument('--only', help='Comma-separated labels from the selected profile, for a focused rerun.')
    parser.add_argument('--timeout', type=int, default=1800, help='Seconds allowed for each command (default: 1800).')
    args = parser.parse_args()
    selected = jobs(args.profile)
    if args.only:
        names = set(args.only.split(','))
        known = {name for name, _ in selected}
        if names-known:
            parser.error('Unknown check names: '+', '.join(sorted(names-known)))
        selected = [(name, cmds) for name, cmds in selected if name in names]
    stamp = datetime.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S')
    output = PACKAGE/'validation'/(stamp+'_'+args.profile)
    output.mkdir(parents=True, exist_ok=True)
    report = {'profile': args.profile, 'started': stamp, 'python': sys.version,
              'platform': platform.platform(), 'checks': [],
              'scope': 'See README.md: quick checks do not repeat every search or full-size computation.'}
    with tempfile.TemporaryDirectory(prefix='gaston_hsilu_reproduce_') as temp:
        work = Path(temp)
        def ignore(directory, names):
            ignored = {n for n in names if n in {'__pycache__','.git','build'}}
            return ignored
        print('Preparing isolated working copy...', flush=True)
        shutil.copytree(PACKAGE/'research', work/'research', ignore=ignore)
        shutil.copytree(PACKAGE/'scripts', work/'scripts')
        shutil.copy2(PACKAGE/'EXPECTED_RESULTS.json', work/'EXPECTED_RESULTS.json')
        env = os.environ.copy()
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        for name, commands in selected:
            start = time.monotonic()
            result = {'name': name, 'status': 'PASS', 'commands': [], 'log': name+'.log'}
            print('RUN', name, flush=True)
            with (output/(name+'.log')).open('w') as log:
                for relative, command in commands:
                    log.write('COMMAND '+json.dumps(command, ensure_ascii=False)+'\n')
                    log.flush()
                    try:
                        completed = subprocess.run(command, cwd=work/'research'/relative,
                                                   stdout=log, stderr=subprocess.STDOUT,
                                                   env=env, timeout=args.timeout)
                        code = completed.returncode
                    except (subprocess.TimeoutExpired, OSError) as error:
                        log.write(str(error)+'\n')
                        code = -1
                    result['commands'].append({'cwd': relative, 'argv': command, 'exit_code': code})
                    if code:
                        result['status'] = 'FAIL'
                        break
            result['seconds'] = round(time.monotonic()-start, 3)
            report['checks'].append(result)
            print(result['status'], name, result['seconds'], 'seconds', flush=True)
        # Retain numerical output from the bounded checks, separately from original records.
        candidates = [
            'research/Gaston_Hsilu_CTP/paper_rewriting_output/ctp_reviewer_20261117/check_results.json',
            'research/Gaston_Hsilu_CTP/paper_rewriting_output/ctp_reviewer_20261117/dl_results.json',
            'research/代数攻击/Hsilu3r_new/validation_current_prefix/prefix_validation.json',
            'research/Cube/results/hsilu_7r_batch_audit.json',
            'research/Cube/gaston_rank_certificate/rerun_r4.json',
            'validation/preimage_accounting.json',
        ]
        # Only copy files produced by checks run in this invocation.
        produced = {'differential': [candidates[0]], 'differential_linear': [candidates[1]],
                    'preimage_structure': [candidates[2]], 'cubes': candidates[3:5],
                    'preimage_accounting': [candidates[5]]}
        for item in report['checks']:
            if item['status'] != 'PASS': continue
            for name in produced.get(item['name'], []):
                if (work/name).is_file():
                    shutil.copy2(work/name, output/Path(name).name)
    report['status'] = 'PASS' if all(x['status']=='PASS' for x in report['checks']) else 'FAIL'
    (output/'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(report['status'], 'Report:', output/'summary.json', flush=True)
    return 0 if report['status']=='PASS' else 1

if __name__ == '__main__':
    raise SystemExit(main())
