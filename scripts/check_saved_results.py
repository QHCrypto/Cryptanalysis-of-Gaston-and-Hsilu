#!/usr/bin/env python3
"""Compare fresh bounded computations with the archived paper records."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1] / 'research'

def cubes():
    directory = ROOT / 'Cube/gaston_rank_certificate'
    fresh = directory / 'rerun_r4.json'
    subprocess.run([str(directory/'certificate'), '4', '0', '32768', str(fresh)], check=True)
    sys.path.insert(0, str(directory))
    from merge import merge
    report = {}
    for rounds, source in [(4, fresh), (5, directory/'complete_r5.json')]:
        value = merge([source])
        saved = json.loads((directory/f'certificate_r{rounds}.json').read_text())
        assert value['coverage_complete']
        for key in ['sums_hex','representative_columns_hex','matrix_rank','covered_points']:
            assert value[key] == saved[key], (rounds, key)
        report[str(rounds)] = {
            'rank': value['matrix_rank'],
            'cube_sums_recomputed_this_run': rounds == 4,
            'rank_recomputed_this_run': True,
            'covered_points_in_record': value['covered_points'],
        }
    print(json.dumps(report, indent=2))

def branches():
    directory = ROOT / '代数攻击/Hsilu3r_new/validation_current_prefix'
    output = subprocess.check_output([str(directory/'enumerate_branches'), str(directory/'branch_system.txt')], text=True)
    fresh = json.loads(output)
    saved = json.loads((directory/'exhaustive_branches.json').read_text())
    for key, value in fresh.items():
        assert value == saved[key], key
    assert fresh['branches_enumerated'] == 2**23
    print(output)

def datf():
    directory = ROOT/'Gaston_Hsilu_CTP/paper_rewriting_output/dl_datf_algorithm5_20261111'
    expected = json.loads((Path(__file__).resolve().parents[1]/'EXPECTED_RESULTS.json').read_text())['datf']
    for target in ['G4a','H4a']:
        result = json.loads((directory/(target+'.json')).read_text())
        assert result['status'] == 'completed', result.get('reason')
        assert result['correlation'] == expected[target], (target, result['correlation'])
        print(target, result['status'], result['correlation'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('check', choices=['cubes','branches','datf'])
    args = parser.parse_args()
    globals()[args.check]()
