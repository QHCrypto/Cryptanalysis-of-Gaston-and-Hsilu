#!/usr/bin/env python3
"""Compare fresh bounded computations with the archived paper records."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

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
    expected = json.loads((Path(__file__).resolve().parents[1]/'EXPECTED_RESULTS.json').read_text())['datf_legacy']
    for target in ['G4a','H4a']:
        result = json.loads((directory/(target+'.json')).read_text())
        assert result['status'] == 'completed', result.get('reason')
        assert result['correlation'] == expected[target], (target, result['correlation'])
        print(target, result['status'], result['correlation'])

def hsilu_search():
    """Rerun the C++ extension search on the selected manuscript core."""
    directory = ROOT/'Hsilu_trail/search_below82_20261117'
    saved = json.loads((directory/'results/best_core.json').read_text())
    with tempfile.TemporaryDirectory(prefix='hsilu_fixed_core_') as temp:
        work = Path(temp)
        (work/'results').mkdir()
        (work/'seeds.txt').write_text(' '.join(saved['a1'])+'\n')
        output = subprocess.check_output([str(directory/'search'), 'seeds'], cwd=work, text=True)
        stats = json.loads(output)
        candidates = [json.loads(line) for line in (work/'results/seeds.jsonl').read_text().splitlines()]
        assert stats['extended_cores'] == 1 and stats['extensions'] == 16384, stats
        assert stats['best'] == 80, stats
        assert any(all(row[key] == saved[key] for key in saved) for row in candidates)
    result = {'status': 'PASS', 'scope': 'Complete extension of the saved core, not a global search',
              'search': stats, 'weights': saved['weights'], 'saved_states_reproduced': True}
    out = ROOT.parent/'validation/hsilu_search.json'
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('check', choices=['cubes','branches','datf','hsilu_search'])
    args = parser.parse_args()
    globals()[args.check]()
