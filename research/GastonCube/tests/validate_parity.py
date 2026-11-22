"""Independent Boolean tests of the report cube, MP parity and rotations."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
sys.path.insert(0,str(ROOT/'src'))
from gaston import rol
from parity_cube import cube_sum, gf2_rank, matrix_from_columns, matrix_from_rows, nonce_word
from run_parity_matrix import summarize


def main():
    folder = ROOT/'results/parity_validation'/(datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8])
    folder.mkdir(parents=True)
    rng = random.Random(2026102264)
    references = {}
    rotation_checks = affine_checks = constants_checks = 0
    for rounds in (1,2,3):
        period = 1 << rounds
        for index in range(1 << (period-1)):
            n = nonce_word(index,rounds)
            assert n.bit_count() % 2 == 0 and rol(n,period) == n
            # A rotated block still has even parity, hence belongs to the same cube.
            assert (rol(n,1) & ((1 << period)-1)).bit_count() % 2 == 0
        base = cube_sum(0,rounds)
        cols = [cube_sum(1,rounds)^base,cube_sum(1 << 64,rounds)^base]
        matrix = matrix_from_columns(cols)
        for k in range(128):
            direct = cube_sum(1 << k,rounds)^base
            assert [(direct >> j) & 1 for j in range(128)] == [matrix[j][k] for j in range(128)]
            rotation_checks += 1
        assert matrix_from_rows([matrix[0],matrix[64]]) == matrix
        for _ in range(3):
            key = rng.getrandbits(128)
            expected = base
            for j,row in enumerate(matrix):
                expected ^= (sum(row[k] for k in range(128) if (key >> k) & 1) & 1) << j
            assert cube_sum(key,rounds) == expected
            affine_checks += 1
        for iv,start in ((0,0),(rng.getrandbits(64),4)):
            changed_base = cube_sum(0,rounds,iv=iv,start=start)
            for a in (0,1):
                assert cube_sum(1 << (64*a),rounds,iv=iv,start=start)^changed_base == cols[a]
                constants_checks += 1
        references[rounds] = dict(matrix=matrix,rank=gf2_rank(matrix),base=base,columns=cols)
    records = []
    def milp(rounds,key,output,split=1,expand=1,checkpoint=None,resume=False,pool=256):
        tag = f'r{rounds}_k{key}_o{output}_split{split}_expand{expand}'+('_resume' if resume else '')
        result = folder/(tag+'.json')
        cmd = [str(ROOT/'gaston_parity_cube'),'--rounds',str(rounds),'--key',str(key),
               '--output-bit',str(output),'--mode','count','--key-split',str(split),
               '--expand-first',str(expand),'--pool',str(pool),'--result',str(result)]
        if checkpoint:
            cmd += ['--checkpoint',str(checkpoint)]
        if resume:
            cmd += ['--resume']
        subprocess.run(cmd,check=True)
        r = json.loads(result.read_text())
        assert r['complete'] and r['cube_kind']=='even_parity_periodic' and r['nonce_indices'] is None
        keys = range(128) if key == 'all' else [key]
        for k in keys:
            assert r['coefficients'][str(k)] == references[rounds]['matrix'][output][k], (tag,k,r)
        records.append(dict(rounds=rounds,key=key,output=output,key_split=split,expand_first=expand,
                            coefficients_checked=len(keys),trails=r['trails_counted'],result=str(result.relative_to(ROOT))))
        return r
    for rounds in (1,2):
        for output in (0,64):
            milp(rounds,'all',output)
    milp(2,'all',0,split=0)
    milp(2,'all',64,expand=0,pool=8192)
    for key in (0,64):
        outputs = {0,64}
        for target_value in (0,1):
            candidates = [j for j in range(128) if references[3]['matrix'][j][key]==target_value]
            if candidates:
                outputs.add(candidates[0])
        for output in sorted(outputs):
            milp(3,key,output,pool=65536)
    # The driver must reconstruct every certified entry with the correct rotation sign.
    driver = folder/'driver'
    subprocess.run([sys.executable,str(ROOT/'src/run_parity_matrix.py'),'--rounds','2',
                    '--pool','7','--directory',str(driver)],check=True)
    matrix = json.loads((driver/'matrix.json').read_text())
    assert matrix['complete'] and matrix['matrix']==references[2]['matrix']
    assert matrix['matrix_rank']==references[2]['rank'] and matrix['constant_term'] is None
    before = milp(2,'all',0,checkpoint=folder/'resume.trails',pool=7)
    after = milp(2,'all',0,checkpoint=folder/'resume.trails',pool=7,resume=True)
    assert before['trails_counted']==after['trails_counted']
    partial = folder/'partial'
    partial.mkdir()
    (partial/'row_0.json').write_bytes((driver/'row_0.json').read_bytes())
    p = summarize(partial,2,64,1)
    assert p['unknown_entries']==8192 and p['matrix_rank'] is None and not p['complete']
    assert p['matrix'][:64]==references[2]['matrix'][:64]
    assert all(x is None for row in p['matrix'][64:] for x in row)
    subprocess.run([sys.executable,str(ROOT/'src/three_cube_reference.py'),'--rounds','3',
                    '--result',str(folder/'three_cube.json')],check=True)
    three = json.loads((folder/'three_cube.json').read_text())
    assert three['matrix']==references[3]['matrix'] and int(three['constant_term_hex'],16)==references[3]['base']
    summary = dict(all_passed=True,word_bits=64,rounds=[1,2,3],
                   independently_evaluated_key_columns=rotation_checks,random_key_affine_checks=affine_checks,
                   iv_round_constant_checks=constants_checks,
                   milp_coefficient_checks=sum(r['coefficients_checked'] for r in records),
                   matrix_driver_passed=True,partial_unknowns_passed=True,checkpoint_resume_passed=True,
                   three_cube_interpolation_passed=True,
                   three_cube_result=str((folder/'three_cube.json').relative_to(ROOT)),
                   matrix_driver_result=str((driver/'matrix.json').relative_to(ROOT)),
                   reference_ranks={str(r):v['rank'] for r,v in references.items()},
                   cpp_sha256=hashlib.sha256((ROOT/'gaston_parity_cube.cpp').read_bytes()).hexdigest(),
                   cpp_core_sha256=hashlib.sha256((ROOT/'gaston_mp.hpp').read_bytes()).hexdigest(),
                   scope='Actual 64-bit words at 1-3 rounds with periodic parity cubes; not a six-round coefficient or rank certificate.',
                   records=records)
    (ROOT/'results/parity_validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},indent=2))


if __name__ == '__main__':
    main()
