"""Independent reference and MP/MILP parity checks; solver logs stay visible."""
import ctypes
from datetime import datetime
import uuid
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / 'results' / 'validation' / (datetime.now().strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8])
sys.path.insert(0, str(ROOT / 'reference'))
from gaston import IV, coefficient, linear, permute


def run(case, tag, pool=256, expand=1, checkpoint=None, resume=False, key_split=1):
    path = VALIDATION / (tag + '.json')
    cmd = [str(ROOT / 'gaston_cube'), '--rounds', str(case['rounds']),
           '--width', str(case.get('width', 64)), '--key', str(case['key']),
           '--nonce', ','.join(map(str, case['nonce'])),
           '--output-bit', str(case['output']), '--mode', 'count',
           '--pool', str(pool), '--expand-first', str(expand), '--key-split', str(key_split),
           '--result', str(path)]
    if checkpoint:
        cmd += ['--checkpoint', str(checkpoint)]
    if resume:
        cmd += ['--resume']
    subprocess.run(cmd, check=True)
    result = json.loads(path.read_text())
    expected = coefficient(case['rounds'], case['key'], case['nonce'],
                           case['output'], case.get('width', 64))
    assert result['complete'] and result['coefficient'] == expected, (case, result, expected)
    return {**case, 'expected': expected, 'trails': result['trails_counted'],
            'classification': result['classification'], 'expand_first': expand, 'key_split': key_split,
            'result': str(path.relative_to(ROOT))}


def main():
    folder = VALIDATION
    folder.mkdir(parents=True)
    official = ROOT.parent / 'gaston_ref/gaston/implementation/reference/gaston.c'
    shared = folder / 'libgaston_reference.so'
    subprocess.run(['gcc', '-shared', '-fPIC', '-O2', str(official), '-o', str(shared)], check=True)
    lib = ctypes.CDLL(str(shared))
    class State(ctypes.Structure):
        _fields_ = [('x', ctypes.c_uint64 * 5)]
    lib.gaston.argtypes = [ctypes.POINTER(State)]
    lib.gaston_round.argtypes = [ctypes.POINTER(State), ctypes.c_uint64]
    rng = random.Random(2026102263)
    for _ in range(64):
        s = [rng.getrandbits(64) for _ in range(5)]
        a = State((ctypes.c_uint64 * 5)(*s))
        lib.gaston(ctypes.byref(a))
        assert list(a.x) == permute(s, 12)
        a = State((ctypes.c_uint64 * 5)(*s))
        for rc in (0x96, 0x87, 0x78, 0x69, 0x5A, 0x4B):
            lib.gaston_round(ctypes.byref(a), rc)
        assert list(a.x) == permute(s, 6)
    # Check independently expanded matrix rows against the reference linear layer.
    from gaston import E, W, T
    for bit in range(320):
        s = [0]*5
        s[bit//64] = 1 << (bit % 64)
        out = [0]*5
        for i in range(5):
            for b in range(64):
                q = (b-W[i]) % 64
                v = (s[i] >> ((q-E[i]) % 64)) & 1
                for j in range(5):
                    for shift in (23, 24, 23+T[j], 41+T[j]):
                        v ^= (s[j] >> ((q-E[j]-shift) % 64)) & 1
                out[i] |= v << b
        assert out == linear(s)
    cases = []
    # Include nonzero examples selected by the independent Boolean reference.
    for rounds in (1, 2):
        found = {0: 0, 1: 0}
        while min(found.values()) < 3:
            case = dict(rounds=rounds, key=rng.randrange(128),
                        nonce=sorted(rng.sample(range(128), (1 << rounds)-1)),
                        output=rng.randrange(128))
            value = coefficient(rounds, case['key'], case['nonce'], case['output'])
            if found[value] < 3:
                cases.append(case)
                found[value] += 1
    cases += [
        dict(rounds=3, key=106, nonce=[12,17,18,26,42,101,112], output=127),
        dict(rounds=3, key=97, nonce=[22,38,50,58,89,93,122], output=14),
        dict(rounds=3, key=35, nonce=[14,27,43,46,76,85,105], output=115),
        dict(rounds=3, key=118, nonce=[5,21,40,45,49,81,125], output=65),
    ]
    records = []
    for i, case in enumerate(cases):
        for key_split in (0, 1):
            tag = f'case_{i:02d}_split_{key_split}'
            records.append(run(case, tag, key_split=key_split))
            if case['rounds'] <= 2:
                records.append(run(case, tag + '_raw', expand=0, key_split=key_split))
    # A many-trail example makes batching/resume and GF(2) cancellation substantive.
    batch_case = cases[12]
    checkpoint = folder / 'resume.trails'
    if checkpoint.exists():
        raise RuntimeError('validation resume checkpoint already exists; preserve it and use a fresh validation directory')
    before = run(batch_case, 'batch', pool=31, checkpoint=checkpoint)
    after = run(batch_case, 'resume', pool=31, checkpoint=checkpoint, resume=True)
    assert before['trails'] == after['trails'] and before['trails'] > 31
    # Old checkpoints must remain reusable with the old circuit, and incompatible
    # split circuits must reject them instead of silently mixing path parities.
    legacy_checkpoint = folder / 'legacy.trails'
    legacy = run(batch_case, 'legacy_checkpoint', checkpoint=legacy_checkpoint, key_split=0)
    assert legacy['trails'] == 169
    mismatch = subprocess.run([str(ROOT/'gaston_cube'), '--rounds','3','--key','106',
        '--nonce','12,17,18,26,42,101,112','--output-bit','127','--mode','count',
        '--key-split','1','--checkpoint',str(legacy_checkpoint),'--resume',
        '--result',str(folder/'mismatch.json')])
    assert mismatch.returncode == 1 and not (folder/'mismatch.json').exists()
    # Constants disappear from the saturated-degree component.
    constant_checks = 0
    for case in cases[12:]:
        expected = coefficient(case['rounds'], case['key'], case['nonce'], case['output'])
        for iv, start in [(0, 0), (rng.getrandbits(64), 4)]:
            assert coefficient(case['rounds'], case['key'], case['nonce'], case['output'], iv=iv, start=start) == expected
            constant_checks += 1
    expected_joint = [coefficient(2,k,[38,46,122],56) for k in range(128)]
    for key_split in (0, 1):
        joint_result = folder / f'all_keys_r2_split_{key_split}.json'
        subprocess.run([str(ROOT/'gaston_cube'), '--rounds','2','--key','all',
                        '--nonce','38,46,122','--output-bit','56','--mode','count',
                        '--key-split',str(key_split),'--result',str(joint_result)], check=True)
        joint = json.loads(joint_result.read_text())
        assert joint['complete'] and [joint['coefficients'][str(k)] for k in range(128)] == expected_joint
    summary = dict(all_passed=True, full_reference_checks=128,
                   linear_basis_checks=320, coefficient_checks=len(records),
                   nonzero_checks=sum(r['expected'] for r in records),
                   cancellation_checks=sum(r['expected']==0 and r['trails']>0 for r in records),
                   constants_checks=constant_checks, batching_resume_passed=True,
                   legacy_checkpoint_rejection_passed=True, all_key_vector_checks=256,
                   all_key_vector_nonzero=sum(expected_joint),
                   scope='Actual 64-bit words, 1-3 rounds; does not certify any six-round coefficient.',
                   cpp_sha256=hashlib.sha256((ROOT/'gaston_cube.cpp').read_bytes()).hexdigest(),
                   cpp_core_sha256=hashlib.sha256((ROOT/'gaston_mp.hpp').read_bytes()).hexdigest(),
                   official_c_sha256=hashlib.sha256(official.read_bytes()).hexdigest(), records=records)
    (ROOT/'results/validation_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='records'}, indent=2))


if __name__ == '__main__':
    main()
