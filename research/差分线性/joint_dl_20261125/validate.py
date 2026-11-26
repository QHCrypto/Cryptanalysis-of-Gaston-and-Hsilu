"""Independent algebra checks and agreement with pre-existing C++ cores."""
import json
import random
import subprocess
from pathlib import Path
import os

os.chdir(Path(__file__).resolve().parent)
from fractions import Fraction
from joint_estimator import *

rng = random.Random(20261125)
checks = {}
def state():
    return tuple(rng.getrandbits(64) for _ in range(5))

# Direct enumeration, independent of the Gaussian elimination implementation.
for trial in range(600):
    n = rng.randrange(1, 11)
    rows = [rng.getrandbits(n) for _ in range(rng.randrange(13))]
    rhs = [rng.randrange(2) for _ in rows]
    s, t = rng.getrandbits(n), rng.randrange(2)
    brute = Fraction(sum(1-2*(((s & x).bit_count() & 1)^t)
                         for x in range(1 << n)
                         if all(((a & x).bit_count() & 1) == b
                                for a,b in zip(rows,rhs))), 1 << n)
    actual, _, _ = affine_character(rows, rhs, s, t)
    assert actual == brute, (rows,rhs,s,t,actual,brute)
checks['gaussian_vs_exhaustive_small_systems'] = 600

# Complete matrix transpose identity for both 320-bit linear maps.
units = [embed(j, 1 << i) for i in range(5) for j in range(64)]
for L,LT in [(lg,lgt),(lh,lht)]:
    forward = [L(e) for e in units]
    backward = [LT(e) for e in units]
    for i in range(320):
        for j in range(320):
            assert dot(forward[i],units[j]) == dot(units[i],backward[j])
checks['linear_transpose_basis_pairs'] = 2*320*320

for _ in range(300):
    x,d,w = state(),state(),state()
    assert add(chi(x),chi(add(x,d))) == add(polar(x,d),chi(d),chi(ZERO))
    assert dot(w,polar(x,d)) == dot(polar_transpose(w,d),x)
checks['quadratic_derivative_and_adjoint_trials'] = 300

# Check each affine image used by the current Markov-prefix calculation
# against direct second-round difference propagation and its adjoint.
for name in ENDPOINTS:
    L, LT = (lg, lgt) if name == 'Gaston' else (lh, lht)
    d0, _ = ENDPOINTS[name]
    j0 = next(j for j in range(64) if column(d0, j))
    first = sorted({sbox(x)^sbox(x^column(d0, j0)) for x in range(32)})
    assert len(first) == 4
    for d1 in first:
        d2 = L(embed(j0, d1))
        off = L(chi(d2))
        for _ in range(8):
            y, mask = state(), state()
            direct = L(add(chi(y), chi(add(y, d2))))
            adj = polar_transpose(LT(mask), d2)
            assert dot(mask, direct) == (dot(adj, y)^dot(mask, off))
checks['current_prefix_affine_trials'] = 2*4*8

subprocess.run(['g++','-O2','-std=c++17','reference_vectors.cpp','-o','reference_vectors'],check=True)
lines = subprocess.check_output(['./reference_vectors'],text=True).splitlines()
for line in lines:
    v = [int(s) for s in line.split()]
    x = tuple(v[:5])
    assert lg(x) == tuple(v[5:10])
    assert lh(x) == tuple(v[10:15])
    assert permutation(x,'Gaston') == tuple(v[15:20])
    assert permutation(x,'Hsilu') == tuple(v[20:25])
checks['existing_cpp_core_vectors_per_target'] = len(lines)
checks['chi_local_autocorrelation_cases'] = 32*32  # Checked on module import.

checks['status'] = 'PASS'
Path('validation.json').write_text(json.dumps(checks, indent=2)+'\n')
print(json.dumps(checks, indent=2))
