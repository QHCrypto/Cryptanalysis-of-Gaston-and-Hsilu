"""Exact highest-degree obstruction to extending the fixed 21-variable path.

Polynomials live in the Boolean ring: monomials are variable-index bitmasks.
Only the homogeneous parts of degrees 1, 2, and 4 are propagated.
The input capacity conditions make the first chi's active part row 0 only.
Consequently these highest-degree coefficients are independent of background
values and round constants, including the optional X2[1][17] condition.
"""
import ctypes
import itertools
import json
import random
import verify_path_algebra as v

U = v.DATA['free_variables_by_lane'][0]

def rol_poly(row, n):
    return [row[(j-n) % 64] for j in range(64)]

def linear_poly(a):
    x = list(a)
    for j in range(5,18):
        left = x[v.PARAM['n1'][j]]
        right = rol_poly(x[v.PARAM['n2'][j]], v.PARAM['m'][j])
        x.append([s ^ t for s,t in zip(left,right)])
    return [rol_poly(x[13+i], v.PARAM['s'][i]) for i in range(5)]

def multiply_top(p, q, degree):
    out = set()
    for a in p:
        for b in q:
            monomial = a | b
            if monomial.bit_count() == degree:
                if monomial in out:
                    out.remove(monomial)
                else:
                    out.add(monomial)
    return out

def chi_top(a, degree):
    return [[multiply_top(a[(i+1)%5][j], a[(i+2)%5][j], degree)
             for j in range(64)] for i in range(5)]

first_chi = [[set() for _ in range(64)] for _ in range(5)]
for k,j in enumerate(U):
    first_chi[0][63-j] = {1 << k}
x2_linear = linear_poly(first_chi)
x3_quadratic = linear_poly(chi_top(x2_linear, 2))
x4_quartic = linear_poly(chi_top(x3_quadratic, 4))
monomials = sorted(set().union(*x4_quartic[0]))
index = {m:k for k,m in enumerate(monomials)}
rows = [sum(1 << index[m] for m in polynomial) for polynomial in x4_quartic[0]]
quartic_rank = v.rank(rows)
degree4_counts = [[len(a[63-j]) for j in range(64)] for a in x4_quartic]

# Independent coefficient check: a 4th finite difference of the C++ reference
# returns exactly that degree-4 ANF coefficient (the total degree is at most 4).
lib = ctypes.CDLL(str(v.HERE / 'verify_reference_wrapper.so'))
lib.verify_round.argtypes = [ctypes.POINTER(ctypes.c_uint64), ctypes.c_uint64]
random.seed(20261013)
umask = sum(1 << (63-j) for j in U)
checks = 0
for sample in range(3):
    background = [random.getrandbits(64) for _ in range(5)]
    background[0] &= ~umask
    background[1] &= ~umask
    background[4] |= umask
    for variable_indices in random.sample(list(itertools.combinations(range(21),4)), 20):
        monomial = sum(1 << k for k in variable_indices)
        numerical = [0]*5
        for assignment in range(16):
            a = background.copy()
            for bit,k in enumerate(variable_indices):
                if assignment >> bit & 1:
                    a[0] ^= 1 << (63-U[k])
            c = (ctypes.c_uint64*5)(*a)
            for rc in (0x69,0x5a,0x4b):
                lib.verify_round(c,rc)
            numerical = [s ^ t for s,t in zip(numerical,c)]
        expected = [sum((monomial in polynomial) << j for j,polynomial in enumerate(row)) for row in x4_quartic]
        assert numerical == expected
        checks += 1

out = {'scope': 'fixed 21 active message positions and the 42 original capacity conditions',
       'path_to_native_bit': '63-j', 'numerical_round_constants': ['0x69','0x5a','0x4b'],
       'quadratic_free_first_block_positions_after_2_rounds':
           [j for j in range(64) if not x3_quadratic[0][63-j]],
       'quartic_free_first_block_positions_after_3_rounds':
           [j for j in range(64) if not x4_quartic[0][63-j]],
       'quartic_coefficient_matrix_rank_first_block': quartic_rank,
       'quartic_coefficient_matrix_monomial_count': len(monomials),
       'dimension_of_first_block_output_masks_cancelling_all_quartics': 64-quartic_rank,
       'degree4_monomial_counts_by_path_row_and_bit': degree4_counts,
       'cpp_reference_coefficient_checks': checks,
       'conclusion_scope': 'No background assignment can remove these highest-degree terms without changing/restricting the active variable space. Does not exclude new three-round paths or approximation methods.'}
(v.HERE / 'three_round_extension_exact_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:w for k,w in out.items() if k!='degree4_monomial_counts_by_path_row_and_bit'},indent=2))
