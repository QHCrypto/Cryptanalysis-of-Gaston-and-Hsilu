"""Exact first-two-round joint difference law plus a Markov tail.

The final two rounds retain joint masks/differences but use fresh uniform
base states. The computed rational is exact for this approximation, NOT an
exact correlation certificate for the fixed four-round permutation.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
import argparse
import json
import math
import time
from core import *


def span_basis(values):
    b = {}
    for v in values:
        while v:
            p = v.bit_length()-1
            if p in b:
                v ^= b[p]
            else:
                b[p] = v
                break
    return list(b.values())


SBOX = [sbox(x) for x in range(32)]
LAT = [[sum(1-2*(((a & x).bit_count() ^ (b & SBOX[x]).bit_count()) & 1)
            for x in range(32)) for a in range(32)] for b in range(32)]
SUPPORT = [[a for a, value in enumerate(row) if value] for row in LAT]
BASES = [span_basis(a ^ support[0] for a in support) for support in SUPPORT]
for b in range(32):
    assert len(SUPPORT[b]) == 1 << len(BASES[b])
    assert len({LAT[b][a]**2 for a in SUPPORT[b]}) == 1
    for d in range(32):
        direct = Fraction(sum(1-2*((b & (SBOX[x]^SBOX[x^d])).bit_count()&1)
                              for x in range(32)), 32)
        form = (0 if any((v & d).bit_count() & 1 for v in BASES[b])
                else 1-2*((SUPPORT[b][0] & d).bit_count() & 1))
        assert direct == form


def affine_character(rows, rhs, sign_row, sign_offset):
    """E[1_{rows*x=rhs} (-1)^(sign_row*x+sign_offset)], uniform x."""
    pivots = {}
    for row, bit in zip(rows, rhs):
        while row:
            p = row.bit_length()-1
            if p in pivots:
                row ^= pivots[p][0]
                bit ^= pivots[p][1]
            else:
                pivots[p] = (row, bit)
                break
        if not row and bit:
            return Fraction(0), len(pivots), 'inconsistent'
    while sign_row:
        p = sign_row.bit_length()-1
        if p not in pivots:
            return Fraction(0), len(pivots), 'balanced_sign'
        row, bit = pivots[p]
        sign_row ^= row
        sign_offset ^= bit
    return Fraction(1-2*sign_offset, 1 << len(pivots)), len(pivots), 'nonzero'


def fraction_record(x):
    return {'numerator': x.numerator, 'denominator': x.denominator,
            'value': float(x), 'log2_abs': math.log2(abs(x)) if x else None}


def estimate(name, prefix='exact', factorize_third=False, target_override=None):
    start = time.perf_counter()
    L, LT = (lg, lgt) if name == 'Gaston' else (lh, lht)
    delta, target = ENDPOINTS[name]
    if target_override is not None:
        target = target_override
    active_input = [j for j in range(64) if column(delta, j)]
    assert len(active_input) == 1
    j0 = active_input[0]
    d0 = column(delta, j0)
    free_mask = ((1 << 320)-1) ^ pack(embed(j0, 31))
    active_out = [j for j in range(64) if column(target, j)]
    last_branches = list(product(*(SUPPORT[column(target, j)] for j in active_out)))
    groups = {}
    for x0 in range(32):
        d1 = embed(j0, SBOX[x0] ^ SBOX[x0 ^ d0])
        groups.setdefault(d1, []).append(x0)
    total = Fraction(0)
    details = []
    status_counts = Counter()
    for choices in last_branches:
        alpha4 = ZERO
        weight4 = Fraction(1)
        for j, a in zip(active_out, choices):
            alpha4 = add(alpha4, embed(j, a))
            weight4 *= Fraction(LAT[column(target,j)][a]**2, 1024)
        b3 = LT(alpha4)
        constraints = [embed(j, v) for j in range(64) for v in BASES[column(b3,j)]]
        origin = ZERO
        for j in range(64):
            origin = add(origin, embed(j, SUPPORT[column(b3,j)][0]))
        branch = Fraction(0)
        term_details = []
        for d1, xs in groups.items():
            d2 = L(d1)
            def adjoint(mask):
                v = polar_transpose(LT(mask), d2)
                if prefix == 'fresh':
                    return pack(v)
                return pack(LT(v)) & free_mask
            rows = [adjoint(k) for k in constraints]
            sign_row = adjoint(origin)
            for x0 in xs:
                u0 = embed(j0, SBOX[x0])
                y2 = ZERO if prefix == 'fresh' else add(L(u0), (RC[9],0,0,0,0))
                offset = L(add(chi(y2), chi(add(y2,d2))))
                rhs = [dot(k, offset) for k in constraints]
                value, rank, reason = affine_character(rows, rhs, sign_row, dot(origin, offset))
                if factorize_third:
                    value = Fraction(1)
                    rank = 0
                    for j in range(64):
                        b = column(b3, j)
                        if not b:
                            continue
                        ks = [embed(j, v) for v in BASES[b]]
                        o = embed(j, SUPPORT[b][0])
                        local, local_rank, _ = affine_character(
                            [adjoint(k) for k in ks], [dot(k, offset) for k in ks],
                            adjoint(o), dot(o, offset))
                        value *= local
                        rank += local_rank
                    reason = 'nonzero' if value else 'factorized_zero'
                status_counts[reason] += 1
                branch += value / 32
                if value:
                    term_details.append({'first_column_value':x0, 'rank':rank,
                                         'conditional':fraction_record(value)})
        contribution = weight4 * branch
        total += contribution
        details.append({'last_chi_input_mask':[f'{x:016x}' for x in alpha4],
                        'last_chi_squared_weight':fraction_record(weight4),
                        'prefix_and_third_round':fraction_record(branch),
                        'contribution':fraction_record(contribution),
                        'nonzero_conditions':term_details})
    scope = ('Exact joint prefix of two chi layers, followed by independent uniform base-state approximation for chi3 and chi4.'
             if prefix == 'exact' else
             'Full-joint differential Markov model: uniform independent base states at each chi.')
    scope += (' Ablation: factorize the chi3 autocorrelation expectation across columns; retain joint last-round masks.'
              if factorize_third else
              ' No cross-column product approximation for differences.')
    return {'target':name, 'delta':[f'{x:016x}' for x in delta],
            'lambda':[f'{x:016x}' for x in target],
            'correlation':fraction_record(total), 'final_mask_branches':len(last_branches),
            'status_counts':dict(status_counts), 'branches':details,
            'seconds':time.perf_counter()-start, 'prefix':prefix,
            'factorize_third':factorize_third,
            'scope':scope}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--target', choices=['Gaston','Hsilu','both'],default='both')
    p.add_argument('--output', required=True)
    p.add_argument('--prefix',choices=['exact','fresh'],default='exact')
    p.add_argument('--factorize-third',action='store_true')
    args = p.parse_args()
    results=[]
    for name in (['Gaston','Hsilu'] if args.target=='both' else [args.target]):
        result=estimate(name,args.prefix,args.factorize_third)
        results.append(result)
        print(name, result['correlation'], result['seconds'], flush=True)
    with open(args.output,'w') as stream:
        json.dump({'results':results},stream,indent=2)
