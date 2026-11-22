"""Hsilu two-chi HNW zero model, exact affine/quadratic first-DDT sum.

Computes E_X prod_j DLCT((L D_delta chi(X))_j,gamma_j).
This is a MODEL expectation, not an exact two-round permutation expectation.
Equivalent algebra to the verified Gaston backend, with Hsilu's L substituted.
"""
import sys
from pathlib import Path
from functools import lru_cache
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'hsilu_dl/tools'))
sys.path.insert(0,str(ROOT/'gaston_dl/tools'))
from hsilu import linear,transpose,chi,embed
from datf_quadratic import uniform_quadratic

def bits(n):
    while n:
        b=n&-n;yield b;n^=b

@lru_cache(maxsize=64)
def affine(delta):
    base=linear(chi(delta))
    d=[(base[i]>>j)&1 for i in range(5) for j in range(64)]
    for k in range(5):
        for j in range(64):
            pattern=(((delta[(k+1)%5]>>j)&1)<<((k-1)%5))^(((delta[(k-1)%5]>>j)&1)<<((k-2)%5))
            if not pattern:continue
            image=linear(embed(j,pattern));var=1<<(k*64+j+1)
            for i in range(5):
                for b in bits(image[i]):d[i*64+b.bit_length()-1]^=var
    return tuple(d)

@lru_cache(maxsize=1024)
def equations(gamma):
    eqs=[]
    for j in range(64):
        for k in range(5):
            a=((k+1)%5)*64+j if gamma[(k-1)%5]>>j&1 else -1
            b=((k-1)%5)*64+j if gamma[(k-2)%5]>>j&1 else -1
            if a!=-1 or b!=-1:eqs.append((a,b))
    return eqs

def two(delta,gamma):
    d=affine(delta);pivots={}
    for a,b in equations(gamma):
        eq=(d[a] if a>=0 else 0)^(d[b] if b>=0 else 0)
        while eq&~1:
            pivot=1<<(eq.bit_length()-1)
            if pivot in pivots:eq^=pivots[pivot]
            else:pivots[pivot]=eq;break
        else:
            if eq==1:return Fraction(0)
    substitutions={}
    def subst(f):
        out=f&1
        for b in bits(f&~1):out^=substitutions[b]
        return out
    for k in range(1,321):
        b=1<<k;substitutions[b]=subst(pivots[b]^b) if b in pivots else b
    def poly(f):return set([0] if f&1 else [])|set(bits(f&~1))
    def multiply(a,b):
        out=set()
        for x in poly(a):
            for y in poly(b):
                t=x|y
                if t in out:out.remove(t)
                else:out.add(t)
        return out
    p=set()
    for i in range(5):
        for bit in bits(gamma[i]):
            j=bit.bit_length()-1
            a=subst(d[i*64+j]);b=subst(d[((i+1)%5)*64+j]);c=subst(d[((i+2)%5)*64+j])
            p.symmetric_difference_update(poly(a^c));p.symmetric_difference_update(multiply(b,c))
    return uniform_quadratic(frozenset(p))/2**len(pivots)
