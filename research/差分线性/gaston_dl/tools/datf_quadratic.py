"""Exact quadratic Bias backend for unchanged ordinary DATF Algorithm 1.

No different transitions, dependency assumptions, partition estimator, or
backtracking are introduced. Only the independent expectation is evaluated by
quadratic elimination when its marginals are uniform, and exact conditioning
otherwise. This backend is separately compared to ordinary Shannon Bias.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache

from datf import OrdinaryDATF, BudgetExceeded, ONE, condition, components, indices, isolated, support


def uniform_quadratic(p):
    rows={};linear=0;constant=int(0 in p)
    for term in p:
        n=term.bit_count()
        if n>2:raise ValueError('non-quadratic polynomial')
        if n==1:
            linear^=term;rows.setdefault(term,0)
        elif n==2:
            a=term&-term;b=term^a
            rows[a]=rows.get(a,0)^b;rows[b]=rows.get(b,0)^a
    weight=0
    while rows:
        edge=next(((a,adj&-adj) for a,adj in rows.items() if adj),None)
        if edge is None:return Fraction(0) if linear else Fraction(-1 if constant else 1,1<<weight)
        a,b=edge;A=rows[a]^b;B=rows[b]^a
        la=bool(linear&a);lb=bool(linear&b)
        constant^=la and lb
        linear^=(B if la else 0)^(A if lb else 0)^(A&B)
        linear&=~(a|b)
        both=a|b
        for v in indices(A|B):
            bit=1<<v
            rows[bit]^=(B if A&bit else 0)^(A if B&bit else 0)
            rows[bit]&=~(both|bit)
        del rows[a];del rows[b];weight+=1
        if any((linear&bit) and not adj for bit,adj in rows.items()):return Fraction(0)
    return Fraction(-1 if constant else 1,1<<weight)


class QuadraticOrdinaryDATF(OrdinaryDATF):
    def independent_correlation(self,p):
        if any(t.bit_count()>2 for t in p):return super().independent_correlation(p)
        probs={v:self.marginal(v) for v in indices(support(p))}

        @lru_cache(maxsize=20000)
        def calc(f):
            self.stats['bias_nodes']+=1
            if self.stats['bias_nodes']>self.max_bias_nodes:
                raise BudgetExceeded(f'exact quadratic Bias exceeded {self.max_bias_nodes} nodes')
            variables=list(indices(support(f)))
            for v in variables:
                if abs(probs[v])==1:return calc(condition(f,v,int(probs[v]==-1)))
            biased=[v for v in variables if probs[v]]
            if not biased:
                self.stats['uniform_quadratic_calls']+=1
                return uniform_quadratic(f)
            lone=isolated(f)
            if lone:
                factor=Fraction(1)
                for v in indices(lone):factor*=probs[v]
                if not factor:return Fraction(0)
                return factor*calc(f-frozenset(1<<v for v in indices(lone)))
            parts=components(f)
            if len(parts)>1:
                result=Fraction(1)
                for part in sorted(parts,key=lambda x:(support(x).bit_count(),len(x))):
                    result*=calc(part)
                    if not result:break
                return result
            counts=Counter(v for t in f for v in indices(t))
            # Once these biased coordinates are fixed, the remaining sum is
            # uniform quadratic and can be eliminated without Shannon trees.
            v=max(biased,key=lambda v:(counts[v],-v))
            return ((1+probs[v])*calc(condition(f,v,0))+(1-probs[v])*calc(condition(f,v,1)))/2
        return calc(p)
