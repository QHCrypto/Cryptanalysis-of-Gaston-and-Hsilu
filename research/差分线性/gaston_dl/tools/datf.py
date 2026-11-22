#!/usr/bin/env python3
"""Ordinary DATF: Boolean ATFs, original isolated-variable bias algorithm.

No enhanced transition threshold, latest-variable backtracking, or input-space
partitioning is used. Correlations (twice biases) are stored as exact Fractions.
Transition rules follow QHCrypto/DATF 4rascon.sage; the estimation loop follows
the original algorithm reproduced as Algorithm 1 in ePrint 2025/1632, Section 2.4.
The algebraic DAG is exact. Independence assumptions enter only its estimation.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache

ZERO = frozenset()
ONE = frozenset((0,))


def variable(index):
    return frozenset((1 << index,))


def add(*polys):
    result = set()
    for p in polys:
        result.symmetric_difference_update(p)
    return frozenset(result)


def mul(a, b):
    if not a or not b:
        return ZERO
    if a == ONE:
        return b
    if b == ONE:
        return a
    if a == b:  # Boolean idempotence applies to every polynomial.
        return a
    result = set()
    for x in a:
        for y in b:
            term = x | y
            if term in result:
                result.remove(term)
            else:
                result.add(term)
    return frozenset(result)


def support(p):
    result = 0
    for term in p:
        result |= term
    return result


def indices(bits):
    while bits:
        bit = bits & -bits
        yield bit.bit_length() - 1
        bits ^= bit


def isolated(p):
    linear = nonlinear = 0
    for term in p:
        if term and term & (term - 1) == 0:
            linear |= term
        elif term:
            nonlinear |= term
    return linear & ~nonlinear


def condition(p, index, value):
    bit = 1 << index
    result = set()
    for term in p:
        if term & bit:
            if not value:
                continue
            term ^= bit
        if term in result:
            result.remove(term)
        else:
            result.add(term)
    return frozenset(result)


def evaluate(p, values):
    return sum((term & values) == term for term in p) & 1


def components(p):
    """Exact separation on disjoint variable sets, as in ordinary Bias()."""
    parent = {}
    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for term in p:
        vs = list(indices(term))
        for v in vs:
            parent.setdefault(v, v)
        for v in vs[1:]:
            parent[find(v)] = find(vs[0])
    groups = defaultdict(set)
    for term in p:
        if term:
            groups[find((term & -term).bit_length() - 1)].add(term)
    result = [frozenset(v) for _, v in sorted(groups.items())]
    if 0 in p:
        if result:
            result[0] = add(result[0], ONE)
        else:
            result = [ONE]
    return result


class BudgetExceeded(RuntimeError):
    pass


class OrdinaryDATF:
    def __init__(self, max_bias_nodes=1000000):
        self.definitions = []
        self.metadata = []
        self.beta_lookup = {}
        self.marginals = {}
        self.estimate_cache = {}
        self.stats = Counter()
        self.zero_witnesses = []
        self.max_bias_nodes = max_bias_nodes

    def fresh(self, expression=None, kind='input', stage=0, row=None, column=None):
        index = len(self.definitions)
        self.definitions.append(expression)
        self.metadata.append({'kind': kind, 'stage': stage, 'row': row, 'column': column})
        return variable(index)

    def transition(self, p, kind, stage, row, column):
        # Original public implementation: keep <=1-variable expressions;
        # reuse beta expressions up to XOR with 1 (no subset matching).
        if kind == 'beta':
            nonconstant = p - ONE
            if nonconstant in self.beta_lookup:
                self.stats['beta_equal_or_complement_reuses'] += 1
                return add(self.beta_lookup[nonconstant], ONE if 0 in p else ZERO)
        if support(p).bit_count() <= 1:
            return p
        out = self.fresh(p, kind, stage, row, column)
        if kind == 'beta':
            self.beta_lookup[p - ONE] = add(out, ONE if 0 in p else ZERO)
        return out

    def marginal(self, index):
        if index not in self.marginals:
            meta = self.metadata[index]
            if meta['kind'] in ('input', 'alpha'):
                # Each alpha slice is a uniform permutation state; DATF
                # additionally treats it as independent of other variables.
                result = Fraction(0)
            else:
                result = self.estimate(self.definitions[index], label=f'var_{index}')
            self.marginals[index] = result
        return self.marginals[index]

    def independent_correlation(self, p):
        """Exact weighted expectation under independent variable marginals.

        Ordinary Separate + Shannon decomposition; no enhanced preprocessing.
        Deterministic marginals and isolated uniform variables are exact shortcuts.
        """
        probs = {v: self.marginal(v) for v in indices(support(p))}

        @lru_cache(maxsize=100000)
        def calc(f):
            self.stats['bias_nodes'] += 1
            if self.stats['bias_nodes'] > self.max_bias_nodes:
                raise BudgetExceeded(f'ordinary Bias exceeded {self.max_bias_nodes} nodes')
            if not f:
                return Fraction(1)
            if f == ONE:
                return Fraction(-1)
            for v in indices(support(f)):
                if abs(probs[v]) == 1:
                    return calc(condition(f, v, int(probs[v] == -1)))
            lone = isolated(f)
            if lone:
                factor = Fraction(1)
                for v in indices(lone):
                    factor *= probs[v]
                if not factor:
                    return Fraction(0)
                return factor * calc(f - frozenset(1 << v for v in indices(lone)))
            parts = components(f)
            if len(parts) > 1:
                result = Fraction(1)
                for part in sorted(parts, key=lambda x: (support(x).bit_count(), len(x))):
                    result *= calc(part)
                    if not result:
                        break
                return result
            counts = Counter(v for term in f for v in indices(term))
            v = max(counts, key=lambda x: (counts[x], -x))
            return ((1 + probs[v]) * calc(condition(f, v, 0))
                    + (1 - probs[v]) * calc(condition(f, v, 1))) / 2

        return calc(p)

    def estimate(self, polynomial, label=None, trace=None):
        """Original Algorithm 1, expressed in normalized correlations.

        Split e=e_nonisolated+e_isolated. Multiply E(-1)^e_nonisolated,
        then substitute all isolated transitional variables simultaneously.
        Repeated factors are assumed independent even if their definitions share
        ancestors. This is deliberately NOT enhanced Algorithm 2.
        """
        if polynomial in self.estimate_cache and trace is None:
            return self.estimate_cache[polynomial]
        e, result, step = polynomial, Fraction(1), 0
        while e:
            lone = isolated(e)
            lone_terms = frozenset(1 << v for v in indices(lone))
            rest = e - lone_terms
            factor = self.independent_correlation(rest)
            self.stats['estimation_steps'] += 1
            record = {'step': step, 'terms': len(e), 'variables': support(e).bit_count(),
                      'isolated_variables': list(indices(lone)), 'factor_terms': len(rest),
                      'factor_variables': support(rest).bit_count(),
                      'factor_numerator': str(factor.numerator), 'factor_denominator': str(factor.denominator)}
            result *= factor
            if not result:
                record['zero_reason'] = 'zero nonisolated-factor correlation under independent marginals'
                if len(self.zero_witnesses) < 32:
                    self.zero_witnesses.append({'label': label, **record})
                if trace is not None:
                    trace.append(record)
                break
            replacement = ZERO
            for v in indices(lone):
                expr = self.definitions[v]
                if expr is None:
                    result = Fraction(0)
                    record['zero_reason'] = 'isolated uniform input variable in an independence-separated factor'
                    record['zero_variable'] = v
                    break
                replacement = add(replacement, expr)
            if trace is not None:
                trace.append(record)
            if not result:
                break
            e = replacement
            step += 1
        self.estimate_cache[polynomial] = result
        return result

    def evaluate_dag(self, input_values):
        assignment = 0
        for v, expr in enumerate(self.definitions):
            value = input_values[v] if expr is None else evaluate(expr, assignment)
            assignment |= value << v
        return assignment

    def serialized_dag(self):
        return [{'id': i, **m, 'expression': None if p is None else
                 [list(indices(t)) for t in sorted(p)]}
                for i, (p, m) in enumerate(zip(self.definitions, self.metadata))]
