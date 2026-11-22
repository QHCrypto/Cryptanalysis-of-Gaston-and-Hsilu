"""JOC DATF Algorithm 5, with exact rational BiasP evaluation.

Input and transitional distributions are stored as correlations E[(-1)^v].
Line 6 evaluates EACH definition under the independent distributions of its
immediate predecessors. It does NOT recursively invoke Algorithm 5's peeling
loop, and it does NOT replace all alpha distributions by uniform distributions.
The DAG uses the original reference code's <=1-variable and equal/complementary
beta simplifications. Uniform/quadratic shortcuts change BiasP's evaluation
order only. Results are model estimates; resource exhaustion is never zero.
"""
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'DATF算法3复核_20261021'))
from algorithm3 import Algorithm3, ResourceLimit
from datf import ZERO, ONE, add, indices, support, isolated, components, condition


class Algorithm5:
    def __init__(self, engine, seconds=120, max_nodes=500000, max_terms=1000000):
        self.engine = engine
        self.started = time.perf_counter()
        self.seconds = seconds
        self.max_nodes = max_nodes
        self.max_terms = max_terms
        self.distributions = {}
        self.cache = {}
        self.stats = Counter()
        self.trace = []
        self.phase = 'transitional_distributions'
        self.uniform = Algorithm3(engine, seconds=seconds, max_nodes=max_nodes,
                                  max_terms=max_terms)

    def check(self, p=ZERO):
        self.stats['peak_polynomial_terms'] = max(self.stats['peak_polynomial_terms'], len(p))
        if len(p) > self.max_terms:
            raise ResourceLimit('Algorithm 5 polynomial term limit')
        if self.stats['weighted_nodes'] + self.uniform.stats['bias_nodes'] > self.max_nodes:
            raise ResourceLimit('Algorithm 5 BiasP node limit')
        if time.perf_counter() - self.started > self.seconds:
            raise ResourceLimit('Algorithm 5 wall time limit')

    def bias_p_correlation(self, f):
        """Exactly 2*BiasP(f), for the currently assigned independent marginals."""
        self.check(f)
        if f in self.cache:
            self.stats['weighted_cache_hits'] += 1
            return self.cache[f]
        self.stats['weighted_nodes'] += 1
        if not f:
            return Fraction(1)
        if f == ONE:
            return Fraction(-1)
        variables = list(indices(support(f)))
        assert all(v in self.distributions for v in variables)
        deterministic = next((v for v in variables if abs(self.distributions[v]) == 1), None)
        if deterministic is not None:
            ans = self.bias_p_correlation(condition(f, deterministic,
                                                   int(self.distributions[deterministic] == -1)))
        else:
            lone = isolated(f)
            if lone:
                factor = Fraction(1)
                for v in indices(lone):
                    factor *= self.distributions[v]
                ans = factor * self.bias_p_correlation(f - frozenset(1 << v for v in indices(lone))) if factor else Fraction(0)
            elif all(self.distributions[v] == 0 for v in variables):
                self.stats['uniform_backend_calls'] += 1
                ans = self.uniform.uniform_correlation(f)
            else:
                parts = components(f)
                if len(parts) > 1:
                    ans = Fraction(1)
                    for part in sorted(parts, key=lambda g: (len(g), support(g).bit_count())):
                        ans *= self.bias_p_correlation(part)
                        if not ans:
                            break
                else:
                    counts = Counter(v for t in f for v in indices(t))
                    v = max(variables, key=lambda i: (self.distributions[i] != 0, counts[i], -i))
                    c = self.distributions[v]
                    ans = ((1+c)*self.bias_p_correlation(condition(f, v, 0))
                           + (1-c)*self.bias_p_correlation(condition(f, v, 1))) / 2
        if len(self.cache) < 100000:
            self.cache[f] = ans
        return ans

    def initialize_distributions(self):
        # Eager topological evaluation implements Algorithm 5, lines 2 and 6.
        for v, (expr, meta) in enumerate(zip(self.engine.definitions, self.engine.metadata)):
            self.check()
            if expr is None:
                assert meta['kind'] == 'input'
                value = Fraction(0)
            else:
                assert not (support(expr) >> v)
                value = self.bias_p_correlation(expr)
            self.distributions[v] = value
            self.stats['distributions_' + meta['kind']] += 1
            self.stats['nonuniform_' + meta['kind']] += int(value != 0)

    def estimate(self, p):
        self.initialize_distributions()
        self.phase = 'isolated_variable_substitution'
        e, result = p, Fraction(1)
        while e:
            self.check(e)
            lone = isolated(e)
            rest = e - frozenset(1 << v for v in indices(lone))
            factor = self.bias_p_correlation(rest)
            result *= factor
            item = dict(terms=len(e), variables=support(e).bit_count(),
                        isolated_variables=list(indices(lone)),
                        nonisolated_factor=str(factor), accumulated=str(result))
            self.trace.append(item)
            if not result:
                break
            replacement = ZERO
            for v in indices(lone):
                expr = self.engine.definitions[v]
                if expr is None:
                    # An isolated original uniform variable has correlation 0.
                    result = Fraction(0)
                    item['zero_reason'] = 'isolated uniform original input'
                    item['zero_variable'] = v
                    break
                replacement = add(replacement, expr)
            if not result:
                break
            e = replacement
        self.phase = 'completed'
        return result

    def statistics(self):
        return dict(self.stats, uniform_backend=dict(self.uniform.stats))
