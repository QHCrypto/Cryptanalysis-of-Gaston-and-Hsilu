"""Original DATF Algorithm 3: eliminate ALL beta, then uniform Bias.

The DAG is constructed by the existing, independently verified permutation code.
No call to OrdinaryDATF.estimate/marginal is made here. Alpha definitions are
deliberately NOT substituted. Values returned are correlations = 2 * paper bias.
Exact quadratic elimination and Shannon decomposition change only the Bias
backend, not the distributional model. A resource limit never returns zero.
"""
import sys
import time
from types import SimpleNamespace
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'gaston_dl/tools'))
from datf import ZERO, ONE, add, indices, support, isolated, condition, components
from datf_quadratic import uniform_quadratic


class ResourceLimit(RuntimeError):
    pass


class Algorithm3:
    def __init__(self, engine, max_terms=1000000, max_products=100000000,
                 max_nodes=200000, seconds=120, progress=None):
        self.engine = engine
        self.max_terms, self.max_products, self.max_nodes = max_terms, max_products, max_nodes
        self.started = time.perf_counter()
        self.seconds = seconds
        self.stats = Counter()
        self.phase = 'beta_expansion'
        self.progress = progress
        self.last_update = 0
        self.expanded_beta = {}
        self.bias_cache = {}
        self.final_polynomial = None
        self.zero_certificate = None
        self.branches = None
        self.beta_mask = sum(1 << i for i,m in enumerate(engine.metadata) if m['kind']=='beta')

    def check(self, terms=0):
        self.stats['peak_polynomial_terms'] = max(self.stats['peak_polynomial_terms'], terms)
        if terms > self.max_terms:
            raise ResourceLimit(f'polynomial terms {terms} exceed {self.max_terms}')
        if self.stats['product_pairs'] > self.max_products:
            raise ResourceLimit(f'monomial product pairs exceed {self.max_products}')
        now = time.perf_counter()
        if now - self.started > self.seconds:
            raise ResourceLimit(f'wall time exceeds {self.seconds} seconds')
        if self.progress and now-self.last_update > 2:
            self.progress(self)
            self.last_update = now

    def multiply(self, a, b):
        if not a or not b: return ZERO
        if a == ONE: return b
        if b == ONE or a == b: return a
        self.stats['product_pairs'] += len(a)*len(b)
        self.check()
        out = set()
        for x in a:
            for y in b:
                t = x | y
                if t in out: out.remove(t)
                else: out.add(t)
            self.check(len(out))
        return frozenset(out)

    def beta(self, v):
        if v not in self.expanded_beta:
            assert self.engine.metadata[v]['kind']=='beta'
            expr = self.engine.definitions[v]
            assert expr is not None
            assert not (support(expr) >> v), 'DAG must be topologically ordered'
            self.expanded_beta[v] = self.expand(expr)
            self.stats['beta_definitions_expanded'] += 1
        return self.expanded_beta[v]

    def expand(self, p):
        self.check(len(p))
        if not (support(p) & self.beta_mask): return p
        out = set()
        for t in sorted(p):
            bs = t & self.beta_mask
            factors = [self.beta(v) for v in indices(bs)]
            factors.sort(key=len)
            prod = frozenset((t & ~self.beta_mask,))
            for f in factors:
                prod = self.multiply(prod, f)
                if not prod: break
            out.symmetric_difference_update(prod)
            self.check(len(out))
        return frozenset(out)

    def uniform_correlation(self, p):
        """Exact expectation for mutually independent uniform formal variables."""
        def calc(f):
            self.check(len(f))
            if f in self.bias_cache:
                self.stats['bias_cache_hits'] += 1
                return self.bias_cache[f]
            self.stats['bias_nodes'] += 1
            if self.stats['bias_nodes'] > self.max_nodes:
                raise ResourceLimit(f'Bias nodes exceed {self.max_nodes}')
            if not f: return Fraction(1)
            if f == ONE: return Fraction(-1)
            if isolated(f):
                self.stats['isolated_uniform_zero_leaves'] += 1
                return Fraction(0)
            if max(t.bit_count() for t in f) <= 2:
                self.stats['quadratic_calls'] += 1
                ans = uniform_quadratic(f)
            else:
                parts = components(f)
                if len(parts)>1:
                    self.stats['separations'] += 1
                    ans=Fraction(1)
                    for part in sorted(parts,key=lambda g:(len(g), support(g).bit_count())):
                        ans *= calc(part)
                        if not ans: break
                else:
                    counts = Counter(v for t in f for v in indices(t))
                    # Exact sum elimination: f=u(v+A)+B, A independent of
                    # u,v. Averaging u imposes v=A, hence E(-1)^f =
                    # (1/2) E(-1)^B[v=A]. No independence approximation.
                    # Restrict to small affine A to avoid polynomial swelling.
                    eliminated = False
                    for u in sorted(counts, key=lambda i:(counts[i],i))[:24]:
                        bit=1<<u
                        derivative=frozenset(t^bit for t in f if t&bit)
                        if len(derivative)>6 or any(t.bit_count()>1 for t in derivative):continue
                        candidates=list(indices(support(derivative)))
                        if not candidates:continue
                        v=min(candidates,key=lambda i:(counts[i],i))
                        A=derivative ^ frozenset((1<<v,))
                        B=frozenset(t for t in f if not t&bit)
                        bv=1<<v
                        B0=frozenset(t for t in B if not t&bv)
                        Bv=frozenset(t^bv for t in B if t&bv)
                        g=add(B0,self.multiply(A,Bv))
                        self.stats['affine_pair_eliminations']+=1
                        ans=calc(g)/2
                        eliminated=True
                        break
                    if eliminated:
                        if len(f)<=2000 and len(self.bias_cache)<20000:self.bias_cache[f]=ans
                        return ans
                    # The choice affects running time only. First-input variables
                    # often cut many products; frequency is a deterministic rule.
                    v = max(counts, key=lambda i:(counts[i], -i))
                    ans=(calc(condition(f,v,0))+calc(condition(f,v,1)))/2
            if len(f)<=2000 and len(self.bias_cache)<20000:
                self.bias_cache[f]=ans
            return ans
        return calc(p)

    def estimate(self, p):
        self.stats['initial_terms']=len(p)
        self.stats['initial_variables']=support(p).bit_count()
        f=self.expand(p)
        assert not (support(f) & self.beta_mask)
        self.final_polynomial=f
        self.stats['expanded_terms']=len(f)
        self.stats['expanded_variables']=support(f).bit_count()
        self.stats['expanded_degree']=max((t.bit_count() for t in f),default=0)
        self.stats['remaining_beta_variables']=0
        self.stats['expanded_input_variables']=sum(self.engine.metadata[v]['kind']=='input' for v in indices(support(f)))
        self.stats['expanded_alpha_variables']=sum(self.engine.metadata[v]['kind']=='alpha' for v in indices(support(f)))
        self.stats['max_component_variables']=max((support(c).bit_count() for c in components(f)),default=0)
        if isolated(f):
            v=next(indices(isolated(f)))
            self.zero_certificate={'kind':'isolated_uniform_variable_after_ALL_beta_expansion','variable':v,'metadata':self.engine.metadata[v]}
        self.phase='uniform_Bias'
        self.check()
        ans=self.uniform_correlation(f)
        self.phase='completed'
        return ans

    def estimate_input_shannon(self, p, max_inputs=12):
        """The SAME Algorithm-3 sum, with exact input Shannon splitting first.

        If F=ExpandBeta(e), substitution commutes with expansion:
        E[F] = 2^-k sum_a E[ExpandBeta(e|inputs=a)]. Alpha variables remain
        free uniform formal variables in every branch. This is not conditional
        sampling of actual permutation states and adds no model assumptions.
        """
        seen=set()
        free=set()
        def collect(q):
            for v in indices(support(q)):
                if self.engine.metadata[v]['kind']=='beta':
                    if v not in seen:
                        seen.add(v);collect(self.engine.definitions[v])
                else:free.add(v)
        collect(p)
        inputs=sorted(v for v in free if self.engine.metadata[v]['kind']=='input')
        if len(inputs)>max_inputs:raise ResourceLimit(f'input Shannon variables {len(inputs)} exceed {max_inputs}')
        self.phase='input_Shannon_then_beta_expansion'
        self.branches=[]
        self.stats['input_shannon_variables']=len(inputs)
        self.stats['input_shannon_total_branches']=1<<len(inputs)
        total=Fraction(0)
        mask=sum(1<<v for v in inputs)
        def subst(q,assignment):
            if q is None:return None
            bad=mask&~assignment
            result=set()
            for t in q:
                if t&bad:continue
                term=t&~mask
                if term in result:result.remove(term)
                else:result.add(term)
            return frozenset(result)
        for b in range(1<<len(inputs)):
            self.check()
            bits=sum(((b>>i)&1)<<v for i,v in enumerate(inputs))
            e=SimpleNamespace(metadata=self.engine.metadata,
                              definitions=[subst(q,bits) for q in self.engine.definitions])
            a=Algorithm3(e,self.max_terms,self.max_products-self.stats['product_pairs'],
                         self.max_nodes-self.stats['bias_nodes'],self.seconds-(time.perf_counter()-self.started))
            try:
                value=a.estimate(subst(p,bits))
                self.branches.append({'assignment':b,'input_variables':inputs,
                                      'numerator':str(value.numerator),'denominator':str(value.denominator),
                                      'statistics':dict(a.stats)})
                total+=value
                self.stats['input_shannon_completed_branches']+=1
            finally:
                for key,val in a.stats.items():
                    if key.startswith(('peak_','expanded_','max_')):
                        self.stats['max_branch_'+key]=max(self.stats['max_branch_'+key],val)
                    elif key in ('bias_nodes','product_pairs','beta_definitions_expanded','quadratic_calls','affine_pair_eliminations','isolated_uniform_zero_leaves'):
                        self.stats[key]+=val
                if a.phase!='completed':self.phase='branch_'+str(b)+'_'+a.phase
        self.phase='completed'
        return total/(1<<len(inputs))
