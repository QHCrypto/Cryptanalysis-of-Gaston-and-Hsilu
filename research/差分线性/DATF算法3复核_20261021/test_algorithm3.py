import random
import unittest
from fractions import Fraction
from algorithm3 import Algorithm3, ResourceLimit, ONE, ZERO, indices, support
from datf import OrdinaryDATF, add, mul, variable, evaluate


def brute(p):
    vs=list(indices(support(p)))
    return Fraction(sum(1-2*evaluate(p,sum(((x>>i)&1)<<v for i,v in enumerate(vs))) for x in range(1<<len(vs))),1<<len(vs))


class TestAlgorithm3(unittest.TestCase):
    def test_uniform_bias_against_enumeration(self):
        rng=random.Random(20261021)
        e=OrdinaryDATF()
        for i in range(10):e.fresh()
        for n in range(1,11):
            for trial in range(35):
                p=frozenset(rng.randrange(1<<n) for _ in range(35))
                a=Algorithm3(e)
                self.assertEqual(a.uniform_correlation(p),brute(p))

    def test_beta_expansion_against_independent_dag_evaluation(self):
        rng=random.Random(116)
        for trial in range(40):
            e=OrdinaryDATF()
            for i in range(4):e.fresh()
            # Alpha has a definition but is a FREE variable in Algorithm 3.
            e.fresh(add(variable(0),variable(1)),kind='alpha')
            for b in range(4):
                n=len(e.definitions)
                p=frozenset(rng.randrange(1<<n) for _ in range(8))
                e.fresh(p,kind='beta')
            p=frozenset(rng.randrange(1<<len(e.definitions)) for _ in range(8))
            a=Algorithm3(e)
            f=a.expand(p)
            self.assertEqual(support(f)&a.beta_mask,0)
            total=0
            for x in range(32):
                values=x
                for i in range(5,len(e.definitions)):
                    values|=evaluate(e.definitions[i],values)<<i
                y=evaluate(p,values)
                self.assertEqual(evaluate(f,x),y)
                total+=1-2*y
            self.assertEqual(a.estimate(p),Fraction(total,32))
            self.assertEqual(Algorithm3(e).estimate_input_shannon(p),Fraction(total,32))

    def test_shared_beta_and_unexpanded_alpha(self):
        e=OrdinaryDATF()
        x=e.fresh();y=e.fresh()
        a=e.fresh(add(x,y),kind='alpha')
        b=e.fresh(add(x,y),kind='beta')
        # Complete beta substitution cancels two dependent terms exactly.
        self.assertEqual(Algorithm3(e).estimate(add(mul(b,x),mul(b,y),b)),1)
        # Alpha must not be expanded, even though its definition is x+y.
        self.assertEqual(Algorithm3(e).estimate(add(a,x,y)),0)

    def test_budget_is_not_zero(self):
        e=OrdinaryDATF();x=e.fresh();y=e.fresh()
        b=e.fresh(add(x,y),kind='beta')
        with self.assertRaises(ResourceLimit):Algorithm3(e,max_terms=1).estimate(b)


if __name__=='__main__':unittest.main()
