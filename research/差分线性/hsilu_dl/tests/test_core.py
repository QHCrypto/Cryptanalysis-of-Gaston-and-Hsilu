import random
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from hsilu import *
from datf import evaluate

rng=random.Random(20261018)
inv=inverse_function(linear);invt=inverse_function(transpose)
for i in range(320):
    a=unpack(1<<i)
    assert inv(linear(a))==a and invt(transpose(a))==a
for R,start in ((2,10),(3,9),(4,8),(3,7)):
    delta=embed(0,3)
    e,a,d=construct(delta,R,start=start)
    for _ in range(4):
        s=tuple(rng.getrandbits(64) for i in range(5))
        assignment=e.evaluate_dag([(s[i]>>j)&1 for i in range(5) for j in range(64)])
        values=lambda t:tuple(sum(evaluate(t[i][j],assignment)<<j for j in range(64)) for i in range(5))
        base=core(s,R,start);pair=core(tuple(u^v for u,v in zip(s,delta)),R,start)
        assert values(a)==base and values(d)==tuple(u^v for u,v in zip(base,pair))
print('PASS 320-basis inverse/inverse-transpose; 2/3/4-chi DAGs; five-round-prefix RC indices')
