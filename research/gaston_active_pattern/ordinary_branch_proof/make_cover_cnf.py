#!/usr/bin/env python3
"""Activity-cover relaxation used to exclude wt_bits(E)>=8."""
from pathlib import Path
import sys
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
pool=IDPool()
a,b,e=[[pool.id((name,j)) for j in range(64)] for name in ("a","b","e")]
east,west=(0,60,22,27,4),(0,56,31,46,43)
clauses=[[-e[j],a[(j+23-east[i])%64],b[(j+23+west[i])%64]]
         for j in range(64) for i in range(5)]
clauses.extend(CardEnc.atmost(a+b,15,vpool=pool,encoding=EncType.seqcounter).clauses)
clauses.extend(CardEnc.atmost([-v for v in e],56,vpool=pool,encoding=EncType.seqcounter).clauses)
clauses.append([e[0]])
canonical = "--canonical" in sys.argv
if canonical:
    # Any >=8-bit E has an 8-bit subset with the same activity cover.
    # Its eight cyclic gaps sum to 64. Rotate a gap >=8 immediately before 0.
    clauses.extend(CardEnc.atmost(e,8,vpool=pool,encoding=EncType.seqcounter).clauses)
    clauses.extend([[-e[j]] for j in range(57,64)])
    # A cover of size <=15 can be padded to exactly 15 vertices.
    clauses.extend(CardEnc.atmost([-v for v in a+b],113,vpool=pool,encoding=EncType.seqcounter).clauses)
name = "cover_canonical" if canonical else "cover_e8_total15"
with Path(name+".cnf").open("w") as f:
    f.write(f"p cnf {pool.top} {len(clauses)}\n")
    for c in clauses:f.write(" ".join(map(str,c))+" 0\n")
print(pool.top,len(clauses))
