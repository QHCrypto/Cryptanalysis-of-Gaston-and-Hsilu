#!/usr/bin/env python3
"""Check the SAT circuit/activity prefix against independent bit destinations."""
import json
from pathlib import Path
import random

ROOT=Path(__file__).resolve().parent
N=(1<<64)-1
EAST,WEST,TWIN=(0,60,22,27,4),(0,56,31,46,43),(25,32,52,60,63)


def main():
    mapping=json.loads((ROOT/"total15.map.json").read_text())
    clauses=[]
    for line in (ROOT/"total15.cnf").read_text().splitlines():
        if not line or line[0] in "cp":continue
        values=list(map(int,line.split()))
        if values[-1]!=0:raise ValueError("unexpected multiline clause")
        clauses.append(values[:-1])
    circuit=clauses[:5184]  # 64 columns × (64 XOR + 12 activity + 5 cover clauses)
    if max(abs(v) for c in circuit for v in c)!=1152:raise ValueError("circuit prefix changed")
    rng=random.Random(20261020)
    inputs=[[0]*5]
    for row in range(5):
        for j in range(64):
            x=[0]*5;x[row]=1<<j;inputs.append(x)
    inputs += [[rng.getrandbits(64) for _ in range(5)] for _ in range(128)]
    for x in inputs:
        y=[0]*5
        p=q=e=0
        for r,w in enumerate(x):
            for j in range(64):
                if not(w>>j&1):continue
                z=(j+EAST[r])%64
                p^=1<<z;q^=1<<((z+TWIN[r])%64)
                y[r]^=1<<((z+WEST[r])%64)
                for d in (0,1,TWIN[r],TWIN[r]+18):
                    e^=1<<((z+d)%64)
                    for rr in range(5):
                        y[rr]^=1<<((z+d+23+WEST[rr])%64)
        assignment={}
        for words,key in ((x,"x"),(y,"y")):
            for row,variables in enumerate(mapping[key]):
                for j,v in enumerate(variables):assignment[v]=(words[row]>>j)&1
        for value,key in ((p,"P"),(q,"Q"),(e,"E")):
            for j,v in enumerate(mapping[key]):assignment[v]=(value>>j)&1
        for words,key in ((x,"active_input"),(y,"active_output")):
            for j,v in enumerate(mapping[key]):
                assignment[v]=int(any(w>>j&1 for w in words))
        aux=961
        for j in range(64):
            sets=(
                [mapping["x"][i][(j-EAST[i])%64] for i in range(3)],
                [mapping["x"][i][(j-EAST[i]-TWIN[i])%64] for i in range(3)],
                [mapping["P"][j],mapping["P"][(j-1)%64],mapping["Q"][j]],
            )
            for vs in sets:
                assignment[aux]=sum(assignment[v] for v in vs)%2;aux+=1
        for index,clause in enumerate(circuit):
            if not any(assignment[abs(v)]==(v>0) for v in clause):
                raise ValueError(f"circuit clause {index} disagrees with direct propagation")
    result={"verification":"PASS","unit_inputs_checked":320,"zero_input_checked":True,
            "random_inputs_checked":128,"circuit_activity_and_cover_clauses_checked_per_input":5184}
    (ROOT/"encoding_verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
