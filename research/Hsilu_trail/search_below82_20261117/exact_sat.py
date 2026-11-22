"""Exact full-DDT three-round Hsilu trail-weight feasibility queries."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'pysat_kissat'))
import hsilu_sat as h
from pysat.card import CardEnc, EncType
DDT=h.exact_table('diff')
W=[5-max(row).bit_length()+1 for row in DDT]
V=[min(W[a] for a in range(32) if DDT[a][b]) for b in range(32)]
FULL=h.SBOX_CNF['diff'][:46]
for a in range(32):
    for b in range(32):
        mask=a|(b<<5)
        accepted=all(any(bool(mask>>(abs(v)-1)&1)==(v>0) for v in row) for row in FULL)
        assert accepted==bool(DDT[a][b])
h.SBOX_CNF['diff']=FULL

class Exact(h.HsiluDifferential):
    def bound(self,bound):
        terms=[]
        for r,(name,table) in enumerate([('y',V),('x',W),('x',W)]):
            state=self.vars(name,r)
            act=self.nactive_sboxes_var(r)
            for z in range(64):
                bits=[row[z] for row in state]
                terms.extend([act[z],act[z]])
                for level in range(3,max(table)+1):
                    flag=self.pool.id()
                    for a in range(32):
                        self._add([(-v if a>>i&1 else v) for i,v in enumerate(bits)]+[flag if table[a]>=level else -flag])
                    terms.append(flag)
        self.cnf.extend(CardEnc.atmost(lits=terms,bound=bound,vpool=self.pool,encoding=EncType.totalizer).clauses)
        self.target=bound

    def result(self):
        data=super().result()
        def cols(s):return [sum((w>>z&1)<<i for i,w in enumerate(s)) for z in range(64)]
        bc=cols(self._values('y',0))
        ac=[min((a for a in range(32) if DDT[a][b]),key=lambda a:(W[a],a)) for b in bc]
        a0=[sum((a>>i&1)<<z for z,a in enumerate(ac)) for i in range(5)]
        weights=[sum(V[b] for b in bc)]+[sum(W[a] for a in cols(self._values('x',r))) for r in (1,2)]
        assert sum(weights)<=self.target
        data.update(exact_round_weights=weights,completed_trail_weight=sum(weights),
                    completed_states={'x0':[f'{x:016x}' for x in a0],**data['states'],
                                      'y2':data['terminal_sbox_completion']['output']},
                    weight_bound=self.target,sbox_scope='full 317-pair exact DDT, all three rounds')
        return data

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--bound',type=int,default=81)
    p.add_argument('--seconds',type=int,default=600)
    p.add_argument('--tag',required=True)
    p.add_argument('--active',type=int,nargs=3)
    p.add_argument('--middle-max',type=int)
    p.add_argument('--fix-core',type=Path)
    p.add_argument('--near',type=Path)
    p.add_argument('--distance',type=int,default=8)
    p.add_argument('--near-state',default='y0',choices=['y0','x1','y1','x2'])
    args=p.parse_args()
    model=Exact(3,topology='reference',column_symmetry=not(args.fix_core or args.near),
                output_prefix=ROOT/'results'/args.tag)
    if args.active:
        for r,k in enumerate(args.active):
            if k>=0:model.set_nactive_sboxes_at_round(r,k)
    if args.middle_max:model.set_max_nactive_sboxes_at_round(1,args.middle_max)
    if args.fix_core:
        data=json.loads(args.fix_core.read_text())
        for key in ['y0','x1','y1','x2']:
            model._fix(key[0],int(key[1]),[int(x,16) for x in data['states'][key]])
    if args.near:
        data=json.loads(args.near.read_text())
        key=args.near_state
        state=[int(x,16) for x in data['states'][key]]
        bits=model.vars(key[0],int(key[1]))
        lits=[-bits[i][z] if state[i]>>z&1 else bits[i][z] for i in range(5) for z in range(64)]
        model.cnf.extend(CardEnc.atmost(lits=lits,bound=args.distance,vpool=model.pool,encoding=EncType.totalizer).clauses)
    model.bound(args.bound)
    print('built',model.pool.top,len(model.cnf.clauses),flush=True)
    model.has_trail(timeout=args.seconds)
    print(model.status,flush=True)
    if model.status=='SAT':print(json.dumps(model.result()),flush=True)

if __name__=='__main__':main()
