"""Fresh legal NIST-XOF prefixes for the rotated 190-bit Hsilu candidate."""
import ctypes as c
import json
from pathlib import Path
import random
import time
import sys
from nist_modes import XOF_IV, prefix_state

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'MitM'))
import validate_hsilu_3r_candidate as v
from hsilu_structural_audit import rank


def rotated_config():
    cfg=json.loads((ROOT/'MitM/results/hsilu_coloring_3r_pair.json').read_text())
    for field in ('red','blue','capacity1zero','capacity4one','matches'):
        cfg[field]=[(j-2)%64 for j in cfg[field]]
    cfg['cancellations']=[[s,64*(j//64)+(j%64-2)%64] for s,j in cfg['cancellations']]
    cfg['coordinate_shift_mod64']=62
    return cfg


def configure(cfg):
    v.CFG=cfg;v.RED=cfg['red'];v.BLUE=cfg['blue']
    v.C1=[j for s,j in cfg['cancellations'] if s==1]
    v.C2=[j for s,j in cfg['cancellations'] if s==3]


def main():
    cfg=rotated_config();configure(cfg)
    (HERE/'rotated_hsilu_3r_config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    U=c.c_uint64;P=c.POINTER(U)
    lib=c.CDLL(str(HERE/'native.so'));f=lib.prefix_scan
    f.argtypes=[U,c.c_int,U,U,U,U,U,P];f.restype=U
    m1=sum(1<<j for j in cfg['capacity1zero']);m4=sum(1<<j for j in cfg['capacity4one'])
    rng=random.Random(20261109);records=[];total=0;begin=time.time();witness=None
    for batch in range(16):
        first=rng.getrandbits(64);start=rng.getrandbits(64)&~((1<<26)-1);found=(U*7)()
        hits=int(f(XOF_IV,3,first,start,1<<26,m1,m4,found));total+=hits
        rec={'first_word_hex':f'{first:016x}','second_word_start_hex':f'{start:016x}','trials':1<<26,'hits':hits}
        if hits:
            words=list(found[:2]);state=prefix_state(words,'hsilu',3);assert state==list(found[2:])
            base=(state[0]&((1<<62)|(1<<63)))^(1<<62)
            rows,basis=v.c1_matrix(base,state[1:]);q=v.quadratic_labels(base,state[1:],basis)
            rec.update(prefix_hex=[f'{x:016x}' for x in words],state_hex=[f'{x:016x}' for x in state],
                       first_rank=rank(rows),fibre_dimension=len(basis),second_polynomial_rank=q['polynomial_rank_mod_constant'],second_quadratic_rank=q['quadratic_part_rank'])
            if rank(rows)==29 and witness is None:witness=rec
        records.append(rec)
        print(f'batch {batch+1}/16, hits {hits}, cumulative {total}, elapsed {time.time()-begin:.1f}s',flush=True)
    result={'IV_hex':f'{XOF_IV:016x}','trials':1<<30,'hits':total,'uniform_model_expected_hits':64,
            'first_rank_29_witness':witness,'batches':records,'elapsed_seconds':time.time()-begin,
            'scope':'Finite actual-prefix sampling; not a proof of asymptotic prefix density or digest coverage.'}
    (HERE/'nist_hsilu_3r_prefixes.json').write_text(json.dumps(result,indent=2)+'\n')
    assert witness is not None

if __name__=='__main__':main()
