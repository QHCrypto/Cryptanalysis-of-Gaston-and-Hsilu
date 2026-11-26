"""Check original HNW column-product and joint-mask extension separately."""
from joint_estimator import *
import random
import subprocess
from pathlib import Path
import os

os.chdir(Path(__file__).resolve().parent)

rng=random.Random(2026112503)
ACT=[[Fraction(sum(1-2*((b&(sbox(x)^sbox(x^d))).bit_count()&1)
                   for x in range(32)),32) for b in range(32)] for d in range(32)]
def direct(name,d,target,split):
    LT=lgt if name=='Gaston' else lht
    columns=[j for j in range(64) if column(target,j)]
    result=Fraction(0)
    for choices in product(*(SUPPORT[column(target,j)] for j in columns)):
        pulled=[];weight=Fraction(1)
        for j,a in zip(columns,choices):
            pulled.append(LT(embed(j,a)))
            weight*=Fraction(LAT[column(target,j)][a]**2,1024)
        if not split:
            pulled=[add(*pulled)]
        for beta in pulled:
            for j in range(64):weight*=ACT[column(d,j)][column(beta,j)]
        result+=weight
    return result

records=[]
for name in ENDPOINTS:
    exe='hnw_fullmask_reference_'+name.lower()
    command=['g++','-O2','-std=c++17','-Wno-return-type','hnw_fullmask_reference.cpp',
             '-I../gaston_dl/include','-I../hsilu_dl/include','-lgmpxx','-lgmp','-o',exe]
    if name=='Gaston':command.insert(1,'-DTARGET_GASTON')
    subprocess.run(command,check=True)
    delta,target=ENDPOINTS[name]
    L=lg if name=='Gaston' else lh
    j0=next(j for j in range(64) if column(delta,j))
    cases=[(ZERO,target)]
    for _ in range(32):
        x0=rng.randrange(32)
        d2=L(embed(j0,sbox(x0)^sbox(x0^column(delta,j0))))
        y=tuple(rng.getrandbits(64) for _ in range(5))
        d3=L(add(chi(y),chi(add(y,d2))))
        cases.append((d3,target))
    cases += [(embed(rng.randrange(64),rng.randrange(1,32)),target) for _ in range(16)]
    data=''.join(' '.join(f'{v:016x}' for v in d+m)+'\n' for d,m in cases)
    values={}
    for split in [True,False]:
        mode='original_column_product' if split else 'joint_mask_extension'
        lines=subprocess.check_output(['./'+exe]+(['--split'] if split else []),input=data,text=True).splitlines()
        assert len(lines)==len(cases)
        values[mode]=[]
        for (d,m),line in zip(cases,lines):
            r=json.loads(line)
            ref=Fraction(int(r['numerator']),int(r['denominator']))
            assert ref==direct(name,d,m,split),(name,mode,d,m,ref)
            values[mode].append(ref)
        records.append({'target':name,'readout':mode,'cases':len(cases),
                        'nonzero_cases':sum(bool(v) for v in values[mode]),'status':'PASS'})
    differences=sum(a!=b for a,b in zip(*values.values()))
    assert differences > 0, (name,'No test distinguishes the two readouts')
    for row in records[-2:]:row['cases_distinguishing_readouts']=differences
Path('hnw_fullmask_validation.json').write_text(json.dumps(records,indent=2)+'\n')
print(records)
