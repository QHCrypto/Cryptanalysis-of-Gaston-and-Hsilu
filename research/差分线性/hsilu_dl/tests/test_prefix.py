"""Check exact two-chi prefix sums without using the C++ derivative formula."""
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from hsilu import *

ROOT=Path(__file__).resolve().parents[1]
def dot(a,b):return sum((x&y).bit_count() for x,y in zip(a,b))%2
def chi5(x):return sum((((x>>i)&1)^((1^((x>>((i+1)%5))&1))&((x>>((i+2)%5))&1)))<<i for i in range(5))
BASIS=[linear(unpack(1<<i)) for i in range(320)]
def independent(col,pattern,gamma,rc):
    total=0;cache={}
    for x in range(32):
        y=chi5(x);b=y^chi5(x^pattern)
        if b not in cache:
            delta=linear(embed(col,b))
            def f(state):
                base=(state[0]^rc,*state[1:]);other=tuple(a^d for a,d in zip(base,delta))
                return dot(gamma,tuple(a^b for a,b in zip(chi(base),chi(other))))
            c=f((0,0,0,0,0))
            # Recover the affine function on the previous output state by
            # evaluating its 320 basis vectors, independently of L^T and the
            # optimized derivative coefficient formula.
            coeff=sum((f(v)^c)<<i for i,v in enumerate(BASIS))
            cache[b]=(c,coeff)
        c,coeff=cache[b]
        column_mask=pack(embed(col,31))
        if coeff&~column_mask:continue
        total+=-1 if c^((coeff&pack(embed(col,y))).bit_count()%2) else 1
    return total

g1=tuple(int(w,16) for w in ('000c800010000124','4000040000120001','4010000000000412','1000044000000001','0088000002482010'))
g2=tuple(int(w,16) for w in ('001c802200000124','2000800000120009','4201000000000492','3000c00000000001','1040000002482430'))
cases=[(0,1,embed(1,1),0x4b),(10,1,g1,0x69),(10,16,g1,0x69),(18,8,g2,0x69),
       (31,7,embed(22,13),0x78),(7,31,embed(0,31),0)]
for col,pattern,gamma,rc in cases:
    actual=int(subprocess.check_output([ROOT/'build/hsilu_dl','prefix-eval',str(col),str(pattern),';'.join(words(gamma)),f'{rc:x}'],text=True))
    assert actual==independent(col,pattern,gamma,rc),(col,pattern)
print('PASS exact two-chi correlation against independent affine basis enumeration')
