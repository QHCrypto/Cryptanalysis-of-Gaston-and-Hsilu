"""Independent truth-table checks and comparison with the author's Ascon list."""
import csv
import json
import re
from pathlib import Path
from prepare import HN1,HN2,HM,HR,multiply

BASE=Path(__file__).resolve().parent

def linear_truth(state,n):
    x=list(state)
    for u,v,s in zip(HN1,HN2,HM):
        x.append([x[u][j]^x[v][(j-s)%n] for j in range(n)])
    return [[x[13+i][(j-HR[i])%n] for j in range(n)] for i in range(5)]

def reduced_check(n,rounds):
    # Bit j in these big integers is one complete (nonce,lambda) assignment.
    size=1<<(2*n)
    allones=(1<<size)-1
    variables=[(allones//((1<<(1<<j))+1))<<(1<<j) for j in range(2*n)]
    state=[[0]*n for _ in range(5)]
    for j in range(n):
        state[1][j]=variables[j]&variables[n+j]
        state[3][j]=variables[j]
    state=linear_truth(state,n)
    for _ in range(1,rounds):
        state=linear_truth([[state[(i+1)%5][j]&state[(i+2)%5][j] for j in range(n)] for i in range(5)],n)
    flat=sum(state,[])
    rows=list(csv.DictReader((BASE/f'results_{n}.tsv').open(),delimiter='\t'))
    details=[]
    mask=(1<<(1<<n))-1
    for row in rows:
        truth=allones
        for label in ('p','q','s','t'):truth &= flat[int(row[label])]
        signature=sum((((truth>>(key*(1<<n)))&mask).bit_count()%2)<<key for key in range(1<<n))
        details.append({'id':int(row['id']),'classification':row['classification'],
                        'nonzero':bool(signature),'coefficient_truth_table_hex':hex(signature)})
    false_exclusions=sum(r['nonzero'] and r['classification']=='excluded' for r in details)
    assert false_exclusions==0
    report={'word_bits':n,'prefix_rounds':rounds,'quartics':len(rows),
            'excluded':sum(r['classification']=='excluded' for r in details),
            'retained':sum(r['classification']=='retained' for r in details),
            'nonzero_coefficients':sum(r['nonzero'] for r in details),
            'false_exclusions':false_exclusions,'details':details}
    (BASE/f'validation_{n}.json').write_text(json.dumps(report,indent=2)+'\n')
    return {k:v for k,v in report.items() if k!='details'}

def ascon_tail_check():
    n=64
    x=[[{1<<(n*i+j)} for j in range(n)] for i in range(5)]
    def sub(x):
        a=[[set(p) for p in w] for w in x]
        for j in range(n):
            a[0][j]^=a[4][j];a[4][j]^=a[3][j];a[2][j]^=a[1][j]
        y=[[a[i][j]^multiply(a[(i+1)%5][j]^{0},a[(i+2)%5][j]) for j in range(n)] for i in range(5)]
        for j in range(n):
            y[1][j]^=y[0][j];y[0][j]^=y[4][j];y[3][j]^=y[2][j];y[2][j]^={0}
        return y
    y=sub(x)
    pairs=[(19,28),(61,39),(1,6),(10,17),(7,41)]
    # Author's coordinate convention: output[j] uses input[j-rotation].
    z=[[y[i][j]^y[i][(j-a)%n]^y[i][(j-b)%n] for j in range(n)] for i,(a,b) in enumerate(pairs)]
    out=sub(z)
    formal={tuple(j for j in range(320) if m>>j&1) for m in out[0][0] if m.bit_count()==4}
    source=BASE/'ascon_author_quartics.txt'
    published={tuple(map(int,re.findall(r'x(\d+)',line))) for line in source.read_text().splitlines() if line.startswith('x')}
    assert formal==published
    result={'formal_quartics_before_prefix_screening':len(formal),
            'author_published_quartics':len(published),'sets_identical':True,
            'meaning':'For this chosen Ascon output, the final published list equals all exact formal quartics. No additional net term deletion is visible relative to this already parity-cancelled ANF.',
            'author_result_path':str(source),
            'author_repository':'https://github.com/keyrecoveryofascon/keyrecoveryofascon',
            'author_commit':'616f5d247bd9decdb39a74f2339f73a115392910'}
    (BASE/'ascon_tail_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    for n,r in [(4,1),(8,2)]:print(reduced_check(n,r))
    print(ascon_tail_check())
