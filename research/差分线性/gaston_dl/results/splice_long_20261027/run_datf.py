import concurrent.futures,json,subprocess,sys
from pathlib import Path
H=Path(__file__).resolve().parent
ROOT=H.parents[2]
jobs=json.loads((H/'manifest.json').read_text())['jobs']
def run(pair):
    i,j=pair;p=H/'raw'/f"{j['id']}_datf.json"
    if not p.exists():
        subprocess.run([sys.executable,str(ROOT/'DATF算法3复核_20261021/worker.py'),
            '--manifest',str(H/'manifest.json'),'--index',str(i),'--output',str(p),
            '--seconds','45'],check=True,stdout=subprocess.DEVNULL,timeout=65)
    x=json.loads(p.read_text())
    if j['rounds']==2 or x['status']!='completed':print(j['id'],'middle',j['rounds'],x['status'],x['correlation'],flush=True)
    return x['status']
if __name__=='__main__':
    from collections import Counter
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        counts=Counter(ex.map(run,enumerate(jobs)))
    print('DATF',dict(counts),flush=True)
