#!/usr/bin/env python3
"""Fresh, fixed-size H5 sampling; preserve the original discovery batch."""
import argparse
import csv
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statistics(row):
    result = {k: (v if k == 'id' else float(v) if k == 'rho' else int(v))
              for k, v in row.items()}
    n, n0, n1 = result['samples'], result['n0'], result['n1']
    assert n0 + n1 == n
    rho = (n0 - n1) / n
    assert rho == result['rho']
    z = 1.959963984540054
    p = n0 / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    radius = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    threshold = math.sqrt(2 * math.log(2 / .05) / n)
    result.update(log2_abs=math.log2(abs(rho)) if rho else None,
                  ci95=[2*(center-radius)-1, 2*(center+radius)-1],
                  null_standard_error=1/math.sqrt(n),
                  null_z_score=rho*math.sqrt(n),
                  hoeffding_threshold=threshold,
                  significant_hoeffding=abs(rho)>threshold,
                  ci95_excludes_zero=abs(center-.5)>radius)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=1 << 28)
    parser.add_argument('--seed', type=int, default=2026102101)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    assert args.samples > 0 and 1 <= args.threads <= 32
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    baseline = ROOT/'results/baseline'
    candidate = next(x for x in json.loads((baseline/'five_selection.json').read_text())
                     if x['id'] == 'R5_pool_000000')
    assert candidate['rounds'] == 5
    assert candidate['input_difference']['bits'] == [[0,57],[2,57],[3,57]]
    assert candidate['core_output_mask']['bits'] == [[0,10],[2,23]]
    with (baseline/'five_samples.csv').open() as f:
        original_rows = list(csv.DictReader(f))
    assert args.seed not in {int(x['seed']) for x in original_rows}
    original = statistics(next(x for x in original_rows if x['id'] == candidate['id']))
    job = f"{candidate['id']} 5 " + ';'.join(candidate['input_difference']['words']) + ' ' + ';'.join(candidate['core_output_mask']['words']) + '\n'
    jobs = out/'jobs.txt'
    jobs.write_text(job)
    command = [str(ROOT/'build/hsilu_dl'), 'verify', str(jobs), str(out/'samples.csv'),
               str(args.samples), str(args.seed), str(args.threads)]
    configuration = {'date':'2026-10-21','id':candidate['id'],'rounds':5,
        'input_difference':candidate['input_difference'],
        'core_output_mask':candidate['core_output_mask'],
        'full_output_mask':candidate['full_output_mask'],
        'core_boundary':'first chi input to last chi output; first RC absorbed by uniform translation; final L omitted',
        'round_constant_indices':[7,8,9,10,11],
        'round_constants_hex':['87','78','69','5a','4b'],
        'samples':args.samples,'seed':args.seed,'threads':args.threads,
        'fixed_design':{'endpoints':1,'alpha':0.05,'stopping':'fixed sample count; no optional stopping',
                        'primary_test':'two-sided Hoeffding threshold for one prespecified endpoint',
                        'confidence_interval':'Wilson 95% for rho=2p-1',
                        'scope':'fresh pseudorandom simulation; statistical interpretation assumes IID',
                        'original_batch_pooled':False},'command':command}
    (out/'configuration.json').write_text(json.dumps(configuration,ensure_ascii=False,indent=2)+'\n')
    selftest = subprocess.run([str(ROOT/'build/hsilu_dl'),'selftest'],check=True,text=True,capture_output=True)
    (out/'selftest.log').write_text(selftest.stdout+selftest.stderr)
    print(selftest.stdout.strip(),flush=True)
    print('START H5 samples',args.samples,'seed',args.seed,'threads',args.threads,flush=True)
    start = time.perf_counter()
    with (out/'sampling.log').open('w') as log:
        subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
    elapsed = time.perf_counter()-start
    with (out/'samples.csv').open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    result = statistics(rows[0])
    assert result['id'] == candidate['id'] and result['rounds'] == 5
    assert result['samples'] == args.samples and result['seed'] == args.seed
    sources = [Path(__file__).resolve(), ROOT/'build/hsilu_dl', ROOT/'src/main.cpp',
               ROOT/'include/hsilu.hpp', ROOT.parents[1]/'hsilu_ref/hsilu.cpp',
               baseline/'five_selection.json',baseline/'five_samples.csv',jobs,
               out/'configuration.json',out/'samples.csv',out/'selftest.log',out/'sampling.log']
    document = {'status':'completed','configuration':configuration,
                'sampling':result,'previous_sampling':original,
                'elapsed_seconds':elapsed,
                'pairs_per_second':args.samples/elapsed,
                'nonzero_signal_confirmed':result['significant_hoeffding'],
                'interpretation':'finite pseudorandom sampling; not an exact correlation or proof of zero',
                'sources':[{'path':str(p),'sha256':digest(p)} for p in sources]}
    (out/'results.json').write_text(json.dumps(document,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'elapsed_seconds':elapsed,**result},ensure_ascii=False),flush=True)


if __name__ == '__main__':
    main()
