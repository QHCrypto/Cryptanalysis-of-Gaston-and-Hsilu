#!/usr/bin/env python3
"""Seal a completed and independently audited paper-tail experiment snapshot."""
import argparse
import hashlib
import json
from pathlib import Path
from datf_paper_parallel import atomic_json


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True)
    args=p.parse_args();root=args.results
    validation=json.loads((root/'validation.json').read_text());assert validation['status']=='PASS'
    assert validation['every_required_rank_prefix_verified']
    receipt=json.loads((root/'reassessment/completed.json').read_text());assert receipt['all_rows_terminal']
    for marker in (root/'reassessment/jobs').glob('*.running.json'):
        eid=marker.name.removesuffix('.running.json')
        assert (marker.parent/(eid+'.summary.json')).exists()
        marker.replace(marker.with_name(eid+'.execution.json'))
    sources=[Path(name) for name in ['Makefile','README.md','include/gaston.hpp','include/datf_middle.hpp',
             'src/datf_paper_shard.cpp','tools/datf_paper_parallel.py','tools/datf_paper_merge.py',
             'tools/datf_paper_reassess.py','tools/datf_paper_report.py','tools/datf_paper_manifest.py','tools/datf_paper_audit.py',
             'tools/datf.py','tools/datf_k06.py','tools/datf_fixed_tail.py','tools/datf_audit.py',
             'tools/datf_quadratic.py','tools/datf_pool_build.py','tools/fixed_tail_search.py',
             'tests/test_datf_quadratic.py','tests/test_datf_paper_full.py']]
    config=json.loads((root/'configuration.json').read_text())
    sources+=[Path(s['path']) for s in config['sources']]
    sources+=list(Path('results/datf_pools/pools').iterdir())
    paths=sources+[f for f in root.rglob('*') if f.is_file() and f.name!='reproducibility.json']
    records=[{'path':str(f.resolve()),'bytes':f.stat().st_size,'sha256':digest(f)} for f in paths]
    atomic_json(root/'reproducibility.json',{'status':'complete_and_validated','scope':'24-process exhaustive matching and ranked full ordinary-DATF reassessment',
                                          'files':records})
    print('SEALED',len(records),'files',sum(r['bytes'] for r in records),'bytes')


if __name__=='__main__':main()
