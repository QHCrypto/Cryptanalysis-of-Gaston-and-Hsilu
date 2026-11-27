"""Recompute G4/H4 with HNW's joint final-layer treatment.

Two-round differential-prefix weighting (HNW26, Eq. 5) is combined with the
joint final-layer treatment in the authors' public PRESENT implementation.
The basic Algorithm 1 column-product readout is exported as a comparison.
The last chi is expanded with squared LAT coefficients. Gaussian elimination
is only the summation backend for the affine prefix difference distribution.
This uses model assumptions and is not an exact fixed-permutation certificate.
"""
from hnw_prefix_comparison import compare, HNW_IMPLEMENTATION
from pathlib import Path
import json
import hashlib

HERE = Path(__file__).resolve().parent

records=[]
comparisons=[]
for name in ['Gaston','Hsilu']:
    r=compare(name,prefix='fresh',split_output=False)
    assert r['correlation']['denominator'] == (256 if name=='Gaston' else 512)
    assert r['correlation']['numerator'] == 1
    records.append(r)
    basic=compare(name,prefix='fresh',split_output=True)
    assert basic['correlation']['numerator'] == 0
    comparisons.append(basic)
    print(name,r['correlation'],flush=True)
doc={'method':'HNW26 two-round differential-prefix weighting with the authors\' joint final-layer treatment',
     'reference':'https://eprint.iacr.org/2026/358',
     'implementation_source':HNW_IMPLEMENTATION,
     'parameters':{'rounds':4,'differential_prefix_rounds':2,'value_masks':'zero',
                   'last_chi':'complete nonzero squared-LAT expansion',
                   'linear_readout':'combine all active final-column masks and apply L^T to the complete mask before taking the tail expectation for each fixed prefix difference'},
     'quantity':'signed correlation, not bias',
     'scope':'Exact rational evaluation of the specified approximation, not an exact true permutation correlation',
     'results':records,
     'basic_column_product_comparison':{'method':'Basic HNW26 Algorithm 1 column-product readout with the same two-round prefix; product E[Z_j] instead of E[product Z_j]',
                                        'results':comparisons},
     'cpp_column_recurrence_validation':json.loads((HERE/'hnw_fullmask_validation.json').read_text()),
     'source_hashes':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                      for p in ['recompute_for_manuscript.py','hnw_prefix_comparison.py',
                                'joint_estimator.py','core.py','hnw_fullmask_reference.cpp']}}
(HERE/'manuscript_hnw_results.json').write_text(json.dumps(doc,indent=2)+'\n')
