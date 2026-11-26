"""Recompute G4/H4 with the original HNW26 column-product readout.

Two-round differential-prefix weighting (HNW26, Eq. 5); columnwise HNW
zero-value-mask propagation in the tail; each final-column mask pulled through L^T
separately. The joint-mask extension is exported only as a comparison.
The last chi is expanded with squared LAT coefficients. Gaussian elimination
is only the summation backend for the affine prefix difference distribution.
This uses model assumptions and is not an exact fixed-permutation certificate.
"""
from hnw_prefix_comparison import compare
from pathlib import Path
import json
import hashlib

HERE = Path(__file__).resolve().parent

records=[]
comparisons=[]
for name in ['Gaston','Hsilu']:
    r=compare(name,prefix='fresh',split_output=True)
    assert r['correlation']['numerator'] == 0
    records.append(r)
    joint=compare(name,prefix='fresh',split_output=False)
    assert joint['correlation']['denominator'] == (256 if name=='Gaston' else 512)
    assert joint['correlation']['numerator'] == 1
    comparisons.append(joint)
    print(name,r['correlation'],flush=True)
doc={'method':'Original HNW26 column-product tail with two-round differential-prefix weighting',
     'reference':'https://eprint.iacr.org/2026/358',
     'parameters':{'rounds':4,'differential_prefix_rounds':2,'value_masks':'zero',
                   'last_chi':'complete nonzero squared-LAT expansion',
                   'linear_readout':'apply L^T separately to each final-column mask; multiply the resulting tail expectations for each fixed prefix difference'},
     'quantity':'signed correlation, not bias',
     'scope':'Exact rational evaluation of the specified approximation, not an exact true permutation correlation',
     'results':records,
     'joint_mask_extension_comparison':{'method':'Joint-mask extension, not original HNW26; compare E[product Z_j] with product E[Z_j]',
                                       'results':comparisons},
     'cpp_column_recurrence_validation':json.loads((HERE/'hnw_fullmask_validation.json').read_text()),
     'source_hashes':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                      for p in ['recompute_for_manuscript.py','hnw_prefix_comparison.py',
                                'joint_estimator.py','core.py','hnw_fullmask_reference.cpp']}}
(HERE/'manuscript_hnw_results.json').write_text(json.dumps(doc,indent=2)+'\n')
