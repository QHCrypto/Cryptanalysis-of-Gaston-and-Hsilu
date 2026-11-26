#!/usr/bin/env python3
"""Validate and export the current seven-row manuscript DL result table.

Four-round results come from fresh two-round-prefix computations. Short pairs
and five-round middle parts are also recomputed with DATF Algorithm 3. Saved
sampling counts are kept separate from exact rational approximation values.
"""
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path
import sys

PACKAGE = Path(__file__).resolve().parents[1]
RESEARCH = PACKAGE/'research'
DL = RESEARCH/'差分线性'
CURRENT = DL/'joint_dl_20261125'
sys.path.insert(0, str(CURRENT))
import datf_two_round_prefix as datf


def load(path):
    return json.loads(path.read_text())


def fraction(value):
    return Fraction(int(value['numerator']), int(value['denominator']))


def record(value):
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator),
            'log2_abs': math.log2(abs(value)) if value else None}


def algorithm3(name, rounds, delta, mask, start):
    delta, mask = (tuple(int(x, 16) for x in words) for words in (delta, mask))
    if name == 'Gaston':
        engine, _, diff, _ = datf.G.construct(delta, rounds, start)
    else:
        engine, _, diff = datf.H.construct(delta, rounds, start=start)
    return datf.Algorithm3(engine, seconds=30).estimate(datf.G.parity_polynomial(diff, mask))


def round_based_middle(name, delta, mask):
    # These archived backends integrate the first-DDT affine image directly,
    # including dense middle masks, without expanding every last-layer LAT term.
    path = (DL/'gaston_dl/results/splice_long_20261027/hnw_middle.py' if name == 'Gaston'
            else DL/'hsilu_dl/results/rounds2to5_20261027/middle2.py')
    spec = importlib.util.spec_from_file_location('current_middle_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.two(tuple(int(x, 16) for x in delta), tuple(int(x, 16) for x in mask))
    return result[0] if name == 'Gaston' else result


def main():
    expected = load(PACKAGE/'EXPECTED_RESULTS.json')['differential_linear']
    short_path = DL/'路径分段汇总_20261027/results.json'
    short = {r['code']: r for r in load(short_path)['rows']}
    gaston_path = DL/'gaston_dl/results/splice_long_20261027/results.json'
    g5 = next(r for r in load(gaston_path)['rows'] if r['id'] == 'D1_M2_D3_chi3_to_4')
    hsilu_path = DL/'hsilu_dl/results/rounds2to5_20261027/results.json'
    h5 = next(r for r in load(hsilu_path)['rows'] if r['code'] == 'H5-new1')
    fresh_datf = {r['target']: r for r in load(CURRENT/'datf_prefix_fullmask_results.json')['results']}
    fresh_hnw = {r['target']: r for r in load(CURRENT/'manuscript_hnw_results.json')['results']}
    output = []
    aliases = {'G2':'G2', 'G3':'G3a', 'G4':'G4a', 'H3':'H3', 'H4':'H4a'}
    for spec in expected:
        ident, name = spec['id'], spec['permutation']
        row = {k: v for k, v in spec.items() if k != 'correlation'}
        if ident in aliases:
            saved = short[aliases[ident]]
            assert (spec['delta'], spec['lambda']) == (saved['delta_core'], saved['lambda_core'])
            row['historical_id'] = aliases[ident]
            sampling = saved['sampling']
            measured = Fraction(sampling['n0']-sampling['n1'], sampling['samples'])
            assert sampling['n0']+sampling['n1'] == sampling['samples']
            assert float(measured) == sampling['rho']
            row['experiment'] = {**sampling, 'exact_sample_correlation': str(measured),
                                 'source': str(short_path.relative_to(PACKAGE)),
                                 'rerun_this_check': False}
            row['data_log2_from_experiment'] = -2*math.log2(abs(measured))
            if ident in ['G4', 'H4']:
                d, h = fresh_datf[name], fresh_hnw[name]
                assert d['delta'] == spec['delta'] and d['lambda'] == spec['lambda']
                assert h['prefix'] == 'fresh' and h['split_final_column_readout'] is False
                assert len(d['first_prefix_branches']) == 4
                assert all(Fraction(b['probability']) == Fraction(1, 4) for b in d['first_prefix_branches'])
                assert d['last_mask_branches'] == 16 and d['all_summed_cases'] == 64
                dc, hc = fraction(d['correlation']), fraction(h['correlation'])
                row['model'] = 'Full two-round differential-prefix weighting; uniform tail base state; full-mask transpose readout'
                row['first_prefix_branches'] = d['first_prefix_branches']
                row['summation_terms'] = d['all_summed_cases']
                row['source'] = str(CURRENT.relative_to(PACKAGE))
            else:
                dc = algorithm3(name, spec['rounds'], spec['delta'], spec['lambda'], 12-spec['rounds'])
                assert dc == fraction(saved['datf_algorithm3']['correlation'])
                hc = fraction(saved['hnw26']['correlation'])
                row['model'] = 'Direct whole-pair evaluation'
                row['round_based_source'] = 'Archived result; its native implementation is tested by the quick profile'
                row['source'] = str(short_path.relative_to(PACKAGE))
        else:
            row['experiment'] = None
            if ident == 'G5':
                candidate = g5['extra']['original_candidate']
                assert (spec['delta'], spec['lambda'], spec['gamma']) == (
                    g5['input_difference'], g5['output_mask'], candidate['middle_output_mask'])
                delta = candidate['middle_input_difference']
                p, q2 = fraction(candidate['differential_probability']), fraction(g5['q_squared'])
                saved_d = fraction(g5['branches'][0]['datf_algorithm3']['correlation'])
                saved_h = fraction(g5['branches'][0]['hnw26']['correlation'])
                source = gaston_path
            else:
                assert (spec['delta'], spec['lambda'], spec['gamma']) == (
                    h5['input_difference'], h5['core_output_mask'], h5['trail']['middle_mask'])
                delta = spec['delta']
                p, q2 = fraction(h5['p']), fraction(h5['q_squared'])
                saved_d, saved_h = fraction(h5['datf_algorithm3']), fraction(h5['hnw26']['zero'])
                source = hsilu_path
            start = 12-spec['rounds']+spec['split'][0]
            middle_d = algorithm3(name, 2, delta, spec['gamma'], start)
            middle_h = round_based_middle(name, delta, spec['gamma'])
            assert (middle_d, middle_h) == (saved_d, saved_h)
            dc, hc = p*middle_d*q2, p*middle_h*q2
            row['factors'] = {'p': str(p), 'middle_datf': str(middle_d),
                              'middle_round_based': str(middle_h), 'q_squared': str(q2)}
            row['data_log2_from_model'] = -2*math.log2(abs(dc))
            row['model'] = 'Selected differential/middle/linear concatenation, not whole-pair evaluation'
            row['source'] = str(source.relative_to(PACKAGE))
        assert dc == hc == Fraction(spec['correlation']), (ident, dc, hc)
        row['datf_algorithm3'] = record(dc)
        row['round_based'] = record(hc)
        output.append(row)
        print(ident, row['evaluation'], 'DATF', dc, 'round-based', hc, flush=True)
    result = {'date': '2026-11-26', 'status': 'PASS', 'quantity': 'signed correlation, not bias',
              'scope': 'Current manuscript Table 7 and Appendix C. Exact rational evaluation of the stated models; sampling counts are historical observations.',
              'legacy': 'Full four-round Algorithm 5 outputs 3/2048 and 0 are retained only by the datf-legacy profile.',
              'rows': output}
    (PACKAGE/'CURRENT_DL_RESULTS.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')


if __name__ == '__main__':
    main()
