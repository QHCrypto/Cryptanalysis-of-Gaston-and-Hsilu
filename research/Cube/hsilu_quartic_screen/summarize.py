"""Merge completed batches, retain unresolved rows, and count filtered tails."""
import csv
import json
from collections import Counter
from pathlib import Path

BASE=Path(__file__).resolve().parent

def summarize():
    data=json.loads((BASE/'terms_64.json').read_text())
    rows={}
    for path in sorted(BASE.glob('results_64*.tsv')):
        for row in csv.DictReader(path.open(),delimiter='\t'):
            index=int(row['id'])
            assert tuple(int(row[k]) for k in ('p','q','s','t'))==tuple(data['classes'][index])
            if len(row)!=10 or row['classification'] not in ('retained','excluded','unresolved'):
                raise ValueError(f'invalid row in {path}')
            if index in rows:
                old=rows[index]
                if old['classification']!='unresolved' and row['classification']!='unresolved':
                    assert old['classification']==row['classification']
                if row['classification']=='unresolved':continue
            rows[index]=row
    byterm={tuple(data['classes'][i]):row['classification'] for i,row in rows.items()}
    mapping={}
    for name,entry in data['mapping'].items():
        counts=Counter(byterm.get(tuple(t),'not_solved') for t in entry['classes'])
        mapping[name]={'formal_quartic_count':entry['count'],
                       'retained_term_count':counts['retained'],
                       'excluded_term_count':counts['excluded'],
                       'unresolved_term_count':counts['unresolved']+counts['not_solved']}
    counts=Counter(row['classification'] for row in rows.values())
    complete=len(rows)==len(data['classes']) and counts['unresolved']==0
    report={
        'complete':complete,
        'model':'Hsilu (IV,K0,K1,N,0), seven rounds, N is a full 64-bit cube; first 128 output bits',
        'screen':'First five rounds, four selected state-coordinate factors, all 64 nonce exponents fixed to one; lambda=1+K1; maximal nonce-degree reduction; COPY/AND/XOR monomial trail rules with valid saturated-degree disjointness cuts',
        'scope':'Trail infeasibility screening, NOT parity enumeration or proof of nonzero coefficients',
        'rotation_classes_total':len(data['classes']),
        'rotation_classes_solved':len(rows),
        'rotation_classes_retained':counts['retained'],
        'rotation_classes_excluded':counts['excluded'],
        'rotation_classes_unresolved':counts['unresolved']+len(data['classes'])-len(rows),
        'tail_counts':mapping,
        'selected_solve_seconds':sum(float(row['seconds']) for row in rows.values()),
        'solves':[dict(row) for _,row in sorted(rows.items())],
        'validation':[json.loads((BASE/f'validation_{n}.json').read_text()) for n in (4,8)],
        'ascon_tail_comparison':json.loads((BASE/'ascon_tail_comparison.json').read_text()),
        'interpretation':('All rotation classes have a feasible trail. This screening proves no additional zero terms; it does not prove all retained superpolies nonzero.'
                          if complete and counts['excluded']==0 else 'Incomplete or filtered result; consult individual statuses.'),
        'effect_on_previous_costs':('No reduction from this screening alone; previous budgets remain unscreened upper-bound schedules, not optimal complexities.'
                                   if complete and counts['excluded']==0 else 'Do not infer a final complexity reduction from a partial screen.')
    }
    destination=BASE.parent/'results'/'hsilu_7r_quartic_screen.json'
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('solves','validation','screen')},indent=2))

if __name__=='__main__':summarize()
