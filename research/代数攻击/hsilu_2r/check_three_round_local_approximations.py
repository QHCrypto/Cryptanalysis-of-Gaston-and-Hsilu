"""Test all four 3/4 local chi-bit approximations before the last mixer.

Scope: fixed 21-variable input space, original 42 capacity conditions,
individual first-block output bits. The probability 3/4 refers to a uniform
local 3-bit chi input; it is NOT asserted for the constrained path distribution.
The feasibility test allows a separate approximation choice for every source
of each output bit, so infeasibility also excludes globally shared choices.
"""
import contextlib
import io
import json

with contextlib.redirect_stdout(io.StringIO()):
    import check_three_round_extension as ext

Q = ext.x3_quadratic
# chi output a+(1+b)c has approximants a, a+c, a+b+1, a+b+c.
truth_table_counts = {}
for alpha in (0,1):
    for beta in (0,1):
        delta = alpha & (1 ^ beta)
        matches = 0
        for a in (0,1):
            for b in (0,1):
                for c in (0,1):
                    exact = a ^ ((1 ^ b) & c)
                    approx = a ^ (alpha & b) ^ (beta & c) ^ delta
                    matches += exact == approx
        assert matches == 6
        truth_table_counts[f'alpha={alpha},beta={beta},delta={delta}'] = matches

symbols = [[{64*i+j} for j in range(64)] for i in range(5)]
first_block_sources = ext.linear_poly(symbols)[0]
monomials = sorted(set().union(*(p for row in Q for p in row)))
index = {m:k for k,m in enumerate(monomials)}
def encode(polynomial):
    return sum(1 << index[m] for m in polynomial)

local_options = {}
for i in range(5):
    for z in range(64):
        local_options[i,z] = [(alpha,beta) for alpha in (0,1) for beta in (0,1)
            if not (Q[i][z] ^ (Q[(i+1)%5][z] if alpha else set())
                    ^ (Q[(i+2)%5][z] if beta else set()))]

records = []
for path_bit in range(64):
    sources = sorted(first_block_sources[63-path_bit])
    target = set()
    columns = []
    for node in sources:
        i,z = divmod(node,64)
        target ^= Q[i][z]
        columns += [encode(Q[(i+1)%5][z]), encode(Q[(i+2)%5][z])]
    rank = ext.v.rank(columns)
    augmented_rank = ext.v.rank(columns+[encode(target)])
    records.append({'output_path_bit': path_bit,
                    'source_Y3_positions_path_indexing': [[node//64,63-node%64] for node in sources],
                    'source_count': len(sources),
                    'source_positions_admitting_locally_affine_approximation':
                        [[node//64,63-node%64] for node in sources if local_options[node//64,node%64]],
                    'choice_matrix_rank': rank, 'augmented_rank': augmented_rank,
                    'affine_output_feasible': rank == augmented_rank})

out = {'scope': 'fixed 21 variables, four standard 3/4 local chi-bit approximations, individual first-block output bits',
       'probability_caveat': '6/8 under uniform local input only; no independence or constrained-distribution assumption verified',
       'local_approximation_truth_table_matches_out_of_8': truth_table_counts,
       'affine_X3_positions_by_row': [[j for j in range(64) if not Q[i][63-j]] for i in range(5)],
       'affine_X3_counts_by_row': [sum(not p for p in row) for row in Q],
       'locally_affine_approximation_counts_by_row':
           [sum(bool(local_options[i,z]) for z in range(64)) for i in range(5)],
       'first_block_linear_layer_fanin_values': sorted(set(map(len,first_block_sources))),
       'feasible_individual_first_block_outputs': [r['output_path_bit'] for r in records if r['affine_output_feasible']],
       'output_records': records}
(ext.v.HERE / 'three_round_local_approximation_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:w for k,w in out.items() if k not in ['output_records','affine_X3_positions_by_row']},indent=2))
print('Example output 7:',json.dumps(records[7]))
