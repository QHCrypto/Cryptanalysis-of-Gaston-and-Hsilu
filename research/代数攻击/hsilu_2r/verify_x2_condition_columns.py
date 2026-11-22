"""Trace chi columns 7 and 17, without imposing X2[1][17]=0."""
import ctypes
import itertools
import json
import random
import verify_path_algebra as v

path_u = v.DATA['free_variables_by_lane'][0]
v.U = [63-j for j in path_u]
v.G = [63-j for j in v.DATA['guessed_by_lane'][0]]
v.S = [63-j for j in v.DATA['selected_outputs_by_lane'][0]]
v.UMASK = sum(1 << j for j in v.U)
columns = [v.linear([1 << j, 0, 0, 0, 0]) for j in v.U]
traces = {}
for j in (7, 17):
    active = [[u for u, col in zip(path_u, columns) if col[i] >> (63-j) & 1] for i in range(5)]
    quadratic = []
    for i in range(5):
        pairs = []
        for k, l in itertools.combinations(range(21), 2):
            a, b = columns[k], columns[l]
            q = (a[(i+1)%5] & b[(i+2)%5]) ^ (b[(i+1)%5] & a[(i+2)%5])
            if q >> (63-j) & 1:
                pairs.append([path_u[k], path_u[l]])
        quadratic.append(pairs)
    assert not any(quadratic)
    traces[str(j)] = {'X2_active_message_positions_by_row': active,
                      'Y2_quadratic_active_pairs_by_row': quadratic}

random.seed(20261013)
lib = ctypes.CDLL(str(v.HERE / 'verify_reference_wrapper.so'))
lib.verify_round.argtypes = [ctypes.POINTER(ctypes.c_uint64), ctypes.c_uint64]
for trial in range(1000):
    state = [random.getrandbits(64) for _ in range(5)]
    state[0] &= ~v.UMASK
    state[1] &= ~v.UMASK
    state[4] |= v.UMASK
    if not (v.rnd(state, 0x5a)[1] >> (63-17) & 1):
        state[0] ^= 1 << (63-24)
    base, cols, qmask, rows = v.affine_data(state)
    if v.rank(rows) != 21:
        continue
    assert not any(qmask >> j & 1 for j in v.S)
    native = (ctypes.c_uint64 * 5)(*state)
    lib.verify_round(native, 0x5a)
    assert list(native) == v.rnd(state, 0x5a)
    assert native[1] >> (63-17) & 1 == 1
    mids = []
    for u in [0, (1<<21)-1, 0x123456]:
        candidate = state.copy()
        candidate[0] ^= sum(((u >> k) & 1) << j for k,j in enumerate(v.U))
        mid = v.rnd(candidate, 0x5a)
        assert mid[1] >> (63-17) & 1 == 1
        mids.append({'active_assignment': u, 'X2_1_17': 1})
    witness = {'native_base_input_words_hex': [f'{w:016x}' for w in state],
               'X2_1_17': 1, 'selected_equation_rank': 21,
               'all_affine_output_path_positions': [j for j in range(64) if not qmask >> (63-j) & 1],
               'active_assignment_checks': mids}
    break
else:
    raise AssertionError('No full-rank witness found')

result = {'path_to_native_bit': '63-j', 'definition': 'Y2 = chi(iota(X2)), before MRM and rho',
          'column_traces': traces, 'condition_equal_one_witness': witness}
(v.HERE / 'x2_condition_column_validation.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
