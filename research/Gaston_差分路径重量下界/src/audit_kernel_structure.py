#!/usr/bin/env python3
"""Exact GF(2) audit and count/bit-weight classification; no MILP solve."""
from datetime import datetime
from functools import reduce
import hashlib
import itertools
import json
from operator import xor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MASK = (1 << 64) - 1
EAST = (0, 60, 22, 27, 4)
WEST = (0, 56, 31, 46, 43)
TWIN = (25, 32, 52, 60, 63)


def rol(word, shift):
    shift %= 64
    return ((word << shift) | (word >> ((64-shift) % 64))) & MASK


def maps(theta_input):
    p = reduce(xor, theta_input)
    q = reduce(xor, (rol(w, t) for w, t in zip(theta_input, TWIN)))
    s = p ^ reduce(xor, (rol(q, j) for j in range(18)))
    e = p ^ rol(p, 1) ^ q ^ rol(q, 18)
    return p, q, s, e


def rank(vectors):
    basis = {}
    for value in vectors:
        while value:
            pivot = value.bit_length()-1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def main():
    unit_states = []
    for bit in range(320):
        state = [0]*5
        state[bit//64] = 1 << (bit % 64)
        unit_states.append(state)
    for a in unit_states:
        p, q, s, e = maps(a)
        assert e == s ^ rol(s, 1)
        # Equality of linear functionals, checked on a complete basis.
        assert (e & 0x5555555555555555).bit_count() % 2 == sum(w.bit_count() for w in a) % 2
    ranks = {"E": rank([maps(a)[3] for a in unit_states]),
             "S": rank([maps(a)[2] for a in unit_states]),
             "PQ": rank([maps(a)[0] | (maps(a)[1] << 64) for a in unit_states])}
    assert ranks == {"E": 63, "S": 64, "PQ": 127}
    # Factorization identity on all 64 basis words for P and Q separately.
    for p, q in [(1 << z, 0) for z in range(64)] + [(0, 1 << z) for z in range(64)]:
        s = p ^ reduce(xor, (rol(q, j) for j in range(18)))
        assert s ^ rol(s, 1) == p ^ rol(p, 1) ^ q ^ rol(q, 18)
    assert rank([(1 << z) ^ rol(1 << z, 1) for z in range(64)]) == 63
    examples = []
    for i, j in itertools.combinations(range(5), 2):
        rho = lambda t: 1 ^ 2 ^ (1 << t) ^ (1 << ((t+18) % 64))
        a = [0]*5
        a[i], a[j] = rho(TWIN[j]), rho(TWIN[i])
        assert maps(a)[3] == 0 and sum(w.bit_count() for w in a) == 8
        examples.append(("row_twin", a))
    for i, j, k in itertools.combinations(range(5), 3):
        a = [0]*5
        a[i] = (1 << TWIN[j]) ^ (1 << TWIN[k])
        a[j] = (1 << TWIN[i]) ^ (1 << TWIN[k])
        a[k] = (1 << TWIN[i]) ^ (1 << TWIN[j])
        assert maps(a)[0] == maps(a)[1] == maps(a)[3] == 0
        assert sum(w.bit_count() for w in a) == 6
        examples.append(("vortex", a))
    examples.append(("constant_coset_one", [MASK, 0, 0, 0, 0]))
    assert maps(examples[-1][1])[2] == MASK and maps(examples[-1][1])[3] == 0
    for _, a in examples:
        b = [rol(w, -s) for w, s in zip(a, EAST)]
        kernel_output = [rol(w, s) for w, s in zip(a, WEST)]
        assert kernel_output == [rol(w, e+wshift) for w, e, wshift in zip(b, EAST, WEST)]

    source = ROOT / "references/activity_profiles.json"
    profile = json.loads(source.read_text())
    forward = profile["forward_lower_by_exact_input_count"]
    reverse = profile["reverse_lower_by_exact_output_count"]

    def edge(a, b):
        return b >= forward[a] and a >= reverse[b] and a+b >= 12 and a <= 5*b and b <= 5*a

    branches = {}
    for middle in (6, 7, 8, 9):
        lower = 1000
        possibilities = []
        for c0, c2 in itertools.product(range(1, 65), repeat=2):
            if not edge(c0, middle) or not edge(middle, c2):
                continue
            min_h1 = 2*((max(c0, middle)+1)//2)
            min_h2 = 2*((max(middle, c2)+1)//2)

            def weight_floor(h1, h2):
                return 2*(c0+middle+c2)+(h1-middle+1)//2+(h2-c2+1)//2

            floor = weight_floor(min_h1, min_h2)
            lower = min(lower, floor)
            if floor > 43:
                continue
            for h1 in range(min_h1, 5*min(c0, middle)+1, 2):
                for h2 in range(min_h2, 5*min(middle, c2)+1, 2):
                    w = weight_floor(h1, h2)
                    if w <= 43:
                        possibilities.append({"counts": [c0, middle, c2],
                                              "internal_bit_weights": [h1, h2],
                                              "weight_lower_bound": w})
        branches[str(middle)] = {"count_and_even_bit_weight_lower_bound": lower,
                                 "remaining_types_under_44": possibilities,
                                 "types_are_necessary_not_sufficient": True}
    assert [branches[str(k)]["count_and_even_bit_weight_lower_bound"] for k in (6, 7, 8, 9)] == [36, 42, 43, 45]
    assert len(branches["8"]["remaining_types_under_44"]) == 1
    assert branches["8"]["remaining_types_under_44"][0] == {
        "counts": [6, 8, 7], "internal_bit_weights": [8, 8], "weight_lower_bound": 43}
    assert len(branches["7"]["remaining_types_under_44"]) == 4
    assert not branches["9"]["remaining_types_under_44"]
    report = {
        "status": "PASS", "checked_local": datetime.now().astimezone().isoformat(),
        "scope": "Exact linear algebra and necessary count/bit-weight types; no sparse-kernel catalog enumeration or new MILP solve",
        "coordinates": "A is theta input, after rho_east; h1/h2 are bit weights at a1/a2 and equal those at b0/b1 under kernel propagation",
        "rank": ranks, "E_kernel_dimension": 257, "PQ_kernel_dimension": 193,
        "E_kernel_state_count": "2^257", "two_S_coset_dimensions": [256, 256],
        "kernel_equivalence": "E=0 iff S=P XOR (XOR_{j=0..17} ROL(Q,j)) is all-zero or all-one",
        "kernel_even_bit_weight": True,
        "kernel_combined_row_rotations": [(e+w) % 64 for e, w in zip(EAST, WEST)],
        "linear_basis_vectors_checked": 320,
        "named_examples_checked": {"row_twins": 10, "vortices": 10, "constant_coset_one": 1},
        "named_examples_exhaust_all_low_weight_kernel_states": False,
        "branch_results": branches,
        "profile_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "historical_single_layer_proofs_replayed": False,
    }
    output = ROOT / "results/kernel_structure_audit_20261024.json"
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({"status": "PASS", "ranks": ranks,
                      "weight_floors": {k: v["count_and_even_bit_weight_lower_bound"] for k,v in branches.items()},
                      "remaining_type_counts": {k: len(v["remaining_types_under_44"]) for k,v in branches.items()},
                      "output": str(output)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
