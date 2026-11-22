"""Recheck full-state first-round affine nonce cubes against reference C.

Run: python3 Cube/src/gaston_first_round_affine_audit.py
The 17-dimensional certificate is an existing workspace result. This audit
does not claim a global maximum among arbitrary linear nonce subspaces.
"""

import ctypes
import hashlib
import json
import random
import subprocess
import tempfile
from pathlib import Path

from linear_cube_structure import BASIS_HEX, bilinear_constraint_rank
from permutations import RC, gaston_linear, permute_window, self_check
from structure import (
    coordinate_graph, maximum_independent_exact, packed, polar, rank,
    state_from_nonce,
)


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT.parent / "gaston_ref/gaston/implementation/reference/gaston.c"
MASK = (1 << 64) - 1


class State(ctypes.Structure):
    _fields_ = [("x", ctypes.c_uint64 * 5)]


def key_direction(bit):
    state = [0] * 5
    state[1 + bit // 64] = 1 << (bit % 64)
    return state


def xor_state(a, b):
    return [x ^ y for x, y in zip(a, b)]


def hex_words(value):
    return [f"{(value >> (64 * i)) & MASK:016x}" for i in range(5)]


def main():
    basis = [int(value, 16) for value in BASIS_HEX]
    nonce_states = [state_from_nonce(1 << j) for j in range(128)]
    nonce_images = [gaston_linear(state) for state in nonce_states]
    key_states = [key_direction(j) for j in range(128)]
    key_images = [gaston_linear(state) for state in key_states]
    basis_states = [state_from_nonce(value) for value in basis]
    basis_images = [gaston_linear(state) for state in basis_states]

    assert rank(basis) == 17
    assert all(polar(a, b) == 0 for i, a in enumerate(basis_images)
               for b in basis_images[:i])
    constraint_rank = bilinear_constraint_rank(basis)
    assert constraint_rank == 111
    graph = coordinate_graph("gaston")
    coordinate_basis, nodes = maximum_independent_exact(graph)
    assert coordinate_basis == [0, 18, 88, 100, 118]

    # Common nullspace of B(L E_N n, L E_K k) for every key direction k.
    mixed_rows = []
    for key_image in key_images:
        columns = [polar(nonce_image, key_image)
                   for nonce_image in nonce_images]
        for bit in range(320):
            row = sum(((value >> bit) & 1) << j
                      for j, value in enumerate(columns))
            if row:
                mixed_rows.append(row)
    mixed_rank = rank(mixed_rows)
    assert mixed_rank == 128

    nonce_witness = polar(nonce_images[0], nonce_images[1])
    key_witness = polar(key_images[0], key_images[1])
    assert nonce_witness and key_witness

    with tempfile.TemporaryDirectory(prefix="gaston_first_round_") as temp:
        library_path = Path(temp) / "reference.so"
        subprocess.run(["cc", "-O2", "-shared", "-fPIC", str(REFERENCE),
                        "-o", str(library_path)], check=True)
        library = ctypes.CDLL(str(library_path))
        library.gaston_round.argtypes = [ctypes.POINTER(State), ctypes.c_uint64]
        library.gaston_round.restype = None

        def reference_round(state, rc):
            value = State((ctypes.c_uint64 * 5)(*state))
            library.gaston_round(ctypes.byref(value), rc)
            return list(value.x)

        def second_difference(base, a, b, rc):
            return (packed(reference_round(base, rc))
                    ^ packed(reference_round(xor_state(base, a), rc))
                    ^ packed(reference_round(xor_state(base, b), rc))
                    ^ packed(reference_round(xor_state(xor_state(base, a), b), rc)))

        zero = [0] * 5
        coordinate_checks = 0
        for i, a in enumerate(nonce_states):
            for j, b in enumerate(nonce_states[:i]):
                assert second_difference(zero, a, b, RC[0]) == polar(
                    nonce_images[i], nonce_images[j])
                coordinate_checks += 1

        rng = random.Random(20261022)
        basis_checks = 0
        random_comparisons = 0
        for start, rc in enumerate(RC):
            base = [rng.getrandbits(64) for _ in range(5)]
            assert reference_round(base, rc) == permute_window(
                base, "gaston", start, 1)
            random_comparisons += 1
            for i, a in enumerate(basis_states):
                for b in basis_states[:i]:
                    assert second_difference(base, a, b, rc) == 0
                    basis_checks += 1
            assert second_difference(base, key_states[0], key_states[1], rc) == key_witness

    report = {
        "date": "2026-10-22",
        "definition": "For every fixed key, all 320 bits after the first initialization round are affine in independent cube coordinates u.",
        "initial_state": ["IV", "K0", "K1", "N0", "N1"],
        "nonce_bit_numbering": "0..63 in N0; 64..127 in N1; bit j has weight 2^j in its word",
        "round": "chi o iota o rho_west o theta o rho_east",
        "criterion": "B(L E_N b_i, L E_N b_j)=0 for all i<j; independent of key, IV, nonce base and round constant",
        "known_answer_checks": self_check(),
        "reference_source": str(REFERENCE),
        "reference_sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
        "reference_checks": {
            "all_coordinate_pair_second_differences": coordinate_checks,
            "basis_pairs_across_all_12_round_constants": basis_checks,
            "random_single_round_comparisons": random_comparisons,
            "key_quadratic_witness_all_12_round_constants": True,
        },
        "coordinate_cube": {
            "maximum_dimension": len(coordinate_basis),
            "example_positions": coordinate_basis,
            "conflict_edges": sum(map(len, graph)) // 2,
            "search_nodes": nodes,
            "proof_method": "exhaustive complement-clique branch and bound",
        },
        "linear_nonce_subspace": {
            "basis_hex": BASIS_HEX,
            "dimension": rank(basis),
            "pair_count": len(basis) * (len(basis) - 1) // 2,
            "all_pair_polars_zero": True,
            "orthogonal_constraint_rank": constraint_rank,
            "orthogonal_complement_dimension": 128 - constraint_rank,
            "maximal_by_inclusion": True,
            "global_maximum_proved": False,
        },
        "joint_key_nonce_affineness": {
            "all_keys_joint_affine": False,
            "K0_bits_0_1_second_difference_words": hex_words(key_witness),
            "mixed_derivative_constraint_rank": mixed_rank,
            "key_independent_nonce_direction_space_dimension": 128 - mixed_rank,
            "interpretation": "Every nonzero nonce direction has a key-dependent full-state derivative; this does not exclude key-independent selected projections.",
        },
        "counterexample": {
            "nonce_bits": [0, 1],
            "second_difference_words": hex_words(nonce_witness),
            "interpretation": "Fixing the key or changing fixed nonce bits cannot remove this full-state quadratic term for these directions.",
        },
        "consequences": {
            "cube_degree_after_r_rounds": "at most min(d, 2^(r-1)) on a first-round affine subspace",
            "five_round_17_dimensional_sum": "identically zero for all keys and all nonce bases by degree",
            "affine_in_key_superpoly": "not implied by first-round affineness in cube variables",
            "full_construction_complexity_claim": False,
        },
    }
    target = ROOT / "results/gaston_first_round_affine_audit_20261022.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "result_path": str(target),
        "reference_checks": report["reference_checks"],
        "coordinate_maximum": len(coordinate_basis),
        "certified_linear_subspace_dimension": rank(basis),
        "orthogonal_constraint_rank": constraint_rank,
        "mixed_key_nonce_constraint_rank": mixed_rank,
        "status": "all exact checks passed",
    }, indent=2))


if __name__ == "__main__":
    main()
