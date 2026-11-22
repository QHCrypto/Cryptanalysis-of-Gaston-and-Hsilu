#!/usr/bin/env python3
"""Exact two-round differential-weight certificate; Python standard library only.

The requested core is A ->chi X ->lambda Y ->chi B, with E(X) != 0.
Weight is -log2(EDP), the sum of local chi DDT weights, not actual DP.
The optional leading lambda completion makes this a full two-round trail
for the official round chi o iota o lambda; it does not change the weight.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
BRANCH = ROOT.parent / "ordinary_branch_proof"
MASK = (1 << 64) - 1
EAST = (0, 60, 22, 27, 4)
WEST = (0, 56, 31, 46, 43)
TWIN = (25, 32, 52, 60, 63)
X = (0, 0x0000008000000100, 0, 0, 1)
Y = (0x0400000000002000, 0x1004000008000020,
     0x0000100002000000, 0x0800010000000000, 0x0100802000000000)
BRANCH_RESULT_SHA256 = "d09b2603ec019fc02a4ed4c91c4da063c259ef309d9877c02978e3ed644a7984"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rol(value, distance):
    distance %= 64
    return ((value << distance) | (value >> ((64 - distance) % 64))) & MASK


def effect(words):
    a = [rol(value, distance) for value, distance in zip(words, EAST)]
    p = q = 0
    for value, distance in zip(a, TWIN):
        p ^= value
        q ^= rol(value, distance)
    return p ^ rol(p, 1) ^ q ^ rol(q, 18)


def linear(words):
    e = rol(effect(words), 23)
    return tuple(rol(rol(value, east) ^ e, west)
                 for value, east, west in zip(words, EAST, WEST))


def linear_by_destinations(words):
    """Separate per-bit implementation, for checking word operations."""
    result = [0] * 5
    for row, value in enumerate(words):
        for bit in range(64):
            if not (value >> bit) & 1:
                continue
            position = (bit + EAST[row]) % 64
            result[row] ^= 1 << ((position + WEST[row]) % 64)
            for shift in (0, 1, TWIN[row], TWIN[row] + 18):
                for output_row in range(5):
                    result[output_row] ^= 1 << (
                        (position + shift + 23 + WEST[output_row]) % 64)
    return tuple(result)


def pack(words):
    return sum(value << (64 * row) for row, value in enumerate(words))


def unpack(value):
    return tuple((value >> (64 * row)) & MASK for row in range(5))


def inverse_image(words):
    """Solve lambda(Z)=words over GF(2); used only to complete the boundary."""
    columns = []
    for index in range(320):
        columns.append(pack(linear_by_destinations(unpack(1 << index))))
    rhs = pack(words)
    rows = [sum(((columns[c] >> r) & 1) << c for c in range(320))
            | (((rhs >> r) & 1) << 320) for r in range(320)]
    for column in range(320):
        pivot = next(r for r in range(column, 320) if rows[r] >> column & 1)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for r in range(320):
            if r != column and rows[r] >> column & 1:
                rows[r] ^= rows[column]
    require(all((row & ((1 << 320) - 1)) == 1 << i
                for i, row in enumerate(rows)), "lambda is not invertible")
    result = unpack(sum((row >> 320) << i for i, row in enumerate(rows)))
    require(linear(result) == words, "inverse boundary completion failed")
    return result


def chi(value):
    return sum((((value >> i) & 1)
                ^ ((1 ^ ((value >> ((i + 1) % 5)) & 1))
                   & ((value >> ((i + 2) % 5)) & 1))) << i for i in range(5))


def column(words, position):
    return sum(((value >> position) & 1) << row
               for row, value in enumerate(words))


def active_columns(words):
    return [j for j in range(64) if column(words, j)]


def describe(words):
    return {
        "words": [f"0x{value:016x}" for value in words],
        "active_columns": active_columns(words),
        "column_weight": len(active_columns(words)),
        "bit_weight": sum(value.bit_count() for value in words),
    }


def identity_chi_transition(words, ddt):
    columns = []
    weight = 0
    for j in range(64):
        value = column(words, j)
        count = ddt[value][value]
        require(count > 0 and count & (count - 1) == 0,
                "invalid or nonintegral-weight chi transition")
        local_weight = 5 - (count.bit_length() - 1)
        weight += local_weight
        if value:
            require(value.bit_count() == 1 and count == 8,
                    "witness active column is not a weight-two identity transition")
            columns.append({
                "column": j, "input_difference": value, "output_difference": value,
                "ddt_count": count, "weight": local_weight,
                "base_values": [v for v in range(32)
                                if chi(v) ^ chi(v ^ value) == value],
            })
    return {"weight": weight, "active_columns": columns}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-branch", action="store_true",
                        help="also check all matching certificates and saved DRAT records")
    parser.add_argument("--full-branch", action="store_true",
                        help="also replay all 1667 previously saved DRAT proofs")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    branch_bytes = (BRANCH / "result.json").read_bytes()
    require(hashlib.sha256(branch_bytes).hexdigest() == BRANCH_RESULT_SHA256,
            "branch-number theorem summary has changed")
    branch = json.loads(branch_bytes)
    require(branch["proof_complete"] is True
            and branch["ordinary_branch_number_outside_theta_kernel"] == 16,
            "missing completed branch-number theorem")
    if args.check_branch or args.full_branch:
        command = [sys.executable, str(BRANCH / "verify_proof.py"),
                   "--jobs", str(args.jobs)]
        if args.full_branch:
            command.append("--full")
        subprocess.run(command, check=True)

    require(len({chi(v) for v in range(32)}) == 32, "chi is not a permutation")
    ddt = [[0] * 32 for _ in range(32)]
    for difference in range(32):
        for value in range(32):
            ddt[difference][chi(value) ^ chi(value ^ difference)] += 1
    require(all(sum(row) == 32 for row in ddt), "wrong DDT row totals")
    require(max(max(row) for row in ddt[1:]) == 8, "wrong differential uniformity")
    require(all(ddt[a][0] == 0 for a in range(1, 32)), "activity not preserved")
    require(all(ddt[0][b] == 0 for b in range(1, 32)), "zero input changed")
    require(all(ddt[1 << i][1 << i] == 8 for i in range(5)),
            "unit difference identity transition failed")

    require(linear(X) == linear_by_destinations(X) == Y, "lambda(X) != Y")
    require(effect(X) == 0x0040000800000000, "middle lambda is not outside kernel")
    require(len(active_columns(X)) == 3 and len(active_columns(Y)) == 13,
            "wrong column counts")
    z = inverse_image(X)
    require(linear_by_destinations(z) == X, "leading lambda check failed")
    require(effect(z) != 0, "leading lambda is not outside kernel")
    first = identity_chi_transition(X, ddt)
    second = identity_chi_transition(Y, ddt)
    require(first["weight"] == 6 and second["weight"] == 26, "wrong trail weight")
    # Since chi is bijective, first-chi input/output have the same activity.
    # Every nonzero DDT entry with nonzero input is <=8/32, hence weight>=2.
    # The certified middle-layer branch bound gives W>=2*(wc(X)+wc(Y))>=32.
    lower = 2 * branch["ordinary_branch_number_outside_theta_kernel"]
    upper = first["weight"] + second["weight"]
    require(lower == upper == 32, "lower and upper bounds do not meet")
    result = {
        "verification": "PASS",
        "round_function": "chi o iota o lambda; constants cancel in differences",
        "requested_core": "A ->chi X ->lambda Y ->chi B",
        "requested_core_witness": "A=X and B=Y",
        "weight_convention": "-log2(EDP), sum of chi DDT restriction weights",
        "minimum_two_round_outside_weight": 32,
        "lower_bound": lower, "upper_bound": upper,
        "kernel_condition": "middle lambda outside E-kernel",
        "witness_also_has_both_round_linear_layers_outside": True,
        "expected_differential_probability": "2^-32",
        "actual_two_round_probability_claimed": False,
        "trail": "Z ->lambda X ->chi X ->lambda Y ->chi Y",
        "Z": describe(z), "X": describe(X), "Y": describe(Y),
        "leading_effect_hex": f"0x{effect(z):016x}",
        "middle_effect_hex": f"0x{effect(X):016x}",
        "chi_1": first, "chi_2": second,
        "ddt_nonzero_input_maximum": 8,
        "branch_theorem_result_sha256": BRANCH_RESULT_SHA256,
        "branch_certificate_checks_this_run": args.check_branch or args.full_branch,
        "branch_drat_replayed_this_run": args.full_branch,
        "lower_bound_argument": "W >= 2*(wc(X)+wc(lambda(X))) >= 2*16 = 32",
    }
    (ROOT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    (ROOT / "chi_ddt.json").write_text(json.dumps(ddt, indent=2) + "\n")
    print("PASS: all 1024 chi input-pair evaluations; nonzero DDT maximum 8/32.")
    print("PASS: Z ->lambda X ->chi X ->lambda Y ->chi Y, both effects nonzero.")
    print(f"Leading effect = {effect(z):016x}; middle effect = {effect(X):016x}.")
    print("PASS: exact two-round minimum weight = 6 + 26 = 32; lower bound = 2*16.")
    if not (args.check_branch or args.full_branch):
        print("The lower bound uses the previously certified branch theorem; "
              "this invocation checked its summary hash, not its proof files.")


if __name__ == "__main__":
    main()
