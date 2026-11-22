#!/usr/bin/env python3
"""Independent certificates and orbit-coverage checks for E != 0 results.

The lower bounds for K <= 5 come from complete unrestricted enumerations,
whose optimal witnesses are checked here to satisfy the additional condition.
For K=6, the per-position lower-bound records use disjoint zero-pair pruning.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import verify as reference

ROOT = Path(__file__).resolve().parent
MASK = (1 << 64) - 1


def rol(x, n):
    n %= 64
    return ((x << n) | (x >> ((64 - n) % 64))) & MASK


def theta_effect(x):
    # Independent bit-destination formula, before the common rotation by 23.
    east, twin = (0, 60, 22, 27, 4), (25, 32, 52, 60, 63)
    e = 0
    for row in range(5):
        for bit in range(64):
            if x >> (row * 64 + bit) & 1:
                for shift in (0, 1, twin[row], twin[row] + 18):
                    e ^= 1 << ((bit + east[row] + shift) % 64)
    return e


def state(data, side):
    return sum(int(w, 16) << (64 * r) for r, w in enumerate(data[side + "_words"]))


def orbit_count(k):
    return sum(math.comb(math.gcd(64, s), k // (64 // math.gcd(64, s)))
               for s in range(64) if k % (64 // math.gcd(64, s)) == 0) // 64


def audit_zero_pairs(inverse, count=26):
    def rank(rows):
        basis = {}
        for value in rows:
            while value:
                bit = value.bit_length() - 1
                if bit in basis:
                    value ^= basis[bit]
                else:
                    basis[bit] = value
                    break
        return len(basis)

    for gap in range(1, 64):
        blocks = [[sum(((inverse[64 * row + col] >> (64 * r + c)) & 1) << (5 * s + r)
                       for s, c in enumerate((0, gap)) for r in range(5))
                   for row in range(5)] for col in range(64)]
        used = set()
        for a in range(64):
            if a in used:
                continue
            for b in range(a + 1, 64):
                if b not in used and rank(blocks[a] + blocks[b]) == 10:
                    used.update((a, b))
                    break
            if len(used) == 2 * count:
                break
        assert len(used) == 2 * count, (gap, len(used))
    print(f"PASS: independent inverse matrix supplies {count} disjoint rank-10 column pairs "
          "for every nonzero output-column gap.", flush=True)


def check_pair(data, forward, k):
    x, y = state(data, "input"), state(data, "output")
    assert reference.apply(forward, x) == y
    assert reference.weight(x) == data["objective"] == data["input_active"]
    assert reference.weight(y) == k == data["output_active"]
    e = theta_effect(x)
    assert e != 0
    # This identity is used to reject E=0 cheaply in the GPU.
    assert (x & MASK) ^ (y & MASK) == rol(e, 23)
    return x, y, e


def validate_six(prefix, forward, inverse):
    path = Path(str(prefix) + ".json")
    data = json.loads(path.read_text())
    assert data["method"] == "exact_enumeration_with_zero_pair_pruning"
    assert data["mode"] == "difference" and data["fixed_side"] == "output"
    assert data["fixed_count"] == 6 and data["theta_effect_constraint"] == "nonzero"
    total = orbit_count(6)
    assert data["position_sets"] == total
    assert 1 <= data["cutoff"] <= 26
    images = [reference.apply(inverse, 1 << j) for j in range(320)]
    seen = set()
    best, tested, witnesses = None, 0, 0
    with Path(str(prefix) + ".positions.tsv").open() as f:
        for offset, row in enumerate(csv.DictReader(f, delimiter="\t")):
            index = int(row["position_index"])
            mask = int(row["position_mask_hex"], 16)
            threshold = int(row["search_below"])
            assert index == data["start_position"] + offset < total
            assert mask.bit_count() == 6 and mask not in seen
            seen.add(mask)
            cols = [j for j in range(64) if mask >> j & 1]
            assert cols[0] == 0 and all(mask <= rol(mask, -j) for j in cols)
            assert 1 <= threshold <= data["cutoff"]
            if data["objective"] is not None:
                assert threshold >= data["objective"]
            evaluated = int(row["tested"])
            assert 0 <= evaluated <= 31 ** 6
            tested += evaluated
            if row["minimum_below_cutoff"] == "null":
                assert row["base31_value_code"] == "null"
                continue
            w, code = int(row["minimum_below_cutoff"]), int(row["base31_value_code"])
            assert 0 < w < threshold and 0 <= code < 31 ** 6
            x = y = 0
            for col in cols:
                value = code % 31 + 1
                code //= 31
                for r in range(5):
                    if value >> r & 1:
                        bit = 64 * r + col
                        y |= 1 << bit
                        x ^= images[bit]
            assert code == 0 and reference.weight(x) == w
            assert theta_effect(x) != 0 and reference.apply(forward, x) == y
            best = w if best is None else min(best, w)
            witnesses += 1
    assert len(seen) == data["completed_position_sets"]
    assert tested == data["tested_candidates_with_duplicates"]
    assert len(seen) * 31 ** 6 == data["covered_output_states"]
    assert data["enumeration_complete"] == (data["start_position"] == 0 and len(seen) == total)
    assert data["proven_optimal"] == (data["enumeration_complete"] and best is not None)
    assert data["objective"] == best
    if best is not None:
        _, _, e = check_pair(data, forward, 6)
        assert int(data["theta_effect_hex"], 16) == e
    print(f"PASS K=6: {len(seen)}/{total} distinct canonical position orbits; "
          f"{witnesses} witnesses; best={best}; proven={data['proven_optimal']}", flush=True)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--six", type=Path, default=ROOT / "result/diff_out6_nonzero")
    parser.add_argument("--summary", type=Path, default=ROOT / "result/nonzero_output_1_to_6.json")
    parser.add_argument("--partial-only", action="store_true",
                        help="Validate only the supplied six-column prefix, including partial searches")
    args = parser.parse_args()
    forward, _ = reference.forward_matrix()
    inverse = reference.inverse(forward)
    audit_zero_pairs(inverse)
    if args.partial_only:
        validate_six(args.six, forward, inverse)
        return
    rows = []
    sources = ["diff_out1_nonzero_witness", "diff_out2", "diff_out3", "diff_out4_gpu", "diff_out5_gpu"]
    for k, name in enumerate(sources, 1):
        path = ROOT / "result" / (name + ".json")
        data = json.loads(path.read_text())
        assert data["mode"] == "difference" and data["fixed_side"] == "output" and data["fixed_count"] == k
        assert data["enumeration_complete"] and data["proven_optimal"]
        if "position_sets" in data:
            expected = orbit_count(k) * 31 ** k
            assert data["completed_position_sets"] == data["position_sets"] == orbit_count(k)
        else:
            expected = math.comb(63, k - 1) * 31 ** k
        assert data["visited"] == expected
        _, _, e = check_pair(data, forward, k)
        rows.append({"output_active": k, "minimum_input_active": data["objective"],
                     "proven_optimal": True, "theta_effect_hex": f"0x{e:016x}",
                     "proof": "complete_unrestricted_enumeration_and_optimal_witness_with_nonzero_effect",
                     "source": str(path), "enumeration_seconds": data["search_seconds"],
                     "input_words": data["input_words"], "output_words": data["output_words"]})
        print(f"PASS K={k}: unrestricted optimum {data['objective']} attained with E=0x{e:016x}", flush=True)
    six = validate_six(args.six, forward, inverse)
    assert six["proven_optimal"], "The six-column search must finish before publishing the complete table"
    rows.append({"output_active": 6, "minimum_input_active": six["objective"],
                 "proven_optimal": True, "theta_effect_hex": six["theta_effect_hex"],
                 "proof": "complete_position_orbits_with_exact_disjoint_zero_pair_pruning",
                 "source": str(args.six) + ".json", "enumeration_seconds": six["search_seconds"],
                 "input_words": six["input_words"], "output_words": six["output_words"]})
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps({"mode": "difference", "condition": "E(rho_east(X)) != 0",
                                      "metric": "nonzero 5-bit S-box columns", "results": rows},
                                     ensure_ascii=False, indent=2) + "\n")
    print(f"Saved validated summary: {args.summary}")


if __name__ == "__main__":
    main()
