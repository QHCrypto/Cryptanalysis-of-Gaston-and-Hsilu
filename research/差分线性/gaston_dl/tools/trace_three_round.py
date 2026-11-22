#!/usr/bin/env python3
"""Exact rational decomposition of the existing three-round approximate model.

The decomposition is exact *within the zero-value-mask model*. It does not
replace its independence assumptions by a proof about the actual permutation.
"""
import csv
import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORD_MASK = (1 << 64) - 1
EAST = (0, 60, 22, 27, 4)
WEST = (0, 56, 31, 46, 43)
THETA = (25, 32, 52, 60, 63)


def rot(x, r):
    return ((x << r) | (x >> ((-r) & 63))) & WORD_MASK


def linear(a):
    a = [rot(x, r) for x, r in zip(a, EAST)]
    p = q = 0
    for x, t in zip(a, THETA):
        p ^= x
        q ^= rot(x, t)
    p ^= rot(p, 1)
    q ^= rot(q, 18)
    p = rot(p ^ q, 23)
    return [rot(x ^ p, r) for x, r in zip(a, WEST)]


def embed(j, a):
    return [((a >> i) & 1) << j for i in range(5)]


def column(a, j):
    return sum(((a[i] >> j) & 1) << i for i in range(5))


def dot(a, b):
    return sum((x & y).bit_count() for x, y in zip(a, b)) % 2


def transpose(m):
    return [sum(dot(linear(embed(j, 1 << i)), m) << j for j in range(64)) for i in range(5)]


def chi(x):
    return sum((((x >> i) & 1) ^ ((1 ^ ((x >> ((i + 1) % 5)) & 1)) &
                                 ((x >> ((i + 2) % 5)) & 1))) << i for i in range(5))


def chi_correlation(u, v):
    return Fraction(sum((-1) ** ((u & x).bit_count() + (v & chi(x)).bit_count()) for x in range(32)), 32)


def dlct(a, v):
    return Fraction(sum((-1) ** ((v & (chi(x) ^ chi(x ^ a))).bit_count()) for x in range(32)), 32)


def state(a):
    return {
        "words": [f"{x:016x}" for x in a],
        "active": [[i, j] for i in range(5) for j in range(64) if (a[i] >> j) & 1],
    }


def main():
    output = ROOT / "results/three_round_trace"
    output.mkdir(parents=True, exist_ok=True)
    counts = {}
    for x in range(32):
        b = chi(x) ^ chi(x ^ 1)
        counts[b] = counts.get(b, 0) + 1
    tails = [(u, chi_correlation(u, 8)) for u in range(32) if chi_correlation(u, 8)]
    terms = []
    total = Fraction(0)
    for b, count in sorted(counts.items()):
        delta = linear(embed(0, b))
        for u, c in tails:
            mask = transpose(embed(19, u))
            connection = Fraction(1)
            local = []
            for j in range(64):
                if column(mask, j):
                    a, v = column(delta, j), column(mask, j)
                    coefficient = dlct(a, v)
                    connection *= coefficient
                    local.append({"column": j, "difference": a, "mask": v, "correlation": str(coefficient)})
            contribution = Fraction(count, 32) * connection * c * c
            total += contribution
            terms.append({
                "beta_pattern": b, "u_pattern": u,
                "first_probability": str(Fraction(count, 32)),
                "middle_connection": str(connection), "tail_correlation": str(c),
                "tail_squared": str(c * c), "contribution": str(contribution),
                "delta_after_chi1": state(embed(0, b)),
                "delta_before_chi2": state(delta),
                "mask_after_chi2": state(mask),
                "mask_before_chi3": state(embed(19, u)),
                "local_connections": local,
            })
    with (ROOT / "results/table13/summary.csv").open() as stream:
        row = next(r for r in csv.DictReader(stream) if r["rounds"] == "3")
    assert total == Fraction(row["model_rho"]) == Fraction(-3, 16)
    nonzero = [t for t in terms if t["contribution"] != "0"]
    assert [(t["beta_pattern"], t["u_pattern"], t["contribution"]) for t in nonzero] == [
        (1, 8, "-1/16"), (1, 9, "-1/16"), (17, 8, "-1/16")]
    trace = {
        "description": "Exact rational decomposition of the approximate HNW-zero-value-mask-rd1-rl0 score",
        "input": state(embed(0, 1)), "output": state(embed(19, 8)),
        "rounds": 3, "omit_first_linear": True,
        "subsequent_constants": ["0x5a", "0x4b"],
        "beta_patterns": sorted(counts), "u_patterns": [u for u, _ in tails],
        "terms": terms, "model_correlation": str(total),
        "experimental_correlation": row["experimental_rho"],
        "samples": int(row["samples"]),
    }
    (output / "trace.json").write_text(json.dumps(trace, ensure_ascii=False, indent=2) + "\n")
    with (output / "contributions.csv").open("w", newline="") as stream:
        names = ["beta_pattern", "u_pattern", "first_probability", "middle_connection", "tail_correlation", "tail_squared", "contribution"]
        writer = csv.DictWriter(stream, names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(terms)
    visual_data = {
        "input": trace["input"]["active"], "output": trace["output"]["active"],
        "branches": [{
            "beta": t["beta_pattern"], "u": t["u_pattern"], "c": t["tail_correlation"],
            "states": [trace["input"]["active"], t["delta_after_chi1"]["active"],
                       t["delta_before_chi2"]["active"], t["mask_after_chi2"]["active"],
                       t["mask_before_chi3"]["active"], trace["output"]["active"]],
            "negative_columns": [p["column"] for p in t["local_connections"] if p["correlation"] == "-1"],
        } for t in nonzero],
    }
    (output / "visual_data.json").write_text(json.dumps(visual_data, separators=(",", ":")) + "\n")
    print(json.dumps(visual_data, separators=(",", ":")))
    print("PASS: all 16 contributions computed exactly as rational numbers; sum=-3/16 matches saved estimator result")


if __name__ == "__main__":
    main()
