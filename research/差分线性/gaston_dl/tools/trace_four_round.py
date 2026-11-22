#!/usr/bin/env python3
"""Extract and check the README's 1+1+2 four-round contribution diagram.

The two-chi prefix is exact. The four-round quantity is a splice score using
one fixed linear trail, not a measured or proven full-permutation correlation.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from fixed_tail_search import ExactPrefix, RC, dot, lat
from linear_seed_scan import chi, column, embed, encode, linear, transpose

ROOT = Path(__file__).resolve().parents[1]


def decode(encoded):
    return tuple(int(word, 16) for word in encoded["words"])


def main():
    source = ROOT / "results/fixed_paper_tail/search.json"
    search = json.loads(source.read_text())
    best = next(row for row in search["best_by_R"] if row["R"] == 4)
    report = next(row for row in search["reports"] if row["id"] == best["tail_id"])
    candidate = next(row for row in report["top_candidates"]
                     if row["input_column"] == best["input_column"]
                     and row["input_pattern"] == best["input_pattern"])
    delta = decode(candidate["input_difference"])
    gamma = decode(report["join_mask"])
    assert delta == embed(1, 0x10)
    assert (report["r_d"], report["r_m"], report["r_l"]) == (1, 1, 2)
    assert report["remaining_constants"] == ["0x69", "0x5a", "0x4b"]
    exact = ExactPrefix(gamma).evaluate(delta, 0x69, details=True)
    assert exact == candidate["prefix"]
    prefix_correlation = Fraction(exact["numerator"], exact["denominator"])
    assert prefix_correlation == Fraction(1, 2)

    tail = []
    previous = gamma
    signed_product = Fraction(1)
    for index, saved in enumerate(report["tail_rounds"]):
        alpha, beta = tuple(saved["input"]), tuple(saved["output"])
        assert transpose(alpha) == previous
        # Check the mask boundary against the forward L, independently of L^T.
        for i in range(5):
            for j in range(64):
                basis = embed(j, 1 << i)
                assert dot(linear(basis), alpha) == dot(basis, previous)
        local = [Fraction(lat(column(alpha, j), column(beta, j)), 32)
                 for j in range(64)]
        corr = math.prod(local)
        assert corr * corr == Fraction(1, 1 << saved["w_squared"])
        rc = RC[report["start_round_index"] + 2 + index]
        phase = (-1) ** ((alpha[0] & rc).bit_count() & 1)
        signed_product *= phase * corr
        tail.append({"chi_number": 3 + index, "input": encode(alpha),
                     "output": encode(beta), "chi_correlation": str(corr),
                     "w_squared": saved["w_squared"],
                     "preceding_constant": f"0x{rc:02x}", "constant_sign": phase})
        previous = beta
    assert encode(previous) == report["output_mask"]
    assert signed_product == Fraction(-1, 1 << 19)
    squared = signed_product * signed_product
    assert squared == Fraction(1, 1 << 38)

    branches = []
    for branch in exact["branches"]:
        beta = decode(branch["first_chi_output_difference"])
        middle_delta = decode(branch["middle_chi_input_difference"])
        assert linear(beta) == middle_delta
        pattern = column(beta, 1)
        ddt_count = sum(chi(x) ^ chi(x ^ 0x10) == pattern for x in range(32))
        assert ddt_count == branch["ddt_count"] == 8
        factors = []
        for j in range(64):
            a, b = column(middle_delta, j), column(gamma, j)
            values = [(-1) ** ((b & (chi(x) ^ chi(x ^ a))).bit_count() & 1)
                      for x in range(32)]
            factors.append(Fraction(sum(values), 32))
        connection = math.prod(factors)
        assert connection == branch["middle_uniform_dlct"] == 1
        # Local derivatives are constant; negative local factors cancel.
        assert all(abs(factor) == 1 for factor in factors)
        contribution = Fraction(branch["exact_signed_count"], branch["denominator"])
        assert contribution == Fraction(1, 4)
        branches.append({"beta_pattern": pattern,
                         "beta": encode(beta), "middle_delta": encode(middle_delta),
                         "ddt_probability": str(Fraction(ddt_count, 32)),
                         "middle_connection": str(connection),
                         "negative_middle_columns": [j for j, factor in enumerate(factors) if factor == -1],
                         "exact_prefix_contribution": str(contribution),
                         "splice_score": str(contribution * squared)})
    assert [row["beta_pattern"] for row in branches] == [0x10, 0x1c]
    score = prefix_correlation * squared
    assert score == Fraction(1, 1 << 39)
    result = {
        "source": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "tail_id": report["id"], "split": [1, 1, 2], "omit_first_linear": True,
        "input": encode(delta), "join_mask": encode(gamma), "branches": branches,
        "tail": tail, "prefix_correlation": str(prefix_correlation),
        "tail_signed_product": str(signed_product), "tail_squared_product": str(squared),
        "splice_score": str(score), "log2_abs_splice_score": -39,
        "limits": "Exact two-chi prefix; fixed single-trail splice score, not the full four-round correlation."
    }
    output = ROOT / "results/four_round_trace"
    output.mkdir(parents=True, exist_ok=True)
    (output / "trace.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "input": result["input"]["bits"],
        "branches": [{"beta": b["beta"]["bits"], "delta": b["middle_delta"]["bits"]} for b in branches],
        "join": result["join_mask"]["bits"],
        "tail": [{"input": t["input"]["bits"], "output": t["output"]["bits"],
                  "chi_correlation": t["chi_correlation"], "constant_sign": t["constant_sign"]} for t in tail],
        "prefix": str(prefix_correlation), "score": str(score)
    }, separators=(",", ":")))
    print("PASS: DDT branches, constant middle derivatives, exact prefix, three L^T boundaries, all tail LAT entries, round-constant signs and splice score.")


if __name__ == "__main__":
    main()
