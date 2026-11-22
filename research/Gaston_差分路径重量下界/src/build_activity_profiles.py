#!/usr/bin/env python3
"""Import completed local activity results, preserve scope, derive count-level cuts."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT.parent
MASK = (1 << 64) - 1


def rol(x, n):
    return ((x << n) | (x >> ((64 - n) & 63))) & MASK


def linear(x):
    a = [rol(w, s) for w, s in zip(x, (0, 60, 22, 27, 4))]
    p = q = 0
    for w, t in zip(a, (25, 32, 52, 60, 63)):
        p ^= w
        q ^= rol(w, t)
    e = p ^ rol(p, 1) ^ q ^ rol(q, 18)
    y = [rol(w ^ rol(e, 23), s) for w, s in zip(a, (0, 56, 31, 46, 43))]
    return y, e


def active(words):
    union = 0
    for w in words:
        union |= w
    return union.bit_count()


def allowed(x, y, nonzero, forward, reverse):
    return (y >= forward[x] and x >= reverse[y] and
            x + y >= 12 + 4 * nonzero and
            (not nonzero or y != 6 or x >= 13) and
            (nonzero or (x <= 5*y and y <= 5*x)))


def main():
    records = []
    def load(path, scope):
        p = WORK / path
        d = json.loads(p.read_text())
        assert d["mode"] == "difference" and d["proven_optimal"] and d["verified"]
        if d.get("method", "").startswith("exact_enumeration"):
            assert d["enumeration_complete"]
        else:
            assert d["status"] == 2 and d["bound"] == d["objective"] and d["gap"] == 0
        assert d.get("theta_effect_constraint", "any") == scope
        x = [int(w, 16) for w in d["input_words"]]
        y = [int(w, 16) for w in d["output_words"]]
        actual, e = linear(x)
        assert actual == y and active(x) == d["input_active"] and active(y) == d["output_active"]
        assert scope != "nonzero" or e
        records.append({"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                        "scope": scope, "fixed_side": d["fixed_side"], "fixed_count": d["fixed_count"],
                        "lower_bound": d["objective"], "evidence": d.get("method", "saved_Gurobi_OPTIMAL"),
                        "witness_rechecked": True})
        return d
    forward = [0] + [max(1, 12-k) for k in range(1, 65)]
    reverse = forward.copy()
    for k in range(1, 9):
        name = "difference_input1.json" if k == 1 else f"diff_in{k}.json"
        d = load("gaston_active_pattern/result/" + name, "any")
        assert d["fixed_side"] == "input" and d["fixed_count"] == k
        forward[k] = d["objective"]
    names = ["diff_out1_nonzero_witness", "diff_out2", "diff_out3", "diff_out4_gpu", "diff_out5_gpu"]
    for k, name in enumerate(names, 1):
        d = load(f"gaston_active_pattern/direct_enumeration/result/{name}.json", "any")
        assert d["fixed_side"] == "output" and d["fixed_count"] == k
        reverse[k] = d["objective"]
    d = load("gaston_active_pattern/direct_enumeration/result/diff_out6_nonzero.json", "nonzero")
    assert d["fixed_side"] == "output" and d["fixed_count"] == 6 and d["objective"] == 13
    certificate = WORK / "gaston_active_pattern/ordinary_branch_proof/result.json"
    cert = json.loads(certificate.read_text())
    assert cert["proof_complete"] and cert["verification"] == "PASS" and cert["lower_bound"] == 16
    baseline = WORK / "Cube/extracted/gaston_definition_for_cube.txt"
    raw_forward, raw_reverse = forward.copy(), reverse.copy()
    pairs = [(a, b) for a in range(1, 65) for b in range(1, 65)
             if b >= forward[a] and a >= reverse[b] and a+b >= 12]
    # Close the two directional tables without changing their allowed count pairs.
    forward = [0] + [min(b for a, b in pairs if a == k) for k in range(1, 65)]
    reverse = [0] + [min(a for a, b in pairs if b == k) for k in range(1, 65)]
    assert set(pairs) == {(a, b) for a in range(1, 65) for b in range(1, 65)
                         if b >= forward[a] and a >= reverse[b] and a+b >= 12}
    middle_total = [0] + [reverse[k]+k+forward[k] for k in range(1, 65)]
    minima = {}
    for n0, n1 in ((0, 0), (0, 1), (1, 0), (1, 1)):
        best, witnesses = 193, []
        for a in range(1, 65):
            for b in range(1, 65):
                if not allowed(a, b, n0, forward, reverse):
                    continue
                for c in range(1, 65):
                    if allowed(b, c, n1, forward, reverse):
                        total = a + b + c
                        if total < best:
                            best, witnesses = total, []
                        if total == best:
                            witnesses.append([a, b, c])
        minima[f"{n0}{n1}"] = {"minimum_count_sum": best, "weight_lower_bound": 2*best,
                                "count_witnesses_not_state_trails": witnesses}
    assert [minima[k]["minimum_count_sum"] for k in ("00", "01", "10", "11")] == [18, 22, 22, 22]
    report = {
        "status": "PASS",
        "source_scope": "Saved completed results; witnesses rechecked, old optimality computations not replayed",
        "forward_lower_by_exact_input_count": forward,
        "reverse_lower_by_exact_output_count": reverse,
        "raw_forward_lower_by_exact_input_count": raw_forward,
        "raw_reverse_lower_by_exact_output_count": raw_reverse,
        "total_count_lower_by_middle_count": middle_total,
        "unconditional_branch_number_lower_bound": 12,
        "outside_branch_number_lower_bound": 16,
        "outside_output6_input_lower_bound": 13,
        "records": records,
        "branch_sources": [
            {"path": str(certificate), "sha256": hashlib.sha256(certificate.read_bytes()).hexdigest(),
             "evidence": "existing checked combinatorial and DRAT certificates; not replayed now"},
            {"path": str(baseline), "sha256": hashlib.sha256(baseline.read_bytes()).hexdigest(),
             "evidence": "original Gaston paper section 7.6/Table 5, unconditional branch number 12"}],
        "three_round_count_minima_by_nonzero_flags": minima,
        "excluded_incomplete_records": ["gaston_active_pattern/result/diff_in6_outker.json",
                                        "gaston_active_pattern/result/diff_out2.json",
                                        "gaston_active_pattern/result/diff_outker.json"],
    }
    destination = ROOT / "references/activity_profiles.json"
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    digest = hashlib.sha256(destination.read_bytes()).hexdigest()
    header = "#pragma once\n#include <array>\n"
    for name, values in (("FORWARD_ACTIVITY_LOWER", forward), ("REVERSE_ACTIVITY_LOWER", reverse),
                         ("MIDDLE_TOTAL_ACTIVITY_LOWER", middle_total)):
        header += "constexpr std::array<int,65> " + name + "{" + ",".join(map(str, values)) + "};\n"
    header += 'constexpr const char* ACTIVITY_PROFILE_SHA256 = "' + digest + '";\n'
    (ROOT / "src/activity_profile_data.hpp").write_text(header)
    print(json.dumps({"forward_1_to_8": forward[1:9], "reverse_1_to_5": reverse[1:6],
                      "three_round_count_minima": minima}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
