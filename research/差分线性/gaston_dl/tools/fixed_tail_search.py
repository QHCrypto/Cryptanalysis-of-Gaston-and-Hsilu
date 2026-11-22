#!/usr/bin/env python3
"""Import verified Gaston trails from ePrint 2024/1867 and search a two-chi prefix.

The prefix correlation is exact. Multiplication by a SINGLE linear trail's
squared correlation is a ranking score, not the full permutation correlation.
Only Python's standard library and pdftotext are required.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
import re
import subprocess
from pathlib import Path

from linear_seed_scan import (MASK, EAST, WEST, TWIN, active, chi, column,
                              embed, encode, linear, rol, rotate, transpose, xor)

RC = (0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5, 0x96, 0x87, 0x78, 0x69, 0x5a, 0x4b)
S = tuple(chi(x) for x in range(32))
ZERO = (0,) * 5


def sxor(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def dot(a, b):
    return xor(x & y for x, y in zip(a, b)).bit_count() & 1


def state_chi(a):
    return tuple(x ^ ((~a[(i + 1) % 5] & MASK) & a[(i + 2) % 5])
                 for i, x in enumerate(a))


def theta_transpose(a):
    p = xor(a)
    return tuple(a[i] ^ rol(p, -23) ^ rol(p, -24)
                 ^ rol(p, -23 - TWIN[i]) ^ rol(p, -41 - TWIN[i]) for i in range(5))


def lat(a, b):
    return sum((-1) ** (((a & x).bit_count() + (b & S[x]).bit_count()) & 1)
               for x in range(32))


def checked_round(a, b):
    values = [lat(column(a, j), column(b, j)) for j in range(64)]
    if not all(values):
        raise ValueError("zero LAT entry in imported trail")
    weight = sum(2 * (5 - (abs(v).bit_length() - 1)) for v in values)
    return {"input": a, "output": b, "w_squared": weight,
            "sign": math.prod(1 if v > 0 else -1 for v in values),
            "active_columns": active(a)}


def import_paper(path):
    text = subprocess.check_output(["pdftotext", "-layout", str(path), "-"], text=True)
    trails = {}
    for key, section, expected in (("D2", "D.2    A 2-round Iterative", [24, 24]),
                                   ("D3", "D.3     A 4-round", [22, 18, 20, 22])):
        pages = [(i + 1, p) for i, p in enumerate(text.split("\f")) if section in p]
        if len(pages) != 1:
            raise ValueError(f"cannot uniquely locate {section}")
        page_no, page = pages[0]
        rows = re.findall(r"^\s*([.1]{64})\s*$", page, re.M)
        if len(rows) != (4 * len(expected) - 2) * 5:
            raise ValueError(f"unexpected mask count in {key}")
        # Printed leftmost bit is 63 (MSB). The opposite convention fails L^T.
        states = [tuple(int(s.replace(".", "0"), 2) for s in rows[i:i + 5])
                  for i in range(0, len(rows), 5)]
        rounds = [checked_round(states[k], states[k + 1])
                  for k in range(0, len(states) - 1, 4)]
        assert [r["w_squared"] for r in rounds] == expected
        for k in range(0, len(states) - 2, 4):
            assert states[k + 1] == tuple(rol(v, -s) for v, s in zip(states[k + 2], EAST))
            assert states[k + 2] == theta_transpose(states[k + 3])
            assert states[k + 3] == tuple(rol(v, -s) for v, s in zip(states[k + 4], WEST))
            assert states[k + 1] == transpose(states[k + 4])
        # Verify each boundary mask against the forward L on all 320 basis bits.
        for r in rounds:
            mask = r["input"]
            for lane, j in itertools.product(range(5), range(64)):
                x = embed(j, 1 << lane)
                assert dot(linear(x), mask) == dot(x, transpose(mask))
        trails[key] = {"physical_pdf_page": page_no, "printed_page": page_no - 1,
                       "printed_rows": rows, "rounds": rounds}
    # D.2 repeats at a single chi boundary; verify rather than infer from its title.
    assert trails["D2"]["rounds"][0] == trails["D2"]["rounds"][1]
    a = trails["D2"]["rounds"][0]
    assert transpose(a["input"]) == a["output"]
    return trails


class ExactPrefix:
    """Exact E (-1)^[gamma . D_Delta (chi o (L+rc) o chi)(X)].

    For fixed beta=D_Delta chi(X), the second derivative is affine in
    Y=L chi(X)+rc: gamma.D_{L beta}chi(Y)=h.Y+k. Pull H=L^T h back.
    Inactive first-chi columns have independent uniform outputs, so a branch
    vanishes unless H is supported inside Delta's active columns. On each
    remaining column, sum (-1)^[H_j.chi(x)] over D_Delta_j chi(x)=beta_j.
    This keeps conditional value dependence and uses integer counts throughout.
    """

    def __init__(self, gamma):
        self.gamma = gamma
        self.ddt = [[0] * 32 for _ in range(32)]
        self.walsh = [[[0] * 32 for _ in range(32)] for _ in range(32)]
        for a, x in itertools.product(range(32), repeat=2):
            b = S[x] ^ S[x ^ a]
            self.ddt[a][b] += 1
            for h in range(32):
                self.walsh[a][b][h] += (-1) ** ((h & S[x]).bit_count() & 1)
        local = {}
        for b, m in itertools.product(range(32), repeat=2):
            f = lambda x: ((m & (S[x] ^ S[x ^ b])).bit_count() & 1)
            k = f(0)
            h = sum((f(1 << i) ^ k) << i for i in range(5))
            assert all(f(x) == (k ^ ((h & x).bit_count() & 1)) for x in range(32))
            local[b, m] = h
        self.entries = []
        for j in range(64):
            row = []
            for b in range(32):
                d = linear(embed(j, b))
                h = [0] * 5
                for k in range(64):
                    v = local[column(d, k), column(gamma, k)]
                    for i in range(5):
                        h[i] |= ((v >> i) & 1) << k
                h = tuple(h)
                row.append((d, h, transpose(h)))
            self.entries.append(row)

    def evaluate(self, delta, rc, details=False):
        positions = [j for j in range(64) if column(delta, j)]
        if not positions:
            raise ValueError("input difference must be nonzero")
        outside = MASK ^ sum(1 << j for j in positions)
        choices = [[b for b in range(1, 32) if self.ddt[column(delta, j)][b]]
                   for j in positions]
        total = zero_value_total = 0
        branches = []
        for bs in itertools.product(*choices):
            H = ZERO
            for j, b in zip(positions, bs):
                H = sxor(H, self.entries[j][b][2])
            if any(v & outside for v in H):
                continue
            count = math.prod(self.walsh[column(delta, j)][b][column(H, j)]
                              for j, b in zip(positions, bs))
            if not count:
                continue
            d = h = beta = ZERO
            for j, b in zip(positions, bs):
                d = sxor(d, self.entries[j][b][0])
                h = sxor(h, self.entries[j][b][1])
                beta = sxor(beta, embed(j, b))
            phase = dot(self.gamma, sxor(state_chi(d), state_chi(ZERO)))
            phase ^= (h[0] & rc).bit_count() & 1
            value = (-1) ** phase * count
            total += value
            if h == ZERO:
                zero_value_total += value
            if details:
                pcount = math.prod(self.ddt[column(delta, j)][b]
                                   for j, b in zip(positions, bs))
                branches.append({"first_chi_output_difference": encode(beta),
                                 "middle_chi_input_difference": encode(d),
                                 "ddt_count": pcount, "denominator": 32 ** len(positions),
                                 "middle_uniform_dlct": (-1) ** phase if h == ZERO else 0,
                                 "exact_signed_count": value})
        return {"numerator": total, "denominator": 32 ** len(positions),
                "zero_value_numerator": zero_value_total, "branches": branches}

    def single_column_scan(self, rc):
        results = []
        disagreements = 0
        for j, a in itertools.product(range(64), range(1, 32)):
            value = self.evaluate(embed(j, a), rc)
            disagreements += value["numerator"] != value["zero_value_numerator"]
            if value["numerator"]:
                results.append((value["numerator"], j, a))
        results.sort(key=lambda x: (-abs(x[0]), x[1], x[2]))
        return results, disagreements

    def two_column_zero_certificate(self):
        # If every output-difference branch vanishes, all 31^2 input patterns
        # at every pair of columns vanish. This test makes no sign assumptions.
        eligible = 0
        witnesses = []
        for j in range(64):
            for k in range(j + 1, 64):
                outside = MASK ^ ((1 << j) | (1 << k))
                for b in range(1, 32):
                    left = self.entries[j][b][2]
                    for c in range(1, 32):
                        right = self.entries[k][c][2]
                        if not any((x ^ y) & outside for x, y in zip(left, right)):
                            eligible += 1
                            if len(witnesses) < 8:
                                witnesses.append([j, b, k, c])
        return {"input_candidates": math.comb(64, 2) * 31 ** 2,
                "eligible_output_difference_branches": eligible,
                "all_two_active_column_prefix_correlations_zero": eligible == 0,
                "witnesses_if_certificate_fails": witnesses,
                "scope": "exactly two active input columns, fixed two-chi prefix and gamma; independent of rc"}


def tails_from(trails):
    tails = []
    for length in range(2, 7):
        tails.append((f"D2_repeat_{length}", "D.2/Table 13", [trails["D2"]["rounds"][0]] * length, 1))
    for start in range(3):
        for length in range(2, 5 - start):
            rounds = trails["D3"]["rounds"][start:start + length]
            tails.append((f"D3_chi{start + 1}_to_{start + length}", "D.3/Table 14", rounds, start + 1))
    return tails


def write_json(path, obj):
    with path.open("x", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
        f.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    trails = import_paper(args.paper)
    imported = {"source": str(args.paper.resolve()),
                "sha256": hashlib.sha256(args.paper.read_bytes()).hexdigest(),
                "printed_bit_order": "left-to-right = bit 63 to bit 0",
                "checks": "all local LAT entries, each rho_east/theta/rho_west transpose, all 320 forward-L basis bits",
                "trails": trails}
    write_json(args.output / "imported_trails.json", imported)
    cache = {}
    reports = []
    for name, source, rounds, offset in tails_from(trails):
        r_l = len(rounds)
        total_rounds = r_l + 2
        start = 12 - total_rounds
        gamma = transpose(rounds[0]["input"])
        if gamma not in cache:
            cache[gamma] = ExactPrefix(gamma)
        search = cache[gamma]
        candidates, disagreements = search.single_column_scan(RC[start + 1])
        weight = sum(r["w_squared"] for r in rounds)
        sign = math.prod(r["sign"] * (-1) ** ((r["input"][0] & RC[start + 2 + i]).bit_count() & 1)
                         for i, r in enumerate(rounds))
        report = {"id": name, "source": source, "source_first_chi": offset,
                  "R": total_rounds, "r_d": 1, "r_m": 1, "r_l": r_l,
                  "start_round_index": start, "first_L_omitted": True,
                  "remaining_constants": [f"0x{x:02x}" for x in RC[start + 1:]],
                  "tail_rounds": rounds, "tail_w_squared": weight,
                  "tail_signed_correlation_product": f"{sign} * 2^-{weight // 2}",
                  "join_mask": encode(gamma), "output_mask": encode(rounds[-1]["output"]),
                  "searched_single_column_inputs": 1984, "nonzero_inputs": len(candidates),
                  "exact_vs_zero_value_disagreements": disagreements,
                  "candidate_scope": "fixed unrotated published tail; all 64*31 single-active-column input differences",
                  "top_candidates": []}
        for numerator, j, a in candidates[:16]:
            delta = embed(j, a)
            v = search.evaluate(delta, RC[start + 1], details=True)
            report["top_candidates"].append({"input_column": j, "input_pattern": a,
                                              "input_difference": encode(delta),
                                              "prefix": v,
                                              "log2_abs_score": math.log2(abs(numerator) / 32) - weight,
                                              "score_sign": 1 if numerator > 0 else -1})
        if name == "D3_chi1_to_4":
            report["two_column_certificate"] = search.two_column_zero_certificate()
        reports.append(report)
        print(name, "R", total_rounds, "w_L_squared", weight, "nonzero", len(candidates),
              "best log2 score", report["top_candidates"][0]["log2_abs_score"] if candidates else None,
              flush=True)
    best = []
    for R in range(4, 9):
        eligible = [r for r in reports if r["R"] == R and r["top_candidates"]]
        winner = max(eligible, key=lambda r: r["top_candidates"][0]["log2_abs_score"])
        c = winner["top_candidates"][0]
        best.append({"R": R, "r_d": 1, "r_m": 1, "r_l": R - 2,
                     "tail_id": winner["id"], "tail_w_squared": winner["tail_w_squared"],
                     "prefix_numerator": c["prefix"]["numerator"], "prefix_denominator": 32,
                     "log2_abs_score": c["log2_abs_score"],
                     "input_column": c["input_column"], "input_pattern": c["input_pattern"],
                     "input_words": ";".join(c["input_difference"]["words"]),
                     "join_words": ";".join(winner["join_mask"]["words"]),
                     "output_words": ";".join(winner["output_mask"]["words"])})
    write_json(args.output / "search.json", {
        "method": "exact two-chi prefix correlation times one published linear trail's squared LAT product",
        "limits": ["The score is not a measurement, bound, or exact full-R correlation.",
                   "It retains one identical pair of tail trails. Other trails and cross terms can change or cancel it.",
                   "Only first L is omitted. The leading L of the tail is included by gamma=L^T(alpha).",
                   "Single-column prefix optimum is exhaustive within the fixed imported tail only.",
                   "No optimality over all linear tails or all input differences is claimed."],
        "reports": reports, "best_by_R": best})
    with (args.output / "summary.csv").open("x", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(best[0]))
        writer.writeheader()
        writer.writerows(best)
    # Sampling jobs validate the strong prefix, not the unobservably small full-R score.
    jobs = []
    seen = set()
    for r in reports:
        c = r["top_candidates"][0] if r["top_candidates"] else None
        delta = tuple(int(x, 16) for x in c["input_difference"]["words"]) if c else embed(0, 1)
        gamma = tuple(int(x, 16) for x in r["join_mask"]["words"])
        identity = (delta, gamma)
        if identity in seen:
            continue
        seen.add(identity)
        jobs.append((r["id"], 2, r["start_round_index"],
                     c["prefix"]["numerator"] if c else 0, 32, delta, gamma))
    with (args.output / "prefix_jobs.txt").open("x") as f:
        for name, n, s, num, den, delta, gamma in jobs:
            f.write(f"{name} {n} {s} {num} {den} " + " ".join(f"{w:016x}" for w in delta + gamma) + "\n")
    print("Saved", args.output.resolve(), flush=True)


if __name__ == "__main__":
    main()
