#!/usr/bin/env python3
"""Two differential chi layers (6 -> 6), one middle chi, fixed paper tails.

Exhaust the 64 rotations of a supplied L witness and every DDT output of chi 2.
A meet-in-the-middle equality tests the exact zero-value-mask connection.
The output is a restricted signed contribution, not a full-prefix correlation.
"""
import argparse
import collections
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

from fixed_tail_search import (MASK, RC, S, ZERO, column, dot, embed, encode,
                              lat, linear, rotate, state_chi, sxor, transpose,
                              write_json, xor)


def flat(state):
    return sum(word << (64 * i) for i, word in enumerate(state))


def rank(rows):
    basis = {}
    for v in rows:
        while v:
            k = v.bit_length() - 1
            if k in basis:
                v ^= basis[k]
            else:
                basis[k] = v
                break
    return len(basis)


def active_mask(state):
    return state[0] | state[1] | state[2] | state[3] | state[4]


def derivative_mask(gamma, delta):
    return tuple((gamma[(i - 1) % 5] & delta[(i + 1) % 5])
                 ^ (gamma[(i - 2) % 5] & delta[(i - 1) % 5]) for i in range(5))


def assemble(positions, patterns):
    return tuple(xor(embed(j, v)[i] for j, v in zip(positions, patterns)) for i in range(5))


def ddt_table():
    ddt = [[0] * 32 for _ in range(32)]
    for a, x in itertools.product(range(32), repeat=2):
        ddt[a][S[x] ^ S[x ^ a]] += 1
    return ddt


def chi_probability(a, b, ddt):
    return math.prod(Fraction(ddt[column(a, j)][column(b, j)], 32) for j in range(64))


def second_round_rank_certificate(first_input, second_input):
    # First DDT event involves only the active first-chi columns. Every other
    # first-chi output column is independent and uniform, even conditioned on it.
    # The second-chi derivative equations are affine in Y=L Z+rc.
    outside = MASK ^ active_mask(first_input)
    all_rows, outside_rows = [], []
    for j in range(64):
        a = column(second_input, j)
        for i in range(5):
            h = (((a >> ((i + 1) % 5)) & 1) << ((i + 2) % 5))
            h ^= (((a >> ((i + 2) % 5)) & 1) << ((i + 1) % 5))
            if not h:
                continue
            local = embed(j, h)
            all_rows.append(flat(local))
            pulled = transpose(local)
            outside_rows.append(flat(tuple(w & outside for w in pulled)))
    r_all, r_outside = rank(all_rows), rank(outside_rows)
    return {"second_derivative_rank": r_all, "rank_on_unconditioned_first_output_columns": r_outside,
            "ddt_product_is_exact_for_this_two_round_trail": r_all == r_outside,
            "reason": "Full rank on untouched first-chi columns makes all feasible second derivative values uniform, conditioned on the first DDT event; constants only change affine offsets."}


def match_outputs(positions, options, table):
    split = len(positions) // 2
    left = collections.defaultdict(list)
    for values in itertools.product(*options[:split]):
        key = xor(table[j][v] for j, v in zip(positions[:split], values))
        left[key].append(values)
    hits = []
    for values in itertools.product(*options[split:]):
        key = xor(table[j][v] for j, v in zip(positions[split:], values))
        hits.extend(a + values for a in left.get(key, []))
    return sorted(hits)


def fraction_record(value):
    return {"numerator": value.numerator, "denominator": value.denominator,
            "log2_abs": math.log2(abs(value)) if value else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tails", type=Path, required=True)
    parser.add_argument("--pattern", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    old = json.loads(args.tails.read_text())
    pattern = json.loads(args.pattern.read_text())
    b1_seed = tuple(int(w, 16) for w in pattern["input_words"])
    a2_seed = tuple(int(w, 16) for w in pattern["output_words"])
    assert linear(b1_seed) == a2_seed
    assert active_mask(b1_seed).bit_count() == active_mask(a2_seed).bit_count() == 6
    ddt = ddt_table()
    first_patterns = []
    first_positions = [j for j in range(64) if column(b1_seed, j)]
    alternatives = []
    for j in first_positions:
        b = column(b1_seed, j)
        best = max(ddt[a][b] for a in range(1, 32))
        choices = [a for a in range(1, 32) if ddt[a][b] == best]
        first_patterns.append(choices[0])
        alternatives.append(choices)
    d0_seed = assemble(first_positions, first_patterns)
    cache, reports, groups = {}, [], {}
    for old_report in old["reports"]:
        gamma = tuple(int(w, 16) for w in old_report["join_mask"]["words"])
        if gamma not in cache:
            table = [[flat(derivative_mask(gamma, linear(embed(j, b))))
                      for b in range(32)] for j in range(64)]
            scan = []
            for shift in range(64):
                d0, b1, a2 = (rotate(s, shift) for s in (d0_seed, b1_seed, a2_seed))
                positions = [j for j in range(64) if column(a2, j)]
                options = [[b for b in range(1, 32) if ddt[column(a2, j)][b]] for j in positions]
                hits = match_outputs(positions, options, table)
                # Independently brute-force the small 4^6 search for every shift.
                brute = [values for values in itertools.product(*options)
                         if derivative_mask(gamma, linear(assemble(positions, values))) == ZERO]
                assert hits == sorted(brute)
                certificate = second_round_rank_certificate(d0, a2)
                assert certificate["ddt_product_is_exact_for_this_two_round_trail"]
                p1 = chi_probability(d0, b1, ddt)
                connections = []
                signed_sum = Fraction(0)
                for index, values in enumerate(hits):
                    b2 = assemble(positions, values)
                    a3 = linear(b2)
                    # Independent local enumeration verifies both feasibility and sign.
                    local = [sum((-1) ** ((column(gamma, j) & (S[x] ^ S[x ^ column(a3, j)])).bit_count() & 1)
                                 for x in range(32)) for j in range(64)]
                    assert all(abs(v) == 32 for v in local)
                    g = math.prod(v // 32 for v in local)
                    assert g == (-1) ** dot(gamma, state_chi(a3))
                    p2 = chi_probability(a2, b2, ddt)
                    contribution = p1 * p2 * g
                    signed_sum += contribution
                    connections.append({"index": index, "b2": encode(b2), "a3": encode(a3),
                                        "p2": fraction_record(p2), "middle_correlation": g,
                                        "two_round_trail_probability": fraction_record(p1 * p2),
                                        "exact_event_contribution_to_three_round_prefix": fraction_record(contribution)})
                scan.append({"rotation": shift, "d0": encode(d0), "b1": encode(b1), "a2": encode(a2),
                             "p1": fraction_record(p1), "rank_certificate": certificate,
                             "second_DDT_output_combinations": math.prod(map(len, options)),
                             "connections": connections, "restricted_signed_prefix_sum": fraction_record(signed_sum)})
            cache[gamma] = scan
            group_id = f"join_{len(groups)}"
            groups[group_id] = {"gamma": encode(gamma), "rotations": scan,
                                "feasible_trails": sum(len(s["connections"]) for s in scan),
                                "noncancelling_rotations": [s["rotation"] for s in scan if s["restricted_signed_prefix_sum"]["numerator"]]}
        group_id = next(k for k, v in groups.items() if v["gamma"]["words"] == encode(gamma)["words"])
        scan = cache[gamma]
        new_R = old_report["R"] + 1
        tail_weight = old_report["tail_w_squared"]
        # Tail constants remain identical: adding a prefix round shifts its start back by one.
        assert RC[12 - new_R + 3:] == RC[old_report["start_round_index"] + 2:]
        previous, signed_tail = gamma, Fraction(1)
        for i, row in enumerate(old_report["tail_rounds"]):
            a, b = tuple(row["input"]), tuple(row["output"])
            assert transpose(a) == previous
            c = math.prod(Fraction(lat(column(a, j), column(b, j)), 32) for j in range(64))
            rc = RC[12 - new_R + 3 + i]
            signed_tail *= c * (-1) ** ((a[0] & rc).bit_count() & 1)
            previous = b
        assert signed_tail ** 2 == Fraction(1, 1 << tail_weight)
        valid = [s for s in scan if s["restricted_signed_prefix_sum"]["numerator"]]
        best = max(valid, key=lambda s: abs(Fraction(s["restricted_signed_prefix_sum"]["numerator"], s["restricted_signed_prefix_sum"]["denominator"]))) if valid else None
        result = {"id": old_report["id"], "R": new_R, "r_d": 2, "r_m": 1, "r_l": old_report["r_l"],
                  "differential_active_chi_counts": [6, 6], "start_round_index": 12 - new_R,
                  "remaining_constants": [f"0x{x:02x}" for x in RC[13 - new_R:]],
                  "connection_group": group_id, "tail_w_squared": tail_weight,
                  "tail_signed_LAT_product": fraction_record(signed_tail),
                  "tail_rounds": old_report["tail_rounds"], "output_mask": old_report["output_mask"],
                  "feasible_trails": groups[group_id]["feasible_trails"],
                  "noncancelling_rotations": groups[group_id]["noncancelling_rotations"],
                  "best_rotation": best["rotation"] if best else None,
                  "best_restricted_splice_score": fraction_record(Fraction(best["restricted_signed_prefix_sum"]["numerator"], best["restricted_signed_prefix_sum"]["denominator"]) * signed_tail ** 2) if best else None}
        reports.append(result)
        print(result["id"], "R", new_R, "trails", result["feasible_trails"],
              "noncancelling shifts", result["noncancelling_rotations"],
              "score", result["best_restricted_splice_score"], flush=True)
    data = {"sources": [{"path": str(p.resolve()), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in (args.tails, args.pattern)],
            "scope": "One supplied 6->6 L witness, all 64 common rotations; one locally optimal first-DDT predecessor, all second-DDT output differences; all 11 existing paper tail candidates.",
            "first_predecessor_choices_per_active_column": alternatives,
            "first_L_omitted": True, "groups": groups, "reports": reports,
            "checks": ["MITM equals independent full 4^6 enumeration for each of 64 rotations and 4 distinct masks",
                       "every second DDT event has full rank 12 on unconditioned first-output columns",
                       "all connecting middle derivatives verified by direct 32-value local enumeration",
                       "fixed tail LAT, transpose and constant signs rechecked"],
            "limits": ["Not a search over all 6->6 propagation patterns.",
                       "Two-round trail probabilities and connecting event signs are exact for this witness by the rank certificate.",
                       "Restricted sums include only the fixed first-DDT output branch and the middle DLCT +/-1 branches.",
                       "The full three-chi prefix correlation is not calculated. A zero restricted sum does not prove zero full-prefix correlation.",
                       "Splice scores retain one fixed tail trail pair; they are not full-R correlation estimates or bounds."]}
    write_json(args.output / "search.json", data)
    # Conditional sampling jobs test selected trails under the exact first-DDT event.
    jobs = []
    for rid, shift in (("D3_chi2_to_3", 24), ("D3_chi2_to_4", 24), ("D3_chi2_to_3", 7)):
        r = next(r for r in reports if r["id"] == rid)
        group = groups[r["connection_group"]]
        s = group["rotations"][shift]
        c = s["connections"][0]
        states = [s["d0"], s["b1"], c["b2"], group["gamma"]]
        jobs.append(f'{rid}_rot{shift} {r["start_round_index"]} {c["middle_correlation"]} '
                    + ' '.join(w for state in states for w in state["words"]))
    (args.output / "conditional_jobs.txt").write_text('\n'.join(jobs) + '\n')
    print("Saved", args.output.resolve(), flush=True)


if __name__ == "__main__":
    main()
