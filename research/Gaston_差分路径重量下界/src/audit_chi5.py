#!/usr/bin/env python3
"""Exhaustive chi5 DDT/rank audit and read-only audit of the existing CNF."""

import hashlib
import json
import runpy
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OLD_CNF = ROOT.parent / "Gaston_trail/pysat_kissat/sbox_cnf.py"


def chi(x):
    return sum(
        (((x >> i) & 1)
         ^ ((1 ^ ((x >> ((i + 1) % 5)) & 1))
            & ((x >> ((i + 2) % 5)) & 1))) << i
        for i in range(5)
    )


def gf2_rank(vectors):
    basis = {}
    for v in vectors:
        while v:
            pivot = v.bit_length() - 1
            if pivot in basis:
                v ^= basis[pivot]
            else:
                basis[pivot] = v
                break
    return len(basis)


def satisfies(clauses, a, b):
    bits = [(a >> i) & 1 for i in range(5)]
    bits += [(b >> i) & 1 for i in range(5)]
    return all(any(bits[abs(v) - 1] == (v > 0) for v in c)
               for c in clauses)


def main():
    table = list(map(chi, range(32)))
    assert len(set(table)) == 32
    ddt = [[0] * 32 for _ in range(32)]
    for a in range(32):
        for x in range(32):
            ddt[a][table[x] ^ table[x ^ a]] += 1

    ranks = []
    for a in range(32):
        constant = table[a] ^ table[0]
        columns = [table[1 << i] ^ table[(1 << i) ^ a] ^ constant
                   for i in range(5)]
        rank = gf2_rank(columns)
        ranks.append(rank)
        for x in range(32):
            derivative = constant
            for i in range(5):
                if (x >> i) & 1:
                    derivative ^= columns[i]
            assert derivative == table[x] ^ table[x ^ a]
        assert all(n == 0 or n == 2 ** (5 - rank) for n in ddt[a])

    reverse_min = [min(ranks[a] for a in range(32) if ddt[a][b])
                   for b in range(32)]
    fuzzy = [0, 2, 3, 3, 4, 4]
    assert all(fuzzy[a.bit_count()] <= ranks[a] for a in range(32))

    source = runpy.run_path(str(OLD_CNF))
    assert source["SBOX_TABLE"] == tuple(table)
    restricted = source["SBOX_CNF"]["diff"]
    assert len(restricted) == 47
    assert restricted[-1] == (-1, -2, -3, -4, -5)
    full = restricted[:-1]
    assert all(satisfies(full, a, b) == bool(ddt[a][b])
               for a in range(32) for b in range(32))
    assert all(not satisfies(restricted, a, b) or bool(ddt[a][b])
               for a in range(32) for b in range(32))
    excluded = [[a, b] for a in range(32) for b in range(32)
                if ddt[a][b] and not satisfies(restricted, a, b)]

    output = {
        "status": "PASS",
        "scope": "Exhaustive local chi5 facts; no new multiround lower bound",
        "bit_convention": "integer bit i is row i, indices modulo 5",
        "chi_table": table,
        "ddt": ddt,
        "forward_weight_by_input": ranks,
        "minimum_weight_by_fixed_output": reverse_min,
        "forward_weight_histogram": dict(Counter(ranks)),
        "reverse_minimum_histogram": dict(Counter(reverse_min)),
        "weights_by_input_hamming_weight": {
            h: sorted({ranks[a] for a in range(32) if a.bit_count() == h})
            for h in range(6)
        },
        "all_legal_outputs_have_weight_equal_to_derivative_rank": True,
        "fuzzy_cost_by_hamming_weight": fuzzy,
        "fuzzy_cost_is_valid_lower_bound": True,
        "supported_pair_count": sum(bool(n) for row in ddt for n in row),
        "old_cnf_path": str(OLD_CNF),
        "old_cnf_sha256": hashlib.sha256(OLD_CNF.read_bytes()).hexdigest(),
        "first_46_clauses_equal_ddt_support": True,
        "pairs_excluded_by_47th_clause": excluded,
    }
    destination = ROOT / "results/chi5_audit.json"
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "output": str(destination),
                      "supported_pair_count": output["supported_pair_count"],
                      "excluded_pair_count": len(excluded),
                      "forward_weight_histogram": dict(Counter(ranks))},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
