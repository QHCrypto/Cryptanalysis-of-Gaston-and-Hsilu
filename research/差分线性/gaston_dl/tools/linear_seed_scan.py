#!/usr/bin/env python3
"""Enumerate 2->2 active-column linear masks of Gaston's L, modulo rotation.

This is a bounded single-layer enumeration, not a multi-round trail search.
For a linear-mask pair (alpha,beta), alpha=L^T(beta).
"""
import argparse
import itertools
import json
import math
from pathlib import Path

MASK = (1 << 64) - 1
EAST = (0, 60, 22, 27, 4)
WEST = (0, 56, 31, 46, 43)
TWIN = (25, 32, 52, 60, 63)


def rol(x, shift):
    shift %= 64
    return ((x << shift) | (x >> ((-shift) % 64))) & MASK


def xor(values):
    value = 0
    for x in values:
        value ^= x
    return value


def linear(x):
    a = tuple(rol(v, s) for v, s in zip(x, EAST))
    p = xor(a)
    q = xor(rol(v, s) for v, s in zip(a, TWIN))
    effect = rol(p ^ rol(p, 1) ^ q ^ rol(q, 18), 23)
    return tuple(rol(v ^ effect, s) for v, s in zip(a, WEST))


def transpose(beta):
    d = tuple(rol(v, -s) for v, s in zip(beta, WEST))
    p = xor(d)
    return tuple(rol(d[i] ^ rol(p, -23) ^ rol(p, -24)
                     ^ rol(p, -23 - TWIN[i]) ^ rol(p, -41 - TWIN[i]), -EAST[i])
                 for i in range(5))


def embed(col, pattern):
    return tuple(((pattern >> i) & 1) << col for i in range(5))


def column(state, j):
    return sum(((v >> j) & 1) << i for i, v in enumerate(state))


def active(state):
    return (state[0] | state[1] | state[2] | state[3] | state[4]).bit_count()


def rotate(state, amount):
    return tuple(rol(x, amount) for x in state)


def chi(x):
    return sum((((x >> i) & 1) ^ ((1 ^ ((x >> ((i + 1) % 5)) & 1))
                                 & ((x >> ((i + 2) % 5)) & 1))) << i for i in range(5))


def encode(state):
    return {
        "words": [f"0x{v:016x}" for v in state],
        "bits": [[i, j] for i, v in enumerate(state) for j in range(64) if v >> j & 1],
        "active_columns": [j for j in range(64) if column(state, j)],
    }


def scan_two_column_seeds():
    """Exhaust all two-column beta, fixing one column at zero by rotation."""
    table = {(j, v): transpose(embed(j, v)) for j in range(64) for v in range(1, 32)}
    seeds = set()
    hits = 0
    for j, p, q in itertools.product(range(1, 64), range(1, 32), range(1, 32)):
        alpha = tuple(a ^ b for a, b in zip(table[0, p], table[j, q]))
        if active(alpha) != 2:
            continue
        beta = tuple(a ^ b for a, b in zip(embed(0, p), embed(j, q)))
        hits += 1
        seeds.add(min((rotate(alpha, s), rotate(beta, s)) for s in range(64)))
    return sorted(seeds), hits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    lat = [[sum((-1) ** (((u & x).bit_count() + (v & chi(x)).bit_count()) % 2)
                for x in range(32)) for v in range(32)] for u in range(32)]
    seeds, hits = scan_two_column_seeds()
    records = []
    for alpha, beta in seeds:
        assert transpose(beta) == alpha
        # Check the transpose relation against the forward linear operation on
        # every input basis vector, independently of the transpose formula.
        for lane, j in itertools.product(range(5), range(64)):
            x = embed(j, 1 << lane)
            lhs = xor(a & b for a, b in zip(linear(x), beta)).bit_count() % 2
            rhs = (alpha[lane] >> j) & 1
            assert lhs == rhs
        before = sum(-math.log2(max(abs(lat[u][column(alpha, j)]) for u in range(32)) / 32)
                     for j in range(64))
        after = sum(-math.log2(max(abs(x) for x in lat[column(beta, j)]) / 32)
                    for j in range(64))
        records.append({"alpha": encode(alpha), "beta": encode(beta),
                        "best_chi_L_chi_trail_w_abs": before + after,
                        "best_chi_L_chi_trail_w_squared": 2 * (before + after)})
    edges = []
    for i, (_, beta) in enumerate(seeds):
        for k, (alpha, _) in enumerate(seeds):
            for shift in range(64):
                target = rotate(alpha, shift)
                if all(lat[column(beta, j)][column(target, j)] for j in range(64)):
                    edges.append({"from": i, "to": k, "shift": shift})
    alpha = embed(0, 1)
    alpha = tuple(a ^ b for a, b in zip(alpha, embed(4, 2)))
    beta = tuple(a ^ b for a, b in zip(embed(0, 1), embed(56, 2)))
    assert transpose(beta) == alpha
    result = {
        "scope": "all single-L linear mask pairs with 2 active columns on each side, modulo common rotation",
        "relation": "alpha = L^T(beta)",
        "bit_index": "(lane,bit), zero-based, LSB=0",
        "enumerated_output_patterns_with_column_zero_active": 63 * 31 * 31,
        "hits_with_column_zero_active": hits,
        "rotation_orbits": len(seeds),
        "direct_chi_compatible_edges_between_orbits": edges,
        "example": {"alpha": encode(alpha), "beta": encode(beta)},
        "seeds": records,
        "limits": [
            "These witnesses attain active-column sum 4; the global lower bound is not proved by this scan alone.",
            "The chi-L-chi weights are products of local LAT correlations, not the full linear hull correlation.",
            "A linear tail following a middle chi includes a leading L: pull the chi-L-chi input mask back through L^T.",
            "Fixed round constants may change signs; this scan evaluates mask feasibility and absolute local weights only."
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({key: result[key] for key in (
        "enumerated_output_patterns_with_column_zero_active", "hits_with_column_zero_active", "rotation_orbits")}, indent=2))
    print("direct chi-compatible edges:", len(edges))
    print("two-chi absolute-correlation weights:", sorted({x['best_chi_L_chi_trail_w_abs'] for x in records}))
    print(args.output.resolve())


if __name__ == "__main__":
    main()
