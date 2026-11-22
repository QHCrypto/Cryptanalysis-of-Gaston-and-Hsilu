#!/usr/bin/env python3
"""Exact SAT decision instances. Requires python-sat only to encode cardinality."""
import argparse
import itertools
import json
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

EAST = (0, 60, 22, 27, 4)
WEST = (0, 56, 31, 46, 43)
TWIN = (25, 32, 52, 60, 63)


def generate(s, p, effect, prefix):
    pool = IDPool()
    clauses = []
    word = lambda name: [pool.id((name, j)) for j in range(64)]
    x = [word(f"x{i}") for i in range(5)]
    y = [word(f"y{i}") for i in range(5)]
    a, b, e = word("P"), word("Q"), word("E")
    ai, ao = word("active_input"), word("active_output")

    def xor(inputs, output):
        inputs = list(inputs)
        while len(inputs) > 3:
            t = pool.id()
            xor(inputs[:3], t)
            inputs = [t] + inputs[3:]
        vs = inputs + [output]
        for values in itertools.product((0, 1), repeat=len(vs)):
            if sum(values) % 2:
                clauses.append([-v if bit else v for v, bit in zip(vs, values)])

    def atmost(lits, bound):
        clauses.extend(CardEnc.atmost(lits=lits, bound=bound, vpool=pool,
                                     encoding=EncType.seqcounter).clauses)

    for j in range(64):
        xor([x[i][(j-EAST[i]) % 64] for i in range(5)], a[j])
        xor([x[i][(j-EAST[i]-TWIN[i]) % 64] for i in range(5)], b[j])
        xor([a[j], a[(j-1) % 64], b[j], b[(j-18) % 64]], e[j])
        for i in range(5):
            xor([x[i][(j-EAST[i]-WEST[i]) % 64], e[(j-WEST[i]-23) % 64]], y[i][j])
        for st, act in ((x, ai), (y, ao)):
            for i in range(5):
                clauses.append([-st[i][j], act[j]])
            clauses.append([-act[j]] + [st[i][j] for i in range(5)])
    atmost(ai, p)
    atmost(ao, s)
    atmost([-v for v in ao], 64-s)
    # Valid redundant effect bounds; nonzero E has even Hamming weight.
    atmost(e, p+s)
    atmost([-v for v in e], 62)
    if effect == "any":
        clauses.append([e[0]])  # simultaneous cyclic rotation representative
    else:
        fixed = int(effect, 16)
        if fixed == 0:
            raise ValueError("fixed E must be nonzero")
        clauses.extend([[v if fixed >> j & 1 else -v] for j, v in enumerate(e)])
    prefix = Path(prefix)
    with prefix.with_suffix(".cnf").open("w") as f:
        f.write(f"c Gaston difference: output={s}, input<={p}, E!=0\n")
        f.write(f"c E search restriction: {effect}; LSB bit numbering\n")
        f.write(f"p cnf {pool.top} {len(clauses)}\n")
        for clause in clauses:
            f.write(" ".join(map(str, clause)) + " 0\n")
    prefix.with_suffix(".map.json").write_text(json.dumps({
        "output_active": s, "input_max": p, "effect_search": effect,
        "x": x, "y": y, "E": e, "P": a, "Q": b,
        "variables": pool.top, "clauses": len(clauses),
    }, indent=2) + "\n")
    print(prefix, pool.top, len(clauses))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=int, required=True)
    parser.add_argument("--input-max", type=int, required=True)
    parser.add_argument("--effect", default="any")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    generate(args.output, args.input_max, args.effect, args.prefix)
