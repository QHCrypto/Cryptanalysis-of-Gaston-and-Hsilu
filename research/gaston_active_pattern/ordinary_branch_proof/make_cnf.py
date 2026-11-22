#!/usr/bin/env python3
"""Exact ordinary branch-number decision CNF; optional search restrictions."""
import argparse
import itertools
import json
from pathlib import Path
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

def generate(args):
    pool, clauses = IDPool(), []
    def word(name):
        return [pool.id((name, j)) for j in range(64)]
    x, y = [word(f"x{i}") for i in range(5)], [word(f"y{i}") for i in range(5)]
    p, q, e = word("P"), word("Q"), word("E")
    ai, ao = word("active_input"), word("active_output")
    def xor(vs, r):
        vs = list(vs)
        while len(vs) > 3:
            t = pool.id()
            xor(vs[:3], t)
            vs = [t]+vs[3:]
        vs.append(r)
        for bits in itertools.product((0, 1), repeat=len(vs)):
            if sum(bits) % 2:
                clauses.append([-v if b else v for v, b in zip(vs, bits)])
    def atmost(vs, bound):
        clauses.extend(CardEnc.atmost(lits=vs, bound=bound, vpool=pool,
                                     encoding=EncType.seqcounter).clauses)
    east, west, twin = (0,60,22,27,4), (0,56,31,46,43), (25,32,52,60,63)
    for j in range(64):
        xor([x[i][(j-east[i]) % 64] for i in range(5)], p[j])
        xor([x[i][(j-east[i]-twin[i]) % 64] for i in range(5)], q[j])
        xor([p[j], p[(j-1) % 64], q[j], q[(j-18) % 64]], e[j])
        for i in range(5):
            xor([x[i][(j-east[i]-west[i]) % 64], e[(j-west[i]-23) % 64]], y[i][j])
            # Redundant propagation clause for the active columns.
            clauses.append([-e[j], ai[(j+23-east[i]) % 64], ao[(j+23+west[i]) % 64]])
        for st, active in ((x, ai), (y, ao)):
            clauses.extend([[-st[i][j], active[j]] for i in range(5)])
            clauses.append([-active[j]]+[st[i][j] for i in range(5)])
    atmost(ai+ao, args.total)
    atmost(e, args.total)
    atmost([-v for v in e], 62)  # E has even weight and is nonzero
    if args.input is not None:
        atmost(ai, args.input)
        atmost([-v for v in ai], 64-args.input)
    if args.output is not None:
        atmost(ao, args.output)
        atmost([-v for v in ao], 64-args.output)
    if args.reduced_domain:
        # Conditional on separately established small-endpoint bounds:
        # any total <= 15 counterexample must have p >= 6 and s >= 7.
        atmost([-v for v in ai], 58)
        atmost([-v for v in ao], 57)
    if args.effect == "any":
        clauses.append([e[0]])  # no additional rotation anchor
    else:
        fixed = int(args.effect, 16)
        if not 0 < fixed < 1 << 64:
            raise ValueError("E must be a nonzero 64-bit hexadecimal integer")
        clauses.extend([[v if fixed >> j & 1 else -v] for j, v in enumerate(e)])
    if args.effect_weight is not None:
        atmost(e, args.effect_weight)
        atmost([-v for v in e], 64-args.effect_weight)
    prefix = Path(args.prefix)
    with prefix.with_suffix(".cnf").open("w") as f:
        f.write(f"c Gaston: E!=0, wc(X)+wc(Y)<={args.total}\n")
        f.write(f"p cnf {pool.top} {len(clauses)}\n")
        for c in clauses:
            f.write(" ".join(map(str, c))+" 0\n")
    prefix.with_suffix(".map.json").write_text(json.dumps({
        "parameters": vars(args), "x":x, "y":y, "E":e, "P":p, "Q":q,
        "active_input":ai, "active_output":ao, "variables":pool.top, "clauses":len(clauses),
    }, indent=2)+"\n")
    print(prefix, pool.top, len(clauses), flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--total", type=int, required=True)
    parser.add_argument("--input", type=int)
    parser.add_argument("--output", type=int)
    parser.add_argument("--effect", default="any")
    parser.add_argument("--effect-weight", type=int)
    parser.add_argument("--reduced-domain", action="store_true")
    parser.add_argument("--prefix", required=True)
    generate(parser.parse_args())
