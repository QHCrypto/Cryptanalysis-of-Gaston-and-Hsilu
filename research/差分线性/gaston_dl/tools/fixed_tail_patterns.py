#!/usr/bin/env python3
"""Try existing six/seven-column L patterns before the full D.3 linear tail."""
import argparse
import hashlib
import json
import math

from pathlib import Path
from fixed_tail_search import ExactPrefix, RC, column, embed, encode, linear, rotate, sxor, write_json, ZERO


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", type=Path, required=True)
    parser.add_argument("--patterns", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = next(r for r in json.loads(args.search.read_text())["reports"] if r["id"] == "D3_chi1_to_4")
    gamma = tuple(int(v, 16) for v in report["join_mask"]["words"])
    engine = ExactPrefix(gamma)
    candidates = []
    sources = []
    for n in (6, 7):
        path = args.patterns / f"diff_in{n}.json"
        pattern = json.loads(path.read_text())
        beta = tuple(int(w, 16) for w in pattern["input_words"])
        out = tuple(int(w, 16) for w in pattern["output_words"])
        assert linear(beta) == out
        sources.append({"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        # The stored pattern is the input of L, i.e. first chi OUTPUT difference.
        # Choose one maximum-probability predecessor independently in each column.
        delta = ZERO
        pcount = 1
        active = 0
        for j in range(64):
            b = column(beta, j)
            if b:
                a = max(range(1, 32), key=lambda a: engine.ddt[a][b])
                delta = sxor(delta, embed(j, a))
                pcount *= engine.ddt[a][b]
                active += 1
        for shift in range(64):
            a = rotate(delta, shift)
            value = engine.evaluate(a, RC[report["start_round_index"] + 1])
            candidates.append({"pattern_columns": n, "shift": shift,
                               "first_chi_output_seed": encode(rotate(beta, shift)),
                               "input_difference": encode(a),
                               "seed_ddt_probability": {"numerator": pcount, "denominator": 32 ** active},
                               "prefix": value,
                               "log2_abs_score": (math.log2(abs(value["numerator"]) / value["denominator"]) - 82)
                               if value["numerator"] else None})
        print(f"{n}-column pattern: 64 rotations evaluated", flush=True)
    write_json(args.output, {"tail": report["id"], "R": 6, "split": [1, 1, 4],
                            "sources": sources, "candidates": candidates,
                            "nonzero_candidates": sum(c["prefix"]["numerator"] != 0 for c in candidates),
                            "scope": "one locally best DDT predecessor per supplied L input pattern and all 64 rotations; not all multi-column differences"})
    print("nonzero:", sum(c["prefix"]["numerator"] != 0 for c in candidates), flush=True)


if __name__ == "__main__":
    main()
