#!/usr/bin/env python3
"""Compare continuous relaxations only; no integer search or solve-time caps."""
import argparse
import json
from pathlib import Path
import gurobipy as gp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original_mps", type=Path)
    parser.add_argument("strengthened_mps", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Output report already exists")
    variants = [
        ("original", args.original_mps, False, False, False),
        ("or_hull_only", args.strengthened_mps, False, True, False),
        ("profiles_without_combined_count_cuts", args.strengthened_mps, True, False, False),
        ("profiles_with_combined_count_cuts", args.strengthened_mps, True, False, True),
        ("all_strengthening", args.strengthened_mps, True, True, True),
    ]
    results = []
    for label, source, profiles, or_hull, combined in variants:
        print("\nLP comparison:", label, flush=True)
        original = gp.read(str(source))
        remove = []
        for row in original.getConstrs():
            name = row.ConstrName
            if name.startswith("profile_") and (
                    not profiles or (not combined and
                                     ("three_round_total" in name or name == "profile_middle_total"))):
                remove.append(row)
            elif name.startswith("chi") and ("_input_bit" in name or "_output_bit" in name) and not or_hull:
                remove.append(row)
        original.remove(remove)
        original.update()
        model = original.relax()
        model.Params.Threads = 16
        model.Params.Method = 1
        model.optimize()
        if model.Status != gp.GRB.OPTIMAL:
            raise RuntimeError(f"LP was not solved to optimality: {label}, {model.Status}")
        results.append({"variant": label, "continuous_relaxation_objective": model.ObjVal,
                        "runtime_seconds": model.Runtime, "constraints": model.NumConstrs,
                        "variables": model.NumVars, "iterations": model.IterCount})
        model.dispose()
        original.dispose()
    report = {"scope": "Continuous LP relaxations before MIP presolve/cuts; no integer search",
              "threads_per_sequential_solve": 16, "results": results}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
