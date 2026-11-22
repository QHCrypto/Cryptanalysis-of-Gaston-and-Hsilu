#!/usr/bin/env python3
"""Extract integer witnesses, then complete independent CNF assignments.

Requires python-sat. The resulting certificates are checked without it by
verify_certificates.py; a SAT assignment proves existence, not optimality.
"""
import hashlib
import json
from pathlib import Path

from pysat.formula import CNF
from pysat.solvers import Solver

ROOT = Path(__file__).resolve().parent
MASK = (1 << 64) - 1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rol(x, n):
    return ((x << n) | (x >> ((64-n) % 64))) & MASK


def cols(words):
    active = 0
    for w in words:
        active |= w
    return [j for j in range(64) if active >> j & 1]


def export(s):
    source = ROOT / f"e2_gap11_out{s}.sol"
    values = {}
    for line in source.read_text().splitlines():
        if line and not line.startswith("#"):
            name, raw = line.split()
            value = float(raw)
            if abs(value-round(value)) > 1e-7:
                raise ValueError(f"nonintegral variable {name}")
            values[name] = round(value)
    states = {side: [sum(values[f"{side}_{i}_{j}"] << j for j in range(64))
                     for i in range(5)] for side in ("input", "output")}
    p, q = 0, 0
    for i, (east, twin) in enumerate(zip((0, 60, 22, 27, 4), (25, 32, 52, 60, 63))):
        a = rol(states["input"][i], east)
        p ^= a
        q ^= rol(a, twin)
    e = p ^ rol(p, 1) ^ q ^ rol(q, 18)
    assert e == 0x801
    assert sum(values[f"effect_{j}"] << j for j in range(64)) == e
    record = {
        "mode": "difference",
        "condition": "E(rho_east(X)) != 0",
        "metric": "nonzero 5-bit S-box columns",
        "claim_output_active": s,
        "claim_input_at_least": 19-s,
        "claim_holds": False,
        "certificate_type": "explicit_counterexample",
        "minimum_input_proven": False,
        "theta_effect_hex": f"0x{e:016x}",
        "theta_effect_weight": e.bit_count(),
        "theta_input_words": [f"0x{rol(w, d):016x}" for w, d in
                              zip(states["input"], (0, 60, 22, 27, 4))],
        "P_hex": f"0x{p:016x}",
        "Q_hex": f"0x{q:016x}",
        "source_solution": source.name,
        "source_solution_sha256": sha(source),
    }
    for side, words in states.items():
        record[side+"_words"] = [f"0x{w:016x}" for w in words]
        record[side+"_columns"] = cols(words)
        record[side+"_active"] = len(cols(words))
        record[side+"_bit_weight"] = sum(w.bit_count() for w in words)
        record[side+"_column_patterns"] = {
            str(j): f"0x{sum(((w>>j)&1)<<i for i,w in enumerate(words)):02x}"
            for j in cols(words)
        }
    assert record["output_active"] == s and record["input_active"] <= 18-s
    record["output_quasi_cost"] = s+max(s, record["input_active"])
    record["input_quasi_cost"] = record["input_active"]+max(s, record["input_active"])
    cnf_path = ROOT / f"unrestricted_out{s}.cnf"
    mapping = json.loads((ROOT / f"unrestricted_out{s}.map.json").read_text())
    assumptions = [var if (states["input"][i] >> j) & 1 else -var
                   for i, row in enumerate(mapping["x"]) for j, var in enumerate(row)]
    cnf = CNF(from_file=str(cnf_path))
    with Solver(name="cadical195", bootstrap_with=cnf.clauses) as solver:
        if not solver.solve(assumptions=assumptions):
            raise ValueError("independent SAT encoding rejected the MILP witness")
        model = solver.get_model()
    assert set(abs(v) for v in model) == set(range(1, cnf.nv+1))
    assignment = ROOT / f"out{s}.sat.model"
    assignment.write_text("s SATISFIABLE\nv " + " ".join(map(str, model)) + " 0\n")
    record["cnf"] = cnf_path.name
    record["cnf_sha256"] = sha(cnf_path)
    record["cnf_model"] = assignment.name
    record["cnf_model_sha256"] = sha(assignment)
    record["runtime_seconds"] = json.loads(
        (ROOT / f"e2_gap11_out{s}.status.json").read_text())["runtime_seconds"]
    (ROOT / f"out{s}.counterexample.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"output {s}, input {record['input_active']}: counterexample and CNF assignment exported")
    return record


if __name__ == "__main__":
    records = [export(s) for s in (7, 8, 9)]
    (ROOT / "results.json").write_text(json.dumps({
        "all_three_claims_refuted": True,
        "minimum_input_counts_computed": False,
        "output_quasi_branch_upper_bound": min(r["output_quasi_cost"] for r in records),
        "input_quasi_branch_upper_bound": min(r["input_quasi_cost"] for r in records),
        "results": records,
    }, indent=2) + "\n")
