#!/usr/bin/env python3
"""Independent exact checks. Python standard library only, no MILP/SAT solver.

Checks explicit bit-destination propagation, a separately implemented word
circuit, GF(2) matrix inversion, E!=0, column counts, every DIMACS clause,
the link between each SAT assignment and the witness, and file hashes.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MASK = (1 << 64) - 1
EAST = (0, 60, 22, 27, 4)
WEST = (0, 56, 31, 46, 43)
TWIN = (25, 32, 52, 60, 63)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rol(x, n):
    n %= 64
    return ((x << n) | (x >> ((64-n) % 64))) & MASK


def pack(words):
    return sum(w << (64*i) for i, w in enumerate(words))


def bit_image(position):
    """Track one input bit through all destinations, with XOR cancellation."""
    i, j = divmod(position, 64)
    z = (j+EAST[i]) % 64
    image = 1 << (64*i+(z+WEST[i]) % 64)
    for d in (0, 1, TWIN[i], TWIN[i]+18):
        for r in range(5):
            image ^= 1 << (64*r+(z+d+23+WEST[r]) % 64)
    return image


def invert(rows):
    augmented = [r | (1 << (320+i)) for i, r in enumerate(rows)]
    for j in range(320):
        pivot = next((i for i in range(j, 320) if augmented[i] >> j & 1), None)
        require(pivot is not None, "singular forward matrix")
        augmented[j], augmented[pivot] = augmented[pivot], augmented[j]
        for i in range(320):
            if i != j and augmented[i] >> j & 1:
                augmented[i] ^= augmented[j]
    require(all((r & ((1 << 320)-1)) == 1 << i for i, r in enumerate(augmented)),
            "matrix reduction failed")
    return [r >> 320 for r in augmented]


def active_columns(words):
    active = 0
    for w in words:
        active |= w
    return [j for j in range(64) if active >> j & 1]


def hash_check(record, key, filekey):
    data = (ROOT / record[filekey]).read_bytes()
    require(hashlib.sha256(data).hexdigest() == record[key], f"hash mismatch {filekey}")


def check_cnf(record, states):
    mapping = json.loads((ROOT / f"unrestricted_out{record['output_active']}.map.json").read_text())
    assignment = {}
    for line in (ROOT / record["cnf_model"]).read_text().splitlines():
        if line.startswith("v "):
            for lit in map(int, line[2:].split()):
                if lit:
                    require(abs(lit) not in assignment, "duplicate assignment variable")
                    assignment[abs(lit)] = lit > 0
    count, pending, variables, declared = 0, [], None, None
    for line in (ROOT / record["cnf"]).read_text().splitlines():
        if not line or line.startswith("c"):
            continue
        if line.startswith("p "):
            _, kind, n, m = line.split()
            require(kind == "cnf", "unexpected formula format")
            variables, declared = int(n), int(m)
            continue
        for lit in map(int, line.split()):
            if lit:
                pending.append(lit)
            else:
                require(any(assignment.get(abs(v)) == (v > 0) for v in pending),
                        f"unsatisfied clause {count+1}")
                count += 1
                pending = []
    require(not pending and count == declared, "incorrect DIMACS clause count")
    require(set(assignment) == set(range(1, variables+1)), "incomplete SAT assignment")
    for side, key in (("input", "x"), ("output", "y")):
        reconstructed = [sum(int(assignment[v]) << j for j, v in enumerate(row))
                         for row in mapping[key]]
        require(reconstructed == states[side], f"CNF/state mismatch: {side}")
    require(sum(int(assignment[v]) << j for j, v in enumerate(mapping["E"])) ==
            int(record["theta_effect_hex"], 16), "CNF/E mismatch")
    return count


def main():
    images = [bit_image(j) for j in range(320)]
    rows = [sum(((images[j] >> i) & 1) << j for j in range(320)) for i in range(320)]
    inverse = invert(rows)
    results = []
    for s in (7, 8, 9):
        record = json.loads((ROOT / f"out{s}.counterexample.json").read_text())
        states = {side: [int(w, 16) for w in record[side+"_words"]]
                  for side in ("input", "output")}
        for side, words in states.items():
            require(len(words) == 5 and all(0 <= w <= MASK for w in words), "invalid state")
            columns = active_columns(words)
            require(columns == record[side+"_columns"], f"{side} column list")
            require(len(columns) == record[side+"_active"], f"{side} column count")
            require(sum(w.bit_count() for w in words) == record[side+"_bit_weight"],
                    f"{side} bit count")
            patterns = {str(j): f"0x{sum(((w>>j)&1)<<i for i,w in enumerate(words)):02x}"
                        for j in columns}
            require(patterns == record[side+"_column_patterns"], f"{side} column patterns")
        x, y = states["input"], states["output"]
        actual = 0
        for j, image in enumerate(images):
            if pack(x) >> j & 1:
                actual ^= image
        require(actual == pack(y), "bit-destination propagation")
        recovered = sum(((r & pack(y)).bit_count() & 1) << i for i, r in enumerate(inverse))
        require(recovered == pack(x), "inverse matrix propagation")
        a = [rol(w, d) for w, d in zip(x, EAST)]
        p, q = 0, 0
        for w, d in zip(a, TWIN):
            p ^= w
            q ^= rol(w, d)
        e = p ^ rol(p, 1) ^ q ^ rol(q, 18)
        require([rol(w ^ rol(e, 23), d) for w, d in zip(a, WEST)] == y,
                "word circuit propagation")
        require(a == [int(w, 16) for w in record["theta_input_words"]], "theta input")
        require(p == int(record["P_hex"], 16) and q == int(record["Q_hex"], 16), "P/Q mismatch")
        require(e != 0 and e == int(record["theta_effect_hex"], 16), "not outside kernel")
        require(e.bit_count() == record["theta_effect_weight"], "E bit count")
        require(record["output_active"] == s, "wrong output count")
        require(record["claim_output_active"] == s and record["claim_input_at_least"] == 19-s,
                "wrong claim")
        require(record["input_active"] < 19-s, "not a counterexample")
        require(record["claim_holds"] is False and record["minimum_input_proven"] is False,
                "incorrect status labels")
        require(record["output_quasi_cost"] == s+max(s, record["input_active"]), "output quasi cost")
        require(record["input_quasi_cost"] == record["input_active"]+max(s, record["input_active"]),
                "input quasi cost")
        for key, filekey in (("source_solution_sha256", "source_solution"),
                             ("cnf_sha256", "cnf"), ("cnf_model_sha256", "cnf_model")):
            hash_check(record, key, filekey)
        clauses = check_cnf(record, states)
        result = {"output_active": s, "input_active": record["input_active"],
                  "theta_effect_hex": f"0x{e:016x}", "claim_refuted": True,
                  "clauses_checked": clauses, "all_checks_passed": True}
        print(f"PASS output={s}, input={record['input_active']}, E=0x{e:016x}, "
              f"all {clauses} CNF clauses satisfied")
        results.append(result)
    aggregate = json.loads((ROOT / "results.json").read_text())
    require(aggregate["results"] == [json.loads((ROOT / f"out{s}.counterexample.json").read_text())
                                     for s in (7, 8, 9)], "aggregate mismatch")
    summary = {"verification": "PASS", "method": "exact_integer_arithmetic",
               "forward_matrix_rank": 320, "requires_solver": False, "results": results}
    (ROOT / "verification.json").write_text(json.dumps(summary, indent=2)+"\n")
    print("PASS: matrix rank 320; forward, inverse, effect, column counts, hashes and SAT certificates.")


if __name__ == "__main__":
    main()
