#!/usr/bin/env python3
"""Independent checks of the exported C++ MILP (requires Python gurobipy).

The reference matrix is assembled by tracking each input bit's destinations,
without calling the C++ word-level implementation. Solver logs remain enabled.
"""
import json
from pathlib import Path
import random
import re
import signal
import subprocess
import sys
import tempfile

import gurobipy as gp

ROOT = Path(__file__).resolve().parent
PROGRAM = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "gaston_active_pattern"
BITS = 320
MASK = (1 << 64) - 1


def reference_matrix():
    east = (0, 60, 22, 27, 4)
    west = (0, 56, 31, 46, 43)
    twin = (25, 32, 52, 60, 63)
    columns = []
    for row in range(5):
        for bit in range(64):
            a = (bit + east[row]) % 64
            image = 1 << (64 * row + (a + west[row]) % 64)
            for shift in (0, 1, twin[row], twin[row] + 18):
                for out_row in range(5):
                    image ^= 1 << (64 * out_row + (a + shift + 23 + west[out_row]) % 64)
            columns.append(image)
    rows = [sum(((columns[c] >> r) & 1) << c for c in range(BITS)) for r in range(BITS)]
    return rows, columns


def inverse(rows):
    a = [row | (1 << (BITS + i)) for i, row in enumerate(rows)]
    for c in range(BITS):
        pivot = next(r for r in range(c, BITS) if (a[r] >> c) & 1)
        a[c], a[pivot] = a[pivot], a[c]
        for r in range(BITS):
            if r != c and (a[r] >> c) & 1:
                a[r] ^= a[c]
    assert all((row & ((1 << BITS) - 1)) == 1 << i for i, row in enumerate(a))
    return [row >> BITS for row in a]


def apply(rows, x):
    return sum(((row & x).bit_count() % 2) << i for i, row in enumerate(rows))


def columns(x):
    result = 0
    for i in range(5):
        result |= (x >> (64 * i)) & MASK
    return result


def weight(x):
    return columns(x).bit_count()


def json_state(data, side):
    return sum(int(word, 16) << (64 * i) for i, word in enumerate(data[side + "_words"]))


def reference_effect(value):
    # Equation (3): XOR destinations of each nonzero input bit directly.
    # Input is BEFORE rho_east; the E word is BEFORE the u=23 shift.
    east, twin = (0, 60, 22, 27, 4), (25, 32, 52, 60, 63)
    result = 0
    for row in range(5):
        word = (value >> (64 * row)) & MASK
        for bit in range(64):
            if word >> bit & 1:
                for shift in (0, 1, twin[row], twin[row] + 18):
                    result ^= 1 << ((bit + east[row] + shift) % 64)
    return result


def check_json(path, mode, side, count, expected, forward, transpose, effect_constraint="any"):
    data = json.loads(path.read_text())
    x, y = json_state(data, "input"), json_state(data, "output")
    assert data["mode"] == mode and data["fixed_side"] == side
    assert data["fixed_count"] == count
    assert data["verified"]
    assert data["proven_optimal"] == (data["status"] == gp.GRB.OPTIMAL)
    if expected is not None:
        assert data["status"] == gp.GRB.OPTIMAL and data["proven_optimal"]
        assert data["objective"] == expected == data["bound"] and data["gap"] == 0
    else:
        assert data["status"] in (gp.GRB.OPTIMAL, gp.GRB.INTERRUPTED)
        # An interrupt during presolve may precede the first finite bound.
        if data["bound"] is not None:
            assert data["bound"] <= data["objective"]
        if data["gap"] is not None:
            assert data["gap"] >= 0
        if data["proven_optimal"]:
            assert data["bound"] == data["objective"] and data["gap"] == 0
    assert apply(forward, x) == y if mode == "difference" else apply(transpose, y) == x
    assert data["total_active"] == weight(x) + weight(y)
    if count is None:
        assert side is None and data["input_bit_weight_min"] == 1
        assert x != 0 and weight(x) + weight(y) == data["objective"]
        assert data["objective_side"] == "input+output"
    else:
        assert weight(x if side == "input" else y) == count
        assert weight(y if side == "input" else x) == expected
        assert data["objective_side"] == ("output" if side == "input" else "input")
    for end, value in (("input", x), ("output", y)):
        assert data[end + "_active"] == weight(value)
        assert data[end + "_bit_weight"] == value.bit_count()
        assert data[end + "_columns"] == [j for j in range(64) if columns(value) >> j & 1]
    if mode == "difference":
        effect = reference_effect(x)
        assert data["theta_effect_constraint"] == effect_constraint
        assert int(data["theta_effect_hex"], 16) == effect
        assert data["theta_effect_weight"] == effect.bit_count()
        assert data["in_theta_kernel"] == (effect == 0)
        assert effect_constraint == "any" or (effect == 0) == (effect_constraint == "zero")
    else:
        assert "theta_effect_constraint" not in data and "theta_effect_hex" not in data


def check_parity_hulls(model):
    # Audit the actual LP, not a second copy of the C++ facet generator.
    # All variables must be binary; each XOR block must contain ALL odd-subset
    # facets, with unit coefficients and the correct right side.
    assert all(var.VType == gp.GRB.BINARY for var in model.getVars())
    groups = {}
    for constraint in model.getConstrs():
        assert "_xor4_facet_" not in constraint.ConstrName
        match = re.fullmatch(r"(.+)_xor([234])_facet_\d+", constraint.ConstrName)
        if not match:
            continue
        name, arity = match.group(1), int(match.group(2))
        row = model.getRow(constraint)
        assert row.size() == arity + 1 and constraint.Sense == gp.GRB.LESS_EQUAL
        coefficients = {row.getVar(i).VarName: row.getCoeff(i) for i in range(row.size())}
        assert set(coefficients.values()) <= {-1.0, 1.0}
        support = frozenset(coefficients)
        positive = frozenset(var for var, coef in coefficients.items() if coef == 1)
        assert len(positive) % 2 == 1 and constraint.RHS == len(positive) - 1
        if name not in groups:
            groups[name] = (arity, support, set())
        n, variables, facets = groups[name]
        assert n == arity and variables == support and positive not in facets
        facets.add(positive)
    assert groups
    checked_arities = set()
    for arity, variables, facets in groups.values():
        assert len(facets) == 2 ** arity
        if arity in checked_arities:
            continue
        checked_arities.add(arity)
        ordered = sorted(variables)
        for assignment in range(1 << len(ordered)):
            active = {v for i, v in enumerate(ordered) if assignment >> i & 1}
            accepted = all(len(active & subset) - len(active - subset) <= len(subset) - 1
                           for subset in facets)
            assert accepted == (assignment.bit_count() % 2 == 0)
    return groups


def kernel_samples():
    east, twin = (0, 60, 22, 27, 4), (25, 32, 52, 60, 63)

    def before_east(rows):
        return sum((((word >> shift) | (word << ((64-shift) % 64))) & MASK) << (64*i)
                   for i, (word, shift) in enumerate(zip(rows, east)))

    def rho(row):
        return 1 ^ (1 << 1) ^ (1 << twin[row]) ^ (1 << ((twin[row] + 18) % 64))

    # A0=rho1, A1=rho0: E=rho0*rho1 + rho1*rho0=0, but P and Q are NONZERO.
    rows = [rho(1), rho(0), 0, 0, 0]
    row_twin = before_east(rows)
    assert rows[0] ^ rows[1]
    q = 0
    for word, shift in zip(rows, twin):
        q ^= ((word << shift) | (word >> (64-shift))) & MASK
    assert q != 0 and reference_effect(row_twin) == 0
    # Feeding the theta-input coordinates directly as X must not pass E=0.
    wrong_coordinates = sum(word << (64*i) for i, word in enumerate(rows))
    assert reference_effect(wrong_coordinates) != 0
    rows = [(1 << twin[1]) ^ (1 << twin[2]), (1 << twin[0]) ^ (1 << twin[2]),
            (1 << twin[0]) ^ (1 << twin[1]), 0, 0]
    vortex = before_east(rows)
    assert reference_effect(vortex) == 0
    rotated = sum(((((row_twin >> (64*i)) & MASK) >> 1)
                   | (((row_twin >> (64*i)) & 1) << 63)) << (64*i) for i in range(5))
    return [0, (1 << BITS)-1, row_twin, vortex, rotated,
            row_twin ^ 1, wrong_coordinates, 1, 1 << 319]


def check_effect_models(tmp, forward, transpose):
    samples = kernel_samples()
    assert any(0 < reference_effect(x).bit_count() < 64 for x in samples)
    for side in ("input", "output"):
        for policy in ("zero", "nonzero"):
            prefix = tmp / f"effect_{side}_{policy}"
            subprocess.run([str(PROGRAM), "--mode", "difference", f"--{side}-active", "8",
                            "--theta-effect", policy, "--no-symmetry", "--export-only",
                            "--prefix", str(prefix)], check=True)
            model = gp.read(str(prefix.with_suffix(".lp")))
            model.Params.Threads = 2
            hulls = check_parity_hulls(model)
            # Every four-input effect XOR must now be one xor3 + one xor2.
            for j in range(64):
                assert hulls[f"effect_eq_{j}_tree_0"][0] == 3
                assert hulls[f"effect_eq_{j}_root"][0] == 2
            if policy == "zero":
                for j in range(64):
                    c = model.getConstrByName(f"theta_effect_zero_{j}")
                    row = model.getRow(c)
                    assert c.Sense == gp.GRB.EQUAL and c.RHS == 0
                    assert row.size() == 1 and row.getCoeff(0) == 1
                    assert row.getVar(0).VarName == f"effect_{j}"
            else:
                c = model.getConstrByName("theta_effect_nonzero")
                row = model.getRow(c)
                assert c.Sense == gp.GRB.GREATER_EQUAL and c.RHS == 1 and row.size() == 64
                assert {row.getVar(i).VarName for i in range(64)} == {f"effect_{j}" for j in range(64)}
                assert all(row.getCoeff(i) == 1 for i in range(64))
            model.remove(model.getConstrByName("fixed_count"))
            model.setObjective(0.0)
            model.update()
            states = {end: [model.getVarByName(f"{end}_{i}_{j}")
                            for i in range(5) for j in range(64)] for end in ("input", "output")}
            effects = [model.getVarByName(f"effect_{j}") for j in range(64)]
            for x in samples:
                y = apply(forward, x)
                for end, state in (("input", x), ("output", y)):
                    for bit, var in enumerate(states[end]):
                        var.LB = var.UB = (state >> bit) & 1
                model.optimize()
                effect = reference_effect(x)
                accepted = (effect == 0) == (policy == "zero")
                assert model.Status == (gp.GRB.OPTIMAL if accepted else gp.GRB.INFEASIBLE)
                if accepted:
                    assert sum(round(v.X) << j for j, v in enumerate(effects)) == effect
            model.dispose()
            # CLI: E=0 accepts the all-zero pair, E!=0 rejects it, for either side.
            prefix = tmp / f"effect_zero_count_{side}_{policy}"
            completed = subprocess.run(
                [str(PROGRAM), "--mode", "difference", f"--{side}-active", "0",
                 "--theta-effect", policy, "--threads", "2", "--prefix", str(prefix)],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            print(completed.stdout, end="", flush=True)
            if policy == "zero":
                check_json(prefix.with_suffix(".json"), "difference", side, 0, 0,
                           forward, transpose, policy)
            else:
                assert "Model is infeasible" in completed.stdout
                assert not prefix.with_suffix(".json").exists()
                assert not prefix.with_suffix(".sol").exists()
    # End-to-end nonzero report with a nonzero optimum.
    prefix = tmp / "effect_nonzero_input1"
    subprocess.run([str(PROGRAM), "--mode", "difference", "--input-active", "1",
                    "--theta-effect", "nonzero", "--threads", "2", "--prefix", str(prefix)], check=True)
    check_json(prefix.with_suffix(".json"), "difference", "input", 1, 18,
               forward, transpose, "nonzero")
    # Explicit any and omitted option must export identical default models.
    for side in ("input", "output"):
        paths = []
        for explicit in (False, True):
            prefix = tmp / f"effect_any_{side}_{explicit}"
            args = [str(PROGRAM), "--mode", "difference", f"--{side}-active", "1",
                    "--export-only", "--prefix", str(prefix)]
            if explicit:
                args += ["--theta-effect", "any"]
            subprocess.run(args, check=True)
            paths.append(prefix.with_suffix(".lp"))
        assert paths[0].read_bytes() == paths[1].read_bytes()
    print("PASS: theta-effect zero/nonzero in both directions; 36 fixed-state cases "
          "(including E=0 with P,Q nonzero); effect coordinates and bit values; "
          "5 CLI optimizations; any/default model identity.", flush=True)


def check_nonzero_default(tmp, forward, transpose, inv_transpose):
    configurations = [("difference", p) for p in ("any", "zero", "nonzero")] + [("linear", None)]
    samples = kernel_samples() + [1 << 1]
    for mode, policy in configurations:
        for symmetry in (False, True):
            prefix = tmp / f"default_{mode}_{policy}_{symmetry}"
            args = [str(PROGRAM), "--mode", mode, "--export-only", "--prefix", str(prefix)]
            if policy is not None:
                args += ["--theta-effect", policy]
            if not symmetry:
                args += ["--no-symmetry"]
            subprocess.run(args, check=True)
            model = gp.read(str(prefix.with_suffix(".lp")))
            model.Params.Threads = 2
            assert model.getConstrByName("fixed_count") is None
            c = model.getConstrByName("input_nonzero")
            row = model.getRow(c)
            assert c.Sense == gp.GRB.GREATER_EQUAL and c.RHS == 1 and row.size() == 320
            assert {row.getVar(i).VarName for i in range(row.size())} == {
                f"input_{i}_{j}" for i in range(5) for j in range(64)}
            assert all(row.getCoeff(i) == 1 for i in range(row.size()))
            c = model.getConstrByName("rotation_representative")
            if symmetry:
                row = model.getRow(c)
                assert c.Sense == gp.GRB.GREATER_EQUAL and c.RHS == 1 and row.size() == 5
                assert {row.getVar(i).VarName for i in range(5)} == {f"input_{i}_0" for i in range(5)}
                assert all(row.getCoeff(i) == 1 for i in range(5))
            else:
                assert c is None
            objective = model.getObjective()
            assert model.ModelSense == gp.GRB.MINIMIZE and objective.size() == 128
            assert {objective.getVar(i).VarName for i in range(128)} == {
                f"active_{end}_{j}" for end in ("input", "output") for j in range(64)}
            assert all(objective.getCoeff(i) == 1 for i in range(128))
            variables = {end: [model.getVarByName(f"{end}_{i}_{j}")
                               for i in range(5) for j in range(64)] for end in ("input", "output")}
            accepted_count = 0
            for x in samples:
                y = apply(forward if mode == "difference" else inv_transpose, x)
                for end, value in (("input", x), ("output", y)):
                    for bit, var in enumerate(variables[end]):
                        var.LB = var.UB = (value >> bit) & 1
                accepted = x != 0 and (not symmetry or bool(columns(x) & 1))
                if policy in ("zero", "nonzero"):
                    accepted = accepted and (reference_effect(x) == 0) == (policy == "zero")
                model.optimize()
                assert model.Status == (gp.GRB.OPTIMAL if accepted else gp.GRB.INFEASIBLE)
                if accepted:
                    assert model.ObjVal == weight(x) + weight(y)
                accepted_count += accepted
            assert 0 < accepted_count < len(samples)
            model.dispose()
    # Exercise the real CLI report/JSON after the initial incumbent is loaded.
    # Interrupt instead of requiring a global branch-number proof in the tests.
    # The kernel-zero case, which may lack an initial seed, is covered above.
    for mode, policy, symmetry in (("difference", "any", True),
                                   ("difference", "any", False),
                                   ("difference", "nonzero", True),
                                   ("linear", None, True)):
        prefix = tmp / f"default_solved_{mode}_{policy}_{symmetry}"
        args = [str(PROGRAM), "--threads", "2", "--prefix", str(prefix)]
        if mode == "linear":
            args += ["--mode", mode]
        if policy not in (None, "any"):
            args += ["--theta-effect", policy]
        if not symmetry:
            args += ["--no-symmetry"]
        interrupted = False
        with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, bufsize=1) as process:
            try:
                for line in process.stdout:
                    print(line, end="", flush=True)
                    if not interrupted and "Loaded user MIP start with objective" in line:
                        process.send_signal(signal.SIGINT)
                        interrupted = True
                assert process.wait() == 0
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
        assert interrupted
        check_json(prefix.with_suffix(".json"), mode, None, None, None,
                   forward, transpose, policy or "any")
        data = json.loads(prefix.with_suffix(".json").read_text())
        assert data["objective"] <= (19 if mode == "difference" else 18)
        if symmetry:
            assert columns(json_state(data, "input")) & 1
    print("PASS: optional counts; 320-bit input >= 1; exact five-bit column-0 symmetry constraint; "
          "8 exported models / 80 fixed-state cases; zero rejected and multi-column inputs accepted; "
          "input+output objective coefficients and values; 4 interrupted CLI reports; "
          "theta-effect combinations; no-symmetry.", flush=True)


def check_exported_model(path, fixed_side, mapping):
    model = gp.read(str(path))
    model.Params.Threads = 2
    hulls = check_parity_hulls(model)
    if any(var.VarName.startswith("inverse_cse_") for var in model.getVars()):
        # Exclude each shared node's own defining block when counting its uses.
        uses = {}
        for name, (_, support, _) in hulls.items():
            for var in support:
                if var.startswith("inverse_cse_") and var != name:
                    uses[var] = uses.get(var, 0) + 1
        assert any(count >= 2 for count in uses.values())
    # Remove only the count and rotation representative: the propagation and
    # activity constraints under test remain exactly those exported by C++.
    model.remove(model.getConstrByName("fixed_count"))
    model.remove(model.getConstrByName("rotation_representative"))
    model.setObjective(0.0)
    model.update()
    other_side = "output" if fixed_side == "input" else "input"
    fixed = [model.getVarByName(f"{fixed_side}_{i}_{j}") for i in range(5) for j in range(64)]
    other = [model.getVarByName(f"{other_side}_{i}_{j}") for i in range(5) for j in range(64)]
    rng = random.Random(20261014)
    samples = [0, (1 << BITS) - 1, 1, 1 << 319,
               sum(1 << j for j in rng.sample(range(BITS), 9)), rng.getrandbits(BITS)]
    for value in samples:
        for c, var in enumerate(fixed):
            var.LB = var.UB = (value >> c) & 1
        model.optimize()
        assert model.Status == gp.GRB.OPTIMAL and model.SolCount > 0
        actual = sum(round(var.X) << c for c, var in enumerate(other))
        assert actual == apply(mapping, value), (path, value)
        for end, state in ((fixed_side, value), (other_side, actual)):
            for j in range(64):
                activity = model.getVarByName(f"active_{end}_{j}")
                assert round(activity.X) == (columns(state) >> j) & 1
    # A wrong output bit and a wrong activity flag must both be rejected.
    expected_bit = (apply(mapping, samples[-1]) & 1)
    other[0].LB = other[0].UB = 1 - expected_bit
    model.optimize()
    assert model.Status == gp.GRB.INFEASIBLE
    other[0].LB, other[0].UB = 0, 1
    activity = model.getVarByName(f"active_{fixed_side}_0")
    activity.LB = activity.UB = 1 - (columns(samples[-1]) & 1)
    model.optimize()
    assert model.Status == gp.GRB.INFEASIBLE
    model.dispose()


def main():
    forward, transpose = reference_matrix()
    inv, inv_transpose = inverse(forward), inverse(transpose)
    for c in range(BITS):
        assert apply(inv, apply(forward, 1 << c)) == 1 << c
        assert apply(inv_transpose, apply(transpose, 1 << c)) == 1 << c
    maps = {
        ("difference", "input"): (forward, 18),
        ("difference", "output"): (inv, 57),
        ("linear", "input"): (inv_transpose, 59),
        ("linear", "output"): (transpose, 17),
    }
    with tempfile.TemporaryDirectory(prefix="gaston_milp_check_") as tmp:
        for (mode, side), (mapping, expected) in maps.items():
            # Shift invariance reduces every one-active-column state to one
            # of these 31 values at column zero, giving an exact reference optimum.
            candidates = [sum(((mask >> i) & 1) << (64 * i) for i in range(5))
                          for mask in range(1, 32)]
            assert min(weight(apply(mapping, value)) for value in candidates) == expected
            for count, objective in ((0, 0), (1, expected)):
                prefix = Path(tmp) / f"{mode}_{side}_{count}"
                subprocess.run([str(PROGRAM), "--mode", mode, f"--{side}-active", str(count),
                                "--threads", "2", "--prefix", str(prefix)], check=True)
                check_json(prefix.with_suffix(".json"), mode, side, count, objective, forward, transpose)
                if count == 1:
                    check_exported_model(prefix.with_suffix(".lp"), side, mapping)
        prefix = Path(tmp) / "export_only"
        subprocess.run([str(PROGRAM), "--mode", "linear", "--input-active", "64",
                        "--no-symmetry", "--export-only", "--prefix", str(prefix)], check=True)
        assert prefix.with_suffix(".lp").is_file()
        assert not prefix.with_suffix(".json").exists()
        check_effect_models(Path(tmp), forward, transpose)
        check_nonzero_default(Path(tmp), forward, transpose, inv_transpose)
    bad_args = [["--input-active", "65"], ["--output-active", "-1"],
                ["--input-active", "2", "--output-active", "3"],
                ["--mode", "unknown", "--input-active", "1"],
                ["--input-active", "1junk"], ["--input-active"],
                ["--input-active", "1", "--export-only"]]
    bad_args += [["--input-active", "1", "--theta-effect", "invalid"],
                 ["--input-active", "1", "--theta-effect"],
                 ["--input-active", "1", "--theta-effect", "zero", "--theta-effect", "nonzero"]]
    for policy in ("any", "zero", "nonzero"):
        bad_args.append(["--mode", "linear", "--input-active", "1", "--theta-effect", policy])
        bad_args.append(["--theta-effect", policy, "--mode", "linear", "--output-active", "1"])
    for args in bad_args:
        result = subprocess.run([str(PROGRAM), *args], capture_output=True, text=True)
        assert result.returncode != 0 and "错误:" in result.stderr
    print("PASS: all-binary models and complete Jeroslow XOR hulls; shared inverse expressions; "
          "320-bit inverses; 4 exhaustive 31-state optima; 8 CLI optimizations; "
          "24 propagation/activity cases; 8 infeasibility checks; theta-effect checks; "
          "optional-count and first-column symmetry checks; "
          f"export-only; {len(bad_args)} argument checks.", flush=True)


if __name__ == "__main__":
    main()
