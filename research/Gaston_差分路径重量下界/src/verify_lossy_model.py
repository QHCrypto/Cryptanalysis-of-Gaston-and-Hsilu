#!/usr/bin/env python3
"""Audit the actual exported MPS: parity facets, official C propagation, lossy DDT."""
import argparse
import ctypes
import hashlib
import itertools
import json
import random
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import gurobipy as gp

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT.parent / "gaston_ref/gaston/implementation/reference"
COST = (0, 2, 3, 3, 4, 4)
MASK = (1 << 64) - 1


def chi(x):
    return sum((((x >> i) & 1) ^
                ((1 ^ ((x >> ((i + 1) % 5)) & 1)) &
                 ((x >> ((i + 2) % 5)) & 1))) << i for i in range(5))


def column(words, z):
    return sum(((w >> z) & 1) << i for i, w in enumerate(words))


def rol(x, n):
    n %= 64
    return ((x << n) | (x >> ((64 - n) & 63))) & MASK


def rank(values):
    basis = {}
    for v in values:
        while v:
            k = v.bit_length() - 1
            if k in basis:
                v ^= basis[k]
            else:
                basis[k] = v
                break
    return len(basis)


def holds(row, values):
    coefficients, sense, rhs = row
    lhs = sum(c * values[v] for v, c in coefficients.items())
    return lhs <= rhs if sense == "<" else lhs >= rhs if sense == ">" else lhs == rhs


def verify_profiles(rows, used, settings):
    profile_path = ROOT / "references/activity_profiles.json"
    assert hashlib.sha256(profile_path.read_bytes()).hexdigest() == settings["activity_profile_sha256"]
    p = json.loads(profile_path.read_text())
    for record in p["records"] + p["branch_sources"]:
        assert hashlib.sha256(Path(record["path"]).read_bytes()).hexdigest() == record["sha256"]
    forward, reverse = p["forward_lower_by_exact_input_count"], p["reverse_lower_by_exact_output_count"]
    def add(*terms):
        out = Counter()
        for scale, term in terms:
            for name, coefficient in term.items():
                out[name] += scale * coefficient
        return {k: v for k, v in out.items() if v}
    def check(name, coefficients, sense, rhs):
        assert rows[name] == (coefficients, sense, rhs), name
        used.add(name)
    counts = [{f"chi{r}_active_{z}": 1 for z in range(64)} for r in range(3)]
    selectors = [{f"count{r}_{z}": 1 for z in range(64)} for r in range(3)]
    total = add(*[(1, c) for c in counts])
    for r in range(3):
        check(f"profile_onehot_{r}", selectors[r], "=", 1)
        check(f"profile_count_{r}",
              add((1, counts[r]), (-1, {f"count{r}_{z}": z+1 for z in range(64)})), "=", 0)
    check("profile_three_round_total", total, ">", 18)
    middle_total = {f"count1_{z}": p["total_count_lower_by_middle_count"][z+1] for z in range(64)}
    check("profile_middle_total", add((1, total), (-1, middle_total)), ">", 0)
    for r in range(2):
        prefix = f"profile_L{r}"
        f = {f"count{r}_{z}": forward[z+1] for z in range(64)}
        b = {f"count{r+1}_{z}": reverse[z+1] for z in range(64)}
        if settings.get("kernel_reduced", False):
            check(f"{prefix}_forward", add((1, counts[r+1]), (-1, f)), ">", 0)
            check(f"{prefix}_reverse", add((1, counts[r]), (-1, b)), ">", 0)
            check(f"{prefix}_branch", add((1, counts[r]), (1, counts[r+1])), ">", 12)
            check(f"{prefix}_kernel_forward", add((1, counts[r+1]), (-5, counts[r])), "<", 0)
            check(f"{prefix}_kernel_reverse", add((1, counts[r]), (-5, counts[r+1])), "<", 0)
            continue
        n = {f"L{r}_nonzero": 1}
        effect = {f"L{r}_e_{z}": 1 for z in range(64)}
        for z in range(64):
            check(f"{prefix}_effect_bit{z}", add((1, {f"L{r}_e_{z}": 1}), (-1, n)), "<", 0)
        check(f"{prefix}_effect_or", add((1, n), (-1, effect)), "<", 0)
        check(f"{prefix}_effect_even_min", add((1, effect), (-2, n)), ">", 0)
        f = {f"count{r}_{z}": forward[z+1] for z in range(64)}
        b = {f"count{r+1}_{z}": reverse[z+1] for z in range(64)}
        check(f"{prefix}_forward", add((1, counts[r+1]), (-1, f)), ">", 0)
        check(f"{prefix}_reverse", add((1, counts[r]), (-1, b)), ">", 0)
        check(f"{prefix}_branch", add((1, counts[r]), (1, counts[r+1]), (-4, n)), ">", 12)
        check(f"{prefix}_outside_output6",
              add((1, counts[r]), (-1, b), (-7, {f"count{r+1}_5": 1}), (-7, n)), ">", -7)
        check(f"{prefix}_kernel_forward", add((1, counts[r+1]), (-5, counts[r]), (-59, n)), "<", 0)
        check(f"{prefix}_kernel_reverse", add((1, counts[r]), (-5, counts[r+1]), (-59, n)), "<", 0)
        check(f"{prefix}_three_round_total", add((1, total), (-4, n)), ">", 18)
    # Independently join the two 64x64 count relations and minimize total count.
    edges = {}
    for nonzero in (0, 1):
        edges[nonzero] = [(a, b) for a in range(1, 65) for b in range(1, 65)
                          if b >= forward[a] and a >= reverse[b] and a+b >= 12+4*nonzero
                          and (not nonzero or b != 6 or a >= 13)
                          and (nonzero or (a <= 5*b and b <= 5*a))]
    minima = {}
    for n0, n1 in itertools.product((0, 1), repeat=2):
        incoming = {b: min(a for a, y in edges[n0] if y == b)
                    for b in {y for _, y in edges[n0]}}
        value = min(incoming[b]+b+c for b, c in edges[n1] if b in incoming)
        minima[f"{n0}{n1}"] = value
    assert minima == {"00": 18, "01": 22, "10": 22, "11": 22}
    raw_f, raw_r = p["raw_forward_lower_by_exact_input_count"], p["raw_reverse_lower_by_exact_output_count"]
    pairs = {(a, b) for a in range(1, 65) for b in range(1, 65)
             if b >= raw_f[a] and a >= raw_r[b] and a+b >= 12}
    assert pairs == {(a, b) for a in range(1, 65) for b in range(1, 65)
                     if b >= forward[a] and a >= reverse[b] and a+b >= 12}
    assert p["total_count_lower_by_middle_count"] == [0] + [
        k + min(a for a, b in pairs if b == k) + min(b for a, b in pairs if a == k)
        for k in range(1, 65)]
    return minima


def verify_kernel_types(rows, used, settings, variables, state_variables):
    audit_path = ROOT / "results/kernel_structure_audit_20261024.json"
    assert hashlib.sha256(audit_path.read_bytes()).hexdigest() == settings["kernel_audit_sha256"]
    audit = json.loads(audit_path.read_text())
    profile = json.loads((ROOT / "references/activity_profiles.json").read_text())
    f, b = profile["forward_lower_by_exact_input_count"], profile["reverse_lower_by_exact_output_count"]
    # Recompute the complete necessary-type list independently from the generated header.
    all_types = []
    minima = {}
    for k in (6, 7, 8, 9):
        possible, minimum = [], 1000
        for c0, c2 in itertools.product(range(1, 65), repeat=2):
            if any(not (y >= f[x] and x >= b[y] and x+y >= 12 and x <= 5*y and y <= 5*x)
                   for x, y in ((c0, k), (k, c2))):
                continue
            lo1, lo2 = 2*((max(c0,k)+1)//2), 2*((max(k,c2)+1)//2)
            floor = 2*(c0+k+c2)+(lo1-k+1)//2+(lo2-c2+1)//2
            minimum = min(minimum, floor)
            if floor > 43:
                continue
            for h1 in range(lo1, 5*min(c0,k)+1, 2):
                for h2 in range(lo2, 5*min(k,c2)+1, 2):
                    floor = 2*(c0+k+c2)+(h1-k+1)//2+(h2-c2+1)//2
                    if floor <= 43:
                        possible.append({"counts": [c0,k,c2], "internal_bit_weights": [h1,h2], "weight_lower_bound": floor})
        assert possible == audit["branch_results"][str(k)]["remaining_types_under_44"]
        minima[k] = minimum
        all_types.extend(possible)
    assert minima == {6:36,7:42,8:43,9:45}
    selected = [(f"kernel_type_{j}", item) for j,item in enumerate(all_types) if item["counts"][1] == settings["middle_active"]]
    assert {v for v in variables if v.startswith("kernel_type_")} == {v for v,_ in selected}
    def check(name, expr, sense, rhs):
        assert rows[name] == ({v:c for v,c in expr.items() if c},sense,rhs), name
        used.add(name)
    check("kernel_type_onehot", {v:1 for v,_ in selected}, "=", 1)
    for r in range(3):
        coeff = {f"chi{r}_active_{z}":1 for z in range(64)}
        coeff.update({v:-item["counts"][r] for v,item in selected})
        check(f"kernel_type_count{r}", coeff, "=", 0)
    for r in range(2):
        coeff = {state_variables[2*(r+1)][i][z]:1 for i in range(5) for z in range(64)}
        coeff.update({v:-item["internal_bit_weights"][r] for v,item in selected})
        check(f"kernel_type_bits{r}", coeff, "=", 0)
        coeff = {f"chi{r+1}_extra_{z}":1 for z in range(64)}
        for v,item in selected:
            coeff[v]=-((item["internal_bit_weights"][r]-item["counts"][r+1]+1)//2)
        check(f"kernel_type_extra{r}", coeff, ">", 0)
    coeff = {v.VarName:v.Obj for v in variables.values() if v.Obj}
    coeff.update({v:-item["weight_lower_bound"] for v,item in selected})
    check("kernel_type_weight_floor", coeff, ">", 0)
    return {"type_count": len(selected), "analytically_excluded_middle9_lower_bound": minima[9], "all_type_lists_recomputed": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    directory = args.directory.resolve()
    metadata = json.loads((directory / "encoding.json").read_text())
    settings = json.loads((directory / "settings.json").read_text())
    reduced = settings.get("kernel_reduced", False)
    east,west = (0,60,22,27,4),(0,56,31,46,43)
    labels = ("a0","b0","a1","b1","a2","b2")
    state_variables = [[[f"{labels[t]}_{i}_{z}" for z in range(64)] for i in range(5)] for t in range(6)]
    if reduced:
        assert settings["kernel_both"] and settings["middle_active"] in (6,7,8) and settings["weight_cap"] == 43
        for t in (2,4):
            state_variables[t] = [[f"b{t//2-1}_{i}_{(z-east[i]-west[i])%64}" for z in range(64)] for i in range(5)]
    if "state_variables" in metadata:
        assert metadata["state_variables"] == state_variables
    model = gp.read(str(directory / "model.mps"))
    variables = {v.VarName: v for v in model.getVars()}
    rows = {}
    for constraint in model.getConstrs():
        expression = model.getRow(constraint)
        rows[constraint.ConstrName] = (
            {expression.getVar(i).VarName: expression.getCoeff(i)
             for i in range(expression.size())},
            constraint.Sense, constraint.RHS)
    assert len(rows) == metadata["constraint_count"]
    assert len(variables) == metadata["variable_count"]
    assert model.ModelSense == gp.GRB.MINIMIZE
    used = set()
    for name, variable in variables.items():
        if "_extra_" in name:
            assert variable.VType == gp.GRB.INTEGER and variable.LB == 0 and variable.UB == 2
            assert variable.Obj == 1
        else:
            assert variable.VType == gp.GRB.BINARY and variable.LB == 0 and variable.UB == 1
            assert variable.Obj == (2 if "_active_" in name else 0)

    # Prove each actual exported gate equals XOR by its complete local truth table.
    gate_outputs = set()
    for gate in metadata["gates"]:
        names = gate["inputs"] + [gate["output"]]
        assert len(set(names)) == len(names) and len(names) in (3, 4)
        assert gate["output"] not in gate_outputs
        gate_outputs.add(gate["output"])
        constraints = [gate["name"] + "_f" + str(mask)
                       for mask in range(1 << len(names)) if mask.bit_count() % 2]
        used.update(constraints)
        local_rows = [rows[name] for name in constraints]
        assert all(set(row[0]) <= set(names) for row in local_rows)
        for values in itertools.product((0, 1), repeat=len(names)):
            accepted = all(holds(row, dict(zip(names, values))) for row in local_rows)
            assert accepted == (sum(values) % 2 == 0), gate["name"]
    if reduced:
        assert len(metadata["gates"]) == 512 and len(metadata["kernel_checks"]) == 128
        for r in range(2):
            for z in range(64):
                check = metadata["kernel_checks"][64*r+z]
                names = [f"L{r}_p_{z}",f"L{r}_p_{(z-1)%64}",f"L{r}_q_{z}",f"L{r}_q_{(z-18)%64}"]
                assert check == {"name":f"L{r}_kernel_{z}","bits":names}
                local = []
                for mask in range(16):
                    if mask.bit_count()%2:
                        name=f"L{r}_kernel_{z}_f{mask}"
                        assert rows[name] == ({v:1 if mask>>j&1 else -1 for j,v in enumerate(names)},"<",mask.bit_count()-1)
                        used.add(name); local.append(rows[name])
                for assignment in itertools.product((0,1),repeat=4):
                    assert all(holds(row,dict(zip(names,assignment))) for row in local) == (sum(assignment)%2==0)
    else:
        assert not metadata.get("kernel_checks")

    # All 192 chi columns must have exactly the intended four support + two cost rows.
    for r in range(3):
        for z in range(64):
            x = {state_variables[2*r][i][z]: 1 for i in range(5)}
            y = {state_variables[2*r+1][i][z]: 1 for i in range(5)}
            active, extra = f"chi{r}_active_{z}", f"chi{r}_extra_{z}"
            expected = {
                "input_min": ({**x, active: -1}, ">", 0),
                "input_max": ({**x, active: -5}, "<", 0),
                "output_min": ({**y, active: -1}, ">", 0),
                "output_max": ({**y, active: -5}, "<", 0),
                "cost_min": ({**x, active: -1, extra: -2}, "<", 0),
                "cost_max": ({extra: 1, active: -2}, "<", 0),
            }
            for suffix, target in expected.items():
                name = f"chi{r}_{z}_{suffix}"
                assert rows[name] == target, name
                used.add(name)
            if settings.get("or_hull", False):
                for side, letter in (("input", "a"), ("output", "b")):
                    for i in range(5):
                        name = f"chi{r}_{z}_{side}_bit{i}"
                        bit = state_variables[2*r+(side=="output")][i][z]
                        assert rows[name] == ({bit: 1, active: -1}, "<", 0)
                        used.add(name)
    expected_nonzero = ({f"chi0_active_{z}": 1 for z in range(64)}, ">", 1)
    assert rows["nonzero_trail"] == expected_nonzero
    used.add("nonzero_trail")
    if settings["column_symmetry"]:
        assert rows["column_shift_representative"] == ({"chi0_active_0": 1}, "=", 1)
        used.add("column_shift_representative")
    count_minima = None
    if settings.get("extra_branch_cuts", False):
        count_minima = verify_profiles(rows, used, settings)
    kernel_type_report = None
    if reduced:
        kernel_type_report = verify_kernel_types(rows,used,settings,variables,state_variables)
    if settings.get("kernel_both", False) and not reduced:
        for r in range(2):
            for z in range(64):
                name = f"subproblem_kernel_L{r}_{z}"
                assert rows[name] == ({f"L{r}_e_{z}": 1}, "=", 0)
                used.add(name)
    if settings.get("middle_active") is not None:
        name = "subproblem_middle_active"
        assert rows[name] == ({f"chi1_active_{z}": 1 for z in range(64)}, "=", settings["middle_active"])
        used.add(name)
    if settings.get("weight_cap") is not None:
        name = "subproblem_weight_cap"
        objective = {v.VarName: v.Obj for v in variables.values() if v.Obj}
        assert rows[name] == (objective, "<", settings["weight_cap"])
        used.add(name)
    assert used == set(rows), "Unaccounted constraints in exported model"

    ddt = [[0] * 32 for _ in range(32)]
    for a in range(32):
        for x in range(32):
            ddt[a][chi(x) ^ chi(x ^ a)] += 1
    chi_rows = [row for name, row in rows.items() if name.startswith("chi0_0_")]
    assert len(chi_rows) == (16 if settings.get("or_hull", False) else 6)
    allowed = 0
    for a in range(32):
        for b in range(32):
            costs = []
            for active, extra in itertools.product(range(2), range(3)):
                values = {f"a0_{i}_0": (a >> i) & 1 for i in range(5)}
                values.update({f"b0_{i}_0": (b >> i) & 1 for i in range(5)})
                values.update(chi0_active_0=active, chi0_extra_0=extra)
                if all(holds(row, values) for row in chi_rows):
                    costs.append(2 * active + extra)
            assert bool(costs) == ((a == 0) == (b == 0))
            if costs:
                allowed += 1
                assert min(costs) == COST[a.bit_count()]
            if ddt[a][b]:
                assert costs and min(costs) <= 5 - (ddt[a][b].bit_length() - 1)

    class CState(ctypes.Structure):
        _fields_ = [("x", ctypes.c_uint64 * 5)]

    inverse_chi = {chi(x): x for x in range(32)}
    with tempfile.TemporaryDirectory(prefix="gaston-lossy-verify-") as temp:
        library_path = Path(temp) / "reference.so"
        subprocess.run(["cc", "-O2", "-shared", "-fPIC", str(REFERENCE / "gaston.c"),
                        "-o", str(library_path)], check=True)
        library = ctypes.CDLL(str(library_path))
        library.gaston_round.argtypes = [ctypes.POINTER(CState), ctypes.c_uint64]
        library.gaston_round.restype = None

        def reference_linear(words):
            s = CState((ctypes.c_uint64 * 5)(*words))
            library.gaston_round(ctypes.byref(s), 0)
            columns = [inverse_chi[column(s.x, z)] for z in range(64)]
            return [sum(((a >> i) & 1) << z for z, a in enumerate(columns))
                    for i in range(5)]

        def encoded_linear(words, r):
            values = {f"b{r}_{i}_{z}": (words[i] >> z) & 1
                      for i in range(5) for z in range(64)}
            for gate in metadata["gates"]:
                if gate["name"].startswith(f"L{r}_"):
                    values[gate["output"]] = sum(values[v] for v in gate["inputs"]) % 2
            if reduced:
                effect = sum((sum(values[v] for v in check["bits"])%2)<<z
                             for z,check in enumerate(metadata["kernel_checks"][64*r:64*r+64]))
                assert (effect&0x5555555555555555).bit_count()%2 == sum(w.bit_count() for w in words)%2
                # The aliases are the true output when effect=0. Verify the exact
                # effect map as well by restoring the deleted correction on arbitrary inputs.
                direct = [sum(values[state_variables[2*(r+1)][i][z]]<<z for z in range(64)) for i in range(5)]
                return [w ^ rol(effect,23+west[i]) for i,w in enumerate(direct)]
            assert sum(values[f"L{r}_e_{z}"] for z in range(64)) % 2 == 0
            return [sum(values[f"a{r+1}_{i}_{z}"] << z for z in range(64))
                    for i in range(5)]

        basis = []
        for bit in range(320):
            words = [0] * 5
            words[bit // 64] = 1 << (bit % 64)
            basis.append(words)
        rng = random.Random(20261023)
        cases = [[0] * 5] + basis + [[rng.getrandbits(64) for _ in range(5)] for _ in range(32)]
        for words in cases:
            expected = reference_linear(words)
            for r in range(2):
                assert encoded_linear(words, r) == expected
            assert reference_linear([rol(w, 17) for w in words]) == [rol(w, 17) for w in expected]
        linear_rank = rank([sum(w << (64*i) for i, w in enumerate(reference_linear(words)))
                            for words in basis])
        assert linear_rank == 320

    report = {
        "status": "PASS",
        "scope": "Exported lossy model encoding; not a solved three-round optimum",
        "subproblem": {k: settings.get(k) for k in ("kernel_both", "middle_active", "weight_cap")},
        "kernel_reduced": reduced,
        "kernel_type_verification": kernel_type_report,
        "kernel_even_bit_weight_identity_verified_on_full_basis": reduced,
        "variables": len(variables),
        "constraints": len(rows),
        "xor_gates_exhaustively_verified": len(metadata["gates"]),
        "chi_columns_checked": 192,
        "chi_input_output_pairs_checked": 1024,
        "lossy_supported_pairs": allowed,
        "true_ddt_supported_pairs": sum(bool(n) for row in ddt for n in row),
        "every_true_ddt_pair_embeds_at_no_greater_cost": True,
        "lossy_minimum_cost_by_input_hw": COST,
        "reference_linear_cases_per_layer": len(cases),
        "reference_linear_layers_checked": 2,
        "linear_rank": linear_rank,
        "initial_linear_layer_can_be_eliminated_for_free_input": True,
        "column_translation_symmetry_checked": True,
        "activity_profiles_enabled": settings.get("extra_branch_cuts", False),
        "activity_or_hull_enabled": settings.get("or_hull", False),
        "count_minima_by_effect_nonzero_flags": count_minima,
        "sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in [directory/"model.mps", directory/"encoding.json",
                             ROOT/"gaston_3round_lossy.cpp", ROOT/"gaston_3round_lossy",
                             REFERENCE/"gaston.c", REFERENCE/"gaston.h",
                             ROOT/"src/activity_profile_data.hpp", ROOT/"references/activity_profiles.json",
                             ROOT/"src/kernel_type_data.hpp", ROOT/"results/kernel_structure_audit_20261024.json"]},
    }
    (directory/"verification.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
