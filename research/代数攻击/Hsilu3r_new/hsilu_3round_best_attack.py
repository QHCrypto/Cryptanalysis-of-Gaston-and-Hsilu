#!/usr/bin/env python3
"""Standalone three-round Hsilu structural linearization MILP.

Hsilu round:
    chi -> theta(MRM) -> rho.

Three-round state chain:
    X1 -> chi -> Y1 -> MRM -> X2
       -> chi -> Y2 -> MRM -> X3
       -> chi -> Y3 -> MRM -> X4

X1[0] contains 64 input variables and X1[1..4] are constant capacity
words.  Conditions are allowed on X1, X2 and X3, the input states of the
three chi layers.  The selected output bits are the structurally affine
X4[0] bits.  The resource constraint is

    free X1[0] variables + weighted conditions <= 64.

This is a structural strategy search.  It does not check equation rank,
concrete condition consistency, or a real preimage.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

try:
    import gurobipy as gp
    from gurobipy import GRB
except ImportError:
    gp = None
    GRB = None


LANES = 5
LANE_BITS = 64
INPUT_BITS = 64
OUTPUT_BITS = 64
STATES = ("X1", "Y1", "X2", "Y2", "X3", "Y3", "X4")
COND_STATES = ("X1", "X2", "X3")
OUTPUT_DIR = Path(__file__).resolve().parent / "work"
Bit = Tuple[str, int, int]
LocalBit = Tuple[int, int]

# Table 9 of ToSC2025_2_06.pdf.
N1 = {5: 0, 6: 3, 7: 1, 8: 0, 9: 5, 10: 5, 11: 7,
      12: 11, 13: 12, 14: 12, 15: 13, 16: 14, 17: 15}
N2 = {5: 4, 6: 4, 7: 2, 8: 1, 9: 2, 10: 3, 11: 6,
      12: 8, 13: 9, 14: 6, 15: 10, 16: 11, 17: 9}
M = {5: 47, 6: 51, 7: 54, 8: 30, 9: 4, 10: 55, 11: 0,
     12: 14, 13: 59, 14: 4, 15: 10, 16: 45, 17: 53}
S = (0, 57, 13, 31, 20)


def rot(bit: int, offset: int) -> int:
    return (bit + offset) % LANE_BITS


def build_mrm_supports() -> Dict[LocalBit, Set[LocalBit]]:
    words: Dict[int, List[Set[LocalBit]]] = {
        lane: [{(lane, bit)} for bit in range(LANE_BITS)]
        for lane in range(LANES)
    }
    for index in range(LANES, LANES + len(N1)):
        words[index] = [set() for _ in range(LANE_BITS)]
        for bit in range(LANE_BITS):
            words[index][bit] = words[N1[index]][bit] ^ words[N2[index]][rot(bit, M[index])]
    supports: Dict[LocalBit, Set[LocalBit]] = {}
    for lane, shift in enumerate(S):
        source_word = LANES + len(N1) - LANES + lane
        for bit in range(LANE_BITS):
            supports[lane, bit] = set(words[source_word][rot(bit, shift)])
    return supports


MRM_SUPPORTS = build_mrm_supports()


def mrm_inputs(src: str, lane: int, bit: int) -> List[Bit]:
    return [(src, source_lane, source_bit)
            for source_lane, source_bit in sorted(MRM_SUPPORTS[lane, bit])]


def add_and(model, out, inputs, name: str) -> None:
    if not inputs:
        model.addConstr(out == 1, name=name + "_empty")
        return
    for index, value in enumerate(inputs):
        model.addConstr(out <= value, name="%s_ub_%d" % (name, index))
    model.addConstr(out >= gp.quicksum(inputs) - len(inputs) + 1, name=name + "_lb")


def add_or(model, out, inputs, name: str) -> None:
    if not inputs:
        model.addConstr(out == 0, name=name + "_empty")
        return
    for index, value in enumerate(inputs):
        model.addConstr(out >= value, name="%s_lb_%d" % (name, index))
    model.addConstr(out <= gp.quicksum(inputs), name=name + "_ub")


def add_or_when_linear(model, out, inputs, linear, name: str) -> None:
    model.addConstr(out <= linear, name=name + "_linear")
    if not inputs:
        model.addConstr(out == 0, name=name + "_empty")
        return
    for index, value in enumerate(inputs):
        model.addConstr(out >= value + linear - 1, name="%s_lb_%d" % (name, index))
    model.addConstr(out <= gp.quicksum(inputs), name=name + "_ub")


def chi_flags(a, b, c, b_value=None, c_value=None):
    b_const = bool(b[0] and not b[1])
    c_const = bool(c[0] and not c[1])
    killed = b_value == 1 or c_value == 0
    product_linear = killed or (b[0] and c[0] and (b_const or c_const))
    linear = bool(a[0]) and product_linear
    participating = linear and (
        bool(a[1])
        or (bool(c[1]) and b_const and b_value != 1)
        or (bool(b[1]) and c_const and c_value != 0)
    )
    return int(linear), int(participating)


def add_chi_bit(model, lin, part, cond0, cond1, src: str, dst: str, lane: int, bit: int):
    a, b, c = [(src, (lane + k) % LANES, bit) for k in range(3)]
    tag = "%s_%d_%d" % (dst, lane, bit)
    any_const, killed, normal, product_linear, c_prop, b_prop = [
        model.addVar(vtype=GRB.BINARY, name=name + "_" + tag)
        for name in ("ANY_CONST", "KILL", "NORMAL", "PRODUCT_LIN", "C_PROP", "B_PROP")
    ]
    b_const, c_const = lin[b] - part[b], lin[c] - part[c]
    add_or(model, any_const, [b_const, c_const], "any_const_" + tag)
    add_or(model, killed, [cond1[b], cond0[c]], "kill_" + tag)
    add_and(model, normal, [lin[b], lin[c], any_const], "normal_" + tag)
    add_or(model, product_linear, [killed, normal], "product_lin_" + tag)
    add_and(model, lin[dst, lane, bit], [lin[a], product_linear], "chi_lin_" + tag)
    add_and(model, c_prop, [part[c], b_const, 1 - cond1[b]], "c_prop_" + tag)
    add_and(model, b_prop, [part[b], c_const, 1 - cond0[c]], "b_prop_" + tag)
    add_or_when_linear(
        model, part[dst, lane, bit], [part[a], c_prop, b_prop],
        lin[dst, lane, bit], "chi_part_" + tag,
    )


def validate_bits(bits: Optional[Sequence[int]], label: str) -> None:
    if bits is not None and (
        len(set(bits)) != len(bits)
        or any(type(value) is not int or not 0 <= value < LANE_BITS for value in bits)
    ):
        raise ValueError(label + " must contain distinct integer positions in [0, 63].")


def parse_list(text: Optional[str]) -> Optional[List[int]]:
    if text is None:
        return None
    value = json.loads(text)
    if not isinstance(value, list) or not all(type(item) is int for item in value):
        raise ValueError("Expected a JSON list of integer positions.")
    validate_bits(value, "positions")
    return value


def build_model(
    resource_bound: int = INPUT_BITS,
    min_conditions: int = 0,
    fixed_guessed: Optional[Sequence[int]] = None,
    fixed_selected: Optional[Sequence[int]] = None,
    condition_weights: Sequence[float] = (1.0, 1.0, 3.4),
    no_conditions: bool = False,
    min_free_vars: int = 0,
    exact_free_vars: Optional[int] = None,
    output_flag: bool = True,
) -> Dict[str, object]:
    if gp is None:
        raise RuntimeError("gurobipy is missing. Use the Gurobi-enabled PyCharm interpreter.")
    if not 0 <= resource_bound <= INPUT_BITS:
        raise ValueError("resource_bound must be between 0 and 64.")
    if len(condition_weights) != 3 or any(
        not math.isfinite(value) or value < 1 for value in condition_weights
    ):
        raise ValueError("Use three finite condition weights, each >= 1.")
    if no_conditions and min_conditions:
        raise ValueError("--no-conditions conflicts with --min-conditions.")
    if not 0 <= min_free_vars <= INPUT_BITS:
        raise ValueError("Invalid min_free_vars.")
    if exact_free_vars is not None and not min_free_vars <= exact_free_vars <= INPUT_BITS:
        raise ValueError("Invalid exact_free_vars.")
    validate_bits(fixed_guessed, "fixed_guessed")
    validate_bits(fixed_selected, "fixed_selected")

    model = gp.Model("hsilu_3round_best_attack")
    model.Params.OutputFlag = int(output_flag)
    lin = model.addVars(STATES, LANES, LANE_BITS, vtype=GRB.BINARY, name="LIN")
    part = model.addVars(STATES, LANES, LANE_BITS, vtype=GRB.BINARY, name="PART")
    cond0 = model.addVars(COND_STATES, LANES, LANE_BITS, vtype=GRB.BINARY, name="COND_EQ_0")
    cond1 = model.addVars(COND_STATES, LANES, LANE_BITS, vtype=GRB.BINARY, name="COND_EQ_1")
    guessed = model.addVars(LANE_BITS, vtype=GRB.BINARY, name="G_INPUT_CONST")
    selected = model.addVars(LANE_BITS, vtype=GRB.BINARY, name="SEL_X4")

    for key in lin:
        model.addConstr(part[key] <= lin[key], name="part_implies_lin_%s_%d_%d" % key)
    for lane in range(LANES):
        for bit in range(LANE_BITS):
            lin["X1", lane, bit].LB = 1
            if lane == 0:
                model.addConstr(part["X1", lane, bit] + guessed[bit] == 1, name="input_%d" % bit)
            else:
                part["X1", lane, bit].UB = 0

    for state in COND_STATES:
        for lane in range(LANES):
            for bit in range(LANE_BITS):
                model.addConstr(
                    cond0[state, lane, bit] + cond1[state, lane, bit]
                    <= lin[state, lane, bit] - part[state, lane, bit],
                    name="condition_on_gray_%s_%d_%d" % (state, lane, bit),
                )
                if state == "X1" and lane == 0:
                    cond0[state, lane, bit].UB = 0
                    cond1[state, lane, bit].UB = 0
                if no_conditions:
                    cond0[state, lane, bit].UB = 0
                    cond1[state, lane, bit].UB = 0

    def chi_layer(src: str, dst: str) -> None:
        for lane in range(LANES):
            for bit in range(LANE_BITS):
                add_chi_bit(model, lin, part, cond0, cond1, src, dst, lane, bit)

    def mrm_layer(src: str, dst: str) -> None:
        for lane in range(LANES):
            for bit in range(LANE_BITS):
                terms = mrm_inputs(src, lane, bit)
                add_and(model, lin[dst, lane, bit], [lin[key] for key in terms],
                        "mrm_lin_%s_%d_%d" % (dst, lane, bit))
                add_or_when_linear(
                    model, part[dst, lane, bit], [part[key] for key in terms],
                    lin[dst, lane, bit], "mrm_part_%s_%d_%d" % (dst, lane, bit),
                )

    chi_layer("X1", "Y1")
    mrm_layer("Y1", "X2")
    chi_layer("X2", "Y2")
    mrm_layer("Y2", "X3")
    chi_layer("X3", "Y3")
    mrm_layer("Y3", "X4")

    free_count = INPUT_BITS - gp.quicksum(guessed[bit] for bit in range(LANE_BITS))
    condition_counts = {
        state: gp.quicksum(
            cond0[state, lane, bit] + cond1[state, lane, bit]
            for lane in range(LANES) for bit in range(LANE_BITS)
        )
        for state in COND_STATES
    }
    weighted = gp.quicksum(
        weight * condition_counts[state]
        for state, weight in zip(COND_STATES, condition_weights)
    )
    for bit in range(LANE_BITS):
        model.addConstr(selected[bit] <= lin["X4", 0, bit], name="select_x4_%d" % bit)
        if fixed_guessed is not None:
            model.addConstr(guessed[bit] == int(bit in set(fixed_guessed)), name="fixed_guess_%d" % bit)
        if fixed_selected is not None:
            model.addConstr(selected[bit] == int(bit in set(fixed_selected)), name="fixed_selected_%d" % bit)
    selected_count = gp.quicksum(selected[bit] for bit in range(LANE_BITS))
    model.addConstr(selected_count == free_count, name="selected_equals_free")
    model.addConstr(free_count + weighted <= resource_bound, name="resource_bound")
    model.addConstr(gp.quicksum(condition_counts.values()) >= min_conditions, name="min_conditions")
    model.addConstr(free_count >= min_free_vars, name="min_free_vars")
    if exact_free_vars is not None:
        model.addConstr(free_count == exact_free_vars, name="exact_free_vars")
    model.ModelSense = GRB.MINIMIZE
    model.setObjectiveN(OUTPUT_BITS - selected_count, 0, priority=2, abstol=1e-9, reltol=0,
                        name="complexity")
    model.setObjectiveN(weighted, 1, priority=1, abstol=1e-9, reltol=0, name="conditions")
    model.update()
    return dict(model=model, lin=lin, part=part, cond0=cond0, cond1=cond1,
                guessed=guessed, selected=selected, condition_weights=tuple(condition_weights),
                resource_bound=resource_bound)


def replay_flags(guessed: Iterable[int], conditions: Mapping[Bit, int]):
    guessed = set(guessed)
    flags = {("X1", lane, bit): (1, int(lane == 0 and bit not in guessed))
             for lane in range(LANES) for bit in range(LANE_BITS)}

    def chi(src: str, dst: str):
        for lane in range(LANES):
            for bit in range(LANE_BITS):
                a, b, c = [(src, (lane + k) % LANES, bit) for k in range(3)]
                flags[dst, lane, bit] = chi_flags(
                    flags[a], flags[b], flags[c], conditions.get(b), conditions.get(c)
                )

    def mrm(src: str, dst: str):
        for lane in range(LANES):
            for bit in range(LANE_BITS):
                terms = mrm_inputs(src, lane, bit)
                linear = all(flags[key][0] for key in terms)
                flags[dst, lane, bit] = (int(linear), int(linear and any(flags[key][1] for key in terms)))

    chi("X1", "Y1")
    mrm("Y1", "X2")
    chi("X2", "Y2")
    mrm("Y2", "X3")
    chi("X3", "Y3")
    mrm("Y3", "X4")
    return flags


def bit(var) -> int:
    return int(var.X > 0.5)


def by_lane(bits):
    result = {}
    for lane, value in sorted(bits):
        result.setdefault(lane, []).append(value)
    return result


def format_state_bits(bits):
    result = {}
    for state, lane, bit_position in sorted(bits):
        result.setdefault(state, {}).setdefault(lane, []).append(bit_position)
    return result


def collect_solution(data):
    model = data["model"]
    if model.SolCount == 0:
        raise ValueError("No incumbent solution.")
    lin, part = data["lin"], data["part"]
    cond0, cond1 = data["cond0"], data["cond1"]
    guessed = [bit_position for bit_position in range(LANE_BITS) if bit(data["guessed"][bit_position])]
    free = [bit_position for bit_position in range(LANE_BITS) if bit_position not in guessed]
    selected = [bit_position for bit_position in range(LANE_BITS) if bit(data["selected"][bit_position])]
    conditions0 = [key for key, value in cond0.items() if bit(value)]
    conditions1 = [key for key, value in cond1.items() if bit(value)]
    conditions = {key: 0 for key in conditions0}
    conditions.update({key: 1 for key in conditions1})
    flags = replay_flags(guessed, conditions)
    for key, expected in flags.items():
        actual = bit(lin[key]), bit(part[key])
        if actual != expected:
            raise RuntimeError("Propagation mismatch at %s" % (key,))
    if len(selected) != len(free):
        raise RuntimeError("Selected output/free-variable mismatch.")
    if any(not flags["X4", 0, bit_position][0] for bit_position in selected):
        raise RuntimeError("Selected X4 bit is not linear.")
    counts = {state: sum(key[0] == state for key in conditions) for state in COND_STATES}
    weighted = sum(weight * counts[state] for state, weight in zip(COND_STATES, data["condition_weights"]))
    if len(free) + weighted > data["resource_bound"] + 1e-6:
        raise RuntimeError("Resource bound exceeded.")
    return {
        "algorithm": "Hsilu", "rounds": 3,
        "state_order": "X1 -> chi -> Y1 -> MRM -> X2 -> chi -> Y2 -> MRM -> X3 -> chi -> Y3 -> MRM -> X4",
        "status": model.Status, "solver_optimal": model.Status == GRB.OPTIMAL,
        "runtime_seconds": model.Runtime, "solution_count": model.SolCount,
        "complexity_exponent": OUTPUT_BITS - len(selected), "objective": OUTPUT_BITS - len(selected),
        "input_bits": INPUT_BITS, "output_bits": OUTPUT_BITS,
        "guessed_count": len(guessed), "free_variable_count": len(free),
        "condition_count": len(conditions), "condition_counts": counts,
        "condition_weights": dict(zip(COND_STATES, data["condition_weights"])),
        "weighted_condition_count": weighted, "resource_usage": len(free) + weighted,
        "resource_bound": data["resource_bound"], "selected_count": len(selected),
        "guessed_by_lane": by_lane((0, value) for value in guessed),
        "free_variables_by_lane": by_lane((0, value) for value in free),
        "selected_outputs_by_lane": by_lane((0, value) for value in selected),
        "final_linearizable_x4_count": sum(flags["X4", 0, bit_position][0] for bit_position in range(LANE_BITS)),
        "conditions_eq_0": format_state_bits(conditions0),
        "conditions_eq_1": format_state_bits(conditions1),
        "mrm_parameters": {"n1": N1, "n2": N2, "m": M, "s": list(S)},
        "validation": {"structural_replay": "passed", "equation_rank": "not_checked",
                       "concrete_condition_consistency": "not_checked",
                       "real_preimage": "not_computed"},
    }


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_svg(path: Path, data, result) -> None:
    lin, part = data["lin"], data["part"]
    cond0, cond1, selected = data["cond0"], data["cond1"], data["selected"]
    cell, left = 8, 86
    panel, top, gap = LANES * cell, 78, 38
    width = left + LANE_BITS * cell + 30
    height = top + len(STATES) * (panel + gap) + 40
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
        % (width, height, width, height),
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Consolas,monospace;font-size:10px;fill:#222}.title{font-size:15px;font-weight:bold}</style>',
        '<text x="20" y="24" class="title">Hsilu 3-round best attack</text>',
        '<rect x="20" y="40" width="10" height="10" fill="#ffd84d"/><text x="35" y="49">linear variable</text>',
        '<rect x="160" y="40" width="10" height="10" fill="#b8b8b8"/><text x="175" y="49">constant</text>',
        '<rect x="270" y="40" width="10" height="10" fill="#d73a49"/><text x="285" y="49">condition = 0</text>',
        '<rect x="410" y="40" width="10" height="10" fill="#2da44e"/><text x="425" y="49">condition = 1</text>',
        '<rect x="550" y="40" width="10" height="10" fill="#8ec5ff"/><text x="565" y="49">selected X4</text>',
    ]
    for state_index, state in enumerate(STATES):
        y0 = top + state_index * (panel + gap)
        lines.append('<text x="20" y="%d" class="title">%s</text>' % (y0 - 12, state))
        for lane in range(LANES):
            row = y0 + lane * cell
            lines.append('<text x="54" y="%d">x%d</text>' % (row + 7, lane))
            for bit_position in range(LANE_BITS):
                key = state, lane, bit_position
                color = "#ffffff" if not bit(lin[key]) else ("#ffd84d" if bit(part[key]) else "#b8b8b8")
                if state == "X4" and lane == 0 and bit(selected[bit_position]):
                    color = "#8ec5ff"
                if state in COND_STATES:
                    if bit(cond0[key]):
                        color = "#d73a49"
                    elif bit(cond1[key]):
                        color = "#2da44e"
                lines.append('<rect x="%d" y="%d" width="7" height="7" fill="%s" stroke="#555" stroke-width="0.25"/>' %
                             (left + bit_position * cell, row, color))
    lines.append('<text x="20" y="%d">complexity %d; free %d; conditions %d</text>' %
                 (height - 12, result["complexity_exponent"], result["free_variable_count"], result["condition_count"]))
    lines.append("</svg>")
    write_text(path, "\n".join(lines))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Standalone three-round Hsilu best-attack MILP.")
    parser.add_argument("--condition-weights", nargs=3, type=float, default=(1.0, 1.0, 5.1),
                        metavar=("W_X1", "W_X2", "W_X3"))
    parser.add_argument("--resource-bound", type=int, default=64)
    parser.add_argument("--no-conditions", action="store_true")
    parser.add_argument("--min-conditions", type=int, default=0)
    parser.add_argument("--min-free-vars", type=int, default=0)
    parser.add_argument("--exact-free-vars", type=int)
    parser.add_argument("--fixed-guessed", type=parse_list)
    parser.add_argument("--fixed-selected", type=parse_list)
    parser.add_argument("--time-limit", type=float)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--mip-focus", type=int, choices=range(4), default=2)
    parser.add_argument("--mip-gap", type=float, default=0.0)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--log-file", type=Path)
    parser.add_argument("--json-output", type=Path, default=OUTPUT_DIR / "hsilu_3round_best_attack_1_5.1.json")
    parser.add_argument("--txt-output", type=Path, default=OUTPUT_DIR / "hsilu_3round_best_attack_1_5.1.txt")
    parser.add_argument("--svg-output", type=Path, default=OUTPUT_DIR / "hsilu_3round_best_attack_1_5.1.svg")
    args = parser.parse_args(argv)
    if gp is None:
        parser.error("gurobipy is missing. Select the Gurobi-enabled interpreter.")
    if args.time_limit is not None and (not math.isfinite(args.time_limit) or args.time_limit <= 0):
        parser.error("--time-limit must be positive.")
    model = None
    try:
        data = build_model(resource_bound=args.resource_bound, min_conditions=args.min_conditions,
                           fixed_guessed=args.fixed_guessed, fixed_selected=args.fixed_selected,
                           condition_weights=args.condition_weights, no_conditions=args.no_conditions,
                           min_free_vars=args.min_free_vars, exact_free_vars=args.exact_free_vars,
                           output_flag=not args.quiet)
        model = data["model"]
        for name, value in (("Threads", args.threads), ("MIPFocus", args.mip_focus),
                            ("Presolve", 2), ("Cuts", 2), ("Heuristics", 0.1),
                            ("MIPGap", args.mip_gap), ("Seed", args.seed),
                            ("TimeLimit", args.time_limit)):
            if value is not None:
                model.setParam(name, value)
        if args.log_file:
            args.log_file.parent.mkdir(parents=True, exist_ok=True)
            model.Params.LogFile = str(args.log_file)
        model.optimize()
        if not model.SolCount:
            print("No feasible incumbent.")
            return 2
        result = collect_solution(data)
        report = "\n".join("%s: %s" % (key, value) for key, value in result.items())
        write_text(args.json_output, json.dumps(result, indent=2) + "\n")
        write_text(args.txt_output, report + "\n")
        write_svg(args.svg_output, data, result)
        print(report)
        print("json output:", args.json_output.resolve())
        print("txt output:", args.txt_output.resolve())
        print("svg output:", args.svg_output.resolve())
        return 0
    except (ValueError, RuntimeError, gp.GurobiError) as exc:
        parser.exit(2, "Hsilu 3-round search failed: %s\n" % exc)
    finally:
        if model is not None:
            model.dispose()


if __name__ == "__main__":
    raise SystemExit(main())
