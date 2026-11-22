#!/usr/bin/env python3
"""Substitute Hsilu X2 conditions into X3 conditions and minimize guesses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import hsilu_3round_best_attack as hsilu

Var = Tuple[int, int]
Monomial = frozenset
Polynomial = Set[Monomial]


def xor_poly(left: Polynomial, right: Polynomial) -> Polynomial:
    return left.symmetric_difference(right)


def multiply_poly(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = set()
    for left_term in left:
        for right_term in right:
            term = left_term | right_term
            if term in result:
                result.remove(term)
            else:
                result.add(term)
    return result


def variable_poly(variable: Var) -> Polynomial:
    return {frozenset((variable,))}


def constant_poly(value: int) -> Polynomial:
    return {frozenset()} if value else set()


def substitute(poly: Polynomial, fixed: Mapping[Var, int]) -> Polynomial:
    result: Polynomial = set()
    for monomial in poly:
        term: Polynomial = {frozenset()}
        for variable in monomial:
            if variable in fixed:
                factor = constant_poly(fixed[variable])
            else:
                factor = variable_poly(variable)
            term = multiply_poly(term, factor)
            if not term:
                break
        result = xor_poly(result, term)
    return result


def collect_fixed_x2(data: Mapping[str, object]) -> Dict[Var, int]:
    fixed: Dict[Var, int] = {}
    for key, value in (("conditions_eq_0", 0), ("conditions_eq_1", 1)):
        states = data[key]
        for lane_text, positions in states.get("X2", {}).items():
            lane = int(lane_text)
            for position in positions:
                variable = (lane, int(position))
                if variable in fixed and fixed[variable] != value:
                    raise ValueError("Contradictory X2 condition at %r" % (variable,))
                fixed[variable] = value
    return fixed


def collect_x3_targets(data: Mapping[str, object]) -> List[Tuple[int, int, int]]:
    targets = []
    for key, value in (("conditions_eq_0", 0), ("conditions_eq_1", 1)):
        states = data[key]
        for lane_text, positions in states.get("X3", {}).items():
            lane = int(lane_text)
            targets.extend((lane, int(position), value) for position in positions)
    return sorted(targets)


def x3_condition_polynomial(
    lane: int,
    bit: int,
    value: int,
    fixed_x2: Mapping[Var, int],
) -> Polynomial:
    # Hsilu chi: Y2[i,j] = X2[i,j] + X2[i+2,j] + X2[i+1,j]X2[i+2,j].
    y2: Dict[Var, Polynomial] = {}
    for source_lane in range(hsilu.LANES):
        for source_bit in range(hsilu.LANE_BITS):
            a = variable_poly((source_lane, source_bit))
            b = variable_poly(((source_lane + 1) % hsilu.LANES, source_bit))
            c = variable_poly(((source_lane + 2) % hsilu.LANES, source_bit))
            expression = xor_poly(xor_poly(a, c), multiply_poly(b, c))
            y2[source_lane, source_bit] = substitute(expression, fixed_x2)

    output: Polynomial = constant_poly(value)
    for source_lane, source_bit in sorted(hsilu.MRM_SUPPORTS[lane, bit]):
        output = xor_poly(output, y2[source_lane, source_bit])
    return output


def format_variable(variable: Var) -> str:
    lane, bit = variable
    return "X2[%d][%02d]" % (lane, bit)


def format_monomial(monomial: Monomial) -> str:
    if not monomial:
        return "1"
    return " * ".join(format_variable(variable) for variable in sorted(monomial))


def format_polynomial(poly: Polynomial) -> str:
    if not poly:
        return "0"
    terms = sorted(poly, key=lambda item: (len(item), tuple(sorted(item))))
    return " + ".join(format_monomial(term) for term in terms)


def solve_minimum_guesses(polynomials: Sequence[Polynomial]):
    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError as exc:
        raise RuntimeError("gurobipy is required to prove the minimum.") from exc

    quadratic_terms = sorted(
        {
            tuple(sorted(monomial))
            for poly in polynomials
            for monomial in poly
            if len(monomial) == 2
        }
    )
    unsupported = sorted(
        {
            monomial
            for poly in polynomials
            for monomial in poly
            if len(monomial) > 2
        },
        key=lambda item: (len(item), tuple(sorted(item))),
    )
    if unsupported:
        raise ValueError("Unexpected degree > 2 terms: %r" % (unsupported,))

    variables = sorted({variable for edge in quadratic_terms for variable in edge})
    model = gp.Model("hsilu_x3_condition_minimum_guess")
    model.Params.OutputFlag = 0
    guess = model.addVars(variables, vtype=GRB.BINARY, name="guess_x2")
    for index, (left, right) in enumerate(quadratic_terms):
        model.addConstr(guess[left] + guess[right] >= 1, name="hit_quad_%d" % index)
    model.setObjective(gp.quicksum(guess[variable] for variable in variables), GRB.MINIMIZE)
    model.optimize()
    if model.Status != GRB.OPTIMAL:
        status = model.Status
        model.dispose()
        raise RuntimeError("Minimum-guess model did not prove optimality (status %d)." % status)
    selected = [variable for variable in variables if guess[variable].X > 0.5]
    minimum = int(round(model.ObjVal))
    model.dispose()
    return minimum, selected, quadratic_terms


def build_report(source: Path) -> str:
    data = json.loads(source.read_text(encoding="utf-8"))
    fixed_x2 = collect_fixed_x2(data)
    targets = collect_x3_targets(data)
    if not targets:
        raise ValueError("No X3 conditions in source result.")

    polynomials = [
        x3_condition_polynomial(lane, bit, value, fixed_x2)
        for lane, bit, value in targets
    ]
    minimum, selected, quadratic_terms = solve_minimum_guesses(polynomials)
    selected_set = set(selected)
    linearized_polynomials = [
        substitute(poly, {variable: 0 for variable in selected_set})
        for poly in polynomials
    ]
    if any(
        len(monomial) > 1
        for polynomial in linearized_polynomials
        for monomial in polynomial
    ):
        raise RuntimeError("Chosen guesses did not linearize every X3 condition.")

    output = [
        "Hsilu 3-round X3 condition linearization over X2",
        "Source result: %s" % source,
        "",
        "State equations:",
        "  Y2[i][j] = X2[i][j] + X2[i+2][j] + X2[i+1][j]*X2[i+2][j]",
        "  X3[lane][bit] = XOR of Y2 bits in the Hsilu MRM support.",
        "All lane indices are modulo 5; all bit indices are modulo 64.",
        "",
        "X2 conditions substituted before linearization:",
    ]
    for (lane, bit), value in sorted(fixed_x2.items()):
        output.append("  %s = %d" % (format_variable((lane, bit)), value))
    output.extend(
        [
            "",
            "Summary:",
            "  X2 condition bits substituted = %d" % len(fixed_x2),
            "  X3 conditions to linearize = %d" % len(targets),
            "  distinct residual quadratic monomials = %d" % len(quadratic_terms),
            "  minimum additional X2 guesses = %d" % minimum,
            "  minimum proven by exact binary optimization = yes",
            "",
            "X3 conditions expanded over X2 after substituting X2 condition values:",
        ]
    )
    for (lane, bit, value), polynomial in zip(targets, polynomials):
        output.append(
            "  X3[%d][%02d] = %d: %s = 0"
            % (lane, bit, value, format_polynomial(polynomial))
        )

    output.extend(["", "One minimum set of additional X2 guesses:"])
    for variable in selected:
        output.append("  %s" % format_variable(variable))
    output.extend(
        [
            "",
            "Linearized equations after the selected X2 guesses:",
        ]
    )
    for (lane, bit, value), linearized in zip(targets, linearized_polynomials):
        output.append(
            "  X3[%d][%02d] = %d: %s = 0"
            % (lane, bit, value, format_polynomial(linearized))
        )

    output.extend(
        [
            "",
            "Interpretation:",
            "  X2 condition bits are treated as already-fixed constants and are not",
            "  counted as additional guesses. Each residual quadratic monomial",
            "  requires guessing at least one of its two X2 variables. The exact",
            "  minimum is therefore the minimum vertex cover of the residual",
            "  quadratic-term graph.",
            "  The displayed guessed values are set to 0 only to show the resulting",
            "  linear form; choosing 1 instead changes coefficients/constants but",
            "  does not change whether the equation is linear.",
        ]
    )
    return "\n".join(output) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        default=Path(__file__).resolve().parent / "work" / "hsilu_3round_best_attack_1_6.3.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent
        / "work"
        / "hsilu_3round_x3_conditions_over_x2_1_6.3.txt",
    )
    args = parser.parse_args()
    report = build_report(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report, encoding="utf-8")
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
