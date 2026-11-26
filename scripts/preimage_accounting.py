#!/usr/bin/env python3
"""Recalculate the three preimage results under Huang's elimination accounting.

This implements the accounting convention used in the submitted manuscript.
Coefficient construction, coordinate reconstruction and candidate enumeration
are not charged in this model. Numerical verification inputs are hashed for provenance.
"""

import hashlib
import json
import math
from pathlib import Path


def E(n, k):
    return sum(i * n - i * (i - 1) // 2 for i in range(1, k))


def result(search_exponent, permutation_bits, recurring_bits, setup_bits=0,
           setups_shared_by=1):
    ts = (recurring_bits + setup_bits / setups_shared_by) / permutation_bits
    time_exponent = search_exponent + math.log2(1 + ts + 2**(64-search_exponent))
    return {
        "search_exponent": search_exponent,
        "permutation_bit_operations": permutation_bits,
        "recurring_elimination_bit_operations_per_candidate": recurring_bits,
        "setup_elimination_bit_operations": setup_bits,
        "candidates_per_setup": setups_shared_by,
        "Ts_elimination_only": ts,
        "Tv": 1,
        "exponent_retaining_1_plus_Ts": search_exponent + math.log2(1 + ts),
        "second_word_tests_exponent": 64,
        "time_formula": "C_s * (1 + T_s) + 2^64",
        "reported_time_exponent": time_exponent,
        "reported_time_exponent_rounded": f"{time_exponent:.2f}",
    }


def main():
    # Reproduce Huang's printed fixed-row and variable-row counts.
    assert E(81, 57) == 100016
    assert E(81, 68) - E(81, 57) == 34386
    assert E(107, 90) == 311055
    assert E(107, 93) - E(107, 90) == 16925
    for n in (2, 17, 22):
        assert E(n, n) == (n**3 - n) // 3

    eta = 8
    g_exponent = 111 - 17 * math.log2(3 / 4)
    g = result(g_exponent, 4096, E(17, 17), E(64, 47), 2**47)
    h2 = result(106, 3712, E(22, 22))
    h3 = result(126, 5568,
                (E(62, 59) - E(62, 55)) / eta + E(2, 2),
                E(62, 55), eta * 2**23)
    root = Path(__file__).resolve().parents[1] / "research/Gaston_Hsilu_CTP"
    verification_inputs = root / "verification_inputs"
    output = {
        "date": "2026-11-26",
        "scope": "Current manuscript Sections 5.1--5.3 and Appendix F; retain solving and verification costs under Huang's precomputation and Gaussian-elimination accounting",
        "warning": (
            "These estimates exclude coefficient formation and candidate handling. "
            "The estimate follows the manuscript's stated accounting model. "
            "The reported time is Cs*(1+Ts)+2^64 within this model; "
            "prefix selection continues to use Huang's max allocation convention."
        ),
        "assumptions": {
            "Gaston_conditional_success": "(3/4)^17",
            "Hsilu2_prefix_density": "2^-41",
            "Hsilu3_mean_backgrounds_per_branch": eta,
            "Hsilu3_prefix_selection": "p_sel * mean_backgrounds_per_prefix >= 1",
            "target_coverage": "uniform target-coverage model of the manuscript",
            "second_word_tests": "2^64, negligible at displayed precision",
        },
        "Gaston_2round": g,
        "Hsilu_2round": h2,
        "Hsilu_3round": h3,
        "Hsilu3_fixed_row_count": E(62, 55),
        "Hsilu3_variable_row_count": E(62, 59) - E(62, 55),
        "Hsilu3_concrete_prefix_eta": 8 + 2**-11,
        "verification_inputs_sha256": {
            name: hashlib.sha256((verification_inputs / name).read_bytes()).hexdigest()
            for name in ("paper_parameters.tex",)
        },
    }
    for name, expected in [("Gaston_2round", "118.54"), ("Hsilu_2round", "106.97"),
                           ("Hsilu_3round", "126.23")]:
        assert output[name]["reported_time_exponent_rounded"] == expected
        assert output[name]["reported_time_exponent"] > output[name]["search_exponent"]
    destination = Path(__file__).resolve().parents[1] / "CURRENT_PREIMAGE_RESULTS.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({name: output[name] for name in
                     ("Gaston_2round", "Hsilu_2round", "Hsilu_3round")}, indent=2))


if __name__ == "__main__":
    main()
