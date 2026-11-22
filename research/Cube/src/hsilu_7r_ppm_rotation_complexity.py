"""Operation counts for the previously discussed 7-round Hsilu PPM schedule.

This is a cost calculation, not a 7-round signature reconstruction or a
certificate of key recovery. The schedule interpolates each coefficient
separately, stores dense coefficient tables for rotation representatives,
and groups partial products by orbit length before rotation symmetrization.
"""

import json
import math
from pathlib import Path


def phi(n):
    return sum(math.gcd(i, n) == 1 for i in range(1, n + 1))


def count(value, unit):
    return {"value": value, "log2": math.log2(value), "unit": unit}


def main():
    subsets = math.comb(64, 32)
    divisors = [d for d in range(1, 33) if 32 % d == 0]
    burnside_numerator = sum(phi(d) * math.comb(64 // d, 32 // d)
                             for d in divisors)
    assert burnside_numerator % 64 == 0
    representatives = burnside_numerator // 64
    interpolation_points = 1 << 32
    prefix_calls = representatives * interpolation_points * (1 << 32)
    pair_operations = 5 * representatives * (1 << 64)
    memory_bits = 320 * representatives * (1 << 32)
    ascon_points = sum(math.comb(32, d) for d in range(16))
    prefix_equivalent = prefix_calls * 6 / 7
    # Analogous to the paper's XOR-only conversion, with AND, NOT,
    # rotations, indexing and data-movement overhead omitted.
    hsilu_xors_per_seven_rounds = 7 * (5 + 13)
    eta = 1 / hsilu_xors_per_seven_rounds
    result = {
        "status": "Conditional cost estimate for a specified schedule; not a verified full key-recovery result",
        "model": "Hsilu initialization prefix (IV,K0,K1,N0,0), 7 rounds, full 64-dimensional N0 cube, first 128 output bits observable",
        "assumptions": [
            "The final 128-bit cube signature leaves sufficiently few K1 candidates; constant candidate count is used for the online estimate.",
            "Uniform intermediate coefficient key-degree bound is 32; each 32-variable coefficient is recovered with 2^32 key settings and a 2^32-point nonce cube per setting.",
            "Rotation representatives are used for coefficient recovery AND product accumulation, with actual orbit lengths respected.",
            "All representative coefficient ANFs are stored densely; this is not a minimum-memory claim.",
            "Main time model uses one XOR-equivalent unit per 64-bit AND-XOR accumulation, as an idealized comparison with the paper; large-table access cost is also parameterized separately."
        ],
        "subset_count": count(subsets, "32-element subsets"),
        "rotation_orbits": count(representatives, "orbits"),
        "coefficient_recovery": count(prefix_calls, "six-round Hsilu prefix evaluations"),
        "coefficient_recovery_seven_round_equivalent": count(prefix_equivalent, "seven-round Hsilu equivalents, linear round scaling"),
        "product_accumulation": count(pair_operations, "64-bit AND-XOR accumulation pairs"),
        "time_parameterization": "T = (6/7)*prefix_calls + eta*pair_operations + lower-order work; eta is seven-round Hsilu equivalents per 64-bit pair",
        "xor_only_eta": eta,
        "product_xor_only_equivalent": count(eta * pair_operations, "seven-round Hsilu equivalents"),
        "main_terms_xor_only_total": count(prefix_equivalent + eta * pair_operations, "seven-round Hsilu equivalents"),
        "sensitivity_two_xor_units_per_pair": count(prefix_equivalent + 2 * eta * pair_operations, "seven-round Hsilu equivalents"),
        "sensitivity_one_permutation_per_pair": count(prefix_equivalent + pair_operations, "seven-round Hsilu equivalents"),
        "data": count(1 << 64, "distinct-nonce queries, one 128-bit output block per query"),
        "raw_observed_data_bits": count(1 << 71, "bits, before streaming XOR accumulation"),
        "dense_coefficient_storage": count(memory_bits, "bits"),
        "dense_coefficient_storage_bytes": count(memory_bits // 8, "bytes"),
        "online_work": "2^64 oracle queries + L*2^64 candidate verification work, where L is the number of K1 candidates; approx 2^65 with L=1 under the usual unit convention",
        "matching_precomputation": "Convert the recovered 64-variable signature ANF to its truth table using a fast Mobius transform, then build a signature-to-key index. Do not evaluate an exponentially large ANF independently for every key.",
        "ascon_paper_2021_194_section_6_2": {
            "reported_time_log2": 123.28,
            "reported_data_log2": 64,
            "reported_memory_bits_log2": 101,
            "coefficient_interpolation_points": ascon_points,
            "coefficient_only_ascon_to_hsilu_operation_count_ratio": (subsets / representatives) * (ascon_points / (1 << 32)),
            "caution": "The two time units are their respective cipher's seven-round permutations, not an assertion of equal wall-clock cost."
        },
        "limitations": [
            "No full seven-round signature, collision distribution, or full key-recovery success rate has been computed.",
            "The degree-32 witness for a six-round intermediate coefficient does not establish that the final seven-round signature is nonconstant or discriminating.",
            "The total neglects coefficient-cube XOR accumulation and setup overhead, like the paper's leading operation estimate; decimals are not hardware predictions.",
            "These are operation counts for the specified interpolation and storage schedule, not optimality claims."
        ]
    }
    destination = Path(__file__).resolve().parents[1] / "results/hsilu_7r_ppm_rotation_complexity.json"
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    for field in ["coefficient_recovery", "product_accumulation", "main_terms_xor_only_total", "data", "dense_coefficient_storage"]:
        print(field, result[field])


if __name__ == "__main__":
    main()
