"""Complete cost schedules and reduced-instance validation for Hsilu PPM.

Two schedules:
 A. Materialize 16+16 -> 32 key ANFs, then 32+32 -> 64, with rotations.
 B. Evaluate the same leading-degree circuit on the three local states shared
    by all (nonce,key) assignments, then contract 3 -> 2 on each coordinate.

The full 64-bit instance is counted, NOT executed. Candidate discrimination
of the actual seven-round signature remains an explicit assumption.
"""

import itertools
import json
import math
import random
import time
from collections import Counter
from pathlib import Path

from permutations import HN1, HN2, HM, HR, RC, permute_window


def rol(x, s, n):
    s %= n
    mask = (1 << n) - 1
    return ((x << s) | (x >> ((n - s) % n))) & mask


def linear(a, n):
    a = list(a)
    for u, v, s in zip(HN1, HN2, HM):
        a.append(a[u] ^ rol(a[v], s, n))
    return [rol(a[13 + i], HR[i], n) for i in range(5)]


def round_full(a, rc, n):
    a = list(a)
    a[0] ^= rc & ((1 << n) - 1)
    return linear([a[i] ^ ((~a[(i + 1) % 5]) & a[(i + 2) % 5])
                   for i in range(5)], n)


def pack(a, n):
    return sum(v << (n * i) for i, v in enumerate(a))


def unpack(a, n, words=5):
    return [(a >> (n * i)) & ((1 << n) - 1) for i in range(words)]


def rotate_packed(a, s, n, words=5):
    return pack([rol(v, s, n) for v in unpack(a, n, words)], n)


def leading_value(x, b, rounds, n):
    # First-round nonce part: (0, (1+K1)*x, 0, x, 0).
    a = linear([0, b, 0, x, 0], n)
    for _ in range(rounds - 1):
        a = linear([a[(i + 1) % 5] & a[(i + 2) % 5]
                    for i in range(5)], n)
    return pack(a, n)


def phi(n):
    return sum(math.gcd(i, n) == 1 for i in range(1, n + 1))


def necklaces(n, alphabet=3):
    # Standard prenecklace recursion. Each output includes its actual period.
    a = [0] * (n + 1)
    def visit(t, period):
        if t > n:
            if n % period == 0:
                yield tuple(a[1:]), period
        else:
            a[t] = a[t - period]
            yield from visit(t + 1, period)
            for value in range(a[t - period] + 1, alphabet):
                a[t] = value
                yield from visit(t + 1, t)
    yield from visit(1, 1)


def ternary_table(n, rounds):
    powers = [3 ** j for j in range(n)]
    table = [None] * (3 ** n)
    calls = 0
    for digits, period in necklaces(n):
        x = sum((value != 0) << j for j, value in enumerate(digits))
        b = sum((value == 2) << j for j, value in enumerate(digits))
        value = leading_value(x, b, rounds, n)
        calls += 1
        index = sum(digits[j] * powers[j] for j in range(n))
        for s in range(period):
            assert table[index] is None
            table[index] = rotate_packed(value, s, n)
            # Rotate a base-3 address in O(1) fixed-width arithmetic. This
            # avoids re-encoding all n digits for every orbit member.
            index = 3 * index - digits[n - 1 - s] * (3 ** n - 1)
    assert all(value is not None for value in table)
    expected = sum(phi(d) * 3 ** (n // d) for d in range(1, n + 1)
                   if n % d == 0) // n
    assert calls == expected
    return table, calls


def contract(table, n):
    # k=0 permits local states 0 or 2; k=1 permits states 0 or 1.
    a = list(table)
    stride = 1
    fibers = 0
    for _ in range(n):
        out = [0] * (len(a) * 2 // 3)
        for group in range(len(a) // (3 * stride)):
            source = group * 3 * stride
            target = group * 2 * stride
            for offset in range(stride):
                v0, v1, v2 = [a[source + j * stride + offset] for j in range(3)]
                out[target + offset] = v0 ^ v2
                out[target + stride + offset] = v0 ^ v1
                fibers += 1
        a = out
        stride *= 2
    assert fibers == 3 ** n - 2 ** n
    return a


def contract_in_place(a, n):
    """3->2 tensor contraction without a second full-size table.

Read each triple before overwriting it. The destination is at or behind
the source, and writes never touch a later unread triple component.
"""
    stride = 1
    for _ in range(n):
        current_size = len(a)
        for group in range(current_size // (3 * stride)):
            source = group * 3 * stride
            target = group * 2 * stride
            for offset in range(stride):
                v0 = a[source + offset]
                v1 = a[source + stride + offset]
                v2 = a[source + 2 * stride + offset]
                a[target + offset] = v0 ^ v2
                a[target + stride + offset] = v0 ^ v1
        del a[current_size * 2 // 3:]
        stride *= 2
    return a


def positions(mask, n):
    return [i for i in range(n) if (mask >> i) & 1]


def deposit(value, pos):
    return sum(((value >> j) & 1) << bit for j, bit in enumerate(pos))


def compress(value, pos):
    return sum(((value >> bit) & 1) << j for j, bit in enumerate(pos))


def fmt(a):
    a = list(a)
    step = 1
    while step < len(a):
        for i in range(len(a)):
            if i & step:
                a[i] ^= a[i ^ step]
        step <<= 1
    return a


def subset_representatives(n, weight):
    reps = {}
    for pos in itertools.combinations(range(n), weight):
        mask = sum(1 << j for j in pos)
        orbit = {rol(mask, s, n) for s in range(n)}
        reps[min(orbit)] = len(orbit)
    return reps


def get_coefficients(tables, mask, n):
    rep = min(rol(mask, s, n) for s in range(n))
    shift = next(s for s in range(n) if rol(rep, s, n) == mask)
    target_pos, source_pos = positions(mask, n), positions(rep, n)
    return [rotate_packed(tables[rep][compress(rol(deposit(i, target_pos), -shift, n),
                                               source_pos)], shift, n)
            for i in range(1 << len(target_pos))]


def quadratic_pair(a, b, n):
    aa, bb = unpack(a, n), unpack(b, n)
    return pack([aa[(i + 1) % 5] & bb[(i + 2) % 5] for i in range(5)], n)


def materialized_ppm(n, rounds):
    quarter, half = n // 4, n // 2
    base = {}
    for rep in subset_representatives(n, quarter):
        pos = positions(rep, n)
        truth = []
        for key in range(1 << quarter):
            keymask = deposit(key, pos)
            value = 0
            for cube in range(1 << quarter):
                x = deposit(cube, pos)
                value ^= leading_value(x, x & ~keymask, rounds - 2, n)
            truth.append(value)
        base[rep] = fmt(truth)
    middle = {}
    middle_reps = subset_representatives(n, half)
    for rep in middle_reps:
        pos = positions(rep, n)
        accum = [0] * (1 << half)
        for subset in itertools.combinations(pos, quarter):
            j = sum(1 << bit for bit in subset)
            jp = positions(rep ^ j, n)
            a = get_coefficients(base, j, n)
            b = get_coefficients(base, rep ^ j, n)
            for ka, va in enumerate(a):
                ma = compress(deposit(ka, subset), pos)
                for kb, vb in enumerate(b):
                    target = ma | compress(deposit(kb, jp), pos)
                    accum[target] ^= quadratic_pair(va, vb, n)
        middle[rep] = [pack(linear(unpack(v, n), n), n) for v in accum]
    grouped = {}
    for rep, period in middle_reps.items():
        accum = grouped.setdefault(period, [0] * (1 << n))
        pos, comp_pos = positions(rep, n), positions(((1 << n) - 1) ^ rep, n)
        a = middle[rep]
        b = get_coefficients(middle, ((1 << n) - 1) ^ rep, n)
        for ka, va in enumerate(a):
            ma = deposit(ka, pos)
            for kb, vb in enumerate(b):
                accum[ma | deposit(kb, comp_pos)] ^= quadratic_pair(va, vb, n)
    anf = [0] * (1 << n)
    for period, accum in grouped.items():
        for keymask, value in enumerate(accum):
            for shift in range(period):
                anf[rol(keymask, shift, n)] ^= rotate_packed(value, shift, n)
    anf = [pack(linear(unpack(v, n), n), n) for v in anf]
    return fmt(anf), {"base_representatives": len(base),
                      "middle_representatives": len(middle),
                      "middle_periods": sorted(grouped)}


def toy_check(n, rounds):
    table, calls = ternary_table(n, rounds)
    truth = contract(table, n)
    assert contract_in_place(list(table), n) == truth
    ppm_truth, counts = materialized_ppm(n, rounds)
    assert truth == ppm_truth
    for k1 in range(1 << n):
        direct = 0
        for x in range(1 << n):
            a = [0x81 & ((1 << n) - 1), 0xA3 & ((1 << n) - 1), k1, x, 0]
            for rc in RC[-rounds:]:
                a = round_full(a, rc, n)
            direct ^= pack(a, n)
        assert truth[k1] == direct
    rate = [v & ((1 << (2 * n)) - 1) for v in truth]
    buckets = Counter(rate)
    # Matching retains all colliding candidates, rather than overwriting them.
    matching = {}
    for key, value in enumerate(rate):
        matching.setdefault(value, []).append(key)
    assert all(key in matching[value] for key, value in enumerate(rate))
    return {"word_bits": n, "rounds": rounds, "keys_exhausted": 1 << n,
            "ternary_points": 3 ** n, "ternary_orbit_evaluations": calls,
            "all_state_bits_match_direct_and_both_schedules": True,
            "in_place_contraction_matches": True,
            "rate_distinct_signatures": len(buckets),
            "rate_max_bucket_size": max(buckets.values()), "ppm_counts": counts}


def actual_width_check():
    # True 64-bit Hsilu words; an 8-dimensional maximal cube at four rounds.
    pos = [1, 4, 13, 20, 31, 40, 52, 63]
    table = []
    for digits in itertools.product(range(3), repeat=8):
        # Reverse product order to make position 0 the least significant trit.
        digits = digits[::-1]
        x = sum((d != 0) << p for d, p in zip(digits, pos))
        b = sum((d == 2) << p for d, p in zip(digits, pos))
        table.append(leading_value(x, b, 4, 64))
    truth = contract(table, 8)
    assert contract_in_place(list(table), 8) == truth
    rng = random.Random(20261020)
    k0, outside = rng.getrandbits(64), rng.getrandbits(64)
    outside &= ~sum(1 << p for p in pos)
    nonces = [deposit(value, pos) for value in range(256)]
    for local_key in range(256):
        k1 = outside | deposit(local_key, pos)
        direct = 0
        for nonce in nonces:
            direct ^= pack(permute_window([0x00001000808C0001, k0, k1, nonce, 0],
                                          "hsilu", 8, 4), 64)
        assert direct == truth[local_key]
    return {"word_bits": 64, "rounds": 4, "cube_positions": pos,
            "local_keys_exhausted": 256, "nonce_points_per_key": 256,
            "all_320_bits_match_reference": True,
            "in_place_contraction_matches": True}


def number(value, unit):
    return {"value": value, "log2": math.log2(value), "unit": unit}


def full_costs():
    n = 64
    r32 = sum(phi(d) * math.comb(n // d, 32 // d)
              for d in [1, 2, 4, 8, 16, 32]) // n
    r16 = sum(phi(d) * math.comb(n // d, 16 // d)
              for d in [1, 2, 4, 8, 16]) // n
    middle = 5 * r32 * math.comb(32, 16) * (1 << 32)
    outer = 5 * r32 * (1 << 64)
    points = 3 ** n
    orbits = sum(phi(d) * 3 ** (n // d) for d in [1, 2, 4, 8, 16, 32, 64]) // n
    fibers = points - (1 << n)
    # Two 64-bit words per grid value. Each 3->2 fiber reads 3 vectors
    # and writes 2 vectors, or 10 word accesses and 4 word XORs.
    accesses = 2 * points + 10 * fibers
    return {
        "common": {
            "data": number(1 << 64, "distinct-nonce queries, each observing 128 bits"),
            "signature_truth_table_storage": number(128 * (1 << 64), "bits"),
            "deterministic_matching_scan": number(1 << 64, "128-bit signature comparisons; two word reads per comparison"),
            "candidate_index_storage": number(192 * (1 << 64), "bits, payload only: 128-bit signature + 64-bit key"),
            "candidate_verification": "L * 2^64 local candidate checks; L is the actual seven-round signature bucket size",
            "full_recovery_condition": "L * 2^64 must not dominate the reported preprocessing time; no seven-round bucket distribution is certified"
        },
        "A_materialized_four_factor_ppm": {
            "base_prefix_calls": number(r16 * (1 << 32), "five-round prefix calls"),
            "middle_16_plus_16": number(middle, "64-bit AND-XOR accumulations"),
            "outer_32_plus_32": number(outer, "64-bit AND-XOR accumulations"),
            "combined_accumulations": number(middle + outer, "64-bit AND-XOR accumulations"),
            "main_memory_accesses_simple_loop": number(4 * (middle + outer), "64-bit accesses: read two factors, read/write accumulator"),
            "time_one_access_one_round": number(4 * (middle + outer) / 7,
                                                "seven-round equivalents, memory-access-dominated paper-style model"),
            "time_one_accumulation_one_round": number((middle + outer) / 7,
                                                     "seven-round equivalents, alternative aggregate primitive model"),
            "time_xor_only": number((middle + outer) / 126,
                                    "seven-round equivalents, AND/memory/indexing overhead ignored"),
            "dominant_storage": number(320 * r32 * (1 << 32), "bits"),
            "additional_lower_order_work": "6 orbit-length accumulators; <=126 rotation symmetrizations; final linear layer; 64-variable fast Mobius transform; matching-table index"
        },
        "B_shared_ternary_batch_evaluation": {
            "points": number(points, "ternary assignments"),
            "rotation_representatives": number(orbits, "ternary necklaces, including short periods"),
            "representative_evaluation_upper_bound": number(orbits, "seven-round Hsilu calls; actual circuit is first linear layer plus six quadratic-only rounds"),
            "grid_fill_writes": number(2 * points, "64-bit writes"),
            "contraction_fibers": number(fibers, "3-vector to 2-vector contractions"),
            "contraction_xors": number(4 * fibers, "64-bit XORs"),
            "contraction_memory_accesses": number(10 * fibers, "64-bit reads/writes"),
            "total_grid_memory_accesses": number(accesses, "64-bit reads/writes"),
            "leading_time_one_access_one_round": number(orbits + accesses / 7,
                                                       "seven-round equivalents; local indexing/rotations and XOR overhead are separate O(3^64) word operations"),
            "leading_time_plus_contraction_xor_cost": number(orbits + accesses / 7 + 4 * fibers / 126,
                                                            "seven-round equivalents; excludes address-generation and rotation arithmetic"),
            "storage_peak_in_place": number(128 * points, "bits; dense two-word ternary table, contraction may overwrite its input"),
            "matching": "After 64 contractions, entries are the actual signature truth table for all 2^64 local key words. The default schedule scans it against the observed signature in 2^64 comparisons and retains every matching key; a hash index is optional, not required.",
            "why_no_missing_key_evaluation_factor": "The 3->2 contraction computes every key's cube sum simultaneously, with no separate evaluation of a large expression at each key.",
            "remaining_word_arithmetic": "O(3^64) for necklace generation, orbit scatter, rotations and fixed-width address calculations; decimals of the leading cost model are not wall-clock predictions"
        },
        "paper_2025_1029_reported": {"time_log2": 95.96, "data_log2": 64,
                                      "memory_bits_log2": 104.85,
                                      "time_unit": "seven-round Ascon permutations",
                                      "status": "Author-reported estimate; the all-key evaluation step is not independently certified here"}
    }


def main():
    start = time.monotonic()
    tests = [toy_check(4, 3), toy_check(8, 4), actual_width_check()]
    report = {
        "status": "Complete counted schedules under explicit cost models, with reduced-instance verification; not a practical execution of the seven-round instance",
        "model": "Hsilu (IV,K0,K1,N0,0), 7 initialization rounds, full 64-dimensional N0 cube, first 128 output bits observable",
        "quartic_screening": {
            "used_to_reduce_these_costs": False,
            "follow_up_report": str(Path(__file__).resolve().parents[1] / "results" / "hsilu_7r_quartic_screen.json"),
            "meaning": "These are explicit unpruned schedules. A separate first-five-round MILP screen audits whether quartic terms can be excluded."
        },
        "tests": tests,
        "costs": full_costs(),
        "proof": [
            "First-round state is B(K)+L(0,(1+K1)*x,0,x,0). A degree-64 nonce monomial after the remaining six quadratic rounds requires 64 nonce-bearing leaves; any first-state constant or later linear/constant branch cannot contribute.",
            "Hence the full cube equals the cube of G(x,b), where b=(1+K1)*x and G is the first linear layer followed by six quadratic-only rounds. G is rotation-covariant.",
            "Each coordinate of (x,b) assumes only (0,0),(1,0),(1,1). At fixed K1 bit 0, cube summation chooses (0,0) and (1,1); at bit 1 it chooses (0,0) and (1,0). Tensoring this local 3->2 XOR map over 64 coordinates gives every key's cube sum exactly.",
            "The contraction schedule is the pointwise-evaluation counterpart of retaining the degree-16 factors and their last-two-round circuit: lower-degree terms disappear under the full 64-dimensional sum, and common subexpressions are evaluated only once per grid point."
        ],
        "author_code_audit": {
            "repository": "https://github.com/keyrecoveryofascon/keyrecoveryofascon",
            "commit": "616f5d247bd9decdb39a74f2339f73a115392910",
            "observed": "The supplied experiment constructs symbolic six-round superpoly expressions; no seven-round all-key matching-table construction was found in the inspected source. This does not prove that no faster algorithm exists."
        },
        "limitations": [
            "No seven-round nonconstant-signature, collision-distribution, or complete recovery-success claim.",
            "The ternary batch algorithm is an additional derivation for this Hsilu transfer, not a result claimed in the Ascon paper and not a claim of novelty.",
            "Both memory estimates are for specified schedules, not minimum memory bounds.",
            "No claim of optimal time complexity or of a flaw in the published Ascon result."
        ],
        "elapsed_seconds": time.monotonic() - start
    }
    path = Path(__file__).resolve().parents[1] / "results/hsilu_7r_batch_audit.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"path": str(path), "tests": tests, "costs": report["costs"],
                      "elapsed_seconds": report["elapsed_seconds"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
