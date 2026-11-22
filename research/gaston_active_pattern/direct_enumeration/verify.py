#!/usr/bin/env python3
"""Independent exact reference; Python standard library only, no Gurobi."""
import json
import math
from pathlib import Path
import signal
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PROGRAM = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "gaston_enumerate"
N = 320
MASK = (1 << 64) - 1


def forward_matrix():
    # Track the destinations of each input bit. Independent of the C++
    # word-level implementation and its std::bitset elimination.
    east, west, twin = (0, 60, 22, 27, 4), (0, 56, 31, 46, 43), (25, 32, 52, 60, 63)
    columns = []
    for row in range(5):
        for bit in range(64):
            position = (bit + east[row]) % 64
            image = 1 << (64 * row + (position + west[row]) % 64)
            for shift in (0, 1, twin[row], twin[row] + 18):
                for out_row in range(5):
                    image ^= 1 << (64 * out_row + (position + shift + 23 + west[out_row]) % 64)
            columns.append(image)
    rows = [sum(((columns[c] >> r) & 1) << c for c in range(N)) for r in range(N)]
    return rows, columns  # The columns also give the rows of the transpose.


def inverse(rows):
    augmented = [row | (1 << (N + i)) for i, row in enumerate(rows)]
    for col in range(N):
        pivot = next(i for i in range(col, N) if augmented[i] >> col & 1)
        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        for i in range(N):
            if i != col and augmented[i] >> col & 1:
                augmented[i] ^= augmented[col]
    assert all((row & ((1 << N) - 1)) == 1 << i for i, row in enumerate(augmented))
    return [row >> N for row in augmented]


def apply(matrix, value):
    return sum(((row & value).bit_count() % 2) << i for i, row in enumerate(matrix))


def active_mask(value):
    return (value | value >> 64 | value >> 128 | value >> 192 | value >> 256) & MASK


def weight(value):
    return active_mask(value).bit_count()


def exact_small_optima(matrix):
    # Enumerate independently with Python integer vectors and XOR, anchoring
    # one nonzero column at zero. Counts 0, 1 and 2 are all covered.
    table = []
    for col in range(64):
        images = []
        for pattern in range(1, 32):
            fixed = sum(((pattern >> row) & 1) << (64 * row + col) for row in range(5))
            images.append(apply(matrix, fixed))
        table.append(images)
    one = min(weight(value) for value in table[0])
    two = min(weight(x ^ y) for col in range(1, 64) for x in table[0] for y in table[col])
    return (0, one, two)


def verify_certificate(data, mapping, side, k, expected=None, complete=True):
    assert data["fixed_side"] == side and data["fixed_count"] == k
    assert data["fixed_metric"] == "active_sbox_columns" and data["verified"]
    states = {}
    for end in ("input", "output"):
        state = sum(int(word, 16) << (64 * i) for i, word in enumerate(data[end + "_words"]))
        states[end] = state
        assert data[end + "_active"] == weight(state)
        assert data[end + "_bit_weight"] == state.bit_count()
        assert data[end + "_columns"] == [j for j in range(64) if active_mask(state) >> j & 1]
    other = "output" if side == "input" else "input"
    assert weight(states[side]) == k
    assert apply(mapping, states[side]) == states[other]
    assert data["objective"] == weight(states[other])
    assert data["enumeration_complete"] == complete and data["proven_optimal"] == complete
    if expected is not None:
        assert data["objective"] == expected
    if data.get("symmetry_method") == "position_rotation_orbits":
        fixed_sets = sum(math.comb(math.gcd(64, s), k // (64 // math.gcd(64, s)))
                         for s in range(64) if k % (64 // math.gcd(64, s)) == 0) // 64
        expected_visits = fixed_sets * 31 ** k
    else:
        expected_visits = (math.comb(63, k - 1) if data["symmetry"] and k else math.comb(64, k)) * 31 ** k
    assert data["candidate_count_estimate"] == expected_visits
    if complete:
        assert data["visited"] == expected_visits
    else:
        assert 0 < data["visited"] < expected_visits


def run(args):
    proc = subprocess.run([str(PROGRAM), *args], text=True, capture_output=True)
    if proc.returncode:
        raise AssertionError(proc.stdout + proc.stderr)


def main():
    forward, transposed = forward_matrix()
    maps = {
        ("difference", "input"): forward,
        ("difference", "output"): inverse(forward),
        ("linear", "input"): inverse(transposed),
        ("linear", "output"): transposed,
    }
    one_known = {("difference", "input"): 18, ("difference", "output"): 57,
                 ("linear", "input"): 59, ("linear", "output"): 17}
    with tempfile.TemporaryDirectory(prefix="gaston_enum_check_") as tmp:
        for (mode, side), mapping in maps.items():
            expected = exact_small_optima(mapping)
            assert expected[1] == one_known[mode, side]
            if (mode, side) == ("difference", "output"):
                assert expected[2] == 48
            for k in range(3):
                prefix = Path(tmp) / "auto_created" / f"{mode}_{side}_{k}"
                run(["--mode", mode, f"--{side}-active", str(k), "--prefix", str(prefix)])
                data = json.loads(prefix.with_suffix(".json").read_text())
                verify_certificate(data, mapping, side, k, expected[k])
                assert data["mode"] == mode
            # Full 1,937,376-state enumeration checks the symmetry reduction.
            prefix = Path(tmp) / f"{mode}_{side}_full"
            run(["--mode", mode, f"--{side}-active", "2", "--no-symmetry", "--prefix", str(prefix)])
            data = json.loads(prefix.with_suffix(".json").read_text())
            verify_certificate(data, mapping, side, 2, expected[2])
            print(f"PASS {mode:10s} fixed {side:6s}: K=0,1,2 -> {expected}; full K=2 agrees", flush=True)
        # Wait for an actual incumbent, then test a graceful user interruption.
        prefix = Path(tmp) / "interrupted"
        proc = subprocess.Popen([str(PROGRAM), "--output-active", "4", "--prefix", str(prefix)],
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            for line in proc.stdout:
                if line.startswith("Improved:"):
                    proc.send_signal(signal.SIGINT)
                    break
            stdout, stderr = proc.communicate(timeout=10)
            assert proc.returncode == 130, stdout + stderr
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
        data = json.loads(prefix.with_suffix(".json").read_text())
        verify_certificate(data, maps["difference", "output"], "output", 4, complete=False)
    for args in [[], ["--output-active", "65"], ["--output-active", "-1"],
                 ["--input-active", "1", "--output-active", "2"],
                 ["--mode", "invalid", "--output-active", "2"],
                 ["--output-active", "2bits"], ["--output-active"], ["--unknown"]]:
        result = subprocess.run([str(PROGRAM), *args], capture_output=True, text=True)
        assert result.returncode != 0 and "错误:" in result.stderr
    print("PASS: independent exact Python reference; 16 complete searches; symmetry/full agreement; "
          "all witness states; SIGINT incumbent reporting; automatic result directories; 8 argument checks.")


if __name__ == "__main__":
    main()
