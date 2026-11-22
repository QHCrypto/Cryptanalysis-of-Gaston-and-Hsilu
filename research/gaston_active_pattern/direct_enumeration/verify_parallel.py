#!/usr/bin/env python3
"""Independent orbit/certificate checks plus CPU/OpenCL cross-validation."""
import argparse
import csv
import json
import math
from pathlib import Path
import signal
import subprocess
import tempfile
import verify as reference

ROOT = Path(__file__).resolve().parent
PROGRAM = ROOT / "gaston_parallel"
MASK = (1 << 64) - 1


def orbit_count(k):
    return sum(math.comb(math.gcd(64, shift), k // (64 // math.gcd(64, shift)))
               for shift in range(64) if k % (64 // math.gcd(64, shift)) == 0) // 64


def validate(prefix, maps, expected=None, check_all_witnesses=True):
    data = json.loads(Path(str(prefix) + ".json").read_text())
    mode, side, k = data["mode"], data["fixed_side"], data["fixed_count"]
    mapping = maps[mode, side]
    reference.verify_certificate(data, mapping, side, k, expected, complete=data["enumeration_complete"])
    total = orbit_count(k) if data["symmetry"] else math.comb(64, k)
    assert data["position_sets"] == total
    assert data["candidate_count_estimate"] == total * 31 ** k
    seen = set()
    best = 65
    # Matrix columns let us independently map each witness by XORing its
    # nonzero input bits, without 320 row dot products per certificate.
    images = [sum(((mapping[r] >> c) & 1) << r for r in range(320)) for c in range(320)]
    last_index = -1
    with Path(str(prefix) + ".positions.tsv").open() as file:
        for record in csv.DictReader(file, delimiter="\t"):
            index = int(record["position_index"])
            mask = int(record["position_mask_hex"], 16)
            minimum = int(record["minimum"])
            code = int(record["base31_value_code"])
            assert last_index < index < total
            last_index = index
            assert mask.bit_count() == k and mask not in seen
            seen.add(mask)
            assert 0 <= code < 31 ** k and 0 <= minimum <= 64
            if data["symmetry"]:
                assert all(mask <= ((mask << s) | (mask >> (64-s))) & MASK for s in range(1, 64))
            if check_all_witnesses:
                image = 0
                for col in range(64):
                    if mask >> col & 1:
                        value = code % 31 + 1
                        code //= 31
                        for row in range(5):
                            if value >> row & 1:
                                image ^= images[64*row+col]
                assert code == 0 and reference.weight(image) == minimum
            best = min(best, minimum)
    assert len(seen) == data["completed_position_sets"]
    assert data["visited"] == len(seen) * 31 ** k
    assert len(seen) == data["cpu_position_sets"] + data["gpu_position_sets"]
    assert best == data["objective"]
    if data["enumeration_complete"]:
        # Every saved set is a distinct canonical orbit representative and
        # the number equals Burnside's count: no position orbit is omitted.
        assert len(seen) == total
    print(f"PASS {prefix}: optimum/candidate={best}, {len(seen)}/{total} position sets, "
          f"complete={data['enumeration_complete']}; independent witnesses and orbit coverage", flush=True)
    return data


def run(prefix, mode, side, k, backend, extra=()):
    args = [str(PROGRAM), "--mode", mode, f"--{side}-active", str(k),
            "--backend", backend, "--threads", "4", "--prefix", str(prefix), *extra]
    proc = subprocess.run(args, capture_output=True, text=True)
    assert proc.returncode == (130 if "--limit-positions" in extra else 0), proc.stdout + proc.stderr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", nargs="*", help="Validate saved result prefixes instead of running small tests")
    args = parser.parse_args()
    forward, transposed = reference.forward_matrix()
    maps = {("difference", "input"): forward, ("difference", "output"): reference.inverse(forward),
            ("linear", "input"): reference.inverse(transposed), ("linear", "output"): transposed}
    if args.results is not None:
        for prefix in args.results:
            validate(prefix, maps)
        return
    with tempfile.TemporaryDirectory(prefix="gaston_parallel_check_") as tmp:
        tmp = Path(tmp)
        for (mode, side), mapping in maps.items():
            expected = reference.exact_small_optima(mapping)
            for k in range(4):
                cpu, gpu = tmp / "cpu", tmp / "gpu"
                run(cpu, mode, side, k, "cpu")
                run(gpu, mode, side, k, "opencl")
                a = validate(cpu, maps, expected[k] if k < 3 else None)
                b = validate(gpu, maps, a["objective"])
                assert Path(str(cpu)+".positions.tsv").read_bytes() == Path(str(gpu)+".positions.tsv").read_bytes()
                if (mode, side, k) == ("difference", "output", 3):
                    assert b["objective"] == 39
            full = tmp / "full"
            run(full, mode, side, 2, "hybrid", ("--no-symmetry",))
            validate(full, maps, expected[2])
        # A real incomplete search must never be labelled optimal.
        partial = tmp / "partial"
        run(partial, "difference", "output", 5, "hybrid", ("--limit-positions", "128"))
        validate(partial, maps)
        interrupted = tmp / "interrupted"
        proc = subprocess.Popen([str(PROGRAM), "--output-active", "5", "--backend", "hybrid",
                                 "--threads", "4", "--prefix", str(interrupted)],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            for line in proc.stdout:
                if line.startswith("Improved:"):
                    proc.send_signal(signal.SIGINT)
                    break
            stdout, stderr = proc.communicate(timeout=20)
            assert proc.returncode == 130, stdout + stderr
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
        validate(interrupted, maps)
    print("PASS: CPU/OpenCL per-position equality in all four directions for K=0..3; "
          "full K=2; periodic position orbits; partial and SIGINT reporting.")


if __name__ == "__main__":
    main()
