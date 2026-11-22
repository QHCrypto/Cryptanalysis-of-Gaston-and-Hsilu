#!/usr/bin/env python3
"""Verify the 6->6 rank certificate with forward L, and conditional sampling."""
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from six_six_prefix_search import (MASK, active_mask, ddt_table, embed, flat,
                                   linear, rank, state_chi, sxor)


def decode(obj):
    return tuple(int(w, 16) for w in obj["words"])


def main():
    result = ROOT / "results/six_six_fixed_tail"
    data = json.loads((result / "search.json").read_text())
    ddt = ddt_table()
    assert max(max(row) for row in ddt[1:]) == 8
    group = next(g for g in data["groups"].values() if g["noncancelling_rotations"] == [24])
    for s in group["rotations"]:
        d0, a2 = decode(s["d0"]), decode(s["a2"])
        inactive = MASK ^ active_mask(d0)
        columns = []
        # Independent of L^T and of the derivative-coefficient formula in search.
        for j in range(64):
            if not (inactive >> j) & 1:
                continue
            for i in range(5):
                y = linear(embed(j, 1 << i))
                linear_part = sxor(sxor(state_chi(y), state_chi(sxor(y, a2))), state_chi(a2))
                columns.append(flat(linear_part))
        assert rank(columns) == s["rank_certificate"]["second_derivative_rank"] == 12
    assert len(group["rotations"][7]["connections"]) == 4
    assert group["rotations"][7]["restricted_signed_prefix_sum"]["numerator"] == 0
    assert group["rotations"][24]["restricted_signed_prefix_sum"] == {
        "numerator": -1, "denominator": 1 << 24, "log2_abs": -24.0}
    with tempfile.TemporaryDirectory(prefix="gaston-six-six-") as tmp:
        outputs = []
        for threads in (1, 4):
            output = Path(tmp) / f"threads{threads}.csv"
            subprocess.run([str(ROOT / "build/conditional_trail_verify"),
                            str(result / "conditional_jobs.txt"), str(output), "65536", "2026101501", str(threads)],
                           check=True, stdout=subprocess.DEVNULL)
            with output.open() as f:
                rows = list(csv.DictReader(f))
            for row in rows:
                row.pop("threads")
                assert row["bad_first"] == row["bad_middle"] == "0"
            outputs.append(rows)
        assert outputs[0] == outputs[1]
    print("PASS 64 independent forward-L rank checks, signed cancellation, and thread-independent conditional sampling")


if __name__ == "__main__":
    main()
