#!/usr/bin/env python3
"""Exhaust the matching-filtered E cases and independently check DRAT proofs."""
import argparse
import concurrent.futures
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import threading
import time

ROOT = Path(__file__).resolve().parent
KISSAT = os.environ.get("KISSAT", str(ROOT / "tools/kissat"))
CHECKER = ROOT / "tools/drat-trim/drat-trim"
lock = threading.Lock()
running = set()
stopped = threading.Event()


def interrupt(signum, frame):
    stopped.set()
    with lock:
        for process in running:
            if process.poll() is None:
                process.send_signal(signal.SIGINT)


def run_process(command, log, stdin=None):
    with log.open("wb") as output:
        process = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT,
                                   stdin=stdin)
        with lock:
            running.add(process)
        try:
            code = process.wait()
        finally:
            with lock:
                running.discard(process)
    return code


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--total", type=int, default=15)
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--min-matching", type=int, default=0)
    parser.add_argument("--max-matching", type=int, default=15)
    args = parser.parse_args()
    signal.signal(signal.SIGINT, interrupt)
    work = ROOT / f"effects_bound{args.total}"
    work.mkdir(exist_ok=True)
    base = (ROOT / f"total{args.total}.cnf").read_text()
    mapping = json.loads((ROOT / f"total{args.total}.map.json").read_text())
    cases = []
    for line in (ROOT / "effects_total15.tsv").read_text().splitlines()[1:]:
        weight, hexword, bound = line.split()
        if args.min_matching <= int(bound) <= min(args.total, args.max_matching):
            cases.append((int(weight), int(hexword, 16), int(bound)))
    cases.sort(key=lambda c: (-c[2], -c[0], c[1]))
    if args.limit:
        cases = cases[:args.limit]
    def solve_case(case):
        if stopped.is_set():
            return None
        weight, e, bound = case
        name = f"e{weight}_{e:016x}"
        prefix = work / name
        status = prefix.with_suffix(".json")
        if status.exists():
            old = json.loads(status.read_text())
            saved_proof=prefix.with_suffix(".drat.gz")
            if old.get("drat_verified") and saved_proof.exists() and \
                    hashlib.sha256(saved_proof.read_bytes()).hexdigest()==old.get("proof_sha256"):
                return old
        formula = base.replace(f"p cnf {mapping['variables']} {mapping['clauses']}\n",
                               f"p cnf {mapping['variables']} {mapping['clauses']+64}\n", 1)
        formula += "".join(f"{v if e >> j & 1 else -v} 0\n" for j,v in enumerate(mapping["E"]))
        cnf = prefix.with_suffix(".cnf")
        proof = prefix.with_suffix(".drat.gz")
        cnf.write_text(formula)
        start = time.monotonic()
        code = run_process([KISSAT, "-f", "--phase=0", str(cnf), str(proof)],
                           prefix.with_suffix(".solver.log"))
        result = {"effect": f"0x{e:016x}", "effect_weight": weight,
                  "matching_lower_bound": bound, "total_max": args.total,
                  "solver_exit_code": code, "solver_seconds": time.monotonic()-start,
                  "cnf_sha256": hashlib.sha256(formula.encode()).hexdigest(),
                  "cnf_generator": "base total bound CNF plus 64 E units",
                  "drat_verified": False}
        if code == 20 and not stopped.is_set():
            # drat-trim accepts a proof on stdin, including binary DRAT.
            decompressor = subprocess.Popen(["gzip", "-dc", str(proof)], stdout=subprocess.PIPE)
            checked = run_process([str(CHECKER), str(cnf), "-i"], prefix.with_suffix(".check.log"),
                                  stdin=decompressor.stdout)
            decompressor.stdout.close()
            unzip_code = decompressor.wait()
            text = prefix.with_suffix(".check.log").read_text()
            result["drat_verified"] = checked == 0 and unzip_code == 0 and "s VERIFIED" in text
            result["checker_exit_code"] = checked
            if not result["drat_verified"]:
                stopped.set()
            else:
                result["proof_sha256"] = hashlib.sha256(proof.read_bytes()).hexdigest()
                # CNF is exactly reproducible from the saved base plus units.
                with gzip.open(prefix.with_suffix(".cnf.gz"), "wb") as f:
                    f.write(formula.encode())
                cnf.unlink()
        elif code == 10:
            stopped.set()
            result["satisfiable"] = True
        status.write_text(json.dumps(result, indent=2)+"\n")
        return result
    def solve(case):
        # Multiple resumptions must never write the same proof concurrently.
        weight, e, bound = case
        claim = work / f"e{weight}_{e:016x}.lock"
        with claim.open("a") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            return solve_case(case)
    start = time.monotonic()
    print(f"Cases={len(cases)}, jobs={args.jobs}, total<={args.total}", flush=True)
    completed = 0
    verified = 0
    found_sat = False
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(solve, case) for case in cases]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result is None:
                continue
            completed += 1
            verified += int(result["drat_verified"])
            if completed % 20 == 0 or not result["drat_verified"] or completed == len(cases):
                print(f"Completed={completed}/{len(cases)}, last={result['effect']}, "
                      f"verified={result['drat_verified']}, elapsed={time.monotonic()-start:.1f}s",
                      flush=True)
            if result.get("satisfiable"):
                found_sat = True
                print("SAT witness in", work, result["effect"], flush=True)
                interrupt(None, None)
    print("Finished. Only records with drat_verified=true prove an E case.", flush=True)
    print(f"Verified={verified}/{len(cases)}", flush=True)
    return 10 if found_sat else (0 if verified == len(cases) else 130 if stopped.is_set() else 1)


if __name__ == "__main__":
    raise SystemExit(main())
