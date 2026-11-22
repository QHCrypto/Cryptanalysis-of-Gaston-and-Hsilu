#!/usr/bin/env python3
"""Verify and run reduced kernel MILPs in order C1=6,7,8; C1=9 is analytically excluded."""
import argparse
from datetime import datetime
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]
ORDER = (6, 7, 8)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def coverage():
    profile = json.loads((ROOT / "references/activity_profiles.json").read_text())
    f = profile["forward_lower_by_exact_input_count"]
    r = profile["reverse_lower_by_exact_output_count"]
    edges = {}
    for n in (0, 1):
        edges[n] = [(a, b) for a in range(1, 65) for b in range(1, 65)
                    if b >= f[a] and a >= r[b] and a+b >= 12+4*n
                    and (not n or b != 6 or a >= 13)
                    and (n or (a <= 5*b and b <= 5*a))]
    branches = {}
    remaining = []
    for n0, n1 in itertools.product((0, 1), repeat=2):
        incoming = {b: min(a for a, y in edges[n0] if y == b)
                    for b in {y for _, y in edges[n0]}}
        lower = {b: 2*(incoming[b]+b+min(c for x, c in edges[n1] if x == b))
                 for b in incoming if any(x == b for x, _ in edges[n1])}
        branches[f"{n0}{n1}"] = lower
        remaining.extend((n0, n1, k) for k, bound in lower.items() if bound < 44)
    assert sorted(remaining) == [(0, 0, k) for k in (6, 7, 8, 9)]
    audit_path = ROOT / "results/kernel_structure_audit_20261024.json"
    audit = json.loads(audit_path.read_text())
    assert audit["status"] == "PASS" and audit["profile_sha256"] == digest(ROOT / "references/activity_profiles.json")
    assert audit["source_sha256"] == digest(ROOT / "src/audit_kernel_structure.py")
    assert audit["branch_results"]["9"]["count_and_even_bit_weight_lower_bound"] == 45
    assert not audit["branch_results"]["9"]["remaining_types_under_44"]
    return {"status": "PASS", "threshold": 44, "integer_weight_cap": 43,
            "argument": "Every true DDT trail embeds in lossy at no greater integer weight; each active chi column costs at least two.",
            "remaining_subproblems": [(0,0,k) for k in ORDER],
            "analytically_excluded_middle9_weight_lower_bound": 45,
            "kernel_audit_sha256": digest(audit_path),
            "count_only_weight_lower_bounds": branches,
            "source_profile_sha256": digest(ROOT / "references/activity_profiles.json"),
            "source_proof_scope": "Uses saved completed single-layer results; historical proofs not replayed"}


def prepare():
    run = ROOT / "results/kernel_campaigns" / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8])
    run.mkdir(parents=True, exist_ok=False)
    save(run / "coverage.json", coverage())
    verification_paths = {}
    for k in ORDER:
        directory = run / f"validation_k{k}"
        subprocess.run([str(ROOT / "gaston_3round_lossy"), "--output", str(directory),
                        "--threads", "16", "--kernel-both", "--middle-active", str(k),
                        "--weight-cap", "43", "--kernel-reduced", "--export-only"], check=True)
        subprocess.run([sys.executable, str(ROOT / "src/verify_lossy_model.py"), str(directory)], check=True)
        verification_paths[str(k)] = str(directory / "verification.json")
    sources = ["gaston_3round_lossy", "gaston_3round_lossy.cpp", "Makefile",
               "src/activity_profile_data.hpp", "src/verify_lossy_model.py",
               "src/kernel_campaign.py", "src/build_activity_profiles.py",
               "references/activity_profiles.json", "src/kernel_type_data.hpp",
               "src/build_kernel_types.py", "src/audit_kernel_structure.py",
               "results/kernel_structure_audit_20261024.json"]
    hashes = {}
    for name in sources:
        source, destination = ROOT / name, run / "source" / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        hashes[name] = digest(destination)
    manifest = {"state": "verified_not_running", "run_directory": str(run), "order": ORDER,
                "kernel_reduced": True, "analytically_excluded": {"9": {"weight_lower_bound": 45}},
                "threads_per_solve": 16, "maximum_concurrent_solves": 1,
                "time_limit": None, "external_timeout": None,
                "verification": verification_paths, "source_sha256": hashes,
                "coverage_sha256": digest(run / "coverage.json")}
    save(run / "manifest.json", manifest)
    save(ROOT / "results/latest_kernel_campaign.json", manifest)
    print("PREPARED_CAMPAIGN=" + str(run), flush=True)


def validate_manifest(run):
    manifest = json.loads((run / "manifest.json").read_text())
    assert digest(run / "coverage.json") == manifest["coverage_sha256"]
    for name, expected in manifest["source_sha256"].items():
        assert digest(run / "source" / name) == expected, name
    for k, path in manifest["verification"].items():
        report = json.loads(Path(path).read_text())
        assert report["status"] == "PASS"
        assert report["subproblem"] == {"kernel_both": True, "middle_active": int(k), "weight_cap": 43}
        assert report["kernel_reduced"] and report["kernel_type_verification"]["type_count"] == {6:62,7:4,8:1}[int(k)]
        for name in ("gaston_3round_lossy", "gaston_3round_lossy.cpp", "src/activity_profile_data.hpp", "references/activity_profiles.json",
                     "src/kernel_type_data.hpp", "results/kernel_structure_audit_20261024.json"):
            original = str(Path(report_path_root(report)) / name)
            assert report["sha256"][original] == manifest["source_sha256"][name]
        assert digest(Path(path).parent / "model.mps") == report["sha256"][str(Path(path).parent / "model.mps")]
    return manifest


def report_path_root(report):
    return next(str(Path(path).parent) for path in report["sha256"] if Path(path).name == "gaston_3round_lossy.cpp")


def worker(run):
    manifest = validate_manifest(run)
    status = {"state": "running", "worker_pid": os.getpid(), "current_k": None,
              "child_pid": None, "completed": {}, "pending": list(manifest["order"]),
              "kernel_reduced": True, "analytically_excluded": manifest["analytically_excluded"]}
    child = None
    stopped = False

    def handle_stop(signum, frame):
        nonlocal stopped
        stopped = True
        if child is not None and child.poll() is None:
            child.send_signal(signal.SIGINT)

    signal.signal(signal.SIGINT, handle_stop)
    signal.signal(signal.SIGTERM, handle_stop)
    try:
        for k in manifest["order"]:
            if stopped:
                break
            status.update(current_k=k, child_pid=None)
            case = run / f"k{k}"
            case.mkdir(exist_ok=False)
            command = ["stdbuf", "-oL", "-eL", str(run / "source/gaston_3round_lossy"),
                       "--output", str(case / "solve"), "--threads", "16", "--kernel-both",
                       "--middle-active", str(k), "--weight-cap", "43", "--kernel-reduced"]
            print("START " + json.dumps(command, ensure_ascii=False), flush=True)
            with (case / "console.log").open("xb") as log:
                child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
                status["child_pid"] = child.pid
                save(run / "status.json", status)
                code = child.wait()
            if code:
                raise RuntimeError(f"Branch k={k} exited with code {code}")
            report = json.loads(Path(manifest["verification"][str(k)]).read_text())
            validated_mps = Path(manifest["verification"][str(k)]).parent / "model.mps"
            assert digest(case / "solve/model.mps") == report["sha256"][str(validated_mps)]
            result = json.loads((case / "solve/result.json").read_text())
            status["completed"][str(k)] = result
            status["pending"].remove(k)
            save(run / "status.json", status)
            print("FINISH k=" + str(k) + " " + json.dumps(result), flush=True)
            if result["status_code"] not in (2, 3):
                stopped = True
                break
        status.update(current_k=None, child_pid=None)
        if stopped:
            status["state"] = "stopped"
        elif all(result["subproblem_infeasible"] for result in status["completed"].values()):
            assert len(status["completed"]) == 3
            status["state"] = "all_three_infeasible_plus_analytic_middle9"
            status["conclusion"] = "Together with coverage and analytic exclusion of middle9, solver results establish lossy and true three-round weight >=44; no independent UNSAT certificate."
        else:
            status["state"] = "completed_with_lossy_feasible_cases"
            status["conclusion"] = "This lossy decomposition does not exclude all weights below44; DDT-check witnesses before treating them as true trails."
    except Exception as error:
        status.update(state="failed", error=str(error))
        raise
    finally:
        save(run / "status.json", status)


def launch(run):
    manifest = validate_manifest(run)
    assert manifest["state"] == "verified_not_running", "Already launched"
    command = [sys.executable, str(run / "source/src/kernel_campaign.py"), "worker", "--run", str(run)]
    with (run / "worker.log").open("xb") as log:
        process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
    manifest.update(state="launched", worker_pid=process.pid,
                    started_local=datetime.now().astimezone().isoformat(), command=command,
                    stop_instruction="Send SIGINT to worker_pid; it interrupts the current solve and does not start another branch.")
    save(run / "manifest.json", manifest)
    save(ROOT / "results/latest_kernel_campaign.json", manifest)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "launch", "worker"))
    parser.add_argument("--run", type=Path)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "launch":
        launch(args.run.resolve())
    else:
        worker(args.run.resolve())
