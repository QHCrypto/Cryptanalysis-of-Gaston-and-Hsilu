#!/usr/bin/env python3
"""Launch the verified C++ model in a persistent session; no timeout or stop condition."""
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import uuid

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verification", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=16, choices=range(1, 17))
    args = parser.parse_args()
    report = json.loads(args.verification.read_text())
    if report["status"] != "PASS":
        raise RuntimeError("Model verification must pass before launch")
    if any(report.get("subproblem", {}).values()) or report.get("kernel_reduced"):
        raise RuntimeError("Use kernel_campaign.py for restricted subproblems; this launcher runs the unrestricted model")
    for path in [ROOT/"gaston_3round_lossy.cpp", ROOT/"gaston_3round_lossy",
                 ROOT/"src/activity_profile_data.hpp", ROOT/"references/activity_profiles.json",
                 ROOT/"src/kernel_type_data.hpp", ROOT/"results/kernel_structure_audit_20261024.json"]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != report["sha256"][str(path)]:
            raise RuntimeError("Source/binary changed after verification: " + str(path))
    run = ROOT / "results/runs" / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8])
    run.mkdir(parents=True, exist_ok=False)
    snapshot = run/"source"
    snapshot.mkdir()
    for source in [ROOT/"gaston_3round_lossy.cpp", ROOT/"Makefile", Path(__file__),
                   ROOT/"src/verify_lossy_model.py", ROOT/"src/activity_profile_data.hpp",
                   ROOT/"src/build_activity_profiles.py", ROOT/"references/activity_profiles.json",
                   ROOT/"src/kernel_type_data.hpp", ROOT/"src/build_kernel_types.py",
                   ROOT/"results/kernel_structure_audit_20261024.json"]:
        destination = snapshot/source.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    shutil.copy2(args.verification, snapshot/"verification.json")
    command = ["stdbuf", "-oL", "-eL", str(ROOT/"gaston_3round_lossy"),
               "--output", str(run/"solve"), "--threads", str(args.threads)]
    # Preserve every Gurobi console line. No Gurobi LogFile/OutputFlag suppression.
    # start_new_session detaches job lifetime from the interactive tool session.
    with (run/"console.log").open("xb") as logfile:
        process = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.DEVNULL,
                                   stdout=logfile, stderr=subprocess.STDOUT,
                                   start_new_session=True)
    metadata = {
        "pid": process.pid,
        "started_local": datetime.now().astimezone().isoformat(),
        "command": command,
        "cwd": str(ROOT),
        "threads": args.threads,
        "time_limit": None,
        "external_timeout": None,
        "run_directory": str(run),
        "console_log": str(run/"console.log"),
        "verified_model_mps_sha256": report["sha256"][str(args.verification.resolve().parent/"model.mps")],
        "source_sha256": report["sha256"][str(ROOT/"gaston_3round_lossy.cpp")],
        "binary_sha256": report["sha256"][str(ROOT/"gaston_3round_lossy")],
        "stop_instruction": "Send SIGINT to the recorded PID to request a graceful stop and save an incumbent.",
    }
    (run/"launch.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2)+"\n")
    (ROOT/"results/latest_run.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
