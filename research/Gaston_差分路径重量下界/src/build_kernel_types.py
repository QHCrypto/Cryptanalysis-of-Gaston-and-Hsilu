#!/usr/bin/env python3
"""Export the audited necessary low-weight count/bit types for the reduced C++ model."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "results/kernel_structure_audit_20261024.json"
audit = json.loads(source.read_text())
assert audit["status"] == "PASS"
assert audit["profile_sha256"] == hashlib.sha256((ROOT / "references/activity_profiles.json").read_bytes()).hexdigest()
types = []
for k in (6, 7, 8):
    types.extend(audit["branch_results"][str(k)]["remaining_types_under_44"])
assert len(types) == 67
header = ["#pragma once", "#include <array>",
          "struct KernelType { int c0,c1,c2,h1,h2,weight_floor; };",
          "constexpr std::array<KernelType,67> KERNEL_TYPES{{"]
for item in types:
    values = item["counts"] + item["internal_bit_weights"] + [item["weight_lower_bound"]]
    header.append("    {" + ",".join(map(str, values)) + "},")
header.extend(["}};", 'constexpr const char* KERNEL_AUDIT_SHA256="' + hashlib.sha256(source.read_bytes()).hexdigest() + '";'])
path = ROOT / "src/kernel_type_data.hpp"
path.write_text("\n".join(header) + "\n")
print(path)
