#!/usr/bin/env python3
"""Check all matching certificates and full rotation-orbit coverage.

Python standard library only. Does not trust the enumerator's matching algorithm.
"""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import struct
import time

ROOT = Path(__file__).resolve().parent
MASK = (1 << 64)-1
EAST, WEST = (0,60,22,27,4), (0,56,31,46,43)


def require(value, message):
    if not value:
        raise ValueError(message)


def orbit_count(k):
    return sum(math.comb(math.gcd(64, t), k//(64//math.gcd(64,t)))
               for t in range(64) if k % (64//math.gcd(64,t)) == 0)//64


def check_matching(e, matching, bits=None):
    if bits is None:
        bits = [j for j in range(64) if e >> j & 1]
    valid_edges = {((j+23-east) % 64) | (((j+23+west) % 64) << 6)
                   for j in bits for east, west in zip(EAST,WEST)}
    aa = bb = 0
    for edge in matching:
        require(edge in valid_edges, "matching edge absent from forced-activity graph")
        a, b = edge & 63, edge >> 6
        require(not (aa >> a & 1) and not (bb >> b & 1),
                "matching edges do not have distinct endpoints")
        aa |= 1 << a
        bb |= 1 << b
    require(aa.bit_count() == bb.bit_count() == 16, "not a size-16 matching")


def main():
    started = time.monotonic()
    seen = set()
    retained = Counter()
    pruned = Counter()
    retained_six = set()
    def register(e):
        k = e.bit_count()
        require(k in (2,4,6) and e & 1, "invalid canonical effect")
        require(e not in seen, "duplicate rotation representative")
        seen.add(e)
        bits = []
        v = e
        while v:
            low = v & -v
            j = low.bit_length()-1
            v ^= low
            require(((e >> j) | (e << (64-j))) & MASK >= e, "noncanonical rotation")
            bits.append(j)
        return k, bits
    for line in (ROOT / "effects_total15.tsv").read_text().splitlines()[1:]:
        weight, hexword, bound = line.split()
        e = int(hexword,16)
        k, _ = register(e)
        require(k == int(weight) and int(bound) <= 15, "invalid retained case")
        retained[k] += 1
        if k == 6:
            retained_six.add(e)
    record = struct.Struct("<Q16H")
    with (ROOT / "pruned_matchings.bin").open("rb") as file:
        require(file.read(8) == b"GSTNMAT1", "invalid certificate header")
        done = 0
        while chunk := file.read(record.size*10000):
            require(len(chunk) % record.size == 0, "truncated matching certificate")
            for e, *matching in record.iter_unpack(chunk):
                k, bits = register(e)
                check_matching(e, matching, bits)
                pruned[k] += 1
                done += 1
            if done % 200000 == 0:
                print(f"Checked {done} size-16 matchings in {time.monotonic()-started:.1f}s",
                      flush=True)
    for k in (2,4,6):
        require(retained[k]+pruned[k] == orbit_count(k),
                f"incomplete orbit coverage at weight {k}")
    extension_record = struct.Struct("<QQ16H")
    extensions = {base:set() for base in retained_six}
    with (ROOT/"extension_matchings.bin").open("rb") as file:
        require(file.read(8) == b"GSTNEXT1", "invalid extension header")
        data = file.read()
        require(len(data)%extension_record.size == 0, "truncated extension certificate")
        for base, e, *matching in extension_record.iter_unpack(data):
            require(base in retained_six and (base & e) == base and e.bit_count() == 8,
                    "invalid extension")
            require(e not in extensions[base], "duplicate extension")
            extensions[base].add(e)
            check_matching(e,matching)
    require(len(retained_six) == 3 and all(len(es)==math.comb(58,2) for es in extensions.values()),
            "incomplete coverage of two-bit extensions")
    result = {
        "verification": "PASS",
        "retained_by_weight": dict(retained),
        "pruned_by_weight": {k:pruned[k] for k in (2,4,6)},
        "orbit_counts": {k:orbit_count(k) for k in (2,4,6)},
        "size16_matchings_checked": sum(pruned.values()),
        "all_orbits_covered_exactly_once": True,
        "all_weight_at_least_eight_excluded": True,
        "eight_bit_extensions_checked": sum(map(len,extensions.values())),
        "extension_certificate_sha256": hashlib.sha256((ROOT/"extension_matchings.bin").read_bytes()).hexdigest(),
        "certificate_sha256": hashlib.sha256((ROOT/"pruned_matchings.bin").read_bytes()).hexdigest(),
        "retained_tsv_sha256": hashlib.sha256((ROOT/"effects_total15.tsv").read_bytes()).hexdigest(),
    }
    (ROOT/"matching_verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2), flush=True)


if __name__ == "__main__":
    main()
