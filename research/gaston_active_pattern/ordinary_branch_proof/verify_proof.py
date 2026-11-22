#!/usr/bin/env python3
"""Verify the complete ordinary-branch proof bundle.

Default: exact matching/witness checks plus hashes and saved DRAT-check results.
--full: additionally rerun drat-trim on all 1667 proofs.
No Gurobi, SAT search, or python-sat dependency.
"""
import argparse
from collections import Counter
import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
MASK=(1<<64)-1
EAST,WEST,TWIN=(0,60,22,27,4),(0,56,31,46,43),(25,32,52,60,63)


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        while chunk:=f.read(1024*1024):h.update(chunk)
    return h.hexdigest()


def rol(w,d):
    d%=64
    return ((w<<d)|(w>>((64-d)%64)))&MASK


def weight(words):
    active=0
    for w in words:active|=w
    return active.bit_count()


def verify_upper():
    d=json.loads((ROOT/"upper_bound16.json").read_text())
    x=[int(w,16) for w in d["input_words"]]
    y=[int(w,16) for w in d["output_words"]]
    require(len(x)==len(y)==5 and all(0<=w<=MASK for w in x+y),"invalid witness words")
    a=[rol(w,t) for w,t in zip(x,EAST)]
    p=q=0
    for w,t in zip(a,TWIN):p^=w;q^=rol(w,t)
    e=p^rol(p,1)^q^rol(q,18)
    require(e!=0 and e==int(d["theta_effect_hex"],16),"upper witness not outside kernel")
    require([rol(w^rol(e,23),t) for w,t in zip(a,WEST)]==y,"word propagation failed")
    # Independent per-bit destinations, including cancellations.
    by_bits=[0]*5
    for row,w in enumerate(x):
        for bit in range(64):
            if not(w>>bit&1):continue
            pos=(bit+EAST[row])%64
            by_bits[row]^=1<<((pos+WEST[row])%64)
            for shift in (0,1,TWIN[row],TWIN[row]+18):
                for r in range(5):
                    by_bits[r]^=1<<((pos+shift+23+WEST[r])%64)
    require(by_bits==y,"bit propagation failed")
    require(weight(x)==d["input_active"]==3 and weight(y)==d["output_active"]==13,
            "wrong witness column counts")
    require(weight(x)+weight(y)==d["total_active"]==16,"wrong upper bound")
    return {"input_active":3,"output_active":13,"theta_effect_hex":f"0x{e:016x}"}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--full",action="store_true",help="recheck every DRAT proof")
    parser.add_argument("--jobs",type=int,default=4)
    args=parser.parse_args()
    upper=verify_upper()
    print("PASS: explicit outside-kernel witness has total activity 16.",flush=True)
    subprocess.run([sys.executable,str(ROOT/"verify_matchings.py")],check=True)
    matching=json.loads((ROOT/"matching_verification.json").read_text())
    require(matching["verification"]=="PASS" and matching["all_weight_at_least_eight_excluded"],
            "matching proof incomplete")
    base=(ROOT/"total15.cnf").read_text()
    mapping=json.loads((ROOT/"total15.map.json").read_text())
    parameters=mapping["parameters"]
    require(parameters["total"]==15 and parameters["input"] is None and parameters["output"] is None
            and parameters["effect"]=="any" and not parameters.get("reduced_domain",False)
            and parameters["effect_weight"] is None,"base model is restricted")
    cases=[]
    for line in (ROOT/"effects_total15.tsv").read_text().splitlines()[1:]:
        w,e,bound=line.split()
        cases.append((int(w),int(e,16),int(bound)))
    require(len(cases)==1667,"wrong case count")
    def check_case(case):
        w,e,bound=case
        prefix=ROOT/"effects_bound15"/f"e{w}_{e:016x}"
        d=json.loads(prefix.with_suffix(".json").read_text())
        require(d.get("drat_verified") is True and d.get("solver_exit_code")==20
                and d.get("checker_exit_code")==0,f"incomplete proof: {prefix.name}")
        require(int(d["effect"],16)==e and d["effect_weight"]==w and d["total_max"]==15,
                "case metadata mismatch")
        formula=base.replace(f"p cnf {mapping['variables']} {mapping['clauses']}\n",
                             f"p cnf {mapping['variables']} {mapping['clauses']+64}\n",1)
        formula+="".join(f"{v if e>>j&1 else -v} 0\n" for j,v in enumerate(mapping["E"]))
        require(hashlib.sha256(formula.encode()).hexdigest()==d["cnf_sha256"],"wrong case formula hash")
        require(gzip.decompress(prefix.with_suffix(".cnf.gz").read_bytes())==formula.encode(),
                "saved CNF differs from the unrestricted model plus E units")
        require(sha(prefix.with_suffix(".drat.gz"))==d["proof_sha256"],"DRAT proof hash mismatch")
        require("s VERIFIED" in prefix.with_suffix(".check.log").read_text(),
                "missing successful checker result")
        if args.full:
            # Each case has its own filenames; this does not alter saved certificates.
            cnf=prefix.with_suffix(".recheck.cnf")
            log=prefix.with_suffix(".recheck.log")
            cnf.write_text(formula)
            unzip=subprocess.Popen(["gzip","-dc",str(prefix.with_suffix(".drat.gz"))],stdout=subprocess.PIPE)
            with log.open("wb") as output:
                process=subprocess.run([str(ROOT/"tools/drat-trim/drat-trim"),str(cnf),"-i"],
                                       stdin=unzip.stdout,stdout=output,stderr=subprocess.STDOUT)
            unzip.stdout.close()
            code=unzip.wait()
            require(process.returncode==0 and code==0 and "s VERIFIED" in log.read_text(),
                    f"DRAT replay failed: {prefix.name}")
            cnf.unlink()
        return w
    counts=Counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1,args.jobs)) as executor:
        for w in executor.map(check_case,cases):counts[w]+=1
    require(counts==Counter({2:32,4:1632,6:3}),"missing cases")
    result={
        "verification":"PASS",
        "ordinary_branch_number_outside_theta_kernel":16,
        "lower_bound":16,"upper_bound":16,
        "upper_witness":upper,
        "weight_at_least_eight_excluded_by_matching_extensions":True,
        "retained_cases_with_checked_unsat_proofs":dict(counts),
        "drat_replayed_in_this_invocation":args.full,
        "initial_drat_checks_completed_for_every_case":True,
        "base_cnf_sha256":sha(ROOT/"total15.cnf"),
        "matching_certificates_sha256":matching["certificate_sha256"],
        "extension_certificates_sha256":matching["extension_certificate_sha256"],
        "statement":"Y=lambda(X), E(rho_east(X)) != 0, metric=nonzero 5-bit columns",
    }
    (ROOT/"proof_verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: all 1667 cases UNSAT with checked proofs. Exact outside-kernel ordinary branch number B=16.")


if __name__=="__main__":
    main()
