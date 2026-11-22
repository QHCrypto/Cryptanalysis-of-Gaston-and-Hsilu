# Code accompanying “Cryptanalysis of Gaston and Hsilu”

This is the lightweight source-code package for the manuscript by Xiaobin Yu,
Zhaoyue Huang, Haoran Li, Jiajun Zhou, and Meicheng Liu. It contains research
implementations, necessary small inputs, selected result records, and a tested
entry point for bounded verification. Correspondence: Meicheng Liu,
liumeicheng@iie.ac.cn.

**Journal**: Cryptographic Theory and Practice.

**Corresponding-author affiliations**: State Key Laboratory of Cyberspace
Security Defense, Institute of Information Engineering, Chinese Academy of
Sciences, Beijing, China; School of Cyber Security, University of Chinese
Academy of Sciences, Beijing, China.

**Public repository**: https://github.com/QHCrypto/Cryptanalysis-of-Gaston-and-Hsilu

**Supplementary material**: Online Resource 1 (`ESM_1.zip`). Source code,
verification scripts, necessary small inputs, and selected computational records
for the differential, differential-linear, preimage, and cube analyses of
Gaston and Hsilu. Machine-readable article information is in
`ARTICLE_METADATA.json`.

The package was assembled from the authors' local research directories on
2026-11-22. The original directory structure is preserved inside `research/`
because several programs share permutation implementations and algebra helpers.
Start with the commands and result map in this README; the historical READMEs
inside individual directories describe a broader research history.

**Large proof certificates, matching-certificate binaries, and regenerable
expanded polynomial/DAG files are omitted**. This package is not a complete
archive for replaying every proof certificate or every exhaustive search.

Chinese guide: [README_zh.md](README_zh.md).

## Quick start

Use Python 3.10 or newer, NumPy, GCC/G++ with C++17 support, `make`, and the GMP
development libraries. The tested environment is recorded under `validation/`.
Gurobi is not needed for either verification command below.

```bash
# Run from this directory.
python3 run_checks.py --profile quick
python3 run_checks.py --profile datf
```

The runner creates a temporary working copy, builds the needed native programs,
and writes logs and a JSON summary under `validation/<timestamp>_<profile>/`.
The archived research inputs remain unchanged. Every selected group must print
`PASS`; the process returns a nonzero exit status if a command fails.

The `quick` profile checks:

1. The current manuscript's permutation parameters, native rounds, and XOF/AEAD
   conventions against the included implementations.
2. The weight-80 Hsilu characteristic, local DDT identities, the Gaston
   weight-32 witness, and selected preimage affine relations.
3. The saved Gaston support-enumeration records and independent nullspace
   calculations used in the weight-44 argument. This does not repeat the
   exhaustive support enumeration or establish the historical activity-profile
   optimality results anew.
4. Differential-linear implementation tests, exact local identities, and a
   seeded Monte Carlo check of the two four-round endpoints, with `2^22`
   samples per permutation.
5. The current three-round Hsilu prefix and linearization conditions, followed
   by all `2^23` branches of the concrete saved algebraic system. This finite
   branch enumeration is not a complete 128-bit preimage search.
6. Fresh four-round Gaston cube sums, rank reconstruction from the saved
   five-round sums, reduced-width quartic checks, and reduced-instance Hsilu
   batch-evaluation checks.
7. The three preimage complexity calculations under the manuscript's Huang
   precomputation and Gaussian-elimination convention.

The separate `datf` profile recomputes DATF Algorithm 5 for the manuscript's
four-round endpoints and compares the exact rational model outputs with the
archived values. The historical endpoint identifiers `G4a` and `H4a` correspond
to the selected Gaston and Hsilu four-round examples. The model gives `3/2048`
and `0`, respectively; the latter is a model output, not a claim that the actual
Hsilu correlation is zero.

To rerun a particular group:

```bash
python3 run_checks.py --profile quick --only preimage_structure,preimage_accounting
```

## Manuscript-to-code map

Paths in this table are relative to `research/`.

| Manuscript component | Main sources and inputs | What is verified or estimated |
|---|---|---|
| Permutations and outer modes | `hsilu_ref/`; `gaston_ref/gaston/implementation/reference/`; `HsiluCube/src/permutations.py`; `nist_compatibility_20261109/` | Native round order, constants, bit conventions, and mode interfaces |
| Gaston differential bounds | `gaston_active_pattern/`; `Gaston_差分路径重量下界/` | Outside-kernel branch number 16, two-round weight 32, and the three-round weight-44 argument; full proof-certificate files are omitted |
| Hsilu weight-80 trail | `Hsilu_trail/search_below82_20261117/` | Explicit trail with round weights 16, 14, 50; this is an upper-bound witness for the minimum trail weight, not a proof of global optimality |
| Differential-linear analysis | `差分线性/gaston_dl/`; `差分线性/hsilu_dl/`; `差分线性/路径分段汇总_20261027/`; `Gaston_Hsilu_CTP/paper_rewriting_output/evaluate_datf_algorithm5.py` | Local exact identities, DATF/HNW model estimates, endpoint selection, and sampled correlations; concatenation estimates retain the manuscript's assumptions |
| Two-round Gaston preimages | `gaston分析/improved_milp.cpp`; `gaston分析/model_common.hpp`; `Gaston_Hsilu_CTP/paper_rewriting_output/ctp_reviewer_20261117/independent_checks.py` | Shared-form strategy, rank 47, and the 17 selected affine relations |
| Two-round Hsilu preimages | `代数攻击/hsilu_2r/`; `Gaston_Hsilu_CTP/paper_rewriting_output/ctp_reviewer_20261117/independent_checks.py` | Current 22-variable construction and selected output relations; the current figure input is `hsilu_no_start.svg` |
| Three-round Hsilu preimages | `代数攻击/Hsilu3r_new/validation_current_prefix/`; `代数攻击/Hsilu3r_new/hsilu_3round_best_attack.py` | Current legal prefix; 32-condition rank; 25 coordinate guesses of rank 23; concrete finite branch and background counts |
| Gaston cube construction | `GastonCube/`; `Cube/gaston_rank_certificate/` | Even-parity cube, rotation reuse, and low-round rank checks. Four-round rank is 128 and saved five-round rank is 125; the six-round rank hypothesis is not verified by these checks |
| Hsilu cube construction | `Cube/hsilu_quartic_screen/`; `Cube/src/hsilu_7r_batch_audit.py` | Four-factor terms, reduced-instance checks, and counted full-size batch schedules. Retained terms alone do not prove a nonzero coefficient |

`research/Gaston_Hsilu_CTP/verification_inputs/paper_parameters.tex` contains the
parameter tables, round constants, mode IVs, and differential-trail witness used
by the verification scripts. These numerical inputs make the code package
self-contained. The complete manuscript and its PDF belong to the separate
submission package.

The current preimage accounting is in `scripts/preimage_accounting.py`:

| Construction | Gaussian-elimination cost per candidate, in permutation units | Reported search exponent |
|---|---:|---:|
| Gaston, 2 rounds | 0.3984375 | 118.05563749, displayed as 118.06 |
| Hsilu, 2 rounds | 0.954202586 | 106 |
| Hsilu, 3 rounds | 0.174075248 | 126 |

These are model estimates under the manuscript's allocation, target-coverage,
and dominant-term conventions. They are not measured full-search running times.
The Hsilu seven-round time `2^102.22` and space `2^108.44` bits are also counted
estimates under the stated table-access and candidate-bucket assumptions. The
full six-round Gaston cube and full seven-round Hsilu recovery are not executed
by the verification profiles.

## Search programs and optional dependencies

The bounded verification entry points above use the Python standard library,
NumPy, native compilers, and GMP. To install the Python dependency in an existing
environment:

```bash
python3 -m pip install -r requirements.txt
```

MILP searches additionally need a local Gurobi installation and a valid license.
The copied C++ Makefiles default to Gurobi 13; override `GUROBI_HOME` for your
installation, and adjust the versioned linker name if using another release.
The three-round Hsilu Python search uses `gurobipy`. Optional SAT searches in
`Hsilu_trail/pysat_kissat/` also require `python-sat` and the configured SAT solver.
No solver binaries, license files, or Python environments are bundled.

Examples of source build entry points:

```bash
make -C research/Hsilu_trail/search_below82_20261117
make -C research/GastonCube
make -C research/Cube/hsilu_quartic_screen
```

The first command builds the standalone C++ trail search; the other two build
Gurobi-dependent programs. The Makefiles and source argument parsers document
the individual search options. Searches may take substantially longer than the
bounded checks, and a resource limit must not be interpreted as infeasibility.

## Provenance, portability, and integrity

- `SOURCE_MANIFEST.json` records each copied research file's original path and
  original SHA-256. Original absolute paths appearing in saved JSON files are
  historical provenance; the tested entry points resolve their inputs within
  this package.
- `PORTABILITY_CHANGES.json` lists changes to copied scripts: relative input
  resolution, current-manuscript parsing, and removal of superseded accounting
  outputs from an old verification helper. Research originals were not changed.
- `OMITTED_LARGE_ARTIFACTS.json` records the large files excluded from this
  lightweight distribution. Sources for the omitted proof computations remain
  included, but proof replay requires the omitted data or its regeneration.
- `SHA256SUMS` hashes the distributed sources, inputs, and documentation;
  generated `validation/` reports are excluded. Check with `sha256sum -c
  SHA256SUMS` from the package root.
- Existing third-party source headers and license files are retained. See
  [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Current result interpretation follows this README and the manuscript. Older
candidate files and directory names are retained where they are inputs to shared
helpers; they should not be substituted for the current constructions above.
