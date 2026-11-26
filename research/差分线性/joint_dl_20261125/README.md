# Current four-round differential-linear evaluation

This directory reproduces the revised manuscript's G4 and H4 estimates with
two-round differential-prefix weighting. The endpoints are specified before
the first chi layer and after the last chi layer; the outer invertible linear
layers are absorbed into the endpoint difference or mask.

| Target | DATF Algorithm 3 correlation | Round-based correlation | Corresponding bias |
|---|---:|---:|---:|
| Gaston G4 | 1/256 | 1/256 | 1/512 |
| Hsilu H4 | 1/512 | 1/512 | 1/1024 |

These are exact rational evaluations of the stated approximation models.
The connection uses a fresh uniform base state. The calculation does not prove
the true correlation of a fixed four-round permutation.

## Complete prefix weighting

For each endpoint, enumerate all four first-round DDT output differences,
each with probability 1/4. Conditional on one such branch, the complete
second-round difference distribution under the local-transition product model
is a uniform affine image. Gaussian elimination integrates this image without
enumerating its individual elements or discarding dependencies between bits.
The image dimensions are 39, 73, 75, 103 for Gaston and 20, 18, 37, 37 for Hsilu.

For every prefix difference, the last chi layer is expanded over its complete
nonzero squared-LAT support. Each full mask is propagated by the transpose of
the mixer before evaluating the preceding chi layer. There are 16 last-mask
combinations and four first-round branches, giving 64 weighted terms per
endpoint. The only nonzero contributions are:

- G4: `(1/4) * (1/16) * (1/4) = 1/256`.
- H4: `(1/4) * (1/16) * (1/8) = 1/512`.

`datf_two_round_prefix.py` computes its local tables with the original DATF
Algorithm 3 backend. Fully expanding its two-round beta definitions and first
averaging the independent final alpha variables gives the same squared-LAT
formula. The formula is checked against direct two-round Algorithm 3 calls.

`recompute_for_manuscript.py` explicitly selects `prefix='fresh'` and
`split_output=False` in the round-based evaluator. Its readout uses the complete
transposed mask. `validate_hnw_fullmask.py` checks this readout against the
existing C++ nonlinear recurrence. Agreement between the two approaches here
does not imply equivalence for arbitrary round counts or model settings.

## Reproduction

From the repository root:

```bash
python3 run_checks.py --profile datf
```

This uses an isolated temporary copy and retains outputs in
`validation/<timestamp>_datf/`. It needs Python 3.10+, a C++17 compiler, and GMP;
Gurobi is not required. The checks cover:

- 600 affine-character calculations against direct enumeration of small systems;
  both complete linear transposes; quadratic derivatives and native C++ cores.
- 2,048 local Algorithm 3 truth-table comparisons and 121 two-round tail
  comparisons per permutation, including nonzero cases and full-state DAG checks.
- 49 C++ round-based tail comparisons per permutation.
- All seven selected manuscript rows, including freshly evaluated short DATF
  pairs and both five-round middle parts. Historical measurements and saved
  short round-based values retain their source records.

The main standalone commands, run from this directory, are:

```bash
python3 validate.py
python3 datf_two_round_prefix.py
python3 validate_hnw_fullmask.py
python3 recompute_for_manuscript.py
```

The last command reads the validation report produced by the preceding command.
The aggregate result builder is `scripts/current_dl_results.py` at the package
root. The runner places its fresh `CURRENT_DL_RESULTS.json` beside the logs.

## Relation to historical records

The original full four-round Algorithm 5 values, `3/2048` for G4 and `0` for
H4, are retained under `--profile datf-legacy`. Earlier search splits describe
endpoint selection; current G4/H4 estimates evaluate the complete endpoints
with the prefix weighting above. The current H4 has one Overall row, with no
additional concatenation row or middle mask. G5 and H5 remain concatenation
estimates with splits `(1,2,2)` and `(0,2,3)`.

当前稿件采用首两轮差分的完整加权求和；各首轮分支下的第二轮差分分布通过仿射像消元求和，
不是逐个遍历多达 `2^103` 个差分。G4、H4 的模型相关度分别为 `2^-8`、`2^-9`。
实验计数沿用原记录，与模型值分开保存；旧 Algorithm 5 结果仅作为历史对照。
