# Third-party material and source attribution

The original Gaston and Hsilu reference implementations are retained under
`research/gaston_ref/` and `research/hsilu_ref/`. Their existing source headers
and available documentation are preserved.

The DATF reference material in
`research/差分线性/gaston_dl/vendor/datf_reference/` retains its upstream
`LICENSE`, `README.md`, and `source.json`.

The HNW26 round-based evaluation follows Hu, Niu, and Wang,
[Round-Based Approximation of (Higher-Order) Differential-Linear Correlation](https://eprint.iacr.org/2026/358).
For G4/H4, the local implementation combines their two-round differential-prefix
weighting with the joint final-layer treatment in their public
[`PRESENT/present.py`](https://github.com/hukaisdu/Round-based-approximation/blob/783e8edc560b848d0f2ba63f11838aecd23da40f/PRESENT/present.py),
functions `getBias_Opt2` and `getBias_Opt4`, commit
`783e8edc560b848d0f2ba63f11838aecd23da40f`. The S-box and mixer are adapted to
Gaston/Hsilu; affine elimination evaluates the prefix sum. The joint final-layer
treatment is not claimed as a new method of this manuscript.

The SAT-related code under `research/Hsilu_trail/pysat_kissat/` retains the
available `LICENSE` and provenance documentation. The DRAT checker source in
`research/gaston_active_pattern/ordinary_branch_proof/tools/drat-trim/` retains
its own license. The large proof files used by that checker are not included.

The Ascon quartic list under `research/Cube/hsilu_quartic_screen/` is accompanied
by a comparison record identifying the source repository and commit:
https://github.com/keyrecoveryofascon/keyrecoveryofascon,
`616f5d247bd9decdb39a74f2339f73a115392910`.

This collection preserves existing notices and does not assign a new blanket
license to the research code or third-party material. Solver dependencies are
obtained separately from their respective providers.
