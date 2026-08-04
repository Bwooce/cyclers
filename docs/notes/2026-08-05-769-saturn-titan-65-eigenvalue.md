# `#769`: Saturn-Titan 6:5 eigenvalue near-miss — follow-up on `#765`'s honest FAIL

**Task:** `#769`, LOW-PRIORITY, exploratory. `#765` (`src/cyclerfinder/search/saturn_titan_resonant_families.py`)
reproduces Vaquero 2013 Table 4.1's four Saturn-Titan resonant/Lyapunov orbits at `C=3.010000`;
three of four pass the `1e-3` relative-error gate on both eigenvalue and period, but the 6:5
resonant orbit's eigenvalue misses by `2.34e-3` (191.1928 recovered vs. 191.641 target) despite
its own IC and period matching to <0.02% relative — as tight as the passing rows. `#765`'s own
docstring proposed a mechanism: a fixed `µ`-precision floor (from Vaquero's own 5-significant-
figure, "≈"-qualified `µ` display) producing a proportionally larger *relative* eigenvalue error
for 6:5 specifically because its eigenvalue magnitude (191.6) is an order of magnitude smaller
than 3:4's (2129.8). `#769` was dispatched to try to close this gap, or else confirm/refute that
mechanism.

## Path 1: search the thesis for a more precise `µ` — clean negative, but confirms the module already uses the best available value

Grepped the OCR text sidecar
(`cyclers_pdf/papers/vaquero-2013-spacecraft-transfer-trajectory-design-resonant-orbits-multibody-environments-purdue-phd.txt`)
for every occurrence of `µ`/"mass parameter"/"mass fraction" and for the digit string `2.36`
(Saturn-Titan `µ`'s leading digits) across the whole document — background/Ch.2 formulas,
Ch.4.3 Saturn-Titan intro, and any appendix. The ONLY place a numeric Saturn-Titan `µ` value
appears anywhere in the thesis is p.132 (in the module's own docstring citation): "the mass
parameter is larger in the Earth-Moon system (µ ≈ 0.0122) than in the Saturn-Titan (µ ≈ 2.3658 ×
10⁻⁴) and Jupiter-Europa (µ ≈ 2.5266 × 10⁻⁵) systems." `VAQUERO_MU = 2.3658e-4` in the module is
already exactly this value. No higher-precision, more-decimal-place, or alternatively-sourced
`µ` exists anywhere else in the thesis (no separate physical-constants table, no appendix). This
path is a clean negative — not because the search was incomplete, but because there is nothing
more precise to find.

## Path 2: diagnostic sensitivity fit — decisive, and it overturns `#765`'s own proposed mechanism

Ran a bisection on `µ` (scratch script, not committed) to find the value that makes 6:5's
recovered eigenvalue hit Table 4.1's printed `191.641` exactly, then checked what that `µ` does
to the other three rows:

| `µ` | 3:4 rel_err | 6:5 rel_err | L1 rel_err | L2 rel_err |
|---|---|---|---|---|
| `2.3658e-4` (thesis, baseline) | 1.07e-6 PASS | 2.34e-3 FAIL | 4.56e-6 PASS | 2.73e-6 PASS |
| `2.36068e-4` (bisected to fix 6:5, Δ=-0.216% rel) | 2.09e-3 **FAIL** | 2.46e-12 PASS | 2.03e-3 **FAIL** | 2.08e-3 **FAIL** |

**No single `µ` value reconciles all four rows** — the `µ` that perfectly fixes 6:5 breaks the
other three, pushing them from near-machine-precision (1e-6-level) matches to ~2e-3, *worse*
than 6:5's own original miss. This alone falsifies a simple "global `µ` imprecision" explanation.

Measuring the actual sensitivity mechanism made this precise. A uniform +0.1% `µ` perturbation
(the same probe `#765`'s own test suite already runs,
`test_eigenvalue_sensitivity_to_mu_is_measured_not_assumed`) shifts each row's eigenvalue by:

| Row | baseline λ | perturbed λ (µ×1.001) | absolute shift | **relative shift** |
|---|---|---|---|---|
| 3:4 | 2129.8077 | 2131.8530 | +2.0453 | **+0.0960%** |
| 6:5 | 191.1928 | 190.9851 | -0.2077 | **-0.1086%** |
| L1 | 1004.7246 | 1003.7882 | -0.9363 | **-0.0932%** |
| L2 | 892.8524 | 891.9966 | -0.8558 | **-0.0959%** |

The **relative** shift is comparable across all four rows (0.093%-0.109%), NOT the **absolute**
shift as `#765`'s docstring originally claimed (absolute shift instead scales with each row's
own eigenvalue magnitude — 6:5's absolute shift, 0.21, is in fact the SMALLEST of the four, not
the largest). Note also 3:4 shifts in the OPPOSITE direction from 6:5/L1/L2. Because relative
`µ`-sensitivity is essentially the same for every row, a `µ`-precision floor would show up as a
comparable RELATIVE eigenvalue error in ALL FOUR rows — not the observed pattern of three rows
at 1e-6-level and one row at 2.3e-3. **`µ` imprecision is excluded as 6:5's root cause.**

### The Jacobi constant `C` is a much better-supported explanation

`C` is the other display-rounded input (Table 4.1 prints `C = 3.010000`, though — unlike `µ` —
*without* Vaquero's own "≈" qualifier). A parallel bisection on `C` (same scratch-script
approach) to hit 6:5's target eigenvalue found:

* Required shift: `C_fit = 3.0099955250`, `ΔC = -4.475e-6` absolute, **`1.487e-6` relative** —
  about 145,000× smaller than the `µ` shift that was needed and failed.
* At `C_fit`, **all four rows land inside the `1e-3` gate simultaneously**:

  | Row | eigenvalue rel_err at `C_fit` | Gate |
  |---|---|---|
  | 3:4 | 2.12e-6 | PASS |
  | 6:5 | 2.39e-12 | PASS |
  | L1 | 6.66e-4 | PASS |
  | L2 | 6.59e-4 | PASS |

* This `ΔC` (1.5e-6 relative) is the same order of magnitude as this module's own
  ALREADY-MEASURED C self-validation residual from `#765`'s own
  `test_table41_ic_reproduces_stated_jacobi_constant` (nondimensionalizing each row's own
  printed IC at this module's derived `l*`/`t*` and the thesis's stated `µ` reproduces
  `C = 3.010000` to 1.05e-6-5.20e-6 relative per row, 6:5's own residual being 1.95e-6) — i.e.
  well within this module's own already-documented `l*`/`t*`/C derivation precision floor, not
  a large or implausible correction.
* A direct sensitivity measurement (fixed IC, `+1e-5` relative `C` perturbation) confirms WHY:
  6:5's eigenvalue is disproportionately `C`-sensitive relative to the other three rows —
  measured coefficients (relative eigenvalue shift per relative `C` shift) of **~1569 for 6:5**
  vs. only **~3 (3:4), ~452 (L1), ~446 (L2)**. A ~350-500× difference in `C`-sensitivity between
  6:5 and the other rows means the SAME tiny `C`/`l*`/`t*`-derivation-scale imprecision that is
  invisible for 3:4/L1/L2 manifests as a materially larger relative eigenvalue error for 6:5
  specifically. This is physically plausible: unstable resonant-orbit eigenvalues are generally
  far more energy-sensitive than mass-ratio-sensitive, and 6:5's own eigenvalue-vs-C curve
  evidently sits in a steeper region than the other three rows' at this energy.

## Verdict: genuinely better-explained honest near-miss, NOT adopted as a fix

This is reported as a **diagnostic finding, not a correction**. Table 4.1 prints `C = 3.010000`
unqualified — unlike `µ`'s own "≈"-qualified display, there is no textual signal in the thesis
that this value is itself a rounded display of a higher-precision internal quantity, so there is
no sourced basis for adopting a different `C`. Per this project's sourced-only discipline,
`VAQUERO_C` and `VAQUERO_MU` are both left EXACTLY as printed; the `1e-3` eigenvalue gate is
NOT weakened; `6:5` continues to honestly report `passed=False` on the eigenvalue criterion,
exactly as `#765` reported it.

What changed as a result of this task: `#765`'s own module docstring "HONEST FINDINGS" item 2
proposed a specific causal mechanism (fixed absolute `µ`-precision floor) for 6:5's near-miss.
This task's own measurement empirically contradicts that mechanism (relative, not absolute,
`µ`-sensitivity is comparable across rows; the required `µ` correction is 3 orders of magnitude
larger than the sensitivity probe and breaks the other three rows). Per this project's own
citation/claim-grounding discipline (never leave a falsified causal claim sitting in a docstring
once directly measured evidence contradicts it), `#765`'s docstring finding 2 was corrected in
place to report both the exclusion of `µ` and the better-supported `C`-sensitivity explanation,
with a pointer to this note. A regression-guard assertion was also added to
`test_eigenvalue_sensitivity_to_mu_is_measured_not_assumed` (checks the four rows' relative
`µ`-sensitivities stay within a generous 5× spread of each other) so a future change that
reintroduces row-specific `µ`-sensitivity asymmetry — which would revive the originally-proposed
but now-excluded mechanism — trips a test rather than silently reappearing only in prose.

## Files touched

* `src/cyclerfinder/search/saturn_titan_resonant_families.py` — docstring-only correction
  (HONEST FINDINGS item 2); no constants, gate tolerances, or seeds changed.
* `tests/search/test_saturn_titan_resonant_families.py` — corrected docstring on
  `test_eigenvalue_sensitivity_to_mu_is_measured_not_assumed` + added a regression-guard
  assertion on cross-row relative-sensitivity spread. All 27 tests still pass (26 pre-existing +
  the strengthened one; no new test count change, same test just asserts more).

## Verification

* `uv run pytest tests/search/test_saturn_titan_resonant_families.py -v`: 27/27 pass.
* `uv run pytest tests/data tests/search -q`: green except two PRE-EXISTING, documented,
  unrelated Mac-local BLAS-sensitivity failures (`test_eggie_ballistic.py::test_gate_b_table4_vinf_reached_but_subsurface`,
  `test_504_pluto_charon_kk_sweep.py::test_504_sweep_33` — see `OUTSTANDING.md`'s many prior
  confirmations of this exact pair as green-on-Linux-CI, e.g. around lines 443-445, 5882-5899,
  9182-9183, 10517-10518, 10845-10846, 11066-11067, 11281-11282, 11379, 11666-11667). Neither
  test touches this module or its dependencies.
* `uv run ruff check` / `ruff format --check` on both touched files: clean.
* `uv run mypy src tests`: `Success: no issues found in 829 source files`.

## Bottom line

6:5's eigenvalue near-miss remains an honest FAIL under `#765`'s `1e-3` gate — this task did not
close it, and per its own LOW-PRIORITY/exploratory framing that is an acceptable outcome. But
the underlying explanation is now considerably better-supported: `µ` imprecision is empirically
EXCLUDED (not just "assumed not to be the whole story"), and a `C`-sensitivity asymmetry
(6:5's eigenvalue ~350-500× more `C`-sensitive than the other three rows, with the reconciling
`ΔC` sitting inside this module's own already-measured C self-validation floor) is a plausible,
evidence-backed diagnostic account of why 6:5 alone shows a large relative eigenvalue error
despite matching its own IC and period as tightly as the passing rows.
