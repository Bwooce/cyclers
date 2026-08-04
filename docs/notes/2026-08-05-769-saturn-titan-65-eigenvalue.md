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

### A single global `C` offset does NOT explain the table either — adversarially reviewed and REFUTED

**This section originally claimed `C` was "a much better-supported explanation" than `µ`. A
Fable adversarial review (2026-08-05, dispatched at the user's explicit request — "review 769
with Fable, this is inconsistent and should be drilled into") independently re-derived every
number below and found that claim overclaimed on three separate counts.** The bottom-line
verdict (6:5 stays an honest FAIL; nothing adopted) is unaffected, but the diagnostic story is
corrected here rather than left standing.

`C` is the other display-rounded input (Table 4.1 prints `C = 3.010000`, though — unlike `µ` —
*without* Vaquero's own "≈" qualifier). A parallel bisection on `C` (same scratch-script
approach) to hit 6:5's target eigenvalue found `C_fit = 3.0099955250`, `ΔC = -4.475e-6` absolute,
`1.487e-6` relative — and at `C_fit` all four rows' eigenvalue errors are: 3:4 `2.12e-6`, 6:5
`2.39e-12`, L1 `6.66e-4`, L2 `6.59e-4` — all technically inside the `1e-3` gate. The original
draft presented this as support for `C`. Three things wrong with that:

1. **"All four inside the gate" is a threshold artifact, not evidence FOR a shared-`C`
   explanation — and by the SAME standard used to exclude `µ`, it's evidence AGAINST one.**
   L1/L2's own eigenvalue errors get **146×/242× WORSE** at `C_fit` (from their baseline
   `4.6e-6`/`2.7e-6` up to `6.66e-4`/`6.59e-4`) — the same qualitative degradation pattern that
   was used to declare `µ` excluded (there, the other three rows degraded ~2000×, to `~2.1e-3`).
   Under a genuine "Vaquero's own true `C` was `C_fit`" hypothesis, L1/L2 should IMPROVE toward
   machine precision, not degrade by two orders of magnitude while staying just inside a loose
   gate.
2. **Inverting each row's own baseline error through its own measured `C`-sensitivity makes this
   precise.** 6:5's baseline error (`-2.339e-3`) implies an effective `ΔC/C` of `-1.49e-6`; L1's
   baseline error (`+4.561e-6`) implies only `-1.0e-8`; L2's (`+2.728e-6`) implies only `-6.1e-9`.
   **L1 and L2 pin Vaquero's effective `C` to the printed `3.010000` to ~1e-8 relative — about
   150× tighter than what 6:5 alone would require.** A single shared `C` offset is excluded by
   this cross-row inconsistency, not supported by it.
3. **The "~350-500× more C-sensitive" claim (used to argue the effect is "invisible for the
   others") is a factual error.** The measured coefficients are `~1569` (6:5), `~452` (L1),
   `~446` (L2), `~3` (3:4, itself ill-determined — a ~45% spread between the `±` probe
   directions). `1569/452 ≈ 3.5×` and `1569/446 ≈ 3.5×` — the "~350-500×" figure only holds
   against 3:4 specifically, not against L1/L2, which is exactly where finding (2) above bites:
   at only a ~3.5× sensitivity ratio, L1/L2 are NOT insensitive enough to hide a `ΔC` large
   enough to explain 6:5's own miss, and indeed they don't (they visibly degrade, just not past
   the gate).
4. **Display-rounding cannot rescue the story either.** A correctly-rounded `C = 3.010000`
   (6 decimal places) bounds any true value within `5e-7` absolute (`1.66e-7` relative) of the
   printed number. At 6:5's own measured `C`-sensitivity (`~1569`), that rounding bound explains
   at most `~2.6e-4` of eigenvalue error — almost an order of magnitude short of the observed
   `2.34e-3` miss. (The thesis also prints the same energy elsewhere as `C = 3.01000`, 5 decimal
   places, Fig. 4.5 caption — consistent with a chosen nominal design energy rather than a
   converged quantity carrying hidden precision, which further weakens the "hidden-digits"
   reading.)
5. **The "plausibility" argument citing this module's own C self-validation residual
   (`test_table41_ic_reproduces_stated_jacobi_constant`, `1.05e-6`-`5.20e-6` relative per row) was
   a category error.** That residual measures THIS codebase's own `l*`/`t*` round-trip
   reproducing the PRINTED `C` — it is evidence about this reproduction pipeline's own numerical
   noise floor, not about whether Vaquero's own internal `C` differed from what she printed. It
   also never enters the eigenvalue corrector at all, which holds `C` exactly at the printed
   `3.010000` throughout every row's convergence.

**What survives**: 6:5's own eigenvalue genuinely is far more `C`-sensitive than the other three
rows' (confirmed smooth and monotone over a wide `C` window — not a fold/bifurcation artifact or
bisection noise), and this row's own numerics are independently verified rock-solid on our own
side (Newton residual `3.6e-12`, eigenvalue stable to `1.8e-10` across a 100× `rtol` sweep,
Barden vs. Floquet agree to `7e-10` — ruling out "our own corrector converges more loosely for
6:5" as a mundane explanation). But no single shared parameter — `µ`, `C`, or a joint `µ`+`C`
offset checked against the full sensitivity matrix — reconciles all four rows simultaneously.
The best-supported honest account is that the miss is **row-specific and source-side**: a
transcription or precision issue specific to 6:5's own printed `191.641` (this same table
already carries one other demonstrated erratum, L2's period, per finding 3 below), whose effect
is then amplified through 6:5's own unusually steep `C`-sensitivity — not a coherent
single-parameter correction of any kind.

## Verdict: genuinely investigated honest near-miss, NOT adopted as a fix

This is reported as a **diagnostic finding, not a correction**. Table 4.1 prints `C = 3.010000`
unqualified — unlike `µ`'s own "≈"-qualified display, there is no textual signal in the thesis
that this value is itself a rounded display of a higher-precision internal quantity, so there is
no sourced basis for adopting a different `C` (and, per the adversarial review above, doing so
would not actually reconcile the table even if it were adopted). Per this project's sourced-only
discipline,
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
place to report the exclusion of `µ`. **A first-pass follow-up `C`-sensitivity explanation was
also drafted at that point, but a subsequent Fable adversarial review (2026-08-05, dispatched at
the user's explicit request) found that draft itself overclaimed — see the "adversarially
reviewed and REFUTED" section above — so the docstring was corrected a second time to report the
row-specific/source-side account instead of a single-parameter `C` explanation.** A
regression-guard assertion was also added to
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
close it, and per its own LOW-PRIORITY/exploratory framing that is an acceptable outcome. `µ`
imprecision is empirically EXCLUDED (not just "assumed not to be the whole story") — a bisected
`µ` that fixes 6:5 breaks the other three rows by ~2000×, and the measured per-row sensitivities
rule out a shared `µ`-precision floor as the cause. A follow-up `C`-sensitivity explanation was
also tried and initially looked promising, but a Fable adversarial review (2026-08-05, dispatched
at the user's explicit request after they judged the first-pass diagnostic "inconsistent") found
it does NOT survive scrutiny either: L1 and L2's own baseline eigenvalue errors pin Vaquero's
effective `C` to the printed value ~150× tighter than 6:5 alone would require, and the originally
-claimed "~350-500× more C-sensitive" figure was a factual error (the true ratio against L1/L2 is
only ~3.5×). No single shared parameter — `µ`, `C`, or a joint offset — reconciles all four rows.
What remains well-supported: 6:5's own eigenvalue genuinely is unusually `C`-sensitive (confirmed
smooth, not a numerical artifact) and this row's own numerics are independently verified
rock-solid on our own side, so the miss is best characterized as a row-specific, source-side
issue (transcription or precision imprecision specific to 6:5's own printed `191.641`, amplified
by its own steep `C`-sensitivity) — genuinely investigated and honestly characterized, but not
reducible to a single tidy explanation the way the first draft of this note claimed.
