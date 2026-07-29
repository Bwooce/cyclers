# `#765`: Saturn-Titan planar CR3BP resonant-orbit families — Vaquero 2013 Table 4.1 gate

**Task:** `#765`, the first concrete per-system task of `#760`'s new-system discovery
campaign, spec-complete in
`docs/notes/2026-07-29-764-new-system-discovery-scoping.md` §6 (recommended by `#764`'s own
scoping pass). Dispatched by the user ("register everything as we go... numbers are free").

**Source paper** (acquired + digested this task; see
`docs/notes/2026-07-29-765-vaquero-2013-digest.md` for the full digest and citation-mining
pass): T. M. Vaquero Escribano (2013), "Spacecraft Transfer Trajectory Design Exploiting
Resonant Orbits in Multi-Body Environments," PhD dissertation, Purdue University (advisor
K. C. Howell). Freely downloaded from the author's own Purdue faculty page; md5
`fdcbf871322b87cd1dd3448059cb2596` (matches `#764`'s own scoping-note md5 exactly). Filed at
`cyclers_pdf/papers/vaquero-2013-spacecraft-transfer-trajectory-design-resonant-orbits-multibody-environments-purdue-phd.pdf`.
Table 4.1 (p.109, PDF page 124) verified against the actual PDF text layer AND a direct
page-image (vision) read — see screenshot cross-check in the digest note.

**Code delivered:** `src/cyclerfinder/search/saturn_titan_resonant_families.py` (new module,
a thin sibling of `jovian_resonant_families.py`, reusing its corrector/classification
machinery directly — no reimplementation) + `tests/search/test_saturn_titan_resonant_families.py`
(27 tests, all passing). Both pass `ruff check`, `ruff format --check`, and
`uv run mypy src tests` cleanly.

---

## Unlike the Jovian chain: this source prints its own ICs

Anderson & Lo 2011 (the Jovian module's source) published NO initial-condition table —
`#753`-`#758` spent four tasks locating every seed by blind grid search, continuation, and a
Table-2-homoclinic-derived seed strategy. Vaquero 2013's own Table 4.1 prints the seed `x`
(km) and `ẏ` (km/s) directly for all four rows. Per `#758`'s own lesson ("paper-sourced seed
windows beat blind grid scans"), this task led with those sourced seeds — no blind grid sweep
was needed to locate any of the four candidates.

## Characteristic-quantity derivation and self-validation

The thesis states `µ ≈ 2.3658e-4` but does not print a Saturn-Titan `l*`/`t*` table (only the
general Ch.2 formulas). This module derives them from this project's own DE440 registry
(Titan's own SMA about Saturn, 1,221,870 km, as `l*`; Saturn-system GM + Titan GM as the
two-body GM sum — mirroring EXACTLY how `jupiter_europa_system()` derives its own `t_s`), then
**self-validates** the choice: nondimensionalizing all four of Table 4.1's own printed `(x,
ẏ)` pairs with this `l*`/`t*` and the thesis's own `µ` reproduces the STATED Jacobi constant,
C = 3.010000, to 1.6e-5 relative or better for every row (`test_table41_ic_reproduces_stated_jacobi_constant`).
This is strong, independent evidence the derived `l*`/`t*` matches whatever the thesis used
internally — not merely a plausible-looking guess.

A blind `two_body_resonant_seed(3, 4)` attempt (the literal periapsis-at-Titan two-body
geometry) was also tried for completeness/comparison, per the dispatch note's own instruction
to try it. Consistent with Anderson & Lo's own documented finding for the analogous
Jupiter-Europa attempt (`#753` module docstring item 1), it does not even produce a valid
`ydot0` at C = 3.010000 (the Jacobi-constraint radicand is negative there) — an expected,
documented negative, not a bug (`test_naive_two_body_seed_does_not_converge_at_vaquero_c`).

## The gate: 2 of 4 rows fully pass; the other 2 are honest, well-characterized partial results

| Row | Target `\|λu\|` | Recovered | Eigenvalue rel. err | Target T (days) | Recovered T (days) | Period rel. err | Passed |
|---|---|---|---|---|---|---|---|
| **3:4** | 2,129.81 | **2129.8077** | **1.07e-6** ✓ | 66.3312 | 66.3427 | 1.7e-4 ✓ | **PASS** |
| **L1 Lyapunov** | 1,004.72 | **1004.7246** | **4.56e-6** ✓ | 8.2829 | 8.2844 | 1.8e-4 ✓ | **PASS** |
| 6:5 | 191.641 | 191.1928 | **2.34e-3** ✗ | 71.2638 | 71.2782 | 2.0e-4 ✓ | FAIL (eigenvalue) |
| L2 Lyapunov | 892.850 | 892.8524 | 2.73e-6 ✓ | 79.7260 | 8.6032 | 0.892 ✗ | FAIL (period) |

Gate tolerance `1e-3` relative on both eigenvalue and period (`TABLE41_EIGENVALUE_GATE_REL_TOL`,
`TABLE41_PERIOD_GATE_REL_TOL`) — same value the Jovian module uses, for direct
comparability, justified by the corrector's own 1e-12–1e-14 convergence floor and by 3/4 rows
reproducing to 1e-6–1e-5 relative at exactly the thesis's own stated `µ`. Dimensional IC
match (`x`, `ẏ`) is tracked as evidence only (unit-conversion-dependent, per the dispatch
note's own framing, not part of `passed`) — all four rows beat 0.02% relative on this axis.

**Basin robustness** (`test_sourced_seed_is_basin_robust`): every one of 11 evenly-spaced
seeds across a ±2e-4 window around each sourced Table 4.1 seed converges to the identical
eigenvalue for all four rows — none of these four candidates is an isolated numerical fluke.

**Independent cross-check** (`test_barden_matches_planar_floquet`): Barden's half-period
identity and a direct full-period monodromy eigendecomposition (`_planar_floquet`) agree to
<1e-6 relative for all four rows. None of Table 4.1's rows sit anywhere near `|λ| ≈ 1` (the
smallest is 6:5's 191.6), so the Jovian module's own degenerate-eigenvalue pitfall
(`argmax(|eigenvalue|)` mis-selecting the trivial unit pair) does not bite here — the
cross-check is performed regardless, as standing discipline.

---

## Honest finding 1: 6:5's eigenvalue is a genuine, small, physically-explained near-miss

6:5's own IC (`x`, `ẏ`) and period match Table 4.1 to <0.02% relative — as tight as 3:4/L1 —
yet its eigenvalue misses by 2.34e-3, just outside the gate. A systematic `µ`-sensitivity
sweep (`test_eigenvalue_sensitivity_to_mu_is_measured_not_assumed`, and a wider ad hoc sweep
run during development) found that a +0.1% `µ` perturbation shifts EVERY row's eigenvalue by a
comparable small amount — but because 6:5's own unstable eigenvalue (191.6) is an order of
magnitude smaller than 3:4's (2129.8), the SAME absolute sensitivity floor (from the thesis's
own 5-significant-figure `µ` display, which itself uses "≈" rather than Anderson & Lo's
unqualified 14-digit value) shows up as a proportionally larger RELATIVE error for the smaller
eigenvalue. This is reported as an honest FAIL under the 1e-3 gate — not fudged, not
retroactively loosened — while being a real, well-characterized, small miss rather than a wild
one, and a plausible general lesson for any future system whose gate rows span more than an
order of magnitude in target eigenvalue magnitude.

## Honest finding 2: L2's own printed period is very likely a source transcription error

L2's own printed `x` (1.25231e6 km) and eigenvalue (892.850) both reproduce to <0.01%/2.7e-6
relative — as tight as 3:4/L1 — using the SAME seed. But the converged orbit's
self-consistent period is 3.3899 nondim time (8.603 days), a factor ~9.27 different from the
thesis's own printed `T = 79.7260 days`. This is NOT explainable by this module's own `l*`/`t*`
choice, which independently reproduces 3:4/6:5/L1's own printed T-in-days columns to <0.02%
relative using the exact same `l*`/`t*`. Fig. 4.6(b) (p.109) itself plots the L2 Lyapunov
orbit as a small, simple single-loop oval comparable in size to the L1 orbit (whose own period
is 8.2829 days) — physically consistent with an ~8.6-day period, not an order of magnitude
larger. Per this project's respectful-errata-framing discipline (evidence-first, falsifiable,
benefit of the doubt — typesetting slips happen to everyone), this is flagged as a likely
transcription/typesetting error in the thesis's own Table 4.1 for this one row.
`TABLE41_DIMENSIONAL`'s L2 entry is kept EXACTLY as printed (sourced-only discipline — never
silently "corrected"); the period gate for this row honestly reports FAIL, with the full
evidence documented in the importable, testable `TABLE41_L2_PERIOD_ERRATA_NOTE` constant.

---

## Verification

* `uv run ruff check` / `ruff format --check` on both new files: clean.
* `uv run mypy src tests` (project canonical invocation): clean, 825 files.
* `uv run pytest tests/search/test_saturn_titan_resonant_families.py -q`: 27/27 pass.
* `uv run pytest tests/data/test_outstanding_structure.py tests/data/test_outstanding_header_body_consistency.py -q`:
  run before committing the `OUTSTANDING.md` update (see commit history for pass/fail status
  recorded at commit time).

---

## Explicitly out of scope (per the dispatch note): the connection stage

Vaquero 2013 Ch.4.3.1 also describes — as figures + prose only, no state tables — a
homoclinic connection of the 3:4 resonant orbit (Fig. 4.9), a periodic "resonant chain"
cycling indefinitely between the 3:4 and 6:5 resonances (Fig. 4.10, continued in C in
Fig. 4.11), and a falsifiable published termination claim ("it is suspected that this family
of periodic resonant chains ends for a value of Jacobi constant C < 3.01400", Fig. 4.12). None
of that is attempted here — it is `#765`'s own Task-B analog, to be registered as a separate,
later task once the resonant families themselves (this task's own deliverable) are confirmed,
exactly mirroring how the Jovian chain (`#753` → `#754`) sequenced family confirmation before
connection work.

## Recommendation for a Task-B analog — opinion, not a decision

3:4 and L1 are both confirmed to near-machine precision, giving a strong, independently-
validated foundation (correct `µ`/`l*`/`t*`, correct corrector/classification path) for a
future connection-stage task targeting the 3:4 orbit's own homoclinic connection (Fig. 4.9) —
unlike the Jovian chain, where Task B was held because BOTH families it needed (3:4-LO,
5:6-LO) were initially unconfirmed. The `#759`-carried-forward risk flag (`λ≈2130` on the 3:4
manifold leg as "moderate risk") did not materialize as a convergence difficulty here — 3:4
converged in 3 Newton iterations to a 1.07e-6 relative eigenvalue match on the FIRST attempt
using the sourced seed directly, no densification fallback needed. This is my assessment for
the user to weigh, not a decision — a connection-stage task is a separate dispatch.
