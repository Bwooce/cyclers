# Continuation driver — results & golden-gate break-point analysis (task #158)

Companion to the spec
`docs/notes/2026-06-07-russell-2004-continuation-deepdive.md`. Records the
outcomes of the circular→ephemeris continuation driver
(`src/cyclerfinder/search/continuation.py`) against the two sourced golden gates
from **Russell 2004 dissertation, Table 5.5 (p.178)**.

## What was built

`continuation_correct` implements the Russell §5.4 model-fidelity homotopy
(deep-dive §1): a free-return closure on the circular-coplanar planet model is
walked out to the true ephemeris by ramping the **solar-system model**, not the
genome. Schedule per rung: a leading **phase** ramp (per-body J2000 mean-longitude
alignment, our addition — see below), then **e**, then **i**, then one final step
to DE440 (`Ephemeris("astropy")`). The large elements `a, Ω, ω, ν` are frozen at
their Table-5.4 J2000 means throughout. The `nstep` ladder `{1, 3, 9, 27, 81}`
(243 skipped for wall-time, recorded in `ContinuationResult.skipped` — no silent
cap) is run rung-by-rung; the lowest-residual final is kept. A full per-step audit
trail (`ContinuationStep`/`ContinuationRung`/`ContinuationResult`) records every
step's λ-trio, model name, residual, `(a, e, t0)` and converged flag.

The inner solve is the existing `search/free_return.free_return_correct`, imported
unchanged. Intermediate models are a thin in-module `_RampedElementsBackend`
injected post-construction into a real `Ephemeris` — **no `core/` file edited**.
At `λ_e = λ_i = λ_p = 0` the ramped backend is bit-identical to
`Ephemeris("circular")` (a fast mechanics gate).

### Why a phase ramp was added to Russell's e→i→ephemeris schedule

Our circular seed lives in a `θ=0-at-J2000` frame; the true ephemeris does not.
The ~100° per-body absolute-longitude offset cannot be absorbed by the corrector's
single `t0` (different bodies need different shifts). Russell folds this into his
Fig.5.4 seeding step ("propagate the simple model until the parent's beginning
phase angle is achieved"); we make it an explicit leading homotopy leg
(`λ_p: 0→1`) so each step stays a small perturbation. This is an
orchestration-level addition, not a deviation from the physics.

### Non-singular-elements note (deep-dive §3 / failure mode 4)

Russell switches his analytic Jacobian to a non-singular `β` element set because
classic-element partials blow up at `e→0, i→0`. Our inner genome is `(a, e, t0)`
with the spacecraft `e` bounded away from zero (`0 < e < 0.95`) and a residual
built from heliocentric longitudes + a Mars-reach margin — it never forms a
classic `(ω, ν)` partial at the circular endpoint. The planet-model singularity
Russell guards against is supplied here as Cartesian state, not differentiated. So
the genome is already non-singular for this homotopy; documented, not re-derived.

## Mechanics gates (fast, no astropy) — ALL GREEN (8)

- `λ=0` reproduces `Ephemeris("circular")` bit-identically (E and M, several t).
- A single e-ramp step moves Mars's state continuously (no jump; in-ecliptic).
- The i-ramp tilts the plane while preserving `|r|` (pure rotation).
- The phase ramp aligns Mars's J2000 longitude with the Table-5.4 mean longitude.
- Schedule shape: `p-ramp(nstep) + e-ramp(nstep) + i-ramp(nstep)`.
- Ladder bookkeeping: requested ladder run, skipped rungs recorded, keep-best.
- Audit trail completeness.

## Golden gates (slow, DE440 via astropy) — Russell Table 5.5, p.178

EXPECTED side is the sourced "total maintenance Δv = 0 m/s / 7 cycles" claim; our
evidence is the continued solution's final (true-ephemeris) step residual being
below the documented closure floor `TOL_KMS = 0.1 km/s` (Russell's "0 m/s" is
below his unprinted SNOPT/post-processing floor, deep-dive §7).

### (a) Aldrin 6.399G1, launch Aug 6 2003 — **PASS** (hard gate)

| quantity | value |
|---|---|
| winning ladder rung | `nstep = 3` |
| final (DE440) step residual | **0.00158 km/s** (≈0 within `TOL_KMS = 0.1`) |
| converged `(a, e)` | 1.5249 AU, 0.3616 |
| launched `t0` | 2003-07-11 (**26 d** from sourced 2003-08-06; < 200 d gate) |
| emerged V_inf | E = 6.08, M = 8.88 km/s (anchors 6.5 / 9.7; within 1.5) |
| 243 skipped & recorded | yes |

The Aldrin cycler closes ballistically end-to-end at the true ephemeris in the
sourced launch window — the headline reproduction of Russell's named, dated,
0-m/s Aug-6-2003 Aldrin result, and a full end-to-end validation of the
continuation driver.

### (b) 4.991gG2 (S1L1), launch Jun 2025 — **strict-xfail break point** (NEVER softened)

This gate fails — but **not** on the driver and **not** on the residual. It fails
on launch-window proximity, and the root cause is the inner single-ellipse genome,
not the continuation infrastructure.

What the driver produces from the (Dec-2026) seed:

| quantity | value |
|---|---|
| winning ladder rung | `nstep = 81` |
| final (DE440) step residual | 0.00153 km/s (ballistic — the corrector DID close) |
| converged `(a, e)` | 1.3147 AU, 0.2637 (≈ sourced 1.30 / 0.257) |
| emerged V_inf | E = 4.93, M = 5.34 km/s (anchors 4.99 / 5.10; within 1.5) |
| launched `t0` | 2026-12-10 (**543 d** from sourced 2025-06-15 — FAILS < 200 d) |

So the driver finds *a* ballistic single-ellipse free-return with the right
elements and V_inf — but in a **different phase basin** ~543 d from Russell's
sourced June-2025 window.

**Seed-sweep diagnostic (the decisive break-point evidence).** Seeding the
corrector AT or near the sourced Jun-2025 window does not recover Russell's
solution — it diverges to the genome's `e→0.95` bound:

| seed | nstep | converged | residual (km/s) | t0 | off (d) | a / e | V_inf E/M |
|---|---|---|---|---|---|---|---|
| 2025-06-15 | 3  | **no**  | 34.57 | 2025-06-15 | 0   | 0.881 / 0.949 | 33.6 / 21.5 |
| 2025-05-01 | 81 | **no**  | 23.22 | 2025-05-01 | 45  | 0.886 / 0.950 | 33.7 / 21.6 |
| 2025-08-01 | 81 | **no**  | 46.10 | 2025-08-01 | 47  | 0.890 / 0.950 | 33.8 / 21.7 |
| 2025-03-01 | 3  | **no**  | 6.13  | 2025-03-01 | 106 | 0.888 / 0.950 | 33.8 / 21.7 |
| 2024-12-01 | 9  | yes     | 0.011 | 2024-12-01 | 196 | 1.370 / 0.351 | 8.4 / 8.1   |

At the sourced window the single-ellipse free-return arc collapses onto the
degenerate high-eccentricity edge (`e→0.95`, residual 23–46 km/s): **the genome
cannot represent Russell's S1L1 closure at June 2025.** The only seeds that close
ballistically land in other basins, and the closest one (2024-12-01) does so with
off-family elements (a=1.37, e=0.35, V_inf 8.4/8.1 — not the sourced 1.30/0.257,
4.99/5.10).

**Root cause (recorded S1L1 multi-arc blocker, MEMORY.md):** S1L1's real closure
is **two generic-return arcs**, not one symmetric free-return ellipse. The
single-ellipse `(a, e, t0)` genome the driver wraps cannot model that, so the
sourced Jun-2025 window is not a residual-zero point for *this genome* even though
the continuation machinery around it is correct. This is consistent with the
Aldrin result: Aldrin is a clean one-synodic symmetric free-return (genome-native)
and closes perfectly; S1L1 is two-synodic multi-arc (genome cannot express it).

**Disposition:** strict-`xfail` (the parametrized `_S1L1_XFAIL` case). It is NOT
softened — the 200-day window assertion stays at 200 days. The day this test
starts passing, the inner genome will have gained multi-arc S1L1 representability
and the xfail marker must be removed (strict mode forces that by failing on an
unexpected pass).

## Conclusion

The continuation driver is complete and validated: Aldrin reproduces Russell's
0-m/s Aug-6-2003 result end-to-end (driver gate green), and the S1L1 failure is a
sourced, documented genome-representability break point (driver-correct,
genome-limited), recorded here and pinned by a strict xfail. The unblock for S1L1
is a multi-arc free-return genome in `search/free_return.py` (out of scope for
this task — the driver will consume it unchanged once it exists).
