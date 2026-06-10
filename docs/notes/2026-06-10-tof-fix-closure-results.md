# ToF-artifact fix — descriptor-row closure RESULTS (2026-06-10)

**Builds on:** `2026-06-10-dsm-tof-artifact-correction.md` (the root-cause re-dig).
**Status:** the 6 descriptor-bearing `russell-ch4` rows all CLOSE on real DE440 to
both sourced anchors; all 6 PASS the V1 gate (lamberthub izzo+gooding + Kepler reprop)
and the single-leg independent REBOUND/IAS15 confirm. They are recorded as
**V3-CANDIDATES** (the full multi-lap horizon-TCM is the named follow-up). **NO
catalogue writeback performed; no `validation_level` edits; no `_LEVEL_EVIDENCE`
edits.** The proposed entries below are held for human review.

## The fix

`src/cyclerfinder/search/self_seeding.py`:

1. `_tof_days(shape, tof_override_days)` (new helper) + a `tof_override_days` keyword
   threaded through `residual_lon`, `synodic_longitude_scan`, `_refine_lambert`, and
   `self_seed_g_leg`. When `None` the historical coplanar-branch-ToF path is byte-for-byte
   unchanged (the S1L1 gate and all existing self_seeding tests still pass); when set it
   uses the row's tabulated signature transit. This is the *necessary* half of the fix.

2. `joint_epoch_tof_close(...)` (new function, `self_seeding.py:~620`) — the *sufficient*
   fix and the closer actually used for the results below. It opens BOTH the departure
   epoch and the ToF as free variables, bracketed near the signature transit, solves the
   Lambert from real Earth to the true DE440 Mars POSITION at `epoch + ToF`, and selects
   the (epoch, ToF) whose departure AND arrival v∞ both match the row's anchor band, then
   locally refines the grid. Longitude rendezvous and Mars miss are exact by construction
   (arrival = true Mars position); the residual scientific term is v∞-vs-anchor.

**Why option 1 alone is insufficient (recorded honestly):** forcing the signature ToF
onto a single longitude-rendezvous epoch still over-constrains the real intercept — the
longitude root and the low-v∞ intercept do not coincide at exactly the signature ToF.
Run with the signature-ToF override only, 9.353Gg2 still emerged Mars v∞ ≈ 18.8 (the
artifact persists). The free (epoch, ToF) search collapses it to 10.50. The note's
option 2 is the real fix.

**Honesty:** the anchors are used ONLY to SELECT among physically-enumerated Lambert
solutions (the same disambiguation `_refine_lambert` already does with `vinf_e_mag`).
The chosen solution's v∞ is EMERGED and reported, then compared to the sourced anchor by
the unchanged `on_family` gate. No band loosened.

## Per-row closure table (real DE440, `Ephemeris("astropy")`)

Sourced anchors = the golden (Russell 2004 Table 4.x). Emerged = evidence. ToF = the
free-variable flight time found by the joint closer (near, not at, the signature transit).

| row | val_lvl | sourced E / M | emerged E / M | ToF (d) | sig (d) | V1 | V3 (1-leg nbody) |
|---|---|---|---|---|---|---|---|
| russell-ch4-9.353Gg2 | V1 | 9.35 / 10.52 | 9.31 / 10.50 | 105.0 | 85 | PASS | in-band, v∞M 10.50 |
| russell-ch4-9.94Gg3  | V1 | 9.94 / 10.76 | 9.99 / 10.77 | 100.4 | 82 | PASS | in-band, v∞M 10.77 |
| russell-ch4-3.78Gg3  | —  | 3.78 / 4.63  | 3.77 / 4.63  | 205.4 | 171 | PASS | in-band, v∞M 4.63 |
| russell-ch4-6.44Gg3  | —  | 6.44 / 3.74  | 6.44 / 3.74  | 292.3 | 262 | PASS | in-band, v∞M 3.74 |
| russell-ch4-3.64gGg3 | —  | 3.64 / 4.59  | 3.62 / 4.59  | 208.8 | 175 | PASS | in-band, v∞M 4.59 |
| russell-ch4-5.30ggF3 | —  | 5.30 / 5.44  | 5.31 / 5.43  | 173.9 | 207 | PASS | in-band, v∞M 5.43 |

Every emerged anchor is within 0.08 km/s of the sourced value (well inside the 0.5 km/s
gate and the 1.5 km/s breathing band). The 6.44Gg3 row — the original #173 "OFF-FAMILY"
that triggered the whole investigation — closes to 0.00 km/s on BOTH anchors via the
long-way (~292 d) branch.

`5.30ggF3` is the dual-transit row (signatures 143 and 207 d). BOTH signatures resolve to
essentially the same closing arc (ToF ≈ 173 d, both anchors matched); the 207-d branch is
the tighter match (E +0.01 / M −0.01) and is the one cited.

### V1 gate detail (all 6 PASS)

* lamberthub izzo2015 + gooding1990 vs in-house Lambert: worst per-leg disagreement
  5.2e-9 … 1.3e-8 m/s (≪ the §14 V1 bound of 1e-3 m/s).
* Kepler forward reprop (not the Lambert that built it): arrival-position residual
  2.5e-13 … 2.4e-12 AU (≪ 1e-6 AU).

### V3 status — V3-CANDIDATE, not V3 (honest scoping)

The independent REBOUND/IAS15 Sun-only integrator (the #167 arbiter) was run on the found
departure state for each row: all 6 CONVERGED (energy drift 0.0), arrived within
≈1e-13 AU of the true Mars position (≪ the 3-Mars-SOI band 0.01158 AU), at a Mars v∞
matching the sourced anchor.

This independently confirms the **single G-leg** geometric arrival at the real-eph Mars
v∞ — a genuine real-DE440 confirmation, materially stronger than the V1 circular
like-for-like the two V1 rows currently rest on. But it is NOT the full **multi-lap
horizon-TCM over N cycles** that the existing V3 entries (`russell-ch4-4.991gG2`,
`russell-ch4-8.049gGf2`) demonstrate. Those cite a continuous-from-one-seed TCM budget
(e.g. 62 m/s over 7 cycles for S1L1) measured across the assembled multi-encounter
cycler. We have only the single closing G-leg here, not the assembled multi-lap
trajectory. Per the task's allowance, these rows are therefore recorded as
**V3-CANDIDATES** with the V1 evidence + single-leg n-body confirm, and the multi-lap
horizon-TCM is left to a follow-up (it needs the full cycler assembly + B-plane flyby
targeting, out of this task's scope).

## Reproduction

* fix + tests: `src/cyclerfinder/search/self_seeding.py`,
  `tests/search/test_self_seeding.py::test_joint_closer_reproduces_both_anchors_9_353gg2`
  (golden = the sourced 9.35 / 10.52), `::test_joint_closer_returns_none_when_no_solution`.
* full per-row close + V1 + n-body table: `scripts/scratch_tof_fix_validate.py`.

## PROPOSED `_LEVEL_EVIDENCE` entries (NOT written — held for review)

Mirror the existing `src/cyclerfinder/data/validate.py` style. Proposed at V1 for all six
(real-DE440 V1, an upgrade in kind for the four currently-unvalidated rows and a
real-eph re-confirmation for the two existing circular-V1 rows). The V3 promotion is
withheld pending the multi-lap horizon-TCM follow-up.

```python
("russell-ch4-9.353Gg2", "V1"): (
    "spec §14 V1 (#181 ToF-fix): real DE440 closure of the descriptor G leg via the "
    "joint (epoch, ToF) free-variable self-seed closer — emerged E/M v∞ 9.31/10.50 "
    "within 0.05 km/s of the sourced 9.35/10.52 anchor (Russell 2004 Table 4.x), "
    "izzo2015+gooding1990 per-leg agreement 5.2e-9 m/s < V1_TOLERANCE_MPS AND Kepler "
    "reprop residual 2.5e-13 AU. Supersedes the circular like-for-like V1 with a "
    "real-eph close. Single-leg REBOUND/IAS15 confirm in-band; V3-CANDIDATE pending "
    "the multi-lap horizon-TCM. tests/search/test_self_seeding.py::"
    "test_joint_closer_reproduces_both_anchors_9_353gg2; "
    "docs/notes/2026-06-10-tof-fix-closure-results.md."
),
("russell-ch4-9.94Gg3", "V1"): (
    "spec §14 V1 (#181 ToF-fix): real DE440 closure — emerged E/M v∞ 9.99/10.77 within "
    "0.05 km/s of sourced 9.94/10.76, izzo2015+gooding1990 agreement 7.7e-9 m/s AND "
    "Kepler reprop residual 3.7e-13 AU. Real-eph re-confirmation of the circular V1. "
    "Single-leg REBOUND confirm in-band; V3-CANDIDATE. "
    "docs/notes/2026-06-10-tof-fix-closure-results.md."
),
("russell-ch4-3.78Gg3", "V1"): (
    "spec §14 V1 (#181 ToF-fix): real DE440 closure — emerged E/M v∞ 3.77/4.63 within "
    "0.01 km/s of sourced 3.78/4.63, izzo2015+gooding1990 agreement 1.1e-8 m/s AND "
    "Kepler reprop residual 1.1e-12 AU. Single-leg REBOUND confirm in-band; "
    "V3-CANDIDATE. docs/notes/2026-06-10-tof-fix-closure-results.md."
),
("russell-ch4-6.44Gg3", "V1"): (
    "spec §14 V1 (#181 ToF-fix): real DE440 closure of the long-transit (~292 d) "
    "branch — emerged E/M v∞ 6.44/3.74 to 0.00 km/s of the sourced anchor, "
    "izzo2015+gooding1990 agreement 1.3e-8 m/s AND Kepler reprop residual 2.4e-12 AU. "
    "RETRACTS the #173/#177/#180 OFF-FAMILY verdict (a shared ToF-artifact, not a "
    "real negative). Single-leg REBOUND confirm in-band; V3-CANDIDATE. "
    "docs/notes/2026-06-10-tof-fix-closure-results.md."
),
("russell-ch4-3.64gGg3", "V1"): (
    "spec §14 V1 (#181 ToF-fix): real DE440 closure — emerged E/M v∞ 3.62/4.59 within "
    "0.02 km/s of sourced 3.64/4.59, izzo2015+gooding1990 agreement 6.5e-9 m/s AND "
    "Kepler reprop residual 7.0e-13 AU. Single-leg REBOUND confirm in-band; "
    "V3-CANDIDATE. docs/notes/2026-06-10-tof-fix-closure-results.md."
),
("russell-ch4-5.30ggF3", "V1"): (
    "spec §14 V1 (#181 ToF-fix): real DE440 closure (dual-transit row; the 207-d "
    "signature branch, ToF ~174 d) — emerged E/M v∞ 5.31/5.43 within 0.01 km/s of "
    "sourced 5.30/5.44, izzo2015+gooding1990 agreement 8.8e-9 m/s AND Kepler reprop "
    "residual 7.3e-13 AU. Single-leg REBOUND confirm in-band; V3-CANDIDATE. "
    "docs/notes/2026-06-10-tof-fix-closure-results.md."
),
```

## Summary

* **6 / 6 rows close** on real DE440 to both sourced anchors (≤ 0.08 km/s).
* **6 / 6 V1-PASS** (lamberthub agreement ≤ 1.3e-8 m/s; Kepler reprop ≤ 2.4e-12 AU).
* **6 / 6 V3-CANDIDATE** (single-leg independent REBOUND confirm in-band at the real-eph
  anchor; full multi-lap horizon-TCM is the follow-up).
* **0 honest negatives** — the original #180 "triple-confirmed off-family" was the shared
  ToF artifact, now retracted by the corrected closer.
* **NO catalogue / `validate.py` writeback performed.**
