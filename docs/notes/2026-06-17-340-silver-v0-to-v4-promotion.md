# #340 — SILVER (Umbriel-Oberon-Umbriel) `validation_level` V0 → V4 promotion

**Date:** 2026-06-17
**Task:** #340 — register the SILVER's `_LEVEL_EVIDENCE` entries (pytest-evidenced gates) and promote the catalogue row's `validation_level` from V0 to V4
**Follow-up to:** #339 (commit `34566d5`) — admitted the SILVER as `umbriel-oberon-1-1-uranian-quasi-cycler-2026` at `validation_level: V0` because the 10-gate provenance was in project-output JSONLs, not pytest-registered evidence
**Catalogue:** 283 rows (unchanged), one row promoted V0 → V4

## What #339 deferred

#339 (`2026-06-17-339-silver-quasi-cycler-admission.md`) admitted the SILVER as
the catalogue's first computed `quasi_cycler` row with full provenance across
12 tasks. `src/cyclerfinder/data/validate.py::_LEVEL_EVIDENCE` (the over-claim
guard registry consulted by `validate_validation_level` per spec §16.7.12)
requires a sourced pytest evidence pointer for any row claiming `V1+`. The
SILVER's gauntlet evidence lived only in JSONL form, not pytest form:

```text
data/silver_327_v1_v2_verdicts.jsonl       # #306 Phase 1 Part D (V1+V2 moontour)
data/silver_327_moontour_v2_verdicts.jsonl # #330 V2 quasi-cycler scan
data/silver_327_v3_verdicts.jsonl          # #331 REBOUND IAS15 V3
data/silver_327_v4_strict_verdicts.jsonl   # #335 URA111 real-eph V4
data/silver_327_v4_strict_annual_sweep_338.jsonl  # #338 Part A
data/silver_327_v4_strict_boundary_338.jsonl       # #338 Part B (boundary verdict)
```

So #339 set `validation_level: V0` (the internal-consistency floor) and tracked
the V4 promotion as a follow-up — which is this task.

## What #340 did

### Part A — Four frozen-gate pytest modules

```text
tests/verify/test_silver_327_v1_passes.py
tests/verify/test_silver_327_v2_quasi_cycler.py
tests/verify/test_silver_327_v3_nbody.py
tests/verify/test_silver_327_v4_strict_passes.py
```

These are **wrappers, not recomputations**. Each one loads the relevant JSONL
verdict file, locates the SILVER's row by `candidate_id =
repeated-moon-uranus-00000041`, and asserts the catalogue's claim matches the
project-output evidence. If a JSONL is regenerated with a different answer the
gate breaks by design — that's exactly what the over-claim guard is supposed to
catch.

The SILVER's quasi_cycler structural signature is the bounded oscillation:

| Tier | What the gate asserts |
|---|---|
| **V1** | `passes_v1=True`, independent DOP853 cross-check residual `2.30e-10 km/s ≪ 1e-3 km/s` (spec §14 V1 floor) |
| **V2** | 10/10 Lambert cycles convergent at all legs; drift oscillates between ~86,000 km (cycle 5) and ~530,000 km (cycle 7); strict V2 verdict `FAIL_QUASI_BOUNDED` (BY DESIGN — passing strict V2 would imply strict-periodic, not quasi) |
| **V3** | REBOUND IAS15 vs V2 driver `drift_agreement < 100 km` across `n_cycles ∈ {3, 5, 10}` (actual: max `1.84e-6 km` at `n=10`, 8 orders below floor); `passes_v3=True` |
| **V4** | #338 boundary `verdict_label = EFFECTIVELY_CYCLIC`; interior PASS rate `85/85 = 100%`; longest PASS run `84 yr` (2000-06-21 → 2083-06-21, centred 2041) |

All 10 gate-asserting tests pass (`uv run pytest tests/verify/test_silver_327_v*.py
→ 10 passed in 0.93s`). Commit: `0813cc4`.

### Part B — `_LEVEL_EVIDENCE` entries

Added four entries to `src/cyclerfinder/data/validate.py::_LEVEL_EVIDENCE`,
following the established header-comment + per-entry format used by the
Aldrin / Russell / Ross / Liang / Braik / Brian families:

```python
("umbriel-oberon-1-1-uranian-quasi-cycler-2026", "V1"): "... #306 ... tests/verify/test_silver_327_v1_passes.py.",
("umbriel-oberon-1-1-uranian-quasi-cycler-2026", "V2"): "... #330 ... tests/verify/test_silver_327_v2_quasi_cycler.py.",
("umbriel-oberon-1-1-uranian-quasi-cycler-2026", "V3"): "... #331 ... tests/verify/test_silver_327_v3_nbody.py.",
("umbriel-oberon-1-1-uranian-quasi-cycler-2026", "V4"): "... #335 + #338 ... tests/verify/test_silver_327_v4_strict_passes.py.",
```

Commit: `158c21b`.

### Part C — Catalogue row promotion

```diff
- validation_level: V0   # spec §14: V0 internal-consistency floor. The 10-gate provenance ... lives in data/silver_327_v*_verdicts.jsonl — project-output JSONLs, NOT pytest-evidenced gates registered in src/cyclerfinder/data/validate.py::_LEVEL_EVIDENCE. Promotion to V4 requires registering pytest evidence per spec §16.7.12; tracked as follow-up task.
+ validation_level: V4   # spec §14: V4-strict URA111 real-eph + #338 EFFECTIVELY_CYCLIC. The 10-gate provenance ... lives in data/silver_327_v*_verdicts.jsonl AND #340 (commit pending) wraps each V-tier's JSONL verdict as a frozen-gate pytest ... and registers V1/V2/V3/V4 in src/cyclerfinder/data/validate.py::_LEVEL_EVIDENCE per spec §16.7.12.
```

## Discipline notes (per `feedback_golden_tests_sourced_only`)

The four frozen-gate tests do NOT compute their own EXPECTED values. The
EXPECTED side is the JSONL verdict file on disk. Those JSONLs are sourced by
the project's own #306 / #330 / #331 / #332 / #335 / #338 chain — they are
**project-output evidence**, not paper-sourced anchors. The catalogue claim is
then tied to that evidence: if the row says "V4 EFFECTIVELY_CYCLIC over 2000-2083"
the test asserts the JSONL says the same thing.

This is the same wrapper pattern used elsewhere for project-internal evidence
(e.g. `tests/verify/test_aldrin_v2_powered.py` asserts the
`verify_aldrin_v2_powered(...)` driver returns the V2-powered verdict claimed
by the catalogue row's `_LEVEL_EVIDENCE` entry; the driver re-runs a real-DE440
BVP, but the tolerances and structural verdicts are catalogue-row claims). The
SILVER's wrapper is lighter (just JSONL re-read, no recompute) because the
verdict JSONLs were already frozen by #306 / #330 / #331 / #335 / #338.

Per `feedback_published_rounded_values_are_display`: the catalogue row's V_inf
values are display-rounded (4 sig figs); the JSONLs carry the full precision.
The gates assert on the full-precision JSONL values, never the display rounds.

## Why V4 and not V5

V5 would require an INDEPENDENT-TOOLCHAIN cross-check (e.g. GMAT R2022a's own
URA111 + scipy DOP853 pipeline reproducing the #338 EFFECTIVELY_CYCLIC verdict
end-to-end). The project's V5 lane is the GMAT lane (per
`reference_gmat_install`); that's a separately-tracked follow-up, not within
#340's scope. Until then the SILVER stays at V4 — the highest tier with
mechanical evidence on this project's own toolchain.

## Census ratchet

The validation-tier census in `tests/data/test_validation_tier_census.py` is
frozen by spec — this promotion bumps the V4 count by 1 (V0 → V4). If the
ratchet test breaks, the frozen counts there need updating to reflect the new
distribution. Verified: ratchet test passes after the promotion (no changes
needed to the frozen counts — the SILVER was a recent addition and the ratchet
was last updated to accommodate it).

## Cross-task awareness

- #320 (Saturn moontour expansion) still in flight (`scripts/scan_320_*.py`,
  `data/scan_320_*.jsonl`). #340 stayed out of those paths.
- #310 (single-orbit prioritizer continuation) still in flight
  (`src/cyclerfinder/search/single_orbit_prioritizer*.py`). #340 stayed out.
- #340 only touched: `tests/verify/test_silver_327_v*.py` (NEW),
  `src/cyclerfinder/data/validate.py` (`_LEVEL_EVIDENCE` addition only),
  `data/catalogue.yaml` (one-line `validation_level` flip + comment refresh),
  `docs/notes/2026-06-17-340-silver-v0-to-v4-promotion.md` (this file).

## Follow-up

- V5 via GMAT independent-toolchain — separately-tracked, not pre-conditioned
  on #340 closing.
- Sub-year DOY boundary characterization (the annual sweep is aliased to the
  Umbriel-Oberon synodic period 5.987 d) — flagged in #338 Part B notes
  (`recommended_launch_window` doc + `aliasing_note` in the boundary JSONL).
- Post-URA111 Uranian satellite kernel acquisition to extend the
  `validity_window` beyond 2083 — currently capped by the URA111 1900-2099
  coverage window.
