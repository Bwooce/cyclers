# Continuation batch — V1→V3 lift results (2026-06-08)

Generalised the #158 single-row continuation driver into a batch tool
(`src/cyclerfinder/search/continuation_batch.py`) with a per-row V3 evidence gate
(`tests/verify/test_continuation_v3_batch.py`, slow / DE440), and ran it over the
**V1→V3 lift set**: the four Russell free-return rows currently at
`validation_level: V1` plus the Aldrin outbound (already V2).

This is the only set the continuation can reach — it walks UP from a
circular-coplanar CLOSURE, so it can only lift rows that already close
circular-coplanar (the #137 set + Aldrin). Rows lacking that closure (the
multi-arc / Lambert-singular set, and S1L1) cannot be seeded here.

## Method (per row)

1. **Seed** (CONSTRAINT): `(a, e)` from the sourced aphelion + outbound transit
   (`free_return.seed_ae_from_aphelion_transit`, the #106 shared physics) for the
   Russell rows; sourced `(a, e)` directly for Aldrin.
2. **Basin search**: the circular phase frame does not coincide with DE440's, so a
   circular-best basin sits ~600 d off the real window. Instead sweep seed dates
   across ±200 d of the sourced launch (10 d steps), run the #158 ladder
   (`ladder=(3,)`, the documented winning rung) at each, and **select** the
   converged basin whose emerged per-body V∞ best matches the INDEPENDENTLY
   sourced anchor (V∞ is the row's physical fingerprint — distinct DE440 basins
   are distinct cyclers), tie-broken toward the window.
3. **V3 half (a) — phase-match**: ballistic close at DE440 (residual < `TOL_KMS`
   0.1 km/s) AND emerged V∞ within 1.5 km/s of the sourced anchor (the #158
   golden band, coplanar-vs-3D) AND `t0` within the §14 window gate (200 d).
4. **V3 half (b) — bounded horizon TCM**: from the continued `t0`, re-phase one
   E–M synodic (2.135 yr) per lap over 3 laps and re-solve the in-family
   maintenance ΔV on DE440 (Aldrin via the #161 `optimise_aldrin_maintenance_dv`
   resolver wrapper; Russell rows via the body-agnostic `optimise_maintenance_dv`
   seeded from the per-lap V∞-matched basin, `closure_body="E"`). Gate: every lap
   converges in-family and stays under the engineering plausibility bar
   (`MAINTENANCE_DV_CONVENTION_KMS` = 3.0 km/s — a CONVENTION, not a sourced
   per-row budget; none is published).

Golden discipline: sourced anchors are EXPECTED; the emerged V∞, phase-matched
`t0` and per-lap ΔV are EVIDENCE. The bar is the project's own convention. No
catalogue writeback is performed by these files (the main session consolidates).

## Per-row outcomes (verbatim, DE440)

### Half (a) — phase-match (continuation closure)

| row | sel basin | resid km/s | emerged V∞ E/M | sourced V∞ E/M | V∞ off E/M | t0 | win off d | half (a) |
|---|---|---|---|---|---|---|---|---|
| russell-ch4-5.30gGf3 | +150 d | 0.00486 | 5.70 / 8.42 | 5.30 / 9.17 | 0.40 / 0.75 | 2003-07-09 | 150 | **PASS** |
| russell-ch4-9.94Gg3  | +180 d | 0.00002 | 10.38 / 10.38 | 9.94 / 10.76 | 0.44 / 0.38 | 2003-08-08 | 180 | **PASS** |
| russell-ch4-5.75ggF3 | +150 d | 0.00513 | 5.76 / 8.56 | 5.75 / 9.36 | 0.01 / 0.81 | 2003-07-09 | 150 | **PASS** |
| russell-ch4-9.353Gg2 | (none) | 20.617 | 4.60 / 8.09 | 9.35 / 10.52 | 4.75 / 2.43 | 2003-02-21 | 200 | **BREAK** |
| aldrin…-outbound     | −20 d  | 0.00098 | 6.72 / 9.09 | 6.50 / 9.70 | 0.22 / 0.61 | 2003-07-17 | 20  | **PASS** |

`9.353Gg2` (the deep-aphelion high-e row) finds NO converged, V∞-matched,
in-window basin in the ±200 d sweep — its DE440 continuation lands off-fingerprint
(V∞ E 4.6 vs sourced 9.35). It breaks at half (a).

### Half (b) — bounded horizon TCM (3 laps, per-cycle ΔV km/s)

| row | lap 0 | lap 1 | lap 2 | all < 3.0 | half (b) |
|---|---|---|---|---|---|
| russell-ch4-5.30gGf3 | 0.000 | 0.000 | **12.980** | no | **BREAK** |
| russell-ch4-9.94Gg3  | 0.000 | 0.000 | **13.182** | no | **BREAK** |
| russell-ch4-5.75ggF3 | 0.000 | 0.000 | 1.214 | yes | **PASS** |
| aldrin…-outbound     | 0.000 | 0.000 | 0.000 | yes | **PASS** |

`5.30gGf3` and `9.94Gg3` close & phase-match (half a) but their 3-lap horizon TCM
exceeds the bar: lap 2 (≈2007-10) re-solves onto a high-e member (e ≈ 0.44),
needing ΔV ≈ 13 km/s. They break at half (b).

### Laps 4–5 (un-run by the gate, NOT silently dropped)

The gate caps at 3 laps (the §14 V3 floor) for the per-row wall budget. The
diagnostic over 5 laps shows the **2009–2012 phasing window** pushes EVERY row
(including Aldrin under some cadences) above the bar at laps 3–4 (high-e ~0.41–0.47
members, ΔV 6–18 km/s). Example (5.75ggF3, the half-b passer at 3 laps): laps 3–4
= 8.14 / 12.13 km/s — it would FAIL at 4–5 laps. This is a real, recorded finding:
the Russell free-return rows are higher-energy than Aldrin and do not maintain
in-family across the full §14 5-lap horizon under synodic re-phasing. The result is
also cadence/start sensitive (the #161 Aldrin test threads between the excursions
from its 2003-07-11 start + 2.135 yr cadence). 3-lap evidence is the floor reported;
the laps-4–5 excursion is the honest ceiling caveat.

## V3 reached (both halves) — recommended writebacks

Two rows clear BOTH V3 halves at the 3-lap floor:

- **aldrin-classic-em-k1-outbound** — corroborates the existing #161 V3 path; it
  is already at V2 and its V3 evidence is the dedicated #161 test. No NEW writeback
  needed beyond what #161 already justifies (this batch re-confirms it).
- **russell-ch4-5.75ggF3** — the genuinely new V1→V3 lift.

### Recommended `validation_level` writeback (for the main session to apply)

> Do NOT apply here. `data/validate.py::_LEVEL_EVIDENCE` must gain the matching
> entry FIRST (the over-claim guard refuses any V1+ row not in the registry), then
> the catalogue row may be tagged. Both edits are the main session's to make.

**russell-ch4-5.75ggF3 → V3**, with the `_LEVEL_EVIDENCE[("russell-ch4-5.75ggF3",
"V3")]` evidence chain:

```
spec §14 V3 (#158 continuation V1->V3 lift): the circular-coplanar free-return
closure (the row's #137 V1 evidence) continued out to the true ephemeris (DE440)
via the #158 homotopy ladder closes BALLISTICALLY (best-final residual 0.0051 km/s
< TOL_KMS) with emerged V_inf 5.76/8.56 km/s matching the sourced Russell Table
4.12 anchor 5.75/9.36 within the 1.5 km/s coplanar-vs-3D band, at t0 2003-07-09
inside the §14 launch-window gate (150 d of the 2003-02-09 priority date); AND the
ephemeris-mode horizon TCM over 3 laps (one E-M synodic per lap) is BOUNDED and
within the engineering plausibility bar (per-cycle dV 0.0/0.0/1.21 km/s, all <
MAINTENANCE_DV_CONVENTION_KMS = 3.0). tests/verify/test_continuation_v3_batch.py::
test_continuation_row_reaches_v3[russell-ch4-5.75ggF3]. NOTE (recorded ceiling):
the 5-lap horizon exceeds the bar at laps 3-4 (2009-2012 phasing window); V3 is
asserted at the §14 3-lap floor.
```

(If the main session prefers to also re-cite Aldrin's continuation V3 corroboration
under a separate evidence key, the chain is the same shape with the −20 d /
0.001 km/s / 6.72-9.09 numbers and the `[aldrin-classic-em-k1-outbound]` test id —
but Aldrin's V3 is already established by #161 and is at V2 in the catalogue, so no
change is required.)

## Break points (recorded, strict-xfail in the gate — never softened)

| row | breaks at | why |
|---|---|---|
| russell-ch4-5.30gGf3 | half (b) horizon TCM | lap 2 high-e (e≈0.44) member, ΔV≈12.98 km/s > 3.0 bar |
| russell-ch4-9.94Gg3  | half (b) horizon TCM | lap 2 high-e (e≈0.44) member, ΔV≈13.18 km/s > 3.0 bar |
| russell-ch4-9.353Gg2 | half (a) phase-match | no converged V∞-matched in-window DE440 basin (continues off-fingerprint) |

These stay at **V1**. The strict-xfail markers in
`test_continuation_v3_batch.py::_EXPECT` record the binding half; the day a marker
flips to PASS, the row gained the missing representability and the marker must be
removed.

## Honest ceiling statement

This lift reaches **only circular-coplanar-closable rows** (the #137 V1 set +
Aldrin) — by construction the continuation cannot seed anything else. Within that
set:

- 1 NEW row reaches V3 at the §14 3-lap floor (**russell-ch4-5.75ggF3**);
- Aldrin re-confirms V3 (already covered by #161);
- 2 rows phase-match but break the horizon-TCM bar (5.30gGf3, 9.94Gg3);
- 1 row breaks the phase-match itself (9.353Gg2);
- and even the passer fails the FULL 5-lap horizon (the 2009–2012 excursion).

The rows that don't close circular-coplanar at all (the multi-arc set, S1L1) are
out of this driver's reach entirely — they need the closer sweep / Appendix C
machinery, not the continuation. The continuation is a precursor to, not a
substitute for, the n-body harness (the V3→V4 step).
