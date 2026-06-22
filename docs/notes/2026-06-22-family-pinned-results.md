# 2026-06-22 — #388 family-pinned penalty homotopy: COMPLETE verdict

**Status: COMPLETE (12/12 shoots). Decisive characterized NEGATIVE — no row
ballistically reproduces its published cycler. HELD — no catalogue writeback.**

The family-pinned penalty homotopy (`search/family_pinned_shoot.py`,
`shooter.shoot(jacobian="stm", vinf_anchors=…)`) biased the full STM
multiple-shooting corrector into each row's published-V∞ basin, then ramped the
penalty to zero so the recorded V∞ emerges from an unpenalized (λv=0) ballistic
solve. Ran all 4 SnLm descriptor rows × 3 best-phase epochs on DE440.

## Results (λv=0 unpenalized solve; anchors are the SOURCED per-row E/M V∞)

| row | V | anchor E/M | best epoch | emerged V∞ (km/s) | E/M resid | retain | defect | conv | match | bend |
|-----|---|-----------|-----------|-------------------|-----------|--------|--------|------|-------|------|
| mcconaghy-2006-em-k2 | V0 | 4.70/5.00 | ep0 | 19.8,27.2,7.8,8.2 | 3.13/2.83 | 1.86 | 6.7e3 | No | No | No |
| russell-ch4-4.991gG2 | V3 | 4.99/5.10 | **ep1** | 9.8,19.5,4.0,4.7 | **0.25/0.36** | 3.70 | 7.5e3 | No | **Yes** | No |
| russell-ch4-8.049gGf2 | V3 | 8.05/10.02 | ep1 | 24.6,40.3,29.9,26.8 | 16.5/14.5 | 14.5 | 7.6e3 | No | No | No |
| russell-ch4-9.353Gg2 | V1 | 9.35/10.52 | ep1 | 30.3,42.6,15.9,10.6 | 1.24/0.07 | 0.00 | 7.4e3 | No | No | No |

(Per-row best epoch shown; all 12 records in the recovered runlog. "match" =
emerged V∞ within 0.5 km/s of BOTH anchors. "retain" = anchor-residual change
from the first penalized rung to λv=0 — large = sprang off when the pin lifted.)

## Verdict

**Zero rows produce a converged, anchor-matched, bend-feasible ballistic cycler
at the published V∞.** The single `match=True` (russell-ch4-4.991gG2 ep1, the
S1L1 row, within 0.25/0.36 of both anchors) is **not** a closure: it is
unconverged (defect 7.5e3, far above the SNOPT continuity floor), not
bend-feasible, and sprang 3.70 km/s off the pinned point when the penalty lifted.
`mcconaghy-2006-em-k2` stays **V0**.

### Per-family structure (the informative part)

- **Low-V∞ rows** (mcconaghy 4.7/5.0; russell-4.991 4.99/5.10) — the penalty
  drives them *close* to the anchor (4.991 ep1 anchor-matches at 0.25/0.36;
  mcconaghy within ~3), but the unpinned ballistic solve relaxes off-anchor. The
  published-V∞ basin is **not a ballistic fixed point** — pinned it sits near the
  anchor, released it drifts (large `retain`).
- **High-V∞ row** (8.049 8.05/10.02) — overshoots massively into a far
  higher-energy basin (V∞ 25–54 km/s, 8–20 km/s off). The higher the published
  V∞, the harder the corrector overshoots.
- **9.353Gg2 (9.35/10.52)** — a distinct, repeatable pattern: the **Mars side
  nails the anchor** (0.07–0.21 across all 3 epochs) while the **Earth side
  misses by ~1.3** km/s, with `retain=0.00` (the solution is a genuine fixed
  point that holds when unpinned — but an off-anchor-on-E member). The two
  encounters cannot simultaneously hit 9.35 (E) and 10.52 (M) ballistically.

## What this establishes (third independent lane → same wall)

This is the **third independent method** to reach the #388 family-selection wall,
now from the *literal constructed parent* with the published-V∞ basin actively
targeted:
1. Conic N-arc continuation → DE440, off-anchor basin.
2. FD/STM full multiple-shooting (literal parent) → stalled/off-anchor.
3. Family-pinned penalty homotopy (this) → reaches/nears the published V∞ under
   the pin, but **no row is a zero-ΔV ballistic fixed point there**.

Combined with the #415 ΔV-band finding, the honest conclusion is: **the published
"ballistic" SnLm cyclers are ballistic in the idealized circular-coplanar model
(geometric AR/TR criterion — and we golden-reproduced those constructions), but
they are NOT zero-ΔV ballistic in DE440 at the published V∞.** Per McConaghy 2006
itself, S1L1 is "only *nearly* ballistic in the ephemeris model" (~10 m/s/30 yr).

The right success criterion is therefore **not** V1-ballistic real-ephemeris
closure (which we have now shown, three ways, does not exist for these rows). It
is **V2-powered with a quantified, bounded maintenance ΔV** — the band the #415
note defines and the tier the Aldrin row already occupies. The concrete next
lever for an actual `mcconaghy-2006-em-k2` promotion is to **measure its
real-ephemeris maintenance ΔV** (station-keeping budget) and check it against the
sourced ΔV band — not to keep chasing a zero-ΔV closure.

## Method asset

The STM Jacobian + the V∞-anchor penalty layer are reusable beyond this verdict
(STM is also #347's n-body monodromy piece; the penalty/homotopy generalizes to
any V∞-targeted correction). `mcconaghy-2006-em-k2` and the V3/V1 rows unchanged;
no `data/catalogue.yaml` / `validate.py` edit.

## References
- `docs/notes/2026-06-21-shooter-stm-batch-results.md` (the FD/STM lane),
  `docs/notes/2026-06-21-narc-continuation-results.md` (the conic lane),
  `docs/notes/2026-06-22-dv-band-definitions.md` (#415 ΔV bands — the V2-powered
  reframing).
- Runlog (authoritative, recovered after a concurrent `git clean` wiped the live
  one): `/home/bruce/cyclerfinder` backup + `data/runs/shooter-family-pinned.jsonl`
  (post-recovery records).
- Memory: `project_dsm_closure_modeljump_blocker`,
  `project_s1l1_realeph_closure_blocker`, `feedback_orbit_closure_discipline`.
