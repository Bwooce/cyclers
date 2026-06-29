# #480 follow-up — EGGIE Stage 2 verdict: right basin, FD-Jacobian closure wall

**Date:** 2026-06-29
**Task:** #480 follow-up 1, Stage 2 (ideal-model n-body close).
**Plan:** `docs/superpowers/plans/2026-06-29-480-eggie-resonant-conic-generator-plan.md`
**Builds on:** Stage 1 (`535d2fb`, resonant-conic generator) + Stage 2 infra
(`38fc0a5`, `nbody/jovian_ideal.py` ideal-model corrector lane).

## One-line verdict

The resonant-conic seed is in the **right basin** (all three V∞ on the Table-4
targets), and the ideal-model multiple-shooting corrector pulls it from ~10 km/s
seed defects down to **~0.1–0.2 km/s leg/wrap velocity continuity** — a large
advance over the M1 off-basin negative (5.9 km/s) — but it **plateaus there and
does NOT reach ballistic closure**. The barrier is the **FD-Jacobian noise floor
on a flyby-sensitive multi-rev tour**, the same wall recorded in
`memory/project_dsm_closure_modeljump_blocker.md`; the known lever (an analytic
STM Jovian Jacobian) is not yet built for the Jovian lane.

## What was tested

`ideal_eggie_shoot` (`nbody/jovian_ideal.py`): the EGGIE resonant-conic guess
(Stage 1) → periapsis-node seed → multiple-shooting in the ideal circular-coplanar
Jovian n-body model (reuses `jovian_defect_residual` / `JovianRestrictedNBody`
with an injected ideal ephemeris; moons = real GM/radius on ideal circular orbits).

## Numbers (`scripts/_stage2_decomp_480.py`, ideal model)

| quantity | patched-conic SEED | after trf (25 iter) |
|---|---|---|
| total residual \|r\| | 8.69e5 | 3.39e2 |
| leg-continuity block | 8.69e5 | 2.82e2 |
| periodicity-wrap block | 5.99e2 | 1.88e2 |
| leg position defects | 2.4e4 – **8.6e5 km** | 2.5 – 54 km |
| leg velocity defects | 0.20 – **10.43 km/s** | **0.10 – 0.19 km/s** |
| wrap velocity | 0.47 km/s | 0.19 km/s |
| correction ΔV | — | 4101 m/s |
| converged (leg < floor) | — | **False** |

Per-encounter V∞ at the plateau (km/s, tgt): Europa 10.27 (9.12), Ganymede 6.91
(7.07), Ganymede 7.07 (7.07), Io 8.60 (8.38), Europa 10.19 (9.12). The two
Ganymede and the Io V∞ stay in band; Europa drifts ~1.1 km/s up as the corrector
trades V∞ accuracy for continuity — a signature of a local plateau, not closure.
lm and trf agree on the plateau (lm ~5.1e2, trf ~3.2–3.4e2); more iterations do
not break it (lm crawled 900 s / nfev≈750 without escaping).

## Why it plateaus (diagnosis)

1. **Patched-conic→n-body flyby gap (seed quality).** The seed places nodes at the
   moon-flyby periapsis (deep in the moon gravity well). Integrating the patched-
   conic node velocity through the real moon gravity deflects the spacecraft by
   tens-to-hundreds of thousands of km within a single leg (worst: Io→Europa,
   8.6e5 km). The seed is far from any continuous n-body trajectory.
2. **FD-Jacobian noise floor (corrector).** `scipy.least_squares` uses a finite-
   difference Jacobian. On a high-eccentricity (e≈0.62) multi-rev arc that grazes
   moon wells, the residual is extremely sensitive to the perijove state, so the
   FD Jacobian is noisy and the Gauss-Newton step stalls at ~0.1–0.2 km/s — far
   above the 1e-6 km/s continuity floor. This is precisely the
   `project_dsm_closure_modeljump_blocker` finding (heliocentric multi-rev cyclers:
   "FD Jacobian dominates"; the fix there was the analytic block-bidiagonal STM
   Jacobian, ~8.6× faster/iter, which made the shooter tractable and closed the
   literal-parent batch).
3. Soft caps: `ideal_eggie_shoot`'s `max_nfev`/`max_wall_sec` are not enforced
   inside the scipy loop (the same soft-cap bug noted in the M1 verdict follow-up
   #2) — `lm` ignores `max_nfev` entirely; use `method="trf"` for a bounded run.

## Status of the #480 reproduction

- **Stage 1 (resonant conic): SUCCESS** — a single conic (a from the 4:5
  resonance, e≈0.62) puts Io/Europa/Ganymede V∞ on the Table-4 targets. The
  off-paper-basin problem of the M1 verdict is SOLVED at the seed level.
- **Stage 2 (n-body close): PARTIAL / characterized negative** — right basin,
  ~0.1–0.2 km/s from closure, plateaued by the FD-Jacobian wall. EGGIE is NOT yet
  reproduced as a ballistic n-body cycler.
- No catalogue change. The golden (`tests/verify/test_ieg_reproduction_golden.py`)
  stays skipped; this note supersedes the M1 verdict's failure mode (it is now an
  in-basin closure wall, not an off-basin relaxation).

## The unblock (next lever, not a re-run)

The corrector, not the seed *family*, is now the wall. To close EGGIE:
1. **Analytic STM Jacobian for the Jovian lane** (port the heliocentric block-
   bidiagonal STM to `jovian_shoot`/`ideal_eggie_shoot`) — the documented fix for
   exactly this FD plateau. Highest-leverage.
2. **Sub-arc the legs** (more shooting nodes per leg) so each defect is small and
   the perijove sensitivity is split across nodes.
3. **Gravity homotopy** (ramp moon GM 0→full): at GM=0 the Lambert seed closes
   exactly (2-body); ramp the flyby gravity on gradually, re-converging.
4. Enforce the `max_nfev`/`max_wall_sec` caps inside the residual callback.
Only after an in-band ballistic *ideal* close does Stage 3 (ideal→real homotopy
+ real-eph confirm) become worthwhile.

Scratch: `scripts/_stage2_*_480.py` (not committed).
