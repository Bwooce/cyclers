# #388 Phase 1b — descriptor-ingestion sweep of the remaining McConaghy-2005 named rows

**2026-06-23.** Following the 4.3.1-5 result, sweep the other "descriptor-published
but not-ingested" McConaghy-Russell-Longuski 2005 Table 2 cyclers, gating each on the
CANONICAL `dsm_descriptor_seed.close_row_dsm` from the start (the 4.3.1-5 lesson: do not
chase a lucky multi-start/epoch hit).

## russell-ocampo-2.5.1+0 — STAYS V0 (high-energy off-anchor collapse, as predicted)

Descriptor (McConaghy-2005 Table 2): `2 g(1-11/14, 11/14 rev, U) f(1:1, 74.919°,
∓144.069°) h(0.5, 0, U, ±15.081°) f(1:1, 74.919°, ±35.931°)` — n=2, 4 legs
(generic/full-rev/half-rev/full-rev). Repeat constraint Σt_f = 1.7857+1+0.5+1 = 4.2857 =
2×(2-1/7) — EXACT. Ingested into `free_return_arcs`.

**Anchor: E 7.8 / M 9.9 km/s — HIGH energy.** Canonical `close_row_dsm`:
- converged=False, max_residual = **32.27 km/s**, anchor_match=False, total DSM 60.8
  km/s, emerged V∞ [7.8, 0.014, 20.78, 0.001] (off-anchor collapse).

This is exactly the [[project_388_wall_energy_selective]] prediction: a HIGH-energy
cycler collapses to the low-energy basin (V∞ → ~0 / garbage), no anchor match. No
promotion; descriptor ingested as sourced data-completeness + this determination.

## Meta-conclusion: the descriptor-seeded close_row_dsm lane promotes NO #365 row

Across the descriptor-bearing #365 set:
- **High-energy** (2.5.1+0 at 7.8/9.9; S1L1 at 4.7/5.0; the V3 regressions) → off-anchor
  collapse (energy-selective wall).
- **Low-energy** (4.3.1-5 at 3.1/2.5, the best case) → anchor recovered but
  epoch/seed-fragile and 164 m/s low_maintenance, not ballistic → no gate.
- **Structural**: these are g/f/h cyclers; the `f` (full-rev, θ=2π) and `h` (half-rev,
  θ=π) legs are singular Lambert cases, so the spec §14 V1 conic izzo+gooding crosscheck
  (the path that earned the generic-pair russell-ch4 rows their V1) does not apply.

So the descriptor-ingestion + close_row_dsm path **characterizes** these rows (and
completes their sourced `free_return_arcs`) but yields **no validation promotion**. The
real promotion lever for the f/h-leg cyclers would be a §14-V1 implementation that
handles full-rev/half-rev legs (a capability build), not more per-row grinding through
the current lane. Per-row grind STOPPED here on that basis.

Descriptor-completeness still worth finishing opportunistically, but it is data hygiene,
not a reproduction advance.
