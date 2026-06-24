# #435 High-e Sun-planet ER3BP seeds — verdict

**Date:** 2026-06-24
**Status:** COMPLETE — capability delivered (CR3BP seed generator at arbitrary μ);
structural NEGATIVE (no novel e>0-only ER3BP family at the probed high-e
Sun-planet seeds). Capability-extends the #432 Earth-Moon-only negative to the
high-departure regime #432 could not reach. No catalogue writeback.

## What was run

#432's ER3BP discovery campaign was seed-limited to Earth-Moon ICs and so could
NOT probe the systems where the ER3BP departs most from the CR3BP — the high-e
Sun-planet systems. #435 closes that gap. `src/cyclerfinder/search/cr3bp_seed_generator.py`
GENERATES rotating-frame CR3BP seeds at arbitrary μ (L1-Lyapunov from the
collinear linearization, refined by the fixed-Jacobi symmetric corrector; DRO
from a retrograde co-orbital guess), and `scripts/run_435_high_e_er3bp.py` feeds
them to the already-built #432 pipeline (`continue_and_monitor` /
`adjudicate_trace`), continuing each to the body's **real** orbital eccentricity.

- Sun-Mercury (e=0.206), Sun-Mars (e=0.093), Sun-Pluto (e=0.249).
- 2 seeds/system (L1-Lyapunov + DRO), n_steps=60.
- Wall: 273 s.

Output: `data/er3bp_discovery_435_highE.jsonl` (6 records).

## Result

| system | e | seed | outcome | e_max | e_star | lit |
|---|---|---|---|---|---|---|
| Sun-Mercury | 0.206 | L1-Lyapunov | survives | 0.206 | None | not-found |
| Sun-Mercury | 0.206 | DRO | survives | 0.206 | None | not-found |
| Sun-Mars | 0.093 | L1-Lyapunov | survives | 0.093 | None | published |
| Sun-Mars | 0.093 | DRO | survives | 0.093 | None | published |
| Sun-Pluto | 0.249 | L1-Lyapunov | survives | 0.249 | None | inconclusive |
| Sun-Pluto | 0.249 | DRO | survives | 0.249 | None | inconclusive |

**Overall: survives=6, dies=0, bifurcates=0.** Every Lyapunov and DRO family
continues smoothly from e=0 to the body's real eccentricity staying in its
stability regime — **no bifurcation** (`e_star=None`), so no e>0-only family
branches off at any of the three high-e Sun-planet systems.

## Reading (honest scope of the negative)

Combined with #432, the (eccentricity-continuation + Floquet-transition) method
now finds **no novel e>0-only ER3BP cycler** across:

- Earth-Moon (e=0.0549, extended to 0.15) — Broucke 7P + 2 Koblick NRHOs (#432);
- the three real high-e Sun-planet systems (Mercury/Mars/Pluto) — L1-Lyapunov +
  DRO (#435).

This is the regime #432's verdict flagged as "the most promising for e>0-only
structure, not reached". It is now reached, and it is **empty under this method**.

The conditions of the negative remain (per #432, two now retired, one still open):

1. ~~Seed-limited to Earth-Moon~~ — **RESOLVED.** High-e Sun-planet systems are
   now probed at real eccentricity.
2. **Seed families are libration (Lyapunov) + co-orbital (DRO)**, which are
   *adjacent to*, not strictly, cycler-class — these are where the published
   ER3BP literature actually lives, so they are the right frontier to probe, but
   a genuinely cycler-class high-e seed (a resonant heliocentric family) is still
   untested.
3. **Method = bifurcation-along-continuation.** A genuine e>0-only family with
   *no CR3BP limit at all* requires **direct e>0 seeding** (task #436); a flagged
   bifurcation would require **branch-switching** to trace (task #436). Neither
   was exercised (there were no bifurcations to switch on).

So the precise claim: *no e>0-only family branches off the L1-Lyapunov or DRO
families of Sun-Mercury, Sun-Mars, or Sun-Pluto up to their real eccentricities.*

## Correctness note caught during harvest (would have under-reported coverage)

The first #435 run reported Mercury & Mars **Lyapunov as `dies@e_max=0.0000`** —
the continuation failed at the *first* e-step (e≈0.002–0.003) with the corrector
plateauing at err ~1e-7, just shy of the 1e-10 tol after 60 iterations. Diagnosis
(not a physical death): the smallest-first amplitude ladder in `lyapunov_seed`
landed a Lyapunov orbit so close to L1 at μ~1e-7 (ydot0≈6.6e-4) that the
(x,ẏ)/(y,ẋ) ER3BP corrector is ill-conditioned. A modestly larger small-amplitude
Lyapunov (amplitude=3e-3 → ydot0≈3e-2, still deep in the linear regime) is
well-conditioned and continues cleanly at Mercury/Mars/Pluto alike. The runner
now requests amplitude=3e-3 (sourced comment in `_generate_seeds`); Pluto, which
already converged at the default, is unaffected. Both previously-inconclusive
seeds now resolve to `survives` — turning 2 inconclusive into 2 genuine negatives.

## Disposition

- No catalogue writeback (no ER3BP V0–V5 gauntlet; matches #432/#430/#307).
- Method-versioned negative registered in `data/empty_regions.jsonl`
  (capability-extends `er3bp-discovery-em-broucke-koblick-2026-06-24`).
- Remaining frontier (task #436): direct-e>0 seeding for families with no CR3BP
  limit, and branch-switching for any future flagged bifurcation. The pipeline is
  ready to consume both.
