# #318 multi-axis joint search — Phase 2 blocker diagnosis + reframe

**Date:** 2026-06-30. Picked up #318 as the top open discovery probe ("strategic keystone —
no published paper has done a joint powered × multi-rev × 3D × epoch-locked sweep"). Phase 1
(`2026-06-16-318-multi-axis-joint-search-phase1.md`) shipped the substrate and concluded
"the next phase must adopt smarter sampling." **That framing is incomplete: the real Phase-2
blocker is a MODEL mismatch, not the sampling strategy.** This note diagnoses it and reframes
the path. No catalogue change.

## The real blocker: the four axes do not share a state model
The Phase-1 substrate composes four axis drivers, but they live in **two incompatible models**:

| Axis | Driver | Model / state |
|---|---|---|
| A powered maintenance | `low_thrust_cycler_search` (#309) | heliocentric patched-conic Lambert cycler, real planet ephemeris |
| B multi-rev Lambert | (same #309 driver, per-leg revs) | same heliocentric Lambert model |
| D epoch-locked window | `close_epoch_locked` (#289) | same heliocentric, real-ephemeris epoch |
| **C 3D / broken-plane** | `correct_general_periodic_3d` (#291) | **CR3BP nondimensional** (x,y,z,ẋ,ẏ,ż,T), single-system rotating frame |

Axes A/B/D are the **heliocentric Lambert cycler**; axis C is a **CR3BP periodic orbit**.
They have no shared state representation — a heliocentric V-E-M cycler is not a CR3BP orbit,
and the CR3BP corrector's IC `(x,y,z,ẋ,ẏ,ż)` in a rotating Earth-Moon (or Sun-planet) frame
does not embed in the Lambert tour. **This is exactly why Phase 1 could only RECORD the z0
amplitude and never drove the 3D corrector inline** — the Phase-1 note attributed this to a
"contract" deferral, but it is structural: you cannot pipe `correct_general_periodic_3d` into
the Lambert sweep because there is nothing to pipe it *to*. The Phase-1 EM probe's finding
that "the 2D Lambert engine converges identically on the 3D-request and planar cells" is the
same fact: the z0 axis is inert in the Lambert model.

## Grounding (read-only, this session)
- `correct_general_periodic_3d` signature confirms a pure CR3BP corrector (`system.mu` + 6D
  nondim state + period) — no heliocentric/ephemeris coupling.
- The substrate runs on heliocentric tours (EM reproduces; a `V-E-M-V` cell emits a
  survivor) but with **default params closes 0 cells** — the A×B×D sweep needs per-tour,
  per-published-family config (ToF guesses/bounds, synodic pair, closure body), exactly as
  the Phase-1 EM probe hand-supplied for the Aldrin tour.
- The catalogue's VEM rows (Jones 2017 family, etc.) carry **no per-leg ToF bounds** and
  have complex multi-encounter topologies (E-M-E-V-V-E, M-E-E-V-E-M) that don't map to the
  joint-search 2-leg interface — so a fresh-surface A×B×D sweep is a real config build, not a
  drop-in run.
- Phase-1's own result already shows the A×B×D Lambert axes are **largely redundant on EM**
  (192 cells → exactly 2 distinct ΔV: 1.2870 and 0.0000 km/s); multi-rev (B) and the recorded
  z0 (C) never opened a new basin.

## Reframe — what Phase 2 actually is
The Phase-1 "smarter sampling" recommendation optimises the wrong thing: Sobol over a
Cartesian product still cannot make the inert 3D axis do anything in the Lambert model. Two
honest Phase-2 options:

1. **Unify the model in real-ephemeris n-body (the right keystone build).** All four axes
   coexist *naturally* in a real-eph n-body trajectory: 3D geometry is intrinsic (no CR3BP
   needed), epoch-locking is the launch date, powered maintenance is the per-cycle ΔV, and
   multi-rev is the Lambert branch of each leg. The project already has the pieces
   (`nbody/shooter.py`, the Jovian/heliocentric n-body lanes, `chain_cycles`). The genuinely-
   novel "joint pocket" — if it exists — lives in the real-eph n-body manifold, not in a
   CR3BP×Lambert Cartesian product. This is a multi-week build (the honest cost).
2. **Drop axis C; run an honest A×B×D Lambert joint sweep across fresh systems** (VEM,
   broader Mars). Bounded, but the EM-redundancy prior + the individually-published axes make
   a clean-negative ("joint pocket empty") the strong expected outcome — valuable as
   empty-region mapping, not a discovery. Requires per-family config construction first.

## Verdict / standing
- **Phase 1's "smarter sampling" framing is superseded:** the binding Phase-2 constraint is
  the A/B/D-vs-C model incompatibility, not the sampler. The "joint 4-axis sweep" cannot be
  jointly 3D in the substrate as built.
- **Recommendation:** if #318 is pursued, do it as a **real-ephemeris n-body joint search**
  (option 1) — that is where all four axes genuinely co-vary and where any novel pocket would
  be. Treat the CR3BP-3D corrector (#291) as a separate lane, not a joint-sweep ingredient.
- No code change, no catalogue change, no novelty claim — a scoping correction on the keystone
  so the next Phase-2 effort builds the right thing instead of a smarter Cartesian product.
