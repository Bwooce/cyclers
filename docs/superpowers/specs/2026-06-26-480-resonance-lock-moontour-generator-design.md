# #480 — Resonance-lock moon-tour generator: Milestone 1 (IEG real-ephemeris reproduction)

**Date:** 2026-06-26
**Task:** #480 (capability lever — resonance-lock / commensurate-phasing multi-rev moon-tour generator)
**This spec covers:** Milestone 1 only — faithful reproduction of the Galilean **IEG** (Io–Europa–Ganymede Laplace-resonance) tour in **real ephemeris**, validated to its published invariants.

## Why this milestone

#473 established that single-window patched-conic closes drift out of phasing, and that even a circular-coplanar *resonance-lock* search re-finds only #339 + two known Galilean tours (IEG, EGC) — 8/10 #468 skeletons are off-basin in that model. The conclusion: **the model is the wall.** The published Galilean/Cassini tours hold together only in the real-ephemeris, resonance-locked regime. So the first worthwhile build is a real-ephemeris reproduction of one published tour, proving the generator's resonance-lock core on the cleanest case.

### Decisions locked during brainstorming (2026-06-26)
1. **Success bar:** faithful reproduction of ≥1 published tour (not discovery, not infra-only).
2. **Match target:** published **invariants** first (achievable, unblocked); bit-exact author match is a **data-gated stretch**.
3. **First target:** Galilean **IEG** (Hernandez-Jones-Jesick 2017) — already bounded in #473's crude model, canonical 4:2:1 Laplace resonance, **ballistic** (no powered DSM).
4. **Approach:** **A** — extend the proven V4-Uranus real-eph pipeline (#335) to Jupiter, warm-started from the #473 seed.

## Architecture — five bounded units (data flows left → right)

### Unit 1 — Jovian ephemeris setup
- Acquire NAIF kernels into `data/kernels/`: `de440.bsp` (planetary) + a Jovian satellite ephemeris (`jup365.bsp` or current NAIF Galilean kernel).
- Register Jupiter + Io/Europa/Ganymede/Callisto in the SPICE-backed system, mirroring #335's URA111/URA107 setup (`v4_uranus.py` is the template).
- **Interface:** given epoch + body, return real state (position/velocity) in the Jupiter-centred frame.
- **Depends on:** existing SPICE loader used by `v4_uranus.py`.

### Unit 2 — IEG seed adapter
- Input: #473's bounded IEG geometry — sequence `(Io, Europa, Ganymede, Io)`, per-leg ToFs, encounter V∞ vectors, phasing (from `scripts/_resonance_lock_moontour_473.py`, the committed bounded flip).
- Output: a machine-loadable initial guess for the real-eph corrector — encounter epochs + V∞ vectors as a flat IC vector.
- **Interface:** `ieg_seed() -> ShooterIC`. Pure, deterministic.
- **Depends on:** the #473 geometry (promote the relevant constants out of the scratch script into a small module/data file; do not import the scratch script).

### Unit 3 — Resonance-locked real-ephemeris corrector (the core new unit)
- Location: `src/cyclerfinder/nbody/moon_tour_shooter.py`, building on `nbody/shooter.py`.
- Multiple-shooting / multi-encounter corrector that converges the IEG tour in real ephemeris (DE440 + Jovian kernel) from the Unit-2 seed, enforcing three constraint groups:
  - (a) **Leg connectivity:** each leg connects moon_i → moon_{i+1}; V∞ continuity / flyby B-plane constraints at each encounter.
  - (b) **Laplace commensurability (the resonance lock):** the encounter phase configuration repeats each cycle (the constraint single-window closes lack).
  - (c) **Bounded inter-cycle drift:** the #473 relative criterion (`max per-cycle drift / rendezvous-moon SMA ≤ ~0.91`, R_REF calibrated from #339).
- Uses the **analytic-STM Jacobian** (already built; the FD Jacobian was the #388-saga bottleneck) for tractability; moon legs are days-scale so each shoot is cheap.
- **Interface:** `converge_ieg(seed: ShooterIC, system) -> ConvergedTour | NonConvergence`.
- **Fallback rung (built in):** if direct convergence stalls, homotopy continuation from circular-coplanar IEG → real-eph (ramp the ephemeris-perturbation strength 0→1).

### Unit 4 — `v4_jupiter` gauntlet
- Location: `src/cyclerfinder/data/validation/v4_jupiter.py`, an analog of `v4_uranus.py`.
- Runs the converged tour through **V3** (independent n-body, REBOUND IAS15 / DE440) and **V4** (high-fidelity real-eph with the Jovian SPICE kernel), emitting per-cycle drift, integrator agreement, the verdict, and the measured `dv_band`.
- **Interface:** `run_v4_jupiter(tour, system, n_cycles) -> V4Verdict`. Reuses the V3/V4 verdict dataclasses + framework.

### Unit 5 — Invariant golden (the success gate)
- Location: `tests/` (e.g. `tests/verify/test_ieg_real_eph_reproduction.py`).
- Compares the converged tour's invariants against **Hernandez-Jones-Jesick 2017** published values: the Laplace resonance ratios (4:2:1 Io:Europa:Ganymede), the encounter V∞ levels, and the tour period.
- **Sourcing discipline (non-negotiable):** the EXPECTED side traces to the paper (sourced golden), never to a value our own code computed. Tolerances documented and justified.

## Data flow
```
#473 IEG geometry
  → Unit 2 seed adapter (epochs + V∞ IC)
  → Unit 3 real-eph corrector (jup365 + de440, analytic STM, resonance + drift constraints)
  → converged IEG tour
  → Unit 4 v4_jupiter gauntlet (V3 + V4, per-cycle drift + dv_band)
  → Unit 5 invariant check vs Hernandez 2017
  → reproduction verdict
```

## Staging
- **M1 (this spec):** IEG *invariant-match* reproduction — V3/V4 pass AND invariants match Hernandez 2017 within documented tolerance.
- **M1-stretch (data-gated):** bit-exact author-trajectory match, ONLY if Hernandez 2017's detailed trajectory data is acquired (arXiv source / tables / contact). Expected to hit the publication gap; tracked as an optional follow-on, not a blocker on M1.
- **Out of scope (separate specs):** the discovery enumerator (Approach B → #467), EGC reproduction, Cassini powered tours (#465/#464 powered legs), cross-system continuation (Approach C), 3D-inclined tours (#451).

## Error handling / risks
- **Shooter non-convergence** — mitigated by the #473 warm start + analytic STM + days-scale legs; fallback is the homotopy rung in Unit 3. If BOTH fail, that is a characterized negative (the IEG tour does not converge in real-eph from this seed) — report honestly, do not force.
- **Kernel acquisition** — NAIF `de440`/`jup365` are public and free; if a needed kernel is unavailable, document and stop (do not fabricate ephemeris).
- **Publication gap on exact-match** — expected; M1's invariant-match is the achievable bar; exact-match is explicitly the gated stretch.
- **Catalogue:** M1 produces at most a *known-reproduction* row (V4 ceiling — a published tour cannot reach V5 by definition). Do NOT self-admit; report the reproduction for separate human admission.

## Testing
- **Unit 5 golden invariant test** (Hernandez 2017 sourced) — the success gate.
- **V3/V4 gauntlet tests** — bounded per-cycle drift, integrator agreement (reuse the V-framework patterns from `test_v4_uranus.py`).
- **Seed-adapter parity** — the #473 geometry round-trips through Unit 2 unchanged.
- **Bounded-drift positive control** — the converged IEG reproduces its bounded relative-drift signature (≤ R_REF).
- Full `tests/data tests/search` ratchet stays green (additive units; no change to existing behaviour).

## Scope boundaries (YAGNI)
- M1 is **IEG only**, **ballistic** (Laplace is ballistic — no powered DSM; powered deferred to Cassini).
- **No discovery enumerator** (known sequence, seeded directly).
- **No cross-system continuation.**
- **No catalogue self-admission** (V4-ceiling known-reproduction → human admission).
