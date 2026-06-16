# 2026-06-16 - #330 moontour V2 Phase 2 verdict on the #327 Umbriel-Oberon SILVER

## Phase 2 scope

#306 Phase 1 (commits `a866190` / `e2f1dac` / `71b44f3`) shipped the V1 + V2
gauntlets for full-3D CR3BP periodic orbits — same-model closure (V1) +
long-span bounded-drift propagation (V2) of a 6D `state0 + period` IC.

The #327 Umbriel-Oberon-Umbriel SILVER (`repeated-moon-uranus-00000041`)
is NOT a periodic CR3BP IC. It is a single-cycle Lambert moontour described
by `(sequence, V_inf-per-encounter, leg-ToFs, rel_offset_deg)`. The Phase-1
V2 driver `run_v2_3d` correctly skipped it as `NOT_RUN_PHASE_1` because no
`state0 + period` exists to propagate.

Phase 2 builds the right V2 for a moontour:
`src/cyclerfinder/data/validation/v2_moontour.py::run_v2_moontour` re-solves
the SILVER's Lambert legs over `n_cycles` consecutive cycles, advancing the
moon ephemerides through their natural Keplerian motion across each cycle.
Verdict gates:
1. every cycle's Lambert converged;
2. per-cycle V_inf-continuity residual <= `V2_MOONTOUR_CLOSURE_FLOOR_KMS`
   (0.05 km/s, matches the #285 / #312 SILVER closure gate);
3. inter-cycle rendezvous drift <= `V2_MOONTOUR_DRIFT_FLOOR_KMS` (50,000 km,
   mirrors the Phase-1 same-model floor).

## Phase 1 V2 vs Phase 2 V2 — the distinction

| Phase | Defining model | V2 input | V2 propagation |
|-------|----------------|----------|----------------|
| Phase 1 (`v2_3d.py`) | full-3D CR3BP rotating frame | `state0` (6D nondim) + `period_nondim` | propagate 6D state over `n_cycles * period` |
| Phase 2 (`v2_moontour.py`) | planet-frame Lambert + circular-coplanar moon ephemeris | `sequence + V_inf-tuple + leg-ToFs + rel_offset_deg` | re-solve Lambert legs x `n_cycles`, moons advance Keplerian |

Both verdicts use the same `n_cycles_min = 3` and the same `drift_floor_kms
= 50,000`. The drift metric is the equivalent in each model: 6D rotating-
frame position deviation at cycle boundaries (Phase 1) vs planet-frame
final-encounter position offset cycle-k vs cycle-0 (Phase 2).

## SILVER V2 verdict at `n_cycles = {3, 5, 10}`

Stored SILVER fields (from `data/silver_327_verified.jsonl`):

```
sequence        = ('Umbriel', 'Oberon', 'Umbriel')
V_inf (km/s)    = (0.9199258810725036, 0.9604309791298091, 0.8946936085078939)
leg ToFs (days) = (14.940560615336594, 14.940560615336594)   # cycle = 29.88 d
rel_offset (deg) = 180.0   # gate-passing basin floor (#327)
phase0 (deg)    = 29.999999999999996
n_rev           = (1, 1)
```

V2-moontour verdict at the three `n_cycles` scales:

| n_cycles | passes_v2 | n_completed | max_drift (km) | max_closure (km/s) |
|----------|-----------|-------------|----------------|--------------------|
| 3        | **False** | 3 / 3       | 5.16e5         | 1.23e-1            |
| 5        | **False** | 5 / 5       | 5.16e5         | 2.09e-1            |
| 10       | **False** | 10 / 10     | 5.30e5         | 3.49e-1            |

Per-cycle trace at `n_cycles = 10`:

```
cycle 0: drift=0.000e+00 km, residual=2.523e-02 km/s, converged_legs=2/2
cycle 1: drift=3.264e+05 km, residual=7.505e-02 km/s, converged_legs=2/2
cycle 2: drift=5.155e+05 km, residual=1.230e-01 km/s, converged_legs=2/2
cycle 3: drift=4.876e+05 km, residual=1.680e-01 km/s, converged_legs=2/2
cycle 4: drift=2.546e+05 km, residual=2.092e-01 km/s, converged_legs=2/2
cycle 5: drift=8.561e+04 km, residual=2.461e-01 km/s, converged_legs=2/2
cycle 6: drift=3.898e+05 km, residual=2.784e-01 km/s, converged_legs=2/2
cycle 7: drift=5.299e+05 km, residual=3.062e-01 km/s, converged_legs=2/2
cycle 8: drift=4.471e+05 km, residual=3.296e-01 km/s, converged_legs=2/2
cycle 9: drift=1.761e+05 km, residual=3.491e-01 km/s, converged_legs=2/2
```

Verdict label: **`FAIL_QUASI_BOUNDED`** — the SILVER fails the strict V2
gates (drift > 50,000 km from cycle 1; closure residual > 0.05 km/s from
cycle 1) but every cycle's Lambert geometry converges over 10 cycles, and
the closure residual stays below 0.5 km/s (within the v4.7 quasi-cycler
admission envelope but FAR above the strict V2 bar).

### Reading the drift trace

The drift OSCILLATES rather than growing monotonically (cycle 5 returns
to 8.6e4 km — close to cycle 0). This is the signature of a
**near-resonant tour** in a phase that does NOT exactly lock. The SILVER's
cycle period (29.88 d) sits at 4.991 × the Umbriel-Oberon synodic period
(5.987 d) — close to 5:1 synodic but not exactly so. The 0.9% offset from
exact 5:1 synodic resonance accumulates as inter-cycle drift; the rendezvous
geometry recovers when accumulated phase mod 2π returns near 0 (around
cycle 5).

The closure residual grows monotonically (0.025 -> 0.349 km/s over 10
cycles) because each cycle the moons are at slightly different relative
phases than cycle 0, so the Lambert solution at the SILVER's fixed
`(rel_offset_deg, phase0_deg, ToFs)` no longer hits the same V_inf-
continuity envelope.

The Lambert legs continue to converge at every cycle — the candidate is
not topologically broken; it is dynamically *near-cycler-like* but not a
genuine cycler at this phasing under a fixed-time-grid V2-ballistic check.

## What 7-of-9 looks like

Gates the SILVER has cleared (from `data/silver_327_verified.jsonl`):

1. **closure** (`closure_pass_two_moon`): true, 0.025 km/s
2. **independent cross-check** (DOP853 at rtol=atol=1e-12): max dr 1.6e-5 km
3. **physical sanity** (5 deg useful-bend floor): true at all 3 encounters
4. **literature fresh** (offline corpus): not-found
5. **ML flagger** (p_fp = 0.59 vs SILVER ceiling 0.75): pass
6. **V1** (Phase 1 — same-model closure): not applicable to moontour (no
   periodic IC); the cycle-0 V_inf-continuity self-consistency residual
   (0.025 km/s) is the closest analogue and matches the SILVER's stored
   value tightly. There is no `passes_v1` boolean for the moontour signature.
7. **V2** (this work — moontour bounded drift): **FAIL** at strict gates;
   `FAIL_QUASI_BOUNDED` at the v4.7 quasi-envelope.

Still gating:

8. **V3** (#331 future task — real-ephemeris 6D nbody / REBOUND)
9. **V4** (#332 future task — HFEM Uranus via GMAT + SPICE Uranus kernels)
10. **Heaton-Longuski / #329 lit-fresh deep dive** (still required pre-admit)

## Phase 3 recommendation

The SILVER's V2 failure is **drift-and-closure-both** but it is NOT
topologically broken (every Lambert converges over 10 cycles, closure < 0.5
km/s). Phase 3 (V3 real-eph 6D nbody / REBOUND) is the right next gate:

* V3 will tell us whether the bounded-drift behavior survives real
  ephemeris (no circular-coplanar idealization). A real-eph V3 PASS would
  promote the candidate to a quasi-cycler per v4.7; a real-eph V3 FAIL
  would retire it.
* V4 (HFEM Uranus via GMAT) is also relevant but more expensive and
  premature — the cleaner question right now is whether real-eph V3 sees
  the same bounded oscillation or rapid divergence. If V3 fails, V4 is
  moot.

Phase 4 (#332) HFEM Uranus is the appropriate next step ONLY if V3 passes
with the v4.7 quasi-cycler bound (closure < 0.5 km/s under real eph over
3 cycles). If V3 fails, the SILVER is honestly retired to the negative-
results registry as a near-resonant tour that does not survive ephemeris
widening.

## Discipline anchors

* NO catalogue writeback. V2 fail / quasi-bounded does NOT admit.
* The V2 verdict is whatever the math says. The SILVER's drift OSCILLATES
  at ~5e5 km (bounded but well above the 50,000 km floor) and closure
  grows monotonically (well above the 0.05 km/s floor). Verdict is FAIL.
* `feedback_orbit_closure_discipline`: clean negative is success.
* Lambert geometry uses the same kernel #327 closed under (self-consistent
  cycle-0 residual reproduces the stored 0.025232 km/s tightly).

## Artifacts

* `src/cyclerfinder/data/validation/v2_moontour.py` — the driver
* `tests/data/test_v2_moontour.py` — 13 tests (Part A)
* `scripts/run_330_silver_moontour_v2.py` — the SILVER scan (Part C)
* `data/silver_327_moontour_v2_verdicts.jsonl` — the verdict JSONL
* this note (Part C)

Run the scan via:

```
uv run python scripts/run_330_silver_moontour_v2.py
```

Verify with:

```
uv run pytest tests/data/test_v2_moontour.py tests/data/test_v1_3d.py \
    tests/data/test_v2_3d.py -x --timeout=180
```
