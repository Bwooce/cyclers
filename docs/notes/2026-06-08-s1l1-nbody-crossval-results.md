# S1L1 independent n-body cross-validation — VERDICT: DRIFT (#165 / #94)

**Date:** 2026-06-08
**Experiment:** the INDEPENDENT n-body cross-check of the #164 two-body S1L1
closure — the decisive gate for whether `russell-ch4-4.991gG2` (#94) earns the
catalogue's first V3. #164 found a CLOSED two-arc geometry with the **two-body
radial-crossing continuation solver** (`search/continuation_chain.py`). This pass
confirms — or refutes — that geometry with an **independent integrator**:
REBOUND / IAS15 over the real DE440 planet ephemeris (`nbody/propagator.py`),
sharing NONE of the continuation solver's machinery.

**Code (this experiment):** new test
`tests/nbody/test_s1l1_nbody_closure.py` (+ a thin in-test construction helper; no
new `nbody/` module needed, none added). Reuses by import (no edit):
`nbody/propagator.py` (RestrictedNBody / RailsEphemerisCache), `nbody/forces.py`,
`core/ephemeris.py` (astropy/DE440), `core/kepler.py`, and the #164
`search/continuation_chain.py` for the initial state ONLY. No edit to any nbody
module, `data/catalogue.yaml`, `docs/spec.md`, or `search/`.

**GOLDEN/HONESTY.** EXPECTED = the row's OWN SOURCED anchors (Russell 2004
Table 4.9): V_inf **E 4.99 / M 5.10** km/s. The n-body miss-distance and V_inf are
EVIDENCE, never imposed. The integrator is independent of the continuation solver —
that independence is the entire point; nothing downstream of the initial state is
fit. **No catalogue writeback** — S1L1 promotion is the main session's call after
this verdict.

---

## VERDICT: DRIFT — S1L1 is two-body closeable but n-body-UNCONFIRMED

Propagating the #164 closed geometry forward in the independent integrator does
**not** reproduce a Mars encounter. Both arcs miss the real DE440 Mars by **~2.6 AU**:

> **The #164 construction is a radial-crossing + V_inf + ToF closure that does NOT
> enforce LONGITUDE RENDEZVOUS with the real ephemeris Mars.** The spacecraft
> crosses Mars's orbital RADIUS at the right time with the right speed, but ~110°
> away in heliocentric longitude from where DE440 Mars actually is at the crossing
> epoch. The "V_inf at Mars = 5.10" that the two-body solver reports is a
> radius-and-speed match, not a real relative velocity at a real encounter.

This is the honest "two-body close but n-body drifts" outcome the task names as a
valid result. S1L1 is **not** independently confirmed; **no V3 writeback is
recommended.**

### Per-encounter table (achieved vs sourced)

Closed geometry re-solved from #164 (winning rung nstep=9; matches the #164 note
exactly): arc1(g) a=1.34479 / e=0.26610 / n_rev=0, arc2(G) a=1.28701 / e=0.24089 /
n_rev=1, t0 = J2000 + 239.42 d.

| arc | Mars crossing epoch | n-body miss | n-body V_inf | sourced V_inf | converged | energy drift | 2-body miss (cross-check) | #164 own V_inf E/M |
|---|---|---|---|---|---|---|---|---|
| arc1 (g) | J2000 + 429.3 d | **2.681 AU** | 38.80 km/s | 5.10 | True | 0.0e+00 | 2.661 AU | 4.988 / 5.103 |
| arc2 (G) | J2000 + 427.8 d | **2.616 AU** | 38.51 km/s | 5.10 | True | 0.0e+00 | 2.593 AU | 4.984 / 5.100 |

(The n-body "V_inf" of ~38 km/s is meaningless as a flyby speed — it is the bare
kinematic velocity difference between the spacecraft and a Mars that is ~110° away,
recorded only to show how far off the encounter is.)

### The tolerance used, and why it cannot be the problem

**Encounter band: miss < 3 Mars SOI ≈ 0.0116 AU**, V_inf within 0.5 km/s of the
anchor (the same band the #164 continuation gate uses for its V_inf half).

`r_SOI(Mars) = a_Mars (m_Mars/m_Sun)^(2/5) ≈ 0.00386 AU ≈ 5.8e5 km`. A real targeted
flyby threads the B-plane well *inside* one SOI; 3 SOI is a deliberately generous
confirmation band. The integrator is trustworthy FAR below this scale — the harness
golden gates pin REBOUND/IAS15 to the analytic two-body solution to **< 1 km** over
120 d (GATE 1), the DE440 planet read to **< 1e-6 km** (GATE 3), and Sun-only energy
drift to **< 1e-10** (GATE 2). The observed ~2.6 AU miss is **~230×** the (already
generous) tolerance and **~4500×** one Mars SOI — it cannot be loosened away, and it
is not integrator noise. The tolerance was NOT touched to chase a CONFIRMED.

### It is structural, not an integrator artifact

A pure two-body (Sun-only) Kepler propagation of the **same absolute initial state**
returns the **same ~2.6 AU miss** (2.661 / 2.593 AU vs the n-body 2.681 / 2.616 AU —
agreement within 0.1 AU; the ~0.02-AU difference is the real Earth+Mars perturbation
over the cruise). The gap is therefore a property of the construction's geometry, not
of REBOUND. An epoch scan over a full Earth–Mars synodic period (>820 d, departure
stepped 5 d) finds the *best achievable* real-Mars miss for this arc shape is still
**0.234 AU** (~38 Mars SOI) — i.e. even the optimally-phased departure of this
geometry never reaches a real Mars flyby, and the #164 closed t0 = 239 d is not even
that best-phased point.

---

## Why the two-body solver closed but the real cycler does not exist

The #164 residual constrains, per arc: emerged V_inf at Earth and Mars vs the sourced
anchors, and the Earth-to-Earth arc ToF vs the descriptor ToF. The Mars term is
evaluated against an **effective planet** placed by `_planet_phase` at a *mean-motion
phase from the shared epoch* — used only to set Mars's eccentric **radius/speed**, NOT
its true DE440 **longitude**. So the solver is free to match `(radius, speed, ToF)`
while the spacecraft's actual crossing longitude and the real Mars longitude are
unrelated. The shared epoch `t0` is a lever over Mars's radius (perihelion-vs-aphelion,
the eccentricity physics that closed the ToF gap in #164), but it carries **no
longitude-rendezvous constraint**. The closure is real in the `(radius, V_inf, ToF)`
space the residual lives in; it is empty in the real 3-D ephemeris.

This is the same lesson the harness has recorded twice before, one level up:
*representable ≠ reachable* (#157), and *radius-and-V_inf-closeable ≠ a real n-body
encounter* (the Phase-C verdict, #133). The continuation closed the ToF axis of the
#163 frontier; the **longitude/phasing axis against the true ephemeris remains open**
and is what this independent check exposes.

---

## V3-writeback recommendation

**Do NOT recommend V3 writeback for `russell-ch4-4.991gG2`.** The #164 closure is a
two-body `(radius, V_inf, ToF)` artifact; the independent DE440 integrator shows no
real Mars encounter (2.6 AU miss, both arcs). The correct status for S1L1 is
**"two-body closeable, n-body-unconfirmed"** — the blocker's caveat #3 is now
discharged with a DRIFT result, not a confirmation.

What a future CONFIRMED would require (for the main session / a later pass): a
construction that adds a **longitude-rendezvous constraint** against the real DE440
Mars at the encounter epoch (i.e. solve for the spacecraft state that actually
intercepts Mars's true position, e.g. a Lambert/multiple-shooting node on the real
ephemeris), then re-run this exact independent gate. Only if that brings the n-body
miss inside the 3-SOI band at V_inf ≈ 5.10 does S1L1 graduate to a V3 evidence chain.
The test (`test_s1l1_nbody_crossvalidation_verdict`) currently PINS the DRIFT verdict;
its asserts must be deliberately inverted if such a construction ever lands inside the
band.

---

## Test results (verbatim)

```
$ uv run pytest tests/nbody/test_s1l1_nbody_closure.py -m "slow or not slow"
2 passed in 5.75s
```

Tests: `test_mars_soi_tolerance_is_grounded_not_inflated` (pins the tolerance
derivation), `test_s1l1_nbody_crossvalidation_verdict` (@slow; pins the DRIFT
verdict — both arcs miss real Mars by > 1 AU, n-body and two-body agree within
0.1 AU, integrator converged with energy drift < 1e-6).

**No harness regression.** The nbody golden gates and the full non-slow nbody suite
stay green on the committed state:

```
$ uv run pytest tests/nbody/test_golden_*.py -m "slow or not slow"   -> 12 passed
$ uv run pytest tests/nbody/ -m "not slow"                           -> 34 passed
```

`ruff check` / `ruff format --check` / `mypy src tests` all clean.

**No catalogue writeback.** The emerged n-body geometry is evidence; the sourced
anchors remain the EXPECTED target. S1L1 promotion (here: NON-promotion) is the main
session's call.
