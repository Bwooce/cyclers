# S1L1 continuous-from-one-seed V4 attempt — VERDICT: PARTIAL (#169)

**Date:** 2026-06-08
**Task:** #169 — the V4 attempt for S1L1. A CONTINUOUS-from-one-seed multi-cycle
propagation that measures the actual maintenance ΔV the #167 confirmation left
unmeasured (it re-anchored v∞ at each App-C node — "Russell's per-leg reproduction
recipe"). This MEASURES; no catalogue writeback (sibling #168 owns the V3 promotion).

**Bottom line:** **PARTIAL.** Propagating ONE continuous trajectory from the first
App-C Earth departure (leg 2, 2026-12-15, v∞ = (−2.278, 5.322, 0.574) km/s) forward
through all 20 App-C nodes WITHOUT re-anchoring v∞, on the independent REBOUND/IAS15
integrator over DE440:

- In **Russell's own model** (Sun-only patched-conic, flybys instantaneous at the
  patch points) the continuous single-seed orbit **holds**: every node lands ≪ SOI
  (Earth miss ≤ 18,434 km, Mars miss ≤ 3,173 km) at the App-C v∞ (|v∞| matched to
  ≤ 2.7 m/s), and the **total maintenance ΔV is 62 m/s over 7 cycles** — UNDER the
  spec §14 V3 budget (120 m/s) and far under the 3.0 km/s/cycle bar (8.8 m/s/cycle).
  This is **CLEAN** and **STRENGTHENS** the V3 #168 is writing.
- But once **Mars's finite continuous gravity** is modelled (Mars-perturbed), the
  legs AFTER each Mars flyby (the M→E returns) **DIVERGE** (>1e8 km, integrator
  non-converged); total continuous ΔV ~122 km/s. The App-C nodes are patched-conic
  states that don't account for the continuous deflection through the Mars
  encounter, so the per-leg re-anchoring across the Mars flyby **is doing real
  work**. This **QUALIFIES** the V3.

So S1L1 is a clean, bounded-maintenance continuous orbit *in the patched-conic
model Russell published it in*, and a powered/maintained orbit once you integrate
the Mars flyby continuously. The 62 m/s number is a genuine V4-grade reproduction
of the trajectory AND its maintenance ΔV by an independent integrator+ephemeris —
the spec §14 V4 definition — with the caveat that canonical V4 names GMAT.

---

## What this measures that #167 did not

#167 reconstructs each Mars-transit leg INDEPENDENTLY from its own App-C v∞ node:
`v_sc = v_planet(DE440) + v∞_appc` at every leg start. That RE-ANCHORS v∞ at each
node, so it verifies Russell's *published per-leg states* but never asks whether a
SINGLE continuous orbit threads them — re-anchoring silently supplies whatever
velocity each leg needs.

This run keeps ONE trajectory. Between nodes the state is propagated with no v∞
re-anchoring. At each node the **flyby patch** convention is applied: the position
is snapped to the real DE440 planet (the few-thousand-km ballistic miss is ≪ SOI —
the standard patched-conic boundary, recorded as evidence) and the velocity is
rotated toward the next leg's App-C v∞ direction (the *free* flyby turn). The
**maintenance ΔV** is the part a ballistic flyby CANNOT supply: the |v∞| magnitude
change, plus any un-bendable shortfall when the required turn exceeds the
safe-periapsis maximum bend (`core.flyby.max_bend`). EXPECTED = the App-C per-leg
v∞ (sourced); the miss / v∞ / ΔV are EVIDENCE.

Helper: `cyclerfinder.search.s1l1_corrected.continuous_chain` (additive).
Test: `tests/nbody/test_s1l1_continuous.py` (1 fast mechanics + 2 @slow gates).

---

## THE GATE — Sun-only continuous (Russell's cruise model) → CLEAN

REBOUND/IAS15 Sun-only ≡ `core.kepler.propagate` to < 1 km / 120 d (GOLDEN GATE 1),
so the Sun-only continuous chain is the independent integrator reproducing the
patched-conic seed exactly. SOI: Earth 924,649 km, Mars 577,239 km. Band = 3·SOI
(same as #165/#167, never loosened).

| leg | body | miss (km) | miss/SOI | v∞ achieved | v∞ App-C | bend° | maxbend° | ΔV (m/s) |
|---|---|---|---|---|---|---|---|---|
| 3  | M | 2,725  | 0.0047 | 5.248 | 5.248 | 10.0 | 34.5 | 0 |
| 4  | E | 7,110  | 0.0077 | 5.768 | 5.767 | 30.9 | 79.9 | 1 |
| 5  | E | 11,209 | 0.0121 | 5.766 | 5.764 | 24.1 | 79.9 | 2 |
| 6  | M | 1,278  | 0.0022 | 7.693 | 7.693 |  5.6 | 18.8 | 0 |
| 7  | E | 16,416 | 0.0178 | 3.773 | 3.771 | 58.5 | 107.7 | 1 |
| 8  | E | 7,391  | 0.0080 | 3.754 | 3.752 | 79.1 | 108.0 | 1 |
| 9  | M | 439    | 0.0008 | 4.657 | 4.657 | 10.5 | 40.8 | 0 |
| 10 | E | 8,976  | 0.0097 | 5.514 | 5.512 | 32.6 | 83.0 | 1 |
| 11 | E | 10,308 | 0.0111 | 5.534 | 5.532 | 35.1 | 82.8 | 2 |
| 12 | M | 3,173  | 0.0055 | 3.198 | 3.198 | 23.2 | 64.2 | 0 |
| 13 | E | 2,475  | 0.0027 | 6.942 | 6.942 | 23.9 | 67.2 | 0 |
| 14 | E | 13,477 | 0.0146 | 6.949 | 6.946 | 18.1 | 67.1 | 3 |
| 15 | M | 2,383  | 0.0041 | 6.262 | 6.263 | 11.9 | 26.4 | 0 |
| 16 | E | 9,625  | 0.0104 | 5.141 | 5.140 | 41.3 | 87.8 | 1 |
| 17 | E | 10,355 | 0.0112 | 5.123 | 5.121 | 33.7 | 88.0 | 2 |
| 18 | M | 879    | 0.0015 | 8.046 | 8.046 |  7.4 | 17.5 | 0 |
| 19 | E | 18,434 | 0.0199 | 4.660 | 4.657 | 45.4 | 94.4 | 2 |
| 20 | E | 9,250  | 0.0100 | 4.648 | 4.646 | 95.0 | 94.5 | 42 |
| 21 | M | 383    | 0.0007 | 3.219 | 3.219 | 16.9 | 63.7 | 0 |
| 22 | E | 480    | 0.0005 | 4.600 | 4.600 |  0.0 | 0.0  | 0 |

**TOTAL maintenance ΔV = 62 m/s over 7 cycles (8.8 m/s/cycle).** Decomposition:
22 m/s of |v∞|-magnitude change (the flyby-uncreditable part) + 40 m/s from a single
un-bendable shortfall at node 20 (a real flyby constraint: 95.0° required vs 94.5°
max bend at Earth safe periapsis — a 0.5° / 40 m/s shortfall, not a tolerance
artifact). Achieved v∞ breathes the real-eph 3.2–8.0 km/s span (matches App-C, NOT
the coplanar 5.10 / Rogers 3.05 anchors). All 20 nodes ≤ 0.02 SOI — far inside band.

Comparison: **62 m/s < 120 m/s V3 budget (52% of budget); 8.8 m/s/cycle ≪ 3.0
km/s/cycle bar (0.3% of bar).** CLEAN.

---

## THE QUALIFIER — Mars-perturbed continuous → DIVERGES on post-flyby legs

Turning on real DE440 Mars as a continuous perturber (the flyby gravity Russell
models as instantaneous, now acting continuously) on the SAME continuous chain:

| leg | body | miss (km) | v∞ in | v∞ App-C | converged | note |
|---|---|---|---|---|---|---|
| 3  | M | 14,145      | 5.796  | 5.248 | yes | G-transit, in-band |
| 4  | E | 361,917,843 | 47.701 | 5.767 | **no** | **M→E return DIVERGES** |
| 5  | E | 11,967      | 5.766  | 5.764 | yes | (re-anchored dep recovers) |
| 6  | M | 6,928       | 8.459  | 7.693 | yes | G-transit, in-band |
| 7  | E | 51,631,865  | 6.510  | 3.771 | yes | **M→E return DIVERGES** |
| 9  | M | 20,236      | 5.091  | 4.657 | yes | G-transit, in-band |
| 10 | E | 231,920,600 | 27.778 | 5.512 | **no** | **M→E return DIVERGES** |
| 12 | M | 37,553      | 3.537  | 3.198 | yes | G-transit, in-band |
| 13 | E | 121,184,501 | 15.396 | 6.942 | **no** | **M→E return DIVERGES** |
| 15 | M | 9,682       | 6.933  | 6.263 | yes | G-transit, in-band |
| 16 | E | 136,608,726 | 19.993 | 5.140 | yes | **M→E return DIVERGES** |
| 18 | M | 6,650       | 8.810  | 8.046 | yes | G-transit, in-band |
| 19 | E | 40,629,048  | 3.305  | 4.657 | yes | **M→E return DIVERGES** |
| 21 | M | 40,866      | 3.529  | 3.219 | yes | G-transit, in-band |

(Earth nodes that are NOT immediately after a Mars node — the g free-return E→E
returns — stay clean, ~10,000 km, since the re-anchored departure resets each
cycle.) **Total continuous ΔV ~122 km/s** — far over both the V3 budget and the bar.

The structure is unambiguous: **the seven Mars (G-transit) arrivals stay in-band**
(6,600–40,900 km, the #167 perturbed numbers), but **every leg leaving a Mars node
diverges** (>1e7 km, several non-converged). The App-C v∞ are patched-conic states;
the continuous deflection through the finite Mars flyby is exactly what the
instantaneous-turn patch (and #167's re-anchoring) papers over. A real continuous
flight needs active maintenance through each Mars flyby — S1L1 is "powered" in
precisely this sense, and the per-leg re-anchoring of #167 was hiding that.

---

## Honesty / why this is not massaged toward CLEAN

- The CLEAN Sun-only result **depends on the patched-conic flyby patch** (snapping
  position to the planet at each node — legitimate because the miss is ≪ SOI, the
  standard conic boundary). A FULLY continuous run with NO position patch diverges
  even Sun-only (the few-thousand-km patched-conic miss compounds over the ~500-day
  legs: 18,434 km → 3.6e8 km by the next node). That is reported plainly: the
  ~10⁴-km node misses ARE a real targeting residual a flyby's B-plane targeting
  nulls for free (a direction adjustment within the SOI, not a ΔV). The 62 m/s is
  the |v∞|-magnitude + un-bendable ΔV, which a flyby genuinely cannot supply.
- An EARLIER probe that let *Lambert* pick each next departure (instead of the
  sourced App-C direction) produced pathological single-rev transfers on the
  resonant g (E→E ~530-d) legs — v∞ blowing to 18 km/s, ~373 km/s total. That
  conflates "what a single-rev Lambert picks" with "the cycler leg"; it is NOT the
  cycler and is NOT used. The faithful measure uses the sourced App-C v∞ direction
  (the actual cycler geometry); only the magnitude residual is charged.
- EXPECTED = App-C per-leg v∞ (sourced, Russell App-C #83). Achieved miss / v∞ / ΔV
  are integrator EVIDENCE. The bands are the SAME 3-SOI bands as #165/#167, never
  loosened. The DIVERGE half is reported as a first-class result, not hidden.

---

## STRENGTHENS or QUALIFIES the V3?

**Both, cleanly separated:**

- **STRENGTHENS** — in the row's defining model (patched-conic, the model the App-C
  states and the catalogue's `delta_v_kms = 0` near-ballistic claim live in), an
  independent integrator+ephemeris reproduces the continuous trajectory AND its
  maintenance ΔV at **62 m/s < 120 m/s**. That is a concrete, bounded horizon-TCM
  number the V3 row previously asserted only as a budget. #168 may cite it as the
  measured continuous horizon-TCM in the patched-conic model.
- **QUALIFIES** — the 0-ΔV / near-ballistic claim is a *patched-conic* property; a
  continuous high-fidelity integration with Mars's flyby gravity acting continuously
  is NOT self-consistent without active per-flyby maintenance (~122 km/s if
  uncorrected). The V3 evidence/caveat should record that the continuous
  self-consistency holds in the patched-conic model only, and that canonical V4 (a
  GMAT/Tudat run that targets each flyby's B-plane continuously) is the finisher
  that would convert the patched-conic 62 m/s into a fully-continuous TCM budget.

---

## Canonical-V4 / GMAT recommendation

Our REBOUND/IAS15 cross-check satisfies the spec §14 V4 *definition* ("independent
codebase + ephemeris reproduces trajectory and maintenance ΔV") for the
patched-conic model, but spec §14 names **GMAT** as the canonical V4 tool. The
finisher is a **GMAT (or Tudat/pykep) continuous run with B-plane-targeted flyby
maintenance**: seed from the same App-C leg-2 state, propagate continuously over
DE440 with full Mars (and Earth) point-mass gravity, and let an impulsive TCM null
the post-flyby targeting each cycle. That run would (a) confirm the patched-conic
62 m/s, and (b) produce the true fully-continuous horizon-TCM for the Mars-perturbed
case — the number the Mars-perturbed DIVERGE here says is non-trivial. Until then
the honest V4-grade claim is "patched-conic continuous, 62 m/s, REBOUND/IAS15
independent; Mars-perturbed continuous needs per-flyby maintenance (GMAT to
quantify)."

---

## Test results (verbatim)

```
$ uv run pytest tests/nbody/test_s1l1_continuous.py::test_continuous_chain_mechanics \
    tests/nbody/test_s1l1_continuous.py::test_continuous_sun_only_bounded_maintenance \
    -m "slow or not slow"
2 passed, 76 warnings in 3.29s

$ uv run pytest tests/nbody/test_s1l1_continuous.py::test_continuous_mars_perturbed_diverges_post_flyby \
    -m "slow or not slow"
1 passed, 19830 warnings in 417.68s (0:06:57)
```

No nbody regression — the #167 test, the S1L1 closure test, the golden gates, and
the S1L1 search test stay green:

```
$ uv run pytest tests/nbody/test_s1l1_corrected_nbody.py tests/nbody/test_s1l1_nbody_closure.py \
    tests/nbody/test_golden_twobody.py tests/nbody/test_golden_conservation.py \
    tests/nbody/test_golden_anchor.py tests/nbody/test_golden_convergence.py \
    tests/search/test_s1l1_corrected.py -m "slow or not slow"
21 passed in 6.77s
```

`ruff check` / `ruff format --check` / `mypy` clean on the two changed files.

**No catalogue writeback.** This note MEASURES; #168 owns the V3 promotion. Finding:
STRENGTHENS the V3 in its patched-conic defining model (bounded 62 m/s continuous
horizon-TCM), QUALIFIES it for high-fidelity continuous flight (per-flyby
maintenance required; GMAT to quantify).
