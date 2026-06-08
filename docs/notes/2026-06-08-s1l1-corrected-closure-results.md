# S1L1 corrected-topology closure — VERDICT: CONFIRMED (#167 / #94)

**Date:** 2026-06-08
**Task:** #167 — build the CORRECTED-topology S1L1 closure and validate it in n-body,
the redemption of #164/#165 now that #166 found both the topology error and the exact
real-eph seed.

**Bottom line:** **CONFIRMED.** The corrected per-cycle sequence
`E → g(Earth-Earth free return, NO Mars) → E flyby → G(Earth-Mars-Earth transit,
true longitude rendezvous with DE440 Mars) → E`, seeded entirely from Russell 2004
Appendix C #83's printed per-leg real-eph state, lands **all 7 Mars encounters inside
the 3-Mars-SOI confirmation band** on an INDEPENDENT REBOUND/IAS15 integrator over the
real DE440 ephemeris, at the published per-leg v∞ to 4 decimals (breathing 3.2–8.0
km/s, NOT the coplanar 5.10/3.05 anchors). The g (E→E) arcs stay far from Mars
(closest 0.67–1.05 AU, aphelion < 1.4 AU). S1L1 / `russell-ch4-4.991gG2` (#94) is
**independently confirmed on DE440**. This inverts the #165 DRIFT — for the corrected
topology, not the broken #164 one.

---

## What changed from #165 (why this confirms where #165 drifted)

| | #164/#165 (DRIFT) | #167 (CONFIRMED) |
|---|---|---|
| topology | BOTH arcs Mars-crossing | g = E→E free return (no Mars); only G = Mars transit |
| Mars target | radius + speed + ToF only | TRUE DE440 Mars position/longitude (App-C seed) |
| seed | two-body continuation FIT | Russell App-C #83 printed real-eph state (sourced) |
| n-body Mars miss | ~2.6 AU (both arcs) | 380–3,200 km (Sun-only) / 6,600–40,900 km (Mars-perturbed) |
| longitude | ~110° off real Mars | rendezvous to < 0.5° (Mars at lon 201.0° on 2027-06-13 for leg 2) |

The #165 DRIFT was structural: a `(radius, V∞, ToF)` closure with no longitude
constraint, on the wrong topology. The corrected build seeds from the real-eph state
that, by construction, places the spacecraft at real Mars's true position — so the
encounter is a genuine longitude rendezvous, not a radius coincidence.

---

## THE GATE — per-encounter Mars miss + v∞ (achieved vs App-C real-eph)

Construction: 21 leg arcs reconstructed from the App-C #83 block
(`src/cyclerfinder/search/s1l1_corrected.py`), Russell's recipe (v_sc = v_planet[DE440]
+ v∞ at each leg start). The 7 Mars-transit (G) outbound legs are 2, 5, 8, 11, 14, 17,
20 (arriving at Mars nodes 3, 6, 9, 12, 15, 18, 21).

**SOI band used:** miss < 3 × r_SOI(Mars) = 3 × 0.00386 AU ≈ **0.0116 AU ≈ 1.73e6 km**
— the SAME band #165 used, kept identical, never loosened.

### (1) Sun-only IAS15 — matches Russell's patched-conic cruise model

| Mars node | cycle | n-body miss (km) | miss (AU) | achieved v∞ | published v∞ | in-band? |
|---|---|---|---|---|---|---|
| M3  | 1 | 2,725 | 1.8e-05 | 5.248 | 5.248 | yes |
| M6  | 2 | 1,279 | 8.6e-06 | 7.693 | 7.693 | yes |
| M9  | 3 |   439 | 2.9e-06 | 4.657 | 4.657 | yes |
| M12 | 4 | 3,173 | 2.1e-05 | 3.198 | 3.198 | yes |
| M15 | 5 | 2,383 | 1.6e-05 | 6.262 | 6.263 | yes |
| M18 | 6 |   879 | 5.9e-06 | 8.046 | 8.046 | yes |
| M21 | 7 |   383 | 2.6e-06 | 3.219 | 3.219 | yes |

All 7 ≪ 1 Mars SOI; v∞ matches the published per-leg value to ≤ 1e-3 km/s. Mean Mars
v∞ 5.47 (App-C avg 5.475); span 3.198–8.046 — the real-eph "breathing", NOT a single
coplanar anchor. Energy drift 0.0; all converged. (REBOUND/IAS15 Sun-only is pinned to
the analytic two-body to < 1 km by GOLDEN GATE 1, so this independent integrator
reproduces the patched-conic seed exactly — same numbers as a pure Kepler propagation.)

### (2) Mars-perturbed IAS15 — real DE440 Mars as a continuous perturber

| Mars node | n-body miss (km) | miss (AU) | achieved v∞ | in-band (< 1.73e6 km)? |
|---|---|---|---|---|
| M3  | 14,145 | 9.5e-05 | 5.796 | yes (0.008 SOI) |
| M6  |  6,928 | 4.6e-05 | 8.459 | yes |
| M9  | 20,236 | 1.4e-04 | 5.091 | yes |
| M12 | 37,553 | 2.5e-04 | 3.537 | yes |
| M15 |  9,682 | 6.5e-05 | 6.933 | yes |
| M18 |  6,650 | 4.4e-05 | 8.810 | yes |
| M21 | 40,866 | 2.7e-04 | 3.529 | yes |

Turning on Mars's continuous gravity (the flyby bend Russell models as instantaneous)
shifts the miss to 6,600–40,900 km — still < 0.1 Mars SOI, comfortably in-band. The
v∞ rises a little (the perturber pulls the spacecraft in, deepening the encounter) but
stays in the real-eph 3.5–8.8 envelope.

### The g (E→E) free-return arcs stay FAR from Mars (the structural fix)

| g leg | cycle | closest approach to DE440 Mars (AU) | arc aphelion (AU) | Earth-return miss (km) |
|---|---|---|---|---|
| 1  | 1 | 1.048 | 1.266 | 11,051 |
| 4  | 2 | 0.869 | 1.283 | 11,209 |
| 7  | 3 | 0.666 | 1.173 |  7,391 |
| 10 | 4 | 0.899 | 1.241 | 10,308 |
| 13 | 5 | 0.967 | 1.329 | 13,477 |
| 16 | 6 | 0.774 | 1.253 | 10,355 |
| 19 | 7 | 0.789 | 1.207 |  9,250 |

Every g arc's closest approach to Mars is 0.67–1.05 AU — orders of magnitude outside
any encounter band — and aphelion stays sub-Mars (< 1.4 AU < Mars's 1.52 AU). Leg 1
matches #166's DE440 probe exactly (closest 1.05 AU, aphelion 1.27 AU). This is the
direct refutation of the #164/#165 both-arcs-Mars-crossing error: forcing the g arc to
cross Mars's radius is geometrically impossible for the real cycler. (The ~7,000–13,000
km Earth-return miss per g leg is the two-body-vs-DE440 deviation that the next leg's
fresh App-C v∞ node re-anchors — Russell's per-leg reproduction recipe.)

---

## On the patched-conic departure artifact (honesty note)

A third run — the #165-style handoff that starts the spacecraft a full Earth-SOI
radius (~925,000 km) DISPLACED along the departure asymptote and then perturbs by E+M
— gives larger misses (0.007–0.033 AU). This is **not** real geometry: the App-C seed
is a departure VELOCITY off Earth, not a state already clear of Earth's gravity well.
Displacing the start point and applying continuous Earth perturbation in the departure
region bends the trajectory; the miss shrinks **monotonically** toward the SOI band as
the start point moves clear of Earth (2,725 → 8.9e6 → 4.2e6 → 1.95e6 → 642,000 km for
0 / 0.5 / 1 / 2 / 5-day clean two-body hops on leg M3). It conflates the seed-handoff
convention with the geometry, so it is documented here but NOT used as the gate. The
faithful independent checks are the Sun-only (Russell's own cruise model) and
Mars-perturbed (real flyby gravity) runs above, both fully in-band.

---

## Golden / honesty

- **EXPECTED = Russell App-C #83 printed per-leg values** (transit time + v∞), verbatim
  in `src/cyclerfinder/search/s1l1_corrected.py:APPC_LEGS` / `APPC_MARS_TRANSIT`. The
  n-body miss / v∞ are EVIDENCE, never imposed. The seed traces entirely to Russell's
  printed state block; nothing in this build is fit.
- The gate targets the **real-eph** numbers, NOT the coplanar 5.10 (Russell Table 4.9)
  or 3.05 (Rogers/CPOM) Mars v∞ — those are idealisations the real cycler does not
  reproduce (#166 §4). The achieved v∞ matches the real-eph per-leg values, which
  breathe 3.2–8.0.
- The 3-SOI band is the SAME one #165 used; it was NOT loosened. A CONFIRMED that needed
  a loosened band would be invalid — here the misses are 50×–4,500× INSIDE the band on
  the Sun-only model and ≥ 40× inside on the Mars-perturbed model.
- **No catalogue writeback.** S1L1 / #94 V3 promotion is the main session's call.

---

## V3-writeback recommendation (CONFIRMED → recommend, for main-session review)

**Recommend the first V3 writeback for `russell-ch4-4.991gG2` (#94)** to the main
session, with this evidence chain:

- `_LEVEL_EVIDENCE` chain:
  1. **Sourced real-eph seed** — Russell 2004 Appendix C #83 per-leg
     `(epoch, time-start, v∞-vector)` block, transcribed verbatim
     (`s1l1_corrected.py:APPC_LEGS`); the published cycler's own reproduction data.
  2. **Topology** — corrected to E → g(E-E free return) → E → G(E-M-E transit) → E
     (#166 §3, Russell §4.8: lowercase = E-E, uppercase = Mars transit; one Mars
     encounter/cycle). The g arcs verified sub-Mars (aphelion < 1.4 AU).
  3. **Independent n-body confirmation (DE440)** — all 7 Mars encounters in the 3-SOI
     band at the published per-leg v∞ (4-decimal), on REBOUND/IAS15 sharing none of any
     solver's machinery; verified on both the Sun-only (patched-conic) and
     Mars-perturbed models. Test `tests/nbody/test_s1l1_corrected_nbody.py`.
  4. **Golden gates intact** — the integrator is pinned to analytic two-body < 1 km
     (GOLDEN GATE 1), DE440 read < 1e-6 km (GATE 3), Sun-only energy drift < 1e-10
     (GATE 2); the confirmation band is grounded in Mars SOI, not loosened.
- **Caveat to record on the row:** the real-eph Mars v∞ is epoch-dependent (3.2–8.0,
  avg 5.48); the catalogue's coplanar anchors 5.10/3.05 are idealisations the real
  cycler does not hold (#166 §4). The V3 evidence should cite the App-C real-eph
  values, not the coplanar anchors.

---

## Test results (verbatim)

```
$ uv run pytest tests/search/test_s1l1_corrected.py
5 passed in 3.57s

$ uv run pytest tests/nbody/test_s1l1_corrected_nbody.py -m "slow or not slow"
2 passed in 5.09s
```

Tests: `tests/search/test_s1l1_corrected.py` (5, construction mechanics, fast —
topology, arc classification, real-Mars intercept + true-longitude rendezvous,
per-leg v∞ breathing, g-arcs-far-from-Mars); `tests/nbody/test_s1l1_corrected_nbody.py`
(2 @slow — the independent CONFIRMED gate + g-arc clearance).

**No harness regression.** The #165 DRIFT test and the nbody golden gates stay green
(the #165 DRIFT of the WRONG-TOPOLOGY #164 construction is real and still valid; it is
NOT inverted — the corrected geometry is confirmed in the new sibling test):

```
$ uv run pytest tests/nbody/test_s1l1_nbody_closure.py tests/nbody/test_golden_*.py -m "slow or not slow"
14 passed in 5.80s
```

`ruff check` / `ruff format --check` / `mypy src tests` all clean.

**No catalogue writeback.** The emerged n-body geometry is evidence; the App-C
real-eph values are the EXPECTED target. S1L1 / #94 V3 promotion is the main session's
call — recommended CONFIRMED.
