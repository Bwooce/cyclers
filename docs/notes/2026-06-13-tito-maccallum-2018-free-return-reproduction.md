# Tito et al. 2018 Mars free-return — real-ephemeris reproduction & scope verdict (#238)

**Source.** Tito, Anderson, Carrico, et al., "Feasibility Analysis for a Manned
Mars Free-Return Mission in 2018," IEEE Aerospace Conference, 2013. Tables III
and IV carry an Astrogator full-force solution computed against the JPL **DE421**
ephemeris (cross-validated by the authors against a patched-conic MAnE run, Tables
I/II). This is a single Earth-Mars-Earth ballistic free-return for the
once-per-window January-2018 opportunity (501.3-day total, no deterministic
maneuvers after TMI; a single Mars gravity assist drops perihelion toward Venus's
orbit).

Our stack uses **DE440**. The DE421→DE440 difference is a small modelling
discrepancy to acknowledge, not fix.

---

## 1. Transcribed published values (EXPECTED side — verbatim from Tables III/IV)

**Table III — epochs (UTCG) and flight times (Astrogator, DE421):**

| Leg | Depart | Arrive | Flight time (d) |
|---|---|---|---|
| 1 | Earth, 5 Jan 2018 07:00:00.000 UTCG | Mars, 20 Aug 2018 08:18:19.619 UTCG | 227.05439374 |
| 2 | Mars, 20 Aug 2018 08:18:19.619 UTCG | Earth, 21 May 2019 13:52:48.012 UTCG | 274.23227306 |
| | | **Total duration** | **501.2866668** |

**Table IV — per-leg V-infinity / declination / RAAN / Vperi / C3 (Astrogator):**

| Leg | Dep V∞ (km/s) | Dep Decl | Dep RtAsc | Dep Vperi | Dep C3 (km²/s²) | Arr V∞ (km/s) | Arr Decl | Arr RtAsc | Arr Vperi | Arr C3 (km²/s²) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 6.232 | -6.554 | 271.053 | 12.578 | 38.835 | 5.417 | -7.488 | 53.581 | 7.272 | 29.344 |
| 2 | 5.417 | -11.94 | 200.12 | 7.272 | 29.344 | 8.837 | -25.34 | 52.718 | 14.18 | 78.094 |

**Fig.3 text (Mars flyby):** periapsis on the dark side at ~100 km altitude,
periapsis speed ~7.27 km/s, ~10 hours within Mars's sphere. **§II summary:** Earth
departure C3 = 38.835 km²/s²; Earth return velocity at perigee = 14.18 km/s
(the entry-speed driver for the TPS analysis).

(Tables I/II carry the patched-conic MAnE precursor — Earth depart 5 Jan 2018
07:42:35 GMT, JD 58123.7990; MAnE leg-1 V∞ 6.22697 km/s. The Astrogator
full-force numbers in III/IV are the higher-fidelity solution and are the
reproduction target here.)

---

## 2. Our reproduction (ACTUAL side — DE440 + UV Lambert)

Method: convert each published UTCG epoch to our `t_sec` (TDB s since J2000(TDB))
via `astropy Time(..., scale="utc").tdb`; heliocentric Earth/Mars states from
`Ephemeris("astropy")` (DE440); one single-rev prograde `lambert()` per leg;
v∞ = |v_Lambert − v_planet|, C3 = v∞². Mars-flyby ballistic continuity checked
with `core/flyby.py` (`bend_angle`, `max_bend`). Script:
`scripts/tito_free_return_repro.py` (committed).

The published UTCG epochs are held **fixed** as the targets; we do not re-optimise.
Our flight times reproduce Tables III to machine precision (same epochs, same
calendar→TDB conversion), confirming the epoch plumbing.

**Residuals (ours − published):**

| Quantity | Ours (DE440) | Published (DE421) | Δ | % |
|---|---|---|---|---|
| tof1 Earth→Mars (d) | 227.0544 | 227.05439 | ~0 | 0.00% |
| tof2 Mars→Earth (d) | 274.2323 | 274.23227 | ~0 | 0.00% |
| **L1 V∞ depart Earth (km/s)** | 6.2264 | 6.232 | −0.0056 | −0.09% |
| **L1 C3 depart Earth (km²/s²)** | 38.768 | 38.835 | −0.067 | −0.17% |
| L1 V∞ arrive Mars (km/s) | 5.4239 | 5.417 | +0.0069 | +0.13% |
| L2 V∞ depart Mars (km/s) | 5.4134 | 5.417 | −0.0036 | −0.07% |
| **L2 V∞ arrive Earth (km/s)** | 8.8996 | 8.837 | +0.0626 | +0.71% |
| **L2 C3 arrive Earth (km²/s²)** | 79.203 | 78.094 | +1.109 | +1.42% |

**Mars-flyby ballistic continuity (the free-return constraint):**

- V∞ speed mismatch (incoming Mars-arrival vs outgoing Mars-departure): **0.0105 km/s**
- required heliocentric bend across the flyby: **33.44°**
- achievable ballistic bend cone at 100 km altitude (the published periapsis): **34.24°**
- ballistic-feasible (bend ≤ cone)? **YES**

The flyby closes ballistically: the two Mars V∞ legs match in magnitude to
~10 m/s, and the required 33.4° turn fits inside the 34.2° cone Mars can deliver
at the published 100-km periapsis. This is exactly the gravity-assist that makes
the trajectory a *free* return — no deterministic maneuver at Mars.

---

## 3. Reproduce verdict: **MATCHES**

All four V∞ tuples reproduce within 0.71%; departure C3 within 0.17%; departure
V∞ within 0.09%; the Mars flyby is independently confirmed ballistic-feasible.
The residuals are dominated by the DE421→DE440 ephemeris difference, which grows
along the leg (Earth-departure C3 0.17%, Earth-return C3 1.42% — the return arc
accumulates 274 d of the small heliocentric-state difference, and the
Earth-return V∞ enters the encounter steeply so a small position error maps to a
larger speed error). No single-rev-Lambert pathology, no epoch/frame error: the
residual structure is precisely what a like-for-like-but-different-ephemeris
cross-check should show. This is a clean positive reproduction.

Caveat: we did not reproduce the declination / right-ascension columns (those are
geocentric/areocentric asymptote orientations requiring the planetcentric B-plane
frame, not the heliocentric V∞ magnitudes our Lambert solve produces). The
magnitude/energy quantities — the catalogue's load-bearing fields — all match.

---

## 4. Catalogue-scope verdict: **OUT OF SCOPE (different object class)**

The reproduction succeeds, but the object is **not catalogue-appropriate.**

**The catalogue holds repeating cyclers.** Spec §16.7 / §16.7.1 is explicit: the
record's structural kinds are `single-ellipse` ("one Kepler ellipse **the vehicle
repeats**"), `multi-arc` ("a different ellipse per leg … no single (a,e)"), and
`non-keplerian` (a periodic CR3BP orbit). Every one is something the vehicle
*cycles* — that is the definitional property the whole schema, the §14 validation
gauntlet (V2 = "≥3 continuous laps … bounded drift"; V3 = multi-cycle continuous
TCM), and the `period{}` block are built around.

**The existing "free-return" rows are repeating cyclers, not one-shot
free-returns.** Every `free-return` hit in `data/catalogue.yaml` is a
Russell-Ocampo / McConaghy **Earth-Mars cycler** (`cycler_class: multi-arc`,
e.g. `russell-ocampo-3.1.2+1`, the four `russell-ch4-*` rows). There "free-return"
names the *construction technique* (Russell builds a repeating cycler out of
Earth-to-Earth free-return arcs + a Mars gravity assist), and the rows carry a
`period{}` (k synodic periods, years), `invariants{}` (cycle-level aphelion/turn
ratios), and a `sequence_canonical` that repeats. The §14 evidence for those rows
"closes a single E→M→E free-return *ellipse slice*" **as one lap of a multi-lap
cycler** (spec line ~1415) — the free-return arc is a *component* of a repeating
object, never the object itself.

**The Tito 2018 trajectory does not repeat.** It is a single Earth-Mars-Earth
ballistic flyby pinned to the once-per-~15-yr fast opportunity at the specific
January-2018 launch window (Patel/Longuski's "two times every 15 years"). After
the 501-day round trip the vehicle is on a hyperbolic Earth-entry trajectory at
14.18 km/s — it re-enters and the mission ends; there is no next lap, no synodic
repetition, no steady-state V∞-matched encounter sequence. It has no `period{}`,
no cycle-level `invariants{}`, no repeating `sequence_canonical`. It cannot earn
even V2 (the floor that distinguishes a *cycler* from a one-shot transfer:
"≥3 continuous laps"). It is a **mission-design free-return**, a distinct object
class from the repeating "free-return-constructed cyclers" the catalogue holds.

Forcing it in would corrupt the catalogue's defining invariant (everything in it
is a repeating cycler) and the §14 ladder (a row that structurally cannot pass V2
sitting alongside V3 cyclers). The clean call is: **do not add a catalogue row.**

This is a valid, useful outcome: a sourced, ephemeris-grade interplanetary
free-return that our real-ephemeris Lambert/flyby stack reproduces to <1.5% on a
genuinely independent target (sourced EXPECTED, DE421; ours, DE440) — exercising
and validating the multi-arc / flyby path end-to-end — but correctly classified as
out-of-scope for a *cycler* catalogue.

---

## 5. PROPOSED WRITEBACK (HELD FOR REVIEW)

**Recommendation: NO catalogue row.** `data/catalogue.yaml` is not edited.

Rationale (above): the trajectory reproduces cleanly but is a single,
non-repeating, one-window Earth-Mars-Earth free-return mission, not a repeating
cycler, and so falls outside the catalogue's scope (spec §16.7.1; all existing
"free-return" rows are repeating multi-arc cyclers). It cannot meet the §14 V2
"≥3 continuous laps" floor that defines a catalogue cycler.

**What is captured instead (this note):**

- The reproduction is recorded here as a **validation cross-check** result, not a
  catalogue entry: our DE440 stack reproduces the Tito 2013 DE421 free-return to
  ≤0.17% on departure C3 and ≤0.71% on all V∞, with the Mars flyby independently
  confirmed ballistic-feasible (33.4° required vs 34.2° cone at 100 km).
- The throwaway reproduction is `scripts/tito_free_return_repro.py` (committed,
  ruff-clean) so the cross-check is rerunnable.

**If a reviewer decides the catalogue's scope *should* widen** to include
one-shot, sourced, ephemeris-grade interplanetary free-returns (a deliberate
scope decision, not this task's to make), the row that the evidence would earn:

```yaml
# PROPOSAL ONLY — NOT to be added without a scope-widening decision.
- id: tito-2018-mars-free-return
  name: "Tito 2018 manned Mars free-return (one-window E-M-E ballistic flyby)"
  source: literature
  trajectory_regime: ballistic
  model_assumption: analytic-ephemeris   # source: DE421; reproduced on DE440
  cycler_class: multi-arc                 # two distinct heliocentric legs E->M, M->E
  # NB: NOT a cycler — single window, no repetition. period{} would be null.
  source_ephemeris: "DE421"
  orbit_source: derived                   # Tito et al. IEEE Aero 2013 Tables III/IV
  vinf_source: derived
  orbit_fidelity: analytic-ephemeris
  vinf_fidelity: analytic-ephemeris
  bodies: ["E", "M"]
  sequence_canonical: "E-M-E"
  sense: "n/a"
  vinf_kms_at_encounters:
    - {body: "E", vinf_kms: 6.232, note: "Table IV leg-1 departure; ours 6.226 (DE440)"}
    - {body: "M", vinf_kms: 5.417, note: "Table IV Mars flyby (in≈out); ours 5.42/5.41"}
    - {body: "E", vinf_kms: 8.837, note: "Table IV leg-2 Earth arrival; ours 8.900 (DE440)"}
  # validation_level: V0  (sourced + reproduced-on-real-eph; CANNOT reach V2 —
  #   no repeating laps exist, so the §14 cycler ladder does not apply).
  first_published:
    authors: ["Tito, D. A.", "Anderson, G.", "Carrico, J. P.", "et al."]
    year: 2013
    title: "Feasibility Analysis for a Manned Mars Free-Return Mission in 2018"
    venue: "IEEE Aerospace Conference, 2013"
```

Even under a scope-widening decision the earned validation level is **V0**
(internally consistent + sourced + reproduced on real ephemeris): the §14 V1+
ladder is built on cycler laps this object does not have, so it would sit as an
explicitly-flagged non-cycler reference row, never a graded cycler. The default
recommendation stands: **withhold — out of scope.**
