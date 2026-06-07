# Jones / Hernandez / Jesick 2017 — METHOD deep-dive (B-plane, SNOPT, broad-search seeding)

Mined 2026-06-07 (Task #142). **Method-section deep-dive** complementing the
existing data note `docs/notes/2026-06-05-jones-aas17-577-vem-mining.md` (which
covered the member tables). This note extracts every implementable detail of the
broad-search seeding, the powered-flyby / B-plane targeting, and the SNOPT
ephemeris optimization — the Phase C shooter implements exactly this.

**Source (cite exactly, no file path):**
Jones, D. R., Hernandez, S., & Jesick, M., "Low Excess Speed Triple Cyclers of
Venus, Earth, and Mars," AAS 17-577, AAS/AIAA Astrodynamics Specialist
Conference, Stevenson WA, Aug 20-24 2017 (NTRS 20190028464). JPL/Caltech.

> Method pages (pp.5-8) render cleanly; every equation transcribed verbatim.

---

## 1. The algorithm in 3 lines

Construct **near-Hohmann seed legs** by Lambert over a (t0, Δt1) grid; grow each
seed into a full itinerary by Lambert broad-search over per-leg flight-time grids
(all revs × fast/slow), filtering on **v∞-continuity** (powered tangential flyby,
B-plane geometry) + **altitude window**; match single-cycle solutions across
opportunities for multi-cycle repeatability; **optimize selected solutions to
fully ballistic in the true ephemeris with a 2-step homotopy + control-point/
break-point model + SNOPT**. Mapping verdict: **MAPS DIRECTLY — this is the
blueprint our Phase C shooter implements** (broad-search seeding + B-plane
targeting + SNOPT). Our `correct.py` vector-residual mode is the in-residual
version of their post-hoc flyby-feasibility filter.

---

## 2. The method stated precisely

### 2.1 Broad-search model (p.5)
- "A **zero-sphere-of-influence patched conic** gravity model ... with the real
  planetary ephemeris, and Lambert's problem is solved to determine legs
  connecting consecutive encounters."
- **Flight time is the primary search variable.** Revolutions 1→max enumerated,
  plus **fast/slow (type 1 / type 2)** arcs. "Lambert's problem admits four
  solutions for a given number of revolutions, but here only **prograde**
  transfers are considered." Integer-π transfers excluded. Cites Russell [12] for
  Lambert.
- Lambert arcs yield incoming/outgoing asymptotes; flybys must have **altitude
  between 100 km and 100,000 km**.

### 2.2 Flyby evaluation — powered tangential maneuver + B-plane (pp.5-6)
After adjacent legs are Lambert-solved, **powered hyperbolic flybys** are computed
to (a) evaluate feasibility and (b) correct any v∞ discontinuity. "velocity
increments below 200 m/sec are permitted since experience has shown these can be
differentially corrected in high-fidelity dynamics to be entirely ballistic."

- **Transfer (turn) angle** (Eq.1): `δ = ∠(v∞⁻, v∞⁺)`.
- **Periapsis radius `r_p`** solved iteratively (Eq.2):
  `asin(μ/(μ + r_p v∞⁻²)) + asin(μ/(μ + r_p v∞⁺²)) = δ`.
  (Subsurface solutions allowed here, removed later by the min-altitude filter.)
- **Periapsis speeds** (Eq.3): `v_p⁻ = sqrt(v∞⁻² + 2μ/r_p)`,
  `v_p⁺ = sqrt(v∞⁺² + 2μ/r_p)`. The tangential maneuver Δv = `v_p⁺ − v_p⁻`.
- **B-plane unit vectors** (Eq.4, body-centered equatorial):
  `Ŝ = v∞⁻/v∞⁻`; `T̂ = (v∞⁻/v∞⁻ × k̂)/‖(v∞⁻/v∞⁻ × k̂)‖`; `R̂ = Ŝ × T̂`,
  with `k̂ = (0,0,1)` (pole).
- **B-plane angle** (Eq.5): `θ_B = atan2(v∞⁺/v∞⁺ · R̂, v∞⁺/v∞⁺ · T̂) − π`
  (atan2 range (π,−π)). "The flyby bends the excess velocity vector such that the
  projection of v∞⁺ onto the B-plane is along the −B vector."
- Two periapsis states formed; propagated fwd/bwd from periapsis to the SOI
  crossing. **Luidens [14]** gives the analytic SOI-propagation time.
> Note: the maneuver is **tangential** (sub-optimal but a fast filter); the
> guess "suffices for filtering poor solutions via constraint evaluation".

### 2.3 Broad-search algorithm (p.6) — the seeding procedure
For a family with repeat period `T = k·T_syn` (`t_f − t_0 = T`):
1. For the initialization year, compute the E-M (or M-E) **Hohmann transfer time
   `t_0*` and flight time `Δt_H`**.
2. Build **near-Hohmann seed legs** via Lambert over grids:
   - `t_0 ∈ {t_0^min, t_0^min+Δt_0, …, t_0^max}` with `t_0^min < t_0* < t_0^max`.
   - `Δt_1 ∈ {Δt_H^min, …, Δt_H^max}` where
     `Δt_H^min = Δt_H + t_0^min − t_0*`, `Δt_H^max = Δt_H + t_0^max − t_0*`.
3. For each seed leg:
   (a) feasible leg-2 options over a grid of `Δt_2`, all revs + fast/slow;
   (b) feasible final-leg options over a grid of `Δt_f` (`t_f` is known);
   (c) for six total flybys, third leg similarly;
   (d) for each feasible leg combination, the **last un-computed leg is assessed
       for feasibility** (only one flight time, but revs + fast/slow still
       enumerated);
   (e) **all fully feasible solutions are saved.**
- **Seed constraints** (p.7): seed-leg epochs/ToF deviate **≤ 50 days** from the
  Hohmann transfer; seed v∞ at Earth/Mars **below v∞^max = 5 km/s**.
- **Interior-flyby feasibility** (Eq.7): `‖v∞(t_i⁺) − v∞(t_i⁻)‖ < Δv∞^max` AND
  `100 km < r_p(t_i) − R_planet < 100,000 km`, with **Δv∞^max between 100 and
  200 m/s**.
- Final filter: ensure **no unintended intermediate flybys** along the trajectory.

### 2.4 Multi-cycle matching (p.7)
Because planetary alignment is not exactly repeatable, single-cycle solutions are
**combinatorially matched** across opportunities (e.g. a 2020 1-synodic set
matched with a 2026 set) rather than broad-searching over 2+ cycles (10-12
encounters). "The matching also permits the mixing of cycler families (e.g. an
EMEEVE followed by a EMEVVE)." Often "a simple shift in the middle (interior leg)
flight times within a given family is sufficient to maintain feasibility into the
next cycle."

### 2.5 True-ephemeris optimization (pp.7-8) — the SNOPT problem
- **Selection criterion**: high-quality approximate solutions = **low excess
  speed + low Δv**; flyby altitudes secondary.
- **Two-step continuation (homotopy)**: Step 1 = Sun + planets, iterated to
  continuity; the converged Step-1 result seeds Step 2 = adds **all planets +
  Earth's moon**.
- **Constraints enforced on the optimization**:
  1. flyby periapsis altitudes **100 km–100,000 km**;
  2. **continuity to 1.0E-3 km position and 1.0E-6 km/s velocity**.
- **Model**: control-point (CP) / break-point (BP) — integration forward+backward
  from each CP, continuity enforced at the BPs (between adjacent CPs). Initial CP
  state from the broad-search solution; **hyperbolic flyby orbits take precedence
  over the Sun-centered Lambert arcs**.
- **Optimizer: SNOPT** (SQP) [15]. "Much effort is taken to set bounds, scaling,
  and step-size control for the state and time parameters to ensure quality
  convergence."

### 2.6 Broad-search results constraints recap (p.8)
- **No 1-synodic (6.4 yr) feasible cyclers** of any family (only single-subsurface
  -flyby near-solutions). **Thousands of 2-synodic (12.8 yr) feasible** cyclers.
- EMEVVE / EMEEVE (six flybys, consecutive Earth/Venus) best overall.
- Matching built 2-cycle (25.6 yr) trajectories, avg transit-leg v∞ < 5 km/s.

---

## 3. Maps to our X / does not map

| Jones method element | Our code / Phase C | Verdict |
|---|---|---|
| Near-Hohmann seed legs via (t0, Δt1) Lambert grid, ≤50 d from Hohmann, seed v∞<5 km/s | Phase C broad-search seeding; `search/scan.py` | **MAPS — this IS the Phase C seeding recipe.** The 50-day Hohmann window + 5 km/s seed cap are directly implementable grid bounds. |
| Broad-search: flight-time-primary, all revs × fast/slow, prograde only, integer-π excluded | `core/lambert.py` (multi-rev, type1/2), `search/scan.py` | **MAPS — our Lambert already enumerates revs + fast/slow.** Confirm we exclude integer-π and restrict prograde in the scan. |
| Powered tangential flyby: δ (Eq.1), r_p iterate (Eq.2), v_p± (Eq.3) | `core/flyby.py` (bend/max_bend) | **PARTIAL — we have δ/turn geometry; the tangential periapsis-speed Δv (Eq.3) is the same physics as Russell's Eq.5.5 powered-SOI Δv.** We do not currently compute the tangential-maneuver Δv as a filter. |
| B-plane frame Ŝ,T̂,R̂ (Eq.4) + θ_B (Eq.5) | (absent) | **DOES NOT MAP — Phase C shooter needs this.** This is the exact B-plane targeting setup the task says Phase C implements; Eqs.4-5 are the implementable formulas (body-centered equatorial, k̂=pole). |
| Interior-flyby feasibility (Eq.7): ‖Δv∞‖<Δv∞^max(100-200 m/s) AND 100<alt<100,000 km | `correct.py` `residual_mode="vector"` (bend hinge) + `bend_feasible` | **MAPS — our vector residual mode is the in-solver version; the 100-200 m/s Δv∞ tolerance + altitude window are directly adoptable thresholds.** |
| SNOPT 2-step homotopy + CP/BP model, continuity 1e-3 km / 1e-6 km/s | our verify/propagate + (no SNOPT) | **PARTIAL — same homotopy/continuity philosophy; we lack a SNOPT driver and a CP/BP multi-shoot.** Russell's Ch.5 multiple-shooting (see russell note) is the more detailed version of the same CP/BP idea. |
| Multi-cycle combinatorial matching (single-cycle sets across opportunities; family mixing) | (absent) | **DOES NOT MAP — design idea.** Cheaper than broad-searching 2+ cycles; relevant if we pursue multi-cycle VEM cyclers. Our catalogue currently treats cyclers per-family, not matched architectures. |
| Hyperbolic flyby orbits take precedence over Sun-centered Lambert arcs (in CP/BP) | (implicit) | **Design note** — when stitching, the flyby state is authoritative, not the heliocentric Lambert leg. |

---

## 4. Candidate test anchors

The method section adds NO new tabulated trajectories beyond Tables 1-4 (already
captured in the data note). The **method constants** are the implementable
anchors (parameter values, not goldens):
- Seed window: **≤50 d** from Hohmann; seed **v∞^max = 5 km/s**.
- Flyby altitude window: **100 km – 100,000 km**.
- Interior-flyby Δv∞ tolerance: **Δv∞^max ∈ [100, 200] m/s**.
- Ephemeris continuity tolerance: **1.0E-3 km position, 1.0E-6 km/s velocity**.
- Homotopy: **Step 1 Sun+planets → Step 2 +all planets+Earth's moon**.
These pin the Phase C shooter's tolerances/bounds with a published provenance.

(Equation cross-check anchors: Eqs.1-5 are standard B-plane / powered-flyby
geometry; a unit test that reproduces θ_B / r_p for a known (v∞⁻, v∞⁺) pair would
have a source-traced EXPECTED only if Jones tabulated a worked example — he does
NOT, so these remain self-consistency checks, not goldens.)

---

## 5. Single most implementable finding (this paper)

**The B-plane targeting setup (Eqs.4-5) + the iterative powered-flyby r_p/Δv
(Eqs.1-3) as the Phase C flyby-feasibility kernel.** These five equations are the
complete, directly-codeable recipe for: given (v∞⁻, v∞⁺), compute the required
turn δ, solve for periapsis r_p, get the tangential Δv (Eq.3), and the B-plane
angle θ_B (Eq.5) for targeting. Combined with the Eq.7 feasibility gate
(Δv∞<100-200 m/s, 100<alt<100,000 km) this is exactly the broad-search filter the
Phase C shooter needs, with published tolerances. (Pair with Russell Eq.5.5 for
the powered-SOI Δv formulation — same physics, two independent sources.)

---

## 6. v4.2 backfill checks

- **center**: heliocentric Lambert legs; flybys body-centered equatorial (B-plane
  frame uses the planet pole k̂=(0,0,1)). No catalogue center ambiguity beyond
  what the data note already flagged.
- **tof_days_bounds**: the method gives ToF *grid* construction (Hohmann ±50 d
  for seeds; per-leg grids of Δt) but the realized per-leg ToFs are in Tables 2-4
  (data note). No new bounds.
- **source_ephemeris**: broad search = "**real planetary ephemeris**" (zero-SOI
  patched conic); optimization = full n-body (Sun+all planets+Earth's moon). The
  data note's `model_assumption: analytic-ephemeris` stands; the optimized cyclers
  are ballistic in **full n-body** (note this for any row promoted from Tables
  2/3 — the optimized version is n-body, not patched-conic). Specific DE version
  not stated in the method section.

---

## 7. Honest "not extractable" list

- No worked numeric example of Eqs.1-5 (no tabulated δ/r_p/θ_B for a specific
  flyby), so those equations yield self-consistency checks only, not goldens.
- The CP/BP control-point counts, SNOPT bounds/scaling values are not given
  numerically ("much effort is taken" — no numbers).
- The Hohmann `t_0*`, `Δt_H` per initialization year are not tabulated.
- Grid step sizes (Δt_0, Δt_2, …) for the broad search are not specified.
