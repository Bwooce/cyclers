# Beeson-Englander-Hughes-Schadegg 2015 — EMTG→GMAT medium-to-high-fidelity tool-chain (PROCESS mine)

Mined 2026-06-07. PROCESS/ALGORITHM mine aimed at ONE consumer: our n-body
validation ladder, specifically the GMAT V4 rung (#131/#136, Phase D of
`docs/superpowers/plans/2026-06-06-nbody-harness.md`). This paper is the
**published precedent** for exactly our planned medium→high-fidelity handoff:
a medium-fidelity global optimizer (EMTG, Sims-Flanagan + MBH + SNOPT) seeds an
auto-generated GMAT high-fidelity script, which re-converges in real dynamics.

**Source (cite exactly — authors/title/AAS number only):**
Beeson, R., Englander, J. A., Hughes, S. P., & Schadegg, M., "An Automatic Medium
to High Fidelity Low-Thrust Global Trajectory Tool-Chain; EMTG-GMAT," AAS 15-278,
AAS/AIAA Space Flight Mechanics Meeting, 2015. EMTG and GMAT are NASA Goddard
Space Flight Center (GSFC) tools; Beeson was then at U. Illinois Urbana-Champaign,
Schadegg at JPL.

> 14 pp., clean digital typeset, all body text + Algorithm 1 + Tables 1-2 read
> unambiguously. Figures are architecture sketches (Figs 1-3), GMAT GUI
> screenshots (Figs 4,5,7-11), and one EMTG scatter (Fig 6) — read for
> qualitative structure, not numbers (the load-bearing numbers are in Tables 1-2).

---

## 0. HEADLINE — the handoff recipe (the core prize), 5 lines

1. **EMTG solves at medium fidelity** (Sims-Flanagan transcription = forward/backward
   half-phases meeting at a match point, bounded impulses at step centers,
   **Kepler/patched-conic propagation**; MBH+SNOPT global search; integer-GA flyby
   sequence) → a converged **control decision vector** (p.2-3).
2. **EMTG auto-writes a GMAT `.script`** that *mimics the EMTG problem formulation*
   (same hierarchical mission/journey/phase/step structure) and **seeds it with a
   trajectory that is "nearly feasible in GMAT's high-fidelity model"** (p.3).
3. **What gets seeded:** the full control decision vector — control (thrust)
   per step, node states, epochs/times — transferred directly (EMTG medium→medium-high
   is a *direct return* of the control vector; the GMAT seed is the same vector
   written as a script) (p.3). **Boundary conditions are MOVED** from body-center
   (patched-conic) to **periapse of the hyperbolic flyby trajectory** (or to the
   parking orbit for departure/arrival), so the flyby body's gravity enters the EOM
   (p.4, Fig.3).
4. **GMAT re-converges in HIGH fidelity:** patched-conic dropped; near-body dynamics
   fully modeled; low-thrust = finite-burn segments; GMAT's NLP solver (VF13ad, an
   SQP) drives match-point defects to ~0 (p.4, p.10). Run **feasibility first**, then
   optionally **optimality** (p.5-6). Hard problems use **homotopy** (1st pass:
   match-point + thrust-magnitude only; 2nd pass: add auxiliary constraints) (p.8).
5. **"Validated" = match-point defects collapse** from patched-conic-seed magnitude
   to numerical zero while TOF/objective stay stable (Tables 1-2). No fixed ΔV/state
   *tolerance threshold* is published — acceptance is "feasible solution found"
   (defects driven below the NLP's feasibility tolerance), i.e. **convergence, not a
   numeric gap bound**. THIS IS THE KEY MAPPING CAVEAT FOR OUR V4 (see §3).

---

## 1. The handoff method (verbatim, with page numbers)

### Why two-step at all (p.1-2)
> "Solving the global optimization, low-thrust, multiple-flyby interplanetary
> trajectory problem with high-fidelity dynamical models requires an unreasonable
> amount of computational resources. A better approach … is a multi-step process
> whereby the solution of the aforementioned problem is solved at a medium-fidelity
> and this solution is used as an initial guess for a higher-fidelity solver." (p.1
> Abstract)

NMDB process: medium-fidelity solution in EMTG → high-fidelity in GMAT (p.2). The
distinguishing feature vs JPL's STOUR→MALTO/Mystic/Copernicus workflow is **the
level of automation**: "In its nominal run case, the proposed work-flow requires
very little human intervention." (p.2) — directly relevant to whether our GMAT
rung should be auto-script-generation (it should; see §2).

### What EMTG hands over (p.3, the seed contents)
> "Solutions with desirable characteristics can then be re-optimized using the
> medium-high fidelity mode in EMTG … The similarity in architecture between
> EMTG's medium and medium-high fidelity modes allows for a seamless transition
> between the two. **The control decision vector from the medium mode can be
> directly returned to EMTG to run in medium-high.** A GMAT script may be
> automatically generated from either EMTG mode. **The GMAT script provides GMAT
> with [a] fully defined low-thrust optimization problem which mimics the problem
> formulation in EMTG. The script also provides a very accurate initial guess with
> a trajectory that is nearly feasible in GMAT's high-fidelity model.** Finally, the
> new optimization problem is solved using one of GMAT's NLP solvers." (p.3)

So the seed is: **(a) the problem definition** (variables, constraints,
force-model, bodies — the *whole* NLP, templated) AND **(b) the initial-guess
trajectory** (control vector + node states + times). Both. Not just a state — the
entire control history and the problem structure.

### The fidelity step that re-converges in GMAT (p.4, Fig.3)
> "The reformulation of the low-thrust problem in GMAT takes on a similar
> formulation as the FBLT mode in EMTG, with the exception that **near body
> dynamics are fully modeled in GMAT.** … steps in GMAT are modeled using finite
> burns, and **the boundary conditions of two phases are moved from the center of
> the body to the periapse of the hyperbolic trajectory, which allows for the flyby
> body's gravity model to be fully included in the equations of motion.** Similarly,
> if a departure/arrival at the body occurs from/to a parking orbit, then the
> initial/final states are moved to that orbit." (p.4)

i.e. the medium→high jump is precisely **dropping the patched-conic approximation
at flybys** and re-converging the match-point defects in full n-body. This is the
SAME operation our shooter/rung does (Phase B/C): take a conic seed, re-converge
in real dynamics, measure the correction.

### Transcription detail (p.2-3, Fig.2) — for our shooter vocabulary
- Mission ⊃ journeys ⊃ phases ⊃ steps (Fig.1a). Phase = path between flybys.
- **MGALT** (multiple-gravity-assist low-thrust, Sims-Flanagan): forward half-phase
  propagates from initial state for ½ flight time, backward half from final state;
  they meet at a **match point** enforcing **state continuity**; thrust = bounded
  impulse at the temporal center of each step; spacecraft propagated by **solving
  Kepler's problem** (patched-conic) (p.3).
- **FBLT** (finite-burn low-thrust): full numerical integration forward/backward
  with piecewise-constant control (p.3). GMAT's reformulation resembles FBLT but
  with full near-body dynamics (p.4).
- This **match-point / forward-backward / continuity-defect** structure is exactly
  our multiple-shooting node/defect structure (plan §3, Tasks C.1-C.3). Their
  "match point" = our defect node; their "match-point defect" = our full-state
  continuity residual.

### GMAT script generation (p.4, Algorithm 1) — see §2.

---

## 2. Automation architecture (informs whether OUR GMAT rung should be script-generation)

**VERDICT: yes — auto-generated script is the published precedent, and it is what
our Phase D already plans.** Confirms `scripts/gmat_v4_aldrin.py` (a generator) is
the right shape.

Mechanism (p.4):
- GMAT is driven via its **custom scripting language** (not its GUI) for batch use.
- EMTG establishes a **GMAT API** ("a disjoint set of the EMTG code base") whose
  class objects **mirror EMTG's hierarchical event types** — `gmat_phase()` mimics an
  MGALT/FBLT phase and owns two `gmat_spacecraft` objects, which own `gmat_thrusters`
  and `gmat_fueltanks`. The API focuses on the prevalent GMAT optimization calls:
  **`Create Variable`, `Vary`, `Calculate`, `NonlinearConstraint`** (p.4).
- GMAT scripts have two sections: **(1) initialization** (data objects, variables,
  outputs) and **(2) the mission sequence** with the optimization sequence embedded
  (p.4).

**Algorithm 1 (GMAT Optimization Script Flow, p.5) — the template skeleton (verbatim
structure):**
```
Create Spacecraft, Tanks, Thrusters
Create Bodies, Propagators, Force Models, Coordinate Systems
Create Plots, Figures, Reports
Create Control Variables and Constants
Provide an Initial Guess          <-- the EMTG seed goes here
Start the Optimization Sequence
while Iterations Not Exceeded or Tolerances Not Satisfied do
  vary mission level parameters, calculate values
  for Journeys in Mission do
    vary journey parameters; calculate values
    for Phases in Journey do
      vary phase parameters
      for Steps in Phase do
        vary control parameters
        calculate and apply thrust magnitude constraint, apply constraint
        propagate spacecraft
      end for
      calculate and apply match-point defect constraints
      enforce continuity between phases
    end for
    calculate and apply journey level constraints
  end for
  apply mission level constraints
  calculate objective function
end while
```
- "Although the algorithm here shows the logic … using nested for-loops, in reality
  the for-loop structure is **completely expanded in script form. This can result in
  very large scripts. The Vesta-Ceres mission … results in a 4000 line script.**" (p.5)
- Templated = the structure (Algorithm 1) + body/force-model/solver boilerplate.
  Computed-and-injected = the **Initial Guess** (EMTG control vector + node states +
  times) and the per-step `Vary`/`Calculate`/`NonlinearConstraint` lines (loops
  fully unrolled).
- Tested against **GMAT R2013a, R2013b, R2014a** (p.10). "If requested by the user,
  EMTG can produce a GMAT script file automatically at the end of an MGALT or FBLT
  simulation. The script will automatically run in GMAT" (p.10).

**Mapping to OUR exporter (Phase D, plan Tasks D.0-D.1):**
- Our `generate_aldrin_script` should template Algorithm 1's two-section structure:
  init (Spacecraft/Tanks/Thrusters/Bodies/Propagators/ForceModel/CoordSys) +
  mission sequence with an `Optimize` block whose **`Provide an Initial Guess`** slot
  is filled with our medium-fidelity (conic/`maintain.py`/Aldrin) solution: node
  states, epochs, and (for the powered cycler) the impulse magnitudes/directions.
- Use the same four GMAT primitives the EMTG API leans on: `Create Variable`,
  `Vary`, `Calculate`, `NonlinearConstraint` (p.4). Our match-point/periodicity
  constraint = their "match-point defect constraints + enforce continuity" line.
- Move boundary conditions to flyby periapse / parking orbit (p.4) so Earth/Mars
  gravity is in the EOM — i.e. seed node states at the encounter, not at body-center.

---

## 3. Fidelity-gap data (golden-eligible cross-check for what "rung gap" magnitude is NORMAL)

These are the published medium(EMTG/patched-conic) → high(GMAT/full-dynamics)
deltas. **They calibrate our SILVER→GMAT expectations** (the SILVER ARTIFACT
verdicts came from exactly such a gap). NOTE the discipline boundary: these are
**EXTERNAL published numbers** about a DIFFERENT problem (low-thrust EMTG→GMAT),
so they are usable as an *order-of-magnitude calibration of what a healthy
patched-conic→n-body handoff looks like*, NOT as a golden EXPECTED for any cycler
row. Do not assert a cycler test against them.

### Table 1 (p.6) — Earth-Mars min-time transfer (impulse+low-thrust mix)
| Parameter | EMTG Initial Guess | GMAT Feasible | GMAT Optimized |
|---|---|---|---|
| Time of Flight (days) | 81.25 | 81.26 | 81.26 |
| ‖Δx‖ (km) [match-point position defect] | **3411.4** | **0.06** | 0.34 |
| ‖Δvx‖ (m/s) [match-point velocity defect] | 0.9 | 5E-7 | 7E-4 |
| Δm (g) [mass defect] | 0.2 | 8E-3 | 3E-2 |

**Reading:** the patched-conic seed carries a **~3400 km** match-point position
defect; full-fidelity GMAT collapses it to **~cm** (0.06 km feasible) while TOF
barely moves (81.25→81.26 d). The "gap" between medium and high fidelity is
~3400 km of position discontinuity that the high-fidelity solver *absorbs*, with a
**negligible TOF change**. This is a *small, healthy* handoff (ROBUST analogue).

### Table 2 (p.8) — Vesta-Ceres (DAWN reproduction, MGALT seed → GMAT)
| Parameter | EMTG Initial Guess | 1st Feasible Pass | 2nd Feasible Pass |
|---|---|---|---|
| TOF (days) | 2921.7 | 2962.8 | 2962.8 |
| Final mass (kg) | 943 | 932.2 | 932.2 |
| sup‖Δxᵢ‖ (km) [worst match-point pos defect] | **12E6 (1.2×10⁷)** | **1E-4** | 3E-2 |
| sup‖Δvxᵢ‖ (km/s) | 0.64 | 2E-10 | 3E-9 |
| sup‖Δmᵢ‖ (kg) | 17.5 | 7E-11 | 3E-10 |

**Reading:** the harder, multi-flyby DAWN-class seed carries a **12 million km**
worst-case match-point position defect; GMAT homotopy drives it to **~10 cm**
(1E-4 km), at the cost of **~41 days TOF growth** (2921.7→2962.8) and
**~11 kg final-mass loss** (943→932.2). The handoff *worked* but VF13ad nearly
failed: "after nearly 200 major iterations, ~18 hours … the simulate failed to
converge" without homotopy (p.8); with homotopy, 1st pass = 4 major iterations /
~40 min (p.8). Body force model: **Earth, Mars, Vesta, Ceres, Sun, and Jupiter**
(p.8) — note Jupiter included (relevant to our Gate-4 Jupiter-sensitivity arm).

### CALIBRATION TAKEAWAYS for our rungs / V4 (the load-bearing cross-check)
1. **A patched-conic seed's pre-correction position discontinuity of 10³-10⁷ km is
   NORMAL and recoverable** in a healthy handoff. Our rung's "terminal closure
   error" before correction being large (km-to-Mkm) is NOT by itself an ARTIFACT
   signal — it is *expected* of a conic seed. The ARTIFACT signal is the
   **correction COST** (their TOF/mass *change*; our node-impulse ΔV), not the raw
   pre-correction defect.
2. **A healthy handoff costs little in the mission metric:** TOF moved +0.01% (Earth-Mars)
   to +1.4% (DAWN), mass −1.2% (DAWN). Their acceptance is *defects→numerical-zero
   with the metric only mildly perturbed*. Translate to our ΔV thresholds: a small
   correction ΔV (plan ROBUST <200 m/s) ↔ their "mild TOF/mass change"; a large one
   ↔ a seed that "lives only in patched-conic land."
3. **No published fixed numeric TOLERANCE.** Their bar is *convergence of the NLP*
   (defects below feasibility tol), not a percentage. This CONFIRMS our plan's
   choice to **declare the V4 5% band ourselves up front** (plan Q5) — the precedent
   gives no external 5% number to borrow, so our self-declared band is the honest
   move, and 2.9138 km/s ±5% is OURS to own (not back-fit, not from this paper).
4. **Run feasibility first, optimality second; use homotopy on hard problems**
   (p.5-6, p.8). Our V4 run-book (Task D.2) should adopt the same two-stage pattern:
   GMAT feasibility pass first (drive match-point/periodicity defects to ~0), then
   an optional optimality pass. Homotopy (relax→add constraints across passes) is
   the published escape when the dense SQP stalls.

---

## 4. MBH specifics (does this paper source `search/mbh.py`'s perturbation default?)

**Partial — it NAMES the components but does NOT give the perturbation distribution
or its tuning.** `search/mbh.py` currently documents its `"cauchy"` default as an
unsourced stand-in pending Englander & Englander 2014 (ISSFD24 S7-3). This paper
does NOT close that gap; it cites the architecture, not the tuning.

What it states (p.2):
> "The EMTG-GMAT tool-chain begins the global optimization search at a
> medium-fidelity using EMTG's stochastic search that combines the Sims-Flanagan
> transcription with **a special variant of monotonic basin-hopping (MBH)** and the
> nonlinear programming solver, Sparse Nonlinear OPTimizer (**SNOPT**) … the
> capability to automatically select the optimal flyby sequence using an **integer
> genetic algorithm (GA)**." (p.2)

- Confirms the **MBH (special variant) + SNOPT + integer-GA** stack and the
  **Sims-Flanagan** transcription — the same architecture our survey Thread 1 and
  `mbh.py` docstring already cite.
- Calls the MBH a "special variant" but **does not specify the perturbation
  distribution (Cauchy/Pareto), the scale schedule, or per-gene scaling** — those
  remain in Englander & Englander 2014 (still on the #116 acquisition list).
- **References worth pulling for the MBH tuning** (this paper's bibliography):
  - Ref [8] Ellison, Englander, Conway, "Robust Global Optimization of Low-Thrust,
    Multiple-Flyby Trajectories," AAS/AIAA 2013 (Hilton Head). ← likely the
    MBH-robustness detail.
  - Ref [11] Englander, Ellison, Conway, "Global Optimization of Low-Thrust,
    Multiple-Flyby Trajectories at Medium and Medium-High Fidelity," AAS/AIAA Space
    Flight Mechanics 2014 (Santa Fe). ← the medium/medium-high split this paper builds on.
  - Ref [9] Englander, Conway, Williams, "Automated Mission Planning via Evolutionary
    Algorithms," JGCD 35(6):1878-1887.
  - Ref [19] Ellison, Englander, Conway, "Analytical Partial Derivative Calculation
    of the Sims-Flanagan Transcription Match Point Constraints," AAS/AIAA 2014. ←
    the analytic Jacobian we currently finite-difference (plan Q1, Risk 2 — our
    Tudat trigger). NOT the MBH tuning, but directly relevant to our STM-noise risk.

**Action for `mbh.py`:** no change to the SPEC CAVEAT — this paper confirms the
*architecture* citation but the `"cauchy"` tuning remains unsourced pending
Englander & Englander 2014. Optionally add Ellison/Englander/Conway 2013 (ref [8])
to the acquisition list as a secondary MBH source.

---

## 5. v4.2 backfill check (confirm low-thrust examples are NOT cycler rows)

**CONFIRMED: nothing here is a catalogue/cycler row.** The two worked examples are
(a) an Earth-Mars minimum-time low-thrust transfer (impulsive departure + LT cruise
+ impulsive Mars rendezvous; unrealistic Isp, illustrative only, p.5) and (b) a
DAWN reproduction (Earth→Mars-flyby→Vesta→Ceres low-thrust rendezvous, p.8). Both
are **one-way low-thrust rendezvous/flyby missions, NOT periodic Earth-Mars
cyclers**. No periodicity, no repeat synodic structure. No v4.2 row backfill applies.

### Ephemeris / force-model provenance flags used in their GMAT runs (record for our V4 run-book)
- **Ephemeris:** SPICE kernels (EMTG ref [12] "SPICE Ephemeris," naif.jpl.nasa.gov;
  p.3 "detailed hardware modeling and ephemeris available in EMTG … including SPICE
  kernels"). GMAT uses its own DE ephemeris. The paper does NOT pin a specific DE
  release (no "DE440/DE430" stated). → For OUR V4 stamp, GMAT's independent DE +
  our shared-DE440 is the intended *independence* (plan §4), consistent with theirs
  (EMTG-SPICE vs GMAT-DE = independent ephemerides). FLAG: paper gives no DE number.
- **Force model (Vesta-Ceres, p.8):** Earth, Mars, Vesta, Ceres, **Sun, and
  Jupiter** point masses. Thrust modeled as constant 0.1 N at Isp 3100 s, 90% duty
  cycle (GMAT lacked accurate power/thrust models at the time, p.8). → Jupiter IS in
  their high-fidelity set; supports keeping Jupiter in our rung body set
  (Sun+E+M+J) and recording the Jupiter on/off sensitivity (plan Gate-4 arm).
- **Solvers:** EMTG = SNOPT (medium); GMAT = VF13ad (dense SQP, ref [17] HSL) or
  Fmincon (ref [16] MATLAB, untested). VF13ad limitations are the paper's main
  caveat: no bound constraints, dense-only, struggles >100 vars, finite-differenced
  partials only (p.4, p.13). → For our V4 we are reproducing one Aldrin solution
  (small problem), so VF13ad-class is adequate; the homotopy/feasibility-first
  pattern (§3.4) is the mitigation if it stalls.

---

## 6. Direct mapping onto OUR n-body harness plan (#131/#136)

| Their step (AAS 15-278) | Our equivalent (plan `2026-06-06-nbody-harness.md`) |
|---|---|
| EMTG medium (Sims-Flanagan + MBH + SNOPT) → control vector | Our conic corrector / `maintain.py` / free-return → `BallisticClosureResult` + node V∞ (`_vinf_nodes`) |
| Auto-write GMAT script mimicking EMTG formulation + seed | `scripts/gmat_v4_aldrin.py` (Task D.0) — confirmed: generator, not hand-script |
| Move BCs to flyby periapse / parking orbit (drop patched-conic) | Seed node states at encounters; GMAT full near-body EOM (Phase D force model) |
| GMAT feasibility pass → optimality pass (+ homotopy if hard) | V4 run-book two-stage procedure (Task D.2) — ADOPT this pattern |
| Acceptance = match-point defects → numerical zero, metric stable | Our V4 = GMAT reproduces 2.9138 km/s ±5% (Q5, self-declared — precedent gives NO external numeric tol) |
| Match-point / fwd-bwd / continuity-defect transcription | Our multiple-shooting node/defect structure (Tasks C.1-C.3) — same skeleton |
| Tables 1-2 medium→high deltas | Calibration of "normal rung gap": raw defect 10³-10⁷ km is normal; COST (ΔV/TOF change) is the ARTIFACT signal (§3) |

**Single most actionable finding:** the ARTIFACT/ROBUST distinction should key on the
**correction COST**, not the raw pre-correction discontinuity — their healthy
handoffs start from 3400 km / 12 Mkm defects and still converge cheaply. Our rung
already does this (it thresholds node-impulse ΔV, not terminal closure), so this
paper *validates the plan's threshold choice* and warns against re-keying on raw
closure error.

---

## Provenance / honesty notes

- All quotes carry page numbers from the AAS 15-278 PDF; Tables 1-2 transcribed
  verbatim.
- Tables 1-2 are EXTERNAL low-thrust numbers about a DIFFERENT problem class; they
  are calibration/order-of-magnitude context only, **never a golden EXPECTED** for a
  cycler row (golden discipline: EXPECTED must be sourced AND about the same thing).
- 2.9138 km/s remains OUR computed Aldrin value; this paper supplies NO Aldrin
  number and NO external tolerance — it *confirms* that self-declaring the V4 band
  (plan Q5) is the correct, honest move.
- MBH perturbation tuning is still UNSOURCED for `search/mbh.py`; this paper names
  the architecture, not the distribution. Englander & Englander 2014 (ISSFD24 S7-3)
  remains the acquisition target; add Ellison/Englander/Conway 2013 (AAS Hilton
  Head) as a secondary.
