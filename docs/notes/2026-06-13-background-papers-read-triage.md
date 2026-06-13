# Read-based relevance triage of background papers (#235)

**Status: READ-BASED TRIAGE.** A prior corpus audit classified these papers
from filenames only. The user correctly objected that title-based dismissal is
wrong — a reusable method (reachability formulation, cost surrogate,
primer-vector identity, missed-thrust margin, sequencing trick, sourced IC) can
hide in a paper whose title sounds off-topic. This note records verdicts after
**actually opening and reading** the abstract, method, and results/conclusion of
each of the 13 never-opened papers.

Relevance lens (what makes a paper ADOPT-CANDIDATE): does it offer something
usable for cycler discovery/search (family continuation, reachable-set/
accessibility, novelty), trajectory closure/correction (STM, shooting, primer
vector, missed-thrust margins), the catalogue (new sourced cycler ICs/rows), or
our validation gauntlet (independent cross-checks, cost/feasibility bounds)? The
project is mostly ballistic/impulsive heliocentric + CR3BP cyclers.

Outcome: **13 read — 2 ADOPT-CANDIDATES, 4 METHOD-NOTED, 7 OUT-OF-SCOPE.**
(Guzman 2002 was already mined separately in
`docs/notes/2026-06-07-guzman-2002-primer-survey-mining.md`; it is logged below
as METHOD-NOTED / already-digested, not re-derived.)

---

## ADOPT-CANDIDATES (turn these into todos)

1. **Zhou, Armellin, Qiao & Li (2025), "Single-Impulse Reachable Set in
   Arbitrary Dynamics Using Polynomials," arXiv:2502.11280.**
   - WHAT: a single-impulse reachable-set (RS) method that works in *arbitrary*
     dynamics (incl. CR3BP) via differential algebra (DA) + automatic domain
     split (ADS); the RS boundary is the envelope of a two-variable Taylor
     polynomial family, reduced to one-variable root-finding plus a local
     polynomial approximation. Validated on two Earth-Moon NRHOs in cislunar
     space (relative error < 0.0658 %, envelope-solve time cut > 84 %).
   - WHY ADOPT: this is a *same-model* (CR3BP) accessibility tool — it answers
     "from this state, with bounded Δv, what is reachable?", which is exactly the
     family-accessibility question. Directly comparable to the Braik-Ross
     accessibility line (#230). DACE/DACEyPy and the `alphashape` Python package
     are named and open-source.
   - PROPOSED FOLLOW-ON TASK: spike a single-impulse CR3BP reachable-set
     computation against our existing Earth-Moon CR3BP propagator (cross-check
     vs the paper's two published NRHO initial states, Tables 2-3, 15-digit) and
     compare the accessibility framing to Braik-Ross #230 before deciding whether
     to adopt DA or roll a cheaper STM-based RS.

2. **Tito, Anderson, Carrico et al. (2013), "Feasibility Analysis for a Manned
   Mars Free-Return Mission in 2018," IEEE Aerospace.**
   - WHAT: carries a fully sourced Earth-Mars-Earth *ballistic free-return*
     trajectory: depart Earth 5 Jan 2018, Mars flyby ~100 km altitude 20 Aug
     2018, return Earth 21 May 2019 (501.5-day total, no deterministic maneuvers
     after TMI), Mars gravity assist drops perihelion near Venus's orbit. Tables
     III/IV give exact UTCG epochs + Julian dates and per-leg V-infinity,
     declination, RAAN, and C3 from a high-fidelity numerical integration with
     the JPL DE421 ephemeris (cross-validated against patched-conic MAnE).
   - WHY ADOPT: a sourced, ephemeris-grade, impulsive interplanetary free-return
     with a flyby — the right *kind* of object for our catalogue/validation, and
     a clean independent cross-check target (sourced expected values, not ones we
     computed). It is single-cycle (not a steady-state cycler) but the closure
     constraints and the published V-infinity/C3 tuples are directly usable.
   - PROPOSED FOLLOW-ON TASK: attempt to reproduce the Tito 2018 free-return as a
     real-ephemeris cross-check (target the published Tables III/IV epochs and
     V-infinity/C3 values) using our ephemeris + Lambert/flyby stack; if it
     reproduces, log it as a sourced V-validation row; if not, it's a clean
     negative that exercises the multi-arc / flyby path (cf. S1L1 multi-arc
     finding).

---

## Per-paper verdicts

### 1. Silvestrini & Lavagna (2022), "Deep Learning and ANN for Spacecraft
Dynamics, Navigation and Control," *Drones* 6:270.
A 39-page review survey of deep-learning / ANN methods for spacecraft GNC
(system identification, optical navigation, control synthesis, hybrid ANN +
classical). Closed-loop on-board GNC scope; no trajectory-design, cycler,
reachability, or sourced-IC content.
**OUT-OF-SCOPE** — general ML-GNC literature survey; nothing for ballistic/
impulsive cycler discovery, closure, or the catalogue.

### 2. Li, Topputo & Baoyin (2019), "Autonomous Time-Optimal Many-Revolution
Orbit Raising for EP GEO Satellites via Neural Networks," arXiv:1909.08768.
Trains a supervised NN (on indirect-method / homotopy-generated trajectories) to
fly time-optimal many-revolution low-thrust LEO-to-GEO orbit raising
autonomously. Earth-orbit, continuous-thrust, closed-loop NN control.
**OUT-OF-SCOPE** — low-thrust GEO orbit-raising + NN guidance; no heliocentric,
impulsive, or cycler content.

### 3. Singh & Junkins (2022), "Stochastic learning and extremal-field map based
autonomous guidance of low-thrust spacecraft," *Sci. Rep.* 12:17774.
Gaussian Process Regression trained on an extremal-field bundle (neighboring
optimal trajectories via "perturbed back-propagation") to predict costates /
primer vector for on-board re-planning of a low-thrust Earth-3671-Dionysus
transfer under off-nominal thrust. The extremal-field-bundle idea and predicting
the primer vector are interesting, but it is low-thrust closed-loop guidance, not
a discovery/closure tool, and the predicted primer values are its own output (not
a sourced cross-check).
**OUT-OF-SCOPE** (marginal method-note on extremal-field bundles).

### 4. Hu, Yang, Li & Baoyin (2024), "Robust Design of Low-Thrust Gravity-Assist
Trajectories via Reinforcement Learning," AIAA J. (G009427).
RL-based multi-phase robust low-thrust gravity-assist trajectory design under
uncertainties (Earth-Earth-Jupiter), MDP per phase, with an analytic
reachability constraint folded into the interior-transfer reward to improve
GA-target reachability. Low-thrust + RL + uncertainty, closed-loop guidance.
**OUT-OF-SCOPE** — the analytic GA-reachability reward term is the only thing of
passing interest; the deliverable is an RL guidance policy, not usable for
ballistic cycler discovery/closure.

### 5. Sinha & Beeson (2025), "Initial Guess Generation for Low-Thrust Trajectory
Design with Robustness to Missed-Thrust-Events," arXiv:2501.06694.
Compares non-conditional vs conditional global-search initial-guess strategies
for robust low-thrust trajectory design, validated on the cislunar Gateway PPE
(min-fuel SEP transfer to the 9:2 NRHO). Entirely low-thrust NLP / global-search
machinery.
**OUT-OF-SCOPE** — the "conditional initial guess seeded from a less-robust
solution" is a generic continuation-style idea but nothing concretely adoptable
for our impulsive cyclers.

### 6. Venigalla, Englander & Scheeres (2020), "Low-Thrust Trajectory
Optimization for Maximum Missed Thrust Recovery Margin," AAS 20-438.
Defines the missed-thrust recovery margin (MTM) — the longest forced coast a
low-thrust *nominal* can absorb at its weakest point and still reach a terminal
manifold once thrusting resumes — and directly optimizes it via a "virtual
swarm" (nominal + many recovery trajectories co-optimized), with a Pareto front
over MTM / delivered mass / arrival date. Sims-Flanagan transcription, SEP duty
cycle: intrinsically low-thrust.
**METHOD-NOTED** — the MTM concept ("longest coast still able to reach the
target") is a clean robustness metric we could conceptually borrow if we ever
cost cycler-maintenance robustness, but the virtual-swarm machinery is
low-thrust-specific. No task.

### 7. Viavattene & Ceriotti (2021), "Artificial Neural Networks for Multiple NEA
Rendezvous Missions with Continuous Thrust," *JSR* 59(2):574-586.
ANN surrogate estimates low-thrust transfer time/cost between near-Earth
asteroids; the surrogate feeds a tree-search sequencer that finds feasible
multi-NEA rendezvous sequences, with each leg verified by an OCP solve.
**METHOD-NOTED** — the decomposition (combinatorial sequence search over a fast
*cost surrogate*, then per-leg verification) is structurally analogous to
sequencing cycler legs/encounters; logged as a pattern. Content itself
(low-thrust NEA tours + ANN cost net) is out-of-scope. No task.

### 8. Zhang, Acciarini, Izzo, Baoyin & Topputo (2026), "Pretrained Approximators
for Low-Thrust Trajectory Cost and Reachability," arXiv:2605.26790.
ML "Value Networks" that predict low-thrust optimal fuel and minimum transfer
time (= reachability, framed as the time-optimal lower bound) across orbital
regimes via a self-similar transformation, trained on a homotopy-ray-generated
dataset; reports a log-log scaling law (no saturation up to ~1e8 samples).
Low-thrust-specific, and the deliverable is a pretrained model gated by a very
large training-data floor (consistent with the prior ML-surrogate trio triage,
where the ozaki blueprint was deferred for the same reason).
**METHOD-NOTED** — the "reachability = time-optimal transfer-time boundary"
framing is a tidy reachability definition; the homotopy-ray sampling (deform
guaranteed-feasible solutions toward the reachability boundary) is a notable
data-generation trick. Both logged; no task (low-thrust + training-floor gated).

### 9. Zhou, Armellin, Qiao & Li (2025), single-impulse reachable set,
arXiv:2502.11280.
See ADOPT-CANDIDATES #1. **ADOPT-CANDIDATE** (DA-based single-impulse reachable
set in CR3BP; same-model accessibility tool; publishes two NRHO ICs).

### 10. Blender & Singh (2025), "Uncertainty-Aware Guidance for Continuous Thrust
Transfers Using Gradient-Boosted Decision Trees," AAS 25-524.
Builds an uncertainty-aware dataset by propagating sigma-point ensembles around a
low-thrust nominal bundle (Earth-3671-Dionysus), then trains LightGBM GBDTs to
map current belief state -> optimal control-to-go corrections for real-time
guidance. Low-thrust closed-loop guidance under state uncertainty.
**OUT-OF-SCOPE** — no cycler / impulsive / discovery content.

### 11. Putnam, Braun, Rohrschneider & Dec (2005), "Entry System Options for
Human Return from the Moon and Mars," AIAA 2005-5915.
Earth-entry (EDL) trade study: ballistic / lifting-capsule / biconic / lifting-
body configurations, direct entry vs aerocapture, thermal / deceleration /
crossrange / entry-corridor analysis; concludes L/D=0.3 suffices for both lunar
and Mars return. Pure atmospheric entry-vehicle design.
**OUT-OF-SCOPE** — EDL / aerocapture vehicle design; no orbital/trajectory/
cycler content.

### 12. Tito, Anderson, Carrico et al. (2013), manned Mars free-return 2018, IEEE
Aerospace.
See ADOPT-CANDIDATES #2. **ADOPT-CANDIDATE** (sourced ephemeris-grade Earth-Mars
ballistic free-return with Mars flyby; published epochs + V-infinity/C3 tuples).

### 13. Guzman, Mailhe, Schiff, Hughes & Folta (2002), "Primer Vector
Optimization: Survey of Theory, New Analysis and Applications," IAC-02-A.6.09.
Primer-vector survey (Lawden COV foundations; STM/adjoint formulation; necessary
conditions; interior/endpoint impulse criteria; the variational adjoint
equation; and — critically — *when the theory fails*: linear-theory breakdown for
eccentric/long arcs and isolated Phi_rv singularities at multiples of the orbit
period, with Stern's true-anomaly singularity condition).
**METHOD-NOTED / already-digested** — this paper was already mined in
`docs/notes/2026-06-07-guzman-2002-primer-survey-mining.md`, which used it to
lift the "PROVISIONAL pending Guzman 2002" caveats on the primer diagnostic
(#144) and primer-refine (#148) lanes. Re-reading confirms that note's use is
faithful (the failure conditions and the multi-rev singularity caveats are the
load-bearing content). No new task.

### 14. Fan, Cheng, Li, Pan, Huo & Qi (2025), "An Initial Trajectory Design for
the Multi-Target Exploration of the Electric Sail," *Aerospace* 12:196 (#234).
Bezier shape-based initial-trajectory design for an electric-sail (E-sail,
infinite-Isp continuous thrust) flying planetary flybys + gravity assists and
then escaping to the solar-system boundary. Heliocentric cylindrical coordinates
are each fit by an 8th-order Bezier curve whose endpoint coefficients are pinned
by boundary conditions; the free interior coefficients, flight times, GA
pericenter altitudes and pre-flyby velocities are optimized by `fmincon`
(interior-point) to minimize total flight time, with the E-sail thrust
coefficient kappa, pitch-angle feasibility and a 70 km minimum GA altitude as
constraints. Three worked cases (Mars-GA 2.62 yr; Jupiter-GA 11.15 yr;
Mars-Jupiter-GA 5.88 yr, the multi-GA shortening flight time 47.3% vs Jupiter
alone), each solved in <100 s. Despite the "multi-target" title, the planet
SEQUENCE is fixed by the analyst per example — there is no combinatorial
sequencing/ordering search (the appealing hook from the relevance lens); the
contribution is the continuous-thrust shape-based + GA-constraint plumbing.
**OUT-OF-SCOPE** — E-sail-specific continuous-thrust propulsion and shape-based
(Bezier) trajectory design, neither of which transfers to our ballistic/
impulsive heliocentric + CR3BP cyclers; and the hoped-for reusable multi-target
sequencing method is absent (sequences are hand-chosen, not searched). The
gravity-assist closure model (eq. 4: v-infinity matching, bend angle delta,
pericenter-altitude floor) is textbook patched-conic, already covered by our
flyby stack. Note the bibliography is a useful index of the shape-based +
multi-GA low-thrust literature (Petropoulos exp-sinusoid [29]; Fan's own MCTS
multi-GA sequencer [22]; Quarta-Mengali E-sail outer-system [10]) if that lane is
ever revisited, but nothing here is adoptable.
