# Triage skim of six marginal acquisitions

**Date:** 2026-06-07
**Method:** TRIAGE SKIM only — abstract, introduction, and method-section
headings read (first ~6 pages + a peek at any algorithm/results section). NOT a
full mine. Relevance judged against the live project frontiers (multi-arc
topology search, basin-selection / MBH, n-body shooting, maintenance-ΔV
optimality / primer vector, reachable-set tooling for TCM windows) per
`docs/notes/2026-06-07-external-algorithms-survey.md`.

Our project is **ballistic** (impulsive, conic + n-body) Earth-Mars cycler
discovery/validation; the only catalogued low-thrust lane is Pony-Express.
Low-thrust-only methods are therefore OUT-OF-SCOPE *as engines*, regardless of
quality; they are only retained if a sub-component (a search architecture, a
diagnostic, a reachability tool) transfers to our impulsive/maintenance stack.

---

## 1. Venigalla, Englander & Scheeres (2020) — Missed-thrust recovery margin

**Citation:** C. Venigalla, J. A. Englander, D. J. Scheeres, "Low-Thrust
Trajectory Optimization for Maximum Missed Thrust Recovery Margin," AAS 20-438,
AAS/AIAA Astrodynamics Specialist Conference, 2020.

**Summary (3 lines):** Defines missed-thrust recovery margin (MTM) as the
longest forced coast at a point that still lets the spacecraft reach its terminal
manifold. Introduces a "virtual swarm" technique that co-optimizes the nominal
trajectory with N recovery trajectories spawned at points along it, maximizing
the worst-case margin (or constraining it while maximizing delivered mass).
Transcription is multi-phase Sims-Flanagan; iteratively adds virtual spacecraft
at the worst margin-violation points; also builds a Pareto front (MTM vs mass vs
arrival date).

**Verdict: OUT-OF-SCOPE.** The whole construct (forced-coast recovery margin,
delivered-mass objective, Sims-Flanagan low-thrust transcription) is intrinsic
to low-thrust missed-thrust robustness and has no impulsive-cycler analogue; our
maintenance ΔV is impulsive and tiny, not a thrust-duty-cycle problem. The
"spawn-N-virtual-trajectories-and-add-at-worst-violation" iteration scheme is a
mildly interesting constraint-refinement pattern but not worth a mine.

---

## 2. Sinha & Beeson (2025) — Initial-guess generation, robust low-thrust

**Citation:** A. Sinha, R. Beeson, "Initial Guess Generation for Low-Thrust
Trajectory Design with Robustness to Missed-Thrust-Events," arXiv:2501.06694
(math.OC), 2025.

**Summary (3 lines):** Compares two initial-guess strategies for global search of
robust low-thrust trajectories: a baseline *non-conditional* global search
(sample from a fixed, global-support distribution) and a *conditional* global
search that seeds the k-robust problem from solutions of a less-robust (k′-robust
or non-robust) problem, giving a sequential ladder of increasingly robust solves.
Validated on the Lunar Gateway PPE; reports statistically significant gains in
convergence rate and solution quality. Frames MBH/global-search basin-shrinkage
in high-dimensional spaces as the core difficulty.

**Verdict: REFERENCE-ONLY.** The conditional/sequential-seeding *idea* (warm-start
a harder global search from an easier solved problem) is conceptually adjacent to
our seed-ladder and to continuation, and the intro is a clean modern survey of
MBH/global-search basin issues citing Englander & Englander, Wales-Doye, Leary,
adaptive-hop distributions. But the concrete method is bound to low-thrust
robust-MTE transcription and the "depth of robustness" ladder has no impulsive
analogue — the transferable kernel is already covered by our MBH (Thread 1) and
continuation (Thread 2) survey threads. Cite if we write up seed-ladder /
warm-start rationale; do not mine.

---

## 3. Zhou, Armellin, Qiao & Li (2025) — Single-impulse reachable set via polynomials

**Citation:** X. Zhou, R. Armellin, D. Qiao, X. Li, "Single-Impulse Reachable Set
in Arbitrary Dynamics Using Polynomials," arXiv:2502.11280 (astro-ph.IM), 2025.

**Summary (3 lines):** Computes the reachable set (RS) after a single velocity
impulse of arbitrary direction (fixed epoch, bounded |Δv|), under *arbitrary*
dynamics — explicitly including CRTBP / three-body, not just two-body. Uses
Differential Algebra (DACE/DACEyPy) + automatic domain splitting to build
high-order Taylor polynomials of the final state in the two impulse-direction
angles, then envelope theory to extract the RS boundary by root-finding a
one-variable polynomial; a local polynomial approximation cuts envelope-solve
time >84%. Demonstrated on two NRHOs with <0.0658% relative error; covers both
state-space and observation-space RS.

**Verdict: FULL-MINE WORTHY** — as candidate tooling for our maintenance/TCM
window analysis under n-body dynamics.

**A full mine should extract:**
- The DA + automatic-domain-split (ADS) recipe for propagating a bounded
  single-impulse Δv (fixed epoch, arbitrary direction) into a reachable-set
  boundary — directly applicable to "where can one bounded maintenance/TCM burn
  put the cycler, and does the next node lie inside that set?" under our n-body
  propagator.
- The envelope-equation → one-variable-polynomial root-finding reduction and the
  high-order local-polynomial approximation (the >84% speedup), since our solves
  are ephemeris/propagation-bound (state() ≈ 87-89%).
- Tooling/licence note: DACE (https://github.com/dacelib/dace) C++ core +
  DACEyPy (https://github.com/giovannipurpura/daceypy) Python wrapper — assess as
  a possible dependency; check licences before adoption.
- The integration setup actually used (RK78, rel/abs tol 1e-12) and the NRHO
  accuracy benchmark as a validation target if we re-implement.
- Caveat to verify on mine: it solves only the *fixed-epoch, single-impulse*
  scenario; multi-burn / free-epoch maintenance windows would need extension, and
  the DA Taylor expansion degrades over long multi-rev arcs (same long-arc
  caveat the primer-vector survey flags).

---

## 4. Takao (2025) — Mission analysis for Saturn Trojan 2019 UO14

**Citation:** Y. Takao, "Mission Analysis for the First-Ever Saturn Trojan
2019 UO14," arXiv:2501.06586 (astro-ph.EP), 2025.

**Summary (3 lines):** Designs flyby and rendezvous missions to a highly inclined
(>32°) Saturn Trojan. The ballistic phase uses a Multiple Powered Gravity-Assist
with one Deep-Space Maneuver (MPGA-1DSM) patched-conic transcription: per-leg
decision vector [V∞, azimuth, elevation, ToF, DSM-timing-fraction η], departure
epoch, with body states from SPICE; Lambert per leg, propagate-then-Lambert for
the DSM, and flyby feasibility via in/out V∞ deflection vs periapsis limit,
permitting powered flybys (Oberth) where the bend deficit allows. Solved by a
meta-heuristic global optimizer; a scaled low-thrust variant is added afterward.

**Verdict: FULL-MINE WORTHY** — for its impulsive multi-arc / multi-flyby
global-search transcription, transferable to our multi-arc cycler topology search.

**A full mine should extract:**
- The exact MPGA-1DSM decision-vector parameterisation (Eqs. 1-2:
  [V∞, α, β, τ, η] per leg + t0) and how it strings legs (Eq. 3 epoch chaining
  off SPICE body states) — this is a concrete, ballistic, ephemeris-based
  multi-arc-with-DSM genome, exactly the topology class our M-ED / S1L1 frontier
  needs (single-ellipse-per-leg is the known blocker; this is the multi-arc
  alternative).
- The powered-flyby model: deflection δ_in/δ_out from eccentricity (Eqs. 8-9),
  max-deflection feasibility test, and the Oberth powered-flyby ΔV formula
  (Eq. 11) — directly comparable to our `hunt_vem_ballistic.py` bend-feasibility
  prune, and lets a flyby absorb a bend deficit instead of being pruned.
- The DSM construction (propagate η·τ ballistically, then Lambert to next body,
  ΔV = velocity mismatch at the DSM point) — a clean recipe for inserting an
  interior impulse into a leg, relevant to multi-arc return modelling.
- Which meta-heuristic global optimizer and decision-vector bounds were used, for
  comparison against our grid + MBH plans.
- Note for the catalogue/scope: the science target is irrelevant (Saturn Trojan);
  the value is purely the impulsive transcription + powered-flyby/DSM machinery.

---

## 5. Hu, Yang, Li & Baoyin (2024) — Robust low-thrust gravity-assist via RL

**Citation:** J. Hu, H. Yang, S. Li, H. Baoyin, "Robust Design of Low-Thrust
Gravity-Assist Trajectories via Reinforcement Learning," Journal of Guidance,
Control, and Dynamics (AIAA DOI 10.2514/1.G009427), 2024.

**Summary (3 lines):** RL-based multi-phase robust guidance for low-thrust gravity
-assist under state/observation/control/missed-thrust uncertainty. Segments legs
into interplanetary-transfer phases (ITPs) and approaching phases (APs), models
each as an MDP, uses a trade-off factor to modify the AP nominal initial state,
and folds an analytical reachability constraint into the ITP reward to improve
AP-target reachability. Demonstrated on an Earth-Earth-Jupiter low-thrust mission.

**Verdict: OUT-OF-SCOPE.** This is a low-thrust-only ML/RL closed-loop *guidance
law* trained as a neural policy under uncertainty — precisely the category the
brief excludes regardless of quality. The "analytical reachability constraint in
the reward" is the only conceptually adjacent fragment, but it is tied to the RL
reward shaping and far weaker than the dedicated RS method in paper #3. No mine.

---

## 6. Blender & Singh (2025) — Uncertainty-aware GBDT guidance, continuous thrust

**Citation:** S. Blender, S. K. Singh, "Uncertainty-Aware Guidance for Continuous
Thrust Transfers Using Gradient-Boosted Decision Trees," AAS 25-524, AAS/AIAA
Astrodynamics Specialist Conference, 2025.

**Summary (3 lines):** Builds an uncertainty-aware dataset by propagating
sigma-point (Unscented Transform) ensembles around an extremal-bundle of nominal
continuous-thrust trajectories, capturing the joint state/control distribution.
Trains LightGBM gradient-boosted decision trees to map current belief state to
optimal control-to-go corrections. Demonstrated on an Earth-to-3671-Dionysus
time-optimal low-thrust transfer; reports reduced tracking error / control effort
at real-time speed vs nominal-only guidance.

**Verdict: OUT-OF-SCOPE.** A low-thrust continuous-thrust real-time ML guidance
law — the exact excluded category. The extremal-bundle + sigma-point
uncertainty-propagation dataset construction (after Izzo et al.) is mildly
interesting as a TCM-covariance idea, but it is indirect-optimal-control /
costate-bundle machinery for continuous thrust with no impulsive-cycler payoff
that the primer-vector (Thread 3) and RS (paper #3) routes don't cover better.
No mine.

---

## Verdict summary

| # | Paper | Verdict |
|---|-------|---------|
| 1 | Venigalla, Englander & Scheeres (2020), MTM virtual swarm, AAS 20-438 | OUT-OF-SCOPE |
| 2 | Sinha & Beeson (2025), conditional initial-guess gen, arXiv:2501.06694 | REFERENCE-ONLY |
| 3 | Zhou, Armellin, Qiao & Li (2025), single-impulse RS via polynomials/DA, arXiv:2502.11280 | FULL-MINE WORTHY (TCM/maintenance reachability tooling, n-body) |
| 4 | Takao (2025), Saturn Trojan MPGA-1DSM mission analysis, arXiv:2501.06586 | FULL-MINE WORTHY (impulsive multi-arc + powered-flyby/DSM transcription) |
| 5 | Hu, Yang, Li & Baoyin (2024), robust low-thrust GA via RL, AIAA G009427 | OUT-OF-SCOPE |
| 6 | Blender & Singh (2025), GBDT continuous-thrust guidance, AAS 25-524 | OUT-OF-SCOPE |
