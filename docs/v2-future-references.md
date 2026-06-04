# v2 future references

This document is the seed bibliography for the **v2 scope expansion** of
cyclerfinder beyond ballistic patched-conic + impulsive flybys (the
current v1 scope per spec §2). The papers here are queued for ingest when /
if the project formally adopts low-thrust / Sims-Flanagan / multi-fidelity
optimisation — at which point individual cyclers from these works become
catalogue entries with `trajectory_regime: low-thrust` or
`model_assumption: analytic-ephemeris` per the schema v2 conventions
(`data/README.md`).

For each entry: full citation, why it matters to the cyclerfinder
architecture, what specifically would change in the codebase to consume it,
and whether a PDF is cached locally.

---

## 1. Yam, Di Lorenzo, Izzo 2010 — Sims-Flanagan low-thrust transcription

| Field | Value |
|---|---|
| Title | "Constrained global optimization of low-thrust interplanetary trajectories" |
| Authors | Chit Hong Yam, Dario Di Lorenzo, Dario Izzo |
| Venue | *IEEE Congress on Evolutionary Computation*, 2010, pp. 1–7 |
| DOI | [10.1109/cec.2010.5586019](https://doi.org/10.1109/cec.2010.5586019) |
| Full text | held offline (cite by DOI `10.1109/cec.2010.5586019`; not stored in repo) |

**v1 catalogue eligibility: OUT OF SCOPE (low-thrust methods paper).** Read
in full 2026-06-01. The paper contains no ballistic cycler trajectory and
no concrete published numbers that fit the v1 circular-coplanar patched-
conic schema. Its two worked examples are one-way **rendezvous** transfers
to Jupiter (the arrival velocity is *constrained to zero relative to the
planet* — "we are considering a rendezvous problem"), not repeating
cyclers, and both require continuous nuclear-electric thrust to close. No
`catalogue.yaml` change is warranted.

**Why it matters.** The Sims-Flanagan transcription is the standard way to
turn a continuous low-thrust trajectory optimisation problem into a finite
nonlinear program (NLP) the optimiser can solve. The paper's own framing:
"Trajectory is divided into legs which begin and end with a planet.
Low-thrust arcs on each leg are modeled as sequences of impulsive maneuvers
∆V, connected by conic arcs." Each leg is propagated forward and backward
to a *matchpoint* (usually halfway), and the position/velocity mismatch
`Smf − Smb` must fall below tolerance to be feasible. This is the
architectural pattern cyclerfinder would adopt to handle low-thrust
cyclers.

**Concrete specifics from the paper text (2026-06-01 read):**

- **Problem class = LT-MGA** (Low-Thrust Multiple Gravity Assist), framed as
  a constrained NLP. The paper's contribution over box-constrained MGA /
  MGA-DSM is handling the "high number of nonlinear constraints" via
  hybridising the global search with a *local* solver, rather than via an
  objective-function penalty term with an unknown weighting factor.
- **NLP dimension** = `(8 + 3N)M` for `M` legs, `N` segments/leg, with
  ~`neq·M` nonlinear constraints (`neq = 7` for a 3-D + mass problem). The
  ∆V-per-segment bound is `∆Vmax = (Tmax/m)(tf − t0)/N` (Eq. 1); the flyby
  turn-angle constraint uses `sin(δ/2) = 1/(1 + rp·V∞²/µ)` (Eq. 3) — the
  same bend formula already noted for the v2 `flyby_mechanics` field in
  `data/README.md`.
- **Two-phase solve.** Phase 1 minimises total ∆V at *constant mass*
  (`min ΣΣ∆Vi`, Eq. 4); Phase 2 re-optimises locally to *maximise final
  mass* `mf` (Eq. 6) propagating mass via the rocket equation
  `mi+1 = mi·exp(−∆Vi/g0·Isp)` (Eq. 5). The local solver is **SNOPT** (SQP).
- **Three global optimisers compared**: Multistart (MS, baseline), Monotonic
  Basin Hopping (BH, a.k.a. Iterated Local Search; MNI=500, perturbation
  r=0.05, time-shift probability p=0.1 by a *synodic period*), and Simulated
  Annealing with Adaptive Neighborhood (SA). **Basin Hopping wins** on both
  test cases (statistically significant by unpaired t-test, both quality and
  number of feasible solutions).
- **Test cases — NOT GTOC problems**, but a NEP mission inspired by the
  cancelled NASA **Jupiter Icy Moons Orbiter** (Table I: initial mass
  20,000 kg, Tmax 2.26 N, Isp 6,000 s, launch window Jan 2020–Jan 2030,
  launch V∞ ≤ 2.0 km/s, min flyby radius 7,000 km; max ToF 10 yr for E-E-J,
  15 yr for E-E-E-J). GTOC/GTOP are only cited as motivating context.
  - **E-E-J** (one Earth flyby): NLP dim 75, 35 nonlinear constraints.
    Phase-1 best total ∆V ≈ 9.558 km/s (BH); Phase-2 best final mass
    ≈ 17,004 kg (BH). Best trajectory: launch 2022-10-15 at V∞ 2 km/s, a
    ~1.5-yr 3:2 resonance loop (2 revs), Earth flyby at 2.9 yr boosting V∞
    to 9.0 km/s, Jupiter rendezvous ~Sept 2029.
  - **E-E-E-J** (two Earth flybys): NLP dim 112, 56 nonlinear constraints.
    Phase-1 best total ∆V ≈ 7.524 km/s (BH); Phase-2 best final mass
    ≈ 17,601 kg (BH). Best trajectory: 1:1 resonance (495 d, V∞ 2.0→5.3
    km/s), then 2:1 resonance (772 d, V∞ →9.0 km/s), total flight ~8.0 yr;
    the extra flyby buys ~600 kg final mass at the cost of ~1 yr ToF.
- **Why none of this is a catalogue entry**: the resonance loops are
  one-shot capture-to-Jupiter sequences, not steady-state repeating cyclers;
  and (per the paper) "unlike the case in the ballistic transfer, in the
  low-thrust transfer case the Earth-Earth flight time does not exactly
  equal some integer multiple of Earth's orbital period" — i.e. these are
  expressly NOT the ballistic resonant geometry the v1 schema encodes.

**Architectural impact on cyclerfinder if adopted:**

- `core/lambert.py` is no longer the right leg constructor for low-thrust
  segments. Replace single-rev Lambert with a shape-based seed (exponential
  sinusoid; see Petropoulos 2003) or directly with Sims-Flanagan as the
  initial leg discretisation.
- `model/cycler.py` `Leg` gains a `segments: tuple[Segment, ...]` field
  carrying the per-segment ΔV vectors. v1 ballistic legs degenerate to
  `len(segments) == 0` (the impulsive ΔV are at encounter boundaries, not
  on the leg).
- `search/optimize.py` adds a low-thrust mode: SLSQP / IPOPT solves the NLP
  with per-segment ΔV as decision variables and a thrust-magnitude bound
  per segment as an inequality constraint.
- `model/score.py` adds `propellant_mass_fraction` or similar so low-thrust
  candidates are comparable to ballistic ones.

## 2. Pascarella, Woollands, Pellegrini, Sanchez-Net, Van Hook 2024 — Solar System Pony Express (**low-thrust**)

| Field | Value |
|---|---|
| Title | "Low-thrust trajectory optimization for the Solar System Pony Express" |
| Authors | Andrea Pascarella, Robyn Woollands, Etienne Pellegrini, Marc Sanchez-Net, Joel Van Hook |
| Venue | *Advances in the Astronautical Sciences*, 2024, pp. 45-61 (conference paper AAS 22-015) |
| DOI | [10.1007/978-3-031-51928-4_4](https://doi.org/10.1007/978-3-031-51928-4_4) |
| Full text | held offline (cite by DOI `10.1007/978-3-031-51928-4_4`; AAS 22-015 open-access preprint from `ai.jpl.nasa.gov/public/documents/papers/AAS-22-015-Paper.pdf`; not stored in repo) |

> **Distinct-paper note.** This is the **low-thrust** Pony Express paper (Pascarella AAS-22-015 / 2024 book chapter). It is DISTINCT from Sanchez Net, Pellegrini, Parker, Vander Hook, Woollands 2022 *Journal of Spacecraft and Rockets* 59(3):861-870 (DOI 10.2514/1.A35091), which is the **near-ballistic** (ΔV≤10 m/s) Pony Express paper and is in-scope for v1. See `data/OUTSTANDING.md` task #38 for the near-ballistic entry.

> **Source-version note.** The 2024 *Advances in the Astronautical Sciences*
> book chapter (DOI `10.1007/978-3-031-51928-4_4`) sits behind Springer's
> paywall (`idp.springer.com` auth redirect — not accessible to the web-fetch
> tool). The cached PDF is the open-access **AAS 22-015** conference paper
> (same title, same five authors — note the PDF byline spells the last author
> "Joshua Vander Hook" vs our citation's "Joel Van Hook"; same JPL group
> supervisor). AAS conference proceedings are what get published as *Advances
> in the Astronautical Sciences*, so the chapter is the proceedings version of
> this paper. All quotes below are verbatim from the AAS 22-015 PDF; treat
> numerics as high-confidence pending verification against the paginated
> chapter.

**Why it matters.** Pascarella et al. design **Earth-Mars cycler orbits for
the SSPE** (Solar System Pony Express) — a JPL NIAC concept (grant
80NM0018D0004) for interplanetary **data mules**: 500 kg ESPA-class smallsats
with optical-laser comms and a NEXT ion engine (Isp 4155 s, thrust 0.235 N),
launched as rideshare with a Mars-bound mission, that fly Earth-Mars cyclers,
retrieve "1-3 petabits of data per flyby" at Mars and downlink at Earth
(>8000 Tbits over 8 flybys in their example). It is the canonical reference
for the **patched-conic → ephemeris transition** — but the actual pipeline is
**five sequential sub-problems**, not three, and the optimiser is **indirect
optimal control** (Pontryagin / primer-vector, MEE state, bang-bang + tanh
continuation, RK 9(8) + Matlab `fsolve` multiple-shooting), **not**
Sims-Flanagan (contrast Paper 1). Verbatim: "the trajectory design is divided
in five sub-problems that are solved sequentially":

1. **Patched-conic** — JPL's STAR software generates "thousands of Earth-Mars
   cycler trajectories with a patched conic model"; filtered on ΔV +
   flyby-altitude requirements.
2. **Two-body impulsive** — a TPBVP in two-body dynamics builds planetocentric
   flyby trajectories matching STAR's in/out velocity vectors "with negligible
   difference."
3. **Medium-fidelity impulsive** — adds solar-system gravity perturbations.
   The documented penalty lives here (verbatim): *"the Sun's gravity
   significantly perturbs the trajectory of the spacecraft during the planetary
   flybys and thus the incoming and outgoing vectors computed by STAR result in
   very large ∆v's for the flyby targeting problem."* Fix: propagate
   backward/forward from each flyby and minimise position/velocity mismatch at a
   heliocentric "breakpoint" midway between flybys.
4. **High-fidelity low-thrust** — step-3 impulsive solutions seed the indirect
   low-thrust TPBVP; impulsive ΔV replaced by optimal thrust/coast arcs. "the
   high-fidelity impulsive solutions are a very good initial guess ... however
   the run-time is significantly longer." They could **not** converge the whole
   STAR set in the ephemeris model — 278 reached step 3, only "several" through
   step 5.
5. **Polynomial thrust-arc fitting** — post-process for precise thruster on/off
   times for ops; negligible fuel difference vs step 4.

**Documented ΔV / propellant numbers (verbatim).** Headline result: a "500 kg
courier spacecraft ... inserted into an Earth-Mars cycler orbit using only
36 kg of propellant, and a further 2 kg of propellant is required to target
eight subsequent flybys over a period of six years." Cost is dominated by
**cycler-orbit injection (COI)**, not maintenance: "most of the propellant
budget is required for COI, and once the spacecraft is on the cycler orbit a
minimal amount of fuel (< 5 kg) is needed to maintain the orbit and target
flybys." Russell-Ocampo 2006 (their ref [13]) is cited for the claim that
"directly obtaining cycler orbits in an ephemeris model is often impractical
if not impossible" — the explicit justification for the graduated
multi-fidelity ladder. Flyby altitude band ≤ 25,000 km, lower bounds 300 km
(Mars) / 1000 km (Earth); engine off during flybys for comms attitude; flyby
target set at 3× Mars SOI so the flyby stays ballistic.

**Architectural impact on cyclerfinder if adopted:**

- This is the **canonical reference for M6b's ephemeris-mode
  optimisation pipeline**. The current M6 slice (`9b2611d`) implemented
  geometric phase-matching only; M6b's full implementation should follow
  the Pascarella graduated multi-fidelity ladder (five steps: patched-conic →
  two-body impulsive → medium-fidelity impulsive → high-fidelity low-thrust →
  polynomial fit) rather than going directly from idealised seed →
  high-fidelity N-body. M6b's v1-scope ballistic mandate maps onto steps 1-3
  (the impulsive ladder); steps 4-5 (low-thrust) are explicitly out of v1
  scope. Note the paper's optimiser is *indirect* optimal control, whereas
  the M6b plan §3.1 chooses a Lambert-chain construction — the alignment
  Pascarella validates is the *fidelity ladder*, not the solver.
- `verify/propagate.py` (M6a, planned) becomes the medium-fidelity
  impulsive stage. `verify_long_term_stability` produces the per-lap TCM
  ΔV that quantifies the patched-conic → medium-fidelity penalty.
- A future `optimize/ephemeris_mode.py` (M6b body) does the
  medium-fidelity → high-fidelity step by re-optimising with a
  multi-segment impulsive sequence, then handing off to the optional GMAT
  V4 gate for the final N-body verification.

## 3. Izzo, Hennes, Simões, Märtens 2015 — GTOC trajectory design

| Field | Value |
|---|---|
| Title | "Designing Complex Interplanetary Trajectories for the Global Trajectory Optimization Competitions" |
| Authors | Dario Izzo, Daniel Hennes, Luís F. Simões, Marcus Märtens |
| Venue | arXiv preprint (ESA Advanced Concepts Team) |
| arXiv | [1511.00821](https://arxiv.org/abs/1511.00821) |
| Full text | held offline (cite by arXiv `1511.00821`; not stored in repo) |

**Why it matters.** Izzo's group at ESA-ACT runs the Global Trajectory
Optimization Competition (GTOC), the premier benchmark for interplanetary
trajectory design. This paper reviews the methods used across the
competitions, including the Sims-Flanagan transcription, the dynamical
systems approach (manifolds), and the algorithmic-evolutionary methods
(differential evolution at scale). It's the standard reference for
**"how does the trajectory-design community actually find the best
trajectories, given a budget of compute?"**.

**Architectural impact on cyclerfinder:**

- Validates the project's current scipy DE choice (Izzo's lab built the
  pygmo / pagmo framework around DE for these competitions; spec §7 lists
  pygmo as an optional v2 extra).
- The paper's archipelago migration topologies + island-based search are
  the natural M7 parallelisation pattern for the §13.6 work queue: each
  cell becomes an island, with periodic migration of best solutions
  between cells to share basin information.
- Catalogue ingest: the GTOC winning trajectories themselves are
  candidates for catalogue entries (with `model_assumption:
  analytic-ephemeris` since GTOC uses mid-fidelity ephemerides).

## 4. Burhani, Fantino, Flores, Sanjurjo-Rivo 2023 — automated inclined low-thrust

| Field | Value |
|---|---|
| Title | "A new automated strategy for optimizing inclined interplanetary low-thrust trajectories" |
| Authors | Burhani M. Burhani, Elena Fantino, Roberto Flores, Manuel Sanjurjo-Rivo |
| Venue | arXiv preprint, 2023 |
| arXiv | [2305.18368](https://arxiv.org/abs/2305.18368) |
| Full text | held offline (cite by arXiv `2305.18368`; not stored in repo) |

**Why it matters.** Same Fantino (UAE University) as the Saturnian CR3BP
work flagged in `data/OUTSTANDING.md` §H — but here from a low-thrust
angle. The paper presents an automated pipeline for optimising
*inclined* interplanetary low-thrust trajectories, i.e. ones that depart
significantly from the ecliptic. Important because most cyclers are
catalogued in the circular-coplanar (i=0) idealisation, but real Mars and
Venus orbital inclinations (1.85° and 3.4°) introduce real out-of-plane
ΔV cost that current cyclerfinder treats as zero.

**Architectural impact on cyclerfinder:**

- M6b's TCM-budget computation should include out-of-plane corrections;
  this paper's automated strategy is the right reference for how to do
  it without over-fitting individual launch dates.
- Once `model_assumption: analytic-ephemeris` entries gain
  `inclination_deg != 0` data, the matcher will need an inclination
  tolerance band — Burhani's empirical numbers help calibrate it.

## 5. Ozimek, Riley, Arrieta 2019 — LInX low-thrust trajectory-optimization tool

| Field | Value |
|---|---|
| Title | "The Low-Thrust Interplanetary eXplorer (LInX): A Medium-Fidelity Algorithm for Multi-Gravity Assist Low-Thrust Trajectory Optimization" |
| Authors | Martin T. Ozimek, James Riley, Arturo Arrieta |
| Venue | *Advances in the Astronautical Sciences*, AAS 19-348, 2019 (JHU APL + Nabla Zero Labs) |
| Full text | held offline (cite by conference number AAS 19-348; not stored in repo) |

**v1 catalogue eligibility: OUT OF SCOPE (low-thrust trajectory-optimization tool paper).** No ballistic
cycler trajectory or patched-conic V∞ anchor is presented. This is a methods contribution describing
the LInX software tool (Sims-Flanagan transcription + SPICE ephemerides + NLP solver), not a cycler
catalogue source.

**Why it matters.** LInX extends the Sims-Flanagan framework (cf. Paper 1, Yam 2010) with SPICE
planetary data and multi-gravity-assist structure in a single medium-fidelity package. It is the
canonical reference for the JHU APL / Nabla Zero Labs low-thrust toolchain and a peer of the JPL
STAR software used in Pascarella (Paper 2). If the project adopts low-thrust v2 modelling and seeks
to benchmark or replicate JHU APL results, this paper is the primary reference for LInX.

**Architectural impact on cyclerfinder if adopted:** same class as Paper 1 (Yam 2010) — provides
the NLP formulation and SPICE integration pattern for a future `search/lt_optimizer.py`
low-thrust leg constructor. Not actionable until v2 scope is formally adopted.

## 6. Hollister & Menning 1969-1971 lineage — pre-Aldrin foundational

| Field | Value |
|---|---|
| Title | Multiple papers, 1969-1971 (see catalogue entry `hollister-menning-1970-ev-periodic` for the canonical 1970 citation) |
| DOIs | 10.2514/3.29664 (Hollister 1969 precursor), 10.2514/3.30134 (Hollister & Menning 1970 E-V), plus Rall & Hollister 1971 Earth-Mars companion |
| Cached PDF | none (all attempts blocked by AIAA paywall) |
| Catalogue | `hollister-menning-1970-ev-periodic` (citation-only — all numerics null) |

**Why it matters.** Pre-dates Aldrin by 15 years. First analytical proof
that a spacecraft can indefinitely cycle between two planets via matched
V∞ swing-bys. The matched-V∞ constraint is the analytical ancestor of
every Tisserand-graph method later used in the catalogue. The papers are
foundational but inaccessible — if a future ingest can obtain them
(institutional access; MIT archives; etc.), the catalogue entry's null
numerics get backfilled and the Earth-Venus cycler family gains
per-member data.

**Architectural impact:** none structurally — the Tisserand machinery
(`search/tisserand.py`) already implements the matched-V∞ contour logic.
The papers are valuable as primary-source attribution and possible
discovery of additional V-E cyclers not yet in the catalogue.

---

## Activation

Each of these papers becomes "real catalogue work" once the project
formally scopes:

- **Sims-Flanagan / low-thrust v2** — Papers 1, 3, 4, 5 cited; new
  `trajectory_regime: low-thrust` entries ingested.
- **M6b ephemeris-mode pipeline** — Paper 2 referenced as the
  architectural template.
- **Hollister & Menning lineage backfill** — Paper 6 numerics ingested.

Until then this document is awareness-only. Updating: when the project
adopts v2, this doc moves to `docs/phases/v2-low-thrust/references.md` (or
similar) and gains an implementation plan alongside it.
