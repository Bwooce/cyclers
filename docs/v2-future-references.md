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
| Cached PDF | not yet downloaded |

**Why it matters.** The Sims-Flanagan transcription is the standard way to
turn a continuous low-thrust trajectory optimisation problem into a finite
nonlinear program (NLP) the optimiser can solve. Each interplanetary leg
is sliced into N segments; the low-thrust history is treated as a sequence
of impulsive ΔV at segment boundaries, connected by Keplerian conic arcs.
This is the architectural pattern cyclerfinder would adopt to handle
low-thrust cyclers.

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

## 2. Pascarella, Woollands, Pellegrini, Sanchez-Net, Van Hook 2024 — Solar System Pony Express

| Field | Value |
|---|---|
| Title | "Low-thrust trajectory optimization for the Solar System Pony Express" |
| Authors | Andrea Pascarella, Robyn Woollands, Etienne Pellegrini, Marc Sanchez-Net, Joel Van Hook |
| Venue | *Advances in the Astronautical Sciences*, 2024, pp. 45-61 |
| DOI | [10.1007/978-3-031-51928-4_4](https://doi.org/10.1007/978-3-031-51928-4_4) |
| Cached PDF | not yet downloaded |

**Why it matters.** Pascarella et al. outline the **patched-conic → medium-fidelity
impulsive → high-fidelity N-body** pipeline for transitioning a candidate
low-thrust cycler trajectory from a circular-coplanar idealisation to a
real-ephemeris flight design. The paper documents the ΔV penalty incurred
by feeding patched-conic outputs *directly* into a high-fidelity optimiser
(spoiler: massive, due to solar gravity perturbations during flybys) and
shows an intermediate medium-fidelity impulsive stage that absorbs most of
the divergence.

**Architectural impact on cyclerfinder if adopted:**

- This is the **canonical reference for M6b's ephemeris-mode
  optimisation pipeline**. The current M6 slice (`9b2611d`) implemented
  geometric phase-matching only; M6b's full implementation should follow
  the Pascarella three-stage pipeline rather than going directly from
  idealised seed → high-fidelity N-body.
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
| Cached PDF | `/tmp/arxiv-1511.00821.pdf` (1.1 MB) — downloaded 2026-06-01 |

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
| Cached PDF | `/tmp/arxiv-2305.18368.pdf` (835 KB) — downloaded 2026-06-01 |

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

## 5. Hollister & Menning 1969-1971 lineage — pre-Aldrin foundational

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

- **Sims-Flanagan / low-thrust v2** — Papers 1, 3, 4 cited; new
  `trajectory_regime: low-thrust` entries ingested.
- **M6b ephemeris-mode pipeline** — Paper 2 referenced as the
  architectural template.
- **Hollister & Menning lineage backfill** — Paper 5 numerics ingested.

Until then this document is awareness-only. Updating: when the project
adopts v2, this doc moves to `docs/phases/v2-low-thrust/references.md` (or
similar) and gains an implementation plan alongside it.
