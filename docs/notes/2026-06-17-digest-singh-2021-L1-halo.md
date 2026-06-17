# Digest: Singh, Anderson, Taheri & Junkins 2021, "Exploiting manifolds of L1 halo orbits for end-to-end Earth-Moon low-thrust trajectory design"

**Task**: #353 (low-thrust digest Agent C, paper 2 of 4).
**Read**: 2026-06-17, all 18 pages.

## Bibliographic

- **Authors**: Sandeep K. Singh (TAMU, PhD candidate), Brian D. Anderson (JPL),
  Ehsan Taheri (Auburn), John L. Junkins (TAMU).
- **Venue**: Acta Astronautica 183 (2021) 255-272.
- **DOI**: 10.1016/j.actaastro.2021.03.017.
- **Received**: 23 May 2020; accepted 21 March 2021; online 3 April 2021.

## What the paper actually is

An end-to-end Earth-Moon mission design study for a low-thrust smallsat (240 kg
wet, 0.096 N max thrust, 1250 s Isp) going from a GTO drop-off to a 1000 km
circular lunar polar orbit (LPO). The vehicle of transport is the **stable-
unstable manifold pair of an Earth-Moon L1 halo orbit**: spacecraft enters via
the stable manifold (Earth side) and exits via the unstable manifold (Moon
side). The paper has two halves:

1. **CR3BP design pass** (Secs 2-6, 8-9): pick the L1 halo, propagate its
   manifolds, find favourable "piercing" points on two cutting planes
   (x = −μ for the Earth leg, x = 1 − μ for the Moon leg), pick IC/FC that
   minimise inclination change at GTO and at lunar capture.
2. **BCP perturbation pass** (Sec 7): re-evaluate the CR3BP-designed manifold
   ICs/FCs in a Bi-Circular Problem (Earth-Moon-Sun) and homotope over solar
   phase angle λ₀ to find solar phasings that preserve manifold geometry.

The mission is split into 7 phases (Table 1, p.257): GTO → nGEO (min-time
spiral, MEE), 30-day coast in nGEO (SSA), nGEO → L1 halo stable manifold
(min-time), coast on manifold, halo waypoint, manifold off → moon side, then
manifold → 1000 km lunar polar orbit (min-time MEE in MCI). 48+ revolutions
in the Earth leg, 3-day powered insertion on the Moon side.

## Method digest

- **Optimisation framework**: indirect minimum-time TPBVP per phase, solved
  piecewise with MATLAB `fsolve` (Levenberg-Marquardt). 368 s total wall time
  on Intel Xeon Gold 6130, 32 GB RAM, for the full end-to-end converge.
- **Coordinates**: MEE (Modified Equinoctial Elements) for the multi-rev
  spiral phases 1A and 2C — exploits MEE's regularity (consistent with paper
  1, Junkins-Taheri 2019). Cartesian synodic CR3BP for the manifold-adjacent
  phases.
- **Smoothing**: continuation on boundary conditions (homotopy via η in Eq 8
  of the paper, η: 0 → 1 walks an easy initial-guess solution to the desired
  final BC). This is **boundary-condition** homotopy, not throttle smoothing —
  different from the HTSL homotopy in paper 1.
- **Manifold construction**: standard linear-eigenvector seed at ε = 1e-6 on
  100 points around the halo orbit (Eq 5, p.258). Stable/unstable manifolds
  cut at planes x = −μ (Earth, IC for stable) and x = 1 − μ (Moon, FC for
  unstable), Sec 6 / Fig 8.
- **BCP integration**: Earth-Moon synodic frame with Sun added (Eq 6 p.260);
  solar phase angle λ₀ parametrises the time dependence. Manifolds re-
  propagated for 29.618 days; the comparison metric is osculating-element
  deviation between CR3BP and BCP propagations (Figs 14, 17).
- **Integrator**: not explicitly stated; standard MATLAB ode45-class tooling.

## **Sourced halo ICs (Table 2, p.257) — KNOWN_CORPUS candidates**

This is the load-bearing data extract. Earth-Moon CR3BP, EM synodic frame,
ND units (DU = 386 274.56245094 km, μ = 0.012150586632602, VU = 1.0281468...
km/s, TU = 375 699.807437... s):

| Jacobi C | x (DU)             | y (DU) | z (DU)            | v_x (VU) | v_y (VU)         | v_z (VU)           |
|----------|--------------------|--------|-------------------|----------|------------------|--------------------|
| 3.128    | 0.990617793994954  | 0      | 0.126355564698763 | 0        | 0.019697666322050 | 0                  |
| 3.143    | 0.989736590417966  | 0      | 0.118345929293312 | 0        | 0.015959748108035 | 0                  |
| 3.158    | 0.989169523068388  | 0      | −0.111126735244195| 0        | 0.013230252533363 | 0                  |

All three orbits have y = v_x = v_z = 0 at the planar crossing (standard
periodic-orbit IC convention for halos with z-symmetry). All three are
Earth-Moon L1 halos, bifurcating from the planar Lyapunov family at L1.

**Selected for the mission**: C = 3.128 (largest amplitude, ~8.4° max inertial
inclination achievable on piercing plane Earth leg, Fig 9).

**Manifold IC/FC adopted in CR3BP** (Table 3, p.263):
- Stable manifold IC (Earth leg, x = −μ piercing): x = −0.0121505866303,
  y = 0.339432988876685, z = 0.0198732432593505,
  v_x = 1.655892806055449, v_y = 0.061882094901634, v_z = 0.274216573159948.
- Unstable manifold FC (Moon leg, x = 1 − μ piercing):
  x = 0.987849413374523, y = 0.001704102052949, z = 8.149831565347900e-4,
  v_x = 1.489681663332230, v_y = 2.153172331531527, v_z = 2.420970920376801.

**BCP-corrected IC/FC at λ₀ = 32°** (Table 5, p.265): different by ~5-8%
in position; the Sun perturbation walks the manifold-piercing condition off
the CR3BP value.

## Relevance to this project

### To #309 (Sims-Flanagan low-thrust substrate) — **moderate, not direct extension**

- The paper is **indirect-method TPBVP**, not SF transcription. It does not
  directly improve `cyclerfinder.core.sims_flanagan`.
- However, the paper's **two-stage decomposition pattern** (CR3BP solution
  first, then BCP re-validation as a post-hoc check) is procedurally what
  the cyclers project does at the corpus level: design in CR3BP, validate in
  ER3BP / BCR4BP. That conceptual pattern is already adopted; nothing to
  port.
- **HTC** (high-thrust correction): paper's solar-phase BCP sweep produces
  a candidate "favourable" phase angle (λ₀ = 32°) at which manifold-piercing
  position/velocity perturbations from CR3BP are nearly minimised (Fig 14
  shows |IΔ| ≈ 1°, |ΔR_p| ≈ 10 km). If we ever build a powered-cycler
  station-keeping pass for an EM L1 halo precursor in BCP, this λ₀-sweep
  recipe could anchor the search initial guess. Currently nothing in our
  substrate does this.

### To #347 (Floquet bifurcation framework Phase 1) — **modest**

- Singh 2021 uses the standard halo-orbit Monodromy/Floquet eigenstructure
  (Eq 4, p.257): λ₁ > 1, λ₂ = 1/λ₁ (the saddle pair, source of stable/
  unstable manifolds); λ₃ = λ₄ = 1 (neutral, periodicity + energy);
  λ₅ = λ̄₆ (complex conjugate pair, quasi-periodic centre). This is
  textbook Floquet-multiplier structure for an L1 halo; no new bifurcation
  insight beyond well-known Doedel/Pernarcic-Howell results.
- **No saddle-centre or period-doubling bifurcation analysis** in this
  paper — Singh et al. select isolated halo orbits and study transport, NOT
  family bifurcations. RTR2026's claim that NRHO bifurcation produces cycler
  families is **not the topic of this paper**.
- Phase 1 of #347 does not need this paper. Direct verdict: **no
  contribution to the Floquet substrate**.

### Catalogue / KNOWN_CORPUS impact — **VALUABLE**

**Three new sourced L1 halo IC candidates** (Table 2 above) are immediately
usable as **precursor_mga or quasi_cycler corpus anchors** in the Earth-Moon
system. These are:

- Full 15-digit synodic-frame state vectors.
- Tied to a specific published source with a DOI.
- All three are EM L1 halos at distinct Jacobi constants spanning a
  1 % family interval (3.128 ↔ 3.158).
- The C = 3.128 IC is mission-selected; the other two are also in the
  "useful for cislunar transport" envelope.

**Recommendation to KNOWN_CORPUS curator**: add three EM L1 halo anchors
under the precursor_mga (or "earth_moon_libration") category, tagged
`source: Singh-Anderson-Taheri-Junkins 2021 Acta Astronautica 183`.

The IC tables also resolve the EM L1 halo family in a 3-point grid that
could **anchor #347's Floquet substrate as a cross-check**: the published
ICs give multipliers we can independently recompute and compare.

## Catalogue impact

This paper does **not directly admit any catalogue rows**:
- The studied trajectory is a single-shot GTO → LPO mission, not a cycler
  or quasi-cycler (no repetitive periodic structure between bodies).
- The L1 halo orbit itself is a single periodic orbit, not a Mars-Earth /
  Earth-asteroid cycler.
- However, the three L1 halo ICs could **support precursor_mga class
  rows** if the project ever admits "precursor: EM L1 halo → LPO transit"
  as a precursor capability. This requires #316 cross-system framework
  treatment — defer to that workstream.

## Action items

1. **Stage three KNOWN_CORPUS additions** (EM L1 halo at C = 3.128, 3.143,
   3.158) for a future literature_check.py refresh — write **as a
   recommendation, not now** (per task instruction: do NOT edit corpus).
2. **Cross-check with #347 Phase 1**: when the Floquet substrate computes
   monodromy/multipliers for EM L1 halo, the three Singh ICs are a sourced
   ground truth for spec-validation.
3. **Note for #309 V2**: paper's λ₀ = 32° BCP-favourable-phase result is a
   non-obvious solar-phase-anchor for any future cislunar low-thrust pass
   our substrate executes in BCR4BP / BCP.
4. **Bibliographic logging**: cite alongside Howell/Barden, Mingotti et al.,
   Pan-Lu-Pan, as the standard reference for **halo-manifold-leveraged
   low-thrust** in Earth-Moon. Note this paper's contribution is the BCP
   perturbation analysis and the published IC table.

## Verdict

**Methodology + sourced ICs paper.** Off-scope for #309's direct SF substrate,
no immediate use for #347's Floquet bifurcation Phase 1 (no bifurcation
analysis here). **HIGH VALUE for KNOWN_CORPUS**: three sourced 15-digit
EM L1 halo ICs with associated stable-unstable manifold piercing-plane
IC/FC pairs at two distinct Sun phasings. Recommend three anchor additions
under earth_moon_libration / precursor_mga category in a future curator pass.
