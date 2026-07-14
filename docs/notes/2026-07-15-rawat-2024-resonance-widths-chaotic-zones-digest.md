# Rawat-Kumar-Rosengren-Ross 2024 — Resonance widths, chaotic zones, cislunar transport (AAS 24-368)

**Date filed/digested:** 2026-07-15 (light digest); **full mining pass:** 2026-07-15 (#597/#15).
**Trigger:** same `ross.aoe.vt.edu/papers/` review as the two companions above. Confirmed genuinely
new via title-page read; full 20-page body read page-by-page for this pass.

**Citation (corrected — the light digest had the wrong venue):**
> Rawat, A., Kumar, B., Rosengren, A. J., and Ross, S. D. (2024). "Resonance Widths, Chaotic Zones,
> and Transport in Cislunar Space." **AAS/AIAA Astrodynamics Specialist Conference**, Broomfield, CO,
> Aug 2024, AAS 24-368. (Not the Space Flight Mechanics Meeting as the light digest said — corrected
> against the actual PDF title page/footer.) No DOI found (AAS conference paper, not indexed in
> CrossRef as of this filing).

Registered as `CorpusAnchor` `rawat-2024-resonance-widths-aas` in
`src/cyclerfinder/search/literature_check.py` (provenance `verified-against-source`).

## Model and method

Planar CR3BP, Earth-Moon, `mu = 1.2150584270571545e-2`, same Jacobi-constant convention as the
exterior-MMR companion (no `mu(1-mu)` shift). **Perigee** Poincaré section (`l=0`), same choice as the
interior-MMR companion (Kumar et al. 2024 IAC) — this paper studies INTERIOR 2:1/3:1, unlike the
apogee-map exterior companion (AAS 25-569).

Formalizes the PIP/BIP/resonance-region vocabulary shared across this paper cluster in its most
careful, definitional form (Definitions 1-5, Figures 3-6): a periodic orbit's "orbit" under the
Poincaré map `P`; a primary intersection point (PIP); a boundary intersection point (BIP) defining a
local boundary between two Poincaré-section regions; the **stable resonance width** = semimajor-axis
span of the outermost closed (quasi-periodic librational) torus around a stable resonant center; the
**chaotic resonance-region width** = semimajor-axis span between the upper and lower BIPs of the
region's unstable-orbit separatrix (strictly larger than the stable width — the paper's central
distinction).

## Semi-analytical baseline (Section VI) — Gallardo's algorithm

Reviews the standard perturbed-Hamiltonian approach to MMR width (Gallardo 2006/2019/2020): expand the
disturbing function `R` in cosines of critical arguments `sigma = k_p*lambda_m - k*lambda + gamma`,
average, canonically transform to a single-harmonic normal form. Explicitly notes the caveats: assumes
fixed coplanar Moon orbit + a fixed nominal semimajor axis during the averaging integral (asteroid-case
justification: slow `(e,i,omega,Omega)` evolution vs. fast `(a,sigma)` oscillation) — an assumption the
paper's own results show breaks down for the highly-perturbed Earth-Moon case.

## Key quantitative findings — CR3BP vs. Gallardo widths

- **Gallardo's semi-analytical widths dramatically UNDERESTIMATE the true 2:1 and 3:1 resonance
  zones** (Figure 8 vs. Figures 9-11): the full CR3BP widths are consistently broader across
  `C in [2.50, 3.42]` (computed in `Delta C = 0.02` steps).
- The **2:1 resonance width does NOT taper toward `e=0`** the way Gallardo's approximation (and
  typical small-body/asteroid CR3BP resonance-width computations) predict; the CR3BP-computed 2:1
  width does not even center on the CR3BP's own 2:1 stable periodic orbit. The 2:1 zone does not show
  the "widest near `e~0.3-0.7`, tapering toward `e=0,1`" shape reported elsewhere in the small-body
  literature.
- **Combined 2:1 + 3:1 resonance regions span ~0.3 NDU (~115,000 km)** in semimajor axis (encompassing
  several higher-order sub-resonances in between) — a materially large fraction of the interior
  cislunar realm, again undermining a "resonances are narrow perturbations" mental model for Earth-
  Moon specifically (vs. e.g. the asteroid belt).
- The **2:1 and 3:1 chaotic-zone BIPs nearly touch** in the gap between them, i.e. they are right at
  the edge of overlapping — consistent with the interior-MMR companion's finding that 2:1<->3:1
  heteroclinics appear as soon as 2:1 unstable orbits exist at all (`C<=3.15`).
- **TLE cross-validation (Figure 12, real spacecraft)**: overlaying historic/current cataloged xGEO
  objects' two-line-element (TLE) time histories onto both Gallardo's widths and the CR3BP widths
  shows the **CR3BP widths correctly capture IBEX (3:1 resonant orbit) and TESS (2:1 resonant orbit)
  and Spektr-R (near-3:1)**, whereas Gallardo's narrower analytic bands do NOT — direct empirical
  validation that the fuller CR3BP treatment, not the classical semi-analytic one, is the right
  physical model for real Earth-Moon resonant spacecraft. This is a materially strong, sourced,
  real-world-validated result (not just an internal-consistency check).

## Heteroclinic connections and transfer TIMES (3:1 <-> 2:1) — the concrete numbers

- Heteroclinic connections between 3:1 and 2:1 unstable resonant orbits exist across
  `C in [~2.50, 3.15]` (matching the interior-MMR companion's finding, here with explicit transfer
  durations computed).
- **Transfer type 1 (short/direct)**: exists for `2.50 <= C <= 3.07`, transfer time **~28-29 days**.
  Ceases to exist for `C >= 3.09` (lower energy => reduced chaos => the direct stable/unstable
  manifold intersection needed for this type disappears).
- **Transfer type 2 (long, via an intermediate 5:2 resonance)**: exists at higher `C` (including
  `C=3.15`, where type 1 has already vanished), transfer time **~56-57 days** — almost exactly double
  type 1's duration, and via a genuinely different resonance-chain path (3:1 -> 5:2 -> 2:1) rather than
  a direct 3:1-2:1 heteroclinic.
- Transfer times are heuristic (defined by when the trajectory's averaged semimajor axis has visibly
  transitioned, since exact fixed-point convergence takes infinite time) but reported as robust across
  the `C` range within each type.
- By CR3BP time-reversal symmetry, a 2:1 -> 3:1 transfer of the same duration exists wherever a
  3:1 -> 2:1 transfer does.
- **4:1 confirmed isolated from 3:1 at every checked `C` (3.10, 2.85)** by a rotationally-invariant-
  circle (RIC) barrier — same finding as the interior-MMR companion, here with an explicit potential-
  RIC location marked on the `C=3.10` Poincaré map between the 3:1 and 4:1 islands.

## L1 tube interaction (Section IX, brief)

Notes a "strong interaction" between the 2:1 unstable periodic orbit's invariant manifolds and the
`L1` Lyapunov-orbit stable/unstable tube's first Poincaré cut, producing a visually "swirling" manifold
pattern from repeated Earth-realm exit / Moon-realm entry / Earth-realm re-emergence cycles (likened to
atomic-physics/chemistry transition-state numerical challenges, ref [34]). Flagged as future work, not
resolved quantitatively in this paper — the exterior-MMR companion (AAS 25-569 above) is the one that
actually quantifies `L1`-`L2` tube heteroclinic transfer times.

## Relevance to this project

- The **strongest single source in this 4-paper batch for `resonance_network.py`'s (#267 Track-B tier
  3) documented data gap**: that module's reproduce-before-trust gate is blocked on the exact
  arXiv:2509.12675/2026-ASR paper's PDF not being in the local mirror, specifically for "the exact
  common Jacobi constant... and the explicit form of the 'generalized distance metric'". This paper
  supplies concrete, sourced, independently-computed 3:1<->2:1 transfer-time numbers (28-29 days direct,
  56-57 days via 5:2) and Jacobi-constant ranges (`C<=3.07` direct-type-exists, `C>=3.09` direct-type-
  vanishes) from the SAME author group's methodologically-identical earlier work. **This is still not
  the same paper as the one the module cites** and should not be substituted for it without an
  explicit decision — but it is a strong, concrete, previously-missing candidate cross-check or
  stand-in source if the exact journal PDF continues to be unobtainable. Flagged as a follow-on
  opportunity for whoever next touches `resonance_network.py`'s xfail gate; not actioned in this pass.
- The Gallardo-vs-CR3BP TLE-validated "resonance widths are broader than semi-analytical prediction"
  finding is a citable, sourced caveat for any future MMR-width claim this project makes about the
  Earth-Moon system.

**Status: full mining pass complete.** No catalogue impact, no writeback (methods/results reference).
`CorpusAnchor` registered.
