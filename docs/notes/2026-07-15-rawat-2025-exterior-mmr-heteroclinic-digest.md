# Rawat-Kumar-Rosengren-Ross 2025 — Exterior MMR regions of influence, Earth-Moon (AAS 25-569)

**Date filed/digested:** 2026-07-15 (light digest); **full mining pass:** 2026-07-15 (#597/#15).
**Trigger:** discovered via a manual review of `https://ross.aoe.vt.edu/papers/`, requested by the
user during the #595 (#498/#499/#503 acquisition-status) staleness audit. Confirmed genuinely new by
reading the actual title page; full 20-page body read page-by-page for this pass.

**Citation:**
> Rawat, A., Kumar, B., Rosengren, A. J., and Ross, S. D. (2025). "Regions of Influence of Exterior
> Mean-Motion Resonances in the Earth–Moon System: Bifurcations, Separatrices, and Heteroclinic
> Pathways." AAS/AIAA Astrodynamics Specialist Conference, AAS 25-569.
> No DOI found (AAS conference paper, not indexed in CrossRef as of this filing).

Registered as `CorpusAnchor` `rawat-2025-exterior-mmr-aas` in
`src/cyclerfinder/search/literature_check.py` (provenance `verified-against-source`).

## Model and method

Planar CR3BP, Earth-Moon, `mu = 1.2150584270571545e-2`. Uses a Jacobi-constant convention WITHOUT the
`mu(1-mu)` additive shift some authors use (so `C(L4/L5) != 3` here) — flagged explicitly as a
footnote; matters if ever cross-comparing `C` values against a source using the other convention.

**Poincaré section: apogee** (osculating mean/true anomaly `l = pi`), not perigee — the opposite
choice from the interior-MMR companion (Kumar et al. 2024 above) and the resonance-widths companion
(Rawat et al. 2024 AAS 24-368 below), both of which use perigee. Rationale given: in the EXTERIOR
region, the spacecraft is least perturbed by the Moon at apogee, giving a cleaner resonant-dynamics
picture there (the reverse of the interior case, where perigee is the "quiet" point). Section
coordinates: `(varpi, a)` where `varpi` is the **synodic longitude of perigee** (angle between perigee
and the Earth-Moon line in the rotating frame) — a cylindrical topology `(varpi, a) in S^1 x I`.

**Resonance naming convention differs from the interior case**: exterior resonances are `k:km` with
`k < km` (spacecraft period `T` longer than Moon's `T_m`; `T/T_m ~ km/k`). A `1:n` exterior resonance
is always a period-1 fixed point of the apogee map (vs. period-`m` for interior `m:n`).

**PIP/BIP/resonance-zone formalism** (shared vocabulary with AAS 24-368 below, presented here in more
detail): a primary intersection point (PIP) is an intersection of stable/unstable manifold branches of
a saddle-type periodic orbit that is "closest" in arc-length; a boundary intersection point (BIP) is
the PIP chosen (shortest arc-length) to define a region's boundary; the resonance zone/region of
influence is the semimajor-axis span between the top and bottom BIPs.

## Key findings — 1:3 and 1:4 exterior resonances

- Unlike interior resonances, purely exterior 1:n resonant orbits do **not** occur as symmetric
  stable/unstable pairs. Instead: **two asymmetric stable islands** (resonant centers away from
  `varpi=0` or `pi`) are born via bifurcation from **symmetric weak unstable** orbits, all embedded
  within the separatrices of a **symmetric strong unstable** orbit. This 3-tier structure (strong
  unstable envelope / weak unstable sub-structure / two asymmetric stable islands) is the paper's
  headline topological finding, distinct from the interior case's simple stable-center + unstable-
  saddle pairing.
- Both the **strong** unstable family (three separate sub-families across differing eccentricity:
  low-e "Family 1" purely exterior, mid-e "Family 2" entering the Earth-Moon gap without entering
  interior, high-e "Family 3" entering interior realm too) and the **weak** unstable family (a single
  continuous family for both 1:3 and 1:4) were computed. 1:3 weak unstable ranges `C in
  [-0.56, 3.37]`; 1:4 weak unstable ranges `C in [-0.35, 3.57]`.
- **Combined 1:3 + 1:4 resonance zones span ~735,511 km (~2 NDU) in semimajor axis** — i.e. these two
  exterior resonances alone dominate a huge fraction of the exterior realm.
- The **weak-unstable quasi-separatrices are always contained inside the strong-unstable separatrix**
  (region of influence of strong unstable is always the outer/dominant one), but for 1:4 specifically
  the weak orbit's quasi-separatrices are **multi-lobed with turnstiles that intersect the strong
  orbit's separatrices at multiple points**, producing complex heteroclinic-like interactions WITHIN
  a single MMR (not just between different MMRs) — a richer structure than 1:3, whose weak
  quasi-separatrices stay single-lobed and distinct from the strong ones (limited internal chaotic
  transport).
- Resonance-zone width/shape is markedly `C`-dependent: lower `C` (higher energy) -> larger, more
  intersecting lobes; higher `C` -> smaller islands, less internal mixing. A `C=3.10..3.15` "1/4-type"
  bifurcation of the **1:2** stable island is separately noted (Greene's-criterion-consistent
  invariant-torus breakup signature), distinct from the 1:3/1:4 main results.
- Standard symmetry-based differential correction is inadequate for the asymmetric stable orbits
  ("the second author very recently developed a suitable methodology" — cited as their own separate,
  in-progress work, ref [21], not detailed here).

## Heteroclinic connections: exterior 1:3 <-> interior 2:1/3:1 via L1-L2 tubes

- Builds on the tube-dynamics machinery of Koon et al. 2000/2001/2003: `L1` and `L2` Lyapunov-orbit
  stable/unstable manifolds each have an interior-realm branch and an exterior/exterior-realm branch;
  intersections of these tubes with resonant-orbit manifold "lobes" enable free transfers.
  Constructs a Poincaré section at `x = 1-mu` (the Moon's location) to track the first crossing of the
  `L2` unstable manifold (forward) and `L1` stable manifold (backward), finding **substantial overlap
  between these tubes at `C=3.10`** — a wide window facilitating exterior<->interior transit via the
  `L1`/`L2` gateways.
- Concretely finds heteroclinic connections between the **1:3 exterior** (both strong and weak
  unstable variants) and the **2:1 interior** unstable resonant orbit, at `C=3.10`: transfer time
  **~9.5 days** (1:3 strong -> 2:1) and **~7 days** (1:3 weak -> 2:1). Similar connections found
  between 1:3 and **3:1 interior** at `C=3.10` (two heteroclinic points for both strong/weak variants).
- **At `C=3.15` the `L1`-`L2` overlap window shrinks drastically**: no heteroclinic connection is found
  between 1:3 (either variant) and 2:1, and the same negative result holds for 1:3-to-3:1 at `C=3.15`.
  This establishes a fairly sharp transition between "interior-exterior transit is easy" (`C~3.10`) and
  "not really open" (`C~3.15`) within a narrow `Delta C ~ 0.05` window.
- Time-reversal symmetry (an explicit modeling assumption of the planar CR3BP equations) means every
  found forward transfer implies the reverse (interior -> exterior) transfer exists at the same `C`.

## Relevance to this project

- Same relevance profile as the interior-MMR companion digest (`kumar-2024-interior-mmr-iac`): direct
  precursor-group literature for #314/#411/#496/#503, and independent (though not identical) supporting
  evidence for `resonance_network.py`'s (#267 Track-B tier 3) documented reproduce-before-trust data
  gap against the exact 2026 ASR/arXiv:2509.12675 paper — this AAS 25-569 conference paper supplies the
  concrete `L1`-`L2` tube transfer-time numbers (7-9.5 days) and the sharp `C=3.10` vs `C=3.15`
  interior-exterior transit threshold that the journal paper's un-mirrored PDF currently cannot supply.
  Flagged as a follow-on opportunity, not actioned here.
- The exterior 1:n asymmetric-bifurcation topology (strong/weak unstable + 2 asymmetric stable
  islands) is a genuinely distinct structural finding from anything currently modeled in this
  project's own genome/search machinery — worth keeping in mind if a future task targets exterior
  cislunar resonant orbits specifically (none currently do).

**Status: full mining pass complete.** No catalogue impact, no writeback (methods/results reference).
`CorpusAnchor` registered.
