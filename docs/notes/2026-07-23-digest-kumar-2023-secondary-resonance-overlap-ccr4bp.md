# Digest: Kumar, Anderson, de la Llave 2023 (AAS 23-397)

**Paper:** "4th Body-Induced Secondary Resonance Overlapping Inside Unstable Resonant Orbit
Families: A Jupiter-Ganymede 4:3 + Europa Case Study"
**Venue:** AAS/AIAA Astrodynamics Specialist Conference, AAS 23-397 (2023)
**Preprint:** arXiv:2309.06073v2 [astro-ph.EP] 16 Sep 2023
**Authors:** Bhanu Kumar (Georgia Tech), Rodney L. Anderson (JPL/Caltech), Rafael de la Llave
(Georgia Tech)
**Filed:** `kumar-anderson-delallave-2023-secondary-resonance-overlap-ganymede-4-3-ccr4bp-AAS-23-397-arxiv-2309.06073.pdf`
**Acquired/digested:** 2026-07-23 (#688). Text-layer PDF, no OCR needed.

**Why in corpus (#688 / #686 Stage B):** this is the paper `#686`'s shortlist item cites for the
claim that the CCR4BP generates GENUINELY 4-body-native objects with NO 3-body analog — the core
"is it really 4-body-native, not a perturbed 3-body orbit" bar the discovery-strategy pass had to
clear. Acquired in the same pass as AAS 21-651 (arXiv 2109.14815) per the `#688` dispatch.

## Setup

- Same Jupiter-Europa-Ganymede planar CCR4BP as AAS 21-651 (concentric circular moons, mass ratios
  µ = m2/(m1+m2), µ3 = m3/(m1+m3); Hamiltonian H_µ3(x,y,px,py,θ3); at µ3=0 reduces to the
  Jupiter-Ganymede PCRTBP).
- Object of study: the family of UNSTABLE Jupiter-Ganymede **4:3** resonant orbits (interior,
  a ≈ 0.83 in Ganymede units), whose stable/unstable manifolds are the most useful for low-time-of-
  flight resonance-hopping transfers.

## Main result (the 4-body-native mechanism)

- Europa's forcing generates SECONDARY RESONANCES between the internal libration frequency of the
  Jupiter-Ganymede 4:3 unstable orbits and Europa's forcing frequency. Despite being high-order,
  four are significant: **11/34, 12/37, 23/71, 25/77**.
- These secondary resonances undergo CHIRIKOV OVERLAP (Europa's forcing m3... i.e. the 4th body is
  strong enough), which is directly computed and confirmed by generating the new orbits living
  inside the secondary-resonance islands.
- The overlap causes "a complete structural change of the higher-energy unstable 4:3 orbits whose
  manifolds are most useful for low-TOF orbit transfers."
- These secondary-resonance objects are ENTIRELY NEW orbit types created by the 4th body — no
  analog exists in any single CRTBP. The authors argue the phenomenon is general (major
  implications for resonant-orbit use in multi-moon tour design).

## Bearing on #686 / #688

- Establishes (grounds `#686`'s §3 claim (2)) that the CCR4BP is 4-body-NATIVE: the secondary
  resonances and their overlap-generated orbits cannot be found by perturbing any 3-body cycler,
  because no single CRTBP contains an object whose frequency locks to a second moon's forcing. This
  is a genuine "no 3-body analog" structure, satisfying the discovery-pass's core bar.
- Reinforces #688's screen-negative interpretation: the transport-relevant Ganymede 4:3 objects are
  interior AND structurally reorganized by Europa's forcing — a regime entirely outside the RS07
  exterior-periapsis Keplerian map's validity. The cheap two-map screen cannot represent them; only
  the full Stage-B CCR4BP + parameterization-method torus/manifold build can.
- CAUTION for Stage B (scoop/complexity risk): near the overlapping secondary resonances the 4:3
  family's structure changes qualitatively, so a naive continuation through that energy band may
  hit the structural change rather than a clean torus. Stage B should measure where the overlap
  band sits before choosing target energies.
