# Digest: Kumar, Anderson, de la Llave, Gunter 2021 (AAS 21-651)

**Paper:** "Computation and Analysis of Jupiter-Europa and Jupiter-Ganymede Resonant Orbits in the
Planar Concentric Circular Restricted 4-Body Problem"
**Venue:** AAS/AIAA Astrodynamics Specialist Conference, AAS 21-651 (2021)
**Preprint:** arXiv:2109.14815v1 [math.DS] 30 Sep 2021
**Authors:** Bhanu Kumar (Georgia Tech), Rodney L. Anderson (JPL/Caltech), Rafael de la Llave
(Georgia Tech), Brian Gunter (Georgia Tech)
**Filed:** `kumar-anderson-delallave-gunter-2021-europa-ganymede-resonant-orbits-ccr4bp-AAS-21-651-arxiv-2109.14815.pdf`
**Acquired/digested:** 2026-07-23 (#688). Text-layer PDF, no OCR needed. md5 distinct from the
adjacent-arXiv-ID SIADS paper 2109.14814 (already in corpus) — genuinely different paper, not a dup.

**Why in corpus (#688 / #686 Stage B):** this is the direct predecessor and model definition for
the Jupiter-Europa-Ganymede CCR4BP that `#686`'s shortlist item (the Laplace-locked repeating
resonant tour) proposes to search. It is the acquisition the `#688` dispatch flagged as PC-support
so Stage B can start immediately on GO. Companion papers: SIADS 24(1):219-258 / arXiv:2109.14814
(the parameterization/GPU-connection machinery, PC1, already digested
`2026-07-03-digest-kumar-anderson-delallave-2025-whiskered-tori-connections.md`) and AAS 23-397 /
arXiv:2309.06073 (secondary-resonance overlap, digested alongside this one, #688).

## The model (planar CCR4BP)

- m1 (Jupiter) + m2 (Europa) + m3 (Ganymede) + massless spacecraft. m2 and m3 revolve in
  CONCENTRIC CIRCLES of radii r12, r13 about m1 at angular velocities Ω2, Ω3, with NO mutual effect
  (m2 does not perturb m3 or vice versa) — this is what makes the model tractable.
- Mass ratios µ = m2/(m1+m2), µ3 = m3/(m1+m3). Physical Galilean masses used.
- KEY STRUCTURE (the tractability hinge, and the structural test Stage B must implement): the
  CCR4BP is simultaneously a TIME-PERIODIC PERTURBATION of BOTH the m1-m2 (Jupiter-Europa) CRTBP
  AND the m1-m3 (Jupiter-Ganymede) CRTBP. At µ3 = 0 it reduces exactly to the Jupiter-Europa
  PCRTBP; at µ̄2 = 0 (Europa mass → 0 in the Ganymede frame) it reduces to the Jupiter-Ganymede
  PCRTBP. (This is precisely the "reduces to each PCRTBP at the other µ→0" structural test the #686
  strategy note prescribes for the Stage-B EOM module, mirroring bcr4bp's µ_S→0 CR3BP test.)
- The 2:1 Laplace commensurability (Europa ≈ 2 revs per 1 Ganymede rev) is what makes the joint
  forcing TIME-PERIODIC and lets the stroboscopic reduction apply. (Our #688 screen independently
  measures the physical period ratio at 2.014 from JPL SMAs.)

## Method

- Unstable PCRTBP resonant periodic orbits persist under the periodic forcing as 2D UNSTABLE
  INVARIANT TORI (quasi-periodic orbits) in the 5D extended phase space (x,y,px,py,θ3). Not periodic
  orbits — tori.
- STROBOSCOPIC MAP F drops 5D→4D and the torus 2D→1D invariant circle K(θ) with F(K(θ)) = K(θ+ω),
  rotation number ω = 2πΩ1/|Ωp−1|.
- Predictor: Poincaré-Lindstedt series. Corrector: the O(N log N) FFT-based quasi-Newton
  parameterization method of their prior work (= the 2109.14814 machinery), which returns the torus
  K(θ) AND the bundle matrix P(θ) whose columns 3 & 4 ARE the linear stable/unstable manifold
  directions with their Floquet multipliers — i.e. whiskers come out of the solve, not via post-hoc
  STM extraction (the structural fix for the #619 failure mechanism).
- Tori are computed in a Jupiter-Europa synodic frame OR a Jupiter-Ganymede rotating frame, related
  by an explicit transformation Φ_θ3; torus SHAPES differ significantly between frames.

## Results (positive controls for Stage B)

1. Successfully computed CCR4BP quasi-periodic analogues of the Jupiter-Europa **3:4** and
   Jupiter-Ganymede **3:2 and 7:5** resonant periodic orbits, with physical Galilean masses.
2. Continuation of the Jupiter-Europa 3:4 torus with Europa mass → 0 → the unstable Jupiter-Europa
   resonant torus "likely approaches a stable Jupiter-Ganymede PCRTBP nonresonant KAM torus lying
   near but not on the 3:2 Jupiter-Ganymede resonance." (One-way linkage, at torus level.)
3. Verified the quasi-Newton corrector on a NEW model (vs the elliptic RTBP of the prior paper) —
   performance "very positive."

## What is explicitly UNCLAIMED (the Stage-B window)

The conclusion names the next step as future work, NOT done here: "searching for heteroclinic
connections in the CCR4BP linking Jupiter-Ganymede resonances to Jupiter-Europa resonant orbits."
The CLOSED, repeating cycle (the cycler-class object) is a further step beyond even that one-way
connection. This is the object `#686`'s shortlist targets and the reason the unclaimed-cycle window
is real (but at scoop risk — these authors own the machinery).

## Bearing on #688's screen result

This paper confirms WHY the cheap exterior-Keplerian-map screen (#688 Stage A) cannot itself find
the tour: the physically relevant Jupiter-Ganymede resonances are INTERIOR (3:2 → a=0.76, 4:3 →
a=0.83 in Ganymede units), the connection objects are torus stable/unstable MANIFOLDS in a 4D
stroboscopic map (not periapsis kicks), and the Europa↔Ganymede link is a heteroclinic transfer,
not a periapsis-pinned migration. The #688 seed geometry (Europa 3:4 ↔ Ganymede 3:2 orbits
radially overlap with a ~30 m/s coplanar speed-match) independently corroborates, at screen
resolution, that this paper's exact target resonance pair is energetically compatible.
