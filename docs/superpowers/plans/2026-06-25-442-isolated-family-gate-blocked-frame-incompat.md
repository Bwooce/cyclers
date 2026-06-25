# #442 µ=0.001 isolated-family capability gate — BLOCKED (frame incompatibility, not digitization precision)

**Date:** 2026-06-25. The bounded capability gate from the #442 gap analysis:
can we digitize + converge ONE Antoniadou & Libert 2018 isolated elliptic
resonant family (the 3/1 MMR, µ=0.001) and confirm it is isolated? Verdict:
**BLOCKED-ON-DIGITIZATION — but the blocker is structural, not pixel error, and
the unblock path is concrete.**

## What was established (a thorough negative, ladder cleared for this approach)
- **Paper + figure found & readable.** `cyclers_pdf/papers/antoniadou-libert-2018-...-arxiv-1805.00288.pdf`. The 3/1 isolated family: Sec 5.3 p.17, config (π,0), "only stable periodic orbits, both bodies highly eccentric." Most transcribable representative: Fig 11(e) header `a1/a2=0.480674, e1=0.659951, ϖ1=M1=0°, ϖ2=M2=180°`, guiding `e2≈0.90`.
- **NO state vectors anywhere in the paper** (graphical DS-maps only) — the digest's "graphical only" is correct.
- **224-cell convention/sign grid** (both P1 apses, ϖ1∈{0,π}, both P2 apses; e1∈[0.55,0.95], e2∈[0.70,0.95], period_f=2π): NO member reached `independent_residual < 1e-8`. P2-pericentre seeds diverge (residual 1e2–1e5); the physically-correct P2-apocentre best case descends to crossing residual **6.0e-4** then the line search **stalls 5–7 orders short** of tol — the independent-closure gate is never even reached.

## The load-bearing finding — a FRAME mismatch
- **Paper frame (Eq.7, p.5):** NON-pulsating rotating frame; the perturber P2
  physically MOVES on the Ox-axis, its position `x'(0)=a2(1±e2)` IS the state
  coordinate encoding e2; the massless P1 is the `(x, ẏ)` perpendicular-crossing
  state; `e1` is P1's osculating eccentricity.
- **Our frame (`src/cyclerfinder/core/er3bp.py`):** Szebehely Ch.10 PULSATING
  (Nechvile) frame; BOTH primaries FIXED at (−µ,0),(1−µ,0); e enters ONLY via
  `scale = 1/(1+e·cos f)`. The particle's e1 is NOT a coordinate — it is implicit
  in the pulsating-frame Cartesian IC.
- These are **different frames.** The best osculating→pulsating map (verified to
  reduce EXACTLY to the validated e=0 seed recipe, max diff 0.00) omits the
  pulsating frame's `dr12/df` radial-velocity coupling at high e and cannot be
  validated without the paper's state vectors. **A tighter figure scan cannot fix
  this** — it is a modeling gap, not a digitization-precision gap.

## Unblock paths (precise)
1. **Build the corrector DIRECTLY in the paper's Eq.7 frame** (moving P2, `x'` as
   a coordinate / continuation variable), with `(e1,e2)` as the natural seed
   coordinates — sidesteps the conversion entirely; is what the paper actually
   uses. **Recommended.** Independently valuable: reproduces ALL Antoniadou
   isolated families, not just the EM target.
2. A faithful Eq.7-rotating → Szebehely-pulsating frame map INCLUDING the
   `dr12/df` velocity coupling, validated against a known orbit.
3. Author/data request to Antoniadou & Libert for the isolated-family state
   vectors (the paper publishes none) — cleanest, but external-dependency.

## Disposition
- The capability gate's TRUE COST is now established: not figure error, but a
  frame-incompatible corrector. The digitization approach is exhausted (224-cell
  ladder + validated circular limit); pursuing it further is wasted.
- **Continue via path 1** (corrector in the paper's frame), bounded by a
  checkable milestone: implement + validate the Eq.7-frame EOM against the e=0
  circular limit BEFORE attempting the isolated family.
- Escalation-ladder methodology recorded in task #446; gap-analysis context in
  registry entry `isolated-er3bp-cycler-novelty-gap-analysis-2026-06-25`.
- NB: this unblocks the µ=0.001 CAPABILITY gate; the downstream EM-cycler step
  remains SPECULATIVE (gap analysis: thin prior + semantic mismatch) with its own
  kill-criterion. Path 1 is justified by its standalone capability value, not the
  EM prize alone.
