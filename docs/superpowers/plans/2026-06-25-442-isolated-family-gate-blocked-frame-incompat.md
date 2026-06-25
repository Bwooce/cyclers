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

## UPDATE 2026-06-25 — Path 1 EOM VALIDATED (frame incompatibility UNBLOCKED)
Built the ER3BP EOM in Antoniadou's Eq.7 frame (`scripts/_scratch_442_paper_frame.py`,
`paper_frame_eom(t, state5, mu, e2)`) and ran the validation gate:
- **(a) e2=0 → CR3BP reduction: PASS, 8.16e-12** (paper-frame EOM reproduces our
  `cr3bp_eom` to machine precision at e2=0).
- **(b1) e2=0 Jacobi conservation: PASS, 3.36e-12.**
- **(b2) known-orbit closure (L4 equilibrium + a converged CR3BP DRO): PASS, 1.03e-14.**
- **OVERALL GATE: PASS** (runs in <1 s — the prior agent hangs were entirely in
  the high-e2 step-3 isolated-family hunt, NOT the EOM).

So the frame-incompatibility blocker is RESOLVED: we now have a validated ER3BP
corrector in the paper's OWN frame, where the published `(e1,e2)` configs are the
natural seed coordinates. Remaining: converge the actual 3/1 (π,0) isolated member
(`build_paper_ic` + `correct_paper_member`, period_f T=2π, gate on independent
residual < 1e-8). This is compute-expensive (e2≈0.90 is the near-collision chaotic
corner) and follows the #446 escalation ladder (apse/sign/θ0/e2-scan variants).
Next: promote the validated EOM to a proper module + test once a member converges
(or the convergence is honestly characterised).

## UPDATE 2026-06-25 (b) — convergence WALL characterised; answers the digitization question
Two corrections + a decisive convergence finding (two independent correctors):
- **PREMISE CORRECTION:** Antoniadou & Libert 2018 "Eq. 7" is the PERIODICITY
  CONDITION, not the EOM. The dynamics are the Lagrangian **Eq. 1**, integrated in
  **TIME t** (not true anomaly f), non-pulsating rotating frame, primaries sliding
  on Ox on their Keplerian ellipse. The validated EOM was built from Eq. 1 (the
  e2=0→CR3BP reduction confirms it term-by-term: Euler θ̈y / Coriolis 2θ̇v /
  centrifugal θ̇²x / two-moving-primary gravity).
- **3/1 (π,0) ISOLATED member does NOT converge.** Cold-shoot from the analytic
  graphical seed (a1=0.4807, e1=0.659951, e2≈0.90, P2 apocentre θ0=π) stalls at
  crossing residual ~0.3–0.5; the concurrent heavier LM/dogbox 13×13-grid corrector
  (`_scratch_442_converge.py`, same frame) reached crossing 2.3e-4 / **independent
  1.1e-2 — still ~6 orders short of the 1e-8 gate**. TWO independent correctors
  confirm: removing the frame incompatibility did NOT close the gate.
- **This IS the digitization-precision answer.** The binding constraint is now the
  graphical-only `(e1,e2)` seed precision into the isolated family's STEEP basin —
  and because the family is ISOLATED (no e=0 / circular limit by definition), there
  is **no continuation fallback** to rescue an imprecise seed. So for *this* family
  the digitization is demonstrably NOT good enough (confirmed by two solvers). The
  earlier "self-certifying via convergence" logic holds — it just certified the
  NEGATIVE: the seed never reached a converging basin.

## The next-step path (higher-probability than cold-shoot)
1. **Continuation from the e=0 CRTBP bifurcation** B^{3/1}_{I,1} (Scheme II, where the
   IC is EXACT), continued in e2 toward 0.90, freeing the period — this lands the
   CONNECTED 3/1 family (guaranteed capability win; reproduces a published family)
   and may give a homotopy bridge toward the isolated branch. NB the ISOLATED family
   has no e=0 limit, so continuation reaches the CONNECTED family; the isolated one
   needs either the homotopy-from-connected bridge or better seed data.
2. **Author/data request** to Antoniadou & Libert for the isolated-family state
   vectors — now JUSTIFIED (the paper publishes none; the graphical seed is proven
   insufficient for the isolated basin).
The validated EOM (Eq. 1 in the natural frame) is the durable capability deliverable
and is reusable for BOTH paths.
