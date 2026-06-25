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

## UPDATE 2026-06-25 (c) — CONNECTED 3/1 (π,0) family CONVERGED e2=0→0.90 (path 1 executed)
Ran the full #446 ladder with trust-region/LM correctors + e2-continuation
(`scripts/_scratch_442_converge.py`, `_scratch_442_converge2.py`,
`_scratch_442_arclength.py`). Decisive results:
- **The 6e-4 / 1e-2 "wall" was a SOLVER artifact, not dynamics.** Replacing the
  2-var half-period (y,vx) damped-Newton with a **4-vector full-period periodicity
  objective** `[x(2π)−x, y(2π), vx(2π), vy(2π)−vy]` under `scipy.least_squares`
  (lm/trf/dogbox), the family converges CLEANLY. From a geometric apse seed at LOW
  e2 (where the basin is wide) the corrector hits indep ~1e-11 immediately.
- **CONNECTED 3/1 (π,0) family fully tracked e2=0 → 0.90**, every member
  independently re-integrated at rtol=1e-13 closing < 2e-10 (gate 1e-8):
  - e2=0.00: x=0.4793554969, vy=0.9624369118, indep=3.1e-12, a1=0.48014, e1≈0
  - e2=0.30: x=0.6909309271, vy=0.5087470978, indep=5.1e-11, a1=0.48062, e1=0.440
  - e2=0.55: x=0.8217405562, vy=0.3046047985, indep=5.8e-11, a1=0.48065, e1=0.713
  - e2=0.85: x=0.9457500071, vy=0.0255730480, indep=3.0e-11, a1=0.48047, e1=0.972
  - e2=0.90: x=0.9341012757, vy=−0.3487825823, indep=1.9e-10, a1=0.48054, e1=0.948
  (IC = [x, y=0, vx=0, vy, θ0=π], T=2π.) osc a1≈0.4806 throughout = the published
  3/1 DS-map value (Fig 3 caption 0.4807). Adaptive natural-parameter AND
  pseudo-arclength continuation agree; arclength reached e2=0.91.
- **It is CONNECTED, not isolated.** Reverse-continuation reaches e2=0 cleanly
  (indep 3e-12, a genuine circular CR3BP 3/1 PO at x=0.4794, vy=0.9624) — i.e. this
  is the Scheme-II BIFURCATING branch, exactly as the digest predicted. Its e1
  GROWS monotonically with e2 (0→0.97), so at e2=0.90 it sits at e1=0.948, NOT the
  paper's ISOLATED representative e1=0.659951.
- **ISOLATED member NOT reached (honest negative, well-characterised).** A coarse
  (x,vy) periodicity-residual map over the whole section at e2=0.90, θ0=π shows the
  ONLY clean a1≈0.48 doubly-symmetric PO is the connected one (3 distinct basins all
  polish to the same x=0.9341, e1=0.948 member). The geometric isolated seed
  (e1=0.659951, M1=0 pericentre → x=0.16, near-collision) sits at residual 0.24 with
  NO low-residual basin within ±0.1×±0.4. The complementary θ0=0 section is
  unphysical (P2 pericentre r=0.1 → all orbits blow up, resid ~3e3). So the isolated
  family is genuinely disconnected and its representative is not a simple
  (x,0,0,vy)|θ0=π perpendicular crossing reachable from the graphical seed — it needs
  the family's own state vectors or a homotopy bridge that this section does not expose.
- **DELIVERABLE:** path 1 produced a reproducible, published-value-matching family
  (capability win). The isolated-family gate remains seed-data-limited exactly as the
  digest §7.2 forecast (direct high-e seeding needs the paper's ICs, which it omits).

## UPDATE 2026-06-25 (d) — EXHAUSTIVE data-availability check (author email VETOED by user)
User vetoed contacting the authors and asked whether the ICs exist anywhere / in software.
Checked every realistic source — the answer is a DEFINITIVE NEGATIVE on data availability:
- **arXiv:1805.00288 ancillary files:** NONE (only PDF + TeX source).
- **arXiv TeX source (downloaded + inspected, /tmp):** `arxiv.tex` + 59 pre-rendered PDF
  figures + `.bbl`. `\begin{tabular}` count = 0; no `.dat`/`.csv`/`.txt` data files. The
  only numeric ICs in the prose are three CIRCULAR-family bifurcation x-values
  (3/2 x≈0.763143, 5/2 x≈0.542884, 4/1 x≈0.39685 — and NONE printed for 3/1) plus
  figure-caption (e1,e2) pairs. NO state-vector tables, NO isolated-family ICs.
- **2019 spatial sister paper:** our digest is explicit — "no IC table anywhere"; Table A1
  is a structural map, not state vectors.
- **Web (thesis / repository / periodic-orbit database / group code):** none found. The
  Voyatzis–Antoniadou (Thessaloniki) school publishes families as graphical characteristic
  curves / DS-maps as a consistent NORM — state vectors are not tabulated anywhere.
- **JPL three-body catalogue:** no match (checked during #444).
CONCLUSION: the isolated-family ICs are genuinely unpublished and unretrievable without the
authors. With email vetoed, the isolated prize is now purely METHOD-gated, not data-gated.

**Bonus — the source TEXT clarifies the 3/1 mechanism (corrects this doc's "no e=0 limit"
framing).** §"3/1 MMR" of arxiv.tex: the circular family is UNSTABLE through the 3/1
neighbourhood; at the two ENDINGS of that unstable segment sit critical orbits (double
eigenvalue −1) from which branches bifurcate at T=2T0=2π (Scheme II): branch **I (stable)**
and branch **II (unstable)**. The published ISOLATED stable family ("Ic stable for
0.75<e1<0.98") is the HIGH-e STABLE SEGMENT OF BRANCH II. #442's converger reached branch I
(stable, low-e, e1 grows monotonically with e2). The isolated prize is therefore reachable
by continuing **branch II** (the unstable circular-bifurcation branch) in e1 to its high-e
stable segment — a method path (er3bp_branching.py + #437 fold-aware pseudo-arclength),
needing NO author data. Recorded as the concrete next-method (was task #455's prize; #455
the EMAIL approach is retired). This is the live, ours-to-build route to the isolated family.

## UPDATE 2026-06-25 (e) — ISOLATED FAMILY SOLVED (#457), method-only, no author data
The UPDATE (d) branch-II method path was EXECUTED in #457 and SUCCEEDED. The published 3/1
isolated stable I_c representative is now REACHED + CLOSED without author data:
- IC (paper-frame): x=0.16151871386593838, y=0, vx=0, vy=3.166889839559741, θ0=π, e2=0.91, T=2π.
- Independent full-period residual 4.3e-12 (gate 1e-8), integrator-invariant (DOP853 1e-12/1e-13/3e-14 + Radau); doubly-symmetric (perpendicular at t=0 AND t=π); STABLE (monodromy on unit circle).
- a1=0.480353 vs published 0.480674 (dev 3e-4); e1=0.659774 vs published 0.659951 (dev 1.8e-4) — the small offsets are because the exact member sits at e2=0.91 vs the figure header's nominal e2≈0.90.
- KEY FIX (broke the #442/#457-interim wall): the plain 2-var full-period corrector floors at ~2e-6
  on a non-doubly-symmetric compromise point; a JOINT doubly-symmetric corrector enforcing
  perpendicularity at BOTH crossings (residual [y(T/2),vx(T/2),y(T),vx(T)]) selects the true member.
- DELIVERABLES: src/cyclerfinder/core/er3bp_paper_frame.py::correct_doubly_symmetric_member +
  osculating_e1; golden test tests/core/test_er3bp_paper_frame.py::test_isolated_ic_stable_member_reached_and_closed
  (sourced expected values). Verified independently: 8/8 tests pass, ruff/mypy clean.
- DISPOSITION: this is a CAPABILITY win + a REPRODUCTION of a published family — it is a µ=0.001
  planetary resonant PO, NOT a cycler and NOT Earth-Moon, so NOT a catalogue row. The original
  #442 "novel Earth-Moon CYCLER prize" remains unachieved, but the capability that blocked it
  (reaching isolated elliptic resonant families) is now SOLVED. So UPDATE (d)'s "author-data-gated"
  is CORRECTED: the isolated family was method-gated, and the method now exists. 4/1 + 5/1 isolated
  members are the natural extension (follow-up task).
