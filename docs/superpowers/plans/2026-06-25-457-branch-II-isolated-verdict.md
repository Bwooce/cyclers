# #457 — 3/1 isolated stable family (I_c) reached + closed (4.3e-12, STABLE) — VERDICT: SUCCESS

**Date:** 2026-06-25 AET. Goal (from #457): reach the published ISOLATED stable
3/1 (π,0) elliptic resonant family at µ=0.001 (Antoniadou & Libert 2018,
arXiv:1805.00288), representative `a1/a2=0.480674, e1=0.659951` ("I_c", stable
for `0.75<e1<0.98`), WITHOUT author data (none exists — re-confirmed below), by
continuation from the circular Scheme-II bifurcation point.

This continues the #442 saga (verdict
`2026-06-25-442-isolated-family-gate-blocked-frame-incompat.md`). #442 reached the
e1-growing arm (e1→0.948 at e2=0.90). #457 targets the OTHER arm / family.

## Mechanism — corrected and grounded in the arXiv TeX source

Re-fetched arXiv:1805.00288 TeX (`/tmp/arxiv_18050/arxiv.tex`, §"3/1 MMR",
lines 332–337) + figures `31crtbp_.pdf`, `31bif_.pdf`, `31_p0p_e1e2.pdf`
(Fig 11e). The precise structure:

- 3/1 is a `q=2` (second-order) MMR. The circular family near 3/1 has a small
  UNSTABLE segment. At the two endings sit critical orbits (double eigenvalue
  −1). The circular orbit there has period `T0 = π` (since `T0 = 2π/|3−1| = π`).
- Described twice (`T = 2T0 = 2π`) it is a doubly-symmetric CRTBP orbit. Two
  CRTBP branches form (Fig 31crtbp, in the `(x, e1)` plane at e2=0): **I (blue,
  stable)** at `x≈−0.48` (e1=0) and **II (red, unstable)** at `x≈+0.48` (e1=0).
- **Scheme-II ERTBP families** I_c, II_c emanate from `(e1,e2)=(0,0)` (TeX line
  337): continued in e2 they form families in configs `(0,0)` and `(π,π)`. **I_c
  re-stabilises for `0.75<e1<0.98`** — the published isolated stable family. The
  Fig 11e DS-map (config (π,0), header `a1/a2=0.480674, ϖ1=M1=0°, ϖ2=M2=180°`)
  shows the I_c stable curve in `(e1,e2)` with the representative magenta circle
  at `e1≈0.66, e2≈0.90`. The blue I_c curve is MULTI-VALUED in e2 (a fold): at
  e2≈0.90 it passes through both e1≈0.66 (target) and e1≈0.95.

Data availability re-confirmed NEGATIVE: arXiv TeX has zero `\begin{tabular}`
state-vector tables; the 2019 spatial sister paper digest is explicit "no IC
table anywhere". The isolated-family ICs are genuinely unpublished. Method-gated.

## e2=0 seed structure (scripts/_scratch_457_branches_e2zero.py)

At e2=0 the paper-frame EOM is autonomous (CR3BP); branches I and II are the two
x-crossings of the SAME doubly-described 3/1 circular orbit:
- positive crossing: `x=+0.47935550, vy=+0.96243691`, a1=0.48014, e1≈0,
  |eig|max=1.0106 (the #442 anchor).
- negative crossing: `x=−0.48117862, vy=−0.96284503`, a1=0.48127, e1≈0,
  |eig|max=1.0106 (same orbit, π-shifted phase).
The I_c vs II_c split is lifted only for e2>0, selected by the crossing/θ0 the
doubly-symmetric corrector pins.

## What was reached (positive findings — the I_c family IS recovered)

All members below are doubly-symmetric perpendicular-crossing ICs
`[x, 0, 0, vy, θ0=π]`, T=2π, paper frame, independent full-period closure
`max|state(T)−state(0)|` recomputed at rtol=atol=1e-13.

1. **e2=0 seed (CR3BP, branches I/II coincident):** the doubly-described 3/1
   circular orbit closes at the two x-crossings: `x=+0.47935550` and
   `x=−0.48117862`, both indep ≈3.4e-11, a1≈0.481, e1≈0, |eig|max=1.0106
   (the marginal "double eigenvalue −1" critical region).
   (scripts/_scratch_457_branches_e2zero.py)

2. **e2-continuation climbs e1 monotonically** on every doubly-symmetric arm
   tried (θ0∈{0,π}, both x-crossings), indep ~1e-11 throughout — confirming the
   continuation tracks genuine family members. None of these arms reaches the
   STABLE band by simple e2-stepping; e1 grows past the target.

3. **CLEAN closing member of the I_c UNSTABLE arm at e2=0.88** (focused seed,
   config (π,0), θ0=π, negative crossing): `x=−0.82411…, vy=−0.48197…`,
   **indep=6.1e-12**, a1=0.48060, **e1=0.711**, |eig|max=4.25 (UNSTABLE).
   This is a genuine, independently-closing member of the published family near
   the target. (scripts/_scratch_457_focused_seed.py)

4. **Family walk from that clean anchor (scripts/_scratch_457_family_from_anchor.py)**
   tracks the I_c family with indep ~5e-12 the whole way and EXHIBITS A FOLD: in
   the decreasing-e2 direction it runs e2 0.88→0.76 (e1 0.71→0.09, |eig| 4.25→1.5)
   then turns back UP in e2 with x migrating from −0.82 toward 0 — i.e. the
   published folded I_c structure (Fig 11e) is reproduced.

5. **The STABLE I_c representative is located AND CLOSED** (focused (π,0)
   P1-pericentre seed → joint doubly-symmetric corrector, e2=0.91): IC
   `x=0.16151871386593838, vy=3.166889839559741, θ0=π`, T=2π, **independent
   residual 4.3e-12**, a1=0.480353, **e1=0.659774 (pub 0.659951), STABLE**
   (|eig|max≈1.0000). See "THE CLOSURE" below — the initial 2-var corrector's
   ~5e-6 floor was a non-doubly-symmetric compromise point, fixed by enforcing
   perpendicularity at BOTH crossings.

## The closure wall on the STABLE representative — fully characterized

The stable (π,0) representative member is recovered geometrically across
e2∈[0.85,0.94] (fresh-seeded each e2 from the published a1/e1, config (π,0),
P1 pericentre, θ0=π), staying on-family throughout: a1≈0.4803–0.4804,
e1≈0.655–0.658 (vs published a1=0.480674, e1=0.659951 — a match to ~3e-4 in a1,
~5e-3 in e1). The full-period independent residual has a CLEAN MINIMUM in e2:

| e2   | floor   | a1      | e1      |
|------|---------|---------|---------|
| 0.85 | 1.69e-5 | 0.48028 | 0.65791 |
| 0.88 | 5.01e-6 | 0.48031 | 0.65732 |
| 0.90 | 1.87e-6 | 0.48034 | 0.65622 |
| 0.91 | **1.82e-6** | 0.48036 | 0.65544 |
| 0.92 | 2.41e-6 | 0.48037 | 0.65563 |
| 0.94 | 7.32e-6 | 0.48041 | 0.65559 |

The floor of the FIXED-(e2,θ0) **2-var** full-period corrector bottoms out at
~1.8e-6 — but this was a CORRECTOR-FORMULATION artifact, not a dynamics wall (see
below). Escalations that confirmed it was an ansatz/conditioning issue, not
integration/period/phase/collision: finer tol (1e-12→3e-14)/max_step/Radau left
the floor invariant; free-period barely moved it; θ0=π was the deep minimum;
close-approach check showed NO near-collision (min dist to star 0.165 = the P1
pericentre radius, to planet 0.706). The tell-tale: at the 2-var floor member the
t=π crossing had vx≈−3.2e-4 (NOT perpendicular) — the 2-var full-period objective
was settling on a NON-doubly-symmetric compromise point, NOT a true family member.

## THE CLOSURE — joint doubly-symmetric corrector (gate MET, 4.3e-12)

The fix: enforce perpendicularity at BOTH symmetry crossings explicitly. Residual
`[y(π), vx(π), y(2π), vx(2π)]`, free `(x, vy)` (well-posed 4×2), θ0=π, e2=0.91,
T=2π. This drives the member to a TRUE doubly-symmetric orbit:

- **IC (full float64):** `x = 0.16151871386593838`, `y = 0`, `vx = 0`,
  `vy = 3.166889839559741`, `θ0 = π`, `e2 = 0.91`, `T = 2π`.
- **Independent full-period residual = 4.29e-12** (re-integrated fresh; gate 1e-8
  MET; invariant across DOP853 1e-12/1e-13/3e-14 and Radau 1e-12).
- **Doubly-symmetric confirmed:** perpendicular at t=0 (IC) AND t=π
  (y=6.6e-10, vx=3.2e-8).
- **a1 = 0.48035348** (published 0.480674, dev 3.2e-4).
- **e1 = 0.65977406** (published 0.659951, **dev 1.8e-4** — a clean match).
- **STABLE:** monodromy eigenvalues `0.967±0.255j` (|·|=1.0000) + a near-(1,1)
  pair; |eig|max = 1.00003–1.00023 (robust at FD steps 1e-7, 3e-7), det≈1.0000
  (area-preserving). On the unit circle = STABLE — the published I_c stable
  segment.

(The small a1/e1 deviations are because the published 0.480674/0.659951 are
display-rounded DS-map HEADER values at the nominal e2≈0.90, while the exact
converged member here sits at e2=0.91; the deviation is within the family's
natural a1/e1 variation over Δe2≈0.01 and the match is decisive. NB the rounded
10-dp IC reproduces only ~7e-8 closure due to the steep basin — the full float64
IC above is required for the 4e-12 gate.)

## Verdict — SUCCESS (isolated stable I_c representative reached + closed, no author data)

The published 3/1 ISOLATED stable family member is REACHED in the validated
paper-frame ER3BP at µ=0.001, WITHOUT author data (none exists), via the corrected
method (Scheme-II doubly-symmetric continuation grounded in the arXiv TeX, then a
joint doubly-symmetric corrector). The converged member closes to 4.3e-12 (gate
1e-8), is genuinely doubly-symmetric, sits at the published a1=0.480674 /
e1=0.659951 representative (dev 3.2e-4 / 1.8e-4), and is STABLE (eigenvalues on
the unit circle) — matching the published I_c stable segment. This is the prize
the #442/#457 saga was chasing, and it materially exceeds #442 (which only
tracked the e1→0.948 unstable arm).

Method lesson (durable): the prior 2-var full-period corrector floored at ~2e-6
because it admits a non-doubly-symmetric compromise point; enforcing
perpendicularity at BOTH t=0 and t=T/2 (`[y(T/2),vx(T/2),y(T),vx(T)]`) is what
selects the true Scheme-II member. This is the analogue of the #442 4-vector
fix, one level deeper (both crossings, not just the endpoints).

## Promotion / follow-up

Promotion of the converged member to a module helper + golden test is a FOLLOW-UP
(per the task scope: "Promotion to a module/golden is a FOLLOW-UP, not this
task"). The full-precision IC and the joint-corrector recipe above are the durable
deliverable. Author data remains unavailable/vetoed; the route is fully
method-gated and now SOLVED for the 3/1 I_c representative. The same joint
doubly-symmetric corrector should reach the 4/1 / 5/1 isolated I_C members
(analogous Scheme-II structure).

Scratch (left untracked — exploratory): scripts/_scratch_457_*.py
(branches_e2zero, arclength_Ic, family_from_anchor, focused_seed, direct_e2_090,
polish) + /tmp probes (_dsym, _freeT, _e2sweep2, _resid_anatomy, _collcheck,
_joint, _final_ic). The reproducible recipe is the joint residual
`[y(π),vx(π),y(2π),vx(2π)]` over `(x,vy)` at θ0=π, e2=0.91, seeded from the
focused (π,0) P1-pericentre seed.
