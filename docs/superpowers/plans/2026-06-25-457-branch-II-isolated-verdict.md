# #457 — 3/1 isolated stable family (I_c) via Scheme-II e2-continuation — VERDICT

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

5. **The STABLE I_c representative basin is located** (focused seed, e2=0.88,
   θ0=π, P1 pericentre, positive crossing): `x≈+0.16271, vy≈+3.14973`,
   a1=0.480313, **e1=0.657318, |eig|max=1.00249 STABLE**
   (eigs 0.980±0.200j, 1.0025, 0.9976) — matching the published representative
   `a1/a2=0.480674, e1=0.659951`. HOWEVER its fixed-(e2,θ0) 2-var
   perpendicular-crossing residual FLOORS at **5.0e-6** and will not descend to
   the 1e-8 gate (verified insensitive to integrator tol/method/max_step, and
   NOT a near-collision — min P1–planet distance 0.68). The non-closure is a
   small residual perpendicularity defect at t=T/2 (y≈5e-5, vx≈5e-4): the exact
   member needs a freed parameter (e2/θ0/period), not just (x,vy).

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

The floor BOTTOMS OUT at ~1.8e-6 (e2≈0.91) and never reaches the 1e-8 gate.
Escalations that did NOT break it (each ruled out a cause):
- finer integrator tol (1e-12→3e-14), max_step, Radau: floor invariant ⇒ NOT an
  integration-accuracy limit;
- free-period (x,vy,T): state residual only 5e-6→2.4e-6, T=6.28301 vs 2π=6.28319
  ⇒ period freedom does not close it;
- θ0 sweep {0, π/2, π, 3π/2}: θ0=π is the deep minimum (others 1e-3) ⇒ phase is
  correct, not the cause;
- doubly-symmetric HALF-period corrector (zero y,vx at t=π): drives half-residual
  to 8.6e-7 but RAISES full-period residual to 1.8e-3 ⇒ the member is NOT exactly
  perpendicular at t=π (t=π state has y≈−6e-6, vx≈−3e-4): it is not a clean
  θ0=π doubly-symmetric perpendicular-crossing PO in this ansatz;
- close-approach check: min dist to star 0.165 (= the P1 pericentre radius), to
  planet 0.706 ⇒ NOT a near-collision / regularization issue.
The irreducible residual is dominated by the y(T) component (1.84e-6); dx, vx(T),
dvy(T) are 1e-8–1.6e-7.

Stability of the best (e2=0.91, full precision) member: monodromy eigenvalues are
a complex pair 0.97±0.23j (|·|=1.000) plus a near-(1,1) pair; |eig|max=1.000–1.007
(FD steps 1e-7, 1e-6), det≈1.0, trace≈3.945 — **consistent with STABLE** (on the
unit circle), matching the published "I_c stable for 0.75<e1<0.98". But because
the member closes only to ~2e-6, the stability is CONSISTENT-WITH-stable, not
certified.

## Verdict — CHARACTERIZED PARTIAL (family reached; exact stable-member closure blocked)

**Reached:** the 3/1 Scheme-II I_c family IS recovered in the validated paper
frame — clean 1e-12 closures on its unstable arm and through a reproduced FOLD
(matching the published folded Fig 11e structure), and the STABLE high-e segment
representative is located at the published config (π,0) with a1≈0.4804 (pub
0.480674) and e1≈0.655 (pub 0.659951), eigenvalues on the unit circle
(consistent with the published stable I_c). This is materially further than #442
(which only tracked the e1→0.948 unstable arm and never located the stable
representative).

**Blocked:** the stable representative does NOT close to the 1e-8 independent gate
— it floors at ~1.8e-6 (clean e2-minimum at 0.91), and the escalation ladder
(integrator tol/method, free-period, θ0-sweep, half-period doubly-symmetric,
collision check) rules out integration error, period, phase, and collision. The
half-period test shows the member is not an exact θ0=π perpendicular-crossing
doubly-symmetric PO. The most likely cause is an ANSATZ mismatch: config (π,0)
for the q=2 / 3/1 resonance uses the resonant-angle pair (θ3,θ1) (arXiv §"3/1",
Eq. 5), whose symmetry line is NOT the plain x-axis perpendicular crossing the
corrector enforces — so the exact symmetric member sits slightly off the
`[x,0,0,vy,θ0=π]` section. This is the same steep-basin / symmetry-ansatz wall
that walled #442's stable member, now precisely localized to the symmetry-section
definition.

**Honest status:** a genuine, well-characterized NEGATIVE on the 1e-8-closure
goal, with a strong POSITIVE on family reachability — the I_c family and its
stable segment are demonstrably reached at the published a1/e1, the published
fold is reproduced, and the only remaining gap is the exact symmetry-section for
the (θ3,θ1)=(π,0) config. NOT promoted to a module/golden (the gate is not met).

## Concrete next method (the unexhausted rung)

Build the corrector on the paper's ACTUAL (θ3,θ1) symmetry condition for the q=2
3/1 resonance (Eq. 5; the doubly-symmetric section for (π,0) is defined by the
resonant angles, not the bare x-axis perpendicular crossing). Re-seed the located
e2≈0.91 member (x≈0.16363, vy≈3.14213) into that section's residual; the ~2e-6
y(T) defect should collapse if the section is the correct one. This needs the
explicit (θ3,θ1) → state map, a bounded follow-up build — NOT a re-run of the
present ansatz. Author data remains unavailable/vetoed; the route stays
method-gated, and this is the next method.

Scratch: scripts/_scratch_457_*.py (branches_e2zero, arclength_Ic,
family_from_anchor, focused_seed, direct_e2_090, polish) + /tmp probes
(_dsym, _freeT, _e2sweep2, _resid_anatomy, _collcheck, _finalstab). Left
untracked (gate not met).
