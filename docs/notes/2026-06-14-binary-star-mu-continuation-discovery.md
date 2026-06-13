# Binary-star mass-parameter (mu) continuation of stable EM cyclers (discovery)

**Date:** 2026-06-14
**Task:** #252 (autonomous run).
**Source motivation:** Roberts-Tsoukkas & Ross 2026 (journal: "Stable Prograde
Earth-Moon Multi-Orbiter Cyclers via Three-Body Dynamics"; journal version of
Ross & Roberts-Tsoukkas 2025, AAS 25-621). Its Figure 3 depicts -- in FIGURES
ONLY, **no printed numbers** -- stable prograde cyclers persisting up into the
binary-star mass-parameter range: a **(1,3) exterior cycler at mu = 0.1**, a
**(3,1) cycler at mu = 0.3**, and a **(1,1) equal-mass cycler at mu = 0.5**, all
drawn stable. Figure 2 separately depicts the (1,1) and (3,3) Earth-Moon
cyclers "near the saddle-center bifurcation of their respective families".

**Writeback: NONE.** The recovered `(mu, C, T, IC, nu)` are OUR OWN computed
values (the paper prints none for any binary-star member) -> they are
**DISCOVERIES** requiring the full V0-V5 gauntlet, not sourced rows. This note
is review-gated evidence; the catalogue is untouched.

## What was built

- `src/cyclerfinder/search/mu_continuation.py` -- pseudo-arclength continuation
  of a symmetric (k1,k2) cycler in the mass parameter mu, plus a fixed-mu
  C-family stable-window scanner.
- `tests/search/test_mu_continuation.py` -- a known small-mu step recovers a
  genuine periodic orbit (residual ~ machine precision, independent-Radau
  re-closure, Jacobi self-consistency, symmetric-IC structure, continuity);
  the fixed-mu C-scan finds the stable subfamily around the held nu=0 member.
- `scripts/cr3bp_mu_continuation.py` -- the discovery driver (emits the runlog
  `out/outcome_log/mu*_discovery.txt`; gitignored `out/`).

## Method

A perpendicular-x-axis-crossing symmetric cycler has IC `(x0,0,0,0,ydot0,0)`
with `ydot0` fixed from the Jacobi constant C (Ross Eq. 9). For a fixed
`(k1,k2)` the single periodicity condition is `xdot(t_half)=0` at the family's
half-period crossing (Ross Eq. 11). The solution set
`r(x0,C,mu) = xdot(t_half) = 0` is a 2-surface in `(x0,C,mu)`; a family at fixed
mu is a curve on it (parameterised by C). We follow a path from the Earth-Moon
mu to the target mu by **pseudo-arclength continuation**: predict along the
local tangent (null vector of dr/dz), correct back onto `r=0` with the
arclength constraint `tan.(z - z_pred)=0` (min-norm Newton, 2 eqns / 3 unknowns),
adaptive step. An exact final landing on `mu = mu_target` is done by a fixed-mu
projection. Each kept member: Barden half-period-monodromy
`nu = 1/2(lambda+1/lambda)` (Ross Eqs. 13-15; |nu|<1 = linearly stable) AND an
independent-Radau full-period re-closure (different integrator than the DOP853
corrector).

Natural-parameter continuation (fix mu, solve x0 OR C alone) was tried first and
**rejected**: with C frozen the corrector traces a C=const slice (drifts off the
member, wildly unstable); with x0 frozen likewise; full-state shooting wandered
to unrelated far orbits. Pseudo-arclength is the only scheme that follows the
actual branch (and can turn folds). [Decisions under ambiguity, below.]

## Results -- family (3,1), mu: 0.01215 -> 0.3

Arclength continuation of the held (3,1) nu=0 midpoint
(`ross-rt-em-cycler-31-2025`), all members genuine periodic orbits
(crossing residual <= 3e-11, independent-Radau Jacobi drift <= 6e-13):

| mu | C | T (TU) | nu | verdict |
|------|------------|-----------|-----------|----------|
| 0.01215 | 3.161784147 | 14.788268 | +0.01545 | STABLE (seed) |
| 0.02520 | 3.158396688 | 15.349447 | +39.398 | unstable |
| 0.09788 | 3.133129905 | 15.831205 | -0.70421 | STABLE |
| 0.18270 | 3.113788765 | 16.246254 | -1.71023 | unstable |
| 0.31381 | 3.094940232 | 16.252893 | -300.10 | unstable |

The branch REACHES mu = 0.3 (genuine periodic orbit at every mu; landed member
crossing residual 3.2e-12, independent-Radau Jacobi drift 4.9e-13), but the
**held nu=0 (3,1) branch oscillates in/out of linear stability** and is strongly
unstable by mu ~ 0.3 (nu = -208 at the landing). Stability is NOT monotone:
stable at mu=0.0121 and mu=0.098, unstable between and beyond.

**Topology check (the decisive figure-match test).** Counting U1-/U2+ crossings
(Ross Def. 1) over one period:

| member | mu | T | (k1,k2)~ | x-range | P2 at | cycler? |
|--------|------|--------|----------|---------|-------|---------|
| (3,1) seed | 0.0122 | 14.788 | (3,1) | [-0.63, +0.99] | +0.99 | YES (reaches P2) |
| (3,1) landed | 0.3 | 16.283 | (3,0) | [-0.74, -0.18] | +0.70 | **NO** (never reaches P2) |

The landed mu=0.3 orbit has **lost the lunar/secondary encounter** (k2 -> 0):
its trajectory is confined to x in [-0.74, -0.18], entirely on the primary side
of P2 (at +0.70). It is a genuine periodic orbit but **no longer a (3,1)
cycler**. A fixed-mu C-family scan at mu=0.3 across C in [3.02, 3.17]
(51 converged members) found **0 linearly stable members** (all nu > +140,
monotone) -- and these all share the same one-sided, non-cycler topology.

**(3,1) verdict: NEGATIVE.** The nu=0-seeded (3,1) branch does not deliver the
paper's stable mu=0.3 cycler; it leaves the cycler regime and the stable window
before mu=0.3. The paper's depicted stable (3,1)@0.3 must live on a different
branch (different x0/C region) not reached from this seed.

## Results -- family (1,1), mu: 0.01215 -> 0.5

Arclength continuation of the held (1,1) nu=0 midpoint
(`ross-rt-em-cycler-11-2025`):

| mu | C | T (TU) | nu | verdict |
|------|------------|-----------|-----------|----------|
| 0.01215 | 3.151175880 | 10.292069 | -0.00334 | STABLE (seed) |
| 0.02039 | 3.152574230 | 15.869088 | +41.272 | unstable |
| 0.09755 | 3.142919760 | 16.054149 | +1.12520 | unstable |
| 0.15416 | 3.106001851 | 16.378558 | +2.68094 | unstable |
| 0.18613 | 3.106931284 | 16.610056 | +21.340 | unstable |
| 0.33551 | 3.133841251 | **2.05729** | -0.51264 | STABLE |
| 0.43946 | 3.131003020 | 2.098012 | -0.55631 | STABLE |
| 0.50000 | 3.129232191 | 2.130037 | -0.58705 | STABLE |

The branch REACHES mu = 0.5 (equal mass) and LANDS on a STABLE periodic orbit
(nu = -0.587; crossing residual 4.6e-12, Radau drift 6.4e-13). A C-family scan
at mu=0.5 found a WIDE stable window: **51/51 scanned members linearly stable**
(nu in [-0.65, -0.52]).

**BUT the topology test fails.** Note the **period collapse** between mu=0.186
(T=16.6) and mu=0.336 (T=2.06): the arclength path passed through a fold and
**branch-switched** onto a short-period orbit. Crossing counts of the landed
mu=0.5 member:

| member | mu | T | (k1,k2)~ | x-range | P2 at | cycler? |
|--------|------|--------|----------|---------|-------|---------|
| (1,1) seed | 0.0122 | 10.292 | (1,1) | [-0.77, +1.05] | +0.99 | YES (reaches P2) |
| (1,1) landed | 0.5 | 2.055 | (0,0) | [-0.696, -0.307] | +0.50 | **NO** (one-primary librational) |

The landed mu=0.5 orbit is a small **librational orbit around the primary P1**
(x in [-0.696, -0.307], never reaching P2 at +0.50); (k1,k2) = (0,0). It is a
genuine, stable periodic orbit at equal mass, but it is **NOT the depicted (1,1)
cycler** -- the continuation fell off the cycler branch at the fold near
mu ~ 0.2 (where the genuine (1,1) cycler branch, period ~16, had already gone
unstable).

## Figure-match verdict

**NEGATIVE for both attempted families (honest result).** The nu=0-seeded
Earth-Moon (1,1) and (3,1) branches were continued as genuine periodic orbits
(machine-precision residuals, independent-Radau closure at every recorded mu)
all the way to the binary-star targets mu=0.3 and mu=0.5 -- but **neither lands
on the stable cycler the paper draws**:

- the (3,1) branch loses the secondary encounter (k2 -> 0) and stays unstable;
- the (1,1) branch branch-switches at a fold onto a one-primary librational
  orbit (loses both encounters), which IS stable but is not a cycler.

The genuine (1,1)/(3,1) *cycler* branches (period ~16 at mid-mu) go LINEARLY
UNSTABLE well before the targets (nu growing through +1 around mu ~ 0.02-0.1).
The paper's claim that a stable subfamily exists at each target mu is not
contradicted -- it asserts a stable *window* exists *somewhere* on the C-family
at each mu; our scans show it is NOT the analytic continuation of the EM nu=0
member, and not in the C-windows we scanned around the continued branch. Finding
it would require seeding the fixed-mu C-scan from the paper's saddle-center
(nu near +1, near C^max) topology at the target mu, or a 2-parameter
(mu, C) sweep -- a larger search than this bounded run.

## Proposed novel-cycler candidate(s) + gauntlet entry point

**No figure-matched stable binary-star CYCLER candidate is proposed** (the
honest negative above). Two by-products are recorded for completeness, NOT as
catalogue candidates:

1. A stable equal-mass (mu=0.5) one-primary librational periodic orbit
   (state0 = [-0.3073442, 0,0,0, -1.8219388, 0], C=3.20423, T=2.0549,
   nu=-0.521) -- genuine + stable but topology (0,0), not a cycler; out of
   catalogue scope (cyclers only).
2. The (3,1) and (1,1) periodic orbits that DO exist at every continued mu (the
   tables above) -- genuine periodic orbits but either unstable or
   non-cycler-topology at the targets.

Were a genuine figure-matched stable binary-star cycler found in future work, its
gauntlet entry point would be V0: a self-computed periodic orbit with an
independent-integrator closure + a Barden stability verdict in a NON-physical mu
regime. It could never be sourced-V1 (the paper prints no numbers), so it would
be a true discovery on the V0-V5 ladder, not a sourced row.

## Decisions under ambiguity

1. **Which families to continue.** The journal targets are (1,3)@0.1,
   (3,1)@0.3, (1,1)@0.5. We hold EM seeds for (1,1),(2,1),(3,1),(3,2),(3,3) but
   NOT (1,3). So the two same-class continuations are (3,1)->0.3 and
   (1,1)->0.5; (1,3)@0.1 has no same-class EM seed and was not attempted
   (would require first constructing a (1,3) EM member -- out of scope here).
2. **Continuation scheme.** Pseudo-arclength in (x0,C,mu); natural-parameter
   schemes rejected (see Method). The symmetric corrector is used throughout to
   preserve the perpendicular-crossing family identity; full-state shooting was
   rejected (wanders off-family).
3. **Member selection across mu.** The arclength path follows the branch through
   the held nu=0 EM member. That branch leaves the stable window (and, at the
   targets, the cycler topology) before reaching mu=0.3/0.5, so a fixed-mu
   C-family scan at the target mu was used to look for the stable subfamily (the
   paper's own method). The scans found no figure-matched stable cycler near the
   continued branch -> reported as a clean negative, NOT forced.
4. **Honest negative is acceptable and is the outcome here.** Both attempted
   families gave genuine periodic orbits at every mu but no stable binary-star
   cycler matching the figures. The exact failure points are localized (fold /
   branch-switch for (1,1) near mu~0.2; secondary-encounter + stability loss for
   (3,1) before mu=0.3). Per task guidance, a clean negative with the failure
   point is a valid result; no discovery was manufactured.
5. **Topology is the binding figure-match test**, not residual or stability
   alone. "It closed and it's stable" was the danger signal that the mu=0.5 land
   was NOT a cycler -- the U1-/U2+ crossing count (Ross Def. 1) caught it.

## Provenance / discipline

Pure PCR3BP (no ephemeris). EXPECTED test quantities are intrinsic correctness
properties (residual~0, Jacobi conserved), never a self-computed value asserted
against itself. The source paper is filed in the private reference repo; no path
to it appears here.
