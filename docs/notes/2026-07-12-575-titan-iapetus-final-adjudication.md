# #575 Titan-Iapetus symmetric-closure thread -- FINAL ADJUDICATION

Date: 2026-07-12. Adjudicator: Opus (culminating trust-bearing call on the
#571->#572->#573->#574->#575 Titan-Iapetus quasi-cycler thread). Inputs read in full:
`data/OUTSTANDING.md` #571-#575, `docs/notes/2026-07-12-{572,573,574-stageA,574-stageB}`,
and the raw probe data `data/probe_575_stage1_inclination_closure.jsonl` +
`data/probe_575_stage1b_omega_only_closure.jsonl`.

## Decision: STOP. Stamp a clean, well-instrumented, method-conditional negative.

This is a **novel** internally-enumerated search (NOT a published-paper reproduction we
are failing to reproduce), so a clean negative is a legitimate terminus per project
policy. The negative is honestly conditional on the search *method* (continuation from
circular seeds); the one method not tried -- a from-scratch node-locked inclined
symmetric construction -- is recorded as the explicit, method-versioned re-open condition
in `data/empty_regions.jsonl`, NOT dispatched as a live task. Rationale below.

Framing throughout: quasi_cycler-CLASS evidence about our own idealized search space
(same standing as #312). No novelty claim. No `catalogue.yaml` edit.

## 1. The near-miss numbers -- robustly infeasible, not a tantalizing near-miss

The coplanar #563-method construction found **9 genuine symmetric closures** that repeat
to machine precision by construction (C2 positive control: 4e-14 to 3e-12 km/s over 3
cycles; negative control #571-branch-1 correctly fails to complete even 1 cycle). That
pipeline demonstrably *can* find genuine repeating cyclers when they exist -- so the
downstream inclination negative is trustworthy, not a broken filter.

Pushing the 9 coplanar seeds to Iapetus's real ~15.5deg inclination was tried **two
ways**, both of which fail decisively and in *opposite* directions -- the signature of a
robust structural non-closure, not a solution one refinement iteration away:

**Stage 1 -- Omega + tof_scale free (continuation, #572 machinery verbatim):** 6/9 seeds
find a single-cycle closure to machine precision (refined residual 3e-10 to 2e-9 km/s).
But the moment the closed cycle is actually *propagated forward* (C2 3-cycle repeat
check), it diverges:

| seed | tof drift off n*T_syn/2 | 3-cycle max residual | 3-cycle max drift |
|---|---|---|---|
| n2 rel0  | 0.23%  | 0.108 km/s | 2,443,740 km |
| n2 rel180| 0.23%  | 0.150 km/s | 2,443,740 km |
| n6       | 0.70%  | 0.280 km/s | 2,443,740 km |
| n8       | 4.64%  | 1.651 km/s |   702,966 km |
| n10 rel0 | 10.5%  | (Lambert infeasible past cycle 0) | -- |
| n10 rel180|11.6%  | (Lambert infeasible past cycle 0) | -- |

A drift of **2.44 million km is ~2x Titan's own orbital radius** (1.22M km) -- the
trajectory wanders completely off, not a small perturbation. The single-cycle
machine-precision closure is misleading: it closes one transfer by letting tof drift, but
the state it returns to does not repeat. The two highest-commensurability seeds (n=10)
lose Lambert feasibility outright past cycle 0 -- the exact #571-branch-1 failure
signature, now reproduced by refinement-induced drift.

**Stage 1b -- Omega only, tof HARD-fixed at the exact commensurate value (zero drift by
construction):** the more-constrained variant. 5/9 seeds find a single-cycle closure, and
in every case Omega lands on the symmetric axis (0.0 / 180.0 / 359.9 / 179.9 deg) -- the
search is trying to align Iapetus's line of nodes with the symmetric-crossing axis. But
even there the single-cycle residual is only ~1e-3 km/s (grid resolution, NOT machine
precision), and the 3-cycle repeat residual grows to **0.006-0.021 km/s (6-21 m/s)**.
None repeat. So periodicity loss is NOT merely the Stage-1 tof-drift mechanism: Iapetus's
real inclination breaks the multi-cycle repeat property in its own right.

**The honest read:** the tightest, most-constrained search leaves a ~6-21 m/s per-cycle
velocity discontinuity that it *cannot* drive to zero; the freer search diverges to
millions of km. Both bracket the true inclined solution and neither converges to it. The
gap to the Uranian precedent is ~12 orders of magnitude: the Umbriel-Oberon coplanar
closures repeat at 4e-15 km/s; the best Titan-Iapetus inclined case is 6e-3 km/s. That is
not "harder to close" -- it is a different object.

## 2. Is continuation from circular seeds fundamentally limited here? Yes.

The coplanar symmetric family is periodic because of a **mirror symmetry** about the
common orbital plane (rel_offset in {0,180deg} + commensurate tof = the classical
perpendicular-crossing condition). Tilting Iapetus by ~15deg breaks that mirror symmetry
*except* in the special configuration where the encounters occur at the mutual line of
nodes (where the two orbital planes intersect). Stage 1b's Omega -> {0,180deg} is exactly
the search reaching for that node alignment -- but the commensurate-timing condition and
the node-crossing condition are generically **incompatible** at Titan-Iapetus's period
ratio (~4.975, non-resonant), which is why the residual never nulls.

A genuine inclined symmetric family would be **node-locked**: a codimension-higher,
structurally *disconnected* family that smooth deformation of a coplanar seed cannot
reach. So the Stage-1/1b negative is formally conditional on **connectivity** -- it
proves the coplanar family does not survive continuation to inclination; it does NOT
prove no disconnected inclined family exists. That residual question is the honest re-open
condition, recorded in the registry.

## 3. Uranian precedent -- an honest, different-regime comparison

This is the load-bearing physical point, and it is NOT "the same correction, harder":

* **Umbriel-Oberon (#312):** both moons orbit within ~0.1-0.15deg of Uranus's
  equatorial/Laplace plane -> **mutual inclination ~0.1-0.2deg**. The circular-COPLANAR
  idealized model is *near-exact* for this pair; the inclination-extension is a negligible
  correction, which is precisely why the coplanar symmetric closures survive it and
  reproduce to machine precision.
* **Titan-Iapetus:** Titan sits ~0.3deg off Saturn's Laplace plane, but Iapetus -- far out
  at a=3.56M km -- orbits near its *local* Laplace plane, tilted ~8deg from Saturn's
  equator and ~15.5deg in the reference frame the probe uses, giving a **mutual
  inclination of ~15deg** between the two orbital planes. The circular-coplanar idealized
  model was never a good starting point for this pair; inclination is a **leading-order
  effect**, not a correction.

So the honest conclusion is the sharper one the dispatch flagged: *Uranian moons are
nearly coplanar with each other, so the idealized model was almost realistic for them;
Titan-Iapetus's moons are NOT nearly coplanar, so the idealized coplanar model was never
going to be a good starting point.* This is a genuinely different dynamical regime, not a
tougher instance of the same one.

## 4. Cost/benefit of a from-scratch inclined search -- not worth dispatching now

The whole #572->#575 probe sequence was explicitly scoped as *"the cheap gate deciding
whether any #552 build is worth scoping"* (#552 = the general 3D/inclined-releg moontour
genome extension, explicitly KILLED on 2026-07-10). A from-scratch node-locked inclined
symmetric construction IS essentially that killed build (or very close to it): it requires
deriving and coding the inclined mirror-symmetry / node-locked closure condition -- genuine
new capability, not reuse of the existing coplanar #563 construction.

The gate has now returned a robust, positive-control-validated negative **twice**:
* #574 Stage B: 0/15 under real eccentricity, root-caused to the #571 periodicity gap.
* #575: 0/6 (and 0/5 in the zero-drift variant) under inclination *alone*, using the
  *correct* symmetric method -- an earlier-stage, more-fundamental breakdown than #574's.

Expected yield of the expensive search is low: even if a marginal node-locked inclined
family exists, it must still survive the eccentricity kill-gate (#574 Stage A) and the
Saturn V2-V4-strict gauntlet (#574 Stage B) that already killed 0/15 of the #571
population, with Titan-Iapetus eccentricity perturbation 7-25x stronger than the
negligible Uranian effect (per #575's own C4 Fable note). And it is not novelty-eligible.
Building the very thing the gate was designed to gate, to chase a low-probability marginal
result with no catalogue standing, would invert the disciplined structure of this thread.

Per the empty-region registry's method-versioned re-open rule, the right disposition is to
**record** the from-scratch node-locked inclined symmetric search as the named, licensed
re-open method (a strictly-more-capable method that would subsume this continuation-based
sweep) -- without dispatching it. If a future, independently-motivated inclined-symmetric
capability is built for another system, it can cheaply re-check this pair then.

## Verdict

**EMPTY** for inclined, repeating Titan-Iapetus symmetric quasi_cyclers reachable by
continuation from the coplanar symmetric family. Nine genuine coplanar symmetric closures
exist in the idealized circular-coplanar model, but none survives extension to Iapetus's
real ~15.5deg inclination as a repeating cycle. Clean negative, thread closed.
Re-open condition: a from-scratch node-locked inclined symmetric construction (disconnected
from the coplanar family; not built here, not dispatched).
