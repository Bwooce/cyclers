# Multi-arc basin-selection — scoped, ranked design (task #162)

**Date:** 2026-06-08
**Type:** design pass (no production code; deliverable = this note + the executable
plan `docs/superpowers/plans/2026-06-08-multiarc-basin-selection.md`).
**Acceptance discipline (binding):** the target for any approach is a row's
**SOURCED** V∞ anchor (6.44Gg3 → E 6.44 / M 3.74; S1L1 `russell-ch4-4.991gG2` →
E 4.99 / M 5.10, its own framing's anchors — never mixed with the 5.65/3.05 or
4.7/5.0 framings, multi-arc-classification §7/§12). Emerged values are EVIDENCE,
never imposed. A clean NEGATIVE (the basin provably does not exist for the
topology) is a valid, publishable outcome; an "inconclusive" outcome is not, and
no approach below is recommended that can only ever return inconclusive.

**No catalogue writeback in scope.**

---

## 0. The elimination chain (so we do not re-propose dead ends)

| pass | genome / method | row | result |
|---|---|---|---|
| #137 | single-ellipse free-return `(a,e,t0)` | 6.44Gg3 | CLOSE-OFF-ANCHOR (emerged 3.01/3.06 vs 6.44/3.74) — "genuinely different family" |
| #145 | MBH over single-ellipse | 6.44Gg3 | Gate-3 negative; 0 hops accepted — single-ellipse basin does not exist here |
| #158 | continuation (circular→DE440) over single-ellipse | S1L1 | strict-xfail; at sourced window the genome collapses to `e→0.95`, residual 23–46 km/s |
| #150 | one-DSM-per-leg, minimal E-M-E (control) | 6.44Gg3 | floor 9.40 km/s, 0/61 accepted |
| #153 | one-DSM-per-leg, full E-M-E-M-E, single-rev | 6.44Gg3 | floor 29.9 km/s — WORSE; mechanism: single-rev Lambert degenerate on the >1-period loop arcs |
| #157 | one-DSM-per-leg, full sequence, **multi-rev** Lambert | 6.44Gg3 / S1L1 | floor 26.9 / 38.9 km/s, 0/61 accepted — representable now, still off-basin |

**The consistent, multiply-confirmed bound:** the topology is now *representable*
(the multi-rev one-DSM-per-leg genome can express two generic loop arcs), but the
search *floors off-basin* — `Σ ΔV_DSM` stays 27–39 km/s with **zero** accepted
hops and emerged encounter V∞ 3–37 km/s, nowhere near the sourced anchors.
"Representable ≠ reachable." The frontier is **seed/basin selection**, not
mechanics.

### 0.1 The precise mechanism (read from the live evaluator)

`evaluate_dsm_chain` (`src/cyclerfinder/search/dsm_leg.py:357`) has TWO structural
properties that jointly explain the off-basin floor:

1. **It FORCES heliocentric continuity at every intermediate flyby**
   (line 448: `v_depart = leg.v_arrive`). So V∞ *magnitude* is continuous across
   the flyby **by construction** — the corrector never has to satisfy a
   `|v∞⁻| = |v∞⁺|` constraint, because the next leg literally departs on the
   arrival velocity. What it pays for that continuity is the **DSM impulse** that
   bends the *next* leg's Lambert arc to reach the next body on time.
2. **The single scalar residual is `Σ ΔV_DSM` (+ terminal arrival)** (line 474).
   There is **no turn-angle / bend-feasibility term at all**. The corrector is
   free to land a "continuous" chain whose intermediate flyby demands an
   *infeasible* turn (the spacecraft arrives and departs at the same heliocentric
   speed but in directions no real periapsis can connect), and pays nothing for
   it — the bend cost is deferred to the catalogue scorer (`core/flyby.py`), which
   never runs inside the search.

This is the exact **"#110 / S1L1 summed-mismatch lets the long leg dominate"**
failure mode, sharpened: the scalar `Σ ΔV_DSM` is dominated by whichever leg's
Lambert is most stressed (the long loop arc), so `least_squares` spends all its
gradient buying down that one leg's impulse — and the geometry it lands in is a
*kinematically* continuous but *dynamically* infeasible chain at high V∞, not the
sourced low-V∞ bend-feasible cycler. MBH then hops, but every hop lands in the
same funnel (0/61 accepted) because the objective has no gradient pointing at the
sourced basin: the sourced anchor is a *bend-feasible, low-V∞* point that the
objective does not reward.

**Corollary that drives the ranking:** the most likely single fix is not "more
search" (MBH already cures selection *within* a transcription, #145) and not "more
revs" (already added, #157) — it is **changing the objective so the sourced basin
is the one with the lowest value**, i.e. adding the flyby V∞-continuity-AND-bend
residual the evaluator currently omits. That is approach #3, and it is the only
approach that directly attacks the diagnosed mechanism.

---

## 1. Approach inventory — scope, experiment, failure mode

### Approach 1 — Branch as a discrete genome coordinate

**Idea.** The per-leg `(n_revs, branch)` is currently chosen *inside* `dsm_leg` by
ΔV-min (greedy/local, #157). Promote it to an outer **discrete** search dimension:
MBH (or a small enumeration) over `(continuous genome) × (branch tuple per leg)`,
each candidate branch tuple corrected independently, keep the best.

**Files/functions.** `dsm_leg` already accepts `rev_branch=(n,"low"|"high")` per
leg (#157, `dsm_leg.py:129`). New work: thread an explicit `rev_branch_per_leg`
through `evaluate_dsm_chain` / `dsm_chain_correct` / `make_dsm_chain_step` (today
they pass only `max_revs` and let each leg re-minimise); add an outer loop in a
probe script that enumerates the small Cartesian product of feasible branches
(e.g. `{single, (1,low), (1,high), (2,low), (2,high)}` on the two loop legs) and
runs the existing MBH per tuple. **No `core/` change.**

**Effort.** LOW (threading + an enumeration loop; the primitive already supports
forced branches).

**Cheapest decisive experiment.** Re-run the #157 6.44Gg3 full-sequence probe, but
instead of one greedy run, enumerate the ~9–25 branch tuples on legs 2 and 4 and
MBH each. Gate below.

**Honest failure mode.** This only changes *which Lambert family* each leg sits in;
it does **not** add a bend-feasibility term. If the off-basin floor is caused by
the missing continuity/bend residual (the diagnosis in §0.1), enumerating branches
will lower the floor modestly (as #157's `(1,low)` already did: 29.9→26.9) but will
**not** cross into the sourced basin. Expected outcome on the diagnosis: NEGATIVE,
floor improves but 0 accepted at anchor. It is a necessary *ingredient* of the
hybrid (#5), not a standalone cure. Risk it returns "inconclusive": LOW — a floor
number + emerged V∞ is always decisive evidence.

### Approach 2 — Symmetric / period-locked arc seeding

**Idea.** Russell's generic-return arcs are (often) *symmetric* (t_in = t_out about
the arc apoapsis) and their ToFs are descriptor-sourced. Seed the two loop arcs at
their true symmetric half-arc geometry (each arc's DSM at η=0.5 = the apoapsis,
ToF split symmetrically), and seed the departure V∞ at the sourced Earth anchor —
mirroring the #137 free-return *radial-crossing* seed idea, but at the multi-arc
(per-arc) level. Today the probes seed `(transit, arc−transit)` asymmetric splits
(#153 table) and η=0.5 free.

**Files/functions.** A seeding helper (probe-side, or a new
`symmetric_arc_seed(...)` in `dsm_leg.py`) that maps each `free_return_arcs[]`
descriptor (ToF in years, ψ angle) to a symmetric per-leg `(tof_days, eta, alpha,
beta)` seed. Reuse `seed_ae_from_aphelion_transit` (`free_return.py:247`) as the
sibling pattern. Feeds the existing `dsm_chain_correct` / MBH unchanged.

**Effort.** LOW–MED (the descriptor→seed map; ψ→α/β geometry needs care).

**Cheapest decisive experiment.** Seed S1L1 4.991gG2 symmetrically from its
descriptor `g(1.4612…)·G(2.8096…)` and run the existing corrector ONCE (no MBH) —
does the symmetric seed land nearer the basin than the #157 asymmetric seed
(floor 38.9)?

**Honest failure mode.** Seeding changes only the *starting point*, not the
objective. If the basin is real but narrow, a better seed helps the local solve
find it — but #135 already showed (single-ellipse) that *the sourced geometry is
not a residual-zero point of the genome*; if the same holds for the multi-arc
objective (likely, given §0.1: the sourced point is bend-feasible-but-not-ΔV-min),
the symmetric seed will *start* near the anchor and the corrector will *walk away*
to the same high-V∞ funnel. Expected on the diagnosis: NEGATIVE — seed lands near
anchor, residual at the seed is large (the decisive datum), corrector diverges.
That "residual at the sourced seed" number is itself the publishable evidence (it
is the multi-arc analogue of the #135 seed-at-truth probe). Risk of inconclusive:
LOW.

### Approach 3 — Explicit flyby V∞-continuity + bend-feasibility residual ★

**Idea.** Add the intermediate-flyby physics the evaluator currently omits as an
**explicit residual term** the corrector drives to zero, so the sourced
bend-feasible low-V∞ chain becomes the *minimum* of the objective rather than an
unrewarded point. Concretely, STOP forcing `v_depart = leg.v_arrive`; instead let
the next leg's departure V∞ direction be a (bounded) genome coordinate, and add per
flyby a residual:

```
r_flyby = w_mag · ( |v∞⁺| − |v∞⁻| )            # magnitude continuity (Eq.5.1)
        + w_bend · flyby_dv(v∞⁻, v∞⁺, μ, rp)   # turn-feasibility (Eq.5.5 surrogate)
```

`core/flyby.py` ALREADY has both pieces ready to wire in: `flyby_dv`
(`flyby.py:356`, the Russell Eq.5.5 powered-SOI surrogate that returns 0.0 when
ballistic-feasible and a positive Δv otherwise) and `dv_powered_flyby_periapsis`
(`flyby.py:117`, Oberth). The residual vector becomes
`[*ΔV_DSM_per_leg, *flyby_dv_per_flyby, arrival]` — a true multi-term root-find,
not a single dominated scalar. This converts the objective so that the sourced
anchor (low V∞, bend-feasible → `flyby_dv≈0`, small DSM) is the global minimum the
optimiser is actually pulled toward.

**Files/functions.** `evaluate_dsm_chain` (`dsm_leg.py:357`): add an optional
`charge_flyby_continuity=True` path that (a) frees the departure-V∞ direction per
intermediate leg (new genome coords, mirroring the leg-0 `alpha/beta`), (b) emits
the per-flyby `flyby_dv` term, (c) returns a *vector* residual. `dsm_chain_correct`
(`dsm_leg.py:614`): change `_res` to return the full residual vector (already uses
`least_squares`, which natively does vector residuals — today it returns a
length-1 array). `dsm_chain_decision_vector` / `_unpack` / `DsmBounds` /
`sequence_keyed_bounds`: extend layout by `2·(n_legs−1)` direction coords. Import
`flyby_dv` from `core/flyby.py` (read-only; no `core/` edit). MBH adapter
unchanged in shape.

**Effort.** MED (it is the one genome-structure change: per-flyby direction coords
+ vector residual + bounds; ~80–120 lines, all in `dsm_leg.py`, all additive
behind a default-preserving flag).

**Cheapest decisive experiment.** ON 6.44Gg3 (the cleanest probe, fastest, and the
one with the sharpest prior negatives), with branches forced to the #157-winning
`(1,low)` on the long leg: run `dsm_chain_correct` with
`charge_flyby_continuity=True` from the descriptor-sourced symmetric seed. The
decisive question: does the residual vector drive to a point where BOTH
`Σ ΔV_DSM` is small AND every `flyby_dv ≈ 0` AND emerged V∞ ≈ 6.44/3.74? Gate
below.

**Honest failure mode — and why it is still decisive.** If, when the bend term is
charged, the corrector CANNOT simultaneously zero the DSM impulses and the
flyby_dv at *any* geometry near the anchor, that is the **strongest possible
negative**: it proves the sourced anchors do not host a bend-feasible ballistic (or
near-ballistic) closure in this topology — the "the basin may not exist" outcome
the chain keeps surfacing (cf. the memory blocker: "the 2-synodic E-M-E-E topology
does not host the published anchors as a feasible closed ballistic cycler"; the
Hughes-2014 corroboration that low-V∞ Mars needs broken-plane the authors exclude).
That is a publishable empty-set result with a *quantified* irreducible
flyby_dv floor — not inconclusive. Conversely if it closes, it is the **first
multi-arc closure** at the sourced anchors. Either way: decisive.

### Approach 4 — Continuation at the multi-arc level

**Idea.** Russell's circular→ephemeris homotopy (#158) reaches sourced basins by
*walking from a circular-coplanar closure*. The #158 blocker for S1L1 was "no
circular-coplanar closure to walk from" — because the *single-ellipse* genome
cannot produce one. So FIRST get the multi-arc genome to produce a
**circular-coplanar closure** for 6.44Gg3 / S1L1 (where #137's single-ellipse
failed), THEN feed that as the continuation seed.

**Files/functions.** Depends entirely on a prior approach succeeding on the
**circular** backend: `evaluate_dsm_chain` already accepts any `ephem`, so run
approach #2/#3 with `Ephemeris("circular")`. If that yields a closed coplanar
multi-arc chain, wire it as the seed to `continuation_correct`
(`continuation.py:335`) — but that corrector's inner solve is
`free_return_correct` (single-ellipse); it would need to accept a multi-arc inner
solve (`dsm_chain_correct`). That is a real change to `continuation.py`'s inner
hook.

**Effort.** HIGH, and **gated on a prerequisite**: continuation has nothing to do
until a circular-coplanar multi-arc closure EXISTS. Today it does not (that is the
whole point of #150/#153/#157). So this approach cannot be the *first* experiment —
it is strictly downstream of #2/#3 producing a coplanar closure.

**Honest failure mode.** If #2/#3 cannot produce a circular-coplanar multi-arc
closure (the diagnosis in §0.1 suggests the coplanar objective floors off-basin
too — #153 Config 2 floored *identically* with ToF frozen, i.e. structurally),
continuation has no seed and the approach is **vacuous** — it would return
"inconclusive" (no seed to walk), which violates the discipline. Therefore
**continuation is NOT a standalone candidate for the first experiment**; it is the
*reward* if #3 closes on the circular backend, at which point #158's machinery
walks it to DE440 for free.

### Approach 5 — Hybrid (symmetric seed + branch coord + explicit continuity residual)

**Idea.** #2 (seed) + #1 (branch) + #3 (objective) together: symmetric
descriptor-seed, branches forced to the per-arc rev count the descriptor implies
(g≈1.5 rev, G≈3 rev), and the explicit flyby-continuity+bend residual, all driven
by MBH.

**Files/functions.** Union of #1/#2/#3 above; no extra new surface.

**Effort.** MED (it is #3's code change + #1's threading + #2's seed helper — but
#3 is the load-bearing piece; #1 and #2 are cheap add-ons).

**Cheapest decisive experiment.** The full hybrid on 6.44Gg3, run *only if the #3
single-experiment shows the bend residual changes the landed basin* (i.e. #3 is
the gate for whether the hybrid is worth assembling).

**Honest failure mode.** Same empty-set possibility as #3 (the hybrid cannot
manufacture a basin that does not exist), but it is the configuration most likely
to *find* the basin if it does exist, because it fixes seed, family, AND objective
simultaneously. It is the natural **second** experiment, conditioned on #3's gate.

---

## 2. Ranking

Ranked by **decisiveness-per-unit-cost against the diagnosed mechanism** (§0.1):

1. **#3 Explicit flyby V∞-continuity + bend residual — MED, ★ recommended first.**
   The only approach that attacks the diagnosed mechanism (the objective does not
   reward the sourced bend-feasible low-V∞ point). Decisive either way: closes =
   first multi-arc closure; cannot close = quantified empty-set with an irreducible
   flyby_dv floor. Reuses `core/flyby.py` machinery already built.
2. **#2 Symmetric/period-locked seeding — LOW.** Cheap, and its "residual at the
   sourced symmetric seed" is a decisive datum (the multi-arc seed-at-truth probe).
   But it cannot fix the objective, so on the diagnosis it most likely returns a
   clean negative. Best run as the *seed* feeding #3, not alone.
3. **#1 Branch as a discrete coordinate — LOW.** Necessary ingredient; alone only
   moves the floor (as #157 already showed). Fold into the hybrid.
4. **#5 Hybrid — MED.** The most likely to *close* if a basin exists; assemble only
   after #3's gate fires positive-or-ambiguous.
5. **#4 Multi-arc continuation — HIGH, downstream only.** Vacuous until a
   circular-coplanar multi-arc closure exists; it is the payoff of #3-on-circular,
   not a first experiment.

---

## 3. Recommended first experiment + decisive gate + cost

**Experiment.** Add the explicit flyby V∞-continuity + bend-feasibility residual
(#3) behind a default-off `charge_flyby_continuity` flag in `evaluate_dsm_chain`,
free the intermediate-leg departure-V∞ direction, return a vector residual, and run
`dsm_chain_correct` on **6.44Gg3** (E-M-E-M-E, long leg branch forced `(1,low)`)
from the descriptor-sourced symmetric seed (a thin slice of #2 to provide the
seed). Single corrector run first (fast, ~20 s); then MBH if the single run shows
movement.

**Sourced acceptance target (EXPECTED side, never imposed):** emerged V∞
E = 6.44, M = 3.74 km/s (catalogue `russell-ch4-6.44Gg3`, Russell 2004 Table 4.13).

**Go/no-go gate (three-way, all decisive):**
- **GO / CLOSE (publishable positive):** `Σ ΔV_DSM < 0.1 km/s` AND every per-flyby
  `flyby_dv < 0.1 km/s` (bend-feasible) AND emerged V∞ within campaign tol
  (E 6.44 ±1.0, M 3.74 ±0.5). → First multi-arc closure; proceed to the hybrid (#5)
  to harden and to S1L1, then continuation (#4) to DE440.
- **NO-GO / EMPTY-SET (publishable negative):** the corrector converges (solver
  success) but the irreducible `min Σ(ΔV_DSM + flyby_dv)` stays bounded away from
  zero (> ~1 km/s) across a seed sweep AND emerged V∞ never approaches the anchors.
  → Record the quantified flyby_dv floor as evidence that the topology does not
  host a (near-)ballistic closure at the sourced anchors (consistent with the
  memory blocker + Hughes-2014). Stop; do not assemble the hybrid.
- **AMBIGUOUS (proceed to hybrid):** the bend residual *changes the landed basin*
  (emerged V∞ moves materially toward the anchors, floor drops below the #157
  26.9 km/s) but does not yet close. → The objective fix is working; assemble the
  full hybrid (#5: symmetric seed + branch enumeration + MBH) as the second
  experiment.

**Cost.** MED. Code: ~80–120 additive lines in `dsm_leg.py` (vector residual,
per-flyby direction coords, bounds), reusing `core/flyby.flyby_dv`. Wall: the
existing probes run in 17–62 s; the first single-run experiment is ~20 s, an MBH
follow-up ≤2 min. Total implementer effort: ~1 focused session under strict TDD.

---

## 4. Honest "could this be an empty set?" assessment

**Yes — and the design is built to PROVE it if so.** Three independent prior
results already point at empty-set for the *sourced anchors in this topology*:

- The memory blocker (2026-06-04, `correct_s1l1_twoarc.py`): closed ballistic
  E-M-E-E cyclers on DE440 exist, but S1L1's 5.65/3.05 is not among them — Mars V∞
  floors ~6.4 in feasible closures; the low Mars V∞ only appears in *discontinuous*
  geometry. "The 2-synodic E-M-E-E topology does not host the published anchors as
  a feasible closed ballistic cycler."
- Hughes-Edelman-Longuski 2014 (#147): pure-ballistic conic chains to low Mars V∞
  require a broken-plane maneuver the authors deliberately exclude — i.e. the
  low-V∞ family is *not* ballistic-representable without DSMs (which is exactly what
  the multi-arc DSM genome adds — so the open question is whether the DSM budget at
  the anchor is *small* or *irreducibly large*).
- #135 seed-at-truth: the sourced geometry was not a residual-zero point of the
  single-ellipse genome (3.2–37.5 km/s at truth). Approach #3's seed-at-symmetric
  probe is the multi-arc analogue and will return the same kind of decisive number.

What approach #3 adds that none of the prior passes had: an objective whose minimum
**is** the bend-feasible low-V∞ point. If that objective still cannot reach zero,
the irreducible `flyby_dv` floor is a *quantified* statement of the empty-set — a
publishable negative ("the sourced anchors require ≥ X km/s of irreducible
flyby/DSM work in this topology; they are not a ballistic multi-arc closure"). The
note that the three S1L1 rows carry three model-dependent V∞ (5.65/3.05, 4.99/5.10,
4.7/5.0) means the *target itself is framing-ambiguous* — which is one more reason
the honest deliverable may be "no ballistic closure at these anchors" rather than a
closure, and #3 is the experiment that decides it cleanly.

---

## 5. Files referenced (all read-only this pass)

- `src/cyclerfinder/search/dsm_leg.py` — the multi-arc genome + chained evaluator
  (`evaluate_dsm_chain:357`, `dsm_chain_correct:614`, `dsm_chain_decision_vector:582`,
  `make_dsm_chain_step:701`); the forced-continuity line is `:448`, the scalar
  objective `:474`.
- `src/cyclerfinder/core/flyby.py` — `flyby_dv:356` (Eq.5.5 surrogate, returns 0
  when ballistic-feasible), `dv_powered_flyby_periapsis:117` (Oberth),
  `is_ballistic_feasible:318`, `max_bend:40` — the continuity/bend pieces #3 wires.
- `src/cyclerfinder/search/free_return.py` — `seed_ae_from_aphelion_transit:247`
  (the #137 seed pattern #2 mirrors).
- `src/cyclerfinder/search/continuation.py` — `continuation_correct:335` (the #158
  homotopy #4 would consume a multi-arc seed into).
- `src/cyclerfinder/search/mbh.py` — the generic basin-hopper (unchanged).
- `src/cyclerfinder/core/lambert.py` — `lambert:511` multi-rev API (already wired).
- Notes: the #150/#153/#157 probes, #145/#156 MBH, #158 continuation, and
  `multi-arc-classification.md` §7/§12; memory `project_s1l1_realeph_closure_blocker.md`.
