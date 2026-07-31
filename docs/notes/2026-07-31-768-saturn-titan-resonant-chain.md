# `#768` — Saturn-Titan periodic 3:4<->6:5 "resonant chain" (Vaquero 2013 Fig. 4.9-4.12)

**Task:** `#768`, the culmination of the Saturn-Titan thread `#764` -> `#765` -> `#767` -> `#768`.
With `#765` confirming the 3:4 resonant orbit to near-machine precision and leaving 6:5 as an
honest, small, well-characterized eigenvalue near-miss, and `#767` building four independent,
Newton-converged homoclinic self-connections of the 3:4 orbit, this task's job (per the dispatch
note) was: (1) make an explicit, evidence-based call on whether 6:5 is ready to use for connection
work; (2) if so, build whatever connection(s) Vaquero 2013 Sec. 4.3.1 actually describes as the
periodic "resonant chain" cycling between the 3:4 and 6:5 resonances (Fig. 4.10, continued in
Jacobi constant in Fig. 4.11); (3) attempt to confirm or refute her own falsifiable claim that the
family "ends for a value of Jacobi constant C < 3.01400" (Fig. 4.12).

**Sources read directly this task**: Vaquero 2013 Sec. 4.3.1 pp.104-118 (the FULL "Design of
Planar Transfers" and "Design of Planar Periodic Resonant Chains" subsections, not just the
summary already excerpted in `#765`'s own note) — `cyclers_pdf/papers/vaquero-2013-...-phd.txt`,
lines ~4010-4335 — read in full this task specifically to resolve the geometry question the
dispatch note flagged ("check the actual geometry... don't assume"). Also re-read
`docs/notes/2026-07-29-765-saturn-titan-resonant-families-vaquero-gate.md`,
`docs/notes/2026-07-31-767-saturn-titan-homoclinic-connection.md`, and the two Jovian-chain
reviewer-verdict precedents (`docs/notes/2026-07-28-755-...md`,
`docs/notes/2026-07-28-758-...md`) and the Jovian heteroclinic-connection precedent
(`docs/notes/2026-07-28-759-...md`) named in the dispatch note.

---

## Verdict (read this first)

1. **Step 1 — 6:5 readiness: proceed, on a narrower and better-supported basis than the dispatch
   note anticipated.** Re-reading the thesis directly (see Step 1 below) shows the "resonant
   chain" needs 6:5's own IC *location* only, not its eigenvalue/manifold structure — and that
   narrower claim is well-supported (new Radau cross-check this task, plus `#765`'s own existing
   evidence). `#765`'s own eigenvalue gate FAIL is **not** overridden or re-litigated here.
2. **Step 2 — a genuine, independently-found homoclinic self-connection of the 3:4 orbit whose
   own Poincare-map crossing sits markedly closer to 6:5's own fixed point than anything `#767`
   itself reported** (self-consistency evidence directly corroborating Vaquero's own qualitative
   Fig. 4.9 description). The further step — correcting this excursion into an exact NEW periodic
   "chain" orbit (Fig. 4.10), per Vaquero's own described single-shooting method — is an **honest
   PARTIAL/NEGATIVE**: a bounded, good-faith STM-based 2-D Newton attempt makes real, ~40x
   monotonic progress (residual 0.253 -> 0.0063) but stalls before reaching machine precision,
   consistent with (not contradicting) `#759`'s own documented Table-3 stall for a comparably
   severe manifold-sensitivity regime.
3. **Step 3 — explicitly out of scope**, per the dispatch note's own escape valve. Registered as
   follow-up `#773`/`#774` below.

No catalogue writeback, literature-novelty check, or "novel" claim is made anywhere in this note —
all three are explicitly out of scope per the dispatch note and flagged as next steps at the end.

---

## Step 1: what does the thesis actually need from 6:5? (re-read before assuming)

The dispatch note assumed the chain might need `Wu(3:4) ∩ Ws(6:5)` and `Wu(6:5) ∩ Ws(3:4)` (two
heteroclinic legs), mirroring the Jovian `#759` construction. **Reading Vaquero 2013 pp.113-117
directly (Sec. 4.3.1, "Design of Planar Periodic Resonant Chains") shows this is NOT what she
does.** Her own words (p.114-115, quoted at length because this is the load-bearing finding):

> "As an application of dynamical systems theory, consider a homoclinic connection that
> asymptotically departs and approaches the 3:4 resonant orbit... A particular subset of the
> manifolds from the 3:4 resonance travel to the interior region, shadowing the 6:5 resonant
> orbit. To highlight the relationship between the two unstable resonant orbits, **the
> intersection on the Poincare map for the homoclinic connection is selected to be near the fixed
> point corresponding to the 6:5 resonant orbit**, as illustrated in Figure 4.9. The path that
> results from propagating this intersection point on the map forward and backward in time
> naturally follows both resonant orbits... To produce a resonant chain, or a cycle between
> resonances, it is necessary to numerically correct this path via a single shooting scheme to
> obtain a periodic orbit that shadows the invariant manifolds associated with the two resonant
> orbits. The resulting trajectory is plotted in Figure 4.10(a) and its periodicity is represented
> by the two perpendicular crossings."

This is a **homoclinic self-connection of the 3:4 orbit alone** (`Wu(3:4) ∩ Ws(3:4)`) — exactly
`#767`'s own machinery — with ONE extra selection criterion (the crossing must sit *near* 6:5's
own fixed point on the map) and ONE extra step this task chain has not yet built (a further
single-shooting periodicity correction of the resulting excursion into a genuinely new, closed
periodic orbit). It is **not** a two-orbit heteroclinic construction at all. This confirms the
dispatch note's own third alternative ("the SAME homoclinic self-connections already built plus
some other mechanism").

**Consequence for the 6:5-readiness question**: this construction only needs 6:5's own **fixed
point location** — `(x0, ydot0)`, a static reference point used purely to *select which of
`#767`'s own homoclinic crossings to use* — never 6:5's own manifolds, eigenvalue, or stability
classification. That is a strictly narrower and better-supported claim than "6:5 is the same real
physical orbit for eigenvalue/manifold purposes," which is what `#765`'s own honest FAIL is about.

### Additional corroborating evidence gathered this task

1. **Independent Radau cross-check of 6:5's own IC** (not run by `#765`; done this task,
   `crosscheck_periodic`-style, `rtol=atol=1e-11`, one full period, DOP853-corrected IC
   re-propagated with the independent implicit-RK Radau integrator):

   | Orbit | Radau closure distance | Jacobi drift over 1 period |
   |---|---|---|
   | 3:4 | `4.719e-10` | `8.882e-16` |
   | 6:5 | `6.814e-11` | `6.128e-14` |

   6:5's own IC is **at least as tightly self-consistent as 3:4's own** under a completely
   independent integrator — direct evidence 6:5's own IC/period is a genuine, well-converged
   periodic solution, not an artifact of the corrector that happens to also miss the eigenvalue
   gate. (3:4's own eigenvalue gate PASSES, so this comparison shows the Radau-closure axis is not
   simply tracking the eigenvalue-gate outcome — it is an independent quality signal.)

2. **`#765`'s own already-existing basin-robustness result stands unchanged**: 11 seeds across a
   +/-2e-4 window around 6:5's own Table 4.1 seed all converge to the identical eigenvalue — not
   an isolated numerical fluke.

3. **The thesis's own text confirms 6:5's IC was itself derived from 3:4's own Poincare map**
   (p.108-109): "The interior resonant orbit in Figure 4.6(a) is computed from an initial guess
   obtained from the Poincare map associated with the 3:4 unstable resonant orbit in Figure 4.5."
   This is exactly consistent with — and a plausible *cause* of — this task's own independent
   numerical finding (Step 2 below) that a genuine homoclinic self-crossing of 3:4 sits
   suspiciously close to 6:5's own fixed point: 6:5's own location and 3:4's own manifold
   structure are not coincidentally near each other, they are directly related by construction in
   the source thesis itself.

### The call

**Judged NOT a genuinely ambiguous close call for this specific, narrower IC-location use** (in
contrast to `#758`'s own 5:6-LO eigenvalue-confirmation call, which WAS genuinely close and needed
an explicit reviewer ruling). The eigenvalue-gate FAIL (`#765`, 2.34e-3) is real, unchanged, and
NOT overridden here — but it is simply not the axis this connection-selection step depends on.
Proceeding to Step 2 using 6:5's own `(x0, ydot0)` as a fixed reference point.

---

## Step 2a: locating the near-6:5 crossing among `#767`'s own homoclinic family

`#765`'s recovered 6:5 candidate gives the target fixed point (on this module's own
`{y=0, ydot>0}` section, unrestricted in `x`):

```
target_65 = (x=0.9347726861768341, xdot=0.0)
```

`#767`'s own originally-reported reflection-symmetric MIRROR pair (`branch_u=branch_s=-1`,
`(k_u,k_s)=(4,5)/(5,4)`, crossing `x=0.844031, xdot=∓0.114813`) sits at Euclidean distance
`0.1463` from this target — already, by construction, the closest of `#767`'s own four originally
reported hits (the on-axis PRIMARY and off-axis TERTIARY hits sit at distance `>2.0`, an order of
magnitude further).

**This task independently re-scanned the same `(branch_u, branch_s)=(-1,-1)` combination** (same
tolerances as `#767`'s own refined pass: `scan_n=12, tol=1e-9, max_iter=60, fd_step=1e-7`, but a
fresh scan-grid start rather than pinning from `#767`'s own light-triage seed) and found a
**different, genuine local root of the same `(k_u,k_s)=(4,5)/(5,4)` residual equation** — a real,
independently-converged mirror pair markedly CLOSER to 6:5's own fixed point:

| branch_u | branch_s | k_u | k_s | residual | crossing `(x, xdot)` | ghost margin | dist to 6:5 |
|---|---|---|---|---|---|---|---|
| -1 | -1 | 4 | 5 | `9.628e-10` | `(0.914073, -0.091737)` | `148x` | `0.09404` |
| -1 | -1 | 5 | 4 | `2.883e-10` | `(0.914073, +0.091737)` | `148x` | `0.09404` |

This is **36% closer to 6:5's own fixed point than `#767`'s own reported MIRROR pair** (`0.094` vs
`0.146`), while carrying the same class of evidence quality: Newton residual `<1e-9`, ghost-guard
margin `148x` the `GHOST_GUARD_DELTA=1e-3` threshold (a real, non-delicate margin), and this
task's own forward/backward re-approach check (`stc.homoclinic_reapproach_check`, `rtol=1e-13,
atol=1e-14` throughout): `backward_distance=1.27e-9`, `forward_distance=0.0053` — tighter than
`#767`'s own PRIMARY hit's own forward re-approach (`0.0112`).

That this task's own fresh grid-scan start converges to a genuinely different (and closer) root
than `#767`'s own pinned-seed refinement, for the nominally SAME `(branch, k)` label, is expected
and not a contradiction: the coarse scan-grid start only guarantees convergence to *some* local
root of the residual equation, and this Poincare-map curve visibly has multiple genuine zero
crossings within the same `(k_u, k_s)` window (a documented feature of this whole task chain's own
manifold-return maps, not specific to this candidate). Both are legitimate, independently-verified
homoclinic self-intersections of the 3:4 orbit.

**This is the strongest available self-consistency corroboration of Vaquero's own qualitative
Fig. 4.9 description** ("the intersection...is selected to be near the fixed point corresponding
to the 6:5 resonant orbit") that this task chain's own machinery can produce, given no digit-grade
Table exists for Fig. 4.9 itself (same honesty caveat as `#767`'s own note).

Elapsed transit times for this candidate (needed for Step 2b below): `t_u = 39.805` nondim,
`t_s = -70.695` nondim (i.e. `|t_s| = 70.695`); total loop time `t_u + |t_s| = 110.4996` nondim
time units (`~4.228` periods of the 3:4 orbit's own `T=26.1408`).

## Step 2b: attempting the periodicity correction (Fig. 4.10) — honest partial/negative

Vaquero's own method: "it is necessary to numerically correct this path via a single shooting
scheme to obtain a periodic orbit that shadows the invariant manifolds associated with the two
resonant orbits... its periodicity is represented by the two perpendicular crossings." This means
the target chain orbit is itself a symmetric periodic orbit (same class as 3:4/6:5 themselves),
just with a much longer period (`~110` nondim vs `3:4`'s own `26.14`).

**Why the existing `correct_connection`/`heteroclinic_cycle` machinery doesn't apply here**: that
machinery corrects a *manifold-phase* parametrization (`tau_u, tau_s`) to match TWO separately-
seeded manifold legs at a section crossing — it is not built for a genuine Poincare-map
fixed-point search (a single trajectory returning to itself after one loop). This task built a
NEW, purpose-specific STM-based 2-D Newton corrector
(`saturn_titan_resonant_connections.attempt_chain_closure`, plus its private helpers
`_chain_ydot`, `_chain_state0`, `_chain_crossings`, `_chain_map_step`):

* **Free parameters**: `(x, xdot)` on the `{y=0}` section (2 unknowns), with `ydot` slaved to the
  Jacobi constant (generalizing
  `cr3bp_periodic.ydot0_from_jacobi`/`correct_symmetric_fixed_jacobi`'s own `xdot=0`-only case to
  general `xdot`).
* **Analytic (STM-based) Jacobian**, not finite-difference: essential here, not a style choice —
  the compounded growth over the `~110` nondim-time / `~4.2`-period loop at `|lambda|~2129.8` is
  `exp(ln(2129.8)/26.14 * 110.5) ~ 1.2e14`, right at the edge of double-precision's dynamic range;
  an FD-based Jacobian at any usable step size would be dominated by cancellation noise at this
  compounded growth (confirmed by inspection, not assumed — this is why `correct_connection`'s own
  FD-based Jacobian, which works fine for `#767`'s own single-crossing (`~1-3` period) connections,
  was not reused for this multi-period problem).
* **Fixed crossing index**, determined once from the seed (nearest `t_target=110.4996`) and held
  fixed across iterations — mirroring `correct_symmetric_fixed_jacobi`'s own `half_crossings`
  discipline, to keep Newton on a single continuous branch.
* **Damped backtracking line search** (accept a step only if the residual norm decreases,
  `max_backtrack` halvings) — this matters: an EARLIER, cruder exploratory version of this
  corrector without backtracking (undamped, `max_step`-clipped only) took a first Newton step that
  left the seed's local linear regime entirely — the qualifying-`{y=0}`-crossing count within the
  same horizon exploded from `16` (at the seed) to `290` after just one step, meaning the "fixed"
  crossing index silently jumped onto an entirely different, much-shorter-period branch (not the
  desired `~110`-time-unit chain). This diagnostic is real and informative (it is direct, measured
  evidence of just how astronomically ill-conditioned this map is at this compounded instability),
  but is NOT the behaviour of the function actually shipped and tested below.

**With backtracking** (the shipped `attempt_chain_closure`, seeded at `node`'s own IC
`(x0=1.0301663, xdot=0)`, `t_target=110.4996`):

| `max_iter` | `max_backtrack` | final residual | converged | notes |
|---|---|---|---|---|
| 1 | 8 | `0.2534298` | No | seed residual only (no step taken) |
| 2 | 3 | `0.0128414` | No | real progress, still exhausted |
| 8 | 6 | `0.0063392` | No | `7` iterations run; line search exhausted (no improving step found) |

**This is a genuine, real, ~40x reduction in residual (`0.253` -> `0.0063`) via a properly damped
Newton scheme with an analytic STM Jacobian — not a wild divergence, and not a forced result** —
but it plateaus and stalls before reaching the `1e-9` convergence tolerance. This mirrors, almost
exactly in character, `#759`'s own documented finding for `Ws(5:6-LO)`'s heteroclinic connection:
*"the residual plateaus around `2.6e-3`-`3.5e-2` without converging... diagnostic of an
ill-conditioned Jacobian in this neighborhood."* Here the mechanism is somewhat different (compound
multi-period growth rather than intra-period fractal manifold sensitivity), but the qualitative
outcome — a genuine, well-instrumented Newton stall, not evidence against the chain's existence —
is the same class of honest limitation.

**This is NOT evidence the periodic chain doesn't exist.** Vaquero's own thesis demonstrates it
does (Fig. 4.10, with a plotted trajectory and an explicit continuation family in Fig. 4.11). It is
evidence that a plain single-shooting 2-D Newton scheme, from a physically-motivated but not
highly-refined seed, is not by itself equal to closing a `~4.2`-period loop at this magnitude of
compounded instability within this task's own effort budget — consistent with (not contradicting)
this whole task chain's repeated experience that severely unstable (`|lambda| > ~1000`) multi-leg
constructions break naive Newton correctors (`#759`).

---

## Step 3: the `C < 3.01400` termination claim — explicitly out of scope

Confirming or refuting Vaquero's own falsifiable claim ("it is suspected that this family of
periodic resonant chains ends for a value of Jacobi constant `C < 3.01400`", Fig. 4.12) requires,
at minimum: (a) a CONVERGED chain orbit at `C=3.01` as a starting point (Step 2b's own honest
stall means this is not yet in hand), then (b) a continuation-in-`C` campaign tracking that family
up toward `3.01400` and checking whether the underlying `Wu(3:4)`/`Ws(3:4)` manifold intersection
near 6:5's own fixed point survives or vanishes at each step — a materially larger undertaking than
this task's own scope, and explicitly gated on Step 2b's own unresolved closure. Per the dispatch
note's own explicit escape valve ("if it's too large for this task's own scope, explicitly say so
and register it as a follow-up"), this is registered as follow-up work below, not attempted in a
rushed/forced form.

---

## Follow-up work registered (not attempted this task)

* **`#773`** — close Step 2b's own honest Newton stall: either (a) a genuinely better seed (e.g.
  digitizing Vaquero's own Fig. 4.10(a) plotted trajectory for a much closer starting `(x,xdot)`),
  or (b) a multiple-shooting corrector (many patch points along the `~110`-nondim-time loop, each
  individually well-conditioned) rather than single-shooting — the standard professional-grade fix
  for exactly this class of severely-unstable multi-revolution closure problem.
* **`#774`** — once `#773` produces a converged chain orbit at `C=3.01`, run the continuation-in-`C`
  campaign to confirm or refute Vaquero's own `C < 3.01400` termination claim (Step 3 above).

## Catalogue/novelty: explicitly out of scope, flagged for the coordinating session

No chain orbit was fully converged this task, so no catalogue-eligibility question actually arises
yet. If/when `#773` produces a genuinely converged periodic chain orbit, per the dispatch note's
own instruction: that object may be independently catalogue-eligible (a genuine published
Saturn-Titan cycler-like object) but is subject to the mandatory `literature_check.py` novelty gate
and a catalogue-writeback design decision, both explicitly out of scope for both this task and
`#773`/`#774` as registered — flagged here as the natural next step once a chain orbit converges.

---

## Code delivered

* `src/cyclerfinder/search/saturn_titan_resonant_connections.py` (extended, not a new module —
  the `#768` work is a direct, small extension of `#767`'s own machinery, per the Step 1 finding
  that this reuses the SAME homoclinic self-connection mechanism):
  `resonant_chain_target_point` (6:5's own IC-only reference point + its honest gate row),
  `ChainProximityCandidate`/`rank_by_proximity_to_65` (re-ranks `#767`'s own `find_homoclinic`
  hits by distance to 6:5's own fixed point instead of residual), `ChainClosureResult`/
  `attempt_chain_closure` (the STM-based 2-D Poincare-map Newton corrector, with its own honest
  non-convergence documented in both the docstring and the test suite) and the private helpers
  `_chain_ydot`, `_chain_state0`, `_chain_crossings`, `_chain_map_step`. No changes to
  `jovian_resonant_connections.py`, `heteroclinic_cycle.py`, `saturn_titan_resonant_families.py`,
  or `cr3bp_periodic.py` (`_ubar_grad_x_at_axis` is imported and reused unchanged).
* `tests/search/test_saturn_titan_resonant_connections.py` (extended, 28 tests total, up from 23):
  6 new tests covering `resonant_chain_target_point`'s own honest gate row, the closer near-6:5
  homoclinic candidate's own full convergence/ghost-margin evidence, `attempt_chain_closure`'s own
  exact reproducible seed residual and its bounded, honest (non-forced) partial progress, and the
  Jacobi-radicand honesty check. None marked `@pytest.mark.slow`.

## Verification

* `uv run ruff check` / `ruff format --check` on both changed files: clean.
* `uv run mypy --strict` on both changed files: clean.
* `uv run mypy src tests` (canonical full invocation): clean, 827 source files.
* `uv run pytest tests/search/test_saturn_titan_resonant_connections.py -q`: 28/28 pass (~112s).
* `uv run pytest tests/search/test_saturn_titan_resonant_connections.py
  tests/search/test_saturn_titan_resonant_families.py -q`: 55/55 pass.
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: run before the `OUTSTANDING.md`
  update commit (see commit history for pass/fail status recorded at commit time).

## Net effect on `#768`

**Step 1: resolved, proceed** (narrower, IC-location-only basis; `#765`'s own eigenvalue gate FAIL
stands unchanged). **Step 2: genuine positive self-consistency corroboration of Vaquero's own
Fig. 4.9 geometry** (a homoclinic self-connection of 3:4 whose crossing sits closer to 6:5's own
fixed point than any previously-found candidate), **plus an honest, well-characterized,
NOT-forced partial/negative on the further periodicity correction into an exact new chain orbit**
(Fig. 4.10) — real Newton progress, genuine stall, not silent divergence. **Step 3: explicitly out
of scope**, registered as `#773`/`#774`. This task chain (`#764`->`#765`->`#767`->`#768`) has now
built the full available self-consistency evidence base for the Saturn-Titan resonant-chain
question that this project's existing single-shooting machinery can reach; closing the loop fully
requires the follow-up work registered above.
