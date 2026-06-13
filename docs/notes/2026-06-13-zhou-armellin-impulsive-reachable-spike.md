# Zhou-Armellin single-impulse reachable-set SPIKE + NRHO cross-check + merge assessment (#239)

**Date:** 2026-06-13
**Method source:** X. Zhou, R. Armellin, D. Qiao, X. Li, *Single-Impulse
Reachable Set in Arbitrary Dynamics Using Polynomials*, arXiv:2502.11280v1
[astro-ph.IM], 16 Feb 2025.
**Mining note (transcription source):**
`docs/notes/2026-06-07-zhou-2025-da-reachable-sets-mining.md`.
**Complements:** the Braik-Ross energy-PRESERVING reachable-set scorer (#236/#247),
`src/cyclerfinder/search/reachable_network.py`,
`docs/notes/2026-06-13-braik-ross-reachable-set-scorer-results.md`,
`docs/notes/2026-06-13-reachable-scorer-ungate.md`.
**Scope:** SPIKE / method prototype. **No catalogue writeback. No validation-level
change.** Sourced / independent-recompute cross-checks only.

---

## TL;DR

- Built `src/cyclerfinder/search/reachable_impulsive.py`: a **single-impulse**
  CR3BP reachable-set prototype — bounded impulse of arbitrary direction on the
  max-magnitude sphere (Zhou Eqs. 4-6), propagated with the existing float
  `core.cr3bp.propagate`, projected onto the plane orthogonal to the nominal
  velocity (Zhou Eq. 8), characterized by convex-hull footprint area/extent.
- It is the **sampling** realization of the method, **not** the differential-
  algebra (DA) realization. The DA/polynomial machinery (Taylor map +
  ADS + envelope root-find + local-poly, Eqs. 7-51) is the >84% CPU optimization
  and needs a DA-evaluable force model (DACE/DACEyPy) our pure-float stack does
  not host (mining note Sec. 3). The spike computes the **same geometry by brute
  force** — the method's geometry without the method's speed.
- **NRHO cross-check (independent recompute, clean PASS).** Recovered the
  Earth-Moon L2 southern **9:2 NRHO** (sourced Gateway/Lee-2019 seed, re-corrected
  under Zhou's μ = 0.0121505839): perilune **3247 km**, apolune/perilune aspect
  **21.9**, T = **1.511 TU** (6.56 d), C_J = **3.0465** — a genuine near-rectilinear
  halo. The independent Monte-Carlo cloud of feasible impulses is **fully
  contained** by the spike's grid footprint, with Zhou's error index
  **P = d_max²/S_RS = 0.000%** on a 13×25 grid — meeting the paper's <0.1%
  acceptance bar (mining note Sec. 4).
- **Energy contrast demonstrated.** A 50 m/s impulse shifts C_J by **−7.2e-3**
  (`test_apply_impulse_changes_jacobi_constant`). This is the defining difference
  from Braik-Ross, whose heading-change maneuver holds C_J fixed by construction.
- **Merge assessment: YES, the two compose into one accessibility framework** —
  Braik-Ross is the *on-manifold* (energy-preserving, screening) layer; Zhou is
  the *cross-manifold* (energy-changing, budgeted) layer. Concrete merge path
  below. Proposed only — `reachable_network.py` was **not** edited (concurrent
  agents / collision avoidance).

---

## What the spike implements

`src/cyclerfinder/search/reachable_impulsive.py`:

| function | role | Zhou ref |
|---|---|---|
| `impulse_vector(dv_mag, α, β)` | impulse on the sphere in the local velocity frame; `(0,0)` = prograde | Eqs. 4-6 |
| `velocity_frame(r, v)` | orthonormal triad `[e_v, e_n, e_h]` from `v`, `h=r×v`, `h×v` | Eq. 8 |
| `apply_impulse(state, dv_local)` | rotate impulse into the rotating frame, add to velocity (energy-changing) | Eq. 4 |
| `reachable_cloud(...)` | sample `(α, β)` grid on the max sphere, propagate each to `t_f`, project to the auxiliary plane | Eqs. 4-11 |
| `monte_carlo_cloud(...)` | independent random-impulse truth cloud (sphere or full ball), same projection | Sec. V (MC validation) |
| `footprint_metrics(...)` | convex-hull area `S_RS`, bounding-box extent, centroid | Eq. 55 `S_RS` |
| `containment_crosscheck(...)` | contained fraction, `d_max`, error index `P = d_max²/S_RS` | Eq. 55 |

**Faithful to the paper:** the impulse model (arbitrary direction, fixed epoch,
max-sphere boundary), the auxiliary-plane projection, and the `P` accuracy gauge.
**Deliberately NOT reproduced (DA-specific, out of spike scope):** the DA Taylor
map, Automatic Domain Splitting, the envelope-equation root-find (Eqs. 34-45),
the local-poly speedup (Eqs. 46-51), and the alpha-shape non-convex boundary. The
spike's convex-hull footprint is an **over-approximation** of the true (possibly
non-convex) RS envelope — adequate for an evidence band, not a tight boundary.

---

## NRHO cross-check (independent recompute — the headline result)

The Zhou method invariant we reproduce (mining note Sec. 4): **the MC cloud lies
inside the RS boundary, with `P = d_max²/S_RS` well under 0.1%.** No published
state vector is asserted (the paper prints none for the RS members); the MC cloud
is the independent truth set. This is an independent-recompute cross-check, not a
circular golden.

### Recovered 9:2 NRHO (sourced seed)

Seed: the widely published EM L2 southern 9:2 NRHO IC (NASA Gateway baseline,
Lee 2019), re-corrected under Zhou's μ with `cr3bp_periodic.correct_periodic`:

| quantity | value |
|---|---|
| period T | 1.511 TU (6.56 d) |
| Jacobi C_J | 3.0465 |
| perilune | 0.00845 LU = **3247 km** |
| apolune | 0.1853 LU |
| aspect (apolune/perilune) | **21.9** (near-rectilinear) |
| z-excursion | southern (z_min −0.182, z_max +0.008) |

These match the canonical 9:2 NRHO (T ≈ 1.51 TU, perilune ~3000 km, C_J ≈ 3.05),
confirming we recovered the right orbit. (Closure residual < 1e-10.)

### Reachable-footprint containment (seed at apolune, dv_max = 10 m/s = Zhou Table 4)

`dv_max` nondimensionalized with VU = 1023.16 m/s (EM length/time units).

| t_f | grid | footprint S_RS [LU²] | extent [km] | MC contained | d_max | **P = d_max²/S_RS** |
|---|---|---|---|---|---|---|
| 0.25 T | 13×25 | 3.98e-5 | 1386 × 1354 | **100%** | 0.0 | **0.000%** |
| 0.25 T (×5 dv = 50 m/s) | 13×25 | 9.95e-4 | 6931 × 6771 | **100%** | 0.0 | **0.000%** |

The footprint scales ~quadratically in area with dv_max (a 5× impulse gives ~25×
area: 3.98e-5 → 9.95e-4), as expected for a linear-dominated reachable map at this
short arc. The 10 m/s single-burn footprint at the NRHO apolune is **~1400 km
across** over a quarter-period — a directly interpretable "what a 10 m/s TCM
buys" number.

### Grid-resolution / nonlinearity stress (seed at apolune, t_f = 0.5 T → reaches perilune)

Zhou's 9:2 case is `t_f = 0.5 T` reaching perilune (highly nonlinear; the regime
where their partial-map inversion *overestimates* the RS and Newton's iteration is
required). The spike has no inversion to fail, but a **coarse grid under-covers
the curved boundary** there:

| grid | footprint S_RS [LU²] | MC contained | d_max | **P** |
|---|---|---|---|---|
| 7×13 (coarse) | 7.45e-5 | 98.5% | 3.43e-4 | **0.158%** |
| 13×25 | 8.24e-5 | 99.75% | 1.90e-5 | 0.00044% |
| 21×41 (fine) | 8.41e-5 | **100%** | 0.0 | **0.000%** |

This reproduces Zhou's **accuracy-vs-resolution tradeoff** (mining note Sec. 4:
tightening the ADS threshold cuts relative error): refining the angle grid
monotonically improves containment and area convergence. The coarse-grid 0.158%
miss near perilune is the spike's analogue of the paper's perilune-nonlinearity
stress — an honest, discriminating cross-check (not a trivially-passing one).

**Cross-check verdict: CLEAN PASS.** The spike reproduces the Zhou method
invariant (MC-cloud containment, `P` < 0.1% at adequate resolution) on a genuine
9:2 NRHO, and exhibits the documented resolution/nonlinearity sensitivity. It does
NOT reproduce any *specific* Zhou table number (different orbit instance, no
published state vector, convex-hull vs alpha-shape boundary) — those are not
golden-eligible for us per the mining note, and the discipline forbids asserting
our own first output as golden.

---

## Does it cross-check against Zhou-Armellin specifically?

Partially, and honestly so. The paper's **method invariants** (mining note Sec.
4) are reproduced: (b) MC cloud fully contained by the footprint with `P` well
under 0.1%; the resolution-error tradeoff. The paper's invariant (a) (local-poly
vs exact-solve envelope agreement ~1e-12) and (c) (inversion-vs-Newton perilune
overestimate) are **not testable here** — they are properties of the DA envelope
layer the spike deliberately omits. The paper's *numeric* table values (12.32 s
CPU, 84.17% reduction, P = 0.0658% over 100 epochs) are DA-implementation and
orbit-instance specific and were never golden-eligible. So: the **geometry and the
acceptance invariant cross-check cleanly**; the **DA speedup and its
failure-modes are out of scope** and flagged as such.

---

## Merge assessment with the Braik-Ross #236 scorer

**Verdict: they merge into one accessibility framework — they are the two
complementary maneuver classes.** This is the cleanest part of the result.

| axis | Braik-Ross (`reachable_network.py`, #236/#247) | Zhou-Armellin (this spike) |
|---|---|---|
| maneuver | energy-**preserving** heading change (pure rotation at fixed speed, Eq. 26 `2v·sin(\|δ\|/2)`) | general single **impulse** (arbitrary direction + magnitude) |
| energy | stays on ONE `C_J` manifold | **changes** `C_J` (cross-manifold) |
| primitive | reachable set per family at fixed energy → N×N proxy graph | reachable footprint per seed/epoch on the auxiliary plane |
| output | family-accessibility ranking (screening) | budgeted reachable footprint ("what ≤X m/s buys from epoch t") |
| status | GATED (C32-dominance gate did not reproduce on the recoverable subset) | spike, cross-checked |

The two answer **disjoint, composable questions**: Braik-Ross answers *"which
families are well-connected at this energy by free heading changes?"* (a cheap,
necessary-not-sufficient pre-screen); Zhou answers *"given a real ΔV budget at a
real epoch, what set of states can I reach, crossing energy levels?"* (a budgeted
footprint). A heading change is the **zero-energy-change special case** of a
general impulse — Braik-Ross is literally the on-manifold restriction of the Zhou
impulse set.

### Concrete merge path (proposed; NOT implemented to avoid collision)

A unified `accessibility` layer with one maneuver abstraction:

1. **Common maneuver interface.** Both already share `core.cr3bp.propagate`, the
   `(r, v)` → velocity-frame triad (Braik-Ross `heading`/`velocity_from_heading`
   ≈ this spike's `velocity_frame`), and a voxel/footprint reachability primitive.
   Define `Maneuver.apply(state) -> state'` with two implementations:
   `HeadingChange(δ)` (energy-preserving; Braik-Ross Eq. 26 cost) and
   `Impulse(dv_mag, α, β)` (energy-changing; this spike). `dv_turn` for the former
   and `‖dv‖` for the latter are both ΔV costs in the same VU units (both modules
   already use VU = 1023.16 m/s), so they are directly comparable on one budget.

2. **Two-tier accessibility.** Use Braik-Ross as the **cheap pre-screen** (which
   family pairs are even plausibly connected at a reference energy) and the Zhou
   impulsive footprint as the **budgeted confirm** on the survivors (does the
   target node lie inside the ≤dv_max single-burn footprint from epoch t?). This
   is exactly the "necessary screen → budgeted feasibility" two-tier the
   Braik-Ross note already frames its proxy ΔV as the first half of.

3. **Shared voxel/hull grid.** Braik-Ross logs reachability into an `(x, y, θ)`
   voxel grid; the Zhou footprint is a `(x_p, z_p)` plane hull. A merged layer
   can voxel-log the impulsive cloud too (3-D `(x, y, z)` or the 2-D plane), so
   "is the target reachable" is one grid-overlap query regardless of maneuver
   class. The Braik-Ross time-reversal mirror (Eq. 14) is maneuver-agnostic and
   reused for the backward set of either class.

4. **Where each is authoritative.** Braik-Ross is planar-only / single-energy /
   screening (its own stated limits); Zhou is 3-D / cross-energy / short-arc
   (its own stated limits — expansion validity degrades over long coasts, mining
   note Sec. 6). The merged framework should route by regime: heading-change
   screening for the *coarse family graph at one energy*; impulsive footprints for
   the *budgeted feasibility of a specific node at a specific epoch* (e.g. a TCM
   on a cycler leg, mining note Sec. 5).

**Why not implemented now:** (i) task scope is a spike + assessment, not a
production merge; (ii) `reachable_network.py` is read-only for this run (concurrent
agents, collision avoidance); (iii) the Braik-Ross scorer is itself **GATED** (its
C32-dominance validation gate has not reproduced — #247), so wiring it into a
production accessibility path would be premature. The merge is *architecturally
clean and low-friction* (shared propagator, shared VU, shared grid primitive) but
should wait until the Braik-Ross gate is resolved or the impulsive layer is
promoted past spike on its own.

---

## Honest limits (carried up front)

- **Sampling, not DA.** No Taylor map, no ADS, no envelope root-find, no
  local-poly speedup, no alpha-shape. The spike is O(n_α·n_β) propagations per
  footprint — fine for a few seeds, not for a swept catalogue. The DA layer is the
  paper's whole performance contribution and is the real adoption cost
  (DACEyPy + a DA-evaluable force model; mining note Sec. 3).
- **Convex-hull over-approximation.** The footprint is the convex hull of the
  sampled cloud, not the true (possibly non-convex) RS envelope. Containment of an
  MC cloud by the hull is a *necessary* coverage check, not a proof the grid
  traces the exact boundary. A real envelope would need the Eq. 34-51 layer.
- **Short arc only.** All cross-checks are `t_f ≤ 0.5 T ≈ 3.3 d`. The same
  long-arc fragility the mining note flags (Sec. 6) applies: our cycler legs are
  100-600 d; the method (and this spike) is credible only on a short final segment
  near a target node — which is where a TCM actually lives.
- **Single impulse, fixed epoch.** Multi-burn / free-epoch TCM windows are a
  genuine gap (mining note Sec. 6), not a wrapper away.
- **No catalogue writeback.** Zhou carries no new sourced cycler tuples and the
  cross-check asserts no golden; nothing was written to `data/catalogue.yaml`.

---

## Files

- `src/cyclerfinder/search/reachable_impulsive.py` — the spike.
- `tests/search/test_reachable_impulsive.py` — 12 fast mechanics tests + 3 `slow`
  NRHO cross-check tests (recovery geometry, MC containment with `P` < 0.1%,
  grid-resolution monotonicity).
- This note.

(Scratch runner used to produce the numbers above lived in gitignored `out/`.)
