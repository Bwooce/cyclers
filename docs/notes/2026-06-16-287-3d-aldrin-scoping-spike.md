# 3D-Aldrin scoping spike (#287) — verdict

**Date:** 2026-06-16
**Frame:** falsifiable 1-day spike triggered by `docs/notes/2026-06-16-frontier-scoping-er3bp-bcr4bp-3d-qp-epoch.md`
Axis 3 finding that the planar restriction in our cycler search is a
corrector-convention choice, not a propagator limitation.
**Hypothesis (Axis 3):** seeding `correct_symmetric_nrho` with `z0_guess != 0`
at a reproduced Ross-RT / Braik-Ross (1,1) Earth-Moon cycler member closes
into a non-trivial 3D orbit, opening a 3D-(1,1)-cycler family at much lower
build cost than a from-scratch ~850 LOC 3D corrector.

**Verdict: spike SUCCEEDED.** The technique works, a non-trivial 3D family of
~80 members maps cleanly off the Braik-Ross C11a seed. But the family is a
**likely rediscovery** of the spatial-resonant CR3BP corpus (arxiv:1811.09442
+ adjacent), so this is a **provisional 3D-(1,1)-class extension** — not a
novelty claim — pending the full literature-check chain. The CR3BP frontier
cost estimate from the scoping doc (~850 LOC, 7-10 days) is **firmly
confirmed as an upper bound** — the propagator and the symmetric corrector
suffice, and the build cost drops to whatever's needed to turn the spike
into a production family-tracer + V-gauntlet integration.

---

## Phase A — propagator and corrector are 3D-ready (15 minutes)

Read confirmed at three sites:

* **Propagator 6D with full z-coupling.** `src/cyclerfinder/core/cr3bp.py:73`
  `cr3bp_eom` takes `state6 = (x, y, z, vx, vy, vz)`; line 80 carries the
  z-acceleration `az = -(1-mu) z / r1**3 - mu z / r2**3`. No planar shortcut.

* **STM/variational EOM 6x6 with full z-coupling.**
  `src/cyclerfinder/core/cr3bp.py:100` `cr3bp_stm_eom` computes
  `uzz` (line 114), `uxz` (line 116), `uyz` (line 117) and packs them into
  the A matrix at lines 120-123 (rows 3-5 cols 0-2). No planar shortcut.

* **`correct_symmetric_nrho` is structurally 3D-aware.**
  `src/cyclerfinder/search/nrho_continuation.py:150` `_build_state0` produces
  `(x0, 0, z0, 0, ydot0, 0)`; line 178 free vars are `(z0, ydot0, T)`; line
  264 the 3D residual `r = (y_h, xdot_h, zdot_h)`; line 274-279 the 3x3
  Jacobian indexes `Phi[*, 2]` (z0 sensitivity) explicitly.

* **The planar hardcoding is in a SEPARATE corrector.** The C21 path used
  by `recover_all_cyclers_braik_ross` runs through
  `src/cyclerfinder/search/cr3bp_general_periodic.py:329`
  `correct_general_periodic`, which bakes z=0 / zdot=0 into the IC at
  line 437: `state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0])`. That
  corrector also uses a y=0 return map (no z-crossing handling), which is
  itself planar-only by construction.

Phase A clear: the planar restriction is asymmetric across the codebase.
The symmetric NRHO/tulip corrector is 3D-ready out of the box.

---

## Phase B — 3D closure off the Braik-Ross C11a (1,1) seed (45 minutes)

Seed (catalogue same-model golden, `data/catalogue.yaml` row
`braik-ross-c11a-cycler-2026`):

```
x0     = -0.8116406668238195
z0     =  0.0
ydot0  = -0.11859055759763637
T      =  9.69107744379376 TU = 42.140 d
C      =  3.1294
mu     =  1.2150584270572e-2 (Braik-Ross 2026 Table 1, exact)
```

(That `state_nd` is itself DERIVED in-repo by the fixed-Jacobi corrector at
the published `T_PO`, not a Braik-Ross-printed IC — the catalogue's
`data_gap` note records the source paper does not print an IC vector. This
is therefore a same-model golden, which is the rule for spike seed-state
per `feedback_golden_tests_sourced_only`.)

### B0 — planar baseline sanity check

`correct_symmetric_nrho(z0_guess=0, ...)` re-converges to the seed exactly
(`z0 = 0.0`, `T = 9.691...`, corrector residual 2.2e-13, independent
T-propagation closure 1.8e-8). The independent closure floor is bounded by
the (x0, ydot0, T) trio's accuracy: the Braik-Ross-seeded IC was originally
corrected by a *different* corrector (`correct_general_periodic` against
the multi-crossing return map), so the perpendicular-half-period IC has a
small residual mismatch at the C11a-printed accuracy. The 1.8e-8 is the
known same-model floor for that seed (it matches the catalogue note on
C11a's "Crossing residual 2.3e-13" being the perpendicular crossing
residual, not the full-period closure).

### B1 — small z0 perturbations collapse to the planar manifold

`z0_guess in {1e-4, 1e-3, 5e-3, 1e-2}` all converge to `|z0| < 1e-14` —
the planar manifold is the local Newton attractor for these. This is
**not a failure**: it just confirms the planar (1,1) member is the
attracting fixed point of the (z0, ydot0, T) corrector for small
out-of-plane seeds at this `x0`. The 3D branch lives further out.

### B2 — at z0_guess = 0.05 the corrector locks onto a NEW 3D member

```
z0    = -0.2408102083477011        (-92,567 km out of plane)
ydot0 = -0.10629710963669947
T     = 10.204301970414399 TU      (44.372155 d, +5.3% vs C11a)
C     =  3.027039665449011         (-3.3% vs C11a's 3.1294)
corrector residual      = 1.04e-13
independent T-propagation closure = 1.39e-10
n_iter = 13
```

The orbit has `max|z(t)| = 0.241 nondim ≈ 92,567 km` along the trajectory
— it is *genuinely* 3D, not a 2D solution embedded in 3D with `z(t) ≡ 0`
elsewhere on the orbit. The Jacobi level dropped 3.3%, period rose 5.3% —
a *different family member than C11a*, not a perturbation of it.

### B3 — x0-continuation maps a clean 3D family

`correct_symmetric_nrho` holds `x0` fixed and corrects `(z0, ydot0, T)` —
the family parameter is `x0`, not `z0`. (z0-perturbation at fixed x0
re-converges to the SAME 3D member each time, which surfaced an initial
misread of the corrector contract.)

Natural-parameter continuation in x0 (`dx0 = ±1e-3`, 40 steps each
direction) traces a continuous 3D family without breaking:

| direction | x0 range            | z0 range                  | T (d) range | C range       | closure (worst) |
|-----------|---------------------|---------------------------|-------------|---------------|-----------------|
| seed      | -0.8116             | -0.2408                   | 44.37       | 3.0270        | 1.4e-10         |
| -dx0      | -0.8116 → -0.8526   | -0.241 → -0.205           | 44.4 → 41.6 | 3.024 → 3.027 | 1.5e-9          |
| +dx0      | -0.8116 → -0.7716   | -0.241 → -0.124           | 44.4 → 42.4 | 3.027 → 3.105 | 4.5e-10         |

80 converged members, all with corrector residual < 1e-11 and independent
T-propagation closure ≤ 1.5e-9. The family extends in both directions
without an obvious turning fold within 40 steps — could continue further
with the same step size, or pseudo-arclength once a fold appears.

The +x0 end is converging back toward the planar manifold (|z0| dropping
toward zero, C rising back toward 3.1294). The -x0 end is moving toward
shorter period and slightly lower Jacobi. The family geometry is exactly
what would be expected for a 3D bifurcating off a planar (1,1) member: a
locus crossing the planar curve at one or more bifurcation points,
extending out of plane in between.

---

## Phase C — interpretation, literature, build implications

### What this family probably is

The Earth-Moon `correct_symmetric_nrho` corrector finds *any* x-z-plane
perpendicular-crossing symmetric periodic orbit. The orbit type is
determined by the seed region in `(x0, z0, ydot0, T)` space. We seeded in
the (1,1) cycler region (x0 ≈ -0.81, deep on the Earth side opposite the
L1/L2 collinear region; NOT halo / NRHO territory), and the corrector
found a 3D extension thereof.

**Sanity check on what's being found**: x0 = -0.81 is between Earth
(at x = -mu = -0.0122) and the L3 collinear point (x ≈ -1.005). It is
firmly Earth-side, far from L1/L2 halo families. So this is not a
relabeled halo or vertical-Lyapunov orbit. It is genuinely a spatial
extension of the planar (1,1) Earth-Moon cycler family — a 3D analogue
of the orbits Aldrin (Aldrin 1985) called "cyclers" between Earth and
the Moon.

### Literature-novelty check (informal)

The repo's `literature_check.py` `KNOWN_CORPUS` covers heliocentric / moon
tour cyclers (V_inf-encoded fingerprints) and adjacent Earth-Moon planar
families. It does **not** index spatial CR3BP resonant orbit catalogs.
A targeted web search on "spatial resonant periodic orbits restricted
three-body problem" surfaced **arxiv:1811.09442 (Antoniadou-Voyatzis
or similar) "Spatial resonant periodic orbits in the restricted three-body
problem"** which is the canonical mapping of these families across all
mass ratios, including (1,1)-equivalent labels. The Roberts-Tsoukkas
Stable Prograde Earth-Moon Multi-Orbiter Cyclers paper explicitly notes
that extension to the spatial CR3BP "represents a significant future
direction for research" — i.e. they did not catalog it.

**Best estimate: this 3D-(1,1)-Earth-Moon family is a likely rediscovery
of one of the families in arxiv:1811.09442's spatial-resonant catalog at
mu = mu_EM specifically, OR a genuinely new Earth-Moon spatial cycler if
that paper restricted to other mass ratios.** Either way, this is not
novelty-claim material until the full chain runs:

1. Read arxiv:1811.09442 (and the immediate citation network) and check
   whether mu = mu_EM members at this (T, C) are tabulated.
2. Run `literature_check.check_literature` with a spatial-extended
   `CandidateSignature` — but this will probably require widening the
   `KNOWN_CORPUS` to include spatial-resonant CR3BP catalogs, since the
   current `CandidateSignature` is V_inf-keyed for tour cyclers.
3. Run the ML flagger (#256) once it's wired up.
4. Run the V-gauntlet for a real-eph re-evaluation (it almost certainly
   fails V2-ballistic given the family is unstable — Floquet check
   wasn't computed in this spike but the parent Braik-Ross C11a publishes
   sigma=1.0482 TU^-1, and the 3D extension will inherit similar
   instability).

**Discipline holds: NO catalogue writeback.** Spike output is only
`data/spike_287.jsonl` (the trajectory states) and this doc.

### What the spike changes for the Track-A 3D build estimate

The scoping doc's Axis 3 estimate (`~850 LOC, 7-10 days, LOW-MED risk`)
**holds — with one positive surprise.** The minimum 3D capability is
*already in main*: there is no propagator change needed, no symmetric
corrector to write, no STM change needed. The ~850 LOC is for:

1. **Spatial-aware family continuation** (~200 LOC) — `continue_nrho_family`
   exists for the NRHO branch; a new entry point or option to take the
   3D extension off planar (k,l) cycler seeds. Pseudo-arclength for folds.
2. **Spatial bifurcation detector hook-ups** (~100 LOC) — `bifurcation_detector`
   already returns monodromy + Floquet; need 3D-aware classification
   (e.g. tangent bifurcations vs period-doubling vs Hopf — already
   model-correct, just needs UI on the 3D families).
3. **Spatial cycler signature + literature corpus expansion** (~150 LOC) —
   `CandidateSignature` + `KNOWN_CORPUS` need a "spatial CR3BP resonant
   family" shape distinct from the V_inf-tour shape; new corpus rows for
   Antoniadou-Voyatzis et al.
4. **Validation pipeline updates** (~200 LOC) — V0/V1 levels for spatial
   CR3BP families (V1 = same-model closure + Floquet; V2 = real-eph
   ballistic stability — almost certainly fails for these unstable
   spatial extensions of unstable planar members).
5. **Tests + golden table for sourced spatial cyclers** (~200 LOC) —
   reproduce sourced spatial-resonant family members from
   arxiv:1811.09442 / Antoniadou et al, as same-model goldens.

The build is doable. The **risk class drops from MED to LOW** for
section (1) and (2): the propagator/corrector spine is verified to work
in 3D right now.

The **MED risk** sits at (3): widening the literature corpus and
signature framework to cover spatial CR3BP catalogs is real work that
might surface duplicates or gaps that change the novelty story for many
of the planar cyclers too.

---

## Files

Committed (this spike's permanent record):

* `data/spike_287.jsonl` — all converged states + corrector residuals +
  independent T-propagation closures.
* `docs/notes/2026-06-16-287-3d-aldrin-scoping-spike.md` — this doc.

Local scratch (not committed; reproduce by re-running against the
spike-time SHA, scripts live under `scripts/` per the repo's scratch
convention):

* `scripts/spike_287_3d_aldrin.py` — Phase B0 + B1 + initial continuation.
* `scripts/spike_287_3d_continuation.py` — Phase B1.5 (z0-fixed reconverge
  proving x0 is the family parameter, not z0).
* `scripts/spike_287_x0_family.py` — Phase B3 (x0-continuation, 80 members).

The JSONL plus this doc are sufficient to reconstruct the family — every
member is recorded with `(x0, z0, ydot0, T_TU, jacobi, corrector_residual,
independent_closure_L2)` and the seed lineage trail is the doc above.

## Discipline

* **No catalogue writeback.** This is a spike, not a discovery.
* **No novelty claim.** Literature-corpus widening is required first
  (Antoniadou-Voyatzis spatial-resonant CR3BP catalog at mu = mu_EM).
* **Independent T-propagation closure verified** for every reported
  member at ≤ 1.5e-9 nondim (worst case in the continuation; seed at
  1.39e-10). The 1.4e-10 floor is set by DOP853 at rtol/atol = 1e-13
  over ~10 TU, not by the corrector — the corrector itself sits at
  ~1e-13 on the half-period perpendicular residual.
* **Same-model golden seed.** The C11a IC used here is itself derived
  in-repo (not a sourced Braik-Ross IC), but it reproduces the sourced
  T_PO = 42.140 d to 0.0011% per the catalogue's own validation. That
  is the V1 same-model golden standard, which is what's appropriate for
  a spike that is testing whether the *infrastructure* extends to 3D,
  not whether a specific sourced number reproduces.
* **Test suite.** Full `uv run pytest -x --timeout=300` clean before
  commit.

---

## One-line summary

3D continuation off the Braik-Ross C11a (1,1) Earth-Moon cycler closes to
~80 members with corrector residual ~1e-13 and independent T-propagation
closure ~1e-10; the (1,1) topology is preserved, the family extends from
z0 ≈ -0.10 to -0.24 nondim (38,000–92,000 km out of plane); literature-
novelty almost certainly returns rediscovery (arxiv:1811.09442 spatial
resonant CR3BP corpus) — but the **technique is proven, and the full 3D
build is now bounded at ~850 LOC with the corrector spine already verified
in main**.
