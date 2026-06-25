# #466 Energy-Moving QP-GMOS Continuation — Design Note

**Date:** 2026-06-26
**Lineage:** #290 (single GMOS torus) → #319/#320 (V0–V2_qp + first sweep, 2 SILVER) →
#333 (iso-energetic pseudo-arclength family) → **#466 (energy-moving walk)**
**Status:** Design + TDD build. Report-only — NO catalogue writeback, NO novelty claim.

## The #333 failure mode

`continue_qp_family_arclength` follows the SVD null tangent of the augmented
residual Jacobian. On the #290 smoke seed (k=4, C_J≈3.1279) that tangent points
almost entirely along the rotation number `rho`: the family is a near-iso-energetic
tube and `C_J` spans only ~6e-7 over 41 members. It therefore never approached the
#320 SILVER Bracket-2 region at C_J≈3.032, and the #290↔#320 connection hypothesis
stayed OPEN (see `2026-06-26-333-qp-gmos-family-harvest.md`).

## What the energy ladder actually looks like

The #299 Neimark-Sacker bracket inventory (`data/family_296_3d_subfamilies_299.jsonl`)
is a clean ladder in `C_J` along the #296 (1,1) 3D vertical parent family:

```
k=6  C≈3.024   k=4  C≈3.032 (SILVER Bracket-2)   k=5  C≈3.057 / 3.084
k=3  C≈3.108   k=4  C≈3.126 / 3.128 (the #333 smoke seed)
```

So C_J≈3.128 (seed) and C_J≈3.032 (SILVER) are TWO brackets on the SAME parent
family, ~0.096 apart in energy. The question "do they connect" is: can a QP-torus
family be continued in energy from one to the other, or does a fold / Arnold-tongue
phase-lock / mode-truncation breach block the descent?

## Design decision: parent-family-driven energy continuation (rho free)

Two options were considered:

1. **Energy-PIN on the torus mean state `c_0`** — add a row
   `C_J(c_0) - C_target = 0` to the GMOS+phase residual and march `C_target`. This was
   built and TIMED (`scripts` diagnostics): a single corrector step is ~75-230s and
   does NOT converge even for a 1e-4 energy step. The energy pin fights the GMOS
   invariance rows for the same 6 DOF of `c_0` (the torus mean state IS the parent
   periodic-orbit state, whose energy is set by the parent family, not by a free
   gradient nudge). Rejected: architecturally wrong + compute-infeasible.

2. **Parent-family-driven continuation (CHOSEN).** The #296 (1,1) 3D vertical parent
   family (`data/family_296_3d_em_11.jsonl`, 265 members) IS already an energy
   continuation: `C_J` is strictly monotone in `step_index`, spanning 2.920…3.148.
   The #333 seed sits at **step 112, C_J≈3.12785**; the #320 SILVER Bracket-2 sits at
   **step 8, C_J≈3.03196** — 105 parent members apart on ONE monotone ladder. So the
   energy-moving QP walk is: step the *parent orbit* down the family member-by-member,
   and at each parent member re-converge the QP-torus with `correct_qp_torus`,
   **warm-started from the previous member's converged Fourier modes / rho / t_strob**.
   The parent orbit's energy IS the continuation parameter; the torus rides it down.
   This sidesteps the energy-pin-vs-GMOS tension entirely (each torus is a clean GMOS
   solve at a fixed parent energy) and reuses the proven `correct_qp_torus` corrector
   with no new ill-conditioned constraint row.

Choice (2) is chosen because the parent family is a *ready-made, monotone, dense*
energy ladder connecting the seed to the SILVER region, and because re-converging the
torus at each parent energy is a clean, fast GMOS solve (the seed build itself is
~40s; warm-started members are faster). It MOVES in energy by construction (the parent
member's `C_J` is the member energy) and the connection question becomes concrete:
does a valid irrational warm-started torus survive all the way from step 112 to step 8?

### Per-member procedure

For each parent member `p` (stepping `step_index` from 112 toward 8):

```
seed_modes = previous member's (coeffs, rho, t_strob)       (warm start; seed uses GMOS bootstrap)
torus_p = correct_qp_torus(parent_state[p], parent_period[p], floquet_pair[p], warm-start=seed_modes)
classify: irrational? V1_qp? tail<guard?  -> EnergyWalkMember at C_J = parent jacobi[p]
```

The parent member already carries its Neimark-Sacker eigen-pair only at the bracket
steps; away from a bracket we re-derive the transverse-mode seed by continuation of the
previous torus's `c_1` (the warm start), correcting onto the new parent orbit.

### Termination / blocker reasons

- `reached_target` — reached the requested `target_step` (e.g. step 8, the SILVER).
- `corrector_fail` — the GMOS corrector cannot produce a valid torus at parent member
  `p`: a **family boundary** (the characterized negative — report the `C_J` where it
  stops, i.e. how far down the ladder the QP structure survives).
- `resonance_lock` — `freq_ratio` crosses a low-order p/q (Arnold tongue): the torus
  degenerates to a phase-locked periodic orbit; report the locking C_J and p:q.
- `mode_truncation_breach` — the N=2 tail `|c_N|/|c_1|` exceeds the guard: the torus
  thickens beyond the truncation's validity; report the C_J.

## Golden discipline

Same as #333: every EXPECTED side asserts topology (irrationality), invariance
(Fourier / off-grid closure), or self-consistency (energy genuinely MOVES this time —
the inverse of the #333 containment assertion: `C_J` must span a real range ≫ the
#333 6e-7 floor). NO frequency or C_J value our own code produced is asserted as a
target. The march TARGETS (C* increments) are inputs, not goldens; the member
energies tracking those targets is a self-consistency check.

**Olikara golden:** NOT in the corpus (`grep -ci olikara docs/notes/CORPUS_INDEX.md`
→ 0; no Olikara/Henderson-Howell QP-tori family table PDF present). Registered as a
standing acquisition gap (per #333 follow-on + `feedback_never_give_up_reproducing_papers`).
The energy walk is NOT blocked on it; it ships with self-consistency goldens.

## No catalogue admission

Any candidate quasi-cycler surfacing from the walk is flagged for human gauntlet
review after `search/literature_check.py`, NEVER self-admitted (it is a QP-torus,
which is never "novel" until it clears the published record).
