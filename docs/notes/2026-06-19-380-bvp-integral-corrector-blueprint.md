# #380 BVP Integral-Constraint Corrector — design blueprint

**Date:** 2026-06-19. **Status:** core module (Steps 1-3) approved for build;
wiring into existing correctors (Steps 4-6) DEFERRED pending the scope decisions
below. Source: code-architect analysis of the existing corrector/STM stack +
the Shepperd-1985 / Bond-Allman-2021 / Gurfil-2007 methodology digests.

## Approach (chosen): augmented-STM quadrature
Append `n_q` integral-accumulator rows and their variational rows to the
existing 42-dim `[X(6), Phi(36)]` ODE, so each integral constraint value
`q_i(T)` and its sensitivity `dq_i/dX0` come out of ONE `solve_ivp` call:
- `dq_i/dt = h_i(X(t), t)`,  `q_i(0)=0`  -> integral value at T.
- `dpsi_i/dt = (dh_i/dX(X,t)) @ Phi(t,t0)` (1x6) -> sensitivity row.
Augmented state dim = `42 + 7*n_q` (49-56 for 1-2 constraints). Residual stacks
point rows (existing) + integral rows `q_i(T) - target`; Jacobian gains rows
`[dq_i/dX0 (free cols) | h_i(X_f,T)*t_scale (T col)]`. Newton + backtracking on
the total residual norm; independent Radau closure cross-check on point rows.

Rejected: collocation (full re-architecture, marginal benefit on smooth
orbits); post-hoc finite-difference quadrature (6x propagation cost/Newton
step). Augmented-STM is the textbook method (Bond-Allman 2021 Ch 14;
Shepperd 1985 quadrature appended to universal-variable propagation).

## Build scope APPROVED NOW (Steps 1-3, additive — new files only)
- CREATE `src/cyclerfinder/search/bvp_integral.py`:
  `IntegralConstraint` (label, integrand `h(t,X,params)`, grad `dh/dX` or None,
  target, weight), `AugmentedArc` (state_f, stm, q_values, dq_dX0),
  `IntegralCorrectorResult`, `propagate_augmented_cr3bp`,
  `propagate_augmented_bcr4bp`, `correct_with_integral_constraints`, and
  factories `cr3bp_jacobi_drift_constraint`, `sun_commensurate_period_constraint`.
- CREATE `tests/search/test_bvp_integral.py` (sourced golden tests):
  1. Jacobi-drift integral = 0 on a converged catalogue periodic orbit
     (analytical conservation; EXPECTED side = 0.0, not our-computed).
  2. BCR4BP Jacobi-drift non-zero but bounded (POL1 / Rosales-Jorba 2023 IC).
  3. Corrector converges with a time-integral (period) constraint on a SOURCED
     L1 Lyapunov orbit (period from catalogue golden, not our-computed).
  4. Augmented Jacobian matches central finite-differences (<1e-6) — consistency.
  5. (Step 6, deferred) arc-duration-sum = Sun-commensurate period T_1.

Seams (for the deferred wiring): `core/cr3bp.py:108-162` `cr3bp_stm_eom`;
`core/bcr4bp.py:286-318` `bcr4bp_stm_eom`; `genome/bcr4bp_genome.py:157-168`
`_propagate_with_stm` + Newton loop `316-392`; `genome/multi_shooting.py:185-267`
`_residual_and_jacobian`; `search/cr3bp_general_periodic_3d.py:124-168`.

## Steps 4-6 — RESOLVED: NOT NEEDED (2026-06-19, routing confirmed in code)

The two consumers #380 was speculatively framed to "gate" do NOT route through
this rotating-frame integral-constraint machinery, so the wiring is unnecessary.
#380 stands **CORE-COMPLETE** (Steps 1-3) as a forward-looking capability.

1. **#378 BCR4BP Phase 3.2 — already handled, no integral rows needed.**
   `genome/bcr4bp_genome.correct_bcr4bp_periodic` already enforces Sun-
   commensurate periodicity by taking `sun_commensurate_n`, correcting toward
   `sun_commensurate_period(omega_sun, n)`, and computing
   `sun_phase_drift = |omega_sun*T - 2*pi*n|` POST-HOC at the converged period
   (bcr4bp_genome.py:425, reported in the result :466). The architect's Risk #5
   is confirmed: the post-hoc drift check suffices; integral rows in the Newton
   system add nothing for Phase 3.2.
2. **#388 SnLm multi-arc — heliocentric, different machinery.**
   `search/dsm_leg.py` is a propagate-then-Lambert HELIOCENTRIC corrector
   (residual = match-point defect in km/s; depends only on core/kepler +
   core/lambert; already uses the Shepperd-STM analytic Jacobian via
   `core/fbs_match_point.match_point_defect_jacobian`). It is NOT a rotating-
   frame CR3BP/BCR4BP corrector, so the augmented-quadrature integral-constraint
   module does not apply. Architect's Risk #4 confirmed.

So #380 = the augmented-quadrature propagator + corrector + time-integral
constraint (Steps 1-3, committed 0b0bdea), available for any FUTURE rotating-
frame integral-constraint need. No wiring lands; the "forced gate" framing was
speculative. Carried-forward design notes (if a future need arises):
- The pure Jacobi-drift constraint has ZERO gradient in CR3BP (conserved) and is
  numerically pathological as an augmented quadrature (removed; #380 commit) —
  monitor Jacobi drift post-hoc, not as a constraint row.
- `IntegralConstraint.weight`: user-supplied scalar (default 1.0) applied to the
  row before stacking, to control `lstsq` conditioning vs the point rows.
