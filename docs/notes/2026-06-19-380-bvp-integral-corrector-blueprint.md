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

## DEFERRED (Steps 4-6) — gated on these open questions
1. **#378 Phase 3.2:** does it strictly need integral rows IN the Newton system,
   or does a post-hoc `sun_phase_drift` check on `correct_bcr4bp_periodic`
   already suffice? (Architect: the latter may suffice — #380 then becomes a
   forward-looking capability, not a hard gate.)
2. **#388 SnLm multi-arc routes through the HELIOCENTRIC DSM corrector**
   (`search/dsm_leg.dsm_chain_correct`, Lambert-based residual) — NOT the
   rotating-frame CR3BP/BCR4BP corrector this module augments. So the
   arc-duration-sum constraint for #388 may need a separate heliocentric variant.
3. The pure Jacobi-drift constraint has ZERO gradient in CR3BP (conserved) -> it
   cannot steer the corrector; the actionable constraint is the time-integral /
   arc-duration-sum (non-trivial gradient).
4. `IntegralConstraint.weight`: user-supplied scalar (default 1.0), applied to
   the row before stacking, to control `lstsq` conditioning vs the point rows.

Revisit Steps 4-6 once #378/#388 routing is confirmed.
