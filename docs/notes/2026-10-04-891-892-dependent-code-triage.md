# #891 / #892 dependent-code triage (2026-10-04)

Scope: everything downstream of the two corrected Sun-Earth-Moon models, `core/bcr4bp.py`
(#891, Sun's sense) and `core/qbcp.py` (#892, parities, alpha_6, two coefficient slips),
corrected in commits `8c5b754b` and `11002122`. Neither core module was edited here.

Labels: COMPUTED (a run I made, numbers given), READ (from a file), INFER (my reasoning, not
checked by a run).

## 1. What failed before any repair

COMPUTED. The whole command of the brief does not fit in 8 minutes on this machine (it was
cut off at 470 s twice), because two non-converging torus correctors ran for minutes (see
section 2). I ran the same file list in groups, per file for `tests/search`. Logs in the
session scratchpad, `fix891/g1.log`, `g2.log`, `s_*.log`. Head `11002122`.

| Test | Failure as first seen |
|---|---|
| `tests/data/test_v3_bcr4bp.py::test_v3_bcr4bp_pol1_integrator_independent_passes` | DOP853 against LSODA 351.5 km over 3 cycles, floor 100 km |
| `tests/genome/test_bct_transfer.py::test_backward_arc_reaches_apoapsis` | apoapsis 5.648 LD, expected 3.9 +- 1.2 |
| `tests/genome/test_bct_transfer.py::test_hiten_signature_band` | apoapsis 5.648 LD > 5.1 |
| `tests/genome/test_bcr4bp_torus.py::test_correct_bcr4bp_torus_convergence` | not converged, residual 1.3e-3 after 960 evaluations, 232 s |
| `tests/search/test_cislunar_bct_search.py::test_search_emits_transfer_capability_records` | no transfer-capability record emitted |
| `tests/search/test_sun_forced_periodic_884.py::test_forced_orbit_closure_independent_integrator_and_reversibility` | reverse continuation lands 4.45e-6 from the start point (y, vx only), bound 1e-6 |
| `tests/search/test_variational_periodic_orbit_qbcp.py::test_low_harmonics_converges_residual_but_not_closure` | the low-harmonic solve now passes closure |
| `tests/search/test_variational_periodic_orbit_qbcp.py::test_positive_control_cold_start_reproduces_qbcp_l1_substitute` | 2.31e-2 from the expected point, bound 1e-8 |
| `tests/search/test_variational_periodic_orbit_qbcp.py::test_planar_symmetry_components_are_near_zero` | symmetric component 3.1e-17, test demanded > 1e-4 |
| `tests/search/test_variational_qbcp_arc.py::test_ghost_solution_rejected_by_independent_closure` | loop defect 3.7e-3, test demanded > 0.1 |
| `tests/search/test_variational_qbcp_torus.py` (whole file) | did not finish in 460 s; 12 F among the first 23 results before the cut |

All other files in the list passed (`tests/core`, the other `tests/data` files,
`tests/genome/test_bcr4bp_continuation.py`, `test_bcr4bp_genome.py`, `test_qbcp_torus.py`,
`tests/search/test_variational_ccr4bp_torus.py`, `test_variational_crnbp_torus.py`,
`test_variational_qp_torus.py`, `tests/scripts/test_run_538_residual_shape.py`).

## 2. The two frame-conversion helpers

### 2.1 `genome/bcr4bp_torus.py`

READ. Before this change `se_to_em_transform` put the Sun at `theta_sun0 + omega_S t`
(prograde), used the factor `1 + omega_S` for the Sun-Earth time scale, and anchored the
map on the Earth at `(-mu, 0, 0)` with a time-varying Earth-Sun distance.
`se_lyapunov_to_bcr4bp_torus_seed` had two further defects independent of the sense:
it transformed sample j at Earth-Moon time `j T_em / n` although the GMOS residual and
`evaluate_bcr4bp_torus` treat every point of the invariant circle as living at the section
time t = 0, so the five samples were rotated against each other (by up to about 1.2 rad in
the old time scale and about 30 rad in the corrected one), and it rebuilt a mu = 0 system
without passing `theta_sun0`.

Derivation (the docstring has it in full). With the Moon's mass zero the Earth-Moon origin
is the Earth+Moon mass point, which is the Sun-Earth secondary of mass
`mu_SE = 1 / (mu_sun + 1)`. The Sun is at angle `theta_sun0 - omega_S t` in the Earth-Moon
frame and at angle pi (seen from the secondary) in the Sun-Earth frame, so positions are
related by the rotation `R(alpha)`, `alpha = theta_sun0 - omega_S t - pi`, and the scale
`a_S`. The Sun-Earth frame turns at the Sun's inertial rate `n_S = 1 - omega_S`, so
Sun-Earth time is `tau = n_S t`. Differentiating `r_EM = a_S R(alpha(t)) p(tau(t))`:
`v_EM = a_S R (n_S v_SE - omega_S J p)`. Both factors of the brief become `1 - omega_S`
and the Sun's velocity changes sign, as expected. The secondary now maps to the origin, not
to the Earth: `mu_SE` is defined with Earth + Moon = 1, so the Sun-Earth secondary is the
Earth-Moon barycentre. At mu = 0 both choices coincide, so the identity test cannot tell
them apart; INFER that the barycentre is right from the definition of `mu_SE`. At mu > 0
there is no exact map (the Moon has no Sun-Earth counterpart); the helper is then a seed.

Identity test, COMPUTED (`tests/genome/test_bcr4bp_torus.py::test_se_to_em_transform_mu0_identity`):
a Sun-Earth state off the x axis and out of plane, about 3.4 Earth-Moon distances from the
Earth, `theta_sun0 = 0.7`, t0 = 1.3, transformed, propagated with `bcr4bp_eom` at mu = 0,
compared with the Sun-Earth CR3BP propagation for `n_S dt` transformed at t0 + dt:

| Sun mass | dt | position miss | velocity miss |
|---|---|---|---|
| module value 328900.5423 | 2 | 2.6e-9 | 2.8e-9 |
| module value | 6 | 3.0e-8 | 3.2e-8 |
| `n_S^2 a_S^3 - 1` (Kepler holds) | 2 | 1.1e-12 | 1.1e-12 |
| `n_S^2 a_S^3 - 1` | 6 | 3.7e-12 | 3.5e-12 |
| old transform (prograde, `1 + omega_S`) | 6 | 1.9e2 | 1.9e2 |
| old transform, `n_S` time scale | 6 | 1.6e1 | 1.3e1 |

The module's Sun mass misses the Kepler relation by 2.33e-8 relative (COMPUTED); with a
Sun mass that satisfies it the identity falls to the integrator floor, so the 3e-8 is that
mismatch and nothing else. The test pins both, plus a negative control with the old time
scale.

Seed quality, COMPUTED (SE L2 Lyapunov at Jacobi 3.0008, `P_SE` = 3.1225875, as the
tests build it): GMOS residual of the seed at mu = 0 was 83 with the old helper and is
1.8e-2 with the new one. The remaining 1.8e-2 is Fourier truncation: the seed residual at
mu = 0 falls 1.8e-2, 3.0e-3, 1.6e-3, 2.5e-4, 3.1e-5 for 5, 7, 9, 11, 15 samples. With 5
samples (2 modes) the corrector stops at a floor of 7.2e-4 (`ftol` termination, 5
evaluations); with 11 samples it converges in 4 evaluations (5 s) to 7.8e-7 at mu = 0 and
3.7e-7 at full mu. At mu = 0 the exact rotation number is known independently:
`2 pi n_S T_s / P_SE` = 1.0222010 (rho / 2 pi = 0.162688); the corrected torus has
1.0222014, within 4e-7. At full mu, 1.0221246. The old helper's seed rho / 2 pi was 0.187034,
which is `frac((1 + omega_S) T_s / P_SE)`, the old time scale.

### 2.2 `genome/qbcp_torus.py`

Same defects, READ, plus: its Sun was a circle starting at angle 0, whereas the QBCP's own
Sun, `(alpha_7, alpha_8)`, starts at angle pi in this module's frame (COMPUTED:
`atan2(alpha_8, alpha_7)` = 3.141593 at t = 0, 2.298 at t = 0.9) and regresses. The new
helper takes the Sun's distance `D` and angle from `(alpha_7, alpha_8)(t)` and their time
derivatives (central difference, step 1e-5), maps the secondary to the origin and uses the
mean time scale `n_S`.

What can be tested exactly, COMPUTED
(`tests/genome/test_qbcp_torus.py::test_se_to_em_transform_maps_sun_and_secondary_exactly`):
the Sun-Earth primary lands on `(alpha_7, alpha_8)` to 6e-14 and on its velocity to 1e-8
(speed 370; the reference derivative is taken analytically from the series in the test),
and the secondary on the origin at rest exactly, at t = 0, 0.9, 3.7.

What cannot: the flow. In the coherent model the Sun's motion is not the uniform circle of
the Sun-Earth CR3BP, so there is no exact identity. Measured bound, COMPUTED
(`test_se_to_em_transform_flow_bound_at_zero_moon_mass`, same state and start time as the
bicircular identity): with the Moon's mass set to zero, QBCP flow against transformed
Sun-Earth flow, position / velocity miss 3.5e-8 / 1.6e-7 after 0.5, 3.5e-7 / 6.2e-7 after
2, 7.5e-7 / 9.6e-7 after 6 time units; at the physical Moon mass 1.5e-5, 9.7e-5 and 7.2e-4
in position. The old time scale misses by more than 1. The docstring calls the helper a
seed generator and states this bound; it is not an identity.

Seed, COMPUTED: GMOS residual of the QBCP seed (full mu) 1.9e-2 with 5 samples, 2.0e-3 with
11; the 11-sample seed converges in 4 evaluations (14 s) to 2.8e-7, rho = 1.022113.
