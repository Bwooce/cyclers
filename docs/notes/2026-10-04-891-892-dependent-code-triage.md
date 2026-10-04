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

### 2.3 A third copy of the QBCP field

COMPUTED. `search/variational_qbcp_torus.py` carries its own vectorised copy of the QBCP
right-hand side and Jacobian (`_qbcp_rhs_grid`, `_qbcp_jacobian_grid`) for the
pseudospectral torus corrector. It still had #592's Sun-only `1/alpha_6` factor, so after
the core correction its own gate tests (`test_qbcp_rhs_grid_matches_pointwise_eom`,
`test_qbcp_jacobian_grid_matches_pointwise_stm_a_matrix`) failed (velocity rows differed by
5e-2). Corrected to alpha_6 on the whole potential; both pass. Everything computed with the
pseudospectral QBCP torus corrector before today (#617, #618, #619, #620, #626, #646, run
through `scripts/run_538_qbcp_cycler.py`, which imports this module, READ) used this copy
as well as the defective core module (INFER from the import; the runs were not re-traced).
No other copy of either field was found in `src/` or `scripts/` except
`scripts/analyze_593_qbcp_l1_substitute_reconciliation.py` (grep, READ; not rerun).

## 3. Triage of every failing test

Classes: (a) structure or identity hard-coding the old convention; (b) regression pin of our
own number; (c) scientific result or positive control that was a self-comparison. All
numbers COMPUTED on 2026-10-04 unless marked.

| Test | Class | What was done | Before | After |
|---|---|---|---|---|
| `genome/test_bcr4bp_torus.py::test_correct_bcr4bp_torus_convergence` | (c), holds | Seed helper fixed (section 2); 11 samples instead of 5; new check of rho against the exact mu = 0 value | not converged, 1.3e-3, 232 s | 7.8e-7, rho within 4e-7 of exact, 5 s |
| `search/test_variational_qbcp_torus.py::test_qbcp_rhs_grid_matches_pointwise_eom`, `..._jacobian_grid_matches_pointwise_stm_a_matrix` | (a) in the code, not the test | grid copy corrected (2.3); tests unchanged | rows off by 5e-2 | pass |
| `search/test_variational_qbcp_torus.py::test_se_l2_positive_control_reproduces_gmos_torus` | (c), holds | builder: 11 samples, BCR4BP bootstrap at the QBCP's Sun phase pi; pins recomputed; added the Sun-Earth three-body rotation number 0.162688 as an independent anchor (rel 2e-3) | rotation number pinned 0.231365, rms 9.405e-5 | 0.162806 (GMOS 0.162809), rms 7.350e-5, closure 3.7e-5 |
| same file, the eleven `se_tori` fixture tests | (b) | rebuilt on the new fixture; all pass except the next row | fixture 240 s, non-converging | fixture 40 to 85 s per worker |
| `search/test_variational_qbcp_torus.py::test_unstable_forward_vs_backward_differ_on_torus` | (b) | bound 0.9 was calibrated on the old torus; measured the angle at 16 phases (2 to 25 deg); now `< 0.99` plus a pin of 0.9567 | 0.9567 > 0.9 | pass |
| `search/test_variational_qbcp_torus.py::test_em_l2_c313_crosses_gmos_plateau`, `test_em_l2_exact_and_lsmr_agree` (both `slow`, not in the gate, did not fail in Task 1) | (c), cannot tell | strict xfail, reason in the file | pins 3.412e-3, 1.9183, plateau 4.771e-1, all old-model | not run (13 to 22 min each) |
| `search/test_variational_qbcp_arc.py::test_ghost_solution_rejected_by_independent_closure` | (b) with a check | ghost defect is now 3.7e-3, not O(1); refinement from the same seed (orders 20, 28, 40, 52: 1.1e-3, 3.7e-3, 1.2e-2, 2.0e-2, residual 1e-15) shows it is still a ghost; test now requires defect > 2e-3, veto, and growth with order | 3.7e-3 < 0.1 | pass |
| `search/test_variational_qbcp_arc.py` SE-L2 builder | (b) | same builder change as above | | pass |
| `search/test_variational_periodic_orbit_qbcp.py::test_positive_control_cold_start_reproduces_qbcp_l1_substitute` | (c), upgraded | the expected state was a multiple-shooting state computed in the defective model; now compared with the PUBLISHED POL1 point | 2.31e-2 from the self-computed state | 1.7e-8 from published POL1 (bound 1e-7) |
| `search/test_variational_periodic_orbit_qbcp.py::test_planar_symmetry_components_are_near_zero` | (a) | asserted y, px > 1e-4, the defective model's asymmetry; now asserts they vanish (< 1e-10) | 3.1e-17 failed `> 1e-4` | pass |
| `search/test_variational_periodic_orbit_qbcp.py::test_low_harmonics_converges_residual_but_not_closure` | (b) with a check | at 8 harmonics the corrected model now closes (5.9e-4); moved to 4 harmonics where residual 1.428e-7 passes `tol` but closure is 0.26, restoring the gate demonstration | old pin: plateau 9.396e-5 at 8 | residual 1.428e-7, closure 0.26 |
| `search/test_sun_forced_periodic_884.py::test_forced_orbit_closure_independent_integrator_and_reversibility` | (a) | the coordinator's reading confirmed: landing point is a phase slip of -1.20e-5 TU at both nodes, 6.9e-13 and 1.8e-12 from the orbit; test now minimises the distance to the orbit over a common phase slip (< 1e-9) | 4.45e-6 from the start point | 1.8e-12 from the orbit |
| `genome/test_bct_transfer.py::test_backward_arc_reaches_apoapsis`, `test_hiten_signature_band`; `search/test_cislunar_bct_search.py::test_search_emits_transfer_capability_records` | (c), holds at another Sun phase | see 3.1 | 5.648 LD at the 70-day window edge | 3.178 LD at 45.9 d, interior |
| `data/test_v3_bcr4bp.py::test_v3_bcr4bp_pol1_integrator_independent_passes` | (c), does not hold | strict xfail; see 3.2 | | |

Two new tests per helper are listed in section 2. Commits: `49df785e`, `7c5e3e55`, `4f65771f`,
`0696575d`, `e175c833`, `9aa470a9`, plus the commit of this note's final version.

### 3.1 The Hiten band

READ: the band (apoapsis 3.9 LD +- 30 %, so 2.7 to 5.1 LD) traces to Belbruno 2004 section
3.4; the Delta-V split in `build_hiten_bct` is the published value carried as a constant, so
`test_hiten_signature_band`'s Delta-V assertion is not a computed check (pre-existing; INFER
that only the apoapsis and the capture energy are computed). COMPUTED: at Sun phase 0 the
backward arc from theta2 = 0.70 climbs to 5.65 LD and is still climbing when the 70-day
window ends, so the reported "apoapsis" is the window edge. In the old-sense model (the
corrected module with omega_S negated, which reproduces `theta0 + omega t`) the same arc
had an interior apoapsis of 3.54 LD at 58.2 d. A scan in the corrected model: over theta2 at
Sun phase 0, one of 24 values (0.785) gives an interior apoapsis in band; over Sun phase at
theta2 = 0.70, phases 0.30 to 0.75 rad give interior apoapses of 3.74 to 2.81 LD at 57 to
40 d, smoothly, and the same at phase + pi (0.524 and 3.665 both give 3.131 LD, the
expected tidal symmetry). The tests now use Sun phase 0.5 rad (3.178 LD at 45.9 d; 3.743 LD
at 54.2 d for theta2 = 0.75) and require the apoapsis to be interior. The band was not
changed. This is a re-selection of a free parameter, disclosed in the test; the claim that
survives is "a backward arc from a W point reaches a genuine apoapsis in the Hiten band for
a range of Sun phases", which is what the test checks.

### 3.2 The V3 "POL1" test (strict xfail)

COMPUTED: the closed orbit's DOP853-LSODA agreement is 0.0034, 1.09 and 351.5 km at cycles 1
to 3; the ratios are 321.7 and 321.9, and the orbit's dominant monodromy eigenvalue is 321.8.
So the disagreement is integrator noise at the 1e-9 level amplified by the orbit's
instability, and the 100 km floor over 3 cycles would need cycle-1 agreement of about 1 m.
I did not loosen the floor or the cycle count. The dominant multiplier in the old-sense
model is 406.6 (COMPUTED), so the old pass was not a matter of a less unstable orbit; INFER
that the old cycle-1 agreement happened to be smaller. Separately, READ: the seed
`_POL1_X = -0.8369...` with `py = -0.8391...` is POL1 in the paper's frame (Earth at +mu).
The bicircular module puts the Earth at -mu, so this seed is 0.82 Earth-Moon distances from
the Earth on the side away from the Moon, not near L1, and the closed orbit (x = -0.856) is
not the L1 dynamical substitute. The same seed is used in `tests/core/test_bcr4bp.py`,
`tests/data/test_v0_bcr4bp.py`, `test_v1_bcr4bp.py`, `test_v2_bcr4bp.py` and
`tests/genome/test_bcr4bp_genome.py`; those pass, but their docstrings call it the L1
substitute. Not changed here (out of the #891 scope); reported.

## 4. Tests left as strict xfail (open scientific questions)

| Test | Open question |
|---|---|
| `tests/data/test_v3_bcr4bp.py::test_v3_bcr4bp_pol1_integrator_independent_passes` | Is a 3-cycle V3 integrator-agreement claim meaningful for an orbit with multiplier 322 per cycle, and what orbit was meant (the seed is not near L1 in this frame)? |
| `tests/search/test_variational_qbcp_torus.py::test_em_l2_c313_crosses_gmos_plateau` (slow) | Does the pseudospectral corrector still cross the EM-L2 GMOS plateau in the corrected QBCP, and at what residual and rotation number? |
| `tests/search/test_variational_qbcp_torus.py::test_em_l2_exact_and_lsmr_agree` (slow) | Same torus: do the two trust-region solvers still agree on one minimum? |

## 5. Exposure report (read-only; no data changed)

READ by a read-only survey (file:line evidence kept in the session scratchpad,
`fix891/task4_exposure.md`); I spot-checked the four rows at lines 25, 39, 60 and 81 and the
line list of the grep. "Integrated" means the method propagated `core/bcr4bp.py` or
`core/qbcp.py`, directly or through a module built on them.

### 5.1 `data/empty_regions.jsonl` entries mentioning the bicircular model (19 lines)

| Line | region_id | Producer | Verdict |
|---|---|---|---|
| 25 | er3bp-discovery-em-broucke-koblick-2026-06-24 | `scripts/run_432_er3bp_discovery.py` | prose (cites a scoping note whose filename contains "bcr4bp"; method is ER3BP) |
| 36 | isolated-er3bp-cycler-novelty-gap-analysis-2026-06-25 | none (#442 adjudication) | prose (same citation) |
| 39 | cislunar-bct-wsb-quasicycler-2026-06-26 | `search/cislunar_bct_search.py` via `genome/bct_transfer.py` and `core/wsb.py` (#378); entry written by `scripts/_apply_378_empty_region.py`; apoapsis spike `scripts/spike_378_bct_apoapsis.py` | INTEGRATED: method-invalid |
| 51 | floquet-branch-c32-b0-em-v4-2026-06-19 | #389 | prose (#425 sentence "CR3BP/BCR4BP have no epoch and no flyby") |
| 52 | floquet-branch-c32-c-3.1774-em-v4-2026-06-19 | #392 | prose (same) |
| 53 | floquet-branch-c11a-b0-em-v4-2026-06-19 | #393 | prose (same) |
| 54 | floquet-branches-em-cycler-nodes-v4-aggregate-2026-06-19 | #392 | prose (same) |
| 55 | closer-sweep-v1-russell-ocampo-3.1.2+1-2026-06-17 | #365 | prose (same) |
| 56 | closer-sweep-v1-russell-ocampo-4.3.1-5-2026-06-17 | #365 | prose (same) |
| 57 | closer-sweep-v1-russell-ocampo-4.5.2-2-2026-06-17 | #365 | prose (same) |
| 58 | closer-sweep-v1-mcconaghy-2006-em-k2-s1l1-2026-06-17 | #365 | prose (same) |
| 59 | cross-system-se-em-l2-patched-cr3bp-2026-06-20 | `scripts/run_405_cross_system_search.py`, `analyze_405_theta_closure.py` | prose (bicircular named as an untried venue) |
| 60 | bcr4bp-phase-b-em-libration-seed-reach-spike-2026-06-20 | `scripts/spike_412_bcr4bp_reach.py` (#412) | INTEGRATED: method-invalid |
| 61 | cross-system-se-em-3d-patched-cr3bp-2026-07-01 | `scripts/run_515_cross_system_3d_search.py` | prose (future-work resweep condition) |
| 62 | cross-system-se-em-3d-asymmetric-patched-cr3bp-2026-07-01 | `scripts/run_517_asymmetric_3d_search.py` | prose (same) |
| 63 | cross-system-se-em-3d-multirev-patched-cr3bp-2026-07-01 | `scripts/run_516_multirev_3d_search.py` | prose (same) |
| 81 | mars-phobos-deimos-symmetric-closure-609-2026-07-16 | `scripts/_apply_609_mars_phobos_deimos_empty_region.py` | prose (literature-gap remark) |
| 83 | cross-system-se-em-l1l1-patched-cr3bp-2026-07-17 | `scripts/run_622_em_l1_se_l1_search.py` | prose (untried venue) |
| 94 | sunmars-bct-wsb-quasicycler-2026-07-22 | `scripts/search_681_sunmars_wsb_chain.py` | prose (`core/sunmars_wsb.py` does not import either model; concept lineage from #378 only) |

Two of the 19 are method-invalid (lines 39 and 60); 17 merely mention the model; none
unknown. For line 39, section 3.1 adds a COMPUTED fact: in the corrected model the same
constructor at Sun phase 0 no longer has an apoapsis inside 70 days, while at other Sun
phases it does. Neither the entry nor `search/cislunar_bct_search.py` sets a Sun phase
(grep, READ), so INFER the #378 sweep ran at the default phase 0 only; a rerun would have to
cover Sun phases.

### 5.2 Task conclusions

| Task | Conclusion (READ) | Depends on the defective model? |
|---|---|---|
| #292 | Bicircular Phase 1: equations, STM, propagator, corrector; CR3BP limit to float precision; "POL1 seed closes to a nearby BCR4BP L1 substitute". | Partly. The CR3BP limit and STM checks are sense-independent. The Sun-on closure is in the old model, and (section 3.2) the seed was not near L1 in this frame. |
| #303 | mu_sun continuation of the L1 Lyapunov (C 3.1294) to the full Sun; 50/50 steps; IC barely moves; probe 0 candidates. | Yes, for every Sun-on number (`data/bcr4bp_l1_family_303.jsonl`, `data/bcr4bp_validation_bridges_303.jsonl`). |
| #304 | Same for the Howell L1 southern halo; z0 -0.0529 to -0.0566; probe 0 candidates. | Yes (`data/bcr4bp_halo_family_304.jsonl`, `..._validation_bridges_304.jsonl`). |
| #334 | System-swap scaling exponent k = 2.89 +- 0.27 across 8 Sun-planet-moon systems; outliers flagged. | Partly: k and the outlier flags come from Sun-on continuations (yes); the constants table and corpus anchors do not. Whether k survives is untested. |
| #412 | EM-L1 seeds keep Earth reach 0.9 to 1.16 LD against a 3.9 LD target; wrong vehicle; scoped negative. | Yes for the numbers (empty-region line 60). The structural argument (an EM-bounded orbit stays EM-scale) may hold; INFER, untested. |
| #533 | Built the coherent QBCP, "verified against circular BCR4BP limits and finite differences". | Yes: the model itself (#892); its checks could not see parity or alpha_6 errors. |
| #538 | SE-EM connection through QBCP tori; the arc #538 to #646 ends in a 166,016 km closure floor, "not a proof of non-existence". | Yes for #538 and #544 (defective QBCP, old seed helper). For #617 to #646, INFER yes: `scripts/run_538_qbcp_cycler.py` imports `search/variational_qbcp_torus.py`, whose grid copy also carried the Sun-only alpha_6 (section 2.3). |
| #544 | EM-L2 torus does not converge in the QBCP (3.4) while SE-L2 does (3.1e-5); diagnosis of a degraded bicircular seed; POL1 gap explained as two coefficient sets. | Yes. Both models and the old seed helper; the POL1-gap explanation is already withdrawn (#892). |
| #593 | Scoped the #592 alpha_6 change: SE-L2 torus moved 0.194; multiple-shooting L1 substitute 4.14e-2 then 1.81e-2 from POL1; kept #592. | Yes, and withdrawn in part: the "fixed" model was still defective and #592 itself was wrong (#892). The 0.84 % size of alpha_6 - 1 is arithmetic on the series and stands. In the corrected model the harmonic-balance L1 orbit lands 1.7e-8 from POL1 (COMPUTED, section 3). |

### 5.3 Scripts that use either model and write results (not rerun)

| Script | Writes |
|---|---|
| `scripts/run_303_bcr4bp_l1_continuation.py` | `data/bcr4bp_l1_family_303.jsonl` |
| `scripts/run_303_catalogue_validation_probe.py` (reads the 303 family) | `data/bcr4bp_validation_bridges_303.jsonl` |
| `scripts/run_304_bcr4bp_halo_continuation.py` | `data/bcr4bp_halo_family_304.jsonl` |
| `scripts/run_304_catalogue_validation_probe.py` (reads the 304 family) | `data/bcr4bp_halo_validation_bridges_304.jsonl` |
| `scripts/scan_313_sun_jupiter_moons.py` | `data/scan_313_sun_jupiter_europa.jsonl`, `..._io.jsonl` (the writer of `data/scan_313_mars_*.jsonl` was not traced) |
| `scripts/scan_334_bcr4bp_system_swap.py` | `data/scan_334_bcr4bp_system_swap.jsonl` |
| `scripts/run_538_qbcp_cycler.py` | `data/runlogs/run_538_qbcp_cycler.jsonl` |
| `scripts/screen_884_sun_forced_em_cyclers.py` | `data/found/884_sun_forced_em_cyclers/` (already marked `INVALID_WRONG_SUN_SENSE.md`) |
| `scripts/_apply_378_empty_region.py` (no model import; records the #378 verdict) | `data/empty_regions.jsonl` line 39 |

Model users that print only (no result file found): `run_305_bcr4bp_gauntlet_probe.py`,
`run_522_coherent_connection.py`, `run_533_qbcp_connection.py`,
`search_coherent_connections.py`, `analyze_593_qbcp_l1_substitute_reconciliation.py` (has
its own copy of the QBCP field), `spike_378_bct_apoapsis.py`, `spike_412_bcr4bp_reach.py`.
Four of them (`run_522`, `run_533`, `run_538`, `search_coherent_connections`) call the
SE-L2 seed helper with 5 samples, which in the corrected models stops at the 7.2e-4
truncation floor (section 2.1); a rerun needs 11.

## 6. Anything suggesting a core model is still wrong

Nothing found. COMPUTED evidence for both: the bicircular identity at mu = 0 holds to the
Kepler miss of the module's Sun mass and to 4e-12 without it; the coherent model at zero Moon
mass agrees with the Sun-Earth three-body flow to 7.5e-7 over 6 time units, its harmonic-
balance L1 orbit lands 1.7e-8 from the published POL1, and its SE-L2 torus rotation number
(0.162806 to 0.162809, two independent correctors) is within 7e-4 of the three-body value.
The one remaining physical-constant inconsistency is the bicircular Sun mass, which misses
`n_S^2 a_S^3 = 1 + m_S` by 2.33e-8 relative (already recorded in the #884 review, 11a); the
QBCP set satisfies it.

## 7. Other things noticed

- The Task 1 command takes 5 min 24 s after the repairs (COMPUTED, 6 workers); the
  slowest item is `test_variational_crnbp_torus.py::test_n1_1_cannot_represent_...` (211 s),
  unrelated to these models. The QBCP torus fixture costs 40 to 85 s per worker.
- `uv run ruff check .` reports four errors, all in another agent's uncommitted #895 files
  (`scripts/screen_895_titania_oberon_realeph.py`,
  `src/cyclerfinder/search/titania_oberon_realeph_895.py`); none in files touched here.
- The bicircular "POL1" seed frame (section 3.2) affects the docstrings of five test files.
