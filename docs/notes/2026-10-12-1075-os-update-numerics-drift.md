# #1075: twelve CI failures after the macOS 27.0.1 update (pre-registration, 2026-10-12)

Owner: numerics-drift-opus. Brief (lead, 2026-10-12 14:30 AEDT): CI run 38180601249 (push 983ba9c1, suite
12:34-13:16 AEDT on this Mac) failed 14 tests; the last green run was ada736b9 (11 Oct evening). The Mac
installed macOS 27.0.1 at 09:38 AEDT on 12 Oct. earthmoon-1054-opus proved the OS cause for the two
`#1054` tests (exact-float goldens; note `docs/notes/2026-10-11-1054-seed-design-fix.md` sec. 10). This
note covers the other twelve. Binding rule: no deselect, skip, xfail, slow mark, underived tolerance, or
re-recorded expected value from the code under test.

## 1. Method (written before any fix)

1. Read each test and the code under test; classify a priori as exact-float golden, tolerance edge,
   convergence of an iterative or chaotic solver (known hit, saddle branch), or time-budgeted.
2. Check that no commit since ada736b9 touches the modules under test.
3. Reproduce each test alone in the main checkout (`pytest -n 0`, one at a time, load recorded).
4. Reproduce in a detached worktree at ada736b9 with its own venv. A failure there today, where it
   passed on 11 Oct, is platform or contention, not code.
5. Separate contention from arithmetic for any time-budgeted test (a quiet-slot re-run, or the budget's
   clock measured).
6. Fix per class. A tolerance edge gets a tolerance derived from the solver's conditioning or residual
   floor. A known hit that no longer converges is re-found under today's arithmetic and its identity is
   checked against the stored values. A src change is stated here first, with its expected effect, and
   kept minimal.

## 2. Step 2: code changes since ada736b9

`git log ada736b9..HEAD -- src/cyclerfinder/search src/cyclerfinder/genome tests/search tests/genome`
lists six commits. They touch only `search/varin_closed_families.py` (new, #975), `search/literature_check.py`
(#1052 Pluto anchors) and their own tests, plus two literature-gate tests (#1066). None of the twelve tests
imports these modules. `uv.lock` and `pyproject.toml` are unchanged; numpy 2.4.6 and scipy 1.17.1 in both
venvs. numpy links the system BLAS/LAPACK (Accelerate), which the OS update replaced.

## 3. Step 1: a priori classification

None of the modules under test reads a clock (grep for perf_counter, time.time, monotonic, deadline,
budget, timeout: no hits). The only time limit is pytest-timeout (600 s, 1200 s on one test), and every
CI failure is an assertion failure, not a timeout. So no test is time-budgeted, and CPU contention cannot
change any of these results: the arithmetic is single-threaded and does not depend on wall time.

| # | Test | Class (a priori) | Mechanism |
|---|---|---|---|
| 1 | genome/test_asymmetric_branch_gram_schmidt::test_phase1_branch_at_saddle_center_survives_gram_schmidt | chaotic solver convergence (saddle branch) | `branch_at_saddle_center` i=124: damped full-asymmetric Newton from parent + eps*v; outcome already documented as BLAS-dependent (#740) |
| 2 | genome/test_asymmetric_branch::test_branch_at_saddle_center_i124_eps_5e_minus_4_converges | same as 1 | same call |
| 3 | genome/test_floquet_phase1_reproduction::test_phase1_end_to_end_reproduction | same as 1 | same call via the detector bracket |
| 4 | search/test_free_return_chain::test_degenerate_reduces_to_single_free_return_ellipse | tolerance edge | turn angle by arccos of a cosine one ulp from 1 |
| 5 | search/test_jovian_resonant_connections::test_find_heteroclinic_returns_known_combo | chaotic solver convergence (known hit) | 8x8 phase grid seeds a 2x2 FD Newton on a manifold section map |
| 6 | search/test_jovian_resonant_connections::test_table3_gate_honest_result | same as 5 | same seed and Newton |
| 7 | search/test_saturn_titan_resonant_connections::test_attempt_chain_closure_natural_multiple_shooting_progress_wrong_energy | chaotic solver walk | asserts properties of a 40-iteration NON-converged multiple-shooting walk |
| 8 | ...::test_find_homoclinic_ranks_by_residual | chaotic solver convergence (known hits) | grid-seeded Newton, tol 1e-7, max_iter 20 |
| 9 | ...::test_find_symmetric_chain_seed_finds_new_closer_near_perpendicular_candidate | chaotic solver convergence (scan) | equal-index scan plus Newton |
| 10-12 | ...::test_known_hit_{converges_with_real_ghost_margin, forward_backward_reapproach, independent_radau_crosscheck}[tertiary] | to be decided by probe | Newton from the stored phases of the tertiary hit |

No module under test uses a process or thread pool either (grep for Pool, concurrent, joblib,
cpu_count, multiprocessing: no hits), so the worker count cannot depend on the machine.

## 4. Step 3: reproduction alone in the main checkout (2026-10-12, pytest -n 0, one at a time)

Log: `scratchpad/numerics-drift-opus/main/summary.txt` (session scratchpad; the values are copied here).

| # | Result | Wall | Load (1-min) | Failure value (same as CI) |
|---|---|---|---|---|
| 4 | FAIL | 5 s | 36 | intermediate_turn_deg 1.2074182697257333e-06 < 1e-06 |
| 1 | FAIL | 233 s | 37 | branch_at_saddle_center returned None |
| 2 | FAIL | 222 s | 37-40 | returned None |
| 3 | FAIL | 246 s | 40-45 | P1.4 returned None |
| 5 | FAIL | 195 s | 40-50 | len(candidates) 0 == 1 |
| 6 | FAIL | 258 s | 40-50 | known-hit combo did not reach tol |
| 7 | FAIL | 547 s | 45 | abs(3.0099370969589163 - 3.01) = 6.2903041083473e-05 > 0.001 |
| 8 | FAIL | 92 s | 45 | len(candidates) 0 >= 1 |

(9-12 and the ada736b9 worktree runs: sec. 5.)

## 5. Step 4: worktree at ada736b9 (own venv, numpy 2.4.6, scipy 1.17.1)

| # | Result at ada736b9 today |
|---|---|
| 4 | FAIL, same value 1.2074182697257333e-06 |

(the others are queued behind the main-checkout runs; filled in when they finish)

## 6. Diagnosis per class (probes `scripts/probe_1075_saddle_branch.py`, `scripts/probe_1075_connection_hits.py`)

### 6.1 Test 4 (free-return degenerate): arccos resolution floor
`_intermediate_turn_geometry` computes the turn as `arccos(v_in . v_out / (|v_in| |v_out|))`. For the
degenerate chain the two arcs differ by 4.3e-12 AU in a and 3.3e-12 in e (a joint least-squares solve), so
the true angle is 6.9e-10 deg (by atan2 of cross and dot). Today the cosine evaluates to 1 - 2.2e-16 (one
ulp below 1, also in pure Python), and arccos of that is sqrt(2 * 2^-52) rad = 1.2074182697257333e-06 deg,
which is exactly the failure value. Before the update the cosine rounded to 1.0 or above (clipped), giving 0.
The test's 1e-6 deg bound sits just below arccos's resolution near 0, so it was always one ulp from failing.

**Fix (src, pre-registered):** compute the turn as `atan2(|v_in x v_out|, v_in . v_out)`, which is well
conditioned at every angle. Expected effect: the degenerate turn becomes about 7e-10 deg; any other turn
changes by at most about 1e-16 / sin(turn) rad, so no `turn < max_turn` gate changes (tests at
test_free_return_chain.py lines 143 and 169, test_continuation_chain.py line 208). The test is unchanged.
`continuation_chain.py` line 613 uses the same arccos pattern; it is not failing and is left alone
(recorded here as a latent edge).

### 6.2 Tests 1-3 (saddle-center branch at i=124, eps=5e-4): the recorded branch was never local
- The eigenvalue picked today is 0.9742446738323984, the documented old-Mac value (0.97424467383240) to
  every digit the #740 docstring prints. The old Mac converged (to (4,3)); today the whole ladder diverges.
  So the platform sensitivity is in the Newton walk, not in the eigendecomposition the #740 docstring blames.
- The seed is not near a solution. Parent i=124 closes over the full period to 1.4e-7. Parent + eps*v
  closes to |R| = 2.30 at eps = 5e-4. The residual is in-plane and scales as eps^2 (|R|/eps^2 = 4.67e7 at
  eps = 1e-5, 3e-5 and 1e-4; it saturates above), the second-order out-of-plane coupling amplified by an
  STM whose largest singular value is 8.8e7. The linear regime needs eps far below 1e-9. Every rung of the
  ladder (1e-4 to 1.5e-3) therefore starts Newton O(0.5-2) from any periodic orbit; the first capped
  Newton step moves y, xdot and ydot by the full 0.1 cap.
- Where Newton ended is therefore a property of an 80-iteration walk. The committed record
  (`data/floquet_phase1_reproduction.jsonl`, Linux/OpenBLAS) is planar (z0 = 1.7e-22, zdot0 = 8.3e-24),
  T = 23.36 against the parent's 17.46, C = 3.7975 against 3.1418. The old-Mac landing was (4,3) at
  C = -0.753. Neither is a member of a family branching at this saddle-center: they are distant periodic
  orbits that a walk happened to reach.
- Consequence: the #347 Phase 1 PASS ("a branched orbit with topology distinct from the parent") is not
  evidence of a branch at the saddle-center. This affects #347 and is reported to the lead for a ruling.
  Re-closing the committed record would not restore a positive control: the record came from the code
  under test, not from a source.

### 6.3 Tests 10-12 (Saturn tertiary known hit): Newton tolerance below the map's noise
Newton from the stored phases (test settings: tol 1e-9, fd_step 1e-7) stops after 7 iterations at
|R| = 3.4e-9, 4.1e-9 from the stored crossing (-1.12647328, 0.13211509): the hit is there. The backtracking
cannot reduce |R| further. The noise of the section-gap map at the hit (13-point line through the phases,
straight-line fit per component, steps 1e-11, 1e-10, 1e-9; `probe_1075_connection_hits.py tertiary_floor`):

| Direction | Fit-residual std (x, xdot) |
|---|---|
| tau_u | (0.9-1.0e-9, 0.6-0.7e-9) |
| tau_s | (3.1-6.9e-8, 2.1-4.6e-8) |

The scatter is the same at all three steps, so it is noise, not curvature. The largest standard deviation
is sigma = 6.9e-8 (x, along tau_s). A residual norm below about sigma cannot be told from zero, so a
Newton tol of 1e-9 (and the test's `residual < 1e-8`) asks for luck.

**Fix (test only, pre-registered):** per-hit Newton tol and residual bound for the tertiary hit =
k * sigma with k = 3 chosen now: 3 * 6.9e-8 = 2.1e-7, used as 2e-7. The primary and mirror hits keep 1e-9 /
1e-8 (they pass, and their floors are not measured). Identity (crossing within 1e-5 of the stored values),
ghost margin (> 200 * GHOST_GUARD_DELTA), Radau cross-check (< 1e-6) and re-approach bounds are unchanged.
Expected: tests 10-12 pass with today's 3.4e-9 landing.

### 6.4 Tests 5-6 (Jovian Table-3 known hit): the hit exists; the grid-seeded walk misses it
The 8x8 phase grid picks the same cell as before, (4,1), with gap 0.175 against 0.185 next; 8 of 64 cells
reach the section. Newton from (4,1) now stalls at |R| = 2.1e-3, 0.026 from the stored crossing; the next
three cells fail too. Re-found leg by leg: each leg's crossing depends on one phase, so two 1-D scans
toward the STORED crossing (160 samples, then three 40-sample refinements) give tau = (14.36867, 7.27806).
Newton (test settings) from there converges in 2 iterations to |R| = 9.5e-8 at (-1.306870677409797,
0.02869078595377826), 3.1e-8 from the stored (-1.3068706902782006, 0.02869075783271665). The map's noise
at the hit is 0.2-2.7e-9 (`jovian_floor`), so tol 1e-7 is honest there.
Fix pending the lead's ruling (sec. 7).

### 6.5 Test 7 (Saturn natural multiple shooting, "wrong energy")
The test asserts that a 40-iteration NON-converged multiple-shooting walk drifts in Jacobi by more than
1e-3. The threshold was read off one walk (#782 note: 3.01 -> 2.976). Today the walk ends at
C = 3.0099370969589163, a drift of 6.29e-5. The claim the test exists for is that
`correct_multiple_shooting` does not hold C, so its result is not at the node's energy.

**Fix (test only, pre-registered):** bound the drift by the project's own same-energy tolerance, the
`jacobi_tol = 1e-6` that `correct_connection` uses to accept two orbits as the same energy. A drift above it
means the result is not at the node's energy by the project's own criterion; today's 6.29e-5 is 63 times
that. The integrator's Jacobi conservation at rtol 1e-12 is far below 1e-6, so the bound cannot be met by
integration error alone. The other assertions (progress below the seed residual, residual > 1e-3) are
unchanged.

### 6.6 Tests 8-9 (Saturn scans)
Per-combination runs behind the two scans (`probe_1075_connection_hits.py saturn_scans`, the tests' settings):

| Test | Combination | Today |
|---|---|---|
| 8 (branch +1, scan_n 8, tol 1e-7, max_iter 20) | (4,4) (4,5) (5,4) | stall at 2.4e-2, 6.5e-2, 6.3e-2 |
| 8 | (5,5), the primary hit | stalls at 1.25e-2 at (-1.6723, 0.0017), near the primary hit (-1.67320084, ~0) |
| 9 (branch -1, scan_n 12, tol 1e-9, fd_step 1e-6) | (4,4), the #782 near-perpendicular hit | stalls at 1.09e-9 against tol 1e-9, xv = (0.9492672648843903, -4.3e-11) |
| 9 | (4,5), (5,4) | converge (5.2e-10, 9.3e-10); not near-perpendicular, so the test rejects them |
| 9 | (5,5) | stalls at 1.2e-2 |

Test 8 is grid-seeded Newton, the same class as the Jovian tests (ruling, sec. 7). The primary hit
still converges from its stored phases (test_known_hit_*[primary] pass).

Test 9 is the tertiary class: the (4,4) hit is found and Newton stalls just above a tol that sits
below the map's noise. Measured at the stall phases (`symmetric_floor`): fit-residual std up to
sigma = 2.5e-9 (xdot, along tau_u; 0.4-2.5e-9 over both phases and all three steps).

**Fix (test only, pre-registered, same rule as sec. 6.3):** pass tol = 3 * sigma = 7.5e-9 to
`find_symmetric_chain_seed` in test 9. The test's own checks are unchanged (dist_to_65 < 0.02,
|xdot| < 1e-6, k_u == k_s, branch -1, residual < 1e-7). Expected: test 9 passes on the (4,4) hit.
The function's src default (1e-9) and the `near65_symmetric_seed_xv` fixture's 1e-9 (fd_step 1e-7;
it converges today) are below the same floor; they are left alone and recorded here as latent.

### 6.7 Speed of the runs
Wall times vary about fivefold with load (the same corrector call took 26 s at 14:40 and 5 s at
15:13). The worktree at ada736b9 and the main checkout give bit-identical values for the saddle
corrector call (`scratchpad/numerics-drift-opus/time_one.py`): lambda, v[2], residual and final x0
agree to every printed digit.

## 7. Open for the lead
1. Tests 1-3: the ruling on sec. 6.2 (narrow the claim and record the #347 consequence, or register a
   local branch method).
2. Tests 5, 6, 8, 9: all are grid-seeded Newton. Options: pin test 6's seed to the re-found phases (the
   identity check against the stored crossing stays); test 5 claims the scan itself finds the hit, so pinning
   narrows it. A project-wide fix is a leg-separable scan in `_scan_starts`: the gap is p_u(tau_u) - p_s(tau_s),
   so sampling each leg in 1-D at n^2 points costs the same 2 n^2 propagations as the n x n grid at n times
   the resolution per phase. It changes every caller's seed, so it needs a ruling.
