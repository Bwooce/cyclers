# Shooter performance pass — parallel FD Jacobian + variational verdict (#159)

Date: 2026-06-07. Rank-1 next rung from the Phase C verdict
(`2026-06-07-nbody-shooter-phase-c-verdict.md` §"Next rungs" item 1). Territory:
`nbody/shooter.py`, `nbody/propagator.py` (additive), `tests/nbody/`.

## The cost model being attacked

The multiple-shooting LM solve (`shooter.py::shoot`) builds its Jacobian by
finite difference over REBOUND. Each Jacobian is `3*n_nodes+1` residual
evaluations; each residual re-propagates the full multi-segment trajectory
serially. The verdict note's cost model: `n_LM_steps × (3·n_nodes+1) ×
n_segments × t_segment` ≈ 0.5–5 CPU-h/member.

Two independent levers were scoped. Lever 1 shipped; lever 2 is a documented
blocker (below) and was NOT forced, per the task's stop rule.

## Lever 1 — parallel finite-difference Jacobian columns (SHIPPED)

The `3*n_nodes+1` FD residual evaluations are independent → evaluated across a
process pool. Wired as opt-in `n_jobs` on `shoot()`:

- `n_jobs == 1` (default): today's serial path EXACTLY. `least_squares` builds
  its own internal-FD Jacobian; no behavioural change, byte-identical, golden
  gates untouched.
- `n_jobs > 1`: an explicit forward-difference Jacobian whose columns are mapped
  over a `ProcessPoolExecutor`, reusing the scan rung's primitives-only worker
  contract (`search/scan.py`): only the picklable `ShootingSeed` (numpy arrays +
  floats) and the perturbed free-vector cross the process boundary; each worker
  constructs its own `Ephemeris(model=...)` once via the pool initialiser — never
  pickling a live `Ephemeris`. The FD step matches scipy's `lm` convention
  (`h = sqrt(eps)·|x|`), so the parallel Jacobian equals the serial one of the
  same scheme to working precision.

### Why processes, not threads (MEASURED, not assumed)

The rails third-body perturbation is installed as a Python `additional_forces`
CFUNCTYPE callback (`propagator.py::_install_rails_forces`), invoked on **every**
IAS15 force evaluation. REBOUND re-acquires the GIL inside that inner loop, so a
thread pool would serialise on the callback — the GIL is held during the
dominant work. Process pool it is.

### Correctness (asserted)

`tests/nbody/test_shooter_parallel_jac.py`:
- `test_parallel_fd_jacobian_equals_serial`: builds the FD Jacobian serially and
  via the pool on a small two-body fixture; asserts
  `allclose(rtol=1e-9, atol=1e-12)`. **Measured: identical to 0.00e+00 relative
  difference** (same arithmetic, different execution mode).
- `test_shoot_n_jobs_reduces_defect`: an end-to-end bounded `shoot(n_jobs=4)`
  drives the defect below the seed (solver health). NB: a full `shoot` with
  `n_jobs>1` swaps scipy's *internal*-FD `lm` for our *explicit*-FD `lm` (a
  different FD scheme from `n_jobs==1`), so the LM iterates legitimately differ
  at FP level between the two `shoot` paths — the byte-identical guarantee is on
  the Jacobian-of-the-same-scheme (test 1), and on the default `n_jobs==1` path
  vs pre-#159.

### Benchmark (MEASURED, not projected)

16-core host (`nproc`=16). Two fixtures, both with DE440 rails perturbers:

Light fixture — 3-node E-M-E, ~4-month legs, 18 FD columns:
```
single residual eval : ~0.06 s
SERIAL  Jacobian      : 1.09 s
PARALLEL n_jobs=4     : 0.89 s  (1.23x)
PARALLEL n_jobs=8     : 0.90 s  (1.22x)
PARALLEL n_jobs=16    : 1.38 s  (0.79x — pool overhead dominates cheap legs)
```
Pool spawn + dispatch overhead swamps the gain when legs are this cheap.

Heavy fixture — 4-node E-V-E-M, multi-year legs (V/E/M perturbers), 24 FD
columns — representative of the 0.5–5 CPU-h members the cost model targets:
```
single residual eval : 5.09 s
SERIAL  Jacobian      : 130.34 s
PARALLEL n_jobs=8     : 22.17 s  (5.88x)   max_rel_diff vs serial = 0.00e+00
PARALLEL n_jobs=16    : 20.11 s  (6.48x)   max_rel_diff vs serial = 0.00e+00
```

**Headline: 5.9–6.5x measured per-Jacobian speedup on the representative heavy
member, with a byte-identical Jacobian.** The speedup is bounded by
`min(n_jobs, 3·n_nodes+1)`; the ~6.5x at 16 workers / 24 columns reflects
per-column propagation-cost imbalance (the columns are not equal-wall) plus the
fixed pool overhead. The lever pays off exactly where it matters — long, hard
legs — and is a small net loss on trivially cheap legs, hence the opt-in default
of `n_jobs=1`.

## Lever 2 — REBOUND variational equations (BLOCKER — NOT shipped)

REBOUND 5.0.0 supports first-order variational particles natively
(`sim.add_variation()`, IAS15-compatible). One augmented propagation would yield
the state-transition matrix and replace the FD sweeps. **It is not cleanly
compatible with our planets-on-rails additional-force setup.** The specific
blockers, established by direct probe (not assumption):

1. **REBOUND does not auto-differentiate `additional_forces` into the
   variational equations.** REBOUND's variational integrator differentiates only
   the forces it computes analytically (its own N-body gravity). Our rails
   third-body perturbation is a black-box Python callback; its Jacobian (the
   tidal matrix `mu_p·(3 d⊗d/|d|^5 − I/|d|^3)`, with `d = r_p − r`) is **never
   applied to the variational particle**. Probe result: with the perturber
   active but the callback touching only the real particle, the STM column
   `∂r_term/∂v0` shifted only ~320 km from the pure-Sun value — that residual is
   the *indirect* effect (the real path bent, so the variational eq sees a
   shifted Sun-gravity gradient); the perturber's *direct* tidal contribution to
   the STM is missing entirely. The STM REBOUND returns under our current setup
   is therefore wrong.

2. **The variational particles are not reachable from inside the public
   `additional_forces` callback.** `sim.particles` inside the callback exposes
   only the real particles (`N=2`); `sim.particles[2]` raises
   `AttributeError: Index 2 ... out of range`. The variational particles live in
   a separate `particles_var` buffer requiring lower-level ctypes access that
   REBOUND neither documents nor tests for the additional-force path (none of
   REBOUND's `test_variational.py` cases use `additional_forces`).

To make lever 2 correct we would have to: add 6 first-order variational
particles per leg (the STM w.r.t. the 6-component spacecraft initial state),
reach them via raw ctypes inside the force callback, and hand-code the rails
tidal matrix applied per variational particle — including the softening-clamp
discontinuity (`d_eff = max(|d|, d_safe)`), where the gradient is non-smooth.
This is exactly the part REBOUND does not support out of the box; getting the
variational-particle wiring wrong produces a *silently* wrong Jacobian (a
plausible-looking STM that omits a physical term — see blocker 1). The risk/
reward did not clear the bar against a shipped 6x from lever 1, so per the task's
explicit stop rule it was not forced. If revisited, the analytic rails Jacobian
itself is the easy part; the REBOUND ctypes variational-particle plumbing and a
golden STM-vs-FD validation gate are the work.

## Gate status

- Full non-slow nbody suite: **33 passed** (was 31; +2 new parallel-jac tests),
  golden gates green, default `n_jobs=1` path byte-identical to pre-#159.
- Parallel Jacobian vs serial: identical to 0.00e+00 relative difference
  (asserted).

## Commit

- `nbody/shooter: parallel finite-difference Jacobian via opt-in n_jobs`
  (+ this note). `propagator.py` needed no change — the rails-cache reuse it
  already exposes is what makes the per-column propagation affordable.
