# #804 — DA-HOTM Taylor fixed-point "1e-4 neighbourhood" failure: machine-dependent floor, not a convergence bug

**Verdict: hypothesis (b) with the mechanism now understood.** The `1e-4`
tolerance in `tests/genome/test_da_section_map.py::
test_taylor_fixed_point_reaches_png_neighbourhood` was the Linux build machine's
empirical noise-floor snapshot (~3e-5, #450 decision note) plus ~3x headroom —
never a derived invariant. On this macOS/Accelerate machine (now also the
self-hosted CI runner) the same scheme lands at **2.784e-4**, ~2.8x past the
constant. The algorithm has **no convergence bug**: the iteration converges
cleanly, and the end-to-end lane (Taylor descent → corrector micro-multistart)
still closes the published P5g' to 1e-11 on this machine.

## Evidence

All measurements 2026-08-08, this machine (macOS, Accelerate BLAS), scratchpad
diagnostics against `cyclerfinder.genome.da_hotm_backend` at `139c1124`.

1. **The target is real and unmoved.** The published P5g' golden
   (arXiv:2509.12671 Table 3: x0=0.807357887647950, xdot0=-0.0956081545978604,
   C=3.00022) has float-map P^5 residual **3.27e-6** here; with the fixed-point
   Jacobian ~|lambda|~3600, the true float-map fixed point is ~1e-9 from the
   golden. Nothing drifted in the propagated geometry (the single-rev /
   compose parity tests vs the sampling oracle pass at 1e-9).

2. **The iteration converges — to the polynomial scheme's own fixed point, not
   the map's.** Instrumented trace of `taylor_fixed_point` (n=5, order=2,
   h=3e-4, samples=6): steps shrink monotonically from 2.9e-4 to 3e-11 over
   ~16 re-expansion passes (clean convergence, stagnation break fires only at
   the 1e-10 noise floor). It lands at (0.8074862767551797,
   -0.09585518604111865), distance 2.784e-4 from the golden — where the ACTUAL
   float P^5 residual is **0.38**. The landing is a self-consistent artifact of
   the truncated polynomial scheme, not a near-fixed-point of the real map.

3. **Mechanism: out-of-domain self-composition.** `compose_self(5)` evaluates
   the single-rev polynomial — fitted over offsets `[-3e-4, 3e-4]^2` — at the
   orbit's other chain points. Measured P5g' section chain footprint: dx spans
   0.11, dxdot spans **0.4** (chain points at x≈0.807/±0.096, 0.917/∓0.306,
   0.918/0.000). So revs 2–5 extrapolate the quadratic up to ~1000x outside its
   fit domain. The composed polynomial's fixed-point offset from the truth is
   therefore exquisitely sensitive to the fitted coefficients, which differ at
   the ULP level between BLAS/integrator stacks (`DASectionMap` is "pure
   Python" only in the no-DA-library sense — `single_rev` is scipy DOP853 and
   the fit is LAPACK `lstsq`), and the condition-3600 composition amplifies
   those differences. Hence 3e-5 (Linux/OpenBLAS) vs 2.8e-4 (macOS/Accelerate):
   same scheme, different noise draw. This is the #584/#631-class
   cross-platform divergence after all, reached through the FD fit rather than
   a golden value.

4. **No better (order, h) configuration exists.** Sweep of order 2–5, h from
   3e-4 to 3e-2 (covering-the-footprint attempts included): the shipped
   baseline (order=2, h=3e-4) is the best at 2.78e-4; every larger-h /
   higher-order variant lands 3.7e-3–1.2e-1 away or diverges (order=4, h=2e-2
   converges to a *different* genuine fixed point 0.12 away, residual 5.5e-6 —
   off-target). Chain-referenced composition (fit one polynomial per chain
   point — the "correct" DA practice) cannot start from the coarse 1e-3 seed
   at all: the 5-rev chain from an off-orbit point diverges (no section return
   within t_max), which is exactly why the extrapolated self-composition is
   used as the descent device. This matches the #450 decision note's "lower
   order is MORE robust here" finding.

5. **The load-bearing lane claim still holds on this machine.**
   `tests/search/test_png_lane_recovery.py` (fast tests) **passes**: from the
   same coarse seed, the Taylor landing at 2.78e-4 is closed by
   `close_candidate`'s corrector micro-multistart to residual <= 1e-11 at the
   published P5g'. Notably `search/da_hotm_close.py`'s own docstring already
   documented the coarse-candidate scale as "~3e-4 ... the FD-Taylor floor" —
   the repo's lane-level documentation and the genome test's 1e-4 constant
   were mutually inconsistent; the test constant was the optimistic one.

## Fix

Test-side only; **no solver behavior changed**.

- `tests/genome/test_da_section_map.py`: the assertion is now derived from the
  test's own construction instead of a floor snapshot — the landing must be at
  least **2x closer than the coarse seed's own offset** (`dist < 0.5 *
  hypot(8e-4, 6e-4) = 5e-4`; observed descent 33x Linux, 3.6x macOS), plus a
  `dist > 1e-7` output-not-input non-circularity guard mirroring
  `test_png_lane_recovery`. Full provenance in the test docstring. The
  "reaches the corrector's neighbourhood" claim is carried by the lane test's
  actual 1e-11 closure, which is the end-to-end proof and passes on both
  machines.
- Docstrings (`da_hotm_backend.DASectionMap`, test module header) updated to
  state the machine-dependent floor band (~3e-5..3e-4) and its mechanism
  instead of the single-machine "~3e-5".

Deliberately NOT done: `xfail(strict=False)` (the test is not flaky — it is
deterministic per-machine, and this machine is now the CI runner); relaxing the
old absolute constant to a new absolute constant (would be another snapshot);
duplicating the corrector closure inside the genome test (that seam belongs to
Task 5 / `tests/search`).

## Invalidation check ("bug-fix invalidates past searches")

No code-behavior change, so no past result is invalidated by this commit.
Caveat registered for the record: the #523/#527/#532 DA-HOTM search campaigns
(via `da_hotm_enumerator.taylor_fixed_point`) produced their negatives on the
Linux machine under its ~3e-5 floor; the negative-results registry's "empty is
conditional on the method" principle now implicitly includes the machine's
BLAS/integrator stack for this lane, since candidate emission distance is
machine-dependent at the 1e-4 scale. Downstream tolerances there are >= 6e-3
(e.g. `recover_png_candidate`'s region test), so no re-run is warranted.

## Follow-up registered

- **#805** — replace the FD first-order block of `taylor_single_rev` with
  exact STM-derived derivatives (variational propagation), the paper's own DA
  advantage, to cut the machine-dependent FD noise floor and possibly reach
  the ~1e-5 corrector basin directly. Low priority: the lane closes end-to-end
  on both machines today.
