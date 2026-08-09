# #805 — STM-exact derivatives for the DA-HOTM Taylor lane: the as-scoped fix is a measured negative; the exact-derivative CHAIN-NEWTON ENDGAME is the fix that works

**Verdict: #805's goal achieved (P5g' landing 2.784e-4 → 4.12e-9 on this
machine, inside the ~1e-5 corrector basin, machine-dependence of the landing
eliminated) — but NOT by the mechanism the task proposed.** Substituting the
exact STM-derived first-order block into `taylor_single_rev`'s polynomial makes
the descent WORSE (measured, decisively); the exact-derivative machinery
instead enters as a multiple-shooting Newton endgame on the section chain
(`section_chain_newton`), which converges from the FD-floor landing to the TRUE
float-map fixed point at ~1e-11 chain residual.

All measurements 2026-08-09, macOS/Accelerate machine, scratchpad diagnostics
against `cyclerfinder.genome.da_hotm_backend`. Baseline reproduced first:
FD-scheme landing for the P5g' golden (arXiv:2509.12671 Table 3, seed 1e-3 off)
= **2.784031e-4**, bit-matching #804's note.

## Part 1 — the as-scoped change (exact affine block in the polynomial): NEGATIVE

Built `single_rev_stm`: the exact 2x2 section-map Jacobian
`d(x',xdot')/d(x,xdot)` assembled from (a) the lift derivative (`ydot0` from
the Jacobi constant: `d ydot0/dx = -(dUbar/dx)/ydot0`, `d ydot0/dxdot =
-xdot/ydot0`), (b) the 6x6 STM from `cr3bp.cr3bp_stm_eom` variational
propagation to the same same-sign y=0 crossing, (c) the first-order
crossing-time correction `delta_t = -(Phi L ds)_y / ydot_f`. Validated against
central FD of the float `single_rev` at three section points: agreement
~1e-8 RELATIVE (the FD probe's own error floor), return parity vs the float
oracle ~5e-11. The machinery is correct.

Wiring the exact constant + linear blocks into `taylor_single_rev` (degree >= 2
still FD-fitted, to the residual after subtracting the exact affine part):

- Standard P5g' config (order=2, h=3e-4): landing **9.2e-4** (worse than
  2.784e-4), float P^5 residual there +inf (compose walls off).
- Full trace without the stagnation break: the re-expansion iteration
  WANDERS (dist 1.0e-3 → 1.8e-3 → 4.3e-3) then falls off the trust region.
  Genuinely non-convergent, not a break-tuning artifact.
- Own (order, h) sweep for the new scheme (order 1-3, h 3e-4..1e-2): best
  9.2e-4, i.e. every configuration is worse than the FD scheme's 2.784e-4.
- 2x2 block decomposition (exact-const-only / exact-lin-only / exact-both vs
  FD-both): 8.8e-4 / 9.9e-4 / 2.8e-3 vs 2.784e-4 — EACH exactness injection
  degrades the descent.

**Mechanism.** The FD least-squares "first-order block" was never a noisy
estimate of the Jacobian: at the golden, the fitted linear coefficients differ
from the exact Jacobian (212.16, 83.16; 541.23, 212.16) by 13.3 / 8.0 / 89.7 /
58.8 — **6-25%**. Under this map's violent curvature (the order-2 polynomial
only tracks the map to ~1e-2 across the tiny [-3e-4, 3e-4]^2 domain), the
regression affine block is a domain-AVERAGED linearization, and the
out-of-domain `compose_self` extrapolation that drives the descent (evaluating
the polynomial at chain points up to ~0.4 away) depends on exactly that
averaging. #804's "(order, h) sweep confirms the shipped baseline is optimal"
finding was a property of the REGRESSION scheme; the exact-Taylor scheme is a
different scheme, and it is uniformly worse at this descent task. The #805
premise ("FD noise in the linear block causes the floor; exact derivatives cut
it") conflated the FD noise floor (machine-dependence, real) with the
truncation-artifact landing offset (dominated by the scheme's own out-of-domain
composition, which exactness makes worse, not better).

`taylor_single_rev` is therefore RESTORED bit-for-bit to the FD fit, with the
finding documented in its docstring.

## Part 2 — the fix that works: exact-derivative section-chain Newton endgame

Direct Newton on `P^5` with the exact chained Jacobian fails (measured: from
the 2.784e-4 FD landing, |r| creeps 0.38 → 0.34 while the iterate wanders to
1.1e-2 and the chain walls off; the residual field is the #804-documented
pathology, and |J_P5 - I| ~ 3.6e3 compounds the instability). But MULTIPLE
SHOOTING over the n chain nodes — solve `P(s_i) = s_{i+1 mod n}` for all n
nodes simultaneously, each 2x2 block Jacobian exact from `single_rev_stm`,
backtracking line search on the chain-residual norm — distributes the
instability across per-rev blocks (|J_rev| ~ 5e2) and converges cleanly:

    from the FD landing (2.784e-4 off, chain residual 0.385):
    |R|: 3.8e-1 → 2.8e-1 → ... → 1.7e-2 → 4.0e-4 → 7.3e-7 → 1.8e-11 (9 iters)
    final node0 distance to published golden: 4.1e-9

which is consistent with #804's item-1 measurement that the true float-map
fixed point is ~1e-9 from the published (rounded) golden. The landing is now
the TRUE fixed point, so the machine-dependence of the landing (3e-5 Linux vs
2.8e-4 macOS — different noise draws of the same FD floor) collapses to
integrator-level agreement whenever the endgame accepts.

Shipped as `DASectionMap.section_chain_newton` + wired into
`taylor_fixed_point` (new `refine: bool = True` parameter) with a STRICT
fallback: any failure (initial chain does not exist — the off-family fail-fast,
measured to fire immediately at 5e-4 and 1e-3 off-family; singular system;
best residual > 1e-9) returns `None` and the caller keeps the FD-floor landing
unchanged. The endgame can only improve a landing, never worsen one.

End-to-end from the standard coarse seed (1e-3 off, the genome-test/lane-test
configuration): landing distance **4.120885e-9** (was 2.784031e-4, ~67,000x),
float P^5 residual at the landing **4.2e-10** (was 0.38 — the pre-#805 landing
was a truncation artifact, not a near-fixed-point; the new landing is the real
thing). Runtime unchanged (41.9 s vs 41.0 s baseline; the polynomial descent
dominates, the endgame adds ~1 s).

## Scope / seam notes

- The corrector (Task 5, `search/da_hotm_close.py`) REMAINS the certified
  closure authority; the endgame hands it a candidate at the radius-0 start.
  The micro-multistart is retained as the robustness net for endgame-declined
  coarse candidates. Docstrings updated (`da_hotm_close`, `da_hotm_enumerator`,
  backend class/methods).
- Higher-order exactness / a fuller DA implementation: still out of scope, and
  Part 1 is direct evidence that partial exactness inside the truncated
  polynomial is actively counterproductive for this descent device — a future
  full-DA lane should replace the whole polynomial stage, not patch blocks of
  the regression fit.

## Tests (tightened, none weakened)

- `tests/genome/test_da_section_map.py`:
  - `test_single_rev_stm_jacobian_matches_central_fd` (NEW): exact Jacobian vs
    central FD probe at 1e-6 relative (probe-limited), return parity 1e-9.
  - `test_section_chain_newton_declines_off_family` (NEW): pins the
    strict-fallback contract (`None` where the chain does not exist).
  - `test_taylor_fixed_point_reaches_png_neighbourhood`: bound TIGHTENED from
    "2x closer than the seed" (#804's honest FD-floor form, 5e-4) to
    **dist < 1e-5** — the corrector's demonstrated reliable basin, i.e. the
    #805 capability claim, with ~2400x observed headroom — plus
    `residual(fp, 5) < 1e-8` (genuine near-fixed-point, machine-independent)
    and a moved-off-the-seed non-circularity guard (replacing `dist > 1e-7`,
    which the endgame's genuine convergence now correctly violates).
- `tests/search/test_png_lane_recovery.py` (fast lane proof): coarse-candidate
  bound TIGHTENED 1e-3 → **1e-5**; non-circularity guard likewise converted to
  moved-off-the-seed form; the 1e-11 closure + published-value asserts
  unchanged. Slow global-sweep test: assertions unchanged (docstring updated),
  verified green in this session (~10 min serial run).

## Invalidation check ("bug-fix invalidates past searches")

No past search result flows through the changed code path. The #527/#532
campaign scripts (`run_527_hilda_dahotm_search.py`, `run_532_*.py`) use
`search/da_hotm_enumeration.py`, which imports `SamplingSectionMap` +
`enumerate_fixed_points` — the float sampling backend, NOT
`DASectionMap.taylor_fixed_point`. (#527's certified orbits closed to ~1e-14
via the corrector; its negative verdict is structural — stable resonant
backbone vs unstable-manifold transport — untouched by candidate precision.)
`taylor_fixed_point` is exercised only by the genome/lane tests. #804's caveat
that the negative-registry entries are implicitly machine-conditional for this
lane is now moot going forward: endgame-accepted landings are
machine-independent.

## Follow-ups registered

- **#812** — point the now-precise Taylor+endgame lane (`taylor_enumerate`
  with `refine=True`) at a fresh multi-rev band as a discovery probe; the
  pre-#805 lane was demonstrably unable to land strongly-unstable basins
  without the corrector multistart, so bands screened only by the sampling
  backend may hide needle-basin families the refined lane can now surface.
  Speculative, low priority.
