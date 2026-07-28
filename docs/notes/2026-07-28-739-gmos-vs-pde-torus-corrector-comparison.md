# #739 -- GMOS (`genome/qp_tori.py`) vs PDE (`search/variational_qp_torus.py`)
head-to-head, same parent orbit, same amplitude

**Task:** `#739`, a small, low-stakes, non-catalogue-affecting engineering-hygiene
follow-up to `#732`'s digest of Baresi, Olikara & Scheeres 2018 ("Fully
Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in
Astrodynamics"). That paper's own head-to-head found their preferred method
(GMOS -- stroboscopic-map shooting) faster, more accurate, and giving free
Floquet stability compared to a PDE(DFT)-class method resembling this
project's own pseudospectral corrector -- but only ever tested **stable**
parent orbits, never the ~1540x-monodromy-amplification regime this
project's PDE corrector (`#612`) was built to survive. `#732`'s digest flagged
"benchmark this project's own PDE corrector against its own GMOS corrector on
a shared low-amplitude torus" as an untested, cheap follow-up. This note is
that follow-up.

## Method

Reused `tests/search/test_variational_qp_torus.py`'s own **L2 positive
control** exactly (same halo builder `_l2_halo_at_315`, same
Neimark-Sacker-eigenvector seed, same GMOS call parameters: `k` from
`_best_k`, `n_trans=4`, `initial_torus_amplitude=5e-4`, `tol=1e-8`,
`max_iter=40`; same PDE bootstrap call: `discover_qp_torus(SYS, gmos, n1=10,
n2=4, max_nfev=300)`). Instrumented both stages with `time.perf_counter()`
instead of asserting on the numbers. System: Earth-Moon CR3BP, `mu =
0.012153643`, EM L2 halo at Jacobi `C=3.15`. Script:
`/private/tmp/.../scratchpad/bench_gmos_vs_pde.py` (one-off, not committed --
see "Why no permanent test" below).

## Premise correction: this "positive control" is NOT a Baresi-style stable
parent orbit

Before trusting the comparison, checked the parent orbit's own monodromy
spectral radius (the same quantity `test_l1_crosses_gmos_amplitude_wall`
already asserts `> 500.0` for the L1 wall case):

```
parent orbit: EM L2 halo at Jacobi=3.15, T_TU=3.410926
parent monodromy spectral radius: 1161.4
```

This is **violently unstable** -- in the same amplification regime as the L1
wall case (`~1540`), not the `|lam|=1` elliptic-center regime
Baresi/Olikara/Scheeres's own two test problems (ECEF circular orbit,
Earth-Moon DRO family) actually assume (their Eq. 1: "one pair of complex
conjugate eigenvalues with unitary modulus"). `#732`'s own digest language
("this project's own GMOS-lineage corrector... already converges cleanly...
documented as working below amplitude ~0.01") was imprecise on this point: GMOS
"works" here not because the parent orbit is dynamically stable, but because
the torus amplitude tested (`5e-4`) is tiny enough that the propagation
amplification hasn't yet degraded the residual past `tol` -- exactly the same
mechanism `#612`'s own docstring documents degrading monotonically with
amplitude on the L1 case. **Neither this project's L1 nor L2 EM halo at
C=3.15 is a genuine positive control for Baresi et al.'s own stable-orbit
comparison scope** -- both are hyperbolic parents; L2 is just being tested at
an amplitude far below where the amplification bites. No genuinely
Floquet-stable (`|lam| ~ 1`) parent-orbit torus case exists anywhere in this
project's test suite (checked: `tests/genome/test_qp_tori*.py`,
`tests/search/test_variational_qp_torus.py`, `test_qp_torus_fixed_jacobi_continuation.py`
-- none use a DRO-class or other genuinely stable family). Building one was
judged out of scope for this task (new periodic-orbit-family machinery, not a
benchmark). This note therefore answers a narrower, still-useful question:
**at the specific small-amplitude, near-Neimark-Sacker-bifurcation point this
project's own regression suite already treats as "the easy case," which
corrector wins?** -- not the literal Baresi stable-orbit case.

## Results (measured, this session, 2026-07-28, 8-core Mac)

```
=== GMOS (genome.qp_tori.correct_qp_torus) ===
wall time:                 95.17 s
converged:                 False   (residual 5.53e-5, above the requested tol=1e-8)
n_iter (nfev):              54
rotation number:            0.021400
invariance_residual (L2):   5.530e-05
independent_closure_resid:  6.369e-05

=== PDE (search.variational_qp_torus.discover_qp_torus), bootstrapped from GMOS above ===
wall time (PDE stage only): 121.11 s
wall time (GMOS+PDE total): 216.28 s
converged:                  False   (residual_rms 3.59e-6, above requested closure_tol)
n_iter (nfev):               139
rotation number:             0.021400
residual_rms (PDE invariance):  3.593e-06
closure_residual (independent): 8.308e-07
```

Rotation number agrees to the printed precision between the two independent
methods (`0.021400` both) -- the cross-check the existing regression test
itself relies on (`test_l2_positive_control_reproduces_gmos_torus` asserts
`res.rotation_number == pytest.approx(gmos_rot, rel=1e-3)`).

## Finding: Baresi/Olikara/Scheeres's own "GMOS is faster and more accurate"
does NOT clearly hold in this codebase's own implementations, at this test
point

1. **Speed:** GMOS alone (95.2 s) was not obviously faster than the PDE stage
   alone (121.1 s) -- both are the same order of magnitude, GMOS somewhat
   ahead but not decisively (their paper's own Table 1 shows GMOS beating
   PDE(DFT) by orders of magnitude on its linear-algebra step specifically).
   Since `discover_qp_torus` in this codebase is a **bootstrap FROM a
   converged (or near-converged) GMOS torus**, not an independent from-scratch
   method, the honest "cost of using the PDE corrector" here is the full
   `GMOS + PDE` pipeline (216.3 s, 2.27x GMOS-alone) -- there is currently no
   way in this codebase to run the pseudospectral PDE corrector without first
   paying for a GMOS solve. That is itself a finding: **for cases where GMOS
   already converges, running the PDE corrector afterward is pure added cost
   with no independent-method benefit** (see recommendation below).
2. **Accuracy:** the PDE corrector's own residual (3.59e-6) and especially its
   **independent closure check (8.31e-7)** were both roughly one to two orders
   of magnitude tighter than GMOS's own residual (5.53e-5) and independent
   closure (6.37e-5) at this specific point -- the opposite of Baresi et al.'s
   own accuracy finding for their PDE(DFT) vs. GMOS comparison.
3. **Neither corrector satisfied its own requested tolerance** at these
   parameters (`tol=1e-8` for GMOS, implicit `1e-5` closure gate for PDE) --
   both landed in a "good enough for the regression test's loose bounds, not
   fully Newton-converged" state, consistent with this exact test file's own
   documented platform-dependent fragility at this near-bifurcation point
   (`#635`'s eigenvector-phase-pinning note). This reduces confidence that the
   measured numbers are a clean apples-to-apples "both fully converged"
   comparison; they are better read as "both correctors' behavior at a shared,
   realistic, moderately-converged stopping point."

**Most likely explanation, checked directly in the source (not speculation):**
this is very plausibly an **implementation-maturity asymmetry, not an
intrinsic per-method result**. `genome/qp_tori.py::_correct_gmos` calls
`scipy.optimize.least_squares(..., method="trf", ...)` **without an analytic
`jac` argument**, so scipy defaults to 2-point **finite-difference**
Jacobians (`diff_step=1e-5`) -- exactly the expensive linear-algebra path
Baresi et al.'s own Table 1 shows GMOS *avoiding* in their implementation
(their GMOS: `0.05 s` of linear-algebra vs. PDE(DFT)'s `49.7 s`, at `N=25`).
By contrast, `search/variational_qp_torus.py` has a **hand-derived analytic
Jacobian** (`_jacobian`), verified against finite differences in
`test_analytic_jacobian_matches_finite_difference` to `<1e-6` relative error
-- the more mature, better-optimized implementation of the two correctors in
this codebase. The paper's comparison assumed both methods were similarly
well-implemented (or FD-based) on their own end; here they are not
symmetrically implemented, which plausibly inverts the paper's own speed
ranking for this codebase specifically.

## Floquet stability "for free" -- confirmed NOT available in this codebase's
GMOS implementation

`genome/qp_tori.py`'s own module docstring is explicit about this already
(re-confirmed by reading the implementation, not just the docstring): the
free-Floquet-stability byproduct Baresi et al. cite as a GMOS advantage is a
feature of the **Olikara 2016 thesis's Gauss-Legendre collocation** variant of
GMOS -- a *different, unimplemented* method. This codebase's actual GMOS
(`correct_qp_torus`/`_correct_gmos`) is the earlier Olikara & Scheeres
2010/2012 **shooting** scheme, solved via a black-box `scipy.optimize.
least_squares` call with no Jacobian introspection for stability information.
No Floquet multipliers of the invariant circle are extracted anywhere in
`qp_tori.py`. So: **Baresi et al.'s "free Floquet stability" advantage does
NOT transfer to this project's own GMOS implementation** -- it would require
implementing the 2016-thesis collocation variant (already flagged
high-priority in `#732`'s own citation-mining pass), not just using the
existing shooting-scheme code differently.

## Net verdict

Baresi/Olikara/Scheeres's own finding (GMOS preferred for stable cases) does
**not** clearly hold when checked directly in this project's own codebase, at
the one shared test point available (`#732`'s own recommended follow-up). Two
important caveats limit how far this generalizes:

* The tested parent orbit is not actually in the paper's own stable-orbit
  scope (spectral radius 1161, not `~1`) -- this project currently has no
  genuinely stable-parent-orbit torus case to test the paper's claim on its
  own terms.
* The result is very plausibly explained by an implementation asymmetry
  (GMOS here uses FD Jacobians via a generic scipy solver; PDE here uses a
  hand-verified analytic Jacobian) rather than a property of the GMOS/PDE
  method classes themselves.

**No change to how the PDE corrector is used for its actual job** (crossing
the unstable-parent-orbit wall GMOS cannot handle) is warranted or proposed --
per `#739`'s own scope, this is purely informational. The one concrete,
actionable takeaway: **for cases where GMOS already converges (this project's
own low-amplitude regime), do not bother also running the PDE corrector
afterward for "extra confidence" -- it costs 2.27x the GMOS-alone wall time
for a method that is currently only reachable by bootstrapping from GMOS's own
output, not an independently-run alternative.** Where accuracy actually
matters at these small amplitudes (independent closure 8e-7 vs. 6e-5 here),
that is available from the PDE stage specifically, at the stated 2.27x total
cost -- a genuine, real trade a future caller can make consciously, not a
reason to prefer one method universally.

## Why no permanent regression test was added

Considered adding a locked-in comparison test but declined: (1) wall-clock
timing assertions are exactly what this same test file's own `#631`/`#632`/`#737`
history shows to be CI-flaky (2-core CI runners vs. 8-core local, `>600s`
timeouts already documented for the existing `slow`-marked tests in this same
file); (2) neither corrector fully converged to its own requested tolerance
at these parameters, so the residual/closure numbers are themselves
platform- and BLAS-dependent (per this file's own `#635` note on
non-bit-reproducible DOP853 stepping) -- pinning them tightly would be
fragile, and pinning them loosely would add ~200s of CI time for a test that
duplicates machinery already exercised by
`test_l2_positive_control_reproduces_gmos_torus`. The numbers above are
reproducible by rerunning that existing test's own setup with timing added
(exact recipe: this note's benchmark script, not committed, available in this
session's scratchpad if needed again).
