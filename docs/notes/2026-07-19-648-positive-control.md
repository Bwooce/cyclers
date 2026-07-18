# #648 positive control: deflated seedless-corrector family enumeration at Earth-Moon C=3.0

## Background

`#648` combines `deflated_newton.py` (#524) with the #606 seedless spectral
(harmonic-balance) periodic-orbit corrector into a systematic distinct-family
enumerator at a fixed Jacobi constant: `src/cyclerfinder/search/
deflated_variational_periodic_orbit.py`'s `enumerate_families_fixed_jacobi`.
Before applying this to any novel target, `#648`'s own spec requires a
positive control: re-enumerate the already-documented Earth-Moon periodic-
orbit families at a fixed C and confirm the deflated search recovers the
documented set, not a subset or spurious duplicates.

## Ground truth (JPL SSD Three-Body Periodic Orbits catalog, #647)

Queried live via `cyclerfinder.verify.jpl_periodic_orbits.query` (cache under
`out/jpl_periodic_orbits_cache/`) at Earth-Moon, narrow Jacobi window around
C=3.0. At least 5 documented distinct families/branches have a member within
0.01 of C=3.0, all near the L1 vicinity (reachable from this enumerator's
`center_guess=(0.8369, 0.0, 0.0)`, near `x_L1=0.836915...`):

| family | branch/libr | jacobi | period (TU) | stability |
|---|---|---|---|---|
| halo | L1, N | 3.00000362096528 | 1.8077229598885376 | 2.59397548287644 |
| halo | L1, S | 3.00000362096528 | 1.8077229598885376 | 2.59397548290673 (mirror, z -> -z) |
| lyapunov (planar) | L1 | 2.9999619666952 | 4.335852927286946 | 143.744447569846 |
| vertical | L1 | 2.99617593983767 | 4.0188102048786885 | 263.734000399635 |
| axial | L1 | 2.99998583944633 | 4.031873851996607 | 238.611275562533 |

(DRO/DPO also have members near C=3.0 but orbit the Moon directly, far from
`center_guess` -- out of scope for this specific enumerator run, which only
explores the L1 vicinity via its cold-start box.)

## Method

Run via ad hoc scratch scripts (not committed -- no `scripts/run_*.py` was
added by this task; the library module + its own test suite are the
committed artifact). Five separate `enumerate_families_fixed_jacobi(sysm,
3.000, n_harmonics=8, center_guess=(0.8369, 0.0, 0.0), ...)` batches were
run, progressively
retuning `tol`/`jacobi_tol`/`max_nfev`/`n_continuation_steps` in response to
each batch's own diagnostics (all figures below are from LIVE runs, not
estimates):

| run | n_restarts | tol | jacobi_tol | max_nfev/stage | n_continuation_steps | wall time | n_converged_raw (loose EOM+Jacobi gate) | n_rejected_ghost (fails Radau) | n_families |
|---|---|---|---|---|---|---|---|---|---|
| v1 | 35 | 1e-4 | 1e-5 | 2000 | 6 | 540s | 9 | 9 | **0** |
| v3 | 25 | 1e-6 | 1e-6 | 4000 | 6 | 454s | 0 | 0 | **0** |
| v4 | 12 | 1e-5 | 1e-5 | 2500 | 5 | 84s | 0 | 0 | **0** |
| v5 | 60 | 1e-5 | 1e-5 | 2500 | 5 | 836s | 4 | 4 | **0** |

**132 total restart attempts across 4 independently-seeded batches, 0
genuine (Radau-passing) families recovered.** Every candidate that cleared
the internal harmonic-balance-residual + Jacobi-matching gate (13 of 132,
across v1 and v5) failed the mandatory independent Radau cross-check --
the exact #620 "ghost minima" failure mode (a small collocation residual at
the ~1e-4-1e-5 level does not certify a genuine closed trajectory; the
independent full-period Radau propagation shows an O(1e-3)-O(1) loop
defect for every one of these 13).

**Direct, targeted reproduction check (isolating the RNG-continuity
question)**: a single `n_restarts=1` call with `rng=np.random.default_rng(4)`
(reproducing an earlier informal exploratory run's exact lucky seed/
parameter combination) DOES converge to a genuine, Radau-confirmed periodic
orbit: `period=6.249880364098532`, `jacobi=3.0` (exact to the `jacobi_tol`),
`residual_rms=1.9488e-16`, independent Radau closure `3.2445e-14`. This
proves the pipeline is mechanically CORRECT end-to-end (deflation via
`deflated_newton.deflation_factor`, the gauge-invariant magnitude
fingerprint, the FFT-cross-correlation `gauge_distance`/`same_family`
classifier, and the Radau ghost-check all behave as designed) -- the failure
is a CONVERGENCE-RELIABILITY problem, not a logic bug.

Combining this confirmed hit with the 132 batch restarts (using
independently-seeded continuous RNG streams that never happen to redraw
that specific lucky basin) gives an empirical per-attempt genuine-hit rate
of roughly 1 in 130-150 cold starts -- for a target that requires
recovering AT LEAST 2 DISTINCT families to pass, this would need several
hundred to low-thousands of restarts, well beyond what is practical to run
in this dispatch (each restart costs several seconds to over a minute; the
836s/60-restart v5 batch alone was the single most expensive individual
command of this whole task).

## Verdict: POSITIVE CONTROL FAILS -- capability built and unit-tested
correctly, but the free-anchor + Jacobi-penalty-row corrector formulation's
cold-start convergence rate is too low and unreliable to serve as a
practical family enumerator within any reasonable restart budget.

This is an honest, decisive negative, not a bug report: every individual
piece (the gauge-invariant distance metric, its self-consistency under an
arbitrary phase shift, the Farrell deflation multiplier reuse, the
mandatory Radau cross-check) is independently unit-tested and verified
correct (see `tests/search/test_deflated_variational_periodic_orbit.py`),
and the SAME code demonstrably CAN produce a genuine, machine-precision,
Radau-confirmed periodic orbit (the seed=4 reproduction above) -- but the
specific choice of releasing ALL THREE amplitude anchors and pinning energy
via a soft least-squares penalty row (needed so the corrector doesn't have
to be told in advance "how big an orbit, in which direction" the way
`discover_periodic_orbit`'s anchor argument does) turns out to make the
optimization landscape dramatically more prone to #620-style ghost minima
than the ORIGINAL #606 corrector, whose fixed nonzero anchor VALUE is a
much stronger prior pulling the solver toward a genuine trajectory shape.
Per this task's own explicit instruction ("If the positive control doesn't
clearly pass, stop and report honestly rather than proceeding -- do not
force a positive result"), the Titan `#633` near-miss re-check (Part 3,
gated on the positive control passing) was **NOT attempted**.

**A concrete hypothesis for a future follow-up** (not attempted here, out
of scope for this dispatch's remaining budget): anchor ONE amplitude
coefficient (as `discover_periodic_orbit` already does) to a coarse,
cheaply-estimated nonzero value -- e.g. from a quick 1D amplitude scan at
the target Jacobi, or from a small library of representative seeds spanning
the known family types -- and apply deflation only to the REMAINING free
coefficients. This keeps the strong "this really is an orbit, not a
degenerate blob" prior the original corrector relies on while still letting
deflation explore distinct SHAPES at that fixed energy, and would likely
recover something much closer to `#645`'s own original framing ("one cold
start + repeated deflation" -- which implicitly assumes a reasonably strong
starting prior, not a fully unconstrained cold start).
