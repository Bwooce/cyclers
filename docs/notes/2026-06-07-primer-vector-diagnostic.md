# Primer-vector optimality diagnostic — results note (2026-06-07)

**Status: DIAGNOSTIC / PROVISIONAL.** Task #144. Read-only first-order
necessary-conditions check on existing impulsive schedules. Implemented in
`src/cyclerfinder/verify/primer.py`; tests in `tests/verify/test_primer.py`.
See the external-algorithms survey Thread 3
(`docs/notes/2026-06-07-external-algorithms-survey.md`) for background.

## What it computes

Classical impulsive primer-vector theory (Lawden 1963; Lion & Handelsman
1968, AIAA J. 6(1):127-132, DOI 10.2514/3.4452; textbook treatment Prussing &
Conway / Prussing 2010 Ch. 2). On a ballistic (coast) arc the primer obeys the
same linear variational equation as a position perturbation in the two-body
field:

    p̈ = G(r) p ,    G(r) = (μ/r³)(3 r̂ r̂ᵀ − I)

`G` is symmetric and traceless by construction. Because `(p, ṗ)` share the
variational ODE with the state perturbation, the primer propagates with the
arc's 6×6 state-transition matrix (STM). The STM is integrated (DOP853, the
reference Kepler arc carried alongside and cross-checked against the project
`core.kepler.propagate` propagator). At each impulse the boundary condition is
`p(t_i) = Δv_i/|Δv_i|`, so `|p(t_i)| = 1` by construction; the interior
coast is a two-point BVP whose unknown `ṗ(t0)` is recovered by inverting the
STM's upper-right 3×3 block `Φ_rv`.

**Verdict per coast:** `max|p|`, time of max, and an enum —
`OPTIMAL_NECESSARY_CONDITIONS_MET` (`|p| ≤ 1 + tol` throughout) or
`IMPROVABLE_ADD_IMPULSE` (`max|p| > 1`, the Lion & Handelsman signature that an
added/relocated midcourse impulse lowers total ΔV). This is a **necessary**
conditions test only: it can refute optimality but never prove it (sufficiency
is not checked).

## Sourced golden gates and one correction to the brief

The task brief proposed two Hohmann golden gates: ratio 2 → OPTIMAL, and
**ratio 20 → interior `max|p| > 1` (IMPROVABLE), citing the ~11.94 Lawden /
Marchal bi-elliptic threshold.** The second gate is physically incorrect as
stated, and the diagnostic confirms it empirically.

**Finding (verified, grid-converged, n_samples up to 2000):** the symmetric
180° two-impulse Hohmann transfer coast has `max|p| = 1.000000` reached **at
the endpoints** for *every* radius ratio tested (2, 11.94, 15.58, 20, 50); the
interior dips strictly below unity. There is **no** interior `|p| > 1` bulge on
the Hohmann transfer arc at any ratio.

This does not contradict the published ~11.94 bi-elliptic threshold. That
threshold is governed by a *different* Lawden necessary condition — the
**endpoint coast-extension condition** (the sign of `d|p|/dt` evaluated on an
arc *extended beyond* the burn, i.e. condition #4 in the survey), not by an
interior bulge on the fixed two-impulse transfer arc. Capturing it requires
propagating a coast extension past the impulse and checking the slope there —
a separate, larger feature (Lawden's add-an-initial/final-coast test). It is
**out of scope** for this read-only interior diagnostic and was deliberately
not implemented; inventing an interior-violation assertion at ratio 20 would
have violated the golden-test discipline (the EXPECTED side must trace to a
source, and no source claims an interior bulge there).

**Gates actually shipped (all sourced or pure construction invariants):**

- `G(r)` symmetry, tracelessness, and the radial/transverse eigenvalues
  `+2μ/r³` / `−μ/r³` (algebraic identities).
- `|p(t_i)| = 1` at both bounding impulses (Lawden BC, by construction).
- Continuity / grid-refinement convergence of `max|p|`.
- The STM reference arc matches `core.kepler.propagate`.
- **Hohmann ratio 2 and ratio 20 → OPTIMAL** (`|p| ≤ 1` on the coast,
  endpoints the maxima) — the *honest* sourced statement for the symmetric
  Hohmann transfer (Prussing & Conway).
- **A fixed-time "long-way" (transfer angle > 180°) Lambert transfer →
  IMPROVABLE** with an interior peak — the classic Lion & Handelsman
  non-optimal case. Only the qualitative `max|p| > 1` necessary-condition
  violation is asserted (no invented magnitude).

## Application — Aldrin E-M-E maintenance schedule

Run via `optimise_aldrin_maintenance_dv` (real ephemeris, `astropy`,
real-window priority 1985-10-28; the V∞-anchor seeding used in
`scripts/maintenance_batch.py`). The impulse schedule is built from the
cycler's heliocentric legs (each leg = a coast; the heliocentric ΔV at each
encounter = the impulse direction).

**Result (grid-converged; gated `@pytest.mark.slow`, runs in ~21 s):**

```
maintenance_dv_kms (our value) = 2.9138   (matches the published in-family value)
coast 0 (E->M, 131.9 d): max|p| = 1.12228 at t/T = 0.355  -> IMPROVABLE_ADD_IMPULSE
coast 1 (M->E, 598.7 d): max|p| = 1.00008 at t ~ endpoint -> IMPROVABLE_ADD_IMPULSE
overall: IMPROVABLE_ADD_IMPULSE
```

**Interpretation (DIAGNOSTIC, PROVISIONAL):**

- **Coast 0 (Earth→Mars) carries a genuine interior bulge**, `max|p| ≈ 1.122`
  peaked ~35% along the 132-day arc and stable under grid refinement
  (300 → 4000 samples, identical to 5 sig figs). By the Lion & Handelsman
  diagnostic this is first-order evidence that a midcourse impulse near
  t ≈ 47 d would reduce the total ΔV — i.e. the current schedule's maneuvers
  are **not** optimally placed/timed on this leg. The recoverable ΔV is *not*
  quantified here (that needs the add-impulse predictor / corrector, which is
  out of scope), so the implication is qualitative only: a fraction of the
  2.9138 km/s is plausibly recoverable, magnitude TBD.
- **Coast 1 (Mars→Earth) is at the noise floor**, `max|p| = 1.00008` with the
  peak at the endpoint (t ≈ 2 d into a 599-day arc). This is a *marginal* touch
  of unity, not a real interior violation; treat coast 1 as effectively
  satisfying the necessary conditions. It tips the overall verdict to
  IMPROVABLE only because of the `> 1.0 + tol` threshold; the substantive
  signal is entirely coast 0.

**Caveat (binding).** The maintenance ΔV (2.9138 km/s) is our own computed
number, not a sourced value — McConaghy 2002 defers it. The long M→E leg
(599 days, multi-revolution-scale) is exactly the regime where the Guzman et
al. 2002 primer survey (IAC-02-A.6.09 / NTRS 20030032208, **not yet
acquired**) flags linearised primer theory as fragile (singularities along
arcs). Until that survey is in hand and the STM propagation is validated
against its cases, **all Aldrin primer results are provisional** and must not
be cited as an optimality claim. The Hohmann and long-way goldens (short,
single-rev, analytic) are the trustworthy validation of the machinery itself.

## What is NOT done (deliberate scope limits)

- No re-optimisation / add-impulse corrector (Jezewski & Rozendaal). The module
  is pure diagnostic.
- No endpoint coast-extension condition (the one that actually detects the
  bi-elliptic threshold). Documented above as a possible follow-up.
- No multi-rev STM validation against Guzman 2002 (PDF not acquired).
