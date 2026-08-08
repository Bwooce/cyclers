# #803 — rails-cache batch-parity test failure: final-knot ULP knife edge (test-design bug)

**Verdict: test-design bug, fixed with an honestly-derived tolerance. No production code
changed; no xfail; no past-search invalidation.**

## Symptom

`tests/nbody/test_propagator_api.py::test_rails_cache_batch_samples_match_per_point`
failed reproducibly in isolation on this Mac (now also the self-hosted CI runner):

```
assert np.float64(2.9802322387695312e-08) <= 1e-09
```

## Root cause

The residual `2.9802322387695312e-08` is exactly `2**-25` — **1 ULP of a double in
`[2^27, 2^28) = [1.342e8, 2.684e8) km`**, i.e. one last-bit of a planetary heliocentric
coordinate. (The dispatch brief's float32-round-trip lead was a red herring: `2**-25` is
the float32 *relative* epsilon scale, but here it appears as an *absolute* 1-ULP error on
a ~1.69e8 km double.) The single failing sample was Mars' x-coordinate
(`169056503.30759174` vs `...17`) at grid index **70 of 70 — the FINAL spline knot**.

Two hypotheses were tested directly:

1. **Batch ephemeris parity (the thing the test says it pins): INTACT.**
   `Ephemeris("astropy").states()` was compared element-by-element against per-point
   `ingest_planet_state()` for all 3 bodies x 71 knots, both `r` and `v`:
   **bit-identical throughout** (`np.array_equal`). The `#692` fix (elementwise
   ICRS→ecliptic rotation instead of BLAS dgemv/dgemm matrix products) holds.

2. **CubicSpline knot reproduction: exact at interior knots, NOT at the final knot.**
   scipy's PPoly interval lookup is left-closed: at any interior knot `t_i` the evaluated
   interval starts at `t_i`, offset 0, and the returned value is the stored constant
   coefficient `c[3,i]` — which `CubicSpline` sets to the sample exactly. Measured
   interior diffs: `0.0` for all bodies. But at the **last** knot, PPoly evaluates the
   *last interval's* cubic at offset `h = 86400 s`. That is a genuine 4-term
   floating-point polynomial evaluation, and it reproduces the endpoint sample only to
   rounding. The exact scipy kernel (`_ppoly.pyx` ascending-power accumulation,
   `res += c[3-k]*z; z *= s` — *not* Horner) was reproduced op-for-op in pure Python and
   yields exactly the observed 1-ULP-low value; a Horner evaluation of the same
   coefficients happens to round exact, confirming the outcome is pure
   evaluation-order/rounding luck.

Whether the final-knot evaluation rounds bit-exact depends on the spline coefficients,
which come from scipy's banded LAPACK solve — **Accelerate** on this machine
(`scipy 1.17.1 show_config`: blas/lapack = Accelerate, aarch64). Venus and Earth happen
to round exact; Mars lands 1 ULP off. The old fixed `1e-9 km` tolerance sits *below* one
ULP of any coordinate in `[1.34e8, 2.68e8) km`, so at the final knot the test was
demanding exact endpoint rounding — a property cubic-spline evaluation never guaranteed.
It passed historically by coefficient-level rounding luck; a scipy/Accelerate build or OS
BLAS drift plausibly flipped a coefficient ULP and broke the knife edge. Same-machine CI
makes "cross-platform divergence" unavailable as an explanation, and none is needed: the
mechanism is fully accounted for without it.

## Fix

`tests/nbody/test_propagator_api.py` (test only — `src/` untouched):

- The I5 batching pin is now asserted **directly and more strongly**: batched
  `ephem.states()` samples must be `np.array_equal` (bitwise) to per-point
  `ingest_planet_state()` for both `r` and `v` at every knot.
- The spline-evaluation check uses an honestly-derived bound:
  `4.0 * np.spacing(np.max(np.abs(ref)))` — 4 ULPs of the largest coordinate
  (~1.2e-7 km at Mars). Headroom of 4 ULPs covers evaluation-order/LAPACK-build
  variation in the final-knot polynomial evaluation; interior knots remain exact in
  practice (measured 0.0) and pass trivially.

Context for the bound: the cache's *actual* accuracy claim is mid-interval interpolation
against DE440, measured at ~0.08 km max (Mars, 1-day grid) — six orders of magnitude
above the new knot bound. Nothing physically meaningful is lost.

## Blast radius / past results

No production code changed — `RailsEphemerisCache` behaves exactly as it always did, and
its consumers (`nbody/shooter.py`, `maintenance_shoot.py`, `jovian.py`, `propagator.py`,
`scripts/reproduce_mcconaghy_table71.py`) consume it inside force evaluations where the
cache's own ~0.08 km interpolation error dominates a final-pad-knot ULP by ~6 orders of
magnitude. No past discovery/search result is affected; no follow-up re-runs needed.

`#804` (`test_da_section_map` 2.8x margin miss) is **not** the same class: that failure
is 2.8x past a 1e-4 convergence bound, far above ULP scale — it needs its own
investigation and stays open.

## Verification

- `uv run pytest tests/nbody -q` — green (exit 0).
- `uv run ruff check .` / `uv run ruff format --check .` — clean.
- `uv run mypy src tests` — clean (839 files).
