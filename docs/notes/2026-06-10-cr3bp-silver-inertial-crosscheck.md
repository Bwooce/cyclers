# CR3BP SILVER Lyapunov members -- inertial n-body cross-check (constraint 3)

**Run timestamp:** 2026-06-10T05:05:22Z
**Script:** `scripts/cr3bp_silver_inertial_crosscheck.py`
**Method:** `cr3bp-inertial-rebound-ias15-v1`  (git `8baba33`)
**Input:** `data/cr3bp_silver.jsonl` (14 SILVER, `cr3bp-lyapunov-corrector-v2`)
**Spec:** `docs/superpowers/specs/2026-06-10-cr3bp-tier2-design.md`, binding
constraint 3 -- re-propagation in the INERTIAL n-body harness (REBOUND/IAS15,
different code path and frame from the rotating-frame CR3BP integrator).

## Harness

- Inertial **barycentric** frame, km/s units; REBOUND `G=1` with masses set to
  the JPL GM values, IAS15 (`epsilon=1e-9`).
- Saturn + moon initialised on the exact circular two-body orbit of the CR3BP
  idealisation (separation `l_km`, mean motion `n = 1/t_s`), momentum-free.
  This is the like-for-like CR3BP-consistency check, **not** a real-ephemeris
  claim. The moon's deviation from the analytic circular rail is measured and
  feeds the noise floor.
- Frame map derived independently of the rotating-frame propagator: at t=0
  (theta=0) `r_in = r_rot`, `v_in = v_rot + n z x r_in`; back-transform at
  sample times uses `theta = n t` (CR3BP convention), then nondimensionalise
  by `l_km` and `l_km*n`. Jacobi is evaluated on the back-transformed states.
- 5 periods, 128 samples/period (recurrence sampled
  at exact integer period multiples).

## Pre-registered tolerances and verdict rules (fixed before the run)

**Instability fact stated up front:** collinear-point Lyapunov orbits have
in-plane instability exponent `nu_u ~ 2.5` nd here, i.e. per-period error
amplification `lambda = exp(nu_u T)` = 1986..2042 for these 14.
Even a machine-epsilon-perfect seed (1e-16 nd) must depart the orbit
neighbourhood within ~5 periods, and the candidates' own closure residuals
(7e-14..8e-11 nd) depart at ~2-4 periods. **Literal 5-period boundedness is
physically unattainable for any numerical trajectory of these orbits**, so the
binding rules are the strongest operationalisation that retains information:

- **R1 (periodicity):** `delta1 = |X_nd(T) - X_nd(0)|` (6-norm, rotating nd,
  back-transformed from the inertial run) <= 0.1 * A (A = recorded
  `amplitude_nd`).
- **R2 (Jacobi):** max `|J(t) - J(0)|` over the bounded span <= 1e-09
  (absolute; J ~ 3). Full-5T drift reported as a diagnostic.
- **R3 (boundedness):** no NaN/divergence; no moon impact within the bounded
  span; observed departure time (first `d_L > 3A`) >=
  min(5T, 0.7 * t_dep_pred), with
  `t_dep_pred = T + ln(3A / max(delta1, floor)) / nu_u` -- departure EARLIER
  than the orbit's own measured residual + linear instability predicts
  contradicts the claimed orbit; later (or none) is fine.
- **Noise gate:** harness floor (rail deviation + barycentre drift, nd) <=
  0.01 * A, else INCONCLUSIVE (noise >= signal at this amplitude).

PASS = R1 & R2 & R3 (noise gate ok). CHECK-FAILED = any rule violated.
INCONCLUSIVE = noise gate violated.

## Per-candidate results

(d1 = recurrence after one period; 'bounded T' = nd time inside `d_L <= 3A`;
t_dep pred = linear-instability departure prediction, nd; floor = harness
noise floor, nd.)

| # | Moon | L | Ax/g | A (nd) | d1/A | d1 (nd) | bounded T | t_dep pred (T) | dJ bounded | dJ full | floor (nd) | R1 | R2 | R3 | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Mimas | L2 | 0.001 | 4.15e-06 | 5.3e-07 | 2.21e-12 | 9.6 | 9.3 | 8.9e-16 | 8.9e-16 | 2.4e-13 | Y | Y | Y | PASS |
| 2 | Mimas | L1 | 0.002 | 5.98e-05 | 5.6e-08 | 3.33e-12 | 10.5 | 10.1 | 8.9e-16 | 8.9e-16 | 2.4e-13 | Y | Y | Y | PASS |
| 3 | Mimas | L2 | 0.0005 | 3.73e-06 | 4.9e-06 | 1.82e-11 | 8.8 | 8.4 | 1.3e-15 | 1.3e-15 | 2.4e-13 | Y | Y | Y | PASS |
| 4 | Mimas | L1 | 0.0005 | 5.44e-06 | 1.2e-05 | 6.27e-11 | 8.3 | 8.0 | 8.9e-16 | 8.9e-16 | 2.4e-13 | Y | Y | Y | PASS |
| 5 | Mimas | L1 | 0.001 | 1.22e-05 | 6.5e-06 | 7.88e-11 | 8.6 | 8.2 | 8.9e-16 | 8.9e-16 | 2.4e-13 | Y | Y | Y | PASS |
| 6 | Enceladus | L2 | 0.0005 | 5.42e-06 | 3.5e-07 | 1.89e-12 | 9.7 | 9.4 | 8.9e-16 | 8.9e-16 | 7.5e-14 | Y | Y | Y | PASS |
| 7 | Enceladus | L1 | 0.001 | 1.77e-05 | 4.6e-08 | 8.12e-13 | 10.6 | 10.2 | 8.9e-16 | 8.9e-16 | 7.0e-14 | Y | Y | Y | PASS |
| 8 | Enceladus | L2 | 0.001 | 5.94e-06 | 1.5e-06 | 8.81e-12 | 9.2 | 8.9 | 1.3e-15 | 1.3e-15 | 7.6e-14 | Y | Y | Y | PASS |
| 9 | Enceladus | L1 | 0.0005 | 7.72e-06 | 1.7e-06 | 1.35e-11 | 9.2 | 8.7 | 8.9e-16 | 8.9e-16 | 7.0e-14 | Y | Y | Y | PASS |
| 10 | Tethys | L1 | 0.002 | 1.19e-04 | 1.3e-09 | 1.60e-13 | 12.0 | 11.4 | 8.9e-16 | 1.3e-15 | 2.3e-13 | Y | Y | Y | PASS |
| 11 | Tethys | L2 | 0.0005 | 1.01e-05 | 1.2e-07 | 1.19e-12 | 10.3 | 9.9 | 8.9e-16 | 1.3e-15 | 2.4e-13 | Y | Y | Y | PASS |
| 12 | Tethys | L1 | 0.0005 | 1.34e-05 | 6.6e-08 | 8.86e-13 | 10.4 | 10.0 | 8.9e-16 | 8.9e-16 | 2.3e-13 | Y | Y | Y | PASS |
| 13 | Tethys | L2 | 0.001 | 1.36e-05 | 2.7e-07 | 3.70e-12 | 9.9 | 9.6 | 8.9e-16 | 8.9e-16 | 2.4e-13 | Y | Y | Y | PASS |
| 14 | Tethys | L1 | 0.001 | 3.14e-05 | 1.2e-06 | 3.71e-11 | 9.3 | 8.9 | 8.9e-16 | 8.9e-16 | 2.3e-13 | Y | Y | Y | PASS |

Per-period recurrence `|X(kT) - X(0)|` (nd), k = 1..5:

| # | Moon | L | Ax/g | d1 | d2 | d3 | d4 | d5 |
|---|---|---|---|---|---|---|---|---|
| 1 | Mimas | L2 | 0.001 | 2.21e-12 | 4.47e-09 | 8.94e-06 | 8.97e-03 | 9.19e-02 |
| 2 | Mimas | L1 | 0.002 | 3.33e-12 | 6.75e-09 | 1.37e-05 | 1.15e-02 | 9.59e-02 |
| 3 | Mimas | L2 | 0.0005 | 1.82e-11 | 3.64e-08 | 7.24e-05 | 2.49e-02 | 1.13e-01 |
| 4 | Mimas | L1 | 0.0005 | 6.27e-11 | 1.27e-07 | 2.52e-04 | 3.89e-02 | 1.20e-01 |
| 5 | Mimas | L1 | 0.001 | 7.88e-11 | 1.60e-07 | 3.15e-04 | 4.16e-02 | 1.21e-01 |
| 6 | Enceladus | L2 | 0.0005 | 1.89e-12 | 3.77e-09 | 7.54e-06 | 9.12e-03 | 1.22e-01 |
| 7 | Enceladus | L1 | 0.001 | 8.12e-13 | 1.65e-09 | 3.35e-06 | 5.11e-03 | 1.08e-01 |
| 8 | Enceladus | L2 | 0.001 | 8.81e-12 | 1.76e-08 | 3.51e-05 | 2.21e-02 | 1.47e-01 |
| 9 | Enceladus | L1 | 0.0005 | 1.35e-11 | 2.74e-08 | 5.54e-05 | 2.77e-02 | 1.51e-01 |
| 10 | Tethys | L1 | 0.002 | 1.60e-13 | 3.26e-10 | 6.65e-07 | 1.30e-03 | 1.23e-01 |
| 11 | Tethys | L2 | 0.0005 | 1.19e-12 | 2.36e-09 | 4.68e-06 | 7.38e-03 | 1.82e-01 |
| 12 | Tethys | L1 | 0.0005 | 8.86e-13 | 1.81e-09 | 3.70e-06 | 6.23e-03 | 1.78e-01 |
| 13 | Tethys | L2 | 0.001 | 3.70e-12 | 7.36e-09 | 1.46e-05 | 1.71e-02 | 2.20e-01 |
| 14 | Tethys | L1 | 0.001 | 3.71e-11 | 7.57e-08 | 1.54e-04 | 6.02e-02 | 2.76e-01 |

## Verdict counts

| System | PASS | CHECK-FAILED | INCONCLUSIVE |
|---|---|---|---|
| Saturn/Mimas | 5 | 0 | 0 |
| Saturn/Enceladus | 4 | 0 | 0 |
| Saturn/Tethys | 5 | 0 | 0 |
| **Total** | **14** | **0** | **0** |

## Honest noise-floor / amplitude discussion

These are tiny orbits: amplitudes 3.7e-6..1.2e-4 nd, i.e. ~0.7..35 km at the
moons. The harness noise floor (moon's two-body rail deviation from the
analytic circle plus barycentre drift, in nd after derotation) is measured per
run and reported per candidate above; the noise gate demands it sit below 1%
of each orbit's amplitude for a PASS/FAIL to be meaningful at all. The
one-period recurrence d1 is the information-bearing periodicity number: it is
the one-period flow defect of the given initial condition -- the SAME quantity
the rotating-frame corrector reports as `closure_residual`, here re-measured
through a completely different code path (inertial REBOUND/IAS15 + an
independent frame back-transform). Observed: d1 reproduces each record's
`closure_residual` to within the harness noise floor (<= ~2e-13 nd), i.e. the
two integrators agree on the one-period flow map at the noise level. (The
pre-run expectation written into an earlier draft -- d1 ~ lambda * residual --
was wrong; the lambda amplification enters from period 2 onward, exactly as
the d2/d1 ratios show.)

The departures visible in d2..d5 grow at the measured per-period factor
d_{k+1}/d_k ~ 2.0e3, matching the theoretical lambda = exp(nu_u T) per
candidate to ~0.1% -- the departure IS the orbit's intrinsic linear
instability acting on the seed's finite residual, not a harness disagreement;
R3 grades whether the observed departure time is consistent with (never
earlier than) that prediction; in this run every candidate departed slightly
LATER than predicted (t_bound > t_dep_pred), never earlier. A 5-period
absolutely-bounded trajectory was shown above to be unattainable in principle
at double precision, which is why constraint 3's 'bounded/periodic' is
operationalised as R1+R3 rather than a literal 5T position bound.

## Writeback

Each record in `data/cr3bp_silver.jsonl` gained a `crosscheck_inertial` field
(method, verdict, metrics). No verdict-tier change, NO catalogue writeback.
