# #567 step (3) diagnosis: post-fix isolated singleton flips are genuine synodic aliasing, not a third artifact

Date: 2026-07-11
Author: diagnostic pass on `data/scan_567_epoch_robustness.jsonl` (scan commit `a777c2d`,
V4-strict fixes commit `6c54bba`)

## Question

After #567(1)+(2) fixed the two diagnosed V4-strict artifacts (Lambert branch-selection
discontinuity; silent planet-crossing misclassification), the epoch-robustness scan still shows
`detect_isolated_singleton_anomalies()` flagging isolated PASS/FAIL flips at wildly
candidate-dependent rates (0% Titania-Oberon daily up to 44% Ariel-Umbriel annual). Per
[[feedback_isolated_sweep_flips_suspect_artifact]], a flip population must be *diagnosed* before it
is believed or discarded. Two live hypotheses:

- **(a)** a third, undiagnosed artifact generator in the same code family;
- **(b)** genuine physical aliasing — the real PASS/FAIL boundary oscillates fast enough that even a
  DAILY grid under-resolves it.

## Verdict: (b), decisively — with a small, separately-explained secondary population

The isolated singletons are **genuine physical aliasing of a real, fast-oscillating
planet-crossing feasibility boundary whose frequency is the moon-pair synodic frequency.** They are
an artifact of the *diagnostic sampling grid*, not of the solver. No third code artifact is
present. Four independent lines of evidence:

### 1. The pass/fail bit is a step function of a continuous geometric quantity

`passes_v4_strict` (v4_uranus_strict.py:874-878) fails a cycle two ways: (i) the aggregate
V4-vs-V3 drift agreement exceeds the 50,000 km floor, or (ii) `planet_crossing_infeasible` —
`_select_leg_transfer` (lines 494-537) finds that every Lambert branch for a leg has osculating
periapsis inside Uranus's `r_eq` (25,559 km). Across the whole scan, essentially all FAILs are
mode (ii): the transfer's periapsis dips below Uranus. Periapsis is a smooth function of the two
moons' positions, hence of launch epoch. The boolean flips exactly when `perijove(epoch)` crosses
`r_eq` — a smooth-crossing step, not a discontinuity.

### 2. FFT of each daily PASS/FAIL series peaks at that pair's synodic period

Dominant spectral period of the daily boolean series vs. the transfer pair's synodic period
(sidereal periods: Ariel 2.520 d, Umbriel 4.144 d, Titania 8.706 d, Oberon 13.463 d):

| pair (first two seq bodies) | synodic P | daily-series dominant FFT P | daily singletons (of ~366) |
|---|---|---|---|
| Titania-Oberon | 24.64 d | 24.4 d | **0** |
| Umbriel-Titania | 7.91 d | 7.96 d | 12-15 |
| Umbriel-Oberon (#312) | 5.99 d | 6.00 d | 35-40 |
| Ariel-Titania | 3.55 d | 3.55 d | 82-83 |
| Ariel-Oberon | 3.10 d | 3.10 d | 105-106 |
| Ariel-Umbriel | 6.43 d | **3.21 d (= P/2)** | 117-120 |

The FFT period matches the synodic period to <2% for five of six pairs. Ariel-Umbriel matches
**half** the synodic period: with `rel_offset_deg = 0` and a symmetric there-and-back
(Ariel-Umbriel-Ariel), the transfer periapsis dips below `r_eq` on *both* the ascending and
descending synodic passes, so the boundary oscillates at 2× synodic frequency. That shortest
effective period (3.2 d) is why Ariel-Umbriel has the highest flip rate.

Singleton count is **monotonic in boundary period**: 24.6 d → 0 flips (band is well-resolved and
FAILs form contiguous runs), down to 3.2 d → 120 flips (Nyquist-violated at daily sampling). This
monotone frequency-dependence is the textbook signature of aliasing; a numerical solver artifact
would not spectrally lock to the astronomical synodic frequency of each specific moon pair.

### 3. Zooming into finer resolution resolves a daily singleton into a contiguous FAIL band

Direct sub-daily probe (4-hour steps) of Ariel-Umbriel across the daily singleton at 2000-01-04
(FAIL flanked by PASS on 01-03 and 01-05 at daily resolution):

```
2000-01-03T00:00  PASS  converged                 v3drift 993
2000-01-03T12:00  PASS  converged                 v3drift 2325
2000-01-03T16:00  FAIL  planet_crossing  perij 17197
2000-01-03T20:00  FAIL  planet_crossing  perij    19
2000-01-04T00:00  FAIL  planet_crossing  perij   912
2000-01-04T08:00  FAIL  planet_crossing  perij  6677
2000-01-04T16:00  FAIL  planet_crossing  perij 17497
2000-01-04T20:00  FAIL  planet_crossing  perij 24619   (approaching r_eq 25559)
2000-01-05T00:00  PASS  converged                 v3drift 4701
2000-01-06T00:00  PASS  converged                 v3drift 3827
```

The "isolated" daily FAIL is really a **contiguous ~1.3-day-wide feasibility gap** (01-03T16:00 →
01-04T20:00). The offending branch's periapsis varies **smoothly and continuously** across it
(down to 19 km near center, back up to ~24,619 km — just below `r_eq` — at the exit edge), and the
V4-vs-V3 drift on both PASS sides is small and smooth (~1,000-4,700 km). The daily grid simply
placed one sample inside a real gap and both neighbors outside it. This is the "zoom into finer
resolution" test prescribed by the feedback rule, and it converts the singleton into a resolved
band — the definitive discriminator: a solver artifact would stay isolated/random at any
resolution; a real under-sampled boundary resolves into a contiguous band, which is what happened.

### 4. Flip rate tracks boundary frequency, not proximity of PASS% to 50%

If flips were generic boundary-proximity noise you would expect the rate to peak where PASS% is
near 50%. It does not. Titania-Oberon at PASS 74% has **0** daily flips; Ariel-Umbriel at PASS 61%
has the most. The controlling variable is the synodic boundary period vs. the sampling step, not
PASS%. (For contrast: an i.i.d. Bernoulli(p) series has expected singleton fraction p(1-p) ≈ 0.19-0.24
here; the well-resolved Titania-Oberon sits far *below* that at 0, the under-resolved Ariel pairs
sit *above* it — i.e. more alternation than random, the sub-Nyquist signature.)

## The secondary population (#312-only converged threshold FAILs) is also genuine, different mechanism

A small, #312-exclusive subpopulation (~1-4%) fails mode (i): all legs converged but the V4-vs-V3
drift-agreement exceeded the 50,000 km floor. Their drift values (50,084 / 52,952 / 54,644 /
57,560 / 58,059 / 71,432 / 84,702 / 90,617 / 97,390 / 103,413 km) span from just-above-floor to
2× the floor — several are far from any knife-edge. These are real high-perturbation epochs where
the true ephemeris drift genuinely diverges from the Kepler-only V3 reference, a threshold on a
different continuous metric. Not an artifact, not knife-edge-only. This is exactly the distinction
`_aggregate_failure_mode`'s docstring preserves ("the FAIL is real dynamical drift, not a solver
artifact").

## Implication for #567 step (4) — read pass-rates as a synodic duty cycle, not a raw window

Because the boundary is a fast periodic (synodic) structure, **the raw PASS% at any fixed grid is
grid-resolution-dependent and is NOT directly a `validity_window` in wall-clock epoch.** A
quasi_cycler here does not have a contiguous multi-day validity window; it has a **feasible synodic
PHASE band that recurs every synodic period** (3-25 days depending on pair; halved for the
symmetric zero-offset pairs). The physically meaningful, grid-independent figures are:

1. the **feasible duty cycle** = fraction of synodic phase for which the cycle is feasible. Because
   a daily step over a year is an irrational fraction of every pair's synodic period, the daily
   samples are approximately phase-uniform, so the **daily PASS% is itself a good estimator of the
   duty cycle** (e.g. Ariel-Umbriel ~61%, #312 ~78%) — *provided it is read as a duty cycle, not as
   epoch-band robustness*;
2. the **synodic boundary period** (per table above), which sets how often a feasible launch
   opportunity recurs and how wide each opportunity is (≈ duty_cycle × synodic_period).

The `detect_isolated_singleton_anomalies` flip% must **NOT** be reported as an
artifact-suspicion/instability metric for these Uranian candidates: it is confirmed expected
physical aliasing and carries no information beyond "the boundary period is short relative to the
grid." Reporting it as instability would wrongly discount otherwise-legitimate duty cycles.

## Concrete next step

No new bugfix is warranted — this is a re-characterization, not a code defect. Step (4) is
**unblocked** and should proceed under the duty-cycle framing above rather than the raw
pass%-and-flip% reading that motivated the block. A dedicated task (**#568**) is allocated to carry
out the corrected step-(4) writeback-readiness characterization (per-candidate feasible duty cycle
+ synodic boundary period, with the flip% explicitly demoted from an artifact flag to an expected
aliasing diagnostic). The catalogue writeback itself remains a separate step after that.
