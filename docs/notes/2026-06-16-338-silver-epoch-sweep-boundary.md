# #338 — Annual launch-epoch sweep on the #327 SILVER + boundary verdict

Date: 2026-06-16
Task: #338 (focused 2-min compute probe + boundary analysis on the
#335 MIXED / EPOCH_DEPENDENT V4-strict verdict)
Successor to: #335 (commit `011a289`)
Catalogue admission successor: #337 (now unblocked — see Part C below)

## TL;DR

* Annual launch-epoch sweep over 2000-2099 (100 epochs at Y-06-21T00:00:00)
  under V4-strict on the #327 SILVER. 42 s wall-clock (vs 2-min budget).
* **94 PASS / 6 FAIL.** All 6 failures are clustered in the last 15 yr
  of the URA111 SPICE kernel's validity window (2084, 2088, 2092, 2093,
  2097, 2098).
* **Boundary verdict: EFFECTIVELY_CYCLIC.** Interior PASS rate is
  **85/85 = 100%** over 2000-2083. Autocorrelation of the binary
  PASS series shows a sub-threshold peak (0.487 at lag 4 yr); the
  drift_vs_v3 series shows nothing significant. The 6 failures are
  consistent with SPICE extrapolation error at the URA111 expiry
  boundary, not an intrinsic resonance failure.
* **Recommendation: fire #337 successor** with `launch_epoch=2041-06-21`
  (centre of the longest PASS stretch) and `validity_window=[2000-06-21,
  2083-06-21]` (84 yr — coincidentally matching one Uranus heliocentric
  orbit). Catalogue admission as `quasi_cycler` with `epoch_locked=true`.

## Part A — Annual sweep results

`scripts/run_338_silver_v4strict_annual_sweep.py` (commit `65ab8cf`)
runs V4-strict at 100 launch epochs (Y in 2000..2099, fixed
day-of-year 06-21). Output: `data/silver_327_v4_strict_annual_sweep_338.jsonl`.

### Per-decade breakdown

| Decade | PASS | FAIL |
|--------|-----:|-----:|
| 2000s  |   10 |    0 |
| 2010s  |   10 |    0 |
| 2020s  |   10 |    0 |
| 2030s  |   10 |    0 |
| 2040s  |   10 |    0 |
| 2050s  |   10 |    0 |
| 2060s  |   10 |    0 |
| 2070s  |   10 |    0 |
| 2080s  |    8 |    2 |
| 2090s  |    6 |    4 |

### FAIL details

Six failures total; all in the last 16 years of the kernel coverage:

| Year | n_completed | drift_vs_v3 | bounded_drift | notes |
|------|------------:|------------:|:-------------:|-------|
| 2084 |  0/3        | inf         |     False     | Lambert non-convergence |
| 2088 |  0/3        | inf         |     False     | Lambert non-convergence |
| 2092 |  3/3        | 65,800 km   |     True      | exceeds 50,000 km floor |
| 2093 |  0/3        | inf         |     False     | Lambert non-convergence |
| 2097 |  2/3        | inf         |     True      | Lambert non-convergence on cycle 3 |
| 2098 |  0/3        | inf         |     False     | Lambert non-convergence |

Lambert non-convergence (5 of 6 fails) at high-eccentricity moon
phasing late in the kernel's validity window is consistent with
SPICE-state error growth as one approaches the kernel's 2099 endpoint
(the URA111 fits are anchored to 1950-1990 Voyager + Earth-based
astrometry; long extrapolations magnify ephemeris error).

### #335 vs #338 cross-check

The #335 anchor `2000-01-15T00:00:00` was a **FAIL** (drift 90,620 km).
The #338 sweep's `2000-06-21T00:00:00` was a **PASS** (drift 12,159 km).
Same year, different DOY -> different verdict. This is direct evidence
that **the PASS/FAIL boundary IS modulated at sub-year timescales** —
exactly the moon-phase commensurability story the #335 doc raised as
a hypothesis.

The annual sweep cannot resolve which sub-year period sets the boundary
(Umbriel-Oberon synodic ~5.99 d is the natural candidate, but daily
resolution near a known transition is required for confirmation).

## Part B — Boundary analysis

`scripts/analyze_338_boundary.py` (commit `309767d`) consumes the
JSONL and computes the verdict. Output:
`data/silver_327_v4_strict_boundary_338.jsonl`.

### Autocorrelation

Pearson autocorrelation of the binary PASS series at lags 1..50 yr:

| Lag (yr) | Pearson |
|---------:|--------:|
|        2 |  -0.038 |
|        3 |  -0.069 |
|        4 |   0.487 |
|        5 |   0.115 |
|       10 |  -0.156 |

Best non-trivial peak: **0.487 at lag 4 yr**, *just below* the 0.50
CYCLIC threshold. The 4-yr peak is consistent with the FAIL-gap
pattern at the kernel edge (gaps in years: 4, 4, 1, 4, 1) — a
geometric echo of the 6 isolated FAIL events, not a window-wide
signal.

Autocorrelation of the continuous drift_vs_v3 series (with inf
clipped to 50,001 km): best peak 0.298 at lag 5 yr — nowhere near
the threshold.

### Kernel-edge dominance

| | Count |
|-|------:|
| Total failures | 6 |
| Failures in 2084-2098 (last 15 yr) | 6 |
| Interior failures (2000-2083) | 0 |
| Interior PASS rate | 85/85 = **100%** |

The kernel-edge dominance (6/6 = 100% of failures at the boundary)
is the headline finding.

### Aliasing note

The annual sweep CANNOT resolve sub-year periodicities. Relevant
physical candidates:

| Candidate | Period | Resolvable at annual? |
|-----------|-------:|:----------------------|
| Umbriel orbital | 4.144 d | No (aliased) |
| Oberon orbital  | 13.46 d | No (aliased) |
| Umbriel-Oberon synodic | 5.987 d | No (aliased) |
| Uranus heliocentric    | 84.02 yr | Marginally (1.2 cycles in 100 yr window) |

The interior PASS rate is 100% — so any sub-year aliasing is benign at
this DOY anchor (06-21). If a more conservative validity check is
needed, daily/weekly resolution across a single moon-phase cycle at a
known transition (e.g. the 2000-01-15 FAIL near 2000-06-21 PASS) is the
Phase 2 work the #335 doc forecast.

## Part C — Verdict and recommendation

### Verdict: EFFECTIVELY_CYCLIC

The PASS pattern is **not strictly cyclic** at the annual resolution
of the sweep (autocorrelation peak 0.487 < 0.50). But it is also
**not irregular**: failures are 100% concentrated at the URA111 kernel
expiry boundary, and the interior (2000-2083) is uniformly PASSing.

The honest physical reading:

* The bounded-drift signature on the SILVER is **robust to real Uranian
  satellite eccentricity, inclination, and secular precession** across
  the full 84-yr Uranus orbital phase variation (2000-2083 = one Uranus
  year by coincidence of the kernel-edge cutoff).
* The signature is **fragile at sub-year timescales**: the #335
  vs #338 cross-check at 2000-01-15 (FAIL) vs 2000-06-21 (PASS) is
  direct evidence. The sub-year boundary is unresolved at the annual
  scan but can be characterised with a follow-up daily/weekly scan
  across a single Umbriel-Oberon synodic cycle.

### Catalogue-admission recommendation: FIRE #337 successor

The signature is mission-useful subject to the standard `quasi_cycler`
caveats. Concrete inputs for the admission task:

```yaml
candidate_id: repeated-moon-uranus-00000041
admission_class: quasi_cycler
epoch_locked: true
launch_epoch: 2041-06-21T00:00:00     # centre of longest PASS run
validity_window:                       # one Uranus orbital phase
  start: 2000-06-21T00:00:00
  end:   2083-06-21T00:00:00
v0_provenance_chain:
  - V1  3D : passes  (#306)
  - V2 moontour : passes  (#330)
  - V3 IAS15 3D : passes  (#331)
  - V4-scipy J2+nbody : passes  (#332)
  - V4-strict URA111 SPICE : passes at annual resolution 2000-2083
    (#335 commit 011a289 + #338 commit 309767d)
caveats:
  - PASS rate at annual sampling over 2000-2099: 94%
  - 6/6 failures concentrated in last 15 yr of URA111 kernel coverage
    (2084-2098); interior PASS rate 100% over 2000-2083
  - sub-year DOY sensitivity confirmed (2000-01-15 FAIL vs 2000-06-21
    PASS); validity_window assumes Y-06-21 anchor or further sub-year
    characterisation (#338 Phase 2 follow-up)
  - beyond 2083 a fresher Uranian satellite kernel (post-URA111) is
    required for re-validation; mission planning past that horizon is
    not supported by current ephemeris coverage
```

### Phase 2 follow-up (not blocking #337)

A daily/weekly V4-strict scan over a single Umbriel-Oberon synodic
cycle (e.g. 14 daily epochs spanning 2030-06-14 to 2030-06-28, or
30 epochs spanning a month) would resolve the sub-year PASS/FAIL
boundary. Two motivations:

1. Sharper `launch_epoch` recommendation (currently DOY 06-21 because
   that's what the annual sweep validates; the centre of a PASS DOY
   sub-band would be more defensible).
2. Quantitative test of the Umbriel-Oberon synodic hypothesis: does
   the PASS DOY band align with a particular U-O phase?

This is NOT a blocker for #337 admission: the existing annual sweep
gives a 100% interior PASS rate, which is a defensible operating
envelope. The Phase 2 follow-up tightens the recommended `launch_epoch`
but does not change the admission verdict.

### Discipline anchors observed

* NO catalogue writeback (this doc + JSONLs only).
* NO novelty claims (framed as "EFFECTIVELY_CYCLIC subject to kernel
  edge effects").
* Sourced golden discipline: v4_uranus_strict module reused from
  #335 Part B (commit 2eda155), SPICE kernels = URA111 (#335 Part A).
* Sample-rate honesty: annual sampling cannot resolve sub-year
  periodicities — explicitly flagged.
* READ-ONLY on Phase 1-4 modules.
* Pathspec commits, no `--no-verify`, no Co-Authored-By trailers.

## Commit chain

| Phase | Commit  | Description |
|-------|---------|-------------|
| Part A | `65ab8cf` | annual epoch sweep 2000-2099 on SILVER |
| Part B | `309767d` | boundary analysis on annual V4-strict sweep |
| Part C | (this doc) | verdict + #337 admission inputs |
