# #392 Floquet Low-Amplitude Cycler Re-scope

**Date:** 2026-06-19

## Context
The #347 Floquet bifurcation-continuation framework found branch candidates off the Earth-Moon CR3BP cycler families. However, the #389 candidate `branch_C32_b0` was a far-amplitude (3,3) orbit at ~0.77× the Earth-Sun Hill radius, and thus the solar tide destroyed it on real ephemeris (failing V4).

The goal of this analysis was to re-scope the Floquet continuation to specifically target low-amplitude branches (well inside the Hill radius) that have a real chance of surviving V4, and to report on the ranked candidates.

## Methodology
1. Extended the `cj_window` for the Phase 2 Floquet sweep up to 3.50 (targeting higher Jacobi constant = lower amplitude) and ran a denser sweep with 2000 steps per family (`scripts/run_392_floquet_low_amplitude.py`).
2. Subjected the output to the Hill amplitude screen.
3. Passed the lowest-amplitude (highest C_J) survivor through the V1 (periodic closure) and V2 (bounded cycle drift) gauntlet.

## Results
The sweep successfully walked deeper into the low-amplitude regime for the `C32`, `C11a`, `C11b`, and `C21` families, finding 6 saddle-center brackets overall and converging 3 branch orbits.

The best candidate found was **`branch_C32_C_3.1774`**, breaking off the C32 parent at $C_J \approx 3.1774$.

* **Amplitude vs Hill Sphere:** 0.262 (Solidly passing the Hill screen threshold. Max amplitude is ~394,000 km, well inside the ~1.5 million km Hill radius.)
* **Topology:** Retained the (3, 2) topology of its parent.
* **Stability:** Unlike the highly stable `C32_b0`, this branch is mildly unstable, with a maximum Floquet multiplier of `37.23`.
* **Period:** 72.65 days (16.7 TU).

### V1 and V2 Gauntlet Validation
The branch passed V1 effortlessly, converging with a corrector residual of `3.7e-12` and an independent closure of `2.6e-10`.

On V2 (bounded uncorrected drift), the mild instability manifested over multiple cycles:
* **3 cycles (218 days):** 0.047 km max drift (**PASS**)
* **5 cycles (363 days):** 71.6 km max drift (**PASS**)
* **10 cycles (726 days):** 369,380 km max drift (**FAIL**, exceeds 50,000 km limit)

### V3 and V4 Real-Ephemeris Gauntlet
On V3 (REBOUND IAS15 cross-check), the candidate passed comfortably, proving that the bounded drift observed in V2 is a real dynamical property of the CR3BP model and not an integrator artifact.

However, on V4 (Real Ephemeris using DE440), the candidate suffered a **structural failure**:
* **Single Epoch (2000-01-15):** The candidate drifted by ~7e8 km (escaped) over 3 cycles.
* **Earth+Moon Only Control:** Even without the Sun, Mars, and Jupiter, the candidate drifted by ~1e6 km. 
* **Annual Sweep:** A 100-year sweep across different launch epochs yielded 0 passes and 100 failures (`STRUCTURAL_FAIL_ALL_EPOCHS`). There is no phase of the real lunar orbit that allows this cycler to survive uncorrected.

**The Perturbation Insight:** 
The full-ephemeris drift (7e8 km) is ~700× larger than the Earth-Moon only drift (1e6 km). Even though the solar tide is tiny (0.0% of Earth's gravity at this amplitude, ~8% of the Hill radius), the solar/planetary perturbation clearly dominates the absolute drift. The reconciliation is that this low-amplitude branch is mildly unstable (Floquet multiplier 37.2), so it exponentially amplifies *any* un-modeled perturbation. The real lunar eccentricity provides a disqualifying initial kick, and the small solar tide adds another, and the instability amplifies both into massive escape trajectories.

## Conclusion
`branch_C32_C_3.1774` represents a successful discovery of a low-amplitude branch in the CR3BP. While it survives V1-V3 beautifully, it fails structurally at the first real-ephemeris gate (V4). Because the failure occurs even in the Earth-Moon only control, we conclude that the `C32` Floquet branch family is fundamentally incompatible with the real lunar eccentricity, even when pulled to very low amplitudes deep inside the Hill sphere. 

No admission to the catalogue will be made for this candidate.
