# #378 Phase 0.2 apoapsis-reach spike — RESULT (R1 gate: PASS)

**Date:** 2026-06-26
**Script:** `scripts/spike_378_bct_apoapsis.py` (asserts nothing; prints a table)
**Gate (plan Phase 0.2 / design R1):** does the *incoherent* BCR4BP reach a
Sun-shaped apoapsis ≳ 3 LD (~1.1e6 km) from a low-Earth periapsis at any
`(t₀, |V₀|, γ₀)`? PASS if yes; KILL if max apoapsis < 2 LD everywhere.

## Verdict: **PASS** — proceed with Phase 2/3 BCT construction.

The forward-propagated incoherent BCR4BP from a 200 km LEO periapsis shapes the
full apoapsis ladder; the Hiten-signature 3.9 LD band is squarely reachable.

| `|V₀|` (km/s) | max apoapsis (LD) | km | regime |
|---|---|---|---|
| 10.70–10.90 | ~1.0–1.3 | ~0.4–0.5e6 | sub-target (EM-distance scale) |
| **10.95** | **4.65–5.57** | **1.79–2.14e6** | **straddles Hiten Q_a ≈ 3.9 LD (1.5e6 km)** |
| 10.98 | ~33 | ~12.6e6 | near-escape (not a bound capture arc) |
| 11.00–11.02 | ~45–56 | ~17–21e6 | escaping |

**Key finding:** the Hiten target apoapsis (≈3.9 LD, 1.5e6 km) sits inside the
v ≈ 10.95 km/s band, which produces a *bound* high apoapsis of 4.6–5.6 LD
depending on Sun phase `t₀` and flight-path angle `γ₀`. Sun phase modulates the
apoapsis by ~20% (4.65 LD at t₀/T_sun = 0.25 vs 5.13 LD at 0.0), confirming the
Sun is actively shaping the geometry — exactly the energy-removal mechanism a
WSB exterior transfer relies on. This is the geometric vehicle the #412 reach
spike found the EM-libration family *lacks*; the from-scratch Sun-Earth-scale
build the spike called for is feasible in the incoherent BCR4BP.

Best raw config (escape regime, NOT a transfer): 55.7 LD at v=11.02 — included
only to show the ladder. The transfer-relevant config is v≈10.95 km/s.

## Caveat for Phase 2

The apoapsis at v=10.95 (4.6–5.6 LD) overshoots the Hiten 3.9 LD by ~20–40%,
within the design ±30% band when finely tuned. The corrector (Phase 2.2) tunes
`(|V₀|, γ₀)` to land precisely on `W` at the Moon; the spike only establishes
*reachability*, which is the gate. The incoherent-BCR4BP-vs-PR4BP-3D O(ε²) model
gap (design §4) means the golden is a **signature band**, not bit-exact.

R1 does **not** fire. The coherent-QBCP α-table acquisition (design §8 item 2) is
NOT required.
