# #588: IC-relation check on the 20 unmatched clusters (pre-adjudication triage)

**Date:** 2026-07-14

Before spending an Opus/Fable adjudication pass, checked whether the 20 clusters with no
same-type/overlapping-range known-family match (see
`docs/notes/2026-07-14-588-gurfil-companion-papers-digest.md`) are nonetheless just further points
on an ALREADY-KNOWN continuous family — several Gurfil-Kasdin families (and Henon's foundational
DRO family f, already anchored in `literature_check.py`) are defined by a characteristic IC
*relation* (e.g. ẏ₀ ≈ -2x₀), not an absolute radius; the paper's own 12 characterization sets each
picked one optimizer-selected point, not a scan along the curve, so a widened search box landing on
a different point of the SAME curve is expected and not novel. This is pure arithmetic on
already-computed `ic_interleaved` vectors (`data/found/583_widened_search/unmatched_20_for_adjudication.json`)
— no propagation, no new compute.

## Result: three clean groups

**Group 1 — 8 DEO clusters (1, 2, 3, 7, 20, 35, 36, 44), ẏ₀/x₀ = -1.94 to -1.99.**
All planar or near-planar, all within 1% of Henon's ẏ₀=-2x₀ relation — the SAME curve Families A, B,
F sit on (and Henon's family f itself, already an anchored corpus item). Radii span 5.4M-40.8M km,
overlapping/extending beyond A's 5.77M-11.7M km band. Our own drift classifier calls these "DEO"
(practically unstable at 5yr) rather than "DRO" — expected: further out on the same curve, weaker
Earth-binding, escapes past the 0.1 AU/5yr threshold sooner. **Read: near-certainly the same
already-known curve sampled beyond where the original 12 sets happened to land — not a new family.**

**Group 2 — 6 of 12 "3D DRO" clusters (39, 40, 42, 43, and borderline 24, 25), ẏ₀/x₀ = -2.0 to -2.5.**
Close to Family J's own defining relation (ẏ₀≈-1.987x₀≈-2x₀, paper's Eq. 8), at radii (1.3M-13.4M km)
spanning both below and around J's 6.9M-8.5M km band. **Read: plausibly the same J-type 3D DRO curve
extended to radii the original set didn't sample; 24/25 (-2.3, -2.5) are a looser fit and worth a
closer look, not a clean dismissal.**

**Group 3 — 6 of 12 "3D DRO" clusters (10, 12, 14, 29, 37, 38), ẏ₀/x₀ wildly divergent (-62.8, +705.8,
-11.5, +27.4, +28.1, -11.6).** These do NOT sit near either known 3D-DRO curve (J's ≈-2, or K's own
≈-5 from its IC). The wild ratios are largely a coordinate artifact — x₀ is very close to zero
(-0.00004 to 0.00294) in all six, making ẏ₀/x₀ numerically unstable rather than meaningfully
different. What IS meaningfully different: **all six sit at 80,274-559,535 km — smaller than the
smallest known 3D DRO family (K, 854,531-939,889 km)**, a radius band the original paper's 12 sets
(and our own #581 reproduction) never sampled for this type. **Read: this is the genuinely
interesting subset** — a possible small-radius 3D DRO regime below the previously-characterized
floor, not obviously explained as a curve extension of J or K.

## Recommendation for the Opus/Fable adjudication

Focus scrutiny on **Group 3** (6 candidates, small-radius 3D DRO below the known floor) — check
internal coherence (do these 6 actually cluster together as one family, or are they 6 unrelated GA
artifacts near a coordinate degeneracy?), stability/periodicity beyond the 1yr window, and a fresh
targeted literature search for "small-radius 3D distant retrograde orbit Sun-Earth" before any
novelty language. Groups 1 and 2 (14 candidates) get a lighter pass — default toward "known-curve
extension, not novel," consistent with this pool's own track record (#577 found 0/36 novelty-clear
on a similarly promising-looking Galilean pool).
