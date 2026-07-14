# #591: ranking validated bounded 3D Sun-Earth structures by out-of-ecliptic excursion

**Date:** 2026-07-14
**Origin:** Fable's salvageable reframe of #589 (see
`docs/notes/2026-07-14-589-fable-design-review.md`) — "maximize out-of-ecliptic excursion SUBJECT TO
remaining bounded/returning," a defined figure of merit over orbits already known to stay bounded,
not a new search. Script: `scripts/analyze_591_bounded_ecliptic_excursion.py` (analysis-only,
correctly outside the `run_*.py` preflight ratchet, matching the `dedup_588_*.py` precedent).

## Method

`characterize()` (in `scripts/run_581_gurfil_reproduction.py`) already computes rmin/rmax (radial
distance from Earth) over 1yr/5yr, but not the out-of-ecliptic excursion (max |z|) alone. Added a
small propagation (same `er3bp_eom`/`geocentric_to_barycentric` machinery, same collision/escape
event convention) over the paper's own 5-year "practically stable" horizon, extracting max |z| in km
and AU. Pool: the 5 known Gurfil-Kasdin (2002) 3D families with nonzero z (J, K, L, M, N from
`TABLE34`) plus #588's cluster 43 (confirmed bounded to 1000+ years by #590's long-horizon check —
the strongest-validated 3D structure in the whole #582-#590 thread).

**Type matters:** L and N are labeled "3D DEO" (Delayed Escape Orbit) by Gurfil-Kasdin's own
convention — i.e. NOT bounded by definition. Confirmed directly: both escape within the 5-year
window itself (t=19.8/28.7 rad ≈ 3.15yr/4.57yr), exactly consistent with their own type label, not a
bug. They are reported separately, not ranked alongside the genuinely bounded types (DRO/ERO).

## Result

**Bounded 3D structures, ranked by max |z| over 5 years:**

| id | type | source | max\|z\| (km) | max\|z\| (AU) |
|---|---|---|---|---|
| 43 | 3D DRO | #588 cluster | **6,604,700** | 0.04415 |
| J | 3D DRO | known family | 5,706,214 | 0.03814 |
| M | 3D ERO | known family | 591,824 | 0.00396 |
| K | 3D DRO | known family | 145,837 | 0.00097 |

**3D DEO (not bounded by type, reported separately):**

| id | type | max\|z\| (km) | escapes within 5yr? |
|---|---|---|---|
| L | 3D DEO | 1,772,342 | yes, at ~3.15yr |
| N | 3D DEO | 1,445,955 | yes, at ~4.57yr |

## Reading

**Cluster 43 — from #588's "not novel, likely Family-J-curve-extension" pool — is actually the
single best out-of-ecliptic excursion among all validated bounded structures**, beating the
published state of the art (Family J) by ~16% (6.60M km vs. 5.71M km), and it does so while being
*far* more durably bounded: #590 confirmed cluster 43 stays bounded to at least 1000 years, while
Family J's own literal published IC escapes at only ~30 years. This doesn't change #588's own
adjudication (cluster 43 is still not claimed as a novel *family* — it's plausibly an unsampled point
on Family J's curve, per #588/#590), but it does mean cluster 43 is a genuinely useful, citable
figure-of-merit result on its own terms: the most Earth-bound-yet-off-ecliptic-reaching validated
orbit currently known in this project's data, better than the published anchor by two independent
metrics (excursion AND longevity).

Families M and K sit far below both (0.004 AU and 0.001 AU) — their own type labels (ERO = "Earth
Return Orbit," near-circular quasi-circular 3D DRO respectively) predict small out-of-plane reach;
this is expected, not a finding.

## Scope note

No new search, no new fitness function, no novelty claim. This is a derived-metric ranking over data
already validated by #581 (known-family reproduction) and #590 (cluster 43's long-horizon stability).
`#591 STATUS: CLOSED.`
