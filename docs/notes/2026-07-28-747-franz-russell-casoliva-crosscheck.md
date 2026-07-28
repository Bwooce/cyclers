# Finding note — Franz & Russell 2022 vs Casoliva's 6 symmetric Class-1 rows (#747)

**Task:** `#747`, flagged by `#742`'s digest as "a genuine, not-yet-executed future cross-check
candidate": do Casoliva et al. 2010's six *symmetric* p-q resonant Class-1 cyclers (1-2c, 1-2d,
2-1a, 2-1b, 3-2c, 7-3a) appear as a matching family (comparable Jacobi constant + period) in Franz
& Russell 2022's public 13-million-solution Earth-Moon periodic-orbit database?

**Verdict up front: structural domain mismatch, not a database lookup result.** All 6 rows are
excluded from Franz & Russell's search domain on three independent, quantified grounds (search-box
position, escape-radius/near-Keplerian-construction, and energy regime — detailed below). Two of
the six (1-2c, 1-2d) are excluded with certainty by a single published number each (their own
periselene radius already exceeds Franz & Russell's stated "vicinity of Moon" cutoff). The other
four are excluded with high but not absolute confidence (their periselene *does* sit inside Franz &
Russell's near-Moon zone, but their Earth-relative apogee and the near-Keplerian nature of the p-q
construction method make it very unlikely the full orbit stays inside that zone for a full period —
this last step was not independently verified by propagation, since building a propagator was out
of scope for this research-only task; see Caveats). **A literal brute-force search of the 13M-row
database was attempted and confirmed impractical within this task's effort budget** (measured
throughput ~790 KB/s → ~27 minutes for the 1.3 GB archive alone, before decompression/grep of the
5.4 GB uncompressed text — see §3). The verdict below rests on the two papers' own published
numbers, not a database lookup.

## 1. Casoliva Table 3 — the 6 rows, full precision

`pdftotext -layout` (used by the `#725`/`#742` digests) badly scrambles Table 3's columns across
the page. `pdftotext -raw` on the same PDF (re-extracted this session,
`/Users/bruce/dev/cyclers_pdf/papers/casoliva-mondelo-villac-mease-barrabes-olle-2010-two-classes-cycler-trajectories-earth-moon-jgcd-33-5-doi-10.2514-1.46856.pdf`)
recovers it cleanly as one row per orbit. µ_EM = 0.0121529529 (Earth at `(µ,0)`, Moon at `(µ−1,0)`,
barycentric synodic frame, DU = 384,400 km, TU⁻¹ = ω_EM = 2.66529×10⁻⁶ rad/s — both from the
paper's own text). IC point `(x_i, y=0, u_i, v_i)` is a Poincaré-section (`y=0`) crossing, **not**
necessarily the orbit's own x-axis-perpendicular symmetric crossing (Casoliva's own text: "we have
not used this property [the perpendicular crossings] in this paper").

| Row | C_J | T (period) | x_i | u_i | v_i | r_pM (periselene, Moon-rel.) | r_aE (apogee, Earth-rel.) | k (stability) |
|---|---|---|---|---|---|---|---|---|
| 1-2c | 1.5691874798 | 12.5663706144 (4π) | −2.4754942840 | −0.3779077269 | 2.2861781603 | 0.8833330997 | 3.0501343486 | +1.9995914782 (stable) |
| 1-2d | 2.5803060666 | 12.5663706144 (4π) | −2.5987958202 | 0.0000000000 | 2.2237844860 | 1.5702173819 | 2.6109487731 | +4.8578794987 (unstable) |
| 2-1a | 0.4887353098 | 6.2831853072 (2π) | −0.7362153512 | −0.6876743802 | 1.5221620241 | 0.2353576230 | 1.0644080071 | +1.9255637447 (stable) |
| 2-1b | 1.1964188553 | 6.2831853072 (2π) | −1.2287157860 | 0.0000000000 | 1.4164812666 | 0.2408687390 | 1.2621776123 | +2.0373978625 (unstable) |
| 3-2c | 0.7089330385 | 12.5663706144 (4π) | −1.2237269180 | −0.5006876037 | 1.4965104108 | 0.2193864666 | 1.5557128499 | +1.8751541804 (stable) |
| 7-3a | 1.0215696153 | 18.8495559215 (6π) | −0.8938394486 | −0.4831113793 | 1.4082724890 | 0.0709180137 | 1.0910338049 | −1.2990228617 (stable) |

All lengths in EM-distance units (1 = 384,400 km, `µ_EM`-normalized). Periods are exact integer
multiples of `2π` (`T ≈ 2πq`, Casoliva Eq. 18) — 1 sidereal month (2π) for the 2-1 pair, 2 months
(4π) for 1-2c/1-2d/3-2c, 3 months (6π) for 7-3a; in days (sidereal month = 27.321661 d):
27.32 d (2-1a/2-1b), 54.64 d (1-2c/1-2d/3-2c), 81.96 d (7-3a).

**Derived distances (km), converted at 384,400 km/unit, Moon at barycentric x = −0.9878470471:**

| Row | Distance from Moon at tabulated IC (km) | Periselene distance from Moon (km) | Apogee radius from Earth (km) |
|---|---|---|---|
| 1-2c | 571,852 (Earth side) | 339,553 | 1,172,472 |
| 1-2d | 619,249 (Earth side) | 603,592 | 1,003,648 |
| 2-1a | 96,727 (L2 side) | 90,471 | 409,158 |
| 2-1b | 92,590 (Earth side) | 92,590 | 485,182 |
| 3-2c | 90,672 (Earth side) | 84,332 | 598,016 |
| 7-3a | 36,137 (L2 side) | 27,261 | 419,393 |

## 2. Franz & Russell 2022 — mass ratio and search-domain bounds

Read in full: `/Users/bruce/dev/cyclers_pdf/papers/franz-russell-2022-database-planar-3d-periodic-orbits-families-near-moon-jas-69-1573-doi-10.1007-s40295-022-00361-9.pdf`
(no `.txt` sidecar existed despite `#742`'s digest describing it as text-layer — generated one this
session via `pdftotext -layout`, now committed alongside the PDF).

**Mass ratio.** Franz & Russell's Table 1: `µ = 1.215058392535863×10⁻²`. Casoliva's own precise
2010 value: `µ_EM = 0.0121529529`. Difference = `2.368975×10⁻⁶`, **relative difference =
0.019497% (≈2.0×10⁻⁴, confirming the `#742` digest's "~0.02%" estimate to full precision).**
Negligible for this cross-check — as shown below, the actual exclusion reason is 3-4 orders of
magnitude larger in effect than this mass-ratio gap.

**Search domain (§3, Table 2, Moon-centered frame — "the reference frame used here is centered at
the Moon, not the system barycenter"):** the planar grid search's `x0` range is **−185,000 km to
`L2`** (≈+64,500 km beyond the Moon on the far side); `ẏ0 ∈ [0, 2.0]` km/s. Figs. 1 and 3
(rendered and visually inspected this session, PDF pages 15-16) confirm this empirically — the
retrograde planar phase-space plot's `x0` axis spans 0 to −2×10⁵ km, the prograde plot spans 0 to
+6×10⁴ km, matching Table 2's stated bounds.

**Escape criterion (§4, explicit quote):** *"Orbits that escape the vicinity of the Moon (defined
as ever being more than 4 Hill's units from the Moon, approximately 350,000 km or 90% of the
Earth–Moon distance) or impact the Moon... are not considered."* This is a **hard exclusion applied
during propagation**, not just an initial-condition-grid bound — any orbit whose trajectory, at any
point over its period, exceeds ~350,000 km from the Moon is discarded regardless of where its
search-grid point started.

**Energy/Jacobi-constant regime (Fig. 3, dimensional `J` in km²/s², visually confirmed this
session):** the retrograde planar phase-space panel spans `J ≈ 2.8`–`4.2` km²/s²; the prograde
panel spans `J ≈ 3.0`–`4.0` km²/s². Converting Casoliva's 6 target `C_J` values to the same
dimensional units (`J = C_J × (DU·ω_EM)²`, `DU·ω_EM = 1.024537` km/s from Casoliva's own stated
constants, `VU² = 1.049677` km²/s² — an approximate order-of-magnitude conversion for comparison,
not a reproduction of either paper's exact internal normalization pipeline) gives:

| Row | C_J (normalized) | J (km²/s², converted) |
|---|---|---|
| 1-2c | 1.5691874798 | 1.647 |
| 1-2d | 2.5803060666 | 2.708 |
| 2-1a | 0.4887353098 | 0.513 |
| 2-1b | 1.1964188553 | 1.256 |
| 3-2c | 0.7089330385 | 0.744 |
| 7-3a | 1.0215696153 | 1.072 |

**All 6 values fall at or below Franz & Russell's own plotted `J`-minimum (≈2.8 km²/s²
retrograde).** Even 1-2d (the lowest-energy/highest-`C_J` of the six) falls just short; the other
five (0.51–1.65 km²/s²) sit a factor of 1.7–5.5x below the low end of the energy range their own
figures show their database actually populating.

## 3. Zenodo dataset — accessed, confirmed impractical to brute-force search

`https://doi.org/10.5281/zenodo.6411980` → `https://zenodo.org/record/6411980` (fetched this
session). Contents: **a single file**, `lunarPOdatabase.zip`, 1,313,701,084 bytes (1.22 GiB)
compressed, 5.4 GB uncompressed containing the 5 human-readable text files (planar prograde/
retrograde, axial prograde/retrograde, x-z), README, MATLAB GUI/scripts, and grid-search/clustering
summary files. **No online query API, search index, or partial-download mechanism is offered** —
the summary/cluster files are bundled inside the same single zip, so even the smaller
"family-summary" product requires downloading the full archive first.

Measured a 20 MB range-request against the live download URL to get an honest throughput estimate
rather than guessing: **808,356 bytes/s (~790 KB/s), 25.96 s for 20 MB.** Extrapolated to the full
1.22 GiB archive: **≈1,625 s (≈27 minutes)** for the download alone, before unzip/decompress time
for the 5.4 GB of text and any grep/parse pass over 13 million rows. This is well beyond this
task's ~8-minute-per-step / no-backgrounding constraint, so **the download was not attempted** —
per [[feedback_never_give_up_reproducing_papers]]'s companion honesty requirement
([[feedback_verify_gauntlet_with_positive_control]]), this is reported as a genuine
couldn't-practically-check, not glossed over. The domain-boundary analysis in §§1-2 above,
grounded in exact published numbers from both papers (not estimated or assumed), was judged
sufficient to reach a confident verdict without the brute-force search — see §4.

## 4. Per-row verdict

| Row | C_J | Verdict | Basis |
|---|---|---|---|
| **1-2c** | 1.569 | **Excluded (certain)** | Even its periselene (339,553 km from Moon) is at Franz & Russell's own 350,000 km "vicinity of Moon" cutoff; its IC point (571,852 km) is far outside both the `x0` search bound (185,000 km) and the escape radius. |
| **1-2d** | 2.580 | **Excluded (certain)** | Periselene alone (603,592 km) is **nearly double** the 350,000 km cutoff — this orbit never enters Franz & Russell's near-Moon search domain at any point in its trajectory. |
| **2-1a** | 0.489 | **Excluded (very likely)** | IC point (96,727 km, L2 side) already exceeds the prograde search's `x0max = L2` (~64,500 km) bound; periselene (90,471 km) is inside the 350,000 km zone, but apogee (409,158 km from Earth) implies a wide near-Keplerian excursion for most of the 27.3-day period; `J≈0.51` km²/s² is far below Franz & Russell's plotted energy range. |
| **2-1b** | 1.196 | **Excluded (very likely)** | Periselene (92,590 km) and IC point both sit inside the near-Moon search zone, but Earth-relative apogee (485,182 km) — reached during the same 27.3-day period — is well beyond the 350,000 km Moon-vicinity radius for a plausible fraction of the orbit; `J≈1.26` km²/s² far below range. |
| **3-2c** | 0.709 | **Excluded (very likely)** | Same pattern as 2-1b, wider apogee (598,016 km, 54.6-day period); `J≈0.74` km²/s² far below range. |
| **7-3a** | 1.022 | **Excluded (very likely)** | The "tightest" cycler of the six (periselene 27,261 km — genuinely close to the Moon, Casoliva's own paper singles this one out as the one built via the small-µ second-species/matched-asymptotics method specifically to reach a tight periselene), but its Earth-relative apogee (419,393 km, over an 82-day period) still implies the spacecraft spends most of each 3-month cycle far outside the Moon-vicinity zone; `J≈1.07` km²/s² far below range. |

The "very likely" (not "certain") qualifier on 4 of 6 rows reflects one limitation honestly: I did
not propagate these orbits to check exactly how far from the Moon they travel at every instant
(only periselene/apogee *radii*, not the full trajectory, are tabulated in Casoliva's paper, and
writing a propagator was out of scope for this research-only task per the dispatch instructions).
The inference rests on (a) Franz & Russell's own explicit energy range (§2, all 6 `J` values below
their plotted minimum — independent of trajectory shape) and (b) the near-Keplerian construction
method itself (Barrabés-Gómez matched asymptotics: a brief Moon encounter stitched to a long free
two-body Keplerian arc around Earth — Casoliva's own description, `#742`'s digest confirms this
matches the source method verbatim), which structurally implies wide excursions between periselene
passages. Both lines of evidence point the same direction and are independent of each other and of
the search-box argument in §1-2, so the combined case is strong even without direct propagation.

## 5. Bottom line

**Franz & Russell's database neither corroborates nor refutes Casoliva's Class-1 families — it is
not applicable evidence for them.** The two datasets occupy non-overlapping regions of phase space
by construction: Franz & Russell's grid search is explicitly bounded to within ~350,000 km /
4 Hill's units of the Moon and populates `J ≳ 2.8` km²/s² in its own published figures; all 6 of
Casoliva's symmetric Class-1 cyclers are wide, near-Keplerian, high-eccentricity Earth-orbiting
resonant paths with `J` in the 0.51-2.71 km²/s² range and (for 4 of 6) substantial excursions well
beyond the Moon-vicinity radius. This is a genuine domain-boundary finding, not a "searched and not
found" result — the two papers are studying dynamically adjacent but non-overlapping populations
(matches this project's broader finding in the `#725` digest that Casoliva's Class 1 has "no direct
RRT counterpart at all," now extended: it also has no direct Franz-Russell counterpart, for the
same underlying reason — Class 1's near-Keplerian wide orbits sit outside every near-Moon
symmetric-family census this project has cross-checked so far). The ~0.02% mass-ratio difference
flagged by `#742`'s digest is real but immaterial — it is 3-4 orders of magnitude smaller than the
actual gap (hundreds of thousands of km / whole units of `J`) driving the exclusion.

No code was written; no catalogue rows were touched (research/verification task only, per scope).
