# Digest — Stone & Miner, Voyager 2 encounters with Uranus (1986) and Neptune (1989)

**Date:** 2026-06-23 AET
**Task:** #429 (flyby-altitude design floors / sourced altitudes for GAP bodies).
**Corpus files digested (private paper corpus):**
- `stone-miner-1986-voyager-2-encounter-uranian-system-science-233-4759.pdf`
  — E. C. Stone & E. D. Miner, "The Voyager 2 Encounter with the Uranian System,"
  *Science* **233** (4759), 39–43, 4 July 1986. (text-layer PDF; 25k words)
- `stone-miner-1989-voyager-2-encounter-neptunian-system-science-246-1417.pdf`
  — E. C. Stone & E. D. Miner, "The Voyager 2 Encounter with the Neptune System,"
  *Science* **246** (4936), 1417–1421, 15 December 1989. (text-layer PDF; 5k words)

Both are the lead **overview/summary** papers in the respective *Science* special
issues — Stone (Voyager Project Scientist) & Miner (Deputy) summarise the
encounter geometry, planetary findings, ring system, satellites and magnetosphere.
They are NOT the navigation/trajectory papers; per-satellite closest-approach
*altitudes* (other than the planet C/A and a couple of mass-determination flyby
distances) are not tabulated here.

---

## Overview

These two papers anchor the **as-flown** (single-trajectory) closest-approach
distances for the two ice giants and for Triton — the only spacecraft data that
exists for these bodies. A flown C/A is a single realised trajectory point, **NOT
a design floor**: it tells us what one mission did, not the minimum-safe altitude
a designer would adopt. They are recorded `type: observed-flown`, `confidence:
sourced` (the C/A *distances from centre* are stated in the paper text; the
altitude is our arithmetic = stated-C/A − body-radius, with the radius source
noted).

The companion **Heaton & Longuski 2003** Galileo-style Uranian-satellite tour
paper (already in corpus) supplies the only **design floor** for the Uranian
moons: a >50 km flyby-altitude constraint (Table 4), with a full 40-flyby tour
table (Table 5) of per-encounter design altitudes. That is digested inline below
since it closes the Uranian-moon GAP that the Stone & Miner overview does not.

---

## Per-body closest-approach / altitude table

### Stone & Miner 1986 (Uranus) — *Science* 233:39–43

| body | C/A distance from centre | body radius (source) | altitude_km | type | citation | confidence |
|---|---|---|---|---|---|---|
| Uranus | 107,000 km | 25,559 km (IAU 2015, `constants.py`) | **81,441** | observed-flown | Science 233:39 abstract ("Closest approach, at 107,000 km from … [Uranus]") | sourced |
| Miranda | 28,260 km (mass-determination Doppler pass; "passed within 28,260 km of the satellite") | 235.8 km (`satellites.py`) | **~28,024** | observed-flown (mass pass, not a targeted CA) | Science 233:411 ("the mass of Miranda was determined … as it passed within 28,260 km of the satellite") | sourced (note: this is the Doppler-mass closest distance, not necessarily the trajectory CA) |

Uranus visible equatorial radius stated in the paper (p.40): 25,600 km
(slightly larger than the IAU 2015 1-bar 25,559 km; altitude using the paper's
own 25,600 km would be 81,400 km — within rounding of 81,441).

**Satellite size table (Stone & Miner 1986, Table 3, p.411)**: the first `(km)`
column gives Miranda 484±10, Ariel 1160±10, Umbriel 1190±20, Titania 1610±10,
Oberon 1550±20, 1985U1 (Puck) 170±10. These are **diameters** — they are ~2× the
modern JPL *radii* in `satellites.py` (Miranda r 235.8, Ariel 578.9, Umbriel
584.7, Titania 788.9, Oberon 761.4; e.g. Miranda diameter 484 → r 242 ≈ JPL
235.8). For any altitude arithmetic use the modern JPL radii from `satellites.py`
(authoritative for code). The 1986 paper does not give per-moon CAs anyway, so no
moon altitude is derived from it.

**GAP not closed by this paper:** Ariel / Umbriel / Titania / Oberon
per-satellite Voyager-2 CA *distances* are NOT in this overview paper (they live
in the navigation paper). Only Uranus-planet C/A and the Miranda mass-pass
distance are stated. The Uranian-moon DESIGN floor is supplied by Heaton &
Longuski 2003 (below), not by Voyager flown values.

### Stone & Miner 1989 (Neptune) — *Science* 246:1417–1421

| body | C/A distance from centre | body radius (source) | altitude_km | type | citation | confidence |
|---|---|---|---|---|---|---|
| Neptune | 29,240 km | 24,764 km (IAU 2015, `constants.py`) / 24,760 km (paper 1-bar radius) | **4,476** (4,480 vs paper radius) | observed-flown | Science 246:1418 ("Closest approach, at 29,240 km from the center of [Neptune]") | sourced |
| Triton | 39,800 km | 1,352.6 km (`satellites.py`; paper diameter 2705±6 km → r 1352.5) | **38,447** | observed-flown | Science 246:1418 ("Triton closest approach occurred at 0910 UTC at a distance of 39,800 km from the [centre]") | sourced |

Neptune 1-bar equatorial radius stated in the paper (Table, p.1418): 24.76 ×10³
km = 24,760 km. Triton **diameter** in the satellite table (Table 3, p.1420):
2705±6 km → radius 1352.5 km (matches `satellites.py` 1352.6).

Proteus (1989N1) satellite-table diameter 400±20 km → radius ~200 km, but its
orbital distance from Neptune (117,600 km) is listed, NOT a flyby CA. **Voyager 2
did not perform a close targeted flyby of Proteus** — no Proteus CA altitude
exists in this paper. Proteus remains a GAP.

### Heaton & Longuski 2003 (Uranian satellite tour) — DESIGN floor, inline

`heaton-longuski-2003-feasibility-galileo-style-tour-uranian-satellites-jsr-doi-10.2514-2.3981.pdf`
(JSR, doi:10.2514/2.3981; AAS 01-468). Galileo-style Ariel-orbiter tour using
Titania/Umbriel/Ariel/Oberon flybys. Already in corpus (OCR'd .txt); Table 4 and
Table 5 vision-read for this digest.

**Table 4 — Guidelines & constraints summary (Ariel orbiter):**
- Arrival V∞ < 1 km/s
- Periapsis constraint R_p > 4 R_U
- Initial inclination < 20 deg
- Nontargeted flybys > 25,000 km
- **Flyby altitude > 50 km**  ← the DESIGN floor for Uranian-moon flybys
- Time of flight < 2.5 years

**Table 5 — Tour U00-01 summary (40 flybys), per-satellite minimum design altitude:**

| moon | min design flyby altitude (km) | range across tour (km) | type | confidence |
|---|---|---|---|---|
| Ariel | 35 (event 12; below the 50 km nominal — design example) | 35 – 651 | design-floor (50 km constraint) / design-CA | sourced |
| Umbriel | 54 (event 14) | 54 – 432 | design-floor (50 km) / design-CA | sourced |
| Titania | 58 (event 3) | 54 – 2189 | design-floor (50 km) / design-CA | sourced |
| Oberon | 109 (event 19) | 109 – 584 | design-floor (50 km) / design-CA | sourced |

(Tisserand graphs Figs. 1–2 use 50-km-altitude flyby tick-marks at Uranus.)
Miranda is NOT used as a flyby body in this tour (innermost; the tour targets the
outer four for an Ariel rendezvous). The Heaton & Longuski **50 km** flyby-altitude
constraint is therefore the sourced design floor for all four major Uranian moons.

### Genova 2016 (PADME) / Conte & Spencer 2018 — Phobos / Deimos, inline

`genova-2016-phobos-deimos-PADME-trajectory-AIAA-2016-5681.pdf` (AIAA 2016-5681):
- PADME design requirement: low-altitude passes of Phobos and Deimos at
  **2–10 km** altitude (mission requirement, p.1 abstract + requirements table
  rows 4: "low altitude passes … between 2 and 10 km"). → design-floor ~2 km.
- Sourced flown cross-reference: Mars Global Surveyor performed 14 Phobos flybys
  within an altitude range of **89–213 km** (p.1, citing prior MGS results). →
  observed-flown range.

`conte-spencer-2018-mission-analysis-earth-to-mars-phobos-dro-...049.pdf`:
- Mars-Phobos Distant Retrograde Orbits with x-amplitude A_x ∈ [15, 300] km;
  DRO solutions pass "with Phobos at altitudes less than 10 km" (p.1). These are
  DRO orbit amplitudes / proximity, NOT targeted flybys — note the distinction.
  → design proximity <10 km (A_x ≥ 15 km the smallest DRO considered).

---

## Discipline notes

- **Flown ≠ floor.** Every Voyager-2 value above is `observed-flown`: a single
  realised CA, not a minimum-safe design altitude. The only sourced *design*
  floors closed here are Heaton & Longuski's 50 km (Uranian moons) and Genova's
  2 km / Conte's <10 km (Phobos/Deimos).
- **Altitude = sourced C/A-distance − body-radius.** The C/A distances are
  literature-stated (`sourced`); the subtraction uses code/IAU radii (radius
  source noted per row). Per the discipline this stays `sourced` (inputs sourced)
  with the radius provenance recorded.
- **GAPs still open after this pass:** Proteus (no flown CA; only orbital
  distance), Pluto moons (Charon/Nix/Hydra), Mercury sourced *design min*,
  small irregulars. Uranian-moon per-Voyager-flown CAs are not in the Stone &
  Miner overview (navigation paper would have them) — but the *design* floor is
  now closed via Heaton & Longuski, which is the cycler-relevant value.
- **No code touched.** This task is data/docs only; `safe_alt_km` floors are NOT
  changed (a separate follow-up after reviewing the design-floor values).
