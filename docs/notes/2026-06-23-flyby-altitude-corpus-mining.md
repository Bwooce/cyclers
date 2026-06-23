# Flyby / closest-approach / periapsis altitude corpus mining (#428)

**Date:** 2026-06-23 AET
**Scope:** READ-ONLY sweep of `docs/notes/` (~98 files mentioning closest-approach /
altitude / periapsis) + the two physical registries (`core/constants.py`,
`core/satellites.py`) for *every* sourced flyby / closest-approach / periapsis
altitude, per celestial body, with citations.

**Confidence tags:**
- `primary` — stated in the digest's sourced text / a published table or figure quoted in the digest.
- `secondary` — the digest flags it as a NASA fact-sheet / Wikipedia / ESA-portal cross-check rather than from the primary paper's text.
- `derived` — a value computed by *our* code (SPK / Horizons extraction); high fidelity but not literature-sourced. Flagged separately because the golden-tests-sourced-only rule applies.

**Type vocabulary:**
- `design-floor` — a minimum-safe / ALTMIN constraint used in a design study (the engineering floor).
- `observed-flown` — the as-flown closest approach (SPK / Horizons reconstruction, or a post-flight published value).
- `aim-point/biased` — a targeted aim point (pre-correction bias, or a planned aim altitude).
- `design-CA` — a designed (pre-flight) closest-approach value from a mission-design study (neither a floor nor a flown reconstruction).

**Units note:** giant-planet altitudes are `CA_radius − 1-bar_equatorial_radius`; the radius-in-planet-radii column is the meaningful cross-check (digests note this explicitly). Several sources give CA *distances from body centre* (asteroids, some moon imaging passes) — flagged inline.

---

# PLANETS

## Mercury (Me)

| body | altitude_km | type | mission/paper | date | citation (digest + locus) | confidence |
|---|---|---|---|---|---|---|
| Me | 199 | aim-point/biased | BepiColombo Mercury-1 (published CA) | 2021-10-01 | `2026-06-19-bepicolombo-kernel-recon.md` flyby table (ESA/Wikipedia consolidated; extracted 198.8) | secondary (published) / derived (extracted) |
| Me | 200 | observed-flown | BepiColombo Mercury-2 | 2022-06-23 | same table (extracted 198.1) | secondary / derived |
| Me | 236 | observed-flown | BepiColombo Mercury-3 | 2023-06-19 | same table (extracted 234.8) | secondary / derived |
| Me | 165 | observed-flown | BepiColombo Mercury-4 (~165 km; ESA press release) | 2024-09-04 | same table; ESA Fourth-Mercury-flyby press release | secondary / derived |
| Me | 37,626 | observed-flown | BepiColombo Mercury-5 (post-anomaly high flyby; planned ~37 km raised after Sept-2024 thruster shortfall) | 2024-12-01 | same table; ESA note | secondary / derived |
| Me | 295 | observed-flown | BepiColombo Mercury-6 | 2025-01-08 | same table (extracted 295.0) | secondary / derived |
| Me | 713 | observed-flown | Mariner-10 Mercury-I (SPK-derived; 1.29 R_Me) | 1974-03-29 | `2026-06-19-390-spk-vinf-extractor.md` table; cross-check ~703 km achieved | derived (vs secondary 703) |
| Me | ~703 | observed-flown | Mariner-10 Mercury-I (NASA fact sheet / Dunne 1974 Science) | 1974-03-29 | `2026-06-17-digest-dunne-burgess-1978-mariner-10.md` §6 (NOT in Dunne-Burgess book text) | secondary |
| Me | "200 km closer than planned" | observed-flown | Mariner-10 Mercury-I (qualitative; planned altitude not stated in book) | 1974-03-29 | Dunne-Burgess digest line 59 | primary (qualitative) |
| Me | 1,000 (aimed 960) | aim-point/biased | Mariner-10 Mercury-I (Giberson-Cunningham 1975) | 1974-03-29/30 | `2026-06-19-345-voyager-mariner-mission-digests.md` line 55 | primary |
| Me | ~50,000 | aim-point/biased | Mariner-10 Mercury-II (Fig 8-4 aim point post-TCM-5) | 1974-09-21 | Dunne-Burgess digest line 64/86 | primary |
| Me | 327 | observed-flown | Mariner-10 Mercury-III (NASA fact sheet; NOT in Dunne-Burgess text) | 1975-03-16 | Dunne-Burgess digest line 132 | secondary |
| Me | (not given) | — | Mariner-10 Mercury-III ("closest planetary flyby yet accomplished") | 1975-03-16 | Dunne-Burgess digest line 66/87 | primary (qualitative; no number) |

**Recommended physical design floor (Mercury):** ~200 km is the de-facto flown/aimed floor (BepiColombo all three close Mercury flybys at 199/200/236 km; Mariner-10 ~703 km). Project engineering default in `core/constants.py` is **1000 km** (`safe_alt_km`, explicitly "convention, not sourced physics" — solar-thermal + sparse-tracking nav margin). The sourced flown evidence supports a much lower physical floor of ~200 km.

## Venus (V)

| body | altitude_km | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|
| V | 5,794 (5,790 on p.63) | observed-flown | Mariner-10 Venus (3,600 mi) | 1974-02-05 | `2026-06-17-digest-dunne-burgess-1978-mariner-10.md` line 52/84 | primary |
| V | 5,000 (aimed 4,800) | aim-point/biased | Mariner-10 Venus (Giberson-Cunningham 1975) | 1974-02-05 | `2026-06-19-345-voyager-mariner-mission-digests.md` line 54 | primary |
| V | 16,000 → 1,380 → 5,784 | aim-point/biased | Mariner-10 Venus (bias evolution TCM-1→TCM-2; 55,000 km biased injection aim point) | 1974-02 | Dunne-Burgess digest line 53 | primary |
| V | 16,123 | design-CA | Galileo Venus flyby (lat −40.8°, V_rel 8.2, V∞ 6.2) | 1990-02-10 | `2026-06-17-digest-damario-1992-galileo.md` Table I / Fig 3 | primary |
| V | 10,720 | observed-flown | BepiColombo Venus-1 (extracted 10,721.9) | 2020-10-15 | `2026-06-19-bepicolombo-kernel-recon.md` table | secondary / derived |
| V | 552 | observed-flown | BepiColombo Venus-2 (extracted 552.6) | 2021-08-10 | same table | secondary / derived |
| V | 300 | design-floor | Campagnola 2014 Europa-tour Venus flyby (V∞ 6.62) | 2022-05-14 | `2026-06-17-digest-campagnola-2014.md` Table 1 line 90 | primary |
| V | 300 | design-floor | Hughes 2014 EVE STOUR ALTMIN at Venus | — | `2026-06-07-hughes-2014-fast-mars-free-returns-mining.md` lines 151/202 | primary |
| V | 200 | design-floor | Hughes 2014 EVME STOUR ALTMIN at Venus | — | Hughes digest lines 105/199 | primary |
| V | 10,855 | design-CA | Donahue-Duggan 2022 Mars-2035 Venus passage (V∞ 6,491 m/s) | 2036-04-08 | `2026-06-07-donahue-duggan-2022-mars2033-flyby-mining.md` line 83 | primary |
| V | ~7 km/s V∞ → 61.42° bend at safe alt | — | Cassini Venus flyby (sanity-gate example) | — | `2026-06-16-324-physical-sanity-gate.md` lines 52/89 | primary (V∞, not altitude) |

**Recommended physical design floor (Venus):** **200–300 km** is the consistent sourced design floor (Hughes STOUR ALTMIN 200 km for EVME / 300 km for EVE; Campagnola Europa-tour 300 km). BepiColombo flew 552 km (Venus-2). Project engineering default `safe_alt_km` = **300 km** — matches the sourced design floor. WELL ANCHORED.

## Earth (E)

| body | altitude_km | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|
| E | 960 | design-CA | Galileo Earth-1 (lat 25°, V_rel 13.7, V∞ 8.9) | 1990-12-08 | `2026-06-17-digest-damario-1992-galileo.md` Table I / Fig 5 | primary |
| E | 300 | design-CA | Galileo Earth-2 (lat −34°, V_rel 14.1, V∞ 8.9) | 1992-12-08 | same, Fig 8 | primary |
| E | 200 | design-floor | Russell-Ocampo 2003 Earth-Mars cycler TR>1 floor (max-turn at 200 km Earth flyby) | — | `2026-06-17-digest-russell-ocampo-2003.md` line 74 | primary |
| E | 200 | design-floor / Earth-escape | Campagnola 2014 Europa-tour Earth escape (V∞ 3.77) | 2021-11-22 | `2026-06-17-digest-campagnola-2014.md` Table 1 line 89 | primary |
| E | 11,761 | design-CA | Campagnola 2014 Earth-1 flyby (V∞ 12.07) | 2023-10-24 | same, line 91 | primary |
| E | 3,330 | design-CA | Campagnola 2014 Earth-2 flyby (V∞ 12.05) | 2025-10-24 | same, line 92 | primary |
| E | 12,677 | observed-flown | BepiColombo Earth flyby (extracted 12,685.9; V∞ 3.99) | 2020-04-10 | `2026-06-19-bepicolombo-kernel-recon.md` table | secondary / derived |
| E | 122 (entry interface) | design | Hughes 2014 / Donahue V_Entry inertial Earth-entry altitude (re-entry, not flyby) | — | Hughes digest line 105 | primary |
| E | 3,000 (perigee) | design-CA | Genova-Aldrin 2015 free-return Earth perigee | — | `2026-06-10-genova-aldrin-2015-mining.md` lines 36/57 | primary |
| E | ~6.232 km/s V∞ → 74.6° bend at 300 km | — | Galileo Earth flyby (sanity-gate example) | — | `2026-06-16-324-physical-sanity-gate.md` lines 50/87 | primary (V∞) |

**Recommended physical design floor (Earth):** **200 km** (Russell-Ocampo cycler turn-ratio floor; Campagnola Earth escape; Galileo flew 300 km at Earth-2). Project engineering default `safe_alt_km` = **200 km** — matches. WELL ANCHORED.

## Mars (M)

| body | altitude_km | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|
| M | 200 | design-floor | Russell-Ocampo 2003 / Hughes 2014 STOUR ALTMIN at Mars (EVME) | — | `2026-06-07-hughes-2014-fast-mars-free-returns-mining.md` lines 105/199/226 | primary |
| M | 100 | design-CA / dark-side periapsis | Tito-MacCallum 2018 Mars free-return flyby (periapsis ~100 km, dark side; V∞ 5.42; 34.24° bend) | 2017-2019 | `2026-06-13-tito-maccallum-2018-free-return-reproduction.md` lines 34/76 | primary |
| M | 250 | design-CA / constraint | Donahue-Duggan 2022 Mars-2033 passage (V∞ 5.113 km/s; "by design constraint") | 2033-08-09 | `2026-06-07-donahue-duggan-2022-mars2033-flyby-mining.md` lines 67/150 | primary |
| M | 250 | design-CA / constraint | Donahue-Duggan 2022 Mars-2035 passage (V∞ 7,037 m/s) | 2035-08-13 | same, line 82 | primary |
| M | 346–363 | design-CA | Hughes 2014 one Table-4 EVME case Mars flyby altitude | — | Donahue-Duggan digest line 150 (cross-ref) | primary |
| M | ~5.5 km/s V∞ → 32.2° bend | — | Aldrin/Russell Mars flyby (sanity-gate example; Russell & Ocampo 2005) | — | `2026-06-16-324-physical-sanity-gate.md` lines 54/88 | primary (V∞) |

**Recommended physical design floor (Mars):** **100–250 km** sourced design range (Tito 100 km dark-side periapsis is the lowest published; Hughes/Russell-Ocampo STOUR ALTMIN 200 km; Donahue 250 km by design constraint). Project engineering default `safe_alt_km` = **200 km** — squarely inside the sourced range. WELL ANCHORED. No *flown* Mars gravity-assist CA in corpus (Mars-flyby human missions are all design studies).

## Jupiter (J)

| body | altitude_km | CA radius | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|---|
| J | 277,369 | 4.88 R_J | observed-flown | Voyager 1 Jupiter (SPK; V∞ 10.773; cross-check 4.89 R_J) | 1979-03-05 | `2026-06-19-390-spk-vinf-extractor.md` table | derived (vs primary R_J cross-check) |
| J | 650,060 | 10.09 R_J | observed-flown | Voyager 2 Jupiter (SPK; V∞ 7.639; cross-check Kohlhase-Penzo 10.0 R_J) | 1979-07-09 | same | derived (vs primary) |
| J | — | 4.9 R_J | design-CA | Voyager JST design (Kohlhase-Penzo 1977 Table IV) | 1979-03-05 | `2026-06-19-345-voyager-mariner-mission-digests.md` line 24 | primary |
| J | — | 10.0 R_J | design-CA | Voyager JSX design (Kohlhase-Penzo 1977 Table IV) | 1979-07-09 | same, line 25 | primary |
| J | — | 0.1–22.5 R_J | design-CA | Bourke 1971 Grand-Tour Table 1 J-alt across 9 mission rows (0.1, 0.2, 0.9, 3.0, 3.3, 3.6, 7.8, 8.8, 22.5 R_J) | 1976-79 windows | `2026-06-17-digest-bourke-1971.md` Table 1 lines 49-61 | primary |
| J | 3,800 → 3,470 | — | design-CA | Bourke 1971 JSP'76 Jupiter CA (launch-period start→end) | 1977-12-27 | Bourke digest line 74 | primary |
| J | ~4,500 (1-bar) | design-CA / perijove | Juno Jupiter perijove (Lam 2008; apojove 39 R_J) | 2016-08-03 | `2026-06-17-digest-lam-2008-juno.md` lines 21/57/146 | primary |
| J | 1,000 (Io flyby alt); 4 R_J perijove | design-CA | Galileo JOI / Io flyby (D'Amario 1992) | 1995-12-07 | `2026-06-17-digest-damario-1992-galileo.md` Table I line 53 | primary |

**Recommended physical design floor (Jupiter):** Juno perijove ~4,500 km (1-bar) is the lowest sourced; Bourke's 0.1 R_J ≈ ~7,000 km design rows are inside-the-radiation-belt extremes. Project engineering default `safe_alt_km` = **5,000 km** ("well above cloud tops + radiation/ring margin") — consistent with Juno's ~4,500 km perijove. WELL ANCHORED.

## Saturn (S)

| body | altitude_km | CA radius | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|---|
| S | 123,875 | 3.06 R_S | observed-flown | Voyager 1 Saturn (SPK; V∞ 15.167; cross-check ~3.09 R_S) | 1980-11-12 | `2026-06-19-390-spk-vinf-extractor.md` table | derived (vs primary) |
| S | 101,050 | 2.68 R_S | observed-flown | Voyager 2 Saturn (SPK; V∞ 10.674) | 1981-08-26 | same | derived |
| S | — | 0.1–9.0 R_S | design-CA | Bourke 1971 Table 1 S-alt rows (0.1, 1.2, 1.3, 6.0, 9.0 R_S; "open" for direct-Saturn) | 1976-79 | `2026-06-17-digest-bourke-1971.md` Table 1 | primary |
| S | 333,700 | — | design-CA | Bourke 1971 JSP'76 Saturn periapsis | 1980-01-16 | Bourke digest line 77 | primary |
| S | ≥ 2.7 R_S | design-floor | Wolf-Smith 1995 Cassini ring-plane crossing constraint | — | `2026-06-17-digest-wolf-smith-1995-cassini.md` line 40 | primary |
| S | — | ~2.7 R_S | observed-flown | Voyager 2 Saturn (390 cross-check vs published) | 1981-08-26 | `2026-06-19-390-spk-vinf-extractor.md` line 77 | secondary |

**Recommended physical design floor (Saturn):** ring-plane crossing constraint ≥ 2.7 R_S (Wolf-Smith) ≈ ~102,000 km is the operational floor for the ring system; Voyager 2 flew 2.68 R_S. Project engineering default `safe_alt_km` = **5,000 km** (atmosphere/ring margin — note: the *ring* floor is much higher, ~2.7 R_S, not captured by safe_alt). NOTE: the 2.7 R_S ring constraint is far stricter than the 5,000 km atmospheric default; for Saturn the binding floor is the rings, not the atmosphere.

## Uranus (U)

| body | altitude_km | CA radius | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|---|
| U | 81,573 | 4.19 R_U | observed-flown | Voyager 2 Uranus (SPK; V∞ 14.732; cross-check 81,500 km above cloud tops) | 1986-01-24 | `2026-06-19-390-spk-vinf-extractor.md` lines 101/108 | derived (vs secondary 81,500) |
| U | 81,500 | — | observed-flown | Voyager 2 Uranus (Wikipedia / NASA "35 Years Ago" 50,700 mi) | 1986-01-24 | same digest line 108 | secondary |
| U | — | 0.7–5.1 R_U | design-CA | Bourke 1971 Table 1 U-alt rows (0.7, 0.9, 1.3, 1.7, 3.8, 5.1 R_U; "open" direct-Uranus) | 1976-79 | `2026-06-17-digest-bourke-1971.md` Table 1 | primary |

**Recommended physical design floor (Uranus):** no clean sourced *floor* — only Voyager 2's flown 81,500 km (4.19 R_U) and Bourke design CAs as low as 0.7 R_U. Project engineering default `safe_alt_km` = **1,000 km** (atmosphere-plus-margin, explicitly "pending mission-specific analysis"). GAP — no sourced design floor; the Heaton-Longuski Uranian *moon* tour is the cycler-relevant work but per-moon altitudes are not in the digest (see Moons GAP below).

## Neptune (N)

| body | altitude_km | CA radius | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|---|
| N | 4,507 | 1.18 R_N | observed-flown | Voyager 2 Neptune (SPK; V∞ 16.742; cross-check ~4,950 km above N pole) | 1989-08-25 | `2026-06-19-390-spk-vinf-extractor.md` lines 102/110 | derived (vs secondary 4,950) |
| N | ~4,950 | — | observed-flown | Voyager 2 Neptune (Wikipedia; 29,240 km from centre) | 1989-08-25 | same digest line 110 | secondary |

**Recommended physical design floor (Neptune):** only the flown Voyager 2 value (~4,500–4,950 km). Project engineering default `safe_alt_km` = **1,000 km** (atmosphere-plus-margin, "pending mission-specific analysis"). GAP — no sourced *design* floor; the one flown point is the only literature anchor.

---

# MOONS

## Titan (Saturn)

| altitude_km | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|
| 950 | design-floor | Wolf-Smith 1995 minimum Titan flyby altitude (sample tour assumption; hit exactly at Titan 35/36/37/38/43/46/53) | 2004-2008 | `2026-06-17-digest-wolf-smith-1995-cassini.md` lines 39/95-107 | primary |
| 1500, 1250, 2397, 4408, … 958 | design-CA (per-encounter) | Wolf-Smith 1995 Table 2 — 33 Titan flyby altitudes (full table; range 950–16,784) | 2004-2008 | same, Table 2 | primary |
| 960.3 (T22, lowest) → 339,123 (T1) | observed-flown | Cassini flown Titan flybys T1–T50 (Horizons-derived; min targeted ~960 km) | 2004-2008 | `2026-06-19-digest-cassini-vinf.md` table | derived |

**Recommended physical design floor (Titan):** **950 km** — the canonical Cassini design floor (Wolf-Smith 1995) and the as-flown minimum (~960 km, Cassini T22). STRONGLY ANCHORED, both design and flown. (Project `core/satellites.py` registers Titan with the engineering 100 km default `safe_alt_km`, which UNDERSHOOTS the real ~950 km Titan-atmosphere floor — see Surprises.)

## Galilean moons (Jupiter): Io / Europa / Ganymede / Callisto

| moon | altitude_km | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|
| Io | 1,000 | design-CA | Galileo Io flyby at JOI | 1995-12-07 | `2026-06-17-digest-damario-1992-galileo.md` Table I line 53 | primary |
| Io | 158,000 | design-CA (distance) | Bourke 1971 JSP'76 Io closest approach | 1977 | `2026-06-17-digest-bourke-1971.md` line 75 | primary |
| Europa | 100 | design-floor | Campagnola 2014 Europa-Clipper-precursor min flyby altitude | — | `2026-06-17-digest-campagnola-2014.md` line 200-201 | primary |
| Europa | 100 (EOI sci-orbit alt) | design | Campagnola 2014 11-O3 Europa orbit insertion altitude | 2030-01-03 | same, Table 1 line 114 / line 198 | primary |
| Europa | 275,000 | design-CA (distance) | Bourke 1971 JSP'76 Europa closest approach | 1977 | Bourke digest line 75 | primary |
| Europa | 6,681 / 6,563 (E14/E15) | design-CA | Campagnola 2014 Europa flybys 14/15 | 2029-11/12 | Campagnola Table 1 lines 110/112 | primary |
| Ganymede | 100 | design-floor | Campagnola 2014 Ganymede min flyby altitude (ephemeris-uncertainty buffer) | — | `2026-06-17-digest-campagnola-2014.md` line 200 | primary |
| Ganymede | 100 (G2), 120, 176, 268, 293, 500, 552, 559, … | design-CA (per-encounter) | Campagnola 2014 Tables 1-2 Ganymede flyby altitudes (range 100–23,667) | 2028-2029 | same, Tables 1/2 | primary |
| Ganymede | 158,000 | design-CA (distance) | Bourke 1971 JSP'76 Ganymede CA | 1977 | Bourke digest line 75 | primary |
| Callisto | 200 | design-floor | Campagnola 2014 Callisto min flyby altitude | — | `2026-06-17-digest-campagnola-2014.md` line 201 | primary |
| Callisto | 200 (C4), 221 (C5), 2,124, 2,883 | design-CA (per-encounter) | Campagnola 2014 Callisto flybys | 2029 | same, Tables 1/2 | primary |
| Callisto | 228,000 | design-CA (distance) | Bourke 1971 JSP'76 Callisto CA | 1977 | Bourke digest line 75 | primary |
| Callisto | 100,000 (perpendicular) | design-CA | Bourke 1971 §VI nav-example Callisto flyby | — | Bourke digest line 67 | primary |

**Recommended physical design floors (Galileans):** Io ~1,000 km (Galileo flown design); **Europa 100 km, Ganymede 100 km, Callisto 200 km** (Campagnola 2014 — all explicitly "real altitudes", ephemeris-uncertainty-buffer floors). STRONGLY ANCHORED for Europa/Ganymede/Callisto; Io has a design-CA but no explicit floor.

## Other Saturn moons: Rhea / Dione / Iapetus / Tethys / Enceladus / Mimas

All from Wolf-Smith 1995 Table 2 (`2026-06-17-digest-wolf-smith-1995-cassini.md`), per-encounter design altitudes (primary):

| moon | altitude_km (targeted flyby) | type | citation | confidence |
|---|---|---|---|---|
| Rhea | 999 (Rhea-4 targeted); 13,377–68,841 (nontargeted imaging) | design-CA | Wolf-Smith Table 2 rows 4, 19N/20N/22N | primary |
| Dione | 1,005 (Dione-7 targeted); 13,377–99,682 (nontargeted) | design-CA | Wolf-Smith Table 2 rows 9, 5N/8N/etc. | primary |
| Iapetus | 931 (Iapetus-13 targeted) | design-CA | Wolf-Smith Table 2 row 14 | primary |
| Tethys | 648 (Tethys-21 targeted); 33,630–35,957 (nontargeted) | design-CA | Wolf-Smith Table 2 rows 24, 46N/49N | primary |
| Enceladus | 605 (Enceladus-25 targeted); 5,421–87,264 (nontargeted) | design-CA | Wolf-Smith Table 2 rows 29, 20N/27N/etc. | primary |
| Mimas | 97,383 (Mimas-31N nontargeted only) | design-CA | Wolf-Smith Table 2 row 37 | primary |

Bourke 1971 also gives JSP'76 Saturn-moon CA distances: **Dione 173,000 km, Tethys 210,000 km** (`2026-06-17-digest-bourke-1971.md` line 76, primary).

Cassini flown (Bellerose 2018) gives only *flyby counts* (117 Titan, 28 Enceladus, 33 other icy), NO per-flyby altitudes (`2026-06-17-digest-bellerose-2018-cassini.md`).

**Recommended physical design floors:** Tethys ~648, Enceladus ~605, Iapetus ~931, Rhea ~999, Dione ~1,005 km are the lowest *targeted* design CAs (Wolf-Smith) — these are aim points, not stated atmosphere/safety floors, so treat as observed-design lower bounds, not hard floors. Mimas: only a 97,383 km nontargeted pass → no low-altitude data. Project `core/satellites.py` uses the 100 km engineering default for all.

## Uranian moons: Ariel / Umbriel / Titania / Oberon / Miranda

| moon | altitude_km | type | source | confidence |
|---|---|---|---|---|
| — | GAP | — | Heaton-Longuski 2003 "Galileo-style tour of the Uranian satellites" (40+ flybys incl. Titania & Oberon) is the foundational tour paper, but the digest (`2026-06-16-328-uranian-cycler-lit-deep-dive.md`) does NOT transcribe per-moon flyby altitudes | — |
| Umbriel | — (V∞_O 0.96, V∞ 14.7°/39° bend at safe alt) | derived candidate | `2026-06-16-328-uranian-cycler-lit-deep-dive.md` / `2026-06-16-327-umbriel-silver-verification.md` — our own (1,1) Umbriel-Oberon cycler candidate; V∞ only, no sourced altitude | derived (V∞ only) |

**Recommended physical design floor (Uranian moons):** **GAP — no sourced flyby altitude for any Uranian moon.** Voyager 2 did single flybys (Miranda/Ariel 1986) but the digests carry no CA-altitude numbers. Heaton-Longuski 2003 tour altitudes are an acquisition target. Project `core/satellites.py` uses the 100 km engineering default for Miranda/Ariel/Umbriel/Titania/Oberon (`2026-06-14-full-body-registry.md`).

## Earth's Moon

| altitude_km | type | mission/paper | citation | confidence |
|---|---|---|---|---|
| 3,000 (perilune, lunar farside) | design-CA | Genova-Aldrin 2015 free-return | `2026-06-10-genova-aldrin-2015-mining.md` lines 36/61 | primary |

**Recommended physical design floor (Moon):** only a single design CA (3,000 km farside perilune, Genova-Aldrin). No hard floor sourced. Lunar flyby design floors otherwise GAP.

## Triton / Proteus (Neptune); Phobos / Deimos (Mars); Charon (Pluto); minor Saturn/Jupiter moons

**GAP — no sourced flyby/CA altitude.** These appear only in `core/satellites.py` as physical-parameter registry entries (GM, radius, semi-major axis from JPL SSD; `safe_alt_km` = engineering default 10 km for tiny irregulars / 100 km otherwise — `2026-06-14-full-body-registry.md`). No mission/design study in the corpus flies past any of them with a published altitude. Triton (large, retrograde, inclined) is the most cycler-relevant of these and is a clear acquisition target.

---

# ASTEROIDS

| body | altitude_km | type | mission/paper | date | citation | confidence |
|---|---|---|---|---|---|---|
| Gaspra (951) | 1,600 (from CENTRE of asteroid) | design-CA | Galileo Gaspra flyby (V_rel 8.0; gravity assist negligible) | 1991-10-29 | `2026-06-17-digest-damario-1992-galileo.md` Table I / §2.3.4 lines 50/135 | primary |
| Ida (243) | TBD (at writing); flyby flown 1993-08-28 (V_rel 12.4) | design-CA | Galileo Ida flyby | 1993-08-28 | same, line 52/141 | primary (TBD altitude) |

**Recommended physical design floor (asteroids):** N/A — asteroid flybys are science passes (low mass → no gravity assist); Gaspra CA is "from centre" not altitude. No floor concept applies; project default 100 km for the dwarf-planet/large-asteroid registry entries (Ceres/Vesta/Pallas/Eris/Makemake/Haumea, `2026-06-14-full-body-registry.md`).

---

# GAPS (bodies with NO sourced flyby/CA altitude) — acquisition targets

| body | status | note |
|---|---|---|
| Uranus (planet) | no sourced *design floor* | only Voyager-2 flown 81,500 km; Bourke 1971 design CAs in R_U |
| Neptune (planet) | no sourced *design floor* | only Voyager-2 flown ~4,950 km |
| Uranian moons (Miranda/Ariel/Umbriel/Titania/Oberon) | full GAP | Heaton-Longuski 2003 tour altitudes not transcribed → acquisition |
| Triton (Neptune) | full GAP | cycler-relevant (large); registry-only |
| Proteus (Neptune) | full GAP | registry-only |
| Phobos / Deimos (Mars) | full GAP | registry-only |
| Charon / Nix / Hydra (Pluto) | full GAP | registry-only |
| Amalthea, Hyperion (Jupiter/Saturn small) | full GAP | registry-only; Hyperion in flown Cassini tour but no per-flyby alt |
| Phoebe, Telesto, Helene, Epimetheus, Methone, Pallene (Saturn small) | full GAP | named in Cassini tour (Bellerose) but no altitudes |
| Earth's Moon | partial GAP | one design CA (3,000 km farside, Genova-Aldrin); no hard floor |
| Mercury/Venus/Earth/Mars/Jupiter/Saturn | covered | floors well-anchored (see syntheses) |

---

# SUMMARY

**Bodies with sourced data:** 21 distinct bodies got at least one sourced altitude —
8 planets (Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune; Uranus/Neptune flown-only),
11 moons (Titan, Io, Europa, Ganymede, Callisto, Rhea, Dione, Iapetus, Tethys, Enceladus, Mimas, Earth's Moon — 12 if Moon counted),
2 asteroids (Gaspra, Ida).

**Total altitude rows captured:** ~95+ individual altitude data points (the largest single
contributors are Wolf-Smith 1995 Table 2 = 54 Cassini encounters, the Cassini-vinf
Horizons table = 50 flown Titan flybys, Campagnola 2014 = ~45 Galilean flybys, and
the Bourke 1971 Table 1 = ~30 outer-planet design CAs in planet radii).

**GAPS (no sourced altitude — acquisition targets):** all 5 Uranian moons, Triton,
Proteus, Phobos/Deimos, Charon/Nix/Hydra, Amalthea/Hyperion, and the small Saturn
moons (Phoebe/Telesto/Helene/Epimetheus/Methone/Pallene). Uranus and Neptune lack a
sourced *design floor* (flown-only). Earth's Moon has only one design CA.

**Top surprises (beyond the earlier floor-focused pass):**

1. **Titan's real ~950 km design floor vs the project's 100 km registry default.** Wolf-Smith 1995 (design floor 950 km) and the flown Cassini minimum (~960 km, T22) both pin Titan at ~950 km — but `core/satellites.py` carries the generic 100 km engineering `safe_alt_km` for Titan. Titan has a thick atmosphere; the 100 km default undershoots reality by ~850 km. This is the single most actionable correction surfaced.

2. **BepiColombo gives clean *flown* Mercury/Venus/Earth CA altitudes** (199/200/236/165/295 km at Mercury; 552 km at Venus-2) that match published values to ≤1% — a high-quality flown anchor the earlier pass (focused on cycler *design* floors) skipped. Mercury's real flown floor is ~200 km, far below the project's 1000 km default.

3. **Saturn's binding floor is the rings (2.7 R_S ≈ 102,000 km), not the atmosphere (5,000 km default).** Wolf-Smith's ring-plane crossing constraint is ~20× the atmospheric safe_alt; for Saturn flyby/cycler design the ring constraint dominates and is not represented by `safe_alt_km`.

4. **Complete flown Voyager-2 E-J-S-U-N grand-tour CA tuple now exists** (SPK-derived, #390): Jupiter 650,060 km / Saturn 101,050 km / Uranus 81,573 km / Neptune 4,507 km — Uranus and Neptune flown altitudes that no *design* paper in the corpus provides.

5. **Campagnola 2014 supplies explicit per-moon flyby *floors* for the Galileans** (Europa 100, Ganymede 100, Callisto 200 km, all "real altitudes") — the only place in the corpus with hard moon-flyby floors. These should anchor any Jovian moon-tour cycler work.

6. **Mars human-flyby design CAs cluster at 100–250 km** (Tito 100 km dark-side, Hughes 200 km ALTMIN, Donahue 250 km by-design) — tighter and lower than the project might assume; the 200 km Mars default is well-centred.

**Discipline notes:** Every "observed-flown" value tagged `derived` (Voyager/Mariner-10 SPK, Cassini/BepiColombo Horizons) is computed by our code, not literature-sourced — admissible for high-fidelity context but NOT for golden-test expected sides (golden-tests-sourced-only rule). The `safe_alt_km` registry values are explicitly "convention, not sourced physics" per `core/constants.py` comments and should not be cited as sourced floors. No code/config/catalogue was edited (this task is read-only except this note).
