# #491 — V0 ingestion staging: Russell-Strange 2009 + Lynam-Longuski 2011 sourced invariants

**Date:** 2026-06-30. The trust-bearing part of #491 — the **sourced golden invariants** for the
~31 documented cyclers, transcribed directly from the papers (read visually; every value traces
to a printed table cell, [[feedback_golden_tests_sourced_only]]). The mechanical follow-on
(construct catalogue rows from this + coordinate the frozen-census ratchets) is deliberately
NOT rushed at session-tail; this note is the verified basis for it.

## Russell-Strange 2009 — Jovian cyclers (Table 3, p.149; descriptors Table 5, p.151)
`(flyby A, target B, ID, syn-period d, V∞_A km/s, V∞_B km/s, n_legs, period d, min flyby alt @B km)`:

| A | B | ID | syn(d) | V∞_A | V∞_B | n_legs | period(d) | minAlt@B(km) |
|---|---|---|---|---|---|---|---|---|
| Europa | Ganymede | 131 | 7.05 | 2.40 | 4.10 | 2 | 21.2 | 1113 |
| Europa | Ganymede | 159 | 7.05 | 2.45 | 4.11 | 3 | 28.2 | 1261 |
| Ganymede | Callisto | 1 | 12.52 | 3.18 | 3.26 | 3 | 37.6 | 247 |
| Ganymede | Callisto | 5 | 12.52 | 3.24 | 3.34 | 2 | 37.6 | 328 |
| Ganymede | Europa | 5 | 7.05 | 1.66 | 2.57 | 1 | 35.3 | 1819 |
| Ganymede | Europa | 43 | 7.05 | 1.87 | 3.89 | 1 | 14.1 | 8861 |
| Ganymede | Europa | 316 | 7.05 | 3.20 | 3.81 | 4 | 49.4 | 1447 |
| Ganymede | Io | 53 | 2.35 | 3.90 | 3.90 | 2 | 21.2 | 518 |
| Ganymede | Io | 185 | 2.35 | 3.97 | 9.90 | 6 | 49.4 | 603 |
| Ganymede | Io | 403 | 2.35 | 4.29 | 4.34 | 2 | 56.4 | 540 |

Leg descriptors (Table 5, the `g/G/h/f(p:q, transit_day, branch)` formal nomenclature) e.g.
GanCal#1 = `G(1.74871,269.53421,U) g(1.50246,540.88534,L) f(2:1,77.40130,0.03291)`;
EurGan#131 = `G(3.95655,704.35739,U) f(2:1,87.95239,90.0)`.

## Russell-Strange 2009 — Saturnian Titan→Enceladus cyclers (Table 4, p.150; descriptors Table 6)
`(ID, V∞_Titan km/s, V∞_Enceladus km/s, n_legs, period d, petal-period yr, min Titan alt km)`. 19 rows:
37: 2.50/5.17/2/94.4/−3.34/1377 · 145: 2.78/7.33/3/94.4/−3.34/1067 · 183: 2.98/4.66/3/64.5/4.16/3304 ·
207: 3.03/5.08/4/63.0/−3.34/3905 · 217: 3.03/5.08/4/63.0/−3.34/10059 · 227: 3.17/5.98/3/80.9/2.90/1183 ·
231: 3.18/6.04/4/97.4/2.41/2140 · 235: 3.18/6.04/5/97.4/2.41/1843 · 314: 3.59/6.88/2/49.5/1.32/1218 ·
370: 3.69/4.15/3/97.4/2.41/1343 · 492: 4.63/5.61/4/97.4/2.41/4085 · 510: 5.03/7.27/2/80.9/2.90/1852 ·
539: 5.14/8.08/4/97.4/2.41/1821 · 552: 5.16/7.69/2/112.4/6.06/3784 · 572: 5.28/4.69/3/112.4/6.06/1798 ·
586: 5.43/6.68/2/112.4/6.06/3874 · 594: 5.55/5.76/3/127.4/−38.20/6348 · 602: 5.55/5.76/4/127.4/−38.20/1469 ·
624: 5.79/7.97/3/95.9/15.77/1961 · 631: 5.91/7.49/3/127.4/−38.20/4243.
(Saturn G-ring hazard out to r~176,000 km; Titan-Enceladus synodic period 1.50 d. Flagship #235:
9 Enceladus + 45 Titan encounters, 2.4 yr; legs `g g f(2:3,…) f(1:2,…)`.)

## Lynam-Longuski 2011 — IEG triple cyclers (the #480 2nd source)
- **Single-period IEG**: E→(perijove)→I→G→(apojove)→E, period 7.05 d, example depart 2024-08-27,
  total ΔV ≈ 11 m/s (patched-conic ephemeris/MALTO); ballistic version needs −175 km Europa flyby.
- **Indefinitely-repeatable GIPEIPE**: Ganymede-Io-(perijove)-Europa-Io-(perijove)-Europa, orbital
  period 3.5 d / repeats every 7.05 d (Laplace period); lowest-ΔV member.

## Ingestion plan (the mechanical follow-on — careful, not rushed)
1. **Class/level:** all are `orbit_class: cycler`, `validation_level: V0` (sourced, published;
   the in-paper real-eph transition is provenance, not our V1). Jovian/Saturnian → `center` =
   Jupiter/Saturn; `trajectory_regime: ballistic` (R-S) / powered (L-L single-period 11 m/s).
2. **Provenance:** `first_published` = Russell-Strange 2009 (DOI 10.2514/1.36610) / Lynam-Longuski
   2011 (DOI 10.1016/j.actaastro.2011.03.011); `vinf_kms_at_encounters` from the tables above;
   `source_ephemeris` = ideal circular-coplanar (+ in-paper real-eph note).
3. **L-L on #480:** add Lynam-Longuski 2011 as a 2nd corroborating source on
   `hernandez-2017-jovian-ieg-triple-family`; the GIPEIPE + single-period as distinct IEG rows.
4. **Ratchets (the care):** ~31 new rows shifts the census → update EVERY frozen-census ratchet
   (`cycler_class_census`, `validation_tier_census`, `rediscovery`, `validate`) + KNOWN_CORPUS
   anchors; run `uv run pytest tests/data tests/search -q` ([[feedback_catalogue_edits_run_all_ratchets]]).
   Do NOT hand-pick a subset of ratchets.
5. **Dedup check:** confirm none already in catalogue (the GanEur/GanCal/GanIo two-moon cyclers are
   new; cross-check vs the Liang CGE + Hernandez IEG rows so the IEG ones aren't double-counted).

## Status
#491 trust-bearing extraction DONE (sourced invariants captured + verified vs the printed tables).
Row construction + ratchet coordination is the bounded mechanical follow-on against this staging note.
