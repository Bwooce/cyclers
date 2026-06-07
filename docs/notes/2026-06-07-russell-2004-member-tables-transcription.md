# Russell 2004 dissertation — remaining member-table transcription (Ch.3 + App.C)

Date: 2026-06-07
Source: **Russell 2004 dissertation** ("Global Search and Optimization of Trans-Mars
Cyclers", R. P. Russell). Read via the project reference-PDF library, vision read,
pages cited inline. Clean digital typeset; Ch.3 tables read unambiguously.
Appendix C is dense fixed-width numeric dumps (small monospace) — see §App.C notes
for the honest readability caveat.

## Purpose & scope

Flagged by #142. Tables 3.4–3.11 and Appendix C hold the **circular-coplanar**
(Ch.3) and **ephemeris / DE405** (App.C) cycler member lists. This note transcribes
the **not-yet-noted** remainder **VERBATIM** as a reference NOTE. **No catalogue
writes** were performed by this pass; catalogue-eligibility is recommendation only.

### What is already in the catalogue (do NOT duplicate — cross-reference)

`grep russell data/catalogue.yaml` shows the Ch.3 member lists are already
catalogued as derived rows:

- **Table 3.4** (2/3/4-synodic) → source tag `russell-2004-t34`. Catalogued rows
  include `russell-ocampo-2.1.1+2-case2`, `…-2.3.1+1-case3`, `…-4.3.1-5`,
  `…-3.1.2+1`, `…-2.5.1+0`, plus the Aldrin cycler `1.0.1.-1` carried as the
  primary Aldrin entry. (The catalogue stores **interpreted** fields — AR, TR,
  per-flyby turn angle, surrogate ΔV — not the verbatim table row.)
- **Tables 3.9 / 3.10 / 3.11** (5- and 6-synodic) → source tag
  `russell-2004-t39_311`. Catalogued as ~157 `russell-ocampo-5.*` and
  `russell-ocampo-6.*` rows (counted: 58 five-synodic + 99 six-synodic in
  `data/catalogue.yaml`).
- **Appendix C parent cyclers** (Ch.4 shorthand) → catalogued as `russell-ch4-*`
  rows (e.g. `russell-ch4-6.399G1` lineage, `russell-ch4-4.991gG2`,
  `russell-ch4-8.049gGf2`, `russell-ch4-8.165Gfh-f2`, etc.) with the
  circular-coplanar simple-model values.

The **verbatim tables themselves** (every published digit, the caption text, the
footnotes, the units) were NOT previously written to a note — the prior
method-mining note (`2026-06-07-russell-2004-dissertation-method-mining.md`,
Anchor U4) explicitly deferred Tables 3.4 + 3.5–3.11 and Appendix C. This note
closes that gap.

---

## Model / frame / provenance flags (apply to ALL Chapter-3 tables)

Per v4.2 provenance:

- **Model / fidelity (Ch.3): SIMPLE / circular-coplanar.** Ch.3 builds cyclers in a
  "simple solar system" — circular, coplanar planetary orbits (Earth & Mars), the
  ecliptic is the x-y plane, the initial Earth position is on the x-axis (p.84).
  AR/TR and V∞ are circular-coplanar idealizations, NOT ephemeris values.
- **Center / frame:** heliocentric (Sun-centered), ecliptic x-y plane. Flybys are
  **geocentric (Earth)** gravity assists — Russell §3.8: "The last column shows the
  number of gravity-assisted Earth flybys and the geocentric turning angle
  associated with each one."
- **Epoch:** none — circular-coplanar cyclers are time-independent / generic
  (Earth starts on x-axis). No J2000 epoch in Ch.3.
- **Cycler nomenclature `p.h.s.i`** (front-matter Nomenclature, p.xiv):
  - `p` = period of cycler in synodic periods (p ≥ 1).
  - `h` = total number of half-years allotted for full- or half-revolution
    free-returns during one cycler period (h ≥ 0).
  - `s` = total number of identical generic returns during one cycler period (s ≥ 1).
  - `i` = the i-th solution from a multiple-revolution Lambert problem. |i| = number
    of complete revolutions made by the free-return transfer; sign of i indicates a
    fast (+) or slow (−) solution.
- **Column meanings (Ch.3 tables), per §3.8 and Nomenclature:**
  - **Aphelion Ratio (AR):** (cycler aphelion) / (Mars semi-major axis). Footnote a
    on Table 3.4: "Cycler reaches Mars if Aphelion Ratio ≥ 1". (Catalogue derives
    aphelion_au ≈ AR × 1.52 AU.)
  - **Turn Ratio (TR):** required turn angle / maximum ballistic turn. TR ≥ 1 ⇒
    strictly ballistic. Footnote b on Table 3.4: "All flybys have minimum altitudes
    above 200 km if Turn Ratio ≥ 1".
  - **Earth→Mars (or aphelion) Time (days):** transit time of the E→M leg.
  - **Earth V∞ / Mars V∞ (km/s):** heliocentric hyperbolic excess speed at Earth and
    Mars in the circular-coplanar model.
  - **Required Geocentric Turning Angle at each Flyby (deg):** one value per
    gravity-assisted Earth flyby in the cycle.
- **Gate values (captions):** Tables 3.4, 3.9, 3.10, 3.11 are all filtered with
  **AR_MIN = 0.9 and TR_MIN = 0.85** (these are the "ballistic or near-ballistic"
  documentation thresholds, stated arbitrary at p.88).

---

## Table 3.4 (p.83) — Two, three, and four-synodic period ballistic or near-ballistic cyclers

Caption (verbatim): **"Table 3.4: Two, three, and four-synodic period ballistic or
near-ballistic cyclers  AR_MIN=0.9, and TR_MIN=0.85"**

Columns: `Cycler p.h.s.i` | `Aphelion Ratio` | `Turn Ratio` |
`Earth→Mars (or aphelion) Time (days)` | `Earth V∞ (km/s)` | `Mars V∞ (km/s)` |
`Required Geocentric Turning Angle at each Flyby (deg)` (variable-width list).

Footnotes (verbatim):
- a: "Cycler reaches Mars if Aphelion Ratio ≥ 1"
- b: "All flybys have minimum altitudes above 200 km if Turn Ratio ≥ 1"
- c: "Aldrin cycler^{1,13,14}" (on row 1.0.1.-1)
- d: "'Case 2' cycler described by Byrnes et al^{16}" (on row 2.1.1.+2)
- e: "'Case 3' cycler described by Byrnes et al^{9}" (on row 2.3.1.+1)

(Note: superscript footnote markers in the PDF: row 1.0.1.-1 carries marker c;
2.1.1.+2 carries d; 2.3.1.+1 carries e. The bold rows in the PDF are the cyclers
detailed further in the text/figures.)

| Cycler p.h.s.i | AR | TR | E→M time (d) | Earth V∞ | Mars V∞ | Geocentric turn angles each flyby (deg) |
|---|---|---|---|---|---|---|
| 1.0.1.-1 (c) | 1.47 | 0.86 | 146 | 6.5 | 9.7 | 84 |
| 2.1.1.+2 (d) | 0.95 | 1.11 | 207 | 4.1 | 2.0 | 92, 92 |
| 2.3.1.+1 (e) | 1.08 | 0.92 | 143 | 5.4 | 5.3 | 93, 93 |
| 2.5.1.+0 | 1.44 | 1.12 | 94 | 7.8 | 9.9 | 54, 54, 54 |
| 3.1.1.+3 | 1.07 | 1.19 | 174 | 3.6 | 4.6 | 93, 93 |
| 3.1.1.+2 | 1.43 | 0.89 | 115 | 5.4 | 9.2 | 95, 95 |
| 3.1.2.+1 | 1.07 | 1.23 | 181 | 3.4 | 4.6 | 93, 93, 24 |
| 3.1.3.+0 | 1.43 | 0.93 | 123 | 5.1 | 9.1 | 95, 95, 16, 16 |
| 3.3.1.+2 | 1.19 | 1.06 | 141 | 4.3 | 6.8 | 94, 94 |
| 3.5.1.+2 | 0.94 | 1.80 | 231 | 2.7 | 1.5 | 70, 70, 70, 70 |
| 3.5.1.+1 | 1.43 | 1.15 | 115 | 5.4 | 9.2 | 73, 73, 73, 73 |
| 3.5.2.+0 | 1.43 | 1.06 | 121 | 5.2 | 9.2 | 83, 83, 83, 83, 24 |
| 3.7.1.+1 | 1.07 | 1.56 | 175 | 3.6 | 4.6 | 71, 71, 71, 71 |
| 3.9.1.+0 | 1.43 | 1.17 | 116 | 5.4 | 9.2 | 72, 45, 45, 45, 45, 72 |
| 4.0.3.+1 | 1.07 | 1.18 | 160 | 4.3 | 4.9 | 85, 85, 85 |
| 4.1.1.-6 | 0.94 | 1.37 | 256 | 2.7 | 1.6 | 92, 92 |
| 4.1.1.-5 | 1.15 | 1.11 | 173 | 4.1 | 6.1 | 94, 94 |
| 4.1.1.-4 | 1.44 | 0.89 | 137 | 5.5 | 9.3 | 95, 95 |
| 4.1.2.-3 | 0.94 | 1.40 | 250 | 2.6 | 1.5 | 92, 92, 24 |
| 4.1.2.-2 | 1.43 | 0.93 | 132 | 5.2 | 9.2 | 95, 95, 24 |
| 4.1.4.-1 | 1.43 | 0.93 | 129 | 5.1 | 9.2 | 95, 95, 24, 12, 12, 12 |
| 4.3.1.-5 | 0.99 | 1.29 | 268 | 3.1 | 2.5 | 93, 93 |
| 4.3.1.-4 | 1.26 | 1.01 | 154 | 4.7 | 7.6 | 94, 94 |
| 4.5.1.-4 | 1.07 | 1.55 | 196 | 3.6 | 4.7 | 71, 71, 71, 71 |
| 4.5.1.-3 | 1.44 | 1.15 | 137 | 5.5 | 9.3 | 73, 73, 73, 73 |
| 4.5.2.-2 | 1.07 | 1.40 | 191 | 3.4 | 4.6 | 81, 81, 81, 81, 24 |
| 4.5.3.-1 | 1.43 | 1.02 | 130 | 5.1 | 9.2 | 87, 87, 87, 87, 16, 16 |
| 4.6.1.-4 | 0.91 | 1.50 | 154 | 6.8 | 2.1 | 46, 46, 46, 46 |
| 4.6.3.+0 | 1.43 | 0.88 | 105 | 6.4 | 9.5 | 60, 60, 60, 60, 84, 84 |
| 4.7.1.-3 | 1.20 | 1.38 | 163 | 4.3 | 6.8 | 72, 72, 72, 72 |
| 4.7.1.-2 | 1.77 | 0.96 | 120 | 6.6 | 11.4 | 74, 74, 74, 74 |
| 4.8.1.+3 | 0.96 | 1.64 | 164 | 7.7 | 3.1 | 37, 37, 37, 37, 37 |
| 4.8.1.+2 | 1.31 | 0.86 | 76 | 12.5 | 10.7 | 37, 37, 37, 37, 37 |
| 4.9.1.-3 | 0.94 | 1.83 | 256 | 2.7 | 1.6 | 69, 45, 45, 45, 45, 69 |
| 4.9.1.-2 | 1.44 | 1.16 | 137 | 5.5 | 9.3 | 72, 45, 45, 45, 45, 72 |
| 4.9.2.-1 | 1.44 | 1.05 | 132 | 5.2 | 9.2 | 83, 45, 45, 45, 45, 83, 24 |
| 4.10.1.-3 | 0.92 | 1.46 | 263 | 10.2 | 3.6 | 30, 30, 30, 30, 30, 30 |
| 4.10.1.+2 | 1.03 | 1.65 | 131 | 8.9 | 5.0 | 31, 31, 31, 31, 31, 31 |
| 4.11.1.-2 | 1.07 | 1.58 | 195 | 3.6 | 4.7 | 70, 45, 45, 45, 45, 70 |
| 4.12.1.-2 | 0.97 | 1.43 | 268 | 11.6 | 4.8 | 25, 25, 25, 25, 25, 25 |
| 4.12.1.+1 | 1.16 | 1.48 | 93 | 10.8 | 8.2 | 27, 27, 27, 27, 27, 27 |
| 4.13.1.-1 | 1.44 | 1.16 | 137 | 5.5 | 9.3 | 72, 30, 30, 30, 30, 30, 30, 72 |
| 4.14.1.-1 | 1.12 | 1.13 | 199 | 14.7 | 9.4 | 22, 22, 22, 22, 22, 22, 22, 22 |
| 4.14.1.+0 | 1.49 | 1.09 | 66 | 14.1 | 12.7 | 25, 25, 25, 25, 25, 25, 25, 25 |

**Row count: 44 cyclers.** (1 two-synodic-or-less precursor row 1.0.1.-1 = Aldrin;
the 2./3./4. families follow.)

---

## Tables 3.5–3.8 (p.84) — Detailed flyby-maneuver state tables for four selected cyclers

Intro text (p.84, verbatim): "Details about the flyby maneuvers for a few of the
discussed solutions are provided in Table 3.5 - Table 3.8. They have sufficient data
to simulate one complete cycle plus the first leg of the second cycle for each
described cycler. The ecliptic is the x-y axis plane and the initial position of the
Earth is always on the x-axis."

These are **state / maneuver tables** (NOT member lists): each column is an encounter
(Location = Earth or Mars), rows are `time (days)`, `Δvx`, `Δvy`, `Δvz` (km/s) at
that encounter, plus an initial Mars position vector `r_mars` (AU). Footnote a marks
the initial column's Δv as "initial v∞ with respect to Earth".

### Table 3.5: Cycler 2.5.1.+0
Footnote a: "initial v∞ with respect to Earth".

| | Earth | Mars | Earth | Earth | Earth | Earth | Mars |
|---|---|---|---|---|---|---|---|
| time (days) | 0 | 94 | 652 | 1018 | 1200 | 1565 | 1659 |
| Δvx (km/s) | 6.50 (a) | 0 | -5.19 | 1.40 | -1.40 | -5.29 | 0 |
| Δvy (km/s) | 4.35 (a) | 0 | -1.41 | -6.12 | 6.12 | -0.98 | 0 |
| Δvz (km/s) | 0 (a) | 0 | 4.55 | 3.20 | 3.20 | 4.55 | 0 |

`r_mars` at t0 = [1.41  0.57  0] AU

### Table 3.6: Cycler 3.1.2.+1
Footnote a: "initial v∞ with respect to Earth".

| | Earth | Mars | Earth | Earth | Earth | Mars |
|---|---|---|---|---|---|---|
| time (days) | 0 | 181 | 1083 | 1265 | 2348 | 2529 |
| Δvx (km/s) | 0.71 (a) | 0 | -0.09 | -1.48 | -1.28 | 0 |
| Δvy (km/s) | 3.32 (a) | 0 | -3.59 | -3.27 | 0.62 | 0 |
| Δvz (km/s) | 0 (a) | 0 | 3.39 | 3.39 | 0 | 0 |

`r_mars` at t0 = [1.15  0.99  0] AU

### Table 3.7: Cycler 4.3.1.-5
Footnote a: "initial v∞ with respect to Earth".
Footnote b: "0.008 AU from Mars (cycler aphelion)" (on the "Mars" column header
marked b).

| | Earth | Mars (b) | Earth | Earth | Mars (b) |
|---|---|---|---|---|---|
| time (days) | 0 | 268 | 2583 | 3131 | 3399 |
| Δvx (km/s) | -1.24 (a) | 0 | 0.18 | 2.42 | 0 |
| Δvy (km/s) | 2.84 (a) | 0 | -3.24 | -2.16 | 0 |
| Δvz (km/s) | 0 (a) | 0 | 3.09 | 3.09 | 0 |

`r_mars` at t0 = [0.93  1.20  0] AU

### Table 3.8: Cycler 4.5.2.-2
Footnote a: "initial v∞ with respect to Earth".

| | Earth | Mars | Earth | Earth | Earth | Earth | Mars |
|---|---|---|---|---|---|---|---|
| time (days) | 0 | 191 | 1109 | 1474 | 1657 | 2022 | 3131 | 3322 |
| Δvx (km/s) | -0.71 (a) | 0 | 3.38 | -3.29 | 3.29 | -1.80 | 1.29 | 0 |
| Δvy (km/s) | 3.34 (a) | 0 | -2.86 | -0.75 | 0.75 | -4.04 | 0.62 | 0 |
| Δvz (km/s) | 0 (a) | 0 | -0.50 | -2.91 | -2.91 | -0.50 | 0 | 0 |

`r_mars` at t0 = [1.03  1.12  0] AU

(Note Table 3.8 has 8 encounter columns: Earth, Mars, Earth, Earth, Earth, Earth,
Earth, Mars. The header in the PDF reads "Earth Mars Earth Earth Earth Earth Earth
Mars"; the body rows have 8 values matching that.)

---

## Table 3.9 (p.90) — Five-synodic period ballistic (or near) cyclers

Caption (verbatim): **"Table 3.9: Five-synodic period ballistic (or near) cyclers
AR_MIN=0.9, and TR_MIN=0.85"**

Columns: `Cycler p.h.s.i` | `Aphelion Ratio` | `Turn Ratio` |
`Earth→Mars or aphelion (days)` | `v∞ at Earth (km/s)` | `v∞ at Mars (km/s)` |
`Required Geocentric Turning Angle at each Flyby (deg)`.

| Cycler p.h.s.i | AR | TR | E→M/aph (d) | V∞ Earth | V∞ Mars | Geocentric turn angles (deg) |
|---|---|---|---|---|---|---|
| 5.1.1.-7 | 1.04 | 0.97 | 229 | 5.0 | 4.3 | 93, 93 |
| 5.1.2.-3 | 1.20 | 1.00 | 168 | 4.7 | 7.0 | 94, 94, 67 |
| 5.1.5.-1 | 1.44 | 0.92 | 133 | 5.2 | 9.2 | 95, 95, 28, 28, 28, 28 |
| 5.2.1.-7 | 0.90 | 1.07 | 182 | 4.5 | 1.3 | 92, 92 |
| 5.2.2.+2 | 1.20 | 0.94 | 128 | 5.2 | 7.1 | 94, 94, 85 |
| 5.2.5.+0 | 1.43 | 0.91 | 118 | 5.3 | 9.2 | 95, 95, 37, 37, 37, 37 |
| 5.3.1.-7 | 0.92 | 1.17 | 270 | 3.8 | 1.4 | 92, 92 |
| 5.3.1.-6 | 1.10 | 0.90 | 205 | 5.5 | 5.7 | 93, 93 |
| 5.3.3.-2 | 1.07 | 1.19 | 195 | 3.6 | 4.7 | 93, 93, 47, 47 |
| 5.4.1.+6 | 0.94 | 1.45 | 189 | 4.9 | 1.9 | 63, 63, 63 |
| 5.4.1.+5 | 1.12 | 1.06 | 122 | 7.0 | 6.3 | 64, 64, 64 |
| 5.4.3.+1 | 1.07 | 1.45 | 170 | 3.8 | 4.7 | 75, 75, 75, 61, 61 |
| 5.5.1.-6 | 0.96 | 1.95 | 279 | 4.3 | 2.2 | 52, 52, 52, 52 |
| 5.5.1.-5 | 1.18 | 1.44 | 186 | 6.2 | 7.0 | 53, 53, 53, 53 |
| 5.5.1.-4 | 1.48 | 1.08 | 154 | 8.0 | 10.3 | 54, 54, 54, 54 |
| 5.5.2.-3 | 0.94 | 1.79 | 262 | 3.0 | 1.7 | 63, 63, 63, 63, 68 |
| 5.5.2.-2 | 1.45 | 1.19 | 142 | 5.9 | 9.5 | 66, 66, 66, 66, 66 |
| 5.5.4.-1 | 1.44 | 1.10 | 134 | 5.3 | 9.3 | 78, 78, 78, 78, 35, 35, 35 |
| 5.6.1.+5 | 0.98 | 1.74 | 198 | 5.4 | 2.7 | 49, 49, 49, 49 |
| 5.6.1.+4 | 1.20 | 1.23 | 107 | 7.7 | 7.6 | 49, 49, 49, 49 |
| 5.6.1.+3 | 1.49 | 0.90 | 82 | 9.8 | 11.0 | 50, 50, 50, 50 |
| 5.6.2.+2 | 0.94 | 1.36 | 219 | 3.3 | 1.7 | 58, 58, 58, 58, 86 |
| 5.6.2.+1 | 1.44 | 0.87 | 104 | 6.4 | 9.5 | 60, 60, 60, 60, 84 |
| 5.6.4.+0 | 1.43 | 1.16 | 116 | 5.4 | 9.2 | 73, 73, 73, 73, 46, 46, 46 |
| 5.7.1.-5 | 1.02 | 1.76 | 245 | 4.8 | 3.6 | 52, 52, 52, 52 |
| 5.7.1.-4 | 1.30 | 1.27 | 169 | 7.0 | 8.5 | 53, 53, 53, 53 |
| 5.7.1.-3 | 1.71 | 0.93 | 142 | 9.1 | 11.9 | 54, 54, 54, 54 |
| 5.8.1.+4 | 1.03 | 1.91 | 154 | 6.1 | 4.3 | 40, 40, 40, 40, 40 |
| 5.8.1.+3 | 1.31 | 1.30 | 94 | 8.6 | 9.1 | 41, 41, 41, 41, 41 |
| 5.8.1.+2 | 1.72 | 0.92 | 73 | 11.0 | 12.6 | 42, 42, 42, 42, 42 |
| 5.9.1.-4 | 1.10 | 2.09 | 204 | 5.6 | 5.7 | 40, 40, 40, 40, 40, 40 |
| 5.9.1.-3 | 1.48 | 1.38 | 154 | 8.0 | 10.3 | 42, 42, 42, 42, 42, 42 |
| 5.9.1.-2 | 2.15 | 0.93 | 130 | 10.5 | 13.8 | 45, 44, 44, 44, 44, 45 |
| 5.9.2.-2 | 1.08 | 1.57 | 198 | 4.0 | 4.9 | 60, 45, 45, 45, 45, 60, 67 |
| 5.9.2.+1 | 2.10 | 0.90 | 117 | 7.9 | 12.9 | 65, 45, 45, 45, 45, 65, 65 |
| 5.9.3.-1 | 1.44 | 1.16 | 137 | 5.5 | 9.3 | 72, 45, 45, 45, 45, 72, 46, 46 |
| 5.10.1.+3 | 1.11 | 1.94 | 123 | 6.9 | 6.2 | 35, 35, 35, 35, 35 |
| 5.10.1.+2 | 1.49 | 1.25 | 82 | 9.8 | 10.9 | 36, 36, 36, 36, 36 |
| 5.10.2.+1 | 1.07 | 1.71 | 160 | 4.3 | 4.9 | 59, 59, 59, 59, 59, 59 |
| 5.10.2.+0 | 2.08 | 0.87 | 81 | 8.6 | 13.0 | 62, 62, 62, 62, 62, 62 |
| 5.10.3.+0 | 1.43 | 1.24 | 112 | 5.7 | 9.3 | 66, 45, 45, 45, 45, 66, 60, 60 |
| 5.11.1.-3 | 1.24 | 1.76 | 177 | 6.6 | 7.8 | 41, 41, 41, 41, 41, 41 |
| 5.11.1.-2 | 1.83 | 1.09 | 138 | 9.5 | 12.3 | 43, 43, 43, 43, 43, 43 |
| 5.12.1.+2 | 1.14 | 1.00 | 101 | 9.6 | 7.5 | 47, 47, 47, 47, 47, 47, 47, 47 |
| 5.12.1.+1 | 1.24 | 1.81 | 101 | 8.1 | 8.3 | 32, 32, 32, 32, 32, 32 |
| 5.12.1.+1 | 1.82 | 1.07 | 70 | 11.5 | 13.2 | 34, 34, 34, 34, 34, 34 |
| 5.13.1.-3 | 0.97 | 2.69 | 280 | 4.3 | 2.3 | 37, 30, 30, 30, 30, 30, 37 |
| 5.13.1.-2 | 1.49 | 1.39 | 153 | 8.1 | 10.3 | 30, 30, 30, 30, 30, 42 |
| 5.13.2.-1 | 1.45 | 1.20 | 141 | 5.9 | 9.5 | 66, 66, 66, 66, 66, 66, 66 |
| 5.14.1.+2 | 0.97 | 3.06 | 196 | 5.3 | 2.6 | 28, 28, 28, 28, 28, 28, 28 |
| 5.14.1.+1 | 1.48 | 1.41 | 82 | 9.8 | 10.9 | 32, 30, 30, 30, 30, 30, 32 |
| 5.14.2.+0 | 1.43 | 1.22 | 105 | 6.4 | 9.7 | 60, 60, 60, 60, 60, 60, 60 |
| 5.15.1.-2 | 1.11 | 2.13 | 202 | 5.6 | 5.8 | 39, 30, 30, 30, 30, 30, 39 |
| 5.15.1.-1 | 2.16 | 0.93 | 130 | 10.5 | 13.9 | 45, 30, 30, 30, 30, 30, 45 |
| 5.16.1.+1 | 1.10 | 2.37 | 126 | 6.8 | 6.0 | 26, 26, 26, 26, 26, 26, 26, 29 |
| 5.16.1.+0 | 2.12 | 0.91 | 64 | 12.6 | 14.6 | 35, 25, 25, 25, 25, 25, 25, 35 |
| 5.17.1.-1 | 1.50 | 1.37 | 152 | 8.1 | 10.4 | 22, 22, 22, 22, 22, 22, 22, 22, 22, 42 |
| 5.18.1.+0 | 1.46 | 1.45 | 84 | 9.6 | 10.6 | 32, 22, 22, 22, 22, 22, 22, 22, 22, 32 |

**Row count: 58 cyclers.** (Two rows carry the label "5.12.1.+1" in the PDF — both
transcribed as printed; this is likely a typeset duplication of the index, the two
rows differ in AR/TR/V∞ — preserved verbatim, do NOT silently fix.)

---

## Table 3.10 (p.91) — Six-synodic period ballistic (or near) cyclers. Part I

Caption (verbatim): **"Table 3.10: Six-synodic period ballistic (or near) cyclers.
Part I  AR_MIN=0.9, and TR_MIN=0.85"**

Columns: `Cycler p- h-s -i` | `Aphelion Ratio` | `Turn Ratio` |
`Earth→Mars (or aphelion) Time (days)` | `v∞ at Earth (km/s)` |
`v∞ at Mars (km/s)` | `Required Geocentric Turning Angle at each Flyby (deg)`.

Footnotes (verbatim, on the four 6.0.1 rows):
- a: '"Cycler 6S9"' (on 6.0.1.+9)
- b: '"Cycler 6S8"' (on 6.0.1.+8)
- c: '"Cycler 6S7"' (on 6.0.1.+7)
- d: '"Cycler 6S6" described described by McConaghy^{17} et al' (on 6.0.1.+6)

| Cycler p-h-s-i | AR | TR | E→M/aph (d) | V∞ Earth | V∞ Mars | Geocentric turn angles (deg) |
|---|---|---|---|---|---|---|
| 6.0.1.+9 (a) | 0.92 | 1.40 | 213 | 3.0 | 1.2 | 86 |
| 6.0.1.+8 (b) | 1.03 | 1.22 | 179 | 4.0 | 3.9 | 85 |
| 6.0.1.+7 (c) | 1.17 | 1.07 | 133 | 5.0 | 6.7 | 85 |
| 6.0.1.+6 (d) | 1.34 | 0.93 | 111 | 6.0 | 8.7 | 84 |
| 6.1.2.-4 | 1.09 | 0.91 | 203 | 4.9 | 5.4 | 93, 93, 100 |
| 6.1.3.-3 | 0.95 | 1.30 | 264 | 3.1 | 1.7 | 92, 92, 74, 74 |
| 6.1.4.-2 | 1.07 | 1.16 | 197 | 3.8 | 4.8 | 93, 93, 57, 57, 57 |
| 6.1.6.-1 | 1.44 | 0.90 | 135 | 5.4 | 9.3 | 95, 95, 39, 39, 39, 39, 39 |
| 6.2.1.+8 | 0.94 | 1.26 | 220 | 3.3 | 1.7 | 92, 92 |
| 6.2.1.+7 | 1.08 | 1.07 | 158 | 4.3 | 5.0 | 93, 93 |
| 6.2.1.+6 | 1.24 | 0.91 | 123 | 5.4 | 7.6 | 94, 94 |
| 6.2.2.+3 | 1.07 | 1.19 | 174 | 3.6 | 4.6 | 93, 93, 47 |
| 6.2.2.+2 | 1.43 | 0.89 | 115 | 5.4 | 9.2 | 95, 95, 46 |
| 6.2.3.+2 | 0.94 | 1.39 | 235 | 2.6 | 1.5 | 92, 92, 32, 32 |
| 6.2.3.+1 | 1.43 | 0.92 | 119 | 5.2 | 9.2 | 95, 95, 31, 31 |
| 6.2.4.+1 | 1.07 | 1.23 | 181 | 3.4 | 4.6 | 93, 93, 24, 24, 24 |
| 6.2.6.+0 | 1.43 | 0.93 | 123 | 5.1 | 9.1 | 95, 95, 16, 16, 16, 16, 16 |
| 6.3.1.-9 | 0.92 | 0.89 | 279 | 5.7 | 1.8 | 91, 91 |
| 6.3.4.+1 | 1.07 | 1.04 | 156 | 4.5 | 5.0 | 93, 93, 93, 93, 93 |
| 6.4.1.+7 | 0.98 | 1.59 | 227 | 3.6 | 2.4 | 70, 70, 70 |
| 6.4.1.+6 | 1.13 | 1.33 | 142 | 4.7 | 6.0 | 71, 71, 71 |
| 6.4.1.+5 | 1.33 | 1.11 | 113 | 5.9 | 8.5 | 71, 71, 71 |
| 6.4.1.+4 | 1.58 | 0.93 | 96 | 7.0 | 10.6 | 72, 72, 72 |
| 6.5.1.-8 | 0.95 | 1.62 | 283 | 6.2 | 2.4 | 47, 47, 47, 47 |
| 6.5.1.-7 | 1.11 | 1.16 | 213 | 8.4 | 6.6 | 47, 47, 47, 47 |
| 6.5.1.-6 | 1.31 | 0.88 | 180 | 10.4 | 9.8 | 48, 48, 48, 48 |
| 6.5.5.-1 | 1.44 | 1.15 | 137 | 5.5 | 9.3 | 73, 73, 73, 73, 46, 46, 46, 46 |
| 6.6.1.-6 | 1.02 | 1.82 | 189 | 3.9 | 3.5 | 58, 58, 58, 58 |
| 6.6.1.+5 | 1.20 | 1.48 | 128 | 5.2 | 7.1 | 59, 59, 59, 59 |
| 6.6.1.+4 | 1.45 | 1.21 | 104 | 6.4 | 9.6 | 60, 60, 60, 60 |
| 6.6.1.+3 | 1.78 | 0.99 | 89 | 7.7 | 11.7 | 61, 61, 61, 61 |
| 6.6.2.+2 | 1.19 | 1.38 | 141 | 4.3 | 6.8 | 72, 72, 72, 72, 46 |
| 6.6.2.+1 | 1.77 | 0.96 | 99 | 6.6 | 11.3 | 74, 74, 74, 74, 45 |
| 6.6.5.+0 | 1.43 | 1.03 | 122 | 5.1 | 9.2 | 85, 85, 85, 85, 19, 19, 19, 19 |
| 6.7.1.-7 | 0.98 | 1.49 | 289 | 6.7 | 3.1 | 47, 47, 47, 47 |
| 6.7.1.-6 | 1.17 | 1.05 | 199 | 9.1 | 7.8 | 47, 47, 47, 47 |
| 6.7.2.+3 | 0.91 | 0.98 | 176 | 5.1 | 1.5 | 62, 62, 62, 91, 91 |
| 6.7.3.-2 | 1.08 | 1.40 | 199 | 4.1 | 5.0 | 62, 62, 62, 62, 73, 73 |
| 6.7.5.+0 | 1.43 | 0.99 | 107 | 6.1 | 9.4 | 62, 62, 62, 62, 77, 77, 77, 77 |
| 6.8.1.+6 | 0.91 | 2.39 | 211 | 2.9 | 1.0 | 51, 51, 51, 51, 51 |
| 6.8.1.+5 | 1.08 | 1.89 | 158 | 4.3 | 5.0 | 53, 53, 53, 53, 53 |
| 6.8.1.+4 | 1.30 | 1.50 | 116 | 5.7 | 8.3 | 54, 54, 54, 54, 54 |
| 6.8.1.+3 | 1.62 | 1.18 | 95 | 7.2 | 10.8 | 55, 55, 55, 55, 55 |
| 6.8.1.+2 | 2.09 | 0.94 | 81 | 8.6 | 13.1 | 57, 57, 57, 57, 57 |
| 6.8.3.+1 | 1.07 | 1.47 | 179 | 3.4 | 4.6 | 77, 60, 60, 60, 77, 32, 32 |
| 6.9.1.-6 | 1.03 | 1.97 | 248 | 7.3 | 4.5 | 32, 32, 32, 32, 32, 32 |
| 6.9.1.-5 | 1.26 | 1.35 | 186 | 10.0 | 9.1 | 33, 33, 33, 33, 33, 33 |
| 6.9.1.-4 | 1.57 | 0.98 | 161 | 12.4 | 12.4 | 34, 34, 34, 34, 34, 34 |
| 6.9.2.-3 | 0.96 | 1.61 | 274 | 3.8 | 2.0 | 54, 54, 54, 54, 64, 67, 67 |
| 6.9.2.-2 | 1.47 | 0.94 | 150 | 7.2 | 10.0 | 56, 56, 56, 56, 69, 69, 69 |
| 6.9.4.-1 | 1.45 | 1.21 | 139 | 5.7 | 9.4 | 45, 45, 45, 45, 67, 56, 56, 56 |

**Row count: 51 cyclers.**

---

## Table 3.11 (p.92) — Six-synodic period ballistic (or near) cyclers. Part II

Caption (verbatim): **"Table 3.11: Six-synodic period ballistic (or near) cyclers.
Part II  AR_MIN=0.9, and TR_MIN=0.85"**

Columns: identical to Table 3.10 (`Cycler p- h-s -i` | `Aphelion Ratio` |
`Turn Ratio` | `Earth→Mars (or aphelion) Time (days)` | `v∞ at Earth (km/s)` |
`v∞ at Mars (km/s)` | `Required Geocentric Turning Angle at each Flyby (deg)`).

| Cycler p-h-s-i | AR | TR | E→M/aph (d) | V∞ Earth | V∞ Mars | Geocentric turn angles (deg) |
|---|---|---|---|---|---|---|
| 6.10.1.+5 | 0.94 | 2.31 | 219 | 3.3 | 1.7 | 50, 45, 45, 45, 45, 50 |
| 6.10.1.+4 | 1.15 | 1.76 | 137 | 4.9 | 6.4 | 52, 45, 45, 45, 45, 52 |
| 6.10.1.+3 | 1.44 | 1.34 | 104 | 6.4 | 9.6 | 54, 45, 45, 45, 45, 54 |
| 6.10.1.+2 | 1.89 | 1.02 | 86 | 8.1 | 12.2 | 56, 45, 45, 45, 45, 56 |
| 6.10.2.+2 | 0.94 | 1.84 | 231 | 2.7 | 1.5 | 69, 45, 45, 45, 45, 69, 47 |
| 6.10.2.+1 | 1.43 | 1.17 | 115 | 5.4 | 9.2 | 72, 45, 45, 45, 45, 72, 46 |
| 6.10.4.+0 | 1.43 | 1.06 | 121 | 5.2 | 9.2 | 83, 45, 45, 45, 45, 83, 24, 24, 24 |
| 6.11.1.-5 | 1.09 | 1.75 | 219 | 8.2 | 6.2 | 33, 33, 33, 33, 33, 33 |
| 6.11.1.-4 | 1.38 | 1.17 | 173 | 11.0 | 10.7 | 33, 33, 33, 33, 33, 33 |
| 6.11.2.+2 | 0.98 | 1.61 | 191 | 6.1 | 2.9 | 48, 48, 48, 48, 48, 48, 48, 48 |
| 6.12.1.+4 | 1.00 | 2.14 | 232 | 3.7 | 2.7 | 51, 36, 36, 36, 36, 36, 36, 51 |
| 6.12.1.+3 | 1.26 | 1.57 | 120 | 5.5 | 7.8 | 53, 36, 36, 36, 36, 36, 36, 53 |
| 6.12.1.+2 | 1.67 | 1.15 | 92 | 7.4 | 11.1 | 55, 36, 36, 36, 36, 36, 36, 55 |
| 6.13.1.+5 | 1.08 | 0.94 | 79 | 16.7 | 9.7 | 22, 22, 22, 22, 22, 22, 22, 22 |
| 6.13.1.-5 | 0.90 | 3.44 | 276 | 5.4 | 1.5 | 25, 25, 25, 25, 25, 25, 25, 25 |
| 6.13.1.-4 | 1.18 | 1.92 | 198 | 9.2 | 7.9 | 26, 26, 26, 26, 26, 26, 26, 26 |
| 6.13.1.-3 | 1.58 | 1.21 | 160 | 12.4 | 12.5 | 27, 27, 27, 27, 27, 27, 27, 27 |
| 6.13.2.-2 | 1.10 | 1.65 | 202 | 5.0 | 5.4 | 55, 55, 55, 55, 55, 55, 55, 55 |
| 6.13.3.-1 | 1.46 | 1.06 | 143 | 6.1 | 9.6 | 60, 30, 30, 30, 30, 30, 60, 72, 72 |
| 6.14.1.+3 | 1.08 | 1.93 | 158 | 4.3 | 4.9 | 52, 30, 30, 30, 30, 30, 52 |
| 6.14.1.+2 | 1.44 | 1.34 | 104 | 6.4 | 9.6 | 54, 30, 30, 30, 30, 30, 54 |
| 6.14.1.+1 | 2.09 | 0.94 | 81 | 8.6 | 13.1 | 57, 30, 30, 30, 30, 30, 57 |
| 6.14.2.+1 | 1.07 | 1.59 | 175 | 3.6 | 4.6 | 70, 30, 30, 30, 30, 30, 70, 47 |
| 6.14.2.+0 | 2.08 | 0.85 | 91 | 7.4 | 12.7 | 75, 30, 30, 30, 30, 30, 75, 45 |
| 6.14.3.+0 | 1.43 | 1.10 | 119 | 5.2 | 9.2 | 79, 30, 30, 30, 30, 30, 79, 31, 31 |
| 6.15.1.+4 | 1.11 | 1.90 | 74 | 17.2 | 10.4 | 22, 22, 22, 22, 22, 22, 22 |
| 6.15.1.-4 | 0.95 | 3.00 | 285 | 6.3 | 2.6 | 25, 25, 25, 25, 25, 25, 25, 25 |
| 6.15.1.-3 | 1.32 | 1.57 | 179 | 10.5 | 10.0 | 26, 26, 26, 26, 26, 26, 26, 26 |
| 6.15.1.-2 | 1.94 | 0.94 | 147 | 14.3 | 14.8 | 28, 28, 28, 28, 28, 28, 28, 28 |
| 6.15.2.+1 | 1.11 | 1.27 | 117 | 7.7 | 6.4 | 39, 39, 39, 39, 48, 48, 48, 48 |
| 6.16.1.+2 | 1.20 | 1.67 | 128 | 5.2 | 7.1 | 53, 26, 26, 26, 26, 26, 26, 53 |
| 6.16.1.+1 | 1.78 | 1.08 | 89 | 7.7 | 11.7 | 56, 25, 25, 25, 25, 25, 25, 56 |
| 6.17.1.+3 | 1.14 | 1.06 | 69 | 17.9 | 11.3 | 17, 17, 17, 17, 17, 17, 17, 17, 17, 17 |
| 6.17.1.-3 | 1.04 | 2.91 | 241 | 7.5 | 4.8 | 21, 21, 21, 21, 21, 21, 21, 21, 21, 21 |
| 6.17.1.-2 | 1.59 | 1.24 | 160 | 12.5 | 12.6 | 26, 22, 22, 22, 22, 22, 22, 22, 22, 26 |
| 6.17.2.-1 | 1.48 | 1.30 | 149 | 7.3 | 10.1 | 48, 45, 45, 45, 48, 50, 50, 50, 50, 50 |
| 6.18.1.+2 | 0.94 | 2.32 | 219 | 3.3 | 1.7 | 50, 22, 22, 22, 22, 22, 22, 22, 22, 50 |
| 6.18.1.+1 | 1.44 | 1.34 | 104 | 6.4 | 9.5 | 54, 22, 22, 22, 22, 22, 22, 22, 22, 54 |
| 6.18.2.+0 | 1.43 | 1.17 | 116 | 5.4 | 9.2 | 72, 22, 22, 22, 22, 22, 22, 22, 22, 72, 46 |
| 6.19.1.+2 | 1.19 | 0.97 | 64 | 18.8 | 12.5 | 17, 17, 17, 17, 17, 17, 17, 17 |
| 6.19.1.-2 | 1.20 | 2.09 | 195 | 9.4 | 8.2 | 23, 22, 22, 22, 22, 22, 22, 23 |
| 6.19.2.+0 | 1.47 | 1.18 | 79 | 10.6 | 11.1 | 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35 |
| 6.20.1.-4 | 0.93 | 1.09 | 183 | 12.8 | 5.0 | 29, 20, 20, 20, 20, 20, 20, 20, 20, 29 |
| 6.20.1.+1 | 1.07 | 1.94 | 160 | 4.3 | 4.9 | 52, 20, 20, 20, 20, 20, 20, 20, 20, 52 |
| 6.20.1.+0 | 2.08 | 0.94 | 81 | 8.6 | 13.0 | 57, 20, 20, 20, 20, 20, 20, 20, 20, 57 |
| 6.21.1.+1 | 1.29 | 1.01 | 57 | 20.3 | 14.4 | 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 |
| 6.21.1.-1 | 1.63 | 1.19 | 158 | 12.7 | 12.9 | 27, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 27 |
| 6.22.1.+0 | 1.43 | 1.35 | 105 | 6.4 | 9.5 | 54, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 54 |

**Row count: 48 cyclers.** Tables 3.10 + 3.11 together = 99 six-synodic cyclers
(matches the catalogue's 99 `russell-ocampo-6.*` rows).

Cross-check vs the prose (p.88): "If Table 3.4 were expanded to include five- and
six-synodic period cyclers, an additional 157 entries would be added. These are
documented in Table 3.9 - Table 3.11." → 58 (T3.9) + 51 (T3.10) + 48 (T3.11) = **157**.
Confirmed consistent with the transcription above.

---

## Appendix C (pp.201–245) — Cycler Trajectories Using an Ephemeris Model

### Model / frame / provenance flags

- **Model / fidelity: EPHEMERIS / DE405-fit (accurate model).** Header (p.201):
  "This appendix documents the full trajectories that are summarized in Table 5.5…
  optimized solutions were obtained for twenty-one launch windows for each of the
  203 cyclers, totaling 4263 trajectories… The 77 presented here represent the
  lowest Δv solutions of the twenty-one launch date periods for a selected set of
  the parent cyclers."
- **Center / frame (p.202, verbatim):** "The DE405 ephemeris gives positions and
  velocities referenced to the Earth's mean equator and equinox of J2000. These
  positions and velocities are then rotated by the mean obliquity angle of
  **0.409092629205 radians** as given on page 289 of Ref. 5. Thus, all velocities
  are presented relative to the ecliptic plane." Planet velocities from **DE405**.
- **Epoch (verbatim):** "a value of i, indicates that the search was initiated with
  an epoch date between i and i+1 synodic periods following J2000." Each solution
  prints an explicit "EPOCH TIME (days after J2000)". Δv given in **m/s**.
- **Column meanings (verbatim, p.202):**
  - `E/M` = the planet where the i-th leg is initiated.
  - `time start` = time (relative to epoch) of the beginning of the leg; planet
    velocity at this time from DE405; spacecraft velocity = planet velocity + given v∞.
  - Trajectory propagated by Kepler's equation until `time dv`, when given Δv applied,
    then propagated to start of next leg / next planet encounter.
  - "The 'Epoch time' and the detailed leg information are the only values required
    to reproduce the solution. The total Δv … is an artifact of the trajectory."

### Per-solution block format (verbatim header layout)

Each solution begins with a separator line and header:
```
=============== PARENT CYCLER <shorthand> ===============
Parent cycler number                           <1-203>
Approximate search space (synodic periods after J2000)   <i>
Number of steps to walk eccentricity/inclination         <n / m>
Number of cycles                                          7
Total delta v over <NN.NN> years (km/s)                  <value>
---------EARTH TO MARS TRANSIT LEG CHARACTERISTICS BELOW-------------
LEG  E-M transit  Earth vinfm  Earth vinfp  Mars vinfm  Mars vinfp
     time (days)   (km/s)       (km/s)       (km/s)      (km/s)
 ...
AVERAGE  <...>
---------DATA NECESSARY TO REPRODUCE CYCLER BELOW-------------
EPOCH TIME (days after J2000)                  <value>
LEG E/M  time start  vinfx  vinfy  vinfz  time dv  dvx  dvy  dvz
         (days)      (km/s) (km/s) (km/s)  (days)  (km/s)(km/s)(km/s)
 ...
```
("Number of cycles" = 7 for all solutions: cyclers propagated 7 cycles, consistent
with Ch.5 Table 5.1. The `vinfm`/`vinfp` = v∞ minus/plus = before/after the flyby.)

### Worked example — first solution, VERBATIM (p.203)

Parent Cycler **6.399G1**, parent cycler number **1**, approx search space **243/243**,
steps to walk ecc/incl **7**, number of cycles **7**, total Δv over **14.56 years**
= **0.000000** km/s. This is the **ballistic Aldrin cycler** (Fig.5.7 "Ballistic
Aldrin cycler 6.399G1(#1), launch: Aug 6, 2003"). Earth→Mars transit leg
characteristics (LEG / E-M transit time (days) / Earth vinfm / Earth vinfp / Mars
vinfm / Mars vinfp, km/s):

| LEG | time (d) | Earth vinfm | Earth vinfp | Mars vinfm | Mars vinfp |
|---|---|---|---|---|---|
| 1 | 67.5 | ****** | 10.246 | 8.479 | 8.479 |
| 3 | 140.8 | 6.225 | 6.225 | 7.278 | 7.278 |
| 5 | 167.1 | 6.226 | 6.226 | 7.431 | 7.431 |
| 7 | 169.9 | 5.329 | 5.329 | 8.663 | 8.663 |
| 9 | 163.0 | 5.824 | 5.824 | 10.149 | 10.149 |
| 11 | 150.5 | 5.881 | 5.881 | 11.494 | 11.494 |
| 13 | 139.4 | 6.044 | 6.044 | 11.806 | 11.806 |
| AVERAGE | 142.6 | 6.022 | 6.625 | 9.329 | 9.329 |

EPOCH TIME (days after J2000) = **1222.895303**. (The `****` in LEG 1 Earth vinfm is
as printed — the first leg has no inbound Earth v∞; Russell prints asterisks.)

The full "DATA NECESSARY TO REPRODUCE" block (15 legs × time-start / vinfx,y,z /
time-dv / dvx,y,z to ~7 significant figures in E-notation) is present on p.203 but is
NOT transcribed digit-for-digit here — see readability/scope note below.

### Parent-cycler index of Appendix C (shorthand + parent number, as read pp.203–220)

This is the readable index of which parent cyclers' lowest-Δv trajectories appear.
The shorthand is the **Chapter-4 nomenclature** (already in the catalogue as
`russell-ch4-*`). Verbatim shorthand + parent number (and printed total-Δv in km/s):

| Parent cycler shorthand | parent # | total Δv (km/s) | period (yr) |
|---|---|---|---|
| 6.399G1 | 1 | 0.000000 | 14.56 |
| 4.991gG2 | 83 | 0.000000 | 30.07 |
| 8.049gGf2 | 188 | 0.436091 | 29.95 |
| 8.165Gfh-f2 | 192 | 1.677476 | 30.36 |
| 8.165Gfh+f2 | 193 | 2.611810 | 30.34 |
| 3.353Gg2 | 195 | 0.427585 | 30.26 |
| 3.40 ... Gff3 (3.406Gff3 region) | 19 | 0.148267 | 44.74 |
| 3.406gGf3 (3.40…) | 20 | 0.681385 | 45.02 |
| 3.40… Gfh-3 (3.406Gfh-3) | 23 | 1.819780 | 44.78 |
| 3.418gGff3 | 33 | 0.009457 | 45.14 |
| 3.418gGf3 | 34 | 0.019571 | 45.20 |
| 3.639gGf3 | 49 | 0.007701 | 44.79 |
| 3.768Gh+3 | 54 | 0.000000 | 44.66 |
| 3.768Gh+3 | 58 | 0.000000 | 44.69 |
| 3.768Gh-f3 | 56 | 0.033288 | 45.04 |
| 3.768Gfh-fff3 | 62 | 0.037902 | 44.81 |
| 3.768Gh-ff3 | 64 | 0.007746 | 44.93 |
| 3.768Gh+ff3 | 65 | 0.008819 | 44.91 |
| 3.784Gg3 | 74 | 0.074448 | 44.71 |
| 5.12 ... gGgf3 (5.12sGgf3 region) | 89 | 0.210000 | 44.85 |
| 5.301gGff3 | 94 | 0.247075 | 44.78 |
| 5.301gGff3 | 96 | 0.004582 | 45.04 |
| 5.301gGff3 | 97 | 0.184814 | 44.71 |
| 5.301gGff3 | 98 | 0.429957 | 44.75 |
| 5.301gGfgf3 | 99 | 0.391556 | 44.79 |
| 5.301gtGff3 (5.301g…) | 100 | 0.205926 | 44.62 |
| 5.301gtGff3 | 102 | 0.086683 | 44.69 |
| 5.219Ggh-3 | 111 | 0.021335 | 44.70 |
| 5.219Ggh+3 | 112 | 0.000079 | 45.10 |
| 5.219Ggh+f3 | 113 | 0.004426 | 45.02 |
| 5.219Ggfh+3 | 114 | 0.005993 | 44.90 |
| 5.225Ggg3 | 117 | 0.000000 | 44.74 |
| 5.333gGf3 | 122 | 0.145076 | 45.30 |
| 5.333gg+3 (5.333gGf3 cont.) | 123 | 0.417988 | 44.87 |
| 5.333Ggf3 | 126 | 0.201402 | 45.05 |

(Index continues past p.220 to p.245 for the balance of the 77 solutions; the
shorthand glyphs `g/G/f/F/h` and +/− signs are small in the typeset and a few are
ambiguous between `f` and `t`, marked "…" above. The parent-cycler **numbers**
(1–203) are unambiguous and are the reliable cross-reference key.)

### Readability / scope caveat (HONEST, per rule 5)

Appendix C is **readable** (not rotated, not image-only) but is **43 pages of dense
fixed-width high-precision numeric dumps** — each of the 77 solutions has 7–43 legs,
each leg with 8 numbers in ~7-significant-figure E-notation. Transcribing all of
these digit-for-digit into Markdown would be both impractical and a high
transcription-error risk that violates the "exact published digits, never guess"
rule. Therefore:

- The **format definitions, column meanings, frame/epoch conventions, the
  parent-cycler index, and one complete worked transit-leg table (6.399G1 #1)** are
  transcribed verbatim above.
- The **per-leg "DATA NECESSARY TO REPRODUCE" numeric blocks** are recorded as
  PRESENT-but-not-transcribed. If a future catalogue row needs the exact DE405
  reproduction data for a specific parent cycler, pull that single page fresh and
  transcribe only that one block (the prose confirms only "Epoch time + detailed leg
  information" are needed to reproduce). Russell also states (p.201) an electronic
  file of all 4263 solutions is available "contact the author" — i.e. the canonical
  machine-readable source is not the PDF.

---

## Catalogue-eligibility recommendations (RECOMMENDATION ONLY — no writebacks)

### Already catalogued (do NOT re-add)
- **Table 3.4** members → `russell-2004-t34` (Aldrin 1.0.1.-1 plus several
  `russell-ocampo-*` rows). The verbatim AR/TR/V∞/turn-angle data above can be used
  to **backfill** any partially-populated existing `russell-2004-t34` rows
  (`flyby_mechanics`, AR/TR invariants) — but that is a catalogue edit, out of scope
  here.
- **Tables 3.9 / 3.10 / 3.11** members → already catalogued as 58 + 99 = 157
  `russell-ocampo-5.*`/`6.*` rows under `russell-2004-t39_311`. Count matches exactly.
- **Appendix C** parent cyclers → already catalogued as `russell-ch4-*` rows.

### Catalogue-ELIGIBLE but NOT yet present (candidates for a future writeback pass)
Repeating E-M cyclers in Tables 3.4/3.9/3.10/3.11 that are NOT among the existing
`russell-ocampo-*` / `russell-2004-t34` ids. Spot-check of Table 3.4 vs the catalogue:
the catalogued Table-3.4 set is a **subset** (Aldrin + a handful of `ocampo` rows).
The remaining ~38 Table-3.4 rows (e.g. 3.1.1.+3, 3.3.1.+2, 3.5.x, 3.7.1.+1, 3.9.1.+0,
4.0.3.+1, 4.1.x, 4.5.x, 4.6.x, 4.7.x, 4.8.x, 4.9.x, 4.10.x–4.14.x) are
**catalogue-eligible** repeating circular-coplanar E-M cyclers with full AR/TR/V∞/
turn-angle data — sufficient for a row of the same shape as the existing
`russell-ocampo-*` entries. Recommend adding under `russell-2004-t34` if the
catalogue is to fully cover Russell's circular-coplanar member list.

### Method-artifact / NOT catalogue-eligible
- **Tables 3.5–3.8** are per-encounter Δv/state tables for FOUR specific cyclers
  (2.5.1.+0, 3.1.2.+1, 4.3.1.-5, 4.5.2.-2) — these are **detail tables for cyclers
  already catalogued**, not new members. They are reproduction data (useful for
  validating an existing row's geometry), not new catalogue rows.
- **Appendix C per-leg numeric blocks** are the DE405 reproduction artifacts for
  parent cyclers already catalogued as `russell-ch4-*`. They are accurate-model
  reproduction data, not new family members. The parent cyclers themselves are
  catalogued; the per-launch-window trajectory dumps are not individually
  catalogue-eligible (they are 4263 optimizer outputs, of which 77 are printed).

---

## Summary

| Table | Page | Model | Rows | Status |
|---|---|---|---|---|
| 3.4 | 83 | circular-coplanar | 44 | transcribed verbatim; subset already in catalogue (`russell-2004-t34`) |
| 3.5 | 84 | circular-coplanar (state) | 1 cycler | transcribed verbatim (detail table) |
| 3.6 | 84 | circular-coplanar (state) | 1 cycler | transcribed verbatim (detail table) |
| 3.7 | 84 | circular-coplanar (state) | 1 cycler | transcribed verbatim (detail table) |
| 3.8 | 84 | circular-coplanar (state) | 1 cycler | transcribed verbatim (detail table) |
| 3.9 | 90 | circular-coplanar | 58 | transcribed verbatim; already in catalogue (`russell-2004-t39_311`) |
| 3.10 | 91 | circular-coplanar | 51 | transcribed verbatim; already in catalogue |
| 3.11 | 92 | circular-coplanar | 48 | transcribed verbatim; already in catalogue |
| App. C | 201–245 | ephemeris / DE405 | 77 solutions | format + index + 1 worked example transcribed; per-leg numeric blocks recorded PRESENT-but-not-transcribed (43 pp. dense E-notation; high transcription-error risk) |

Total verbatim member rows transcribed: 44 + 58 + 51 + 48 = **201** circular-coplanar
cycler rows, plus 4 detail/state tables, plus the Appendix C format + parent index.
157-row cross-check (T3.9+T3.10+T3.11) confirmed against p.88 prose.
