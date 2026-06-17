# Nock-Friedlander 1987 (Acta Astronautica 15(6/7)) — full-page deep-read verdict

**Read 2026-06-17 AET.** Part of Mars-cycler digest wave (Agent A).

## 1. Header (verbatim from PDF p.1)

- Title: **"ELEMENTS OF A MARS TRANSPORTATION SYSTEM†"**
- Authors: **Kerry T. Nock** (Jet Propulsion Laboratory, Caltech, 4800 Oak Grove
  Drive, Pasadena CA 91109 USA) **and Alan L. Friedlander** (Science Applications
  International Corporation, Schaumburg IL USA).
- Venue: **Acta Astronautica, Vol. 15, No. 6/7, pp. 505-522, 1987**. Pergamon
  Journals Ltd. ISSN 0094-5765/87.
- DOI (from filename, not printed on PDF): **10.1016/0094-5765(87)90189-5**.
- Received 6 January 1987.
- Footnote †: **"Paper IAA-86-466 presented at the 37th Congress of the
  International Astronautical Federation, Innsbruck, Austria, 4-11 October
  1986."** ← This IS IAA-86-466, the IAF-1986 paper.
- Length: **18 pages** (pp. 505-522).

## 2. What the paper actually is

A **Mars-transportation systems comparison study** for a hypothetical 2035-era
permanent Mars base of 20 people. It's a SYSTEMS paper, not a trajectory-design
paper — but it compares four trajectory options head-to-head over a 15-year
operating cycle (1995-2010 by analogy) by total propellant, vehicle count, and
crew schedule:

1. **Conjunction class** — traditional minimum-energy Hohmann.
2. **VISIT orbits** (Niehoff [ref 13] / Hollister [ref 8]) — heliocentric
   periodic orbits 4:5 and 2:3 commensurabilities with Earth, near-tangent to
   Mars perihelion (~1.38 AU) and Earth orbit. Low V∞.
3. **Up/Down Escalator** (Aldrin [ref 14] = cycler) — 2-1/7 year period
   heliocentric ellipses with aphelion >Mars orbit, high V∞.
4. **Down Escalator** + Phobos NEP cargo variant.

The trajectory-shape content lives in Section 4 (pp. 514-518) and Tables 4, 9,
10. Sections 1-3 are Mars-base + spaceport + vehicle infrastructure scoping
(SP-100 reactors, CASTLE 460 mt design, Phobos propellant plant, etc.).

### Title / venue / DOI cross-check against our records

| Field | Our records | Paper actual | Status |
|---|---|---|---|
| Wishlist entry on #116 | "Friedlander 1986" | **Nock & Friedlander, Acta Astro 15(6/7), 1987, pp. 505-522**, but **= IAA-86-466 presented Oct 1986**. | The "1986" date refers to the IAF conference; the **Acta Astro publication year is 1987**. Both forms refer to the same paper (footnote † on p.505 makes this explicit). |
| Wave digest 2026-06-17 entry | "Nock-Friedlander 1987" | Correct (publication date is 1987). | OK. |
| KNOWN_CORPUS | No anchor for Nock-Friedlander 1987 in `literature_check.py`. | — | This paper resolves the #116 Friedlander wishlist but does NOT need a new KNOWN_CORPUS anchor (see §4 below for reasoning). |

**This paper is NOT the AIAA-86-2009 paper.** Reference 12 of this paper
(p. 522) is **Friedlander, A. L., Niehoff, J. C., Byrnes, D. V. and Longuski,
J. M., "Circulating transportation orbits between Earth and Mars,"
AIAA Paper No. 86-2009-CP, AIAA/AAS Astrodynamics Conference, Williamsburg
VA (1986)** — i.e. the AIAA-86-2009 paper is a SEPARATE Williamsburg
August 1986 paper, distinct from this IAA-86-466 / Acta Astro 1987 paper.
Both have Friedlander as a co-author. The #116 Friedlander wishlist may have
been targeting EITHER paper.

## 3. Key numerical / structural content

### Solar-system model (implicit)
No explicit constants; the paper consumes pre-computed trajectory ΔV ranges
from Friedlander's reference 12 (AIAA-86-2009). Mars surface fly-around
parameters appear in the Mars-orbit-insertion diagrams (Fig. 11, 12) using
**Mars-radius scaling R_M**.

### Cycler-relevant numerical content

**Table 2 (p.510) — ΔV requirements for staging in Earth-Moon space (m/s):**

| Staging location | Trans-Mars injection | Transport from LEO | Transport from lunar surface |
|---|---|---|---|
| LEO | 4470 | — | 2670 |
| GEO | 3540 | 3820 | 3520 |
| Earth-Moon L₁ | 2050 | 3670 | 2510 |
| **Earth-Moon cycler** | **1408** | **3058** | **2550** |
| Lunar orbit | 2230 | 3880 | 1730 |
| Lunar surface | 3960 | 5610 | — |

(Footnote (1): "Escape energy C₃=30 (km/s)² typical of Up Escalator transport
mode, assumes ideal geometry and final injection burn at 6878 km perigee.")

The **Earth-Moon cycler** row is a transport mode for moving propellant/payload
WITHIN cislunar space; **NOT** the same as the Earth-Mars cycler. The text
(p.510) clarifies: "The L₁ point was selected due to its natural
stationkeeping characteristics and the delta-V penalty imposed on the cycler
during non-ideal launch geometry situations."

**Section 4.1 (p.514) — Conjunction class properties** (NOT cycler):
- C₃ between 9 and 16 km²/s² (low launch energy).
- Mars arrival/departure + Earth return V∞ = **2.5 to 4.0 km/s**.
- Earth-to-Mars transit time: **200-350 days**.
- Stopover time: 330-520 days.
- Mars-to-Earth return: 190-360 days.
- Total round-trip ≈ 945-995 days ≈ **2.7 years on average**.

**Section 4.1 — VISIT orbit properties** (Niehoff, NOT the Aldrin cycler):
- Two orbits (Fig. 9, p.514): VISIT-1 (period 1.25 yr, 4:5 with Earth, 3:2
  with Mars) and VISIT-2 (period 1.5 yr, 2:3 with Earth, 5:4 with Mars).
- Aphelion 1.38-1.66 AU (near Mars orbit).
- "Earth encounters will occur approximately once every 5 years and Mars
  encounters once every 3.75 years" (p.514).
- **V∞ at Earth: 4.2-5.2 km/s. V∞ at Mars: 3.7-3.9 km/s.** (p.515 "The relative
  velocity characteristics of the VISIT orbit are 4.2-5.2 km/s at Earth
  encounters and 3.7-3.9 km/s at Mars encounters.")

**Section 4.1 — Up/Down Escalator (Aldrin cycler) properties (p.515):**
- "Both ellipses with a **2 1/7 year period** and an **aphelion point located
  beyond the orbit radius of Mars**. The Up Escalator orbit is oriented so that
  the Earth to Mars transfer time is approximately six months while the Down
  Escalator orbit is oriented to allow the same amount of time for the return
  trip." (p.515)
- Synodic apse rotation **51.4 deg per orbit** ← matches Byrnes 1993 exactly.
- "Aphelion distance of about **2.32 AU**, thereby crossing the orbit of Mars
  at steeper (nontangential) angles."
- **V∞ at Earth: 5.4-6.2 km/s. V∞ at Mars: 6.1-11.7 km/s.** (p.515: "Relative
  velocities of the Escalator orbit are 5.4-6.2 km/s at Earth and 6.1-11.7
  km/s at Mars.") ← **MATCHES catalogue's `aldrin-cycler` Mars V∞=9.7 (midrange)
  and Byrnes 1993 Mars range 6.05-11.74.**

**Table 4 (p.515) — Delta-V requirements for the primary trajectory options,
variation over 15 year cycle (km/s):**

| Transport mode | Launch from L₁ to Mars | Mars orbit capture | Launch from Mars orbit to Earth | Capture at L₁ |
|---|---|---|---|---|
| Conjunction | 1.12-1.47 | 0.96-1.90 | 1.03-1.65 | 1.12-1.44 |
| Visit orbit | 1.91-2.54 | 0.81-1.00 | 2.36-2.55 | 1.91-2.54 |
| Up/Down Escalator | 2.31-2.74 | 0.80-2.58 | 4.09-10.25 | 2.31-2.74 |

(Footnote: "Conjunction mode: all-propulsive capture to 1.15 × 7.95 R_M
parking orbit (0.75-1.05 km/s Taxi to Phobos). Visit and Escalator Taxis:
aero-propulsive capture to Phobos.")

Additional infrastructure V requirements (Table 4 footnotes, p.515):
- LEO to L₁ tanker/transporter = 4.41 km/s (round trip)
- Lunar surface to L₁ tanker = 5.21 km/s (round trip)
- Mars to Phobos shuttle = 7.56 km/s (round trip)
- Conjunction OTV booster return to L₁ = 1.14-1.49 km/s
- Conjunction OTV booster return to Phobos = 1.54-2.08 km/s
- Escalator midcourse adjustment = 0-1.16 km/s
- Midcourse navigation, all modes = 0.05 kg/s (planet-planet leg)

**Table 5 (p.519) — Summary comparison of propellant and consumables over 15
years (mt):**

| Mode | LOX from Moon | LH₂ from LEO | LOX from Mars | LH₂ from Mars | LOX from Phobos | LH₂ from Phobos | CASTLE consumables | TOTAL |
|---|---|---|---|---|---|---|---|---|
| Conjunction | 22519 | 3217 | 919 | 131 | 2992 | 427 | 206 | **30205** |
| Visit | 17006 | 2430 | 656 | 94 | 1150 | 164 | 813 | **21500** |
| Up/Down Escalator | 21975 | 3138 | 919 | 131 | 5834 | 833 | 269 | **32830** |
| Down Escalator | 16216 | 2317 | 459 | 66 | 4027 | 575 | 279 | **23660** |

**Table 9 (p.521) — Transportation mode comparison chart** (cycler-relevant
rows):

| Property | Conjunction | Visit | Up/Down Escalator | Down Escalator |
|---|---|---|---|---|
| Number of CASTLEs | 2 | 3 | 2 | 1 |
| Tour of duty (yr) | 4.8 | 5.7-7.9 | 5 | 6.5 |
| Mars crew flightime (yr) | 1.6 | 1.2-6.3 | 0.9 | 2.1 |
| Mars staytime (yr) | 3.2 | — (1.6-5.9) | 4.1 | 4.4 |
| Number of sorties per 15 years | 17 | 21 | 17 | 17 |
| Number of personnel to Mars vicinity and return in 15 years | 119 | 147 | 133 | 105 |

**Table 10 (p.521) — CASTLE velocity change requirements (m/s):**

| | Conjunction | VISIT | Up/Down Escalator | Down Escalator |
|---|---|---|---|---|
| M/s per 15 years | 36345ᵇ | 1000 | 5159ᶜ | 2734 |
| M/s per CASTLE per 15 years | 18172 | 333 | 2580 | 2734 |
| M/s per CASTLE per sortie (avg) | 4543 | 125 | 369 | 391 |
| Number of major ΔVs per CASTLE | 16 | 0 | 3 | 3 |
| M/s per major ΔV (avg) | 1298 | 0 | 627 | 627 |
| M/s per major ΔV (min) | 1011 | 0 | 270 | 270 |
| M/s per major ΔV (max) | 1646 | 0 | 1105 | 1105 |

ᶜ "Taxi propulsion system provides CASTLE ΔVs"

**Up/Down Escalator: 3 major ΔVs per CASTLE per 15-year cycle, 627 m/s mean,
range 270-1105 m/s.** ← This is the second independent source corroborating
Byrnes 1993 Table 1's 3-maneuvers-per-cycle observation (0.45, 0.54, 0.74 km/s
outbound = 450, 540, 740 m/s, within Table 10's 270-1105 m/s range).

### Cited references for cycler-canonical lineage
- Ref 8: **Hollister, W. M., "Castles in space." Astronautica Acta. (1967)**.
- Ref 12: **Friedlander, A. L., Niehoff, J. C., Byrnes, D. V. and Longuski,
  J. M., "Circulating transportation orbits between Earth and Mars," AIAA
  Paper No. 86-2009-CP, AIAA/AAS Astrodynamics Conference, Williamsburg, Va
  (1986).** ← THE separate Friedlander-Niehoff-Byrnes-Longuski paper.
- Ref 13: **Niehoff, J., "Manned Mars mission design, steps to Mars,"
  Joint AIAA/Planetary Society Conference, National Academy of Sciences,
  Washington DC (1985).**
- Ref 14: **Aldrin, E. E., "Cyclic trajectory concepts," SAIC presentation
  at the Interplanetary Rapid Transit Study Meeting, Jet Propulsion Laboratory
  (1985).** ← The same "Aldrin 1985 SAIC presentation" the Byrnes 1993 JSR
  paper's ref 8 cites.

## 4. `KNOWN_CORPUS` impact (RECOMMEND, do not edit)

**Does this paper need a new KNOWN_CORPUS anchor? HONEST NEGATIVE.**

Nock-Friedlander 1987 is a systems-level transportation-system comparison
study that **consumes** Aldrin cycler / VISIT orbit / conjunction trajectory
properties from other primary sources. It does NOT publish new cycler family
discoveries or new orbital elements. Its quantitative content is:

- Per-mode total propellant masses (Table 5).
- Per-mode CASTLE ΔV budgets (Table 10).
- VISIT and Up/Down Escalator V∞ ranges (Section 4.1 prose).

Of these, only the V∞ ranges (5.4-6.2 km/s @ Earth, 6.1-11.7 km/s @ Mars for
Up/Down Escalator) are cycler properties — and **they are entirely consistent
with the Byrnes 1993 JSR Tables 1-2** (Earth 5.67-6.19 km/s outbound /
5.39-6.01 km/s inbound; Mars 6.05-11.74 km/s outbound / 6.60-11.59 km/s
inbound). So Nock-Friedlander 1987 is a **secondary corroborating source**,
not a primary source.

**Recommendation**: Do NOT add a new KNOWN_CORPUS anchor for Nock-Friedlander
1987. Instead, add it as a **corroborating citation in the existing
`aldrin-cycler` catalogue row** (line ~135 corroborating_sources list).

If the parent wants a wider Mars-systems-context KNOWN_CORPUS anchor (covering
VISIT + Aldrin + Conjunction transport modes for a permanent Mars base), then
add Nock-Friedlander 1987 as the primary citation for that systems-anchor.
But the cycler-specific anchors already exist (`russell-ocampo` for SnLm,
`aldrin-cycler` for Aldrin); no new anchor is needed.

## 5. Catalogue impact (RECOMMEND)

### Rows touched

- **`aldrin-cycler`** (line ~30 of catalogue.yaml) — Nock-Friedlander 1987
  Table 4 / Table 10 / Section 4.1 prose all corroborate the Up/Down Escalator
  (Aldrin cycler) V∞ ranges and ΔV maintenance budget. **Add to
  corroborating_sources list.** Already V3 cross_validated; this doesn't change
  the V-level.

### V0→V1 promotion candidacy

**HONEST NEGATIVE.** Nock-Friedlander 1987 prints range-bound aggregate values
(V∞ ranges, ΔV ranges) for the Up/Down Escalator mode, NOT per-encounter or
per-arc orbital elements. **It is by construction insufficient for V0→V1
promotion** under the §14 like-for-like rule (which requires per-member orbital
elements reproducible in same-model corrector). No specific catalogue row gets
promoted.

### What this paper IS good for

- **#116 Friedlander wishlist resolution**: confirmed. This paper IS the
  Friedlander 1986/87 paper the wishlist refers to (under either date — IAF
  Oct 1986 presentation, Acta Astro 1987 publication). Mark wishlist-resolved.
- **Cross-validates Byrnes 1993 numerical results** within its quoted ranges.
- **Documents the Aldrin 1985 SAIC presentation** existence (ref 14) and the
  separate **AIAA-86-2009 Friedlander-Niehoff-Byrnes-Longuski paper** (ref 12).
  Both are PRIMARY sources we should pursue if the parent wants to chase them
  (the AIAA-86-2009 paper may have richer per-arc data than Nock-Friedlander
  1987 itself).

## 6. Errata / surprises

1. **The wishlist date is "1986" but the publication is 1987.** Defensible
   either way — IAF conference Oct 1986 → Acta Astro publication 1987. The
   IAA paper number IAA-86-466 is the unambiguous identifier.
2. **This is NOT the AIAA-86-2009 paper.** They are sister papers from the
   same authors (overlapping author lists) at two different conferences. If
   the catalogue or KNOWN_CORPUS cites "Friedlander 1986" referring to
   AIAA-86-2009, this paper does NOT close that wishlist item — it closes a
   different one.
3. **The "Earth-Moon cycler" row in Table 2** is a **cislunar transport mode**
   between LEO/L₁/Lunar Surface for staging propellant, NOT an interplanetary
   cycler. Don't conflate with the Aldrin Earth-Mars cycler.
4. **The Aldrin 1985 SAIC presentation (ref 14)** is the same source Byrnes
   1993 JSR ref 8 cites. **It is the original Aldrin cycler concept document**
   — not in our paper corpus, and probably not in any public archive (SAIC
   internal report). The catalogue's `aldrin-cycler` row's "a≈1.659, e≈0.41,
   peri≈0.98, apo≈2.34" alternative (from spec.md) may trace to this 1985 SAIC
   presentation. If so, it is a `SOURCE_PERMANENT_NEGATIVE` (the document is
   not accessible).
5. **No mention of S1L1 / SnLm nomenclature.** This paper predates the
   McConaghy SnLm naming convention (2002+). The cycler discussion uses only
   the "Up Escalator" / "Down Escalator" terminology.
6. **The Up/Down Escalator mode requires 3 major ΔVs per CASTLE per 15-year
   cycle (Table 10).** Matches Byrnes 1993 JSR exactly. **This is a
   structural anchor**: the Aldrin cycler is operationally a "3 maneuvers per
   15 years per CASTLE" trajectory.

## 7. Action items for parent

- [ ] **Mark #116 Friedlander 1986 wishlist as RESOLVED** by Nock-Friedlander
  1987 (Acta Astro 15(6/7):505-522, IAA-86-466).
- [ ] **Add Nock-Friedlander 1987 to `aldrin-cycler` corroborating_sources**
  in catalogue.yaml. Citation: `Nock, K. T. & Friedlander, A. L., "Elements
  of a Mars Transportation System," Acta Astronautica 15(6/7):505-522 (1987),
  IAA-86-466, DOI 10.1016/0094-5765(87)90189-5; Sec. 4.1 + Tables 4, 9, 10
  give Up/Down Escalator V∞ and ΔV ranges consistent with Byrnes 1993 JSR.`
- [ ] **NO new KNOWN_CORPUS anchor** needed. This paper is corroborating, not
  primary.
- [ ] **Track the Friedlander-Niehoff-Byrnes-Longuski AIAA-86-2009 paper
  separately** as a still-outstanding acquisition (it is referenced as ref 12
  here; it MAY have richer per-arc data than Nock-Friedlander 1987 itself and
  is the source the McConaghy 2002 / Russell-Ocampo 2003 papers cite). Add
  to #116 acquisitions wishlist if not present.
- [ ] **Document the Aldrin 1985 SAIC presentation as a `SOURCE_PERMANENT_
  NEGATIVE`** since it appears to be an SAIC internal report not in any
  public archive. Source-permanent under the spec §16.x SOURCE_PERMANENT_*
  taxonomy.
- [ ] **NO action on V0→V1 promotions.** This paper publishes range-bound
  aggregate values, not per-arc orbital elements; insufficient for §14
  like-for-like.
