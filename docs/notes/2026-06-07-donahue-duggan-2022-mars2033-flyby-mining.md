# Donahue & Duggan 2022 (Boeing) — Mars 2033 Human Flyby — data mining

Mined 2026-06-07. A **DATA-focused** mine: sourced 2033/2035 free-return anchors
(dates, C3, V∞, TOF, flyby altitude, ΔV budget) for cross-checking our
launch-window / free-return machinery, plus a topology verdict against the
Hughes 2014 EVME finding (`docs/notes/2026-06-07-hughes-2014-fast-mars-free-returns-mining.md`).

**Source (cite exactly — authors/title/IAC number only):**
Donahue, B., & Duggan, M., "Mars 2033 Human Flyby Mission," IAC-22,B3,8,x70674,
73rd International Astronautical Congress (IAC), Paris, France, 18–22 September
2022. Boeing Exploration Launch Systems (Huntsville, AL) + Boeing Space Systems
(Houston, TX). © 2022 IAF.

> 17 pp., clean digital typeset. Trajectory numbers live in **Figure 4** (p.10,
> "2033 Mars Flyby Trajectory Diagram", annotated), **Figure 18** (p.17, "2035
> Opportunity Mars Flyby Trajectory Diagram", annotated), and the **ΔV-budget
> Tables 1-A / 1-B** (p.11). Body text pp.2-7 gives the topology, constraints,
> and a few numbers that disagree slightly with the figures (logged in §6).
> Most other figures (Figs 7-9, 16-17) are renders; Figs 10-15 are launch-
> manifest mass bars (mission-architecture, not trajectory).

---

## 0. HEADLINE VERDICT

**TOPOLOGY: a single-arc, MARS-ONLY ballistic free return with an inbound
powered Deep-Space Burn (DSB) — NOT a Venus gravity assist.** This is the
key difference from Hughes 2014's EVME family:

- The **2033** baseline is **Earth → Mars flyby → (DSB near Venus's orbital
  radius) → Earth**. Mars is the only gravity-assist body. "Venus" appears
  *only* as a heliocentric distance marker for where the inbound perihelion /
  DSB occurs ("reaches perihelion (Venus AU=0.7, Fig 4)", p.5). **There is no
  Venus flyby in 2033.** The outbound Earth→Mars leg and the Mars flyby itself
  are purely **ballistic** ("No propulsive maneuver is required, as the vehicle
  is on a free return trajectory", p.5); the only deterministic deep-space burn
  is the **inbound DSB** added to pull the Earth-entry speed down from >14 km/s
  to 12.5 km/s.
- The **2035** alternate (Fig 18, p.17) **DOES include a Venus passage on the
  inbound leg** — "the 2035 passage employs a close Venus passage on the inbound
  leg that eliminates the need for any propulsive maneuvers (such as a DSB)
  following Earth departure" (p.7). So 2035 is the closer analogue to the Hughes
  EVME free-return topology (Mars flyby then Venus flyby, fully ballistic).

**Mars-encounter V∞ (the number our pipeline cares about): 5.113 km/s (2033),
7.037 km/s (2035)** — both **above** the Jones Mars-~3 km/s cycler class and
squarely **inside the Hughes 2014 "5-7 km/s class" for ballistic Mars free
returns** (Hughes Fig 1 Tisserand: Mars ~5 km/s). This is independent third-party
corroboration that *ballistic* Mars free returns sit at ~5-7 km/s Mars V∞, and
that getting lower needs added energy management (Hughes: broken-plane DSM;
here: an inbound DSB, but to lower *entry* speed, not Mars V∞).

**The mission is NOT a cycler** (one-shot human flyby) → **NOT catalogue-
eligible.** Confirmed in §4.

---

## 1. GOLDEN-ELIGIBLE ANCHORS (verbatim, with page numbers)

### Anchor A — 2033 Mars Flyby (Figure 4 annotations, p.10) — PRIMARY

Free-return, optimized, opening day of a 21-day launch period.

| Event | Date | Quantity (verbatim from Fig 4) |
|---|---|---|
| Earth departure | Dec 2032 (8/12 → 12/8/2032) | C3 = 36.35 km²/s²; TMI ΔV = 1,643 m/s; Time = 0 days |
| Mars passage | 8/9/2033 | Altitude 250 km; **V∞ = 5.113 km/s**; No propulsive burn; Time = 244 days |
| Inbound perihelion / DSB | 2/7/2034 | ΔV = 1,465 m/s; Distance = Venus AU; Time = 426 days |
| Earth arrival | 5/22/2034 | V_entry = 12,500 m/s (ref) / 13,900 m/s (option); Time = **531 days** |

Sidebar (Fig 4, p.10): "Earth Arrival Velocity 12.5 km/s Reference, 13.9 km/s
Option"; "return hyperbolic excess speed is 3.69 km/s" (p.3).
Abstract/Summary trip time = "530 days" (p.1, p.7) ≈ Fig-4 531 days.

### Anchor B — 2035 Mars Flyby (Figure 18 annotations, p.17) — BACKUP OPPORTUNITY

The "missed-2033" fallback; includes an inbound Venus passage (ballistic, no DSB).

| Event | Date | Quantity (verbatim from Fig 18) |
|---|---|---|
| Earth departure | 1/23/2035 | C3 = 66.18 km²/s²; Time = 0 days |
| Pass Mars | 8/13/2035 | V∞ = 7,037 m/s; Altitude 250 km; Time = 202 days |
| Pass Venus | 4/8/2036 | V∞ = 6,491 m/s; Altitude 10,855 km; Time = 441 days |
| Arrive Earth | 8/18/2036 | V∞ = 5,416 m/s; V_entry = 12,345 m/s; Time = **573 days** |

Body text cross-checks (p.7): Mars passage 250 km, V∞ "7,037 m/s"; Venus closest
approach 10,855 km; Earth arrival V∞ 5,416 m/s → entry 12,345 m/s ("highest
Earth Entry Capsule speed for this mission"); "natural arrival entry speed at
Earth is less than 12.5 km/s for all launch dates" → permits direct entry, **no
DSB**. Earth-departure TMI ΔV from EM-L2 = 2,641 m/s opening day. Total mission
trip time 573 days, "varies only 1 day over the entire launch period." Total
mission ΔV for 2035 < total for 2033.

### Anchor C — 2033 ΔV budgets (Tables 1-A / 1-B, p.11, verbatim)

The full maneuver ledger (note the dates are mistyped in the table as "1032"
for "2032" — clearly the same Dec-2032/2033/2034 timeline as Fig 4).

**Table 1-A — Habitat Expended (reference):**

| Maneuver date | ΔV (m/s) | Event |
|---|---|---|
| 11/1/[20]32 | 9 | EM-L2 orbit departure |
| 11/3/[20]32 | 41 | Mid-course correction |
| 11/10/[20]32 | 202 | Powered lunar swingby |
| 11/15/[20]32 | 13 | First perigee ΔV |
| 12/4/[20]32 | 4 | Third apogee ΔV |
| 12/8/[20]32 | 1,649 | TMI at fourth perigee |
| 8/9/2033 | 0 | Pass Mars |
| 2/7/2033 [sic; =2/7/2034] | 1,461 | Deep space burn |
| 5/22/2034 | 0 | Earth entry |
| **Total** | **3,379** | Total Mission ΔV |

**Table 1-B — Habitat Recovered (alternative, captures Hab back to EM-L2):**

| Maneuver date | ΔV (m/s) | Event |
|---|---|---|
| 10/29/[20]32 | 9 | EM-L2 orbit departure |
| 10/31/[20]32 | 41 | Mid-course correction |
| 11/7/[20]32 | 202 | Powered lunar swingby |
| 11/12/[20]32 | 13 | First perigee ΔV |
| 12/1/[20]32 | 4 | Third apogee ΔV |
| 12/5/[20]32 | 1,674 | TMI at fourth perigee |
| 8/8/2033 | 0 | Pass Mars |
| 2/12/2033 [sic; =2/12/2034] | 2,115 | Deep space burn |
| 6/13/2034 | 707 | Earth orbit insertion |
| 6/17/2034 | 4 | First apogee ΔV |
| 7/6/2034 | 13 | Fourth perigee ΔV |
| 7/11/2034 | 202 | Powered lunar swingby |
| 7/18/2034 | 41 | Final mid-course correction |
| 7/20/2034 | 9 | EM-L2 orbit insertion |
| **Total** | **5,034** | Total Mission ΔV |

Key invariant in BOTH ledgers: **Pass Mars ΔV = 0** (ballistic flyby). The only
large deep-space burn is the inbound DSB (1,461 / 2,115 m/s), plus the
Earth-departure TMI (~1,649-1,674 m/s) done in cis-lunar space, not Mars.

---

## 2. TOPOLOGY vs HUGHES 2014 (the requested comparison)

| Aspect | Hughes 2014 EVME (AIAA 2014-4109) | Donahue-Duggan 2022 (this paper) |
|---|---|---|
| GA bodies | **Venus then Mars** (Earth-Venus-Mars-Earth) | **2033: Mars only.** 2035: Mars then Venus |
| Venus role | Powered/unpowered Venus flyby (energy bend) | 2033: NONE — just the heliocentric radius of inbound perihelion. 2035: an actual inbound Venus flyby |
| Deterministic ΔV | **ZERO** (pure free return; lower-V∞ gated behind broken-plane DSM they exclude) | 2033: an **inbound DSB** (1,461-2,115 m/s) — but to cut *Earth-entry* speed, not to enable the return. 2035: ZERO DSB (Venus flyby does it) |
| Mars-encounter V∞ | ~5-7 km/s (Fig 1 Tisserand; not tabulated) | **5.113 km/s (2033), 7.037 km/s (2035)** — explicit, annotated |
| Earth-arrival V∞ | 6.34-9 km/s | 2035: 5.416 km/s (2033 not given as V∞; V_entry 12.5/13.9 km/s) |
| Total TOF | 499-582 d (EVME) | 531 d (2033), 573 d (2035) |
| Mars flyby altitude | 346-363 km (one Table-4 case) | **250 km** (both 2033 & 2035, by design constraint) |
| Search method | STOUR grid (patched-conic, analytic eph) | MAnE™ optimizer + SpaceFlightSolutions software (p.3) |

**Verdict:** Donahue-Duggan's Mars-encounter V∞ values (5.11 / 7.04 km/s)
**land exactly in Hughes 2014's ballistic "5-7 km/s class"** — independent
confirmation that *ballistic* Mars free returns do not reach the Jones
Mars-~3 km/s regime. The 2033 baseline is *not even* a Venus-GA trajectory; it
is a **plain Mars free return** (Patel-Longuski-Sims class, no intermediate
flyby), with a powered inbound DSB added solely to bring Earth-entry speed from
>14 km/s into the Orion TPS limit (12.5 km/s). Where Hughes parks lower V∞
behind an excluded broken-plane DSM, Boeing simply accepts the higher entry
speed and burns it off inbound — a *powered* topology, NOT a pure free return.
The 2035 case (Mars→Venus, zero DSB) is the true ballistic-free-return analogue
of Hughes EVME, and its higher Mars V∞ (7.04) is consistent with Hughes's range.

---

## 3. MAPS-TO-OUR-X verdicts

| Element | Our code / frontier | Verdict |
|---|---|---|
| 2033 single-arc Mars free return, ballistic Mars flyby (ΔV=0 at Mars) + inbound DSB | single-ellipse-per-leg ballistic corrector (`correct.py`); free-return genome (#137) | **MAPS on the ballistic legs**, NOT on the DSB. The outbound E→M and the Mars flyby are single-ellipse ballistic — representable. The inbound DSB makes the *full* return a 2-arc (powered) trajectory — needs the multi-arc / DSM path, mirroring the S1L1 multi-arc finding and Hughes's broken-plane gate. |
| Mars V∞ = 5.113 km/s (2033), 7.037 km/s (2035), annotated | Jones Mars-~3 km/s target; #137 E-M sub-arc ~2.81; #110 floors | **DOES NOT REACH Jones regime** — corroborates Hughes: ballistic Mars free returns are 5-7 km/s class. A *sourced* Mars-V∞ anchor (Hughes had none tabulated — this paper fills that gap). |
| Dec-2032 Earth departure (C3 36.35), 250-km Mars flyby 8/9/2033 | our 10-yr horizon launch-window scans | **NEIGHBOURHOOD ANCHOR** — 2032-2033 is inside our scan horizon. A check target: does our window finder surface a Dec-2032 / Aug-2033 Mars free-return arc near C3≈36 with Mars-pass alt≈250 km, V∞≈5.1? See §5 uncertainty before treating as golden. |
| 2035 Mars→Venus ballistic free return (zero DSB), Venus alt 10,855 km | Hughes EVME family; our VEM hunt | **MAPS to Hughes EVME topology** (Mars-then-Venus ballistic). Note: order is M-before-V here, opposite to Hughes's V-before-M EVME; closer to Hughes's "EMVE" path which Hughes found gave only TOF>800 d — yet Boeing's MVE-return closes at 573 d. Worth flagging as a counterexample to Hughes's EMVE pessimism. |
| MAnE™ optimizer + SpaceFlightSolutions; "trajectory is optimized" (p.3) | our optimizer / verify chain | **DIFFERENT TOOL** — not STOUR, not Copernicus. See §4 provenance. |

---

## 4. v4.2 BACKFILL CHECKS + catalogue-eligibility ruling

**CATALOGUE-ELIGIBLE? NO — confirmed.** A "2033 crewed flyby mission" is a
**one-shot human Mars flyby free return**, explicitly framed as "a precursor to
a later Mars surface mission" (p.1). It does NOT re-encounter a body to sustain
a repeating orbit. **Not a cycler → not a catalogue row.** Valid only as a
**cross-check / sourced-anchor** for a free-return pipeline.

Provenance flags recorded anyway (per the v4.2 checklist) in case any anchor is
ever staged as a pipeline golden:

- **center**: heliocentric (Sun-centered) for the transfer arcs; Mars-centered
  (alt 250 km, unpowered) for the flyby; the DSB is a heliocentric deep-space
  maneuver near r ≈ Venus's orbit (~0.7 AU). EM-L2-centered departure sequence
  for the TMI stack. No catalogue `center` field applies (not a row); if staged,
  legs are `center: "Sun"`.
- **tof_days_bounds**: figures give **per-event cumulative times** (Mars 244 d,
  DSB 426 d, Earth 531 d for 2033; Mars 202 d, Venus 441 d, Earth 573 d for
  2035) — so per-leg ToFs ARE derivable here (2033: E→M 244 d, M→DSB 182 d,
  DSB→E 105 d; 2035: E→M 202 d, M→V 239 d, V→E 132 d). These are READ FROM the
  annotated figure, so flag as **SOURCED (figure-annotated)**, not derived.
  Total TOF: 531 d (2033) / 573 d (2035). 2035 total "varies only 1 day over the
  entire launch period" (p.7) → a genuine narrow `tof_days_bounds` ≈ [573, 574].
- **source_ephemeris**: **NOT a named DE kernel.** The trajectory was "optimized
  with the MAnE™ software from SpaceFlightSolutions" (p.3). So any anchor carries
  `source_ephemeris: MAnE (SpaceFlightSolutions), DE version UNSPECIFIED`.
  **Flag UNSPECIFIED.** (This is a distinct provenance from our existing STOUR
  and Copernicus sources — first MAnE-sourced anchor candidate.)

---

## 5. HONEST UNCERTAINTY / "not golden-grade as-is"

- **Numbers come from annotated figures + a budget table, not a clean data
  table.** Fig 4 / Fig 18 annotations are legible and self-consistent, but they
  are figure callouts, not a tabulated trajectory listing. Treat as **good
  sourced anchors, second-tier to a true data table.**
- **Internal disagreements (logged, not resolved):**
  - DSB ΔV: Fig 4 = **1,465 m/s**; body text p.4 = "1,465 m/s" and p.4 also
    "1,461 m/s"; Table 1-A = **1,461 m/s**. (~4 m/s spread — rounding.)
  - TMI ΔV: Fig 4 = **1,643 m/s**; Table 1-A = **1,649 m/s**; Table 1-B
    (Hab-recovered) = **1,674 m/s**; body text p.5 = "C3=36.2 km²/s², 1,947 m/s"
    for the boosted TMI (different reference point — full TMI from the elliptic
    orbit vs. the perigee burn). Do NOT average; pick the figure that matches
    your reference frame.
  - C3: Fig 4 = **36.35**; text p.5 = **36.2**; text p.3 = **36.35**. Use 36.35.
  - Trip time: abstract/summary "530 days"; Fig 4 "531 days". Use 531 (figure).
- **2033 Mars V∞ = 5.113 km/s** appears in BOTH Fig 4 ("Vinf = 5.113 m/s" [sic,
  km/s]) and body text p.5 ("Mars arrival V-infinity is 5.1 km/s"). Consistent →
  **the most defensible single golden number in the paper.**
- **No outbound DSM, no powered Mars flyby** — confirmed by ΔV=0 at "Pass Mars"
  in both Tables 1-A/1-B. High confidence the Mars flyby is ballistic.
- **2035 C3 = 66.18 km²/s²** is very high (vs 36.35 for 2033) — consistent with
  the text noting 2035 needs no DSB but pays in departure energy; plausible but
  worth a sanity check if used.

---

## 6. CROSS-REFERENCES WORTH ACQUIRING (their citations for trajectory design)

From the References (p.8):

- **[4] Tito, Dennis; MacCallum, Taber; Carrico, John, "Inspiration Mars:
  Feasibility Analysis for a Manned Mars Free-Return Mission in 2008" [sic; 2018],
  May 8, 2013.** — The direct heritage of this whole free-return concept; same
  Inspiration-Mars lineage as Hughes 2014 (which shares Carrico/Loucks/Tito).
  **Acquire — ties the two papers together; the IM baseline both back up.**
- **[3] Z. R. Putnam, R. D. Braun, et al., "Entry System Options for Human Return
  from the Moon and Mars," AIAA 5915, 2005.** — Source of the 12.5 km/s entry
  limit / aerodynamic corridor (Fig 6, Table 2) that *bounds* why the inbound DSB
  exists. **Acquire if we model the entry-speed constraint.**
- **[2] Dunham, David, et al., "Earth-Moon Halo Orbit — Gateway or Tollbooth?,"
  AAS 19-756, 2019.** — Source of the EM-L2 departure sequence (5-maneuver,
  61,347 km) feeding TMI. Departure-side only; lower priority for trajectory.
- **MAnE™ / SpaceFlightSolutions** — the optimizer (no citation given, just named
  p.3). Note for provenance: not STOUR, not Copernicus, not GMAT.

NOT cited here but the obvious upstream for the *no-Venus* 2033 topology:
Patel, Longuski & Sims, "Mars Free Return Trajectories," JSR 35(3), 2002 (the
no-intermediate-flyby Mars free return) — already flagged as acquisition-worthy
in the Hughes 2014 note.

---

## 7. SINGLE MOST DECISIVE FINDING

**Boeing's 2033 baseline is a single-arc Mars-only ballistic free return (Mars
flyby ΔV = 0, V∞ = 5.113 km/s, alt 250 km, 12/8/2032 → 8/9/2033 → 5/22/2034,
531 d, C3 36.35) with an inbound powered DSB (1,461-2,115 m/s) used only to
trim Earth-entry speed — NOT a Venus gravity assist.** Its Mars V∞ (5.11 km/s),
and the 2035 backup's (7.04 km/s), land squarely in Hughes 2014's ballistic
"5-7 km/s class," giving us — for the first time — *explicit, sourced* Mars-
encounter V∞ anchors for a ballistic Mars free return (Hughes tabulated none).
This independently corroborates that ballistic Mars free returns do not reach the
Jones Mars-~3 km/s cycler regime; lower Mars V∞ remains gated behind added
energy management (Hughes: broken-plane DSM; Boeing: an inbound DSB), i.e. a
multi-arc / powered topology our single-ellipse corrector cannot represent.
