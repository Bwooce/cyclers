# Data Mining: Agrawal 2022, Landau-Longuski 2009, Howe 2025

Extracted 2026-06-04. Method: pdftotext where text layer present; pdftoppm (PNG/JPEG) +
Read for image-only or scrambled-text PDFs.

---

## Full Citations

**[AGRAWAL-2022]**
Rachana Agrawal, "Design and Analysis of an Orbital Logistics Architecture for Sustainable
Human Exploration of Mars," PhD Dissertation, Purdue University, 2022. Advisors:
James M. Longuski, Sarag J. Saikia. 144 pp.

**[LANDAU-2009]**
D. F. Landau and J. M. Longuski, "Comparative Assessment of Human-Mars-Mission
Technologies and Architectures," Advances in the Astronautical Sciences (AAS),
Vol. 126, 2007 (?). [NOTE: PDF metadata shows creation date "Wed May 31 05:49:52
2006 AEST"; venue header read from p.1 reads "COMPARATIVE ASSESSMENT OF
HUMAN-MARS-MISSION TECHNOLOGIES AND ARCHITECTURES" — 32 pp.
Reference [43] in Howe 2025 cites this as Landau & Longuski 2009. Exact journal volume
and issue not printed in the scanned PDF. Authorship and institution confirmed: D.F.
Landau and J.M. Longuski, School of Aeronautics and Astronautics, Purdue University,
West Lafayette, IN 47907-2023.]

**[HOWE-2025]**
A. Scott Howe, John Blincow, Theodore W. Hall, Colin Leonard, "Tackling a Mars Cycler
Design Head-On," 54th International Conference on Environmental Systems (ICES),
Paper ICES-2025-555, 13–17 July 2025, Prague, Czech Republic. 15 pp.

---

## Paper 1: Agrawal 2022

### What the paper is

A Purdue PhD dissertation on an **orbital logistics node** ("Mars Spacedock") for sustained
human Mars exploration. The core content is: Spacedock conceptual design; stationing orbit
optimisation (LMO/HMO/HEMO candidates); EDL/landing-site accessibility from orbit;
station-keeping over 15 years; propellant comparison with/without refuelling node.

**This paper does NOT contain cycler trajectory data.** Cyclers are mentioned only
in passing:

- p. 72: "In the long term, Earth-Mars cyclers can be added to the MTV architecture. Since
  minimum two cycler vehicles are required for crewed missions..."
- pp. 129–130 (Future Work): cycler trajectories listed as a topic for future investigation.
  Text: "Earth-Mars cycler has been proposed as alternative concepts for long term human
  Mars exploration [81][82]. The Mars Spacedock would provide even more advantages
  especially with the habitat since the cycler will not stop in Mars orbit and crew will
  travel in a taxi vehicle from the cycler to the Mars orbit. A cycler will also impose more
  constraints on the arrival and departure trajectory."

### Trajectory data for direct missions (not cyclers)

Extracted from Chapter 3 (pp. 68–70, trajectory constraint tables):

| Mission type | Constraint | Value |
|---|---|---|
| Long-stay crewed (E→M one-way) | Max ToF | < 240 days |
| Long-stay crewed (E→M one-way) | Max V∞ at Mars | < 6 km/s |
| Short-stay crewed (round trip) | Total round-trip ToF | < 800 days |
| Short-stay crewed (surface stay) | Surface duration | 30–60 days |
| Cargo, ballistic | — | standard Hohmann-like |
| Cargo, hybrid (E→M) | Max ToF | < 600 days |
| Cargo, hybrid (E→M) | Max cumulative SEP ΔV | < 10 km/s |
| Cargo, hybrid (E→M) | Max V∞ at Mars | < 4.5 km/s |

Synodic periods analysed: 7 opportunities from 2030 to 2044 (p. 75).

### Cycler references cited

- [18] Byrnes, Longuski, Aldrin, "Cycler orbit between Earth and Mars," JSR 30(3):334–336, 1993.
- [81] B. Aldrin, "Cyclic trajectory concepts," SAIC presentation, JPL, 1985.
- [82] S. Pelle et al., "Earth-Mars cyclers for a sustainable human exploration of Mars,"
  Acta Astronautica, vol. 154, pp. 286–294, 2019.

### S1L1 relevance

None. No orbital elements, V∞ by encounter, or per-arc ToFs for any cycler are given.

---

## Paper 2: Landau-Longuski 2009

### What the paper is

A comparative architecture study of six mission classes for human Mars missions, evaluated
primarily by **Injected Mass to Low-Earth Orbit (IMLEO)** as a function of one-way ToF
(120–270 days). The six architectures (Table 1, p.3) are:

1. Direct
2. Semi-Direct
3. Stop-Over
4. Mars-to-Earth (M-E) Semi-Cycler
5. Earth-to-Mars (E-M) Semi-Cycler
6. Cycler

**This paper does NOT give per-encounter V∞ values or orbital elements for any
cycler trajectory.** IMLEO is modelled as a function of ToF using the Zola [54] low-thrust
method (ref. [54] = Zola, "A Method for Approximating Propellant Requirements of
Low-Thrust Trajectories," NASA TN-D-3400, 1966).

### Mission assumptions (p.7)

| Parameter | Value |
|---|---|
| One-way ToF range (parametric study) | 120–270 days |
| Nominal one-way ToF | 210 days |
| Crew size | 4 |
| Taxi capsule mass | 1.5 mt/person |
| Transfer vehicle cabin | 6 mt/person |
| Consumables | 20 kg/person/day |
| Mars stay (long-stay baseline) | ~500 days |
| ΔV to reorient parking orbit at Earth | 350 m/s |
| ΔV to reorient parking orbit at Mars | 180 m/s |

### Architecture ranking (Table 13, p.25)

For the "advanced" scenario (electric propulsion + Mars propellant tanker technology):

| Rank | Architecture |
|---|---|
| 1 | E-M Semi-Cycler |
| 2 | Cycler |
| 3 | Semi-Direct |
| 4 | Stop-Over |
| 5 | M-E Semi-Cycler |
| 6 | Direct |

For "current" technology (chemical, no tanker), ranks shift; direct and stop-over improve
relative to cycler options.

### Key findings on cycler IMLEO (qualitative, from figures pp. 14–24)

- Cycler IMLEO is competitive with (and often lower than) direct/stop-over for ToFs ≥ 180 days.
- E-M semi-cycler generally yields lower IMLEO than full cycler, because crew only cycles
  in one direction (Earth→Mars) and returns by taxi on a standard trajectory.
- IMLEO sensitivity to ToF is steeper for direct missions than for cycler-based missions.
- No tabulated IMLEO numbers for specific cycler orbits are printed in the paper; results
  are presented in graphs only.

### Semi-cycler architecture definition (Table 1, p.3)

- **E-M Semi-Cycler**: cycler vehicle makes repeated Earth-to-Mars transfers; crew
  returns to Earth directly (not on the cycler).
- **M-E Semi-Cycler**: cycler vehicle makes repeated Mars-to-Earth transfers; crew
  departs Earth on standard trajectory to meet cycler near Mars.

### References of interest extracted from pp. 31–32

Selected cycler-relevant entries from the reference list:

- [40] Niehoff, Friedlander, McAdams, "Earth-Mars Transport Cycler Concepts," IAF Paper
  91-438, IAC Montreal, Oct 1991.
- [41] Penzo, Nock, "Earth-Mars Transportation Using Stop-Over Cyclers," AIAA 2002-4424.
- [42] Aldrin, Byrnes, Jones, Davis, "Evolutionary Space Transportation Plan for Mars Cycling
  Concepts," AIAA 2001-4677, Aug 2001.
- [43] Landau, Longuski, "Mars Exploration via Earth-Mars Semi-Cyclers," AAS 05-269,
  Aug 2005.
- [45] McConaghy, Landau, Yam, Longuski, "A Notable Two-Synodic-Period Earth-Mars
  Cycler," JSR Vol. 43, No. 2, pp. 456–465, Mar–Apr 2006.
  [This is the canonical S1L1 reference.]
- [50] Landau, Longuski, "A Reassessment of Trajectory Options for Human Missions to
  Mars," AIAA 2004-5095, Aug 2004.
- [51] Landau, Longuski, "Method for Parking-Orbit Reorientation for Human Missions to
  Mars," JSR Vol. 42, No. 3, pp. 517–522, May–Jun 2005.
- [52] Friedlander, Niehoff, Byrnes, Longuski, "Circulating Transportation Orbits Between
  Earth and Mars," AIAA 1986-2009, Aug 1986. [The original cycler paper.]
- [53] Penzo, Nock, "Hyperbolic Rendezvous for Earth-Mars Cycler Missions," AAS 02-162,
  pp. 763–772, Jan 2002.
- [54] Zola, "A Method for Approximating Propellant Requirements of Low-Thrust
  Trajectories," NASA TN-D-3400, 1966. [IMLEO model used throughout this paper.]

### S1L1 relevance

Reference [45] (McConaghy et al. 2006) is the canonical S1L1 source; it is cited but not
summarised — no orbital elements or per-leg V∞ values are reproduced in this paper.

---

## Paper 3: Howe 2025

### What the paper is

A design paper for a **physical Mars Cycler habitat structure**: conceptual design, mass
breakdown, assembly concept, rotating artificial-gravity torus, propulsion requirements.
The trajectory basis is the up/down escalator cycler from Rauwolf, Friedlander, Nock 2002,
cited throughout.

### Trajectory data (p.4, citing Rauwolf, Friedlander, Nock 2002)

From Figure 3 caption and surrounding text (p.4):

| Trajectory | Departure day | Arrival day | ToF (days) |
|---|---|---|---|
| Up escalator (E→M) | Day 1 | Day 151 | 151 |
| Down escalator (E→M) | Day 38 | Day 229 | 191 [see note] |

**NOTE on down escalator ToF**: The text states "the down escalator Cycler starts on Day 38
and arrives at Mars on Day 229, for a total duration of 170 days" (p.4). Arithmetic
229 − 38 = 191 days; the paper prints "170 days." This is an internal inconsistency in Howe
2025 — report the paper's stated value (170 days) and flag the discrepancy. The day numbers
(38 and 229) are explicitly printed; "170 days" is also explicitly printed. **The correct ToF
from day numbers is 191 days; 170 days appears to be an error in Howe 2025.**

After Mars encounter, up escalator takes > 6 months before next useful Earth encounter (p.4).

**Source for these numbers**: Rauwolf, G. A., Friedlander, A., and Nock, K., AIAA paper
number not printed in Howe 2025 but the citation is "[2] Rauwolf, Friedlander, Nock 2002."

### V∞ constraint (p.5)

- Maximum V∞ for rendezvous: **< 6,000 m/s** (citing Rauwolf et al. 2002) — stated as
  design constraint, not a computed value for a specific encounter.

### ΔV budget (pp. 4–5)

| Maneuver | Value | Notes |
|---|---|---|
| ΔV from Earth-Mars L1 to cycler trajectory | 2,500 m/s | p.4 |
| ΔV to achieve cycler orbit (LEO injection) | 4,739 m/s | p.5, computed from Earth GM, R_Earth = 6,378,000 m, LEO alt = 300,000 m |
| Total ΔV for chemical propulsion to establish cycler | ~5,000 m/s | p.7 |

### Launch mass (p.7)

| Configuration | Launch requirement |
|---|---|
| Minimal cycler (dry ~467,990 kg) | 2 × Starship 2 (150 t capacity) + 3 × Starship 3 propellant |
| Cruise-ship cycler (~1,000 persons) | 15,000 t propellant = 83 × Starship 3 launches; 50–75 Starship launches for propellant alone |

### Structural dimensions (pp. 5, 12)

| Feature | Minimal cycler | Cruise-ship cycler |
|---|---|---|
| Outer (non-rotating) torus major circumference | 200 m (24-cassette ring) | larger |
| Inner (pressurized, rotating) torus minor circumference | 20 m (22 panels) | larger |
| Inner torus rotation rate | 4.62 rpm | 5.79 rpm (counter-rotating flywheel) |

### S1L1 relevance

None directly. The escalator cycler cited (Rauwolf 2002) is a different family from S1L1
(McConaghy 2006). The V∞ < 6 km/s design constraint is consistent with S1L1 encounter
speeds but does not resolve S1L1 orbital elements or real-ephemeris values.

---

## Section: S1L1 Resolution

**None of the three papers resolves or confirms the S1L1 real-ephemeris target numbers.**

- Agrawal 2022: no cycler trajectory data at all.
- Landau-Longuski 2009: cites McConaghy et al. 2006 (JSR 43(2):456–465) as ref [45] but
  does not reproduce its orbital elements, V∞ values, or ToFs. The paper's IMLEO
  comparisons treat cycler performance parametrically.
- Howe 2025: uses escalator cycler (Rauwolf 2002), not S1L1. No S1L1 data present.

The primary source for S1L1 orbital parameters remains McConaghy et al. 2006 (JSR).

---

## Candidate Catalogue Rows

**From Howe 2025, citing Rauwolf, Friedlander, Nock 2002:**

These are escalator cyclers (distinct from S1L1):

```
up_escalator_E_to_M:
  E_departure_day: 1          # relative day, single synodic cycle
  M_arrival_day:  151
  tof_days:       151
  source: Howe-2025 p.4, citing Rauwolf-Friedlander-Nock-2002
  caution: patched-conic/impulsive model; real-ephemeris values may differ

down_escalator_E_to_M:
  E_departure_day: 38
  M_arrival_day:   229
  tof_days_stated: 170        # value printed in Howe-2025
  tof_days_computed: 191      # 229-38; Howe-2025 internal inconsistency
  source: Howe-2025 p.4, citing Rauwolf-Friedlander-Nock-2002
  caution: day numbers and stated ToF disagree; do not use without checking primary source
```

No V∞ values at specific encounters are given for either escalator leg.

**No candidate catalogue rows from Agrawal 2022 or Landau-Longuski 2009.**

---

## Extraction Failures and Caveats

| Issue | Detail |
|---|---|
| Howe 2025 text layer | Scrambled/garbage characters in pdftotext; all 15 pages read as PNG via pdftoppm |
| Landau-Longuski 2009 text layer | No text layer (scanned document); all pages read as PNG via pdftoppm |
| Landau-Longuski 2009 pp.31–32 | PNG rendering at 200/150/100 DPI failed API image-size limit; read at 80 DPI as JPEG — legible, references extracted successfully |
| Landau-Longuski 2009 exact venue | Publication year 2009 is from Howe-2025 citation; PDF metadata shows 2006 creation date; precise journal volume/issue not printed in the scanned document |
| IMLEO figure values | Landau-Longuski 2009 IMLEO comparisons are in graphs, not tables; no tabulated numbers could be extracted without digitising the figures |
| Rauwolf 2002 primary source | Not in the PDF set; escalator values above are second-hand from Howe 2025 |

---

## Summary

| Paper | Biggest finding | Cycler orbital data? | S1L1 resolved? | Catalogue rows |
|---|---|---|---|---|
| Agrawal 2022 | Mars Spacedock logistics node; ~60 Mg NRHO mass saving with refuel node; cyclers only mentioned in future-work | No | No | 0 |
| Landau-Longuski 2009 | E-M semi-cycler ranks #1 for IMLEO under advanced propulsion; cites McConaghy 2006 (S1L1 source) without reproducing data | No (IMLEO graphs only) | No | 0 |
| Howe 2025 | Up escalator ToF 151 d, down escalator ~170/191 d (inconsistency); V∞ < 6 km/s design constraint; cycler structural design details | ToF from Rauwolf 2002 only; no V∞ by encounter | No | 2 (with caveats) |
