# Digest: Lam, Buffington & Campagnola 2018, "A Robust Mission Tour for NASA's Planned Europa Clipper Mission" (#960 batch 28; #943)

T. Lam, B. Buffington and S. Campagnola (JPL), AIAA 2018-0202 (SciTech 2018), doi 10.2514/6.2018-0202.
- Filed as `cyclers_pdf/papers/lam-buffington-campagnola-2018-robust-mission-tour-europa-clipper-aiaa-2018-0202-doi-10.2514-6.2018-0202.pdf`.
  21 pp., text layer, md5 cfe25e50603bc353b6e7156f5dc7868a. Supplied by the owner.
- I read the tour descriptions and the full flyby tables (Tables 2-3, all rows), and skimmed the
  requirements section.

## 0. Verdict (sent to twobody-gen2-opus)

**No collision with `#943` gc-1, gc-2 or ge-1..3. There is no repeating G-C or G-E segment.**
- Tours 17F12 (46 Europa, 4 Ganymede and 9 Callisto flybys over 3.7 yr) and 17F13 (54 Europa, 4
  Ganymede and 9 Callisto) have the same structure: **G^4 E^n C^8 E^n C**.
- **Ganymede is used only in the pump-down:**
  - G00 at JOI (v_inf 8.45 km/s);
  - 01G01-06G04 (v_inf 6.77-7.03 km/s);
  - then a single G -> E handoff (07E01 at 4.04 km/s).
- **Callisto is used in one contiguous "switch-flip" block of 8 flybys:** in 17F12, 33C01-45C08 at
  v_inf 2.73-3.63 km/s, entered from 32E24 and left to 47E25 by E->C and C->E pi-transfers. In 17F13 the
  block is 39C01-56C08. Disposal is a Callisto impact (C09).
- This answers the 21F31 gap: these tours state the transition types, and **none is a G-C alternation.**

## 1. Content (READ)

- Tour phases: capture and pump-down (robust to JOI anomalies, after Scott et al. AAS 17-437), COT-1 and
  COT-2 (crank-over-the-top Europa sequences on the anti-Jovian hemisphere), a non-resonant transfer, a
  petal rotation, the illuminated-hemisphere transition, the Callisto switch-flip, then COT-3 to COT-5.
- 17F13 variants let COT-3 repeat COT-1 or COT-2 groundtracks (Tables 4-5, with dV costs).
- Table 6 compares the tours from JUP230 to 17F13.
- Radiation dose to the last Europa flyby is given (17F13: 2.98 Mrad).

## 2. Citation mining

- Held: Buffington, Campagnola & Petropoulos 2012 [4]; Lam, Arrieta-Camacho & Buffington 2015 [5];
  Anderson, Buffington & Campagnola 2014/2018 [15] (the JGCD form is held); Lantukh & Russell 2012
  [17].
- Already on the wanted list (Clipper background row): Buffington 2014 [6]; Uphoff, Roberts & Friedman
  1976 [11].
- Not held, low priority, not added: Johannesen & D'Amario 1999 [1], Buffington, Strange & Campagnola
  ISSFD 2012 [12], Scott et al. AAS 17-437 [13], Smith AAS 98-106 [16], and project and radiation
  reports.
