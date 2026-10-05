# Digest: Lam, Arrieta-Camacho & Buffington 2015, "The Europa Mission: Multiple Europa Flyby Trajectory Design Trades and Challenges" (#960)

T. Lam, J. J. Arrieta-Camacho & B. B. Buffington (JPL), "The Europa Mission: Multiple Europa Flyby
Trajectory Design Trades and Challenges", AAS 15-657, AAS/AIAA Astrodynamics Specialist Conference, 2015.
No DOI.
- Filed as `cyclers_pdf/papers/lam-arrieta-camacho-buffington-2015-europa-mission-multiple-europa-flyby-trajectory-design-trades-AAS-15-657.pdf`.
  20 pages, text layer, md5 c746297062f7b9bcbd2698e514512d78. Former wanted-list row 3.
- I read all pages from the text layer. The text layer splits the tour tables into columns, so I read
  Tables A1 and A2 (pp.18-19) on the page images. Table 2 (p.8) was read from the text layer only.

## 0. Verdict

The paper gives four full Europa Clipper candidate tours (13F7-A21, 14F8-S22, 15F9-A22, 15F10-S22) with
every flyby. None has a repeating Ganymede-Callisto or Ganymede-Europa path. So there is **no literal
collision for `#943` X1**, and both verdicts stand:
- **X1 Ganymede-Callisto stays PARTIAL** (from the Campagnola 2019 GCGC cycler).
  - The only G-C alternation is the 13F7-A21 pump-down, G3-C1-G4-C2 (Table A1).
  - Its transfers are non-resonant (m:n of 0.9:0.8, 3.8:1.1, 1:1.1), and v_inf falls at each moon:
    G 6.84 to 5.81 km/s, C 5.92 to 4.67 km/s. So it is energy pumping, not a cycler.
  - It is an earlier published G-C-G-C flyby chain. A G-C cycler result should cite it as a precursor
    (INFERRED).
- **X1 Ganymede-Europa stays OPEN.** No tour has a G-E chain. Ganymede appears only in the pump-down and
  at end of mission.

**What "15F09" was.** Campagnola et al. 2019 say their GCGC cycler "is inspired by the Clipper tour 15F09".
By the naming rule (p.5: YY T N - V L), 15F9-A22 is the ninth flyby tour, designed in 2015 (INFERRED to be
the same tour).
- Its 180-deg "switch flip" uses Callisto only: E23, then C1-C7, then E24, with a Callisto pi-transfer
  C4-C5 (Table A2).
- There is no Ganymede in that rotation, so 15F09 inspired the fast rotation, not the G-C pattern. This
  closes the open follow-up in the Campagnola 2019 digest.

## 1. Content (READ)

- **Stages (pp.3-4):** pump-down (Ganymede and/or Callisto), COT-1, COT-2, petal rotation, switch flip
  ("resonant and pi-transfers of Europa and Callisto"), COT-3, COT-4, and impact on Ganymede or Callisto.
  "Flybys of Ganymede and Callisto are only used for gravity assist purposes."
- **Naming (p.5 footnote):** YYTN-VL.
  - YY: design year. T: type (F = flyby). N: sequence number.
  - V: launch vehicle (A = Atlas V 551, S = SLS, D = Delta IV H). L: launch year.
- **Baseline 14F8-S22 (pp.5-8, Table 2):**
  - SLS direct launch, June 2022. Arrival 24 Dec 2024 or 5 Mar 2025.
  - Tour: 45 E, 5 G and 7 C flybys over 3.68 yr; maximum inclination 21.2 deg; deterministic ΔV
    174 m/s; TID 3.14 Mrad.
  - Petal rotation uses alternating 4:1+ and 5:1- non-resonant transfers (Europa periods 14.6 and
    17.6 d).
  - Switch flip (text layer): E27, then C1 (8 Jun 2027) to C7 (30 Nov 2027), then E28. Callisto v_inf
    2.92-2.98 km/s. C3 to C4 is half a Callisto period ("0.5:0.5", 8.3 d).
  - There is also one non-targeted Ganymede passage, "32G (NT)", at 31,470 km during the petal rotation.
- **Table 3 (p.13), tour comparison:**

| Tour | Launch | Cruise | Tour | E/G/C flybys | Tour ΔV | Total ΔV | i_max | TID |
|---|---|---|---|---|---|---|---|---|
| 13F7-A21 | Nov 2021 | 6.4 yr | 3.5 yr | 45/5/9 | 164 m/s | 1218 m/s | 20.1 deg | 2.82 Mrad |
| 14F8-S22 | Jun 2022 | 2.7 yr | 3.7 yr | 45/5/7 | 174 m/s | 1181 m/s | 21.2 deg | 3.14 Mrad |
| 15F9-A22 | Jun 2022 | 7.6 yr | 3.3 yr | 40/5/8 | 152 m/s | 1092 m/s | 19.8 deg | 2.74 Mrad |
| 15F10-S22 | Jun 2022 | 2.7 yr | 3.4 yr | 42/4/8 | 118.2 m/s | 1093 m/s | 21.2 deg | 2.82 Mrad |

  Count note: Table A2 lists six Ganymede flybys for 15F9-A22 (G0-G5) against 5 in Table 3. The
  Europa and Callisto counts match the appendix tables.
- **13F7-A21 pump-down (Table A1, page image):**

| Flyby | In/Out | Date (ET) | Alt, km | v_inf, km/s | m | n | Next enc., d |
|---|---|---|---|---|---|---|---|
| 0G0 | I | 03-Apr-2028 11:58 | 500 | 8.36 | 28.2 | 1 | 202.2 |
| 1G1 | O | 22-Oct-2028 17:15 | 100 | 6.93 | 8 | 1 | 57.2 |
| 2G2 | O | 18-Dec-2028 23:07 | 100 | 6.98 | 4 | 1 | 28.6 |
| 3G3 | O | 16-Jan-2029 13:49 | 1035 | 6.84 | 0.9 | 0.8 | 14.9 |
| 4C1 | I | 31-Jan-2029 12:30 | 912 | 5.92 | 3.8 | 1.1 | 27.3 |
| 5G4 | O | 27-Feb-2029 18:43 | 524 | 5.81 | 1 | 1.1 | 16.5 |
| 6C2 | O | 16-Mar-2029 06:22 | 2636 | 4.67 | 2.7 | 0.8 | 9.7 |

  m = moon orbits, n = spacecraft orbits. Non-integer values mark non-resonant transfers.
- **15F9-A22 switch flip (Table A2, page image):**
  - 32E23 (O, 01-Feb-2032).
  - 37C1 (I, 14-Mar-2032, 300 km, v_inf 2.77), 2:3, 33.3 d.
  - 40C2 (I, 272 km), 1:1, 16.6 d.
  - 41C3 (I, 50 km), 1:1, 16.7 d.
  - 42C4 (I, 1826 km), 0.5:0.5, 8.3 d: the Callisto pi-transfer.
  - 42C5 (O, 554 km), 1:1, 16.6 d.
  - 43C6 (O, 154 km), 2:3, 33.3 d.
  - 46C7 (O, 17-Jul-2032, 1653 km, v_inf 2.75), 1:1.5, 15.2 d.
  - 48E24 (I, 01-Aug-2032).
  - Check: the pi-transfer's 8.3 d is half the Callisto period (16.689 d / 2 = 8.34 d).
- **Printed date slip (Table A2, INFERRED):** flybys 26E19-29E22 are printed in 2032 (10-Nov, 25-Nov,
  13-Dec, 27-Dec-2032). The next-encounter column puts them in 2031: 24-Oct-2031 + 17.6 d = 10-Nov-2031,
  and 29E22 + 35.6 d = 32E23 on 01-Feb-2032. Use the 2031 dates.
- **Robustness (pp.13-15):**
  - Decision-tree tours, designed in advance.
  - "Reverse-crank": reverse the next flyby's crank to repeat a missed groundtrack, at about 25-30 m/s
    (Table 4).
  - Solar-conjunction rules: no maneuvers below SEP 5 deg; resume at about SEP 12 deg. These force empty
    orbits, for example the 11:3 resonance after 14F8 E12.

## 2. Gate relevance

- **`#943` X1 Ganymede-Callisto:** PARTIAL stands; no collision. 13F7-A21 G3-C1-G4-C2 is the earliest
  published G-C-G-C flyby chain known to us. It is non-resonant and pumping, not periodic.
- **`#943` X1 Ganymede-Europa:** OPEN stands. There are no G-E chains.
- **Cite with:** Campagnola et al. 2019 (GCGC; held), Buffington, Campagnola & Petropoulos 2012 (11-F5,
  the Ganymede switch-flip; filed today), and Anderson, Campagnola & Buffington 2018 (petal rotation;
  held).

## 3. Positive controls

- Tables 2 and A1-A3: complete flyby lists (dates to the minute, altitudes, B-plane angles, v_inf,
  inclinations, m:n, next-encounter times) for four tours. Read values on the page images; the text
  layer splits the columns.
- The 15F9-A22 Callisto pi-transfer (C4-C5, 8.3 d) is a direct check of the half-period rule.

## 4. Citation mining (references 1-13, p.16-17)

Held:
- [5] Buffington, Campagnola & Petropoulos 2012 (filed today).
- [12] Anderson, Buffington & Campagnola, AIAA 2014-4350: held as the 2018 JGCD version
  (anderson-campagnola-buffington-2018-...).

Not held:
1. Buffington, B., Strange, N. & Campagnola, S. (2012), "Global Moon Coverage Via Hyperbolic Flybys". The
   COT reference. No venue or DOI given.
2. Buffington, B. (2014), "Trajectory Design for the Europa Clipper Mission Concept", AIAA 2014-4105.
   Already in the wanted list (Clipper and Jovian background row).
3. Wu, X. et al. (2001), "Probing Europa's hidden ocean from tidal effects on orbital dynamics", GRL
   28(11):2245-2248, doi 10.1029/2000GL012814 (printed). Petal-rotation science background.
4. Johannesen & D'Amario 1999 (AAS 99-360); Clark 2007 (Europa Explorer report); Garrett 2002 and Divine
   & Garrett 1983 (radiation models); Pappalardo 2010; NRC 1999. Background.
