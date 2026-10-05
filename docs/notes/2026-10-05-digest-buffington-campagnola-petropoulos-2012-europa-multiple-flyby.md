# Digest: Buffington, Campagnola & Petropoulos 2012, "Europa Multiple-Flyby Trajectory Design" (#960)

B. Buffington (JPL), S. Campagnola (JAXA/ISAS) & A. Petropoulos (JPL), "Europa Multiple-Flyby Trajectory
Design", AIAA 2012-5069, AIAA/AAS Astrodynamics Specialist Conference, 2012, doi 10.2514/6.2012-5069
(wanted-list CONFIRMED).
- Filed as `cyclers_pdf/papers/buffington-campagnola-petropoulos-2012-europa-multiple-flyby-trajectory-design-aiaa-2012-5069-doi-10.2514-6.2012-5069.pdf`.
  20 pages. The upload (md5 34f489483cde5cb94af44500c713a235) was image-only. I OCR'd it with
  `ocrmypdf --skip-text` (10.8k words); the filed PDF carries that layer. Formerly part of wanted-list row 58.
- The OCR of the tables is unreliable. Every number below was read on the page images (pp.1, 5-10, 12,
  15, 16, 20).

## 0. Verdict

The 11-F5 tour is a one-shot Jovian tour: 34 Europa and 9 Ganymede flybys over 2.4 years, ending in a
Ganymede impact. Its one Ganymede-Europa segment is the "switch-flip": E23, then a Europa-to-Ganymede
pi-transfer to G24, a 3.5-d Ganymede pi-transfer to G25, a 1:1 Ganymede transfer to G26, then a
non-resonant transfer to E27. It happens once and moves the Europa flybys 180 deg around Jupiter.
- **`#943` X1 Ganymede-Europa:** no literal collision. There is no repeating G-E path.
- **`#943` X1 Ganymede-Callisto:** no collision. There are no Callisto flybys at all.

It is the source of the crank-over-the-top (COT) definitions and of a full, published flyby list
(Table 4) with dates, altitudes, v_inf and resonances.

## 1. Content (READ)

- **Baseline (pp.1, 4-7):**
  - Atlas V 551. Launch 21 Nov 2021 (21-d period, C3 15.0 km^2/s^2). VEEGA: Venus 14 May 2022; Earth
    24 Oct 2023 and 20 Oct 2025.
  - G0 at 500 km on 3 Apr 2028. JOI of 857 m/s at 12.8 Rj on 4 Apr 2028 (Table 1).
  - Tour: maximum inclination 14.9 deg (abstract: 15 deg), deterministic ΔV 157 m/s after PJR, TID
    2.0 Mrad (abstract: 2.06).
- **Phases (Table 3, p.7):**
  - Pump-down: G1-G4, Apr 2028-Feb 2029.
  - COT-1: 7 Europa flybys at the ascending node, anti-Jovian hemisphere.
  - Non-resonant transfer.
  - COT-2: 6 flybys at the descending node.
  - Pump-down / crank-up, then the switch-flip (Jan-Feb 2030).
  - COT-3: 8 flybys.
  - Non-resonant transfer.
  - COT-4: 6 flybys.
  - Ganymede impact on 25 Aug 2030 (G42).
- **COT definition (p.10, sec. III.F):**
  - Start from an equatorial orbit, crank the inclination up to i_max, then back to the equator, using
    resonant transfers.
  - Starting from an inbound flyby, the sequence changes the flybys to outbound "over the top", and the
    reverse.
  - With the same resonance throughout (pure cranking), the closest approaches lie near the prime or
    180-deg meridians.
  - i_max depends on the orbit period and v_inf. When the orbit period exceeds the moon's, i_max is
    reached with the moon at the spacecraft's periapsis.
- **Pi-transfer (p.12 footnote):** a non-resonant, typically inclined transfer whose two flybys are n*pi
  apart in true anomaly (n odd), on opposite sides of Jupiter.
- **Lighting change options (p.12):**
  1. non-resonant Callisto and/or Ganymede transfers (longest, lowest dose);
  2. non-resonant Europa transfers only (highest dose);
  3. the switch-flip (fastest; used).
- **Table 4 (pp.8-9) Ganymede segment, read on the image:**

| Flyby | In/Out | Date | Alt, km | v_inf, km/s | Inc, deg | Resonance m:n | Period, d |
|---|---|---|---|---|---|---|---|
| Europa23 | I | 28 Dec 2029 08:03:48 | 805.1 | 3.89 | 13.87 | NR | 5.22 |
| Ganymede24 | O | 26 Jan 2030 00:58:15 | 1346.7 | 2.78 | 14.86 | pi-transfer | 7.15 |
| Ganymede25 | O | 29 Jan 2030 13:39:45 | 123.1 | 2.79 | 11.93 | 1:1 | 7.11 |
| Ganymede26 | O | 05 Feb 2030 16:18:38 | 1584.7 | 2.75 | 10.2 | NR | 5.39 |
| Europa27 | I | 14 Feb 2030 04:41:04 | 100 | 3.51 | 10.81 | 5:3 | 5.92 |

  The pump-down Ganymede v_inf are G0 7.382, G1 6.34, G2 6.42, G3 6.37 and G4 6.40 km/s. Europa v_inf
  falls from about 3.9 km/s (COT-1/2) to about 3.5 km/s (COT-3/4).
- **Navigation (p.15):**
  - First G and E flybys at 500 and 724 km; ratchet down to 100 km and 25 km.
  - COT-1 alternates 4:1 (14.2 d) and 7:2 (24.88 d); COT-2 is five back-to-back 4:1; COT-3 alternates
    3:1 (10.65 d) and 5:2 (25.44 d); COT-4 is five back-to-back 3:1.
  - The 3.5-d G-G pi-transfer is ballistic. If it proves too aggressive, a 3-, 5- or 7-pi-transfer
    (10.5, 14 or 17.5 d) can replace it.
- **ΔV (Table 5, p.16):** CBE 1311 m/s, MEV 1675 m/s. Tour deterministic 157/200; JOI 857/900;
  PJR 114/135.

## 2. Gate relevance

- **`#943` X1 Ganymede-Europa:** no collision.
  - The switch-flip is a published ballistic E-G-G-G-E chain. A G-E cycler result should cite it as the
    nearest operational relative (INFERRED).
  - The 3.5-d Ganymede pi-transfer is half the Ganymede period (7.155 d / 2 = 3.58 d), as a pi-transfer
    should be.
- **`#943` X1 Ganymede-Callisto:** not touched (no Callisto flybys).
- **Cite with:** McElrath, Campagnola & Strange 2012 (filed today), Campagnola et al. 2014 and 2019
  (held), and Lam, Buffington & Campagnola 2018 (in the wanted list, Clipper and Jovian background row).

## 3. Positive controls

- Table 4: the full 11-F5 flyby list (dates to the second, altitudes, B-plane angles, v_inf,
  inclinations, perijove and apojove, resonances). Usable as a one-shot real-ephemeris tour check. Read
  every value on the page image; the OCR layer garbles the table.

## 4. Citation mining (references 1-13, p.20)

Held:
- none.

Not held:
1. Kloster, K. W., Petropoulos, A. E. & Longuski, J. M. (2011), "Europa Orbiter Tour Design with Io
   Gravity Assists", Acta Astronautica 68:931-946, doi 10.1016/j.actaastro.2010.08.041 (printed).
2. Kahn, M., Campagnola, S. & Croon, M. (2004), "End-to-End Mission Analysis for a Low-Cost,
   Two-Spacecraft Mission to Europa", Adv. Astronaut. Sci. 119:463-472. No DOI.
3. Boutonnet, A. & Schoenmaekers, J. (2012), "Mission Analysis for the JUICE Mission", AAS 12-207. No
   DOI.
4. Johannesen & D'Amario (1999), AAS 99-360; Petropoulos, Longuski & Bonfiglio (2000), JSR 37(6).
   Both are also cited by McElrath 2012.
5. JPL study reports (Europa Explorer 2007, JEO 2008, Europa Study 2012) and science background
   (Khurana 1998, Anderson 1998, Schubert 2009, Decadal Survey 2011). Background.
