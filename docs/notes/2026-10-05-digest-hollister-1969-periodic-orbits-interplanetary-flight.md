# Digest: Hollister 1969, "Periodic Orbits for Interplanetary Flight" (#960)

W. M. Hollister (MIT), "Periodic Orbits for Interplanetary Flight", J. Spacecraft and Rockets 6(4):366-369
(April 1969), doi 10.2514/3.29664.
- Crossref-confirmed in batch 1. Presented as AAS 68-102, September 1968.
- Filed as `cyclers_pdf/papers/hollister-1969-periodic-orbits-interplanetary-flight-JSR-6-4-366-doi-10.2514-3.29664.pdf`.
  4 pages, text layer, md5 32233e45afefa46b04b5a78a0f9e6ac5.
- I read all pages from the text layer. Tables 1-2 (pp.368-369) were read as page images.

## 0. Verdict

This is the original published record of the Earth-Venus cyclers ("periodic orbits connecting Earth and
Venus").
- It defines them: free-fall, flybys at both ends, "back and forth between Earth and Venus forever".
- It finds three orbits, I, II and III, first in the circular-coplanar model (3.2-yr repeat) and then in
  the real inclined-elliptic model (16-yr repeat).
- Menning 1968 and Hollister & Menning 1970 extend these to 15 orbits. Orbits I, II and III are Menning's
  1H, 2H and 3H.
- There is NO Earth-Mars or Venus-Mars periodic orbit. Mars is discussed only as future work (p.369). The
  reason given: the limiting turn at Mars "is half as large as it is at Earth or Venus" (Fig. 4), and a
  Lambert solution "might call for the spaceship to go beneath the surface during a Mars flyby".
- Gate for `#942`:
  - R1(b): published prior art at the origin. Same itineraries as the 15 H&M rows.
  - R1(a) and (c): no collision. These are stated as open ("unending orbits ... which periodically visit
    Mars as well" are suggested but not found).

## 1. Content (READ)

- Repeat geometry (p.366):
  - Earth 32 revolutions, Venus 52, Mars 17 in 32 yr.
  - Earth-Venus absolute orientation repeats every 8 yr (8 and 13 revolutions), with 5 alignments.
  - Synodic periods: Earth-Venus 1.6 yr, Venus-Mars 32/35 yr, Earth-Mars 32/15 yr. In the circular model
    the pattern repeats every 6.4 yr.
- Direct-return types (p.367):
  - "full-revolution return": the same period as the planet, a double infinity of them. The
    "half-revolution return" is a special inclined case.
  - "symmetric return": coplanar, about 1.41 solar revolutions (after Ross; Fig. 3).
  - Both keep the same |V_inf| at arrival as at departure.
- **Circular-coplanar case (p.367):** periodic if it repeats after a multiple of 1.6 yr. The combinations
  found total 4.2 solar revolutions in 3.2 yr.
  - **Orbit I:** one full-revolution return at Earth, a transfer to Venus, two full-revolution returns at
    Venus, and a transfer back. Each transfer is 0.485 yr and 0.6 revolution. V_inf is 0.126 EMOS (3.8
    km/s) at Venus and 0.107 EMOS (3.2 km/s) at Earth.
  - **Orbit II:** Earth full-revolution; at Venus one symmetric and one full-revolution return.
  - **Orbit III:** Earth symmetric; Venus two full-revolution returns.
  - Flybys must stay above 1.1 planet radii (Fig. 4: turn limits at 1.1 Earth/Venus radii and 1.3 Mars
    radii).
- **General inclined-elliptic case (p.368):**
  - 16-yr repeat with 10 Earth-Venus transfers.
  - A 10-dimensional iteration (10 Lambert solutions) on the transfer endpoint times t_i to null the
    inbound-outbound speed differences V_i. Steepest descent worked; Newton also worked but converged less
    well.
  - For orbits II and III the symmetric-return duration was held constant and equal speeds were assumed
    there.

**Table 1, periodic orbit I, general case (READ p.368, image).** Columns: event (LV = launch, AR =
arrival), body, JD - 2440000, V_inf (EMOS), angle (deg; in the orbital plane, clockwise from the
circumferential direction), elevation (deg; positive above the orbital plane).

| Event | JD-2440000 | V_inf (EMOS) | Ang | Elev |
|---|---|---|---|---|
| LV E | 0806 | 0.154 | 163 | -52 |
| AR V | 0971 | 0.178 | 38 | 55 |
| LV V | 1421 | 0.178 | 335 | -60 |
| AR E | 1592 | 0.154 | 201 | 50 |
| LV E | 1957 | 0.154 | 143 | -44 |
| AR V | 2125 | 0.203 | 8 | 67 |
| LV V | 2575 | 0.203 | 294 | 14 |
| AR E | 2797 | 0.191 | 240 | 12 |
| LV E | 3163 | 0.191 | 186 | 58 |
| AR V | 3316 | 0.194 | 36 | -62 |
| LV V | 3765 | 0.194 | 343 | 65 |
| AR E | 3935 | 0.156 | 215 | 46 |
| LV E | 4300 | 0.156 | 154 | 50 |
| AR V | 4471 | 0.192 | 17 | -64 |
| LV V | 4921 | 0.192 | 319 | 58 |
| AR E | 5077 | 0.174 | 183 | -56 |
| LV E | 5442 | 0.174 | 123 | -12 |
| AR V | 5664 | 0.223 | 69 | -15 |
| LV V | 6114 | 0.223 | 8 | -69 |
| AR E | 6285 | 0.154 | 215 | 47 |
| LV E | 6650 | 0.154 | 163 | -52 |

"Repeating after 16 yr."

Table 1 is the only published set of v_inf DIRECTIONS for an Earth-Venus cycler. The directions are needed
to reconstruct the full-revolution return legs.

**Table 2 (READ p.369, image).** V_inf in EMOS x 10^3, lowest/average/highest:

| Orbit | Earth | Venus |
|---|---|---|
| I | 154/165/191 | 178/198/223 |
| II | 136/145/159 | 232/243/249 |
| III | 173/196/214 | 144/174/232 |

**Applications (p.369):**
- Each orbit is two synodic periods, so two spacecraft are needed for every opportunity, giving a Venus
  flyby every 6.4 months on average.
- "a minimum of twenty combinations to investigate", which Menning later counted as at least 1024.
- Suggested uses: shelter, communications link, rescue station, navigation beacon.

## 2. Cross-checks against held sources (COMPUTED, comparison)

- Orbit I dates (806 ... 6650) equal Menning's 1H dates (thesis p.A-2: 441 806 971 1196 1421 1592 ...). Table
  1 omits the full-revolution intermediate dates.
- Menning's 1H closes at 6285 and Hollister lists 6285 (AR E) then 6650 (LV E). These agree, since 6285 + 365
  = 6650 is the next full-revolution return.
- Menning orbit 1 (the converged rigorous version) has V_r 0.155 at 441/806, against Hollister's 0.154 at
  806. Menning orbit 1 is the re-converged solution, so the last-digit difference is expected.
- The H&M 1970 Table 3 orbit 1 matches Menning orbit 1.

## 3. Positive controls

- Circular-coplanar orbit I: 0.485 yr and 0.6 revolution transfers; V_inf 0.126 EMOS at Venus and 0.107 at
  Earth. A clean analytic control for a circular-coplanar Earth-Venus generator (Earth period 1 yr, Venus
  8/13 yr).
- Table 1: the general-case orbit I with v_inf directions, the only directional control.
- Table 2: V_inf ranges for I, II and III.

## 4. Citation mining (references 1-8, p.369)

Held:
- [1] Breakwell & Perko 1965: its results are in the held Breakwell & Perko 1974; AIAA 65-689 itself is not
  held.

Not held, in priority order:
1. [4] Hollister, W. M. (1967), "Periodic Orbits Connecting Earth and Venus", ZAMM 47, Sonderheft. DOI not
   checked (GAMM proceedings). The earliest announcement.
2. [6] Ross, S. (1963), "A Systematic Approach to the Study of Nonstop Interplanetary Round Trips", AAS 9th
   Annual Meeting. No DOI. Source of the symmetric-return idea and the E-V-M-E round trips.
3. [2] Gillespie & Ross 1967, JSR 4(2):170-175. Also cited by VanderVeen and Menning.
4. [8] Hollister, W. M. & Prussing, J. E. (1966), "Optimum Transfer to Mars via Venus", Astronautica Acta
   12(2). DOI not checked.
5. [3] Lockheed 1962, 3-17-62-1; [5] Deerwester 1966, JSR 3(10):1564; [7] Casal & Ross 1965 (AAS
   symposium).
