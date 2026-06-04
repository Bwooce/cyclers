# S1L1 Target Topology Mining

Sources mined:
- **Patel 2019**: Bhumika Patel, "Earth – Mars Cycler Vehicle Conceptual Design," M.S. Thesis, Florida Institute of Technology, December 2019. 73 pp. (FIT ETD repository)
- **Sanchez Net 2022**: Marc Sanchez Net et al., "Cycler Orbits and Solar System Pony Express," *Journal of Spacecraft and Rockets*, Vol. 59, No. 3, May–June 2022, pp. 861–870. DOI: 10.2514/1.A35091

---

## 1. S1L1 Target Resolution

### 1.1 What source does Patel 2019 use for S1L1?

Patel's reference [15] is:

> T. T. McConaghy, D. Landau, C. Hong Yam and J. M. Longuski, "Notable Two-Synodic-Period Earth-Mars Cycler," *Journal of Spacecraft and Rockets*, vol. 43, no. 2, pp. 456–465, March–April 2006.

This is cited explicitly as "McConaghy et al. [15]" throughout Chapter 1 and is the sole source for every S1L1 trajectory number in the thesis. No other trajectory paper (Russell 2004, Rogers 2012) is cited for the orbital parameters. [Patel p. 3–5, ref list p. 62–63]

### 1.2 Circular-coplanar or real-ephemeris?

Patel Section 1.3 (Limitations, p. 6) states explicitly:

> "It was assumed that the Earth and Mars have **circular coplanar orbits** and that the Earth-Mars synodic period (τ) is exactly 2 1/7 years."

The thesis operates entirely in the circular-coplanar model. No real-ephemeris computation is performed or cited. All values below derive from McConaghy 2006 in the coplanar model.

### 1.3 Table of S1L1 values from Patel (all sourced from McConaghy et al. [15])

| Parameter | Value stated in Patel | Page / Table |
|---|---|---|
| Semi-major axis | 1.945 × 10⁸ km (1.30 AU) | p. 5, Table 1 |
| Eccentricity | 0.257 | p. 5, Table 1 |
| Aphelion radius | 2.454 × 10⁸ km (1.64 AU) | p. 5, Table 1 |
| Perihelion radius | 1.451 × 10⁸ km (0.97 AU) | p. 5, Table 1 |
| Earth–Mars transfer time | ~154 days | p. 5, Table 1; confirmed p. 3, 6 |
| Repeat time (full cycle) | 4 2/7 years | p. 5, Table 1 |
| v∞ at Earth flyby | **3.657 km/s** | p. 5, Table 1 (row labelled "v∞ Flyby") |
| v∞ at Earth (launch) | **~2.492 km/s** (as "launch v∞") | p. 3 (text only, from Rogers [14] comparison) |
| v∞ at Earth (flyby, from Rogers comparison) | **~3.657 km/s** | p. 3 (text) |
| Direct Δv LEO→S1L1 (one impulsive burn) | 3.796 km/s | p. 5, Table 1; p. 7, ConOps |
| Return time if Mars arrival aborted | 2.8276 years (~1032 days) | p. 6 |
| Design lifetime (full reset period) | 30 years (~7 cycles of 4 2/7 yr) | p. 4 |

**Note on the v∞ row labelled "v∞ Flyby" (3.657 km/s):** Patel's Table 1 has a single v∞ row with value 3.657 km/s, labelled "v∞ Flyby." Context (p. 3) clarifies this is the flyby v∞ at Earth. No separate Mars v∞ figure appears anywhere in Patel. The paper gives no V∞ at Mars encounter.

**Additional v∞ values on p. 3 (from Rogers et al. [14] comparison, not McConaghy):**
- Aldrin cycler launch v∞: 3.449 km/s; flyby v∞: 6.546 km/s
- S1L1 launch v∞: 2.492 km/s; flyby v∞: 3.657 km/s

These are cited as Rogers [14] results, not Patel's own calculation.

### 1.4 Topology as described in Patel

Patel pp. 3–4 and Figure 1 (p. 4) describe the S1L1 topology as follows (all attributed to McConaghy [15]):

- The cycler **repeats every two synodic periods** (4 2/7 years = two × 2 1/7 yr).
- Within one full cycle there is a **short period leg** ending at t = 2.8276 years, when the cycler encounters Earth again after two solar revolutions (crossing Mars orbit four times).
- After that Earth flyby the cycler completes **one and a half more solar revolutions** and returns to Earth at t = 4 2/7 years (the **long period leg**).
- Figure 1 (p. 4) is a 2D diagram labelled "Ballistic S1L1 Cycler Trajectory [15]" showing two arcs with Earth encounters at t = 0, t = 2.8276 yr, t = 4 2/7 yr.
- Section 2.1 ConOps (p. 8) names the two legs:
  - **"Earth–Earth loop"** (short leg, only Earth encounters, used for resupply)
  - **"Earth–Mars loop"** (long leg, carries crew to Mars)
- The ConOps describes a **one-way crew trip of 154 days** from Earth to Mars on the long leg.
- **The cycler encounters Mars during the short leg** ("the cycler makes two solar revolutions while on short period leg which means it crosses Mars orbit four times," p. 3), but there is **no crew-carrying return leg directly from Mars back to Earth**; after the Mars encounter the cycler swings back toward Earth on the short leg. The Mars crossing on the short leg is not described as a usable transfer leg in Patel.
- Only Earth provides gravity-assist maneuvers; no Δv maneuver at Mars. [p. 6, Limitations]
- The first transfer opportunity (t = 0, Earth→Mars) and the fourth Mars crossing (t = 2.8276 yr, Mars→Earth direction for the cycler) both have a transfer time of **153.15 days** [p. 3, attributed to McConaghy [15]].

**Encounter sequence as Patel describes it:** E(t=0) → [Mars crossing on outbound arc, ~154 days] → E(t=2.8276 yr, short-leg Earth flyby) → E(t=4 2/7 yr, long-leg Earth return). Mars is crossed during the short-period leg. The crew transfer to Mars uses the arc from E(t=0) to Mars; crew return from Mars uses the arc from Mars back to E on the same short leg at t=2.8276 yr. This gives a roughly **E → M → E → E** structure per 4 2/7-year cycle, with both E→M and M→E legs of ~154 days.

### 1.5 Reconciliation with the three catalogue V∞ pairs

| Source | V∞ Earth (km/s) | V∞ Mars (km/s) | Model |
|---|---|---|---|
| Patel 2019 (= McConaghy 2006) | 3.657 (flyby) / 2.492 (launch) | **Not stated** | Circular coplanar |
| Russell 2004 (not in Patel) | 4.99 | 5.10 | Circular coplanar |
| McConaghy 2006 (= Patel [15]) | ~4.71 (E flyby, from p. 3 text) | — | Circular coplanar |
| Our catalogue "spec §9" | 5.65 | 3.05 | Unknown |
| McConaghy 2006 real-eph (catalogue "4.7/5.0") | 4.7 | 5.0 | Real ephemeris |

**Critical observations:**

1. Patel's Table 1 gives only a **single v∞ value of 3.657 km/s**, which is the Earth flyby v∞, not the Mars encounter v∞. Mars V∞ is absent from Patel entirely.

2. Patel p. 3 quotes McConaghy [15] as finding that "the incoming and outgoing v∞ at Earth (~4.71 km/s) are **identical**." This refers to the gravity-assist flyby constraint, not to a departure v∞. This 4.71 km/s figure is consistent with the "4.7/5.0" catalogue entry being the real-ephemeris McConaghy result, but Patel does not make this distinction.

3. The 3.657 km/s figure in Table 1 and the 4.71 km/s in the text are **inconsistent with each other** within Patel itself. The 3.657 km/s appears to correspond to the coplanar circular-orbit flyby v∞ for the first leg of the S1L1, while the 4.71 km/s is the symmetric Earth flyby v∞ for the gravity-assist maneuver. Patel does not reconcile these.

4. **The 5.65/3.05 pair (catalogue "spec §9") does not appear anywhere in Patel.** Its origin is not traceable to Patel or McConaghy 2006 from the text of this thesis.

5. **The 4.99/5.10 pair (Russell 2004 coplanar) does not appear in Patel.** Russell 2004 is cited as ref [12] but only for the general method of constructing cyclers, not for any specific trajectory numbers.

6. **The 4.7/5.0 pair** is most consistent with the McConaghy 2006 real-ephemeris results mentioned in the wider literature, but Patel 2019 never cites a real-ephemeris v∞ at all — it operates only in the coplanar model.

**Summary for target resolution:** Patel takes all S1L1 trajectory data from McConaghy 2006 in the **circular-coplanar** model. The Earth flyby v∞ used for vehicle design is **3.657 km/s** (Table 1, p. 5). No Mars V∞ is given. The topology is **E → M → E (short loop) → E (long loop)** with ~154-day E→M and M→E legs. The thesis does not contain real-ephemeris V∞ values and does not use or cite the 5.65/3.05, 4.99/5.10, or 4.7/5.0 pairs by those values.

---

## 2. EEM/EM Cyclers from Sanchez Net 2022

### 2.1 Definitions

Sanchez Net defines (Nomenclature, p. 861):
- **EM** = Earth–Mars sequence of visits (single-ellipse style, one Mars visit per cycle before returning to Earth)
- **EEM** = Earth–Earth–Mars (–Earth) sequence of visits (multi-arc, with an intermediate Earth flyby before the Mars visit)

The paper notes (p. 863, Sec. III.A): "A single cycle is either an Earth–Mars(–Earth) (EM) trajectory or an Earth–Earth–Mars(–Earth) (EEM) trajectory. There are many more combinations that could be explored."

### 2.2 Search constraints (p. 863)

The Star patched-conic enumeration uses:
- Launch date: 2030–2034
- Total time of flight: 15 years
- Number of cycles per trajectory: exactly 3
- Maximum Δv per flyby: **≤10 m/s** (near-ballistic filter; this is the unoptimized patched-conic value)
- **Batch 1**: every Mars–Earth (ME) leg shorter than **12 months** (to limit data-loss risk)
- **Batch 2**: ME legs shorter than **18 months**
- Over 7,900 trajectories returned before Set Cover selection

### 2.3 Example cycler orbits (Fig. 2, p. 862)

Figure 2 shows two annotated example cycler orbits with full event lists. These are the most complete citable cycler solutions in the paper.

**Cycler 1** (Fig. 2a caption — "1 Mars visit every 2 synodic periods"):

| Event # | Body | Date | v∞ (km/s) | Altitude (km) | Δv (m/s) |
|---|---|---|---|---|---|
| 1 | Earth | 05/11/2032 | 17.6 | — | — |
| 2 | Earth | 11/13/2032 | 4.283 | 6200 | — |
| 3 | Mars | 07/09/2035 | 6.207 | 65000 | — |
| 4 | Earth | 11/22/2035 | 3.605 | — | 5 m/s |
| 5 | Earth | 05/20/2036 | 3.709 | 14100 | — |
| 6 | Mars | 10/25/2041 | 5.210 | 2040 km | — |
| 7 | Earth | 07/07/2042 | 5.225 | 6350 | — |
| 8 | Earth | 07/18/2044 | 5.234 | 14000 | — |
| 9 | Mars | 02/06/2046 | 7.276 | 14800 | — |
| 10 | Earth | 06/13/2049 | 7.285 | — | — |

Notes: The caption reads "C∞ = 17.6 km²/s², Dec. = 43.7°" at first Earth (launch conditions). Δv = 5 m/s listed at event 4. This is an **EEM** structure: Earth (launch) → Earth flyby → Mars → Earth flyby → Earth flyby → Mars → … The intermediate Earth-to-Earth arc between events 1 and 2 is ~6 months (May to Nov 2032).

**Cycler 2** (Fig. 2b caption — "8 years between Mars 1 and 2"):

| Event # | Body | Date | v∞ (km/s) | Altitude (km) | Δv (m/s) |
|---|---|---|---|---|---|
| 1 | Earth | 01/02/2034 | 10.7 | — | — |
| 2 | Mars | 12/19/2037 | 6.466 | 10400 | — |
| 3 | Earth | 10/10/2040 | 3.090 | 370 | 7 m/s |
| 4 | Earth | — | 5.285 | 7560 | — |
| 5 | Mars | 11/02/2045 | 6.871 | — | 6 m/s |
| 6 | Earth | 08/24/2048 | 5.721 | — | Dec. = 20.8° |
| 7 | Earth | 09/19/2048 | — | — | — |

Notes: C∞ = 10.7 km²/s² at first Earth (launch). Δv = 7 m/s at event 3, 6 m/s at event 5. This is an **EM** structure at macro scale (Earth → Mars → Earth…) but with intermediate Earth flybys visible. ME transit for first Mars visit is very long (~33 months, 2037 to 2040).

**Important caveat:** Fig. 2 annotations are dense and some values are partially readable from the paper's small print. The values above are read directly from the figure captions as visible in the PDF. Treat individual flyby altitudes as approximate; exact values may require higher-resolution source.

### 2.4 Fleet / network results (Tables 3–5, p. 867)

The paper does not provide individual cycler orbital elements or V∞ tables for the full selected network. Instead it gives **downlink schedules** (date + cycler ID number + ME transit days) for representative networks:

**Table 3: Three-cycler network / nine downlinks** (p. 867)
Cycler IDs used: 51, 1, 1141. ME transit durations range 143–296 days.

| Downlink date | Cycler No. | ME transit (days) |
|---|---|---|
| 29 Oct. 2033 | 51 | 174.16 |
| 25 Nov. 2035 | 1 | 233.62 |
| 19 Dec. 2037 | 1141 | 165.09 |
| 24 April 2038 | 51 | 220.83 |
| 29 May 2040 | 1 | 332.28 |
| 30 May 2042 | 51 | 212.25 |
| 9 June 2044 | 1 | 278.25 |
| 24 Aug. 2046 | 1141 | 295.64 |
| 19 Sept. 2048 | 1141 | 143.25 |

**Table 4: Five-cycler network / 15 downlinks** (p. 867)
Cycler IDs: 50, 84, 397, 1066, 1141. ME transits range 105–365 days.

**Table 5: Six-cycler network / 18 downlinks** (p. 867)
Cycler IDs: 53, 84, 362, 529, 612, 1102. ME transits range 126–358 days.

The paper does not tabulate V∞ values, orbital elements, or leg ToFs for these numbered cyclers; they are identified only by index within the Star database. Their individual parameters are not recoverable from this paper alone.

### 2.5 ME transit duration statistics

From Fig. 3 (p. 863) and associated text:
- Full data set: ME transit durations roughly **11–15 years** total flight time per trajectory (3 cycles)
- ME leg individual durations: **batch 1** constrains each ME leg < 12 months; **batch 2** constrains < 18 months
- The selected network downlinks (Tables 3–5) show ME transit legs of **105–366 days** (approx. 3.5–12 months), consistent with batch 1 constraint

### 2.6 ΔV per flyby

All cyclers in the selected set have unoptimized Δv ≤ 10 m/s per flyby in the patched-conic model (p. 863, constraint 5). The paper notes the real optimized Δv will be higher once full-model perturbations are applied. Selected network cyclers need "at most 60 m/s once injected into their trajectories" (p. 863, second column). The specific per-flyby Δv values for the named network cyclers (IDs 51, 1, 1141, etc.) are not tabulated.

### 2.7 Fleet phasing / cadence

- 6 couriers can achieve ~yearly downlinks (f ≈ 1 per Earth-year) per Fig. 4 (p. 864)
- Phase offsets between couriers are implicit in the Set Cover solution but no explicit phase angle table is given
- Launch epoch spans 2030–2034 per enumeration constraint (p. 863)
- The 6-courier fleet spacing is set by the Set Cover optimizer, not a uniform phase formula; actual dates appear only in Table 5 (p. 867)

### 2.8 SEP insertion trajectory (Sec. VI, p. 868)

The paper separately analyzes a solar-electric propulsion (SEP) trajectory to insert a Data Mule Spacecraft (DMS) from Earth into the Earth–Mars cycler orbit (EMCO). Example: depart Earth 12 July 2035, arrive at EMCO 12 July 2037 (2 years), using NEXT ion engine (134 mN, Isp = 4310 s), wet mass 250 kg, fuel 41.1 kg. This is a one-off insertion example, not a citable cycler row.

### 2.9 Candidate catalogue rows from Sanchez Net

The only fully-specified cycler orbits with dateable events and partial V∞ data are **Cycler 1** and **Cycler 2** from Fig. 2 (p. 862). They are the best candidates for catalogue rows.

| Field | Cycler 1 (Fig. 2a) | Cycler 2 (Fig. 2b) |
|---|---|---|
| Sequence type | EEM (Earth–Earth–Mars per macro-arc) | EM (Earth–Mars per macro-arc) |
| Launch epoch | 05/11/2032 | 01/02/2034 |
| Launch C∞ | 17.6 km²/s² (v∞ = 4.20 km/s) | 10.7 km²/s² (v∞ = 3.27 km/s) |
| First Mars v∞ | 6.207 km/s (07/09/2035) | 6.466 km/s (12/19/2037) |
| First ME transit | event 3→4: ~136 days (Jul→Nov 2035) | event 2→3: ~1026 days (Dec 2037→Oct 2040) |
| Max Δv per flyby (unoptimized) | 5 m/s (event 4) | 7 m/s (event 3), 6 m/s (event 5) |
| Period (approx.) | ~4.3 yr between Mars visits 1 and 2 (2035→2041) | ~8 yr between Mars visits 1 and 2 ("8 years between Mars 1 and 2") |
| Source | Sanchez Net 2022, Fig. 2a, p. 862 | Sanchez Net 2022, Fig. 2b, p. 862 |
| Model | Patched conics (Star tool), unoptimized | Patched conics (Star tool), unoptimized |

**Data gaps for catalogue rows:**
- Orbital elements (a, e, i) not given for any of these cyclers
- V∞ at Earth for return legs (events 4, 5, 7 etc. in Cycler 1) are partially listed in Fig. 2 but not all are legible at paper resolution
- No epoch in J2000 or MJD form
- The numbered cyclers in Tables 3–5 (IDs 51, 84, 397, etc.) have no published parameters other than downlink date and ME transit duration; they cannot be used as catalogue rows without access to the Star database output

---

## 3. Recommendations

*(Clearly separated from sourced facts above — these are interpretive conclusions, not paper citations.)*

### 3a. V∞ target and leg topology for S1L1 differential corrector

Patel 2019 is a vehicle design thesis that takes its trajectory entirely from McConaghy 2006 in the **circular-coplanar** model. It provides no real-ephemeris V∞. The only V∞ values it gives are Earth-referenced (3.657 km/s flyby, 2.492 km/s launch departure), not Mars V∞.

For the differential corrector on real ephemeris:
- **Do not use Patel as the target.** Patel provides no Mars V∞ and operates in a different dynamical model.
- The **4.7/5.0 pair** (McConaghy 2006 real-ephemeris, catalogue entry) remains the best-sourced real-ephemeris target, since McConaghy is exactly Patel's ref [15] and the 4.7 figure at Earth is consistent with Patel's text reference to ~4.71 km/s.
- The **5.65/3.05 pair** has no traceable source in either paper mined here. It should remain flagged as unverified until the original §9 specification source is identified.
- The **4.99/5.10 pair** (Russell 2004) is a different coplanar solution and not what Patel or McConaghy describe; it is a separate catalogue entry.
- The topology to target is confirmed as **E → M on a ~154-day arc**, with the cycler encountering Earth first (at t=0), Mars ~154 days later, then returning to Earth on a short arc at t≈2.83 yr, then completing the long arc back to Earth at t≈4.29 yr. This is the **E-M-E-E** structure in the cycler literature.
- The closure failure (real-eph Mars V∞ ≥ 6.4 km/s) is consistent with the known difficulty of the real-ephemeris problem; the corrector may be landing on a different family. Resonance-anchored construction (as noted in project memory) is more likely to succeed than free optimisation from the coplanar seed.

### 3b. EEM/EM cyclers for catalogue

Cycler 1 and Cycler 2 from Sanchez Net Fig. 2 are citable but incomplete for full catalogue rows:
- **Cycler 1** (EEM, launch 2032) has short ME legs (~136 days) and low Δv, making it operationally attractive. Recommend adding as a candidate row with the data in §2.9 above, flagged as patched-conic only.
- **Cycler 2** (EM, launch 2034) has a 1026-day first ME transit, violating the <12-month bound Sanchez Net themselves impose. The second Mars visit has no legible ME transit in Fig. 2b. Do not add to catalogue without further verification.
- The numbered network cyclers (IDs 51, 84, etc.) cannot be added to the catalogue from this paper; the underlying Star database output is unpublished.
- Missing data for any catalogue row: a/e/i, real-ephemeris validation, Mars V∞, full per-leg ToF table, epoch in standard format.

---

## 4. Extraction Caveats

1. **Patel gives no Mars V∞.** The thesis design is sized on Earth flyby v∞ (3.657 km/s) and Earth-to-Mars transfer time (154 days). Mars V∞ is not tabulated or discussed as a design driver. The absence is real, not an extraction failure.

2. **Patel's Table 1 v∞ of 3.657 km/s vs. the 4.71 km/s in body text are both attributed to McConaghy [15]** but represent different quantities (possibly leg 1 vs. the symmetric flyby condition). Patel does not reconcile these. This is an ambiguity in the source document, not a reading error.

3. **Sanchez Net Fig. 2 annotations are at small scale.** Some flyby altitude and v∞ values for intermediate events in Cycler 1 and Cycler 2 were partially legible. Values marked with approximate signs above reflect genuine resolution limits in the printed figure. The event dates, primary body identifications, and primary v∞ values are clearly legible.

4. **Sanchez Net does not give orbital elements or individual leg times for the network cyclers** in Tables 3–5. These tables give only downlink date, cycler index, and ME transit duration. No V∞ or topology data can be extracted for these orbits from this paper.

5. **No additional trajectory-parameter pages in Patel were found beyond pages 3–8.** The remainder of the thesis (pp. 23–59) covers vehicle subsystem design (ECLSS, power, thermal, GNC, structures, propulsion) and does not revisit or refine the trajectory numbers. Pages 58–61 summarise mass/power totals only. The trajectory content is entirely in Ch. 1 (pp. 1–6) and Sec. 2.1 (pp. 7–9).

6. **No V∞ at Mars appears anywhere in Patel 2019.** This is the key finding: the thesis that is cited as the design source for the S1L1 vehicle does not specify Mars encounter V∞ at all.
