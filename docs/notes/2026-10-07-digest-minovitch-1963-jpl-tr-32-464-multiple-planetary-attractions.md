# Digest: Minovitch 1963, "The Determination and Characteristics of Ballistic Interplanetary Trajectories Under the Influence of Multiple Planetary Attractions" (JPL TR 32-464) (#960 batch 31)

M. A. Minovitch, Jet Propulsion Laboratory Technical Report No. 32-464, 31 October 1963. Prepared under NASA
contract NAS 7-100. 60 printed pages (viii + 60). No DOI and no NTRS record.
- **Source:** the HathiTrust public-domain scan of the Cornell Engineering Library copy (Google-digitised),
  hdl.handle.net/2027/coo.31924098468527. The file is `133c6079-JPL_TR_32-464_complete.pdf`, 72 scan pages, 18 MB,
  md5 `be5c46f703fd5f46b5080858fd8bca18`.
  - Library stamps on the copy: "Engineering Library DEC 17 1963" (front cover) and "Cornell University DEC 6 1963 Library" (back cover).
  - Scan page = printed page + 10, for printed pp. 1-60 (scan pp. 11-70).
  - Pages below are cited as "p. N (scan M)".
- **Proposed corpus filename:**
  `cyclers_pdf/papers/minovitch-1963-ballistic-interplanetary-trajectories-multiple-planetary-attractions-jpl-tr-32-464-hathitrust-coo-31924098468527.pdf`
- **Companion files (proposed), all built in this batch:**
  - `...-merged.txt`: the reconciled full text, page-marked, with every table cell.
  - `...-ocr-disagreement-report.txt`: every disagreement, the readings, the decided value and how it was decided.
  - `...-ocr-stats.md`: agreement statistics, including the Google layer.
- **How I read it:**
  - Three OCR witnesses, by the same method as the VanderVeen 1969 and Niehoff 1965 reports:
    - G: the Google/HathiTrust text layer;
    - T: tesseract at 400 dpi;
    - V: macOS Vision.
  - **Every table cell (4,236 cells in 24 tables) was read on the 400 dpi page image.** Tables 1-4 were read with no
    OCR draft in view.
  - Every display equation was read on the image, row by row, and so was every disputed or numeric token in the text.
  - The prose was proof-read page by page on column crops of the page images.
  - Arithmetic checks:
    - `arith.py`: 811 checks. They cover the flyby conic at every swingby (VACA and DA from HEV, DOCA and the
      Appendix B constants), the leg-time sums, the calendar sums, and the 1.84 launch-energy ratio on p. 17.
    - `digest_checks.py`: the trajectories in the text.
    - `vv_check.py`: the VanderVeen 1969 cross-check.
  - Every disagreement was re-read at full zoom. The ones that survive are listed in sec. 6 as source misprints.
- **Wanted list:** this is **row 1**. It can be closed when the PDF is filed.
  - The row notes that the report is "also cited as TR 32-468" (Niehoff 1965). That is wrong: the report number
    on the cover, the title page and every running head of this scan is 32-464.

## 0. Verdict

**What it is.** The 1963 JPL report of Minovitch's multi-planet patched-conic method, with the first published
numbers.
- It covers the whole 1965-1975 decade:
  - Earth-Venus-Mercury tours (one for each of the six Earth-Venus launch periods);
  - Earth-Venus-Mars tours;
  - manned Earth-Venus-Earth, Earth-Mars-Earth and **Earth-Venus-Mars-Earth** round trips;
  - and one 5-year **"space bus" trajectory, E-V-M-E-M-E-V-E**, offered as one member of an
    "interplanetary transportation network".
- It is the report that Menning, VanderVeen and Minovitch 1972 cite as TR 32-464. Niehoff 1965 cites it with the same
  title and date, but as TR 32-468.

**Periodic orbits and cyclers.**
- The report never uses the words "periodic" or "cycler", and it does not compute a repeating orbit.
- It does state two ideas that come near to one:
  - **A repeated E-V-M-E tour**, as a trajectory type (p. 15, scan 25; quoted in sec. 3).
  - **A transportation network of large space buses on long multi-planet free-fall trajectories that are "used as
    often as desired"**, with supplies renewed at each Earth rendezvous and crews carried by excursion modules
    (pp. 50-51, scan 60-61; quoted in sec. 3).
- Its one worked example (Table 23) is a **single, non-repeating** 1,811-day sequence. The report gives no closure
  condition, no period, and no statement that the bus returns to the same state.
- So it is a precursor in concept (a reusable interplanetary "space bus"), not a cycler. The prior-art verdict for
  #942 belongs to the other agent; this digest only lists the trajectories (sec. 2).

**What it gives the project:**
- (1) The 1963 numbers for every E-V-Me, E-V-M, E-V-E, E-M-E and E-V-M-E trajectory, image-read with the
  internal checks done (sec. 2). These are positive controls for patched-conic multi-flyby code. Each row has
  dates, HEV, flight times, B-plane components and the flyby DOCA, VACA and DA, which can be reproduced from
  Appendix B's constants.
- (2) The method's printed form (sec. 1), which matches the held TM 312-130 (1961).
- (3) The earliest dated "transportation network" wording (1963), five years before Minovitch's 1968 AAS paper
  "Gravity thrust and interplanetary transportation networks" (not held).

**Catalogue implication (PROPOSAL only).**
- No new cycler row: nothing here is periodic.
- If the project keeps a `precursor_mga` or `mga_tour` class for dated historical tours, these could become V0
  rows sourced to this report:
  - the Table 23 space-bus sequence;
  - the Aug 12, 1970 E-V-M-E manned flyby (Table 16);
  - the June 4, 1972 E-V-M-E trajectory "A" (text, p. 49).
  They would need `inserts_into` to be waived or left empty, since they insert into no cycler. Owner decision.
- A cheaper first step: file the report and point the existing wanted-list row 1 notes to it.

## 1. Method (Parts II-III, pp. 1-14, scan 11-24)

The model is a three-dimensional patched conic. The only assumption is that "at any instant, one, and only one,
gravitating body influences the vehicle's motion" (Abstract, p. viii, scan 10). Every equation below was read on
the image.

Part II, conic tools:
- The e and h vectors (Eqs. 1-7):
  - R × V = h (2);
  - V × h = μ(R̂ + e) (3), so e = (1/μ) V × h − R̂ (4);
  - R = l/(1 + e cos θ), with l = h²/μ (5-7).
- e and h from two position vectors (Eqs. 8-12):
  - h = ± (R₁ × R₂ / |R₁ × R₂|) (aμ|1 − e²|)^(1/2) (8);
  - e = αR₁ + βR₂ (9), with α and β from the 2 × 2 system (10-12).
- The velocity, V = (μ/h²) h × (R̂ + e) (13), and the energy equation V² = μ(2/R ∓ 1/a) (14).
- The classical elements from e and h, and back (15-23).
- **Lambert's theorem in Lagrange's form, (24)-(29).** These are six branches, for transfer angles in
  [0,180), [180,360) and [360,540) deg, each with F, F* not separated or separated by PQ. The 360-540 deg branch
  ("Type III") is one extra revolution. Table 15 uses Type III Mars-Earth legs. Nothing goes beyond 540 deg.
  - The report requires all transfer angles ≤ 540 deg (p. 10).
  - The eccentricity follows from (30)-(31).
- The launch-energy minimum is a Hohmann-like transfer, with θ = 0 and a = s/2 (p. 7).
- The asymptotes TA₁, TA₂ and the joining times TM₁, TM₂, TM₃ of the six T(a) branches are given (Fig. 5).
- Kepler's second law holds for every conic (Eq. 34). The time from perihelion is given as an integral (35).

Part III, the swingby:
- The sphere of influence is ρ* = (m/M)^(2/5) R (p. 9).
- **Fundamental equation (40):**
  V²(T₂*) − V²(T₁*) = 2 V₂ • [V(T₂*) − V(T₁*)].
  This is the heliocentric energy change at the flyby when |V′| in = |V′| out (38). That is the same
  V-infinity-magnitude closing condition as TM 312-130 Eq. (26).
- **Algorithm (p. 11-12):**
  - The first leg P₁P₂ comes from Lambert's theorem, with T₁ and T₂ given.
  - The second leg P₂P₃ comes from Lambert's theorem with a trial T₃.
  - T₃ is stepped until (40) holds. "In general, there is an infinite set of values of T₃ ... we shall choose that
    solution which gives T₃ − T₂ the smallest value" (p. 12), subject to d > 0 (Fig. 10: T₃′, T₃″, T₃‴, T₃⁗).
- **The flyby hyperbola:**
  - a₂ from (41);
  - e₂ = [2V′₁V′₂ / (V′₁V′₂ − V′₁ • V′₂)]^(1/2) (43);
  - d = a₂(e₂ − 1) − radius (44);
  - V′CA (45), h₂ (46), the vectors e₂ and h₂ (47-48);
  - the time in the sphere of influence, 2ΔT, from (49);
  - the entry and exit points (50).
- **Chaining:** "instead of terminating the mission at P₃ ... we simply take P₂P₃ as an initial condition for
  proceeding to P₄" (p. 14). So any P₁-P₂-...-Pₙ tour is a sequential date march with one free parameter per
  added flyby.
- **Planet ephemeris (Appendix A, pp. 58-59):**
  - an osculating ellipse is fitted through ephemeris positions from HMSO "Planetary Coordinates 1960-1980"
    (Ref. 7), spaced t₃ − t₀ = 10k days (30, 50, 70 and 90 days for Mercury to Mars);
  - interpolation uses a Taylor series in δ to the third term, or the fourth for Mercury (A-3).
- **Constants (Appendix B, p. 60):**
  - au = 1.495990 × 10⁸ km;
  - μSun = 2.9591221 × 10⁻⁴ au³/day²;
  - obliquity 23° 26′ 44.84″ (1950);
  - planet radii 2330.0, 6100.0, 6378.2 and 3415.0 km. Venus includes the cloud layer, so DOCA at Venus is measured to
    the cloud top;
  - μ = 4.835167 × 10⁻¹¹, 7.241303 × 10⁻¹⁰, 8.887552 × 10⁻¹⁰ and 9.582649 × 10⁻¹¹ au³/day²;
  - the (m/M)^(2/5) column is consistent with these μ to 0.03-0.16%.
- **Accuracy (p. 52):** JPL's precision integrator, started from these patched-conic solutions, showed "very rapid
  convergence".

**Comparison with the held TM 312-130 (23 Aug 1961, digest 2026-10-05-966):**
- The algorithm is the same:
  - the same model assumptions;
  - the same ρ*;
  - the same energy-match unknown (TM Eq. 26 = TR Eq. 40);
  - the same choice of the next encounter date as the unknown;
  - the same d > 0 gate;
  - the same hyperbola formulas (TM 30-31 = TR 41, 43; TM's d = TR 44).
- What the TR adds:
  - the e- and h-vector machinery for three-dimensional conics (Part II);
  - Lambert in six branches up to 540 deg;
  - the time-in-sphere formula (49);
  - the ephemeris method (Appendix A);
  - and, above all, **numbers**: TM 312-130 had none.
- What the TR limits: TM 312-130 allowed any number k of complete solar circuits on a leg (adding kP). The TR caps
  every leg at 540 deg, that is, k ≤ 1 (Type III).
- The 1961 example E-V-M-E-S-P-J-E is not computed here. The computed tours stop at Mars.
- The Abstract says the conic study was done "at the Jet Propulsion Laboratory during the summer of 1961" (p. viii).
  The E-V-M-E discovery is dated "the early part of June 1962" (p. 40).

## 2. Trajectories with numbers

The values are as printed and image-read. The full tables are in merged.txt, and as one file per table in
`tables/tNN.txt`. Notation:
- HEVₖ is the hyperbolic excess speed (km/s);
- Tₖ,ₖ₊₁ is the leg time (d);
- θ is the heliocentric transfer angle (deg);
- DOCA is the altitude at closest approach (km, measured to Venus's cloud top);
- VACA is the speed at closest approach;
- DA is the turn angle;
- TFT is the total flight time.
Dates are at 1200 GMT (p. 21).

### 2.1 Direct transfers (Tables 1-3, pp. 17-19; Tables 18-19, p. 47)
- **Table 1**, Earth-Venus: 6 Type I and 6 Type II launch periods, 1965-1973.
  - The midpoints are about 19.2 months apart. Arithmetic: 19.0-19.6 months for Type I.
  - Check: (3.65/2.69)² = 1.84, matching the text on p. 17.
- **Table 2**, Earth-Mercury: 11 Type I and 11 Type II periods, 1965-1974. Minimum HEV₁ is 6.43 km/s.
- **Table 3**, Earth-Mars: 4 + 4 periods.
  - The minimum is HEV₁ = 2.81 km/s (Type I, 1971).
  - The periods are about 780 days apart (p. 19).
- **Tables 18-19**: Venus-Earth (1970-1974) and Mars-Earth (1970-1973) optimum return periods.

### 2.2 Earth-Venus-Mercury (Tables 4-10)
**All six Earth-Venus launch periods of the decade also allow E-V-Me** (p. 20). These rows match VanderVeen 1969
Table 2 (held) to ±1 day and ±1 deg on all six launch rows (`vv_check_out.txt`). The 1969 total, which VanderVeen
gives as 297 d, exposes the Table 7 misprint (sec. 6).

| Table (page) | Launch dates | Rows | Representative row (as printed) |
|---|---|---|---|
| 4 (p. 22) | Nov 28, 1965 - Jan 5, 1966 | 20 | Dec 18, 1965: HEV₁ 3.97, T₁₂ 170.17, DOCA 1,560, DA 56.57, T₂₃ 105.62, HEV₃ 9.83, TFT 275.43 |
| 5 (p. 24) | June 7, 1967, Venus date stepped by 0.01 d | 26 | DOCA rises from 6.7 to 200.9 km, then falls to 31.0 km; T₁₂ 107.32 to 107.58 |
| 6 (p. 25) | June 5 - July 3, 1967 | 15 | June 19: HEV₁ 3.70, T₁₂ 96.28, DOCA 311, DA 65.49, T₂₃ 71.49, TFT 167.77 |
| 7 (p. 27) | Jan 3 - Feb 2, 1969 | 16 | Jan 23: HEV₁ 3.96, T₁₂ 189.38, DOCA 645, T₂₃ 107.67, TFT 279.05 (misprint for 297.05) |
| 8 (p. 29) | July 25 - Sept 13, 1970 | 26 | Aug 18: HEV₁ 3.50, T₁₂ 101.22, DOCA 3,768, T₂₃ 59.00, HEV₃ 11.91, TFT 160.22 |
| 9 (p. 31) | Mar 18 - Apr 21, 1972 | 18 | Apr 1: HEV₁ 4.03, T₁₂ 196.58, DOCA 521, T₂₃ 85.00, TFT 281.58 |
| 10 (p. 33) | Oct 21 - Nov 16, 1973 | 14 | Nov 4: HEV₁ 4.25, T₁₂ 93.79, DOCA 2,983, T₂₃ 58.00, TFT 151.79 |

### 2.3 Earth-Venus-Mars (Tables 11-13)

| Table (page) | Launch dates | Rows | Notes |
|---|---|---|---|
| 11 (p. 35) | Dec 22, 1968 - Jan 21, 1969 | 16 | TFT about 500 d; Venus-Mars leg Type II (T₂₃ about 390-405 d) |
| 12 (p. 37) | July 23 - Aug 28, 1970 | 19 | Aug 12: HEV₁ 3.26, T₁₂ 129.28, DOCA 3,850, T₂₃ 180.00, HEV₃ 6.75, TFT 309.28 |
| 13 (p. 38) | May 11 - June 6, 1972 | 14 | the table has no HEV₂ column; TFT 283-348 d |

### 2.4 Round trips (manned reconnaissance), sec. IV.B
- **Earth-Venus-Earth (Table 14, p. 41).** Three free returns in about one year:

  | Launch | HEV₁ | T₁₂ | DOCA | T₂₃ | HEV₃ | TFT | Earth return |
  |---|---|---|---|---|---|---|---|
  | 8/20/70 | 2.92 | 114.00 | 725.5 | 250.96 | 7.13 | 364.96 | Aug 19, 1971 |
  | 4/3/72 | 3.69 | 114.00 | 3,983 | 260.95 | 8.18 | 374.95 | Apr 12, 1973 |
  | 11/4/73 | 3.76 | 110.00 | 5,343.5 | 269.46 | 7.90 | 385.46 (sum 379.46) | Nov 24, 1974 |

  The 1970 case returns to Earth 364.96 d after launch, which is 0.999 yr. This is a one-year Earth-Venus-Earth
  free return (`digest_checks_out.txt`). The report does not comment on the one-year period or suggest repeating it.
- **Earth-Mars-Earth (Table 15, p. 41).**
  - 6/8/71: TFT 1,111.83 d.
  - 8/20/73: TFT 1,027.7 d.
  Both have Type III Mars-Earth legs. The 1973 case spends 5.46 d inside Mars's sphere, at HEV₂ 2.53 km/s (p. 40).
- **Earth-Venus-Mars-Earth.**
  - **1970, Table 16 (p. 42):** launches July 15 - Aug 28, 1970, TFT 607.82-650.65 d.
    - Aug 12, 1970: HEV₁ 3.26; Venus DOCA 3,848 km after 129.28 d; Mars DOCA 6,590 km; T₃₄ 312.36 d;
      Earth-return HEV₄ 9.34; TFT 621.63 d.
    - Text (p. 40): this needs "only about one-half of the launch energies and flight times required by the best
      Earth-Mars-Earth reconnaissance trajectories".
  - **1972, Table 17 (p. 44):** launches May 13 - June 8, 1972, TFT 459.07-481.47 d. May 27 is the figured case.
  - **Trajectory A (text, p. 49; not in the tables):** E-V-M-E launched June 4, 1972.
    - Encounters Nov 19, 1972, May 23, 1973 and Oct 17, 1973.
    - HEV₁ 4.33, DOCA₂ 9,164, HEV₃ 5.97, DOCA₃ 609, HEV₄ 9.51.
    - T₁₂ 167.56, T₂₃ 185.44, T₃₄ 146.67, TFT 499.67 (the dates and the sum check).
  - **Trajectory B (text, p. 49):** E-V-M launched May 31, 1972.
    - Encounters Nov 17, 1972 and May 12, 1973.
    - HEV₁ 4.27, DOCA₂ 9,223, HEV₃ 6.03, T₁₂ 170.00, T₂₃ 175.65, TFT 345.64.
    - B lands a crew on Mars on May 12, 1973. A passes Mars 11 days later, on May 23, and collects the crew, who fly up
      in the Mars excursion module ("the Saturn 5 possibility", pp. 45-50).

### 2.5 Manned landing profiles (Tables 20-22, p. 48)
These are E-V(stay)-E, E-M(stay)-E, and E-V-M(stay)-E with conjunction-style stays (ΔT, days on the planet).
- Table 20, Venus stays:
  - Aug 18, 1970: 450 d total, stay 32 d;
  - Mar 26, 1972: 416 d total, stay 20 d;
  - Nov 10, 1973: 404 d total, stay 14 d.
- Table 21, Mars stays:
  - May 19, 1971: 408 d, stay 9 d;
  - July 27, 1973: 442 d, stay 9 d.
- Table 22, via Venus:
  - Aug 12, 1970: 574 d, stay 19 d;
  - June 4, 1972: 586 d, stay 61 d. The dates and the durations of this row disagree by 16 d (sec. 6).
- All calendar sums check, except in Table 22 row 2.
- The Table 22 rows are the E-V-M legs of Table 12 (Aug 12, 1970: HEV₁ 3.26, TFT 309.28) and Table 13 (June 4, 1972:
  HEV₁ 4.32, TFT 335.13, HEV₃ 6.19).

### 2.6 The space-bus network trajectory (Table 23, p. 51)
E-V-M-E-M-E-V-E:

| i | Planet | Date | HEV | Leg into i (d) | θ (deg) | DOCA (km) | VACA | TISI | DA |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Earth | 8/14/70 | 3.28 | — | — | — | — | — | — |
| 2 | Venus | 12/20/70 | 5.44 | 128.20 | 151.25 | 3,817 | 9.75 | 2.47 | 63.44 |
| 3 | Mars | 6/17/71 | 6.74 | 179.30 | 171.65 | 6,838 | 7.34 | 1.85 | 9.69 |
| 4 | Earth | 4/25/72 | 9.34 | 312.22 | 290.81 | 8,089 | 11.93 | 2.26 | 27.79 |
| 5 | Mars | 12/11/73 | 13.46 | 595.03 | 203.10 | 592 | 14.23 | 0.97 | 6.41 |
| 6 | Earth | 5/31/74 | 11.81 | 171.58 | 191.74 | 6,249 | 14.23 | 1.82 | 21.27 |
| 7 | Venus | 3/26/75 | 10.83 | 298.37 | 196.13 | 9,455 | 12.61 | 1.29 | 17.39 |
| 8 | Earth | 7/30/75 | 7.91 | 126.26 | 220.57 | — | — | — | — |

- The total is 1,811 d (4.96 yr). The Earth-to-Earth intervals are 620, 766 and 425 d.
- **The first three legs are almost the Table 16 Aug 14, 1970 E-V-M-E row.** Compare (Table 23 vs Table 16):
  - T₁₂ 128.20 vs 128.16;
  - Venus HEV 5.44 vs 5.43;
  - T₂₃ 179.30 vs 180.00;
  - Mars HEV 6.74 vs 6.70;
  - T₃₄ 312.22 vs 311.66;
  - Earth HEV 9.34 vs 9.30.
  So the space bus is the 1970 reconnaissance trajectory, extended by a second Mars, Earth and Venus sequence.
- The arrival HEV changes from encounter to encounter, which is consistent with a non-periodic tour.
- The VACA and DA at all six flybys follow from HEV and DOCA by the conic formulas. The worst case is 0.06 deg at
  Venus i = 2.
- Two notes on the printed table:
  - The leg-time row is labelled "Tᵢ, ᵢ₊₁" but holds the leg ending at column i.
  - The DOCA unit is printed "km/sec".

## 3. Statements on repeating or periodic trajectories (quoted; image-checked)

- **p. 15 (scan 25), example 5** of the trajectory types considered:
  > "Trajectories of a vehicle that is launched from Earth at time T₁ and that at time T₂ makes a closest approach to
  > Venus which causes the vehicle to intercept Mars; the gravitational influence of Mars causes the vehicle to return
  > to Earth where the Earth's gravitational influence causes the vehicle to repeat the same flight; the vehicle is
  > sent to Venus such that the Venusian influence sends it to Mars whereupon the Martian gravitational influence
  > causes the vehicle to return to Earth (Earth-Venus-Mars-Earth-Venus-Mars-Earth)."

  It is listed as a type. No example of type 5 is computed anywhere in the report.
- **p. 15:** the number of sequences P₁-P₂-P₃ is "9³⁻¹ or 81!", and "in general" the number with n − 1 encounters is
  printed "9ⁿ". The text's own examples imply 9ⁿ⁻¹ (example 4, E-V-M-E, "is only one of the 729 different types").
- **pp. 50-51 (scan 60-61), sec. IV.C, "Interplanetary Transportation Networks to Support Manned Bases on Venus and
  Mars":**
  > "This problem of economics can be conveniently solved by constructing a long-lasting interplanetary
  > transportation network designed for the sole purpose of transporting personnel from one planet to another."

  > "Preliminary calculations have shown that if the planets Pᵢ are restricted to Mercury, Venus, Earth, and Mars,
  > where P₁ = Earth and Pᵢ ≠ Pᵢ₊₁ for i = 1, 2, • • • , n, it is possible to find sequences P₁-P₂-• • •-Pₙ such that the
  > flight times Tᵢ₊₁ - Tᵢ are comparable to those required for optimum Pᵢ-Pᵢ₊₁ transfers. Moreover, many of these
  > multiplanet trajectories were found to have very low launch energies."

  > "Each vehicle will carry extra provisions and life-support equipment to last until it makes its first Earth
  > rendezvous, whereupon its supply can be replenished to last until it makes its second Earth rendezvous, etc."

  > "One notices that all the vehicles involved in the network (which could be made up of several dozen space buses
  > all on their own different interplanetary trajectories) can be used as often as desired." (p. 51)

  The space buses are toroidal, "20 to 60 persons", with an outside diameter of "perhaps 200 to 300 ft", for
  artificial gravity (p. 50). Excursion modules and tankers carry the crews between bus and planet.
- **There is no statement** of periodicity, of a repeating geometry, or of a bus returning to its initial state. "As
  often as desired" refers to the long multi-planet path with repeated Earth passes, and the worked example
  (Table 23) is one finite sequence.

## 4. Other content worth a line

- Mercury missions: the E-V-Me tours need "less than one-third of the launch energies" of direct Earth-Mercury
  trajectories in 1965-1966 (p. 21), and arrive with "less than one-half" the approach energy (p. 21).
- 1969 Mars: E-V-M arrives with about nine times the approach energy of direct flights (p. 34). "The direct-flight
  trajectories of 1969 should be employed."
- Saturn 5 feasibility arithmetic for the Mars landing (pp. 47-50):
  - re-entry speed √(9.34² + 11.00²) = 14.43 km/s;
  - a 62,000 lb fly-by vehicle against 78,300 lb of Saturn 5 capacity on that trajectory;
  - a 145,000 lb B vehicle, so two Saturn 5s.
  - The rocket-equation results agree with the printed values to 1-2%. These are rough calculations, rounded in the
    source.

## 5. Citation mining (references of TR 32-464, p. 57)

| # | Reference | Held? | Wanted list |
|---|---|---|---|
| 1 | Wintner, A., Analytical Foundations of Celestial Mechanics, Princeton Univ. Press, 1947 | not held | not listed (textbook) |
| 2 | Battin, R. H., "The Determination of Round-Trip Planetary Reconnaissance Trajectories," J. Aero/Space Sci. 26(9), Sept 1959 | not held | row 42 (Battin 1959) |
| 3 | Clarke, V. C., Jr., A Summary of the Characteristics of Ballistic Interplanetary Trajectories, 1962-1977, JPL TR 32-209, 15 Jan 1962 | not held | not listed. It is the source of the Earth-Venus launch periods (p. 16). Low priority |
| 4 | Makemson, M. W., Baker, R. M. L., Jr., and Westrom, G. B., "Analysis and Standardization of Astrodynamic Constants," J. Astronaut. Sci. 8(1), Spring 1961 | not held | not listed (constants) |
| 5 | Hammock, D., and Jackson, B., "Vehicle Design for Mars Landing and Return to Mars Orbit," AAS Symposium on the Exploration of Mars, Denver, 6-7 June 1963 | not held | not listed |
| 6 | Dixon, F., and Stimpson, L., "A Systems Approach to Vehicle Design for Earth Re-entry From an Interplanetary Mission," AAS Symposium on the Exploration of Mars, Denver, 6-7 June 1963 | not held | not listed |
| 7 | Planetary Coordinates for the Years 1960-1980, HMSO (H. M. Almanac Office), London, 1958 | not held | not listed (ephemeris) |
| 8 | Clarke, V. C., Jr., Constants and Related Data Used in Trajectory Calculations at the JPL, JPL TR 32-273, 1 May 1962 | not held | not listed. The Niehoff 1965 digest also notes it as not held |

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i CORPUS_INDEX.md` for each author: Wintner, Battin,
Clarke, Makemson, Hammock, Dixon, Stimpson, and "planetary coordinates". None is held.

The report cites none of the cycler literature, which comes later. It does not cite TM 312-130 either.

## 6. Suspected misprints (printed values kept; each re-read at full zoom on the image)

Table cells:
- **Table 4, Dec 18, 1965:** TFT 275.43, but T₁₂ + T₂₃ = 275.79. VanderVeen 1969 gives 276.
- **Table 4, Dec 4:** TFT 290.91, but the sum is 290.97.
- **Table 5, last row:** VACA 12.17. The conic from DOCA 31.0 and HEV₂ 6.39 gives 12.115, and the printed DA (68.74)
  matches the conic. So VACA is probably 12.12.
- **Table 7, Jan 23, 1969:** TFT 279.05, but the sum is 297.05 (transposition). VanderVeen gives 297.
- **Table 11, Jan 19, 1969:** TFT 492.42, but the sum is 497.42.
- **Table 13, May 17, 1972:** T₂₃ 155.50 should be 115.50. TFT − T₁₂ = 115.50, and the T₂₃ trend agrees.
- **Table 14, Aug 20, 1970:** DA 76.16. The conic from VACA 11.39 and DOCA 725.5 gives about 70.8. **Unresolved.**
- **Table 14, Nov 4, 1973:** TFT 385.46, but the sum is 379.46.
- **Table 15, Aug 20, 1973:** TFT 1,027.7, but the sum is 1,026.72.
- **Table 16, July 29, 1970:** T₂₃ 210.85. From TFT, T₂₃ = 201.85 (transposition). Table 12's July 29 T₂₃ (201.95) supports this.
- **Table 16, Aug 2, 1970:** TFT 630.77, but the leg sum is 627.77. T₃₄ = 315.08 would fit. Table 12's Aug 2 TFT (315.69)
  equals T₁₂ + T₂₃, so T₃₄ is the misprinted cell.
- **Table 22, row 2 (June 4, 1972):** the launch date is right, because it is the Table 13 June 4 row. But June 4 + 335 d is
  May 5, 1973, only 45 d before the printed Mars departure of June 19 (ΔT is printed 61). The printed total, 586 =
  335 + 61 + 190, disagrees with the printed dates (570 d). Either ΔT and the total are misprinted, or the two later dates
  are 16 d early. **Unresolved.**

Text:
- **Nomenclature (p. 56):** "β ... defined by Eq. (11)". β is defined in Eq. (10).
- **p. 15:** the count of sequences is printed "9ⁿ". The worked cases (9³⁻¹ = 81; 729 for P₁-P₂-P₃-P₄) imply 9ⁿ⁻¹.
- **Report number:** Niehoff 1965 cites this report as "TR 32-468". The scan reads 32-464 throughout. Minovitch 1972,
  VanderVeen 1969 and Golubev 2014 all use 32-464.

## 7. Files (in `dg31-minovitch-ocr/`)

- `merged.txt`, `ocr-disagreement-report.txt` and `ocr-stats.md`: the deliverables.
- Witnesses:
  - G per page: `G/pNN.txt` (layout) and `Gb/pNN.html` (bbox);
  - T per page: `T/pNN.txt` and `.tsv`;
  - V per page: `V/pNN.txt` and `.tsv`;
  - `vocr.swift` and the binary `vocr`.
- Tables:
  - `tables/tNN.txt`: the image transcriptions;
  - `tables/score_tNN.tsv`: the per-cell witness scores;
  - `tables/score_summary.txt`.
- Body text:
  - `body/pNN.tsv`: the aligned slots;
  - `dec/pNN.tsv`: the image decisions;
  - `merged_body.json`, `report_body.txt`, `stats.json`.
- Checks:
  - `arith.py` / `arith_out.txt`;
  - `digest_checks.py` / `digest_checks_out.txt`;
  - `vv_check.py` / `vv_check_out.txt`;
  - `blind.py` / `blind_read.txt` / `blind_key.tsv`;
  - `gsample.py` / `gsample_read.txt` / `gsample_key.json` / `gsample_score.py` / `gsample_out.txt`.
- Scripts: `wit.py`, `align_body.py`, `finalize.py`, `score_tables.py`, `build_outputs.py`, plus the helpers
  `lines.py`, `need.py`, `cols.py`, `crop.py`, `sheet.py`, `stack.py`, `draft.py`, `mk.py`, `fin.py`, `tocdraft.py`
  and `zoom.py`.
- The page renders (`png/`, `view/`) are deleted.

## 7. #942 R1 prior-art verdict (reported to the lead and twobody-gen2-opus first, 2026-10-07)
- **ev-A (E-V, k=2, 4.89/10.36 km/s): no collision.** The only Earth-Venus returns are the three one-shot
  ~1-yr E-V-E free returns of Table 14 (p.41), HEV1 2.92-3.76 out, 7.13-8.18 back. Asymmetric, not repeated,
  no k=2 structure. Not even related prior art beyond the generic E-V-E free-return idea.
- **ev-C (E-V, k=2, 9.07/13.17, Earth-massless, Venus-hosted): no collision.** No trajectory returns to
  Venus, no Venus-Venus leg. Table 14 as above.
- **vm2-1 (V-M, k=3, 1001.770 d): no collision; weak related prior art.** Venus-Mars passages appear only
  inside one-shot E-V-M-E chains (Tables 16-17, Table 23 legs 2->3), with Mars bending (DA 9.69 deg in
  Table 23; Mars DA 7.68-16.20 deg in Table 16). There is no Mars-to-Venus leg and no Venus-to-Venus return.
  Venus HEV in Table 16 (5.37-6.11 km/s) brackets vm2-1's 6.086 km/s, but that is a one-way E-V-M-E.
- **Concept priority (for the record, not a collision):** p.15 example 5 is a 1963 statement of a
  repeating E-V-M-E-V-M-E free-fall itinerary, and p.50-51 is a reusable "space bus" network. These predate
  Hollister 1969 / Rall 1969 as ideas but give no periodic solution. This fits the #966 gallery digest
  (TM 312-130, 1961: the method; TR 32-464, 1963: first numbers) and is consistent with Menning's citation of
  TR 32-464 for E-V-M round trips. It does not touch the Hollister & Menning 1970 Table 3 rows.

Checked on the images by a separate fast pass and by the filer (p.15 example 5 and Table 23, p.51).

*Filed as `cyclers_pdf/papers/minovitch-1963-ballistic-interplanetary-trajectories-multiple-planetary-attractions-jpl-tr-32-464-hathitrust-coo-31924098468527.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
