# Digest: Niehoff 1965, "An Analysis of Gravity Assisted Trajectories in the Ecliptic Plane", IITRI Report T-12 (#960 batch 24; R1 history)

J. Niehoff, IIT Research Institute, Astro Sciences Center, Report No. T-12, 25 May 1965 (approved by
C. A. Stone). Contract NASr-65(06) for NASA Lunar and Planetary Programs. N65-34460, NTRS 19650024859.
- Filed as `cyclers_pdf/papers/niehoff-1965-analysis-gravity-assisted-trajectories-ecliptic-plane-iitri-report-t-12-n65-34460-ntrs-19650024859.pdf`.
  - This is an `ocrmypdf --redo-ocr` copy: the original page images are kept and the text layer is
    replaced with tesseract.
  - The source scan is md5 efad7a65e31b78d332d9513141f0f454, 79 pp.: front matter, 27 pp. of text,
    Appendices A-C, and 24 figures.
  - Companions: `.txt` (three-witness reconciled text) and `-ocr-disagreement-report.txt`.
- I read secs. 1 and 5-6, the Summary, the references and Appendix A from the reconciled text. Tables
  1-3 and all 24 figures were read on the 400 dpi page images (see sec. 4).

## 0. Verdict

**No repeating, periodic, free-return or cycler trajectory anywhere in the report.**
- It is a one-way, single-assist (Earth - P2 - P3) parametric study in a planar circular solar system
  (Mercury eccentric), with patched conics.
- Earth-Venus, Earth-Mars and Venus-Mars round trips are not discussed. Ross 1963 and Hollister are not
  cited.
- **R1 history:** this is the 1965 root of Niehoff's mission-analysis line, not of his cycler work.
  - The later cycler papers (AAS 86-172 "Pathways to Mars", VISIT; wanted rank 11) have no antecedent
    here.
  - The only periodicity remark is about launch opportunities: Jupiter-assisted solar probes recur
    "about every 13 months" (Jupiter's synodic period) with "only minor changes in trajectory
    requirements" (p.25). The Earth-Mars synodic period (780 d) makes Mars-assisted Jupiter missions
    rare (p.20).
- **Minovitch citation (for the lead):** the reference list prints "Minovitch, Michael A., 'The
  Determination and Characteristics of Ballistic Interplanetary Trajectories Under the Influence of
  Multiple Planetary Attractions', **JPL TR No. 32-468**, October 31, 1963" (read on the page image).
  - VanderVeen 1969 and Golubev 2014 give the same title and date as **TR 32-464**.
  - No NTRS id or N-number is given.
  - **If an NTRS search for 32-464 fails, try 32-468.** One of the two is a typo; two sources against
    one favour 32-464.
- Niehoff used Minovitch's 1970 E-V-Me launch window (7/25/70 to 9/13/70) as a check. His 2-D minimum-
  time curve is within 1 % of ideal velocity of Minovitch's 3-D window at the closest point (Fig. 7).

## 1. Content (READ)

- **Method (secs. 2-4, Appendices B-C):** a 2-D gravity-assist model with circles of influence
  (Table 2 radii).
  - The analytic maximum velocity change is dV_max = sqrt(K/R0), the circular speed at the surface.
    The maximum energy change follows from it (App. C).
  - The numerical approach sweeps ideal velocity (Earth launch dV including a 36,178 ft/s escape
    characteristic velocity and 4000 ft/s losses, App. A), miss distance and the injection angle.
  - Results are presented as "data" graphs (trip time vs miss distance at fixed ideal velocity) and
    "minimum time" graphs.
- **Sec. 5.1, Earth-Venus-Mercury** (Mercury at perihelion 0.31 AU and aphelion 0.47 AU):
  - Venus fly-by beats direct below 52,200 ft/s (0.31 AU). At 0.47 AU the break-even is 46,250 ft/s and
    95 d.
  - Examples: direct to 0.31 AU needs 51,400 ft/s and 96 d; via Venus the same trip time needs 48,700
    ft/s; a 120 d trip needs 45,000 ft/s.
  - 0.47 AU via Venus in 170 d needs 41,500 ft/s. That gives Atlas-Agena 580 lb or Atlas-Centaur 1900 lb
    (Summary: injected weight 400 to 1200 lb for the same 115 d trip).
  - Fig. 8 shows the August 1970 mission (sec. 4 below). Guidance: 150 lb of propellant on a 1300 lb
    spacecraft (Cutting & Sturms 1964).
- **Sec. 5.2, Earth-Venus-Jupiter:** the least favourable combination. Post-Venus aphelia stay below
  3 AU for 48,000-54,000 ft/s, against aphelion above 10 AU direct at 54,000 ft/s.
- **Sec. 5.3, Earth-Mars-Jupiter:** little gain except below 51,000 ft/s. The next opportunity is 1984
  (Fig. 11). The Mars approach speed is 14.9 km/s, more than twice a direct Mars fly-by.
- **Sec. 5.4, Earth-Jupiter-Saturn:**
  - The minimum-time miss distance is several Jupiter radii (geometric trade-off, Fig. 14).
  - A minimum-energy trajectory to Jupiter is extended to about 370 AU aphelion by the assist.
  - At 3.5 yr the ideal velocity falls from 55,000 (direct) to 52,000 ft/s.
  - 1977 launch (also 1976 and 1978). Equivalent dV at Jupiter is 18.7 km/s, against a 21 km/s maximum
    at 4 R_J.
- **Sec. 5.5, outer regions:** 20 AU at 56,000 ft/s takes 5 yr instead of 13 yr. 50 AU takes 11.5 yr,
  and a direct flight cannot reach 50 AU at that speed. Thrusted plus gravity-assist hybrids are
  recommended.
- **Sec. 5.6, solar probe via Jupiter:** about 55,000 ft/s for 0.1 AU or Sun impact, against about
  70,000 ft/s direct to 0.1 AU and almost 100,000 ft/s to 0.005 AU. Opportunities recur about every
  13 months.

## 2. Tables (READ on the page images; arithmetic checks in the disagreement report)

- **Table 1, planet data** (Clarke 1962 K; Explanatory Supplement 1961; Ehricke 1960). Transcribed in
  full in the `.txt`.
  - Every perihelion and aphelion equals a(1 -/+ e).
  - Mean orbital velocity equals a^-1/2, except **Mars: printed 0.8068546 against 0.8101241**, a source
    discrepancy.
- **Table 2, sphere-of-influence radii:** e.g. Jupiter 0.3216 AU, Venus 0.00412 AU, Earth 0.00618 AU,
  Mars 0.00378 AU. All km values equal AU x 1.49599e8.
- **Table 3, ordered maxima** (dV_max km/s; dE_max km^2/s^2):
  - Jupiter 42.6 / 583.7; Saturn 25.7 / 261.7; Neptune 16.8 / 91.9; Uranus 15.1 / 107.6;
    Earth 7.9 / 239.4; Venus 7.2 / 255.2; Pluto 6.9 / 42.1; Mars 3.6 / 95.5; Mercury 3.0 / 173.7.
  - 3B also gives the heliocentric peri/aphelion before and after the maximum-energy assist.
  - **The Summary says 42.5 km/s and 584 km^2/s^2** (rounding of 42.58 and 583.7).

## 3. Source inconsistencies (all confirmed on the page images)

- (a) Text p.21: "minimum energy trajectory to Jupiter (dV = 51,200 ft/sec)". Figs. 13, 15, 17, 18 and
  20 all say 50,200 ft/sec.
- (b) Summary 42.5 / 584 against Table 3 42.6 / 583.7.
- (c) Table 1 Mars mean orbital velocity.
- (d) Minovitch TR 32-468 against 32-464 elsewhere.

## 4. Mission illustrations (figure data boxes, read on the page images; tesseract and Vision agree)

| Fig. | Mission | Launch | Ideal dV ft/s | Key data |
|---|---|---|---|---|
| 8 | E-V-Me | Aug 1970 | 41,950 | 100 d to Venus, V approach 7.1 km/s, closest 1.66 R_V, total 160 d, Me approach 9.2 km/s, max comm 1.15 AU, equiv. dV at Venus 5.48 km/s |
| 11 | E-Ma-J | 2 Mar 1984 | 50,000 | Mars approach 14.9 km/s, closest 2 R_M, J approach 7.6 km/s, 706 d, max comm 6.2 AU, equiv. dV at Mars 1.6 km/s |
| 16 | E-J-S | Sep 1977 (also Jul-Aug 1976) | 54,000 | 502 d to Jupiter, J approach 12.7 km/s, closest 4 R_J, total 1072 d, S approach 17.8 km/s, max comm 10.5 AU, equiv. dV 18.7 km/s |
| 22 | E-J-solar probe | Jan 1970 to Aug 31 1977, eight dates about 13 months apart | 54,000 | J approach 12.7 km/s, closest 5.3 R_J, perihelion 0.02 AU, 298 km/s there, 3 yr, max comm 5.5 AU, equiv. dV 17.5 km/s |

All 24 figures' curve labels and axes are listed in the `.txt`.

## 5. Citation mining (pp.28-29)

| Item | Status |
|---|---|
| Clarke 1962, JPL TR 32-273 (constants) | Not held; constants only, not listed |
| Cutting & Sturms 1964, JPL TM 312-505 | Not held; already in wanted rank 71 (from VanderVeen) |
| Dobson, Huff & Zimmerman 1962, NASA TN D-1106 | Not held; osculating elements, not listed |
| Ehricke 1960, Space Flight I | Not held; textbook, not listed |
| Hohmann 1925 | Not held; history, not listed |
| Hunter 1964, "Future Unmanned Exploration of the Solar System", Astronautics and Aeronautics, May 1964, p.16 | Not held. Jupiter-assisted solar probe priority. Added to the wanted list |
| **Minovitch 1963, "JPL TR No. 32-468"** | Not held; wanted rank 1 (as 32-464). Report-number discrepancy recorded |
| Narin 1964, IITRI Report T-6, accessible-regions method | Not held. Added to the wanted list with Hunter (IITRI series context) |
| Stearns 1963, Navigation and Guidance in Space | Not held; textbook, not listed |
| Witting 1965, private communication | n/a |

Catalogue implications: none. The report is one-way, single-assist mission analysis.
