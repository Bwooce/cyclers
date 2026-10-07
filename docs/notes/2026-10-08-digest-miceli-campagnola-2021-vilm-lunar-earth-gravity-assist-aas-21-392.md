# Digest: Miceli & Campagnola 2021, "Effect of V-infinity Leveraging with Lunar-Earth Gravity Assist on Interplanetary Trajectories"

G. E. Miceli (ISAE-SUPAERO, Toulouse) & S. Campagnola (JPL), AAS/AIAA Space Flight Mechanics Meeting, virtual,
1-4 Feb 2021, **paper AAS 21-392** (the preprint prints "AAS XX-XXX"; the number and venue come from a web search that
found NTRS 20220003813 and the ResearchGate record, which both give AAS 21-392; I could not find a Crossref DOI for it).
JPL clearance CL#21-0437. No DOI found.

- Source file: upload `f1d5ddb2-CL21_0437.pdf`, 20 pp, md5 f195b852bbb3b525d3c82b106d3444c5. Digital (publisher
  LaTeX text layer, good). No OCR needed.
- Proposed filename: `miceli-campagnola-2021-v-infinity-leveraging-lunar-earth-gravity-assist-interplanetary-trajectories-AAS-21-392.pdf`
- Compared with nothing held (not in corpus; not in CORPUS_INDEX; the index hits for "Miceli" are the Miceli-Bosanac
  Neptune papers, a different author pair).
- How I read it: text layer read in full; equations 2-9 and 10-12 and Tables 1-3 read on page images (pp.4, 5, 6, 18, 19,
  110-130 dpi); tables compared cell by cell with the text layer (`witness-comparison.tsv`: 17 rows x 4 columns = 68
  cells, 68 agree); then every table row checked by arithmetic (`checks.py`, `checks.out`, `orbit_check.py`). Figures
  3-31 are curves only (no numbers); I looked at no figure for a value.

## 0. Verdict

A short patched-conic parametric study. It reuses the Sims-Longuski-Staugler 1997 V-infinity-leveraging (VILM) set-up
(held), keeps the Moon phase free, and adds a lunar flyby before or after the Earth perigee. It prints **three small
tables** (17 rows) for a Jupiter-bound mission, and nothing else numeric except the headline percentages.

What it gives the project:

1. **A VILM cost-axis anchor for #983.** Table rows give C3, delta-V_launch+VILM and rp_post for six (E)n:m families.
   The paper states its own cost model (eqs 7-9), and **my check shows the printed numbers obey it to 1e-4 km/s** for
   all 17 rows (section 3). That makes the tables a usable golden for a re-implementation of eqs 7-9, with the caveat
   that the tables do not print the leveraging radius, so the orbit itself must be rebuilt (section 4).
2. **A size for the lunar effect.** Of the six family-and-transfer cases in the tables, five have a plain-EGA baseline. The Moon-assisted
   Earth flyby saves 0.12 to 0.60 percent of delta-V_tot in four of them, but 6.35 percent (and 21.06 percent of the leveraging delta-V) for 3:2+ Long. That is the number
   for the Earth-Moon lane (#1000): the effect is large only when the plain Earth flyby is weak.
3. **Honest limits**, all from the paper itself: planar, circular, coplanar Earth and Moon orbits; Moon phase free
   (no phasing constraint); interior-leveraging curves partly interpolated and partly wrong (Fig. 20, p.13).

PROPOSALS only (nothing done): (a) use Tables 1-3 as a read-only regression for any VILM-plus-lunar-flyby cost
function in the X1 tooling, tolerance 1e-4 km/s on delta-V_tot given C3 and delta-V_lev; (b) do **not** use rp_post or
the percentages as goldens without the orbit definitions; (c) no catalogue row follows from this paper.

## 1. Formulation (read on page images, pp.2-6)

Patched conics, ecliptic plane, Earth and Moon on circular coplanar orbits, Moon phase free (pp.2, 6).

- Notation (eq. 1, p.4): `(EI)n : m sigma_k`, E = +1 exterior, I = -1 interior; n secondary-body revolutions; m
  spacecraft revolutions; sigma = +1 Long (H+) or -1 Short (H-); k full revolutions of the non-tangent arc. Special case
  of the Campagnola-Strange-Russell 2010 general form, with the first arc tangent (departure at an apse, r_va1 = 1 AU).
- Families (p.4): exterior 1:1(+ only), 4:3, 3:2, 2:1, 3:1, 4:1 (each with +/-1); interior 4:5, 3:4, 2:3, 3:5, 4:7, 1:2.
  Same six per type as Sims 1997. Short-transfer exterior has five families (1:1- not computable, p.8).
- Phasing (eqs 2-5): `d_theta1 + d_theta2 = 2 pi (k1 + k2 + 1 - n)`, with d_theta(ra, rp) from eq. 3; with r_va1 = 1 and
  sigma_1 = 0 it reduces to eq. 5 with the single unknown r_va2, solved by MATLAB trust-region-dogleg (p.5).
- Cost (eqs 7-9, p.5): `dV_tot = dV_launch + dV_VILM`; `dV_launch = | sqrt(v_inf^2 + 2 v_c^2) - v_c |` with v_c the
  circular speed at **185 km** altitude (p.5); `dV_VILM = | v_le2 - v_le1 |` by vis-viva at the leveraging radius.
- Lunar flyby (eqs 10-12, p.6): `v_sc/Earth = v_moon + v_inf,moon`; the Moon rotates v_inf,moon by delta_m. The paper
  prints `delta_m = 2 sin( mu_m / (mu_m + r_p |v_inf,m|^2) )` (eq. 11). **The usual form has arcsin**
  (`delta = 2 arcsin(1/(1 + r_p v_inf^2/mu))`), so eq. 11 as printed has lost the inverse sine (typesetting slip; I
  read it on the page image at 130 dpi). Minimum Moon-flyby altitude 300 km (cited to Schoenmaekers, IAC-14); Earth
  flyby perigee altitude from 200 km to the Moon's SOI; direct and retrograde hyperbolae; "LEGA" = Moon met before
  Earth perigee, "ELGA" = after. For each Earth hyperbola and each delta_m in [-max, +max] the best final apse is kept
  (max aphelion exterior; min perihelion interior).

## 2. Tables 1-3 (pp.18-19), transcribed

All 17 rows, as printed (comma decimals in the paper; I write points). Columns: rp_post (AU, perihelion radius "before
the GA (Long) or after the GA (Short)", p.17), DV_tot (km/s), DV_lev (km/s), C3 (km2/s2). Short = family "-", Long =
family "+". Two witnesses: image reading and text layer, 68 of 68 cells agree. The EGA row of 3:2- is "/" (the plain
Earth flyby cannot reach Jupiter there, p.17). Machine-readable copy: `miceli-campagnola-2021-jupiter-vilm-lega-tables.yaml`.

| family | row | rp_post | DV_tot | DV_lev | C3 |
|---|---|---|---|---|---|
| 2:1- | LEGA | 0.99914 | 4.8606 | 0.50389 | 26.1524 |
| 2:1- | ELGA | 0.99916 | 4.8606 | 0.50389 | 26.1524 |
| 2:1- | EGA | 0.99709 | 4.89 | 0.53221 | 26.1793 |
| 2:1+ | LEGA | 0.91701 | 4.8811 | 0.46682 | 27.5539 |
| 2:1+ | ELGA | 0.91636 | 4.8853 | 0.47057 | 27.5665 |
| 2:1+ | EGA | 0.91311 | 4.9067 | 0.48937 | 27.6296 |
| 3:1- | LEGA | 0.999999 | 5.391 | 0.16409 | 48.0549 |
| 3:1- | ELGA | 0.999999 | 5.391 | 0.16409 | 48.0549 |
| 3:1- | EGA | 0.999999 | 5.3977 | 0.17058 | 48.0617 |
| 3:1+ | LEGA | 0.96415 | 5.4079 | 0.1596 | 48.6127 |
| 3:1+ | ELGA | 0.96415 | 5.4079 | 0.1596 | 48.6127 |
| 3:1+ | EGA | 0.96283 | 5.4146 | 0.16557 | 48.6319 |
| 3:2- | LEGA | 0.97138 | 5.1813 | 1.3309 | 14.1044 |
| 3:2- | ELGA | 0.9754 | 5.1792 | 1.329 | 14.1004 |
| 3:2- | EGA | / | / | / | / |
| 3:2+ | LEGA | 0.83321 | 5.0074 | 1.1099 | 15.2026 |
| 3:2+ | ELGA | 0.83266 | 5.0118 | 1.1138 | 15.2159 |
| 3:2+ | EGA | 0.79164 | 5.3468 | 1.406 | 16.2177 |

Headline numbers in the text, checked against the table (ours, `checks.out`):

| statement (page) | printed | recomputed from the table |
|---|---|---|
| 3:2+ DV_tot saving, LEGA vs EGA (abstract p.1; p.19) | 6.4 percent | **6.35 percent** (6.27 for ELGA) |
| same, p.18 | "6.5 percent on DV_tot" | 6.35 percent |
| 3:2+ DV_lev saving (abstract; p.18) | 21.06 percent | **21.06 percent** (ELGA 20.78) |
| 2:1- DV_VILM saving (p.18) | 5.3 percent | **5.32 percent** |
| conclusion p.19: "DV_VILM ... cut of 6.4 percent at the most" | 6.4 percent on DV_VILM | wrong quantity: 6.35 is DV_tot; DV_VILM is 21.06 |

So the paper states 6.4 (abstract, p.19), 6.5 (p.18) and attaches the 6.4 to two different quantities. The table gives
6.35 percent on DV_tot and 21.06 percent on DV_lev. Use the table values.

Other table facts (ours): 3:1 rows are the same for LEGA and ELGA (the Moon does not matter at high v_inf, p.17);
3:2- has a solution only with the lunar flyby (EGA "/"); in Long transfers LEGA beats ELGA, in Short transfers
ELGA beats LEGA for 2:1 and 3:2 (consistent with p.17 and the conclusion).

## 3. Checks (our arithmetic; scripts and output kept)

`checks.py` / `checks.out`:

- **Cost model closes on every row.** With v_c from mu_Earth = 398600.4418 km3/s2 and r_c = 6378.137 + 185 km
  (v_c = 7.79315 km/s), `dV_launch = sqrt(C3 + 2 v_c^2) - v_c`, then `dV_launch + DV_lev` equals the printed DV_tot
  with a largest difference of 1e-4 km/s over all 17 rows (14 of 17 rows differ by under 0.00005; three rows, 3:2- LEGA, 3:2+ LEGA and 3:2+ EGA, differ by 0.0001, which is
  rounding of the printed digits). This is **not circular**: C3, DV_lev and DV_tot are all printed; I only supplied
  the Earth constants and the 185 km radius from p.5. It confirms the tables' digits in C3, DV_lev and DV_tot are
  mutually consistent, which is a third witness for those three columns. It cannot check rp_post.
- Launch orbit from C3 (`orbit_check.py`, `orbit_check.out`): treat v_inf = sqrt(C3) as tangent at 1 AU
  (Earth circular, v_E = 29.7847 km/s). The implied heliocentric periods are 2.01-2.07 yr (2:1), 3.00-3.04 yr (3:1),
  1.60-1.67 yr (3:2). So the printed launch orbits are **near** the stated resonances but not exactly on them (3:2 is
  about 1.6 yr, not 1.5). The paper says r_le1 "varies to span a range of v_inf within a family" (p.5), which fits.
  The 3:2+ EGA row (C3 16.2177) implies a 1.667 yr orbit, the LEGA row (15.2026) 1.633 yr. Observation only; I do not
  claim the paper is wrong.
- Eq. 3 and eq. 5 were read on page images; I did not re-derive or integrate them. I did not run a phasing solve.

## 4. What is NOT in the paper (so what a golden cannot do)

- No r_le1, r_va2, rp (Earth flyby altitude), perilune altitude, Moon phase or flyby date for any row. The tables
  give only outputs. A re-implementation can match C3 to DV_tot (section 3) but cannot match rp_post without
  searching for r_le1.
- No inclination, no 3D, no ephemeris. "Effect of LEGA on orbit inclination" is future work (p.19).
- No full-ephemeris or flown comparison. The only external link is the JUICE LEGA idea (Schoenmaekers, not held).
- Figures 3-31 hold the family curves (final aphelion or perihelion versus DV_tot, DV_launch and DV_VILM) but no
  digits; I took no value from any figure.
- Interior leveraging: initial parts of the Long-transfer curves were filled by interp1 (spline or pchip), Fig. 20
  "shows incorrect results" (p.13), so no interior number is safe.

## 5. Citation mining (16 references)

Held = `ls papers | grep -i` AND CORPUS_INDEX. Wanted-list rows: none of the not-held items has a row today (the only
related row is #48). Proposed new rows use numbers from 64.

| ref | work | status |
|---|---|---|
| 1 | Kondratyuk, NASA TT F-9285 (1965), cited as "Mel'kmov" | not held; low |
| 2 | Crocco 1956, VII IAC Rome | **held** (crocco-1956-...-proc-vii-iac-rome-227-252.pdf) |
| 3 | Dowling, Kosmann, Minovitch, Ridenoure 1990, IAF-90 | not held; no row; Minovitch lane is covered elsewhere |
| 4 | Williams 1990, Purdue MS thesis | not held; no row; proposed row 64 (low) |
| 5, 6 | Sims, Longuski & Staugler 1997, JGCD 20(3):409, doi 10.2514/2.4064 | **held** (sims-longuski-staugler-1997-...-checks files present); [5] and [6] are the same paper listed twice |
| 7 | Campagnola & Russell 2010, JGCD 33:463, doi 10.2514/1.44258 | **journal version not held**; the AAS 09-224 Part A version is held (campagnola-russell-2009-endgame-partA-...). Proposed row 65 |
| 8 | Johannesen & D'Amario 1999, AAS 99-360 | not held; no row; proposed row 66 (low) |
| 9 | Ross & Scheeres 2007, SIADS 6(3):576 | **held** (ross-scheeres-2007-...-SIADS-6-3.pdf) |
| 10 | Boutonnet, De Pascale, Canalias 2008, Laplace mission, IAC 08-C1.6 | not held; no row; low |
| 11 | McAdams et al. 2006, MESSENGER, JSR 43(5):1054, doi 10.2514/1.18178 | not held; no row; low (only the McAdams-2011 Uranus AAS 11-188 paper is held, different) |
| 12 | Jehn, Campagnola, Garcia, Kemble 2004, ISSFD-18, ESA SP-548 | not held; no row; low |
| 13 | Langevin 2000, Acta Astronaut. 47:443 | not held; no row; low |
| 14 | **Schoenmaekers 2014, "Improved interplanetary transfers with Lunar-Earth Gravity Assists", IAC-14-C1.9.12** | **not held, no row; the key reference for the Earth-Moon lane (JUICE LEGA, 300 km Moon altitude rule). Proposed row 67, priority high for #1000** |
| 15 | Strange, Campagnola & Russell 2009, AAS 09-435 | not held; **already row 48** |
| 16 | Campagnola, Strange & Russell 2010, CMDA 108:165, doi 10.1007/s10569-010-9295-1 | **journal version not held**; the AAS 10-164 conference version is held (campagnola-strange-russell-2010-fast-tour-design-non-tangent-...-AAS-10-164.pdf). Proposed row 68 (journal version, to compare with AAS) |

Proposed wanted rows (text for the list):
- 64 | Williams, S. N. (1990), "Automated design of multiple encounter gravity-assist trajectories", MS thesis, Purdue | cited by Miceli-Campagnola 2021 as origin of VILM | low
- 65 | Campagnola, S. & Russell, R. P. (2010), "Endgame Problem Part 1: V-Infinity-Leveraging Technique and the Leveraging Graph", JGCD 33(2):463-475, doi 10.2514/1.44258 | journal version of held AAS 09-224 | medium (X1 VILM)
- 66 | Johannesen & D'Amario (1999), "Europa Orbiter Mission Trajectory Design", AAS 99-360 | low
- 67 | Schoenmaekers, J. (2014), "Improved interplanetary transfers with Lunar-Earth Gravity Assists", IAC-14-C1.9.12 | the only published LEGA design source cited; JUICE | high (#1000)
- 68 | Campagnola, Strange & Russell (2010), CMDA 108:165-186, doi 10.1007/s10569-010-9295-1 | journal version of held AAS 10-164 | medium (X1, non-tangent VILM; the phasing equations eq. 2-5 of this digest come from it)

*Filed as `cyclers_pdf/papers/miceli-campagnola-2021-v-infinity-leveraging-lunar-earth-gravity-assist-interplanetary-trajectories-aas-21-392-jpl-cl-21-0437.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/miceli-campagnola-2021-jupiter-vilm-lega-tables.yaml`.*
