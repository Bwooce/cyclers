# Digest: Patel, Longuski & Sims 1998, "Mars Free Return Trajectories" (#960)

M. R. Patel, J. M. Longuski & J. A. Sims, "Mars Free Return Trajectories", J. Spacecraft and Rockets
35(3):350-354 (May-June 1998), doi 10.2514/2.3333 (wanted-list CONFIRMED; printed on every page).
- Filed as `cyclers_pdf/papers/patel-longuski-sims-1998-mars-free-return-trajectories-jsr-35-3-350-doi-10.2514-2.3333.pdf`.
  5 pages, text layer, md5 228278211701a354f78585932e426977.
- I read all pages from the text layer. Tables 2 and 3 are scrambled there, so I read them on the page
  images (pp.353-354).

## 0. Verdict

A survey of all Earth-Mars-Earth free returns with launch dates 1995-2020, launch v_inf 4-8 km/s and
time of flight under 4 years, made with a patched-conic automated search (STOUR lineage). It also links
them to Henon's consecutive-collision orbits. It is a usable cross-check for the R1 one-body (Earth
free-return) generator:
- the three time-of-flight bands (about 1.5, 2 and 3 yr);
- Henon's Table 3 collision orbits (a, e, period);
- the 15-yr (7 synodic period) repeat.
No cycler is computed; the escalator data are quoted from Byrnes, Longuski & Aldrin 1993 (held).

## 1. Content (READ)

- **Repeat (p.351):** the Mars period is about 1 7/8 yr, so the Earth-Mars synodic period is about 2 1/7
  yr, and the inertial positions repeat in seven synodic periods (about 15 yr).
- **Free-return families (pp.351-352, Figs. 1-6):**
  - TOF bands near 1.5, 2 and 3 yr.
  - The 2-yr band has a second window offset by about 0.6 yr (Mars met before or after aphelion; the
    escalator geometry).
  - Fast free returns (TOF about 1.4 yr) in 2000 and 2002 repeat in 2015 and 2017.
  - Encounter cases 1-8 (E1/E2, M1/M2: before or after periapsis).
  - Resonant periods (Wolf 1991): n x orbital period = m x Earth period, j x orbital period = k x Mars
    period. VISIT-I [4 5 3 2] and VISIT-II [2 3 5 4].
- **Table 2 (p.353), up escalator from Byrnes et al. 1993 (page image):**
  - Earth-1 19 Nov 1996, 6.19 km/s.
  - Mars-2 1 May 1997, 10.69.
  - Earth-3 1 Jan 1999, 5.94.
  - Mars-4 28 May 1999, 11.74.
  - Earth-5 8 Feb 2001, 5.67.
  - ... Earth-15 13 Nov 2011, 5.81.
  - Maneuvers of 0.45-0.74 km/s.
  - The down-escalator columns were not re-read on the image; use Byrnes et al. 1993 (held) for
    either.
- **Table 3 (p.354, "Excerpt from Henon's tables", page image):**

| tau/pi, yr | eta/pi | a, AU | e | Period, yr |
|---|---|---|---|---|
| 2.00000 | 1.00000 | 1.58740 | 0.37004 | 2.00141 |
| 3.00000 | 2.00000 | 1.31037 | 0.23686 | 1.50003 |
| 3.00000 | 1.00000 | 2.08008 | 0.51925 | 3.00083 |

  - Check: a^1.5 gives 2.0000, 1.5000 and 3.0000 yr. The printed periods differ by up to 0.07 %
    (probably Henon's mass-ratio units, INFERRED; Henon 1968 is held for a check).
- **Collision orbits (pp.353-354):** Henon's timing equation (1) and a, e from Eqs. (2)-(3). Free returns
  to Mars need the collision orbit's aphelion beyond Mars's orbit. Large Mars-flyby-altitude members are
  Henon's collision orbits; fast ones are strongly bent by Mars.

## 2. Gate relevance

- **`#942` R1 one-body generator:** a cross-check on the family structure (TOF bands, 15-yr repeat) and
  on Henon collision-orbit a and e.
- **R1(a):** the paper notes that Mars gravity strongly perturbs the fast free returns (Fig. 3,
  opposition class). These are one-shot trajectories, not cyclers.

## 3. Positive controls

- Table 3 (three collision orbits).
- The up-escalator Table 2 dates and speeds (sourced to Byrnes et al. 1993; prefer the held original).

## 4. Citation mining (references 1-26)

Held (by filename check): Byrnes, Longuski & Aldrin 1993 [10]; Friedlander et al. 1986 [8]; Henon 1968
[24]. Not held: Breakwell & Perko 1965 [14] (already in the wanted list).

Not held (selected):
1. [23] Wolf, A. A. (1991), "Free Return Trajectories for Mars Missions", AAS 91-123. The resonant
   free-return framework. Not in the wanted list; worth adding at low priority.
2. [25] Howell, K. C. (1987), "Consecutive Collision Orbits in the Limiting Case mu = 0 of the Elliptic
   Restricted Problem", Celest. Mech. 40:393-407. Not held (filename check).
3. [26] Prado & Broucke (1994), JGCD 17(5):1075-1081. Henon's transfer problem via Lambert. Not held.
4. [7] Niehoff 1986 (AAS 86-172) is already in the wanted list. [11]-[22] are STOUR-lineage papers and
   theses; background.
