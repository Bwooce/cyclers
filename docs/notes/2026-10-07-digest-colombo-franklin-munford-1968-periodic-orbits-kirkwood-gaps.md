# Digest: Colombo, Franklin and Munford 1968, "On a Family of Periodic Orbits of the Restricted Three-Body Problem and the Question of the Gaps in the Asteroid Belt and in Saturn's Rings" (#960 batch 32)

G. Colombo, F. A. Franklin and C. M. Munford (Smithsonian Astrophysical Observatory), The Astronomical Journal 73(2):111-123
(March 1968), received 5 September 1967, final 8 January 1968. DOI 10.1086/110607, ADS 1968AJ.....73..111C. 13 pp.
- Source: `f36dcd45-1968AJ.....73..111C` (a PDF without extension), md5 fc46ec15f27fd5c89d30aece52cce3c1 (matches the brief).
  ADS scan with an OCR layer. The layer is good for the tables (digits agree with the images) and poor for equations.
- Proposed corpus filename:
  `colombo-franklin-munford-1968-family-periodic-orbits-restricted-three-body-problem-gaps-asteroid-belt-saturn-rings-aj-73-111-ads-1968AJ-73-111C.pdf`
  with companion `colombo-franklin-munford-1968-tables.yaml` (this batch).
- How I read it: all 13 page images at 110 dpi; the four tables (pp.114, 115, 117) at 170 dpi; p.2 and p.12 crops at 200 dpi.
  Tables I-IV parsed from the text layer by script (`parse_cfm.py`) and checked against the images and by the arithmetic below.

## 0. Verdict

**A table of symmetric periodic orbits of the planar Sun-Jupiter restricted problem (mu_J), one per initial radius, for sidereal
periods near 1/2, 2/3, 3/4 and 4/5 of Jupiter's.** It is not a pure-math paper. It gives 154 reproducible orbits, and the first
numerical map of the four "unconnected groups" behind the Kirkwood-gap argument.
- Orbit class: asteroid between Sun and Jupiter on the Sun-Jupiter line at t = 0 (inferior conjunction), velocity perpendicular,
  direct; at perihelion or aphelion at t = 0 and again at superior conjunction (T_s/2). These are the "Schwarzschild type" symmetric orbits.
- Relation to Bruno-Varin: the held Bruno-Varin JAMM digest (line 168) says families i_2 and i_3 "were computed at mu_J by Colombo
  and Franklin 1968". I read Tables II and III (2/3 and 3/4 groups) as those. I did not check Bruno-Varin's text for the match.
  **Unresolved:** whether Table I (1/2 group) is their i_1 (they say i_1 is closed at mu_M).
- **Proposal only:** file the scan and `-tables.yaml` under `data/sources/`. No catalogue row (heliocentric asteroid orbits, not cyclers).
  Useful as a positive control for any mu = 9.5e-4 symmetric-orbit corrector (sec. 3).

## 1. Model and conventions (pp.112, checked on image)

- G = 1, Sun mass 1, Jupiter mass **m = 0.000954786** (Jupiter-Sun ratio), circular, unit separation. Rotating frame rate `sqrt(1+m)`.
  Jupiter at xi = -0.99904612 at t = 0 (equals -1/(1+m)); the Sun at +m/(1+m). Period `T_J = 2 pi / sqrt(1+m) = 6.28019`.
  In standard form this is mu = m/(1+m) = 9.5388e-4 and lengths the same; velocities divide by `sqrt(1+m)`.
- Runge-Kutta, double precision, seven figures kept; Jacobi constant held to that accuracy. Floquet analysis used doubled precision.
- For each xi0 they solve for etadot0 so that the orbit crosses the xi-axis perpendicularly after half a synodic period.
  Sidereal period from `1/T* - 1/T_J = 1/T_s`. Elements are osculating about the Sun.
- Columns (all four tables): xi0, etadot0, T_s, e(0), a(0), e(T_s/2), a(T_s/2). Tables I, II, III, IV have 25, 33, 57, 39 rows.
  Table V (3 rows, p.117): the 2/3 orbit with Jupiter circular, at aphelion, at perihelion at t = 0 (Jupiter eccentricity 0.05).
- T_s is printed as a closed form for seven commensurate rows (e.g. `3 pi / 2(1+m)^(1/2)`); the grouping is ambiguous on the page.
  The reading `3pi/(2 sqrt(1+m))` fits the ordering and the labels 3/7, 2/5, 1/3, 1/4 in Fig. 1.

## 2. Results printed

- Four groups, unconnected (Fig. 2): "Hecuba" (a about 0.630), "Hilda" (0.763), "Thule" (about 0.82-0.83), "4/5" (about 0.87).
- **Gaps**: no periodic orbits for a between 0.630 and 0.645 (stable Hilda minimum near 0.645). At aJ = 5.2 AU that is
  3.28 to 3.36 AU (my arithmetic: 0.630 x 5.2 = 3.28, 0.645 x 5.2 = 3.35). Cluster at 0.76314 = (2/3)^(2/3), i.e. 3.97 AU. The abstract
  gives exactly those two numbers. The Hilda/Thule gap in a: largest Hilda a 0.765, smallest Thule a 0.770 (for e < 0.3).
- "No gaps are present around period ratios of 3/5, 3/7, 2/5, 1/3 or 1/4" in this orbit class (p.118); they suspect Jupiter's eccentricity is needed
  for those. The gaps at n/(n+1) are tied to Jupiter's circular motion. They apply the same argument to Saturn's rings (Cassini division at 1/2 of Mimas's period).
- **Stability (pp.121-122).** All Hecuba orbits stable; the lower-left near-horizontal Hilda branch unstable (critical turning
  point "near xi0 = -0.660, etadot0 = -0.54, e about 0.03"; Table II shows the turn at xi0 = -0.667, between rows 22 and 23).
  "All other periodic orbits in Fig. 2 are quite certainly stable" (not every one was checked). Floquet exponents for T* = 1/3, 2/5, 3/7, 3/5 of T_J
  are printed as purely imaginary, `+-i 2 pi/(2 T_s)`, `2 pi/(3 T_s)`, `2 pi/(4 T_s)`, `2 pi/(2 T_s)`, "to at least 14 figures".

## 3. My checks (scripts and outputs kept; none of these numbers is in the paper)

- **State vs elements** (`check_cfm_elements.py`): a(0), e(0) from the printed (xi0, etadot0) with the Sun-relative inertial velocity
  agree for all 154 rows to 2e-5 (the 5-digit rounding). A different frame rate (omega = 1) or GM = 1+m fails at 3e-4 to 1e-2.
  So the frame reading of sec. 1 is confirmed and no misprint of xi0, etadot0, e(0), a(0) exists.
- **Integration** (`integrate_cfm.py`, DOP853, rtol 1e-12, `integrate_cfm_all.out`): from each printed state to T_s/2,
  a(T_s/2) and e(T_s/2) agree with the print to 1.2e-5 (Table I), 7.7e-6 and 1.6e-5 (II), 6.2e-6 and 6.6e-6 (IV).
  Table III: 1.2e-4 and 1.4e-4, worst in the strongly unstable rows 35-53, where the 5-digit initial state limits the orbit
  (the true y = 0 crossing is 2e-3 off T_s/2 in III-44). **Only the state columns have a tight witness for Table III rows 35-53.**
  The Jacobi constant is conserved to 1e-11 on every orbit.
- **Symbolic-period rows close to 1e-7** (closure after one T_s): I-20 1.0e-7, I-21 1.8e-7, I-23 6.9e-8, I-25 1.9e-8, II-6 6.1e-7,
  II-18 2.1e-7. These use the printed 8-digit xi0 and so confirm those rows.
- **Jacobi constants (paper units / standard mu = m/(1+m) form, no mu(1-mu) constant):**
  I-5 2.791743 / 2.789080; I-20 3.270409 / 3.267290; I-23 3.470077 / 3.466767; II-6 2.870008 / 2.867270; IV-8 3.003499 / 3.000634.
  Ranges (standard form): Table I 2.306-3.779, II 2.341-3.165, III 2.767-3.049, IV 2.917-3.018.
- **Table V** row 1 equals Table II row 6 (0.76393, 0.45373, 0.76321, 0.45491), as it must.
- **Stability recomputed** (`floquet_cfm.py`, `stability_cfm.py`; positive control: at m = 0 the multipliers equal the closed forms
  exactly, -1, -1, +-i, e^(+-2 pi i/3); finite-difference monodromy agrees with the variational one to 1e-6 relative).
  - At m = 9.5e-4 the closed forms are only first order: I-20 phase 0.50107 pi (not 0.5 pi), I-21 0.66746 pi (not 0.66667 pi).
    **The "14 figures" statement is not reproduced.** For T* = 1/3 (row I-23) the multipliers are -0.994994 and -1.005031
    (trace -2.000025, marginally hyperbolic, sensitive to m: at m = 5e-4 the trace is -2.000007). For T* = 3/5 (row II-18)
    they are -0.935 and -1.069, clearly unstable, which agrees with the paper's own statement that part of the Hilda family is unstable.
  - Trace test on every row (stable if |trace| < 2): Table I stable except I-23; Table II unstable at rows 1, 18, 24-33;
    Table III unstable at rows 2-3, 22, 27-56; Table IV unstable at rows 1-4, 12-15, 17-38 (traces up to 1e3, orbits that pass
    within about 1.5 Hill radii of Jupiter). The two isolated rows 1 and 18 in Table II, and the short runs, sit at trace -2.
    **This contradicts the paper's "all other orbits quite certainly stable" for the Thule and 4/5 outer branches.** The paper
    checked only sample orbits; I have not found a cause on my side but I did not check against a second implementation.

## 4. What it means for the project

- A ready set of orbits at mu about 9.5e-4 for testing a symmetric-orbit corrector or continuation (use the rows with 8-digit xi0).
- A check of the p/(p+1) families in Bruno-Varin and Henon's generating families at a real mass ratio.
- It is not an ejection-collision or torus paper. Only link to #899: none.

## 5. Citation mining

- Colombo and Franklin 1968 (this paper) and Franklin and Colombo 1970 (Icarus 12:338-347): not held; both are wanted-list row 56.
- Kotoulas and Voyatzis 2004, Bray and Goudas 1967: not held; row 56. Voyatzis et al. 2005: not held; row 56.
- Schwarzschild 1903 (AN 160:385), Dziewulski 1909 (AN 183:65), Hagihara 1957 (Stability in Celestial Mechanics): not held, not listed.
- Brown, Goddard and Kane 1967 (ApJS 14:57, asteroid distribution), Message 1966 (IAU Symp. 25:197): not held, not listed.
- Poincare 1899 (Methodes nouvelles): not held. Related held works: Henon 1997 generating families (HELD), Hitzl-Henon 1977 (HELD),
  Hadjidemetriou 1975 (HELD, continuation and stability).

*Filed as `cyclers_pdf/papers/colombo-franklin-munford-1968-family-periodic-orbits-restricted-three-body-gaps-asteroid-belt-saturn-rings-aj-73-111-doi-10.1086-110607.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. The table transcription is `data/sources/colombo-franklin-munford-1968-tables.yaml`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
