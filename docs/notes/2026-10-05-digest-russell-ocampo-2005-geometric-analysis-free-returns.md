# Digest: Russell & Ocampo 2005, "Geometric Analysis of Free-Return Trajectories Following a Gravity-Assisted Flyby" (#960)

R. P. Russell & C. A. Ocampo (University of Texas at Austin), "Geometric Analysis of Free-Return
Trajectories Following a Gravity-Assisted Flyby", J. Spacecraft and Rockets 42(1):138-152
(January-February 2005), doi 10.2514/1.5571 (wanted-list CONFIRMED; printed on every page).
- Presented as AAS 03-508 (Big Sky, August 2003).
- Filed as `cyclers_pdf/papers/russell-ocampo-2005-geometric-analysis-free-return-trajectories-gravity-assisted-flyby-jsr-42-1-138-doi-10.2514-1.5571.pdf`.
  15 pages, text layer, md5 720dd05609113414512b1dfde3540297.
- I read pp.1-3, 5-6 and 12-14 from the text layer, and checked Eqs. (17)-(19) on the page image (p.143).
  I skimmed the rest.

## 0. Verdict

This is the method paper for the v_inf-globe geometry behind the Russell-Strange and R1 generators. It
gives the full set of free returns after a flyby:
- **full-rev** (even-n pi) returns lie on spheres in velocity space;
- **half-rev** (odd-n pi) returns lie on circles (the "equal gamma cone");
- **generic** returns are found numerically, as "dots" on the v_inf sphere.
It is a method check for twobody-gen-opus; it computes no new two-working-body cycler.

## 1. Content (READ)

- **Full-rev return (Eq. 2, p.140):** M 2 pi sqrt(a_B^3/mu) = N 2 pi sqrt(a_F^3/mu), so
  a_F = a_B (M/N)^(2/3). Here N is the spacecraft revolutions and M the body revolutions.
- **Full-rev speed (Eqs. 12-13, p.142):** v_F = sqrt(2 mu/r - mu (N/M)^(2/3)/a_B). The direction is free,
  so the locus of departure velocities is a sphere of radius v_F centred at the base of v_B.
- **Full-rev circle (Eq. 17, p.143, page image):** z_F = (v_F^2 - v_inf^2 - v_B^2)/(2 v_B). z is along
  v_B, measured from the tip of v_B; the circle is where the v_inf sphere meets the full-rev sphere.
  - My check: |v_B + v_inf| = v_F with |v_inf| = v_inf gives exactly this.
- **Half-rev return (Eqs. 18-19, p.143, page image):**
  - v_Hr^2 = mu [2/(r1 + r2) - 1/a]
  - v_Htheta^2 = 2 mu r2/(r1^2 + r1 r2)
  - Rotating v_H about r_B gives a cone whose base is the half-rev circle. "Fast" transfers have negative
    initial radial speed and "slow" ones positive.
- **Probable print slip (p.143):** the Fig. 13 example (a = 1.18 AU, r1 = r2 = 1 AU, mu = 1) prints
  "v_Hr = 0.153 AU/TU and v_Htheta = 1 AU/TU". By Eq. (18), v_Hr^2 = 1 - 1/1.18 = 0.153, so
  v_Hr = 0.391 AU/TU. The printed number looks like v_Hr^2 (INFERRED). v_Htheta = 1 checks.
- **Lambert view (pp.139-141, Table 1):**
  - The vacant-focus sphere of radius 2a - r for r1 = r2.
  - Four transfers per point: direct and retrograde, fast and slow.
  - The rectilinear-ellipse solution is impractical.
- **Generic returns (pp.146-150):**
  - A grid search over the time of flight (about half-hour steps; over 180,000 solutions in 6 yr) and
    the angle psi, then 1-D refinement to machine precision.
  - The dots lie on non-uniform arcs in the ecliptic.
- **Cycler applications (pp.150-151, Table 4):** replacing chains of 1 pi and 2 pi transfers (from
  Russell & Ocampo 2004 JGCD) by single n pi transfers. Two cases:
  - 4.11.1-2: seven Earth free returns become one 5.5-yr half-rev return. The turn needed rises to
    107 deg, still above 200 km altitude.
  - 4.13.1-1: Mars v_inf drops from 9.3 to 5.6 km/s (a 40 % decrease) for 3 more transit days.
  - Table 4 also lists 4.11.1-4 as a new cycler (v_inf 3.8 km/s at Earth, 4.8 km/s at Mars, 168-d
    transit).

## 2. Gate relevance

- **`#942` R1 generator (method check):** Eqs. (2), (13), (17)-(19) are the closed forms for the full-rev
  and half-rev loci that a v_inf-globe generator should reproduce. The generic-return search is the
  numeric part.
- **R1(a)/(c):** no new two-planet cycler; Earth is the only host of free returns in the examples.

## 3. Positive controls

- Fig. 10 setup: mu = 1, r = a_B = 1 AU, N = 7, M = 4, v_inf = 0.5 AU/TU. Compute z_F from Eq. (17).
- Table 4 entries (sourced to Russell & Ocampo 2004 and this paper).

## 4. Citation mining (references 1-15)

Held: Patel et al. 1998 [1] (filed today); Hollister 1969 [4]; Menning 1968 [8]; Rall 1969 [9]; Byrnes,
McConaghy & Longuski 2002 [10].

Not held:
1. [11] Russell & Ocampo (2004), JGCD 27(3):321-335. Already in the wanted list (content held in the
   Russell 2004 dissertation).
2. [3] Wolf (1991), AAS 91-123 (also cited by Patel 1998).
3. [7] Uphoff, Roberts & Friedman (1976). Already in the wanted list. [6] Uphoff & Crouch (1993), JAS
   41(2):189-205, lunar cyclers: not held, not in the wanted list (low priority).
4. [2] Miele, Wang & Mancuso (2000), JAS 48(2-3); [5] Uphoff (1989); [14] Prussing (2000) and [15] Shen
   & Tsiotras (2003) are already in the wanted list (methods row).
