# Digest: Fornari & Pontani 2020, "A global exploration method to identify families of cycling Earth-Mars trajectories" (#960, #942)

E. Fornari and M. Pontani (Sapienza University of Rome), Aerotecnica Missili & Spazio 99:187-194 (2020),
doi 10.1007/s42496-020-00050-6. Received 19 Dec 2019, accepted 13 May 2020. 8 pp.
- File: supplied file `5c389825-fornari2020.pdf`, md5 8843391bdc5808b9e7ffc76ee3e90aeb, publisher text layer.
- Filed as `cyclers_pdf/papers/fornari-pontani-2020-global-exploration-families-cycling-earth-mars-trajectories-aerotec-missili-spazio-99-187-doi-10.1007-s42496-020-00050-6.pdf`.
- How I read it: the full text layer, then pp. 3-4 (Tables 1 and 2, Figs. 1-2, the Mars-flyby statements) on
  110-dpi page images. Both tables agree with the text layer cell for cell. The p.188 Mars statement and the
  p.191 S1L1 paragraph were read in the text layer and match the rendered pages.
- Tables transcribed: `data/sources/fornari-pontani-2020-cycler-families-tables.yaml`.
- Wanted-list row 64 (added 2026-10-07, "MUST ACQUIRE before any em row"): received.

## 0. Verdict

**A one-working-body census method in the Russell-Ocampo / Russell-Strange class. Mars is a massless target.
No collision with the `#942` em-1 to em-5 candidates.**
- **Mars does not bend the path.** Family I: "Mars flyby can be proven to have negligible effect on the
  spacecraft trajectory [1], therefore the heliocentric path in the interval [t0, t3] is a single elliptic arc"
  (p.188). For both tables, the printed Mars flyby altitude is the one that gives a 2-deg deflection: "the
  minimum flyby altitude that allows neglecting the effect of Mars flyby on the heliocentric ellipse". The
  effect "was assumed to be negligible, and this was verified for all of the cases" (pp.189-190).
- **Only Earth bends.** Family I is one ellipse per repetition k Tsyn, closed by one Earth flyby. Family II is
  two ellipses joined by one Earth slingshot. The governing equations (1)-(12) contain no Mars turn.
- **So it cannot produce the em candidates.** em-1/2/3 need two 38.4-deg Mars turns (a Mars generic return).
  em-4/5 need a Mars 1:1 return with 5.6-5.9-deg turns. Neither fits Family I or II.
- **Numbers do not collide either.** The only k = 3 row (Family II, l1 2, l2 3, m 1, s 4) has E/M v_inf
  7.537/6.382 km/s, against em-1/2 5.333/4.713 and em-3 4.684/4.539. No Table 2 row is near the R-O 2.5.1+0
  neighbour (k2 7.593/9.865).
- Reported to the lead and twobody-gen2-opus first (2026-10-07).
- **What it gives the project:**
  - A sourced, independent statement of the Aldrin cycler (k1 m1 l2: a = 2.3928e8 km, e = 0.3904, t1 = 145.33 d,
    v_inf 6.393/9.677 km/s) and of S1L1 (k = 2: a1 = 1.9506e8 km, e1 = 0.2554, a2 = 1.5682e8 km, e2 = 0.1609,
    v_inf 4.726/5.035 km/s).
  - S1L1 arc 1 then has perihelion 0.971 AU and aphelion 1.637 AU (our arithmetic). That matches the catalogue
    `s1l1-2syn-em-cpom` elements 0.97/1.64 (data/catalogue.yaml lines 518-519, Rogers 2012). The v_inf pair
    matches McConaghy's S1L1-B literal pair 4.7/5.0 quoted in the same row (lines 510-514).
  - A census count: 3339 Family I solutions for k = 1-150 and 23 Family II solutions for k = 1-10, under their
    filters (a <= 2.2 AU, flyby altitude below the SOI). Only the "remarkable" rows are printed.
- **Catalogue implication (PROPOSAL only):** `s1l1-2syn-em-cpom` could cite this paper as an independent
  circular-coplanar source for its elements and the 4.7/5.0 km/s pair. No new row: every printed solution is a
  one-working-body cycler of the kind the R-O/R-S census already covers, and no row is listed with a
  repeat structure we lack. A future check could test whether the Family I k = 8 and k = 15 rows (very low
  v_inf, e.g. k15 m22: 3.163/3.656 km/s with t1 = 211.94 d) are in the catalogue's R-O ladder.

## 1. Method

- Circular coplanar Earth and Mars, Keplerian arcs, Tsyn = 779.87 d (the paper's value; circular-orbit
  periods give 779.95 d).
- **Family I** (eqs. 1-4): integers k (synodic periods per repetition), m (cycler revolutions), l (Earth
  revolutions). Earth-side timing fixes the injection true anomaly f0 (eq. 3); then one nonlinear equation in e
  is solved by regula falsi on [e_Hohmann, 1]. The Mars crossing gives f1, t1 and the Mars phase. No Lambert
  solver is used.
- **Family II** (eqs. 5-12): integers k, l1, l2, m, s. Everything reduces to one equation in e1, with the Earth
  flyby feasibility condition eq. (12) (equal v_inf magnitude before and after).
- Filters: a <= 3.291e8 km (2.2 AU); perihelion below 1 AU and aphelion above Mars (Family I, and arc 1 of
  Family II); flyby altitude below the SOI radius.
- **Appendix:** when the required Earth turn needs a periapsis below the surface or atmosphere, two symmetric
  impulses at the SOI raise the perigee without changing the asymptotes. Aldrin: 179 m/s to reach 400 km.
  The k5 l1 7 l2 2 m4 s3 Family II row (90 km altitude): 137 m/s.

## 2. Checks (our arithmetic)

Re-solving eqs. (1)-(4) with Tsyn = 779.87 d, AU = 1.495978707e8 km, R_M = 1.523679 AU and
mu_Sun = 1.32712440018e11 km^3/s^2:

| row | a (ours / printed) | e | t1 d | v_inf E | v_inf M |
|---|---|---|---|---|---|
| k1 m1 l2 | 2.3928e8 / 2.3928e8 | 0.3904 / 0.3904 | 145.33 / 145.33 | 6.394 / 6.393 | 9.677 / 9.677 |
| k8 m12 l17 | 1.8876e8 / 1.8876e8 | 0.2119 / 0.2119 | 252.00 / 252.01 | 3.252 / 3.250 | 2.830 / 2.829 |
| k15 m16 l32 | 2.3749e8 / 2.3749e8 | 0.3707 / 0.3707 | 131.35 / 131.27 | 5.132 / 5.129 | 9.163 / 9.162 |
| k15 m22 l32 | 1.9206e8 / 1.9206e8 | 0.2216 / 0.2216 | 212.02 / 211.94 | 3.165 / 3.163 | 3.657 / 3.656 |
| k30 m43 l64 | 1.9503e8 / 1.9503e8 | 0.2350 / 0.2349 | 199.91 / 199.77 | 3.438 / 3.432 | 4.388 / 4.386 |

- The small t1 and v_inf offsets (up to 0.14 d and 0.006 km/s) come from the planet constants, which the paper
  does not print. With Tsyn from circular periods (779.95 d) the offsets are larger (up to 0.9 d), so the
  paper's 779.87 d is the right reading.
- **Two printed inconsistencies (kept as printed):**
  - p.188 says the Aldrin case has "k and l both set to 1". Table 1 prints k = 1, m = 1, l = 2, and l = 2 is
    what reproduces (one 779.87-d repetition spans two full Earth revolutions). The text means k = m = 1.
  - S1L1: the p.191 text gives t1 = 152.65 d; Table 2 prints 153.06 d. Not reconciled.
- Script: `fp_check.py`, output `fp_check.out`, filed beside the PDF.

## 3. Citation mining (8 references)

| Ref | Item | Status |
|---|---|---|
| [1] | Byrnes, Longuski & Aldrin 1993, JSR 30:334-336 | HELD (`byrnes-longuski-aldrin-1993-...`) |
| [2] | McConaghy, Landau, Yam & Longuski 2006, JSR 43:456-465 | HELD (`mcconaghy-landau-yam-2006-...`) |
| [3] | McConaghy, Longuski & Byrnes 2002, AIAA 2002-4420 | HELD (`mcconaghy-longuski-byrnes-2002-...`) |
| [4] | Pontani & Conway 2018, JGCD 41(2):360-376 | HELD (`pontani-conway-2018-...`) |
| [5] | Pontani 2019, "Optimal low-thrust hyperbolic rendezvous for Earth-Mars missions", Acta Astronautica 162:608-619 | not held; new candidate (rendezvous operations, low priority; not added) |
| [6] | Prussing & Conway 2012, Orbital Mechanics, 2nd edn (ch. 5) | not held; textbook |
| [7] | Pontani 2018, Spaceflight mechanics lecture notes | not held; not public |
| [8] | Pontani 2009, "Simple method to determine globally optimal orbital transfers", JGCD 32(3):899-914 | not held; method background, not added |

No reference is a missed cycler source.
