# Digest: Hatten & Russell 2015, "Comparison of three Stark problem solution techniques for the bounded case" (#960 batch 37, #999, #864)

N. Hatten and R. P. Russell (University of Texas at Austin), Celestial Mechanics and Dynamical Astronomy
121(1):39-60 (2015), doi 10.1007/s10569-014-9586-z (Crossref: CMDA 121(1), pp.39-60, 2015-01; title and authors
match the title page). Received 31 Mar 2014, revised 15 Jul 2014, accepted 11 Sep 2014, online 7 Oct 2014 (text layer).
The supplied file name says "2014" (online year); the issue year is 2015.
- Supplied file `c3189fcf-hatten2014.pdf`, 22 pp. = journal pp.39-60 (PDF page n = p.38 + n). md5
  3d64c6ad9f96322b5fbc79fcfc8a2f89. Springer text layer.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/hatten-russell-2015-comparison-three-stark-problem-solution-techniques-bounded-case-cmda-121-39-doi-10.1007-s10569-014-9586-z.pdf`
- **How I read it:** the whole text layer. On page images (130-200 dpi): p.45 (sec. 5.1 test definition),
  p.50 (Tables 2-3, rows to the theta_1 line), p.51 (Table 4 and the transfer definition), p.52 (the 24 SPO
  sentence), p.54 (conclusion), p.55 (eqs 11-25, Lantoine bounded-case summary). Figures 2, 3 and 10 are read
  only qualitatively.
- **Checks:** its eq. 25 (the separation constant c) is conserved and equals the c of Lantoine eq. 8 (see
  `check_stark.out` sec. B and the Lantoine digest). Online Resource 1 (their Fortran Lantoine and Biscani
  routines) was not supplied.
- Wanted list: not listed (supplied by the owner).

## 0. Verdict

**A head-to-head, same-language (Fortran, gfortran 4.8.2, p.45 text layer) comparison of three Stark solvers for bound
orbits: Lantoine-Russell (Jacobi elliptic), Biscani-Izzo (Weierstrass elliptic) and Pellegrini-Russell-
Vittaldev (F and G Taylor series). Result: the series is fastest for short segments; the closed forms win only
for long arcs. It publishes test definitions and timings, not state values.**
- **Single segment (sec. 5.1):** mu = 1 LU^3/TU^2, eps = 1e-3 LU/TU^2, a0 = 1 LU, i0 = 28.5 deg,
  Omega0 = 5 deg, omega = 10 deg, e0 from 0 to 0.9, start at periapsis, truth = quad-precision LSODE with
  relative tolerance 1e-25 (p.45, image). Series orders 5-18.
  - Lantoine is faster than Biscani in their implementations (Tables 2-3: per-step elliptic work 4.15e-6 s vs
    9.15e-6 s (the Biscani sum from the text layer); the cost driver is Pi (third kind) for Lantoine and theta_1 for Biscani).
  - "The Pellegrini method is found to be the fastest of the three for series orders up to approximately 16"
    (p.54, image). For relative position tolerance 1e-12 the series saves time "for propagations of up to
    approximately 60 deg of eccentric anomaly if the series is expanded in eccentric anomaly, even for
    significantly eccentric orbits" (p.54, image).
  - The series propagated in physical time t diverges at e = 0.9 with Delta M = 5 deg (p.45, text layer). Steps of equal
    eccentric anomaly (Sundman alpha = 1) are the robust choice.
  - The closed forms lose accuracy through transcendental-function evaluation and through the t(tau) inversion
    (Table 5, p.53, text layer). Hatten caps the elliptic modulus m at 1 minus two machine epsilons when roundoff pushes it
    above 1 (error ~1e-15, p.42, text layer).
- **Many segments (sec. 5.2):** a Q-law low-thrust transfer, constant mass, no coasts, thrust 0.01 m/s^2,
  1 LU = 6378 km, initial a = 7000 km, e = 0.3, i = Omega = omega = 1 deg, nu = 0; target e = 0.7, i = 25 deg
  (Table 4, p.51, image); r_p,min = 6578 km; tolerances e 0.01, i 1 deg; 8 to 96 segments per orbit (SPO).
  Analytic cost grows linearly with SPO; series cost is nearly flat because the needed order drops as SPO
  rises (Fig. 10). "The utility of the analytical solutions as orbit propagation tools is superseded by the
  Pellegrini solution if at least ~24 SPO are deemed necessary" (p.52, image).
- **Implementation facts worth keeping:** the 3D Lantoine formulas do not reduce to the 2D ones for planar
  initial conditions, so a full solver needs both (p.42); for bounded motion only xi1-eta2 (2D) and (xiI, eta)
  (3D) are needed (p.42); boundedness is decidable analytically from the initial state (footnote 1, p.40).
  Biscani-Izzo publish no 2D formulas (p.43).
- **Correction to Lantoine 2011:** Hatten's eq. 25, c = x'(x y' - y x') + mu y/r - (1/2) eps x^2 (field along +y,
  p.55, image), is the conserved form; Lantoine's printed eq. 9 is not conserved and equals -c - eps x^2 (our check; see the
  Lantoine digest for the two one-step repairs).
- **What it gives the project:** the decision data for `#999`'s propagator choice (see `stark-plan.md`):
  short constant-thrust segments favour a Taylor series in eccentric anomaly; long single arcs favour the
  closed form. Two reusable test definitions (sec. 5.1 orbit grid; Table 4 transfer), but no published end
  states, so no state goldens.
- **Catalogue implication (PROPOSAL only):** none.

## 1. Content

- **Sec. 1 (pp.39-41).** Sims-Flanagan (Kepler arcs plus impulses) vs Stark (Kepler plus constant thrust per
  segment) vs full integration (Fig. 1). The Stark model is exact for piecewise-constant thrust; earlier
  speed-ups vs numerical integration: Lantoine & Russell 2009, Lantoine 2010 thesis.
- **Sec. 2 (pp.41-42), Lantoine method.** r'' = -mu r/r^3 + eps k-hat (eq. 1); dt = 2r dtau (eq. 2), i.e.
  Sundman c_s = 2, alpha_s = 1 (tau proportional to eccentric anomaly). Elliptic functions and integrals by
  Fukushima's routines.
- **Sec. 3 (pp.42-43), Biscani method.** Weierstrass p, p', p^-1, sigma, zeta; no case split, one time
  variable, 3D direct; complex intermediates; Python scripts (github bluescarni/stark_weierstrass) re-coded in
  Fortran.
- **Sec. 4 (pp.43-44), Pellegrini method.** r = F r0 + G v0 + H p (eq. 3); any Sundman variable; no
  elliptic functions; no case split; 2D and 3D in one formula; truncation error estimate available; drawback:
  symbolic coefficients and large files; most efficient at orders 12-18 (p.44, text layer); impractical beyond about 25.
- **Sec. 5.1 (pp.45-50).** Fig. 2 (error vs CPU for Delta M = 5 deg), Fig. 3 (CPU vs Delta M and Delta E,
  5-90 deg, minimum order for rtol 1e-12, eq. 5), Figs 4-7 (count of elliptic evaluations), Tables 2-3
  (per-function timings over 1e6 trials).
- **Sec. 5.2 (pp.50-53).** Q-law transfer, Figs 8-10.
- **Sec. 6, Table 5 (p.53).** Summary: speed (Lantoine medium, Biscani slow, Pellegrini medium to fast by
  order); code size (medium, small, large); independent variable (tau, tau2; tau; any Sundman variable);
  error sources.
- **Appendices (pp.54-58).** Lantoine bounded 2D (eqs 6-34) and 3D (eqs 35-44) relations; Biscani (eqs
  45-61); Pellegrini recursion (eqs 62-68).

## 2. Golden candidates (published values only)

- None are state values. Usable published definitions:
  - **G-H1 orbit grid:** sec. 5.1 elements (above), eps = 1e-3, Delta M = 5 deg and Delta M or Delta E up to
    90 deg; the truth must be our own high-precision integration. Published claim to reproduce: at rtol 1e-12 the series in
    eccentric anomaly is cheaper than the closed forms for Delta E up to about 60 deg (Fig. 3, p.54).
  - **G-H2 transfer:** Table 4 plus the sec. 5.2 rules (Q-law nominal gains; not reproducible exactly without
    Petropoulos 2005's gains). The end state is not published.
- Timings (Tables 2-3, Figs 2-3, 10) depend on hardware and implementation; not goldens.

## 3. Citation mining

| Cited work | Held? | Wanted list |
|---|---|---|
| Biscani & Izzo 2014, MNRAS 439:810-822, doi 10.1093/mnras/stt2501 (Weierstrass solution) | not held | not listed; **candidate row** (the third method; also its Python code) |
| Biscani & Izzo 2013, Python code (github bluescarni/stark_weierstrass) | code, not held | - |
| Lantoine & Russell 2011, CMDA 109:333 | supplied this batch | - |
| Pellegrini, Russell & Vittaldev 2014, CMDA 118:355 | supplied this batch | - |
| Lantoine & Russell 2009, ISSFD 21 | not held | not listed; candidate row (see Lantoine digest) |
| Lantoine 2010, PhD thesis, Georgia Tech, pp.125-181 | not held | not listed; optional |
| Petropoulos 2005, Q-law refinements, AAS/AIAA | not held | not listed; needed only to repeat G-H2 |
| Sims & Flanagan 1999, AAS/AIAA | not held (our `sims_flanagan.py` cites it via Yam 2010) | not listed |
| Sims et al. 2006 (MALTO); Yam & Longuski 2006 (GALLOP) | not held | not listed |
| Fukushima 2012, 2013a, 2013b; Fukushima & Ishizaki 1994 (elliptic-function algorithms) | not held | not listed; needed only if we code the closed form |
| Fenton & Gardiner-Garden 1982 (theta functions) | not held | not listed |
| Nacozy 1976, Celest. Mech. 13:495 | not held | not listed |
| Radhakrishnan & Hindmarsh 1993 (LSODE) | not held | not listed |
| Bate, Mueller & White 1971 | not held | not listed |
| Beletsky 2001; Isayev & Kunitsyn 1972; Kirchgraber 1971; Rufer 1976 | not held | not listed (see Lantoine digest) |
| Sundman 1912, Acta Math. 36:105 | not held | not listed |

Held status checked with `ls cyclers_pdf/papers | grep -i` and `grep -i CORPUS_INDEX.md` for each author; no
hit's file name was the cited work.

*Filed as `cyclers_pdf/papers/hatten-russell-2015-comparison-three-stark-problem-solution-techniques-bounded-case-cmda-121-39-doi-10.1007-s10569-014-9586-z.pdf`.*
