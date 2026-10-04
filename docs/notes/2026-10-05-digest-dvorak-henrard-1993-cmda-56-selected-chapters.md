# Digest: Dvorak & Henrard (eds.) 1993, Qualitative and Quantitative Behaviour of Planetary Systems (CMDA 56, nos. 1-2), selected chapters

Date: 2026-10-05 (Sydney). Reading, plus throwaway numerical checks run outside the repository (section 12) and not
committed. No source file was edited. Tasks touched: `#924` (chaos indicators), `#899`, `#905`, `#884`, `#929`, `#912`,
`#931`, `#933`.

Source: R. Dvorak and J. Henrard (eds.), "Qualitative and Quantitative Behaviour of Planetary Systems", Proceedings of the
Third Alexander von Humboldt Colloquium on Celestial Mechanics (Ramsau, Austria, 29 March to 4 April 1992), Kluwer Academic
Publishers 1993; reprinted from Celestial Mechanics and Dynamical Astronomy volume 56, nos. 1-2 (1993); ISBN
978-94-010-4898-9 (print), 978-94-011-2030-2 (eBook); book DOI 10.1007/978-94-011-2030-2. Filed in the private paper corpus
as
`dvorak-henrard-eds-1993-qualitative-quantitative-behaviour-planetary-systems-humboldt-colloquium-springer-doi-10.1007-978-94-011-2030-2.pdf`
(393 PDF pages, text layer present; printed pages 1 to 396 plus the front matter and a participants photograph). Equations
and table digits quoted below were read from page images at 120 to 300 dpi; where OCR and image disagreed the image was
used and the disagreement is stated.

Evidence tags: READ (p.N) is read at printed page N. COMPUTED is my own calculation (section 12). INFERRED is my
reasoning. "PDF p.N" is the page of the PDF file (it differs from the printed page by an offset that changes through the
volume: +9 near the start, +2 around p.200, 0 at p.307 to 324, -2 at p.373).

## 0. Findings in one place

1. Volume coverage and the Olle-Simo question are in sections 1 and 2. **No chapter by Olle and/or Simo is in this volume**,
   and the volume has 393 printed pages, so a chapter at pp.513-526 cannot be in it.
2. **Hadjidemetriou (pp.201-219)** prints 43 table rows (39 distinct; the e' = 0 rows repeat) of initial conditions of
   symmetric periodic orbits of the planar elliptic restricted problem (Sun-Jupiter, mu = 0.00095387535) at the 2:1 and 4:1
   resonances as a function of the primaries' eccentricity e'. I integrated every printed row and then Newton-corrected
   (x0, ydot0) onto the exact symmetric orbit (perpendicular crossing at t = pi), once the rotating-frame velocity uses the
   instantaneous angular velocity of the Sun-Jupiter line at t = 0 (section 5.3). **Reproduced to the printed digits (a
   correction of at most 2e-5 in x0 and ydot0): 19 of 43 rows** (all of 4:1 I_c, the low-e' rows of 4:1 II_c, I_e and II_e,
   and two high-e' rows); **a further 16 rows have an exact symmetric orbit within 1e-3**, and 7 within 6e-3 (the
   high-eccentricity, very sensitive orbits); **one printed e' is a misprint**: p.211 row "e' = 0.010, x0 = 0.186246"
   has no orbit there and corresponds to e' = 0.100 (correction 5.6e-5 and 6.0e-4, the size of the neighbouring 2:1 rows).
   The printed e and a at e' = 0 agree with the osculating values at t = 0 to the last digit (section 5.3), which
   independently confirms the frame convention. A ready control set for `#912` and `#931` in the planar elliptic problem
   (section 5.4).
3. **Hagel & Trenkler (pp.81-98)**, a Sitnikov-problem chapter outside the requested list, prints Table I (monodromy matrix
   of the linearised Sitnikov problem against e). The project's own `core.er3bp` state-transition matrix reproduces it
   (mu = 0.5, state at the origin; section 8.2): rows at e = +-0.20, +-0.40, +-0.60 agree to 3e-4 absolute and e = +-0.80 to 3e-3 absolute (the paper's RK4 step is
   2 pi/200), apart from the printed slips listed in section 8.2;
   the rows at e = +-0.99 do not (printed values are from a fixed-step RK4 at 2 pi/200 and are not converged). Their
   eq. (26) coefficient 24/121 disagrees with numerics at small e; the second-order coefficient derived independently is
   **21/124** (probable digit transposition).
4. **Valsecchi et al. (pp.373-380)**: the printed Delaunay series for the rates of the lunar node and perigee (eqs. 8 and 9)
   and the printed Saros solution reproduce to the printed digits: (e, i) = (0.0623, 5.57 deg) for the 223:239:3 relation is
   computed as (0.06228, 5.571 deg), and the two shifted-nu cases (0.0871, 5.41 deg; 0.0137, 5.72 deg) are computed as
   (0.08707, 5.412 deg; 0.01367, 5.724 deg). OCR read 0.0871 as "0.OS71"; the page image gives 0.0871.
5. **Froeschle, Froeschle & Lohinger (pp.307-314)** and **Lohinger, Froeschle & Dvorak (pp.315-322)** print almost no
   numbers (the results are figures), but the setup is fully specified and I reproduced the qualitative claims (section 3
   and 4): standard map GLI statistics, and the 3D restricted-problem transition from regular to chaotic at d0 = 0.255.
   **Not the same as FLI**: the generalised indicators are moments of the distribution of the local stretch, a different
   object from the fast Lyapunov indicator already digested.
6. **Yoshida (pp.27-43)** is a review with few numbers; the checkable ones all pass (section 7).
7. **Dvorak (pp.71-80)** prints Table I (classification of 11 eccentricities against 60 starting values of T) and Table II;
   the printed equations (1) to (7) are internally consistent and identical to the project's `core.er3bp` z-equation at
   mu = 0.5; a spot check with a Lyapunov-style indicator agrees with 15 of 16 sampled cells (section 6.3).

## 1. Table of contents of the volume (READ, pp.v-vii)

Preface p.ix. Session titles and chairmen as printed. Pages are printed pages.

| Session | Chapter | Authors | Pages |
| --- | --- | --- | --- |
| Planetary theories (H. Eichhorn) | Review of Planetary and Satellite Theories | P.K. Seidelmann | 1 |
| | The Mystery of Pluto's Mass - The Ring Hypothesis | C. Marchal | 13 |
| | Recent progress in the theory and application of symplectic integrators | H. Yoshida | 27 |
| | Stable Planetary Orbits Around One Component in Nearby Binary Stars (II) | D. Benest | 45 |
| | Stability of Outer Planetary Orbits Around Binary Stars: a Comparison of Hill's and Laplace Stability Criteria | A. Kubala, D. Black, V. Szebehely | 51 |
| | An Excess Motion of the Ascending Node of Mercury in the Observations Used by Le Verrier | T. Inoue | 69 |
| Sitnikov problem (Yi-Sui Sun) | Numerical Results to the Sitnikov-Problem | R. Dvorak | 71 |
| | A Computer Aided Analysis of the Sitnikov Problem | J. Hagel, Th. Trenkler | 81 |
| | The Original Sitnikov Article - New Insights | K. Wodnar | 99 |
| Asteroids (J. Henrard) | Proper Elements: What are They? | A. Lemaitre | 103 |
| | The High Eccentricity Libration of the Hildas II: Synthetic-Theory Approach | T. Michtchenko, S. Ferraz-Mello | 121 |
| | A Stability Study of Asteroid Families Near the 3:1 and 5:2 Resonance with Jupiter | G. Hahn, C.-I. Lagerkvist, B.A. Lindblad | 131 |
| | Chaotic Behaviour of Trajectories for the Asteroidal Resonances | M. Sidlichovsky | 143 |
| | Low-Eccentricity Motion of Asteroids Near the 2/1 Jovian Resonance | J. Schubart | 153 |
| | Numerical Experiments in the 3/1 and 1/6 Overlapping Resonance Region | Ch. Froeschle, H. Scholl | 163 |
| Resonance problems (S. Ferraz-Mello) | An Introduction to Hamiltonian Dynamical Systems and Practical Perturbation Methods: New Insight by Successive Elimination of Perturbation Harmonics | A. Morbidelli | 177 |
| | Frequency Analysis of a Dynamical System | J. Laskar | 191 |
| | An Application of KAM Theory to the Planetary Three Body Problem | Ph. Robutel | 197 |
| | Resonant Motion in the Restricted Three Body Problem | J. Hadjidemetriou | 201 |
| | A Fourth-Order Solution of the Ideal Resonance Problem | B. Erdi, J. Kovacs | 221 |
| General dynamical systems (P.K. Seidelmann) | Wavelet Analysis and Applications to Some Dynamical Systems | Ph. Bendjoya, E. Slezak | 231 |
| | Wavelet Analysis of the Standard Map: Structure and Scaling | A.D. Gilbert, Cl. Froeschle, U. Frisch | 263 |
| | Asteroids: 2/1 Resonance and High Eccentricity | M. Moons, A. Morbidelli | 273 |
| | On the Second Order Long-Period Motion of Hyperion | P.J. Message | 277 |
| | The Problem of Critical Inclination Combined with a Resonance in Mean Motion in Artificial Satellite Theory | F. Delhaise, J. Henrard | 285 |
| Chaos and stability (P.J. Message) | Meteorites from the Asteroid 6 Hebe | P. Farinella, Ch. Froeschle, R. Gonczi | 287 |
| | Generalized Lyapunov-Characteristic Indicators and Corresponding Kolmogorov like Entropy of the Standard Mapping | Cl. Froeschle, Ch. Froeschle, E. Lohinger | 307 |
| | Generalized Lyapunov Exponents Indicators in Hamiltonian Dynamics: An Application to a Double Star System | E. Lohinger, Cl. Froeschle, R. Dvorak | 315 |
| | Asteroid 522 Helga is Chaotic and Stable | A. Milani, A.M. Nobili | 323 |
| Miscellaneous (V. Szebehely) | Classical Periodic Orbits and Quantum Mechanical Eigenvalues and Eigenfunctions | G. Contopoulos | 325 |
| | Generalized Least-Squares Adjustments: A Timely but Much Ignored Tool | H. Eichhorn | 337 |
| | Least Squares Parameter Estimation in Chaotic Differential Equations | J. Kallrath, J.P. Schloder, H.G. Bock | 353 |
| | The arrangement in Mean Elements Space of the Periodic Orbits Close to that of the Moon | G. Valsecchi, E. Perozzi, A.E. Roy, B.A. Steves | 373 |
| | Resonance Trapping of Circumstellar Dust Particles by an Alleged Planet | H. Scholl, F. Roques, B. Sicardy | 381 |
| | Orbital Resonances and Poynting-Robertson Drag Confining Dust Particles in beta-Pic Disk | D. Lazzaro, B. Sicardy, F. Roques, H. Scholl | 395 |

Pagination check (READ, chapter footers): the Scholl-Roques-Sicardy chapter is pp.381-393 and the last chapter (Lazzaro et
al.) is pp.395-396 ("56: 395-396, 1993" on its first page), so the last printed page of the volume is 396 and the PDF's
393 pages run to it (the printed numbering skips p.394 between the last two chapters; the extract does not show it). The
volume carries a group photograph of the participants (PDF p.8), not a list of names, so it cannot be searched (section 2).

Coverage of this digest: read in full (text and the key page images): Froeschle-Froeschle-Lohinger, Lohinger-Froeschle-Dvorak,
Hadjidemetriou, Valsecchi et al., Yoshida, Dvorak, Hagel-Trenkler (sections 2 to 3 and the tables; the first-integral
construction of sections 4 to 5 not rederived), Laskar, Milani-Nobili, Benest. Read at the level of method and tables:
Kallrath-Schloder-Bock. Read at abstract or first-page level only: the other chapters listed in section 11.4. Not opened:
Seidelmann, Inoue, Lazzaro et al.

## 2. Olle and Simo: not in this volume (READ, contents pp.v-vii; COMPUTED text search)

A citation claimed "Olle & Simo, Bifurcations in the RTBP with equal masses, pp.513-526" in this volume. Checks:

- The contents page (section 1) lists no chapter by Olle, Simo, or both, and no chapter on bifurcations in the restricted
  three-body problem with equal masses.
- A case-insensitive text search of the whole extracted text (about 18,800 lines) for "Olle" and "Simo" (whole word)
  returns nothing; the only matches for "Simo" are inside "Simon" (Bretagnon and Simon planetary theories, in Seidelmann's
  reference list).
- The volume ends at printed p.396 (the last chapter is pp.395-396), so pp.513-526 do not exist in it.
- The only equal-mass chapters are the Sitnikov chapters (Dvorak, Hagel and Trenkler, Wodnar), and Benest's mass-ratio
  survey; none is by Olle or Simo.

Conclusion: the claimed attribution is wrong for this volume. The citation should be resolved under `#909` (where the
Olle and Simo 1990 and 1993 papers are tracked), not against this volume.

## 3. Froeschle, Froeschle and Lohinger, "Generalized Lyapunov characteristic indicators and corresponding Kolmogorov like entropy of the standard mapping", CMDA 56:307-314 (PDF pp.307-314)

READ in full. Short (8 pages, figures 1 to 5 only, no tables).

### 3.1 Method and equations (READ pp.308 to 313)

The indicator is the usual Gram-Schmidt Lyapunov machinery (Benettin et al. 1980; p.308), with the **distribution** of the
local stretches kept instead of only its mean. Iterate the map and its tangent map together with a unit initial tangent
vector, renormalising every step (period tau = 1; "j = 1, ... n", p.308); the local stretch is ln(alpha_j), alpha_j the norm
of the evolved vector before renormalisation. Then

    chi(X0, n) = (1/(n tau)) sum_{j=1..n} ln(alpha_j)        (p.308; the first moment)
    m = E(X),  sigma = sqrt(E[(X - m)^2])                    (p.309; X = ln alpha_j)
    gamma_1 = mu_3 / sigma^3,   gamma_2 = mu_4 / sigma^4 - 3,   mu_p = E[(X - m)^p]    (p.309, image-checked)

called the generalised Lyapunov indicators (GLI) when estimated from a finite number n of iterations. The authors call
m, sigma, gamma_1, gamma_2 "generalized Lyapunov exponents" in the limit; gamma_1 and gamma_2 vanish for a normal
distribution, and the authors say gamma_1 = 0 for a symmetric unimodal one (p.309).

The map (p.309, image-checked): **x1 = x0 + a sin(x0 + y0), y1 = x0 + y0 (mod 2 pi)** (the same form as eq. (2) of the
Froeschle-Lega-Gonczi 1997 digest). Pesin (p.312): rho(P) = chi_1(P) for the 2D map; the Kolmogorov entropy is
h = integral over M of rho(P) dmu. The Monte Carlo estimators (p.313) are

    h_k = (1/(N M)) sum_{j=1..M} sum_{i=1..N} ln(alpha_ij)    (k = 1 to 4; no exponent printed)

The printed formula shows h_k as the double sum of ln(alpha_ij) with no visible exponent k (p.313 image); the text says h_1
estimates the Kolmogorov entropy, h_2 "a measure of the dispersion", h_3 and h_4 "the mean over the phase space of the
asymmetry and the flatness parameters". The formula as printed cannot be the four different quantities; INFERRED: the
exponent or a centring was lost. Do not code h_2 to h_4 from this page; use m, sigma, gamma_1, gamma_2 per orbit and
average. M = 50 random orbits and N = 2000 iterations per value of a, a step of -0.2 in a (Fig. 4 caption, p.312).

### 3.2 Printed observations (READ)

- Regular orbit a = -1.3, (x0, y0) = (1, 0) (invariant curve, Fig. 1, p.310): m converges to 0; sigma converges fast to a
  constant; gamma_1 and gamma_2 converge more slowly "to a non-zero value"; the distribution for 1e5 iterations is bimodal
  (Fig. 5a), "neither symmetric unimodal nor normal" (p.309).
- Chaotic orbit a = -1.3, (2, 0) (Fig. 2, p.310): m does not go to zero and "a good approximation to a constant value is not
  reached before 1e4 iterations"; sigma reaches its limit within 1e2 iterations; m and gamma_2 oscillate (cantori); the
  convergence of gamma_1 is "quite fast (~ 10^2 iterations)" (exponent read from the page image), and "the small value
  does not reflect a symmetric distribution as seen on Fig. 5b".
- Transversal scan (Fig. 3, p.311): y0 = 0, x0 from 0 to pi, 10000 iterations. **The step is stated twice and differently**:
  0.02 in the text (p.310) and 0.05 in the Fig. 3 caption (p.311). In the ordered region sigma decreases "to ~ 0.32 when x
  goes to 1" (p.312); sigma is "the smallest and more discriminating variation".
- Fig. 5 (p.313): four distributions of ln(alpha): (a) a = -1.3, (1, 0), invariant curve; (b) a = -1.3, (2, 0), chaotic;
  (c) a = -0.1, (1, 0), "invariant curve in the circulation case"; (d) a = -10, (1, 0), "strong chaotic orbit".
- Fig. 4 (p.312): h_1, h_3, h_4 increase with |a|; h_2 reaches a constant for |a| > 5.
- Conclusion (p.314): m and especially sigma "seem to provide good hints of the global stochasticity"; further work
  "tools of artificial intelligence" for bimodal distributions. The paper itself calls this preliminary.

### 3.3 Reproduction (COMPUTED, section 12.1)

Standard map as printed, tangent vector renormalised each iteration, n up to 1e6:

| Case (a, x0, y0 = 0) | n | m | sigma | gamma_1 | gamma_2 |
| --- | --- | --- | --- | --- | --- |
| -1.3, 1 (regular) | 1e4 / 1e6 | 7.6e-4 / 1.2e-5 | 0.3357 / 0.3355 | 0.509 / 0.512 | -1.232 / -1.230 |
| -1.3, 2 (chaotic) | 1e2 / 1e4 / 1e6 | 0.359 / 0.196 / 0.225 | 0.542 / 0.553 / 0.557 | -0.67 / -0.74 / -0.68 | -0.45 / -0.44 / -0.52 |
| -0.1, 1 (circulation) | 1e6 | 1.4e-5 | 0.0815 | -0.072 | -1.503 |
| -10, 1 (strong chaos) | 1e6 | 1.620 | 0.771 | -1.158 | 1.157 |

All printed qualitative statements reproduce: m to zero for the regular orbit with sigma constant (0.3355; the paper's
"~ 0.32" is a figure reading at x near 1 and is 5 percent from my 0.3355 at exactly x = 1); m settling to about 0.22 for
the chaotic orbit with sigma settled by 1e2 iterations; gamma_2 = -1.23 for the regular orbit (a flat, bimodal
distribution); a strongly skewed distribution at a = -10 whose mean 1.62 agrees with the large-|a| estimate ln(|a|/2) =
1.609 (my external knowledge, not from this volume). The chaotic-sea Lyapunov exponent about 0.22 at |a| = 1.3 agrees
with the Laskar-form map in section 9 (0.2197 at a = +1.3 from (2, 0.5) and 0.2265 at a = -1.3 from (1, 0.3)). The unexplained feature: the paper's "gamma_1
small" for the chaotic orbit (Fig. 2c) is not what I compute (-0.68); the vertical axes of the figure are illegible so I do
not call this a disagreement.

### 3.4 Comparison with the held FLI digest

`docs/notes/2026-10-04-digest-froeschle-lega-gonczi-1997-fli.md` defines psi_1 to psi_3 as inverse norms of the tangent
vectors and uses the time to threshold. The 1993 GLI uses the **same tangent integration** but keeps the distribution of
per-step stretches (m, sigma, gamma_1, gamma_2) and has no threshold-time. It costs nothing extra over FLI, and the extra
moments carry information about cantori (slow convergence, oscillating m and gamma_2) that a single time-to-threshold
loses. The paper offers **no calibrated threshold** and no regular-versus-chaotic decision rule: sigma is "more
discriminating" in a figure only. For `#924` it is therefore a diagnostic to log alongside FLI, not a replacement; the
FLI digest's finding that every model needs a regular and a hyperbolic calibration case applies unchanged.

## 4. Lohinger, Froeschle and Dvorak, "Generalized Lyapunov exponents indicators in Hamiltonian dynamics: an application to a double star system", CMDA 56:315-322 (PDF pp.315-322)

READ in full. Eight pages, figures 1 to 4, no tables.

### 4.1 Model and initial conditions (READ p.316, image-checked)

Three-dimensional circular restricted problem, equal primaries mu = m2/(m1+m2) = 0.5, units as usual (total mass 1,
separation 1, angular velocity 1), rotating frame (xi, eta, zeta), m1 at (-mu, 0, 0), m2 at (1-mu, 0, 0). Equations of motion
as printed (standard; eq. for the third component has no Coriolis term). Satellite-type start from Gonczi and Froeschle
(1981):

    xi = 0.5 - d0,  eta = 0,  zeta = 0
    xi_dot = 0.235749 d0^(-1/2),  eta_dot = -0.640312 d0^(-1/2) + d0,  zeta_dot = xi_dot     (p.316)

d0 = r2(t = 0) is the distance to m2 and the parameter. Sampled values: d0 = 0.10, 0.13, 0.22, 0.24, 0.255, 0.27, 0.31,
0.35, 0.40, 0.50 (p.316; one printed "0040" is OCR for 0.40).

### 4.2 Method (READ pp.317 to 318)

Bulirsch-Stoer integration with a "conventional sequence" n = 2, 4, 6, 8, 12, ...; Gram-Schmidt renormalisation of 3
tangent vectors every time unit for t = 20000 time units (20000 revolutions of the primaries); 42 differential equations
(6 for the orbit plus 6 variational equations for each of the 6 vectors); only the first three local numbers ln(alpha^i)
are examined because chi_i = -chi_(N-i+1). Accuracy parameter epsilon = 1e-9 and 1e-13; "we could only produce the correct
value of the Jacobi constant ... if we chose epsilon = 1e-13 or smaller" (p.318).

### 4.3 Printed results (READ pp.318 to 321)

- Mean (the LCN): near zero for d0 < 0.25, "strong increase, especially of ln alpha^1 and ln alpha^2 ... indicates
  stochasticity for d0 >= 0.255" (p.318). ln alpha^1 and ln alpha^2 are robust to epsilon even in the chaotic zone; ln
  alpha^3 is extremely sensitive to epsilon (p.318 image-checked; OCR had dropped the superscripts).
- sigma is "more or less constant in the chaotic zone, and its value increases with ln alpha^1"; in the regular zone it
  reflects the geometry of the invariant manifold (p.318).
- gamma_1 is "very close to zero, except for ln alpha^1"; gamma_2 behaves like sigma "but of smaller extent"; "in all cases
  we are quite far from a normal distribution" (p.321).
- Distributions (Figs. 2 to 4, x-axis [-6, 8], 100 bins, epsilon = 1e-13): in the regular zone (0.1 <= d0 < 0.25) the
  modal value moves with d0 and the shape changes; in the stochastic zone (d0 > 0.25) the modal value is near zero and the
  shape is nearly the same for every d0, "a uniform chaos" (p.321). No numbers other than the transition (0.25 to 0.255).
- Planned: 1e5 revolutions, other mass ratios, cometary dynamics.

### 4.4 Reproduction (COMPUTED, section 12.2)

My integrator (adaptive Dormand-Prince 5(4) in numba, 24 equations, three tangent vectors, renormalised each time unit,
t = 20000, relative tolerance 1e-11 and 1e-9) gives, for the first vector ln(alpha^1):

| d0 | mean (tol 1e-11) | sigma | mean (tol 1e-9) |
| --- | --- | --- | --- |
| 0.100 | 8.6e-4 | 0.803 | 9.0e-4 |
| 0.130 | 8.1e-4 | 0.711 | 8.8e-4 |
| 0.220 | 7.6e-4 | 0.359 | 8.0e-4 |
| 0.240 | 7.0e-4 | 0.710 | 7.3e-4 |
| **0.255** | **0.113** | 1.303 | 0.120 |
| 0.270 | 0.103 | 1.173 | 0.103 |
| 0.310 | 0.146 | 1.220 | 0.113 |
| 0.350 | 0.141 | 1.195 | 0.148 |
| 0.400 | 0.148 | 1.150 | 0.185 |
| 0.500 | 0.148 | 1.216 | 0.152 |

The transition falls between d0 = 0.24 and 0.255, exactly as printed; regular means are 7e-4 to 9e-4 (of the order of the
finite-time shear floor ln(n)/n, not zero); chaotic means are 0.10 to 0.19 and vary by tens of percent between tolerances
(a chaotic orbit's finite-time mean is not a number to pin); sigma is about 1.15 to 1.33 for every chaotic case ("more or
less constant"); the ln alpha^3 mean differs by a factor 1.5 to 7 between the two tolerances (0.0005 to 0.0013 at 1e-11 against 0.0012 to
0.0033 at 1e-9; d0 = 0.40 is the 7) while the ln alpha^1 mean differs by 0 to 29 percent, the same qualitative sensitivity as
printed. The orbits for d0 >= 0.255 pass
within 0.002 to 0.006 of both primaries (a DOP853 check at t = 2000: r2 min = 0.002 to 0.004, r1 min = 0.003 to 0.006): **a
fixed-step integrator at 5e-4 gave garbage for those cases** (means of 0.002 with kurtosis 1e4), so any reproduction needs
an adaptive or regularised integrator. The Jacobi constant (DOP853, 1e-12, t = 2000) varied by 4e-9 to 6e-8 for d0 = 0.10 to 0.255 and 0.5 and by 5.5e-7 for d0 = 0.31.

### 4.5 Use for `#924`

A published, parameter-controlled regular-to-chaotic transition in the **project's own model class** (3D CR3BP, here
mu = 0.5), with the exact initial conditions printed, so it is a ready calibration for the CR3BP leg of `#924`: regular
cases d0 = 0.10, 0.13, 0.22, 0.24 and chaotic cases d0 = 0.255 to 0.50. The paper gives no threshold; the numbers above
(regular mean <= 1e-3 at t = 2e4, chaotic mean >= 0.1) are my results, not the paper's. The close approaches mean it also
exercises an integrator through repeated primary encounters, a stress test relevant to `#929`.

## 5. Hadjidemetriou, "Resonant motion in the restricted three body problem", CMDA 56:201-219 (PDF pp.203-221)

READ in full; tables on pp.211, 216 and 217 and all equations image-checked. Compare with the held Hadjidemetriou digests
(1975, 1975b, Hadjidemetriou-Christides 1975): those are the general three-body problem; this chapter is the **restricted**
problem at three resonances, circular then elliptic.

### 5.1 Content (READ)

Planar restricted problem, Sun-Jupiter, mu = 0.00095387535 (p.205). Rotating frame xOy with origin at the barycentre and the
x-axis along the Sun-Jupiter line; units: total mass 1, Jupiter's orbit radius 1, G = 1, so Jupiter's period is 2 pi and the
frame's angular velocity n' = 1 in the circular case.

Unperturbed (mu = 0) circular orbits: period T = 2 pi/(n - n') (eq. 1), initial conditions x0 = +-a, y0 = 0, xdot0 = 0,
ydot0 = +-(-a + a^(-1/2)) (eq. 2). Resonant (elliptic) families bifurcate at n : n' = p : q with a = n^(-2/3) = (p/q)^(2/3),
x0 = +-a (1 - e), y0 = 0, ydot0 = +-[-x0 + a^(-1/2) (1 + e)^(1/2) / (1 - e)^(1/2)] (eq. 3); e > 0 pericentre, e < 0
apocentre. Two phase branches (type I, II) per resonance (p.205).

Findings per resonance (READ pp.205 to 209): continuation of the circular family is impossible at 2:1 (the family breaks
down), possible at 3:1 with a small unstable interval AB generated at the resonance, possible at 4:1 with no instability
generated. 2:1: families I (stable) and II (unstable) of the second kind, II stable beyond x0 > 1.0325 after a complicated
structure near a collision orbit (p.206). 3:1: family I unstable, II stable. 4:1: I stable, II unstable. Classification of
all low-order resonances: (i+1)/i behaves like 2:1, (i+2)/i like 3:1, all others like 4:1 (p.218).

Elliptic problem (READ pp.209 to 217): a non-uniformly rotating frame; the problem is a periodic time-dependent two-degree
system, period 2 pi; **symmetric periodic orbits start perpendicularly on the x-axis with Jupiter at perihelion or aphelion**;
the period is 2 pi or a multiple; families are curves in (x0, ydot0, e'), bifurcating from circular-problem orbits whose
period is (p/q) pi (p.210). Stability remark (p.217): in the circular case one pair of multipliers is the unit pair from
the energy integral; when e' is switched on the pair is free to leave the unit circle, which is how instability appears;
every elliptic family that starts at non-zero e does so from the stable branch of the circular problem.

Averaged Hamiltonian at 3:1 (eq. 5, p.213, image-checked): H = H0(N, S) + mu H1(N, S, sigma) + mu e' H2(S, sigma, N, nu) with
H0 = -2(1-mu)^2/(N-S)^2, H1 = 2 F S - b (S/N) cos 2 sigma, H2 = sqrt(2S)[G cos(sigma+nu) + D cos(sigma-nu)] + 2 mu e' K cos 2 nu;
S = sqrt((1-mu) a)(1 - sqrt(1-e^2)), sigma = (3 lambda' - lambda)/2 - omega, N = sqrt((1-mu)a)(3 - sqrt(1-e^2)),
nu = -(3 lambda' - lambda)/2 + omega'; coefficients **b = 2.392398, F = -0.205070, G = 0.198705, D = 2.656407, K = -0.181477,
e' = 0.048**. The paper's correction term that supplies the missing high-eccentricity fixed points (eq. 6, p.214):
**H_c = -3.4 mu (1 - 6.4 e' cos(sigma + nu)) S^3**. These are as printed; I did not test them.

### 5.2 Printed tables (READ, image-checked at 130 dpi)

Isolated bifurcation orbits: p.210 2:1 orbit A at x = 0.1657759, ydot0 = 3.057715, e = 0.72 (the table gives x0 = 0.165776
and e = 0.735 for the same point: the two eccentricities are not the same number, INFERRED: 0.72 is the
Morbidelli-Giorgilli-style mean value and 0.735 the osculating value at t = 0); p.212 3:1 orbits (x0, ydot0) = (0.479420,
0.962393) period pi, and (-0.866333, 0.384417), e = 0.798, period 2 pi; p.215 4:1 orbits (0.395733, 1.190494) period
2 pi/3 (taken three times), (0.299295, 1.733614) e = 0.243, and (0.077836, 4.700496) e = 0.801 ("4.700496" here against
"4.700478" in the p.217 tables: one of them is a slip).

Family tables (x0, ydot0, e, a; columns as printed, e' is the primaries' eccentricity; e < 0 means apocentre):

2:1 resonance (p.211), "Jupiter at perihelion" family I_e and "Jupiter at aphelion" family II_e:

| family | e' | x0 | ydot0 | e | a |
| --- | --- | --- | --- | --- | --- |
| I_e | 0 | 0.165776 | 3.057715 | 0.735 | 0.6295 |
| I_e | 0.020 | 0.172130 | 2.975535 | 0.725 | 0.6301 |
| I_e | 0.048 | 0.187429 | 2.796163 | 0.701 | 0.6308 |
| I_e | "0.010" | 0.186246 | 2.785537 | 0.704 | 0.6322 |
| II_e | 0 | 0.165776 | 3.057715 | 0.735 | 0.6295 |
| II_e | 0.020 | 0.136000 | 3.473983 | 0.782 | 0.6286 |
| II_e | 0.048 | 0.109003 | 3.971428 | 0.824 | 0.6263 |
| II_e | 0.075 | 0.077703 | 4.809806 | 0.873 | 0.6202 |

4:1 resonance (pp.216 to 217), families I_c (peri), II_c (aph) from the first-kind orbit; I_e, II_e from the orbit at
e = 0.243; III_e (peri), IV_e (aph) from the orbit at e = 0.801:

| family | e' | x0 | ydot0 | e | a |
| --- | --- | --- | --- | --- | --- |
| I_c | 0 | 0.395733 | 1.190494 | 0.0003 | 0.3972 |
| I_c | 0.05 | 0.378593 | 1.237525 | 0.0447 | 0.3973 |
| I_c | 0.10 | 0.354497 | 1.325581 | 0.1057 | 0.3975 |
| I_c | 0.20 | 0.300954 | 1.564030 | 0.2414 | 0.3980 |
| I_c | 0.35 | 0.225103 | 2.016289 | 0.4346 | 0.3998 |
| II_c | 0.01 | 0.397015 | 1.191898 | -0.0021 | 0.3937 |
| II_c | 0.02 | 0.396386 | 1.202629 | -0.0006 | 0.3941 |
| II_c | 0.03 | 0.390295 | 1.240289 | 0.0147 | 0.3946 |
| II_c | 0.05 | 0.360970 | 1.404728 | 0.0875 | 0.3966 |
| II_c | 0.10 | 0.292971 | 1.826506 | 0.2585 | 0.3964 |
| II_c | 0.16 | 0.228092 | 2.321809 | 0.4213 | 0.3959 |
| I_e | 0 | 0.299295 | 1.733614 | 0.243 | 0.3967 |
| I_e | 0.048 | 0.217642 | 2.332853 | 0.450 | 0.3974 |
| I_e | 0.080 | 0.146939 | 3.142349 | 0.628 | 0.3979 |
| I_e | 0.090 | 0.106069 | 3.892089 | 0.732 | 0.3992 |
| I_e | 0.100 | 0.052047 | 5.872047 | 0.870 | 0.4079 |
| II_e | 0.01 | 0.317271 | 1.627262 | 0.1986 | 0.3971 |
| II_e | 0.02 | 0.336945 | 1.517788 | 0.1490 | 0.3971 |
| II_e | 0.03 | 0.361906 | 1.386701 | 0.0452 | 0.3967 |
| II_e | 0.04 | 0.385989 | 1.268952 | 0.0245 | 0.3967 |
| II_e | 0.05 | 0.410069 | 1.158271 | -0.0360 | 0.3966 |
| III_e | 0 | 0.077836 | 4.700478 | 0.801 | 0.3967 |
| III_e | 0.05 | 0.080606 | 4.599587 | 0.796 | 0.3990 |
| III_e | 0.08 | 0.078027 | 4.682457 | 0.803 | 0.4007 |
| III_e | 0.10 | 0.071419 | 4.923792 | 0.820 | 0.4027 |
| III_e | 0.12 | 0.06056432 | 5.40045369 | 0.849 | 0.4067 |
| III_e | 0.15 | 0.03636179 | 7.10683320 | 0.914 | 0.4325 |
| IV_e | 0.010 | 0.083159 | 4.525573 | 0.788 | 0.3963 |
| IV_e | 0.020 | 0.213797 | 2.398451 | 0.459 | 0.3966 |
| IV_e | 0.030 | 0.239965 | 2.176001 | 0.392 | 0.3966 |
| IV_e | 0.035 | 0.170342 | 2.864093 | 0.568 | 0.3964 |
| IV_e | 0.041 | 0.023039 | 8.956967 | 0.936 | 0.3764 |

(Also printed at e' = 0: II_c = I_c row, II_e = I_e row, IV_e = III_e row; omitted above as duplicates. Two entries of the
p.216 table, II_c e' = 0.01 and e' = 0.02, have eccentricities -0.0021 and -0.0006, and I_c e' = 0 has e = 0.0003; these are
as printed and the negative sign means apocentre per eq. 3.)

Printed stability labels (READ): 2:1 I_e stable, II_e unstable (pp.206, 211: the text on p.212 gives e = 0.71 stable and
0.76 unstable for the two e' = 0.048 orbits P1 and P2, taken from Morbidelli and Giorgilli; "the numerical computations we
made for the stability are not very accurate"); 3:1 I_e stable and the other three unstable (p.212); 4:1 I_c and II_c unstable,
I_e and II_e stable, stability of III_e and IV_e "could not be obtained with good accuracy" (p.215). The p.218 summary says
the 4:1 labels as above.

Flag, p.217: the IV_e family is **not smooth in e'**: x0 = 0.083159 at e' = 0.010 jumps to 0.213797 at 0.020, then 0.239965,
0.170342, 0.023039; the printed e = 0.788, 0.459, 0.392, 0.568, 0.936 are likewise non-monotonic. My closure test (below)
shows that the e' = 0.020, 0.030 and 0.035 rows close to 1.4e-4, 1.0e-3 and 5.5e-3, while e' = 0.010 and 0.041 close only to 4.7e-2 and 1.5e-2 (high-eccentricity, sensitive orbits).

### 5.3 Reproduction of the tables (COMPUTED, section 12.3)

Method: planar problem, the Sun at -mu r and Jupiter at (1 - mu) r from the barycentre with r(t) the Kepler distance
(a = 1, GM = 1, Jupiter at perihelion or aphelion at t = 0 on the +x axis); asteroid started at (x0, 0) with the inertial
velocity (0, ydot0 + omega x0), **where omega is the instantaneous angular velocity of the Sun-Jupiter line at t = 0**
(sqrt(1+e')/(1-e')^(3/2) at perihelion, sqrt(1-e')/(1+e')^(3/2) at aphelion). Integrate DOP853 (1e-12) over 2 pi and
compare the rotating-frame state at t = 2 pi with the start. The paper does not say what omega it uses for ydot0; with
omega = 1 (the mean angular velocity) the elliptic rows do **not** close (residuals 0.3 to 1.2), so the printed ydot0 is
measured in a frame rotating at the instantaneous rate. The circular row (e' = 0) is unaffected by the choice.

Residual |(dx, y, dxdot, dydot)| at t = 2 pi, by family (all rows listed in section 5.2 were run):

- 4:1 I_c: 1.2e-5, 2.6e-5, 1.8e-4, 2.8e-4, 1.0e-3 (e' = 0 to 0.35).
- 4:1 II_c: 1.2e-5, 1.5e-4, 3.3e-5, 2.4e-6, 5.7e-4, 2.0e-3, 2.2e-3 (e' = 0 to 0.16).
- 4:1 I_e: 1.4e-4, 6.6e-4, 2.6e-3, 9.9e-4, 6.9e-3. 4:1 II_e: 1.4e-4, 1.3e-4, 2.5e-4, 8.3e-5, 9.6e-5, 1.2e-3.
- 4:1 III_e: 9.3e-4, 5.7e-2, 3.8e-2, 1.6e-2, 2.3e-3, 1.7e-2. 4:1 IV_e: 9.3e-4, 4.7e-2, 1.4e-4, 1.0e-3, 5.5e-3, 1.5e-2.
- 2:1 I_e: 9.9e-4, 8.8e-4, 2.9e-3 and the e' = "0.010" row 2.7 (fails). 2:1 II_e: 9.9e-4, 7.6e-3, 3.8e-3, 2.9e-2.

**A small closure residual does not by itself show that a row is periodic** (the orbits are unstable and a 6-digit initial
condition is amplified), so I also solved for the nearest exact symmetric orbit: Newton iteration (SciPy fsolve, 1e-11
integration) on the two half-period symmetry conditions y(pi) = 0 and xdot(pi) = 0 for (x0, ydot0) at fixed e', started at the
printed row, with the maximum of |Delta x0|, |Delta ydot0| as the measure of how far the printed row is from an exact
orbit. Result for the 43 rows (excluding the e' = "0.010" row, counted separately):

| max(|Delta x0|, |Delta ydot0|) | rows | which |
| --- | --- | --- |
| <= 5e-6 | 15 | 4:1 I_c e' = 0, 0.05, 0.10, 0.20, 0.35 (4e-7 or less); 4:1 II_c e' = 0, 0.01, 0.02, 0.03; 4:1 I_e e' = 0.048, 0.080; 4:1 II_e e' = 0.01, 0.02, 0.03; 2:1 II_e e' = 0.075 (3e-6) |
| 5e-6 to 2e-5 | 4 | 4:1 I_e and II_e at e' = 0 (5.1e-6), 4:1 I_e e' = 0.090 (1.5e-5), 4:1 III_e e' = 0.15 (9.5e-6) |
| 2e-5 to 1e-3 | 16 | the 2:1 rows at e' = 0 (3.9e-4), I_e 0.048 (7.9e-4), II_e 0.02, 0.048 (9e-5, 3.7e-4); 4:1 II_c 0.16, I_e 0.100, II_e 0.04, III_e 0 and 0.10 and 0.12, IV_e 0, 0.010, 0.020, 0.035, 0.041 |
| 1e-3 to 6e-3 | 7 | 2:1 I_e e' = 0.02 (1.8e-3); 4:1 II_c e' = 0.05 (5.0e-3) and 0.10 (2.1e-3); II_e e' = 0.05 (6.1e-3); III_e e' = 0.05 (4.1e-3) and 0.08 (1.4e-3); IV_e e' = 0.03 (1.4e-3) |

So 19 of the 43 rows reproduce to the printed digits (<= 2e-5), and every other row has an exact symmetric periodic orbit within
6e-3 in (x0, ydot0); the rows with the large corrections are the high-eccentricity, strongly unstable orbits that the paper
itself says "are very sensitive to the initial conditions" (p.212) and whose stability it could not compute. At e' = 0
the 2:1 pair differs from the exact orbit by 3.2e-5 in x0 and 3.9e-4 in ydot0 (exact 0.165808, 3.057323 against the printed
0.165776, 3.057715), which is a property of the paper's table, not of my frame convention (the circular problem has no
convention to choose). **The e' = "0.010" row of 2:1 I_e (p.211)** has no orbit anywhere near (Newton diverges to
(-0.054, 5.22)); scanning e' at the printed (x0, ydot0), the closure residual is small only at e' = 0.100 (8.6e-4; none of
the other e' tried, 0.040 to 0.125 in steps of 0.005), and a Newton correction there is (5.6e-5, -6.0e-4), the same size as
the neighbouring 2:1 rows. The printed e' = 0.010 is therefore a misprint for 0.100: the corrected row is e' = 0.100, x0 =
0.186246, ydot0 = 2.785537 (the family then reads x0 = 0.165776, 0.172130, 0.187429, 0.186246 at e' = 0, 0.02, 0.048, 0.100,
not monotone in x0, with e' increasing in the table's order).

**The two printed values of the 4:1 III_e circular-problem orbit are not settled**: for x0 = 0.077836 the closure
residual at t = 2 pi is 9.3e-4 for ydot0 = 4.700478 and a Newton solve from either value lands on a different nearby orbit
(4.700535 from 4.700478; 4.700174 from 4.700496) with residuals of 5e-12 and 3.5e-9, i.e. a flat valley in which the 6-digit
difference of 1.8e-5 cannot be resolved. Both are consistent with the paper; use 4.700478 (the family tables) and record the
text value 4.700496 (p.215) as a variant.

**Frame convention confirmed independently of closure (COMPUTED).** With the heliocentric velocity of the asteroid (the Sun
moves at -mu times Jupiter's velocity, so the barycentric velocity ydot0 + omega x0 is corrected by +mu omega r_J), GM = 1 - mu
and r = x0 + mu r_J, the osculating elements at t = 0 are: 2:1 I_e at e' = 0: a = 0.62953, e = 0.73515 (printed 0.6295,
0.735); 4:1 I_e e' = 0: 0.39673, 0.24319 (printed 0.3967, 0.243); 4:1 III_e e' = 0: 0.39671, 0.80139 (printed 0.3967, 0.801);
4:1 I_c e' = 0: a = 0.39679, e = 0.00027 (printed 0.3972, 0.0003, so a differs in the third digit here). So **at e' = 0
the printed a and e are the osculating elements at t = 0**, to the last printed digit. For e' > 0, with omega = the
instantaneous angular velocity and the Sun's velocity -mu r_J omega, the osculating eccentricity at t = 0 is 0.7251 (2:1
I_e e' = 0.020; printed 0.725), 0.7009 (e' = 0.048; printed 0.701), 0.7824 (2:1 II_e e' = 0.020; printed 0.782), 0.4491
(4:1 I_e e' = 0.048; printed 0.450), 0.0436 (4:1 I_c e' = 0.050; printed 0.0447), 0.0877 (4:1 II_c e' = 0.050; printed
0.0875): the printed eccentricity matches the osculating value at t = 0 to 1e-3 or better (0.0011 for the 0.0447 row). The
printed a for e' > 0 is not the osculating value (it varies with e' as 0.6295, 0.6301, 0.6308 while the osculating a at
t = 0 stays at 0.6295 to 0.6296), so it is some other average; INFERRED, the paper does not say.

### 5.4 Techniques applicable to the project's problems

- **`#912` (elliptic-problem periodic-orbit code lacks what the published method needs).** This chapter is the planar
  elliptic benchmark that `#912` lacks: 43 printed rows (39 distinct), six digits, both apse phases, e' to 0.35, 2:1 and 4:1,
  two high-eccentricity branches. They are in a **non-pulsating, instantaneously rotating** frame with physical (not
  pulsating) distances; the project's `core/er3bp.py` is pulsating, so the conversion is: distance unit Jupiter's semi-major
  axis, state at the apse; pulsating position = physical position / r(t0) with r(t0) = 1 -+ e', and the velocity
  conversion needs r(t0), f-dot(t0) and e'. Use my section 5.3 construction as the reference and a pulsating-frame run as
  the system under test. Test: closure after 2 pi, perpendicular crossing at t = pi (the symmetry), and the printed
  stability labels (monodromy of the 2 pi map).
- **`#931` (Hadjidemetriou follow-ups).** The stability statement on p.217 is a ready test for the classifier: at e' = 0
  the monodromy of a circular-problem orbit has a unit pair (from the energy integral) and one more pair; for e' > 0 the
  unit pair splits, and a family that starts stable at e' = 0 can lose stability through that pair. A correct classifier on
  the elliptic monodromy must therefore never use "a unit pair is present" as a test once e' > 0.
- **`#905`/`#884` (Earth-Moon cyclers in the Sun-perturbed model).** The Sun's eccentric forcing is the same structure
  (a 2 pi-periodic time-dependent perturbation of a circular-problem family); the chapter's principle that **continuation
  of a circular-problem periodic orbit into the elliptic problem is possible only at points whose period is (p/q) pi** is the
  elliptic analogue of the commensurability condition used in `#905`, and its statement that such families bifurcate only
  from the stable branch of the second kind restricts which Earth-Moon members to try first. INFERRED mapping.
- **`#933` (Broucke 1969).** The chapter cites Broucke (1968, 1969) for the theory (p.210) and its periodic orbits are in the
  same symmetry class (start perpendicular to the x-axis with the primary at an apse; period 2 pi or a multiple). The two
  sources use different frames and mass ratios (Broucke's tables are in the pulsating frame, this chapter's mu = 0.00095
  rows in a physical rotating frame), so no printed row of one is an input to the other; the family structure (which
  circular-problem points are bifurcation points, p.210) is the shared content. INFERRED, not compared.
- **`#899`.** The 4:1 and 2:1 first-kind and second-kind families of the circular problem near mu = 0.00095 (a mass ratio
  far smaller than Earth-Moon) are not a one-moon cycler test; use the book only for the tabulated families' continuation
  structure (which branches exist, which are stable).
- The averaged-model remark (p.203 and pp.213 to 215): an averaged Hamiltonian that omits the high-eccentricity resonant
  fixed points gives qualitatively wrong long-term evolution (Fig. 10); the authors fix this by adding a correction term
  found from the periodic orbits (eq. 6). Not used in this project (no averaged model).

## 6. Dvorak, "Numerical results to the Sitnikov-problem", CMDA 56:71-80 (PDF pp.78-87)

READ in full; Table I and II image-checked at 300 dpi.

### 6.1 Model (READ pp.72 to 73, image-checked)

Sitnikov problem: massless body on the axis through the barycentre perpendicular to the plane of two equal primaries on
Keplerian ellipses of eccentricity e. Equations (units: primaries of mass 1/2 each, relative orbit semi-major axis 1,
period 2 pi; INFERRED from the printed forms and confirmed by the computation below):

    z'' + z/(r^2 + z^2)^(3/2) = 0          (1)
    r(v) = (1 - e^2) / (2 (1 + e cos v))   (2)   (distance of one primary to the barycentre, v true anomaly)
    T := z/(2r) = z (1 + e cos v)/(1 - e^2)                                  (3)
    T' = dT/dv = z' sqrt(1-e^2)/(1 + e cos v) - z e sin v/(1 - e^2)          (4)
    z = T (1 - e^2)/(1 + e cos v)                                             (5)
    z' = (1/sqrt(1-e^2)) [T e sin v + T' (1 + e cos v)]                      (6)
    T'' + [ (1/sqrt(1/4 + T^2))^3 + e cos v ] / (1 + e cos v) * T = 0         (7)

(Wodnar's form; T is half the tangent of the angle subtended by the body from one primary.) The Poincare section is T
against T' at every pericentre passage of the primaries (v = 2 pi k); integration by a variable-step Lie-series method
(Hanslmeier and Dvorak 1984) over 1000 primary revolutions, 5000 where a fractal island structure appeared; initial
conditions v = 0, T' = 0, T varied, e between 0.33 and 0.66 (p.73).

COMPUTED (section 12.4): integrating eq. (1) with Kepler r(v) and integrating eq. (7) from the same start (z0 = T0 (1-e) by
eq. (5), z-dot = 0 by eq. (6)) give the same z(v) to 10 digits at v = pi, 2 pi, 6 pi for (e, T0) = (0.33, 1.0) and (0.66,
0.5). **The printed equations (1) to (7) are mutually consistent.** The project's `core.er3bp.er3bp_eom` z-equation,
`zdoubleprime = -scale * (e cos f * z + grav_z)` with scale = 1/(1 + e cos f), at x = y = 0 and mu = 0.5 (so r1 = r2 =
sqrt(1/4 + z^2)), is **identically eq. (7)** with T = z (pulsating z, unit primary separation): so every Sitnikov number in
this volume is a direct test of `core.er3bp` and of the Floquet and Poincare-section machinery built on it.

### 6.2 Printed tables (READ pp.77 to 78; image-checked)

Legend of Table I (p.77): "0" orbit on an invariant curve; "a" invariant curve with some accumulation of points, close to a
separatrix; "n" (a digit) number of islands of the orbit; "s" on or very close to a separatrix; "k" chaotic character; "*"
10 or more islands. Each cell block gives T from the header value in steps of 0.01 (ten characters per block); the
printed glyph for "0" is a lowercase "o"; I normalise it to 0. The text says the T interval is 0.4 to 1.1 with step 0.01
(p.77) but the table shows only 0.40 to 0.99 (six blocks of ten).

| e | T = 0.40 to 0.49 | 0.50 to 0.59 | 0.60 to 0.69 | 0.70 to 0.79 | 0.80 to 0.89 | 0.90 to 0.99 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.33 | 0000000000 | 1000000a00 | 0a00222220 | 0a00111110 | 0000000000 | 044000669k |
| 0.37 | 0000000000 | 0010000000 | 0000000000 | 0000a22222 | 2000000000 | 0044?5?0kk |
| 0.41 | 0000000000 | 0001000000 | 0000000000 | 0000002222 | 2220000000 | 00440006kk |
| 0.44 | 0000000000 | 0000100000 | 0000004000 | 0a00002222 | 2222200000 | 300444s6kk |
| 0.47 | 0000000000 | 0000010000 | 0000000400 | 000000a222 | 2222220000 | 030a44skkk |
| 0.51 | 0000000000 | 0000001100 | 0000000000 | 0000000022 | 22222222s0 | 003*k44kkk |
| 0.55 | 0000000000 | 0000a00111 | 0000000000 | 000a0000s* | 2222222222 | 2d8k3kk4kk |
| 0.58 | 0000000000 | 0000000001 | 1000000000 | 000000000s | 2222222222 | 2sssk3kkkk |
| 0.61 | 0000a00000 | a000a00000 | 1110000006 | 0004000000 | k2*2222222 | 222kkkkkkk |
| 0.64 | 0000000000 | 0000000000 | 0111100000 | 0000440030 | kkkkk22222 | 22222kkkkk |
| 0.66 | 0000000000 | 0000000000 | 0111110000 | 0000044403 | kkkkkk2222 | 22220*kkkk |

Flagged symbols: in e = 0.37, block 0.90 the printed string is "oo44?5?okk" with two question
marks (illegible or undefined in the legend); in e = 0.55 block 0.90 the string "2d8k3kk4kk" contains "d", which is not
in the legend; both are as printed.

Table II, "The last invariant curve" (p.78, image-checked):

| e | 0.33 | 0.37 | 0.41 | 0.44 | 0.47 | 0.51 | 0.55 | 0.61 | 0.64 | 0.66 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T | 0.99 | 0.98 | 0.98 | 0.97 | 0.92 | 0.91 | 0.92 | 0.80 | 0.80 | 0.80 |

**Table II is not consistent with Table I as a "last invariant curve"**: at e = 0.33 Table I has an island glyph "9" at
T = 0.98 and "k" (chaotic) at 0.99, so 0.99 is the first chaotic T rather than the last invariant curve; at e = 0.61, 0.64
and 0.66 Table II gives 0.80, which is the first "k" of the 0.8 block in Table I, yet Table I shows island blocks ("2")
and more regular cells at larger T. The text (p.78) calls the row "the last island ... from here on global chaos can
arise", which differs from the caption. Table II has no e = 0.58 entry although Table I has the row. Do not use Table II as
a threshold without recomputing it.

Other printed observations (READ pp.73 to 79): (i) the one-to-three period splitting of a stable periodic orbit with e at
T0 = 1.13 (Fig. 1, six values of e, not printed); (ii) at T0 = 1.05, e = 0.23 shows an invariant curve, e = 0.28 a chain of
7 islands, e = 0.33 17 islands (Fig. 3, p.73, p.75); (iii) invariant curves reach larger T for smaller e; the onset of global
chaos moves to smaller T as e grows; (iv) the first island moves outward with e (p.78); (v) Fig. 6 is described in the text
(p.78 to 79) as "e = 0.33" and captioned e = 0.37 (p.77): a mismatch in the source; (vi) the complexity is "already
present for very small eccentricities, e ~ 0.0001" (Liu and Sun 1991; p.71). The conclusion calls the problem "generic"
(p.78) and says Lyapunov exponents were not computed because they are "expensive" (p.72), "kept for a future project".

### 6.3 Spot check of Table I (COMPUTED, section 12.4)

A single-orbit finite-time Lyapunov exponent per primaries' revolution (eq. (7) with its variational equation, 300
revolutions, renormalised each revolution; regular values 0.011 to 0.030, chaotic 0.07 to 0.27 on this scale; the thresholds
are mine) against the printed cell at the same (e, T):

| e | T0 | my exponent | printed cell | agree |
| --- | --- | --- | --- | --- |
| 0.33 | 0.45, 0.65, 0.85, 0.95 | 0.027 each | 0, 2 (islands), 0, 0 | yes |
| 0.33 | 0.99 | 0.273 | k | yes |
| 0.51 | 0.55 | 0.027 | 0 | yes |
| 0.51 | 0.85 | 0.019 | 2 | yes |
| 0.51 | 0.94 | 0.151 | k | yes |
| 0.51 | 0.97 | 0.072 | k | yes (weak) |
| 0.66 | 0.45 | 0.026 | 0 | yes |
| 0.66 | 0.62 | 0.021 | 1 (island chain) | yes |
| 0.66 | 0.78 | 0.028 | 0 | yes |
| 0.66 | 0.80 | 0.030 | k | **no** (inconclusive: a thin layer; 300 revolutions) |
| 0.66 | 0.82 | 0.107 | k | yes |
| 0.66 | 0.86 | 0.019 | 2 | yes |
| 0.66 | 0.95 | 0.011 | * (10 or more islands) | yes (regular) |

15 of 16 sampled cells agree; the one that does not is the first "k" of a block and a single short orbit. This is a sampled
spot check, not a reproduction of the table; a full reproduction needs a 60 by 11 grid at 1000 revolutions and a
classification rule for islands and separatrix layers, which the paper does not give.

### 6.4 Techniques applicable

- A **sourced elliptic restricted control with a known chaotic boundary** in the project's own equation set
  (`core.er3bp` at mu = 0.5, x = y = 0): positives (e = 0.33, T0 = 0.99; e = 0.51, T0 = 0.94; e = 0.66, T0 = 0.82) and
  negatives (every "0" cell), for the `#924` calibration on the elliptic model, with the caveat that Table I is graphical
  classification by eye, no thresholds, and some glyphs are unclear (above).
- For `#912`/`#933` it is **not** a periodic-orbit control (the z-axis orbits of the Sitnikov problem are not in the
  Earth-Moon plane). The nearest periodic content is Hagel-Trenkler's linear monodromy (section 8) and the period-1 to
  period-3 splitting statement above.

## 7. Yoshida, "Recent progress in the theory and application of symplectic integrators", CMDA 56:27-43 (PDF pp.36-52)

READ in full; equations (7), (8), (28), (37), (40), (41), (43) to (46), (61) to (63) image-checked. A review; it prints no
tables and only qualitative integration results.

### 7.1 Content (READ)

Motivation: Euler multiplies the oscillator energy by (1 + tau^2) per step (eq. 7), classical RK4 damps it by
(1 - tau^6/72 + ...) (eq. 8). For autonomous Hamiltonian systems a scheme cannot conserve both energy and symplectic
structure except by reparametrising the exact flow (Ge and Marsden 1988; p.28).

Implicit schemes (pp.29 to 31): generating-function methods (Feng and Qin; Channell and Scovel up to order 6); implicit
Runge-Kutta schemes with the Sanz-Serna/Lasagni condition (18) (conserve quadratic integrals); implicit midpoint (eq. 21,
order 2) and the 2-stage Gauss-Legendre scheme (eq. 22, order 4); the s-stage Gauss-Legendre method has order 2 s.

Explicit schemes for H = T(p) + V(q) (pp.31 to 34): composition S_T(c tau) S_V(d tau) with coefficients solving the order
conditions. Printed: order 2 leap-frog c = (1/2, 1/2), d = (1, 0); **order 3 (Ruth 1983, eq. 28): c = (7/24, 3/4, -1/24),
d = (2/3, -2/3, 1)**; **order 4 (Forest-Ruth, eq. 37): c1 = c4 = 1/(2(2 - 2^(1/3))), c2 = c3 = (1 - 2^(1/3))/(2(2 - 2^(1/3))),
d1 = d3 = 1/(2 - 2^(1/3)), d2 = -2^(1/3)/(2 - 2^(1/3)), d4 = 0**, needing 3 force evaluations per step against 4 for RK4;
Yoshida's composition S_4(tau) = S_2(x1 tau) S_2(x0 tau) S_2(x1 tau) with **x0 = -2^(1/3)/(2 - 2^(1/3)), x1 = 1/(2 -
2^(1/3))** from x0 + 2 x1 = 1, x0^3 + 2 x1^3 = 0 (eqs. 38 to 41); the order-6 scheme is a composition with three sets of
(w0, w1, w2, w3) "obtained numerically" and five for order 8 (eq. 42, values not printed); Suzuki: no solution with all
positive coefficients for n >= 3.

Energy behaviour (pp.34 to 36): the symplectic Euler map for the oscillator is q' = q + tau p, p' = p - tau q' (eq. 43) with
the conserved quantity (p^2 + q^2)/2 + (tau/2) p q (eq. 44), so the orbit through (1, 0) lies on q^2 + p^2 + tau p q = 1;
Theorem 1 (pp.35 to 36): the symplectic map exactly describes the evolution under a formal "modified Hamiltonian" H~ = H + tau
H_1 + tau^2 H_2 + ..., H_1 = (1/2) H_p H_q, H_2 = (1/12)(H_pp H_q^2 + H_qq H_p^2), H_3 = (1/12) H_pp H_qq H_p H_q (eq. 46),
by the Baker-Campbell-Hausdorff formula (eqs. 47 to 53); for an order-n scheme the energy error stays of order tau^n; the
series is not guaranteed to converge for nonlinear systems.

Applications (pp.37 to 38): Kinoshita, Yoshida and Nakai on Kepler e = 0.1: RK4 better over 10 periods, SI4 better after 2000
periods (RK4 errors in a and e grow secularly; the mean anomaly error grows quadratically under RK4, linearly under SI4;
the argument-of-pericentre error grows linearly in both and is larger for SI4); Wisdom and Holman 1991, second-order scheme,
tau = 1 year, to 1e9 years (printed "109 years", OCR; INFERRED 10^9); Gladman and Duncan 1990, fourth order.

Variable step (pp.38 to 39): with the step set by a dt = r ds (eq. 59) a symplectic integrator loses its bounded energy error
(Fig. 1: Kepler a = 1, e = 0.5, constant tau = 0.05, 10000 iterations = 80 periods; the variable-step RK4 error drops by one
order while SI4 develops a secular growth), "no advantage to use a variable step when integrating by known symplectic
methods". Large step (pp.39 to 41): the dual first-order map p' = p + tau sin q, q' = q + tau p' (eq. 61) is the standard
map p' = P + k sin Q, Q' = Q + P' with tau p = P, q = Q, tau^2 = k (eqs. 61 and 62); its modified Hamiltonian to third order
is H~ = p^2/2 + cos q + (tau/2) p sin q + (tau^2/12)(sin^2 q - p^2 cos q) - (tau^3/12) p sin q cos q + o(tau^4) (eq. 63); for
tau = 0.2 the contours of H~ match the iterates, at 0.7 a visible discrepancy, at 1.0 the series is not the conserved
quantity except near the elliptic fixed point (q, p) = (pi, 0) and most initial conditions are chaotic.

### 7.2 Checks (COMPUTED, section 12.5)

- Eq. (8): RK4 on the oscillator multiplies p^2 + q^2 per step by 0.99999998613 at tau = 0.1, 0.99999911556 at 0.2,
  0.99978977 at 0.5, against 1 - tau^6/72 = 0.99999998611, 0.99999911111, 0.99978299: the printed leading term is right.
- Eq. (63): I iterated eq. (61) for 20000 steps from (q, p) = (1.0, 0): the spread of the printed H~ is 1.19e-4 at tau = 0.2 and
  4.6e-7 at tau = 0.05, a ratio 259 against tau^4 scaling 256, while the unmodified H spreads by 0.25 and 0.063 (ratio 4, O(tau)).
  So the signs and coefficients of eq. (63) are right (including the sign flip of the odd terms against eq. 46 for the dual
  ordering).
- Eq. (37) numerical values (COMPUTED from the printed closed form): x1 = 1.3512071920, x0 = -1.7024143839; c1 = c4 =
  0.6756035960, c2 = c3 = -0.1756035960, d1 = d3 = 1.3512071920, d2 = -1.7024143839.
- The Fig. 2 values tau = 0.2, 0.7, 1.0 correspond to k = tau^2 = 0.04, 0.49, 1.00 (eq. 62); the critical parameter of the
  standard map is k_c = 0.9716 (Laskar's chapter, section 9), so the statement that at tau = 1.0 "most initial conditions"
  are chaotic is consistent with k = 1.0 just above k_c (INFERRED, cross-chapter).

### 7.3 Techniques applicable (for `#929` and `#924`)

- **`#929` (integrator controls through one moon flyby).** Reversibility and energy diagnostics for a symmetric
  integrator are standard; the chapter adds a concrete test and two warnings. Test: iterate the dual map (61) with the
  modified Hamiltonian (63) as a conserved quantity (O(tau^4) scaling above), a six-line positive control for a symplectic
  routine; **warning 1**: a symplectic scheme with a variable step has no bounded energy error (p.38), so the project's
  adaptive DOP853 and any variable-step symplectic variant cannot be compared on energy drift; **warning 2**: the benefit
  is for separable or Kepler-plus-perturbation Hamiltonians with a constant step and does not extend to time-dependent
  systems (BCR4BP, ER3BP), where there is no conserved energy at all; the project's models are of that kind, so the
  "symplectic integrator" lane of `#929` reduces to the reversibility control, as the OUTSTANDING entry already concludes.
- **`#924`.** Eq. (61) to (62) gives the exact map between the pendulum Hamiltonian and the standard map used for the FLI
  and GLI controls (tau^2 = k): the standard-map control and a pendulum-with-symplectic-Euler control are the same system.
- The page-37 statement that SI4 is worse over short times and better over long ones argues for choosing an integrator by
  the **length** of the run; the project's cycler arcs (a few revolutions) are in the regime where RK-type methods are
  competitive. INFERRED.

## 8. Hagel and Trenkler, "A computer aided analysis of the Sitnikov problem", CMDA 56:81-98 (PDF pp.88-105)

Outside the requested list; included because it prints sourced numbers for the Sitnikov equation of section 6. READ the
abstract, sections 2 to 3, Tables I, II, V and the conclusion; sections 4 to 5 (the construction of the approximate first
integral by Dirac-delta pulse discretisation, p.89 to p.97) were read at the level of the equations printed and the
Tables, not rederived.

### 8.1 Equations (READ pp.82 to 85, image-checked)

Same equation in Wodnar's form (eq. 6 of the chapter: T'' + [e cos phi + (1/4 + T^2)^(-3/2)]/(1 + e cos phi) T = 0 = Dvorak's
eq. 7). P(T) = (1/4 + T^2)^(-3/2) has the Taylor series 8 - 48 T^2 + 240 T^4 - 1120 T^6 + 5040 T^8 - 22176 T^10 + 96096 T^12 + ...
(eq. 9; I checked every coefficient against 8 (1 + 4 T^2)^(-3/2)), with radius of convergence 0.5; Chebyshev (Remez) form
to order 12 with coefficients "rounded to integer numbers" (eq. 10): **P(T) = 8 - 47 T^2 + 203 T^4 - 616 T^6 + 1168 T^8 - 1206
T^10 + 512 T^12**, valid for |T| <= 0.8 (maximum error on that interval 0.046, COMPUTED, about 3 percent of P(0.8) = 1.19; the
polynomial is wildly wrong at T = 0.9: 2.63 against 0.92). Linearised problem (eq. 12): T'' + (8 + e cos phi)/(1 + e cos phi) T = 0,
a Hill equation; monodromy R = M(2 pi); stable iff |Tr R| <= 2; characteristic exponent cos mu = Tr R/2; average frequency
Q = mu/(2 pi) = (1/2 pi) integral of d phi/w^2(phi) (eq. 23).

### 8.2 Table I as printed (p.86, image-checked) and my recomputation

Linearised Sitnikov problem: matrix R = [[r1, r2], [r3, r4]], amplitude function w(0), w(pi)/w(0), Tr R, Q. Negative e means
the primaries start at their greatest separation (phi = 0 at apocentre).

| e | r1 | r2 | r3 | r4 | w(0) | w(pi)/w(0) | Tr R | Q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -0.99 | 0.9325 | -0.0135 | 9.6365 | 0.9325 | 0.1936 | 3.5687 | 1.8650 | 4.9381 |
| -0.80 | -0.5498 | 0.1373 | -5.0817 | -0.5438 | 0.4054 | 1.6663 | 1.0995 | 3.3426 |
| -0.60 | 0.9563 | 0.0675 | -1.2706 | 0.9563 | 0.4798 | 1.3737 | 1.9125 | 3.0472 |
| -0.40 | 0.8534 | -0.1455 | 1.8667 | 0.8534 | 0.5297 | 1.2120 | 1.6117 | 2.8992 |
| -0.20 | 0.5777 | -0.2606 | 2.5563 | 0.5777 | 0.5651 | 1.0960 | 1.1554 | 2.8480 |
| 0.00 | cos(2 pi sqrt8) | sin(2 pi sqrt8)/sqrt8 | -8 r2 | cos(2 pi sqrt8) | 8^(-1/4) | 1.0000 | 0.9461 | sqrt8 |
| 0.20 | 0.5777 | -0.3131 | 2.1280 | 0.5777 | 0.6193 | 0.9124 | 1.1554 | 2.8480 |
| 0.40 | 0.8534 | -0.2139 | 1.2703 | 0.8534 | 0.6406 | 0.8249 | 1.7068 | 2.9127 |
| 0.60 | 0.9563 | 0.1271 | -0.6733 | 0.9563 | 0.6591 | 0.7279 | 1.9125 | 3.0472 |
| 0.80 | -0.5498 | 0.3812 | -1.8303 | -0.5450 | 0.6756 | 0.6002 | -1.0995 | 3.3426 |
| 0.99 | 0.9325 | -0.1007 | 0.4458 | 0.9325 | 0.6894 | 0.2801 | 1.8650 | 4.9381 |

COMPUTED (a converged integration of the linearised equation, 1e-13, and separately the project's `core.er3bp` STM; sections
12.6 and 12.7): Tr R at e = 0 is 2 cos(2 pi sqrt8) = 0.94614 (printed 0.9461, the row is a closed form); at e = +-0.20,
+-0.40, +-0.60 the printed entries match to 3e-4 absolute and at e = +-0.80 to 3e-3 absolute (largest differences: r2 0.0675
against 0.0673 and r3 -1.2706 against -1.2703 at e = -0.60; r1 -0.5498 against -0.5494 and r3 -5.0817 against -5.0842 at
e = -0.80), a pattern consistent with the paper's fixed-step RK4 whose error grows with |e|; the exceptions are the
following printed slips:

- Row e = -0.40 is not the mirror of e = +0.40: printed Tr R = 1.6117, Q = 2.8992, w(0) = 0.5297, but the paper's own text
  says Q is exactly symmetric in e (p.86), and the converged values are Tr R = 1.7068, Q = 2.9127, w(0) = 0.5284 (the +0.40
  row matches); the other printed entries of the row (r1, r2, r3) match. The entries Tr R, Q, w(0) of the e = -0.40 row are
  slips.
- Rows e = +-0.80: printed r4 = -0.5438 and -0.5450 against r1 = -0.5498 and, with det R = 1 and equal diagonals, my r4 =
  r1 = -0.5494; printed Tr R is +1.0995 at e = -0.80 and -1.0995 at +0.80 against my -1.0989 for both (Tr R must be
  symmetric in e: the printed sign at -0.80 is a slip).
- Rows e = +-0.99 do **not** match: printed r1 = 0.9325, Tr R = 1.8650, Q = 4.9381, w(0) = 0.1936 and 0.6894; my values
  r1 = 0.9745, Tr R = 1.9490, Q = 4.9640, w(0) = 0.1926 and 0.6896. The paper's method (RK4 with step 2 pi/200, p.86) is not
  converged near |e| = 1; my Q = 4.9640 agrees with the **Table V entry at e = 0.99, T0 = 0: 4.96398**, so the paper's own
  other table has the converged value.
- w(pi)/w(0) at e = -0.99: printed 3.5687, mine 3.581; r3 at e = -0.99: printed 9.6365, mine 6.0509.

Project code (CHECKED, `core.er3bp.propagate_er3bp` with `with_stm=True`, mu = 0.5, state (0,0,0,0,0,0), integrated over
one period from f = 0 or f = pi, the entries STM[2,2], STM[2,5], STM[5,2], STM[5,5]): (r1, r2, r3, r4) = (0.5777, -0.3131,
2.1280, 0.5777) at e = +0.20 and (0.5777, -0.2606, 2.5564, 0.5777) at e = -0.20 (printed 2.5563); +-0.40:
(0.8534, -0.2139, 1.2704, 0.8534) and (0.8534, -0.1455, 1.8669, 0.8534); +-0.60: (0.9563, 0.1270, -0.6730, 0.9563) and
(0.9563, 0.0673, -1.2703, 0.9563); +-0.80: (-0.5494, 0.3813, -1.8308, -0.5494) and (-0.5494, 0.1373, -5.0842, -0.5494);
+-0.99: (0.9745, -0.1067, 0.4719, 0.9745) and (0.9745, -0.0083, 6.0509, 0.9745). Apart from the paper's slips above, this
reproduces the printed table for |e| <= 0.80. **A test-ready, sourced, independent check of the project's STM in the
eccentric problem**, with the e = +-0.99 and the slip entries excluded and the tolerances above (3e-4, and 3e-3 at |e| = 0.8).

### 8.3 Equation (26) is probably misprinted

Eq. (26), p.87: Q(e) = sqrt8 [1 + (24/121) e^2 + O(e^4)], "reproduces the Q values in Table I for |e| < 0.6" (Hagel 1992).
At e = 0.20 it gives 2.8509 against Table I's 2.8480; at e = 0.40 2.9182 against 2.9127. The second-order coefficient c in
Q = sqrt8 (1 + c e^2): converged numerics give c = 0.16939 (e = 0.02), 0.16959 (0.05), 0.17030 (0.10), tending to **0.16935
= 21/124**; an independent perturbation derivation (g = 8 - 7 e cos phi + 7 e^2 cos^2 phi - ...; shift of Omega^2 =
7/2 e^2 - (49/4)(2/31) e^2 = 2.7097 e^2, giving c = 2.7097/16 = 21/124; derived by hand from the equation as printed, INFERRED, and consistent with the numerics) agrees. 24/121 and 21/124 differ
by a transposition of digits; **eq. (26) as printed does not reproduce Table I** (for e = 0.2, 21/124 gives 2.8476 against
2.8480). Do not code eq. (26) from this page.

### 8.4 Table V (amplitude-dependent frequency) and Table II (integral coefficients)

Table V (p.98, image-checked): Q(e, T0) for T(0) = T0, T'(0) = 0 from eqs. (61) to (63) (the truncated first integral).

| T(0) | e = 0 | 0.2 | 0.4 | 0.6 | 0.8 | 0.99 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.00 | sqrt8 | 2.84802 | 2.91273 | 3.04723 | 3.34258 | 4.96398 |
| 0.10 | 2.77598 | 2.80674 | 2.86849 | 3.00643 | 3.29348 | 4.89341 |
| 0.20 | 2.67461 | 2.69700 | 2.76219 | 2.89060 | 3.17505 | 4.73265 |
| 0.30 | 2.51413 | 2.53423 | 2.60368 | 2.72519 | 2.99919 | 4.60023 |
| 0.40 | 2.30821 | 2.34725 | 2.42091 | 2.52520 | 2.81878 | 4.50321 |

The T(0) = 0 row is the exact linear Q and reproduces (2.84802, 2.91273, 3.04723, 3.34258, 4.96398 against my 2.8480,
2.9127, 3.0472, 3.3426, 4.9640). **The nonlinear rows are approximations that deteriorate with amplitude**: at e = 0 the exact
frequency of T'' + 8 T (1 + 4 T^2)^(-3/2) = 0 is, by quadrature and by integration (COMPUTED), 2.76655, 2.59971, 2.37077,
2.12203 for T0 = 0.1, 0.2, 0.3, 0.4; the printed 2.77598, 2.67461, 2.51413, 2.30821 are higher by 0.34, 2.9, 6.0 and 8.8
percent. Direct integration of eq. (6) from T(0) = T0, T'(0) = 0 at e = 0.2, 0.4, 0.6, 0.8, 0.99 (200 revolutions, mean
frequency from zero crossings) gives, for T0 = 0.10, 0.20, 0.30, 0.40: e = 0.2: 2.79117, 2.63647, 2.42220, 2.18655; e = 0.4:
2.85948, 2.71414, 2.51092, 2.28451; e = 0.6: 2.99684, 2.85855, 2.66362, 2.44413; e = 0.8: 3.29455, 3.16186, 2.97303,
2.75752; e = 0.99: 4.91856, 4.79055, 4.60241, 4.38914. Against the printed Table V these are lower for e <= 0.8 (by
0.03 percent at e = 0.8, T0 = 0.1, up to 8.8 percent at e = 0 and T0 = 0.4) and mixed at e = 0.99 (+0.5, +1.2, +0.05 and
-2.5 percent for T0 = 0.1 to 0.4). The authors themselves describe the
procedure as "semi convergent" (p.95). **Use Table V only at T0 = 0 (exact) or as a check on order of magnitude**, not
as a sourced nonlinear frequency.

Table II (p.94, integral coefficients): the printed e = 0 column (y^2: 8.000, v^2: 1.000, y^4: -8.310, y^6: 8.463) is
reproduced by transforming the Chebyshev potential of eq. (10): 2H = T'^2 + 8 T^2 - (47/2) T^4 + (203/3) T^6, T = 8^(-1/4) y,
dividing by w0^2 gives y^4 coefficient -(47/2)/sqrt8 = -8.3085 and y^6 coefficient (203/3)/8 = 8.458 (COMPUTED; printed
-8.310 and 8.463, a 0.02 to 0.06 percent difference, attributable to the unrounded Chebyshev coefficients; INFERRED).
Tables III and IV (residual variations of the integral in percent against order and amplitude, e = 0.6 and 0.99, pp.96
to 97): the "semi-convergence" (best order 6 to 8; order 10 and 12 get worse, up to 145 percent at T0 = 0.4) is the content;
no further check made.

### 8.5 Techniques applicable

- **A converged, independent check of the project's eccentric-problem STM** (section 8.2), in the one sub-problem where
  everything closes in elementary form (e = 0: Tr R = 2 cos(2 pi sqrt8)). It tests the pulsating-frame equations and the
  `er3bp_stm_eom` variational block for the z degree of freedom only; it does not test the planar blocks.
- The paper's statement that a Hill equation with this g(phi) has **no instability islands for |e| < 0.99999** is atypical
  of Mathieu-type equations (p.87): for any classifier built under `#931` on the elliptic monodromy, this is a case where
  "|Tr R| < 2 for all e" is the right answer and a classifier that flags instability is wrong.
- For `#924`: finite-time exponents on this eq. (7) have the exact linear monodromy as ground truth at small T.

## 9. Laskar, "Frequency analysis of a dynamical system", CMDA 56:191-196 (PDF pp.194-199)

Outside the requested list (it bears on chaos indicators for `#924`). READ in full; Fig. 1 and the maps image-checked.

Method (READ pp.191 to 192): approximate a numerically obtained complex function f(t) on [-T, T] by a quasi-periodic sum
sum a_k exp(i sigma_k t): find sigma_1 as the maximum of |<f(t), exp(i sigma t)>|, with scalar product (1/2T) integral f(t) g(t)-bar
chi(t) dt and the Hanning weight **chi(t) = 1 + cos(pi t/T)** (the printed normalisation is (1/2T) integral chi dt = 1);
subtract, repeat, with Gram-Schmidt of the exp(i sigma_k t) set. For a regular orbit the frequency is the same on every
window [t, t+T]; for a chaotic one it drifts (diffusion).

Standard map form here (p.192): **x' = x - a sin y, y' = x' + y (mod 2 pi)**, a different convention from the 1993 GLI
chapter (x1 = x0 + a sin(x0 + y0), y1 = x0 + y0) and from Yoshida (p' = p + k sin q, q' = q + p'); I checked numerically
that the chaotic-sea Lyapunov exponent at |a| = 1.3 is about 0.22 in both forms
(0.2197 at a = +1.3 from (2, 0.5) and 0.2265 at a = -1.3 from (1, 0.3) for the Laskar form; 0.2254 at a = -1.3 from
(2, 0) for the GLI form). The three forms are different parametrisations of the same family; initial conditions and
the sign of the parameter must be converted before reusing any number, and I have not constructed the conversion.

Printed numbers (READ p.193, Fig. 1 and its caption; p.193 text):

- Criterion: if x < x' on the line y = 0 have frequencies nu(x) > nu(x'), there is no KAM curve of irrational rotation number
  between nu(x') and nu(x).
- Fig. 1: golden rotation number nu_g = (3 - sqrt5)/2; x origin x0 = 4.176550; frequency and x units 1e-6; panels (a) a =
  0.9717, (b) a = 0.9718, (c) a = 0.9710, (d) a = 0.9720. T = 12516 for the panel b; "Fig. 1b shows the disappearance of the
  golden curve for a = 0.9718", "very close to and compatible with the value **a_c = 0.971635 derived by Greene (1979)**".
- 4D map (p.194, image-checked): x1' = x1 + a1 sin(x1 + y1) + b sin((x1 + y1 + x2 + y2)/2), y1' = x1 + y1, x2' = x2 + a2
  sin(x2 + y2) + b sin((x1 + y1 + x2 + y2)/2), y2' = x2 + y2 (mod 2 pi), the Froeschle (1972) map; parameter values a1 = a2
  = -1.3, b = 0.01, initial conditions on x1 = x2 = 0, **516 iterations** per orbit, frequency plane (nu1, nu2) of Fig. 2
  (axes 0.12 to 0.18); regular zones are smooth non-distorted curves, the top of the nu1 = 1/6 vertical resonance
  shows regular spacing in nu2 and erratic in nu1 ("Arnold diffusion is probably possible"), a chaotic outer zone at small nu.

Numerics (COMPUTED, Laskar-form map at a = +1.3, two starts: LCE 0.2197 from (2, 0.5) and 6e-6 (regular) from (1, 0.3)).

Techniques applicable (`#924`): frequency analysis is the independent confirming method named in the Froeschle-Lega-Gonczi
1997 digest (cost 2e4 iterations against FLI's 2e3 at a thin layer); for a periodic-orbit family in the project it needs a
quasi-periodic signal, which a cycler arc (a few revolutions) does not provide, so it is practical only for long orbits in
the circular problems. A **statement of its own limitation** (p.191): over a time shorter than the divergence time a
chaotic orbit looks regular; this is the same short-arc caveat as the FLI digest. The standard-map a_c = 0.971635 is a
sourced published number (Greene) that any standard-map calibration can use: below it the golden curve exists.

## 10. Valsecchi, Perozzi, Roy and Steves, "The arrangement in mean elements space of the periodic orbits close to that of the Moon", CMDA 56:373-380 (PDF pp.371-378)

READ in full; the equations (1) to (10) and Fig. 2 image-checked at 130 to 170 dpi. This is not a chapter on Earth-Moon
CR3BP cyclers: the model is a **geocentric lunar orbit perturbed by the Sun on a circular orbit** ("restricted circular
3-dimensional 3-body problem" with the Earth as the central body), with the Moon's mean elements as variables.

### 10.1 Content (READ)

Eight symmetric periodic orbits of 223 synodic months (the Saros) were found earlier (Perozzi et al. 1991; Valsecchi et al.
1991, 1992), each built from a pair of mirror configurations; the restricted circular 3D problem has 16 mirror
configurations, hence 8 distinct periodic orbits for a given set of periods (p.373). Periods (p.374): synodic month
T_Dlambda (period of the angle lambda - lambda'), anomalistic month T_M (period of the mean anomaly M), nodical month T_theta
(period of theta = omega + M). **Observed means T_Dlambda = 29.530589 d, T_M = 27.554550 d, T_theta = 27.212221 d;
223 T_Dlambda ~ 239 T_M ~ 242 T_theta (eq. 1), and with T_omega the period of the argument of perigee, 223 T_Dlambda ~ 239 T_M ~ 3 T_omega
(eq. 2), exact for the periodic orbits.** General condition (eq. 3): N_Dlambda T_Dlambda = N_M T_M = N_omega T_omega with integers.
With n-bar the observed mean mean-motion and n' the Sun's: T_Dlambda = 2 pi/(n-bar - n') (eq. 4); n-bar = n (1 - nu^2) (eq. 5),
nu = n'/n-bar; T_M = 2 pi/(n-bar - varpi-dot) (eq. 6), T_omega = 2 pi/(varpi-dot - Omega-dot) (eq. 7), the secular rates from
Delaunay (1872).

The two series (pp.375 to 376, image-checked; written here for e' = 0 and gamma = sin(i/2), e the mean eccentricity; with
e' the Sun's eccentricity, set to 0 in the paper, the e'-terms are printed in the source and omitted here):

    Omega-dot = -n-bar [ (3/4 - (3/2) g^2 + (3/2) e^2 + (51/8) g^2 e^2 - (21/64) e^4) nu^2
        - (9/32 - (27/16) g^2 - (189/32) e^2 + (27/16) g^4 + (567/16) g^2 e^2 - (675/256) e^4) nu^3
        - (273/128 - (843/128) g^2 - (2739/128) e^2) nu^4 - (9797/2048 - (7185/1024) g^2 - (165411/2048) e^2) nu^5
        - (199273/24576) nu^6 - (6657733/589824) nu^7 + ((45/32) nu^2 + (1935/512) nu^3) alpha^2 ]            (8)
    varpi-dot = n-bar [ (3/4 - 6 g^2 - (3/8) e^2 - (45/4) g^4 + (69/8) g^2 e^2 - (3/32) e^4) nu^2
        + (225/32 - (189/8) g^2 - (675/64) e^2 + (1107/16) g^4 + (81/32) g^2 e^2) nu^3
        + (4071/128 - (3963/32) g^2 - (31605/512) e^2) nu^4 + (265493/2048 - (335403/512) g^2 - (1483665/4096) e^2) nu^5
        + (12822631/24576 - (25291729/16384) e^2) nu^6 + (1273925965/589824 + (352038885/1179648) e^2) nu^7
        + (71028685589/7077888) nu^8 + (32145882707741/679477248) nu^9 + ((45/32) nu^2 + (7425/512) nu^3) alpha^2 ]   (9)

    alpha = a-bar/a' = [ mu nu^2 (1 - nu^2)^2 / (1 + mu) ]^(1/3)         (10)   (mu the mass ratio of the Earth-Moon system to the Sun)

The e'-dependent terms (printed, e'^2 and e'^4 terms in the same brackets) were read and are omitted from the display; the
check below uses e' = 0 as the paper does. Procedure: for fixed n-bar and n' find (e-bar, i-bar) satisfying eq. (3) for
integer triples with a two-dimensional root finder (Press et al. 1986), restricted circular problem, gamma from the
observed i-bar; Fig. 1 plots (e, i) for N_omega <= 10, 20, 40, 80 with the Moon's actual mean elements **e-bar = 0.0549,
i-bar = 5.133 deg** marked; Fig. 2 (p.378) for N_omega <= 115; "city map" structure with "squares" at points of small N_omega;
the Saros orbit (N_Dlambda, N_M, N_omega) = (223, 239, 3) at **e-bar = 0.0623, i-bar = 5.57 deg** is "the centre of the most
conspicuous square"; every plotted point is really a set of 8 orbits, verified by sampling; the largest set found has
**N_Dlambda = 8534, N_M = 9146, N_omega = 115** (the Hipparchus cycle, Steves 1990) (p.377). Effect of the lunar mean motion
(p.377 to 378, Fig. 3 for N_omega <= 50): for nu = 0.9995 nu0 the Saros square is at **e-bar = 0.0871, i-bar = 5.41 deg**; for nu =
1.0005 nu0 at **e-bar = 0.0137, i-bar = 5.72 deg**; for nu = 0.999 nu0 it has "still to enter from the right", for nu = 1.001 nu0
"just been absorbed" from the left; increasing nu raises the inclination slowly and lowers the eccentricity much faster.
Remark (p.378): if a substantial fraction of the Fig. 2 orbits turned out unstable, the figure would be an Arnold web for
the lunar problem; stability was not computed.

OCR note: OCR gave "0.OS71" for the 0.9995 nu0 eccentricity; the page image reads **0.0871** (and the computation below
confirms it).

### 10.2 Reproduction (COMPUTED, section 12.8)

With my transcription of eqs. (8), (9), (10), the observed T_Dlambda, and **two constants not printed in the paper**: the
Sun's mean motion from a sidereal year of 365.256363 d, and mu = 1/328900.56 (Sun to Earth-Moon mass ratio; external
constants, flagged): n-bar = 2 pi/T_Dlambda + n', nu0 = 0.0748013, alpha = 0.0025623.

- At the Moon's mean elements (0.0549, 5.133 deg): T_M = 27.55436 d (observed 27.554550), T_theta = 27.21228 d (observed
  27.212221), T_omega = 2191.9 d.
- Solving eq. (3) for (223, 239, 3): **(e, i) = (0.062277, 5.5709 deg)** against printed (0.0623, 5.57 deg).
- Shifted nu with n-bar fixed (T_Dlambda recomputed): nu = 0.9990 nu0 gives (0.10628, 5.248 deg) (not printed); 0.9995 nu0 gives
  **(0.08707, 5.412 deg)** against printed (0.0871, 5.41 deg); 1.0005 nu0 gives **(0.013672, 5.7235 deg)** against printed
  (0.0137, 5.72 deg); 1.001 nu0 has no root, matching "just been absorbed".
- Arithmetic of eq. (1): 223 T_Dlambda = 6585.32 d, 239 T_M = 6585.54 d, 242 T_theta = 6585.36 d; T_omega from the observed
  months is 1/(1/T_M - 1/T_theta) = 2190.35 d, 3 T_omega = 6571.05 d, 0.2 percent from the others (eq. 2 is "approximately",
  and holds exactly only for the periodic orbits).

The series coefficients beyond nu^6 contribute 0.7 percent (nu^7), 0.2 percent (nu^8) and 0.08 percent (nu^9) of the leading
term of varpi-dot at nu0; the agreement above to four digits in (e, i) means my transcription of the used coefficients is
right to that level but does not independently verify the highest-order coefficients.

### 10.3 Techniques applicable (`#899`, `#884`/`#905`)

- **For `#884`/`#905`**: the paper's model (Moon as a massless body around the Earth, Sun on a circular orbit) is the
  geocentric analogue of the BCR4BP's Sun forcing, in the regime where the Moon's orbit is nearly Keplerian (not the
  Earth-Moon rotating CR3BP in which the project's cyclers live). The **method** transfers as a predictor: find which
  periodic orbits of the forced problem exist by requiring the three secular periods (here synodic, anomalistic, nodical) to
  be commensurable, using closed-form secular rates, then confirm and refine by shooting; the mirror theorem gives 8
  orbits per commensurability (a count to test against the number of orbits found by the project's symmetric shooting).
  The paper's Fig. 2 shows how dense such orbits are in (e, i) space, which is the context for "which orbits survive the
  Sun". INFERRED, not tested.
- **For `#899`**: a published check that a generator of periodic orbits close to the Moon's own orbit reproduces
  eight orbits at (e, i) = (0.0623, 5.57 deg) for 223 T_Dlambda; the paper does **not** print initial conditions of the 8
  orbits (they are in the unpublished report and the submitted A&A paper, Valsecchi et al. 1991, 1992), so it is a model
  and a mean-element test, not an initial-condition control.
- **For `#924`**: the Arnold-web remark (p.378) is the proposal to compute the stability (Floquet and a chaos indicator) of the
  orbits in Fig. 2; nobody in this volume does it.

## 11. Other chapters that bear on the project

### 11.1 Milani and Nobili, "Asteroid 522 Helga is chaotic and stable", CMDA 56:323-324 (PDF pp.323-324)

READ in full (2 pages). A short note (full details in Milani and Nobili 1992, Nature): a variational vector linearised about
the reference orbit with random initial direction, gamma(t) = log(|v(t)|/|v0|); the maximum Lyapunov exponent is the slope of a
linear fit; **Helga: slope 1.45e-4 per year, Lyapunov time about 6900 years** (1/1.45e-4 = 6897, COMPUTED), over 7 My with
Jupiter to Neptune perturbing and initial conditions referred to the barycentre of the inner Solar System. Elements: a 3.628
to 3.632 AU, e 0 to 0.094, i < 3.9 deg; could in principle approach Jupiter to 0.9 AU but **never closer than 1.31 AU**;
the argument 7 lambda - 12 lambda_J completes about 830 revolutions in 7 My and alternates libration and circulation
(7:12 resonance); either the longitudes of perihelion librate or both eccentricities are near their minima; proper
eccentricity e_p = 0.039 with standard deviation 0.0015 (so the protection mechanism is **not** the 7:12 resonance, which
only causes the chaos). Conclusion: "the first clear-cut example in which the Lyapunov exponent has very little
significance"; "stable chaos" = small average change of proper elements over one Lyapunov time. For `#924`: a published
warning that a large indicator is not a stability verdict; the FLI digest already records this (its section 1). Also:
Mercury's eccentricity changes by 0.05 over about 40 Lyapunov times; Pluto computed for 50 Lyapunov times shows no
significant instability.

### 11.2 Benest, "Stable planetary orbits around one component in nearby binary stars II", CMDA 56:45-50 (PDF pp.53-58)

READ in full. Plane elliptic restricted problem, S-type orbits (around one star); the mass parameter mu and the binary's e are
the only parameters; frame rotating-pulsating with the primary at X = 0 and the other star at X = -1; the initial plane
(X0, v0) with X0 > 0, Y0 = 0, U0 = 0 (the symmetric-start argument of Benest 1978); an orbit is "stable" if there is neither
escape nor collision during **100 binary revolutions**; the grid is 0.05 by 0.05 (2420 orbits for Fig. 1). Systems studied:
alpha Cen (e = 0.52, mu = 0.45 or 0.55), Sirius (e = 0.592, mu = 1/3, 2/3), eta CrB (e = 0.28, mu = 0.45, 0.55; the new case);
result: stable orbits to between half and three quarters of the periastron separation and nearly circular orbits (e < 0.1) to
a quarter. Lyapunov exponents "for all the orbits" planned (p.49). Relevance: the method (a systematic scan of a
two-parameter initial plane with a fixed-time survival criterion in the elliptic restricted problem) is the style of the
project's `#908` capture sweeps; no numbers to reproduce.

### 11.3 Kallrath, Schloder and Bock, "Least squares parameter estimation in chaotic differential equations", CMDA 56:353-371 (PDF pp.352-370)

READ at the level of the method summary, the model, Tables I to III and the conclusion text; the optimisation algebra
(sections 2 to 4) was not rederived. The method fits parameters and initial values of an ODE model to a (possibly noisy)
chaotic time series by a **boundary-value-problem approach**: multiple shooting (50 to 100 nodes) or collocation, with the
node values as extra unknowns, continuity as constraints, and a structure-exploiting generalised Gauss-Newton solver
(Bock 1981). Their stated reason (p.354): for systems with positive Lyapunov exponents "the link of an initial value
problem solver" (single shooting, guess-integrate-correct) is not advisable because errors propagate exponentially, whereas
the multiple-shooting formulation tolerates initial trajectories that are discontinuous and far from the data (Fig. 1 and
Fig. 4: an initial guess a = 9.2, b = 1, c = 2.1 is an escape orbit, yet the original parameters (1, 1, -1) are recovered).
The test system is the Henon-Heiles system, generalised to ODEs x1' = x3, x2' = x4, x3' = -a x1 - 2 x1 x2, x4' = -b x2 -
x1^2 - c x2^2 (the printed eq. 5.5 reads "- b x2 - x1^2 - c x2^2" for the last line, consistent with a potential), with
(a, b, c) = (1, 1, -1) the 1964 system; the section of Fig. 2 is x1 = 0, x3 >= 0 and the text says that between E- = 1/12 and
E+ = 1/8 one passes from regular to ergodic behaviour (p.361).

Printed numbers that are usable (READ pp.361 to 365):

| Case | Energy E | Initial values (x1, x2, x3, x4) | Parameters (a, b, c) | Largest Lyapunov exponent |
| --- | --- | --- | --- | --- |
| Table I | 0.125 | (0, 0, sqrt(0.1275) = 0.35707, -0.35) | (1, 1, -1) | lambda_1 = 0.044; lambda_2 and lambda_3 of opposite sign, magnitude of order 1e-7; lambda_4 = -lambda_1 "with five digits accuracy" |
| Table II | 0.12905 | (0, 0, 0.3, -0.41) | (1, 1, -1) | lambda_1 = 0.053 |
| Table III | 0.1849 | (0, 0, 0.43, 0.43) | 1.3 (1, 1, -1), i.e. (1.3, 1.3, -1.3) | lambda_1 = 0.037 |

The energies are the kinetic energies at the origin (E = (x3^2 + x4^2)/2 = 0.125, 0.12905, 0.1849, COMPUTED, matching the printed
values and eq. 5.7). The exponents are in the time unit of the ODE and are printed with two significant figures; the
integration length used for them is not stated in the pages I read. Recovered-parameter accuracy from exact data: absolute
error below 5e-5 (Table I discussion); Fisher factor 3.75 to convert the printed standard errors to confidence intervals;
the integrator error bound 1e-13 and the estimation accuracy 1e-5 (p.361).

COMPUTED (section 12.9): the largest Lyapunov exponent of the same system by renormalised tangent vectors (DOP853, 1e-11,
unit renormalisation interval), finite-time values at t = 2000, 4000, ..., 10000: Table I orbit 0.052, 0.035, 0.033, 0.038,
0.042; Table II orbit 0.018, 0.038, 0.048, 0.052, 0.056 (still drifting); Table III orbit 0.030, 0.039, 0.038, 0.040, 0.040.
So they agree with the printed 0.044, 0.053 and 0.037 to within the finite-time scatter (about 0.005 to 0.01), which is
the right level of agreement for a two-figure printed number: **use them as bands (0.03 to 0.06), not as values to 1e-3**.

Techniques applicable:
- **`#924`**: three Henon-Heiles chaotic orbits with printed largest exponents and exact initial values, in the
  potential that the project's `#924` entry already uses (Henon and Heiles 1964), to check the sign and rough size of any
  chaos indicator, with the Table I note that the Lyapunov spectrum of a Hamiltonian flow comes out as
  (+lambda_1, ~0, ~0, -lambda_1) in this four-dimensional phase space, a ready consistency test for a tangent-vector code
  (the two middle values must be near 1e-7, and the first and last must cancel to five digits).
- **`#905`**: the method's central point is the project's own situation: Leiva and Briozzo's single-shooting failure over
  five or more Sun periods (as scoped in the OUTSTANDING entry) is the failure mode this chapter says single shooting has
  on any chaotic system, and the cure it prints is multiple shooting with the node values as unknowns. This chapter does
  not look for periodic orbits (it fits trajectories to data), so it supports the diagnosis, not a specific algorithm.
  INFERRED mapping.
- Eichhorn's companion chapter (pp.337-351) is the review of generalised least-squares adjustment that this one extends;
  not read beyond the abstract.

### 11.4 Chapters read at abstract or first page only

- Kubala, Black and Szebehely (pp.51-68): Hill's analytic stability criterion governs the closest stable outer orbit for
  0 <= mu <= 0.15, Laplace's (Graziani-Black numerical) criterion for mu > 0.15, the latter good to a few percent (abstract).
- Wodnar (pp.99-101): a three-page note re-examining Sitnikov's 1960 oscillatory-motion proof (equal primaries, restricted
  problem); sequence-shift chaos on infinitely many symbols for most e in (0, 1) near the escape energy; Alekseev's
  generalisation; refers to Dvorak for numerics. No numbers.
- Sidlichovsky (pp.143-152): surfaces of section of third- and fourth-order asteroidal resonances from the truncated planar
  elliptic Hamiltonian, compared with 4th- and 6th-order symplectic integrators on the full Hamiltonian: good agreement for
  small eccentricity, differences from e about 0.3 (abstract). A volume-internal example of the Yoshida schemes in use.
- Robutel (pp.197-199): KAM (Arnold 1963) applied to the planar planetary problem; no numerics (abstract).
- Gilbert, Froeschle and Frisch (pp.263-272): wavelet analysis of points of the standard map's large chaotic region; k = 1.3
  shows voids from islands, k = 10 chaos fills the space (abstract). The standard-map parameters k = 1.3 and 10 recur in
  this volume's chaos chapters.
- Contopoulos (pp.325-336): periodic orbits of two-degree Hamiltonian systems against quantum eigenfunctions (abstract).
- Bendjoya and Slezak (pp.231-262): wavelet transform review (abstract). Marchal (pp.13-26), Hahn et al. (pp.131-142): not
  relevant.
- First page or abstract read, nothing relevant to cycler, second-species or restricted-problem periodic orbits or to chaos
  indicators: Froeschle and Scholl (pp.163-176, numerical experiments in the 3/1 and nu_6 overlapping region, 1 Myr
  integrations in a Sun-Jupiter-Saturn model); Moons and Morbidelli (pp.273-276, action-angle theory of the planar 2/1
  resonance, "very preliminary"); Michtchenko and Ferraz-Mello (pp.121-130, Hildas by numerical integration and Fourier
  analysis, 50000 years); Schubart (pp.153-162, asteroid 903 Nealley at the 2/1 resonance, 110000 yr in a four-body model);
  Erdi and Kovacs (pp.221-230, fourth-order solution of the ideal resonance problem H = B(y) - eps^2 A(y) cos x); Lemaitre
  (pp.103-120, proper elements); Message (pp.277-284, Hyperion); Delhaise and Henrard (pp.285-286, critical inclination of
  Molniya and Tundra orbits); Farinella et al. (pp.287-306, fragments of 6 Hebe); Morbidelli (pp.177-190, introduction to
  KAM, Nekhoroshev and successive elimination of harmonics; a review that might serve as a reference for an averaged-model
  treatment, not read). Not opened at all: Seidelmann (pp.1-12), Inoue (pp.69-70), Lazzaro et al. (pp.395-396).

## 12. What I computed (throwaway; run outside the repository; not committed)

All with Python 3 and SciPy, plus numba (via the project's uv environment) for the double-star case. Scripts are in the
session scratch area and are not kept; the method of each is stated so it can be redone.

12.1 Standard map statistics (section 3.3): tangent vector renormalised each step, the paper's map, sums over n iterations
of ln(norm); moments from the n samples.
12.2 Double-star (section 4.4): CR3BP, 3D, mu = 0.5, the printed initial conditions, Dormand-Prince 5(4) at relative
tolerance 1e-9 and 1e-11 with absolute 1e-12, 3 tangent vectors with Gram-Schmidt at every integer time, t = 20000; first
attempt with fixed-step RK4 at h = 5e-4 failed for d0 >= 0.255 (close encounters at 0.002 to 0.006).
12.3 Hadjidemetriou rows (section 5.3): inertial two-primary integration (Sun at -mu r(t), Jupiter at (1-mu) r(t) with r from
Kepler's equation for a = 1, n = 1), start at perihelion or aphelion, rotating-frame closure after 2 pi, both omega
conventions; scan in e' for the failing row.
12.4 Sitnikov (section 6): eq. (1) and (7) integrated side by side; single-orbit FTLE per revolution with the variational
equation of (7), DOP853 1e-10, 300 revolutions.
12.5 Yoshida (section 7.2): RK4 energy factor of the oscillator; iteration of eq. (61) with eq. (63).
12.6 Hagel-Trenkler Table I (section 8.2): converged integration of the linearised Hill equation (1e-13), amplitude function
from the Courant-Snyder form with beta(0) = r2/sin(mu), psi = integral of d phi/beta; second-order coefficient by finite
difference in e and by perturbation.
12.7 Project code: `cyclerfinder.core.er3bp.propagate_er3bp(with_stm=True)` at mu = 0.5, e in {+-0.2, ..., +-0.99}.
12.8 Valsecchi (section 10.2): eqs. (8) to (10) transcribed from the page image and evaluated at e' = 0, two-dimensional root
find with SciPy fsolve, scan over nu.
12.9 Kallrath et al. (section 11.3): Henon-Heiles tangent-vector Lyapunov exponent, DOP853 1e-11, t = 10000.

## 13. Test-ready numbers (for the tasks named)

Per the project's no-circular-goldens rule the two tables are separate: Table A holds expected values **printed in the
volume** (with page), the only legitimate goldens; Table B holds values from my own runs, which are **not goldens**: they
are sanity bands that show the printed statement is reachable and must not be written into a test as expected values.

### Table A: printed expected values

| Task | What | Source | Printed expected value and test condition |
| --- | --- | --- | --- |
| `#912`, `#931` | Planar elliptic restricted problem, Sun-Jupiter mu = 0.00095387535: initial conditions (x0, ydot0, e', phase) of symmetric periodic orbits, 43 rows (39 distinct), tables in section 5.2 | Hadjidemetriou pp.211, 216, 217 | each row is a periodic orbit of period 2 pi starting perpendicular to the x-axis with Jupiter at the printed apse, in a frame rotating at the instantaneous angular velocity of the Sun-Jupiter line (the printed ydot0 is relative to that frame). Tolerance for a test: 2e-5 in (x0, ydot0) for the 19 rows listed in section 5.3 (a Newton-corrected exact orbit is the reference); 6e-3 for the remaining rows; the 2:1 row printed e' = 0.010 is run at e' = 0.100 |
| `#912` | Bifurcation orbits of the circular problem: 2:1 orbit A x = 0.1657759, ydot0 = 3.057715, e = 0.72, period 2 pi (p.210); 3:1 (0.479420, 0.962393) period pi and (-0.866333, 0.384417) e = 0.798 period 2 pi (p.212); 4:1 (0.395733, 1.190494) period 2 pi/3, (0.299295, 1.733614) e = 0.243, (0.077836, 4.700496 on p.215 and 4.700478 on p.217) e = 0.801 | Hadjidemetriou pp.210, 212, 215 | periodic in the circular problem; the printed e at e' = 0 is the osculating eccentricity at t = 0 (0.735 for the 2:1 orbit, not 0.72) |
| `#931` | Stability by family: 2:1 I_e stable, II_e unstable; 3:1 I_e stable and the other three unstable; 4:1 I_c and II_c unstable, I_e and II_e stable; III_e and IV_e not determined. Unit pair of the circular monodromy splits for e' > 0 | Hadjidemetriou pp.206, 212, 215, 217 | a classifier on the elliptic monodromy returns these labels (2:1 labels are for the circular-problem families; the paper warns its own stability computations for the 2:1 orbits were not accurate) |
| `#924`, `#912` | Linearised Sitnikov monodromy R = M(2 pi), mu = 0.5, state at the origin of the pulsating frame, rows e = +-0.20, +-0.40, +-0.60, +-0.80 of Table I (section 8.2); Tr R at e = 0 is 2 cos(2 pi sqrt8) printed as 0.9461 | Hagel-Trenkler p.86 | agreement 3e-4 absolute for |e| <= 0.6 and 3e-3 at |e| = 0.8, excluding the printed slips of section 8.2 (row e = -0.40 columns Tr R, Q, w(0); the sign of Tr R and r4 at e = +-0.80); the e = +-0.99 rows are excluded (not converged in the source); |Tr R| < 2 for all |e| <= 0.99 |
| `#924` | Exact linear frequency Q(e) at T(0) = 0: 2.84802, 2.91273, 3.04723, 3.34258, 4.96398 at e = 0.2, 0.4, 0.6, 0.8, 0.99; sqrt8 at e = 0 | Hagel-Trenkler p.98 (Table V, first row) | the printed first row of Table V |
| `#924` | Sitnikov Poincare classification, cells of Table I | Dvorak p.78 | chaotic "k" at (e, T) = (0.33, 0.99), (0.51, 0.94), (0.66, 0.82); invariant curves "0" at the cells marked 0, islands "2" at (0.33, 0.65), (0.51, 0.85); no thresholds printed; one printed cell (0.66, 0.80) "k" is thin and was not reproduced by a short run (Table B) |
| `#924` | 3D restricted problem, mu = 0.5, printed initial conditions of section 4.1 | Lohinger et al. pp.316, 318 | regular (means near zero) for d0 < 0.25, chaotic for d0 >= 0.255; sampled d0 = 0.10, 0.13, 0.22, 0.24 regular and 0.255, 0.27, 0.31, 0.35, 0.40, 0.50 chaotic. No number is printed for the exponents |
| `#924` | Standard map x1 = x0 + a sin(x0 + y0), y1 = x0 + y0: a = -1.3, (x0, y0) = (1, 0) invariant curve (m tends to 0, sigma constant, distribution bimodal at 1e5 iterations); (2, 0) chaotic (m non-zero, settles after 1e4 iterations, sigma after 1e2); a = -0.1 circulation; a = -10 strong chaos | Froeschle et al. pp.309 to 313 | the qualitative statements only; sigma tending to about 0.32 as x0 tends to 1 is a figure reading, not a printed number |
| `#924` | Henon-Heiles largest Lyapunov exponents: E = 0.125, (0, 0, 0.35707, -0.35): 0.044; E = 0.12905, (0, 0, 0.3, -0.41): 0.053; E = 0.1849, (0, 0, 0.43, 0.43), parameters (1.3, 1.3, -1.3): 0.037; lambda_2 and lambda_3 of order 1e-7 and lambda_4 = -lambda_1 to five digits (Table I) | Kallrath et al. pp.363 to 365 | two printed significant figures, finite-time values; a test should use a band (the scatter in my runs is about 0.01) |
| `#924` | Standard-map critical parameter a_c = 0.971635 (Greene); golden curve gone at a = 0.9718, T = 12516 | Laskar p.193 | Greene's printed value |
| `#929` | RK4 on the unit oscillator multiplies p^2 + q^2 per step by 1 - tau^6/72 + ...; Euler by 1 + tau^2; symplectic Euler orbit through (1, 0) lies on q^2 + p^2 + tau p q = 1; eq. (63) is conserved under eq. (61) with error O(tau^4); Ruth order-3 and Forest-Ruth order-4 coefficients (eqs. 28 and 37), x0 + 2 x1 = 1 and x0^3 + 2 x1^3 = 0 | Yoshida pp.28, 32, 33, 35, 41 | the printed formulas |
| `#899`, `#884`, `#905` | Observed mean periods T_Dlambda = 29.530589 d, T_M = 27.554550 d, T_theta = 27.212221 d; 223 T_Dlambda ~ 239 T_M ~ 242 T_theta ~ 3 T_omega (exact for the periodic orbits); Saros orbit (223, 239, 3) at (e, i) = (0.0623, 5.57 deg); shifted nu = 0.9995 nu0: (0.0871, 5.41 deg); 1.0005 nu0: (0.0137, 5.72 deg); the Moon's mean elements (0.0549, 5.133 deg); largest set N = (8534, 9146, 115) | Valsecchi et al. pp.374 to 378 | the printed values, with the series of eqs. (8) and (9); the series needs two external constants not printed in the paper (the Sun's mean motion and mu), so a test must declare them |

### Table B: sanity values from my runs (not goldens)

| What | My value | Note |
| --- | --- | --- |
| GLI statistics, standard map, n = 1e6 | a = -1.3, (1, 0): m = 1.2e-5, sigma = 0.3355, gamma_1 = 0.512, gamma_2 = -1.23; a = -1.3, (2, 0): m about 0.22, sigma = 0.557, gamma_1 = -0.68, gamma_2 = -0.52; a = -0.1: sigma = 0.0815; a = -10: m = 1.62, sigma = 0.771 | section 3.3; chaotic values scatter by 10 percent between n = 1e4 and 1e6 |
| Lohinger double star, first-vector mean of ln alpha at t = 2e4 | regular d0 <= 0.24: 7e-4 to 9e-4; chaotic d0 >= 0.255: 0.10 to 0.19 | section 4.4; depends on tolerance by up to 29 percent for chaotic cases; the adaptive integrator is required |
| Sitnikov finite-time exponent per revolution, 300 revolutions | regular cells 0.011 to 0.030; chaotic cells 0.07 to 0.27 | section 6.3; thresholds are mine; 15 of 16 sampled cells agree with Table A cells |
| Henon-Heiles finite-time exponent at t = 1e4 | 0.042, 0.056, 0.040 for the three cases | section 11.3 |
| Saros solution with my transcription of eqs. (8) to (10) | (0.062277, 5.5709 deg); 0.9995 nu0: (0.08707, 5.412 deg); 1.0005 nu0: (0.013672, 5.7235 deg); external constants used: sidereal year 365.256363 d, mu = 1/328900.56 | section 10.2 |
| Linearised Sitnikov monodromy, converged | agrees with the `core.er3bp` STM; the values at e = +-0.99 differ from the printed ones | section 8.2; **not to be used in place of the printed rows**; if a converged reference is wanted for e = +-0.99 it must be independently sourced |

## 14. Follow-ups and open items

1. `#912`/`#931`: build the 43-row closure test (section 5.3 recipe, pulsating-frame conversion for the project's frame; assert
   closure to the Newton-correction table, not to 1e-6) and record the corrected e' = 0.100 row; the 4:1 bifurcation point
   ydot0 = 4.700496 (p.215) against 4.700478 (p.217) is not resolvable (flat valley) and the 2:1 e = 0.72 (p.210) against
   0.735 (p.211) is an unresolved printed difference (0.735 is the osculating value at t = 0).
2. `#924`: add the double-star regular and chaotic pair (4.4) and the Sitnikov Table I cells (6.3) as calibration cases for
   the project's CR3BP and ER3BP chaos indicators; neither source prints a threshold.
3. The Hagel-Trenkler printed slips (e = -0.40 row, e = +-0.80 row, e = +-0.99 rows, eq. 26) should be recorded under `#910`
   (printed items not reproduced).
4. Dvorak Table II is inconsistent with Table I; recompute before any use.
5. Not done: a full Table I reproduction of Dvorak (60 by 11 grid at 1000 revolutions with an island classifier); the
   project's own `core.cr3bp` STM through the double-star orbits (my reproduction used a private integrator); the
   Valsecchi orbits' initial conditions are in unpublished reports and the submitted A&A paper (acquire if `#899` wants a
   Moon-orbit control).
6. Chapters not read beyond an abstract and possibly of interest: Morbidelli (pp.177-190, perturbation methods and KAM),
   Eichhorn (pp.337-351, generalised least squares, the companion of section 11.3), Gilbert et al. (pp.263-272, standard
   map at k = 1.3 and 10).
