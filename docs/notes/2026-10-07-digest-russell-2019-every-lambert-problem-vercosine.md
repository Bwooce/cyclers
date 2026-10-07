# Digest: Russell 2019, "On the Solution to Every Lambert Problem" (#960 batch 37)

R. P. Russell (University of Texas at Austin), "On the Solution to Every Lambert Problem", Celestial Mechanics and
Dynamical Astronomy 131(11), article 50 (2019), doi 10.1007/s10569-019-9927-z (Crossref: CMDA vol. 131, issue 11,
article-number 50, online 2019-10-30, print 2019-11; the Zenodo v2 records cite it as "Vol. 131, No. 50, pp. 1-33").
Received 19 October 2018, accepted 24 September 2019 (p.1, image).
- Supplied file `bed06c77-19preprint_CelMech_interpolatedVercosineLambert.pdf`, 37 pp., md5
  cc60905cbc5f2a698119a01636ee7485. LaTeX author preprint (pdfTeX, created 2019-10-15). Red running header on every
  page: "draft version Oct. 14, 2019, to appear in Celestial Mechanics and Dynamical Astronomy". The p.1 footnote
  (image) says: "This manuscript is a post-peer-review, pre-copyedit version of an article published in Celestial
  Mechanics and Dynamical Astronomy. The final authenticated version is available online at
  http://dx.doi.org/10.1007/s10569-019-9927-z".
- **Preprint vs journal:** this is the accepted manuscript before copy-editing. The journal version has 33 pages
  (Zenodo citation), this one has 37, so the page numbers differ. I cannot tell whether the content differs: I did
  not see the journal version. All page numbers below are preprint pages (PDF page = printed page).
- **Proposed corpus filename:**
  `cyclers_pdf/papers/russell-2019-solution-every-lambert-problem-cmda-131-50-doi-10.1007-s10569-019-9927-z-accepted-manuscript.pdf`
- **How I read it:** whole text from the text layer, then every equation, algorithm and number quoted below on
  200 dpi page images: pp.1, 4 (footnote i), 5 (Eq. 1), 6, 7 (Fig. 2 labels), 8 (Table 1), 10 (Eqs 8-15), 12
  (Eqs 22-29), 13 (Eqs 30-35), 14 (Eqs 36-38), 15 (Eq. 39, Fig. 5 legend), 16 (Fig. 6 legend, Fig. 7 labels,
  Eqs 42-43), 17 (Eq. 46, Alg. 1), 18 (Eq. 47), 19 (Alg. 2), 20 (Eqs 49-51), 23 (Eqs 53-54, Alg. 3), 24 (Alg. 4),
  26 (Alg. 5), 27 (Alg. 6), 28 (Table 2), 30 (Gooding times), 31 (Eq. 55, test spec), 32 (Table 3).
  Checks: `check_russell2019.py` -> `check_russell2019.out` (both in this folder).
- **Code consumer:** `src/cyclerfinder/core/lambert.py`.
- Wanted list: not listed. Successor paper (Russell 2022 JGCD, below) also not listed and not held.

## 0. Verdict

**A survey of the whole Lambert solution space plus an improved vercosine formulation (iteration variable k)
and a precomputed 2D biquintic-spline initial guess (ivLam). The formulation is fully printed; I transcribed it
and it reproduces every legend value of the two worked geometries and agrees with lambert.py to 1e-13. There is no
numerical stress-case table: the only published input/output cases are two planar geometries (Figs 1/5/6).**

- **What it is.** Sec. 2-3: the Lambert parameter space reduces, by scaling, to two continuous parameters per
  (direction d, signed revolution count N~), and to a rectangle after a shift by T_min (Eq. 49). Sec. 4: the
  vercosine time equation (Eq. 27), singularity-free at theta = pi for the root solve, with series patches at the
  parabola and at k = 0 and a large-k series for short long-way hyperbolas. Sec. 5-6: biquintic interpolation of
  the solution k(tau, Gamma) for N = 0..100, then 0 or 1 unguarded third-order Householder step.
- **What it gives the project.**
  1. A candidate THIRD independent Lambert formulation (the cross-check today uses lamberthub's Izzo 2015 and
     Gooding 1990). The formulation needs only Eqs 25, 27-35 and Algs 2-4; the spline files only supply the seed,
     so a test oracle can be built from the paper with no download and no GPL code.
  2. Published goldens (legend values, 5 significant digits): S, tau, eps_min and the parabolic time Tp for two
     geometries, both directions (sec. 3 below). All 14 agree with lambert.py to the last printed digit.
  3. A published solution count: geometry A, T* = 30 TU gives "3 total solutions ... for each direction"
     (p.6, image). lambert.py returns 3 for each direction.
  4. A sourced meaning for lambert.py's multi-rev branch labels. In all 8 checked cases lambert.py `low` = Russell's
     long-period solution (N~ > 0, larger a) and `high` = short-period (N~ < 0) (definition p.6, image). The
     lambert.py comment calls the low/high mapping "empirical".
  5. A finding about an existing dev dependency: `lamberthub.arora2013` (the 2013 vercosine solver, zero-rev only)
     is already installed. At the published parabolic times it disagrees with both lambert.py and my transcription
     by up to 6.3e-6 (geometry B direct), while those two agree to 1e-15. Not investigated further.
- **Capability difference (not a bug).** lambert.py raises `LambertGeometryError` for |dnu - pi| < 1e-9. In the
  vercosine form the root solve is regular at tau = 0 (theta = pi); only the velocities are singular (g -> 0,
  Alg. 2, p.19). Russell's code adds eps_tiny = 1e-150 to g and warns that the exact half-rev velocities are
  corrupted (p.19). So neither code gives exact-pi velocities; that is a property of the problem.
- **Figure-label mismatches (reported, not goldens).** The trajectory panel titles of Figs 2 and 7 print T values.
  Five of the 16 do not match the true values that three independent routes agree on (lambert.py, my vercosine
  transcription, and the Eq. 47 root), each to 6 or more digits. See sec. 3. Every printed T_min is at or above the
  true minimum, by 0.013 to 0.019 TU. That fits markers placed on a sampled curve near a flat minimum, but this is
  a hypothesis. Use only the legend values as goldens.
- **Reference implementation (PROPOSAL context only).** ivLam, Fortran + MATLAB (sec. 4 below). Not downloaded.
- **Catalogue implication:** none. **Code PROPOSAL:** a test-only vercosine oracle (sec. 5).

## 1. The vercosine formulation (all read on the page images)

Inputs: r1, r2 (vectors), T*, d = +1 short way (0 <= theta <= pi) or -1 long way (pi < theta < 2 pi), N~ (signed
revolution count). From direct/retrograde s: d = sgn(s) sgn(r1x r2y - r2x r1y) (Eq. 1, p.5).

- **Free variable** (Eq. 8, p.10): k^2 = vercos(x / sqrt(a)) = 1 + cos(x / sqrt(a)), x the Sundman universal
  variable (dt/dx = r / sqrt(mu), Eq. 9). Inverse (Eq. 14):
  - dE = 2 pi - arccos(k^2 - 1) + 2 N pi, if -sqrt2 <= k < 0;
  - dE = arccos(k^2 - 1) + 2 N pi, if 0 <= k <= sqrt2;
  - dF = arcosh(k^2 - 1), if k > sqrt2.
  Parabola at k = sqrt2; ellipse -sqrt2 < k < sqrt2; hyperbola k > sqrt2 (p.11).
- **Geometry** (Eqs 25-26, 28, p.12): c_theta = r1.r2 / (r1 r2), s_theta = d |r1 x r2| / (r1 r2);
  tau = d sqrt(r1 r2 (1 + c_theta)) / (r1 + r2), or (s_theta / (r1 + r2)) sqrt(r1 r2 / (1 - c_theta)) when
  (1 + c_theta) <= eps_tau (suggested 1e-2); -sqrt2/2 < tau < sqrt2/2. S = sqrt((r1 + r2)^3 / mu).
- **Semi-major axis** (Eq. 23): a = (1 - k tau)(r2 + r1) / (2 - k^2).
- **Time equation** (Eq. 27): T = S sqrt(1 - k tau) [tau + (1 - k tau) W], with (Eq. 29) M = 1 / (2 - k^2) and
  W = dE sqrt(M^3) - k M (ellipse), W = -dF sqrt(-M^3) - k M (k > sqrt2). N enters only through W (p.14).
  Normalised: T~ = T / S = sqrt(p) Z, p = 1 - k tau, Z = tau + p W (Eqs 43, 46).
- **Practical W** (Eq. 30, p.13): five bands, W_e- (k <= -eps0), W_e0 series (|k| <= eps0), W_e+ (N != 0 or
  eps0 < k <= sqrt2 - eps_p), W_p series (|k - sqrt2| <= eps_p, N = 0), W_h (k > sqrt2 + eps_p).
  eps_p = eps0 = 0.02 for double precision (0.0002 for quad), series to 8th order (p.13).
  - W_e0 (Eq. 32): Q/4 - k + 3Qk^2/16 - 2k^3/3 + 15Qk^4/128 - 2k^5/5 + 35Qk^6/512 - 8k^7/35 + 315Qk^8/8192,
    Q = sqrt2 (2N + 1) pi.
  - W_p (Eq. 34): sqrt2/3 - P/5 + 2 sqrt2 P^2/35 - 2P^3/63 + 2 sqrt2 P^4/231 - 2P^5/429 + 8 sqrt2 P^6/6435
    - 8P^7/12155 + 8 sqrt2 P^8/46189, P = k - sqrt2.
  - Check (`.out`, "Series vs closed form"): both series match the closed form to 1e-15 at the band edges
    (N = 0, 1, 5). So the coefficients as I read them are right.
- **Derivatives** (Eq. 36, p.14): W' = (3 W k - 2) M, W'' = (5 W' k + 3 W) M, W''' = (7 W'' k + 8 W') M away from
  the parabola; Eq. 37 gives the W_p series derivatives (I checked by hand that they are the term-by-term
  derivatives of Eq. 34). The W' recursion matches a central difference (`.out`).
- **Parabolic time** (Eq. 39, p.15, image): Tp = S sqrt(1 - tau sqrt2)(tau + sqrt2)/3. (The text layer drops the
  bracket and reads as tau + sqrt2/3; the image is right, and the printed form equals Eq. 27 at k = sqrt2.)
- **Minimum time per revolution** (Eq. 47, p.18): (3 k tau^2 - 3 tau) W + (2 k^2 tau^2 - 4 k tau + 2) W' - tau^2 = 0;
  its root k_Tmin depends only on tau (and N), so T~_min(tau) is a 1D function (Eq. 48), which ivLam interpolates.
  My Eq. 47 roots match a direct minimisation of Eq. 27 to 1e-8 in k in all 8 cases (`.out`).
- **Root solve** (Alg. 3, p.23): f~ = sqrt(p)[tau + p W] - T~*, with closed-form f~', f~'', f~''' (checked against
  finite differences of Eq. 43 at 4 points, `.out`). Alg. 4 (p.24): third-order Householder,
  q_f = -1/f', d1 = f q_f, q_h = d1^2 q_f, d2 = q_h f''/2, d3 = q_h d1 f'''/6 + d2 f f'' q_f^2,
  alpha <- alpha + d1 + d2 + d3.
- **Precision patches.** Short-flight hyperbolas: (a) tau > 0: iterate on p instead of k when N~ = 0 and
  p < eps_p ~ 0.1 (Alg. 6, p.27); (b) tau < 0, k > eps_k ~ 1000: Alg. 1 (p.17) gives Z as a series in kappa = 1/k
  (t1 = -ln2 + 2 ln kappa + 2, t2 = -3 ln2 + 6 ln kappa + 5, Z = t2 kappa^5 - t2 tau kappa^4 + t1 kappa^3
  - t1 tau kappa^2 + kappa). The series matches the closed Z to 13 digits at k = 1000 and 1e-10 at k = 200 (`.out`).
- **Velocities** (Alg. 2, p.19): p_a = p (r1 + r2), f = 1 - p_a/r1, g = S tau sqrt(p), g_dot = 1 - p_a/r2,
  v1 = (r2 - f r1)/(g + eps_tiny), v2 = (g_dot r2 - r1)/(g + eps_tiny), eps_tiny = 1e-150.
- **Initial guess.** The 2019 guess is the spline k(tau, Gamma_j), Gamma_0 = T~*, Gamma_i = T~* - T~_min,i
  (Eq. 49, p.20; Alg. 5, p.26), domain tau in [-sqrt(1/2) + 1e-7, sqrt(1/2) - 1e-7], Gamma_0 in [1e-5, 1e3],
  Gamma_i in [1e-8, 1e3] (Table 2 note, p.28). The older analytic guesses are in Arora & Russell 2013 (p.17; not
  held). At Gamma_j = 1e3 the minimum eccentricity over all geometries is 0.98, 0.96, 0.89, 0.74 for N = 0, 4, 16,
  64 (p.20).
- **Stress cases (Table 1, p.8, qualitative only):** #1 T too small, long way, N = 0; #2 T too small, short way,
  N = 0; #3 T too large, all N; #4 T near T_min,i, N > 0; #5 theta ~ pi (velocity divisor); #6 r1 ~ r2 with theta ~ 0
  or 2 pi, the "only physical singularity". No numbers.

## 2. Accuracy and speed claims (published, all read on images)

- Table 2 (p.28): four coefficient sets, ID 4-7, target rms k error 1e-4 to 1e-7; 16.0 to 65.3 MB for the N = 0
  file, 2.7 to 37.7 MB per N~ != 0 file; ~3,640 patches per MB.
- Velocity test (p.31): r2 = 10^random(-2, 2), theta = random(0, 2 pi), T* = Tp 10^random(-3, 3), mu = 1,
  r1 = [1, 0, 0], 10 million problems per N~, theta kept 0.5 deg away from n pi; error metric Eq. 55 against a
  quad-precision truth.
- Results (p.32): ID 7 + 1 iteration matches truth to ~13 digits rms (~11 worst) for low N, like Gooding; Gooding is
  ~1 digit better at N~ = 0 and loses accuracy at very high N, mainly short-period. ID 7 with 0 iterations: ~7 digits
  rms (~4 worst) zero-rev, ~8 (~6 worst) multi-rev.
- Table 3 (p.32, N = 0..10): low fidelity ID 7, 0 iter, error 3e-7, 3.4x, 920 MB; medium ID 4, 1 iter, 2e-9, 3.0x,
  69 MB; high ID 7, 1 iter, 8e-15, 2.6x, 920 MB. Speedups are relative to Gooding's Fortran.
- Absolute anchor (p.30): Gooding needs 0.84 us (N = 0, one solution) and 2.05 us (N > 0, both solutions) on a
  Xeon E5-2680v3 at 2.50 GHz, Intel Fortran 16.0, O3. Overall 2 to 5 times faster than Gooding (p.31).
- Parallel (p.32): O(1e8) solutions per second with OpenMP on 24 cores.
- Best-case root-solve errors (Fig. 12, p.25) and Figs 14-17 are colour maps only; no tabulated cells.

## 3. Published numbers and the checks

Geometries (Fig. 1 caption p.6; Fig. 6 caption p.16): mu = 1 LU^3/TU^2, r1 = [1, 0, 0].
A: r2 = [-2, 2, 0], T* = 30 TU. B: r2 = [2, 1, 0], T* = 14 TU (from the Fig. 7 labels).
"Direct" = s = +1, "retrograde" = s = -1. Both geometries have r1 x r2 along +z, so direct is short way here.

**Golden candidates (legend values, Figs 5 and 6, pp.15-16, images; formulation-independent):**

| Geometry | S | tau (direct / retro) | eps_min | Tp direct | Tp retro |
|---|---|---|---|---|---|
| A | 7.4908 | +0.23774 / -0.23774 | -0.26903 | 3.3606 | 3.3957 |
| B | 5.8214 | +0.63601 / -0.63601 | -0.43008 | 1.2615 | 2.0812 |

eps_min and both A values of Tp also appear in Fig. 1 (p.6). All 14 values MATCH to the last printed digit
(`.out`, SUMMARY): S and tau from Eqs 25/28; eps_min = -mu/s (s the semi-perimeter); Tp both from Eq. 39 and from
lambert.py (the T at which its N = 0 solution has zero energy; the two routes agree to 1e-8 relative).

**Solution count (p.6, image):** geometry A, T* = 30: N_max+ = N_max- = 1, "3 total solutions exist for each
direction". lambert.py with max_revs = 5 returns exactly 3 per direction.

**Figure-panel labels (Fig. 2 p.7, Fig. 7 p.16): loose checks only.**

| Label | Printed | lambert.py / transcription | Verdict |
|---|---|---|---|
| A direct T(eps_min), pt 1 | 7.94 | 7.9419 | match |
| A retro T(eps_min), pt 8 | 7.98 | 7.9773 | match |
| A direct T_min1, pt 3 | 23.14 | 23.1227 | MISMATCH 0.017 |
| A retro T_min1, pt 10 | 23.16 | 23.1581 | match |
| A direct T_min2, pt 6 | 39.36 | 39.3458 | MISMATCH 0.014 |
| A retro T_min2, pt 13 | 39.4 | 39.3811 | match at 1 decimal; 0.019 if "39.4" means 39.40 |
| B direct T(eps_min), pt 1 | 3.46 | 3.4663 | MISMATCH 0.006 |
| B retro T(eps_min), pt 8 | 4.4 | 4.4098 | match at 1 decimal |
| B direct T_min1, pt 3 | 10.98 | 10.9787 | match |
| B retro T_min1, pt 10 | 11.92 | 11.9188 | match |
| B direct T_min2, pt 6 | 19.02 | 19.0040 | MISMATCH 0.016 |
| B retro T_min2, pt 13 | 19.96 | 19.9462 | MISMATCH 0.014 |
| Tp panels 7/14 | 3.36, 3.4, 1.26, 2.08 | as legend | match |

I suspected my side first: T_min from lambert.py's golden-section search, from my vercosine Eq. 27 minimisation,
and from the Eq. 47 root all agree to 1e-6 TU or better. The legend values of the same figures match exactly.
So the panel labels are the outliers. Do not use them as goldens.

**Other checks (`.out`):**
- My transcription (Algs 3, 4, 2 with a bracketed seed) vs lambert.py at T* for N~ = 0, +1, -1, both geometries,
  both directions: velocity agreement 2e-16 to 4e-13. This is "transcription vs our code", NOT published values
  (the paper prints no velocities). Same for a km/s Earth case (mu = 398600.4418): 1.3e-15 km/s.
- Branch mapping: LP (N~ > 0, larger a) = lambert.py `low` in all 8 cases; SP = `high`. The k positions agree with
  the figure markers (e.g. A direct LP k = 0.90, SP k = -0.16; Fig. 5 points 4 and 5).
- lamberthub.arora2013 vs lambert.py: 0 to 6.3e-6 near the parabola, 1.5e-8 or better at T*. My transcription
  sides with lambert.py at the parabola (1e-15).

## 4. Reference implementation (checked 2026-10-07, read-only `curl -I` / Zenodo API; nothing downloaded)

- Footnote i (p.4, image): "Ryan P. Russell. (2019, October 14). ivLam (Version 1.06). Zenodo.
  http://doi.org/10.5281/zenodo.3479924; Also see http://russell.ae.utexas.edu/index_files/lambert.html".
- doi 10.5281/zenodo.3479924 resolves (HTTP 200). Files: code `ivLamV1p06_737712p63970.zip` (138 kB, Fortran and
  MATLAB), README.txt, and coefficient zips of 66 MB (0-rev hi-res), 931 MB (0-10 rev), 7.7 GB (0-100 rev hi-res),
  538 MB (0-100 rev low-res).
- **Licence mismatch:** the Zenodo metadata for v1.06 says GPL-2.0; its README says "covered under the gpl-3.0.txt
  license file inside the code package".
- Later versions (same concept record 10.5281/zenodo.3479923): ivLam2 2.41 (2021-09-13, record 5196639, GPL-2.0
  metadata, 2.0 MB zip) and ivLam2 v2.50 (2025-12-30, record 18102116, GPL-3.0-or-later, 2.0 MB zip). v2 uses one
  ~1 MB coefficient file, "essentially no limits on flight time or number of revolutions", and adds first- and
  second-order sensitivities; it accompanies Russell 2022, JGCD 45(2):196-212, doi 10.2514/1.G006089 (Zenodo text).
- The author page http://russell.ae.utexas.edu/index_files/lambert.html returns 404; https gave no response.
- No Python version is listed. GPL code must not be vendored into this project without a licence decision.

## 5. PROPOSAL: a third-formulation cross-check for lambert.py (not implemented)

- Write `tests/core/_vercosine_oracle.py` (test-only) from Eqs 25, 27-36 and Algs 2-4 of this paper: about 80
  lines, no coefficient files (seed from a bracketed scalar root, then Householder polish). Bracketing (p.15, p.18):
  N = 0 on (-sqrt2, 1/tau) for tau > 0 (or an open upper bound for tau <= 0); N > 0 split at the Eq. 47 root,
  LP = the side with the larger a. No GPL code involved.
- Golden tests (published values only, preprint pages):
  1. Geometry A direct: S = 7.4908, tau = 0.23774, eps_min = -0.26903, Tp = 3.3606 (Fig. 5 legend, p.15).
  2. Geometry A retrograde: tau = -0.23774, Tp = 3.3957; 3 solutions at T* = 30 (p.6).
  3. Geometry B direct: S = 5.8214, tau = 0.63601, eps_min = -0.43008, Tp = 1.2615 (Fig. 6 legend, p.16).
  4. Geometry B retrograde: tau = -0.63601, Tp = 2.0812 (p.16).
  Assert to half a unit of the last printed digit. Tp is tested through lambert.py by root-solving for zero energy.
- Cross-check-only corners (no published outputs; oracle vs lambert.py vs lamberthub): the Table 2 domain edges
  tau = +-(sqrt(1/2) - 1e-7), Gamma_0 = 1e-5 and 1e3, Gamma_i = 1e-8 (near T_min) and 1e3, N up to 100; plus
  theta = pi +- 1e-6 (the oracle's k stays regular, the velocities do not) and the p.31 random distribution with
  theta kept 0.5 deg from n pi.
- Add a test that pins `low` = long period (larger a) and `high` = short period, citing p.6.
- Look at `lamberthub.arora2013`'s near-parabola error before relying on it as a fourth opinion.

## 6. Citation mining

Held-status checked per work with `ls cyclers_pdf/papers | grep` and `grep CORPUS_INDEX.md` (2026-10-07); wanted
list `2026-10-05-960-wanted-papers.md` checked by name.

| Cited work | Relevance | Held? | Wanted row |
|---|---|---|---|
| Russell & Ocampo 2005, JSR 42(1):138, doi 10.2514/1.5571 (half-rev velocities; r1 = r2 resonant solutions) | medium | HELD (`russell-ocampo-2005-geometric-analysis-free-return-...`) | - |
| Stiefel & Scheifele 1971 | low | HELD (`stiefel-scheifele-1971-...grundlehren-174...`) | - |
| Arora & Russell 2013, AAS 13-728 (original vercosine, initial guesses) | HIGH for the oracle seed | not held | not listed |
| Arora, Russell, Strange & Ottesen 2015, JGCD 38(9):1563 (Lambert partials) | medium | not held | not listed |
| Russell 2022, JGCD 45(2):196, doi 10.2514/1.G006089 (ivLam2, successor; from Zenodo, not cited in the paper) | medium | not held | not listed |
| Gooding 1990, CMDA 48:145 | medium (in use via lamberthub) | not held (the only "gooding" file is Zhang et al. 2024) | not listed |
| Izzo 2015, CMDA 121:1 | medium (in use via lamberthub) | not held | not listed |
| Battin 1977 AIAA J; Battin 1999 textbook | low | not held | Battin 1999 in row 42 |
| Vallado 2013, 4th ed. (lambert.py's Algorithm 5.2 source) | medium | not held (only the 1991 USAFA TR) | not listed |
| Bate, Mueller & White 1971 | low | not held | not listed |
| Bombardelli, Gonzalo & Roa 2018, JGCD 41(3):792 (approximate multi-rev) | low | not held (index hit is a 2025 ERTBP paper) | not listed |
| Pellegrini, Russell & Vittaldev 2014, CMDA 118:355 (F and G Taylor series) | low | not held (held Pellegrini file is the 2016 JGCD STM paper) | not listed |
| Lantoine, Russell & Dargent 2012, ACM TOMS 38(3) (multicomplex derivatives) | low | not held (held Lantoine files are other papers) | not listed |
| Shen & Tsiotras 2003, JGCD 26:50 (two-impulse multi-rev rendezvous) | low | not held | row 42 lists Shen & Tsiotras AAS 03-568 (a related item, not this paper) |
| Ochoa & Prussing 1992 (multi-rev Lagrange) | low | not held (Prussing 2000 is held, a different paper) | not listed |
| Longuski & Williams 1991, CMDA 52:207 | low | not held (other Longuski papers are held) | not listed |
| Torre, Flores & Fantino 2018, Acta Astronaut. 153:26; Torre & Fantino review | low | not held | not listed |
| Ahn, Bang & Lee 2015 JGCD; He, Li & Han 2010 JGCD; Avanzini 2008 JGCD; Der 2011; Kriz 1976; Jezewski 1976; Lancaster & Blanchard 1969 / Lancaster et al. 1966; Thorne 2004; Sun 1971; Nelson & Zarchan 1992; Mahajan & Vadali 2018; Klumpp 1991; Peterson et al. 2010; Escobal 1965; Herrick 1971; Herrick & Liu 1959; Gauss 1857 | low (Lambert method history) | none held | none listed |
| Colasurdo et al. 2014 (GTOC6); Petropoulos et al. 2007 (GTOC1); Healy 2014, Healy et al. 2019; Ottesen & Russell 2017 | low (applications) | none held | none listed |
| W. S. Russell 1995, Appl. Numer. Math. 17:129 (biquintic splines); Psiaki et al. 2019; Junkins et al. 1973; Fornberg 1981 | low (interpolation) | none held | none listed |

If the oracle in sec. 5 is built, Arora & Russell 2013 (AAS 13-728, no DOI) is the one worth acquiring, for its
initial-guess scheme; it is not needed for a bracketed oracle.

## 7. Unresolved

- Journal-vs-preprint content differences: not determinable without the journal PDF.
- Why five figure-panel T labels are 0.006 to 0.019 TU off the true values (hypothesis: markers on a sampled
  curve).
- Why lamberthub.arora2013 loses ~6 digits near the parabola for geometry B direct (not investigated).

*Filed as `cyclers_pdf/papers/russell-2019-on-the-solution-to-every-lambert-problem-cmda-131-50-doi-10.1007-s10569-019-9927-z-author-preprint.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the pre-batch-37 numbering; the list was renumbered in batch 37.*
