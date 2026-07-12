# Digest: Richardson 1980 — Analytic construction of periodic orbits about the collinear points

Single-paper digest. Read 13/13 pages of the user-supplied PDF on 2026-07-12 AET.
Filed to the private corpus as
`richardson-1980-analytic-construction-periodic-orbits-collinear-points-celest-mech-22-241-doi-10.1007-BF01229511.pdf`.

## 1. Header

- **Title (verbatim)**: *Analytic Construction of Periodic Orbits about the
  Collinear Points*
- **Author**: David L. Richardson, Dept. of Aerospace Engineering, University
  of Cincinnati
- **Venue**: *Celestial Mechanics* **22** (1980) 241-253
- **DOI**: 10.1007/BF01229511 (confirmed via CrossRef bibliographic search
  2026-07-12; matches the citing pattern of the project's existing
  `howell-1984-...-doi-10.1007-BF01358403.pdf` filing exactly — same journal,
  same Springer/Reidel retro-DOI numbering scheme)
- **Received / accepted**: 28 March 1979 / 30 August 1979
- **Publisher**: D. Reidel Publishing Co.
- **Length**: 13 pages (pp. 241-253, incl. Appendix I + references)

**Companion paper found, NOT yet in the corpus** (flagging per "if you find
more references you need just ask"): Richardson, D. L. (1980), "A Note on a
Lagrangian Formulation for Motion about the Collinear Points," *Celestial
Mechanics* **22**, 231-236, DOI 10.1007/BF01229509 (same issue, immediately
preceding this paper — confirmed via the same CrossRef query). This is the
paper this digest's paper cites as "(Richardson, 1980)" for the Lagrangian
`L` (Eq. 5) the equations of motion (Eq. 9a-c) are derived from. Not blocking
— this paper reproduces the final equations of motion (Eqs. 9a-c) and the
`c_n` coefficient definitions (Eqs. 8a-b) explicitly, so the derivation is
usable without the companion. Worth acquiring only if a future need arises to
audit the Lagrangian derivation step itself.

## 2. What the paper actually is

**This is THE foundational "Richardson third-order" analytic halo-orbit
construction** — the single most widely-cited closed-form approximation for
CR3BP halo-type periodic orbits about the collinear libration points (L1,
L2, L3), used ubiquitously in the field as an initial-guess generator that a
numerical differential-correction scheme (e.g. Howell 1984's method, already
digested — [[project reference below]]) then refines to full periodicity.
Where Howell (1984) is a *numerical* survey (differential correction +
continuation, tabulating discrete converged orbits across mass ratio), this
paper is the complementary *analytic* half: a closed-form third-order
Lindstedt-Poincaré-style series solution `x(τ), y(τ), z(τ)` parameterized
directly by the desired out-of-plane amplitude `A_z`, requiring no shooting
or iteration to get an initial guess.

## 3. Model and method

### 3.1 Setup (p. 241-244, §1-2)
Standard CR3BP, rotating frame centered at a collinear point (L1, L2, or L3),
x-axis along the syzygy line away from the larger primary, z completes the
right-handed triad (Fig. 1, p. 242). Linearized equations of motion (Eq. 1a-c)
have the classic saddle×center×center structure: in-plane (x,y) has 2 real +
2 imaginary roots, out-of-plane z is simple harmonic. If the in-plane
frequency λ and out-of-plane frequency ν happen to be *equal* (achieved by
choosing amplitudes of "sufficient magnitude" that the *nonlinear*
contributions shift the two linear frequencies together), the linearized
motion becomes genuinely periodic (not quasi-periodic) — this is the
halo-orbit resonance condition, Eq. (4a-c).

### 3.2 Derivation (p. 244-247, §3-4)
The Lagrangian `L` (Eq. 5, from the companion paper) is expanded in Legendre
polynomials of the primaries' perturbing potential (Eq. 7), giving the full
nonlinear equations of motion (Eq. 9a-c) with `c_n` coefficients (Eq. 8a-b)
that depend only on the mass ratio `μ` and the dimensionless `γ_L` (the
libration point's distance ratio to its nearby primary). A **Lindstedt-
Poincaré** successive-approximation procedure is applied through **third
order** in the amplitudes `A_x, A_z`: a frequency correction `ω = 1 + Σω_n`
removes secular terms order by order (`ω_1 = 0` identically; `ω_2` given by
Eq. 17), and — because the z-equation's secular terms can't be removed by a
frequency correction alone — an **amplitude-constraint relationship**
`l1·A_x² + l2·A_z² + Δ = 0` (Eq. 18) plus a **phase-constraint**
`ψ = φ + nπ/2` (Eq. 19, n=1 or 3) close the system. This exact bifurcation
mechanism (switching n=1 vs n=3 via the switch function `δ_n = 2-n`, Eq. 21)
produces the paper's two solution branches, **Class I** (n=1) and **Class
II** (n=3) — the same "two mirror-symmetric branches" convention used
throughout the halo-orbit literature (matches Farquhar & Kamel 1973's
terminology, cited p. 243).

### 3.3 Final closed-form solution (Eq. 20a-c, p. 247)
```
x(τ1) = a21·Ax² + a22·Az² - Ax·cos(τ1) + (a23·Ax² - a24·Az²)·cos(2τ1)
        + (a31·Ax³ - a32·Ax·Az²)·cos(3τ1)
y(τ1) = k·Ax·sin(τ1) + (b21·Ax² - b22·Az²)·sin(2τ1)
        + (b31·Ax³ - b32·Ax·Az²)·sin(3τ1)
z(τ1) = δn·Az·cos(τ1) + δn·d21·Ax·Az·(cos(2τ1) - 3)
        + δn·(d32·Az·Ax² - d31·Ax³)·cos(3τ1)
```
with `τ1 = λτ + φ`. **Every coefficient** (`a21..a32`, `b21..b32`, `d21,
d31, d32`, `l1, l2, s1, s2, k, λ`) is given **in closed form as a function of
`c2, c3, c4` alone** (Appendix I, p. 250-252) — these `c_n` in turn depend
only on `μ` (Eq. 8a-b). This means the whole solution is directly
**algorithmically reproducible from μ and a chosen A_z** with zero numerical
iteration — a genuinely different construction path from every
differential-correction-based halo generator already in the codebase
(`src/cyclerfinder/search/cr3bp_seed_generator.py::lyapunov_seed_3d`,
`src/cyclerfinder/genome/known_corpus_3d.py`'s Howell-1984-anchored family).

## 4. Sourced numeric golden values (Table I, p. 253) — Sun-Earth system

This is the **fully-tabulated, directly reusable golden data** for this
digest. All values below trace VERBATIM to Table I, p. 253 (read via Claude
vision on the table page, not OCR/digitized from a plot — every entry is a
typeset table cell, so this is a genuine sourced value, not a figure-derived
one):

**System constants**: `n1 = 1.99099e-7 rad/sec`, `μ = 3.04036e-6` (Sun-Earth,
*including* the Moon's mass per the paper's own footnote — "the actual value
of μ includes the mass of the Moon"), `a1 = 1.49598e8 km`, `A_z = 1.25e5 km`
(the specific amplitude Table I's derived constants are evaluated at).

**CORRECTION (2026-07-12, #580):** the table below was originally
transcribed via a single Claude-vision pass and contained several errors,
caught while building task #580's golden test by (a) re-extracting the
page's embedded PDF text layer directly (`pdftotext -layout`, exact,
non-OCR) and (b) cross-checking against an independent third-party
implementation of the same paper's Appendix I
(`jacobwilliams/Fortran-Astrodynamics-Toolkit`'s `halo_orbit_module.f90`).
Confirmed errors in the original transcription: `λ_L3` was mistranscribed as
`3.22729` (actually `1.00000`; that value is literally `k_L1`, apparently
copy-pasted into the wrong cell), `k_L1` had a digit transposition
(`3.22729` vs the correct `3.22927`), `c3_L2` was missing its minus sign
(the table's L2 column carries `(-1)^n` through `c_n`, so odd-`n` `c3`, `s1`,
`l1`, `a1`, `d31` are negative at L2), `s2`'s exponent was `e1` instead of
`e-1` for L1/L2, `l1_L3` had a digit transposition (`-1.57177e-5` vs the
correct `-1.57717e-5`), and the `a24`/`b21`/`b32`/`d21`/`d31` rows were
missing a `e-1`/`e-2` scale factor for one or more columns. The table below
is the corrected, doubly-verified transcription (max relative deviation
between this table and this implementation's live output —
`search/cr3bp_seed_generator.py::richardson_halo_coefficients`, see #580 —
is ~4e-6 across all 84 cells, consistent with Table I's 6-significant-figure
publication precision).

| constant | L1 (period 177.704 d) | L2 (period 180.145 d) | L3 (period 365.255 d) |
|---|---|---|---|
| γ_L | 1.00109e-2 | 1.00782e-2 | 9.99998e-1 |
| λ | 2.08645 | 2.05701 | 1.00000 |
| k | 3.22927 | 3.18723 | 2.00000 |
| Δ | 2.92214e-1 | 2.90785e-1 | 2.66029e-6 |
| c2 | 4.06107 | 3.94052 | 1.00000 |
| c3 | 3.02001 | -2.97984 | 1.00000 |
| c4 | 3.03054 | 2.97026 | 1.00000 |
| s1 | -8.24661e-1 | -7.44452e-1 | -1.59141e-6 |
| s2 | 1.21099e-1 | 1.25047e-1 | 6.29433e-6 |
| l1 | -1.59656e1 | -1.48288e1 | -1.57717e-5 |
| l2 | 1.74090 | 1.67369 | 1.40258e-5 |
| a1 (coef) | -8.78563 | -8.52882 | -1.25889e-5 |
| a2 (coef) | 6.86546e-1 | 6.15466e-1 | 1.43702e-6 |
| d1 (coef) | 3.11184e2 | 2.93192e2 | 1.20002e1 |
| d2 (coef) | 1.58787e3 | 1.49800e3 | 7.20008e1 |
| a21 | 2.09270 | -2.05304 | 5.00000e-1 |
| a22 | 2.48298e-1 | -2.51646e-1 | 2.50000e-1 |
| a23 | -9.05965e-1 | 8.96284e-1 | -4.99999e-1 |
| a24 | -1.04464e-1 | 1.06600e-1 | -2.49999e-1 |
| a31 | 7.93820e-1 | 7.80646e-1 | 3.75000e-1 |
| a32 | 8.26854e-2 | 8.36960e-2 | 1.25000e-1 |
| b21 | -4.92446e-1 | 4.91357e-1 | -2.50000e-1 |
| b22 | 6.07465e-2 | -6.27190e-2 | 2.49998e-1 |
| b31 | 8.85701e-1 | 8.55305e-1 | 2.91666e-1 |
| b32 | 2.30198e-2 | 2.04354e-2 | -1.24999e-1 |
| d21 | -3.46865e-1 | 3.52118e-1 | -4.99999e-1 |
| d31 | 1.90439e-2 | 1.88290e-2 | 2.53854e-7 |
| d32 | 3.98095e-1 | 3.94028e-1 | 3.74999e-1 |

Also quantitative in prose: `A_x minimum` for the bifurcation (`A_z=0`
boundary, Eq. 22) is "about 14% of the normalized distance r1" for Sun-Earth
L1/L2, i.e. **~200,000 km minimum halo x-amplitude** (p. 247). The
third-order solution's fidelity was cross-checked (p. 248-249) against a
differential-corrections numerical solution for L1 halos at `A_z = 110,000
km`: **maximum state-variable discrepancy < 3%**, "consistent with the order
of magnitude of the truncation error in the third-order development" — this
is Richardson's OWN reported accuracy bound, useful as a sourced tolerance
if this construction is ever adopted as a seed generator.

## 5. Figures (qualitative, not digitized — see fidelity guard)

Figs. 2-4 (pp. 248-250) plot Class I / Class II L1, L2, L3 orbit projections
(x-y, x-z, y-z planes) at `A_z = 1.25e5 km` (L1/L2) and `A_z = 1.25e7 km`
(L3, note the two-order-of-magnitude-larger amplitude — L3 halos are much
larger for the Sun-Earth system). These confirm the qualitative "north/south"
mirror-branch halo shape familiar from every subsequent halo paper; no
numeric values were extracted from them (Table I already gives the exact
coefficients that generate these curves analytically — reading the plots
would be redundant digitization of data already available in closed form).

## 6. Relevance to the cyclerfinder codebase (context for the Fable review)

**UPDATE (2026-07-12, #580): this gap is now filled.**
`search/cr3bp_seed_generator.py::richardson_halo_coefficients` /
`richardson_halo_ic` / `richardson_halo_seed` implement exactly the
construction described below, golden-tested against this digest's (now
corrected, see §4) Table I. The framing below is left as historical context
for why the task was scoped this way; see `data/OUTSTANDING.md`'s #580
resolution note for the implementation summary.

- The codebase already has an extensive CR3BP halo/vertical-Lyapunov capability:
  `search/cr3bp_seed_generator.py::lyapunov_seed_3d` (numeric seed generation)
  and `genome/known_corpus_3d.py` (Howell 1984 + Folta-Bosanac-Guzzetti-Howell
  2015 KNOWN_CORPUS anchors for the literature-novelty gate on any halo/NRHO
  rediscovery). No existing code path implements Richardson's analytic
  third-order construction, and no `KNOWN_CORPUS` anchor currently cites this
  paper specifically (only Howell 1984 and the 2015 cislunar reference catalog
  are anchored for the `{"halo"}`/`{"halo","nrho"}` topology labels).
- Richardson's construction is a **different capability axis** from what's
  there: a *closed-form, iteration-free* seed (μ, A_z) → (x0,y0,z0,vx0,vy0,vz0)
  map, vs. the existing numeric linearized-eigenvector + differential-
  correction approach. Whether this is worth building is a genuine judgment
  call (faster/more-robust seeding vs. redundant with an already-working
  numeric path) — this is exactly the question posed to the Fable review.
- The paper is Sun-Earth-specific in its published Table I numerics, but the
  closed-form coefficients (Appendix I) are valid for **any** collinear point
  of **any** mass ratio μ given the corresponding `c2,c3,c4` — i.e. it is a
  general-mu tool, same generality class as Howell 1984, just analytic
  instead of numeric.

## 7. References cited (p. 253)

- Broucke, R. (1968, 1969), JPL Reports 32-1168, 32-1360.
- Dasenbrock, R. (1973), U.S. Naval Research Lab Report 7564U — the algebraic
  manipulation program used to automate the successive-approximation algebra.
- Deprit, A. (1965), Icarus 4, 242.
- Farquhar, R. & Kamel, A. (1973), Celestial Mechanics 7, 458 — the
  Lindstedt-Poincaré quasi-periodic L2 analytic precedent this paper extends
  to full third-order periodicity.
- Farquhar, R., Muhonen, D., Newman, C., Heuberger, H. (1979), AAS/AIAA Paper
  79-126 — ISEE-3 mission design context motivating this analytic method.
- Moulton, F. (1920), *Periodic Orbits*, Carnegie Institute Publication 161.
- Pedersen, P. (1933, 1935), Monthly Notices Roy. Astron. Soc. 94:167, 95:482.
- Plummer, H. (1901, 1903a, 1903b), Monthly Notices Roy. Astron. Soc. 62:6,
  63:436, 64:98.
- Richardson, D. (1979), AAS/AIAA Paper 79-127 — "Halo-Orbit Formulation for
  the ISEE-3 Mission" (also cited by Howell 1984's own reference list — see
  [[2026-06-25-digest-howell-1984-halo-orbits]]).
- Richardson, D. (1980), "A Note on a Lagrangian Formulation for Motion about
  the Collinear Points," Celestial Mechanics — the companion paper flagged
  in §1 above, not yet in the corpus.
- Szebehely, V. (1967), *Theory of Orbits* — already in the corpus, digested
  2026-06-19.

End of digest.
