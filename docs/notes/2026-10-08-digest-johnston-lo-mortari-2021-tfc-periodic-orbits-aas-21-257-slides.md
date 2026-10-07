# Digest: Johnston, Lo & Mortari 2021, "A Functional Interpolation Approach to Compute Periodic Orbits in the CR3BP" (conference slides, AAS 21-257)

H. Johnston (Texas A&M), M. W. Lo (JPL) & D. Mortari (Texas A&M), 31st AAS/AIAA Space Flight Mechanics Meeting,
virtual, 1-4 Feb 2021. Every slide footer prints **AAS 21-257**. **This file is the slide deck (17 slides), not the
paper.** Slide 1 carries a "(c) 2019" line. The matching full paper is J. Mathematics 9(11):1210 (27 May 2021),
**doi 10.3390/math9111210** (Crossref: title, the same three authors, container "Mathematics", issued 2021-05-27). I
could not open the MDPI page (HTTP 403), so I have **not read the journal paper** and say nothing about its content.

- Source file: upload `44c3f900-CL21_0389.pdf` (JPL clearance CL#21-0389), 17 pp (slides, 960 x 540 pt), md5
  59c5299c923a03ea2b34901d6b819335. The text layer is poor (pdftotext gives about 10.7 kB with the slide text mixed up).
- OCR copy made in my folder (`ocrmypdf --force-ocr -l eng`, PDF/A, 17 pp): `johnston-ocr.pdf`, **md5
  1328194884a33ff33d04087f8dd340cf**. This is the copy to file. Render check, pages 3, 12 and 15 at 100 dpi, original
  versus copy: same pixel size (1334 x 750), mean absolute grey difference 8.7, 12.0 and 9.4 of 255 (page rasterised
  at 400 dpi; font anti-aliasing). Page 12 of the original viewed in full: plots, tables and arrows all present.
- Proposed filename: `johnston-lo-mortari-2021-functional-interpolation-periodic-orbits-cr3bp-slides-AAS-21-257.pdf`
- How I read it: every slide viewed as an image (slides 5, 6, 8-16 at 80 dpi, slide 12 plot zoomed at 250 dpi);
  numbers quoted below come from those images. No table of orbit data exists (section 2).

## 0. Verdict

A 17-slide talk, 10 of them content. It shows that the **theory of functional connections (TFC)** can replace the
classical differential corrector for CR3BP periodic orbits, on **Earth-Moon L1 and L2 planar Lyapunov orbits** and
**northern-bifurcation halo orbits**, with residual and run-time scatter plots. The slides state the claim in one line
(slide 16): "comparable results in terms of speed and accuracy to the differential corrector method" and "easier
implementation through the TFC toolbox".

What it gives the project:

1. **For #948 (alternative corrector):** the formulation and the solver settings (section 1) are enough to understand
   the method, but **not enough to reproduce it**: the basis function type is not on any slide.
2. **For #970's control as a test case:** **no.** The deck prints **no orbit initial conditions, no periods, no
   per-orbit Jacobi constants** (only a Jacobi-constant axis from 2.90 to 3.20 on the plots). There is nothing to
   transcribe and nothing to check by integration. The brief asked for a transcription and a DOP853 check; I
   report instead that the source does not contain the data. The journal paper (doi above) is the place to look; it
   is on the wanted list proposal below.
3. **An accuracy and speed claim to test,** not to trust: TFC residuals reach about 1e-15 to 1e-12 against about 1e-14
   to 1e-12 for the differential corrector (L1 Lyapunov), at 0.25 to 1.7 s per orbit (slide 12). The deck compares
   against its own differential corrector, in Python (numpy, JAX), so the speed claim is not portable.

PROPOSALS only: (a) file the OCR copy as a slides-only item; (b) add the Mathematics paper to the wanted list (row 69
below), read it before any #948 design decision; (c) do not cite the slides for orbit data.

## 1. Formulation (slides 3-11)

- TFC (slides 4-6): "functional interpolation" in additive form. A constrained expression `y(t, g(t))` satisfies the
  constraints for any free function g(t); the ODE then becomes an unconstrained problem in g. Slide 5 example:
  `y' = 2y, y(0) = 3` gives `y(t, g) = g(t) - g(0) + 3` (this satisfies y(0) = 3 for any g; the exact answer is
  3 e^(2t)). Method (slide 6): define the free function g, discretise the domain, solve the algebraic system by least
  squares.
- Dynamics (slides 8-9): barycentric rotating frame, m1 at {-mu, 0, 0}, m2 at {1-mu, 0, 0}. **Slide 8 prints
  `mu = m1/(m1 + m2)`; with m1 at -mu the standard definition is `m2/(m1 + m2)`** (slide typo; the numbers on slides
  12-15 are consistent with m2/(m1+m2); see section 3). Slide 9: `Omega = (x^2 + y^2)/2 + (1-mu)/R1 + mu/R2 +
  (1-mu) mu / 2` and `C = (x^2 + y^2) + 2(1-mu)/R1 + 2 mu/R2 + (1-mu) mu - (xdot + ydot + zdot)`. **The velocity term
  is printed without squares** (should be xdot^2 + ydot^2 + zdot^2). Note the `+ (1-mu) mu` constant: this Jacobi
  constant is 0.012003 higher than the form used in, e.g., Blanchard et al. 2023 (`C = 2 Omega - v^2` with no
  constant). Periodicity: `r(0) = r(T)`, `v(0) = v(T)`.
- Time transformation (slide 11, for Lyapunov orbits): `dt/dzeta = R2`, mapping `t in [0,T] -> zeta in [0, zeta_f] ->
  tau in [-1, +1]`, so `rdot_i = (b^2 / R2) r_i'(tau)`. The slide's plot shows "scaled" versus "non-scaled"
  residuals: scaled stays at about 1e-14 to 1e-12 over the whole range (L1), non-scaled starts near 1e-2 at C = 2.90 and
  falls to 1e-14 only above C of about 3.15.
- Implementation (slide 10): public `tfc` Python toolbox on GitHub; JAX for automatic differentiation of Jacobians.

## 2. Numerical tests (slides 12-15); no orbit table

Earth-Moon: **m1 = 5.9724e24 kg, m2 = 7.346e22 kg** (slides 12, 14, 15), giving mu = 0.0121505 (ours, `checks.out`).
TFC settings as printed, in `johnston-lo-mortari-2021-tfc-settings.yaml`:

| test | slide | N points | m basis terms | epsilon | max iterations | least squares |
|---|---|---|---|---|---|---|
| Lyapunov L1 | 12 | 140 | 130 | 2.22e-16 | 20 | numpy.pinv() |
| Lyapunov L2 | 14 | 140 | 130 | 2.22e-16 | 20 | numpy.pinv() |
| Halo (northern bifurcation, L1 and L2) | 15 | 200 | 190 | 2.22e-16 | 20 | numpy.pinv() |

Table digits: these are slide-image readings, confirmed against the OCR text of the same slides (two witnesses:
image and OCR; both give 140/130, 200/190, "2.22 x 10^-16" and "20"; the OCR of slide 14 garbled 10^-16, read on the
image).

Computation-time breakdown (slide 13, L1 Lyapunov, as printed): loss function 5.69 ms (std 3.46), Jacobian 91.1 ms
(std 55.8), least squares 329.8 ms (std 203.8). Least squares is 78 percent of the sum of the three means (ours:
329.8 / 426.59). The slide-13 text layer gave 50.8 for the Jacobian std; the image shows 55.8; **image reading kept**.

Result statements, read on the plots only (no digits printed, so approximate):
- L1 Lyapunov: Jacobi constant 2.90 to 3.20 (about 30 orbits, 0.01 steps). TFC residuals fall from about 2e-12 to
  about 1e-15 as C rises; differential corrector scatters at 1e-14 to 1e-12. Time per orbit 0.25 to 1.7 s.
- L2 Lyapunov (slide 14): the **differential corrector diverged** for about nine orbits at the low-C end (C about 2.90
  to 2.98), and TFC hit its maximum-iteration limit (flagged on the slide) for about a dozen orbits at low C (about 2.91
  to 3.04). So the "comparable" claim fails at the low-C end of L2 for the differential corrector, and TFC pays in
  iterations there. Read from a 80 dpi image; counts are approximate.
- Halo (slide 15, "northern bifurcation"): C axis 3.025 to 3.20; residuals about 1e-14 to 1e-15 for most orbits, a few
  low-C orbits at 1e-10 to 1e-12; run times about 1.5 to 9 s (L1 and L2).

## 3. Checks (our arithmetic; `checks.py`, `checks.out`)

- With m2/(m1+m2) and the slide-9 Jacobi form including `+ mu (1-mu)`: **C(L1) = 3.200343, C(L2) = 3.184162**
  (collinear points by root finding of dOmega/dx, ours). The vertical dashed lines on the slide-12 residual plot
  (250 dpi zoom, gridline 3.00 at x = 470 px, 3.20 at about 1067 px) sit at about 3.2005 and 3.1847, so **the slides'
  "E(L1)" and "E(L2)" lines are at C(L1) and C(L2)** under this convention (agreement about 0.0002 to 0.0005, within
  my pixel reading). With the constant omitted the values would be 3.18834 and 3.17216, which do **not** match the
  lines. So the deck's Jacobi axis includes the `+ mu (1-mu)` constant. C(L3) (for context) is 3.012147 without the
  constant.
- Consequence: Lyapunov orbits plotted up to C about 3.20 start essentially at L1 (small amplitude); the low end
  C = 2.90 lies well below C(L3) + 0.012 = 3.0241, so those orbits are large and enclose the Moon (slide 12 shows big
  loops around the Moon).
- Not checked: no orbit can be integrated (no ICs). TFC accuracy claims are untested by us.

## 4. Citation mining

The deck has **no reference slide** ("Ref [6]" is cited on slide 11 for the time transformation, but no list is
printed). Nothing to mine. Related items:

| item | status |
|---|---|
| Johnston, Lo & Mortari 2021, Mathematics 9(11):1210, doi 10.3390/math9111210 (the full paper) | **not held, no row.** Proposed row 69, priority medium-high for #948 |
| Mortari 2017, "The theory of connections: connecting points", Mathematics 5(4):57 | not held; no row; proposed row 70 (low; the TFC foundation) |
| Johnston & Mortari 2018 / Johnston, Schiassi, Furfaro, Mortari 2020 TFC papers | not named on any slide; not requested |
| Other held corrector papers | Baresi-Olikara-Scheeres 2018, Olikara 2016 thesis are held; not cited by the deck |

Proposed wanted rows:
- 69 | Johnston, H., Lo, M. W. & Mortari, D. (2021), "A Functional Interpolation Approach to Compute Periodic Orbits in the Circular-Restricted Three-Body Problem", Mathematics 9(11):1210, doi 10.3390/math9111210 (open access) | the full paper behind the held AAS 21-257 slides; likely holds the orbit tables the slides lack | medium-high (#948)
- 70 | Mortari, D. (2017), "The Theory of Connections: Connecting Points", Mathematics 5(4):57 | TFC foundation | low

*Filed as `cyclers_pdf/papers/johnston-lo-mortari-2021-functional-interpolation-periodic-orbits-cr3bp-aas-21-257-slides-jpl-cl-21-0389.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*
