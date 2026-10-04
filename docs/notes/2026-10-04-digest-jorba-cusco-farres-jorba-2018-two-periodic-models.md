# Digest: Jorba-Cusco, Farres & Jorba (2018), "Two Periodic Models for the Earth-Moon System"

Frontiers in Applied Mathematics and Statistics 4:32 (2018), DOI 10.3389/fams.2018.00032, open access (CC BY), 14 pages.
Filed in the private paper corpus as
jorba-cusco-farres-jorba-2018-two-periodic-models-earth-moon-system-bcp-qbcp-frontiers-ams-4-32-doi-10.3389-fams.2018.00032.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with page and section or equation) or INFERRED
(our reading, or a comparison with project code). Page numbers are the journal's printed page numbers (1 to 14).
Context: tasks #891 and #892 corrected the project's two Sun-Earth-Moon models (`core/bcr4bp.py`, `core/qbcp.py`); this
paper is the source of the coefficient tables in `core/qbcp.py`.

## 0. What the paper is

READ (abstract, p1): two alternatives to the Earth-Moon RTBP for a massless particle, the Bicircular Problem (BCP) and
the Quasi-Bicircular Problem (QBCP), both periodically time dependent. The paper concludes (p1, p12) that the BCP is more
adequate near the triangular points and the QBCP near the collinear points. Methods: stroboscopic map, fixed points as
periodic orbits with the period of the Sun, parameterisation method for invariant manifolds (section 3.2), parallel
shooting. Computations used a Taylor integrator with demanded accuracy 1e-16 (section 7, p12).

## 1. The bicircular model as printed (section 4, p5, eq. 3)

READ. Units are those of the Earth-Moon RTBP (section 4, p5: "It is usual to take the units and the synodic coordinates
of the Earth-Moon RTBP"). RTBP Hamiltonian (eq. 1, p2):
`H_RTBP = 1/2 (px^2+py^2+pz^2) - x py + y px - (1-mu)/r_PE - mu/r_PM`, with `r_PE^2 = (x-mu)^2 + y^2 + z^2` and
`r_PM^2 = (x-mu+1)^2 + y^2 + z^2` (p2). So Earth (mass 1-mu) is at x = mu and Moon (mass mu) at x = mu-1.

Eq. 3 (p5), transcribed:

    H = 1/2 (px^2+py^2+pz^2) - x py + y px - (1-mu)/r_PE - mu/r_PM
        - (m_S/a_S^2) (y sin(theta) - x cos(theta)) - m_S/r_PS

Definitions, quoted from p5: "Here mu, r_PE and r_PM denote the same quantities as in (1). Moreover, m_S denotes the mass
of Sun, a_S the averaged semi-major axis of Sun, theta = omega_S t, omega_S is the frequency of Sun in this system of
reference, T_S = 2 pi / omega_S is its period and finally, r_PS^2 = (x - a_S cos theta)^2 + (y - a_S sin theta)^2 + z^2."

Sign by sign, what the page shows:

| Item | As printed |
| --- | --- |
| indirect term | `- (m_S/a_S^2)(y sin(theta) - x cos(theta))`, equivalently `+ (m_S/a_S^2)(x cos(theta) - y sin(theta))` |
| r_PS^2 | `(x - a_S cos(theta))^2 + (y - a_S sin(theta))^2 + z^2` |

The project's belief is stated correctly: the printed r_PS has `(y - a_S sin theta)`, which places the Sun at
(a_S cos theta, +a_S sin theta), while the indirect term as printed is the linearisation of `-m_S/r_PS` for a Sun at
(a_S cos theta, -a_S sin theta). INFERRED, not resolved here: expanding `-m_S/r_PS` about the origin for a Sun at
(a_S cos theta, s a_S sin theta) gives the cancelling term `+(m_S/a_S^2)(x cos theta + s y sin theta)`; the printed
indirect term matches s = -1, the printed r_PS matches s = +1.

Other statements on the page and its neighbours that bear on the sense of the Sun's motion (READ):
- Figure 2 (p5) draws the Sun S in the lower right (x > 0, y < 0) with theta marked by an arrow running clockwise from the
  +x axis towards S. INFERRED: this agrees with a Sun at (a_S cos theta, -a_S sin theta).
- Taylor expansion of the Sun's potential (p6), printed: `(1/a_S) (1 + (x cos theta - y sin theta)/a_S)`. This again
  carries `x cos theta - y sin theta`, the Sun at (a_S cos theta, -a_S sin theta). INFERRED.
- Section 4.1 (p7): "We have named T(x, y, theta) = -x cos theta + y sin theta" and the order-2 term
  `H_S^2(theta,x,y,z) = (1/a_S^3) ( (3/2) T(x,y,theta)^2 - (1/2)(x^2+y^2) )`, invariant under "(x, y, theta) ->
  (x, -y, -theta)" (the paper calls this the symmetry of the Coriolis-cancelled term). The order-3 term H_S^3 (p7) is
  printed in terms of rho^2 = x^2 + y^2 and T, with a polynomial in T that is not even; it is not transcribed here. INFERRED:
  T as named is minus the combination x cos theta - y sin theta of the p6 expansion, which is only an overall sign
  convention for T (T enters H_S^2 squared); it does not by itself fix the Sun's sense.
- The truncated model at linear order is `H_BCP^{<2} = H_RTBP - m_S/a_S` (p6): "the Coriolis term and the truncated Sun's
  potential cancel out and the dynamics is the one of the RTBP"; the BCP is a perturbation "with size O(m_S/a_S^3) ~ 0.0056".
- The BCP is described (section 1, p2) as "a restricted four body problem" with Earth and Moon "along a circular orbit
  around their common center of masses" and Sun and C_EM "in another circular orbit around C_SEM"; "the motion assumed for
  the primaries does not verify Newton's laws" (p2).

Constants (Table 3, p10, READ):

| mu | a_S | m_S | omega_S |
| --- | --- | --- | --- |
| 0.012150581623433623 | 388.81114302335106 | 328900.54999999906 | 0.92519598551829646 |

Table 3 serves both models (section 7: "Table 3 and Table 4 contains the values of the parametres used"). INFERRED: the
comment above `_QBCP_MU_EM` in `core/qbcp.py` line 185 attributes Table 3 to "Gimeno-Jorba (2018)"; the authors of this
paper are Jorba-Cusco, Farres and Jorba.

## 2. The quasi-bicircular model as printed (section 5, p9 to p10, eqs. 4 and 5)

READ. Frame (section 5, p9): the QBCP solution is planar and computed in the Jacobi frame; to study the Earth-Moon
vicinity "one has to perform three different transformations. First, one has to use a translation to move the origin from
the global barycenter to Earth's and Moon's center of masses. Second, one has to use a rotating (synodic) frame to keep
Earth and Moon fixed on the horizontal axis. Third, the unit of length is scaled so the distance between Earth and Moon is
equal to one." Earth and Moon sit on the x axis; the distances in eq. 4 are `r_pe^2 = (x-mu)^2 + y^2 + z^2`,
`r_pm^2 = (x-mu+1)^2 + y^2 + z^2`, `r_ps^2 = (x-alpha_7)^2 + (y-alpha_8)^2 + z^2` (p10). So Earth is at x = mu and Moon
at x = mu-1, as in the RTBP and the BCP of the same paper.

Eq. 4 (p9), transcribed:

    H = 1/2 alpha_1 (px^2+py^2+pz^2) + alpha_2 (px x + py y + pz z) + alpha_3 (px y - py x)
        + alpha_4 x + alpha_5 y - alpha_6 ( (1-mu)/r_pe + mu/r_pm + m_S/r_ps )

alpha_6 multiplies all three potential terms (Earth, Moon, Sun). No other alpha multiplies a potential term.

Eq. 5 (p10): `alpha_i(theta) = a_0^i + sum_{k>=0} a_k^i cos(k theta) + sum_{k>=0} b_k^i sin(k theta)`, with
"theta = omega_S t" and "omega_S the frequency of Sun", for i = 1..8, alpha_i: T -> R.

The parity sentence (p10), verbatim: "Moreover, alpha_i is odd for i = 1, 3, 4, 6, 7 and even for i = 2, 5, 8." INFERRED:
by the paper's own symmetry and by Table 4 this wording uses odd and even the other way round from the usual meaning for
functions of theta: i = 1, 3, 4, 6, 7 are the cosine series (even functions of theta) and i = 2, 5, 8 the sine series
(odd functions). Andreu (1998), as quoted in the module docstring of `core/qbcp.py`, says "The functions alpha_1, alpha_3,
alpha_4, alpha_6, alpha_7 are even. The other ones, alpha_2, alpha_5, alpha_8 are odd."

Symmetry (p10), verbatim: "it is easy to see that the Hamiltonian function (4) has the symmetry
(theta, x, y, z, xdot, ydot, zdot) -> (-theta, x, -y, z, -xdot, ydot, -zdot)." READ.

Meaning of the alpha_i (p10, READ, paraphrase with quoted fragments):
1. "(alpha_7, alpha_8, 0) is the position of Sun in the plane of motion of the primaries."
2. "alpha_1, alpha_2, alpha_3 and alpha_6 capture the fact that the distance between Earth and Moon is not constant."
3. "alpha_4 and alpha_5 take into account the Coriolis effect due to the rotating frame of reference."

Source of the values (p10): "one can only have a numerical approximation of these functions. In this case, we take
advantage on the computations done in [15] and take the same values for the Fourier coefficients of the periodic functions
alpha_i's." READ. The paper therefore does not compute the coefficients; they are Andreu's.

## 3. Table 4 (p11): headers, k = 0 entries, comparison with `core/qbcp.py`

READ. Caption: "Coefficients of the functions alpha_j, j = 1, ..., 8, in (5)." Footnote: "Due to the symmetries of the
model, each alpha_j only contains either sin or cos terms, so we only list either the a_k or b_k coefficients."

Column headers exactly as printed (the text layer and the rendered page agree):

| alpha_1 | alpha_2 | alpha_3 | alpha_4 | alpha_5 | alpha_6 | alpha_7 | alpha_8 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| a_k | a_k | b_k | a_k | b_k | a_k | a_k | b_k |

k = 0 entry of each column, as printed:

| alpha | k = 0 | k range |
| --- | --- | --- |
| 1 | 1.001841608924835e+00 | 0 to 12 |
| 2 | 0.e0 | 0 to 12 |
| 3 | 9.999999999999983e-01 | 0 to 12 |
| 4 | -9.755242327484885e-04 | 0 to 11 |
| 5 | 0.e0 | 0 to 11 |
| 6 | 1.000907457708158e+00 | 0 to 12 |
| 7 | -6.314069568006227e-02 | 0 to 13 |
| 8 | 0.e0 | 0 to 13 |

Observations (INFERRED):
- The headers label alpha_2 as a_k (cosine) and alpha_3 as b_k (sine). The paper's own k = 0 entries contradict that: a
  cosine series has a non-zero constant a_0, a sine series has none, and the table prints 0 for alpha_2 and 0.99999... for
  alpha_3. The parity sentence above, the symmetry of eq. 4 and Andreu's table all make alpha_2 a sine series and
  alpha_3 a cosine series. So the two headers look transposed in print. The other six headers agree with the parity
  statement read in the usual sense.
- alpha_1 at k = 5 is printed `-38.068581391005552e-08` (the mantissa has two digits before the point, unlike every other
  entry in the table). Read literally this is -3.8068581391005552e-07. Andreu's table (the cited source, text layer of the
  filed thesis: "8:06858139100555e 08" in the k = 5 row) has -8.06858139100555e-08; the neighbouring entries (k = 4:
  1.176e-04, k = 6: 9.843e-07) fit -8.07e-08, not -3.81e-07. The two numbers share the digits 06858139100555 and differ in
  the leading digits, consistent with a typesetting slip in the leading digits. This note does not decide which is right
  beyond that observation; `core/qbcp.py` now carries Andreu's value (#892).

Comparison with `_COEFFS_ALPHA1` to `_COEFFS_ALPHA8` in `src/cyclerfinder/core/qbcp.py` (done by a script that parsed the
text layer of Table 4 and compared against the lists with the stated sign flips: alpha_4, alpha_5, alpha_7, alpha_8
multiplied by -1; tolerance 1e-9 relative):
- Magnitudes: every printed entry (k ranges above, 8 columns) matches the code in magnitude except alpha_1 at k = 5 (above).
  No other magnitude differs. READ plus INFERRED (script comparison).
- Signs, INFERRED, not a magnitude difference: the code does not apply the reflection uniformly at the tail. The entries
  alpha_7 at k = 12 and 13, and alpha_8 at k = 11, 12 and 13, are in the code with the same sign as printed, whereas the
  stated rule (multiply alpha_7 and alpha_8 by -1) would flip them. All other entries of alpha_4, alpha_5, alpha_7, alpha_8
  carry the flip. The affected terms are at most 1.9e-8 in absolute coefficient (alpha_8 k = 11) and 1.6e-10 at k = 13, against
  leading coefficients of order 4e2, so the effect is far below the other terms; worth fixing for consistency but unlikely to be visible. I did not check whether the code's listing of alpha_7 and alpha_8 at
  these k comes from the printed table or from Andreu's.
- Length: the paper prints alpha_4 and alpha_5 to k = 11 (12 entries), alpha_6 to k = 12, alpha_7 and alpha_8 to k = 13; the
  code lists have the same lengths.

## 4. Printed periodic orbits and candidates for positive controls

Nothing in the paper prints a full initial state (x, y, z, px, py, pz at a stated phase theta) for any orbit. The data
available are as follows.

### 4.1 BCP, triangular points (section 4.1, 4.2, Figures 3 to 5)

READ (p7): "each triangular point is replaced by three periodic orbits with the same period as Sun. One small and unstable
(the actual replacement of L4) and two which are stable." Named PO1 (saddle x center x center, the replacement of L4),
PO2 and PO3 (totally elliptic) (p7). By symmetry the dynamics near L5 is the same (p7). Figure 3 left (p6) is a
continuation diagram in the Sun mass parameter epsilon (vertical axis, 0 to 1) against x (horizontal axis, -0.8 to 0); the
orbits sit near x between -0.8 and 0 with y near 0.8 to 0.95 (Figure 4, p6). Figure 4 shows the stroboscopic map near the
triangular points; the three fixed points are marked as crosses, PO1 near (x, y) = (-0.5, 0.87) (read off the plot,
INFERRED). No numerical coordinates and no eigenvalues are printed for PO1, PO2, PO3.
Reproducibility: label only. Not reproducible from the page. The continuation from the RTBP (epsilon from 0 to 1) is a
procedure the project could run, and the figures give coarse targets (about 2 significant digits).

Vertical families VF1, VF2, VF3 of 2D tori (Figure 3 right, p6; section 4.3, p8): printed as curves of frequency against
pz. Figure 5 (p7): effective stability regions for tori of VF3 at pz = 0.5 and pz = 0.8, plotted in (alpha, r). Setup
printed in section 4.3 (p8): polar grid with h_r = 0.001 and h_alpha = 0.0002, 15000 Moon revolutions, termination on
collision or y < 0. This is a procedure with parameters, not a printed state; reproducing the figure is a long
computation.

### 4.2 BCP, L2 weakness (section 4.4, Figure 6, p8 to p9)

READ: the continuation of the L2 point to the BCP "reaches a turning point and it never reaches the homotopy level of the
BCP"; "the translunar dynamical structure is lost in the BCP." Figure 6 left: epsilon against x over x from -1.19 to -1.11
with L2 marked near x = -1.15. Figure 6 right: a large periodic orbit near L2 in the BCP (planar, x from about -1.19 to
-1.11, y from -0.15 to 0.1). Not reproducible from the page; useful only as a qualitative negative control (the BCP has no
translunar equivalent).

### 4.3 QBCP, dynamical equivalents of L1, L2, L3 (section 5.1, Figure 7, Table 1, p9)

READ (p9 to p10): "the orbits replacing L1 and L2 are small, their maximal distance to the corresponding equilibrium point
is of order O(10^-6)"; "linear normal behavior ... saddle x center x center"; "the unstable direction of L1 (of order 10^8)
and of L2 (of order 10^6) are large and this implies huge propagation of error near these orbits"; "the dynamical
equivalent of L3 has a very weak unstable direction" (p10 to p11). Figure 7 axis ranges give the plotted positions
(INFERRED, read off axis labels): L1 orbit at x about -0.836915 (axis ticks -8.36916e-01 to -8.36914e-01, y within about
3e-6), L2 orbit at x about -1.155682 (ticks -1.155683 to -1.155682 e0, y within about 4e-6), L3 orbit at x about 1.0050
(ticks 1.0048 to 1.0054, y within about 1e-3). The L1 and L2 numbers are the three-body points to the printed digits.

Table 1 (p9), eigenvalues of the stroboscopic map (one period of the Sun), "We only put three for each orbit. The rest are
given by their inverses due to the symplectic character of the stroboscopic map.":

| index | L1 real | L1 imag | L2 real | L2 imag | L3 real | L3 imag |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 460182151.57 | 0 | 2397196.84 | 0 | 3.370855 | 0 |
| 2 | -0.987151 | 0.159784 | 0.995818 | 0.0913562 | 0.863840 | -0.503764 |
| 3 | -0.963639 | 0.267205 | 0.917527 | 0.3976716 | 0.841148 | 0.5408042 |

Reproducibility: eigenvalues are printed to 7 or more significant digits for the three orbits, so they are a quantitative
check, but no state is printed; the orbit has to be found by the project's own corrector (a T_S-periodic fixed point near
the three-body L1, L2, L3). INFERRED: the L1 and L2 unstable multipliers (4.6e8 and 2.4e6) are very sensitive to the
integrator tolerance and to the model; matching them to a few digits is a strong test, matching them to the order of
magnitude is a weak one. The paper does not state the starting phase of the Poincare section (theta = 0 is implied by
"the section is temporal", section 3.2, p4; INFERRED).

### 4.4 QBCP, resonant orbits (section 5.2, Table 2, Figures 8 and 9, p9 to p12)

READ: Table 2 (p9) lists continuations of low order resonant orbits from the RTBP to the QBCP, taken originally from
Andreu [15]. Columns: RTBP label (012, 014, 018, 01C, 01E, 022, 026, 02A, 02E, 026), resonance (1:2, 1:1, 1:1, 1:3, 1:3,
1:2, 1:6, 1:2, 1:3, 1:4), number of bifurcating orbits (2 or 4), and the QBCP labels (12, 13; 14 to 17; 18, 19, 1A, 1B;
1C, 1D; 1E, 1F; 22 to 25; 26 to 29; 2A to 2D; 2E, 2F; 2G, 2H). Colour code (p11): blue = saddle x center x center, green
= saddle x saddle x center, cyan = totally hyperbolic, yellow = totally elliptic, red = continuation does not reach the
homotopy level of the QBCP. The label's second digit refers to the Lyapunov or Halo family of L1 or L2 (p11). Figure 9
(p10): the resonant orbit 2G in the (y, z) projection (about 0.04 in x extent, z up to 0.2) and its manifolds of order 64.
Only labels, resonance orders and linear stability type are printed. Not reproducible from the page: no state, no period
value (periods are the Sun period times the resonance), no eigenvalues. Andreu's thesis is the source for the underlying
numbers (see the existing digest of that thesis).

### 4.5 Order-64 manifold approximations (section 5.3, Figures 8 and 9)

READ (p11 to p12): manifolds of the L1, L2, L3 orbits and of 2G computed to order 64 by the parameterisation method; the
L3 manifold "passes very close to the triangular points"; L1 used 128-bit extended precision with single shooting, L2
multiple shooting with two sections. Graphical only.

## 5. What the paper cites (reference list, pp12 to 14)

READ (numbers as in the paper):
- QBCP coefficients: [15] M. A. Andreu, "The Quasi-Bicircular Problem", PhD thesis, University of Barcelona (1998). The
  authors "take the same values for the Fourier coefficients" from it (p10), and Table 2 is "originally in Andreu [15]"
  (p11). Related: [16] Andreu and Simo, "The quasi-bicircular problem for the Earth-Moon-Sun parameters" (2000),
  online at http://www.maia.ub.es/dsg/2000/index.html; [17] M. A. Andreu, "Dynamics in the center manifold around L2
  in the quasi-bicircular problem", Celestial Mech. 84:105-33 (2002), doi 10.1023/A:1019979414586; [18] Bihan, Masdemont,
  Gomez, Lizy-Destrez, "Invariant manifolds of a non-autonomous quasi-bicircular problem computed via the
  parameterization method", Nonlinearity 30:3040 (2017), doi 10.1088/1361-6544/aa7737. The QBCP "was introduced by
  C. Simo" (p2, no reference attached to that sentence). Continuation-based alternatives for the quasi-bicircular
  solution (p2): [20] Gabern, PhD thesis, Barcelona (2003); [21] Gabern and Jorba, Discrete Contin. Dyn. Syst. B 1:143-82
  (2001); [22] Gabern, Jorba and Robutel, Discrete Contin. Dyn. Syst. B 4:843-54 (2004).
- Bicircular model: "The BCP is a restricted four body problem [9, 10]": [9] Huang, "Very Restricted Four-Body Problem",
  NASA TN D-501 (1960); [10] Cronin, Richards and Russell, "Some periodic solutions of a four-body problem", Icarus
  3:423-8 (1964). Utilised in other cases [11] Barrabes, Gomez, Mondelo, Olle, "Pseudo-heteroclinic connections between
  bicircular restricted four-body problems", MNRAS 462:740-50 (2016). Derivation of the equations of motion: [12] Gomez,
  Jorba, Masdemont, Simo, "Study Refinement of Semi-Analytical Halo Orbit Theory", ESOC contract 8625/89/D/MD(SC), final
  report (1991). The replacement of each triangular point by three periodic orbits and the use of the BCP near the
  triangular points: [51] Simo, Gomez, Jorba, Masdemont, "The Bicircular model near the triangular libration points of
  the RTBP", in From Newton to Chaos, Plenum (1995), pp343-70. Vertical families in the BCP: [13] Castella and Jorba,
  Celestial Mech. 76:35-54 (2000). Regions of stable motion in the real Earth-Moon system: [14] Jorba, Astron. Astrophys.
  364:327-38 (2000).
- Other earlier work on coherence of solutions close to bicircular: [19] Siegel and Moser, Lectures on Celestial Mechanics
  (1971).

## 6. What this means for the project

Statements are INFERRED unless marked.

Can be used as tests:
- The printed constants (Table 3) and the full Table 4: they define the QBCP the project implements. The magnitude check in
  section 3 above is complete: only alpha_1 at k = 5 differs, and the sign inconsistencies of the tail entries are listed.
- Table 1 eigenvalues (7 or more digits) for the T_S-periodic orbits replacing L1, L2, L3 in the QBCP: a quantitative
  positive control for a corrected `core/qbcp.py` once the orbits are located by the project's own corrector. It tests the
  model (alpha_6 multiplying the whole potential, the table orientation) and the corrector. Note the L1 and L2 unstable
  multipliers are huge and tolerance-sensitive; treat the leading digits with care.
- The size statement (maximal distance from the three-body point of order 1e-6 for L1 and L2) and the stability type
  (saddle x center x center for all three, with the L3 unstable direction weak): a cheap qualitative check. The module
  docstring of `core/qbcp.py` already uses the 3e-6 statement.
- Qualitative BCP checks: three T_S-periodic orbits replace L4 (one saddle x center x center, two totally elliptic), and
  the continuation in epsilon (Sun mass multiplier) of the L2 equivalent turns back before epsilon = 1, so the BCP has no
  translunar equivalent. The project's `core/bcr4bp.py` should show both. These are controls on topology and stability type,
  not on numbers.
- The BCP sign question: the page contains an internal difference between the printed r_PS (y - a_S sin theta) and the
  indirect term, the Taylor expansion and Figure 2 (all consistent with a Sun at (a_S cos theta, -a_S sin theta)).
  Anything that decides between them has to be computed, not read: a candidate check is that the L4 equivalent of the
  corrected model has the stability type and the epsilon-continuation of Figures 3 and 4 for the sign adopted, since the
  indirect term and r_PS must have the same sense for the Coriolis cancellation of section 4 (p6) to hold.

Cannot be used as tests:
- PO1, PO2, PO3, the vertical tori VF1 to VF3 and the Figure 5 stability regions: no state or eigenvalue printed, only
  figures and a procedure.
- Table 2 resonant orbits (including 2G): labels and stability type only; the numbers are in Andreu's thesis.
- Any quantity needing an initial state: none is printed anywhere in the paper.

Contradictions with project docstrings found while reading (INFERRED, factual):
1. `core/qbcp.py` line 185 labels Table 3 "Gimeno-Jorba (2018)"; the paper's authors are Jorba-Cusco, Farres and Jorba.
2. The module docstring quotes Andreu's parity statement ("even" for alpha_1, 3, 4, 6, 7). The 2018 paper's own sentence
   says "odd" for the same set (usual meaning reversed); the docstring's claim that the cosine/sine assignment follows
   "Andreu 1998" is right, and it correctly avoids relying on the 2018 wording.
3. The comment block above `_COEFFS_ALPHA1` states that alpha_4, alpha_5, alpha_7, alpha_8 are reflected (multiplied by
   -1); the code does so for all entries except alpha_7 k = 12, 13 and alpha_8 k = 11, 12, 13, which carry the printed sign.
4. Nothing else in the docstrings of `core/qbcp.py` or the quoted description of `core/bcr4bp.py` conflicts with the page.
   In particular the docstring's statement that alpha_6 multiplies the whole Newtonian potential, and that the 2018 table
   prints alpha_1 at k = 5 as `-38.068581391005552e-08` and labels alpha_2 as a_k and alpha_3 as b_k, are all confirmed
   by the page.
