# Digest: Rosales, Jorba & Jorba-Cusco (2023), "Invariant manifolds near L1 and L2 in the quasi-bicircular problem"

Celestial Mechanics and Dynamical Astronomy 135:15 (2023), DOI 10.1007/s10569-023-10129-4, 47 pages (received 13 July 2022,
revised 1 December 2022, accepted 13 February 2023, published online 11 March 2023; open access).
Filed in the private paper corpus as
rosales-jorba-jorba-cusco-2023-invariant-manifolds-near-l1-l2-quasi-bicircular-problem-cmda-135-15-doi-10.1007-s10569-023-10129-4.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with printed page and section, equation, table or figure)
or INFERRED (our reading, a derivation from printed quantities, or a comparison with project code). Page numbers are the
journal's printed "Page n of 47". All numbers were taken from the PDF text layer and checked against the page images for the
tables that matter most (Tables 2 to 5, the normal-frequency table and Table 12); Tables 14 and 15 were extracted by script
from the text layer (count check: 77 and 68 coefficients, equal to the number of printed numbers on the pages).
Values read off a plotted figure are marked "graph read" and are good to a few units in the last plotted digit only.
Context: tasks #892 (quasi-bicircular module `core/qbcp.py`, corrected 2026-10-04) and #893 (model identity tests).

## 0. What the paper is

READ (abstract, p1): "The quasi-bicircular problem (QBCP) is a periodic time-dependent perturbation of the Earth-Moon
restricted three-body problem (RTBP) that accounts for the effect of the Sun. It is based on using a periodic solution of the
Earth-Moon-Sun three-body problem to write the equations of motion of the infinitesimal particle. The paper focuses on the
dynamics near the L1 and L2 points of the Earth-Moon system in the QBCP. By means of a periodic time-dependent reduction to
the center manifold, we show the existence of two families of quasi-periodic Lyapunov orbits around L1 (resp. L2) with two
basic frequencies. The first of these two families is contained in the Earth-Moon plane and undergoes an out-of-plane
(quasi-periodic) pitchfork bifurcation giving rise to a family of quasi-periodic Halo orbits. This analysis is complemented
with the continuation of families of 2D tori. In particular, the planar and vertical Lyapunov families are continued, and
their stability analyzed. Finally, examples of invariant manifolds associated with invariant 2D tori around the L2 that pass
close to the Earth are shown. This phenomenon is not observed in the RTBP and opens the room to direct transfers from the
Earth to the Earth-Moon L2 region."

READ: contents (p2): 1 Introduction (1.1 RTBP, 1.2 BCP, 1.3 QBCP, 1.3.1 dynamical substitutes of the collinear points);
2 Center manifold around L1 and L2 (2.1 testing the software, 2.2 L1, 2.3 L2); 3 Families of 2D invariant tori (3.1 L1,
3.2 L2); 4 Transfers in the QBCP; 5 Conclusions; Appendix A (coefficients of the center manifolds).

READ: the paper does NOT print the Fourier coefficients of the alpha functions. Table 3 caption (p8): "The values of the
Fourier coefficients for the functions alpha_i can be found in Andreu (1998) and are available upon request to the
corresponding author, J.J.R". It also prints no initial condition of any torus, invariant curve or periodic orbit other than
Table 4. Its numerical content that a dynamical module can test is: Table 3 (constants), Table 4 (the two substitute orbits at
t = 0), Table 5 (their monodromy eigenvalues), the gamma and normal-frequency tables (p10), and rotation numbers of tori given
in the text (section 4 below). Everything else is the centre-manifold coefficients (Tables 14 and 15) and figures.

## 1. The model exactly as printed

### 1.1 The RTBP (section 1.1, p4)

READ. Eq. (1): `H_RTBP = (1/2)(P_X^2 + P_Y^2 + P_Z^2) + Y P_X - X P_Y - (1 - mu)/R_PE - mu/R_PM`, with
`R_PE^2 = (X - mu)^2 + Y^2 + Z^2`, `R_PM^2 = (X - mu + 1)^2 + Y^2 + Z^2`, and the momenta "P_X = X' - Y, P_Y = Y' + X and
P_Z = Z'" (primes are time derivatives; the paper prints dots). Text: "the primary of mass 1 - mu (resp. mu) is at x = mu
(resp. x = mu - 1)". Table 1 (p4) lists mu for three systems: Sun-Earth 3.04042339E-6, Sun-Jupiter 9.54791915E-4,
Earth-Moon 1.21505816E-2.

### 1.2 The BCP (section 1.2, pp4-6)

READ. Eq. on p5: `H_BCP = (1/2)(P_X^2 + P_Y^2 + P_Z^2) + Y P_X - X P_Y - (1 - mu)/R_PE - mu/R_PM - m_S/R_PS -
(m_S/a_S^2)(Y sin(theta) - X cos(theta))`, with `R_PS^2 = (X - X_S)^2 + (Y - Y_S)^2 + Z^2`. Text (p5): "the parameters m_S and
a_S are the mass of the Sun and its distance to the Earth-Moon barycenter, respectively. The frequency of the Sun around the
Earth-Moon barycenter is omega_s, and theta = omega_s t, (X_S, Y_S) = (a_S cos(theta), -a_S sin(theta)) is the Sun position
vector". Fig. 1 (p5) shows the Moon M at negative x, L1 between E and M, L2 beyond M, and S at the lower right with theta
measured clockwise from +x. "Note that in this reference system the Sun moves around the origin in a circular motion."
H_BCP = H_RTBP + H_S with `H_S = -m_S/R_PS - (m_S/a_S^2)(Y sin(theta) - X cos(theta))` (p6), and the homotopy
`H^epsilon = H_RTBP + epsilon H_S` (eq. 2, p6).
READ: the paper says the BCP "is, therefore, not coherent in the sense that the motion of Earth, Moon and Sun does not verify
Newton's laws". "Remark 1 ... Notice from Table 2 that [omega_s] is a little bit smaller than 1; therefore, the period T_s of
the vectorfield (for both the BCP and QBCP) is about 30 days."
INFERRED: Sun at angle -theta, so it moves clockwise in this frame; the indirect term equals `+(m_S/a_S^2)(X cos(theta) -
Y sin(theta))`, which cancels the first-order term of -m_S/R_PS about the barycentre, so the printed BCP is internally
consistent (same finding as for the 2020 paper, see the digest of Jorba, Jorba-Cusco and Rosales 2020).

Table 2 (p5), "Parameters of the BCP", as printed:

| Symbol | Value |
| --- | --- |
| mu | 0.012150581623433 |
| m_s | 328900.5499999991 |
| omega_s | 0.925195985518289 |
| a_s | 388.8111430233511 |

### 1.3 The QBCP (section 1.3, p7, eq. 3, and p8)

READ. Construction (p7, quoted): "The first step is to compute a quasi-bicircular solution that models the motion of the Sun,
the Earth, and the Moon under each other's gravitational influence. This is accomplished by expressing the three-body problem
in the Jacobi formulation. Then, an approximation to the Jacobi decomposition of the three-body problem is obtained as Fourier
series, solving for the coefficients. The details are in Andreu (1998). With this solution, the origin of the (inertial)
reference frame is translated from the center of masses of the Sun, Earth, and Moon to the Earth-Moon barycenter. Then, the
reference frame is rotated such that the x-axis contains both the Earth and the Moon. A third change is a time-dependent
transformation that keeps the Earth and the Moon fixed on the x-axis. This defines a pulsating reference frame with period
equal to one revolution of the Earth and the Moon around their common barycenter. Also, the unit of distance is scaled such
that the distance between the Earth and the Moon is equal to one, the time is scaled such that one revolution of the pulsating
reference frame is equal to 2 pi, and the unit of mass is scaled such that m_E + m_M = 1".

READ. Eq. (3), p7, transcribed in full:

    H_QBCP = (1/2) alpha_1 (P_X^2 + P_Y^2 + P_Z^2) + alpha_2 (P_X X + P_Y Y + P_Z Z) + alpha_3 (P_X Y - P_Y X)
             + alpha_4 X + alpha_5 Y - alpha_6 ( (1 - mu)/R_PE - mu/R_PM - m_S/R_PS )

Checked on a 220 dpi crop of p7: the second line is printed exactly `+ alpha_4 X + alpha_5 Y - alpha_6 ( (1 - mu)/R_PE - mu/R_PM
- m_S/R_PS )`, that is a plus sign on the Earth term inside the bracket and minus signs on the Moon and Sun terms. As typeset
that is not the Newtonian potential (the Moon and the Sun would repel); eq. (1) and (2) of the same paper print all three with
minus signs (`- (1 - mu)/R_PE - mu/R_PM - m_S/R_PS`). Stated factually: the signs inside the bracket of eq. (3) are a typesetting
slip. INFERRED: the intended expression is `- alpha_6 [ (1 - mu)/R_PE + mu/R_PM + m_S/R_PS ]`, the form in the 2018 paper's eq. 4
and in the project's `core/qbcp.py`. What the equation fixes unambiguously is which terms alpha_6 multiplies: all three
gravitational terms, as a factor.

READ. Definitions under eq. (3): `R_PE^2 = (X - mu)^2 + Y^2 + Z^2`, `R_PM^2 = (X - mu + 1)^2 + Y^2 + Z^2`,
`R_PS^2 = (X - alpha_7)^2 + (Y - alpha_8)^2 + Z^2`. So the Sun's position in the pulsating frame is (alpha_7(theta),
alpha_8(theta), 0) and the Earth is at X = mu, the Moon at X = mu - 1, both fixed.

READ (slip). The prose on p7 says "the Earth is located at (mu, 0, 0) and the Moon at (1 - mu, 0, 0)". The formula for R_PM
two lines below, section 1.1 on p4 and Fig. 1 all put the Moon at X = mu - 1. The prose is a typesetting slip; the formulas are
used here.

READ. Eq. (4), p8: "The coefficients alpha_i, i = 1, ..., 8 are 2 pi-periodic real functions of the form
`alpha_i(theta) = a_0^i + sum_{k>=0} a_k^i cos(k theta) + sum_{k>=0} b_k^i sin(k theta)`" (the sums are printed with k >= 0;
the constant a_0 is separate, so the cosine sum effectively starts at k = 1; this is as printed). The variable theta is the
same angle as in the BCP, theta = omega_s t (INFERRED: the paper does not restate it for the QBCP, but the Hamiltonian is
written as a function of theta and the period T_s = 2 pi/omega_s is stated for the QBCP in eq. 5 and the text of 1.3.1).
READ. How many terms: not stated in this paper. The 2018 paper (Jorba-Cusco, Farres and Jorba) table runs to k = 13 (see
`core/qbcp.py`).

READ. Parity, p8: "A property of the coefficients alpha_i, i = 1, ..., 8, is that they are odd functions for i = 1, 3, 4, 7
and even for the rest. These properties imply that the following symmetry holds:
`H_QBCP(theta, X, Y, Z, P_X, P_Y, P_Z) = H_QBCP(-theta, X, -Y, Z, -P_X, P_Y, -P_Z)`".
INFERRED (our check, term by term): the printed symmetry holds exactly when alpha_1, alpha_3, alpha_4, alpha_6 and alpha_7 are
EVEN in theta (cosine series) and alpha_2, alpha_5, alpha_8 are ODD (sine series): under the map, alpha_2(P_X X + P_Y Y + P_Z Z)
changes sign (needs alpha_2 odd), alpha_3(P_X Y - P_Y X) is unchanged (needs alpha_3 even), alpha_4 X unchanged (even), alpha_5 Y
changes sign (odd), and (Y - alpha_8)^2 needs alpha_8 odd. That is the parity that Andreu (1998) prints and that the project
now uses. The sentence on p8 states the opposite words ("odd for 1, 3, 4, 7 and even for the rest"); the symmetry formula printed
directly after it, and Table 4 (x and p_y given with y = p_x = 0 at t = 0, a point on the reversing symmetry), agree with the
corrected parity, so the sentence appears to be a wording slip (an even/odd swap and the omission of i = 6). It is flagged here
only so that nobody takes the sentence as a second source for the sine/cosine assignment.

READ. Physical interpretation, p8: "alpha_1(theta), alpha_2(theta), alpha_3(theta) and alpha_6(theta) capture instantaneous
distance between the Earth and the Moon; alpha_4(theta) and alpha_5(theta) are the instantaneous Coriolis effect due to the
rotating reference frame; alpha_7(theta) and alpha_8(theta) capture the instantaneous position of the Sun within its plane of
motion."

READ. Momenta for the QBCP: the paper defines P_X = X' - Y etc. only for the RTBP and BCP. For the QBCP the P are the momenta
conjugate to (X, Y, Z) of eq. (3); no velocity relation is printed. INFERRED from eq. (3): `X' = dH/dP_X = alpha_1 P_X +
alpha_2 X + alpha_3 Y` and `Y' = alpha_1 P_Y + alpha_2 Y - alpha_3 X`, `Z' = alpha_1 P_Z + alpha_2 Z`. See section 6, item C1,
for a consequence for Table 4.

READ. Sun's sense in the QBCP: not stated and not testable from this paper alone (alpha_7, alpha_8 are not printed). It follows
from the sign of the first sine coefficient of alpha_8 in Andreu's table. INFERRED: the BCP of this paper has the Sun moving
clockwise (section 1.2), and the QBCP is the coherent version of the same system, so the same sense is expected.

### 1.4 Constants

Table 3 (p8), "specific values of the parameters used in this work", as printed:

| Symbol | Value |
| --- | --- |
| mu | 0.012150581600000 |
| m_s | 328900.5423094043 |
| omega_s | 0.925195985520347 |
| a_s | 388.8111430233511 |

READ: these are the QBCP constants. Compare with Table 2 (the BCP constants, section 1.2 above): mu differs at the 11th decimal
(0.0121505816 against 0.012150581623433), m_s differs by 7.7e-3 (328900.5423094043 against 328900.5499999991, relative 2.3e-8),
omega_s differs by 2.06e-12 (0.925195985520347 against 0.925195985518289), a_s is identical to all printed digits.
INFERRED (arithmetic): `T_s = 2 pi/omega_s = 6.791193871907917` time units with the Table 3 value (6.791193871923023 with the
Table 2 value). One time unit is about 4.348 days, so T_s is about 29.5 days, which is the "about 30 days" of Remark 1;
`6 T_s` is about 177 days, which the paper rounds to "approximately 180 days" (section 4, p31).
INFERRED, and important for the project: the mu of Table 3 is a rounding. The gamma values printed on p10 (section 3 below) are
reproduced to 1e-16 by the Earth-Moon RTBP only with mu = 0.012150581623433623, and are off by 1.1e-10 with mu = 0.0121505816
(our computation, see section 6, item C5). So the centre-manifold computation used the full mu, and Table 3's mu is a truncation
of it. The project's `core/qbcp.py` default constants (`_QBCP_MU_EM`, `_QBCP_MU_S`, `_QBCP_A_S`, `_QBCP_OMEGA_S`) are the BCP
constants of the 2018 paper (mu = 0.012150581623433623, m_S = 328900.54999999906, omega_S = 0.92519598551829646, a_S =
388.81114302335106), not the m_s and omega_s of Table 3; the older comment block in `core/qbcp.py` and the constants in
`core/bcr4bp.py` lines 11 to 14 quote Table 3 values. Whether the Fourier table was generated with the Table 3 or the Table 2
constants is not stated in this paper. The effect on the published orbits is below the project's current 2e-8 agreement
(INFERRED, not measured), but it is a loose end worth closing by trying both.

## 2. Tables, transcribed digit by digit

### 2.1 Table 4 (p9): initial conditions at t = 0 of the dynamical substitutes

READ. Caption: "Initial conditions at t = 0 of the dynamical substitutes for L_i, i = 1, 2. (L1, top and L2, bottom)". Footnote:
"Note that y = z = p_x = p_z = 0". Frame: the paper's frame (Earth at X = mu, Moon at X = mu - 1), canonical momentum p_y of
the QBCP Hamiltonian (eq. 3), pulsating coordinates, t = 0 meaning theta = 0.

| j | x | p_y |
| --- | --- | --- |
| 1 (L1, POL1) | -0.8369141677649317 | -0.8391311559808445 |
| 2 (L2, POL2) | -1.1556836078332600 | -1.1587306159501061 |

Checked against the page image and the text layer: identical. Both match the values in the second-hand digest
(`2026-06-14-andreu-quasi-bicircular-digest.md`) digit for digit.

### 2.2 Table 5 (p10): monodromy matrix eigenvalues of the dynamical substitutes

READ. Caption: "Monodromy matrix eigenvalues lambda_{i,j}, i = 1, 2, j = 1, 2, 3 of the dynamical substitutes for L_i, i = 1, 2.
(L1, top and L2, bottom)". Only three eigenvalues per orbit are listed; the other three are the inverses (Hamiltonian
system, as B states for the BCP). Period of the monodromy matrix: T_s. Printed as absolute value and argument (radians).

| i | j | abs(lambda_{i,j}) | arg(lambda_{i,j}) |
| --- | --- | --- | --- |
| 1 (L1) | 1 | 460182151.5759 | 0.000000000000 |
| 1 (L1) | 2 | 1.000000000000 | 2.871101174766 |
| 1 (L1) | 3 | 1.000000000000 | 2.981120162511 |
| 2 (L2) | 1 | 2397196.843443 | 0.000000000000 |
| 2 (L2) | 2 | 1.000000000000 | 0.408977840813 |
| 2 (L2) | 3 | 1.000000000000 | 0.091483781904 |

Stability (p9): "in the QBCP there are no changes of stability, and throughout the continuation process the stability type of
the periodic orbits is saddle x center x center for all values of epsilon in [0, 1]."

### 2.3 gamma values (p10)

READ. The scaling of step N = 1 makes "the unit of distance equal to the distance between the libration point studied and the
closest primary. We call this distance gamma_i":

| i | gamma_i |
| --- | --- |
| 1 | 0.1509342729900642 |
| 2 | 0.1678327317370704 |

INFERRED: these are the RTBP values (distance from the Moon to the RTBP L1 and L2), not the distances of the substitute orbits
(see section 6, item C5).

### 2.4 Normal frequencies (p10, step N = 2)

READ. "The normal frequencies chosen in each case are", with "kappa_1 corresponds to the hyperbolic part, and omega_1 and
omega_2 to the elliptical parts":

| Case | kappa_1 | omega_1 | omega_2 |
| --- | --- | --- | --- |
| POL1 | 2.93720564115629 | 2.27316022488810 | 2.33661946019073 |
| POL2 | 2.16306748237037 | 1.79017018257069 | 1.86386291350378 |

"Note that, for each case, these normal frequencies are very similar to their associated equilibrium points counterparts in the
RTBP." The paper adds the vector `omega = (kappa_1, i omega_1, i omega_2)`. READ (p17): the frequencies used by Andreu (2002)
for POL2 are `omega~_1 = 1.34709425E-02`, `omega~_2 = 2.16306748E+00` (hyperbolic), `omega~_3 = -6.02217885E-02`, "related to
the ones used here by `omega~_1 = omega_1 - 2 omega_s`, `omega~_3 = omega_2 - 2 omega_s`".
INFERRED (our arithmetic, Table 3 omega_s): `omega_1 - 2 omega_s = -0.0602218` and `omega_2 - 2 omega_s = +0.0134709`, so the
printed numbers agree with the pairing `omega~_1 <-> omega_2` and `omega~_3 <-> omega_1`; the subscripts in the printed relation
are interchanged relative to the numbers. A slip of indices only, the numerical values are consistent.

### 2.5 Radius of convergence and accuracy tables (pp11 to 17)

READ. Smallest small divisor delta_D (p11): "For the center manifold reduction computation around L1, the smallest value for
delta_D was delta_D ~ 0.011, and for the L2 case, delta_D ~ 0.013." Coefficients computed "up to degree N = 16" (POL1 and
POL2); N = 12 used for the L2 tests.

Table 6 (p12), radius of convergence r_n = 1/(||H_n||_1)^(1/n) for the reduced Hamiltonian:

| n | POL1 r_n | POL2 r_n |
| --- | --- | --- |
| 6 | 9.813101e-01 | 8.199574e-01 |
| 8 | 9.913491e-01 | 8.108276e-01 |
| 10 | 9.909848e-01 | 7.983601e-01 |
| 12 | 9.838444e-01 | 7.106946e-01 |
| 14 | 9.708615e-01 | 5.779491e-01 |
| 16 | 9.609837e-01 | 5.137823e-01 |

Table 7 (p13), POL1, N = 16, error ||v0 - v01||_2 against lambda0 (initial point (lambda0, lambda0, lambda0, lambda0)/2 in
centre-manifold coordinates, integration from t = 0 to t_f = 1):

| lambda0 | error | lambda0 | error |
| --- | --- | --- | --- |
| 0.125 | 2.532617e-10 | 0.250 | 3.989719e-08 |
| 0.150 | 3.631822e-10 | 0.275 | 1.817547e-07 |
| 0.175 | 5.019000e-10 | 0.300 | 7.241818e-07 |
| 0.200 | 1.267081e-09 | 0.325 | 2.579780e-06 |
| 0.225 | 7.452637e-09 | 0.350 | 8.355658e-06 |

Table 8 (p13), POL1, estimated truncation order N_j:

| j | lambda_j | lambda_{j+1} | N_j |
| --- | --- | --- | --- |
| 0 | 0.125 | 0.150 | 1.97717 |
| 1 | 0.150 | 0.175 | 2.09857 |
| 2 | 0.175 | 0.200 | 6.93523 |
| 3 | 0.200 | 0.225 | 15.04336 |
| 4 | 0.225 | 0.250 | 15.92378 |
| 5 | 0.250 | 0.275 | 15.90966 |
| 6 | 0.275 | 0.300 | 15.88740 |
| 7 | 0.300 | 0.325 | 15.87174 |
| 8 | 0.325 | 0.350 | 15.85841 |

Table 9 (p17), POL2, N = 12:

| lambda0 | error | lambda0 | error |
| --- | --- | --- | --- |
| 0.100 | 2.226642e-12 | 0.225 | 3.051407e-09 |
| 0.125 | 3.706322e-12 | 0.250 | 1.095514e-08 |
| 0.150 | 2.248650e-11 | 0.275 | 3.497710e-08 |
| 0.175 | 1.457249e-10 | 0.300 | 1.014818e-07 |
| 0.200 | 7.336179e-10 | 0.325 | 2.719555e-07 |

Table 10 (p17), POL2, N = 12 (the last row carries the index "9" as printed; the sequence of j would give 8):

| j | lambda_j | lambda_{j+1} | N_j |
| --- | --- | --- | --- |
| 0 | 0.100 | 0.125 | 2.28349 |
| 1 | 0.125 | 0.150 | 9.88844 |
| 2 | 0.150 | 0.175 | 12.12324 |
| 3 | 0.175 | 0.200 | 12.10403 |
| 4 | 0.200 | 0.225 | 12.10166 |
| 5 | 0.225 | 0.250 | 12.13174 |
| 6 | 0.250 | 0.275 | 12.18007 |
| 7 | 0.275 | 0.300 | 12.24192 |
| 9 | 0.300 | 0.325 | 12.31541 |

(Slip, stated factually: the text on p17 refers to "Table 7" and "Table 8" in the L2 discussion where Tables 9 and 10 are meant.)

### 2.6 Table 11 (p32): rotation numbers of the L2-Halo representatives plotted in Fig. 24

| Orbit | Rotation number rho |
| --- | --- |
| Blue | -0.0480876152458433 |
| Red | 3.6403791158911880 |
| Green | 1.0224171606049586 |

### 2.7 Table 12 (p37): the three QBCP Halo invariant curves used for transfers

| Invariant curve | Rotation number | lambda_u |
| --- | --- | --- |
| ICQ1 | 3.239814740891185 | 1269.060394604636 |
| ICQ2 | 1.022417160604956 | 58362.76296971765 |
| ICQ3 | 0.517157160604977 | 206452.6867125494 |

READ: curves of the stroboscopic map (period T_s), rotation number rho = 2 pi omega_1/omega_s (definition in section 3.1 of the 2021 paper), lambda_u the unstable
eigenvalue of the invariant curve. INFERRED: ICQ1's rotation number coincides to all 16 printed digits with the first row of
Table 3 of Rosales, Jorba and Jorba-Cusco (2021) (the BCP paper), a Halo orbit of the RTBP of energy -1.510315749412583; ICQ2
is the "green" orbit of Table 11 above (1.022417160604956 against 1.0224171606049586). So the QBCP tori were seeded from the same RTBP Halo orbits as the BCP
tori, and the rotation number was held fixed while the Sun's effect was switched on.

### 2.8 Table 13 (p39): transfer costs to the QBCP Halo orbits

READ. Caption: "Transfer cost to QBCP Halo orbits". Integration time 6 T_s; target a 200 km altitude Earth parking orbit
(R_E = 6400 km, p31). Columns: invariant curve, side of the unstable manifold (+ is the side between the Halo orbit and the
Moon, p34), cost function J1 (minimum delta-v), J2 (minimum time), J3 (minimum norm of delta-v and time), total delta-v in km/s,
time in days, latitude of the intersection with the LEO sphere in degrees.

| Curve | Side | Cost | delta-v (km/s) | t (days) | Latitude (deg) |
| --- | --- | --- | --- | --- | --- |
| ICQ1 | + | J1 | 3.2386 | 134.2429 | 10.710279 |
| ICQ1 | - | J1 | 3.2003 | 137.4482 | 6.415619 |
| ICQ1 | + | J2 | 3.8470 | 131.3539 | -18.440223 |
| ICQ1 | - | J2 | 3.3394 | 118.9735 | -2.317154 |
| ICQ1 | + | J3 | 3.8470 | 131.3539 | -18.440223 |
| ICQ1 | - | J3 | 3.3394 | 118.9735 | -2.317154 |
| ICQ2 | + | J1 | 3.2271 | 159.5806 | 18.505784 |
| ICQ2 | - | J1 | 3.1517 | 125.3764 | -13.777695 |
| ICQ2 | + | J2 | 6.3825 | 121.0911 | -54.610093 |
| ICQ2 | - | J2 | 3.2460 | 107.9764 | -4.959981 |
| ICQ2 | + | J3 | 3.7862 | 121.6507 | -21.937209 |
| ICQ2 | - | J3 | 3.2460 | 107.9764 | -4.959981 |
| ICQ3 | + | J1 | 3.1581 | 127.7909 | -5.262186 |
| ICQ3 | - | J1 | 3.1587 | 132.4915 | 5.678865 |
| ICQ3 | + | J2 | 3.7272 | 115.9231 | -19.960734 |
| ICQ3 | - | J2 | 3.2713 | 104.0634 | -6.622813 |
| ICQ3 | + | J3 | 3.7272 | 115.9231 | -19.960734 |
| ICQ3 | - | J3 | 3.1586 | 132.4914 | 5.678865 |

READ (text, p36): cheapest total delta-v is {ICQ2, -, J1} at 3.1517 km/s, 125.4 days; shortest is {ICQ3, -, J2}, about 104 days
at about 3.3 km/s; overall range 3.1517 km/s to slightly more than 13 km/s (Fig. 34). Distances to the invariant curve used:
2.5e-7 (about 100 m) for ICQ1, 7.5e-7 (about 290 m) for ICQ2, 7e-7 (about 270 m) for ICQ3 (p34); N = M = 1000 initial
conditions per fundamental cylinder.

### 2.9 Tables 14 and 15 (Appendix A, pp42 to 45): Hamiltonian reduced to the centre manifold up to order 6

READ. Caption: "Hamiltonian reduced to the central manifold up to order 6 around POL1" (Table 14) and "... around POL2"
(Table 15). Form (eq. 7, p12): `H_k = sum_{k1+k2+k3+k4 = k} a(k1,k2,k3,k4) Q1^k1 P1^k2 Q2^k3 P2^k4`. The variables are the
scaled, Floquet- and Lie-normalised coordinates of section 2 (unit distance = gamma_i, normal frequencies of 2.4), not
synodic coordinates. The paper prints two columns of (k1 k2 k3 k4, a) per row; the lists below are sorted by degree and then
by exponents, which changes nothing but the order. Printed digits kept (13 decimals in mantissa).

Table 14, POL1 (77 coefficients):

| k1 | k2 | k3 | k4 | a(k1,k2,k3,k4) |
| --- | --- | --- | --- | --- |
| 2 | 0 | 0 | 0 | 1.1365801124440E+00 |
| 0 | 2 | 0 | 0 | 1.1365801124440E+00 |
| 0 | 0 | 2 | 0 | 1.1683097300953E+00 |
| 0 | 0 | 0 | 2 | 1.1683097300953E+00 |
| 2 | 0 | 1 | 0 | -4.2742797554386E-01 |
| 1 | 1 | 0 | 1 | -1.2254290645138E-04 |
| 0 | 2 | 1 | 0 | -5.3891327233143E-05 |
| 0 | 0 | 3 | 0 | 2.5523418206125E-02 |
| 0 | 0 | 1 | 2 | -4.9529829287648E-01 |
| 4 | 0 | 0 | 0 | -1.0387633163417E-01 |
| 2 | 2 | 0 | 0 | 8.5654706992094E-02 |
| 2 | 0 | 2 | 0 | 2.1622139838010E-01 |
| 2 | 0 | 1 | 1 | -9.1489731924294E-08 |
| 2 | 0 | 0 | 2 | -2.4182953302687E-01 |
| 1 | 1 | 2 | 0 | -1.4863019899213E-09 |
| 1 | 1 | 1 | 1 | -3.2495127186968E-02 |
| 1 | 1 | 0 | 2 | 1.4246394326115E-09 |
| 0 | 4 | 0 | 0 | 1.0812958900733E-05 |
| 0 | 2 | 2 | 0 | -1.5360957052390E-02 |
| 0 | 2 | 1 | 1 | -3.2599524388118E-08 |
| 0 | 2 | 0 | 2 | 9.9396670609705E-02 |
| 0 | 0 | 4 | 0 | -1.5779796388201E-02 |
| 0 | 0 | 3 | 1 | -2.1506277895067E-08 |
| 0 | 0 | 2 | 2 | 2.8794821677007E-01 |
| 0 | 0 | 1 | 3 | -9.9796480381400E-08 |
| 0 | 0 | 0 | 4 | -1.4074479895471E-01 |
| 4 | 0 | 1 | 0 | 3.7745746907786E-02 |
| 4 | 0 | 0 | 1 | -6.0213290782196E-08 |
| 3 | 1 | 1 | 0 | -3.1696042934014E-09 |
| 3 | 1 | 0 | 1 | -6.3675915101523E-02 |
| 2 | 2 | 1 | 0 | -1.2726077950140E-01 |
| 2 | 2 | 0 | 1 | 8.9009096319913E-09 |
| 2 | 0 | 3 | 0 | -1.1083737547103E-01 |
| 2 | 0 | 0 | 3 | -1.2846150897253E-07 |
| 1 | 3 | 0 | 1 | 1.7507059394155E-02 |
| 1 | 1 | 3 | 0 | -6.0652645741202E-09 |
| 1 | 1 | 0 | 3 | -7.4012409047879E-02 |
| 0 | 4 | 1 | 0 | 1.0507633803701E-02 |
| 0 | 2 | 3 | 0 | 2.2665985616829E-02 |
| 0 | 2 | 0 | 3 | 9.6176086630088E-09 |
| 0 | 0 | 5 | 0 | 1.1494979183962E-02 |
| 0 | 0 | 2 | 3 | 8.0037887245600E-08 |
| 0 | 0 | 1 | 4 | 1.4436762488537E-01 |
| 0 | 0 | 0 | 5 | -6.7934034132082E-08 |
| 6 | 0 | 0 | 0 | 6.2094210958681E-03 |
| 5 | 1 | 0 | 0 | -6.7815393271166E-09 |
| 4 | 2 | 0 | 0 | -2.0086057404615E-02 |
| 4 | 0 | 2 | 0 | -2.1866965033703E-03 |
| 4 | 0 | 1 | 1 | 7.1268148429479E-08 |
| 4 | 0 | 0 | 2 | 2.0982296568260E-02 |
| 3 | 1 | 2 | 0 | 1.6754220796143E-09 |
| 3 | 1 | 1 | 1 | 8.6670022069764E-02 |
| 3 | 1 | 0 | 2 | -3.3030448427256E-08 |
| 2 | 4 | 0 | 0 | 2.7250801950251E-02 |
| 2 | 2 | 2 | 0 | 1.0778375887518E-01 |
| 2 | 2 | 1 | 1 | -4.8046664060818E-08 |
| 2 | 2 | 0 | 2 | -3.1233429812347E-02 |
| 2 | 0 | 4 | 0 | 5.0908816751363E-02 |
| 2 | 0 | 3 | 1 | -7.0795221695100E-08 |
| 2 | 0 | 2 | 2 | -8.3158755543959E-02 |
| 2 | 0 | 1 | 3 | 1.7764130617398E-07 |
| 1 | 3 | 2 | 0 | 1.0666533414263E-08 |
| 1 | 3 | 1 | 1 | -3.7167629573763E-02 |
| 1 | 3 | 0 | 2 | -4.9317299850824E-09 |
| 1 | 1 | 4 | 0 | 1.0801455668335E-08 |
| 1 | 1 | 3 | 1 | -1.0655124578491E-01 |
| 1 | 1 | 2 | 2 | -1.4711027029686E-09 |
| 0 | 6 | 0 | 0 | -1.2378414626373E-03 |
| 0 | 4 | 2 | 0 | -8.5673296189717E-03 |
| 0 | 4 | 1 | 1 | -1.6472545989214E-08 |
| 0 | 4 | 0 | 2 | 2.0641390341789E-02 |
| 0 | 2 | 4 | 0 | -1.2873305122172E-02 |
| 0 | 2 | 3 | 1 | -3.9935005395564E-08 |
| 0 | 2 | 2 | 2 | 1.0852170163338E-01 |
| 0 | 0 | 6 | 0 | -5.5676966532490E-03 |
| 0 | 0 | 5 | 1 | -2.4022343097574E-08 |
| 0 | 0 | 4 | 2 | 1.2851003627049E-01 |

Table 15, POL2 (68 coefficients):

| k1 | k2 | k3 | k4 | a(k1,k2,k3,k4) |
| --- | --- | --- | --- | --- |
| 2 | 0 | 0 | 0 | 8.9508509128534E-01 |
| 0 | 2 | 0 | 0 | 8.9508509128534E-01 |
| 0 | 0 | 2 | 0 | 9.3193145675189E-01 |
| 0 | 0 | 0 | 2 | 9.3193145675189E-01 |
| 2 | 0 | 1 | 0 | 6.5589636328480E-05 |
| 1 | 1 | 0 | 1 | -1.4657320225294E-04 |
| 0 | 2 | 1 | 0 | 6.4841149489243E-01 |
| 0 | 0 | 3 | 0 | -6.4947365185738E-02 |
| 0 | 0 | 1 | 2 | 8.3042596977058E-01 |
| 4 | 0 | 0 | 0 | 1.6691540956563E-05 |
| 2 | 2 | 0 | 0 | 1.6501717240559E-01 |
| 2 | 0 | 2 | 0 | -4.9579201703060E-02 |
| 2 | 0 | 0 | 2 | 2.1143854714294E-01 |
| 1 | 1 | 1 | 1 | 1.0973656675138E-01 |
| 0 | 4 | 0 | 0 | -1.8016477271676E-02 |
| 0 | 2 | 2 | 0 | 3.5651315214778E-01 |
| 0 | 2 | 0 | 2 | -4.7292944632242E-02 |
| 0 | 0 | 4 | 0 | -4.1231015606744E-02 |
| 0 | 0 | 2 | 2 | 5.9236862155832E-01 |
| 0 | 0 | 0 | 4 | -3.1058453169198E-02 |
| 4 | 0 | 1 | 0 | -4.3777802018475E-02 |
| 3 | 1 | 0 | 1 | 7.4667616272107E-02 |
| 2 | 2 | 1 | 0 | 2.8508460013478E-01 |
| 2 | 2 | 0 | 1 | -1.0319019428507E-09 |
| 2 | 0 | 3 | 0 | -7.6670187245196E-02 |
| 2 | 0 | 1 | 2 | 2.8479184552457E-01 |
| 1 | 3 | 0 | 1 | -1.3880815534462E-01 |
| 1 | 1 | 2 | 1 | 4.1875686746481E-01 |
| 1 | 1 | 1 | 2 | -1.6240816892809E-09 |
| 1 | 1 | 0 | 3 | -1.7844450052689E-01 |
| 0 | 4 | 1 | 0 | -8.3453433400644E-03 |
| 0 | 2 | 3 | 0 | 1.9426009938526E-01 |
| 0 | 2 | 2 | 1 | -1.5463320553586E-09 |
| 0 | 2 | 1 | 2 | -2.4495800601342E-01 |
| 0 | 0 | 5 | 0 | -3.1013224023379E-02 |
| 0 | 0 | 3 | 2 | 5.8051203522045E-01 |
| 0 | 0 | 2 | 3 | -3.0394344483381E-09 |
| 0 | 0 | 1 | 4 | -3.0140880721764E-01 |
| 6 | 0 | 0 | 0 | -6.4307281988146E-03 |
| 4 | 2 | 0 | 0 | 8.1725097260177E-02 |
| 4 | 0 | 2 | 0 | -2.7581282162579E-02 |
| 4 | 0 | 0 | 2 | 4.1938736410204E-02 |
| 3 | 1 | 1 | 1 | 1.9906619251976E-01 |
| 2 | 4 | 0 | 0 | -4.2728780806097E-03 |
| 2 | 2 | 2 | 0 | 3.0570142682015E-01 |
| 2 | 2 | 1 | 1 | -1.7463324777085E-09 |
| 2 | 2 | 0 | 2 | -4.3943735317670E-02 |
| 2 | 0 | 4 | 0 | -3.0925874741699E-02 |
| 2 | 0 | 2 | 2 | 3.5211138710868E-01 |
| 2 | 0 | 1 | 3 | -1.9587319456264E-09 |
| 2 | 0 | 0 | 4 | -6.3932143025888E-03 |
| 1 | 3 | 1 | 1 | -1.5531042568898E-01 |
| 1 | 1 | 3 | 1 | 4.7597114070644E-01 |
| 1 | 1 | 2 | 2 | -3.0449382625031E-09 |
| 1 | 1 | 1 | 3 | -2.9968970918954E-01 |
| 0 | 6 | 0 | 0 | -1.3308183673882E-02 |
| 0 | 4 | 2 | 0 | 4.1375312168077E-02 |
| 0 | 4 | 0 | 2 | -7.3614758826606E-02 |
| 0 | 2 | 4 | 0 | 9.6429491577036E-02 |
| 0 | 2 | 3 | 1 | -3.0413181902431E-09 |
| 0 | 2 | 2 | 2 | -1.0297348505084E-01 |
| 0 | 2 | 0 | 4 | -1.2340310730044E-01 |
| 0 | 0 | 6 | 0 | -1.0289658815507E-02 |
| 0 | 0 | 4 | 2 | 4.5913199291929E-01 |
| 0 | 0 | 3 | 3 | -5.6536567341678E-09 |
| 0 | 0 | 2 | 4 | -2.7360468675825E-01 |
| 0 | 0 | 1 | 5 | -1.2601762375995E-09 |
| 0 | 0 | 0 | 6 | -6.5234840557094E-02 |

INFERRED (our arithmetic): the degree-2 coefficients are half the normal frequencies of 2.4: for POL1,
a(2,0,0,0) = a(0,2,0,0) = 1.1365801124440 = omega_1/2 (1.13658011244405) and a(0,0,2,0) = a(0,0,0,2) = 1.1683097300953 = omega_2/2
(1.1683097300953); for POL2, 0.89508509128534 = omega_1/2 (0.895085091285345) and 0.93193145675189 = omega_2/2 (0.93193145675189).
These are a consistency check on the printed tables, not an independent test of any module.

## 3. The dynamical substitutes of the collinear points

READ (section 1.3.1, p8): "In the QBCP, the collinear points in the RTBP are replaced by small periodic orbits with the same
period as the perturbation, T_s = 2 pi/omega_s. These orbits are computed by continuation from the RTBP to the QBCP... we
consider the family of Hamiltonians H^epsilon = H_RTBP + epsilon (H_QBCP - H_RTBP), epsilon in [0, 1] (eq. 5)... we start the
continuation scheme from a collinear equilibrium point (L_i, i = 1, 2) and epsilon = 0, then the value of epsilon is increased
until it reaches epsilon = 1... For each value of epsilon in [0, 1], there is a T_s-periodic orbit." Note the homotopy here
interpolates between the whole QBCP and the RTBP, not (as for the BCP) the Sun's mass.

READ: "In all two cases, there is a direct connection between the starting point and the final periodic orbit. We recall that,
in the BCP, where L2 is connected with a 1:2 resonant planar Lyapunov orbit (see Jorba-Cusco et al. 2018). We remark that this
does not happen for the QBCP." Stability: saddle x center x center for all epsilon (p9). Names: POL1 and POL2 (p9).

READ (Fig. 2, p9): first column, x(t = 0) against epsilon: for L1 the curve runs from x about -0.83691 (epsilon = 0, the RTBP
point) out to x about -0.83693 near epsilon = 0.5 and back to the printed Table 4 value at epsilon = 1 (graph read; the
plotted x-axis is -0.83693 to -0.83691). Second column, the orbits at epsilon = 1 in the (x, y) plane: for L1 they span about
-0.836916 to -0.836914 in x and about +-3e-6 in y; for L2 about -1.155683 to -1.155681 in x and about +-4e-6 in y (graph read).
INFERRED from Table 4 and the gamma values: the t = 0 position of POL1 is 9.8e-7 and of POL2 1.46e-6 from the RTBP point
(section 6, item C5), consistent with "of order 1e-6". The 2021 paper states of the QBCP: "L2 is replaced by a periodic orbit
that is small in the sense that its maximal distance to L2 is of the order of 1e-6, and it has the same stability type of the L2
point".

Eigenvalues: Table 5 (2.2). Hyperbolic multipliers 4.6e8 (L1) and 2.4e6 (L2); the elliptic pair on the unit circle.
INFERRED (our arithmetic, section 6, item C3): the hyperbolic exponent ln(abs(lambda_1))/T_s reproduces the printed kappa_1 to
2e-14 (2.93720564115631 against 2.93720564115629 for POL1, 2.16306748237038 against 2.16306748237037 for POL2), and the
printed arguments equal omega_j T_s modulo 2 pi to 12 digits up to sign: POL1 omega_1 T_s = 15.4374717891 = 2.871101174766 +
4 pi, omega_2 T_s = 15.8684357590 = -2.981120162511 + 6 pi; POL2 omega_1 T_s = 12.1573927735 = -0.408977840813 + 4 pi,
omega_2 T_s = 12.6578543963 = +0.091483781904 + 4 pi. So Tables 5 and the normal-frequency table are one set of data.

## 4. Families of tori and periodic orbits computed

READ (section 3, p20): the families are computed by numerical continuation of 2D invariant tori (invariant curves of the
stroboscopic map at time T_s), parametrised by the rotation number rho, "with the algorithms described in Jorba (2001) and
Rosales et al. (2021a)". The paper also reduces to the centre manifold up to degree 16 (POL1) or 12 (POL2) and shows
Poincare sections (Figs. 4, 5, 7, 8) at energy levels (e.g. h = 0.2, 0.4, 0.7, 0.9 for POL1, p15) of the reduced Hamiltonian.
Values are not tabulated, and the centre-manifold energy is not the synodic Hamiltonian value.

L1 (section 3.1, pp21 to 22):
- Vertical family (quasi-periodic vertical Lyapunov), born from POL1 along the vertical direction. Tori shown: rho =
  2.8710835247657562 (Fig. 10, "very small, and close to the periodic orbit that replaces L1"), rho = 1.7158771247657665
  (Fig. 11), rho = 1.0158771247657681 (Fig. 12, large). Stability: one real pair, the largest eigenvalue "of the order of 1e8
  and decreases with the rotation number until a value of the order of 1e6"; the other pair complex of norm 1 (partially
  elliptic). "No bifurcations were identified", although Fig. 4 implies one; "the step size ... probably jumped over the
  bifurcation." Fig. 9 has "a sharp turn" between x = 0.13 and x = 0.14 on the third (vertical) coordinate axis, which the authors
  say "reminds [of] a pitchfork bifurcation obtained by symmetry breaking" and that they were unable to verify.
- Horizontal family L1-HLy (quasi-periodic planar Lyapunov), born from POL1 along the planar frequency. A bifurcation on it
  (identified in the stability analysis: the last eigenvalue pair, complex of norm 1 at first, becomes real) gives the family
  L1-QV (Fig. 13, 14). Graph read: the bifurcation is at rho about 3.25 on Fig. 14, and in Fig. 13 the L1-QV branch leaves L1-HLy
  near x about -0.825; L1-HLy runs from rho about 3 (x about -0.835) to rho about 6 (x about -0.80).
- L1-Halo: Halo orbits of the RTBP continued to the QBCP and then continued in the QBCP (purple in Fig. 13). "Numerical evidences
  suggest that these two families [L1-Halo and L1-QV] are not connected." Compared pair: L1-Halo rho = 3.4622727594120977 and
  L1-QV rho = 3.4623791625106679 (Fig. 15): "Both orbits are different in size and position." The L1-QV representative is "a
  Halo-like orbit". L1-Halo is mostly elliptic; L1-QV undergoes a bifurcation from elliptic to hyperbolic (Fig. 16, 17).
  (The text calls the QV family "L1-Q1" once, p22; a slip of the label.)

L2 (section 3.2, pp22 to 24, 29):
- Vertical family from POL2: tori rho = -0.4089841068128386 (Fig. 19, close to the reference periodic orbit),
  -0.8717553068128412 (Fig. 20), -1.0173803068128409 (Fig. 21). Largest real eigenvalue from order 1e6 down to 1e5; the other
  pair complex of norm 1 until the end of the family (rho about -1.0179) where "it seems that the two eigenvalues become real".
- Horizontal family L2-HLy: bifurcation (Fig. 23, graph read at rho about -0.06) gives L2-QV (rho = -0.0721362180958642 and
  -0.2449362180958645, Figs. 26, 27, "not Halo-like"). "In Andreu (1998) three other small bifurcations were found. These were
  not noticed here, probably because the step size ... was not small enough."
- L2-Halo: RTBP Halo orbits continued to the QBCP. Rotation numbers of three representatives (Table 11, section 2.6).
  "the family L2-QV and L2-Halo are not connected", but "the L2-Halo family connects to another family of 2D tori resonant with
  the frequency of the Sun" near (x, rho) about (-1.12, -0.05) (Fig. 22, right, graph read). "This connection was conjectured in
  Andreu (1998), and the numerical evidence provided here seems to prove it." A member of the resonant family: rho =
  -0.0774976152458405 (Fig. 25).
- Stability: L2-Halo mostly elliptic with "small pockets of real eigenvalues"; the other two eigenvalues real, "between 1e2 and
  1e6" for L2-Halo and "between 1e5 and 1e6" for L2-QV and the resonant family; L2-QV tori all real eigenvalues (Fig. 29).

READ (conclusion, p38): "the family of out-of-plane orbits born from the bifurcation seemed not to be the RTBP Halo
counterparts in the QBCP. The RTBP Halo orbits do survive in the QBCP, but do not seem to be connected to the quasi-periodic
planar Lyapunov family."

## 5. Statements comparing the BCP and the QBCP, and objects absent from one

All READ, quoted:
- p3 (introduction): "Notice that the center manifold of L2 in the QBCP has also been analyzed in Le Bihan et al. (2017b) by means
  of the parameterization method. And in Andreu (1998, 2002), using the same approach as in this work". Also p3: "we use the
  center manifold approach to study L2 as the QBCP has a similar qualitative behavior to the RTBP (while the BCP has not).
  Notice that this behavior is also observed in the high-fidelity model used in Lian (2013)."
- p6, known facts on the BCP: "When it comes to L2, in Jorba-Cusco et al. (2018) it is shown that there is no dynamical equivalent
  of L2 in the BCP. Indeed, the dynamical equivalent of L2 merges with a 1:2 resonant horizontal Lyapunov orbit. However, at some
  distance of L2 the model displays common features with the RTBP."
- p6: the BCP has Type I and Type II Halo families of tori near L2 (Rosales et al. 2021a); "Type II is a family of quasi-Halo orbits
  which is in 1:2 resonance with the Sun". Whether the QBCP has a Type II counterpart is not stated here.
- p9: in the QBCP the substitute of L2 connects directly to the starting point: "We remark that this does not happen for the QBCP."
- 2021 paper, p9 of that paper (see its digest): "The comparison between the BCP and QBCP illustrates this phenomena."
- Abstract and section 4: manifolds of L2 tori that pass close to the Earth "[are] not observed in the RTBP"; they were found in the
  BCP too (Rosales et al. 2021b), and section 4 repeats the analysis: "the behavior of the cases studied in the QBCP are pretty
  similar to their counterparts in the BCP".
- p6: the BCP (Jorba et al. 2020) "horizontal family of Lyapunov of invariant tori undergoes a 1:1 resonance and bifurcate producing
  a Halo family of invariant tori"; this paper finds for the QBCP that the planar Lyapunov family bifurcates to an out-of-plane
  family (L1-QV, L2-QV) that is not the continuation of the RTBP Halo family.

## 6. Positive controls for the project

Frame recipe (INFERRED from the printed frame and the project's convention): the paper's frame has Earth at X = mu and the Moon at
X = mu - 1, with the Sun at angle 0 at t = 0 (theta = 0). The project's frame is rotated by pi about the z axis: Earth at -mu,
Moon at 1 - mu; positions x, y, velocities x', y' and momenta p_x, p_y change sign; z and p_z do not; the Sun starts at angle pi.
The rotation preserves the sense of the Sun's motion, so the Sun is still clockwise. Monodromy eigenvalues and normal frequencies are
unchanged by the rotation. Time zero is the same (theta = 0). For Table 4:

| Orbit | project-frame state at t = 0 (x, y, z, p_x, p_y, p_z) |
| --- | --- |
| POL1 | (+0.8369141677649317, 0, 0, 0, +0.8391311559808445, 0) |
| POL2 | (+1.1556836078332600, 0, 0, 0, +1.1587306159501061, 0) |

Already tested by the project (`tests/core/test_qbcp.py`, end of file): Table 4 states (`_POL1_PUBLISHED`, `_POL2_PUBLISHED`;
position stays within 1e-5 of the point for a quarter period each way; multiple-shooting orbit matches the printed x and p_y to
1e-7; code docstring: 2e-8) and the largest multiplier to 1e-7 against the 2018 paper's table (4.60182151e8, 2.39719684e6).
New candidates, in order of value:

C1. Velocity at t = 0 from the printed Hamiltonian (new; also a correction to a project note). INFERRED from eq. (3): at theta = 0
with y = p_x = 0, alpha_2(0) = 0 (sine series), so `x' = 0` and `y' = alpha_1(0) p_y - alpha_3(0) x`. With the project's tables
(`core/qbcp.py`, alpha_1(0) = 1.0169220392610525, alpha_3(0) = 1.0196084040004443, both computed by us as sums of the cosine
coefficients), the paper-frame velocity is y' = -6.2475e-6 for POL1 and +6.0180e-6 for POL2 (project frame: opposite signs),
which is the size expected for an orbit of extent 1e-6 over a time of order one. The relation `v_y = p_y - x` (as used in the
project note on Andreu's canonical-momentum states, which is the RTBP relation) would give -2.2170e-3 for POL1 and -3.0470e-3 for
POL2, three orders of magnitude too large, so it must not be used with Table 4 in the QBCP. `core/qbcp.py` takes (x, p) as the
state, so the module is not affected; any helper that converts Table 4 to velocities is.

C2. Full Table 5 (new beyond the largest multiplier). Recipe: integrate the state variational equations of the project's QBCP over
one period T_s = 2 pi/omega_s (6.7911938719 time units) from the project-frame state above, take the eigenvalues of the 6 x 6
monodromy matrix. Compare as sets: abs 460182151.5759 and 2397196.843443; unit-modulus pairs with arguments +-2.871101174766 and
+-2.981120162511 (L1), +-0.408977840813 and +-0.091483781904 (L2). Agreement the digits allow: the printed moduli have 13
significant digits; the hyperbolic multiplier amplifies a 1e-16 rounding of the printed state by 5e-8 relative at L1 (4.6e8 times
1e-16, INFERRED) but the STM eigenvalue itself depends smoothly on the state, so 1e-9 relative is a fair target and 1e-6 a safe
first tolerance; the arguments should reach 1e-9 absolute (INFERRED; tighten after the first measurement). The project test uses
only the 2018 paper's 9-digit largest multiplier.

C3. kappa_1 and omega_j (new, derived). `kappa_1 = ln(abs(lambda_1))/T_s` reproduces 2.93720564115629 and 2.16306748237037; the
elliptic exponents `omega_j` equal the arguments over T_s plus integer multiples of omega_s (printed omega_1, omega_2 are
particular lifts, section 3 above; the lifts for Andreu's frequencies are omega - 2 omega_s). Not independent of C2 beyond the lift
convention; useful because the printed lifts are the ones the centre-manifold tables use.

C4. Two constants sets. Run C1 and C2 with (a) the project's current 2018 BCP constants and (b) Table 3's m_s = 328900.5423094043,
omega_s = 0.925195985520347 (and mu = 0.012150581623433623, see 1.4). If (b) reproduces the 12-digit arguments better than (a),
the Fourier table belongs with the Table 3 constants. This is a way to close the loose end in 1.4 and costs two runs.

C5. gamma values and substitute offsets (new, project-independent of the QBCP). The RTBP L1 and L2 distances to the Moon from the
collinear quintics with mu = 0.012150581623433623 are 0.15093427299006432 and 0.1678327317370705, against the printed
0.1509342729900642 and 0.1678327317370704 (agreement 1e-16); with mu = 0.0121505816 they are 0.15093427289819 and 0.16783273162351
(off by 1.1e-10). The substitute orbits are at distance 0.1509352506350683 (POL1) and 0.16783418943326 (POL2) from the Moon at
t = 0, so 9.78e-7 and 1.458e-6 from the RTBP points (our arithmetic, mu = 0.0121505816 and the paper-frame Moon at X = mu - 1; the
choice of mu changes these by 1e-11). Consistent with the project's "within 3e-6" check.

C6. Rotation numbers of tori (new, needs a torus solver; the project has none in the QBCP). The first torus of each vertical family
has a rotation number within 2e-5 of the corresponding printed multiplier argument: L1 2.8710835247657562 against 2.871101174766
(difference 1.8e-5), L2 -0.4089841068128386 against -0.408977840813 (6.3e-6) (INFERRED: in the limit of vanishing amplitude the
rotation number tends to the argument). Further printed rotation numbers: sections 2.6, 2.7, 4. The unstable eigenvalues of the
three Halo curves (Table 12) and the transfer costs (Table 13) require the invariant curves themselves, not given.

C7. The centre-manifold coefficients (Tables 14 and 15): not a realistic control; they depend on the normalisation (scaling by
gamma_i, Floquet gauge, Lie transformation, choice of omega lifts) and need an implementation of the whole normal-form algorithm.
The degree-2 coefficients are the exceptions (omega_j/2, section 2.9).

Not testable: the transfer figures; the Sun's sense in the QBCP (no alpha table printed here); the number of Fourier terms.

## 7. Corrections to the existing second-hand digest (`2026-06-14-andreu-quasi-bicircular-digest.md`)

1. The Hamiltonian block prints `- (1-mu)/R_PE - mu/R_PM ... - m_S/(alpha_6 R_PS)`, i.e. alpha_6 only on the Sun, as a divisor.
   Source (eq. 3, p7): alpha_6 multiplies the whole Newtonian potential, as a factor. (This is the defect corrected in
   `core/qbcp.py` as item 2 of its 2026-10-04 note.) The gloss "alpha_6 scales the Sun distance" is also wrong: the paper says
   alpha_1, alpha_2, alpha_3 and alpha_6 capture the instantaneous Earth-Moon distance.
2. The digest says the Rosales 2023 Hamiltonian has the Sun term only as `m_S/(alpha_6 R_PS)` and cites the Frontiers paper for
   alpha_7, alpha_8. Source: alpha_7 and alpha_8 appear in R_PS of eq. (3) of this paper itself.
3. "the alpha_i are odd/even under (theta, x, y, z) -> (-theta, x, -y, z)": the source prints `(theta, X, Y, Z, P_X, P_Y, P_Z) ->
   (-theta, X, -Y, Z, -P_X, P_Y, -P_Z)`, the momenta flipping as well, and states the parity in words that disagree with the
   formula (see 1.3).
4. "The Sun's perturbation is O(epsilon^2) in size (Coriolis and linear-order Sun terms cancel)": not found in this paper; not
   a statement of the source. The indirect term cancels the first-order term of the Sun's potential about the barycentre in the
   BCP (section 1.2), which is a different statement. Treat as unsupported.
5. The digest's table header "Gimeno 2018" is a misattribution of the paper by Jorba-Cusco, Farres and Jorba (2018) (this paper's
   reference list: "Jorba-Cusco, M., Farres, A., Jorba, A.: Two periodic models for the Earth-Moon system. Front. Appl. Math.
   Stat. 4, 32 (2018)"). The 2018 column of that table is the BCP constants; Table 3 here is the QBCP set (1.4).
6. Table 4 values, the Table 3 constants and the period (6.7912, about 30 days) agree with the source to all printed digits.
7. The project note on Andreu's canonical-momentum states ("convert via vy = py - x") should be revised: that relation is the RTBP
   one; the QBCP relation is `y' = alpha_1 p_y - alpha_3 x` at the symmetry point (C1 above).
8. The digest says halo ICs are not tabulated: confirmed for this paper; the rotation numbers it prints are in 2.6, 2.7, 4.

## 8. References cited by the paper that bear on this project (full citations as printed)

Quasi-bicircular model and coefficients:
- Andreu, M.A.: The quasi-bicircular problem. Ph.D. thesis, University of Barcelona (1998). (Source of the Fourier coefficients.)
- Andreu, M.A.: Dynamics in the center manifold around L2 in the quasi-bicircular problem. Celest. Mech. Dyn. Astron. 84(2),
  105-133 (2002).
- Gabern, F., Jorba, A.: A restricted four-body model for the dynamics near the Lagrangian points of the Sun-Jupiter system.
  Discrete Contin. Dyn. Syst. Ser. B 1(2), 143-182 (2001).
- Jorba-Cusco, M., Farres, A., Jorba, A.: Two periodic models for the Earth-Moon system. Front. Appl. Math. Stat. 4, 32 (2018).
- Le Bihan, B., Masdemont, J.J., Gomez, G., Lizy-Destrez, S.: Systematic study of the dynamics about and between the libration
  points of the Sun-Earth-Moon system. In: International Symposium on Space Flight Dynamics (ISSFD), pp. 1-10, Matsuyama, JP
  (2017a).
- Le Bihan, B., Masdemont, J.J., Gomez, G., Lizy-Destrez, S.: Invariant manifolds of a non-autonomous quasi-bicircular problem
  computed via the parameterization method. Nonlinearity 30(8), 3040 (2017b).

Bicircular model:
- Huang, S.S.: Very restricted four-body problem. Technical note TN D-501, Goddard Space Flight Center, NASA (1960).
- Cronin, J., Richards, P.B., Russell, L.H.: Some periodic solutions of a four-body problem. Icarus 3, 423-428 (1964).
- Gomez, G., Jorba, A., Masdemont, J., Simo, C.: Study of Poincare maps for orbits near Lagrangian points. ESOC contract
  9711/91/D/IM(SC), final report, European Space Agency, 1993. Reprinted as Dynamics and mission design near libration points.
  Vol. IV, Advanced methods for triangular points, volume 5 of World Scientific Monograph Series in Mathematics (2001). (The text
  cites "Gomez et al. (2001)" for the derivation of the BCP equations; the reprint of this report is the 2001 item with that
  author list, INFERRED.)
- Simo, C., Gomez, G., Jorba, A., Masdemont, J.: The bicircular model near the triangular libration points of the RTBP. In: Roy,
  A.E., Steves, B.A. (eds.) From Newton to Chaos, pp. 343-370. Plenum Press, New York (1995).
- Scheeres, D.J.: The restricted Hill four-body problem with applications to the Earth-Moon-Sun system. Celest. Mech. Dyn.
  Astron. 70(2), 75-98 (1998).
- Jorba, A., Jorba-Cusco, M., Rosales, J.J.: The vicinity of the Earth-Moon L1 point in the bicircular problem. Celest. Mech. Dyn.
  Astron. 132(2), 11 (2020).
- Rosales, J.J., Jorba, A., Jorba-Cusco, M.: Families of Halo-like invariant tori around L2 in the Earth-Moon bicircular problem.
  Celest. Mech. Dyn. Astron. 133, 03 (2021a). (The article number is printed 16 in the paper itself; "03" is as in this reference list.)
- Rosales, J.J., Jorba, A., Jorba-Cusco, M.: Transfers from the Earth to L2 Halo orbits in the Earth-Moon bicircular problem.
  Celest. Mech. Dyn. Astron. 133, 12 (2021b). (Not yet in the corpus as far as this digest knows; it holds the BCP version of
  section 4.)
- Jorba, A., Nicolas, B.: Transport and invariant manifolds near L3 in the Earth-Moon bicircular model. Commun. Nonlinear Sci.
  Numer. Simul. 89, 105327 (2020); and: Using invariant manifolds to capture an asteroid near the L3 point of the Earth-Moon
  Bicircular model. Commun. Nonlinear Sci. Numer. Simul. 102, 105948 (2021).

Earlier periodic-orbit and torus continuation in these models, and methods:
- Gomez, G., Mondelo, J.: The dynamics around the collinear equilibrium points of the RTBP. Physica D 157(4), 283-321 (2001).
- Jorba, A., Masdemont, J.: Dynamics in the center manifold of the collinear points of the restricted three body problem. Physica D
  132, 189-213 (1999).
- Jorba, A.: A methodology for the numerical computation of normal forms, centre manifolds and first integrals of Hamiltonian
  systems. Exp. Math. 8(2), 155-195 (1999).
- Jorba, A.: Numerical computation of the normal behaviour of invariant curves of n-dimensional maps. Nonlinearity 14(5), 943-976
  (2001).
- Jorba, A., Villanueva, J.: On the persistence of lower dimensional invariant tori under quasi-periodic perturbations. J. Nonlinear
  Sci. 7, 427-473 (1997).
- Castella, E., Jorba, A.: On the vertical families of two-dimensional tori near the triangular points of the bicircular problem.
  Celest. Mech. Dyn. Astron. 76(1), 35-54 (2000).
- Gabern, F., Jorba, A., Robutel, P.: On the accuracy of restricted three-body models for the Trojan motion. Discrete Contin. Dyn.
  Syst. Ser. B 11(4), 843-854 (2004).
- Lian, Y., Gomez, G., Masdemont, J., Tang, G.: A note on the dynamics around the L1,2 Lagrange points of the Earth-Moon system in a
  complete solar system model. Celest. Mech. Dyn. Astron. (2013). https://doi.org/10.1007/s10569-012-9459-2
- Szebehely, V.: Theory of Orbits. Academic Press, London (1967); Broucke, R.A.: Periodic orbits in the restricted three-body problem
  with Earth-Moon masses. JPL technical report (1968).
