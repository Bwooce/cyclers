# Digest: Rosales, Jorba & Jorba-Cusco (2021), "Families of Halo-like invariant tori around L2 in the Earth-Moon Bicircular Problem"

Celestial Mechanics and Dynamical Astronomy 133:16 (2021), DOI 10.1007/s10569-021-10012-0, 30 pages (received 23 February 2020,
revised 20 February 2021, accepted 27 February 2021, published online 27 March 2021). Topical collection "Toward the Moon and
Beyond".
Filed in the private paper corpus as
rosales-jorba-jorba-cusco-2021-families-halo-like-invariant-tori-l2-earth-moon-bicircular-problem-cmda-133-16-doi-10.1007-s10569-021-10012-0.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with printed page and section, equation, table or figure)
or INFERRED (our reading, a derivation from printed quantities, or a comparison with project code). Page numbers are the
journal's printed "Page n of 30". Numbers were taken from the PDF text layer and checked against the page images for Tables 1 to 3.
Values read off a plotted figure are marked "graph read" and are good to a few units in the last plotted digit only.
Context: tasks #891 (bicircular module `core/bcr4bp.py`, Sun sense corrected 2026-10-04), #892 and #893. Companion digests:
Rosales, Jorba and Jorba-Cusco (2023) on the quasi-bicircular model, and Jorba, Jorba-Cusco and Rosales (2020) on L1 in the BCP.

## 0. What the paper is

READ (abstract, p1): "The Bicircular Problem (BCP) is a periodic time dependent perturbation of the Earth-Moon Restricted
Three-Body Problem that includes the direct gravitational effect of the Sun. In this paper we use the BCP to study the existence
of Halo-like orbits around L2 in the Earth-Moon system taking into account the perturbation of the Sun. By means of computing
families of 2D invariant tori, we show that there are at least two different families of Halo-like quasi-periodic orbits around
L2."

READ: the paper has three tables only: Table 1 (the BCP constants), Table 2 (the six eigenvalues of the periodic orbit that
replaces L2 at epsilon = 1) and Table 3 (appendix: energy and rotation number of 21 RTBP Halo orbits). It prints no initial
condition of any periodic orbit or torus. The numbers a dynamical module can test are in Table 2 (the one that matters), Table 3,
and the rotation numbers and energies quoted in text and captions (section 2.2 below).

## 1. The model exactly as printed

### 1.1 Hamiltonian, frame and Sun (section 2, p3, eq. 1; Fig. 1, p4)

READ. Quoted setup (p3): "In the BCP, the dynamics of the Earth, Moon and Sun are simplified considering that the three bodies
orbit in the same plane. Also, it is considered that the Earth and the Moon follow a circular orbit around their barycenter (as in
the RTBP), and that B is orbiting around the S-E/M barycenter. Note that this model is not coherent, in the sense that the motion
of the three massive bodies is not described by the Newton's equations of motion... In our case, the RTBP is the Earth-Moon system
and the perturbing body is the Sun. It is then natural to use the units and reference frame of the Earth-Moon RTBP, so that the Sun
is moving around in a circular orbit." Momenta: "px = x' - y, py = y' + x, pz = z'" (dots in the paper).

Eq. (1), transcribed (lowercase variables as in the paper):

    H_BCP = (1/2)(px^2 + py^2 + pz^2) + y px - x py - (1 - mu)/r_PE - mu/r_PM - m_S/r_PS - (m_S/a_S^2)(y sin(theta) - x cos(theta))

with "r_PE^2 = (x - mu)^2 + y^2 + z^2, r_PM^2 = (x - mu + 1)^2 + y^2 + z^2, r_PS^2 = (x - x_S)^2 + (y - y_S)^2 + z^2, x_S = a_S
cos(theta), y_S = -a_S sin(theta), and theta = omega_S t with omega_S being the frequency of the Sun around the Earth-Moon
barycenter" (pp3 to 4). "Note that in this reference system the Sun moves around the origin in a circular motion (see Fig. 1)."
The derivation is referred to Gomez et al. (1993), Chapter 3.

READ: H_BCP(X, P_X, theta) = H_RTBP(X, P_X) + H_S(X, P_X, theta), with H_RTBP the autonomous part and H_S the Sun's contribution (p4).

INFERRED, from the formulas (same finding as for the 2023 and 2020 papers):
- Earth at x = mu, Moon at x = mu - 1. Fig. 1 (p4) shows M to the left of E, L1 between them, L2 beyond M (the sketch is schematic); S at the lower right with the angle theta drawn from the +x axis downward.
- The Sun is at (a_S cos(theta), -a_S sin(theta)), so it moves clockwise; at t = 0 it is on the +x axis, on the side away from the
  Moon. The indirect term `-(m_S/a_S^2)(y sin(theta) - x cos(theta))` equals `+(m_S/a_S^2)(x cos(theta) - y sin(theta))` and
  cancels the first-order term of `-m_S/r_PS` about the origin, so the printed Hamiltonian is internally consistent.
- In the project's frame (Earth at -mu, Moon at 1 - mu) the whole picture is rotated by pi: x, y, p_x, p_y and velocities change
  sign, the Sun starts at angle pi and is still clockwise.

### 1.2 Table 1 (p4): constants, as printed

| Symbol | Value |
| --- | --- |
| mu | 0.012150581623433 |
| m_s | 328900.55 |
| omega_s | 0.925195985518289 |
| a_s | 388.811143023351121 |

READ: a_s is printed with 15 decimals (beyond double precision; digits after the 16th significant one carry no information).
INFERRED: these agree to every printed digit with the constants in Table 2 of the 2023 paper (mu = 0.012150581623433, m_s =
328900.5499999991, omega_s = 0.925195985518289, a_s = 388.8111430233511) and with those used by `tests/core/test_bcr4bp_sun_sense.py`
(`_JJR_2020_*`, from the 2020 paper). Period of the Sun, `T = 2 pi/omega_S` (p5): 6.791193871923023 time units (our arithmetic).

### 1.3 Homotopy (section 2, p5, eq. 2)

READ. `H^epsilon = H_RTBP + epsilon H_S`, "with |epsilon| <= 1. Note that for epsilon = 0, H^0 = H_RTBP, and for epsilon = 1,
H^1 = H_BCP." The parameter multiplies "the mass of the Sun", and negative epsilon is used "which has no physical sense" (p6).
Period: "T = 2 pi/omega_S. More generally, a periodic orbit of the RTBP whose period is (p/q)T becomes a periodic orbit of period pT
once epsilon is set to be different from zero. This is a consequence of the Implicit Function Theorem."

## 2. Tables and printed numerical values, transcribed digit by digit

### 2.1 Table 2 (p7): eigenvalues of the monodromy matrix of the periodic orbit at epsilon = 1 (displayed in Fig. 3)

READ. Caption: "Eigenvalues of the monodromy matrix related to the periodic orbit displayed in Fig. 3". Footnote: "Due to the
Hamiltonian nature of the system, the other three eigenvalues are lambda_i^{-1}, i = 1, 2, 3. Also, note that due to the
non-autonomous character of the BCP Hamiltonian, there is no double eigenvalue 1".

| i | Re(lambda_i) | Im(lambda_i) |
| --- | --- | --- |
| 1 | 776607.104649077169597 | 0.0000000000000000 |
| 2 | 1.660211640235458 | 0.0000000000000000 |
| 3 | 0.865694004478591 | - 0.500573561636870 |

INFERRED (our arithmetic): abs(lambda_3) = 1.0000000000000 to the printed digits (so the pair lambda_3, conjugate is on the unit
circle, an elliptic direction; the printed digits are self-consistent), arg(lambda_3) = -0.5242611942478 rad, which over the
period T = 6.7912 is an exponent of -0.0772 per time unit modulo omega_S. lambda_1 and lambda_2 are real and larger than 1
(two saddles). The digits of lambda_1 beyond the 16th significant one (776607.1046490772) are not meaningful.
READ (p7): "After this point, the stability of the periodic orbits is of the type saddle x saddle x center until epsilon = 1. The
eigenvalues lambda_i, i = 1,...,6 of the monodromy matrix associated to the periodic orbit in the BCP are captured in Table 2.
Notice that, the final orbit has not the same stability character as L2."
READ (p9): "We recall from Table 2 that the largest eigenvalue of the monodromy matrix of the periodic orbit around L2 found in the
BCP is of order of 1e6."
Frame and epoch: the paper's frame (Earth at x = mu), Sun at angle theta = 0 at t = 0 (INFERRED: the stroboscopic map and the
"initial condition ... at time zero" of Fig. 2 use t = 0 with theta = omega_S t). The monodromy matrix is that of the period T of
the Sun, from the "initial condition corresponding to epsilon = 1 in Fig. 2" (Fig. 3 caption); eigenvalues do not depend on the
starting phase of the orbit along itself, but the orbit does depend on the Sun's phase at t = 0, so a project test must use the
orbit that the continuation (section 3 below) leads to, not an arbitrary phase.

### 2.2 Table 3 (appendix, p26): RTBP Halo orbits continued to the BCP to produce the Type I family

READ. Caption: "List of Halo orbits in the RTBP that have been continued to the BCP in order to obtain the Type I Halo family
displayed on Fig. 5". Note: "Each orbit can be identified by its energy or its rotation number. We display both". Rows in the
order printed:

| # | Rotation number | Energy |
| --- | --- | --- |
| 1 | 3.239814740891185 | -1.510315749412583 |
| 2 | 2.784894528517858 | -1.512100136308846 |
| 3 | 2.675226847819367 | -1.512638705372481 |
| 4 | 2.567114959430050 | -1.513218320582068 |
| 5 | 2.303428352712991 | -1.514862882292650 |
| 6 | 2.251789879351074 | -1.515228024256718 |
| 7 | 2.109091801535878 | -1.516320632450459 |
| 8 | 2.048773355923625 | -1.516822276667714 |
| 9 | 1.900124462997684 | -1.518170910364937 |
| 10 | 1.851244512820770 | -1.518652245225281 |
| 11 | 1.658983813333735 | -1.520751836889357 |
| 12 | 1.471783276849562 | -1.523162066038365 |
| 13 | 1.425750459444525 | -1.523818776410092 |
| 14 | 1.380018549762754 | -1.524498837931237 |
| 15 | 1.200041740490371 | -1.527471722488730 |
| 16 | 1.111784105690475 | -1.529124549475959 |
| 17 | 0.913023767640371 | -1.533412636021507 |
| 18 | 0.938620819394460 | -1.532811535496076 |
| 19 | 0.853672999241473 | -1.534868225111571 |
| 20 | 0.769787160604950 | -1.537084025294045 |
| 21 | 0.645906459334169 | -1.540737852533387 |

(Count: the printed list has 21 rows, not 20; the numbering is ours.) Observation, stated factually: rows 17 and 18 are out of
order against the otherwise monotone decrease of both columns (0.913... with energy -1.5334... precedes 0.938... with energy
-1.5328...); the pairing within each row is consistent (energy and rotation number both rise from row 17 to row 18), so only the
order of the two rows is exchanged, INFERRED. All rotation numbers are below pi. These are "Halo orbits in the RTBP"; their
energies are values of H_RTBP of section 1.1.
INFERRED (derivation): with the momenta of section 1.1, `H_RTBP = (1/2) v^2 - (1/2)(x^2 + y^2) - (1 - mu)/r_PE - mu/r_PM`, which is
`-C/2` for the Jacobi constant `C = x^2 + y^2 + 2(1 - mu)/r_PE + 2 mu/r_PM - v^2`; so the energies above correspond to C between
3.0206 (row 1: C = 3.020631498825166) and 3.0815 (row 21: 3.081475705066774). That range is the near-rectilinear and lunar-near
part of the L2 Halo family, not the small Halos next to L2 (C near 3.15).

### 2.3 Numerical values quoted in text and captions

READ (all with printed digits):
- Fig. 4 (p14) caption: "Transition from Halo orbit with energy -1.5244988379312372 in the RTBP (green) to a torus in the BCP (red). The
  torus is the dynamical equivalent (in the BCP) to the periodic orbit in the RTBP and has rotation number 1.3800185497627542."
  This is row 14 of Table 3 (energy and rotation number agree to every printed digit).
- Rotation numbers of tori, section 4.1 (p17): Type I Fig. 4, rho = 1.380018549762754; Type I Fig. 9, rho = 2.675226847819367
  ("close to the resonance value of rho = 6 pi/7 = 2.6927937..."; this is row 3 of Table 3); Type I Fig. 10, near the resonance,
  rho = 2.692464347819371; Type II Fig. 11, rho = 3.116137168026786; Type II Fig. 12, near the resonance rho = pi,
  rho = 3.130357871578353.
- Type II tori continued back to the RTBP (Fig. 13, p23): rho = 0.739476685309787 and rho = 0.858771705123796; their RTBP
  centre-manifold energy levels (of the reduced Hamiltonian, not the RTBP value): -1.56702372645620 and -1.56376188467649 as printed in
  the caption "Left, RTBP energy level: -1.56702372645620. Right, RTBP energy level: -1.56376188467649".
- Appendix (p25): planar tori H1 rho = 0.522687812628674 and H2 rho = 0.258684108104417 (Fig. 16); V1 rho = 0.651014628070470
  (Fig. 17); V2 rho = 0.585297052915989 (Fig. 18).
- Largest eigenvalue of the tori (p22): Type I "ranges from 2300 to 318600 (approximately)", Type II "from 23400 to 794260
  (approximately)". Family H1: "one real eigenvalue of the order of 1e6 (and its inverse)" (p16).
- Numerical method (pp10 to 12): r = 4 shooting sections, Newton tolerance 1e-6 for plots and 1e-10 for stability, Fourier degree
  N from 5 up to 252, step doubling of the continuation increment when fewer than 6 iterations are needed and halving when more.
- Resonance values on the axes of Fig. 5: pi/4, pi/3, pi/2, 2 pi/3, pi. Fig. 6a text: gap "corresponding to the 1:3 resonance, i.e.,
  rotation number close to 2 pi/3"; the plotted rotation-number range of that panel is 2.58 to 2.74 (graph read) and section 4.1
  puts the same gap at 6 pi/7. The caption and the text of section 4.1 do not agree on which resonance the gap is; the plotted range
  agrees with 6 pi/7.

## 3. The replacement of L2 in the BCP, and its stability (section 2, pp4 to 7)

READ (p4): "L2 is an equilibrium point of the RTBP. However, as opposed to the RTBP, the BCP is not an autonomous system... so that
the L2 point is not an equilibrium point anymore."
READ (pp5 to 7), the continuation in epsilon (multiple shooting, "the total number of sections used is four"): "Starting from L2, and
moving to the left the parameter increases until it hits a local maximum, and then decreases to cross the horizontal line and become
negative. The point on epsilon = 0 corresponds to a planar Lyapunov orbit whose period is half the one of the Sun, so we can see it
as a closed trajectory traveling twice around L2 in a single period of the Sun (i.e., a 1:2 resonant planar Lyapunov orbit). Moving
from L2 to the right, although the parameter epsilon becomes negative (which has no physical sense), it decreases until it hits a
turning point, and then increases to become positive and reach epsilon = 1, that is, the BCP. Again, in this case, the crossing point
with epsilon = 0 corresponds to the previous 1:2 resonant planar Lyapunov orbit... When epsilon > 0, this 1:2 resonant orbit becomes
a periodic orbit that nearly travels the same trajectory twice before closing the loop. This behavior is maintained until epsilon = 1,
see Fig. 3."
READ (p6): "the value of epsilon is not small enough and there is no natural dynamical substitute of the L2 point in the BCP. In other
words, there is no direct connection between L2 and a periodic orbit in the BCP."
READ (Fig. 2, p5; graph read, low confidence): the plot is x at t = 0 (horizontal, -1.18 to -1.10) against epsilon (0 to 1); the
curve is S shaped. Left end near x = -1.18 at epsilon about 0; local maximum epsilon about 0.07 near x = -1.17; passes epsilon = 0
near the Earth/Moon L2 point (-1.1557, drawn as a blue dot); local minimum epsilon slightly below 0 near x = -1.14; crosses 0 again
near x about -1.13; reaches epsilon = 1 at the right end near x about -1.10 (the plot edge). Colours: green = saddle x center x
center, red = saddle x saddle x center. The orbit at epsilon = 1 is plotted in Fig. 3 (x about -1.19 to -1.09, y about +-0.12, loops
around L2 twice, the Moon at x about -0.988 drawn for reference).
READ (stability sequence, p7): "It is observed that the periodic orbits alternate between the types saddle x center x center (green
regions) and saddle x saddle x center (red regions). Starting from L2, the linear stability is of the type saddle x center x
center. Moving to the left, epsilon increases and the periodic orbits keep this linear stability type until they hit the local
maximum. In this turning point, the linear stability becomes of the type saddle x saddle x center until another bifurcation point at
resonant 1:2 planar Lyapunov (epsilon = 0)... A similar pattern but with different sign for epsilon is observed when moving to the
right... Finally, this resonant planar Lyapunov orbit is continued until the last bifurcation point. This is a pitchfork
bifurcation, and it is where the 1:2 resonant (with the Sun) Halo orbit in the RTBP ends (this is shown in Andreu (1998), where a
bifurcation diagram and a suitable analysis are provided). The implication is that the 1:2 resonant Halo orbit in the RTBP does not
reach the BCP. As we will see in Sect. 4, this is not the case for all Halo orbits, and there is a dense set of Halo orbits that
survive the perturbation of the Sun as modeled in the BCP."
READ: values of epsilon at the bifurcations and at the turning points are not printed.
INFERRED: the period of the epsilon = 0 orbit is T/2 = 3.3956 time units. The planar normal frequency of POL2 in the 2023 paper,
1.86386291350378, corresponds to a period of 3.371, so the orbit is a modest-amplitude planar Lyapunov orbit, consistent with a
half-width of about 0.025 in x in Fig. 2 and Fig. 3 (graph read).

## 4. Families of tori and periodic orbits computed (sections 3 and 4, pp8 to 23)

READ. Method (section 3): invariant curves of the stroboscopic map F (flow over T = 2 pi/omega_S), `W(theta + rho) = F(W(theta))`,
`rho = 2 pi omega_1/omega_S` (eq. 3, p9), truncated Fourier series in theta (eq. 5), r-invariant curves with r = 4 for stability of
the unstable region (eqs. 7 to 12), stability by the generalised eigenproblem of Jorba (2001), continuation parameterised by rho.
The centre-manifold route was tried first and abandoned for L2 because its radius of convergence "was very small" (p8).
Families are Cantorian (gaps at resonances); passing a gap is done by going back to the RTBP by lowering the Sun's mass (p13).

READ. Six families (Fig. 5, p15; x is the x component of the invariant curve at theta = 0, vertical axis the rotation number):
- H1 and H2: planar quasi-periodic Lyapunov tori. H1 starts from the periodic orbit that replaces L2 (section 3); "most of the
  tori are hyperbolic"; always an eigenvalue 1 of multiplicity two plus a real eigenvalue of order 1e6 and its inverse; the remaining
  pair meets the unit circle at two points (Fig. 7), each giving a small interval of partially elliptic tori (Fig. 8). Graph read
  from the axes of Fig. 8: the two bifurcations are at x about -1.0862 and x about -1.0696. "Looking at the family H1 in Fig. 5,
  from left to right, the first bifurcation gives rise to the Type II Halo family, while the second one to the V1 family." The values
  of rho at the bifurcations are not printed. H2 was reached from V2 (below).
- Type I Halo: obtained by continuing RTBP Halo orbits (Table 3) from epsilon = 0 to 1 and then along rho. "Type I family is to be
  understood as the dynamical equivalent in the BCP to the classical Halo family of the RTBP." Stability: "mostly behave like their
  counterparts in the RTBP", saddle x center x center, with a real pair, a unit-modulus pair for almost every torus, and the double
  unit eigenvalue.
- Type II Halo: born from the first bifurcation of H1; "comes from a quasi-Halo orbit of the RTBP which has one frequency in resonance
  with the frequency of the Sun" (Fig. 13); stability saddle x saddle x center ("The other pair ... is also real and positive").
  "Type II family is less stable" (larger eigenvalue). Both Type I and Type II tori keep line of sight to the Earth for some members
  (Fig. 15, p24; the Moon is inside the loop as seen from the Earth).
- V1 and V2: vertical families that fall behind the Moon, "not Halo-like". V2 was found after passing a small resonance on V1 by the
  RTBP detour, and "eventually, the V2 branch met a planar quasi-periodic Lyapunov orbit of a new family, called H2". "A complete
  study of the H and V families is left for another work."
- The extent of each family in x is given only by the axes of Figs. 5 and 6 (graph read: x from about -1.10 to -1.02 for the
  families of Fig. 5, rotation number from about 0.2 to about 3.2); no table of it is printed.

## 5. Statements comparing the BCP and the QBCP, and objects absent from one

All READ, quoted:
- p2: "Focusing on the BCP and QBCP, it is interesting to mention that despite modeling the same system, there are qualitative
  differences between these two models around L2 (see Jorba-Cusco et al. 2018 for a discussion)."
- p2: "As opposed to the BCP, the QBCP is coherent and the motion of the three primaries is a solution of a three-body problem. The QBCP
  can also be formulated as a time-periodic perturbation of the RTBP. Hence, from a formulation point of view, the motion of the
  primaries is the only difference between the two models (see Andreu 1998 for the details on the QBCP derivation)."
- p3: "Notice that, due to the absence of a natural replacement of L2, the properties of some of these families change near the
  coordinates of the translunar point (which is no longer an equilibrium point in the BCP). In particular, we report the existence of a
  family of Halo-like orbits that does not come from the original Halo family in the RTBP."
- p7: "The comparison between the BCP and QBCP illustrates this phenomena. In the QBCP, L2 is replaced by a periodic orbit that is small
  in the sense that its maximal distance to L2 is of the order of 1e-6, and it has the same stability type of the L2 point. See (Andreu
  2002; Jorba-Cusco et al. 2018) and references therein for the details."
- p24: "the QBCP also models the direct gravitational effect of the Sun but, as of today, only the quasi-periodic counterparts of the
  Halo orbits (Type I family) have been computed (see Andreu 1998). The existence of Type II Halo-like orbits provides mission analysts
  with new potential candidates." (Later work, the 2023 paper, computes a QBCP Halo family and a QV family, not a Type II.)
- p7: the 1:2 resonant Halo orbit of the RTBP "does not reach the BCP"; the Type I tori "survive the perturbation of the Sun as modeled
  in the BCP".

## 6. Positive controls for the project

Frame recipe (INFERRED): the paper's frame has Earth at x = mu, Moon at x = mu - 1, Sun at (a_S cos(theta), -a_S sin(theta)) with
theta = omega_S t. The project's frame is rotated by pi: x, y, x', y', p_x, p_y change sign, z and p_z do not, Sun at angle pi at
t = 0 (`theta_sun0 = pi` in `core/bcr4bp.py`, as used in `tests/core/test_bcr4bp_sun_sense.py`), still clockwise. Eigenvalues and
rotation numbers are unchanged by the rotation. Already tested by the project: the 2020 paper's replacement orbit of L1
(`test_l1_replacement_orbit_matches_jorba_2020`) and the Sun's sense and phase conventions. Nothing in this paper is tested yet.

D1. Table 2, the periodic orbit that replaces L2 in the BCP (NEW, and the most valuable control for the bicircular module at L2).
Recipe: in `core/bcr4bp.py` (project frame) with the constants of section 1.2, take the RTBP planar Lyapunov orbit around L2 whose
period is exactly T/2 = pi/omega_S = 3.3956 time units (epsilon = 0, traversed twice in T), and continue it in the Sun's mass
(`mu_sun = epsilon * m_S`) from epsilon = 0 to 1 by multiple shooting over the Sun's period T = 6.7912 with the starting phase
theta_sun0 = pi, as the project did for L1 (4 or more segments; the multiplier is 7.8e5 over one period, so use at least four).
There are two crossings of epsilon = 0 in Fig. 2 (x about -1.18 and about -1.13 in the paper frame, i.e. x about +1.18 and +1.13 in
the project frame, graph read; the same orbit at two phases). Start from the right-hand one (paper frame x about -1.13) and increase
epsilon: on the right of Fig. 2 epsilon then rises to 1 (graph read). The left-hand crossing leads along a branch that folds at
epsilon about 0.07, so ordinary continuation in epsilon from there stalls; a pseudo-arclength continuation, as the paper used, is the
safe choice. Targets: real multipliers
776607.1046490772 and 1.660211640235458 (and their inverses), a unit-modulus pair 0.865694004478591 +- 0.500573561636870 i
(arg +-0.5242611942478 rad); the orbit loops twice around L2 in T (Fig. 3). Agreement the digits allow: 15 significant digits are
printed for the two real multipliers and the pair; the achievable agreement is set by conditioning (a 7.8e5 multiplier amplifies
errors of the solved initial state), so 1e-8 relative is a fair target (INFERRED). This control discriminates the Sun's sense and
phase strongly because the orbit is large (about 0.05 in x and 0.12 in y) and the Sun term is not a small perturbation here (the
paper: "the value of epsilon is not small enough"). The Sun's phase at t = 0 is part of the orbit's definition: the test should run
the continuation for theta_sun0 = pi, and, as a negative control, for the reversed sense (the multipliers should then differ).
If the epsilon = 1 orbit reached is a different one (different loop count, different stability type), compare the stability type
first: the printed type is saddle x saddle x center.

D2. Table 3, RTBP Halo energy against rotation number (NEW; independent of the Sun's sense; tests only the CR3BP, the energy
convention and the Sun's period). Convention (INFERRED): `C = -2 H_RTBP` (section 2.2); for a Halo orbit of RTBP period P the
rotation number in the stroboscopic map of period T is `rho = 2 pi T/P` reduced modulo 2 pi and possibly mirrored (rho and -rho
are the same curve up to orientation; the paper's tabulated values are all below pi). Test: for each row, find the L2 Halo orbit of
Jacobi constant C = -2 * energy (rows 1, 14 and 21: C = 3.020631498825166, 3.048997675862474, 3.081475705066774), compute P, and
check that some integer k and sign give `rho = +-(2 pi T/P) + 2 pi k`. With T = 6.7912 and Halo periods of a few time units or less
(from general knowledge of the family, not from this paper) the integer k is several, so it cannot be fixed in advance; the test is
the pattern over the 21 rows (a single smooth monotone P(C) must reproduce all rotation numbers with k changing only where rho wraps).
Agreement allowed by the digits: 15 digits printed for both columns; if the hypothesis is right, the match should reach 1e-9 or
better with the paper's T = 6.791193871923023, otherwise the hypothesis (a torus rather than RTBP rotation number) is wrong and the
result is the first thing to report. Status: hypothesis, not yet checked.
Consistency evidence (INFERRED): the rotation number of the QBCP curve ICQ1 in the 2023 paper, 3.239814740891185, equals row 1 to all
printed digits, although the QBCP and the BCP constants differ in the 12th significant digit of omega_S; that is what one expects if
rho was held fixed while the Sun's mass was switched on, so the printed rho is a property of the RTBP Halo orbit and the period
used by the authors.

D3. Cross-check between Fig. 4 and Table 3 (READ, no computation): energy -1.5244988379312372 and rho = 1.3800185497627542 (Fig. 4)
equal row 14 of Table 3. Not a control for a module; confirms that Table 3's rotation number is the one of the torus continued from
that Halo orbit.

D4. Ranges of the largest eigenvalue of tori (Type I 2300 to 318600, Type II 23400 to 794260) and the 1e6 for H1: require tori;
not realistic for the project without an invariant-curve solver. The multiplier of D1 (7.77e5) lies inside the Type II range, as one
would expect of a nearby object (INFERRED observation, no claim).

Not testable here: bifurcation values (epsilon at the turning points, rho at the H1 bifurcations: not printed), the Fig. 5 plot,
the Type II families' initial states.

## 7. Corrections to the existing second-hand digest (`2026-06-14-andreu-quasi-bicircular-digest.md`)

That digest does not cite this paper. Its statement that halo initial conditions are not tabulated in the open literature it used is
consistent with this paper: Table 3 gives energies and rotation numbers, not states. The Hamiltonian, symmetry, "O(epsilon^2)",
attribution and canonical-momentum issues found by reading the 2023 paper are listed in section 7 of that paper's digest. One point
relevant here: that digest treats the QBCP constants as the BCP's; this paper's Table 1 holds the BCP set (mu = 0.012150581623433,
m_s = 328900.55, omega_s = 0.925195985518289), which is the BCP set, whereas the constants quoted in the
`core/bcr4bp.py` header (lines 11 to 14) and named `_ANDREU_*` (lines 101 to 110) are the QBCP set of Table 3 of the 2023 paper (mu =
0.0121505816, m_S = 328900.5423094043, omega_S = 0.925195985520347, a_S = 388.8111430233511). INFERRED from reading those lines; the
BCR4BP test file uses the BCP set, so this concerns defaults and documentation only.

## 8. References cited by the paper that bear on this project (full citations as printed)

Bicircular model and earlier work in it:
- Huang, S.: Very restricted four-body problem. Technical note TN D-501, Goddard Space Flight Center, NASA (1960).
  https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19890068606.pdf
- Cronin, J., Richards, P., Russell, L.: Some periodic solutions of a four-body problem. Icarus 3, 423-428 (1964).
- Gomez, G., Jorba, A., Masdemont, J., Simo, C.: Study of Poincare maps for orbits near Lagrangian points. ESOC contract 9711/91/D/IM(SC),
  final report, European Space Agency, 1993. Reprinted as Dynamics and mission design near libration points. Vol. IV, Advanced methods for
  triangular points, volume 5 of World Scientific Monograph Series in Mathematics (2001). (Derivation of the BCP equations, Chapter 3.)
- Gomez, G., Llibre, J., Martinez, R., Simo, C.: Study on orbits near the triangular libration points in the perturbed Restricted
  Three-Body Problem. ESOC contract 6139/84/D/JS(SC), final report, European Space Agency, 1987. Reprinted as Dynamics and mission design
  near libration points. Vol. II, Fundamentals: the case of triangular libration points, volume 3 (2001).
- Simo, C., Gomez, G., Jorba, A., Masdemont, J.: The bicircular model near the triangular libration points of the RTBP. In: Roy, A.,
  Steves, B. (eds.) From Newton to Chaos, pp. 343-370. Plenum Press, New York (1995).
- Castella, E., Jorba, A.: On the vertical families of two-dimensional tori near the triangular points of the bicircular problem. Celestial
  Mech. 76(1), 35-54 (2000). https://doi.org/10.1023/A:1008321605028
- Castella, E.: Sobre la dinamica prop dels punts de Lagrange del sistema Terra-Lluna. PhD thesis, Univ. Barcelona (2003).
- Jorba, A., Jorba-Cusco, M., Rosales, J.J.: The vicinity of the Earth-Moon L1 point in the Bicircular Problem. Celestial Mech. 132(2),
  11 (2020). https://doi.org/10.1007/s10569-019-9940-2
- Jorba, A., Nicolas, B.: Transport and invariant manifolds near L3 in the Earth-Moon Bicircular model. Commun. Nonlinear Sci. Numer.
  Simul. 89, 105327 (2020). https://doi.org/10.1016/j.cnsns.2020.105327
- Jorba-Cusco, M., Farres, A., Jorba, A.: Two periodic models for the Earth-Moon system. Front. Appl. Math. Stat. 4, 32 (2018).
  https://doi.org/10.3389/fams.2018.00032
- Scheeres, D.J.: The restricted Hill four-body problem with applications to the Earth-Moon-Sun system. Celestial Mech. 70(2), 75-98
  (1998). https://doi.org/10.1023/A:1026498608950
- Rosales, J., Jorba, A., Jorba-Cusco, M.: The effect of the Sun on direct transfers from Earth to translunar Halo orbits. In preparation
  (2020). (Published as Rosales, Jorba and Jorba-Cusco, "Transfers from the Earth to L2 Halo orbits in the Earth-Moon bicircular problem",
  Celest. Mech. Dyn. Astron. 133:12 (2021), per the 2023 paper's reference list; INFERRED link.)

Quasi-bicircular coefficients:
- Andreu, M.A.: The Quasi-Bicircular Problem. PhD thesis, Univ. Barcelona (1998).
- Andreu, M.A.: Dynamics in the center manifold around L2 in the Quasi-Bicircular Problem. Celestial Mech. 84(2), 105-133 (2002).
  https://doi.org/10.1023/A:1019979414586
- Le Bihan, B., Masdemont, J., Gomez, G., Lizy-Destrez, S.: Invariant manifolds of a non-autonomous quasi-bicircular problem computed
  via the parameterization method. Nonlinearity 30, 3040-3075 (2017).

Periodic orbits and Halo families in the RTBP, methods:
- Breakwell, J., Brown, J.: The 'Halo' family of 3-dimensional periodic orbits in the Earth-Moon restricted 3-body problem. Celestial
  Mech. 20(4), 389-404 (1979).
- Gomez, G., Mondelo, J.: The dynamics around the collinear equilibrium points of the RTBP. Phys. D 157(4), 283-321 (2001). (Multiple
  shooting for periodic orbits.)
- Jorba, A., Masdemont, J.: Dynamics in the center manifold of the collinear points of the restricted three body problem. Physica D 132,
  189-213 (1999). https://doi.org/10.1016/S0167-2789(99)00042-1
- Jorba, A.: A methodology for the numerical computation of normal forms, centre manifolds and first integrals of Hamiltonian systems.
  Exp. Math. 8(2), 155-195 (1999).
- Jorba, A.: Numerical computation of the normal behaviour of invariant curves of n-dimensional maps. Nonlinearity 14(5), 943-976 (2001).
  https://doi.org/10.1088/0951-7715/14/5/303
- Jorba, A., Villanueva, J.: On the persistence of lower dimensional invariant tori under quasi-periodic perturbations. J. Nonlinear Sci.
  7(5), 427-473 (1997). https://doi.org/10.1007/s003329900036
- Jorba, A., Olmedo, E.: On the computation of reducible invariant tori on a parallel computer. SIAM J. Appl. Dyn. Syst. 8(4), 1382-1404 (2009).
- Gonzalez, J., Mireles James, J.: High-order parameterization of stable/unstable manifolds for long periodic orbits of maps. SIAM J. Appl.
  Dyn. Syst. 16 (2016). https://doi.org/10.1137/16M1090041
- Gomez, G., Llibre, J., Martinez, R., Simo, C.: Station keeping of libration point orbits. ESOC contract 5648/83/D/JS(SC), final report,
  European Space Agency (1985). Reprinted as Dynamics and mission design near libration points. Vol. I (2001).
- Gomez, G., Jorba, A., Masdemont, J., Simo, C.: Study refinement of semi-analytical Halo orbit theory. ESOC contract 8625/89/D/MD(SC),
  final report, European Space Agency (1991). Reprinted as Dynamics and mission design near libration points. Vol. III (2001).
- Farres, A., Jorba, A.: On the high order approximation of the centre manifold for ODEs. Discrete Contin. Dyn. Syst. Ser. B 14(3), 977-1000
  (2010). https://doi.org/10.3934/dcdsb.2010.14.977
- Duarte, G.: On the Dynamics Around the Collinear Points in the Sun-Jupiter System. PhD thesis, Univ. Barcelona (2020).
- Stoer, J., Bulirsch, R.: Introduction to Numerical Analysis, Texts in Applied Mathematics 12, Springer (2002); Seydel, R.: Practical
  Bifurcation and Stability Analysis, Springer (2009).
