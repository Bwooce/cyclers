# Digest: Boudad, Howell and Davis (2020), "Dynamics of synodic resonant near rectilinear halo orbits in the bicircular four-body problem"

K.K. Boudad, K.C. Howell (Purdue University), D.C. Davis (a.i. solutions), Advances in Space Research 66:2194-2214 (2020),
DOI 10.1016/j.asr.2020.07.044, 21 pages. Received 23 June 2020, accepted 30 July 2020.
Filed in the private paper corpus as
boudad-howell-davis-2020-dynamics-synodic-resonant-NRHO-bicircular-four-body-problem-asr-66-2194-doi-10.1016-j.asr.2020.07.044.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with page, section, equation, figure or table number) or
INFERRED (our reading or arithmetic). Page numbers are the journal's printed numbers (2194 to 2214). Values were read from page
images; values that exist only in a figure are marked "read off the figure" and are approximate.
Context: tasks #884 and #891. The paper is the published BCR4BP continuation in Sun mass (epsilon from 0 to 1) of
commensurate-period three-body orbits, with several counterparts per orbit and branches that stop short of full Sun mass.

## 0. What the paper is, and a limit on its use as a numerical control

READ (abstract, p2194): NRHOs of the Earth-Moon CR3BP that are synodic resonant with the Sun are transitioned to the BCR4BP;
geometry, perilune and apolune radii are "generally preserved" (conclusion, p2214); stability, eclipse avoidance and energy are
examined.

IMPORTANT, READ: the paper prints NO initial conditions, NO state vector, NO exact period, and NO table of perilune radii or Sun
angles. Its only tables are Table 1 (eigenvalues of the instantaneous equilibrium points in the Earth-Moon frame), Table 2 (the
same in the Sun-B1 frame) and Table 3 (Lyapunov exponents of four NRHOs). The numerical constants of the model (mass ratio, Sun
mass, Sun distance, units) are not printed either; only the Sun's angular rate, 0.9253 (four digits), appears (p2196). The
other printed numbers are in prose and figures (section 5 below). So there is no closure test to run against this paper;
section 6 below lists what can be checked instead. Constants are presumably in Boudad (2018, M.S. thesis, Purdue) and Boudad,
Davis and Howell (2019, SciTech Forum); neither is held (INFERRED).

## 1. The model as printed

### 1.1 CR3BP (section 2.1, p2195, eqs. 1 to 3) - READ

    xddot = 2 ydot + dU*/dx,  yddot = -2 xdot + dU*/dy,  zddot = dU*/dz                 (1)
    U* = 1/2 (x^2 + y^2) + mu / r_(e-sc) + (1 - mu) / r_(m-sc)                           (2)
    C = 2 U* - sqrt(xdot^2 + ydot^2 + zdot^2)                                            (3)

with `mu = m_m / (m_e + m_m)` ("is the mass parameter for the Earth-Moon CR3BP system"), characteristic length the Earth-Moon
distance, mass the sum of the primaries, time such that the gravitational constant is 1.
Slip to note, factual: with `mu` defined as the Moon's mass fraction, eq. (2) as printed puts `mu` on the Earth distance
`r_(e-sc)` and `1 - mu` on the Moon distance; the usual weights are the other way round, and eq. (10) of the same paper (Sun-B1
frame) has `(1 - mu)` with `r_(e-sc)` and `mu` with `r_(m-sc)`. INFERRED: a typesetting transposition in eq. (2); it has no
effect on the paper's results.
Primary positions are not printed in this section (INFERRED: Earth at -mu and Moon at 1 - mu, the standard placement; Fig. 1
shows the Earth at the left and the Moon at the right of the barycenter B1).

### 1.2 BCR4BP, Earth-Moon rotating frame (section 2.2, p2195-2196, eqs. 4 to 7) - READ

Quote (p2195): "The BCR4BP is not a coherent model: the perturbing acceleration from the Sun does not influence the motion of
the Earth and the Moon, thus, the motion of the Moon is not a solution to the Sun-Earth CR3BP. Coherent bicircular models have
been investigated previously (Andreu, 1998), but are not necessary in this analysis."

    xddot = 2 ydot + dY*/dx,  yddot = -2 xdot + dY*/dy,  zddot = dY*/dz                 (4)
    Y* = U* + mu_s / r_(s-sc) - mu_s / a_s^3 (x_s x + y_s y + z_s z)                     (5)

"mu_s = m_s / (m_e + m_m) is the nondimensional mass of the Sun and a_s = r_s / r_(e-m) is the nondimensional distance between
the Earth-Moon barycenter and the Sun. Note that this adaptation of the BCR4BP assumes the Sun moves in the Earth-Moon plane of
motion." The Sun position in the rotating frame (eq. 6, p2196):

    [x_s; y_s; z_s] = a_s [cos(theta); sin(theta); 0] = a_s [cos(-omega t + theta_0); sin(-omega t + theta_0); 0]       (6)

"where the Sun angle theta is measured from the rotating x axis to the Sun position vector as defined in Fig. 1, and omega =
0.9253 is the magnitude of the nondimensional angular velocity of the Sun as viewed in the Earth-Moon rotating frame. This
angular velocity is computed as the difference between the nondimensional mean motion of the Sun in the inertial frame centered
at the Earth-Moon barycenter, that is, n_S = sqrt((1 + mu_s) / a_s^3), and the nondimensional mean motion of the Earth-Moon
system with respect to the same observer, n, that is, the value one."

Sun's sense, READ, one sentence for the hand-back: the Sun is at `a_s (cos(theta), sin(theta), 0)` with `theta = -omega t +
theta_0` and `omega = 0.9253 > 0`, so the angle DECREASES with time and the Sun moves clockwise in the Earth-Moon rotating frame
(the same sense as the corrected project module, which uses `theta_S0 - omega_S t`). INFERRED: `1 - n_S` with
`n_S = sqrt((1 + mu_s)/a_s^3)` is about 1 - 0.0748 = 0.9252 for project values, consistent with the printed 0.9253.

Indirect term (READ, eq. 5): `- mu_s / a_s^3 (x_s x + y_s y + z_s z)`; its gradient is `- mu_s r_s / a_s^3`, the usual indirect
acceleration for a Sun at (x_s, y_s, z_s); consistent with the Sun position in eq. (6).

Hamiltonian-like quantity (eq. 7): `H(theta) = 2 Y* - sqrt(xdot^2 + ydot^2 + zdot^2)`, "a scaled version of the Hamiltonian
value ... consistent with the Jacobi constant in the CR3BP". Because time appears explicitly, there is no integral of motion
(p2196).

Continuation parameter (section 3.2, p2201, eq. 15) - READ:

    Gamma* = U* + eps mu_s / r_(s-sc) - eps mu_s / a_s^3 (x_s x + y_s y + z_s z)         (15)

"such that epsilon is the coefficient that scales the mass of the Sun, i.e., mu_s. Note that for epsilon = 0, the updated
pseudo-potential Gamma* is equivalent to the CR3BP pseudo-potential, introduced in Eq. (2). Conversely, for epsilon = 1, Eq. (15)
is equivalent to the BCR4BP pseudo-potential, in Eq. (5)." Epsilon scales the direct and the indirect term alike. The plots
show "Sun mass [%]" (100 times epsilon) on the vertical axis.

Sun angle at the reference epoch: `theta_0` is the Sun angle at t = 0 (eq. 6). In the solar mass exclusion plots (section 4.5)
the angle is quoted at each perilune, in the Earth-Moon frame: "Sun angle at perilune [deg]", `theta = -omega t + theta_0`
(p2207). Reference epoch (section 4.5, p2210, READ): "Periodic orbits in the BCR4BP model are defined for a specific epoch."
The epoch is then characterised in the plots by the Sun angle at each perilune (INFERRED from Figs. 24, 25, 27, 28).

### 1.3 BCR4BP, Sun-B1 rotating frame (section 2.2, p2196-2198, eqs. 9 to 14) - READ

Sun and B1 (Earth-Moon barycenter) fixed on the x axis, the barycenter B2 of the Sun-B1 pair at the origin; Earth and Moon on
circles about B1.

    xddot = 2 ydot + dY*/dx, ... (9)
    Y* = 1/2 (x^2 + y^2) + (1 - 1/(mu_s+1))/r_(s-sc) + (1/(mu_s+1)) (1 - mu)/r_(e-sc) + (1/(mu_s+1)) mu/r_(m-sc)       (10)

Moon angle (printed): "theta_ = pi - theta = omega_ t + theta_0_ is the Moon angle, and omega_ = omega / (1 - omega) is the
nondimensional angular rate of the Earth and the Moon in their motion around their common barycenter B1." The Moon moves
counter-clockwise in the Sun-B1 frame; the Earth-Moon-frame Sun angle theta and the Moon angle are linked by
`theta_ = pi - theta`. The printed Earth and Moon positions (eq. 11) were not read reliably (as rendered, the y entry of the
first vector appears to repeat the x entry); they are not relied on. The two formulations "describe the same dynamical system"
and a periodic orbit in one is periodic in the other (p2198, p2202).

## 2. Synodic resonance, family selection and the number of counterparts

### 2.1 Definition (sections 3, 3.1; p2199-2200) - READ

"Periodic orbits in the CR3BP are solutions that precisely repeat in all six position and velocity states ... Periodic solutions
in the BCR4BP require an additional condition, i.e., the Sun location must also be commensurate with the periodic cycle over
which the states repeat. ... Thus, the period for any periodic orbit in the BCR4BP is a multiple of the Earth-Moon-Sun period,
that is, approximately 29.5 days. As a consequence, all periodic solutions in the BCR4BP are isolated synodic resonant orbits."
(p2199-2200.) "Across a family of orbits, the ratio of the orbital period to the synodic period, that is, P_orbital / P_synodic,
is computed. Synodic resonant orbits are characterized by a rational quotient, denoted the resonance ratio. This ratio is
represented as P : Q, where P is the number of orbital periods and Q is equal to the number of lunar synodic periods." (p2200;
the plots use T_synodic / Period.) The paper also says: "Recall that the resonance ratio is inversely proportional to the
orbital period."
Examples (p2200 and p2203): 3:1 means the orbital period is exactly one third of the synodic period; 9:2 means nine revolutions
in exactly two synodic periods (59 days, p2203). The BCR4BP counterparts are called "T_syn-periodic" for P:1 and "2 T_syn-periodic"
for 9:2.

### 2.2 How members are selected (section 3.1, p2200; section 4.1, p2202) - READ

Initial guess: stack P revolutions of the CR3BP orbit to fill Q synodic periods; discretise into patch points; parallel shooting
(Keller 1976) in the BCR4BP with epsilon stepped from 0 to 1 (Fig. 9: 25 equal steps). Selection by simple resonance ratios P:Q
from the L2 halo family, Fig. 8b: 3:1, 4:1, 9:2 NRHOs and 5:1 NRHO+ (NRHO+ = the part of the family beyond the NRHO boundary,
defined by the first stability change, with perilune above the lunar surface for some members; section 2.3, p2199). 6:1 NRHO+
(perilune about 500 km, below the lunar radius of about 1700 km) is excluded.

### 2.3 Sun phase and the number of distinct counterparts (section 4.5, p2210-2213, Figs. 24 to 29) - READ

The Sun phase at the start is not fixed by symmetry; it is scanned (the "solar mass exclusion plots"). Quote (p2210): "To explore
the impact of this initial epoch selection on the final BCR4BP synodic resonant orbit, solar mass exclusion plots are
exploited. Solar mass exclusion plots relate the initial epoch to the convergence of a BCR4BP periodic orbit." Horizontal axis:
Sun angle at each perilune from 0 to 360 degrees; vertical axis: percentage of Sun mass reached. "blue dots correspond to initial
conditions that reach 100% of the actual solar mass, that is, initial conditions that successfully yield a periodic solution in
the Earth-Moon-Sun BCR4BP. Red dots correspond to solutions that only converge for assumed masses smaller than 100% of the true
solar mass; that is, red dots are associated with initial Sun angles that do not produce a periodic solution in the Earth-Moon-Sun
BCR4BP." (p2211).

Counts:
- 3:1 (Fig. 24, p2211): "Three nearly vertical blue lines are apparent ... Note that the three lines, associated with Sun
  angles approximately equal to 0, 120 and 240 degrees, all correspond to the same orbit. The T_syn-periodic orbit in the BCR4BP
  corresponding to the 3:1 NRHO in the CR3BP includes three lobes around the Moon, and, thus, three distinct perilunes ... Each
  vertical line thus corresponds to the Sun angle at one perilune passage. Note that the lines are separated by approximately 120
  degrees, or, equivalently, approximately one third of the synodic period." One counterpart. "Most of the initial conditions,
  marked in blue, immediately converge to the BCR4BP periodic orbit when the Sun is introduced, i.e., for an assumed mass greater
  than or equal to 1% of the actual solar mass. The remaining initial conditions, denoted by red dots, yield a periodic orbit
  offset by 60 degrees. Note that the red lines marking the Sun angles at perilune approximately equal to 60, 180, and 300
  degrees only exist for assumed solar mass values between 1% and 17% of the true solar mass; convergence for assumed mass greater
  than 17% is not achieved. The CR3BP synodic resonant 3:1 NRHO has only one BCR4BP counterpart, plotted in Fig. 12(b)."
- 4:1 (Fig. 25-26, p2211-2212): "Two T_syn-periodic counterparts in the Earth-Moon-Sun BCR4BP of the 4:1 synodic resonant CR3BP
  NRHO are obtained. One counterpart, labeled orbit A ... presents the following Sun angle at perilune combination: 45, 135, 225
  and 315 degrees. The second analog's apolunes occur for Sun angles equal to 0, 90, 180 and 270 degrees, and is indicated by the
  purple dot in Fig. 25." Eight blue lines separated by about 45 degrees. (INFERRED: the text says "apolunes" for B where the
  figure caption of Fig. 25 refers to Sun angle at perilune; read it as the Sun angles at perilune of B being offset 45 degrees
  from A.) Geometry: A four perilunes at about the same radius (about 5940 km), B two groups (about 4970 km and about 7100
  km).
- 9:2 (Figs. 27-29, p2212-2213): "The number of BCR4BP counterparts to a synodic resonant CR3BP orbit is not limited to two."
  Nine lobes give nine perilune Sun angles; repeating vertical pattern every 40 degrees (two consecutive perilunes are about 80
  degrees apart, but since the nine perilunes take two synodic periods "the smallest angular difference between two non-consecutive
  perilunes is 40 degrees"). Zoom 0 to 40 degrees (Fig. 28); four converged counterparts shown (Fig. 29), at Sun angles at perilune
  read off Fig. 28 of about 0, 5, 16 and 33 degrees. Quote (p2213): "Numerous counterparts of the CR3BP 9:2 NRHO are available in
  the BCR4BP."
- Why several counterparts: no symmetry argument is made. The only explanations are geometric (the Sun's tidal pull stretches
  alternate lobes in opposite directions; p2203: "Tidal effects from the Sun introduce similar effects in two opposite quadrants
  representing the Sun angle theta between 0 and 90 degrees, and the quadrants between -180 and -90 degrees ... Conversely, the
  tidal effects are also similar for Sun angle quadrants defined between 0 and -90 degrees, as well as between 90 and 180 degrees.
  Thus, the first and the third revolutions of the BCR4BP T_syn-periodic orbit are 'pulled' by the Sun's gravitational
  acceleration in one direction; the second and the fourth revolutions are impacted by the solar acceleration in the opposite
  direction") and numerical (the counterparts are found by scanning the starting Sun angle). Concluding remarks (p2214): "While
  certain synodic resonant CR3BP NRHOs present a unique BCR4BP counterpart, others offer two or more counterparts when
  transitioned to the BCR4BP."

## 3. Orbit families treated

READ and plain: Earth-Moon CR3BP near rectilinear halo orbits (the L2 halo family, NRHO and NRHO+ subsets) only: 3:1, 4:1, 9:2
NRHOs and the 5:1 NRHO+. In section 3.3 of the paper (p2201-2202) the L1 and L2 1:1 Lyapunov orbits are transitioned as an illustration
(Figs. 10-11). NO Earth-Moon cycler, no transfer orbit and no resonant orbit with lunar flybys (in the cycler or free-return
sense) is treated; the NRHOs have close lunar passes of 1500 to 15000 km perilune but are Moon-centred orbits, not cyclers.
INFERRED: the project's expectation is confirmed.

Orbit data stated in prose (READ; km as printed; these are the only numbers describing the orbits):

| Orbit | CR3BP period | CR3BP perilune | CR3BP apolune | BCR4BP counterpart |
| --- | --- | --- | --- | --- |
| 3:1 NRHO (p2202) | "exactly one third of the synodic period, that is, approximately 9.79 days" | about 15,000 km | about 84,500 km | T_syn-periodic; three lobes; perilunes 14,300 to 15,200 km; apolunes 82,300 to 89,400 km |
| 4:1 NRHO (p2203) | one quarter of the synodic resonance, "7.34 days" | about 5,600 km | "about 75,335 km" | T_syn-periodic; four lobes; apolunes 74,800 to 75,500 km; all four perilunes about 5,900 km (A: about 5,940 km; B: about 4,970 km pair and about 7,100 km pair, p2212) |
| 9:2 NRHO (p2203) | "6.53 days"; nine revolutions in two synodic periods, 59 days | about 3,100 km | about 71,000 km | 2 T_syn-periodic; nine lobes; perilunes 3,100 to 3,900 km; apolunes 69,900 to 71,700 km |
| 5:1 NRHO+ (p2204) | about one fifth of the synodic period | 1,650 km | 66,000 km | T_syn-periodic; five lobes; perilunes 1,500 to 3,000 km; apolunes 64,500 to 68,200 km |

INFERRED arithmetic note, respectfully: the printed day-values imply a synodic period of about 29.4 days (3 x 9.79 = 29.37;
4 x 7.34 = 29.36; 4.5 x 6.53 = 29.39), whereas the paper elsewhere uses about 29.5 days and the mean synodic month is 29.53 days
(which would give 9.84, 7.38 and 6.56 days). The three printed day-values are each about 0.5 percent short of that. Do not use
them as quantitative targets; use the ratio.
Section 4.4 also states: the NRHO subset has perilune radii "between 1,800 and 17,300 km"; NRHO+ orbits have perilune "less than
1,800 km".

## 4. What happens along the continuation (stability, folds, bifurcations)

READ, section 4.4 and 4.5 (p2207-2211), Fig. 23.
- Stability definition: Lyapunov exponent `phi_i = Re(ln(lambda_i) / T)` of the monodromy matrix, T the period (eq. 14, p2198); two
  of six are always zero; stable if all six are zero. "One advantage ... is that they are not influenced by the orbital period."
- 4:1 and 9:2: unstable in the CR3BP with a small unstable mode, and "generally consistent" in the BCR4BP (Table 3, section 5).
  "Both orbits are linearly unstable, but the magnitude of the Lyapunov exponent associated with the unstable mode is small."
  Lyapunov-exponent curves in Fig. 23b, c are nearly vertical lines over the continuation, i.e. no bifurcation along epsilon.
- 3:1: linearly stable in the CR3BP, linearly unstable in the BCR4BP. "By examination of the evolution of the Lyapunov exponents
  over the continuation in Sun mass ... a bifurcation occurs in the continuation process when the assumed mass is approximately
  15% of the true solar mass, and the T_syn-periodic in the BCR4BP that corresponds to the 3:1 NRHO in the CR3BP is unstable in the
  linear sense." (p2209-2210, Fig. 23a: two lines leave phi = 0 near 15 to 17 percent and reach about plus or minus 0.62 at
  100 percent, read off the figure.) Also (p2209-2210): "No discontinuity is observed but the stability
  properties for the periodic orbits are different between the CR3BP and BCR4BP, a bifurcation occurs along the evolution of the
  family in Sun mass."
- 5:1 NRHO+: stable in the CR3BP, weakly unstable in the BCR4BP (0.0322); "a bifurcation occurs at an assumed mass around 90% of
  the true solar mass" (p2210, Fig. 23d: four of six exponents stay at zero, two leave at about 90 percent).
- Caution stated by the authors (p2209): "The stability differences for certain synodic resonant periodic orbits signal a need
  for further verification that the BCR4BP periodic solutions are the true analogs of the CR3BP synodic resonant NRHOs."
  The discontinuities-hypothesis in the same passage is hedged ("suggest that the type of solutions ... might have also shifted
  during the corrections and continuation process and that the periodic orbit is not the true counterpart of the CR3BP orbit"),
  and is followed by "No discontinuity is observed". Conclusion on the check (p2210): "it is determined that the NRHOs, as constructed
  in the BCR4BP, are the counterparts of the CR3BP orbits, but the linear stability characteristics are not necessarily the same
  in both models."
- Period-multiplying bifurcations: none discussed (INFERRED from reading all 21 pages; "bifurcation" is used only for the
  eigenvalue leaving the unit circle in Fig. 23 and for the halo-family bifurcation from the Lyapunov family). Energy-related
  observation (section 4.3, p2207-2208): the Hamiltonian value along the 3:1 BCR4BP orbit "dips below the line associated with
  the Hamiltonian value of E_3(theta) for certain values of theta" (p2208).

Folds and branches that stop short of full Sun mass: the paper does not use the word "fold". The relevant READ statements are
about the solar mass exclusion plots.
- 3:1 (p2211): the offset-by-60-degrees family (perilune Sun angles 60, 180, 300 degrees) "only exist[s] for assumed solar mass
  values between 1% and 17% of the true solar mass; convergence for assumed mass greater than 17% is not achieved."
- 9:2 (p2213): "the orange and red lines that experience similar evolution to the blue lines as the assumed mass of the Sun is
  increased, but do not reach 100% of the true solar mass. Further experience is required to establish whether adjustments to the
  numerical corrections and continuation scheme could extend the convergence or whether these sets of initial conditions do not
  yield a periodic orbit in the Earth-Moon-Sun BCR4BP."
- INFERRED: these are branches (fixed starting Sun angle) that converge from epsilon = 0 up to an intermediate Sun mass and then
  fail; the paper uses a natural-parameter-style continuation (epsilon steps, parallel shooting), not pseudo-arclength, so a fold
  in epsilon appears only as lost convergence and the branch is not followed round the fold. It does not report whether the
  branch turns back to epsilon = 0 (that is what Oshima 2022 shows with pseudo-arclength continuation, see the other digest).
  Note too the observation that the family of Sun-angle starting values moves with epsilon: "Some perilunes are shifted toward
  smaller Sun angles as the mass of the Sun is increased while others are shifted to larger Sun angles, creating the woven pattern
  in Fig. 28" (p2213), so the Sun phase of an equivalent is not constant along the continuation.

## 5. Tables, transcribed digit by digit

No table of initial conditions, periods, perilune radii or Sun angles exists in the paper. All three tables are transcribed.
All digits legible; none marked unclear.

Table 1 (p2197), "Eigenvalues associated with the equilibrium points in the CR3BP and in the BCR4BP", Earth-Moon rotating frame,
nondimensional (eigenvalues of the 3 x 3 matrix A of eq. 8; the printed pairs are the real saddle pair and two imaginary pairs):

| Point | lambda_i in the CR3BP | lambda_i range in the BCR4BP |
| --- | --- | --- |
| L1 / E1(theta) | +/-2.932; +/-2.334 i; +/-2.269 i | +/-2.916 to +/-2.940; +/-2.324 i to +/-2.337 i; +/-2.258 i to +/-2.276 i |
| L2 / E2(theta) | +/-2.159; +/-1.863 i; +/-1.786 i | +/-2.137 to +/-2.200; +/-1.847 i to +/-1.887 i; +/-1.776 i to +/-1.811 i |

Table 2 (p2198), "Eigenvalues associated with the equilibrium point in the CR3BP and in the BCR4BP", Sun-B1 rotating frame,
nondimensional (the point labels are underlined in the paper):

| Point | lambda_i in the CR3BP | lambda_i range in the BCR4BP |
| --- | --- | --- |
| L1 / E1(theta) | +/-2.533; +/-2.087 i; +/-2.015 i | +/-2.531 to +/-2.538; +/-2.084 i to +/-2.090 i; +/-2.015 i to +/-2.018 i |
| L2 / E2(theta) | +/-2.484; +/-2.057 i; +/-1.985 i | +/-2.483 to +/-2.489; +/-2.055 i to +/-2.060 i; +/-1.985 i to +/-1.988 i |

Table 3 (p2209), "Lyapunov exponents for the sampled synodic resonant NRHOs as computed in the CR3BP and the BCR4BP", columns
headed "CR3BP, 1 period" and "BCR4BP, 1 period". Units: the exponent of eq. (14) (per nondimensional time, INFERRED). Each
orbit lists three of the six exponents (the others are their negatives, by the reciprocal-pair property, p2198):

| Orbit | CR3BP, 1 period | BCR4BP, 1 period |
| --- | --- | --- |
| 3:1 NRHO | 0; 0; 0 | 0; 0; +/-0.6223 |
| 4:1 NRHO | 0; 0; +/-0.6277 | 0; 0; +/-0.6364 |
| 9:2 NRHO | 0; 0; +/-0.5157 | 0; 0; +/-0.5208 |
| 5:1 NRHO | 0; 0; 0 | 0; 0; +/-0.0322 |

(The 5:1 row is the NRHO+ of the text; the table caption's label is "5:1 NRHO".)
Other printed numbers (prose): distances from L1 to the instantaneous equilibrium E1 range 200 to 650 km and L2 to E2 400 to 1700
km in the Earth-Moon frame (p2196, Fig. 2); in the Sun-B1 frame 600 to 2500 km for both (p2197). Jacobi constants and Hamiltonian
values are shown only in Figs. 20 and 21 (read off Fig. 21: the four resonant orbits lie near C between 3.02 and 3.06; not
tabulated).

## 6. Positive controls for the project

Honest answer: this paper cannot supply a digit-level closure test, because it prints no state, period or model constants (section
0). What can be tested, in order of usefulness (all INFERRED recipes; the printed targets are READ):

1. CR3BP half (needs `core/cr3bp.py` and an L2 halo family continuation): find the L2 halo orbits with period equal to one third,
   one quarter, two ninths and one fifth of the lunar synodic period `2 pi / omega` = 6.7912 TU (project constants; 1.6978 TU
   for 4:1). Check perilune and apolune (3:1: about 15,000 and 84,500 km; 4:1: about 5,600 and 75,335 km; 9:2: about 3,100 and 71,000
   km; 5:1 NRHO+: 1,650 and 66,000 km) and the CR3BP Lyapunov exponents of Table 3 (4:1: 0.6277; 9:2: 0.5157; 3:1 and 5:1:
   all zero). Exponent is `Re(ln(lambda)) / T` with T the orbit's own period; the tolerance is the printed four digits. The
   period of the 9:2 orbit is 2 T_syn / 9. These involve no Sun and so are independent of the sense of the Sun's motion.
2. BCR4BP half, sense-sensitive: stack P revolutions, continue in epsilon with the Sun at `theta = theta_0 - omega t` (the project's
   `theta_S0 - omega_S t`, no sign change: the paper's Sun angle decreases with time), find the T_syn-periodic orbit at epsilon = 1
   and compare:
   - 3:1: a single counterpart; perilune Sun angles at 0, 120, 240 degrees (all one orbit); perilunes 14,300 to 15,200 km, apolunes
     82,300 to 89,400 km; Lyapunov exponent pair +/-0.6223, with the stability-change at about 15 to 17 percent of the Sun mass; a
     second branch with perilune Sun angles 60, 180, 300 degrees that exists only for epsilon between 0.01 and 0.17.
   - 4:1: two counterparts, A with perilune Sun angles 45, 135, 225, 315 degrees (perilunes about 5,940 km) and B with Sun angles
     0, 90, 180, 270 degrees (perilunes about 4,970 km and 7,100 km); apolunes 74,800 to 75,500 km; Lyapunov exponent +/-0.6364 for
     one or both (the paper does not say which).
   - 9:2: Sun angles at perilune repeat every 40 degrees; stability exponent +/-0.5208.
   Sun-angle convention: measured from the Earth-Moon rotating +x axis (toward the Moon) to the Sun, as in the project (READ,
   p2196; the +x direction toward the Moon is INFERRED from Fig. 1).
3. Sense check: a 4:1 BCR4BP orbit built with the reversed Sun gives different perilune Sun-angle patterns, not published; the
   printed mirror of the counterpart A/B offset of 45 degrees is a useful consistency check only if the Sun's sense is right
   (INFERRED, not a printed test).
Best candidate (hand-back): the 4:1 NRHO, because both CR3BP (0.6277) and BCR4BP (0.6364) exponents and two counterparts with
specific Sun angles are printed. Printed numbers: 4:1, perilune about 5,600 km, apolune about 75,335 km, one quarter of the
synodic period; BCR4BP A perilune Sun angles 45, 135, 225, 315 degrees (about 5,940 km), B 0, 90, 180, 270 degrees (about 4,970
km and 7,100 km); exponents 0.6277 (CR3BP) and 0.6364 (BCR4BP). Expected agreement with project constants (unknown model
constants in the paper): the exponents to about 1e-3 if the project reproduces the paper's model (INFERRED).

## 7. What this means for #884

All READ, with the limits noted.
- Method: published. "In this investigation, each BCR4BP periodic orbit is constructed using this method, in which a synodic
  resonant CR3BP orbit serves as the initial guess for a continuation process with the Sun's mass as the continuation parameter."
  (p2201.) This is the project's method (commensurate-period three-body orbit, epsilon from 0 to 1).
- Several equivalents per three-body orbit: published, found by scanning the Sun angle at perilune; counts 1 (3:1), 2 (4:1), many
  (9:2). The reason is not given as a symmetry argument (in contrast to Oshima 2022, who uses the time-reversal and mirror
  symmetries).
- Branches that do not reach full Sun mass: published as observation (3:1 offset branch up to 17 percent; 9:2 branches that "do not
  reach 100%"), with the open question whether the continuation or the orbit is at fault. The return to epsilon = 0 is not shown
  here; it is shown by Oshima (2022) for 12:11 retrograde orbits.
- Stability inherited from the parent: published as a mixed result: 4:1 and 9:2 (unstable parents) keep a small unstable mode;
  3:1 and 5:1 (stable parents) become unstable at about 15 and about 90 percent of the Sun mass. So "only forced orbits whose CR3BP
  parent is stable are stable" is not what this paper says; here stable parents lose stability. Check against the #884 literature
  note's wording.
- Not in the paper: Earth-Moon cyclers, the sense-of-motion question, the model constants, any digit-level target.
- For the novelty question: the mechanism is in print for NRHOs and 1:1 Lyapunov orbits (Boudad et al. 2020) and retrograde orbits
  (Oshima 2022); the project's application to cycler families is not covered by either paper.

## 8. References cited (relevant to periodic orbits of the bicircular or quasi-bicircular problem)

Bicircular, quasi-bicircular and four-body periodic-orbit work cited by Boudad et al. (p2195, p2201, p2214):
- Andreu, M.A., 1998. The Quasi-Bicircular Problem. Ph.D. Dissertation. Universitat de Barcelona, Barcelona, Spain.
- Gomez, G., Llibre, J., Martinez, R., Simo, C., 2001. Dynamics And Mission Design Near Libration Points - Vol II: Fundamentals: The Case Of Triangular Libration Points. World Scientific Monograph Series In Mathematics, World Scientific Publishing Company. doi:10.1142/4402. (As printed, the author list reads Gomez, Llibre, Martinez, Simo; Oshima cites the same book with Gomez, Jorba, Masdemont, Simo.)
- Jorba-Cusco, M., Farres, A., Jorba, A., 2018. Two periodic models for the earth-moon system. Front. Appl. Mathe. Stat. 4, 32. (Held and digested.)
- Boudad, K.K., 2018. Disposal Dynamics From The Vicinity Of Near Rectilinear Halo Orbits In The Earth-Moon-Sun System. M.S. Thesis. Purdue University, West Lafayette, Indiana.
- Boudad, K.K., Davis, D.C., Howell, K.C., 2019. Near rectilinear halo orbits in cislunar space within the context of the bicircular four-body problem. In: 2nd IAA/AAS SciTech Forum, Moscow, Russia.
- Bosanac, N., 2012. Exploring the Influence of a Three-Body Interaction Added to the Gravitational Potential Function in the Circular Restricted Three-Body Problem: a Numerical Frequency Analysis M.S. Thesis. Purdue University. (Cited for the stability index.)
- Broucke, R.A., 1969. Stability of periodic orbits in the elliptic, restricted three-body problem. AIAA J. doi:10.2514/3.5267.
- Keller, H.B., 1976. Solution of two point boundary value problems. In: AAS/AIAA Astrodynamics Specialist Conference. Society for Industrial and Applied Mathematics. (Parallel shooting.)

Other references in the paper (halo and NRHO background): Davis (2011); Davis, Phillips, Howell, Vutukuri, McCarthy (2017);
Grebow (2006); Hambleton (2017); Henon (1997), "Generating Families in the Restricted Three-Body Problem", Springer; Howell (1998),
J. Astronaut. Sci. 49; Huang (1960), NASA technical note; Lee (2019), Gateway Destination Orbit Model; Meyer and Hall (2013);
Ortiz Longo and Rickman (1995); Roy (2004), Orbital Motion, fourth ed.; Szebehely (1967); Warner (2018); Williams, Lee, Whitley,
Bokelmann, Davis, Berry (2017); Yakubovich and Starzhinskii (1975); Zimovan, Howell, Davis (2017) and Zimovan-Spreen, Howell, Davis
(2020, Celest. Mech. Dyn. Astron.); Zimovan-Spreen and Howell (2019). Full page ranges not transcribed for these (not needed for
the bicircular question).
