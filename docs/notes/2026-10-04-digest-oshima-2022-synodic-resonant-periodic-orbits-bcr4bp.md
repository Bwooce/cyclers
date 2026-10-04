# Digest: Oshima (2022), "Multiple families of synodic resonant periodic orbits in the bicircular restricted four-body problem"

K. Oshima, Advances in Space Research 70:1325-1335 (2022), DOI 10.1016/j.asr.2022.06.009, 11 pages, Hiroshima Institute of
Technology. Received 17 March 2022, accepted 5 June 2022.
Filed in the private paper corpus as
oshima-2022-multiple-families-synodic-resonant-periodic-orbits-bicircular-restricted-four-body-problem-asr-70-1325-doi-10.1016-j.asr.2022.06.009.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with page, section, equation, figure or table number) or
INFERRED (our reading or arithmetic). Page numbers are the journal's printed numbers (1325 to 1335). All numbers were read
from page images; where a figure is the only source the value is marked "read off the figure" and is approximate.
Context: tasks #884 and #891 (the project's bicircular model had the Sun going round the wrong way; corrected in
`core/bcr4bp.py`). This paper is the published work on the mechanism the project met: continuing a three-body periodic orbit
into the bicircular problem at a Sun-commensurate period, several equivalents per orbit, branches that fold back.

## 0. What the paper is

READ (abstract, p1325). "The present paper deals with two mechanisms of the generation of multiple families of synodic
resonant periodic orbits in the bicircular restricted four-body problem through numerical examples adopting planar and
three-dimensional retrograde periodic orbits around the Earth." Part one (sections 3 and 4) is planar: period-multiplying
bifurcations of the Earth-Moon CR3BP orbit interplay with a 12:11 synodic resonant orbit and give a second 12:11 orbit in the
BCR4BP. Part two (section 5) is spatial: a doubly symmetric 1:1 orbit gives four families.

## 1. The model as printed

### 1.1 Equations of motion (section 2.1, p1326, eqs. 1 to 5) - READ

Non-dimensional, Earth-Moon rotating frame (p1326: "Non-dimensional equations of motion for the CR3BP and BCR4BP in the
Earth-Moon rotating frame"):

    xddot - 2 ydot = -dU/dx,   yddot + 2 xdot = -dU/dy,   zddot = -dU/dz          (1)

with a continuous parameter epsilon unifying the two models:

    U = U_3BP - eps m_S / r_3 + eps m_S / a_S^2 (x cos(theta_S) + y sin(theta_S))   (2)

    U_3BP = -1/2 (x^2 + y^2) - (1 - mu)/r_1 - mu/r_2 - 1/2 mu (1 - mu)             (3)

    r_1 = sqrt((x + mu)^2 + y^2 + z^2)
    r_2 = sqrt((x - 1 + mu)^2 + y^2 + z^2)
    r_3 = sqrt((x - a_S cos(theta_S))^2 + (y - a_S sin(theta_S))^2 + z^2)          (4)

    theta_S = theta_S0 + omega_S t                                                  (5)

"and theta_S0 is the solar phase angle at initial time t = 0." Overbars on U in the printed equations denote this potential
(the sign convention is that of a negative pseudo-potential, so that the accelerations are minus its gradient).
The rendering of the third term of eq. (2) was read at page-image resolution; its sign and the cos/sin pairing are as given
above (READ at limited resolution).

- Frame and primaries (READ, eqs. 3 and 4): Earth at (-mu, 0, 0), Moon at (1 - mu, 0, 0). This is the same placement as
  `core/bcr4bp.py`.
- Epsilon (READ, p1326): "epsilon = 1 expresses the BCR4BP dynamics incorporating the full mass of the Sun m_S and epsilon = 0
  corresponds to the CR3BP dynamics ignoring the solar gravity." Epsilon scales the Sun mass in both the direct and the indirect
  term.
- Indirect term (READ, eq. 2): `+ eps m_S / a_S^2 (x cos(theta_S) + y sin(theta_S))`. INFERRED: its gradient gives the
  acceleration `-eps m_S (cos(theta_S), sin(theta_S)) / a_S^2 = -eps m_S r_Sun / a_S^3` for a Sun at
  `a_S (cos(theta_S), sin(theta_S))`, which is the cancelling term of the direct attraction; direct and indirect terms are
  mutually consistent in the printed equations (no sign slip seen), and identical in form to the project's
  `_sun_acceleration`.
- The Sun is not on the Earth-Moon axis by construction; the Sun moves in the Earth-Moon plane (p1326, section 2: "the Earth-Moon
  barycenter and the Sun move on circular orbits around their common barycenter on the same orbital plane of the Earth and
  Moon").

### 1.2 The sense of the Sun's motion - READ

Eq. (4) puts the Sun at `(a_S cos(theta_S), a_S sin(theta_S), 0)`. Eq. (5) is `theta_S = theta_S0 + omega_S t`, and Table 1
prints `omega_S = -0.925195985 1/TU`, a NEGATIVE number. So `theta_S` DECREASES with time: the Sun moves clockwise as seen from
+z in the Earth-Moon rotating frame (retrograde against the primaries' own counter-clockwise inertial motion of the frame).
Sentence for the hand-back: Oshima's Sun is at angle theta_S = theta_S0 + omega_S t with omega_S = -0.925195985, so the angle
decreases and the Sun goes clockwise, the same sense as the corrected `core/bcr4bp.py` (theta_S = theta_S0 - omega_sun t with
omega_sun = +0.925195985). No sign flip and no angle shift is needed to transfer a state.
The period `|2 pi / omega_S|` is quoted as "approximately 29.5 days" (p1327, section 3).

### 1.3 Constants - Table 1 (p1326), READ, "Parameters used for the CR3BP and BCR4BP (Topputo, 2013)"

| Parameter | Value | Unit |
| --- | --- | --- |
| Distance unit (DU) | 384405 | km |
| Time unit (TU) | 4.34811305 | day |
| Velocity unit (VU) | 1.02323281 | km/s |
| Mass parameter mu | 0.0121506683 | - |
| Sun's mass m_S | 328900.541 | - |
| Sun's orbital radius a_S | 388.811143 | DU |
| Sun's angular velocity omega_S | -0.925195985 | 1/TU |
| Earth's radius | 6378 | km |
| Moon's radius | 1738 | km |

INFERRED comparison with the project constants (`core/bcr4bp.py`, Andreu set):

| Constant | Oshima (Topputo) | `andreu_default()` | Difference |
| --- | --- | --- | --- |
| mu | 0.0121506683 | 0.0121505816 | 8.7e-8 absolute (7e-6 relative) |
| m_S | 328900.541 | 328900.5423094043 | 1.3e-3 absolute (4e-9 relative) |
| a_S | 388.811143 | 388.8111430233511 | 2e-8 |
| omega (magnitude) | 0.925195985 | 0.925195985520347 | 5.2e-10 |

The mass parameter differs in the 7th significant digit. That matters for any closure test against Tables 2 to 5 (section 6
below): build a `BCR4BPSystem` with Oshima's printed constants, not `andreu_default()`.

## 2. Synodic resonance, family selection and the symmetry argument

### 2.1 Definition (section 3, p1327) - READ

"Periodic orbits in the BCR4BP must be in resonance with the orbital motion of the Sun in the Earth-Moon rotating frame. In
other words, the period of a synodic resonant periodic orbit is

    T_4BP = N x |2 pi / omega_S|                                    (10)

with an integer N and |2 pi / omega_S| approx 29.5 days. An M:N synodic resonant periodic orbit exhibits an M-revolutional
geometry while the Sun revolves N times in the Earth-Moon rotating frame. Thus, the period of a single-revolutional orbit in the
CR3BP that can be expanded into the M-revolutional M:N synodic resonant periodic orbit is

    T_3BP = (N/M) x |2 pi / omega_S|                                (11)"

So M is the number of revolutions of the orbit geometry and N the number of Sun synodic revolutions (months) in one period; the
CR3BP parent has period (N/M) synodic months. INFERRED: same direction as Boudad et al.'s P:Q (orbital periods : synodic
periods), with M = P and N = Q.

### 2.2 Symmetries (section 2.2, p1326-1327, eqs. 6 to 9) - READ

The equations of motion are invariant under

    s1: (x,y,z,vx,vy,vz,t,theta_S) -> (x,-y,-z,-vx,vy,vz,-t,-theta_S)       (6)
    s2: (x,y,z,vx,vy,vz,t,theta_S) -> (x,-y,z,-vx,vy,-vz,-t,-theta_S)       (7)
    s3: (x,y,z,vx,vy,vz,t,theta_S) -> (x,y,-z,vx,vy,-vz,t,theta_S)         (8)
    s_p (planar): (x,y,vx,vy,t,theta_S) -> (x,-y,-vx,vy,-t,-theta_S)       (9)

"These properties reflect symmetries with respect to the x-axis, xz-plane, and xy-plane, respectively." s1 and s2 are
time-reversal symmetries (the Sun angle also changes sign), s3 is not. "Ignoring theta_S reduces to the symmetries in the
CR3BP (Russell, 2006)." Symmetric periodic orbits must start at `y0 = vx0 = z0 = 0` or `y0 = vx0 = vz0 = 0` "with theta_S0 = 0
or theta_S0 = pi" (p1327); in the planar problem `y0 = vx0 = 0` with theta_S0 = 0 or pi.

### 2.3 How the families are labelled and why there are several (p1327-1328) - READ

- Odd M (Fig. 2a): "The difference in theta_S0 distinguishes the S_0 and S_pi families" (initial point on the x axis, Sun at
  phase 0 or pi).
- Even M (Fig. 2b): "the difference in states separated by a half period distinguishes the T_0 and T_1/2 families." The paper
  adopts T_0 and T_1/2 instead of the "left" and "right" families of Campagnola et al. (2008), because with the same x0 at a
  half period apart "the definition of the left and right families is not applicable".
- Why: "Although S_0 and S_pi families or T_0 and T_1/2 families are identical in the Earth-Moon CR3BP, their dynamical
  differences arise once the solar gravity appears." (p1328). INFERRED: the CR3BP orbit has several symmetric starting points
  (equivalent in the autonomous problem) that become inequivalent when the Sun phase is attached; each gives a distinct
  continuation. This is the same mechanism as the project's "equivalents at symmetric Sun phases".
- Spatial doubly symmetric orbit (section 5, p1330): "Since the number of revolution of the 1:1 synodic resonant orbit is odd,
  Fig. 2(a) indicates that S_0 and S_pi families ... would emerge in the BCR4BP once I set initial time t = 0 and theta_S0 = 0 or
  theta_S0 = pi on a point satisfying y = vx = 0 and start the continuation procedure. However, the doubly symmetric orbit in
  Fig. 7 provides two candidates for the special points y = vx = z = 0 and y = vx = vz = 0. Therefore, S_0 and S_pi families
  can emanate from each of the cases ... and in total four distinct families of synodic resonant periodic orbits may emerge in
  the BCR4BP from a doubly symmetric orbit in the CR3BP." (p1330-1331). The four are named z0-S0, z0-Spi (start `y0 = vx0 = z0 =
  0`) and vz0-S0, vz0-Spi (start `y0 = vx0 = vz0 = 0`) (p1331).
- Period-multiplying route (section 4, p1329-1330): an m-revolutional period-m CR3BP orbit also gives an M-revolutional M:N orbit
  when M/m is an integer: "if m is a factor of M, period-m orbits can be hopeful candidates for generating M:N synodic resonant
  periodic orbits as well as misleading mazes connected with an original single-revolutional orbit."
- Sun phase: the Sun angle at the start is fixed to 0 or pi by the symmetry (not scanned). Other Sun phases are not explored in
  this paper (INFERRED: no phase scan; contrast Boudad et al., who scan the Sun angle at perilune).

## 3. Orbit families treated

READ and INFERRED, plainly: NO Earth-Moon cycler, no transfer orbit, and no resonant orbit with lunar flybys is treated.

- Planar part: retrograde periodic orbits (RPOs) around the Earth, 12:11 synodic resonant, orbits that were taken from
  Oshima (2022a, ASR 69:2210). Initial conditions of these are NOT printed in this paper (only the plot of x0 against epsilon,
  Fig. 3, and the surface-of-section plots, Figs. 5 and 6). The one stated value of the starting point in the figures: read off
  Fig. 3 the T_0 branch starts near x0 = -1.45 and the T_1/2 branch near x0 = -0.55 at epsilon = 0 (approximate).
- Spatial part: a doubly symmetric, linearly stable, 1:1 synodic resonant three-dimensional RPO (Oshima 2022b, "3D stable and
  weakly unstable periodic orbits around the Earth near the retrograde co-orbital resonance with the Moon"). INFERRED from the
  tables: positions at y = 0 are at about 1.09 to 1.12 DU from the Earth with |z| about 0.18 to 0.20 DU, so the orbit stays
  about 0.2 DU (about 80000 km) or more from the Moon at those crossings; no lunar flyby is discussed anywhere in the paper.

## 4. What happens along the continuation

### 4.1 12:11 planar family (section 4, p1328-1330, Figs. 3 to 6) - READ

- Method: pseudo-arclength continuation in epsilon (Keller 1977; details in Oshima 2022a). "However, the algorithm does not
  guarantee the convergence into periodic orbits in the BCR4BP (epsilon = 1)." (p1328)
- Fig. 3 caption: "Continuation pathways of T_0 and T_1/2 families of 12:11 synodic resonant planar RPOs around the Earth
  computed in Oshima (2022a)."
- Decisive statement (p1328): "The T_1/2 family successfully reaches epsilon = 1, but the T_0 family encounters a fold point
  near epsilon = 1 and eventually returns back to epsilon = 0. Such a spontaneous evolution of epsilon is owing to the
  pseudo-arclength continuation procedure, but it nevertheless fails to reach epsilon = 1."
- Extension beyond the usual stopping rule (p1328): "it is not mandatory to stop the continuation at epsilon = 0 or epsilon = 1
  in Fig. 3. Extending the continuation curve of the T_0 family in Fig. 3 would lead to negative epsilon corresponding to a
  negative mass of the Sun, which is not physically meaningful at first glance. ... However, there could be the possibility that
  extended continuation curves make it back through the unphysical realms to 0 <= epsilon <= 1 and reach epsilon = 1 generating
  other solutions converged into the BCR4BP." The extension is run "until it returns back to the initial point, from which the
  continuation starts, and thus a further extension would trace the identical closed curve."
- Fig. 4 result (p1328-1329): both extended curves go to substantially negative epsilon (read off the figure: T_0 down to about
  -8, T_1/2 down to about -3.5) "but they eventually return back to the initial point (i) and produce the closed curves." Labels:
  lower-case roman numerals i to v are solutions at epsilon = 0, upper-case I and II are solutions at epsilon = 1. "(ii) and (I)
  are identical to the end points of the continuation pathways of the T_0 and T_1/2 families in Fig. 3, respectively, whereas (II)
  is a newly found periodic orbit in the BCR4BP. The extended pathway of the T_1/2 family beyond (I) instantly reaches (II)
  crossing epsilon = 1 again." (Fig. 4 text; "(ii)" in the first clause is as printed.)
- "It is noticeable that the extended curve of the T_0 family crosses epsilon = 0 multiple times indicated by (i)-(iv), but is
  not able to reach epsilon = 1. The curve of the T_1/2 family also crosses epsilon = 0 indicated by (v). Indeed, the newly
  found solution (II) in the BCR4BP is straightforwardly accessible from (v), which is a periodic orbit in the CR3BP." (p1329)
- The CR3BP solutions (i) to (v) (Fig. 6, p1330): "Although (i) is the original single-revolutional orbit, the others are
  multi-revolutional higher-period orbits that are originated from the single-revolutional family via period-multiplying
  bifurcations (Lara et al., 2007; Zimovan et al., 2020). The number of intersections on the surface of section reveals that
  (ii) is a period-2 orbit, (iii) is a period-12 orbit, (iv) is a period-2 orbit, and (v) is a period-12 orbit. Moreover, (ii)
  and (iv) are identical whereas (iii) and (v) are distinct, which are understandable because (ii) and (iv) are on the same
  continuation curve in Fig. 4 whereas (iii) and (v) are not. See Pushparaj et al. (2021) for the existence of two distinct
  families for a multi-revolutional orbit."
- Stability colour: max over the four Lyapunov exponents `phi_k = Re(ln(lambda_k) / T_4BP)` of the monodromy matrix, with
  `T_4BP = 11 x |2 pi / omega_S|` (eq. 12, p1328). Zero is linear stability. The extended curves are coloured by this in Fig. 4;
  the colour scale runs 0 to 0.6. No numeric stability values for individual 12:11 members are printed in the text (only the
  colour map).
- Conclusion (section 6, p1334): "M-revolutional synodic resonant periodic orbits in the BCR4BP may emerge from m-revolutional
  period-m orbits in the CR3BP if m is a factor of M. The newly found 12:11 synodic resonant planar RPO in the BCR4BP naturally
  connects with one of the period-12 orbits in the CR3BP."

### 4.2 1:1 spatial families (section 5, p1330-1333, Figs. 7 to 12) - READ

- All four families converge straightforwardly to epsilon = 1 "without the interference from higher-period orbits as expected"
  (p1331).
- Fig. 8 (p1331): the absolute values of the six monodromy eigenvalues along epsilon. For (a) z0-S0 and z0-Spi, a pair of
  real eigenvalues leaves the unit circle at epsilon > 0, reaching about 1.06 and 0.94 at epsilon = 1 (read off the figure,
  approximate); for (b) vz0-S0 and vz0-Spi, all six stay at modulus 1 over the whole path.
- Statement (p1331): "A periodic orbit is linearly stable if and only if each of |lambda| equals unity and thus the z0-S0 and
  z0-Spi families are weakly unstable whereas the vz0-S0 and vz0-Spi families are linearly stable in the BCR4BP. Note that the
  instabilities of the z0-S0 and z0-Spi families are peculiar to the presence of the solar gravity as they only appear at
  epsilon > 0 and the original orbit in the CR3BP (epsilon = 0) is linearly stable."
- Differences between the S0 and Spi members (p1331, Fig. 10): "The evolutions of z are in opposite phase between the S_0 and
  S_pi families in terms of theta_S. Moreover, the out-of-plane amplitude |z| of the vz0-S0 and vz0-Spi families is larger than
  that of the z0-S0 and z0-Spi families." "Note that the converged orbits in the BCR4BP slightly violate the doubly symmetric
  property of the original orbit in the CR3BP."
- No folds and no bifurcations are reported for the 1:1 families.
- Long-term test (p1332-1333, Table 6 and Figs. 11, 12): displacements of eight sizes added to the apolune state of #2 in Table
  4 and propagated three years in the BCR4BP. READ: "The cases of (1) and (2) stay in the vicinity of the RPO. The larger
  displacements in the cases of (3) and (4) still generate regular behaviors. The stability boundary appears to exist between
  the conditions of (5)-(6) and (7)-(8) indicating that the corresponding stability region is substantially large."

Table 6 (p1333), displacements in the Earth-Moon rotating frame (READ; units as printed):

| Case | (1) | (2) | (3) | (4) | (5) | (6) | (7) | (8) | Unit |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dx | +1000 | -1000 | +5000 | -5000 | +10000 | -10000 | +20000 | -20000 | km |
| dy | +1000 | -1000 | +5000 | -5000 | +10000 | -10000 | +20000 | -20000 | km |
| dz | +1000 | -1000 | +5000 | -5000 | +10000 | -10000 | +20000 | -20000 | km |
| dvx | +5 | -5 | +25 | -25 | +50 | -50 | +100 | -100 | m/s |
| dvy | +5 | -5 | +25 | -25 | +50 | -50 | +100 | -100 | m/s |
| dvz | +5 | -5 | +25 | -25 | +50 | -50 | +100 | -100 | m/s |

## 5. The tables (all transcribed digit by digit)

Tables 2 to 5 are on p1332. Captions are quoted: "Position, velocity, and theta_S of the [family] at y = 0. The modulo
operation expresses theta_S between 0 and 2 pi." The captions give no units. INFERRED: non-dimensional (DU, DU/TU, rad) in the
Earth-Moon rotating frame of Table 1, evaluated on the BCR4BP periodic orbit (epsilon = 1). `y = 0` at all four rows, so the
state is `(x, 0, z, vx, vy, vz)`; the y column is omitted in the table. All digits legible; none marked unclear.

Table 2, z0-S0 family:

| # | x | z | vx | vy | vz | theta_S |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | -1.107328855 | 0 | 0 | 2.056256309 | -0.156230339 | 0 |
| 2 | 1.111203054 | -0.180245301 | -0.000134973 | -2.062045643 | -0.000001655 | 4.712253115 |
| 3 | -1.106981252 | 0 | 0 | 2.056236796 | 0.156277994 | 3.141592654 |
| 4 | 1.111203054 | 0.180245301 | 0.000134973 | -2.062045643 | -0.000001655 | 1.570932192 |

Table 3, z0-Spi family:

| # | x | z | vx | vy | vz | theta_S |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | -1.106981252 | 0 | 0 | 2.056236796 | -0.156277994 | 3.141592654 |
| 2 | 1.111203054 | -0.180245301 | 0.000134973 | -2.062045643 | 0.000001655 | 1.570932192 |
| 3 | -1.107328855 | 0 | 0 | 2.056256309 | 0.156230339 | 0 |
| 4 | 1.111203054 | 0.180245301 | -0.000134973 | -2.062045643 | 0.000001655 | 4.712253115 |

Table 4, vz0-S0 family:

| # | x | z | vx | vy | vz | theta_S |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.090174251 | -0.204803847 | 0 | -2.061909684 | 0 | 0 |
| 2 | -1.120233045 | 0.000079580 | -0.000178477 | 2.042532822 | 0.177628656 | 4.712573766 |
| 3 | 1.090649738 | 0.204909100 | 0 | -2.061914819 | 0 | 3.141592654 |
| 4 | -1.120233045 | 0.000079580 | 0.000178477 | 2.042532822 | -0.177628656 | 1.570611541 |

Table 5, vz0-Spi family:

| # | x | z | vx | vy | vz | theta_S |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.090649738 | -0.204909100 | 0 | -2.061914819 | 0 | 3.141592654 |
| 2 | -1.120233045 | -0.000079580 | 0.000178477 | 2.042532822 | 0.177628656 | 1.570611541 |
| 3 | 1.090174251 | 0.204803847 | 0 | -2.061909684 | 0 | 0 |
| 4 | -1.120233045 | -0.000079580 | -0.000178477 | 2.042532822 | -0.177628656 | 4.712573766 |

Text accompanying the tables (p1332): "The z0-S0 and z0-Spi families satisfy only y = vx = z = 0 whereas the vz0-S0 and vz0-Spi
families intersect only with y = vx = vz = 0." "Note also that S_0 and S_pi reflect the symmetric property with respect to the
xy-plane in Eq. (8), e.g., flipping the sign of the out-of-plane components z and vz of #2 of the z0-S0 family coincides with
#4 of the z0-Spi family."

INFERRED consistency checks on the printed numbers (arithmetic only, not entries):
- Tables 2 and 3: theta_S of rows #2 and #4 sum to 2 pi (4.712253115 + 1.570932192 = 6.283185307); Table 4 and 5 likewise
  (4.712573766 + 1.570611541 = 6.283185307).
- Each Table 3 row is the Table 2 row at the same theta_S with z and vz changed in sign (the s3 mirror), for example Table 3 #1
  against Table 2 #3, and Table 3 #2 against Table 2 #4. Table 5 against Table 4 the same way (Table 5 #1 against Table 4 #3).
  So the four families come as two mirror pairs; the independent information is in Tables 2 and 4.
- The theta_S values at the second and fourth crossings are not exactly 3 pi / 2 and pi / 2 (4.712253115 against 4.712388980),
  consistent with the printed remark that the BCR4BP orbits "slightly violate" the double symmetry.
- The z0 rows have vx = 0 and z = 0 at the symmetric crossings (#1, #3); the vz0 rows have vx = 0 and vz = 0 (#1, #3), as the
  symmetry argument requires.

## 6. Positive controls for the project

All use `core/bcr4bp.py` (`BCR4BPSystem`, `propagate_bcr4bp(system, state6, t, t0=...)`), state ordering (x, y, z, vx, vy, vz)
(INFERRED from the signature; check the docstring before use).

Common recipe (INFERRED):
1. Build `BCR4BPSystem(mu=0.0121506683, mu_sun=328900.541, a_sun_nondim=388.811143, omega_sun_nondim=0.925195985,
   theta_sun0=theta_k)`, with the printed Table 1 constants (not `andreu_default()`: mu differs by 8.7e-8, which alone moves
   the orbit by about 1e-6 to 1e-5 over a period; this size is an estimate).
2. State `(x, 0, z, vx, vy, vz)` from the table row k, `theta_sun0 = theta_S` of that row, start time t0 = 0.
3. No sign flip and no angle shift: the frame (Earth at -mu, Moon at 1 - mu), the velocities (rotating-frame xdot, ydot, zdot),
   and the Sun's sense (omega_S = -0.925195985 with theta + omega t, against the project's theta - omega t with +omega) all
   agree, so Oshima's `theta_S` is the project's `theta_sun0` as printed. READ for Oshima's side, INFERRED for the identity.
4. Period: T = 2 pi / 0.925195985 = 6.791193875727 TU (29.5289 d with TU = 4.34811305 d) - one synodic month, since this is a
   1:1 orbit (N = 1). The project's own `omega` would give T 3.8e-9 TU longer; use the printed omega for the test.

Candidate A (best): Table 4 row #1, the vz0-S0 orbit. It is linearly stable (all |lambda| = 1, Fig. 8b), so the printed 9-digit
rounding cannot be amplified by the dynamics.
- Initial state: `(x, y, z, vx, vy, vz) = (1.090174251, 0, -0.204803847, 0, -2.061909684, 0)`, `theta_sun0 = 0`.
- Propagate T = 6.791193875727. Expected closure (INFERRED, not printed): position and velocity miss the start by of order
  1e-8 or less (printed digits are rounded to 1e-9), plus the integrator error (set rtol = atol = 1e-12 or tighter).
- Intermediate crossings (INFERRED, derived from the printed theta_S assuming theta falls at |omega| and one time axis):
  the time from row #1 to row k is `t_k = (2 pi - theta_k) / 0.925195985`, giving #2 at t = 1.697598743 TU, #3 at 3.395596937
  TU, #4 at 5.093595133 TU. At each, y should be 0 (crossing) and the state should reproduce the printed row to about 1e-8 (the
  printed rows are on the same orbit). Crossing times are not printed in the paper; they follow from the printed Sun angles.
- Negative controls: the same start with the Sun's sense reversed (the project's pre-#891 sense) should fail to close by a gross
  margin, and with the Sun off (mu_sun = 0) the orbit should not close either (the periodic orbit is a Sun-forced one at
  epsilon = 1). The size of these misses is not printed.

Candidate B: Table 2 row #1, the z0-S0 orbit: `(-1.107328855, 0, 0, 0, 2.056256309, -0.156230339)`, `theta_sun0 = 0`, T as
above. Weakly unstable (a real pair of about 1.06 and 0.94 per period, read off Fig. 8a), so rounding of 1e-9 grows to about
1e-9 times 1.06, still of order 1e-8 after one period. Crossing times from row #1: #2 at 1.697945319 TU, #3 at 3.395596937
TU, #4 at 5.093248557 TU (INFERRED).

Candidate C (second epoch test): start from Table 4 row #3 with `theta_sun0 = 3.141592654` (Sun on the opposite side), state
`(1.090649738, 0, 0.204909100, 0, -2.061914819, 0)`; closure after T tests that the project reproduces the same orbit when the
Sun phase is shifted by pi.

Not available: the 12:11 planar orbits (Figs. 3 to 6) have no tabulated initial conditions in this paper, so nothing there can be
tested digit by digit. The 12:11 fold (T_0 branch fold near epsilon = 1) can only be tested qualitatively after obtaining the
initial conditions from Oshima (2022a), which is not held.

## 7. What this means for #884

READ with quotes; all are statements about retrograde orbits around the Earth, not cyclers.

- Several equivalents per three-body orbit from symmetry: published. Planar, two families (S_0/S_pi for odd M, T_0/T_1/2 for even
  M); spatial doubly symmetric, four families (section 2.3 above). The argument that these are identical in the CR3BP and differ
  once the Sun is present is quoted at p1328.
- Folds returning to zero Sun mass: published, for the 12:11 T_0 family: "the T_0 family encounters a fold point near epsilon
  = 1 and eventually returns back to epsilon = 0", and the extended curve closes on itself after crossing epsilon = 0 several
  times. Extending past epsilon = 0 and 1 (negative and super-unit epsilon) is also published, with a new epsilon = 1 solution
  (II) found by crossing epsilon = 1 a second time.
- Higher-period parents as sources of additional epsilon = 1 solutions: published (period-2 and period-12 orbits, the "mazes").
  The project's commensurate members may be seen the same way: a period-m parent with m dividing M gives M:N solutions.
- Stability inherited from the parent: NOT a general statement here. READ: the stable 1:1 parent gives two linearly stable
  families (vz0) and two weakly unstable ones (z0) in the BCR4BP, the instability "only appear[ing] at epsilon > 0". So
  stability is not always inherited, and Oshima gives no rule for which of the symmetric starting points stays stable.
- Which branches reach physical Sun mass depends on the starting point in the symmetry set: T_1/2 reaches epsilon = 1, T_0 does
  not, directly.
- Not in this paper: any Earth-Moon cycler, any Melnikov or splitting estimate, any quantitative fold location beyond "near
  epsilon = 1", and any treatment of the Sun's sense (the paper uses the physical sense; there is nothing to compare with the
  project's pre-#891 error). A published control for the sense is available through Table 4 (section 6).
- Novelty bearing: the mechanism (continuation in Sun mass from a commensurate-period three-body orbit, symmetric-phase
  equivalents, folds back to zero) is in print for these orbit families. See also Boudad, Howell and Davis (2020), digest filed
  alongside.

## 8. Slips noted (factual, minor)

- Reference list: the Boudad, Howell and Davis (2020) entry gives the pages as "2194-2294"; the paper runs 2194-2214.
- Table 1 prints `omega_S` as -0.925195985 (negative) while the text of section 3 writes `|omega_S|` and "approximately 29.5
  days"; these are consistent (the modulus gives 29.5289 d with TU = 4.34811305 d). No error.

## 9. References cited by Oshima (full list, p1334-1335)

Earlier periodic-orbit work on the bicircular or related problem, as cited in the introduction and sections 2, 3: Gomez, Jorba,
Masdemont and Simo (2001) and Boudad, Howell and Davis (2020) for the continuation from the CR3BP to the BCR4BP; Simo, Gomez, Jorba
and Masdemont (1995) for the BCR4BP itself; Oshima (2022a) for the symmetric-orbit method; Campagnola, Lo and Newton (2008) for the
symmetric ER3BP method; Broucke (1969), Macris, Katsiaris and Goudas (1975), Ichtiaroglou and Voyatzis (1990), Palacian, Yanguas,
Fernandez and Nicotra (2006), Peng and Xu (2015), Ferrari and Lavagna (2018), Fitzgerald and Ross (2022), Voyatzis, Tsiganis and
Gaitanas (2018) for the ER3BP.

- Belbruno, E., Miller, J., 1993. Sun-perturbed Earth-to-Moon transfers with ballistic capture. J. Guid. Control Dyn. 16, 770-775. doi:10.2514/3.21079.
- Boudad, K.K., Howell, K.C., Davis, D.C., 2020. Dynamics of synodic resonant near rectilinear halo orbits in the bicircular four-body problem. Adv. Space Res. 66, 2194-2294 [sic, 2194-2214]. doi:10.1016/j.asr.2020.07.044.
- Broucke, R.A., 1968. Periodic orbits in the restricted three-body problem with earth-moon masses. JPL technical report 32-1168, Pasadena.
- Broucke, R.A., 1969. Stability of periodic orbits in the elliptic, restricted three-body problem. AIAA J. 7, 1003-1009. doi:10.2514/3.5267.
- Campagnola, S., Lo, M., Newton, P., 2008. Subregions of motion and elliptic halo orbits in the elliptic restricted three-body problem. 18th AAS/AIAA Space Flight Mechanics Meeting, AAS 08-200, Galveston, Texas.
- Ferrari, F., Lavagna, M., 2018. Periodic motion around libration points in the elliptic restricted three-body problem. Nonlinear Dyn. 93, 453-462. doi:10.1007/s11071-018-4203-4.
- Fitzgerald, J., Ross, S.D., 2022. Geometry of transit orbits in the periodically-perturbed restricted three-body problem. Adv. Space Res. doi:10.1016/j.asr.2022.04.029.
- Gomez, G., Jorba, A., Masdemont, J., Simo, C., 2001. Dynamics and Mission Design Near Libration Points, Vol II: Fundamentals: The Cases of Triangular Libration Points. World Scientific, Singapore.
- Howell, K.C., Breakwell, J.V., 1984. Almost rectilinear halo orbits. Celest. Mech. 32, 29-52. doi:10.1007/BF01358402.
- Ichtiaroglou, S., Voyatzis, G., 1990. On the effect of the eccentricity of a planetary orbit on the stability of satellite orbits. J. Astrophys. Astr. 11, 11-22. doi:10.1007/BF02728017.
- Keller, H.B., 1977. Numerical solution of bifurcation and nonlinear eigenvalue problems. In: Rabinowitz, P. (Ed.), Applications of Bifurcation Theory. Academic Press, New York.
- Koon, W.S., Lo, M.W., Marsden, J.E., Ross, S.D., 2011. Dynamical Systems, the Three-Body Problem and Space Mission Design. Marsden Books, Wellington.
- Lara, M., Russell, R., Villac, B.F., 2007. Classification of the distant stability regions at Europa. J. Guid. Control Dyn. 30, 409-418. doi:10.2514/1.22372.
- Macris, G., Katsiaris, G.A., Goudas, C.L., 1975. Doubly-symmetric motions in the elliptic problem. Astrophys. Space Sci. 33, 333-340. doi:10.1007/BF00640102.
- Morais, M.H.M., Namouni, F., 2019. Periodic orbits of the retrograde coorbital problem. MNRAS 490, 3799-3805. doi:10.1093/mnras/stz2868.
- Oshima, K., 2021. Capture and escape analyses on planar retrograde periodic orbit around the Earth. Adv. Space Res. 68, 3891-3902. doi:10.1016/j.asr.2021.07.012.
- Oshima, K., 2022a. Continuation and stationkeeping analyses on planar retrograde periodic orbits around the Earth. Adv. Space Res. 69, 2210-2222. doi:10.1016/j.asr.2021.12.020.
- Oshima, K., 2022b. 3D stable and weakly unstable periodic orbits around the Earth near the retrograde co-orbital resonance with the Moon. Astrophys. Space Sci. 367, 42. doi:10.1007/s10509-022-04071-4.
- Palacian, J.F., Yanguas, P., Fernandez, S., Nicotra, M.A., 2006. Searching for periodic orbits of the spatial elliptic restricted three-body problem by double averaging. Physica D 213, 15-24. doi:10.1016/j.physd.2005.10.009.
- Peng, H., Xu, S., 2015. Stability of two groups of multi-revolution elliptic halo orbits in the elliptic restricted three-body problem. Celest. Mech. Dyn. Astr. 123, 279-303. doi:10.1007/s10569-015-9635-2.
- Pushparaj, N., Baresi, N., Ichinomiya, K., Kawakatsu, Y., 2021. Transfers around Phobos via bifurcated retrograde orbits: Applications to Martian Moons eXploration mission. Acta Astronaut. 181, 70-80. doi:10.1016/j.actaastro.2021.01.016.
- Russell, R.P., 2006. Global search for planar and three-dimensional periodic orbits near Europa. J. Astronaut. Sci. 54, 199-226. doi:10.1007/BF03256483.
- Simo, C., Gomez, G., Jorba, A., Masdemont, J., 1995. The bicircular model near the triangular libration points of the RTBP. In: Roy, A.E., Steves, B.A. (Eds.), From Newton to Chaos. Springer, Boston.
- Sliz-Balogh, J., Barta, A., Horvath, G., 2018. Celestial mechanics and polarization optics of the Kordylewski dust cloud in the Earth-Moon Lagrange point L5-I. Three-dimensional celestial mechanical modelling of dust cloud formation. MNRAS 480, 5550-5559. doi:10.1093/mnras/sty2049.
- Szebehely, V., 1967. Theory of Orbits: The Restricted Problem of Three Bodies. Academic Press Inc, New York.
- Topputo, F., 2013. On optimal two-impulse Earth-Moon transfers in a four-body model. Celest. Mech. Dyn. Astr. 117, 279-313. doi:10.1007/s10569-013-9513-8.
- Uesugi, K., 1996. Results of the MUSES-A "HITEN" mission. Adv. Space Res. 18, 69-72. doi:10.1016/0273-1177(96)00090-7.
- Voyatzis, G., Tsiganis, K., Gaitanas, M., 2018. The rectilinear three-body problem as a basis for studying highly eccentric systems. Celest. Mech. Dyn. Astr. 130, 3. doi:10.1007/s10569-017-9796-2.
- Whitley, R., Martinez, R., 2016. Options for staging orbits in cislunar space. 2016 IEEE Aerospace Conference, pp. 1-9, Big Sky, Mar. 2016.
- Zimovan, E.M., Howell, K.C., Davis, D.C., 2020. Near rectilinear halo orbits and nearby higher-period dynamical structures: orbital stability and resonance properties. Celest. Mech. Dyn. Astr. 132, 28. doi:10.1007/s10569-020-09968-2.

(Some DOIs and page ranges above were read at page-image resolution; treat a single digit as unverified until checked against the
publisher record.)
