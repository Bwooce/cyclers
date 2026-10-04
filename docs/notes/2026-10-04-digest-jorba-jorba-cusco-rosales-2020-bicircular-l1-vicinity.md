# Digest: Jorba, Jorba-Cusco & Rosales (2020), "The vicinity of the Earth-Moon L1 point in the bicircular problem"

Celestial Mechanics and Dynamical Astronomy 132:11 (2020), DOI 10.1007/s10569-019-9940-2, 25 pages (received 17 April 2019,
accepted 27 November 2019, published 7 February 2020).
Filed in the private paper corpus as
jorba-jorba-cusco-rosales-2020-vicinity-earth-moon-l1-point-bicircular-problem-cmda-132-11-doi-10.1007-s10569-019-9940-2.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with page and section or equation) or INFERRED (our
reading, or a comparison with project code). Page numbers are the journal's printed "Page n of 25". Numbers were transcribed
from the page and cross-checked against the PDF text layer; printed digits are given in full. Values read off a plotted
figure are marked "graph read" and are good to a few units in the last plotted digit only.
Context: tasks #891, #892 and #884. The bicircular module `core/bcr4bp.py` was corrected on 2026-10-04 (Sun sense) and this
paper is a published positive control for it (the ledger entry for #891 records the first use of this paper).

## 0. What the paper is

READ (abstract, p1): "The bicircular model is a periodic time-dependent perturbation of the Earth-Moon restricted three-body
problem that includes the direct gravitational effect of the Sun on the infinitesimal particle. In this paper, we focus on
the dynamics in the neighbourhood of the L1 point of the Earth-Moon system. By means of a periodic time-dependent reduction
to the centre manifold, we show the existence of two families of quasi-periodic Lyapunov orbits, one planar and one
vertical. The planar Lyapunov family undergoes a (quasi-periodic) pitchfork bifurcation giving rise to two families of
quasi-periodic halo orbits. Between them, there is a family of Lissajous quasi-periodic orbits, with three basic
frequencies."

READ: the paper does not give any table of initial conditions, any table of frequencies other than the two numbers in
section 3.2, and no bifurcation value in physical units. Its only tables are Table 1 (constants). Everything else is in text
and Figures 2 to 7. In particular the pitchfork bifurcation is located only by comparing four energy levels (Fig. 4).

## 1. The bicircular Hamiltonian as printed

### 1.1 Section 1, p3 (no homotopy parameter)

READ. Momenta defined on p2: "Defining the momenta as px = x' - y, py = y' + x and pz = z'" (primes for time derivatives).
Hamiltonian, transcribed from p3:

    H_BCP = (1/2)(px^2 + py^2 + pz^2) + y px - x py
            - (1 - mu)/r_PE - mu/r_PM - m_S/r_PS - (m_S/a_S^2)(y sin(theta) - x cos(theta))

with, quoted from p3: "r_PE^2 = (x - mu)^2 + y^2 + z^2, r_PM^2 = (x - mu + 1)^2 + y^2 + z^2, r_PS^2 = (x - x_S)^2 +
(y - y_S)^2 + z^2, x_S = a_S cos(theta), y_S = -a_S sin(theta) and theta = omega_S t."

### 1.2 Section 3, p11 (with homotopy parameter epsilon)

READ. Same Hamiltonian with epsilon multiplying the Sun terms: `- epsilon m_S/r_PS - (epsilon m_S/a_S^2)(y sin(theta) -
x cos(theta))`, same definitions "x_S = a_S cos(theta), y_S = -a_S sin(theta), theta = omega_S t, and omega_S is the mean
angular velocity of the Sun in these synodic coordinates". Epsilon = 0 is the RTBP and epsilon = 1 the BCP (p12).

### 1.3 Frame, Sun position and sense

READ (the printed Hamiltonian is internally consistent here, unlike the 2018 paper's eq. 3 where the printed r_PS and
indirect term disagree on the sign of y_S):
- Earth is at x = mu and the Moon at x = mu - 1, from `r_PE^2 = (x - mu)^2 + ...` and `r_PM^2 = (x - mu + 1)^2 + ...`. INFERRED
  from those expressions. Fig. 1 (p3) agrees: E at the origin of the sketch, M at negative x, L1 between them, L2 beyond M.
  This is the mirror image (rotation by pi about z) of the project's convention (Earth at -mu, Moon at 1 - mu).
- Sun at (a_S cos(theta), -a_S sin(theta)) with theta = omega_S t and omega_S > 0, so the Sun moves clockwise in the synodic
  frame. READ: Fig. 1 draws S at lower right with the arrow for theta running clockwise from +x.
- The printed indirect term is consistent with this Sun. INFERRED (our check): the first-order expansion of -m_S/r_PS about
  the origin for a Sun at (a_S cos(theta), -a_S sin(theta)) is `-m_S/a_S - (m_S/a_S^2)(x cos(theta) - y sin(theta))`; the
  printed term `-(m_S/a_S^2)(y sin(theta) - x cos(theta)) = +(m_S/a_S^2)(x cos(theta) - y sin(theta))` cancels it, which is
  the requirement for the Sun's acceleration of the Earth-Moon barycentre.
- Theta = 0 puts the Sun on the +x axis of this paper's frame, which is the side away from the Moon (Moon at negative x).
  INFERRED. In the project's frame (Earth at -mu, Moon at +1 - mu) the same instant has the Sun at (-a_S, 0), that is
  `theta_sun0 = pi` in the convention `theta_S = theta_S0 - omega_S t` of `core/bcr4bp.py`.

### 1.4 Table 1 (p3), constants, all digits as printed

READ. Caption: "Floating point values used for the different constants of the BCP".

| Symbol | Value as printed |
| --- | --- |
| mu | 0.012150581623433623 |
| m_S | 328900.54999999906 |
| omega_S | 0.925195985518289646 |
| a_S | 388.81114302335106 |

Units are those of the Earth-Moon RTBP (section 1, p2: "it is usual to take the same units and reference frame as in the
RTBP"): length = Earth-Moon distance, time unit = 1/n of the Earth-Moon mean motion, Earth + Moon mass = 1. The Sun's period
in the frame is `T_S = 2 pi / omega_S` (section 3, p12). INFERRED (arithmetic): `T_S = 6.791193871923018` time units, and
`1 - omega_S = 0.0748040144817` is the Sun's inertial mean motion in these units, which corresponds to a period of about
365.2 days with a time unit of about 4.3484 days.

Comparison with project code (INFERRED, our check):
- `core/qbcp.py` `_QBCP_MU_EM`, `_QBCP_MU_S`, `_QBCP_A_S` match Table 1 digit for digit. `_QBCP_OMEGA_S = 0.92519598551829646`
  has 17 significant digits whereas Table 1 prints `0.925195985518289646` (18 significant digits, the other three constants
  having 17). The two differ by 6.8e-15 and carry an apparent extra digit "8" in the printed value; this is immaterial
  numerically. Which is the double actually used by the authors cannot be told from the page.
- `core/bcr4bp.py` uses the Andreu-digest constants, not Table 1: mu = 0.0121505816 (Table 1: 0.012150581623433623),
  m_S = 328900.5423094043 (Table 1: 328900.54999999906, a relative difference of 2.3e-8), omega_S = 0.925195985520347 (Table 1:
  0.925195985518289646, a difference of 2.1e-12), a_S the same to printed precision. Any comparison at the 1e-8 level and
  below with this paper must use Table 1 and not the module defaults. The ledger's published positive control agrees
  with the paper's second frequency to 15 digits, so the test evidently already does this or is insensitive to it; check
  which before quoting a tolerance.

### 1.5 The "indirect effect"

READ, p2 (section 1): "Sun's gravitational acceleration upon the particle is one of the most relevant forces ignored by the
RTBP... This effect is called direct effect of Sun's gravity. There is, however, another effect of Sun's gravity on the
particle, the indirect one: the gravity of Sun changes the motion of Earth and Moon; therefore, the motion of the test
particle suffers a small deviation according to the new trajectories of Earth and Moon. This effect is especially important
near Earth and Moon. This work does not consider the indirect effect."
READ, p3: "Note that this motion does not follow Newton's laws, since we are not taking into account the effect of Sun on the
motion of Earth and Moon." READ, p20 (section 5): L2 was not studied "because the motion near this point is severely affected
by the indirect effect of Sun's gravity on the particle. A better suited model such as the quasi-bicircular problem (a
coherent version of the BCP) should be used to investigate this point."
INFERRED: the paper's "indirect effect" is the Sun's perturbation of the Earth-Moon orbit (the part the QBCP adds), not the
`m_S/a_S^2` acceleration-of-the-origin term in the Hamiltonian, which the paper keeps. The two uses of the word "indirect"
(this one and the project's "indirect term" in `core/bcr4bp.py`) are different things; the module's incoherent model matches
the paper's BCP, which includes the origin-acceleration term and omits the Sun's effect on the primaries' orbit.

## 2. Section 3.1 and 3.2: the periodic orbit that replaces L1

### 2.1 What is printed (READ)

p12, section 3.1: "Due to the periodic perturbation due to Sun, the Lagrangian points are no longer equilibria; they are
replaced by periodic orbits with the same period as Sun's (T_S = 2 pi/omega_S). We name these replacements as dynamical
equivalents of the Lagrangian points."
p12: "Due to the high instability of L1, we have applied a multiple shooting technique combined with a continuation method to
go from epsilon = 0 to epsilon = 1. When epsilon reaches 1, the replacement is a small unstable periodic orbit with the same
normal behaviour as L1 (four elliptic directions and two hyperbolic ones). The size of the orbit is around 10^-3, and its
(x, y) projection revolves L1 twice in T_S units of time (Fig. 2). The linear normal behaviour is of type
saddle x centre x centre."
p12: "The unstable eigenvalue of the monodromy matrix is large, around 10^8." The explanation printed: "the eigenvalue of the
dynamical replacement is, at first order, the exponential eigenvalue of L1 in the RTBP multiplied by the period of Sun."
p4 (section 1.1): "the monodromy matrix around the periodic orbit that replaces L1 in the BCP has an hyperbolic eigenvalue
close to 4.287 x 10^8."
p13 (section 3.2): eigenvalue set `(lambda_h, lambda_h^-1, lambda_e1, lambda_e1^-1, lambda_e2, lambda_e2^-1)` with
`|lambda_h| >> 1` and `|lambda_e1| = |lambda_e2| = 1`.

### 2.2 Normalised logarithms and frequencies (READ, p13, section 3.2)

"The normalized logarithms are given by alpha_1 = log lambda_h, omega_1,2 = log lambda_e1,2. Here, omega_1,2 are to be
understood as a complex logarithms of the elliptic eigenvalues. As mentioned, any combination +-(omega_i + k omega_S) for
i = 1, 2 and k in Z is also an admissible choice. In Remark 2.3, we discuss the optimal choice for these logarithms in terms
of the decay of the Fourier series representing each entry of the Floquet change. In particular, we are interested in the
change of variables which is as close as possible to constant coefficients. This is obtained, in this case, selecting
omega_1 = 2.32981963603288 and omega_2 = 2.26695149158478, which are close to the frequencies related to the equilibrium
point L1 in the RTBP (2.33438585628816 and 2.2688310655411, respectively)."

Definition (READ, p7, footnote 1 and step 3): "By normalized, we mean that the logarithms are divided by the period T." The
eigenvalue of the monodromy matrix is `lambda = exp(i omega T)` with T = T_S here; alpha is defined by `lambda_h = exp(alpha T)`.
Remark 2.3 (p7 to p8): the optimal choice of the multiple of omega_S "is the one that makes the dominant Fourier
coefficients be at the beginning of the Fourier series"; "In problems which are a perturbation of an autonomous one, we know
in advance that the logarithms are to be chosen as close as possible to the frequencies of the dynamical equivalent of the
periodic orbit in the autonomous system."

Real Floquet normal form of the second-order part (p13): `H_2 = alpha_1 x px + omega_1 (y^2 + py^2)/2 + omega_2 (z^2 + pz^2)/2`
(the printed text of the last term carries `z^3 + pz^3` in the typesetting of the display on p13; the second-order meaning
is clear from the line above it on the same page, which has `z^2 + pz^2`). Reported as a typesetting slip, not a defect.
So omega_1 belongs to the in-plane elliptic direction (y, py) and omega_2 to the vertical direction (z, pz). INFERRED from the
form of H_2 and from the comparison numbers: 2.33438585628816 is the planar elliptic frequency of L1 in the RTBP and
2.2688310655411 the vertical one.

INFERRED (our check, short script, Table 1 mu, high-precision root of the L1 collinear equation): L1 at x = 0.836915145386502
from the Earth-Moon barycentre side (project frame); planar elliptic frequency 2.33438585398 and vertical frequency
2.26883106319. These agree with the printed RTBP values to eight to nine significant digits but differ from them by about
2.3e-9 and 2.4e-9 respectively, in the same direction. The cause is not known (rounding of mu or a loose L1 root in the
authors' program are both possible); it does not matter for the use below, but a test that quotes the printed RTBP
values should allow 1e-8.

### 2.3 Cross-check with the project's ledger entry (INFERRED, our comparison)

The `#891` ledger entry records "4.287389e8, 2.3298196303 and 2.266951491584771" from the corrected module. Against the printed
omega_1 = 2.32981963603288 the entry's 2.3298196303 differs by 5.7e-9 (the digits run 2.329819630 against 2.329819636), while
omega_2 agrees to all printed digits. Either the ledger entry has a transcription slip in the ninth decimal of omega_1
(2.3298196303 against 2.32981963603), or the module's omega_1 differs from the paper's at 2.5e-9 relative, which would be
worth a look because omega_2 is matched so closely. The permanent test should be inspected to settle which.

### 2.4 Other numbers in section 3 and the figure

- Fig. 2a (p12, graph read): continuation of the replacement orbit in epsilon, x against epsilon, a nearly straight line from
  about x = -0.8369 at epsilon = 0 to about x = -0.8377 at epsilon = 1. The epsilon = 0 value is the RTBP L1, which in this
  paper's frame (Moon at x = mu - 1) is at x = -0.836915 (INFERRED from our recomputation above), so the plot's start is
  consistent. Axis tick labels run from -0.8377 to -0.8369.
- Fig. 2b (p12, graph read, tick labels small and not fully certain): the (x, y) projection is one closed curve, roughly an
  ellipse, with y between about -0.006 and +0.006 and x spanning most of an axis labelled from about -0.8378 to -0.8358.
  The text says the size is "around 10^-3"; the plotted extent is a few times 1e-3 in x and about 1.1e-2 in y. READ, both
  statements as printed; they are not contradictory if "around 10^-3" is an order-of-magnitude statement. It is the
  geometry to compare with, not the text.
- Section 3.3 (p13): the Hamiltonian is scaled by gamma, "the distance between Moon and the average of the periodic
  orbit", so that this distance is one. Section 4.1 (p15): "this unit corresponds to the average distance from the periodic
  orbit that replaces L1 to the Moon, which is 57,355 km." INFERRED: 57,355 km over 384,400 km is 0.14921 Earth-Moon
  distances, against 0.15093 for the RTBP L1 (58,019 km), so the orbit's time-average is about 1.1 percent closer to the Moon
  than the RTBP point (the plotted shift of x from -0.8369 to -0.8377 in Fig. 2a has the same sign). The Earth-Moon
  distance in km used for the conversion is not printed.

## 3. Sections 3.3 to 4: normal form, centre manifold, the numbers that exist

### 3.1 Method in a few lines (READ, sections 2 to 4)

1. Autonomise the periodic Hamiltonian (add the angle theta and its action I_theta), translate the periodic orbit to the
   origin and scale by gamma (section 3.3).
2. Symplectic Floquet change (Theorem 2.1, p6) turns the second-order terms into constant-coefficient real Floquet normal form
   `H_2` above. The periodic coefficients of the change are stored as Fourier series; the monodromy matrix is computed in
   extended (mpfr) precision (appendix A.2, p23).
3. Lie-series normal form (sections 2.2 to 2.4) removes both the hyperbolic part (set to zero) and the time dependence, up to
   order 12 for the figures (Fig. 4 caption: "The expansion used for the Hamiltonian is of order 12"; Fig. 3 runs orders 4
   to 16). The result is an autonomous 2-degree-of-freedom Hamiltonian in (q1, p1, q2, p2), so its energy `h` is conserved and
   Poincare sections at fixed `h` are area-preserving maps.
4. Return to synodic coordinates by composing the changes of variable; each periodic orbit of the reduced Hamiltonian gives a
   quasi-periodic orbit of the BCP, with the Sun's frequency added to the frequencies of the reduced orbit.
What this means for numbers: `h` is the energy of the reduced Hamiltonian in normal-form units (the "normalized energy"); it
has no simple counterpart in the BCP's own Hamiltonian, so the printed `h` values (0.1, 0.2, 0.4, 0.5, 0.7, 0.9) cannot be
reproduced without rebuilding the normal form. INFERRED.

### 3.2 Energy levels, sections and what is shown (READ, section 4.1, p16 to p18)

- Normal-form coordinates (p16): "we name q1, p1, q2 and p2 the coordinates in the normal form". Section `h = {q2 = 0}`
  "corresponds, at first order, to fix z = 0 in the synodical coordinates. We will name this section as the horizontal one."
  Section `v = {q1 = 0}` is the vertical one (p17).
- Figs. 4 and 5 (pp16, 18): sections at normalised energy `h` = 0.2, 0.5, 0.7, 0.9, expansion of order 12. Axes: Fig. 4
  (horizontal section) q1 horizontal and q3 vertical (label q3 as printed; the text names the coordinates q1, p1, q2, p2, so
  q3 is probably p1 or q2 under a different naming, INFERRED, reported as a labelling slip); Fig. 5 (vertical section) q2 and q4.
- Statements (p17), quoted: "The periodic orbit replacing L1 is at the origin in the centre manifold coordinates; it is
  totally elliptic and has zero energy." "In Fig. 4, the outer limit of the plots corresponds to a planar Lyapunov orbit, and
  the fixed point at the centre to a vertical Lyapunov orbit." "for the energy level h = 0.2 (the plot (a)) the invariant
  curves are Lissajous orbits. When the energy level is increased to h = 0.5 (plot (b)), the planar family of Lyapunov orbit
  goes through a pitchfork bifurcation and the halo family of periodic orbits appears." So the halo bifurcation lies between
  `h` = 0.2 and `h` = 0.5. The graph reads: Fig. 4a (h = 0.2) shows one nest of curves about the centre, Fig. 4b (h = 0.5) two
  extra islands left and right at about q1 = +-0.6, Fig. 4c and 4d show the islands growing with a hyperbolic region between.
  No exact `h` for the bifurcation is printed.
- p17 to p18: "Going back to the original synodical coordinates, the family of Lyapunov periodic orbits becomes a (Cantor) family of
  quasi-periodic solutions (2D tori) with two basic frequencies, by adding the frequency of the Sun to the frequencies of the
  family. Of course, the periodic orbits whose frequency is (close to be) in resonance with that of the Sun are destroyed and
  this gives the Cantorian structure to the family." Same for halos. "Therefore, in synodical coordinates, the family of
  quasi-periodic Lyapunov orbits undergoes a pitchfork bifurcation giving rise to quasi-periodic halo orbits."
- Families, how many (READ): planar quasi-periodic Lyapunov (one family, the outer limit of Fig. 4); vertical quasi-periodic
  Lyapunov (one family); two quasi-periodic halo families (north and south; "two families of quasi-periodic halo orbits" in
  the abstract; the two island centres in Fig. 4b to 4d); Lissajous orbits with three frequencies for `h` below the
  bifurcation. Parameterisation: all by the normalised energy `h` of the reduced Hamiltonian.
- Fig. 6 (p19): a quasi-periodic halo orbit at `h = 0.4` (graph read: x from about -0.86 to -0.82, y between about -0.08 and
  +0.08 in the x-y projection; in the y-z panel the vertical axis is labelled z but runs from about -1 to -0.7, centred near
  the x of L1, so it is probably x mislabelled, INFERRED, reported as a labelling slip). It is computed as a fixed point of
  the section map by Newton's method, then sent to synodic coordinates.
- Fig. 7 (p19): a Lissajous orbit at `h = 0.1` (graph read: x from about -0.843 to -0.832, y within about +-0.015, z within
  about +-0.03).
- No rotation numbers, no frequency table and no stability data for the tori are printed.

### 3.3 Radius of convergence (READ, section 4.1, p14 to p15, Fig. 3)

Estimator printed: `r_n^-1 = (||H_n||_1)^(1/n)`, `||H_n||_1 = sum over |k| = n of |h_k|`, `3 <= n <= N`. Graph read, Fig. 3: for
orders 4 to 9 the autonomous (time-dependence removed) and non-autonomous curves agree, rising from about 0.77 at order 4 to a
peak of about 1.10 at order 6, then about 0.93 at order 9, peaking again at about 1.05 at order 10; after order 10 the
autonomised curve falls to about 0.56 at order 13, 0.72 at order 14, 0.55 at order 15 and 0.42 at order 16, while the
non-autonomous curve stays between about 0.87 and 0.98. The text says "It stays around 0.9 at order 16." Units: 57,355 km.

## 4. Statements about the model's applicability worth remembering (READ)

- p20: "the structure of the phase space is similar to the one observed in the RTBP. In particular, the bifurcation that gives rise
  to the halo orbit (well known in the RTBP) has its quasi-periodic counter part in the BCP."
- p20: L2 not studied because of the indirect effect; the QBCP "should be used to investigate this point."
- The validation of the software (A.3, p23 to p24): the difference between the reduced flow sent to synodic coordinates and
  direct numerical integration of the BCP is checked to behave as `h^(m+1)` for truncation order m, "passed successfully in
  all the cases" (no numbers printed).

## 5. Positive controls for the project

The only high-precision published numbers are those of section 3.2 (the second frequency to 15 digits), and the project has
already used them. Beyond that, the paper supports mostly structural controls. Each recipe below is for `core/bcr4bp.py`
(project frame: Earth at -mu, Moon at 1 - mu, Sun clockwise) and uses Table 1 constants, not the module defaults.

1. Unstable eigenvalue 4.287 x 10^8 (printed to four digits only). Recipe: refine the T_S-periodic orbit near L1 by
   multiple shooting with epsilon continuation from the CR3BP point; monodromy over one T_S; take the largest eigenvalue.
   Tolerance: the printed value supports about 1e-4 relative. (Used already: 4.287389e8.) The cheap consistency check on the
   same orbit is `ln(lambda_h)/T_S`: for lambda_h = 4.287389e8 and T_S = 6.791193871923 this is 2.92678, against 2.93206 for the
   RTBP hyperbolic exponent at L1 (our computation); the paper's "at first order" statement allows a fraction of a percent.
   INFERRED.
2. omega_2 = 2.26695149158478 (matched to 15 digits already) and omega_1 = 2.32981963603288 (the ledger entry shows
   2.3298196303; see section 2.3 above). Recipe: eigenvalues of the monodromy matrix, `omega = Im(log lambda)/T_S`, then
   choose the branch `+-(omega + k omega_S)` closest to the RTBP frequencies 2.33438585628816 and 2.2688310655411.
   Because the freedom is exactly `omega_S` (0.9252), any frequency within 0.5 of the RTBP value is unambiguous.
3. Sense discriminator. The paper's two frequencies differ from the old-sense computation (2.33046 and 2.26711 per the
   ledger) at the 1e-4 level, so these numbers pin the Sun's sense; they would not distinguish a wrong phase theta_S0
   (the T_S-periodic orbit is the same curve for any phase, only the time origin moves). INFERRED.
4. Two revolutions about L1 per T_S: the (x, y) projection of the replacement orbit makes two loops. Recipe: count winding of
   (x - x_L1, y) over one period. Qualitative, robust, and also sense-dependent.
5. Fig. 2a start and end: x(epsilon = 0) = -0.8369 (RTBP L1, exact to the digits shown) and x(epsilon = 1) about -0.8377
   at the point of the orbit at t = 0 (graph read, good to about 1e-4). Recipe: in the project frame x = +0.8369 and about
   +0.8377 with the Sun at (-a_S, 0) at t = 0 (theta_sun0 = pi); the sign of the shift (towards the Moon) is the useful
   content. INFERRED that the plotted point is the t = 0 point of the orbit; a plotted maximum or the orbit's mean would
   change the comparison, so use the sign and about 1e-4 only.
6. Orbit extent: x about 2e-3 wide, y about 1.1e-2 tall (Fig. 2b, graph read; the figure's tick labels are small, so check against
   the image before using). Recipe: extent of the refined orbit in each coordinate. Loose.
7. Mean distance to the Moon of the replacement orbit: 57,355 km, to be compared with the time average of the orbit's
   position (the paper says gamma is the distance between the Moon and the average of the periodic orbit). Needs the
   Earth-Moon distance in km used by the authors (not printed; 384,400 km gives 0.14921). Loose, about 1e-3 at best.
8. Table 1 constants as a self-consistency check on `core/qbcp.py` (they match) and on `core/bcr4bp.py` (they do not;
   the module uses Andreu's values, listed in section 1.4).
9. RTBP L1 frequencies 2.33438585628816 and 2.2688310655411 as controls on `core/cr3bp.py` linearisation: our recomputation gives
   2.33438585398 and 2.26883106319, so the controls hold only to about 1e-8.
10. The `epsilon = 0` end of the continuation: the periodic orbit at epsilon = 0 is the RTBP point L1 (a constant orbit),
    so the continuation in epsilon from the CR3BP fixed point is a model-independent test of the homotopy code.
Not testable without rebuilding a Lie-series normal form: the energy levels `h`, Figs. 3 to 7, the bifurcation `h`. A cheap
substitute for the bifurcation in the synodic frame would be to look for the halo family of quasi-periodic orbits by
continuation, but the paper prints no number to compare against.

## 6. References worth chasing from this paper

- Boudad, K., Howell, K., Davis, D.: Near rectilinear Halo orbits in cislunar space within the context of the Bicircular
  Four-Body Problem. IAA-AAS-SciTech-039 (2019). (Cited on p4 for resonant periodic orbits near L1 in the BCP.)
- Andreu, M.A.: The quasi-bicircular problem. Ph.D. thesis, Univ. Barcelona (1998); Andreu 2002, Celest. Mech. 84(2), 105-133,
  "Dynamics in the center manifold around L2 in the quasi-bicircular problem"; Andreu and Simo, "The quasi-bicircular
  problem for the Earth-Moon-Sun parameters", preprint (2000).
- Bihan, Masdemont, Gomez, Lizy-Destrez, Nonlinearity 30(8), 3040 (2017), invariant manifolds of the QBCP by the
  parameterisation method.
- Castella and Jorba, Celest. Mech. 76(1), 35-54 (2000), vertical families of 2D tori near the triangular points of the BCP.
- Simo, Gomez, Jorba, Masdemont (1995), "The Bicircular model near the triangular libration points of the RTBP", in From
  Newton to Chaos, pp. 343-370.
- Huang (1960), NASA TN D-501, and Cronin, Richards and Russell, Icarus 3, 423-428 (1964), for the origin of the BCP.
- Jorba (1999), Exp. Math. 8(2), 155-195, the normal-form software this paper adapts.

## 7. What this means for the project

- READ: this paper confirms, in its own words and in a self-consistent printed Hamiltonian, the Sun's sense (clockwise in the
  synodic frame), the presence of the indirect origin-acceleration term, and the Table 1 constants. It is therefore an
  independent source for the #891 correction in addition to the 2018 paper.
- INFERRED: the BCP is the right model near L1 and the QBCP near L2 (the paper's own conclusion, p20 and, from the 2018
  paper, p12). For #884 (Sun-forced Earth-Moon cyclers) this means the incoherent `core/bcr4bp.py` is an adequate L1-region
  model but a coarse one near the Moon's far side.
- No Earth-Moon cycler, resonant orbit with lunar flybys or transfer orbit appears in this paper. It concerns only the
  neighbourhood of the L1 replacement orbit (size of order 1e-3 to 1e-2 length units). Nothing here can be used to test a
  cycler.
- Wrongly-read numbers to avoid: the figures' energy `h` is not the BCP Hamiltonian value; the paper's `m_S` is not the
  value in `core/bcr4bp.py`; the paper's Moon is at negative x.
