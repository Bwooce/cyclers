# Digest: Henon and Heiles 1964, the applicability of the third integral of motion

Date: 2026-10-05 (Sydney). Reading and reasoning, plus one throwaway numerical check (section 6) run outside the repository
and not committed.

Source: M. Henon and C. Heiles, "The applicability of the third integral of motion: some numerical experiments",
Astronomical Journal 69:73-79 (February 1964), DOI 10.1086/109234, received 7 August 1963. Filed in the private paper corpus
as `henon-heiles-1964-applicability-third-integral-of-motion-aj-69-73-doi-10.1086-109234-ads-scan.pdf` (7 pages, ADS scan
with no text layer). I rendered all 7 pages at 130 dpi and read every page image, including Figs 1-12.

Evidence tags: READ (p.N) is read at printed page N. COMPUTED is my calculation. INFERRED is my reasoning. "Read off the
figure" is an approximate reading of a plotted point or curve, not a printed number.

Related: `docs/notes/2026-10-04-digest-cincotta-simo-1999-conditional-entropy.md` (uses this model at h = 0.118), the MEGNO
digest `docs/notes/2026-10-04-digest-cincotta-simo-2000-megno.md` and the FLI digest
`docs/notes/2026-10-04-digest-froeschle-lega-gonczi-1997-fli.md`; the consumer is `#924` in `data/OUTSTANDING.md`.

## 1. What the paper is

A numerical experiment on whether a third isolating integral exists for an axisymmetric galactic potential U_0(R, z). The
authors reduce the problem to motion in a plane in an arbitrary potential U(x, y) (a system with two degrees of freedom),
choose one simple potential, and map the surface of section (y, y-dot) at x = 0 for several energies. The result that made
the paper famous (READ pp.76-77, 79): at low energy every orbit lies on an invariant curve (a third integral exists
in practice); above a critical energy about 0.11 an "ergodic" region appears in which points fill an area at random, and its
share of the allowed area grows rapidly with energy until, near the energy of escape 1/6, it is almost all of it. It is the
origin of the Henon-Heiles model and of the standard surface-of-section picture of mixed phase space.

## 2. The potential, exactly as printed

READ (p.75, eq. 11): the potential chosen for study is

    U(x, y) = (1/2) ( x^2 + y^2 + 2 x^2 y - (2/3) y^3 )                                       (11)

The equations of motion (p.75, eq. 12) are

    x-ddot = -dU/dx = -x - 2 x y,       y-ddot = -dU/dy = -y - x^2 + y^2                      (12)

The energy is E = U(x, y) + (x-dot^2 + y-dot^2)/2 with unit mass (eqs 6, 7, p.74).

Confirmation of the standard form inferred in the Cincotta and Simo digest: that digest recorded phi = (q1^2 + q2^2)/2 +
q1^2 q2 - q2^3/3 as the standard Henon-Heiles potential, not READ. Expanding eq. (11): U = x^2/2 + y^2/2 + x^2 y - y^3/3.
This is identical, with (q1, q2) = (x, y). The inferred form is CONFIRMED by the source; no correction is needed. The
equations of motion agree term by term with the ones I used (COMPUTED: -dU/dx = -x - 2xy, -dU/dy = -y - x^2 + y^2).

Geometry stated by the paper (p.75): "Near the center [equipotential lines] tend to be circles; farther out they become
elongated in three directions. The particular equipotential U = 1/6 consists of three straight lines, forming an
equilateral triangle" (Fig. 2). Fig. 2 labels U = 0.0100, 0.0417, 0.0833, 0.1250, 0.1667. COMPUTED: the three saddles of U
are at (0, 1) and (+/- sqrt(3)/2, -1/2), each with U = 1/6 exactly, so the escape energy E = 1/6 (printed, p.77:
"E = 1/6 is the energy of escape") is the saddle energy. For E > 1/6 the equipotentials open and the star can escape
(p.77).

## 3. The surface-of-section method (pp.73-75)

Reduction (pp.73-74, eqs 1-10). The axisymmetric problem with the angular momentum constant C_2 is equivalent to planar
motion in U = U_0 + C_2^2/(2 R^2) (eq. 4); the phase space (x, y, x-dot, y-dot) has four dimensions and there is one known
integral, the energy (eq. 6). The section is the plane x = 0 crossed with x-dot > 0: the successive points P_1, P_2, ... of
the trajectory in the (y, y-dot) plane satisfy x = 0, x-dot > 0 (eq. 9, Fig. 1). Given E, a point (y, y-dot) fixes x-dot from
eq. (7), so the passage P_i to P_(i+1) is a mapping, area-preserving (p.75, citing Birkhoff 1927 and Moser 1962). The allowed
region of the section is U(0, y) + y-dot^2/2 <= E (eq. 10). If a second isolating integral exists the points lie on a curve;
if none, they fill an area. The criterion: "compute a number of points P_i, plot them in the (y, y-dot) plane and see whether
they lie on a curve or not" (p.74).

COMPUTED: the section's y range at E = 0.118 (y-dot = 0) is from -0.4284 to +0.6425, from the roots of U(0, y) = E. This
matches the C&S Fig. 1a abscissa range (0.3 to 0.62 shown). Integration (p.75): Runge-Kutta; two independent computations by
the two authors on different machines (CDC 1604 and IBM 7090, Adams and Runge-Kutta); the energy decreased "very slightly
(< |0.00003| for 150 orbits)". Whether that is per orbit or total is as printed; I read it as the drift over the 150 orbits
(INFERRED).

A single trajectory on a curve rotates around it with a constant angle alpha between O P_i and O P_(i+1) (p.75): for the Fig. 3
orbit at E = 0.08333, "its approximate value is alpha = 0.1143 (taking one revolution as the unit)". alpha is generally
irrational; if alpha = p/q the orbit is periodic. Chains of islands (q islands around a stable periodic orbit) are described
on p.76: invariant points at the middle of four small loops in Fig. 4 correspond to stable periodic orbits, and the "three
intersections of curves are also invariant points" (unstable periodic orbits).

## 4. Printed results

### 4.1 Energies studied

READ (pp.75-76, 77): E = 1/12 = 0.08333 (Fig. 3 one orbit, Fig. 4 the complete picture), E = 0.12500 (Fig. 5), E = 1/6 =
0.16667 (Fig. 6, the escape energy). Figure 7 covers a range of energies from about 0.01 to about 0.17 (p.76). The
mapping experiment of section 4 uses a = 1.6 (Fig. 8).

### 4.2 Qualitative statements (all READ)

- E = 0.08333 (p.75, p.76; Fig. 4): "the complete picture in the (y, y-dot) plane": every orbit shown lies on a closed curve.
  "These curves form a one-parameter family which fills completely the available area, defined by (10)." The boundary of the
  area is "almost identical with the outer curve on Fig. 4" (p.76). The picture shows a central region of closed curves, four
  invariant points in the middle of four small loops (stable periodic orbits) and three intersection points (unstable
  periodic orbits). The area is "completely covered with curves" (p.77).
- E = 0.12500 (Fig. 5, p.76): "We still have a set of closed curves around each stable invariant point. But these curves no
  longer fill the whole area." The isolated points all belong to one trajectory; "It is clearly impossible to draw any curve
  through them. They seem to be distributed at random, in an area left free between the closed curves." The change "seems to
  occur abruptly across some dividing line in the plane." A chain of five small loops on the right belongs to one trajectory
  jumping from loop to loop (a "chain of islands"). Open circles in the middle (an "eight-shaped" trajectory) are an
  intermediate kind between closed curves and ergodic behaviour (p.77).
- E = 1/6 = 0.16667 (Fig. 6, p.76): "All the isolated points correspond to one trajectory, and it is apparent that this
  'ergodic' trajectory covers almost the whole area." Two of the sets of closed curves of Fig. 5 (those on the y-dot axis) have
  disappeared, presumably because their central invariant point became unstable; the other two degenerated into chains of two
  islands. "The outer line on Fig. 6 is the limit given by (10)."
- Three hypotheses (p.76): (1) there is an infinite number of islands and chains; (2) the set of islands is dense everywhere;
  (3) the islands do not cover the area since they become very small; there is a "sea" between them in which the ergodic
  trajectory is dense.
- The critical energy (p.77): "Up to a critical energy (about E = 0.11) the curves cover the whole area; there is no ergodic
  orbit. For higher energies the area covered by curves shrinks very rapidly. Thus the situation could be very roughly
  described by saying that the second integral exists for orbits below a 'critical energy,' and does not exist for orbits
  above that energy." "In the present case the critical energy is less than the energy of escape."
- Conclusions (p.78): "If the energy is small, it seems that a third isolating integral always exists ... If the energy is
  higher than the critical energy, there are an infinite number of separated regions in the phase space where such a third
  integral still seems to exist. The space left free between these regions is the 'ergodic region' ... The proportion of
  allowable phase occupied by this ergodic region increases very rapidly and tends to be the whole space." (pp.78-79).

### 4.3 The area-fraction method and Fig. 7 (p.77)

Method (READ p.77): to decide whether a point P_1 lies on a curve or in an ergodic orbit, a second point P_1' is taken "very
close to P_1 (usually at a distance 10^-7)"; then "usually 25" successive transforms of both are computed. On a curve the
distance P_i P_i' increases only slowly (about linearly with i); in the ergodic region roughly exponentially. The quantity

    mu = sum_(i=1..25) (distance P_i P_i')^2                                                  (13)

is computed, and P_1 is in the ergodic region if mu > mu_c, on a curve if mu < mu_c. The values of mu "covered a very wide
range, from about 10^-12 to 10^+1 [illegible exponent: printed digits read as 10^-12 to 10^+1]" and the criterion is very
sensitive, the exact mu_c unimportant: "Here mu_c approx 10^-4" (p.77). So this is itself a finite-time divergence indicator,
a close cousin of FLI.

Fig. 7 (p.76 caption): "Relative area covered by the curves as a function of energy". Axes: relative area 0 to 1.0,
energy 0 to 0.18. READ qualitative: relative area 1.0 for energies up to the critical energy about 0.11, then a rapid
fall. Dots read off the figure, approximate (about +/- 0.03 in area and +/- 0.003 in energy; the figure has no numeric
labels on the points):

| Energy | Relative area covered by curves (read off Fig. 7) |
| --- | --- |
| 0.01 to about 0.11 (about 8 dots) | 1.0 |
| about 0.118 | about 0.9 |
| 0.125 | about 0.7 |
| about 0.14 | about 0.5 |
| about 0.15 | about 0.2 |
| about 0.16 | about 0.17 |
| about 0.167 (escape) | about 0.1 |

The caption, the plateau at 1.0 and the sharp fall are printed; the dots are my readings. No numeric table is printed, and
the paper states no fraction at any energy in the text (READ absence). Note the dot at 0.125 (about 0.7) is consistent with
the paper's description of Fig. 5 ("no longer fill the whole area" but closed curves still occupy much of it) and the dot
at 0.1667 (about 0.1) with Fig. 6 ("covers almost the whole area" by the ergodic trajectory); both are qualitative.

### 4.4 The mapping experiment (section 4, pp.77-78, Figs 8-12)

The area-preserving mapping (p.78, eq. 14):

    X_(i+1) = X_i + a (Y_i - Y_i^3),     Y_(i+1) = Y_i - a (X_(i+1) - X_(i+1)^3)              (14)

with a constant (a = 1.6 in Fig. 8). A cubic mapping suggested by Dr. Kruskal as a quicker substitute, a factor about 1000
cheaper than integrating orbits (p.77). Fig. 8: a central region of simple closed curves around the stable invariant point
X = Y = 0, "a chain of six islands (instead of five)" and an outer "ergodic" region; up to 10^5 points were computed for some
curves. Figs 9 to 12: initial points on a grid (0.02 for Fig. 9, 0.002 for Figs 10 and 12, 0.0002 for Fig. 11), 1000
iterations, a point marked "nonergodic" if all 1000 iterates stay in the vicinity of the origin (X^2 + Y^2 < 100; p.78).
Enlargements: area A (ten times), then D (another ten times) reveal a multitude of small islands, supporting the dense-islands
hypothesis; area C of Fig. 9 contains no dots at grid 0.002 (p.78). The density of islands falls rapidly with distance
from the central region. Not used by the project; recorded for completeness.

## 5. Stated limits (READ)

- The potential is chosen "arbitrary function of R and z, not necessarily representing an actual galactic potential" (p.73);
  the authors say it is "probable that the potential (11) is typical of the general case" (p.75), a conjecture. "The ultimate
  answer ... should rest on rigorous mathematical proofs, not on numerical experiments" (p.79).
- The critical energy is potential-dependent and the authors say further potentials with other angular momentum and higher
  energies are needed (p.77); computations with U = (x^2 + y^2 - x^2 y^2)/2 (not shown) and with Ollongren's approximation to
  the Galactic potential (1962) indicate the opposite situation (critical energy above the escape energy, p.77).
- Whether the curves are exactly or only approximately invariant, and whether the ergodic set is truly a single connected
  region, is left open (p.79).
- Everything is a finite-length computation by 1960s machines; curves are those of a sample of initial conditions.

## 6. Printed numbers usable as sourced tests, and my one check

Sourced (the expected side is a printed value or statement, not computed by the project):

1. The potential (eq. 11) and equations of motion (eq. 12); the escape energy E = 1/6 (p.77), the saddles at U = 1/6.
2. Equipotentials of Fig. 2: U = 0.0100, 0.0417, 0.0833, 0.1250, 0.1667 (labelled, p.74). The U = 1/6 contour is the triangle.
3. Energies of the sections: E = 1/12 (0.08333), 0.12500, 1/6 (0.16667) (pp.75-76).
4. At E = 1/12 every orbit lies on a closed curve (the area is completely covered by curves; pp.76-77).
5. Rotation angle alpha = 0.1143 (one revolution as the unit) for the particular orbit of Fig. 3 at E = 0.08333 (p.75).
   The initial conditions of that orbit are not printed (READ absence; the points are numbered 1 to 28 on Fig. 3 within
   y in about 0 to 0.5 and y-dot in about -0.2 to 0.2, read off the figure), so alpha cannot be tested directly.
6. Critical energy about 0.11 (p.77). Energy of the "abrupt" appearance of a sea is between 0.0833 and 0.125 by the figures.
7. Fig. 7 relative-area plateau 1.0 up to about 0.11 and a rapid fall to about 0.1 at E = 1/6 (p.76, p.77; the point
   values are read off the figure, section 4.3).
8. The criterion mu = sum_(i<=25) distance^2, delta = 1e-7, mu_c about 1e-4 (p.77).
9. The map (14) with a = 1.6: closed curves around the stable invariant point X = Y = 0, a chain of six islands and an outer
   ergodic region (Fig. 8, p.78).

Not printed, so no golden: any initial condition of any orbit; any numerical y-value of an invariant point; any measured
fraction as a number; the integration tolerances.

COMPUTED check (throwaway script in the session scratch area, not committed): at E = 0.118 with the printed potential,
section x = 0, x-dot > 0, a Newton solve for a fixed point of the one-crossing map near (y, y-dot) = (0.305, 0) gives
(0.29546, 0) with a one-crossing map determinant 1.0000 (area-preserving, as the paper states), trace -0.608 (so the orbit is
elliptic, stable). The section's y range at E = 0.118 is -0.4284 to +0.6425. Two cross-checks against the other digest: (a)
the printed C&S label "q2 about 0.305 corresponds to the stable 1-periodic orbit" is close but not identical to my fixed
point 0.2955 (a 3 percent difference; C&S say "about" and give a fixed-point scan only to two or three digits; my number is
for the exact one-crossing map). I record the discrepancy and do not resolve it. (b) The five C&S orbit labels were
tested in the C&S digest's check with this potential, and J matched the labels; that check is now sourced to this paper's
eq. 11 rather than an inferred form.

## 7. Techniques applicable to the project's problems

### `#924`: a sourced regular/chaotic classification test for an FLI or MEGNO module

The expected side must come from the literature, not from a project run. The two printed sources combine as follows.

Test model: eq. 11 and 12 of this paper (sourced). Section: x = 0, x-dot > 0 (this paper, eq. 9; the same section as C&S
Fig. 1a).

Test A (this paper only, qualitative and strong, E = 1/12 = 0.08333): any initial condition on the section lies on an
invariant curve, so a correct indicator must classify EVERY sampled initial condition as regular (FLI-type indicators:
no early threshold crossing; MEGNO: Jbar about 2 within its band). Rule: no chaotic classification anywhere in the allowed
area at this energy (p.76-77, Fig. 4). The paper says "completely covered with curves"; the small chaotic layers at
unstable periodic orbits are not resolved at the paper's resolution and a real indicator may flag a thin layer near the
three intersection points. Use a tolerance on the flagged FRACTION (the paper implies it is at most a few percent; this
is my reading, INFERRED), not zero.

Test B (this paper only, E = 1/6 = 0.16667): the ergodic trajectory covers almost the whole area (Fig. 6, p.76); a correct
indicator classifies the vast majority (the paper's Fig. 7 reading about 0.1 regular) as chaotic.

Test C (this paper, monotonicity): the regular fraction of the section against energy has a plateau at 1.0 up to about
0.11, then falls monotonically to about 0.1 at E = 1/6 (Fig. 7). A module that gets the plateau, the knee near 0.11 and a
monotone fall passes; a tolerance on the knee location of a few hundredths is appropriate given the figure reading.

Test D (this paper and C&S, E = 0.118): the five labelled orbits of C&S (p.204 there), at h = 0.118, y = q2_0 on the line
y-dot = p2 = 0, x = 0, x-dot > 0: q2_0 = 0.305 (near the stable 1-periodic orbit) regular; 0.5 regular; 0.5085 regular
(near the separatrix of the 5-periodic island); 0.509 chaotic (stochastic layer); 0.6 chaotic (stochastic sea). The
independence argument: C&S's labels are printed from their own computation (J^(2) and the Poincare map), not ours, and the
Henon-Heiles potential is now sourced here, so both the model and the expected classes trace to published sources. The C&S
section energy 0.118 is just above this paper's critical energy of about 0.11 and the area covered at 0.118 reads about 0.9
off Fig. 7 (approximate, section 4.3), so a mix with a regular majority and a chaotic sea is what both papers
imply. Caution: orbit 0.5085 sits next to 0.509 across a layer of width of order 1e-3 in q2; an indicator at finite T may
misclassify one of the two. That is a sensitivity test, so give the pair a separate tolerance, and report the five
verdicts.

Positive-control discipline (the project's rule): run the indicator first on a known regular case (0.305 or any E = 1/12
point) and a known chaotic case (0.6), and judge by the same criterion the candidates are held to, before reading a 0/N from
any sweep.

Independence check on C&S's own labels: the printed J^(2) criterion is "about 2 versus much larger" and my throwaway
check in the C&S digest section 7.3 reproduces the labels with the paper's own formula. Both are for the same potential
form, now confirmed.

### Other tasks

- `#908` (capture-sweep sampling): this paper's mu criterion (distance of two close points after 25 map iterates, mu_c about
  1e-4) is a time-to-divergence indicator, a cousin of FLI, on a Poincare section; it shows the same method as the FLI
  digest on a map. Nothing in it applies directly to a non-autonomous capture sweep.
- `#890`, `#895`, `#916`, `#917`: no technique beyond classification practice; this model is autonomous and smooth, and the
  paper contains nothing about moons, flybys or time-dependent forcing.
- Broader value: the paper is the sourced origin of the "mixed phase space" picture (curves, islands, a sea) that the project's
  language about quasi-cyclers, tori and stability maps borrows. It gives an independent, cheap positive control for ANY
  chaos indicator before it is used in an orbital model.

## 8. Follow-ups (no task numbers registered)

1. In the `#924` plan, replace "the logarithmic potential band" as the only control with a second one: Henon-Heiles at
   E = 1/12 (all regular), E = 1/6 (nearly all chaotic) and the five C&S labels at 0.118, all with expected sides from the two
   papers. A test of the fixed-point and chaos fractions with documented tolerances; check what the project asserts before
   writing it.
2. Resolve the 0.305 versus 0.2955 difference: check whether C&S's "stable 1-periodic orbit" refers to a different point
   (for example the center of the regular island measured in q2 on their section) or is a loose label; read their Fig. 1a
   again at higher resolution.
3. Digitise Fig. 7 (7 to 9 dots) at higher resolution if a sourced numeric area fraction is wanted; there is no printed
   table, so the check is a monotone-knee check unless this is done.
4. The mapping (14) is a cheap extra control (a = 1.6, a chain of six islands around X = Y = 0, an outer ergodic region); the
   standard-map control in the FLI digest is the better-sourced one. Optional.
5. Acquire nothing further for this task. Optional context: Contopoulos 1958/1963 and Ollongren 1962 are cited here for the
   same finding in other potentials (p.76).

## 9. Transcription flags

- p.77: "from about 10^-12 to 10^+1" for the range of mu: the exponents are small print in the scan; read as printed, doubtful.
- p.77: mu_c approx 10^-4: legible.
- Fig. 7 dot values: read off the figure; no labels.
- Fig. 3 and Fig. 4 axis ranges: read off the figure.
- p.75: the energy-drift figure "(< |0.00003| for 150 orbits)": legible but the precise meaning is my reading.
