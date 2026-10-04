# Digest: Hadjidemetriou and Christides 1975, "Families of periodic orbits in the planar three-body problem"

Date: 2026-10-05 (Sydney; file name carries the 2026-10-04 series date). Reading, reasoning and independent re-computations; no project code was changed.

Source: J. D. Hadjidemetriou and Th. Christides, "Families of periodic orbits in the planar three-body problem", Celestial
Mechanics 12:175-187 (1975), DOI 10.1007/BF01230210, received 18 March 1974 (University of Thessaloniki). Filed in the
private paper corpus as
`hadjidemetriou-christides-1975-families-periodic-orbits-planar-three-body-problem-celest-mech-12-175-doi-10.1007-BF01230210.pdf`.
A 13-page scan (PDF page n is journal page 174 + n). Every page was read as an image at 130 dpi, and Table I (p.181) again at
220 dpi; no digit was illegible. This is the mass-continuation paper that the 1975b digest said was missing.

Evidence tags: READ (p.N) is read at printed page N. COMPUTED is my own arithmetic or integration (scratch scripts, not
committed; recipe in section 6). INFERRED is my reasoning, not stated by the author.

Companions: `docs/notes/2026-10-04-digest-hadjidemetriou-1975-restricted-to-general-continuation.md` (the existence theorem,
Celest. Mech. 12:155, cited here as "Hadjidemetriou 1974 / 1975, this issue p.155") and
`docs/notes/2026-10-04-digest-hadjidemetriou-1975b-stability-periodic-orbits-three-body.md` (stability, Table I of 26
equal-mass members).

## 0. What the paper is

READ (abstract p.175): "A periodic orbit of the restricted circular three-body problem, selected arbitrarily, is used to
generate a family of periodic motions in the general three-body problem in a rotating frame of reference, by varying the
mass m3 of the third body. This family is continued numerically up to a maximum value of the mass of the originally small
body, which corresponds to a mass ratio m1:m2:m3 ~ 5:5:3. From that point on the family continues for decreasing masses m3
until this mass becomes again equal to zero. It turns out that this final orbit of the family is a periodic orbit of the
elliptic restricted three body problem."

One orbit, one path, one table: 19 points of one characteristic curve. The numbers are real and checkable. The scientific
claims are (i) the continuation works up to m3 = 0.2245 (fixed ratio m1/m2 = 1, total mass 1) and then turns back, (ii) the
return branch ends at m3 = 0 in an elliptic-restricted-problem periodic orbit with primaries' eccentricity 0.292, (iii)
families for fixed masses follow by repeating this for all members of a restricted family.

Naming (same as the 1975b paper, opposite to the 1975 existence paper): P1 and P2 are the primaries (masses m1, m2), P3 is
the third body that is massless in the restricted limit. mu = m2/(m1 + m2) (eq. 2); P1 at x1 = mu r, P2 at x2 = -(1 - mu) r
in the rotating frame G1xy (eq. 1), r the P1-P2 distance. Here m1 = m2, so mu = 1/2 throughout.

## 1. Equations as printed (READ pp.176-178)

Frame (Fig. 1, p.176): G is the centre of mass of the system (fixed in the inertial frame), G1 the centre of mass of P1 and
P2; the rotating frame G1xy has its origin at G1 and its x axis always through P1 and P2, at angle theta to a fixed inertial
direction. Generalised coordinates (x3, y3) of P3 in the frame, r, and theta (four degrees of freedom). READ eq. (1):
x1 = mu r, x2 = -(1 - mu) r, G1G = (m3/m) r-vector (r-vector the radius vector of P3), with m = m1 + m2 + m3 (eq. 2).
Lagrangian (eq. 3), identical to eq. (1) of the 1975b paper after renaming:

    L = (1/2)(m1 + m2) { mu(1 - mu)(r'^2 + r^2 th'^2) + (m3/m)[ x3'^2 + y3'^2 + th'^2 (x3^2 + y3^2) + 2 th'(x3 y3' - x3' y3) ] } - V,
    V = -G m1 m3/r13 - G m2 m3/r23 - G m1 m2/r12   (eq. 4).

Both dL/dt = 0 and dL/dtheta = 0, so the energy and the angular momentum are integrals. Equations of motion for G = 1
(eq. 5):

    r'' - r th'^2 = -(m1 + m2)/r^2 - m3 A,
    x3'' - 2 th' y3' - th'^2 x3 - th'' y3 = m B,
    y3'' + 2 th' x3' - th'^2 y3 + th'' x3 = m C,
    th'' + 2 r' th'/r = -(m3 / (mu(1 - mu) r^2)) [ C x3 - B y3 ],

with A = -(x3 - mu r)/r13^3 + (x3 + (1 - mu) r)/r23^3 (eq. 6), B = -(1 - mu)(x3 - mu r)/r13^3 - mu (x3 + (1 - mu) r)/r23^3
(eq. 7), C = -(1 - mu) y3/r13^3 - mu y3/r23^3 (eq. 8). The authors integrate the whole system (5), with theta kept as a
variable, rather than the angular-momentum-reduced Routhian (p.178). Transcription note: I have copied eq. (5) as printed;
I have not re-derived the sign of the m3 A and m B terms, and my verification below does not use them (it uses an inertial
Newtonian integration), so a typesetting slip in (5) to (8), if there were one, would not show up in my numbers.

## 2. The method: continuation in m3, step by step (READ pp.178-181)

### 2.1 Symmetry and periodicity conditions (eq. 9-14)

The transformation r -> r, x3 -> x3, y3 -> -y3, theta -> -theta, t -> -t leaves (5) unchanged (eq. 9). Initial conditions
(eq. 10), at t = 0: r = r0, r' = 0, x3 = x30, y3 = 0, x3' = 0, y3' = y3'0, theta' = theta'0. At t = T the orbit again crosses
y3 = 0 with r' = x3' = 0 (eq. 11); by (9) the motion in the rotating frame is then periodic with period 2T, with P3 starting
perpendicularly from the x axis at the instant P1 and P2 are at rest on it (so the primaries' separation is at an extremum
whenever P3 crosses the axis). Normalisation (p.178): G = 1, m1 + m2 + m3 = 1, theta'0 = 1. Hence the angular momentum
varies along the family. Start: m3 << 1 and (eq. 12) r = r0 (about 1), x30, y3'0 near a restricted-problem periodic orbit.
Periodicity conditions (eq. 13): r'(r0, x30, y3'0, T) = 0, x3'(r0, x30, y3'0, T) = 0, with T the first (in general the nth)
crossing of the axis, defined by y3(r0, x30, y3'0, T) = 0 (eq. 14).

### 2.2 What is held and what is varied (READ p.179, p.185-186)

For fixed masses there are three free parameters (r0, x30, y3'0) and two conditions, so one parameter family. Two ways of
proceeding (p.179): fix m3 and vary x30 (a family of periodic orbits for fixed masses m1, m2, m3), or fix x30 and vary m3
(a family for fixed ratio m1/m2 and fixed x30, "the continuation of the periodic orbit of the restricted problem which
corresponds to those values of m1/m2 and x30 and a certain y3'0, or the Jacobi constant"). The paper takes the second:
x30 and m1/m2 fixed, m3 varied, with m1 + m2 + m3 = 1 held so that m1 = m2 = (1 - m3)/2.

In the discussion (pp.185-186): in general x30 = f(m3) is a free choice; each f gives a different "characteristic curve", a
different (m3)max and a different final orbit. "What in fact happens" is that with m3 as an extra parameter the family is
two-parametric: the initial conditions (r0, x30, y3'0) lie on a two-dimensional "characteristic surface" in the
three-dimensional space, the lines m3 = constant being the monoparametric families for fixed masses; a given continuation
is a path on this surface and "the particular shape of a characteristic curve does not have any real significance; the
complete picture can only be revealed from the study of the topology of the whole characteristic surface".

### 2.3 The corrector (eq. 15, READ p.179)

Newton iteration on the two unknowns (r0, y3'0) with x30 and m3 given: integrate (5) to the first y3 = 0 at t = T. If
Delta r0, Delta y3'0 are the corrections that give the periodic motion, whose period is 2(T + Delta T), then (eq. 15)

    A1 Dr0 + A2 Dy3'0 = r'(r0 + Dr0, x30, y3'0 + Dy3'0, T + DT) - r'(r0, x30, y3'0, T)
    A3 Dr0 + A4 Dy3'0 = x3'(r0 + Dr0, x30, y3'0 + Dy3'0, T + DT) - x3'(r0, x30, y3'0, T),

with r'(T + DT) = 0 and x3'(T + DT) = 0 imposed. The 2 x 2 matrix coefficients are the derivatives of the solution at t = T
with respect to the initial values corrected for the moving crossing time:

    A1 = dr'/dr0 - [r''(T)/y3'(T)] dy3/dr0,   A2 = dr'/dy3'0 - [r''(T)/y3'(T)] dy3/dy3'0,
    A3 = dx3'/dr0 - [x3''(T)/y3'(T)] dy3/dr0, A4 = dx3'/dy3'0 - [x3''(T)/y3'(T)] dy3/dy3'0

(the printed primes and dots are time derivatives; "r-double-dot(T)" is my reading of the printed second derivatives, which
the page shows as dots). These come from two extra integrations (each stopped at y3 = 0) with increments of 10^-6 in r0 and
then y3'0. Practical economies, READ p.179: once |r'| + |x3'| < 10^-4 the same A_i are reused in all later iterations;
"in practice three to five iterations were enough to minimize r' and x3' up to ten decimal places"; the scheme is
"similar in principle" to Szebehely and Peters (1967). Predictor, READ p.181: in part EF of the characteristic curves the
approximate initial conditions have to be right to better than 10^-3, "otherwise a periodic motion belonging to a different
family was discovered": a warning that the corrector can hop families.

Integrator (READ p.180): series expansion (Deprit and Price 1965), 11 decimal places retained; "retained its accuracy even when
two of the bodies approached to a distance of the order of 10^-5, so it was not necessary ... to regularize". (For the
project's `#928`, this is a counter-example of sorts: a high-order Taylor integrator carried a close pass at 1e-5 without
regularisation in 1974.)

## 3. Results as printed (READ pp.181-185)

Start (p.181): the restricted circular orbit with x30 = 0.181, y3'0 = -0.947 110 34 for m1/m2 = 1, from Szebehely (1967),
"belongs to the family (g)". Fix x30 = 0.181, keep m1 + m2 + m3 = 1 and m1/m2 = 1, theta'0 = 1, G = 1. "It was found that this
family of periodic motions is continued up to a maximum value of the mass m3 equal ... to m3 = 0.2245, which corresponds to
a ratio of the masses equal to about m1:m2:m3 = 5:5:3." (COMPUTED: m1 = m2 = (1 - 0.2245)/2 = 0.38775; 5:5:3 would be
m3 = 3/13 = 0.2308, so the ratio is approximate, as the text says.) Beyond that point the family turns back and continues
for decreasing m3 to m3 = 0 (points A, B, C, D, E, F of Figs 2, 3, 5). The final orbit is a periodic orbit of the elliptic
restricted problem with eccentricity e = 0.292 (abstract, p.176, p.185); it is "similar to the collision orbits of the
elliptic restricted three body problem studied in detail by Broucke (1969)", but the family has no collision orbit; at F
the minimum P1-P3 distance is "of the order of 0.001" (p.185).

### Table I (READ p.181), the 19 members of the characteristic curve

Columns: m3, x10 (the coordinate of P1 at t = 0, so r0 = 2 x10 for m1 = m2), y3'0, x1(T), x3(T) (positions at the half
period, the perpendicular crossing), the period 2T, and theta/2pi, the rotation of the frame over one period divided by
2 pi. x30 = 0.181 and theta'0 = 1 for every row. The rows are given in order along the curve; the m3 values repeat on the
way up (rows 2 to 12) and the way down (rows 13 to 19). The row numbers are mine.

| row | m3 | x10 | y3'0 | x1(T) | x3(T) | period 2T | theta/2pi | point |
|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 0.50000000 | -0.94711034 | 0.50000000 | 0.72101839 | 1.8484 | 0.29419126 | A |
| 2 | 0.0010 | 0.50092191 | -0.94511399 | 0.50033711 | 0.72165922 | 1.8560 | 0.29490893 | |
| 3 | 0.0100 | 0.50909746 | -0.92781324 | 0.50315431 | 0.72712559 | 1.9240 | 0.30138460 | |
| 4 | 0.0500 | 0.54307785 | -0.86342556 | 0.51142589 | 0.74554248 | 2.2244 | 0.33111567 | |
| 5 | 0.1000 | 0.58106608 | -0.80535147 | 0.51335472 | 0.75619699 | 2.6078 | 0.37251151 | |
| 6 | 0.1200 | 0.59504112 | -0.78781819 | 0.51168416 | 0.75629375 | 2.7668 | 0.39105923 | B |
| 7 | 0.1500 | 0.61470336 | -0.76738429 | 0.50631794 | 0.75064525 | 3.0140 | 0.42200972 | |
| 8 | 0.1800 | 0.63254710 | -0.75535104 | 0.49643753 | 0.73460805 | 3.2760 | 0.45888914 | |
| 9 | 0.2000 | 0.64291315 | -0.75497081 | 0.48543791 | 0.71289250 | 3.4654 | 0.49000811 | |
| 10 | 0.2200 | 0.64988645 | -0.77337080 | 0.46325010 | 0.66276240 | 3.6880 | 0.53839883 | |
| 11 | 0.2240 | 0.64910244 | -0.78970281 | 0.45126636 | 0.63398810 | 3.7526 | 0.55965075 | C |
| 12 | 0.2245 | 0.64854088 | -0.79446754 | 0.44814903 | 0.62643112 | 3.7646 | 0.56477324 | max m3 |
| 13 | 0.2240 | 0.64284614 | -0.82658400 | 0.42939154 | 0.58102449 | 3.8070 | 0.59285687 | |
| 14 | 0.2200 | 0.63554886 | -0.85813511 | 0.41310999 | 0.54248007 | 3.8116 | 0.61439245 | D |
| 15 | 0.1800 | 0.58857781 | -1.02869781 | 0.34032721 | 0.39120985 | 3.6656 | 0.70005821 | |
| 16 | 0.1200 | 0.53513036 | -1.21749666 | 0.28459823 | 0.30014806 | 3.5334 | 0.78750832 | E |
| 17 | 0.0500 | 0.48265003 | -1.41273701 | 0.25296896 | 0.25639133 | 3.5264 | 0.89757983 | |
| 18 | 0.0100 | 0.45321081 | -1.52856693 | 0.24485183 | 0.24601948 | 3.5784 | 0.97702084 | |
| 19 | 0 | 0.44559306 | -1.55950866 | 0.24406810 | 0.24495936 | 3.5986 | 1.00000000 | F |

(The PDF prints digits in groups of two or three; I have joined them. The "point" column is my assignment from the Fig. 4
caption, p.183: (a) m3 = 0, point A; (b) m3 = 0.12, point B; (c) m3 = 0.2240, point C; (d) m3 = 0.22, point D; (e) m3 = 0.12,
point E; (f) m3 = 0, point F. Ambiguity: the table does not label the points, and each of 0.22, 0.12 appears twice, so
"B = ascending 0.12 and E = descending 0.12" and "C = ascending 0.2240, D = descending 0.22" are my reading of Figs 2 and 3
and of the caption, not printed in the table. In particular the caption's "(c) m3 = 0.2240 point C" and "(d) m3 = 0.22
point D" put C before the maximum 0.2245 and D after it.)

Characteristic curve shapes (Figs 2 and 3, pp.180): x10 against m3 rises from 0.5 to a maximum near 0.65 at C and falls to
0.4456 at F; y3'0 against m3 rises from -0.947 to a maximum near -0.755 around m3 = 0.2, then falls to -1.5595 at F.
Fig. 5 (p.184): theta/2pi against m3 runs from 0.2942 at A up to about 0.565 at the maximum and back through 2/3 to 1 at F.

Orbit shapes (Fig. 4a-f, pp.182-183), READ: "all the periodic orbits are simple periodic orbits around the primary P1. Up to
the orbit indicated by C ... they are simple orbits and from there on a loop is formed before the orbit crosses the x axis."
Fig. 4a is the circular restricted orbit (a near-circle around P1, between x about 0.2 and 0.7), 4f the elliptic one
(two lobes with a figure-of-eight shape, almost a collision pass with P1).

Inertial-frame periodicity (READ p.185): "we have simple periodic motion of the three bodies in the inertial system for
values of theta/2pi equal to the simple ratios 1/3, 1/2, 2/3, 1" (and infinitely many with other rational values; isolated
for fixed masses). Fig. 6 (p.184) is an example for theta/2pi = 1/2: the caption prints m1 = 0.3975, m2 = 0.3975,
m3 = 0.2050 (m1:m2:m3 ~ 2:2:1), x10 = 0.644 35, y3'0 = -0.758 77, all three orbits symmetric about both axes, the numbers
1 to 4 marking times t = 0, 0.5, 1.0, 1.78. See section 5 for what this does and does not reproduce.

Terminal elliptic orbit (READ p.185): "As m3 tends to zero the periodic motion of the three bodies reduces to a periodic orbit
of the elliptic restricted three-body problem. This is easily seen as r0 not equal to 1 at m3 = 0 while theta'0 = 1. The eccentricity
of the two primaries is equal to e = 0.292." Why it is elliptic: with theta'0 = 1 and G(m1 + m2) = 1 the primaries are on a
circle only if r0 = 1; at row 19, r0 = 2 x10 = 0.89118612.

Discussion claims (pp.185-187): the family is not only an isolated point; every member of a restricted family can be
continued, so a family of restricted orbits continues to a family of the general problem for each m3; the continuation could
instead lead to m3 = 1 (two massless particles in Kepler orbits around P3), or back to the circular restricted problem or to
the starting orbit, depending on f(m3); an orbit of the elliptic restricted problem can be continued to the general case in
the same way; starting from a general-problem orbit one may be unable to extend it beyond a minimum value of a mass, "which
would imply that this periodic orbit of the general problem is of a different nature than the periodic orbits of the
restricted problem since no connection can be established between them by a continuous variation of the parameters"; the
Szebehely and Standish orbits, isolated in the inertial frame, might be continued to m3 = 0 in the rotating frame (stated
as a possibility, p.186); Henon (1974) "in agreement, in principle" (p.187). Normalisation remark (p.187): theta'0 = 1 makes the
angular momentum vary along the family; the alternative is angular momentum equal to unity, together with G = 1 and
m1 + m2 + m3 = 1; both fail for rectilinear motion.

The paper prints no energy, angular momentum, Jacobi constant, stability index or eigenvalue anywhere. The only
stability-related statement is in the Introduction and abstract of the companion 1975b paper.

## 4. Independent re-computation (COMPUTED)

An inertial Newtonian three-body integrator (not the paper's equations (5) to (8); DOP853, rtol 1e-13, atol 1e-14), G = 1,
masses (m1, m2, m3) = ((1 - m3)/2, (1 - m3)/2, m3), set up at t = 0 as in the 1975b recipe: in the frame rotating with the P1-P2
line, P1 = (x10, 0), P2 = (-x10, 0), P3 = (0.181, 0), rotating-frame velocities (0, 0), (0, 0), (0, y3'0), inertial velocity
= rotating velocity + theta'0 (= 1) times the 90-degree rotation of the position; centre-of-mass subtraction. T is the first
return of the rotating-frame y of P3 to zero; theta/2pi is the unwrapped rotation of the P1-P2 line over 2T divided by 2 pi.

All 19 rows of Table I reproduce from the printed initial conditions:
- 2T equals the printed period to within 1e-4 (the printed precision) in 18 rows. Row 15 (m3 = 0.18, descending) is the
  exception: printed 3.6656, computed 3.66967 (difference 4.1e-3). Everything else about row 15 matches (x1(T), x3(T) and
  theta/2pi, the last computed over my 2T), so the printed 3.6656 is very probably a misprint for 3.6697 (INFERRED, a
  digit slip); the neighbouring periods 3.8116 and 3.5334 and a smooth curve agree with 3.6697. Treat 3.6656 as a known
  erratum and test row 15's period at 3.6697, or skip that column for row 15;
- x1(T) and x3(T) agree with the printed values (eight decimals) to within 2.1e-8 in every row;
- x1'(T) and x3'(T) are zero to 3.2e-6 or better in every row (3.2e-6 at row 19, 1e-7 or better in rows 1 to 16);
- theta/2pi over 2T agrees with the printed column to within 7.1e-8 in every row (rows 1 and 19 to the printed digits);
- the full state at 2T returns to the initial rotating-frame state to 1e-9 (rows 1 to 14), 4e-7 (row 15) and up to 4e-6 (rows
  16 to 19, where P1 and P3 pass within 3e-3 to 9e-4 of each other and the printed 8-digit initial conditions limit the closure).
This is a clean positive control: the table is internally consistent, the conventions (G = 1, total mass 1, theta'0 = 1, x30 =
0.181, P1 at +x10) are confirmed, and the 1975b-style inertial integrator is validated against a second, independent source.

Quantities not printed in the paper (COMPUTED, for use as outputs of a test code, not as sourced expected values):

| row | m3 | E (total energy about the centre of mass, G = 1) | p (angular momentum) |
|---|---|---|---|
| 1 | 0 | -0.125000 | 0.250000 |
| 5 | 0.1 | -0.176301 | 0.293703 |
| 9 | 0.2 | -0.203034 | 0.314048 |
| 12 | 0.2245 | -0.211188 | 0.306848 |
| 14 | 0.22 | -0.216625 | 0.294028 |
| 16 | 0.12 | -0.220993 | 0.232190 |
| 19 | 0 | -0.181248 | 0.198553 |

Two of these have closed forms and so give independent analytic checks of the setup (COMPUTED): row 1 is two equal masses 0.5
on a circle of radius 1 about their centre of mass with unit angular velocity and a massless P3: E = -m1 m2/(2 r) = -0.125, p =
0.25 x 1 x 1 = 0.25, as computed. Row 19: the primaries are a Kepler pair, with h = r0^2 theta'0 = 0.79416, semi-latus rectum
h^2/(m1 + m2) = 0.63069, r0 = 0.89118612 at an extremum so r0 = p/(1 - e) = apoapsis, giving e = 0.292209 (printed 0.292),
semi-major axis a = 0.689661, period 2 pi a^(3/2) = 3.598600 (printed 2T = 3.5986), E = -m1 m2/(2a) = -0.18125, p = 0.25 x
sqrt(0.63069) = 0.19855, all as computed numerically. The minimum P1-P3 distance over the period at row 19 is 8.9e-4
(printed "of the order of 0.001") and at rows 17 and 18 3.4e-3 and 1.2e-3. The primaries' separation at row 19 ranges 0.4881 to 0.8912
(= a(1 -/+ e), as expected).

## 5. Test-ready numbers and what does not reproduce

Sourced test material (every value printed, page and table):
1. Table I (p.181), all 19 rows: each row is a pair (m3, x10, y3'0) from which the period 2T, x1(T), x3(T) and theta/2pi are
   expected outputs (with x30 = 0.181, theta'0 = 1, G = 1, m1 = m2 = (1 - m3)/2). Tolerances that pass with my integrator: 2T
   to 1e-4 (the printed precision) in every row except row 15, whose printed 3.6656 is a probable misprint (computed
   3.6697); x1(T), x3(T) to 5e-8; theta/2pi to 1e-7; closure at 2T not better than 1e-5 for rows 16 to 19.
2. Szebehely restricted orbit (p.181): x30 = 0.181, y3'0 = -0.947 110 34, m1/m2 = 1 is a periodic orbit of the planar CR3BP
   with mu = 1/2; row 1 is its half-period state. Check against the project's own code (done, section 7).
3. m3 maximum: 0.2245 (p.181), at which x10 = 0.64854088, y3'0 = -0.79446754 (row 12); the maximum m3 itself is read from
   the table's turn, not separately computed.
   Erratum candidate: row 15 period (printed 3.6656, computed 3.6697).
4. Terminal orbit: e = 0.292 (pp.176, 185), and the printed period 3.5986 and theta/2pi = 1 reproduce analytically as above.
5. 5:5:3 ratio (abstract, p.181): m3 = 0.2245 gives m1:m2:m3 = 0.38775:0.38775:0.2245 = 5:5:2.89.

Does not reproduce, or is not usable as a number:
- Fig. 6 (p.184), x10 = 0.64435, y3'0 = -0.75877, m3 = 0.2050. I solved for the periodic orbit at x30 = 0.181 and m3 = 0.2050:
  the corrector gives x10 = 0.645171, y3'0 = -0.756661 (T = 1.7581, theta/2pi = 0.49942, closure 2e-10). Taking theta/2pi = 1/2
  exactly as the unknown too: m3 = 0.205293, x10 = 0.645297, y3'0 = -0.756795, T = 1.7596. The printed Fig. 6 values differ from
  these by about 1e-3 in x10 and 2e-3 in y3'0, and the caption's last time label 1.78 is not T. Run directly from the printed
  5-digit values at m3 = 0.2050 the orbit closes only to 5e-3 and x1'(T) = -2.9e-3. INFERRED: the figure's initial conditions were read
  from the characteristic curves of Figs 2 and 3 (linear interpolation of Table I between m3 = 0.20 and 0.22 gives x10 = 0.6448
  at 0.2053), not computed to the printed digits; treat the Fig. 6 numbers as illustrative, with a tolerance of at least 3e-3, and
  use the theta/2pi = 1/2 property (m3 = 0.2053) rather than the caption's digits.
- E, p, Jacobi constant, stability: not printed.

## 6. Recipe for the re-computation

As in section 4 and the recipe of the 1975b digest (section 7 there). For this paper's rows the same integrator needs no
re-solve of theta'0 (theta'0 = 1 is fixed and x30 = 0.181 is fixed; the periodicity is a property of the printed (x10, y3'0) pair),
so no root finder is needed for the table check, only an event-located first return of y3 to zero, and the unwrapped
rotation angle of the P1-P2 line. The Fig. 6 re-solve used a Newton solve in (x10, y3'0) (and m3 for the 1/2-rotation case)
on x1'(T) = x3'(T) = 0. Time for the table: seconds per row.

## 7. Does it connect to the 1975b Table I equal-mass family?

Short answer: not from anything printed in either paper. Facts (READ unless tagged):
- Same convention: both use G = 1, total mass 1, theta'0 = 1, the frame G1xy with the x axis through P1 and P2 (the binary), P3
  the body that is massless in the restricted limit. So E, p, T and x, y, y' are directly comparable in form; my E and p
  (total energy and angular momentum about the centre of mass) reproduce the 1975b table, and I use the same here.
- The 1975b paper (p.267, section 7) says one member of its equal-mass family "has been obtained by continuing numerically a
  periodic orbit of the restricted circular three-body problem by varying the mass m3" and refers to this paper for the
  procedure. But this paper's continuation fixes x30 = 0.181 and stops (turns back) at m3 = 0.2245, where m1 = m2 = 0.38775;
  the equal-mass case needs m3 = 1/3. And the 1975b Table I has x30 from 0.310 to 1.022, whereas this one has x30 = 0.181.
  So the single path printed here does not reach the 1975b family.
- The paper itself says (p.186) that a different choice f(m3) = x30(m3) gives different (m3)max and different end orbits
  (INFERRED consequence): the 1975b member must have come from a different path on the characteristic surface (a different
  f(m3), a different start orbit, or a continuation at fixed m3 in x30, p.179), none of which is printed.
- Orbit class differs in the printed figures (INFERRED from the pictures): here every member is a simple loop around P1 (Fig. 4)
  in the rotating frame, with the primaries oscillating by at most 0.49 to 0.90 in separation; the 1975b description is a binary
  with the third body circling it at a distance that varies along the family. I did not compare the two classes beyond that.
- Numbers do not overlap: this family has E from -0.125 to -0.225 and p from 0.199 to 0.315 at total mass 1; 1975b Table I
  has E from -0.811 to -0.123 and p from 0.341 to 0.363. A shared member would need equal masses, which this path never
  reaches, so there is no common row to compare.
Conclusion: the two papers are about the same method and normalisation, and 1975b's family is built by the method described
here, but the actual 1975b starting orbit and continuation path are not printed in either paper. Hadjidemetriou's 1974/75
existence digest therefore stays the existence statement, this paper the worked example for ratio m1/m2 = 1 only to m3 = 0.2245.

## 8. Reconciliation with the project's code

- No general three-body code in `src/` (verified again: no `hadjidemetriou`, `isoenergetic`, `routh` hits; the `nbody/`
  package is a restricted REBOUND wrapper; `core/*` are restricted models). The paper's setting (all three bodies massive,
  primaries recoiling, frame tied to the P1-P2 line) has no project model. That is the `#931` gap, and this paper
  adds a second, independent source for the control integrator.
- One existing module reproduces part of the table directly: `core/cr3bp.py::cr3bp_eom` at mu = 0.5. COMPUTED (DOP853, rtol 1e-13)
  from the printed row-1 state (x, y, vx, vy) = (0.181, 0, 0, -0.94711034) in the project's rotating frame (the frame of
  `cr3bp_eom`, origin at the barycentre, which equals G1 for equal primaries; unit separation and unit angular velocity
  agree with theta'0 = 1, G(m1 + m2) = 1): the first y = 0 return is at T = 0.9242291 (printed 2T/2 = 0.9242), with
  x = 0.721018392 (printed x3(T) = 0.72101839) and vx = -1.3e-9, and the orbit returns to the start at 2T. Jacobi constant of this
  orbit, `cr3bp.jacobi_constant`: 3.738968 (COMPUTED, not printed). So the project's restricted model reproduces the printed
  starting orbit to eight digits: a sourced positive control for `core.cr3bp` with mu = 1/2, new to the repository (grep for
  0.94711034 in `src` and `tests`: no hits).
- Not covered by any module: rows 2 to 19 (finite m3, recoil).

## 9. Techniques applicable to the project's problems

Entries read from `data/OUTSTANDING.md` first: `#931` (a) to (c), `#890` and `#895` (Titania-Oberon, near lines 591 and 1406),
`#899` (near 1247).

### 9.1 `#931` (a), the equal-mass test integrator

1. Add this paper's Table I as a second sourced control set to the same test-only module, with m1 = m2 = (1 - m3)/2 and m3 as
   a parameter (the equal-mass 1/3 case of 1975b is the special case m1 = m2 = m3). Test: for every row, the first return of
   the rotating-frame y of P3 to zero occurs at half the printed period, x1(T) and x3(T) match the printed values, and the
   unwrapped rotation over the period matches the printed theta/2pi, tolerances as in section 5. Control: the m3 = 0 rows
   (1 and 19) have the closed forms of section 4 (E, p, e = 0.2922, period 3.5986), so the same code is checked against
   analytics at both ends and against print in between.
2. Positive control for the corrector, not just the integrator: a Newton corrector on (r0 or x10, y3'0) at fixed x30 = 0.181
   and fixed m3, started from a perturbed row, must reproduce the printed (x10, y3'0) (the paper's eq. 15 corrector, or any
   equivalent one). Starting 1e-3 off should converge in about 3 to 5 iterations (the printed count).
3. Mass-continuation test: step m3 from row 1 upward with a predictor (tangent or secant); the corrector must return the
   printed (x10, y3'0) at m3 = 0.001, 0.01, 0.05, 0.1, 0.12, 0.15, 0.18, 0.2, 0.22, 0.224, and must detect that no solution
   continues past m3 = 0.2245 on this branch (the fold, where the path turns); passing through the fold needs
   pseudo-arclength in (m3, x10, y3'0), not parameter continuation. That is a sharper check of the project's continuation
   machinery than any restricted-problem fold it has, and it is the same structure as a mu-continuation (`search/mu_continuation.py`)
   whose fold behaviour is currently only checked on restricted problems (not read here; I did not open that module).
4. Warn on branch-hopping: the paper says (p.181) that in EF the predictor must be within 1e-3 or a different family is found;
   a continuation test should assert that the corrected solution lies within a stated distance of the secant prediction.

### 9.2 `#890` and `#895` (Titania-Oberon; methods only)

The parameter the paper continues in is a body's mass with the others fixed in ratio, from zero to finite; `#890` raised the
moons' masses with a massless spacecraft (companion digest, section 5.1), so this paper is the other direction. What carries
over as a method:
1. The continuation target is a periodic orbit in the frame tied to two bodies (here the P1-P2 line), with symmetric
   perpendicular crossings; for Titania and Oberon the analogue is the frame of the two moons and a spacecraft orbit crossing
   the axis, but the moons' mutual force and Uranus's mass would be part of a four-body generalisation that this paper does not
   treat (three bodies only).
2. The characteristic-surface view (p.186): the solution set is a surface in (r0, x30, y3'0) with m3 as a coordinate; a
   continuation is a path on it; paths differ and the shape of one path "does not have any real significance". For `#890`'s
   continuation in moon-mass scale, record the whole (scale, x, vy, T) data, not a single curve, and test whether another
   path (another f) reaches the same endpoint at scale 1; that tests whether the reported orbit is on a connected component
   or a path artefact. Control: this paper's own claim that a different f(m3) changes (m3)max, so the test must be expected to
   differ at the fold unless the path is reparametrised.
3. A test the paper makes possible: the integrator-and-corrector control of 9.1 is a sourced check of the finite-mass
   recoil dynamics, which neither `#890` (prescribed moons) nor `#895` (kernel positions) contains. If the moons' reaction to
   each other (not to the spacecraft) ever matters at the level of the Oberon sensitivity reported in `#895` (1 m/s moves the flyby
   periapsis by 1,200 to 1,400 km), a finite-mass three-body Titania-Oberon-Uranus orbit with the spacecraft massless is
   a restricted-problem inside a finite-mass one, which is exactly the setting where this paper gives checks: massless limit
   reproduces the restricted orbit (row 1 at eight digits) and the recoil of order m3 grows smoothly (the x10 and 2T columns).
   The moons' masses are tiny compared with Uranus's, far from this paper's range (m3 up to 0.2245 of the total), so no
   number from the table carries over; only the methods and the validated integrator do.
4. Close approaches: the paper integrated a 9e-4 P1-P3 pass (relative to the unit primary separation) without
   regularisation using a high-order series method and 11 digits retained. A Titania or Oberon flyby at 1,000 to 2,000 km
   altitude is a pass at a few 1e-3 of the moon's distance from Uranus (moon orbital radii are of order 4e5 km), the same order
   as that pass; this is mild support for the `#895` choice of an unregularised integrator (DOP853, IAS15), with the caveat that
   the comparison of scales is mine and the paper's pass is relative to the primaries' separation, not to a planet.
   The near-collision seeds of `#899` (1e-5 scale) are beyond what is shown here.

### 9.3 `#899` (continuation of second-species seeds in mass)

The paper is the printed precedent for continuation in the third mass from the restricted orbit to a finite mass, with the
corrector and predictor stated, and it reports two behaviours the Casoliva et al. summary also mentions, in a different
setting: (a) the continuation need not be monotonic in the mass; it reached m3 = 0.2245 and turned back, so a
mass-continuation code that assumes a monotonic parameter will stop at the fold (this is what a fixed-period continuation
that "in most cases ends on lunar impact" can also do; the mechanism differs, a fold here and an impact there); (b) the final
orbit can be a periodic orbit of a different restricted problem (here the elliptic one), "linked" to the circular one through
orbits of the general problem. For `#899`: record at each continuation step the min distance to the primaries (as the paper does
at F, 9e-4) and distinguish termination by fold (m3 turning), by near-collision (distance falling) and by family hop
(predictor miss over 1e-3). Control: this paper's curve, where all three occur or are excluded as stated: a fold at C, a
near-collision at F (8.9e-4), and the warning about hops in EF. The paper does not address second-species orbits or
close passes of the Earth-Moon type, and its smallest P1-P3 distance (9e-4 of the primary separation, of the order 1e-3) is two
orders of magnitude above `#899`'s 1e-5-scale passes, so it tests the continuation logic, not the seeding.

## 10. Recommended follow-ups (no task numbers registered)

1. Add Table I of this paper to the `#931` (a) test module: all 19 rows (section 5 tolerances), the analytic ends (rows 1 and 19:
   E = -1/8, p = 1/4; e = 0.2922, a = 0.68966, period 3.5986), and `core.cr3bp` at mu = 1/2 from row 1 (T = 0.92423, x = 0.721018392,
   return to the start at 2T).
2. A pseudo-arclength continuation test on the (m3, x10, y3'0) curve through the fold at m3 = 0.2245, with the printed rows as
   the expected output, as a control for mass-continuation in the project (including `#899`).
3. A second characteristic curve to resolve the 1975b ambiguity: continue the x30 = 0.181 orbit's neighbours in x30 at fixed m3
   = 1/3 (the other option on p.179) to see whether a path reaches an orbit whose (x10, x30, y3'0) matches a 1975b Table I
   row; the 1975b rows have x30 from 0.31 upward, so a different restricted start orbit is the likelier candidate. This is a
   search, not a reproduction, and may fail.
4. Treat the Fig. 6 digits (x10 = 0.64435, y3'0 = -0.75877) as unreproduced published values: the exact orbit at the caption's mass has
   x10 = 0.645171, y3'0 = -0.756661. Under the project's "never give up reproducing papers" ladder the next rung is whether Fig. 6 uses a
   slightly different x30 or a different (m1, m3) than printed; no evidence either way. Low priority.
5. Obtain Szebehely and Peters 1967 (Astron. J. 72:1187) for the 2 x 2 corrector it is "similar in principle" to, and Henon 1974
   (Celest. Mech. 10:375).

## 11. Summary

- Method: start from a circular restricted orbit (mu = 1/2, x30 = 0.181, y3'0 = -0.94711034, Szebehely 1967 family (g)), keep x30
  and m1/m2 = 1 fixed, set m1 + m2 + m3 = 1, G = 1, theta'0 = 1, raise m3. Unknowns (r0, y3'0); conditions r' = x3' = 0 at the first
  axis crossing T; a 2 x 2 Newton corrector with a moving-crossing-time correction (eq. 15), increments 1e-6, three to five
  iterations to 1e-10; series-expansion integrator, 11 digits, no regularisation needed at 1e-5 close passes.
- Results: continuation reaches m3 = 0.2245 (about 5:5:3), turns, and returns to m3 = 0 as an elliptic-restricted-problem orbit,
  e = 0.292, with a 9e-4 minimum P1-P3 distance. 19 printed rows (Table I, p.181). The family is a path on a two-dimensional
  characteristic surface; other paths give other (m3)max and end orbits.
- Reproduced: all 19 rows (2T to the printed 1e-4 except row 15, whose printed period 3.6656 looks like a misprint for
  3.6697; x1(T), x3(T) to 2e-8; theta/2pi to 7e-8), the analytic ends, and the starting
  orbit by the project's own `core.cr3bp` at mu = 1/2 (to 8 digits). Not reproduced: the Fig. 6 initial conditions (differ by about 1e-3).
- Connection to 1975b Table I: not printed in either paper; the path here never reaches equal masses (stops at m3 = 0.2245).
- Project: no general three-body code; this table is a second sourced control set for `#931` (a), and the continuation through the
  fold is a control for mass-continuation work (`#899`).
