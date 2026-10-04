# Digest: Hadjidemetriou 1975, "The stability of periodic orbits in the three-body problem"

Date: 2026-10-04 (Sydney). Reading, reasoning and small independent re-computations; no project code was changed.

Source: J. D. Hadjidemetriou, "The stability of periodic orbits in the three-body problem", Celestial Mechanics 12:255-276
(1975), DOI 10.1007/BF01228563, received 11 June 1974 (University of Thessaloniki). Filed in the private paper corpus as
`hadjidemetriou-1975b-stability-periodic-orbits-three-body-problem-celest-mech-12-255-doi-10.1007-BF01228563.pdf`. A 22-page
scan with an OCR text layer (PDF page n is journal page 254 + n). I read every page as an image at 130 dpi, and Table I
(journal p.268, printed rotated) at 220 dpi after rotating it, so every digit of the table was read from the image.
No digit was illegible.

Evidence tags: READ (p.N) is read at printed page N. COMPUTED is my own arithmetic or integration of 2026-10-04 (the
scratch scripts are not committed; the recipe is in section 7). INFERRED is my own reasoning, not stated by the author.

Companion: `docs/notes/2026-10-04-digest-hadjidemetriou-1975-restricted-to-general-continuation.md` (the existence theorem,
Celest. Mech. 12:155). This paper is the stability half: it gives the numerical method and a worked family.

## 0. What the paper is, and a naming trap

A method paper plus one worked family. The method (READ, abstract p.255) is "linear and involves the computation of a
4 x 4 variational matrix by integrating numerically the differential equations for time intervals of the order of a
period"; for a symmetric orbit the matrix "can be computed by integrating for half the period only". The worked family
is equal-mass three-body periodic orbits in a rotating frame (two bodies in a binary, the third circling the binary in
the same direction), which the author continued from the restricted problem by raising the third mass
(Hadjidemetriou and Christides 1975, not held). The nonlinear check is by Poincare-section iterates (more than 1000
intersections in some cases).

NAMING TRAP, opposite to the companion paper. Here (READ p.255, p.256 eq. 1-3) the binary is P1, P2 (masses m1, m2,
q = m1/m2), the rotating frame G1xy has its origin at the centre of mass G1 of P1 and P2 with the x axis always through P1
and P2, and P3 (mass m3, coordinates x3, y3 in that frame) is the third body, the one that is massless in the restricted
limit. In the companion paper P2 was the small body. Do not carry equations between the two papers without renaming.

## 1. The reduction, as printed

Lagrangian, READ p.256 eq. (1)-(3):

    L = (1/2)(m1 + m2) { q (x1'^2 + x1^2 th'^2) + (m3/m) [ x3'^2 + y3'^2 + th'^2 (x3^2 + y3^2) + 2 th' (x3 y3' - x3' y3) ] } - V
    V = -G m1 m3 / r13 - G m2 m3 / r23 - G m1 m2 / r12,   m = m1 + m2 + m3,   q = m1/m2

(primes are time derivatives; th is the angle of the G1x axis against a fixed inertial direction; the centre of mass of
the whole system is at rest in the inertial frame; x1 is the abscissa of P1; the abscissa of P2 follows from G1 being
the centre of mass of P1 and P2.)

The angle th is ignorable, so the angular momentum p = dL/d th' is an integral, and (eq. 4)

    th' = [ p/(m1 + m2) - (m3/m)(x3 y3' - x3' y3) ] / [ q x1^2 + (m3/m)(x3^2 + y3^2) ].

Eliminating th gives a three-degree-of-freedom problem in (x1, x3, y3) with the Routhian (eq. 5)

    R = (1/2)(m1 + m2) { q x1'^2 + (m3/m)(x3'^2 + y3'^2) - [ p/(m1+m2) - (m3/m)(x3 y3' - x3' y3) ]^2 / [ q x1^2 + (m3/m)(x3^2 + y3^2) ] } - V,

with an energy integral E = f(x1, x3, y3, x1', x3', y3', p) (eq. 6), a hypersurface in the six-dimensional phase space
(x1, x3, y3, x1', x3', y3'). So the problem is: 8 phase variables including (th, p), minus 2 (th ignorable, p fixed) gives 6,
minus 1 (energy) gives 5, minus 1 (section) gives the 4-dimensional map below.

What the reduced frame does and does not carry (INFERRED from the printed equations): the frame is not uniformly rotating,
its origin G1 is not at rest, and P1 and P2 move on the axis. For m3 -> 0 it collapses to the restricted problem in the
frame rotating with the binary.

## 2. The method, step by step

### 2.1 Isoenergetic Poincare map (READ pp.256-258, eq. 7-13)

A symmetric periodic orbit starts on the section y3 = 0 with x1' = x3' = 0 (eq. 7):
(x1, x3, y3, x1', x3', y3') = (x100, x300, 0, 0, 0, y300). For a neighbouring orbit with the SAME p and the SAME E
(eq. 8): x10 = x100 + xi10, x30 = x300 + xi20, y30 = 0, x1' = xi30, x3' = xi40, with y3'0 fixed by the energy integral.
The next crossing of y3 = 0 in the same direction as initially (eq. 9) defines the four-dimensional map T of the
variables (x1, x3, x1', x3') (eq. 10-11), Q1 = T Q0; fixed points are periodic orbits and the n-th crossing is T^n Q0.
Linearised (eq. 12): xi = A xi0 with A_ij = partial g_i / partial (initial variable j), a 4 x 4 matrix; its determinant is
the Jacobian of T (eq. 13). The linear stability of the orbit, for isoenergetic displacements, depends on the eigenvalues
of A.

Where the trivial unit eigenvalues went (INFERRED, the paper never says so in words): the flow direction is removed by
using a section, the energy direction by holding E fixed and eliminating y3' through eq. (6), and the rotation by the
Routh reduction with p held fixed. Consequently A has no forced unit pair, unlike a flow monodromy matrix, and all four
of its eigenvalues are stability information. The three-dimensional bookkeeping is the restricted problem's: a flow
monodromy of the planar circular restricted problem is 4 x 4 with eigenvalues {1, 1, lambda, 1/lambda}, and the project's
`barden_stability` discards the pair nearest 1; Hadjidemetriou's A is the full-problem analogue with two nontrivial pairs.

### 2.2 Properties of A (READ pp.258-262, eq. 14-38)

1. det A = 1 (eq. 18), proved by passing to canonical variables (q1, q2, p1, p2), where the map is volume preserving
   (Liouville; Siegel and Moser 1971), and noting that the Jacobians relating (x1, x3, x1', x3') to the canonical
   variables at the start and end are equal because start and end points coincide for a periodic orbit. A second proof,
   with no canonical variables, is given in section 6 (eq. 66-67).
2. Reversing symmetry (eq. 19): the transformation x1 -> x1, x3 -> x3, x1' -> -x1', x3' -> -x3', t -> -t leaves the
   system unchanged. Reading T backwards from the point with reversed velocities gives (eq. 21-23)
   A^-1 = J A J with J = diag(1, 1, -1, -1) (eq. 22). With A written in 2 x 2 blocks (eq. 24), A^-1 = [[A1, -A2],[-A3, A4]]
   (eq. 25), and hence (eq. 26) A1^2 - A2 A3 = I, A1 A2 = A2 A4, A3 A1 = A4 A3, A4^2 - A3 A2 = I.
3. Scalar consequences: det A1 = det A4 and trace A1 = trace A4 (eq. 27-29, if det A2 is not zero), so
   a11 + a22 = a33 + a44; and with the permutation matrices K (eq. 30) and M (eq. 34): det[[a11, a13],[a31, a33]] =
   det[[a22, a24],[a42, a44]] (eq. 33) and det[[a11, a14],[a41, a44]] = det[[a22, a23],[a32, a33]] (eq. 36).
4. a_ii = A_ii, the cofactors (eq. 37), because A^-1 = J A J has the same diagonal as A and det A = 1; so trace A is
   the sum of the four 3 x 3 principal minors (eq. 38).

### 2.3 Characteristic polynomial and the stability criterion (READ pp.261-262, eq. 39-49)

    lambda^4 + alpha lambda^3 + beta lambda^2 + gamma lambda + delta = 0   (eq. 39)
    alpha = -trace A,  beta = sum of the six principal 2 x 2 minors of A (eq. 41),  gamma = -(A11 + A22 + A33 + A44),  delta = det A.

By (18) and (38): gamma = alpha and delta = 1, so the equation is reciprocal (eq. 44):
lambda^4 + alpha lambda^3 + beta lambda^2 + alpha lambda + 1 = 0, with roots in reciprocal pairs lambda1 lambda2 = 1,
lambda3 lambda4 = 1 (eq. 45). Factorised (eq. 46-48):

    (lambda^2 + b1 lambda + 1)(lambda^2 + b2 lambda + 1) = 0,
    b1 = (alpha + sqrt(Delta))/2,  b2 = (alpha - sqrt(Delta))/2,  Delta = alpha^2 - 4(beta - 2).

Stability, READ p.262 eq. (49): the orbit is linearly stable when all four roots are complex conjugate on the unit circle,
which happens exactly when

    Delta > 0,   |b1| < 2,   |b2| < 2.

"In all other cases the motion is unstable." The author notes the same quartic arises for the elliptic restricted problem
(Broucke 1969, NASA TR 32-1360, seven stability regions, only one stable) and for the three-dimensional restricted
problem (Bray and Goudas 1967). He states, READ p.262: the condition "Delta < 0" for stability in Bray and Goudas "is not
correct as, in general, Delta < 0 corresponds to unstable motion". That is a published erratum worth knowing if the
project ever relies on Bray and Goudas for a stability condition.

DERIVED by me (COMPUTED, algebra), consistent with eq. (49): each factor lambda^2 + b lambda + 1 has lambda + 1/lambda =
-b, so the one-pair stability index nu = (lambda + 1/lambda)/2 used by `barden_stability` and in Ross's papers is
nu = -b/2, and |b| < 2 is |nu| < 1. The boundaries |b| = 2 are lambda = +1 (b = -2) and lambda = -1 (b = +2); in the
(alpha, beta) plane lambda = 1 is a root iff beta = -2 alpha - 2, lambda = -1 is a root iff beta = 2 alpha - 2, and
Delta = 0 is the parabola beta = alpha^2/4 + 2. Crossing Delta = 0 with |b| < 2 on both sides is a collision of two unit-circle pairs
(a Krein-type event if the pairs have opposite signature, a complex quartet off the unit circle afterwards); crossing
|b| = 2 is the usual tangent or period-doubling event of one pair. This is all derived from the printed quartic, not
printed by the author.

### 2.4 Numerical computation of A (READ pp.263-267, eq. 50-53)

A_ij is computed by finite differences of the map: integrate four neighbouring orbits, each with one initial value
(x10, x30, x1'0, x3'0) increased by a small increment, until y3 returns to zero with the same sign of y3', then take the
ratio of final to initial displacement. READ p.263: "an increment of the order of 10^-5 in the initial conditions ...
gave satisfactory results (the accuracy of the integration was of the order of 10^-11)"; "the elements of A were obtained
with an accuracy of at least four significant figures"; "a smaller or a larger increment resulted in poorer accuracy."
The identities (18), (28), (29), (33), (36), (37) are used to check the matrix. Using them, alpha and beta need only the
elements with i, j = 1, 2, 3 (so three extra integrations suffice):

    alpha = -2 (a11 + a22)  (eq. 51),   beta = 2 { det[[a11,a12],[a21,a22]] + det[[a11,a13],[a31,a33]] + det[[a22,a23],[a32,a33]] }  (eq. 52),
    A44 = a11 + a22 - a33  (eq. 53, an accuracy check).

(The paper labels the reference "the coefficients alpha and beta in (43)"; the definitions are eqs. 40 and 41. Minor
internal cross-reference slip, harmless.)

### 2.5 Half-period computation for symmetric orbits (READ pp.264-267, eq. 54-67)

Let Q0 be on the section at t = 0, Q' the intersection half a period later (opposite direction), Q the next
same-direction intersection (eq. 54-58): xi' = B1 xi0 (start to half period) and xi = B2 xi' (half period to full period).
Then A = B2 B1 (eq. 60). The reversing symmetry (19) gives B2^-1 = J B1 J, i.e. B2 = J B1^-1 J (eq. 63-65), hence

    A = J B1^-1 J B1   (eq. 66),

which needs integration over HALF the period only. det J = 1 gives det A = 1 again, and det B1 = 1 (eq. 67) is proved with
the canonical momenta p1 = q (m1 + m2) x1', p2 = (m3/m)(m1 + m2) x3' at y3 = 0 (eq. 71-72). It is the accuracy check for
B1. Important caveat, READ pp.266-267: "the properties (18), (28), (29), (33) and (37) can no longer be used to check the
accuracy in the computation of the matrix A ... any matrix which is expressed in the form (66) has a determinant equal to
unity and obeys the property (23) ... all the above properties of A will be verified identically no matter how large an
error exists in the computation of B1." That is, once A is built from the half-period formula, determinant and
reciprocity are automatic and say nothing about accuracy.

### 2.6 Stability for non-isoenergetic displacements (READ pp.274-275, section 9)

Linear isoenergetic stability implies stability for general displacements, argued as follows: a displacement that changes
the energy to E' can be read as an isoenergetic displacement from a neighbouring member of the same one-parameter family
(energy varying along the family) that has energy E', and that neighbour is itself isoenergetically stable.
INFERRED gap: a general displacement changes both E and p, two parameters, while the argument uses a family with one
parameter; the argument is explicit about energy only (p varies along the family in Table I, but a displacement that
changes E and p independently is not obviously a displacement from a family member). The nonlinear tests of section 8
all keep E and p equal to the periodic orbit's. Section 9 also says (p.275) that disintegration "may take place after
several hundred intersections (periods) and for this reason great care should be taken in deciding whether or not a
certain motion is stable or not".

### 2.7 How stability changes along the continuation from the restricted to the general problem

The paper does not do the mass continuation itself (it cites Hadjidemetriou and Christides 1975, not held) and states
only (READ p.255, p.267): the orbits were obtained "by varying the mass m3", and "one member of this family has been
obtained by continuing numerically a periodic orbit of the restricted circular three-body problem by varying the mass m3"
(section 7, p.267). It does not say which member of Table I, nor print the continuation path, nor print the restricted
orbit. What it does give is the stability variation along the characteristic curve at fixed equal masses (section 3 below):
stable branch, a Delta = 0 transition to complex-quartet instability, then a real-hyperbolic branch. From the companion
digest the first-order inheritance statement is: the primaries' motion is stable and the small body's stability is that of
the restricted orbit; larger mass needs numerical continuation, with no guarantee the character persists.

## 3. Table I and the family (READ p.268, Table I; Fig. 1 p.269; Fig. 2 p.270)

Normalisation (READ p.267): m1 = m2 = m3 = 1/3, G = 1, and the initial angular velocity of the rotating frame th'0 = 1.
Consequently p varies along the family. The table lists, for 26 symmetric members, the initial conditions on the section
x1' = x3' = 0, y3 = 0 (x10, x30, y3'0 = ydot30), the half-period tau/2, the energy E, the angular momentum p, the
parameters b1 and b2 of eq. (46)-(47), and S (stable) or U (unstable). READ: "Where no values are given for b1 and b2, it
implies that Delta < 0." Fig. 1 shows the characteristic curves (x30 against x10, curve I, and x30 against y3'0, curve II),
with letters A, B, C marking ends and the orbits a to f and 1 to 6 marked on curve I. Fig. 2 shows orbits 1 to 6 in the
rotating frame; the initial conditions of numbers 1 to 6 are "shown in Figure 1" (graph only, not tabulated).
Where the printed table does not pair the numbered orbits with rows, I do not guess; the mapping of orbits 1-6 and a-f to
rows is not printed.

| row | x10 | x30 | y3'0 | tau/2 | E | p | b1 | b2 | S/U |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.09554886 | 1.02242545 | -2.94637505 | 0.076672 | -0.81130995 | 0.36301930 | -1.9977 | -1.9977 | S |
| 2 | 1.00261239 | 0.90242545 | -2.47957163 | 0.123789 | -0.61131972 | 0.35387510 | -1.9918 | -1.9919 | S |
| 3 | 0.87680981 | 0.70242545 | -1.78210271 | 0.292671 | -0.38509299 | 0.34399856 | -1.9266 | -1.9280 | S |
| 4 | 0.80204765 | 0.52242545 | -1.27163802 | 0.640001 | -0.27372385 | 0.34187413 | -1.5065 | -1.5310 | S |
| 5 | 0.77189287 | 0.40242545 | -0.99810693 | 1.092636 | -0.22929139 | 0.34394187 | -0.3685 | -0.6434 | S |
| 6 | 0.76284895 | 0.34242545 | -0.88988142 | 1.504882 | -0.21035787 | 0.34630060 | 1.0163 | 0.0881 | S |
| 7 | 0.76191915 | 0.33242545 | -0.87614909 | 1.610221 | -0.20679462 | 0.34684774 | 1.3755 | 0.2333 | S |
| 8 | 0.76130801 | 0.32242545 | -0.86550686 | 1.743966 | -0.20270082 | 0.34748142 | 1.7757 | 0.4164 | S |
| 9 | 0.76139866 | 0.31242545 | -0.86250490 | 1.950458 | -0.19705229 | 0.34829444 | 1.9276 | 0.9085 | S |
| 10 | 0.76153702 | 0.31122545 | -0.86360433 | 1.989432 | -0.19605400 | 0.34842257 | 1.7211 | 1.1766 | S |
| 11 | 0.76159848 | 0.31082545 | -0.86414896 | 2.004173 | -0.19568097 | 0.34846882 | - | - | U |
| 12 | 0.76167048 | 0.31042545 | -0.86481637 | 2.020124 | -0.19528000 | 0.34851749 | - | - | U |
| 13 | 0.76858510 | 0.31848942 | -0.94718871 | 2.669244 | -0.18052934 | 0.34931891 | - | - | U |
| 14 | 0.77206853 | 0.32848942 | -0.98868862 | 2.875031 | -0.17627753 | 0.34920022 | - | - | U |
| 15 | 0.77273038 | 0.33048942 | -0.99639009 | 2.911089 | -0.17554935 | 0.34916969 | - | - | U |
| 16 | 0.77338585 | 0.33248942 | -1.00395435 | 2.946047 | -0.17484813 | 0.34913823 | -5.0820 | -10.5790 | U |
| 17 | 0.77403622 | 0.33448942 | -1.01139699 | 2.980052 | -0.17417047 | 0.34910617 | -4.5745 | -12.7540 | U |
| 18 | 0.78234159 | 0.36048942 | -1.10090160 | 3.371221 | -0.16669338 | 0.34872540 | -3.3105 | -37.0384 | U |
| 19 | 0.80079171 | 0.41548942 | -1.26881323 | 4.092251 | -0.15444017 | 0.34872335 | -2.9674 | -83.1946 | U |
| 20 | 0.80437480 | 0.42548942 | -1.29768483 | 4.219792 | -0.15247090 | 0.34887701 | -2.9380 | -89.8951 | U |
| 21 | 0.81065194 | 0.44250230 | -1.34605758 | 4.437199 | -0.14924209 | 0.34925431 | -2.9011 | -99.8154 | U |
| 22 | 0.82145194 | 0.47040678 | -1.42383773 | 4.797651 | -0.14422363 | 0.35018878 | -2.8835 | -111.5570 | U |
| 23 | 0.83105194 | 0.49394628 | -1.48839480 | 5.107752 | -0.14021381 | 0.35127507 | -2.9235 | -117.0173 | U |
| 24 | 0.84305194 | 0.52195540 | -1.56444275 | 5.486094 | -0.13566873 | 0.35290618 | -3.0580 | -118.3948 | U |
| 25 | 0.86427822 | 0.56736052 | -1.68973493 | 6.139448 | -0.12860411 | 0.35635235 | -3.5960 | -110.4138 | U |
| 26 | 0.88185998 | 0.60436052 | -1.78720860 | 6.671667 | -0.12347567 | 0.35959209 | -4.4477 | -97.1528 | U |

(The PDF prints digits in groups, e.g. "1.095 548 86"; I have joined them. The row numbers are mine, counted from the top.
Rows 1 to 10 are S, rows 11 to 26 are U. Some rows of the printed table have identical-looking x30 prefixes by design: the
continuation steps x30 by round increments.)

Reading the family, READ p.267: toward end A the period and the energy decrease, the system tends to a binary and the
third body moves around it at larger distance (wording slip below); toward end C the period and the energy increase. "All orbits represented
by the part AB are stable and from there on become unstable. The instability (measured by the magnitude of the eigenvalues
of the matrix A) increases rapidly and the part around the point C represents highly unstable orbits." Wording slip in the
paper (flagged, not a reading problem): in section 7 (p.267) and in orbit a (p.268) the text calls the binary "P1 and P3"
and the circling body "P2", while the setup of sections 1 to 6 has the binary P1, P2 and the circling body P3 (and Fig. 2's
caption says "the motion of P2 on the x-axis ... is symmetric to that of P1", consistent with the binary being P1, P2).
For equal masses the labels are interchangeable physically; use the sections 1-6 convention for any code.

Structure of the stability changes along Table I (COMPUTED from the printed b values, consistent with eq. 49):
- Rows 1 to 4: b1, b2 both real and near -2 (|b| < 2, Delta > 0 small, the pairs close to lambda = +1): stable but
  approaching the |b| = 2 boundary, the near-binary limit at end A.
- Rows 5 to 10: |b| < 2 and Delta > 0; b1 and b2 approach each other from rows 9 (1.9276, 0.9085) to 10 (1.7211, 1.1766).
- Row 11: b values missing, so Delta < 0, and the transition S to U is a Delta = 0 event (two unit-circle pairs
  colliding and leaving the circle as a complex quartet), NOT a |b| = 2 event. My re-computation (section 7) finds for row 11
  Delta = -0.0645 and eigenvalue moduli about 1.096 and 0.912, consistent.
- Rows 16 to 26: Delta > 0 again with both b real and less than -2: all four eigenvalues real and positive (reciprocal
  pairs), one pair of modulus about 10 to 118 (b2 reaching about -118) - "highly unstable" in the author's words.
Row 2: b1 and b2 are within 1e-4 of each other (Delta ~ 0): see the ill-conditioning remark in section 5.

## 4. Nonlinear checks (READ pp.269-275, section 8)

Orbits a to f are members on curve I of Fig. 1, a to e linearly stable and f unstable. For orbit a the author prints the
full 11-digit initial condition (eq. 73, READ p.271):

    x1 = 0.940 821 152 56,  x3 = 0.812 425 450 00,  x1' = 0,  x3' = 0,  y3 = 0,  y3' = -2.150 687 295 28,  th' = 1,

(equal masses 1/3, G = 1), with perturbations Dx1 = Dx3 = Dx1' = Dx3' = eps (eq. 74), eps = 0.02, 0.04, 0.06, 0.10.
READ: the intersections of the perturbed orbits with y3 = 0 in the same direction lie on smooth closed curves (invariant
curves, Fig. 3a x1-x1' plane and Fig. 3b x3-x3' plane), almost ellipses for small eps, distorted for larger; "about 40 to
45 intersections (i.e. revolutions) were needed to obtain a complete invariant curve" (p.272). For the perturbation (75),
Dx1 = Dx1' = 0.05, Dx3 = Dx3' = 0, about 1200 revolutions were computed; the points form ovals that are almost closed:
spiralling out in x1-x1' and spiralling in then out (reversal after about 500 periods, transition at the 9th and 10th
oval) in x3-x3', with a slow rotation (Fig. 4a-b, p.272); the author concludes linear stability implies stability "in
general" for this case. Orbit b (perturbation 0.01 on all four): stable, points slightly diffused, several hundred
intersections. Orbit c: same, more diffusion. Orbit d: diffusion profound, "complete stability". Orbit e: linearly stable
"but on the verge of instability", eps = 0.01 gives about 1000 intersections, diffusion complete, bounded, maximum
deviations observed +-0.05 (x1), +-0.03 (x3), +-0.02 (x1'), +-0.05 (x3'); eps = 0.10 disintegrates (a binary forms
between P1 and P3 and P2 escapes) after about 300 intersections. Orbit f (unstable, eps = 0.01): disintegrates, a binary
forms between P2 and P3 and P1 escapes after about 300 intersections. These are qualitative; the quantitative items usable
as tests are the numbers just listed (counts, bounds) and eq. (73).

Discussion, READ p.275: the family describes triple systems with a close binary and a distant third body, stable when the
third body is farther than the binary's dimensions and unstable when the distance is comparable; "interplay" (Szebehely)
is considered as a perturbed motion of a stable periodic orbit; the existence of bounded points around a stable periodic
orbit "suggest[s] the existence of an additional integral of motion, valid near the periodic orbit". Henon 1974b (private
communication, cited) found a retrograde-third-body equal-mass family stable over a large part.

## 5. Test-ready numbers and how I checked them

All of Table I is printed and usable. I independently re-derived and re-integrated several rows (COMPUTED, section 7 gives
the recipe; the code is an inertial Newtonian equal-mass integrator, not any project model, and does not use the
paper's reduced equations, so the check is independent of a transcription of eq. 1 to 6).

(a) Conventions confirmed against the printed digits. With G = 1, m = 1/3 each, th'0 = 1, P1 = (x10, 0), P2 = (-x10, 0),
P3 = (x30, 0) in the frame rotating with the P1-P2 line, rotating-frame velocities P1, P2 at rest and P3 = (0, y3'0), the
inertial total energy E and the angular momentum p about the system centre of mass of the printed initial conditions are:

| row | E computed | E printed | p computed | p printed |
|---|---|---|---|---|
| 1 | -0.81130984 | -0.81130995 | 0.36301931 | 0.36301930 |
| 2 | -0.61131968 | -0.61131972 | 0.35387510 | 0.35387510 |
| 3 | -0.38509299 | -0.38509299 | 0.34399857 | 0.34399856 |
| 10 | -0.19605401 | -0.19605400 | 0.34842256 | 0.34842257 |
| 24 | -0.13566874 | -0.13566873 | 0.35290618 | 0.35290618 |

(Rows 22 and the other stability rows were run for the matrix in (c) only.) So E is the total energy (G = 1) in the centre-of-mass frame and p the total
angular momentum about the centre of mass, with the P1-P2 line as the frame: eq. (1)-(6) as printed, in these units.
Differences are at the 1e-7 level, the printed 8-digit initial conditions' rounding. (Row 1 E differs by 1.1e-7; I did
not chase it.)

(b) Half period and symmetry. For rows 1, 2, 3, 10, 24 the first return to y3 = 0 (rotating-frame y of P3) in the
integration occurs at t = 0.076672, 0.123789, 0.292671, 1.989432, 5.486094, equal to the printed tau/2 to all printed
digits, and at that instant x1' = x3' = 0 to 1e-9 (the perpendicular-crossing condition). So tau/2 is the time to the
first return of y3 to zero and Table I is internally consistent.

(c) The b parameters. Computing A as a central finite-difference Jacobian (step 1e-5, as the paper) of the isoenergetic
map built from the printed initial conditions, with E and p held at the printed values and the angular velocity and y3'
re-solved at each perturbed point, then alpha = -trace A, beta = sum of principal 2 x 2 minors, Delta, b1, b2:

| row | printed b1, b2 | computed b1, b2 | computed det A |
|---|---|---|---|
| 5 | -0.3685, -0.6434 | -0.36843, -0.64328 | 0.9999994 |
| 6 | 1.0163, 0.0881 | 1.01636, 0.08813 | 0.99998 |
| 9 | 1.9276, 0.9085 | 1.92757, 0.90876 | 1.00001 |
| 10 | 1.7211, 1.1766 | 1.72109, 1.17667 | 0.99983 |
| 11 | none (Delta < 0) | Delta = -0.0645 | 1.00002 |
| 16 | -5.0820, -10.5790 | -5.08171, -10.57983 | 0.99913 |
| 22 | -2.8835, -111.5570 | -2.88360, -111.55060 | 0.99635 |
| 2 | -1.9918, -1.9919 | alpha = -3.98363 (printed b1 + b2 = -3.9837) | 0.99996 |

So Table I's stability parameters are reproduced to 4 significant figures for the modest rows, and to 3 for rows 16 and 22
(the unstable rows with a multiplier of order 1e2: the finite-difference determinant departs from 1 by 1e-3 to 4e-3, the
size of the difference). All sign, S/U and Delta < 0 labels checked agree. Row 2 is a warning: Delta is about 1.6e-4, so
the pair of numbers (b1, b2) splits as +-sqrt(Delta)/2 around alpha/2 and the split is ill-conditioned (computed
-1.9854 and -1.9982 against the printed -1.9918 and -1.9919) while their sum alpha is reproduced to 7e-5. A test on rows 1
or 2 must pin alpha (and beta), not b1 and b2 separately, or use a tolerance of about 1e-2 on each.

(d) Orbit a, eq. (73). COMPUTED, not printed: E = -0.49376285, p = 0.34848785, full period 0.3624990 (half period 0.1812495:
the first symmetric return is at x1 = 0.87662330, x3 = 1.00501900 with x1' = x3' = 0, and the second return, one full
period, is back at the printed x1 = 0.94082115, x3 = 0.81242545 to 1e-9, repeating for at least five periods). The printed
11-digit state is therefore a genuine symmetric periodic orbit of the equal-mass problem; and it lies between Table I rows 2
and 3 in E, tau/2, p, x10 and y3'0 (x30 = 0.81242545 is on the same x30 grid with step 0.1, 0.9024 and 0.7024 at rows 2 and
3). The map's eigenvalues computed for orbit a: alpha = -3.95656, beta = 5.91361, Delta = -5.5e-5 (zero within the
finite-difference noise), eigenvalues about 0.98919 +- 0.14657 i and 0.98909 +- 0.14741 i, both of modulus 0.99999 and
argument about 0.147 rad, i.e. 2 pi / 0.147 = 42.7 intersections per circuit: consistent with the printed "about 40 to
45 intersections" per invariant curve (p.272). Caution: the printed criterion Delta > 0 would label this orbit unstable
on the numerical sign of Delta, though it is stable; Delta is within noise of zero because the two pair angles are nearly
equal (0.1466 and 0.1474). Use modulus of the eigenvalues (all |lambda| = 1 within tolerance), not the sign of Delta,
in code.

## 6. Reconciliation against the project's code

Searched `src/cyclerfinder` for monodromy, floquet, stability_index, eigenstructure, barden, isoenergetic, routh and
hadjidemetriou. Hits: `search/cr3bp_periodic.py::barden_stability`, `search/er3bp_floquet.py::{er3bp_monodromy,
floquet_classify}`, `search/er3bp_periodic.py::monodromy_eigenstructure`, `search/sun_forced_periodic_884.py::monodromy`,
`search/jovian_resonant_families.py` (with its own `_planar_floquet` full-period cross-check), `search/binary_star_search.py`,
`core/cr3bp.py::propagate(with_stm=True)`. No isoenergetic section map, no Routhian reduction, and no `hadjidemetriou`
anywhere in `src/` or `tests/`.

Models: `core/cr3bp.py`, `er3bp.py` (+ `er3bp_geocentric.py`, `er3bp_paper_frame.py`), `bcr4bp.py`, `qbcp.py`,
`ccr4bp*.py`, `crnbp.py` are all restricted-problem models: the third body is massless and every other body's motion is
prescribed (circular, elliptic or ephemeris). `nbody/` is a REBOUND restricted integrator ("RestrictedNBody", massless
spacecraft, planets on rails from DE440). There is no finite-mass three-body (or general N-body) Newtonian code in `src/`,
so the paper's own setting (all three bodies massive, recoil, angular-momentum and energy reduction) has no project
model. This is the main reconciliation fact: the paper's numbers cannot be reproduced by any current module.

| Paper | Project | Status |
|---|---|---|
| 4 x 4 isoenergetic map A, eigenvalues of the section map | Flow monodromy 6 x 6 (`core/cr3bp.propagate`) or planar 4 x 4 sub-block, with the trivial pair discarded by hand | Same information; the paper's A has no trivial pair. Difference is bookkeeping, not physics. |
| Half-period formula A = J B1^-1 J B1 (eq. 66) | `barden_stability`: M = G Phi(T/2)^-1 G Phi(T/2) with G = diag(1, -1, -1, 1) on (x, y, vx, vy) | Same construction. G is the reflection in the project's variables; the paper's J = diag(1, 1, -1, -1) is the same time-reversal in (x1, x3, x1', x3'). As the paper warns, det M = 1 and reciprocity hold identically, so they are no check; the project already has an independent full-period cross-check in `jovian_resonant_families._planar_floquet` (agreement 3.0e-7). Any new use of `barden_stability` for a strongly unstable orbit should keep such a check. |
| One-pair index nu = (lambda + 1/lambda)/2, stable iff abs(nu) < 1 | `barden_stability` returns nu and lambda | Same quantity, nu = -b/2. Valid only when the map has ONE nontrivial pair (planar circular restricted). |
| Two pairs, criterion (49): Delta > 0, abs(b1) < 2, abs(b2) < 2 | `floquet_classify`: unstable iff max abs(lambda) > 1 + 1e-3 | Equivalent in exact arithmetic; the eigenvalue-modulus form is the numerically safer one (section 5(d): Delta's sign is unreliable near a double angle). `floquet_classify`'s "stable" branch (all abs(lambda) < 1 - tol) cannot occur for a symplectic spectrum; "marginal" is the real stable case. Not a defect, but the tag names are misleading. |
| Complex-quartet instability (Delta < 0, abs(lambda) off the unit circle, Table I rows 11 to 15) | `floquet_classify` catches it by modulus; `monodromy_eigenstructure` takes the largest modulus as a saddle and a unit-circle complex pair as a centre | `monodromy_eigenstructure` would, on a quartet (modulus about 1.1 for both members of one complex pair), classify the larger as r and then look for a centre among moduli within 0.5 of 1: it returns a number rather than raising, but the labels "saddle" and "centre" would be wrong. (This is `data/OUTSTANDING.md` item `#912`'s note that it "would raise on three real pairs"; the quartet is a second case.) Check before relying on it for elliptic-problem classification. |
| The criterion needs det A = 1 and A^-1 = J A J (symmetric orbit) | none in `src/`, no test of symplectic or reversing identities on monodromies | Cheap controls missing, see the next section. |

The project's own restricted models and the paper's A: the paper's section map is for the finite-mass problem. In the project's
models the same statements hold with the restricted flow's own reversing symmetry (CR3BP planar: (t, x, y, vx, vy) ->
(-t, x, -y, -vx, vy)), used already in `search/cr3bp_periodic.py` and `search/two_moon_periodic_890.py` (docstring of the
latter). In a time-periodic restricted model (ER3BP, BCR4BP, QBCP, prescribed-moon four-body) there is NO energy integral
and NO unit-eigenvalue pair guaranteed, so the paper's isoenergetic reduction does not apply; the 4 x 4 planar monodromy
of such a model has the same quartic form (44) because it is symplectic (Broucke 1969), so (49) is the right classifier
there (that is the case the paper itself notes, p.262).

Companion-paper reconciliation: the companion digest's first-order statement "stability inherited from the restricted
orbit, primaries stable" is consistent with this paper's section 7: the stable region of the equal-mass family (rows 1 to
10) is the one with the third body well outside the binary.

## 7. Recipe for the independent re-computation (so it can be repeated or turned into a test)

Inertial Cartesian Newtonian integration, G = 1, masses (1/3, 1/3, 1/3), DOP853 with rtol 1e-12 or tighter (1e-13 and
atol 1e-14 for the stability rows). At t = 0 in the frame rotating with the P1-P2 line (angle th = 0): P1 = (x10, 0),
P2 = (-x10, 0), P3 = (x30, 0), rotating-frame velocities (0, 0), (0, 0), (0, y3'0); inertial velocity = rotating velocity +
th'0 times the 90-degree rotation of the position vector with th'0 = 1; then subtract the centre of mass position and
velocity. E = total energy; p = sum m_i (x_i vy_i - y_i vx_i). At later times recover the rotating frame by the angle of
(P1 - P2), the frame origin at the centre of mass of P1 and P2, and the "section" y3 = 0 as the rotating-frame y of P3.
For A: perturb one of (x10, x30, x1'0, x3'0) by +-1e-5 (central differences), re-solve th'0 and y3'0 by a root finder so that
E and p equal the periodic orbit's (the paper's isoenergetic condition; branch chosen near th' = 1 and the printed y3'0),
integrate to the first same-direction return of y3 = 0 near tau, and read (x1, x3, x1', x3') in the rotating frame, with
the velocities from a central difference of the rotating-frame positions (1e-7 step). Then alpha, beta, Delta, b1, b2 as in
eq. (41), (47), (48). This took tens of seconds per row.

## 8. Techniques applicable to the project's problems

Entries read from `data/OUTSTANDING.md` before writing: `#890` (lines near 591, RESULT 2026-10-04), `#895` (near 1406),
`#899` (near 1247, corrections and scoping), `#912` (near 1173), `#925` (near 1125), `#916` and `#917` (near 1016), `#924`
(near 1060), `#928` (near 1078). The `#897` synthesis note does not name this paper or its companion (searched), so no
earlier synthesis considered Hadjidemetriou's continuation; the companion digest and the Musielak-Quarles review did, for
`#890`, and the companion digest withdrew that application (wrong parameter: it continues the particle's mass, whereas
`#890` raised the moons' masses with a massless spacecraft).

### 8.1 `#890` and `#895` (Titania-Oberon; methods only)

What matters here is the moons' mutual gravity and recoil, which both project models drop in `#890` (prescribed circles,
no moon-moon force) and keep only through DE-kernel positions in `#895` (prescribed, not responding to the spacecraft; the
spacecraft mass is negligible, about 1e-21 of a moon's, so recoil is irrelevant there; moon-moon force enters through
the kernel). The paper's contribution is therefore not a model but a method: (i) classify the stability of a periodic
orbit of a system with several massive bodies from a four-dimensional section map, no trivial pair, using (49); (ii) the
half-period formula; (iii) the warning that identities obtained from the half-period formula cannot validate it.

Concrete tests:
1. Reported multiplier 8.4e5 for the `#890` orbit (largest Floquet multiplier). Independent control: integrate the full
   cycle (not the half) with the STM, check det(M) = 1 within 1e-6 (relative, with the multiplier near 1e6 the
   finite-difference route is hopeless; use the variational STM), check the multiplier appears with its reciprocal (the
   pairing eq. 45) and that M^-1 = G M G with the model's reflection. If the reported multiplier came from a half-period
   Barden construction, these three identities hold automatically and the check is vacuous; then the independent control
   must be a full-period STM or a full-period finite difference. State which was used in the `#890` and `#895` notes.
2. Control for any new full-period monodromy code: Table I rows 5, 6, 10 (stable) and 16, 22 (unstable) reproduce to 4
   (3 for the unstable) significant figures; see section 5(c). This does not validate a restricted model, only the
   method code; it is a sourced control for the method, which `#890` lacks.
3. If a future `#890` follow-up wants to treat Oberon and Titania as massive and mutually interacting (a genuine
   four-body periodic orbit with Uranus as the centre), the reduction is Hadjidemetriou's: remove the rotation by the Routh
   reduction at fixed angular momentum, fix E, take a section; the paper treats N = 3
   only, so the reduced map's dimension for four bodies is not derived there and I do not state it.

### 8.2 `#899` (second-species seeds and continuation through close approaches)

The paper's family structure is the useful part: along a monoparametric family at fixed masses the stability type changes
only at (a) Delta = 0 (complex quartet, Table I row 10 to 11) and (b) |b| = 2 (one pair crossing lambda = +1 or -1,
visible where b passes through -2, which is a tangent or period-doubling bifurcation). For a mass or Jacobi-constant
continuation of a Casoliva-type orbit, record (alpha, beta, Delta, b1, b2) of the planar monodromy at every step. A change
of regime at a step is then identified rather than guessed, and a failed continuation ending at lunar impact (Casoliva et
al.) can be compared with whether it ended at |b| = 2 (a bifurcation to follow) or at Delta = 0 (a quartet, no new
family) or at neither (an impact). Control: Table I's own S to U transition at rows 10 to 11 where the program must flag
Delta = 0 and not |b| = 2. The paper does NOT cover close approaches (its section map is smooth only away from
collision), so it adds nothing to the second-species seeding itself.

### 8.3 `#912` and `#925` (elliptic-problem corrector; controls that do not reproduce)

For `#912`, the segmented monodromy product needs a classifier that is valid for two or three reciprocal pairs including
complex quartets. The planar elliptic problem has a 4 x 4 symplectic monodromy (no energy integral), so (44)-(49) apply
exactly (this is the case Broucke 1969 treated, which the author cites): classify by (alpha, beta, Delta, b1, b2), with
the modulus test as the decision and b1, b2 as the labels. Test: build a symplectic 4 x 4 matrix from a prescribed
spectrum (two circle pairs; a circle pair and a real pair; a complex quartet; two real pairs) and require the classifier
to return the right regime in each; Table I rows 1 to 4, 6 and 11, 16, supply the sourced regimes (quartet, real-real,
unit-circle-circle) with printed b values, but only as numbers (not matrices), so the sourced control is the full-problem
reproduction of section 5(c) through a test-only integrator, and the matrix-level test is a synthetic one.

For `#925` (two Mako-Salamon unstable bands not reproduced, Fig. 7 [178, 215] deg, Fig. 8 [119, 259] deg): add one
diagnostic before attributing the miss to the model or to the paper. At the printed unstable true anomalies, compute
(alpha, beta, Delta, b1, b2) of the planar monodromy, not only the largest modulus. If the printed unstable bands are
complex-quartet instabilities (Delta < 0) of modest growth (modulus about 1.01 to 1.1), a one-pair stability index or a
real-eigenvalue-only test would call them stable while a modulus test would not; `floquet_classify` is a modulus test, so
for it the miss would not be a classifier artefact, but a test written with `barden_stability`'s single index would be.
This is a hypothesis (INFERRED), not a finding: whichever classifier the xfail tests use, the (alpha, beta, Delta) values at
the missing bands are the cheapest discriminating evidence, and the paper's Bray-Goudas erratum (sign of Delta, p.262) is
the precedent for a published literature slip of exactly this kind.

### 8.4 `#916` and `#917` (invariant curves of the one-period map; Casoliva rows in the elliptic problem)

The paper's section 8 is a worked example of the object `#916` wants: the one-period map in the neighbourhood of a stable
periodic orbit, iterated about 40 to 45 times per closed invariant curve, with a measured relationship between the map's
eigenvalue angle and the iteration count. Positive control for any such code: orbit a (eq. 73, 11 digits printed):
start from the printed state and perturbation (74) with eps = 0.02 and expect closed invariant curves in the (x1, x1') and
(x3, x3') projections for eps = 0.02, 0.04, 0.06, 0.10 (Fig. 3), and a rotation of about 2 pi / 0.147 = 42.7 iterations per
circuit from the eigenvalue angle (section 5(d), computed), against the printed 40 to 45. Failure control: perturbation
(75) (0.05 on x1 and x1') is expected to show nearly closed ovals that spiral and slowly rotate over about 1200
revolutions (Fig. 4); the same code must not report that case as perfectly closed, nor as escaping. Second control
(near the stability boundary): orbit e with eps = 0.01 bounded diffusion with the printed maximum deviations (+-0.05,
+-0.03, +-0.02, +-0.05), eps = 0.10 escape after about 300 intersections. The paper does not tabulate orbits b to f, so
only orbit a is fully specified; e's state is not printed (only its position on Fig. 1), which means the e checks can be
done only as "a table member between rows 9 and 10 with Delta near 0".

For `#917` (the elliptic problem with the fold): the paper supplies the invariant-curve picture only for the finite-mass
circular family; it says nothing about elliptic perturbation or folds.

### 8.5 `#924` (FLI, MEGNO as a screen)

The paper's evidence is an independent, if informal, positive control for a chaos indicator: orbits a to d are
regular (closed invariant curves or diffuse-but-bounded ovals), e is bounded with complete diffusion (mixed), and f and
the unstable rows disintegrate. An FLI or MEGNO screen run on Table I orbits with the printed perturbations should rank
a < b < c < d < e < f in increasing indicator; the printed text gives this order qualitatively (diffusion "more evident",
"profound", "complete"). This is a rank-order control, not a number-for-number reproduction, and rests on a
finite-mass code the project does not have; the cheapest source of a numeric control remains the standard map for FLI and
the logarithmic potential for MEGNO already named in `#924`.

### 8.6 `#928` (regularised propagator and transition matrix for close passes)

The paper contributes nothing to regularisation. It does bear on the validation of a regularised STM: the controls the
paper uses (det A = 1, eq. 18; reciprocity; A^-1 = J A J for a symmetric orbit; accuracy 1e-11 integration with 1e-5
increments giving 4 significant figures) are the standard checks, and Table I is a sourced check of the final classifier
(section 5(c)), not of a regularised STM. Note the paper's own caveat (p.263) that finite-difference increments both too small
and too large degrade A: a regularised STM that is analytic would remove that dependence, which is an argument for the
variational STM that `#928` plans.

## 9. Recommended follow-ups (no task numbers registered)

1. A test-only finite-mass equal-mass three-body Newtonian module (about 60 lines) and a test file that pins Table I rows 1,
   2, 3, 10 and 24 for E, p and tau/2 and the perpendicular crossing, rows 5, 6, 10, 16, 22 for (b1, b2) within 2e-3 absolute (3e-3
   relative for the unstable rows), row 11 for Delta < 0, rows 1 and 2 pinned on alpha only, and eq. (73) for symmetry
   (return to the printed state after a computed period of about 0.362499). Expected values are the printed ones; my
   re-integration is only evidence that they are reproducible, not the oracle (no circular golden).
2. A classifier function for a 4 x 4 symplectic monodromy that returns (alpha, beta, Delta, b1, b2, regime) and decides
   by eigenvalue modulus, with synthetic-spectrum tests for each of the four regimes; use it in the elliptic problem
   (`#912`, `#925`) and keep `monodromy_eigenstructure` away from complex quartets until it is checked on one.
3. Add two checks to every reported large multiplier (the `#890` 8.4e5, jovian families): full-period STM determinant
   and reciprocity, and state whether the reported value came from the half-period construction.
4. Compute (alpha, beta, Delta) at the Mako-Salamon missing unstable bands before attributing the miss to model or paper.
5. Acquire Hadjidemetriou and Christides 1975 (the mass continuation of Table I's first member; not held) and Henon 1974a,b,
   and Broucke 1969 (NASA TR 32-1360; the full regime catalogue for the quartic), then decide whether the restricted to
   general continuation is worth a project lane; the companion digest concludes it is not needed at the project's
   spacecraft mass ratios, so this is background, not a task.
6. If `#890` or `#895` ever treats the moons as responding bodies, derive the reduced map dimension for four bodies from
   the Routh reduction; the paper stops at three.

## 10. Summary

- Method: a 4 x 4 isoenergetic Poincare map A on the section y3 = 0 of the Routh-reduced three-body problem (angular
  momentum and energy fixed), computed by finite differences of full-period integrations, or, for symmetric orbits, by
  A = J B1^-1 J B1 from half a period. det A = 1 and A^-1 = J A J give a reciprocal quartic; stability iff Delta > 0 and
  |b1|, |b2| < 2 (all four eigenvalues on the unit circle). The half-period route makes those identities vacuous as
  accuracy checks.
- Results: a 26-member equal-mass family (Table I, 8-digit initial conditions, tau/2, E, p, b1, b2, S/U), stable for the
  third body outside the binary, a Delta = 0 transition to a complex-quartet instability, then strongly real-hyperbolic;
  nonlinear checks by Poincare-section iteration (40 to 45 points per invariant curve, 1200 revolutions on one, disintegration
  after about 300 intersections for the unstable and verge orbits); linear stability argued to imply stability for general
  displacements (with the gap noted in section 2.6).
- Independent re-computation (COMPUTED): E, p, tau/2 of rows 1, 2, 3, 10, 24 and the b parameters of rows 5, 6, 9, 10, 16, 22
  reproduce the printed values (4 significant figures, 3 for the unstable rows), row 11 is Delta < 0, and orbit a (eq. 73) is
  periodic with period 0.362499 and eigenvalue angle 0.147 rad (42.7 iterations per circuit).
- Not a model the project has: no finite-mass three-body code in `src/`; the paper's numbers need a new test-only integrator.

## Note 2026-10-05: the mass-continuation paper is now held and digested

Hadjidemetriou and Christides 1975 (Celest. Mech. 12:175, DOI 10.1007/BF01230210) is digested in
`docs/notes/2026-10-04-digest-hadjidemetriou-christides-1975-families-planar-three-body.md`. It corrects three statements
above. (1) It is NOT the source of this paper's equal-mass family: its printed path (x30 = 0.181, m1 = m2, total mass 1)
turns back at m3 = 0.2245 and never reaches m3 = 1/3, and its x30 is 0.181 against 0.31 to 1.02 in Table I here; the
starting orbit and path of the 1975b member remain unprinted. (2) The convention check in section 5(a) is confirmed from a
second source: all 19 rows of that paper's Table I reproduce with the same inertial integrator and normalisation (G = 1,
total mass 1, theta'0 = 1) to 2e-8 in x1(T), x3(T) and 7e-8 in theta/2pi, with one probable misprint (row 15 period). (3) Section 6 says no current module can reproduce the paper's numbers; that holds for finite third mass, but the project's CR3BP at mu = 1/2
reproduces that paper's m3 = 0 starting orbit to eight digits (T = 0.9242291, x = 0.721018392 from x30 = 0.181,
y3'0 = -0.94711034), so it controls the restricted end of a continuation; no module controls finite m3.

## Note 2026-10-05: Henon 1974a is now held and digested

Henon, "Families of periodic orbits in the three-body problem", Celest. Mech. 10:375 (DOI 10.1007/BF01586865) is digested in
`docs/notes/2026-10-04-digest-henon-1974a-families-periodic-orbits-three-body.md`. Relevant here: its two families have masses 3:4:5
(E = -47/288 in its normalisation), so none of its orbits is, or connects through anything printed to, an orbit of this paper's
Table I (equal masses) or of Hadjidemetriou and Christides (m1 = m2, m3 varied). Henon states (p.385) that the Hadjidemetriou-Christides sequence is
a section of a three-parameter family (m2/m1 and x30 constant), and its dimensional count (periodic orbits at fixed masses form
one-parameter families) is the counting behind the continuation used here. An independent re-integration reproduces all 13 non-collision rows
of its Tables I and II (E, A, phi to 1e-9, 5e-9, 8e-7), a third sourced control set for `#931` (a); the closure test there needs the rotation by -phi.
