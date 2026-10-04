# Digest: Hadjidemetriou 1975, continuation of periodic orbits from the restricted to the general three-body problem

Date: 2026-10-04 (Sydney). Reading and reasoning only; no code was run for this digest. No existing file was edited.

Source: J. D. Hadjidemetriou, "The continuation of periodic orbits from the restricted to the general three-body
problem", Celestial Mechanics 12:155-174 (1975), DOI 10.1007/BF01230209, received 21 December 1973. Filed in the private
paper corpus as `hadjidemetriou-1975-continuation-periodic-orbits-restricted-to-general-three-body-problem-celest-mech-12-155-doi-10.1007-BF01230209.pdf`.
The file is a 20-page scan of the whole paper (pp.155-174). The text layer is jumbled, so I read every page as an image,
including the formulas.

Evidence tags: READ (p.N) is read in the paper at printed page N. COMPUTED is a small calculation of mine. INFERRED is
my own reasoning, not stated by the author.

## 0. The answer to the question that was asked

The review digest (`docs/notes/2026-10-04-digest-musielak-quarles-2014-three-body-problem-review.md`, section 2.4) calls
this paper "the existence theory for the step the project takes by hand in `#890`". That claim is wrong, and the reason
is the question of which mass is continued.

The paper's parameter is the mass of the THIRD body, that is, the particle that is massless in the restricted problem. It
is option (i), not option (ii). The two large bodies are kept at fixed masses and the particle's own mass m2 is raised
from zero to a small positive number; the large bodies then recoil and no longer move on exact circles. The `#890`
continuation did the opposite: the spacecraft stayed massless throughout and the two moons' masses were raised from 1
percent to 100 percent of their physical values. The theorem does not cover that.

Quotations, READ:

- Abstract (p.155): "It is proved that a symmetric periodic orbit of the circular planar restricted three-body problem
  can be continued analytically, when the mass of the third body is small but not negligible, to a periodic motion of the
  general three-body problem in a rotating frame of reference whose origin coincides with the center of mass of the two
  bodies with large masses and its x axis always contains these bodies."
- Section 8 (p.169): "We assume that the masses m1, m3 are fixed and we consider m2 as a parameter."
- Section 7 (p.167): "We allow now the mass m2 of P2 to be finite, but small, so that second order terms in m2 to be
  neglected."

A naming trap, READ (pp.155-156): in this paper the bodies are P1, P2, P3 with masses m1, m2, m3, and P2 (mass m2) is the
SMALL third body, while P1 and P3 are the two primaries. The mass-ratio parameter of the primaries is mu = m3/(m1 + m3)
(eq. 10, p.157) and q = m3/m1 (eq. 6). A reader used to the convention in which "m2" is the secondary will misread every
equation. The small mass is m2 throughout.

## 1. Setting, exactly as printed

Planar problem of three bodies, centre of mass at rest in an inertial frame (p.155-156). G1 is the centre of mass of
P1 and P3. The rotating frame G1xy has its origin at G1 and its x axis always through P1 and P3 (Figs 1, 2, p.156-157).
Generalised coordinates: r = distance P1P3, theta = polar angle of P3, and the Cartesian (x2, y2) of P2 in the rotating
frame (eq. 9). Positions of the primaries on the x axis: x1 = -mu r, x3 = (1 - mu) r (eq. 14).

The angle theta is ignorable, and the angular momentum integral p_theta (eq. 12-13, pp.157-158) removes it, leaving a
three-degree-of-freedom problem in (r, x2, y2). Unlike the restricted problem, the frame does not rotate uniformly, the
origin G1 is not at rest, and the primaries' positions on the axis are not fixed (p.158). Units: k^2 = 1 and
m1 + m3 = 1 (p.158).

The equations of motion are eqs (15)-(17) (p.158), with the Lagrangian rather than canonical form chosen "because in the
latter case the small mass appears, formally, in the denominators in the expressions for the momenta" (p.155). Expanding
to first order in m2 (and in rdot, which is taken to be of order m2 because the primaries start on a circle) gives
eqs (29)-(31) (p.159). At m2 = 0 they reduce to r = r0 with r0 = p^2 and the restricted problem in the frame rotating at
omega = r0^(-3/2) (eqs 33-35, p.160). The unperturbed primaries are on a circular orbit; that is a hypothesis (p.159,
p.169: "r0 = 1, rdot0 = 0").

## 2. The method, step by step

The proof has two layers. The first (sections 3 to 7) is a direct construction at first order in m2. The second
(section 8) is a standard implicit-function argument that needs no symmetry. They agree.

### 2.1 Linearisation about the restricted orbit (sections 4 to 6, pp.160-167)

The solution is written r = r0 + r1(t), x2 = x2(t) + xi1(t), y2 = y2(t) + xi2(t) with the corrections of order m2
(eq. 36, p.160), and the variational equations (37)-(38) follow. The restricted-problem orbit x2(t), y2(t) is assumed
symmetric about the x axis: x2(t) = x2(-t) and y2(t) = -y2(-t) (p.161). Then x2 has a cosine Fourier series and y2 a sine
series, and every forcing term that appears has a definite parity (eq. 39).

Primaries (section 5, pp.161-162). The equation for the separation correction decouples:
r1'' + omega^2 r1 = m2 f(t), with f an even function of period T (eq. 40). Its solution (eq. 46) contains a
homogeneous part at the primaries' own frequency omega, and a particular part with Fourier coefficients
a_n = beta_n / (omega^2 - omega0^2 n^2), omega0 = 2 pi / T (eq. 45), "provided 2 pi n / T is not equal to omega" (p.161).
The choice r1(0) = m2 sum beta_n/(omega^2 - omega0^2 n^2), rdot1(0) = 0 (eq. 48) removes the homogeneous part, so that
r = r0 + m2 phi(t) (eq. 49): the distance between the primaries is periodic with the same period T as the particle's
orbit.

Small body (section 6, pp.162-167). The four-component first-order system (eq. 50-52) is mapped by a linear
transformation S1 (eq. 54-55) to the Hamiltonian form of the restricted problem's variational equations (eq. 56-57). Its
four basic solutions (eq. 59) are two with period T (one is the flow direction, the other carries a secular term t) and a
complex pair exp(+/- i a t) times T-periodic functions, where +/- i a are the nonzero characteristic exponents, taken to
be of the stable type for this construction (p.164). The forced solution is built by variation of constants following the
author's 1973 paper (his "Paper A", cited as J. Mecanique 12:1, not held). The constant-of-motion bookkeeping gives
eq. (65), and the freedom in the initial constants is used to cancel every secular term: lambda2 = -m2 k and
lambda3 = i m2 N30* (eq. 66). A condition is printed here: the T-periodicity of the functions N3 requires that
"2 pi / (a T) is not a rational number" (p.165). The result is a T-periodic orbit of P2 (eq. 67-68, p.166).

The free constant lambda1 is shown to be only a phase shift along the orbit (eq. 69-70, p.166-167). Setting lambda1 = 0
makes the right-hand sides of (70) even (x) and odd (y), so the continued orbit "is symmetric with respect to the x axis
of the rotating frame" (p.167).

### 2.2 The families statement (section 7, pp.167-169)

READ, p.167 (italics in the original): "the motion of P1, P3 on the rotating x axis is periodic with period T and the
motion of the body P2 with the small mass in the rotating frame of reference G1xy, whose x axis contains always the
bodies P1 and P3, is also periodic with the same period T. This orbit is symmetric with respect to the x axis of the
rotating frame, if we choose lambda1 = 0".

If lambda1, lambda2, lambda3 are left free the motion of P2 is quasi-periodic with periods T and 2 pi / a (p.168).

For a fixed mass ratio q the symmetric orbit belongs to a one-parameter family, a curve in the (x0, C) plane (initial
crossing point and Jacobi constant). Since the correction functions depend only on x0 and C, the continued solution is
x2' = x2'(x0, C, m2; t) (eq. 72-73), and the family relation x0 = x0(C) (eq. 74) carries over. READ, p.168 (italics):
"a family of symmetric periodic orbits of the restricted three-body problem can be continued analytically, for a fixed
value of the small mass m2, to a family of symmetric periodic orbits in the general three-body problem."

Two further statements, READ (p.169). First, the argument was given for a stable restricted orbit; "The same result can
be obtained if we start from an unstable periodic orbit of the restricted problem. The reasoning is completely analogous
to the one used above and we shall not repeat it here." So the unstable case, which is the cycler case, is asserted, not
proved in the text. Second, "The method used in this paper holds only for small values of m2. The existence of periodic
orbits in the general problem of three bodies for larger values of m2 can only be established by numerical continuation
of the above periodic orbits."

### 2.3 The implicit-function argument (section 8, pp.169-171)

This section is the general theorem, and it states that "No assumption will be made here for the symmetry of the orbit"
(p.169). The periodicity conditions are the six equations r(...) - r0 = 0 and five similar ones, in the six variables
(r, rdot, x2, y2, xdot2, ydot2) at fixed m1, m3, with m2 as the parameter. The Jacobian of these conditions at m2 = 0 is

    J(T) = [[A1(T), 0], [A2(T), Delta(T)]] - I6,

where A1(T) = [[cos wT, sin(wT)/w], [-w sin wT, cos wT]] is the primaries' radial oscillation (epicycle at frequency w),
Delta(T) is the 4 x 4 monodromy matrix of the restricted orbit, and the zero block is 2 x 4 (p.169-170). Then
det J(T) = 2 (1 - cos wT) det(Delta(T) - I4), "always equal to zero" because Delta has a unit eigenvalue (p.170). My check
of the first factor: det(A1 - I2) = (cos wT - 1)^2 + sin^2 wT = 2 - 2 cos wT (COMPUTED, matches).

The standard cure is used: the energy integral removes one condition. Keep y20 fixed (a phase condition), vary the
other five initial values, and delete the last row and fourth column to get J*. Then
det J*(T) = 2 (1 - cos wT) det Delta*(T), and "det Delta*(T) is not zero, provided the periodic orbit of the restricted
problem under consideration belongs to a monoparametric family of periodic orbits, with the period varying along the
members of the family" (p.170). The sixth periodicity condition is recovered because dE/d(ydot2) is nonzero once
m2 is nonzero, although it vanishes at m2 = 0 only through the common factor m2 (p.171).

The theorem, READ (p.171, italics): "any periodic orbit of the restricted three-body problem which belongs to a
monoparametric family of periodic orbits can be continued to a periodic motion of the general problem, with the same
period, provided its period is not equal to 2 k pi (since in the usual normalized units w = 1). The explicit form of the
periodic solution for m2 not zero is given, to a linear approximation, by Equations (49) and (70)."

Stated hypotheses, collected (all READ):

1. Planar problem; primaries circular at m2 = 0 (p.159, p.169).
2. m2 small, first order in m2 (p.159, p.167; "second and higher order terms" are dropped, p.155).
3. The restricted orbit lies in a monoparametric family with the period varying along it (det Delta* not zero, p.170).
4. w T is not a multiple of 2 pi, equivalently 2 pi n / T is not equal to w for any n (pp.161, 170-171).
5. For the explicit construction: the nonzero exponents are stable, and 2 pi / (a T) is irrational (p.165). For the
   stability discussion: T, 2 pi / w and 2 pi / a are incommensurable (p.172).

## 3. Stability result (section 9, pp.171-173)

READ, p.172 (italics): "The motion of the two primaries P1, P3 on the rotating x axis is always stable. The motion of
the body P2 with the small mass is stable or unstable, depending on whether or not the nonzero characteristic exponents
of the reference periodic orbit of the restricted problem are of the stable or unstable type, respectively."

The primaries' separation is bounded for every choice of orbit and initial data (eq. 46), given that 2 pi / w and T are not
commensurable (p.171). The secular terms k t in the general solution (eq. 75) are shown to be only a rescaling of the
phase along the orbit, x2((1 + k p) t) (eq. 76), of the same nature as the secular terms of the restricted variational
equations, so they do not cause instability. Instability can arise only from a real characteristic exponent of the
reference orbit.

Scope, READ (p.172-173): this is first-order orbital stability, and it is stability in the rotating frame G1xy. If P2's
motion is unstable, the primaries "will deviate much from their original orbit, with respect to an inertial frame, as can
be seen from (32a)", where GG1 = -(m2/m) r2 shows the displacement of G1 from the system's centre of mass is driven by
P2's position.

## 4. Numerical example

None. READ absence: the paper has no numerical orbit, table, figure of an orbit, or printed initial condition. The only
figures are the two coordinate diagrams. Consequently no printed number can be turned into a check. The only checkable
statements are structural identities, listed in section 6.

## 5. What this does and does not cover for the project

### 5.1 For `#890` (two-moon periodic lane, `search/two_moon_periodic_890.py`)

It does not cover it. Reasons, in order of how firmly I can state them.

1. Wrong parameter (READ for the paper, READ in the module docstring for `#890`). The paper continues the particle's
   mass. The module continues "a moon-mass scale" `lam` that multiplies both moon GMs, with the planet GM adjusted so the
   total stays fixed, and the spacecraft is a massless test particle at every `lam`.
2. The `#890` model has no recoil (READ, module docstring). The moons sit on prescribed circles (Titania as the second
   primary of the rotating frame, Oberon on a concentric circle, no moon-moon force), so the primaries do not respond to
   the spacecraft at all. The paper's central object, the motion r(t) of the primaries driven by the particle, does not
   exist in that model.
3. The generating orbit is different (INFERRED). At `lam = 0` the `#890` orbit is a Kepler ellipse about the planet. In
   the frame of the moons it is a periodic orbit of a time-periodic system, but a Kepler orbit is a degenerate generating
   orbit: all orbits of the two-body problem are periodic, so the monodromy matrix is far from having a simple unit
   multiplier. Poincare's continuation then needs a selection of the generating orbit by the first-order perturbation
   (the project started at 1 percent of the moon masses, not at zero). Hadjidemetriou's hypothesis 3 is a nondegeneracy
   condition on a restricted-problem orbit that is already a member of a one-parameter family with varying period. That
   is not the situation at `lam = 0`.
4. The paper does not discuss degenerate or bifurcation generating orbits at all (READ absence), so the review digest's
   remark that a degenerate orbit "needs the second-species treatment" comes from the Gomez-Olle and Barrabes-Gomez
   papers, not from this one.

What would be the existence statement for the `#890` step, INFERRED: the system with prescribed circular moons is
time-periodic, so a periodic orbit is a fixed point of the stroboscopic map over the common period, and Poincare's
continuation in the parameter `lam` holds at any `lam` where the monodromy matrix does not have 1 as a multiplier. A
symmetric orbit stays symmetric. The practical check is the one the project already has: the corrector's reduced
Jacobian stays nonsingular along the continuation. The `#890` result was found numerically and was independently
re-integrated (see `data/OUTSTANDING.md`, `#890`), so no existence theorem was needed for it, but if one were wanted for a
write-up the references to follow are the continuation chapters of the two-body-seeded work (Casoliva et al. 2010 and the
Gomez-Olle digests) and not this paper.

Correction to carry back to the review digest: its sentence that this is "the existence theory behind the `#890`
two-moon lane" should be withdrawn. Its description of the theorem (continuation of a symmetric restricted orbit to one
where the small body has finite mass) is correct; its application to `#890` is not.

### 5.2 For the catalogue's restricted-problem cyclers in general

The statement is relevant to one thing only: a real spacecraft has nonzero mass. The theorem says that, to first order,
a nondegenerate symmetric periodic orbit of the circular planar restricted problem persists when the particle has small
mass, with a periodic motion of the primaries of size of order m2 (eq. 49). A 10-tonne spacecraft against the Earth-Moon
system is a mass ratio of order 10^-21 (COMPUTED from 1e4 kg over 6e24 kg), so the recoil correction is far below any
other modelling error, and for the project's purposes it does not matter.
The stability result, INFERRED to carry over: a hyperbolic cycler stays hyperbolic and a stable one stays stable at the
same order, with the multipliers unchanged to first order.

Two small notes. (a) The hypothesis w T not equal to 2 k pi excludes orbits whose rotating-frame period is an integer
number of primary revolutions. The reason is visible in eq. (45): the correction to the primaries' separation has the
small divisor omega^2 - omega0^2 n^2, which vanishes when T = 2 pi n / omega, so the primaries' epicyclic oscillation is
forced at its own frequency. This is the generic situation for an orbit that is periodic both in the rotating frame and in
the inertial frame (resonant with the primaries). The paper does not say whether the exclusion is a real obstruction or an
artefact of keeping the primaries on a circle and the period fixed (INFERRED; the natural remedy would be to let the
primaries' eccentricity vary). Near such periods the first-order correction gets large, so the "small mass" range shrinks.
This has no practical effect for a massless-limit spacecraft; it matters only if someone wants the finite-mass theorem for
the project's resonant orbits. (b) The paper is planar only and restricted to circular primaries, so it says nothing
about the project's three-dimensional, bicircular or elliptic lanes.

### 5.3 Is the symmetric-crossing argument the same one the project's symmetric shooting relies on?

The symmetry is the same one; the argument is not. READ: the symmetry is x2(t) = x2(-t), y2(t) = -y2(-t) (p.161), which is
the reflection (t, x, y) to (-t, x, -y) that the project's symmetric shooters use, written in `#890`'s docstring as
(tau, x, y, vx, vy) to (-tau, x, -y, -vx, vy). The paper's route, though, is Fourier parity: it shows that the first-order
solution built from even and odd forcing terms is itself even and odd (pp.162-167, Appendix p.173). It does not state or use
the other statement that the project's shooters depend on, namely that an orbit crossing the symmetry axis perpendicularly
twice is symmetric and periodic (the mirror theorem, standard in the texts the paper cites; not in this paper), and it does
not formulate a boundary-value problem from the crossing. So the project's shooting is a different, equivalent
formulation and the paper does not validate it (INFERRED). `#890`'s extra condition, that the perturber lies on the axis at
tau = T, comes from the periodic forcing and is not covered.

One related point in the paper's favour: section 8 makes no symmetry assumption. Its argument (nondegenerate monodromy
after removing the unit multiplier and the integral) supports continuation of asymmetric periodic orbits too, which is the
case the project's general corrector (`search/cr3bp_general_periodic.py`) addresses. INFERRED: the same removal of a phase
condition and one integral-redundant equation is done in practice by the project's fixed-Jacobi correctors. The condition
that the period varies along the family corresponds to a nonzero derivative of the period along the characteristic curve,
which is where a symmetric shooter that fixes the period loses its Jacobian, though I have not checked the project's
correctors for this.

## 6. Structural checks that could be made from the printed formulas

These involve no printed numerical orbit, only identities.

1. The primaries' radial block: det(A1(T) - I2) = 2 (1 - cos wT) (COMPUTED above).
2. For any restricted-problem periodic orbit, the 4 x 4 monodromy matrix Delta(T) must have a unit eigenvalue and
   det(Delta - I4) = 0 (p.170). The project's monodromy code should show a unit eigenvalue to rounding error.
3. det Delta*(T) not zero needs a family with varying period; at a turning point of the period along a family it vanishes.
4. First-order recoil of the primaries (eq. 49): the separation correction r - r0 = m2 phi(t) has period T and is
   bounded when 2 pi n / T differs from w; its size grows without bound as T approaches 2 pi k.

None of these is a golden value in the project's sense. They are consistency conditions on a project computation.

## 7. What is not held and may be worth obtaining

The companion paper on the stability of the continued orbits and the numerical families built from this theorem are cited
in the review digest and are not held (D): Hadjidemetriou, "The stability of periodic orbits in the three-body problem",
Celest. Mech. 12:255-276 (1975); Hadjidemetriou and Christides 1975; Bozis and Hadjidemetriou 1976; and Hadjidemetriou's
1973 Paper A (J. Mecanique 12:1) on the basic solutions of the variational equations. The paper's own reference list
(p.174) is Hadjidemetriou 1967, 1971, 1973; Siegel and Moser 1971; Standish 1969; Szebehely 1967, 1969, 1971; Szebehely
and Feagin 1973; Whittaker 1960.

## 8. Summary for the coordinator

- Which mass: the third body's (option (i)). The review digest's application of the theorem to `#890` is wrong; `#890` is
  option (ii) in a model without recoil.
- Method: first-order perturbation in m2 of the restricted orbit, with symmetry used through Fourier parity, plus a
  general implicit-function argument that removes the unit multiplier through the energy integral and a phase
  condition. Key hypotheses: a one-parameter family with varying period, and T not equal to 2 k pi.
- Results: continuation of one orbit and of a whole family at fixed small m2; the primaries stay in stable periodic
  motion; the particle's orbit is as stable or unstable as the original (first order, rotating frame). The unstable case
  is asserted, not proved. Larger m2 needs numerical continuation (stated by the author).
- No numerical example.
- For the project the paper is background for the mass of the spacecraft, which is irrelevant at the project's mass
  ratios. It is not the existence theory for moon-mass continuation.
