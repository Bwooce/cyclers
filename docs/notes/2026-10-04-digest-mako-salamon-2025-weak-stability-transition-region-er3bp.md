# Digest: Makó & Salamon (2025), "A note on the weak stability transition region in the planar elliptic restricted three-body problem"

Celestial Mechanics and Dynamical Astronomy 137:34 (2025), DOI 10.1007/s10569-025-10264-0, 21 pages (article pp1-19, appendix
pp16-19, references pp20-21). Zoltán Makó and Júlia Salamon, Department of Economic Sciences, Sapientia Hungarian University of
Transylvania, Miercurea Ciuc, Romania. Received 29 July 2025, revised 16 September 2025, accepted 14 October 2025, published online
30 October 2025. Open access (CC BY 4.0). Filed in the private paper corpus as
mako-salamon-2025-weak-stability-transition-region-planar-elliptic-restricted-three-body-problem-cmda-137-34-doi-10.1007-s10569-025-10264-0.pdf
(usable text layer: `pdftotext` gives about 8,870 words).

Digested 2026-10-04 from all 21 pages. Each statement is marked READ (seen on the page, page and section given), COMPUTED (our
arithmetic on printed numbers) or INFERRED (our reading). Page numbers are the journal's printed "Page n of 21". The paper contains
NO numbered tables; its numerical content is in the text of Sections 5 and 6, in Figures 2 to 8, and in the appendix equations. Figure
values are not transcribed beyond what the text prints.

## 0. What the paper is

READ (abstract, p1): "The paper investigates the weak stability transition region within the framework of planar elliptic restricted
three-body problem, considering the initial velocity as an independent variable. The boundary curves of the studied region are
defined, and a specific application to the Sun-Earth system is presented. Near the lower boundary of the transition region, the
influence of the true anomaly on the orbital stability is also analyzed."

In one paragraph: earlier weak-stability-boundary (WSB) work fixes the true anomaly of the primaries, the particle's initial
eccentricity about the secondary and the initial direction, and then makes the speed a function of position (eq. 1). This paper
instead scans the speed on a grid at each position. For each position the weakly stable speeds then form a finite or countable union
of open intervals (Proposition 1, a Cantor-type set), and the "transition region" is the band of positions and speeds in which stable
and unstable initial conditions alternate. Four boundary curves of that region are defined (Definition 5) and computed for the
Sun-Earth planar elliptic problem on a 1000 x 1000 grid, and the dependence on the true anomaly of Earth's orbit is shown by sweeping
f0 for two launch points (Figures 7 and 8). It is a short numerical note: no tables, no comparison to a published numerical value,
no capture-transfer design.

## 1. The model exactly as printed

### 1.1 Frame, variables, units (Section 2 p4; Appendix A.1, pp16-17; Fig. 9 p16)

READ. Planar elliptic restricted three-body problem (PER3BP). Massless P3 moves in the field of P1 (Sun, mass m1) and P2 (Earth, mass
m2) which move on Keplerian ellipses about their barycentre O with the same eccentricity e. P3 is restricted to the primaries'
orbital plane. Nondimensionalisation (p4): "normalized units such as the Earth's orbital period P = 2 pi and a(1-e) = 1, where a is
the semimajor axis of the Earth's orbit and e is its orbital eccentricity, then from the Kepler's third law we obtain G.m2 =
mu/(1-e)^3." The mass ratio mu = m2/(m1+m2) (INFERRED from the use of 1-mu for P1 in Fig. 9; the paper writes only "mu is the mass
ratio"). Independent variable: the true anomaly f of the primaries' relative orbit (p4),

    df/dt = sqrt(G(m1+m2)) / (a^(3/2) (1-e^2)^(3/2)) (1+e cos f)^2 = 1 / (1-e^2)^(3/2) (1+e cos f)^2

(Section 2, p4, unnumbered; the second form is in the units P = 2 pi, a = 1). Primary separation (A1, p17):
R = a(1-e^2)/(1+e cos f).

Rotating-pulsating frame O xi-tilde eta-tilde (a tilde marks the pulsating variables; Fig. 9 p16 shows P1 at (-mu, 0) on the left and
P2 at (1-mu, 0) on the right, so P1 is the Sun at -mu and P2 the Earth at 1-mu). Equations of motion (A2, p17), "Szebehely (1967)":

    xi-tilde'' - 2 eta-tilde' = d omega / d xi-tilde,    eta-tilde'' + 2 xi-tilde' = d omega / d eta-tilde,

with primes d/df, omega(f) = (1 + e cos f)^(-1) Omega(xi-tilde, eta-tilde) and

    Omega = (1/2)(xi-tilde^2 + eta-tilde^2) + (1-mu)/r1-tilde + mu/r2-tilde + (1/2) mu (1-mu).

"The planar circular restricted three-body problem (PCR3BP) can be considered as a special case of the elliptic problem, by setting
e = 0." READ. (INFERRED: this is the standard planar ER3BP, the same as `core/er3bp.py` restricted to z = 0, up to a constant in
Omega that does not affect the equations.)

### 1.2 Polar form about P2 (Appendix A.2 and A.3, pp17-18)

READ. Change of variables, with r2-tilde = |P2 P3| and theta-tilde the polar angle of P3 about P2 in the rotating-pulsating frame:

    xi-tilde = r2-tilde cos theta-tilde + 1 - mu,   eta-tilde = r2-tilde sin theta-tilde,
    xi-tilde' - eta-tilde = p_r2 cos theta-tilde - (p_theta / r2-tilde) sin theta-tilde,
    eta-tilde' + xi-tilde = p_r2 sin theta-tilde + (p_theta / r2-tilde) cos theta-tilde + 1 - mu.

Hamiltonian (p17-18), as printed:

    H = (1/2)(p_r2^2 + p_theta^2 / r2^2) - p_theta
        + 1/(1+e cos f) ( e cos f r2^2/2 - 1/r2 )
        - (1-mu)/(1+e cos f) ( r2 cos theta - 1/r2 + 1/sqrt(r2^2 + 2 r2 cos theta + 1) ),

(tildes dropped here for legibility; the paper carries them) and equations of motion (A3, p18):

    r2' = p_r2,      theta' = p_theta / r2^2 - 1,
    p_r2' = p_theta^2 / r2^3 - 1/(1+e cos f) ( e r2 cos f + 1/r2^2 )
            + (1-mu)/(1+e cos f) ( cos theta + 1/r2^2 - (r2 + cos theta) / (r2^2 + 2 r2 cos theta + 1)^(3/2) ),
    p_theta' = (1-mu) r sin theta / (1+e cos f) ( 1/(r2^2 + 2 r2 cos theta + 1)^(3/2) - 1 ).

COMPUTED check: the coefficient of 1/r2 in H is 1 rather than mu, but the (1-mu)/r2 term from the P1 bracket combines with it to give
-mu/r2, so the equations are consistent with a Newtonian attraction mu/r2 from P2. In the last line "r sin theta" is printed without
the subscript and tilde; this is a typesetting slip, r2-tilde and theta-tilde are meant.

### 1.3 Fixed (non-rotating) frame P2xy used to define capture (A4, p18)

READ. "To investigate the weak capture, we need a new, fixed (nonrotating) reference frame P2xy." Relations:

    theta = theta-tilde + f,      r2 = (1+e)/(1+e cos f) r2-tilde,
    theta-dot = 1/(1-e^2)^(3/2) (1+e cos f)^2 p_theta / r2-tilde^2   (printed with a second equal form using df/dt),
    r2-dot = (1+e)/(1-e^2)^(3/2) [ (1+e cos f) p_r2 + e r2-tilde sin f ],
    v2^2 = r2-dot^2 + r2^2 theta-dot^2
         = (1+e cos f)^2 / (2 (1+e)(1-e)^3) [ (p_r2 + D(f) r2-tilde)^2 + p_theta^2 / r2-tilde^2 ],   D(f) = e sin f / (1 + e cos f).

COMPUTED check of the v2^2 line: from the printed r2-dot and theta-dot, r2-dot^2 = (1+e cos f)^2 / ((1+e)(1-e)^3) (p_r2 + D r2-tilde)^2,
so v2^2 should carry the prefactor (1+e cos f)^2 / ((1+e)(1-e)^3) with no factor 2 in the denominator. The factor 2 as printed in the
v2^2 line is the one that belongs to E2 = v2^2/2 in (A5), which is internally consistent (COMPUTED: the mu term of A5 reduces to
G m2 / r2 with r2 = (1+e) r2-tilde/(1+e cos f)). So the v2^2 line appears to be the A5 prefactor carried over; use the component
equations and A5. (A typesetting slip, no consequence for the figures.)

Kepler energy about P2 (A5, p18):

    E2 = v2^2/2 - G m2 / r2
       = (1+e cos f)^2 / (2 (1+e)(1-e)^3) ( (p_r2 + D(f) r2-tilde)^2 + p_theta^2 / r2-tilde^2 - (1/(1+e cos f)) (2 mu / r2-tilde) ).

READ (p19): "The equations of motion (A3) and (A4), as well as the expression for the Kepler energy (A5) in the case of e = 0, coincide
with the formulas used by Topputo and Belbruno (2009) within the framework of the PCR3BP model."

Orthogonality (A6, A7, p19): cos(r2, v2) = (p_r2 - r2-tilde D(f)) / sqrt((p_r2 + D(f) r2-tilde)^2 ... ) as printed; r2 is perpendicular to v2
when `p_r2 = e r2-tilde sin f / (1 + e cos f)`. Direct motion in P2xy: `x y-dot - y x-dot > 0`, equivalent to `p_theta > 0`.

Hill sphere (A.6, p19): R_P2-bar(f) = R (m2/(3 m1))^(1/3) = a(1-e^2)/(1+e cos f) (mu/(3(1-mu)))^(1/3); its minimum is R_P2 =
a(1-e)(mu/(3(1-mu)))^(1/3), and in the normalised unit a(1-e) = 1, R_P2 = (mu/(3(1-mu)))^(1/3) (A8). READ.

### 1.4 System, numerical values (Section 5, p8)

READ: "For the Sun-Earth system, the mass ratio is mu = 3.003158242 * 10^-6 and the eccentricity of the elliptical orbit of Earth is
e = 0.0167." (First digits of the mass ratio are clear on the page; the figure `3.003158242e-6` is read from the printed page.)
Earth radius R_E = 6378 km (p8, Remark 1 text: "the distance between the massless particle and the Earth becomes equal to the Earth's
radius at some time, ||P3P2|| = R_E = 6378 km"). COMPUTED: with a = 149,597,870.7 km this gives a Hill radius a(1-e)(mu/(3(1-mu)))^(1/3)
= 1.4715e6 km, consistent with the R_P2 about 1.5e6 km marked in Figure 3. The model has no Moon: the Earth is a single point mass
(READ: "the Earth's orbit", "the Sun-Earth system"); see section 6 for a wording oddity regarding "the Moon's orbit".

## 2. The definitions

### 2.1 Two-body stability (Section 2, p3)

READ. Kepler energy E = v^2/2 - GM/r = h. h<0: bound, "both the starting point and orbit are considered stable"; h>0: hyperbolic,
unstable; h=0: parabolic, escape velocity v_e(r) = sqrt(2GM/r), "classified as unstable". Weak stability is "the notion of weak
stability (Fig. 1) introduced by Belbruno (1987)". "Weak capture (also known as ballistic capture)" is "the event where the Kepler
energy of the massless particle (e.g., a spacecraft or asteroid) relative to one of the primaries changes its sign from positive to
negative" (p1-2).

### 2.2 Definition 1: weak stability (p3-4)

READ, verbatim in the decisive part: "For a fixed value of the true anomaly f = f0, we consider the half-line l(f0, alpha) starting
from P2 [...] and making an angle alpha in [0, 360 deg) with the axis P1P2. Let P3 be the initial position of the massless particle,
taken at the periapsis of an osculating ellipse around P2, with semimajor axis on l(f0, alpha). The initial eccentricity e3 of P3 is
fixed. The initial velocity of P3 is perpendicular to the line l(f0, alpha). The initial distance between P2 and P3 in the fixed
coordinate system P2xy is r2 and the semimajor axis is a3 = r2/(1-e3). The modulus of the sidereal initial velocity relative to the
P2-centred reference frame is

    v2^2(r2, e3) = G m2 (1+e3)/r2 = mu (1+e3) / ((1-e)^3 r2) in [ v_c^2(r2), v_e^2(r2) ]            (1)

where mu is the mass ratio, v_c(r2) = v2(r2, 0) is the circular initial velocity and v_e(r2) = v2(r2, 1) is the escape initial velocity
with respect to primary P2 in the context of the two-body problem. The motion of a particle is said to be weakly stable relative to
P2, under three-body problem dynamics, if after leaving l(f0, alpha) it makes a full cycle around P2 without going near to P1 or
crashing into P2 and returns to a point on l(f0, alpha) with negative Kepler energy with respect to P2 in fixed P2xy system. Otherwise,
the motion will be weakly unstable."

So (READ): the section is the half-line l(f0, alpha), the cycle is "a full cycle around P2" returning to "a point on l(f0, alpha)", the
test is the sign of the Kepler energy about P2 at that return, and the failure modes are "near to P1" and "crashing into P2". "Near
to P1" is not quantified in the paper. The number of revolutions in this definition is exactly one (no n-revolution generalisation
appears in the paper; INFERRED from Definitions 1-3 and the code pseudocode of p10, which test one return).

### 2.3 Definition 2 (García-Gómez WSB, circular problem, pp4-5) and Definition 3 (the ER3BP extension, p5)

READ. Definition 2: for fixed (alpha, e3) in [0, 360 deg) x [0, 1), finitely many points 0 = r1* < r2* < ... < r2n* such that
`r2 in S*(alpha, e3) = union_{k=1..n} (r*_{2k-1}, r*_{2k})` gives a weakly stable orbit, otherwise weakly unstable, and the WSB is
`W_{G-G}(e3) = { (r2 cos alpha, r2 sin alpha) : alpha in [0, 360 deg) and r2 in dS*(alpha, e3) }`. Belbruno's original definition
assumed a single finite r2* = r2*(alpha, e3) > 0 (p4); García-Gómez "identified shortcomings in Belbruno's definition" and showed the
transition is a Cantor-type set (p4).

Definition 3 (p5): for fixed (f0, alpha, e3) in [0,360 deg) x [0,360 deg) x [0,1),

    S*(f0, alpha, e3) = { r2 > 0 : the orbit of a test particle with initial condition (r2, v2(r2, e3)) is weakly stable and the angular
                          velocity theta-dot(T) > 0 at the first return time T to half-line l(f0, alpha) },

and `W(f0, e3) = {(r2 cos alpha, r2 sin alpha) : alpha in [0,360 deg) and r2 in dS*(f0, alpha, e3)}`. The paper states that this
"extends García and Gómez's definition of WSB to the PER3BP model, taking into consideration the observations of Cecaroni et al.
(2012)" (the reference list spells the name Ceccaroni). The lower boundary of the WSB "starts at the initial point r1* where the first
weakly unstable orbit arises as we move away from P2".

INFERRED note on a possible slip: the condition theta-dot(T) > 0 at return in Definition 3 and in Proposition 1 selects counter-clockwise
(direct) return in P2xy. Yet the paper plots retrograde cases (Figures 5, 6; Remark 4) under the same set names. The proof of
Proposition 1 uses theta(f*) = 2 pi + alpha, which is likewise the direct case. A retrograde analogue would need theta(f*) = -2 pi +
alpha and theta-dot < 0; the text does not say how the retrograde sets were defined, and the pseudocode (p10) uses the generic code
values 1/0/-1 without the sign test.

### 2.4 Section 4: v2 as an independent variable (pp5-8)

READ. "Let f0 and e3 be given, and alpha, r2 be the variable parameters. In the procedure of constructing the WSB, first, we calculate
v2(r2, e3), the magnitude of the velocity and then we examine the weak stability of trajectories for initial values (r2, v2(r2, e3)).
In this case, the independent variable is r2. The question is: what are the properties of WSTR if v2 is an independent variable and
belongs to the interval [v2^min(r2), v2^max(r2)]?"

Lower bound on the velocity inside the Hill sphere (eq. 3, p5), "for r2 <= R_P2":

    v2^2 >= mu (1+e3) / (1-e)^3 * ( 3(1-mu)/mu )^(1/3),                                                   (3)

with the sentence "For a given r2, the smallest value of v2^2 is obtained by the formula (3)". INFERRED/COMPUTED: the right side has no
r2 in it; at e3 = 0 it equals v_c^2 evaluated at r2 = R_P2 (the Hill radius), so (3) is the circular speed at the Hill radius, a
lower bound over the whole Hill sphere, not a function of r2. Equation (3) is not used in the Section 5 computation, which scans
v2 in [0.9 v_c, 1.2 v_e] (p8).

Figure 2 (p6): initial conditions in the (r2, v2) plane for a fixed f0 and alpha; green = weakly stable, red = weakly unstable, yellow =
collision (r2 < R_E reached); two black curves are the 2BP circular (C) and escape (E) speeds; a vertical segment at a fixed r2 is
magnified and "contains multiple weakly stable and unstable subsegments". Figure 2 is a plot, no numbers other than axes (about 4 to
13 x 10^5 km on r2 and 0.5 to 1.2 km/s on v2).

**Proposition 1 (p6-7), verbatim statement:** "For fixed f0, alpha in [0, 360 deg) and distance r2 = ||P2P3|| the
S*(alpha, f0, r2) = { v2 > 0 : the orbit of test particle with initial condition (r2, v2), with r2 . v2 = 0 is weakly stable and the
angular velocity theta-dot(T) > 0 at the first return time T to half-line l(alpha, f0) } is an open set, and there exists a finite or
countably many points v*_k(r2) with v*_1(r2) = 0, where v*_i(r2) < v*_{i+1}(r2) for all i >= 1, such that
S*(alpha, f0, r2) = union_{k>=1} (v*_{2k-1}(r2), v*_{2k}(r2))." The proof (p6-7) takes a stable initial velocity, uses the implicit
function theorem on theta(f; v2) = 2 pi + alpha (angular component of A3) with theta-prime > 0, gets an open neighbourhood V of v2 on
which the return time g(v) is differentiable, notes the Kepler energy E2(g(v)) is continuous so V2 = {v in V : E2(g(v)) < 0} is open,
and writes an open subset of R as a countable union of open intervals. Closing remark (p7): "This also shows that S*(alpha, f0, r2)
has a Cantor-set-type structure." INFERRED: what is proved is openness (hence a countable union of open intervals); a Cantor
structure is the standard observation of the earlier García-Gómez work and is asserted rather than established here.

**Definition 4 (p7).** "There exist distances r2 from P2, such that: Property A: the maximum velocity for which the orbit remains
weakly stable is strictly greater than the minimum velocity at which the orbit becomes weakly unstable. Property B: weakly unstable
trajectories can occur when the initial velocity v2(r2) is lower than the circular velocity v_c(r2) = sqrt(mu/r2)."
(In the printed formula v_c = sqrt(mu/r2) the factor (1-e)^(-3/2) of eq. (1) is not shown; INFERRED that it is the normalised-unit
shorthand.)

**Definition 5 (p7), the four boundary curves in P2xy for fixed f0:**
1. LB(f0), the lower boundary: `{(r2* cos alpha, r2* sin alpha) : alpha in [0,360 deg) and r2* is the smallest distance from P2 on half-line l(f0, alpha) where Property A appears}`; it is "the bifurcation point r2* along the half-line from which the separatrix curve begins to transform into a bounded region".
2. CB(f0): smallest distance on l(f0, alpha) where Property B appears (weakly unstable orbits at v2 < v_c).
3. HB, the Hill circle around P2 of radius R_P2.
4. UB(f0), the upper boundary: the distance where Property A disappears, "from where the region begins to transform back into a separatrix curve".

**Definition 6 (p8), the weak stability transition region (WSTR):**

    WSTR(f0, alpha) = union over r2 in [r2^l, r2^u] of {r2} x [ v2*(r2), v2^max(r2) ],

where `(r2^l cos alpha, r2^l sin alpha) in LB(f0)`, `(r2^u cos alpha, r2^u sin alpha) in UB(f0)`, and `v2^max = sup S*(alpha, f0, r2)`.
READ. The lower velocity end v2*(r2) is the first transition point v*_1... as printed (the text reuses the star notation, v2*(r2)
being the first point where stability is lost; INFERRED).

READ, Remark 1 (p8): "The lower boundary of WSB corresponding to e3 = 0 is the smallest distance r2* from P2, where we obtain a weakly
unstable orbit when v2^2(r2, 0) = v_c^2(r2). In fact, CB(f0) estimates the lower bound of W(f0, 0)."

READ, Remark 2 (p8): decreasing R_E shrinks the yellow (colliding) zone, "the yellow region becomes point-like when R_E = 0. However,
this conclusion is not entirely rigorous and relies on semi-analytic and numerical results, as discussed by Belbruno (2024)."

## 3. What is computed and found

### 3.1 Setup (Section 5, pp8-10)

READ (p8). Sun-Earth, mu = 3.003158242e-6, e = 0.0167. "The distance r2, measured from the center of the Earth, is varied in the
interval [2 R_E, 1.5 R_P2]. This interval is divided into M = 1000 parts, which means that the step size is step_r2 = 14587 km. The
initial velocity v2 is varied in the interval [0.9 v_c, 1.2 v_e]. This interval is divided into N = 1000 parts with step_v2 = 0.0063
km/s = 22.6 km/h step." Stability is tested with `test_ode78`; on a stability change the boundary is refined by bisection with
sensitivity parameters eps_v2 = step_v2/10, eps_1 = eps_v2/2 (lower-boundary detection), eps_2 = step_r2/50 (p10). Integrator
(p10): "the seventh-eighth-order Runge-Kutta-Fehlberg integration scheme (ode78 in MATLAB) [...] with a relative error tolerance of
10^-10". The p10-11 pseudocode `Compute_S*` (v2 grid, classify code 1 weakly stable / 0 weakly unstable / -1 colliding with the Earth,
r2 < R_E; bisection on each code change) and `LB` (scan r2, call Compute_S* on consecutive r2 values, bisect on r2 where the first
boundary point moves) is reproduced in the paper and is a usable recipe.

COMPUTED discrepancy (stated as an observation, not a criticism): with R_P2 about 1.47e6 km (section 1.4) the interval [2R_E, 1.5R_P2]
spans about 2.19e6 km, so 1000 parts would be about 2195 km, not 14587 km; 14587 km x 1000 = 1.46e7 km would correspond to a
range about 6.6 times longer. Figure 3 (p9) plots r2 to about 2.5e6 km, which agrees with the shorter range. The printed step size
may therefore refer to a different grid than the interval quoted; it does not affect the qualitative results.

### 3.2 Figures 3 and 4: WSTR(0 deg, 45 deg) and WSTR(0 deg, alpha) (pp8-12)

READ. Figure 3 (p9) shows the Sun-Earth structure for f0 = 0 deg and alpha = 45 deg; first vertical line "the location of the lower
boundary point r2^l of the WSTR corresponding to the angle alpha = 45 deg of the boundary curve LB(0 deg). The property A appears for
the first time when r2* = r2^l"; r2^c "is the distance from where there exists initial velocity v2(r2^c) less than the circular
velocity v_c(r2^c) which leads to weakly unstable trajectory"; the last vertical line is R_P2 (the second panel is a magnified view
between r2^c and R_P2), and r2^u is marked beyond R_P2 (top-axis ticks 0.5, 1.5, 2, 2.5 x 10^6 km). INFERRED from the figure: r2^c
is about 0.8e6 km and R_P2 about 1.5e6 km; these are read off a plot and are not printed values.

READ. Figure 4 (p12) repeats WSTR(0 deg, alpha) for alpha = 0, 45, 90, 135, 180, 225, 270, 315 deg (direct initial velocities); the
dot marks the lower boundary point r2^l of LB(0 deg) in each panel. No numeric r2^l values are printed.

### 3.3 Figures 5 and 6: the three boundary curves in 3D and in the plane (pp11-14)

READ. Figure 5 (p13): 3D boundary curves for f0 in {0 deg, 180 deg}, alpha in [0,360 deg), direct and retrograde; third axis is the
lowest initial velocity at which the orbit is weakly unstable, v_I(r2*, alpha, f0); red = first boundary (LB), black = second (CB),
blue = fourth (UB). Also drawn is the circular velocity surface (CVS), eq. (6), p11: `x = r2 cos alpha, y = r2 sin alpha, z = sqrt(mu/r2)`.
Remark 3 (p11): "The red dots are located above the CVS, demonstrating that the lowest velocity v_I(r2, alpha, f0) at which Property
A is always between the circular velocity and the escape velocity. The black dots lie approximately on the CVS, which is expected
since the second boundary curve corresponds to Property B. The blue dots are situated below the CVS, demonstrating that the fourth
boundary curve is realized when the lowest initial velocity is less than the circular velocity. Consequently, at larger distances,
weakly unstable orbits can occur even with lower velocities."

READ. Figure 6 (p14): the boundary curves in P2xy, f0 in {0 deg, 180 deg}, direct and retrograde; LB red, CB black, UB blue; inner
circle is the Moon's orbit, outer circle the Hill circle. Remark 4 (p12): "in the case of direct motion, the dispersion for f0 =
180 deg is less than that for f0 = 0 deg. In the case of retrograde motion, the boundary curves become more diffuse (exhibiting
greater dispersion), the second boundary curve extends beyond the Hill circle, and the fourth boundary curve approaches it. In the
case of direct motion, the second boundary curve encloses a smaller area compared to the retrograde case. This reconfirms the
property that the relative frequency of the weakly stable points for direct motion is less than that of the retrograde motion. The
points of the LB(f0) are particularly significant for designing ballistic escape trajectories in the Sun-Earth system." (p12), and
"the lower boundary LB(f0) can be approximated by a Bernoulli lemniscate (Makó and Salamon 2014). According to Fig. 6, the farthest
lower boundary points from the Earth are near to the alpha = 45 deg" (p12-13).

### 3.4 Figures 7 and 8: dependence on the Earth's true anomaly f0 (pp13-15)

READ (p13). Figure 7: initial conditions "alpha = 45 deg, r2 = 381 829 km and v2 = 1.20446 km/s; is launched from a stable point of
WSTR for f0 in [0 deg, 360 deg) near to the first boundary curve; has a direct motion." Printed result: "If f0 in [0 deg, 177 deg] at
the starting moment then the orbit will be weakly stable, but if f0 in [178 deg, 215 deg] then the orbit will be weak unstable, and
for f0 in [216 deg, 360 deg) we get weak stable orbits again. This phenomenon is a consequence of the elliptical motion of the
Earth, because if we consider the Earth's orbit as a circular orbit, then we would obtain the weakly stable orbit as in first
picture of Fig. 7 for all f0 in [0 deg, 360 deg)." Figure 7 panels (p14) show f0 = 0, 90, 177, 178, 180, 215, 216, 270 deg
(trajectories in P2xy, star marks the start, Hill circle drawn).

READ (p13-15). Remark on the Moon's orbit: the Moon's semimajor axis is 384,400 km and its mean orbital velocity 1.022 km/s,
"Therefore, beyond the Moon swing, it has to be secured energy to obtain velocity about v_p = 0.20226 km/s = 202.26 m/s = 728.136
km/h to reach the escape velocity (order of magnitude of a passenger airplane speed)."
COMPUTED: 1.20446 - 1.022 = 0.18246 km/s, while 1.22426 - 1.022 = 0.20226 km/s exactly; and 202.26 m/s x 3.6 = 728.136 km/h as
printed. So the printed 0.20226 km/s is consistent with a start speed of 1.22426 km/s, not 1.20446 km/s: the two digit pairs 44 and
42 look transposed in one of the two places (INFERRED; the paper does not say which is intended).

READ (p15). Figure 8 initial conditions: "alpha = 135 deg, r2 = 12 754 km, and v2 = 7.63689 km/s; is launched from a stable point of
WSTR for f0 in [0 deg, 360 deg) near to the first boundary curve; has a direct motion." "The circular velocity for this distance is
v_c = 5.496891 km/s and the escape velocity is v_e = 7.773778 km/s." Result: "when the initial true anomaly f0 in [0 deg, 118 deg] at
the initial time, the resulting orbit is weakly stable. If f0 in [119 deg, 259 deg], the orbit will be weakly unstable, and for f0 in
[260 deg, 360 deg), we get weakly stable orbits again. [...] the weak instability of the orbits results from the fact that the Kepler
energy becoming positive and the test particle leaves the Hill circle. If f0 in [119 deg, 191 deg], [...] If f0 in [192 deg, 259 deg],
loops are formed on the orbits. Assuming circular motion of the Earth, all orbits remain weakly stable, as illustrated in the first
graph of Fig. 8." (The text gives the intervals as printed: "If f0 in [119, 191], the weak instability ... leaves the Hill circle"
and "If f0 in [192, 259], loops are formed", the exact sentence layout on p15 is split across the two statements; the intervals are
read correctly but the grouping of the explanation to each interval is INFERRED.) Figure 8 panels (p15) show f0 = 0, 90, 118,
119, 180, 259, 260, 270 deg.

COMPUTED check of the units: r2 = 12,754 km is 2 R_E to within the paper's own value 6377 km (the quoted R_E is 6378 km, giving
12,756 km), the lowest r2 on the grid. With GM_Earth = 398,600.44 km^3/s^2, sqrt(GM/r2) = 5.59044 km/s and sqrt(2 GM/r2) = 7.90601
km/s. The paper's printed v_c = 5.496891 and v_e = 7.773778 km/s are smaller by the factor 0.98327 (v_c) and 0.98327 (v_e), which equals
1 - e = 0.9833 for e = 0.0167. So the km/s speeds in the paper are the standard two-body speeds multiplied by (1 - e), to about
3e-5 relative (5.49708 computed versus 5.496891 printed). INFERRED: this is the velocity scale implied by the normalisation a(1-e) = 1
with period 2 pi, and it is not stated anywhere in the paper. Any reproduction in km/s has to apply this factor to match the printed
v2 values; this also means "km/s" in this paper are not the heliocentric-inertial km/s of the Earth-relative velocity.
The same factor applied to the Figure 7 start gives v_c = 1.00466 km/s at 381,829 km (COMPUTED), so 1.20446 (or 1.22426) is between
v_c and v_e = 1.42081 km/s there (both COMPUTED with the (1-e) factor).

Findings and the stated interpretation (Section 6, p16): the weak stability of orbits near LB depends on the Earth's true anomaly,
"in line with the conclusions of Hyeraci and Topputo (2013)" (cited as Hyeraci and Topputo, Celest. Mech. Dyn. Astron. 116, 175-193);
the shape of the first boundary curve "raises the question of whether it is possible to express this curve analytically in terms of
ER3BP or CR3BP". A design claim (p15-16): "weakly unstable points are located very close to the Earth, when the angle between the
Sun-Earth direction and the Earth-test particle direction is around 135 deg. [...] Starting from these points, escape trajectories
can be planned even if the initial velocity is lower than the escape velocity (see Fig. 8)."

### 3.5 Comparison with earlier definitions, as the paper itself frames it

READ (Sections 1 and 3, pp2-5). Belbruno (1987, 2004): a single finite r2*(alpha, e3), speed from eq. (1), PCR3BP. García and Gómez
(2007): Cantor-set transition, Definition 2, still circular and speed from eq. (1). Cecaroni (Ceccaroni) et al. (2012): analytic
definition of the WSB and the topology of the weakly stable set in the PCR3BP. Belbruno (2024): the WSB is a union of infinitely many
Cantor sets of hyperbolic points when the particle makes infinitely many cycles about the secondary. Hyeraci and Topputo (2013), Makó
(2014): true-anomaly dependence in the ER3BP. Érdi et al. (2009): stable regions about L4 in the ER3BP. The paper's stated
difference (p2): earlier work treats f0, e3 and the direction as constants and "the independent variable is r2, and the dependent
variable is v2"; here "v2 appears as an independent variable". Applications cited: Sousa Silva and Terra (2012); Vetrisano et al.
(2015); Topputo (2013); Campana and Topputo (2025); Qi and Xu (2014); Topputo and Belbruno (2015); Oshima et al. (2017);
Sousa-Silva et al. (2018); Dutt et al. (2018); Wang et al. (2025). No numerical comparison with any of them is made.

## 4. Method, as far as needed to read the numbers

READ/INFERRED. For each (r2, alpha, f0): integrate the planar ER3BP in the pulsating-rotating polar form A3 from the section
l(f0, alpha), initial state periapsis-like (r2 perpendicular to v2, A7), for 1000 speeds between 0.9 v_c and 1.2 v_e; code each speed as
weakly stable, weakly unstable or colliding (r2 < R_E), refine every code change by bisection to eps_v2. The WSTR boundary points are
then found with a second nested bisection in r2. Integration tolerance 1e-10 relative. Collision radius R_E = 6378 km; "near P1"
undefined. The "full cycle" is detected through theta in the fixed frame returning to 2 pi + alpha (proof of Proposition 1).

## 5. Positive controls for the project

READ and COMPUTED: there is no printed boundary point of the form "mass ratio, eccentricity, true anomaly, direction, boundary r2 and
v2" in the paper. The following are the only printed numbers that `core/er3bp.py` plus a small capture predicate could be tested
against; none has been run (task instruction: no computations).

| Candidate | Printed inputs | Printed outcome | Page |
|---|---|---|---|
| Fig. 7 sweep | mu = 3.003158242e-6, e = 0.0167, alpha = 45 deg, r2 = 381 829 km, v2 = 1.20446 km/s (see slip in section 3.4), direct, f0 integer degrees | weakly stable for f0 = 0 to 177 deg, weakly unstable for 178 to 215 deg, weakly stable for 216 to 359 deg | 13 |
| Fig. 8 sweep | same mu and e, alpha = 135 deg, r2 = 12 754 km, v2 = 7.63689 km/s, direct | weakly stable for f0 = 0 to 118 deg, unstable 119 to 259 deg, stable 260 to 359 deg | 15 |
| Circular control | same, e = 0 | weakly stable for every f0 (statement, Fig. 7 and 8 first panels) | 13, 15 |
| Unit check | r2 = 12 754 km | v_c = 5.496891 km/s, v_e = 7.773778 km/s (equal to the two-body values times 1 - e, section 3.4) | 15 |

Recipe: (1) set mu, e as above; (2) choose f0; place P3 at distance r2 from P2 in the direction angle alpha measured from the P1P2 axis
in the P2xy frame at that f0, velocity perpendicular to the radius, counter-clockwise (direct), speed v2 converted with the (1 - e)
factor to the paper's km/s scale (or equivalently scale the printed km/s up by 1/(1-e) to the standard two-body scale); (3) convert to
the pulsating-rotating coordinates through A4; (4) integrate `propagate_er3bp` in f from f0 forward until theta in P2xy has advanced by
2 pi (return to the half-line), stopping if r2 < R_E (colliding) or P3 gets "near P1"; (5) classify weakly stable if E2 (A5) < 0 at
the return. The step-like threshold at integer degrees (177/178, 215/216, 118/119, 259/260) makes these good bracketing tests: they
are crisp, not a fit of a smooth curve.

Caveats that must be settled before using this as a control: (a) "near P1" is undefined and the collision radius R_E = 6378 km versus
the 6377 km implied by r2 = 12,754 km; (b) the printed 1.20446 versus 1.22426 ambiguity; (c) the origin of the (1 - e) velocity scale
is inferred, not stated (if the wrong scale is used the classification will fail and the cause would be the scale, not the code);
(d) the paper's sign for "direct" (p_theta > 0 in the fixed frame, A.5); (e) whether f0 is the true anomaly of the Earth about the
Sun at the section crossing, with f0 = 0 at perihelion (INFERRED: yes, Fig. 9 and A1). Also: (f) the project's `core/er3bp.py` has
an identity test only against its own physics (see `#893`: item (a) not yet done for this module); a boundary point reproduced from
this paper would be a (b)-type control of the independent-published-numbers kind, though a weak one because only a classification
(not a trajectory) is printed. Per the project's rule "no circular goldens", the expected values here are the printed intervals,
not anything our code computed.

## 6. What this means for the project

**Scope (stated plainly).** Nothing in the paper falls in the catalogue's classes (cycler, quasi_cycler, precursor_mga, mga_tour).
It contains no periodic orbit, no transfer, no flyby sequence and no delta-v. It is a definitional and numerical-structure note on the
weak stability boundary. Its only possible bearing is on the capture-substrate work (`#378`, `#681`).

**Ballistic capture and weak stability work.** The project's code already uses a related but different definition in three places.

| Item | Mako-Salamon 2025 | Project code |
|---|---|---|
| Model | planar ER3BP, Sun-Earth, polar pulsating-rotating form | `core/wsb.py`: Earth-Moon BCR4BP (`core/bcr4bp.py`, the model corrected under `#891`); `core/sunmars_wsb.py`: planar Sun-Mars, Mars on a Keplerian ellipse (e = 0.093418671), heliocentric inertial frame, physical units |
| Section | half-line l(f0, alpha) from P2, particle at periapsis of an osculating ellipse | `wsb.py`: periapsis surface `r-dot_23 = 0` (Belbruno eq 3.9), W given by (r_2, theta_2, e_2); `sunmars_wsb.py`: Mars-relative periapsis predicate |
| Speed | independent variable; WSTR is a region in (r2, v2); e3 dropped | dependent: speed from the osculating eccentricity e_2 (Belbruno eq 3.6/3.29 analytic W; `capture_periapsis_state` in `sunmars_wsb.py`), as in earlier work |
| Revolutions | one full cycle about P2, return to the half-line, energy sign test | `wsb.stability_class`: "propagate one Moon-revolution from a periapsis state" with osculating-period horizon, five labels (stable, unstable, capture, escape, primary_interchange); `n_rev` argument allows more. `sunmars_wsb.py` follows Topputo-Belbruno 2015 n-stability (digest `2026-07-22-digest-topputo-belbruno-2015.md`) |
| Direction | direct only in the definitions (theta-dot > 0); retrograde shown | both branches |
| True anomaly | f0 is an explicit parameter and the main variable of Figs 7-8 | Earth-Moon BCR4BP: Sun phase, no Earth eccentricity; Sun-Mars: f0 is a parameter of `mars_state` and `sunmars_eom` |
| "Near P1", collisions | collision r2 < R_E coded -1; "near P1" unquantified | `escape_radius = 1.5` LD and Earth-relative distance used for primary interchange (wsb.py) |

COMPUTED/INFERRED comparison: the paper's Definition 1/3 is the single-revolution case of the project's `stability_class` with the
same energy-sign test, so the two agree on what "weakly stable" means; they differ in (i) the independent variable (the paper scans
speed, the project fixes e_2 and derives speed), (ii) the model (ER3BP versus BCR4BP and the project's Sun-Mars ER3BP-equivalent),
(iii) the return section (half-line at fixed alpha versus the periapsis condition). Consequence for `#378` and `#681`: the project's
searches sample the (r_2, theta_2, e_2) periapsis surface. The paper shows that at fixed position the stable speeds form disjoint
intervals and that stability can occur above and below the circular speed (Property B: weakly unstable below v_c; Properties A and B
together imply the stable set is not a simple shell between v_c and v_e). The e_2-parametrised surface therefore samples only the
band v_c to v_e at each position (e_2 in [0,1) maps to v in [v_c, v_e)) and misses any weakly stable or unstable structure at lower
speed; it would be insufficient to claim that no weakly stable orbit exists at lower speeds. This does not invalidate the clean
negatives of `#378` and `#681`, whose object is a repeating capture-escape chain, but it is a stated limit of the sampling: "empty"
is conditional on the e_2 parametrisation (see the negative-results registry rule). INFERRED, not tested.

The paper's Sun-Earth case is the Earth analogue of what `#681` did for Mars (strong sensitivity of weak stability to the primary's
true anomaly, shown in Figs 7-8); it supports keeping f0 a swept parameter in any further sweeps. It does not alter the project's
catalogue; no catalogue row cites it.

Relevance to the corrected bicircular code: this paper's model is the ER3BP, not the BCR4BP, so the `#891` defect has no bearing on it.
It could, however, serve as an independent published check on `core/er3bp.py` itself (`#893` item (b)), subject to the caveats in
section 5.

## 7. References cited that the project might want

Corpus status checked against `docs/notes/CORPUS_INDEX.md` on 2026-10-04 (index rows grepped by author and title words). "Held"
means a row exists; "not held" means no row.

| Reference | Held? |
|---|---|
| Belbruno, E.: Lunar capture orbits, a method for constructing. AIAA 87-1054, Proc. Int. Propul. Conf. (1987) | not held |
| Belbruno, E.: Capture Dynamics and Chaotic Motions in Celestial Mechanics. Princeton University Press (2004) | held (digest 2026-06-17-digest-belbruno-2004.md) |
| Belbruno, E.: Cantor set structure of the weak stability boundary for infinitely many cycles in the restricted three-body problem. Celest. Mech. Dyn. Astron. 136, 53 (2024) | not held; this is the most relevant for the definition (union of Cantor sets) |
| Belbruno, E., Topputo, F., Gidea, M.: Resonance transitions associated to weak capture in the restricted three-body problem. Adv. Space Res. 42, 1330-1351 (2008) | not held |
| Belbruno, E., Gidea, M., Topputo, F.: Weak stability boundary and invariant manifolds. SIAM J. Appl. Dyn. Syst. 3 (sic: 9), 1061-1089 (2010) | not held (the paper's reference list prints volume 3) |
| Belbruno, E., Gidea, M., Topputo, F.: Geometry of weak stability boundaries. Qual. Theory Dyn. Syst. 12, 53-66 (2013) | not held |
| Campana, F., Topputo, F.: Ephemeris refinement of low-energy Earth-Moon transfers via an astrodynamic model chain. Acta Astronaut. 232, 271-282 (2025) | not held |
| Ceccaroni, M., Biggs, J., Biasco, L.: Analytic estimates and topological properties of the weak stability boundary. Celest. Mech. Dyn. Astron. 114, 1-24 (2012) | not held; the analytic definition that Definition 3 extends |
| Dutt, P., Anilkumar, A.K., George, R.K.: Design and analysis of Weak Stability Boundary trajectories to Moon. Astrophys. Space Sci. 363, 161-172 (2018) | not held |
| Érdi, B., Forgács-Dajka, E., Nagy, I., Rajnai, R.: A parametric study of stability and resonances around L4 in the elliptic restricted three-body problem. Celest. Mech. Dyn. Astron. 104, 145-158 (2009) | not held |
| García, F., Gómez, G.: A note on weak stability boundaries. Celest. Mech. Dyn. Astron. 97, 87-100 (2007) | not held; defines Definition 2 |
| Hamilton, D.P., Burns, J.A.: Orbital stability zones about asteroids. Icarus 96, 43-64 (1992) | not held |
| Hyeraci, N., Topputo, F.: The role of true anomaly in ballistic capture. Celest. Mech. Dyn. Astron. 116, 175-193 (2013) | not held; the closest earlier result on true-anomaly dependence |
| Makó, Z.: Connection between Hill stability and weak stability in the elliptic restricted three-body problem. Celest. Mech. Dyn. Astron. 120, 233-248 (2014) | not held |
| Makó, Z., Salamon, J.: Weak stability transition region near the orbit of the Moon. Multi-scale (time and mass) dynamics of space objects, Proc. IAU Symp. 364, Vol 15, 184-190, Cambridge Univ. Press (2021) | not held |
| Oshima, K., Topputo, F., Campagnola, S., Yanao, T.: Analysis of medium-energy transfers to the Moon. Celest. Mech. Dyn. Astron. 127, 285-300 (2017) | not held (the corpus holds other Oshima papers, not this one) |
| Qi, Y., Xu, S.: Lunar capture in the planar restricted three-body problem. Celest. Mech. Dyn. Astron. 120, 401-422 (2014) | not held |
| Romagnoli, D., Circi, C.: Earth-Moon weak stability boundaries in the restricted three and four body problem. Celest. Mech. Dyn. Astron. 103, 79-103 (2009) | not held |
| Sousa Silva, P.A., Terra, M.O.: Applicability and dynamical characterization of the associated sets of the algorithmic weak stability boundary in the lunar sphere of influence. Celest. Mech. Dyn. Astron. 113, 141-168 (2012) | not held |
| Sousa-Silva, P.A., Terra, M.O., Ceriotti, M.: Fast Earth-Moon transfers with ballistic capture. Astrophys. Space Sci. 363, 210-221 (2018) | not held |
| Szebehely, V.: Theory of Orbits. Academic Press, New York (1967) | held (digested, per the index text on the Szebehely gap) |
| Topputo, F., Belbruno, E.: Computation of weak stability boundaries: Sun-Jupiter system. Celest. Mech. Dyn. Astron. 105, 3-17 (2009) | not held; the paper says its A3-A5 reduce to this at e = 0 |
| Topputo, F.: On optimal two-impulse Earth-Moon transfers in a four-body model. Celest. Mech. Dyn. Astron. 117, 279-313 (2013) | not held |
| Topputo, F., Belbruno, E.: Earth-Mars transfers with ballistic capture. Celest. Mech. Dyn. Astron. 121, 329-346 (2015) | held (digest 2026-07-22-digest-topputo-belbruno-2015.md) |
| Vetrisano, M., Van der Weg, W., Vasile, M.: Navigating to the Moon along low-energy transfers. Celest. Mech. Dyn. Astron. 114, 25-53 (2012) | not held (the text cites it as Vetrisano et al. 2015; the reference list gives 2012) |
| Wang, M., Yin, W., Zhang, C., Zhang, H.: Mechanism and characteristics analysis of weak stability boundary transfers to the 2:1 distant retrograde orbit. Adv. Space Res. 75, 2221-2250 (2025) | not held |

Priority for acquisition if the capture work is extended: Hyeraci and Topputo 2013 (the true-anomaly dependence), Ceccaroni et al.
2012 and Belbruno 2024 (the definition), García and Gómez 2007 (the Cantor structure), and Topputo and Belbruno 2009 (the PCR3BP
polar equations A3-A5 reduce to).

## 8. Minor observations (factual, for completeness)

- Citation year slips in the text: "Cecaroni" for Ceccaroni (p1, p5); "Belbruno et al. (2010)" in text and reference list differ in the
  SIAM volume number (3 in the list; the journal volume is 9); "Vetrisano et al. (2015)" in the text and 2012 in the list. These do
  not affect the results.
- Definition 4 prints v_c = sqrt(mu/r2) and the Section 5 v_c is a factor (1 - e) below the two-body value; both can be reconciled by
  the normalisation of section 1.1 but the paper does not spell it out.
- Section 6 says the first boundary curve "lies close to the Moon's orbit" and the Figure 6 caption draws "the orbit of the Moon", yet
  the model has no Moon; the Moon's orbit is a length marker only (semimajor axis 384,400 km).
- The data availability statement reads "No datasets were generated or analyzed during the current study." There are no tables and no
  code link.

## Addendum 2026-10-04 (`#896`): corrections found when the printed numbers were made tests

Test: `tests/core/test_er3bp_mako_salamon_2025.py` (commits `b752b68b`, `fd2eeabf`); the
non-reproduction of the unstable bands is registered as `#925`. Points of this digest found against
the PDF:

- Section 1.3 copies A7 (p19), p_r2 = + e r2~ sin f / (1 + e cos f), without comment. A4 (p18) gives
  r2-dot proportional to (1 + e cos f) p_r2 + e r2~ sin f, which vanishes only for the opposite sign,
  so A7 and the A6 numerator have a sign slip. Using A7 as printed shifts the boundary speed by less
  than 5e-4.
- Section 4 presents only the fixed return half-line; the proof of Proposition 1 (p6) reads as a
  rotating one (theta-dot = (df/dt)(theta' + 1) holds for the rotating-pulsating angle). Both readings
  give the same classifications.
- Not in the digest: p3 speaks of "threshold values of the mean anomaly" while the sweeps are in the
  true anomaly f0 (differing by up to 1.9 deg at e = 0.0167).
- Sections 3.4 and 5: the 1.22426 km/s variant for Figure 7 escapes at every f0 in the model; the
  printed 1.20446 km/s sits on the model's boundary, so the p15 v_p remark is the likelier slip.
- Not in the digest: at the printed Figure 8 speed (0.982 v_e) the two-body apoapsis is about
  3.5e5 km, but the plotted orbits reach several million km, which the model gives only near
  0.998 v_e; the printed speed is probably not the one Figure 8 was computed with.
- Section 1.4 attributes R_E = 6378 km to Remark 1; it is in the main text between Remarks 1 and 2
  (p8).
