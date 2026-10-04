# Digest: Gomez and Olle 1986, "A note on the elliptic restricted three-body problem" (mu = 0)

Date 2026-10-04. Purpose: the mu = 0 groundwork of the elliptic second-species theory: symmetric periodic orbits of the rotating
pulsating two-body problem, arc solutions (orbits with consecutive collisions, OCC) and double-collision orbits of the planar
elliptic problem as the primaries' eccentricity e_p runs from 0 to 1. It is the elliptic generalisation of Henon 1968 (digest
`docs/notes/2026-10-04-digest-henon-1968-consecutive-collision-orbits.md`) and the base of the 1991 papers (digest
`docs/notes/2026-10-04-digest-gomez-olle-1991-second-species-circular-elliptic-I-II.md`, where it is cited as "[2]" in Part II and
"[1]" in Part I). Written for #899 (an elliptic enumerator alongside the Henon one), #912, #925, #930 and #906.

Citation: G. Gomez and M. Olle, "A note on the elliptic restricted three-body problem", Celestial Mechanics 39:33-55 (1986), DOI
10.1007/BF01232287. 23 pages (journal pp. 33-55; PDF page n is journal page n+32). Received 20 May 1985, accepted 27 June 1986.
Filed in the private paper corpus as
gomez-olle-1986-note-elliptic-restricted-three-body-problem-mu-0-arcs-celest-mech-39-33-doi-10.1007-BF01232287.pdf.
Text layer: `pdftotext <file> - | wc -w` gives 6879 words (usable for prose; every equation and figure label below was read from the page
images, all 23 pages).

Evidence labels: READ = from the printed page; COMPUTED = my own arithmetic on 2026-10-04 (double precision, scratch scripts not kept);
INFERRED = my reading of what the print implies. Journal page numbers. The paper prints NO numerical tables: its numbers are formulas,
a handful of values in the text and figures, and plotted curves. Nothing below was invented; plot readings are marked approximate.

## 0. Summary

1. Part 1 (pp. 34-38). Symmetric periodic orbits of the rotating pulsating two-body problem (a mu = 0 problem: the secondary is absent
   from the dynamics): they are Kepler ellipses with rational mean motion n = q/p; the synodic period is 2 pi p; the perpendicular axis
   crossings are exactly the passages through the pericentre or apocentre (Proposition 1); the admissible initial epochs are listed
   (Proposition 2); there are no circular orbits except the point and the unit circle with n = +-1 and e = e_p (Proposition 3, Corollary:
   the elliptic problem at mu = 0 has no circular orbits in the synodic system). Characteristic curves are given in the (I0, x) plane,
   with I0 = h - c (sidereal energy minus angular momentum): eq. (7) for circles, eq. (10) for ellipses (themselves ellipses).
2. Part 2 (pp. 38-49). Orbits with consecutive collisions: the single implicit equation (17) in (tau, eta) with the parabolic (19) and
   hyperbolic (20) variants, the 8 sign configurations (eps, eps', eps'', eps_p), the families A_i, B_i, C_ij following Henon, their
   diagrams at e_p = 0.5 and 0.98 (Figs. 6-9), the existence condition (21) of the double-point families C_ij, and orbit shapes (Figs. 10-13).
   For e_p = 0 and eps_p = +1 the equations reduce to Henon's (p. 42).
3. Part 3 (pp. 49-55). Double-collision orbits (rectilinear sidereal orbits colliding with both primaries): equations (22)-(27), the
   circular-case theorem (Pinol), a numerical Proposition (p. 54) giving the number of double-collision orbits as a function of the semi-major
   axis a, the sign pair (eps_p, eps'') and four tangency values a^{++}, a^{+-}, a^{-+}, a^{--} that are plotted against e_p (Fig. 20).
4. What it is not: no mu > 0 content, no corrector, no stability, no tables of orbits, no existence proof of the numerically found families.
   Scope statement (p. 34): "At present, we are studying the extension, for mu /= 0, of second species solutions of the elliptic restricted
   problem, which is, in fact, the interesting case for practical purposes (i.e. the rendez-vous problem) in Space Dynamics." (That is the
   1991 pair.)
5. My verification (section 6): the equations reproduce Henon's parabolic orbit at e_p = 0 (tau/pi = 0.163926, matching the value
   0.16393 in Henon's Table 2 first row); Eq. (21)'s printed numbers are right (1.539, 8, 16, 4.617, 32, 7.695, 1.015); the first
   double-point families at e_p = 0.98 (C67,68 and C68,69) are right; the four printed tangency values of Figs. 16-19 are reproduced by
   solving system (26) (0.3448, 0.7105, 0.5847, 0.4636) PROVIDED (26b) carries the factor eps_p that the print omits. Five printed defects
   found (section 7), one of which (the parabolic tau/pi at e_p = 0.5) differs from my solution in the third decimal.

## 1. Setting and the rotating pulsating two-body problem (pp. 34-38, READ)

Frame (p. 34): "In the whole paper, we assume that one primary is located at the origin and the other one describes a direct ellipse with
mean angular velocity equal to unity." The infinitesimal body is P3; the primary at the origin is P1; the other primary P2 does not interact
with P3 in Part 1 ("assuming that the other primary P2 is entirely absent, in the sense that it does not interact with the third body P3",
p. 34; mu = 0).

Eq. (1) (p. 34), rotating pulsating (synodic) system, prime = d/df, f the true anomaly of the primaries' orbit, e their eccentricity:

x'' - 2 y' = x/(1 + e cos f) (1 - 1/r^3),   y'' + 2 x' = y/(1 + e cos f) (1 - 1/r^3),   r^2 = x^2 + y^2.

(READ; COMPUTED comparison with the code: `core.er3bp.er3bp_eom` at mu = 0 and z = 0 reads xdd = 2 y' + (x - x/r^3)/(1 + e cos f),
ydd = -2 x' + (y - y/r^3)/(1 + e cos f), which is eq. (1) exactly; the project's frame has the big primary at (-mu, 0) = origin at mu = 0.)

Symmetry (p. 34): (x(f), y(f), x'(f), y'(f)) -> (x(-f), -y(-f), -x'(-f), y'(-f)) leaves the equations invariant, "so a periodic orbit will be
symmetric if it intersects twice the synodical axis x perpendicularly."

Proposition 1 (p. 34): "The points where a rotated pulsated ellipse intersects perpendicularly the synodical axis x, correspond to passages of
P3 by the pericenter or the apocenter of its side real orbit." Proof: perpendicular crossing means dr/df = 0 for r = r_sid/r_p, and P2 must be at
pericentre or apocentre of its own orbit for the synodic orbit to be symmetric (p. 35).

Proposition 2 (p. 35): let M = n (t - t0) and M_p = t be the mean anomalies of P3 and P2, n = q/p rational. If the initial epoch satisfies
(i) n t0 = j pi, j in Z; or (ii) n t0 = n pi, 2 n pi, ..., (p-1) n pi; or (iii) n t0 = (n+1) pi, (2n+1) pi, ..., ((p-1) n + 1) pi, then the
synodic periodic orbit is symmetric. Proof (i): first orthogonal crossing at t = t0, second at t = T_syn/2 = pi p; at that instant E_p = pi p,
f_p = pi p, E = (q - j) pi. In (ii) the first crossing is the first passage of P3 through its pericentre; in (iii) through its apocentre.

Proposition 3 (p. 35): for t0 = 0 and e_p /= 0, an orbit of the rotating pulsating two-body problem is a point or a circle iff n = +-1 and
e = e_p (e the sidereal eccentricity of P3). The proof (p. 36) splits n = q/p into odd/even parity cases via eq. (2) 1 - e cos E =
K-bar (1 - e_p cos E_p) and the values (3) E_p = 0, pi p; E = 0, pi q. [The proof's cases (ii) and (iii) are both headed "q even, p odd" and
"q odd, p odd"; the second duplicates case (i), so one of the three case labels is a misprint (INFERRED); the conclusion is unaffected.]
Corollary (p. 36): "The plane elliptic restricted three body problem for mu = 0, has no circular orbits in the synodical system."
The critical points at mu = 0 of the vector field (1) are (cos theta, sin theta, 0, 0), theta in R (mod 2 pi) (p. 36).

### 1.1 Characteristic curves (pp. 36-38, READ)

I0 = h - c (eq. 4): sidereal energy h minus angular momentum c; x the abscissa of a perpendicular crossing. For a sidereal orbit of semi-major
axis a and semi-minor b: +-b = I0 sqrt(a) + 1/(2 sqrt(a)) (eq. 5; from c = +-a sqrt(1-e^2), h = -1/(2a); the sign is direct/retrograde).
(i) Circles (b = a): |x| = a/r_p (eq. 6) with r_p = 1 + e_p or 1 - e_p; for x > 0, I0 = +- sqrt(x r_p) - 1/(2 x r_p) (eq. 7), "four equations for a
fixed e_p" (two per value of r_p); Fig. 1 plots them at e_p = 0.3 (curves cross I0 = -3/2 and the x-axis at 1/r_p). Crossings with x < 0 take -x
in eq. (7). Only a dense subset is periodic: |n| = a^(-3/2) = q/p and T_syn = 2 pi p (cites Wintner [20]).
(ii) Ellipses: with n = a^(-3/2) rational, eq. (8): (I0 + 1/(2a))^2/a + e^2 = 1; the crossings (eq. 9) |x_1| = a(1-e)/r_p, |x_2| = a(1+e)/r_p,
r_p = 1 +- e_p; eliminating e (eq. 10): (I0 + 1/(2a))^2/a + (|x| - a/r_p)^2/(a^2/r_p^2) = 1, "ellipses centered at (-1/2a, a/r_p) and with
semiaxis sqrt(a), a/r_p". Fig. 2: both curves (r_p = 1 +- e_p) for x > 0, e_p = 0.3, orbit labels (1)-(4) with different a (or N = |n|) [values of a
for the labels not printed].

## 2. Orbits with consecutive collisions (pp. 38-49)

### 2.1 Definition and setting (pp. 38-40, READ)

P2 has mass zero in the dynamics but collisions P3-P2 are possible (the point x = 1, y = 0 of the synodic orbit). "Such points of collision break
a solution of the problem of two bodies into pieces, and each such piece is an independent solution of the elliptic restricted problem of three
bodies for mu = 0." Pieces that begin and end at a collision are Poincare's "orbits with consecutive collisions". "Our results are an extension
of Henon's paper [6], where the primaries describe circular orbits, while we consider the eccentricity of the elliptical orbit of the primaries."
Method: the sidereal plane; P3 elliptic, hyperbolic or parabolic.

Elliptic case (p. 39-40). P and Q the first and second collision points (Fig. 3), 2 tau the angle between them, collisions at f_p = tau and f_p = -tau
(f_p the true anomaly of the primaries' orbit). Case (i): P = Q, 2 tau = 2 pi p; the ellipse is described q times in 2 tau, period 2 pi p/q,
semi-major axis (p/q)^(2/3); "an elliptic orbit which collides with P2 at some point and with a rational mean motion, will be a solution with
consecutive collisions." Case (ii): P /= Q, P2 at t = 0 at a point R on the bisectrix of P P1 Q (pericentre or apocentre of P2's orbit); eight
initial configurations.

Coordinates (eqs. 11, 12; READ): P2: X = eps_p r_p cos f_p, Y = eps_p r_p sin f_p, t = E_p - eps_p e_p sin E_p. P3: X = eps a (eps'' cos E - e),
Y = eps eps' a sqrt(1 - e^2) eps'' sin E, t = a^(3/2) (E - eps'' e sin E). Signs: eps = +1/-1 if the pericentre of P3 is positive/negative; eps' = +1/-1
direct/retrograde; eps'' = +1/-1 if P3 begins at pericentre/apocentre; eps_p = +1/-1 if P2 begins at pericentre/apocentre ("defined as in Henon's
paper" for the first three; eps_p is new). "If eps'' = -1 or eps_p = -1, then the eccentric and true anomalies are measured from the apocenter."

### 2.2 The equations (eqs. 13-17, READ)

With K_p = ((1 + e_p)/(1 - e_p))^(1/2), tan(f/2) = sqrt((1+e)/(1-e)) tan(E/2), and the collision at Q at f_p = tau, E = eta:

(13a) eps_p r_p cos tau = eps a (eps'' cos eta - e)
(13b) eps_p r_p sin tau = eps eps' a sqrt(1 - e^2) eps'' sin eta
(13c) 2 arctan(K_p^(-eps_p) tan(tau/2)) - eps_p e_p sin[2 arctan(K_p^(-eps_p) tan(tau/2))] = a^(3/2) (eta - eps'' sin eta).

[(13c) as printed has "eps'' sin eta" without the factor e on the right; from (12) it must read a^(3/2) (eta - eps'' e sin eta); the printed (17) below
is nevertheless the correct elimination. INFERRED misprint.] Three equations, four unknowns (tau, eta, a, e): an infinity of solutions.
(14) r_p = a (1 - eps'' e cos eta). (15) a = r_p (1 - eps_p eps eps'' cos tau cos eta)/sin^2 eta, e = (eps'' cos eta - eps_p eps cos tau)/(1 - eps_p eps eps''
cos eta cos tau). (16) r_p = (1 - e_p^2)/(1 + eps_p e_p cos tau). Substituting (15), (16) into (13c) gives the single implicit equation (17):

( 2 arctan(k_p^(-eps_p) tan(tau/2)) - eps_p e_p sin[ 2 arctan(k_p^(-eps_p) tan(tau/2)) ] ) |sin eta|^3
 = ( (1 - e_p^2)/(1 + eps_p e_p cos tau) )^(3/2) (1 - eps_p eps eps'' cos eta cos tau)^(1/2) x
   [ eta (1 - eps_p eps eps'' cos eta cos tau) - sin eta (cos eta - eps_p eps eps'' cos tau) ],   k_p = (1-e_p)/(1+e_p))^(1/2) [printed k_p in (17), K_p in (13)].

"an implicit equation that relates the angles tau and eta only." The elements follow from (15) once (tau, eta, signs) are fixed. (The 1991 Part II
eq. (1) is this equation written with the same symbols; the e_p = 0 form is Henon's eq. 30.)

Parabolic and hyperbolic cases (pp. 41-42, eqs. 18-20, READ). Hyperbolic: X = eps a (e - cosh F), Y = eps eps' a sqrt(e^2 - 1) sinh F, t = a^(3/2)
(e sinh F - F) (18). Parabolic: X = eps (p/2)(1 - s^2), Y = eps eps' p s, t = p^(3/2) (s/2 + s^2/6) [sic; the correct cubic is s^3/6: COMPUTED,
dt = r^2 d theta/sqrt(p) with r = p(1+s^2)/2, tan(theta/2) = s gives t = p^(3/2)(s/2 + s^3/6); the printed s^2 is a misprint]. Timing equations:

(19, parabolic) 2 arctan(K_p^(-eps_p) tan(tau/2)) - eps_p e_p sin[...] = (1/3) ((1 - e_p^2)/(1 + eps_p e_p cos tau))^(3/2) sqrt(1 - eps_p eps cos tau) (2 + eps_p eps cos tau)

(20, hyperbolic) ( 2 arctan(K_p^(-eps_p) tan(tau/2)) - eps_p e_p sin[...] ) |sinh eta|^3 = ( (1 - e_p^2)/(1 + eps_p e_p cos tau) )^(1/2) (1 - eps_p eps cos tau cosh eta)^(1/2) x
 | sinh eta (cosh eta - eps_p eps cos tau) - eta (1 - eps_p eps cos tau cosh eta) |.

[(20) as printed has the prefactor power 1/2 where (17) has 3/2 for the same factor; the 3/2 should carry over (the right side is r_p^(3/2)-scaled). I did
not test (20) numerically: unresolved, possible misprint.] "For the particular case e_p = 0 and eps_p = +1 (for e_p = 0 and eps_p = -1, we have the same
situation if we turn the whole figure by an angle pi), Equations (17), (19) and (20) become those ones obtained by Henon [6], in the circular case."

### 2.3 Solutions found (pp. 42-49, READ)

Method (p. 42): "solved numerically (by a method of continuation to calculate roots of a system of equations) Equations (17), (19) and (20) for different
values of e_p, between 0 and 1." No step sizes or tolerances printed.

Parabolic case (p. 43): for fixed e_p equation (19) has a unique solution: for eps_p = +1 and eps = eps' = -1, "tau/pi = 0.2318 if e_p = 0.5" (Fig. 5 caption
shows tau/pi = 0.2317); for eps = +1 tau = 0 (P and Q coincide when the motion begins). For eps_p = -1 one solution for eps = +1, eps' = -1, "tau/pi = 0.1055 if
e_p = 0.5", and tau = 0 for eps = -1. [COMPUTED: eq. (19) at e_p = 0.5 gives 0.231149 for eps_p = +1 and 0.105578 for eps_p = -1; at e_p = 0 it gives 0.163926, the
parabolic value of Henon's Table 2. See section 6: the printed 0.2318 and 0.2317 are not reproduced.]

Hyperbolic case: (20) gives a simple family for eps = eps' = -1 if eps_p = +1 and for eps = +1, eps' = -1 if eps_p = -1; Fig. 5 (e_p = 0.5): the (tau/pi, eta/pi)
diagram with orbit shapes at tau/pi = 0.03, 0.1, 0.15, 0.2317 for eps_p = +1; "For eta = 0 and the corresponding value of tau, we obtain exactly the
unique solution of the parabolic case."

Elliptic case: "an infinity of families of solutions", named A_i, B_i, C_ij following Henon; "The description of these families is similar to the circular
case." "Due to the evolution of these curves, which is not continuous when e_p varies from 0 to 1, we show the different diagrams of solutions for e_p = 0.5
and e_p = 0.98" (Figs. 6-9: Fig. 6 eps_p = +1, e_p = 0.5, tau < 5.5 pi, eta < 7 pi; Fig. 7 eps_p = -1, e_p = 0.5; Fig. 8 eps_p = +1, e_p = 0.98;
Fig. 9 eps_p = -1, e_p = 0.98; sign triplets (eps, eps', eps'') on the curves). Orbit shapes: A0 (Fig. 10, eps_p = +1, half orbit f_p from -tau to 0, panels
tau/pi = 0.231, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4); A1 (Fig. 11, panels tau/pi = 4, 2, 1.5, 2.5, 3.5, 1.75, 1.75, 3,
3, 1.5, 2, 3.5, 2.5, 1.2102, 2.25, 4); A2 (Fig. 12, panels 4, 2, 1, 2.5, 3.5, 1.5, 1.25, 3, 3, 1, 1.5, 3.5, 2.5, 0.967, 2, 4); C34 (Fig. 13, panels 3, 3, 2.95,
3.02, 3.03, 2.97, 3, 3.1, 3, 2.98, 3). [The panel values are the plotted tau/pi labels; the correspondence of each panel to its (a, e) is not printed.]

Qualitative statement (p. 44, READ): "A qualitative difference in these diagrams, when varying e_p, is the appearance and disappearance of families C_ij."
Every C_ij contains a double point (i pi, j pi) (i, j in Z) "it corresponds to an ellipse, described by P3, which is tangent to that one described by P2
(the points P and Q coincide). The parameters a, e of this ellipse verify:

a = (i/j)^(2/3),   e = | 1 - (j/i)^(2/3) (1 + (-1)^(i+1) eps_p e_p) |."

Existence of C_ij (eq. 21, READ): for j > i, once i is fixed,   j <= ( 2/(1 + (-1)^(i+1) eps_p e_p) )^(3/2) i.   "For the particular case e_p = 0, (21) becomes the
condition given by Henon." At e_p = 0.5 (eps_p = +1): (2/(1+e_p))^(3/2) = 1.539 and (2/(1-e_p))^(3/2) = 8; so "if i/pi = 1 then j/pi <= 1.539; i/pi = 2: j/pi <= 16;
i/pi = 3: <= 4.617; i/pi = 4: <= 32; i/pi = 5: <= 7.695, etc." (the notation i/pi, j/pi is as printed; i, j are integers). "Therefore, comparing with the diagram for e_p = 0
we can conclude that, when e_p = 0.5 and eps_p = +1, families C12, C35, C36, ... disappear (if i/pi is odd) and new families C26, C27, ..., C2,16, etc appear (if i/pi
is even) (see Figures 6 and 4 of Henon's paper)." At e_p = 0.98: (2/(1+e_p))^(3/2) = 1.015 and (2/(1-e_p))^(3/2) = 1.000 [printed; the true value is 1000, a lost
thousands separator: COMPUTED], "so the first family C_ij, with i/pi odd, will be C67,68 (if eps_p = +1) and C68,69, with i/pi even, if eps_p = -1."

## 3. Double-collision orbits (pp. 49-55)

Setting (p. 49, READ): P3 collides with both primaries; the number of such orbits is characterised "according to the value of the sidereal energy h = -1/(2a)".
P3 on a rectilinear sidereal orbit (e = 1) at angle alpha to the X axis, beginning at pericentre or apocentre (eps_p = 1 only for these equations in the text;
four configurations); P3: X = a(1 - eps'' cos E) cos alpha, Y = a(1 - eps'' cos E) sin alpha, t = a^(3/2)(E - eps'' sin E); at the collision with P2, f_p = tau,
E_p = sigma, E = eta, hence alpha = tau and

(22a) eps_p r_p cos tau = a (1 - eps'' cos eta) cos tau,  (22b) eps_p r_p sin tau = a (1 - eps'' cos eta) sin tau,  (22c) sigma - eps_p e_p sin sigma = a^(3/2)(eta - eps'' sin eta).
(23) r_p = a (1 - eps'' cos eta). (24) ( 2 arctan(K_p^(-eps_p) tan(tau/2)) - eps_p e_p sin[...] ) (1 - eps'' cos eta)^(3/2) = (1 - e_p^2)(eta - eps'' sin eta)/(1 + eps_p e_p cos tau)^(3/2).
For e_p = 0 (25): tau (1 - eps'' cos eta)^(3/2) = eta - eps'' sin eta. The diagram of solutions (Fig. 15, e_p = 0.5, eps_p = +1 and -1) restricts eta < 2 pi if eps'' = +1 and
eta < 3 pi if eps'' = -1 ("we require both collisions before P3 collides with P1 again").

Circular-case theorem (p. 51, READ, quoted): "THEOREM. If h < -1, there exist no double collision orbits. If h = -1, there exist two double collision orbits. If -1 < h < 0
there exist two or four double collision orbits (for more details see Pinol [18])." "In general, there will be an even number, because, if tau, eta is a solution, so will
-tau, -eta (the symmetrical solution respect the synodical x-axis)."

Elliptic case (pp. 51-55). Fix a (equivalently h); the system (26): (a) sigma - eps_p e_p sin sigma = a^(3/2)(eta - eps'' sin eta); (b) a(1 - eps'' cos eta) = 1 - e_p cos sigma.
[(b) as printed omits eps_p; COMPUTED: with the factor, a(1 - eps'' cos eta) = 1 - eps_p e_p cos sigma, the printed tangency values of Figs. 18-19 are reproduced, and the
alternative form (27) a(1 + eps_p e_p cos tau)(1 - eps'' cos eta) = 1 - e_p^2 carries eps_p. See section 6.] The equivalent form (27) pairs (13c-type time equation) with that relation.
"The curve (26a) eta = eta(sigma) passes by (0,0) and is monotone increasing, while (26b) is 2 pi periodic in eta and sigma. Therefore, for eps'' = +1, we shall have 0, 1 or 2 intersection
points ... If eps'' = -1, then eta varies from 0 to 3 pi, and there are 0, 1, 2 or 3 possible intersection points." Figs. 16-19: solutions of (26) at a = 0.3, 0.345, 0.75 (eps_p = eps'' = +1);
0.65, 0.71, 0.75 (eps_p = +1, eps'' = -1); 0.4, 0.585, 0.75 (eps_p = -1, eps'' = +1); 0.4, 0.461, 0.55, 0.75 (eps_p = -1, eps'' = -1). [The e_p of Figs. 16-19 is not stated; Fig. 20 is a function of e_p.]

PROPOSITION (p. 54, READ, the table transcribed): the number of intersection points of (26) as a function of the semi-major axis a.

| eps_p | semi-major axis a | eps'' = +1 | eps'' = -1 |
|---|---|---|---|
| +1 | 0 < a < (1 - e_p)/2 | 0 | 0 |
| +1 | a = (1 - e_p)/2 | 0 | 1 |
| +1 | (1 - e_p)/2 < a < a^{++} | 0 | 1 |
| +1 | a = a^{++} | 1 | 1 |
| +1 | a^{++} < a < a^{+-} | 2 | 1 |
| +1 | a = a^{+-} | 2 | 2 |
| +1 | a^{+-} < a | 2 | 3 |
| -1 | 0 < a < a^{--} | 0 | 0 |
| -1 | a = a^{--} | 0 | 1 |
| -1 | a^{--} < a < a^{-+} | 0 | 2 |
| -1 | a = a^{-+} | 1 | 2 |
| -1 | a^{-+} < a < (1 + e_p)/2 | 2 | 2 |
| -1 | (1 + e_p)/2 <= a | 2 | 3 |

"where a^{++}, a^{+-}, a^{-+}, a^{--} are those values of a for which there is a tangential intersection. Figure 20 shows that a^{++} < a^{+-}, for eps_p = +1, and a^{--} < a^{-+} for eps_p = -1; and the
two signs are those of eps_p and eps''." "For each intersection point, we obtain two double collision orbits, plus the two symmetrical ones. There is an exceptional case: eps'' = -1, eta < pi,
a^(3/2) + (eps_p e_p sin sigma - sigma)/pi not in N [and] a^(3/2) (3 pi/2 - eta) not in N [the second line is clipped across the page break and read with low confidence]; for that intersection point,
we obtain one double collision orbit (plus its symmetrical one)." The statement is "numerical evidence" ("we have numerical evidence of the following results", p. 52): no proof.

Fig. 20 (p. 55; READ approximate): the four tangency curves a(e_p) all start at a = 0.5 at e_p = 0; "+-" rises to about 0.77 near e_p = 0.9; "-+" rises to about 0.64; "--" stays near 0.48 to 0.5;
"++" falls to about 0.02 at e_p near 0.98. [Digitisation from the figure, not exact.]

## 4. Reduction to Henon 1968 (e_p = 0)

Printed (p. 42, p. 44, p. 50): at e_p = 0 and eps_p = +1, eqs. (17), (19), (20) are Henon's. At e_p = 0, eps_p = -1 the situation is the same turned by pi. Eq. (21) at e_p = 0 gives Henon's existence
condition for C_ij. Eq. (25) is the e_p = 0 double-collision equation. Henon's notation (digest sec. 1.1: eps, eps', eps'' = sigma0, sigma1, sigma2; equations numbered 3, 4, 25, 30) maps one-to-one onto this
paper's (eps, eps', eps''); the new signs are eps_p, e_p. COMPUTED checks at e_p = 0:
- Eq. (19), eps_p = +1, eps = -1: tau/pi = 0.163926; Henon Table 2 first row (parabolic orbit) prints 0.16393. Eq. (19) with eps_p = -1, eps = +1 gives the same value, as the "turn by pi" statement requires.
- Eq. (15) with (tau, eta) = (0.17 pi, 0.16734 pi), eps = -1, eps'' = +1: a = 6.9272, e = 0.98922; Henon Table 2 row prints 6.92689 and 0.98922 (a is sensitive to the five-digit eta: relative difference 5e-5).
- Eq. (17) residual at that row: left 0.0674949, right 0.0674948.

## 5. What is not in this paper (READ and INFERRED)

No mu > 0 dynamics, no second-species periodic orbit, no corrector, no stability or Floquet result, no table of (tau, eta, a, e, x, C), no numerical value of any element at any e_p except the few quoted in
the text, no proof that the numerically found families persist for e_p near 1 ("not continuous" is stated), no discussion of hyperbolic or parabolic families beyond the equation and Fig. 5. Hitzl's stability work
on the circular mu = 0 problem is cited (p. 33) as a reference only. The 1991 papers (digest) give the O(mu) matching and the mu = 1e-6 families; the e_p dependence of the OCC enumeration is only here.

## 6. Verification I performed (COMPUTED, 2026-10-04)

| item | result |
|---|---|
| Eq. (19) at e_p = 0, eps_p = +1, eps = -1 | tau/pi = 0.163926 (Henon Table 2: 0.16393) |
| Independent solution of the parabolic collision problem from first principles (three equations X, Y, time with P2's Kepler equation, P3's parabola, unknowns tau, p, s) | e_p = 0: 0.163926; e_p = 0.5, eps_p = +1: 0.231149 (same as eq. 19) |
| Eq. (19), e_p = 0.5, eps_p = +1, eps = -1 | 0.231149; PRINTED 0.2318 (text), 0.2317 (Fig. 5 caption): difference 6e-4 to 7e-4, not reproduced |
| Eq. (19), e_p = 0.5, eps_p = -1, eps = +1 | 0.105578; printed 0.1055 (agrees to the printed digits) |
| (2/(1+e_p))^(3/2), (2/(1-e_p))^(3/2) at e_p = 0.5 | 1.5396 and 8 (printed 1.539 and 8) |
| Eq. (21) bounds at e_p = 0.5: i = 1, 2, 3, 4, 5 | 1.539, 16, 4.617, 32, 7.695 (printed, all reproduced) |
| e_p = 0.98 | (2/1.98)^(3/2) = 1.0153 (printed 1.015); (2/0.02)^(3/2) = 1000 (printed "1.000": lost separator) |
| first C_ij, e_p = 0.98: i odd, eps_p = +1: smallest i with a j > i satisfying j <= 1.0153 i | i = 67, j = 68 (C67,68); eps_p = -1, i even: i = 68, j = 69 (C68,69): both as printed |
| Fig. 7 (eps_p = -1, e_p = 0.5), eq. (21): i = 1 gives j <= 8 (C12 to C17 drawn), i = 2 gives j <= 3.08 (only C23 drawn) | consistent with the figure |
| Tangency of system (26), e_p = 0.5, solving G1 = G2 = det = 0 for (a, sigma, eta), with (26b) written a(1 - eps'' cos eta) = 1 - eps_p e_p cos sigma | a^{++} = 0.344799 (Fig. 16 plots 0.345); a^{+-} = 0.710487 (Fig. 17: 0.71); a^{-+} = 0.584733 (Fig. 18: 0.585); a^{--} = 0.463637 (Fig. 19: 0.461, the plotted value) ; ordering a^{++} < a^{+-}, a^{--} < a^{-+} as printed; (1 - e_p)/2 = 0.25 and (1 + e_p)/2 = 0.75 appear as the degenerate solutions |
| With (26b) as printed (no eps_p) | a^{-+} and a^{--} do not come out near 0.585 and 0.461 (0.259/0.823 and 0.312/0.548/0.644), so the print omits eps_p (INFERRED misprint) |
| 1991 starting orbits (A0 tau = pi, A1 tau = 3 pi, B1 tau = 2 pi, e_p = 0, eps_p = +1, eps eps'' = -1, -1, +1) at eta/pi = 0.53286, 0.68000, 0.63192 (derived from a0 = |x|/2 and cos eta = 1/a0 - 1, so not independent of the t1 = k pi check) | eq. (17) residual left - right = -1.7e-5, -2.0e-5, -3.8e-5 against sides of size 3.1, 5.7, 4.8 |
| code identity: `core.er3bp.er3bp_eom` at mu = 0 | equals eq. (1) (read from source; no run) |

Not tested: eq. (20) (hyperbolic) and (24); the full set of family diagrams; the double-collision counts of the Proposition table (my own brute-force counting was unreliable and is not reported).

## 7. Printed defects and ambiguities (kept as printed in the transcription above)

1. p. 43 and Fig. 5 caption: parabolic tau/pi at e_p = 0.5, eps_p = +1: printed 0.2318 (text) and 0.2317 (caption); the equation and an independent solution give 0.23115.
2. p. 42: parabolic time law printed with s^2/6; correct s^3/6.
3. p. 41 eq. (13c): right side misses the factor e in the eps'' e sin eta term (eq. 12 and eq. 17 are consistent with it).
4. p. 42 eq. (20): prefactor power 1/2 where the analogous (17) has 3/2 (not tested).
5. p. 51 eq. (26b): missing eps_p (shown above); p. 44: "1.000" for 1000.
6. p. 36 Proposition 3 proof: case labels (ii) and (iii) both read as q, p parity combinations that duplicate (i).
7. The notation "i/pi", "j/pi" in the C_ij condition (p. 44) means the integers i, j of the double point (i pi, j pi).
8. e_p of Figs. 16-19 is not stated; the agreement of the four plotted a values with the e_p = 0.5 computation above suggests e_p = 0.5 (INFERRED).

## 8. Test-ready numbers (sourced to this paper's print; each test must label the quantity as printed or computed)

| id | quantity | value | source | status |
|---|---|---|---|---|
| T1 | eq. (21) bounds, e_p = 0.5, eps_p = +1, i = 1..5 | 1.539, 16, 4.617, 32, 7.695 | p. 44, printed | reproduced (COMPUTED) |
| T2 | (2/(1+e_p))^(3/2) and (2/(1-e_p))^(3/2) at e_p = 0.5 | 1.539, 8 | p. 44 | reproduced |
| T3 | same at e_p = 0.98 | 1.015, 1000 (printed 1.000) | p. 44 | 1.015 reproduced; second printed digits are a typesetting loss |
| T4 | first C_ij at e_p = 0.98 | C67,68 (eps_p = +1), C68,69 (eps_p = -1) | p. 49 | reproduced |
| T5 | families gained and lost at e_p = 0.5, eps_p = +1 | C12, C35, C36 lost; C26 to C2,16 gained | p. 44 | consistent with eq. (21) |
| T6 | double-point ellipse of C_ij | a = (i/j)^(2/3); e = abs(1 - (j/i)^(2/3)(1 + (-1)^(i+1) eps_p e_p)) | p. 44 | derived from the tangency condition; not otherwise tested |
| T7 | parabolic tau/pi, e_p = 0.5 | 0.2318 (eps_p = +1), 0.1055 (eps_p = -1) | p. 43 | NOT reproduced for +1 (0.23115); reproduced for -1 (0.10558). Use as a strict expected failure or record the discrepancy; do not copy |
| T8 | e_p = 0 limit of eqs. (17), (19) | Henon's Tables 1-9 (digest of Henon 1968) | pp. 42, 44 | the project's Henon tests apply with e_p = 0 |
| T9 | tangency a values (plotted) | 0.345, 0.71, 0.585, 0.461 | Figs. 16-19 | 0.3448, 0.7105, 0.5847, 0.4636 computed (Fig. 19 differs by 0.003 from the plotted 0.461) |
| T10 | double-collision count (Proposition) | table in sec. 3 | p. 54 | partly tested; boundaries (1 +- e_p)/2 appear in the tangency solve |
| T11 | circular theorem | h < -1: none; h = -1: two; -1 < h < 0: two or four | p. 51 | cited from Pinol, not tested |
| T12 | eq. (1) at mu = 0 | equals `er3bp_eom` at mu = 0 | p. 34 | READ vs source |

## 9. Techniques applicable to the project's problems

### 9a. #899, an elliptic-problem enumerator beside the Henon one

- The enumeration of generating arcs at e_p > 0 is eq. (17) with the signs (eps, eps', eps'', eps_p), eq. (15) for (a, e), eq. (16) for the primaries' distance at the collisions, and the parabolic and hyperbolic
  variants (19), (20). The Henon module (task step 1 of #899: `second_species_arcs`, timing equation with the e = 1 factorisation, controls from Henon Tables 1-9 and the parabolic value 0.163926) is the e_p = 0 slice: make
  e_p and eps_p parameters (default 0 and +1) and the Henon tables become the e_p = 0 regression tests. T7 and T1-T5 are the new tests; note that (17) has the same structure for any e_p, so the
  e = 1 factorisation of Henon's equation should be re-derived for e_p /= 0 (not done here; at e = 1, (15) gives e = 1 iff cos tau = -eps_p eps for eps'' = -1, i.e. tau = k pi: INFERRED from (15), consistent with 1991's
  t1 = k pi).
- Solve it as the paper does: a continuation in e_p from 0 (roots of a system by continuation), because the family structure is "not continuous when e_p varies from 0 to 1" (p. 44): families C_ij appear and
  disappear; do not assume a root follows the same label. Use (21) as the counting oracle: for given (i, eps_p, e_p) it predicts exactly which C_ij exist; a missing or extra family is a bug.
- The number of generating arcs at a given energy: the double-collision Proposition (p. 54) is the only counting result, and only for double collisions. For single-collision OCC there is no count; the enumerator
  must produce one and the C_ij oracle checks its C_ij part.
- Period class: an OCC chain in the elliptic problem is periodic only when the collisions occur at f_p = tau = k pi (pericentre or apocentre of P2) and the sidereal period commensurable (Proposition 2: n = q/p,
  T_syn = 2 pi p). The enumerator should flag, for each (tau, eta) root, whether tau is a multiple of pi (symmetric periodic candidate) or not (only an arc).
- Note the sign eps_p: apocentre starts of P2 are a separate set of eight configurations; the e_p = 0, eps_p = -1 set is the e_p = 0, eps_p = +1 set turned by pi (p. 42), so at e_p = 0 do not double count.

### 9b. #912, the elliptic periodic-orbit code lacks the published method

- Proposition 2 is the published statement of WHICH initial epochs make a symmetric synodic orbit: pericentre (case ii), apocentre (case iii) or the j pi epochs (case i). This is exactly the "starting true anomaly
  argument (pericentre or apocentre group)" that #912 says the project lacks; use it as the specification for the start-condition argument and for the half-period (T_syn/2 = pi p) crossing. The p of
  Proposition 2 is the denominator of n = q/p, so periods are 2 pi p, matching the 1991 result (period 2 k pi) and the project's catalogue rows with periods that are multiples of 2 pi.
- Eq. (1) at mu = 0 identical to `er3bp_eom`: the base two-body (mu = 0) limit is a free positive control for any new elliptic corrector or segment product: with mu = 0 every orbit of the rotating pulsating
  problem is a Kepler orbit in disguise, so the monodromy over p revolutions is known in closed form from the two-body STM (rotation/pulsation composed with the Kepler flow). An elliptic multi-segment corrector
  that cannot return the Kepler orbit at mu = 0 has an error; this is a model-independent test of the product of segment matrices that #912 asks for.
- Characteristic curves (eqs. 7, 10) give the (I0, x) location of every symmetric periodic mu = 0 orbit: use them to seed or to test the choice of x0 for a pericentre-start orbit at mu = 0 (x0 = a(1 -+ e)/r_p).

### 9c. #925, the Mako & Salamon elliptic stability bands (still unreproduced)

This paper has no stability content, so it cannot settle the bands. What it offers: (i) the mu = 0 limit of the elliptic stability problem is two-body and exactly solvable (as above): compute the Mako-Salamon
case's monodromy at mu = 0 analytically and as a limit of `core.er3bp` to separate a model defect from a paper slip (the same limit procedure #925's "independent formulation" step asks for, but with a closed-form
answer); (ii) eq. (1) agreeing with `er3bp_eom` removes one place where a defect could hide at mu = 0. It does not help with the mu > 0 part. Do not infer anything about the Mako-Salamon bands from this paper.

### 9d. #930, the elliptic corrector's silent acceptance

- The mu = 0 family gives a SOURCED closed-form target for the corrector's closure test: at mu = 0 a symmetric orbit of Proposition 2 must close exactly after T_syn = 2 pi p, so a corrector (and the
  `genome.er3bp_continuation` caller) can be unit-tested on the mu = 0 case, where a failing full-period closure must raise, with an exact expected result and no integrator-resolution excuse (no 1e-7 periapsis at
  mu = 0). A defect of the kind #930 describes (a full-period closure only logged as a warning) would be exposed by a deliberately perturbed mu = 0 start whose Kepler orbit has the wrong period.
- The paper's own symmetry (eq. above) gives the right residual: perpendicular axis crossings at f = 0 and f = pi p are event-located conditions, not fixed-time residuals (what #930 recommends); Proposition 1
  says the crossings ARE the pericentre/apocentre passages, so locating them by event (dr/df = 0) is the paper's definition.

### 9e. #906, hardening the demanded-turn gate

- The elliptic generalisation of the equal-speed node condition: at a collision the incoming and outgoing relative speeds are equal for e_p = 0; for e_p > 0 the paper's (13)-(17) place the collision at f_p = +-tau with
  P2's velocity different at the two collisions (different r_p), so the demanded-turn gate, if ever applied to an elliptic-problem chain, must evaluate relative velocities at the two different phases
  (the primaries' velocity and distance differ at f_p = +tau and -tau only through cos tau, i.e. the equal-distance symmetric case still gives equal |relative speed| only when the two collision points are mirror images:
  they are, by the symmetry in section 1: INFERRED). Treat this as a flag for #906(a) (inputs that are not body-relative velocities) in an elliptic setting; no test is available from this paper.
- (b) Tangency and resonance cases: the double point (i pi, j pi) of C_ij is an ellipse tangent to P2's orbit (P = Q coincide): the demanded turn is zero there. This is the elliptic version of the Hitzl-Henon
  near-resonant seeds that #906's amendment says must be returned as "indeterminate", not rejected. The explicit location a = (i/j)^(2/3), e = |1 - (j/i)^(2/3)(1 + (-1)^(i+1) eps_p e_p)| lets the gate recognise a
  near-tangent seed at any e_p (the nearest low-order p/q in the amendment becomes the nearest (i, j) of the C_ij double points at the current e_p).

## 10. The #896 item (j) finding (the p163 starting orbits of the 1991 digest) and whether this paper bears on it

The #896 note (`docs/notes/2026-10-04-896-published-checks-added.md`, item 8, "Gomez & Olle 1991 II, p163 starting orbits") found that the printed (x, ydot) of the three starting orbits does not give the
axis crossing at k pi (the crossing comes 1.35e-6 to 1.46e-6 late with xdot -1.1 to -1.26) and the passage is at 1e-5 to 3.5e-5 instead of about 1e-7; that the computed members whose periapsis lies on the axis at
t = k pi differ by (dx, dydot) = (+5.89e-6, -5.89e-6), (+9.08e-6, -9.08e-6), (-1.537e-5, +1.538e-5) with x + ydot equal to the printed sum; and (INFERRED there) that the printed parameter is right and the printed x is off.

Does this paper bear on it? Only weakly, and in one direction. (a) It is a mu = 0 paper: it contains no O(mu) initial-condition offsets, no K0, K1, a_pq, and no mu > 0 orbit, so it cannot confirm or explain the
size or sign of (dx, dydot). (b) It does confirm the mu = 0 generating arcs the starting orbits come from: the rectilinear ellipse with collision at tau = k pi is the e = 1 end of the A0, A1, B1 families
(from eq. (15), e = 1 iff cos tau = -eps_p eps at eps'' = -1) and eq. (17) is satisfied by the generating (tau, eta) of the three printed orbits to relative 4e-6 to 8e-6 (section 6; this check uses a0 = |x|/2 from the printed x,
so it is not independent of the printed x to that order). A discrepancy of dx of order 1e-5 in a quantity of size 2 to 4 is 2.5e-6 to 7e-6 relative, below what eq. (17) as a mu = 0 test can resolve. So the
1986 paper neither supports nor contradicts the #896 inference that the printed x is off by O(1e-5); that inference rests on the project's own integration of the mu = 1e-6 model. If anything the O(mu) shift is the
whole story: x, ydot printed to 11 digits are generating-orbit values plus an offset the paper does not state how it computed (Theorem 10's K0, K1 route needs a21, a24), and an error of order mu in a printed x
would be a first-order (O(mu)) correction omitted or mis-signed, not a transcription issue (INFERRED; untested).

## 11. Follow-ups (not registered as tasks)

1. Implement the e_p-parametrised enumerator (9a) and make eqs. (15), (17), (19), (21) tests: T1-T5 as sourced tests, T7 as a recorded non-reproduction.
2. A mu = 0 elliptic monodromy test (9b, 9c, 9d): closed-form Kepler flow in the rotating pulsating frame; check `er3bp_monodromy` and the corrector at mu = 0.
3. Decide eq. (20) (hyperbolic) by a first-principles solution like the one in section 6 and record the result in the Henon/Gomez-Olle test notes.
4. Check the e = 1 factorisation of eq. (17) for e_p /= 0 (needed for the rectilinear seeds that the 1991 theory uses).
5. Re-verify the tangency values of Figs. 16-19 and the counts of the Proposition table with a proper intersection counter (my brute-force count was unreliable).

## 12. References (as printed, p. 55) and corpus status

Held = a file is in the corpus directory (checked by listing); "digested" = a digest exists in `docs/notes/`.
- [1] A. D. Brjuno, 1978, Celest. Mech. 18, 9. Held and digested (`docs/notes/2026-10-04-digest-brjuno-1978-periodic-solutions-arcs-mu-0.md`); the paper's own reference for characteristic curves of mu = 0 periodic orbits of first and second kind.
- [2] G. Gomez and J. Llibre, 1981, Celest. Mech. 24, 335. Not held.
- [3] P. Guillaume, 1969, Astron. Astrophys. 3, 57. Not held (the held Guillaume 1973 paper, Celest. Mech. 8, 199, is a different work).
- [4] P. Guillaume, 1975, Celest. Mech. 11, 213 (Linear analysis of one type of second-species solutions). Held and digested (`docs/notes/2026-10-04-digest-guillaume-1975a-linear-analysis-second-species.md`).
- [5] P. Guillaume, 1975, Celest. Mech. 11, 449 (Extension of Breakwell-Perko's matching theory). Held and digested (`docs/notes/2026-10-04-digest-guillaume-1975-extension-breakwell-perko-matching.md`).
- [6] M. Henon, 1966 [as printed; the digest records 1968], Bull. Astron. Paris 1, 377. Held and digested (Henon 1968).
- [7] J. Henrard, 1980, Celest. Mech. 21, 83. Not held.
- [8] D. L. Hitzl, 1977, AIAA J. 15, 1410. Not held.
- [9] D. L. Hitzl and M. Henon, 1977, Celest. Mech. 15, 421. Held and digested (`docs/notes/2026-10-04-digest-hitzl-henon-1977-critical-generating-orbits.md`).
- [10] D. L. Hitzl and M. Henon, 1977, Acta Astronautica 4, 1019. Not held.
- [11] J. Llibre, 1982, Celest. Mech. 28, 83. Not held.
- [12] L. M. Perko, 1964, AIAA J. 2, 2187. Not held.
- [13] L. M. Perko, Ph. D. Dissertation, Stanford University, 1965. Not held.
- [14] L. M. Perko, 1974, SIAM J. Appl. Math. 27, 200. Not held (the 1991 Part I cites the same page range as volume 41, 1981: one of the two is a slip).
- [15] L. M. Perko, 1981, SIAM J. Appl. Math. 41, 181. Held and digested (Perko 1981b, first-second species bifurcation).
- [16] L. M. Perko, 1981, Celest. Mech. 34, 155 [as printed; the held paper is Celest. Mech. 24, 155 (O(mu^nu), 0 < nu < 1 near-moon passage), so the volume 34 is a slip, INFERRED from the page number and the 1991 lists]. Held (file perko-1981-...-celest-mech-24-155); digest not checked.
- [17] L. M. Perko, 1983, Celest. Mech. 30, 115. Not held.
- [18] C. Pinol, Master Thesis, Universitat Autonoma de Barcelona, 1983. Not held.
- [19] V. Szebehely, 1967, Theory of Orbits, Academic Press. Held.
- [20] A. Wintner, 1947, The Analytical Foundations of Celestial Mechanics, Princeton University Press. Not held.
Reference status was assembled from a directory listing and the digests named (re-checked after the corpus grew today); I did not open the held files here. Also held and relevant but not cited by this paper: Perko 1976 (CM 14, 395), Perko 1977 (CM 16, 275), Guillaume 1973 (CM 8, 199), Bruno 1981 (CM 24, 255).
