# Digest: Font, Nunes and Simo 2009, a numerical study of the orbits of second species of the planar circular RTBP

Date 2026-10-04. Purpose: establish what is known, numerically and analytically, for ONE small secondary about
orbits with consecutive close encounters, so that the #890 result (a periodic orbit of the planar concentric
circular restricted four-body problem, Uranus plus Titania plus Oberon, one flyby of each moon per cycle; see
`docs/notes/2026-10-04-890-titania-oberon-candidate.md` and `docs/notes/2026-10-04-890-literature-check.md`)
is worded correctly and the one-moon results can serve as checks. Companion digests:
`docs/notes/2026-10-04-digest-font-nunes-simo-2002-consecutive-quasi-collisions.md` (the analytic base paper, to be
read with this one) and `docs/notes/2026-10-04-digest-bolotin-mackay-2006-nonplanar-second-species.md`.

Citation: J. Font, A. Nunes and C. Simo, "A numerical study of the orbits of second species of the planar
circular RTBP", Celestial Mechanics and Dynamical Astronomy 103:143-162 (2009), DOI 10.1007/s10569-008-9176-z.
20 pages (journal pages 143 to 162; PDF page n is journal page n+142). Received 8 August 2008, accepted 17
November 2008, published online 25 December 2008.
Filed in the private paper corpus as
font-nunes-simo-2009-numerical-study-orbits-second-species-planar-circular-RTBP-cmda-103-143-doi-10.1007-s10569-008-9176-z.pdf.

Evidence labels: READ means read from the printed page in this session (all 20 pages as page images; the table
digits re-read at 400 dpi). INFERRED means my reading of what the print implies. Page numbers are journal pages.
Anything marked [unclear] could not be read with confidence.

## 0. Abstract (p. 143, READ, quoted in full)

"We present a numerical study of the set of orbits of the planar circular restricted three body problem which
undergo consecutive close encounters with the small primary, or orbits of second species. The value of the
Jacobi constant is fixed, and we restrict the study to consecutive close encounters which occur within a
maximal time interval. With these restrictions, the full set of orbits of second species is found numerically
from the intersections of the stable and unstable manifolds of the collision singularity on the surface of
section that corresponds to passage through the pericentre. A 'skeleton' of this set of curves can be computed
from the solutions of the two-body problem. The set of intersection points found in this limit corresponds to
the S-arcs and T-arcs of Henon's classification which verify the energy and time constraints, and can be used
to construct an alphabet to describe the orbits of second species. We give numerical evidence for the existence
of a shift on this alphabet that describes all the orbits with infinitely many close encounters with the small
primary, and sketch a proof of the symbolic dynamics. In particular, we find periodic orbits that combine
S-type and T-type quasi-homoclinic arcs."

## 1. Setting (pp. 144-146, READ)

Model (p. 145): the planar circular restricted three-body problem with "masses of the two primaries ... m_1 = 1 -
mu, m_2 = mu, mu in [0, 1/2], the angular velocity of their motion around the fixed centre of mass is 1 and
their positions in the plane are given in synodic coordinates (x, y) by (mu, 0) and (mu - 1, 0)", equations
(1.1), Jacobi function (1.2) identical to the 2002 paper: C_J = 2 Omega - (xdot^2 + ydot^2). "[They] are singular
at the points (x, y) = (mu, 0) and (x, y) = (mu - 1, 0), which correspond to collision of the massless body Z
with the big primary E and with the small primary M, respectively." The frame is the synodic frame; the
quantities are in units where the primaries' separation is 1 and their period is 2 pi.

Goal and restrictions (p. 145), quoted: "Our main goal is to study numerically the set of all the orbits that,
for a given value of the Jacobi constant C_J and of the mass parameter mu, undergo consecutive close encounters
with the small primary M. We restrict the study to orbits such that these consecutive encounters occur within
time intervals smaller than a multiple q of 2 pi, and such that the number of passages through the pericentre
between close encounters does not exceed p. We have considered C_J = 2.8 and q = p = 4 throughout, a choice for
which the alphabet of the symbolic dynamics is rich but still compatible with the concise presentation of a
complete numerical study. We have taken mu = 1e-4 in most of the examples, and the dependence on mu is
illustrated by comparison with mu = 1e-6."

Close encounter and in/out maps: as in the 2002 paper. "Let the set of initial conditions exiting the circle C of
radius mu^alpha around the small primary M be parameterized by an angle phi over C and an angle psi that
determines the direction of the velocity" (p. 145 to 146); "The value of alpha is set to 0.4 through the paper"
(p. 146). So the encounter distance is mu^(2/5): 0.0251 at mu = 1e-4, 0.0040 at mu = 1e-6 (INFERRED arithmetic,
in units of the primaries' separation). The OUT-map and IN-map are defined as in the 2002 paper (p. 156 to 157):
"The OUT-map sends a point on a p-q resonant strip to the position and velocity angles (phi, psi) of the orbit
with those initial conditions at the time of its first intersection with the circle C; the IN-map sends the
position and velocity angles (phi, psi) of an orbit entering the disk D around the small primary M to the
position and velocity angles of that orbit at the time that it crosses again the circle C leaving the disk D.
The return map on the strips is the composition of the IN- and the OUT-maps."

Regularisation: "numerically using the Levi-Civita transformation that regularizes Eq. (1.1) with respect to
collision with the small primary (Stiefel and Scheifele 1971)" (p. 150). The angle theta in the regularised
velocity plane relates to the synodic exit direction by psi = 2 theta (p. 151). Integration time cut-off 8 pi for
the manifold curves (p. 151).

Henon's classification (p. 144), quoted: "These limit orbits are formed by arcs, which are pieces of Keplerian
ellipses that begin and end in a collision, and an arc is called of type S (resp. T) if it begins and ends at
different points (resp. at the same point) on the ellipse. In the rotating frame of reference of the restricted
three body problem orbits formed by S arcs are symmetric with respect to the line joining the two primaries,
while orbits that contain T arcs are in general not symmetric." (The first clause is from my reading of the
page image at 100 dpi, "formed by arcs, which are pieces of"; verify before quoting the first clause.)

## 2. Historical and context statements (pp. 144-145, READ, quoted)

- "The results of Perko (1976), Guillaume (1975a,b), Bruno (1981) for the circular restricted three body problem,
  extended to the elliptic case in Gomez and Olle (1991), use perturbative calculations to construct different
  classes of symmetric periodic orbits of second species. Levi-Civita regularization and a more geometric
  approach were introduced in Henrard (1980), and following these ideas a complete proof of the existence of
  symmetric periodic orbits with two collisions per period was given in Marco and Niederman (1995). All these
  symmetric periodic orbits converge, when the mass parameter mu -> 0, either to sequences of arcs of Henon type
  S, or to two arcs of type T each of which is the symmetric image of the other."
- "The existence of a large set of periodic and chaotic orbits of second species that converge to sequences of
  arcs of type T was proved in Bolotin and Mackay (2000) for the circular problem, and a similar result was later
  reported in Font et al. (2002), where a numerical study of these in general asymmetric periodic orbits was also
  presented. The proof of Bolotin and Mackay (2000) was extended to the elliptic problem in Bolotin (2005, 2006)."
- "All the results mentioned above deal with planar orbits. The existence of nonplanar periodic and chaotic
  orbits of second species was proved only recently (Bolotin and Mackay 2006), by extending the variational
  method developed in Bolotin and Mackay (2000) to sequences of segments of noncoplanar Keplerian orbits that
  collide with the small primary at opposite ends of a straight line through the big primary."
- "In Bolotin and Mackay (2006), the possibility of enlarging the shift of planar orbits by combining arcs of
  type S and arcs of type T is also mentioned. It is formally stated as a conjecture in Henon (1997) that every
  infinite periodic sequence of S- and T- arcs which does not contain two identical T- arcs in succession is the
  limit when mu -> 0 of a family of periodic orbits of the restricted three body problem."
- "Here we give numerical evidence in support of this conjecture, and a description of the shift that corresponds
  to a large set of planar orbits of second species, that is, orbits either periodic or chaotic with infinitely
  many close encounters with the small primary, for a given value of the mass parameter mu."

(The paper makes no statement about the "small angle changes" restriction of the 2002 paper in these terms. What
it adds is the S-arcs, which begin and end in collision, and so include encounters in which the exit direction is
unrelated to the entry direction; see section 4.)

## 3. The two-body skeleton and the homoclinic families (Section 2, pp. 145-154, READ)

Unperturbed problem (mu = 0, p. 146): initial conditions (x, y) = (-1, 0), (xdot, ydot) = (v^s cos psi,
v^s sin psi) at t = 0, v^s = sqrt(3 - C_J), "psi in I = [psi*, psi**] = [arcsin((2 - C_J)/(2 sqrt(3 - C_J))),
pi - arcsin((2 - C_J)/(2 sqrt(3 - C_J)))], the interval that corresponds to elliptic orbits." Figures 1 and 2 plot
the first pericentre position and the time to it as functions of psi (C_J = 2.8, so v^s = sqrt(0.2) = 0.447,
INFERRED arithmetic). A self-intersection of the synodic first-pericentre curve (Fig. 1(a)) gives a pair of exit
angles psi_1, psi_2 = pi - psi_1 whose two-body orbits have the same pericentre in the synodic frame; there are
three such pairs (p. 147). The relation printed on p. 148:
"phi_0(psi_1) - t_1 = phi_0(psi_2) - t_2 - 2(k + 1) pi, k = 1, 2, 3, ... or, given that psi_2 = pi - psi_1,
phi_0(psi_1) = 2 pi - phi_0(psi_2), 2 phi_0(psi_2) = t_2 - t_1 + 2(k + 2) pi, k = 1, 2, 3, ...", where phi_0 is the
argument of the pericentre in sidereal coordinates and t_1, t_2 the times to reach it. Using the symmetry
(x, y, t) -> (x, -y, -t) for mu >= 0, the "1-homoclinic" arcs (exactly one intermediate pericentre passage
between two passages through (-1, 0)) come from the intersection of the pericentre curve with its reflection
(Fig. 5).

Counts (pp. 149-154):
- Eleven 1-homoclinic arcs with time of travel to or from the pericentre less than 8 pi, in two families: "One
  family of (p, q)-resonant orbits, in the notation of Font et al. (2002), or T-arcs of Henon (1997), that begin
  and end at the same point, in fixed coordinates, q being the exact number of full turns of the rotating
  reference frame, and p the number of full turns of the massless body in its elliptic orbit. For all the
  1-homoclinic orbits, p = 1, and q goes from 1 to 4 due to the time cut-off that we introduced. For each (p, q)
  pair we have the + and the - orbits, according to whether the orbit starts going inwards or outwards the unit
  circle. These 8 orbits are not symmetric ... The other family has 3 symmetric homoclinic orbits (labels 5, 7
  and 8), which can be obtained from the basic orbit of Fig. 6c by adding half turns ... The symmetric homoclinic
  orbits correspond to the S-arcs of Henon's notation."
- "All these primary 1-homoclinic orbits are associated with transversal intersections of the pericentre curve
  with its image by the y -> -y symmetry. These intersections correspond to exit angles that do not involve
  intermediate collisions with the small primary, and so they will persist as transversal intersections of the
  stable and unstable manifolds of collision for the perturbed problem." (p. 150)
- For C_J = 2.8 and p = q = 4 (p. 154): "the set of basic homoclinic arcs of the rotating two-body problem has
  41 orbits, 18 resonant orbits or T-arcs (R, p, q, s), where p and q are the integers that denote the number of
  full turns of the massless body and of the reference frame, respectively, that take place along the orbit and
  s = +1 (resp. s = -1) for the ingoing (resp. outgoing) orbits; and 23 symmetric orbits or S-arcs (S, p, q, s)."

Tables of the admissible orbits (p. 154 to 155, READ; counts agree with the 18 and 23 stated in the text):

Table 1 "Admissible resonant orbits (R)":

| p | q | s |
|---|---|---|
| 1 | 1,2,3,4 | +-1 |
| 2 | 1,3 | +-1 |
| 3 | 2,4 | +-1 |
| 4 | 3 | +-1 |

Table 2 "Admissible symmetric orbits (S)" (rows with a blank p continue the p above):

| p | q | s |
|---|---|---|
| 0 | 1,2,3 | -1 |
| 1 | 1,2,3 | -1 |
| 2 | 0,1,2,3 | +1 |
| (2) | 1,2,3 | -1 |
| 3 | 1,2,3 | +1 |
| (3) | 2,3 | -1 |
| 4 | 1,2,3 | +1 |
| (4) | 2,3 | -1 |

The meaning of p and s for the S-arcs is not separately restated beyond the sentence quoted above (so p, q, s
for S-arcs follow the same general description; INFERRED that p is then a count of pericentre passages along the
symmetric arc).

Perturbed manifolds (pp. 150-154, mu = 1e-4, C_J = 2.8): Figures 7 and 8 show the curve W1 of first pericentre
positions of orbits leaving collision, and the first pericentre distance as a function of psi. Quoted (p. 153):
"For the perturbed problem, in a very narrow range of exit angles close to psi_1 there are quasi-collision orbits
that undergo large deflections and leave again a small neighbourhood of the point (-1, 0) along all directions,
after which they follow closely a two-body orbit again until they reach the next pericentre. ... Deflected
orbits will leave collision along all directions, and psi = pi/2 is approximately the exit angle that separates
these orbits from the ones that are deflected without passing through a pericentre." And (p. 154): "The
intervals of values of psi for which these close encounters generate significant deflections are of order
mu^(2 alpha) (Font et al. 2002)." Discontinuity points A_i, A'_i of the first-pericentre curve are the places
where a deflected orbit gains an extra pericentre passage.

## 4. The return map, symbolic dynamics and the theorem (Sections 3 and 4, pp. 154-162, READ)

Homoclinic strips (p. 154 to 155): each of the 41 basic homoclinic orbits corresponds to a strip in the (phi,
psi) torus of initial conditions that leave C and return to C and in between stay close to the basic orbit; "We
shall call it the (X, p, q, s)-homoclinic strip." Figure 9 shows the 41 strips for C_J = 2.8, p, q <= 4 at
mu = 1e-6 (a), the remaining strips at mu = 1e-4 (b) and mu = 1e-3 (c). Lemma 1 (p. 155): for given C_J a
homoclinic strip is a set of (phi, psi) in the annulus {|phi - psi| < pi/2} such that psi = psi~(X, p, q, s, C_J)
+ mu^alpha (B(X, p, q, s, C_J, phi) + xi C(X, p, q, s, C_J, phi)), xi in [-xi_m, xi_m] subset [-1, 1]. "According
to the preceding result, the width of all the homoclinic strips is of order O(mu^alpha) ... mu can be taken small
enough so that all the homoclinic strips of each class are isolated ... provided that the psi~ ... are all
different for p, q < N." (p. 156)

Lemma 2 (p. 158), quoted: "For a given value of the Jacobi constant C_J, the return map on a symmetric (S, p,
q, s) homoclinic strip is of the form T o F, where F ... and T(phi, psi) = (phi + Delta(p, q, s, C_J), psi +
Delta(p, q, s, C_J)) is a translation along the line phi = psi on the torus and F can be approximated by the
expressions given by Lemmas 5.1 and 5.2 of Font et al. 2002 in their regions of validity."

Alphabet (p. 158 to 159): "For the selected value of C_J and for mu = 1e-4, the shift that describes the dynamics
of the planar orbits of second species built with arcs that verify the constraints imposed on time and number of
intermediate pericentre passages is a large subshift of the full shift on the alphabet of 41 symbols that
corresponds to the admissible two-body problem homoclinic arcs. Every transition (R, p, q, s) -> (R, p, q, s) is
forbidden, and a few others may have to be excluded. For instance, in Fig. 11a, we see that for mu = 1e-4 the
image of (R, 2, 1, -1) fails to intersect the strip (S, 1, 3, -1), but the intersection subsists for mu = 1e-6."
"The symbolic dynamics on this large subshift of the full shift implies, in particular, the existence of
periodic orbits that combine arcs of types S and T." "The approach we have followed here and the numerical
results support a complete description of the orbits of second species for the chosen parameter values. For
smaller values of mu, the behaviour of the return map is closer to the estimates given in Font et al. (2002) ...
In particular, the gap between the two horizontal branches of the image of each strip closes down to zero,
because the size of the projection on the psi axis of these horizontal branches tends to zero as mu^alpha.
Hence, as mu -> 0, the subshift of the symbolic dynamics is enlarged, in the sense that transition matrix gains
more non-zero elements" (p. 159 to 161).

Theorem 1 (p. 161), quoted exactly: "For a given upper bound N on p and q, mu can be chosen small enough so that,
for C_J outside a finite union of small intervals, the set of all the orbits of the second species is described
by a subshift on the alphabet of all the admissible S- and T-arcs of the two body problem, for which the only
forbidden transition is the one that concatenates two identical resonant homoclinic arcs."

Remark after Theorem 1 (p. 161 to 162): "Note that, for a fixed C_J, when mu -> 0 the value of N can increase and,
hence, the set of admissible basic homoclinic orbits contains more elements and the alphabet for the subshift
has more symbols. On the other hand, in contrast with similar results on embedding of a subshift in Celestial
Mechanics (see, e.g. Simo and Martinez 1988), instead of having passages near a few homoclinic orbits which
require an additional label to denote some number of revolutions between two successive passages, in the
present problem the number of basic homoclinic orbits is large and increases when mu -> 0."

Proof sketch (p. 162): by Lemma 2 and the 2002 results, "it remains to show that for a given upper bound N on p
and q and mu small enough, the image by the return map of a symmetric homoclinic strip intersects all the
homoclinic strips, and the image of a resonant strip intersects all the other homoclinic strips. The latter
assertion was proved in Font et al. (2002) and it follows from the behaviour of the IN-map, provided that
resonant orbits with intermediate close encounters can be avoided. As argued in Bolotin and Mackay (2000), Font
et al. (2002), this can be achieved for mu small enough by choosing C_J outside a finite union of small
intervals ... the image of any strip intersects any horizontal strip except those inside a neighbourhood of
order mu^(1-2 alpha) of the image of the strip by the OUT-map of the strip. Therefore, the former assertion, that
the image by the return map of a symmetric homoclinic strip intersects all the homoclinic strips, requires that
the projections on the psi axis of the image of the symmetric strip by the OUT-map and of the remaining
homoclinic strips have non-overlapping neighbourhoods of order mu^(1-2 alpha). Again, this can be avoided by
tuning C_J." The paper calls this "a sketch of the proof".

## 5. Numerical data (READ)

### 5.1 What was computed, for which mass ratios and Jacobi constants

- Jacobi constant: C_J = 2.8 throughout (p. 145: "We have considered C_J = 2.8 and q = p = 4 throughout").
  No other value is used. (C_J = 3 is the threshold of the Hill-type region; 2.8 gives v^s = sqrt(3 - C_J) =
  0.447 at distance 1 from the large primary at mu = 0, INFERRED arithmetic.)
- Mass parameter: mu = 1e-4 (most examples, all periodic orbits), mu = 1e-6 (comparison, Figures 9a, 11b, 13), and
  mu = 1e-3 (Figure 9c only). The maximum mu in the paper is 1e-3. The Earth-Moon value (about 0.0122) and any
  value above 1e-3 are not reached. In the 2002 companion paper the numerics reach 1.5e-3 (merging of strips).
- Encounter circle radius mu^alpha with alpha = 0.4 (p. 146).
- Time and passage bounds: consecutive encounters within q <= 4 turns of the frame and p <= 4 pericentre
  passages. Integration of the manifolds to time 8 pi (p. 151).
- Method: the strips are computed by the numerical procedure of the 2002 paper (Levi-Civita regularised
  integration; shown in Fig. 9 as "found numerically with the method used in Font et al. 2002"); periodic orbits
  are found by "the symbolic dynamics" (return-map fixed points on the strips, p. 158); no continuation in mu of
  a periodic orbit is described (READ, by absence: the paper does not follow any periodic orbit in mu; only the
  strips are compared at mu = 1e-6, 1e-4, 1e-3).

### 5.2 Comparison with the collision-orbit (two-body) limit as mu grows

- Fig. 7 (mu = 1e-4) against Fig. 1 (mu = 0): "The perturbed curve almost overlaps the two-body problem curve, and
  has three additional double branches branching out (for mu > 0) from the self-intersection points." (p. 153)
- Fig. 9: the homoclinic strips are thin and isolated at mu = 1e-6, thicker and overlapping at 1e-4, and heavily
  deformed at 1e-3 (figure only; no number is printed).
- Fig. 11: at mu = 1e-4 the image of (R,2,1,-1) misses the strip (S,1,3,-1), at 1e-6 it crosses it (the example
  quoted in section 4).
- Fig. 13: the strip (S,2,1,-1) and its return-map image at mu = 1e-6 and 1e-4: the image has two nearly parallel
  horizontal branches whose gap closes as mu -> 0.
- No position error, no table of orbit versus patched-conic distance, and no comparison of a periodic orbit at
  mu = 1e-3 with its collision-orbit limit is printed. The statement used is qualitative: transitions are lost as
  mu grows, and the shift enlarges as mu -> 0.

### 5.3 Printed periodic-orbit data (all mu = 1e-4, C_J = 2.8)

Tables 3 and 4, transcribed digit by digit from the 400 dpi image (no entry computed; all digits as printed;
(phi, psi) are the initial angles on the circle of radius mu^0.4 as in section 1). The sequences are the
concatenations of homoclinic arcs in time order; the "Period" is the period of the periodic orbit in the
paper's time units (INFERRED: nondimensional time in which the primaries' period is 2 pi).

Table 3 "Examples of periodic orbits for mu = 1e-4 and C_J = 2.8":

| | Figure 12a | Figure 12c | Figure 12e |
|---|---|---|---|
| (X,p,q,s) | (S,1,1,-1) | (S,1,2,-1), (S,3,1,1) | (S,2,1,1), (R,3,2,-1) |
| phi | 2.685003594282268 | 3.0141823898790555595 | 2.5457836425596942403 |
| psi | 2.639352410113041 | 3.0316423709511530133 | 2.4204845762785093917 |
| Period | 7.933918152222289 | 27.003867331650326275 | 22.310316405350279103 |
| Stability parameter | 0.2639321981E+05 | 0.1946971773492E+08 | -0.1201986472185E+09 |

Table 4 "Examples of periodic orbits for mu = 1e-4 and C_J = 2.8":

| | Figure 12g | Figure 12i |
|---|---|---|
| (X,p,q,s) | (R,2,1,1), (S,2,3,1), (R,1,3,1) | (R,1,2,-1), (S,1,2,-1), (R,1,2,1) |
| phi | 0.98857090907035220406767336927899 | 3.27824276059703769776947746511 |
| psi | 0.954760770266642116216411782118 | 3.21924695055804868387324783072 |
| Period | 45.8085897638589254831031186215 | 40.8236298146304851361874699996 |
| Stability parameter | 0.2867317722679597242E+14 | 0.9608334138976017932E+11 |

Notes on the tables: the stability parameter is not defined in the text I read (pp. 143 to 162); its magnitude
(1e4 to 1e14) is the printed expected sensitivity. Figures 12b, d, f, h, j are close-ups of the orbit passing
through the circle C of radius mu^alpha (the arcs are numbered 0, 1, ... in time order in the panels e to j,
Figure 12 caption). Figure 12 has five orbits (a/b, c/d, e/f, g/h, i/j) and Tables 3 and 4 give one column per
orbit; the digit strings of the (phi, psi) of Table 4 run to 30 significant decimals, beyond double precision
(INFERRED: the authors used extended precision; the 2002 paper says "We use extended precision to compute the
periodic orbits and the stability parameter").

From the 2002 companion (printed in that digest): the periodic orbits of 2002 Figures 8, 9, 10 at the same mu and
C_J, with (phi, psi) and stability parameters but no periods.

## 6. Statements on more than one small body, elliptic problem, time dependence, spacecraft

- More than one small body: none (READ; the entire paper is the three-body problem).
- Elliptic problem: only through the citations of the elliptic extension "Bolotin (2005, 2006)" and "Gomez and
  Olle (1991)... extended to the elliptic case" (p. 144, quoted in section 2); nothing is computed.
- Time-dependent problem: none beyond that.
- Spacecraft applications: none in this paper. The only related sentence is in the 2002 paper's conclusion
  (successive flybys of the small primary by a spacecraft "can help put the spacecraft on a desired nominal
  orbit around the largest primary at low cost"). The 2009 paper contains no application sentence (READ: the
  paper ends with the proof sketch and references).

## 7. What this does and does not cover for #890

(a) One small secondary. The complete family of two-body skeleton arcs (41 at C_J = 2.8 with p, q <= 4: 18 resonant
T-arcs and 23 symmetric S-arcs) and a numerical symbolic dynamics on them at mu = 1e-4, with a theorem (Theorem 1,
stated for mu small enough, C_J outside a finite union of small intervals, a stated limit on p, q <= N), a proof
sketch only. Computed up to mu = 1e-3 (strips shown, none of the periodic orbits); periodic orbits only at
mu = 1e-4. Nothing at mu of order 1e-2.

(b) Two secondaries with different periods: nothing is said.

(c) Reproducible in `src/cyclerfinder/core/cr3bp.py` as checks of the project's flyby-chain continuation
(INFERRED recipes, not tested): (i) the five periodic orbits of Tables 3 and 4: set mu = 1e-4, C_J = 2.8,
alpha = 0.4; state x = mu - 1 + mu^alpha cos(phi), y = mu^alpha sin(phi), speed from the Jacobi constant,
direction psi (convention INFERRED from the leaving condition |phi - psi| < pi/2, to be tested on the shortest
orbit first, Figure 12a: (S,1,1,-1), one arc, period 7.933918152222289); the pass is closure at the printed
period; the printed digits go beyond double precision, so closure to the integrator's precision is the test,
and the orbits are very unstable (1e4 to 1e14), so the shortest one first. (ii) The mu = 0 skeleton: enumerate
the two-body homoclinic arcs at C_J = 2.8 with p, q <= 4 and check the counts against 41 = 18 + 23 and Tables 1
and 2; this uses only Kepler solutions and tests the project's patched-conic enumeration logic, which is the
closest analogue of the #888 enumeration. (iii) A continuation of the strips in mu from 1e-6 to 1e-3 and the
loss of the (R,2,1,-1) to (S,1,3,-1) transition between 1e-6 and 1e-4 (Figure 11) as a check on where
continued orbits stop connecting. Not reproducible from the print: strip positions at 1e-6 or 1e-3 (figures
only), and any periodic orbit away from mu = 1e-4.

(d) Wording supported: "For the planar circular restricted three-body problem with one small secondary, Font,
Nunes and Simo (2009) computed, at C_J = 2.8 and mu = 1e-4 (comparison at 1e-6 and 1e-3), the full set of
orbits with consecutive close encounters of bounded time and pericentre count, identified them with S- and T-arcs
of the two-body problem, gave numerical evidence for the symbolic dynamics, a proof sketch, and printed five
periodic orbits that combine S and T arcs." Supported also: encounters of the one-moon problem include
arbitrary exit directions after a near-collision passage (S-arcs), so large turns at a close flyby are part of
the one-moon second species picture at this paper's level (it is the S-arcs that give them; the 2002 theorem
handles only the T-type resonant arcs, with O(mu^alpha) deflection). Not supported: any statement about two
secondaries, any mu above 1e-3 (the largest printed), any claim about the real Uranus system, and any claim that
a periodic orbit was continued in mu. The #890 object (two moons, mu of order 4e-5 each INFERRED from the note's
GM values and a Uranus GM of about 5.8e6 from memory) lies inside the mass-ratio range these papers cover for ONE
moon (1e-6 to 1e-3), but with a second moon the setting is outside them.

## 8. References cited (p. 162, READ, full list as printed)

Bolotin, S.V.: Second species periodic orbits of the elliptic 3 body problem. Celest. Mech. Dyn. Astron. 93,
343-371 (2005) [printed pages 343-371; other papers give 343-373].
Bolotin, S.V.: Shadowing chains of collision orbits. Discrete Contin. Dyn. Syst. 14, 235-260 (2006).
Bolotin, S.V., Mackay, R.S.: Nonplanar second species periodic and chaotic trajectories for the circular
restricted three body problem. Celest. Mech. Dyn. Astron. 94, 433-449 (2006).
Bolotin, S.V., Mackay, R.S.: Periodic and chaotic trajectories of the second species for the n-centre problem.
Celest. Mech. Dyn. Astron. 77, 49-75 (2000).
Bruno, A.D.: On periodic flybys to the Moon. Celest. Mech. 24, 255-268 (1981).
Font, J., Nunes, A., Simo, C.: Consecutive quasi-collisions in the planar circular RTBP. Nonlinearity 15,
115-142 (2002).
Gomez, G., Olle, M.: Second species solutions in the circular and elliptic restricted three body problem,
I and II. Celest. Mech. Dyn. Astron. 52, 107-146 and 147-166 (1991).
Guillaume, P.: Linear analysis of one type of second species solution. Celest. Mech. 11, 213-254 (1975a).
Guillaume, P.: An extension of Breakweel-Perko's matching theory. Celest. Mech. 11, 449-467 (1975b).
Henon, M.: Sur les orbites interplanetaires qui rencontrent deux fois la Terre. Bull. Astron. 3, 377-393 (1968).
Henon, M.: Generating Families in the Restricted Three-Body Problem. Springer-Verlag, Berlin Heidelberg (1997).
Henrard, J.: On Poincare's second species solutions. Celest. Mech. 21, 83-97 (1980).
Marco, J.P., Niederman, L.: Sur la construction des solutions de seconde espece dans le probleme plan
restreint des trois corps. Ann. Inst. H. Poincare Phys. Theor. 62, 211-249 (1995).
Perko, L.: Second species solutions with an O(mu) near-moon passage. Celest. Mech. 14, 395-427 (1976).
Poincare, H.: Les Methodes Nouvelles de la Mecanique Celeste. Tome I. Gauthier-Villars (1892).
Poincare, H.: Les Methodes Nouvelles de la Mecanique Celeste. Tome III. Gauthier-Villars (1899).
Simo, C., Martinez, R.: Qualitative study of the planar isosceles three body problem. Celest. Mech. 41,
179-251 (1988).
Stiefel, E.L., Scheifele, G.: Linear and Regular Celestial Mechanics. Springer-Verlag, Berlin Heidelberg (1971).
Szebehely, V.: Theory of Orbits. Academic Press, New York (1967).

Note: Henon 1968 is listed in this paper as "Bull. Astron. 3, 377-393 (1968)" and in the 2002 paper as "Bull.
Astronomique Paris 1, 377-402 (1966)"; not resolved here. The two Henon papers on interplanetary orbits are
presumably the same work cited with different metadata (INFERRED).

## Addendum 2026-10-04 (`#896`): corrections found when the printed numbers were made tests

Test: `tests/core/test_cr3bp_font_nunes_simo_second_species.py` (commits `a18cf2bc`, `553daa8b`).
All five orbits of Tables 3 and 4 close in `core.cr3bp`.

- Section 5.3, Table 4, Fig. 12g phi: the digest gives 0.98857090907035220406767336927899 (32
  decimals). The table (p155, re-read at 400 dpi) prints 0.988570907035220406767336927899 (30
  decimals, like its psi); the digest has an inserted "09". The digested value misses the fixed point
  of the return map by 2.0e-9; the printed value closes to 1e-13. The other entries of Tables 3 and 4
  are correct.
- Section 7(c), the recipe for a state on the encounter circle, is in the PAPER's frame: big primary
  at (mu, 0), small primary at (mu - 1, 0) (eq. (1.1), p145). `core.cr3bp` has them at (-mu, 0) and
  (1 - mu, 0) with the same Coriolis sign, so the paper's frame is the project's turned by pi; every
  angle shifts by pi. Used literally in `core.cr3bp` the recipe never returns to the circle (a control
  in the test).
- The Jacobi constant of the paper includes the constant mu(1 - mu) (eq. (1.2), p145; Omega has
  + mu(1 - mu)/2). `core.cr3bp.jacobi_constant` omits it, so C_J = 2.8 is 2.79990001 in the project.
  Using 2.8 directly makes orbit 12a miss by 0.33.
- The convention is printed, not inferred: psi is "the angle between the x-axis and the synodic
  velocity" (p151; also (xdot, ydot) = (v cos psi, v sin psi) on p146); phi is counter-clockwise from
  +x on the circle of radius mu^(2/5) about the small primary, and the points are exits from the disk
  (p145).
- The "stability parameter" is not defined in either paper; it equals the trace of the linearised
  return map to 7e-10 .. 3.8e-7 (INFERRED definition).
