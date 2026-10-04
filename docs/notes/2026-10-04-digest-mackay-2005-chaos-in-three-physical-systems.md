# Digest: MacKay 2005, chaos in three physical systems (second-species subshifts, Anosov linkage, ergodic pumping)

Date 2026-10-04. Purpose: this is the paper that the 2006 Bolotin-MacKay paper cites for "a summary, some minor
additions and a correction" of the 2000 second-species paper, and that the project's digest of the 2000 paper
(`docs/notes/2026-10-04-digest-bolotin-mackay-2000-second-species-n-centre.md`, section 5, item 5) could not read.
It is now read. Only Section 1 (about three and a half pages of fourteen) concerns the project: the statement of
the planar circular restricted result, a short sketch of the proof, one footnote that contains the correction of
the direction-change lemma, and one paragraph on instability and "controllability" for flyby missions. Sections 2
and 3 concern unrelated systems and are summarised only for completeness. Related digests, not repeated here:
the 2000 and 2006 Bolotin-MacKay digests, `docs/notes/2026-10-04-digest-dimare-2010-quasi-collision-3-centre.md`,
`docs/notes/2026-10-04-digest-marco-niederman-1995-seconde-espece-plan-restreint.md`, and the companion digest of
Bolotin and Negrini (2001), `docs/notes/2026-10-04-digest-bolotin-negrini-2001-regularization-entropy-n-center.md`.

Citation: R. S. MacKay, "Chaos in three physical systems", in Equadiff 2003 (World Scientific, 2005, pp. 59-72, DOI
10.1142/9789812702067_0006. The held file is the author's PostScript converted to PDF (running head "December 6,
2003 14:51 WSPC/Trim Size: 9in x 6in for Proceedings 3chaos"), 14 pages; PDF page n is journal page n + 58 (so
journal pages 59 to 72). Filed in the private paper corpus as
mackay-2005-chaos-in-three-physical-systems-equadiff-2003-doi-10.1142-9789812702067_0006-author-postscript-converted.pdf.
Text layer: present but defective, about 6,250 words by `pdftotext | wc -w`; the conversion drops the letters c
and the Greek letters (for example "ex ept", "dire tion", "Ja obi", and every Omega and epsilon vanishes), so
every statement below was read from the page images, and the text layer was used only to confirm that the page
images and the text agree where the text survives. The text layer must not be used for the formulas.

Evidence labels: READ means read from the printed page in this session. INFERRED means my reading of what the
printed text implies, or a comparison with project work. All 14 pages were read as page images. Page numbers
below are the paper's own page numbers 1 to 14 (journal page minus 58). Quotation marks mean verbatim;
"not equal" is written in words.

## 0. Abstract (p. 1, READ)

"Examples of various types of chaos are described in some physical systems: (1) Subshifts of second species
orbits in the circular restricted three body problem (with S.V. Bolotin); (2) Anosov energy levels for the
frictionless dynamics of a mechanical linkage, and a normally hyperbolic Anosov submanifold for weak friction and
suitable feedback control (with T.J. Hunt); (3) "Ergodic pumping": a proposed mechanism for the power stroke of
myosin, and a design principle for nanobiotechnology (with D.J.C. MacKay)."

## 1. Section 1, second-species subshifts in the three-body problem (pp. 1-4, READ)

### 1.1 The history paragraph (p. 1)

"The idea that a deterministic system could exhibit chaotic motion arose during Poincare's work on the three body
problem, but he came across it in a near integrable part of the phase space which was hard to analyse because the
resulting chaos is exponentially weak in the relevant parameter. In a less well known part of his work^1, Poincare
studied orbits of the three body problem with close encounters. 100 years later we proved that this provides a much
more prominent source of chaos^5, as I will review here."

"Poincare's suggestion was that the three body problem of celestial mechanics contains "second species" periodic
orbits, i.e. ones which, as the two smaller masses tend to zero, tend to sequences of arcs of Kepler orbits about
the heavy mass, joined at collisions between the smaller bodies. It is agreed^2 that he did not have a valid proof
and that his claim is false in the generality he stated. Numerical computations^3, however, suggest many such
orbits exist, in particular in the planar circular restricted case, where I shall refer to the bodies as the Sun,
Jupiter and an Asteroid. Marco and Niederman^4 provided a proof of existence of a second species periodic orbit,
more precisely of a one-parameter family with two close encounters of the Asteroid with Jupiter per period for
small enough mass of Jupiter. Bolotin and I^5 extended this dramatically, proving existence of many one-parameter
families of chaotic subshifts of second species orbits in the same case."

Reference numbering as printed: 1 Poincare 1899 chapter XXXII; 2 Levy's notes p. 629 of Poincare Oeuvres VII
(1952); 3 Henon, Generating families in the restricted three body problem (1997); 4 Marco and Niederman, Ann.
IHP Phys. Theor. 62 (1995) 211-249; 5 Bolotin and MacKay, CMDA 77 (2000) 49-75.

### 1.2 Theorem 1.1, quoted in full (p. 2)

"Theorem 1.1.^5 In the planar circular restricted three body problem with masses 1 - epsilon, epsilon, 0, for all
values C in the interval (-sqrt(8), 3) there are arbitrarily large subsets T of rationals and epsilon_0(T) > 0 such
that for any sequence (m_n/k_n), n in Z, in T and epsilon in (0, epsilon_0) there is a unique trajectory of Jacobi
constant C near to a chain of collision trajectories formed by transforming to the rotating frame Kepler ellipses
of frequencies m_n/k_n traversed m_n times, and it converges to the chain as epsilon -> 0."

Differences from the 2000 paper's Theorem 2.1 as digested (READ for both): the 2000 statement has a dense subset
S_C of rationals in the allowed frequency set A_C and then "for all finite subsets T subset S_C"; this statement
has "arbitrarily large subsets T of rationals". The set T is chosen by the authors (INFERRED: so that no early collision
occurs, as in the 2000 Lemma 3.1), so it is not all rationals in A_C. The statement here has no exception clause;
the exception is in footnote a (section 1.5 below). The statement says nothing about uniformity of epsilon_0
beyond its dependence on T, nothing about closest approach (the 2000 theorem's inequality (1.5), distance of at
least c epsilon from collision) and no value for epsilon_0.

### 1.3 The remark on instability and flyby missions (p. 2)

"The result is hyperbolic horseshoes, which occupy a set of measure zero and are unstable, so their physical
relevance is not immediate. Horseshoes like these are, however, the principle of solar system exploration
missions like Voyager (though our result is much more restricted because our horseshoes make close encounters
with only one planet and assume it is in a circular orbit). The angle change in a close encounter is a smooth
function of impact parameter with slope of order 1/epsilon, so these orbits have Lyapunov exponent of the order
of log 1/epsilon (in units of the frequency of Jupiter). This is a much stronger instability, and hence offers
greater controllability, than schemes using unstable periodic orbits around the colinear Lagrange points, for
example, and Sitnikov's horseshoe."

The word "strong controllability" is not used on this page; that phrase is from the 2006 paper (p. 445 there, per
the 2006 digest). The statement about spacecraft is the sentence above: flyby tours are "the principle" of such
horseshoes, restricted to one planet in a circular orbit, and the strong instability "offers greater
controllability". The mechanism is stated: the deflection angle depends smoothly on the impact parameter with
slope of order 1/epsilon, which is the origin of the log(1/epsilon) exponent. No numerical example, no mass ratio
and no cost estimate is given.

### 1.4 The proof sketch (pp. 2-3)

"Here is the idea of our proof. In the frame rotating with Jupiter about the Sun, the Jacobi invariant

  C = -|q'|^2 + |q|^2 + 2(1 - epsilon)/|q| + 2 epsilon/|q - j| - 2 epsilon x   (1)

for the Asteroid is preserved, where q = (x, y) is its position and Jupiter is at j = (1, 0)." (q' is the
rotating-frame velocity; I transcribe the first term as printed, with the sign convention that makes C equal to
minus twice the energy; this matches the 2000 paper's Section 2 as digested.)

"Collisions with Jupiter can be regularised by taking a square root coordinate change about j and rescaling time,
to convert the orbits of the Asteroid to those of energy epsilon of a new C-dependent Hamiltonian which is smooth
at j and has a saddle over it at energy 0. The saddle has rotationally symmetric quadratic part. Each segment of
Kepler orbit for the Asteroid with epsilon = 0 that happens to start and end on Jupiter gives a homoclinic orbit
to the saddle. In particular, this holds for each rational frequency m/k for which there is no intermediate
collision during k revolutions of Jupiter's orbit and m revolutions of the Asteroid's. Figure 1 shows two such
ellipses. Each ellipse produces two homoclinic orbits of the regularised system because of the square-root
cover. A large subset of the homoclinic orbits are non-degenerate (including all those which make a complete
number of revolutions), i.e. they are non-degenerate critical points of the corresponding Jacobi action
functional. Non-degenerate homoclinic orbits to a saddle with rotationally symmetric quadratic part can be glued
together to produce nearby true orbits for the regularised system for epsilon small and positive, if the angle
changes are less than pi/2. We proved this by taking a small ball around the saddle and defining a discrete-time
action functional for successive entries and exits of the ball using the Jacobi action and then using the
implicit function theorem to continue the sequence of collision orbits with respect to epsilon (this is a
perturbation from an "anti-integrable limit"). The resulting orbits of the regularised system convert to true
orbits of the Asteroid in the rotating frame if the angle changes of the sequence of Kepler orbits there are
never 0 or pi (the first keeps the angle changes for the regularised system less than pi/2 if the appropriate
branch of the cover is taken, and the second keeps the Asteroid away from collisions with Jupiter).^a"

Reading notes (INFERRED, the second only from the printed sentence): (1) "the angle changes of the sequence of
Kepler orbits there" are the angle changes in the rotating frame ("true orbits of the Asteroid in the rotating
frame"), so the condition of the theorem is stated in the rotating frame already in this text; (2) the square-root
(Levi-Civita type) cover halves angles at the saddle, which is why a physical angle change strictly between 0
and pi corresponds to a regularised angle change below pi/2: the physical reversal (pi) is the excluded boundary
and the physical no-turn (0) is excluded by the first clause of the parenthesis; this is consistent with, but not
spelled out in, the printed parenthesis.

Figure 1 (p. 3), caption verbatim: "The two orbits of an Asteroid with frequency 1 and given Jacobi constant C in
(1, 3), starting and ending at a zero mass Jupiter, in the inertial frame." The picture: a dashed unit circle
(Jupiter's orbit) and two ellipses about the Sun that share one point on the unit circle (I read the picture as two
mirror-image ellipses, as the caption's "two orbits ... with frequency 1 and given Jacobi constant" implies; I
did not rely on finer detail of the drawing). It is the only figure of Section 1, and it is drawn in the INERTIAL frame,
the frame in which the original direction-change lemma was (wrongly) proved.

There is no graph of collision arcs in this paper, no vertex and edge language, no mention of positive
topological entropy, and no statement of the topological Markov chain. The structure is "any sequence in T", with
T a set of frequencies. For the graph formulation (vertices K, edges G, branched subgraph) see the 2000 digest,
section 1.3.

### 1.5 THE CORRECTION: footnote a (p. 3), quoted in full

"^a I take the opportunity to correct an omission from Ref. 5, which does not, however, affect the truth of the
main result. In Lemma 3.3 we claimed that for all sequences Omega_n of rationals in a set S_C there is an 'orbite
a chocs' consisting of ellipses of frequency Omega_n and Jacobi constant C connecting collisions at Jupiter which
leave each collision in neither the same nor the opposite direction as they arrive. We proved this in the
inertial frame, but neglected to check the result holds in the rotating frame, which is where we need it to prove
our Theorem 2.1. In fact, it holds in the rotating frame for all such sequences except those containing two
consecutive Omega_n satisfying Omega^{2/3} + 2 Omega = C. The obstacle in that case is that the two ellipses with
such a frequency have angular momentum +1, so they transform to orbits in the rotating frame which both leave and
arrive at Jupiter tangent to the radial direction. For all other frequencies the transformation of the two
ellipses to the rotating frame results in two orbits that leave Jupiter in different non-radial directions,
because the angular momentum is (C - Omega^{2/3})/2 Omega."

The corrected condition, therefore (READ): the direction-change condition (arrival and departure velocities
neither the same nor opposite) must be checked in the ROTATING frame, the frame in which the small body (Jupiter)
is stationary, that is, with velocities relative to Jupiter. In that frame it holds for every sequence of
frequencies in S_C except those that contain two consecutive entries Omega_n satisfying Omega^{2/3} + 2 Omega =
C. The 2006 paper (Section 4.5, p. 444, per the 2006 digest) states the same condition in the same words: "the
velocity in the rotating frame just after each collision [must be] neither parallel nor opposite to the velocity
just before the collision". MacKay says the omission does not affect the truth of the main result; the planar
theorem therefore stands, with its hypothesis on the sequences corrected. What the corrected planar theorem says
is Theorem 1.1 above together with this footnote. There is no separately restated corrected theorem in the paper.

My own checks of the footnote (INFERRED, to be treated as my reading and not as the paper's text):

1. Why only the case of angular momentum +1 fails. With the 2000 paper's relation C = a^{-1} + 2 G (G the
   signed inertial angular momentum, prograde positive, mean motion Omega = a^{-3/2} so that a^{-1} =
   Omega^{2/3}), at the collision point |q| = 1 the velocity relative to Jupiter is (v_r, v_t - 1) in radial and
   tangential components, because Jupiter's own velocity is purely tangential with speed 1. This is radial exactly
   when v_t = 1, that is when G = 1. The two ellipses of the same (a, C) leaving the same point are mirror images
   with the same G and opposite v_r, so for G = 1 both have purely radial relative velocity. Two consecutive
   arcs of that type then have arrival and departure velocities in the rotating frame that are parallel or
   antiparallel, which is the excluded case; for G not equal to 1 the two departure velocities (plus or minus
   v_r, G - 1) are different and non-radial, so at least one of the two choices differs from the arrival direction
   and from its reverse. That agrees with the footnote's reasoning.
2. A discrepancy to resolve before using the printed formula. By the relation just used, G = 1 is equivalent to C
   = Omega^{2/3} + 2, with no factor Omega on the second term, and G = (C - Omega^{2/3})/2, again with no factor
   Omega in the denominator. The footnote prints "Omega^{2/3} + 2 Omega = C" and "(C - Omega^{2/3})/2 Omega". The
   two agree when Omega = 1 (the case drawn in Figure 1) and differ otherwise. I read the page image twice and
   the printed form has the factor Omega in both places. Either the footnote measures angular momentum in a
   different normalisation (for instance per unit of the mean motion), or this is a slip in the printed formula;
   I cannot decide which from the paper. Consequence for the project: do not implement the printed exceptional
   condition. Implement the physical statement (the relative velocity at the encounter is radial, equivalently
   v_t equals the moon's orbital speed) and test it directly on the computed vectors. Also note that Omega^{2/3} +
   2 Omega is increasing in Omega (for either form), so for a given C at most one frequency is exceptional.
3. Since at most one frequency is exceptional for given C, and T is any large set of rationals chosen by the
   authors, avoiding it costs nothing, which is consistent with the author's statement that the correction does
   not affect the truth of the main result (INFERRED).

### 1.6 Extensions and open questions (p. 4, READ)

The whole of the statement of extensions in this paper is one sentence at the end of Section 1 (p. 4): "Some
extensions may be possible to the spatial, elliptic, unrestricted and N-body cases." Spatial is done in the 2006
paper (per the 2006 digest). Elliptic means the elliptic restricted problem (time-periodic); "unrestricted" means
the full three-body problem in which both small masses are positive; "N-body" means more than three bodies. The
paper lists no other open question for the second-species problem, no list of open problems, no mention of two
small masses with different periods, no mention of a second Jupiter, and no mention of a numerical value of
epsilon_0. This exhausts what the paper says on extensions (searched across all 14 pages).

## 2. The other two systems (pp. 4-14, READ, brief)

Section 2, Anosov behaviour in the triple linkage (pp. 4-7; with T. J. Hunt): the triple linkage of Thurston and
Weeks (three disks pivoted at points fixed in a roughly equilateral triangle, each connected by a rod from an
eccentric point to a common floating pivot) has a configuration space that is a surface of genus 3. Theorem 2.1
(p. 5, citing Hunt and MacKay, Nonlinearity 16 (2003) 1499-1510): "For a non-empty open set of lengths and mass
distribution, the geodesic flow of the triple linkage is Anosov." The proof uses D_6 symmetry, the limit of
small rod length, and the sectional curvature formula (2) (negative or zero on the Schwarz P-surface, so Anosov by
structural stability). Corollary 2.1: mixing, exponential decay of correlations of Holder functions, diffusive
motion in Abelian covers. Weak friction with feedback torques Gamma_j = -gamma(K - E) theta_j-dot (gamma greater
than sqrt(I/2E)) gives a normally hyperbolic Anosov submanifold. Section 2 also discusses why smooth physical
Anosov examples are rare.

Section 3, ergodic pumping (pp. 7-13; with D. J. C. MacKay): a proposed mechanism for the power stroke of myosin,
in which entropy increase of fast variables does mechanical work on slow ones. It derives an effective
Fokker-Planck equation (6) for slow variables with diffusion tensor D, shows that free energy decreases along
paths except where the diffusive source dominates, and presents a toy Metropolis simulation of a piston in a
wedge (Figure 4) delivering about 16 kT work from an entropy increase of about 24 against a free energy budget of
about 25 kT per ATP in muscle. This part is a speculative biophysical model and is not related to the project.

Nothing in Sections 2 and 3 bears on orbital mechanics.

## 3. Open questions listed for the second-species problem (READ)

None, beyond the one sentence on extensions quoted in section 1.6. The paper is a summary talk; it has no
section of open problems. The open questions for second-species theory that are printed elsewhere are in the 2006
paper (infinitely many arcs and unbounded orbits, class 2 arcs, combined planar and nonplanar subshifts; see the
2006 digest) and are not repeated here.

## 4. What this paper does and does not cover for the two-moon object (printed text only)

Covered: the planar circular restricted problem with ONE small secondary on a circular orbit; Kepler arcs about the
primary beginning and ending at the secondary, of rational frequency m/k with whole revolutions; a common Jacobi
constant C in (-sqrt(8), 3); direction change in the rotating frame; existence and uniqueness of nearby orbits
for small epsilon; Lyapunov exponent of order log(1/epsilon), with the mechanism stated as angle change with slope
1/epsilon in impact parameter.

Not covered (all hypotheses are printed in Theorem 1.1 and the proof sketch): any second small body; different
periods; any non-circular orbit of the secondary; any time-dependent problem; any energy other than a single
Jacobi constant. The single sentence on extensions names "elliptic, unrestricted and N-body" and not a second
small body on a different circular orbit. The restricted problem here is autonomous in the frame rotating with
the one secondary; with two secondaries of different periods there is no such frame (INFERRED, as in the 2000
digest). MacKay himself notes the one-planet restriction when he speaks of Voyager: "our horseshoes make close
encounters with only one planet and assume it is in a circular orbit" (p. 2).

Safe wording (INFERRED from READ): "MacKay (2005) summarises the Bolotin-MacKay theorem for one small secondary on
a circular orbit, corrects the direction-change lemma of the 2000 paper (the condition must be checked with
velocities relative to the secondary, in the rotating frame, and fails only for a sequence containing two
consecutive arcs of one exceptional frequency), and says that extensions to the spatial, elliptic, unrestricted
and N-body cases may be possible. It describes the mechanism as a close-encounter deflection with slope of order
1/epsilon in impact parameter, giving Lyapunov exponents of order log(1/epsilon). It does not treat two
secondaries with different periods."

## 5. Techniques applicable to the project's problems

Written for someone who will build a search. READ means the idea is stated in the papers; INFERRED means it is my
construction from them and has not been tested. Where I refer to project code I have read it. The earlier digests
give the printed theorems and I do not repeat them.

### (a) The graph of collision arcs as an enumeration method for flyby chains

The construction (READ in the 2000 paper section 1.2 and 1.3, and in this paper section 1.4; the transcription
to the project's setting is INFERRED):

- Vertices: zero-mass-limit collision arcs. For one small body: Kepler arcs about the primary that start and end
  exactly on the small body, at one common Jacobi constant C, whole numbers of revolutions (class 1 in the 2006
  paper's terms), with no intermediate collision ("early collision", 2000 Lemma 3.1). For the project: Kepler
  (Lambert) arcs about Uranus between the position of a moon at one epoch and the position of a moon at a later
  epoch, with a chosen revolution count and branch.
- Edges: arc k to arc l allowed when the arrival body of k is the departure body of l (beta_k = alpha_l), the
  velocity is continuous in the sense below, and the direction changes. Positive entropy follows if the graph has a
  connected branched subgraph. A periodic chain, a cycle in the graph, gives a periodic orbit by uniqueness.
- Admissibility tests, and how to test each numerically for Kepler arcs between moon encounters:
  1. Same energy. For one secondary this is equality of the Jacobi constant C of all arcs; for Kepler arcs this
     means the same C = a^{-1} + 2 G (in the moon's units, G signed) for each arc; test the closed form, or
     compute C from the state relative to the moon. For two moons (INFERRED) there is no common conserved C; the
     patched-conic analogue of the same-energy condition is that at each flyby the magnitude of the velocity
     relative to the moon is the same before and after (the unpowered flyby conserves it); the project already
     has this as the magnitude-mismatch part of its gate. Each moon then has its own approximate invariant (a
     Tisserand-type parameter with respect to that moon) conserved across that moon's flybys but changed by the
     other moon's flybys.
  2. Nondegeneracy. The 2000 paper gives five equivalent forms. The most direct numerical one (INFERRED
     transcription for the planar case): consider the map from (departure direction angle lambda at fixed C, time of
     flight tau) to the arrival position in the plane; the arc is nondegenerate if the 2 by 2 Jacobian of
     position with respect to (lambda, tau) is nonsingular at the arc, equivalently the derivative of the arrival
     point with respect to lambda is not parallel to the arrival velocity (MacKay: "fortunate that we want
     nondegeneracy for fixed Jacobi constant rather than fixed time", 2000 digest section 2.2). In the
     time-dependent two-moon case one would instead ask that the Jacobian of the vector of flyby-condition
     residuals with respect to the vector of event epochs be nonsingular (INFERRED; see (d)).
  3. Early collision. Along each arc compute the minimum distance to every moon at the moon's position at that
     time; reject arcs that pass within the moon's sphere of influence or collision radius at an intermediate
     time that is not intended as an encounter.
  4. Direction change in the correct frame (the point of the correction). At each encounter form the arrival
     velocity and the departure velocity RELATIVE TO THE MOON (v_spacecraft minus v_moon at that epoch), compute
     the angle between them, and require it to be bounded away from 0 and from 180 degrees by a tolerance. Never
     use the planet-centred inertial velocities.
- Enumeration: build the arc list for each C (or each admissible set of flyby-speed pairs), build the directed
  graph, enumerate simple cycles up to a length cap, and take each cycle as a seed for the continuation in the
  moon mass (the project's mass ladder, which is the Bradley-Russell continuation applied to the shadowing
  orbit; the theory guarantees uniqueness and a distance of at least c epsilon from collision only for epsilon
  below an unquantified epsilon_0).
- Dense candidates: the theory chooses a dense but special set of rational frequencies for which no early
  collision occurs; in practice the project would generate arcs from a grid of flight times and keep those that
  pass the early-collision test.

Condition (3) and the "no early collision" requirement are the places where the theory is most easily
under-applied: the 2000 paper says additional collisions "really occur for many other collision trajectories as
the Jacobi constant is varied", so a search that skips test 3 will emit chains that contain unintended
encounters (READ, 2000 digest section 3(e)).

### (b) The direction-change condition against the demanded-turn gate (`src/cyclerfinder/verify/turn_gate.py`)

What the gate computes (READ in the code): `demanded_turn` is the angle between the incoming and outgoing
V-infinity vectors, `arccos(v_in . v_out / (|v_in| |v_out|))` in [0, pi], and `available_bend` is `2 asin(1 / (1 +
r_p v^2 / GM))` at `r_p = R + alt_floor` and `v = min(|v_in|, |v_out|)`; the flyby is feasible when the demanded
turn does not exceed the available bend; `required_alt_km` is the periapsis altitude at which an unpowered flyby
turns by exactly the demanded angle, `r_p = GM/v^2 (1/sin(demanded/2) - 1)` minus the radius.

How the two relate (INFERRED):
- The direction-change condition is the statement that the demanded turn is not 0 and not pi (measured in the
  frame relative to the body). The gate's demanded turn is the same quantity. So the theory's condition is a
  condition on the one number the gate already computes.
- The theory is a small-mass statement. The turn needed for a given periapsis scales as r_p = GM/v^2 (1/sin(delta/
  2) - 1), so the periapsis is proportional to the moon's mass (GM) at a given turn delta and relative speed v. The
  theorem's closest-approach bound, the distance at least c epsilon from collision (2000 Theorem 1.1, inequality
  (1.5)), is the same scaling: the closest approach is of order epsilon with a constant depending on the arcs. At
  vanishing mass any demanded turn delta strictly between 0 and pi is attainable at a vanishing periapsis, which is
  why the theory has no upper bound on the turn other than pi, and no finite radius; the 2000 paper's Remark 1
  (digest section 2.3) accommodates a finite body by requiring its radius to be less than c epsilon.
- The gate adds what the theory cannot: at the physical mass and the physical radius the periapsis needed must
  exceed R + alt_floor, i.e. the demanded turn must not exceed the available bend (ratio at most 1). A chain that
  satisfies the theory's direction-change condition at every encounter can still fail the gate when the moon is
  too light or the demanded turn too large (a demanded turn near 180 degrees needs a required altitude near minus
  the radius). Conversely, passing the gate does not give the theory's uniqueness or hyperbolicity, which need
  epsilon below an unquantified epsilon_0.
- Differences in the degenerate limits: the theory excludes both 0 and pi. The gate excludes pi in effect (the
  required altitude for a turn of pi is minus the body radius, infeasible) but it ACCEPTS a demanded turn of 0
  (ratio 0, `turn_feasible` true). A zero demanded turn means the periapsis is infinite, that is no encounter at
  all, so such an edge is not a flyby; the gate would pass it and the theory would not. The tolerance in the gate
  for "zero" (`turn_tol_rad`) is round-off only, so near-zero demanded turns are accepted. If the gate is used as
  an admissibility test in a search, the lower bound on the demanded turn (a minimum turn, or a maximum periapsis)
  must be added separately.
- The halving at the regularised saddle (section 1.4 above) suggests another way to state the same thing: in the
  moon-centred regularised variables the admissible angle changes are those below 90 degrees.

### (c) What hyperbolicity and the Lyapunov-exponent statements imply for navigation cost and for continuation

READ: the angle change in a close encounter is a smooth function of impact parameter with slope of order
1/epsilon, so Lyapunov exponents are of order log(1/epsilon) in units of the secondary's frequency (this paper
p. 2; 2006 Theorem 2.2, which the 2006 paper says "can be deduced from the proof ... but was not proved" in
2000). The same text says this gives "greater controllability".

INFERRED consequences for the project:
- Error growth. A position or velocity error at a flyby, expressed as a change of impact parameter, becomes an
  angle error 1/epsilon times larger after the flyby. For Titania and Oberon, epsilon (mass ratio to Uranus) is of
  order 4e-5, so log(1/epsilon) is about 10.1 and 10.2 (the adversarial review `docs/notes/2026-10-04-890-
  adversarial-review.md` section 8 reports exactly these numbers, and a leading multiplier of 8.381e5 per cycle,
  6.8 in logarithm per flyby, the same order). Navigation: an impact-parameter error is multiplied by roughly
  e^6.8, about 9e2, at each flyby, so a free-flying cycle needs a correction at or between flybys and only a
  few flybys can be flown before the error is of the order of the periapsis itself. A statement of the
  total cost would need the cost per unit of impact parameter, which the paper does not give.
- Controllability. The same slope is what makes a small manoeuvre effective: a change of the impact parameter by
  a small fraction of the periapsis rotates the outgoing direction by a large angle. This is why the flyby orbit
  is "controllable" in the sense of steering to a chosen next encounter with a small manoeuvre, which is the
  standard mission-design experience (this paper's Voyager remark). The effect is on the manoeuvre needed to
  steer, not the cost to hold a given orbit (where the instability works against).
- Conditioning of continuation. A single-shooting continuation over several flybys has a Jacobian whose condition
  number grows like the product of the per-flyby multipliers (here about 1e6 per cycle, so about 1e12 or more
  over two cycles), which is why long single shoots fail and why multiple shooting with a node at each
  encounter (as the project did) works. Parameterising the nodes by the b-plane impact parameter and
  V-infinity direction rather than by position and velocity keeps the encounter-to-encounter map well scaled,
  and a regularised integration (Levi-Civita in the plane) removes the step-size collapse near periapsis.
  Exploiting a symmetry of the orbit, as the project's symmetric multiple shooting does, halves the unknowns.
  These recommendations are mine, not the paper's.
- Existence versus detection. Uniform hyperbolicity of the invariant set (2006) means any such periodic orbit
  has real expanding and contracting multipliers; the project's #890 orbit has both multiplier pairs real (review
  section 5 as quoted in section 8), consistent with this (INFERRED, a consistency remark, not verification).

### (d) What would have to be proved or computed to extend the construction to two moons with different periods

All INFERRED; the printed papers do not address two secondaries, and none of this is established.

1. The setting. Moon A and moon B on concentric circular orbits with periods T_A and T_B about the primary. If
   T_A : T_B = p : q is rational the system is time-periodic with period lcm; otherwise it is quasi-periodic. The
   single-moon theorem needs a time-independent Lagrangian and an energy level. The nearest published
   generalisation, Bolotin (2005) on the elliptic restricted problem (time-periodic, one small body), is not held
   by the project (the 2000 and 2006 digests list it; the 2006 paper attributes the elliptic extension to
   Bolotin). That is the paper to read first.
2. The vertices. With no energy to fix, the zero-mass arcs are Kepler arcs between space-time events: arrival at
   a moon at epoch t_i, departure at the next moon at epoch t_{i+1}. For N flybys in the planar problem the
   unknowns are the N event epochs (plus discrete choices: revolution counts, branch), and the N equations are
   the equal-magnitude conditions of the relative velocities at each flyby. Counting N equations in N unknowns
   makes isolated solutions generic; this is the zero-mass closure the project already solves (the symmetric
   closures of `turn_gate_closures.symmetric_closure` have two unknowns, the flight time and the relative phase,
   and two magnitude equations).
3. Nondegeneracy. The determinant of the Jacobian of the magnitude-mismatch vector with respect to the epoch vector
   must be nonzero; this plays the role of the 2000 paper's nondegeneracy. It is a computation, and for the #890
   chain it can be done now.
4. Direction change at each event, in the frame relative to the moon at that epoch; computed as in (a) test 4.
5. Gluing. The theorem would need an implicit-function argument in the extended phase space (position, velocity and
   epoch) near the "anti-integrable limit" epsilon = 0, with the local encounter analysis done in the field of a
   MOVING moon (the transit time minus log epsilon plus a bounded term, and the closest approach of order
   epsilon) and uniform control of the time-dependent forcing between encounters. This is the part that is a
   proof, not a computation, and whether the estimates survive time dependence is not stated anywhere I have read.
6. Hyperbolicity. The Lyapunov-exponent statement for the two-moon case would follow only from the same gluing
   argument; the project's computed multipliers (both real, order 1e6 per cycle) are a numerical datum, not a
   proof.
7. Periodicity condition. A periodic orbit of the time-periodic problem needs the total duration of the chain
   to be a common multiple of the relevant periods (or the phase to repeat); for #890 the cycle time is 123.162
   days. Which commensurabilities admit chains is a computation over the graph.
8. Mass range. The theory is a small-mass statement with no numerical epsilon_0; the question whether the
   physical Titania and Oberon masses are inside the range of validity is answered only by computing (the
   project's continuation in mass, with the gate at the end), not by the theory.

### (e) The convention that caused the 2000 error, and which frame the project's gate measures

The error (READ, section 1.5): the 2000 paper proved the direction-change lemma for ellipses in the INERTIAL frame
(about the Sun; its Figure 3 and the Figure 1 here are drawn in that frame), but the theorem needs the condition
for the velocities in the ROTATING frame, the frame in which the secondary is at rest and in which the theorem is
formulated. Two arcs that leave a collision point in different, non-opposite directions in the inertial frame can
leave in parallel or antiparallel directions relative to the secondary: the transformation subtracts the
secondary's velocity, which is the same vector for both, but it is not a rotation and does not preserve angles
between the two arcs' velocities. In the inertial frame the two mirror ellipses always leave in different,
non-opposite directions, which is why the lemma looked true; relative to Jupiter they both leave radially when the
angular momentum is 1.

What the project's gate measures (READ in `src/cyclerfinder/verify/turn_gate.py` and
`src/cyclerfinder/verify/turn_gate_closures.py`):
- `turn_gate.evaluate_encounter` takes `Encounter.vinf_in` and `vinf_out`, documented as "body-relative velocities
  (km/s) in any common frame", and computes the angle between them with `bend_angle`. The docstring says the
  demanded turn is "frame-invariant, any common frame works (inertial, rotating, body-local), provided both
  vectors are expressed in it". That is correct for rotations of the frame, but it is true only because the
  inputs are already V-infinity (velocity relative to the body); the translation that makes the vectors
  body-relative is the part that matters, and the gate does not check that its inputs have been made body-relative.
- In `turn_gate_closures.symmetric_closure`, which builds the two-moon closures, the vectors are `out0 = leg0.v1 -
  w0`, `in1 = leg0.v2 - w1`, `out1 = leg1.v1 - w1`, `in2 = leg1.v2 - w2`, `out2 = leg2.v1 - w2`, where `leg.v1` and
  `leg.v2` are the Lambert arc's planet-centred inertial velocities at its ends and `w0`, `w1`, `w2` are the moon's
  inertial velocities at the same epoch and body. So each encounter compares the arriving and departing velocities
  relative to the SAME moon at the SAME epoch. For a moon on a circular orbit, the rotating-frame velocity at the
  moon's own position equals the inertial velocity minus omega cross r, and omega cross r at the moon's position
  is the moon's inertial velocity, so these V-infinity vectors are identical to the velocities in the frame rotating
  with the moon (INFERRED, simple kinematics for a circular orbit). The module header states "All frames are
  planet-centred (or Sun-centred) inertial, coplanar", which describes the frame in which the Lambert velocities
  are computed, and the subtraction of the body velocity is what supplies the body-relative frame.
- Conclusion: the project's turn gate measures the turn in the body-relative (rotating-with-the-body) frame, the
  frame the corrected direction-change condition requires, and does not repeat the 2000 paper's inertial-frame
  error for the encounters it evaluates.
- Cautions found in the code, none of which is a bug today: (i) the gate does not verify its inputs are
  body-relative (a caller passing planet-centred inertial velocities would reproduce the 2000 error silently);
  (ii) the wrap encounter in `symmetric_closure` compares `in2` and
  `out2` at the same epoch and body, which is correct in any frame; (iii) `turn_gate.encounters_from_vinf_nodes`
  with `wrap_rotation_rad` compares a V-infinity vector at one epoch with another rotated about +z by the home
  body's advance over one period, and `to_body_local` is documented for comparing vectors at different epochs
  "when the chain repeats in the body's rotating frame ... the case for circular-coplanar ideal models". Those
  rotations of body-relative vectors are valid only if the whole configuration repeats up to that rotation. With
  two moons of different periods the advance of the home moon over one cycle and the phase of the other moon
  differ, so for a two-moon chain the wrap comparison must be made at the same epoch (as `symmetric_closure` does)
  or the repeat condition must be checked first (INFERRED).
- A positive control for any future frame check (INFERRED, in the spirit of the project's rule about known-good
  cases): build a collision chain whose inertial-frame directions differ but whose body-relative directions are
  parallel (the exceptional case of section 1.5, G = 1, V-infinity purely radial) and confirm that the gate
  reports a demanded turn of 0 or 180 degrees for it, not a nonzero turn; and conversely a pair of ellipses with
  G not equal to 1 should report a turn strictly between 0 and 180 degrees.
