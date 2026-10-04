# Digest: Henon 2001, "Generating Families in the Restricted Three-Body Problem II. Quantitative Study of Bifurcations", part A (front matter, chapters 11 to 16: general equations and type 1)

M. Henon, *Generating Families in the Restricted Three-Body Problem. II. Quantitative Study of Bifurcations*, Lecture Notes in Physics Monographs m65, Springer, Berlin, 2001, ISBN 3-540-41733-8, DOI 10.1007/3-540-44712-1, xii + 308 pages in the file. Filed in the private paper corpus as
`henon-2001-generating-families-restricted-three-body-problem-II-quantitative-study-bifurcations-lnp-m65-springer-doi-10.1007-3-540-44712-1.pdf` (308 PDF pages; text layer present but its mathematics is garbled in places).
This note is part A. Part B (chapters 17 to 23, type 2 and the conclusions) is `docs/notes/2026-10-04-digest-henon-2001-generating-families-II-part-b-type-2.md`; I write nothing about those chapters beyond cross-references.
Volume I is `docs/notes/2026-10-04-digest-henon-1997-generating-families.md`. Related digests read for this note: `docs/notes/2026-10-04-digest-hitzl-henon-1977-critical-generating-orbits.md`, `docs/notes/2026-10-04-digest-henon-1968-consecutive-collision-orbits.md`, `docs/notes/2026-10-04-digest-bruno-1981-periodic-flybys-of-the-moon.md`, `docs/notes/2026-10-04-digest-brjuno-1978-periodic-solutions-arcs-mu-0.md`, `docs/notes/2026-10-04-digest-guillaume-1973-periodic-symmetric-solutions-small-mu.md`, `docs/notes/2026-10-04-digest-guillaume-1975-extension-breakwell-perko-matching.md`, `docs/notes/2026-10-04-digest-guillaume-1975a-linear-analysis-second-species.md`, `docs/notes/2026-10-04-digest-perko-1976-second-species-O-mu-near-moon.md`.
Evidence tags: READ (printed page) means read in the text layer and, for every formula and table I rely on, on the page image; COMPUTED means my own arithmetic on 2026-10-04 (scratch scripts, not kept); DERIVED means my own algebra from printed relations (not printed in the book); INFERRED means my reading across sources.

## 0. Page numbers and what was read

Chapters 11 to 23 continue the numbering of volume I, and the page numbering restarts at 1 with chapter 11. For the whole of chapters 11 to 16 and the appendix 16.8 the offset is constant: PDF page = printed page + 13 (PDF 14 is printed 1, PDF 52 is printed 39, PDF 159 is printed 146, checked at both ends and in the middle). The preface is PDF 6 (printed V), contents PDF 8 to 13. Every page of PDF 6 and 14 to 159 was read in the text layer. Page images were read for: (12.31) (printed 21), (12.78) (printed 30), (12.106) to (12.110) and Fig. 12.3 (printed 34 to 35), (12.111) to (12.115) (printed 35 to 36), (13.41) to (13.50) (printed 47 to 50), (14.12) to (14.16) (printed 82), (14.18) to (14.21) (printed 86), Fig. 14.5 (printed 87), Tables 13.1 to 13.2 (printed 73), 13.7 to 13.10 (printed 77 to 78), and Table 15.2 (printed 112). Tables 13.3 to 13.6, 14.1 to 14.5 and 15.3 to 15.7 were read in the text layer only, where the column layout is scrambled; their content is described but not transcribed (section 8.3).

## 1. What the book is and how it is organised

Volume I (digest above) defined generating orbits (limits of periodic orbits as mu tends to 0), the arc families (T, S, and the first species ellipses), the three bifurcation types and the invariants (symmetry, side of passage) that settle many junctions. Where several families meet at a bifurcation orbit the invariants fail. Volume II (READ, preface p.V) replaces the qualitative argument by an asymptotic study of the families at small mu near a bifurcation orbit: "in almost all cases, the first-order asymptotic approximation of the families in the neighbourhood of the bifurcation can be derived", which "allows, in particular, a quantitative comparison with numerically found families". Chapter 11 gives the machinery; chapters 12 to 16 treat type 1 (the generating orbit is a Keplerian ellipse of rational period that crosses the Moon's circle transversally at two distinct points); chapters 17 to 23 treat type 2 (tangent to the circle; part B); type 3 (the retrograde unit circle) "had not yet been completed at the time of writing" (preface). The author states the work is "sometimes lacking in mathematical rigor" and rests his confidence on agreement with volume I, with numerical computation, internal consistency and intuition (preface).

The type 1 programme has three layers, all of which are in this part:
1. Chapter 12: the general derivation, valid for every n, in regimes of the distance nu of the Jacobi constant from the bifurcation value, ending in an explicit finite system (12.114) for nu = 1/2.
2. Chapters 13 and 14: solving (12.114) for partial (13) and total (14) bifurcations; the output is the junction tables (Tables 13.1 to 13.10, 14.1 to 14.5).
3. Chapters 15 and 16: two independent re-derivations (Newton polyhedra after Bruno 1998/2000; and an algebraic proof of the general structure for all n).

## 2. Chapter 11: notation and general equations (pp.1 to 16)

### 2.1 The O notation (11.2, pp.1 to 4; READ)

Following Graham, Knuth and Patashnik 1989 sect. 9.2. O[g(x1..xn)] is the set of functions f with |f| <= C|g| for |x_j| <= eps_j (Definition 11.2.1); g is generally a product of powers x1^q1 ... xn^qn with real, possibly negative exponents; C and eps_j may depend on other finite variables y. An expression with O terms is the set obtained by letting each O range over its domain (Definition 11.2.2). The "=" sign in an equation containing O means "is a subset of" (Definition 11.2.3), so the relation is not symmetric. O cannot be nested (g must not itself contain O). Rules used throughout (printed (11.3) to (11.14)):
- g = O(g); h O(g) = O(g) for a function h of the other variables, and |h| g = O(g).
- Substitution: only a left-hand O-expression may be substituted into a right-hand one.
- O(g1 + g2) = O(g1) + O(g2) (reverse false); O(g1) O(g2) = O(g1 g2) with the reverse true; larger exponents are absorbed: O(x^q) is contained in O(x^q') for q >= q'.
- Elimination of a composite term: O(g1) + O(g1^lambda g2^(1-lambda)) + O(g2) = O(g1) + O(g2) for 0 <= lambda <= 1.
- Truncation of a convergent series: f = sum_{j<=q} h_j x^j + O(x^(q+1)).
- Shorthand: O(g1) + O(g2) + ... is written O(g1, g2, ...).
Logarithmic factors are deliberately neglected in orders of magnitude (p.8, citing Perko 1976b p.399); an estimate such as Delta p = O(mu ln eps) is written O(mu).

### 2.2 Set-up and numbering (11.3.1, pp.4 to 6)

- A second species bifurcation orbit Q is a closed chain of basic arcs joined at collisions with the Moon M2. Partial bifurcation of order n: collisions numbered 0 to n (both ends included), basic arcs 1 to n; internal collisions C = {1, ..., n-1}. Total bifurcation of order n: collisions 0 to n-1 with collision 0 at the origin, arcs 1 to n, arc n joining collision n-1 to collision 0, C = {0, ..., n-1} (arcs are numbered from the origin; the book writes the index of C as "i taken modulo n" in chapter 14). The arc set A = {1, ..., n} in both cases.
- Coordinates: the fixed-direction system (X, Y) with origin in M2, related to the synodic (x, y) by (11.17): X = (x - 1) cos t - y sin t, Y = (x - 1) sin t + y cos t. Equations of motion (11.18), (11.19), and in vector form (11.31): p'' = -(1 - mu) (p - p_M)/|p - p_M|^3 ... with p_M the (known) position of the Sun-side primary seen from M2, (11.32) p_M = (-cos t, -sin t) up to the book's sign convention; the mu term is -mu p/|p|^3. [The exact printed form of (11.31) is partly illegible in the text layer; the structure is standard and is not used below.]
- For mu small but non-zero the collisions become encounters, "following the tradition of stellar dynamics", with minimum approach distance eps_i at encounter i. The relative velocity modulus at a collision is the same, v, in the rotating and the fixed frames (11.20, 11.21), finite for types 1, 2 and 3 (11.22), and the same at every collision of the orbit.
- Deflection by encounter i (Mihalas and Routly 1968 p.174), (11.23): theta_i = 2 mu / (v^2 eps_i) (the small-angle form; READ p.6).

### 2.3 The exponent nu_i and its constraint (11.3.1, p.6)

Fundamental assumption (11.27): eps_i = O(mu^(nu_i)) as mu tends to 0, so theta_i = O(mu^(1 - nu_i)) (11.28). Because eps_i must tend to 0 (11.25) and theta_i must tend to 0 near a bifurcation orbit (11.26), 0 < nu_i < 1 (11.29), the same constraint as Guillaume 1971 p.98. This nu_i is the exponent of the encounter distance, and it is NOT the nu of section 2.6 (below). The two are related in section 3.5.

### 2.4 Intermediate arcs and their accuracy (11.3.2 to 11.3.4, pp.6 to 11)

Between encounters the real orbit is replaced by the osculating Keplerian ellipse (about the Sun, mass 1) taken at a point far from M2 (|p(t_0i)| = Theta(1)): the intermediate orbit i, p(t) = p_i(t) + Delta p_i(t) (11.30). Integrating the perturbation equation (11.36), Delta p_i'' = O(mu) p/|p|^3 + O(Delta p_i) ... and estimating the dominant term near the encounter (the logarithmic factor neglected), the printed results are:
- Delta p_i'(t) = O(mu / (t_i - t)) and Delta p_i(t) = O(mu ln eps) = O(mu) (11.38 to 11.41).
- At t = t_i - eps_i: Delta p_i' = O(mu / eps_i), Delta p_i = O(mu) (11.42); inside the encounter Delta p_i' = O(mu / eps_i) and Delta p_i = O(mu) (11.44 to 11.46).
- **Proposition 11.3.1:** the intermediate orbit i approximates the true orbit to O(mu) over the whole piece between encounters i-1 and i, and also during those encounters.
So the true periodic orbit is, to O(mu), a sequence of Keplerian arcs joined near (not at) the Moon, within distance Theta(eps_i) of M2 with nu_i in (0, 1); the joining point has latitude because the two successive arcs differ by an angle of order theta = mu^(1 - nu_i) over a span mu^(nu_i), a displacement O(mu). The only case in the whole book where O(mu) is not enough is the bifurcation 2P1 (p.7; section 23.2 in part B).

The velocity error is refined to second order (11.3.4): approximating the encounter by uniform straight-line motion on the osculating tangent (impact parameter d_b = O(eps_i), unit tangent vector, polar angle phi with phi_oi = -pi/2), the printed results are (11.61) and (11.62):
Delta p_i'(t_i) = (mu / (v_b d_b)) [ i_b cos phi_b - j_b (1 + sin phi_b) ] + O(mu) + O(mu^2/eps_i^2)
(and the mirror expression with (1 - sin phi_a) for the arc after the encounter). The subscripts b and a are before and after.

### 2.5 Matching relations at an internal encounter (11.3.5, 11.3.6; pp.11 to 13)

For i in C (READ):
- First matching relation (11.65): p_{i+1}(t_i) = p_i(t_i) + O(mu).
- Second matching relation (11.70): p'_{i+1}(t_i) - p'_i(t_i) = -(2 mu / (v_b d_b)) j_b + O(mu) + O(mu^2/d_b^2), with j_b the unit vector orthogonal to the incoming tangent, directed from M2 toward the tangent, d_b the impact parameter and v_b the speed (the first term is the hyperbolic-flyby kick; it is the Guillaume 1971 relation; it is the same in a rotating frame centred on M2).
- The second term of (11.70) is significant when mu/d_b > mu and mu/d_b > mu^2/d_b^2, which gives exactly the two conditions (11.25) and (11.26); "a confirmation that (11.70) is indeed the correct relation for bifurcations".
- At the ends of a partial bifurcation (11.71): the bifurcating arc starts and ends in true collisions, p_1(t_0) = 0 and p_n(t_n) = 0 to O(mu).
- For mu = 0 (11.3.6): eps_i is zero at a node and non-zero at an antinode, theta_i the reverse, so eps_i theta_i = 0 (11.72) replaces (11.23); the matching equations still hold, with (11.70) multiplied by d_b: (11.79) d_b [p'_{i+1}(t_i) - p'_i(t_i)] = -(2 mu / v_b) j_b + O(mu d_b) + O((mu^2/d_b) ...), which is meaningful at a node (d_b = 0 and mu = 0 give 0 = 0) and is the form actually used in chapter 12.

### 2.6 The general method and the parameter nu (11.4, pp.14 to 16)

Orbits near a bifurcation depend on two small parameters, mu (distance from the two-body problem) and Delta C (distance from the bifurcation orbit; the curves C = const are approximately parallel lines near Q). The strategy: start from the asymptotic branches of the characteristics (the parts far from Q approaching the mu = 0 branches, volume I Fig. 1.1) and work inward until the branches are joined two by two. Because the interest is in mu tending to 0, define the exponent (11.80):

**Delta C = O(mu^nu)**, equivalently ln|Delta C| = nu O(ln mu) (11.81).

nu = 0 is the mu = 0, Delta C not zero situation (the full lines of volume I Fig. 1.1); increasing nu goes closer to the bifurcation; large nu is the vicinity of Q. Each regime of nu (an open interval or an isolated value) is treated by a nine step procedure (printed p.15 to 16):
1. extrapolate the orders of magnitude from the previous regime;
2. choose new variables of order 1 (intuition required);
3. rescale to simplify;
4. collect dominant terms on the left, divide so they are O(1), the right-hand side must be o(1);
5. find the condition under which the right-hand side is o(1) (this fixes an upper limit of nu), and set the right-hand side to zero (the asymptotic equations);
6. solve them; with several solutions select by continuity from the previous regime;
7. compute the Jacobian determinant |J| (here "Jacobian" means the determinant) of the asymptotic system at the asymptotic solution; if non-zero, the implicit function theorem gives a true solution near it for mu > 0, at a distance given by the largest right-hand member;
8. optionally refine individual error estimates;
9. return to the physical variables.
**Regimes for type 1 (11.82): nu = 0; 0 < nu < 1/2; nu = 1/2; nu > 1/2.** For type 2 (11.83, part B): nu = 0, 0 < nu < 1/3, nu = 1/3, 1/3 < nu < 1/2, nu = 1/2 (and nu > 1/2 does not exist; part B).

## 3. Chapter 12: quantitative study of type 1, general n (pp.17 to 38)

### 3.1 Fundamental equations (12.1, pp.17 to 23)

Definitions: t_i is chosen as the time of intersection of the true orbit with the orbit of M2; y_i is the oriented arc length along M2's orbit from M2 to M3 at that time (the lead of M3 over M2; Fig. 12.1). Because type 1 has non-zero radial velocity at the crossing, the intermediate orbit's crossing time differs from t_i by O(mu) (12.1 to 12.4). The intermediate orbit is a point (A_i, Z_i) of the (A, Z) plane (volume I chapter 4; A = a^(3/2), Z from the orbital elements), A_i = A + Delta A_i, Z_i = Z + Delta Z_i, C_i = C + Delta c_i (12.5). The bifurcation orbit alternates first basic arcs PQ and second basic arcs QP (volume I sect. 6.2.1.2); s_i = -1 if the arc is a first basic arc (the encounter at P), +1 for a second (at Q) (12.12).

Arc relation (12.17), merged from (12.8) and (12.11) (READ pp.18 to 19):

y_i - y_{i-1} = -2 pi s_i (dZ/dC)_A Delta C + 3 pi s_i sqrt(a) [bracket in (dZ/dA)_C, beta_0, J and s_i] Delta a_i + O(mu) + O(Delta C^2) + O(Delta a_i^2), for i in A

(the bracket is illegible in the text layer; its content is the constants G3 and K below, which is the form actually used in (12.32a)). Delta A_i = (3/2) sqrt(a) Delta a_i + O(Delta a_i^2) (12.16). Additional relation for two consecutive arcs (12.20, READ p.20), with no O(Delta C^2) term:
y_{i+1} - y_{i-1} = -G3 (1 + K s_i) Delta a_i - G3 (1 - K s_i) Delta a_{i+1} + O(mu) + O(Delta a_i^2) + O(Delta a_{i+1}^2), i in C (12.33).

Encounter relation (12.29, READ p.21; "identical with the equation obtained by Guillaume 1971, p.83"):
**y_i (Delta a_{i+1} - Delta a_i) = -4 mu a^2 / v + O(mu y_i) + O(mu Delta a_i) + O(mu Delta a_{i+1}), i in C.**
For a partial bifurcation y_0 = O(mu), y_n = O(mu) (12.30).

Constants (12.31, READ on the page image):
- G1 = 4 a^2 / v (always positive);
- G2 = 2 pi (dZ/dC)_A;
- G3 = (3 pi J sqrt(a)) / 2 (always positive; J is the integer of the supporting ellipse, volume I);
- K = 1 + (2/J) [ (dZ/dA)_C - beta_0 ] (reproduces volume I (8.14); beta_0 is the arc's slope parameter of volume I eq. 6.17).
G2 and K are functions of the partial derivatives of Z(A, C) (volume I (4.28), Fig. 4.15), "expected to be arbitrary, ordinary real numbers"; hence two conjectures: **Conjecture 12.1.1: K never takes an integer value** (supported by numerical computation of K for many bifurcations, volume I sect. 8.2.1) and **Conjecture 12.1.2: G2 never vanishes.**

The fundamental equations for type 1 (12.32), with s_i = (-1)^i (12.35; the starting point is taken in P, which costs nothing because changing the sign of K and of all s_i maps solutions to solutions, so only the start in P is studied but all real K are covered):
- (a) y_i - y_{i-1} = -G2 s_i Delta C - G3 (1 + K s_i) Delta a_i + O(mu) + O(Delta C^2) + O(Delta a_i^2), i in A;
- (b) y_i (Delta a_{i+1} - Delta a_i) = -G1 mu + O(mu y_i) + O(mu Delta a_i) + O(mu Delta a_{i+1}), i in C;
- (c, d) y_0 = O(mu), y_n = O(mu) for a partial bifurcation.
Count: a partial bifurcation has 2n + 1 equations for 2n + 2 variables, a one-parameter (the usual family); a total one has 2n equations for 2n + 1 variables.

### 3.2 Exclusion of successive identical T-arcs (12.2, pp.23 to 24)

The formalism applies also to an ordinary generating orbit (one family only). For a sequence of identical T-arcs (collisions only in P, or only in Q; total or partial T-sequence) the arc relation becomes y_i - y_{i-1} = -3 pi J sqrt(a) Delta a_i + ... (12.38). Multiplying by Delta a_i and summing with the encounter relation gives 0 = -(3 pi J sqrt(a)) sum (Delta a_i)^2 - (4 mu a^2/v) h' + smaller terms (12.40, 12.42), impossible for a total T-sequence and for a partial sequence of h >= 2 arcs. This **proves volume I Proposition 4.3.2** (an ordinary generating orbit of the second species cannot contain two identical T-arcs of type 1 in succession; a single T-arc is allowed). The book notes this proposition had "already been used many times".

### 3.3 nu = 0 (12.3, pp.25 to 27)

The equations (12.43), (12.44) reduce to the known mu = 0 solutions (Delta C finite):
- First species ellipse (two basic arcs): Delta a_1 = Delta a_2 = 0 (exact invariance of the period), and, using the symmetry y_1 = -y_0, **y_1 = (G2/2) Delta C + O(Delta C^2)** (12.50).
- T-arc: Delta a = 0 on both basic arcs, and **y_{i+1} = G2 s_i Delta C + O(Delta C^2)** at the antinode (12.57).
- S-arc of m basic arcs from node i to node i + m: all Delta a_{i+j} equal to a common value Delta a, **Delta a = G2 s_i Delta C / ((m s_i - K) G3)** (12.61, from 12.60), and y_{i+j} = (m-j)... : **y_{i+j} = [(m - j)/(m s_i - K)] G2 Delta C for j odd and -(j/(m s_i - K)) G2 Delta C for j even** (12.62; 12.110 below).
Nodes have y = 0 at this order (collisions), antinodes y of order Delta C.

### 3.4 0 < nu < 1/2 (12.4, pp.27 to 34)

Variables y_i = Delta C y_i*, Delta a_i = Delta C x_i* (12.63) and a second rescaling (12.64); the right-hand members are O(mu Delta C^-1) + O(Delta C) and O(mu Delta C^-2), which are o(1) exactly when **0 < nu < 1/2** (12.67). The asymptotic equations (12.68), (12.69) are linear, and solved for each species:
- First species orbit: Jacobian |J| = 2 (12.74), so for 0 < nu < 1/2 a solution exists with error O(Delta C) + O(mu Delta C^-2), and
  **y_0 = -(G2/2) Delta C [1 + ...], y_1 = (G2/2) Delta C [1 + ...], Delta a_1 = (K + 1) (G1/G2) mu Delta C^-1 [1 + ...], Delta a_2 = (K - 1)(G1/G2) mu Delta C^-1 [1 + ...]** (12.78, READ on image; the bracket is 1 + O(Delta C) + O(mu Delta C^-2)).
- Second species orbit or bifurcating arc: each T- or S-arc decouples. A T-arc: node y's of order mu/Delta C, antinode Y* = 1 (12.81). An S-arc of m basic arcs: X = 1/(m - K s_i) with alternating signs, Y* = (m - j)/(m - K s_i) (odd j) and j/(m - K s_i) (even j) (12.86, 12.87). The Jacobian is a product of per-arc determinants: -2 Y*_{i+1} for a T-arc (12.89) and a product of Y*'s for an S-arc (12.90), non-zero.

The refined results for the physical variables (12.106) to (12.110), READ on the page image (printed p.34), with the bracket [1 + O(Delta C) + O(mu Delta C^-2)] on every line, are:
- TS node (T-arc then S-arc, the S-arc having m_a basic arcs): **y_i = -(m_a s_i - K) (G1 G3 / G2) mu Delta C^-1**.
- ST node (S-arc of m_b arcs then T-arc): **y_i = -(m_b s_i + K) (G1 G3 / G2) mu Delta C^-1**.
- SS node: **y_i = -[(m_b s_i + K)(m_a s_i - K)] / [(m_b + m_a) s_i] (G1 G3 / G2) mu Delta C^-1**.
- T-arc antinodes: Delta a_{i+1} = [(K + s_i) - g'(m_b s_i + K) + g''(m_a s_i - K)] (G1/(2 G2)) mu Delta C^-1, Delta a_{i+2} = [(K - s_i) - g'(m_b s_i + K) + g''(m_a s_i - K)] (G1/(2 G2)) mu Delta C^-1, y_{i+1} = G2 s_i Delta C, where g' = 0 or 1 according as node i is an end of the bifurcating arc or a junction with an S-arc, and g'' likewise for node i + 2 (12.103 to 12.105, 12.109).
- S-arc antinodes: Delta a_{i+j} = [1/(m s_i - K)] (G2/G3) Delta C, y_{i+j} as in (12.62) (12.110).
- Fig. 12.3 (READ on the image) plots these orders of magnitude on log-log axes between mu^(1/2) and 1: node y goes from mu^(1/2) at Delta C = mu^(1/2) down to mu at Delta C = 1 (slope -1 in Delta C); antinode y goes from mu^(1/2) up to 1 (slope +1); Delta a at a T-arc or first species orbit goes from mu^(1/2) down to mu; Delta a at an S-arc is of order Delta C. All of them equal mu^(1/2) at Delta C = mu^(1/2).

Side of passage (12.4.3): y_i has the sign of s_i sigma cos(phi), which, with the definition of the sign of G2 Delta C, recovers the side-of-passage rules of volume I sect. 8.2.1 and 8.3.1 (an independent confirmation of volume I's chapter 8).

### 3.5 nu = 1/2 (12.5, pp.35 to 36): the exact transition

At Delta C = O(mu^(1/2)) all y_i and Delta a_i are of order mu^(1/2); node and antinode, T- and S-arcs, and the first and second species fuse. Variables (12.111), (12.112): y_i = mu^(1/2) Y_i*, Delta a_i = mu^(1/2) X_i*, Delta C = mu^(1/2) W*, then Y_i*, X_i*, W* are rescaled by sqrt(G1 G3) and G2; the printed result (12.115, READ on the image):
**y_i = s_i sqrt(G1 G3) mu^(1/2) Y_i + O(mu), Delta a_i = -s_i sqrt(G1/G3) mu^(1/2) X_i + O(mu), Delta C = -(1/G2) sqrt(G1 G3) mu^(1/2) W.**
(The printed (12.115) omits the factor sign(cos phi) that (12.112) and (12.117) include; the magnitudes are unaffected. Printed slip, flagged.)
The asymptotic system (12.114) is
- Y_i + Y_{i-1} - W - (1 + K s_i) X_i = 0, i in A;
- Y_i (X_{i+1} + X_i) + 1 = 0, i in C;
- Y_0 = 0, Y_n = 0 (partial bifurcation).
It "cannot be solved explicitly in general" (step 6; chapters 13, 14 study it). Its Jacobian "is non-zero in general (but it can vanish in isolated points on the characteristics)"; the distance of the true solution from the asymptotic one is O(mu^(1/2)) in these variables.

### 3.6 nu > 1/2 (12.6, pp.36 to 37)

The nu = 1/2 system is solved in chapters 13 and 14 for all branches, which are joined two by two, so the programme "is completed". For completeness the regime Delta C << mu^(1/2) corresponds to W tending to 0; the characteristics in the (W, Y_1) plane sometimes cross the axis W = 0 (Figs 13.1 to 13.7, 14.1 to 14.4). The asymptotic equations (12.119) are those of (12.114) without the W-dependence in the error term; solutions are a subset of those of (12.114) with W = 0. They are isolated, and the book says nothing further: no new branches exist for nu > 1/2.

### 3.7 Summary of the type 1 scaling (what chapter 12 delivers)

For type 1 the following hold, all with the factors G1 = 4a^2/v, G3 = 3 pi J sqrt(a)/2, G2 = 2 pi (dZ/dC)_A and K = 1 + (2/J)((dZ/dA)_C - beta_0), which are constants of the bifurcation orbit:

| Delta C regime | node lead y | antinode lead y | Delta a |
|---|---|---|---|
| nu = 0 (Delta C finite) | Theta(mu / Delta C) -> 0 (collision at mu = 0) | G2 s Delta C (T-arc) or order Delta C (S-arc) | Theta(mu / Delta C) at a T-arc, Theta(Delta C) in an S-arc |
| 0 < nu < 1/2 | -(...)(G1 G3/G2) mu / Delta C, eqs 12.106 to 12.108 | order Delta C | as above |
| nu = 1/2 | sqrt(G1 G3 mu) |Y| = Theta(mu^(1/2)) | Theta(mu^(1/2)) | Theta(mu^(1/2)) |
| nu > 1/2 | Theta(mu^(1/2)) | Theta(mu^(1/2)) | Theta(mu^(1/2)) |

Connecting to the encounter distance exponent of section 2.3 (INFERRED; the book does not write it): since the lead y_i is measured along M2's orbit and the impact parameter is d = |y_i| |sin phi_i| (volume II (12.28) p.21, with phi_i the angle of the relative velocity), eps_i = Theta(|y_i|) when the radial velocity is non-zero (as it is for type 1). Hence the exponent nu_i of (11.27) is nu_i = 1 - nu at a node, nu_i = nu at an antinode, for 0 < nu < 1/2, and nu_i = 1/2 at both when nu >= 1/2. This is the dictionary used in section 9.

## 4. Chapter 13: partial bifurcation of type 1 (pp.39 to 78)

### 4.1 The system and its symmetries (13.1, pp.39 to 40)

(12.114) with s_i = (-1)^i: Y_i + Y_{i-1} - W - (1 + K s_i) X_i = 0 (i = 1..n), Y_i (X_{i+1} + X_i) + 1 = 0 (i = 1..n-1), Y_0 = Y_n = 0 (13.1). 2n + 1 equations in the 2n + 2 variables W, Y_0..Y_n, X_1..X_n: one-parameter families, the ordinary families of periodic orbits for a given mu (the map from (W, Y, X) to (Delta C, y, Delta a) is one to one at fixed mu). Properties:
1. Y_i (i = 1..n-1) never vanishes (13.1b), so each Y_i keeps its sign along a family: **Broucke's principle of invariance of the side of passage**, recovered.
2. Symmetry Sigma (13.3): (Y_i, X_i, s_i) -> (Y_{n-i}, X_{n+1-i}, s_{n+1-i}) maps solutions to solutions (the fundamental symmetry of the restricted problem).
3. Symmetry Sigma' (13.4): (Y_i, X_i, W) -> (-Y_i, -X_i, -W). Hence one studies Y_1 > 0 (13.5); the figures show only the half plane Y_1 > 0.
4. Eliminating X_i gives the three-term relation (13.6):
   **(1 + K s_i) Y_{i+1} + (1 - K s_i) Y_{i-1} + 2 Y_i + (1 - K^2)/Y_i - 2W = 0, i = 1..n-1.**
5. With xi_k = Y_{2k+1}, eta_k = Y_{2k+2} (13.7), (13.6) is a plane mapping (13.8), area preserving; Henon's unpublished numerical explorations show the mixture of regular and chaotic orbits of a non-integrable system, so (13.1) is conjectured not to be solvable explicitly for general n.
6. The sign s of the branch (volume I definition (8.22)) is the sign of W (13.9).

### 4.2 Asymptotic branches (13.1.1, 13.1.2, 13.1.5; pp.40 to 47)

For |W| large (equivalent to nu < 1/2, chapter 12.4) with X_i = W X_i-bar etc. (13.10), the printed asymptotics (13.11) to (13.14) are: first species orbit, X_1, X_2 = O(1/W), Y_0 = -W/2... (as (12.78)); node Y_i = O(1/W); T-arc Y_{i+1} = W + O(1/W); S-arc X_{i+j} = (-1)^j W/(m - K s_i), Y_{i+j} = ((m - j)/(m - K s_i)) W (j odd) and ... (j even) (the 12.109, 12.110 forms with Delta C replaced by W). Always X_i = O(W) (13.15) and Y_i at an antinode is O(W) (13.17).
Variational equations (13.18) to (13.30): dY_i = -dY_{i-1} + (1 + K s_i) dX_i, dX_{i+1} = -dX_i + dY_i / Y_i^2 (the Y_i^2 denominator comes from differentiating the encounter equation Y_i (X_{i+1} + X_i) + 1 = 0, DERIVED). Across a T-arc the variations are dX_{i+2} = O(W^-2) dY_i - dX_{i+1}[1 + O(W^-2)], dY_{i+2} = dY_i[1 + O(W^-2)] - 2 dX_{i+1}[1 + O(W^-2)] (13.21, 13.22); across a S-arc dY_{i+m} = -(-1)^m dY_i[...] + (m - K s_i) dX_{i+1}[...] (13.24). Define u_alpha = -2 for a T-arc and u_alpha = 2... (printed (13.26), the garbled symbol is the quantity m_alpha + K s for an S-arc; the sign rule (13.27) is u_alpha > 0 for a normal S-arc, < 0 for an abnormal S-arc or a T-arc). Iterating from the origin gives (13.29), and the orders of magnitude (13.30): "The variations are strongly amplified after each node". The Jacobian for |W| tending to infinity is dY_n/dX_1 = u_1 u_2 ... u_N / (Y_{i_1}^2 ... Y_{i_{N-1}}^2) = Theta(W^(2N-2)) (13.40): it never vanishes in that limit.

### 4.3 Jacobian and stability (13.1.3, 13.1.4; pp.43 to 46)

The 2n x 2n Jacobian of (13.31) reduces to the n-variable chain (13.33). **Proposition 13.1.1: in a partial bifurcation of type 1 the Jacobian vanishes if and only if the bifurcating arc is critical**, "perturbations in the direction of the departure velocity have no effect on the impact parameter at the next encounter" (Hitzl and Henon 1977b p.1029), written dY_n/dX_1 = |J| = 0 (13.34 to 13.37). A stable interval appears when the stability index jumps from +infinity to -infinity (Henon and Guyot 1970 p.364). A critical bifurcating arc is either an extremum of the characteristic in W (the Jacobian can be made non-zero by reparametrising with X_1; not a true singularity), or a saddle of the surface Y_n(W, X_1), i.e. an intersection of two families at that point (a true singularity: for mu > 0 the four branches can be joined in different ways, as in a generating-orbit bifurcation of volume I Fig. 1.1). The intersections of this paragraph are at the "finer level of description" of the bifurcation (they shrink to the bifurcation point in (Delta C, y_i)), and should not be confused with bifurcations proper.
First derivatives (13.39): dY_1/dX_1 = 1 - K, dY_2/dX_1 = (1 - K^2)... the printed expressions are partly garbled in the text layer and are not used.

### 4.4 Small n: closed forms (13.2, pp.47 to 53)

- **n = 2, partial (1P2)** (13.41): **2 Y_1 + (1 - K^2)/Y_1 - 2 W = 0** (READ on the image), a hyperbola with asymptotes Y_1 = 0 and Y_1 = W; one of three shapes for K < -1, |K| < 1, K > +1 (Fig. 13.1). Branch +2 is the T-arc (asymptotic to Y_1 = W), branch +-11 the two basic S-arcs (asymptotic to Y_1 = 0, Y_1 = (1 - K^2)/(2W) + ...). For |K| < 1 the characteristic has an extremum in W, the critical arc, **Y_1 = sqrt((1 - K^2)/2), W = sqrt(2 (1 - K^2))** (13.42); the Jacobian is generally non-zero. This W is the smallest |W| reached on the family, so Delta C cannot approach zero closer than |Delta C|_min = sqrt(2(1 - K^2) G1 G3 mu)/G2 (DERIVED from 13.42 and 12.115).
- **n = 3 (1P3)** (13.43): the two equations subtract to (Y_1 - Y_2)(Y_1 Y_2 + K - 1) = 0 (13.44), so two families: Y_1 = Y_2 (symmetric bifurcating arc), 2W = (3 - K) Y_1 + (1 - K^2)/Y_1 (13.45); and Y_1 Y_2 = 1 - K (asymmetric), W = Y_1 + (1 - K)/Y_1 (13.46). Hyperbolas with asymptotes Y_1 = 0, Y_1 = W, Y_1 = 2W/(3 - K). For K < 1 they intersect at **Y_1 = sqrt(1 - K), W = 2 sqrt(1 - K)** (13.47), an extremum in W of the second family. Critical arcs: first family, (3 - K) Y_1^4 - 4(1 - K) Y_1^2 + (1 - K^2)(1 - K) = 0 (13.48), solutions (13.49) (the intersection, K < 1) and **Y_1 = sqrt((1 - K^2)/(3 - K)), W = sign(3 - K) sqrt((1 - K^2)(3 - K))** (13.50; exists for -1 < K < 1 or K > 3); second family Y_1^4 - 2(1 - K) Y_1^2 + (1 - K)^2 = 0 (13.51), double root (13.49). The junctions agree with volume I sect. 7.3.1.1 and 8.4.1 when the complement is symmetric (Restriction 7.3.1); for K < 1 and an unsymmetric complement the quantitative approach cannot decide, "likely a real fact": the junction then depends on higher-order perturbations at the two ends.
- **n = 4**: the three coupled equations (13.52) are not solvable in closed form; computed numerically for K = -4, -2, 0, 2, 4 (Figs 13.3 to 13.7, with the n = 3 and n = 2 characteristics dashed and dotted). The figure is qualitatively the same inside each of K < -3, -3 < K < -1, -1 < K < +1, +1 < K < +3, +3 < K (proof in 13.3). The junctions for n = 4 are fully determined for all K, where volume I's invariants left some open.

### 4.5 The positional method (13.3, pp.51 to 73)

Observation: characteristics of different n never intersect in the (W, Y_1) plane, since Y_n = 0 for a branch of order n < n' while Y_n not zero on a branch of order n' (13.1b), and Y_n is a single-valued function of (W, Y_1). So the characteristics of orders up to n - 1 divide the (W, Y_1) plane into regions, and
- **Proposition 13.3.1:** two branches of order n can be joined only if they lie in the same region.
- **Proposition 13.3.2:** the same must hold in the (W, Y_{n-1}) plane, which by the symmetry Sigma reduces to the (W, Y_1) plane with the same K for odd n and the opposite K for even n; two branches can be joined only if their symmetric branches lie in the same region.
The position of a branch in the plane is found from the large-|W| ordering. Two branches with common first arcs U_1..U_alpha (arcs of the bifurcating arc as T or S with m basic arcs) have dY_{i_alpha}/dX_1 = Theta(W^(2(alpha-1))), with sign determined by the number of T-arcs and abnormal S-arcs (13.53 to 13.57); the ordering at the end of the common part is read from (13.62) (T then S: Y_{i_alpha} = (m_{alpha+1} - K s) W^-1 [1 + O(W^-2)], increasing in m), (13.65), (13.66) (S then S), (13.67) (S then T); normal and abnormal arcs (volume I Definition 8.2.1) are the cases where m + K s is positive or negative; the case where the bifurcating arc ends is represented as m = 0 and is positioned left of the normal S-arcs and right of the abnormal ones (sequences 13.63 to 13.64, 13.68 to 13.70). Branches with different first arcs are ordered by Y_1 (13.71 to 13.75: 1,3,5,...;2 for K < 1; for 1 < K < 3 the order is 1;2;...7,5,3, for 3 < K < 5 it is 3,1;2;...9,7,5, and so on). Hierarchy ("packets", 13.3.2.5): dY_1 = O(W^(1 - 2 alpha)) (13.77; for alpha = 0 this is O(W), the first-arc separation): branches with a given first arc lie within O(W^-1) of each other (first-order packet), those sharing two arcs within O(W^-3), and so on. This is also the source of the numerical difficulty: branches with long common prefixes are exponentially close.
Worked example 1 < K < 3, n <= 5 (Figs 13.8 to 13.9). The "trident" configuration: four branches in one region, two symmetric (H_1, H_2), two asymmetric and exchanged by symmetry (H_3, H_4), resolved by the Chapter 7 argument (the H_1H_2 family is symmetric, the H_3H_4 family meets it at a symmetric orbit that is an extremum in W). In a partial type 1 bifurcation tridents exist only for odd n (volume I sect. 7.3.1.1).

### 4.6 Results (13.3.3, pp.60 to 78)

Junctions are given as figures (13.10 for n <= 3, 13.11 for n <= 4, 13.12 to 13.14 for n <= 5, 13.15 to 13.18 for n <= 6) and tables (Tables 13.1 to 13.10, same format as volume I Tables 8.4 to 8.11). The hand calculation went to n = 6 and was then automated, which verified it and extended it to n = 7 for -3 < K < -1 and 1 < K < 3, where all junctions are determined (printed in Tables 13.4 and 13.7). For n = 8 and the same two K ranges two groups of 4 branches cannot be resolved by the positional method. Printed numerical exceptions at n = 6 (READ p.67):
- For -1 < K < +1 a group of four branches +2112, +213, +33, +312 is not resolved by the positional method (and for +3 < K < +5, K > +5: -33, -321, -11121, -1113; for -5 < K < -3 and K < -5: +33, +3111, +12111, +123). The reason is that the junctions between the four change inside the K interval.
- **Change of junctions at K = -5.669369... (and at K = +5.669369... by Proposition 13.3.2, p.67).** Fig. 13.19 (printed 72) shows the characteristics at K = -5.66 and K = -5.68 (the W axis 8 to 9, Y_1 axis 1 to 2.5), with branches +33, +3111, +12111, +123. For K < -5.669369 the junctions are as in Table 13.1, for -5.669369 < K < -5 as in Table 13.2, otherwise as Table 13.3.
- **Change at K = 0 for -1 < K < +1** (Fig. 13.20 at K = -0.01, 0, +0.01; axis ranges as printed in the figure, not read): at K = 0 the symmetry of Sigma for even n inverts the sequences of Y_i and X_i; the branches +2112 and +33 are invariant, +213 and +312 exchanged, so at K = 0 there is a trident with all four branches meeting at a common symmetric solution; for K not zero they pair off. The junctions are as in Table 13.5 for -1 < K < 0 and Table 13.6 (Table 13.5 with two changes) for 0 < K < 1. No other change was found inside -1 < K < +1.
Positional method fails exactly where the two sets of junctions differ at the same normal/abnormal pattern, and the figures are numerical (Henon's program with a relaxation method for the asymptotic T/S decomposition, then shooting).

## 5. Chapter 14: total bifurcation of type 1 (pp.79 to 91)

System (14.1): the same equations with i modulo n, 2n equations in 2n + 1 variables W, Y_0..Y_{n-1}, X_1..X_n; one studies K >= 0 and Y_0 > 0 by the two isomorphisms (volume I sect. 8.5.1), with the origin in P. The properties are as in chapter 13 with (14.2) adding the periodicity conditions f_{2n+1} = Y_0 - Y_n = 0, f_{2n+2} = X_1 - X_{n+1} = 0.
- **Jacobian and stability (14.1.1, 14.1.2; pp.79 to 81):** for a total bifurcation the stability index z can be computed since the whole orbit is known, z = (1/2) trace of the 2 x 2 monodromy-like matrix (14.5, 14.6), and **(14.11): z = 1 - |J|/2. Proposition 14.1.1: the Jacobian vanishes if and only if the orbit is a critical orbit of the first kind (stability index z = 1).** For |W| tending to infinity (14.12): z = u_1 u_2 ... u_N / (2 Y_{i_1}^2 ... Y_{i_N}^2)[1 + O(W^-2)], z = Theta(W^(2N)), strongly unstable, sign(z) = (-1)^zeta where zeta is the number of T-arcs and abnormal S-arcs (14.13).
- **n = 2 (1T2)** (14.14, 14.15): the two equations give Y_0 = Y_1 and **W = 2 Y_0 + (1 - K^2)/(2 Y_0), X_1 = -(1 + K)/(2 Y_0)** (READ on image). Hyperbola with asymptotes Y_0 = 0 and Y_0 = W/2 (Fig. 14.1); branch +E (first species orbit) is asymptotic to Y_0 = W/2, branches +-11 to Y_0 = 0. For 0 <= K < 1 an extremum in W at **Y_0 = sqrt(1 - K^2)/2, W = 2 sqrt(1 - K^2)** (14.16). Stability index (14.17) is partly legible only (not transcribed); on the +E branch Y_0 tends to infinity and z tends to 1 from below: the first species orbit is stable, in agreement with Bruno 1976 and 1994 chapter VII). The 1T2 bifurcation was already studied by Guillaume 1971 pp.112 to 119.
- **n = 4 (1T4)** (Figs 14.2 to 14.4 for K = 0.3, 2.2, 4): the figure is qualitatively constant inside 0 <= K < 1, 1 < K < 3, 3 < K, so the junctions are determined. Special orbits (14.18) to (14.21), READ on the image:
  - Omega_2: W = (3 + K) sqrt((1 - K)/2), Y_0 = sqrt((1 - K)/2), X_1 = -(1 + K)/sqrt(2 (1 - K)) (exists only for 0 <= K < 1);
  - Omega_3: W = (3 - K) sqrt((1 + K)/2), Y_0 = sqrt((1 + K)/2), X_1 = -sqrt((1 + K)/2);
  Omega_2 and Omega_3 are intersections of an n = 4 family with the n = 2 family described twice, also an extremum in W of the n = 4 family; z = +1 for the n = 4 family, z = -1 for the n = 2 family there.
  - Omega_4: W = (3 + K)/sqrt(2), Y_0 = (1 + sqrt(K))/sqrt(2), X_1 = -(1 - sqrt(K))/(sqrt(2) (1 + sqrt(K)));
  - Omega_5: W = sign(1 - K) (3 + K)/sqrt(2), Y_0 = |1 - sqrt(K)|/sqrt(2), X_1 = -(1 + sqrt(K))/(sqrt(2)|1 - sqrt(K)|) (the printed extra factor structure of X_1 is read from the image);
  Omega_4 and Omega_5 exist for all K, are intersections of two n = 4 families (one symmetric, one asymmetric: another trident), also an extremum in W for one of them, and have z = +1 for both.
  I COMPUTED the equations (14.1) for n = 4 at K = 0.5, propagating from the printed (W, Y_0, X_1) and testing the periodicity conditions Y_4 = Y_0, X_5 = X_1: all four orbits close to better than 1e-13 (Omega_2: Y_i = 0.5 for all i; Omega_3: Y_i = 0.866025; Omega_4: Y = 1.207107, 1.207107, 0.207107, 0.207107; Omega_5: Y = 0.207107, 0.207107, 1.207107, 1.207107). So the printed forms are internally correct at that K.
- **n = 6** (14.2.4): numerical, from printouts of the branches; the figure is qualitatively constant in 0 <= K < 1, 1 < K < 3, 3 < K < 5; for K > 5 some junctions change at **K = 5.2612...** (the text truncates the digits; Fig. 14.5 on the image shows K = 5.25 and K = 5.27 for branches +-1131-, +-111111-, +21111, +231, W axis 5.8 to 6, Y_0 axis 2 to 3). (The text layer prints the caption as "K5.25 and K = -5.27"; the image shows both are positive. Printed typo in the text layer only, none on the page.)
- **Method (14.2.2):** at large W relaxation using the T/S decomposition (arc solved as a linear system, nodes from the encounter equation), followed downward in |W| until the decomposition ceases to converge, then trial-and-error shooting with Newton-Raphson on (Y_0, X_1) with the 2 x 2 monodromy matrix (14.5) giving z. No equivalent of the positional method was found for total bifurcations.
- **Conclusions (14.3):** the type 1 study is complete; all junctions are determined for n <= 6 (Tables 13.1 to 13.10 and 14.1 to 14.5), agree with volume I wherever both exist, and all the cases left undecided in volume I chapter 8 up to n = 6 are now solved. "Nothing prevents in principle the solution of higher values of n, using numerical computation. However, the amount of work grows exponentially."

## 6. Chapter 15: the Newton approach (pp.93 to 129)

Prompted by Bruno (1998, 2000) (the book "Local Methods in Nonlinear Differential Equations" chapter 1 and the nonlinear-equation chapter 2); the chapter was added after the monograph was otherwise complete. The aim is a systematic derivation of the chapter 12 results from the Newton polyhedra of the equations (12.32), and the author says at the outset that it "is not applicable to the general case" and "in practice it can only be used for small values of n". The variables are x_1 = mu, x_2 = Delta C, x_3 = Delta a_1, x_4 = y_1, x_5 = Delta a_2, ..., x_{2n+1} = Delta a_n (15.2), n_B = 2n + 1 variables and m_B = 2n - 1 equations after eliminating y_0 and y_n (15.4), so m_B = n_B - 2 and the solutions lie on two-dimensional manifolds (one for mu, one for the family parameter); the commonest case in Bruno's book has m_B = n_B - 1 (curves). Solutions are sought as powers x_i = b_i tau^(p_i) (1 + o(1)) with tau tending to infinity and all p_i < 0 (the cone of the problem, 15.15, 15.16) and
**nu = p_2/p_1** (15.14), tying the exponent of chapter 12 to the exponents of the Newton method.
- **Encounter equation:** support of 3 points, Newton polyhedron a triangle (1 face of dimension 2, three of dimension 1; 4 normal cones labelled a, c, d, b, 15.25 to 15.30). **Arc equation (general):** 4 points, a tetrahedron, 11 faces (1 + 4 + 6) with 11 normal cones (labels a, d, c, e, b, h, k, j, g, f, ...; 15.31 to 15.36); initial and final arcs have 3 points and 4 normal cones each (15.37 to 15.47). **Additional relations (12.33):** 4-point tetrahedra again (labels A to K, 15.48 to 15.51), with special first, last and n = 2 cases (15.52 to 15.66).
- Intersections with the cone of the problem are non-empty for all faces (15.4). A "coherent boundary subset" is a combination of one normal cone per equation with non-empty intersection (the cone of truncation); the number of possible combinations is 4^(n+1) 11^(n-2) (15.68), that is 64, 2816, 123904 for n = 2, 3, 4, which "quickly becomes impracticable". The intersection is found by the Motzkin-Burger algorithm (Theorem 15.5.1 of the book, from Bruno 1.4 Theorem 4.1; Theorem 15.5.2 for equalities), with strict inequalities replaced by non-strict ones and the parasitic solutions then removed (15.5.2).
- **Table 15.1 (printed p.111, READ), number of valid combinations for 1Pn and the computing time on an HP 720 workstation:**

| n | valid combinations | seconds |
|---|---|---|
| 2 | 12 | 0 |
| 3 | 39 | 2 |
| 4 | 138 | 24 |
| 5 | 505 | 1143 |
| 6 | 1920 | 69600 |

- **1P2** (15.5.4, 15.6, 15.7): n_B = 5, m_B = 3; 64 possible, 12 valid combinations (Table 15.2, READ on the image, p.112; the entries are the vectors N^i, "-2 -1 -1 -1 -1" being (-2, -1, -1, -1, -1), which is present in every row; columns are Case, i, N^i, dim Pi, d):
  aaa: N1 = (-2,-1,-1,-1,-1), dim 1, d 4.
  acc: N1 = (0,0,0,0,-1), N2 = (-2,-1,-1,-1,-1), dim 2, d 3.
  aba: N1 = (-1,0,0,0,0), N2 = (-2,...), dim 2, d 3.
  dad: N1 = (-1,0,0,-1,0), N2 = (-2,...), dim 2, d 3.
  dbd: N1 = (-1,0,0,0,0), N2 = (-1,0,0,-1,0), N3 = (-2,...), dim 3, d 2.
  cac: N1 = (-1,0,-1,0,-1), N2 = (-2,...), dim 2, d 3.
  ccc: N1 = (0,0,0,0,-1), N2 = (-1,0,-1,0,-1), N3 = (-2,...), dim 3, d 2.
  cda: N1 = (0,0,-1,0,0), N2 = (-2,...), dim 2, d 3.
  cdc: N1 = (0,0,-1,0,0), N2 = (-1,0,-1,0,-1), N3 = (-2,...), dim 3, d 2.
  cbc: N1 = (-1,0,0,0,0), N2 = (-1,0,-1,0,-1), N3 = (-2,...), dim 3, d 2.
  bab: N1 = (0,-1,0,0,0), N2 = (-2,...), dim 2, d 3.
  bbb: N1 = (-1,0,0,0,0), N2 = (0,-1,0,0,0), N3 = (-2,...), dim 3, d 2.
  The truncated systems are Table 15.3. Four cases (cac, ccc, cdc, cbc) are degenerate (first and third equations identical); the remedy (Bruno 2.6 Remark 6.2) replaces the last cone c by the cone A for the additional relation (15.8), giving the two valid cases caA and cbA (Tables 15.4, 15.5).
- **Power transformations (15.7).** For each of the 10 resulting cases a monomial change of variables x_i = prod w_j^(theta_ij) reduces the truncated system to d variables. Results (all READ, pp.115 to 126):
  - aaa (d = 4): a one-parameter family w_1 = (w_3/G2) + (1 - K^2) G1 G3/(2 G2 w_3), etc., giving the equations of the nu = 1/2 characteristic: **W = Y_1 + (1 - K^2)/(2 Y_1)** (15.107, reading the text layer; the same equation as (13.41), with X_1 and X_2 equal to the values obtained from (13.1a), up to the sign conventions of (12.112)). With nu = 1/2 (15.106). Two points of the characteristic at W = 0 for |K| > 1 are not covered because all variables must be non-zero; they come from bab.
  - acc, aba, cda, dbd, cbA, bbb: no solution (the system forces a variable to zero, contradicting (15.10), or is inconsistent).
  - dad (d = 3): Jacobian determinant -2 G2 G3 (non-zero, no critical point), solution w_i = psi_i (15.117, 15.118); with the vectors N1 = (-1,0,0,-1,0), N2 = (-2,-1,-1,-1,-1) and parasitic elimination lambda_1, lambda_2 > 0, nu = lambda_2/(lambda_1 + 2 lambda_2), so **0 < nu < 1/2** and the physical results (15.123, 15.136) with error terms equal to (12.78)-type and the S-arc of two basic arcs of (12.110) and (12.108): Delta a_1 = [G2/(G3 (1 - K))] Delta C [1 + ...], y_1 = [G1 G3 (1 - K^2)/(2 G2)] mu Delta C^-1 [1 + ...], Delta a_2 = -[G2/(G3 (1 + K))] Delta C [1 + ...]. The error terms of chapter 12 are recovered exactly.
  - caA (d = 3): the T-arc of (12.109), 0 < nu < 1/2: Delta a_1 = ((K + 1) G1/(2 G2)) mu Delta C^-1, y_1 = G2 Delta C, Delta a_2 = ((K - 1) G1/(2 G2)) mu Delta C^-1.
  - bab (d = 3): two solutions if |K| > 1, w_3 = ±sqrt(G1 (K^2 - 1)/(2 G3)) etc.; N1 = (0,-1,0,0,0), nu = (lambda_1 + lambda_2)/(2 lambda_2) so **nu > 1/2** (15.164, 15.165); it is the case nu > 1/2 of chapter 12.6, and the physical form is (15.160, 15.166) Delta a_1 = ±(sqrt(...)) mu^(1/2) etc. This is the solution that carries the two intersections of the 1P2 characteristic with W = 0 for |K| > 1.
  So of the 12 combinations two give the first-species-or-T-or-S-arc regimes of 0 < nu < 1/2 (dad, caA), one the nu = 1/2 characteristic (aaa), one the nu > 1/2 points (bab), and the rest are empty.
- **Total bifurcation (15.8):** variables n_B = 2n + 2, m_B = 2n (15.177), 44^n possible combinations; Table 15.6: valid combinations for 1Tn.

| n | valid combinations | seconds |
|---|---|---|
| 2 | 32 | 0 |
| 3 | 85 | 10 |
| 4 | 300 | 487 |
| 5 | 1095 | 26079 |

  1T2 gives 11 x 4 x 11 x 4 = 1936 possible and 32 valid combinations (Table 15.7 lists them); in the case kbkb there are t = 5 vectors but the rank is 4 (N1 - N2 - N4 + N5 = 0). The vector (-2,-1,-1,-1,-1,-1) is always present. The 32 cases are not worked out.
- **Conclusions (15.9):** the approach confirms chapters 12 to 14 "at least in the simplest case 1P2" in a more rigorous way, but it must be done n by n, by enumeration, and the cost grows exponentially, so it is limited to 1P2. The objective of general results valid for any n was reached in chapter 12 instead.

## 7. Chapter 16: proving general results (pp.131 to 148)

A third approach: keep the exponent ansatz x_i = b_i tau^(p_i) (16.5) but argue algebraically for all n. Four cases of dominant terms in the encounter equation (16.8, labels a to d, corresponding to the normal cones of 15.30). Results (all READ):
- **Proposition 16.3.1:** p_{2i+2} <= max(p_2, p_1/2) for every encounter i in C (the lead y_i cannot be larger than max(Delta C, mu^(1/2))). Proof by contradiction: if violated at the first i, the neighbouring arc and encounter equations force all following exponents equal and coefficients alternating, and the chain ends in the final arc where the arc equation cannot balance (partial bifurcation); for a total bifurcation, summation of the additional relation over every other i gives G3 n b = 0, contradicting (15.10).
- **Proposition 16.3.2:** p_{2i+1} <= max(p_2, p_1/2) for every arc i. Together they say everything depends on which of p_2 and p_1/2 is larger, i.e. on nu = p_2/p_1 against 1/2.
- **p_2 = p_1/2 (nu = 1/2)** (16.4, Proposition 16.4.1: p_{2i+2} = p_1/2 for any i): re-derives (12.114) exactly with the changes of variable (16.29), (16.30).
- **p_2 < p_1/2 (nu > 1/2)** (16.5): re-derives (12.119).
- **p_2 > p_1/2 (nu < 1/2)** (16.6): "nodes*" are encounters with p_{2i+2} < p_2 (coefficient zero), "antinodes*" those with p_{2i+2} = p_2, "arcs*" the parts between nodes*. It is shown that a total bifurcation without nodes* is a first species orbit with y_i = -(-1)^i (G2/2) Delta C and Delta a_i = [G1/((K - (-1)^i) G2)] mu Delta C^-1 (16.52, 16.53, the leading terms of 12.78). An arc* is an S-arc* of m basic arcs (m odd, 16.59, 16.61, reproducing 12.110) or, for m = 2, a T-arc* (16.62); even m > 2 is impossible. Nodes* reproduce the three node formulas (12.106) to (12.108) (16.70 to 16.72). T-arcs* are refined to (12.109) (16.82), with the particular cases g' = 1, m_b = 1 or g'' = 1, m_a = 1 where the coefficient of Delta a_{i+3} (or Delta a_{i+5}) vanishes and "cannot be determined without further computation".
- **Proposition 16.6.1: there cannot be two T-arcs* in succession (no TT node*).** Proved in the appendix 16.8 (pp.143 to 148) for a partial T-sequence and a total T-sequence, in three cases p* greater than, equal to, and less than p_1/2 (p* = max of the exponents of the two arcs of a T-arc*): in each case sums of the arc and encounter equations give a positive-definite sum equal to zero (for example 2 G3 sum b^2 = 0, 16.105, 16.127), or in the third case a single-T-arc partial sequence. This is the perturbed (mu > 0) analogue of Proposition 4.3.2 and is more involved than the mu = 0 proof of chapter 12.2.
- **Conclusions (16.7):** all basic results of chapter 12 are rederived "in a completely independent way" using only the fundamental equations and without recourse to the qualitative analysis of volume I; in particular the existence and properties of nodes, antinodes, S-arcs and T-arcs are derived. The method "is an ad hoc method", strongly dependent on the specific form of the equations, not generalisable by default.

## 8. Printed numbers, formulas and tables usable as sourced tests

### 8.1 Closed-form values (all READ on the page image; test the formula at any K in the stated range)

| Item | Printed form | Page, equation | Note |
|---|---|---|---|
| Type 1 constants | G1 = 4a^2/v, G2 = 2 pi (dZ/dC)_A, G3 = 3 pi J sqrt(a)/2, K = 1 + (2/J)[(dZ/dA)_C - beta_0] | 21, (12.31) | G1, G3 > 0; Z from volume I (4.28) |
| 1P2 characteristic | 2 Y_1 + (1 - K^2)/Y_1 - 2W = 0 | 47, (13.41) | check: at K = 0.5 the extremum has Y_1 = 0.612372, W = 1.224745 (COMPUTED from 13.42) |
| 1P2 extremum | Y_1 = sqrt((1 - K^2)/2), W = sqrt(2 (1 - K^2)), |K| < 1 | 48, (13.42) | COMPUTED: the printed pair satisfies (13.41) to 4e-16 at K = 0, 0.5, -0.5 |
| 1P3 symmetric family | 2W = (3 - K) Y_1 + (1 - K^2)/Y_1 | 48, (13.45) | |
| 1P3 asymmetric family | W = Y_1 + (1 - K)/Y_1, Y_1 Y_2 = 1 - K | 48, (13.46) | |
| 1P3 intersection | Y_1 = sqrt(1 - K), W = 2 sqrt(1 - K), K < 1 | 48, (13.47) | COMPUTED: lies on both (13.45) and (13.46) (residual 0 at K = 0) |
| 1P3 critical arc | Y_1 = sqrt((1 - K^2)/(3 - K)), W = sign(3 - K) sqrt((1 - K^2)(3 - K)), -1 < K < 1 or K > 3 | 50, (13.50) | COMPUTED: at K = 0.5, Y_1 = 0.547723, W = 1.369306, lies on (13.45) |
| 1T2 characteristic | W = 2 Y_0 + (1 - K^2)/(2 Y_0), X_1 = -(1 + K)/(2 Y_0), Y_0 = Y_1 | 82, (14.15) | |
| 1T2 extremum | Y_0 = sqrt(1 - K^2)/2, W = 2 sqrt(1 - K^2), 0 <= K < 1 | 82, (14.16) | COMPUTED: minimum of (14.15) |
| 1T4 orbits Omega_2 to Omega_5 | see 5 | 86, (14.18) to (14.21) | COMPUTED at K = 0.5: all close under (14.1) to 1e-13 (section 5) |
| Stability index, total type 1 | z = 1 - |J|/2 | 81, (14.11) | |
| Large-W stability index | z = u_1...u_N / (2 Y_{i_1}^2 ... Y_{i_N}^2) [1 + O(W^-2)] | 82, (14.12) | order Theta(W^(2N)) |
| Number of possible combinations | 4^(n+1) 11^(n-2): 64, 2816, 123904 for n = 2, 3, 4 | 107, (15.68) | COMPUTED: arithmetic matches |
| Possible combinations, 1Tn | 44^n (1936 for n = 2) | 127, (15.182) | |
| Newton method counts | Tables 15.1 and 15.6 (section 6) | 111, 127 | counts and HP 720 seconds |

### 8.2 K breakpoints and boundary values (READ)

| Value | Meaning | Page |
|---|---|---|
| K = -5.669369... (and +5.669369... by Proposition 13.3.2) | 1P6 junctions change here | 67, 73 (Tables 13.1, 13.2, 13.9, 13.10) |
| K = 0 | 1P6 junctions of +2112, +33, +213, +312 change here | 67 |
| K = 5.2612... | 1T6 junctions change here (digits truncated by the author as printed) | 86 |
| K = integer | excluded by Conjecture 12.1.1 | 22 |
| K in (-3, -1), (1, 3) | the two ranges taken to n = 7 | 73 |
| |K| < 1: W_min = sqrt(2(1-K^2)) for 1P2 | smallest |W| on the family | 48 |

The number K = 5.669369... is a candidate positive control for any independent implementation of (13.1) at n = 6: the junctions of the four branches +33, +3111, +12111, +123 switch there. I did not recompute it (follow-up).

### 8.3 Junction tables

Tables 13.1 to 13.10 (printed pp.73 to 78) and 14.1 to 14.5 (pp.88 to 91) list, for each n and K interval, the families as groups of branches joined together; each branch is named by its sign (that of W) and the sequence of numbers of basic arcs in each T or S arc (a "2" is a T-arc, "1", "3", "5" are S-arcs with that number of basic arcs, "11" two S-arcs of one basic arc, etc.). The format is that of volume I Tables 8.4 to 8.11 and 8.14 to 8.17, and the book states they agree with volume I wherever both exist. Transcribed from the page image here, as examples:
- Table 13.1 (K < -5.669369...): group 1P6+++++A: -15, -1311; group 1P6++---A: +33, +3111, +123, +12111 (others as Table 13.3).
- Table 13.2 (-5.669369... < K < -5): group 1P6+++++A: -15, -1311 (others as Table 13.3).
- Table 13.9 (5 < K < 5.669369...): group 1P5++++S: -5, -11111; group 1P6+++++A: -51, -1131 (others as Table 13.8).
- Table 13.10 (K > 5.669369...): 1P5++++S: -5, -11111; 1P6+++++A: -51, -1131; 1P6++++--A: -33, -1113, -321, -11121 (others as Table 13.8).
- Table 13.8 (3 < K < 5), read from the 90 dpi image (families as printed, one family per class label): 1P2+A: +2, -11. 1P3++S: -3, -111. 1P3+-A: +21, -12. 1P4+++A: -31, -1111. 1P4+--A: +211, -112, -13, -121. 1P5++++S: +5, -11111. 1P5++++A: -311, -113. 1P5+++-A: -32, -1112. 1P5+--+S: +212, -131. 1P5+--+A: -1211, -1121. 1P5+---A: +23, +2111. 1P6+++++A: +51, -1131, -3111, -111111. 1P6+++--A: -33, -321, -1113, -11121, -312, -11112. 1P6+--++A: +213, +2121, +2112, -11211, -1311, -12111. 1P6+----A: +231, +21111, -15, -123, -132, -1212. The reading is of a small page image; a transposed digit is possible, so re-read the page before using any single entry as a test (the first and last rows of each family are the least certain). The underlines printed in the branch labels (which mark something about the arc sequence) are not explained in the text I read and are omitted.
- The other tables (13.3 to 13.7, 14.1 to 14.5) are for the junction groups of 1P2 to 1P7 and 1T2 to 1T6 at the K ranges named in their captions: Table 13.3 for -5 < K < -3, 13.4 for -3 < K < -1 (n <= 7), 13.5 for -1 < K < 0, 13.6 for 0 < K < 1, 13.7 for 1 < K < 3 (n <= 7), 13.8 for 3 < K < 5; 14.1 for 0 <= K < 1, 14.2 for 1 < K < 3, 14.3 for 3 < K < 5, 14.4 for 5 < K < 5.26..., 14.5 for 5.26... < K (n = 2 to 6). Not transcribed (the text-layer layout is scrambled and a transcription would need image reading of each page: follow-up).

### 8.4 What the book does NOT print in this part

No table of numerical orbits, no value of mu, no comparison of the asymptotic formulas with computed orbits at a stated mu, and no value of G1, G2, G3 or K for any actual bifurcation orbit (K values are in volume I sect. 8.2.1 and Hitzl and Henon 1977). The "quantitative comparison with numerically found families" promised in the preface is, in chapters 11 to 16, a qualitative agreement of the junctions with the volume I invariants and with numerically computed characteristics in the scaled (W, Y) plane (figures 13.3 to 13.7, 13.19, 13.20, 14.2 to 14.5), not a comparison at a stated mu. Part B should be checked for any such comparison.

## 9. #906: how the turn and periapsis scale with mu near a type 1 bifurcation

The question: the project's demanded-turn gate must treat as "indeterminate" those encounters whose demanded turn is small near a low-order resonance (volume I and Hitzl-Henon: 64 of 160 non-degenerate critical generating orbits demand under 2 degrees). What does volume II say about the scaling of the turn and the periapsis with mu?

**Answer (READ for the relations used; DERIVED for the combination):** near a type 1 bifurcation orbit, every encounter has a periapsis distance and a turn controlled by the lead y_i and the exponent nu = log|Delta C|/log mu of the distance to the bifurcation, with these scalings (mu small, impact parameter d_i = |y_i| |sin phi_i|, theta_i = 2 mu/(v^2 d_i) from 11.23):

| regime | node: periapsis d, turn theta | antinode: periapsis d, turn theta |
|---|---|---|
| Delta C finite (nu = 0) | d = Theta(mu/Delta C), theta = Theta(Delta C) from y_i (12.106 to 12.108) | d = Theta(Delta C), theta = Theta(mu/Delta C) = Theta(mu) |
| 0 < nu < 1/2 | d = Theta(mu^(1 - nu)), theta = Theta(mu^nu) | d = Theta(mu^nu), theta = Theta(mu^(1 - nu)) |
| nu = 1/2 (|Delta C| about sqrt(mu)) | d = Theta(mu^(1/2)), theta = Theta(mu^(1/2)) | the same |
| nu > 1/2 | d = Theta(mu^(1/2)), theta = Theta(mu^(1/2)) | the same |

Explicitly, using (12.106) and (11.23), at a node with a T-S junction and a finite impact angle, to leading order,

theta_node = 2 G2 |Delta C| / ( v^2 |m_a s_i - K| G1 G3 |sin phi| ) = G2 |Delta C| / ( 3 pi J a^(5/2) v |m_a s_i - K| |sin phi| )

(DERIVED: this is the printed y_i substituted in the printed (11.23) with d = |y| |sin phi|; mu cancels, so at a node the turn is independent of mu and proportional to the distance Delta C from the bifurcation orbit; it is the demanded turn, not a bound, and it applies while it is small, so the small-angle form is valid). At nu = 1/2, the leading statement from (12.115) is d_i = sqrt(G1 G3 mu) |Y_i| |sin phi_i| and theta_i = 2 sqrt(mu) / ( v^2 sqrt(G1 G3) |Y_i| |sin phi_i| ), with Y_i = Theta(1) a root of (13.1) (for example Y_1 = sqrt((1 - K^2)/2) at the 1P2 extremum, giving the smallest periapsis reached on that family, d_1,min = sqrt(G1 G3 mu (1 - K^2)/2) |sin phi_1|, and the largest turn there). The consequences for the gate:

1. **The turn is small near a bifurcation, and the smallness is a power of mu**: theta = Theta(mu^nu) at a node (which is where the large encounters are) and the whole family stays in theta = O(mu^(1/2)) once |Delta C| is below the order of sqrt(mu) (the family cannot be followed to |Delta C| < sqrt(2 (1 - K^2) G1 G3 mu)/|G2| for 1P2 and |K| < 1, W_min above). At the Earth-Moon mass mu = 0.0121529529 (the value in the Hitzl and Henon digest), sqrt(mu) = 0.110240 and mu |ln mu| = 0.053597 (COMPUTED); both are order 0.1, i.e. the asymptotic theory is at best marginal for the Earth-Moon system.
2. **No turn threshold is natural: the turn tends to zero continuously along the family as Delta C tends to the bifurcation.** A rule "turn near zero is not an encounter" would discard exactly the generating-orbit seeds near a bifurcation; this is the physical reason for the Hitzl and Henon #906 amendment (return "indeterminate", never "reject").
3. **A scale-free criterion exists**: classify a node by Delta C / sqrt(mu) (regimes nu < 1/2, = 1/2, > 1/2: the transition is the same in every type 1 bifurcation), and report the first-order periapsis and turn of the nearest bifurcation (the formulas above) when it is within 1 order of magnitude of sqrt(mu). The same quantity r_p / sqrt(mu) already appears as the second validity number of the Hitzl and Henon note; this section shows that r_p/sqrt(mu) = Theta(1) is the nu = 1/2 transition itself. So that number is not arbitrary.
4. **The first-order theory breaks down by mu^(1/2) in the leads** (the O(mu) corrections to (12.115) are relative O(mu^(1/2))): at the Earth-Moon mass sqrt(mu) = 0.11, so about 11 percent relative accuracy at best (INFERRED), which supports the earlier remark that the first-order theory is outside its validity for most Table I orbits.
5. **Relation to Perko's nu.** Perko 1977 and 1981 (digest above) use the exponent nu_P of the near-Moon distance, Delta_1 = O(mu^(nu_P)), with deflection O(mu^(1 - nu_P)); this is the exponent nu_i of (11.27) here, not Henon's nu (which is the exponent of Delta C). The dictionary of section 3.7 gives: Perko's case i (0 < nu_P < 1/2) is Henon's antinode with nu_P = nu_H; Perko's case ii (1/2 < nu_P < 1) is Henon's node with nu_P = 1 - nu_H; Perko's nu_P = 1/2 is Henon's nu_H = 1/2. Perko's variation dv_1 = O(mu^(nu_P)) (case i) and O(mu^(1 - nu_P)) (case ii, "pointed out ... by Dr Michel Henon") agree with the lines of the table: Delta C = mu^(nu_H) at the antinode and Delta a = O(mu/Delta C) at the node (INFERRED; checked qualitatively against the Perko digest, not numerically).

What the book does not give: no printed bound on the turn at a stated mu, no worked evaluation of G1, G2, G3, K for any real orbit (K is in volume I and in the Hitzl and Henon tables), so the coefficient of the turn law needs (dZ/dC)_A and (dZ/dA)_C from the project's own (A, Z) surface (volume I (4.28)), which `src/` does not contain.

## 10. Do the chapters supply the coefficients of Guillaume 1973's hyperbola and cubic?

The question (from the Guillaume 1973 and 1975 digests): the 1973 conference paper writes (3'') (alpha dx0 + beta dydot0)(alpha' dx0 + beta' dydot0) = gamma mu and (4'') dydot0 (a dydot0 + b dxdot0)(a dydot0 - b dxdot0) = delta mu without the Greek and Latin coefficients, and neither Guillaume 1975 paper prints them.

**Answer:** only partly, and only for the hyperbola of a type 1 bifurcation, and not in Guillaume's variables.
- The type 1 characteristic in the plane of the lead and the Jacobi-constant displacement is the hyperbola (13.41), in physical variables (DERIVED from 13.41 and 12.115 with s_1 = -1, signs of cos(phi) absorbed): y_1 (y_1 - G2 Delta C) = (K^2 - 1) G1 G3 mu / 2 + O(mu^(3/2)), a product of two linear forms (the asymptotes y_1 = 0 and y_1 = G2 Delta C) equal to a constant times mu. So the structure of Guillaume's (3'') is reproduced with an explicit constant, gamma_eff = (K^2 - 1) G1 G3 / 2 in the variables (y_1, Delta C), and with explicit expressions for the asymptote slopes (G2). The error term order is the same as the book's nu = 1/2 error O(mu^(1/2)) relative.
- The coefficients are expressed in terms of G1, G2, G3, K of (12.31), which need (dZ/dC)_A, (dZ/dA)_C and beta_0 of the (A, Z) surface of volume I. Henon 2001 does not tabulate them for any orbit. A numerical evaluation requires the project's own computation of Z(A, C) (volume I (4.28)).
- Guillaume's variables are the displacements of (x0, ydot0) of the initial state at the orthogonal crossing, not (y_1, Delta C). The map from (y_1, Delta C) to (dx0, dydot0) requires the orbit's Jacobian, not given. So Guillaume's alpha, beta, alpha', beta', gamma are not delivered.
- Guillaume 1973 example 2 family I (his "(3'')") is a type 1 bifurcation (two collisions per period on a T-arc, INFERRED; in the Guillaume 1973 digest). Its cubic (4'') is a different bifurcation: three basic characteristics meeting (one T-conic family and two C24 arcs), not a type 1 bifurcation. The book's type 2 chapters (nu = 1/3 transition, part B) are the place to look for a cubic; those chapters are not covered here, and the cross-reference is the part B digest. Not claimed.
- Guillaume's example 1 (retrograde unit circle at C = -1) is type 3 and is not treated in the book at all ("its analysis had not yet been completed", preface).

## 11. Relation to the other held sources

- Volume I: the type 1 results of this part are, per the book, consistent with every case of volume I chapter 8 up to the order where both exist, and settle all cases left open up to n = 6 (14.3). The symmetry, Restriction 7.3.1, Broucke's principle and the tridents of volume I chapter 7 are re-derived as consequences (section 4.1 and 4.5 above).
- Hitzl and Henon 1977: the "critical bifurcating arc" (13.1.4) is the arc at an extremum of C along its family, which is the arc family critical member tabulated there. Proposition 13.1.1 ties its location to a vanishing Jacobian and an infinite-to-minus-infinity stability index jump (Henon and Guyot 1970).
- Henon 1968 and Bruno 1981: no new formula; the collision velocity v used here is the volume I (8.3) quantity (the Bruno W), not new.
- Brjuno 1978: the mu = 0 background (arcs T_N and S, intersection types); this part only uses its results through volume I.
- Perko 1976b, 1977 and 1981 and Guillaume 1971, 1975: the intermediate-arc idea is Guillaume 1971 p.72 and 1975b p.452 and Perko 1977a p.277 (printed p.7); the encounter relation is Guillaume 1971 p.83; the nu = 1/2 hyperbola of Perko 1981 eq. 9 is the type 1 case n = 2 of (13.41) (INFERRED, structure only). Chapter 12 generalises these to symmetric or asymmetric orbits and to partial or total bifurcations (p.1).
- Bruno 1998 and 2000: chapter 15 uses Bruno's Newton-polyhedron method (chapters 1 and 2 of that book; I did not hold it, so the Theorem and chapter citations are as printed by Henon).

## 12. Reconciliation with project code (grep of `src/cyclerfinder` and `scripts`, 2026-10-04)

- `src/cyclerfinder/verify/turn_gate.py` is the demanded-turn gate (#888): it computes the demanded turn between incoming and outgoing V-infinity, the available bend `2 asin(1/(1 + r_p v^2/GM))`, the ratio, and the required altitude. It has no scaling in mu, no reference to sqrt(mu), no "indeterminate" verdict for small turns, and does not use a bifurcation or second-species orbit as input. Its relation r_p = (GM/v^2)(1/sin(turn/2) - 1) is the exact hyperbola relation, consistent with (11.23) to first order in the turn (theta = 2 mu/(v^2 d) for small theta; the exact tan(theta/2) = mu/(v^2 d) is Mihalas and Routly's, not printed in the book). No contradiction. The amendment in #906(b) (indeterminate near resonance, report first-order periapsis, never reject) is supported, and its two validity numbers, mu |ln mu|/v^3 and r_p/sqrt(mu), are respectively the logarithmic-factor neglect of (11.40) and the nu = 1/2 transition of (12.115) (section 9).
- `src/cyclerfinder/search/bifurcation_detector.py`, `mu_continuation.py`: Floquet-multiplier root-of-unity detection and pseudo-arclength continuation in mu for symmetric periodic orbits at finite mu. This part's bifurcations are bifurcations of the generating (mu = 0) families, a different object: the in-code detector will not see a junction of Henon's type 1, because at finite mu the branches are joined two by two and the "intersection" is a fold or a stability index jump (Proposition 14.1.1, z = 1 - |J|/2). The Floquet-index crossing +1 found by `detect_saddle_center_bracket` is the same event as the vanishing Jacobian at an extremum in W (section 4.3), INFERRED.
- Second-species code: `search/earth_moon_class1_resonant_connections.py`, `earth_moon_resonant_families.py` and similar state that they do NOT implement a second-species differential correction (module docstrings); `genome/multi_shooting.py` records that no regularised transition matrix exists (task #928). There is no implementation of the (A, Z) plane, the G1 to G3, K constants, the lead y_i, or the scaled system (12.114). Nothing in the repository contradicts this part; nothing implements it.

## 13. Techniques applicable to the project's problems

Everything below is INFERRED use; nothing was built or run.

**#899 (second species seeds and continuation in mu to the Earth-Moon value).**
1. A reduced model for the small-mu end of the continuation. For a type 1 bifurcation the junction structure at small mu is the finite system (13.1), whose solution is a one-parameter curve in (W, Y_1) at fixed K for each n. A continuation driver following a family from mu = 1e-6 (the Casoliva seed mass) upward can use the scaled variables of (12.115) as a predictor: Delta C = -sqrt(G1 G3 mu) W/G2, y_i = s_i sqrt(G1 G3 mu) Y_i, Delta a_i = -s_i sqrt(G1/G3 mu) X_i. Test: at fixed K and n, a family followed in mu at constant W should keep its scaled variables (Y, X) approximately constant; deviation grows as mu^(1/2) (the stated error order).
2. Branch selection at a fold. The "policeman" problem of Guillaume 1973 is exactly the junction problem of chapters 13 and 14: Tables 13.1 to 13.10 and 14.1 to 14.5 give, for each (n, K interval), which branches join, and give the K values where the rule changes (5.669369 and 0 at n = 6, 5.2612 for 1T6). A driver that has computed K from the (A, Z) surface of the seed can use these to pick the branch it will land on after a fold at small mu and flag the ambiguity where K is within a breakpoint. Positive control: the symmetric-complement case with K < 1 for n = 3 (Fig. 13.10), where volume I's Restriction 7.3.1 and this part's positional method must give the same junctions.
3. Distance to a fold. The smallest |Delta C| at which a type 1 family can exist is |Delta C|_min = sqrt(2 (1 - K^2) G1 G3 mu)/|G2| for 1P2 (|K| < 1), with similar closed forms for 1P3 (13.47, 13.50); a continuation step control in mu can use these values to size its steps (a family with a fold at W_min will fold at |Delta C| proportional to sqrt(mu): a step that holds Delta C fixed while mu increases will cross the fold). This supports the Casoliva observation that fixed-period continuation in mu "in most cases" ends on lunar impact: at a node the periapsis falls as Delta C grows (d = Theta(mu/Delta C)), so a continuation that raises Delta C at fixed mu drives the node toward the Moon, which is a reading of section 9 and not a statement of the book.
4. K as a classifier. Because the junction pattern depends on K only through the integer part of |K| and the finitely many break values, a seed catalogue can store K per bifurcation (not per orbit), ordered by the junction regime.

**#906 (turn gate).** Section 9: the turn and the periapsis scale as mu^nu and mu^(1-nu) at nodes and as mu^(1-nu) and mu^nu at antinodes; nu = 1/2 is the unique scale-free regime; the gate should compute the nearest bifurcation orbit's first-order d and theta and report "indeterminate: near bifurcation, nu = log(|Delta C|)/log(mu) in (0, 1/2]" with those numbers, instead of a threshold. The exact thresholds on Delta C/sqrt(mu) are not delivered by the book; a defensible choice is the nu = 1/2 transition |Delta C| ~ sqrt(mu) (a labelled convention, not a sourced value).

**#928 (regularised propagator and transition matrix).** The intermediate-arc approximation (Proposition 11.3.1) is an alternative to integrating through the encounter at all: the orbit is O(mu) close to a sequence of Keplerian arcs joined at the encounters with the matching relations (11.65), (11.70). A continuation that only needs the O(mu) state at the next encounter can avoid a regularised propagator for the near-collision pass, and the first-order transition matrix across an encounter is available in closed form (the variational equations (13.18) to (13.24) amplify perturbations by Theta(W^(2 alpha - 2)) after alpha nodes: the book's estimate of the condition number of the monodromy matrix; they also say Delta p = O(mu) so the first-order map has error O(mu), compared to Heggie's regularised accuracy). This is a comparison point for the regularised transition matrix test, not a substitute. Positive control for any three-body version of the matrix: the Jacobian-versus-stability identities (13.1.1 and 14.1.1, z = 1 - |J|/2 for total bifurcations).

**#925 (elliptic-problem controls).** No direct content: volume II treats only the circular restricted problem. The matching and intermediate-arc scheme would extend to the elliptic problem with the Moon's motion a parameter (the book's (11.19) uses the circular motion of M2 explicitly in p_M); that extension is not in the book. Nothing here helps decide the Mako and Salamon xfails or the Neelakantan M4N2 miss; the only transferable item is the reminder that a published numerical figure should be checked against an equation set that the author prints (as I did for 14.18 to 14.21).

## 14. Defects, uncertainties and things not checked

- (12.115) omits sign(cos phi) which (12.112) and (12.117) include (printed slip; magnitude unaffected).
- The text layer of the page for Fig. 14.5 prints "K5.25 and K = -5.27"; the image shows both K positive.
- The 1T2 stability index (14.17) and the first derivatives (13.39) are partly illegible in the text layer; not transcribed. The packet exponent (13.77) was inferred from (13.54) and (13.56) (O(W^(1-2 alpha))), not read from the image.
- Transcribed from page images: Tables 13.1, 13.2, 13.8, 13.9, 13.10 (13.8 at 90 dpi, with the caveat given there). Tables 13.3 to 13.7 and 14.1 to 14.5 were not transcribed.
- The node formulas (12.106) to (12.108) were read from the page image and their common factor G1 G3/G2 is clean. A scaled intermediate result in the text layer, (12.96), appears to carry a factor 1/2 relative to the rescaled formula; I judge this to be the book's intermediate variable Y* before the rescaling (12.64) but did not redo the algebra.
- All DERIVED statements (the turn at a node, the physical hyperbola, |Delta C|_min) use the book's d = |y| |sin phi| (12.28) and (11.23), and the sign conventions of (12.112); I did not verify them numerically against a computed orbit. The constant gamma_eff in section 10 is up to a sign convention.
- I did not run any code and did not recompute the K = 5.669369 and 5.2612 break values.

## 15. Follow-ups (no task numbers registered)

1. Independent control for the junction rule: implement (13.1) and (14.1), reproduce 13.41, 13.45, 13.46, (14.18) to (14.21) and the n = 6 break values K = 5.669369 and K = 5.2612, then use them as an offline classifier for #899. The 14.18 to 14.21 residual check done here is already a start.
2. Transcribe Tables 13.3 to 13.7 and 14.1 to 14.5 from page images (each is a short junction list) into a machine-readable file with the K ranges, and check them against volume I Tables 8.4 to 8.17.
3. Compute G1, G2, G3, K for the project's type 1 bifurcation orbits (the first row of volume I Table 6.2, (I, J, L) = (2, 1, 0), C = -0.406767, and others) from the (A, Z) surface, and compare the sqrt(mu) law of section 9 against the minimum periapsis of a computed family at two or three masses. This is the numerical check that the book's chapter 11 to 16 do not give.
4. Read part B for any numerical comparison at a stated mu and for the type 2 cubic (Guillaume 1973 (4'')).
5. Decide the #906 gate wording around |Delta C|/sqrt(mu) (section 9) once item 3 gives the coefficient.
6. Look for volume I Tables 6.2 and 6.3 K values in the digest (volume I sect. 8.2.1) and pair them with the junction tables.
