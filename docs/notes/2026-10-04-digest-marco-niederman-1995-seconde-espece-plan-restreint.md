# Digest: Marco and Niederman 1995, construction of second species solutions in the planar restricted three-body problem

Date 2026-10-04. Purpose: establish exactly what this paper proves, so the project can word the #890 result (a periodic
orbit of the planar concentric circular restricted four-body problem, Uranus plus Titania plus Oberon, one close flyby
of each moon per cycle, found by continuing a patched-conic flyby chain in the moon mass; see
`docs/notes/2026-10-04-890-titania-oberon-candidate.md`) correctly against published second species theory. Related
digests, not repeated here: `docs/notes/2026-10-04-digest-bolotin-mackay-2000-second-species-n-centre.md`,
`docs/notes/2026-10-04-digest-bolotin-mackay-2006-nonplanar-second-species.md`,
`docs/notes/2026-10-04-digest-font-nunes-simo-2002-consecutive-quasi-collisions.md`,
`docs/notes/2026-10-04-digest-font-nunes-simo-2009-second-species-numerical-study.md`,
`docs/notes/2026-10-04-digest-bradley-russell-2014-patched-conics-to-full-gravity-continuation.md`.

Citation: Jean-Pierre Marco and Laurent Niederman, "Sur la construction des solutions de seconde espece dans le
probleme plan restreint des trois corps", Annales de l'Institut Henri Poincare, section A (Physique theorique),
62(3):211-249 (1995). Manuscript received 16 September 1993, revised 3 January 1995. Filed in the private paper corpus
as marco-niederman-1995-construction-solutions-seconde-espece-probleme-plan-restreint-trois-corps-ann-ihp-a-62-3-211-numdam-french.pdf.
The file has a usable text layer (about 10,600 words from pdftotext) but all displayed formulas are missing from it, so
every formula below was read from page images.

Evidence labels: READ (printed page, journal page numbers given), TRANSLATED (my English rendering of the French
original; all translations are mine), INFERRED (my reading of what the printed text implies, or a comparison with
project work). All 40 PDF pages were read: pages 211-249 as page images, with the text layer used as a cross-check.
PDF page n is journal page n+209 (PDF page 1 is the NUMDAM cover sheet). The scan is clean; no page was unreliable.
The paper has no table, and its figures (Figures 3, 4, 5) carry
no numbers beyond the energy thresholds listed in section 4.

## 0. Summary

READ (abstract, p. 211, French): "Apres une regularisation de Levi-Civita, les solutions de seconde espece apparaissent
comme bifurcations d'orbites homoclines relatives a un point fixe hyperbolique. L'idee est de J. Henrard, nous
ameliorons sa demonstration en utilisant un simple theoreme d'approximation autour du point fixe, et nous donnons dans
ce contexte de maniere analytique les conditions necessaires d'existence des solutions." TRANSLATED: after a
Levi-Civita regularisation, second species solutions appear as bifurcations of homoclinic orbits of a hyperbolic fixed
point; the idea is Henrard's; the authors improve his proof with a simple approximation theorem near the fixed point and
give, in this setting, analytical necessary conditions for existence of the solutions. (The French "conditions
necessaires" is printed in the abstract; the body proves conditions that are sufficient, Theorem 2, and shows by
density that they can be met, Theorem 1. INFERRED: the abstract word is loose.)

The paper is a construction and a density result, not a closed catalogue. It proves: (Theorem 2) if two symmetric
Kepler "double collision" arcs of the same energy, with consecutive time intervals, satisfy three explicit
non-degeneracy conditions, then for every mass parameter mu below an unspecified mu_1 there is a transverse
intersection of two surfaces in the energy manifold, hence a periodic orbit near the pair; (Theorem 1) such good pairs
exist arbitrarily near any prescribed pair of angular momenta in the admissible interval. No numerical value of mu_1,
no distance of closest approach, no uniqueness statement and no numerical example appear.

## 1. Setting (section 1, pp. 213-221)

### 1.1 Model, frames, units (READ, pp. 213-216)

Planar circular restricted problem: "une particule P3 de masse nulle situee dans le plan de deux points P2 et P1, de
masses respectives mu in ]0,1[ et 1 - mu, en mouvement circulaire autour de leur centre de masse C. On posera toujours
nu = 1 - mu." (p. 214, TRANSLATED: a massless particle P3 in the plane of two points P2 and P1 of masses mu and 1 - mu,
in circular motion about their centre of mass C; nu = 1 - mu.) So P2 is the small secondary (mass mu), P1 the primary
(mass 1 - mu). The units are implicit in the equations: separation of P1 and P2 equal to 1, angular velocity 1 (the
paper says in 2.2.1 "un point P2 en orbite circulaire de rayon 1 autour de P1 avec une vitesse angulaire unite", p. 224),
so total mass 1 with G = 1 (INFERRED from equation (1)).

Fixed frame (1.1.1, eq. (1), p. 214): d^2X/dt^2 = -nu (X + mu e^{it})/|X + mu e^{it}|^3 - mu (X - nu e^{it})/|X - nu e^{it}|^3,
X in C the position of P3. Rotating frame (1.1.2): X = xi e^{it}, Y = eta e^{it} (eq. (2)), giving (eq. (3))
dxi/dt = eta - i xi, deta/dt = -i eta - nu (xi + mu)/|xi + mu|^3 - mu (xi - nu)/|xi - nu|^3. Hamiltonian (eq. (4), p. 215):
R_mu(xi, eta) = |eta|^2 + i (eta conj(xi) - conj(eta) xi) - 2 nu/|xi + mu| - 2 mu/|xi - nu| + 2 nu + nu^2, with the symplectic
form Omega = d eta ^ d conj(xi) + d conj(eta) ^ d xi (eq. (5)). The constant is "choisie pour simplifier l'etude des
regularisations". Hill regions: Figure 2 (p. 215), critical values h_0 < h_1 < h_2 < h_3 = 0; the paper always takes
the energy h above h_3 = 0 (and in fact h > 0), so the particle can pass between the neighbourhoods of P1 and P2.
Quote (p. 215): "La nature meme des solutions de seconde espece conduit a ne considerer que des energies superieures a
la premiere valeur critique h_0, rendant ainsi possibles les transitions de la particule P3 entre P1 et P2. Dans le but
de simplifier l'etude, on supposera meme toujours l'energie h superieure a la derniere valeur critique h_3 = 0."
(TRANSLATED: second species solutions require energies above the first critical value so that P3 can move between P1
and P2; to simplify, the energy is always taken above the last critical value, 0.)

Frame centred on P2 (1.1.4, eq. (6)): M = xi - nu, N = eta - i nu; Hamiltonian (eq. (7), p. 216):
H_mu(M, N) = |N|^2 + i (conj(M) N - M conj(N)) - 2 mu/|M| - nu (2/|1 + M| - 2 + (M + conj(M))), symplectic form
dM ^ d conj(N) + d conj(M) ^ dN (eq. (8)). INFERRED from eqs. (3) and (6): at M = 0 (P3 at P2), N = eta - i nu =
d xi/dt + i xi - i nu = d xi/dt, so N is the velocity of P3 relative to P2 in the rotating frame. The paper calls N_s and
N_u "les vitesses" (p. 229-230) with |N|^2 = h at M = 0 (p. 230, eq. after (23)).

### 1.2 The limit problem mu = 0 (p. 216-219, READ)

At mu = 0, H_0 (eq. (9)) is a Kepler problem about P1 in the rotating frame, and (eq. (10), (11), (11')):
H_0 = |eta|^2 + i (eta conj(xi) - conj(eta) xi) - 2/|xi| + 3 = F_0 - 2 sigma, with sigma = (i/2)(eta conj(xi) - conj(eta) xi)
the angular momentum and F_0 = |eta|^2 - 2/|xi| + 3 the fixed-frame energy (with the offset 3). Solutions are called
elliptic if F_0 < 3, parabolic if F_0 = 3, hyperbolic if F_0 > 3; direct if sigma > 0, retrograde otherwise (p. 217).
For h > 0 the energy manifold H = H_0^{-1}(h) splits into H_+ (hyperbolic, cylinders), H_0 (parabolic), H_- (elliptic,
invariant 2-tori) (1.2.2). The velocity circle C_0 over P2 plays the central role: "les orbites donnant naissance aux
solutions de seconde espece le rencontrent deux fois" (p. 217, TRANSLATED: the orbits giving rise to second species
solutions meet it twice). The elliptic singular domain D_e is the union of tori of elliptic orbits meeting C_0. Section
1.2.4 and Figure 3 (pp. 217-219) give, as a function of h, the limiting angular momenta sigma_-(h) = 1 - sqrt(h) and
sigma_+(h) = 1 + sqrt(h) for 0 < h < 3 - 2 sqrt(2), and sigma_+(h) = (3 - h)/2 for larger h up to 3 + 2 sqrt(2); D_e is
empty for h >= 3 + 2 sqrt(2) (p. 218). Hence section 2 takes h in ]0, 3 + 2 sqrt(2)[ (p. 222). READ (the thresholds are
as printed; the sentence at h = 3 - 2 sqrt(2) to 1 has the sigma_+ value (3 - h)/2 for the cylinder H_0).

### 1.3 Levi-Civita regularisation (pp. 220-221, READ)

At fixed energy h > 0, rho(z, w) = (z^2/h, sqrt(h) w / conj(z)) (eq. (12)), with rho*(Omega) = (2/sqrt(h)) Omega, and the
regularised Hamiltonian L_mu = (|z|^2/h)(H_mu - h) o rho (eq. (13)), which extends analytically to C^2:
L_mu(z, w) = |w|^2 - |z|^2 - 2 mu + i (|z|^2/h^{3/2})(conj(z) w - z conj(w)) - (nu/h^3) f(z) (eq. (14)),
f(z) = -|z|^6 + (3/4)|z|^2 (z^2 + conj(z)^2)^2 + O_8(z) (eq. (15)). C_mu is the intersection of L_mu = 0 with z = 0, the
circle corresponding to collision. rho restricted to the complement is a two-sheeted cover of H_mu, and the vector field
is rescaled by the factor xi(z, w) = |z|^2/h^{3/2}.

Local structure (1.3.2): the origin is a hyperbolic fixed point of X_{L_mu} with linear part dz/dt = w, dw/dt = z
(eq. (16)), stable and unstable subspaces w = -z and w = z, signature 2-2 (p. 213). Zero level (eq. (17)):
|w|^2 - |z|^2 = 2 mu + O_6(|z|, |w|). For mu > 0 the zero level meets a small band |z|^2 <= alpha in a full solid torus
with core the circle C_mu; as mu tends to 0 the circle collapses to the origin, and at mu = 0 the level set near the
origin is a cone, "n'est donc plus une sous-variete" (p. 221). Quote (p. 221): "Le cas limite mu = 0 est different. Le
cercle C_0 des vitesses au-dessus du point P2, initialement transverse au champ, est tout entier envoye par rho sur le
point fixe O du systeme de Levi-Civita. [...] Le champ initial est donc singularise lorsque mu = 0." TRANSLATED: the limit
case is different; the circle of velocities over P2, initially transverse to the flow, is sent entirely to the fixed point
O; the original field is therefore singularised at mu = 0. This is the "singular perturbation" the title of the method
refers to. The paper notes that Conley's thesis studies the case h < 0 (particle trapped near P2), which gives an elliptic
fixed point, rho(z, w) = (z^2/h, sqrt(-h) w/conj(z)) (p. 221); that case is not used.

## 2. Definition of second species solutions and the density theorem (section 2, pp. 222-227)

### 2.1 Generating orbits, definition (READ, pp. 222-223)

Fix h in ]0, 3 + 2 sqrt(2)[ and two angular momenta sigma_u, sigma_s in ]sigma_-(h), sigma_+(h)[. Let phi_u, phi_s be
two solutions of the singular limit problem P_0 (the Kepler problem in the rotating frame with the velocity circle over
P2 removed), of angular momenta sigma_u, sigma_s, each asymptotic in both time directions to the singularity circle:
"On suppose de plus que leurs intervalles de definition sont consecutifs, i.e. de la forme Dom(phi_s) = ]-t_c, t_c[ et
Dom(phi_u) = ]t_c, t_c + l[. Dans le plan de configuration, leurs trajectoires sont donc deux arcs d'ellipses (vues en
repere tournant), d'origines et d'extremites au point P2." (p. 222, TRANSLATED: their intervals of definition are
consecutive; in the configuration plane their trajectories are two elliptic arcs, in the rotating frame, starting and
ending at P2.) So each generating arc is a full Keplerian revolution about P1 that leaves P2 and returns to P2 (a double
collision orbit, a loop of the rotating-frame ellipse). Figure 4 (p. 222) shows two doubly singular trajectories of the
same rotating-frame ellipse; Figure 5 (p. 223) shows the arcs phi_s, phi_u and a second species solution (qualitative).

Definition (p. 223): "On note phi la fonction reunion de phi_s et phi_u [...]. On appelle classiquement famille de
solutions de seconde espece associee a phi_u et phi_s une famille (psi_mu), definie pour mu > 0 et voisin de 0, de
fonctions periodiques solutions pour chaque mu du probleme P_mu, qui converge vers phi sur Dom(phi_s) U Dom(phi_u)
lorsque mu -> 0. S'il existe une telle famille, les solutions limites phi_u et phi_s sont dites generatrices."
TRANSLATED: a family of second species solutions associated with phi_u and phi_s is a family (psi_mu), defined for
mu > 0 near 0, of functions that are periodic solutions of P_mu for each mu and converge to phi on Dom(phi_s) U
Dom(phi_u) as mu -> 0; if such a family exists the limit solutions are called generating. The authors then widen the
definition to allow members that meet the collision circle (p. 223): "il sera commode de generaliser un peu la
definition precedente des solutions de seconde espece en prenant en compte les (psi_mu) possedant eventuellement des
singularites de collision P2 P3" (TRANSLATED: it is convenient to allow members with collision singularities).

Consequence of the definition (INFERRED from the periodicity lemma, p. 228): the period is 2 tau with tau = t_c + l/2
and the limit orbit has exactly two collision passages per period, at times t_c and t_c + l, equivalent modulo the period
to -t_c and t_c. So the generating orbit is two Keplerian loops, each beginning and ending at P2, glued at P2 with a
change of velocity direction. "Two collisions per period" in Bolotin and MacKay's description of this paper corresponds to
this reading.

### 2.2 Theorem 1 (density), p. 223-224 (READ)

"THEOREME 1. - Soient sigma_1 < sigma_2 deux moments cinetiques dans l'intervalle ]sigma_-, sigma_+[. Pour tout eta > 0,
il existe des solutions generatrices phi_u et phi_s du probleme P_0, de moments cinetiques respectifs m_1 et m_2
verifiant |m_1 - sigma_1| < eta et |m_2 - sigma_2| < eta."
TRANSLATED: let sigma_1 < sigma_2 be two angular momenta in the interval ]sigma_-, sigma_+[; for every eta > 0 there exist
generating solutions phi_u and phi_s of P_0 with angular momenta m_1, m_2 within eta of sigma_1, sigma_2. (The printed
text says "respectivement m_1 et m_2"; which of phi_u, phi_s takes which is as printed. Energy h is fixed throughout,
section 2 opening.) This is a density theorem: it does not assert that a given pair generates, only that generating pairs
are dense in the (sigma_s, sigma_u) square.

The authors also state the intermediate corollary (p. 227): "Pour une energie h in ]0, 3 + 2 sqrt(2)[, le point fixe du
probleme Q_0 possede une infinite d'orbites homoclines symetriques." (TRANSLATED: for such an energy the fixed point of
the regularised problem has infinitely many symmetric homoclinic orbits.) Proof route, section 2.2.1-2.2.2 (READ, pp.
224-227), after Henon [14] and Bruno [13], [15]: work in the fixed frame; a Keplerian ellipse (a, e) about P1 that meets
the circular orbit of P2 (radius 1, angular velocity 1) twice. Restrict to the symmetric case, time origin midway
between the two collisions (collisions at -t_0, +t_0), P3 starting at pericentre direct, excentric anomaly E from t = 0:
X1 = a (e - cos E), X2 = -a sqrt(1 - e^2) sin E (eq. (18)), t = a^{3/2} (E - e sin E) (eq. (19)). Collision at t_0 gives
cos t_0 = a (e - cos E_0), sin t_0 = -a sqrt(1 - e^2) sin E_0, t_0 = a^{3/2} (E_0 - e sin E_0) (eq. (20)). Eliminating,
cos E_0 = (a - 1)/(a e), cos t_0 = (1 - a(1 - e^2))/e (eq. (21)), and the third equation becomes
2 pi (k a^{3/2} - m) + a^{3/2} [ arccos((a-1)/(a e)) - sqrt(e^2 - ((a-1)/a)^2) ] + arccos((1 - a(1 - e^2))/e) = 0
with k, m integers (p. 226). Writing the bracket expression as G(a, e), a torus (a, e) contains a double collision orbit
iff G(a, e) = 2 pi (k a^{3/2} - m). The case where the orbital period T of P3 is a rational multiple of 2 pi (a^{3/2}
rational) gives infinitely many collisions at one fixed-frame point and is excluded: "On se limite donc au cas ou le
rapport T/2 pi est irrationnel" (p. 224). Density (p. 226-227): at fixed energy h = -1/a + 2 sqrt(a (1 - e^2)) (the
paper's offset relation, as printed) the set {k a^{3/2} - m} is dense in R, and the map psi: (a', e') -> (h(a', e'),
G(a', e') - 2 pi k a'^{3/2}) is a local diffeomorphism near (a, e) for large |k|, so its preimage of (h, 2 pi m) gives
nearby double collision tori of the same energy. The proof of the density corollary is a few lines; the invertibility of
psi for large |k| is asserted ("On verifie alors facilement", p. 227), not computed.

## 3. Method of construction (section 3, pp. 227-233)

### 3.1 Symmetry and the periodicity lemma (READ, pp. 227-228)

The involution S(M, N) = (conj(M), -conj(N)) leaves H_mu invariant for every mu and S*(Omega) = -Omega, so
phi~(t) = S(phi(-t)) is a solution when phi is. A solution is symmetric if its orbit is S-invariant. The set I_mu =
{(M, N) in H_mu : Im(M) = Re(N) = 0} is a one-dimensional S-invariant submanifold with two components: P3 on the axis
through P1 and P2 with velocity perpendicular to it ("a croisement orthogonal", p. 228). The lemma, attributed to
Poincare: "LEMME. - Symetrie et periodicite. Soit phi une solution de P_mu pour laquelle il existe deux instants t_1 et
t_2 dans Dom(phi) tels que phi(t_1) in I_mu et phi(t_2) in I_mu, avec t_1 != t_2. Alors phi est periodique, de periode
T = 2(t_2 - t_1), et symetrique, plus precisement S(phi(T - t)) = phi(t), pour tout t." (p. 228, TRANSLATED: if a solution
meets I_mu at two distinct times, it is periodic of period T = 2(t_2 - t_1) and symmetric.) This is the one place where
the symmetry is used in the existence proof (and it is also what makes the generating solutions symmetric).

### 3.2 The gluing scheme (READ, pp. 228-229)

Start: phi_s, phi_u of energy h with Dom as in 2.1, joined by p o phi_s(t_c) = p o phi_u(t_c) = P2, symmetric with
phi_s(0) in I_0 and phi_u(tau) in I_0, tau = t_c + l/2; A_s = phi_s(0), A_u = phi_u(tau), n_s = p o phi_s(0), n_u = p o
phi_u(tau) the orthogonal axis crossings (p the projection to configuration space). For small mu > 0 let A_s(mu), A_u(mu)
in I_mu be the points with p(A_s(mu)) = n_s, p(A_u(mu)) = n_u. Choose small intervals I_s(mu), I_u(mu) in I_mu around
them and define the two surfaces in the three-dimensional energy manifold H_mu:
V_s(mu) = Phi_mu([0, t_c + eps[, I_s(mu)) and V_u(mu) = Phi_mu(]-tau - eps, 0], I_u(mu)) (p. 229),
obtained by flowing the arcs forward and backward. "Le lemme de periodicite precedent entraine qu'une condition
suffisante pour que les solutions phi_s et phi_u soient generatrices est que les varietes V_s(mu) et V_u(mu) aient une
intersection non vide pour tout mu > 0, l'orbite associee a un point de cette intersection etant periodique." (p. 229;
TRANSLATED: a sufficient condition for the solutions to be generating is that V_s(mu) and V_u(mu) intersect for every
mu > 0, the orbit through an intersection point being periodic.) INFERRED: a point of the intersection lies on an orbit
that left I_mu along the s-side and arrives on I_mu along the u-side, so it crosses I_mu at two distinct times and the
lemma applies; the intersection is then the orbit segment, a curve. The search is restricted to transverse
intersections: "La difficulte, traduisant la nature singuliere de la perturbation provient du fait que cette intersection
disparait pour mu = 0." (p. 229; TRANSLATED: the difficulty, reflecting the singular nature of the perturbation, is that
the intersection disappears at mu = 0.)

Outline of the argument (p. 229, READ): work near the origin of the regularised systems Q_mu; build isolating blocks
(Conley and Easton [16]) B_eps, balls of radius proportional to eps about 0 in a norm adapted to the problem; lift V_s,
V_u by the Levi-Civita map to W_s(mu), W_u(mu); at mu = 0 find the (analytic) intersections of W_s(0), W_u(0) with the
entry set d+B_eps and exit set d-B_eps; these depend analytically on mu, giving first-order approximations for mu != 0;
for small eps the transition maps f_{mu eps} of the blocks are C^1-approximated, uniformly in mu, by their linearisation
at the origin, which does not depend on mu; this gives a transversality condition for [f_{mu eps}(W_s(mu) meet d+B_eps)]
and (W_u(mu) meet d-B_eps) and hence for V_s, V_u.

### 3.3 Local analysis at the origin (READ, pp. 229-233)

Near the collision time, phi_s(t) = (N_s (t - t_c) + O_2, N_s + O_1) and likewise for phi_u (eq. (23)); |N_s|^2 = |N_u|^2 = h.
Lifts gamma_s = (z_s, w_s), gamma_u with rho o gamma = phi, with z_s(t) = v_s sqrt(t_c - t) + O_{3/2}, w_s = -v_s sqrt(t_c - t)
+ ..., z_u = v_u sqrt(t - t_c) + ..., w_u = v_u sqrt(t - t_c) + ... (eq. (24), (24')), v_s = i sqrt(h) r_s, v_u = sqrt(h) r_u,
with r_s, r_u square roots of N_s, N_u. Condition (25) on a complex basis (e_1, e_2, f_1, f_2), choice (28)
e_1 = v_s, f_1 = v_u, e_2 = i v_u, f_2 = i v_s, and the requirement that v_s, v_u be non-perpendicular, "ce qui equivaut a
N_s et N_u non colineaires" (p. 231). In the new coordinates the Hamiltonian is
l_mu(s~, u~) = 2 (s~_1 u~_1 + s~_2 u~_2) + mu/alpha + (sum s_i u_j) O_2 + O_6 (eq. (27)), alpha = (e_1|f_1), and after
straightening the invariant manifolds, l_mu(s, u) = 2 (s_1 u_1 + s_2 u_2) + mu/alpha + ... (eq. (27')), with the arcs phi_s,
phi_u carried by the axes O s_1, O u_1. The isolating block is B_eps = B(eps h^{3/4}) for the sup norm
||(s, u)|| = sup(|s|, |u|), with entry set C_s(eps) = {|s| = eps h^{3/4}, |u| <= eps h^{3/4}} and exit set C_u(eps)
(p. 232). The linearised flow is ds/dt = -s, du/dt = u (independent of mu) with transition map f_eps(s, u) =
((|u|/(eps h^{3/4})) s, (eps h^{3/4}/|u|) u) (eq. (32)); in the entry and exit angle coordinates (x, y, theta) it is
x_u = x_s, y_u = -y_s, e^{i theta_u} = u/|u| (eq. (33), (34)). The approximation theorem used is quoted (p. 233): "l'application
de transition linearisee f_eps approche l'application de transition exacte f_{eps mu} en topologie C^1 a l'ordre eps^2"
(TRANSLATED: the linearised transition map approximates the exact one in the C^1 topology to order eps^2), with a
reference to Wiggins [17], section 3.2, propositions 3.2.5 and 3.2.6. The authors say this replaces Henrard's use of a
Hartman-type conjugacy, which does not preserve energy; see section 5.

## 4. Existence theorem and its conditions (section 4, pp. 233-249)

### 4.1 The three lemmas (READ, pp. 233-238)

The semi-major axis a parametrises I_mu near phi_s(0) and phi_u(0): analytic parametrisations i_s^mu(a), i_u^mu(a) with
i_s^0(a_s) = phi_s(0), i_u^0(a_u) = phi_u(0) (p. 233). Lemma a (p. 234): for eps and mu small, the implicit relation
Phi~_mu(t, i_s(a)) in C_s(eps) defines a time tau_s^{mu eps}(a), the map q_s^{mu eps}(a) = Phi~_mu(tau_s, i_s(a)) is an
analytic diffeomorphism of an interval I_{delta_s} = ]a_s - delta_s, a_s + delta_s[ onto its image in C_s(eps), and the
images are the arcs W_s(mu) meet d B_eps (same for the u side). Lemma b (p. 234, formulas (35), (36)): first-order Taylor
expansion of h_s(t, a) = Phi_0(t, i_s^0(a)) at (t_c, a_s):
M(t, a) = N_s (t - t_c) + Delta M_s (a - a_s) + O_2, N(t, a) = N_s + Delta N_s (a - a_s) + O_2, with
Delta M_s = 1/a_s + (N_s + i) [ ((1 - a_s sqrt(sigma_s)) / (2 a_s^2 e_s)) sin t_c - 3 t_c/(2 a_s) ]
+ (a_s sigma_s^2 - sigma_s)/(2 a_s^3 e_s) (a_s - i 2 e_s sin t_c/(1 - e_s^2)) e^{-i t_c}, and a similar expression for
Delta N_s, sigma_s = sqrt(a_s (1 - e_s^2)) (as printed; the same sigma_s symbol is used for the angular momentum in
section 2, which is the usual a_s (1 - e_s^2) under the root; I did not re-derive these two expressions and transcribe
them as printed, and the exact grouping of terms in Delta M_s is hard to read in the scan). The proof (p. 235) uses the
fixed-frame excentric-anomaly parametrisation X = a (e + cos E) + i a sqrt(1 - e^2) sin E, t = a^{3/2}(E + e sin E)
(note the sign convention here differs from eq. (18)-(19); I transcribe as printed). Lemma c (pp. 235-238, formulas
(37)-(39)): in the block coordinates, q_s^0(a_s) = (0, 0, 0), dx_s/da = O(1), dy_s/da = -Im(Delta M_s conj(N_s))/(2 h^{5/2} eps^2)
+ O(1), dtheta_s/da = (1/eps) ds_2/da + O(eps), and the analogous u-side formulas.

### 4.2 Theorem 2 (existence of generating solutions), pp. 238-239 (READ)

Definition (p. 238): V_s(0) (resp. V_u(0)) is "p-reguliere en P2" when the restriction of the projection p to V_s(0) is
a local diffeomorphism near the lift of P2. By Lemma b, a necessary and sufficient condition is Im(Delta M_s conj(N_s)) != 0
(resp. Im(Delta M_u conj(N_u)) != 0).

"THEOREME 2. - On suppose que les varietes V_s(0) et V_u(0) sont p-regulieres en P2, et que les vitesses N_s et N_u ne
sont pas colineaires. Il existe un reel mu_1 > 0 tel que, pour tout mu < mu_1, les varietes V_s(mu) et V_u(mu) aient une
intersection transverse dans la variete H_mu. Les solutions phi_s et phi_u sont donc alors generatrices."
TRANSLATED: assume V_s(0), V_u(0) are p-regular at P2 and the velocities N_s, N_u are not collinear. There is a real
mu_1 > 0 such that for every mu < mu_1 the surfaces V_s(mu), V_u(mu) intersect transversally in H_mu; the solutions
phi_s, phi_u are therefore generating.

Hypotheses, stated in one place: (i) phi_s, phi_u are symmetric double collision Kepler arcs of the same energy h in
]0, 3 + 2 sqrt(2)[ with consecutive time intervals (section 2.1, 3.2); (ii) the orbital period ratio is irrational
(section 2.2.1); (iii) Im(Delta M_s conj(N_s)) != 0 and Im(Delta M_u conj(N_u)) != 0; (iv) N_s, N_u not collinear. The
proof also uses a choice of square roots with alpha = (v_s|v_u) < 0 (p. 246), which the authors say is possible whenever
N_s, N_u are not collinear. No numerical value or estimate of mu_1 is given anywhere; it depends on the generating pair
and is "sufficiently small" in the printed text only through the choice of eps with mu = eps^{11/2}. Uniqueness is not
stated; the proof produces a transverse intersection point with a_s, a_u ranging over intervals "dont les longueurs
respectives sont de l'ordre de eps^{17/3} et eps^{5/2}" (p. 248, as printed; these exponents appear inconsistent with
the delta_s = eps^{7/3} chosen on p. 240 and I could not reconcile them from the printed text; TRANSLATED: the intervals
of semi-major axis have lengths of order eps^{17/3} and eps^{5/2}). The convergence of the periodic solutions to the
generating pair is stated to follow immediately: W_s(mu) meet W_u(mu) meets d B_eps in points converging to the
intersections of phi_s, phi_u with d B_eps (p. 248).

### 4.3 Proof of Theorem 2 in six steps (READ, pp. 239-248)

Step a (pp. 239-240): domains of the arcs. Gronwall type estimate with a Lipschitz constant C_0 of the Kepler field and
analyticity of i_s^0 give delta_s = eps^{7/3}; then in Levi-Civita coordinates, p(C_eps) image: (M, N)(Phi_0(tau_s^0(a),
i_s^0(a))) = (-eps^2 N_s, N_s) + O(eps^{7/3}), so (s, u)(q_s^0(a)) = (eps, 0) + O(eps^{4/3}). Same for u with
delta_u = eps^{7/3}. Step b (pp. 240-241): derivative estimates for the unperturbed arc, d(s,u)/da (q_s^0(a)) = d(s,u)/da
(q_s^0(a_s)) + O(eps^{4/3}), via a Gronwall bound on Z_0. Step c (pp. 241-242): gap between the perturbed (mu != 0)
and unperturbed arcs: outside a disc of radius eps about P2 the gap is O(mu/eps^2) over a time of order 1; between the
disc and C_s(eps) (travel time of order eps) it is O(mu/eps^4) * eps = O(mu/eps^3); analytic closeness of i_s^mu to i_s^0
gives O(mu). Step d (p. 242): the choice that links the two small parameters: "Nous allons fixer ici mu = eps^{11/2}, ce
qui garantira la majoration suivante || Phi_mu - Phi_0 || = O(eps^{5/2})" (TRANSLATED: we fix mu = eps^{11/2}, which
guarantees the bound O(eps^{5/2}) on the gap). The perturbed arc then crosses rho(C_s(eps)) with
|M| = eps^2 |N_s| + O(eps^{7/3}) + O(eps^{5/2}) > 0, and (s, u)(q_s^mu(a)) = (s, u)(q_s^0(a)) + O(eps^{3/2}) (eq. (40)),
with |tau_s^mu(a) - tau_s^0(a)| = O(eps^{5/2}). Step e (pp. 243-245): the same for the a-derivatives (Gronwall again),
giving d(s,u)/da (q_s^mu(a)) = d(s,u)/da (q_s^0(a)) + O(eps^{-1/2}) (eq. (40')); equivalent estimates near i_u.
Step f (pp. 245-248): straighten the level set L_mu meet C_s(eps) by the analytic maps T_mu^s, T_mu^u (eps-close to the
identity) so that x~_s = -mu/(2 alpha eps^2 h^{3/2}); using the hypothesis Im(Delta M_s conj(N_s)) != 0, u_2(q_s^mu(a)) varies
over an interval of length order eps^{4/3} as a varies, hence a point m_mu^(s) exists with u_2 = 0, with (s_1, s_2, u_1) = (eps,
0, 0) + O(eps^{4/3}) (eq. (41)) and u_1(m_mu^(s)) = -mu/(2 alpha eps) + O(mu/eps^{2/3}) = -eps^{9/2}/(2 alpha) + O(eps^{29/6})
(eq. (42)); with alpha < 0 this u_1 is positive, with argument zero. Image under the linearised transition map: y'_u = O(eps^{2/3}),
theta'_u = 0 and, via the C^1 approximation, theta'_u = O(eps^2) at the exact image; derivative with respect to a of order
1/eps^{11/2}-type quantities as printed (dy'_u/da = O(1/eps^2), dtheta'_u/da = O(1/eps^{11/2}) in eq. after (42)). The u-side arc is
treated symmetrically (eqs. (43), (43')), and the tangents are shown to be nearly along the theta_u axis (the image of
W_s(mu) meet C_s(eps)) and the y_u axis (W_u(mu) meet C_u(eps)), which gives transversality. I did not re-derive these
estimates; they are transcribed as printed, and several of the later exponents (for example the exponent of the a-derivative
of theta') are not individually legible with full certainty in the scan beyond the figures given above.

### 4.4 Proof of Theorem 1 (p. 248, READ)

"On fixe toujours les solutions phi_s et phi_u de moments cinetique sigma_s et sigma_u. Remarquons d'abord que la fonction
moment cinetique, restreinte a un voisinage des conditions initiales (correspondant a phi_s et phi_u) dans I_0, est
strictement monotone. Le theoreme de densite 2.2.2 entraine donc la densite des conditions initiales sur I_0
correspondant a des solutions a double collision. L'idee est maintenant tres simple, il suffit de verifier que les fonctions
Im(Delta M_s conj(N_s)) et Im(Delta M_u conj(N_u)), restreintes a la sous-variete I_0, sont analytiques et non constantes,
ce qui est clair au moins dans un voisinage des conditions initiales considerees. Le principe des zeros isoles entraine
alors directement l'existence de conditions initiales arbitrairement proches de celles de depart, donnant lieu a des
solutions generatrices, et le theoreme est demontre."
TRANSLATED: angular momentum is strictly monotone on I_0 near the initial data, so double collision data are dense; it
suffices to check that Im(Delta M_s conj(N_s)) and Im(Delta M_u conj(N_u)) are analytic and non-constant on I_0, which is
"clear at least in a neighbourhood of the initial data"; the isolated zeros principle then gives nearby data where both
are nonzero. INFERRED: the non-constancy is asserted, not computed, and non-collinearity of N_s, N_u is not discussed in
this proof (it is presumably an open condition satisfied away from collinear pairs). Note also that the proof of
Theorem 1 does not discuss the choice of root with alpha < 0 or the sign conditions on the double-collision pairing of
k, m for two distinct arcs with consecutive intervals. These are small gaps in the printed argument, not errors I can
demonstrate.

### 4.5 Closing remark (pp. 248-249, READ)

"Remarque. - La transversalite observee dans la demonstration precedente, associee a la forme de l'application f_{eps mu},
ne contredit pas l'integrabilite eventuelle du probleme pour mu != 0. On aurait en effet une situation analogue pour le
probleme des deux centres fixes. On peut simplement affirmer qu'il n'existe pas de famille g_mu d'integrales premieres du
probleme, dependant differentiablement de mu, qui converge pour mu -> 0 vers le moment cinetique. L'etude detaillee et
geometrique de ce type de probleme donnera lieu a un prochain travail." TRANSLATED: the transversality does not contradict
possible integrability for mu != 0 (the same would happen for two fixed centres); one can only say that there is no
differentiable-in-mu family of first integrals converging to the angular momentum as mu -> 0; a detailed geometric study
will be the subject of a later work.

## 5. What the paper says about earlier work and what remains open (READ)

The paper has no separate open-problems section. Everything is in the introduction (pp. 212-213) and in the remarks above.

* Poincare (p. 212): "Dans le premier tome des Methodes Nouvelles de la Mecanique Celeste, H. Poincare met en evidence
l'importance des solutions periodiques [...] Depuis Poincare, les solutions periodiques ainsi construites sont dites de
premiere espece. Dans le troisieme tome [1], Poincare s'apercoit que d'autres solutions periodiques du probleme general
s'obtiennent en perturbant des solutions avec singularites du probleme a deux masses nulles - les singularites
correspondant aux collisions de ces deux particules [...] Il baptise solutions de seconde espece ces nouvelles solutions
periodiques. La demonstration donnee par Poincare de l'existence des solutions de seconde espece est tres courte, et
parait assez difficile a justifier completement. En particulier, la geometrie du probleme n'est pas eclaircie."
TRANSLATED: Poincare's proof is very short and seems quite hard to justify fully; in particular the geometry of the
problem is not clarified. (Reference [1] cites Methodes Nouvelles, tome 3, sections 385-391.) Note that Poincare's
original problem has two small masses; the paper deliberately treats "un probleme plus simple, directement lie au
probleme initial de Poincare, et etudie aussi par Hill, celui de la construction des solutions de seconde espece dans le
probleme restreint" (p. 212; TRANSLATED: a simpler problem, directly tied to Poincare's and also studied by Hill, the
restricted one).
* Henrard [8] (1980): the idea of the construction is his (abstract, p. 213). Critique, p. 213: "Il faut noter que la
demonstration de Henrard est basee sur l'utilisation d'un theoreme de conjugaison (a la Hartman) au voisinage du point
d'equilibre, difficile a utiliser correctement dans un cadre symplectique. En particulier, cette conjugaison ne conserve
pas l'energie et le procede de raccordement de Henrard porte sur des segments de solutions d'energies differentes, ce qui
est incorrect. Pour eviter ce probleme, nous utilisons ici un simple theoreme d'approximation C^1 autour du point fixe, qui
nous permet de tenir explicitement compte de la conservation de l'energie." TRANSLATED: Henrard's proof uses a
Hartman-type conjugacy near the equilibrium, hard to use correctly symplectically; the conjugacy does not preserve energy
and Henrard's gluing involves solution segments of different energies, which is incorrect; the authors instead use a C^1
approximation theorem that keeps the energy conservation explicit. (A direct statement that an earlier published proof is
flawed; for the project this matters only as a caution that the energy bookkeeping in the matching is decisive.)
* Numerical analytic precedents (p. 213): "Nous donnons de plus les expressions analytiques explicites des conditions
d'existences des solutions de seconde espece en termes de conditions initiales pour les orbites homoclines limites,
(alors que les etudes precedentes de ce probleme etaient seulement numeriques [9])." TRANSLATED: they give explicit
analytic existence conditions in terms of the initial data of the limit homoclinic orbits, whereas previous studies were
only numerical [9] (Henon and Hitzl 1977, critical generating orbits).
* Others cited as "nombreuses etudes analytiques (cf. [2]-[7])" (p. 212): Perko 1976, 1977, 1981 (near-Moon passage at
O(mu), O(mu^nu) with 1/3 < nu < 1 and 0 < nu < 1), Guillaume 1973 and 1975, Gomez and Olle 1991 (circular and elliptic
restricted problem). No individual comment is made on any of them. Henon [14] (orbits meeting the planet twice, 1968) and
Bruno [13], [15] (1990 book; "On periodic flybys to the Moon", 1981) are credited only for the double-collision
generating ellipses used in section 2.2.1 ("d'apres une etude de Henon [14] et Bruno [15], [13]", p. 224). Conley [12]
(1963) and Chenciner [11] are cited for the h < 0 regularisation and for the regularisation reference.
* Open or deferred: only the closing remark above ("L'etude detaillee et geometrique [...] prochain travail"). The paper
says nothing about more than two collisions per period, the elliptic restricted problem, the spatial problem, the full
three-body problem with two small masses, two secondaries, or time-dependent problems. It states only the planar
circular restricted case with one small secondary. (READ: searched the whole text; no such sentence.)

## 6. What this does and does not cover for #890

### (a) What is established for one small secondary

READ. For the planar circular restricted three-body problem with one secondary of mass mu, at a fixed energy h in
]0, 3 + 2 sqrt(2)[ (paper's constant), the following holds. If two symmetric Kepler loops about the primary, each starting
and ending at the secondary's position, with consecutive intervals, an irrational period ratio, the same energy, and
satisfying Im(Delta M conj(N)) != 0 on each loop and non-collinear collision velocities N_s, N_u, are given, then for all
mu below an unspecified mu_1 the matching surfaces intersect transversally and there is a symmetric periodic solution
near the pair, with two collision-type passages per period and period about 2 t_c + l (INFERRED from the lemma). Generating
pairs satisfying the conditions are dense in the pair of angular momenta (Theorem 1, with a proof that asserts the key
non-constancy). Not printed: any numerical mu_1, a minimum distance (the only printed length is the block entry distance
eps^2 |N_s| = mu^{4/11} sqrt(h) with mu = eps^{11/2}, INFERRED from p. 242), uniqueness, stability, the number of
families, or any numerical example.

### (b) Same idea as the project's continuation

INFERRED, with the paper's own limit object being the same as the project's zero-mass object. In Marco and Niederman the
mu = 0 limit is two Kepler arcs joined at the secondary's position with a change of velocity direction (|N_s| = |N_u| =
sqrt(h), N_s not collinear with N_u). That is exactly a patched-conic chain whose flyby has been shrunk to a point mass
pivot: a patched-conic flyby of given v-infinity and turn angle has periapsis r_p = (GM/v_inf^2)(1/sin(delta/2) - 1)
(standard hyperbolic formula, not from the paper), which tends to zero with GM at fixed turn angle. So the generating
orbit of this paper is the project's zero-mass seed. The method differs in what it attaches to that seed. Marco and
Niederman (like Bolotin and MacKay) take mu to zero asymptotically, regularise the collision, treat the saddle at the
Levi-Civita fixed point by a linearised transition map and obtain an existence statement for mu < mu_1, with the
perturbed orbit's closest approach of order mu in magnitude (INFERRED: from eq. (27') the saddle level set s.u ~ -mu/alpha
puts |z| of order sqrt(mu), hence |M| = |z|^2/h of order mu/h; this is my estimate, not printed). The project instead
continues in mu from a small value to the physical masses with each flyby seeded by the hyperbola that provides the
patched-conic turn, by multiple shooting on the full model, which is the approach in the Bradley and Russell digest
(continuation of a patched-conic trajectory to full gravity by parametrising masses, flyby altitudes and sphere sizes). So
the three share the same start (a patched or collision chain) and the same direction (continue in mass), but only the
project and Bradley and Russell compute the continuation numerically up to finite physical masses; Marco and Niederman
prove only an asymptotic existence result at small mu and give no scheme to follow it to a finite mass. The
transverse-intersection formulation (two surfaces V_s, V_u in the energy manifold) is in principle a numerical matching
scheme, but the paper does not use it numerically.

### (c) Does anything depend on autonomy and a Jacobi integral

INFERRED, steps named. Yes, in four places, so the proof does not carry over to two moons with different periods without
new work:
1. Fixed energy: the Levi-Civita regularisation is done on a fixed level H_mu = h (eq. (12), (13), time rescaling by
|z|^2/h^{3/2}), which needs an autonomous Hamiltonian. With two moons on circles of different periods there is no
rotating frame in which both are stationary, so there is no Jacobi integral and no fixed energy level.
2. Energy bookkeeping in the matching: the paper's central correction of Henrard is that both arcs have the same energy h.
For two moons the matching would be on a time-dependent flow; the sections C_s(eps), C_u(eps) and the surfaces V_s, V_u
would have to be built in the extended phase space, with the moon phases as extra coordinates.
3. Time-reversal symmetry S: the periodicity lemma (p. 228) uses the reversing involution S of an autonomous system and
orthogonal crossings of the axis through P1 and P2; with two moons of incommensurate periods the symmetry holds only at
instants when both moons are on a common axis, which does not recur, so an orbit symmetric in the paper's sense cannot be
periodic unless the periods are commensurate (the project's orbit is periodic in the full model, so it would need a
commensurate period relation; see the #890 note).
4. Local saddle analysis: this step is local to one secondary and does not itself use the Jacobi integral beyond fixing the
level L_mu = 0 (eq. (17)); I expect it to carry over to each flyby separately (INFERRED), the open issue being the matching
of the two local analyses through a time-dependent Keplerian transfer.
The irrationality condition (ii) and the density argument use only the one-secondary double collision geometry and are
also specific to one moon.

### (d) Wording

Supported by the printed text: "second species solutions of the planar circular restricted three-body problem with one
small secondary were shown to exist near symmetric Kepler double-collision generating pairs for sufficiently small mass
parameter (Marco and Niederman 1995, Theorem 2; such pairs are dense, Theorem 1)"; "the zero-mass limit of the #890 chain
is a pair of Kepler arcs patched at the moon's position, the same generating object used in the one-moon second species
constructions"; "the proof gives no numerical bound on the mass parameter, so it does not establish a periodic orbit at the
Uranian moon masses". Not supported: that the paper proves the existence of a periodic orbit for a given generating orbit
or for all mu (only for mu < mu_1, unspecified, under three explicit conditions that must be checked); that it covers
two secondaries, elliptic or spatial motion, three or more collisions per period, stability, or uniqueness; that its
construction is the same as the project's finite-mass continuation (it is an asymptotic existence argument, not a
continuation to finite mass); and that it says anything on non-autonomous (two-period) problems. The #890 two-moon periodic
orbit should therefore be described as an instance of a patched-conic chain continued to physical masses by shooting,
motivated by but not covered by the one-moon second species theory of Poincare, Henrard-type constructions and Marco and
Niederman; Bolotin and MacKay 2000, not this paper, is the one that states uniform statements for chains of collision arcs
(see their digest).

## 7. References on second species solutions as printed (p. 249)

[1] H. Poincare, Les Methodes Nouvelles de la Mecanique Celeste, 1987, Blanchard, tome 3, Sections 385-391.
[2] L. M. Perko, Second species periodic solutions with an O(mu) near-Moon passage, Cel. Mech., vol. 14, 1976, p. 395-427.
[3] L. M. Perko, Second species solutions with an O(mu^nu), 1/3 < nu < 1, near-Moon passage, Cel. Mech., vol. 16, 1977,
    p. 275-290.
[4] L. M. Perko, Second species solutions with an O(mu^nu), 0 < nu < 1 near-Moon passage, Cel. Mech., vol. 24, 1981,
    p. 155-171.
[5] P. Guillaume, Families of symmetric periodic orbits of the restricted three-body problem, when the perturbing mass is
    small, Astronomy and Astrophysics, vol. 3, 1973, p. 57-76.
[6] P. Guillaume, Linear analysis of one type of second species solutions, Cel. Mech., vol. 11, 1975, p. 449-467.
[7] G. Gomez et M. Olle, Second species solutions in the circular and elliptic restricted three body problem, I et II,
    Celestial Mechanics, vol. 52, 1991, p. 107-166.
[8] J. Henrard, On Poincare second species solutions, Cel. Mech., vol. 21, 1980, p. 83-97.
[9] M. Henon et D. L. Hitzl, Critical generating orbits for second species solutions of the restricted problem, Cel. Mech.,
    vol. 15, 1977, p. 421-452.
[10] T. Levi-Civita, Sur la resolution qualitative du probleme restreint des trois corps, Opere matematishe, vol. 2,
    Bologna, 1956.
[11] A. Chenciner, Le probleme de la lune et les systemes dynamiques, (I), Preprint Universite Paris-VII, 1989.
[12] C. Conley, On some new long periodic solutions of the plane restricted three-body problem, Comm. on Pure and Applied
     Math. XVI, 1963, p. 449-467.
[13] A. D. Bruno, Le probleme restreint des trois corps, 1990, Nauka, Moscou.
[14] M. Henon, Sur les orbites interplanetaires qui rencontrent deux fois la terre, Bulletin Astronomique, Serie 3, 1968,
     p. 377-393.
[15] A. D. Bruno, On periodic flybys to the Moon, Cel. Mech., vol. 24, 1981, p. 255-268.
[16] C. Conley and R. Easton, Isolated invariant sets and isolating blocs, Trans. A.M.S., vol. 158, 1971, p. 35-61.
[17] S. Wiggins, Global Bifurcations and Chaos, 1988, Springer.
[18] M. Chaperon, Geometrie differentielle et singularites des systemes dynamiques, Asterisque, 1986.
Bibliographic spellings are as printed (including "Cel. Mech." and "Opere matematishe").

## 8. Reliability notes

All formulas were read from page images. Items not re-derived and uncertain in transcription: the exact grouping of terms
in Delta M_s and Delta N_s (eq. (36), p. 234); the sign convention in the Lemma b proof (p. 235) which differs from eq.
(18)-(19) (printed as such); the final-step exponents on p. 248 (eps^{17/3}, eps^{5/2}), which do not obviously agree with
delta_s = eps^{7/3} on p. 240; and the derivative orders in the equation after (42). None of these affect the statements of
Theorems 1 and 2, which are legible.
