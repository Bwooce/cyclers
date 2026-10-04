# Digest: Burdet 1967, "Regularization of the two body problem"

C. A. Burdet, "Regularization of the Two Body Problem", Zeitschrift fuer angewandte Mathematik und Physik (ZAMP) 18:434-438 (1967), Kurze Mitteilungen (Brief Reports), DOI 10.1007/BF01601283 (received 4 February 1967; Institute for Applied Mathematics, Eidgenoessische Technische Hochschule, Zurich). Journal pages 434 to 437 are the paper and the references run onto the top of p.438; the rest of p.438 is the start of an unrelated note (H. Brunner, "Stabilization of Optimal Difference Operators"), which I did not use. Filed in the private paper corpus as `burdet-1967-regularization-two-body-problem-zamp-18-434-doi-10.1007-BF01601283.pdf` (5 PDF pages = journal pages 434 to 438; a scan with a poor OCR layer).
Which pages I read: all five were read as 150 dpi page images; the text layer was not used. Every equation (1) to (14) is on pp.434 to 437 and was read on an image. No symbol was illegible.
Evidence tags: READ (journal page) is what the printed page says; COMPUTED is my own hand algebra of 2026-10-04 (each marked, with what it checked); INFERRED is my reading across sources; RECALLED is from my memory of the literature, not from any document held, and is not to be relied on without a source.
Companions: `docs/notes/2026-10-04-digest-heggie-1974-global-regularisation.md`, `docs/notes/2026-06-19-digest-bond-allman-2021-modern-astrodynamics.md` (its Chapter 9 develops this method into the full Sperling-Burdet element system with numerical tables), `docs/notes/2026-10-04-896-published-checks-added.md`.

## 0. What the paper is

A four-page communication announcing a new regularisation of the two-body problem, and of the perturbed two-body problem, which keeps the ORDINARY Cartesian position vector as the dependent variable and changes only the independent variable. The author says (p.434) that a thesis is forthcoming with the theory of the general two-body problem and "various integration methods based on our regularization"; the paper itself gives no numerical example, no table and no number to reproduce. It treats elliptic orbits only (p.435: "To simplify we shall here consider orbits of elliptic type only"; the parabolic and hyperbolic cases are deferred to "the announced global study of Kepler motion"). The impulse is credited to Sperling (ARS J., May 1961, pp.660-661), who obtained regular third-order equations for unperturbed Kepler motion by differentiating the classical equations. The motivation given is the Kustaanheimo-Stiefel (KS) result (J. reine angew. Math. 218:204, 1965) of regular linear independent second order equations in 4 dimensions; Burdet claims "similar properties" with 3 equations in 3 dimensions.

Why it matters to this project: it is the origin of the Sperling-Burdet formulation (elements with variable coefficients, equations regular at r = 0, a Kepler orbit propagated EXACTLY by constant elements), and it is the only regularisation in the held corpus that needs no coordinate map, so it works identically in 2 and 3 dimensions and for a velocity-dependent perturbation. The project's close-pass problems (#899, #924, #928) are perturbed two-body motion about a moon.

## 1. The regularisation, written out (READ pp.434 to 437)

Setting. Cartesian axes centred on the attracting body; position x_i(t), i = 1, 2, 3; r^2 = sum_k x_k^2; gravitational parameter M, set to 1 "without loss of generality". Unperturbed motion (1): d^2x_i/dt^2 + (M/r^3) x_i = 0. Burdet names two defects: (a) singular at r = 0, which prohibits integrating ejection-collision trajectories; (b) nonlinear and coupled, which is "responsible for poor numerical stability".

### 1.1 Two integrals of motion (READ p.435)

Equations (2) and (3):
- 1/a = 2/r - v^2 (the vis-viva energy integral, a the semi-major axis);
- A_i = (v^2 - 1/r) x_i - u dx_i/dt, with u = sum_k x_k dx_k/dt (that is, r dr/dt) and v^2 = sum_k (dx_k/dt)^2.
A_i is "one of the so-called Laplace vectors" (the eccentricity vector with this sign and scale, |A| = e, see (8) and the standard form below). Both are constants of pure Kepler motion; proof is by (10d) with P = 0.

Substituting (2) and (3) into (1) gives (4):
d^2x_i/dt^2 + (u/r^2) dx_i/dt + (1/(r^2 a)) x_i + (1/r^2) A_i = 0.
COMPUTED check: with a and A constant, (u/r^2) v + x (2/r - v^2)/r^2 + ((v^2 - 1/r) x - u v)/r^2 = x/r^3, so (4) is (1). The point is that the nonlinear -x/r^3 has been rewritten as terms that are LINEAR in x and dx/dt with coefficients that are constants of the motion (but still singular as 1/r^2 in (4)).

### 1.2 The independent variable (READ p.435)

Equation (5): dt = sqrt(a) r dE, so t = integral of sqrt(a) r dE (E "is the eccentric anomaly up to the additive constant of integration"). Multiplying (4) by a r^2 gives the boxed pair (6a, 6b):

 d^2x_i/dE^2 + x_i + a A_i = 0, i = 1, 2, 3,
 dt/dE = sqrt(a) r.

COMPUTED check of (6a): d/dE = sqrt(a) r d/dt, so d^2x/dE^2 = a (u v + r^2 d^2x/dt^2), and (4) gives r^2 d^2x/dt^2 = -u v - x/a - A, so d^2x/dE^2 = -x - aA. This confirms (6a) and shows where the energy enters: through the factor a in the time scale and in the forcing a A.

Burdet's statement of the result: three regular, linear, independent second order equations with constant coefficients for the position vector, because a and A are elements "they can be computed from the initial conditions with the help of (2) and (3)". Note what is regular and what is not: x_i(E) is regular in E through r = 0 (the equation has no r in it); the time equation (6b) is the only remaining place where r appears, and r = |x| is smooth in E.

### 1.3 Solution and element map (READ pp.435 to 436)

General solution (7a, 7b): x_i(E) = [alpha_i cos E + beta_i sin E - A_i] a; dx_i/dE = [-alpha_i sin E + beta_i cos E] a, with six arbitrary constants alpha_i, beta_i (the initial data). Setting E = 0 and r = r0 at t = t0: alpha_i = A_i + (1/a) x_i(t0), beta_i = (r0/sqrt(a)) dx_i/dt at t0. COMPUTED check: (7a) at E = 0 gives x0 = a (alpha - A); (7b) at E = 0 gives dx/dE = a beta = sqrt(a) r0 v0; both agree with the printed alpha and beta.

Correspondence of A with alpha, beta (8a to 8c, "without going into details of a transformation to principal axis"):
- A_i = e [alpha_i cos E_p + beta_i sin E_p],
- e^2 = 2 - sum alpha_k^2 - sum beta_k^2,
- E_p = (1/2) arctan[ 2 sum alpha_k beta_k / (sum alpha_k^2 - sum beta_k^2) ], the eccentric anomaly at t = t0.
Standard form, when the start is at pericentre (E_p = 0): the vector alpha points to pericentre and has unit length, beta is perpendicular to alpha in the orbital plane, and with p = a (1 - e^2): sum alpha_k^2 = 1, sum beta_k^2 = 1 - e^2, sum alpha_k beta_k = 0, A_i = e alpha_i. Check, COMPUTED: 1 + (1 - e^2) = 2 - e^2, consistent with the printed e^2 formula.
Unperturbed Kepler's equation follows by integrating (6b) with r = |x| from (7a) at standard form: r = a (1 - e cos E) and t = a^(3/2) (E - e sin E) (COMPUTED from (5) and (7a); not printed).

### 1.4 Perturbed motion (READ pp.436 to 437)

Perturbing force P_i = P_i(x_j, dx_k/dt, t), so velocity dependence is allowed (p.436). Equation (9): d^2x_i/dt^2 + (M/r^3) x_i = P_i. Now a and A_i vary with E, a = a(E), A_i = A_i(E), and (6) becomes the boxed system (10):

 (10a) d^2x_i/dE^2 + x_i + a(E) A_i(E) = a(E) r^2 P_i,
 (10b) dt/dE = sqrt(a(E)) r,
 (10c) d/dE (1/a) = -2 f1,
 (10d) dA_i/dE = 2 f1 x_i - f2 P_i - f3 dx_i/dE,

with f1 = sum_k P_k dx_k/dE ("work done by the perturbing force P_i during the motion", per unit E), f2 = sum_k x_k dx_k/dE, f3 = sum_k P_k x_k. Relations (10c) and (10d) "can be established by a straightforward differentiation of (2) resp. (3), replacing d^2x_i/dt^2 by its value in (9)" (p.437, citing Stiefel's chapter in Dynamics of Rockets and Satellites, ed. Groves, North-Holland).
COMPUTED checks, by hand, 2026-10-04: (10c) follows from d(2/r - v^2)/dt = -2 v.P; (10d) follows from differentiating (3) with dv/dt = P - x/r^3, giving dA/dt = 2 (v.P) x - (x.P) v - u P, and multiplying by dt/dE = sqrt(a) r (using f2 = sqrt(a) r u). Both agree with the printed forms. (10a) was checked in the form d^2x/dE^2 = a (u v + r^2 d^2x/dt^2) above.

Order count (READ p.437): (10) is "a system of 8 simultaneous differential equations for the quantities x_i(E), t(E), a(E), A_i(E)" of total order 11 (3 second-order in x, plus t, a and three A: 6 + 1 + 1 + 3 = 11, COMPUTED). Key statement: "r = 0 is no longer a singularity", and (10a) is linear in x. The perturbation P enters on the right as r^2 P (sqrt of a absorbed); for a velocity-dependent P the velocity dx/dt is dx/dE divided by sqrt(a) r, so a force with a term linear in velocity (a Coriolis term) gives r^2 P ~ r (dx/dE)/sqrt(a), still regular (INFERRED, my reading of the printed dependence of P on dx_k/dt; the paper does not discuss it).

### 1.5 Variation of constants (READ p.437)

Because (10a) is linear in x, Burdet introduces variable natural elements: (11a, 11b) x_i(E) = [alpha_i(E) cos E + beta_i(E) sin E - A_i(E)] a(E) and dx_i/dE = [-alpha_i(E) sin E + beta_i(E) cos E] a(E), with (12) alpha_i(E) = (alpha_i)_{t=t0} + Delta alpha_i(E), beta_i(E) = (beta_i)_{t=t0} + Delta beta_i(E). From (11) and (10a) the perturbations of the elements are (13), boxed:

 d alpha_i/dE = d(Delta alpha_i)/dE = 1g_i cos E - 2g_i sin E,
 d beta_i/dE = d(Delta beta_i)/dE = 2g_i cos E + 1g_i sin E,

with 1g_i = -f2 P_i - f3 dx_i/dE and 2g_i = r^2 P_i - f1 dx_i/dE. (The left superscripts 1 and 2 on g are as printed; the right-hand side of the second line is read on the image as above.) I did not re-derive (13); I checked only that its structure is the standard variation-of-constants form (cos/sin rotation of a two-component forcing), and note the factors of a(E) are not displayed on the page. Equations (13), (10b), (10c), (10d) "build a system of first order differential equations for the functions alpha_i(E), beta_i(E), t(E), A_i(E), a(E) which characterizes completely any perturbed two body problem", regular with respect to the inclination and to r = 0, each containing the perturbation P_i as a factor (except (10b)), so that when P = 0 the elements stay constant. "Expressions (8) still hold when the motion is perturbed and can be advantageously used at every step to eliminate 3 quadratures because (10d) is then no longer needed."

Time element, equation (14) (READ p.437): in practice (10b) is replaced by an equation for the physical-time perturbation Delta t(E):
 d(Delta t)/dE = sqrt(a(E)) r - sqrt(a(0)) [1 - e(0) cos(E + E_p(0))].
(The subtracted term is the unperturbed d t/dE for the initial elements; this removes the secular E-linear part so that the integrated quantity is small when the perturbation is small. INFERRED from the form; the paper does not say this in words.)

### 1.6 Stability remarks (READ p.437, section 4)

- Kepler motion is only orbitally stable (the geometry is stable, the position along the orbit is not).
- The theory uses "all the stability allowed" because the geometry of the orbit is computed regularly and strictly stably (constant-coefficient linear equations (6)), which gives "a very reasonable error propagation" for alpha, beta, a, A; the position-along-the-orbit instability appears only in the integration of the physical time t(E).
- "Except for orbits with very small eccentricity, error propagation was found smaller when the vector A_i(E) was not integrated (with (10d)), but computed directly from alpha_i(E) and beta_i(E) with help of (8). This seems to be due to the reduction of the number of required numerical quadratures." This is an observation reported without data.

## 2. What is not in the paper

- No parabolic or hyperbolic case; no variable for them. The variable E needs a > 0 (a real, positive; (5) takes sqrt(a)). COMPUTED, not printed (see section 4.2): a hyperbolic analogue exists with dt = sqrt(-a) r dF and the sign of the oscillator term reversed.
- No numerical example, table, step count or error figure. The thesis was to carry those; the downstream source in the corpus is Bond and Allman 2021 Chapter 9 (digest `docs/notes/2026-06-19-digest-bond-allman-2021-modern-astrodynamics.md`, which records an e = 0.95 oblate Earth plus Moon case and a radial thrust case against Tsien's analytic solution). Those numbers belong to that book, not to this paper.
- No rectilinear (collision) example worked out; "ejection-collision" is the stated motivation only.
- No spatial (3D) specialisation beyond the vector equations, which are 3D as written.
- No perturbation example, no Hamiltonian or canonical statement, no variational equation.

## 3. Comparison

### 3.1 With Levi-Civita and KS

Levi-Civita (planar, z = w^2 on the complex plane) and KS (3D, the 4-vector map printed in the Heggie digest, section 1 step 3) change the DEPENDENT variables so that, together with dt = r ds, the Kepler problem becomes a set of harmonic oscillators of a single frequency (the frequency is set by the energy; fixed energy makes it exactly linear, with variable energy a coefficient). Their independent variable is the Sundman time with a constant scale. Burdet instead keeps x and changes independent variable AND adds a, A as redundant variables: the equation is linear in x only with the constant coefficient 1, and the energy a and Laplace vector A appear as constants-to-be-varied. The price is the 11th-order system with redundant variables (a, A) and the need for the standard element relations; the gain is that no coordinate map is needed and 2D, 3D and velocity-dependent forces are one formula. KS needs a bilinear constraint and a gauge variable; Burdet has none (INFERRED from the two sets of equations, which I have compared by their printed forms).

Burdet's equation (6a), x'' + x + a A = 0, is a 3-component unit-frequency oscillator in E with a constant forcing; KS gives a 4-component oscillator after a quadratic map. Both use the eccentric-anomaly-type fictitious time for an ellipse; I did not compare the two derivations beyond that (INFERRED).

### 3.2 With Sundman

Sundman's regularisation is dt = r ds (or r1 r2 ds in the restricted problem). Burdet's (5) is Sundman's with the energy-dependent scale sqrt(a): dt = sqrt(a) r dE. For a constant a the two time variables differ only by a constant factor. The difference is not the time but the equation: with Sundman time alone the project's `core/cr3bp_regularized.py` still integrates the Cartesian 6-state with the singular force law, whereas Burdet's formulation replaces the force law by the linear form (10a) whose coefficients are slowly varying elements. COMPUTED, 2026-10-04: with dt = r ds (scale 1 instead of sqrt(a)), the same algebra gives d^2x/ds^2 + w(s) x + A(s) = r^2 P with w = 1/a = 2/r - v^2, dw/ds = -2 f1 and dA/ds the same form as (10d) in s. This form is valid for any sign of the energy (elliptic, parabolic, hyperbolic) because it uses 1/a rather than a. It is my derivation from the printed (2), (3), (9); not printed in the paper.

### 3.3 With Heggie 1974

Heggie (digest `2026-10-04-digest-heggie-1974-global-regularisation.md`) regularises a many-body problem with KS maps on each pair separation, in Hamiltonian form with an enlarged phase space, and gives a restricted-problem limit with the primaries on any Kepler orbit. Burdet is a single-centre perturbed-two-body regularisation in non-canonical form. For the project's single close pass of a moon (a perturbed two-body problem about the moon, the Earth and the rotating frame being a perturbation) Burdet's form is the cheaper one: three scalar equations, no map, no gauge. For a pass that is close to BOTH primaries or for a time-dependent (elliptic) two-primary geometry, Heggie's (39) and (46) are the stated route. Both authors use a redundant extended state; Heggie's appendix gives a theorem for it (any such enlargement preserves the solution if the initial data satisfy the constraints), Burdet gives no such statement and relies on explicit relations (2), (3), (8). Neither paper gives a numerical test of the near-collision regime in the restricted problem.

## 4. Printed numbers

None. The paper contains no numerical result. Structural statements usable as exact checks: the eleventh-order count (p.437), the standard-form identities sum alpha^2 = 1, sum beta^2 = 1 - e^2, sum alpha.beta = 0, p = a (1 - e^2), A = e alpha (p.436), and the vis-viva and Laplace definitions (2), (3).

## 5. Reconciliation with project code

Searched on 2026-10-04: `src/cyclerfinder/core/cr3bp_regularized.py`, `genome/multi_shooting.py` header, the #670 scripts, and the notes.
- `core/cr3bp_regularized.py`: Sundman only, dt/ds = r1 r2, r1 or r2, augmented state (x, y, z, vx, vy, vz, t) integrated in s, Cartesian and singular force law (its docstring cites Szebehely 1967 section 3.7 and Stiefel and Scheifele 1971). It does not carry any element, and has no a(E)/A(E) variable or linear form. Nothing in it contradicts Burdet; Burdet's (5) is the same type of time transformation with a different scale.
- #670 (`scripts/certify_670_wz_oterma_regularized.py`, `scripts/_validated_taylor_integrator.py`): planar Levi-Civita, z = w^2 about the secondary, dt = r2 dtau on the extended-phase-space Hamiltonian Gamma = r2 (H - h), interval arithmetic. This is the only coordinate-regularised propagator in the repository and is a proof tool, not a `src/` float propagator. Burdet's formulation is not implemented anywhere in the repository (grep for Burdet, Sperling and "Laplace vector" in `src/` found no implementation; the only repository mentions are the notes digests listed above).
- `genome/multi_shooting.py`: header records that a regularised STM is not implemented and uses a dual integration; Burdet's element form is a candidate (section 6 below) but nothing here is built.
- Constants: none of the paper's constants is in the code. No reconciliation of a numeric value is needed because the paper has none.
- Citation reconciliation: the repository's `docs/notes/2026-06-19-digest-bond-allman-2021-modern-astrodynamics.md` cites "Burdet (1968)" for the Sperling-Burdet system; this paper is Burdet 1967 (ZAMP 18:434), the announcement. The 1968 work I RECALL as C. A. Burdet, "Le mouvement keplerien et les oscillateurs harmoniques", J. reine angew. Math. 238 (1969), with a 1968 ZAMP/Bull. note on a related theme; I have not seen these and the dates are from memory (RECALLED). Do not cite either until the corpus holds them.

## 6. Techniques applicable to the project's problems

All uses are INFERRED; nothing was built or run. Throughout, the attracting centre is the moon being passed and everything else (the other primary, the rotating frame) is the perturbation P.

### 6.1 #928 (a regularised propagator and transition matrix for close passes)

Which formulation to build. Three candidates: Levi-Civita (planar; already derived in #670 in intervals), KS (3D; the Heggie digest section 1 has the map), and the Burdet/Sperling element form (this paper, any dimension). Recommendation:
1. Planar first: Levi-Civita, as #928 already says, because the #670 derivation exists and the Hamiltonian form gives a symplectic check on the transition matrix. In the rotating frame the perturbation set is the other primary's attraction plus the Coriolis and centrifugal terms; Levi-Civita handles it through the energy (Jacobi) integral as #670 does.
2. 3D: choose between KS (4-vector, a gauge relation to maintain, Heggie's restricted-problem form available) and Burdet (3-vector, no gauge, the element system). Burdet's form is the shorter to implement as a first prototype because it needs only (6), (10), the element initialisation (2), (3) or (8), and a time equation; its weak point is the elliptic-only variable E: a pass of a moon from the heliocentric or planetocentric side is HYPERBOLIC with respect to the moon (positive two-body energy about the moon). So the version to build is the energy-independent form of section 3.2, d^2x/ds^2 + w x + A = r^2 P with dt = r ds, w = 1/a, and the rates dw/ds = -2 P.(dx/ds), dA/ds = 2 (P.dx/ds) x - (x.P) dx/ds - (x.dx/ds) P (a COMPUTED derivation of this note, unverified numerically). That form is valid for hyperbolic and parabolic passes. It must be tested on the controls below before it is trusted; the paper's stability remark (compute A from alpha, beta instead of integrating it) was only observed for the elliptic form.
3. Controls (analytic, in order):
 (i) A Kepler ellipse in the regularised variables is reproduced exactly: with P = 0, the elements are constant and x(E) = [alpha cos E + beta sin E - A] a. The test is that numerical integration of (6a) or (10) with P = 0 holds alpha, beta, a, A constant to the integrator tolerance, and that physical time from (6b) equals Kepler's equation t = a^(3/2) (E - e sin E), against `core/kepler.py` to 1e-12 (the same control as for the Levi-Civita version in the Heggie digest section 8).
 (ii) Radial fall (the collision orbit): e = 1 gives x(E) = -a (1 - cos E) alpha_hat and r = a (1 - cos E), t = a^(3/2) (E - sin E), the cycloid; pass through r = 0 at E = 0 and 2 pi must be smooth (x a double zero, the orbit staying on one side of the centre, which is the usual reflecting continuation of Levi-Civita and KS) and the time equation exact to the integrator tolerance. COMPUTED from (7a) and (8a, 8c) in standard form with sum beta^2 = 0; not printed in the paper. A hyperbolic radial case needs the section 3.2 form, whose closed form (cosh) is a COMPUTED extension and must be checked against the radial Kepler problem before use.
 (iii) A near-collision case at e = 1 - 1e-8 against `core/kepler.py` through the pericentre passage; the unregularised propagator fails or needs a tiny step there, so this is also the demonstration of the regularisation.
 (iv) A perturbed control: add a small constant perturbation P (a uniform acceleration) and compare with a plain high-accuracy integration of (9) at a pass distance above 1e-2, then at 1e-3 and 1e-4 (the same ladder as #928).
 (v) Determinant and finite-difference Jacobian of the transition matrix (below).
Transition matrix in the Burdet form: the state is (x, dx/ds, w, A, t) with algebraic constraints (2), (3); the variational equations of (10a), (10c), (10d) are linear and regular, with P-dependent terms O(P), so the matrix is regular through r = 0. The paper does not give them. Two cautions: the state is redundant, so the transition matrix of the redundant state is larger than 6 x 6 and must be projected onto the 6 physical components (or integrate only the 6 + t independent ones with A computed from (8)); and a symplectic check is not available because the formulation is not canonical (Heggie's and Levi-Civita's are). So choose Levi-Civita for the first float transition matrix and use Burdet as an independent cross-check with a different structure, which is a stronger test than a second Sundman integrator.

### 6.2 #899 (continuation through near-collision seeds)

Seeds from the second-species theory (Guillaume, Henon, Perko, Casoliva) start at mu small, with the passes at distances of order mu. A regularised propagator is what makes the corrector survive these. The Burdet form is attractive here because the perturbation is genuinely small near the moon (the other primary's tidal term is small at the moon) so the elements change slowly: an integrator step in E or s is nearly free, and the main error is in the time equation, as the paper states. Use: propagate the near-moon arc in the regularised variables with a pass radius below a threshold and switch to the plain integrator outside; the junction needs only x, v and t, which are available from (11a, 11b) and (10b). Gate: results agree with plain DOP853 on passes above 1e-3 and with the Levi-Civita propagator (when built) on passes below 1e-4; record any disagreement, since the `#896` tolerances were set by two integrators of the same structure.

### 6.3 #924 (variational vectors through a close pass: FLI and MEGNO)

The FLI and MEGNO screens integrate a variational vector along the orbit; through a pass the plain variational vector grows as 1/r^3 and the indicator is dominated by the pass. In Burdet's variables the variational equations are regular, so the screen can be computed in regularised time without the pass dominating. The caution in section 6.1 applies: the redundant state must be projected; and the indicators are defined in physical time, so the conversion through dt = sqrt(a) r dE or dt = r ds must be applied before averaging, as a MEGNO average in fictitious time weights the pass by its (small) duration differently. Control: the same finite-difference Jacobian check, plus the Kepler ellipse where the exact variational solution is known from the constant elements (the transition matrix of a Kepler orbit is obtained by differentiating (7) with respect to the six constants, COMPUTED as a statement of structure, not carried out).

### 6.4 #929 (integrator controls through one moon flyby)

The paper does not bear on the controls (a) to (e) of #929 except as one more independent integrator: a regularised propagator of different structure is the right third member of the three-way comparison in (b), and the analytic-against-finite-difference transition matrix check in (e) is the control of section 6.1 (v). The Burdet form adds a check that none of the others has: A (the Laplace vector) and a computed from the integrated state through (2), (3) must agree with their integrated values (10c), (10d) (the paper notes the integrated and direct values differ in error propagation); the difference is an error monitor that costs nothing, and is the natural analogue of an energy drift monitor in a time-dependent model (#929 says there is no invariant to watch in most models). INFERRED.

## 7. Recommended follow-ups (not registered)

1. Implement Burdet's element system in its energy-independent (dt = r ds, w = 1/a) form as a float prototype in a scratch module, run controls (i) to (iv) of section 6.1, and decide whether it is worth carrying into `src/` as a cross-check on the Levi-Civita propagator of #928. Verify the section 3.2 and 6.1 derivation (the 1/a form and the rate equations) by symbolic differentiation first; it is hand algebra here.
2. Derive and test the Burdet variational equations (the 6 + t independent components) and compare with the Levi-Civita transition matrix on a near-collision orbit.
3. Obtain the follow-up papers the author announces (the thesis and the 1968 and 1969 articles, RECALLED titles in section 5) and Stiefel's chapter (reference [3]) and Sperling 1961 (reference [2]); none is held. Bond and Allman 2021 Chapter 9 is held and carries the complete Sperling-Burdet system and numerical tables; the digest of 2026-06-19 summarises it but the Tables 9.2 and 9.4 values could be used as sourced controls for an implementation of the formulation, after reading the book pages against the code (not done here).
4. Record the unresolved elliptic-only limitation (paper p.435) in the section of #928 that chooses the formulation: for a hyperbolic moon pass the paper's E variable does not apply as printed.
5. When the turn gate (#906) is next touched, take the integrated deflection from the regularised propagator (as in the Heggie digest) and compare with the first-order relation; Burdet's elements give the deflection directly as the change of the Laplace vector direction through the pass (INFERRED, since A is the eccentricity vector of the instantaneous osculating orbit).

References printed on p.438 (READ): [1] P. Kustaanheimo and E. Stiefel, "Perturbation Theory of Kepler Motion Based on Spinor Regularisation", J. reine angew. Math. 218:204 (1965); [2] H. Sperling, "Computation of Keplerian Conic Sections", ARS J., May 1961, pp.660-661; [3] E. Stiefel, "Many-body Problem and Interplanetary Flight", in Dynamics of Rockets and Satellites, G. V. Groves ed., North-Holland. None of the three is held in the corpus.

## Note of 2026-10-05: the full paper (Burdet 1968) has been read

The full paper behind this announcement, Burdet 1968, ZAMP 19:345-368, is digested in `docs/notes/2026-10-04-digest-burdet-1968-theory-kepler-motion-perturbed-two-body.md`. Consequences for this note:
- The elliptic-only limit of the printed 1967 form (section 2 and section 6.1 above) is removed in 1968: the general form uses dt = (1/kappa) r ds with omega^2 = 1/a, valid for elliptic, parabolic and hyperbolic motion. Its equations (111), (201), (202), (203) are the energy-independent form that sections 3.2 and 6.1 of this note derived by hand; the derivation is confirmed term for term and numerically (including a perturbed hyperbolic-energy close pass), so the "check it symbolically" caveat is closed.
- Sign and scale conventions differ: alpha(1968) = -alpha(1967) and beta(1968) = beta(1967)/omega with E = omega s; the vector A is the same.
- The 1968 element equations as printed have three errors (found by numerical differentiation of the exact elements): the alpha term of gamma' has the wrong sign in (206c), (311) and Appendix II (6), and (206d) lacks the tilde on B3.
- Section 5's recollection of a 1968 title was imprecise: the 1968 paper is the one digested above; "Le mouvement keplerien et les oscillateurs harmoniques" (J. reine angew. Math. 238, 1969) is still only a recollection.
- The 1968 paper prints one usable benchmark (an exact L4 perturbed two-body problem with three printed orbits and periods) and no tables; see its section 3.
