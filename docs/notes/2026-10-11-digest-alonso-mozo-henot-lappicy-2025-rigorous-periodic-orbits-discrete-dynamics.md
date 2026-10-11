# Digest: Alonso Mozo, Hénot and Lappicy 2025, "A Practical Guide to Rigorously Locate Periodic Orbits in Discrete Dynamics" (arXiv:2510.06373v1) (#960, owner-requested)

Lucia Alonso Mozo (Universidad Complutense de Madrid; now Google Zurich), Olivier Hénot (National Taiwan University) and Phillipo Lappicy
(Universidad Complutense de Madrid). arXiv:2510.06373v1 [math.DS], submitted 7 October 2025, DOI 10.48550/arXiv.2510.06373, MSC2020 37E05, 37G15,
37N25, 65G30, 92D25. Not a wanted-list row; owner-requested. A method paper in computer-assisted proof for discrete maps, not celestial mechanics.

- **Filed as** `cyclers_pdf/papers/alonso-mozo-henot-lappicy-2025-practical-guide-rigorously-locate-periodic-orbits-discrete-dynamics-arxiv-2510.06373v1.pdf`.
  Numbers: `data/sources/alonso-mozo-henot-lappicy-2025-numbers.yaml`.
  - arXiv PDF as received from the lead, md5 47f9aee19a5e864e6ee137d54e357fe7, 17 pages, 6,410,193 bytes, text layer present, no OCR.
  - Reproduction of the worked example: `...-check_logistic_example.py` and `.out`, beside the PDF in the private papers directory.
- **How it was read:** the whole text layer (all 17 pages: Secs. 1-3 and the reference list). No tables exist. I did not read the figures beyond their captions.
  The worked example (Sec. 2.1) was reproduced in exact rational arithmetic as a positive control (sec. 3 below). The Julia code repository was not fetched.

## 0. Verdict

**A short, self-contained recipe for proving that a numerically found periodic orbit of a map really exists, is unique in an explicit ball, and has a
rigorous error bound. It is the radii-polynomial (Newton-Kantorovich type) contraction test on a quasi-Newton operator, with interval arithmetic for the
bounds. The demonstrations are two one-dimensional maps. It does not touch ODEs, Hamiltonian systems or celestial mechanics.**
- Method in three steps: (1) rewrite a p-periodic orbit as a zero of F in R^{np} by treating all p iterates as unknowns (Eq. 3); (2) find an approximate zero
  xbar and an approximate inverse A of DF(xbar) by ordinary floating point; (3) bound Y, Z1, Z2 rigorously by interval arithmetic and check that
  P(r) = Y + r (Z1 - 1) + r^2 Z2 / 2 <= 0 with Z1 + r Z2 < 1 for some r.
- A passing check gives existence of a unique zero of F in the ball B_r(xbar), and so a rigorous a posteriori error bound r on the numerical orbit. Stability
  follows from interval enclosure of the eigenvalues of Df^p (Gershgorin disks if n > 1). Period-doubling curves in a parameter plane are followed with a
  uniform version of the same contraction (Theorem 2) and Chebyshev interpolation.
- The examples: the logistic map (8906 unstable and 372 stable orbits of periods 2 to 80 proven) and a one-dimensional predator-prey map (820,921 unstable and
  143,477 stable orbits of periods 2 to 10, plus period-doubling curves for kappa in [-31, -13] with error at most 1e-4 per window).

## 1. The method stated precisely

- **Setting:** x_k = f(x_{k-1}) with f a diffeomorphism of R^n (a discrete model, or "the Poincare map of a differential equation", p.2; no ODE example is given).
- **Operator:** T(x) = x - A F(x), A linear and injective with A approx DF(xbar)^-1 (Eqs. 4-5). The paper notes A = DF(xbar)^-1 would give Z1 = 0 but costs an
  interval matrix inverse; a floating-point inverse wrapped in intervals is enough.
- **Theorem 1** (Banach space X, F of class C^1, A injective linear). Suppose Y, Z1, Z2 = Z2(R) satisfy
  ||A F(xbar)|| <= Y, ||I - A DF(xbar)|| <= Z1, ||A (DF(x) - DF(xbar))|| <= Z2 ||x - xbar|| for x in B_R(xbar).
  If some r in [0, R] gives P(r) <= 0 and Z1 + r Z2 < 1, a unique zero of F lies in B_r(xbar). Proof deferred to Breden's thesis (Ref. 2, Thm 1.2.1).
  - Reading: Y is the residual of xbar, Z1 is the quality of A, Z2 is a Lipschitz constant of A DF, R is an a priori radius.
  - The admissible radii form [r-, r+] with r+- = (1 - Z1 +- sqrt((1 - Z1)^2 - 2 Y Z2)) / Z2. The smallest is the tightest error bound, the largest is the
    biggest ball where uniqueness is guaranteed.
  - If Z1 < 1, ADF(xbar) is invertible, so A is injective. This is a finite-dimensional check.
- **Norms:** the 1-norm on R^{np}, and the maximum absolute column sum for operators (Eq. 8).
- **A genuine period p** (not a divisor): the p balls B_r(xbar_k) must be disjoint (Eq. 9).
- **Interval arithmetic (Sec. 1.1):** intervals [a, b] with floating-point endpoints, outward rounding of +, -, x, /, and elementary functions. The paper's tooling
  is Julia with IntervalArithmetic.jl; INTLAB (Matlab) and MPFI (C++) are named. IEEE 1788-2015 is cited.
- **What "rigorous" means here:** the proof is conditional on the correctness of the interval-arithmetic implementation and on F and DF being coded as intended.
  The conclusion is a theorem about the exact F, not a numerical tolerance. The certificate says nothing about whether F is the right equation.
- **Nonlinearities (Sec. 3.1-3.2):** where f is not polynomial, a Taylor expansion of order N with a rigorously bounded remainder (here N = 10, using
  sup |h'''| = 1/8 and sup |h''| = sqrt(3)/18 for the sigmoid). The authors say automatic differentiation or Poisson summation would control this better, and
  that a Banach-algebra setting on Chebyshev coefficients (l^1) is the longer-term route (Refs. 3, 13).
- **Limits stated by the authors:** the non-degeneracy of the eigenvalue crossing -1 is not proved, so the period-doubling curves are candidates. At a fold the
  implicit function theorem fails and Theorem 2 does not apply; pseudo-arclength continuation is not done (Refs. 35, 3).

## 2. The examples (all numbers in the numbers file)

- **Logistic map, mu = 3.2, period 2:** xbar = (0.51, 0.79) gives Y in [0.012775, 0.0127751], Z1 in [0, 6.66e-16], Z2 in [20.7422, 20.7423] and r* = 2^-6, and
  the eigenvalue lies in [-0.0644712, 0.314457], so the orbit is stable. With a 16-digit xbar the error is about 1e-12.
- **Logistic map sweep:** 8906 unstable and 372 stable orbits for p = 2 to 80, consistent with the period-doubling cascade (Fig. 1).
- **Predator-prey map** f(x) = beta + x - kappa / (1 + e^x): periodic orbits of periods 2 to 10 over a grid in (beta, kappa), then period-doubling curves ("shrimp"
  structures) by uniform contraction with K = 24 Chebyshev nodes and Taylor order 10.

## 3. Our checks (all computed by us; not from the authors' Julia code)

- **Positive control on the worked example** (`check_logistic_example.py`). Exact rational arithmetic stands in for interval arithmetic.
  - It reproduces every printed enclosure: Y = 0.0127750744, Z1 = 3.3e-16 (below the printed upper bound 6.7e-16), Z2 = 20.742247, the roots r- = 0.015158 and
    r+ = 0.081264 with r* = 2^-6 inside, P(r*) = -3.2e-4, and the eigenvalue range [-0.0644712, 0.3144566] over the ball of radius r_m. The true period-2 points at
    mu = 3.2 are 0.5130445 and 0.7994555 and the true eigenvalue is 0.16.
  - **One print slip found.** The A of Eq. (16) and the DF of Eq. (17) are the matrices for the rows in the order (f(x_0) - x_1, f(x_1) - x_0). Eq. (13)-(14) as
    written have the other order. Used with Eq. (14) as printed, Eq. (16) gives Z1 = 4.16 and the contraction fails. Used with Eq. (17), Z1 = 1.9e-12. The printed
    enclosures of Eq. (19) reproduce from the code-snippet order, so no result is affected. It looks like a typesetting slip in a worked example, the kind that
    happens to everyone; it is worth knowing only for a reader who retypes the matrices.
- **Counts:** 8906 + 372 = 9278 and 820,921 + 143,477 = 964,398, above the abstract's lower bounds of 8,000 and 800,000.
- The "positive control" matters here: a certificate routine that returned "proved" on everything would also pass the paper's statements. The control above
  shows our reading of Theorem 1 gives the paper's numbers and, by the row-swap, fails when it should.

## 4. Relevance to this project

**Short version: a rigorous existence certificate is stronger than our V1 and V2 re-fly checks for the one property it covers (a true periodic orbit of the
row's defining model exists near the numerical one, with an explicit error radius). The method as published cannot be applied to our CR3BP orbits without
a validated ODE integrator and a phase or energy constraint, and for the patched-conic date maps it is close to what the project already built in #663.**

- **What it would certify, and what it would not.**
  - It would certify that F = 0 has a unique solution in a ball, for the F we write. V1 (spec sec. 14: solver cross-check and Kepler re-propagation) and V2
    (multi-lap bounded drift) are numerical agreements; they do not exclude a nearby spurious closure or a missed root. A certificate does.
  - It says nothing about model fidelity (V3, V4: ephemeris, GMAT) or novelty (V5). The spec's own note that spectral CR3BP stability does not imply V4 survival
    applies in full: a proven periodic orbit of the CR3BP is still an orbit of the CR3BP.
  - **It certifies only the F you wrote.** The project's orbit-closure discipline (memory note "orbit-closure-discipline"; the S1L1 closure notes
    `2026-06-08-s1l1-corrected-closure-results.md` and `2026-07-12-572-closure-adjudication.md`) records closures declared on a subset of the binding constraints
    (#164: V-infinity, time of flight and radius with no longitude rendezvous; the orbit missed real Mars by 2.6 AU). A rigorous proof of a zero of a
    residual that omits a constraint is a rigorous proof of the wrong thing. Any use would need the residual reviewed against the full constraint list first.
    The rule "all binding constraints in the residual" is the gate; the certificate is not a substitute.
- **CR3BP resonant families (Poincare-map or multiple-shooting orbits).**
  - **Fit of the formulation:** Eq. (3), all iterates as unknowns, is multiple shooting. The project already uses it for strongly unstable orbits (#1000, and the
    arc-by-arc closure of #1047 rows A, B, D with |lambda| about 110 to 720). Short arcs keep Z2 manageable, which suggests the contraction test is feasible
    in principle for those rows. I have not tested this.
  - **Two gaps the paper does not address (our analysis):**
    1. *A validated ODE integrator.* f is a flow or return map, so Y needs interval-valued flow arcs and Z2 needs bounds on the second variation. Nothing in the
       paper does this; the authors cite Tucker 2002 only among the achievements of the field (Sec. 1) and do not mention CAPD. The #653 scoping note records the same blocker:
       the codebase has no rigorous ODE integrator, and the Wilczak-Zgliczynski proofs rest on the CAPD C++ package.
    2. *Non-isolated roots.* A periodic orbit of an autonomous Hamiltonian flow comes in a one-parameter family and has a unit eigenvalue pair, so DF is singular
       unless the Jacobi constant or a coordinate is fixed (a section plus a fixed-energy constraint). Their examples are one-dimensional dissipative maps, where this
       does not arise. This is the same reason our own fixed-h checks (the 2005 Voyatzis-Kotoulas control) were needed. Also |lambda| = 1 is "inconclusive" for stability.
- **Patched-conic date maps (the two-working-body generator and the #563/#600 symmetric closure).**
  - Here f is built from Lambert solves, a transcendental map, not an ODE. #663 already built a certified interval-Newton/Krawczyk closure test on an augmented
    system (the Lambert universal variable promoted to an unknown), using `mpmath.iv` (task #610's optional `interval` extra). Krawczyk and the radii-polynomial
    contraction are close relatives: both give existence and uniqueness in a box. The paper's version adds an explicit a posteriori radius range [r-, r+]
    and a route to parameter continuation and to infinite dimensions. That comparison is my reading, not the authors' claim.
  - Open obstacle already recorded in the project: bounding a multi-revolution, branch-selecting Lambert solve rigorously ("still unsolved", #610 and #625). This
    paper does not touch it. The 2026-10 topology-gap finding (VEM legs need 3 to 8 revolutions) makes branch selection the real difficulty.
- **Cost to adopt.**
  - Python side: `mpmath.iv` is already an optional extra (`interval` in pyproject). It is slow and has no matrix or ball-arithmetic tools. A ball-arithmetic library
    for Python (arb through python-flint, if available for 3.11) would be a candidate; I have not checked that it installs or that it covers the functions needed.
  - The authors' tooling is Julia and no Julia is installed on this machine (checked: no `julia` on PATH). The Julia route means a second language for one tool.
  - For ODE orbits, the minimum is a validated integrator (CAPD bindings, or a Taylor-model code), which #653 estimated at weeks.
  - For algebraic or Lambert maps, the existing Krawczyk routine could be reused with the contraction test layered on it, a much smaller job.
- **Where it bears on the project's gates.** Spec sec. 14: it could sit above V1 as an optional "certified" annotation on the idealized-model closure of a row,
  never as a replacement for V2 to V4 and never as novelty evidence. Closure notes and the orbit-closure memory: it would be an independent cross-check only if run
  on a residual built independently of the finder's, with an independent integrator or arithmetic, per the rule that the checker takes constants from the reference model.
- **Catalogue implication:** none. No row is a cycler. Treated as method background; no recommendation to adopt, only a scoped note of what adoption would need.

## 5. Citation mining (against CORPUS_INDEX.md and the wanted list, 2026-10-11; proposals only, no row edited)

- **Held:** Wilczak and Zgliczynski, heteroclinic connections in the PCR3BP, Parts I and II (the CAPD-based proofs this paper's reference list leans on);
  Capinski, Gidea and de la Llave 2016; Mireles James and Murray 2017 and Gonzalez and Mireles James 2017 (parameterization methods). The paper cites none of these.
- **Proposed rows (celestial-mechanics or periodic-orbit relevance, neither held nor wanted):**
  1. Calleja, Garcia-Azpeitia, Hénot, Lessard and Mireles James (2024), "From the Lagrange triangle to the figure eight choreography: proof of Marchal's conjecture",
     arXiv:2406.17564. The one rigorous periodic-orbit result for an N-body problem in the list, co-authored by Hénot, one of this paper's authors. Medium priority.
  2. Hungria, Lessard and Mireles James (2016), "Rigorous numerics for analytic solutions of differential equations: the radii polynomial approach",
     Math. Comp. 85:1427-1459. The standard radii-polynomial source; the paper's Theorem 1 is its finite-dimensional form. Medium.
  3. van den Berg, Lessard and Mischaikow (2010), "Global smooth solution curves using rigorous branch following", Math. Comp. 79:1565-1584. Rigorous continuation of
     solution branches, the route to family-level certificates. Medium.
  4. Tucker (2002), "A rigorous ODE solver and Smale's 14th problem", Found. Comput. Math. 2:53-117. The validated ODE solver whose absence in the project is the
     blocker in sec. 4. Low to medium.
  5. Breden (2025), "Computer-assisted proofs for differential equations and dynamical systems", thesis, Institut Polytechnique de Paris (tel-05142432v1), the proof of
     Theorem 1. Low; open access, probably easy to fetch.
- **Not proposed:** Mireles James and Mischaikow 2015 and van den Berg and Lessard 2015 (surveys; useful only if the method is adopted), Wilczak 2003 (Kuramoto-Sivashinsky),
  Zgliczynski 1997, Jaquette 2019, Lanford 1982, Mischaikow and Mrozek 1995, Nakao-Plum-Watanabe 2019, Gomez-Serrano 2019, Moore 1966, Tucker 2011 (book), IEEE 1788-2015,
  Church et al. 2025 (arXiv:2504.03058, predator-prey), Kryzhevich et al. 2021, and the logistic-map, bifurcation and "shrimp" references (Lorenz, May, Li and Yorke,
  Devaney, Field and Golubitsky, Sander and Yorke, Gallas, Stoop et al., Facanha et al.). No celestial-mechanics content.

*Filed as `cyclers_pdf/papers/alonso-mozo-henot-lappicy-2025-practical-guide-rigorously-locate-periodic-orbits-discrete-dynamics-arxiv-2510.06373v1.pdf`. Numbers: `data/sources/alonso-mozo-henot-lappicy-2025-numbers.yaml`.*
