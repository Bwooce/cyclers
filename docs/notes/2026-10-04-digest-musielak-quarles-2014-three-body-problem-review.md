# Digest: Musielak & Quarles 2014, "The three-body problem" (review)

Date: 2026-10-04 (Sydney). Reading only; no code or computation. Source: Z. E. Musielak and B. Quarles, "The three-body
problem", Reports on Progress in Physics 77:065901 (2014), DOI 10.1088/0034-4885/77/6/065901; the held copy is the arXiv
version 1508.02312v1, 68 pages, filed in the private paper corpus as
`musielak-quarles-2014-three-body-problem-review-rep-prog-phys-77-065901-arxiv-1508.02312v1.pdf`. Page numbers below are
the arXiv page numbers printed in the running head. Tags: READ (p.N) means read in the paper's text layer; INFERRED
means my reasoning, not a statement of the paper.

## 1. Scope and level

READ (pp.1-6, 52-68, and the sections cited below). A broad survey for astrophysicists and astrodynamicists: general
three-body problem (sections 2 to 5), circular restricted problem (6), elliptic restricted problem (7), applications of
the restricted problems with weight on exoplanets (8), the relativistic problem (9), and a closing summary (10). About
half the text is history from Newton and Poincare to the present; the mathematics is stated, not derived, and each
method gets one to three paragraphs plus pointers. The review is not written for trajectory designers: there is no
treatment of flybys, gravity assists, cyclers or mission design beyond one paragraph (section 8.1, pp.39-40). Its
value to this project is as a map of where to look for chaos diagnostics and for second-species and continuation
literature, not as a source of procedure. The reference list (pp.50-68) is the most reusable part.

## 2. Techniques applicable to the project's problems

The project's synthesis inventory is `docs/notes/2026-10-04-897-technique-synthesis-papers-to-problems.md` (T1 to T19).
Everything below is either absent from it or present only in passing. A search of `src/cyclerfinder` for Lyapunov,
MEGNO, FLI finds no chaos-indicator code (INFERRED from a text search; one eigenvalue-based classifier comment in
`ml/seed_generation.py:599` is not a Lyapunov-exponent computation).

### 2.1 Chaos indicators: maximum Lyapunov exponent, FLI, MEGNO (sections 4.3 to 4.5, pp.23-25). New to the inventory.

What the review says, READ:

- Maximum Lyapunov exponent (4.3, p.23): "A positive Lyapnuov exponent lambda denotes that a trajectory ... is diverging
  from the neighboring orbit"; definition lambda = lim 1/(t - t0) ln(delta(t)/delta(t0)) (eq. 25). First applied to the
  general problem by Benettin et al. 1976, 1980 and to the restricted problem by Jeffreys and Yi 1983 (p.23). The
  Jacobian of the acceleration field (eq. 26) is what the tangent-space integration needs. (The review's statement that a
  Hamiltonian system "has 3 positive and 3 negative exponents" in the CR3BP is loose: there are zero exponents from the
  integrals, INFERRED.)
- FLI (4.4, p.24): the stated motivation is that exponents cost "extra computational time" and have "normalization and
  dimensionality" problems (Froeschle, Lega, Gonczi 1997); the indicator is "FLI(t) = sup_j ||v_j(t)||" over the
  tangent vectors (eq. 27). "The slope of the time series of FLI are used as indicators of chaos. If the slope is steep
  and positive, then chaos is inferred ... The slope may also remain flat ... which is the stable case." The review says
  chaos can be detected "within 100,000 year timescales; note that previous methods required much longer timescales"
  (p.24).
- MEGNO (4.5, pp.24-25): Y(t) = (2/t) integral_0^t (delta-dot/delta) s ds, with time average <Y>(t) = (1/t) integral Y
  (eqs. 28, 29); two extra differential equations (eq. 30). Stated advantage: "<Y>(t) converges faster to its limit value
  ... the time-weighting factor amplifies the presence of stochastic behavior, which allows the early detection ... of
  chaos" (p.25). Primary sources: Cincotta and Simo 1999, 2000; applications Gozdziewski et al. 2001, 2002; Goździewski
  et al. 2013 and Migaszewski et al. 2012 use it for periodic-orbit and resonance searches (p.22).
- Time-frequency diagnostics (6.7, p.36): "Vela-Arevalo and Mardsen [2004] demonstrated that a method based on
  time-varying frequencies is a good diagnostic tool to distinguish chaotic trajectories from regular ones."

What the review does not give (READ, absence): no threshold values, no statement of how the indicator behaves on a
regular orbit versus a chaotic one beyond "flat" versus "steep", no treatment of renormalisation, and no discussion of
non-autonomous or finite-duration problems. The asymptotic MEGNO values (<Y> tending to 2 for regular quasi-periodic
motion, growing linearly with t for chaotic) are in Cincotta and Simo's papers, not in this review; I state them from
memory of that literature and they must be checked at the source before use (INFERRED).

Where it applies to the project (all INFERRED):

1. A classification tool for the question the project keeps meeting: is a computed object a regular quasi-periodic
   motion near a torus, or a chaotic one that happens to stay bounded for the integration time. The `#895` real-ephemeris
   arcs and the quasi-cycler rows are the cases. The project's present substitute is bounded drift-oscillation versus
   divergence (the `#339` criterion), which is a trajectory-level test; FLI and MEGNO would add a phase-space-neighbourhood
   test using the variational equations the project already propagates (the STM lane in `nbody/shooter.py`, the
   variational Jacobian in the `#895` shooter).
2. A caution that cuts against over-use. A cycler is hyperbolic by construction: each close flyby amplifies a
   perturbation. A positive finite-time exponent on a multi-flyby arc is therefore expected and does not disqualify it.
   The informative quantity is how the growth compares across a family or across epochs, and whether it is exponential or
   linear over many cycles. Short arcs (the `#895` chain is about seven encounters over a few tens of days) do not give
   enough revolutions for FLI or MEGNO to separate regular from chaotic; the review's own framing is thousands to
   millions of revolutions of a nearly Keplerian system. The tools suit stability maps over a grid of initial conditions
   (the weak-stability and capture sweeps `#378`, `#681`, `#908`, which synthesis item T14 already covers by a different
   criterion) better than single multi-flyby arcs.
3. A close encounter makes the tangent vector grow by large factors in one step, so a Lyapunov or MEGNO integral across
   a flyby is dominated by that one passage. Using the regularised variables of 2.2 across the passage, or computing the
   indicator between encounters on the return map, would be needed; the review says nothing on this and I know of no
   statement in the held papers that does.
4. A cheap cross-check for the Floquet-multiplier approach already used for periodic orbits: for a periodic orbit the
   largest multiplier gives the exponent exactly, so the indicator is a positive control (a known-unstable orbit must show
   exponential FLI growth, a known-stable one a flat slope) before it is applied to anything uncertain.

### 2.2 Collision regularisations (section 3.4, pp.16-17; 3.5, pp.17-19). Partly in the inventory (T17).

What the review says, READ: the section is historical, not procedural. Levi-Civita 1903 treated the restricted problem
and Bisconcini 1906 the general one; both discussed "regularization that allows extending possible solutions beyond a
singularity" (p.16). "Regularization means that the motion is extended beyond a singularity ... through an elastic bounce
and without any loss or gain of energy" (p.16). Siegel 1941 showed an analytic solution through triple collision cannot be
found "for practically all masses" (p.16), extended by McGehee 1974 (collision manifold). "There are two currently used
regularization schemes, the so-called K-S regularization scheme developed by Kustaanheimo and Stiefel [1965], and the
B-H regularization scheme developed by Burdet [1967] and Heggie [1973, 1976]. These two analytical schemes describe close
two-body encounters and they are typically used to supplement numerical simulations [Valtonen and Karttunen, 2006]"
(p.16). The review never writes down the Levi-Civita or K-S transformation, so it does not supply a method; it supplies the
primary citations. Sundman (3.5, pp.17-19): the independent variable u = integral dt/t (eq. 16) removes a double-collision
singularity, with the coordinates expandable in powers of (t - t1)^(1/2); the complete series solution "requires millions
of terms to find the motion of one body for insignificantly short durations of time", "no practical applications" (p.18).
That is a statement that Sundman's series should not be pursued.

Where it applies (INFERRED): the project already has `core/cr3bp_regularized.py` and the inventory entry T17. The new
information is the pointer to the Burdet-Heggie "time-transformation only" family (Heggie 1973, 1976), which regularises
by a change of independent variable and needs no coordinate change; that is the lighter route for a numerical integrator
that must pass a flyby at periapsis distances of tens of kilometres, and for the unsoftened-moon propagator the `#895` and
`#900` lanes need. For a flyby of a small body inside a restricted model the Levi-Civita form in the project is the
right tool and the review adds nothing to it. The Sundman time transformation is the general-problem analogue and could be
used as a step-size control (a fictitious-time step) in the full n-body shooter, where the project reports compute
infeasibility (single shoot over 400 s).

### 2.3 Numerical search for periodic orbits and resonances (section 4.2, pp.21-22; 6.5, pp.31-33). Mostly new.

READ: (a) Resonance location from perturbation theory: "the nominal resonance location, Ra ... defined in terms of two
integers (k, l), where l is the order the resonance to form a ratio of k/(k + l)", Ra = (k/(k+l))^(2/3) (mp/(mp+ms))^(1/3)... as
printed eq. 23 (the printed exponent layout is garbled in the text layer; the quantity is the Kepler-third-law radius for
the period ratio) (p.22). (b) The resonant argument phi = j1 lambda' + j2 lambda + j3 varpi' + j4 varpi + j5 Omega' + j6
Omega (eq. 24) with the d'Alembert rules: "a zero sum of the coefficients and an even sum for the j5 and j6 coefficients
due to symmetry conditions"; "a numerical search for resonances can be performed in a tractable way using the period
ratio P'/P as an initial guess" (p.22). (c) Chaos indicators double as periodic-orbit finders (Gozdziewski et al. 2013,
Migaszewski et al. 2012, p.22). (d) Periodic orbit existence in the CR3BP: Poincare's method of continuation from mu = 0,
with his three kinds (first kind from circular two-body orbits, second kind from elliptical, third kind inclined) (6.5,
p.32); Henon 1965, 1974 as the systematic classification; Sinclair 1970 for families near commensurability; Henrard 1970,
2002 for infinitely many non-symmetric families at L4; "All orbits discovered in the above studies were symmetric
periodic orbits, which are much more easy to generate than non-symmetric periodic orbits" (p.33).

Where it applies (INFERRED): (a) and (b) give the cycler project a cheap prefilter for resonant seeds in the restricted
and elliptic problems: the order-(k, l) location from the period ratio, with the d'Alembert rule as a test that a
proposed resonance can occur at the given order in the perturbation expansion. This is background the project already
uses in the form of synodic-period ratios. The point most useful to the project is the statement that symmetric orbits are
the easy ones to find: the project's searches (symmetric closures) are by that statement biased away from asymmetric
periodic orbits, which is a recorded limitation of any "no X found" result and belongs in the negative-results registry
(project memory: empty regions are conditional on method).

### 2.4 The continuation theory behind the project's seed-and-continue lanes (section 3.2, pp.12-14). New citations.

READ, p.13-14: "Poincare established that there were periodic orbits for all sufficiently small values of mu", with the
Hessian (non-degeneracy) conditions Det(d2 Psi / d omega1 d omega2) not equal to zero and the same for F0, and with the
frequency detuning eta_i' = eta_i (1 + epsilon) giving period T' = T/(1 + epsilon). "Hadjidemetriou [1975a] ... proved that
any symmetric periodic orbit of the CR3BP could be continued analytically to a periodic orbit in the planar general
three-body problem" (p.13); the families were built by Bozis and Hadjidemetriou 1976 and checked for stability by
Hadjidemetriou and Christides 1975; Katopodis 1979 extended this to three dimensions; and Markellos 1981 showed the
resulting general-problem orbits "form continuous mono-parametric families for given masses" (p.13).

Where it applies (INFERRED): this is the existence theory for the step the project takes by hand in `#890` and in the
two-moon periodic lane (`search/two_moon_periodic_890.py`): continue a symmetric periodic orbit of the restricted problem
to one in which the small body has finite mass. Hadjidemetriou 1975a is the primary source that says the symmetric
restricted orbit does continue, and that the continued orbit is again symmetric; it is not held. The synthesis inventory
cites Casoliva et al. for mass continuation (T6) but not this theorem. It also states the failure mode precisely: the
continuation works from a non-degenerate orbit, so a degenerate (bifurcation) generating orbit needs the second-species
treatment (2.5).

### 2.5 Second-species orbits (section 7.2, pp.37-38; 6.6, p.34). Pointers only; sources already held.

READ, p.37: "Periodic orbits showing chains of collision orbits of the original Kepler problem were first studied by
Poincare [1892], who named them the second species solutions. Searches for such orbits in the 3D ER3BP were performed by
Gomez and Olle [1991a], Gomez and Olle [1991b], Bertotti [1991], Bolotin and MacKay [2000], Bolotin [2005] and Palacian and
Yanguas [2006]. In Bertotti's work, a variational method designed to find such periodic orbits was developed. Bolotin
[2005] considered the planar ER3BP with a small mass ratio and eccentricity. Then he proved the existence of many second
species periodic orbits." "The fact that there are ejection-collision orbits in the 3D ER3BP was demonstrated by Llibre and
Pinol [1990]" (p.38). Circular Hill problem (6.6, p.34): "Applications of the circular Hill problem to the Sun-Jupiter-
asteroid system was considered by Bolotin and MacKay [2000], who proved the existence of an infinite number of periodic
and chaotic (almost collision) orbits".

Where it applies (INFERRED): Gomez and Olle 1991 I and II, Bolotin and MacKay 2000 and 2006 are held and digested. Not
held: Bertotti 1991 (variational method), Bolotin 2005 (planar elliptic problem), Palacian and Yanguas 2006, Llibre and
Pinol 1990. Bolotin 2005 is the one that bears on the project's elliptic work: it is the elliptic-problem existence
statement for second-species orbits, and the project's elliptic lane (`core.er3bp`) is where an eccentric primary pair
would test it. Bertotti's variational method is the same family as the discrete-action method of synthesis item T4
(Bolotin and MacKay), so it is a second source for that technique rather than a new one.

### 2.6 The elliptic restricted problem (section 7, pp.36-39). Mostly background; two points.

READ: there is no Jacobi constant in the ER3BP (p.36). Contopoulos 1967 found two integrals for the planar problem in
rotating coordinates with the x-axis through the primaries, and Sarris 1982 three integrals in 3D for small primary
eccentricity and small distance from one primary, which "depended periodically on time" (p.37). Five libration points exist,
three collinear and two triangular, which "pulsate together with their own coordinate systems" (Szebehely 1967) and whose
coordinates "depend explicitly on time" (Choudhry 1977) (p.38). Lie-series integration in rotating-pulsating coordinates
with the true anomaly as independent variable (Delva 1984, p.38); elliptic Hill problem, whose first periodic families were
"all of them highly unstable" (Ichtiaroglou 1981) and which Voyatzis, Gkolias and Varvoglis 2012 mapped with stability for
a satellite around a planet (p.38). Hadjidemetriou 1992 and Haghighipour et al. 2003 found stable 3:1 and 1:2 planar
families (p.38).

Where it applies (INFERRED): (a) The absence of a Jacobi integral is the reason the project's ER3BP lane has no energy
gate and must use a stroboscopic or true-anomaly map; the review states it but adds no method. (b) Delva's Lie-series
integrator is a high-order alternative integrator for the elliptic problem, comparable to the Taylor-method package the
review cites (Jorba and Zou 2005, p.21), which is the most useful integrator pointer for the project's close-passage
work: a variable-order Taylor integrator with automatic order and step selection is the standard tool for the
bicircular and elliptic problems in the Barcelona school whose papers the project already holds. (c) The elliptic Hill
problem is the natural model for a small body near a planet on an eccentric orbit (a moon-tour limit) and the Voyatzis
et al. 2012 families with stability are a ready control population; not held.

### 2.7 Earth-Moon and Sun-Earth-Moon applications (sections 8.1, 8.2, pp.39-40). Little content.

READ, p.39: "by decoupling the Jovian moon n-body system into several three-body systems, it is possible to design a
spacecraft orbit which follows a prescribed itinerary in its visit to Jupiter and its moons. In fact, the Jupiter-
Ganymede-Europa-spacecraft system is approximated as two coupled planar three-body systems [Gomez et al., 2001]."
"For a spacecraft in transit from Earth to the outer planets, its flight may involve a gravity assist, which yields a
four-body problem; once far away from the Earth, it reduces to a three-body problem" (p.39). The Sun-Earth-Moon section
(8.2, p.40) is the Hill-stability criterion (limit about 3.2 Hill radii, Gladman 1993; limits of the criterion, Cuntz and
Yeager 2009) and says nothing on Earth-Moon cyclers or on the bicircular and quasi-bicircular models. The review does not
cite the Barcelona-school bicircular or quasi-bicircular papers, Koon-Lo-Marsden-Ross, or any cycler paper.

Where it applies (INFERRED): nothing new for Blocker 3 (do Earth-Moon families survive the Sun and eccentricity). The one
transferable item is the Hill-radius stability estimate (3.2 R_H, with R_H = (m/3M)^(1/3) a) as a screening rule for which
moon-system flybys sit inside the region where a three-body treatment holds, with the review's own caveat that the
criterion has limits (p.40). Gomez et al. 2001 (the coupled three-body method) is not held under that citation; the
project holds the 2004 spatial-problem paper by the same group (`gomez-koon-lo-marsden-masdemont-ross-2004-...`).

### 2.8 Numerical integrators (section 4.1, pp.19-21). Background, two points.

READ: symplectic maps (Wisdom and Holman 1992) rely on "a low occurrence of events where the orbital velocity could
quickly become large"; "When considering simulations where close encounters are highly probable, an adaptive step-size
method is preferred. Otherwise large errors in the total energy can accrue" (p.21). Leapfrog is stable for a constant step
"less than 2/omega" (p.20). Taylor methods with automatic order selection, TAYLOR by Jorba and Zou 2005 (p.21). Hybrid
symplectic integrators that handle close encounters (Chambers 1999; MERCURY) are cited (p.19, reference list).

Where it applies (INFERRED): confirms that fixed-step symplectic integrators are the wrong tool for a flyby-dominated arc,
which the project already knows. The hybrid close-encounter integrator (Chambers 1999) is the only entry that targets a
mixed regime of long Keplerian arcs punctuated by close passages, which is exactly a cycler's structure; the project's
REBOUND-based lane could use a hybrid scheme of this kind (REBOUND ships one, MERCURY-style) for the multi-year
heliocentric arcs where the full shooter is too slow. The project's memory notes a REBOUND variational-particle gotcha for
custom forces; that stays.

## 3. Statements on cyclers, repeated close encounters, second-species orbits and continuation

The words "cycler", "free return" and "gravity assist" occur only in the places quoted here (text search of the full
review, READ; the one gravity-assist mention is the p.39 quotation in 2.7). Verbatim, the statements that bear on the
project:

- Second species, p.37: "Periodic orbits showing chains of collision orbits of the original Kepler problem were first
  studied by Poincare [1892], who named them the second species solutions." (Note the review's wording: Poincare's
  "second species" are chains of collision orbits, so the review defines them in terms of collisions, not of near-passages
  at finite mass.)
- Continuation, p.13: "Poincare established that there were periodic orbits for all sufficiently small values of mu."
- Continuation to finite mass, p.13: Hadjidemetriou "proved that any symmetric periodic orbit of the CR3BP could be
  continued analytically to a periodic orbit in the planar general three-body problem."
- Families rather than isolated orbits, p.13: the general-problem periodic orbits "form continuous mono-parametric families
  for given masses of the three bodies" (Markellos 1981).
- Existence limit, p.33: "the existence of periodic orbits near the triangular libration points is limited to the mass
  ratios 0 <= mu < mu0, where mu0 is Routh's critical mass ratio" (the one statement in the review of an explicit mass-ratio
  bound on an orbit family).
- Symmetric bias, p.33: "All orbits discovered in the above studies were symmetric periodic orbits, which are much more
  easy to generate than non-symmetric periodic orbits."
- Close triple encounters, p.17: "after a close triple encounter one of the bodies of the three-body problem would
  escape from the system" (McGehee 1974 in the collinear case; Waldvogel 1976 planar), with sufficient escape conditions
  from Standish 1971 and Anosova and Orlov 1994, and the conditions for return of the ejected body from Standish 1972.
  This is a statement about triple collisions, not the finite-mass close passages of cyclers; its relevance is to the
  four-body (two massive moons) lane, where a close three-body encounter is the failure mode (INFERRED).
- Multi-flyby tours, p.39: the Jovian system "is approximated as two coupled planar three-body systems [Gomez et al.,
  2001]."
- Chaos from homoclinic structure, p.35: Poincare "realized that there could be an infinite number of them leading to a
  homoclinic tangle, which appeared to be the first mathematical manifestation of chaos in the CR3BP."

## 4. Primary sources the review cites that the project should obtain

Checked against `docs/notes/CORPUS_INDEX.md` by text search on the authors' names and subjects. "Held" means a file
for that exact paper is in the index.

| Source | Why | Held |
|---|---|---|
| P. M. Cincotta and C. Simo, "Simple tools to study global dynamics in non-axisymmetric galactic potentials - I", A&A Suppl. 147:205-228 (2000); and "Conditional entropy", Celest. Mech. Dyn. Astron. 73:195-209 (1999) | defines MEGNO and its asymptotic values; the source needed before any indicator is implemented | not held |
| C. Froeschle, E. Lega, R. Gonczi, "Fast Lyapunov indicators. Application to asteroidal motion", Celest. Mech. Dyn. Astron. 67:41-62 (1997) | defines FLI, its normalisation and the regular versus chaotic slope criterion | not held |
| R. Gonczi and C. Froeschle, "The Lyapunov characteristic exponents as indicators of stochasticity in the restricted three-body problem", Celest. Mech. 25:271-280 (1981) | Lyapunov exponents in the restricted problem, the nearest analogue of the project's case | not held |
| G. Benettin, L. Galgani, A. Giorgilli, J.-M. Strelcyn, "Lyapunov characteristic exponents for smooth dynamical systems and for Hamiltonian systems; a method for computing all of them", Meccanica 15:9-30 (1980) | the standard computation (tangent-space integration with renormalisation) | not held |
| J. D. Hadjidemetriou, "The continuation of periodic orbits from the restricted to the general three-body problem", Celest. Mech. 12:155-174 (1975); and "The stability of periodic orbits in the three-body problem", Celest. Mech. 12:255-276 (1975) | the continuation theorem and stability treatment for finite-mass moons (`#890` lane) | not held |
| S. Bolotin, "Second species periodic orbits of the elliptic 3 body problem", Celest. Mech. Dyn. Astron. 93:343-371 (2005) | second-species existence in the elliptic problem | not held |
| M. L. Bertotti, "Periodic solutions for the elliptic planar restricted three-body problem: a variational approach", in Predictability, Stability and Chaos in N-Body Dynamical Systems, NATO ASI Series, pp.467-473 (1991) | variational search for second-species orbits | not held |
| D. C. Heggie, "Regularization using a time-transformation only", in Recent Advances in Dynamical Astronomy, Astrophys. Space Sci. Libr. 39:34 (1973); and "Redundant variables for 'global' regularization of the three-body problem", Celest. Mech. 14:69-71 (1976); C. A. Burdet, "Regularization of the two body problem", ZAMP 18:434-438 (1967) | the time-transformation regularisation, the lighter alternative to coordinate regularisation for a numerical integrator | not held |
| P. Kustaanheimo and E. Stiefel, "Perturbation theory of Keplerian motion based on spinor regularization", J. reine angew. Math. 218:204 (1965); T. Levi-Civita, "Traiettorie singolari ed urti nel problema ristretto dei tre corpi", Ann. Mat. Pura Appl. (3) 9:1-32 (1903) | the two standard regularisations; the review cites them without writing them down | not held under these entries; Levi-Civita regularisation is used in the held Font-Nunes-Simo and Marco-Niederman papers |
| G. Gomez, W. S. Koon, M. W. Lo, J. E. Marsden, J. Masdemont, S. D. Ross, "Invariant manifolds, the spatial three-body problem and space mission design", AAS/AIAA Astrodynamics Specialist Meeting, Proc. Vol. 109, pp.3-22 (2001) | the coupled-three-body Jovian tour design cited at p.39 | not held (the same group's 2004 Nonlinearity paper is held) |
| L. V. Vela-Arevalo and J. E. Marsden, "Time-frequency analysis of the restricted three-body problem: transport and resonance transitions", Class. Quantum Grav. 21:351-375 (2004) | a frequency-analysis regular-versus-chaotic diagnostic for the restricted problem | not held |
| J. E. Chambers, "A hybrid symplectic integrator that permits close encounters between massive bodies", MNRAS 304:793-799 (1999) | the mixed Keplerian-plus-close-encounter integrator | not held (title read in the review's reference list; the volume and pages are from my memory and must be checked) |
| A. Jorba and M. Zou, "A software package for the numerical integration of ODEs by means of high-order Taylor methods", Exp. Math. 14(1):99-117 (2005) | the Taylor integrator the Barcelona-school papers in the corpus use | not held |
| J. Henrard and J. F. Navarro, "Families of periodic orbits emanating from homoclinic orbits in the restricted problem of three bodies", Celest. Mech. Dyn. Astron. 89:285-304 (2004) | a fast multiple-Poincare-section method for quasi-periodic orbits, "robust near chaotic regions" (p.33) | not held |
| M. Henon, "Exploration numerique du probleme restreint", Ann. d'Astrophysique 28:499-511 and 992-1007 (1965); and the 1974 follow-up, Celest. Mech. 10:375-388 | the systematic family classification | not held (the project's Henon material is by citation in other digests) |
| Delva 1984; Palacian and Yanguas 2006; Voyatzis, Gkolias, Varvoglis 2012 (CMDA 113:125-139); Llibre and Pinol 1990 | elliptic-problem integrator, second-species and elliptic-Hill sources | not held |

Already held and cited by the review, no action: Gomez and Olle 1991 parts I and II; Bolotin and MacKay 2000 and 2006;
Font, Nunes and Simo 2002 (and 2009).

The two sources to obtain first are Cincotta and Simo 2000 and Hadjidemetriou 1975a: the first is a precondition for any
indicator work, the second is the continuation theorem the `#890` lane leans on.

## 5. Verdict

A short, mostly historical survey that adds, for this project, one method it lacks (the chaos indicators, with citations but
without the thresholds needed to use them, and with the caveat that a hyperbolic cycler is expected to show positive
exponents) and a good reference list; it contains nothing on cyclers or flybys, and its collision-regularisation and
Earth-Moon sections are too thin to act on.
