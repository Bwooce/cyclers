# Frontier scoping beyond #281–#287 — ER3BP, BCR4BP, 3D, QP-tori, epoch-aware

**Date:** 2026-06-16
**Frame:** beyond the current discovery envelope. The current envelope is the
planar, autonomous, mass-parameter-only CR3BP plus its prioritizer stack
(`five_tier_prioritizer.py` + neural prefilter). Campaigns #283-#287
(`docs/notes/2026-06-16-281-discovery-campaign-roadmap.md`) all sit inside that
envelope; the roadmap planner explicitly **rejected** BCR4BP, low-thrust, and
halo-NRHO heteroclinic-search as "multi-week Track-A builds." The user has now
asked: **what's the frontier BEYOND #281-#287, and what would it cost?**

This note scopes five candidate axes. Per-axis 7-field template (definition,
existence prior, representable in current code, build cost, first deliverable
IC, risk, what it unlocks). Rank-ordered summary table at the end. **No
catalogue change, no source code change.** Planning doc only.

Per CLAUDE-md rules: catalogue scope is currently *repeating cyclers only*
(`project_catalogue_scope_cyclers_only`). Axis 5 (epoch-aware / quasi-cyclers
of opportunity / precursor MGA tours) is genuinely outside that scope as
stated; the scope tradeoff is surfaced explicitly, not pre-decided.

---

## What's IN-envelope today (do not re-propose)

* **Planar CR3BP propagator** (`src/cyclerfinder/core/cr3bp.py:73 cr3bp_eom`).
  Already 6D state `(x, y, z, vx, vy, vz)`; the EOM has full z-coupling
  (`az = -(1-mu) z / r1^3 - mu z / r2^3`). "Planar" today is a *corrector
  convention*, not a propagator limitation.
* **Symmetric NRHO/tulip corrector**
  (`src/cyclerfinder/search/nrho_continuation.py:178 correct_symmetric_nrho`).
  Free vars `(z0, ydot0, T)` at the fixed-x0 perpendicular x-z-plane crossing
  IC `(x0, 0, z0, 0, ydot0, 0)`. **Already 3D**: residual is
  `(y, xdot, zdot)` at the half-period perpendicular crossing.
* **Asymmetric planar corrector**
  (`src/cyclerfinder/search/cr3bp_general_periodic.py:329
  correct_general_periodic`). Drops xdot0=0; planar (z=0).
* **Tulip multi-shooting**
  (`src/cyclerfinder/genome/multi_shooting.py multi_shoot_periodic`).
  Sundman-regularised; tulip-symmetric IC; arbitrary `n_segments`.
* **Repeated-moon Lambert genome**
  (`src/cyclerfinder/search/moon_cycler_genome.py`).
* **Real-ephemeris ICs**
  (`src/cyclerfinder/core/ephemeris.py Ephemeris("astropy")`,
  JPL DE440 via astropy). Today this is a *validation oracle* (V3); the V4 GMAT
  lane (`reference_gmat_install`) is the high-fidelity check.
* **5-tier prioritizer stack** (Braik-Ross + Zhou-Armellin + resonance-network
  + FTLE + lobe-overlap, composed by
  `src/cyclerfinder/search/two_tier_prioritizer.py`). EM-VU hardcoded.
* **Binary-star search** (`src/cyclerfinder/search/binary_star_search.py:143
  figure_seeded_search`), reaches mu up to 0.5.

What is NOT in-envelope: time-periodic Hamiltonians (BCR4BP/QBCP);
non-circular primaries (ER3BP); 3D family-coverage (vertical Lyapunov, halo,
inclined; the *propagator* is 3D — the *family search* is not pointed at z!=0);
quasi-periodic invariant 2-tori; real-ephemeris discovery (epoch as a free
variable).

---

## Axis 1 — ER3BP (Elliptic Restricted 3-Body)

### (a) Definition
Drop the assumption that the primaries are on circular orbits — they orbit
their barycentre on a Keplerian ellipse of eccentricity `e`. The Hamiltonian
becomes explicitly time-periodic (period = 2π in true anomaly) and the Jacobi
integral is replaced by a *pulsating* Jacobi-like integral that conserves only
on average. Definitive treatment: Szebehely (1967, *Theory of Orbits in the
Restricted Problem of Three Bodies*) §10.4 ("The Elliptic Restricted Problem"),
which is in our corpus (`szebehely-1967-theory-of-orbits...`).

### (b) Existence prior
**MEDIUM-HIGH**. ER3BP periodic-orbit literature is well-established:
* Broucke (1969), "Stability of Periodic Orbits in the Elliptic, Restricted
  Three-Body Problem" — *AIAA Journal* 7:1003 — continues CR3BP planar
  symmetric families into `e>0`.
* Hou & Liu (2011), "Continuation of periodic orbits in the elliptic
  restricted three-body problem" — gives systematic Lindstedt-Poincaré
  continuation from CR3BP roots into `e>0`.
* Peng & Xu (2015), "Stability of two groups of multi-revolution elliptic
  halo orbits in the elliptic restricted three-body problem" — *CeMDA* 123:279,
  publishes ICs for ER3BP halos at Earth-Moon `e=0.0549`.

**Cycler-class specifically:** thin. The ER3BP literature is mostly libration-
point families and DROs, not Aldrin-style flyby-cyclers. The Hou-Liu
continuation framework would carry CR3BP cyclers into `e>0` as a *systematic*
extension; we have no published evidence of ER3BP-only cycler families that
don't continue back to CR3BP. The honest answer to "are there ER3BP cyclers
that DO NOT continue to CR3BP as e→0?" is **unknown** — no published
counter-example, but the question has not been systematically answered in the
cycler corpus.

### (c) Representable in current code?
**NO.** The CR3BP EOM
(`src/cyclerfinder/core/cr3bp.py:73 cr3bp_eom`) hard-codes circular primaries
through the rotating-frame Coriolis terms (`+2 vy`, `-2 vx`) and the constant
primary positions at `(-mu, 0, 0)` and `(1-mu, 0, 0)`. The ER3BP requires
either (i) a pulsating-rotating frame with explicit true-anomaly dependence in
the EOM coefficients, or (ii) inertial integration with two primaries on a
prescribed ellipse. The corrector
(`cr3bp_general_periodic.py:329 correct_general_periodic`) uses the
y=0 return map and Jacobi-from-`x0,xdot0` algebra; both assumptions break in
the ER3BP (Jacobi is not conserved; the IC must include the primaries' true
anomaly).

### (d) Build cost
* Pulsating-rotating EOM module (`er3bp.py`, ~250 LOC) — vector field with `f`
  (true anomaly) as a free parameter; pulsating-Jacobi diagnostic.
* True-anomaly-of-IC parametrisation in the corrector (the second free phase
  alongside x0); ~150 LOC patch to `cr3bp_general_periodic.py` *or* a new
  `er3bp_periodic.py` (preferred — keeps the CR3BP path clean).
* Family continuation in `e` (natural-parameter starting from a CR3BP root at
  `e=0` per Hou-Liu); ~150 LOC.
* Closure test extension (period closes in true-anomaly units; period
  multiplier `k` against `T_primary`); ~50 LOC.
* Scorer adaptation: Jacobi-based scorers (`reachable_network`, FTLE) need a
  *pulsating-Jacobi* substitute or to fall back to V∞ only; ~100 LOC patch.
* Tests + goldens against Hou-Liu and Peng-Xu published ICs; ~200 LOC.

**Total: ~900 LOC, 2-3 weeks of focused work.** Continuation is the load-
bearing piece — the rest is plumbing.

### (e) First deliverable IC
**ER3BP Earth-Moon halo from Peng-Xu (2015) Table 1** (need to digitise; not
in our PDF corpus, but referenced in the open-access *CeMDA* abstract).
Substitute: continue the Earth-Moon L2 Southern 9:2 NRHO (already in
`tulip.py:KOBLICK_2023_TABLE4` at Np=1) in `e` from 0 to 0.0549 via Hou-Liu
natural-parameter; use the Np=1 NRHO IC `(x0, 0, z0, 0, ydot0, 0)` at `e=0` as
the continuation seed with `f0=0` (perihelion). This is the **e=0 root**, not
the ER3BP novel IC — but the existence of the continuation *to* `e=0.0549` is
itself the deliverable for the first build.

For a *novel* search seed: Earth-Mars Aldrin S1L1 IC at `(e_E=0.0167,
e_M=0.0934)` — but that lives in the multi-body solar-system case, which is
strictly beyond ER3BP (ER3BP is two primaries, not Sun-Earth-Mars). The
correct first IC is cislunar.

### (f) Risk
* **Cyclers don't appear ER3BP-natively.** ER3BP literature is dominated by
  libration families; no published ER3BP-only cycler exists today. The most
  likely outcome of a cycler-search in ER3BP is "every member continues back to
  a CR3BP member" — a *refined* catalogue, not a new species.
* **Pulsating-Jacobi scorer regression.** The 5-tier prioritizer is built on
  conserved-Jacobi. Adapting it costs build time with no guarantee the new
  scorers retain their cislunar validation.
* **Small `e` for Earth-Moon (0.0549) means small perturbation** — most CR3BP
  families survive with O(e) corrections. The genuinely interesting `e` for
  cyclers is Jupiter-Europa (`e=0.0094`, even smaller) or Pluto-Charon
  (`e≈0`, no effect). Mars-Phobos (`e=0.0151`) is too tight a system for
  cyclers. There is no body pair with both a published cycler context AND
  large `e`.

### (g) What it unlocks
* Refined family curves at real `e` for cislunar discovery — likely closer to
  the real-ephemeris (V3) results than CR3BP, but strictly inferior to V3.
* The genuinely new species — ER3BP-only families that don't continue back —
  would be a real discovery, but the existence prior is weak. This is the
  "no-prior speculative type" tier; **must be backed by an explicit dynamical
  hypothesis** to justify the build.

---

## Axis 2 — BCR4BP / QBCP (Bicircular Restricted 4-Body)

### (a) Definition
Quasi-Bicircular Problem: the Earth-Moon CR3BP made explicitly periodically
time-dependent through the Sun, where the Sun-Earth-Moon primaries follow a
*coherent* (genuine, periodic, 3-body) trajectory close to bicircular.
Definitive treatment: Andreu (1998), *The Quasi-Bicircular Problem*, PhD
thesis, Univ. de Barcelona, advisor C. Simó. Modern open-access
restatements: Gimeno-Jorba (2018) *Frontiers in Applied Mathematics and
Statistics* 4:32; Rosales-Jorba-Jorba (2023) *CeMDA* 135:15. All three are
addressed in `docs/notes/2026-06-14-andreu-quasi-bicircular-digest.md` which
ALREADY digests the Hamiltonian, α-Fourier parameters, and POL1/POL2 dynamical-
substitute goldens.

### (b) Existence prior
**HIGH for refined libration families; LOW for Aldrin-style cycler families.**
* Rosales-Jorba-Jorba 2023 *CeMDA* 135:15 — explicit QBCP invariant manifolds
  near L1/L2 with full ICs (Table 4) — POL1 `x=-0.8369141677649317, py=
  -0.8391311559808445`; POL2 `x=-1.1556836078332600, py=-1.1587306159501061`
  (digested already).
* Gimeno-Jorba 2018 — same model, parameter table to order k=13.
* IAGSP (intermediate auxiliary geometric Sun perturbation) family of Howell et
  al. — referenced by Rosales-Jorba as a Sun-coherent libration-family
  extension. No sourced cycler ICs.
* **No Sun-Earth-Moon QBCP CYCLER catalogue exists in the open literature
  per the Andreu digest.** Andreu's payoff is in transfer-orbit structure
  (Sun-synodic-resonant transfers; the "fast periodic transfer orbits extended
  from EM-RTBP to the Sun-Earth-Moon QBCP" paper, CeMDA 2008), not cyclers.

### (c) Representable in current code?
**NO.** The QBCP needs the α_i(ϑ) Fourier-series Hamiltonian as a new module;
the CR3BP propagator has no time-periodic forcing. No existing file is
extensible — this is a new `qbcp.py` analogous to `cr3bp.py`.

### (d) Build cost (smallest viable BCR4BP genome)
The Andreu digest already lists the cost; this scoping just confirms the
breakdown:
* α_i(ϑ) Fourier coefficient table (Gimeno-Jorba 2018 Table 4, k=13, 8
  functions — ~104 sourced floats), typed once: ~150 LOC.
* QBCP propagator: H_QBCP Hamilton's equations + α_i(ϑ) evaluator:
  ~200-300 LOC.
* Symmetric-corrector extension (the dynamical-substitute periodic orbits
  cross y=0 at ϑ=0; the corrector becomes 4 free vars including ϑ_0): ~200
  LOC.
* Closure test (T = T_s = 2π/ω_s = 6.7912 EM TU — period of the Sun's synodic
  forcing): ~50 LOC.
* Scorer adaptation: NONE. Defer prioritizer to CR3BP-screened candidates
  THEN refine in QBCP — keeps the build narrow.
* Tests + POL1/POL2 reproduce gate: ~150 LOC.

**Total smallest viable: ~800 LOC, 2-3 weeks.** That meets the user's
"smallest possible BCR4BP genome (propagator + symmetric corrector + closure
test) that we could build in <2 weeks" question — the answer is **just barely
feasible at 2 weeks** for the model + POL1/POL2 reproduction; cycler-class
discovery on top is another build.

### (e) First deliverable IC
**POL1 substitute (Rosales-Jorba 2023 Table 4):**
```
state6 = (-0.8369141677649317, 0.0, 0.0, 0.0, -0.8391311559808445, 0.0)
T_target = 2π / ω_s = 2π / 0.9251959855... ≈ 6.7912 EM TU
ϑ_0 = 0  (Sun phase at perihelion)
mu = 0.012150581600000, m_S = 328900.5423094043, a_S = 388.8111430233511
```
This is THE sourced QBCP IC; the goldens are digested in the Andreu note.
Integrating the QBCP vector field one period from POL1 and recovering closure
is the build's first deliverable. **NO published QBCP cycler IC exists** — so
the first cycler-class deliverable is the second build, not the first.

### (f) Risk
* **No cycler families published.** This is the dominant risk: QBCP is most
  productive for transfer orbits and libration-family substitutes, neither of
  which is a cycler. The build pays off if-and-only-if synodically-locked
  Aldrin-class cyclers exist that have no CR3BP analogue — and this is an
  open question (the Andreu digest §3 marks it "narrow watch-item").
* **V3 (DE440 ephemeris) dominates fidelity.** The Andreu digest already
  documents: QBCP sits strictly below V3, so as a validation rung it is
  redundant. Its value is methodological (deterministic, paper-reproducible).
* **Sun phase ϑ_0 is a fourth free variable.** This breaks the (k1, k2)
  winding-class taxonomy at the corrector level; the dedupe pipeline needs
  re-thought.

### (g) What it unlocks
* Cislunar cyclers whose existence is *gated* by the ~29.5-d Sun-synodic
  forcing (resonant with T_s). If discovered, this is a genuine new class.
* POL1/POL2 substitute orbits as sourced regression-rung between CR3BP and
  V3 — useful when a real-ephemeris cislunar discovery's CR3BP origin is
  opaque.
* The dynamical bridge between cislunar (Earth-Moon CR3BP) and
  interplanetary (Sun-influenced) trajectory designs.

---

## Axis 3 — 3D / broken-plane cyclers

### (a) Definition
Lift the z=0, ż=0 restriction *in the discovery search* — search for periodic
orbits with non-trivial out-of-plane motion. CR3BP families: vertical Lyapunov
orbits (planar Lyapunov bifurcation in z), halo orbits (planar/3D bifurcation),
NRHOs (halo continuation to near-rectilinear regime), inclined repeated-
encounter cyclers (3D analogues of the Liang CGE class). Definitive treatment:
Howell (1984), "Three-dimensional, periodic, 'halo' orbits," *CeMDA* 32:53;
Doedel-Romanov-Paffenroth-Keller-Dichmann-Galán-Vanderbauwhede (2007),
"Elemental periodic orbits associated with the libration points in the
circular restricted 3-body problem," *Int. J. Bifurcation and Chaos* 17:2625
(comprehensive 3D family atlas).

### (b) Existence prior
**HIGH for cislunar libration-family extensions; MEDIUM for 3D cyclers
specifically.**
* Howell 1984 (above) — original halo.
* Folta-Bosanac-Cox-Howell (2017), "Trajectory design tools for libration and
  cislunar environments," *CeMDA* 129:283 — current standard reference.
* Ross-Roberts-Tsoukkas 2025 AAS-25-621 (in our corpus,
  `ross-roberts-tsoukkas-2025-stable-ballistic-earth-moon-cyclers-AAS-25-621.pdf`)
  — already covered 5 stable Earth-Moon cycler families. Per memory
  `feedback_published_rounded_values_are_display`, the C21 family has
  asymmetric (xdot0!=0) members at C=3.129389531...; the 3D extension to
  halo/NRHO-rooted cyclers is partially in their analyses but not exhaustively
  searched.
* Roberts-Tsoukkas-Ross 2026 journal version (in our corpus,
  `roberts-tsoukkas-ross-2026-stable-prograde-em-cyclers-journal.pdf`) —
  binary-mu prograde stable cyclers, partially 3D.

**Cycler-specific 3D ICs**: a small handful are published. The Liang 2024 JGCD
Callisto-Ganymede-Europa triple cyclers (in corpus,
`liang-2024-callisto-ganymede-europa-triple-cyclers-JGCD.pdf`) are *planar
multi-rev*; the 3D extension is open. Howell-Marchand-Lo (2001) "Temporary
satellite capture of short-period Jupiter family comets" gives 3D heteroclinic
chains — not cyclers, but the methodology transfers.

### (c) Representable in current code?
**MOSTLY YES — much smaller build than ER3BP/BCR4BP.** Critical finding from
reading `cr3bp.py`: the propagator `cr3bp_eom` (line 73) is *already* 6D and
has full z-coupling. The variational STM `cr3bp_stm_eom` (line 100) computes
all 9 partials including `uzz, uxz, uyz`. The NRHO corrector
(`nrho_continuation.py:178 correct_symmetric_nrho`) is *already 3D* — free
vars `(z0, ydot0, T)`, residual `(y, xdot, zdot)` at the half-period x-z-
plane crossing.

**The constraint is family-coverage, NOT propagator dimensionality.** The
`correct_symmetric_nrho` corrector handles the Howell-symmetric x-z plane
crossing IC. The *asymmetric* corrector
(`cr3bp_general_periodic.py:329 correct_general_periodic`) is *planar*
(z=0 hardcoded) — extending it to 3D asymmetric is ~150 LOC of generalisation.
Discovery drivers currently call symmetric/planar correctors; pointing them at
z!=0 ICs is mostly driver work, not solver work.

### (d) Build cost
* Asymmetric 3D corrector (extend `correct_general_periodic` to 6D state, drop
  z=0): ~150 LOC.
* Vertical Lyapunov + halo continuation drivers from L1/L2 (modify
  `nrho_continuation.continue_nrho_family` to walk along z0 and ydot0 rather
  than just x0): ~200 LOC.
* Halo-bifurcation detector (existing `bifurcation_detector.py` is 2D Floquet
  — extend to 3D unit-circle eigenvalue tracking): ~150 LOC.
* 3D-aware family-switching corrector: ~100 LOC.
* Scorer adaptation: minimal — most scorers operate on Jacobi + V∞, which are
  3D-natural. The 2D winding classifier
  (`binary_star_search.winding_topology`) needs a 3D analogue (vertical
  rotation count); ~100 LOC.
* Tests + goldens against Howell 1984 + Folta-Bosanac-Cox-Howell halo ICs:
  ~150 LOC.

**Total: ~850 LOC, 1.5-2 weeks. Lowest-friction frontier axis by far.** The
propagator and tools are already 3D-capable; the build is mostly converting
drivers and family-coverage.

### (e) First deliverable IC
**Earth-Moon L1 Northern Halo Az=8000 km (Folta-Bosanac-Cox-Howell 2017
standard archive value; not in our PDF corpus, but the JPL three-body periodic-
orbit catalog
https://ssd.jpl.nasa.gov/tools/periodic_orbits.html publishes the IC):**
```
state6 ≈ (0.8233, 0.0, 0.0181, 0.0, 0.1263, 0.0)
T ≈ 2.7430 TU  (Earth-Moon)
mu = 1.215058560962404e-2
```
*Verified path:* JPL's catalog is the routine reference; an exact-digit
sourcing requires either (i) one query of the JPL three-body catalog API or
(ii) reading the L1 halo IC from Folta-Bosanac-Cox-Howell 2017 §4. Either
satisfies the orbit-closure-discipline rule.

For a *novel* cycler: continue an Earth-Moon Ross-RT (1,1) cycler (already
reproduced at `cr3bp_general_periodic` per #249) into z0!=0 via a 3D
asymmetric corrector. The (1,1) family's planar root is sourced; the 3D
extension is open.

### (f) Risk
* **Halo families are EXHAUSTIVELY MAPPED.** Howell 1984 + the JPL catalog +
  Folta-Bosanac-Cox-Howell 2017 cover the Earth-Moon libration halos
  comprehensively. The discovery payoff is NOT in halos themselves, but in
  3D Aldrin-class cyclers (Earth-flyby chains lifted out of plane) — a
  smaller, less-populated literature.
* **Out-of-plane Aldrin cyclers have a thin existence prior.** Patel 2019
  (corpus, `patel-2019-earth-mars-cycler-vehicle-conceptual-design`) treats
  vehicle-design at planar Aldrin cyclers; no 3D Earth-Mars cycler family
  has, to our knowledge, been published. The search-space is large and the
  prior is thin.
* **3D bifurcation detection is messier than 2D.** Floquet multipliers
  bifurcate in 6 (vs 4) dimensions; the family-switching corrector needs more
  care.

### (g) What it unlocks
* Vertical-Lyapunov and halo families at every CR3BP body system the
  catalogue already touches (Earth-Moon, Mars-moons, Saturnian moons,
  Jovian moons, Pluto-Charon). Each new system inherits the same 3D
  family-coverage.
* 3D Aldrin-class cyclers (Earth-Mars, Earth-Venus) — would be a genuine
  novel-class discovery if any exists.
* Halo-rooted heteroclinic-network cyclers — connects to the #267
  resonance-network scorer's natural targets (the resonant manifold tier was
  built as a scorer; converting to a *search driver* in 3D is one path to
  novel-cycler discovery the planner already deferred).

---

## Axis 4 — Quasi-periodic invariant 2-tori

### (a) Definition
Drop strict periodicity. Search for *invariant 2-tori* — trajectories that
densely fill a 2-dimensional torus in phase space, parameterised by two
incommensurate frequencies (longitudinal and transverse). A trajectory on
such a torus *never closes* but stays bounded and structured. Closure
condition is replaced by *invariance*: KAM-style, the torus is invariant
under the flow, and the trajectory traces a circle on its Poincaré section.
Definitive treatment: Olikara (2016), "Computation of quasi-periodic
invariant tori in the restricted three-body problem," PhD thesis, Purdue,
advisor K. Howell; Olikara-Scheeres (2012), "Numerical method for computing
quasi-periodic orbits and their stability in the restricted three-body
problem," *Adv. Astronaut. Sci.* 145:1085 (the **GMOS** algorithm —
"Generalized Method Of Sections").

### (b) Existence prior
**MEDIUM. Olikara is the canonical source.**
* Olikara 2016 thesis (above) — full QP-torus methodology + ICs at L1/L2/L3
  Earth-Moon halo-extension tori.
* Olikara-Howell (2010), "Quasi-periodic trajectories around the libration
  points in the restricted three-body problem," AAS/AIAA Space Flight
  Mechanics Meeting AAS-10-167 — the original QP-torus continuation.
* Baresi-Scheeres (2017), "Bounded relative motion under stochastic
  perturbations on a torus," *CeMDA* 128:493 — QP-torus stability methodology.

**Cycler-class on QP tori:** essentially none. QP tori in the cycler literature
are mostly stability-analysis structures around periodic cyclers (the orbit
plus its quasi-periodic neighbourhood). The user's premise — "cyclers that
almost close — repeating-up-to-a-rotation — are a real class our closure
filter discards by construction" — is correct in principle: a closed cycler
with one Floquet multiplier at `e^{i theta}` on the unit circle is
*surrounded* by an invariant torus. The QP-torus version would be the
quasi-cycler living on that torus.

### (c) Representable in current code?
**NO.** The closure gate (#259) and corrector residuals require literal state
return. The torus discovery condition is fundamentally different: parameterise
the torus by `(theta_1, theta_2)`, represent it as a truncated Fourier series
`X(theta_1, theta_2)`, and the invariance condition is that the flow advances
`theta_i` linearly in time. Our entire pipeline assumes single-trajectory
closure.

### (d) Build cost
* Fourier-coefficient torus parameterisation:
  ~200 LOC (one body of indexed Fourier modes plus FFT-based stroboscopic
  Poincaré map).
* GMOS (Generalized Method of Sections) algorithm — Olikara-Scheeres 2012 —
  for invariance-residual solve: ~400-500 LOC. This is the load-bearing
  piece; non-trivial nonlinear solver design.
* Closure-test replacement: invariance residual norm rather than
  state-return norm: ~50 LOC.
* Continuation in (longitudinal, transverse) frequencies: ~200 LOC.
* Stability via Floquet-Bloch of the torus map: ~150 LOC.
* Scorer adaptation: the prioritizer assumes periodic orbits; FTLE and
  reachable-set scorers translate but need a torus-aware integration: ~250
  LOC patch.
* Tests + goldens against Olikara 2016 thesis Tables: ~300 LOC.

**Total: ~1500-2000 LOC, 4-6 weeks (the largest of the five axes).** GMOS
is the build that makes-or-breaks the time estimate; an undergraduate-quality
implementation might fit; a publication-quality one is a thesis chapter.

### (e) First deliverable IC
**Olikara-Scheeres (2012) Earth-Moon L1 halo-extension torus, longitudinal
freq `ω_1 = 2.842 rad/TU`, transverse freq `ω_2 = 0.512 rad/TU` (digitised
estimates — source has full digits).** The IC is not a 6-vector but a
Fourier-series — *one circle* on the torus's Poincaré section. The minimum
representation is N_theta1=64 Fourier modes around that circle, each a
6-vector — so the "IC" is 384 floats sourced from Olikara's thesis.

This is qualitatively different from every other axis: there is no
"6-vector first IC." Sourcing the torus from Olikara requires either the
thesis (not in our PDF corpus) or the 2012 paper (open access, available).

### (f) Risk
* **No cycler-class QP-torus prior.** The QP-torus literature is libration-
  family-focused. Cycler-extension is conceptually clean but unexplored.
* **Build cost is real.** GMOS is the hardest algorithm of the five axes.
* **Scope tension with `project_catalogue_scope_cyclers_only`.** A
  QP-cycler does not "repeat" in any literal sense — it lives on a torus.
  Either the catalogue scope expands or QP-cyclers are a separate deliverable
  class.
* **Validation chain.** V0-V5 gauntlet assumes single-trajectory state-return.
  Adapting the gauntlet to torus-aware closure is non-trivial.

### (g) What it unlocks
* The entire QP-orbit universe — Olikara's published catalog of cislunar
  QP tori becomes searchable terrain.
* Cyclers-as-tori, where a planar cycler with one Floquet multiplier on
  the unit circle generates a 1-parameter family of QP cyclers. The known
  unstable Ross-RT Earth-Moon members have unit-modulus Floquet multipliers
  for some Jacobi values — those are the natural seeds.
* Stochastic-perturbation robustness analysis (Baresi-Scheeres methodology)
  on existing cyclers — a different and valuable use case.

---

## Axis 5 — Epoch-locked / quasi-cycler / precursor-MGA genome (user-explicit)

### (a) Definition
Two sub-axes treated as one:

**(5a) Cycler precursor tours.** The MGA chain that *inserts* a spacecraft
into a steady-state cycler, with the launch epoch as a free variable.
Aldrin's Earth-Mars cycler requires both the cycler AND the insertion
trajectory; the insertion is epoch-locked through the real planetary
ephemeris. Without the precursor, the cycler catalogue describes orbits no
spacecraft can actually reach.

**(5b) Quasi-cyclers / cyclers of opportunity.** Orbits that satisfy
cycler-closure for 3-5 returns inside a 10-15 year alignment window of the
real ephemeris, then drift out. Long enough to be useful, structured enough
to search as closure problems (not arbitrary flybys), but not strictly
periodic.

Definitive treatment for (5a): McConaghy (2004), "Design and optimization of
interplanetary spacecraft trajectories," PhD thesis, Purdue, advisor J.
Longuski (in our corpus, `mcconaghy-2004-design-optimization-...`); Ceriotti
(2010), "Global optimisation of multiple gravity assist trajectories," PhD
thesis, Glasgow (in our corpus, `ceriotti-2010-global-optimisation-mga-glasgow-phd`).

Definitive treatment for (5b): Russell-Strange (2009),
"Cycler trajectories in planetary moon systems," *J. Spacecraft and Rockets*
46:439; Petropoulos-Longuski (2000), "Shape-based algorithm for the automated
design of low-thrust, gravity assist trajectories," *JSR* 37 (in spirit —
sub-publication-class repeated tours).

### (b) Existence prior
**STRONG for both sub-axes.**

For (5a) precursor tours (real, mission-validated trajectories):
* Galileo VEEGA — Diehl-Belbruno-Roberts 1986 IAF.
* Cassini VVEJGA — D'Amario-Bright-Wolf 1992 *AAS* 92-110.
* Voyager Grand Tour — Flandro 1966 *Acta Astronautica* 12:329.
* Petropoulos-Longuski 2000 pump tours (above).
* Strange-Russell 2007, "Modeling Tisserand-leveraging cycler trajectories for
  Mars cycling," *AAS-07-126*.
* Vasile-Campagnola 2009 (in our corpus,
  `vasile-campagnola-2009-lowthrust-mga-europa-JBIS-arxiv-1105.1823.pdf`).
* GTOC winners (Vavrina, Strange) — multiple competition publications.
* Heaton-Strange-Longuski (2002), "Optimal low-thrust gravity-assist
  trajectories to Jupiter using nuclear electric propulsion," *J. Guidance,
  Control, and Dynamics* 25:1147.
* Campagnola-Russell 2009 endgame Parts A&B (in our corpus,
  `campagnola-russell-2009-endgame-part*`).

For (5b) quasi-cyclers:
* Hughes-Edelman-Longuski 2014 (in our corpus,
  `hughes-edelman-longuski-2014-fast-mars-free-returns-venus-ga-AIAA-2014-4109.pdf`)
  — fast Mars free returns via Venus GA. These ARE the
  cycler-of-opportunity prototype: structured Mars-return geometry, but
  epoch-locked through Venus phasing.
* Jones-Hernandez-Jesick 2017 (corpus, `jones-...-vem-triple-cyclers-AAS-17-577.pdf`)
  — Venus-Earth-Mars triple cyclers — partially epoch-locked.
* Genova-Aldrin 2015 (corpus, two Earth-Moon-cycler papers) — Earth-Moon
  cyclers with real-ephemeris insertion epochs.
* Tito-MacCallum-Carrico 2013 (corpus,
  `tito-maccallum-carrico-2013-feasibility-manned-mars-free-return-2018-IEEE`)
  — the 2018 Mars free-return mission. Per
  `project_catalogue_scope_cyclers_only`, withheld as validation artifact;
  but this IS the epoch-locked-discovery model.

The user's framing — "orbits we can only find at one point in time, they're
not a static cycler so much as a 10-year chaining of unique opportunities" —
is **precisely the established mission-design class**. The literature is
mature; the gap is that our catalogue is built around the abstracted (static)
cycler, not the deployed (epoch-locked) trajectory.

### (c) Representable in current code?
**PARTIALLY.** The components exist:
* `src/cyclerfinder/core/ephemeris.py Ephemeris("astropy")` — JPL DE440 real
  planetary positions/velocities at an arbitrary epoch. Today this is a
  validation oracle (V3 lane).
* `src/cyclerfinder/core/lambert.py` — multi-rev Lambert solver.
* `src/cyclerfinder/core/flyby.py` — patched-conic flyby.
* `src/cyclerfinder/core/fbs_match_point.py` — finite-burn-shooting match
  point. FBS analytic gradients exist for the optimization stage.
* `src/cyclerfinder/search/free_return_chain.py` — multi-arc chained
  free-return closure (already used for S1L1 / Mars free-return reproduction).
* `src/cyclerfinder/search/multiarc_closure.py` — chained closure residual.

What's MISSING:
* **Epoch as a free variable in the closure-search driver.** Today
  free_return_chain assumes a fixed epoch (V0/V3 testing convention). For
  discovery, the launch epoch `t0` must be a search dimension — sweeping
  10-15 year windows.
* **Closure window relaxation.** A 3-5-return quasi-cycler closes within
  some tolerance for 3-5 returns then drifts; the existing closure gate
  (#259) is strict.
* **MGA combinatorics over real ephemeris.** Tisserand graph at fixed t0,
  flyby-sequence enumeration. We have `tisserand.py` but it's body-pair
  static, not real-ephemeris driven.

This is NOT a from-scratch genome; it is the existing pipeline pointed at a
different question, with epoch promoted from constant to search variable.
**Substantially smaller build than ER3BP/BCR4BP/QP-tori.**

### (d) Build cost
* Epoch-window search driver (sweeps `t0` over a 10-15-yr window, calls the
  existing closure for each `t0`): ~200 LOC.
* Closure-window relaxation (3-5-return closure inside tolerance, allow
  drift outside): ~150 LOC patch to `multiarc_closure.py`.
* Real-ephemeris Tisserand pre-screen: ~150 LOC patch to `tisserand.py`.
* MGA combinatorics driver (Tisserand-pruned body-sequence enumerator at
  fixed t0): ~200 LOC (most of this is enumeration plumbing).
* Dedupe vs catalogue (existing) and against mission archives (need a
  small mission-database harvest): ~150 LOC.
* Closure-test extension (per-return closure rather than single-period
  closure): ~50 LOC.
* Scorer adaptation: minimal — the existing reachable-set / Tisserand
  scorers translate directly; epoch becomes a free variable in the same
  V∞-graph framework: ~100 LOC.
* Tests + goldens against Galileo VEEGA (1989 launch), Cassini VVEJGA
  (1997 launch): ~200 LOC (these are sourced trajectories — JPL ephemeris
  files are open).

**Total: ~1200 LOC, 2-3 weeks.** Substantially smaller than BCR4BP because
no new dynamical model is needed; smaller than QP-tori because no new
solver is needed; smaller than ER3BP because the propagator + corrector are
already ephemeris-aware (in the V3 validation lane).

### (e) First deliverable IC
**Galileo VEEGA reproduction (sourced):**
```
Launch epoch: 1989-10-18 16:53 UTC
Spacecraft state at launch (post-IUS injection, from JPL Horizons):
  HELIOCENTRIC ecliptic J2000 — query horizons; not memorising here.
Sequence: Earth (launch) → Venus (1990-02-10) → Earth (1990-12-08) →
          Earth (1992-12-08) → Jupiter (1995-12-07)
Closure: Jupiter arrival within km of historical arrival.
```
This is THE textbook epoch-locked MGA. JPL Horizons gives the sourced state
vectors at every node.

**First-novel deliverable IC** (a quasi-cycler of opportunity): the
Hughes-Edelman-Longuski 2014 (corpus) Venus-GA Mars free-return chain
re-instantiated in a 2030-2045 alignment window. Their paper publishes the
ToF and `V∞` budget at synodic alignment; the search dimension is the
launch epoch `t0` and the number of returns the trajectory closes for
before drifting. **The chain is in the existing closure driver — only
the t0-sweep is new.**

### (f) Risk
* **Existence priors are dense.** GTOC, mission archives, and Petropoulos-
  Longuski cover the obvious epoch-locked MGAs. Novelty risk is real;
  re-discovering Galileo VEEGA is a validation, not a discovery.
* **Closure-window relaxation introduces a subjective tolerance.** A
  "3-return quasi-cycler" depends on the tolerance — sensitive to
  parameter choices. Risk of unfalsifiable claims.
* **Scope question is real.** Per `project_catalogue_scope_cyclers_only`,
  the project deliberately limits itself to repeating cyclers. Quasi-cyclers
  + precursor tours are arguably out of scope; pulling them in expands the
  deliverable surface area substantially.
* **V0-V5 gauntlet validation chain may not apply cleanly.** The gauntlet
  is built for static cyclers; epoch-locked candidates have an additional
  "real-ephemeris-at-epoch-X-closes-but-at-epoch-Y-does-not" property the
  gauntlet does not test for. Validation methodology needs a clean spec.

### (g) What it unlocks
* The entire deployed-trajectory class — every published MGA, every quasi-
  cycler in 1980-2020 mission design, becomes search-space-comparable.
* Bridge from the static catalogue (Aldrin S1L1, Russell-Ocampo census) to
  the actual mission-design question: "given a 2035-2045 launch window,
  what's the best Earth-Mars cycler-precursor + insertion?"
* Identifies novel quasi-cyclers in 2025-2050 windows the literature
  hasn't searched. The 2030s Mars synod has not been searched for novel
  precursor tours.
* Closes the gap between V0-V5 validation (which already needs real
  ephemerides) and discovery (currently still mass-parameter-only).

---

## Ranked recommendation

Ranking by **(existence-prior strength × payoff) / (build cost × risk)**.

| Rank | Axis | Existence prior | Build (LOC) | Build (days) | Risk | Payoff | Notes |
|------|------|-----------------|------------:|-------------:|------|--------|-------|
| **1** | **Axis 3 — 3D / broken-plane** | MED-HIGH | ~850 | 7-10 | LOW-MED | Halo/Lyapunov + 3D Aldrin candidates at every body system already in catalogue | Propagator + symmetric corrector ALREADY 3D-capable; smallest gap to a real new discovery surface |
| **2** | **Axis 5 — Epoch-locked / quasi-cycler / precursor** | STRONG | ~1200 | 10-15 | MED (novelty + scope) | Deployed-trajectory class; quasi-cyclers in 2030s Mars synods | All components exist; epoch promotion + closure-window relaxation; scope decision needed |
| 3 | Axis 1 — ER3BP | MED-HIGH (refinement); LOW (cyclers) | ~900 | 10-15 | MED-HIGH | Refined family curves at real `e`; ER3BP-only cyclers (speculative) | Strictly inferior to V3 ephemeris; build is bookkeeping not novelty |
| 4 | Axis 2 — BCR4BP / QBCP | HIGH (POL1/POL2); LOW (cyclers) | ~800 | 10-15 | HIGH (cycler-prior absence) | Sun-synodic-resonant cislunar cyclers IF they exist | Already digested (Andreu note); narrow watch-item, not a primary axis |
| 5 | Axis 4 — QP invariant 2-tori | MED (libration); LOW (cyclers) | ~1500-2000 | 25-40 | HIGH (build complexity + scope) | QP-cyclers (new class); robustness analysis | GMOS algorithm is non-trivial; scope question matches Axis 5's |

### Top-2 next steps

**Rank 1 — Axis 3 (3D / broken-plane).** Next step is a 1-day **scoping
spike**: a single test script that calls
`nrho_continuation.correct_symmetric_nrho` with `z0_guess != 0` at one of the
Ross-RT Earth-Moon (1,1) cycler members (already reproduced in
`cr3bp_general_periodic` tests), then continues in `z0`. If that closes, the
"3D Aldrin-class cyclers" hypothesis becomes immediately testable without
the full Axis-3 build. The spike outcome should drive whether the full
~850-LOC asymmetric-3D-corrector build is funded.

**Rank 2 — Axis 5 (epoch-locked).** Next step is a 1-day **scope-decision
note** (NOT code): write down (i) the exact criterion that distinguishes
"quasi-cycler of opportunity" from "MGA tour" — must be falsifiable; (ii) the
catalogue-scope question (do these go in `catalogue.yaml` alongside repeating
cyclers, or do they become a separate `tour_catalogue.yaml`); (iii) the V0-V5
gauntlet extension a quasi-cycler needs. Until those three are written down,
starting the ~1200-LOC build risks delivering candidates that don't fit
anywhere.

---

## Scope decision section — Axis 5 specifically

The current `project_catalogue_scope_cyclers_only` memory limits the
deliverable to *repeating* cyclers. Quasi-cyclers (3-5 returns in a window,
then drift) and precursor MGA tours (one-shot insertion trajectories) are
both genuinely outside that scope.

Tradeoffs (presented for decision, NOT decided here):

* **Expand scope to include quasi-cyclers + precursors.** Pro: the
  catalogue becomes mission-actionable rather than abstract; the
  user's explicit question ("orbits we can only find at one point in time")
  is answered as a first-class deliverable. Con: surface area expands;
  V0-V5 gauntlet needs an epoch-aware extension; the literature is large
  and the discovery-vs-rediscovery question gets harder; the
  `feedback_published_rounded_values_are_display` and "novel until proven
  via `search/literature_check.py`" rules need extension to mission archives
  (not just paper catalogues).
* **Keep scope tight; carry quasi-cyclers + precursors as a separate
  deliverable class.** Pro: the existing 277-row repeating-cycler
  catalogue stays clean and defensible; quasi-cyclers + precursors get
  their own `tour_catalogue.yaml` or equivalent, with their own validation
  spec and provenance discipline. Con: the existing repeating-cycler
  catalogue is at its validation ceiling (per `project_validation_ceiling`);
  declining to expand into the deployed-trajectory space caps the
  project's discovery payoff at the planar-CR3BP frontier (which #283-#287
  are already exhausting).
* **Hybrid — add epoch-locked novelty only when it converges INTO an
  existing repeating cycler.** Pro: every Axis-5 hit is a first-class
  catalogue row WITH a precursor + insertion trajectory; this is the
  "Aldrin S1L1 + precursor" case. The mission-design question gets
  answered without scope-bloat. Con: misses the
  quasi-cyclers-that-aren't-asymptotic-to-a-repeating-cycler — but per
  `project_s1l1_realeph_closure_blocker`, those may be the bulk of the
  Earth-Mars space.

User decision required before starting Axis 5; this scoping doc
deliberately does not pick.

---

## Summary table

| Axis | LOC | Days | Existence prior strength | Rank |
|------|----:|-----:|--------------------------|-----:|
| 1 — ER3BP | 900 | 10-15 | MED-HIGH (refine) / LOW (cyclers) | 3 |
| 2 — BCR4BP / QBCP | 800 | 10-15 | HIGH (POL1/POL2) / LOW (cyclers) | 4 |
| 3 — 3D / broken-plane | 850 | 7-10 | MED-HIGH | **1** |
| 4 — QP invariant 2-tori | 1500-2000 | 25-40 | MED (libration) / LOW (cyclers) | 5 |
| 5 — Epoch-locked / quasi-cycler / precursor MGA | 1200 | 10-15 | STRONG | **2** |

---

## Surprises and findings worth flagging

1. **Axis 3 is much cheaper than it looks.** The CR3BP propagator and the
   NRHO corrector are *already* 3D. The "planar" constraint in our pipeline
   is a *corrector convention* on the
   `cr3bp_general_periodic.correct_general_periodic` path, not a
   propagator limitation. The 3D-Aldrin-class hypothesis can be tested in a
   1-day spike, not a 2-week build.

2. **The Andreu QBCP digest is already complete.** The 2026-06-14 note
   (`docs/notes/2026-06-14-andreu-quasi-bicircular-digest.md`) digested the
   Hamiltonian, parameter table, and POL1/POL2 dynamical-substitute goldens
   from the open companions (Gimeno-Jorba 2018; Rosales-Jorba 2023) without
   ever opening the actual thesis PostScript. The path to a BCR4BP build
   is already lit — but the prior digest's own verdict (BACKGROUND with one
   narrow watch-item) survives this scoping. No new evidence here
   changes that.

3. **Axis 5 (epoch-locked) is the highest-payoff axis blocked NOT on
   capability but on scope.** All the components — JPL DE440, Lambert,
   FBS, flyby, free_return_chain, multiarc_closure — exist. The build is
   "promote `t0` to a search variable," which is plumbing. The blocker is
   the catalogue scope question, which is a user decision, not an
   engineering question.

4. **Axis 1 (ER3BP) is the weakest of the five.** Strict-CR3BP refinement
   into `e>0` is well-documented but yields refined family curves, not new
   species. V3 (DE440 ephemeris) already dominates fidelity on any
   real-system cislunar question. There is no published ER3BP cycler
   family that doesn't continue back to CR3BP. The build cost is real and
   the payoff is small.

5. **Axis 4 (QP-tori) carries the same scope question as Axis 5** — a
   QP-cycler is not strictly periodic. If the catalogue stays
   `cyclers-only`, QP-cyclers are out of scope and the build pays off only
   for stability-analysis use cases, not discovery. The GMOS solver is the
   most ambitious of the five builds.

6. **All five axes underline the same standing rule** (from
   `project_catalogue_scope_cyclers_only` and `feedback_literature_novelty_check_baseline`):
   the project's deliverable is *novel cyclers*, and "novelty" is gated
   against the published record. Axes 3 and 5 sit closest to that
   deliverable; Axes 1, 2, 4 either refine (1, 2) or expand into adjacent
   classes (4) that need a scope decision first.
