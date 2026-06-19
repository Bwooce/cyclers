# Szebehely 1967 — Theory of Orbits: The Restricted Problem of Three Bodies

**Deep-read verdict note, 2026-06-19 AET.** The long-overdue foundational
CR3BP textbook digest. Textbook-summary scope discipline: this is a
342-PDF-page (661-book-page) image-only scan (NO text layer), so the
Read tool's page-image OCR path was used selectively. Deep-read of the
foundational chapters the project genome implicitly rests on — Ch1
(EOM + Jacobi integral), Ch4 (libration-point location + zero-velocity
curves / Hill regions), Ch5 (libration-point stability), Ch8
(periodic-orbit taxonomy), Ch9 (Copenhagen family classification +
Jacobi-constant convention zoo). Section-header sampling of the rest.
PDF↔book page map: book_page ≈ 2×PDF_page − 16 (the Read tool renders
two facing book pages per PDF "page").

## Header

- **Title:** *Theory of Orbits: The Restricted Problem of Three Bodies*
- **Author:** Victor G. Szebehely (Yale University Observatory at time
  of writing)
- **Publisher:** Academic Press, New York and London
- **Year:** 1967 (1st edition)
- **Pages:** 661 book pages (xvi prefatory + ~647 body + author/subject
  indices); 10 chapters + per-chapter appendices.
- **Format:** The canonical graduate-level / research monograph on the
  circular restricted three-body problem. Comprehensive: history,
  formulation, regularization, zero-velocity theory, libration points,
  stability, Hamiltonian/canonical methods, periodic orbits, numerical
  family catalogues (the Copenhagen school), and modifications
  (3D / elliptic / Hill's problem).

This is THE foundational reference for the entire project genome
(`core/cr3bp.py`, `search/cr3bp_periodic.py`, the L-point solver, and
the Lyapunov / halo / NRHO / tulip families). It has been in the corpus
since an early batch but was never digested — closed here.

## Full table of contents (chapter level, book pages)

| Ch | Title | book p | Project relevance |
|----|-------|--------|-------------------|
| — | Introduction (history: Euler→Lagrange→Jacobi→Hill→Poincaré→Birkhoff→Strömgren) | 1 | context |
| 1 | Description of the Restricted Problem | 7 | **CORE — EOM + Jacobi integral** |
| 2 | (streamline analogies / Hamiltonian flow) | ~45 | low (hydrodynamic analogy) |
| 3 | Regularization | 70 | moderate (Levi-Civita, Birkhoff, Thiele-Burrau, Lemaître) |
| 4 | Totality of Solutions | 126 | **CORE — L-point location + zero-velocity/Hill regions** |
| 5 | Motion near the Equilibrium Points | 231 | **CORE — linear + nonlinear stability** |
| 6 | Hamiltonian Dynamics in the Extended Phase Space | 319 | moderate |
| 7 | Canonical Transformations of the Restricted Problem | 343 | moderate (Delaunay, regularization with canonical vars) |
| 8 | Periodic Orbits | 381 | **CORE — periodic-orbit taxonomy** |
| 9 | Numerical Explorations | 443 | **CORE — Copenhagen family catalogue** |
| 10 | Modifications of the Restricted Problem | 556 | high (3D / elliptic / Hill's problem — frontier) |

Each of Ch1-Ch9 has its own Notes + References section; Ch4 has four
appendices (collinear-point characteristics + Jacobian values), Ch5 has
three (characteristic-equation roots at collinear + triangular points).

---

## Chapter 1 — Description of the Restricted Problem (book p7-41)
### THE canonical equations-of-motion + Jacobi-integral derivation

The project's `core/cr3bp.py` rotating-frame equations and Jacobi
constant are exactly this chapter's results.

- **§1.2 (p8) — sidereal formulation:** two primaries m1 > m2 in
  circular orbit about their barycenter; the third (massless) body's
  inertial EOM d²X/dt*² = ∂F/∂X with force function
  F = k²(m1/R1 + m2/R2) (eq 5). Kepler's third law k²(m1+m2) = n²l²
  (eq 2). The time-dependence of the primaries' positions makes the
  sidereal force function explicitly time-dependent.
- **§1.3 (p10) — energy is NOT conserved in the restricted problem**
  (Part C, eqs 19-23): the third body's energy h3 is not constant
  because the model neglects m3's back-reaction on the primaries; the
  "violation of energy conservation" is the formal signature of the
  restriction (m3 = 0). Important conceptual anchor.
- **§1.4 (p12) — synodic (rotating) formulation + the JACOBI
  INTEGRAL.** The rotation transformation (eq 24-25) to barycentric
  rotating coordinates; the complex-variable derivation (eqs 26-34)
  yielding the rotating EOM
  ẍ − 2nẏ = ∂F*/∂x, ÿ + 2nẋ = ∂F*/∂y (eqs 32-34).
  In the rotating frame the force function F* has NO explicit time
  dependence — this is what admits the Jacobi integral, "the only
  known integral of the restricted problem."
- **§1.8-1.9 (p28-31) — derivation as the degenerate n-body case** and
  the **classification of modifications** (elliptic / 3D / mass-ratio
  μ / time-varying-mass / non-central-force variants). Defines μ = m2/M
  and notes μ = 1/2 is the **Copenhagen problem**; small-μ is
  Poincaré's perturbation regime.
- **§1.10 (p31) — order-of-magnitude justification** that the
  restricted-problem model is valid to ~1% for the Sun-Earth-Moon
  system and to O(10⁻¹⁶) for a lunar space probe (the probe's
  back-reaction is negligible). This is the rigorous warrant for using
  CR3BP at all.
- **§1.10C / §1.11 (p36) — Tisserand's criterion** (eqs 83-86): the
  sidereal Jacobi integral 1/a + 2[a(1−e²)]^(1/2) = C̄ for comet
  identification, and the explicit caution (referencing Ch9) that
  **"the Jacobian integral, even in its exact form, is NOT a good
  measure of orbit accuracy"** — its engineering use for swing-by
  trajectory evaluation is noted but bounded.

### Historical-errata note (project-relevant: respectful-errata
discipline)

Szebehely's own §1.11 Notes (p37) catch **published sign/typesetting
errors in the founding papers**: Poincaré's *Méthodes* Vol 1 p22 prints
r2² where r1² is meant; G.D. Birkhoff's 1915 Palermo paper has the text
mass of S = μ contradicting Fig 1 (where S is clearly the larger mass)
and the r1²/r2² definitions in text vs figure are inconsistent. This is
a 1967 precedent for the project's "typesetting slips happen to
everyone; evidence-first respectful errata" discipline — even Poincaré
and Birkhoff had transcription errors in their CR3BP papers.

---

## Chapter 4 — Totality of Solutions (book p126-211)
### Libration-point location + zero-velocity curves / Hill regions

This chapter is the direct theoretical basis for the project's L-point
solver and zero-velocity / Hill-region tooling.

### §4.4 (p134) — Computation of the collinear libration points

The exact recipe the project's collinear-point solver implements:

- **Bounding regions** (eqs 16a-16c): L1 ∈ [μ−2, μ−1], L2 ∈ [μ−1, μ],
  L3 ∈ [μ, μ+1]; μ-independent bounds −1.3 ≤ x1 ≤ −1, −1 ≤ x2 ≤ 0,
  1 ≤ x3 ≤ 1.2.
- **The three collinear QUINTICS** (one per region), e.g. for L1
  (eq 19): ξ⁵ + (3−μ)ξ⁴ + (3−2μ)ξ³ − μξ² − 2μξ − μ = 0; for L2 and L3
  similarly (eqs 26, 32). Descartes' sign rule gives exactly one
  positive root for 0 < μ ≤ 1/2.
- **Iteration recipe** (eq 20-21): ξ³ = μ(1+ξ)²/[3−2μ+ξ(3−μ+ξ)],
  with the recommended starting value **ξ = [μ/3(1−μ)]^(1/3)** —
  precisely the standard cubic-root seed for L1/L2 used across the
  field. Series solutions in ν = [μ/3(1−μ)]^(1/3) given to O(ν⁷)
  (eqs 22, 27, 34).
- **§4.5 + Appendices I-III (p138, 214-225)** tabulate L1/L2/L3
  characteristics over μ ranges, with **Appendix ID specifically for
  the Earth-Moon μ ≈ 0.012** (21 equidistant entries 0.011 ≤ μ ≤
  0.013) and Appendix IC for Jupiter μ ≈ 0.00095. The "critical" entry
  marks the Routh stability boundary.

### §4.7 (p159-165) — Regions of motion / zero-velocity curves

The canonical Hill-region treatment the project's zero-velocity tooling
computes:

- **The Ω function** (eq 95):
  Ω = ½[(1−μ)r1² + μr2²] + (1−μ)/r1 + μ/r2,
  with the Jacobi integral v² = 2Ω − C.
- **Zero-velocity curves** defined by 2Ω − C = 0 (the equipotential /
  Hill-region boundaries); motion is possible only where Ω > C/2
  (v² ≥ 0).
- **The 3D Ω(x,y) surface picture** (p165): two infinite peaks at m1
  and m2, three saddle "passage ways" at L1/L2/L3, two minima at L4/L5
  with Ω = 3/2 (C = 3). The C-ordered topology of the Hill regions
  (outer/inner ovals for large C; necks opening at L1/L2 as C
  decreases) — this is the Belbruno-2004-Ch3 Hill-region figure in its
  original textbook form. The §4.7.1 two-body warm-up (eqs 74-93,
  circles of zero velocity) builds the intuition.

---

## Chapter 5 — Motion near the Equilibrium Points (book p231-308)
### Linear + nonlinear stability of the libration points

The theoretical basis for the project's Lyapunov / halo family
generation (linearized periodic orbits → analytic continuation) and for
the L4/L5 stability boundary.

### §5.3-5.4 — collinear vs triangular linear stability

- **Collinear points are ALWAYS UNSTABLE** (§5.5.2, p266): the
  characteristic equation has a positive real root for ANY μ, so the
  general linearized solution is unbounded. BUT particular solutions
  exist — **infinitesimal (linearized) periodic orbits** corresponding
  to the imaginary-root pair. Their period **depends only on μ and is
  independent of the (infinitesimal) orbit size** — "a characteristically
  linear phenomenon." When higher-order terms are included, the mean
  motion changes with orbit size (eq 82a, via Horn's theorem analytic
  continuation). **This is the foundational existence statement for the
  project's Lyapunov family** (the collinear-point linear periodic orbit
  continued to finite amplitude).
- **Triangular points L4/L5** (§5.4): linearly stable iff μ < the Routh
  critical value **μ0 = 0.0385208965...** (the root of 27μ(1−μ) = 1).
  At μ = μ0 the characteristic roots become the double pair
  λ = ±i/2^(1/2), giving SECULAR terms (eq 74) → unstable. For μ < μ0
  the motion is bounded (two real frequencies → the long/short-period
  librations). Earth-Moon μ ≈ 0.0123 < μ0, so L4/L5 are stable — the
  basis for the project's Trojan / tadpole / horseshoe-adjacent work.
  Ω*_xx(L4) = 3/4, Ω*_yy(L4) = 9/4, Ω*_xy(L4) = (3·3^(1/2)/2)(μ−1/2)
  (the linearization coefficients, p264).

### §5.5 — Nonlinear phenomena (Horn's theorem continuation)

The Taylor expansion of Ω about a libration point gives
ξ̈ − 2η̇ − Ω_xx ξ − Ω_xy η = X(ξ,η;μ), etc. (p266), and the existence
of finite-amplitude periodic orbits as **analytic continuation** of the
infinitesimal linear orbits follows from **Horn's theorem**. This is
the textbook justification for the project's
"linearize-then-differential-correct-then-continue" Lyapunov/halo
generation pipeline.

---

## Chapter 8 — Periodic Orbits (book p381-440)
### The periodic-orbit taxonomy the genome's nomenclature extends

§8.2 (p384) is the foundational definitions chapter for the project's
resonant-cycler / quasi-periodic vocabulary:

- **Periodic in sidereal vs synodic systems** (p384): T_sid = 2π/n;
  the synodic path closes after |p| sidereal periods when n = p/q
  rational, giving **T_syn = 2πq**. Circular orbits are periodic in
  both frames; elliptic orbits need not be — this distinction is the
  origin of **Poincaré's periodic orbits of the FIRST and SECOND
  kind**. (Directly relevant to the project's S/L resonant-interval
  nomenclature; cf. memory: S1L1 nomenclature.)
- **Almost-periodic (Bohr) and quasi-periodic (Bohl)** functions
  (p385): quasi-periodic = finite set of basic incommensurable
  frequencies ωᵢ, x = X(ω1 t, ..., ωn t). **This is the formal
  definition of the QP-tori the project's frontier work targets** (the
  ER3BP/BCR4BP/QP frontier per memory: speculative-high-effort).
- **Small divisors / commensurability / resonance** (p386, item D):
  Σωᵢkᵢ = 0; the "great inequality" 5T_Jupiter − 2T_Saturn example;
  resonance = the small-divisor phenomenon. The substrate for resonant
  cyclers.
- **The torus concept** (p387, items E-F): periodic motion on a torus
  when frequency ratios are rational; for irrational ratios the orbit
  is dense and uniformly distributed (ergodic on the torus); Poincaré
  recurrence with recurrence time T* = nᵢ·2π/ωᵢ. This is the
  dynamical-systems backbone for invariant-tori cycler families.

§8.3-8.8 (skimmed): surface-of-section / torus representation, analytic
continuation, the small-mass-parameter expansion (Poincaré's first/second
kind), **Whittaker's criterion for periodic-orbit existence** (p428),
characteristic exponents (p430 — the project's Floquet-multiplier
substrate). Szebehely explicitly states (p384) that a full periodic-orbit
treatment is "far beyond the scope" and that Ch8 is "an outline of
another volume" — the analytic periodic-orbit theory is deliberately
compressed here, with the NUMERICAL families deferred to Ch9.

---

## Chapter 9 — Numerical Explorations (book p443-550)
### The Copenhagen family catalogue + the Jacobi-constant convention zoo

This is the historical origin of periodic-orbit FAMILY classification —
the lineage the project's catalogue family-structure descends from.

- **§9.2 / Table II (p454) — the Jacobi-constant convention ZOO.** Nine
  historical conventions (Standard, Birkhoff, Wintner, Darwin,
  Strömgren, Moulton, Broucke, Charlier, Rabe) each with a different
  mass-parameter normalization and a different Jacobi-constant
  definition, with explicit conversion relations. **The
  Copenhagen/Strömgren constant C_K = 4C̄** (eq 11), Darwin's
  C_D = 11C̄ at μ = 10/11, etc. **This is directly load-bearing for the
  project's published-value cross-referencing** (memory: published
  rounded values are display; the C21 family-extent saga): a published
  Jacobi constant is meaningless without knowing WHICH convention and
  WHICH μ-normalization the author used. Any literature_check / mining
  pass that compares a published C against a project-computed C MUST
  first identify the convention (Table II is the Rosetta stone).
- **§9.4 (p455) — the Copenhagen category.** Elis Strömgren's
  classification of periodic-orbit families, based on the **seven
  special points** (5 equilibria L1-L5 + the 2 primary positions),
  computed at μ = 1/2 (unit mass ratio, y-axis symmetric). Families
  are named by the special point they encircle and their direct/
  retrograde sense, e.g. "retrograde periodic orbits around L3" (direct
  don't exist), the families around the primaries, asymptotic orbits to
  L4/L5. Table III (p457) tabulates the seven points' coordinates +
  three Jacobi conventions for μ = 1/2. **This seven-point,
  named-family taxonomy is the ancestor of the project's family
  classification.**
- §9.5-9.7 (skimmed): periodic lunar orbits, motion around the
  triangular points, lunar trajectories. §9.10 stability of the
  numerical families.

---

## Chapters skimmed / sampled (honest relevance assessment)

- **Ch2 (streamline analogies, ~p45):** Hamiltonian-flow /
  hydrodynamic-streamline analogy for the restricted problem.
  Mathematically elegant, low direct project relevance.
- **Ch3 (Regularization, p70):** Levi-Civita local regularization,
  Birkhoff global regularization, Thiele-Burrau and Lemaître
  transformations for collision orbits. **Moderate relevance** if the
  project ever needs close-approach / collision-orbit regularization
  (e.g., very-low-periapsis flyby legs); the project's current
  integrators do not regularize. Flag as a reference if numerical
  trouble near primaries ever arises.
- **Ch6-7 (Hamiltonian / canonical transformations, p319-379):**
  extended phase space, canonical transformation sidereal↔synodic,
  Delaunay elements, **regularization with canonical variables**
  (§7.8). Moderate relevance for any future
  Hamiltonian/symplectic-integrator or action-angle work; pairs with
  the Gurfil-Ch3 generating-function machinery (#380).
- **Ch10 (Modifications, p556):** the **three-dimensional restricted
  problem** (§10.2, p557 — the basis for halo / NRHO / 3D families),
  the **elliptic restricted problem** (§10.3, p587 — the ER3BP
  frontier), **Hill's problem** (§10.4, p602). **HIGH relevance to the
  project's frontier work** (memory: speculative-high-effort —
  ER3BP/3D are the multi-week frontier). Not deep-read here (out of the
  "foundational substrate" scope of this digest), but flagged as the
  next Szebehely chapter to deep-read if/when the ER3BP or 3D frontier
  fires.

---

## What the genome SHOULD cite but (likely) doesn't

The project genome implements all of the following Szebehely-1967
results but, being foundational "everyone-knows-it" CR3BP theory, almost
certainly cites none of them. Recommend a single methodology-bibliography
anchor (NOT a catalogue/KNOWN_CORPUS family entry — this is substrate
theory):

1. **Rotating-frame EOM + Jacobi integral** (§1.4, eqs 32-34) —
   `core/cr3bp.py`.
2. **Collinear-point quintics + the [μ/3(1−μ)]^(1/3) iteration seed**
   (§4.4, eqs 19-22) — the L-point solver.
3. **Routh critical mass μ0 = 0.0385208965** for L4/L5 linear stability
   (§5.4) — any Trojan / triangular-point work.
4. **Ω(x,y) + zero-velocity-curve / Hill-region topology** (§4.7,
   eq 95) — the zero-velocity tooling.
5. **Collinear-point linear periodic orbit (period independent of
   amplitude) + Horn's-theorem continuation** (§5.5.2) — the Lyapunov /
   halo generation pipeline.
6. **Periodic-orbit-of-first/second-kind + quasi-periodic-torus
   definitions** (§8.2) — the resonant-cycler / QP-frontier vocabulary.
7. **The Jacobi-constant convention table** (§9.2, Table II) — MUST be
   consulted by any literature_check / mining pass that cross-references
   a published Jacobi constant (different schools use C, C̄, C_K, C_D…).

## KNOWN_CORPUS impact

**No catalogue writeback. No cycler-family anchor.** Szebehely 1967 is
substrate CR3BP theory, not a cycler-family source — there is no named
cycler member in the book (cyclers as a mission concept postdate it).

**Methodology-bibliography entry RECOMMENDED.** Szebehely 1967 should be
the canonical foundational citation for the project's CR3BP genome.
Citation:

  Szebehely, V.G. (1967). *Theory of Orbits: The Restricted Problem of
  Three Bodies*. Academic Press, New York and London.

Specific sections to cite are listed in "What the genome SHOULD cite"
above.

## Action items for the parent

1. **Lodge the foundational methodology anchor.** Szebehely 1967 is now
   the canonical CR3BP-genome reference. Recommend adding it to the
   methodology bibliography (alongside Belbruno 2004, Hintz 2023) with
   the section-specific citations above. No catalogue edit.
2. **Jacobi-constant-convention hazard flagged.** §9.2 Table II is the
   convention Rosetta stone. Any mining / literature_check that compares
   a published Jacobi C against a project-computed C should first
   identify the author's convention + μ-normalization. This is a
   concrete reinforcement of the "published rounded values are display"
   memory — the convention difference (C vs C_K = 4C̄) can be an
   order-1 discrepancy, not a rounding one.
3. **Historical-errata precedent recorded.** Szebehely catches
   Poincaré and Birkhoff sign/typeset errors in their founding CR3BP
   papers (§1.11). Reinforces the project's respectful-errata
   discipline: even the founders had transcription slips.
4. **Ch10 is the next deep-read when the frontier fires.** §10.2 (3D
   restricted), §10.3 (elliptic restricted / ER3BP), §10.4 (Hill's
   problem) are the foundational theory for the project's frontier
   (3D halo/NRHO families + ER3BP/BCR4BP). Deliberately NOT deep-read
   here (this digest is scoped to the planar-CR3BP substrate); flag for
   a follow-up digest if/when #347-style frontier work graduates.
5. **Honest negative on cycler corpus.** Zero cycler families. This is
   the bedrock theory the whole genome stands on, not a corpus
   contribution. Its value is (a) a canonical foundational citation and
   (b) the Jacobi-convention table for cross-referencing.
