# #733 — Olikara 2016 PhD thesis + Haro et al. 2016 parameterization-method
# book: two more user-supplied top-5 items from the `#730` backlog

**Task:** `#733`, items #2 and #3 of `#730`'s own §2 top-5 cluster, PDFs
supplied directly by the user (page-1-verified by the coordinating session
before dispatch), found in `~/Downloads`. Sibling task `#732` concurrently
processed the other three top-5 items (Blazevski & Ocampo 2012, Negri & Prado
2020, Baresi/Olikara/Scheeres 2018) — its own digest
(`2026-07-27-732-blazevski-negri-baresi-foundational-papers-digest.md`) is
read and cross-referenced below where relevant (the Baresi/Olikara/Scheeres
2018 paper is itself a head-to-head GMOS-vs-PDE comparison directly adjacent
to this task's own Olikara-thesis cross-check).

## 0. Acquisition, filing, text-layer check

Both filed in the private `cyclers_pdf` corpus
(`/Users/bruce/dev/cyclers_pdf/papers/`, separate repo, committed there as
`1d796f3`, never committed to the public `cyclers` repo):

| Paper | Filename | Pages | Text layer |
|---|---|---|---|
| Olikara 2016 PhD thesis | `olikara-2016-computation-quasi-periodic-tori-heteroclinic-connections-collocation-phd-thesis-colorado.pdf` | 113 | Native, clean (`pdftotext -layout` extraction verified readable including all equations; `pdfinfo` shows a digitally-produced PDF 1.3, PyPDF2-processed) |
| Haro, Canadell, Figueras, Luque & Mondelo 2016 book | `haro-canadell-figueras-luque-mondelo-2016-parameterization-method-invariant-manifolds-book-doi-10.1007-978-3-319-29662-3.pdf` | 280 | Native, clean (Adobe InDesign-produced PDF 1.4, embedded fonts; `pdftotext -layout` extraction verified) |

No OCR needed for either — both are digitally-typeset originals with
embedded, extractable text layers. Read via the extracted text, applying the
"chapter-summary scope" convention: full TOC read for both; for the thesis,
Chapters 1 (Introduction) and 2 (Computing Tori in Hamiltonian Systems, the
core algorithm chapter) were deep-read in full, Chapter 3's quasi-periodic
example (§3.1.3) was deep-read for its numeric positive control, Chapters
4-6 (heteroclinic connections, ephemeris transitioning, conclusion) were
sampled via the author's own Chapter 1 §1.3 "Contributions" synopsis of each;
for the book, Chapter 1 (historical overview, deep-read in part) and the
opening of Chapter 3 ("skew-product systems / fiberwise hyperbolic invariant
graphs," the whiskered-quasi-periodic-tori chapter) and Chapter 5 ("A
Newton-like Method for Computing Normally Hyperbolic Invariant Tori," deep-read
in substantial part) were the load-bearing chapters for this project's own
whiskered/hyperbolic-torus work; Chapters 2 and 4 (fixed-point manifolds,
KAM theory) were sampled via TOC + targeted grep, not deep-read line-by-line.

---

## 1. Olikara 2016 — "Computation of Quasi-Periodic Tori and Heteroclinic
Connections in Astrodynamics Using Collocation Techniques," PhD thesis,
University of Colorado Boulder (advisor D. J. Scheeres)

### 1.1 What the thesis actually does

Six chapters: **Ch. 1** introduces quasi-periodic tori (invariance equation,
reducibility, Hamiltonian setting) and surveys prior computational approaches
(§1.2, see below); **Ch. 2** is the thesis's own core algorithmic
contribution — a general-purpose **Gauss–Legendre collocation** scheme for
computing 2D quasi-periodic tori (autonomous, periodically-forced, and
quasi-periodically-forced Hamiltonian systems) plus their linear stability;
**Ch. 3** demonstrates the method on CR3BP/ER3BP/HR4BP examples; **Ch. 4**
computes heteroclinic connections between tori in the spatial CR3BP (indirect,
extending Calleja et al. 2012's periodic-orbit approach, and direct, a BVP for
the connection family); **Ch. 5** transitions CR3BP quasi-periodic
orbits/connections to an ephemeris model via a closest-trajectory
optimization; **Ch. 6** concludes.

### 1.2 The central, checkable fact: the thesis's own explicit shooting-vs-
collocation framing, and where `genome/qp_tori.py` actually sits

The thesis's own §1.2.1 (p. 12) states the shooting/collocation distinction
in exactly the terms this project's own code needs to be checked against:

> "With a shooting method, we essentially follow the flow of the vector field
> and then iterate to satisfy continuity and boundary constraints. On the
> other hand, with a collocation method we start with a continuous
> trajectory, and then iterate such that it matches the vector field and
> boundary constraints... Before convergence, a shooting method exactly
> follows the vector field (up to the order of the numerical integrator) but
> is not periodic. A collocation trajectory, on the other hand, is exactly
> periodic but does not match the vector field at the collocation points
> until convergence." (Fig. 1.5 illustrates both, side by side.)

And in §1.2.2 (p. 13), the thesis's own literature survey places its direct
predecessor explicitly:

> "A general purpose shooting-based scheme for computing tori parameterized
> by states is presented in **Olikara & Scheeres [2012]**."

And in §1.3 "Contributions" (p. 15), the thesis states its own relationship
to that predecessor in one sentence:

> "Many concepts follow the approach of **Olikara & Scheeres [2012]** but
> rather than using a shooting method, a **Gauss–Legendre collocation**
> method is used."

**This is the exact fact needed to answer the mandatory cross-check.**
`genome/qp_tori.py`'s own `_gmos_residual` propagates each Fourier-mode sample
point forward by the **entire** stroboscopic time `t_strob` via a single
nonlinear numerical integration (`cr3bp.propagate(system, u_samples[j],
t_strob, with_stm=False)`), then FFT-matches the propagated points against
`c_n * exp(i n rho)` — i.e., each Newton/least-squares iterate **exactly
follows the vector field but is not yet invariant until convergence**. This
is, by the thesis's own definition and its own figure, precisely a
**shooting** method, not a collocation method: there is no piecewise
polynomial mesh over the stroboscopic-time direction, no Gauss–Legendre
collocation points, no sparse block-structured Jacobian (§2.3.2/Fig. 2.1),
and no continuity constraints between mesh intervals (Eq. 2.34). The
module's own residual equation, `_gmos_residual`'s
`FFT[phi_{t_strob}(u_j)]_n - c_n * exp(i n rho) = 0`, is algebraically
identical to the thesis's own mode-space invariance condition
`X(0) = R(-ρ)X(1)` (Eq. 2.29, itself derived from the general invariance
equation 2.1a/2.1b and the earlier flow-based equation 1.4) — but the thesis
uses that *same* invariance equation as the boundary condition for its
collocation scheme, whereas `qp_tori.py` uses it directly as a shooting
residual, exactly the pattern the thesis's own Olikara & Scheeres [2012]
citation describes.

**Direct answer: `qp_tori.py` does NOT implement Olikara 2016's thesis
method.** It implements the *earlier* Olikara & Scheeres 2010/2012 GMOS
**shooting** algorithm — the thesis's own explicitly-cited predecessor, the
one this thesis was written to supersede with collocation. This is not a bug
or an error; it is a genuine, sourced, previously-undocumented gap between
what this project's module docstring cites (both Olikara & Scheeres 2010
*and* Olikara 2016 together, as if interchangeable lineage) and what is
actually implemented (only the 2010/2012 shooting predecessor — none of the
2016 thesis's own collocation contribution).

Additional specific matches confirming the *shooting*-method identification
(not merely a family resemblance):

- **Seed formula.** Thesis Eq. 2.3 (§2.1.4, p. 23): for a Neimark-Sacker
  eigenvector `y` of the monodromy at eigenvalue `e^{i2πρ}`,
  `û2(0;θ2) = Re[e^{i2πθ2} y] = cos(2πθ2) Re[y] − sin(2πθ2) Im[y]`. This is
  **term-for-term identical** to `qp_tori.py`'s own
  `_seed_invariant_circle`: `u_seed(theta) ~= s_parent + amplitude *
  (Re(v) cos theta − Im(v) sin theta)`.
- **Invariance condition.** Thesis Eq. 1.4 (`w(θ2+ρ) = φ(w(θ2))`, the flow
  form) and Eq. 2.29 (`X(0) = R(−ρ)X(1)`, the mode-space form after Fourier
  discretization) are the *same* equation `qp_tori.py`'s `_gmos_residual`
  solves — both parameterize the invariant circle by truncated Fourier
  modes of the state and match the stroboscopic-map image to a
  rotation-shifted copy of the circle.
- **What is missing.** The thesis's *own* contribution — §2.3.2's
  Gauss-Legendre-collocation piecewise-polynomial representation of each
  torus-surface trajectory (Eq. 2.31-2.34), the associated sparse Jacobian
  (§2.3.3, Fig. 2.1), the error-control/mesh-refinement scheme (§2.3.4,
  Russell & Christiansen 1978), and the free-by-product Floquet stability
  extraction from that same collocation Jacobian (§2.4) — has **no
  counterpart anywhere in `qp_tori.py`**. The module computes no torus
  stability/Floquet information at all, consistent with the parallel gap
  `#732`'s own Baresi/Olikara/Scheeres-2018 digest found in this project's
  *other* torus correctors (`variational_qp_torus.py` /
  `variational_crnbp_torus.py`, a PDE(DFT)-class method with the same
  stability-computation gap, per that digest's §3.2).

**Consistency with `#732`'s own finding, not a contradiction.** `#732`'s
Baresi/Olikara/Scheeres-2018 cross-check independently concluded
`genome/qp_tori.py` is "GMOS-lineage" (i.e. the Olikara & Scheeres
2010/2012 shooting method) and separately flagged the thesis's own
collocation contribution as unimplemented and higher-priority to acquire.
This task's own direct read of the thesis text **confirms and sharpens**
that finding with the thesis's own explicit self-description (quoted above)
rather than merely inferring it from a citation elsewhere.

**Is this a problem?** Not urgently. `#732`'s own analysis (citing this
project's `variational_qp_torus.py` docstring) established that GMOS-class
shooting methods converge cleanly for this project's actual use case
(low-amplitude tori near a mildly-unstable-to-stable parent orbit) and the
Baresi et al. 2018 paper's own head-to-head comparison found GMOS *faster and
more accurate* than any PDE-class alternative for such cases. The gap that
matters is: (1) no Floquet/stability output from `qp_tori.py` (a genuine,
reusable enhancement this thesis's Ch. 2.4 directly specifies how to add,
would require switching to collocation to get "for free" per the thesis's own
claim, or could be added to the existing shooting scheme by finite-differencing
the one-period map directly — a cheaper partial fix); and (2) the
thesis's own collocation scheme's chief practical advantage — much better
conditioning/sparsity when scaling to larger Fourier order `N` or when the
one-period shooting propagation itself becomes numerically fragile (the same
"stiff hyperbolic parent orbit" pathology `#611`/`#612` document for the PDE
correctors) — is not available as a fallback if `qp_tori.py`'s shooting
approach ever hits that wall for its own use cases. Neither is currently
observed to be blocking any live search.

### 1.3 Positive-control numeric example (Ch. 3, §3.1.3)

Earth-Moon L2 quasi-halo torus family emanating from an L2 halo orbit, all
members at fixed Jacobi constant `C = 3.132` (`H = -1.566`, Fig. 3.8). Each
torus computed (unparallelized) in ~2 s using Fourier order `M2 = 25`
(`N2 = 51` states around the circle), `m = 6` Gauss-Legendre collocation
points per mesh interval, global error tolerance `1e-10` along each
trajectory (Jacobi constant preserved to `1e-15` along each mesh trajectory).
A companion L1 planar-Lyapunov-seeded Lissajous-torus family and both
families' `z=0`-plane center-manifold crossings at the same `C = 3.132` are
also shown (Fig. 3.13). This is a clean, fully-sourced numeric example
(energy level, Fourier order, mesh order, tolerances all stated) reusable as
a positive control for any future GMOS-vs-collocation benchmark this project
might run on its own EM L1/L2 quasi-halo tori.

### 1.4 Mandatory citation-mining pass (Ch. 1 §1.2 + bibliography)

Read the full Introduction/Computational-Approaches survey (§1.1-1.2) and
cross-checked all cited names against `CORPUS_INDEX.md`:

**Already in corpus (false gaps, confirmed, not re-acquired):** Andreu 1999
(`andreu-1998-quasi-bicircular-problem-phd-thesis.pdf` — the thesis's own
1999 citation year is off by one from the corpus's 1998, same document);
Koon, Lo, Marsden & Ross 2000, "Heteroclinic connections between periodic
orbits and resonance transitions in celestial mechanics," *Chaos* 10
(`koon-lo-marsden-ross-2000-heteroclinic-connections-resonance-transitions-chaos-10-2.pdf`);
Gómez, Koon, Lo, Marsden, Masdemont & Ross 2004, "Connecting orbits and
invariant manifolds in the spatial restricted three-body problem,"
*Nonlinearity* 17
(`gomez-koon-lo-marsden-masdemont-ross-2004-connecting-orbits-spatial-rtbp-nonlinearity-17.pdf`).

**Genuinely new candidates, not previously flagged by `#730`/`#732`:**

- **Kolemen, E., Kasdin, N. J. & Gurfil, P. [2012], "Multiple Poincaré
  sections method for finding the quasiperiodic orbits of the restricted
  three-body problem," *CMDA* 112(1):47-74.** Already independently flagged
  by `#732`'s own Baresi-Olikara-Scheeres-2018 citation-mining pass (the
  "KKG" method source) — **no new action**, this is a second, independent
  hit reinforcing the existing flag.
- **Schilder, F., Osinga, H. M. & Vogt, W. [2005], "Continuation of
  quasi-periodic invariant tori," *SIAM J. Appl. Dyn. Syst.* 4(3):459-488.**
  Same — already flagged by `#732` (the PDE(CD) method source, also directly
  cited in this project's own `variational_qp_torus.py` docstring). **No new
  action.**
- **Calleja, R. C., Doedel, E. J., Humphries, A. R., Lemus-Rodríguez, A. &
  Oldeman, E. B. [2012], untitled in the thesis's own in-text citation but
  described as computing "a hyperbolic periodic orbit and indirectly [finding]
  a connection to a quasi-periodic orbit" via a collocation-computed unstable
  manifold and continuation on the time-to-section.** **Not in corpus** —
  a genuinely new, distinct paper from the already-flagged Calleja &
  de la Llave 2010 (`#730` §2 item 8) and Calleja/del-Castillo-Negrete/
  Martínez-del-Río/Olvera 2021 (`#730` §2 item 7). Directly on-point for this
  project's own periodic-to-quasi-periodic and quasi-periodic-to-quasi-periodic
  connection work (`#701`-`#708` Umbriel-Titania torus-homoclinic connection).
  **Medium-high priority new gap.**
- **Canalias, E., Delshams, A., Masdemont, J. J. & Roldán, P. [2006], "The
  scattering map in the planar restricted three body problem."** **Not in
  corpus** (distinct from the already-acquired `canalias-2007-thesis`, which
  is Canalias's own PhD thesis, a different, later document). The "scattering
  map" technique is directly on-point for characterizing torus-to-torus
  heteroclinic connection geometry. **Medium priority new gap.**
- **Delshams, A., Masdemont, J. J. & Roldán, P. [2008], "Computing the
  scattering map in the spatial Hill's problem."** Companion spatial
  extension of the item above. **Not in corpus. Medium priority, companion
  gap.**
- **Arona, L. & Masdemont, J. J. [2007], "Computation of heteroclinic orbits
  between normally hyperbolic invariant tori" (spatial Hill R3BP).** **Not
  in corpus.** Directly on-point title for this project's own
  torus-to-torus heteroclinic connection work. **Medium-high priority new
  gap.**
- **Jorba, À. & Olmedo, E. [2009], "On the computation of reducible
  invariant tori in a parallel computer," *SIAM J. Appl. Dyn. Syst.*
  8(4):1382-1404.** **Not in corpus.** A combined torus+stability method
  specific to non-autonomous (known-forcing-frequency) systems — directly
  relevant to this project's own CCR4BP/CRNBP (non-autonomous, time-periodic)
  torus correctors and their currently-missing stability output (§1.2 above).
  **Medium priority new gap.**
- **Baresi, N., Olikara, Z. P. & Scheeres, D. J. [2016], "Survey of numerical
  methods for computing quasi-periodic..." (title truncated in the thesis's
  own in-text citation).** Likely an earlier conference-proceedings version
  of the already-acquired 2018 *JAS* paper (same three authors, overlapping
  scope) — **not independently verified as a distinct document this pass**;
  flagged low priority (probable duplicate/precursor of an already-acquired
  paper, per the same "check before acquiring" caveat `#730`'s own §1 already
  applies to several other probable-duplicate entries).
- **Lower-priority classical CR3BP background, confirmed absent but
  orthogonal to this project's CCR4BP/CRNBP/torus-heteroclinic search
  domain, not pursued further:** Farquhar & Kamel 1973 (Earth-Moon L2
  Poincaré-Lindstedt expansion); Richardson & Cary 1975 (Sun-Earth L1/L2
  multiple-time-scales expansion); Howell & Pernicka 1988 (early
  quasi-periodic shooting corrector); Gómez, Masdemont & Simó 1998 (quasihalo
  Poincaré-Lindstedt expansions); Jorba & Masdemont 1999 (center-manifold
  reduction); Muñoz-Almaraz, Freire, Galán, Doedel & Vanderbauwhede 2003
  (periodic-orbit continuation in Hamiltonian systems); Mohn & Kevorkian 1967
  ("Some limiting cases of the restricted four-body problem" — a historical
  curiosity given this project's own N=4/N=5 model work, but a 1967
  limiting-case analysis, not a numerical method); Castella & Jorba 2000
  (four-body quasi-periodic orbits).

---

## 2. Haro, Canadell, Figueras, Luque & Mondelo 2016 — *The Parameterization
Method for Invariant Manifolds: From Rigorous Results to Effective
Computations*, Applied Mathematical Sciences vol. 195, Springer

### 2.1 What the book actually does

A 5-chapter monograph unifying the parameterization method across four
contexts: **Ch. 1** overview (historical survey + the general invariance-equation
framework for both fixed points and tori); **Ch. 2** invariant manifolds of
fixed points of vector fields (power-series/automatic-differentiation
algorithms, 3 worked examples including the Earth-Moon L1 spatial CR3BP
center manifold); **Ch. 3** "The Parameterization Method for Quasi-Periodic
Systems: From Rigorous Results to Validated Numerics" — fiberwise hyperbolic
invariant graphs (whiskered tori) of **skew-product** (quasi-periodically
forced) systems, directly built on Haro & de la Llave's own 2006ab/2007 papers
(`HdlL06b`, `HdlL06c`, `HdlL07` in the book's own citation keys), including
two full validation theorems and computer-assisted-proof examples; **Ch. 4**
the parameterization method in KAM theory (existence/persistence of
quasi-periodic tori, built on de la Llave, González, Jorba & Villanueva 2005,
`dlLGJV05`); **Ch. 5** "A Newton-like Method for Computing Normally
Hyperbolic Invariant Tori" — the most general-purpose numerical chapter,
specifying one Newton-like step to correct BOTH the torus parameterization
`K` (solving the invariance equation `F∘K = K∘f`, Eq. 5.1/5.5) AND its stable/
unstable normal bundles `N^S, N^U` (solving the bundle invariance equation
5.9), via an adapted frame `P = (L, N)` and block-diagonalized linearized
dynamics `Λ = diag(Λ_L, Λ_S, Λ_U)` (Eq. 5.7-5.11).

### 2.2 Cross-check: does the book's theoretical framework match how the
already-acquired Kumar-lineage papers describe using it?

**Yes — confirmed consistent, no discrepancy found.** The book's own
Chapter 3 (whiskered tori of quasi-periodic maps) is explicitly built on
exactly the three papers (`HdlL06b`/`HdlL06c`/`HdlL07` = Haro & de la Llave
2006/2006/2007) that `#730`'s own §5 cluster (items 31-36) and `#728`'s
whiskered-tori digest already independently identified as the direct
theoretical ancestors of the already-acquired Kumar, Anderson & de la Llave
2021/2022 CMDA/CNSNS papers — this project's citation trail (Kumar papers →
Haro & de la Llave 2006/2007 → this book's own Ch. 3) is now confirmed from
the book's own text, not merely inferred through an intermediate citation.
Likewise, the book's Chapter 4 (KAM theory) is explicitly built on de la
Llave, González, Jorba & Villanueva 2005 (`dlLGJV05`, already flagged as
`#730` §5 item 33, "KAM theory without action-angle variables" — itself
cited by the already-acquired 2022 CMDA whiskered-tori paper for its
"center-bundle/symplectic-conjugate concept"). The book's Chapter 5 Newton-like
method — correcting the torus AND its stable/unstable bundles in a single
generalized scheme via the adapted-frame/block-triangular-linearization
construction (Eq. 5.7-5.11) — is the fully worked-out, textbook-level version
of the same "whiskered tori" object class (a torus plus its 1D or higher
stable/unstable manifolds, corrected together) the Kumar-lineage papers'
own titles describe computing ("Rapid and accurate methods for computing
whiskered tori...", per `#730`'s own §1 false-gap table). No discrepancy
between the book's own theoretical framework and how those already-acquired
papers describe using it was found — this is a genuine, sourced
confirmation, not a contradiction to flag.

One nuance worth recording: the book's own Chapter 5 is scoped to
**diffeomorphisms** (discrete maps), explicitly noting a flow-based version
"can be developed... or one can use an appropriate Poincaré map leading to
the FE [functional equation] approach" (§5.1, p. 190) — i.e. the book itself
treats the continuous-flow case as a reduction to the discrete-map case via a
Poincaré section, not as an independently-worked-out chapter. This matches
this project's own practice (e.g. `genome/qp_tori.py`'s stroboscopic-map
reduction of the CR3BP flow to a 1-period return map) but is worth noting as
a structural simplification the book itself flags, not something unique to
this project's implementation.

### 2.3 Mandatory citation-mining pass (Ch. 1 historical overview +
Ch. 3 introduction + reference list)

The book's own ~700-entry reference list is overwhelmingly pure
dynamical-systems/KAM-theory mathematics (functional analysis, computer
algebra, ergodic theory) not topically overlapping this project's applied
CCR4BP/CRNBP discovery-search domain — most of the specific
astrodynamics-adjacent lineage this book itself draws on (Haro & de la Llave
2006/2006/2007, de la Llave/González/Jorba/Villanueva 2005, Cabré/Fontich/de
la Llave 2003ab/2005, Fontich/de la Llave/Sire 2009, Calleja/de la Llave
2010, Capiński/Gidea/de la Llave 2016, Zhang/de la Llave 2018) is **already**
independently catalogued as a flagged, not-yet-acquired gap by `#730`'s own
§5 cluster (items 31-36) — this pass **confirms and reinforces** those
existing flags (all of the above citation keys appear repeatedly throughout
the book's own Ch. 1/3/4/5, as expected since Haro is himself a co-author of
that lineage) rather than surfacing new ones. Two items sampled from the
book's own reference list that are new, but low priority, background-only,
and not pursued further (consistent with `#730`'s own §9 "background/textbook
tail" treatment):

- **Sánchez, J., Net, M. & Simó, C. [2010], "Computation of invariant tori by
  Newton-Krylov methods in large-scale dissipative systems," *Physica D*
  239(3-4):123-133.** Method-adjacent (Newton-Krylov large-system torus
  solver) but scoped to **dissipative**, not Hamiltonian, systems — orthogonal
  to this project's conservative CR3BP/CCR4BP/CRNBP focus. Low priority.
- **Schilder, F., Vogt, W., Schreiber, S. & Osinga, H. M. [2006], "Fourier
  methods for quasi-periodic oscillations," *IJNME* 67(5):629-671.** A
  companion/extension of the already-flagged Schilder, Osinga & Vogt 2005
  (`#730`/`#732`'s existing PDE(CD) gap) — same author lineage, same class
  of gap, not independently re-flagged as distinct.

No genuinely new, high-priority, directly-on-point astrodynamics gap was
found in this book's own citation list beyond what `#730`/`#728` already
catalogued — consistent with this book being a pure-mathematics monograph
whose astrodynamics-specific worked example (Ch. 2's Earth-Moon L1 center
manifold) draws on the same already-well-covered KLMR/Gómez/Jorba/Masdemont
lineage extensively documented elsewhere in this project's corpus.

---

## 3. Summary for the coordinating session

Both papers filed (`cyclers_pdf` commit `1d796f3`) + digested + citation-mined
per `[[feedback_corpus_document_policy]]`. Direct answers to the two mandatory
cross-check questions:

1. **Does `genome/qp_tori.py`'s implementation faithfully match Olikara's
   2016 thesis method? NO — specifically and sourcedly no.** `qp_tori.py`
   implements the *earlier* Olikara & Scheeres 2010/2012 GMOS **shooting**
   method (full-period nonlinear propagation of each Fourier-mode sample +
   FFT-matched residual), which the thesis's own §1.2.2/§1.3 text explicitly
   identifies as its own predecessor and explicitly states it *replaces* with
   a Gauss-Legendre **collocation** scheme (§2.3.2, Eq. 2.31-2.34) — the
   thesis's actual, novel contribution. The seed formula (thesis Eq. 2.3) and
   invariance condition (thesis Eq. 1.4/2.1/2.29) match `qp_tori.py` exactly
   in form, confirming the GMOS-lineage identification precisely rather than
   loosely; the collocation machinery (piecewise-polynomial mesh, sparse
   block Jacobian, free Floquet-stability byproduct) has no counterpart
   anywhere in this project's code. This reinforces and sharpens (does not
   contradict) `#732`'s own independent Baresi/Olikara/Scheeres-2018-based
   finding that `qp_tori.py` is GMOS-lineage. Not urgent to fix — GMOS
   shooting is credibly the faster/more-accurate choice for this project's
   actual (low-amplitude, mildly unstable) use cases per Baresi et al. 2018's
   own benchmark — but the module's docstring citing "Olikara 2016" as
   co-equal lineage alongside "Olikara & Scheeres 2010" is imprecise: only
   the latter's method is actually implemented.
2. **Haro et al. 2016 cross-check against the Kumar-lineage papers: no
   discrepancy found — confirmed consistent.** The book's own Ch. 3
   (whiskered tori of quasi-periodic maps) and Ch. 5 (Newton-like NHIT
   method) are built on exactly the citation lineage (`HdlL06b/06c/07`,
   `dlLGJV05`) that `#730`/`#728` already independently traced as the
   theoretical ancestry of the already-acquired Kumar, Anderson & de la Llave
   2021/2022 papers — now confirmed directly from the book's own text. One
   structural nuance recorded: the book's own Ch. 5 Newton-like method is
   scoped to diffeomorphisms, treating the continuous-flow case only via a
   Poincaré-map reduction — a simplification the book itself states, not a
   project-specific gap.

Citation-mining across both surfaced five new gaps beyond what `#730`/`#732`
already flagged, all from the Olikara thesis's own literature survey:
**Calleja, Doedel, Humphries, Lemus-Rodríguez & Oldeman 2012** and **Arona &
Masdemont 2007** (medium-high priority — both directly on-point for this
project's own torus-to-torus heteroclinic-connection work); **Canalias,
Delshams, Masdemont & Roldán 2006** and its companion **Delshams, Masdemont &
Roldán 2008** (medium priority, scattering-map connection-geometry technique);
and **Jorba & Olmedo 2009** (medium priority, combined torus+stability method
for non-autonomous systems, relevant to this project's own currently-missing
CCR4BP/CRNBP torus stability output). None acquired here, per this task's
scope (digest + flag, not acquire).
