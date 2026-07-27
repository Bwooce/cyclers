# Digest: Kumar, Anderson & de la Llave 2021 (Celestial Mechanics and Dynamical Astronomy)

**Paper:** "Rapid and Accurate Methods for Computing Whiskered Tori and their Manifolds in
Periodically Perturbed Planar Circular Restricted 3-Body Problems"
**Venue:** Celestial Mechanics and Dynamical Astronomy, DOI `10.1007/s10569-021-10057-1`.
**Preprint:** arXiv:2105.11100v4 (5 Oct 2021 revision).
**Authors:** Bhanu Kumar (Georgia Tech), Rodney L. Anderson (JPL/Caltech), Rafael de la Llave
(Georgia Tech).
**Filed:** `kumar-anderson-delallave-2021-whiskered-tori-manifolds-cmda-arxiv-2105.11100.pdf`
(private `cyclers_pdf` repo).
**Acquired/digested:** 2026-07-27 (`#728`). Native text-layer PDF (LaTeX, embedded Type-1 fonts
confirmed via `pdffonts`); `pdftotext -layout` extracts cleanly. No OCR needed.

**Why in corpus:** this is the METHOD paper underlying the already-corpused Kumar/Anderson/de la
Llave application papers. In particular, `qp_torus_manifold.py`'s own module docstring already
cites "Kumar/Anderson/de la Llave 2025" for the claim "most unstable periodic orbits persist as
whiskered tori" — that 2025 SIAM ADS paper (arXiv:2109.14814, already digested
`2026-07-03-digest-kumar-anderson-delallave-2025-whiskered-tori-connections.md`) explicitly cites
**this** paper as "their prior work" for the quasi-Newton torus/bundle solver and Fourier-Taylor
manifold parameterization it reuses (confirmed directly: this paper's own References list its
CNSNS companion piece — task #728's *second* paper, below — as "Kumar et al. 2021a", and the 2025
SIAM ADS paper's method section credits this CMDA paper's quasi-Newton scheme as "ref 25"). So this
paper is the foundational algorithmic source the corpus had been citing one hop removed without
ever actually holding.

## 1. What the paper actually presents

A **method paper**, not an application/discovery paper. It develops, in the **periodically
perturbed planar CRTBP (PCRTBP)** — demonstrated concretely on the planar elliptic RTBP (PERTBP,
i.e. Jupiter-Europa with the real orbital eccentricity as forcing) — a five-part toolchain:

1. **A quasi-Newton method that solves SIMULTANEOUSLY for an invariant torus (as a 1D invariant
   circle `K(θ)` of the stroboscopic map `F`) and its center/stable/unstable bundle directions**
   `P(θ)`, `Λ(θ)` (§4). The near-diagonal structure of `Λ(θ)` (Eq. 12) decouples the linear system
   at each correction step into scalar cohomological equations solved by FFT/fixed-point contraction
   (§4.5-4.6), giving `O(N log N)` time / `O(N)` storage per step — versus `O(N^3)` time / `O(N^2)`
   storage for prior single-shooting/collocation methods that solve for `K` alone and get stability
   info as an afterthought (§4.5, Remark 4; explicitly contrasted with Farrés et al. 2017 and
   Olikara 2016's collocation approach, §1, §4.5). This is the key novelty versus the pure-torus
   parameterization methods of Haro & de la Llave (2006) and Haro et al. (2016): adding a genuine
   **center bundle** (the symplectic-conjugate direction, §4.3, §4.9) for UNSTABLE (not just
   elliptic) resonant tori, which the corpus's own torus-relevance codepaths care about.
2. **Continuation by both perturbation strength `ε` and rotation number `ω`** (§4.9-4.10), including
   a "gap-jumping" procedure (§4.10, Fig. 1) for crossing rotation-number values where the
   quasi-Newton method locally diverges (small-divisor breakdown).
3. **Fourier-Taylor parameterizations of the 2D stable/unstable manifolds** `W(θ,s) = K(θ) +
   Σ_{k≥1} W_k(θ) s^k` solving the invariance equation `F(W(θ,s)) = W(θ+ω, λs)` (§5), computed
   order-by-order via automatic differentiation + jet transport (§5.2) — the SAME recursive
   parameterization-method machinery as the companion CNSNS paper (paper 2 below), but extended
   from fixed points (1D manifolds) to 1D tori (2D manifolds).
4. **A Levi-Civita-style regularization of the equations of motion, DERIVED FRESH for the PERTBP**
   (§6.2) — the standard PCRTBP Levi-Civita regularization relies on the Hamiltonian being constant
   along trajectories, which fails once the periodic perturbation is added; the paper works out the
   canonical-transformation + time-rescaling derivation from scratch for the extended-phase-space
   Hamiltonian (Eq. 100-107), needed because globalized torus manifolds were observed passing
   arbitrarily close to (or through) the Europa singularity (§6.2, Fig. 7).
5. **Mesh-based globalization/visualization** of the manifolds beyond their fundamental domain of
   convergence (§6.1), by repeatedly applying `F`/`F^-1` to fundamental-domain points.

## 2. Key numerical results (sourced, with citations)

- **Fundamental-domain improvement**: comparing the manifold's fundamental domain of validity
  (`E_tol = 1e-5`/`1e-6`) at degree `d=1` (linear/eigenvector approximation) vs `d=5`
  (Fourier-Taylor): across 5 tested 3:4/5:6 PERTBP torus manifolds, degree-5 domains were **50-200x
  larger** than degree-1 — linear domains `~1e-4` at best (generally `~1e-5`), degree-5 domains
  `~1e-3` to `0.01` (§5.4, p. 30).
- **Physical demonstration system**: Jupiter-Europa PERTBP, `mu` = Jupiter-Europa mass ratio,
  eccentricity `ε = 0.0094` (real Europa eccentricity, §4.12). Both 3:4 and 5:6 resonant tori
  computed and continued by `ω` (Figs. 3-4); solution tolerance `1e-7` on the invariance equations.
- **Large-perturbation robustness test**: the SAME Jupiter-Europa mass ratio continued to `ε =
  0.206` — larger than the Sun-Mercury eccentricity, "one of the most eccentric two-body orbits of
  any pair of large solar system bodies" (§4.12, Fig. 5) — method remained robust (`N=1024`,
  `Δε=0.0005` continuation step).
- **Runtime**: on a 2017-era quad-core i7 laptop (Julia implementation), full continuation to
  `ε=0.206` took ~230 s; continuation to the physical `ε=0.0094` took <10 s (§4.12, p. 24). Manifold
  parameterization (C implementation, GSL RK8(9) integrator) took <10 s for `s^5`-order expansions
  on the same hardware (§5.3); globalized-mesh computation (Julia, `N=1024`, `L=101`, `k_max=15`
  for the 3:4 manifold) took ~250 s (§6.2.2, p. 35).
- **Regularization validation**: manifold meshes computed via the new regularized equations matched
  those from the unregularized equations where both were computable, but the regularized equations
  produced no integrator-divergence warnings while the unregularized ones did (§6.2.2, p. 35),
  confirming correctness while resolving the near-Europa numerical failures.

All values above carry their own page/section citation; none are read off a plot (no digitization
was needed — every number quoted here is stated in the body text or a labeled table/equation).

## 3. Relation to the companion CNSNS paper (this task's second paper, arXiv:2109.14800)

This CMDA paper's own References list cites the CNSNS paper as "Kumar B, Anderson RL, de la Llave
R (2021a)" (p. 38) — confirming the two are a matched pair: the CNSNS paper (submitted first,
2019 IAC presentation, published early 2021) computes 1D stable/unstable manifolds of ORDINARY
hyperbolic **fixed points** of the stroboscopic map (i.e., plain PCRTBP resonant periodic orbits,
no perturbation) via the parameterization method + jet transport, then finds heteroclinic
connections via Poincaré-section bisection. This CMDA paper extends the SAME parameterization-method
and jet-transport toolchain from 0D fixed points / 1D manifolds to **1D invariant tori / 2D
manifolds** in a PERIODICALLY PERTURBED PCRTBP, adding the bundle-computation quasi-Newton method
(entirely absent from the CNSNS paper, which only needs monodromy-matrix eigenvectors since there
is no torus, only a periodic orbit) and the PERTBP-specific Levi-Civita regularization (also absent
from CNSNS, which needed no regularization for its Jupiter-Europa/Earth-Moon test cases). In short:
CNSNS solves the unperturbed 3-body resonant-orbit-manifold problem; this CMDA paper solves the
perturbed 3-body resonant-**torus**-manifold problem that generically replaces it once a periodic
forcing (third body, eccentricity, etc.) is switched on.

## 4. Relevance assessment to `cyclerfinder`'s own codebase

**Grep results** (`src/cyclerfinder/`):
- `parameteriz*`: present in many files (`genome/qp_tori.py`, `genome/qbcp_torus.py`,
  `genome/bcr4bp_torus.py`, `genome/cross_system_cycle.py`, `genome/da_hotm_backend.py`, others) —
  but these all use the **GMOS (Gomez-Mondelo-Olikara-Scheeres)** collocation-based invariant-circle
  representation (`genome/qp_tori.py`'s own module docstring: "Following Olikara-Scheeres 2010 and
  Olikara 2016 (Purdue PhD)..."), i.e. exactly the O(N^3)-per-step collocation method this paper
  explicitly sets out to improve on and contrasts itself against (§1, §4.5 Remark 4). No file in the
  codebase implements this paper's O(N log N) quasi-Newton simultaneous-bundle solve.
- `whiskered`: two files, `search/ccr4bp_whisker.py` and `genome/qp_torus_manifold.py`. Neither
  implements this paper's method; `ccr4bp_whisker.py` (`#691`) extracts one-period-STM
  manifold-tangent EIGENVECTORS via a segment-anchored discrete-QR/covariant-Lyapunov-vector (CLV)
  technique (Benettin 1980/Ginelli et al. 2007/Dieci-Van Vleck) — a LINEAR local-direction extractor,
  not a Fourier-Taylor high-order manifold parameterization. `genome/qp_torus_manifold.py`'s own
  docstring explicitly cites "Kumar/Anderson/de la Llave 2025" for the "most unstable periodic
  orbits persist as whiskered tori" fact but then computes manifold directions via plain one-period
  STM eigenvectors (linear local approximation) — precisely the accuracy limitation (§1: "these
  linear stable/unstable directions are directly used as approximate local stable/unstable
  manifolds... neglecting higher order terms and thus losing accuracy") that this paper's §5
  Fourier-Taylor parameterization was built to fix.
- `manifold_expansion` (as a literal phrase/identifier): no hits — no Fourier-Taylor or
  high-order-Taylor manifold expansion machinery exists anywhere in the codebase.
- `Fourier-Taylor` / `jet transport` / `quasi-Newton` (as phrases): zero hits anywhere in
  `src/cyclerfinder/`. The order-by-order automatic-differentiation + jet-transport machinery this
  paper (and its CNSNS companion) rely on for computing `E_k(θ)` (§5.2) is simply not present.
- `stroboscopic`: present in many `search/`/`genome/` files, confirming the codebase already reasons
  in stroboscopic-map terms for QBCP/BCR4BP/CCR4BP tori — the right conceptual framework for this
  paper's methods to plug into, if ever adopted.

**Verdict: grounding + a real, currently-unexploited reuse opportunity, not an active dependency.**
`cyclerfinder`'s torus infrastructure (GMOS collocation for the invariant circle, linear
STM-eigenvector local manifold directions) is methodologically DISTINCT from — and, per this
paper's own §1/§4.5 argument and its own quoted 50-200x fundamental-domain improvement (§3 above),
strictly less accurate/efficient far from the torus than — what this paper develops. There is no
code today implementing the quasi-Newton simultaneous-bundle solve or the Fourier-Taylor
high-order manifold expansion. The concrete, honest reuse opportunity: `genome/qp_torus_manifold.py`
currently returns ONLY a linear (STM-eigenvector) local stable/unstable direction per torus point,
which is exactly the limitation this paper's §5 method targets — a future task replacing that
linear direction with a genuine Fourier-Taylor expansion (as this paper does, or per its CNSNS
companion's simpler fixed-point case) would be a direct, sourced adoption of this paper's algorithm,
not a novel invention. No positive-control opportunity was found (this paper's own numerical results
are all in the periodically-perturbed PCRTBP/PERTBP, not any system currently gauntletted in the
catalogue). This is a documented "background + latent reuse opportunity" verdict, not "no relevance."

## 5. Citation-mining pass

Checked `docs/notes/CORPUS_INDEX.md` (full-text grep for author surnames) before flagging. **None
of the following foundational parameterization-method papers — which this paper (and its CNSNS
companion) directly build on and cite repeatedly — are in the corpus:**

- **Cabré, Fontich, de la Llave (2005)**, "The parameterization method for invariant manifolds III:
  overview and applications," J. Differential Equations 218(2):444-515 — foundational parameterization-
  method reference, cited by BOTH papers in this task.
- **Haro, Canadell, Figueras, Luque, Mondelo (2016)**, *The Parameterization Method for Invariant
  Manifolds: From Rigorous Results to Effective Computations*, Applied Mathematical Sciences vol.
  195, Springer — the standard reference book for this entire method family; cited repeatedly by
  both papers (and by the already-corpused 2025 SIAM ADS paper).
- **Haro & de la Llave (2006, 2007)**, "A parameterization method for the computation of invariant
  tori and their whiskers in quasi-periodic maps" (parts I numerical algorithms / II
  explorations/breakdown), Discrete & Continuous Dynamical Systems-B / SIAM J. Appl. Dyn. Syst. —
  this paper's O(N log N) method directly extends theirs (§7 Conclusion: "extends the O(N log N)
  method of Haro and de la Llave (2006) and Haro et al. (2016) to unstable tori with center
  directions").
- **Fontich, de la Llave, Sire (2009)**, "Construction of invariant whiskered tori by a
  parameterization method. Part I," J. Differential Equations 246(8):3136-3213 — source of the
  "vanishing lemma" this paper invokes (§4.5.2) and the a-posteriori KAM framework it uses for its
  convergence argument (§4.7).
- **Huguet, de la Llave, Sire (2012)**, "Computation of whiskered invariant tori and their
  associated manifolds: New fast algorithms," Discrete & Continuous Dynamical Systems-A
  32(4):1309-1353.
- **de la Llave, González, Jorba, Villanueva (2005)**, "KAM theory without action-angle variables,"
  Nonlinearity 18(2):855-895 — cited for the center-bundle/symplectic-conjugate concept (§4.3).
- **Capiński, Gidea, de la Llave (2016)**, "Arnold diffusion in the planar elliptic restricted
  three-body problem," Nonlinearity 30(1):329 — cited for the KAM persistence argument underlying
  §3's NHIM/torus existence claims.
- **Olikara (2016)**, PhD thesis (Purdue) — the collocation method this paper directly benchmarks
  against; note `genome/qp_tori.py` ALREADY relies on this exact thesis's method (via Olikara-Scheeres
  2010) without the thesis itself being in the private paper corpus — a pre-existing, independently
  motivated gap this digest pass simply reconfirms.
- **Zhang & de la Llave (2018)**, "Transition state theory with quasi-periodic forcing," CNSNS
  62:229-243 — closely related prior use of a similar Fourier-Taylor manifold algorithm in a
  lower-dimensional setting (§5, explicitly compared).

This is a CLUSTER gap, not isolated: essentially the entire foundational parameterization-method
literature that both papers in this task (and the already-corpused 2025 SIAM ADS and 2023 Acta
Astro papers) repeatedly cite as their method's theoretical basis is absent from the private corpus.
Flagging for the coordinating session's prioritization, not acquiring myself per this task's scope.
