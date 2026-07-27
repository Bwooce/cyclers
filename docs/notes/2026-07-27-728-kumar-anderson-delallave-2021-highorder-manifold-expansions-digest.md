# Digest: Kumar, Anderson & de la Llave 2021 (Communications in Nonlinear Science and Numerical Simulation)

**Paper:** "High-order resonant orbit manifold expansions for mission design in the planar circular
restricted 3-body problem"
**Venue:** Communications in Nonlinear Science and Numerical Simulation 97:105691, DOI
`10.1016/j.cnsns.2021.105691`. Also presented at the 70th International Astronautical Congress,
21-25 October 2019, Washington, D.C. (per Acknowledgements).
**Preprint:** arXiv:2109.14800v1 (30 Sep 2021).
**Authors:** Bhanu Kumar (Georgia Tech), Rodney L. Anderson (JPL/Caltech), Rafael de la Llave
(Georgia Tech).
**Filed:** `kumar-anderson-delallave-2021-highorder-resonant-manifold-expansions-cnsns-arxiv-2109.14800.pdf`
(private `cyclers_pdf` repo).
**Acquired/digested:** 2026-07-27 (`#728`). Native text-layer PDF (LaTeX, embedded Type-1 fonts
confirmed via `pdffonts`); `pdftotext -layout` extracts cleanly. No OCR needed.

**Why in corpus:** companion paper to `#728`'s other paper (CMDA, arXiv:2105.11100, digested
`2026-07-27-728-kumar-anderson-delallave-2021-whiskered-tori-manifolds-digest.md`, same-day). That
CMDA paper's own References list cites this CNSNS paper as "Kumar B, Anderson RL, de la Llave R
(2021a)" — this is the EARLIER, simpler-scope member of the pair: it develops the parameterization
method + jet-transport toolchain for ORDINARY (unperturbed) PCRTBP resonant-orbit manifolds, which
the CMDA paper then generalizes to periodically-perturbed PCRTBP tori.

## 1. What the paper actually presents

A **method + demonstration paper** (NOT a periodically-perturbed / CCR4BP / N-body paper — this is
plain, single-perturbation-free PCRTBP, i.e. the ordinary 3-body problem with no third forcing
body). Three-part toolchain:

1. **Melnikov-based persistence analysis** (§2) to determine, ANALYTICALLY (to leading order in
   `µ`), which Keplerian `n:m` resonant periodic orbits survive as `µ` is turned on from 0 to a
   physical mass ratio. Uses synodic Delaunay coordinates (§1.2) and proves the Melnikov function
   `M(g_i)` is odd, hence has zeros at `g_i = 0` and (empirically, "strong numerical evidence" but
   not analytically proven for general `n,m`, §2 p.5) `g_i = π/n` — one elliptic, one hyperbolic
   orbit per resonance, per the Poincaré-Birkhoff fixed point theorem.
2. **Numerical continuation** (§3) of the Melnikov-identified Keplerian orbits from `µ=0` to
   physical `µ` (MATLAB `fsolve` continuation on the periodic-orbit fixed-point equation, varying
   both the state and period `T_sc`; §3, Eq. 19-20). Successfully continued **100 different
   Keplerian resonant orbits to the Jupiter-Europa mass ratio** and **32 to the Earth-Moon mass
   ratio** (§3, p. 6).
3. **High-order Taylor parameterization of 1D stable/unstable manifolds of the fixed points** of
   the period map `F` (i.e. plain PCRTBP resonant periodic orbits — a fixed point of `F`, not a
   torus; §4), via the SAME order-by-order automatic-differentiation + jet-transport recursion later
   reused (and extended to tori) in the CMDA companion paper. Manifolds computed to **degree 25-50**
   Taylor polynomials (§4, §6) — one order of magnitude higher polynomial degree than the CMDA
   companion's degree-5 torus manifolds, since fixed-point (not torus) manifolds are cheaper per
   order.
4. **Heteroclinic-connection search** (§5) via a Poincaré-section (`y=0, x<0`) reduction: propagate
   the manifold's fundamental domain to the section, globalize by repeated forward/backward
   Poincaré mapping (`W_p(λs) = P_+(W_p(s))`, Eq. 37-38), then find intersections of the resulting
   stable/unstable curves via a segment-based bisection/refinement algorithm (§5, Fig. 6 documents a
   spurious-intersection false positive caught and removed by refinement — a "the current segment
   endpoint pair, once refined, may no longer actually intersect" discipline directly relevant to
   this project's own orbit-closure verification norms).

## 2. Key numerical results (sourced, with citations)

- **Mass ratios used**: Jupiter-Europa `µ_E = 2.5266448850435028e-5`; Earth-Moon `µ_M =
  1.2150584270571545e-2` (§3, p. 6).
- **Continuation yield**: 100 Keplerian resonant orbits continued to `µ_E`, 32 to `µ_M` (§3, p. 6).
- **Fundamental-domain improvement**: over 60 resonant-orbit stable manifolds tested, degree-25
  Taylor parameterizations gave fundamental domains (`E_tol=1e-5`/`1e-6`) "orders of magnitude"
  larger than degree-1 (linear/eigenvector) approximations — linear domains `~1e-4` at best
  (generally `~1e-5`), degree-25 domains typically `~0.1` or even `~1` (§4.4, p. 10). The Conclusion
  states this as roughly a **"1000x improvement"** in domain of accuracy (§7, p. 15).
- **Worked example — Jupiter-Europa 3:4 <-> 5:6 resonance transfer at Jacobi constant `C=3.0024`**
  (§6, Table 1, p. 13): initial conditions, periods `T_sc`, and monodromy eigenvalues `λ_s, λ_u` for
  both periodic orbits are given to full published precision, e.g. 5:6 orbit `x=-1.231240907544348,
  T_sc=38.328135171743014, λ_s=0.001256465177783, λ_u=795.8835769446018`; 3:4 orbit
  `x=-1.391929713356257, T_sc=25.338526603095760, λ_s=0.011341070996024, λ_u=88.175093899915780`.
  These are genuine SOURCED state-vector/period/eigenvalue values (published table, not a
  digitization) — a rare case of a paper in this method family printing full-precision IC data
  directly.
- **Heteroclinic connections found**: 6 candidate Poincaré-section segment intersections detected
  between the 3:4 unstable and 5:6 stable manifolds at `C=3.0024`; 3 confirmed as genuine after
  bisection refinement (3 were spurious/false positives from unrefined segment endpoints, §6, Fig.
  6-7). The 3 genuine connection points' coordinates and `s_s, s_u` parameter values are given in
  Table 2 (p. 13), e.g. connection 3: `x=-1.1110838, y=5.780044e-15, ẋ=-0.10187786, ẏ=0.14762036,
  s_s=14.24735921, s_u=-3874.28227`.
- **Runtime**: "a few minutes on a laptop for a given pair of resonances" (§7 Conclusion, p. 15); no
  more granular hardware/timing breakdown given (contrast with the CMDA companion's more detailed
  per-step timings).

All values above carry their own page/section/table citation; none are digitized off a figure — the
IC/eigenvalue/connection-point tables (Tables 1-2, p. 13) are printed to full precision in the body
text itself.

## 3. Relation to the companion CMDA paper (this task's other paper, arXiv:2105.11100)

See the companion digest (`2026-07-27-728-kumar-anderson-delallave-2021-whiskered-tori-manifolds-digest.md`,
§3) for the full comparison. In short: this CNSNS paper is the simpler, earlier-conceived member of
the pair — plain (unperturbed) PCRTBP, fixed-point manifolds (1D), no bundle/center-direction
solve needed (only ordinary monodromy-matrix stable/unstable eigenvectors, since there is no torus).
The CMDA paper generalizes every piece of this paper's machinery (parameterization method, automatic
differentiation, jet transport, fundamental-domain/globalization concepts) from 0D fixed points to
1D invariant tori in a PERIODICALLY PERTURBED PCRTBP, additionally requiring the quasi-Newton
bundle solver and PERTBP-specific Levi-Civita regularization that this paper does not need.

## 4. Relevance assessment to `cyclerfinder`'s own codebase

Grep results are identical in substance to the companion CMDA digest's §4 (same codebase, same
grep pass, run once for both papers of this task): `parameteriz*` hits exist throughout
`src/cyclerfinder/genome/` and `src/cyclerfinder/search/`, but implement the GMOS
(Gomez-Mondelo-Olikara-Scheeres) collocation method (Olikara-Scheeres 2010 / Olikara 2016 PhD) for
tori, not this paper's parameterization-method-plus-jet-transport high-order Taylor expansion for
fixed-point manifolds. `Fourier-Taylor`, `jet transport`, and `quasi-Newton` (as literal phrases)
return zero hits anywhere in `src/cyclerfinder/`. `manifold_expansion` returns zero hits.

This paper's SPECIFIC applicability (unperturbed PCRTBP resonant-orbit manifolds — as opposed to
the CMDA companion's periodically-perturbed torus case) is narrower for `cyclerfinder`'s purposes:
the codebase's own resonant-orbit search machinery (`search/halo_family_at_jacobi.py`,
`genome/qp_tori.py`, etc.) already operates in the plain PCRTBP/CR3BP setting this paper targets,
and its manifold computations there currently rely on the standard LINEAR
eigenvector-of-monodromy-matrix local approximation (the exact limitation this paper's abstract
calls out: "most methods currently used in mission design rely on using eigenvectors of the
linearized dynamics as local approximations... not accurate except very close to the base invariant
object"). No code implements this paper's degree-25/50 Taylor-polynomial manifold expansion or its
Melnikov-based resonant-orbit-persistence search (§2) — the latter is a genuinely different,
currently-absent SEARCH STRATEGY (analytically narrowing down candidate resonant orbits via a
zero of an integral, rather than the codebase's typical numerical sweep/continuation approach) that
could in principle accelerate `cyclerfinder`'s own resonant-orbit discovery work if ever adopted.

**Verdict: pure background grounding, with one narrow, concrete reuse idea flagged (Melnikov
pre-screening of candidate resonant orbits before numerical continuation) rather than an active
dependency.** No code in this repository currently implements or depends on this paper's specific
methods. The heteroclinic-connection bisection/refinement procedure (§5) is conceptually close to
(but numerically distinct from — 1D fixed-point manifolds on a 2D Poincaré section, vs this
project's typical full-state or torus-level closure residuals) this project's own connection-finding
code; no direct code-reuse opportunity beyond the general "verify refined intersections, don't trust
the first Poincaré-segment crossing" discipline the paper's own Fig. 6 example illustrates.

## 5. Citation-mining pass

Checked `docs/notes/CORPUS_INDEX.md` (full-text grep for author surnames) before flagging. This
paper's reference list overlaps heavily with its CMDA companion's (Cabré-Fontich-de la Llave 2005,
Haro et al. 2016 book — both already flagged in the companion digest, not re-duplicated here). Two
ADDITIONAL candidates from this paper's own reference list not already covered:

- **Mireles James & Murray (2017)**, "Chebyshev-Taylor parameterization of stable/unstable manifolds
  for periodic orbits: Implementation and applications," International Journal of Bifurcation and
  Chaos 27(14):1730050 — cited (§6, p. 14) as a closely related PRIOR/COMPETING high-order manifold
  parameterization method (2D Chebyshev-Taylor manifolds vs this paper's 1D-manifold-plus-Poincaré-
  section approach); not in corpus and worth a look if the corpus ever wants a second, independent
  high-order-manifold-parameterization reference for cross-checking Kumar/Anderson/de la Llave's own
  approach.
- **Pérez-Palau, Masdemont, Gómez (2015)**, "Tools to detect structures in dynamical systems using
  jet transport," Celestial Mechanics and Dynamical Astronomy 123(3):239-262, DOI
  `10.1007/s10569-015-9634-3` — the primary jet-transport reference this paper (and its CMDA
  companion) cite for the automatic-differentiation/jet-transport machinery itself (§4.2); not in
  corpus.

Anderson & Lo (2010, 2010 JAS), Vaquero & Howell (2014), and the Anderson/Campagnola/Koh/McElrath/
Woollands (2019) Europa Lander endgame study are also cited here but were already flagged (the
latter explicitly) as uncorpused in the `#727` digest of the Acta Astronautica 2023 Kumar et al.
paper — not re-flagged as new here, just cross-confirmed as recurring across three different
Kumar/Anderson/de la Llave papers now.
