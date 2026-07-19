# Litteri, Gil, Vasile, Rodriguez-Fernandez & Camacho 2026 — VAE generation of CR3BP periodic orbits (CMDA 138:25)

**Date filed/digested:** 2026-07-19 (task `#637`).
**Trigger:** `#608`'s generative-ML-seed-model POC cited this paper as an external validating
anchor, but only had the paywalled abstract + the group's earlier open conference version
(arXiv:2408.03691) at the time. This task acquired the actual journal version.

**Citation:**
> Litteri, W., Gil, A. F., Vasile, M., Rodriguez-Fernandez, V., and Camacho, D. (2026).
> "Generation of periodic orbits in the restricted three-body problem with a variational
> autoencoder." *Celestial Mechanics and Dynamical Astronomy*, 138:25.
> DOI: 10.1007/s10569-026-10299-x. Received 20 Jan 2026 / Accepted 14 May 2026 / Published
> online 4 June 2026. **Open Access, CC BY 4.0** (confirmed directly on the Springer article
> page — not paywalled as of this filing; #608's "paywalled" characterization is now stale,
> see cross-check section below).

Filed as `litteri-2026-generation-periodic-orbits-variational-autoencoder-cmda-138-25-doi-10.1007-s10569-026-10299-x.pdf`
(45 pp., native text layer, no OCR needed) in the private corpus.

**Companion file also filed this pass:** the earlier open conference version this project's
`#608` actually worked from, Gil, Litteri, Rodriguez-Fernandez, Camacho & Vasile (2024),
"Generative Design of Periodic Orbits in the Restricted Three-Body Problem," arXiv:2408.03691
(6-7 pp. conference paper, SAIS 2024) — filed as
`gil-2024-generative-design-periodic-orbits-restricted-threebody-arxiv-2408.03691.pdf`. This
was NOT previously in the corpus despite being `#608`'s actual working source; filing it closes
that gap and lets the specific numbers `#608` cited be checked against their real origin (below).
Note the author ORDER differs between the two papers: conference = Gil, Litteri,
Rodriguez-Fernandez, Camacho, Vasile; journal = Litteri, Gil, Vasile, Rodriguez-Fernandez,
Camacho (Litteri promoted to first author, Vasile moved up) — `#608`'s citation order matches
the JOURNAL paper correctly.

## Model and method (journal version)

Planar+3D CR3BP (Eq. 1-3, standard normalized form), Earth-Moon (`mu=1.21505856e-2`, Table 1).
Working dataset: NASA JPL's public periodic-orbit database
(`ssd.jpl.nasa.gov/tools/periodic_orbits.html`), `M=44,112` initial states across 40 orbital
families (Lyapunov/Halo/Axial/Vertical/Butterfly/Dragonfly/DPO/DRO/resonant, Table 2); 33 of
the 40 families retained (7 dropped for numerical robustness — Vertical/Axial at L2/L3,
Long-Period+Vertical at L4/L5, westward LPO). Each orbit integrated to a `N=300`-node time
series (Julia Vern7, 1e-13 tol), later downsampled to 100 nodes for the main experiments; 250
samples/family (200 train / 30 val / 20 test).

**Architecture**: 1D-CNN VAE (not the naive fully-connected VAE of the 2024 conference
precursor), symmetric encoder/decoder, ReLU, beta-VAE loss (Eq. 7-9) with reparameterization
trick. Latent dimension swept `l={2,3,6,12}`; `l=2` used for most reported results (worse
reconstruction loss than `l=6`, but chosen for visualization clarity — explicitly stated as
carrying over to higher dimensions). Compared against an InceptionTime-based alternative
architecture, which overfits/is unstable by comparison (reconstruction loss an order of
magnitude worse, 0.372 vs 0.0194).

**Orbit classification**: random forest (scikit-learn, 500 trees) on FFT-coefficient + period +
Jacobi features, 100% test precision except L1-Lyapunov (95%).

**Refinement**: damped Newton-Raphson multiple shooting (Eq. 13-17) with an analytic Jacobian
(citing Spreen 2021), tolerance 1e-6, `Nmax=50` (undamped, coarse-grid) or `Nmax=100` (damped,
dense-grid).

**Novelty metric (the paper's real methodological contribution vs. the 2024 conference
precursor)**: NOT a raw convergence-rate statistic. A weighted sum `D_nov = D1+D2+D3` (Eq. 19)
of (1) global Mahalanobis distance to the whole training distribution, (2) minimum Mahalanobis
distance to any individual family's distribution, (3) classifier-confusion (`1 - max class
probability`). The paper explicitly motivates this by first showing (Sect 2.5.2) that plain
similarity/convergence is "not sufficient to define a robust novelty measure" — direct
methodological evidence this final version deliberately supersedes a cruder novelty notion.

## Key quantitative results

- **Coarse-grid search** (100 uniformly-spaced latent points, undamped refinement, full
  300-node/100-node discretisation): **43/100 converged** to a physically meaningful orbit
  (Sect. 3.3.2). This is general convergence, NOT a novelty-gated count.
- **Dense-grid search** (150 points/latent-dim, 2D latent -> 22,500 points, damped refinement,
  100 nodes, `Nmax=100`): **10,420/22,500 (46.3%) converged** to a physical periodic orbit
  (Sect. 3.3.3). Also general convergence, not novelty-gated — most of these land inside known
  family basins; novelty is a SEPARATE post-hoc ranking (the 100 highest-`D_nov`-scoring
  converged orbits are highlighted, four labelled "a"-"d" as most novel and shown in detail,
  Fig. 24-25).
- **Manifold continuation in latent space** (Sect. 3.3.1, Appendix A) is shown to be
  qualitatively equivalent to natural-parameter continuation in physical space but without
  needing incremental stepping — sampling far from a known family manifold directly, refinement
  still converges (demonstrated on L3-Axial -> L4-Axial/North-Dragonfly transitions, and on a
  specific case, "orbit 38," recovered independently via both latent continuation and physical
  NPC with matching branches).
- **Out-of-distribution / novel-family discovery mechanism** (Appendix A): shows explicitly how
  varying discretisation node count (100 vs. 11 nodes) on the SAME generated latent sample can
  converge to geometrically different physical orbits — one matching a known family, one
  belonging to a genuinely new 3D family not present in the 40-family training set (orbits
  "A"-"F," IC table given, Table 5) — establishing these aren't refinement artifacts via
  independent Jacobi-constant continuation branches.
- Prior companion (Vasile et al. 2025, cited but not yet acquired by this corpus) reports a
  family-holdout rediscovery experiment: entire families removed from training were
  regenerated by latent-space exploration, when their underlying features overlapped the
  retained distribution.

## Cross-check against `#608`'s citation (this task's mandated step 4)

`#608`'s bullet states: *"...externally de-risked by Litteri, Gil, Vasile,
Rodriguez-Fernandez & Camacho, CMDA 138:25 (paywalled; worked from the abstract + the group's
earlier open conference version, arXiv:2408.03691, which gives the real architecture: CNN-VAE,
2D latent, 44,112 ICs across 40 Earth-Moon families each as a 100-node time-series, refined by
multiple-shooting, 46% of 100 latent samples converged to a genuinely new orbit)."*

Checked directly against BOTH papers now in hand:

1. **The "46% of 100... converged to a genuinely new orbit" figure is accurately sourced to the
   2024 CONFERENCE paper, not the journal paper** — confirmed verbatim: "*46 out of the 100
   generated orbits were sufficiently accurate guesses for the refinement algorithm to
   successfully converge... Interestingly, all of the final refined orbits were new, meaning
   none of them were present in the training data*" (arXiv:2408.03691, Sect. 3.3). `#608`'s
   attribution to the conference version, not the journal, was correct and is now independently
   verified against the actual source rather than an inherited paraphrase — this is a genuine
   confirmation, not just a pass-through.
2. **However, the JOURNAL paper's own novelty definition is materially stricter than the
   conference paper's, and this nuance was not carried into `#608`'s summary.** The conference
   paper's "all refined orbits were new" claim uses the weakest possible novelty bar
   (not-bit-identical-to-a-training-sample — true of virtually any refined continuous-latent
   sample). The journal paper explicitly identifies this style of claim as methodologically
   insufficient (Sect. 2.5.2, quoted above) and replaces it with the multi-component
   Mahalanobis/classification-confusion `D_nov` metric, under which "novel" becomes a
   comparative RANKING (top-100 of 10,420 converged orbits) rather than a coarse pass/fail —
   and the journal paper never reports a "% converged that are novel" statistic at all. So
   `#608`'s "46%... genuinely new orbit" number, while correctly attributed to its actual
   conference-paper source, characterizes an EARLIER, methodologically superseded definition of
   novelty that the authors' own peer-reviewed follow-up explicitly moved past. Framed as "46%
   converged to a genuinely new orbit" without that caveat somewhat overstates how rigorously
   "novel" was being defined in that specific number.
3. **The paywalled characterization is now stale**: the journal version is Open Access (CC BY
   4.0), confirmed directly on Springer's article page — `#608` was correct that it was
   paywalled/inaccessible AT THE TIME (2026-07-16), but that is no longer the live status as of
   this filing.
4. **The `44,112`/`40 families`/`100-node` numbers are confirmed accurate and consistent across
   BOTH the conference and journal papers** — not conference-only as `#608`'s phrasing might
   suggest, but this is a non-issue since both papers report identical dataset provenance
   (same NASA JPL source).
5. **No correction to `#608`'s underlying engineering decision is warranted.** `#608`'s own
   bounded numpy/scipy POC (PCA + k-means + per-cluster Gaussian, operating on `(state0,
   period, Jacobi)` rather than the 100-node time series) remains a reasonable "linear-Gaussian
   analog" characterization independent of this correction — the journal paper's architecture
   (1D-CNN VAE) and dataset scale match what `#608` already described.

**Net verdict**: `#608`'s specific "46%" citation is verified GENUINE and correctly sourced
(not fabricated, not misattributed) — a positive outcome for `[[feedback_ground_citations_
against_content]]`. The one substantive correction is a precision/nuance gap (the number
reflects a superseded, looser novelty definition than the one the same authors later adopted in
peer review), plus the now-stale "paywalled" status. Both are minor annotation-level fixes to
`#608`'s bullet, not a retraction.

## Relevance to this project

- Confirms `#608`'s POC is philosophically aligned with a real, peer-reviewed precedent: a
  learned generative model over CR3BP periodic-orbit ICs, refined by a Newton/multiple-shooting
  corrector, DOES measurably increase physically-valid convergence and can cross basin
  boundaries into orbits (or novel families) not directly present in the training set — the
  same qualitative claim `#608` made for this project's own 54,165-orbit corpus (49% vs. 4%
  uniform-baseline convergence).
- The journal paper's rigorous `D_nov` novelty metric (Mahalanobis-distance-based,
  family-relative + global + classification-confusion) is a stronger, directly reusable design
  than `#608`'s own physical-sanity-filter-based novelty proxy, and a natural candidate
  enhancement if `#608`'s generative POC (or its `#614` follow-up) is ever advanced past
  bounded-POC status.
- Confirms the conference-paper novelty claim ("all refined orbits were new") is real and
  correctly attributed — not a paraphrase artifact.

**Status: digest complete.** No catalogue impact (methods/ML reference, not a specific citable
cycler), no code changes. `CorpusAnchor` registration not applicable (this is a background-ML
methods paper, same class as the existing `## ML / surrogate / GNC background` section of
`CORPUS_INDEX.md`, not a `literature_check.py` novelty-gate anchor).
