# Mining digest — Ross & Roberts-Tsoukkas (2026), "Stable Families of Ballistic Prograde Cyclers in the RTBP"

**Citation:** Shane D. Ross & Michael Roberts-Tsoukkas, *"Stable Families of Ballistic Prograde
Cyclers in the Restricted Three-Body Problem,"* arXiv:2606.29189v1 [nlin.CD], 28 Jun 2026 (dated
June 30 2026). Virginia Tech, Aerospace & Ocean Engineering. PRL-format letter, 5 pp.
**Filed:** `cyclers_pdf/papers/ross-roberts-tsoukkas-2026-stable-ballistic-prograde-cyclers-rtbp-arxiv-2606.29189.pdf`
(text-layer, no OCR). **Mined for:** #494 (binary (k₁,k₂)-cycler μ-family). User-acquired 2026-06-30.

## What it is
The formal μ-family generalization of AAS 25-621 (ref [30]; the source of the catalogue's 5 EM
`ross-rt-em-cycler-*` rows, V2). Reports **stable, ballistic, prograde (k₁,k₂)-cyclers** — periodic
orbits that alternately undergo temporary capture about each primary — and demonstrates them across
**>2 orders of magnitude in mass ratio, μ=0.001 (Sun-Jupiter) → 0.5 (equal mass)**. Central claims:
1. Construction: cyclers lie in the intersection of the **stable/unstable manifold tubes of the L1
   Lyapunov orbit**; symmetric (ẋ=0) seeds from the discrete intersection Γ = P(S^{u1}_{k1}) ∩ S^{u2}_{k2},
   refined by **differential correction at fixed C** (perpendicular crossing ẋ=0 at half period), then
   **pseudo-arclength continuation in the (x₀, C) plane**. This is EXACTLY the repo's existing lane
   (`cr3bp_periodic.correct_symmetric_fixed_jacobi` + `ydot0_from_jacobi` + `mu_continuation`).
2. Every family is **born in a saddle-center bifurcation** of the return map at its maximal Jacobi
   constant C^bif_{(k1,k2)} < C(k1,k2), creating a planar-stable + a hyperbolic branch simultaneously.
   **Conjecture (universality):** saddle-center birth is universal ⇒ every cycler family has a stable
   subfamily.
3. Stability splits into **planar (sₚ) and vertical (sᵥ)** indices via the 6×6 spatial monodromy
   (sₚ = ½(λₚ+1/λₚ), sᵥ = ½(λᵥ+1/λᵥ)); stable ⇔ |sₚ|<1 AND |sᵥ|<1. Vertical instability arises only
   via isolated parametric resonance of the Hill equation z̈ + Ū_zz z = 0 (Ū_zz>0). Across 9 EM
   families (1≤k₁,k₂≤3), |sᵥ|>1 occurs almost only on the planar-UNSTABLE branch; every family has a
   doubly-stable subfamily.

## Conventions — IDENTICAL to the repo (golden directly usable)
- Primaries m1=1−μ at x=−μ, m2=μ at x=1−μ; μ∈(0,½]. Eq. 1 EOM = our `cr3bp_eom`.
- **Jacobi (Eq. 2): C = −2Ū − (ẋ²+ẏ²), Ū = −½(x²+y²) − (1−μ)/r1 − μ/r2.** Expand: −2Ū =
  (x²+y²)+2(1−μ)/r1+2μ/r2, so C = (x²+y²)+2(1−μ)/r1+2μ/r2 − v² — **byte-for-byte our
  `core/cr3bp.py::jacobi_constant`** (verified 2026-06-30).
- (k₁,k₂)-cycler: k₁ prograde circuits about m1, k₂ about m2 per period; k₁,k₂ = #crossings of U1⁻,
  U2⁺. Reversibility sₓ:(x,y,ẋ,ẏ,t)↦(x,−y,−ẋ,ẏ,−t). Symmetric IC: y₀=ẋ₀=0, ẏ₀ from C (= our
  `ydot0_from_jacobi`, Ross Eq. 9).
- L1-threshold energy: cyclers exist for C < C1 (transport through the L1 neck open).

## TABLE I — the sourced golden (verbatim, 15–16 digits)
Representative **fully-stable** symmetric (k₁,k₂)-cyclers. Each orbit: perpendicular x-axis crossing
x₀, with y₀=ẋ₀=0, ẏ₀ determined by C. T = period. All rows satisfy max(|sₚ|,|sᵥ|) < 1.

| μ | (k₁,k₂) | x₀ | C | T | sₚ | sᵥ |
|---|---|---|---|---|---|---|
| 0.001 | (1,1) | -0.647047499999966 | 3.031605708907296 | 14.774502790974823 | 0.4121 | -0.2943 |
| 0.012150584270572 | (1,1) | -0.768217354461248 | 3.151175879917331 | 10.291893641936499 | 0.8210 | 0.6358 |
| 0.012150584270572 | (3,3) | -0.322477620583087 | 3.183379082910527 | 19.503763587070285 | 0.9855 | 0.6207 |
| 0.1 | (3,2) | -0.694376003123377 | 3.573367616904619 | 12.295263874014290 | 0.5686 | 0.9175 |
| 0.3 | (3,1) | -0.804725783387797 | 3.701958166478617 | 9.094576400494693 | 0.0294 | 0.8307 |
| 0.5 | (1,1) | -0.519689929077496 | 3.628400000000000 | 8.792013561462247 | 0.9376 | 0.2130 |

(Stored as `data/golden/ross_rt_2026_cycler_families.yaml`.)

## #494 implications — acquisition risk RESOLVED, μ-gap goldens in hand
- The #494 design's flagged "honest risk" (ICs might be figures-only) is **resolved**: Table I
  tabulates `(μ, x₀, C, T)` directly. The recover-from-(μ,C,T) machinery that did the EM slice
  reproduces these directly (no digitization).
- **μ=0.1 (3,2)** is the sourced anchor adjacent to **Pluto-Charon (μ=0.10851)** — the Phase-3
  instantiation = continue the μ=0.1 (3,2) family to μ=0.10851 (or direct-search seeded at the
  μ=0.1 IC). a_bin=19,596 km, Charon GM≈106 (massive enough for the capture leg — the #489/#492
  small-moon objection does NOT apply to the (k₁,k₂) Charon-capture cycler).
- Pluto-Charon is **NOT in Table I or anywhere in the paper** → a Pluto-Charon (k₁,k₂)-cycler is a
  genuinely fresh REAL-SYSTEM instantiation (modest novelty, the #312-Uranus framing), not a
  reproduction of a tabulated row. The paper's Outlook states it as open: *"allowing the third body
  to possess non-negligible mass raises the possibility of stable ballistic planets that alternately
  orbit the members of a binary star system. Whether such cycling planets exist is an open
  question."* — the explicit #315 binary-star lead.
- The EM (1,1) Table-I row (C=3.151176, T=10.291894) is an **independent cross-check** of the
  catalogued `ross-rt-em-cycler-11` golden + the #494 Phase-0 positive control (the catalogued rows
  trace to AAS 25-621; agreement here is a second sourced confirmation, not circular).

## Stability-index mapping (for the gauntlet)
Paper's sₚ/sᵥ = ½(λ+1/λ) are the **full-period** planar/vertical indices. The repo's
`barden_stability` returns ν = ½(λ+1/λ) on the **half-period** STM (reciprocal-pair convention).
Map carefully: |sₚ|<1 ⇔ planar-elliptic; cross-check the repo's planar ν against Table I's sₚ at the
EM rows during Phase 0 (the vertical sᵥ needs the 6×6 spatial monodromy — the repo's planar-CR3BP ν
covers sₚ only; sᵥ requires the spatial variational integration, a Phase-2/3 add for the
doubly-stable verdict).

## Provenance / cross-refs
- ref [30] = AAS 25-621 (Ross & Roberts-Tsoukkas 2025) = catalogue `ross-rt-em-cycler-*` source.
- ref [36] = Braik & Ross 2026 orbital networks (arXiv:2605.31543) = filed/mined (#249 C11a/C11b/C21/C32 V1).
- refs [25],[26] = Koon-Lo-Marsden-Ross (filed: shoot-the-moon + KLMR book) — the manifold-tube basis.
- Distinct from the filed `roberts-tsoukkas-ross-2026-stable-prograde-em-cyclers-journal.pdf`
  (a MULTI-ORBITER EM paper) and the VSGC student summary — different papers, same author group.
