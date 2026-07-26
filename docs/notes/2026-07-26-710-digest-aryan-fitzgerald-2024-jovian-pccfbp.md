# #710 — Digest: Aryan & Fitzgerald 2024, AAS 24-103, "Four Body Invariant Structures And Chaos Analysis For Jovian Multi-Moon Ballistic Transfers"

**Citation.** Suryansh Aryan, Riley M. Fitzgerald, "Four Body Invariant Structures And Chaos
Analysis For Jovian Multi-Moon Ballistic Transfers," (Preprint) AAS 24-103, Virginia Tech Dept. of
Aerospace and Ocean Engineering. No DOI (confirmed — an unpublished AAS conference preprint, not
later journal-published). ResearchGate lists the conference as "AAS Astrodynamics Specialist
Conference 2024, Broomfield, CO, August 2024" (per its own metadata page) — this differs from the
`#710` dispatch bullet's "34th AAS/AIAA Space Flight Mechanics Meeting, Orlando FL, January 2024"
attribution; flagging the discrepancy rather than silently picking one (the in-paper preprint header
itself only says "AAS 24-103" with no venue/date, so this is unresolved from primary evidence in
this document alone — worth a one-line correction in `data/OUTSTANDING.md` if the coordinating
session wants to chase it, not done here per the "don't edit OUTSTANDING.md" instruction).

**Acquisition.** No DOI, no arXiv posting, no author-hosted preprint found (Virginia Tech faculty
page for Riley Fitzgerald and VTechWorks both checked, no hit). ResearchGate
(`researchgate.net/publication/383155359`) hosts the abstract and a "Public Full-text" author-
uploaded copy, but its own download link is Cloudflare-challenge-gated and never resolved after two
independent wait cycles (confirmed via live browser: "Verification successful. Waiting for
www.researchgate.net to respond" hangs indefinitely) — consistent with `#693`/`#700`'s prior
"ResearchGate 403" findings. Read the full inline "Public Full-text" preview directly in-browser
(all 22 pages, scrolled and read page-by-page) to confirm scope and extract the abstract + model
formulation + Table 1/2 numeric values before the user supplied the actual manuscript PDF directly
(`AAS_Paper_manuscript_final.pdf`, confirmed page 1 = "(Preprint) AAS 24-103" + matching title/
authors/affiliation). That user-supplied file is the one filed and digested here — a native
LaTeX/pdfTeX PDF (`pdftotext` extracts cleanly, 9,479 words, 22 pages; no OCR needed).

**Filed** at `cyclers_pdf/papers/aryan-fitzgerald-2024-four-body-invariant-structures-jovian-multi-moon-ballistic-transfers-AAS-24-103.pdf`
(private corpus repo, commit `72d5e21`).

## 1. The PCCFBP / PCC4BP model formulation

Same model class the paper calls interchangeably "PCCFBP" (abstract) and "PCC4BP" (body text) —
the Planar Concentric Circular Restricted Four-Body Problem. Three large bodies: Jupiter (`M`),
a "second primary" moon `m2` (Europa or Callisto, one model per system), and Ganymede (`m3`) as
the periodic PERTURBER (not a co-equal primary) on a fixed-frequency circular orbit. Massless
spacecraft, planar (`x`-`y`), all four bodies coplanar.

- Rotating (synodic) frame is the ordinary Jupiter–(second primary) CR3BP frame: `+x` from Jupiter
  to the second primary, origin at the Jupiter–second-primary barycenter. Ganymede's position is
  `r31 = (-mu1 + |r31| cos(Theta3), |r31| sin(Theta3))` with perturber phase
  `Theta3(t) = Theta3,0 + (Omega3 - Omega2) t` — i.e. Ganymede's synodic angle advances at the
  DIFFERENCE of its own and the second primary's mean motions, exactly the same "outer-moon
  synodic-rate" construction this project's own `core/ccr4bp.py` uses (`omega_gan = n3/n2 - 1`).
- Mass ratios: `mu1 = m2/(M+m2)` (the base CR3BP ratio for Europa or Callisto), `mu2 = m3/(M+m2)`
  (Ganymede's mass expressed relative to the Jupiter+second-primary system, NOT Ganymede's own
  `m3/(M+m3)` CR3BP ratio — a deliberately different normalization from this project's own
  `mu_gan = m3/(m1+m2)` convention in `core/ccr4bp.py`, which is actually the SAME quantity under a
  relabeled `m1`/`m2`, confirming the two projects independently arrived at an identical
  normalization choice).
- Time-dependent Hamiltonian (their Eq. 2, canonical `x,y,px,py,Theta3` — note `px = ẋ - y`,
  `py = ẏ + x`, standard CR3BP convention):
  `H = (px²+py²)/2 + px·y - py·x - (1-mu1)/r1 - mu1/r2 - mu2/r3 + (mu2/r31²)(x·cosΘ3 + y·sinΘ3)`.
- Non-canonical EOM (their Eq. 3, the one directly comparable to this project's `ccr4bp_eom`):
  standard CR3BP gravity terms from Jupiter and the second primary, PLUS a Ganymede term
  `-mu2(x-x3)/r3³ - mu2·cosΘ3/r31²` (and the `y` analogue) — i.e. direct + INDIRECT Ganymede
  acceleration, matching this project's own documented "direct + indirect Sun/perturber term"
  structure in `core/bcr4bp.py`/`core/ccr4bp.py`. `Θ̇3 = Ω3 - Ω2` is carried as a 5th phase-tracking
  state (their state vector is `{x,y,vx,vy,Θ3}`, 5D augmented — same shape as this project's own
  "theta1 tracked as an explicit phase" torus-corrector convention).
- **Numeric mass parameters actually used** (their Tables 1–2, both at `Θ3,0 = 0`):

  | System | `mu1` (2nd primary) | `mu2` (Ganymede) | Jacobi const. `C` |
  |---|---|---|---|
  | Jupiter–Europa–Ganymede | `2.5263e-05` | `7.7893e-05` | `3.0034` |
  | Jupiter–Callisto–Ganymede | `5.6623e-05` | `7.7890e-05` | `3.0044` |

  Note `mu2` (Ganymede's ratio) is ALMOST but not exactly identical between the two tables
  (`7.7893e-05` vs `7.7890e-05`) — a small residual from the different `(M+m2)` normalizer between
  the Europa- and Callisto-based systems, not a typo (the paper's own Figure 5 continuation study
  independently confirms `mu2 = 7.7892769e-05` as the "full" Ganymede value used for the
  Europa-Ganymede case specifically).

## 2. Torus computation method

Gauss-Legendre COLLOCATION (not GMOS, not the shooting-method Kumar 2021 uses) representing the
torus as a piecewise Lagrange-polynomial-in-time cylinder function `u(τ,θ2)`, discretized in the
free torus angle via a Discrete Fourier Transform (`N2` Fourier modes) and in time via `N` collocation
intervals of polynomial degree `m` (their Eq. 9–14). Two torus angles: `θ1 ≡ Θ3` (LOCKED to
Ganymede's forcing clock, not free — the same "theta1 locked, theta2 free" asymmetry this project's
own `search/variational_ccr4bp_torus.py` docstring documents explicitly as inherited from `#617`'s
QBCP corrector) and `θ2` (the free internal/orbital angle, rotation number `ρ = ω2/ω1`). An
"unfolding parameter" `Λ` augments the vector field (their Eq. 4) to build a well-posed BVP — this
is functionally the KAM-theory gauge/anchor-row machinery, though the paper frames it via Jorba's
normal-form literature rather than this project's own "three gauge rows (phase, amplitude,
rotation-number pin)" framing; same underlying purpose (removing solution non-uniqueness from the
non-autonomous system's lack of a global integral of motion).

Continuation: pseudo-arclength stepping (`Δs = 1e-8`) in TWO directions — the usual family
continuation (constant `I1`, varying rotation number `ρ`) AND a perturbation-strength continuation
directly varying `mu2` from `0` (pure CR3BP Lyapunov orbit) up to the full physical Ganymede value,
used explicitly to verify the family connects smoothly back to the known CR3BP limit (their Figure
5). Manifold/Floquet extraction reuses the SAME collocation Jacobian (no separate STM integration
needed) — their Eq. 15–18, standard Floquet-multiplier eigenspace of the monodromy-equivalent
matrix, propagated via 6th-order Runge-Kutta for the actual manifold trajectories.

## 3. Torus family results (Jupiter-Europa-Ganymede vs Jupiter-Callisto-Ganymede)

Both systems: 2:1 resonant CR3BP planar Lyapunov periodic orbit as the initial guess (`C=3.0034`
Europa case, `C=3.0044` Callisto case, both `ω=0.5`), converged L1 and L2 quasi-periodic tori,
accuracy `1e-6` to `1e-12`. Synodic period `T = 2π/|Θ3−Θ2| ≈ 7.0517 days` for the Europa case (the
Ganymede-relative-to-Europa synodic period). Full parameter table:

| Param | Europa L1 | Europa L2 | Callisto L1 | Callisto L2 |
|---|---|---|---|---|
| `ρ` (rotation no.) | 0.1595 | 0.20863 | 0.8964 | -4.0310 |
| `Λ` | 1e-10 | 1e-10 | 2e-09 | 2e-09 |
| `I1` (torus action) | 0.0027 | 0.004 | 0.0091 | 0.0156 |
| `N2` (Fourier modes) | 41 | 35 | 30 | 30 |
| `N` (collocation intervals) | 50 | 50 | 50 | 50 |
| `m` (polynomial degree) | 8 | 8 | 8 | 8 |

L1 family computed at constant `I1=0.02`; L2 family at constant `I1=0.0042` (zoomed). Continuation
in `μ2 = {0, 3.89464e-5, 7.7892769e-5}` shows invariant circles shrink toward the libration point as
Ganymede's perturbation strength increases (system stays bounded). Continuation in Ganymede's
initial phase `Θ3,0 = {0, π/2, π, 3π/2}` shows circles growing for `Θ3,0 ∈ [0,π/2)` then shrinking
back for `[π/2, 3π/2)` — four extrema total, global max in `[0,π/2)`, global min in `[3π/2,2π)`,
persisting across `I1` values. Largest unstable Floquet multiplier DECREASES monotonically with
increasing `μ2` (Ganymede has a net STABILIZING effect on the manifold direction) — this is a
genuine, reproducible-in-principle numeric trend claim, not just a qualitative statement.

## 4. Manifold-transit / heteroclinic connection results — the second positive-control candidate

This is the most load-bearing section for this project's own `#694`-built manifold-globalization +
heteroclinic-search pipeline (`search/ccr4bp_manifold_globalize.py`,
`search/ccr4bp_heteroclinic_search.py`).

- **Within-system L1↔L2 connections (CR3BP vs PCCFBP comparison, their Fig. 9–10).** For matched
  initial Hamiltonian values (to `1e-6`), the Europa system LOSES previously-existing CRTBP
  L1-unstable→L2-stable heteroclinic connections once Ganymede's perturbation is added (connections
  that existed at low energy in CR3BP cease to exist in PCCFBP). The Callisto system shows the
  OPPOSITE effect: connections that did NOT exist in CR3BP appear once Ganymede is added. Both
  effects are attributed to a shift in the Poincaré-map intersection region (`y = 1-μ2` plane) away
  from collision with the moon surface.
- **Cross-system Europa↔Callisto transit via a Ganymede flyby (their Fig. 11–13) — this DIRECTLY
  answers `#700`'s open question and CONFIRMS its finding rather than correcting it.** A transit
  trajectory was found from the Callisto L1 UNSTABLE manifold to the Europa L2 STABLE manifold, via
  an intermediate Ganymede flyby — i.e. a genuine two-hop, Ganymede-mediated chain, NOT a direct
  Europa-Callisto pair skipping Ganymede. Matched initial energy `HC ≈ HE ≈ -1.505`. Transfer
  times: **~60 days** for the Callisto L1 unstable manifold to reach Ganymede rendezvous, **~74
  days** for the L2 stable manifold from Ganymede to reach Europa — both figures match `#700`'s own
  reported numbers exactly (that prior pass inferred them from web-search summaries without reading
  the source directly; this digest is the first direct-from-source confirmation). Resonance
  structure: the Callisto-side unstable manifold transits a 2:3 MMR with Ganymede before the flyby,
  switching to 5:4 MMR in Ganymede's interior realm afterward (and continues to "skim" Europa's
  orbit for years with multiple encounters); the Europa-side stable manifold (backward-propagated)
  starts in 5:4 MMR with Ganymede in the interior realm, transitions to a non-resonant orbit
  escaping toward Callisto's exterior realm. Post-flyby energy: `HE ≈ -1.644` (Europa-side, more
  bound) vs `HC ≈ -1.75` (Callisto-side, slingshot to higher-potential exterior). A common-frame
  Poincaré map (Ganymede synodic frame, `y=0` plane) shows explicit phase-space intersections
  between the Europa-L2 and Callisto-L1 manifold bundles near Ganymede's orbit (`x ≈ 0.999922`) for
  specific eigendirections (`λ_stable = 0.0649`, `λ_unstable = 10.7185`), used to mark two candidate
  pre-/post-flyby connection states.
- **This does NOT reopen `#693`'s disqualification of the Ganymede-Callisto CCR4BP candidate** —
  that verdict rested on this paper covering the SAME model class for the SAME pair (Jupiter-
  Callisto-Ganymede, a genuine base+perturber system, not a "skip-Ganymede" pair), which this direct
  read confirms is exactly right: Table 2's `mu1=5.6623e-05` IS Callisto/(Jupiter+Callisto) with
  Ganymede as perturber — precisely `#693`'s disqualified candidate's model.
- **LCE (finite-time Lyapunov exponent) comparison (their Fig. 12c):** CRTBP and PCCFBP propagated
  manifolds track closely near the endpoint moons, diverge (LCE spikes) during the Ganymede
  interaction segment, and settle to a higher post-flyby divergence rate in PCCFBP than in CRTBP.
  Europa-side stable-manifold LCE magnitude is lower than Callisto-side unstable-manifold LCE, even
  though Europa's own leading Lyapunov exponent is larger than Callisto's — flagged in their own
  text as a slightly counterintuitive result, not further explained.

## 5. Direct comparison against this project's own CCR4BP pipeline (`#689`-`#694`)

Read `src/cyclerfinder/core/ccr4bp.py`, `src/cyclerfinder/search/variational_ccr4bp_torus.py`,
`src/cyclerfinder/search/ccr4bp_manifold_globalize.py`, and
`src/cyclerfinder/search/ccr4bp_heteroclinic_search.py` directly for this comparison (not from
memory).

**Model.** Structurally identical formulation: same EOM shape (base CR3BP + direct/indirect
Ganymede term), same synodic-frame convention (base = Jupiter–second-primary, Ganymede on a
concentric circle with linearly-advancing synodic phase), same mass-parameter normalization
(`mu2`/`mu_gan` both defined relative to `M + m_(second primary)`, not Ganymede's own barycenter).
This project's `core/ccr4bp.py` docstring cites the two Kumar papers as its model source, not Aryan
& Fitzgerald (unsurprising — Kumar 2021/2023 are this project's actual EOM derivation source and
predate this paper's 2024 date), but the EOM this paper derives independently (their Eq. 2–3,
citing Blazevski & Ocampo 2012 for the underlying Hamiltonian) is the SAME equation set. This is a
genuine third independent literature confirmation of the model (after Kumar 2021's own derivation
and this project's own from-scratch implementation) — worth noting as model-correctness triangulation
even though it wasn't this project's build-time source.

**Torus corrector.** Different numerical method (Gauss-Legendre collocation with DFT circle
discretization vs. this project's own real-tensor-product Fourier-in-both-angles pseudospectral
corrector, itself an "EOM swap" of the `#617` QBCP corrector) but the SAME underlying BVP structure
(non-autonomous, one locked angle + one free angle, phase + pseudo-arclength gauge constraints,
Newton solve). This project's corrector is MORE general in one respect (both torus angles are
Fourier-discretized; the paper collocates only in time and Fourier-discretizes only the free angle)
and the paper's is arguably simpler to implement (collocation avoids the 2D tensor-product Jacobian
this project's own corrector documents as its main complexity driver). Neither is validated against
the other in this pass — a genuine cross-implementation numeric comparison (same `mu1`/`mu2`/`C`
inputs, compare converged `ρ`/`I1`) would be a well-scoped, cheap follow-up task if a second
positive control beyond Kumar 2021 is wanted for the torus corrector specifically. NOTE: Kumar 2021
remains the stronger torus positive control for THIS project's own solver (`#690`) since Kumar
report actual physical resonant tori (3:4 Europa, 3:2/7:5 Ganymede) with directly reproducible
numbers already used in `#690`'s positive control; this paper's Table 1–2 values (rotation numbers,
`I1`) are a candidate SECOND, independent check but were not built against a resonance target the
way Kumar's were — same `C`, same 2:1 resonant Lyapunov seed, but the specific `ρ`/`I1` outputs are
this paper's own numerical result, not an independently-sourced target this project could grade
against without re-deriving their exact seed.

**Manifold globalization + heteroclinic search.** This is the closer methodological match to
`#694`'s own work, and the genuinely useful second positive control the coordinator flagged:
- Both use Floquet-multiplier eigenspace extraction off the converged torus to get manifold
  departure directions (this project: `ccr4bp_whisker.py`'s segmented-CLV eigenvector; paper: Eq.
  17-18's direct Floquet-matrix eigenvector from the collocation Jacobian — same underlying
  mathematical object, different extraction machinery).
- Both propagate the true nonlinear flow forward/backward from a small offset along the eigenvector
  to globalize a manifold tube (this project's `eps=1e-6` nondimensional vs. the paper's `ε` stated
  only as "exponentially small," not numerically specified — a genuine gap if an exact reproduction
  were attempted).
  Both explicitly justify keeping flow TIME continuous (not stroboscopically frozen) when searching
  for intersections — this project's `ccr4bp_manifold_globalize.py` docstring gives a full dimension-
  count argument for this; the paper's Poincaré-map approach (fixed `y=1-μ2` or `y=0` plane) is the
  same idea from the opposite direction (a codimension-1 SECTION rather than this project's own
  continuous-time 4-unknown least-squares match), i.e. the paper's method is closer to a classical
  Poincaré-section heteroclinic search while this project's `#694` module is a direct nonlinear
  least-squares tube-intersection refinement — different search STRATEGY, same underlying manifold
  objects.
- **The paper's genuinely reproducible target for a positive control**: the Callisto-L1-unstable →
  Ganymede-flyby → Europa-L2-stable transit (Section 4 above) is a concrete, numbered result (energy
  levels `HC≈HE≈-1.505`, transfer times 60d/74d, specific resonance transitions) that this project's
  own `#694` pipeline could in principle attempt to reproduce as a SECOND positive control beyond
  Kumar 2021 — but note this is a THREE-body-chain result (Callisto→Ganymede→Europa via an
  intermediate flyby), not a direct two-torus intersection the way `#694`'s own module is scoped
  (`ManifoldTube` pairs within ONE `CCR4BPSystem`); reproducing it would require chaining two
  separate `CCR4BPSystem` instances (Jupiter-Callisto-Ganymede then Jupiter-Europa-Ganymede) through
  a common Ganymede-frame transformation — a real scope expansion, not a drop-in rerun. Flagging as
  a candidate future task, not attempting it here (this task is acquire+digest only).

**Net assessment**: no prior-art conflict with anything already built (`#689`-`#694` all predate
this digest and were built from the Kumar papers, independently arriving at the same model); this
paper is a genuine, useful SECOND external validation source for the model equations and a
candidate second positive control for the manifold/heteroclinic search specifically, not yet
exploited.

## 6. Does this change the `#706`/`#708` Uranus Umbriel-Titania literature-clearance verdict?

**No — this paper has ZERO Uranian content and does not touch `#706`'s verdict.** `#706` cleared
`#701`'s Uranus Umbriel-Titania CCR4BP torus-homoclinic connection against a real-ephemeris-
survivable, epoch-recurring-connection claim specifically for the URANIAN system (Uranus + Umbriel +
Titania, with Titania as the perturber). This paper is exclusively Jovian — Jupiter + Europa/Callisto
+ Ganymede — and its own text never mentions Uranus, Umbriel, Titania, or any non-Jovian system
(confirmed by direct full-text read, not just its abstract). It IS the same PCCFBP/CCR4BP method
FAMILY applied to a different, disjoint physical system, so it is methodologically adjacent but not
prior art for the Uranian finding — the two results (this paper's Jovian tori/manifolds and `#701`'s
Uranian torus-homoclinic connection) are independent applications of the same technique to different
planets, exactly the relationship `#706` already anticipated when it distinguished "is there prior
CCR4BP/torus work on THIS pair" from "is the method itself novel" (the method is not novel — Kumar
2021/2023, McCarthy & Howell 2022 (BCR4BP), and now this paper are all prior CCR4BP/torus-manifold
work — but none of them touch Uranus). `#706`'s own six `KNOWN_CORPUS` Uranian anchors are unchanged
by this filing; this paper adds a Jovian corpus entry, not a Uranian one.

## 7. References worth noting (their bibliography)

- Blazevski & Ocampo 2012, *Physica D* 241(13) — the CCR4BP Hamiltonian source they cite (ref [9]),
  a DIFFERENT source paper from this project's own Kumar-derived EOM (worth checking in a future
  pass whether Blazevski & Ocampo's own formulation is already in this project's corpus — a quick
  grep of `CORPUS_INDEX.md` and `KNOWN_CORPUS` at digest time found no existing entry for it).
- McCarthy & Howell 2022, *Advances in Space Research* 71 — GMOS-algorithm BCR4BP quasi-periodic
  tori for ballistic lunar transfers, cited as the closest sibling four-body-torus paper (different
  system: cislunar BCR4BP, not Jovian CCR4BP) — not yet checked against this project's own
  `core/bcr4bp.py` corpus coverage in this pass.
- Olikara 2016 PhD thesis (Colorado Boulder) — the collocation-technique source (ref [12]), likely
  the same lineage as `search/variational_qbcp_torus.py`'s own collocation-adjacent literature if a
  future digest wants to trace it.
