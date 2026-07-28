# Scoping `#752`: classical-CR3BP resonant-orbit invariant-manifold tour-design machinery — GO/NO-GO

**Task:** `#752` (research/scoping only — no code, no catalogue changes, no runs). Question: how
large a lift would Anderson/Lo-lineage classical planar-CR3BP resonant-orbit
invariant-manifold / homoclinic-heteroclinic machinery actually be, and is it worth reopening
capability-building (`[[project_capability_frontier_complete]]`) for?

**Sources read directly this pass** (not just digests): Anderson & Lo 2011 (JAS 58:167, full
method + continuation + Tables 1-3 sections of the text layer), Anderson & Lo 2010 (JGCD 33:1899,
via `#745`'s digest plus targeted text checks), Barrabés-Mondelo-Ollé 2009 (Nonlinearity 22:2901,
pp. 2903-2907 — the §2.2 zero-finding formulation, Eq. (3)-(6), and the multiple-shooting
augmentation). Project code read directly: `search/ccr4bp_heteroclinic_search.py` (all 637 lines),
`search/resonance_network.py`, `genome/heteroclinic_cycle.py`, `search/cr3bp_periodic.py`
(`correct_symmetric_fixed_jacobi`), `core/flyby.py`, `search/resonant_conic.py`,
`search/ieg_seed.py` (docstrings + structure), plus a very-thorough read-only sub-agent inventory
of the full `search/`+`genome/` CR3BP periodic-orbit/manifold stack, with every load-bearing
claim below re-verified first-hand against the files.

---

## Headline verdict: **GO — but as a small ASSEMBLY task, not a capability-building program.**

The single most important finding of this pass: **the premise that "cyclerfinder has NO
implementation of this machinery" (`#742`/`#745`/`#749`, three independent confirmations) is
substantially overstated.** All three digests grepped for "homoclinic"/"resonant" in the
Jovian-tour modules (`core/flyby.py`, `search/resonant_conic.py`, `search/ieg_seed.py` — which
are indeed patched-conic/Lambert throughout, confirmed again here) and missed that the exact
machinery class lives elsewhere in the tree under Track-B/Earth-Moon names:

- **`search/resonance_network.py` (`#267`, Kumar/Rawat/Rosengren/Ross 2025 cislunar
  resonant-transport reproduction)** already computes, in the planar CR3BP: unstable
  mean-motion-RESONANT periodic orbits at fixed Jacobi constant (3:1, 4:1, 2:1 Earth-Moon
  members, seeded from the paper's own Table 6 sourced ICs via
  `correct_symmetric_fixed_jacobi`), their full-period planar monodromy + largest
  Floquet eigenpair (`_planar_floquet`), manifold globalization by ε-offset along the
  eigenvector with forward/backward integration (`compute_floquet_manifold`), a
  perigee/perilune Poincaré section, and a manifold-overlap accessibility metric between two
  different resonant orbits' sections (`perigee_overlap`, `ResonanceNetworkScorer`) — with the
  standing Radau-vs-DOP853 independent-integrator cross-check discipline built in.
- **`genome/heteroclinic_cycle.py` (`#314`, validated against Wilczak-Zgliczyński's
  computer-assisted Sun-Jupiter L1↔L2 proof)** already has the genuine transversal
  **manifold-to-manifold Newton intersection solver**: `correct_connection` does a 2×2 Newton on
  the two manifold departure phases (τ_u, τ_s) driving the {y=0}-section gap
  Wu(A)∩Ws(B) → 0 between two different equal-energy unstable orbits, with STM transport of the
  Floquet eigenvector to arbitrary phase (`_seed_on_manifold`), k-th-crossing section logic
  (`_section_crossing`), coarse phase-scan seeding, and an independent-integrator
  `crosscheck_cycle`. Its docstring even says it "replicates the focused Floquet
  manifold-seeding pattern from search/resonance_network".
- `correct_symmetric_fixed_jacobi` (`search/cr3bp_periodic.py`) already takes
  `half_crossings: int | None` — exactly the "first x-axis crossing is not necessarily the
  desired intersection" modification Anderson & Lo 2011 describe as their own change to the
  Howell-Breakwell symmetric corrector for resonant orbits. Asymmetric general fixed-C and
  free-C correctors (`cr3bp_general_periodic{,_free_c}.py`, analytic STM Jacobians), natural
  Jacobi continuation (`cr3bp_continuation.py`), fold-crossing pseudo-arclength continuation
  (`cr3bp_jacobi_arclength.py`), multiple shooting with `monodromy`/`floquet_multipliers`
  helpers (`cr3bp_multiple_shooting.py`, `bifurcation_detector.py`), and seedless
  harmonic-balance enumeration (`{deflated_,}variational_periodic_orbit.py`) all exist.

What is genuinely missing is narrow: (i) no **planet-moon application** of any of it — the
resonant-orbit work is all Earth-Moon cislunar (`_RESONANT_SEEDS` is an Earth-Moon table); (ii)
no **Newton-corrected resonant-to-resonant connection point** — the resonant side has only the
overlap *score* (`perigee_overlap`) and lobe-area geometry (`lobe_overlap_scorer.py`), while the
Newton solver is wired to `LyapunovNode` libration orbits only; (iii) no **homoclinic**
(self-connection, A=B) mode anywhere (`#749`'s finding, still true); (iv) no **BMO-style
continuation of connection FAMILIES in energy** (the full Eq.-(6) system with eigenvector
conditions folded in — `#746` already established there is nothing to port to, and nothing here
changes that; it stays out of MVP scope).

The correction to the record matters beyond this task: `#745`'s digest sentence "no code here
computes a classical-CR3BP resonant-periodic-orbit family, its monodromy matrix, or its
Poincaré-section manifold intersections" is **wrong as literally stated** (resonance_network
does all three, with the intersection as an overlap metric rather than a corrected point), and
the `#751`/`#752` dispatch bullet inherited that framing.

## Q1 — Minimum viable version, stage by stage

Anderson & Lo 2011's pipeline, read from the paper directly: grid-search x-axis ICs near the
two-body p:q resonance → symmetric single-shooting corrector (k-th crossing) → linear
extrapolation continuation in C + secant to land on a target C → monodromy eigenvectors →
~1e-6 eigenvector offset → globalize by integration → one-sided {y=0} Poincaré section →
intersect the two manifolds' section curves → interpolate nearest points and integrate to get
the homoclinic/heteroclinic connection.

| Stage | Status in repo | Delta needed |
|---|---|---|
| (a) unstable resonant PO at fixed C, planar CR3BP | **EXISTS** (`correct_symmetric_fixed_jacobi` + `half_crossings`; `recover_resonant_family` pattern; natural + pseudo-arclength continuation; `barden_stability`/`floquet_multipliers`) | Jupiter-Europa (µ=2.5266448850435e-5) seed grid for the 3:4 / 5:6 exterior families — a small two-body-ellipse-IC sweep + corrector loop; the paper documents the exact procedure. Kumar 2021 (AAS 21-651, in corpus) can supply an independent 3:4 seed cross-check (`#750`'s lineage question, bonus synergy). |
| (b) manifold globalization off monodromy eigenvectors | **EXISTS twice** (`compute_floquet_manifold` for resonant POs; `_seed_on_manifold` + STM phase transport for Lyapunov POs) | Essentially none — parameter plumbing only. |
| (c) manifold-to-manifold intersection on a Poincaré section | **PARTIAL** — Newton solver exists (`correct_connection`) but only for `LyapunovNode`s on {y=0}; resonant side has overlap scores only | The one real piece of new code: retarget `correct_connection`'s node abstraction to resonant members (both node types already carry `state0`/`period`/`jacobi`/Floquet eigenpair), add the homoclinic A=B mode, keep the exact-Jacobi equality as a hard gate (autonomous system — *stronger* than the CCR4BP's approximate quasi-Jacobi check), port the ghost-guard discipline (trivial-near-departure exclusion + independent-integrator re-check, both already implemented in the two donor modules). |

**Digit-grade positive controls exist and are already in hand** — rare and decisive for this
project's gauntlet discipline (`feedback_verify_gauntlet_with_positive_control`,
`feedback_golden_tests_sourced_only`): Anderson & Lo 2011 states µ = 2.5266448850435e-5 and
C_flyby = 2.99163956830415 explicitly; Table 1 gives max monodromy eigenvalues per family
(3:4-LO 1036.116088, 5:6-LI 1.000008, 5:6-LO 4445.387515, 5:6-NO 28178.258323); Table 2 gives
the homoclinic intersection state (x=-1.28427733, y=0.0, ẋ=0.00000009, ẏ=0.46372205); Table 3
the heteroclinic one (x=-1.43029175, y=0.0, ẋ=0.00018678, ẏ=0.67262261). All read from the
text layer this pass. Both stage gates (family eigenvalues; connection states) are therefore
sourced-anchor testable, not self-referential.

## Q2 — Is `ccr4bp_heteroclinic_search.py`'s approach the thing to generalize? **No — wrong donor.**

Read in full. Its architecture (coarse KD-tree near-coincidence scan over discretized tubes →
continuous 4-unknown `least_squares` refine on (θ2_u, t_u, θ2_s, t_s) → ghost guard) *is*
conceptually transferable, and its guard discipline (trivial-departure exclusion,
independent-integrator consistency, the `#701`/`#702` seed-anchored `ref_vec` lesson) should be
ported as *discipline*. But its complexity exists to handle problems the autonomous CR3BP does
not have: segmented-CLV direction extraction with no sign-continuity (a torus has no monodromy
matrix — the classical PO does, so the eigenvector is globally well-defined per orbit and STM
transport is exact), a 2-phase × 2-time unknown space (the PO case reduces to the 2-unknown
(τ_u, τ_s) Newton at fixed section + fixed C that `correct_connection` already implements), and
an approximate quasi-Jacobi plausibility check standing in for a conserved quantity the CCR4BP
lacks (the CR3BP has the exact Jacobi constant — a hard equality gate). `#746`'s "architecturally
unrelated, nothing to port" verdict w.r.t. BMO continuation is confirmed and extends here:
the right donor for the resonant-connection MVP is `genome/heteroclinic_cycle.correct_connection`
(same model, same object class, already validated against a computer-assisted proof), with
`resonance_network` supplying the resonant nodes.

## Q3 — What it concretely unlocks

1. **A catalogue-eligible new object class: ballistic resonance-transition cyclers in
   planet-moon systems.** Anderson & Lo 2011's constructed trajectory is *periodic* — it cycles
   3:4↔5:6 via two Europa flybys indefinitely, ballistically, closed to 0.091 km / 5.3e-6 km/s —
   i.e. a genuine Jupiter-Europa cycler in the CR3BP, squarely in the catalogue's `cycler`/
   `quasi_cycler` classes. No such object class exists in the catalogue today (383 rows; grep
   confirms zero Jovian 3:4/5:6 resonant rows, `#745`).
2. **The novel-discovery play, matching all three prior novel writebacks' pattern** (published
   method → unstudied system: `#312` Uranus, CCR4BP Umbriel-Titania, N=5 CRNBP torus): the
   published record covers Jupiter-Europa (Anderson-Lo 2010/2011), Uranus-Oberon *heteroclinic
   connections* (Anderson-Kumar AAS 24-288, digested `#728` — PCRTBP heteroclinics at
   C=3.0028-3.0072, so plain Oberon connections are NOT novel), and Earth-Moon cislunar chains
   (Kumar 2025/2026). Unstudied, per the corpus: **Saturn-Titan, Saturn-Enceladus,
   Uranus-Ariel/Umbriel/Titania, Neptune-Triton (retrograde — `#599` capability), Pluto-Charon**
   resonance-transition periodic cyclers, and even at Jupiter-Europa the *ballistic
   resonance-cycling periodic orbit as a catalogued cycler family* (vs the paper's single
   worked example). Every candidate must still clear `search/literature_check.py`
   (`feedback_literature_novelty_check_baseline`) — the Oberon example shows the prior-art
   surface here is actively moving.
3. **Secondary unlocks**, real but not the justification: reproducing Casoliva Class 2
   (He1-4/Hm1-2 Earth-Moon homoclinic-shadowing cyclers — published, digested `#725`, zero
   catalogue rows today) becomes a homoclinic-mode application of the same solver; the
   `resonant_conic.py` EGGIE Stage-3 gap ("final ballistic closure left to the multiple-shooting
   corrector") gains the CR3BP-manifold seeding route Anderson & Lo 2011's own "Potential
   Techniques" section recommends; and `#750`'s two-hop lineage question gets its digit-grade
   3:4-family cross-check essentially free from stage (a)'s gate.

The honest counterweight: the MVP's own two positive-control gates *are* reproductions of
published results, and step 2's novelty is prospective — the same "expect census, novel hits
are rare" prior applies as everywhere else. But the cost side (below) is far below the
`[[feedback_speculative_high_effort_required]]` threshold where that debate matters.

## Q4 — Effort, calibrated against `#714`→`#736`

The N=5 CRNBP torus arc was ~23 task numbers over ~2-3 days wall-clock and required **new
dynamics** (N=5 EOM+STM `#717`, sign discrepancy resolved from source), **new object machinery**
(torus continuation `#720`), a phase-default bug cycle (`#721`/`#723`), and a multi-round
novelty fight (`#721`/`#722`/`#724`). This build needs **no new dynamics** (planar CR3BP core is
the oldest, best-validated model in the repo), **no new object class machinery** (PO monodromy
manifolds already exist in two validated modules), and has **digit-grade sourced gates at both
stages** (the N=5 arc had none and paid for it in adversarial-review rounds).

Estimate: **MVP = 2 dispatchable tasks, ~1-2 days wall-clock at this project's cadence —
clearly SMALLER than the `#714`→`#736` arc** (roughly a third to a half of it, including
verification overhead). A follow-on discovery campaign over new moon systems is separable and
incremental (per-system tasks, each gated by the same machinery + literature check). Main
technical risks: (i) extreme instability (max eigenvalue ~2.8e4 per period) makes long manifold
integrations sensitive — mitigated by the existing independent-integrator discipline and by
keeping globalization arcs short (the paper's own approach); (ii) Anderson & Lo publish no IC
tables for the resonant families, so stage (a)'s families must be regenerated by the documented
grid+continuation procedure with only Table-1 eigenvalues + family topology as anchors —
acceptable, and the Kumar 2021 3:4 table is an independent partial anchor.

## Recommendation: **GO.** Concrete first task (dispatchable as written)

**Task A — "Jupiter-Europa 3:4/5:6 unstable resonant families + Anderson-Lo Table-1 gate"**
(spec-complete, Sonnet-tier per `[[feedback_subagent_model_tiering]]`): in a new
`search/` module, build a resonant-orbit seed sweep for the Jupiter-Europa planar CR3BP at
µ = 2.5266448850435e-5 (two-body p:q ellipse x-axis ICs, p:q ∈ {3:4, 5:6}, both interior/outer
branch signs), converge members via the existing `correct_symmetric_fixed_jacobi`
(`half_crossings` set per resonance), continue each family in C via the existing natural/
pseudo-arclength tools to C = 2.99163956830415, and classify via the existing
`_planar_floquet`/`barden_stability` path. **Gate (sourced, no self-reference):** recover ≥ the
four Table-1 families and match max monodromy eigenvalues (1036.116088 / 1.000008 / 4445.387515 /
28178.258323) to a stated tolerance justified by the corrector's own convergence (expect
~1e-3-relative given the paper's 1e-11 shooting floor and OCR-limited digits); negative result =
family-by-family honest report, not a silent pass.
**Task B** (dispatch only after A passes): generalize `heteroclinic_cycle.correct_connection` to
a resonant-member node type + homoclinic A=B mode + one-sided {y=0} section; gate on Tables 2/3
connection states. Then (separately authorized) the new-system discovery campaign.

If instead the user weighs the pivot strictly (`[[project_capability_frontier_complete]]`): the
NO-GO case would rest on "MVP gates only reproduce published results" — but at ~2 tasks of
assembly on existing validated parts, with a new catalogue-eligible object class and five
unstudied moon systems behind it, that case is weaker than it was for the N=5 build the project
already chose to do at larger cost. The pivot rule was meant to stop multi-week speculative
builds, not two-task assemblies with sourced gates.
