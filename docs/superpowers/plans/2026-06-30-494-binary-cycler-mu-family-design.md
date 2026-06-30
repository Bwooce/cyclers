# #494 — Binary / circumbinary (k₁,k₂)-cycler μ-family — design

**Date:** 2026-06-30. **Status:** design (build not started). **Supersedes the framing of**
#315 (binary-star μ-gap sweep), and folds in the #252/#255 binary-star negatives.
**Task:** #494 (allocated; next-unused per [[project_task_numbering_convention]]).

## TL;DR — the honest reframe

This was queued as "build a net-new circumbinary (P-type) genome to find novel Pluto-Charon
cyclers." The grounding (B-sourcing agent 2026-06-30 + repo survey) shows that framing is **wrong
on two counts**, both in our favour:

1. **It is NOT net-new dynamics.** The CR3BP solver is already μ-agnostic
   (`core/cr3bp.py::cr3bp_eom(t, state, mu)`; `search/binary_star_search.py::_system(mu)`
   hand-builds a system at *any* μ; `search/mu_continuation.py` already continues (k₁,k₂)
   cyclers in μ). The binary (k₁,k₂)-cycler is a periodic orbit alternately captured about each
   primary — the **same object** the project's #252/#255/#315 thread was chasing, and the
   **same object** Ross & Roberts-Tsoukkas study.

2. **It is primarily REPRODUCTION + census extension, not novel discovery.** The project has
   **already reproduced** the Earth-Moon slice of this family: the 5 stable prograde EM
   (1,1)/(3,1) cyclers (`ross-rt-em-cycler-*`, catalogue lines ~48645+) are **catalogued at V2**
   (#229), recovered from the published `(μ, C, T)` via a 1-D symmetric-orbit solve
   (`docs/notes/2026-06-11-ross-roberts-tsoukkas-2025-mining.md`). The NEW published source —
   **Ross & Roberts-Tsoukkas 2026, arXiv:2606.29189, "Stable Families of Ballistic Prograde
   Cyclers in the RTBP"** (the journal extension of AAS 25-621) — extends that family across
   **μ = 0.001 (Sun-Jupiter) → 0.01215 (EM, done) → 0.5 (equal-mass binary)**, i.e. exactly the
   #315 μ-gap.

So #494 = **mine the now-published μ-family goldens and reproduce the binary (k₁,k₂)-cycler
across the μ range using the EXISTING machinery, then instantiate at Pluto-Charon (μ=0.1085).**
This closes the three-task-old #315/#252/#255 thread as a *positive* result and gives the
catalogue its binary-cycler class — modest novelty (a fresh real-system instantiation, the
#312-Uranus framing), not a novel-discovery claim.

## Two regimes — keep them distinct

| Regime | What it is | Feasible cycler? | Golden |
|---|---|---|---|
| **(k₁,k₂) capture cycler** (the target) | Periodic orbit alternately captured about Pluto AND Charon; Charon (GM ≈ 106 km³/s²) supplies the energy change. The Persephone-tour regime. μ=0.1085. | **Yes** — Charon is massive enough (this is why #489/#492's small-moon objection does NOT apply here). | Ross-RT 2026 `(μ,C,T)` representatives |
| **P-type exterior circumbinary** (among small moons) | Test particle exterior to the binary, orbiting both, flying by Styx/Nix/Kerberos/Hydra. | **No** — small moons too tiny for assists (#489/#492, confirmed). | Holman-Wiegert a_crit (stability bound only) |

The catalogue target is the **(k₁,k₂) capture cycler**. The P-type exterior regime is used only
as a *validation* anchor (the Holman-Wiegert critical-stability radius is a clean sourced golden
for the exterior dynamics of the same μ-agnostic solver) — not as a cycler source.

## Goldens (sourced only — never circular; [[feedback_golden_tests_sourced_only]])

1. **PRIMARY — Ross & Roberts-Tsoukkas 2026** (arXiv:2606.29189; journal ext. of AAS 25-621).
   Stable ballistic prograde (k₁,k₂)-cyclers, families traced in the (x₀, C) plane by
   pseudo-arclength continuation, born at a saddle-center bifurcation at maximal C; symmetric
   members satisfy ẋ=0 at the perpendicular +x̂ crossing → recoverable from `(μ, C, T)` by the
   existing 1-D solve. **ACQUISITION REQUIRED** (paper is 2 days old; NOT yet in corpus — the
   filed Ross items are AAS 25-621 (EM only) + Braik-Ross 2605.31543 (orbital networks) +
   Roberts-Tsoukkas 2026 multi-orbiter, all distinct). Mine its per-μ representative table the
   way `2026-06-11-ross-roberts-tsoukkas-2025-mining.md` mined the EM slice. **Honest risk:** the
   B-sourcing agent read the arXiv HTML and saw (x₀,C) in *figures*; the 2025 AAS version printed
   15-16-digit `(μ,C,T)` per stable representative, so the journal version probably tabulates them
   too — confirm on acquisition. If only figures, this degrades to a digitize-the-figure golden
   (still sourced) per [[feedback_never_give_up_reproducing_papers]].

2. **SECONDARY — Holman & Wiegert 1999** (AJ 117, 621; DOI 10.1086/300695). P-type critical
   semimajor axis, e=0 reduction `a_c/a_bin = 1.60 + 4.12μ − 5.09μ²`. **Verified in-repo
   2026-06-30:** μ=0.1085 → 1.9871 → **a_crit = 38,939 km** (falls just inside Styx at 42,656 km
   — Kenyon-Bromley's "Styx near innermost stable orbit" corroborates); μ=0.5 → 2.39 (reproduces
   HW99's own quoted value). Closed-form, ready now; gates the exterior-stability behaviour of the
   solver at the Pluto-Charon μ.

3. **TERTIARY — Jbara 2025** (arXiv:2510.13479, non-refereed; use only as an *independently
   recomputable* cross-check): Pluto-Charon CR3BP `C_L1 = 3.6210` at μ≈0.109. L-point Jacobi
   values are deterministic in μ → recompute ours and cross-check. Optional: Langford & Weiss 2023
   (AJ 165, 140) Zenodo ICs (10.5281/zenodo.7532982) for a circumbinary periodic-orbit corrector
   golden if a P-type periodic orbit is ever wanted.

## Reuse map (do NOT rebuild)

- **Solver/corrector:** `core/cr3bp.py` — `cr3bp_eom`, `cr3bp_stm_eom`, `propagate(..., with_stm,
  stm_mode)`, `jacobi_constant`, `jacobi_gap_dv_min`; `search/cr3bp_periodic.py` —
  `correct_symmetric_fixed_jacobi(...)`, `ydot0_from_jacobi(...)` (Ross Eq. 9),
  `barden_stability(...)` (|ν|<1 ⇒ linearly stable), `crosscheck_periodic(...)`.
- **μ machinery:** `search/mu_continuation.py` — `continue_in_mu(...)`, `scan_c_family_at_mu(...)`,
  `MuMember/MuBranch/MuStopReason`. `search/binary_star_search.py` — `figure_seeded_search(mu, ...)`,
  `winding_topology(...)` (the U1-/U2+ **figure-match gate** — the binding test per #252, NOT
  residual/stability), `topology_3d(...)`, `collinear_lpoints(mu)`, `_system(mu)`.
- **Gauntlet entry (from the #252 note):** a self-computed periodic orbit + independent-integrator
  closure + Barden verdict in a non-physical μ is a true **V0** discovery; **V1** if the member is
  recovered from a *sourced* `(μ,C,T)` (the EM five are V2 via same-model corrector + 100-period
  REBOUND/IAS15 bounded-band + Barden — `2026-06-13-ross-v2-longspan-evidence.md`). Real-ephemeris
  Pluto-Charon validation (SPICE) is the V3/V4 lever, optional.
- **Data:** `core/satellites.py` Pluto block — **Charon (GM 106.1, a 19,600) ✓; Nix ✓; Hydra ✓;
  Styx + Kerberos MISSING.** Add Styx/Kerberos from Gakis-Gourgouliatos 2022 (Styx a=42,656 km,
  P=20.16 d; Kerberos a=57,783 km, P=32.17 d) — needed only for the exterior-regime sanity, not
  the (k₁,k₂) cycler. System GM 975.5 in-repo matches Brozović (verify Charon 106.1 vs 105.9).

## Phased build (each phase positive-controlled)

- **Phase 0 — positive control (no new data needed; start here).** Re-reproduce the 5 catalogued
  EM V2 members from `(μ=0.01215, C, T)` through `scan_c_family_at_mu` + `correct_symmetric_fixed_jacobi`,
  and confirm `winding_topology` returns the published (1,1)/(3,1). This proves the lane still works
  before trusting any μ-extension all-negative ([[feedback_verify_gauntlet_with_positive_control]]).
- **Phase 1 — acquire + mine arXiv:2606.29189** → per-μ `(μ,C,T)` representative table (or digitized
  figure). File + digest + CORPUS_INDEX per [[feedback_corpus_document_policy]].
- **Phase 2 — μ-family reproduction.** Recover the sourced representatives at μ ∈ {0.001 SJ, 0.01215
  EM✓, …, 0.5 equal-mass}; winding-topology figure-match each; Barden stability. Sourced ⇒ **V1**
  rows. **Closes #315/#252/#255** (the μ-continuation failed because it branch-switched at a fold;
  the direct fixed-μ recover-from-(μ,C,T) is the lever #252 named).
- **Phase 3 — Pluto-Charon instantiation (μ=0.1085).** Direct fixed-μ search for the (k₁,k₂) member
  at Pluto-Charon; cross-check Holman-Wiegert a_crit (exterior) + Jbara C_L1 (Jacobi). If Ross-RT
  2026's continuous family covers μ≈0.109, the member is **V0/V1-known** (reproduction); the
  Pluto-Charon-specific real-system entry is the modest-novelty deliverable ("first Pluto-Charon
  cycler" = fresh instantiation of a published family — the #312-Uranus framing). **No novel-
  discovery claim** beyond what literature_check clears.
- **Phase 4 — gauntlet + catalogue.** V0→V1 per the entry above; flag any candidate for human
  adjudication; run ALL data+search ratchets before any catalogue edit
  ([[feedback_catalogue_edits_run_all_ratchets]]).

## Honest novelty standing

Primary value = **closing a long-open thread (#315/#252/#255) as a positive reproduction + giving
the catalogue its binary (k₁,k₂)-cycler μ-family class**, anchored to a brand-new published source.
Pluto-Charon adds a fresh real-system instantiation (modest novelty). This is consistent with the
#492 lesson: genuinely-novel cyclers are scarce; faithful sourced census extension is real
deliverable value. Acquisition of arXiv:2606.29189 is the one external dependency; Phase 0 can
start immediately without it.

## Corpus acquisition targets (user-fetched PDFs, I file/digest/index)

- **Ross & Roberts-Tsoukkas 2026** — arXiv:2606.29189 (PRIMARY golden; blocking Phase 1+).
- Holman & Wiegert 1999 (AJ 117, 621) — secondary golden (closed-form already verified; PDF for the record).
- Gakis & Gourgouliatos 2022 (arXiv:2202.13319) — Pluto-moon elements (Styx/Kerberos registry).
- Kenyon & Bromley 2019 (arXiv:1810.01277) — Pluto-Charon stability map corroboration.
