# #347 Phase 0 — Floquet bifurcation framework design doc

**Status:** Phase 0 design — no code. Awaiting user review.
**Date:** 2026-06-17.
**Predecessors:** #284 / #343 (asymmetric scan Phase 1+2 clean negatives, commit `1f93128`).
**Successors (proposed):** Phase 1 reproduction substrate (5 days); Phase 2 discovery sweep; Phase 3 µ-scaling; Phase 4 gauntlet adaptation.

---

## Source-availability finding (read this first)

The #347 brief named **Braik & Ross 2026, "Orbital Networks in the Three-Body Problem"** (`papers/braik-ross-2026-orbital-networks-three-body-problem-arxiv-2605.31543.pdf`) as the method reference. After reading the paper end-to-end against the brief's expected content:

- **Braik-Ross 2026 is NOT a Floquet-bifurcation-continuation paper.** It is a *reachable-set / network accessibility* analysis. The 13 representative periodic-orbit families it uses (Table 2, p.12) — including (1,1)a-cycler, (1,1)b-cycler, (2,1)-cycler, (3,2)-cycler — are taken as INPUT from external catalogs ([Jet Propulsion Laboratory 2025], [Rawat et al. 2026], [Ross and Roberts-Tsoukkas 2025]), not computed via bifurcation continuation in the paper itself.
- "Floquet" appears in the paper exactly twice — in equation (20) on p.12, defining a per-family *instability rate* σ = ln(|λmax|)/T from the largest non-trivial monodromy multiplier. This is a **stability diagnostic**, not a discovery method.
- The paper cites Doedel et al. 2007 and Hénon 1997 for "geometry, stability, and bifurcation structure" of cislunar families (p.4) and Guzzetti et al. 2016 / Leiva-Briozzo 2006 / Zimovan 2017 for "systematic coverage of the Earth-Moon family landscape." **None of those papers are present in `/home/bruce/dev/cyclers_pdf/papers/`** (corpus survey, 79 PDFs, 2026-06-17).

The fallback the brief offered — `cuevas-del-valle-2023-optimal-floquet-stationkeeping-relative-dynamics-three-body-aerospace.pdf` — was checked too. It uses Floquet theory for **stationkeeping control** (a hybrid SDRE-Koopman scheme on top of the linearized Floquet basis), not for **family discovery**. Wrong tool for the same reason.

**However**, the *correct* reference IS in our corpus, and the brief's named #347 outcome is still reachable through it: **Roberts-Tsoukkas & Ross 2026, "Stable Prograde Earth-Moon Cyclers"** (`papers/roberts-tsoukkas-ross-2026-stable-prograde-em-cyclers-journal.pdf`) is the upstream paper that produced the (1,1)a/(1,1)b/(3,2) representatives Braik-Ross consumed. It documents:

- The **saddle-center bifurcation** as the canonical mechanism by which symmetric cycler families are created (Fig. 4, p.5, for the (3,2) Earth-Moon family).
- Two **single-parameter continuation** schemes (vary x0 hold C, vary C hold x0; p.4).
- An explicit open problem (verbatim, p.6): *"it is known that asymmetric cyclers exist, but due to their numerical complexity, much less is known about these families and it is not so clear that the same results will hold."*

This is exactly the open frontier #284/#343 ran into. The rest of this doc is grounded in Roberts-Tsoukkas-Ross 2026 (cited as **RTR2026** below), plus the project's own pre-built Floquet substrate.

---

## Section 1 — Background & problem statement

The #343 / #284 asymmetric scan ran 1,944 (k1, k2) seed cells across 9 bands and converged 12 orbits — **all** of which had topology (2, 0) or (3, 0), i.e. symmetric. The corrector's basin of attraction is symmetric-only; no amount of seed-grid expansion will reach the asymmetric branch. This is structural, not coverage-limited.

In the project's `winding_topology` nomenclature (`src/cyclerfinder/search/binary_star_search.py:63-69`), (k1, k2) is a pair of integer winding numbers about the two primaries — k1 = signed windings about P1 (Earth in EM), k2 about P2 (Moon). The cyclers RTR2026 documented at length are all **symmetric** in the time-reversal sense R(x,y,ẋ,ẏ) = (x,-y,-ẋ,ẏ) (Braik-Ross p.9, eq. 12), which constrains them to perpendicular x-axis crossings — the corrector exploits ẋ=0 at t=T/2 to halve the residual. Symmetric cyclers cover (m, 0), (0, n), and any (k1, k2) where k1=k2 (the diagonal). The Umbriel-Oberon (1,1) admission (#339) is on this diagonal.

The genuinely open class — and the one #284 was trying to reach — is **asymmetric (k1, k2) with k1 ≠ k2**, including (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2). RTR2026 explicitly notes these exist but were not pursued in their study due to numerical complexity. Braik-Ross's Table 2 row "C32: (3,2)-cycler" suggests (3,2) is in fact reachable, but its provenance traces back to JPL/Rawat catalogs — i.e. it was computed somewhere else, by methods we don't have in our corpus.

The structural answer to the #343 wall is **branch-switching at a bifurcation**: when a one-parameter family of symmetric orbits passes a parameter value where a non-trivial Floquet multiplier touches +1 (saddle-center / pitchfork) or a primitive k-th root of unity (period-multiplying), a new family branches off in the direction of the marginal eigenvector. For the symmetry-breaking (pitchfork) case at λ=+1, the branching family is asymmetric. This is the classical Hamiltonian-with-symmetry bifurcation theory of [Golubitsky-Stewart-Marsden] (cited by RTR2026, ref. [26]).

---

## Section 2 — Method digest (RTR2026 + Hamiltonian bifurcation theory)

Because the brief's named paper (Braik-Ross 2026) does not contain the method, this section pulls from **RTR2026** (the correct upstream) plus the **project's own pre-existing Floquet substrate** (`tulip.py`, `bifurcation_detector.py`, `family_switch.py`).

**The saddle-center bifurcation (RTR2026, "Continuation Methods" + "Theoretical Results" p.5).** Verbatim from RTR2026 p.5: *"In all symmetric cycler-families observed, there exists a maximum value of C for which the family exists. At this maximum value of C, a saddle-center bifurcation occurs ... The saddle-center bifurcation, always results in two branches of periodic orbits near the bifurcation point — a stable branch and an unstable branch."* Figure 4 (RTR2026 p.5) illustrates this for the (3,2) family on a (C, x0) plot, showing the two branches meeting at a critical C*.

The saddle-center is the bifurcation at λ=+1 with **codimension 1**: the trivial Floquet multiplier pair (1, 1) is augmented by a second eigenvalue arriving at +1 along the family. The classical pitchfork — the symmetry-breaking sibling — also has a non-trivial multiplier crossing +1, but the branching family BREAKS the time-reversal symmetry. Both are detected by the same condition: a non-trivial Floquet multiplier on the unit circle at +1.

**Continuation schemes (RTR2026 p.4).** Verbatim: *"Once a (k1, k2) cycler is found, nearby cyclers in the same family can be found trivially through continuation methods. Two single-parameter based methods were employed in the present study. In the first method, C0 is incremented by a small δC and held constant, and x0 is varied in order to locate the new periodic orbit. Conversely, the second method involves incrementing x0 by a small δx0 and varying C0 until a new periodic orbit is found. If the orbit is assumed to be symmetric, the condition that ẋ = 0 for a periodic orbit can be used as a targeting condition."*

This is natural-parameter continuation, NOT pseudo-arclength. It fails at folds (where the family turns back in the continuation parameter), which is exactly where the saddle-center sits. **Pseudo-arclength** is the standard upgrade — instead of fixing one of (C, x0), we follow arclength s along the family and treat both (C, x0) plus the period T as free, closing with an arclength constraint. The project's `cr3bp_general_periodic_3d.py` (#291) header already uses pseudo-arclength in z0 for the 3D family tracer; the pattern is in the repo.

**Asymmetric-branch corrector.** The symmetric corrector enforces ẋ=0 at T/2 — three residuals on three unknowns (z0, ẏ0, T) for the (z, y, ẋ, ż)|T/2 = 0 closure. For an asymmetric orbit, this is wrong: the orbit does NOT cross the x-axis perpendicularly at T/2. The correct residual is **full-period closure**: state(T) - state(0) = 0, six residuals on seven unknowns (six IC components + T), closed minimum-norm. The project's `cr3bp_general_periodic_3d.py` already has this mode (Phase 1 mode 2: "Full 3D asymmetric (broken-plane) — free vars (x0, y0, z0, ẋ0, ẏ0, ż0, T), residual (x-x0, ..., ż-ż0) at T"). For the planar case this collapses to four residuals on five unknowns; same structure.

**Branching off the symmetric family (Hamiltonian-with-symmetry theory).** At the bifurcation point on the parent (symmetric) family, the monodromy has a non-trivial +1 eigenvalue. The corresponding right-eigenvector ξ ∈ R^6 points along the branching family in state space. The standard recipe (analogous to `genome/family_switch.py:64-125` for period-multiplying):

1. Read parent monodromy M, eigendecompose.
2. Pick the eigenvalue closest to +1 that is NOT one of the trivial pair (energy + time-translation; these are always (1, 1) for an autonomous Hamiltonian system).
3. Extract the matching real eigenvector ξ.
4. Perturb parent IC s₀ → s₀ + ε·ξ for small ε (~1e-3 LU).
5. Hand to the **asymmetric** corrector (full-period closure), NOT the symmetric one. *This is the step `family_switch.py` does wrong for our use case — it always re-corrects via `correct_symmetric_nrho`, which is why the asymmetric branch is invisible to it.*

**Termination / validation.** RTR2026 does not explicitly cover validation. The project's own discipline applies: independent topology cross-check (winding_topology on the converged orbit must yield the target (k1, k2)); independent Jacobi cross-check via Radau-integrated state(T); Floquet stability check via `bifurcation_detector.monodromy + floquet_multipliers`.

---

## Section 3 — Integration with existing project machinery

**Re-usable as-is:**

- `src/cyclerfinder/search/bifurcation_detector.py:117-161` — `monodromy(system, state0, period)`. Already integrates state+STM via the existing DOP853 propagator at rtol=atol=1e-12. **No changes needed.**
- `src/cyclerfinder/search/bifurcation_detector.py:164-175` — `floquet_multipliers(monodromy)`. Eigendecomposition with reciprocal-pair-aware sorting. **No changes needed.**
- `src/cyclerfinder/search/cr3bp_general_periodic_3d.py` — single-shooting Newton corrector with configurable free-vars + residuals on full 6D state. Already supports both symmetric (z0, ẏ0, T) and asymmetric (full 6D + T) modes. **This is the asymmetric corrector Phase 1 needs**; do not reinvent.
- `src/cyclerfinder/search/binary_star_search.py:63-95` — `winding_topology(mu, state0, period)`. Independent topology classifier. **Phase 1 reproduction gate.**
- `src/cyclerfinder/genome/tulip.py` — sourced-IC + reproduce-gate pattern. **Borrow the dataclass structure for asymmetric IC tables.**

**Gaps that need filling:**

- `bifurcation_detector.detect_period_multiplying` excludes k=1 (line 195: *"k = 1 is excluded externally"*). The **pitchfork / saddle-center detector** is a new function: scan adjacent family members for a non-trivial multiplier crossing +1 (sign-flip on `dist(λ, +1) - tol` after excluding the two trivial unit multipliers). This is ~30 LOC alongside `detect_period_multiplying`.
- `family_switch.switch_family` hard-codes the symmetric corrector. The **asymmetric branch corrector** is a new entry point that takes the bifurcation eigenvector + a perturbation amplitude and hands to `cr3bp_general_periodic_3d` in full-asymmetric mode. ~80 LOC.
- **Pseudo-arclength continuation in the (C, x0) plane** is not in the repo for cyclers (`cr3bp_general_periodic_3d` has it for z0). New module or extension. ~150 LOC.
- **Sourced asymmetric IC anchor.** RTR2026 documents existence but publishes no asymmetric ICs. The (1,1)a/(1,1)b/(2,1)/(3,2) ICs cited in Braik-Ross Table 2 trace to upstream catalogs we don't hold. **The Phase 1 anchor must be a symmetric IC we already have**, with the bifurcation hop discovering the asymmetric branch as the reproduction test.

**Modules to create:**

- `src/cyclerfinder/search/bifurcation_detector.py` — add `detect_saddle_center` / `detect_pitchfork` (one function with a `bifurcation_type` enum).
- `src/cyclerfinder/genome/asymmetric_branch.py` — new module mirroring `family_switch.py` but routing through `cr3bp_general_periodic_3d` asymmetric mode.
- `src/cyclerfinder/search/cr3bp_jacobi_arclength.py` — already exists (per `ls` output above); audit whether it covers the (C, x0) pseudo-arclength scheme RTR2026 describes.

---

## Section 4 — Phase 1 work plan

**Goal:** A reproduction substrate that takes a published symmetric cycler IC, walks its family in C via natural-parameter continuation, detects the saddle-center bifurcation, branches off via the marginal eigenvector, corrects via the asymmetric corrector, and confirms the branched orbit has winding topology distinct from the parent (or the same with stability flipped). **NOT a discovery sweep.** Discovery is Phase 2.

**Anchor IC.** The (3,2) Earth-Moon symmetric cycler. Provenance: Braik-Ross 2026 Table 2 row C32 (P=78.613 days, σd=0.1583 day⁻¹, CJ=3.1294) — these are reproducible numerical anchors. RTR2026 Figure 4 illustrates the (3,2) saddle-center directly. **Phase 1 success = either reproduce Braik-Ross's (3,2)-cycler period to ≤ 1% AND find the saddle-center bifurcation along its family OR fail and report which step broke.**

**Task list (each ≤ 1 day):**

- **P1.1 (≤ 0.5d)** — Reproduce-gate the (3,2) Earth-Moon symmetric anchor. New test in `tests/genome/test_asymmetric_anchor.py`. Use the existing `cr3bp_general_periodic_3d` symmetric mode to converge a (3,2) cycler near CJ=3.1294. Cross-check period against Braik-Ross 78.613 days. **Go/no-go:** period within 1%, Jacobi within 1e-9, winding_topology = (3, 2). No-go → STOP, raise scope: need a sourced IC, not a seed.
- **P1.2 (≤ 0.5d)** — Natural-parameter continuation in CJ along the symmetric (3,2) family. Use the existing `search/continuation.py` or `cr3bp_continuation.py` (audit which is appropriate). Walk from the anchor toward the critical CJ* (RTR2026 indicates CJ* > 3.1294 for (3,2)). Log monodromy + Floquet multipliers at each member. ~10 family members. **Go/no-go:** continuation runs to a fold (residual diverges or step adapts to zero) OR explicitly reports "no fold within range" with a sensible upper bound. No-go on step-adaptation failure → broaden step control; this is mechanical, not a method failure.
- **P1.3 (≤ 1d)** — Add `detect_saddle_center` to `bifurcation_detector.py`. Scan the P1.2 family for the first non-trivial multiplier reaching +1 (within tol=1e-3). Returns a `BifurcationPoint` with `k=1` (re-use the existing dataclass with a `bifurcation_type` discriminator added). **Go/no-go:** detector flags exactly one saddle-center on the (3,2) family within the P1.2 range. No-go → check whether the trivial pair was correctly excluded (the two energy/time-translation eigenvalues are always at +1; the detector must skip them).
- **P1.4 (≤ 1d)** — New module `genome/asymmetric_branch.py`. Mirrors `family_switch.py:switch_family` but: (a) picks the non-trivial +1 eigenvector at the bifurcation; (b) perturbs IC by ε·ξ with default ε=1e-3 (test 1e-4, 1e-3, 1e-2); (c) hands to `cr3bp_general_periodic_3d` in **full-asymmetric mode** (free vars: x0, y0, z0=0, ẋ0, ẏ0, ż0=0, T; residual: full-state closure at T). **Go/no-go:** corrector converges to residual < 1e-10 from the perturbed seed AND the resulting orbit has winding_topology ≠ parent's. No-go on convergence → multi-shoot the asymmetric corrector (single-shooting may be unstable along the marginal direction); multi-shooting infrastructure already exists in `multi_shooting.py`.
- **P1.5 (≤ 1d)** — Reproduction-of-record test. New test `tests/genome/test_asymmetric_reproduction.py`. End-to-end: anchor → continue → detect → branch → corrector → topology check. Output a JSONL row: `(parent_state0, parent_T, parent_k1, parent_k2, branch_state0, branch_T, branch_k1, branch_k2, residual, max_floquet_mag)`. **Go/no-go (Phase 1 exit):** one JSONL row written with a branched orbit whose (k1, k2) is either (a) the same but different stability character, OR (b) genuinely (k1', k2') ≠ (3, 2). Either is a successful demonstration that the asymmetric corrector escapes the symmetric basin.

**Total Phase 1 effort:** 4 days nominal, 5 days with buffer.

**Phase 1 deliverable:** one passing test, one JSONL row, one design-doc update with the actual numerics found. **No catalogue writeback.**

---

## Section 5 — Risks + open questions

- **Is the saddle-center C* sourceable from RTR2026?** Their Figure 4 shows the (3,2) bifurcation qualitatively but does NOT publish the numerical C*. Likely we have to search for it via P1.2 continuation. **Confidence in finding it via continuation: high** — RTR2026 Fig. 4 is on the same family Braik-Ross Table 2 uses; the bifurcation is at the top of the family's C range, finite-distance from CJ=3.1294.
- **Does the project's DOP853 integrator hold corrector tolerance at the bifurcation?** Near a saddle-center, the monodromy has multipliers near (+1, +1, +1, +1, …) — clustered eigenvalues are ill-conditioned. The existing `monodromy` runs at rtol=atol=1e-12 (the project's standard). **Unknown** whether the eigenvalue separation will survive at the cluster point. Mitigation: cross-check the multiplier via the trace/determinant invariants (the characteristic polynomial coefficients should match a known structure for a Hamiltonian monodromy). Add this as a check in P1.3.
- **Will the asymmetric branch survive physical-sanity / V0-V1 gates?** The asymmetric branch may have ΔV-at-perilune outside cycler regime, or violate the lit-fresh corpus check. **This is the gauntlet question, deferred to Phase 4.** Phase 1 only demands corrector convergence + topology change.
- **Could the bifurcation actually be a Neimark-Sacker (multiplier crossing the unit circle off-axis)?** RTR2026 only claims saddle-center; the bifurcation_detector handles k=2 already (period-doubling at λ=-1) but Neimark-Sacker (complex pair leaving the unit circle) is a separate mechanism — it spawns invariant tori, not new periodic families. **Mitigation:** P1.3 explicitly checks the multiplier is real (+1), not complex.
- **Effort confidence: medium.** Phase 1 leans heavily on infrastructure that already exists, but the asymmetric corrector has never been exercised on a bifurcation branch in this codebase. The #287 spike (cited in `cr3bp_general_periodic_3d.py:30`) closed a 3D (1,1) Braik-Ross member to 1e-13 residual via the asymmetric mode, which is encouraging. The unknown is the eigenvector-step seeding — if ε is too small the corrector lands back on the parent, if too large it diverges. P1.4 budgets one day for ε-sweeping; could blow out to 1.5 days.

---

## Section 6 — Phase 2-4 sketch

- **Phase 2 — Discovery sweep on Earth-Moon.** Run P1's substrate across the symmetric cycler families Braik-Ross Table 2 enumerates: (1,1)a, (1,1)b, (2,1), (3,2) on Earth-Moon. For each, walk its C-range, detect saddle-center, branch via asymmetric corrector, log the branched (k1', k2'). Replaces the #284/#343 seed-grid sweep with a Floquet-anchored sweep. Expected outcome: a handful of genuinely (k1, k2) ≠ (k2, k1) cyclers entering the V0/V1 pipeline. Effort: ~2 weeks. Anchor source: RTR2026.

- **Phase 3 — µ-scaling.** Re-run Phase 2 substrate with the existing `mu_continuation.py` to Saturn (Enceladus-Mimas, Titan-Hyperion), Uranus (already touched by #285/#341), and binary-star µ values RTR2026 reaches (µ=0.1, 0.3). Each system needs a sourced symmetric anchor as the seed; absent one we run the Phase 2 substrate from a known symmetric cycler discovered in that system. Effort: ~1 week per system, mostly substrate-stable.

- **Phase 4 — V0-V5 gauntlet adaptation for the asymmetric class.** The existing gauntlet (#319 for QP-tori, #340 for symmetric cyclers) needs an asymmetric-aware physical-sanity gate (the V∞-vs-escape check pre-uses symmetric closure), an asymmetric-aware lit-fresh check (the corpus-anchor matcher currently keys on (k1, k2) and symmetric topology), and an asymmetric ML-flagger calibration. Effort: ~1 week.

---

## Decision gate for the user

**Approve Phase 1 (5-day reproduction substrate)?** Yes/no/redirect.

- If **yes**, Phase 1 fires per Section 4.
- If **redirect**: most likely alternative is to acquire one of Doedel 2007, Hénon 1997, Guzzetti 2016, Leiva-Briozzo 2006, or Zimovan 2017 (Braik-Ross's cited sources for the bifurcation method) before Phase 1, which would shift Section 2's method digest from RTR2026 + project-internal substrate to a directly sourced bifurcation-continuation recipe. Effort impact: +0.5 to +1 week (acquisition + digest) but stronger sourcing for the Phase 2 sweep.
- If **decline**: the #343 wall remains the rate-limiting step on genuine asymmetric (k1, k2) cyclers; the catalogue's open class at V0 is not closeable without it.
