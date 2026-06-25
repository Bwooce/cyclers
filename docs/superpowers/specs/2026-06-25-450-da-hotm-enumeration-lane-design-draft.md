# Global multi-rev DA/HOTM fixed-point enumeration lane — DESIGN DRAFT (#450)

**Date:** 2026-06-25
**Status:** DESIGN-DRAFT — for user review. No production code written by this doc.
**Issue:** #450. Ranked the **#1 capability lever** by the discovery-strategy
prioritization draft (`2026-06-25-discovery-strategy-prioritization-design-draft.md`,
§6 lever 2) because it strictly *subsumes* single-family Jacobi continuation and
re-opens the largest method-limited dead region we have a published validation
target for.

This draft answers: **what is the smallest, well-bounded capability that lets us
GLOBALLY enumerate multi-revolution CR3BP periodic-orbit fixed points over a
target energy slice — reaching the families that seed-local continuation
structurally cannot — and how does it plug into the existing genome/gauntlet
pipeline without duplicating it?**

---

## 0. TL;DR

- **Goal:** a global, domain-covering, multi-revolution fixed-point *enumerator*
  for the planar (then spatial) CR3BP, operating on a Poincaré-section return map
  reduced to `(dx0, dxdot0)` by the Jacobi + section constraints (the
  "Png'-class" method, arXiv:2509.12671). It produces candidate fixed-point ICs
  at a target Jacobi `C`; the **existing** `correct_general_periodic` corrector
  refines and certifies each one. The enumerator is a *seed source*, not a new
  corrector — it plugs into the seam our continuation lane already uses.
- **Re-opens:** `cr3bp-em-cj3.00-dro-lyapunov-band-newfamily-2026-06-13` (the EM
  C≈3.0 DRO/Lyapunov band, explicitly recorded EMPTY *conditional on
  single-family continuation*, with the verbatim re-open condition "RE-OPEN this
  region when a DA/HOTM-class GLOBAL multi-rev discovery lane ships"), the three
  `cr3bp-continuation-saturn-*-lyapunov-2026-06-12` Saturn-moon bands, and the
  planar-CR3BP-envelope dead axis (Dead Axis D2) — all flagged in the negative
  registry with the SAME re-open key: a global multi-rev enumeration method.
- **Validation/golden target (the lane-works proof):** recover the published
  **Png' hybrid DRO/Lyapunov-bridge family** of arXiv:2509.12671 — specifically
  re-derive the lowest member **P5g'** (and the P7g'/P9g' members) *as outputs of
  our enumerator*, not handed in as ICs, then close them with
  `correct_general_periodic` to the residual we already achieved (3.45e-12,
  T to 11 sig figs). This is the decisive test that the global sweep surfaces a
  family seed-local continuation provably misses.
- **CRITICAL HONESTY — buildability.** The METHOD is fully digested
  (`2026-06-13-high-order-transfer-map-2509.12671-mining.md`), the target Png'
  ICs are transcribed AND already re-closed under our corrector, and a related DA
  reachable-set spike already exists (`reachable_impulsive.py`). BUT the paper's
  performance contribution — the **differential-algebra (DA) Taylor-map layer**
  (DACEyPy/DACE + Automatic Domain Splitting + a convex polynomial optimizer,
  MOSEK) — is **NOT installable in this environment today** (verified: no
  `daceypy`, `dace`, or `pyaudi` module; no DA/MOSEK in `pyproject.toml`). So the
  honest confidence split is:
  - **The GEOMETRY and the validation are buildable from digested material now**
    via a brute-force/sampling realization (exactly the path the Zhou spike took)
    — this recovers Png' and proves the lane, at the cost of the paper's speed.
  - **The full DA/ADS/SCOP acceleration is ACQUISITION-GATED.** The first plan
    task is therefore a backend decision/acquisition (DACEyPy vendoring + MOSEK
    licence, OR a pure-Python truncated-Taylor-map fallback), NOT method
    invention. We do not need to invent DA math — the paper specifies it — but we
    cannot `import` it yet.

---

## 1. Why the current lanes can't reach this (the method-capability gap)

The negative-registry entry `cr3bp-em-cj3.00-dro-lyapunov-band-newfamily-2026-06-13`
states the gap precisely (verbatim interpretation): our #218/#219 sweep surfaced
**0 distinct new families** in the EM C≈3.0 band, but an independent group found a
previously-undocumented hybrid family (Png') in *exactly* this region. The reason
our lane misses it is structural, not a bug:

| Property of Png' | Why our lane can't reach it |
|---|---|
| **Non-symmetric** (`xdot0 ≠ 0`) | `correct_symmetric_fixed_jacobi` only finds perpendicular-crossing orbits; the symmetric-seed continuation campaign never enumerates the asymmetric class. |
| **Multi-revolution g'-bifurcation** (P5g' is a 5-rev bridge, T=11.175 TU) | Natural-parameter continuation steps *along one family* from a sourced seed; it never *enumerates fixed points over a 2-D Poincaré domain*, so a family with no sourced seed and no continuation path from one is invisible. |
| **DRO/Lyapunov hybrid bridge** | It lives between two known families; continuation from either parent stays on that parent. |

The capability we lack is a **global, domain-covering, multi-revolution
fixed-point enumerator** over a Jacobi slice. This is precisely what the
DA/HOTM method provides and what re-opens the region under the
capability-subsumption rule (a global enumerator strictly subsumes seed-local
continuation).

**Regions re-opened (all carry the same registry re-open key):**

1. `cr3bp-em-cj3.00-dro-lyapunov-band-newfamily-2026-06-13` — EM C≈3.0
   DRO/Lyapunov band (the headline; has the Png' validation target *in* it).
2. `cr3bp-continuation-saturn-mimas-lyapunov-2026-06-12`,
   `…-enceladus-…`, `…-tethys-…` — three Saturn-moon Lyapunov bands, each
   EMPTY only "under this method+version" (`cr3bp-jacobi-continuation-v1`).
3. **Dead Axis D2** (planar CR3BP envelope) — the strategy draft names DA/HOTM
   Png'-class enumeration as the *explicit* re-open lever (§4 D2, §1.3 table).

---

## 2. What a "fixed point" is here, and how multi-rev is parameterized

**Fixed point.** On the Poincaré section `Σ = {y=0, ẏ≥0, 0<x<1−µ}`, define the
single-revolution return map `P : Σ → Σ` (first return to `Σ`). A periodic orbit
of period-n is a **fixed point of the n-th iterate**: `Pⁿ(s) = s`. Because the
Jacobi integral `C(s)=const` removes one DOF and the section `y=0` removes
another, a planar section state is fully described by `(x, ẋ)` — `ẏ` is recovered
from `C` (sign `ẏ≥0`), `y=0`. So the planar map is a function of just
`(dx0, dxdot0)` at fixed `C`; the spatial map adds `(dz0, dzdot0)` (4-D).
This is the same reduction `core/cr3bp.jacobi_constant` +
`cr3bp_periodic.ydot0_from_jacobi` already implement.

**Multi-rev parameterization (the scalability source).** The paper builds ONE
single-revolution map and **composes** it: `Pⁿ = P ∘ P ∘ … ∘ P` (n times). In the
DA realization the single-rev map is a Taylor polynomial in `(dx0,dxdot0)` about a
section reference; composing polynomials is cheap relative to re-propagating, so
arbitrary `n` is reached without new integrations. The revolution count `n` is the
**outer loop parameter** of the enumeration: for each `n ∈ {1,…,N}` and each
target `C` in a band, enumerate the fixed points of `Pⁿ`.

**Why "global", not seed-local.** The enumeration covers the *whole* admissible
`(x, ẋ)` sub-rectangle of `Σ` at the given `C` (via Automatic Domain Splitting in
the DA realization, or a brute-force grid in the sampling realization), so it
finds every fixed point in the domain regardless of whether a sourced seed or a
continuation path to it exists. That is the subsumption.

---

## 3. How candidates are filtered BEFORE the expensive corrector

The enumerator is deliberately a *coarse, cheap* fixed-point locator; the
expensive, high-precision certification stays in the existing corrector. The
filter cascade (cheapest first — the C21 lesson, strategy §5):

1. **Section-residual gate.** A grid/sub-domain cell is a fixed-point candidate
   only if `‖Pⁿ(s) − s‖` on `Σ` drops below a coarse tolerance (e.g. 1e-4). In the
   DA realization this is the polynomial root-isolation (ADS + the
   reduced-Poincaré-orbit convex optimization, RPO/SCOP). In the sampling
   realization it is a contour/zero-crossing detector on the brute-force map
   residual.
2. **Distinct-candidate dedup.** Cluster surviving cells (a family produces a
   continuum of near-residual cells); keep one representative IC per cluster.
3. **Reproduction/known-family screen** — reuse `spatial_novelty_prefilter` and
   the catalogue/JPL-oracle signature dedup so a recovered DRO or Lyapunov member
   is routed to the *reproduction ledger*, not the discovery path. (P5g' must
   survive this — it is genuinely outside JPL's single-period parametrization, as
   the mining note's JPL triangulation already confirmed.)
4. **Then, and only then,** hand the deduped, novel-screened ICs to
   `correct_general_periodic` for the 1e-12 certification + STM/Floquet stability.

This keeps the per-candidate cost of the expensive corrector bounded to the
handful of distinct families, not the thousands of grid cells.

---

## 4. Architecture, new files, and the plug-in seam

The lane is one new genome module + one new search driver + thin glue. It reuses
the propagator, the corrector, the prefilter, and the lit-check unchanged.

```
                 ┌─────────────────────────────────────────────┐
 target (C, n)──▶│  genome/da_hotm_enumerator.py  (NEW)         │
                 │   - build single-rev section map (DA or      │
                 │     sampling backend)                        │
                 │   - compose to Pⁿ                            │
                 │   - locate fixed points over the (x,ẋ) domain│
                 │   - emit coarse candidate ICs at fixed C     │
                 └───────────────┬─────────────────────────────┘
                                 │ candidate (x0, xdot0, C, n)
                                 ▼
                 ┌─────────────────────────────────────────────┐
                 │  search/da_hotm_enumeration.py (NEW driver)  │
                 │   filter cascade §3:                         │
                 │   section-residual → dedup → prefilter/JPL   │
                 └───────────────┬─────────────────────────────┘
                                 │ deduped novel-screened ICs
                                 ▼
   EXISTING  search/cr3bp_general_periodic.py :: correct_general_periodic
             (1e-12 closure + STM Jacobian)  ── unchanged ──
                                 │ GeneralPeriodicOrbit + monodromy
                                 ▼
   EXISTING  novelty pipeline: spatial_novelty_prefilter →
             literature_check.check_literature (V0) → V1/V2 gauntlet
                                 │
                                 ▼
                 reproduction ledger (known) | discovery candidate (novel)
```

**New files:**

| File | Responsibility | Reuses |
|---|---|---|
| `src/cyclerfinder/genome/da_hotm_backend.py` | Abstract `SectionMap` interface + two backends: `SamplingSectionMap` (brute-force, buildable now) and `DASectionMap` (DACEyPy/Taylor, acquisition-gated). One swappable seam. | `core/cr3bp.propagate`, `jacobi_constant`, `ydot0_from_jacobi` |
| `src/cyclerfinder/genome/da_hotm_enumerator.py` | Compose `Pⁿ`, enumerate fixed points over the `(x,ẋ)` domain at target `C`, emit coarse candidate ICs. Backend-agnostic. | `da_hotm_backend` |
| `src/cyclerfinder/search/da_hotm_enumeration.py` | Driver: loop `(C, n)`, run filter cascade §3, dedup, route to corrector, classify novel vs reproduction. | `correct_general_periodic`, `spatial_novelty_prefilter`, `jpl_periodic_orbits`, `literature_check` |
| `data/golden/png_hybrid_family.yaml` | Sourced golden: the five arXiv:2509.12671 Png' ICs (P5g'/P7g'-I/P7g'-II/P7g'/P9g' at C=3.00022; P3g' at C=3.020052) + published periods. Golden EXPECTED side traces only to the paper (per `feedback_golden_tests_sourced_only`). | — |

**The plug-in seam.** `correct_general_periodic(system, x0, xdot0, c_target, …)`
is *already* the asymmetric, fixed-Jacobi, y=0-return-map corrector the Png' class
needs (it is non-symmetric and STM-Jacobian'd; the mining note confirms P5g'
closes through `cr3bp_periodic.correct_periodic` and the general corrector handles
the asymmetric class). The enumerator therefore produces exactly the
`(x0, xdot0, c_target)` tuple this corrector consumes — **no corrector change,
no continuation-driver change**. The lane is additive, mirroring how
`asymmetric_novel_scan_parallel.py` and `multi_axis_search.py` already feed this
corrector.

**Catalogue-scope boundary (carried from the mining note, load-bearing).** Png'
is a CR3BP **periodic orbit, NOT a cycler** (no Earth–Earth resonant transfer
leg; it is a DRO/Lyapunov bridge in the rotating frame). It is a **cross-check /
validation target only — never a catalogue row**. The lane's *catalogue* payoff
is only realized if a future target band yields a fixed point that *also* carries
a cycler transfer structure; the lane itself is a capability that makes such a
find *possible*, with Png' as the proof it works.

---

## 5. Data flow

1. **Input:** `(system, C_band, n_range, domain_box)` — e.g. EM, `C∈[3.000,3.021]`,
   `n∈{1,3,5,7,9}`, `(x,ẋ)` box covering the DRO/Lyapunov sub-rectangle.
2. **Backend builds the single-rev section map** about a reference on `Σ`.
3. **Enumerator composes `Pⁿ`** and locates fixed points → coarse candidate ICs.
4. **Driver filter cascade** (§3): section-residual → dedup → prefilter/JPL/
   catalogue dedup.
5. **`correct_general_periodic`** certifies survivors to 1e-12 + monodromy/Floquet.
6. **Novelty pipeline:** `spatial_novelty_prefilter` → `check_literature` (V0,
   necessary-not-sufficient) → V1/V2 if a cycler structure is present.
7. **Output:** classified ledger — `{reproduction | known-family | novel-PO |
   novel-cycler-candidate}` with provenance, residual, period, stability.
   Empty bands re-stamp the corresponding `empty_regions.jsonl` entry with the
   new method+version (capability-subsumption record), NOT silent deletion.

---

## 6. Validation strategy (the sourced/golden proof the lane works)

**Primary golden (decisive).** Recover the **Png' family as an OUTPUT of the
enumerator**, then close with `correct_general_periodic`:

- **Lane-recovery test (the new capability proof):** run the enumerator on EM,
  `C=3.00022`, `n=5`, over the published domain box; assert it *emits a coarse
  candidate within tolerance of P5g'* `(x0≈0.8074, ẋ0≈−0.0956, T≈11.175)`. This
  is the test seed-local continuation cannot pass — it certifies the *global
  enumeration* capability, not just the corrector.
- **Closure golden (re-uses an already-proven result):** the emitted candidate,
  fed to `correct_general_periodic`, must re-close to residual ≤1e-11 and
  `T=11.17510869…` to ~11 sig figs (we already achieved 3.45e-12 / 11 sig figs in
  the mining note — this is a regression lock, EXPECTED side sourced to the paper +
  the prior reproduction).
- **Family breadth:** repeat for P7g'-I, P7g'-II, P7g', P9g' (n=7,9, same C) and
  P3g' (n=3, C=3.020052) from `data/golden/png_hybrid_family.yaml`. Each is a
  sourced cross-check target the mining note explicitly stages.
- **Base-family triangulation (independence):** the enumerator must ALSO recover
  the n=1 DRO and L1/L2 Lyapunov members at C≈3.0002, which JPL independently
  confirms (mining note's JPL triangulation) — proving the enumerator finds the
  known families *and* the hybrid, i.e. it is genuinely global, not Png'-tuned.

**Why this is golden-clean.** The EXPECTED side traces to arXiv:2509.12671's
printed tables + JPL's live oracle — never to a value our own enumerator computed
(satisfies `feedback_golden_tests_sourced_only`). The lane "works" iff it
*surfaces* a family our continuation lane provably misses AND re-closes it to the
published state — a falsifiable, sourced, non-circular bar.

**Honesty on the golden's reach.** Passing it proves the *enumeration capability*
and *re-opens the dead region*; it does NOT, by itself, yield a catalogue cycler
(Png' is a PO, not a cycler). That is the correct, disciplined outcome: a
capability win with a published validation target, exactly the lever the strategy
draft ranked #1.

---

## 7. Risks + kill-criteria

| Risk | Mitigation / kill-criterion |
|---|---|
| **DA backend un-installable** (verified today: no daceypy/dace/pyaudi/MOSEK). | Plan Task 1 is a backend decision. The **sampling backend recovers Png' now** (the Zhou spike already did DA-method geometry by brute force). KILL the DA-acceleration sub-goal only if neither DACEyPy vendoring nor a pure-Python Taylor fallback is viable AND the sampling backend is too slow for the target bands — but the *capability + validation* ships regardless on sampling. |
| **Sampling backend too slow to sweep real bands** (O(n_grid · n_rev) propagations). | Acceptable for the *validation* (a handful of `(C,n)` points). For production sweeps, gate the DA acceleration behind a measured speedup target (paper: >84% CPU reduction). If sampling can't sweep a band in <1 day on the available cores, the production-sweep goal is DA-gated; the validation goal is not. |
| **Enumerator emits Png' only because the domain box was tuned to it** (circularity). | The base-family triangulation test (must also recover DRO + L1/L2 Lyapunov from the *same* untuned box) is the guard. KILL the "global" claim if the enumerator only finds targets it was boxed around. |
| **Recovered fixed points are all reproductions** (no genuinely new family in a *new* band). | That is a legitimate NEGATIVE for that band — re-stamp `empty_regions.jsonl` with the subsuming method+version. Not a lane failure; the lane's job is enumeration, emptiness is a finding. |
| **Floquet instability makes members un-exploitable** (P5g' max\|λ\|≈3600). | Expected for bridge orbits; recorded as a stability attribute, not a gate. Only matters if/when a member is a *cycler* candidate needing station-keeping (then M7/dv-band applies). |
| **Scope creep into a cycler claim.** | Hard boundary: a Png'-class PO never gets a catalogue row on PO-ness alone. The lit-novelty gate + the cycler-structure requirement (Earth–Earth resonant transfer leg) stay mandatory. |

---

## 8. PREREQUISITE corpus / environment acquisitions (honest flags)

1. **DA backend (environment, not corpus) — REQUIRED for the acceleration, NOT
   for the validation.** No `daceypy`/`dace`/`pyaudi` and no MOSEK in the env
   today (verified). Decision needed: (a) vendor DACEyPy + obtain a MOSEK
   academic licence (matches the paper exactly), OR (b) a pure-Python truncated
   Taylor-map + a non-commercial polynomial root-finder (no MOSEK), OR (c) ship
   on the sampling backend and defer DA. **The method is fully specified by the
   digested paper — this is an install/licence decision, not method invention.**
2. **No new corpus acquisition is strictly required for the lane.** The canonical
   reference arXiv:2509.12671 is **already digested** (CORPUS_INDEX line 154,
   `2026-06-13-high-order-transfer-map-2509.12671-mining.md`), the Png' ICs are
   transcribed, and one (P5g') is already re-closed under our corrector. The
   companion DA reachable-set paper (arXiv:2502.11280, Zhou-Armellin) is also
   digested with a working DA-geometry-by-sampling spike (`reachable_impulsive.py`)
   that is a direct precedent for the sampling backend.
3. **Optional precision upgrade:** if the printed Png' table digits prove
   precision-limiting for the lane-recovery tolerance, the paper's authors or the
   arXiv source code (DACEyPy is the cited library) could supply higher-precision
   ICs — but the existing 11-sig-fig re-closure of P5g' suggests the printed
   tables are already adequate.

**Net:** the lane is **buildable and validatable from currently-digested material
on the sampling backend now.** Only the paper's *performance* layer (DA Taylor
map) is environment-gated, and that gate blocks production-scale sweeps, not the
capability or its proof.

---

## 9. Confidence

- **Method buildable from digested material:** **HIGH** for the geometry +
  validation (sampling backend; Png' already re-closed; Zhou spike is the
  precedent). **MEDIUM** for the DA-accelerated production sweep (install/licence
  decision pending; method itself fully specified, not invented).
- **Re-opens a real dead region with a sourced validation target:** **HIGH** —
  the registry entry names this exact lever as its re-open key and contains the
  validation target *in the region*.
- **Yields a catalogue CYCLER:** **LOW and honestly so** — Png' is a PO not a
  cycler; the catalogue payoff is conditional on a future band yielding a fixed
  point that also carries cycler transfer structure. The lane makes that
  *searchable*; it does not promise it.

---

## 10. Open questions for the user

1. **DA backend choice (Task 1):** vendor DACEyPy + MOSEK (paper-exact, fastest,
   licence cost) vs pure-Python Taylor fallback vs ship-on-sampling-defer-DA.
   Recommendation: ship the **sampling backend + full validation first** (proves
   the lane, re-opens the region), then decide DA acceleration on a measured
   sampling-sweep cost.
2. **Production sweep scope after validation:** which bands first? Recommend the
   EM C≈3.0 band (has the target) then the three Saturn-moon Lyapunov bands (same
   re-open key, cheap).
3. **Catalogue policy for a novel non-cycler PO:** if the lane finds a *new* CR3BP
   PO family (not in JPL/literature) that is NOT a cycler, does it earn any
   record? (Recommendation: a negative-registry/repro-ledger note, not a
   catalogue row, consistent with the Png' boundary.)
