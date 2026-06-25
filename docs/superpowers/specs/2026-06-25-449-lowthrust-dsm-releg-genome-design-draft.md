# Low-thrust / DSM "releg" genome — re-open the ice/gas-giant moon-tour dead regions (DESIGN DRAFT, #449)

**Date:** 2026-06-25
**Status:** DESIGN-DRAFT — for user review. No production code written by this doc.
**Issue:** #449. Ranked **Region F / capability-lever #1** by the
discovery-strategy prioritization draft
(`2026-06-25-discovery-strategy-prioritization-design-draft.md`, §2 row F, §6
lever 1): *"A low-thrust / DSM releg genome — re-opens the entire ice-giant +
gas-giant moon-tour space … the structural emptiness is ballistic-only; the
fresh-but-infeasible Amalthea-class cycles could be rescued by a powered releg.
**Highest reach.**"*

This draft answers ONE question: **what is the smallest, well-bounded capability
that lets us re-solve the TRANSFER LEGS ("relegs") of an existing moon-tour
skeleton with a low-thrust or deep-space-maneuver arc — instead of a pure
ballistic Lambert leg — so that the V∞-continuity defect that kills the ballistic
tour becomes a BUDGETED ΔV the powered leg can absorb, and how does it plug into
the existing moon-tour genome / gauntlet / lit-check without rebuilding the
optimiser?**

---

## 0. TL;DR

- **Goal:** a *leg-swap* genome. A moon-tour cycle today is closed by re-solving a
  ballistic Lambert between each consecutive moon and requiring the V∞ to be
  *continuous* across each flyby (the unbridged |ΔV∞| is the closure residual,
  gated at 0.05 km/s). The releg genome replaces the per-leg `lambert(...)` call
  with a **powered leg solver** (a one-DSM-per-leg impulsive arc, OR a
  Sims-Flanagan low-thrust arc) that *delivers* a budgeted ΔV to make the leg
  close, and scores the tour on **total delivered ΔV per cycle** against the
  powered dv-band, not on a near-zero ballistic residual.
- **Re-opens (exact registry re-open keys):** the ice/gas-giant moon-tour
  negatives whose `interpretation` field names *low-thrust* or *DSM/powered*
  relegs as the re-open condition — see §1.2. Headline keys:
  `uranus-neptune-regular-moon-endgame-vilm-2026-06-23`,
  `jupiter-galilean-amalthea-repeated-moon-quasicycler-2026-06-24`,
  `jovian-IEG-vilm-2026-06-09`, `jovian-perm-vilm-2026-06-09`,
  `saturnian-titan-vilm-2026-06-09`, `saturnian-titan-endgame-vilm-2026-06-10`,
  and the five `repeated-moon-{jupiter,saturn,uranus,neptune,mars}-sweep`
  entries.
- **Reuse over rebuild — the optimiser ALREADY EXISTS.** The DSM leg solver
  (`search/dsm_leg.py`, `search/mga_dsm_placement.py`, #307) and the low-thrust
  Sims-Flanagan leg solver (`core/sims_flanagan.py` + `search/lowthrust.py`,
  #309) are both built, tested, and golden-anchored. The powered-releg *cost
  model* for the inter-moon transfer is **also already built and golden-validated**:
  `search/vilm.py::vilm_dv_min` reproduces Campagnola-Russell Part-1 Tables 1/2,
  and `europa_endgame_dv()` reproduces the Europa 154/147 m/s endgame. The new
  work is the **releg SWAP at one call site + a moon-tour-specific driver that
  loops the existing leg solver over a tour skeleton**, NOT a new transcription.
- **The exact swap seam (load-bearing):**
  `search/discovery_campaign.py::RepeatedMoonTarget._close_one_phasing`,
  **line 494** (`sols = lambert(r_a, r_b, tofs[k] * DAY_S, mu=mu, …)`) on the
  DISCOVERY side, and the mirror in
  `data/validation/v2_moontour.py::_cycle_residual`, **line 275** on the
  VALIDATION side. Both currently call `core.lambert.lambert`. The releg genome
  injects an alternative leg solver at exactly these two sites.
- **Validation / golden target (sourced, non-circular):** the
  Campagnola-Russell *Endgame Part-1* V∞-leveraging ΔV (already digested,
  `2026-06-05-endgame-tisserand-mining.md`, and already golden-validated in
  `vilm.py`). The releg genome's powered-leg ΔV for a moon-to-moon transfer must
  reproduce the published Table 1 (no-GA) / Table 2 (with-GA) ΔV_min to inside
  the 10% linked-conic-vs-CR3BP band — **AND** reproduce the *structural*
  Uranus/Neptune disjoint-contour finding (the VILM `min_vinf` floor forbids a
  bridge), proving the genome is honest about which dead regions stay dead.
- **#450 dependency finding: INDEPENDENT.** #449 relegs EXISTING tour skeletons
  (the closed-but-too-expensive ballistic cycles already enumerated in
  `repeated-moon-*-sweep` and the VILM negatives). It does NOT need #450's
  DA/HOTM enumerator as a seed source — #450 is a *CR3BP periodic-orbit*
  enumerator (its own draft, §4, is explicit that Png′ "is a PO, NOT a cycler"),
  not a moon-tour skeleton generator. The two are orthogonal capability levers.
  The historical "blocked by #450" tag is **incorrect** and should be dropped.
- **CRITICAL HONESTY — buildability:** **HIGH** confidence the genome is
  buildable from currently-digested material. Every primitive (DSM leg, SF leg,
  VILM cost model) exists, is tested, and traces to a digested+golden source. No
  new corpus acquisition is *required*. One OPTIONAL acquisition (a published
  *low-thrust* outer-moon tour with per-leg ΔV — Vasile-Campagnola 2009 DFET
  Europa SOT, partially transcribed) would strengthen the *low-thrust* branch's
  golden beyond the impulsive-VILM branch; see §8.

---

## 1. Which dead regions this re-opens — and why ballistic relegs fail there

### 1.1 Why a ballistic releg cannot close these tours

A moon-tour cycle is a chain `M0 → M1 → … → M0` where each leg is a ballistic
transfer between two moons of the same primary, and the only "free" energy is the
gravity assist at each intermediate flyby. The closure condition the current
genome enforces (`_close_one_phasing`, lines 503-528) is **V∞-continuity**: the
hyperbolic excess speed arriving at a flyby moon must equal the speed departing
it (a flyby rotates V∞ but cannot change its magnitude). The unbridged magnitude
defect `|V∞_in − V∞_out|` is the closure residual, gated at 0.05 km/s.

For the ice/gas-giant systems this residual is **structurally large** for two
distinct reasons recorded in the registry:

1. **High-V∞ basin (Jovian / Saturnian).** The no-leveraging ballistic corrector
   closes a *family*, but only at V∞ far above the feasible floor — gaps of
   **7.8–21.2 km/s** vs the 6.0 km/s floor (`jovian-*-vilm`, `saturnian-*-vilm`
   reverification). A flyby cannot remove that magnitude defect; a ballistic leg
   has no knob to spend.
2. **Disjoint Tisserand contours (Uranian / Neptunian — STRUCTURAL).** The
   regular moons are too widely-spaced and low-mass for *any* V∞ in 4–15 km/s to
   put consecutive moons on a shared resonance contour
   (`uranus-neptune-regular-moon-endgame-vilm-2026-06-23`: "every consecutive-moon
   leg unlinkable (disjoint Tisserand/resonance contours)"). Even the leveraging
   endgame (which walks V∞ *within* a contour) "cannot bridge disjoint contours."

A **powered releg** introduces the missing knob: a deterministic ΔV (impulsive
DSM, or distributed low-thrust) inside the leg that *changes* V∞ — exactly the
V∞-leveraging maneuver (VILM) physics, but generalised to an arbitrary in-leg
impulse rather than the apsidal-VILM special case. Reason (1) is then a *budget*
question (is the required ΔV within the powered dv-band?), not a feasibility wall.
Reason (2) is the honest caveat: a powered leg can climb V∞ to *reach* a contour,
but if two moons' contours are disjoint at *every* V∞ the leg has to cross a
forbidden region, and the cost may be prohibitive — the genome must report that
truthfully (see §6 kill-criterion and §7).

### 1.2 Exact re-open keys (from `data/empty_regions.jsonl`)

The registry has no dedicated `reopen` key; the re-open condition lives in each
entry's `verdict` / `interpretation`. The moon-tour negatives naming a
low-thrust/DSM releg are:

| `region_id` (verbatim) | Why ballistic fails | Re-open phrase (verbatim from `interpretation`/note) |
|---|---|---|
| `jupiter-galilean-amalthea-repeated-moon-quasicycler-2026-06-24` | fresh (Amalthea) cycles physically infeasible ballistically | *"Re-sweep warranted only if a NEW capability subsumes this: 3D/inclined relegs (Amalthea is inclined), **low-thrust relegs (could rescue the fresh-but-infeasible Amalthea cycles)**"* |
| `uranus-neptune-regular-moon-endgame-vilm-2026-06-23` | disjoint Tisserand contours (STRUCTURAL) | strategy-draft re-open key D3: *"Re-open only with a low-thrust or 3D/inclined releg genome (Region F)"* — note the structural caveat below |
| `jovian-IEG-vilm-2026-06-09` | high-V∞ basin, gap 7.8–21.2 km/s | strategy-draft §1.3 row 1 re-open: *"3D/inclined relegs, low-thrust"* |
| `jovian-perm-vilm-2026-06-09` | high-V∞ basin | same |
| `saturnian-titan-vilm-2026-06-09` | high-V∞ basin, gap ~13 km/s | same |
| `saturnian-titan-endgame-vilm-2026-06-10` | high-V∞ basin, gap ~6.8 km/s | same |
| `repeated-moon-jupiter-sweep` | empty as far as ballistic genome reached | strategy-draft §1.3 row 3 re-open: *"capability-subsuming genome"* |
| `repeated-moon-saturn-sweep` | empty | same |
| `repeated-moon-uranus-sweep` | empty | same |
| `repeated-moon-neptune-sweep` | empty | same |
| `repeated-moon-mars-sweep` | empty | same |

(The strategy-draft §1.3 negative-registry table is the canonical cross-index;
the `repeated-moon-*-sweep` re-open key "capability-subsuming genome" is satisfied
by a powered releg under the capability-subsumption rule because a powered leg
strictly subsumes the ballistic leg — the ballistic leg is the zero-ΔV limit.)

**Honest structural caveat on Uranus/Neptune.** The disjoint-contour negative is
the *one* moon-tour dead region where a low-thrust releg is **not guaranteed** to
help: if the contours are disjoint at every V∞, the powered leg must brute-force
across a region the leveraging physics cannot reach cheaply. This region is
re-opened in the sense that it *gets a powered re-test*, but the prior carries a
strong expectation of a *second clean negative* (now at the higher powered
dv-band). The genome re-stamps it with the subsuming method either way (a powered
emptiness is a stronger, more final negative than a ballistic one).

---

## 2. What a "releg" is, precisely

Define the moon-tour leg solver as a function with a fixed CONTRACT (this is the
seam — §4):

```
releg(r_a, v_a_moon, r_b, v_b_moon, tof_s, mu, n_rev) -> RelegResult
    inputs : planet-frame departure moon state (r_a, v_a_moon),
             planet-frame arrival moon state  (r_b, v_b_moon),
             leg time-of-flight, primary mu, revolution count
    outputs: vinf_out      (departure V∞ magnitude at moon A, km/s)
             vinf_in       (arrival   V∞ magnitude at moon B, km/s)
             dv_kms        (DELIVERED ΔV inside the leg, km/s; 0.0 for ballistic)
             feasible      (bool)
```

The **ballistic releg** (today, lines 494 / 275) is the `dv_kms == 0` special
case: it calls `lambert(...)`, picks the lowest-energy branch, and the
V∞-continuity defect across the flyby is the *uncompensated* error.

The **powered relegs** spend a budgeted `dv_kms`:

- **DSM releg (impulsive, the primary build).** One deep-space maneuver splits the
  leg into two Lambert sub-arcs joined at a free in-leg point — *exactly* the
  Vasile-Campagnola BS3 "one-DSM-per-leg" model (digest
  `2026-06-07-vasile-campagnola-dfet-method-mining.md`, §2.5: *"Each phase split
  into two Lambert subarcs joined by a deep-space maneuver"*). This is **already
  implemented**: `search/dsm_leg.py::dsm_leg(r0, v0, tof, eta, target_r, mu, …,
  max_revs)` propagates the front arc to the DSM point and solves a multi-rev
  Lambert from the DSM to the target, returning `dv_dsm_kms`. The releg wraps it:
  optimise the DSM fraction `eta` (and revolution branch) to *retarget the
  arrival V∞ to match the next leg's required departure V∞*, spending `dv_dsm_kms`.
  The leveraging primitive's analytic floor (`vilm.py::vilm_dv_min`) is the
  golden lower bound on this DSM ΔV.
- **Low-thrust releg (Sims-Flanagan, the secondary build).** The leg is a
  Sims-Flanagan N-segment arc with boundary states pinned to the two moon
  encounter states and the deliverable ΔV distributed across the thrust train —
  *already implemented*: `core/sims_flanagan.py::SimsFlanaganLeg` +
  `search/lowthrust.py::solve_leg_min_dv(leg)`. The SF leg's zero-thrust limit is
  the ballistic leg (digest cross-ref: *"Our SF leg's all-zero schedule =
  ballistic"*), so it strictly subsumes it. This branch is gated behind the DSM
  branch because SF is more expensive to converge and its golden is weaker (§8).

**Scoring change.** The ballistic genome minimises a near-zero residual. The
releg genome minimises **Σ dv_kms over the cycle's legs** (total delivered ΔV per
cycle) subject to: every leg feasible, V∞-continuity satisfied *after* the
powered retarget, and the closed-cycle wrap continuous. A tour is a hit if its
per-cycle ΔV lands inside the **powered dv-band** (`[300 m/s, 3.5 km/s × 7]` per
7 cycles, `verify/dv_band_acceptance.py`), not below the 0.05 km/s ballistic gate.

---

## 3. Architecture, new files, and the plug-in seam

The genome is **one new leg-solver module (the swap target) + one new moon-tour
driver + thin glue + one golden**. It reuses the DSM solver, the SF solver, the
VILM cost model, the moon registry, the dv-band gate, the gauntlet, and the
lit-check unchanged.

```
   tour skeleton (sequence, per-leg tof, n_rev, moon ephemerides)
   [from repeated-moon-*-sweep / VILM-negative cycles — EXISTING enumerations]
                              │
                              ▼
        ┌───────────────────────────────────────────────────────────┐
        │  search/releg_solver.py  (NEW — the swap target)           │
        │   Releg protocol + 3 backends:                             │
        │     BallisticReleg  -> core.lambert.lambert (dv=0; today)  │
        │     DsmReleg        -> search.dsm_leg.dsm_leg (#307)       │
        │     LowThrustReleg  -> core.sims_flanagan + search.lowthrust(#309)│
        │   Each retargets arrival V∞ to the next leg's needed V∞;   │
        │   returns RelegResult(vinf_out, vinf_in, dv_kms, feasible) │
        └───────────────────────────┬───────────────────────────────┘
                                    │ per-leg RelegResult
                                    ▼
        ┌───────────────────────────────────────────────────────────┐
        │  search/releg_moontour.py  (NEW driver)                    │
        │   - loop the tour's legs with the chosen Releg backend     │
        │   - enforce V∞-continuity AFTER powered retarget           │
        │   - sum delivered ΔV per cycle                             │
        │   - VILM-floor prefilter (cheap): skip legs whose          │
        │     vilm_dv_min already exceeds the powered dv-band         │
        └───────────────────────────┬───────────────────────────────┘
                                    │ cycle ΔV, per-leg ΔV, vinf chain
                                    ▼
   EXISTING  verify/dv_band_acceptance.py  (powered window gate)
             search/literature_check.py (V0, necessary-not-sufficient)
             data/validation/v2_moontour.py (V2 over n_cycles — releg-aware)
                                    │
                                    ▼
        SILVER holding tier (novel, powered) | reproduction | empty-restamp
```

**New files:**

| File | Responsibility | Reuses |
|---|---|---|
| `src/cyclerfinder/search/releg_solver.py` | The `Releg` protocol + `BallisticReleg` / `DsmReleg` / `LowThrustReleg` backends. One swappable seam. Each backend retargets arrival V∞ to a requested value and reports delivered ΔV. | `core.lambert.lambert`, `search.dsm_leg.dsm_leg`, `core.sims_flanagan`, `search.lowthrust.solve_leg_min_dv`, `search.vilm.min_vinf_for_vilm` |
| `src/cyclerfinder/search/releg_moontour.py` | Driver: take a tour skeleton + a `Releg` backend, close the cycle with powered legs, sum per-cycle ΔV, run the VILM-floor prefilter, classify. | `releg_solver`, `search.vilm.vilm_dv_min/vilm_dv_floor`, `core.satellites` |
| `data/golden/campagnola_endgame_releg.yaml` | Sourced golden: the Campagnola-Russell Part-1 Table 1/2 ΔV_min targets + the Europa 154/147 m/s endgame + the Uranus/Neptune disjoint-contour structural assertion. EXPECTED side traces only to the paper (`feedback_golden_tests_sourced_only`). | — |

**The plug-in seam (two call sites, identical signature).** The discovery-side
ballistic leg is `search/discovery_campaign.py:494` inside
`RepeatedMoonTarget._close_one_phasing`; the validation-side mirror is
`data/validation/v2_moontour.py:275` inside `_cycle_residual`. Both currently
call `lambert(r_a, r_b, tof*DAY_S, mu=mu, max_revs=nrev)` and derive V∞ from the
returned branch. The minimal-footprint change: both methods already accept (or
trivially can accept) an injected `lambert`/leg-solver callable
(`_close_one_phasing` already takes `lambert: Any` as a parameter, line 463 — the
seam is *already* dependency-injected on the discovery side). The releg genome
passes a `Releg` backend whose `BallisticReleg` is call-compatible with today's
behaviour (regression-preserving) and whose `DsmReleg`/`LowThrustReleg` add the
powered path. **No optimiser is written; no corrector is changed.**

**Catalogue-scope boundary.** A powered moon-tour that closes inside the powered
dv-band is a `quasi_cycler` or `mga_tour` class candidate (per
`project_catalogue_scope_expanded_2026-06-15`), NOT an automatic catalogue row —
it must still clear V0 lit-check + the powered V2/V3 gauntlet + the mandatory
literature-novelty gate. The Callisto-Ganymede-Callisto known-attribution V0 note
(`jupiter-galilean-amalthea-…` result) is the boundary exemplar: a *known*
powered tour is a V0-known admission, a *novel* one is a SILVER holding-tier
candidate pending human promotion.

---

## 4. Data flow

1. **Input:** a tour skeleton `(primary, sequence, per-leg tof, n_rev,
   moon-ephemeris phasing)` taken from the EXISTING enumerations that produced
   the negatives — `repeated-moon-*-sweep` closed-but-too-expensive cycles and
   the VILM-negative high-V∞ families. (#449 does **not** generate skeletons; it
   relegs given ones — see §5 #450 finding.)
2. **VILM-floor prefilter (cheapest first, the C21 lesson).** For each leg,
   compute `vilm.vilm_dv_min(moon_a, moon_b)`. If the *analytic* leveraging floor
   already exceeds the powered dv-band, the leg is unbridgeable even by an optimal
   powered leg — skip the whole skeleton (this is the Uranus/Neptune disjoint
   case: `min_vinf_for_vilm` shows no efficient V∞ exists). This avoids running
   the expensive DSM/SF solve on structurally-dead legs.
3. **Per-leg powered solve.** For surviving legs, the chosen `Releg` backend
   retargets the arrival V∞ to the next leg's required departure V∞ (so the flyby
   continuity is satisfied *after* the maneuver), returning delivered ΔV.
4. **Cycle close + wrap.** Enforce V∞-continuity at every interior flyby and the
   closed-cycle anchor wrap (same definition as today's `_close_one_phasing`),
   now post-retarget. Sum delivered ΔV.
5. **Powered dv-band gate** (`verify/dv_band_acceptance.py`): the per-7-cycle ΔV
   must land inside `[300 m/s, 3.5 km/s × 7]`.
6. **Novelty pipeline:** `spatial`/signature dedup → `literature_check`
   (V0, necessary-not-sufficient) → V2 moontour (releg-aware) → SILVER holding.
7. **Output:** classified ledger — `{reproduction | known-powered-tour |
   novel-powered-tour-candidate | empty}`. An empty band **re-stamps** the
   corresponding `empty_regions.jsonl` entry with the subsuming powered method +
   version (capability-subsumption record), never silent deletion.

---

## 5. #450 dependency finding (explicit)

**Finding: #449 is INDEPENDENT of #450. The historical "blocked by #450" tag is
incorrect and should be removed.**

Evidence:

- **#450 is a CR3BP *periodic-orbit* enumerator, not a moon-tour skeleton
  source.** Its own design draft
  (`2026-06-25-450-da-hotm-enumeration-lane-design-draft.md`, §4
  catalogue-scope-boundary, verbatim): *"Png′ is a CR3BP periodic orbit, NOT a
  cycler (no Earth–Earth resonant transfer leg…). It is a cross-check /
  validation target only — never a catalogue row."* #450 emits
  `(x0, xdot0, C, n)` fixed-point ICs for the `correct_general_periodic`
  corrector — a *rotating-frame periodic orbit*, a completely different object
  from a planet-frame moon-tour leg chain.
- **#449's seeds already exist.** The skeletons #449 relegs are the
  closed-but-too-expensive ballistic cycles ALREADY enumerated by the
  repeated-moon sweeps (`repeated-moon-*-sweep`: e.g. 20,488 closed Saturn cycles)
  and the VILM negatives. #449 needs no new global enumerator — it re-scores
  existing tour geometry with a powered leg.
- **The two levers are orthogonal in the strategy draft** (§6): lever 1 is "a
  low-thrust / DSM releg genome" (this, #449); lever 2 is "a global multi-rev
  fixed-point enumeration lane (DA/HOTM)" (#450). They re-open *different* dead
  regions (moon-tours vs the EM C≈3.0 CR3BP-PO band) and share no data path.

**Conclusion:** #449 can be built and validated with #450 absent. The only thing
that would couple them is a *future* desire to releg a moon-tour skeleton that
itself came from a CR3BP-PO source — not in scope and not a dependency.

---

## 6. Validation strategy (the sourced/golden proof the releg works)

**Primary golden (decisive, already in-repo & golden-validated).** The powered
inter-moon releg ΔV must reproduce the Campagnola-Russell *Endgame Part-1*
V∞-leveraging cost — the published, digested
(`2026-06-05-endgame-tisserand-mining.md`) reference for powered moon-to-moon
transfer ΔV:

- **Leveraging-floor golden (regression lock).** For a moon-A → moon-B powered
  releg with the leg's required ΔV∞, the `DsmReleg` delivered ΔV must be **≥** the
  analytic VILM floor `vilm.vilm_dv_min(A, B)` and reproduce Part-1 Table 1
  (no-GA) / Table 2 (with-GA) ΔV_min to inside the 10% linked-conic-vs-CR3BP band
  — *this golden already passes* in `vilm.py` (docstring: *"GOLDEN: reproduces
  Part-1 Table 1 (no-GA) and Table 2 (with-GA) ΔV_min to well inside the 10%
  band"*). The releg's job is to *deliver* a ΔV consistent with that floor; the
  test asserts `dv_kms ≥ vilm_dv_min − tol` and within the band.
- **Europa-endgame golden.** A Ganymede/Callisto → Europa powered releg's
  delivered ΔV must bracket the published Europa endgame **154 m/s (discrete
  3-VILM) / 147 m/s (CR3BP re-optimised)** with the continuous floor
  (`vilm.europa_endgame_dv()` returns the Eq.13 floor + 46-day phasing) as the
  lower bound. EXPECTED side traces to the Part-1 A6 mining note.
- **Structural-emptiness golden (the honesty test).** The driver, run on a
  Uranian Ariel→Umbriel leg, must report **unbridgeable** because
  `min_vinf_for_vilm` shows no efficient V∞ exists / the contours are disjoint —
  reproducing `uranus-neptune-regular-moon-endgame-vilm-2026-06-23` as a powered
  re-test that *also* finds it empty. This proves the genome does not fabricate a
  bridge where the physics forbids one.
- **Positive-control golden.** Run on a Jovian Io-Europa-Ganymede or Saturnian
  Dione-Rhea-Titan skeleton (the registry's own positive controls: "both link at
  vinf=4"), the DSM releg must close the cycle inside the powered dv-band — the
  *capability proof* that a powered leg rescues a ballistically-too-expensive
  cycle.

**Why golden-clean.** Every EXPECTED value traces to Campagnola-Russell's printed
tables / the digested mining note — never to a number the releg solver itself
computed (satisfies `feedback_golden_tests_sourced_only`). The genome "works" iff
it (a) delivers ΔV consistent with the published leveraging floor on a feasible
leg, (b) closes a positive-control cycle inside the powered band, and (c) honestly
reports the disjoint-contour legs as unbridgeable — a falsifiable, sourced,
non-circular bar.

**Honesty on the golden's reach.** Passing it proves the *powered-releg
capability* and re-opens the gas-giant high-V∞ moon-tour regions with a budgeted
ΔV. It does NOT, by itself, yield a *novel* catalogue cycler — the only
literature-fresh outer-moon cycles found ballistically (the Amalthea class) are
the speculative prize, and whether a powered leg makes them land inside the band
is the open empirical question the campaign (a separate issue) answers.

---

## 7. Risks + kill-criteria

| Risk | Mitigation / kill-criterion |
|---|---|
| **DSM releg just rediscovers the ballistic high-V∞ family with a token ΔV** (no real V∞ reduction). | The releg must *retarget* arrival V∞ to the next leg's needed value; assert the post-maneuver V∞-continuity residual drops below the ballistic gate AND the delivered ΔV ≥ the VILM floor. KILL the "powered helps" claim for a region if the only closures need ΔV above the 3.5 km/s/cycle powered ceiling (then it is a *powered* empty — re-stamp the registry). |
| **Uranus/Neptune stays empty even powered** (disjoint contours, expected). | This is a legitimate, stronger NEGATIVE — re-stamp `uranus-neptune-regular-moon-endgame-vilm` with the powered method+version. NOT a lane failure; emptiness is a finding. The VILM-floor prefilter (§4.2) makes this cheap. |
| **SF low-thrust branch too slow / non-convergent** on multi-leg tours. | Gate the SF branch behind the DSM branch. The DSM branch (impulsive, one Lambert pair per leg) carries the *capability proof and the golden*; SF is an optional refinement whose convergence is a known cost (it is the #309 leg solver, already characterised). If SF cannot close a positive-control cycle in bounded time, ship on DSM only and defer SF. |
| **Powered tour "closes" but is not a cycler** (no repeat structure). | The V2-moontour gate (≥3 cycles, drift-bounded) + the closed-cycle wrap continuity stay mandatory; a single powered leg-close is necessary-not-sufficient. |
| **Scope creep into a novelty claim from a known powered tour.** | Hard boundary: the lit-novelty gate (§16.5) + the dv-band classification stay mandatory. Callisto-Ganymede-Callisto is a V0-known admission, not a novelty. |
| **The two swap sites drift out of sync** (discovery vs validation leg models differ). | The `Releg` protocol is the single source of truth; both sites consume the same `BallisticReleg`/`DsmReleg`. A regression test asserts `BallisticReleg` reproduces today's `lambert`-path V∞ exactly (bit-for-bit on a fixed skeleton). |

---

## 8. Prerequisite corpus / environment acquisitions (honest flags)

1. **No acquisition is REQUIRED.** Every primitive is built, tested, and
   golden-anchored to digested material:
   - DSM leg: `search/dsm_leg.py` (#307, Takao 2025 transcription) + tests.
   - SF low-thrust leg: `core/sims_flanagan.py` + `search/lowthrust.py` (#309,
     Yam/Ozimek transcription) + tests.
   - Powered-releg cost model + golden: `search/vilm.py::vilm_dv_min` /
     `europa_endgame_dv` (Campagnola-Russell Endgame Part-1, digested
     `2026-06-05-endgame-tisserand-mining.md`, golden-validated to <10% band).
   - Moon registry, dv-band gate, gauntlet, lit-check: all present.
2. **OPTIONAL strengthener for the LOW-THRUST branch only.** The impulsive-VILM
   golden (above) anchors the *DSM* branch tightly. The *SF low-thrust* branch's
   golden is currently weaker — the best digested low-thrust outer-moon tour with
   per-leg ΔV is **Vasile-Campagnola 2009 "Low-Thrust MGA to Europa" (DFET)**,
   PARTIALLY transcribed (`2026-06-05-vasile-tables-retranscription.md`, Tables
   5/6 Ganymede/Europa SOT). Its DFET transcription *does not map* to our SF leg
   (digest §2.6: *"DOES NOT MAP — a DIFFERENT transcription"*), so it is a *target
   number to bracket*, not a state-level golden. If the SF branch is built and a
   tighter low-thrust golden is wanted, completing the Vasile Tables 7/8 ΔV
   transcription (already started) is the cheapest acquisition — but it is
   OPTIONAL and the DSM branch ships without it.
3. **No environment/install gate** (unlike #450's DA/MOSEK gate). The releg
   genome is pure-Python on the existing scipy stack.

**Net:** the releg genome is **buildable and validatable from currently-digested
material now**, on the DSM (impulsive) branch, with the SF (low-thrust) branch as
an optional second backend whose tight golden is the one (optional) acquisition.

---

## 9. Confidence

- **Genome buildable from digested material:** **HIGH** for the DSM branch (every
  primitive built + tested; the cost-model golden already passes in `vilm.py`;
  the swap seam is already dependency-injected on the discovery side).
  **MEDIUM-HIGH** for the SF low-thrust branch (the #309 leg solver exists but its
  outer-moon-tour golden is bracket-only, not state-level).
- **Re-opens real dead regions with a sourced validation target:** **HIGH** — the
  registry entries name low-thrust/DSM relegs as their explicit re-open key, and
  the Campagnola-Russell leveraging golden is digested + already validated.
- **Independent of #450:** **HIGH** — #450 is a CR3BP-PO enumerator by its own
  draft; no shared data path; the "blocked by #450" tag is wrong.
- **Yields a NOVEL catalogue cycler:** **LOW-MEDIUM and honestly so** — the
  high-V∞ gas-giant regions become *budget* questions a powered leg can answer,
  and the positive controls will close; but the only *literature-fresh* outer-moon
  cycles (Amalthea class) landing inside the powered band is the speculative prize
  the follow-on campaign tests. Uranus/Neptune is expected to stay empty (a
  stronger powered negative). The capability win is near-certain; the novel-cycler
  payoff is the open bet.

---

## 10. Open questions for the user

1. **Build order:** DSM branch first (carries the golden + capability proof),
   defer SF? Recommendation: **yes** — DSM is fully golden-anchored and pure
   impulsive; SF is the optional refinement.
2. **Skeleton source:** releg the existing `repeated-moon-*-sweep`
   closed-but-too-expensive cycles, or re-enumerate skeletons with a powered-aware
   prefilter? Recommendation: **releg existing** (cheaper, and the negatives are
   exactly those skeletons).
3. **Catalogue policy for a powered tour:** a novel powered moon-tour that closes
   inside the dv-band is a `quasi_cycler`/`mga_tour` SILVER candidate — confirm it
   follows the same V0-known vs SILVER-novel split as the Callisto-Ganymede-Callisto
   exemplar.
4. **Drop the "blocked by #450" tag on #449** (§5 finding) — confirm.
