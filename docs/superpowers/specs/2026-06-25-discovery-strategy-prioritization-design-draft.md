# Discovery-strategy prioritization — where to look for novel cyclers FIRST (DESIGN DRAFT)

**Date:** 2026-06-25
**Status:** DESIGN-DRAFT — for user review. No code, no catalogue writeback, no campaign launched by this doc.
**Scope:** A prioritized, falsifiable map of `(search-DIMENSION × BODY-SYSTEM × ENERGY/RESONANCE)` regions, scored on `P(novel) · P(usable cycler)` against capability-readiness and cost. Grounds every claim in the negative-results registry (`data/empty_regions.jsonl`, 36 entries) and the catalogue's real body coverage.

This draft answers ONE question: **given a finite search budget, which region has the highest probability of yielding a GENUINELY NOVEL, USABLE cycler?** It deliberately rejects re-mining and names the dead axes so budget is not wasted.

---

## 0. TL;DR — the recommendation up front

1. **#1 search target: the ISOLATED-ER3BP elliptic-family cycler-usability probe, run as the gated phased plan already adjudicated in `isolated-er3bp-cycler-novelty-gap-analysis-2026-06-25`.** It is the *only* region in the whole map where the DYNAMICAL family is provably *unmapped* for our µ regime (Antoniadou & Libert 2018 compute isolated high-`e` families with **no circular limit**, but nobody has done so at a cycler-relevant µ or asked whether they carry an Earth↔Moon transfer leg). Every other "open" region is either published-class, structurally empty, or a refinement of a known family. **But it is a GATED GO, not a green light:** the cycler prize is speculative (semantic-mismatch risk), so it runs as a bounded capability gate FIRST (reproduce a µ=0.001 isolated family) and a separate speculative cycler pass SECOND, with a hard kill-criterion.

2. **#2 (parallel, cheap): cycler-USABILITY characterization of the already-found STABLE 3D C21/C32 families** (`3d-broken-plane-em-4cycler-lift-2026-06-25`). These are NOT novel *dynamical families* (axis closed — see Dead Axis D1), but their *transport utility* (V∞, encounter geometry as low-thrust/ballistic transfer infrastructure) is uncharacterized. This is a usability question on known structure, not a discovery — lower novelty ceiling but near-zero marginal cost (families already computed, gauntlet partly built).

3. **#3: a genuinely-new µ/resonance combination absent from the planetary-systems literature** — specifically a high-`e` Sun-planet ER3BP DRO/Lyapunov-rooted cycler at **Mercury (e=0.206) or Mars (e=0.093)**, the regime the `er3bp-discovery-em-broucke-koblick` campaign explicitly SKIPPED (no CR3BP seed; ICs were not fabricated). This is the one ER3BP corner where the model departs most from CR3BP and where the negative is only *conditional*.

**Main dead-axis warnings (do NOT spend budget here):**
- **3D lifts of known planar cyclers (Aldrin / Braik-Ross C11/C21/C32) are PUBLISHED-CLASS** (Antoniadou-Voyatzis-Libert spatial-resonant lineage). Axis closed for novelty (`3d-extension-of-known-planar-cyclers-AXIS-published-class-2026-06-25`).
- **The planar CR3BP envelope is mined** — single-family Jacobi continuation surfaces only reproductions / same-family points (9 seeds, 0 distinct new families).
- **Ballistic moon-tours in the ice-giant and most gas-giant systems are STRUCTURALLY empty** — Uranian/Neptunian regular moons have disjoint Tisserand contours (unlinkable at any V∞); Jovian/Saturnian feasible closures are published resonant tours.
- **ER3BP-only families branching off Earth-Moon CR3BP roots do NOT exist up to e=0.15** (3 mechanisms exhausted; the families just persist).
- **Classic Earth-Mars precursor-MGA space is literature-dense** — the global engine reaches it but finds 0 literature-fresh insertions.

---

## 1. The grounding facts (what is MINED vs OPEN)

### 1.1 Catalogue body coverage (from `data/catalogue.yaml`, 319 rows)

| Body-set | Rows | Source cluster | Density |
|---|---:|---|---|
| Earth–Mars (`E,M`) | 269 | Russell 2004 (~222), Rall 1970, Rogers 2012, McConaghy, Aldrin | DENSE (census), but seed-level — most are summary-table family members, not reproducible ICs |
| Earth–Venus (`E,V`) | 15 | Hollister-Menning 1970 | THIN (bare canonical orbits) |
| Earth–Moon (`E,Moon`) | 12 | Ross-Braik / Roberts-Tsoukkas 2025-26 | THIN but recent; the only V2 cluster |
| Earth–Mars–Venus (`E,M,V`) | 4 | Jones 2017 VEM-triple | THIN |
| Jovian moons (CGE / IEG) | ~7 | Liang 2024, Russell-Strange 2009 | family-seed placeholder, null-numeric |
| Saturnian / Uranian / Neptune moons | ~5 | Russell-Strange, Heaton-Longuski | family-seed placeholder |
| Spacecraft grand-tours | ~8 | Voyager/Cassini/Galileo/Juno | historical, not cyclers |

Validation tiers: V4:1, V3:2, V2:7, V1:22, V0/unspecified: rest. **The catalogue is at its data-limited validation ceiling** (`project_validation_ceiling`): the ~200 ocampo V0 rows are V0 by a *publication gap* (summary-only `n.m.k` format), not laziness. Past the ceiling is **NEW-INPUT-gated**, not iteration-gated.

**Implication for discovery:** Earth-Mars is the most-published interplanetary cycler system on Earth; novelty there is structurally improbable. Beyond-E-M moon systems are thin in the catalogue precisely because the ballistic genome found them empty/structural (see 1.3). The OPEN territory is therefore *new dimensions* (ER3BP isolated families, 3D, QP) and *new µ/resonance combinations*, NOT new rows in the dense E-M heliocentric envelope.

### 1.2 Capability readiness (genomes + gauntlets, as code TODAY)

| Genome | Module(s) | State |
|---|---|---|
| Planar CR3BP correctors | `search/cr3bp_general_periodic.py`, `nrho_continuation.py`, `cr3bp_jacobi_arclength.py` | COMPLETE |
| 3D CR3BP family tracer | `search/cr3bp_general_periodic_3d.py`, `cr3bp_3d_family_tracer.py` | COMPLETE (Phase 1+2) |
| ER3BP | `genome/er3bp_periodic.py`, `er3bp_continuation.py`, `er3bp_branching.py`; `search/er3bp_discovery.py`, `er3bp_direct_seeding.py`, `er3bp_isolated_seeds.py` | COMPLETE; **isolated-seed lane is the freshest, least-exercised path** |
| BCR4BP / QBCP | `genome/bcr4bp_genome.py`, `bcr4bp_continuation.py`, `bcr4bp_systems.py` | COMPLETE (libration families only) |
| QP tori (GMOS) | `genome/qp_tori.py` (+ `qp_tori_continuation.py`) | Single-torus COMPLETE; pseudo-arclength family continuation IN-PROGRESS (#333 draft) |
| Epoch-aware MGA / precursor | `genome/epoch_aware_genome.py`, `search/global_precursor_engine.py`, `free_return_chain.py`, `multiarc_closure.py` | COMPLETE |
| Repeated-moon / VILM / endgame | `search/moon_cycler_genome.py`, `vilm.py`, `endgame_graph.py` | COMPLETE |

| Gauntlet tier | Planar/heliocentric | 3D CR3BP-PO | QP tori | BCR4BP |
|---|---|---|---|---|
| V0 (novelty: lit-check + spatial prefilter) | YES | YES (`spatial_novelty_prefilter`) | partial (no QP fingerprint) | n/a |
| V1 | YES | YES (`v1_3d.py`) | YES (`v1_qp.py`) | MISSING (only in-corrector dual-closure) |
| V2 | YES | YES (`v2_3d.py`) | YES (`v2_qp.py`) | MISSING |
| V3 | YES (dv-band) | PARTIAL (`v3_3d_periodic.py` exists; moontour `v3_3d` does not apply to a PO — #306 gap) | MISSING | MISSING |
| V4/V5 | partial (Hill-screen pre-V4; GMAT lane uncommitted) | as planar | MISSING | MISSING |

**Readiness summary:** the **3D and ER3BP lanes are the most gauntlet-ready discovery surfaces**; QP and BCR4BP need gauntlet builds (`#305`, `#333` drafts) before any candidate could be promoted. This makes QP/BCR4BP poor *first* targets regardless of dynamical promise — a find there cannot be validated yet.

### 1.3 Negative registry — what is RULED OUT (anti-catalogue, 36 entries)

| Ruled-out region | Verdict | Re-open condition |
|---|---|---|
| Jovian / Saturnian ballistic moon-tours (VILM, IEG/perm) | EMPTY (high-V∞ basin, gap 7-21 km/s) | 3D/inclined relegs, low-thrust |
| Uranian + Neptunian regular-moon tours | EMPTY (STRUCTURAL — disjoint Tisserand contours) | none ballistic; low-thrust only |
| Repeated-moon multi-rev (Jup/Sat/Ura/Mars/Nep sweeps) | empty | capability-subsuming genome |
| Jupiter Galilean+Amalthea quasi-cyclers | NEGATIVE (fresh⇒infeasible; feasible⇒published) | 3D/inclined or low-thrust relegs |
| CR3BP single-family Jacobi continuation (Ross 1,1/2,1/3,x; Saturn moons; Arenstorf) | EMPTY (reproductions / same-family only) | DA/HOTM global multi-rev enumeration (Png′) |
| EM C_J≈3.00 DRO/Lyapunov new-family | EMPTY (method-limit, NOT void — Png′ exists) | global enumeration lane |
| Binary-star µ-continuation (EM→µ0.3, µ0.5) | EMPTY (family folds before target µ) | — |
| Precursor-MGA into Aldrin/S1L1 (local + global engine) | NOVELTY EMPTY (0 literature-fresh) | Earth-launch enforced, wider epoch search |
| ER3BP off-EM-CR3BP-roots up to e=0.15 | NEGATIVE (conditional — families persist, no branch) | **high-e Sun-planet seeds (SKIPPED)**, branch-switching |
| ER3BP high-e Sun-planet Lyapunov/DRO continuation | NEGATIVE (conditional — 6/6 survive, no e>0-only) | direct-e>0 seeding, robust tracker |
| ER3BP direct-e>0 blind grid | NEGATIVE (methodological — discriminator unreliable) | fold-aware pseudo-arclength tracker |
| 3D lift of known planar cyclers (C11/C21/C32) | **AXIS CLOSED for novelty** (published Antoniadou lineage) | — (closed) |
| Stable 3D C21 (2,0,10) candidate | NOT novelty-claimable (known spatial 2:1-resonant class) | — (usability remains open, see Region B) |
| Isolated-ER3BP cycler novelty-gap | **GATED GO** (genuine dynamical gap; speculative cycler prize) | runs as phased plan — Region A |

**The registry's own redirect** (`3d-extension...AXIS-published-class`, verbatim): the genuine novel-cycler frontier is *(a)* ISOLATED elliptic families with no circular/planar limit; *(b)* cycler-USABILITY where the dynamical family is known but transport utility is not; *(c)* genuinely new µ/resonance combinations absent from the planetary-systems literature. **The three top regions below map exactly onto (a), (b), (c).**

---

## 2. The prioritization matrix

Scoring, each 1 (low) – 5 (high). `P(novel)` = not in registry/corpus. `P(usable)` = low-V∞ / good encounter geometry vs pure dynamical PO. `Ready` = genome+gauntlet exist. `Cost` is inverted in the score (cheaper = better). **Priority = P(novel) · P(usable) · Ready / Cost**, then judgment.

| # | Region (DIM × SYSTEM × ENERGY/RESONANCE) | P(novel) | P(usable) | Ready | Cost | Verdict |
|---|---|:--:|:--:|:--:|:--:|---|
| **A** | **ER3BP ISOLATED elliptic family × cycler-relevant µ × high-e {3/2,5/2,3/1,4/1,5/1}** | **5** | 2 | 4 | high | **#1 — only provably-unmapped dynamical region; gated phased plan** |
| **B** | **3D stable C21/C32 USABILITY × Earth–Moon × 2:1/3:2 resonance** | 2 (family known) | **4** | 4 | **low** | **#2 — usability of known structure; cheap, parallel** |
| **C** | **ER3BP DRO/Lyapunov-rooted × high-e Sun-planet (Mercury 0.206, Mars 0.093) × resonant** | 4 | 3 | 4 | med | **#3 — the SKIPPED conditional corner** |
| D | QP-torus cyclers × Earth–Moon × Floquet-unit-circle Braik-Ross members | 4 | 2 | 1 (no V3+_qp gauntlet) | high | DEFER — dynamically interesting, un-validatable today |
| E | BCR4BP synodic-locked cyclers × cislunar × Sun-resonant | 3 | 2 | 1 (no gauntlet) | high | DEFER — no cycler prior; libration-only genome |
| F | Low-thrust moon-tours × ice/gas-giant × resonant (re-open structural negatives) | 4 | 4 | 1 (no low-thrust genome) | very high | NEW-INPUT-gated — needs a low-thrust releg genome first |
| G | Global multi-rev (DA/HOTM) CR3BP enumeration × EM × C≈3.0 band | 4 | 2 | 1 (no DA lane) | very high | NEW-INPUT-gated — would re-open Png′ region |
| — | 3D lift of known planar cyclers | 1 | 3 | 5 | low | **DEAD — published-class** |
| — | Planar CR3BP single-family continuation | 1 | 3 | 5 | low | **DEAD — mined** |
| — | Ballistic ice-giant / gas-giant moon-tours | 1 | 2 | 5 | low | **DEAD — structural** |
| — | Classic E-Mars precursor-MGA | 1 | 5 | 4 | med | **DEAD — literature-dense, 0 fresh** |
| — | ER3BP off-EM-CR3BP-roots (e≤0.15) | 1 | 2 | 4 | med | **DEAD (conditional) — families persist** |

---

## 3. Top regions to search FIRST

### Region A (#1) — ER3BP isolated elliptic-family cycler probe

**Why likely novel.** The registry entry `isolated-er3bp-cycler-novelty-gap-analysis-2026-06-25` adjudicates this as a **GATED GO**: the dynamical gap is *genuine* — isolated elliptic resonant families at µ_EM for {3/2, 5/2, 3/1, 4/1, 5/1} are computed by **no published work**. Antoniadou & Libert (2018, DOI 10.1007/s10569-018-9834-8) prove these families *exist* (high-`e`, stable, no circular limit) but only at µ=0.001 (exoplanetary) and only as inner-body-vs-distant-perturber MMRs — never at a cycler µ, never asked about a transfer leg. Reproducing them at a cycler-relevant µ is a legitimate CAPABILITY win, NOT a C21-style rediscovery. This is the one region where the dynamical family itself is unmapped, not merely its utility.

**Why it could be a usable cycler.** Speculative, and this is the honest weak point. The isolated families are *stable* and *highly eccentric* — eccentricity is exactly what a cycler's Earth↔Moon (or Earth↔Mars) transfer leg needs. IF an isolated high-`e` member at the right µ carries an encounter geometry linking the two primaries, it would be a structurally new cycler class with no CR3BP ancestor.

**Capability that runs it.** `genome/er3bp_isolated_seeds.py` + `er3bp_periodic.py` + the fold-aware `er3bp_continuation.py` pseudo-arclength tracker (the `er3bp-direct-e0-blind-grid` negative proved the *blind-grid secant* discriminator is unreliable — the fold-aware tracker is the fix). Validation: ER3BP candidates lift into the 3D/ heliocentric gauntlet via V0 lit-check (Antoniadou corpus already anchored in `known_corpus_3d`) + V1/V2.

**Falsifiable kill-criterion (verbatim from the adjudication).** Phase the work:
1. **Bounded capability gate:** reproduce ONE µ=0.001 isolated family (digitize Antoniadou 2018 high-`(e1,e2)` config — the precision-limiting step) and confirm the isolated signature via reverse-continuation *death before e=0* (no circular limit). If it cannot be reproduced/confirmed isolated, STOP — the discriminator is not trustworthy.
2. **Speculative cycler pass (separate):** direct high-`e` EM seeding at the digitized config. **KILL:** if it dies into a distant-perturber resonant PO with **no Earth↔Moon transfer leg**, that is a clean negative — **close the axis**.

This is the C21 lesson applied *at scale* before committing weeks: a pre-build novelty check on a possibly-published axis.

### Region B (#2) — Cycler-USABILITY of the stable 3D C21/C32 families

**Why likely novel — and the honest caveat.** The *dynamical family* is NOT novel (Dead Axis D1 — `3d-extension...AXIS-published-class`). What is uncharacterized is the **transport utility** of the already-computed STABLE out-of-plane C21 (2,0,10) and C32 members (`3d-broken-plane-em-4cycler-lift-2026-06-25`, 2260 closure-verified members). The registry's own redirect names "cycler-USABILITY in regions where the dynamical family is known but its transport utility is not characterized" as a genuine frontier. So the novelty here is *engineering* (a usable stable spatial cycler), not *dynamical discovery*.

**Why it could be a usable cycler.** Stable 3D cyclers are *rare* — most cyclers are strongly unstable (the Ross EM members run ~2000×/period). A stable out-of-plane member with bounded station-keeping is operationally valuable even if its family is published. The question: does the stable C21 family have low maintenance-ΔV and feasible Earth/Moon encounter geometry?

**Capability that runs it.** Families already exist (`cr3bp_3d_family_tracer`). Needs: the **3D V3-periodic gauntlet** (`v3_3d_periodic.py` exists; the #306 draft completes the periodic-orbit IAS15 cross-check) + the M7 maintenance-ΔV / dv-band acceptance on the 3D members. Near-zero marginal discovery cost.

**Falsifiable kill-criterion.** Run the 3D V3 + dv-band gauntlet on the stable C21/C32 members. **KILL:** if maintenance-ΔV exceeds the powered ceiling (3.5 km/s/cycle) OR encounter geometry has no low-V∞ Earth/Moon transfer window, the families are pure dynamical PObjects with no transport utility — log as characterized, do not promote.

### Region C (#3) — High-`e` Sun-planet ER3BP cyclers (the SKIPPED corner)

**Why likely novel.** The `er3bp-discovery-em-broucke-koblick` negative is **explicitly conditional**: it covered only Earth-Moon seeds and **SKIPPED the high-`e` Sun-planet systems (Mercury 0.206, Mars 0.093, Pluto 0.249) where the ER3BP departs MOST from the CR3BP** — because no CR3BP seed IC existed and ICs were (correctly) not fabricated. The negative's own re-sweep condition is "high-e Sun-planet seeds (generate CR3BP Lyapunov/DRO at those µ first)." This is unsearched, not ruled out.

**Why it could be a usable cycler.** Mars (e=0.093) is the canonical cycler destination; an ER3BP Mars cycler that exploits real Mars eccentricity could be lower-V∞ than its circular-model approximation (the eccentricity is a free energy source at the right phasing).

**Capability that runs it.** Generate CR3BP Lyapunov/DRO seeds at Sun-Mars / Sun-Mercury µ (existing planar correctors), then continue in `e` via `er3bp_continuation.py` to the real planetary eccentricity, watching for bifurcation (`er3bp_branching.py`). V0 lit-check + V1/V2.

**Falsifiable kill-criterion.** Continue 3 Sun-Mars and 3 Sun-Mercury CR3BP families to their real `e`. **KILL:** if all families persist with no bifurcation and no V∞ reduction vs the circular model (exactly the EM outcome), the ER3BP-departure hypothesis is falsified for cyclers — close the high-e ER3BP corner. (Note: this is the weaker of the three #1-3 — it may well reproduce the EM persistence result. Run it AFTER A and B.)

---

## 4. DEAD AXES — do NOT spend budget here

- **D1. 3D lifts of known planar cyclers (Aldrin / Braik-Ross C11/C21/C32).** Lifting a known planar root out of plane lands by construction on its published spatial-bifurcation family (Antoniadou-Voyatzis-Libert lineage; the vertical-critical mechanism). `3d-extension-of-known-planar-cyclers-AXIS-published-class-2026-06-25`: **axis closed for novelty.** Valuable only as reproduction/gauntlet exemplars. *(The spatial-novelty prefilter exists precisely to route these away cheaply — see §5.)*
- **D2. Planar CR3BP envelope.** Single-family Jacobi continuation from sourced symmetric seeds surfaces only reproductions / same-family points (9 seeds, 0 distinct families). The envelope is mined. Re-open ONLY when a global multi-rev enumeration lane (DA/HOTM, Png′-class) ships — capability-subsumption, not re-run.
- **D3. Ballistic moon-tours in ice-giant / most gas-giant systems.** Uranian + Neptunian regular moons have **disjoint Tisserand contours — unlinkable at any V∞ (4-15 km/s)**; Jovian/Saturnian feasible closures are published resonant tours, fresh ones are physically infeasible (high-V∞ basin, gap 7-21 km/s). Structural negatives. Re-open only with a low-thrust or 3D/inclined releg genome (Region F).
- **D4. ER3BP-only families off Earth-Moon CR3BP roots (e≤0.15).** All three mechanisms (continuation-from-CR3BP, high-e-seed, direct-e>0) exhausted; families persist, no branch-off. The blind-grid discriminator is unreliable. (Region A/C differ: A is *isolated* families with no CR3BP root; C is the *skipped high-e Sun-planet* µ.)
- **D5. Classic Earth-Mars precursor-MGA.** The global DE440 + per-leg-DSM engine reaches geometries the local optimizer could not, but finds **0 literature-fresh insertions** — the classic E-M precursor space is well-studied. Re-open only with Earth-launch-enforced enumeration or a 2030s-synod-specific epoch sweep (an engineering, not discovery, payoff).

---

## 5. The standing discovery PIPELINE (pre-filter runs FIRST)

The C21 lesson (`spatial_novelty_prefilter.py` docstring): the C21 candidate cost a full gauntlet (Floquet + V1 + V2 + JPL full-family sweep) BEFORE a cheap literature check closed it as published-class. **The expensive verification ran first and the cheap decisive check ran last.** Invert that. Every campaign composes the existing capabilities as:

```
0. NOVELTY PRE-FILTER (cheapest, FIRST)
   - genome/spatial_novelty_prefilter.py  (is the planar/known root a published class? → route to reproduction ledger, NOT discovery)
   - check data/empty_regions.jsonl: is this region a method-versioned negative? (capability-subsumption rule — only re-run if a strictly-more-capable method)
   - check catalogue.yaml signature dedup
        │  (kill published-class & re-mined regions before any propagation)
        ▼
1. GENOME (the right model for the region)
   - ER3BP isolated (A) / 3D tracer (B) / ER3BP continuation (C) / ...
   - produces closure-verified candidates with Floquet stability
        ▼
2. GAUNTLET V0→V5 (model-specific ladder)
   - V0 lit-check (search/literature_check.py + known_corpus_3d) — necessary-not-sufficient
   - V1 (code-path agreement) → V2 (multi-lap / real-eph closure) → V3 (dv-band maintenance-ΔV)
   - SILVER holding tier (novel, un-sourced) caps here pending human review
        ▼
3. LITERATURE NOVELTY CONFIRMATION (mandatory gate, §16.5)
   - full publication-record check; "not-found" clears, never certifies
        ▼
4. CATALOGUE (human-promoted only)
   - V0 = sourced; V1+ = independently reproduced; provenance-tagged
```

**Three standing rules baked in:** (1) pre-filter is the FIRST expense, not the last (C21 lesson); (2) the negative registry is consulted before every sweep and only a capability-subsuming method re-opens a region; (3) no row is "novel" until `literature_check` clears it against the published record — not-found is necessary-not-sufficient.

---

## 6. New-input levers (what OPENS currently-closed regions)

Per the validation-ceiling logic, the closed regions are NEW-INPUT-gated, not iteration-gated. The levers, ranked by reach:

1. **A low-thrust / DSM releg genome** — re-opens the entire ice-giant + gas-giant moon-tour space (D3, Region F): the structural emptiness is *ballistic-only*; the fresh-but-infeasible Amalthea-class cycles could be rescued by a powered releg. Highest reach.
2. **A global multi-rev fixed-point enumeration lane (DA / HOTM, arXiv:2509.12671 Png′-class)** — strictly subsumes single-family continuation (D2, Region G); re-opens the EM C≈3.0 band and every CR3BP system. Already has a published validation target (Png′).
3. **Figure digitization of Antoniadou 2018 high-`(e1,e2)` configs** — the precision-limiting input for Region A; cheapest high-value acquisition.
4. **CR3BP Lyapunov/DRO seeds at Sun-Mars / Sun-Mercury µ** — un-blocks Region C (the skipped corner); a generation step, not an acquisition.
5. **3D-inclined releg genome** — re-opens inclined moon-tours (Amalthea is inclined).
6. **The QP V3+_qp and BCR4BP V0-V5 gauntlets** (`#333`, `#305` drafts) — until these ship, QP and BCR4BP finds cannot be promoted; building them is the gate that makes Regions D/E searchable at all.
7. **JPL 3-Body Periodic Orbit Catalog (#116)** + Russell&Ocampo 2006 / McConaghy 2006 full text — turns sparse 3D lit-checks from false-negative into real adjudication; uncertain payoff (likely same parent-only coverage).

---

## 7. Recommended #1 search target + rationale

**RECOMMENDED #1: Region A — the ER3BP isolated-family cycler probe, run as the gated phased plan.**

Rationale, in one line each:
- It is the **only region in the entire map where the DYNAMICAL family is provably unmapped** for a cycler-relevant µ — every other open region refines a known family, characterizes known structure, or re-tests a conditional negative.
- The existence prior is **published and strong** (Antoniadou & Libert 2018 prove isolated high-`e` no-circular-limit families exist for all five MMRs) — so the *capability* win is guaranteed even if the cycler prize fails.
- The capability is **the freshest, least-exercised lane we have** (`er3bp_isolated_seeds.py`) and the gauntlet (ER3BP→3D/heliocentric V0-V2, Antoniadou corpus already anchored) is ready.
- It is **disciplined by design**: a bounded capability gate first, a speculative cycler pass second, a hard falsifiable kill-criterion (death into a distant-perturber PO with no Earth↔Moon leg) — it cannot become an unfalsifiable money-pit.
- It runs **in parallel with the near-zero-cost Region B** (usability of the stable 3D C21/C32 families), so the budget buys two disjoint shots: one high-novelty/speculative-usability (A) and one low-novelty/high-usability (B). Region C is the third, lower-priority shot, run only if A and C-prerequisite seeds are cheap.

**The single most important caveat to carry into review:** Region A's cycler prize is speculative because of the *semantic mismatch* — Antoniadou's MMRs are inner-body-vs-distant-perturber with **no transfer leg between two primaries**, which is the defining feature of a cycler. The bounded capability gate is worth running regardless (it reproduces an unmapped family); the cycler pass must be held to its kill-criterion ruthlessly. If the high-`e` EM member dies into a distant-perturber PO, that is a *clean, publishable negative*, and the axis closes — which is itself a successful outcome under this project's discipline.

---

## 8. Open questions for the user (decisions this draft does NOT make)

1. **Budget split A vs B vs C.** Recommendation: A (capability gate) + B (cheap, parallel) immediately; C after seed-generation cost is scoped.
2. **Does a stable-but-published 3D cycler (Region B) earn a catalogue row?** It is not novel dynamically; it would be a V0 known-attribution admission with a usability characterization — needs the same scope decision as the Jovian CGE rows.
3. **Low-thrust genome (lever 1) vs DA/HOTM lane (lever 2) as the next capability build** — both re-open large dead regions; lever 1 has wider reach, lever 2 has a ready validation target.
4. **Region A digitization precision** — the Antoniadou figure digitization is the accuracy bottleneck; is figure-digitized precision acceptable for a capability-gate reproduction, or is a data request to the authors warranted?
