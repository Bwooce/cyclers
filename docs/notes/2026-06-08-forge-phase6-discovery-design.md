# Forge Phase 6 — forward NOVELTY-DISCOVERY campaign (design note)

**Date:** 2026-06-08
**Scope:** docs-only design pass. No `src/` or `tests/` production file touched.
Deliverable = this design note + the executable first-campaign plan
(`docs/superpowers/plans/2026-06-08-forge-phase6-first-campaign.md`).
**Status of the Forge:** Phases 0–5 are COMPLETE and SHIPPED
(`docs/superpowers/plans/2026-06-03-the-forge-pipeline.md`, "Completion notes —
Phases 4 + 5"). Phase 6 is the forward extension: finding cyclers that are NOT in
the literature, as distinct from everything to date, which has *rediscovered or
validated KNOWN* cyclers.

---

## 0. What is genuinely new in Phase 6 (and what is not)

Everything before Phase 6 has been *rediscovery / validation* — close a cell,
match its signature, and confirm it lands on a sourced catalogue/literature row
(the Aldrin/Russell/Jones/McConaghy anchors). Phase 6 inverts the goal: close a
cell, **prove the signature misses every catalogue row AND the literature**, and
hold it as a candidate discovery.

The machinery that does this *already exists* — Phase 4/5 built
`discover_novel` + the adversarial panel + the human-review queue and ran a first
campaign (the E-M multi-arc Sanchez window: 12 closed families, 10 REJECTED,
2 SILVER). Phase 6 is therefore **NOT new infrastructure**. It is:

1. **A family choice** — point the existing loop at the single most
   under-explored family with the best novelty odds and lowest cost.
2. **A novelty-gate hardening** — the catalogue dedup (`canonical_signature` /
   `signature_distance` / `match`) is wired; the *literature* half of the V5 gate
   (spec §14: "misses catalogue AND literature") is currently a manual human
   step. Phase 6 makes the literature check an explicit, recorded queue step, not
   an implicit one.
3. **A first-class negative-result format** — "region X swept, no ballistic
   closure to V∞ floor Y" must be captured (the Hughes/Jones empty-set lesson),
   never silently dropped. Today the empty-set lessons live only in
   `data/OUTSTANDING.md` prose; Phase 6 gives them a structured report.

This note designs all three against the named modules and recommends the first
campaign. The honest headline (stated up front, governs §7): the classic
heliocentric E-M / VEM spaces are well-mined and the documented sweeps there
returned **empty** of bend-feasible novelty; the best remaining odds are in the
**planet-centric moon-system space**, whose substrate Tier-1 (#76) just shipped
but has never had a *forward novelty sweep* run over it.

---

## 1. The two-part problem: FINDING that closes + PROVING novelty

### 1a. FINDING — construction-first closure (already built)

The finding half is `discover_novel`
(`src/cyclerfinder/data/discover_novel.py`): drive the real-ephemeris ballistic
corrector (`search/correct.py:ballistic_correct` via `search/scan.py:scan_parallel`)
over a topology × epoch grid, keep only *closed* results
(`closed_only=True`), bridge each to a full `Cycler`
(`cycler_from_closure` → `construct_cycler`), and evaluate it. A closure is only
the *entry ticket*; physical admissibility (converged AND bend-feasible AND
within the V∞ cap) is enforced in `evaluate_closure` — a bend-infeasible closure
is auto-`falsified` and routes to REJECTED, never SILVER. This is the
orbit-closure-discipline rule baked in: "a flyby that cannot deliver the required
turn is an impossible bend."

Phase 6 reuses this verbatim. The only finding-side change is the
**topology/centre the loop scans** (the family choice, §5) and wiring the
**incremental-pruning gate** (§2) so the scan does not burn budget on
sequences whose first leg already fails feasibility.

### 1b. PROVING novelty — the dedup gate (catalogue) + literature check

The proving half is two independent checks, both of which must return "miss":

- **Catalogue miss** — `canonical_signature` reduces the closed cycler to the
  spec §16.2 identity object (body set, canonical sequence, sense, `period_k`,
  binned V∞ multiset, binned leg-element multiset, `model_assumption`).
  `match(signature, catalog)` (`data/catalog.py:967`) filters the catalogue to
  the candidate's **`(model_assumption, primary, bodies, k)` bucket** and returns
  `known` (exact hash hit), `probable-match-NEEDS-HUMAN` (weighted-L1
  `signature_distance` below `TAU_NEAR`), or `novel` (no hit). Only `novel` —
  *and no unresolved `probable`* — can become a discovery.
- **Literature miss** — the catalogue is the *transcribed* literature
  (Aldrin/Byrnes, Russell&Ocampo, McConaghy, Niehoff, Jones; seed list spec
  §14 "Seed catalogue with citations"). A signature that misses the catalogue has
  missed *the literature we have ingested*. The residual risk is a published
  cycler we never transcribed. Phase 6 closes this with an explicit
  **literature-review queue step** (§3c): the human-review queue entry carries a
  `literature_checked` field that a reviewer must fill before promotion — the V5
  gate's "documented literature review returned nothing" (spec §613 hard rule),
  made a recorded artefact rather than an implicit assumption.

The firewall, restated from spec §14 / §448 / the trust ladder: **a re-derived
known cycler is NOT a discovery.** This is enforced structurally —
`classify_candidate_verdict` passes `known_id` through to `run_gauntlet` whenever
`match_outcome != "novel"`, so a matched candidate can only reach GOLD via an
*independent source* (which the discovery loop never fabricates), and a novel
candidate has no source by definition → it caps at SILVER pending human review.
GOLD is structurally impossible for the discovery loop to self-assign.

---

## 2. The search pipeline (wired to the named modules)

The pipeline is five stages; the corrector/gauntlet/dedup stages are shipped, the
**enumerate** and **incremental-prune** stages are the Phase 6 wiring:

1. **Enumerate sequences** — pick the under-explored family (§5: planet-centric
   moon systems, now body-agnostic post-#76). Topologies are `TopologySpec`
   objects (sequence + per-leg rev/branch + period + ToF seed + slack leg), the
   same shape `em_multiarc_topologies()` returns. For moon systems the centre is a
   primary (`Ephemeris(model="circular", center="Jupiter")`,
   `mu_central=PRIMARIES["Jupiter"]`).
2. **Tisserand + bend + VILM prune (incremental, Ceriotti-style)** — gate each
   leg by the *previous* leg's feasibility box before the next leg's grid is
   built, per Ceriotti 2010 Ch.3 (incremental pruning = ~90% lossless volume
   collapse; the partial objective is a Bellman lower bound, so pruning never
   discards the global optimum). The per-leg gate is NOT the ΔV — it is the
   geometric/physics feasibility we already compute: `tisserand.linkable` /
   `linkable_3d` contour intersection (centre-aware post-#76), `_max_bend_deg`
   bend feasibility, and — the moon-system-specific piece — the
   **VILM admissible ΔV-floor** (`search/vilm.py:vilm_dv_floor` /
   `vilm_dv_min`): prune a moon-pair leg whose V∞-leveraging ΔV-floor already
   exceeds the budget. This is the cheapest-decisive prune for the moon space
   because the #76 honest-risk finding is precisely that the *no-leveraging*
   corrector lands bend-infeasible at ~10 km/s; VILM is the leveraging physics
   that tells us which legs *can* be bend-feasible.
3. **Scan + correct + continue** — `scan_parallel(grid, closed_only=True)` over
   the surviving (pruned) grid; each closed result flows to the per-candidate
   evaluation.
4. **Per-candidate independent cross-check (the orbit-closure discipline)** —
   the finding-solver's residual is NEVER sufficient. Each survivor is re-closed
   by an *independent* path: Axis-A `crosscheck_code_paths` (in-house Lambert vs
   lamberthub izzo+gooding; forward Kepler re-propagation; the resonance-vs-
   optimiser path is dropped for multi-arc candidates per the cross-fidelity
   rule), and — for any SILVER survivor — the n-body rung
   (REBOUND/IAS15 over the shared DE440 BSP, the SILVER-ARTIFACT template, §4).
5. **Gauntlet → tier** — `run_gauntlet` combines Axes A–D + supersession into
   GOLD / SILVER / BRONZE / REJECTED. V3+ is the credibility floor; only
   SILVER-novel + catalogue/literature miss routes to the human V5 gate.

**Five-line summary (for the report):**
1. Enumerate moon-system topologies (centre = primary, flyby = its moons).
2. Incrementally prune leg-by-leg on Tisserand `linkable` + `_max_bend_deg` +
   the VILM ΔV-floor (Ceriotti box-gate; ~90% lossless collapse).
3. `scan_parallel(closed_only=True)` over the survivors; bridge each closure to a
   `Cycler`.
4. Dedup-gate (`canonical_signature` → `match`) for catalogue miss; independent
   cross-check (Axis-A + n-body rung) for the orbit-closure discipline.
5. Gauntlet → SILVER-novel routes to the human queue with a literature-check
   field; everything else (REJECTED / known / probable) is recorded, not promoted.

---

## 3. The novelty gate wired into the discovery loop

### 3a. How `canonical_signature` / `signature_distance` flag a catalogue miss

`canonical_signature(cycler, model_assumption=..., period_k=...)` produces the
identity object; `match` buckets by `(model_assumption, bodies, k)` (and, post-#76,
`primary` via `signature_bucket_key`) and returns one of three states. The
moon-system bucket is critical here: a **Jovicentric V∞ never compares to a
heliocentric V∞** (the (model_assumption, primary) pre-filter, #76 Phase 6) — so a
moon-tour candidate is matched only against the two Jovian rows + the Saturnian
row, all of which are **family-seed null-numeric records with no sourced
Jovicentric V∞ multiset**. That means a *closed* Jovicentric cycler with a real
V∞ multiset will almost certainly return `novel` (there is nothing populated in
its bucket to match). This is exactly why the moon space has the best novelty odds
(§5) — and exactly why the novelty gate must be paired with a *strong* literature
check (§3c), because "novel" here partly reflects a sparsely-populated bucket, not
necessarily a genuinely unpublished trajectory.

### 3b. The supersession-aware match (R1 delta 3)

`evaluate_closure` extracts the matched row's `superseded_by` chain
(`_superseded_chain`) and surfaces it in the verdict. A candidate matching a row
that carries `superseded_by` (e.g. the UNREALIZED `vem-emeeve-3syn`) must NOT
report a clean "known" against an invalidated premise — the supersession chain
rides the verdict so the human sees the candidate matched a *retracted* row.

### 3c. The literature-check step (the Phase 6 addition)

The V5 gate (spec §613) requires four conditions for a novelty claim:
**V5 AND no exact match AND no unresolved probable-match AND a documented
literature review returned nothing.** The first three are machine-checkable and
wired. The fourth is human and is currently *implicit*. Phase 6 makes it explicit:

- The human-review queue entry (`data/review_queue.jsonl` via
  `data/review_queue.py:append_review_entry`) gains a recorded
  `literature_check` block: `{checked: bool, reviewer, date, sources_searched:
  [...], result: "no-match" | "matches <cite>"}`. (Schema addition in the queue
  payload — the queue is a non-catalogue artefact, so this is additive and does
  not touch catalogue schema.)
- A SILVER-novel candidate is **not eligible for promotion** until
  `literature_check.checked == true`. If the review finds a match, the candidate
  is downgraded to `known-reproduction` and the citation attached (spec §612
  retroactive-correction path) — the candidate becomes a *catalogue ingest*, not
  a discovery.

### 3d. The firewall — re-derived known is NOT a discovery

Restated as the binding invariant for Phase 6: the discovery loop NEVER writes a
catalogue row (`is_catalogue_source()` is `False` by contract) and NEVER
auto-promotes past SILVER. The only paths out of SILVER are (a) human V5
promotion after a clean literature check → a new `first_published` = us, or (b)
human downgrade to `known-reproduction` on a literature hit. Both are human
decisions on the queue; the machine holds at SILVER.

---

## 4. The SILVER-ARTIFACT independent-verification template

The template is the n-body SILVER-rungs result
(`docs/notes/2026-06-06-nbody-silver-rungs.md`) — the pipeline's first honest
false-positive kill. Both held SILVER candidates (E-M-E-E, E∞~9.7 / M∞~12-13)
graded **ARTIFACT**: the conic seed lives only in patched-conic land; an
*independent* integrator (REBOUND/IAS15 over the shared DE440 BSP) diverged on it,
and the divergence — not the ΔV magnitude — is the load-bearing verdict.

The template, binding for every Phase 6 SILVER survivor before it reaches a human:

1. **Aggressive search, ruthless independent verification.** The finding-solver's
   residual is never a V-level. Re-close with an integrator/method that did not
   build the candidate.
2. **All binding constraints in the residual.** V∞, ToF, radius, longitude/phase
   rendezvous with the true ephemeris body, bend feasibility — "it closed!" is the
   danger signal; immediately ask "closed on WHICH constraints, and which did I
   omit?" (the #164 subset-closure self-deception).
3. **Golden target in the SAME MODEL as the test.** A real-eph candidate compares
   to real-eph values, never to coplanar/spec idealisations (the S1L1
   5.65-vs-4.99 episode).
4. **A clean negative (DIVERGE / ARTIFACT / EMPTY-SET / quantified-partial) is a
   SUCCESS.** Never loosen a tolerance to manufacture CONFIRMED.
5. **Hold the writeback until independent confirmation.** Promotion is a human
   decision on the queue, never the finding-solver's say-so. The n-body rung
   records `promoted=False` and `independence: "shared-DE440 cross-check (NOT a V4
   stamp)"` — only the GMAT lane (independent codebase + ephemeris) is V4.

For Phase 6, the rung's body set must be justified by the same evidence rule:
include a perturber iff it moves the rung metric by more than the consumer's
tolerance (the #131 Jupiter sensitivity test: Jupiter moved the metric ~0.95 km/s
> 0.2 km/s ROBUST threshold, so it earned inclusion). For moon-system candidates
the analogous question is which other moons / the Sun must be in the rung.

---

## 5. Recommended FIRST campaign — the planet-centric moon-system VILM space

### Why this family (the single best novelty odds + lowest cost)

The candidate families and their state:

| Family | Substrate | Documented result | Novelty odds |
|---|---|---|---|
| Classic heliocentric E-M | shipped | well-mined; Phase 4 first run found 2 SILVER (graded ARTIFACT by n-body) | LOW — mined |
| VEM higher-beat | shipped | #110 dense scan: 0 bend-feasible (floors 17.9/18.5 vs sourced 2.4-7.0) | NIL — empirically refuted |
| Inclined / 3D VEM | shipped (M-3D) | #120: inclination moved floors -0.37/-0.10 km/s, still 0 bend-feasible | NIL — refuted |
| Vector-residual / B-plane VEM | shipped | #122: in-residual hinge collapses to ~0 closures; the one survivor sits ~21 km/s | NIL — refuted |
| **Planet-centric moon systems** | **shipped #76, NEVER swept forward** | I-E-G CLOSES about Jupiter; bend-infeasible in the *no-leveraging* model; **VILM is the missing leveraging physics** | **BEST — substrate built, search never run** |

The moon-system space is the clear pick:

- **The substrate just shipped and has never had a forward novelty sweep.** #76
  delivered the `SATELLITES`/`PRIMARIES` registry, the planet-centred ephemeris,
  the centre-agnostic corrector (`mu_central`), centre-aware Tisserand, the VILM
  module, and the `(model_assumption, primary)` dedup bucket. Every piece the
  pipeline needs exists; the loop has simply never been pointed at it.
- **The bucket is near-empty, so closures will read `novel`.** The only Jovian
  catalogue rows are family-seed null-numeric records (no sourced V∞ multiset) —
  a closed Jovicentric cycler with a real multiset has almost nothing to match
  against (§3a).
- **The cost is low and the first sweep is decisive.** The #76 honest-risk
  finding already tells us the no-leveraging corrector lands bend-infeasible at
  ~10 km/s. The cheapest-decisive first sweep is therefore: **VILM-route the
  topology FIRST** (use `vilm_dv_floor` to keep only moon-pair legs whose
  leveraging ΔV-floor is within budget), then run the corrector only on the
  VILM-admissible legs. This directly tests the one open question — does the VILM
  leveraging layer + Laplace-resonance phasing produce a *bend-feasible* Jovian
  tour the no-leveraging #76 closure could not? — at the cost of one moderate
  parallel scan, not a blind dense grid.

### The cheapest-decisive first sweep

Galilean Jovian system, the empirically-closing **Io-Europa-Ganymede** chain
(#76 Phase 3 closes it; the open question is bend feasibility), VILM-gated:

1. Enumerate I-E-G(-I) topologies at Laplace-resonance ToF seeds (Io 1.769 d,
   Europa 3.551 d, Ganymede 7.155 d synodic spacings).
2. Prune each moon-pair leg on `vilm_dv_floor` ≤ budget AND `linkable`
   (Jovicentric) AND `_max_bend_deg` feasibility — incrementally, leg by leg.
3. `scan_parallel(closed_only=True)` over the survivors, `center="Jupiter"`,
   `mu_central=PRIMARIES["Jupiter"]`.
4. Dedup-gate each closure (almost all `novel` by §3a); independent cross-check
   each SILVER (n-body rung about Jupiter, body set justified by sensitivity).
5. If empty: emit the empty-region report (§6) — "Jovian I-E-G VILM space swept,
   no bend-feasible ballistic closure below V∞ floor Y" — and that is a publishable
   negative (the #76 honest-risk finding made rigorous and bounded).

### Honest expected yield

**Most likely: mostly negatives, possibly an empty set — and that is a success.**
The same pattern that produced the VEM empty-sets (#110/#120/#122) and the E-M
SILVER-ARTIFACTs (#131) is the base rate here. Three honest outcomes, in
decreasing likelihood:

- **EMPTY (most likely):** VILM gating prunes the bend-infeasible legs but no
  *closed* chain is simultaneously bend-feasible AND within budget — i.e. the #76
  honest-risk finding generalises (the no-leveraging family is the only one that
  closes). This is the Hughes/Jones empty-set lesson and must be captured as a
  first-class result (§6), not buried in OUTSTANDING prose.
- **SILVER-then-ARTIFACT:** a chain closes and reads novel, but the n-body rung
  diverges (the conic seed is patched-conic-only, exactly the #131 outcome). The
  candidate is recorded REJECTED-style for the human; no promotion.
- **SILVER-survives (least likely, the prize):** a chain closes, reads novel,
  survives the n-body rung AND the adversarial panel, and passes a clean
  literature check → the project's first genuine machine-confirmed novelty
  candidate, held at SILVER pending human V5. We should *expect not to get one on
  the first sweep* and say so plainly.

The value of the campaign is NOT a guaranteed discovery — it is (a) the first
forward sweep of a freshly-built space, (b) a rigorous bounded negative if it is
empty, and (c) proof the Phase 6 wiring (VILM prune + dedup + literature gate +
empty-region report) works end to end on a non-heliocentric centre.

---

## 6. "Region proven empty" as a first-class result

The empty-set lesson (Hughes/Jones: a thorough sweep that finds nothing is a
publishable scientific result, not a failed run) must be captured structurally.
Today the empty-set findings (#110/#120/#122) live as prose in
`data/OUTSTANDING.md`. Phase 6 gives them a machine-readable report.

### The empty-region report format

A JSON artefact (`data/empty_regions.jsonl`, one record per swept region),
emitted by the campaign orchestrator whenever a region yields zero promotable
candidates. Each record carries enough to make the negative *reproducible and
bounded*:

```json
{
  "region_id": "jovian-IEG-vilm-2026-06-08",
  "family": "planet-centric moon system (Jupiter)",
  "centre": "Jupiter",
  "topologies": [{"sequence": ["Io","Europa","Ganymede","Io"],
                  "per_leg_revs": [0,0,0], "period_k": 1}],
  "method_capability": {"genome": "single-ellipse free-return",
                        "corrector": "ballistic_correct (no-leveraging)",
                        "capability_tags": ["ballistic", "patched-conic",
                                            "single-arc", "coplanar"],
                        "git_sha": "..."},
  "search_extent": {"n_epochs": 256, "span_days": 64.0,
                    "n_topologies": 11, "points_total": 2816,
                    "ephem_model": "circular", "center": "Jupiter"},
  "prune_gates": ["vilm_dv_floor<=budget", "linkable(Jovicentric)",
                  "max_bend_deg feasibility"],
  "result": {"closed": 597, "distinct_families": 312,
             "bend_feasible": 0,
             "best_max_vinf_kms": 10.4,
             "vinf_floor_target_kms": 6.0,
             "gap_kms": 4.4},
  "verdict": "EMPTY — no bend-feasible ballistic closure below the V_inf floor",
  "interpretation": "no-leveraging corrector closes a higher-V_inf family; VILM "
                    "gating did not surface a bend-feasible Laplace-resonant tour",
  "source_anchors": "none populated in (circular-coplanar, Jupiter) bucket; "
                    "Russell-Strange/Hernandez rows are family-seed null-numeric",
  "run": {"date": "2026-06-08", "host": "...", "cores": 16,
          "git_sha": "...", "wall_s": ...}
}
```

The required fields (the bar for a negative to count): **search extent**
(how much space was actually covered — epochs × topologies × points), **prune
gates** (what was pruned and why, so the negative is not an artefact of
over-pruning), **best achieved vs target** (the V∞ floor gap — the quantified
near-miss), **interpretation** (the sharpened hypothesis, the #110-style "a
denser-scan-still-fails result is real science"), and — load-bearing for the
re-sweep gate (§6a/§6b) — the **method-capability descriptor** (what the method
could *reach*, not just which SHA produced it). A negative is only first-class
if it is *bounded* (states what was NOT covered), *reproducible* (carries the
grid + git SHA), and *method-versioned* (carries the capability envelope). An
empty region with no search-extent record is a silently-dropped negative; an
empty region with no method-capability descriptor is an *unconditional* "empty"
claim — and "empty" is never unconditional (§6a). Both are exactly what this
format forbids.

The mirror of the human-review queue: SILVER survivors go to
`data/review_queue.jsonl`; empty regions go to `data/empty_regions.jsonl`. Both
are non-catalogue artefacts; neither auto-promotes anything.

### 6a. Method-capability versioning — "empty" is conditional on the method

A bare git SHA does NOT express what a method can *reach*: two SHAs of the
single-ellipse genome are equivalent in capability, while the multi-arc genome at
*any* SHA reaches strictly more of the trajectory space. So every
`empty_regions.jsonl` record carries an explicit **method-capability descriptor**
alongside the run/SHA fields — the genome/corrector plus an ordered (partially-
ordered) **capability tag set** describing the envelope it can represent:

- `single-ellipse free-return` — one ballistic arc, no DSM, patched-conic.
- `one-DSM-per-leg` — adds one deep-space manoeuvre DOF per leg.
- `two-arc free-return chain` (multi-arc, #163) — two generic-return arcs;
  reaches geometries a single ellipse cannot represent.
- `n-body shooter` — full-force integrated, supersedes patched-conic seeds.
- `low-thrust` — continuous control, supersedes ballistic/impulsive.

The descriptor is the load-bearing field: it is what lets a *later* campaign
decide whether a recorded "empty" still binds. The crux, restated as the binding
invariant: **"region X is empty" always means "empty as far as method M could
reach"** — never an absolute claim (the same sourced-floor discipline: the floor
is OUR computation, recorded as method-conditional evidence, never an absolute).

### 6b. The re-sweep gate — capability-subsumption, not region-match

Before a campaign sweeps a region it queries `empty_regions.jsonl` and applies the
re-sweep gate: **SKIP the region ONLY IF a prior sweep's method-capability
SUBSUMES the proposed method's.** Region-match alone is NOT the skip criterion — a
prior empty result over the identical box does *not* license a skip if the new
method can reach ground the old one could not.

**The subsumption partial order (`⊐` = "strictly more capable than").** Methods
form a *partial* order — not a total one — over capability:

- multi-arc `⊐` single-ellipse (more arcs reach more geometry; #163 ⊐ #137),
- n-body `⊐` patched-conic (full force supersedes the conic seed),
- powered (DSM / low-thrust) `⊐` ballistic (added control DOF),
- one-DSM-per-leg `⊐` single-ellipse free-return (added per-leg DOF).

`A` **subsumes** `B` operationally iff every capability tag of `B` is reached by
`A` — i.e. `B`'s envelope is contained in `A`'s under the partial order. Then:

- **Proposed method ⊑ prior method (prior subsumes proposed): SKIP.** A weaker or
  equal method re-running a region a stronger method already emptied learns
  nothing new — that is the compute-saving half of the gate.
- **Proposed method ⊐ prior method (strictly more capable): RE-SWEEP.** A new
  genome, an added DOF, higher fidelity, or a different regime JUSTIFIES
  re-sweeping *exactly* the regions the old method could not represent — that is
  the option-preserving half.
- **Incomparable methods (neither subsumes the other): RE-SWEEP.** E.g. a
  broken-plane method vs a coplanar multi-arc reach *different* ground; neither
  envelope contains the other, so the old "empty" does not bind the new method.
  Incomparable never skips.

**Rationale (why both halves matter).** "Empty" is always conditional on the
method: the multi-arc genome (#163) reopened 6.44Gg3/S1L1 ground the single-
ellipse genome (#137) had closed as "empty". A SHA-only or region-match-only gate
would have *permanently* foreclosed that re-discovery. The capability-subsumption
gate is precisely what prevents BOTH re-burning compute on proven-empty regions
(weaker/equal method ⇒ skip) AND losing the option to revisit when a new
capability arrives (more-capable/incomparable method ⇒ re-sweep). It is the
direct analogue, on the negative-result side, of the catalogue's
supersession-aware match (§3b) on the positive side.

---

## 7. Self-review

- **Builds on Phases 0–5, not new infra.** The finding loop, adversarial panel,
  human queue, dedup, and n-body rung are all shipped; Phase 6 is a family choice
  + two wiring additions (incremental VILM prune; literature-check field) + one
  new artefact (the empty-region report, method-versioned with a
  capability-subsumption re-sweep gate, §6a/§6b).
- **Golden-clean.** No computed value is asserted as a golden EXPECTED. The dedup
  is a classification; the VILM gates are sourced (Endgame Part-1 Tables); a novel
  candidate caps at SILVER by construction; GOLD is structurally unreachable by
  the loop.
- **Honest about odds.** The recommended family is chosen *because* its substrate
  is new and its bucket near-empty — but that same near-emptiness is flagged as a
  reason the literature check must be strong, and the expected yield is stated as
  "mostly negatives, possibly empty, and that is a success."
- **Risk.** The biggest risk is a false "novel" from a sparsely-populated bucket
  (§3a). Mitigation: the explicit literature-check gate (§3c) + the n-body
  ARTIFACT rung (§4) — a candidate must survive BOTH an independent integrator and
  a documented literature review before a human can promote it. The second risk is
  that the moon space is *also* empty of bend-feasible novelty (the #76 honest-risk
  finding generalises); mitigation is the empty-region report (§6), which makes
  that outcome a rigorous bounded negative rather than a dead end.
</content>
</invoke>
