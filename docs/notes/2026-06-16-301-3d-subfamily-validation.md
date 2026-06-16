# 3D sub-family validation — lit-check + ML flagger (#301 / #291 Phase 4)

Date: 2026-06-16
Issue: #301 (Phase 4 of #291); follow-up to #299 (commit f8c3fde)
Inputs:
- `data/family_296_3d_subfamilies_299.jsonl` (#299 Phase 3 output, 4 accepted
  Neimark-Sacker sub-families: k ∈ {3, 4, 5, 6}, 145 members total)
- `src/cyclerfinder/search/literature_check.py` (minimum-necessary edit: added
  optional `period_band_tu` field to `CorpusAnchor` + `CandidateSignature`;
  Antoniadou-Voyatzis 2018 anchor declares `period_band_tu=(0.0, 15.0)`)
- `src/cyclerfinder/ml/falsepos_flagger.py` (#256 / #275 trained on the labeled
  corpus seed; AUC(train)=0.986, AUC(LOO)=0.983)
Outputs:
- `data/sub_families_301.jsonl` — 1 header + 145 member rows
- This doc — verdict + #302 recommendation

## Discipline framing

- NO catalogue writeback. Admission is gated by V0-V5 gauntlet. The current
  gauntlet is planar-only and cannot run on 3D periodic orbits. Catalogue
  admission is impossible until #302 adapts the gauntlet.
- The independent-closure residuals come from #299's JSONL (already verified
  by Phase 2's CR3BP corrector + Radau independent-closure gate at 1e-6).
  Per `feedback_orbit_closure_discipline` we did NOT re-integrate; we record
  the existing residuals.
- "Honest negatives are correct outcomes" — and that is the verdict here.

## Per-sub-family verdict

| k | n_members | T_TU range | T_days range | lit_status | anchor | n_lit_fresh | n_ml_fresh | n_novelty_claimable | doubly-hyperbolic |
|---|-----------|------------|--------------|------------|--------|-------------|------------|---------------------|-------------------|
| 3 | 22 | [20.28, 20.39] | [88.19, 88.65] | published (0.95) | Earth-Moon CR3BP family network | 0 | 22 | **0** | 22 / 22 |
| 4 | 41 | [35.84, 36.03] | [155.86, 156.66] | published (0.95) | Earth-Moon CR3BP family network | 0 | 41 | **0** | 41 / 41 |
| 5 | 41 | [29.93, 30.13] | [130.15, 131.01] | published (0.95) | Earth-Moon CR3BP family network | 0 | 41 | **0** | 0 / 41 |
| 6 | 41 | [43.49, 43.61] | [189.13, 189.62] | published (0.95) | Earth-Moon CR3BP family network | 0 | 41 | **0** | 0 / 41 |

All 145 members are classified as **likely-rediscovery** (lit_check status =
`published`, confidence 0.95) of the Earth-Moon CR3BP family-network corpus
(Braik-Ross 2026 / Roberts-Tsoukkas-Ross 2026 / Kumar-Rosengren-Ross 2026 —
the matcher picks the first overlapping anchor with score ≥ 0.70; manual
inspection shows multiple anchors score 0.95 simultaneously).

The Antoniadou-Voyatzis 2018 anchor was CORRECTLY EXCLUDED for all four
sub-families by the new `period_band_tu` filter — the AV catalogue covers
low-integer p:q spatial resonant orbits at T_TU ≲ 15; our sub-families sit
at T_TU 20-44, period-multiplied by the k ∈ {3, 4, 5, 6} Neimark-Sacker
bifurcations. That filter does what the task asked of it.

But the broader Earth-Moon CR3BP corpus (Braik-Ross orbital networks,
Roberts-Tsoukkas-Ross stable cyclers, Kumar et al. resonant transport)
all carry "family network" / "period-multiplying" / "heteroclinic pathway"
language without period-scope restrictions. At the matcher's structural
granularity our sub-families overlap these published frameworks: the
literature CAN claim to cover Neimark-Sacker-born period-multiplied
derivatives of a base family as part of the broader cislunar family network.

This is the honest call. **Not a novelty result.**

## k=3 and k=4 doubly-hyperbolic structural signature

Both the k=3 and k=4 sub-families show the **doubly-hyperbolic 3D
heteroclinic web signature** in 100% of members:

- k=3 first member: Floquet moduli sorted = [261.3, 30.2, ~1, ~1, 0.033, 0.0038]
  — two distinct reciprocal pairs strictly off the unit circle plus the
  trivial pair at 1.
- k=4 first member: Floquet moduli sorted = [248.3, 85.5, ~1, ~1, 0.012, 0.0040]
  — two distinct reciprocal pairs strictly off the unit circle plus the
  trivial pair at 1.

This is the signature that a non-trivial 3D heteroclinic web is anchored on
these orbits: TWO independent unstable / stable manifold pairs at every
member. (The task description named only k=4 as doubly-unstable; the data
shows k=3 carries the same structural signature.)

k=5 and k=6 are singly-hyperbolic + center: one reciprocal pair off the unit
circle, the remaining four multipliers on the unit circle. Structurally less
rich; standard 3D Lyapunov-vertical-style continuations of the parent.

## ML false-positive flagger results

The flagger was trained on the labeled corpus seed (38 samples; 19 FP / 19 TR;
AUC(train) = 0.986, AUC(LOO) = 0.983 — well above the 0.75
reproduce-before-trust gate).

All 145 sub-family members scored `p_fp < 0.5` — the flagger does not suspect
false positives. But this is a **WEAK signal in this context**: the flagger's
feature set was hand-crafted for heliocentric mission-design SILVERs (V_inf
per encounter, bend feasibility, encounter resonance ratios, epoch-artifact
SHA checks). For a 3D CR3BP periodic orbit, only `max_residual_kms`,
`period_days`, `model_assumption`, `closure_method_version`, and
`closure_date` map naturally; the V_inf / bend / encounter-resonance features
fall back to NaN and the flagger imputes with training medians.

The flagger is correct that no past-bug signature fires (closures are clean,
post-fix dates, CR3BP model declared) — but the score does not constitute
strong evidence of a real cycler in the 3D-periodic-orbit class. A purpose-
built 3D-periodic-orbit flagger (Floquet-spectrum sanity, doubly-hyperbolic
gating, energy-band consistency vs the base family) would be more
informative; not built here, deferred to #302+.

## What this means for novelty

- **0 / 145 members are novelty-claimable** under the current corpus +
  matcher.
- The four sub-families are structurally consistent with the Earth-Moon
  CR3BP family-network framework already published.
- The Neimark-Sacker birth mechanism is itself a generic CR3BP feature
  (period-multiplying bifurcation off a Lyapunov family); no individual
  k-row constitutes a discovery beyond what Braik-Ross / RTR / Kumar et
  al. framework-style papers have covered.
- A clean novelty claim from this lineage would require either (a) a
  specific feature of one sub-family that the published frameworks
  CANNOT cover (none surfaces here), or (b) a tighter literature scan
  against the AV/RTR/Braik-Ross PDFs themselves to confirm no specific
  catalogued row lands at our (T, C, state) point — that's a future
  deeper-search task, not the gating concern here.

## Recommendation for #302 (3D V0-V5 gauntlet adaptation)

Even though no sub-family member is novelty-claimable, the **k=4
sub-family** is the highest-priority structural target for the future
3D-gauntlet trial in Phase 5, on these grounds:

1. **3D heteroclinic web signature is strongest at k=4** (Floquet ratio
   248.3 / 85.5 ≈ 2.9 — two distinct unstable directions with similar
   strength, the canonical doubly-unstable feature).
2. **k=4 sits at a NICE period regime**: T_days ≈ 156, which is a
   Sun-frame ~5.1-month quasi-period — useful for any cycler-as-staging
   application.
3. **k=4 is the densest sub-family** (41 members; same size as k=5/6 but
   with the doubly-unstable signature).

### Concrete IC for the priority candidate

Selected as the lowest-`p_fp`, doubly-hyperbolic, closure-OK member of k=4
(JSONL row member_step_index = -19, arc_length = -0.095):

```text
state_nd = [-0.7793493004626235, 0.0, -0.2592654091352357,
             0.0, -0.35160377428751494, 0.0]
T_TU     = 35.8476
T_days   = 155.8787
jacobi   = 2.9370348
parent   = family_296_3d_em_11 (Antoniadou-Voyatzis 2018 spatial CR3BP
           branch; #299 Phase 2 step_index ~ 112 territory)
floquet  = [-251.348, 85.770, ~1, ~1, 0.01166, -0.003979]  (doubly-unstable)
closure  = 9.45e-11  (Radau independent gate, < 1e-9)
ml p_fp  = 0.374  (NaN-imputed for non-mapped features)
```

The full priority block is in
`data/sub_families_301.jsonl` → header.highest_priority_for_302.

### What #302 needs to land before this IC can be admitted to the catalogue

The current `cyclerfinder.gauntlet.*` V0-V5 layers assume the candidate
carries:
- A heliocentric / patched-conic trajectory with body encounters,
- Per-encounter V_inf values, bend feasibility, and encounter epochs,
- Resonance / multi-rev labels in the SnLm sense,
- Closure residuals expressible in km.

None of these apply natively to a 3D CR3BP periodic orbit. The 3D-gauntlet
adaptation needs to (at minimum):
- V0: structural sanity for a 3D periodic orbit — Floquet sextuple sanity
  (trivial pair at 1, the remaining four reciprocal-pair-structured),
  Jacobi constant within the family's range, closure residual ≤ 1e-9 in
  nondim.
- V1: independent re-integration in a DIFFERENT 3D integrator (not the
  Phase 1 Radau used to build the family) — same closure gate, same
  Floquet spectrum within tolerance.
- V2: source-cross-check — confirm the orbit IS in the AV-2018 / Braik-Ross
  / RTR catalogues at the family-row level (a real PDF dig, not just the
  matcher).
- V3+: scientific contribution gate — does this 3D orbit add anything
  beyond the published framework? At present the answer for k=4 is "no
  obvious increment beyond what the family-network framework covers,"
  which is exactly why no #301 member is novelty-claimable. V3 is the
  right place to enforce that.

## Files

- Edited (minimum-necessary):
  `src/cyclerfinder/search/literature_check.py` — added optional
  `period_band_tu` to `CorpusAnchor` + `CandidateSignature`; populated AV-
  2018 anchor's band; per-anchor disjoint-band filter in
  `_candidate_anchors`. Backward-compatible (None preserves old behaviour).
  All 42 existing tests still pass.
- New: `scripts/run_301_subfamily_validation.py` — runner.
- New: `data/sub_families_301.jsonl` — 1 header + 145 member rows.
- New: this doc.
- Read-only (per discipline): `data/family_296_3d_subfamilies_299.jsonl`,
  `src/cyclerfinder/core/cr3bp_general_periodic_3d.py`,
  `src/cyclerfinder/core/cr3bp_3d_family_tracer.py`.

## Discipline / negative-results registry

This row should land in the persistent negative-results registry
(`feedback_negative_results_registry`) once the registry exists, keyed as:

```text
method_version: literature_check.py + 3D-Earth-Moon-CR3BP-corpus + AV-2018-period-banded
input_window:   T_TU ∈ [20, 44], k ∈ {3, 4, 5, 6}, n_members = 145
verdict:        0 / 145 novelty-claimable; all flagged as likely-rediscovery
                of the Earth-Moon CR3BP family-network corpus.
re-sweep when:  a richer corpus or a per-family-row matcher subsumes this
                (e.g. AV-2018 family-row-level scan, Braik-Ross node-set scan).
```
