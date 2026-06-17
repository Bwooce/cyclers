# 2026-06-17 — #344 Phase 2 Stage A verdict — Saturn Titan-Rhea-Titan (1,1)

## What Stage A is

The 10-gate gauntlet for the Saturn Titan-Rhea-Titan (1, 1) SILVER
(residual 0.0102 km/s at the ps=96 phase grid, identified in #344
Phase 1 Part A.2) became lit-fresh-candidate after #346 tightened the
Davis-Phillips-McCarthy 2018 anchor's `body_set` from
`{Titan, Enceladus, Rhea, Dione}` to `{Titan, Enceladus}` (commit
`dabf4a6`).

Stage A is the EXISTENCE-CONFIRMATION pre-flight: four sub-gates that
must all clear before the V1+ gauntlet (3D corrector, moontour,
REBOUND, GMAT — Stages B-E) is fired. This is the same shape as the
#327 verification ladder (Umbriel-Oberon SILVER → Stage B → catalogue
admission at #339).

If all four sub-gates clear → proceed to Stage B (V1 3D corrector).
If any sub-gate fails → HALT and report the blocker honestly.

## Approach

`scripts/verify_344_titan_rhea_silver.py` wraps the production
primitives (no modification of `discovery_campaign.py`,
`literature_check.py`, `physical_sanity.py`, or
`falsepos_flagger.py`) and runs the four Stage A sub-gates on the
SILVER's stored ps=96 IC. The closure machinery (`_close_one` +
`_sweep_one_cycle`) is the same code path as `scan_344_…_finer.py`
Part A.2, which itself is the verbatim port from
`scan_320_epoch_aware_moon_systems.py`. Output:
`data/silver_344_verified.jsonl` (7 rows).

## Findings

### (1) IC verification

**Sub-gate target:** residual < 0.013 km/s AND stored point reproduces
to 1e-9 km/s.

The stored ps=96 record in `data/scan_344_saturn_robustness.jsonl`:

| Field            | Stored value                                    |
|------------------|-------------------------------------------------|
| residual         | 0.010188096573990224 km/s                       |
| phase0           | 273.749… deg                                    |
| rel_offset       | 288.75 deg                                      |
| tof_scale        | 2.0                                             |
| V_inf (per enc)  | (1.7375, 1.6463, 1.7273) km/s                   |
| ToF per leg      | 16.977266 days                                  |

Point reproduction at the stored phasing returns residual
**0.010188096573990224 km/s** (matches stored to machine precision).
Full ps=96 sweep (96 × 96 grid × 4 tof_scales = 36 864 cells)
reproduces the basin-floor best at the same cell (phase0 = 273.75°,
rel_off = 288.75°, tof_scale = 2.0, residual 0.01018810 km/s,
V_inf [1.7375, 1.6463, 1.7273] km/s) — confirming both the stored row
AND the basin floor are stable under the same `_sweep_one_cycle`
convention used in #344 Phase 1 A.2.

**Verdict:** PASS. The IC is byte-reproducible from the stored record
and the basin floor is in the same cell.

### (2) Lit-fresh anchor count post-#346

**Sub-gate target:** anchor count = 0 (strict gate per task brief);
> 0 → HALT.

`_candidate_anchors` over the post-#346 `KNOWN_CORPUS`
(`literature_check.py` at HEAD `1f93128` / corpus state `dabf4a6`)
applied to `CandidateSignature(primary='Saturn',
sequence=('Titan','Rhea','Titan'), period_k=2,
vinf_per_encounter_kms=(1.7375, 1.6463, 1.7273), n_rev=(1,1))` returns
**1 anchor**:

| Anchor                                            | `body_set`                                       | Matches because                  |
|---------------------------------------------------|--------------------------------------------------|----------------------------------|
| Cassini-Huygens Saturn-Titan satellite tour design | `{Titan, Enceladus, Rhea, Dione, Iapetus}`       | `seq_set <= anchor.body_set` ✓   |

The pre-#346 reported count was 2 (Davis-Phillips-McCarthy 2018 +
Cassini-Huygens). After #346 tightened Davis's `body_set` to
`{Titan, Enceladus}`, the `seq_set <= anchor.body_set` subset test
fails for Davis (Rhea ∉ {Titan, Enceladus}), so Davis no longer
matches. Δ count: **2 → 1**.

The remaining anchor (Cassini-Huygens, line 887 of
`literature_check.py`) keeps `body_set = {Titan, Enceladus, Rhea,
Dione, Iapetus}` because the published Cassini tour did include Rhea
flybys (4 Rhea flybys over the mission). Body-set subset is the
discipline's necessary-not-sufficient gate
(`feedback_literature_novelty_check_baseline`); the anchor fires
on the structural-fingerprint test even though Cassini-Huygens is
not specifically a (1, 1) Titan-Rhea-Titan repeating ballistic cycle
(it is an MGA tour with Titan as the pump moon and Rhea as a one-off
target).

**Verdict:** FAIL under the Stage A strict zero-anchor gate.

### (3) Physical-sanity gate (#324)

**Sub-gate target:** max bend > 5° at every encounter.

`candidate_passes_physical_gate(sequence, V_inf_tuple,
min_useful_bend_deg=5.0)` at the verified IC:

| Encounter   | V_inf (km/s) | Max bend (deg) | Useful? |
|-------------|-------------:|---------------:|:-------:|
| Titan #1    |        1.7375 |          49.91 |   yes   |
| Rhea        |        1.6463 |           7.07 |   yes   |
| Titan #2    |        1.7273 |          50.27 |   yes   |

All three encounters clear the 5° floor. Rhea is the tightest
(7.07°, 1.4× the floor — Rhea's GM 153.94 km³/s² is two orders of
magnitude lower than Titan's 8978.14 km³/s² and the V_inf is only
marginally lower, so the Rhea bend is naturally the gate-critical
encounter). The Titan flybys clear the floor 10× over.

The numbers match the #344 Phase 1 doc qualitative summary (Titan
~50°, Rhea ~6.8° at V_inf ≈ 1.68 km/s); the small shift in Rhea bend
(6.8° → 7.07°) is because the verified ps=96 V_inf at Rhea is
1.6463 km/s (deeper in the basin), and bend angle is monotonic in
1/V_inf² for a fixed body.

**Verdict:** PASS.

### (4) ML flagger (#256)

**Sub-gate target:** `p_fp ≤ 0.75` (production `P_FP_SILVER_MAX`
from #274) → classification "real".

`FalsePosFlagger` (logistic regression over the labeled-corpus seed,
trained from `build_training_set()`) on the production feature
schema (mirrored from `saturn_uranus_campaign.score_candidate`):

```
p_fp = 0.604123
threshold (P_FP_SILVER_MAX) = 0.75
classification = real
```

`0.604 < 0.75`. This is qualitatively similar to the #327 Umbriel
SILVER's `p_fp = 0.5918` — both are "real but borderline" by the
flagger's seed corpus (the corpus contains few moon-system SILVER
labels, so the prior is conservative).

**Verdict:** PASS.

## Sub-gate summary

| Sub-gate              | Target          | Measured                | Verdict |
|-----------------------|-----------------|-------------------------|:-------:|
| (1) IC residual       | < 0.013 km/s    | 0.01018810 km/s         | PASS    |
| (2) Lit-fresh anchors | == 0            | **1 (Cassini-Huygens)** | **FAIL**|
| (3) Physical sanity   | > 5° all encs   | 49.9° / 7.07° / 50.3°   | PASS    |
| (4) ML flagger p_fp   | ≤ 0.75 ("real") | 0.604 ("real")          | PASS    |

## Verdict

**HALT — Stage A sub-gate (2) fails.**

The candidate is NOT lit-fresh under the strict zero-anchor gate
post-#346. The Cassini-Huygens Saturn-Titan satellite-tour-design
anchor still matches because its `body_set` includes Rhea, and the
`seq_set <= anchor.body_set` subset test fires on `{Titan, Rhea}`.

Per the task brief disposition: this is a #346-followon, not for
this Stage A task to fix. The Cassini-Huygens anchor's `body_set`
is conservatively broad (the entire body set the Cassini mission
visited), but the published archetype is a Titan-pumped MGA tour,
not a Titan-Rhea-Titan (1, 1) ballistic repeating cycler. A future
ticket should deep-read Strange et al. JGCD 2010-2017 and Goodson
et al. JGCD 2008 to determine whether Rhea should remain in the
anchor's `body_set` for the structural-fingerprint test, or whether
the anchor should be split (e.g. one anchor for the Titan-pump MGA
archetype with body_set ⊇ {Titan} only, a separate anchor for the
once-off Rhea encounters during the Cassini mission, etc.).

## Path forward

* **NO Stage B (V1 3D corrector).** The lit-fresh prerequisite is
  unmet; firing the V1+ gauntlet would force admission of a
  candidate the matcher flags as anchored, contradicting the
  `feedback_literature_novelty_check_baseline` discipline.
* **NO catalogue writeback.** Stage A is pre-gauntlet; admission is
  Stage E (with V4 evidence proven) only.
* **#346-followon ticket** to deep-read the Cassini Saturn tour
  design literature and decide whether the Cassini-Huygens anchor's
  `body_set` ⊇ Rhea is sound or whether it should be tightened. If
  tightened to `body_set = {Titan}` (Titan-pump archetype only), the
  Stage A (2) sub-gate flips to PASS (anchor count → 0) and Stage B
  becomes the unblocking step.
* **Empirical results in this Stage A remain valid as evidence**:
  the IC is byte-reproducible, the basin floor is at the same cell
  on the ps=96 grid, the physical-sanity gate clears with margin
  (the Rhea bend at 7.07° is the gate-critical encounter for future
  fidelity steps), and the ML flagger classifies the row as "real".
  Sub-gates (1), (3), and (4) are stored in `silver_344_verified.jsonl`
  for re-use if (2) reopens later.

## Anti-claim ledger

What this Stage A verification does NOT establish:

* That the candidate is novel — sub-gate (2) is the discipline's
  gate-keeper for the novelty claim, and it fails. The 0.0102 km/s
  closure remains a body-set-anchored finding under
  Cassini-Huygens.
* That the candidate cannot survive the V1+ gauntlet — Stages B-E
  were not run. The math suggests the candidate has the structural
  characteristics of a real SILVER (deep basin, clean physical
  sanity, ML "real"), but the lit-fresh blocker is upstream of those
  steps in the discipline.
* That sub-gate (2)'s failure means the underlying physical
  trajectory is rediscovery. The anchor matches on body_set only
  (necessary-not-sufficient); the topology / V_inf / period match
  to the actual Cassini-Huygens tour is NOT established. The honest
  reading is "lit-fresh prerequisite unmet" not "rediscovery".

## File index

* Script: `scripts/verify_344_titan_rhea_silver.py`
* Output: `data/silver_344_verified.jsonl`
* Phase 1 doc: `docs/notes/2026-06-17-344-saturn-titan-rhea-extended-sweep.md`
* Source robustness JSONL:
  `data/scan_344_saturn_robustness.jsonl` (ps=96 row)
* Reference template:
  `docs/notes/2026-06-16-327-umbriel-silver-verification.md`
* Corpus state: `src/cyclerfinder/search/literature_check.py`
  at commit `dabf4a6` (Davis-Phillips-McCarthy 2018 body_set
  `{Titan, Enceladus}`; Cassini-Huygens body_set
  `{Titan, Enceladus, Rhea, Dione, Iapetus}`).
