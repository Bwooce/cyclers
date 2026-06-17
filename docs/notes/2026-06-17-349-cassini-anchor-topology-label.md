# #349 — Cassini-Huygens anchor scope investigation + `topology_label` schema

**Date:** 2026-06-17 AET
**Task:** #349 — resolve the #344 Phase 2 Stage A HALT verdict (Cassini-Huygens
anchor still matched the Saturn Titan-Rhea-Titan (1,1) candidate after #346's
Davis-2018 anchor tightening). Choose between path (a) source-acquire + anchor
body_set tightening OR path (b) `topology_label` schema extension.

**Verdict:** Schema extension landed (path b). User acquired the three relevant
PDFs (Strange-Russell-Buffington AAS-07-277, Yam-Davis-Longuski JSR 2009,
Valerino SpaceOps 2014), and the deep-read confirms Cassini's tour topology is
structurally distinct from a `(k1, k2)` repeated-moon cycler. The two changes
compose: schema gives us a per-anchor topology field, the deep-read gives us
the right labels.

## Decisive evidence from the new PDFs

**Strange, N., Russell, R., Buffington, B. (2007). "Mapping the V-infinity
Globe." AAS 07-277, JPL/Caltech.** Abstract opening sentence:

> This paper presents a graphical method for the design of transfers between
> the same gravity-assist body (i.e., **same-body transfers**). This graphical
> method collapses a large and complex space of possible trajectories to a map
> on which a tour designer may use intuition and experience to design a
> gravity-assist tour. **This method was used with great success in the Cassini
> extended mission design.**

The V-infinity-globe method targets **same-body transfers** (Titan→Titan
pumping), not repeated-moon (k1, k2) cycler topology. This IS Cassini's tour
topology.

**Yam, Davis, Longuski, Howell, Buffington (2009). "Saturn Impact Trajectories
for Cassini End-of-Mission." JSR, DOI 10.2514/1.38760.** Abstract:

> To impact Saturn with short-period orbits, a series of **successive Titan
> flybys** are required to increase inclination and decrease periapsis, while
> simultaneously avoiding the rings and mitigating ΔV expenditures. To ensure
> that the spacecraft is not prematurely damaged by material in the rings,
> **Tisserand graphs** are employed to determine when the ring-plane crossing
> distance is within the F–G ring gap...

Confirms: Cassini end-of-mission is **successive Titan** flybys + Tisserand
graphs. Same-body Titan pumping, NOT Titan-Rhea-Titan repetition. Rhea / Dione
/ Enceladus / Iapetus were single-visit science targets reached opportunistically
during the Titan pump phase, not repeating tour members in the cycler sense.

**Valerino (2014). "Updating the Reference Trajectory for the Cassini Solstice
Mission." SpaceOps 2014, DOI 10.2514/6.2014-6.2014-1880.** Confirms the
Titan-flyby reference trajectory was the tour's backbone during Solstice phase.

## Schema extension (the structural fix)

`CorpusAnchor` and `CandidateSignature` both gain an optional
`topology_label: frozenset[str] = frozenset()`. The matcher
(`_candidate_anchors`) gets one extra filter:

```python
if (
    sig.topology_label
    and anchor.topology_label
    and not (sig.topology_label & anchor.topology_label)
):
    continue
```

**Backward compatibility:** empty topology set on either side falls through
to body-set-only matching — the historical behaviour. All 41 existing anchors
that don't yet declare a topology continue to match exactly as before. New
behaviour only activates when BOTH sides opt in.

**Standard labels (documented in the dataclass docstring):**

* `"repeated-moon"` — Aldrin / (k1, k2) repeating-encounter cyclers
* `"pump-tour"` — V∞-leveraging same-body pump (Cassini Titan-pump)
* `"mga-tour"` — non-repeating multi-flyby tours (Galileo VEEGA, Heaton-Longuski U00-01)
* `"tulip"` — Sundman / petal Np-petal periodic orbits
* `"halo"`, `"nrho"`, `"resonant"`, `"binary-coorbital"`

## Anchor annotations landed in this commit

* **Cassini-Huygens Saturn-Titan satellite tour design**:
  `topology_label = frozenset({"pump-tour", "mga-tour"})`. Citation updated
  to lead with Strange-Russell-Buffington AAS-07-277 (the canonical V∞-globe
  paper); doi populated to 10.2514/1.38760 (Yam et al. 2009 JSR is the
  end-of-mission reference). authors tuple widened to include the AAS-07-277
  + JSR-2009 + SpaceOps-2014 author rolls.
* **Davis-Phillips-McCarthy Saturnian Ocean Worlds orbiters (2018)**:
  `topology_label = frozenset({"tulip", "halo", "nrho"})`. Belt-and-suspenders
  with the #346 body_set tightening — a candidate with
  `topology={"repeated-moon"}` is now excluded on both axes.

No other anchors annotated in this pass (`topology_label` defaults to
`frozenset()`). Future Phase 2 work may backfill labels on other anchors
when their topology is decisively sourceable.

## Stage A re-run after #349 lands

Driver `scripts/verify_344_titan_rhea_silver.py` now emits the candidate
signature with `topology_label = frozenset({"repeated-moon"})`. Re-run
result (artifact: `data/silver_344_verified.jsonl`):

| # | Sub-gate | Pre-#349 | Post-#349 |
|---|---|---|---|
| 1 | IC residual | 0.01018810 km/s — PASS | 0.01018810 km/s — PASS |
| 2 | Lit-fresh anchor count | **1** — FAIL | **0** — PASS |
| 3 | Physical sanity | Titan 49.91° / Rhea 7.07° / Titan 50.27° — PASS | same — PASS |
| 4 | ML flagger | p_fp=0.604 "real" — PASS | same — PASS |
| Verdict | | **HALT** | **PASS_PROCEED_TO_STAGE_B** |

`#344` is now eligible to enter the V1 → V2 → V3 → V4 ladder (modeled on the
#327 → #340 sequence that admitted the Umbriel-Oberon SILVER as the catalogue's
first computed `quasi_cycler` at #339).

## Tests + backward-compat verification

* `tests/search/test_literature_check.py`: 18/18 pass — all existing
  signature/anchor fixtures continue to match exactly as before (none of them
  declare `topology_label`, so they fall through the new filter).
* `uv run ruff check`, `uv run ruff format --check`, `uv run mypy` — all
  pass (pre-commit verified the same).

## Phase 2 follow-ons (not blocking)

* **Topology backfill across `KNOWN_CORPUS`:** the 41 anchors not yet
  annotated could benefit from explicit topology labels (the Aldrin family
  is clearly `repeated-moon`; Howard et al. Persephone Pluto is `mga-tour`
  + `halo`; the Pluto-Charon CR3BP tadpole/horseshoe anchor is
  `binary-coorbital`; etc.). Each backfill needs a sourced citation
  supporting the claim, not just inference. Track as a separate batch task
  after #344 Stage B-E.
* **Topology emission across discovery paths:** the production scan paths
  (`search/discovery_campaign.py`, `scripts/scan_320_*`, etc.) should emit
  topology labels on their `CandidateSignature` objects so the matcher's
  topology filter fires in routine scans, not just in the #344 Stage A
  driver. Trackable as a follow-on cleanup.

## Files in this commit

* `src/cyclerfinder/search/literature_check.py` — schema + 2 anchor
  annotations + matcher filter (the load-bearing change)
* `scripts/verify_344_titan_rhea_silver.py` — Stage A driver emits
  `topology_label={"repeated-moon"}` on the candidate
* `data/silver_344_verified.jsonl` — refreshed Stage A artifact with the
  new PASS verdict
* `docs/notes/2026-06-17-349-cassini-anchor-topology-label.md` — this note
