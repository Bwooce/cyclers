# #468 at-scale multi-rev leveraging discovery campaign — VERDICT

**Date:** 2026-06-26
**Capability spent:** `MultiRevLeveragingReleg` (#465) via
`search.releg_moontour.close_powered_cycle`
**Driver script:** `scripts/campaign_468_multirev_tour.py`
**Ledger:** `data/admission_proposals_468.jsonl` (10 proposed rows, human-gated)
**Summary:** `out/campaign_468_summary.json`
**Status:** SHIPPED — capability confirmed at scale; **HONEST OUTCOME (b): every
in-band tour is a REPRODUCTION of a published VILM/leveraging moon tour. NO novel
catalogue row.** Four Uranian/Neptunian skeletons re-stamped as stronger
powered-empty. No catalogue.yaml row self-admitted.

---

## 1. The campaign

Spent the #465 multi-rev leveraging chain across the repeated-moon endgame
skeletons of every in-band system. For each skeleton the driver swept ToF-scale
∈ {1.0, 1.2, 1.5, 1.8, 2.0} × phase-offset ∈ {0.0, 0.5, 1.0} (the #465 lesson:
the high-V∞ finite-chain stall is a phasing/ToF artifact, so sweep to find a
reachable phasing) and kept the cheapest IN-BAND (< 3.5 km/s/cycle) feasible
close. Each in-band closure was classified NOVEL vs REPRODUCTION via
`search.literature_check.check_literature` against the curated published-record
corpus.

**14 skeletons probed: 10 closed in-band, 4 structurally empty.**

| System | Sequence | ΔV/cycle (km/s) | V∞ target | In-band | Novel/Repro |
|---|---|---|---|---|---|
| Jupiter | Io-Europa-Ganymede-Io | 0.607 | 6.0 | YES | REPRODUCTION |
| Jupiter | Europa-Ganymede-Callisto-Europa | 0.274 | 4.0 | YES | REPRODUCTION |
| Jupiter | Io-Europa-Io | 0.124 | 5.0 | YES | REPRODUCTION |
| Jupiter | Ganymede-Callisto-Ganymede | 0.329 | 8.0 | YES | REPRODUCTION |
| Jupiter | Callisto-Ganymede-Europa-Callisto | 0.224 | 4.0 | YES | REPRODUCTION |
| Saturn | Titan-Rhea-Dione-Titan | 0.214 | 5.0 | YES | REPRODUCTION |
| Saturn | Rhea-Dione-Tethys-Rhea | 0.166 | 4.0 | YES | REPRODUCTION |
| Saturn | Dione-Tethys-Enceladus-Dione | 0.459 | 4.0 | YES | REPRODUCTION |
| Saturn | Rhea-Dione-Rhea | 0.047 | 6.0 | YES | REPRODUCTION |
| Saturn | Tethys-Enceladus-Tethys | 0.051 | 5.0 | YES | REPRODUCTION |
| Uranus | Ariel-Umbriel-Ariel | inf | — | EMPTY | structural |
| Uranus | Titania-Oberon-Titania | inf | — | EMPTY | structural |
| Uranus | Ariel-Umbriel-Titania-Ariel | inf | — | EMPTY | structural |
| Neptune | Proteus-Triton-Proteus | inf | — | EMPTY | structural |

All 10 in-band closures sit deep in the `powered_dsm` band, ~5–70× under the
3.5 km/s/cycle ceiling — the #465 multi-VILM-minimum economics hold at scale.

---

## 2. Novel-vs-reproduction split: 0 novel / 10 reproduction (Outcome b)

**Every in-band tour matched the published record. No lit-fresh candidate.**

- **Jovian (5 tours):** matched the **Hernandez-Jones-Jesick** "One Class of
  Io-Europa-Ganymede Triple Cyclers" (AAS/AIAA 2017) for the Io/Europa/Ganymede
  sequences and the **Liang-Yang-Li-Bai-Qin** Callisto-Ganymede-Europa triple
  cyclers (JGCD Engineering Note 2024, DOI 10.2514/1.G008387) for the
  Callisto/Ganymede/Europa sequences. The Galilean endgame is the canonical
  VILM-leveraging region (Campagnola-Russell "The Endgame Problem" Part-1/2);
  the #465 decomposition golden already anchors to its Table 1.
- **Saturnian (5 tours):** matched the **Cassini-Huygens Saturn-Titan satellite
  tour design** corpus anchor (`body_set` = {Dione, Enceladus, Iapetus, Mimas,
  Rhea, Tethys, Titan}; citation chains Wolf-Smith 1995, Strange-Russell-
  Buffington AAS 07-277, Yam et al. JSR 2009, Valerino SpaceOps 2014). A live
  WebSearch corroborated this directly: Campagnola/Strange/Russell explicitly
  publish the "generalized v∞ leveraging tour of Saturn's low-mass satellites
  Rhea, Dione, Tethys and Enceladus" (~0.5 km/s for an Enceladus orbiter) and
  the "Titan → Rhea → Dione → Tethys" endgame — exactly the icy-moon sequences
  the chain closed.

### Honesty note — a corrected false-novel (the load-bearing caveat)

The FIRST campaign run reported the 5 Saturnian tours as `not-found` → NOVEL.
That was a **false-novel artifact of a weak injected offline search**, NOT a real
novelty signal. The curated `_candidate_anchors` ALREADY returns the
Cassini-Huygens anchor for every Saturnian icy-moon sequence (its body set covers
them all and the `pump-tour` topology label intersects), but the offline
corpus_search hit was titled "tour" not "cycler", so the structural scorer's
mandatory `cycler/cyclic` floor was not met and the anchor never corroborated. A
live WebSearch immediately surfaced the published Saturnian leveraging-tour
record. Once the corpus hit carried the (correct) repeating-cycler framing, all 5
classified REPRODUCTION at confidence 0.95. **Lesson (memory updated): a
`not-found` is only as honest as the injected search — ground a mature-field
not-found with a live WebSearch and verify the curated anchor's body_set/topology
actually overlaps before believing novelty. A leveraging "tour" that revisits a
moon IS a (loose) cycler.**

---

## 3. Uranus / Neptune — stronger powered-empty re-stamps (Outcome 3)

The 4 Uranian/Neptunian skeletons are prefiltered EMPTY before any chain solve:
every adjacent-moon leg is unlinkable (disjoint Tisserand/resonance contours at
all probed V∞ up to 15 km/s). A powered multi-rev CHAIN that ALSO cannot bridge
disjoint contours is a STRONGER negative than the prior single-DSM/ballistic one
(`multi-rev-leveraging` ⊐ `one-dsm-per-leg` in the method-capability partial
order). Re-stamped into `data/empty_regions.jsonl` with
`multirev_leveraging_method_capability`:

- `uranus-ariel-umbriel-ariel-multirev-leveraging-2026-06-26`
- `uranus-titania-oberon-titania-multirev-leveraging-2026-06-26`
- `uranus-ariel-umbriel-titania-ariel-multirev-leveraging-2026-06-26`
- `neptune-proteus-triton-proteus-multirev-leveraging-2026-06-26`

This subsumes / corroborates the prior
`uranus-neptune-regular-moon-endgame-vilm-2026-06-23` negative with a more
capable method (chaining walks V∞ *within* a contour; it cannot jump disjoint
ones — the genome does not fabricate a bridge).

---

## 4. No self-admit — the admission-proposal ledger

Per discipline NO catalogue.yaml row is written by this campaign. Instead all 10
in-band tours (every one a reproduction) are PREPARED as proposed catalogue-row
records in `data/admission_proposals_468.jsonl` (id, class `mga_tour`, sequence,
ΔV/cycle, dv_band, target V∞, novel-vs-reproduction + the matched published
source/DOI, validation evidence). Each is necessary-not-sufficient: a closed
cycle + lit-check only — a row would still need to clear the V2 moontour gauntlet,
and since all 10 are reproductions the honest expectation is V0-known (a
published rediscovery), not a novel SILVER. The human (via the main agent)
reviews + admits; this campaign only proposes.

---

## 5. Outcome summary

- **(a) lit-fresh novel candidate:** NONE. The Galilean + Saturnian leveraging
  endgames are thoroughly mapped by the published record.
- **(b) all in-band tours reproductions:** YES — the realised outcome. The
  capability is confirmed at scale (10 in-band closures across 2 systems), but it
  re-derives known VILM/leveraging tours.
- **(c) finite-chain ceiling blocks systems:** The Uranian/Neptunian systems are
  blocked, but by the *disjoint-contour prefilter* (structural unlinkability),
  not the finite-chain V∞ ≲ 2·V_M reach ceiling — the chain never gets to solve
  them. Within Jupiter/Saturn the ToF/phasing sweep found reachable phasings for
  every probed skeleton (the #465 reach ceiling did not bite at scale).

**Net:** the multi-rev leveraging capability spends cleanly across the in-band
systems and confirms the #465 economics, but finds no novel catalogue row — the
in-band moon-tour endgame is a published-mapped region. A genuine, defensible
negative-for-novelty (capability-positive).

---

## 6. What shipped

| Artifact | Path |
|---|---|
| Campaign driver | `scripts/campaign_468_multirev_tour.py` |
| Admission-proposal ledger (10 rows) | `data/admission_proposals_468.jsonl` |
| JSON summary | `out/campaign_468_summary.json` |
| Empty-region re-stamps (4) | `data/empty_regions.jsonl` (appended) |
| Verdict | this doc |

Ratchet status: `tests/data` + `tests/search` GREEN (the re-stamps are
`validate_empty_region`-valid; the ledger is a scratch results file, not a
frozen-census input). No `catalogue.yaml` / `catalogue.schema.json` touched
(left to #305's mid-schema-bump).
