# Sourced ΔV-band definitions for the cycler taxonomy (#415)

**Date:** 2026-06-22. **Status:** literature synthesis + proposal. Review-only —
no `catalogue.yaml` / schema / `validate.py` edits (those are census-ratchet
sensitive and concurrently owned; this note recommends, a later coordinated task
applies). Every threshold below is tagged **[sourced]** (traces to a cited paper)
or **[project-convention]** (no published threshold exists; we chose it).

Companion: `docs/notes/2026-06-22-415-dv-band-definitions-task-brief.md` (the
task), and the raw extraction it was built from (corpus pass, 16 files).

## TL;DR

1. The literature's **"ballistic" is a geometric property of the idealized
   circular-coplanar model** (AR ≥ 1 *and* TR ≥ 1), **not** a real-ephemeris
   zero-ΔV claim. The two diverge — sometimes in opposite directions (Aldrin is
   geometrically *powered*, TR = 0.86, yet has a real-ephemeris **0 m/s** launch
   window; S1L1 is geometrically *ballistic* yet needs ~10 m/s in the ephemeris).
2. The only **sourced quantitative ΔV tiers** are Russell 2004's:
   **< 1 / < 10 / < 300 m/s total deterministic maneuver over a 7-cycle
   real-ephemeris propagation** (Table 5.5/5.6). Everything else in the corpus is
   either a geometric ratio, a one-off mission number, or a *computed surrogate*.
3. ΔV figures across sources are **not comparable** without fixing two things:
   the **basis** (per-orbit / per-cycle / per-synodic / total-over-N-cycles /
   averaged-over-launch-windows) and the **launch-window policy** (best vs
   average diverge ~100× for the same cycler). Any band MUST state both.
4. **Establishment ΔV ≠ maintenance ΔV** — the corpus routinely reports both;
   they must never be folded into one band. The catalogue already separates them
   (`v_infinity_leveraging_dv_kms` vs `delta_v_kms`); keep that.

## Two orthogonal axes (the key structural finding)

"Ballistic vs powered" is **not one scale** — it's two independent axes, and
conflating them is the source of most of the corpus's apparent contradictions:

- **Axis A — geometric feasibility (idealized model).** Can the cycle close with
  unpowered flybys in the circular-coplanar model? Sourced criterion: **AR ≥ 1**
  (reaches Mars) **and TR ≥ 1** (all flyby turns achievable at ≥ 200 km Earth
  altitude) ⇒ "strictly ballistic" [sourced: Russell-Ocampo 2003 p.13; Russell
  2004 §3.8 footnotes a/b; McConaghy 2002 Table 4 fn e: "required turn < max
  possible turn"]. This axis is **launch-window-blind** and **real-ephemeris-blind**.
- **Axis B — real-ephemeris deterministic ΔV.** What does it actually cost to fly
  the cycle in the true (eccentric, inclined, epoch-dependent) ephemeris? This is
  the engineering quantity, and it is what #388 measures.

A cycler can be (A-ballistic, B-cheap) [S1L1 at its window], (A-powered,
B-cheap-at-some-window) [Aldrin Aug-2003], or (A-ballistic, B-needs-ΔV) [S1L1
inbound]. **The catalogue must record both axes, not collapse them.**

## Proposed bands

Bands are defined on **Axis B** (real-ephemeris deterministic ΔV), with **Axis A**
kept as a separate descriptor. Axis-B basis is fixed below.

| Band | Definition (Axis B) | Threshold | Tag |
|------|--------------------|-----------|-----|
| **strictly-ballistic** | total deterministic maneuver over a 7-cycle real-ephemeris propagation, at the best launch window | **< 1 m/s / 7 cycles** | [sourced] Russell 2004 Table 5.5 ("completely ballistic"; 9 parents) |
| **essentially-ballistic** | same basis | **< 10 m/s / 7 cycles** | [sourced] Russell 2004 ("essentially ballistic for all launch dates"; 39 parents); corroborated by McConaghy 2006 S1L1 ~10 m/s/30 yr |
| **low-maintenance cycler** (≈ "near-ballistic") | same basis | **< 300 m/s / 7 cycles** | [sourced] Russell 2004 (widest low-maneuver net; 74 parents) |
| **powered cycler (DSM)** | discrete deterministic maneuvers above the low-maintenance ceiling, cycle still repeats | **≥ 300 m/s / 7 cycles**, impulsive | [project-convention] (boundary = top of Russell's net; McConaghy 2006 DSM trades 0.11–0.81 km/s/synodic sit here) |
| **low-thrust / SEP** | continuous propulsion over transit legs; Lambert inapplicable | n/a (regime, not magnitude) | [sourced] data/README §trajectory_regime; out of v1 scope |

Geometric descriptor (Axis A), recorded alongside, **not** as the band:
`AR ≥ 1 & TR ≥ 1` ⇒ geometric-ballistic; `AR ≥ 0.9 & TR ≥ 0.85` ⇒
geometric-near-ballistic [sourced: Russell 2004 dissertation cutoffs — see
conflict C1].

### Mandatory basis + window policy (resolves the comparison problem)

- **Basis:** Russell's native **"total over 7 cycles"** is the anchor because it
  is the only sourced quantitative tiering. When a per-synodic figure is needed,
  convert *explicitly* using the row's own cycle period ↔ synodic-period ratio —
  **do not assume 1 cycle = 1 synodic** (Aldrin's cycle ≈ several synodic
  periods; a naive division mis-scales the band). Flag the conversion per row.
- **Launch-window policy:** report the **best-window** value (matches Russell
  Table 5.5) AND the **21-window average** as a robustness flag — they diverge
  ~100× (Aldrin: 0 m/s best vs 3297 m/s avg). A cycler that is only cheap at
  isolated windows ("spiky") is materially different from one cheap everywhere
  ("essentially ballistic"); the band should carry that flag.

## Conflicts found (carried from the extraction — surfaced, not smoothed)

- **C1 — near-ballistic cutoff differs between Russell sources.** 2003 paper:
  TR_MIN = **0.9**; 2004 dissertation: TR_MIN = **0.85** (AR_MIN = 0.9 both).
  Both authors call it "arbitrary." Catalogue adopts **0.85**. *Resolution:* use
  0.85 (dissertation, the catalogue's existing choice) and cite it as such; note
  it is an arbitrary documentation filter, not a physical boundary.
- **C2 — geometric "ballistic" ≠ real-ephemeris zero-ΔV** (the Axis-A/B split
  above). *Resolution:* two-axis encoding; never infer Axis B from Axis A.
- **C3 — ΔV basis inconsistent** (per-orbit / per-cycle / per-synodic /
  total-7-cycles / avg-over-windows). *Resolution:* fixed basis + window policy
  above.
- **C4 — same cycler varies ~100× by maneuver model** (our GMAT runs: S1L1 62 m/s
  single-seed vs 7.29 km/s node-anchored per-flyby; Aldrin 0.175 vs 2.91 km/s
  with/without Oberth). These are *our computed* values, not literature, but they
  prove the maneuver model (Oberth credit, single-seed vs node-anchored) must be
  recorded with any Axis-B number. *Resolution:* band figures must state the
  maneuver model; prefer the optimized (Oberth-crediting, free-node) value.
- **C5 — establishment vs maintenance conflated in sources.** *Resolution:*
  separate axes; bands are maintenance-only.
- **C6 — free-return ≠ cycler.** Tito/Okutsu "ballistic free returns" have zero
  post-injection ΔV but do not repeat (can't earn ≥ V2). *Resolution:* the bands
  apply to *repeating* cyclers; free-returns are a different object class.
- **C7 — catalogue `delta_v_kms` is 3-state** (null / 0.0 / positive; absence ≠
  ballistic). *Resolution:* band logic must respect null ≠ 0.0.

## Encoding recommendation (proposal — not applied here)

1. Add an explicit **Axis-A geometric descriptor** field (or reuse AR/TR
   invariants) distinct from the Axis-B ΔV band, so `trajectory_regime` stops
   carrying two meanings.
2. Define the Axis-B band on the **< 1 / < 10 / < 300 m/s per-7-cycle**, best-window
   basis above, with the maneuver model + window policy recorded per row.
3. Map to validation tiers: **strictly/essentially-ballistic → V1-ballistic
   eligible**; **low-maintenance/powered → V2-powered** (the Aldrin tier, with a
   *quantified, bounded* maintenance ΔV); **low-thrust → out of v1**.
4. Keep `v_infinity_leveraging_dv_kms` (establishment) strictly separate from
   `delta_v_kms` / `maintenance_dv_kms_per_synodic` (maintenance).

## Rows that look mis-binned under the sourced definition (list only — no edits)

- **`mcconaghy-2006-em-k2`** — tagged `trajectory_regime: ballistic`, but
  McConaghy 2006 p.461 states S1L1 is "**only NEARLY ballistic in the ephemeris
  model**" (~10 m/s/30 yr outbound; inbound needs a modest DSM), and the #388
  family-pinned run shows it does not close zero-ΔV in DE440. Under the bands it
  is **essentially-ballistic (Axis B), geometric-ballistic (Axis A)** — a
  candidate for **V2-powered with a quantified maintenance ΔV**, not strict
  V1-ballistic. (This is the direct #388 payoff: the promotion target is V2-powered.)
- **`aldrin-cycler`** — `trajectory_regime: powered` (Axis A: TR = 0.86) is
  correct *geometrically*, but Russell 2004 finds a **0 m/s/7-cycle** real-eph
  window (Aug 6 2003) — so its Axis-B band is window-dependent ("spiky":
  essentially-ballistic at repeat windows, ~3.3 km/s average). The row should
  carry the best-window Axis-B value + the spiky flag, not just the powered tag.
- The **`russell-ch4-*`** V3 rows carry geometric AR/TR (Axis A) but no recorded
  Axis-B real-ephemeris ΔV; #388 shows they relax off-anchor / high-energy. They
  need an Axis-B measurement before any band assignment.

## Honesty notes

- The only sourced ΔV *thresholds* are Russell 2004's < 1 / < 10 / < 300 m/s
  per-7-cycle tiers; the "powered" lower bound (≥ 300 m/s) is **project-convention**
  (top of Russell's net), not a published cutoff.
- Several headline ΔV numbers in the corpus are **computed surrogates**
  (catalogue `delta_v_kms` M2-flyby surrogate; `maintenance_dv_kms_per_synodic`
  turn-deficit over-estimate) or **our GMAT reproductions**, NOT source-attested.
  They must never become golden-test targets; only the geometry (AR/TR, turn
  angles) and the sources' own ΔV figures are golden.
- Rauwolf's 1500 m/s/15-yr SEP-maintenance figure (via Howe 2025) traces to an
  **uncaptured** paper (Rauwolf-Friedlander-Nock 2002, AIAA 2002-5046) — treat as
  unverified until that source is acquired.

## Multi-sourcing status — the ΔV thresholds are under-cited (acquisition needed)

Honest assessment of the evidence base, by axis:
- **Axis A (geometric ballistic criterion, AR ≥ 1 / TR ≥ 1):** **well multi-sourced**
  — Russell-Ocampo 2003, Russell 2004, and McConaghy 2002 independently state the
  same turn-angle/aphelion criterion. Solid.
- **Axis B (real-ephemeris ΔV *magnitude* thresholds, < 1 / < 10 / < 300 m/s
  per-7-cycle):** rests on **one primary source** (Russell 2004 Table 5.5/5.6).
  This is the weak point — single-primary-sourced numeric cutoffs.

A 2026-06-22 web pass for independent corroboration found (leads, NOT yet pinned
to citable primaries — do not treat as golden):
- The broader Mars-cycler literature cites **annual station-keeping of ~20–100
  m/s/yr** to realign cyclers (matches Russell's tens-of-m/s near-ballistic tier
  in magnitude) — but the figure surfaced via a secondary aggregator, not a
  primary; its source must be traced before use.
- McConaghy 2006 (the S1L1 / `mcconaghy-2006-em-k2` paper) is independently
  confirmed: 153-day legs, V∞ 4.7/5.0 km/s, "total required ΔV very small, though
  not zero" — corroborates *essentially-ballistic*, not strict-ballistic.

**Acquisition targets to give Axis B ≥ 3 independent anchors** (currently 1):
1. **Genova-Aldrin, "Cycler Orbits and Solar System Pony Express," JSR (2021),
   DOI 10.2514/1.A35091** — reports impulsive + SEP cycler-maintenance ΔV
   sequences (paywalled; acquire via the private corpus / institutional access).
2. **The primary source of the "~20–100 m/s/yr station-keeping" figure** — trace
   from the aggregator to a citable cycler review / station-keeping paper.
3. **Rauwolf-Friedlander-Nock 2002 (AIAA 2002-5046)** — the 1500 m/s/15-yr SEP
   maintenance figure (flagged uncaptured above; #384-adjacent acquisition).

Until ≥ 2 independent primaries corroborate the Axis-B cutoffs, the < 1 / < 10 /
< 300 m/s tiers should be cited as **"Russell 2004 (sole primary); magnitude
consistent with the broader literature (~tens of m/s station-keeping), pending
multi-source confirmation."** The band *structure* (two axes + the geometric
criterion) is firm; the *numeric Axis-B cutoffs* are provisional on acquisition.

## References (all [sourced] items above)
- Russell-Ocampo 2003 (`digest-russell-ocampo-2003.md`): AR/TR criteria, TR_MIN=0.9.
- Russell 2004 dissertation (`russell-2004-member-tables-transcription.md`,
  `russell-2004-continuation-deepdive.md`): TR_MIN=0.85; Table 5.5/5.6 the
  < 1/< 10/< 300 m/s/7-cycle tiers; Aldrin 0-m/s window; S1L1 0-m/s window.
- McConaghy 2002 (`digest-mcconaghy-2002.md`): required-turn < max-turn; n=7 ballistic.
- McConaghy 2006 (`digest-mcconaghy-2006.md`): S1L1 "nearly ballistic ~10 m/s",
  DSM per-synodic trades.
- Byrnes-Longuski-Aldrin 1993 (`digest-byrnes-longuski-aldrin-1993.md`):
  Aldrin 1.73/2.04 km/s/15-yr; ~230 m/s/orbit idealized.
- Genova-Aldrin 2015, Rogers 2015, Tito 2018, Okutsu — corroborating / boundary.
- `data/README.md` §trajectory_regime, §delta_v_kms (field semantics, 3-state).
