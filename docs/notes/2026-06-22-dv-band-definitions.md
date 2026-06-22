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
| **strictly-ballistic** | total deterministic maneuver over a 7-cycle real-ephemeris propagation, at the best launch window | **< 1 m/s / 7 cycles** | [sourced] Russell-Ocampo 2006 JGCD abstract ("nine parent cyclers … less than 1 m/s over seven full cycles") = Russell 2004 Table 5.5 (9 parents) |
| **essentially-ballistic** | same basis | **< 10 m/s / 7 cycles** | [sourced] Russell-Ocampo 2006 JGCD ("39 … parent cyclers … less than 10 … m/s") = Russell 2004 (39 parents); magnitude corroborated by McConaghy 2006 S1L1 ~10 m/s/30 yr |
| **low-maintenance cycler** (≈ "near-ballistic") | same basis | **< 300 m/s / 7 cycles** | [sourced] Russell-Ocampo 2006 JGCD ("74 parent cyclers … less than … 300 m/s") = Russell 2004 (74 parents) |
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

## Multi-sourcing status (updated 2026-06-22 after a corpus-wide ΔV extraction)

A 2026-06-22 pass read the **local primary PDFs** (`cyclers_pdf/papers/`, not just
the digests) + the user-supplied Rauwolf 2002. Evidence base, by claim:

- **Axis A (geometric ballistic criterion, AR ≥ 1 / TR ≥ 1):** **well
  multi-sourced** — Russell-Ocampo 2003, Russell 2004, McConaghy 2002 agree. Solid.
- **The ballistic *concept* (zero deterministic ΔV is achievable for a good
  cycler):** **multi-sourced** — Russell 2004, McConaghy 2006 (S1L1), Patel 2019
  ("no Δv maneuvers necessary … 'Ballistic S1L1 Cycler'"), Landau-Longuski 2006
  ("ballistic … DSM ΔV = 0 versions exist"). Robust, qualitatively.
- **Powered-cycler maintenance ΔV magnitude (~0.5–2 km/s per synodic / per 15 yr):**
  **now firmly multi-sourced (5+ primaries):** Rauwolf 2002 (Aldrin 1.561/1.605
  km/s per 15 yr, 3 of 7 orbits — `digest-rauwolf-friedlander-nock-2002.md`),
  Byrnes-Longuski-Aldrin 1993 (1.73/2.04 km/s/15 yr), Chen 2002 (Aldrin "≈ 0.54
  km/s per synodic"), Chen-McConaghy-Landau-Longuski-Aldrin 2005 (0.58–1.05 km/s
  per synodic per vehicle), Landau-Longuski 2006 (0.78 km/s avg cycler DSM). All
  agree at the ~km/s-per-15-yr / hundreds-of-m/s-to-~1-km/s-per-synodic scale.
- **Near-ballistic *m/s cutoffs* (< 1 / < 10 / < 300 m/s per-7-cycle):** now in a
  **peer-reviewed JGCD primary** — **Russell-Ocampo 2006 JGCD 29(2)** states the
  tiers in its abstract verbatim ("nine parent cyclers … less than 1 m/s over
  seven full cycles … 39 and 74 parent cyclers … less than 10 and 300 m/s,
  respectively"), with the **same 9/39/74 parent counts** as Russell 2004 Table
  5.5/5.6. CAVEAT: this is the *published version of the same author's
  dissertation work*, not an independent re-derivation — so the cutoffs are now
  **publication-grade (JGCD, peer-reviewed)** but still rest on the **single
  Russell line of work** in this exact basis. (This corrects the earlier draft,
  which wrongly called Russell-Ocampo 2006 "NOT in the corpus" — it was filed,
  misnamed `russell-2006-systematic-method-design`, and already digested.) An
  independent second-author same-basis source is still not in hand:
  - **Pascarella/Pony-Express AAS-22-015** (the paper I'd mis-listed as the
    paywalled JSR) corroborates the *magnitude* — real-ephemeris maintenance is
    cheap ("< 5 kg … ~2 kg of propellant" over 8 flybys / ~6 yr; ≈ 175 m/s/6 yr
    by back-calc) — but reports **propellant mass, not ΔV in m/s**, so it cannot
    be band-mapped to the cutoffs.
  - McConaghy 2006 S1L1 "~10 m/s / 30 yr" corroborates the *essentially-ballistic*
    tier magnitude. Also consistent, also not the same basis.

  - **NEW (2026-06-22, #416 corpus-wide + Friedlander acquisition):** the ballistic
    *floor* and the **~200–300 m/s "differentially-correctable-to-ballistic"
    boundary** are now corroborated by **genuinely independent** sources, in three
    systems: **Friedlander-Niehoff-Byrnes-Longuski 1986** (VISIT "free orbit",
    pre-Russell — the earliest primary), **Jones-Hernandez-Jesick 2017** (VEM triple,
    JPL/Jesick, real-eph: "velocity increments below 200 m/s … differentially
    corrected … to be entirely ballistic"), and **Liang 2024** (Callisto-Ganymede-
    Europa, real-eph: max 1.04e-7 m/s/cycle). See `2026-06-22-416-cycler-dv-universe-sweep.md`.

**Net (#416 CLOSED, 2026-06-22):** the powered band and the ballistic concept are
genuinely multi-sourced; the **near-ballistic band is now well-grounded** — its
*floor* and its *~200–300 m/s boundary* are corroborated across independent author
groups and three systems (Friedlander 1986, Jones 2017 VEM, Liang 2024 CGE,
McConaghy 2004/2006). The **exact < 1 / < 10 / < 300 m/s tiered census (9/39/74
parents)** remains **Russell-Ocampo's own published quantification** (JGCD 29(2)
2006 = Russell 2004 diss.) — no other paper reproduces that specific parent-count
census in the same basis, but it is no longer a *fragile* single source. Cite the
cutoffs as: **"Russell-Ocampo 2006 JGCD 29(2) census (9/39/74 parents at
<1/<10/<300 m/s over 7 cycles, real-eph best window); ballistic floor + ~200–300
m/s correctable-to-ballistic boundary independently corroborated (Friedlander 1986;
Jones 2017 VEM; Liang 2024 CGE; McConaghy 2004/2006)."**

**Acquisition status: no gap remaining.** Friedlander 1986 (the last genuine cycler-
ΔV primary not in the corpus, the #116/#416 wishlist item) was user-supplied and is
now filed + digested (`2026-06-22-digest-friedlander-niehoff-byrnes-longuski-1986.md`).
Russell-Ocampo 2006 JGCD is in-corpus (was misfiled). Pony Express is **not** a
cutoff corroboration: Sanchez-Net 2022 is patched-conic and defers real-eph ΔV as
"large"; Pascarella 2022 is low-thrust fuel-optimal (~163 m/s, one candidate,
magnitude only). Rauwolf 2002 is a powered-band anchor (Aldrin), not near-ballistic.

## References (all [sourced] items above)
- Russell-Ocampo 2003 (`digest-russell-ocampo-2003.md`): AR/TR criteria, TR_MIN=0.9.
- Russell 2004 dissertation (`russell-2004-member-tables-transcription.md`,
  `russell-2004-continuation-deepdive.md`): TR_MIN=0.85; Table 5.5/5.6 the
  < 1/< 10/< 300 m/s/7-cycle tiers; Aldrin 0-m/s window; S1L1 0-m/s window.
- Russell-Ocampo 2006 JGCD 29(2), DOI 10.2514/1.13652
  (`2026-06-20-digest-new-papers.md` §3): peer-reviewed publication of the above
  tiers (9/39/74 parents at <1/<10/<300 m/s over 7 cycles); 203 parent cyclers.
- McConaghy 2002 (`digest-mcconaghy-2002.md`): required-turn < max-turn; n=7 ballistic.
- McConaghy 2006 (`digest-mcconaghy-2006.md`): S1L1 "nearly ballistic ~10 m/s",
  DSM per-synodic trades.
- Byrnes-Longuski-Aldrin 1993 (`digest-byrnes-longuski-aldrin-1993.md`):
  Aldrin 1.73/2.04 km/s/15-yr; ~230 m/s/orbit idealized.
- Genova-Aldrin 2015, Rogers 2015, Tito 2018, Okutsu — corroborating / boundary.
- `data/README.md` §trajectory_regime, §delta_v_kms (field semantics, 3-state).
