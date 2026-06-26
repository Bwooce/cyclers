# #470 moon-tour reproduction admission wave — V2 verdict: 0/10 admitted (clean negative)

**Date:** 2026-06-26
**Task:** #470 (admit the 10 in-band moon-tour reproductions from #468 as
known-reproduction catalogue rows)
**Outcome:** **0 of 10 admitted.** Every proposal CLOSES IN-BAND (the #468
single-window powered leveraging close) but **FAILS the V2 moontour gauntlet**
(spec §14 ≥3-cycle bounded-drift). Per the admission procedure (admit only
V2-passers; report V2-failers honestly), nothing is written to `catalogue.yaml`.
This is the repo's "clean negative is success" discipline.

## What was tested

The 10 in-band proposals in `data/admission_proposals_468.jsonl`:

* **Jovian (5):** Io-Europa-Ganymede-Io, Europa-Ganymede-Callisto-Europa,
  Io-Europa-Io, Ganymede-Callisto-Ganymede, Callisto-Ganymede-Europa-Callisto
  (reproduce Hernandez/Jones/Jesick 2017 IEG triple cyclers + Liang et al. 2024
  CGE triple cyclers, JGCD DOI 10.2514/1.G008387).

  > **Citation grounding (#483/#485 audit, 2026-06-26).** The Jovian reference
  > here is the JOVIAN-MOON paper *"One Class of Io-Europa-Ganymede Triple
  > Cyclers"* by Hernandez, Jones & Jesick (AAS/AIAA Astrodynamics Specialist
  > Conf., Stevenson WA, Aug 2017; Adv. Astronaut. Sci. 162, pp. 973-984). It is
  > **NOT** the same authors' heliocentric paper AAS 17-577 *"Low Excess Speed
  > Triple Cyclers of Venus, Earth, and Mars"* — a `{Venus, Earth, Mars}`
  > interplanetary work. The two share the phrase "triple cycler" (a
  > concept-collision across systems) and the same JPL author team; the original
  > admission citation was attached by that concept match rather than by opening
  > the source, the exact failure class in memory
  > `feedback_ground_citations_against_content`. Ground-truthed 2026-06-26: the
  > on-disk PDF `cyclers_pdf/papers/jones-hernandez-jesick-2017-low-excess-speed-
  > vem-triple-cyclers-AAS-17-577.pdf` is VEM-only (title page), and the distinct
  > Jovian IEG paper is independently confirmed (Semantic Scholar
  > `7e1de63096852b5422107ffc23a9312ea3de54f3`). The citation as printed above is
  > body/system-CORRECT (Jovian → Jovian work); this note records that it is now
  > GROUNDED, not concept-matched. The `data/admission_proposals_468.jsonl`
  > `matched_source` strings and the catalogue row
  > `hernandez-2017-jovian-ieg-triple-family` were independently verified to cite
  > the Jovian work, not AAS 17-577. Full sweep: `docs/notes/2026-06-26-citation-
  > audit.md`.
* **Saturnian (5):** Titan-Rhea-Dione-Titan, Rhea-Dione-Tethys-Rhea,
  Dione-Tethys-Enceladus-Dione, Rhea-Dione-Rhea, Tethys-Enceladus-Tethys
  (reproduce the Cassini icy-moon leveraging tour, Wolf & Smith 1995 /
  Campagnola-Strange-Russell endgame).

For each, the WINNING #468 sweep config (cheapest in-band powered close — the
exact `(ToF-scale, phasing)` `campaign_468_multirev_tour._run_skeleton`
selected) was reconstructed and fed to the **official**
`cyclerfinder.data.validation.v2_moontour.run_v2_moontour` with the
`MultiRevLeveragingReleg` powered backend over 3 cycles.

## Result — independently cross-checked, 0/10 pass

`run_v2_moontour` verdicts (drift floor 50,000 km; closure floor 0.05 km/s):

| sequence | n_cycles | max_drift_km | max_resid_km/s | passes_v2 |
|---|---|---|---|---|
| Io-Europa-Ganymede-Io | 3 | 832,718 | 19.40 | **fail** |
| Europa-Ganymede-Callisto-Europa | 3 | 1,168,959 | 1.04 | **fail** |
| Io-Europa-Io | 3 | 802,265 | 11.77 | **fail** |
| Ganymede-Callisto-Ganymede | 3 | 2,140,770 | 0.00 | **fail** |
| Callisto-Ganymede-Europa-Callisto | 3 | 3,593,638 | 1.90 | **fail** |
| Titan-Rhea-Dione-Titan | 3 | 1,446,990 | 8.54 | **fail** |
| Rhea-Dione-Tethys-Rhea | 3 | 826,693 | 10.41 | **fail** |
| Dione-Tethys-Enceladus-Dione | 3 | 700,460 | 12.71 | **fail** |
| Rhea-Dione-Rhea | 3 | 916,116 | 0.00 | **fail** |
| Tethys-Enceladus-Tethys | 3 | 164,478 | 1.50 | **fail** |

Every tour blows the 50,000 km inter-cycle rendezvous-drift floor by **3×–72×**
(smallest is Tethys-Enceladus at 164,478 km; largest is Liang CGE at 3.6 M km).
The closure-residual floor is also blown on most tours where the ballistic
continuity is measured.

A private mirror driver (`scripts/_v2_admit_470.py`, campaign phasing advanced
by natural mean motion each cycle) reproduces the same 0/10 verdict, confirming
the negative is the gauntlet's, not a phasing-convention artifact.

## Interpretation

The #468 in-band closes are **single-window, epoch-locked patched-conic
geometries**: the campaign found ONE phasing where the powered multi-rev
leveraging cycle closes in-band (< 3.5 km/s/cycle). But the moons drift out of
that phasing within 1–2 cycles — the tours are **not in a resonance lock that
re-phases the encounter geometry each cycle**. They therefore do not qualify as
repeating cyclers / quasi-cyclers on the spec §14 ≥3-cycle bounded-drift
horizon: closed-in-band ≠ V2-stable.

This matches the existing taxonomy: a one-window leveraging tour is at best an
`mga_tour`-shaped artifact, not a cycler; and an `mga_tour` is the non-repeating
class the proposals were originally (correctly, as it turns out) labelled. The
"reclassify to cycler/quasi_cycler" step was contingent on V2-repeating
behaviour, which the gauntlet falsifies.

## Action

* **Admitted: 0.** No `catalogue.yaml` rows added; no census ratchets touched.
* Full suite `uv run pytest tests/data tests/search` green pre- and post-work
  (1538 passed, 5 xfailed [pre-existing #54/Hiraiwa], 1 xpassed [Kumar]).
* The #468 ledger and proposals stand as an honest record that these published
  tours close in-band under our powered leveraging method but are **not V2-stable
  multi-cycle cyclers in the circular-coplanar moon-ephemeris model** — a
  model-limitation negative (the published tours live in the real-ephemeris,
  resonance-locked regime; our single-window conic close cannot hold phasing).

Driver: `scripts/_v2_admit_470.py`,
cross-check: `scripts/_v2_official_check_470.py`.
