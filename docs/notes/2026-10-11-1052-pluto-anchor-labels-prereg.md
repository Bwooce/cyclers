# #1052 (Pluto part): pre-registration draft of the anchor labels and three controls (2026-10-11)

Status: DRAFT, written before any change and before any rerun. No edit of `literature_check.py` (held by
twobody-gen2-opus for #1056); this note is the plan for when it is free. Method as for the Earth-Moon part:
`docs/notes/2026-10-11-1047-earth-moon-families-literal-table-and-v3-gate.md` sec. 6-7. A label is the
trajectory type the paper itself computes or designs, in the gate vocabulary; the grounding is a digest section,
never memory. Where no digest exists, the grounding is the lowest level available and is said to be so.

All baselines below were run on the committed `literature_check.py` (`git show HEAD:...`, last commit f03cc22f,
saved to a scratch file and imported from there, so that the other agent's uncommitted edits do not matter).

## 1. Anchors in scope

Anchors whose primary is Pluto (the gate's `KNOWN_CORPUS`, 81 anchors): five. No other anchor reaches Pluto's
moons. The Russell-Strange anchors (Jovian, Saturnian, heliocentric) are the architecture precedent: the gate
returns "known-architecture-new-system" for a one-working-node signature at a system they do not cover. The
Ross and Roberts-Tsoukkas anchor is Earth-Moon only (sec. 4, decision D2).

| # | Anchor (key) | Current label | Provenance | `system` | Digest |
|---|---|---|---|---|---|
| 1 | Howett et al. 2021, Persephone (`howett-2021-persephone`) | none (deliberately empty, #881) | verified-against-source | pluto-charon | `2026-10-03-digest-howett-2021-persephone.md` |
| 2 | Stern, Tapley, Finley & Scherrer 2020, JSR (name "Stern/SwRI Pluto Game-Changer") | {mga-tour} | inherited-unverified | empty | `2026-06-23-digest-stern-2020-pluto-orbiter.md` |
| 3 | Showalter & Hamilton 2015, Nature 522:45 | {resonant} | inherited-unverified | empty | none |
| 4 | Brozovic et al. 2015, Icarus 246:317 (and AJ 163:241, 2022) | {ephemeris} | inherited-unverified | empty | none |
| 5 | "Tadpole/horseshoe at L3/L4/L5 (2025)", arXiv:2510.13479 | {binary-coorbital} | inherited-unverified | empty | none |

Working-body scope (`working_bodies_scope`) is unset on all five.

## 2. Proposed labels and grounding

| # | Proposed topology label | Working-body scope | Grounding |
|---|---|---|---|
| 1 Howett | {repeated-moon} (decision D1; the alternative is to leave it empty) | unset | Digest, "The trajectory content (Sec. 4.2)": "the science orbit consists of multiple periodic orbits designed in the Pluto/Charon restricted three-body dynamics model"; the orbits "enabling close encounters with Pluto and all its satellites"; "repeating ground tracks"; four orbits (Figs. 13-14); minor-satellite encounter velocity "<300 m/s". "What it does NOT contain": no initial conditions, no gravity-assist tour, no patched-conic or Lambert construction. The paper does not say there is a satellite encounter each period; "repeated-moon" is the nearest label in the vocabulary (periodic orbits with repeated satellite encounters), not a stated property. |
| 2 Stern | {mga-tour} unchanged; provenance to verified-against-source after a check of the title page against the digest | "one" | Digest, "What it is": "a Pluto-system orbital tour that maneuvers almost entirely on Charon gravity assists"; "Charon = primary GA maneuvering body, but NO numeric periapsis altitude here"; spec 2 (satellites observed "at radial distance <10,000 km") is observation, not a flyby. So Charon is the working body and the small moons are observed, not used: "one". The Finley et al. companion paper ("An Orbital Tour of Pluto and Its Moons") has no DOI and was never found published (`2026-06-15-279-finley-paper-disposition.md`), so it is not an anchor. |
| 3 Showalter-Hamilton | {ephemeris} (changed from {resonant}) | unset | No digest. Title (Crossref, 2026-10-11): "Resonant interactions and chaotic rotation of Pluto's small moons", Showalter and Hamilton, Nature 522, 45-49 (2015), DOI 10.1038/nature14469 (matches the anchor). The 2026-06-15 Pluto review (class 3) records the Styx-Nix-Hydra Laplace-like resonance as "a property of the moons' own orbits, not a spacecraft trajectory". That is the Brozovic rule (#881): a source about the bodies' own motion never matches a trajectory topology. Title-level plus review-note grounding only. |
| 4 Brozovic | {ephemeris} unchanged | unset | No digest. Title and DOI (Crossref, 2026-10-11): Brozovic, Showalter, Jacobson and Buie, "The orbits and masses of satellites of Pluto", Icarus 246, 317-329 (2015), DOI 10.1016/j.icarus.2014.03.015 (matches). The anchor's second reference, "AJ 163:241 (2022)", was not verified. The JPL SSD Pluto pages (accessed 2026-10-11) cite Brozovic and Jacobson (2024), "Post-New Horizons orbits and masses for the satellites of Pluto", AJ 167:256; the anchor's citation should add it once its DOI is checked. |
| 5 Tadpole (arXiv:2510.13479) | {binary-coorbital} unchanged | unset | arXiv abstract read 2026-10-11: "Chaotic Dynamics and Zero-Velocity Structures in the Pluto-Charon CR3BP", single author Abdul Wahab Jbara (independent researcher), submitted 2025-10-15, v2 2025-12-10; "compares the tadpole and horseshoe orbits of a massless particle with known low-mu orbits", RK4 simulations, "chaotic and unpredictable trajectories", "cannot host long-lived Trojan companions". The label is right. The anchor's citation is wrong: its title string ("Tadpole and horseshoe orbits in the Pluto-Charon CR3BP at L3/L4/L5") is not the arXiv title, and its authors field is a placeholder ("Pluto-Charon", "Charon"). Proposed correction: authors ("Jbara",), the real title, provenance verified-against-source at abstract level. Not a periodic-orbit paper. |

Wanted-list rows proposed (none of the four is in CORPUS_INDEX.md; grep for Showalter, Brozovi, Jbara, the
DOIs and the arXiv number found nothing):
1. Showalter, M. R. and Hamilton, D. P., Nature 522:45-49 (2015), DOI 10.1038/nature14469 (Crossref-confirmed 2026-10-11).
2. Brozovic, M., Showalter, M. R., Jacobson, R. A. and Buie, M. W., Icarus 246:317-329 (2015), DOI 10.1016/j.icarus.2014.03.015 (Crossref-confirmed 2026-10-11).
3. Brozovic, M. and Jacobson, R. A., AJ 167:256 (2024), DOI not checked (named by the JPL page).
4. Jbara, A. W., arXiv:2510.13479 (2025), free to fetch.

## 3. Decisions for the lead or owner

D1, Howett. The digest records the empty label as deliberate ("a narrower one would under-cover"). Today the
empty anchor swallows every Pluto signature, whatever its topology (baseline below: pump-tour, mga-tour and
resonant signatures all read "published" citing Howett). Choices:
- (a) {repeated-moon} (proposed). Keeps the class-level prior-art reading for any repeated-moon signature,
  which is what the digest wants (a human compares with Figs. 13-14), and stops Howett citing for tours and
  resonances it does not cover. Consequence for #1036: the 224 #998 candidates (signature repeated-moon, one
  working body) still read "published" through Howett, as today. This pass changes nothing for them.
- (b) leave empty: no change at all; every Pluto signature keeps reading "published".
- (c) a new anchor-only label for restricted-problem periodic orbits (like "transfer" and "model-transition"
  in #1047): a repeated-moon patched-conic signature then no longer matches Howett and would fall to
  "known-architecture-new-system" (the Russell-Strange architecture). This is the reading that contradicts the
  digest's class-level statement (consequence 2: "any Pluto-system periodic or repeating-encounter candidate
  should keep reading published"). Not recommended without an owner ruling.

D2, Ross and Roberts-Tsoukkas 2026 (arXiv:2606.29189). The digest (`2026-06-30-digest-ross-roberts-tsoukkas-2026-stable-ballistic-cyclers.md`)
says the mu-family spans 0.001 to 0.5, mu = 0.1 is adjacent to Pluto-Charon (0.10851), and "Pluto-Charon is NOT
in Table I or anywhere in the paper". No anchor represents it, and a two-primary CR3BP cycler is not a
patched-conic Charon-block signature. Recommendation: no new anchor now.

D3, Showalter-Hamilton {resonant} to {ephemeris}. Not the "usual" extra anchor edit: today a Pluto signature
labelled {resonant} can cite a paper about natural moon resonance. The change is the same as #881 made for
Brozovic. It rests on title-level grounding, so it is flagged for a digest (wanted row 1).

## 4. Baseline (committed gate, before any change)

Signature: primary Pluto, sequence (Charon, Nix), period_k 4, V_inf (0.178, 0.130) km/s, working_bodies "one".

| Signature | Result |
|---|---|
| topology {repeated-moon}, working_bodies one (a #998 type) | published 0.95, cites Howett; in-footprint anchors: Howett only |
| same, working_bodies two | published, Howett |
| sequence (Nix, Hydra), {repeated-moon} | published, Howett |
| topology {pump-tour} at Pluto | published, Howett (Howett swallows it) |
| topology {mga-tour} at Pluto | published, cites Howett (in-footprint: Howett and Stern; Howett is first in corpus order) |
| topology {resonant}, sequence (Nix, Hydra) | in-footprint: Howett and Showalter-Hamilton; published |
| topology {resonant} at Pluto | published, Howett |
| unlabelled at Pluto | published, Howett |
| the same {repeated-moon} signature with primary Jupiter (wrong system) | known-architecture-new-system, citing Russell-Strange |
| same with primary Saturn | known-architecture-new-system, citing Russell-Strange |

## 5. Controls (run first, after the edit; a failed control stops the run)

| # | Control | Expected after the proposed edit (D1 (a), D3) | Reason |
|---|---|---|---|
| C1 | positive: the #998 type above, {repeated-moon}, one working body | published, citing Howett | Howett would carry {repeated-moon}; it is the only Pluto anchor with that label (the others are mga-tour, ephemeris, ephemeris, binary-coorbital). Status is judged, the citation is recorded and not gated (the #1062 lesson). |
| C2 | negative, wrong system: the same signature with primary Jupiter | not published (known-architecture-new-system, as in the baseline) | No Jupiter anchor has Charon or Nix in its body set. |
| C3 | negative, topology: the same signature with {pump-tour} at Pluto | not published (not-found) | No Pluto anchor carries pump-tour or an empty label after the edit; baseline is "published", so this control fails before the edit and must pass after it. |

Additional checks, recorded not gated:
- {mga-tour} at Pluto, sequence (Charon, Nix): published, citing Stern (Howett no longer in the footprint); with
  working_bodies "two" Stern's scope "one" excludes it, giving not-found.
- {resonant} at Pluto, sequence (Charon, Nix): not-found (baseline: published via Howett only). Sequence (Nix, Hydra): not-found
  (baseline: published, Howett and Showalter in footprint); after D3 Showalter carries {ephemeris} and never matches a
  trajectory topology.
- Every other control of the gate (`scripts/litcheck_942_943_scope.py`: 26 controls, the negative and the five
  candidates) must give output identical to its before-run; none is at Pluto.
- If C1 fails (for instance if the footprint rule drops Howett for the label change), report it; do not move a
  label to make it pass.

## 6. What this pass does not do

- It does not clear or condemn any #998 candidate: under D1 (a) the 224 still read "published" (class level).
  The comparison with Howett's Figs. 13-14 stays a human step (#1036).
- No digest of Showalter-Hamilton, Brozovic or Jbara exists; those three labels rest on titles, abstracts and the
  2026-06-15 review note. The wanted-list rows are the remedy.
- `literature_check.py` is untouched in this step; the two proposals that need code are the label edits (sec. 2)
  and the Jbara citation correction.

## 7. Lead rulings (2026-10-11) and the record that results

Rulings after sec. 1-6 were sent:
- D1: option (b). Howett's topology label stays EMPTY. A label must come from the digest's text, and the digest left
  it empty on purpose. If someone later re-reads Howett and cites the section that describes repeated-moon
  cycler trajectories, (a) can be applied under its own pre-registration.
- D2: agreed, no new anchor for Ross and Roberts-Tsoukkas.
- D3: agreed, Showalter-Hamilton {resonant} to {ephemeris}. Stern gets working_bodies_scope "one" as proposed.
- The tadpole anchor's citation is a data error: correct its title and authors to the arXiv:2510.13479 record
  (Jbara, "Chaotic Dynamics and Zero-Velocity Structures in the Pluto-Charon CR3BP"), keep its label {binary-coorbital}.
  Recorded here on 2026-10-11. The anchor itself is edited only after rule U lands.
- The Brozovic anchor's provenance field marks the AJ 163:241 (2022) reference unverified.
- Wanted rows: sent by message to corpus-file-opus on 2026-10-11 (Showalter-Hamilton, Brozovic 2015,
  Brozovic-Jacobson 2024, Jbara, in that priority order).
- Rule U (twobody-gen2-opus): an untopologied signature can read "published" only through an anchor carrying the
  repeated-moon label; every other anchor's hit becomes inconclusive naming the anchor. Under rule U the
  untopologied Pluto case resolves without relabelling Howett: an empty-label Howett caps an untopologied Pluto
  signature at inconclusive, naming Howett, which is the honest reading until the digest supports a label.
- The anchor edits (Showalter-Hamilton, Stern scope, tadpole citation, Brozovic provenance) are applied ONLY after
  the lead says rule U's commit has landed; then C1-C3 and the gate's 26 controls and the #998 litcheck are rerun,
  everything byte-identical except C3 and the Howett citations.

Sec. 2 row 1 (Howett {repeated-moon}) and the D1 option (a) text of sec. 3 are therefore NOT adopted; the
expected outcomes of sec. 5 for C1 and C3 change as follows.

## 8. Expected gate outcomes for a #998-type signature, before and after rule U and the edits (check for gen2)

Signature: primary Pluto, sequence (Charon, Nix), period_k 4, V_inf (0.178, 0.130) km/s, working_bodies "one".
"Before" was run on the committed gate (f03cc22f, sec. 4); "after" is expected, not run.

| Signature | Before rule U (measured) | After rule U, Howett left empty (expected) |
|---|---|---|
| labelled {repeated-moon}, one working body | published, cites Howett | published, cites Howett (an empty-label anchor is unrestricted for a LABELLED signature; rule U concerns untopologied signatures) |
| UNLABELLED (no topology) | published, cites Howett | inconclusive, naming Howett (Howett carries no repeated-moon label) |
| labelled {pump-tour} at Pluto (C3) | published, cites Howett | published, cites Howett: UNCHANGED, because Howett's empty label is unrestricted; the pre-registered C3 outcome "not-found" (sec. 5) no longer applies under D1 (b) and is withdrawn |
| labelled {mga-tour}, sequence (Charon, Nix) | published, cites Howett (Stern also in footprint) | published, citing Howett or Stern by corpus order (unchanged in status) |
| labelled {resonant}, sequence (Nix, Hydra) | published (Howett and Showalter in footprint) | published via Howett (unchanged in status), no longer able to cite Showalter once it carries {ephemeris} |
| wrong system, primary Jupiter (C2) | known-architecture-new-system | unchanged |

C3 therefore needs a replacement negative that does not depend on Howett: use {pump-tour} at a system whose anchors
are all labelled (a negative that must stay "not published"), or the untopologied Pluto signature, which must
change from "published" to "inconclusive" under rule U. Both are for gen2 to run with rule U; none is run here.
Everything else of the gate's 26 controls and the #998 litcheck must be byte-identical after the anchor edits.
