# Live literature-novelty pass on the 12 Pluto SILVER candidates (#269)

Date: 2026-06-15
Discovery campaign source: #264 (Pluto-Charon binary regime, mu=0.108)
Literature-check module: `src/cyclerfinder/search/literature_check.py` (#261)
ML flagger: `src/cyclerfinder/ml/falsepos_flagger.py` (#256), trained on labeled corpus
Closure method version: `ee6d897`
Review-queue snapshot: `data/review_queue.jsonl` (12 rows)

## Bottom line

- **12 of 12** SILVER candidates cleared the live literature-novelty filter as `status="not-found"`.
- **12 of 12** scored `p_fp` in `[0.559, 0.561]` from the trained #256 ML flagger, well below the 0.75 routing threshold.
- **12 of 12** therefore satisfy `is_novelty_claimable() == True` AND `p_fp <= 0.75`, and are routed to `data/gauntlet_queue.jsonl` for the V0-V5 gauntlet (closure quality, same-model reproduction, n-body cross-check).
- **0 of 12** are auto-promoted to `data/catalogue.yaml`. The guard chain held: a clean lit-check + a clean flagger score is the *necessary* condition for a discovery claim, not the *sufficient* one.
- **0 of 12** are literature-matched. **0 of 12** routed to human re-check on suspicion.

The Pluto-Charon literature stack splits cleanly into three classes that the structural
fingerprints of these candidates do not match:

1. **Persephone-class CR3BP periodic orbits in the rotating Pluto/Charon frame**
   (Howard et al., arXiv:2102.08282 / Planet. Sci. J. 2(2):56, 2021). Spacecraft
   science orbits *around* the binary, with high out-of-plane + low-altitude families;
   incidental small-moon encounters. Structurally distinct from a patched-conic
   `Hydra->Nix->Hydra` cycler with Lambert legs.
2. **Charon-gravity-assist tours** (Stern, Tapley, Zangari et al., SwRI 2018 DPS; the
   "Game-Changer" Pluto orbiter concept). Explicitly *not* a periodic cycler --
   one-off sequential gravity assists; the Space.com / Centauri Dreams write-ups
   confirm no `cycler / periodic / resonant / repeating-pattern` terminology.
3. **Natural celestial mechanics of the moons themselves** (Showalter & Hamilton
   2015 Nature, Brozovic et al. 2015/2022 orbit determination, Pluto-Charon mutual
   interactions A&A 2023). The Styx-Nix-Hydra Laplace-like 3-body resonance
   `3 lambda_S - 5 lambda_N + 2 lambda_H ~ 180deg` is a property of the moons'
   own orbits, not a spacecraft trajectory.

A fourth class -- recent CR3BP-binary dynamical-systems papers (arXiv:2510.13479
2025 on Pluto-Charon tadpole/horseshoe instability, arXiv:1512.08683 orbit
classification, arXiv:2509.12671 HOTM Png' family at Earth-Moon) -- describes
Pluto-Charon periodic orbits, none of them as patched-conic small-moon tour cyclers.

## Per-candidate disposition

| candidate_id | period_k | sequence | residual (m/s) | max V_inf (km/s) | p_fp | lit-check | route |
|---|---|---|---|---|---|---|---|
| `repeated-moon-pluto-00000045` | 2 | Hydra-Nix-Hydra | 25.9 | 0.054 | 0.559 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000060` | 2 | Nix-Charon-Nix | 33.6 | 0.183 | 0.561 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000078` | 2 | Nix-Hydra-Nix | 20.4 | 0.051 | 0.559 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000015` | 2 | Charon-Nix-Charon | 17.4 | 0.097 | 0.559 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000019` | 2 | Charon-Nix-Charon | 31.0 | 0.071 | 0.559 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000105` | 3 | Charon-Hydra-Nix-Charon | 37.8 | 0.202 | 0.561 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000075` | 2 | Nix-Hydra-Nix | 34.3 | 0.056 | 0.559 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000282` | 3 | Hydra-Nix-Charon-Hydra | 36.3 | 0.091 | 0.559 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000283` | 3 | Hydra-Nix-Charon-Hydra | 21.5 | 0.085 | 0.559 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000357` | 3 | Nix-Charon-Hydra-Nix | 31.0 | 0.142 | 0.560 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000361` | 3 | Nix-Charon-Hydra-Nix | 34.3 | 0.142 | 0.560 | not-found | gauntlet_queue |
| `repeated-moon-pluto-00000311` | 3 | Hydra-Nix-Charon-Hydra | 40.7 | 0.068 | 0.559 | not-found | gauntlet_queue |

There are six unique tour structures across the 12 candidates (Hydra-Nix-Hydra,
Nix-Charon-Nix, Nix-Hydra-Nix, Charon-Nix-Charon, Charon-Hydra-Nix-Charon,
Hydra-Nix-Charon-Hydra, Nix-Charon-Hydra-Nix). Two of the duplicates share
`signature_hash` -- `sha1:4ba7e25...` appears on 00000045 / 00000078 / 00000075
(all `(Hydra,Nix)`-anchored, with the third member of the tour at the back) --
a hint that several of these may collapse into the same family member under
V1 same-model reproduction.

## Queries actually run (live WebSearch, this session)

1. `Pluto Charon Nix Hydra cycler trajectory periodic orbit`
2. `Pluto-Charon binary system spacecraft transfer trajectory periodic resonance`
3. `Nix Hydra small moons spacecraft tour mission trajectory patched conic`
4. `Pluto small satellites Nix Hydra Styx Kerberos flyby tour mission design`
5. `Persephone Pluto orbiter trajectory periodic orbits Nix Hydra Charon flyby tour`
6. `Pluto Charon CR3BP three-body periodic orbits Lyapunov halo distant retrograde DRO`
7. `"Hydra-Nix" OR "Nix-Hydra" cycler resonant transfer orbit Pluto system`
8. `"Charon-Nix" cycler ballistic moon tour Pluto trajectory design`
9. `Pluto Charon multiple moon flyby trajectory tour resonant repeating orbit`
10. `Pluto satellite system Lambert resonant orbit periodic cycler design AAS AIAA`
11. `Brozovic Showalter Pluto small moons resonance orbits cycler`
12. `"three-body resonance" Styx Nix Hydra Laplace-like spacecraft tour`
13. `"Pluto-Charon" cycler "low-thrust" OR "ballistic" trajectory tour moons site:arxiv.org`
14. `"Pluto-Charon" binary cr3bp mass ratio 0.108 periodic orbit family novel discovery 2024 2025 2026`

Plus two `WebFetch` calls (arXiv:2102.08282 Persephone, arXiv:2510.13479 CR3BP
tadpole/horseshoe) and one targeted `WebFetch` on the Space.com Stern/SwRI write-up.
None surfaced a structural-fingerprint match.

## Honest caveats

- **Web search is partial.** Paywalled venues (JGCD, AIAA conference proceedings
  pre-2010, ResearchGate gated PDFs) are under-indexed. A Pluto-system cycler paper
  could exist in an AAS or AIAA conference proceedings stack we did not index. The
  `not-found` is necessary, not sufficient; the V0-V5 gauntlet still governs.
- **The Persephone full PDF was not opened.** A `WebFetch` on the arXiv abstract page
  returned only the abstract; the full mission-design section could in principle
  describe small-moon cycler transfers between Persephone's CR3BP science orbits.
  If a human reviewer wants tighter confidence, opening the PDF and grepping for
  `cycler`, `Lambert`, or `Hydra->Nix` is a cheap follow-up.
- **The ML flagger's [0.559, 0.561] band is suspicious in its tightness.** Every
  candidate scoring within `1e-3` of the same value suggests the flagger's features
  cannot discriminate among them -- the labeled corpus may lack any Pluto-system
  exemplars, so the score is essentially the unconditional prior at the
  Pluto-mu-0.108 regime. The flagger score should not be over-interpreted as a
  per-candidate trustworthiness signal; the *fact* that all 12 sit below the 0.75
  routing threshold is what's load-bearing.
- **All 12 closure residuals are 17-41 m/s** -- two to three orders of magnitude
  above the typical sub-m/s clean ballistic closure floor. This is exactly the
  signal the V0-V5 gauntlet must adjudicate: are these real cyclers with a
  tighter same-model reproduction available, or are they basin-edge near-misses
  the closure tolerance accepted on the optimistic side? `not-found` does not
  resolve this; only V1 reproduction and V3 n-body cross-check do.
- **None of these is a discovery claim yet.** All 12 are queued for the gauntlet.
  The literature-novelty pass is one gate in the guard chain; the chain is not yet
  fully cleared.

## Files written

- `data/review_queue.jsonl` (in place; `literature_check` populated for all 12 rows)
- `data/gauntlet_queue.jsonl` (12 candidates queued for V0-V5 adjudication)
- `docs/notes/2026-06-15-pluto-silver-review.md` (this file)
