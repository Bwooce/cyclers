# Discovery campaign roadmap #283..#287 (post-#281/#282)

**Date:** 2026-06-16
**Frame:** discovery-program spec
(`docs/notes/2026-06-13-discovery-program-spec.md`). Three tracks:
A = richer genome, B = prioritizer, C = daemon.
**Already queued / running:** #281 multi-mu tulip Np=2 sweep across 14 moon
systems (RUNNING; 14 rows already in `data/tulip_sweep_281.jsonl`); #282
prioritized re-run of the #264 repeated-moon search with the full 5-tier
scorer stack (QUEUED).

This note proposes the next 5 campaigns. The default impulse is to chase
breadth (more moon systems, more body sets). The spec's yield rule pushes
back: highest-yield scouts have an existence prior, or a representable
topology our current genome covers but hasn't pointed at, or a cheap screen
that prunes the combinatorics. Five campaigns survive that filter; the rest
are documented as **do-not-do** with the reason so we don't re-litigate.

---

## Track-A capability snapshot (what's representable today)

Before recommending campaigns I checked the actual genome reach against
each candidate. Loaded findings (cite, do not guess):

* **Tulip genome (`src/cyclerfinder/genome/tulip.py`)** — Phase 5 operational
  (commit `54c69f0`). Direct-seed path `find_tulip_at_system(np_target=2)`
  lands Np=2 at any CR3BP system where the Koblick AMOSTECH seed corrects.
  Phase 4 multi-shooting (`src/cyclerfinder/genome/multi_shooting.py`)
  added the `n_segments` parameter. Higher-Np reproduction is *partial*:
  per the Phase 4 commit body (commit `42d60be`), Earth-Moon Np=3 closes
  via `n_segments=3` from the impactor-branch seed — but Koblick Table 4
  Np=3..7 have **negative `r_min_km`** (lunar impactors; `tulip.py:217-252`).
  Phase 4 file gates these as `xfail` honestly. Np ≥ 3 needs continuation
  OFF the impactor branch before the result is physically meaningful.
* **Repeated-moon multi-rev (`src/cyclerfinder/search/moon_cycler_genome.py`,
  `src/cyclerfinder/search/discovery_campaign.py`)** — Liang CGE class
  reproduced; daemon-driven; #264 swept Jupiter/Saturn/Uranus/Mars/Neptune
  empty, Pluto produced 12 SILVER. Body-agnostic via the
  `_registry_moons_for(primary)` helper.
* **Asymmetric (general) CR3BP corrector**
  (`src/cyclerfinder/search/cr3bp_general_periodic.py:329`,
  `correct_general_periodic`). Built for #249 C21 reproduction; only used
  in `tests/search/test_cr3bp_general_periodic.py` as a reproduction
  oracle. Not used in any novel-search driver yet.
* **Binary-star figure-seeded search**
  (`src/cyclerfinder/search/binary_star_search.py:143`, `figure_seeded_search`).
  Built (#255) — clean negative on the Roberts-Tsoukkas-Ross figures. Has the
  winding-(k1,k2) classifier and reaches mu up to 0.5.
* **5-tier Track-B prioritizer**: reachable-network (Braik-Ross),
  reachable-impulsive (Zhou-Armellin), resonance-network
  (Kumar/Rawat/Rosengren/Ross), FTLE (Canales-Howell), lobe-overlap
  (Hiraiwa-Bando), composed by
  `src/cyclerfinder/search/two_tier_prioritizer.py`. Tier-0 neural
  prefilter (`neural_reach_prefilter.py`) from Zhang-Topputo. The
  resonance-network is a scoring API (`ResonanceNetworkScorer.score_pair`)
  not a search driver.
* **BCR4BP / QBCP genome** — NOT representable. Spec digest
  (`docs/notes/2026-06-14-andreu-quasi-bicircular-digest.md`)
  recommends DEFER: no integrator, no Fourier series typed in, no goldens.
* **Sims-Flanagan low-thrust** — single-leg machinery exists
  (`src/cyclerfinder/core/sims_flanagan.py`, `search/lowthrust.py`) but
  there is no low-thrust *cycler* genome (closed periodic chain). The
  existing optimizer is single-leg heliocentric, not cycler-shaped.

---

## Campaign #283 — Higher-Np tulip continuation off the impactor branch (Earth-Moon Np=3)

1. **What is it?** Take the Earth-Moon Np=3 mathematical root that #268
   already closes from the Koblick Table-4 impactor IC (negative
   `r_min_km`, i.e. orbit passes through the lunar interior) and continue
   it in (x0 or z0) at fixed Np=3 topology until `r_min_km > 0`. Phase 4
   demonstrated the bifurcation root exists; the missing step is moving off
   the impactor sample point along the family curve to a physically valid
   member.
2. **Existence prior** — STRONG. Koblick 2023 AMOSTECH Table 4 publishes
   the full 15-petal family curve with explicit ICs for Np=1..15 (the IC
   is one sample point; the family is one-parameter). The impactor sign is
   incidental to the sample, not the family. Davis-Phillips-McCarthy
   (referenced in some surveys) is also believed to discuss higher-Np
   tulip-like families but no copy is on hand; treat the Koblick anchor
   as sufficient by itself.
3. **Representable?** YES with one small build. `tulip.py` has the
   topology classifier (`petal_count`), the multi-shooting corrector
   (`multi_shoot_periodic`), and the family-switch primitives
   (`family_switch.py`). What it lacks is a **continuation driver that
   holds Np constant** while sweeping the second free parameter of the
   one-parameter family. ~150-300 LOC: the natural-parameter wrapper
   already exists for Np=1->2 in `nrho_continuation.py`; the new code is
   a fixed-Np continuation in (z0, ydot0) that re-asserts the petal count
   via `petal_count` every step. <1 week including tests.
4. **Plausibly non-empty?** YES — Koblick's family curve is sourced;
   `r_min_km < 0` is local to the published IC.
5. **Cost** — minutes to hours per Np value. The continuation driver
   re-runs `multi_shoot_periodic` per step; the per-step cost is well
   under a second at Np=3 (Phase 4 evidence).
6. **Payoff** — 1-5 new physical Np=3 Earth-Moon tulips (and Np=4-6 if
   the continuation extends), each a clean SILVER candidate for the
   gauntlet. If they pass, they are the first non-Np=2 published tulips
   beyond the original Koblick table samples.
7. **Risk** — moderate. The family curve at Np=3 may stay below the
   lunar surface for the entire branch (rather than just at the sample),
   in which case the result is a clean negative ("Koblick's Np≥3 ICs are
   not the physical branch, and our continuation can't find one"). That
   is still a useful negative — it would close the higher-Np tulip
   question at Earth-Moon.
8. **Dependencies** — Track-A build (the fixed-Np continuation driver).
   No daemon needed; runs in-process.

---

## Campaign #284 — Asymmetric-corrector novel cycler scan at Earth-Moon (energy-sweep + winding-class enumeration)

1. **What is it?** Repurpose the `correct_general_periodic` asymmetric
   corrector (built for the #249 C21 reproduction at C_J=3.129389531...)
   as a *novel-search* driver: enumerate (k1, k2) winding classes 1..4 on
   each axis at a coarse grid of Jacobi C in [2.95, 3.20], shoot the
   corrector from random asymmetric ICs along the L1/L2 vertical line,
   classify each closed orbit by `winding_topology`, dedupe against the
   Ross-RT / Braik-Ross / Liang catalogue entries.
2. **Existence prior** — MEDIUM. Braik-Ross 2026 sweep their own
   Earth-Moon (k1, k2) grid and report 4 cycler members at one common
   energy plus the prior families. Asymmetric cyclers are explicitly
   what the C21 capability addition opened. There is no published
   guarantee of NEW members, but Ross's papers leave winding cells they
   never probed (k1>3, mixed signs, off-common-energy).
3. **Representable?** YES, no new code beyond a thin driver. The
   corrector exists, the winding classifier exists
   (`binary_star_search.winding_topology`), the literature check
   (`literature_check.py`) deduplicates against the catalogue and the
   pre-registered Braik-Ross / Ross-RT corpus. A ~200 LOC driver script
   wires the grid and emits SILVER candidates.
4. **Plausibly non-empty?** Cheap pre-screen: the Tisserand graph at
   Earth-Moon already validates accessibility, and the 5-tier prioritizer
   (#282) will score continuation neighbourhoods. If #282 surfaces any
   high-score (k1, k2) cell not yet represented, that cell is the first
   target for this scan.
5. **Cost** — fits in a day at 16 cores, asymmetric corrector wall-clock
   per shot ~ 0.5-2 s, 5000-shot grid manageable.
6. **Payoff** — 1-3 new Earth-Moon cycler members at unexplored winding
   classes. Each clears literature-novelty (#261) and routes to V0-V5.
7. **Risk** — winding-cell coverage may already be Braik-Ross's faithful
   negative ("they searched, they reported all the (k1, k2) ≤ (3, 2)
   members at C=3.1294"). If true, #284 returns a refined empty-region
   record. Even that is useful: it ratchets the empty-region registry
   onto the asymmetric corrector capability.
8. **Dependencies** — none beyond the existing corrector. Benefits from
   #282 completing first (Track-B output suggests targets); not blocked
   by it.

---

## Campaign #285 — Liang/CGE-style repeated-moon scan at Saturn (Titan-Rhea-Dione-Tethys) and Uranus regulars

1. **What is it?** Re-point the #254 repeated-moon multi-rev daemon
   genome at Saturn (Titan/Rhea/Dione/Tethys/Enceladus) and the five
   Uranian regulars (Miranda/Ariel/Umbriel/Titania/Oberon), using the
   5-tier prioritizer (#282 stack) to rank candidate sequences before
   spending closure compute. The point is that #264's Saturn+Uranus
   passes used the *unprioritized* enumeration; the 5-tier scorer is
   built; it has not yet been pointed at these moon systems.
2. **Existence prior** — MEDIUM. Liang reproduces CGE on Jupiter; the
   Takubo 2210.14996 tour-design paper exercises Titan-Enceladus
   resonance hopping (not as a cycler, but as evidence the V∞ graph
   admits chained transfers). Uranus has no direct cycler prior in our
   corpus but is a clean Tisserand-screening test (the screen
   self-prunes; if all V∞ graph edges are absent the campaign self-
   terminates cheaply).
3. **Representable?** YES, both the genome (#254) and the prioritizer
   (5-tier composition via `two_tier_prioritizer.py`, plus #267 manifold
   tier) are operational. No new code; spec is a config + a per-system
   `RepeatedMoonTarget` instantiation.
4. **Plausibly non-empty?** Cheap pre-screen: Tisserand admit/reject on
   each (primary, moon-pair, V∞) candidate, then prioritizer score.
   #264 hit only Pluto with 12 SILVER (review-gated); the unprioritized
   Saturn/Uranus sweeps returned empty. Re-running with the prioritizer
   targets *only* the high-score cells, so the per-cell budget is
   re-allocatable to deeper phasing sampling (current default 24 phase
   samples / max-rev 3).
5. **Cost** — ~1-2 days at 4 workers each (matching the #264 budget),
   restricted to high-score cells.
6. **Payoff** — best case: 1-3 new Saturn CGE-analogue triple cyclers
   (Titan/Rhea/Dione or Titan/Rhea/Enceladus). Realistic: a clean
   prioritized negative that closes the question at the current genome.
7. **Risk** — moderate. Saturn V∞ graph at the relevant Tisserand
   parameters may be sparse (Titan dominates the system; the inner moons
   sit at very different V∞ regimes than Titan); the prioritizer may
   honestly report no high-score cells. The campaign reduces to a
   ratchet on the empty-region registry. Mars/Neptune already swept
   empty in #264 unprioritized; the prioritizer may now rescue cells
   the unprioritized sweep missed.
8. **Dependencies** — #282 (prioritized #264 re-run) completing first
   to confirm the 5-tier-vs-unprioritized delta is real. If #282
   surfaces new SILVER at Jupiter/Mars/Neptune (vs the original empty
   verdicts), that re-validates the scorer composition before spending
   Saturn/Uranus budget.

---

## Campaign #286 — Direct fixed-mu binary-regime cycler scan at Pluto-Charon (mu = 0.108)

1. **What is it?** Re-purpose the `binary_star_search.figure_seeded_search`
   harness (built #255, applied as a clean negative on Roberts-Tsoukkas-Ross
   figure-read seeds) for direct fixed-mu *novel* search at Pluto-Charon's
   actual mu (~0.10876), enumerating winding (k1, k2) in 1..3 across a
   coarse Jacobi grid, with the corrector seeded from the L4/L5
   tadpole region rather than EM continuation. The motivation: at mu=0.108
   we are firmly in the binary regime (#255 found EM continuation does
   NOT reach it; the cyclers there are intrinsically binary-regime
   topology). The #281 tulip sweep already landed an Np=2 at Pluto-Charon
   (line 14 of `data/tulip_sweep_281.jsonl`, T=4.13 TU, J=3.048) — that's
   a *tulip*, NOT a Roberts-Tsoukkas multi-orbiter cycler. The two are
   different topological classes.
2. **Existence prior** — MEDIUM. Roberts-Tsoukkas-Ross 2026 figure-shows
   stable prograde cyclers at mu=0.1, 0.3, 0.5 (no numeric tables; the
   #255 figure-seeded search was a clean negative on those exact seeds).
   Pluto-Charon mu is published; the EM continuation explicitly failed
   to reach it (note: `2026-06-14-binary-star-mu-continuation-discovery.md`),
   so any cycler discovered here is a genuine new family member.
3. **Representable?** YES. The `figure_seeded_search` driver already
   handles arbitrary mu and winding classification. No new code; ~80
   LOC for a Pluto-specific seed-grid driver.
4. **Plausibly non-empty?** Tisserand isn't the right screen at mu=0.108
   (the patched-conic linearisation breaks); use the winding-class
   classifier directly. #275 trained the false-positive flagger on
   Pluto-class exemplars, so survivors get cheap routing.
5. **Cost** — fits in a day. The fixed-mu corrector is per-shot fast;
   the grid is small (Jacobi × winding).
6. **Payoff** — first true Pluto-Charon ballistic cycler (the existing
   12-row Pluto SILVER are repeated-moon Lambert tours from #264, a
   DIFFERENT topology). Direct existence-of-cycler-at-high-mu evidence.
7. **Risk** — Pluto-Charon may genuinely have no ballistic cyclers in
   the winding-class range we can search. mu=0.108 is in a regime where
   the libration-point geometry shifts (#252 showed the EM (1,1)/(3,1)
   members go linearly unstable before reaching such mu). The honest
   negative is still useful — it ratchets the binary-regime empty-region
   record off the figure-seeded sample onto a real grid scan.
8. **Dependencies** — none. Independent of #283/#284/#285.

---

## Campaign #287 — VEM heliocentric prioritized re-scan (5-tier scorer applied to the Earth-Mars-Venus space)

1. **What is it?** Apply the 5-tier prioritizer to Venus-Earth-Mars
   heliocentric trajectory enumeration, using the existing free-return
   chain and dsm_leg correctors as the closure layer. The motivation:
   the 5-tier scorer stack has NEVER been pointed at heliocentric space;
   #281/#282 are cislunar. The Jones AAS 17-577 VEM triple-cycler catalog
   sits in our literature corpus as a known anchor. The historical dense
   VEM scan (#110, 2816 points/row) used neither the prioritizer nor the
   asymmetric corrector and returned essentially empty (`OUTSTANDING.md`
   line ~835: EMEVVE outbound 0/2816, MEEVEM inbound 2/2816).
2. **Existence prior** — STRONG for reproduction (Jones VEM cyclers are
   published); UNCERTAIN for new discovery (Jones sweep is exhaustive in
   its parameter window).
3. **Representable?** Partial. The closure primitives exist (`free_return_chain`,
   `dsm_leg`, `multiarc_closure`). The 5-tier scorers are **EM-VU
   hardcoded**: `reachable_network.py:518` pins `VU_MS` to the Earth-Moon
   velocity unit; the impulsive and lobe scorers similarly assume EM
   normalisation. Repointing requires a scale-agnostic refactor: ~300-500
   LOC across all 5 scorers, primarily parameterising VU/LU/TU through.
   Multi-week build. **THIS IS A TRACK-A SHIM, NOT A SEARCH** — call it
   out honestly.
4. **Plausibly non-empty?** Heliocentric is the most populated cycler
   regime (Aldrin, S1L1, the Russell-Ocampo census of 200+ V0 rows).
   The space is sourced-dense; the question is whether new SILVER
   exists below the publication threshold.
5. **Cost** — 2-3 weeks build, then days of compute.
6. **Payoff** — potentially the largest novel-discovery yield of any
   campaign (Earth-Mars is where the cycler-application paying audience
   lives), and re-validation of the Jones catalog at full prioritizer
   strength. Worst case: re-confirms the published anchors.
7. **Risk** — the build cost is real, and the prioritizer's published
   validations are all cislunar. There's a non-trivial risk the EM-tuned
   scorers don't generalise (the Braik-Ross heading fan was validated
   only at C_J=3.1294 EM). Build deferred to *after* #283-#286 ship
   their results; if those underperform, the heliocentric pivot is the
   next natural move and the time spent rounding out cislunar discovery
   isn't lost.
8. **Dependencies** — Track-A scorer-scale shim (the only build item
   here, but a real one). Pairs naturally with #240 (KKT/surrogate
   amplifier) once Track-B has a heliocentric tier.

---

## Don't-do (and why)

* **BCR4BP / QBCP genome (campaign-spec candidate 5).** Andreu (1998) digest
  recommends DEFER; the QBCP is a *refinement* of CR3BP, not an opening
  of new families; the cost is full α_i Fourier series + integrator +
  frame, ≥3 weeks; the project has full DE440 ephemeris (V3) which sits
  strictly above QBCP in fidelity. The watch-item is unchanged: revisit
  only if cislunar candidates close in ephemeris but not CR3BP and the
  origin is opaque. No campaign now.
* **Davis-Phillips-McCarthy "Titan Np=6" hunt (candidate 2 in the
  prompt).** No copy of the paper is in the corpus; the reference
  appears to be a speculative ask. Without source the campaign violates
  the existence-prior yield rule. The substance — higher-Np tulips at
  Saturn-Titan — is partially covered by #281 (Np=2 done) and #283
  (Np≥3 continuation at Earth-Moon; Saturn-Titan can follow once the
  Earth-Moon path is proven). Re-list if the paper is sourced.
* **Halo-NRHO-DRO heteroclinic-connection campaign (candidate 9 in the
  prompt).** The #267 resonance-network module is a SCORER (a
  prioritizer tier), not a search driver. Turning it into a discovery
  driver is a Track-A capability project (search the network for
  heteroclinic chains with closure residual under budget), not a
  campaign. Tracked separately as a future capability item, not as #283-7.
* **Low-thrust cycler genome (candidate 10).** The single-leg
  Sims-Flanagan machinery exists, but a *cycler* (closed periodic
  chain, multiple flybys, repeating) needs a multi-leg chain
  formulation, low-thrust flyby-continuity, and the M5-style optimizer
  wired into a closure metric. Estimated 1500-2500 LOC + several weeks
  (large enough to be its own quarter). The spec's
  recommendation matches: defer until A+B discovery cyclers are found.
* **Asymmetric corrector at Saturn/Uranus moon systems.** The
  asymmetric corrector is at Earth-Moon mu; porting it is mostly
  re-normalising constants but adds risk without an existence prior at
  those systems. Defer until #284 either succeeds (then port) or fails
  cleanly (don't bother).
* **More moon-system tulip scans.** #281 is sweeping 14 systems and 4
  already converged direct-seed. Extending the body list (Mars-Deimos,
  Saturn-Mimas/Tethys/Dione, Jupiter-Amalthea, Neptune-Proteus, Pluto-
  Nix/Hydra) is mechanically possible but the satellites GMs there are
  ≤10 km^3/s^2 — the registry already notes the Tisserand screen
  self-prunes them. No payoff.

---

## Summary table

| Rank | Task # | Name | Cost | Payoff | Risk | Dependency |
|------|--------|------|------|--------|------|------------|
| 1 | #283 | Higher-Np tulip continuation off impactor branch (EM Np=3) | <1 week build + hours compute | 1-5 new physical Np≥3 EM tulips | family curve may stay sub-lunar | Track-A: fixed-Np continuation driver (~200 LOC) |
| 2 | #284 | Asymmetric-corrector novel cycler scan at EM (winding × Jacobi grid) | fits in a day | 1-3 new asymmetric EM cyclers | Braik-Ross may already cover the cells | none beyond existing corrector; benefits from #282 |
| 3 | #285 | Prioritized repeated-moon scan at Saturn/Uranus regulars | 1-2 days compute | 1-3 new Saturn triple cyclers | Tisserand screen may self-prune | #282 closes first (validates prioritizer delta) |
| 4 | #286 | Direct fixed-mu binary-regime scan at Pluto-Charon (mu=0.108) | fits in a day | first Pluto-Charon ballistic cycler | mu=0.108 may have none in winding range | none |
| 5 | #287 | VEM heliocentric prioritized re-scan (5-tier ported) | 2-3 weeks build + days compute | new VEM cyclers / Jones re-validation | scorer port risk; published anchors may dominate | Track-A: scale-agnostic scorer port |

Rank 1-4 are scout-and-search (genome / driver work measured in days).
Rank 5 is the largest payoff opportunity but is a multi-week build and
should follow at least #284 succeeding (proves the asymmetric corrector
finds novel cells) before the heliocentric port is funded.

---

## Observations on #281/#282 priority

* **#281 already covers more than was advertised.** 14 rows in
  `data/tulip_sweep_281.jsonl`; 5 converged direct-seed (EM, Jupiter-
  Ganymede, Saturn-Titan, Neptune-Triton, Pluto-Charon). The "easy"
  multi-mu tulip space is essentially mapped. Continuing #281 to
  Np=3..6 at the 5 systems where Np=2 converged is mechanically the
  next step — but that IS campaign #283, sourced at Earth-Moon first.
  **Recommendation:** treat #281 as effectively done at Np=2;
  immediately fold its 5-system Np=2 success into the gauntlet queue
  (literature-novelty → V0-V5) rather than letting it run further.
* **#282 (prioritized #264 re-run) is correctly second priority.** It
  validates the 5-tier prioritizer delta before #285 spends Saturn/Uranus
  budget on the assumption the prioritizer helps. Do not start #285
  before #282 returns its first prioritizer-vs-unprioritized comparison.
* **No re-shuffle of #281/#282 themselves.** The order is right; the
  observation is just that #281's tail (further moon systems beyond the
  14) is already returning fast-fail `seed_no_converge` and is unlikely
  to surface anything before the next campaign.
