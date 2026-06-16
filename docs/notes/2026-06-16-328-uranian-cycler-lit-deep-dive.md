# #328 — Uranian cycler literature deep-dive

**Date**: 2026-06-16
**Task**: #328 (wider literature pass for #327 verified SILVER)
**Scope**: WebSearch + WebFetch over published Uranian dynamics / mission-design literature; triage hits against the verified SILVER candidate; expand `KNOWN_CORPUS` deterministically with anchors uncovered.
**Discipline**: research-only — no catalogue writeback, no novelty claim. A `clean literature-fresh` verdict here is NECESSARY-NOT-SUFFICIENT for novelty; the #306 3D V0-V5 gauntlet still has to pass.

## The candidate being checked

Per `#327` verification (commit `3500c71`):

- **Primary**: Uranus
- **Tour**: Umbriel → Oberon → Umbriel (a (k1,k2) = (1,1) repeated-encounter cycler — one Umbriel-side leg, one Oberon-side leg per cycle)
- **V∞ tuple (km/s)**: (0.92, 0.96, 0.89) — comparable across the three legs (ballistic-low regime by Uranian standards)
- **ToF**: ~14.94 d/leg
- **Bend**: max 14.7° / 39° (within ballistic-flyby limits)
- **Closure residual**: 0.025 km/s (passes SILVER gate)
- **DOP853 cross-check**: 2.7e-11 nondim residual
- **Offline `literature_check.py`**: not-found (KNOWN_CORPUS had no Uranian anchors)
- **ML p_fp**: 0.591

This task is the wider published-literature pass per `feedback_literature_novelty_check_baseline`: offline not-found is necessary-not-sufficient.

## Search methodology

Strategy: **most-specific (the exact (Umbriel, Oberon) (1,1) cycler topology) → moon-pair adjacencies (Titania-Oberon, Ariel-Umbriel) → general Uranian moon-tour mission design → general Uranian-system CR3BP/CR4BP dynamics → historical (Voyager 2)**.

WebSearch was primary; WebFetch was attempted for arXiv abstract pages and `archive.org`/`ntrs.nasa.gov` records. AIAA `arc.aiaa.org` URLs 403'd as expected; institutional repositories (Purdue, ResearchGate, ScienceDirect, arXiv html mirrors) and ADS/IOPscience were tried instead.

### Query trail (in execution order)

1. `"Uranus" cycler trajectory repeated encounter` — established the gap (no dedicated Uranus cycler trajectory in the surfaced corpus).
2. `Titania Oberon resonance orbital dynamics` — surfaced the candidate's moon pair in current resonance / dynamics literature.
3. `Uranus mission trajectory design satellite tour gravity assist` — surfaced UOP / decadal mission corpus.
4. `Uranus orbiter probe Decadal mission trajectory Titania Oberon tour` — confirmed UOP tour design literature.
5. `"Trajectory Options for a Uranus Orbiter and Probe" satellite tour` — 2024 Acta Astronautica paper.
6. `Umbriel Oberon Titania flyby tour V-infinity moon` — no direct cycler hit (only dynamics / mission-design surveys).
7. `"polar Uranus orbiter" satellite tour Heaton Strange Longuski` — surfaced Heaton-Longuski 2003 JSR (the canonical Galileo-style Uranian tour anchor).
8. `Uranus moon tour Tisserand graph V-infinity leveraging satellite` — confirmed Tisserand/VILT generalises to Uranus.
9. `Heaton Longuski "Galileo-style" Uranus tour Journal Spacecraft Rockets` — pinned the JSR Vol.40 No.4 (2003) reference.
10. `"Uranian system" CR3BP periodic orbit Umbriel Titania Oberon three-body` — surfaced Kumar's Uranus-Oberon PCR3BP / Uranus-Titania-Oberon CCR4BP work.
11. `Kumar "Uranus" Oberon resonant invariant manifold heteroclinic periodic` — confirmed Kumar 2025 multi-shooting paper covers Uranus-Oberon 3:4/4:5/5:6 exterior + 4:3/5:4/6:5 interior MMRs with heteroclinic transitions.
12. `Uranian "Lambert problem" intermoon transfer cycler "Umbriel" "Oberon" trajectory` — surfaced Canales/Howell MMAT method (arXiv:2110.03683) applied to Titania-Oberon halo-to-halo transfer.
13. `"Canales" "Howell" Uranian system 2BP-CR3BP patched moon transfer` — confirmed MMAT 2BP-CR3BP patching applied to Uranus moons.
14. `"Mutation Operator for Resonance Flybys" "Moon Tour Design" Uranus authors` — Choi/Abdelkhalik/He JGCD 2024 (Europa-Clipper-like / Saturnian, not Uranian; noise).
15. `"Uranus" 1:1 co-orbital trojan satellite trajectory dynamics` — Trojan dynamics (not cycler).
16. `Voyager 2 Uranus trajectory Miranda Ariel flyby gravity assist 1986` — historical reference, single flyby.
17. `"QUEST" Uranus orbiter New Frontiers tour design satellite` — Jarmak et al. Acta Astronautica 2020 mission concept (polar orbit, no satellite tour).
18. `Russell Buffington Strange "V-infinity globe" Uranus inclination reduction` — confirmed V∞-globe methodology applied to Uranus (inclination-reduction sequences in mission-design context).
19. `arxiv "Uranian" tour design ballistic 2024 2025 cycler trajectory Russell` — no novel cycler hit; surfaced 2025 Flagship community-input poll (arXiv:2505.05514).
20. `"Uranian satellite" 1:1 resonance cycler orbit ballistic moon-moon transfer` — explicit `1:1 cycler` query at Uranus: no hit. Confirmed (negative-result) that the candidate's topology is not directly in the surfaced record.

## Per-paper triage

| # | Paper | Classification | Notes |
|---|---|---|---|
| 1 | Heaton & Longuski, **"The Feasibility of a Galileo-Style Tour of the Uranian Satellites,"** *J. Spacecraft & Rockets* 40(4):591-596 (2003), DOI 10.2514/2.3981 (NTRS 20020021945; AIAA 2001-3859) | **TOUR-TYPE** | 811-day three-phase Galileo-style one-shot ballistic tour with 40+ flybys including Titania & Oberon. STOUR (patched-conic). Foundational anchor for any Uranian moon-tour work but NOT a cycler (no repeated periodic pattern; designed for end-of-mission insertion). Tour pair-equivalence noted: "Titania-Oberon ≈ Ganymede-Callisto; Ariel-Umbriel ≈ Io-Europa." |
| 2 | Sims, Finlayson, Rinderle, Vavrina & Kawalkowski et al., **"Conceptual mission design of a polar Uranus orbiter and satellite tour"** (Acta Astronautica / IEEE Aerospace 2014, ResearchGate 265491803) | **MISSION CONCEPT** / TOUR-TYPE | Baseline 424-day, 619 m/s tour with TWO targeted flybys each of the five major moons. Two-flyby-per-moon is the (1,1) count by VISITS but it is a one-shot insertion tour (Hohmann-like inter-moon legs after orbit-insertion), NOT a periodic ballistic cycler. |
| 3 | Jarmak, Brinckerhoff et al., **"QUEST: A New Frontiers Uranus orbiter mission concept study,"** *Acta Astronautica* 170:6-26 (2020), DOI 10.1016/j.actaastro.2020.01.030 (ADS 2020AcAau.170....6J) | **MISSION CONCEPT** | Polar orbit; one-year primary; no satellite tour (focuses on Uranus interior / magnetosphere). NOT cycler-relevant. |
| 4 | (Uranus Orbiter & Probe / UOP teams) **"Trajectory Options for a Uranus Orbiter and Probe"** (2024 ResearchGate 386106230) | **TOUR-TYPE** | UOP trade-space study. Two targeted flybys per major moon over multi-year science phase. Inclination-reduction sequence via repeated Titania flybys. Inclination-reduction phase IS effectively (1,k) at Titania, but it's a planned insertion-and-reduction sequence — bound to the post-orbit-insertion energy regime, NOT a free-floating periodic cycler. |
| 5 | UOP Decadal mission concept design study (2023) — Cohen, Simon et al. (NASA GSFC / APL); Aerocapture variant Saikia et al. *Acta Astronautica* (2022) DOI 10.1016/j.actaastro.2022.10.026 (ScienceDirect S0094576522005422) | **MISSION CONCEPT** | Flagship decadal architecture; tour sub-detail varies but always insertion-based, not periodic. |
| 6 | UOP "Mission Challenges and Concept Updates Since the OWL Decadal Survey" *Planet. Sci. J.* 6:ae680c (2025), DOI 10.3847/PSJ/ae680c | **MISSION CONCEPT** | Confirms repeated Titania flybys for inclination reduction → equatorial tour of all five moons in 4.5 yr science phase. Same insertion-tour topology. |
| 7 | Flagship Science-Driven Tour Design Community Input Poll (Simon et al., arXiv:2505.05514, LPSC 1207, 2025) | **MISSION CONCEPT** | Community poll on tour priorities; no specific cycler / repeated-encounter topology beyond UOP's 2/3-flyby-per-moon baseline. |
| 8 | Kumar, **"Multi-shooting parameterization methods for invariant manifolds and heteroclinics of 2-DOF Hamiltonian Poincaré maps, with applications to celestial resonant dynamics,"** arXiv:2509.03655 (2025, math.DS) | **STRUCTURAL ADJACENT** (resonant orbit dynamics, NOT cycler-class) | **Section 6.2: Study of MMR overlap in the Uranus-Oberon PCRTBP.** Computes 3:4, 4:5, 5:6 exterior + 4:3, 5:4, 6:5 interior MMR unstable periodic orbits at OBERON, plus stable/unstable manifolds and heteroclinic transitions between them at varying Jacobi energy. The (k1,k2) = (1,1) Umbriel-Oberon two-body cycler is NOT in scope — Kumar's work is single-moon (Oberon) interior/exterior MMRs at SC-Uranus-Oberon PCR3BP. STRUCTURAL ADJACENT: same primary, related dynamics framework, no overlap with the candidate's two-moon-pair (1,1) topology. |
| 9 | Kumar (companion work), **"4th Body-Induced Secondary Resonance Overlapping Inside Unstable Resonant Orbit Families: a Jupiter-Ganymede 4:3 + Europa Case Study,"** arXiv:2309.06073; extensions to Uranus-Titania-Oberon CCR4BP (6:5 Oberon MMR + Titania perturbation, secondary resonances at 25/69, 21/58, 17/47, 30/83, 13/36, 22/61, 9/25) | **STRUCTURAL ADJACENT** | Jupiter case study is the demonstration system; Uranus-Titania-Oberon CCR4BP extension appears in companion search excerpts (Kumar's broader corpus). The CCR4BP secondary resonance ratios are **inside-the-6:5-Oberon-family** sub-harmonics (9/25, 13/36, etc.) — these are TORUS BIFURCATIONS in the perturbed 4-body problem, NOT moon-to-moon (1,1) repeated-encounter cycler topology. STRUCTURAL ADJACENT, not direct. |
| 10 | Canales, Howell, Fantino, **"Transfer design between neighborhoods of planetary moons in the circular restricted three-body problem: the moon-to-moon analytical transfer method,"** *Celestial Mechanics and Dynamical Astronomy* 133:36 (2021), ADS 2021CeMDA.133...36C; arXiv:2110.03683 | **STRUCTURAL ADJACENT** (moon-to-moon transfer but ONE-SHOT halo-to-halo, NOT cycler) | Case study: spacecraft from L2 northern halo orbit in Uranus-Titania system → L1 southern halo orbit in Uranus-Oberon system via unstable manifold from Titania halo + stable manifold into Oberon halo. **2BP-CR3BP patched analytical method.** Same moons as the candidate's neighborhood (Titania/Oberon), but it's a one-shot manifold-mediated transfer between halo orbits — NOT a repeated-encounter cycler at Umbriel-Oberon. STRUCTURAL ADJACENT. |
| 11 | Canales, Howell et al., AAS 21-625 / AAS 21-234 "Using Finite-Time Lyapunov Exponent Maps for moon transfer design"; *Transfers Between Moons with Escape and Capture Patterns via Lyapunov Exponent Maps* (JGCD 2023, arXiv:2308.10029) | **STRUCTURAL ADJACENT** | Same MMAT framework, extends transfer modes (capture / escape / landing) with FTLE maps. Same one-shot-transfer-not-cycler classification. |
| 12 | Choi, Abdelkhalik & He, **"Mutation Operator for Resonance Flybys in Moon Tour Design Optimization,"** *JGCD* 47(11):2287 (2024), DOI 10.2514/1.G007708 | **NO MATCH** | Methodology paper; applies to Europa Clipper / Saturnian system. Not Uranus. |
| 13 | Various Uranian moon dynamics (resonance history / tidal evolution): Crida-Charnoz 2020 arXiv:2005.12887, Caldas et al. arXiv:2509.24631 (Ariel-Umbriel 2:1 capture/escape), arXiv:2403.17896 / 2403.17897 (Ariel-Umbriel 5:3 MMR evolution), Cuk et al. (Quillen et al.) MNRAS 445:3959 (resonant chains in inner moons) | **DYNAMICS** | Long-term tidal-resonance evolution of the Uranian moon system. Uranian moons currently NOT in active mean motion resonance with each other (Cuk: contrast to Jovian Laplace). NOT spacecraft cycler design. |
| 14 | Ribeiro-de-Sousa et al. (and others) **"Mapping Long-Term Natural Orbits about Titania,"** arXiv:2203.14445 (2022) | **DYNAMICS** | Long-term orbit-stability around Titania (orbiter context). Not cycler. |
| 15 | Strange, Russell & Buffington, **"Mapping the V-infinity globe,"** AAS 07-277 (2007) | **DYNAMICS / METHODOLOGY** | The V∞-globe formalism that Heaton-Longuski and UOP teams apply to Uranian inclination-reduction sequences. Already in KNOWN_CORPUS (Strange-Campagnola-Russell anchor + Strange-Russell pump-tour anchor). |
| 16 | de la Fuente Marcos, **"Systematic survey of the dynamics of Uranus Trojans,"** A&A 632:A19 (2020), arXiv:1912.10273; **"Asteroid 2014 YX49: a large transient Trojan of Uranus,"** MNRAS 467:1561 (2017) | **DYNAMICS** | 1:1 mean motion resonance for Trojans at L4/L5 around Uranus — tadpole regime only (no horseshoes). Natural-body dynamics; NOT a spacecraft cycler. |
| 17 | Voyager 2 historical: Stone-Miner JGR 92:14873 (1987); various review articles | **DYNAMICS / HISTORICAL** | Single 1986 flyby. Not cycler-relevant. |

### Triage tally

- **DIRECT MATCH** (paper publishes the Umbriel-Oberon (1,1) cycler family): **0**
- **STRUCTURAL ADJACENT** (Uranian primary + related dynamics or moons): **4** (Kumar 2509.03655 Uranus-Oberon PCRTBP; Kumar/companion Uranus-Titania-Oberon CCR4BP; Canales-Howell 2110.03683 Titania-Oberon halo-halo MMAT; Canales-Howell-Fantino 2308.10029 escape/capture)
- **TOUR-TYPE** (Uranus moon TOUR, one-shot): **4** (Heaton-Longuski 2003; Sims et al. 2014; UOP Trajectory Options 2024; UOP Mission Challenges 2025)
- **DYNAMICS** (resonance/tidal evolution/Trojan/orbiter-stability): **6** (Trojans; Caldas; Crida-Charnoz; Cuk-Quillen MNRAS; long-term Titania orbits; obliquity/libration)
- **MISSION CONCEPT** (Decadal whitepaper class): **3** (QUEST; UOP Decadal; UOP Aerocapture; Flagship community poll)
- **NO MATCH**: rest

## Cross-check of the candidate against the closest adjacencies

### Kumar (arXiv:2509.03655) — Uranus-Oberon PCRTBP

Kumar's resonant periodic orbits at Oberon are **SC-Uranus-Oberon CR3BP family** — the spacecraft sits in a 3:4 / 4:5 / 5:6 / 4:3 / 5:4 / 6:5 mean motion resonance with **Oberon's orbital period**, periodically encountering Oberon at one (or few) longitudes per cycle. The candidate is a two-moon (Umbriel ↔ Oberon) repeated-encounter Lambert-style cycler — fundamentally different topology. Per-encounter V∞ at Oberon in Kumar's family is set by Jacobi constant ∈ specific range; we don't have the exact numeric to cross-check against the candidate's V∞_O = 0.96 km/s without pulling figures from the full PDF. But the **topology mismatch is decisive**: single-moon-period MMR ≠ moon-to-moon (k1,k2) cycler.

### Canales-Howell (arXiv:2110.03683) — Titania-Oberon halo manifold transfer

This is the closest moon-pair adjacency. But: (a) the moons are Titania-Oberon, not Umbriel-Oberon; (b) the trajectory is a ONE-SHOT manifold-mediated transfer between Uranus-Titania-L2 halo and Uranus-Oberon-L1 halo; (c) it requires the spacecraft to insert into a halo orbit at each end — no free-flight repeated flyby. Different topology. Not the candidate.

### Sims et al. (2014) / UOP Trajectory Options (2024) — orbiter tour with 2 flybys per moon

The "two flybys of each of Titania, Oberon, Umbriel, Ariel, Miranda" requirement looks superficially like the candidate's (1,1) pattern, but: (a) the orbiter is Uranus-centered (Uranus is the primary), so each "flyby" is from a planet-centered orbit, not a free-flight ballistic chain that closes back on itself; (b) the V∞ regime is set by post-orbit-insertion energetics (typically < 0.5 km/s at the moons after inclination-reduction; cf. the 619 m/s tour ΔV); (c) it's an open chain, not a periodic cycler. No direct match.

## Verdict

**CLEAN LITERATURE-FRESH** — no published paper in the surfaced corpus reports the Umbriel-Oberon (1,1) repeated-encounter ballistic cycler at V∞ ≈ 0.9 km/s as a periodic free-flight trajectory. The closest adjacencies are:

- **STRUCTURAL ADJACENT (single-moon MMR)**: Kumar 2025 (arXiv:2509.03655) — covers Uranus-Oberon PCR3BP resonant orbits at Oberon's period; complementary dynamics framework but NOT moon-pair cycler.
- **STRUCTURAL ADJACENT (moon-pair one-shot transfer)**: Canales-Howell 2021 (arXiv:2110.03683) — covers Titania-Oberon halo-to-halo manifold transfer; different moons (Titania vs Umbriel), one-shot not cyclic.
- **TOUR-TYPE (one-shot)**: Heaton-Longuski 2003 JSR + Sims 2014 + UOP 2024 — one-shot insertion tours with multiple per-moon flybys; not periodic free-flight cyclers.

Discipline reminder: this verdict is **necessary-not-sufficient** for novelty. The #306 3D V0-V5 gauntlet (parallel agent) and the human reviewer still govern. A "literature-fresh" candidate at SILVER closure with a non-trivial ML p_fp (0.591) and a topology that has demonstrable structural neighbors in the literature is still subject to V0-V5 invariants.

### Confidence in the verdict

- WebSearch breadth: 20 queries covering targeted (Umbriel-Oberon), adjacent (Titania-Oberon, Ariel-Umbriel), general (moon-tour, Tisserand, V∞-globe), and natural-dynamics (Trojan, resonance) angles.
- WebFetch coverage: arXiv abstracts loaded cleanly for the structural-adjacent papers (Kumar 2509.03655, Canales-Howell 2110.03683). AIAA `arc.aiaa.org` 403'd as expected, but the JSR / Acta Astronautica abstracts and key technical claims were recovered via search-result snippets and institutional repository abstracts (Purdue, NTRS, ResearchGate, ADS).
- Critical paywalled venue: full text of Heaton-Longuski 2003 JSR (DOI 10.2514/2.3981) only via arc.aiaa.org (403) — claims about tour structure were recovered from secondary sources (NTRS 20020021945 record + spaceflight history blogs that cite the paper). I am confident from those secondary citations that the paper is a Galileo-style one-shot tour (40+ flybys in 811 days from a single Uranus-arrival), NOT a cycler. **Phase 2 recommendation**: a human reviewer with AIAA member access should pull the JSR PDF to verify (a) no late-mission cycler-class trajectory is described, and (b) per-leg V∞ values for the Titania/Oberon legs to check whether the candidate sits in the same regime.

## Phase 2 recommendation

The verdict is **clean literature-fresh** but with one identified paywalled gap:

- **Heaton & Longuski JSR 40(4):591-596 (2003) DOI 10.2514/2.3981** — full PDF on AIAA arc.aiaa.org (403 for me). Human access (AIAA membership) recommended to confirm the candidate's Umbriel-Oberon (1,1) topology is not buried in a late-mission tour phase. Analogous to the #279 Finley disposition framework.

If that PDF (and the Sims 2014 paper full PDF, also typically paywalled) confirm tour-only / no cycler-class repeated-encounter trajectory at Umbriel-Oberon, the verdict upgrades to **strong literature-fresh**. If either paper describes a quasi-cycler Umbriel-Oberon repeated-encounter sub-pattern within the tour, the verdict downgrades to **REDISCOVERY** (in which case re-classify as `precursor_mga` per the #294 four-class taxonomy).

## KNOWN_CORPUS expansion

Six new anchors added to `src/cyclerfinder/search/literature_check.py::KNOWN_CORPUS` covering the Uranian dynamics + mission-design literature uncovered here. Each cites by DOI / arXiv / public venue per `feedback_golden_tests_sourced_only`.

1. **Heaton-Longuski 2003 Galileo-style Uranian tour** (JSR 40(4):591-596, DOI 10.2514/2.3981) — TOUR-TYPE foundational anchor.
2. **Sims et al. polar Uranus orbiter + satellite tour** (2014, ResearchGate 265491803) — TOUR-TYPE conceptual mission design.
3. **Kumar Uranus-Oberon PCRTBP MMR study** (arXiv:2509.03655) — STRUCTURAL ADJACENT single-moon MMR dynamics.
4. **Canales-Howell MMAT method (Titania-Oberon case study)** (Celest. Mech. Dyn. Astron. 133:36, arXiv:2110.03683) — STRUCTURAL ADJACENT moon-pair one-shot transfer.
5. **Jarmak QUEST Uranus orbiter concept** (Acta Astron. 170:6, DOI 10.1016/j.actaastro.2020.01.030) — MISSION CONCEPT.
6. **UOP Decadal Mission Concept + Trajectory Options** (UOP 2024 + Saikia aerocapture 2022 + UOP 2025 PSJ DOI 10.3847/PSJ/ae680c) — MISSION CONCEPT / TOUR-TYPE current-Decadal anchor.

These let the matcher surface a published-anchor warning whenever a future candidate's structural fingerprint overlaps (primary = Uranus + body_set ⊆ {Miranda, Ariel, Umbriel, Titania, Oberon}).

---

## Files written

- This note: `/home/bruce/dev/cyclers/docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md`
- `KNOWN_CORPUS` expansion: `/home/bruce/dev/cyclers/src/cyclerfinder/search/literature_check.py` (six new `CorpusAnchor` entries at the end of the tuple)

## Sources cited

- Heaton & Longuski 2003 JSR DOI 10.2514/2.3981 (NTRS 20020021945; archive.org/details/nasa_techdoc_20020021945)
- Sims et al. 2014 conceptual polar Uranus orbiter tour (ResearchGate 265491803)
- Jarmak et al. 2020 QUEST Acta Astronautica DOI 10.1016/j.actaastro.2020.01.030
- UOP "Trajectory Options" 2024 (ResearchGate 386106230)
- UOP "Mission Challenges and Concept Updates Since the OWL Decadal" 2025 PSJ DOI 10.3847/PSJ/ae680c
- UOP Aerocapture (Saikia et al. 2022) Acta Astronautica DOI 10.1016/j.actaastro.2022.10.026
- Simon et al. Flagship Science-Driven Tour Design 2025 (arXiv:2505.05514, LPSC 1207)
- Kumar 2025 "Multi-shooting parameterization methods..." (arXiv:2509.03655)
- Kumar / Kumar-Anderson-Jorba "4th Body-Induced Secondary Resonance Overlapping..." Jupiter case study (arXiv:2309.06073) + Uranian extensions
- Canales, Howell & Fantino 2021 MMAT "Transfer design between neighborhoods of planetary moons in the CR3BP" Celest. Mech. Dyn. Astron. 133:36 (arXiv:2110.03683)
- Canales-Howell-Fantino 2023 JGCD "Transfers Between Moons with Escape and Capture Patterns via Lyapunov Exponent Maps" (arXiv:2308.10029)
- Choi-Abdelkhalik-He JGCD 2024 DOI 10.2514/1.G007708 (Europa Clipper / Saturnian; noise)
- Strange, Russell & Buffington 2007 "Mapping the V-infinity globe" AAS 07-277
- Voyager 2 Uranus encounter: Stone-Miner JGR 92:14873 (1987)
- Cuk et al. "Resonant chains and three-body resonances in the closely packed inner Uranian satellite system" MNRAS 445:3959
- de la Fuente Marcos "Systematic survey of the dynamics of Uranus Trojans" A&A 632:A19 (arXiv:1912.10273)
- Caldas et al. arXiv:2509.24631 (Ariel-Umbriel 2:1)
- Crida-Charnoz et al. arXiv:2005.12887 (Uranian dynamical history)
- Ribeiro-de-Sousa et al. arXiv:2203.14445 (Long-term natural orbits about Titania)
