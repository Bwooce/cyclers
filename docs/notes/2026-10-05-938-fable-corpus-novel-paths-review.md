# #938 Corpus review: untried routes to novel cycler and quasi-cycler orbits

**Date:** 2026-10-05 (AEDT). **Task:** `#938` (analysis only; no code, no catalogue edit, no dispatch).
**Author:** corpus-review-fable (Fable, sections 0-4 first draft, stalled 17:26 AEDT) and
corpus-review-fable-2 (Fable, from 19:00 AEDT: sections 3a, 5, 6, the matrix CSV, the source checks
and the corrections marked "[F2]"), with six Sonnet extract readers, four Sonnet technique-inventory
readers and one Sonnet source-check reader. Companion file:
`docs/notes/2026-10-05-938-technique-case-matrix.csv` (32 techniques x 16 cases).
**Question:** which routes to NOVEL cycler-class orbits are (a) grounded in a held paper, (b) not
already tried or on the `#864` do-not-do list, (c) runnable with the present code or a bounded build,
and (d) likely to produce something outside the published record.

## 0. How to read this note

- Evidence tags. `[P p.N]` a reader opened the PDF at page N (PDF page unless "journal p."). `[D]`
  read in a project digest under `docs/notes/`, not at the source. `[C]` my own small calculation.
  `[I]` my inference. Every quotation below carries one of these.
- Papers are cited by author, year and title, with the corpus path as `cyclers_pdf/papers/<file>`.
  The corpus is a private repository and is never named here.
- Policy labels follow `docs/spec.md` 16.4 and `docs/notes/2026-09-07-875-novelty-policy-decision.md`:
  (i) author-excluded class, (ii) known architecture at a never-treated body set, (iii) theorem-generic
  object not in a published family. "Literal collision first" applies to every route.
- Costs are GUESSES unless a measured unit cost is cited (`#864` rule 8). Probabilities are my judgment.
- Where a route touches a system the `#864` sec. 8 list says not to re-sweep for novelty (Earth-Moon,
  Earth-Mars, Sun-planet CR3BP), the entry says so. The lead decides.

## 1. Method and a finding about the corpus itself

Six Sonnet readers were given disjoint themes (moon systems and tours; heliocentric cyclers and
global search; Earth-Moon cyclers, resonance networks and Sun-forced continuation; four-body models,
tori and connections; CR3BP and ER3BP periodic-orbit theory; methods, flown tours and textbooks).
Each read the project digest first, then the PDF's abstract, conclusions and future-work pages, and
returned per-paper extracts of: method, body sets and parameter ranges actually covered, quoted scope
exclusions and future work, and existence claims made without computation. Their extracts (about
450 kB) are in the session scratch area, not the repository. I read all six in full, cross-checked
the candidate routes against `data/OUTSTANDING.md` (tasks `#864`-`#937`), `data/empty_regions.jsonl`
(107 stamps), `data/catalogue.yaml` (392 rows; `our_status` is set on 24 rows, all
`known-class-member`; no `candidate-novel` row exists today), `src/cyclerfinder/search/*` and the
`#897` technique synthesis (`docs/notes/2026-10-04-897-technique-synthesis-papers-to-problems.md`,
proposals P1-P9), and kept the ranking.

**Corpus state.** At 17:20 AEDT the private corpus clone on this Mac was 87 commits behind its
origin and 71 files listed in `CORPUS_INDEX.md` were absent locally (290 files on disk), so the
first-draft claims about Henon, Bolotin, Gomez-Olle, Perko, Bradley-Russell, Landau, Ellison, McAdams,
Strange, Kumar Titan-Rhea, Oshima, Brown, the Jorba school, Leiva-Briozzo and Bruno 1981 were `[D]`.
**[F2]** The clone was updated before 19:00 AEDT (431 files on disk; every key the readers searched
for is present). Section 6 lists which `[D]` claims were then checked at the PDF; the tags in
sections 2-4 are left as written, with corrections marked "[F2]" where the source said something
else.

## 2. What the corpus says is NOT done (the raw material)

Compressed from the six extracts. Each line is a stated gap in a held paper, with where it is said.

**Heliocentric.**
- Russell & Strange 2009 (`russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-doi-10.2514-1.36610.pdf`):
  "Currently, we restrict the ideal model problem to be free-return cyclers only" (one body massless),
  and they cite their ref. [18] as "a recent cycler effort" that allows "flybys at both bodies in the
  ideal model cycler problem" `[P p.2]`. Ref. [18] is Pisarevsky, Kogan & Guelman 2008, "Interplanetary
  Periodic Trajectories in Two-Planet Systems", JGCD 31(3):729-739, doi 10.2514/1.30046. **It is not in
  the corpus** (cited in `docs/notes/2026-06-05-jones-aas17-577-vem-mining.md` as Jones's ref. 4 only).
- Russell & Ocampo 2006 (`russell-ocampo-2006-optimization-broad-class-ephemeris-model-earth-mars-cyclers-JGCD-29.pdf`):
  "the basic problem could be reformulated to define new classes of cyclers, such as those that
  incorporate Venus and or Mars flybys in the original idealized model" `[P p.13]`. Russell 2004
  dissertation: "The problem is immediately simplified by excluding flybys at other planets such as
  Venus" `[P p.185]`; Table 6.1 `[P p.188]` lists every conditional in the "global" claim (one time of
  flight per n-pi return, no back-to-back half-year returns, Earth v-inf below 10 km/s, at most three
  synodic periods).
- McConaghy, Longuski & Byrnes 2002 (`mcconaghy-longuski-byrnes-2002-analysis-broad-class-earth-mars-cycler-trajectories-AIAA-2002-4420.pdf`)
  p.8 "Other Possible Extensions": resonant (integer-year) legs, multiple intermediate gravity assists,
  Venus, Mars "or other planets" as intermediate gravity assists, repeat times other than two synodic
  periods, one-year and half-year "backflips" `[P p.8]`. McConaghy 2004 PhD p.165: "allow for Mars
  gravity assists"; p.167 Venus assists with the 6.4-year composite synodic period `[P]`. McConaghy et
  al. 2004 (JSR) p.5: "gravity-assist maneuvers at Mars, Venus, or the moon" `[P]`.
- Hollister & Rall 1970 (`hollister-rall-1970-periodic-orbits-NASA-CR.pdf`): Earth-Venus-Mars periodic
  orbits with direct returns at Venus and Earth "has not been carried out, but it is the next logical
  step" `[P thesis p.101]`; Mars-Venus periodic orbits were attempted and "this investigation simply
  failed to find any" `[P p.86]`, blamed on missing knowledge of returns "which traverse the Sun a
  different number of times than the planet" `[P p.85-86, rec. 2 p.126-128]` (the Russell n-pi return
  catalogue now supplies exactly that).
- Jones, Hernandez & Jesick 2017 (`jones-hernandez-jesick-2017-low-excess-speed-vem-triple-cyclers-AAS-17-577.pdf`):
  scope "1 or 2 synodic periods in a cycle, and a maximum of six flybys per cycle" `[P p.3]`; the
  one-synodic class had only subsurface near-misses, "It is possible, albeit unverified, that
  reasonably small deep-space maneuvers could enable one synodic period triple cyclers" `[P p.8]`.
  Already the target of `#866`/`#867`.
- No held paper names Mercury, Jupiter or asteroids as a body set for a repeating ballistic cycler.
  Mariner 10 (two 2:1 Mercury returns held by five trajectory corrections) is the only repeating
  encounter precedent for a non-Earth/Venus/Mars body `[P Dunne & Burgess 1978 ch.2; Giberson &
  Cunningham 1975 p.5]`.

**Planet-moon systems.**
- Pairs and triples with a computed REPEATING ballistic trajectory in any held paper: Jupiter
  Ganymede->Io, Ganymede->Europa, Ganymede->Callisto, Europa->Ganymede (Russell & Strange 2009
  Table 1) and the triples Io-Europa-Ganymede (Lynam & Longuski 2011, Hernandez et al. 2017) and
  Callisto-Ganymede-Europa (Liang et al. 2024); Saturn Titan->Enceladus only. Uranus, Neptune and
  Pluto: none `[P/D, reader table]`. Liang et al. 2024 exclude Io because "Io has the nearest orbit
  around Jupiter, which leads to a severe radiation environment" `[P p.2]`, so no held paper treats a
  Callisto-Ganymede-Io or Callisto-Europa-Io triple.
- Hernandez et al. 2017 (`hernandez-jones-jesick-2017-one-class-io-europa-ganymede-triple-cyclers-AAS-17-608.pdf`):
  "cyclers may exist which do not conform to the assumption that the orbit must intersect all three
  flyby bodies ... will be investigated in future work" `[P p.4]`; constant-energy families only
  `[P p.11]`.
- Anderson & Kumar 2024 (`anderson-kumar-2024-oberon-mmr-unstable-orbit-survey-aas-24-288.pdf`): "we
  do not yet look at the heteroclinics between mean motion resonances in the CCR4BP" `[P p.12]`;
  "Titania, Umbriel, and Ariel ... for which this study should be repeated" `[P p.19]`.
- Kumar, Anderson & de la Llave 2025 SIADS (`kumar-anderson-delallave-2025-gpu-connections-tori-perturbed-crtbp-siam-ads-24-219-arxiv-2109.14814.pdf`):
  "Another further step in the program would be to develop algorithms chaining the computed
  heteroclinic connections together" `[P p.36]`. No held paper closes a connection chain into a
  cycle in any four-body model; the only rigorous chain-closure is Wilczak & Zgliczynski Part II
  (planar CR3BP, Sun-Jupiter) `[P Thm 5.10 p.15]`.
- Blazevski & Ocampo 2012 (`blazevski-ocampo-2012-periodic-orbits-ccr4bp-invariant-manifolds-physica-d-241-1158-doi-10.1016-j.physd.2012.03.008.pdf`):
  "It would be interesting to see if there are periodic orbits in the five-body system under the
  assumption that omega_I = 2 omega_E" `[P journal p.1166-67]`. Baresi, Owen & Scheeres 2023
  (`baresi-owen-scheeres-2023-exploiting-laplace-resonance-tri-circular-problem-AAS-23-201.pdf`) build
  that Laplace-locked, time-periodic 5-body model but compute only libration-point substitutes and
  their tori; no resonant orbit, no fixed point of the stroboscopic map that meets two moons `[P p.13]`.
- Liang et al. 2026 ice-giant review lists "cycler orbits between satellites" as an open exploration
  `[P p.32, sec. 6; checked F2]`; Simon et al. 2026 (UOP) report that the baseline tour achieves
  "a minimum of 30 days between encounters" for downlink volume `[P sec. 3.3; F2: an achieved
  property of the baseline, not a stated requirement]`; Strange, Landau & Longuski 2013: "a v_inf of
  2 km/s is the lowest v_inf at either body to allow periapses near the rings" `[P p.13; F2: a
  v-infinity floor at Titania or Oberon, not a delta-v floor; ring limits are periapsis altitude
  >= 5,000 km and ring-plane crossing >= 51,140 km]`.

**Earth-Moon and periodic-orbit theory.**
- Ross & Roberts-Tsoukkas 2026 (`ross-roberts-tsoukkas-2026-stable-ballistic-prograde-cyclers-rtbp-arxiv-2606.29189.pdf`)
  Outlook: "the construction extends naturally to asymmetric, spatial, and exterior cyclers"; the
  persistence conjecture under eccentricity and the Sun `[P Outlook]`. Casoliva et al. 2010: "Further
  research is needed to organize and assess cyclers in the three-dimensional problem and to account
  for solar perturbations" `[P p.1640]`.
- Zhou, Anoe, Armellin, Qiao & Li 2025 (`2025-fixed-points-three-body-high-order-transfer-map-arxiv-2509.12671.pdf`)
  Table 6 lists near-periodic "Cycler2-like, P2g'-like, Cycler3-like, Cycler4-like-I and -II,
  Cycler5-like" initial conditions that "shuttle between both sides of the Earth", resemble the Ross &
  Roberts-Tsoukkas cyclers (their ref. [65]), and are "not strictly periodic, as they cannot be
  refined through differential correction" `[P Sec. V.B]`. The project's June mining note
  (`docs/notes/2026-06-13-high-order-transfer-map-2509.12671-mining.md`) does not mention Table 6.
- Hitzl & Henon 1977b announced and never delivered (in the corpus) the Earth-Moon stability study of
  second-species orbits: "the case mu* = 1/82.30 = 0.01215067 ... will be examined in detail. Of
  particular interest are the distances of closest approach to both the Moon and the Earth"
  `[P PDF p.21 = journal p.1039; checked F2]`.
  Bruno 1981 Table III ("Places of the families S in Figure 2 acceptable for flight round the Moon")
  lists 11 orbits in 10 families with e = 1 and a < 0.725 `[P PDF p.12 = journal p.266; checked
  F2]`; two of them (C12 0.55756, C25 0.52411) violate Bruno's other condition a > 0.57045 (eq. 16)
  and are listed with that caveat [F2];
  `#899` step 1 reproduces them (`second_species_arcs.e1_arcs`) and `#899` step 2 records that
  continuing the 1-2 family through Earth collision "needs Earth regularisation (`#928` regularises the
  Moon only)". Genova & Aldrin 2015 state their 3-petal cycler "was not shown to exist in the
  restricted three-body problem" `[P p.3]`; the row `genova-aldrin-2015-em-3petal-cycler` has no
  validation level.
- Font, Nunes & Simo 2002 leave untreated the case where quasi-collisions occur "near different
  locations along the orbit of the smallest primary" and say "The same ideas can be used" `[D p.141-142]`;
  Barrabes & Gomez 2003 compute no spatial exact orbit; Henon 1968 case (2) (rotation of the arc about
  the chord PQ) is never tabulated `[D]`. Bolotin 2006 Remark 1 allows several moving singularities
  with "proofs ... essentially the same", unproven `[D p.238]`.
- Franz & Russell 2022 and Restrepo & Russell 2018 both exclude cycler geometry by construction (near-
  Moon only; planar symmetric only), so absence from those databases is not evidence of novelty `[P]`.

## 3. Ranked routes

Aim stated by the lead: 8 to 15. I give 10 and park the rest in section 4. Rank is by expected
novelty-bearing yield per unit cost, after the policy and do-not-do checks.

### R1. Two-working-body heliocentric generator: Earth-Venus with Venus returns, Earth-Mars with Mars turns, Venus-Mars

- **Sources.** Russell & Strange 2009 `[P p.2]` (ideal model restricted to free returns; ref. [18]
  allows flybys at both bodies); Russell & Ocampo 2006 `[P p.13]` and Russell 2004 `[P p.185, 188]`
  (Venus and Mars flybys excluded from the idealised model); McConaghy et al. 2002 `[P p.8]`, 2004
  `[P p.5]`, McConaghy 2004 PhD `[P p.165, 167]`; Hollister & Rall 1970 `[P p.84-86, 98-101, 126-128]`
  (Earth-Venus-Mars and Mars-Venus attempts; recommendation 2); Hollister & Menning 1970
  (`hollister-menning-1970-periodic-swingby-earth-venus-JSR-7-10.pdf`) `[P p.1193-1198]` (15 Earth-Venus
  orbits with direct returns at BOTH planets).
- **Mechanism.** Russell's generating model (circular-coplanar, Earth-to-Earth n-pi and generic
  returns, target massless) is in `search/generic_return.py::RussellModel` with
  `periods_yr = {"E": 1.0, "M": 1.875}` and `search/cycler_search.py`. Three cells are unbuilt:
  (a) returns hosted at BOTH bodies (Venus bends 84 degrees at 5 km/s and 300 km, Earth 90, Mars 37,
  Mercury 29 `[C]`), so Earth-Mars cyclers in which Mars takes part of the turn, the exact
  heliocentric analogue of the two-massive-moon chain the project built for Titania-Oberon (`#890`);
  (b) Earth-Venus with the Russell n-pi return catalogue at Venus (Rall's missing ingredient);
  (c) Venus-Mars with Venus returns and Mars massless (Rall's 1969 failure, never revisited).
- **Why unpublished.** Every held source states the exclusion. Hollister used only full- and
  symmetric-return types at Venus and a 3.2-year basic period `[P p.1193]`.
- **Done near it.** `#388` lane (Earth-only returns, Mars massless, 4 rows, decisive negative on
  ephemeris closure); `#867` (Jones-style ephemeris grid, VEM, not dispatched); `#913`/`#915` (Appendix
  C reproduction and the Russell-Ocampo optimiser as published). The 15 Hollister-Menning rows are
  catalogued at validation null. No registry stamp covers Earth-Venus or Mars-as-working-body.
- **Hard gate first.** Acquire and read Pisarevsky, Kogan & Guelman 2008. It is the only published
  two-planet periodic-trajectory method with flybys at both bodies that the corpus knows of (by
  citation). If it already enumerates Earth-Mars with Mars turns or Earth-Venus with Venus returns,
  cells (a) and (b) become reproductions; cell (c) and the n-pi-at-Venus extension probably survive.
  **[F2] Second gate:** Russell & Strange 2009 p.5 states that the ideal-model enumeration was also
  run for Earth-Mars, Earth-Venus, Venus-Mars and Venus-Mercury "heliocentric calibration" sets, with
  the results in their ref. [26] (Russell & Strange, "Planetary moon cycler trajectories", AAS 07-118,
  2007; not held; free at the JPL open repository). That paper may already hold the ONE-working-body
  Earth-Venus and Venus-Mars catalogues, which would make cell (b)'s one-body part and R5's ideal
  model reproductions; the two-working-body cells (a), (c) are not covered by the R-S architecture.
- **First experiment.** Add Venus to `RussellModel` (period 0.61520 yr) and allow a return leg to be
  hosted at either body; enumerate Earth-Venus cyclers with up to three Earth-Venus synodic periods
  (584 d each). Positive control: recover Hollister's three circular-coplanar 3.2-year orbits (one Earth
  return, transfer, two Venus returns, transfer back `[P p.1193]`) and their 15 real-ephemeris
  descendants' turn angles and passing radii (Table 3). Then Earth-Mars with Mars turns on the Russell
  parents of `#913`. Then Venus-Mars. Turn gate `#888`/`#937` at every encounter.
- **Cost.** 3 to 5 agent-days build (GUESS), hours of compute. Ceiling V3 (REBOUND + DE440, `#866` lane).
- **Policy.** (a) and (b): (i) author-excluded classes, attribute Russell & Ocampo and Rall; (c): a
  published negative attempt, candidate-novel if it closes. Literal collision: Hollister rows, Jones
  rows, Pisarevsky 2008 (unread).
- **Likelihood.** 45 percent of at least one ideal-model member not in the held record; 25 percent at
  V3. Risk: the two-body turn budget at Mars is small, so (a) may only recolour Russell parents; (c)
  may repeat Rall's negative (which would then be a registry stamp with a published antecedent).

### R2. Io-containing Jovian triples with Liang's alternating-double-cycler construction

- **Sources.** Liang et al. 2024 (`liang-2024-callisto-ganymede-europa-triple-cyclers-JGCD.pdf`) `[P p.2-3, 20]`;
  Lynam & Longuski 2011 and Hernandez et al. 2017 (Io-Europa-Ganymede only).
- **Mechanism.** Liang's method needs two synodic periods sharing a hub moon to be near a small-integer
  ratio and absorbs the leftover by alternating two double cyclers. Io is excluded by the authors for
  radiation, a practicality ground. The hub ratios `[C]` from the registry periods: hub Callisto with
  Ganymede and Io 6.328 against 19/3 = 6.333 (0.08 percent); hub Ganymede with Callisto and Io 5.328
  against 16/3 (0.09 percent); hub Callisto with Europa and Io 2.280 against 16/7 (0.27 percent). All
  three are closer than Liang's 7:4 (1.5 percent), because Io-Europa-Ganymede are Laplace-locked.
- **Done near it.** `#526` built `genome/alternating_double_cycler.py` (`analyze_near_resonance`,
  `build_alternating_double_cycler_seed`) with Liang's own 7/4 as positive control and a Saturn demo;
  no Jovian triple other than Callisto-Ganymede-Europa was run. `#887` did the triple arithmetic for
  Uranus only. `#791` (blind Galilean sequence enumeration) is shelved; this is a specific construction,
  not a blind sweep, and the `#864` sec. 8 line against `#563`-class symmetric sweeps does not cover it.
  `#577`'s 0/36 ruling is about Russell-Strange two-moon closures.
- **First experiment.** Run `analyze_near_resonance` on the three Io-containing triples (minutes);
  build seeds with the `#526` operator for the best hub; close with the moon-cycler genome and the
  `#888` turn gate; positive control = `reproduce_before_search_gate()` on Liang's members (exists).
  Then JUP365-class ephemeris continuation (`nbody/jovian.py`) for ten cycles, as Liang did.
- **Cost.** 1 to 3 agent-days (GUESS); minutes to hours of compute. Ceiling V2-V3 at Jupiter.
- **Policy.** (i) author-excluded class (radiation), attribute Liang; literal collision against the
  Lynam/Hernandez Io-Europa-Ganymede rows (different body set, but check itineraries).
- **Likelihood.** 40 percent of an ideal-model closure; novelty high if it closes. Risk: Io's
  radiation makes the row mission-irrelevant (the catalogue records it anyway); W12 (Galilean hub
  energy) means several itineraries may close and must be deduplicated by signature.

### R3. Refine Zhou et al. 2025 Table 6 "cycler-like" fixed points with a multiple-shooting corrector

- **Source.** Zhou et al. 2025 `[P Sec. V.B, Table 6]`: five near-periodic Earth-Moon planar orbits with
  printed initial conditions (for example n = 2: x0 = 0.693257903603195, xdot0 = -0.0211935974664688,
  ydot0 = 0.585989070384689) that the authors could not refine with their single-section differential
  corrector and labelled invalid, while noting they resemble the Ross & Roberts-Tsoukkas cyclers.
- **Mechanism.** The project's asymmetric multiple shooting (`search/second_species_continuation.py`,
  `search/cr3bp_multiple_shooting.py`, the `#905` driver) can take non-perpendicular, multi-revolution
  seeds that a y = 0 single-crossing corrector cannot. If a seed closes it is either a Ross &
  Roberts-Tsoukkas or Casoliva or Leiva-Briozzo member (known-reproduction, and a correction to the
  paper's "invalid" label), or an orbit in no published family.
- **Done near it.** `#231`/`#450` mined the paper's method and adopted its Png' golden; the Table 6
  rows were never examined (grep of the mining note). `#905` found 70 Sun-forced counterparts of
  Earth-Moon families and is the nearest live lane.
- **First experiment.** Refine the five rows at the paper's C_J and mu with multiple shooting over the
  stated revolution count; on closure, match against `ross-rt-em-cycler-*`, `casoliva-*` and the
  Leiva-Briozzo atlas by signature. Positive control: the same table's valid Png' row and one Ross &
  Roberts-Tsoukkas row re-closed by the same code.
- **Cost.** 0.5 agent-day; minutes of compute. Ceiling V1-V2 (CR3BP).
- **Policy.** If not in any published family: (iii)-like, "first computed", attribute Zhou et al. for
  the seed. Earth-Moon: `#864` sec. 8 says do not re-sweep for novelty; this is not a sweep but a
  five-seed test of a printed claim, at near-zero cost.
- **Likelihood.** 60 percent that at least one row closes; 30 percent that one is outside the known
  families.

### R4. Earth-grazing second-species cyclers: Bruno's e = 1 arcs and Gomez-Olle double-collision orbits continued to the Earth-Moon mass

- **Sources.** Bruno 1981 Table III (11 Earth-Moon arcs with e = 1, a < 0.725, passing near both
  primaries at mu = 0) `[D]`; Hitzl & Henon 1977b's announced Earth-Moon stability study `[D p.1039]`;
  Gomez & Olle 1986 double-collision orbits (collide with both primaries, elliptic problem, mu = 0)
  `[D]`; Genova & Aldrin 2015 `[P p.3]` (3-petal 3:1 cycler "was not shown to exist in the restricted
  three-body problem"); Casoliva et al. 2010 reject Table 3 members that fly "through Earth" `[D]`.
- **Mechanism.** `#899` already reproduces the arcs (`second_species_arcs.e1_arcs`, Bruno Table III to
  5e-6) and continues single-arc and two-arc chains in mu, but its Levi-Civita propagator regularises
  the Moon only, and the 1-2 family "heads into Earth collision ... passing needs Earth regularisation"
  (`#899` step 2). The missing piece is a both-primary regularisation (Waldvogel 1967 B3, or Birkhoff,
  both digested) so that the e = 1 arcs and the near-Earth branches can be followed to mu = 0.01215.
- **Why unpublished.** Hitzl & Henon promised it and no later paper in the corpus delivers it; Bruno
  1981 stops at the mu = 0 formula; Casoliva discarded the Earth-grazing end; Genova & Aldrin make a
  non-existence statement without a search.
- **First experiment.** Build the Earth-side regularisation as an extension of `#928`/`#899`'s
  propagator (controls: radial fall, Kepler ellipse; **[F2]** Broucke's 7P/7A are families of the
  ELLIPTIC problem continued in e from the circular family C at Earth-Moon mu, TR 32-1360 report
  p.39 `[P]`, so they are an ER3BP control, not a circular one; the circular controls are `#899`'s
  reproduced Casoliva rows and Gomez & Olle 1991 II at mu = 1e-6); continue the 11 Bruno arcs and
  the 1-2c/d/e reverse paths; report
  perigee and perilune of every closed orbit; test for a 3:1 member with perigee and perilune near
  3000 km (Genova & Aldrin's geometry). Positive control: `#899`'s reproduced 2-1a, 3-2c, 7-3b/c.
- **Cost.** 3 to 5 agent-days (GUESS); CPU-hours. Ceiling V1-V2.
- **Policy.** (iii) theorem-generic (Perko existence) for members not in Casoliva's table; a confirmed
  CR3BP 3:1 Earth-grazing member would also be a sourced correction to a published claim. Earth-Moon:
  `#864` sec. 8 tension, but this is the stated continuation of `#899`, an owner-approved lane.
- **Likelihood.** 35 percent of a cycler-class orbit outside Casoliva's table; 50 percent of settling
  the Genova & Aldrin statement either way.

### R5. Venus-Mercury (and Earth-Venus-Mercury) ballistic cyclers on the real ephemeris

- **Sources.** No held paper names the body set. Precedent: Mariner 10's two 2:1 Mercury returns held
  by five trajectory corrections `[P Dunne & Burgess 1978 ch.2; Giberson & Cunningham 1975 p.5]`;
  Hughes et al. 2014 call the Hollister-Menning orbits "Earth-Venus cycler trajectories" `[P p.2]`;
  catalogue rows `mariner-10-venus-mercury` and `bepicolombo-earth-venus-mercury` are `mga_tour`.
- **Mechanism.** Venus-Mercury synodic period 144.6 d; Hohmann Venus to Mercury gives v-inf 5.8 km/s
  at Venus and 6.8 km/s at Mercury; Venus bends 84 degrees at 5 km/s, Mercury 29 degrees at 5 km/s and
  13 at 8 `[C]`. Mercury's eccentricity 0.206 makes a circular-coplanar generator poor; the Jones-style
  real-ephemeris Lambert grid of `#867` (multi-revolution legs, v-inf continuity, altitude window,
  integer synodic closure) is the right tool, with Venus as the working body and Mercury returns of
  the Mariner type as the resonant legs.
- **Done near it.** Nothing in the project. The registry has no Mercury stamp; `literature_check.py`
  has no Mercury cycler anchor (the gate would return "not found", necessary not sufficient).
  **[F2] Hard gate:** Russell & Strange 2009 p.5 says a Venus-Mercury ideal-model free-return search
  was run (results in the unheld AAS 07-118). Read it before dispatch; the real-ephemeris part of R5
  (Mercury e = 0.206, i = 7 deg) is outside that paper's circular-coplanar model either way.
- **First experiment.** Web literature pass first (MESSENGER and BepiColombo resonant-return design
  papers; Yen 1985 "reverse V-EGA" named by Barrabes & Gomez 2003 `[P p.145]`): is any repeating
  Venus-Mercury or Mercury-only ballistic sequence in print. Then the `#867` enumerator with body set
  {V, Me} and {E, V, Me}, cycle = k times 144.6 d, over three Venus-Mercury opportunities. Positive
  control: reproduce Mariner 10's V-Me-Me-Me geometry (dates 1974-02-05, 1974-03-29, 1974-09-21,
  1975-03-16) as a near-ballistic chain with its known small corrections.
- **Cost.** 1 to 2 agent-days on top of `#867` (GUESS); CPU-hours. Ceiling V3.
- **Policy.** No source to attribute; candidate-novel only after the literature gate clears with a
  Mercury-aware anchor set. Heliocentric multi-planet is the roadmap's own lane (`#864` item 3), not
  the Sun-planet CR3BP regime of sec. 8.
- **Likelihood.** 25 percent of a ballistic multi-cycle member (Mercury's eccentricity and small turn
  argue against; Venus's large turn argues for).

### R6. Earth-NEA one-working-node cyclers: a resonant near-Earth asteroid as the massless passive target

- **Sources.** Russell & Strange 2009 architecture (working body with free returns, target
  "considered massless") `[P p.2-4]`; Ozaki et al. 2022 (`ozaki-2022-neural-network-surrogate-global-cycler-search-arxiv-2111-11858.pdf`)
  build Earth-asteroid-Earth blocks on Russell free returns but target a different asteroid every leg
  so that no periodic closure is needed `[P p.1-4]`; Adamo 2025 (`adamo-2025-spanning-earth-mars-chasm-synodic-resonant-waypoints-AIAA-houston-LnL.pdf`)
  lists 493 catalogued NEAs between Earth and Mars and uses one as a powered loiter waypoint `[P p.9, 11]`;
  de la Fuente Marcos 2018 (Arjuna 1:1 co-orbitals) `[D]`.
- **Mechanism.** The heliocentric analogue of Titan-Enceladus: Earth is the only body that bends, the
  asteroid is met on a transit leg, and the cycle closes after an integer number of Earth-asteroid
  synodic periods. The bend gate is irrelevant at the target. Candidates are NEAs whose period is
  within a fraction of a percent of p/q Earth years (q at most 3) and whose orbit crosses the cycler's.
- **Done near it.** `#308` (`search/asteroid_leveraging.py`) treated NEAs as FLYBY nodes and found,
  correctly, that they cannot bend; the passive-target role was not posed. No registry stamp.
- **First experiment.** JPL SBDB query for NEAs with |P/P_E - p/q| < 0.5 percent, q <= 3, perihelion
  below 1 au and aphelion above 1 au; for each, enumerate Russell free returns at Earth
  (`search/free_return.py`, `bodies` parameter) with the asteroid as the target; phase-match on the
  real asteroid ephemeris over 2030-2060. Positive control: reproduce one Ozaki Earth-asteroid-Earth
  block (DESTINY+ Phaethon) as a single leg.
- **Cost.** 2 to 3 agent-days (GUESS); CPU-hours. Ceiling V3 (DE440 plus SBDB elements).
- **Policy.** (ii) known architecture at a never-treated body set, attribute Russell & Strange. OWNER
  QUESTION: is an asteroid an admissible cycler endpoint under the spec (the catalogue has no such
  row; the `#864` mission-scope decision did not cover it).
- **Likelihood.** 50 percent that a geometric closure exists for some NEA; 30 percent that it survives
  the real ephemeris for three cycles. Novelty is policy-dependent.

### R7. Exact periodic triple cycler in the Laplace-locked tri-circular (5-body) model

- **Sources.** Blazevski & Ocampo 2012 `[P journal p.1166-67]`; Baresi, Owen & Scheeres 2023
  `[P p.12-13]` (the Jupiter-Io-Europa-Ganymede geometry "repeats after one full revolution of
  Ganymede"); Hernandez et al. 2017 Tables (EGGIE 4-synodic, 0.70 m/s; EGIEIE) `[P p.7-11]`; Bradley &
  Russell 2014 continuation `[D]`.
- **Mechanism.** In the 4:2:1 locked model the vector field is periodic with Ganymede's period, so a
  triple cycler is a fixed point of the stroboscopic map and can be corrected as a true periodic orbit
  with all three inner moons massive. Seed: Hernandez's ballistic ideal-model member; continuation in
  the moons' mass (kappa) with the `two_moon_periodic_890.py` machinery extended to three moons, or the
  `core/crnbp.py` Laplace-locked field.
- **Done near it.** `#714`-`#736` N = 5 CRNBP torus at Europa 3:4 (a torus, not an encounter-bearing
  periodic orbit); `#890` (two-moon periodic orbit at Uranus by mass continuation); `#480` Laplace
  triple-cycler construction bug and fix. No exact 5-body periodic triple cycler exists in the corpus.
- **First experiment.** Rebuild Hernandez's EGGIE at conic level (gate exists in
  `search/moon_cycler_genome.py`), then continue in kappa with the stroboscopic closure condition.
  Positive control: Blazevski & Ocampo's m2-Lyapunov periodic orbit in the CCR4BP `[P p.1158]`
  (periodic in both rotating frames).
- **Cost.** 3 to 5 agent-days (GUESS); CPU-hours.
- **Policy.** The itinerary is published, so the result is a V-tier lift of
  `hernandez-2017-jovian-ieg-triple-family` (validation null today) and a model-native object, not a
  new row; it is also the validation lane R2's Io-containing closures would need. Ranked here for
  that reason.
- **Likelihood.** 50 percent convergence; 0 novelty on its own.

### R8. Lunar gravity assists inside heliocentric cyclers

- **Sources.** McConaghy et al. 2004 `[P p.5]` ("gravity-assist maneuvers at Mars, Venus, or the moon");
  McConaghy et al. 2002 `[P p.8]` (backflips); McConaghy 2004 PhD `[P p.162]` (tool blind spot at
  180 + n 360 degree Earth-to-Earth transfers). Russell & Ocampo 2006: 9, 39 and 74 of 203 parents
  under 1, 10 and 300 m/s per seven cycles `[P results]`.
- **Mechanism.** A lunar swingby at an Earth passage changes the effective Earth turn by up to a few
  degrees and the v-inf by up to about 1 km/s at no propellant cost, so near-ballistic Russell-Ocampo
  descendants (the 10-300 m/s tier) may become ballistic. No held paper computes a cycler with a lunar
  assist; the project has no lunar-GA node in its heliocentric legs (grep: `backflip` absent from
  `src/`).
- **Done near it.** `#913`/`#915` (Appendix C and the optimiser as published) will produce the parent
  set and costs; `#388` lane negatives are all Earth-only.
- **First experiment.** After `#913`: take the parents with 1-100 m/s printed cost, add a Moon node at
  each Earth encounter (patched conic, DE440 Moon), and re-close ballistically. Positive control: the
  Galileo VEEGA Earth-1 2:1 return `[D damario-1992]` reproduced with and without the Moon to show the
  node's effect is modelled.
- **Cost.** 3 to 5 agent-days (GUESS). Ceiling V3.
- **Policy.** (i), attribute McConaghy; many rows would be V-tier lifts of catalogued parents rather
  than new rows. Earth-Mars: `#864` sec. 8 says do not re-sweep for novelty; the value here is mostly
  validation.
- **Likelihood.** 30 percent of a ballistic descendant not in print; 60 percent of V-tier lifts.

### R9. Exterior-realm cyclers (through L2), the named future work of Ross & Roberts-Tsoukkas

- **Sources.** Ross & Roberts-Tsoukkas 2025 `[P pp.14-15]` and 2026 `[P Outlook]`; Rawat et al. 2025
  exterior MMR `[P Discussion]`; the mu = 0.1 (1,3) exterior cycler is Fig. 3 of the 7-page VSGC
  student paper (`roberts-tsoukkas-2026-vsgc-multiorbiter-cyclers-student-summary.pdf`), figure
  only, no numbers `[P PDF p.5; checked F2]`. **[F2]** The file named `...-journal.pdf` is
  byte-identical to that student paper (md5 f5851cdf...), so no journal version is held; the arXiv
  2606.29189 version prints Table I with numeric members (mu = 0.1 (3,2): x0 = -0.694376003123377,
  C = 3.573367616904619, T = 12.295263874014290) and shows NO exterior cycler. The positive control
  for the first experiment is therefore the (3,2) Table I row, and the (1,3) figure is a shape target.
- **Mechanism.** The project's (k1,k2) construction (`#315`, `#504`, `#549`, `#656`, `#657`) follows
  the L1-tube interior recipe. Exterior members (L2 tubes, apogee beyond the Moon) are a distinct
  class; at Pluto-Charon (mu = 0.109) the mu = 0.1 figure makes them known-class-member; at Earth-Moon
  they are stated-not-done.
- **First experiment.** Rediscover the mu = 0.1 (1,3) figure orbit with the existing corrector at
  mu = 0.1 (shape match only); then sweep (k1,k2) up to 3 at Earth-Moon mu through L2.
- **Cost.** 2 to 3 agent-days (GUESS). Ceiling V2.
- **Policy.** Not-enumerated class; Earth-Moon sec. 8 tension; likely `known-class-member` once the
  journal figure is counted as a published member of the class.
- **Likelihood.** 70 percent existence; 20 percent novelty.

### R10. Spatial second-species seeds corrected at a physical mass ratio

- **Sources.** Barrabes & Gomez 2002, 2003 (`barrabes-gomez-2003-three-dimensional-pq-resonant-orbits-second-species-solutions-cmda-85-2-doi-10.1023-A1022098510161.pdf`):
  spatial family only for phi_0 = +-pi/2 (vertical flyby), C_J = 2 + (p/q)^(2/3), p < q, matching to
  O(mu^(1-alpha)), "no numerical correction of a seed" `[P/D]`; Henon 1968 case (2) `[D]`; Casoliva 2010
  future work `[P p.1640]`.
- **Mechanism.** Add the out-of-plane seed to `#899`'s generator and correct at Earth-Moon or
  Jupiter-Ganymede mass with the spatial transition matrix (`second_species_lc.py` has planar and
  vertical STMs).
- **Done near it.** `#438`/`#444`: 3D lifts of KNOWN planar cyclers ruled `known-class-member` via the
  vertical-bifurcation mechanism (spec 16.4). A vertical-flyby second-species orbit may or may not
  connect to a planar family through a vertical-critical orbit; if it does, the same ruling applies.
- **First experiment.** Correct the 1-2 and 3-5 spatial seeds of Barrabes & Gomez 2003 Fig. 4 at
  mu = 1e-6, then continue in mu; at the target mu compute the family's connection to the planar
  p-q family (vertical stability index crossing).
- **Cost.** 3 to 4 agent-days on top of `#899` (GUESS).
- **Policy.** Likely `known-class-member` under `#444` unless disconnected from the planar family; the
  value is the first 3D second-species cycler at a physical mass, a census item with a sourced seed.
- **Likelihood.** 50 percent of a closed spatial orbit; 15 percent novelty.

## 3a. Cross-paper transfer analysis: technique x case matrix [F2]

**Method.** Four Sonnet readers rebuilt, from the six cluster extracts, the digests and (where the
extract lacked them) the PDFs, an inventory of techniques (what each computes, needs and ASSUMES) and
of orbital cases (system, model, mass ratio, eccentricity, energy, resonance). I merged their 95
technique blocks into 32 technique rows and their 110 case lines into 16 case columns, then marked
every cell: `P` published in that case (cite), `D #NNN` done by this project, `E` empty with the
technique's assumptions holding, `e` empty but low value or assumptions only partly met, `X`
assumptions fail, `-` not meaningful. The full matrix is
`docs/notes/2026-10-05-938-technique-case-matrix.csv` (one row per technique; the `note` column
carries the per-row judgment). Assumption checks are the readers' and mine from the stated hypotheses;
mass ratios and eccentricities quoted in the CSV are standard values, not read from the papers unless
the inventory says so.

**Case columns.** C01 Earth-Moon CR3BP (mu 0.01215); C02 Earth-Moon with the Sun or lunar e
(BCR4BP/QBCP/HR4BP/ER3BP); C03 Sun-Earth-Mars circular-coplanar one-working-body patched conic;
C04 the same on the real ephemeris; C05 Sun-Venus-Earth(-Mars) with Venus as a working body; C06
Sun-Venus-Mercury (e 0.206, i 7 deg); C07 Sun-Earth with a resonant NEA as massless target; C08
Sun-Jupiter CR3BP/ER3BP (mu 9.5e-4); C09 Jupiter-Europa or Jupiter-Ganymede single-moon CR3BP/ER3BP;
C10 Jupiter two-moon CCR4BP and the Laplace-locked N = 5 field; C11 Jovian multi-moon patched-conic
cyclers; C12 Saturn-Titan (mu 2.4e-4, e 0.029) with Rhea or Enceladus; C13 Uranus-Titania/Oberon
(mu 3.5-4e-5) and Umbriel; C14 Neptune-Triton (mu 2.1e-4, e about 0); C15 Pluto-Charon (mu 0.109);
C16 binaries and star-planet systems (mu 0.1-0.5; 1e-3).

**Non-empty rows and columns.** Columns C01, C02 and C10/C11 are the densest in `P` (Earth-Moon
theory and forced models; Jovian tours and tori); C12, C13, C14 carry `P` only in the resonant-family,
tour-design and torus rows (T09, T10, T17, T18, T21, T25, T28); C05, C06, C07 have almost no `P` at
all beyond Hollister-Menning, Jones and the flown Mariner 10. Rows with the most `E` cells whose
assumptions hold are T13 (Font-Nunes-Simo complete enumeration), T12 (Bolotin shadowing), T02
(two-working-body corrector), T22 (DA transfer-map enumeration), T15 (homoclinic-shadowing cyclers),
T31 (spatial second-species seeds) and T23 (Keplerian map). Rows that are design aids or produce
non-cycler objects (T09, T24, T26, T27, T28, T32) are kept in the CSV for completeness and ranked
nowhere.

**Identifications found while building the matrix** (each removes a candidate or changes a rank):
- Jones 2017's two-synodic VEM classes (EMEVE etc.) ARE Liang's alternating-double-cycler
  construction with hub Earth: 4 S(E,V) = 2336 d, 3 S(E,M) = 2340 d, 7 S(V,M) = 2337 d, a 0.2 percent
  mismatch against Liang's 1.5 percent `[C, reader-helio-moons]`. So "Liang at the heliocentric
  case" is published, not a transfer; R2's Io-containing triples remain the open Liang cell.
- Russell's n-pi and generic Earth returns are Henon's mu = 0 consecutive-collision arcs (T10 at
  C03 is `P`-equivalent). The second-species literature and the heliocentric cycler literature are
  the same object at two mass ratios, which is why T12 and T13 transfer to the heliocentric case.
- Casoliva's Class-2 cyclers (symmetric periodic orbits shadowing the L1 Lyapunov homoclinic) are
  exactly what `#868`/`#782` closed at Neptune-Triton. The published technique and the project's
  result were never linked; Titan, Ganymede and Titania are the direct repeats (within `#872`'s
  pilot scope).
- Casoliva's Class-1 generator at a planet-moon pair produces the symmetric p:q resonant families
  that Restrepo & Russell 2018 tabulate for 24 systems and Anderson & Lo, Vaquero, Anderson & Kumar
  and Spear study as resonant flyby orbits (`P/e` cells in T10). The unpublished residue at moons is
  the asymmetric members, the tight members below the grids' cut-offs and the spatial vertical-flyby
  branch of Barrabes & Gomez 2003 (R10), not the symmetric families.
- Russell & Strange 2009 p.5 states that Earth-Venus, Venus-Mars and Venus-Mercury ideal-model
  searches were run, results in the unheld AAS 07-118 (T01 at C05, C06 marked `P?`). This gates R1(b)
  and R5 (section 3, "[F2]" lines).

**Ranked transfer cells** (a technique proven in one case applied to a quite different case where
nobody has used it and its assumptions hold; `X-` prefix distinguishes them from the R routes).

- **X1. Hollister's two-working-body date-residual corrector at Jupiter (T02 -> C11).** Hollister &
  Menning found 15 Earth-Venus periodic orbits with direct returns at BOTH planets by iterating N
  encounter dates on a v-infinity-mismatch residual `[P pp.1193-1198]`; Rall's own limit is that the
  returning bodies must be massive enough to bend `[P p.124]`. Russell & Strange restrict the moon
  problem to a massless target and only assert that a massive-target cycler "is possible as long as
  target flyby altitudes are sufficiently high" `[P p.4]`. Ganymede-Callisto (mu 7.8e-5 and 5.7e-5)
  and Ganymede-Europa satisfy Rall's condition; both moons bend tens of degrees at the R-S excess
  speeds `[C]`. Project state: `#577` ruled 0/36 on SYMMETRIC two-leg two-moon closures at Jupiter
  (sec. 8 line); `#890`/`#895` did a two-moon chain at Uranus by mass continuation, not this
  corrector. Positive controls: Hollister-Menning's 15 rows (same corrector, heliocentric) and the
  R-S Ganymede-Callisto rows in the one-body limit (the `#888` gate already reproduces their
  altitudes). Policy: (i) author-excluded class (R-S's stated model restriction), attribute Russell &
  Strange and Hollister. Cost 3 to 5 agent-days (GUESS), shared generator with R1. Likelihood 40
  percent of an ideal-model closure not in the held record; radiation is not an objection at
  Ganymede-Callisto (Europa Clipper flies them). Ceiling V2-V3 (JUP365 continuation, T07).
- **X2. Font, Nunes & Simo's complete second-species enumeration at Europa or Ganymede (T13 ->
  C09).** Their 2009 method finds, at one Jacobi constant, THE FULL SET of second-species periodic
  orbits up to a maximal time as intersections of the stable and unstable manifolds of the collision
  singularity on the pericentre section, with an S/T symbolic label on each `[P fns2009 abstract;
  X]`; validated for mu from 1e-6 to 1e-3 (full set at 1e-4) and never at Earth-Moon mu, which the
  authors say is too large `[R 2002 pp.139-141]`. Europa (2.5e-5), Ganymede (7.8e-5), the Uranian
  moons (4e-5) and Titan (2.4e-4) are INSIDE the validated range. Nobody has run it at a physical
  moon. Output: every flyby-chain periodic orbit at that moon and energy, including the asymmetric
  S-arc chains that no symmetric grid database (Restrepo & Russell, Franz & Russell) can contain.
  Positive control: their 41-orbit table at C = 2.8, mu = 1e-4 (`#896` already reproduces 8 printed
  orbits), then the symmetric subset must match Anderson & Lo's Europa 3:4 and 5:6 families.
  Policy: members outside the published symmetric families are (iii) theorem-generic objects not in
  a published family, candidate-novel with attribution to Font, Nunes & Simo; the `#855` periselene
  rule decides cycler class per member (the orbits pass within mu^0.4 of the moon by construction).
  Cost 5 to 8 agent-days (GUESS; the collision-manifold machinery is the new part; `second_species_lc.py`
  supplies the regularised flow). Likelihood 85 percent of a complete table; 40 percent that it
  contains cycler-class members outside the published families. Ceiling V1-V2 (CR3BP; JUP365 for V3).
- **X3. Bolotin's elliptic shadowing theory plus the ER3BP corrector at Saturn-Titan (T12 + T21 ->
  C12).** Bolotin 2005 proves, for the planar elliptic restricted problem with FIXED eccentricity and
  small mass ratio, that chains of Kepler collision arcs with a nondegenerate summed action are
  shadowed by periodic orbits `[P p.4 of Bolotin 2006 for the regime statement; D bolotin-2005]`; it
  is stated for Sun-Jupiter-type bodies and computed nowhere. Titan (e 0.029, mu 2.4e-4, e/mu about
  120) is the cleanest physical case in the solar system for the "eps fixed, mu small" regime; the
  Earth-Moon case (e 0.055) is excluded from novelty sweeps by `#864` sec. 8 and Titania (e 0.0011,
  about 0.03 mu) is in the wrong regime. The corrector side is published too: Broucke 1969 (ER3BP
  symmetric families at Earth-Moon mu, 7P/7A `[P report p.39]`), Gomez & Olle 1991 II (elliptic
  second-species families at mu = 1e-6), Peng & Xu and Martinez-Cacho (Titan ER3BP QSOs). The project
  holds `core/er3bp.py` with a published control, `#917` (Casoliva rows in e) is registered, and
  `#899`'s arc chains give the seeds. Positive controls: Gomez & Olle II's printed elliptic families
  at 1e-6; Martinez-Cacho's Titan triplets. Pass: an ER3BP periodic orbit of period 2 k pi at Titan's
  e whose arcs form a prescribed S/T word with at least two distinct Titan passages. Policy: (iii)
  theorem-generic; a two-arc member is a Titan-only resonant-hopping cycler as an exact ER3BP object
  ("first second-species periodic orbits of the elliptic problem at a physical moon"). Cost 4 to 6
  agent-days on top of `#899` and `#917` (GUESS). Likelihood 55 percent of convergence; novelty
  moderate (single-moon, Titan-only objects are mission-relevant for Cassini-class tours).
- **X4. The DA transfer-map enumerator at Pluto-Charon and Titan (T22 -> C15, C12).** Zhou et al.'s
  method is complete within a domain and admits asymmetric fixed points `[P Sec. I, VI]`; the
  project's `#450` enumerator is already stamped at Earth-Moon and Sun-Jupiter (Hilda 3:2). Every
  Pluto-Charon negative (`#504`, `#549`, `#656`: "exactly one cycler family") came from the SYMMETRIC
  (k1,k2) construction, so "empty" there is conditional on symmetry. Running the enumerator over the
  cycler domain at mu 0.109 is hours of compute. Positive control: Ross-RT's (3,2) Table I row at
  mu = 0.1. Policy: the owner closed the Pluto-Charon lane (`#864` sec. 8: no resonant-manifold lane;
  parked here as the cheapest reopening). Titan is the un-gated alternative (`#633`'s 0/16,375 was
  also symmetric). Likelihood 30 percent of an asymmetric cycler-class fixed point at either.
- **X5. Casoliva Class-2 homoclinic-shadowing cyclers at Titan and Ganymede (T15 -> C12, C09).**
  The `#868` recipe (seed near a base orbit's perpendicular point, first y = 0 crossing with small
  xdot, one member per seed) is this technique; it has run at Neptune-Triton only. Controls: the two
  closed Neptune-Triton orbits and Casoliva's He1. Already inside `#872`'s "G1 corrector + symmetric
  pilot" scope by my reading; listed so the link to the published class is recorded. Theorem-generic,
  (iii). Cost 1 to 2 agent-days per moon once `#872` exists.
- **X6. Bolotin-MacKay nondegeneracy as a persistence predictor for the 203 Russell parents (T12 ->
  C03; a higher-order link).** `#388`'s wall is family selection: ideal parents do not continue to
  finite mass. Bolotin & MacKay's condition for a chain to be shadowed is a nonsingular tridiagonal
  Hessian of the summed fixed-energy Lambert action `[R 2000 pp.62, 66-67]`, with Remark 1 of Bolotin
  2006 extending the field to "several planets" `[P p.4]`. Computing that Hessian for each parent is
  cheap (hours) and would sort the parents into shadowable and degenerate BEFORE any ephemeris
  attempt; the `#897` finding that whole-revolution returns are degenerate at fixed time and regular
  at fixed energy is the same mathematics. Hypothesis-level; a clean negative (no correlation with
  Russell & Ocampo's published cost tiers 9/39/74) would itself be informative. No novelty on its own.
- **X7. Wilczak-Zgliczynski covering relations to certify the Neptune-Triton homoclinic-accumulating
  orbits (T20 -> C14).** The `#636` machinery exists (planar CR3BP, validated integration); it has
  only been pointed at Sun-Jupiter controls and the `#646` SE-EM negative. Applying it to the two
  `#868` orbits would turn "first computed" into "proved to exist", which is the wording a preprint
  needs. Publication value, no new row. Cost 3 to 5 agent-days (GUESS).
- **X8. Kumar's stroboscopic invariant-circle solver at Sun-Earth forced by Venus (T17 -> C05).**
  Reader-helio-moons' hint 4: the CCR4BP needs no commensurability (Kumar 2021 treats Europa-Ganymede
  as an approximate 2:1 with irrational ratio `[P Kumar 2021 text]`), so Sun-Earth with Venus as the
  concentric circular perturber (e 0.007) is a valid heliocentric CCR4BP and the quasi-periodic
  counterparts of Earth's resonant free-return orbits (VISIT 4:5, 2:3) could be computed with
  `search/pertbp_strob_889.py`. A dynamical-object bet (`#864` item 13), not a cycler; ranked last.

**Rejected transfer cells (assumptions fail or cell already published).** Lynam's phase-angle
method anywhere but Jupiter (needs an exact Laplace lock); Kumar's tori at Neptune (no circular
coplanar second moon; `#864` sec. 8); RRT (k1,k2) at Titan-class mu (`#627`/`#633` negatives; the
reader's prediction is the same for Ganymede and Triton); the QBCP Fourier model at any moon system
(its tables exist only for Sun-Earth-Moon); linking numbers in any forced model (authors say it does
not apply); Tisserand graphs at Venus-Mercury (e, i); Casoliva seeds at Pluto-Charon (mu too large,
`#504`); Liang's construction at Venus-Earth-Mars (it is Jones 2017).

**Links that remove a named obstacle.** (1) `#388` family-selection wall: T06 "ask for a cost, not a
closure" (`#897` P5) plus X6's nondegeneracy sort plus T08's grid enumeration (`#897` P9). (2) S1L1
off-basin (`[[project_s1l1_realeph_closure_blocker]]`): S1L1 is a two-arc chain; Bradley & Russell's
kappa continuation with the pseudo-arclength step its authors named but did not use (T07) and the
fixed-energy arc solver (T12 remark) are the two published pieces the June attempts lacked. (3)
`#577`/`#563` symmetric two-moon walls: X1's asymmetric N-date corrector is the published alternative
to symmetric closure. (4) The Pluto-Charon "exactly one family" result: X4's asymmetric enumeration is
the method under which that negative was never tested.

**Re-ranked list (R routes and X cells together).** Rank by expected novelty-bearing yield per unit
cost after the gates named in each entry; the first two share a generator.

| Rank | Item | Gate before dispatch | Cost (GUESS) | P(novel row) |
|---|---|---|---|---|
| 1 | R1 two-working-body heliocentric generator (Venus returns; Mars turns; Venus-Mars) | Pisarevsky 2008 and AAS 07-118 read | 3-5 d | 0.25-0.45 |
| 2 | X1 the same corrector at Jupiter Ganymede-Callisto, Ganymede-Europa | none beyond X1's controls | +2 d on R1 | 0.4 |
| 3 | X2 Font-Nunes-Simo complete enumeration at Europa or Ganymede | `#896` controls pass | 5-8 d | 0.4 |
| 4 | R2 Io-containing Jovian triples (Liang construction) | none | 1-3 d | 0.4 |
| 5 | X3 elliptic second-species cyclers at Titan (Bolotin 2005 + ER3BP corrector) | `#899` step 3, `#917` | 4-6 d | 0.3 |
| 6 | R3 Zhou Table 6 refinement by multiple shooting | none | 0.5 d | 0.3 |
| 7 | R4 Earth-grazing second species (both-primary regularisation) | owner: Earth-Moon tension | 3-5 d | 0.35 |
| 8 | X4 DA enumerator at Titan (Pluto-Charon if the owner reopens it) | owner for Pluto | 1-2 d | 0.3 |
| 9 | R6 Earth-NEA passive-target cyclers | owner: admissible endpoint | 2-3 d | 0.3 (policy) |
| 10 | R5 Venus-Mercury on the real ephemeris | AAS 07-118 read; `#867` built | 1-2 d on `#867` | 0.25 |
| 11 | X5 homoclinic-shadowing cyclers at Titan and Ganymede | `#872` | 1-2 d per moon | 0.2 (thin) |
| 12 | R7 exact periodic triple cycler in the Laplace-locked N = 5 field | none | 3-5 d | 0 (V-tier lift) |
| 13 | X6 Bolotin-MacKay nondegeneracy sort of the 203 Russell parents | none | 0.5-1 d | 0 (diagnostic) |
| 14 | R8 lunar gravity assists inside heliocentric cyclers | `#913` | 3-5 d | 0.3 of a row; mostly validation |
| 15 | R9 exterior-realm (k1,k2) cyclers at Earth-Moon | owner: Earth-Moon tension | 2-3 d | 0.2 |
| 16 | R10 spatial second-species seeds at a physical mass | `#899` | 3-4 d | 0.15 |
| 17 | X7 covering-relation proof of the Neptune-Triton orbits | `#868` adjudicated | 3-5 d | 0 (publication) |
| 18 | X8 Venus-forced Earth resonant tori | none | 2-3 d | 0 (dynamical object) |

## 4. Ideas considered and rejected or parked, with the reason

- **Variable-energy and not-all-three-bodies Galilean triple cyclers** (Hernandez 2017 future work
  `[P p.4, 11]`). This is the shelved `#791` blind enumeration under another name; W12 says energy
  cannot prune the Galilean alphabet. Rejected as a duplicate; R2 is the specific construction instead.
- **Two-moon chains with same-moon resonant returns at Uranus** (`#897` P8). Already registered as the
  two-moon part of `#899` (ledger line "P8 is `#899`"). Not re-proposed.
- **Four-body connection chains closed into a cycle** (Kumar SIADS `[P p.36]`, Olikara thesis
  `[P p.97]`). Covered by `#872` (G1 corrector) and `#886`. Not re-proposed.
- **Four-body tori and connections at untreated Uranian or Saturnian base moons** (Anderson & Kumar
  `[P p.19]`, Kumar 2025 Titan-Rhea `[D p.14]`). Roadmap item 12 and `#886` scope (3). Not re-proposed.
- **Secondary-resonance stroboscopic periodic orbits of the Oberon interior 4:3 family that cross
  Titania's orbit** (Kumar 2023 AAS `[P p.18-19]`, Anderson & Kumar `[P p.17]`). Kumar's own body set
  and "work in progress" by the authors; `#864` sec. 8 excludes Titania-Oberon CCR4BP.
- **Pluto-Charon one-working-node with Styx, Nix, Kerberos or Hydra as passive targets** (period
  ratios to Charon 3.16, 3.89, 5.04, 5.98 `[C]`; Howett et al. 2021 show figure-only periodic orbits
  with "close encounters with Pluto and all its satellites" `[D]`). Class-level prior art with no
  numbers; `#320` already returned "Pluto Hydra-Nix V0-known"; `#864` sec. 8 says no Pluto-Charon lane.
  Parked: if the owner reopens Pluto, this is the cheapest cell (the Russell-Strange genome plus the
  Pluto-Charon CR3BP already in `verify/pluto_charon_realeph.py`).
- **Jones one-synodic VEM class rescued by small deep-space manoeuvres** (`[P p.8]`). Powered; a
  Merrill-class forced trajectory under the spec, and `#864` sec. 8 excludes powered "novel cycler"
  sweeps. Record it as a `#867` sub-cell only if the owner's powered-cycler admissibility rule allows.
- **Earth-Moon exterior Franz-Russell / Restrepo-Russell complement grid search** (asymmetric,
  Earth-circulating). A re-sweep of Earth-Moon (`#864` sec. 8, W2); the databases' exclusions are
  recorded so that absence from them is not misread as novelty.
- **Stability of Casoliva's Class-2 L1-homoclinic-shadowing orbits; Broucke E1/F and Kumar-Moreno
  Fig. 12/13 families as cyclers; Leiva-Briozzo atlas label mapping.** Published families; census only;
  useful inside `#905`'s deduplication, not a discovery route.
- **Wittal 2022 Earth-to-NRHO 5-petal pair.** Figure-only in an unstated model; reproduction with no
  sourced initial condition; Earth-Moon.
- **Belbruno Theorem 3.60 periodic symbol sequences (repeating pseudoballistic captures).** `#378`,
  `#681`, `#908` cover the capture sweeps; `#891` stamped `#378` method-invalid.
- **Mars-Venus periodic orbits as a standalone route.** Folded into R1 cell (c) because the generator
  is the same.
- **Simon et al. 2026 "30 days between encounters" and Strange 2013's 2 km/s ring floor.** Not routes;
  record them as mission-relevance filters on any Uranian candidate. **[F2]** Both are weaker than
  the first draft said: the 30 d is an achieved property of the UOP baseline tour, not a requirement,
  and the 2 km/s is a v-infinity floor for periapses near the rings, not a delta-v. The six catalogued
  Uranian rows (7.75-15 d legs, 0.9-2.2 km/s) sit at or below both numbers and would need a
  ring-clearance check (periapsis >= 5,000 km above the rings, crossing >= 51,140 km) before any
  mission-relevance claim.

## 5. Corpus gaps (papers to acquire), in priority order [F2]

Checked against `CORPUS_INDEX.md`, the digest bodies and the 431 filenames on disk (keys grepped
2026-10-05 19:35 AEDT); none of the items below is held. DOIs marked CONFIRMED were resolved through
the Crossref API to the stated title and authors; UNCONFIRMED means the record was found only by web
search or has no DOI. Items the first draft listed that ARE held were removed: Bolotin 2005 (CMDA
93:343, held and digested), Lantoine & Russell 2011, Guillaume 1973/1975, Rhouma & Chicone 2000,
Jorba & Villanueva 1997, Sanaga & Howell 2025, Peng & Xu 2015/2017, Miceli et al. AIAA 2024-1280,
Campagnola et al. 2014 Acta Astronautica.

1. **[Received 2026-10-05, filed under `#960`.]** **Russell, R. P. & Strange, N. J. (2007)**, "Planetary moon cycler trajectories", AAS 07-118,
   AAS/AIAA Space Flight Mechanics Meeting, Sedona; Advances in the Astronautical Sciences 127. No
   DOI (conference); free copy at the JPL open repository, handle 2014/40318 (UNCONFIRMED by DOI;
   located by web search). Unlocks: it holds the Earth-Venus, Venus-Mars and Venus-Mercury
   ideal-model free-return searches that Russell & Strange 2009 p.5 cite as ref. [26]; gates R1(b),
   R5 and the literal-collision check for X1. [`#960` correction: it holds NO Earth-Venus set (Table 1);
   see `2026-10-05-digest-russell-strange-2007-aas-07-118-planetary-moon-cyclers.md`.]
2. **[Received 2026-10-05, filed under `#960`.]** **Pisarevsky, D. M., Kogan, A. & Guelman, M. (2008)**, "Interplanetary Periodic Trajectories in
   Two-Planet Systems", J. Guidance, Control, and Dynamics 31(3):729-739, doi 10.2514/1.30046
   (CONFIRMED). Unlocks: the only published two-planet periodic method with flybys at BOTH bodies;
   gates R1 cells (a), (b) and X1.
3. **[Received 2026-10-05, filed under `#960`.]** **Perko, L. M. (1974)**, "Periodic Orbits in the Restricted Three-Body Problem: Existence and
   Asymptotic Approximation", SIAM J. Applied Mathematics 27(1):200-237 [`#960` correction: issue 1, July 1974, per the PDF and Crossref; was printed here as 27(2)], doi 10.1137/0127016
   (CONFIRMED). Unlocks: the original second-species existence theorem that every held Perko paper
   cites; needed to state X2/X3/R4 novelty claims as "theorem-generic" with the theorem in hand.
4. **[Received 2026-10-05, filed under `#960`.]** **Campagnola, S., Buffington, B. B., Lam, T., Petropoulos, A. E. & Pellegrini, E. (2019)**,
   "Tour Design Techniques for the Europa Clipper Mission", JGCD 42(12):2615-2626, doi
   10.2514/1.G004309 (CONFIRMED). Unlocks: literal-collision check for X1, R2 and R7 (Ganymede-
   Callisto repeated legs, the "Callisto pi-transfer").
5. **[Received 2026-10-05, filed under `#960`.]** **Liang, Y., Xu, M., Peng, K. & Xu, S. (2020)**, "A cislunar in-orbit infrastructure based on p:q
   resonant cycler orbits", Acta Astronautica 170:539-551, doi 10.1016/j.actaastro.2020.02.029
   (CONFIRMED). Unlocks: Earth-Moon p:q cycler prior art for R3, R4 and R9 literal collision
   (quasi-periodic p:q cyclers refined in the bicircular model).
6. **Binder, D. & Arnas, D. (2024)**, "Reliable and Repeatable Transit Through Cislunar Space Using
   2:1 Resonant Spatial Orbits", JGCD 47(9):1973-1979, doi 10.2514/1.G007800 (CONFIRMED); arXiv
   2304.13584. Unlocks: spatial 2:1 Earth-Moon cycler prior art (R10, R3 collision check).
7. **Kevorkian, J. & Lancaster, J. E. (1968)**, "An Asymptotic Solution for a Class of Periodic
   Orbits of the Restricted Three-Body Problem", Astronomical Journal 73:791-806, doi 10.1086/110701
   (CONFIRMED). Unlocks: the earliest matched-asymptotic Earth-Moon periodic flyby orbits; R4 and X2
   attribution.
8. **[Received 2026-10-05, filed under `#960`.]** **Hitzl, D. L. (1977)**, "Generating Orbits for Stable Close Encounter Periodic Solutions of the
   Restricted Problem", AIAA Journal 15(10):1410-1418, doi 10.2514/3.60808 (CONFIRMED). Unlocks: the
   stability side of the Hitzl-Henon programme that 1977b announced; R4's "announced, never delivered"
   claim must be re-checked against it.
9. **[Received 2026-10-05, filed under `#960`.]** **Arenstorf, R. F. (1963)**, "Existence of Periodic Solutions Passing Near Both Masses of the
   Restricted Three-Body Problem", AIAA Journal 1(1):238-240, doi 10.2514/3.1516 (CONFIRMED).
   Unlocks: the existence theorem behind Genova & Aldrin's non-existence remark (R4).
10. **Oshima, K. (2022)**, "Continuation and stationkeeping analyses on planar retrograde periodic
    orbits around the Earth", Advances in Space Research 69(5):2210-2222, doi
    10.1016/j.asr.2021.12.020 (CONFIRMED). Unlocks: the 12:11 initial conditions the held Oshima
    papers lack (`#905`/P3 controls).
11. **Bruno, A. D. & Varin, V. P. (2006)**, "On families of periodic solutions of the restricted
    three-body problem", Celest. Mech. Dyn. Astron. 95:27-54, doi 10.1007/s10569-006-9021-1
    (CONFIRMED). Unlocks: symmetric families for all mu in [0, 1/2], the cross-check Leiva & Briozzo
    name; X4's Pluto-Charon and R9 controls.
12. **Henon, M. & Guyot, M. (1970)**, "Stability of Periodic Orbits in the Restricted Problem", in
    Giacaglia (ed.), Periodic Orbits, Stability and Resonances, Reidel, pp. 349-374, doi
    10.1007/978-94-010-3323-7_33 (CONFIRMED). Unlocks: critical (stability-boundary) orbits of
    families f, g, h, i, l, m for all mu, the limits mu -> 0 and mu -> 1; controls for X4 and R9.
13. **Huang, S.-S. (1962)**, "Preliminary study of orbits of interest for moon probes", Astronomical
    Journal 67:304-310, doi 10.1086/108730 (CONFIRMED); **Huang, S.-S. & Wade, C. Jr. (1963)**,
    "Preliminary study of periodic orbits of interest for moon probes. II", AJ 68:388-391, doi
    10.1086/108988 (CONFIRMED). Unlocks: the 1960s Earth-Moon periodic flyby lineage (R4 collision).
14. **Schwaniger, A. J. (1963)**, "Trajectories in the Earth-Moon Space with Symmetrical Free-Return
    Properties", NASA TN D-1833; **Hoelker, R. F. & Winston, B. P. (1968)**, "A Comparison of a Class
    of Earth-Moon Orbits with a Class of Rotating Kepler Orbits", NASA TN D-4903 (both NTRS, no DOI;
    UNCONFIRMED). **Newton, R. R. (1959)**, "Periodic orbits of a planetoid passing close to two
    gravitating masses", Smithsonian Contributions to Astrophysics 3:69-78 (no DOI found;
    UNCONFIRMED). Unlocks: the remaining primer-cited Earth-Moon cycler ancestors (R4).
15. **[Held as the 2015 journal version: IEEE AES Magazine 30(7), doi 10.1109/MAES.2015.140119,
    digested 2026-10-03; the 2026-10-05 upload was an md5 duplicate and was not filed.]**
    **Campagnola, S., Boutonnet, A., Martens, W. & Masters, A. (2014)**, "Mission Design for the
    Exploration of Neptune and Triton", 24th ISSFD, paper S6-1; free PDF at issfd.org
    (ISSFD_2014/ISSFD24_Paper_S6-1_campagnola.pdf; no DOI; UNCONFIRMED by DOI). Unlocks: the
    Neptune-Triton tour prior art the `#864` review named as unread before any `#868` writeback.
16. **Yen, C.-W. L. (1989)**, "Ballistic Mercury orbiter mission via Venus and Mercury gravity
    assists", J. Astronautical Sciences 37(4):417-432 (ADS 1989JAnSc..37..417Y; no DOI;
    UNCONFIRMED). Unlocks: the reverse Delta-V-EGA and near-resonant Mercury returns, the Mercury
    anchor the literature gate lacks for R5.
17. **Finley, T., Barth, E., Howett, C., Zangari, A., Tapley, M., Scherrer, J. & Stern, A.**, "An
    Orbital Tour of Pluto and Its Moons", J. Spacecraft and Rockets, cited by Stern et al. 2020 as
    "to be published"; no Crossref record found on 2026-10-05 (possibly never published;
    UNCONFIRMED). Unlocks: the only numerical Charon-assist tour; needed only if Pluto reopens (X4).
18. **[Held as journal versions: AAS 03-509 = McConaghy, Landau, Yam & Longuski 2006, JSR 43(2),
    doi 10.2514/1.15215; AAS 03-510 = Chen et al. 2005, JSR 42(5). The 2026-06-17 Hintz digest records
    that the catalogue covers both; the 2026-10-05 upload of the 2006 paper was an md5 duplicate and was
    not filed.]** **McConaghy, T. T., Yam, C. H., Landau, D. F. & Longuski, J. M. (2003)**, "Two-Synodic-Period
    Earth-Mars Cyclers with Intermediate Earth Encounter", AAS 03-509; **Chen, K. J. et al. (2003)**,
    AAS 03-510 (conference, no DOI; UNCONFIRMED). Unlocks: check whether they add members to the held
    AIAA 2002-4420/4422 content (R1(a), R8).
19. **Menning, M. D. (1968)**, MS thesis (MIT), Earth-Venus periodic orbits; **VanderVeen, A. A.
    (1969)**, E-V-M-V-E flyby trajectories; **Minovitch, M. A. (1965)**, "Utilizing Large Planetary
    Perturbations for the Design of Deep-Space, Solar-Probe, and Out-of-Ecliptic Trajectories", JPL
    TR 32-849 (NTRS 19660005935; the first draft's "1967" was wrong). No DOIs; UNCONFIRMED except
    Minovitch's NTRS id. Unlocks: the pre-Rall Earth-Venus-Mars record for R1.
20. **Olle, M. (1989)**, doctoral thesis (Univ. de Barcelona), evolution of all second-species
    families; **Llibre & Pinol**, collision orbits ("to appear" in Llibre 1982); the shortened CMDA
    version of **Leiva & Briozzo 2006** (only the arXiv preprint is held). No DOIs; UNCONFIRMED.
    Unlocks: R4 and X2 family bookkeeping.
21. Still wanted from `#909` and relevant here: Kumar, ASC 2026 Paper 1023 (Titan-Rhea untargeted
    moons); Komachi, ASC 2026 Paper 688 (EML2-SEL2 cycler in the BCR4BP); Kumar & Anderson 2026
    ISSFD (UNCONFIRMED, not searched).

## 6. What was verified at the source versus what was inferred [F2]

**Verified by a reader opening the PDF in the first pass** (tagged `[P]` in sections 2-4): Russell &
Strange 2009 p.2 restriction and ref. [18] identity (re-read by the first author in the text layer);
Russell & Ocampo 2006 p.13 and Russell 2004 pp.185-191; McConaghy 2002 p.8, 2004 p.5, PhD pp.162-167;
Hollister & Rall pp.71-101, 124-128; Hollister & Menning 1970; Jones 2017 pp.1-11; Liang 2024 pp.2-3,
20; Hernandez 2017 pp.2-11; Anderson & Kumar 2024 pp.12, 19; Kumar SIADS p.36; Blazevski & Ocampo
journal pp.1158, 1165-67; Baresi-Owen-Scheeres p.13; Zhou et al. 2025 Sec. V.B and VI; Ross &
Roberts-Tsoukkas 2025 conclusion and 2026 Outlook; Casoliva 2010 p.1640; Genova & Aldrin p.3; Barrabes
& Gomez 2002 p.406 and 2003 p.145; Ozaki 2022 pp.1-4, 20; Adamo 2025 slides; Wolf & Smith 1995 Table
2; Dunne & Burgess 1978 ch.2 and Giberson & Cunningham 1975.

**Source checks run in the second pass (reader-sourcecheck, 18 claims, PDFs now on disk):**
- CONFIRMED as written: Hitzl & Henon 1977b p.1039 quotation (C2); Gomez & Olle 1986 double-collision
  orbits in the planar ELLIPTIC problem at mu = 0 (C3; also: the circular limit has none for h < -1,
  two at h = -1, two or four for -1 < h < 0, and continuation to mu > 0 is stated as under study);
  Rosengren et al. 2026 primer pp.51-52 name every early-lineage reference the draft listed, plus
  Arenstorf 1963 and Hitzl 1977 (C4; no DOIs printed; references transcribed into section 5); Henon
  1968 case (2) rotation about PQ not tabulated (C6; the rotation is of the plane pi3 about the line
  Delta through M1, journal pp.378-379 and 390); Bolotin 2006 Remark 1 p.238 (C7; wording is "several
  planets", not "moving singularities"); Font, Nunes & Simo 2002 pp.141-142 (C8); Kumar 2025
  Titan-Rhea p.14 future work (C10: Dione, Tethys, Enceladus "remain to be carried out"); Liang 2026
  review sec. 6 (C11); Howett 2021 figure-only orbits (C14); D'Amario 1992 Galileo 2-year Earth-Earth
  leg (C15; the paper never writes "2:1"); Perko: the held papers are 1967, 1976a/b, 1977, 1981a/b and
  Breakwell-Perko 1974, all circular planar small mu; the 1974 SIAM existence theorem is NOT held (C17).
- CORRECTED (and propagated into sections 2-4 as "[F2]"): Bruno 1981 Table III count 11 holds but two
  orbits violate his eq. 16 bound (C1); the `...-journal.pdf` of Roberts-Tsoukkas & Ross is the
  student summary, the (1,3) exterior cycler is figure-only there and absent from the arXiv version,
  which prints Table I numbers (C5); Broucke's 7P/7A are elliptic-problem families in TR 32-1360
  (1969), not in TR 32-1168 (1968), and are continued in e from the circular family C (C9); Simon et
  al. 2026's 30 d is achieved, not required (C12); Strange 2013's 2 km/s is a v-infinity, not a
  delta-v (C13); Leiva & Briozzo 2005 is two QBCP orbits from one family (not an atlas), the 2006
  preprint is the atlas, and the 2008 paper prints Table 1 numeric states (x = -0.836915310 section)
  for its RTBP orbits plus QBCP results (C16); no held Oshima paper prints 12:11 initial conditions
  (C18: they are in the unheld ASR 69:2210).
- Still `[D]` (not re-checked at the PDF; not load-bearing for a top-ranked item): Henon 1997/2001
  books, Gomez & Olle 1991 I/II page cites, Bradley & Russell 2014 numerical results, Rosengren
  primer body text beyond pp.51-52 and the references, McAdams 2011, Landau 2023/2025, Ellison 2025,
  Brown 2024/2025, the Jorba school, Leiva & Briozzo 2006 preprint body, Kumar Titan-Rhea body beyond
  p.14, Miceli 2024.

**Technique-inventory evidence.** The four inventory readers tagged each block `[P]` (PDF opened this
run), `[X]` (taken from the first-pass extract, itself page-cited) or `[D]` (digest). Roughly a third
of the 95 blocks are `[P]`; the rest rest on the extracts. Load-bearing `assumes` lines for X1-X8
were checked by me against the inventories and, for Kumar 2021 (single-frequency forcing without
commensurability), Bolotin 2006 p.4 (regime statement), Font-Nunes-Simo 2009 (fixed C, alpha = 0.4,
mu range) and Blazevski & Ocampo pp.1-3, against the PDF text by a reader this run.

**Calculations (`[C]`)**: turn capacities from GM, radius and a 200-300 km altitude floor; synodic
periods from sidereal periods; Jovian hub ratios; VEM composite period (2336/2340/2337 d); Hohmann
v-infinity Venus-Mercury; e/mu ratios for the Bolotin regime. None was checked by a second agent.

**Checked in the repository**: `RussellModel.periods_yr`; `free_return.py` `bodies`;
`alternating_double_cycler.py` and the `#526` note; `asteroid_leveraging.py` (`#308`); `#899` step-1
and step-2 ledger text; `#905` result text; `#887` stage-one note; `#627`/`#629`/`#633` (RRT at Titan
negatives, symmetric); `#504`/`#549`/`#656` (Pluto-Charon, symmetric); `#450`/`#523`/`#527`/`#532`
(DA-HOTM enumerator and stamps); `#636`/`#646` (covering-relation machinery); `#868`/`#782` recipe;
`#890`/`#895` results (no row); `#889` positive control; `#917`/`#928`/`#930`-`#937` registrations;
registry stamps grepped for kk/titan/ganymede/triton/charon/venus/mercury/second/casoliva/hotm/strob
(no Earth-Venus, Mercury, Earth-NEA, Europa second-species or Titan ER3BP stamp exists).

**Not done**: no code run beyond `grep`, `ls`, `wc`, `pdftotext` (readers) and one 30-line arithmetic
script (first pass); no literature gate run; no probability calibrated against a measured unit cost;
the matrix cells marked `E` were not checked against `literature_check.py` anchors.

## 7. Papercuts

First pass (corpus-review-fable):
- Teammates cannot name sub-agents (`name` parameter refused: "team roster is flat"); six dispatches
  had to be re-sent without names.
- All six sub-agents ended "without delivering a report through SubagentHandback"; their output
  survived only because the brief told them to write the extract to a scratch file incrementally.
- The private corpus clone was 87 commits behind origin; 71 indexed PDFs were missing locally, so a
  third of the extraction was digest-only.
- `CORPUS_INDEX.md` abbreviates some filenames with "..." so a filename-to-disk check cannot be
  scripted exactly.
- Short search keys (strange-2013, landau-2023, zhang-2025, stern-2020, kumar-2025) do not match the
  actual filenames; one key (zhang-2025) matches an unrelated paper.
- Only about 30 of 290 PDFs had `.txt` sidecars; readers regenerated text with `pdftotext` in scratch.
- The `Read` tool caps a file at 25k tokens; each 60-80 kB extract needed two reads.
- zsh: unquoted `--include=*.py` fails ("no matches found"); macOS `sed -i` needs a suffix.

Second pass (corpus-review-fable-2):
- All five Sonnet readers again ended without a SubagentHandback report (one reported "agent that
  spawned you is no longer running"); every result was recovered from the readers' own scratch
  files. Incremental scratch writing is mandatory, not optional.
- `roberts-tsoukkas-ross-2026-stable-prograde-em-cyclers-journal.pdf` is byte-identical to the VSGC
  student summary; the filename says "journal" and misled the first pass.
- `hitzl-henon-1977b-...` is the Acta Astronautica stability paper, not the "Critical generating
  orbits" paper (CM 15:421, a separate file); the first draft cited the wrong one by title.
- Broucke's 7P/7A live in TR 32-1360 (1969, elliptic), not TR 32-1168 (1968); the digest's family
  labels were read as circular-problem controls.
- The Crossref API returned HTTP 429 after about 12 calls in one minute; DOI confirmation had to be
  spread over three rounds (Yen 1989 and Newton 1959 still have no DOI record).
- `ugrep` on this Mac rejects a moderately long alternation with "exceeds complexity limits"; the
  OUTSTANDING status scan had to be split.
- Two bullet formats coexist in `data/OUTSTANDING.md` (`- **\`#NNN\`**` and `- \`#NNN\` —`), so a
  status grep anchored on one format silently misses half the tasks.
- The epub-only papers (Simon 2026, Howett 2021) have no page numbers; citations are by section.
- Bruno 1981 is an image-only ADS scan; the reader read it page by page.
