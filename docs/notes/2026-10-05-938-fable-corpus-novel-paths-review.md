# #938 Corpus review: untried routes to novel cycler and quasi-cycler orbits

**Date:** 2026-10-05 (AEDT). **Task:** `#938` (analysis only; no code, no catalogue edit, no dispatch).
**Author:** corpus-review-fable (Fable), with six Sonnet readers over disjoint corpus themes.
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

**Corpus state (reported to the lead at 17:20 AEDT).** The private corpus clone on this Mac is 87
commits behind its origin (`git fetch` run, no pull). 71 files listed in `CORPUS_INDEX.md` are not in
the local `papers/` directory (290 files on disk), including almost the whole 2026-10-03/04/05 digest
wave (Henon, Bolotin, Gomez-Olle, Perko, Bradley-Russell, Landau, Ellison, McAdams, Strange, Kumar
Titan-Rhea, Oshima, Brown, the Jorba school, Leiva-Briozzo 2005 and 2008, Bruno 1981). Every claim
about those papers below is `[D]`. A fast-forward pull would restore them; I did not touch that
repository.

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
  `[D sec. 6, pp.30-32]`; Simon et al. 2026 (UOP) require "a minimum of 30 days between encounters"
  `[D sec. 3.3]`; Strange, Landau & Longuski 2013 set a ring-avoidance floor of about 2 km/s for
  same-moon resonant sequences `[D]`.

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
  particular interest are the distances of closest approach to both the Moon and the Earth" `[D p.1039]`.
  Bruno 1981 gives 11 Earth-Moon e = 1 arcs that pass near BOTH primaries at mu = 0 `[D Table III]`;
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
  propagator (controls: radial fall, Kepler ellipse, Broucke 1969 7P/7A rows at Earth-Moon mu which
  approach the larger primary `[D]`); continue the 11 Bruno arcs and the 1-2c/d/e reverse paths; report
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
- **Done near it.** Nothing. The registry has no Mercury stamp; `literature_check.py` has no Mercury
  cycler anchor (the gate would return "not found", necessary not sufficient).
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
  exterior MMR `[P Discussion]`; Roberts-Tsoukkas & Ross 2026 journal shows an mu = 0.1 (1,3) exterior
  cycler as a figure with no numbers `[D]`.
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
  record them as mission-relevance filters on any Uranian candidate (the six catalogued rows run
  7.75-15 d legs at 0.9-2.2 km/s and would fail both).

## 5. Corpus gaps (papers to acquire), in priority order

1. **Pisarevsky, Kogan & Guelman 2008**, "Interplanetary Periodic Trajectories in Two-Planet Systems",
   JGCD 31(3):729-739, doi 10.2514/1.30046. Russell & Strange's ref. [18] (flybys at both bodies in the
   ideal model), Jones's ref. [4]. Gates R1. Not held, not in `#909`.
2. **Finley et al.**, "An Orbital Tour of Pluto and Its Moons" (JSR), Stern et al. 2020 ref. [5]
   `[P p.4]`: the only numerical source for the Charon gravity-assist tour. Needed if Pluto reopens.
3. **Campagnola et al. 2019**, JGCD 42(12), Europa Clipper tour using Ganymede-Callisto cyclers and the
   "Callisto pi-transfer" (named in Yang et al. 2023 `[D p.5]`). Bears on R2 and R7 literal collision.
4. **Liang et al. 2020**, Acta Astronautica 170:539 (p:q resonant cycler orbits, Earth-Moon),
   **Binder & Arnas 2024**, **Kevorkian & Lancaster 1968**, **Schwaniger 1963**, **Hoelker & Winston
   1968**, **Huang 1962/63**, **Newton 1959**: the early and recent Earth-Moon cycler prior art named
   by the Rosengren et al. 2026 primer `[D pp.51-52]`. Bear on R3, R4, R9 literal collision.
5. **Bruno & Varin 2006** (symmetric periodic orbits for 0 <= mu <= 1/2), named by Leiva & Briozzo
   2006 as the unexploited cross-check `[D]`.
6. **Yen 1985** ("reverse V-EGA" Mercury approach), MESSENGER and BepiColombo resonant-return design
   papers (McAdams et al.; Jehn et al.): the literature gate for R5 has no Mercury anchor today.
7. **Menning 1968** MS thesis, **VanderVeen 1969** (E-V-M-V-E flyby trajectories), **Minovitch 1967**
   (interplanetary transportation network): the pre-Rall Earth-Venus-Mars record, named by Hollister &
   Menning `[P p.1193]` and Rall `[P p.71]`. Bear on R1 and `#867`.
8. **McConaghy, Yam, Landau & Longuski AAS 03-509** and **Chen et al. AAS 03-510**, cited by Hintz 2023
   `[P p.401]`; check whether they duplicate the held AIAA 2002-4420/4422 content.
9. **Olle 1989 thesis** (evolution of all second-species families, Gomez & Olle II p.153 `[D]`),
   **Llibre & Pinol** (collision orbits, "to appear" in Llibre 1982 `[D]`), **Oshima 2022a** (12:11
   initial conditions), **Henon & Guyot 1970**, any Hitzl-Henon follow-up on Earth-Moon stability
   windows (verify whether it was ever published).
10. Already on `#909`'s wanted list and still relevant here: Kumar ASC 2026 Paper 1023 (Titan-Rhea
    untargeted moons), Komachi ASC 2026 Paper 688 (EML2-SEL2 cycler in the BCR4BP), Kumar & Anderson
    2026 ISSFD.

Also: restore the 71 indexed-but-missing PDFs on this machine by pulling the private corpus clone
(section 1).

## 6. What I verified at the source versus what I inferred

- Verified by a reader opening the PDF (tagged `[P]` above): Russell & Strange 2009 p.2 restriction and
  ref. [18] identity (I re-ran `pdftotext` on the held PDF and read the reference entry myself);
  Russell & Ocampo 2006 p.13 and Russell 2004 pp.185-191; McConaghy 2002 p.8, 2004 p.5, PhD pp.162-167;
  Hollister & Rall pp.71-101, 124-128; Hollister & Menning 1970; Jones 2017 pp.1-11; Liang 2024 pp.2-3,
  20; Hernandez 2017 pp.2-11; Anderson & Kumar 2024 pp.12, 19; Kumar SIADS p.36; Blazevski & Ocampo
  journal pp.1158, 1165-67; Baresi-Owen-Scheeres p.13; Zhou et al. 2025 Sec. V.B and VI; Ross &
  Roberts-Tsoukkas 2025 conclusion and 2026 Outlook; Casoliva 2010 p.1640; Genova & Aldrin p.3;
  Barrabes & Gomez 2002 p.406 and 2003 p.145; Ozaki 2022 pp.1-4, 20; Adamo 2025 slides; Wolf & Smith
  1995 Table 2; Dunne & Burgess 1978 ch.2 and Giberson & Cunningham 1975 (Mariner 10).
- Digest-only (`[D]`): everything about Bruno 1981, Hitzl & Henon 1977b, Gomez & Olle 1986/1991, Henon
  1968/1997/2001, Font-Nunes-Simo, Bolotin, Perko, Bradley & Russell 2014, the UOP-era Uranian papers,
  Kumar Titan-Rhea, Howett 2021, Rosengren 2026 primer, Leiva & Briozzo 2005/2008, Oshima, Brown,
  Jorba school, Miceli 2024. These PDFs are indexed but absent locally (section 1).
- My own calculations (`[C]`): turn capacities from GM, radius and a 200-300 km altitude floor;
  synodic periods from sidereal periods; Jovian hub-synodic ratios; Pluto moon period ratios;
  Hohmann v-inf Venus-Mercury. None was checked by a second agent.
- Checked in the repository: `RussellModel.periods_yr`, `free_return.py` `bodies` parameter,
  `alternating_double_cycler.py` API and the `#526` note, `asteroid_leveraging.py` docstring (`#308`),
  the `#899` step-1/step-2 ledger text and scan-script arguments, the `#231` mining note (no Table 6),
  `#905` result text, `#887` stage-one note (Uranus only), registry stamps (no Venus, Mercury, Earth-Venus
  or Earth-NEA stamp), catalogue rows for Hollister, Jones, Hernandez, Genova, Casoliva, Ross-RT,
  Mariner 10 and BepiColombo.
- Not done: no code was run beyond `grep`, `ls`, `pdftotext` and a 30-line arithmetic script; no
  literature gate was run; no probability here is calibrated against a measured unit cost.

## 7. Papercuts

- Teammates cannot name sub-agents (`name` parameter refused: "team roster is flat"); six dispatches
  had to be re-sent without names.
- All six sub-agents ended "without delivering a report through SubagentHandback"; their output
  survived only because the brief told them to write the extract to a scratch file incrementally.
- The private corpus clone is 87 commits behind origin; 71 indexed PDFs are missing locally, so a
  third of the extraction was digest-only.
- `CORPUS_INDEX.md` abbreviates some filenames with "..." so a filename-to-disk check cannot be scripted
  exactly (71 is a lower bound on exact matches, an upper bound on true absences is 87 files).
- Short search keys (strange-2013, landau-2023, zhang-2025, stern-2020, kumar-2025) do not match the
  actual filenames; one key (zhang-2025) matches an unrelated paper.
- Only about 30 of 290 PDFs have `.txt` sidecars; readers regenerated text with `pdftotext` in scratch.
- The `Read` tool caps a file at 25k tokens; each 60-80 kB extract needed two reads.
- zsh: unquoted `--include=*.py` fails ("no matches found"); macOS `sed -i` needs a suffix.
