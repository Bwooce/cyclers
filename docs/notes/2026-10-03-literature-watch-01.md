# Literature watch 01 (2026-10-03)

First dated watch pass. Web research only; no code, catalogue or test was touched. Facts only: this note does not judge whether any project result is novel.

Method. Searches were run with WebSearch, the arXiv query API (export.arxiv.org, sorted by submission date, newest first; the newest entries returned are dated 2026-09-29), the Crossref API, publisher landing pages and the two AAS/AIAA Astrodynamics Specialist Conference (ASC) programmes. Held status was checked by grepping `docs/notes/CORPUS_INDEX.md` for the arXiv number or DOI, the first-author surname and a distinctive title word; the anchors in `src/cyclerfinder/search/literature_check.py` were also listed. Paywalled pages (Springer, AIAA ARC) redirected to a login or consent handshake and were not pursued.

Protected areas (numbering from the task brief): 1 = Uranian two-moon quasi-cyclers; 2 = Neptune-Triton 4:5 homoclinic and symmetric periodic orbits; 3 = four-body torus connections at moons (CCR4BP, linking number); 4 = VEM triple and Earth-Mars cyclers; 5 = Jovian and Saturnian moon cyclers; 6 = Earth-Moon cycler-type orbits and resonance networks.

Observation: query 1 (a plain web search) returned the project's own GitHub repository (Bwooce/cyclers) and cyclers.space as results; those are not external literature and are not counted below.

## 1. Queries run and what each returned

| # | Channel | Query (abridged) | Returned |
|---|---|---|---|
| 1 | WebSearch | Uranus moon cycler Ariel Umbriel Titania Oberon repeating | Uranian dynamical-history papers (2005.12887, 2403.17896, 2509.24631, 2511.18776), project GitHub; no external cycler paper |
| 2 | WebSearch | Uranus Orbiter and Probe satellite tour design 2026 arXiv | 2505.05514 (Simon et al. poll), PSJ ae680c (held), 2310.14514, 1907.02963 |
| 3 | WebSearch | Neptune-Triton 4:5 resonant homoclinic periodic orbits | Motion-primitive papers (held), CMDA retrograde 4/5 Neptune MMR (small bodies, 10.1007/s10569-022-10106-3), Barrabes-Mondelo-Olle heteroclinic classics; nothing at Triton 4:5 |
| 4 | WebSearch | Triton tour design resonant orbits Neptune 2026 Bosanac | Miceli-Bosanac JAS 2026 (held), AAS 24-161, SciTech 2024-1280, Campagnola ISSFD 2014 |
| 5 | WebSearch | Kumar Anderson de la Llave CCR4BP whiskered tori heteroclinic | Acta Astronautica 2023, 2109.14815, 2309.06073, 2105.11100, 2601.00149 (all held) |
| 6 | WebSearch | linking number knot theory heteroclinic quasi-periodic tori four-body moons | Owen-Baresi AAS 23-110 and ISSFD 2024 (both held), 2507.06123 (not held) |
| 7 | WebSearch | Venus-Earth-Mars triple cycler enumeration 2026 | Jones/Hernandez/Jesick 2017 (held), selenianboondocks EVMVE-2034 blog (April 2026) |
| 8 | WebSearch | Earth-Mars cycler new class 2026 Longuski | 2000s Longuski-group papers (held); nothing dated 2025-2026 |
| 9 | WebSearch | Jupiter Saturn moon cycler double triple Europa Ganymede Callisto 2026 | Liang-Yang-Li JGCD 2025 (held), Laplace-resonant triple cyclers 2012, arXiv 2601.00786 (disk formation, not trajectory) |
| 10 | WebSearch | Russell Strange cycler Titan Enceladus ballistic free-return | Russell-Strange 2009 (held), Pony Express JSR (held), 10.1007/s10569-024-10228-w |
| 11 | WebSearch | Earth-Moon cycler periodic orbits resonance network heteroclinic Ross 2026 | 2605.31543 (held), 2509.12675 (held), AAS 25-621 (held), 2606.26367 (not held) |
| 12 | WebSearch | Antoniadou resonant periodic orbits restricted problem 2026 | 2018-2019 papers (held anchor); arXiv listing later showed 2609.24405 (white-dwarf dust, unrelated) |
| 13 | WebSearch | Uranian satellites trajectory design three-body Titania Oberon Umbriel transfer | 2110.03683 (held), Uranian-system dynamical-evolution papers (2403.17896/7, 2509.24631, 2005.12887) |
| 14 | WebSearch | AAS/AIAA ASC 2025 Boston programme, moon tour Uranus Neptune | Programme landing page; full PDF later fetched (query 29) |
| 15 | WebSearch | AAS/AIAA SFMM 2026 programme resonant orbits moons | ASC 2026 programme PDF (fetched, query 28), SFMM 2026 call for papers |
| 16 | WebSearch | Bhanu Kumar Rice Anderson torus connections four-body 2026 | Kumar publication page (fetched, query 27) |
| 17 | WebSearch | arXiv 2026 Uranus moons periodic orbits patched three-body ballistic flybys | 2509.12671 (held), 2603.21750 (Uranian moons fragility, dynamical), 2603.07085 (held) |
| 18 | WebSearch | Neptune Triton heteroclinic homoclinic resonant CR3BP arXiv 2026 | Computer-assisted proofs math/0201278, math/0401146 (old); nothing 2026 |
| 19 | WebSearch | Spear Bosanac Stuart Neptune-Triton resonant orbit tour | Miceli-Bosanac-Stuart SciTech 2024 (held), JAS 2026 (held) |
| 20 | WebSearch | cycler trajectory arXiv 2026 astro-ph.EP cyclers | 2606.29189 (held), 2605.31543 (held), 2606.26367, 2111.11858 (asteroid-flyby cyclers, 2022) |
| 21 | WebSearch | JGCD 2026 cycler trajectories moon system | Russell-Strange 2012, Casoliva 2010 (held), Roberts-Tsoukkas journal draft (held anchor) |
| 22 | WebSearch | CMDA 2026 heteroclinic invariant tori restricted four-body | 10.1007/s10569-024-10203-5 (Hill R4BP tori, 2024), 2507.06123, older four-body papers |
| 23 | WebSearch | Acta Astronautica 2026 Uranus satellite tour gravity assist | 2603.07085 (held), 2609.02189 (Halley, unrelated), Heaton 2003 (held) |
| 24 | WebSearch | Hernandez Jones Jesick cycler 2025 2026 Jovian Saturnian | Only 2017 AAS paper and JGCD 2025 Liang (both held) |
| 25 | WebSearch | "Neural Keplerian Maps" resonance hopping Earth-Moon | CMDA 10.1007/s10569-026-10318-x (Keplerian-map transport barriers; Springer redirected, abstract not read), 2505.10138, ESA ACT inter-moon page |
| 26 | WebSearch | IAC 2025 Sydney astrodynamics cycler moons Uranus tour | JAS 10.1007/s40295-025-00565-9 (Saturnian pump-down tours, directed graphs), 2505.05514; IAC archive has no title-level index reachable |
| 27 | WebFetch | bhanukumar314.github.io | 13 publications listed (see section 4); only Uranus item is the 2024 Oberon survey (held) |
| 28 | WebFetch | ASC 2026 programme PDF (space-flight.org), text-extracted | Titles only (no abstracts); relevant titles in section 2 |
| 29 | WebFetch | ASC 2025 full programme PDF, text-extracted | Abstracts; AAS 25-668, 25-669, 25-621, 25-569, 25-677, 25-523, 25-759; nothing on Uranus or Neptune moons |
| 30 | WebSearch | AIAA SciTech 2026 Orlando SFMM resonant moon tour Neptune Uranus | Programme PDF link only (scitech.aiaa.org ST26 programme, not fetched); 2606.08485 |
| 31 | WebSearch | Astrodynamics journal 2026 cycler / resonant / moon tour | 10.1007/s40295-025-00509-3 (Earth-Moon resonant orbits for cislunar access), 2606.08485, 2504.12470, Springer chapter 10.1007/978-3-032-20347-2_1 |
| 32 | WebSearch | exact title "Ballistic Cycler Trajectory between EML2 and SEL2 in the BCR4BP" | No paper found; only general BCR4BP literature |
| 33 | WebSearch | exact title "Physics-Guided Diffusion for Pseudospectral Moon-to-Moon Tour Design" | No paper found |
| 34 | arXiv API | cycler AND (moon OR Uranus OR Jupiter OR Saturn), newest first | 2 hits: 2606.29189, 2605.31543 (both held) |
| 35 | arXiv API | (Uranus OR Triton OR Neptune) AND (trajectory OR orbits OR resonant) | Dominated by exoplanet and chemistry hits; trajectory-relevant: 2609.21714 (Naiad survival, dynamical evolution), 2609.31580 (Miranda inclination), 2511.08698 (Uranus gravity-field orbits) |
| 36 | arXiv API | (torus OR tori) AND (heteroclinic OR homoclinic) AND restricted | Only older items; 2109.14814 (held) is the newest relevant |
| 37 | arXiv API | Venus AND Mars AND cycler | 0 results |
| 38 | arXiv API | Triton AND (three-body OR trajectory OR spacecraft) | All nuclear or software hits; no Neptune-moon trajectory paper |
| 39 | arXiv API | Uranus AND (gravity assist OR flyby OR satellite tour OR three-body) | 2505.05514, 2511.08698; no Uranian-moon tour or cycler paper |
| 40 | arXiv API | "restricted three-body" AND (resonant OR heteroclinic OR cycler OR periodic orbits), 30 newest | Listed in section 2 and 3 (2604.00679, 2606.08485, 2602.16354, 2512.03849, 2507.07940, 2507.04739 and others) |
| 41 | arXiv API | "four-body" AND tori AND (moon OR Jupiter OR Ganymede) | 1 hit, 2010 (irrelevant) |
| 42 | arXiv API | (Ganymede OR Europa OR Callisto OR Titan OR Enceladus) AND trajectory terms, astro-ph.EP | 2607.03505, 2608.04226, 2603.07085 (held); no cycler |
| 43 | arXiv API | parameterization method / whiskered / invariant tori AND heteroclinic AND three/four-body | 2509.03655 (held), 2507.06123 (not held), 2607.28452 and 2607.28472 (non-autonomous KAM theory, pure maths) |
| 44 | arXiv API | gravity assist / resonance hopping / moon tour AND moons | 2603.07085 (held), 2308.10029 (Canales-Howell-Fantino 2023), 2210.14996 (Takubo-Landau-Anderson, Saturn tours) |
| 45 | arXiv API | authors Bosanac, Baresi, Anderson, de la Llave + orbit terms | No usable results (name collision with unrelated computing papers) |
| 46 | Crossref | "cycler trajectory", 2026 journal articles | No relevant hits (query returned unrelated engineering papers) |
| 47 | Crossref | Uranus satellite tour OR Uranian moons trajectory, since 2025-06 | Barnes and do Vale Pereira (Aerospace 2025), Zhang-Li-Baoyin TAES 2025 (held), Ivanyukhin (Jupiter low-thrust tour, Cosmic Research 2025), Gomes-Keizer (Uranian dynamics III), Clement et al. (Icarus 2026) |
| 48 | Crossref | Neptune Triton periodic orbits resonant three-body, since 2025 | Generic CR3BP periodic-orbit items (section 3); nothing at Neptune-Triton |
| 49 | Crossref | CCR4BP invariant tori connections moons | HTTP 429 (rate limit); no result |
| 50 | WebSearch | Advances in Space Research 2026 Earth-Mars cycler / semicycler | Only 2000s-2015 Earth-Mars cycler literature (held) |
| 51 | WebSearch | JAS 2026 Neptune Triton OR Uranus OR Titania OR Oberon | Miceli-Bosanac JAS 2026 (held); no Uranian-moon item |
| 52 | WebSearch | CNSNS 2026 heteroclinic tori three-body spacecraft transfer | Bonasera-Bosanac 2023 JGCD tori transitions, 10.2514/1.G009219 (Gateway tori orbit capacity), ASR 10.1016/S0094-5765(26)00327-9 review |
| 53 | WebSearch | Kumar "Titan and Rhea" AAS 26 abstract | Abstract not indexed; title and author confirmed only from the programme |
| 54 | WebSearch | Callisto Ganymede Europa cycler Yang Liang 2025 2026 | Liang-Yang-Li JGCD 2025 (held) only |
| 55 | WebSearch | Mars Earth Venus ballistic cycler 2026 Aldrin new family | Wikipedia, project GitHub, Sorensen blog on VanderVeen 1969 EVMVE |
| 56 | WebFetch | selenianboondocks.com 2026/04 EVMVE-2034 | Blog by Kirk Sorensen (27 Apr 2026) on a one-off 750-day EVMVE ballistic flyby after VanderVeen 1969; stated by the post not to be a cycler |
| 57 | WebFetch | arXiv abstract pages 2606.29189, 2606.26367, 2505.10138, 2507.06123, 2505.05514, 2604.00679, 2606.08485 | Abstracts read (see section 2) |
| 58 | WebFetch | Crossref records for Barnes, Nagai, Gao et al. | Metadata only; abstracts absent |

58 rows; of these, queries 1-26, 29-33, 50-55 are web searches (about 38) and 34-48, 57 are database or listing queries. More than 25 differently worded queries were run.

## 2. NEW items (not in CORPUS_INDEX.md)

None of the items below is rated COLLISION RISK. Ratings are of closeness to a protected result, as the brief defines them.

### 2.1 SHOULD HOLD

1. Rosengren A. J., Rawat A., Kumar B., Ross S. D., "The Astrodynamics Primer on Cislunar and Translunar Space", arXiv:2606.26367, submitted 24 Jun 2026 (77 pp, over 300 references). Abstract: a unified spatial description of the Earth-Moon region: inner cislunar zone dominated by secular effects, an outer zone structured by lunar resonances, circumlunar space organised by gateway geometry, and translunar space. Area 6. Survey-level; consolidates the Kumar-Rawat-Rosengren-Ross resonance papers that the corpus holds individually.
2. Rawat A., Kumar B., Rosengren A. J., Ross S. D., "Cislunar Mean-Motion Resonances: Definitions, Widths, and Comparisons with Resonant Satellites", arXiv:2505.10138, 15 May 2025; accepted JGCD, DOI 10.2514/1.G009336. Abstract: defines the resonance zone through the separatrix of unstable resonant periodic orbits for the 2:1 and 3:1 resonances and finds regions of influence broader than semi-analytical predictions. Area 6. The corpus holds the neighbouring 2024-2025 Rawat/Kumar papers but not this one.
3. Fernandez-Mora A., Haro A., de la Llave R., Mondelo J.-M., "Simultaneous computation of whiskered tori and their whiskers in Hamiltonian systems using flow maps", arXiv:2507.06123, 8 Jul 2025. Abstract: a flow-map algorithm that computes partially hyperbolic invariant tori and high-order expansions of their stable and unstable manifolds together, tested in the CR3BP. Area 3 (method; de la Llave co-author; CR3BP only, no moon or four-body system).
4. Guido A. F., Efthymiopoulos C., "Arches of chaos, heteroclinic connections of first-order MMRs and the chaotic transport of small bodies in the Sun-Jupiter system", arXiv:2604.00679, 1 Apr 2026. Abstract: examines "the heteroclinic connections between stable and unstable manifolds of unstable periodic orbits" of mean-motion resonances in the planar Sun-Jupiter restricted three-body problem, and reports heteroclinic links between L3 short-period orbits and resonant orbits and between interior and exterior resonance structures. Areas 2 and 6 (method: resonant-orbit heteroclinics; Sun-Jupiter mass ratio, small-body application).
5. Park B., Howell K. C., "Linking Averaged and Unaveraged Three-Body Dynamics Near Smaller Primaries: Symmetric Periodic Orbits", arXiv:2606.08485, 7 Jun 2026. Abstract: links equilibria of averaged models to symmetric periodic orbits in the Hill and circular restricted problems and builds bifurcation diagrams described as "a comprehensive atlas of the symmetric periodic orbit web". Areas 2 and 6 (symmetric periodic orbit families near the smaller primary; no named Neptune-Triton or Uranian system in the abstract).
6. Simon A., Cohen I., Hedman M., Hofstadter M., Mandt K., Nimmo F., "Uranus Flagship Science-Driven Tour Design: Community Input Poll", arXiv:2505.05514, 7 May 2025. Abstract: a community poll to set tour-parameter priorities for the Uranus Orbiter and Probe, with data on Zenodo. Area 1 (tour requirements; no trajectory result). The corpus holds the 2026 Simon et al. PSJ concept update, not this poll.
7. Kumar B., "The Effect of Untargeted Moons on Resonant Orbit Structure Between Titan and Rhea", AAS/AIAA ASC 2026 (Whistler, 26-30 Jul 2026), Paper 1023. Title and author only; no abstract or preprint located. Area 3 (same author as the held Oberon CCR4BP survey; Saturn, not Uranus; whether it uses a four-body model is not stated in the title).
8. Bonasera S., Bosanac N., "Computing Natural Transitions Between Tori Near Resonances in the Earth-Moon System", JGCD 2023, DOI 10.2514/1.G006941. Abstract: Poincare maps plus manifold learning to find natural transitions between families of invariant 2-tori near resonances in the Earth-Moon CR3BP. Area 3 and 6 (torus-to-torus transfers; three-body only).

### 2.2 NOTE ONLY

- Komachi S. (ANU), "Ballistic Cycler Trajectory between EML2 and SEL2 in the Bicircular Restricted Four-Body Problem", ASC 2026 Paper 688. Title only; no abstract located. Area 6 (a cycler in a four-body model, Earth-Moon/Sun-Earth; not a moon-system cycler).
- Okada H. (Univ. Tokyo), "Physics-Guided Diffusion for Pseudospectral Moon-to-Moon Tour Design in the CR3BP", ASC 2026 Paper 854. Title only. Areas 1 and 5 (moon-to-moon tour method; system not stated).
- Gomez R. (Purdue), "Anatomical Analysis of Orbit-Chain Trajectories in the Earth-Moon System", ASC 2026 Paper 874. Title only. Area 6.
- Henry D. (Univ. Minnesota), "Heteroclinic Connections in the Elliptic Restricted Three-Body Problem", ASC 2026 Paper 936. Title only. Area 6 (method).
- Takao Y. (Yokohama Nat. Univ.), "Asteroid Flyby Cycler Trajectories via Solar Sailing", ASC 2026 Paper 1012. Title only. Heliocentric flyby cycler; outside the six areas.
- Grossi G. (Politecnico di Milano), "Quasi-Periodic Invariant Manifolds in the In-Plane Quasi-Hill Restricted Four-Body Problem", ASC 2026, Paper 622. Title only. Area 3 (four-body tori; Hill-type model, no moon system stated).
- Sailor N., Santacesaria M. et al., AAS 25-630, "Leveraging Low-Thrust Optimization for Tori Intersection Transfers in the Earth-Moon System" (ASC 2025); Ito S., Baresi N., Izzo D., AAS 25-523, "Neural Keplerian Maps for Resonance Hopping Transport in the Earth-Moon System" (ASC 2025; abstract not read; Baresi is on the protected-author list, three-body Earth-Moon only); "Identifying cislunar transport barriers with Keplerian maps", CMDA, DOI 10.1007/s10569-026-10318-x (abstract not read). Area 6, methods.
- Rodriguez A. et al., AAS 25-759, BEACON mission concept using a "three-petal Earth-Moon cycler orbit" (ASC 2025). Area 6, mission concept.
- Englander J., Ellison D. et al., AAS 25-669, "Robust Access to Uranus via Solar Electric Propulsion" (ASC 2025). Area 1 (interplanetary leg only).
- Barnes D., do Vale Pereira P., "Preliminary Proof of the Feasibility of a Novel Mission Concept and Spacecraft Trajectory for Exploring Uranus with Small Satellites", Aerospace 12(12):1069, 30 Nov 2025, DOI 10.3390/aerospace12121069. Jupiter-Uranus gravity assist, CubeSat constellation; no moon tour. Area 1 (adjacent).
- Mankovich C. R. et al., arXiv:2511.08698, 11 Nov 2025: Uranus gravity-field science from close orbits of a UOP. Area 1 (adjacent; planet, not moons).
- Ivanyukhin A. V., "An Analysis of Low-Thrust Satellite Tour Strategy in the Jupiter System", Cosmic Research, Oct 2025, DOI 10.1134/s0010952525601318. Area 5 (low-thrust tour, not a cycler).
- Gomes S. R. A., Keizer T., "Dynamical evolution of the Uranian satellite system III. The passage through the 7/4 MMR between Miranda and Ariel", Icarus, May 2026, DOI 10.1016/j.icarus.2026.116974; Clement M. S. et al., "The fragility of the Uranian moons during the giant planet instability", Icarus, Jul 2026, DOI 10.1016/j.icarus.2026.117056 (arXiv:2603.21750); Agrusa H. et al., arXiv:2609.21714 (Naiad); El Moutamid and Cuk, arXiv:2609.31580 (Miranda inclination); Lari and Rossi, arXiv:2607.03505 (Ganymede-Callisto 7:3). Natural-system dynamics, no trajectory design.
- Aydin C., arXiv:2602.16354 (comet-type periodic motions, Earth-Moon CR3BP) and arXiv:2508.17286 (DRO vertical self-resonant bifurcations); Joung C., Koh D., van Koert O., arXiv:2512.03849 (highly inclined near halo bifurcations, includes Saturn-Enceladus mass ratio); Gao C. et al., Nonlinear Dynamics 114(10), May 2026, DOI 10.1007/s11071-026-12604-7 (cislunar resonant periodic orbits in the ER3BP); Nagai Y., "Rose-like periodic orbits in the restricted three-body problem", CMDA 138(5):52, DOI 10.1007/s10569-026-10327-w (abstract not available). Periodic-orbit family work, Areas 2 and 6 (method, adjacent).
- Sorensen K., blog post "Manned Circumnavigation of the Inner Solar System", 27 Apr 2026: EVMVE one-off ballistic flyby (departure 4 Aug 2034) after VanderVeen 1969; the post states it is not a cycler. Area 4 (adjacent; blog, not a publication).
- Orbit Capacity of Higher-Dimensional Tori Around Gateway, JGCD, DOI 10.2514/1.G009219 (online 9 Feb 2026); The Linear Parametrization of Two Dimensional Tori Families in the Restricted Three Body Problem, JAS 2026 (title seen in a search snippet only); 2504.12470 (frequency-domain differential corrector for quasi-periodic trajectories); Springer chapter 10.1007/978-3-032-20347-2_1 (Multi-Body Dynamics and Periodic Orbits); JAS 10.1007/s40295-025-00565-9 (Saturnian pump-down tours, directed graphs); JAS 10.1007/s40295-025-00509-3 (Earth-Moon resonant orbits for cislunar access). Titles and landing pages only; Springer and AIAA pages were not opened past their access redirect.

### 2.3 COLLISION RISK

None found. Specifically, no item located in this pass: (a) describes cyclers, repeating or periodic moon-to-moon trajectories or resonant inter-moon transfers among Ariel, Umbriel, Titania and Oberon; (b) treats homoclinic or heteroclinic connections or symmetric periodic orbits near the 4:5 resonant orbit of the Neptune-Triton problem; (c) computes connections between tori in a four-body model at Uranus beyond the held Oberon survey and Titania work; (d) enumerates new VEM or Earth-Mars cycler classes. This is conditional on the channels reached (section 5): the ASC 2026 programme was read at title level only, and the Kumar Titan-Rhea paper (Paper 1023) has no abstract available.

## 3. Items checked and already held

2606.29189 (Ross, Roberts-Tsoukkas, stable ballistic prograde cyclers; its abstract states families across mass ratios "from the Sun-Jupiter regime to the equal-mass limit"); 2605.31543 (Braik-Ross); 2509.12675 and its ASR DOI (Kumar-Rawat-Rosengren-Ross); AAS 25-621 (Ross-Roberts-Tsoukkas); AAS 25-569 and 25-677 (Kumar et al.); 2601.00149 and 2509.03655 (Kumar); 2109.14814, 2109.14815, 2309.06073, 2105.11100 (Kumar-Anderson-de la Llave line); Kumar-Anderson Oberon survey AAS 24-288; Owen-Baresi journal and AAS 23-110 knot-theory papers; Owen-Baresi-Scheeres tri-circular ISSFD 2024; Miceli-Bosanac motion primitives JAS 2026 and Spear 2021; Canales-Howell-Fantino 2110.03683; Pozzi et al. 2603.07085; Simon et al. PSJ ae680c; Ellison et al. AAS 25-668; Zhang-Li-Baoyin TAES 2025 (10.1109/taes.2025.3618517); Liang-Yang-Li JGCD 2025; Hernandez-Jones-Jesick 2017; Russell-Strange 2009 and the Pony Express papers; 2509.12671 (Zhou et al. fixed points); Roberts-Tsoukkas and Casoliva anchors; Antoniadou-Libert 2019; Trajectory Options for a Uranus Orbiter and Probe (ResearchGate listing).

## 4. Author pages followed

bhanukumar314.github.io lists 13 items; all were already held except 2505.10138 / 10.2514/1.G009336 (section 2.1, item 2) and the ASC 2026 Titan-Rhea paper (item 7, not yet on the page). It lists no Uranian work after the 2024 Oberon survey. Other author pages (Bosanac, Anderson, Baresi, Owen, Longuski) were not located as separate listings; the arXiv author query for them returned only unrelated computing papers (query 45).

## 5. Gaps (not accessed or not completed)

- AIAA SciTech 2026 and the 2026 Space Flight Mechanics Meeting: only the call for papers and a programme PDF link were reached; no paper-title list was read.
- IAC 2025 and IAC 2026: the paper archive (iafastro.directory) offered no title-level index through the tools used; no IAC astrodynamics session list was read. IAC 2026 has not been searched separately.
- ASC 2026 programme: titles only; abstracts were not available. Kumar Paper 1023, Komachi Paper 688, Okada Paper 854, Gomez Paper 874, Henry Paper 936 and the Quasi-Hill R4BP paper could not be read.
- ASC 2025 programme was read in full text for the abstract lines matched by keyword only (Uranus, Neptune, Triton, cycler, torus, tori, four-body, heteroclinic, moon names); a Neptune or Uranus abstract using none of those words would be missed.
- Springer (JAS, CMDA, Astrodynamics), AIAA ARC (JGCD), ScienceDirect (Acta Astronautica, ASR, CNSNS): landing pages redirect to a login or consent handshake; no workaround was attempted. Online-first lists of these journals were therefore examined only through Crossref and search snippets. Crossref returned HTTP 429 once (query 49) and was not retried.
- Journal abstracts not read: 10.1007/s10569-026-10318-x, 10.1007/s40295-025-00565-9, 10.1007/s40295-025-00509-3, 10.1007/s10569-026-10327-w, 10.1007/s11071-026-12604-7, 10.1109/taes.2025.3618517 (held, not re-read).
- Google Scholar, ADS, NTRS and the Zenodo record for the Uranus poll were not queried.
- The arXiv API was used instead of the astro-ph.EP, math.DS and physics.space-ph listing pages; physics.space-ph was not queried by category.

## Coordinator's check (2026-10-03)

Checked against `CORPUS_INDEX.md` after the pass:

- Two "should hold" items are ALREADY HELD and are not new: Rawat, Kumar, Rosengren & Ross,
  "Cislunar Mean-Motion Resonances" (held as the JGCD 49(4) version) and Fernandez-Mora, Haro,
  de la Llave & Mondelo (whiskered tori by flow maps).
- Not held, and worth obtaining in this order: Kumar, "The Effect of Untargeted Moons on Resonant
  Orbit Structure Between Titan and Rhea" (ASC 2026 Paper 1023; title only so far; a four-body
  moon-system study by the author of the Oberon work, so it bears on `#882`); Guido & Efthymiopoulos
  (arXiv:2604.00679); Park & Howell (arXiv:2606.08485); Rosengren et al. primer (arXiv:2606.26367);
  Bonasera & Bosanac (DOI 10.2514/1.G006941); Simon et al. community poll (arXiv:2505.05514).
- One ASC 2026 title is a CYCLER and should be read when a paper exists: Komachi, "Ballistic Cycler
  Trajectory between EML2 and SEL2 in the Bicircular Restricted Four-Body Problem" (Paper 688). It
  is an Earth-Moon/Sun-Earth object, not a moon-system cycler, so it does not bear on a catalogued
  row, but it is in the catalogue's scope.
- The pass could not read publisher landing pages or the conference abstracts, so "no collision
  found" is conditional on titles and arXiv abstracts only.

Next pass due about 2026-10-17.
