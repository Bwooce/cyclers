# Digest — Yang, Hu, Bai & Li (2023), "Review of Trajectory Design and Optimization for Jovian System Exploration" (Space: Science & Technology 3, 0036)

**Digested:** 2026-10-03 (text-layer PDF, no OCR needed; read in full, all 15 pages).
**Why digested:** an adjudication of the project's 36 two-moon Galilean "symmetric closure" trajectories cited this review as surveying "the whole double-cycler field" before the paper was in hand, and project notes call it "the 2024 Jovian review". This digest lets those two claims be checked against the text. It does not judge whether any project trajectory is new.

## Citation (as printed)
Hongwei Yang, Jincheng Hu, Xiaoli Bai, and Shuang Li. "Review of Trajectory Design and Optimization for Jovian System Exploration." *Space: Science & Technology* 2023;3: Article 0036. https://doi.org/10.34133/space.0036. REVIEW ARTICLE. Affiliations: College of Astronautics, Nanjing University of Aeronautics and Astronautics (Yang, Hu, Li); Rutgers, The State University of New Jersey (Bai).
**Submitted 6 February 2023; Accepted 5 May 2023; Published 23 May 2023** (first page). Copyright 2023, CC BY 4.0. Journal ISSN 2692-7659 (last page). The year printed everywhere is 2023, not 2024.
Filed in the private paper corpus as `yang-hu-bai-li-2023-review-trajectory-design-optimization-jovian-system-exploration-space-sci-technol-3-0036-doi-10.34133-space.0036.pdf` (md5 `002b083a6699cd0174465ef9c2af3116`, 15 pages). Page numbers below are the printed journal page numbers (PDF page 1 = printed page 1).

## What it is
A literature review. Abstract: "This review provides a systematic summarization of the past and state-of-art methodologies for 4 main exploration phases, including Jupiter capture, the tour of the Galilean moons, Jupiter global mapping, and orbiting around and landing on a target moon." Its stated contribution (p.2): "to provide a systematic summarization of the past and state-of-the-art methodologies of trajectory design and optimization for the Jovian system exploration". It is a methods review, not a cycler survey; cyclers are one bullet inside one subsection.

Sections (with printed page): Introduction (1); Jupiter Capture Trajectories (2) — Satellite-aided capture; Jupiter capture using electrodynamic tether propulsion; Solar electric propulsion trajectories for Jupiter capture; Cloudtops arrivals; Jupiter capture trajectory design and optimization; Jupiter capture navigation analysis; Tour of Galilean moons (4) — Two-body techniques for trajectory design [1. Special flyby sequences (resonant hopping, petal rotation, COT, switch-flip, **Cyclers**), 2. V-infinity leveraging transfer, 3. Graphic methods (Tisserand graph; (V-Infinity, Resonance) graph)] (5), Three-body techniques for trajectory design [1. T-P graph, Flyby map, TILT; 2. Invariant manifolds patching; 3. AI-based gravity-assist models] (6-7), High-fidelity trajectory design (7), Multiple-flyby trajectory optimization [flyby sequence optimization; impulsive and continuous optimization] (7-8), Robust trajectory design and optimization (8); Jupiter Global Mapping Trajectories (8); Moon Orbiter and Lander Trajectories (9) — Orbits around Galilean moons; Orbit capture at Galilean moons; Landing on Galilean moons; Conclusions and Future Development (10). Table 1 (missions), Table 2 (multiple-satellite-aided captures), Figs 1-7 (schematics).

## CYCLERS — every passage

The word "cycler" occurs in exactly two places in the body text, both on p.5, both in Section "Tour of Galilean moons", subsection "Two-body techniques for trajectory design", item "1. Special flyby sequences". It also occurs in two reference titles ([31], [32]) and one cited reference [30]. No other occurrence (checked by text search of the whole PDF, including the abstract, Introduction, Table 2, Conclusions).

**Passage A (p.5, left column), a Galileo-era-to-Europa-Clipper tour-design usage, not a free-standing cycler study:**
"For instance, Campagnola et al. [7] utilized Ganymede–Callisto cyclers, a Callisto π-transfer sequence with 1:1 resonant transfers and Callisto petal rotation to achieve fast and low-total-ionizing-dose rotations in the tour design of Europa Clipper. The special flyby sequences briefly introduced herein include resonant hopping, petal rotation, crank-over-the-top (COT) sequences, switch-flip, and Cyclers."
- Class: two-moon (Ganymede-Callisto), used as a building block in a Europa Clipper tour.
- Cited: [7] only. No cycler details given (no period, V-infinity or flyby count).

**Passage B (p.5, right column), the "Cyclers" bullet, quoted in full:**
"Cycler trajectories allow a spacecraft to visit planetary moons repeatedly using little or no fuel. Extensive cyclers exist in the Jupiter-moons system. Russell and Strange [30] first investigated the double-cyclers for exploring the Galilean moons. Specifically, the cyclers of Ganymede–Io, Ganymede–Europa, Ganymede–Callisto, and Europa–Ganymede are searched. The equations for calculation of nonresonant free-return are derived from two-body dynamics. Then, a broad search in an ideal model consisting of circular and coplanar celestial body orbits is performed to find hundreds of ideal model ballistic cycler geometries [30]. In addition to the double-cyclers, triple-cyclers also exist due to the Laplace resonance of Io, Europa, and Ganymede. Three types of Laplace-resonance triple-cyclers are constructed by Lynam and Longuski [31] according to the orbital period of the cycler compared with the period of the Laplace resonance, including half-period, single-period, and multiple-period cyclers. Due to the large phase space in searching for the triple-cyclers, Hernandez et al. [32] further developed a good initial guess strategy. In their initial guess search, a desired trajectory is approximated by a two-body conic and classified by 3 geometrical parameters, which are period, eccentricity, and orientation with respect to an inertial frame. This approach reduces the search complexity through constraining the energy of the cycler families."

Breakdown of Passage B:
| Class | Bodies named | What is said | Reference |
|---|---|---|---|
| double-cyclers (two-moon) | Ganymede-Io, Ganymede-Europa, Ganymede-Callisto, Europa-Ganymede (the four pairs; the paper states "cyclers of ... are searched") | free-return equations from two-body dynamics; broad search in circular-coplanar ideal model; "hundreds of ideal model ballistic cycler geometries" | [30] |
| triple-cyclers (three-moon) | Io, Europa, Ganymede (Laplace resonance) | three types by cycler period versus Laplace-resonance period: half-period, single-period, multiple-period | [31] |
| triple-cyclers, search strategy | triple-cyclers (moons implicit from [32]'s title: Io-Europa-Ganymede) | good-initial-guess strategy; conic classified by period, eccentricity, orientation; energy constraint on cycler families | [32] |

(One apparent slip in the paper's wording: "Hernandez et al. [32]" while the printed reference lists Hernandez, Jones, Jesick.)

## Explicit answers

**(a) Does it describe or name "double cyclers" or two-moon cyclers, and for which pairs?** Yes, in one sentence of attribution to Russell and Strange [30]: "Russell and Strange [30] first investigated the double-cyclers for exploring the Galilean moons. Specifically, the cyclers of Ganymede–Io, Ganymede–Europa, Ganymede–Callisto, and Europa–Ganymede are searched." Plus Ganymede-Callisto cyclers used in Europa Clipper tour design via [7] (p.5). It adds nothing beyond [30]'s own pair list and gives no results from them.

**(b) Any cycler involving Io and Callisto together, or Europa and Callisto together?** Not stated. No cycler is named for Io-Callisto or Europa-Callisto, and no statement that such cyclers do or do not exist. Callisto appears in cyclers only in Ganymede-Callisto (Passages A and B). Caution against a false hit: Table 2 (p.3) and the capture text list the moon combinations IC and EC (and CEI, CGI, CGE, CGEI), but those are multiple-satellite-aided Jupiter-capture sequences (Table 2 caption: "Summary of multiple-satellite-aided captures in Refs. [12,14]"), not cyclers. Quote from the table: "D | IE, IG, IC, EG, EC, GC | 24 | ...". The text says of Callisto: "Callisto is not in Laplace resonance but near-resonances can be found".

**(c) Does it claim to survey all cycler work, or only mention some?** It does not claim to survey cycler work as a whole; it gives a single short paragraph (Passage B) citing three works ([30], [31], [32]) plus a use in [7]. The review's own scope claim is methods across four mission phases (abstract), and its own admission that coverage of earlier reviews is partial is about other reviews ("only a small part ... are described in Refs [8-10]", p.2). Nothing in it says the cycler list is complete or exhaustive. The sentence "Extensive cyclers exist in the Jupiter-moons system" is the only statement of extent, and it is a statement about existence, not about survey coverage.

**(d) Any table or numbers for cyclers (periods, V-infinity, flyby altitudes, flyby counts) usable as reproduction targets?** No. There is no cycler table or figure. The only number attached to cyclers is "hundreds of ideal model ballistic cycler geometries" (from [30]). Nothing else is printed (no periods, V-infinity, altitudes or flyby counts); for Russell-Strange numbers use the project's own digest of [30]. Nothing to transcribe. Table 2 is about capture sequences (counts D=24, T=56, Q=16 sequences), not cyclers.

## Other trajectory classes surveyed (anchoring references)
- **Satellite-aided Jupiter capture** (pp.2-4, Table 2): single/double/triple/quadruple-satellite-aided capture. Cline [13]; Lynam et al. [14] (24 double sequences, 56 triple, 16 quadruple; Io-Ganymede lowest delta-V, Ganymede-Callisto best for perijove > 8 RJ); broad-search by Lynam [12,16-18]; Lynam and Longuski [15,19]; Macdonald and McInnes [20]; Scott et al. [11]; Galileo saved "about 175 m/s [11,12]".
- **Electrodynamic tether capture**: Lantoine et al. [21] (MAGNETOUR), Schadegg et al. [22].
- **SEP capture**: Strange et al. [23], McElrath et al. [24], Landau et al. [25], Patrick and Lynam [28]; **cloudtops arrivals** McElrath et al. [24] ("over 500 m/s can be saved").
- **Capture optimisation and navigation**: Labroquere et al. [26], Izzo et al. [27], Lynam and Longuski [19].
- **Patched-conic tour techniques** (pp.5-6): resonant hopping, petal rotation (Anderson et al. [29]), COT, switch-flip, all via Campagnola et al. [7]; **V-infinity leveraging transfers** (Campagnola and Russell [33], Campagnola, Strange and Russell [34], Lantukh et al. [35,36]); **Tisserand graph** (Strange and Longuski [37]); **(V-infinity, resonance) graph** ([7]).
- **Three-body techniques** (pp.6-7): Tisserand-Poincare graph (Campagnola and Russell [38]); flyby map ([39,40] on Keplerian map [41]; spatial [42]); Tisserand-leveraging transfers [43]; invariant-manifold patching (Fantino and Castelli [44]; moon-to-moon analytical transfer, Canales et al. [45-48]; resonant-orbit manifolds, Anderson [49]; concentric circular restricted four-body problem and whiskered tori, Kumar et al. [50,51]); AI/neural-network gravity-assist models ([52,53]).
- **High-fidelity conversion**: continuation method Bradley and Russell [54]; [55,56].
- **Multiple-flyby optimisation** (pp.7-8): enumeration/Bellman [57]; Monte Carlo [58]; branch and bound [33]; GTOC6 tree search [59,60]; Keplerian-map approach [61]; impulsive [62,39]; continuous, GTOC6 winner [63]; [21]; two-loop [64].
- **Robust design** (p.8): Lam et al. [67], Park and Scheeres [68], PCE methods [69,70], Greco et al. [71].
- **Jupiter global mapping** (p.8-9): 3-D Tisserand graph [72]; repeating-ground-track orbits, GTOC10 [73,74]; Jiang et al. [75]; J2-Lambert DNN [76].
- **Orbits, capture and landing at moons** (pp.9-10): Europa science and frozen orbits [77-81], Ganymede orbits [82,83], low-energy and distant orbits [84-92], tether orbits [93], approach from resonance [94-96], temporary and tight/loose capture at Europa [97-101], landing [102-106].

## What it does NOT contain
- No cycler table, figure, period, V-infinity, flyby altitude or flyby count.
- No mention of any Io-Callisto or Europa-Callisto cycler (see 4b); no three-moon cycler involving Callisto.
- No statement that the cycler literature is complete; no statement that no other cyclers exist.
- No cycler results of its own; no new numerical computations anywhere (it is a review).
- No description of how [30]'s cyclers were found beyond two sentences, and no symmetric-closure or two-moon-closure terminology ("closure" does not occur anywhere in the text; "symmetric" occurs once, in "axisymmetric and doubly symmetric periodic orbits" near Europa, p.9, ref. [91], unrelated to cyclers).
- No cyclers at Saturn or other planets, apart from reference [86] Saturnian orbiters, which is not about cyclers.

## Reference list: every cycler-related entry (as printed, reference number as in the paper)
- **[7]** Campagnola S, Buffington BB, Lam T, Petropoulos AE, Pellegrini E. Tour design techniques for the Europa Clipper mission. J Guid Control Dyn. 2019;42(12):2615-2626. (Ganymede-Callisto cyclers used in the Europa Clipper tour, p.5.)
- **[30]** Russell RP, Strange NJ. Cycler trajectories in planetary moon systems. J Guid Control Dyn. 2009;32(1):143-157. (Double-cyclers: Ganymede-Io, Ganymede-Europa, Ganymede-Callisto, Europa-Ganymede.)
- **[31]** Lynam AE, Longuski JM. Laplace-resonant triple-cyclers for missions to Jupiter. Acta Astronaut. 2011;69(3-4):158-167.
- **[32]** Hernandez S, Jones DR, Jesick M. Families of Io-Europa-Ganymede triple cyclers. Paper presented at: 2017 AAS/AIAA Astrodynamics Specialist Conference; 2017 Aug 20-24; Stevenson, WA.

No reference in the list has the word "Liang" as an author (text search of the whole document); the paper does not cite any "Liang et al." triple-cycler work. Reference entries near the cycler topic that are about capture, not cyclers: [12], [14]-[19] (Lynam multiple-satellite-aided capture series).

## Status
Filed; to be indexed in CORPUS_INDEX by the caller (this digest does not edit it). Net effect for the project's claims, from the text alone: the paper names the four [30] pairs and the Io-Europa-Ganymede Laplace triple-cyclers [31,32] in a single paragraph, states no Io-Callisto or Europa-Callisto cycler, makes no claim to survey the cycler field, and was published in 2023.
