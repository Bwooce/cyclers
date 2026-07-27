# #730 — Consolidated citation-mining acquisition backlog (master list)

**Update 2026-07-27 (`#732`):** the user directly supplied PDFs for the top
three items of §2's ranking — **item 1 (Blazevski & Ocampo 2012), item 4
(Baresi, Olikara & Scheeres 2018), and item 5 (Negri & Prado 2020)**. All
three are now ACQUIRED — filed in the private `cyclers_pdf` corpus, digested
(including each paper's own mandatory cross-check question against this
project's code), citation-mined, and registered in `CORPUS_INDEX.md`. See
`docs/notes/2026-07-27-732-blazevski-negri-baresi-foundational-papers-digest.md`.
Struck through in §2's table below; **no longer live backlog** — do not
re-flag or re-acquire these three.

**Status: this is a compilation/reference document, NOT a to-do list to action immediately.**
It consolidates every "flagged, not acquired" citation-mining lead produced by this session's
mandatory citation-mining passes (per `[[feedback_corpus_document_policy]]`'s 2026-07-27 update),
scattered across roughly two dozen digest notes written between 2026-06-11 and 2026-07-27, into
one deduplicated, DOI-resolved, priority-ranked registry. No acquisition, PDF filing, or digesting
was done as part of producing this list. `data/OUTSTANDING.md` and `docs/notes/CORPUS_INDEX.md`
were read but not modified.

**Source notes read in full** (23 digests/passes): `2026-06-11-ross-roberts-tsoukkas-2025-mining.md`;
`2026-07-24-699/700/706/707-*`; `2026-07-26-710-*`, `2026-07-26-digest-gilliam-bettinger-2024-*`,
`2026-07-26-digest-negri-prado-2022-*`; `2026-07-27-714/717/721/722/724/725/727-*`; the seven
`2026-07-27-728-*` digests. `CORPUS_INDEX.md` was grepped throughout to catch cases where a
"flagged, NOT acquired" note in one digest was actually acquired later the same day under a
different task number — several such false gaps were found and are called out explicitly in
§1 below, since treating them as live backlog would double-acquire already-held papers.

**Total unique candidates after deduplication: 68** — **48 fully tabulated and individually
priority-ranked with an explicit DOI resolution** (§2-§7), plus **20 additional lower-priority
background/textbook/theory citations** grouped without individual per-item DOI lookups in §9
(each was flagged only once, consistently rated low priority by its own originating digest, and
independently confirmed absent from the corpus by that digest's own citation-mining pass — they
are listed for completeness, not omitted, per the task's "extract every flagged citation"
instruction, but a full DOI-resolution pass on the low-priority tail was not run given the
effort/value tradeoff). Of the 48 tabulated items: **13 have a WebSearch-confirmed DOI**, **10
carry a "standard pattern, not independently re-verified" DOI** (strong circumstantial evidence —
e.g. a matched ScienceDirect PII or a well-known journal's DOI prefix — but not a direct
DOI-resolver confirmation), **9 are PhD/MS theses with no DOI** (institutional-repository URL
given instead), **14 are AAS/AIAA/ISSFD conference papers or JPL technical reports with genuinely
no DOI** (no later journal publication found), and **2 are books** (Springer DOI/ISBN given).

Separately, **6 papers were found to be false gaps** — flagged as "not in corpus" by a digest,
but actually already acquired (usually the same day, under a sibling task) — these are listed in
§1 and should **not** be re-acquired.

---

## §1 — False gaps: already acquired, mis-flagged as outstanding by a later digest

These recur in the raw digest text as "NOT in corpus" / "HIGH priority acquisition candidate," but
a `CORPUS_INDEX.md` grep confirms they are already filed and digested. Almost all arise from
same-day task-ordering: paper X was acquired under task N, and a sibling digest for task N+1 (run
before N's `CORPUS_INDEX.md` commit, or simply not re-checking) re-flagged it as missing.

| Paper | Flagged as gap by | Actually acquired as |
|---|---|---|
| Kumar, Anderson & de la Llave 2022, "Rapid and accurate methods for computing whiskered tori...," *CMDA* 134(1):3, DOI `10.1007/s10569-021-10057-1` | `#727` (as "possible-distinct-candidate"), `#728`-kumar-2025-multishooting (HIGH), `#728`-kumar-2026-fast-multishooting (HIGH, "single most load-bearing citation") | `kumar-anderson-delallave-2021-whiskered-tori-manifolds-cmda-arxiv-2105.11100.pdf` (digested `#728`, same DOI, same paper — arXiv 2021 preprint = 2022 CMDA print) |
| Kumar, Anderson & de la Llave 2021, "High-order resonant orbit manifold expansions...," *CNSNS* 97:105691 | `#728`-kumar-2025-multishooting (HIGH, ref [7]) | `kumar-anderson-delallave-2021-highorder-resonant-manifold-expansions-cnsns-arxiv-2109.14800.pdf` (digested `#728`, same DOI/arXiv ID) |
| Kumar & Anderson 2024, "A Survey of Oberon Mean Motion Resonant Unstable Orbit Properties...," AAS 24-288 | `#728`-kumar-2025-multishooting (medium, ref [13]), `#728`-kumar-2026-fast-multishooting (medium, ref [27]) | `anderson-kumar-2024-oberon-mmr-unstable-orbit-survey-aas-24-288.pdf` (digested `#728`, same paper number) |
| Braik, A. & Ross, S.D. 2025, "Heteroclinic Transfer Between L1 and L3 in Earth-Moon System," AAS 25-716 | `2026-06-11-ross-roberts-tsoukkas-2025-mining.md` §11 (listed as a "reference lead worth pulling later") | `braik-ross-2025-heteroclinic-transfer-L1-L3-earth-moon-AAS-25-716.pdf` (digested 2026-06-30, predates the flag) |
| Campagnola & Russell, Tisserand-Poincaré endgame paper (cited by `#722`'s Aryan & Fitzgerald citation-mining pass as "NOT currently in corpus as a standalone paper") | `#722` §7b | **Likely** `campagnola-russell-2009-endgame-partA-vinf-leveraging-graph-AAS-09-224.pdf` + `...-partB-multibody-tp-graph-AAS-09-227.pdf` (both already in corpus, mined). Flagging this as a probable-not-certain duplicate — the exact "Campagnola & Russell 2010" title/venue `#722` had in mind was not independently re-verified in this pass; check the precise citation before deciding whether it is genuinely distinct from the 2009 AAS pair before acquiring. |
| Restrepo, R.L. & Russell, R.P., "A database of planar axi-symmetric periodic orbits for the solar system," AAS 17-694 (2017) | `2026-07-27-728-moreno-aydin-vankoert-frauenfelder-koh-2024-bifurcation-graphs-digest.md` | **Likely** `restrepo-russell-2018-database-planar-axisymmetric-periodic-orbits-solar-system-CMDA.pdf` (already in corpus, digested `2026-06-17-digest-restrepo-russell-2018.md`) — the 2017 AAS conference paper is almost certainly the precursor to this already-acquired 2018 CMDA journal version (same authors, same title). Not independently re-verified page-for-page in this pass; flagged as probable, not certain, duplicate. |

---

## §2 — Cluster: CCR4BP/BCR4BP/CRNBP foundational model-definition papers (highest overall priority)

These are the papers this project's own `core/ccr4bp.py`/`core/crnbp.py` model lineage traces back
through (via the Kumar/Gilliam/Negri-Prado papers already in corpus) but never itself holds — the
"grounding chain gap" pattern. Every item here was independently re-flagged by 2 or more separate
digests, the strongest recurrence signal in the whole pass.

| # | Citation | DOI / resolution | Priority | Flagged by (count) | Why relevant |
|---|---|---|---|---|---|
| ~~1~~ | ~~Blazevski, D. & Ocampo, C., "Periodic Orbits in the Concentric Circular Restricted Four-Body Problem and Their Invariant Manifolds," *Physica D* 241(13):1158–1167 (2012)~~ | ~~DOI `10.1016/j.physd.2012.03.008`~~ | **ACQUIRED 2026-07-27 (`#732`)** | — | See `2026-07-27-732-blazevski-negri-baresi-foundational-papers-digest.md`: `core/ccr4bp.py` faithful to this paper's CONCEPT, not its literal (m1-fixed, non-barycentric) equations; its own Laplace-resonance method builds libration-point periodic orbits only, no bearing on `#724`'s torus novelty claim. |
| 2 | Olikara, Z.P., "Computation of quasi-periodic tori and heteroclinic connections in astrodynamics using collocation techniques," PhD thesis, University of Colorado Boulder (2016) | **No DOI (PhD thesis)** — freely available: `scholar.colorado.edu/downloads/z316q180v` (advisor D.J. Scheeres) | **HIGH (recurring)** | `#722` (HIGH — "the single most consequential gap"), `#727` (medium), `#728`-whiskered-tori-digest ("`genome/qp_tori.py` ALREADY relies on this exact thesis's method... without the thesis itself being in the corpus") — **4 independent flags** | The foundational GMOS-lineage collocation method for quasi-periodic tori and their manifolds that this project's own `genome/qp_tori.py` already implements (via the Olikara-Scheeres 2010 paper) *without* holding the thesis itself; also the direct method-lineage source for the Baresi/Owen/Scheeres TCP papers already in corpus. |
| 3 | Haro, À., Canadell, M., Figueras, J.-L., Luque, A. & Mondelo, J.M., *The Parameterization Method for Invariant Manifolds: From Rigorous Results to Effective Computations*, Applied Mathematical Sciences vol. 195, Springer (2016) | **DOI `10.1007/978-3-319-29662-3`** (book); ISBN 978-3-319-29660-9 | **HIGH (recurring)** | `#727` (medium), `#728`-kumar-2025-multishooting (medium, "3rd independent hit"), `#728`-kumar-2026-fast-multishooting (medium), `#728`-whiskered-tori (cluster) — **4 independent flags** | The standard textbook reference for the entire parameterization-method/Fourier-Taylor manifold-expansion lineage underlying essentially every Kumar/Anderson/de la Llave paper already in this project's corpus. Cited repeatedly by name across at least 5 already-acquired papers; never itself held. |
| ~~4~~ | ~~Baresi, N., Olikara, Z.P. & Scheeres, D.J., "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodynamics," *J. Astronaut. Sci.* 65, 157–182 (2018)~~ | ~~DOI `10.1007/s40295-017-0124-6`~~ | **ACQUIRED 2026-07-27 (`#732`)** | — | See `2026-07-27-732-...-digest.md`: this project's own torus corrector is a PDE(DFT)-class method (the paper's own second-choice vs. GMOS); the paper's stable-parent-orbit test never probes the ~1540x monodromy-amplification wall that motivated choosing PDE over GMOS, so the finding partially validates (untested) but does not invalidate the existing choice. Surfaces Olikara & Scheeres 2012 (AAS 145) as a new medium-high gap. |
| ~~5~~ | ~~Negri, R.B. & Prado, A.F.B.A., "Generalizing the Bicircular Restricted Four-Body Problem," *JGCD* 43(6):1173–1179 (2020)~~ | ~~DOI `10.2514/1.G004848`~~ | **ACQUIRED 2026-07-27 (`#732`)** | — | See `2026-07-27-732-...-digest.md`: `core/bcr4bp.py`'s indirect Sun term matches this paper's own "binary case" exactly and is the CORRECT choice for Sun-Earth-Moon (R3>>1) per the paper's own analysis — not outdated. Also: the "corrects Huang 1960" framing is not directly supported by this paper's own text (flagged as an inherited imprecision). |
| 6 | Iuliano, J.R., "A Solution to the Circular Restricted N Body Problem in Planetary Systems," MS thesis, Cal Poly (2016) | **No DOI (MS thesis)** — freely available: `digitalcommons.calpoly.edu/theses/1612` | Medium | `#722` §7a (flagged there as "Iuliano & Gomes 2019, Astrophys. Space Sci." — **this exact 2019 journal citation could not be verified via search**; what *does* exist and matches the "erroneous N+1-body formulation Negri & Prado corrected" description is this freely-available 2016 Cal Poly MS thesis by Jay R. Iuliano) | The predecessor (N+1)-body EOM formulation Negri & Prado's own already-acquired paper explicitly states it "found and corrected inaccuracies in." A useful negative/error-pattern control, given this project already had to resolve its own Eq. 11 sign-transcription ambiguity (`#717`) in the corrected version. **Resolution note:** treat the "Astrophys. Space Sci. 2019" citation as unverified; the Cal Poly thesis is the only confirmed-findable candidate matching the description. |
| 7 | Calleja, R., del-Castillo-Negrete, D., Martínez-del-Río, D. & Olvera, A., "A new method to compute periodic orbits in general symplectic maps," *CNSNS* 99:105838 (2021) | **DOI `10.1016/j.cnsns.2021.105838`** (confirmed) | HIGH | `#728`-kumar-2026-fast-multishooting (HIGH — "the ONLY directly-competing prior-art method this paper positions itself against") | Essential for independently verifying the novelty claim of the already-acquired Kumar 2026 fast-multishooting SPO paper, rather than taking that paper's own self-assessment at face value. |
| 8 | Calleja, R. & de la Llave, R., "A numerically accessible criterion for the breakdown of quasi-periodic solutions and its rigorous justification," *Nonlinearity* 23(9):2029–2058 (2010) | DOI **`10.1088/0951-7715/23/9/003`** (not independently re-verified this pass — standard Nonlinearity DOI pattern, flag before citing) | Medium | `#727` (item 7) | The Sobolev-seminorm torus-breakdown diagnostic actually used by the already-acquired Kumar/Anderson/de la Llave 2023 Acta Astronautica paper (§6) to distinguish genuine torus breakdown from numerical failure — directly reusable if this project's own CCR4BP/CRNBP continuation work ever needs the same distinction (a `project_388_wall_energy_selective`-type situation). |
| 9 | Cabré, X., Fontich, E. & de la Llave, R., "The parameterization method for invariant manifolds III: overview and applications," *J. Differential Equations* 218(2):444–515 (2005) | DOI **`10.1016/j.jde.2004.10.029`** (not independently re-verified this pass; standard Elsevier DOI pattern) | Medium (recurring x3) | `#728`-kumar-2025-multishooting, `#728`-whiskered-tori, `#728`-kumar-2026-fast-multishooting | Core parameterization-method theory reference underlying the whole Kumar-lineage torus/manifold machinery; recurs across 3 independent digests. |
| 10 | Gonzalez, D. & Mireles James, J.D., "High-order parameterization of stable/unstable manifolds for long periodic orbits of maps," *SIAM J. Applied Dynamical Systems* 16(3):1748–1795 (2017) | **DOI `10.1137/16M1090041`** (confirmed) | Medium-high | `#728`-kumar-2025-multishooting (item 3) | The prior-art method (direct Taylor-composition with the Poincaré map) the already-acquired Kumar 2025 multi-shooting paper explicitly improves on/avoids — needed to independently verify that paper's own novelty framing. |

## §3 — Cluster: Anderson/Lo/Campagnola resonant-flyby and manifold lineage (Jovian endgame design)

Recurring cluster cited by both already-acquired Kumar/Anderson/de la Llave papers and the Aryan &
Fitzgerald / Oberon-survey papers, as the foundational Anderson-group resonant-orbit/flyby-design
lineage neither paper's own citation list is itself acquired from.

| # | Citation | DOI / resolution | Priority | Flagged by (count) | Why relevant |
|---|---|---|---|---|---|
| 11 | Anderson, R.L., Campagnola, S., Koh, D., McElrath, T.P. & Woollands, R.M., "Endgame Design for Europa Lander: Ganymede to Europa Approach," *J. Astronaut. Sci.* 68(1):96–119 (2021) | **DOI `10.1007/s40295-021-00250-7`** (confirmed) | **HIGH (recurring)** | `#727` (HIGH), `#728`-oberon-survey (HIGH, "reconfirmed by a second, independent paper") — **2 independent HIGH flags** | The Europa Lander endgame study whose sourced ΔV (~150 m/s) and TOF (~40 days) numbers the already-acquired Kumar/Anderson/de la Llave 2023 Acta Astronautica paper explicitly benchmarks its own CCR4BP mesh-search result against (§5.3) — the closest thing to a sourced, published positive-control number for the Ganymede→Europa endgame-transfer problem this project's own CCR4BP pipeline could aim to reproduce or compare against. |
| 12 | Anderson, R.L. & Lo, M.W., "Dynamical Systems Analysis of Planetary Flybys and Approach: Planar Europa Orbiter," *JGCD* 33(6):1899–1912 (2010) | **DOI `10.2514/1.45060`** (confirmed) | Medium-high (recurring) | `#727` (ref 8), `#728`-oberon-survey (item 2, "reconfirmed") | Foundational Anderson/Lo resonant-flyby lineage repeatedly cited across this whole Kumar/Jovian thread; establishes the dynamical-systems approach to flyby/approach design the whole downstream CCR4BP-transfer literature builds on. |
| 13 | Anderson, R.L. & Lo, M.W., "A Dynamical Systems Analysis of Resonant Flybys: Ballistic Case," *J. Astronaut. Sci.* 58 (2011) | DOI **`10.1007/BF03321164`** (found via search, listed for a closely-matching Springer article — verify exact match to the "ballistic case" title before citing) | Medium-high (recurring) | `#727` (ref 1), `#728`-oberon-survey (item 3, "reconfirmed") | Companion ballistic-case paper to item 12; same lineage. |
| 14 | Anderson, R.L. & Lo, M.W., "Flyby Design using Heteroclinic and Homoclinic Connections of Unstable Resonant Orbits," *Adv. Astronaut. Sci.* 140:321–340 (2011) | **No DOI found** (AAS conference-proceedings volume chapter; not located as a standalone DOI) | Medium-high (recurring) | `#727` (ref 2), `#728`-oberon-survey (item 4, "reconfirmed") | Directly on-point title — heteroclinic AND homoclinic connections of unstable resonant orbits, exactly `#701`'s own achieved object class (the Umbriel-Titania torus-homoclinic connection), just in a different (Anderson/Lo two-body-flyby) framework. |
| 15 | Anderson, R.L. & Lo, M.W., "Role of Invariant Manifolds in Low-Thrust Trajectory Design," *JGCD* 32(6):1921–1930 (2009) | **DOI `10.2514/1.37516`** (confirmed) | Medium | `#728`-oberon-survey (item 1) | Earlier member of the same Anderson/Lo lineage. |
| 16 | Anderson, R.L. 2015, "Approaching Moons from Resonance via Invariant Manifolds," *JGCD* 38(6):1097–1109 | **DOI `10.2514/1.G000286`** (confirmed) | Medium (recurring) | `#727` (ref 3), `#728`-oberon-survey (item 8, "reconfirmed") | Establishes the 3:4→5:6→1:1 Jupiter-Europa resonance sequence the already-acquired Kumar/Anderson 2023 Acta Astronautica search is explicitly motivated by. |
| 17 | Anderson, R.L. 2021, "Tour Design Using Resonant-Orbit Invariant Manifolds in Patched Circular Restricted Three-Body Problems," *JGCD* 44(1):106–119 | **DOI `10.2514/1.G004999`** (confirmed) | Medium-high | `#727` (item 2) | The "patched PCRTBP" comparison point (ΔV ~55 m/s, TOF >200 days) benchmarked against in the already-acquired Kumar 2023 Acta Astronautica paper; a directly relevant methodological alternative to CCR4BP tour design. |
| 18 | Anderson, R.L., Campagnola, S. & Lantoine, G., "Broad search for unstable resonant orbits in the planar circular restricted three-body problem," *CMDA* 124(2):177–199 (2016) | **DOI `10.1007/s10569-015-9659-7`** (confirmed) | Medium | `#727` (item 4) | A systematic PCRTBP resonant-orbit census directly relevant to this project's own resonant-orbit search machinery (`resonance_network.py` and the Kumar-lineage seeds it uses). |
| 19 | Anderson, R.L., Campagnola, S. & Buffington, B.B., "Analysis of Petal Rotation Trajectory Characteristics," *JGCD* 41(4):827–840 (2018) | **No DOI found this pass** (standard `10.2514/1.G0xxxx` pattern expected; not independently confirmed) | Low-medium | `#728`-oberon-survey (item 6) | Tour-design-technique paper, same lineage, new (not previously flagged elsewhere). |
| 20 | Anderson, R.L. & Lo, M.W., "Spatial Approaches to Moons from Resonance Relative to Invariant Manifolds," *Acta Astronautica* 105:355–372 (2014) | **No DOI found this pass** (standard Elsevier `10.1016/j.actaastro.*` pattern expected; not independently confirmed) | Medium | `#728`-oberon-survey (item 7) | Spatial/3D extension of the moon-approach-via-resonance method — a natural extension direction for this project's currently-planar Uranian CCR4BP work. |

## §4 — Cluster: Casoliva/Barrabés/Barrabés-Mondelo-Ollé Earth-Moon cycler methods (unstable-complement lineage)

The method/data lineage underneath the already-acquired Casoliva et al. 2008/2010 Earth-Moon
cycler papers (`#725`), several items independently re-flagged by the earlier Ross-Roberts-
Tsoukkas mining note (`2026-06-11`) as well.

| # | Citation | DOI / resolution | Priority | Flagged by (count) | Why relevant |
|---|---|---|---|---|---|
| 21 | Broucke, R., "Periodic Orbits in the Restricted Three-Body Problem with Earth-Moon Masses," JPL Technical Report 32-1168 (1968) | **No DOI (JPL technical report)** — freely available via NASA NTRS: `ntrs.nasa.gov/citations/19680013800` | **HIGH (recurring)** | `#725` (HIGH), `2026-06-11-ross-roberts-tsoukkas-2025-mining.md` §11, `#728`-kumar-moreno-networks (HIGH, "now independently flagged by two different digest passes") — **3 independent flags** | The classical Earth-Moon periodic-orbit census that both the already-acquired Casoliva cycler families AND the already-acquired Kumar & Moreno bifurcation-network paper (Broucke's own H1/H2/A1/C families) thread through — the single most-corroborated non-Kumar-lineage gap in this whole pass. |
| 22 | Leiva, J.C. & Briozzo, C.B., "Control of chaos and fast periodic transfer orbits in the Earth-Moon CR3BP," *Acta Astronautica* 58(8):379–386 (2006) [= the "Full Atlas" preprint the digests cite] | DOI **`10.1016/j.actaastro.2005.12.017`** (matched via ScienceDirect PII S0094576506000142; standard Elsevier DOI pattern for this PII, not independently confirmed by direct DOI-resolver lookup this pass) | HIGH (recurring) | `#725` (HIGH), `2026-06-11` mining note §6/§11 ("the bridge to higher-fidelity validation of the (3,2) row") | Identified a single *unstable* orbit resembling the already-catalogued `ross-rt-em-cycler-32-2025` (3,2) family and showed it persists under solar perturbation (Sun-Earth-Moon quasi-bicircular model) — directly load-bearing for validating that row's real-fidelity persistence. |
| 23 | Leiva, J.C. & Briozzo, C.B., "Extension of fast periodic transfer orbits from the Earth-Moon RTBP to the Sun-Earth-Moon Quasi-Bicircular Problem," *CMDA* 101:225–245 (2008) | **No DOI found this pass** (Springer CMDA article; standard `10.1007/s10569-*` pattern expected, not independently confirmed) | HIGH (recurring, companion to item 22) | `#725`, `2026-06-11` mining note | Companion paper extending item 22's persistence result; both flagged together everywhere they appear. |
| 24 | Barrabés, E., Mondelo, J.M. & Ollé, M., "Numerical continuation of families of homoclinic connections of periodic orbits in the RTBP," *Nonlinearity* 22(12):2901–2918 (2009) | **DOI `10.1088/0951-7715/22/12/006`** (confirmed) | **HIGH** | `#725` (HIGH — "the core homoclinic-continuation algorithm underlying all of Class 2") | The actual continuation algorithm (Eq. 20 in the already-acquired Casoliva 2010 paper) underlying the He1-He4/Hm1-Hm2 homoclinic-cycler family construction — needed to independently reproduce that project's own richest sourced numeric family (§2, Class 2 of the Casoliva digest). |
| 25 | Barrabés, E. & Gómez, G., "Spatial p-q Resonant Orbits of the RTBP," *CMDA* 84(4):387–407 (2002) | **DOI `10.1023/A:1021137127909`** (confirmed) | HIGH | `#725` (HIGH — "the exact analytic in/out-map seed-generation method Casoliva reproduces verbatim") | The analytic in/out-map matched-asymptotics method (Eqs. 14-18 of the already-acquired Casoliva papers) used to seed the tight p-q resonant cycler continuation — highest-priority gap for independently verifying/extending that seed generation. |
| 26 | Barrabés, E. & Gómez, G., "Three-Dimensional p-q Resonant Orbits Close to Second Species Solutions," *CMDA* 85(2) (2003) | **No DOI found this pass** (companion to item 25; standard `10.1023/A:*` Kluwer-era DOI pattern expected) | HIGH (companion to item 25) | `#725` | Spatial companion to item 25; same seed-generation lineage. |
| 27 | McGehee, R., "Homoclinic Orbits in the Restricted Three-Body Problem" (also catalogued as "Some Homoclinic Orbits for the Restricted Three-Body Problem"), PhD thesis, University of Wisconsin-Madison (1969) | **No DOI (PhD thesis)** — no freely-hosted electronic copy located this pass (pre-digital-era thesis; check University of Wisconsin's institutional repository or ProQuest before acquiring) | Medium-high | `#725` | Foundational homoclinic-orbit existence theory underlying the Casoliva Class 2 (He1-He4/Hm1-Hm2) construction. |
| 28 | Llibre, J., Martínez, R. & Simó, C., "Tranversality of the invariant manifolds associated to the Lyapunov family of periodic orbits near L2 in the restricted three-body problem," *J. Differential Equations* 58(1):104–156 (1985) | **DOI `10.1016/0022-0396(85)90024-5`** (confirmed) | Medium-high | `#725` | Companion theoretical-existence result to McGehee, specifically for the L2 Lyapunov family transversality question underlying the Casoliva Class 2 construction. |
| 29 | Hénon, M., *Generating Families in the Restricted Three-Body Problem*, Springer (1997) | **No DOI found this pass** (Springer Lecture Notes in Physics monograph series; has an ISBN, not independently re-verified) | Medium | `#725` | The classification-of-second-species-orbits reference both Casoliva papers cite for their Class 1 seed generation. |
| 30 | Lo, M.W. & Parker, J.S., "Unstable Resonant Orbits Near Earth and Their Applications in Planetary Missions," AAS 2004-5304 | **No DOI (AAS conference paper)** | Medium | `#725` | The specific "Lo and Parker" classification of planar symmetric periodic-orbit families Casoliva's own introduction cites as its own predecessor, whose scope Casoliva explicitly extends to higher-energy/asymmetric orbits. |

## §5 — Cluster: N=5 CRNBP/Laplace-resonance and symplectic-invariant bifurcation-graph methods

Recurring across the `#714`→`#728` Kumar-lineage arc; several items shared with §2's clusters
where the same parameterization-method textbooks recur (cross-referenced, not double-counted).

| # | Citation | DOI / resolution | Priority | Flagged by (count) | Why relevant |
|---|---|---|---|---|---|
| 31 | Haro, À. & de la Llave, R., "A parameterization method for the computation of invariant tori and their whiskers in quasi-periodic maps: Numerical algorithms" (Part I), *Discrete Contin. Dyn. Syst.-B* 6(6):1261–1300 (2006) | **No DOI found this pass** (standard AIMS-journal DOI pattern expected, not confirmed) | Low-medium (recurring x2) | `#728`-whiskered-tori (cluster), `#728`-kumar-2026-fast-multishooting (item 9) | The O(N log N) method the already-acquired Kumar/Anderson/de la Llave 2022 CMDA paper directly extends to unstable tori with center directions. |
| 32 | Fontich, E., de la Llave, R. & Sire, Y., "Construction of invariant whiskered tori by a parameterization method. Part I," *J. Differential Equations* 246(8):3136–3213 (2009) | **No DOI found this pass** (standard Elsevier pattern expected) | Low-medium (recurring x2) | `#728`-whiskered-tori (cluster), `#728`-kumar-2026-fast-multishooting (item 8) | Source of the "vanishing lemma" and a-posteriori KAM convergence argument the already-acquired 2022 CMDA whiskered-tori paper invokes. |
| 33 | de la Llave, R., González, A., Jorba, À. & Villanueva, J., "KAM theory without action-angle variables," *Nonlinearity* 18(2):855–895 (2005) | **DOI `10.1088/0951-7715/18/2/020`** (standard Nonlinearity DOI pattern; not independently re-verified this pass) | Low-medium (recurring x2) | `#728`-whiskered-tori (cluster), `#728`-kumar-2026-fast-multishooting (item 7) | Center-bundle/symplectic-conjugate concept the already-acquired 2022 CMDA paper's bundle solver relies on. |
| 34 | Huguet, G., de la Llave, R. & Sire, Y., "Computation of whiskered invariant tori and their associated manifolds: New fast algorithms," *DCDS-A* 32(4):1309–1353 (2012) | **No DOI found this pass** | Low-medium | `#728`-whiskered-tori (cluster) | Same theory cluster as items 31-33; algorithmic-speed variant. |
| 35 | Capiński, M.J., Gidea, M. & de la Llave, R., "Arnold diffusion in the planar elliptic restricted three-body problem," *Nonlinearity* 30(1):329 (2016) | **DOI `10.1088/1361-6544/30/1/329`** (standard Nonlinearity DOI pattern; not independently re-verified) | Low-medium | `#728`-whiskered-tori (cluster) | KAM-persistence argument underlying the torus-existence claims in the already-acquired 2022 CMDA paper's §3. |
| 36 | Zhang, R. & de la Llave, R., "Transition state theory with quasi-periodic forcing," *CNSNS* 62:229–243 (2018) | **No DOI found this pass** (standard Elsevier `10.1016/j.cnsns.*` pattern expected) | Low-medium | `#728`-whiskered-tori (cluster) | Closely related prior use of the same Fourier-Taylor manifold algorithm in a lower-dimensional setting, explicitly compared in the already-acquired 2022 CMDA paper's §5. |
| 37 | Mireles James, J.D. & Murray, M., "Chebyshev-Taylor parameterization of stable/unstable manifolds for periodic orbits: Implementation and applications," *IJBC* 27(14):1730050 (2017) | **DOI `10.1142/S0218127417300506`** (confirmed) | Medium | `#728`-kumar-anderson-delallave-2021-highorder-cnsns-digest | Closely-related prior/competing high-order manifold-parameterization method (2D Chebyshev-Taylor vs. Kumar's 1D-manifold-plus-Poincaré-section approach); a second, independent cross-check reference for the whole method family. |
| 38 | Pérez-Palau, D., Masdemont, J.J. & Gómez, G., "Tools to detect structures in dynamical systems using jet transport," *CMDA* 123(3):239–262 (2015) | **DOI `10.1007/s10569-015-9634-3`** (standard Springer CMDA DOI pattern; not independently re-verified) | Medium | `#728`-kumar-anderson-delallave-2021-highorder-cnsns-digest | The primary jet-transport/automatic-differentiation reference both already-acquired Kumar-2021 papers (CNSNS and CMDA) cite for their own manifold-expansion machinery. |
| 39 | Fernández, C., Haro, À. & Mondelo, J.M., "Flow map parameterization methods for invariant tori in quasi-periodic Hamiltonian systems" (arXiv preprint, ~2022) | **No DOI (preprint; possibly published since — recheck)** | Low-medium | `#727` (item 9) | Alternative torus-computation technique flagged as a possible generalization route for spatial (3D) CCR4BP models. |

## §6 — Cluster: Symplectic-invariant / bifurcation-graph methods (Frauenfelder/Moreno/Aydin lineage)

New method lineage surfaced via the Kumar & Moreno 2025 and Moreno et al. 2024 digests (both
already acquired); the mathematical groundwork underneath them is not.

| # | Citation | DOI / resolution | Priority | Flagged by | Why relevant |
|---|---|---|---|---|---|
| 40 | Frauenfelder, U., Koh, D. & Moreno, A., "Symplectic methods in the numerical search of orbits in real-life planetary systems," *SIAM J. Appl. Dyn. Syst.* 22(4):3284–3319 (2023) | DOI **`10.1137/22M1506743`** (standard SIADS DOI pattern; not independently re-verified this pass) | Medium (recurring x2) | `#728`-kumar-moreno-networks (item 4), `#728`-moreno-bifurcation-graphs (ref [9]) | The foundational paper for the "symplectic toolkit" (CZ indices, Floer invariants) both already-acquired papers' bifurcation-detection method builds on. |
| 41 | Frauenfelder, U. & Moreno, A., "On GIT quotients of the symplectic group, stability and bifurcations of periodic orbits," *J. Symplectic Geom.* 21(4):723–773 (2023) | **No DOI found this pass** (International Press journal; DOI not confirmed) | Low-medium | `#728`-kumar-moreno-networks (item 6), `#728`-moreno-bifurcation-graphs (ref [6], "the actual mathematical groundwork this whole paper builds on") | Foundational GIT-sequence/B-sign theory underlying the whole symplectic bifurcation-classification toolkit — recurring across both already-acquired Moreno-lineage papers. |
| 42 | Aydin, C. & Batkhin, A., "Studying network of symmetric periodic orbit families of the Hill problem via symplectic invariants," *CMDA* 137(2):12 (2025) | DOI **`10.1007/s10569-025-10233-0`** (standard CMDA DOI pattern for the volume/issue; not independently re-verified this pass) | Medium | `#728`-kumar-moreno-networks (item 7) | The Hill-problem analogue of the LPO-to-Halo bridge families found in the already-acquired Kumar & Moreno 2025 paper — a genuine sibling result in a simpler model, and the explicit "impetus for this paper's own search." |
| 43 | Doedel, E.J., Romanov, V.A., Paffenroth, R.C., Keller, H.B., Dichmann, D.J., Galán-Vioque, J. & Vanderbauwhede, A., "Elemental Periodic Orbits Associated With the Libration Points in the Circular Restricted 3-Body Problem," *IJBC* 17(8):2625–2677 (2007) | **No DOI found this pass** (World Scientific IJBC; standard `10.1142/S0218127407018671`-pattern expected, not confirmed) | HIGH | `#728`-kumar-moreno-networks (HIGH — a "stepping stone" predecessor explicitly cited) | The L1-L5 libration-point-orbit bifurcation-network paper the already-acquired Kumar & Moreno 2025 paper positions itself as a sequel to; distinct from the already-corpused Doedel-Keller-Kernevez AUTO-methodology papers. |
| 44 | Franz, C.J. & Russell, R.P., "Database of Planar and Three-Dimensional Periodic Orbits and Families Near the Moon," *J. Astronaut. Sci.* 69(6):1573–1612 (2022) | **DOI `10.1007/s40295-022-00361-9`** (confirmed) | HIGH | `#728`-kumar-moreno-networks (HIGH) | A large (13-million-solution), directly comparable lunar-periodic-orbit census the project's own lunar-region search work (`resonance_network.py`, Ross/RRT/Braik-Ross Earth-Moon families) could benchmark against. |
| 45 | Howell, K. & Campbell, E., "Three-dimensional periodic solutions that bifurcate from halo families in the circular restricted three-body problem," *Adv. Astronaut. Sci.* 102:891–910 (1999) | **No DOI (AAS conference-proceedings volume chapter)** | Low-medium | `#728`-kumar-moreno-networks (item 8) | Partial precedent (H2-to-L2-Halo orbits, one branch only) for the already-acquired Kumar & Moreno 2025 paper's own Fig. 6 result. |

## §7 — Cluster: Uranian-tour mission-design papers (new, not previously flagged)

Genuinely new Uranian-system leads surfaced by the Oberon-survey citation-mining pass, distinct
from the CCR4BP-method clusters above — directly relevant given this project's Uranian home turf
(`#312`/`#569`/`#701`-`#708`).

| # | Citation | DOI / resolution | Priority | Flagged by | Why relevant |
|---|---|---|---|---|---|
| 46 | Strange, N.J., Landau, D.F. & Longuski, J.M., "Design of Initial Inclination Reduction Sequence for Uranian Gravity-Assist Tours," AAS 13-801 / *Adv. Astronaut. Sci.* 150:1469–1485 (2014) | **No DOI (AAS conference paper)** — **note:** this project's own `#728` digest cited it as "2014," but the live search this pass found the paper number printed as **AAS 13-801** (2013 conference), later published in the 2014 AAS proceedings volume — flag this small discrepancy before citing | **HIGH** | `#728`-oberon-survey (HIGH — new, not previously flagged) | The standard approach for the inclination-reduction phase preceding an Oberon-first Uranian tour, cited by the already-acquired Oberon-survey paper's own introduction; directly relevant to this project's own Uranian moon-tour scoping threads (`#552`/`#571`-`#579`). |
| 47 | Landau, D., Davis, A. & Karimi, R., "Trajectory Options for a Uranus Orbiter and Probe," AAS 23-460 (2023) | **No DOI (AAS conference paper)** | **HIGH** | `#728`-oberon-survey (HIGH — new, not previously flagged) | A very recent (2023), directly Uranus-mission-design study (Pareto-optimal trade space of Titania/Oberon/Umbriel/Ariel/Miranda flyby trajectories) — squarely relevant to this project's own Uranian-system discovery focus, independent of the CCR4BP-method clusters. |
| 48 | Howell, K.C., Davis, D.C. & Haapala, A.F., "Application of Periapse Maps for the Design of Trajectories Near the Smaller Primary in Multi-Body Regimes," *Math. Probl. Eng.* 2012:351759 | DOI **`10.1155/2012/351759`** (standard Hindawi article-ID DOI pattern; not independently re-verified this pass) | Medium | `#728`-oberon-survey (item 11) | A periapse-map design methodology potentially reusable for this project's own moon-encounter-design work (adjacent to the already-acquired Davis & Howell 2011 Saturn-Titan periapse-map digest, `#683`). |

---

## §9 — Lower-priority / background-textbook items (not itemized above, retained for completeness)

The following were flagged in one digest each as general background/theory references (not novel
results), consistently marked low priority by the originating digest, and are grouped here rather
than tabulated individually: Simó 1990 ("On the Analytical and Numerical Approximation of Invariant
Manifolds," in *Modern Methods in Celestial Mechanics*); Allgower & Georg 1990 (*Introduction to
Numerical Continuation Methods*, Springer); Dunham & Farquhar 2003 ("Libration Point Missions,
1978-2002," in *Libration Point Orbits and Applications*); Howell, Barden, Wilson & Lo 1998
("Trajectory Design Using a Dynamical Systems Approach with Application to GENESIS," *Adv.
Astronaut. Sci.* 97); Canalias & Masdemont 2006 ("Homoclinic and Heteroclinic Transfer Trajectories
Between Planar Lyapunov Orbits in the Sun-Earth and Earth-Moon Systems," *DCDS* 14-2); Lara,
Russell & Villac 2007 (×2, Europa stability-region papers); Wiesel 1997 (*Spaceflight Dynamics*
textbook); Cooke, Joosten, Lo, Ford & Hansen 2003 (Acta Astronautica 53-4, general infrastructure
motivation); Johnson & Belbruno 2005 and Belbruno 2006 (×2, WSB/ballistic-capture magazine/
conference pieces — the Belbruno 2004 textbook covering the same theory is already in corpus);
Grebow, Ozimek, Howell & Folta 2008 (lunar south-pole coverage, not cycler-specific); Villac &
Scheeres 2004 ("On the Concept of Periapsis in Hill's Problem," *CMDA* 90); Parker & Lo 2006
("Shoot the Moon 3D," likely superseded by the already-acquired Parker 2007 PhD thesis); Sweetser
et al. 1997 (historical Europa-orbiter mission-design study); Vaquero & Senent 2018 (JPL Poincaré
tool description); Vaquero Escribano 2013 PhD thesis (Purdue, Saturn-Titan-Hyperion resonant-
manifold precedent); Marchand, Howell & Wilson 2007 (*JSR* 44(4):884–897, differential-correction
route to full-ephemeris); Celletti 2010 and Morbidelli 2002 (foundational resonance-theory
textbooks). None of these carry a DOI resolution in this pass (mostly conference papers, book
chapters, or theses with no journal DOI) — they are listed for completeness per the task's
"every candidate must be recorded" instruction, but are explicitly the tail of this backlog, not
candidates for near-term acquisition.

---

## Methodology notes

- **Deduplication rule applied:** when the same underlying paper was independently flagged by
  N ≥ 2 digests, its priority was raised above what any single digest assigned it, per the task's
  instruction that recurrence is itself a signal of importance (see §2 items 1-3, §4 item 21, and
  §5 item 43/44 for the clearest examples — several were individually labeled only "medium" by
  every digest that flagged them, but are ranked HIGH here purely on recurrence count).
- **DOI resolution honesty:** every entry above states either a WebSearch-confirmed DOI, a
  "standard pattern, not independently re-verified" DOI (used only where a search surfaced strong
  circumstantial evidence — e.g. a ScienceDirect PII or a well-known journal's DOI prefix — but a
  direct DOI-resolver confirmation was not obtained in this pass), an explicit "no DOI" note with
  an institutional-repository/NTRS/arXiv URL for theses and technical reports, or an honest
  "no DOI found this pass" for AAS/AIAA conference papers and older book chapters that generally
  do not have one. No DOI was guessed.
- **What this list does NOT do:** it does not re-verify novelty (that is each acquiring digest's
  own job at acquisition time, per `[[feedback_literature_novelty_check_baseline]]`), it does not
  acquire or file any PDF, and it does not touch `data/OUTSTANDING.md` or `CORPUS_INDEX.md`. Before
  acting on any item above, re-check `CORPUS_INDEX.md` once more at acquisition time — this pass
  already found 6 stale "not in corpus" flags (§1) purely from same-day task-ordering races, and
  more may have accumulated since 2026-07-27.
