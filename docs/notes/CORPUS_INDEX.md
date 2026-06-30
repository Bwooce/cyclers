# Corpus Master Index (CORPUS_INDEX.md)

**LIVING DOCUMENT.** One row per file in the private corpus
`cyclers_pdf/papers/`. This is the discoverability layer **and** the anti-slip
ledger mandated by `docs/notes/corpus-document-policy.md` §3 — its absence is
precisely what let Szebehely 1967 (the foundational CR3BP textbook the genome
rests on) sit undigested for weeks.

**Maintenance rule (policy §3):** this index is updated **in the same commit**
as any future filing or digest. Filing a PDF is not "done" until its digest
note (or mined-by pointer) **and** its index line both exist — *"Filed,
digested, indexed."* Re-run the coverage sweep (task #397 method, below) after
every large acquisition wave.

**Built:** 2026-06-19 (task #397), sweeping all 130 files then in `papers/`.
Supersedes the one-shot `2026-06-13-paper-corpus-digest-audit.md` (which
covered ~67 files and is now a historical snapshot).

## Status vocabulary

Per file: a **digest/mined status** and an **OCR status**.

- `digested` — has a dedicated digest note `docs/notes/*digest*.md` (or a
  per-paper section in a multi-paper digest), full verdict + sourced citations.
- `mined` — referenced by a `docs/notes/*mining*.md` method/data note. Its
  tables/methods power real work even without a standalone digest note.
- `mined-by-catalogue` — cited by a `data/catalogue.yaml` row
  (first_published / corroborating / orbit_source, DOI- or table-traced).
- `mined-by-KNOWN_CORPUS` — a `literature_check.py` `KNOWN_CORPUS` anchor
  (DOI/author/citation), used for the novelty gate.
- `triaged` — read and given an explicit in-/out-of-scope verdict in a triage
  note (background ML/GNC/EDL sweep). A verdict counts as processed.
- `undigested-unmined` — none of the above. The Szebehely-class gap.

(A file often carries several at once; the row lists the strongest pointer
plus any others.)

OCR status (from `verify/ocr.has_text_layer`, 10.0 chars/page threshold,
measured 2026-06-19):

- `text-layer` — native/searchable text. Just `pdftotext`.
- `image-only` — scanned page-images; needs `ocrmypdf` before any digest can
  `pdftotext` it. Vision-read precision pages per policy §1 hybrid rule.
- `text-layer (thin)` — passes the 10 char/page floor but carries only a
  sparse/partial OCR layer (≈20–250 chars/page); usable for navigation, but a
  precision read should treat it like image-only and vision-read the page.

## Genuine gap list — see §"Gaps" at the bottom

**ZERO genuine `undigested-unmined` files remain** (2026-06-19). The last
three — shepperd-1985 (#402), bond-allman-2021 (#396), willis-2008 — were
digested 2026-06-19; Szebehely was digested earlier (#394/#400). All other
prior gaps (Fan 2025, Li-Topputo 2019, Silvestrini 2022, Singh-Junkins 2022)
were closed by the 2026-06-13 background-papers triage. Every `papers/` file
is now digested, mined, or triaged.

---

## Foundational theory / textbooks / methods (the highest-slip-risk class)

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| szebehely-1967-theory-of-orbits-restricted-problem-three-bodies-book.pdf | 2026-06-19-digest-szebehely-1967.md; -mining; -goldens | The foundational CR3BP textbook (Jacobi integral, periodic-orbit families, stability) — the genome's implicit base | digested + mined | image-only → OCR'd .txt committed (2026-06-22; 217676 words; #400 cache never persisted) |
| gurfil-ed-2007-modern-astrodynamics-elsevier-astrodynamics-series-vol-1-textbook.pdf | 2026-06-19-digest-gurfil-2007-modern-astrodynamics.md | Methods reference (Elsevier Astrodynamics vol. 1); perturbation / optimization methods | digested | text-layer |
| hintz-2023-orbital-mechanics-astrodynamics-techniques-tools-space-missions-springer-textbook.pdf | 2026-06-17-digest-hintz-2023.md | Springer astrodynamics textbook; techniques/tools reference | digested | text-layer |
| belbruno-2004-capture-dynamics-chaotic-motions-...-low-energy-transfers-princetonUP-textbook.pdf | 2026-06-17-digest-belbruno-2004.md | Low-energy/ballistic-capture dynamics textbook (WSB transfers) | digested | text-layer |
| parker-2007-low-energy-ballistic-lunar-transfers-phd-thesis-cu-boulder.pdf | 2026-06-17-digest-parker-2007.md | PhD thesis; low-energy ballistic lunar transfer design | digested | text-layer |
| vallado-1991-methods-astrodynamics-computer-approach-USAFA-TR-91-6.pdf | 2026-06-10-vallado-1991-tr916-mining.md | Methods of astrodynamics, computer approach (USAFA TR) | mined | text-layer |
| russell-2004-dissertation.pdf | 2026-06-07-russell-2004-dissertation-method-mining.md (+continuation/member-tables); catalogue orbit_source | Periodic-orbit continuation dissertation; major catalogue orbit source | mined + mined-by-catalogue | text-layer |
| mcconaghy-2004-design-optimization-interplanetary-spacecraft-trajectories-purdue-phd.pdf | 2026-06-17-digest-mcconaghy-2004.md; 2026-06-10-mcconaghy-2004-dissertation-mining.md | Interplanetary trajectory design/optimization PhD; cycler tables | digested + mined | text-layer |
| ceriotti-2010-global-optimisation-multiple-gravity-assist-glasgow-phd.pdf | 2026-06-07-ceriotti-2010-mga-global-opt-mining.md | Global optimisation of MGA trajectories PhD | mined | text-layer |
| agrawal-2022-orbital-logistics-architecture-sustainable-mars-exploration-purdue-phd.pdf | 2026-06-04-agrawal-landau-howe-mining.md | Orbital logistics architecture for sustainable Mars (PhD) | mined | text-layer |
| doedel-keller-kernevez-1991-numerical-analysis-bifurcation-problems-I-...IJBC-1(3).pdf | 2026-06-17-digest-doedel-keller-kernevez-1991.md | Bifurcation numerics part I (finite dimensions) — AUTO foundations | digested | text-layer (thin) |
| doedel-keller-kernevez-1991-numerical-analysis-bifurcation-problems-II-...IJBC-1(4).pdf | 2026-06-17-digest-doedel-keller-kernevez-1991.md | Bifurcation numerics part II (infinite dimensions) | digested | text-layer (thin) |
| doedel-paffenroth-keller-2003-computation-periodic-solutions-conservative-systems-3-body-IJBC.pdf | 2026-06-17-digest-doedel-2003.md | Periodic-solution computation in conservative 3-body systems | digested | text-layer |
| bond-allman-2021-modern-astrodynamics-fundamentals-perturbation-methods-princetonUP-...textbook.pdf | 2026-06-19-digest-bond-allman-2021-modern-astrodynamics.md | Bond/Allman two-body + perturbation textbook; Sundman/Sperling-Burdet regularization (Ch 9, feeds tulip Sundman), universal-variable f/g + Lambert (Ch 5-6), VOP element rates (Ch 8), metre-level L4/L5 perturbed-two-body case (Table 9.3) usable as CR3BP golden | digested | text-layer |
| willis-2008-book-review-modern-astrodynamics-gurfil-ed-asr-doi-10.1016-j.asr.2007.07.047.pdf | 2026-06-19-digest-willis-2008-book-review.md | 2-page book review of the Gurfil-ed 2007 volume; corroborating metadata only (defer to the #395 Gurfil digest), zero catalogue impact | digested | text-layer |
| canalias-2007-thesis-mission-design.pdf | 2026-06-20-digest-canalias-2007-se-em-manifolds.md; canalias_se_em_connection.yaml | SE-EM manifold connections, heteroclinic bifurcations PhD | digested + mined | text-layer |
| koon-lo-marsden-ross-2000-shoot-the-moon.pdf | 2026-06-21-digest-koon-2000-shoot-the-moon.md | Shoot the Moon - Hiten trajectory via invariant manifolds | digested | text-layer |
| koon-lo-marsden-ross-2006-dynamical-systems-mission-design.pdf | 2026-06-21-digest-klmr-2006-book.md | KLMR Dynamical Systems and Space Mission Design book | digested | text-layer |

## Earth-Mars cyclers (core domain)

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| friedlander-niehoff-byrnes-longuski-1986-circulating-transportation-orbits-earth-mars.pdf | 2026-06-22-digest-friedlander-niehoff-byrnes-longuski-1986.md; 2026-06-22-dv-band-definitions.md | Canonical early cycler paper: VISIT "free orbit" (near-ballistic, V∞ 4.2–4.5 km/s) + Up/Down Escalator (230 m/s/orbit → ~1.6 km/s/15yr). The last #116/#416 acquisition gap (filed 2026-06-22) | digested | text-layer |
| byrnes-longuski-aldrin-1993-cycler-orbit-earth-mars-jsr-doi-10.2514-3.25519.pdf | 2026-06-17-digest-byrnes-longuski-aldrin-1993.md; KNOWN_CORPUS; catalogue | The Aldrin cycler paper | digested + mined-by-catalogue + KNOWN_CORPUS | text-layer |
| rauwolf-friedlander-nock-2002-mars-cycler-low-thrust-AIAA-2002-5046.pdf | 2026-06-22-digest-rauwolf-friedlander-nock-2002.md | SEP Aldrin cycler architecture; maintenance ΔV 1.56/1.61 km/s/15yr (3 of 7 orbits) + SEP exec; #415/#416 ΔV-band powered-class anchor | digested | text-layer |
| mcconaghy-longuski-byrnes-2002-analysis-broad-class-earth-mars-cycler-...AIAA-2002-4420.pdf | 2026-06-17-digest-mcconaghy-2002.md; catalogue | Broad-class E-M cycler analysis (conference) | digested + mined-by-catalogue | text-layer |
| byrnes-mcconaghy-longuski-2002-various-two-synodic-period-earth-mars-cycler-trajectories-AIAA-AAS-monterey.pdf | 2026-06-19-digest-byrnes-mcconaghy-longuski-2002-two-synodic-cyclers.md | Byrnes-led sibling to AIAA 2002-4420 (#384); three 2-synodic E-M cycler Cases (1/2/3) circular-coplanar + Case 3 real-eph V_inf table; Case 3 ≈ S1L1-B near-twin; precedence/envelope source for the s1l1 V_inf | digested | text-layer |
| mcconaghy-longuski-byrnes-2004-analysis-class-earth-mars-cycler-...jsr-doi-10.2514-1.11939.pdf | 2026-06-17-digest-mcconaghy-2004.md; catalogue | Broad-class E-M cycler analysis (journal) | digested + mined-by-catalogue | text-layer |
| mcconaghy-russell-longuski-2005-standard-nomenclature-earth-mars-cycler-...jsr-doi-10.2514-1.8123.pdf | 2026-06-17-digest-mcconaghy-2005.md; catalogue | Standard cycler nomenclature (S/L convention) | digested + mined-by-catalogue | text-layer |
| mcconaghy-landau-yam-2006-notable-two-synodic-period-earth-mars-cycler-...jsr-doi-10.2514-1.15215.pdf | 2026-06-17-digest-mcconaghy-2006.md; catalogue | Two-synodic-period (S1L1) cycler | digested + mined-by-catalogue | text-layer |
| rogers-2012-vinf-leveraging-cyclers-AIAA-2012-4746.pdf | catalogue orbit_source rogers-2012-t1 (source_quotes) | V∞-leveraging cyclers; major establishment-DV source | mined-by-catalogue | text-layer |
| rogers-hughes-longuski-2015-establishing-cycler-trajectories-earth-mars-...j.actaastro.2015.03.002.pdf | 2026-06-17-digest-rogers-2015.md; catalogue | Establishing E-M cycler trajectories | digested + mined-by-catalogue | text-layer |
| russell-ocampo-2003-systematic-method-earth-mars-cyclers-direct-return-...AAS-03-145.pdf | 2026-06-17-digest-russell-ocampo-2003.md; KNOWN_CORPUS | Systematic E-M cycler + direct-return method | digested + KNOWN_CORPUS | text-layer |
| luidens-1964-mars-nonstop-round-trip-trajectories-AIAA-journal-2-2.pdf | 2026-06-17-digest-luidens-1964.md | Mars non-stop round-trip free-return trajectories | digested | text-layer |
| gravier-marchal-culp-1972-optimal-trajectories-earth-mars-true-planetary-orbits-jota.pdf | 2026-06-17-digest-gravier-1972.md | Optimal E-M trajectories in true planetary orbits | digested | text-layer |
| pontani-conway-2018-optimal-trajectories-hyperbolic-rendezvous-earth-mars-cycling-...jgcd.pdf | 2026-06-17-digest-pontani-2018.md | Hyperbolic rendezvous with cycling spacecraft | digested | text-layer |
| jesick-2019-mars-trojan-orbits-continuous-earth-mars-communication-jas.pdf | 2026-06-17-digest-jesick-2019.md | Mars-Trojan orbits for continuous E-M comms | digested | text-layer |
| patel-2019-earth-mars-cycler-vehicle-conceptual-design-FIT-etd.pdf | s1l1-target-topology-mining.md; 2026-06-05-v42-backfill-sweep.md §7 | E-M cycler vehicle conceptual design (thesis) | mined | text-layer |
| adamo-2025-spanning-earth-mars-chasm-synodic-resonant-waypoints-AIAA-houston-LnL.pdf | 2026-06-17-digest-adamo-2025.md | Synodic-resonant waypoints spanning the E-M chasm | digested | text-layer |
| howe-2025-tackling-mars-cycler-design-head-on-ICES-2025-555.pdf | 2026-06-17-digest-howe-2025.md; 2026-06-04-agrawal-landau-howe-mining.md | Tackling Mars cycler design head-on | digested + mined | text-layer |
| howe-blincow-hall-2025-tackling-mars-cycler-design-head-on-ICES-2025-555.pdf | (duplicate of howe-2025 above, full author list) | Same paper, full-author filename | digested + mined | text-layer |
| landau-longuski-2006-human-mars-trajectories-pt1-impulsive-JSR.pdf | 2026-06-04-agrawal-landau-howe-mining.md; catalogue | Human Mars trajectories pt1 (impulsive) | mined + mined-by-catalogue | text-layer |
| landau-longuski-2009-comparative-assessment-human-mars-technologies-architectures.pdf | 2026-06-04-agrawal-landau-howe-mining.md | Comparative assessment of human Mars architectures | mined | text-layer |
| nock-friedlander-1987-elements-mars-transportation-system-...0094-5765(87)90189-5.pdf | 2026-06-17-digest-nock-friedlander-1987.md | Elements of a Mars transportation system | digested | text-layer |
| vasile-summerer-depascale-2005-earth-mars-evolutionary-branching-ActaAstro.pdf | 2026-06-05-v42-backfill-sweep.md; 2026-06-07-external-algorithms-survey.md | E-M evolutionary-branching global search | mined | text-layer |
| chen-2002-earth-mars-cyclers.pdf | 2026-06-20-digest-new-papers.md | Earth-Mars cyclers (Master's thesis) | digested | text-layer |
| chen-2012.pdf | 2026-06-20-digest-new-papers.md | K.J. Chen 2012 work on Earth-Mars trajectories | digested | text-layer |
| chen-mcconaghy-okutsu-longuski-2002-low-thrust-aldrin-cycler-AIAA-2002-4421.pdf | 2026-06-20-digest-new-papers.md | "A Low-Thrust Version of the Aldrin Cycler" (AIAA 2002-4421). NOTE: was misnamed `chen-russell-ocampo-2002-...-free-return` (wrong authors + not free-return); renamed 2026-06-22 | digested | text-layer |
| chen-mcconaghy-landau-longuski-aldrin-2005-powered-em-cycler-three-synodic-JSR-42-5.pdf | 2026-06-20-digest-new-papers.md | "Powered Earth-Mars Cycler with Three-Synodic-Period Repeat Time" (JSR 42(5)). NOTE: was misnamed `chen-russell-ocampo-2005-...-multiple-impulse-free-return`; renamed 2026-06-22 | digested | text-layer |
| machado-2020.pdf | 2026-06-20-digest-new-papers.md | Parametric design of E-M cycler crew transfer vehicle | digested | text-layer |
| russell-ocampo-2006-optimization-broad-class-ephemeris-model-earth-mars-cyclers-JGCD-29.pdf | 2026-06-20-digest-new-papers.md; 2026-06-22-dv-band-definitions.md | Optimization of a Broad Class of Ephemeris-Model E-M Cyclers (JGCD 29(2), DOI 10.2514/1.13652) — the peer-reviewed source for the <1/<10/<300 m/s-per-7-cycle tiers (9/39/74 parent cyclers). NOTE: previously misfiled as `russell-2006-systematic-method-design-earth-mars-cyclers.pdf`; the actual "systematic method" paper is the 2003 AAS-03-145 (separate row above). | digested | text-layer |

## Mars free-return / human-flyby missions

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| okutsu-longuski-2002-mars-free-returns-via-gravity-assist-venus-JSR-39-1.pdf | 2026-06-07-okutsu-tito-free-returns-mining.md | Mars free returns via Venus GA | mined | image-only → OCR'd .txt (2026-06-22) |
| hughes-edelman-longuski-2014-fast-mars-free-returns-venus-ga-AIAA-2014-4109.pdf | 2026-06-07-hughes-2014-fast-mars-free-returns-mining.md | Fast Mars free returns via Venus GA | mined | text-layer |
| tito-maccallum-carrico-2013-feasibility-manned-mars-free-return-2018-IEEE-aerospace.pdf | 2026-06-07-okutsu-tito-free-returns-mining.md; 2026-06-13-tito-maccallum-2018-free-return-reproduction.md | Inspiration Mars 2018 manned free-return feasibility | mined | text-layer |
| donahue-duggan-2022-boeing-mars-2033-human-flyby-IAC-22-B3-8-x70674.pdf | 2026-06-07-donahue-duggan-2022-mars2033-flyby-mining.md | Boeing Mars-2033 human flyby | mined | text-layer |
| conte-spencer-2018-mission-analysis-earth-to-mars-phobos-dro-...j.actaastro.2018.06.049.pdf | 2026-06-17-digest-conte-spencer-2018.md | E-to-Mars/Phobos DRO mission analysis | digested | text-layer |
| kakoi-howell-folta-2014-access-mars-from-earth-moon-libration-orbits-...j.actaastro.2014.06.010.pdf | 2026-06-17-digest-kakoi-2014.md | Accessing Mars from E-M libration orbits | digested | text-layer |
| putnam-braun-2005-entry-system-options-human-return-moon-mars-AIAA-2005-5915.pdf | 2026-06-13-background-papers-read-triage.md #11 — OUT-OF-SCOPE (EDL) | Entry-system options for human Moon/Mars return | triaged | text-layer |

## Earth-Moon / cislunar cyclers & dynamics

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| genova-aldrin-2015-earth-moon-cycler-AAS-15.pdf | 2026-06-10-genova-aldrin-2015-mining.md | Genova & Aldrin, "A Free-Return Earth-Moon Cycler Orbit" (AAS 2015); real-eph maintenance 20–62 m/s/cycle (~39 m/s/month). NOTE: byte-identical dup `...-free-return-...pdf` removed 2026-06-22 | mined | text-layer |
| ross-roberts-tsoukkas-2025-stable-ballistic-earth-moon-cyclers-AAS-25-621.pdf | 2026-06-11-ross-roberts-tsoukkas-2025-mining.md; 2026-06-12-ross-adoption-results.md; catalogue | Stable ballistic E-M cyclers (5 families) | mined + mined-by-catalogue | text-layer |
| roberts-tsoukkas-ross-2026-stable-prograde-em-cyclers-journal.pdf | 2026-06-13-roberts-tsoukkas-2026-multi-orbiter-journal-mining.md; KNOWN_CORPUS | Stable prograde E-M cyclers (journal) | mined + KNOWN_CORPUS | text-layer |
| roberts-tsoukkas-2026-vsgc-multiorbiter-cyclers-student-summary.pdf | 2026-06-13-roberts-tsoukkas-2026-multi-orbiter-mining.md | VSGC multi-orbiter cyclers (student summary, companion) | mined | text-layer |
| wittal-2022-earth-moon-cycler-lunar-logistics-IAC-22-C1.6.6.pdf | 2026-06-11-wittal-2022-iac-mining.md | E-M cycler for lunar logistics | mined | text-layer |
| wittal-smith-cassell-2021-robotic-lunar-gateway-payload-return-AAS-21-724.pdf | 2026-06-14-wittal-2021-aas-21724-digest.md | Robotic lunar Gateway payload return | digested | text-layer |
| AAS-22-015-pascarella-pony-express.pdf | OUTSTANDING H.3 (#38); 2026-06-11-forward-citation-sweep.md | Pony Express low-thrust E-M cycler | mined-by-catalogue (Pony Express rows) | text-layer |
| sanchez-net-2022-cycler-orbits-solar-system-pony-express-JSR.pdf | catalogue sanchez-net-2022-* rows; OUTSTANDING H.4 (#38) | Cycler orbits across solar system (Pony Express) | mined-by-catalogue | text-layer |
| merrill-2025-low-thrust-forced-periodic-em-cr3bp-arxiv-2502.05140.pdf | 2026-06-13-merrill-2502-05140-lowthrust-forced-periodic-mining.md | Low-thrust forced periodic orbits in E-M CR3BP | mined | text-layer |
| hiraiwa-2026-lobe-dynamics-cislunar-transfers-arxiv-2602.17444.pdf | 2026-06-07-hiraiwa-lobe-dynamics-method-mining.md; KNOWN_CORPUS | Lobe-dynamics cislunar transfers | mined + KNOWN_CORPUS | text-layer |
| singh-anderson-taheri-2021-exploiting-manifolds-L1-halo-em-low-thrust-...j.actaastro.2021.03.017.pdf | 2026-06-17-digest-singh-2021-L1-halo.md; KNOWN_CORPUS | Exploiting L1-halo manifolds for low-thrust E-M | digested + KNOWN_CORPUS | text-layer |
| singh-anderson-taheri-2021-low-thrust-transfers-southern-L2-NRHO-invariant-manifolds-jota.pdf | 2026-06-17-digest-singh-2021-NRHO.md; KNOWN_CORPUS | Low-thrust transfers to southern L2 NRHO | digested + KNOWN_CORPUS | text-layer |
| cuevas-del-valle-2023-optimal-floquet-stationkeeping-relative-dynamics-three-body-aerospace.pdf | 2026-06-11-cuevas-del-valle-2023-floquet-mining.md | Optimal Floquet stationkeeping in CR3BP | mined | text-layer |
| cuevas-del-valle-2026-fuel-optimal-rendezvous-CR3BP-MPC-proximal-EuroGNC.pdf | 2026-06-10-cuevas-del-valle-2026-cr3bp-mpc-mining.md; 2026-06-13-cuevas-2026-l1-halo-seed-run.md | Fuel-optimal CR3BP rendezvous MPC | mined | text-layer |

## CR3BP / periodic-orbit families / three-body networks

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| restrepo-russell-2018-database-planar-axisymmetric-periodic-orbits-solar-system-CMDA.pdf | 2026-06-17-digest-restrepo-russell-2018.md | Database of planar axisymmetric periodic orbits | digested | text-layer |
| russell-2012-survey-spacecraft-trajectory-design-strongly-perturbed-environments-jgcd.pdf | 2026-06-17-digest-russell-2012.md | Survey of trajectory design in strongly-perturbed regimes | digested | text-layer |
| braik-ross-2026-orbital-networks-three-body-problem-arxiv-2605.31543.pdf | 2026-06-13-braik-ross-2026-orbital-networks-mining.md; KNOWN_CORPUS | Orbital networks in the three-body problem | mined + KNOWN_CORPUS | text-layer |
| 2025-fixed-points-three-body-high-order-transfer-map-arxiv-2509.12671.pdf | 2026-06-13-high-order-transfer-map-2509.12671-mining.md | Fixed points in 3BP via high-order transfer map | mined | text-layer |
| fu-2026-datamining-escape-families-arxiv-2601.11881.pdf | 2026-06-13-fu-2601-11881-datamining-escape-family-mining.md | Data-mining escape families | mined | text-layer |
| tagliaferri-2024-mbh-manifold-transfers-arxiv-2405.18916.pdf | 2026-06-13-tagliaferri-2405-18916-global-opt-manifolds-mining.md | MBH global-opt with manifold transfers | mined | text-layer |
| andreu-1998-quasi-bicircular-problem-phd-thesis.pdf | 2026-06-14-andreu-quasi-bicircular-digest.md | Quasi-bicircular problem PhD thesis (POL1/POL2 ICs) | digested | text-layer |
| andreu-1998-quasi-bicircular-problem-phd-thesis.ps.gz | (PostScript source of the PDF above) | Original PS of the andreu thesis — same content | digested (via PDF) | n/a (compressed PS) |
| de-la-fuente-marcos-2018-geometric-characterization-arjuna-orbital-domain-arxiv-1410.4104v2.pdf | 2026-06-17-digest-fuente-marcos-2018.md | Geometric characterization of the Arjuna orbital domain | digested | text-layer |
| Wilczak-Zgliczynski-2002-Heteroclinic-Part-I.pdf | (duplicate of wilczak-zgliczynski-2003) | Same paper, different date/version | digested + KNOWN_CORPUS | text-layer |
| wilczak-zgliczynski-2003-heteroclinic-connections-part-I-math-0201278.pdf | 2026-06-19-digest-wilczak-zgliczynski-oterma-heteroclinic.md; KNOWN_CORPUS | W-Z Part I: L1<->L2 heteroclinic cycle proof | digested + KNOWN_CORPUS | text-layer |
| Wilczak-Zgliczynski-2005-Heteroclinic-Part-II.pdf | (duplicate of wilczak-zgliczynski-2006) | Same paper, different date/version | digested + KNOWN_CORPUS | text-layer |
| wilczak-zgliczynski-2006-heteroclinic-connections-part-II-math-0401146.pdf | 2026-06-19-digest-wilczak-zgliczynski-oterma-heteroclinic.md; KNOWN_CORPUS | W-Z Part II: rigorous numerics for the proof | digested + KNOWN_CORPUS | text-layer |
| canalias-2007-thesis-mission-design.pdf | 2026-06-20-digest-canalias-2007-se-em-manifolds.md | SE-EM manifold connections and heteroclinic bifurcations (PhD thesis) | digested | text-layer |

## Outer-planet / moon tours & endgame (Galilean, Saturnian, Uranian)

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| niehoff-1970-touring-galilean-satellites-AIAA-paper-70-1070.pdf | 2026-06-17-digest-niehoff-1970.md; KNOWN_CORPUS | Touring the Galilean satellites | digested + KNOWN_CORPUS | text-layer |
| campagnola-russell-2009-endgame-partA-vinf-leveraging-graph-AAS-09-224.pdf | 2026-06-05-endgame-tisserand-mining.md | Endgame Part A (V∞-leveraging graph) | mined | text-layer |
| campagnola-russell-2009-endgame-partB-multibody-tp-graph-AAS-09-227.pdf | 2026-06-05-endgame-tisserand-mining.md | Endgame Part B (multibody Tisserand-Poincaré graph) | mined | text-layer |
| campagnola-buffington-petropoulos-2014-jovian-tour-design-europa-orbiter-lander-...actaastro.pdf | 2026-06-17-digest-campagnola-2014.md; KNOWN_CORPUS | Jovian tour design (Europa orbiter/lander) | digested + KNOWN_CORPUS | text-layer |
| strange-russell-buffington-2007-mapping-v-infinity-globe-AAS-07-277.pdf | KNOWN_CORPUS (Mapping the V∞ globe) | Mapping the V∞ globe (multi-body resonance hopping) | mined-by-KNOWN_CORPUS | text-layer |
| heaton-longuski-2003-feasibility-galileo-style-tour-uranian-satellites-jsr-doi-10.2514-2.3981.pdf | catalogue Heaton-Longuski U00-01 (Table 5 transcribed); KNOWN_CORPUS | Galileo-style Uranian satellite tour feasibility | mined-by-catalogue + KNOWN_CORPUS | image-only → OCR'd .txt (2026-06-22) |
| davis-phillips-mccarthy-2018-saturnian-ocean-worlds-poincare-maps-...j.actaastro.2017.11.004.pdf | 2026-06-17-346-davis-2018-deep-read.md; KNOWN_CORPUS | Saturnian ocean-worlds tour via Poincaré maps | digested + KNOWN_CORPUS | text-layer |
| takao-2025-mission-analysis-first-saturn-trojan-2019-uo14-arxiv-2501.06586.pdf | 2026-06-07-takao-2025-mpga-1dsm-mining.md | First Saturn-Trojan (2019 UO14) MPGA+1DSM mission | mined | text-layer |
| liang-2024-callisto-ganymede-europa-triple-cyclers-JGCD.pdf | 2026-06-11-liang-2024-cge-triple-cyclers-mining.md; KNOWN_CORPUS; catalogue | Callisto-Ganymede-Europa triple cyclers | mined + mined-by-catalogue + KNOWN_CORPUS | text-layer |
| jones-hernandez-jesick-2017-low-excess-speed-vem-triple-cyclers-AAS-17-577.pdf | 2026-06-05-jones-aas17-577-vem-mining.md; KNOWN_CORPUS; catalogue | Low-excess-speed VEM triple cyclers | mined + mined-by-catalogue + KNOWN_CORPUS | text-layer |
| hernandez-jones-jesick-2017-one-class-io-europa-ganymede-triple-cyclers-AAS-17-608.pdf | 2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md; catalogue | One-class Io-Europa-Ganymede triple cyclers (AAS 17-608; Jovian sibling of VEM 17-577) | digested + catalogued | text-layer |
| russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-doi-10.2514-1.36610.pdf | 2026-06-30-digest-russell-strange-2009-planetary-moon-cyclers.md | THE canonical moon-cycler census (free-return cyclers; Jovian Ganymede/Europa-flyby pairs + Saturnian Titan→Enceladus ONLY — corrects the #320 "any two moons" over-claim) | digested | text-layer |
| lynam-longuski-2011-laplace-resonant-triple-cyclers-jupiter-acta-astronautica-69-doi-10.1016-j.actaastro.2011.03.011.pdf | 2026-06-30-digest-lynam-longuski-2011-laplace-resonant-triple-cyclers.md | Laplace-resonant IEG triple-cyclers (the #480 prior; Hernandez 2017 ref [7]; 2nd independent IEG source) | digested | text-layer |
| vasile-campagnola-2009-lowthrust-mga-europa-JBIS-arxiv-1105.1823.pdf | 2026-06-07-vasile-campagnola-dfet-method-mining.md; 2026-06-05-vasile-tables-retranscription.md | Low-thrust MGA to Europa (DFET) | mined | text-layer |
| genova-2016-phobos-deimos-PADME-trajectory-AIAA-2016-5681.pdf | 2026-06-05-v42-backfill-sweep.md §7 (PADME) | Phobos/Deimos PADME trajectory | mined | text-layer |

## Historic / reference flight missions (Voyager, Mariner, Galileo, Cassini, Juno)

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| bourke-friedman-penzo-1971-design-grand-tour-missions-AIAA-71-187.pdf | 2026-06-17-digest-bourke-1971.md | Design of Grand Tour missions | digested | text-layer |
| kohlhase-penzo-1977-voyager-mission-description-space-science-reviews-21-77-101.pdf | 2026-06-19-345-voyager-mariner-mission-digests.md; catalogue | Voyager mission description | digested + mined-by-catalogue | text-layer |
| mckinley-vanallen-1976-mariner-jupiter-saturn-1977-navigation-strategy-jsr-...57113.pdf | 2026-06-19-345-voyager-mariner-mission-digests.md; catalogue | Mariner Jupiter/Saturn 1977 navigation strategy | digested + mined-by-catalogue | text-layer |
| dunne-burgess-1978-voyage-mariner-10-mission-venus-mercury-NASA-SP-424.pdf | 2026-06-17-digest-dunne-burgess-1978-mariner-10.md | Mariner 10 mission (Venus/Mercury) NASA SP-424 | digested | text-layer |
| giberson-cunningham-1975-mariner-10-mission-venus-mercury-...0094-5765(75)90012-0.pdf | 2026-06-19-345-voyager-mariner-mission-digests.md | Mariner 10 mission to Venus and Mercury | digested | text-layer |
| johnson-yeates-young-1992-galileo-mission-overview-space-science-reviews-vol-60.pdf | 2026-06-17-digest-johnson-yeates-young-1992.md | Galileo mission overview | digested | text-layer |
| damario-bright-wolf-1992-galileo-trajectory-design-space-science-reviews-...bf00216849.pdf | 2026-06-17-digest-damario-1992-galileo.md; catalogue | Galileo trajectory design (VEEGA) | digested + mined-by-catalogue | text-layer |
| diehl-kaplan-penzo-1983-satellite-tour-design-galileo-mission-AIAA-83-0101.pdf | 2026-06-17-digest-diehl-1983.md; catalogue | Galileo satellite-tour design | digested + mined-by-catalogue | text-layer |
| damario-byrnes-1983-interplanetary-trajectory-design-galileo-mission-AIAA-83-0099.pdf | 2026-06-19-digest-damario-byrnes-1983-galileo-interplanetary.md | D'Amario-Byrnes 1986 DIRECT E→J Galileo interplanetary leg (pre-Challenger, never flown; companion to Diehl-Kaplan-Penzo 0101); Jupiter VHP 5.8/6.1 km/s, no GA flybys; Galileo errata corroborator, not catalogue-admissible (#384) | digested | text-layer |
| young-1998-galileo-probe-mission-jupiter-science-overview-jgr-doi-10.1029-98JE01051.pdf | 2026-06-17-digest-damario-1992-galileo.md (probe context); catalogue | Galileo probe mission science overview | digested + mined-by-catalogue | text-layer |
| young-2000-correction-galileo-probe-mission-jupiter-science-overview-jgr-...2000JE001251.pdf | (errata to young-1998 above) | Correction to the Galileo probe overview | digested (via 1998) | image-only → OCR'd .txt (2026-06-22; 105-word erratum) |
| wolf-smith-1995-design-cassini-tour-trajectory-saturnian-system-...0967-0661(95)00172-7.pdf | 2026-06-17-digest-wolf-smith-1995-cassini.md; catalogue | Cassini tour trajectory design | digested + mined-by-catalogue | text-layer |
| bellerose-roth-wagner-2018-cassini-reconstructing-thirteen-years-...2018-2646.pdf | 2026-06-17-digest-bellerose-2018-cassini.md | Cassini: reconstructing 13 years of gravity assists | digested | text-layer |
| valerino-2014-updating-reference-trajectory-cassini-solstice-mission-spaceops-...2014-1880.pdf | 2026-06-17-digest-bellerose-2018-cassini.md (Solstice context); KNOWN_CORPUS | Updating Cassini Solstice reference trajectory | digested + KNOWN_CORPUS | text-layer |
| yam-davis-longuski-2009-saturn-impact-trajectories-cassini-end-of-mission-jsr-...38760.pdf | 2026-06-04-agrawal-landau-howe-mining.md; catalogue | Cassini end-of-mission Saturn-impact trajectories | mined + mined-by-catalogue | text-layer |
| lam-johannesen-kowalkowski-2008-planetary-protection-trajectory-analysis-juno-...2008-7368.pdf | 2026-06-17-digest-lam-2008-juno.md | Juno planetary-protection trajectory analysis | digested | text-layer |
| stone-miner-1986-voyager-2-encounter-uranian-system-science-233-4759.pdf | 2026-06-23-digest-stone-miner-voyager2-uranus-neptune.md | Voyager 2 Uranus encounter overview (Science 233:39); Uranus flown C/A 107,000 km from centre = 81,441 km alt; Miranda mass-pass 28,260 km (#429) | digested | text-layer |
| stone-miner-1989-voyager-2-encounter-neptunian-system-science-246-1417.pdf | 2026-06-23-digest-stone-miner-voyager2-uranus-neptune.md | Voyager 2 Neptune encounter overview (Science 246:1417); Neptune flown C/A 29,240 km from centre = 4,476 km alt; Triton C/A 39,800 km = 38,447 km alt (#429) | digested | text-layer |
| stern-tapley-finley-scherrer-2020-pluto-orbiter-kuiper-belt-explorer-jsr-a34658.pdf | 2026-06-23-digest-stern-2020-pluto-orbiter.md | Pluto Orbiter-KB Explorer mission-design study (JSR A34658); Charon-GA Pluto-system tour; Pluto design floor 100 km SOURCED (atmospheric periapse 100-500 km); Charon periapsis deferred to Finley [5] (#429) | digested | text-layer |
| harch-bhaskaran-2017-accommodating-navigation-uncertainties-pluto-encounter-sequence-design-space-operations-springer.pdf | 2026-06-23-digest-harch-2016-pluto-encounter-sequence.md | New Horizons Pluto-flyby encounter-sequence nav-uncertainty methodology (SpaceOps 2016, DOI 10.2514/6.2016-2623); flown Pluto C/A 13,695 km→12,507 km alt, Charon C/A 29,432 km→28,826 km alt (observed-flown, #429) | digested | text-layer |
| hollister-menning-1970-periodic-swingby-earth-venus-JSR-7-10.pdf | catalogue hollister-menning-1970-ev-orbit-01..15 family | Earth-Venus periodic swingby (E-V cycler family) | mined-by-catalogue | image-only → OCR'd .txt (2026-06-22) |
| hollister-rall-1970-periodic-orbits-NASA-CR.pdf | 2026-06-07-hollister-rall-1970-periodic-orbits-mining.md; -appendices-transcription | Periodic orbits NASA-CR (appendix tables transcribed) | mined | text-layer |

## Methods: optimization, STM, primer-vector, low-thrust, reachable sets

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| guzman-mailhe-schiff-hughes-folta-2002-primer-vector-optimization-survey-IAC-02-A.6.09.pdf | 2026-06-07-guzman-2002-primer-survey-mining.md | Primer-vector optimization survey | mined | text-layer |
| damario-byrnes-stanford-1981-new-method-optimizing-multiple-flyby-trajectories-jgc-doi-10.2514-3.56115.pdf | 2026-06-19-digest-damario-byrnes-stanford-1981-moses-multiconic.md | Canonical JPL multiple-flyby optimizer (AIAA 80-1676R): bounds-constrained Newton, analytic gradient+Hessian, weighted SOS-of-ΔV cost, multiconic (1STEP/MULCON) + analytic STM; corrector/optimizer-lane provenance (#380/#347) (#384) | digested | text-layer |
| ellison-2018-analytic-gradient-bounded-impulse-trajectory-two-sided-shooting-JGCD.pdf | 2026-06-10-ellison-2018-analytic-gradients-mining.md | Analytic-gradient bounded-impulse two-sided shooting | mined | text-layer |
| iorfida-2016-geometric-perpendicular-thrust-trajectory-optimization-JGCD.pdf | 2026-06-10-iorfida-2016-perpendicular-thrust-mining.md | Geometric perpendicular-thrust optimization | mined | text-layer |
| shakouri-2019-shape-based-multiple-impulse-coplanar-maneuvers-arxiv.pdf | 2026-06-10-shakouri-2019-shape-based-mining.md | Shape-based multiple-impulse coplanar maneuvers | mined | text-layer |
| junkins-taheri-2019-alternative-state-vector-choices-low-thrust-...jgcd-doi-10.2514-1.G003686.pdf | 2026-06-17-digest-junkins-taheri-2019.md | Alternative state-vector choices for low-thrust opt | digested | text-layer |
| woollands-taheri-junkins-2019-efficient-computation-optimal-low-thrust-gravity-perturbed-jas.pdf | 2026-06-17-digest-woollands-2019.md | Efficient low-thrust optimal control (gravity-perturbed) | digested | text-layer |
| pellegrini-russell-2016-computation-accuracy-trajectory-state-transition-matrices-jgcd.pdf | 2026-06-17-digest-pellegrini-russell-2016.md | Computation/accuracy of trajectory STMs | digested | text-layer |
| saloglu-2023-infinitely-many-optimal-iso-impulse-trajectories-JGCD.pdf | 2026-06-10-saloglu-2023-iso-impulse-mining.md | Infinitely many iso-impulse trajectories | mined | text-layer |
| saloglu-2025-iso-impulse-3d-classification-feasibility-arxiv.pdf | 2026-06-10-saloglu-2025-iso-impulse-3d-mining.md | Iso-impulse 3D classification/feasibility (preprint) | mined | text-layer |
| saloglu-taheri-2025-iso-impulse-3d-JAS-vor-s40295-025-00528-0.pdf | 2026-06-10-saloglu-2025-iso-impulse-3d-mining.md (VoR of the preprint above) | Iso-impulse 3D (JAS version-of-record) | mined | text-layer |
| zhou-armellin-2025-single-impulse-reachable-set-polynomials-arxiv-2502.11280.pdf | 2026-06-07-zhou-2025-da-reachable-sets-mining.md | Single-impulse reachable-set polynomials (DA) | mined | text-layer |
| beeson-englander-hughes-2015-emtg-gmat-lowthrust-tool-chain-AAS-15-278.pdf | 2026-06-07-beeson-2015-emtg-gmat-toolchain-mining.md | EMTG/GMAT low-thrust toolchain | mined | text-layer |
| englander-englander-2014-tuning-monotonic-basin-hopping-ISSFD24-S7-3.pdf | 2026-06-07-englander-2014-mbh-tuning-mining.md | Tuning monotonic basin hopping (MBH) | mined | text-layer |
| ozimek-2019-linx-lowthrust-mga-trajectory-optimization-AAS-19-348.pdf | 2026-06-07-ozimek-linx-aas19-348-mining.md | LinX low-thrust MGA optimization | mined | text-layer |
| shepperd-1985-universal-keplerian-state-transition-matrix-celest-mech-35.pdf | 2026-06-19-digest-shepperd-1985-universal-keplerian-stm.md | Closed-form universal-variables (Goodyear/Sundman/Stumpff) two-body STM valid for all conics; Kepler solve via one transcendental (q-variable, Newton dt/du=4(1-q)r), STM from 9 M-coefficients + W,U; canonical analytic STM for conic-leg correctors (cf. Pellegrini-Russell #372) | digested | image-only → OCR'd .txt (2026-06-22; was thin) |
| montenbruck-markgraf-2004-gps-sensor-impact-point-prediction-sounding-rockets-jsr-...1962.pdf | 2026-06-17-digest-montenbruck-2004.md | GPS impact-point prediction for sounding rockets | digested | text-layer |
| rinker-jacobson-wood-1976-statistical-analysis-trim-maneuvers-low-thrust-interplanetary-navigation-jsr.pdf | 2026-06-19-345-voyager-mariner-mission-digests.md — off-scope for #345 | Statistical analysis of nav trim maneuvers | triaged | text-layer |

## ML / surrogate / GNC background (triaged sweep)

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| zhang-2024-neural-angle-only-od-earth-moon-libration-remotesensing-16-03287.pdf | 2026-06-07-zhang-2024-neural-od-mining.md | Neural angle-only OD at E-M libration | mined | text-layer |
| zhang-acciarini-2026-pretrained-approximators-lowthrust-cost-reachability-arxiv-2605.26790.pdf | 2026-06-07-ml-surrogate-investigation.md Paper 1; background-triage #8 | Pretrained approximators for low-thrust cost/reachability | mined/triaged | text-layer |
| zhang-2026-neural-porkchop-lowthrust-asteroid-rendezvous-astronautics-01-00006.pdf | 2026-06-07-ml-surrogate-investigation.md Paper 2 | Neural porkchop low-thrust asteroid rendezvous | mined | text-layer |
| ozaki-2022-neural-network-surrogate-global-cycler-search-arxiv-2111.11858.pdf | 2026-06-11-ml-surrogate-trio-triage.md — blueprint DEFERRED | NN surrogate for global cycler search | triaged | text-layer |
| leifsson-2022-global-surrogate-modeling-nn-uncertainty-ICCS.pdf | 2026-06-11-ml-surrogate-trio-triage.md — background-only | Global surrogate modeling, NN uncertainty | triaged | text-layer |
| wu-2024-physics-informed-ml-review-condition-monitoring-ESWA.pdf | 2026-06-11-ml-surrogate-trio-triage.md — background-only | Physics-informed ML review (condition-monitoring) | triaged | text-layer |
| viavattene-ceriotti-2021-neural-multiple-nea-rendezvous-continuous-thrust-JSR.pdf | 2026-06-07-ml-surrogate-investigation.md Paper 3; background-triage #7 — OUT-OF-SCOPE | ANN multiple-NEA rendezvous (continuous thrust) | triaged | text-layer |
| silvestrini-lavagna-2022-deep-learning-ann-spacecraft-gnc-drones-6-270.pdf | 2026-06-13-background-papers-read-triage.md #1 — OUT-OF-SCOPE | Deep-learning ANN for spacecraft GNC | triaged | text-layer |
| li-topputo-baoyin-2019-neural-time-optimal-orbit-raising-ep-geo-arxiv-1909.08768.pdf | 2026-06-13-background-papers-read-triage.md #2 — OUT-OF-SCOPE | Neural time-optimal EP-to-GEO orbit raising | triaged | text-layer |
| singh-junkins-2022-stochastic-learning-extremal-field-lowthrust-guidance-sci-rep.pdf | 2026-06-13-background-papers-read-triage.md #3 — OUT-OF-SCOPE | Stochastic learning of extremal field, low-thrust | triaged | text-layer |
| hu-yang-li-2024-robust-lowthrust-gravity-assist-rl-AIAA-G009427.pdf | 2026-06-07-marginal-papers-triage.md #5 — OUT-OF-SCOPE | Robust low-thrust GA via RL | triaged | text-layer |
| blender-singh-2025-uncertainty-aware-guidance-gbdt-continuous-thrust-AAS-25-524.pdf | 2026-06-07-marginal-papers-triage.md #6 — OUT-OF-SCOPE | Uncertainty-aware GBDT continuous-thrust guidance | triaged | text-layer |
| sinha-beeson-2025-initial-guess-lowthrust-missed-thrust-robust-arxiv-2501.06694.pdf | 2026-06-07-marginal-papers-triage.md #2 — REFERENCE-ONLY | Initial-guess generation, robust low-thrust | triaged | text-layer |
| venigalla-englander-scheeres-2020-lowthrust-missed-thrust-recovery-margin-AAS-20-438.pdf | 2026-06-07-marginal-papers-triage.md #1 — OUT-OF-SCOPE | Low-thrust missed-thrust recovery margin | triaged | text-layer |
| fan-2025-electric-sail-multi-target-trajectory-design-aerospace-12-00196.pdf | 2026-06-13-background-papers-read-triage.md #14 — OUT-OF-SCOPE (E-sail) | E-sail multi-target trajectory design | triaged | text-layer |

## Continuous-thrust / electric-sail / misc applied

| File | Pointer | One-line summary | Status | OCR |
|---|---|---|---|---|
| rickman-intro-orbital-mechanics-spacecraft-attitudes-thermal-engineers-NASA-NESC-slides.pdf | 2026-06-10-rickman-nesc-slides-triage.md | Intro orbital-mechanics lecture slides (NESC) | triaged | text-layer |
| acton-1996-ancillary-data-services-NASA-NAIF-SPICE-planet-space-sci-...95)00107-7.pdf | 2026-06-17-digest-acton-1996.md; catalogue | SPICE/NAIF ancillary data services | digested + mined-by-catalogue | text-layer |
| ccsds-2023-orbit-data-messages-502.0-B-3-blue-book.pdf | 2026-06-05-ccsds-odm-502-mining.md | CCSDS Orbit Data Messages standard | mined | text-layer |
| vergaaij-2018-time-optimal-solar-sail-heteroclinic-connections.pdf | 2026-06-20-digest-new-papers.md | Time-optimal solar sail heteroclinic connections | triaged | text-layer |
| antoniadou-libert-2019-spatial-resonant-periodic-orbits-RTBP-mnras-sty3195-arxiv-1811.09442.pdf | 2026-06-25-digest-antoniadou-libert-2019-spatial-resonant.md; known_corpus_3d | Spatial resonant periodic orbits in the RTBP (mu=0.001 planetary MMRs; spatial families born at vertical-critical-orbits). TAXONOMY anchor for the #434/#438 3D-novelty gate — no Earth-Moon ICs. Corrects the prior known_corpus_3d "Antoniadou-Voyatzis" mislabel (co-author is Libert) | digested | text-layer |
| antoniadou-libert-2018-origin-continuation-resonant-periodic-orbits-circular-elliptic-RTBP-cmda-s10569-018-9834-8-arxiv-1805.00288.pdf | 2026-06-25-digest-antoniadou-libert-2018-circ-ellip-continuation.md; ER3BP erratum | Circular+elliptic resonant PO continuation. KEY: documents ISOLATED high-e families with NO circular limit → reframes the #432/#435/#436 e>0 negatives as method-blind-spot (docs/superpowers/plans/2026-06-25-er3bp-negatives-erratum-antoniadou-libert.md) | digested | text-layer |
| howell-1984-three-dimensional-periodic-halo-orbits-celest-mech-32-53-doi-10.1007-BF01358403.pdf | 2026-06-25-digest-howell-1984-halo-orbits.md | Canonical CR3BP halo-family survey (L1/L2/L3, mu in (0,1)); 3 golden IC tables at mu=0.04/0.96 (NOT Earth-Moon — those trace to Breakwell & Brown 1979). Halo taxonomy anchor for #434/#438 | digested | text-layer |
| peng-xu-2017-continuation-periodic-orbits-sun-mercury-ERTBP-cnsns-47-1.pdf | 2026-06-25-digest-peng-2017-sun-mercury-ERTBP.md; ER3BP erratum | Peng-Bai-Xu 2017 Sun-Mercury ERTBP multi-rev resonant-Halo continuation; folds/disconnected branches/isolated loops in e (narrows #435 to its low-rev seed scope; #437 same-model golden Tables 2/3) | digested | text-layer |
| 2025-planar-retrograde-periodic-orbits-ERTBP-acta-astronautica-229-430-S0094576525000086.pdf | 2026-06-25-digest-planar-retrograde-ERTBP-2025.md; ER3BP erratum | Martinez-Cacho-Gil-Bombardelli-Baresi 2025 planar retrograde QSO=DRO families in the ERTBP; all circular-rooted (corroborates the #435 DRO negative; #437 fold golden at e=0.0324) | digested | text-layer |

## Gaps — `undigested-unmined` + image-only OCR backlog

### Genuine gaps: NONE remaining (all 3 closed 2026-06-19)

The last three `undigested-unmined` files were all digested on 2026-06-19,
bringing the corpus to **zero genuine gaps** (every `papers/` file is now
digested, mined, or triaged):

1. **shepperd-1985-...celest-mech-35.pdf** → digested
   `2026-06-19-digest-shepperd-1985-universal-keplerian-stm.md` (#402; full
   16-page read; universal-variable two-body STM). The thin text-layer was
   handled by reading the rendered page images directly.
2. **bond-allman-2021-...textbook.pdf** → digested
   `2026-06-19-digest-bond-allman-2021-modern-astrodynamics.md` (#396;
   chapter-summary scope; Sundman/Sperling-Burdet regularization Ch 9 feeds
   the tulip Sundman work, Table 9.3 metre-level L4/L5 case usable as a
   CR3BP golden).
3. **willis-2008-...asr.2007.07.047.pdf** → digested
   `2026-06-19-digest-willis-2008-book-review.md` (2-page review of the Gurfil
   2007 volume; corroborating metadata only, zero catalogue impact — defer to
   the #395 Gurfil digest).

### Image-only OCR backlog — CLEARED 2026-06-22 (0 remaining)

A corpus-wide `pdftotext` sweep (2026-06-22) found 6 image-only / near-zero-text
PDFs with no committed `.txt`. **All 6 were OCR'd (`ocrmypdf --force-ocr`) to
committed `.txt` sidecars** — the corpus is now fully text-searchable, and the
text survives a fresh clone / `git clean` (the prior #400 OCRs went to a
throwaway astropy cache that never persisted — so the index's "OCR done" claims
were false-confident; now made real):

- **szebehely-1967-...book.pdf** — 661pp textbook → 217 676 words.
- **heaton-longuski-2003-...uranian-satellites-...pdf** — 5 475 words.
- **okutsu-longuski-2002-mars-free-returns-...pdf** — 6 062 words.
- **hollister-menning-1970-periodic-swingby-earth-venus-...pdf** — 5 331 words.
- **shepperd-1985-...keplerian-stm-...pdf** — 4 671 words (was a thin/garbled layer).
- **young-2000-correction-...galileo-probe-...pdf** — 105 words (short erratum).

Also committed 8 previously-uncommitted text-extract sidecars (chen ×4,
klmr-2006, machado-2020, shoot-the-moon-2000, vergaaij-2018). Per policy §1 the
OCR text is for **navigation**; precision values (equations, table cells) must
still be vision-read from the page images. No image-only black box remains.

### Recent Acquisitions (Cassini)

- **peralta-flanagan-1995-cassini-interplanetary-trajectory-design-cep.pdf** | undigested | Cassini interplanetary leg | - | text-layer |
- **wolf-smith-1995-design-cassini-tour-trajectory-saturnian-system-control-eng-practice-doi-10.1016-0967-0661(95)00172-7.pdf** | digested | Cassini tour trajectory | KNOWN_CORPUS | text-layer |
- **wolf-1996-touring-saturnian-system-spie.pdf** | 2026-06-20-digest-wolf-1996.md | Cassini touring | - | text-layer |
- **canales-howell-2023-arxiv-2308.10029.pdf** | undigested | Transport corridor reproduction | - | text-layer |
- **kumar-2025-arxiv-2509.12675.pdf** | undigested | Resonant orbits (r41u, 3:1 to 2:1 heteroclinic) | - | text-layer |
