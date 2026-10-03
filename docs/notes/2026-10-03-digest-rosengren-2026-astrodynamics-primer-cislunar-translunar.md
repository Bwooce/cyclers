# Digest — Rosengren, Rawat, Kumar & Ross (2026), "The Astrodynamics Primer on Cislunar and Translunar Space" (arXiv:2606.26367v1)

**Digested:** 2026-10-03 (text-layer PDF, 77 pages; text read in full, Tables 1, 2 and 4 checked
against the rendered pages). **Citation:** Aaron J. Rosengren, Anjali Rawat, Bhanu Kumar, Shane
D. Ross, "The Astrodynamics Primer on Cislunar and Translunar Space," arXiv:2606.26367v1
[astro-ph.EP], 24 June 2026 (preprint). Filed in the private paper corpus as
`rosengren-rawat-kumar-ross-2026-astrodynamics-primer-cislunar-translunar-space-arxiv-2606.26367v1.pdf`
(77 pages, md5 `452de66aa449f75482eb127fe8b6ebf7`). **Corpus index description (one line):**
Long review partitioning Earth-Moon space into terrestrial / cislunar / circumlunar /
translunar regimes (Laplace radius, lunar and solar mean-motion resonance ladders, gateway and
Hill boundaries) with MEGNO and fate-class maps; one short section (6.4) on Earth-Moon cyclers.

## What it is
A review and synthesis (no new trajectory-design method). Its "spatiography" patches two
restricted problems (Earth-Moon, and Sun-(Earth+Moon)) and combines orbit-averaged
perturbation theory, Gallardo-style semi-analytical resonance atlases, restricted-problem
periodic-orbit and manifold structure, and numerical MEGNO and fate-class maps. Its argument is
that "cislunar" is used loosely for regions that are dynamically distinct.

## Section list
1 Introduction; 2 Historical Context and Persistent Misconceptions (2.1 lexicon, 2.2 positional
meaning of cislunar, 2.3 sociolinguistic drift, 2.4 forgotten Space-Age taxonomy, 2.5
Spatiography and the patched-CR3BP viewpoint, 2.6 This Primer's usage); 3 The Curated Cislunar
and Translunar Catalog (3.1 perturbation hierarchy, 3.2 Hill-region topology, 3.3 catalog in
the synodic frame); 4 Perturbative Treatment of Distant Geocentric Orbits (4.1 perturbed
Hamiltonian, with disturbing-function subsections; 4.2 von Zeipel-Lidov-Kozai; 4.3 lunisolar
secular resonances; 4.4 lunar and solar mean-motion resonances; 4.5 secular and resonant
structures in the catalog); 5 Spatiography (5.1 Laplace radius; 5.2 inner cislunar zone; 5.3
outer cislunar zone, lunar MMRs; 5.4 spheres of influence, Hill scales and Jacobi gateways;
5.5 circumlunar space; 5.6 lunisolar tidal parity; 5.7 translunar realm and exterior
resonances); 6 Restricted Multi-Body Dynamical Systems Theory (6.1 phenomenological models; 6.2
resonant dynamics and onset of chaos; 6.3 stability and bifurcations of periodic-orbit
families; **6.4 Earth-Moon Cyclers and Collision Orbits**); 7 Astro-Cartographies (7.1 MEGNO
with REBOUND; 7.2 fate classes; 7.3 spatiographic dynamical maps, Figs. 13-16); 8 Conclusions;
Appendix A (Table 3, etymological timeline of spatial terms); Appendix B (Table 4, integration
spans); References. Sixteen figures, four tables.

## Regime definitions and boundary values (as printed)
Units: a/a_Moon is the geocentric semi-major axis over the lunar value; T is the Earth-centred
Keplerian period. Constants (Sec. 5.1): R_Earth = 6378.1363 km, J2 = 1.08263552549e-3 (GGM02);
lunar mean elements from Simon et al. 1994: a_Moon = 383397.7725 km, e_Moon = 0.055545526,
I_Moon = 5.15668983 deg; Sun: a = 1.0000010178 au, e = 0.0167086342.
- **Terrestrial to cislunar:** geolunar Laplace radius r_L, where the combined lunisolar torque
  equals the oblateness-driven precession (Eq. 98); r_L ~ 7.7 R_Earth ~ 0.13 a_Moon (T = 1.24 d).
- **Secularly dominated zone:** 0.13 to 0.34 a_Moon (von Zeipel-Lidov-Kozai cycles, lunar
  nodal and apsidal secular resonances).
- **Cislunar resonant zone:** interior lunar mean-motion resonances from 5:1 (0.34 a_Moon)
  to 5:4 (0.86); "low-order" means |k - k_Moon| <= 4 (Sec. 5.3).
- **Circumlunar space:** Earth-Moon L1 at 0.84 a_Moon to L2 at 1.16 (inner and outer
  zero-velocity surfaces); local lunar vicinity treated as an enclave (Table 2).
- **Translunar:** exterior lunar resonances 4:5 (1.16) to 1:5 (2.92), first solar resonance 5:1
  at 1.93, lunisolar tidal parity 1.17 a_Moon (34.6 d), Laplace patched-conic SOI 2.41, Earth Hill
  sphere 3.90 (210.88 d).

**Table 1 (p. 26), transcribed** (resonance label k:k_b with k for the satellite and k_b for
the perturber; "Moon" = lunar commensurability, "Sun" = solar):

| Zone / feature | a/a_Moon | T [days] | Description as printed |
|---|---|---|---|
| Cislunar lower bound: geolunar Laplace radius r_L | 0.13 | 1.24 | where lunisolar torques ~ Earth's oblateness |
| Secularly dominated zone | 0.13-0.34 | 1.24-5.47 | orbits of OGO, HEOS, Vela, IMP, Prognoz, GEOTAIL, INTERBALL, CXO, XMM-Newton, Cluster, INTEGRAL, &c. |
| 5:1 Moon | 0.34 | 5.47 | innermost low-order lunar MMR |
| 4:1 Moon | 0.40 | 6.84 | |
| 3:1 Moon | 0.48 | 9.11 | orbit of IBEX and Tiandu-1 |
| 5:2 Moon | 0.54 | 10.94 | |
| 2:1 Moon | 0.63 | 13.67 | orbit of TESS |
| 5:3 Moon | 0.71 | 16.41 | |
| 3:2 Moon | 0.76 | 18.23 | orbit of DRO-B |
| 4:3 Moon | 0.83 | 20.51 | |
| 5:4 Moon | 0.86 | 21.88 | |
| Earth-Moon L1 | 0.84 | 20.94 | inner zero-velocity surface |
| Moon's orbit (1:1) | 1.00 | 27.34 | lunar semi-major axis |
| Earth-Moon L2 | 1.16 | 34.13 | outer zero-velocity surface |
| 4:5 Moon | 1.16 | 34.18 | |
| 3:4 Moon | 1.21 | 36.46 | |
| 2:3 Moon | 1.31 | 41.02 | |
| 3:5 Moon | 1.41 | 45.57 | |
| 1:2 Moon | 1.59 | 54.69 | |
| 2:5 Moon | 1.84 | 68.36 | |
| 5:1 Sun | 1.93 | 73.05 | innermost low-order solar MMR |
| 1:3 Moon | 2.08 | 82.00 | |
| 4:1 Sun | 2.23 | 91.31 | |
| 1:4 Moon | 2.52 | 109.38 | |
| 3:1 Sun | 2.71 | 121.75 | |
| 1:5 Moon | 2.92 | 136.72 | outermost low-order lunar MMR |
| 5:2 Sun | 3.06 | 146.10 | |
| 2:1 Sun | 3.55 | 182.63 | outermost low-order solar MMR |
| Lunisolar tidal parity | 1.17 | 34.6 | lunar-internal equals solar-external quadrupole tide (Eq. 127) |
| Laplace's patched-conic SOI | 2.41 | 102.41 | (r_SOI) = a (mu_E / mu_Sun)^(2/5) |
| Earth's Hill sphere | 3.90 | 210.88 | (r_H) = a (mu_E / 3 mu_Sun)^(1/3) |

**Table 2 (p. 39), circumlunar zone in selenocentric distance rho / R_Moon:** low lunar orbit
(h = 100 km) 1.06, T 0.08 d; selenoterrestrial Laplace radius 2.21, 0.25 d; exterior terrestrial
resonances 8:1 to 7:4 at 12.73, 13.92, 15.43, 17.42, 18.69, 20.22, 22.10, 22.83, 24.49, 26.49,
27.65, 28.96, 29.67, 32.09, 34.42, 35.08 (periods 3.42 to 15.63 d; the 8:1 is "innermost
7th-order terrestrial commensurability", 7:4 "outermost terrestrial resonance interior to
(rho_H)"); Chebotarev sphere 24.47 (9.11 d); Battin SOI Earthward 29.93 (12.32 d); Earth-Moon L1
33.31 (14.46 d); Moon's Hill sphere 35.32 (15.79 d, text gives 35.32 R_Moon = 61364 km); Battin
SOI anti-Earthward 36.95 (16.90 d); Earth-Moon L2 37.04 (16.95 d); Laplace patched-conic SOI
37.99 (17.61 d).

**Table 4 (p. 66), map domains** (a/a_Moon range; span in years; revolutions): Secularly
dominated cislunar 0.13-0.35, 19, 5419-1227; Cislunar resonant 0.33-0.89, 19, 1340-303;
Circumlunar gateway 0.84-1.16, 19, 330-203; Inner translunar resonant 1.08-2.03, 38, 453-176;
Outer translunar resonant 1.91-3.34, 57, 289-125; Translunar fringe 3.03-3.90, 57, 144-99.
Table 3 is an etymological timeline of spatial terms (no dynamical values).

Numerical maps (Sec. 7.3): elliptic point-mass Earth-Moon model (EM) and an Earth-Moon-Sun
point-mass model (EMS), both initialised from JPL Horizons at 2027 August 2, 10:06:37 UTC;
Omega = 311.07, omega = 355.84, M = 0 degrees, inclination in the Moon's orbital plane;
REBOUND with IAS15; MEGNO (Cincotta and Simo) plus fate classes. The EM model is explicitly
"rather than in the autonomous circular restricted three-body problem".

## Mean-motion resonances, resonant periodic orbits, heteroclinic pathways (Earth-Moon)
- **Interior lunar resonance ladder:** 5:1, 4:1, 3:1, 5:2, 2:1, 5:3, 3:2, 4:3, 5:4 (Table 1);
  Gallardo-atlas widths (Fig. 8, coplanar, a 2 R_Hill cutoff); "broadest and dynamically most
  conspicuous" are the 3:1 and 2:1 families, citing Rawat et al. 2026. IBEX near 3:1, TESS in
  2:1, Spektr-R neighbouring higher-order structure. The 4:3 and 5:4 are "the last clearly
  marked interior lunar resonances before the circumlunar transition".
- **Exterior (translunar):** Poincare maps at osculating apogee expose exterior 1:n lunar
  resonances with stable asymmetric libration zones coexisting with weak and strong symmetric
  unstable zones; separatrices from the unstable resonant orbits, and the Earth-Moon L1-L2 tube
  geometry, "mediate transitions between exterior and interior resonant domains" (Rawat et al.
  2025). Chirikov overlap of exterior lunar and interior solar resonances is used for the onset
  of large-scale chaos in the translunar region.
- **Heteroclinic pathways, verbatim (Sec. 6.2):** "Together with the heteroclinic connections
  among the interior 4:1, 3:1, and 2:1 resonances and the lunar (L1) region, these structures
  define a network of ballistic pathways by which the semimajor axis of a spacecraft may change
  substantially under lunar perturbations alone (Belbruno et al 2008; Liang et al 2017; Lei and
  Xu 2018; Peng et al 2024; Kumar et al 2026; Pan et al 2026)." No connection counts, Jacobi
  constants or transfer times are given in the Primer; it defers to the cited papers.
- **Periodic-orbit atlas (Sec. 6.3):** Lyapunov, vertical, halo (NRHO), DRO, L3 and resonant
  symmetric and asymmetric families, cited to Broucke 1968, Howell 2001, Doedel et al. 2007,
  and **Bonasera and Bosanac 2023** (the only mention of that paper in the body); stability
  index nu = (lambda + 1/lambda) / 2 (Eq. 141); resonant orbits described as finite-mu
  continuations of subharmonic Keplerian orbits. Elliptic problem: the Jacobi integral is lost
  and only commensurable circular-problem members continue at fixed eccentricity.
- No resonant-orbit period, Jacobi-constant or initial-condition tables are printed.

## Cyclers, repeated encounters, gravity assists (Sec. 6.4 "Earth-Moon Cyclers and Collision Orbits", pp. 51-52)
Counts, body text only (references excluded): "cycler" 21 occurrences (case-insensitive, includes
"cyclers" and "cycler-like" and the section 6.4 heading; they fall on pp. 2, 51 and 52);
"periodic orbit" 19 (plus "periodic-orbit" 15); "heteroclinic" 4; "homoclinic" 3; "resonan"
330; "flyby" 7; "gravity assist" 1; "network" 2 (the passage above, and "a wider network of stable islands, unstable periodic or quasi-periodic invariant objects, tube dynamics, whiskered tori, and chaotic transport" in Sec. 6.1 describing the CR3BP). Reference list
only: cycler 4, periodic orbit 19, heteroclinic 5, homoclinic 1, resonan 49, gravity assist 2.
Key passages, quoted:
- Definition: "repeated-encounter periodic orbits: trajectories that return successively to the
  neighborhoods of both primaries. In the Earth-Moon problem these are naturally interpreted as
  cyclers when they provide recurring access to the terrestrial and lunar vicinities without
  remaining permanently attached to either body."
- Early literature: Newton 1959 (periodic orbits passing close to two masses); Huang 1962 and
  Huang and Wade 1963 (two planar families of Earth-Moon periodic orbits enclosing both bodies);
  Arenstorf 1963 (existence proof); Kevorkian and Lancaster 1968 (matched asymptotics, "prescribed
  synodic commensurability"); Hoelker and Winston 1968. "This early literature already contains
  the essential idea of an Earth-Moon cycler."
- Free return: Schwaniger 1963 classified symmetric free returns into circumlunar (far-side
  periselenum) and cislunar (near-side) cases, "noting that one such free-return case closes as
  a periodic trajectory". The Apollo figure-eight is "not as a cycler in the repeated-access
  sense, but as the one-pass mission analogue". The only "gravity assist" in the body is
  "the lunar encounter itself supplies the passive gravity assist onto a transearth return".
  Genova and Aldrin 2015 embed Apollo-like figure-eight segments in "monthly Earth-Moon
  cyclers, including Arenstorf-type four-leaf-clover or big-loop geometries".
- Modern two-class distinction: "high-energy, near-Keplerian cyclers" organised by "p:q
  resonant circumterrestrial orbits corrected for lunar flybys, together with their
  three-dimensional CR3BP generalizations" (cites Binder and Arnas 2024), and "lower-energy
  gateway-mediated cyclers" organised by "the homoclinic and manifold structure associated
  with the Earth-Moon collinear gateways" (Davidson 1964). Casoliva et al. 2010 is cited for
  the two-class split in the planar CR3BP (resonant high-energy versus lower-energy
  L1-associated homoclinic-type cyclers, the latter "an outgrowth of unstable homoclinic
  structure rather than ... a stable family"). **Ross and Roberts-Tsoukkas 2025**: "stable,
  low-energy prograde Earth-Moon cyclers can also exist in this gateway-mediated class".
  Liang et al. 2020 (p:q resonant cycler orbits for a cislunar infrastructure) "repeatable
  access patterns, regular phasing opportunities, and distributed coverage of geolunar space".
- Collision lineage: Henon 1968, Hitzl 1977, Prado and Broucke 1994 and 1996; "Periodic
  collision orbits may therefore be read as the singular skeleton of repeated flyby families";
  "cyclers, near-collision families, and resonant transport are different faces of the same
  global Earth-Moon phase-space architecture".
- Cartography caveat: a stable cycler family "would appear only indirectly" as coherent
  low-MEGNO pockets, "but such signatures are not by themselves cycler-family
  identifications; rather, they mark regions where periodic-orbit continuation or targeted
  recurrence analysis would be dynamically well motivated."
- Braik and Ross 2026 ("Orbital networks in the three-body problem") is cited as motivation for
  repeated-encounter orbits as candidates for communications, navigation, logistics and
  staging architectures.
- The Primer gives no cycler orbit, period, V-infinity, encounter sequence or Jacobi constant.
  "Lunar gravity assist" as a topic is represented only by Negri et al. 2019 and Ross and
  Scheeres 2007 in the references.

## Which of the group's earlier papers it summarises, and what the corpus index holds
(Checked read-only against `docs/notes/CORPUS_INDEX.md` by file-name and title search.)

| Cited group paper (as in Primer references) | In CORPUS_INDEX |
|---|---|
| Rawat, Kumar, Rosengren, Ross (2026) "Cislunar mean-motion resonances: definitions, widths, and comparisons with resonant satellites," JGCD 49:4 | held (rawat-...-2026-...-JGCD-49-4; digest 2026-06-30) |
| Rawat, Kumar, Rosengren, Ross (2025) "Regions of influence of exterior mean-motion resonances in the Earth-Moon system," AAS 25-569 | held (digest 2026-07-15) |
| Kumar, Rawat, Rosengren, Ross (2026) "Cislunar resonant transport and heteroclinic pathways: From 3:1 to 2:2 to L1," Adv. Space Res. 77:3815-3843 | **not held under this title**; the index holds a related 2024 paper by the same four authors, "Interior 4:1/3:1/2:1 MMR bifurcations and heteroclinic connections" (IAC-24-C1.9.5); whether it is a precursor of the 2026 article is not stated in the Primer |
| Kumar and Anderson (2026) "Resonances, invariant manifolds, and low-thrust lunar transfers: the case of ESA's SMART-1," 30th ISSFD | **not held** |
| Ross and Roberts-Tsoukkas (2025) "Stable, low-energy prograde Earth-Moon cycler orbits," AAS 25-621 | held (ross-roberts-tsoukkas-2025-...-AAS-25-621) |
| Braik and Ross (2026) "Orbital networks in the three-body problem," arXiv:2605.31543 | held |
| Rosengren et al. 2019, 2020 (dynamical cartography of Earth-satellite orbits; lunar multipoles on elongated orbits), Amato et al. 2020, and the 2014 Laplace-plane papers | **not found** in the index |

Held in the corpus but **not cited** in this Primer's reference list (searched): Rosengren,
Ross, Kumar, Rawat 2024 (xGEO resonant structure, AMOS) and Rawat et al. 2024 (AAS 24-368).
Other cited works the index does hold: Casoliva et al. 2010 (and the 2008 AIAA precursor),
Genova and Aldrin 2015 (index title differs: "A Free-Return Earth-Moon Cycler Orbit"), Ross
and Scheeres 2007, Broucke 1968, Doedel et al. 2007, Koon-Lo-Marsden-Ross 2000 (heteroclinic
connections and resonance transitions). Not found in the index by name: Belbruno-Topputo-Gidea
2008, Topputo-Belbruno-Gidea 2008, Lei and Xu 2018, Liang et al. 2017 and 2020, Binder and
Arnas 2024, Pan et al. 2026, Peng et al. 2024, Li et al. 2026, Henry and Scheeres 2023,
McCarthy and Howell 2023, Davidson 1964, Arenstorf 1963, Huang 1960-1969, Schwaniger 1963,
Hitzl 1977, Prado and Broucke 1994 and 1996, Bonasera and Bosanac 2023 (now on disk; digested
separately 2026-10-03), Bosanac 2026 (clustering).

## What it does NOT contain
- No cycler orbit, periodic-orbit family member, initial condition, period or V-infinity table;
  no encounter sequence or lunar-gravity-assist itinerary; no resonance-hopping tour design.
- No heteroclinic or homoclinic connection computed or tabulated (counts and Jacobi levels are
  deferred to Rawat et al. and Kumar et al.); no connection between quasi-periodic tori
  (tori appear only as MEGNO-map signatures and in lists of references).
- No Jupiter, Saturn, Uranus or other planetary-moon system; Earth-Moon and Sun-Earth only.
- No new periodic-orbit continuation; the numerical results are (a, e) MEGNO and fate maps
  in elliptic Earth-Moon and Earth-Moon-Sun point-mass models at a single epoch.
- No statement on novelty of any cycler family.

## Bibliography relevant to resonances, cyclers, connections (entries as printed in the Primer)
- Arenstorf RF (1963) Existence of periodic solutions passing near both masses of the restricted three-body problem. AIAA Journal 1:238-240
- Belbruno EA, Miller JK (1993) Sun-perturbed Earth-to-Moon transfers with ballistic capture. Journal of Guidance, Control, and Dynamics 16:770-775
- Belbruno EA, Topputo F, Gidea M (2008) Resonance transitions associated to weak capture in the restricted three-body problem. Advances in Space Research 42:1330-1351
- Binder D, Arnas D (2024) Reliable and repeatable transit through cislunar space using 2:1 resonant spatial orbits. Journal of Guidance, Control, and Dynamics 47:1973-1979
- Bonasera S, Bosanac N (2023) Computing natural transitions between tori near resonances in the Earth-Moon system. Journal of Guidance, Control, and Dynamics 46:443-454
- Bosanac N (2026) Clustering natural trajectories in the Earth-Moon circular restricted three-body problem. The Journal of the Astronautical Sciences 73:2 (41 pp.)
- Braik A, Ross SD (2026) Orbital networks in the three-body problem, arXiv:2605.31543
- Broucke RA (1968) Periodic Orbits in the Restricted Three-Body Problem with Earth-Moon Masses. JPL Technical Report 32-1168, NASA
- Broucke RA (1969) Periodic Orbits in the Elliptic Restricted Three-Body Problem. JPL Technical Report 32-1360, NASA
- Casoliva J, Mondelo JM, Villac BF, et al (2010) Two classes of cycler trajectories in the Earth-Moon system. Journal of Guidance, Control, and Dynamics 33:1623-1640
- Davidson MC (1964) Numerical examples of transition orbits in the restricted three body problem. Astronautica Acta 10:309-313
- Doedel EJ, Romanov VA, Paffenroth RC, et al (2007) Elemental periodic orbits associated with the libration points in the circular restricted 3-body problem. International Journal of Bifurcation and Chaos 17:2625-2677
- Dutt P, Sharma RK (2010) Analysis of periodic and quasi-periodic orbits in the Earth-Moon system. Journal of Guidance, Control, and Dynamics 33:1010-1017
- Ferrari F, Lavagna M (2018) Periodic motion around libration points in the elliptic restricted three-body problem. Nonlinear Dynamics 93:453-462
- Gallardo T (2006) Atlas of the mean motion resonances in the Solar System. Icarus 184:29-38; (2019) Strength, stability and three dimensional structure of mean-motion resonances in the Solar System. Icarus 317:121-134; (2020) Three-dimensional structure of mean-motion resonances beyond Neptune. Celestial Mechanics and Dynamical Astronomy 132:9
- Gallardo T, Beauge C, Giuppone CA (2021) Semianalytical model for planetary resonances: Application to planets around single and binary stars. Astronomy and Astrophysics 646:A148
- Genova AL, Aldrin B (2015) Circumlunar free-return cycler orbits for a manned Earth-Moon space station. AAS/AIAA Astrodynamics Specialist Conference, Vail, CO, Paper AAS 15-794
- Henry DB, Scheeres DJ (2023) Quasi-periodic orbit transfer design via whisker intersection sets. Journal of Guidance, Control, and Dynamics 46:1929-1944
- Henon M (1968) Sur les orbites interplanetaires qui rencontrent deux fois la terre. Bulletin astronomique 3:377-402
- Hitzl DL (1977) Generating orbits for stable close encounter periodic solutions of the restricted problem. AIAA Journal 15:1410-1418
- Hoelker RF, Winston BP (1968) A Comparison of a Class of Earth-Moon Orbits with a Class of Rotating Kepler Orbits. NASA TN D-4903
- Howell KC (2001) Families of orbits in the vicinity of the collinear libration points. The Journal of the Astronautical Sciences 49:107-125
- Huang SS (1962) Preliminary study of orbits of interest for Moon probes. The Astronomical Journal 67:304-310; Huang SS, Wade C Jr (1963) ... II. The Astronomical Journal 68:388-391
- Kevorkian J, Lancaster JE (1968) An asymptotic solution for a class of periodic orbits of the restricted three-body problem. The Astronomical Journal 73:791-806
- Koon WS, Lo MW, Marsden JE, Ross SD (2000) Heteroclinic connections between periodic orbits and resonance transitions in celestial mechanics. Chaos 10:427-469
- Koon WS, Lo MW, Marsden JE, Ross SD (2022) Dynamical Systems, the Three-Body Problem and Space Mission Design. Marsden Books
- Kumar B, Anderson RL (2026) Resonances, invariant manifolds, and low-thrust lunar transfers: the case of ESA's SMART-1. 30th International Symposium on Space Flight Dynamics, Toulouse, France
- Kumar B, Rawat A, Rosengren AJ, Ross SD (2026) Cislunar resonant transport and heteroclinic pathways: From 3:1 to 2:2 to L1. Advances in Space Research 77:3815-3843
- Lei H, Xu B (2018) Resonance transition periodic orbits in the circular restricted three-body problem. Astrophysics and Space Science 363:70
- Li R, Masdemont JJ, Zhu Z, Gao C (2026) The Sun-Earth heteroclinics in restricted four-body nonautonomous models. Communications in Nonlinear Science and Numerical Simulation 152:109,356
- Liang Y, Xu M, Xu S (2017) The cislunar polygonal-like periodic orbit: Construction, transition and its application. Acta Astronautica 133:282-301
- Liang Y, Xu M, Peng K, Xu S (2020) A cislunar in-orbit infrastructure based on p:q resonant cycler orbits. Acta Astronautica 170:539-551
- McCarthy B, Howell K (2023) Construction of heteroclinic connections between quasi-periodic orbits in the three-body problem. The Journal of the Astronautical Sciences 70, 24
- Negri RB, Sukhanov A, de Almeida Prado AFB (2019) Lunar gravity assists using patched-conics approximation, three and four body problems. Advances in Space Research 64:42-63
- Newton RR (1959) Periodic orbits of a planetoid passing close to two gravitating masses. Smithsonian Contribution to Astrophysics 3:69-78
- Pan S, Urashi T, Bando M, Yoshimura Y, Chen H, Hanada T (2026) Data-driven prediction of chaotic transition in periapsis Poincare maps. Nonlinear Dynamics 114:575
- Prado AFBA, Broucke R (1994) Study of Henon's orbit transfer problem using the Lambert algorithm. JGCD 17:1075-1081; (1996) Transfer obits in the Earth-Moon system using a regularized model. JGCD 19:929-933
- Rawat A, Kumar B, Rosengren AJ, Ross SD (2025) Regions of influence of exterior mean-motion resonances in the Earth-Moon system: Bifurcations, separatrices, and heteroclinic pathways. AAS/AIAA Astrodynamics Specialist Conference, Boston, paper AAS 25-569
- Rawat A, Kumar B, Rosengren AJ, Ross SD (2026) Cislunar mean-motion resonances: Definitions, widths, and comparisons with resonant satellites. Journal of Guidance Control and Dynamics 49:4 (15 pp.)
- Ross SD, Roberts-Tsoukkas M (2025) Stable, low-energy prograde Earth-Moon cycler orbits. AAS/AIAA Astrodynamics Specialist Conference, Boston, paper AAS 25-621
- Ross SD, Scheeres DJ (2007) Multiple gravity assists, capture, and escape in the restricted three-body problem. SIAM Journal on Applied Dynamical Systems 6:576-596
- Schwaniger AJ (1963) Trajectories in the Earth-Moon Space with Symmetrical Free-Return Properties. NASA TN D-1833
- Topputo F, Belbruno E, Gidea M (2008) Resonant motion, ballistic escape, and their applications in astrodynamics. Advances in Space Research 42:1318-1329
- Vaquero M, Howell KC (2014a) Design of transfer trajectories between resonant orbits in the Earth-Moon restricted problem. Acta Astronautica 94:302-317; (2014b) Leveraging resonant-orbit manifolds to design transfers between libration-point orbits. JGCD 37:1143-1157
- Zimovan-Spreen EM, Howell KC, Davis DC (2020) Near rectilinear halo orbits and nearby higher-period dynamical structures: orbital stability and resonance properties. Celestial Mechanics and Dynamical Astronomy 132:28

(Selected, not complete: the Primer has several hundred references on perturbation theory,
catalogue objects and terminology which are omitted here. Diacritics were dropped in
transcription. The text cites "Peng et al 2024" in the heteroclinic passage; the reference list
has Peng C, Zhang Y, He S (2024) "3:1/3:2 resonant orbits touring L3-L5 in cislunar space,"
Advances in Space Research 73:2499-2514.)
