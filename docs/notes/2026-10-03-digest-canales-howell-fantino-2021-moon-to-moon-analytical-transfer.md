# Digest — Canales, Howell, Fantino (2021), "Transfer design between neighborhoods of planetary moons in the circular restricted three-body problem: The Moon-to-Moon Analytical Transfer Method" (arXiv v1)

**Digested:** 2026-10-03 (text-layer PDF, 47 pages, read in full from `pdftotext -layout`; Figure 39 checked on the
rendered page; the tables are text-layer tables). **Version:** the file is the arXiv v1 preprint. A later review
cites the published version as Celestial Mechanics and Dynamical Astronomy 133(8):36 (2021); the published
version was not checked and this paper's text does not itself print that citation. **Purpose:** record exactly
which moon pairs, orbits and transfers this paper treats, because the project's records say it covers
Ganymede-Europa and Titania-Oberon. Facts only; no novelty verdict is written here.

## Citation
David Canales, Kathleen C. Howell, Elena Fantino, "Transfer design between neighborhoods of planetary moons in the
circular restricted three-body problem: The Moon-to-Moon Analytical Transfer Method".

- Printed identifier: "arXiv:2110.03683v1 [astro-ph.EP] 7 Oct 2021". The front page also prints the Springer
  template placeholders "Noname manuscript No. (will be inserted by the editor)" and "Received: date / Accepted:
  date".
- Affiliations as printed: Canales and Howell, School of Aeronautics and Astronautics, Purdue University, West
  Lafayette, IN 47907; Fantino, Aerospace Engineering Department, Khalifa University of Science and Technology,
  Abu Dhabi, United Arab Emirates.
- Keywords as printed: Multi-body dynamical systems; Circular restricted three-body problem (CR3BP); Libration
  point orbits; Spacecraft trajectory design; Moon-to-Moon Transfers; Moon tour design.
- Filed in the private paper corpus as
  `canales-howell-fantino-2021-transfer-design-neighborhoods-planetary-moons-cr3bp-moon-to-moon-analytical-transfer-arxiv-2110.03683v1.pdf`,
  md5 `bfc826dc65dfc116759686c06c9ffac0`. (A different Canales-Howell 2023 arXiv paper is a separate file in the
  corpus and is not this one.)

## What the paper is
Abstract (quoted): "Given the interest in future space missions devoted to the exploration of key moons in the
solar system and that may involve libration point orbits, an efficient design strategy for transfers between
moons is introduced that leverages the dynamics in these multi-body systems. The moon-to-moon analytical transfer
(MMAT) method is introduced, comprised of a general methodology for transfer design between the vicinities of the
moons in any given system within the context of the circular restricted three-body problem, useful regardless of
the orbital planes in which the moons reside. A simplified model enables analytical constraints to efficiently
determine the feasibility of a transfer between two different moons moving in the vicinity of a common planet. In
particular, connections between the periodic orbits of such two different moons are achieved. The strategy is
applicable for any type of direct transfers that satisfy the analytical constraints. Case studies are presented
for the Jovian and Uranian systems. The transition of the transfers into higher-fidelity ephemeris models
confirms the validity of the MMAT method as a fast tool to provide possible transfer options between two
consecutive moons."

Sections: 1 Introduction; 2 Dynamical Models (CR3BP; 2BP-CR3BP patched model; coupled spatial CR3BP; higher-
fidelity ephemeris model); 3 Moon-to-moon transfers relying solely on the coupled spatial CR3BP; 4 The MMAT
Method (4.1 coplanar moons, Theorem 1; 4.2 non-coplanar moons, Theorem 2; 4.3 comparison; 4.4 Titania to Oberon;
4.5 epoch dependence along revolutions; 4.6 dependence on the sphere-of-influence definition); 5 Transition to a
higher-fidelity ephemeris model; 6 Discussion and Concluding Remarks; Appendices A-E. Tables 1-5, Figures 1-43.

## Systems and moon pairs treated (verified from the paper; complete list)
1. **Ganymede to Europa, Jovian system** (Sections 3, 4.1, 4.2, 4.3, 4.6, 5): L1 Lyapunov orbit in the
   Jupiter-Ganymede CR3BP (departure) to L2 Lyapunov orbit in the Jupiter-Europa CR3BP (arrival). Treated three
   ways: moons coplanar, moons in their true planes, and a demonstration with the moons placed in arbitrary
   planes (Omega_JG = 100 deg, i_JG = 60 deg, Omega_JE = 200 deg, "i_JG = 20 deg" as printed, evidently the Europa
   inclination; "uniquely presented to validate the method").
2. **Titania to Oberon, Uranian system** (Sections 4.4, 4.5, 4.6, 5): L2 northern halo orbit in the
   Uranus-Titania CR3BP (departure) to L1 southern halo orbit in the Uranus-Oberon CR3BP (arrival).

No other moon pair is a transfer endpoint. Direction is one way in both cases: departure moon (Ganymede, Titania)
to arrival moon (Europa, Oberon). The paper states that, being geometric, the constraints apply "regardless of
the journey being outward or inward". Table 1 gives orbital data for exactly four moons (Europa, Ganymede,
Titania, Oberon).

### Table 1, exact transcription
Caption as printed: "Orbital data for Europa, Ganymede, Titania and Oberon obtained from the SPICE database and
referred to the Ecliptic J2000.0 reference frame (Acton et al., 2017). Last accessed 08/05/2020."

| Moon | Semi-major axis [10^5 km] | Orbital period [day] | CR3BP mass ratio [10^-5] | Eccentricity [10^-3] | Inclination [degree] | Longitude ascending node [degree] |
|---|---|---|---|---|---|---|
| Europa | 6.713 | 3.554 | 2.528 | 9.170 | 2.150 | 331.361 |
| Ganymede | 10.706 | 7.158 | 7.804 | 2.542 | 2.208 | 340.274 |
| Titania | 4.363 | 8.708 | 3.917 | 1.871 | 97.829 | 167.627 |
| Oberon | 5.836 | 13.471 | 3.544 | 1.170 | 97.853 | 167.720 |

### Umbriel, Ariel, Miranda
Umbriel: 1 occurrence; Ariel: 1 occurrence; Miranda: 0. Both hits are one sentence in Section 5: "The higher-
fidelity ephemeris model for this example is defined with Uranus as the central body, and the Sun, Ariel,
Umbriel, Titania and Oberon as the perturbing bodies." They are perturbers only, not transfer endpoints. For the
Jovian ephemeris example the perturbers are "the Sun, Io, Europa, Ganymede, Callisto and Saturn".

## Models
- CR3BP (circular, planar or spatial) near each moon; both moons assumed on circular orbits about the planet;
  moon planes fixed by the Ecliptic J2000.0 elements at the epoch.
- "2BP-CR3BP patched model": CR3BP arcs near a moon are propagated to the edge of a moon "sphere of influence"
  and then replaced by an osculating Keplerian conic about the planet; the departure and arrival conics are
  intersected analytically. The sphere of influence is defined by the acceleration ratio a_SoI = a_m/a_p, with
  a_SoI = 5e-4 used (R_SoI approximately 3e5 km for Ganymede, 1.2e5 km Europa, 9.54e4 km Titania, 1.23e5 km Oberon,
  as printed).
- "Coupled spatial CR3BP": "two CR3BPs with a common primary ... It is not a dynamical four-body problem, but
  rather a kinematical blending." Differential corrections by multiple shooting (the tau-alpha method of Haapala
  and Howell, 2015).
- Ephemeris model: SPICE ephemeris, Ecliptic J2000.0 planet-centred; "three Delta-v's are allowed across the
  entire transfer at: (1) instant 0 to depart the departure periodic orbit; (2) instant 2 to transfer from the
  departure arc towards the arrival arc; (3) instant 4 to ensure that the s/c arrives on the specified periodic
  orbit."
- Theorem 1 (coplanar): condition 2 a_d a_a (1 + e_d e_a) >= b_d^2 + b_a^2 >= 2 a_d a_a (1 - e_d e_a) (Eq. 8)
  for a tangential intersection of two confocal conics; Theorem 2 (non-coplanar): a_a(1 - e_a) <= a_d(1 -
  e_d^2)/(1 + e_d cos(theta_d,Int + n pi)) <= a_a(1 + e_a), n = 0, 1 (Eq. 21). Rephasing of the arrival moon by
  Eq. 20 gives the required arrival-moon phase; one impulsive maneuver Delta-v_tot at the intersection.
- Manifold arcs: departure arcs are unstable manifolds, arrival arcs stable manifolds, of the L1/L2 periodic
  orbits; "The CR3BP arcs in this investigation are manifold trajectories, but it could also be applied to other
  types of arcs, such as transit orbits."

## Departure and arrival orbits (Jacobi constants as printed; periods are not printed)
| Transfer | Departure orbit | JC departure | Arrival orbit | JC arrival |
|---|---|---|---|---|
| Ganymede to Europa | L1 Lyapunov, Jupiter-Ganymede | JC_d = 3.0061 | L2 Lyapunov, Jupiter-Europa | JC_a = 3.0024 (Sections 3, 4; Figs. 18-38) |
| Titania to Oberon | L2 northern halo, Uranus-Titania | JC_d = 3.0035 | L1 southern halo, Uranus-Oberon | JC_a = 3.003 |

Inconsistency inside the paper, as printed: the captions of Figures 7, 9 and 11 give the Jupiter-Europa L2
Lyapunov orbit as "JC_a = 3.0028", whereas the text of Section 3 and every later figure caption give "JC_a =
3.0024". The paper does not comment. The orbital periods of the Lyapunov and halo orbits are not stated, and the
halo orbits' amplitude and stability index are not stated; only Jacobi constants identify them.

## Transfer results, exact transcriptions

### Table 2 ("Comparison of resulting transfers assuming moons reside in coplanar orbits and in their true orbital planes in the coupled CR3BP")
| Case | Delta-v_tot [km/s] | t_tot [days] |
|---|---|---|
| Coplanar moon orbits | 0.9456 | 9.470 |
| True moon orbits: Minimum-Delta-v_tot transfer | 0.9422 | 9.473 |
| True moon orbits: Minimum-t_tot transfer | 0.9428 | 9.472 |

### Table 3 ("Summary of the resulting transfers between halo orbits of the Titania and Oberon vicinities depending on Titania's a_SoI (coupled spatial CR3BP)")
| Case | Delta-v_tot [m/s] | t_tot [days] |
|---|---|---|
| a_SoI = 5e-4 for Titania and Oberon | 45.7 | 28.44 |
| Titania a_SoI = 5e-3; Oberon a_SoI = 5e-4 | 66.7 | 17.72 |

### Table 4 ("Comparison between coupled spatial CR3BP and higher-fidelity ephemeris model for selected transfers")
| Ganymede to Europa transfer | CR3BP Delta-v_tot [km/s] | CR3BP t_tot [days] | Ephemeris Delta-v_tot [km/s] | Ephemeris t_tot [days] |
|---|---|---|---|---|
| Minimum-Delta-v_tot transfer | 0.9422 | 9.473 | 0.966 | 9.47 |
| Maximum-t_tot transfer | 1.08 | 12.13 | 1.02 | 12.13 |

| Titania to Oberon transfer | CR3BP [m/s] | CR3BP [days] | Ephemeris [m/s] | Ephemeris [days] |
|---|---|---|---|---|
| Minimum-Delta-v_tot transfer | 45.7 | 28.44 | 42.1 | 28.44 |

### Table 5 ("Summary of the sample applications for transfers between Lyapunov orbits of Ganymede and Europa and between halo orbits of Titania and Oberon")
| Transfer between Lyapunov orbits of Ganymede and Europa | CR3BP Delta-v_tot [km/s] | CR3BP t_tot [days] | Ephemeris Delta-v_tot [km/s] | Ephemeris t_tot [days] |
|---|---|---|---|---|
| Minimum-Delta-v_tot transfer | 0.9422 | 9.473 | 0.966 | 9.47 |
| Maximum-t_tot transfer | 1.08 | 12.13 | 1.02 | 12.13 |
| Minimum-t_tot transfer | 0.9428 | 9.472 | - | - |

| Transfer between halo orbits of Titania and Oberon | CR3BP [m/s] | CR3BP [days] | Ephemeris [m/s] | Ephemeris [days] |
|---|---|---|---|---|
| Minimum-Delta-v_tot transfer | 45.7 | 28.44 | 42.1 | 28.44 |
| Minimum-Delta-v_tot transfer with longer t_tot | 50 | 39.6 | - | - |
| Minimum-Delta-v_tot transfer reducing Titania's SoI | 66.7 | 17.72 | - | - |

### Figure-caption values (not in the tables)
- Fig. 13 (Section 3.1, coplanar, converged by Poincare-section and nearest-neighbour search): Delta-v_tot =
  1.0192 km/s, t_tot = 28.31 days, Ganymede L1 Lyapunov to Europa L2 Lyapunov, "coupled planar CR3BP".
- Fig. 19(b) (coplanar patched model): 0.9433 km/s, 9.47 days; Fig. 20 (coupled planar CR3BP): 0.9456 km/s, 9.47
  days. Fig. 24 (patched model, true planes): (a) minimum Delta-v 0.9448 km/s, 9.473 days; (b) minimum t_tot
  0.9455 km/s, 9.471 days. Fig. 25 (coupled spatial CR3BP): (a) 0.9422 km/s, 9.473 days; (b) 0.9428 km/s, 9.472
  days. Text: "transfers between the libration points of both Ganymede and Europa in less than 10 days using a
  single Delta-v_tot equal to 942 m/s"; "about 5.6 days is required for the s/c to leave the Ganymede SoI and 2.9
  days to arrive to the L2 Lyapunov orbit in the Europa vicinity from the Europa SoI".
- Fig. 26 (arbitrary-plane demonstration, not the real moons): Delta-v_tot = 12.98 km/s, t_tot = 11.83 days.
- Fig. 31 (Titania to Oberon, patched model): 62.6 m/s, 28.5 days. Fig. 32 (coupled spatial CR3BP): 45.7 m/s,
  28.44 days; text "about 7 days for the s/c to leave the Titania SoI, and 10.8 days measured from the Oberon SoI
  to arrive into the L1 southern halo orbit". Fig. 33 (an extra revolution of the arrival conic): 50 m/s, 39.6
  days. Fig. 36: 45.7 m/s, 28.44 days and 66.7 m/s, 17.72 days. Fig. 40 (ephemeris): 42.1 m/s, 28.44 days.
- Fig. 37 (ephemeris, Ganymede to Europa): 0.966 km/s, 9.47 days at epoch JD2459323.365 ("April 18th, 2021").
  Fig. 38: (a) CR3BP 1.08 km/s, 12.13 days; (b) ephemeris 1.02 km/s, 12.13 days at epoch JD2459138.765 ("April
  16th, 2021" as printed; reader's arithmetic: that JD is about 185 days before JD2459323.365, which is mid-
  October 2020, so the printed date does not match the printed Julian date). Fig. 40 epoch JD2459766.565 ("July
  6th, 2022").
- Text summary: "transfers between halo orbits near both Titania and Oberon are identified that occur in less
  than a month and using a total single Delta-v_tot of less than 50 m/s"; "transfers between halo orbits near
  both Titania and Oberon that require about 2 weeks with a Delta-v_tot as little as 66.7 m/s".

### Phase conditions
- Section 3 demonstration only: Ganymede at an orbital phase of 45 deg relative to Europa at t0; Poincare
  sections at 90 deg (no intersection) and 270 deg.
- The method's phase variable is the departure angle from the departure moon, theta_0^Gan (Ganymede) or
  theta_0^Tit (Titania), measured from the moon's ascending-node direction in its orbit; for each value, Theorem 2
  either fails (no direct transfer) or gives a unique arrival-moon phase theta_4 and a unique single Delta-v_tot
  and t_tot. "For a transfer between two moons involving a single maneuver, the cost and time are entirely
  dependent upon the initial epoch, represented in this example as theta_0^Gan, with some values of the angle
  providing no access whatsoever." Coplanar case: the Delta-v_tot and t_tot "remain the same regardless of the
  angle of departure". Numerical values of the optimal phases are not given in the text; they appear in Figures
  23, 39 and 42, which are plots. Read by eye from the left panel of Figure 39 (approximate): four
  minimum-Delta-v departure angles near 35, 130, 215 and 310 deg with Delta-v about 0.94 km/s, and
  infeasible (yellow) bands near 0-15, 155-190 and 335-350 deg.
- Epoch dependence: an extra revolution on the arrival conic changes the relative Titania-Oberon phase at nearly
  the same Delta-v at about 10 days more time of flight (Fig. 33, Table 5 row "longer t_tot"); the sphere-of-
  influence ratio a_SoI is varied from 5e-5 to 1e-2 for Ganymede (Fig. 35).

## One-way or repeating?
All transfers are single one-way direct transfers with one impulsive Delta-v_tot (three in the ephemeris model:
leave the departure orbit, switch arcs, enter the arrival orbit). No return leg, no second transfer and no
repeating pattern is constructed. The nearest repeating-type content is Figure 39, "Event recurrence of the
minimum-Delta-v_tot configurations over the span of 5 years": the text says "A representation of minimum-
Delta-v_tot configurations that recur over the span of 5 years appears in Fig. 39." This plots the dates over
which the Ganymede-Europa phasing that allows the minimum-cost one-way transfer recurs (access / no-access vs
Julian date, 2.459e6 to 2.461e6); it is a set of departure windows for the same one-way transfer, not a
trajectory that returns to the moons. The introduction mentions in general terms that tour design is hard for
"leveraging libration point orbits, captures or even returning to a previously visited moon, since the gravity
assist sequence is affected"; that is a statement about other methods, not a result.

## Term search (grep -a, lines joined, ligatures normalised)
cycler 0; cycling 0; periodic 36 (periodic orbits, periodic orbit families, "planar and spatial periodic orbits",
reference titles); repeating 0; free-return 0; resonant 2 ("combine the resonant gravity assist technique with
manifolds that originate from Lyapunov periodic orbits", attributed to Lantoine et al. 2011, and the title of that
reference "resonant hopping transfers"); resonance 1 (title of Lynam et al., "Laplace resonance"); heteroclinic
1; homoclinic 0; quasi-periodic 2 (one body sentence, one reference title); torus 0; four-body 1; Umbriel 1;
Titania 39; Oberon 40; Ariel 1; Miranda 0; Triton 0; Proteus 0; Nereid 0. (Counts include captions, headers and
references.)

Sentences for the requested terms:
- cycler, cycling, repeating, free-return: none, those words do not occur.
- heteroclinic: "Note that if a crossing also occurs at the y-xdot section, a heteroclinic connection is produced
  as explained in Haapala and Howell (2015)." (Section 3.1, the Ganymede-Europa coplanar demonstration; the
  converged solution there has a Delta-v, so it is not a heteroclinic connection.)
- quasi-periodic: "The motion is frequently categorized by different types of families of periodic and quasi-
  periodic orbits according to their geometry and stability" (Section 2.1, background; no quasi-periodic orbit is
  computed).
- four-body: "It is not a dynamical four-body problem, but rather a kinematical blending." (Section 2.3).
- resonant: the Lantoine et al. sentence above (introduction, background).

## Does the paper contain a trajectory that returns to the same two moons, or a periodic or quasi-periodic orbit that encounters two moons
No. Each solution is a one-way transfer between a periodic orbit near one moon and a periodic orbit near another
moon (Lyapunov or halo orbits around L1 or L2 of one planet-moon CR3BP each); the periodic orbits are each
confined to one moon's vicinity, and no periodic or quasi-periodic orbit is constructed that passes both moons.
The connecting arc is a manifold-to-manifold transfer with an impulsive maneuver.

## What it does NOT contain
- Any transfer back from the arrival moon, or any sequence Moon A to Moon B to Moon A.
- Any moon pair other than Ganymede-Europa and Titania-Oberon as endpoints; Umbriel and Ariel only as ephemeris
  perturbers; Miranda absent; no Saturnian or Neptunian system.
- Any quasi-periodic or periodic orbit that encounters two moons, any resonant (m:n) cycler or tour, any
  gravity-assist sequence, any V-infinity leveraging, any Tisserand-graph result of its own.
- Orbit periods for the Lyapunov and halo orbits, stability indices, or the halo amplitude; the transfers are
  identified only by Jacobi constant.
- Phase angles for the optimum configurations stated in the text (they appear only as plots).
- Statements of an end-to-end moon tour; the paper says "Additional insights are then possible when designing
  tours", and the single-transfer result is offered as a building block.

## Relevant bibliography (works on in-system transfers and moon-tour methods, as printed)
Campagnola S, Skerritt P, Russell RP (2012) Flybys in the planar, circular, restricted, three-body problem.
Celestial Mechanics and Dynamical Astronomy 113(3):343-368, DOI 10.1007/s10569-012-9427-x
Colasurdo G, Zavoli A, Longo A, Casalino L, Simeoni F (2014) Tour of jupiter galilean moons: Winning solution of
gtoc6. Acta Astronautica 102:190-199, DOI 10.1016/j.actaastro.2014.06.003
Diehl R, Kaplan D, Penzo P (1983) Satellite tour design for the Galileo mission, AIAA, 21st Aerospace Sciences
Meeting. DOI 10.2514/6.1983-101
Fantino E, Castelli R (2017) Efficient design of direct low-energy transfers in multi-moon systems. Celestial
Mechanics and Dynamical Astronomy 127, DOI 10.1007/s10569-016-9733-9
Fantino E, Flores RM, Al-Khateeb AN (2018) Efficient two-body approximations of impulsive transfers between halo
orbits. 69th International Astronautical Congress (IAC), Bremen, Germany
Fantino E, Salazar F, Alessi EM (2020) Design and performance of low-energy orbits for the exploration of
enceladus. Communications in Nonlinear Science and Numerical Simulation 90:105393, DOI
10.1016/j.cnsns.2020.105393
Grover P, Ross S (2009) Designing trajectories in a planet-moon environment using the controlled keplerian map.
Journal of Guidance, Control, and Dynamics 32:437-444, DOI 10.2514/1.38320
Haapala A, Howell K (2015) A framework for constructing transfers linking periodic libration point orbits in the
spatial circular restricted three-body problem. International Journal of Bifurcation and Chaos 26:1630013, DOI
10.1142/S0218127416300135
Kakoi M, Howell K, Folta D (2014) Access to mars from earth-moon libration point orbits: Manifold and direct
options. Acta Astronautica 102, DOI 10.1016/j.actaastro.2014.06.010
Koon W, Lo M, Marsden J, Ross S (2001) Constructing a Low Energy Transfer Between Jovian Moons, vol 292, pp
129-145. DOI 10.1090/conm/292/04919
Lantoine G, Russell RP, Campagnola S (2011) Optimization of low-energy resonant hopping transfers between
planetary moons. Acta Astronautica 68(7):1361-1378, DOI 10.1016/j.actaastro.2010.09.021
Lynam AE, Kloster KW, Longuski JM (2011) Multiple-satellite-aided capture trajectories at jupiter using the
laplace resonance. Celestial Mechanics and Dynamical Astronomy 109(1):59-84, DOI 10.1007/s10569-010-9307-1
Short CR, Blazevski D, Howell KC, Haller G (2015) Stretching in phase space and applications in general
nonautonomous multi-body problems. Celestial Mechanics and Dynamical Astronomy 122(3):213-238, DOI
10.1007/s10569-015-9617-4
Strange N, Goodson T, Hahn Y (2002) Cassini Tour Redesign for the Huygens Mission, AIAA/AAS Astrodynamics
Specialist Conference and Exhibit. DOI 10.2514/6.2002-4720
Topputo F, Vasile M, Bernelli-Zazzera F (2005) Low energy interplanetary transfers exploiting invariant manifolds
of the restricted three-body problem. Journal of the Astronautical Sciences 53:353-372
Viale A (2016) Low-energy tour of the galilean moons. Master's thesis, Universita degli Studi di Padova
Also cited as mission references: Plaut et al. (2014) JUICE; Phillips and Pappalardo (2014) Europa Clipper;
Lorenz et al. (2018) Dragonfly. The reference list has no Uranus-orbiter mission study and no Neptune or Triton
work.

## Suggested topology label
**halo** (libration-point orbit transfers). Justification from the text: the abstract states "connections between
the periodic orbits of such two different moons are achieved" and the Uranian case is "Transfer from an L2
northern halo orbit of the Uranus-Titania system, to an L1 southern halo orbit of the Uranus-Oberon system".

## Provenance of this digest
Full text read from `pdftotext -layout` (47 pages). Tables 1-5 were transcribed from the text layer and compared
to the printed values in the text. Not checked: the published CMDA version, the arXiv later versions, and the
contents of cited works.
