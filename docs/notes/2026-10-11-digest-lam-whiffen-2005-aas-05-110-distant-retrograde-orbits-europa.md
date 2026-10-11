# Digest: Lam and Whiffen 2005, "Exploration of Distant Retrograde Orbits Around Europa", AAS 05-110 (#960 wanted-list row 25, formerly 26)

Try Lam and Gregory J. Whiffen (Outer Planet / Small Body Flight Dynamics Section, JPL), AAS 05-110, AAS/AIAA Space Flight Mechanics Meeting, Copper Mountain, Colorado, January 2005.
Author order: the first page lists Lam first (Lam*, Whiffen†). The wanted list gave the NTRS order "Whiffen & Lam" (NTRS 20210001733); the paper's own order is Lam and Whiffen, and
that is the order to cite. The wanted-list row has been corrected (see the wanted list change note dated 2026-10-11).
Filed as `cyclers_pdf/papers/lam-whiffen-2005-exploration-distant-retrograde-orbits-europa-aas-05-110.pdf`
(md5 `8dfae1b0aee142e411e6073f0b0b8410`, 19 pp., text layer present, a macOS Preview re-save of the paper; no DOI; original upload not modified).
No copy existed before: CORPUS_INDEX.md and `docs/notes/*.md` were grepped for "Whiffen", "05-110" and the title on 2026-10-11 (hits are other Whiffen citations only).

Status: READ in full (text layer, every equation, every caption, the reference list). Page renders viewed for pp.5, 10, 12, 13, 14 and 18 (Figs. 5, 11, 13 second, 14, 15, 16, 22).
The other figures are described from captions and running text and are marked CAPTION. The paper has NO numbered tables. All printed numbers are in
`data/sources/lam-whiffen-2005-aas-05-110-numbers.yaml` (text layer plus page render; the Fig. 22 annotation block is image text with one witness). COMPUTED = our arithmetic
2026-10-11; INFERRED = our reading. Page numbers are the paper's (title page = 1; PDF page = printed page).
Related digests, not repeated: `2026-10-06-digest-henon-1969-hill-case-periodic-orbits-stability.md` (families a, c, f, g, g', the source of the family nomenclature),
`2026-10-11-digest-henon-1974-vertical-stability-II-hill-case.md`, `2026-10-11-digest-benest-1976-libration-retrograde-satellites-celest-mech.md`,
`2026-10-08-digest-russell-2005-global-search-periodic-orbits-near-europa-aas-05-290.md` (same year, same lab, periodic orbits near Europa).

## 0. Verdict

A JPL mission-design paper that argues distant retrograde orbits (DROs, Henon's family f) about Europa are good end-of-mission quarantine orbits and good stable transfer paths
between low Europa orbits and distant ones. It is a Europa-centred, single-moon study with Jupiter as the perturber. There is no cycler, no multi-moon sequence, no flyby of
another moon and no resonance with a second moon. The only "resonance" language is the labelling of instability necks of family f by resonance ratios (5:2, 3:2, 4:3) with Europa's
orbital period and the period-three family g3 (sec. 3). See "Relevance to this project".

## 1. Abstract and introduction (pp.1-3)

- Claim (abstract): DROs "can be stable for a very long time"; they exist because of the larger perturbation from Jupiter; they suit proposed missions such as JIMO (Jupiter Icy Moons Orbiter);
  "ideal quarantine orbits" owing to long-term stability and existence over a large range of distances; continuous families of DROs give "instantaneously stable transfer paths" for escape and capture.
- Motivation (p.1): at low altitude Europa's field causes a secular eccentricity drift and surface impact; at high altitude Jupiter's third-body perturbation causes collision or escape within a few
  revolutions. DROs have stability "hundreds of years" (p.1). Compared with impacting Jupiter or escaping the Jovian system as end-of-mission options.
- Description (pp.1-2): in the inertial frame a DRO is retrograde with small loops that alternate between retrograde and direct motion, so classical elements do not apply. In the rotating frame (axis along the
  Jupiter-Europa line) it is a retrograde "imperfect ellipse" with its semi-minor axis along the Jupiter-Europa line, about 2:1 axis ratio for large orbits, near circular at about half the distance to L1 (p.2-3).
  Fig. 1 CAPTION: a 35,000 km DRO, inertial (a) and L1-centred rotating (b). Fig. 2: the same in the rotating frame.
- Stability background (p.3): Henon and Benest showed family f (simple-period planar DROs) is stable in Hill's approximation except at two members where the period-three family g3 crosses f twice
  ("instability necks", near the 3:1 resonance in Hill's case). Benest later showed stability regions exist for all mass ratios but pockets of instability appear (refs 3, 5, 6; Broucke, ref 4). Fig. 3 CAPTION: a period-three orbit of family g3.
- Tools (p.3): nearly all analysis done numerically with Mystic (JPL's low-thrust optimiser, Static/Dynamic Optimal Control, ref 7), with Hill's approximation for comparison.

## 2. Equation of motion and orbit characterisation (pp.3-8)

- Eq. (1): mu = M1 / (M1 + M2), M1 Europa, M2 Jupiter; Hill's approximation takes mu to 0. Printed mu(Jupiter-Europa) = 0.0000252803. COMPUTED: 2.52799e-5 from GM values (agrees).
- Eq. (2): Hill's equations centred on Europa: xi'' = 2 eta' + 3 xi - xi/r^3, eta'' = -2 xi' - eta/r^3, xi along Europa to Jupiter.
- Eq. (3): Jacobi constant Gamma = 3 xi^2 + 2/r - xi'^2 - eta'^2 (the text layer drops a glyph; standard Hill form). Larger orbits have Gamma toward minus infinity (higher energy).
- Plots without units are Hill-model plots; plots with units are Mystic (ephemeris) plots. Ephemeris orbits "can never be exactly periodic"; near-periodic ones are called periodic (p.4).
- Fig. 4 CAPTION (near-periodic DROs, 20 d, Mystic, bodies Sun, Jupiter, Io, Europa, Ganymede, Callisto): (a) 10,000 to 50,000 km in 5,000 km steps, (b) 50,000 to 150,000 km in 25,000 km steps.
  Large DROs distort and the tips bend toward Jupiter (p.5).
- Definitions (p.5): "simple-periodic" = symmetric about the rotating X-Z plane, crossing the rotating X axis twice per period, with no period-two or larger orbits within the family; N-periodic orbits bifurcate from
  family f with increasing energy. Henon's nine natural families, five simple-periodic symmetric: a and c (L2 and L1 Lyapunov), f (retrograde about M2), g (direct about M2), g' (branches from g). Hill-case
  orbits become circular (first kind), elliptic (second kind) or consecutive-collision (second species) in the limit (ref 10); the close approaches to M2 serve as starting states for transfers.
- Fig. 5 (viewed): rotating-frame period versus DRO X size, 0 to 100,000 km; the period rises from about 0.5 d at small size toward Europa's 3.55 d. No numbers printed.
- Eqs. (4) to (6): the asymptotic family f solution xi = K1 cos t, eta = -2 K1 sin t, Gamma = -K1^2, an exact ellipse centred on M2 with 2:1 axes and period 2 pi (Europa's period about Jupiter),
  retrograde; the abscissa and initial velocity for given Gamma: xi0 = -(-Gamma)^(1/2), V_eta0 = -2 xi0. This is the relative motion of two Keplerian orbits about Jupiter with equal periods (p.6).
- Family a (L2 Lyapunov) asymptotically collides with Europa as Gamma goes to minus infinity (Fig. 6 CAPTION, p.6); families a and c are prograde near M2 and retrograde far from it, which makes them transfer references.
- Non-periodic orbits (pp.6-8): in the ephemeris or the perturbed CR3BP, DROs are quasi-periodic (a librating ellipse oscillating along Y). Eq. (7), xi'^2 <= 3 xi^2 + 2/xi - Gamma, bounds the Hill Poincare section at Y = 0.
  Fig. 7 CAPTION: the Gamma = 1 section, a triangular stable region around a simple-periodic retrograde point, period-3 g3 at the three corners, period-4 as four bands; all orbits in the triangle are stable and
  only retrograde orbits remain at this energy. Fig. 8 CAPTION: ephemeris Poincare data for 35,000 km DROs with varied velocity (X = 35,000 km corresponds to x = 1.76579 Hill units; the outermost curve is
  X = 35,000 km, Vy = -0.816 km/s). Fig. 9: the same in Hill's model, with a longer "tail" because orbits close to Europa persist in Hill's case but escape in the ephemeris.
  COMPUTED: 35,000 km / 1.76579 implies a Hill length unit of 19,821 km; our constants give 19,696 km (x = 1.7770). The printed value is 0.6 percent lower; the paper does not state its constants.

## 3. Stability of DROs (pp.8-12)

- Method (p.8-9): because the ephemeris has no constant Jacobi integral, stability is shown with Benest's method (ref 6): the Y velocity at each positive X-axis crossing of the rotating X-Z plane,
  against X. Fig. 10 (CAPTION) and Fig. 11 (viewed): the dotted curve is the near-periodic DRO family; the two thick curves are the minimum and maximum Y velocities, at each X, that do not escape or collide within
  200 days ("stable" = does not escape or hit Europa for at least 200 d); velocities are inertial.
- Reading (p.9): except at three points on the curve the DROs are stable for at least 200 d. Two necks, where the bounds merge, are put "near the 5:2 and 4:3 resonances" and attributed to the unstable
  triple-periodic family g3 crossing f; the resonance where the necks occur differs between ephemeris and Hill cases. A third neck near 3:2 "seems to appear only for Europa". The instabilities appear for all Galilean moons.
  Each neck is a region where DROs should not be used as quarantine orbits.
  INFERRED inconsistency: the Fig. 11 render marks vertical lines at 5:2 (about 13,000 km), 3:2 (about 21,000 km) and 4:3 (about 24,500 km), and the text of p.9 names 5:2 and 4:3 as the g3 necks and
  3:2 as the Europa-only neck. The intro (p.3) says the Hill-case necks are near 3:1. The paper does not define the ratio (which period over which). Not resolved here. The neck near 13,000 to 13,500 km
  agrees with the L1/L2 distance (13,656 km COMPUTED), which Fig. 11 also marks near 13,800 km.
- The DRO Vy requirement exceeds that of a two-body circular orbit and rises linearly beyond about twice the L1 distance (p.9); the circular velocity falls with distance. Fig. 10-11 give initial conditions (Y = 0, Vx = 0)
  for any DRO size after rotating the X axis onto the Jupiter-Europa line (p.10).
- Out-of-plane DROs (pp.10-12): Henon showed vertical stability (refs 13, 14). Adding a Z velocity gives Lissajous-like orbits about L1 or L2. Fig. 12 and the first Fig. 13 (p.11): the 35,000 km DRO, Vy = -0.767, Vz = -0.55 km/s, propagated
  200 d. The second Fig. 13 (p.12, viewed): stable region in the (Vz, Vy) plane for 11,000 km (a wide U-shape between Vy about -0.50 and -0.60 km/s) and 13,200 km (a narrow strip near the neck). The paper
  has two figures numbered 13.
- Fig. 14 "Red Sea Plot" (viewed, p.12): orbital lifetime (colour scale 0 to 100 d) versus DRO X size (0 to 70,000 km, periapsis axis) and inclination (110 to 180 deg), by rotating the velocity of planar DROs out of the plane.
  A sharp "cliff" separates the stable blue region from the unstable red; the cliff persists under 1,000-year propagations. A stable region exists from very large to small radii at inclinations from 180 degrees
  down to about 135 degrees depending on distance, which the authors read as a continuous family of three-dimensional DROs in the ephemeris and "stable pathways accessible to low-thrust vehicles". The plot is conservative
  because the velocity magnitude was not varied. A thrust loss of "as little as two days" can mean impact or escape (p.12).

## 4. Transfers to and from DROs (pp.13-18)

- Impulsive transfers (pp.13-14): Ocampo and Rosborough (refs 11, 12) "Earth-Return periodic Orbits" (ERO) classes carried to Europa. Class A uses a family-a Lyapunov orbit tangent to the DRO at the X-axis crossing, a Hohmann-like
  two-impulse transfer with 2 X-axis crossings. Class B has the asymptotic form of family g' (4 crossings), Class C 8 crossings. Class A and C differ in Delta-V by an amount going to zero for large DROs while C takes nearly
  three times as long; B is least preferable (Ocampo, p.14). Numbers (data file): Fig. 15, 600 km altitude parking orbit at 0 deg to a 20,000 km DRO, Class A, 0.596 km/s; Fig. 16a, to a 6,000 km DRO, Class B, 0.343 km/s;
  Fig. 16b, to a 15,000 km DRO, Class C, 0.566 km/s (the caption's second "Class B" is a slip for C). Fig. 15 (viewed): ellipse about 7.5e4 km tall, delta-V1 at the low orbit, delta-V2 at the far crossing.
- Low-thrust (pp.14-17): in the N-body ephemeris many periodic orbits cease to exist, seen as an oscillation of the Jacobi constant (Fig. 17 CAPTION, a 35,000 km near-planar DRO; the short-period oscillation is half the DRO period).
  Spiral-out along the velocity vector from a small retrograde circular orbit (radius under 6,000 km, inclination 180 deg, acceleration under 0.5 mm/s^2) stays in the DRO family (Fig. 18: 100 km altitude start, 90 d burn, 0.1915 mm/s^2,
  concentric DROs up to about 0.4 mm/s^2), so thrust can stop anywhere and leave a DRO; the transfer is not optimal but safe. Above 0.5 mm/s^2 the method fails. Spiral-outs starting beyond 6,000 km or with plane changes need
  manoeuvres. Out-of-plane spiral-outs work down to about 145 deg inclination; Fig. 19 (155 deg) shows a profile shifted off the DRO curve so that stopping between 15,000 and 23,000 km escapes; 140 deg never re-enters.
- Optimised transfers with plane change (p.17, Mystic): Fig. 20 from a 100 km altitude 110 deg orbit (axis label 1,665 km) to a 35,000 km planar DRO, accelerations 0.4069, 0.6103, 1.2207 mm/s^2, Delta-V roughly 1.3 to 1.8 km/s versus
  coast fraction; Fig. 21 from the 35,000 km DRO to a 5,000 km 125 deg orbit, accelerations 0.15 to 0.60 mm/s^2, Delta-V roughly 0.6 to 1.1 km/s. Transfers with many coast arcs can be very unstable; for inclinations 110 to about 135 deg
  the Red Sea plot has very short lifetimes with no obvious way around. Using the DRO profile one can move safely between 5,000 km orbits at 180 to about 135 deg and very distant orbits.
  COMPUTED: 100 km altitude plus Europa's radius is about 1,661 km, consistent with the axis label.
- DRO-type escapes and captures (pp.17-18): a DRO acts as a "coiled spring": a spacecraft leaving a large DRO and swinging by Europa gains a larger departure V-infinity (Whiffen, ref 18, AAS 03-544); the reverse is a very efficient capture,
  compared with capture along an L2 halo's stable manifold, at the price of long spiral-down time. Fig. 22 (viewed, one witness): capture (green, from J. Kangas) and a 35 d spiral-in (blue) from a 30,000 km DRO at 175 deg to a 5,000 km orbit at 125 deg:
  coast 5.909 d, Delta-V 770.5463 m/s, acceleration about 0.3 mm/s^2, mass 18932.5114 kg falling to 18473.2966 kg, final semi-major axis 5001.92 km, inclination 127.34 deg, eccentricity 0.0003196, Europa energy -0.3200 km^2/s^2,
  dates 2 Nov to 7 Dec 2023. COMPUTED: the mass and Delta-V imply an exhaust speed of about 31,380 m/s, a specific impulse of about 3,200 s.

## 5. Conclusion (p.18)

DROs are "among the most stable periodic orbits in any system with third-body perturbations"; they help safe transfer and quarantine design; the study extends to other multi-body systems; out-of-plane DROs and their use in safe transfers
need further work. Acknowledgments name Henon, Ocampo (thesis copy), Sims, Ryan Russell (stability of DROs), and Julie Kangas (Ganymede-to-Europa transfer).

## 6. Equation roles

| Eq. | Role |
|---|---|
| 1 | mass ratio mu, Hill limit mu to 0, mu of Jupiter-Europa printed |
| 2 | Hill's equations about Europa |
| 3 | Hill Jacobi constant |
| 4 | asymptotic family f, an exact 2:1 ellipse of period 2 pi |
| 5, 6 | initial abscissa and velocity for given Gamma on family f |
| 7 | Poincare-section bound on xi' at Y = 0 |

## 7. Figures (one line each)

1 and 2 the 35,000 km DRO (CAPTION). 3 a g3 period-three orbit (CAPTION). 4 DRO families 10,000 to 150,000 km (CAPTION). 5 period versus size (viewed). 6 L2 Lyapunov family a (CAPTION). 7 Hill Poincare section Gamma = 1 (CAPTION). 8 and 9 ephemeris and Hill
sections of 35,000 km DROs (CAPTION). 10 and 11 stability diagram and enlargement with resonances (11 viewed). 12 and 13 out-of-plane stability (13 second viewed). 14 Red Sea plot (viewed). 15 and 16 impulsive transfers (viewed).
17 Jacobi-constant oscillation (CAPTION). 18 and 19 spiral-outs (CAPTION). 20 and 21 optimised transfers (CAPTION and axis text). 22 DRO-type capture plus spiral-in (viewed).

## 8. Relevance to this project

**DROs as quarantine and transfer orbits versus the catalogue's cycler classes.** The paper's DROs are single-body satellites of Europa in the Hill and ephemeris dynamics. They have no encounter sequence, no V-infinity at a
second body, no repetition across bodies and no ballistic flyby chain; the closest the paper comes to transfer-orbit content is the Class A to C impulsive transfers and the spiral-out profile, which are stationary-orbit transfers.
They map to none of the catalogue classes (cycler, quasi_cycler, precursor_mga, mga_tour). They fit the project's existing treatment of raw CR3BP periodic orbits as not cyclers
(literature_check.py header note on halo, DRO, Lyapunov and resonant CR3BP periodic orbits).

**Europa resonant or repeated-moon content.** None in the repeated-moon sense. The resonance content is: (a) the period-three family g3 and instability necks of family f labelled 5:2, 3:2, 4:3 (p.9, Fig. 11), whose ratio
definition is not printed and whose assignment between text and figure is unclear (sec. 3); (b) the L2 and L1 Lyapunov and g' families as transfer references. No Europa-Ganymede or other multi-moon resonance, no Laplace resonance.
Mystic includes the Galilean moons as perturbers (p.4, Fig. 4 caption) but the orbits are about Europa alone.

**Labels the content supports** (label vocabulary read from literature_check.py on 2026-10-11: repeated-moon, mga-tour, resonant, transfer, halo, tulip, nrho, quasi-satellite, retrograde-satellite, pump-tour, model-transition, ephemeris, planar, axisymmetric):
- `retrograde-satellite`: supported (the subject; DROs about Europa, family f).
- `transfer`: supported as secondary content (secs. 4, impulsive Class A to C and low-thrust spirals, DRO-type escape and capture).
- `resonant`: weakly supported only for the instability-neck labelling and period-three family; the paper has no resonant-orbit family study. I would not add it.
- NOT supported: `repeated-moon`, `mga-tour`, `pump-tour`, any cycler or quasi-cycler label.

**Working-body scope.** Europa and Jupiter in the Hill and ephemeris-model analyses; the ephemeris propagation also contains Io, Ganymede, Callisto and the Sun as perturbers. Working bodies (the moons the spacecraft orbits or encounters): Europa only.
For the literature gate: primary Jupiter, body_set {Europa}, single moon. A candidate involving two or more moons is not covered by it.

**Recommendation.** There is no gate anchor for this paper in literature_check.py (grep 2026-10-11). If one is added (for example as prior art for Europa DRO or quarantine candidates under #956 R9 or #949 X4 Jovian families),
use primary "Jupiter", body_set {"Europa"}, topology_label {"retrograde-satellite", "transfer"}, authors ("Lam", "Whiffen"), the citation "Lam, T. & Whiffen, G. J., AAS 05-110, Space Flight Mechanics Meeting, Copper Mountain, Jan 2005"
and no DOI; it should not act as prior art for any cycler, repeated-moon or multi-moon signature. The existing Henon family-f anchor (Sun-Earth) is a different system. The earlier wanted-list line "Family f (DRO) members at Europa;
#956 R9 / X4 family-identity context" is accurate: the paper gives Europa DRO sizes from 10,000 to 150,000 km and stability ranges, which can serve as a family-identity control for a Europa DRO search.

Other uses: (a) Fig. 11's necks give three Europa DRO sizes where a DRO search should expect instability (about 13,000 to 13,500, 21,000 and 24,500 km, approximate from the render, resonance names unreliable). (b) The Hill unit mismatch of 0.6 percent
(sec. 2) is a caution when comparing the paper's x = 1.76579 with our own Hill-frame values. (c) Table-free paper: no Jacobi constants, periods or monodromy values are printed for any DRO (only the two Hill-unit numbers
and Fig. 5 as a curve), so it cannot serve as a numeric golden for a DRO reproduction under the project's sourced-values rule beyond the transfer Delta-Vs of Figs. 15 and 16 and Fig. 22's annotation block.

## 9. Internal inconsistencies and cautions

- Two figures numbered 13 (pp.11 and 12).
- Fig. 16b caption prints "Class C ... Class B impulsive Delta-V = 0.566 km/sec"; the second "Class B" is a slip.
- Resonance names of the necks: text (5:2 and 4:3 for g3, 3:2 Europa-only), the figure's three vertical lines (5:2, 3:2, 4:3) and the introduction (3:1 in Hill's case) are not reconciled, and no ratio definition is given.
- The "100 km altitude" science orbit versus the "1,665 km" axis label of Fig. 20 are consistent only if the label is the orbit radius.
- x = 1.76579 Hill units for 35,000 km versus 1.7770 from our constants.
- The Fig. 22 dates (2023) are an example trajectory, not a mission plan.
- "DROs are among the most stable periodic orbits in any system with third-body perturbations" is a claim of the authors, supported here by Europa runs of 200 d to 1,000 y and Hill's model, not proved.

## 10. Citation mining

All 18 references examined. Held means a file in `cyclers_pdf/papers/` plus an index row (grep and filenames, 2026-10-11).

| Ref | Item | Status |
|---|---|---|
| 1 | Henon 1969 A&A 1:223 | held (digest) |
| 2 | Henon 1970 A&A 9:24 | held (digest) |
| 3 | Benest 1976 A&A 54:563 (IV, three-dimensional stability of plane periodic orbits) | NOT held; not on the wanted list |
| 4 | Broucke 1968 JPL TR 32-1168 | held |
| 5 | Benest 1974 A&A 32:39 (Effects of the mass ratio ..., I) | NOT held; not on the list |
| 6 | Benest 1975 A&A 45:353 (II) | NOT held; not on the list |
| 7 | Whiffen and Sims, AAS 01-209 (Mystic) | NOT held; not on the list |
| 8 | Szebehely 1967 book | held |
| 9 | Scheeres, Guman, Villac, JGCD 24(2):778-787 (2001), stability of planetary satellite orbiters, Europa | NOT held; not on the list |
| 10 | Henon 1997 Generating Families (LNP m52) | held |
| 11 | Ocampo and Rosborough AAS 93-180 (transfer trajectories for distant retrograde orbiters of the Earth) | NOT held; wanted-list rows mention Ocampo only for other papers |
| 12 | Ocampo and Rosborough AAS 93-681 (optimal low-thrust transfers between a class of restricted three-body trajectories) | NOT held |
| 13 | Henon 1973 A&A 28:415 | held |
| 14 | Henon 1974 A&A 30:317 | held |
| 15 | Kechichian et al. 2002, 7th Intl Conf on Libration Point Orbits (DRO and halo picosatellites, solar surveillance) | NOT held; low relevance |
| 16 | Benest 1976 Celest. Mech. 13:203 | held (ADS scan; see the third item of this queue) |
| 17 | Benest 1976 A&A 53:231 (III) | NOT held; not on the list |
| 18 | Whiffen 2003 AAS 03-544 (Jupiter Galilean moon orbiter trajectory; DRO-type escape) | NOT held; not on the list |

### Proposed wanted-list rows (not added; the lead decides)

1. Benest, D. (1974, 1975, 1976a, 1976b), "Effects of the mass ratio on the existence of retrograde satellites in the circular plane restricted problem", I A&A 32:39-46 (1974), II A&A 45:353-363 (1975), III A&A 53:231-236 (1976), IV "Three-dimensional stability of plane periodic orbits" A&A 54:563-568 (1976). Why: the four-part study behind the family-f stability regions this paper reads (refs 3, 5, 6, 17); same acquisition route as the held Benest 1976 Celest. Mech. paper (ADS scans). No identifiers checked here.
2. Ocampo, C. A. & Rosborough, G. W. (1993, 1994), "Transfer trajectories for distant retrograde orbiters of the Earth", AAS 93-180, and "Optimal low-thrust transfers between a class of restricted three-body trajectories", AAS 93-681 (Adv. Astronaut. Sci. 85 Part II, 1994). Why: the Class A, B, C ERO nomenclature and transfer costs the paper quotes (p.13-14); also Ocampo's thesis as a possible source. Low priority for cyclers.
3. Whiffen, G. J. (2003), "An investigation of a Jupiter Galilean moon orbiter trajectory", AAS 03-544, and Whiffen and Sims (2001), AAS 01-209. Why: the DRO-type escape and capture idea (ref 18) and the Mystic algorithm; low priority.
4. Scheeres, D. J., Guman, M. D. & Villac, B. F. (2001), "Stability analysis of planetary satellite orbiters: application to the Europa Orbiter", JGCD 24(2):778-787. Why: Europa stability of low orbits; low priority.
5. Kechichian et al. (2002), libration point orbits conference: not proposed (solar surveillance; no bearing on cyclers).
