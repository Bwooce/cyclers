# Digest: Carruba, Caritá, Aljbaae, Araujo & Domingos (2026), "Co-orbital asteroids of terrestrial planets affected by the von Zeipel-Lidov-Kozai mechanism"

Celestial Mechanics and Dynamical Astronomy 138:8 (2026), DOI 10.1007/s10569-026-10275-5, 22 pages (article pp1-16, appendix and
tables pp17-20, declarations pp21, references pp21-22). V. Carruba (UNESP Guaratinguetá and Laboratório Interinstitucional de
e-Astronomia, Rio de Janeiro), G. Caritá (INPE), S. Aljbaae (INPE and Universidad de Atacama), R. A. N. Araujo (UNESP), R. C. Domingos
(UNESP São João da Boa Vista). Received 22 September 2025, revised 31 December 2025, accepted 5 January 2026, published online
13 February 2026. Open access (CC BY 4.0). Filed in the private paper corpus as
carruba-carita-aljbaae-araujo-domingos-2026-co-orbital-asteroids-terrestrial-planets-von-zeipel-lidov-kozai-cmda-138-8-doi-10.1007-s10569-026-10275-5.pdf
(usable text layer: `pdftotext` gives about 9,790 words).

Digested 2026-10-04 from all 22 pages. Each statement is marked READ (seen on the page, with page, section, equation or table),
COMPUTED (our arithmetic on printed numbers) or INFERRED (our reading). Page numbers are the journal's printed "Page n of 22". Tables
2 to 6 are transcribed from the PDF text layer and checked against the rendered pages; column units are as printed, and where a
caption does not state a unit the unit given in brackets is INFERRED and marked. No table entry is computed or inferred.

## 0. What the paper is

READ (abstract, p1): "We investigate the dynamical behavior of co-orbital asteroids of Venus, Earth, and Mars potentially affected by
the von Zeipel-Lidov-Kozai (ZLK) mechanism. Semi-analytical models of the ZLK dynamics for NEOs in the terrestrial planet regions
predict that at high eccentricities and inclinations, corresponding to i_max = arccos[sqrt(1 - e^2) cos(inc)] > 30 deg, equilibrium
points should appear at 90 deg and 270 deg. At low eccentricities and inclinations in the vicinity of Earth and Venus, perturbations
from these two planets are dominant, with two equilibrium points occurring at 0 deg and 180 deg. Using numerical simulations performed
with the SWIFT and REBOUND integrators, we classify the Kozai states of known co-orbitals and assess their long-term stability. Our
results indicate that no current Venus co-orbital is a robust ZLK resonator, while several Earth co-orbitals and one Mars co-orbital
(2017 XG62) exhibit libration around one of the four ZLK equilibrium points. Libration around the 0 deg and 180 deg equilibrium points
protects asteroids from close encounters with their perturbing planet, Earth, enhancing long-term stability on timescales of about 10^5
years. 2017 XG62 is the only known co-orbital asteroid of a terrestrial planet that is currently librating around the high-inclination
90 deg equilibrium point. This configuration prevents this asteroid from experiencing close encounters with other terrestrial planets."

In one paragraph: a classification study of natural objects. It integrates the known co-orbitals of Venus (22 objects), Earth (69
objects, 54 from a published catalogue plus 15 temporary ones) and Mars (48 objects) for 78,000 years with two N-body integrators,
labels each object by whether and about which equilibrium point its argument of perihelion omega librates, repeats the analysis on
729 clones per object to account for orbit uncertainty, and counts close encounters over 10^5 years for the librating Earth and Mars
objects against a control sample of low-i_max Earth co-orbitals. It is an observational-dynamics paper about asteroids; it has no
spacecraft, trajectory design, delta-v or periodic-orbit family in it.

## 1. The model and the numerical set-up

### 1.1 Co-orbital types, as defined (Section 1, p2)

READ. A co-orbital asteroid "is a small celestial body whose orbital motion is synchronized with a planet due to the existence of
stable equilibrium points associated with a 1:1 mean-motion resonance"; the critical angle sigma is "the difference between the mean
longitudes of the small body and the planet". Five types: "Tadpoles (T): These objects librate (oscillate) around the known L4 and L5
Lagrangian equilibrium points at +/- 60 deg"; "Quasi-satellites (QS): Also known as retrograde satellites. The critical angle librates
around 0 deg"; "Horseshoes (HS): These orbits wrap around both L4 and L5, as well as the unstable L3 equilibrium point (at sigma about
180 deg)"; "H-QS compound orbits, which are mergers of quasi-satellites and horseshoes"; "T-QS compound orbits, formed by combining
quasi-satellites and tadpoles". Also "sticking" orbits (those that evolved on just one side of the resonance and for which the change
in semi-major axis was at least one-fifth of the resonance width), as defined by Pan and Gallardo (2025) (p2).

READ (p2-3). Populations quoted from the literature: Venus 20 known (Carruba et al. 2024; Pan and Gallardo 2025); Earth 54 identified
plus 15 possible temporary co-orbitals; Mars 48, "including the well-known Eureka and 1999 UJ7"; 2023 FW14 "confirmed as the second
Mars L4 tadpole, showing stability for at least 10,000 years"; "A new Mars quasi-satellite, 2021 FV1". Mars tadpoles show an
L4/L5 asymmetry (more at L4), unlike Jupiter and Neptune; possible origins named: rotational fission differences or impact ejecta
from Mars (Polishook et al. 2017). INFERRED/COMPUTED note: the Venus count in the text (20) differs from the 22 listed in Table 4,
presumably because two were added after the cited catalogues; the paper does not comment.

### 1.2 The ZLK mechanism and the semi-analytical model (Section 2, pp3-5)

READ. The ZLK mechanism "describes the coupled oscillations of eccentricity and inclination in the orbit of a body that is subject
to the gravitational perturbation of an inclined massive object" (p3); eccentricity and inclination exchange while the semi-major axis
stays nearly constant. Michel and Thomas (1996) showed for near-Earth asteroids with a < 2 au that two additional equilibrium points
near 0 deg and 180 deg appear besides 90 deg and 270 deg because of Venus and Earth perturbations (p3).

The secular model is the one-degree-of-freedom Hamiltonian in the Delaunay variable g = omega for the generalised circular inclined
problem; "all planets are assumed to move on circular orbits in the invariant plane of the Solar System, in which both H and L are
preserved" (p4), with L = sqrt(a), H = sqrt(a(1 - e^2)) cos(inc). The conserved quantity is sqrt(1 - e^2) cos(inc) and

    i_max = arccos[ sqrt(1 - e^2) cos(inc) ],                                                                  (1)   p4

"which is the value of inclination when e = 0". Michel and Thomas (1996) computed "a semi-analytical averaged Hamiltonian for the
ZLK mechanism while accounting for the contribution of each planet." Figure 1 (p5, adapted with permission from Michel and Thomas
1996 Fig. 1, "AI-enhanced from") shows level curves of the averaged Hamiltonian in the (e cos omega, e sin omega) plane for a =
0.98 au (just outside the 1:1 resonance) and six values of i_max; the panel labels legibly show 1, 2, 20 and 30 deg for panels (a) to (d), and the labels of (e) and (f) are
too small to transcribe reliably. READ (p4) description: "For i_max = 2 deg, two small islands of libration appear at omega = 0 deg and omega = 180 deg. As
i_max increases, the phase space becomes more stretched, and the libration islands expand. At i_max = 30 deg, the equilibrium point
e = 0 is unstable and a separatrix curve passes through it. This curve defines two new islands of libration in omega, with stable
equilibrium points at omega = 90 deg and omega = 270 deg. As i_max grows, the centers of the islands at omega = 0 deg and 180 deg
become unstable equilibrium points. Finally, at very high i_max values, these equilibrium points disappear." For Venus the dynamics
are comparable; for Mars "no libration regime appears in the (e, omega) plane for 0 < i_max < 30 deg. For i_max > 30 deg, Jupiter's
gravity alone drives the resonant behavior" (p4). The caption states the Earth and Venus orbit-crossing requirements are shown by
the dashed and solid lines.

READ (p3-4). Nesvorný et al. (2002) studied secular dynamics inside co-orbital motion and found that the Kozai-libration islands
at omega found by Kozai (1985) "were likely computational artifacts due to Kozai's approximation"; "no omega libration exists at the
asymmetric stationary point L4 and L5, even at high eccentricities and inclinations" (p3). The paper's inference (p4): "if the model
of Nesvorný et al. (2002) holds in the inclined case, the ZLK equilibrium points at 90 deg and 270 deg should disappear for TL4 and
TL5 orbits. For other co-orbital configurations, equilibrium points should still occur at multiples of 90 deg, although the Hamiltonian
contour levels will differ from the secular case." The Michel-Thomas model "is exclusively valid for objects outside of mean-motion
resonances" and is used only for comparison, not quantitatively (p4). There is no new equation, Hamiltonian or semi-analytical model
in this paper beyond (1) and the Hill-radius expression (2); the "semi-analytical" part is a review of Michel and Thomas (1996).

### 1.3 Numerical integrations (Section 3, p6)

READ. Two packages, to avoid reliance on one integrator: "SWIFT (Levison and Duncan 1994) and REBOUND (Rein and Liu 2012)". SWIFT is
described as "based on the Wisdom-Holman symplectic mapping method", and its Bulirsch-Stoer option is used: "We performed numerical
integrations using the SWIFT-BS code, a Bulirsch-Stoer integrator from the SWIFT package (Levison and Duncan 2013), with a tolerance
of 10^-8." REBOUND uses "REBOUND's Bulirsch-Stoer integrator to maintain consistency with the SWIFT simulations". Set-up: "Each
asteroid's orbit was integrated for 78,000 years with a 1-day time step, under the gravitational influence of all the planets and the
Moon." The span "was chosen to ensure that even asteroids with very long ZLK periods (e.g., 2008 WM64) complete at least one full Kozai
cycle." Elements "obtained from the JPL Horizons program at the Julian Date (JD) of 2460800.5, corresponding to midnight of May 5th,
2025" and "transformed to the Solar System's invariant plane reference frame". Masses (Appendix, p17), in solar-mass units:

| Body | Mass [solar masses] |
|---|---|
| Mercury | 1.6601208254808336e-07 |
| Venus | 2.447838287784771e-06 |
| Earth | 3.0034896154502038e-06 |
| Moon | 3.694303350091508e-08 |
| Mars | 3.2271559174983593e-07 |
| Jupiter | 9.545942479871165e-04 |
| Saturn | 2.8581500081698117e-04 |
| Uranus | 4.365793681773934e-05 |
| Neptune | 5.1503084159015005e-05 |

READ (p6-7). Classification by the omega time series (Table 1): objects whose omega circulates 0 to 360 deg are label 0; librating
about 0 deg label 5 (alternating) or 10 (librating); 180 deg labels 15 and 20; 90 deg labels 25 and 30; 270 deg labels 35 and 40.
"Overall, labels that are multiples of 10 correspond to stable libration states, while labels ending in 5 indicate intermittent
(temporary) libration." "Asteroids for which there was disagreement between the two classifications were considered dubious. In this
work, we classify them as alternating between states, and we further investigate them only in case one of the two classifications is a
pure libration, and neglect the other cases." Figures 2 and 3 (pp7-8) show example (e cos omega, e sin omega), (omega, e) and
(omega, t) plots for one object in each of the 5, 10, 15, 20, 25, 30, 35 and 0 states.

**Table 1 (p9): the type of resonant behavior for the ZLK mechanism for terrestrial planets** (caption: "We report the equilibrium
point (None for circulating orbits), the status (A for alternating, L for librating), the color used in our figures, and the index
used to identify each orbit type").

| Eq. point | Status | Colour | Index |
|---|---|---|---|
| None | None | Black | 0 |
| 0 deg | A | Cyan | 5 |
| 0 deg | L | Blue | 10 |
| 180 deg | A | Bluish green | 15 |
| 180 deg | L | Green | 20 |
| 90 deg | A | Orange | 25 |
| 90 deg | L | Red | 30 |
| 270 deg | A | Vermillion | 35 |
| 270 deg | L | Reddish purple | 40 |

### 1.4 Statistical clone analysis (Section 4, pp11-13)

READ (p11-12). Orbit uncertainty is measured with the Orbital Condition Code (OCC), "a scale of 0 to 9: 0 denotes a well-established
orbit (with < 1 arcsec of longitude drift per decade), while 9 indicates an extremely uncertain orbit". Covariance matrices were not
available for every asteroid, so: "Using JPL's Horizons system (Small-Body Database, accessed August 20, 2024), we first obtained each
asteroid's orbital elements along with their stated uncertainties [...] For each orbital element we used three values: the nominal
value, as well as the nominal value plus its 1 sigma error and minus its 1 sigma error. This yields 3^6 = 729 possible sets of orbital
elements (a, e, i, Omega, omega, M) for each asteroid. We integrated each set of 729 clone orbits under the gravitational influence
of all eight planets using the SWIFT-BS integrator with the same settings described in Sect. 3. At this stage, for simplicity and
consistency, we restricted our analysis to the SWIFT-BS results and discontinued use of the REBOUND integrator." For each ensemble the
mean omega <omega> and the oscillation amplitude Amp. (maximum minus minimum omega over the 78,000 years) are computed, with their
standard deviations across the 729 clones as uncertainties; for objects near the 0 deg equilibrium (labels 5 and 10) omega(t) is first
offset by +180 deg before averaging and the 180 deg removed afterwards. Figure 7 (p12) shows histograms of <omega> and Amp. for the
clones of 2017 XG62 (all clones remain in libration near omega = 90 deg) and of 2002 AA29 (clones alternate between libration and
circulation, <omega> near 180 deg, Amp. often near 360 deg).

### 1.5 Close-encounter test (Section 5, pp13-16)

READ (p15). Encounters are counted when the test asteroid comes within a planet's Hill radius,

    R_H = a_p (1 - e_p) cbrt( M_p / (3 (M_p + M_sun)) ),                                                      (2)   p15

"For reference, R_H for Venus, Earth, the Moon, and Mars are 0.0067, 0.0098, 0.0004, and 0.0073 au, respectively." The Moon is included
"because recent work suggests some Earth co-orbital asteroids may have originated from a lunar impact (Sfair et al. 2025)"; "We also
applied the low-velocity encounter detection techniques from Carruba et al. (2020)". Each simulation is "extended to 10^5 years using the
same scheme as in Sect. 3. Extending beyond 10^5 years was deemed not meaningful, given the highly chaotic nature of many terrestrial
co-orbitals." The comparison sample is the Earth co-orbitals with i_max < 9.5 deg (2003 YN107, 2006 JY26, 2009 SH2, 2013 BS45,
2014 QD364, 2019 XH2), chosen "to ensure our non-resonant sample lies in the core region where libration around the 0 deg/180 deg
equilibria is theoretically possible".

## 2. What is found

### 2.1 Venus (Section 3.1, pp8-9; Fig. 4 p9; Table 4)

READ. Observational bias against Venus co-orbitals with eccentricity below 0.38 ("the value for which the apocenter of the asteroid is
equal to the pericenter of Earth"); all known Venus co-orbitals have i_max larger than 20 deg (Table 4). "While some asteroids appear to
be in alternating ZLK states according to the REBOUND simulations, this is not confirmed by the SWIFT simulations, which may be caused
by sensitivity near resonance borders. No asteroid was found in a pure librating state. Based on this analysis, we conclude that there
are no convincing candidates for the ZLK mechanism among the current population of co-orbital asteroids of Venus." In Table 4 the
nonzero labels are 322756 (R = 35), 2014 JU15 (R = 35) and 2022 BL5 (R = 15), all from REBOUND (S = 0 for all 22).

### 2.2 Earth (Section 3.2, pp9-10; Fig. 5 p10; Table 5)

READ. 54 known (Pan and Gallardo 2025) plus 15 possible temporary co-orbitals from other sources (de la Fuente Marcos and de la
Fuente Marcos papers; Kaplan and Cengiz 2020; Dvorak et al. 2012; Carruba et al. 2025 for 2025 FP4; personal communications for 2025 SC
and 2025 RB7). READ (p10): "Several asteroids co-orbital to Earth (such as tadpoles, horseshoes, or quasi-satellites) are known to
exhibit ZLK behavior: for example, asteroid 2015 SO2, an Earth's horseshoe librator, is subject to ZLK mechanism, with omega
oscillating around 270 deg", and three quasi-satellites (164207 2004 GU9, 277810 2006 FV35, 2013 LX28) show "ZLK-like" dynamics around
270, 180 and 0 deg. "Results from Dvorak et al. (2012) show that above about 40 deg the ZLK mechanism greatly amplifies e and i, leading to
orbital instability." Results (p10): 2015 SO2 gets label 35, 2004 GU9 label 35, 2006 FV35 label 15, 2013 LX28 label 5; "we also identify 17
other possible resonators" (preliminary classification, Table 5 and p10). Differences from the secular non-resonant model: "2013 LX28
has a high value of i_max, larger than 50 deg, for which we would not expect libration around 0 deg. Conversely, 2002 AA29 (label
25) and 2001 GO2, 2004 GU9, 2015 SO2, and 2023 GC2 (label 35) all have values of i_max lower than 30 deg, for which we would not expect
equilibrium points at 90 deg." (The text there puts 2002 AA29 in the label-35 group while giving its label as 25; Table 2 has it at 25.
Typesetting slip.)

### 2.3 Mars (Section 3.3, p11; Fig. 6; Table 6)

READ. 48 known Mars co-orbitals; "The stable Trojan region is concentrated at inclinations of about 15 deg - 30 deg (Scholl et al. 2005).
Above about 35 deg, the ZLK mechanism becomes dominant and rapidly removes objects (increasing to planet-crossing orbits). Thus,
most of the known Mars tadpoles reside in long-term stable orbits [...] because they are outside the destabilizing ZLK zone. 2017
XG62 is a recently discovered co-orbital asteroid that may escape this rule, being at high inclination (inc = 42.148 deg)." "We
found a new asteroid librating around 90 deg (label 30), 2017 XG62, not previously reported in the literature."

### 2.4 Clone statistics: the 21 candidates (Section 4, p13; Table 2)

READ (p13). "We identify the following resonant asteroids: Label 10: 0 deg equilibrium point: 2012 FC71 and 2021 VU12. Label 20: 180 deg
equilibrium point: 2008 WM64, 2019 NC1, 2020 DX1, and 2022 UO10. Label 30: 90 deg equilibrium point: 2017 XG62. Label 40: 270 deg
equilibrium point: 2016 CA138." Further: "2008 WM64 is the librating asteroid with label 20 with the highest value of i_max (34.020 deg).
For such a high value of i_max, the equilibrium point at 180 deg is still present, but the size of the libration island is rather
small (see Fig. 1). [...] the asteroid is the only one to conclude a single libration cycle" over the 78,000 years. 2017 XG62 is "the
only Martian co-orbital asteroid found to be in a ZLK state" and "the only one identified to librate around 90 deg, in a ZLK regime
dominated by Jupiter"; it is classified a QS in Pan and Gallardo (2025). "2016 CA138 is the only object that may oscillate around
270 deg. While this classification is not certain, since about 40% of the clones of the asteroids did not remain in such a
configuration, its value of i_max is compatible with a high-inclination ZLK dynamics. 2016 CA138 is currently in an HS-QS orbit,
according to Pan and Gallardo (2025)." "In terms of locations of equilibrium points, our findings are qualitatively compatible with
the semi-analytical model of Michel and Thomas (1996). The disagreements with the model from the preliminary classification discussed
in Sect. 3.2 are not confirmed by our statistical studies."

### 2.5 Role of the mechanism: protection from close encounters (Section 5, pp13-16)

READ (p15-16). "As shown in Table 3, none of the asteroids in the 10 or 20 states experienced close encounters with Earth, whereas the
asteroids not affected by the ZLK mechanism did. A similar protective effect is hinted for the one object tentatively in the 40 state
(270 deg libration), as it likewise showed no Earth encounters. The ZLK mechanism does not shield asteroids from Moon encounters;
however, such encounters are extremely infrequent in our simulations. Indeed, only two asteroids (2003 YN107 and 2008 WM64)
experienced even a single close encounter with the Moon over 10^5 years, and in both cases the approach distances were not very
small." For Mars (p16): 2017 XG62 (i_max = 42.148 [text; Table 6 prints 46.727 in the i_max column, see section 3 below]) "experienced no
close encounters with Earth or Venus. None of these asteroids had any close encounters with Mars during the simulation"; the
high-i_max comparison objects 2019 BG3 and 2020 LE1 are not in the ZLK state and did have encounters (Table 3).

Alignment explanation (p3): librating around 0 deg or 180 deg aligns "the line of nodes and line of apsides. This causes the
intersections of the object's orbit with the orbits of the planets to occur either at its periapsis or apoapsis, maximizing the distance
to the planet."

READ (conclusion, p16): eight robust ZLK resonators ("objects librating around omega = 0 deg, 180 deg, 90 deg, and, in one case,
potentially 270 deg"); none among the Venus co-orbitals ("This absence could be due to observational biases [...] any that librate
around 90 deg or 270 deg could be intrinsically rare"); "several Earth co-orbitals librating around the 0 deg and 180 deg points -
a configuration that shields them from close Earth encounters and stabilizes their orbits over about 10^5 years"; all found to be
"sticking" 1:1 orbits per Pan and Gallardo (2025) "meaning that the secular model of Michel and Thomas (1996) should be applicable to
their dynamics"; "We did not identify any TL4 or TL5 co-orbital asteroids currently librating in ZLK states for terrestrial planets,
consistently with the results of Nesvorný et al. (2002)." For Mars: 2017 XG62, "which is currently on a transitional HS/QS
co-orbital (Pan and Gallardo 2025)", is "the only known terrestrial-planet co-orbital presently librating around the 90 deg Kozai
equilibrium. This high-inclination ZLK state insulates 2017 XG62 from close encounters with the other terrestrial planets."

COMPUTED note on the numbers quoted in the text: the conclusion and Section 5 call the Mars object's i_max 42.148 where Table 2, Table
3 and Table 6 print 46.727 in the i_max column; 42.148 is the printed inclination (Table 6, "inc" column) of 2017 XG62, so the text
appears to use inc where i_max is meant. The text and Table 2 are consistent for every other object checked.

## 3. Every table, transcribed

The tables are reproduced in the printed order and with the printed column headings. Table 2 columns (p14): "#", "P.L." (preliminary
label from Table 1), "Ast. ID", "OCC" (orbital condition code), "i_max [degrees]", "P" (the planet the asteroid is co-orbital to: E
Earth, M Mars), "<omega> [degrees]", its error, "Amp. [degrees]" and its error, and "F.L." (final label). In Tables 4, 5, 6 the
captions say "We report the number, asteroid ID, a, e, inc, omega, i_max and the orbital status label obtained from a SWIFT (S) and a
REBOUND (R) integration" (Table 4 and 6 captions; the Table 5 caption says "the orbital status label"); the omega column is NOT
printed in the tables. Units are not given in the captions: a is INFERRED to be in au (values 0.72 to 1.28, with the planets' semi-major
axes 0.723, 1.000, 1.524), inc and i_max in degrees (in Table 2 the i_max header is "[degrees]").
The S and R columns hold the Table 1 index labels (0, 5, 10, ... 40).

**Table 2 (p14)**, caption: List of the co-orbital asteroids of terrestrial planets that are likely in ZLK states. F.L. is the final label after the clone analysis. Entries marked "40?" are as printed.

| # | P.L. | Ast. ID | OCC | i_max [deg] | P | <w> [deg] | err(<w>) [deg] | Amp. [deg] | err(Amp) [deg] | F.L. |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 5 | 2013LX28 | 1 | 54.997 | E | 22 | 20 | 346 | 51 | 5 |
| 2 | 5 | 2014EK24 | 0 | 6.2643 | E | 167 | 170 | 355 | 29 | 5 |
| 3 | 5 | 2020PN1 | 0 | 8.7497 | E | 171 | 163 | 351 | 31 | 5 |
| 4 | 10 | 2012FC71 | 5 | 7.069 | E | 56 | 130 | 100 | 3 | 10 |
| 5 | 10 | 2021GN1 | 8 | 28.874 | E | 12 | 1 | 360 | 1 | 5 |
| 6 | 10 | 2021VU12 | 1 | 19.434 | E | 351 | 22 | 110 | 17 | 10 |
| 7 | 10 | 2022UP20 | 2 | 18.859 | E | 259 | 147 | 247 | 105 | 5 |
| 8 | 15 | 2006FV35 | 1 | 23.235 | E | 178 | 15 | 354 | 29 | 15 |
| 9 | 15 | 2019GM1 | 0 | 7.9056 | E | 178 | 15 | 365 | 24 | 15 |
| 10 | 15 | 2022UY | 1 | 16.878 | E | 179 | 15 | 202 | 115 | 15 |
| 11 | 20 | 2008WM64 | 0 | 34.020 | E | 147 | 1 | 160 | 2 | 20 |
| 12 | 20 | 2019NC1 | 0 | 9.9915 | E | 182 | 1 | 110 | 1 | 20 |
| 13 | 20 | 2020DX1 | 7 | 13.857 | E | 177 | 1 | 136 | 1 | 20 |
| 14 | 20 | 2022UO10 | 7 | 18.324 | E | 186 | 21 | 225 | 90 | 20 |
| 15 | 25 | 2002AA29 | 0 | 10.773 | E | 182 | 18 | 354 | 24 | 25 |
| 16 | 30 | 2017XG62 | 0 | 46.727 | M | 90 | 1 | 11 | 3 | 30 |
| 17 | 35 | 2001GO2 | 7 | 10.720 | E | 182 | 13 | 360 | 3 | 35 |
| 18 | 35 | 2004GU9 | 0 | 15.696 | E | 197 | 27 | 350 | 42 | 35 |
| 19 | 35 | 2015SO2 | 1 | 11.068 | E | 200 | 30 | 350 | 42 | 35 |
| 20 | 35 | 2016CA138 | 0 | 27.848 | E | 250 | 37 | 182 | 94 | 40? |
| 21 | 35 | 2023GC2 | 1 | 10.921 | E | 196 | 26 | 355 | 22 | 35 |

**Table 3 (p15)**, caption: planet to which the asteroid is co-orbital (E Earth, M Mars), its identification, i_max, final label and the number of encounters with Earth, Moon and Venus during the simulations. "The last three rows refer to the simulations with Martian co-orbital asteroids."

| Planet | Ast. ID | i_max [deg] | F.L. | # enc Earth | # enc Moon | # enc Venus |
|---|---|---|---|---|---|---|
| E | 2003YN107 | 4.3939 | 0 | 2008 | 1 | 0 |
| E | 2006JY26 | 4.9745 | 0 | 1638 | 0 | 0 |
| E | 2009SH2 | 8.6893 | 0 | 489 | 0 | 0 |
| E | 2013BS45 | 4.866 | 0 | 944 | 0 | 0 |
| E | 2014QD364 | 4.6603 | 0 | 1545 | 0 | 0 |
| E | 2019XH2 | 9.4888 | 0 | 183 | 0 | 0 |
| E | 2012FC71 | 7.069 | 10 | 0 | 0 | 0 |
| E | 2021VU12 | 19.434 | 10 | 0 | 0 | 0 |
| E | 2008WM64 | 34.020 | 20 | 0 | 1 | 0 |
| E | 2019NC1 | 9.9915 | 20 | 0 | 0 | 0 |
| E | 2020DX1 | 13.857 | 20 | 0 | 0 | 0 |
| E | 2022UO10 | 18.324 | 20 | 0 | 0 | 0 |
| E | 2016CA138 | 27.848 | 40? | 0 | 0 | 0 |
| M | 2019BG3 | 53.56 | 0 | 5 | 0 | 1 |
| M | 2020LE1 | 40.526 | 0 | 6 | 0 | 8 |
| M | 2017XG62 | 46.727 | 30 | 0 | 0 | 0 |

**Table 4 (p17)**, caption: List of the currently known co-orbital asteroids of Venus

| # | Ast. ID | a [au] | e | inc [deg] | i_max [deg] | S | R |
|---|---|---|---|---|---|---|---|
| 1 | 322756 | 0.72455 | 0.38268 | 8.1313 | 23.852 | 0 | 35 |
| 2 | 524522 | 0.72359 | 0.41018 | 9.0368 | 25.752 | 0 | 0 |
| 3 | 2012XE133 | 0.7232 | 0.43262 | 6.7272 | 26.444 | 0 | 0 |
| 4 | 2013ND15 | 0.7236 | 0.61184 | 4.7986 | 37.982 | 0 | 0 |
| 5 | 2014JU15 | 0.72605 | 0.46649 | 9.1023 | 29.145 | 0 | 35 |
| 6 | 2015WZ12 | 0.72161 | 0.4127 | 3.6358 | 24.628 | 0 | 0 |
| 7 | 2019HH3 | 0.72153 | 0.51333 | 13.561 | 33.461 | 0 | 0 |
| 8 | 2019XQ | 0.71786 | 0.42729 | 7.0382 | 26.194 | 0 | 0 |
| 9 | 2020CL1 | 0.72354 | 0.42003 | 4.0941 | 25.151 | 0 | 0 |
| 10 | 2020HY6 | 0.72447 | 0.39191 | 45.482 | 49.831 | 0 | 0 |
| 11 | 2020QU5 | 0.73794 | 0.43131 | 9.3735 | 27.107 | 0 | 0 |
| 12 | 2020SB | 0.72231 | 0.58177 | 0.769 | 35.582 | 0 | 0 |
| 13 | 2020SV5 | 0.72563 | 0.36158 | 5.7504 | 21.929 | 0 | 0 |
| 14 | 2021VH1 | 0.72186 | 0.53163 | 14.452 | 34.898 | 0 | 0 |
| 15 | 2021XA1 | 0.71496 | 0.40121 | 20.624 | 30.987 | 0 | 0 |
| 16 | 2021XO3 | 0.72295 | 0.40951 | 6.2924 | 24.932 | 0 | 0 |
| 17 | 2022BL5 | 0.72475 | 0.48359 | 4.9457 | 29.304 | 0 | 15 |
| 18 | 2022CD | 0.72354 | 0.53038 | 12.21 | 34.047 | 0 | 0 |
| 19 | 2022UA13 | 0.72438 | 0.68265 | 6.6723 | 43.465 | 0 | 0 |
| 20 | 2023BB1 | 0.72374 | 0.58544 | 8.2493 | 36.647 | 0 | 0 |
| 21 | 2023QS7 | 0.7227 | 0.58983 | 10.8 | 37.512 | 0 | 0 |
| 22 | 2024AF6 | 0.72557 | 0.40722 | 15.055 | 28.117 | 0 | 0 |

**Table 5 (pp18-19)**, caption: List of the currently known Earth co-orbital asteroids. "The last 15 entries are asteroids in temporary co-orbital status that were not reported by Pan and Gallardo (2025)."

| # | Ast. ID | a [au] | e | inc [deg] | i_max [deg] | S | R |
|---|---|---|---|---|---|---|---|
| 1 | 3753 | 0.99771 | 0.51486 | 19.806 | 36.238 | 0 | 0 |
| 2 | 2001GO2 | 1.0067 | 0.16818 | 4.6238 | 10.72 | 35 | 35 |
| 3 | 2002AA29 | 0.99252 | 0.012986 | 10.748 | 10.773 | 25 | 25 |
| 4 | 2004GU9 | 1.0013 | 0.13612 | 13.65 | 15.696 | 35 | 40 |
| 5 | 2005QQ87 | 1.0005 | 0.30205 | 33.964 | 37.754 | 0 | 0 |
| 6 | 2006FV35 | 1.0013 | 0.37752 | 7.104 | 23.235 | 15 | 15 |
| 7 | 2008WM64 | 1.0048 | 0.10678 | 33.528 | 34.02 | 20 | 20 |
| 8 | 2009HE60 | 0.99585 | 0.26472 | 1.5833 | 15.43 | 0 | 0 |
| 9 | 2010NY65 | 0.99959 | 0.3697 | 11.67 | 24.502 | 0 | 0 |
| 10 | 2010SO16 | 1.0029 | 0.075391 | 14.52 | 15.136 | 0 | 0 |
| 11 | 2010TK7 | 0.99943 | 0.19045 | 20.894 | 23.488 | 0 | 0 |
| 12 | 2013LX28 | 1.0016 | 0.45207 | 49.978 | 54.997 | 5 | 5 |
| 13 | 2014EK24 | 1.0071 | 0.070189 | 4.8042 | 6.2643 | 5 | 5 |
| 14 | 2014OL339 | 0.99924 | 0.46077 | 10.187 | 29.129 | 0 | 0 |
| 15 | 2015SO2 | 0.99641 | 0.10857 | 9.1643 | 11.068 | 35 | 35 |
| 16 | 2015XX169 | 1.0023 | 0.18474 | 7.6119 | 13.062 | 0 | 0 |
| 17 | 2015YA | 0.99559 | 0.27944 | 1.6182 | 16.305 | 0 | 0 |
| 18 | 2016CA138 | 1.0006 | 0.047868 | 27.723 | 27.848 | 40 | 35 |
| 19 | 2016CO246 | 1 | 0.12483 | 6.3917 | 9.5949 | 0 | 0 |
| 20 | 2016FU12 | 1.0023 | 0.16563 | 2.0736 | 9.7545 | 0 | 0 |
| 21 | 469219 | 1.0011 | 0.10385 | 7.7759 | 9.7866 | 25 | 0 |
| 22 | 2016JP | 0.99456 | 0.38315 | 11.366 | 25.1 | 0 | 0 |
| 23 | 2017SL16 | 1.0013 | 0.15372 | 8.639 | 12.338 | 0 | 0 |
| 24 | 2017XQ60 | 1.001 | 0.21389 | 27.202 | 29.678 | 0 | 0 |
| 25 | 2018AN2 | 1.001 | 0.1542 | 22.085 | 23.717 | 0 | 0 |
| 26 | 2018XW2 | 0.99886 | 0.302 | 19.715 | 26.178 | 0 | 0 |
| 27 | 2019GM1 | 1.0049 | 0.067671 | 6.893 | 7.9056 | 15 | 15 |
| 28 | 2019HS2 | 0.99666 | 0.21395 | 19.597 | 23.036 | 0 | 0 |
| 29 | 2019NC1 | 0.99364 | 0.12713 | 6.8366 | 9.9915 | 20 | 20 |
| 30 | 2019SB6 | 0.99759 | 0.266 | 7.2072 | 16.99 | 0 | 0 |
| 31 | 2019VL5 | 1.0029 | 0.27916 | 1.5448 | 16.281 | 0 | 0 |
| 32 | 2019XH2 | 1.0074 | 0.15289 | 3.5765 | 9.4888 | 0 | 0 |
| 33 | 2019YB4 | 1.0047 | 0.19279 | 0.42825 | 11.124 | 0 | 0 |
| 34 | 2020CX1 | 0.99576 | 0.16321 | 12.748 | 15.789 | 0 | 0 |
| 35 | 2020DX1 | 0.99273 | 0.14109 | 11.274 | 13.857 | 20 | 20 |
| 36 | 2020HE5 | 1.0056 | 0.07706 | 7.1635 | 8.4111 | 0 | 15 |
| 37 | 2020PN1 | 0.99504 | 0.12783 | 4.7689 | 8.7497 | 10 | 5 |
| 38 | 2020PP1 | 1.0027 | 0.073168 | 5.8072 | 7.1602 | 0 | 25 |
| 39 | 2020XL5 | 1.0007 | 0.38724 | 13.846 | 26.467 | 0 | 0 |
| 40 | 2021BA | 0.99524 | 0.23037 | 12.495 | 18.185 | 0 | 0 |
| 41 | 2021GN1 | 0.99783 | 0.19041 | 26.874 | 28.874 | 10 | 10 |
| 42 | 2021VU12 | 1.006 | 0.14717 | 17.56 | 19.434 | 10 | 10 |
| 43 | 2021XS4 | 0.99673 | 0.18796 | 14.687 | 18.179 | 0 | 0 |
| 44 | 2022UO10 | 1.0049 | 0.11191 | 17.197 | 18.324 | 20 | 20 |
| 45 | 2022UY | 0.99749 | 0.11392 | 15.593 | 16.878 | 20 | 15 |
| 46 | 2022UP20 | 1.0029 | 0.13562 | 17.227 | 18.859 | 10 | 10 |
| 47 | 2022VR1 | 0.99515 | 0.17184 | 5.4144 | 11.266 | 0 | 0 |
| 48 | 2022YG | 0.99533 | 0.19575 | 2.332 | 11.524 | 0 | 0 |
| 49 | 2023FW13 | 1.0021 | 0.17785 | 2.7424 | 10.601 | 25 | 35 |
| 50 | 2023GC2 | 0.99397 | 0.14847 | 6.8353 | 10.921 | 35 | 35 |
| 51 | 2023QR1 | 1.0076 | 0.15135 | 4.7758 | 9.9203 | 0 | 0 |
| 52 | 2023TG14 | 1.0068 | 0.21455 | 4.4077 | 13.139 | 0 | 0 |
| 53 | 2024AV2 | 1.0042 | 0.25316 | 5.7005 | 15.711 | 0 | 0 |
| 54 | 2024JR16 | 1.003 | 0.35125 | 5.9303 | 21.366 | 0 | 0 |
| 55 | 54509 | 1.0062 | 0.23019 | 1.5993 | 13.402 | 0 | 0 |
| 56 | 1991VG | 1.0285 | 0.051551 | 1.4366 | 3.2854 | 0 | 0 |
| 57 | 2003YN107 | 0.98866 | 0.013932 | 4.321 | 4.3939 | 0 | 0 |
| 58 | 2006JY26 | 1.0103 | 0.083023 | 1.4391 | 4.9745 | 0 | 0 |
| 59 | 2009SH2 | 0.9913 | 0.094254 | 6.8112 | 8.6893 | 0 | 0 |
| 60 | 2012FC71 | 0.98777 | 0.088206 | 4.9424 | 7.069 | 10 | 10 |
| 61 | 2013BS45 | 0.99189 | 0.083755 | 0.77252 | 4.866 | 0 | 0 |
| 62 | 2014QD364 | 0.98612 | 0.041482 | 4.0094 | 4.6603 | 0 | 0 |
| 63 | 2014UR | 0.99561 | 0.016368 | 8.2389 | 8.2917 | 15 | 5 |
| 64 | 2015YQ1 | 1.0036 | 0.40401 | 2.4844 | 23.951 | 0 | 0 |
| 65 | 2020TK7 | 0.97509 | 0.31539 | 3.7726 | 18.754 | 0 | 0 |
| 66 | 2025FP4 | 1.0013 | 0.26861 | 11.002 | 18.995 | 0 | 0 |
| 67 | 2025PN7 | 1.0030 | 0.10750 | 1.9795 | 6.480 | 0 | 15 |
| 68 | 2025SC | 0.99201 | 0.09305 | 3.727 | 6.5081 | 0 | 15 |
| 69 | 2025RB7 | 0.99848 | 0.29283 | 16.676 | 23.659 | 0 | 0 |

**Table 6 (pp19-20)**, caption: List of the currently known co-orbital asteroids of Mars

| # | Ast. ID | a [au] | e | inc [deg] | i_max [deg] | S | R |
|---|---|---|---|---|---|---|---|
| 1 | 5261 | 1.5235 | 0.064844 | 20.282 | 20.606 | 0 | 0 |
| 2 | 1998VF31 | 1.5241 | 0.10031 | 31.298 | 31.77 | 0 | 0 |
| 3 | 1999ND43 | 1.5234 | 0.31433 | 5.5495 | 19.115 | 0 | 0 |
| 4 | 1999UJ7 | 1.5244 | 0.039262 | 16.749 | 16.895 | 0 | 0 |
| 5 | 2000XH47 | 1.5256 | 0.22535 | 16.944 | 21.254 | 0 | 0 |
| 6 | 2001DH47 | 1.5238 | 0.034643 | 24.401 | 24.476 | 0 | 0 |
| 7 | 2001FG24 | 1.5201 | 0.13248 | 19.195 | 20.597 | 0 | 0 |
| 8 | 2006XY2 | 1.5255 | 0.248 | 13.976 | 19.934 | 0 | 0 |
| 9 | 2007NS2 | 1.5237 | 0.05398 | 18.62 | 18.866 | 0 | 0 |
| 10 | 2007UR2 | 1.5254 | 0.061041 | 8.2637 | 8.9694 | 0 | 0 |
| 11 | 2009SE | 1.5244 | 0.065095 | 20.627 | 20.948 | 0 | 0 |
| 12 | 2010RL82 | 1.5211 | 0.099414 | 19.886 | 20.656 | 0 | 0 |
| 13 | 2010XH11 | 1.5256 | 0.24539 | 16.179 | 21.404 | 0 | 0 |
| 14 | 2011SC191 | 1.5238 | 0.044072 | 18.746 | 18.909 | 0 | 0 |
| 15 | 2011SL25 | 1.5238 | 0.11451 | 21.497 | 22.434 | 0 | 0 |
| 16 | 2011SP189 | 1.5237 | 0.040317 | 19.9 | 20.028 | 0 | 0 |
| 17 | 2011UB256 | 1.5236 | 0.070992 | 24.303 | 24.622 | 0 | 0 |
| 18 | 2011UN63 | 1.5237 | 0.064696 | 20.364 | 20.685 | 0 | 0 |
| 19 | 2012QR50 | 1.5263 | 0.17209 | 22.922 | 24.865 | 0 | 0 |
| 20 | 2013PM43 | 1.5276 | 0.13885 | 9.4212 | 12.324 | 0 | 0 |
| 21 | 2014HM187 | 1.5258 | 0.18972 | 0.60929 | 10.953 | 0 | 0 |
| 22 | 2014JU24 | 1.5212 | 0.23659 | 7.159 | 15.413 | 0 | 0 |
| 23 | 2014SA224 | 1.5217 | 0.28202 | 8.7836 | 18.53 | 0 | 0 |
| 24 | 2015CW12 | 1.5242 | 0.31403 | 5.8019 | 19.17 | 0 | 0 |
| 25 | 2015TL144 | 1.5239 | 0.078139 | 19.608 | 20.094 | 0 | 0 |
| 26 | 2016AA165 | 1.5229 | 0.089633 | 18.721 | 19.39 | 0 | 0 |
| 27 | 2016CP31 | 1.5236 | 0.058725 | 23.131 | 23.361 | 0 | 0 |
| 28 | 2016NZ55 | 1.525 | 0.25527 | 28.649 | 31.951 | 0 | 0 |
| 29 | 2016QY10 | 1.5228 | 0.25437 | 9.8374 | 17.657 | 0 | 0 |
| 30 | 2017QW35 | 1.5224 | 0.09739 | 28.295 | 28.797 | 0 | 0 |
| 31 | 2017XG62 | 1.5237 | 0.38103 | 42.148 | 46.727 | 30 | 30 |
| 32 | 2018EC4 | 1.5235 | 0.060547 | 21.837 | 22.098 | 0 | 0 |
| 33 | 2018FC4 | 1.5238 | 0.017115 | 22.146 | 22.166 | 0 | 0 |
| 34 | 2018FM29 | 1.5238 | 0.047251 | 21.5 | 21.662 | 0 | 0 |
| 35 | 2019BG3 | 1.5178 | 0.79844 | 9.4031 | 53.56 | 0 | 0 |
| 36 | 2019KF1 | 1.5223 | 0.17118 | 14.576 | 17.536 | 0 | 0 |
| 37 | 2020JO1 | 1.5245 | 0.29341 | 16.194 | 23.356 | 0 | 0 |
| 38 | 2020LE1 | 1.5239 | 0.6316 | 11.357 | 40.526 | 0 | 0 |
| 39 | 2020VT1 | 1.523 | 0.16698 | 18.719 | 20.964 | 0 | 0 |
| 40 | 2021FV1 | 1.5247 | 0.28428 | 3.7686 | 16.928 | 0 | 0 |
| 41 | 2021JK4 | 1.5256 | 0.18604 | 28.895 | 30.659 | 0 | 0 |
| 42 | 2021TB3 | 1.5256 | 0.23096 | 18.595 | 22.754 | 0 | 0 |
| 43 | 2021WX6 | 1.5233 | 0.37585 | 20.531 | 29.794 | 0 | 0 |
| 44 | 2022OG2 | 1.5218 | 0.21493 | 9.3525 | 15.496 | 0 | 0 |
| 45 | 2023FW14 | 1.5237 | 0.15809 | 13.273 | 16.044 | 0 | 0 |
| 46 | 2023QS3 | 1.525 | 0.32284 | 8.9637 | 20.789 | 0 | 0 |
| 47 | 2023RJ6 | 1.5245 | 0.30981 | 10.118 | 20.608 | 0 | 0 |
| 48 | 2023WV1 | 1.5238 | 0.29263 | 26.194 | 30.904 | 0 | 0 |

Observations on the tables (COMPUTED from the printed entries, no entry altered): Table 2 has 21 rows: P.L. 5 (3 objects), 10 (4), 15 (3),
20 (4), 25 (1), 30 (1), 35 (5); the F.L. column differs from P.L. in rows 5, 7 (to 5), 20 (35 to "40?"), and the F.L. 10 group is
2012 FC71 and 2021 VU12 only, matching the text of section 2.4. Table 3 has 16 rows (6 non-resonant Earth, 7 Earth ZLK, 3 Mars); the
six non-resonant Earth controls have Earth encounters 2008, 1638, 489, 944, 1545 and 183 (sum 6807), and all ZLK-state Earth
objects have 0. Table 5 and Table 2 use different label columns (Table 5 holds the S and R preliminary labels; Table 2 the clone-based
F.L.), so e.g. 2022 UP20 appears as 10/10 in Table 5 and F.L. 5 in Table 2, as printed.

## 4. Method in a few lines

READ/INFERRED. Integrate each known co-orbital for 78,000 years with a Bulirsch-Stoer N-body code (SWIFT-BS tolerance 1e-8, 1-day
step, all eight planets and the Moon, ecliptic elements at JD 2460800.5 rotated to the invariant plane) and a second integrator
(REBOUND BS), compute omega(t), assign Table 1 labels, then repeat with 729 clones from the 3 x 3 x 3 x 3 x 3 x 3 grid of
(a, e, i, Omega, omega, M) at nominal and plus or minus 1 sigma, summarise by <omega> and Amp., and finally integrate 10^5 years to
count Hill-sphere encounters. The semi-analytical content is only the averaged-Hamiltonian level curves of Michel and Thomas (1996)
for a non-resonant asteroid at a = 0.98 au (Fig. 1), used qualitatively. No co-orbital resonance Hamiltonian with a ZLK term is
derived or integrated in this paper. The cut that separates the four equilibrium regimes is i_max = 30 deg (eq. 1).

## 5. Positive controls for the project

There is none. The paper prints no trajectory, no periodic orbit, no transfer, no delta-v and no mass-ratio-model result that
`core/er3bp.py`, `core/wsb.py` or any project module computes. The only reproducible quantities are the Hill radii quoted on p15 (Venus
0.0067 au, Earth 0.0098 au, Moon 0.0004 au, Mars 0.0073 au) and the ZLK threshold i_max = 30 deg from eq. (1); both are textbook
closed forms, and the planetary masses of the Appendix are printed to 16 digits. The Earth Hill radius from eq. (2) with a_p = 1 au and
e_p = 0.0167 would be about 0.0098 au (COMPUTED, 1 au x 0.9833 x (3.0035e-6/(3 x 1.0000030))^(1/3) = 0.00983), a consistency check on
the printed number, not a project positive control.

## 6. What this means for the project

**Scope (stated plainly).** Out of scope for the catalogue. Nothing in the paper is a cycler, quasi_cycler, precursor_mga or mga_tour;
the objects are natural minor planets (Venus, Earth and Mars co-orbitals) and no spacecraft appears. It does not bear on
transfer design.

**The one project connection.** `docs/notes/2026-06-17-digest-fuente-marcos-2018.md` filed the Earth 1:1 co-orbital (Arjuna) population as
target data for a possible asteroid-leveraging cycler search (task `#308`, never reactivated) and as the dynamical-class definition of
a hypothetical Arjuna `quasi_cycler`. This paper adds, for that same population: (a) a larger and current list of 69 Earth co-orbitals
with a, e, inc and i_max (Table 5; also 22 Venus and 48 Mars co-orbitals, Tables 4 and 6) at the epoch the elements were drawn (JD
2460800.5 for the integrations; the tables themselves do not state their epoch, INFERRED to be the same); (b) a dynamical statement
that Earth co-orbitals librating about omega = 0 or 180 deg (2012 FC71, 2021 VU12, 2008 WM64, 2019 NC1, 2020 DX1, 2022 UO10) had no
close Earth encounters in 10^5 years, unlike low-i_max controls (183 to 2008 encounters); that is, for a target-leveraging scheme these
are the long-lived, predictable, low-encounter bodies among the Earth co-orbitals, while the ones the Arjuna digest called cheapest
(2006 RH120, 2009 BD) are not in this paper's Table 5 (INFERRED from their absence; not checked designation by designation). (c) The
paper's own caveat that many co-orbitals are chaotic on 150-year scales, so any use as a mission target is epoch-parametric. These
remain background for `#308` only; they do not change any catalogue row or ledger item and nothing needs to be done.

**Mechanism relevance.** The ZLK mechanism here is a secular coupling of e and inc in the heliocentric problem (planets as perturbers)
that matters at i_max above about 30 deg. The project's search and validation machinery is planar or low-inclination n-body and
restricted-three-body; none of the project's cycler families is analysed with secular averaged Hamiltonians, and the paper offers no
design tool. It is not relevant to quasi-cyclers or to ballistic transfer design.

**Possible side note, no action:** 469219 in Table 5 (row 21, a = 1.0011, e = 0.10385, inc = 7.7759 deg, i_max = 9.7866 deg, S = 25,
R = 0) is the minor-planet number of an Earth quasi-satellite that has been named as a spacecraft rendezvous target elsewhere (an
external fact, INFERRED from the number, not stated in the paper).

## 7. References cited that the project might want

Corpus status checked against `docs/notes/CORPUS_INDEX.md` on 2026-10-04 (rows grepped by author and title words). "Held" means a row
exists; "not held" means none.

| Reference | Held? | Why it might matter |
|---|---|---|
| de la Fuente Marcos, C., de la Fuente Marcos, R.: Meet Arjuna 2025 PN7, the newest quasi-satellite of Earth. Res. Notes Am. Astron. Soc. 9(9), 235 (2025) | not held | currency of the quasi-satellite list |
| de la Fuente Marcos, C., de la Fuente Marcos, R.: Dynamical evolution of near-earth asteroid 1991 VG. Mon. Not. R. Astron. Soc. 473(3), 2939-2948 (2018) | not held | Arjuna/minimoon dynamics |
| de la Fuente Marcos, C., de la Fuente Marcos, R.: Infrequent visitors of the Kozai kind: the dynamical lives of FC71, 2014 EK24, 2014 QD364, and 2014 UR? Astron. Astrophys. 580, A109 (2015) (the reference list prints 2012 as year) | not held | the earlier Earth-co-orbital Kozai paper |
| de la Fuente Marcos, C., de la Fuente Marcos, R.: From horseshoe to quasi-satellite and back again: the curious dynamics of earth co-orbital asteroid 2015 SO2. Astrophys. Space Sci. 361, 16 (2016) | not held | |
| de la Fuente Marcos, C., de la Fuente Marcos, R.: Asteroid 2014 OL339: yet another earth quasi-satellite. Mon. Not. R. Astron. Soc. 445(3), 2985-2994 (2014) | not held | |
| de la Fuente Marcos, C., de la Fuente Marcos, R.: A resonant family of dynamically cold small bodies in the near-earth asteroid belt. Mon. Not. R. Astron. Soc. 434, L1-L5 (2013) | not held | the original Arjuna-class definition paper; already on the acquisition wishlist in `2026-06-17-digest-fuente-marcos-2018.md` |
| Pan, T., Gallardo, T.: An attempt to build a dynamical catalog of present-day solar system co-orbitals. Celest. Mech. Dyn. Astron. 137(1), 2 (2025) | not held | the 54 plus 48 plus 20-object co-orbital catalogue and the "sticking" definition used throughout |
| Michel, P., Thomas, F.: The Kozai resonance for near-earth asteroids with semimajor axes smaller than 2 au. Astron. Astrophys. 307, 310-318 (1996) | not held | source of the Fig. 1 semi-analytical model |
| Nesvorný, D., Thomas, F., Ferraz-Mello, S., Morbidelli, A.: A perturbative treatment of the co-orbital motion. Celest. Mech. Dyn. Astron. 82(4), 323-361 (2002) | not held | analytical co-orbital model including ZLK |
| Thomas, F., Morbidelli, A.: The Kozai resonance in the outer solar system and the dynamics of long-period comets. Celest. Mech. Dyn. Astron. 64(3), 209-229 (1996) | not held | |
| Christou, A.A.: A numerical survey of transient co-orbitals of the terrestrial planets. Icarus 144(1), 1-20 (2000) | not held | |
| Morais, M.H.M., Morbidelli, A.: The population of near earth asteroids in coorbital motion with Venus. Icarus 185(1), 29-38 (2006) | not held | |
| Namouni, F.: Secular interactions of coorbiting objects. Icarus 137(2), 293-314 (1999) | not held | secular theory of co-orbitals |
| Sfair, R., Gomes, L.C., Winter, O.C., Moraes, R.A., Borderes-Motta, G., Schafer, C.M.: The Moon as a possible source for Earth's co-orbital bodies. Open J. Astrophys. 8 (September 2025), doi:10.33232/001c.144365 | not held | lunar-origin hypothesis for Earth co-orbitals |
| Carruba, V., Sfair, R., Winter, O.C., Mourao, D.C., Di Ruzza, S., Aljbaae, S., Carita, G., Domingos, R.C., Alves, A.A.: The invisible threat: assessing the collisional hazard posed by the undiscovered Venus co-orbital asteroids. Astron. Astrophys. 699, A86 (2025) | not held | |
| Carruba, V., Spoto, F., Barletta, W., Aljbaae, S., Martins, B.: The population of rotational fission clusters inside collisional asteroid families. Nat. Astron. 4, 83 (2020) | not held | source of the low-velocity encounter detection method |
| Levison, H.F., Duncan, M.J.: The long-term dynamical behavior of short-period comets. Icarus 108, 18-36 (1994); and SWIFT: a solar system integration software package, Astrophysics Source Code Library, ascl:1303.001 (2013) | not held | integrator references |
| Rein, H., Liu, S.-F.: REBOUND: an open-source multi-purpose N-body code for collisional dynamics. Astron. Astrophys. 537, A128 (2012) | not held | the project already uses REBOUND (see memory of the REBOUND variational-equation gotcha); the paper is the standard citation |
| Kozai, Y.: Secular perturbations of asteroids with high inclination and eccentricity. Astron. J. 67, 591-598 (1962); Kozai, Y.: Secular perturbations of resonant asteroids. Celest. Mech. Dyn. Astron. 36, 47-69 (1985); Lidov, M.L.: The evolution of orbits of artificial satellites of planets under the action of gravitational perturbations of external bodies. Planet. Space Sci. 9(10), 719-759 (1962); von Zeipel, H.: Astron. Nachr. 183, 345 (1910) | not held | the primary ZLK sources |

None of these is recommended for acquisition by this digest: the paper is out of scope for the catalogue. If task `#308` is ever
reactivated, Pan and Gallardo (2025) and the de la Fuente Marcos 2013 paper are the two to fetch first.

## 8. Minor observations (factual, for completeness)

- The abstract and Section 2 refer to "Michel and Thomas (1996)" for a model valid outside the 1:1 resonance and the authors state it
  is not adequate for quantitative analysis of co-orbitals (p4); the paper's conclusion on qualitative agreement is accordingly
  qualitative.
- The "Time domain astronomy; Time series analysis" keywords are the publisher's classification; they have no bearing on the content.
- The code is public (GitHub repository named on p21: `valeriocarruba/ZLK-mechanism-for-co-orbital-asteroids-of-terrestrial-planets-`)
  and the data are available "upon reasonable request" (p21).
- The text count of Venus co-orbitals (20, p2) and of the Table 4 rows (22) differ, and the text and Table 2 differ for 2002 AA29's label
  group on p10 and for the 2017 XG62 i_max (42.148 versus 46.727), as noted in sections 2.2 and 2.5; none changes a conclusion.
