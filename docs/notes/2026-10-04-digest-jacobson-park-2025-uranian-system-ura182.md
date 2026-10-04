# Digest: Jacobson & Park (2025), "The Orbits of Uranus, Its Satellites and Rings, the Gravity Field of the Uranian System, and the Orientation of the Poles of Uranus and Its Satellites" (the URA182 solution)

The Astronomical Journal 169:65 (17pp), 2025 February, DOI 10.3847/1538-3881/ad99d1 (open access; received 2024 August 15, accepted
2024 November 24, published 2025 January 13). Authors Robert A. Jacobson and Ryan S. Park, Jet Propulsion Laboratory.
Filed in the private paper corpus as
jacobson-park-2025-orbits-uranus-satellites-rings-gravity-field-poles-URA182-aj-169-65-doi-10.3847-1538-3881-ad99d1.pdf

Digested 2026-10-04 from all 17 pages (page images read; table digits cross-checked against the PDF text layer). Each statement is marked
READ (seen on the page; page, section, equation or table given), COMPUTED (our arithmetic, shown) or INFERRED (our reading or a comparison
with project code). Page numbers are the journal's printed page numbers. This note extends, and does not repeat, the digest of the
previous solution: `docs/notes/2026-10-04-digest-jacobson-2014-uranian-satellites-gravity-field.md` (URA111), cited below as "the 2014 digest".
Why the project holds this paper: tasks #890, #894 and #895 propagate Uranian moons against the URA111 kernel with constants whose sources
were unknown; this paper is the URA182 successor and the registry's semi-major axes turn out to be its values.

## 0. Headline findings for the project

1. READ (Table 6, p7) and COMPUTED: the registry's five Uranian semi-major axes (`core/satellites.py`: Miranda 129846, Ariel 190929,
   Umbriel 265986, Titania 436298, Oberon 583511 km) are digit for digit this paper's Table 6 "Satellite Equatorial Geometric Orbital
   Elements" values. Difference zero for all five. The 2014 paper's Table 2 mean axes differ from them by -12, -1, +4, +16, +62 km
   (the 2014 digest's section 6). The registry's moon GMs and its Uranus system GM are digit for digit Jacobson 2014 Table 12 (the 2014
   digest, headline 1), not this paper's Table 2. COMPUTED: the registry therefore mixes URA182 geometry with URA111 masses. The
   coordinator's reading is confirmed digit by digit (section 4).
2. READ (Table 2, p3): URA182 GMs are Titania 222.80 +/- 1.45, Oberon 214.21 +/- 2.14, Ariel 83.43 +/- 0.56, Umbriel 85.40 +/- 0.68,
   Miranda 4.11 +/- 0.11, system 5794560.69 +/- 2.18, Uranus 5793950.61 +/- 2.16 km^3 s^-2. COMPUTED: Titania changed by -4.10 (-1.81
   percent), Oberon by +8.91 (+4.34 percent) relative to URA111: the coordinator's arithmetic is correct (section 7).
3. READ (Table 6, p7) and COMPUTED: Kepler's third law on a geometric axis does not reproduce the printed period. With the registry's
   model (system GM minus the other moon, URA111 GMs) the periods are 8.706399 d (Titania) and 13.465984 d (Oberon) against
   the printed 8.705869 d and 13.463237 d: too long by 6.1e-5 and 2.0e-4 (relative). Planet oblateness plus the inner moons acting as
   a ring, evaluated with the paper's own "augmented J2" idea, account for only 8 to 14 percent of the gap (section 5). The paper does
   not say what accounts for the rest and prints no mean-motion relation (INFERRED reading in section 5).
4. READ (Table 6 vs the 2014 digest's Table 2) and COMPUTED: the six-decimal periods are unchanged between the two solutions
   (Titania 8.705869 vs 360/41.3514187 = 8.705868; Oberon 13.463237 both; Ariel, Umbriel, Miranda likewise to the sixth decimal). The
   mean-motion information in the two solutions is the same at this precision; the semi-major axes are not (they were redefined as geometric elements).
5. READ (Table 12, p10): the prime-meridian rates give ten-digit mean rates: Ariel 142.83570520, Umbriel 86.86887774, Titania
   41.35141483, Oberon 26.73948040, Miranda 254.69071249, Puck 472.54457439 deg/day. These are the best-printed rates in the paper
   (section 3.4 and 5).
6. READ (p4, section 3.1; Table 3): URA182 zonal harmonics at the 25,559 km reference radius are J2 = 3508.967 (COR) or 3509.709 (COO)
   x 10^-6, J4 = -35.793 or -35.044 x 10^-6, J6 = 0.575 x 10^-6 (set from Neuenschwander and Helled 2022). The project's J2 of 3509.291e-6
   is the French et al. 2024 value, which is neither URA111's (3510.7) nor URA182's. The satellite orbits were integrated with the COR set.
7. READ (p15, concluding remarks): "The planet ephemeris designation is DE442 and that of the satellite ephemeris is URA182." (the abstract and section 1
   say DE440 was updated for this work; the paper does not reconcile DE440 with DE442).

## 1. What the paper does

READ (abstract and section 1, p1): an update of Jacobson (2014). It redetermines "the Uranus and satellite masses, satellite orbits, and
the orbit of Uranus", obtains "a value for the Uranus tidal dissipation factor", produces "an independent determination of the Uranian
ring orbits, Uranus pole direction, and Uranus gravity harmonics" and "new expressions for the orientations of the satellites". "We
processed all of the data from our previous work plus the new and rereduced astrometry and the ring occultations. We extended our data
arc forward to 2016 and backward to 1847. Our new orbit for Uranus is an update of that in ephemeris DE440 incorporating the occultation
data and Gaia astrometry."

Motivation (p1): new astrometry since 2014; some of the original astrometry rereduced against Gaia DR3; and French et al. (2024)
redetermined ring orbits, pole and harmonics from the full 1977 to 2006 ring occultation set. The paper updates the planetary
ephemeris to DE440 (Park et al. 2021) to support the satellite ephemerides.

Data (section 2.2, p2-3; Tables 13 to 18, p11 to p15):
- Earth-based and HST planet and satellite astrometry (including the visual micrometer observations before 1911 omitted in 2014; Table 14
  begins with Herschel 1787 to 1832); Earth-based transits; mutual events; satellite stellar occultations.
- Uranian ring stellar occultations from Earth and from the Voyager 2 photopolarimeter; Voyager 2 radio occultations; the set was extended
  to match French et al. (2024), 1977 to 2006. Eight of the nine main rings (the gamma ring omitted); 566 observations remain after
  dropping 83 (section 3.5, p8; Figure 7, p15).
- Voyager 2 Doppler, range, delta-DOR and optical navigation imaging.
- New: Sheshan and Yunnan CCD and photographic astrometry; Pulkovo CCD; the 29-year LNA/MCTI set rereduced against Gaia DR3 (Camargo et
  al. 2022); Gaia astrometry of the satellites themselves (Tanga et al. 2023, Ariel, Umbriel, Titania, Oberon, 2014 August to 2016
  December); HST positions from the MAST archive (all five major satellites and Puck); Oberon Carlsberg Meridian Circle positions
  1992 to 1995; 58 astrometric positions of Uranus from ring occultation geometry.
- READ (p13): "we extended the observational data forward to the most recent astrometry available and back to the earliest astrometry of
  Ariel and Umbriel (we omitted earlier observations of Titania and Oberon as they are of insufficient accuracy to warrant the
  computational effort required to process them)."

Dynamical model (section 2.1, p1-2). READ:
- "The dynamical model for the motion of Uranus is the same that used in the development of JPL planetary ephemeris DE440 (R. S. Park et
  al. 2021). C. F. Peters (1981) provides the basic model for the satellite dynamics that was used by R. A. Jacobson (2014)."
- New terms: the zero GM of Puck replaced with a value derived from its size and an assumed density (1 g cm^-3, p3); the dissipation
  effect of the tides raised on Uranus by the satellites; "a post-Newtonian general relativistic correction, and the Lense-Thirring effect".
- The solar GM is augmented with the GMs of the inner planets, the Earth and the Moon ("bodies with orbits interior to Jupiter's orbit
  have a negligible direct effect on the orbits of outer planet satellites, while adding their GMs to that of the Sun approximates their
  average effect and simplifies the orbit integration"). Table 1 (p2): Jovian system GM 126,712,761.1898 (JUP387), Saturnian system GM
  37,940,584.9205 (Jacobson 2022), Neptunian system GM 6,836,527.1006 (Jacobson 2009), augmented Sun GM 132,713,233,433.8835 (Park et al. 2021), all km^3 s^-2.
- Uranus gravity "is represented by the second, fourth, and sixth zonal harmonics of its gravitational potential, but instead of the
  previous value of zero for the sixth harmonic, we adopted the value from B. A. Neuenschwander & R. Helled (2022)."
- The Uranus pole is "obtained from the numerical integration of its rotational equations with torques from the Sun and the five major
  satellites"; the old rotational model of Jacobson (2014) was replaced by the more general one used for Saturn (Jacobson 2022).
- Rings: the processing precessing-Keplerian-ellipse model is replaced by the extended model of French et al. (2024), which allows normal
  modes of radial distortion (Table 9).
- READ (p2, section 2.1): the Voyager 2 spacecraft dynamics are unchanged. The paper does not mention ring mass as a force term (the 2014 paper ignored it; INFERRED that this still holds).
  (The integration frame, ICRF centred on the Uranian system barycentre, is stated in the 2014 paper and not repeated here; INFERRED unchanged.)

Estimated parameters (p2): epoch state of the Uranian system barycentre; epoch state of each satellite; ring elements and normal modes; the
GMs of the Uranian system and satellites; J2 and J4 of Uranus (J6 fixed from Neuenschwander and Helled 2022, with 20 percent uncertainty
treated as consider, p4); pole RA and Dec at the integration epoch; the tidal quality factor Q. Spacecraft parameters (epoch state,
manoeuvres, nongravitational accelerations, RTG thermal radiation corrections).

Fit method (section 2.3, p3): weighted least squares in square-root-information form (Bierman 1977), singular-value decomposition of the
combined array. Twelve observation-model parameter sets listed (observer biases, Sytinskaja position-angle bias, absolute-astrometry
opposition biases, Descamps orientation and scale, HST phase-angle biases, double star sigma Sgr corrections, stellar proper-motion
corrections and timing offsets, Voyager camera sample/line biases, imaging phase-angle biases, camera pointing, range biases,
ionosphere and troposphere delay corrections).

Formal-error scale (p3-4): "the computation of our current formal errors includes a scale factor (Lawson & Hanson 1974) zeta^2 = sos/(m - n)",
sum of squares, number of data points, number of estimated parameters; "For our final solution is zeta = 0.4652. This scale factor accounts
for most of the reduction in our GM uncertainties, although the larger astrometric data set contributes somewhat to that reduction."
(INFERRED: the quoted formal uncertainties in Tables 2, 3, 5 are scaled down by 0.4652; see section 7.)

Ephemerides produced (p15): "Ephemerides for Uranus and the satellites are available electronically from the On-Line Solar System Data
Service at the Jet Propulsion Laboratory (http://ssd.jpl.nasa.gov) and from NASA's Navigation and Ancillary Information Facility
(http://naif.jpl.nasa.gov). The planet ephemeris designation is DE442 and that of the satellite ephemeris is URA182."
Time span: the satellite-orbit fit interval for the mean elements and tidal fit is 1100 yr (1550 to 2650, p7); the pole integration
was over 1000 yr for the pole series and 300 centuries for the verification (p4); the rotational element fit is over 1900 to 2100
(p11); the data arc is 1847 to 2016 (abstract). Orbit uncertainty statement: the 2025 to 2075 period (Table 7).

Accuracy claims (p7-8): Table 7 uncertainties "are larger than the formal ones quoted by R. A. Jacobson (2014) and represent a more
accurate appraisal of the orbit uncertainty"; Uranus orbit "uncertainty is optimistic because it assumes the data errors are well
characterized and not systematic"; "we recommend adopting a realistic uncertainty equal to triple the formal uncertainty" for the Uranus orbit (p8).

Other statements of note:
- READ (p3): "There are no resonances in the Uranian satellite system to enhance the determination of the satellite GMs ... The near-commensurability
  among Miranda, Ariel, and Umbriel (Section 3.4) does help by adding a constraint on the product of the Ariel and Umbriel GMs. During
  the Voyager 2-Uranus encounter, continuous tracking data were acquired only during the Uranus and Miranda close flybys. Consequently,
  the GMs of only those two bodies are determined directly from the Doppler data."
- READ (p3): "Our current analysis produced no significant change in the Miranda, Ariel, and Umbriel GMs, but that Titania decreased whereas
  Oberon increased. We attribute the changes in the latter two GMs primarily to the rereduced and additional astrometry."
- READ (p14): "That astrometry also introduced small but significant changes in the orbital periods of the satellites, and the revised
  Uranus pole direction altered the orientation of the satellites' orbital planes."
- READ (p7): the tidal acceleration of Ariel and Umbriel from Emelyanov and Nikonchuk (2013) is near this paper's; Titania, Oberon, Miranda are not.

## 2. Tidal dissipation (section 3.3, p5-7; Table 5; Equation 1)

READ: tidal quality factor Q from the phase lag epsilon = 2 (Delta t / r0) |omega p x r0 - r0-dot|, Q = cot epsilon. Love number k2 = 0.300 adopted
(same as Nimmo 2023). The authors could not obtain independent lags per satellite; they imposed the same Q for all satellites and found
Q = 678 +/- 231 (p6). Kaula (1964) approximation (Equation 1, p6): n-dot/n = -(9/2) (k2/Q) (m/m0) (R0/a)^5 n, with n the satellite mean motion,
m its mass, m0 and R0 the planet's mass and radius; with k2/Q = 4.421936e-4 (estimated GMs, and a and n from Table 6).

### Table 5 (p7): Tidal accelerations

| Satellite | approximate n-dot/n (yr^-1) | estimated n-dot/n (yr^-1) | a-dot (cm yr^-1) | Emelyanov & Nikonchuk (2013) n-dot/n (yr^-1) |
| --- | --- | --- | --- | --- |
| Miranda | -6.76e-10 | (-6.77 +/- 2.29)e-10 | 5.86 +/- 1.98 | -1.24e-8 |
| Ariel | -1.12e-9 | (-1.12 +/- 0.38)e-9 | 14.27 +/- 4.86 | -1.89e-9 |
| Umbriel | -1.33e-10 | (-1.33 +/- 0.45)e-10 | 2.37 +/- 0.80 | -4.12e-10 |
| Titania | -1.39e-11 | (-1.35 +/- 0.47)e-11 | 0.39 +/- 0.14 | -2.93e-9 |
| Oberon | -2.02e-12 | (-3.60 +/- 0.69)e-12 | 0.14 +/- 0.03 | -5.07e-9 |

Note (READ): "The quoted uncertainties are the formal 1 sigma uncertainties." The second column is obtained by fitting a quadratic
polynomial to the angular coordinate over the 1100 yr integration (p6-7). READ: the approximate expression "provides a good match to our
estimated values for all but Oberon".
INFERRED relevance: tidal secular acceleration of Titania is -1.35e-11 per year of n; over a century that is a fractional mean-motion
change of 1.4e-9, negligible for the project's circular two-moon model.

## 3. Definitions that matter for using the numbers

### 3.1 Geometric elements, epoch, reference plane (section 3.4, p7; Table 6 caption)

READ (p7): "The precise orbits of the satellites are computed using numerical integration. The orbits can be roughly represented by
geometric elements (S. Renner & B. Sicardy 2006) referred to the Uranian equatorial plane. The elements for the current orbits appear in
Table 6: a SMAA, e eccentricity, and i inclination to the equator. The table also provides the orbital period and the precession periods
of the longitude of periapsis and the node. The elements are the average values of the osculating elements and longitude rates over an 1100 yr
integration (1550-2650)." (The sentence breaks across the two columns of p7; the second half is in the right column.)

- Epoch: none. READ: Table 6 has no epoch row and no longitude columns; the elements are averages over 1550 to 2650, not values at an epoch.
  (Contrast the 2014 digest's Table 2, which was "Mean Equatorial Orbital Elements at 2000 January 1.5" with a longitude epoch.)
- Reference plane: the Uranian equatorial plane (READ). The pole orientation used is the paper's own (Table 4, section 3.2), not the 2014
  pole, so inclinations differ from 2014 even where the orbits are unchanged (INFERRED: e.g. Titania i 0.114 deg now vs 0.1129 then;
  Ariel 0.026 vs 0.0167).
- Geometric vs osculating vs mean Keplerian: the paper gives no definition of its own and refers to Renner and Sicardy (2006)
  (not held). The paper's phrase "geometric elements ... roughly represent the orbits" and "average values of the osculating elements"
  is all it says. INFERRED from the title and use of that reference: geometric elements describe the shape of the actual precessing
  orbit in the oblate field (radius r = a(1 - e cos M) with an epicyclic anomaly) rather than the osculating Kepler ellipse of the
  instantaneous state, so the geometric a is a mean physical orbit size, not a quantity tied to the mean motion by n^2 a^3 = GM. This is
  an inference from the reference title ("Use of geometric elements in numerical simulations"), not a statement in the paper. The ring
  elements of Tables 8 and 9 are also geometric (the paper calls the ring orbits "geometric ring midlines", p4).

### 3.2 The mean-motion relation: the paper prints none

READ: no equation in the paper relates the mean motion to the semi-major axis, the GM and the oblateness. A search of all 17 pages found no
epicyclic-frequency formula and no statement of "n^2 a^3 = GM". The nearest content is the pole-precession expression of Ward (1975) and
Ward and Hamilton (2004) on p4, unnumbered:

    psi-dot = -(3/2) (n_sun^2 / s) [(J2 + q)/(gamma + l)] cos(epsilon),
    q = (1/2) sum_k (mu_k/mu_0) (a_k/R)^2,     l = sum_k (n_k/s) (mu_k/mu_0) (a_k/R)^2,

where n_sun is the planet's heliocentric mean motion, s its rotation rate, J2 its degree-2 zonal harmonic, gamma its normalised polar
moment of inertia, R its equatorial radius, mu_0 and mu_k G times the planet's and satellite k's mass, n_k and a_k satellite k's mean
motion and semi-major axis. READ (p4): "However, we assumed an 'augmented' J2 = J2 + q and an 'augmented' gamma = gamma + l." That is the
paper's statement that equatorial satellites act as a ring adding q to J2. With the URA182 values (COMPUTED: q = 0.5 * sum (GM_k / GM_U)(a_k/R)^2 over the five moons with GM_U = 5793950.61, R = 25559): q = 0.016447, which is 4.69 times
J2 (0.003509). Used to compute the pole precession, not the mean motions; the same ring idea is applied to the moons' mean motions in section 5.
READ (p4): the test integration "applying only the direct solar torque on Uranus", with the augmented J2 and gamma, "obtained a value of
psi-dot in close agreement with the theoretical prediction"; the Nettelmann et al. (2013) parameters give psi-dot = 0.00649 arcsec per
year, "yielding a precession period of 200 million years".

The only mean-motion statement is Equation (1), the tidal n-dot/n, and the commensurability relations on p7:
    lambda-dot_5 - 3 lambda-dot_1 + 2 lambda-dot_2 = -0.0785 deg/day   (5 Miranda, 1 Ariel, 2 Umbriel)
    lambda-dot_1 - lambda-dot_2 - 2 lambda-dot_3 + lambda-dot_4 = 0.0038 deg/day   (3 Titania, 4 Oberon; "not associated with a resonance; it
    is merely an interesting distribution of orbital periods"; the paper also cites Harris 1949).
COMPUTED from Table 6 periods (rate = 360/P): the first is 254.69073 - 3(142.83566) + 2(86.86888) = -0.07850 deg/day (matches -0.0785). The
second is 142.83566 - 86.86888 - 2(41.35141) + 26.73948 = +0.00344 deg/day, which does not match the printed 0.0038 (difference 0.00036 deg/day,
about 10 percent). The period digits are not the cause (their rounding changes the sum by at most about 1e-6). Unresolved; the printed 0.0038
may come from unrounded rates or be a slip. It does not affect the project's use of the rates.

### 3.3 Frame and pole (section 3.2, p4-6; Table 4; Figures 1 to 3)

READ (p4-5): the pole is represented by a six-term Fourier series fitted to the numerically integrated rigid-body pole under torques from
the Sun and five satellites (integration over 1000 yr with J2, J4, J6 and gamma = 0.2224). "The first five periodic terms in the series are
related to the nodal precessions of the Miranda, Ariel, Umbriel, Titania, and Oberon orbits. The period of the sixth term is half the
orbital period of Uranus." The series "match the integration with an error of about 3 microdegrees over the 1000 yr." They could not
estimate the precession or gamma from the data ("The 14 yr of ring occultations dominate the determination of the pole orientation").
Final J2000 pole (not the leading series terms), formal 1-sigma (p5): alpha = 77.311258 +/- 0.000091 deg, delta = 15.171977 +/- 0.000511 deg,
correlation -0.5293; error ellipse SMAA 0.00051332 deg, SMIA 0.00007657 deg, rotation angle 95.486510 deg (Figure 3, p6, shows the
URA182, URA111 and French et al. poles on 1986 January 19 12:00:00; the URA111 pole lies well outside the URA182 and French ellipses).
IAU relation (READ p5): alpha_IAU = alpha + 180 deg, delta_IAU = -delta, W_IAU = 180 deg - W. The IAU "north pole" of Uranus is the pole
of rotation on the north side of the invariable plane; this paper's pole is the spin angular momentum. At J2000 the prime meridian W = 156.19 deg
(p5), in agreement with the IAU definition.
Two conventions: the satellites' prime meridians are measured "from the node of the satellite's equator on the ICRF reference plane" (p10).

## 4. Constants audit (extending the 2014 digest's)

Values in km^3 s^-2 for GM, km for lengths. "Reg." is the project's `core/satellites.py` / `v4_uranus.py` value. "2014" is Jacobson 2014
(Table 12 for GM and J, Table 2 for mean elements, as transcribed in the 2014 digest; the 2014 Uranus-only GM is from this paper's Table 2).
"URA182" is this paper (Table 2 for GM, Table 6 for a, Table 3 for J, with the J2 reference radius 25,559 km).

### 4.1 Masses

| Quantity | Registry | 2014 (URA111) | URA182 (Table 2) | Reg. - 2014 | Reg. - URA182 | Registry's source |
| --- | --- | --- | --- | --- | --- | --- |
| Uranus system GM | 5794556.4 (`PRIMARIES["Uranus"]` = 5.7945564e6) | 5794556.4 +/- 4.3 | 5794560.69 +/- 2.18 | 0.0 | -4.29 | the 2014 value (comment says "JPL DE440 planetary constants, astro_par") |
| Uranus alone GM | not carried | 5793951.3 +/- 4.4 | 5793950.61 +/- 2.16 | | | |
| Miranda GM | 4.3 | 4.3 +/- 0.2 | 4.11 +/- 0.11 | 0.0 | +0.19 | 2014 |
| Ariel GM | 83.5 | 83.5 +/- 1.4 | 83.43 +/- 0.56 | 0.0 | +0.07 | 2014 |
| Umbriel GM | 85.1 | 85.1 +/- 1.9 | 85.40 +/- 0.68 | 0.0 | -0.30 | 2014 |
| Titania GM | 226.9 | 226.9 +/- 4.1 | 222.80 +/- 1.45 | 0.0 | +4.10 | 2014 |
| Oberon GM | 205.3 | 205.3 +/- 5.8 | 214.21 +/- 2.14 | 0.0 | -8.91 | 2014 |
| Puck GM | not carried | 0 (massless) | 0.1275 +/- 0.0425 | | | |

COMPUTED check of the paper's own bookkeeping (Uranus alone = system - sum of satellite GMs): URA182: 5794560.69 - (83.43 + 85.40 + 222.80 +
214.21 + 4.11 + 0.1275) = 5794560.69 - 610.0775 = 5793950.6125, equals the printed 5793950.61. URA111: 5794556.4 - (83.5 + 85.1 + 226.9 +
205.3 + 4.3) = 5794556.4 - 605.1 = 5793951.3, equals the printed 5793951.3 (the 2014 paper's Table 12 printed only the system GM; the Uranus-alone figure is this paper's Table 2 middle column).
Total moon mass: URA111 605.1, URA182 610.08 including Puck (609.95 without).

### 4.2 Orbital elements

| Moon | Reg. a | URA182 Table 6 a | Reg. - URA182 | 2014 Table 2 mean a | Reg. - 2014 |
| --- | --- | --- | --- | --- | --- |
| Miranda | 129846 | 129,846 | 0 | 129858 | -12 |
| Ariel | 190929 | 190,929 | 0 | 190930 | -1 |
| Umbriel | 265986 | 265,986 | 0 | 265982 | +4 |
| Titania | 436298 | 436,298 | 0 | 436282 | +16 |
| Oberon | 583511 | 583,511 | 0 | 583449 | +62 |

Registry mean motions are not carried as data: `SatelliteData.mean_motion_deg_day` is derived at import by Kepler's third law from `sma_km` and the
primary's GM (`mean_motion_deg_day_about`), and `TwoMoonModel.n_base` uses `sqrt((gm_planet + gm_base)/length_km^3)` with
`gm_planet = gm_sys - gm_base - gm_pert`.

### 4.3 Gravity field

| Quantity | Project | 2014 | URA182 COR | URA182 COO | French et al. 2024 | Project - 2014 | Project - COR | Project - COO |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| J2 (x10^6) | 3509.291 (`URANUS_J2`) | 3510.7 +/- 0.7 | 3508.967 +/- 0.064 | 3509.709 +/- 0.064 | 3509.291 +/- 0.412 | -1.409 | +0.324 | -0.418 |
| J4 (x10^6) | not carried | -34.2 +/- 1.3 | -35.793 +/- 0.140 | -35.044 +/- 0.140 | -35.522 +/- 0.466 | | | |
| J6 (x10^6) | not carried | 0.0 +/- 1.0 (set to zero) | 0.575 +/- 0.115 | 0.575 +/- 0.115 | 0.500 | | | |
| reference radius | 25559.0 (`URANUS_R_EQ_KM`) | 25559 | 25,559 | 25,559 | (see Table 3 note) | 0 | 0 | 0 |

READ (Table 3 note, p4): "The reference radius for the zonal harmonics is 25,559 km." Footnotes: (a) formal 1 sigma uncertainties; (b) the French et
al. column carries realistic 1 sigma uncertainties; (c) J6 from Neuenschwander and Helled (2022), uncertainty 20 percent. The project's J2 equals the French
et al. (2024) value digit for digit; the coordinator's code comment already says so. COMPUTED relative to this paper's formal 0.064: the project is
0.324/0.064 = 5.1 formal sigma from COR and 6.5 from COO, but 0.8 of French's realistic 0.412 from COR.

### 4.4 Radii and flyby floors

The moon mean radii (Miranda 235.8, Ariel 578.9, Umbriel 584.7, Titania 788.9, Oberon 761.4 km) and the flyby altitude floors (100 km Miranda, 50 km
others) are not in this paper (READ: no radii or shapes appear in any of the 17 pages, apart from the rings). Source unknown after this paper; the
registry comment attributes the radii to JPL SSD phys_par and the floors to Heaton-Longuski 2003 Table 4.

### 4.5 Verdict on the mixed-source reading

CONFIRMED digit by digit. Semi-major axes: five of five equal this paper's Table 6 (differences 0, 0, 0, 0, 0) and none equals the 2014 Table 2 (differences -12 to +62
km). GMs: six of six (system plus five moons) equal Jacobson 2014 Table 12 "Current Results" and none equals this paper's Table 2 (differences up to 8.91 km^3 s^-2
for Oberon and 4.29 for the system). INFERRED explanation: the registry's comments say the GMs are from JPL SSD "phys_par, ref URA111" and the elements from SSD's
mean-elements page (accessed 2026-06-14); so the SSD satellite physical-parameter page evidently still lists the URA111 masses while its mean-elements page lists the
URA182 geometric elements. That reading cannot be checked offline and the registry should record the two sources separately.

## 5. Mean motions for a circular two-moon model

All arithmetic is COMPUTED with `python3` (Kepler's third law, P = 2 pi sqrt(a^3 / GM)).

### 5.1 Kepler's third law against the printed period

Printed (Table 6): Titania 8.705869 d, Oberon 13.463237 d.

(a) Registry GM and a, in the model's own form (central GM = system GM minus the other moon; Titania: 5794556.4 - 205.3 = 5794351.1, a = 436298;
    Oberon: 5794556.4 - 226.9 = 5794329.5, a = 583511):

| Moon | P (a) | P(a) - printed (days) | relative |
| --- | --- | --- | --- |
| Titania | 8.706399 | +0.000530 | +6.09e-5 |
| Oberon | 13.465984 | +0.002747 | +2.04e-4 |

(b) URA182 GMs with the same form (Titania: 5794560.69 - 214.21 = 5794346.48; Oberon: 5794560.69 - 222.80 = 5794337.89), same a:

| Moon | P (b) | P(b) - printed | relative |
| --- | --- | --- | --- |
| Titania | 8.706402 | +0.000533 | +6.13e-5 |
| Oberon | 13.465975 | +0.002738 | +2.03e-4 |

(The system-GM-only variants give 8.706245 / 8.706241 (Titania, 2014 / 2025 system GM) and 13.465721 / 13.465716 (Oberon); same gap.)
Changing from URA111 to URA182 masses moves the Kepler periods by only 3e-6 and 9e-6 d: the mass change is not the issue. The gap is in the semi-major
axis (a shift of 1e-4 in a gives 1.5e-4 in P). The registry's synodic period from (a) is 1/(1/8.706399 - 1/13.465984) = 24.632445 d, against
1/(1/8.705869 - 1/13.463237) = 24.637400 d from the printed periods: 0.00496 d (7.1 minutes) per synodic period, 2.0e-4 relative.

(c) The paper prints no oblateness-corrected mean-motion relation (section 3.2). INFERRED standard circular-orbit relation for the radial force on a
moon at radius a in the field of an oblate planet plus interior moons treated as rings plus exterior moons treated as rings (the same ring idea as
the paper's "augmented J2", with q R^2 GM_U = (1/2) sum GM_k a_k^2):

    n^2 a^3 = [GM_U + sum_(interior) GM_k + GM_m] + (3/2)[J2 R^2 GM_U + (1/2) sum_(interior) GM_k a_k^2] / a^2 - (1/2) sum_(exterior) GM_k (a/a_k)^3

(the last term is the outward in-plane pull of an exterior ring, whose radial acceleration at r << a_k is +G M r/(2 a_k^3); the first bracket is the
monopole, with the moon's own GM added for relative motion; J4 and higher are negligible (J4 (R/a)^4 is about 1e-9 of the monopole). Inputs: URA182
GMs, GM_U = 5793950.61, R = 25559, J2 = 3508.967e-6 (COR; COO changes the result by 1e-7), geometric a for all moons.)

| Moon | monopole only (U + interior + self) P (d) | + J2 + interior rings + exterior rings P (d) | printed P | remaining gap (relative) |
| --- | --- | --- | --- | --- |
| Titania | 8.706402 | 8.706330 | 8.705869 | +5.3e-5 |
| Oberon | 13.465716 | 13.465516 | 13.463237 | +1.69e-4 |

The J2 + ring terms shorten Titania's period by 7.2e-5 d (8.3e-6 relative; the interior rings are 34 percent of the planet's quadrupole at Titania
and 194 percent at Oberon) and Oberon's by 2.0e-4 d (1.5e-5 relative). Of the gap between the monopole value and the printed period (5.3e-4 d and 2.48e-3 d),
this explains 14 percent for Titania and 8 percent for Oberon. The remainder does not come from planet oblateness or from the inner moons acting as a ring.

Which effect accounts for the rest? The paper does not let one tell (READ: no discussion). Two COMPUTED diagnostics narrow it:
1. Using the 2014 Table 2 mean axes (436282, 583449 km) in the same full relation, the model-versus-printed residual in n^2 is only -4e-6 (Titania) and +2e-5 (Oberon),
   and +6e-6 (Ariel), +9e-6 (Umbriel). With this paper's geometric axes the residual is +1.06e-4 (Titania), +3.39e-4 (Oberon), -9e-6 (Ariel), +5.4e-5 (Umbriel). So the 2014
   "mean" axes are consistent with the printed rates through Kepler's law plus oblateness; this paper's geometric axes are not, by up to 3e-4. (Miranda is the exception: 2014
   mean a gives +2.3e-4, the geometric a -5e-5, which may reflect its near-resonant perturbation; INFERRED.)
2. The geometric axes that would reproduce the printed periods in the full relation are 436282.6 km (Titania, 15.4 km below the printed 436298) and 583445.2 km (Oberon,
   65.8 km below 583511); with plain Kepler on the system GM (5794560.69) they are 436285.6 km (12.4 km below) and 583439.4 km (71.6 km below). The full-relation values are within 4 km of the 2014 mean axes (436282, 583449), the plain-Kepler values within 4 and 10 km.
INFERRED conclusion: the geometric a of Table 6 is a different average from the mean-motion axis (it is averaged over 1100 yr of osculating-derived geometric elements, which with
eccentricities of 0.0016 and the mutual perturbations of Titania, Oberon, Umbriel and Ariel need not satisfy n^2 a^3 = GM), and the difference is 12 to 16 km for Titania and 66 to 72 km for Oberon. The
two effects named in the request (planet oblateness; the inner moons acting as a ring) are real but contribute at most about a seventh of the gap. The paper does not state the
cause and this note does not claim one.

Scale for the model: the orbital eccentricity amplitudes e a are 698 km (Titania, e = 0.0016) and 934 km (Oberon, e = 0.0016) (Umbriel 1037 km, Ariel 267 km, Miranda 182 km), an
order of magnitude larger than the 12 to 72 km axis ambiguity. A circular model cannot represent either better than the eccentricity allows.

### 5.2 Recommendation (INFERRED)

For a circular two-moon model in which both the timing (periods, hence the synodic phase) and the distance (hence the flyby geometry) matter:
1. Take the angular rates from the printed rates, not from Kepler's law: Table 12 prime-meridian rates (ten digits) Titania 41.35141483 and Oberon 26.73948040 deg/day, or
   Table 6's periods 8.705869 d and 13.463237 d (six decimals). The two sources differ by 1.4e-6 d for Oberon (13.463238 d from the W rate vs 13.463237 d printed) and by 1e-6 d for
   Titania; over a century that phase difference is about 0.1 degree for Oberon. The W rate being the orbital angular rate is INFERRED from synchronous rotation (p10); it agrees with
   the 2014 mean longitude rates to 1e-7 relative (Titania 41.3514187 vs 41.35141483; Oberon 26.7394835 vs 26.73948040).
2. Take the distances from the geometric Table 6 axes, which are what the published orbits and the URA182 kernel reproduce on average (INFERRED; the paper says "average values" over 1100 yr).
3. Do not force n^2 a^3 = GM for the moons. If the model's nondimensionalisation (`TwoMoonModel`, `CCR4BPSystem`) requires it, the effective central mass for the frame is
   GM_eff = n^2 a^3 = 5795056.3 km^3 s^-2 (Titania) and 5796694.5 (Oberon), which exceeds the system GM by 8.6e-5 and 3.7e-4: that single number cannot serve both moons (they differ by
   2.8e-4 relative), so the circular two-moon CCR4BP, with its single planet GM, cannot reproduce both moons' periods and distances at once. The faithful alternative is to prescribe
   each moon kinematically as a circle of radius a_geo at its printed rate (as an ephemeris) while the spacecraft feels GM_U plus the moon point masses; the moons need not be Kepler
   solutions of the model's own attraction. The spacecraft's own Kepler problem then uses GM_U (5793950.61) or the system GM, whichever the model's accounting dictates. (Size of
   the choice: GM_U vs system GM is 1.05e-4.)
4. If Kepler-derived rates must be kept, use the axes that reproduce the periods (Titania 436285.6 km, Oberon 583439.4 km with the system GM), accepting a 12 and 72 km distance error,
   both far below the eccentricity amplitudes.

## 6. Which gravity field goes with which ephemeris

READ (p4, section 3.1): the harmonics are determined chiefly from ring occultation precession rates. "We initially determined our zonal gravity harmonics from the precession rates of the
geometric ring midlines (COR). The radii of those midlines are the radii in the ring orbit model that are estimated in the fit to the observations (see Section 3.5). We subsequently estimated the
gravity harmonics that produce the precession rates at the COO using the ring radius corrections from Table 16 of R. G. French et al. (2024). Table 3 provides the zonal gravity harmonics. Our
uncertainties are the formal uncertainties produced by the data fit." "The values from French et al. (2024) fall between our estimates found using the COR and those found using the COO. We made
no attempt to quantify our actual harmonic uncertainties, but as the determination of the gravity harmonics is dominated by the ring occultations, we expect that our knowledge of them is no
better than that postulated by French et al. (2024)." (COR is the centre of the geometric ring midline; COO the ring "center of opacity".)

READ (p4): "The satellite orbits and Voyager 2 trajectory are integrated with the harmonics found from the COR precession rates, i.e., those found from the fit to the observations. The integrations
are actually insensitive to the differences between the COR and COO harmonics." So the URA182 satellite orbits were produced with the COR set: J2 = 3508.967, J4 = -35.793, J6 = 0.575 (x10^-6) at 25,559 km.
The difference between the COR and COO sets (J2: 0.742e-6, J4: 0.749e-6) is the quantified inconsistency tolerance the authors accept. French et al. adopted the COO-based set with realistic uncertainties.

How the URA182 values differ from URA111's (COMPUTED from Table 3): J2 is 1.733e-6 lower (COR) or 0.991e-6 lower (COO): 2.5 and 1.4 of the 2014 formal sigma 0.7 (and 27 formal 2025 sigmas for COR); J4 is 1.593e-6 lower (COR) or 0.844e-6 lower
(COO), 1.2 and 0.6 of the 2014 sigma 1.3; J6 changes from fixed 0 to 0.575. READ (p4): J6 "is indeterminate with the current set of occultation data"; the adopted value is a
literature value (Neuenschwander and Helled 2022) with 20 percent assumed uncertainty, treated as consider.

INFERRED consistent sets (an ephemeris and the gravity field it was fitted with):
- URA111 (the kernel the project propagates against, `v4_uranus_strict.py`): the 2014 paper's Table 12 "Current Results": system GM 5794556.4 (Uranus alone 5793951.3), moon GMs 4.3, 83.5, 85.1,
  226.9, 205.3 (Miranda, Ariel, Umbriel, Titania, Oberon), Puck massless (zero), J2 = 3510.7e-6, J4 = -34.2e-6, J6 = 0, all at 25,559 km. (The 2014 digest established these.) The mean elements and rates of the
  2014 paper's Table 2 (a 129858, 190930, 265982, 436282, 583449 km; rates in the 2014 digest) go with this solution, since its mean axes satisfy Kepler plus oblateness (section 5.1).
- URA182: GMs of this paper's Table 2 (system 5794560.69, Uranus 5793950.61, Miranda 4.11, Ariel 83.43, Umbriel 85.40, Titania 222.80, Oberon 214.21, Puck 0.1275) with the COR harmonics J2 = 3508.967e-6, J4 =
  -35.793e-6, J6 = 0.575e-6 at 25,559 km; the COO set (3509.709, -35.044, 0.575) is the within-model alternative the authors say the orbits are insensitive to. French et al.'s (3509.291, -35.522, 0.500) lies between the
  two but is not the set the URA182 orbits were integrated with.
- The project's present mixture (URA111 GMs, URA182 geometric axes, French J2) is none of these. INFERRED consequence: with J2 = 3509.291e-6 propagating against URA111, the J2 difference from the URA111 value
  is 1.4e-6 (0.04 percent of J2, 2.0 sigma of the 2014 value), a small effect; the larger inconsistency is the semi-major axes and the Kepler-derived rates (section 5).
- One more gap in consistency: URA111 treats Puck as massless; URA182 gives it 0.1275 km^3 s^-2, which changes nothing at the project's level.

## 7. How much the moons' masses changed between the solutions

COMPUTED from Table 2 (arithmetic: change = new - old; percent = 100 * change / old; sigma = change / formal sigma):

| Body | URA111 GM +/- 1 sigma | URA182 GM +/- 1 sigma | change | percent | in 2014 sigma | in 2025 sigma | in quadrature sigma |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Miranda | 4.3 +/- 0.2 | 4.11 +/- 0.11 | -0.19 | -4.42 | -0.95 | -1.7 | -0.83 |
| Ariel | 83.5 +/- 1.4 | 83.43 +/- 0.56 | -0.07 | -0.08 | -0.05 | -0.1 | -0.05 |
| Umbriel | 85.1 +/- 1.9 | 85.40 +/- 0.68 | +0.30 | +0.35 | +0.16 | +0.4 | +0.15 |
| Titania | 226.9 +/- 4.1 | 222.80 +/- 1.45 | -4.10 | -1.81 | -1.00 | -2.8 | -0.94 |
| Oberon | 205.3 +/- 5.8 | 214.21 +/- 2.14 | +8.91 | +4.34 | +1.54 | +4.2 | +1.44 |
| system | 5794556.4 +/- 4.3 | 5794560.69 +/- 2.18 | +4.29 | +7.4e-5 | +1.0 | +2.0 | +0.89 |

The coordinator's -1.8 percent and +4.3 percent are confirmed (-1.81 and +4.34). The two Titania/Oberon changes are inside the 2014 uncertainties (1.0 and 1.5 sigma) but 2.8 and 4.2 sigma in this paper's
own scaled formal errors. The 2025 formal errors are scaled by zeta = 0.4652 (section 1); the unscaled equivalents are about 1.45/0.4652 = 3.1 (Titania) and 2.14/0.4652 = 4.6 (Oberon), under which the changes
are 1.3 and 1.9 sigma (INFERRED: assumes the scale factor multiplies the formal 1-sigma uncertainties uniformly, which the paper states for the GM uncertainties only). READ: the authors attribute the Titania and Oberon changes "primarily to the rereduced and additional astrometry"
(p3). The sum of Titania and Oberon changes is +4.81 (COMPUTED: -4.10 + 8.91), so most of the Oberon increase is not compensated by Titania; the sum of all five moon GMs changed by +4.85
(COMPUTED: 83.43 + 85.40 + 222.80 + 214.21 + 4.11 = 609.95 against 605.1; the system GM increased by 4.29 and Uranus alone decreased by 0.69).
Mass ratio Oberon/Titania: 205.3/226.9 = 0.905 (URA111) vs 214.21/222.80 = 0.961 (URA182). For the project's two-moon flyby work the two moons' masses are the load-bearing numbers: a
flyby bend angle scales with the moon's GM, so a 1.8 percent smaller Titania and a 4.3 percent larger Oberon change their gravity-assist capability by those fractions (INFERRED, to first order at fixed periapsis and v-infinity).
The paper states no preferred mass for the project's use; the URA182 values carry the smaller formal uncertainty and the more complete data.

## 8. Tables, transcribed (the observation-residual tables summarised)

Tables 1, 2, 3, 5 and the Uranus Table 4 values appear in sections 1 to 6 above. The remaining tables follow. Table numbers are this paper's (not the 2014 paper's).

### Table 1 (p2): Gravitational constants (km^3 s^-2)
Jovian system GM 126,712,761.1898 (JUP387); Saturnian system GM 37,940,584.9205 (R. A. Jacobson 2022); Neptunian system GM 6,836,527.1006 (R. A. Jacobson 2009); Sun GM (augmented) 132,713,233,433.8835 (R. S. Park et al. 2021).

### Table 2 (p3): Uranus gravity parameters, GM (km^3 s^-2), formal 1-sigma

| Body | R. A. Jacobson et al. (1992) | R. A. Jacobson (2014) | Current |
| --- | --- | --- | --- |
| System | 5794548.6 +/- 1.5 | 5794556.4 +/- 4.3 | 5794560.69 +/- 2.18 |
| Uranus | 5793939.3 +/- 2.8 | 5793951.3 +/- 4.4 | 5793950.61 +/- 2.16 |
| Ariel | 90.3 +/- 2.3 | 83.5 +/- 1.4 | 83.43 +/- 0.56 |
| Umbriel | 78.2 +/- 2.3 | 85.1 +/- 1.9 | 85.40 +/- 0.68 |
| Titania | 235.3 +/- 1.9 | 226.9 +/- 4.1 | 222.80 +/- 1.45 |
| Oberon | 201.1 +/- 1.8 | 205.3 +/- 5.8 | 214.21 +/- 2.14 |
| Miranda | 4.4 +/- 0.2 | 4.3 +/- 0.2 | 4.11 +/- 0.11 |
| Puck | | | 0.1275 +/- 0.0425 |

### Table 3 (p4): Uranus zonal gravity harmonics (x10^6), reference radius 25,559 km

| Parameter | R. A. Jacobson (2014) (a) | Current COR (a) | Current COO (a) | R. G. French et al. (2024) (b) |
| --- | --- | --- | --- | --- |
| J2 | 3510.7 +/- 0.7 | 3508.967 +/- 0.064 | 3509.709 +/- 0.064 | 3509.291 +/- 0.412 |
| J4 | -34.2 +/- 1.3 | -35.793 +/- 0.140 | -35.044 +/- 0.140 | -35.522 +/- 0.466 |
| J6 | 0.0 +/- 1.0 | 0.575 +/- 0.115 (c) | 0.575 +/- 0.115 (c) | 0.500 |

Notes (READ): (a) formal 1 sigma; (b) realistic 1 sigma; (c) from Neuenschwander and Helled (2022). The French J6 carries no uncertainty; they assumed 0.50e-6 (p4).

### Table 4 (p6): Uranus body orientation angles (degrees; d = days from J2000, T = Julian centuries from J2000)

    alpha = 77.311864 + 0.0001858313 T - 0.000176 sin U1 + 0.000009 sin U2 - 0.000090 sin U3 + 0.000903 sin U4 + 0.000098 sin U5 + 0.000002 sin U6
    delta = 15.171402 + 0.0000182771 T - 0.000170 cos U1 + 0.000009 cos U2 - 0.000087 cos U3 + 0.000871 cos U4 + 0.000094 cos U5 - 0.000014 cos U6
    W     = 156.189841 + 501.1600928000 d + 0.000046 sin U1 - 0.000002 sin U2 + 0.000024 sin U3 - 0.000236 sin U4 - 0.000026 sin U5

    U1 = 79.140037 + 2023.9699523 T     U2 = 290.329944 + 623.6558178 T     U3 = 6.577578 + 276.1324960 T
    U4 = 325.874083 + 26.3248132 T      U5 = 104.911743 + 185.3598278 T    U6 = 293.114711 + 856.9614776 T

(The signs of the sine and cosine coefficients are as in the page image and the PDF text layer. The constant terms differ from the final estimated pole of p5 because the series terms at J2000 do not sum to zero.)

### Table 5 (p7): see section 2.

### Table 6 (p7): Satellite equatorial geometric orbital elements (averages over 1550 to 2650; no epoch)

| Element | Ariel | Umbriel | Titania | Oberon | Miranda | Puck |
| --- | --- | --- | --- | --- | --- | --- |
| a (km) | 190,929 | 265,986 | 436,298 | 583,511 | 129,846 | 86,004 |
| e | 0.0014 | 0.0039 | 0.0016 | 0.0016 | 0.0014 | 0.0001 |
| i (deg) | 0.026 | 0.083 | 0.114 | 0.125 | 4.421 | 0.327 |
| P (days) | 2.520379 | 4.144177 | 8.705869 | 13.463237 | 1.413479 | 0.761833 |
| P_varpi (yr) | 57.833 | 126.794 | 895.800 | 894.289 | 17.969 | 4.452 |
| P_Omega (yr) | 57.770 | 129.945 | 1644.649 | 192.798 | 17.787 | 4.454 |

(Uncertainties are not printed in this table. The 2014 digest's Table 2 gave Puck a = 86005, e = 0.00019, i = 0.3562, and the rates in deg/day or deg/yr; COMPUTED apsidal and nodal periods from those rates, 360/rate in years: varpi Ariel 57.78,
Umbriel 126.64, Titania 360.25, Oberon 1343.28, Miranda 17.96; node 59.17, 136.80, 17308, 855.11 and 17.78. The apsidal periods of Ariel, Umbriel and Miranda agree with Table 6 to 0.1 percent; their nodal periods differ by 2.4, 5 and 0.0 percent. Titania's apsidal and nodal periods
differ from Table 6 by factors of 2.5 and 10.5 and Oberon's by 1.5 and 4.4, which is not a printing issue: their eccentricity and inclination are small and their apsides and nodes poorly defined. They have no bearing on the project.)

### Table 7 (p7): Satellite orbit uncertainties for 2025 to 2075 (km; T-dot in km yr^-1). R radial, T tangential, N normal; "formal ... larger than ... Jacobson (2014)".

| Satellite | R | T | N | T-dot | Satellite | R | T | N | T-dot |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ariel | 19 | 600 | 60 | 6.2 | Oberon | 44 | 200 | 69 | 1.6 |
| Umbriel | 29 | 200 | 60 | 1.7 | Miranda | 11 | 800 | 290 | 8.7 |
| Titania | 44 | 200 | 68 | 1.5 | Puck | 11 | 175 | 31 | 1.6 |

READ (p7-8): "The rather large uncertainties in the Ariel and Miranda T positions and growth rates reflect their higher sensitivity to possible errors in the tidal forces (Section 3.3). Miranda's large N uncertainty reflects its significant orbital inclination."

### Table 8 (p9): Ring centerline elements referred to the Uranus equator; epoch for the longitudes 1986 Jan 19 12:00:00 (TDT)
Angles in degrees, rates in deg/day, e multiplied by 10^3, formal 1-sigma.

| Element | Ring 6 | Ring 5 | Ring 4 |
| --- | --- | --- | --- |
| a (km) | 41,837.1026 +/- 0.0924 | 42,234.9324 +/- 0.0894 | 42,571.1339 +/- 0.0895 |
| e (x10^3) | 1.01532 +/- 0.00122 | 1.89922 +/- 0.00113 | 1.06437 +/- 0.00118 |
| varpi (deg) | 181.692972 +/- 0.072831 | 176.829324 +/- 0.031647 | 255.806004 +/- 0.057418 |
| i (deg) | 0.061264 +/- 0.000079 | 0.055570 +/- 0.000085 | 0.031643 +/- 0.000060 |
| Omega (deg) | 89.955906 +/- 0.125412 | 296.293251 +/- 0.125345 | 335.952351 +/- 0.209739 |
| varpi-dot (deg/day) | 2.76202261 +/- 0.00001340 | 2.67157993 +/- 0.00000895 | 2.59807961 +/- 0.00001060 |
| Omega-dot (deg/day) | -2.75656742 +/- 0.00001336 | -2.66640225 +/- 0.00000892 | -2.59312168 +/- 0.00001056 |

| Element | alpha | beta | eta |
| --- | --- | --- | --- |
| a (km) | 44,718.3831 +/- 0.0852 | 45,661.0117 +/- 0.0846 | 47,176.0119 +/- 0.0844 |
| e (x10^3) | 0.75832 +/- 0.00092 | 0.43997 +/- 0.00103 | |
| varpi (deg) | 206.165484 +/- 0.070977 | 318.051627 +/- 0.119635 | |
| i (deg) | 0.015639 +/- 0.000073 | 0.004934 +/- 0.000055 | |
| Omega (deg) | 204.216686 +/- 0.208135 | 227.971322 +/- 0.759842 | |
| varpi-dot (deg/day) | 2.18521243 +/- 0.00000732 | 2.03070087 +/- 0.00000676 | |
| Omega-dot (deg/day) | -2.18143486 +/- 0.00000728 | -2.02733449 +/- 0.00000673 | |

| Element | delta | epsilon |
| --- | --- | --- |
| a (km) | 48,300.2297 +/- 0.0809 | 51,149.2877 +/- 0.0778 |
| e (x10^3) | | 7.93512 +/- 0.00083 |
| varpi (deg) | | 307.074205 +/- 0.006027 |
| varpi-dot (deg/day) | | 1.36326174 +/- 0.00000293 |

(Blank cells are blank in the paper: eta, delta, epsilon rings have no entries for the elements not estimated.)

### Table 9 (p9): Ring normal modes

| | Mode | Value | | Mode | Value |
| --- | --- | --- | --- | --- | --- |
| delta ring amplitude (km) | 02 | 3.1889038 +/- 0.0402299 | epsilon ring amplitude (km) | 24 | 0.4378138 +/- 0.0428276 |
| phase (deg) | 02 | 170.2540431 +/- 0.4507573 | phase (deg) | 24 | 2.7556718 +/- 0.2382426 |
| frequency (deg/day) | 02 | 562.5163165 +/- 0.0013998 | frequency (deg/day) | 24 | 1074.5227896 +/- 0.0000973 |
| amplitude (km) | 23 | 0.3641671 +/- 0.0425702 | amplitude (km) | 14 | 0.3827226 +/- 0.0490925 |
| phase (deg) | 23 | 6.7289722 +/- 0.2977699 | phase (deg) | 14 | 13.9984600 +/- 0.4313380 |
| frequency (deg/day) | 23 | 1074.5230423 +/- 0.0001222 | frequency (deg/day) | 14 | 956.4179480 +/- 0.0001676 |
| eta ring amplitude (km) | 03 | 0.5914081 +/- 0.0466169 | | | |
| phase (deg) | 03 | 75.7593278 +/- 1.6156016 | | | |
| frequency (deg/day) | 03 | 776.5841246 +/- 0.0010459 | | | |

### Table 10 (p9): Ring centerline radius comparison (km)

| Ring | R. G. French et al. (2024) | URA182 | Diff. |
| --- | --- | --- | --- |
| 6 | 41,837.092 +/- 0.096 | 41,837.103 +/- 0.092 | 0.011 |
| 5 | 42,234.893 +/- 0.091 | 42,234.932 +/- 0.089 | 0.039 |
| 4 | 42,571.124 +/- 0.091 | 42,571.134 +/- 0.090 | 0.010 |
| alpha | 44,718.473 +/- 0.086 | 44,718.383 +/- 0.085 | -0.090 |
| beta | 45,661.056 +/- 0.087 | 45,661.012 +/- 0.085 | -0.044 |
| eta | 47,176.009 +/- 0.088 | 47,176.012 +/- 0.084 | 0.003 |
| delta | 48,300.227 +/- 0.082 | 48,300.230 +/- 0.081 | 0.003 |
| epsilon | 51,149.279 +/- 0.081 | 51,149.289 +/- 0.078 | 0.010 |

READ (p8): "For the most part agreement is well within 1 sigma. The exception is the alpha ring, however, there is agreement within the overlapping formal uncertainties."

### Table 11 (p10): Voyager 2 navigation performance, Uranus B-plane

| B.R (km) | B.T (km) | SMAA (km) | SMIA (km) | theta (deg) | TCA | sigma TCA (s) | Source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 25,384 | 128,680 | 1 | 1 | 70 | 1986 Jan 24 17:59:46.5 | 0.07 | A. H. Taylor et al. (1986) |
| 25,384.9 | 128,678.8 | 0.9 | 0.1 | 101 | 1986 Jan 24 17:59:46.544 | 0.004 | R. A. Jacobson & B. Rush (2007) |
| 25,385.5 | 128,678.5 | 0.6 | 0.05 | 101 | 1986 Jan 24 17:59:46.531 | 0.002 | R. A. Jacobson (2014) |
| 25,386.1 | 128,678.4 | 0.89 | 0.14 | 99 | 1986 Jan 24 17:59:46.533 | 0.004 | Current |

### Table 12 (p10): Satellite body orientation angles (degrees; d = days from J2000, T = centuries from J2000)
READ from the page and the text layer; the interleaving of coefficient signs was checked against the page image.

    Ariel   alpha = 77.312 - 0.013 sin U1 - 0.013 sin U2 + 0.020 sin U3 - 0.006 sin U4 - 0.010 sin U5 - 0.002 sin U6
            delta = 15.169 - 0.012 cos U1 - 0.012 cos U2 + 0.019 cos U3 - 0.007 cos U4 - 0.009 cos U5 - 0.003 cos U6
            W     = 20.210 + 142.83570520 d + 0.003 sin U1 + 0.004 sin U2 + 0.008 sin U3 - 5.187 sin U4 - 0.001 sin U5 - 0.003 sin U6
                    - 0.001 sin U9 - 0.012 sin U10 - 0.098 sin U11 + 0.006 sin U12
    Umbriel alpha = 77.313 - 0.001 sin U1 + 0.003 sin U2 + 0.081 sin U3 - 0.023 sin U4 - 0.037 sin U5 - 0.001 sin U6
            delta = 15.162 - 0.001 cos U1 + 0.002 cos U2 + 0.078 cos U3 - 0.022 cos U4 - 0.035 cos U5 - 0.003 cos U6
            W     = 71.044 + 86.86887774 d - 0.001 sin U2 - 0.020 sin U3 - 0.359 sin U4 + 0.010 sin U5 - 0.006 sin U6 + 0.004 sin U10 + 0.033 sin U11
    Titania alpha = 77.316 - 0.010 sin U3 - 0.088 sin U4 - 0.077 sin U5 - 0.001 sin U6
            delta = 15.140 - 0.010 cos U3 - 0.085 cos U4 - 0.074 cos U5 - 0.006 cos U6
            W     = 101.566 + 41.35141483 d + 0.003 sin U3 + 0.031 sin U4 + 0.021 sin U5 - 0.012 sin U6
    Oberon  alpha = 77.316 + 0.002 sin U3 - 0.106 sin U4 + 0.060 sin U5 - 0.002 sin U6
            delta = 15.132 + 0.002 cos U3 - 0.103 cos U4 + 0.058 cos U5 - 0.009 cos U6
            W     = 172.662 + 26.73948040 d - 0.001 sin U3 + 0.216 sin U4 - 0.013 sin U5 - 0.018 sin U6
    Miranda alpha = 77.312 + 4.590 sin U1 - 0.001 sin U2 + 0.002 sin U3 - 0.001 sin U5 + 0.048 sin U7 - 0.003 sin U8
            delta = 15.148 + 4.426 cos U1 - 0.001 cos U2 + 0.002 cos U3 - 0.001 cos U5 + 0.023 cos U7 - 0.001 cos U8
            W     = 145.922 + 254.69071249 d - 1.200 sin U1 + 0.005 sin U2 + 0.013 sin U3 - 4.959 sin U4 - 0.022 sin U5 - 0.002 sin U6 - 0.098 sin U7
                    + 0.002 sin U8 + 0.018 sin U9 + 0.176 sin U10 + 1.435 sin U11
    Puck    alpha = 77.312 + 0.006 sin U1 + 0.001 sin U4 - 0.339 sin U13
            delta = 15.171 + 0.005 cos U1 + 0.001 cos U4 - 0.327 cos U13
            W     = 86.237 + 472.54457439 d - 0.001 sin U1 + 0.002 sin U4 - 0.001 sin U6 + 0.089 sin U13

    U1-U6 as in Table 4;
    U7  = 158.107946 + 4047.4585218 T     U11 = 317.981514 + 2868.9479866 T
    U8  = 55.379520 + 6063.0529503 T      U12 = 173.777592 + 428.5000000 T
    U9  = 53.633870 + 8606.1928945 T      U13 = 143.956927 + 8082.8458396 T
    U10 = 276.017756 + 5737.9905163 T

(Miranda's last three W terms (U9, U10, U11) and Ariel's U9 to U12 terms were read from the page image because the text layer drops lines there; unclear digits: none.) READ (p11): the series were fitted to the satellite orbits over 1900 to 2100; "Our series differ significantly from the currently
available IAU series (B. A. Archinal et al. 2018) that are based on the Voyager 2 pre-Uranus-encounter satellite ephemerides. We recommend that our series be considered as a possible replacement." U1 to U6 are the nodal precessions of the five satellite orbits and half Uranus's orbital period; U7 is twice
U1 (approximately; the printed rates differ), U8 triple U1 (approximately), U9 related to a Miranda-Puck interaction, U11 related to the near-commensurability of Miranda, Ariel and Umbriel, U10 about twice U11, U13 Puck's nodal precession rate, the source of U12 unknown (all READ, p11).
The table says the satellites are "assumed to be in synchronous rotation, where their prime meridians are defined to point toward Uranus at their periapsis and apoapsis" (READ p10; so W increases at roughly the orbital mean rate).

### Tables 13 to 18 (p11 to p15): observation residuals, summarised (not transcribed)
- Table 13 (p11): satellite transits (Carlsberg 1999, 79 each of alpha and delta, rms 0.351 and 0.299 arcsec, 1992 to 1995; Arlot et al. 2008, 227 each, 0.153 and 0.185, 1997 to 2005), mutual events (3 rows: Christou 2009, Mallama 2009, Arlot 2013), Titania and Umbriel stellar occultation astrometry (Herald et al. 2020; 2 and 1 observations). 7 rows.
- Table 14 (p12-13): residual statistics of every Earth-based astrometric data set of the major satellites: type, number, rms (arcsec), dates, site, source; filar micrometer 1787 to 1949, photographic 1913 to 1990, CCD 1981 to 2016, with Gaia (757 along-scan, along-cross) in 2014 to 2016; about 74 rows (counted from the text layer, approximate). Footnotes a (not used in the fit: Herschel 1835, Lamont 1840) and b (used to update both the satellite and Uranus orbits: Camargo 2022, Stone and Harris 2000, Stone 2000, 2001, 2005, Owen 2009, Monet 2007, Harris 2016, Zhang 2022, Xie 2019).
- Table 15 (p14): Uranus observation residuals (transit, CCD, stellar occultation 58 positions from French 2024 with rms 0.005 and 0.004 arcsec, stellar relative astrometry); 23 rows (counted from the page image).
- Table 16 (p14): Puck observation residuals; 5 rows (Pascu 1998, Descamps 2002, Showalter 2008, Veiga and Bourget 2006, Showalter 2017).
- Table 17 (p15): Voyager imaging residuals in pixels (sample, line): Miranda 116 (0.249, 0.229), Ariel 109 (0.399, 0.317), Umbriel 103 (0.265, 0.331), Titania 64 (0.244, 0.341), Oberon 64 (0.167, 0.346), Puck 49 (0.189, 0.136), stars 1406 (0.269, 0.272). 7 rows. (Transcribed because short.)
- Table 18 (p15): ring radius residuals by event (km), 32 events (Voyager RSS, sigma Sgr, beta Per and Earth-based events U0 to U0602); rms between 0.1324 km (U17B) and 0.5180 km (U144). READ: 566 ring observations, rms 0.3045 km (p11; Figure 7, p15).
- Figures 1 and 2 (p5): pole right ascension (77.31 to 77.315 deg) and declination (15.17 to 15.172 deg) over 1600 to 2600; Figure 3 (p6): pole error ellipses (URA182, URA111, French); Figure 4 (p8): Uranus orbit formal uncertainties 1960 to 2060 (R.A. and decl. up to 100 and 40 mas, range up to about 300 km); Figure 5 (p14): Voyager Doppler residuals during the encounter (1070 points, rms 0.351437 mm/s); Figure 6 (p14): optical residuals (1010 points, rms 0.289510 pixels); Figure 7 (p15): ring radius residuals (566 points, rms 0.304517 km).

## 9. References the project might want

Corpus status checked against `docs/notes/CORPUS_INDEX.md` on 2026-10-04 (only Jacobson 2014 and Stone and Miner 1986 of the Uranus literature are indexed; none of the others is).

| Reference | Why it matters | In corpus |
| --- | --- | --- |
| Jacobson, R. A. 2014, AJ 148:76, DOI 10.1088/0004-6256/148/5/76 | URA111 solution; the registry's GMs | yes (digest 2026-10-04) |
| French, R. G., McGhee-French, C. A., Gordon, M. K., Nicholson, P. D., Longaretti, P.-Y., et al. 2024, Icarus 411, 115957 (arXiv:2401.04634) "The Uranus system from occultation observations (1977-2006): Rings, pole direction, gravity field, and masses of Cressida, Cordelia, and Ophelia" | the project's J2 source; Table 16 (COO corrections); realistic harmonic uncertainties | not in CORPUS_INDEX as of this digest; an uncommitted digest `2026-10-04-digest-french-2024-uranus-system-occultations-gravity-field.md` was present in the working tree (another agent's) |
| French, R. G., McGhee-French, C. A., Nicholson, P. D., Longaretti, P.-Y., et al. 2023, Icarus 395, 115474 (ring occultation data) | the occultation data description | no |
| Renner, S. & Sicardy, B. 2006, Celest. Mech. Dyn. Astron. 94, 237 "Use of the geometric elements in numerical simulations" | the definition of the geometric elements of Table 6 and the relation to the mean motion; read this to resolve the axis question of section 5 | no |
| Neuenschwander, B. A. & Helled, R. 2022, MNRAS 512, 3124 | the J6 value adopted | no |
| Park, R. S., Folkner, W. M., Williams, J. G., Boggs, D. H. 2021, AJ 161, 105 (DE440) | planetary ephemeris and GMs | no |
| Archinal, B. A., et al. 2018, Celest. Mech. Dyn. Astron. 130, 22 (IAU rotational elements) | the IAU series this paper proposes to replace | no |
| Jacobson, R. A. et al. 1992, AJ 103, 2068; Jacobson & Rush 2007 (AAS/AIAA Astrodynamics Specialist Conf.) | the earlier GM solution; the Voyager encounter reconstruction | no |
| Ward, W. R. & Hamilton, D. P. 2004, AJ 128, 2501; Ward 1975, AJ 80, 64 | the "augmented J2" precession expression | no |
| Peters, C. F. 1981 (basis of the satellite equations of motion; full citation in the paper's reference list, not transcribed here) | equations of motion | no |
| Greenberg, R. 1975, MNRAS 173, 121; 1976, Icarus 29, 427 | the Miranda-Ariel-Umbriel near-resonance | no |
| Laskar, J. & Jacobson, R. A. 1987, A&A 188, 212 (GUST86) | analytic satellite theory | no |
| Emelyanov, N. V. & Nikonchuk, D. V. 2013, MNRAS 436, 3668 | independent satellite ephemerides and tidal accelerations | no |
| Tanga, P., Pauwels, T., Mignard, F., et al. 2023, A&A 674, A12 | Gaia astrometry of the satellites | no |
| Baland, R.-M., Filice, V., Le Maistre, S., et al. 2025, Icarus 426, 116371 | satellite orientation and interior models | no |
| Nimmo, F. 2023, PSJ 4, 241; Cuk, M., El Moutamid, M. & Tiscareno, M. S. 2020, PSJ 1, 22 | Uranus tidal dissipation | no |
| Heaton, A. F. & Longuski, J. M. 2003 (flyby floors) | the registry's 50 km floor | yes (heaton-longuski-2003 file is indexed) |

## 10. Still unknown after this paper

READ and COMPUTED limits of this paper for the project:
1. Moon mean radii (235.8, 578.9, 584.7, 788.9, 761.4 km) and flyby floors (100, 50, 50, 50, 50 km): not in this paper (no satellite radii or shapes appear in it). Source unknown here; left to the
   registry comments (JPL SSD phys_par; Heaton-Longuski 2003 Table 4).
2. The definition of the geometric semi-major axis and why it departs from the mean-motion axis by 12 to 16 km (Titania) and 66 to 72 km (Oberon): not stated; Renner and Sicardy (2006), not held, would settle it.
3. The Uranus system and moon GMs on SSD's pages: this paper does not say which solution SSD lists; the inference in section 4.5 is unverified.
4. Whether the printed 0.0038 deg/day Titania-Oberon commensurability is a slip (computed 0.00344 from Table 6).
5. The project's Uranus J2 and J4 for URA182 and URA111: sourced here (section 6), but which set the project should adopt is a decision for the coordinator.
