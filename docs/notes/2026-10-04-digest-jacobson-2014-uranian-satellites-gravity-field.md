# Digest: Jacobson (2014), "The orbits of the Uranian satellites and rings, the gravity field of the Uranian system, and the orientation of the pole of Uranus"

The Astronomical Journal 148:76 (2014 November), DOI 10.1088/0004-6256/148/5/76, 13 pages (author affiliation: Jet Propulsion
Laboratory; received 2014 March 14, accepted 2014 July 2, published 2014 September 19).
Filed in the private paper corpus as
jacobson-2014-orbits-uranian-satellites-rings-gravity-field-uranian-system-pole-aj-148-76-doi-10.1088-0004-6256-148-5-76.pdf

Digested 2026-10-04 from all 13 pages. Each statement is marked READ (seen on the page; page and section or table given),
COMPUTED (our arithmetic, shown) or INFERRED (our reading or a comparison with project code). Page numbers are the journal's
printed page numbers. Tables were transcribed from the page images; digits that were hard to resolve are marked "unclear".
Why the project holds this paper: it was cited, unread, for the Uranus J2, the reference radius and (by inference) the
Uranian GM values; tasks #890, #894, #895 depend on them.

## 0. Headline findings for the project

1. READ (Table 12, p9): the project's Uranus system GM (5,794,556.4) and all five moon GMs (Miranda 4.3, Ariel 83.5, Umbriel
   85.1, Titania 226.9, Oberon 205.3 km^3 s^-2) are digit for digit the "Current Results" column of Table 12. COMPUTED: difference zero.
2. READ (Table 12 notes, p9): "The reference radius for the Uranus zonal harmonics is 25559 km." The project's
   `URANUS_R_EQ_KM = 25559.0` and its stated source (Jacobson 2014) are correct on that point. The paper calls it the
   "reference radius for the zonal harmonics", not the physical equatorial radius (INFERRED: the two are the same number only by convention).
3. READ (Table 12, p9): Current Results J2 = 3510.7 +/- 0.7 (x 10^-6), J4 = -34.2 +/- 1.3 (x 10^-6), J6 = 0.0 +/- 1.0 (x 10^-6, "Not
   estimated"). The project's `URANUS_J2 = 3509.291e-6` (French et al. 2024) differs by 1.409e-6, about 0.04 percent, 2.0 sigma of Jacobson's
   own 0.7 (COMPUTED). The older 3.34343e-3 that the project carried until 2026-10-04 does not appear anywhere in this paper (READ: Table 12 is
   the only table with J2 values; the paper's "Table 4" is the Titania occultation star catalogue, p6, not a gravity table).
4. READ (Table 2, p5): the paper's mean semi-major axes differ from the project's `SATELLITES` values by 1 to 62 km (section 6). The
   project's values do not come from this paper; they are the JPL Satellite Mean Elements values by the registry comments.
5. READ (Table 2, p5) and COMPUTED: the paper's mean longitude rates give Titania 8.705868 d and Oberon 13.463237 d, exactly the URA111
   fitted periods the #890 review reported. The project should use these (section 7).
6. READ (p10): "The ephemeris designation is URA111." This paper is the source of the URA111 satellite ephemeris.

## 1. What the paper does

READ (abstract, p1; section 1, p1-2): a revision and extension of Jacobson and Rush (2008). It determines the orbits of the five main
Uranian satellites plus Puck (new; "the largest of the 10 small satellites discovered by Voyager 2", section 1, p2), the orbits of the
Uranian rings, the gravity field of the Uranian system (system GM, satellite GMs, J2, J4, with J6 not estimated) and the orientation of
the Uranus pole, in a combined weighted least-squares fit.

Data (section 3, p2-3):
- Earth-based astrometry of the satellites 1911 to 2013 (filar micrometer, photographic, CCD, transit, mutual events, Hubble Space Telescope
  relative positions of Puck and Miranda or Ariel, Titania stellar occultation timing). Observation statistics are Tables 3 and 5.
- Ring stellar occultations 1977 to 1992 (Earth-based, Voyager 2 photopolarimeter, Voyager 2 radio occultation), from French et al. (1986,
  1987, 1988, 1996), Elliot et al. (1987), Millis et al. (1987) and a 1982 Palomar set.
- Voyager 2 imaging (1985 August 22 to 1986 January 31) in original pixel-line form and Voyager 2 Doppler, range and delta-DOR, 1985 August 22 to 1986 February 13.

Dynamical model (section 2, p2): gravitational forces from the Sun, the solar-system planets and the Uranian satellites (Puck massless);
planet gravity as a spherical harmonic expansion retaining J2, J4, J6; ring mass ignored. READ quote (p2): "JPL planetary ephemeris DE430
(Folkner et al. 2014) provides the positions and GMs of the Sun and planets." Equations of motion of the satellites from Peters (1981), in
Cartesian coordinates centred on the Uranian system barycentre and referenced to the ICRF; satellite orbits integrated with a variable-order,
variable-step-size Gauss-Jackson method (Jackson 1924). The spacecraft is integrated with JPL's Orbit Determination Program.

The ephemeris the paper produces: READ (p10, concluding remarks): "Ephemerides for all of the satellites are available electronically from the
JPL Horizons online solar system data and ephemeris computation service and from NASA's Navigation and Ancillary Information Facility. The
ephemeris designation is URA111." Relation to the kernel: INFERRED (the paper says only that URA111 is this ephemeris; the project's
`v4_uranus_strict.py` uses ura111.bsp, so the kernel is this paper's fit). The time span of the fitted orbits and of the mean-element fit is
1900 to 2100 (READ, section 4.2, p4: "obtained by fitting precessing ellipses to the integrated orbits over 1900-2100"); the orbit uncertainty
statement covers 2000 to 2050 (Table 7).

Accuracy claims: READ (section 4.2, p6, quoting): "we have found that actual orbit errors are rarely twice our formal errors; we gauge actual
errors by the size of the orbit change required to fit newly acquired data." Table 7 gives the maximum formal uncertainties in 2000-2050
(radial, tangential, normal, and tangential error growth rate); for example Titania 35, 160, 50 km and 3 km/yr (section 5 below).

Other statements of note:
- READ (p4): "Our semimajor axes, eccentricities, and inclinations are comparable those in the GUST86 theory (Laskar & Jacobson 1987)."
- READ (p4): the Miranda (V), Ariel (I), Umbriel (II) near-commensurability: `lambda_V' - 3 lambda_I' + 2 lambda_II' = -0.0785 deg/day`; the
  Greenberg (1976) perturbation on Miranda's longitude has amplitude about 3500 km and period about 12.55 yr.
- READ (section 4.5, p8): "There are no resonances in the Uranian satellite system to enhance the determination of the satellite GMs."
- READ (p8): J6 could not be estimated: "we set its value to 0 in our dynamical model but took into account its a priori uncertainty when
  computing the statistics for the other parameters (i.e., it was treated as a consider parameter...)". "We use the system GM rather than the
  Uranus GM as our gravity parameter; given the magnitude of the system GM the system is predominately the planet."
- READ (p8): the Jacobson et al. (1992) analysis "erroneously included a considerable amount of simultaneous two-way and three-way Voyager 2
  Doppler"; the extraneous three-way data was removed in this fit.

## 2. Parameters estimated (section 4.1, p3)

READ: satellites: epoch position and velocity of each satellite; rings: elements of each ring. Common gravitational parameters: the GMs of the
Uranian system and the satellites; J2, J4, J6 of Uranus; the Right Ascension and Declination of Uranus's pole at epoch J2000. Plus observation-model
biases (twelve listed items), spacecraft parameters. Table 13 gives parameter counts per data subset.

## 3. Reference plane and frame

READ (Table 2 note, p5): "The longitudes are measured from the intersection of the Uranus equator and the ICRF reference plane." The Table 2
elements are therefore Uranus-equatorial elements with the node reference on the ICRF plane. READ (section 2, p2): the satellite orbits are integrated
in Cartesian coordinates "centered at the Uranian system barycenter and referenced to the ICRF" and Table 1 state vectors are in that frame
(INFERRED: Table 1 does not name a frame in its caption; the Ariel z of -2109 km out of about 190,000 km shows it is not Uranus-equatorial, consistent with ICRF). Ring elements (Table 8) are referred to the Uranus equator; ring longitudes' epoch is 1977 March 10 20:00:00 UTC (Table 8 note).

## 4. Pole model (section 2, p2; Appendix B, p11-12; Figures 2 and 3, p8, p12)

READ: the pole is defined by ICRF right ascension and declination "specified in terms of trigonometric series. The series were developed by
fitting the angles to the numerically integrated polar motion of a rigid body Uranus with torques from the Sun and the satellites." So the pole is
NOT constant: it is a time series, but the variation is small. Appendix B integrates equations B6 and B7 (Uranus spin angular momentum torqued by the Sun and five satellites) over 1600 to 2600. READ (p12): "We set J2 = 3510.7 x 10^-6, gamma = 0.2269 (Hubbard & Marley 1989), and s = 501.1600928 deg day^-1 (Desch et al. 1986). The initial pole orientation angles (at J2000) were alpha = 77.309980 and delta = 15.172395."
(gamma is the axial moment of inertia factor, s the spin rate.) The integrated change is Figures 2 (right ascension) and 3 (declination): amplitude a few 10^-3 deg over 1600 to 2600.

The fitted representation, READ p12, T in Julian centuries from J2000 (degrees):

    alpha = 77.309980 + 0.000173 T + 0.000885 sin S1 + 0.000180 sin S2 + 0.000098 sin S3 + 0.000075 sin S4 - 0.000818 sin S5
    delta = 15.172395 + 0.000019 T + 0.000851 cos S1 + 0.000173 cos S2 + 0.000094 cos S3 + 0.000072 cos S4 + 0.000818 cos S5

    S1 = 328.616724 + 26.9601 T       S2 = 259.275089 + 2024.7285 T      S3 = 102.827444 + 182.8030 T
    S4 = 185.361668 + 276.4108 T      S5 = 137.359959

(digits of the coefficient lists are from the page image; the S5 coefficient 0.000818 appears with opposite sign in alpha and delta as read; unclear whether the printed alpha and delta coefficients of S2 are 0.000180 and 0.000173 as transcribed.) READ quote: "the added S5 term forces the series to sum to 0 at epoch. The series can be simplified by adding the S5 term to the constant term." The paper also states "d is days from J2000" but no equation uses d (READ). Estimated pole at J2000 (Table 12 Current Results): alpha 77.310 +/- 0.002, delta 15.172 +/- 0.002 deg.
READ (p10): "The Uranus pole precession rate is about 1.3 mas yr^-1 which is too small to be detected with the current data set".

## 5. Tables, transcribed

### Table 1 (p4): State Vectors at 1985 August 1 (TDT). Position (km), velocity (km s^-1). "16 digits are quoted to facilitate future orbit integrations" (p4).

| Satellite | x (km) | y (km) | z (km) | vx (km/s) | vy (km/s) | vz (km/s) |
| --- | --- | --- | --- | --- | --- | --- |
| Ariel | -185785.2177189803 | 42477.81018746200 | -2109.273462727150 | -0.384730129923274 | -1.393752472818678 | 5.325004225204424 |
| Umbriel | -176566.9475784755 | 89016.12833946147 | -176154.9418970623 | -3.350588413391897 | -0.153568184837806 | 3.273855499527411 |
| Titania | -221240.1941919138 | 145452.9878060127 | -346697.1461249496 | -3.049048602775958 | 0.138409610142017 | 1.991437563896877 |
| Oberon | -155108.4287158760 | 181606.6634411168 | -532879.3651011410 | -2.962407899779837 | 0.385864135361727 | 0.993694238058708 |
| Miranda | -127430.9607930668 | 23792.64617013941 | -3464.554580724168 | -0.422514450329333 | -1.271890082631948 | 6.552338419694388 |
| Puck | -24369.49145882789 | 27011.79870872380 | -77882.20359704649 | -7.667367460395494 | 1.014093590954378 | 2.753303665833902 |

(Arithmetic caution: the Umbriel z in the image reads -176154.94, consistent in magnitude with its radius about 266,000 km given x, y; unclear digits in none of the other rows, but all 16-digit strings were read from a page image and should be checked against the kernel before use as initial conditions.)

### Table 2 (p5): Mean Equatorial Orbital Elements at 2000 January 1.5 (TDT)

| Element | Ariel | Umbriel | Titania | Oberon | Miranda | Puck |
| --- | --- | --- | --- | --- | --- | --- |
| a (km) | 190930. | 265982. | 436282. | 583449. | 129858. | 86005. |
| e | 0.00122 | 0.00394 | 0.00123 | 0.00140 | 0.00135 | 0.00019 |
| i (deg) | 0.0167 | 0.0796 | 0.1129 | 0.1478 | 4.4072 | 0.3562 |
| lambda (deg) | 203.0922 | 251.2562 | 281.5845 | 352.5701 | 328.6961 | 266.0933 |
| varpi (deg) | 83.3294 | 352.9610 | 228.3816 | 212.8435 | 256.3256 | 98.9184 |
| Omega (deg) | 289.7415 | 195.4845 | 26.4013 | 30.4840 | 100.7027 | 216.0360 |
| lambda-dot (deg/day) | 142.8356506 | 86.8688753 | 41.3514187 | 26.7394835 | 254.6906573 | 472.5445452 |
| varpi-dot (deg/yr) | 6.2308 | 2.8428 | 0.9993 | 0.2680 | 20.0409 | 80.8938 |
| Omega-dot (deg/yr) | -6.0839 | -2.6316 | 0.0208 | -0.4210 | -20.2470 | -80.8624 |

Notes (READ): "The longitudes are measured from the intersection of the Uranus equator and the ICRF reference plane." Text (p4): the elements are "obtained by fitting precessing ellipses to the integrated orbits over 1900-2100"; a is the semi-major axis, e eccentricity, i inclination, lambda mean longitude, varpi longitude of periapsis, Omega longitude of ascending node.

### Table 3 (p5-6): Statistics of the Residuals of the Astrometric Observations
Columns: source, dates, type, number, rms (arcsec), second type, number, rms (arcsec). READ from the page image; this table is observational bookkeeping and is transcribed at lower confidence (unclear digits possible in the rms columns).

Filar micrometer

| Source | Dates | Type | No. | rms | Type | No. | rms |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Aitken (1912) | 1911 | rho dtheta | 23 | 0.279 | d rho | 23 | 0.314 |
| Barnard (1912) | 1911 | rho dtheta | 23 | 0.440 | d rho | 46 | 0.328 |
| Barnard (1915) | 1913 | rho dtheta | 12 | 0.309 | d rho | 12 | 0.239 |
| Harris (1949) | 1913-1948 | da cos d | 62 | 0.173 | d delta | 61 | 0.182 |
| Aitken (1914) | 1914 | rho dtheta | 24 | 0.252 | d rho | 24 | 0.327 |
| Barnard (1916) | 1915 | rho dtheta | 18 | 0.298 | d rho | 18 | 0.292 |
| Barnard (1919) | 1916-1918 | rho dtheta | 45 | 0.342 | d rho | 42 | 0.324 |
| Barnard (1927) | 1919-1922 | rho dtheta | 53 | 0.353 | d rho | 54 | 0.287 |
| Hall (1921) | 1920 | rho dtheta | 21 | 0.321 | d rho | 21 | 0.347 |
| Hall & Bower (1923) | 1922 | rho dtheta | 11 | 0.126 | d rho | 11 | 0.231 |
| Struve (1928) | 1927-1928 | rho dtheta | 41 | 0.333 | d rho | 41 | 0.237 |
| Steavenson (1948) | 1947-1948 | rho dtheta | 23 | 0.441 | | | |
| Steavenson (1964) | 1949 | rho dtheta | 7 | 0.202 | | | |

Photographic

| Source | Dates | Type | No. | rms | Type | No. | rms |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Nicholson (1915) | 1914 | da cos d | 19 | 0.301 | d delta | 19 | 0.433 |
| Sytinskaja (1930) | 1926 | rho dtheta | 62 | 0.530 | d rho | 62 | 0.438 |
| van Biesbroeck (1970) | 1948-1964 | da cos d | 757 | 0.158 | d delta | 757 | 0.149 |
| Whitaker & Greenberg (1973) | 1960-1973 | rho dtheta | 16 | 0.179 | | | |
| Tomita & Soma (1979) | 1964-1977 | rho dtheta | 95 | 0.363 | d rho | 123 | 0.443 |
| van Biesbroeck (1970) | 1966 | da cos d | 43 | 0.202 | d delta | 43 | 0.172 |
| van Biesbroeck et al. (1976) | 1966 | da cos d | 55 | 0.148 | d delta | 55 | 0.123 |
| Soulie (1968) | 1966-1967 | da cos d | 4 | 0.398 | d delta | 4 | 0.377 |
| Soulie (1972) | 1968-1969 | da cos d | 57 | 0.412 | d delta | 57 | 0.388 |
| Soulie (1975) | 1970-1971 | da cos d | 30 | 0.497 | d delta | 30 | 0.334 |
| Soulie (1978) | 1973-1974 | da cos d | 67 | 0.423 | d delta | 67 | 0.444 |
| Mulholland (1985) | 1974-1982 | da cos d | 215 | 0.085 | d delta | 215 | 0.083 |
| Walker et al. (1978) | 1975-1977 | da cos d | 51 | 0.050 | d delta | 51 | 0.042 |
| Veillet (1983b) (Pic du Midi) | 1977-1982 | da cos d | 427 | 0.090 | d delta | 427 | 0.103 |
| Veillet (1983b) (OHP) | 1977-1982 | da cos d | 58 | 0.157 | d delta | 58 | 0.171 |
| Veillet (1983b) (ESO) | 1977-1982 | da cos d | 448 | 0.057 | d delta | 448 | 0.049 |
| Veillet (1983b) (CFH) | 1977-1983 | da cos d | 146 | 0.053 | d delta | 146 | 0.097 |
| Harrington & Walker (1984) | 1979-1983 | da cos d | 283 | 0.056 | d delta | 283 | 0.055 |
| C. Veillet (1985, private comm.) | 1982-1984 | da cos d | 568 | 0.068 | d delta | 568 | 0.063 |
| Debehogne et al. (1981) | 1980 | da cos d | 22 | 0.583 | d delta | 23 | 0.452 |
| Veiga et al. (2003) | 1982-1988 | da cos d | 1395 | 0.095 | d delta | 1397 | 0.077 |
| Standish (1996) | 1983-1986 | da cos d | 481 | 0.210 | d delta | 481 | 0.190 |
| R. S. Harrington (1985, private comm.) | 1985 | da cos d | 41 | 0.025 | d delta | 41 | 0.036 |
| Walker & Harrington (1988) | 1985-1986 | da cos d | 104 | 0.038 | d delta | 104 | 0.041 |
| Chanturiya et al. (2002) | 1987-1994 | da cos d | 91 | 0.557 | d delta | 91 | 0.576 |
| Kulyk et al. (1990) | 1990 | da cos d | 40 | 0.291 | d delta | 40 | 0.239 |

CCD

| Source | Dates | Type | No. | rms | Type | No. | rms |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Pascu et al. (1987) | 1981-1985 | da cos d | 76 | 0.137 | d delta | 76 | 0.108 |
| Veiga et al. (2003) | 1989-1998 | da cos d | 5221 | 0.080 | d delta | 5221 | 0.089 |
| Jones et al. (1998) | 1990-1991 | da cos d | 494 | 0.038 | d delta | 489 | 0.031 |
| Shen et al. (2002) | 1995-1997 | da cos d | 442 | 0.054 | d delta | 442 | 0.056 |
| Stone & Harris (2000) | 1998 | alpha | 100 | 0.155 | delta | 102 | 0.137 |
| Qiao et al. (2013) | 1998-2007 | alpha | 2358 | 0.081 | delta | 2358 | 0.097 |
| Stone (2000) | 1999 | alpha | 117 | 0.134 | delta | 118 | 0.131 |
| Descamps et al. (2002) | 1999-2000 | alpha | 60 | 0.003 | delta | 60 | 0.013 |
| W. M. Owen (2009, private comm.) | 1999-2009 | alpha | 186 | 0.062 | delta | 184 | 0.110 |
| Stone (2001) | 2000-2001 | alpha | 245 | 0.143 | delta | 248 | 0.133 |
| B. Gladman (2001, private comm.) | 2001 | da cos d | 6 | 0.341 | d delta | 6 | 0.271 |
| Stone (2005) | 2002-2005 | alpha | 230 | 0.118 | delta | 232 | 0.132 |
| Garradd & McNaught (2003) | 2003 | alpha | 6 | 0.079 | delta | 6 | 0.116 |
| M. R. Showalter (2007, private comm.) | 2003-2006 | da cos d | 180 | 0.006 | d delta | 179 | 0.017 |
| Veiga & Bourget (2006) | 2004 | alpha | 287 | 0.072 | delta | 287 | 0.078 |
| Izmailov et al. (2007) | 2005 | alpha | 20 | 0.068 | delta | 20 | 0.145 |
| Monet (2007) | 2005-2006 | alpha | 71 | 0.146 | delta | 71 | 0.152 |
| Khovritchev (2009) | 2007 | da cos d | 22 | 0.073 | d delta | 22 | 0.038 |
| Harris (2013) | 2007-2013 | alpha | 337 | 0.139 | delta | 337 | 0.167 |

Transit: Arlot et al. (2008), 1997-2005, alpha, 227, 0.153; delta, 227, 0.185.
Mutual events: Christou et al. (2009), 2007, da cos d, 1, 0.001; d delta, 2, 0.008. Mallama et al. (2009), 2007, da cos d, 2, 0.038; d delta, 2, 0.009. Arlot et al. (2013), 2007-2008, da cos d, 34, 0.010; d delta, 34, 0.017.
Titania stellar occultation timing: Widemann et al. (2009), 2001, 2003, number 58, rms 0.150 (unit not legible on the image; unclear).
Type indicators (READ, p4): rho dtheta is planet-relative separation and position angle (position-angle residuals scaled by the separation); da cos d and d delta are differential RA (scaled by cos declination) and declination; alpha, delta are RA and declination of the satellite.

### Table 4 (p6): Titania Occultation Star Catalog Positions and Corrections

| No. | Cat. R.A. (deg) | Cat. Decl. (deg) | d alpha cos delta (mas) | d delta (mas) |
| --- | --- | --- | --- | --- |
| Hipp. 106829 | 324.55809702 | -14.91006013 | 14.3 +/- 3.0 | 8.0 +/- 3.7 |
| TYC 5806-696-1 | 333.97730960 | -11.61565510 | 38.5 +/- 3.2 | 18.4 +/- 4.1 |

### Table 5 (p6): Statistics of the Residuals for the Puck Astrometric Observations

| Source | Dates | Type | No. | rms | Type | No. | rms |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Pascu et al. (1998) | 1994 | rho dtheta | 31 | 0.012 | d rho | 32 | 0.006 |
| Descamps et al. (2002) | 1999-2000 | alpha | 30 | 0.032 | delta | 30 | 0.143 |
| M. R. Showalter (2007, private comm.) | 2003-2006 | da cos d | 152 | 0.004 | d delta | 152 | 0.013 |
| Veiga & Bourget (2006) | 2004 | alpha | 135 | 0.102 | delta | 135 | 0.091 |

### Table 6 (p6): Voyager Imaging Residuals rms

| Object | No. | Sample | Line | Object | No. | Sample | Line |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ariel | 109 | 0.263 | 0.178 | Oberon | 64 | 0.167 | 0.336 |
| Umbriel | 103 | 0.215 | 0.278 | Miranda | 116 | 0.241 | 0.221 |
| Titania | 65 | 0.192 | 0.267 | Puck | 49 | 0.177 | 0.145 |
| Stars | 879 | 0.253 | 0.250 | | | | |

### Table 7 (p6): Satellite Orbit Uncertainties for the Period 2000-2050 (maximum formal)

| Satellite | R (km) | T (km) | N (km) | T-dot (km/yr) |
| --- | --- | --- | --- | --- |
| Ariel | 15 | 150 | 30 | 3 |
| Umbriel | 20 | 150 | 50 | 3 |
| Titania | 35 | 160 | 50 | 3 |
| Oberon | 30 | 160 | 60 | 3 |
| Miranda | 6 | 250 | 25-125 | 5 |
| Puck | 6 | 150 | 30 | 3 |

(R radial, T tangential, N normal; READ p6-7: orbital period errors make the tangential error grow at the rate T-dot; Miranda's inclination makes the precession-rate error add about 2 km/yr to N.)

### Table 8 (p7): Ring Elements Referred to the Uranus Equator, Number of Observations (N), rms of the Residuals
Epoch for the longitudes: 1977 March 10 20:00:00 (UTC). Rates in deg/yr. e is multiplied by 10^3.

| Element | 6 | 5 | 4 | alpha |
| --- | --- | --- | --- | --- |
| a (km) | 41837.27 +/- 0.30 | 42234.85 +/- 0.20 | 42570.99 +/- 0.21 | 44718.43 +/- 0.22 |
| e (x10^3) | 1.012 +/- 0.008 | 1.898 +/- 0.005 | 1.060 +/- 0.007 | 0.761 +/- 0.006 |
| i (deg) | 0.062 +/- 0.003 | 0.055 +/- 0.002 | 0.032 +/- 0.001 | 0.015 +/- 0.001 |
| varpi (deg) | 242.50 +/- 0.72 | 170.04 +/- 0.38 | 126.98 +/- 0.39 | 333.33 +/- 0.44 |
| Omega (deg) | 11.49 +/- 1.41 | 286.40 +/- 0.87 | 89.81 +/- 2.45 | 61.06 +/- 5.14 |
| varpi-dot (deg/yr) | 1008.767 +/- 0.052 | 975.767 +/- 0.046 | 948.935 +/- 0.046 | 798.166 +/- 0.029 |
| Omega-dot (deg/yr) | -1006.775 +/- 0.052 | -973.875 +/- 0.046 | -947.124 +/- 0.044 | -796.786 +/- 0.029 |
| N | 30 | 39 | 39 | 44 |
| rms (km) | 0.33 | 0.24 | 0.33 | 0.36 |

| Element | beta | eta | epsilon | lambda |
| --- | --- | --- | --- | --- |
| a (km) | 45661.05 +/- 0.18 | 47176.02 +/- 0.27 | 51149.21 +/- 0.30 | 50024.16 +/- 0.96 |
| e (x10^3) | 0.441 +/- 0.004 | 0.003 +/- 0.009 | 7.930 +/- 0.007 | 0.0 |
| i (deg) | 0.005 +/- 0.001 | 0.001 +/- 0.001 | 0.001 +/- 0.001 | 0.0 |
| varpi (deg) | 224.83 +/- 0.65 | 321.24 +/- 103.53 | 214.591 +/- 0.154 | |
| Omega (deg) | 311.60 +/- 12.54 | 186.11 +/- 54.12 | 171.69 +/- 104.374 | |
| varpi-dot (deg/yr) | 741.742 +/- 0.024 | 661.372 +/- 0.023 | 497.941 +/- 0.018 | |
| Omega-dot (deg/yr) | -740.513 +/- 0.024 | -660.345 +/- 0.023 | -497.284 +/- 0.018 | |
| N | 41 | 35 | 44 | 4 |
| rms (km) | 0.21 | 0.50 | 0.60 | 0.80 |

(The 4-ring column's varpi-dot uncertainty and the epsilon Omega uncertainty 104.374 are unclear on the image.) Ring plane model: Appendix A.

### Table 9 (p7): Observatory Time Offsets (s)

| Station | Offset (s) | Station | Offset (s) |
| --- | --- | --- | --- |
| European Southern Obs. (a) | -0.072 +/- 0.063 | European Southern Obs. (b) | -0.127 +/- 0.024 |
| European Southern Obs. (c) | 0.062 +/- 0.029 | Las Campanas Obs. | -0.030 +/- 0.022 |
| Pic du Midi | 3.685 +/- 0.033 | Tenerife Ingress | -0.070 +/- 0.035 |
| Tenerife Egress | 0.368 +/- 0.035 | DSS43 | -0.011 +/- 0.023 |
| PPS (beta Per) | -0.054 +/- 0.053 | PPS (sigma Sgr) | 0.450 +/- 0.526 |

Notes: (a) 1980 August 15, (b) 1982 April 22, (c) 1992 July 11.

### Table 10 (p8): Star Catalog Positions and Corrections

| Star | Cat. | Cat. R.A. (deg) | Cat. Decl. (deg) | d alpha cos d (mas) | d delta (mas) |
| --- | --- | --- | --- | --- | --- |
| S77 | Hipp. | 219.54921290 | -14.95473933 | 7.308 +/- 0.054 | -5.913 +/- 0.131 |
| KM11 | UCAC2 | 233.40997800 | -18.90128950 | -294.172 +/- 0.151 | 27.147 +/- 0.091 |
| KM12 | UCAC2 | 229.54172530 | -17.99479950 | -196.954 +/- 0.065 | 125.594 +/- 0.211 |
| KME13 | Hipp. | 237.10643560 | -19.77402446 | -3.149 +/- 0.042 | 22.473 +/- 0.119 |
| KME14 | Hipp. | 242.14934740 | -20.80743248 | 8.837 +/- 0.053 | -21.867 +/- 0.086 |
| KME15 | UCAC2 | 241.79331860 | -20.74504340 | 120.387 +/- 0.039 | 23.800 +/- 0.095 |
| KME16 | UCAC2 | 240.36678420 | -20.48854530 | 1.543 +/- 0.052 | 97.397 +/- 0.078 |
| KME17B | Hipp. | 247.63035990 | -21.74199001 | -35.879 +/- 0.047 | 4.279 +/- 0.054 |
| U23 | UCAC2 | 256.37848620 | -22.87389030 | 97.926 +/- 0.038 | 33.601 +/- 0.058 |
| U25 | UCAC2 | 255.59000530 | -22.80714590 | 64.089 +/- 0.041 | -291.815 +/- 0.047 |
| U28 | UCAC2 | 261.49121500 | -23.29306250 | 105.295 +/- 0.031 | 3.573 +/- 0.057 |
| U103 | UCAC2 | 287.39833590 | -22.91141370 | -6.625 +/- 0.075 | 13.173 +/- 0.044 |
| sigma Sgr | Hipp. | 283.81633188 | -26.29651960 | 0.0 | 0.0 |
| beta Per | Hipp. | 47.04231934 | +40.95561187 | 0.0 | 0.0 |

### Table 11 (p8): Voyager 2 Navigation Performance, Uranus B-plane

| Source | B.R (km) | B.T (km) | SMAA (km) | SMIA (km) | theta (deg) | TCA (H:M:S) | sigma TCA (s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Taylor et al. (1986) | 25384 | 128680 | 1 | 1 | 70 | 17:59:46.5 | 0.07 |
| Jacobson & Rush (2008) | 25384.9 | 128678.8 | 0.9 | 0.1 | 101 | 17:59:46.544 | 0.004 |
| This work | 25385.5 | 128678.5 | 0.6 | 0.05 | 101 | 17:59:46.531 | 0.002 |

### Table 12 (p9): Uranus Gravity Parameters and Pole. READ note: "The reference radius for the Uranus zonal harmonics is 25559 km. (a) From French et al. (1988). (b) Not estimated (see text)."
GM in km^3 s^-2; J2, J4, J6 multiplied by 10^6; angles in degrees.

First block:

| Parameter | Jacobson et al. (1992) | Current Results | Astrometry Only | Voyager Only |
| --- | --- | --- | --- | --- |
| GM System | 5794548.6 +/- 1.5 | 5794556.4 +/- 4.3 | 5795232.8 +/- 625.1 | 5794566.1 +/- 9.3 |
| GM Ariel | 90.3 +/- 2.3 | 83.5 +/- 1.4 | 68.2 +/- 10.9 | 132.8 +/- 10.4 |
| GM Umbriel | 78.2 +/- 2.3 | 85.1 +/- 1.9 | 99.4 +/- 11.2 | 49.4 +/- 17.8 |
| GM Titania | 235.3 +/- 1.9 | 226.9 +/- 4.1 | 224.4 +/- 7.0 | 214.2 +/- 9.2 |
| GM Oberon | 201.1 +/- 1.8 | 205.3 +/- 5.8 | 193.1 +/- 10.1 | 201.1 +/- 9.6 |
| GM Miranda | 4.4 +/- 0.2 | 4.3 +/- 0.2 | 4.0 +/- 0.7 | 5.0 +/- 0.3 |
| J2 (x10^6) | 3513.2 +/- 0.3 (a) | 3510.7 +/- 0.7 | 3606.2 +/- 83.1 | 3595.0 +/- 41.0 |
| J4 (x10^6) | -31.8 +/- 0.4 (a) | -34.2 +/- 1.3 | 716.0 +/- 797.3 | -239.1 +/- 177.1 |
| J6 (x10^6) | | 0.0 +/- 1.0 (b) | 0.0 +/- 1.0 (b) | 0.0 +/- 1.0 (b) |
| alpha (deg) | 77.311 +/- 0.010 (a) | 77.310 +/- 0.002 | 77.326 +/- 0.018 | 77.691 +/- 0.246 |
| delta (deg) | 15.174 +/- 0.010 (a) | 15.172 +/- 0.002 | 15.157 +/- 0.017 | 15.002 +/- 0.718 |

Second block:

| Parameter | Astrometry & Voyager | Rings Only | Astrometry & Rings | Voyager & Rings |
| --- | --- | --- | --- | --- |
| GM System | 5794564.0 +/- 5.6 | 5794519.0 +/- 3460.7 | 5795208.1 +/- 623.8 | 5794574.6 +/- 8.7 |
| GM Ariel | 78.8 +/- 3.3 | | 77.6 +/- 4.3 | 93.5 +/- 6.7 |
| GM Umbriel | 89.6 +/- 3.4 | | 89.8 +/- 5.2 | 111.1 +/- 13.2 |
| GM Titania | 225.6 +/- 4.5 | | 224.6 +/- 7.0 | 220.1 +/- 6.5 |
| GM Oberon | 200.0 +/- 6.3 | | 194.0 +/- 10.1 | 209.6 +/- 9.3 |
| GM Miranda | 4.1 +/- 0.2 | | 4.5 +/- 0.5 | 4.4 +/- 0.3 |
| J2 (x10^6) | 3520.3 +/- 11.2 | 3510.5 +/- 1.3 | 3510.3 +/- 0.7 | 3510.5 +/- 0.7 |
| J4 (x10^6) | -47.3 +/- 62.2 | -34.4 +/- 1.3 | -34.4 +/- 1.3 | -34.5 +/- 0.7 |
| J6 (x10^6) | 0.0 +/- 1.0 (b) | 0.0 +/- 1.0 (b) | 0.0 +/- 1.0 (b) | 0.0 +/- 1.0 (b) |
| alpha (deg) | 77.322 +/- 0.011 | 77.310 +/- 0.003 | 77.310 +/- 0.002 | 77.309 +/- 0.002 |
| delta (deg) | 15.174 +/- 0.011 | 15.172 +/- 0.002 | 15.172 +/- 0.002 | 15.172 +/- 0.003 |

The system GM printed in the table is in km^3 s^-2 (the table header reads "GM (km^3 s^-2)"). Unclear digits: the "Voyager Only" and "Voyager & Rings" J4 uncertainties (177.1 and 0.7 as read; the Voyager & Rings J4 uncertainty may be 1.3).

### Table 13 (p10): Parameter Sets and the Number of Parameters Contained in the Set That Were Fit to Each Data Subset

| Parameter set | Astrometry only | Voyager only | Astrometry and Voyager |
| --- | --- | --- | --- |
| Gravity and pole | 10 | 10 | 10 |
| Satellite states | 36 | 36 | 36 |
| Micrometer observation biases | 22 | | 22 |
| Absolute astrometry biases | 72 | | 72 |
| CCD observation orientation and scale | 198 | | 198 |
| Star position corrections | 4 | | 4 |
| Spacecraft state | | 6 | 6 |
| RTG acceleration | | 2 | 2 |
| Spacecraft non-gravitational accel. | | 542 | 542 |
| Trajectory correction maneuver | | 3 | 3 |
| Impulsive maneuvers | | 32 | 32 |
| Tracking data ionosphere delays | | 36 | 36 |
| Range biases | | 3 | 3 |
| Quasar location | | 2 | 2 |
| Spacecraft camera pointing angles | | 936 | 936 |
| Imaging data phase corrections | | 5 | 5 |

| Parameter set | Rings only | Astrometry & Rings | Voyager & Rings |
| --- | --- | --- | --- |
| Gravity and pole | 10 | 10 | 10 |
| Satellite states | | 36 | 36 |
| Micrometer observation biases | | 22 | |
| Absolute astrometry biases | | 72 | |
| CCD observation orientation and scale | | 198 | |
| Star position corrections | 24 | 28 | 24 |
| Ring elements | 36 | 36 | 36 |
| Observatory time offsets | 10 | 10 | 10 |
| Spacecraft state | | | 6 |
| RTG acceleration | | | 2 |
| Spacecraft non-gravitational accel. | | | 542 |
| Trajectory correction maneuver | | | 3 |
| Impulsive maneuvers | | | 32 |
| Tracking data ionosphere delays | | | 36 |
| Range biases | | | 3 |
| Quasar location | | | 2 |
| Spacecraft camera pointing angles | | | 936 |
| Imaging data phase corrections | | | 5 |

(Table 13 has no all-data column: the all-data fit is the "Current Results" column of Table 12, INFERRED from the table-12 text and the absence of that set here.)

## 6. Which solution is final

READ (section 4.5, p8): "Table 12 contains the Uranian system gravity parameters and their formal 1 sigma errors found by Jacobson et al. (1992) and those obtained from the current analysis... The remaining columns in Table 12 show the gravity results obtained from fits to various subsets of our data." So the second column, "Current Results", is the fit to all the data (the paper never labels any column "recommended" or "final"; INFERRED from this sentence, from the correlation statement "When all of the data are used the correlation between the harmonics becomes 0.978 and between the Ariel and Umbriel GMs is 0.810" (p9), and from the J2 = 3510.7 used in the Appendix B pole integration (p12), which equals this column). Concluding remarks (p10), READ: "Our gravity field parameter determination required the combination of the satellite astrometry, ring occultations, and radiometric tracking of the Voyager 2 spacecraft during its passage through the Uranian system."

Recommended (our reading of the paper's preferred column): J2 = (3510.7 +/- 0.7) x 10^-6, J4 = (-34.2 +/- 1.3) x 10^-6, J6 set to 0 (not estimated; a priori sigma 1.0 x 10^-6), reference radius 25,559 km, pole alpha = 77.310 +/- 0.002 deg, delta = 15.172 +/- 0.002 deg (ICRF, epoch J2000), system GM 5,794,556.4 +/- 4.3 km^3 s^-2.

The Jacobson et al. 1992 column's J2 and J4 are marked (a) "From French et al. (1988)" in the note (READ), not estimated by Jacobson et al. 1992.

## 7. Constants audit: project values against this paper

Project sources: `src/cyclerfinder/core/satellites.py` (`PRIMARIES["Uranus"]` line 70; `SATELLITES` lines 251 to 259) and `src/cyclerfinder/data/validation/v4_uranus.py` (`URANUS_J2` line 131, `URANUS_R_EQ_KM` line 148). All "difference" entries are COMPUTED (project minus paper).

### 7.1 Gravity parameters

| Quantity | Project value | Paper (Current Results, Table 12) | Difference | Paper is the plausible origin |
| --- | --- | --- | --- | --- |
| Uranus system GM (km^3/s^2) | 5.7945564e6 = 5794556.4 | 5794556.4 +/- 4.3 | 0 | Yes: identical to 8 digits. Project comment says "JPL DE440 planetary constants"; the DE440 value is the same Jacobson 2014 URA111 number (INFERRED: identical digits) |
| Miranda GM | 4.3 | 4.3 +/- 0.2 | 0 | Yes (project comment: JPL SSD phys_par, ref URA111) |
| Ariel GM | 83.5 | 83.5 +/- 1.4 | 0 | Yes |
| Umbriel GM | 85.1 | 85.1 +/- 1.9 | 0 | Yes |
| Titania GM | 226.9 | 226.9 +/- 4.1 | 0 | Yes |
| Oberon GM | 205.3 | 205.3 +/- 5.8 | 0 | Yes |
| `URANUS_J2` | 3509.291e-6 | 3510.7e-6 +/- 0.7 | -1.409e-6 (-0.040 percent; 2.0 sigma of the paper's 0.7, 1.7 sigma combined with French et al. 2024's 0.412) | No: from French et al. 2024 (arXiv:2401.04634) per the code docstring; not checkable here |
| `URANUS_R_EQ_KM` | 25559.0 | "reference radius for the Uranus zonal harmonics is 25559 km" | 0 | Yes as a J2 reference radius; the paper does not give a separate physical equatorial radius (INFERRED: it is not the physical radius, which is 25,559 km as the IAU 2015 nominal value anyway) |
| J4 | not carried | -34.2e-6 +/- 1.3e-6 | | n/a |
| J6 | not carried | 0 (not estimated) | | n/a |
| Retired J2 3.34343e-3 (paired with 25,559 km until 2026-10-04) | removed | not printed anywhere in the paper | | No. COMPUTED: 3.34343e-3 x R_old^2 = 3510.7e-6 x 25559^2 gives R_old = 25559 x sqrt(3510.7 / 3343.43) = 26,190 km; neither that J2 nor a 26,190 km radius appears in this paper |

### 7.2 Semi-major axes: what the paper prints and the project values

READ: Table 2 gives "Mean Equatorial Orbital Elements at 2000 January 1.5 (TDT)", from "fitting precessing ellipses to the integrated orbits over 1900-2100" (p4). These are mean (not osculating) elements, in the Uranus equator frame, with longitudes from the ICRF-equator node. The paper gives no other semi-major axes. The Table 1 state vectors (1985 August 1) are osculating, not elements; semi-major axis could be computed from them but the paper does not print them.

| Satellite | Project a (km) | Paper Table 2 a (km) | Difference (km) | Relative | Match? |
| --- | --- | --- | --- | --- | --- |
| Miranda | 129846 | 129858 | -12 | -9.2e-5 | no |
| Ariel | 190929 | 190930 | -1 | -5.2e-6 | within 1 km (the paper prints integer km) |
| Umbriel | 265986 | 265982 | +4 | +1.5e-5 | no |
| Titania | 436298 | 436282 | +16 | +3.7e-5 | no |
| Oberon | 583511 | 583449 | +62 | +1.06e-4 | no |

INFERRED: the project values (the JPL Satellite Mean Elements page, per the registry comment) are a different mean-element fit (a different epoch or element definition) of the same URA111 ephemeris, not Table 2. The paper cannot be their source. The #890 review's "fitted a" values of about 436,28x and 583,451 km are consistent with Table 2 (583,449 here; the review's 583,451 is a 2 km difference, from a different averaging).

### 7.3 Radii and flyby floors

Not in this paper (READ: no satellite radius appears in any table or section; the only radius is Uranus's zonal reference radius). The project's mean radii (Miranda 235.8, Ariel 578.9, Umbriel 584.7, Titania 788.9, Oberon 761.4 km) and the 50 and 100 km flyby floors come from other sources (the registry comments cite JPL SSD phys_par and Heaton-Longuski 2003); their source is not determined by this paper and no guess is made.

## 8. Mean motions and periods (the #890 question)

READ: Table 2 mean longitude rates, deg/day. COMPUTED: period = 360 / rate; synodic = 360 / (rate difference).

| Satellite | lambda-dot (deg/day) | Period 360/rate (d) | Project registry Kepler period (see below) |
| --- | --- | --- | --- |
| Miranda | 254.6906573 | 1.413479 | |
| Ariel | 142.8356506 | 2.520379 | |
| Umbriel | 86.8688753 | 4.144177 | |
| Titania | 41.3514187 | 8.705868 | 8.706399 (#890 review) |
| Oberon | 26.7394835 | 13.463237 | 13.465984 (#890 review) |
| Puck | 472.5445452 | 0.761833 | |

COMPUTED: the paper's Titania and Oberon periods are 8.705868 d and 13.463237 d, identical to the six printed digits of the URA111 fits the #890 review gave. Titania-Oberon synodic period from the paper: 360 / (41.3514187 - 26.7394835) = 360 / 14.6119352 = 24.6374 d (the review's "real" value 24.63740 d; the model's was 24.63245 d).

COMPUTED check of Kepler with Table 2 a and the paper's own GMs (n = sqrt((GM_sys + GM_sat) / a^3)): Titania 8.705595 d, Oberon 13.463336 d, so the paper's tabulated rates are not exactly Keplerian from the tabulated a (relative difference -3.1e-5 for Titania, +7.4e-6 for Oberon; Miranda +1.6e-4); the difference is the J2 and mutual-perturbation modification of the mean motion that "mean elements" absorb (INFERRED). With the project's registry a and GM (GM sum, Kepler): Titania 8.706074 d, Oberon 13.465482 d. I could not reproduce the review's 8.706399 d and 13.465984 d from the registry numbers with either the sum or the planet GM alone (planet-only gives 8.706245 d and 13.465721 d); the review's method for those two figures should be checked (COMPUTED here; not resolved).

Recommendation (INFERRED): for any model meant to represent the real Titania-Oberon geometry, use the paper's lambda-dot (periods 8.705868 d and 13.463237 d) to set the moons' angular rates, because the registry Kepler periods are 6e-5 (Titania) and 2e-4 (Oberon) too long relative to them (as the #890 review found); take a and the rate from the same element set, or accept that a Keplerian model with Table 2 a gives 8.705595 d and 13.463336 d, which are within 3.1e-5 and 7.4e-6 of the fitted rates.

## 9. Ring model (Appendix A, p10-11), brief

READ: ring plane coordinate system in equinoctial form with p = tan(I/2) sin(Omega), q = tan(I/2) cos(Omega); ring radius `r = a (1 - h^2 - k^2) / (1 + k cos(lambda) + h sin(lambda))`, with h = e sin(varpi), k = e cos(varpi), node and periapsis precessing linearly; secular precession rates from Borderies-Rappaport and Longaretti (1994) augmented by external-satellite terms from Brouwer and Clemence (1961). The rate expressions contain J2, J4, J6, J2^2, J2 J4, J2^3 terms in (R/a), with R the planet's equatorial radius (the zonal reference radius, 25,559 km per Table 12) and mu_0 the planet GM (equation set of Appendix A, p11; not transcribed term by term).

## 10. Appendix B pole torque equations (p11), brief

READ: Equations B1 to B7 give the torque of the Sun and five satellites on the oblate Uranus; the long-term averaged form is B3; with H0 = m0 R^2 gamma s k (B4) the pole rates are `alpha-dot cos(delta) = -(3/2)(J2/(gamma s)) [ ... ]` (B6) and `delta-dot = (3/2)(J2/(gamma s)) [ ... ]` (B7). Constants: J2 = 3510.7e-6, gamma = 0.2269, s = 501.1600928 deg/day.

## 11. Statements the project should not carry forward without care

- The project docstring for the retired J2 cited "Table 4" of this paper. READ: Table 4 is the Titania occultation star catalogue. The J2 is in Table 12.
- The J2 value 3510.7e-6 is paired with R = 25,559 km, as READ in the Table 12 note. Using it with another radius changes J2 R^2 and is the failure mode #894 describes.
- Table 12's "Current Results" J2 uncertainty (0.7e-6) and French et al. 2024's (0.412e-6) are both small; the 1.409e-6 difference is 2.0 sigma on the paper's own error, so the two J2 values are consistent and there is no reason to prefer the older on accuracy grounds (INFERRED). The project's choice of French 2024 for J2 is not contradicted.
- No J2 or J4 value in this paper includes a statement about tides or about whether satellites are treated as point masses beyond "the GMs of the Uranian system and the satellites" being fit (READ p3); the system GM is the planet plus satellites (READ p8, "given the magnitude of the system GM the system is predominately the planet"). INFERRED: using the system GM at the centre together with explicit moon GMs double counts the moon masses at the level of about 1.0e-4 of the system GM (COMPUTED: sum of five moon GMs = 4.3 + 83.5 + 85.1 + 226.9 + 205.3 = 605.1 km^3/s^2; 605.1 / 5794556.4 = 1.04e-4); this is the project's modelling choice, not a statement in the paper.
