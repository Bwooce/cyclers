# Digest: French, Hedman, Nicholson, Longaretti and McGhee-French (2024), "The Uranus system from occultation observations (1977-2006): Rings, pole direction, gravity field, and masses of Cressida, Cordelia, and Ophelia"

Icarus 411:115957 (2024), DOI 10.1016/j.icarus.2024.115957. The held file is the HAL accepted manuscript hal-04778362
(draft version January 8, 2024; submitted to HAL 12 November 2024; CC BY 4.0), 95 pages, text layer present.
Filed in the private paper corpus as
french-hedman-nicholson-longaretti-mcghee-french-2024-uranus-system-occultation-observations-1977-2006-rings-pole-gravity-field-icarus-411-115957-doi-10.1016-j.icarus.2024.115957-hal-accepted-manuscript.pdf

Digested 2026-10-04 from all 95 pages (five reads of at most 20 pages), with the tables cross-checked digit by digit against the
text layer (`pdftotext -layout`). Each statement is marked READ (seen on the page; page, section, table or equation given),
COMPUTED (our arithmetic, shown) or INFERRED (our reading, or a comparison with project code). Page numbers are the manuscript's
printed page numbers (printed page = PDF page minus 1; the cover page is unnumbered). Why the project holds this paper: its
abstract supplied the Uranus J2 now in `data/validation/v4_uranus.py`, read before the paper was held; tasks #890, #894 and #895 depend on the
Uranian constants.

## 0. Headline findings for the project

1. READ (Table 17, p59): the adopted solution is Fit 15, J2 = (3509.291 +/- 0.412) x 10^-6, J4 = (-35.522 +/- 0.466) x 10^-6, J6 fixed at
   0.5 x 10^-6, rho(J2, J4) = 0.9861, reference radius R = 25559 km. The project's `URANUS_J2 = 3509.291e-6` and `URANUS_R_EQ_KM = 25559.0`
   match this exactly. J4 and J6 are not carried by the project.
2. COMPUTED: the J2 the project carried until 2026-10-04, 3.34343e-3 with 25,559 km, is the French et al. (1988) value expressed at that paper's own
   reference radius. READ (p59): "French et al. (1988) ... Fit 1 in Table 17, converted from their assumed reference radius of R = 26200 km to
   R = 25559 km", and Table 17 Fit 1 is J2 = 3513.23 +/- 0.34 (x 10^-6, at 25559 km). 3513.23 x (25559 / 26200)^2 = 3513.23 x 0.951667 = 3343.43
   (x 10^-6), which is 3.34343e-3. So the retired constant was the 1988 J2 at R = 26200 km, filed against the wrong radius, not a Jacobson (2014) number at all.
   This also settles the "26,190 km" radius the Jacobson 2014 digest backed out (3510.7 / 3343.43 gives 26,190 km; 26,200 is the printed radius
   and 3513.23 the corresponding J2).
3. READ (Table 17 lower block, p59): "GM_U 5793950.300 km^3 s^-2, Jacobson (2023)". This is the PLANET's GM, not the system's. COMPUTED evidence: the
   paper's Table 2 mass ratios M_sat/M_Ur are the moon GMs divided by 5793950.3, to every printed digit (x 10^-6: Ariel 82.30 / 5793950.3 = 14.2045, printed 14.204;
   Umbriel 14.8431, printed 14.843; Titania 39.8001, printed 39.800; Oberon 35.8305, printed 35.830; Miranda 0.72489, printed 0.7249). Dividing instead by the Jacobson (2014)
   system GM 5794556.4 gives 14.2030, 14.8415, 39.7960, 35.8267 and 0.72482, which fail the printed digits for Ariel, Umbriel, Titania and Oberon. So the divisor is a planet-only mass.
   The project's `PRIMARIES["Uranus"] = 5.7945564e6` is the Jacobson (2014) SYSTEM GM (Table 12 of that paper); the difference from GM_U is
   5794556.4 - 5793950.3 = 606.1 km^3/s^2 (COMPUTED), against 605.1 for the sum of the five moon GMs of Jacobson (2014) and 610.7 for the five of
   Jacobson (2023) (COMPUTED).
4. READ (Table 2, p9): the newer Jacobson (2023) moon GMs differ from Jacobson (2014) and from the project by up to 3.7 km^3/s^2 (Titania 230.60 vs
   226.9) and the mean semi-major axes by 1 to 30 km (section 6 audit). The project's five moon GMs match Jacobson (2014) only.
5. READ (abstract; section 7.3.2; Table 17): the adopted J2 is 1.174e-6 below Jacobson (2023) (3510.465) and 1.409e-6 below Jacobson (2014)
   (3510.7). The authors attribute the whole difference to systematic effects that earlier fits neglected (section 3 below); the J2 depends on the
   assumed satellite masses (it moves by 1.432e-6 on adding Cordelia and Ophelia alone) and on J6 (dJ2/dJ6 = +0.39909).
6. READ (abstract; Table 13): the adopted pole is alpha = 77.311327 +/- 0.000141 deg, delta = 15.172795 +/- 0.000618 deg at epoch TDB 1986 Jan 19 12:00,
   fixed (precession not modelled in the adopted fit). Jacobson (2014) pole (77.310, 15.172) is 5.4 arcsec away (COMPUTED, section 5).
7. READ (p7, section 3.1.1; p38, section 4.8): this paper's ephemerides are the ura178 series, a JPL update of the ura111 series (the project's
   kernel) that used the ring occultation data. The authors state ura111 and de440 carry "small but measurable systematic errors in the Uranus
   ephemeris" and that an ura111 drift "amounted to several hundred km in the sky plane by the time of the final ring occultation in 2006".
8. COMPUTED (section 7): the J2 difference between any two of the candidate values moves the acceleration at Titania by at most 2.2e-13 km/s^2
   (2.2e-10 m/s^2), a naive constant-acceleration displacement of 12.5 km in 123 days but about 0.14 km as a mean-motion (phase) effect. The
   system-GM double count in the V4-strict lane (INFERRED) is four orders of magnitude larger.
9. FOUND WHILE CHECKING THE CORPUS: a file for Jacobson and Park (2025), "Orbits of Uranus satellites, rings, gravity field, poles", URA182, AJ 169:65,
   DOI 10.3847/1538-3881-ad99d1 is present in the private paper corpus (file timestamp 2026-10-04 21:18) and has no row in `docs/notes/CORPUS_INDEX.md`
   (grep for "jacobson-park" and "URA182": zero hits). It is presumably the published successor to the Jacobson (2023) solution this paper quotes,
   and would supersede the second-hand 2023 numbers below. I did not open it.

## 1. What the paper does

READ (abstract, p1-2; sections 1 and 2, p2-3): "From an analysis of 31 Earth-based stellar occultations and three Voyager 2 occultations spanning
1977-2006 (French et al. 2023a), we determine the keplerian orbital elements of the centerlines (COR) of the nine main Uranian rings to high
accuracy, with typical RMS residuals of 0.2 - 0.4 km and 1-sigma formal errors in a, ae, and a sin i of order 0.1 km, registered on an absolute
radius scale accurate to 0.2 km at the 2-sigma level." This is "Paper 2"; the data are the Paper 1 set (French et al. 2023a, Icarus 395:115474).

Data (section 2, p3-6):
- 31 Earth-based stellar occultations by the narrow rings, 1977 to 2006, mostly infrared K band (2.2 micron), most observed 1980 to 1990 when Uranus crossed
  the Milky Way and the view was nearly pole-on. Ring midtimes and edge times come from square-well model fits to each ring profile (Figures 1 to 3).
  Table 8 (p31-34) lists 30 distinct events from U0 (1977-03-10, KAO) to U0602 (2006-09-20).
- Three Voyager 2 occultations: the radio science occultation (RSS) and the PPS and UVS stellar occultations of sigma Sgr and beta Per.
- Accuracy of the square-well model: READ (p5-6 and Appendix A, Table A1, p78) the mean width difference to the diffraction-corrected Voyager profiles is
  dW = 0.002 km with sigma 0.27 km; the midline accuracy "is reduced by a factor of 1/sqrt(2) to about 0.2 km".

Ring orbit model (section 3.2, p8-10): the RINGFIT code. READ eq. (1): `r = a (1 - e^2) / (1 + e cos f)` with `f = lambda - varpi = lambda - varpi_0 - varpi_dot (t - t_0)`, plus
inclination i, node Omega_0 and regression Omega_dot; the zero-point of inertial longitude is the ascending node of Uranus's equator on Earth's equator of J2000, with the
pole taken in the direction of positive angular momentum, "180 deg from the IAU definition of the Uranus north pole" (p8 and Table 5 note (b)). Optional
secular rates are given by eqs. (2) and (3) (J2, J4, J6 and the satellite terms with Laplace coefficients; the planet term in eq. (2) is
`sqrt(GM/a^3) {(3/2) J2 (R/a)^2 (1 + e^2 - 2 sin^2 i) - (15/4) J4 (R/a)^4 + [(27/64) J2^3 - (45/32) J2 J4 + (105/16) J6] (R/a)^6 + ...}`). Normal modes
are added by eqs. (5) to (8): `Delta r = -A_m cos(m theta)`, `theta = lambda - Omega_P (t - t_0) - delta_m`, pattern speed near `Omega_P = [(m - 1) n + varpi_dot_sec] / m` (eq. 7).

Fit and what is determined (sections 4 to 7):
- Normal modes of the ring centerlines and edges (section 4.1; Tables 5 and 6; Figures 5 to 19): the known gamma m = 0 and delta m = 2 modes, two new gamma outer-Lindblad
  modes (m = -1, -2), a possible gamma m = 3 inner-Lindblad mode, and five satellite-forced modes (Cressida-eta 3:2, Ophelia-gamma 6:5, Cordelia-delta 23:22, Ophelia-epsilon
  14:13 outer edge, Cordelia-epsilon 25:24 inner edge).
- Orbit-fit statistics (section 4.4, p31): 651 data points in the COR fit (RMS per degree of freedom 0.348 km); 1174 in the IER/OER fit (0.655 km).
- The pole direction and radius scale (section 5; Tables 12 and 13), the widths, shapes and masses of the rings (section 6; Tables 14 and 15), the gravity field J2 and J4
  (section 7; Tables 16 and 17), anomalous precession of the alpha, beta and gamma rings (section 8; Table 18), and the masses and densities of Cressida, Cordelia and Ophelia
  (section 9; Tables 3, 19 and 20).
- No fit of satellite orbits is made. READ (section 3.3, p11): the satellite and planetary ephemerides are the JPL ura178 series; this paper does not fit the Uranian
  system's GM, the major-moon GMs or the moon orbits. Those come from Jacobson (2023) (Table 2; GM_U in Table 17).

## 2. Tables, transcribed

Tables that are purely observational bookkeeping are summarised in one sentence each, with row counts, as agreed: Table 1 (p8, observatory and telescope coordinates; 20 rows from
CAL, 807, C60, ESO, ES2, ES1, IRT, LAS, LAV, 688, 711, 414, TEE, 675, 586, PI1, SAA, 413, ANU, UKI); Table 7 (p31, fitted station offset times; 20 rows, events U12 to Vgr2 beta Per,
values from -8.769 +/- 0.262 s for U36A IRTF to +3.717 +/- 0.009 s for U14 Pic du Midi); Table 8 (p31-34, COR RMS residuals by event and observatory; 30 events, per-event RMS 0.028 to 0.476 km);
Table 10 (p37, eight candidate lambda-ring events); Table 12 (p39, fitted corrections to proper motions, star positions and planet ephemeris offsets; 30 star rows, sky-plane
ephemeris offsets f_0, g_0 all under about 200 km, mostly under 100 km); Table 11 (p38, skyplane positions of three secondary stars, not transcribed); Table A1 (p78, 18 ring-width
comparison rows). Tables 5, 6 and 14 to 16, 18 and 20 are ring-orbit and satellite-mean-motion results and are transcribed because they carry the radii, extents and mean motions.

### Table 2 (p9): Major satellite orbits and masses from Jacobson (2023)

| Satellite | a (km) | GM_sat (km^3 s^-2) | M_sat / M_Ur x 10^-6 | M_sat (x 10^20 kg) |
| --- | --- | --- | --- | --- |
| Ariel | 190928. | 82.30 +/- 1.20 | 14.204 +/- 0.207 | 12.331 +/- 0.180 |
| Umbriel | 265981. | 86.00 +/- 1.50 | 14.843 +/- 0.259 | 12.885 +/- 0.225 |
| Titania | 436283. | 230.60 +/- 3.40 | 39.800 +/- 0.587 | 34.550 +/- 0.509 |
| Oberon | 583447. | 207.60 +/- 5.00 | 35.830 +/- 0.863 | 31.104 +/- 0.749 |
| Miranda | 129828. | 4.20 +/- 0.20 | 0.7249 +/- 0.0345 | 0.6293 +/- 0.0300 |
| Puck | 86004. | 0.1275 +/- 0.0425 | 0.0220 +/- 0.0073 | 0.0191 +/- 0.0064 |

READ (p9 and footnote 3): "The orbital semimajor axes a and masses used to compute the secular precession rate due to the major satellites are given in Table 2, from Jacobson (2023)." The
paper does not define "a" (mean or osculating; epoch; element convention) and gives no epoch for Table 2; INFERRED: they are mean semi-major axes of the ura178 fit, as in Jacobson (2014)
Table 2, but that is not stated here. The footnote: "Since satellite masses are variously quoted in the literature in units of GM_sat, fractional planet mass M_sat/M_Ur, and kg, we list all three."
The table order (Ariel, Umbriel, Titania, Oberon, Miranda, Puck) is the printed order. Puck's GM (0.1275) is a fitted value here; in Jacobson (2014) Puck was massless.

### Table 3 (p10): Minor satellite orbits, dimensions, densities and masses

GM_sat in 10^-3 km^3 s^-2, M_sat/M_Ur in 10^-10, M_sat in 10^16 kg, V in 10^5 km^3, rho in g cm^-3. "a" is from Jacobson (1998); dimensions from Karkoschka (2001a,b); square brackets mean assumed.
READ note: "Except as noted, masses and uncertainties computed from satellite shapes from Table V of Karkoschka (2001b) and assumed density rho = 0.90 gm cm^-3." Rows marked (b) are the masses inferred in this paper
from normal-mode amplitudes.

| Satellite | a (km) | sqrt(AB) (km) | B/A | A x B (km x km) | V (x10^5 km^3) | sigma(V)/V | rho (g cm^-3) | GM_sat (x10^-3 km^3 s^-2) | M_sat/M_Ur (x10^-10) | M_sat (x10^16 kg) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cordelia | 49752.000 | 21 +/- 3 | 0.7 +/- 0.2 | 25x18 | 0.339 | 0.349 | [0.90] | 2.04 +/- 0.71 | 3.52 +/- 1.23 | 3.05 +/- 1.06 |
| Cordelia (b) | | | | | | | 1.79 +0.97/-0.49 | 4.06 +/- 0.38 | 7.00 +/- 0.66 | 6.08 +/- 0.57 |
| Ophelia | 53764.000 | 23 +/- 4 | 0.7 +/- 0.3 | 27x19 | 0.408 | 0.504 | [0.90] | 2.45 +/- 1.24 | 4.23 +/- 2.13 | 3.67 +/- 1.85 |
| Ophelia (b) | | | | | | | 0.87 +0.89/-0.30 | 2.38 +/- 0.22 | 4.11 +/- 0.37 | 3.57 +/- 0.32 |
| Bianca | 59165.000 | 27 +/- 2 | 0.7 +/- 0.2 | 32x23 | 0.709 | 0.299 | [0.90] | 4.26 +/- 1.27 | 7.35 +/- 2.20 | 6.38 +/- 1.91 |
| Cressida | 61767.000 | 41 +/- 2 | 0.8 +/- 0.3 | 46x37 | 2.638 | 0.380 | [0.90] | 15.85 +/- 6.02 | 27.35 +/- 10.40 | 23.74 +/- 9.03 |
| Cressida (b) | | | | | | | 0.70 +0.44/-0.21 | 12.27 +/- 1.41 | 21.18 +/- 2.44 | 18.39 +/- 2.12 |
| Desdemona | 62659.000 | 35 +/- 4 | 0.6 +/- 0.2 | 45x27 | 1.374 | 0.375 | [0.90] | 8.25 +/- 3.09 | 14.25 +/- 5.34 | 12.37 +/- 4.63 |
| Juliet | 64358.000 | 53 +/- 4 | 0.5 +/- 0.1 | 75x37 | 4.301 | 0.230 | [0.90] | 25.83 +/- 5.95 | 44.59 +/- 10.26 | 38.71 +/- 8.91 |
| Portia | 66097.000 | 70 +/- 4 | 0.8 +/- 0.1 | 78x63 | 12.968 | 0.148 | [0.90] | 77.90 +/- 11.55 | 134.44 +/- 19.93 | 116.71 +/- 17.30 |
| Rosalind | 69927.000 | 36 +/- 6 | 1.0 +/- 0.2 | 36x36 | 1.954 | 0.314 | [0.90] | 11.74 +/- 3.69 | 20.26 +/- 6.36 | 17.59 +/- 5.52 |
| Belinda | 75255.000 | 45 +/- 8 | 0.5 +/- 0.1 | 64x32 | 2.745 | 0.327 | [0.90] | 16.49 +/- 5.39 | 28.46 +/- 9.30 | 24.71 +/- 8.07 |
| Puck | 86004.000 | 81 +/- 2 | 1.0 +/- 0.1 | 81x81 | 22.261 | 0.078 | [0.90] | 133.72 +/- 10.47 | 230.79 +/- 18.06 | 200.35 +/- 15.68 |

READ (section 9.2, p70): the masses of Cressida, Cordelia and Ophelia are the "(b)" rows, from the amplitudes of the forced modes, eq. (40) `A_m = 2 alpha a^2 (M_sat/M_Ur) |f_d| / (3 (m - 1) |Delta a_P|)` (Murray and Dermott 1999
eq. 10.22), valid "strictly" only in the evanescent region; "probably best thought of as an approximate scaling estimate". The same masses are tabulated in Table 19 (section 2 below). The paper's own masses of these three
moons are therefore: Cressida GM = 12.27 +/- 1.41 x 10^-3, Cordelia 4.06 +/- 0.38 x 10^-3, Ophelia 2.38 +/- 0.22 x 10^-3 km^3 s^-2 (kg: 18.39, 6.08, 3.57 x 10^16). The (b) values are READ; the units conversion from the
M_sat/M_Ur column printed in Table 19 is consistent: Cressida 21.18e-10 x 5793950.3 = 12.27e-3 (COMPUTED).

### Table 4 (p12): Spice kernels

| File name | Description |
| --- | --- |
| pleph.ura178.bsp | Planetary ephemerides, constrained by ring data in Paper 1 |
| ura178.bsp | Major Uranus satellite ephemerides, constrained by ring data in Paper 1 |
| vgr2.ura178.bsp | Voyager 2 trajectory at Uranus, constrained by ring data in Paper 1 |
| ura115.bsp | Minor Uranus satellite ephemerides |
| earthstns_itrf93_040916.bsp | Geocentric coordinates of DSN groundstations |
| earth_720101_070426.bpc | Earth rotation model |
| ObsCodes_Uranus_20220212.spk | Custom geocentric observer coordinates |
| pg3f0000r.bsp | HST ephemeris for U137 |
| pg490000r.bsp | HST ephemeris for U138 |
| urkao_v1.bsp | Custom ephemeris for KAO (U0) |
| naif0012.tls | leap seconds |

READ (Table 4 note): "All kernels listed are available from ftp://ssd.jpl.nasa.gov/pub/eph or as part of the Uranus ring occultation support bundle on PDS. See Appendix A of French et al. (2023b)." Underscores
in the file names were lost in the text layer; they are restored here from the page image. The text of section 3.3 (p11) calls the planetary kernel `peph.ura178.bsp`, the table `pleph.ura178.bsp`: a
typographic inconsistency inside the paper (READ, both; unresolved). READ (section 3.3): "we used the recent JPL ura178.bsp Uranus satellite ephemerides, the peph.ura178.bsp planetary ephemerides, and the associated Voyager 2
ephemeris vgr2.ura178.bsp at Uranus. These ephemerides were updated from the ura111.bsp, ura116.bsp, vgr2.ura111.bsp, and de440.bsp series (Jacobson 2014; Park et al. 2021), using constraints provided by the Uranus ring occultation data in Paper 1,
excluding the gamma and lambda rings (Jacobson 2023)." The minor-satellite kernel is ura115 (Table 4). The Table 19 note (b) says the epoch mean longitudes of the satellites come "from the ura111 ephemeris" while Fig. 37 and p69 compare with "ura115" values for Ophelia and Cordelia: another internal inconsistency (READ, unresolved).

### Table 5 (p26-28): Uranian ring keplerian orbital elements

Epoch TDB 1986 Jan 19 12:00 (note b). Columns: a; ae and a sin i; varpi_0 and Omega_0; varpi_dot and Omega_dot (deg/day); Delta varpi_dot and Delta Omega_dot (the difference between observed and predicted rates for the fitted a, with
the predicted rates including the planet field and satellites; deg/day); Delta a_varpi_dot and Delta a_Omega_dot (km; the corresponding radial offset, eqs. 9 and 10). N is the number of ring event times and RMS is the post-fit RMS residual (km), printed under each a. Eccentricity and inclination in square brackets
were held fixed during the orbit determination (note c). Each ring feature is two lines in the source (a/ae/varpi/varpi-dot/Delta varpi-dot/Delta a_varpi-dot on the first; N/RMS/a sin i/Omega/Omega-dot/Delta Omega-dot/Delta a_Omega-dot on the second); they are shown here as two rows.
"COR: center of ring. IER: inner edge of ring, OER: outer edge of ring."

| Ring | Feat. | a (km) | ae (km) | varpi_0 (deg) | varpi_dot (deg/d) | Delta varpi_dot (deg/d) | Delta a_varpi_dot (km) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| | | N, RMS (km) | a sin i (km) | Omega_0 (deg) | Omega_dot (deg/d) | Delta Omega_dot (deg/d) | Delta a_Omega_dot (km) |
| 6 | IER | 41835.920 +/- 0.102 | 42.068 +/- 0.157 | 181.733 +/- 0.213 | 2.7620198 +/- 0.0001046 | -0.0003140 | 1.352 +/- 0.450 |
| | | 45, 0.830 | 44.747 +/- 0.320 | 90.465 +/- 0.498 | -2.7569020 +/- 0.0001882 | -0.0000256 | -0.110 +/- 0.813 |
| 6 | COR | 41837.092 +/- 0.096 | 42.499 +/- 0.077 | 181.657 +/- 0.112 | 2.7619579 +/- 0.0000516 | -0.0001040 | 0.448 +/- 0.222 |
| | | 50, 0.263 | 44.643 +/- 0.104 | 89.727 +/- 0.272 | -2.7565287 +/- 0.0000502 | 0.0000766 | 0.331 +/- 0.217 |
| 6 | OER | 41838.237 +/- 0.102 | 42.930 +/- 0.158 | 181.526 +/- 0.209 | 2.7619237 +/- 0.0001028 | 0.0001277 | -0.550 +/- 0.442 |
| | | 45, 0.603 | 44.446 +/- 0.323 | 88.850 +/- 0.507 | -2.7560947 +/- 0.0001898 | 0.0002455 | 1.061 +/- 0.820 |
| 5 | IER | 42233.577 +/- 0.089 | 80.052 +/- 0.148 | 176.823 +/- 0.089 | 2.6715990 +/- 0.0000451 | -0.0003241 | 1.456 +/- 0.203 |
| | | 59, 0.603 | 40.572 +/- 0.270 | 297.500 +/- 0.482 | -2.6667586 +/- 0.0001679 | -0.0000154 | -0.069 +/- 0.757 |
| 5 | COR | 42234.893 +/- 0.091 | 80.237 +/- 0.076 | 176.819 +/- 0.048 | 2.6715756 +/- 0.0000227 | -0.0000548 | 0.246 +/- 0.102 |
| | | 67, 0.212 | 40.951 +/- 0.114 | 296.068 +/- 0.259 | -2.6663516 +/- 0.0000445 | 0.0000998 | 0.450 +/- 0.200 |
| 5 | OER | 42236.261 +/- 0.089 | 80.463 +/- 0.148 | 176.851 +/- 0.089 | 2.6715450 +/- 0.0000449 | 0.0002192 | -0.985 +/- 0.202 |
| | | 59, 0.641 | 41.331 +/- 0.270 | 294.460 +/- 0.475 | -2.6658983 +/- 0.0001677 | 0.0002496 | 1.125 +/- 0.756 |
| 4 | IER | 42569.511 +/- 0.089 | 44.931 +/- 0.146 | 256.398 +/- 0.183 | 2.5978170 +/- 0.0000927 | -0.0006562 | 3.057 +/- 0.431 |
| | | 57, 0.744 | 23.183 +/- 0.313 | 339.064 +/- 1.489 | -2.5940762 +/- 0.0004674 | -0.0005631 | -2.630 +/- 2.183 |
| 4 | COR | 42571.124 +/- 0.091 | 45.347 +/- 0.076 | 255.880 +/- 0.093 | 2.5980293 +/- 0.0000433 | -0.0000977 | 0.455 +/- 0.202 |
| | | 63, 0.264 | 23.413 +/- 0.076 | 336.828 +/- 0.671 | -2.5932795 +/- 0.0001121 | -0.0001117 | -0.522 +/- 0.524 |
| 4 | OER | 42572.743 +/- 0.089 | 45.742 +/- 0.143 | 255.294 +/- 0.181 | 2.5981594 +/- 0.0000908 | 0.0003798 | -1.769 +/- 0.423 |
| | | 57, 0.727 | 24.329 +/- 0.308 | 333.891 +/- 1.361 | -2.5925059 +/- 0.0004191 | 0.0003157 | 1.475 +/- 1.958 |
| alpha | IER | 44714.884 +/- 0.083 | 32.512 +/- 0.118 | 205.885 +/- 0.213 | 2.1851508 +/- 0.0000949 | -0.0007266 | 4.229 +/- 0.552 |
| | | 73, 0.542 | 11.853 +/- 0.276 | 206.062 +/- 1.839 | -2.1831510 +/- 0.0005365 | -0.0010535 | -6.145 +/- 3.130 |
| alpha | COR | 44718.473 +/- 0.086 | 33.916 +/- 0.062 | 206.074 +/- 0.112 | 2.1855003 +/- 0.0000510 | 0.0002392 | -1.392 +/- 0.297 |
| | | 81, 0.249 | 12.005 +/- 0.089 | 203.837 +/- 0.894 | -2.1813455 +/- 0.0001311 | 0.0001373 | 0.801 +/- 0.765 |
| alpha | OER | 44722.161 +/- 0.084 | 35.335 +/- 0.118 | 206.243 +/- 0.195 | 2.1857991 +/- 0.0000898 | 0.0011714 | -6.819 +/- 0.523 |
| | | 73, 0.543 | 11.676 +/- 0.289 | 202.596 +/- 1.749 | -2.1799630 +/- 0.0005158 | 0.0008880 | 5.186 +/- 3.012 |
| beta | IER | 45656.770 +/- 0.087 | 18.573 +/- 0.128 | 319.589 +/- 0.435 | 2.0309488 +/- 0.0001626 | -0.0004852 | 3.103 +/- 1.040 |
| | | 71, 0.574 | 3.799 +/- 0.194 | 224.376 +/- 8.293 | -2.0272655 +/- 0.0019207 | 0.0007997 | 5.128 +/- 12.315 |
| beta | COR | 45661.056 +/- 0.087 | 20.106 +/- 0.071 | 317.925 +/- 0.222 | 2.0308229 +/- 0.0000857 | 0.0000590 | -0.377 +/- 0.548 |
| | | 78, 0.266 | 4.004 +/- 0.070 | 232.878 +/- 4.169 | -2.0279790 +/- 0.0006048 | -0.0005821 | -3.733 +/- 3.880 |
| beta | OER | 45665.286 +/- 0.089 | 21.745 +/- 0.131 | 316.622 +/- 0.383 | 2.0308319 +/- 0.0001456 | 0.0007290 | -4.664 +/- 0.932 |
| | | 71, 0.468 | 3.851 +/- 0.206 | 239.903 +/- 7.558 | -2.0287602 +/- 0.0016596 | -0.0020226 | -12.970 +/- 10.649 |
| eta | IER | 47174.853 +/- 0.093 | [0.00] | | | | |
| | | 56, 0.679 | [0.00] | | | | |
| eta | COR | 47176.009 +/- 0.088 | [0.00] | | | | |
| | | 60, 0.297 | [0.00] | | | | |
| eta | OER | 47177.080 +/- 0.093 | [0.00] | | | | |
| | | 55, 0.525 | [0.00] | | | | |

gamma ring, including the m = 3 normal mode:

| Ring | Feat. | a (km) | ae (km) | varpi_0 (deg) | varpi_dot (deg/d) | Delta varpi_dot (deg/d) | Delta a_varpi_dot (km) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gamma | IER | 47624.606 +/- 0.103 | 5.016 +/- 0.122 | 47.274 +/- 1.570 | 1.7516631 +/- 0.0006611 | 0.0001409 | -1.091 +/- 5.119 |
| | | 76, 0.536 | [0.00] | | | | |
| gamma | COR | 47626.170 +/- 0.089 | 5.306 +/- 0.071 | 48.720 +/- 0.858 | 1.7514650 +/- 0.0004210 | 0.0001448 | -1.121 +/- 3.260 |
| | | 83, 0.401 | [0.00] | | | | |
| gamma | OER | 47627.865 +/- 0.100 | 5.470 +/- 0.130 | 47.145 +/- 1.458 | 1.7524301 +/- 0.0005903 | 0.0013286 | -10.286 +/- 4.572 |
| | | 76, 0.692 | [0.00] | | | | |

gamma ring, excluding the m = 3 normal mode:

| Ring | Feat. | a (km) | ae (km) | varpi_0 (deg) | varpi_dot (deg/d) | Delta varpi_dot (deg/d) | Delta a_varpi_dot (km) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gamma | IER | 47624.606 +/- 0.129 | 5.016 +/- 0.153 | 47.274 +/- 1.964 | 1.7516631 +/- 0.0008271 | 0.0001409 | -1.091 +/- 6.405 |
| | | 76, 0.536 | [0.00] | | | | |
| gamma | COR | 47626.289 +/- 0.094 | 5.314 +/- 0.074 | 47.573 +/- 0.881 | 1.7524471 +/- 0.0003578 | 0.0011422 | -8.842 +/- 2.771 |
| | | 83, 0.490 | [0.00] | | | | |
| gamma | OER | 47627.987 +/- 0.122 | 5.524 +/- 0.157 | 46.394 +/- 1.775 | 1.7536847 +/- 0.0006681 | 0.0025990 | -20.112 +/- 5.175 |
| | | 76, 0.895 | [0.00] | | | | |

| Ring | Feat. | a (km) | N, RMS (km) | ae, a sin i (km) |
| --- | --- | --- | --- | --- |
| delta | IER | 48297.775 +/- 0.082 | 73, 0.608 | [0.00], [0.00] |
| delta | COR | 48300.227 +/- 0.082 | 80, 0.273 | [0.00], [0.00] |
| delta | OER | 48302.752 +/- 0.083 | 72, 0.683 | [0.00], [0.00] |
| lambda | COR | 50026.557 +/- 1.314 | 8, 3.477 | [0.00], [0.00] |

| Ring | Feat. | a (km) | ae (km) | varpi_0 (deg) | varpi_dot (deg/d) | Delta varpi_dot (deg/d) | Delta a_varpi_dot (km) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| epsilon | IER | 51120.014 +/- 0.077 | 386.530 +/- 0.107 | 307.089 +/- 0.016 | 1.3632690 +/- 0.0000070 | -0.0028286 | 30.128 +/- 0.075 |
| | | 79, 0.583 | [0.00] | | | | |
| epsilon | COR | 51149.279 +/- 0.081 | 405.894 +/- 0.062 | 307.076 +/- 0.009 | 1.3632575 +/- 0.0000042 | -0.0001081 | 1.154 +/- 0.045 |
| | | 89, 0.396 | [0.00] | | | | |
| epsilon | OER | 51178.588 +/- 0.079 | 425.242 +/- 0.111 | 307.061 +/- 0.015 | 1.3632618 +/- 0.0000070 | 0.0026240 | -28.026 +/- 0.075 |
| | | 77, 0.560 | [0.00] | | | | |

The eta, delta (all three features) and lambda rows carry no eccentricity, inclination, apse or node because e and a sin i were held at zero. The gamma COR "excluding m = 3" entry is the alternate fit; the
nominal gamma fit is the "including m = 3" one (READ, Table 5 and p20). The IER/OER Delta Omega_dot and Delta a_Omega_dot entries sit in the second line of each feature, as printed.

### Table 6 (p29-30): Uranian ring normal modes

A_m amplitude (km), delta_m phase (deg), Omega_P pattern speed (deg/d), Delta Omega_P (observed minus predicted, deg/d), Delta a_P (km). Note (b): longitudes reduced to the minimum value of lambda mod 360/|m|; epoch TDB 1986 Jan 19 12:00.

| Ring | Feat. | m | A_m (km) | delta_m (deg) | Omega_P (deg/d) | Delta Omega_P (deg/d) | Delta a_P (km) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| eta | IER | 3 | 0.452 +/- 0.122 | 77.208 +/- 5.798 | 776.585539 +/- 0.002490 | 0.086574 | -3.499 +/- 0.101 |
| eta | COR | 3 | 0.600 +/- 0.069 | 75.901 +/- 2.359 | 776.584048 +/- 0.001164 | 0.113697 | -4.595 +/- 0.047 |
| eta | OER | 3 | 0.677 +/- 0.124 | 75.946 +/- 3.819 | 776.584059 +/- 0.001818 | 0.140201 | -5.666 +/- 0.073 |

gamma ring, including m = 3:

| Feat. | m | A_m (km) | delta_m (deg) | Omega_P (deg/d) | Delta Omega_P (deg/d) | Delta a_P (km) |
| --- | --- | --- | --- | --- | --- | --- |
| IER | 0 | 6.137 +/- 0.138 | 309.462 +/- 1.178 | 1145.576306 +/- 0.000538 | -0.050871 | 1.411 +/- 0.015 |
| IER | 6 | 0.590 +/- 0.110 | 33.405 +/- 2.208 | 956.418079 +/- 0.000877 | -0.022759 | 0.754 +/- 0.029 |
| IER | -2 | 1.099 +/- 0.146 | 35.175 +/- 2.940 | 1720.120475 +/- 0.001212 | -0.071814 | 1.325 +/- 0.022 |
| IER | -1 | 1.656 +/- 0.120 | 31.309 +/- 4.318 | 2292.907872 +/- 0.002002 | -0.098006 | 1.357 +/- 0.028 |
| COR | 0 | 5.509 +/- 0.076 | 311.695 +/- 0.750 | 1145.576845 +/- 0.000352 | 0.006033 | -0.167 +/- 0.010 |
| COR | 3 | 0.577 +/- 0.073 | 54.989 +/- 2.254 | 765.399832 +/- 0.001152 | -0.065364 | 2.706 +/- 0.048 |
| COR | 6 | 0.637 +/- 0.063 | 32.518 +/- 1.171 | 956.419622 +/- 0.000474 | 0.025958 | -0.861 +/- 0.016 |
| COR | -2 | 0.690 +/- 0.079 | 26.843 +/- 2.762 | 1720.117977 +/- 0.001217 | 0.010437 | -0.193 +/- 0.022 |
| COR | -1 | 1.822 +/- 0.067 | 25.852 +/- 2.169 | 2292.904186 +/- 0.001054 | 0.011240 | -0.156 +/- 0.015 |
| OER | 0 | 4.702 +/- 0.128 | 313.066 +/- 1.574 | 1145.576675 +/- 0.000698 | 0.066941 | -1.857 +/- 0.019 |
| OER | 3 | 0.935 +/- 0.129 | 56.867 +/- 2.520 | 765.397669 +/- 0.001028 | -0.026590 | 1.101 +/- 0.043 |
| OER | 6 | 0.892 +/- 0.117 | 28.027 +/- 1.363 | 956.419529 +/- 0.000560 | 0.076982 | -2.552 +/- 0.019 |
| OER | -1 | 1.955 +/- 0.122 | 25.258 +/- 3.578 | 2292.903174 +/- 0.001739 | 0.132603 | -1.836 +/- 0.024 |

gamma ring, excluding m = 3:

| Feat. | m | A_m (km) | delta_m (deg) | Omega_P (deg/d) | Delta Omega_P (deg/d) | Delta a_P (km) |
| --- | --- | --- | --- | --- | --- | --- |
| IER | 0 | 6.137 +/- 0.172 | 309.462 +/- 1.474 | 1145.576306 +/- 0.000673 | -0.050871 | 1.411 +/- 0.019 |
| IER | 6 | 0.590 +/- 0.138 | 33.405 +/- 2.763 | 956.418079 +/- 0.001097 | -0.022759 | 0.754 +/- 0.036 |
| IER | -2 | 1.099 +/- 0.183 | 35.175 +/- 3.678 | 1720.120475 +/- 0.001517 | -0.071814 | 1.325 +/- 0.028 |
| IER | -1 | 1.656 +/- 0.150 | 31.309 +/- 5.403 | 2292.907872 +/- 0.002505 | -0.098006 | 1.357 +/- 0.035 |
| COR | 0 | 5.450 +/- 0.079 | 311.233 +/- 0.762 | 1145.576581 +/- 0.000376 | 0.010044 | -0.279 +/- 0.010 |
| COR | 6 | 0.698 +/- 0.066 | 30.504 +/- 1.075 | 956.419653 +/- 0.000442 | 0.029566 | -0.980 +/- 0.015 |
| COR | -2 | 0.635 +/- 0.083 | 27.010 +/- 3.104 | 1720.120139 +/- 0.001183 | 0.019028 | -0.351 +/- 0.022 |
| COR | -1 | 1.807 +/- 0.070 | 26.105 +/- 2.325 | 2292.903886 +/- 0.001134 | 0.019506 | -0.270 +/- 0.016 |
| OER | 0 | 4.724 +/- 0.160 | 313.481 +/- 1.841 | 1145.575881 +/- 0.000844 | 0.070553 | -1.958 +/- 0.023 |
| OER | 6 | 1.057 +/- 0.140 | 26.147 +/- 1.358 | 956.419580 +/- 0.000563 | 0.080719 | -2.676 +/- 0.019 |
| OER | -1 | 1.966 +/- 0.151 | 22.643 +/- 4.393 | 2292.903576 +/- 0.002122 | 0.141832 | -1.964 +/- 0.029 |

delta and epsilon rings:

| Ring | Feat. | m | A_m (km) | delta_m (deg) | Omega_P (deg/d) | Delta Omega_P (deg/d) | Delta a_P (km) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| delta | IER | 2 | 2.349 +/- 0.105 | 171.686 +/- 1.527 | 562.516712 +/- 0.000580 | -0.042324 | 2.415 +/- 0.033 |
| delta | COR | 2 | 3.169 +/- 0.059 | 170.383 +/- 0.666 | 562.516205 +/- 0.000258 | 0.000146 | -0.008 +/- 0.015 |
| delta | COR | 23 | 0.339 +/- 0.063 | 6.844 +/- 0.472 | 1074.523021 +/- 0.000189 | -0.072580 | 2.173 +/- 0.006 |
| delta | OER | 2 | 4.035 +/- 0.105 | 171.248 +/- 0.885 | 562.516313 +/- 0.000351 | 0.044481 | -2.539 +/- 0.020 |
| epsilon | IER | -24 | 1.011 +/- 0.107 | 2.916 +/- 0.259 | 1074.522889 +/- 0.000099 | -0.034140 | 1.082 +/- 0.003 |
| epsilon | COR | 14 | 0.383 +/- 0.071 | 14.105 +/- 0.628 | 956.418015 +/- 0.000269 | -0.798192 | 28.425 +/- 0.010 |
| epsilon | COR | -24 | 0.443 +/- 0.061 | 2.828 +/- 0.346 | 1074.522703 +/- 0.000142 | 0.888501 | -28.177 +/- 0.005 |
| epsilon | OER | 14 | 0.590 +/- 0.130 | 14.904 +/- 0.704 | 956.418119 +/- 0.000298 | 0.024894 | -0.887 +/- 0.011 |

READ: no normal modes were detected for rings 6, 5, 4, alpha or beta (abstract; section 4.1).

### Table 9 (p35): 2-sigma eccentricity and inclination limits

| Ring | ae (km) | a sin i (km) |
| --- | --- | --- |
| eta | 0.126 | 0.290 |
| gamma | - | 0.286 |
| delta | 0.088 | 0.284 |
| epsilon | - | 0.158 |

### Table 13 (p41): Uranus pole direction and ring plane radius scale (every fit)

alpha_P and delta_P are the ICRF/J2000 right ascension and declination of the pole (the positive-angular-momentum direction, section 3.2). rho is the correlation coefficient rho(alpha_P, delta_P).
Delta alpha_P and Delta delta_P are differences from Fit 1. <Delta a> is the mean offset of the fitted ring semi-major axes from Fit 1 and sigma(Delta a) its standard deviation.

| Fit | alpha_P (deg) | delta_P (deg) | rho | Delta alpha_P (deg) | Delta delta_P (deg) | <Delta a> (km) | sigma(Delta a) (km) | Description |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 77.311327 +/- 0.000128 | 15.172795 +/- 0.000472 | -0.22 | - | - | - | - | nominal pole direction |
| 2 | 77.311327 +/- 0.000141 | 15.172795 +/- 0.000618 | -0.42 | - | - | - | - | w/ Vgr trajectory errors |
| 3 | 77.311210 +/- 0.000128 | 15.172762 +/- 0.000471 | -0.22 | -0.000117 | -0.000033 | -0.001 | 0.003 | ura178 pole |
| 4 | 77.311246 +/- 0.000274 | 15.173393 +/- 0.002013 | -0.88 | -0.000081 | 0.000597 | 0.278 | 0.027 | Earth-based only |
| 5 | 77.311200 +/- 0.000400 | 15.172400 +/- 0.001700 | -0.39 | -0.000127 | -0.000395 | - | - | Jacobson (2023) |
| 6 | 77.310000 +/- 0.003000 | 15.172000 +/- 0.002000 | [1.0] | -0.001327 | -0.000795 | -0.015 | 0.091 | Jacobson (2014) |
| 7 | 77.310877 +/- 0.003400 | 15.174564 +/- 0.003300 | [1.0] | -0.000450 | 0.001768 | -0.055 | 0.092 | French et al. (1991) |

READ: the Fit 6 alpha error is printed 0.003000, the paper's text on p40 says "(alpha = 77.310 +/- 0.002 deg and delta = 15.172 +/- 0.002 deg)" and Jacobson (2014) Table 12 prints +/- 0.002: an internal discrepancy of the paper in the alpha error of the Jacobson (2014) entry. The rho of 1.0 for Fits 6 and 7 is an assumed
value in square brackets ("assumed correlation rho = 0" in the text of p40, the table prints [1.0]; another internal inconsistency, unresolved).

### Table 14 (p45): Width-radius results

| Ring | a (km) | mean width (km) | m | A_m (km) | a delta e_m (km) | delta varpi_0m (deg) | q_em | q_varpi m | dW_m/dr | delta Omega_P (deg/yr) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6 | 41837.092 | 2.316 | 1 | 42.499 | 0.862 | -0.207 | 0.372 | 0.077 | 0.020 | -0.035 |
| 5 | 42234.893 | 2.684 | 1 | 80.237 | 0.412 | 0.028 | 0.153 | -0.004 | 0.005 | -0.020 |
| 4 | 42571.124 | 3.231 | 1 | 45.347 | 0.810 | -1.103 | 0.251 | 0.277 | 0.018 | 0.125 |
| alpha | 44718.473 | 7.277 | 1 | 33.916 | 2.823 | 0.359 | 0.388 | -0.139 | 0.083 | 0.237 |
| beta | 45661.056 | 8.516 | 1 | 20.106 | 3.172 | -2.967 | 0.372 | 1.105 | 0.158 | -0.043 |
| eta | 47176.009 | 2.227 | 3 | 0.600 | 0.224 | -1.262 | 0.101 | 0.127 | 0.374 | -0.540 |
| gamma | 47626.170 | 3.258 | 1 | 5.306 | 0.454 | -0.129 | 0.139 | 0.018 | 0.086 | 0.280 |
| gamma | | | 0 | 5.509 | -1.435 | 3.605 | -0.440 | 1.587 | -0.261 | 0.135 |
| gamma | | | 6 | 0.637 | 0.302 | -5.378 | 0.093 | 0.499 | 0.475 | 0.530 |
| gamma | | | -1 | 1.822 | 0.299 | -6.051 | 0.092 | 0.556 | 0.164 | - |
| gamma | | | -2 | 0.690 | -1.099 | - | -0.337 | - | -1.592 | 4.742 |
| gamma | | | 3 | 0.577 | 0.935 | - | 0.287 | - | 1.622 | -0.790 |
| delta | 48300.227 | 4.977 | 2 | 3.169 | 1.686 | -0.438 | 0.339 | 0.148 | 0.532 | -0.146 |
| delta | | | 23 | 0.339 | - | - | - | - | - | - |
| epsilon | 51149.279 | 58.574 | 1 | 405.894 | 38.712 | -0.028 | 0.661 | -0.019 | 0.095 | -0.003 |

### Table 15 (p53): Ring masses, surface densities, libration periods

| Ring | M_SSG (x10^14 kg) | Sigma_SSG (g cm^-2) | P_lib (yr) | M_CSG (x10^14 kg) | Sigma_CSG (g cm^-2) |
| --- | --- | --- | --- | --- | --- |
| 6 | 0.04 | 0.58 | 37.6 | 6.78 | 111.37 |
| 5 | 0.23 | 3.20 | 8.5 | 7.51 | 105.47 |
| 4 | 0.11 | 1.23 | 25.6 | 8.13 | 94.05 |
| alpha | 0.19 | 0.93 | 66.3 | 12.46 | 60.96 |
| beta | 0.15 | 0.60 | 118.1 | 13.52 | 55.33 |
| delta (d) | 3.84 | 25.42 | 3.0 | | |
| epsilon | 32.70 | 17.37 | 17.1 | 69.89 | 37.13 |

Notes (READ): (a) N = 500 streamlines in the standard self-gravity model, collisional effects neglected; (b) from eq. (32); (c) N = 20,000 streamlines in the collisional self-gravity model of Chiang and Goldreich (2000) for dispersion velocity c_b = 1 cm/s (mass and surface density scale as c_b^1.5);
(d) the delta ring value is from the two-streamline solution of eq. (35) for m = 2. The gamma ring is excluded "because of the confounding and uncertain influence of its many normal modes" (p53). Ring masses do not enter the project's Uranian models.

### Table 16 (p57): Ring center of opacity

| Ring | a_COR (km) | Delta a_COO (km) | Delta a_COO / sigma | a_COO (km) |
| --- | --- | --- | --- | --- |
| 6 | 41837.092 +/- 0.096 | 0.040 +/- 0.193 | 0.208 | 41837.132 +/- 0.215 |
| 5 | 42234.893 +/- 0.091 | 0.169 +/- 0.100 | 1.686 | 42235.061 +/- 0.135 |
| 4 | 42571.124 +/- 0.091 | 0.404 +/- 0.087 | 4.620 | 42571.528 +/- 0.126 |
| alpha | 44718.473 +/- 0.086 | 0.172 +/- 0.132 | 1.309 | 44718.645 +/- 0.157 |
| beta | 45661.056 +/- 0.087 | 0.322 +/- 0.274 | 1.175 | 45661.378 +/- 0.288 |
| gamma | 47626.170 +/- 0.089 | 0.119 +/- 0.125 | 0.954 | 47626.289 +/- 0.153 |
| epsilon | 51149.279 +/- 0.081 | 1.134 +/- 0.488 | 2.324 | 51150.414 +/- 0.495 |

### Table 17 (p59): Uranus gravity parameters (every fit with its note)

J2, J4 and J6 are x 10^-6, at the reference radius R = 25559 km. Square brackets mean held fixed. The note column is verbatim. AUTOMP is the JPL satellite-precession routine named in the notes (not otherwise defined in the paper).

| Fit | J2 x10^-6 | J4 x10^-6 | J6 x10^-6 | rho(J2, J4) | Note |
| --- | --- | --- | --- | --- | --- |
| 1 | 3513.23 +/- 0.34 | -30.32 +/- 4.73 | - | - | French et al. (1988) (converted to R = 25559 km) |
| 2 | 3510.7 +/- 0.7 | -34.2 +/- 1.3 | 0.0 +/- 1.0 | 0.9784 | Jacobson (2014) adopted solution |
| 3 | 3510.5 +/- 1.3 | -34.4 +/- 1.3 | 0.0 +/- 1.0 | 0.9784 | Jacobson (2014) "rings only" solution |
| 4 | 3510.465 +/- 0.058 | -34.145 +/- 0.082 | 0.58 +/- 0.12 | 0.9887 | Jacobson (2023) (ura178 ephemeris) |
| 5 | 3510.464 +/- 0.066 | -34.158 +/- 0.092 | [0.58] | 0.9384 | COR - 6 satellites (AUTOMP) |
| 6 | 3510.474 +/- 0.056 | -34.135 +/- 0.081 | [0.58] | 0.9853 | COR wtd fit to apse/node rates (AUTOMP) |
| 7 | 3511.175 +/- 0.065 | -33.492 +/- 0.092 | [0.50] | 0.9384 | COR + COO offsets (Table 16) - 6 satellites (AUTOMP) |
| 8 | 3510.975 +/- 0.065 | -34.030 +/- 0.092 | [0.0] | 0.9384 | Fit 7 but J6 = 0 x 10^-6 |
| 9 | 3511.374 +/- 0.065 | -32.954 +/- 0.092 | [1.0] | 0.9384 | Fit 7 but J6 = 1 x 10^-6 |
| 10 | 3511.201 +/- 0.065 | -33.528 +/- 0.092 | [0.50] | 0.9384 | Fit 7 with Delta a = +0.2 km for all rings |
| 11 | 3510.723 +/- 0.065 | -33.957 +/- 0.092 | [0.50] | 0.9384 | Fit 7 including all but Cordelia and Ophelia |
| 12 | 3509.291 +/- 0.067 | -35.522 +/- 0.094 | [0.50] | 0.9383 | Fit 7 incl. all sats, basis for Monte Carlo Fit 15 |
| 13 | 3513.217 +/- 0.065 | -31.669 +/- 0.091 | [0.50] | 0.9384 | Fit 7 excluding all satellites |
| 14 | 3509.337 +/- 0.065 | -35.438 +/- 0.091 | [0.50] | 0.9387 | COO, all satellites, rings 6,5,4, epsilon for J_n |
| **15** | **3509.291 +/- 0.412** | **-35.522 +/- 0.466** | **[0.50]** | **0.9861** | **Adopted solution: Fit 12 w/ composite error estimates** |
| 16 | 3509.291 +/- 0.385 | -35.522 +/- 0.433 | [0.50] | 0.9957 | Monte Carlo errors from COO fits to apse/node rates only |
| 17 | 3509.291 +/- 0.134 | -35.522 +/- 0.147 | [0.50] | 0.9997 | Monte Carlo errors from satellite mass uncertainties only |
| 18 | 3509.291 +/- 0.026 | -35.522 +/- 0.036 | [0.50] | 0.1412 | Monte Carlo errors from radius scale uncertainty only |

Quantities listed below the fits:

| Quantity | Value | Note |
| --- | --- | --- |
| GM_U | 5793950.300 km^3 s^-2 | Jacobson (2023) |
| R | 25559 km | Reference radius for J_n |
| dJ4/dJ2 | +1.134 | slope of error ellipse of adopted solution |
| dJ4/dJ2 | +2.69540 | slope of line connecting COO J6 = (0, 1) x 10^-6 |
| dJ2/dJ6 | +0.39909 | over range J6 = (0, 1) x 10^-6 |
| dJ4/dJ6 | +1.07570 | over range J6 = (0, 1) x 10^-6 |
| dJ2/dDelta a | +0.130 x 10^-6 km^-1 | from Fits 7 and 10 |
| dJ4/dDelta a | -0.180 x 10^-6 km^-1 | from Fits 7 and 10 |

Cross-checks (COMPUTED): Fit 2 equals Jacobson (2014) Table 12 "Current Results" (3510.7 +/- 0.7; -34.2 +/- 1.3) and Fit 3 equals its "Rings Only" column (3510.5 +/- 1.3; -34.4 +/- 1.3), so the Jacobson 2014 digest
transcription is confirmed by an independent reading. The J6 column of Fit 2 is "0.0 +/- 1.0" here, as the Jacobson (2014) Table 12 "(b) Not estimated". The numbers 3510.465 and 0.58 for Fit 4 are Jacobson (2023) as quoted by this paper; the paper that holds them is not held.
Appendix C eq. (C4) and (C7) print the adopted solution as 3509.281 and -35.533 (x10^-6), 0.010 and 0.011 away from Table 17's 3509.291 and -35.522 (READ, p85-86; unresolved: either a rounding of an earlier iteration or a typographical slip; Table 17, the abstract and the conclusions all say 3509.291 and -35.522).

### Table 18 (p63): Limits on anomalous precession rates (10^-5 deg/day)

| Ring | Delta varpi_dot | Description |
| --- | --- | --- |
| 6 | -9.07 +/- 5.31 | No detected anomalous precession |
| 5 | -1.04 +/- 2.30 | No detected anomalous precession |
| 4 | +0.02 +/- 4.48 | No detected anomalous precession |
| alpha | 27.66 +/- 4.99 | Consistent with Chancia and Hedman (2016) wake model |
| beta | 12.79 +/- 8.74 | Consistent with Chancia and Hedman (2016) wake model |
| gamma | 31.21 +/- 43.84 | m = 3 mode included in orbit model |
| gamma | 123.16 +/- 36.54 | m = 3 mode not included in orbit model |
| epsilon | - | Not determinable from orbit model |

### Table 19 (p68): Satellite-driven normal modes (the masses of Cressida, Cordelia and Ophelia from the mode amplitudes)

| Satellite | Ring | a (km) | Delta a_P (km) | m | A_m (km) | m_sat/M_Ur (x10^-10) | Omega_P (deg/d) | delta'_m (deg) | lambda_sat (deg) | Delta lambda (deg) | apse |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cressida | eta (COR) | 47176.009 | -4.595 | 3 | 0.600 +/- 0.069 | 21.18 +/- 2.44 | 776.584048 +/- 0.001164 | 15.90 +/- 2.36 | 17.44 | -1.54 +/- 2.36 | apoapse |
| Ophelia | gamma (IER) | 47624.606 | 0.754 | 6 | 0.590 +/- 0.110 | | 956.418079 +/- 0.000877 | 303.40 +/- 2.21 | 298.08 | 5.33 +/- 2.21 | apoapse |
| Ophelia | gamma (COR) | 47626.170 | -0.861 | 6 | 0.637 +/- 0.063 | 4.23 +/- 0.42 | 956.419622 +/- 0.000474 | 302.52 +/- 1.17 | 298.08 | 4.44 +/- 1.17 | apoapse |
| Ophelia | gamma (OER) | 47627.865 | -2.552 | 6 | 0.892 +/- 0.117 | | 956.419529 +/- 0.000560 | 298.03 +/- 1.36 | 298.08 | -0.05 +/- 1.36 | apoapse |
| Ophelia | epsilon (COR) | 51149.279 | 28.425 | 14 | 0.383 +/- 0.071 | | 956.418015 +/- 0.000269 | 296.96 +/- 0.63 | 298.08 | -1.12 +/- 0.63 | periapse |
| Ophelia | epsilon (OER) | 51178.588 | -0.887 | 14 | 0.590 +/- 0.130 | 3.67 +/- 0.81 | 956.418119 +/- 0.000298 | 297.76 +/- 0.70 | 298.08 | -0.32 +/- 0.70 | periapse |
| Cordelia | delta (COR) | 48300.227 | 2.173 | 23 | 0.339 +/- 0.063 | 5.58 +/- 1.09 | 1074.523021 +/- 0.000189 | 69.45 +/- 0.47 | 70.00 | -0.55 +/- 0.47 | periapse |
| Cordelia | epsilon (IER) | 51120.014 | 1.082 | -24 | 1.011 +/- 0.107 | 7.83 +/- 0.83 | 1074.522889 +/- 0.000099 | 70.42 +/- 0.26 | 70.00 | 0.41 +/- 0.26 | apoapse |
| Cordelia | epsilon (COR) | 51149.279 | -28.177 | -24 | 0.443 +/- 0.061 | | 1074.522703 +/- 0.000142 | 70.33 +/- 0.35 | 70.00 | 0.33 +/- 0.35 | apoapse |

READ notes: the epoch is 1986 Jan 19 12:00 (TDB); lambda_sat is the mean longitude of the satellite at the epoch "from the ura111 ephemeris" (note b); the Delta lambda error bars "do not include the uncertainties in the satellite ephemerides (see Table 4, Jacobson (1998))".
Table 19 gives one mass ratio per resonance (on the first line of each group for Ophelia, a single entry for the others); the final masses are the Table 3 "(b)" rows. COMPUTED check: the Table 3 (b) mass ratios are the inverse-variance weighted means of the Table 19 entries. Cordelia: (5.58 / 1.09^2 + 7.83 / 0.83^2) / (1 / 1.09^2 + 1 / 0.83^2) = 16.063 / 2.2933 = 7.004, sigma = 1 / sqrt(2.2933) = 0.660, printed 7.00 +/- 0.66. Ophelia: (4.23 / 0.42^2 + 3.67 / 0.81^2) / (1 / 0.42^2 + 1 / 0.81^2) = 29.57 / 7.193 = 4.111, sigma 0.373, printed 4.11 +/- 0.37.

### Table 20 (p68): Satellite mean motions (deg/day)

| Satellite | Source | Omega_P (deg/d) |
| --- | --- | --- |
| Cressida | This work | 776.584048 +/- 0.001164 |
| Cressida | Jacobson (1998) | 776.582414 +/- 0.000022 |
| Cressida | Showalter and Lissauer (2006) | 776.582789 +/- 0.000059 |
| Cressida | Robert French (pers. comm.) | 776.5826393 +/- 0.0000053 |
| Ophelia | This work | 956.418404 +/- 0.000171 |
| Ophelia | Jacobson (1998) | 956.42833 +/- 0.00908 |
| Ophelia | Robert French (pers. comm.) | 956.418293 +/- 0.000037 |
| Cordelia | This work | 1074.522858 +/- 0.000075 |
| Cordelia | Jacobson (1998) | 1074.51832 +/- 0.00187 |

COMPUTED periods (360 / rate): Cressida 0.463569 d (from "This work"), Ophelia 0.376404 d, Cordelia 0.335032 d. READ (p69): the weighted mean pattern speed of Ophelia's modes is 0.01 deg/day below the ura115 value, "just ~0.65 sigma from the unpublished value 956.418293".
The measured mean motions of Ophelia and Cordelia differ from the Jacobson (1998) values by -0.0099 and +0.0045 deg/day (COMPUTED: 956.418404 - 956.42833 = -0.00993; 1074.522858 - 1074.51832 = +0.00454), as the figure on p70 shows.

## 3. The adopted solution and why

READ (section 7.3.2, p62): "Our adopted solution for the Uranus gravity field is from Fit 15: J2 = (3509.291 +/- 0.412) x 10^-6, J4 = (-35.522 +/- 0.466) x 10^-6, and rho(J2, J4) = 0.9861. These represent the best-fitting values of
J2 and J4 from Fit 12, which used the entire set of ring observations, corrected for the estimated differences between COO and COR values, including the contribution to ring precession of both major and minor satellites, and assume
J6 = 0.5 x 10^-6. The formal errors sigma(J2) and sigma(J4) and correlation rho(J2, J4) for Fit 12 account for random errors of the observations alone, whereas the Monte Carlo results summarized in Fit 15 account as well for the estimated uncertainties
in Delta a_COO, Delta a, and the satellite contributions to the ring precession rates." Pole, radius scale: the adopted pole is Fits 1 and 2 (section 5.1, p40; "we adopt a fixed pole direction of Fits 1 and 2 for our preferred orbit solution, with the error bars from Fit 2"), and the radius scale
is accurate "to ~0.2 km at the 2-sigma level" (p43).

How it was fitted (READ, section 7.2, p59-62): the ring apse and node rates of Table 5 are fitted for J2 and J4 with J6 held at an assumed value. J6 cannot be fitted: "in a separate fit, not included in Table 17, we allowed J6 to be an additional free parameter, but the resulting formal
errors in all three gravity parameters were very large, the retrieved value of J6 was implausibly large, and the correlations between the parameters were nearly singular" (p60). The radial span of the data is under 10,000 km (R/a from 0.61 to 0.50), so J2, J4 and J6 are "strongly correlated"
(p55). Fit 12 is Fit 7 (COR plus the COO offsets of Table 16, J6 = 0.5e-6) with all satellites, major and minor, included in the precession. Error estimate (section 7.3.1, p62): 50,000 Monte Carlo fits that randomly displace each ring's semi-major axis by the COO error, the radius scale by 0.2 km (sigma),
and the satellite masses by their (Jacobson 2023 for the major, Table 3 for the minor moons including the paper's own masses for Cressida, Cordelia and Ophelia) uncertainties; the dominant error is the COO location (Fit 16, 0.385) then the satellite masses (Fit 17, 0.134), then the radius scale (Fit 18, 0.026). The ellipse
enclosing 1-sigma (e^-1/2 = 0.6065 of samples) sets sigma and rho.

Why it differs from the earlier fits (READ, abstract, p1, quoted): "This result differs significantly from both earlier and more recent results (Jacobson 2014, 2023), owing to our inclusion of previously neglected systematic effects, such as the offset of semimajor axes of the geometric
ring centerlines from their estimated dynamical centers of mass and the significant contributions of Cordelia and Ophelia to the precession rate of the epsilon ring." Conclusions (p75): "This result is significantly displaced from previous results (Jacobson 2014), owing to the inclusion of previously neglected systematic effects. The quoted errors account for systematic
differences (and their uncertainties) between the fitted semimajor axes of the ring centerlines (COR) and their estimated centers of opacity (COO), a surrogate for the center of mass (COM) semimajor axes that should be used when computing the radial gradient in the apse and node rates. The minor satellites Cordelia and Ophelia contribute significantly to the
precession rate of the epsilon ring, affecting the final solution for J2 and J4." Section 7.3.2 adds: "A striking characteristic of the results in Table 17 ... is the small size of formal error ellipses in J2 and J4 from the individual fits, compared to the much larger differences in J2 and J4 that stem from a range of assumptions about the value of J6, the satellites included in the
calculation of secular precession, the choice of COR or COO for ring semimajor axes, the absolute radius scale, and realistic uncertainties in these quantities."

Decomposition of the J2 differences (COMPUTED from Table 17; the paper does not print this accounting, and the steps are not strictly additive since J6 differs between Fits 5 and 7):

| Step | Fits | J2 before | J2 after | Change (x10^-6) |
| --- | --- | --- | --- | --- |
| Jacobson (2023) COR fit, 6 satellites, J6 = 0.58 | 4 vs 5 | 3510.465 | 3510.464 | -0.001 |
| Use COO instead of COR semi-major axes (and J6 0.58 to 0.50) | 5 to 7 | 3510.464 | 3511.175 | +0.711 |
| Add all satellites except Cordelia and Ophelia | 7 to 11 | 3511.175 | 3510.723 | -0.452 |
| Add Cordelia and Ophelia | 11 to 12 | 3510.723 | 3509.291 | -1.432 |
| Net, Fit 5 to Fit 12 | | 3510.464 | 3509.291 | -1.173 |
| Include no satellites at all (Fit 13 vs Fit 12) | 12 vs 13 | 3509.291 | 3513.217 | +3.926 |
| J6 from 0 to 1 x10^-6 (Fit 8 to Fit 9) | | 3510.975 | 3511.374 | +0.399 (= dJ2/dJ6) |
| Radius scale +0.2 km (Fit 7 to 10) | | 3511.175 | 3511.201 | +0.026 (= 0.2 x 0.130) |

So the move from Jacobson (2014) or (2023) to French et al. is mostly the minor-satellite masses (-1.884 from all added satellites, of which Cordelia and Ophelia alone -1.432), partly offset by the COO correction (+0.711). READ (section 7.2.6, p61): "Previous published determinations of J2 and J4
have taken into account the secular precession due to Ariel, Umbriel, Titania, Oberon, Miranda, and Puck, but have ignored the contributions of the smaller moonlets". Fit 12's J2 uses the paper's own Cordelia and Ophelia masses (Table 3, rows b).

Dependence on the assumed satellite masses and on J6, as asked:
- Masses: READ Fit 17 (satellite mass uncertainties alone) gives sigma(J2) = 0.134e-6, sigma(J4) = 0.147e-6, rho = 0.9997; so the mass uncertainties are the second-largest error source. The shift from "no satellites" to "all satellites" is 3.926e-6 in J2 (COMPUTED above),
  about 9.5 sigma of the adopted 0.412, i.e. the satellite contribution to the ring precession is not negligible and the satellite masses are part of the J2 definition. J2 as adopted is therefore conditional on the Jacobson (2023) Table 2 GMs and on the Table 3 minor-satellite GMs; INFERRED: a propagator that
  pairs this J2 with the Jacobson (2014) moon GMs (as the project does) is mixing sets. The moon-GM differences between 2014 and 2023 (at most 1.1 sigma each, section 5.1) are of the size of the mass uncertainties that Fit 17 propagates, so the mismatch maps to a J2 uncertainty of order 0.13e-6, about one tenth of the 1.4e-6 gap to Jacobson (2014); it is second-order.
- J6: READ eq. (C7), (C8) and Table 17: `J2(J6', Delta a) = 3509.281e-6 + 0.39909 (J6' - 0.5e-6) + 0.130e-6 Delta a(km)`, `J4(J6', Delta a) = -35.533e-6 + 1.07570 (J6' - 0.5e-6) - 0.180e-6 Delta a(km)`. COMPUTED consequence: the adopted J2 and J4 shift
  by +0.399e-6 and +1.0757e-6 per 1e-6 of J6. For J6 = 0 (as in Jacobson 2014) the model value would be J2 about 3509.09 (from 3509.291 - 0.5 x 0.39909 = 3509.091) and J4 about -36.06 (from -35.522 - 0.5 x 1.0757 = -36.060): J6 explains only about 0.2e-6 of the 1.4e-6 difference from Jacobson (2014).
  READ (section 10.3, p73): Neuenschwander and Helled (2022) require knowledge of J6 to 10 percent; "values J6 = (59.9-69.04) x 10^-8 can only be explained by deep winds" and "J6 = (46.12-53.76) x 10^-8 do not allow for deep winds"; the paper's adopted 0.5e-6 sits in that interior-model range. Jacobson (2023) set
  J6 = (0.58 +/- 0.12) x 10^-6 "based on a suite of Uranus interior models by Neuenschwander and Helled (2022)" (p59).

Comparison with the earlier determinations (READ, section 7.2.1, p59):
- French et al. (1988): "determined J2 and J4 from orbit fits to ring occultation observations from 1977-1986, with the results listed as Fit 1 ... converted from their assumed reference radius of R = 26200 km to R = 25559 km." (Fit 1: 3513.23 +/- 0.34, -30.32 +/- 4.73; the new adopted J2 is 3.94e-6 lower (COMPUTED), 10 times the paper's own uncertainty.)
- Jacobson (2014): "solved for J2 and J4 from a comprehensive analysis of all available Voyager navigation data and selected ring occultation measurements, and in a separate solution from the ring data alone" (Fits 2 and 3). "Jacobson (2014) assumed that J6 = 0, but the quoted uncertainties in J2 and J4 incorporated an a priori uncertainty sigma(J6) = 1.0 x 10^-6, as well as
  a priori uncertainties on other parameters of the fit including the absolute radius scale and the pole direction." The paper's explanation of the Jacobson (2014) pole offset (p40) is "due primarily to differences between the updated ura178 ephemeris and the earlier spacecraft ephemeris, as well as to the use by Jacobson (2014) of keplerian ring orbital elements that modeled the eta ring as inclined and eccentric,
  excluded the gamma and delta rings, and neglected the contributions of normal modes to the eta and epsilon ring shapes."
- Jacobson (2023): "estimated J2 and J4 as part of a global solution for the Uranus ephemeris and system properties, based on a comprehensive analysis of historical astrometry, navigation data, and ring occultation observations. He set J6 = (0.58 +/- 0.12) x 10^-6 ... and the satellite contribution to the precession of the rings was included for Ariel,
  Umbriel, Titania, Oberon, Miranda, and Puck only, ignoring the minor satellites." (p59). The authors' own COR fit (Fit 5, 6 satellites, J6 = 0.58) "are very similar to the Jacobson (2023) result (Fit 4), showing that these two independent solutions based largely on the same set of observations yield nearly identical results" (p60). Hence the Jacobson (2023) value
  is reproduced by this paper's code under Jacobson's assumptions; the 1.17e-6 difference is entirely this paper's added physics, a statement the authors make implicitly (COMPUTED in the table above).

How the authors present the status of the result: they do not claim Jacobson's J2 wrong; they present theirs as the better treatment and argue the pair (J2, J4, J6) is only meaningful jointly: "Although we cannot set useful independent limits on J6, we obtain strong joint constraints on combinations of J2, J4, and J6 that are consistent with our measurements" (abstract), with the Mahalanobis recipe of Appendix C
(r <= 1 for 1-sigma; for CCP = 30, 10, 1 percent, r = 1.55, 2.15, 3.03; eigenvalues lambda_1 = 3.84252e-13, lambda_2 = 2.64830e-15, theta = 48.5686 deg). COMPUTED test of Jacobson's two published pairs against that recipe (a short script applying eqs. C7 to C9 and C16 to C19 with the printed constants; the paper does not apply the recipe to them). Each pair is first compared with the J2, J4 the adopted solution would give at that pair's own J6 (eqs. C7, C8, with Delta a = 0):
Jacobson (2014), J6 = 0: model J2 = 3509.281 + 0.39909 x (0 - 0.5) = 3509.081, J4 = -35.533 + 1.0757 x (0 - 0.5) = -36.071; x = 3510.7 - 3509.081 = 1.619, y = -34.2 + 36.071 = 1.871 (x10^-6); with theta = 48.5686 deg, x' = 2.474, y' = 0.024; sqrt(lambda_1) = 0.6199, sqrt(lambda_2) = 0.0515 (x10^-6);
r = sqrt((2.474 / 0.6199)^2 + (0.024 / 0.0515)^2) = 4.0, CCP = exp(-4.0^2 / 2) = 3e-4 (the paper's value for 4 sigma is 0.0003). Jacobson (2023), J6 = 0.58: model J2 = 3509.313, J4 = -35.447; x = 1.152, y = 1.302; x' = 1.738, y' = -0.002; r = 2.8, CCP = 0.02.
Using Table 17's 3509.291 and -35.522 instead of the C-constants as the base changes r by under 0.03. INFERRED: by the authors' own criterion (r <= 1 within the 1-sigma ellipse; r = 1.55, 2.15, 3.03 for CCP of 30, 10, 1 percent) both Jacobson fields lie outside the adopted ellipse, Jacobson (2014) by far more than Jacobson (2023). The ellipse is wide (sigma(J2) = 0.412e-6), but it is narrow across the J2-J4 correlation line,
so the discrepancy is mostly a miss in the (J2, J4) combination, not in J2 alone.

## 4. The pole

READ (abstract, p1): "The Uranus pole direction at epoch TDB 1986 Jan 19 12:00 is alpha_P = 77.311327 +/- 0.000141 deg and delta_P = 15.172795 +/- 0.000618 deg. The slight pole precession predicted by Jacobson (2023) is not detectable in our orbit fits, and the absolute radius scale is not strongly correlated with the pole direction."
Frame (READ p8, p41, Fig. 24, "Uranus pole (J2000 frame)"): J2000 (ICRF), the direction of positive angular momentum, "180 deg from the IAU definition of the Uranus north pole" (so the IAU north-pole values near alpha_0 = 257.31, delta_0 = -15.18 are the antipode of this one; INFERRED: the IAU values are not in this paper).
The adopted pole is Fits 1 and 2 (Fit 2 has Voyager 2 trajectory uncertainties added by Monte Carlo displacement of the spacecraft position, which raised sigma(alpha) from 0.000128 to 0.000141 deg and sigma(delta) from 0.000472 to 0.000618 deg).

Precession (READ, section 5.1 Fit 3, p39-40): "Jacobson (2023) developed a trigonometric series representation of these contributions to the pole direction over time, derived from numerical integrations of the ura178 satellite ephemerides. Over the course of the 29-year interval of the ring occultation observations considered here, the corresponding predicted direction of
the pole varied periodically by +/- 0.0002 deg in alpha_P and +/- 0.00015 deg in delta_P, following the dominant short-period term in the series with a period of 17.79 years (6494 days)." Fit 3 (fixed trigonometric series, pole fitted at the reference epoch) gives (77.311210, 15.172762), within 1.17e-4 and 3.3e-5 deg of Fit 1, "virtually identical" ring semi-major axes, so "predicted pole precession has little effect on the derived ring system geometry,
as was also found by Jacobson (2014, 2023). ... Based on this insensitivity to the inclusion of pole precession, we adopt a fixed pole direction." So: precession is NOT modelled in the adopted pole; it is predicted to be tiny (a few 1e-4 deg) over the data span and 1600-2600 per Jacobson 2023 (the series itself is not printed in this paper; the Jacobson 2014 digest holds the 2014 series).

Comparison (READ Table 13 and section 5.1): Jacobson (2023) pole (Fit 5) 77.311200 +/- 0.000400, 15.172400 +/- 0.001700, "determined as part of the development of the ura178 series ... somewhat larger error ellipse (displaced slightly from our adopted Fit 1 pole)"; Jacobson (2014) (Fit 6) 77.310 +/- 0.003, 15.172 +/- 0.002, its offset from Fit 1 "due primarily to differences ..."
(quoted in section 3); French et al. (1991) (Fit 7, B1950 77.5969 +/- 0.0034, 15.1117 +/- 0.0033 precessed to J2000) 77.310877, 15.174564. Fit 4 (Earth-based occultations only) shows "the importance of the unique geometry of the Voyager occultations": rho = -0.88, sigma(delta) = 0.002013 deg, and a radius-scale offset <Delta a> = +0.278 km.
COMPUTED separations from the adopted pole (Table 13): Jacobson (2014): Delta alpha = 0.001327 deg x cos(15.1728 deg) = 0.001281 deg, Delta delta = 0.000795 deg, total sqrt(0.001281^2 + 0.000795^2) = 0.001507 deg = 5.43 arcsec = 2.63e-5 rad (11.5 km at Titania's distance 436,283 km, 15.4 km at Oberon's 583,447 km); Jacobson (2023): Delta alpha cos(delta) = 1.226e-4 deg, Delta delta = 3.95e-4 deg,
total 4.14e-4 deg = 1.49 arcsec = 7.2e-6 rad (3.2 km at Titania).

## 5. Constants audit

Project sources: `src/cyclerfinder/core/satellites.py` (`PRIMARIES["Uranus"]` line 70; `SATELLITES` lines 251 to 259) and `src/cyclerfinder/data/validation/v4_uranus.py` (`URANUS_J2` line 131, `URANUS_R_EQ_KM` line 148). (Note: `satellites.py` has uncommitted edits in the working tree
by another agent at the time of writing; the values below are those committed at HEAD and visible in the file at the time I read it.) Column (i) is Jacobson (2014) as transcribed in `docs/notes/2026-10-04-digest-jacobson-2014-uranian-satellites-gravity-field.md` (Tables 2 and 12); (ii) is Jacobson (2023) as quoted in THIS paper (Tables 2 and 17; the 2023 paper itself is not held);
(iii) is this paper's own adopted value. All differences are COMPUTED (project minus source) from the printed digits.

### 5.1 GM values (km^3 s^-2)

| Quantity | Project | (i) Jacobson 2014 | (ii) Jacobson 2023 (as quoted) | (iii) French 2024 own | Project minus (i) | Project minus (ii) |
| --- | --- | --- | --- | --- | --- | --- |
| Uranus GM | 5794556.4 ("system GM") | 5794556.4 +/- 4.3 (system) | GM_U 5793950.300 (planet; Table 17) | none fitted; adopts (ii) | 0 | +606.1 (but different objects: system vs planet) |
| Miranda | 4.3 | 4.3 +/- 0.2 | 4.20 +/- 0.20 | none | 0 | +0.10 (0.5 sigma) |
| Ariel | 83.5 | 83.5 +/- 1.4 | 82.30 +/- 1.20 | none | 0 | +1.20 (1.0 sigma) |
| Umbriel | 85.1 | 85.1 +/- 1.9 | 86.00 +/- 1.50 | none | 0 | -0.90 (0.6 sigma) |
| Titania | 226.9 | 226.9 +/- 4.1 | 230.60 +/- 3.40 | none | 0 | -3.70 (1.1 sigma) |
| Oberon | 205.3 | 205.3 +/- 5.8 | 207.60 +/- 5.00 | none | 0 | -2.30 (0.5 sigma) |
| Sum of five moons | 605.1 | 605.1 | 610.70 | | 0 | -5.60 |

Which GM is printed here: READ Table 17: "GM_U 5793950.300 km^3 s^-2". Eq. (2) (p8) defines GM in `sqrt(GM/a^3)` for the planet term and the satellites enter separately via the (m_j / M) terms with M the planet mass, and section 6 uses "M_Ur ... the planet mass" and `n ~ sqrt(GM_Ur / a^3)` (eq. 29). COMPUTED, and conclusive for the question asked: GM_U is the planet's GM, because the five M_sat/M_Ur ratios of Table 2 are
GM_sat / 5793950.3 to every printed digit (x 10^-6: Ariel 14.2045, Umbriel 14.8431, Titania 39.8001, Oberon 35.8305, Miranda 0.72489, Puck 0.022006 against the printed 14.204, 14.843, 39.800, 35.830, 0.7249, 0.0220), whereas dividing by the 2014 system value 5794556.4 gives Ariel 14.2030, Umbriel 14.8415, Titania 39.7960, Oberon 35.8267,
which miss the printed digits in the third or fourth place. The Puck ratio is not a test (0.1275 / 5793950.3 = 0.022006 against 0.0220 either way). Eq. (29) and section 6 also use "M_Ur" as the planet mass. The paper never prints a 2023 system GM. COMPUTED illustration, not a printed value: GM_U + the five moons = 5793950.300 + 610.70 = 5794561.0; adding Puck's 0.1275 and the minor moons' Table 3 sum (0.1648, assumed densities) gives about 5794561.3, 4.9 above the project's and Jacobson (2014)'s system value, 1.1 sigma of the latter's 4.3. With the 2014 values, GM_sys - sum of the five = 5794556.4 - 605.1 = 5793951.3 (planet only), 1.0 above Jacobson (2023)'s GM_U.

### 5.2 Semi-major axes (km)

| Satellite | Project | (i) Jacobson 2014 Table 2 | (ii) Jacobson 2023 (as quoted, Table 2) | Project minus (i) | Project minus (ii) | (ii) minus (i) |
| --- | --- | --- | --- | --- | --- | --- |
| Miranda | 129846 | 129858 | 129828 | -12 | +18 | -30 |
| Ariel | 190929 | 190930 | 190928 | -1 | +1 | -2 |
| Umbriel | 265986 | 265982 | 265981 | +4 | +5 | -1 |
| Titania | 436298 | 436282 | 436283 | +16 | +15 | +1 |
| Oberon | 583511 | 583449 | 583447 | +62 | +64 | -2 |

READ: both published columns are integer-kilometre mean semi-major axes of a precessing-ellipse fit ("mean equatorial orbital elements" in the 2014 paper). COMPUTED: the project's values match none of the sources: Ariel is within 1 km of both (a coincidence of integer rounding at that level); the Titania and Oberon values are 15 and 64 km above both. INFERRED, as in the earlier digest: the project's values come from the
JPL Satellite Mean Elements page, a different mean-element definition (epoch, averaging) of the URA111 ephemeris, so they are not in error by a mistake but cannot be audited against these papers. READ (this paper): neither Table 2 nor any other table gives an epoch or definition. The project's semi-major axes are therefore STILL of unknown source relative to this paper. Mean motions: the Jacobson (2014) Table 2 rates (periods 8.705868 d, 13.463237 d) remain the right
ones for URA111 (earlier digest, section 8); this paper prints no satellite mean motions for the five major moons.

### 5.3 Zonal harmonics and reference radius

| Quantity | Project | (i) Jacobson 2014 | (ii) Jacobson 2023 (as quoted) | (iii) French 2024 adopted | Project minus (i) | Project minus (ii) | Project minus (iii) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| J2 x10^-6 | 3509.291 | 3510.7 +/- 0.7 | 3510.465 +/- 0.058 | 3509.291 +/- 0.412 | -1.409 (2.0 sigma of 0.7) | -1.174 (20 sigma of 0.058; 2.8 sigma of French 0.412) | 0 |
| J4 x10^-6 | not carried | -34.2 +/- 1.3 | -34.145 +/- 0.082 | -35.522 +/- 0.466 | - | - | - |
| J6 x10^-6 | not carried | 0 (not estimated; sigma 1.0) | 0.58 +/- 0.12 | 0.5 (fixed) | - | - | - |
| Reference radius (km) | 25559.0 | 25559 | 25559 (Table 17 "Reference radius for J_n") | 25559 | 0 | 0 | 0 |
| Pole alpha, delta | not examined | 77.310 +/- 0.002, 15.172 +/- 0.002 | 77.3112 +/- 0.0004, 15.1724 +/- 0.0017 | 77.311327 +/- 0.000141, 15.172795 +/- 0.000618 | | | |
| Retired J2 3.34343e-3 with 25559 km | removed 2026-10-04 | not printed there | not printed | not printed | | | |

The retired J2 is traceable to this paper. COMPUTED: Table 17 Fit 1 gives French et al. (1988) as J2 = 3513.23e-6 at R = 25559 km, "converted from their assumed reference radius of R = 26200 km". J2 R^2 is the invariant, so the same field at R = 26200 km is 3513.23 x (25559 / 26200)^2 = 3513.23 x 0.951667 = 3343.43 (x10^-6) = 3.34343e-3, the old constant to all six digits.
(The earlier digest backed out about 26,190 km from Jacobson's J2 and the old value; 26,200 km is the printed radius, and the match is to 3513.23, not to Jacobson's 3510.7, so the origin is French et al. 1988, not Jacobson 2014. The code's attribution to "Jacobson 2014 Table 4" was wrong on both counts.) The corresponding 1988 J4 at R = 26200 would be -30.32 x (25559 / 26200)^4 = -27.46e-6 (COMPUTED); it was never carried.
The reference radius 25559 km is a normalisation of the harmonics in all four sources; this paper does not give a separate physical equatorial radius of Uranus. The project's `URANUS_R_EQ_KM` docstring cites "Jacobson 2014 ... also IAU 2015 nominal value (Mamajek et al.)"; this paper supports only the first as a J2 reference radius.

### 5.4 Which project values match which source

- Match Jacobson (2014) exactly: the system GM (5794556.4) and all five moon GMs.
- Match this paper (French 2024) exactly: `URANUS_J2` (3509.291e-6) and the reference radius 25559.
- Match Jacobson (2023) (as quoted here) on nothing: GMs differ by 0.1 to 3.7 km^3/s^2 (all within 1.1 sigma); semi-major axes by 1 to 64 km; J2 by 1.174e-6.
- Match no source: all five semi-major axes (a different mean-element set), per the earlier digest and confirmed against the second source.
- Still of unknown source after this paper: the five moon mean radii (Miranda 235.8, Ariel 578.9, Umbriel 584.7, Titania 788.9, Oberon 761.4 km), the 50 km and 100 km flyby-altitude floors (the registry attributes the 50 km values to Heaton-Longuski 2003 Table 4 and the 100 km for Miranda to convention), the semi-major axes (the JPL Satellite Mean Elements, per the registry comment, not auditable here), and the DE440 attribution of the system GM (this paper does not
  mention DE440 GM values; the number equals Jacobson 2014's system GM digit for digit). This paper contains no major-moon radius at all (READ: Table 3 gives dimensions only for the ten minor satellites and Puck; the major moons' radii are not tabulated), so it cannot source or contradict the registry's radii or floors.

## 6. Which J2, J4, J6 and pole should a trajectory propagator use (item 6)

### 6.1 What the paper says about consistency between a gravity field and the ephemeris fitted with it

READ, bearing on consistency:
1. Section 3.3 (p11): the ura178 satellite ephemerides are an update of ura111 that used the ring occultation data (not the gamma and lambda rings); section 7.2.1 (p59) calls Jacobson (2023) "a global solution for the Uranus ephemeris and system properties". Jacobson (2023)'s J2, J4, J6 (Table 17 Fit 4) are thus the field of the ura178 ephemeris; Jacobson (2014)'s (Fit 2) are the field of ura111.
2. Section 3.1.1 (p7): "With the release of the Gaia EDR3 and DR3 catalogs ... the situation has reversed, with small but measurable systematic errors in the Uranus ephemeris in ura111.bsp and de440.bsp exceeding the star position uncertainties. To reduce these systematic errors, the ura178 series of ephemerides ... was developed". Section 4.8 (p38): "the substantial systematic drift in time of the ura111 Uranus ephemeris that
   amounted to several hundred km in the sky plane by the time of the final ring occultation in 2006 has been effectively eliminated in the ura178 ephemeris" (this is the planet's heliocentric ephemeris, relevant to Earth-pointing geometry, not to the moon-relative geometry of a Uranus-centred tour).
3. Section 5.1 (p40): the Jacobson (2014) pole offset from the modern pole is "due primarily to differences between the updated ura178 ephemeris and the earlier spacecraft ephemeris", i.e. the fitted pole and the ephemeris travel together, and the pole offset is not random.
4. Section 7.3.2: the adopted French J2 and J4 differ from Jacobson's because of ring-physics systematics (COO versus COR, minor-satellite masses), not because of new data; this paper's J2 is a ring-precession field that has not been fed back into any satellite-orbit fit, so it is not self-consistent with any JPL ephemeris in the sense that the ephemeris's satellite-orbit fit used a different J2.
5. The satellite-precession terms in the paper's own J2 fit use the Jacobson (2023) masses (p9, footnote 3), so French's J2 is conditional on the 2023 moon GMs (Fit 17), as noted in section 3.

### 6.2 Recommendation (INFERRED)

(a) With the URA111 moon ephemeris: use the Jacobson (2014) set as one self-consistent bundle: system GM 5794556.4 (or, if the moons are also modelled as separate bodies, the planet GM of the same solution; see 6.4), moon GMs from Table 12 of that paper (the project's present values), J2 = 3510.7e-6, J4 = -34.2e-6, J6 = 0, R = 25559 km, pole alpha = 77.310, delta = 15.172 deg. Reason: those values were fitted together with URA111 (Jacobson 2014 section 2
and Table 12; this paper's Fit 2 reproduces Table 12 exactly). The French (2024) J2 is a newer and, by its authors' account, better ring-physics estimate, but its gain does not reach the ephemeris; mixing it with URA111 is not wrong at any level the data can see (6.3) but it is a mixed set. If the project prefers the French field for scientific reasons, take it together with the French pole (77.311327, 15.172795) and note the mix.

(b) With a newer JPL kernel (ura178 or later): the self-consistent field is Jacobson (2023): J2 = 3510.465e-6, J4 = -34.145e-6, J6 = 0.58e-6, R = 25559 km, GM_U = 5793950.300 km^3/s^2 (a PLANET GM, to be used together with the moon GMs of Table 2 of this paper), pole alpha = 77.3112, delta = 15.1724 deg (+/- 0.0004, 0.0017). These are second-hand, from this paper's Table 17 (Fit 4), Table 13 (Fit 5) and Table 2; the primary source (Jacobson 2023, DPS meeting abstract)
is not held. The newer JPL solution with the same family of numbers is Jacobson and Park (2025, URA182), whose PDF appears in the private corpus and is not indexed (headline finding 9); it should be read and, if its J2, J4, J6, pole and GMs match, those should replace the second-hand 2023 numbers. None of the three solutions' J2 values is "the" correct one: they differ in modelling choices, and the paper's Fit 2 to Fit 15 range is 1.4e-6.

(c) Pole: for any propagator the pole matters at the level computed below; for URA111 use the 2014 pole, for a newer kernel use the 2023 or 2025 pole. The precession terms (about 1e-4 deg over decades) are negligible for a propagation of a few hundred days.

### 6.3 Size of the J2 choice (COMPUTED; arithmetic shown)

Method. In the equatorial plane at radius r, the J2, J4 and J6 terms change the radial acceleration magnitude by `g = (GM / r^2) [1 + (3/2) J2 x^2 - (15/8) J4 x^4 + (35/16) J6 x^6]` with x = R / r, from the Legendre values P2(0) = -1/2, P4(0) = 3/8, P6(0) = -5/16. Take R = 25559 km, GM = GM_U = 5793950.3 km^3/s^2 (the planet GM; using the system GM changes everything below by 1e-4 relative).
Titania: r = 436,283 km (Table 2 of this paper), x = 25559 / 436283 = 0.058582, x^2 = 3.43203e-3, point-mass g = GM / r^2 = 5793950.3 / (436283)^2 = 3.04395e-5 km/s^2. Oberon: r = 583,447 km, x^2 = 1.91904e-3, g = 1.70205e-5 km/s^2. T = 123 d = 123 x 86400 = 1.06272e7 s; T^2 = 1.12937e14 s^2.

Candidate values: French (3509.291e-6, -35.522e-6, 0.5e-6), Jacobson 2023 (3510.465e-6, -34.145e-6, 0.58e-6), Jacobson 2014 (3510.7e-6, -34.2e-6, 0).

J2 term alone, Delta g = (GM / r^2) x (3/2) x Delta J2 x x^2:

| Difference | Delta J2 | Titania Delta g (km/s^2) | Titania (m/s^2) | Oberon Delta g (km/s^2) | Oberon (m/s^2) |
| --- | --- | --- | --- | --- | --- |
| Jacobson 2014 minus French | +1.409e-6 | 2.208e-13 | 2.21e-10 | 6.90e-14 | 6.90e-11 |
| Jacobson 2023 minus French | +1.174e-6 | 1.840e-13 | 1.84e-10 | 5.75e-14 | 5.75e-11 |
| Jacobson 2014 minus Jacobson 2023 | +0.235e-6 | 3.68e-14 | 3.68e-11 | 1.15e-14 | 1.15e-11 |
| French 1-sigma (0.412e-6) | 0.412e-6 | 6.46e-14 | 6.46e-11 | 2.02e-14 | 2.02e-11 |
| Retired project value (3343.43e-6) minus French | -165.86e-6 | -2.60e-11 | -2.60e-8 | -8.13e-12 | -8.13e-9 |

Example line: Titania, Jacobson 2014 minus French: 3.04395e-5 x 1.5 x 1.409e-6 x 3.43203e-3 = 3.04395e-5 x 7.2536e-9 = 2.208e-13 km/s^2. Including the J4 and J6 terms with their own values changes the full-field difference to 2.199e-13 (Titania, J2014 minus French) and 6.89e-14 (Oberon), i.e. under 0.5 percent from the J2-only number, because the J4 term is down by a factor of order x^2 and has the opposite sign (J4 is more negative in French's solution: the strong rho = 0.986 correlation pairs a lower J2 with a more negative J4). So at moon distances the differences between candidate (J2, J4, J6) sets are essentially their J2 difference.

Displacement over 123 days. Two bounds:
- Naive upper bound, the whole Delta g acting as a constant unopposed acceleration: d = (1/2) Delta g T^2. Titania: 0.5 x 2.208e-13 x 1.12937e14 = 12.5 km (Jacobson 2014 minus French), 10.4 km (Jacobson 2023 minus French), 2.1 km (Jacobson 2014 minus Jacobson 2023). Oberon: 0.5 x 6.90e-14 x 1.12937e14 = 3.9 km, 3.2 km and 0.65 km. For scale, the retired value: Titania -1468 km, Oberon -459 km (0.5 x 2.60e-11 x 1.12937e14 and 0.5 x 8.13e-12 x 1.12937e14), which is what the #894 correction removed.
- Mean-motion (phase) estimate, appropriate to a body in a near-circular orbit at fixed radius where a small change in the radial force is absorbed as a change in angular rate: Delta n / n = (1/2) Delta g / g, drift = r (Delta n) T = r x (1/2)(Delta g / g) x n x T. Titania: Delta g / g = 2.208e-13 / 3.04395e-5 = 7.25e-9; n T = (2 pi / 8.7063 d) x 123 d = 88.77 rad; drift = 436,283 km x 0.5 x 7.25e-9 x 88.77 = 0.140 km (Jacobson 2014 minus French), 0.116 km (Jacobson 2023 minus French). Oberon: Delta g / g = 4.05e-9, n T = (2 pi / 13.464) x 123 = 57.40 rad, drift = 583,447 x 0.5 x 4.05e-9 x 57.40 = 0.068 km (0.056 km, Jacobson 2023 minus French).
The truth for a spacecraft leg lies between the two (a re-targeted leg with fixed endpoints absorbs most of the effect). INFERRED: at the 1 km level the J2 choice among the three candidates does not matter for a self-consistent re-targeted tour (phase bound 0.14 km) and could matter, up to 12 km, only for an open-loop propagation of a fixed initial state for the full 123 days.

Pole tilt. If the J2 axis is wrong by the Jacobson (2014) minus French angle of 1.507e-3 deg = 2.63e-5 rad (section 4), the equatorial plane at Titania's distance is displaced by 436,283 x 2.63e-5 = 11.5 km, and the out-of-plane J2 acceleration of a particle in the true equator is about 3 J2 x^2 g (z / r) with z / r = 2.63e-5: 3 x 3.5105e-3 x 3.43203e-3 x 3.04395e-5 x 2.63e-5 = 2.9e-14 km/s^2, a naive displacement 0.5 x 2.9e-14 x 1.12937e14 = 1.6 km over 123 days. The Jacobson (2023) minus French angle of 7.2e-6 rad gives a 3.2 km plane displacement and 0.4 km naive displacement (scaling). So the pole choice is also about the 1 km level at worst.

Against this, the ephemeris errors dominate: READ (Jacobson 2014, Table 7, per the earlier digest) the URA111 orbit uncertainties for 2000 to 2050 are 35 km radial, 160 km tangential and 50 km normal for Titania, growing 3 km/yr along track; the choice between these J2 values is two orders of magnitude below that.

### 6.4 A larger effect the arithmetic exposes (INFERRED, not run)

`src/cyclerfinder/data/validation/v4_uranus_strict.py:769` sets `mu_primary = PRIMARIES[primary]`, the system GM 5794556.4, and the same module adds all five major moons as third-body perturbers. INFERRED: the moons' mass is then counted twice, once inside the system GM and once as separate perturbers; the excess central GM is the sum of the five moon GMs, 605.1 km^3/s^2 (Jacobson 2014) = 605.1 / 5794556.4 = 1.044e-4
of the central GM. COMPUTED: 1.044e-4 x 3.044e-5 km/s^2 = 3.18e-9 km/s^2 at Titania's distance and a naive displacement 0.5 x 3.18e-9 x 1.12937e14 = 1.8e5 km over 123 days; as a phase effect (n changes by half of 1.044e-4) the drift is 436,283 km x 0.5 x 1.044e-4 x 88.77 = 2.0e3 km. Both are four orders of magnitude above the largest J2 choice (3.18e-9 / 2.21e-13 = 1.4e4) and above the 160 km tangential uncertainty of URA111. This paper supplies the cure
as a number: the planet-only GM (Jacobson 2023 GM_U = 5793950.300 for the 2023 moon set, or 5794556.4 - 605.1 = 5793951.3 for the 2014 set, COMPUTED) used with explicit moon perturbers. The Jacobson 2014 paper itself says "We use the system GM rather than the Uranus GM as our gravity parameter" (earlier digest, section 1). The coordinator should check whether the lane really double counts before relying on this paragraph; I read only
the cited lines.

## 7. Other things a moon-tour trajectory designer would want from this paper

- Ring radii and extents (hazards; READ Table 5, COMPUTED extents with R = 25559 km). The nine main rings plus lambda lie between 41,794 km and 51,604 km from the planet's centre. Innermost: ring 6 IER periapse a(1 - e) = 41835.920 - 42.068 = 41,793.9 km (1.635 R_U). Outermost: epsilon OER apoapse 51178.588 + 425.242 = 51,603.8 km (2.019 R_U). The epsilon ring (width "from 19.9 km at periapse to 97.3 km at apoapse", p49; Table 14 mean width 58.574 km) has COR periapse 51149.279 - 405.894 = 50,743.4 km and apoapse 51,555.2 km, IER periapse 50,733.5 km. Lambda: a = 50026.557 +/- 1.314 km, "its semimajor axis should probably be regarded as uncertain at the level of a few km" (p37), and it is a dusty ring. The paper
  covers narrow, dense rings only; the dusty rings and broad dust sheets (the Voyager "mu" and "nu" outer rings and the rings between) are mentioned in the introduction (p2) but no radii are given, so they are not a source here.
- Ring-plane inclination constraint (READ Table 5): the ring planes coincide with the planet's equatorial plane to a few tens of km. a sin i (the out-of-plane displacement amplitude at each ring's own radius) is 44.6 (ring 6), 40.9 (5), 23.4 (4), 12.0 (alpha), 4.0 (beta) km; eta, gamma, delta and epsilon held at zero with 2-sigma limits of 0.290, 0.286, 0.284 and 0.158 km (Table 9). Out-of-plane distance of any ring particle from the equatorial plane is therefore under about 45 km. The plane's orientation is the pole of section 4. The paper
  contains no ring-plane-crossing times for a spacecraft; the 2007 Earth ring-plane crossing is mentioned only in the references (de Pater et al. 2013). READ the geometry statements of section 2.1: from Earth the 1977-2006 view was nearly pole-on, so the ring "inclinations and the planet's pole direction" were initially poorly determined, which the Voyager chords fixed.
- Inner moons (READ Table 3, section 9): the orbital radii and masses of the ten Voyager-era inner satellites and Puck, section 2 above. The three whose masses are measured from the rings are Cressida (GM = 12.27 +/- 1.41 x 10^-3 km^3/s^2 at a = 61,767 km, mean motion 776.584048 deg/day, sqrt(AB) 41 km), Cordelia (4.06 +/- 0.38 x 10^-3 at 49,752 km, 1074.522858 deg/day, 21 km) and Ophelia (2.38 +/- 0.22 x 10^-3 at 53,764 km, 956.418404 deg/day, 23 km); the radii are the Table 3 sqrt(AB) values. Cordelia and Ophelia are the epsilon ring's shepherds ("the epsilon ring is securely established as being confined by shepherd satellites", p72) and sit 1,850 km inside and 2,160 km outside the epsilon ring's outer extent
  (49,752 and 53,764 against 51,604 km at apoapse OER: COMPUTED differences 51,603.8 - 49,752 = 1,852 km and 53,764 - 51,603.8 = 2,160 km, approximate since their a is a mean value). Puck has a = 86,004 km, GM = 0.1275 +/- 0.0425 km^3/s^2 (Table 2) or 133.72e-3 (Table 3, assumed density) and a measured mass ratio of 0.0220e-6 of the planet. A capture trajectory that reaches periapsis inside 51,600 km crosses the ring system; one that stays outside 53,800 km clears Ophelia's orbit, and the tour moons (Miranda to Oberon, 129,828 to 583,447 km) are all outside.
- Satellite orbit accuracy and normal-mode pattern speeds: Table 20 gives the three inner moons' mean motions with 1e-4 to 2e-3 deg/day accuracy, "ura115" minor-satellite kernel (Table 4); nothing in this paper changes the five major moons' orbits (ura178, Jacobson 2023).
- Voyager-era ring data tell nothing about the moons' absolute orbits beyond what Jacobson already used; this paper does not fit moon orbits.
- A sentence for the V4 lane: the J2 and J4 here are ring-precession fields over 41,000 to 52,000 km; extrapolation to 130,000 to 583,000 km (the moon tour) is by power law in R/r and is exact only for the harmonics, not for any unmodelled structure; J6 and higher terms are negligible there (COMPUTED: the J6 term at Titania relative to the point mass is (35/16) x 0.5e-6 x x^6 = 2.1875 x 0.5e-6 x 4.04e-8 = 4.4e-14).

## 8. References the project might want

Corpus status from `docs/notes/CORPUS_INDEX.md` and the private paper corpus directory listing at the time of writing (grep of the index for each name).

| Reference | Full citation | Held? |
| --- | --- | --- |
| Jacobson (2023) | R. A. Jacobson, "Update of the orbits of the regular Uranian satellites, the gravity field of the Uranian system, and the Uranus pole direction", DPS Meeting 2023, San Antonio, TX (a conference abstract; as cited on p92) | No |
| Jacobson and Park (2025) | R. A. Jacobson and R. S. Park, "Orbits of Uranus satellites, rings, gravity field, poles", AJ 169:65, DOI 10.3847/1538-3881/ad99d1 (URA182) (title from the corpus file name, not read) | In the private paper corpus (file timestamp 2026-10-04 21:18); NOT in `CORPUS_INDEX.md` |
| Jacobson (2014) | R. A. Jacobson, "The orbits of the Uranian satellites and rings, the gravity field of the Uranian system, and the orientation of the pole of Uranus", AJ 148:76, DOI 10.1088/0004-6256/148/5/76 | Yes, digested 2026-10-04 |
| French et al. (2023a) | R. G. French and 27 colleagues, "Uranus ring occultation observations: 1977-2006", Icarus 395:115474, DOI 10.1016/j.icarus.2023.115474 (Paper 1; the data of this paper) | No |
| French et al. (1988) | R. G. French and 10 colleagues, "Uranian ring orbits from Earth-based and Voyager occultation observations", Icarus 73:349-378 (J2 = 3513.23e-6 at 25559 km as converted here; 3343.43e-6 at 26200 km; the source of the retired project J2) | No |
| French et al. (1991) | R. G. French, P. D. Nicholson, C. C. Porco, E. A. Marouf, "Dynamics and structure of the Uranian rings", in Uranus, pp. 327-409 (the 1991 pole, Fit 7) | No |
| Jacobson (1998) | R. A. Jacobson, "The orbits of the inner Uranian satellites from Hubble Space Telescope and Voyager 2 observations", AJ 115:1195-1199 (Cordelia, Ophelia, Cressida mean motions, Table 20) | No |
| Chancia et al. (2017) | R. O. Chancia, M. M. Hedman, R. G. French, "Weighing Uranus' moon Cressida with the eta ring", AJ 154:153, DOI 10.3847/1538-3881/aa880e | No |
| Neuenschwander and Helled (2022) | B. A. Neuenschwander, R. Helled, "Empirical structure models of Uranus and Neptune", MNRAS 512:3124-3136, DOI 10.1093/mnras/stac628 (interior models behind the J6 = 0.58e-6 value) | No |
| Park et al. (2021) | R. S. Park, W. M. Folkner, J. G. Williams, D. H. Boggs, "The JPL planetary and lunar ephemerides DE440 and DE441", AJ 161, DOI 10.3847/1538-3881/abd414 | No (the index has no DE440 paper; grep for "de440" finds only unrelated rows) |
| Karkoschka (2001a, b) | E. Karkoschka, "Comprehensive photometry of the rings and 16 satellites of Uranus with the Hubble Space Telescope", Icarus 151:51-68; "Voyager's eleventh discovery of a satellite of Uranus and photometry and the first size measurements of nine satellites", Icarus 151:69-77 (the minor-satellite dimensions of Table 3) | No |
| Showalter and Lissauer (2006) | M. R. Showalter, J. J. Lissauer, "The second ring-moon system of Uranus: Discovery and dynamics", Science 311:973-977 (Cressida mean motion in Table 20) | No |
| Laskar and Jacobson (1987) | GUST86 satellite theory (cited in the Jacobson 2014 digest, not in this paper) | No |

Kernel names this paper uses (Table 4): ura178.bsp, pleph.ura178.bsp (text: peph.ura178.bsp), vgr2.ura178.bsp, ura115.bsp, naif0012.tls, and the older series it was updated from: ura111.bsp, ura116.bsp, vgr2.ura111.bsp, de440.bsp. The project's kernel is ura111.bsp (from `v4_uranus_strict.py`, per its docstring); `naif0012.tls` is already vendored in the repository
(see the NAIF LSK memory note). The paper states the kernels are "available from ftp://ssd.jpl.nasa.gov/pub/eph" (note: the ura178 kernels are not yet in the project's GMAT-bundled kernel set per `v4_uranus_strict.py`; whether they are public was not checked).

## 9. Proposed corpus-index row

| french-hedman-nicholson-longaretti-mcghee-french-2024-uranus-system-occultation-observations-1977-2006-rings-pole-gravity-field-icarus-411-115957-doi-10.1016-j.icarus.2024.115957-hal-accepted-manuscript.pdf | 2026-10-04-digest-french-2024-uranus-system-occultations-gravity-field.md | Uranus gravity field J2 = 3509.291e-6, J4 = -35.522e-6 (J6 fixed 0.5e-6, R = 25559 km, rho 0.9861), pole alpha 77.311327 / delta 15.172795 (epoch TDB 1986 Jan 19), nine-ring orbits, masses of Cressida, Cordelia, Ophelia; quotes Jacobson (2023) GM_U = 5793950.300 (planet), moon GMs and semi-major axes; Table 17 (every fit), Tables 2, 3, 5, 13, 15 transcribed; origin of the retired 3.34343e-3 J2 identified as French et al. 1988 at R = 26200 km | digested | text-layer |

## 10. Unresolved points

- The 4-digit evidence for planet versus system GM_U is suggestive from the mass ratios but not conclusive on its own; section 6 of the paper uses `M_Ur` "the planet mass" with the same symbol. A reading of Jacobson (2023) or of Jacobson and Park (2025) would settle it.
- Jacobson (2023)'s Table 2 semi-major axes, and whether they are mean elements with Jacobson (2014)'s definition, are not stated here.
- Internal inconsistencies of this manuscript (all READ): Fit 6 alpha error 0.003 against Jacobson (2014)'s 0.002; Fits 6 and 7 correlation printed [1.0], text says 0; eq. (C4) and (C7) constants 3509.281 and -35.533 against Table 17's 3509.291 and -35.522; the planetary kernel name in section 3.3 versus Table 4; the Table 19 epoch longitude source, ura111 in the note and ura115 in the text. None changes a number the project uses.
- Appendix C's Mahalanobis values (r about 4.0 for Jacobson 2014 and 2.8 for Jacobson 2023) are our computation from the printed constants (script, not by hand); the paper does not apply its recipe to Jacobson's values.
