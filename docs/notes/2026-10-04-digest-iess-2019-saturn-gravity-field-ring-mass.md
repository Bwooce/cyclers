# Digest: Iess et al. (2019), "Measurement and implications of Saturn's gravity field and ring mass"

Science 364 (2019), DOI 10.1126/science.aat2965 (authors L. Iess, B. Militzer, Y. Kaspi, P. Nicholson, D. Durante, P. Racioppa,
A. Anabtawi, E. Galanti, W. Hubbard, M. J. Mariani, P. Tortora, S. Wahl, M. Zannoni). The held file is the "first release" of
17 January 2019 (received 12 February 2018, accepted 19 December 2018), 15 pages, whose footer reads "Page numbers not final at time
of first release"; the page numbers below are the first-release page numbers 1 to 15, so they will not match the printed issue pages.
Filed in the private paper corpus as
iess-et-al-2019-measurement-implications-saturn-gravity-field-ring-mass-science-364-doi-10.1126-science.aat2965.pdf

The supplementary materials (Materials and Methods, Figs. S1 to S6, Tables S1 and S2) are NOT in the held file (page 10 only lists
them). Everything below about Tables S1 and S2 is what the main text says about them.

Digested 2026-10-04 from all 15 pages. Statements are marked READ (seen on the page, page and section or table given), COMPUTED
(our arithmetic, shown) or INFERRED (our reading or a comparison with project code). Digits were read from page images.

## 0. Headline findings for the project

1. READ (Table 1, p13): J2 = 16290.573 +/- 0.028 (x10^6), un-normalised, reference radius 60,330 km. The project's `SATURN_J2 = 16290.573e-6` and
   `SATURN_J2_REF_RADIUS_KM = 60330.0` match the paper digit for digit. The #894 correction (do not pair this J2 with 60,268 km) is correct.
2. READ (Table 1 caption, p13): "The J2 value includes a constant tidal term owing to the average tidal perturbation from the satellites."
3. READ (Table 1 caption, p13): "The associated uncertainties are recommended values to be used for analysis and interpretation. For the zonal harmonics
   they correspond to 3 times the formal uncertainties." So 0.028 is a 3-sigma-formal figure (formal 1 sigma about 0.0093, COMPUTED 0.028/3).
4. The paper does not print Saturn's GM, nor any moon GM except Mimas (GM 2.5026 km^3/s^2, p2). `PRIMARIES["Saturn"]` is not checkable here.
5. COMPUTED: at Titan's orbit the J4 equatorial acceleration is 1.75e-4 of the J2 acceleration and 1.0e-8 of the central acceleration; at Iapetus 2.06e-5
   of J2 (section 6). J4 and higher are negligible for the project's Saturn lanes at the precision they work to.
6. The paper's total ring GM is 1.02 +/- 0.41 km^3/s^2 on p3 and 1.02 +/- 0.43 on p7 (an internal inconsistency in the first-release text, READ);
   the abstract and Table 1 give 0.41 +/- 0.13 Mimas masses.

## 1. What the paper does

READ (abstract, p1; section "Cassini gravity measurements", p1-2): during the Grand Finale, Cassini dove between Saturn and its innermost ring "at altitudes
2600-3900 km above the cloud tops". A radio link with Earth was monitored for six crossings to determine Saturn's gravitational field and the ring mass.
READ (abstract): "We find that Saturn's gravity deviates from theoretical expectations and requires differential rotation of the atmosphere extending to a depth
of at least 9000 km. The total mass of the rings is (1.54 +/- 0.49)x10^19 kg (0.41 +/- 0.13 times that of the moon Mimas), indicating that the rings may have
formed 10^7-10^8 years ago."

Data (p1-2): coherent microwave link, range-rate (Doppler) at X-band (7.2 GHz uplink, 8.4 GHz downlink) with auxiliary Ka-band (32.5 GHz). In April 2017
Cassini entered a series of inclined, highly eccentric orbits; of the 22 Grand Finale orbits (Rev 271 to 293) six were selected, five (Revs 273, 274, 278, 280,
284) provided useful data (Rev 275 lost to a station configuration error). Doppler counts at 30 s, 24 to 36 hours of data about each closest approach, DSN
complexes at Goldstone, Madrid, Canberra plus two ESA ESTRACK antennas (Malargue, New Norcia). Two-way X-band is 93 percent of the data set. Elevation
below 15 degrees discarded. Doppler noise rms 0.020 to 0.088 mm/s at 30 s. The Doppler signatures of the weakest measurable harmonics (J3, J10) and of the ring
are 40 to 200 times larger than the average Doppler noise.

Dynamical model (p2): the JPL MONTE code (previously used for the Titan and Enceladus gravity fields). Estimated: Saturn's GM, zonal harmonics J2 to J20,
tesseral degree-2 terms (C21, C22, S21, S22, "to account for possible non-principal axis rotation"), the masses of the A, B and C rings. The zonal truncation was set
to twice the degree of the highest harmonic whose central value is above its uncertainty (degree 10). Rings assumed coplanar with Saturn's equator, constant surface
density each; Saturn's spin axis position and precession rate were taken from ring occultations (French et al. 2017) as a prior, "about ten times more accurate than that
obtained from our orbital fitting". READ quote (p2): accelerations accounted for: "the point-mass gravitational accelerations from Saturn and its satellites (including
the ring moons), computed from the JPL planetary and satellite ephemerides DE430, SAT389 and SAT393 (18), the acceleration from the Sun, the planets and satellites of
the Solar System, and Saturn's tidal response to its satellites (12)". Multi-arc, weighted least-squares filter; global parameters (GM, J2 to J20, C/S terms, ring masses) common to all arcs;
local parameters per arc (position and velocity at the start, a priori sigma 100 km and 1 m/s), plus random stochastic accelerations (a priori 4e-7 m/s^2 over 10 minute
intervals within +/-1 hour of pericentre). Consider parameters: Love number k22, pole direction, RTG thermal and solar radiation pressure accelerations. The "A+B+C" ring mass sum was constrained to the satellite-ephemeris value.
Mimas: READ (p2): "Mimas has a GM of 2.5026 km^3 s^-2".

Rotation periods considered (p3): 10h32m45s, 10h39m22s (System III), 10h45m45s, 10h47m06s.

Result statements (p3-4): the Grand Finale data are "consistent with" the earlier estimates of J2, J4, J6 from moon orbit perturbations and Cassini itself (reference 18); the new data add J8 and J10 and the odd harmonics J3 and J5, "the only odd harmonics whose values are larger than the associated uncertainties". The baseline solution (Table 1) uses empirical random accelerations; a tesseral-field alternative and a normal-mode alternative agree with it (fig. S4, table S2).

## 2. Tables, transcribed

### Table 1 (p13). READ caption, verbatim: "Measured gravity harmonic coefficients of Saturn (un-normalized; reference radius 60330 km) and total ring mass (in units of Mimas' mass). The J2 value includes a constant tidal term owing to the average tidal perturbation from the satellites. The associated uncertainties are recommended values to be used for analysis and interpretation. For the zonal harmonics they correspond to 3 times the formal uncertainties. The solution for the total ring mass (A+B+C) is stable independently of the adopted dynamical model (table S2) and the uncertainty reported is the 1 sigma formal uncertainty. See table S2 for our total ring mass estimates for several models of the unknown accelerations."

| Quantity | Value | Uncertainty |
| --- | --- | --- |
| J2 (x10^6) | 16290.573 | 0.028 |
| J3 (x10^6) | 0.059 | 0.023 |
| J4 (x10^6) | -935.314 | 0.037 |
| J5 (x10^6) | -0.224 | 0.054 |
| J6 (x10^6) | 86.340 | 0.087 |
| J7 (x10^6) | 0.108 | 0.122 |
| J8 (x10^6) | -14.624 | 0.205 |
| J9 (x10^6) | 0.369 | 0.260 |
| J10 (x10^6) | 4.672 | 0.420 |
| J11 (x10^6) | -0.317 | 0.458 |
| J12 (x10^6) | -0.997 | 0.672 |
| Ring mass (M_M, Mimas masses) | 0.41 | 0.13 |

(The first-release text prints "x10^6" for each harmonic; the J-column unit marker on J12 is read as x10^6 as well. Saturn's GM is not in Table 1.)

### Table 2 (p13). READ caption, verbatim: "Comparison of observed and calculated gravitational harmonics (un-normalized; reference radius 60330 km). Where two values are given they denote the minimum and maximum values from the suite of models. The physical models in column 3 match the observed J2 and J4 in Table 1, over a parameter space considering ranges of S_met, Y_mol, Z_mol, r_c and rotation periods from 10h32m44s to 10h47m06s. For the same span of rotation periods, column 4 reports a wider range from models that match only J2 and allow for density modifications assuming r_c = 0.2. For J6-J10, the discrepancy between measurements and uniform rotation models is large for all models that assume uniform rotation. Column 5 shows a representative model with DR on cylinders and a deep rotation period of 10h39m22s that matches measurements from J2 to J10." Values x10^6.

| | Measurements | Physical models, uniform rotation (min, max) | Uniform rotation with modified density profiles (min, max) | Physical model with differential rotation |
| --- | --- | --- | --- | --- |
| J2 | 16290.573 +/- 0.028 | 16290.57 | 16290.57 | 16290.573 |
| J4 | -935.314 +/- 0.037 | -935.31 | -990.12 to -902.93 | -935.312 |
| J6 | 86.340 +/- 0.087 | 80.74 to 81.76 | 75.69 to 90.42 | 86.343 |
| J8 | -14.624 +/- 0.205 | -8.96 to -8.70 | -10.26 to -7.97 | -14.616 |
| J10 | 4.672 +/- 0.420 | 1.08 to 1.13 | 0.97 to 1.33 | 4.677 |

### Table 3 (p14). READ caption: "Contribution to the higher gravity harmonics dJ8 and dJ10 resulting from differential rotation and thermal-wind optimization. The deviation (Column 1) is the difference between the measured J8 and J10 (Table 1) and the average of the computed values from the 11 CMS models with uniform rotation (Table 2). Two optimizations are shown: one without latitudinal truncation of the zonal flow, resulting in the reconstructed zonal wind profile shown in Fig. 4A and with a flow depth of 9363 km (Column 2), and the second with the flows truncated at latitude 60 deg (Fig. 4B) and a flow depth of 8832 km (Column 3). Columns 4 and 5 show the deviations calculated with the thermal-gravity equation (48) for similar wind profiles. The solutions from thermal wind are closer to the measurement because the optimization was done using the thermal wind method, but the thermal-gravity solutions also match the observations within 10%." Values x10^6.

| | Deviation | Thermal-wind solution | Thermal-wind solution truncated at latitude 60 deg | Thermal-gravity solution | Thermal-gravity solution truncated at 60 deg |
| --- | --- | --- | --- | --- | --- |
| dJ8 | -5.600 +/- 0.205 | -5.624 | -5.533 | -5.758 | -5.759 |
| dJ10 | 3.528 +/- 0.659 | 3.570 | 3.660 | 3.974 | 4.037 |

### Supplementary tables
Tables S1 and S2 and Figs. S1 to S6 are not in the held file. Text-only fragments from the main text: Table S1 holds the orbit table and Doppler rms by arc (p1-2 references "table S1"); Table S2 holds the ring mass estimates and a priori values (p2 "see table S2": initial A, B, C ring masses from ring occultation data with 100 percent a priori uncertainty; the B ring had a priori uncertainty of 10 Mimas masses).

### Figures (not tables; for orientation only)
Fig. 1 (p11): zonal harmonics J2 to J12 against degree, uniform rotation predictions diverging from measurements above degree 8. Fig. 2 (p11): differential rotation profiles from CMS models against observed cloud-level profiles. Fig. 3 (p12): composition parameters of interior models (core mass 15 to 18 Earth masses). Fig. 4 (p12): observed and reconstructed wind profiles.

## 3. Interpretation, in brief (numbers quoted from the text)

- Uniform rotation is ruled out (p4-5): models matching J2 and J4 cannot reproduce J6 to J10; Table 2 shows the observed J8 = -14.624 against -8.96 to -8.70 and J10 = 4.672 against 1.08 to 1.13. READ (p5): "The inability to reproduce the unusually large values of J6, J8, and J10 leads us to conclude that models with uniform rotation are ruled out by the gravity data."
- Differential rotation on cylinders (CMS+DR, p5-6): "By assuming Saturn's equatorial region rotates approximately 4% faster than the deep interior, agreement between models and data for all coefficients J2-J10 can be achieved". Core masses between 15.0 and 18.2 Earth masses fit best; the mass of heavy elements in the envelope is 1.3 to 4.8 Earth masses. A minimum rotation region near l = 0.83 (about 10,000 km from the surface at the equator).
- Differential rotation with finite depth (thermal-wind method, p6-7): flow depth 9363 +/- 357 km (full wind profile) and 8832 +/- 295 km (truncated at latitude 60 deg). READ (p7): "Regardless of the exact meridional profile, all vertical flow profiles are constrained to contain a very deep flow of about 9000 km. This depth, corresponding to 15% of Saturn's radius ... suggesting that the flow should extend down to the levels of magnetic dissipation." Extending the cloud-level flow into the interior gives dJ8 about -1.5e-6 and dJ10 about 1e-6, "at least a factor of two too small" (p6); the observed deviations are dJ8 = -5.600e-6 and dJ10 = 3.528e-6 (Table 3).
- Ring mass (p3, p7): the data constrain the sum of the A, B, C ring masses; the individual masses are poorly determined. Total ring GM = 1.02 +/- 0.41 km^3/s^2 (p3) "equivalent to 0.41 +/- 0.13 Mimas masses". On p7: B ring GM 0.58 +/- 0.48 km^3/s^2, A 0.38, C 0.06, total 1.02 +/- 0.43 km^3/s^2 (the p3 and p7 uncertainties differ, 0.41 against 0.43, as printed). READ (abstract): (1.54 +/- 0.49)x10^19 kg. COMPUTED consistency: 1.02 / 2.5026 = 0.408 Mimas masses. Voyager value 2.8x10^19 kg (0.75 Mimas masses) is larger; a combined density-wave estimate before this paper was GM = 1.01 or 0.40 Mimas masses (p7). Numerical self-gravity wake simulations suggest an upper limit of 9.7x10^19 kg or about 2.5 Mimas masses (p7).
- Ring age (p7-8): from the low mass, evolutionary ages around 10^8 yr (A and B rings, "assuming the Voyager measurements of ring masses and interplanetary impact fluxes"); pre-Cassini estimates of 80-150 Myr (A) and 30-100 Myr (B); the revised B ring mass would increase the latter by about 25 percent. READ (p8): "On balance, we favor a scenario whereby the present rings of Saturn are relatively young, at least compared to the planet itself, although they may have evolved substantially in the past 10^7-10^8 years and were perhaps once more massive than they are today." Viscous spreading models approach an asymptotic mass of about 1.5x10^19 kg or 0.40 Mimas masses after 5 Gyr (p8).

## 4. Constants audit against this paper

Project sources: `src/cyclerfinder/data/validation/v4_saturn.py` lines 84 to 110 and `src/cyclerfinder/core/satellites.py` (`PRIMARIES["Saturn"]` line 67; Saturnian entries lines 219 to 231 and 301 to 305). Differences COMPUTED (project minus paper).

| Quantity | Project value | Paper | Difference | Notes |
| --- | --- | --- | --- | --- |
| `SATURN_J2` | 16290.573e-6 | 16290.573 +/- 0.028 (x10^-6), Table 1 | 0 | Identical to all printed digits; the 0.028 is 3 times the formal uncertainty |
| `SATURN_J2_REF_RADIUS_KM` | 60330.0 | "reference radius 60330 km" (Tables 1 and 2 captions) | 0 | READ |
| `SATURN_R_EQ_KM` | 60268.0 | not printed | n/a | The paper's only radius is the 60,330 km reference radius; 60,268 km is the project's IAU/JPL value, not in this paper (and the 60,330 km is a reference radius, not stated to be a physical radius) |
| `SATURN_J2_AT_R_EQ` | J2 x (60330/60268)^2 | derived | COMPUTED: 16290.573e-6 x (60330/60268)^2 = 16324.108e-6 | Ratio (60330/60268)^2 = 1.0020585 (0.206 percent, matching the #894 value 0.21 percent). The arithmetic is right; the paper contains nothing against it |
| `PRIMARIES["Saturn"]` GM | 3.7931207e7 km^3/s^2 | not printed in the main text (GM is estimated, p2, value in the supplement or reference 18, not held) | n/a | Source remains the project's stated JPL DE440 value; not checkable here |
| Mimas GM | 2.503 | 2.5026 (p2) | +0.0004 (0.016 percent) | Paper value is quoted as the ring mass unit; the project's 2.503 is consistent to its rounding |
| Mimas radius, a | 198.2 km, 185540 km | not printed | n/a | |
| Enceladus, Tethys, Dione, Rhea, Titan, Iapetus, Hyperion (GM, radius, a, flyby floors) | various | not printed | n/a | None of these moon constants appears in the held file. The paper names SAT389 and SAT393 as the satellite ephemerides used (p2) and DE430; it does not tabulate them |
| J4, J6 and higher | not carried | J4 = -935.314e-6, J6 = 86.340e-6 (Table 1) | n/a | See section 6 |
| Ring mass | not carried | GM 1.02 km^3/s^2 | n/a | |

## 5. The tidal term in J2 and what it means for propagation

READ (Table 1 caption, p13): "The J2 value includes a constant tidal term owing to the average tidal perturbation from the satellites." READ (p2): the paper's own dynamical model contained both the explicit "point-mass gravitational accelerations from Saturn and its satellites (including the ring moons)" and "Saturn's tidal response to its satellites (12)", with Love number k22 as a consider parameter, and J2 as an estimated global parameter. So the published J2 is the number that this paper's fit used together with explicit moon point masses.

INFERRED: using the Table 1 J2 together with explicit moon point masses in a propagator (the project's v4_saturn lane, J2 plus eight moons as third bodies) is the same arrangement the paper's own fit used, so there is no inconsistency in the pairing in kind. The "constant tidal term" means the J2 is the observed effective J2 (the permanent shape plus the time-averaged tidal bulge raised by the moons), not a rotation-only J2; adding the moons as point masses does not double count it except at the level of the paper's separate tidal-response (k22) acceleration, which the project does not model. The project's omission of that tidal-response acceleration is a model simplification; its size is not computed here and the paper does not give it. Whether the omission matters for the project's accuracy needs a separate estimate (INFERRED, not computed).

## 6. Do J4 and higher terms matter at the project's Saturn distances

Setup (COMPUTED, standard zonal potential). For a particle in Saturn's equatorial plane at radius r, the radial acceleration from the zonal field is `g = -(mu / r^2) [1 - sum_n (n+1) J_n (R/r)^n P_n(0)]` with R = 60,330 km (the reference radius of Table 1), P_2(0) = -1/2, P_4(0) = 3/8, P_6(0) = -5/16. Hence the J2 term is `+(3/2) J2 (R/r)^2`, the J4 term `-(15/8) J4 (R/r)^4` (positive because J4 is negative), the J6 term `+7 x (5/16) J6 (R/r)^6`, in units of mu/r^2. Ratios of the J4 to the J2 acceleration:
`a4 / a2 = (5/4) |J4/J2| (R/r)^2 = 1.25 x 0.057413 x (R/r)^2 = 0.071766 (R/r)^2`, with J4/J2 = 935.314 / 16290.573.

| Orbit | r (km) | R/r | a2 / central | a4 / central | a6 / central | a4 / a2 | a6 / a2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Titan (registry a) | 1,221,870 | 0.049377 | 5.96e-5 | 1.04e-8 | 2.7e-12 | 1.75e-4 | 4.6e-8 |
| Iapetus (registry a) | 3,561,700 | 0.016938 | 7.01e-6 | 1.44e-10 | 4.5e-15 | 2.06e-5 | 6.4e-10 |
| Rhea (registry a, for reference) | 527,070 | 0.114464 | 3.20e-4 | 3.01e-7 | 4.2e-10 | 9.4e-4 | 1.3e-6 |

Conclusion (INFERRED from the COMPUTED ratios): across Titan to Iapetus the J4 acceleration is 2e-5 to 2e-4 of the J2 acceleration and 1e-10 to 1e-8 of the central attraction; J6 and higher are smaller still by a further factor of 1e3 or more. The J2 uncertainty itself (0.028e-6 on 16290.573e-6, 1.7e-6 relative) is far smaller than the J4 effect, so J4 would matter before the J2 uncertainty does, but at 1.75e-4 of J2 at Titan it is the same size as a 0.009 percent error in J2 R^2, well below the 0.206 percent radius-pairing error corrected in #894. Whether it is below the project's other model omissions (moons as third bodies, Saturn's tidal response) was not computed here. At Rhea-Tethys distances (about 5e5 km or less) J4 reaches 1e-3 of J2.

## 7. Items needing care

- The first-release text has two values for the ring mass uncertainty (0.41 and 0.43 km^3/s^2; READ, p3 and p7).
- The citation in `v4_saturn.py`, "Science 364(6445), 1052-1056", cannot be checked from the held file (first-release pagination). The DOI is correct as printed on p1.
- The J2 uncertainty is 3 times the formal value; do not treat 0.028 as a 1-sigma figure.
- Saturn's GM, the Saturnian moon constants and the Saturn ring masses by ring (the individual A, B, C values are text-only) are not in the held file.
