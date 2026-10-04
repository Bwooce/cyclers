# Digest: Casoliva, Mondelo, Villac, Mease, Barrabes & Olle 2008, "Families of cycler trajectories in the Earth-Moon system" (AIAA 2008-6434), seeds and method (#920)

AIAA/AAS Astrodynamics Specialist Conference, Honolulu, 18-21 Aug 2008. Filed in the private paper corpus as
`casoliva-mondelo-villac-mease-barrabes-olle-2008-families-cycler-trajectories-earth-moon-AIAA-2008-6434.pdf` (20 pages; a `.txt` text layer is alongside; `pdftotext | wc -w` = 10021).
Status: READ (all 20 page images; the Table 2 and Table 3 digits were also compared with the `.txt` layer, one digit differs between my first image reading and the text layer, see Table 2 note). COMPUTED = my arithmetic 2026-10-04.
Page references are "p.N of 20" of the AIAA PDF. The journal version (JGCD 33(5):1623-1640, 2010) is digested in `docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md`; this note adds what #725 did not extract:
the seed table's exact content, the seed method as printed in 2008, what the printed numbers verify, and the build recipe. For the seed formulas' origin see `docs/notes/2026-10-04-digest-barrabes-gomez-2002-spatial-pq-resonant-orbits.md`
and `docs/notes/2026-10-04-digest-barrabes-gomez-2003-pq-resonant-second-species.md`. To read the 2010 section IV.C and its Table 3 I also read journal pages 1626-1630 of the JGCD paper (PDF pages 4-8 of the 2010 file); quotes marked [2010].

## 1. Model and conventions (pp.2-4)
- Planar circular restricted three-body problem (PCR3BP). Earth at (mu, 0), Moon at (mu - 1, 0) (Moon on the left), "mu = 0.01215 is the Moon-Earth mass ratio" (p.3) for the Earth-Moon system. Units: separation 1, period 2 pi.
- Eq. 1-2: dX/dt = f(X), X = [x y u v]^T, f = [u, v, 2v + dOmega/dx, -2u + dOmega/dy]^T, Omega = (x^2 + y^2)/2 + (1 - mu)/r_1 + mu/r_2, r_1 = sqrt((x - mu)^2 + y^2), r_2 = sqrt((x - mu + 1)^2 + y^2).
- Eq. 3: C_J(X) = 2 Omega - ||V||^2, V = (u, v)^T. Eq. 4-5: Hamiltonian form with p_x = u - y, p_y = v + x, H = (p_x^2 + p_y^2)/2 - x p_y + y p_x - (1 - mu)/r_1 - mu/r_2; "C = -2H".
- Eq. 6-8: variational equation; monodromy eigenvalues {lambda, 1/lambda, 1, 1}; stability index k = lambda + 1/lambda, stable iff |k| <= 2.
- II.C continuation (p.4): unknowns mu, h, T, X_0 = (x_0, p_x, y, p_y) (planar), equations H(x_0) - h = 0, g(x_0) = 0 (Poincare section), phi_T(x_0) - x_0 = 0: six equations, seven unknowns; fixing one of mu, h, T follows a family with that fixed.
  Verbatim: "The Stromgren-Wintner's natural termination principle states the possibility of analytic continuation of any one-parameter family of periodic orbits in the restricted three body problem as long as the solutions stay clear of a collision singularity and their periods remain bounded."
- Criteria for a cycler (p.4): apogee radius, perigee radius, periselene radius, velocities at perigee and periselene, period, stability index k, energy (Jacobi constant), backside-of-the-Moon coverage.

## 2. Resonance constraints (section IV.A, pp.5-6) and the labelling conflict
Eq. 9, verbatim: "We say that the spacecraft is in a p-q resonant orbit when p T_M = q T_s/c where p and q are relative prime integers, T_M and T_s/c represent the period of the Moon and the spacecraft when defined in a sidereal frame."
Prose: "such resonance relations imply that the spacecraft will complete q times its (inertial) elliptic orbit while the Moon will have completed p revolutions around the Earth". Eq. 10: a_s/c = (p/q)^(2/3), E_s/c = (mu_E/4)(q/p)^(2/3) (mu_E = 1).
CONFLICT (COMPUTED, section 6.3): in Eq. 9, Eq. 10, Table 1 and Figure 1 the roles are p = Moon revolutions, q = spacecraft revolutions. But the designations of the Table 2 orbits and Eq. 15-16 (copied from Barrabes & Gomez) use the OPPOSITE roles:
p = spacecraft revolutions, q = Moon revolutions, a = (q/p)^(2/3). The 2010 JGCD paper uses the Barrabes & Gomez roles throughout (its Table 1, rows q and columns p, a_s = (q/p)^(2/3); [2010] p.1626 Eq. 10 prints a_s = mu_E^(1/3)(q/p)^(2/3)).
2008 Table 1 is the TRANSPOSE of 2010 Table 1. Rule for the project: use the Barrabes & Gomez / 2010 roles; never use the 2008 Eq. 9 reading for anything that touches the designations.
Consequence for a sentence in the 2008 conclusions (p.19): "As a rule of thumb, a resonant ellipse with semi-major axis less that 0.6 the Earth-Moon distance may not exist in the Earth-Moon system. This includes the 2:5, 3:7 or 3:8 resonances." In 2008 labelling "3:7"
means a = (3/7)^(2/3) = 0.568, the orbit that 2010 reaches and catalogues as 7-3 (a = (3/7)^(2/3), designation 7-3a to 7-3c; "tight 7-3 cycler"). So the 2008 conclusion that it may not exist was reversed by the 2010 seed-and-continuation method ([2010] p.1630: "a 7-3 resonant orbit for small mu can be continued
to a 7-3 resonant orbit for the Earth-Moon mass ratio. Tight cyclers cannot be computed via the differential correction of an elliptical orbit.").

### Table 1 (p.5), as printed (rows p, columns q; spacecraft semi-major axis = (p/q)^(2/3) in the 2008 labelling; grey cells are forbidden resonances)
| p \ q | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.000 | 0.630 | 0.481 (grey) | | | | | | |
| 2 | 1.587 | 1.000 | 0.763 | 0.630 | 0.543 | 0.481 (grey) | | | |
| 3 | 2.080 | 1.310 | 1.000 | 0.825 | 0.711 | 0.630 | 0.568 | 0.520 | 0.481 (grey) |
| 4 | 2.520 | 1.587 | 1.211 | 1.000 | 0.862 | 0.763 | 0.689 | 0.630 | 0.582 |
| 5 | 2.924 | 1.842 | 1.406 | 1.160 | 1.000 | 0.886 | 0.800 | 0.731 | 0.676 |
COMPUTED check: (p/q)^(2/3) for all printed cells agrees to the printed three decimals (the cell p = 4, q = 9 is 0.5824, printed 0.582; the others are rounded). The three grey cells are exactly the ones with q > 2.83 p ((1,3), (2,6), (3,9)).
Footnote: "The value of semi-major axis for non-relative prime integers has been given for ease of reference ... a 6-3 resonance is really a 2-1 resonance."
Eq. 11 (p.6), verbatim: "Taking epsilon ~ 0.02289, that is, allowing ~300 km altitude at closest approach from the Earth and ~70 km around the Moon, the physical constraint on the resonances is then q < p (2/(1 + epsilon))^(3/2), that is, q <~ 2.8284 p."
COMPUTED: (2/1.02289)^1.5 = 2.7340, while 2^1.5 = 2.8284 is the epsilon = 0 value. The printed "2.8284" is the epsilon = 0 limit; the epsilon value given would yield 2.734. The 2010 paper prints "p <= 2.7372911990 q" with epsilon = 0.0220747554 ([2010] p.1626),
consistent with its different epsilon (computed 2.7373 when (2/1.0220747554)^1.5 is evaluated; COMPUTED). Use the 2010 numbers.
Eq. 12-13 (p.6): apogee r_a = 2 (p/q)^(2/3) - r_p, r_a <= (p/q)^(2/3) - R_E with "R_E ... on the order of 7000 km"; e between (R_M/2)(p/q)^(2/3) and 1 - (R_E/2)(p/q)^(2/3) (2008 labelling, so under the Barrabes & Gomez labelling swap p and q).

## 3. Two-body seeds and their failure (section IV.B, pp.7-8)
Method: ignore the Moon, take the Kepler ellipse of the resonant semi-major axis, differentially correct at fixed period T at the physical mass; examples 1-2, 2-3, 3-4 (Figure 2); continuation in C_J (Figure 3) at fixed mu.
Verbatim limits (p.7): "Note that while the above approach is sufficient for cyclers with p > q, it may not provide a sufficiently good approximation for resonances with smaller semi-major axis (p < q). This is especially true for resonances with small semi-major axis such as the 3 : 7 and 3 : 8 resonances, for which the apogee is near the Moon.
From our experience, these approximations fails when the semi-major axis is on the order of 0.6 or below. This approach does not allow also a fine control on the longitude of perigee as this quantity generally drifts during the correction process." (2008 labelling "3:7": a = 0.568.)
Alternative (p.8): "start from the two-body approximation and continue the solution as a function of the mass parameter, mu. This approach leads to better success when the semi-major axis is still sufficiently large ... fails for limiting resonances, such as the 3:8 resonances for which the apogee must lie very close to the Moon."

## 4. The seed method (section IV.C, pp.8-10): what is printed
Premise (p.8): "Trajectories that encircle both the Earth and Moon are part of families which pass through the Moon and are known as second species solution in the terminology of Poincare. Second species solutions are thus periodic orbits generated from segments of two-body orbits as the mass parameter mu is let to tend to zero. Henon provides a systematic classification of such arcs ... for actual continuation, one needs numerical values of initial conditions to generate these periodic orbit families. For these, we use the work of Barrabes and Gomez."
Eq. 14 (p.9), initial conditions on a circle C of radius mu^alpha about the Moon:
  x_i = mu - 1 + mu^alpha cos(theta),  y_i = mu^alpha sin(theta),  u_i = v cos(psi),  v_i = v sin(psi)                                       (14)
(the symbol v is both the speed and the second velocity component's name; v_i is the y-velocity.)
Eq. 15: for p, q in N with p/q <= 2 sqrt 2, "there is a family of p-q resonant orbits with the Jacobi constant C_J within the interval"
  (p/q)^(2/3) - 2 sqrt( 2 - (p/q)^(2/3) ) < C_J < (p/q)^(2/3) + 2 sqrt( 2 - (p/q)^(2/3) )                                                     (15)
Eq. 16: "For each value of C_J within the above set, there are two admissible values of the polar angle of the velocity psi given by the relation"
  sin(psi) = ( 2 - C_J + (p/q)^(2/3) ) / ( 2 sqrt(3 - C_J) )                                                                                (16)
Eq. 17: "for every admissible pair (C_J, psi), there exists one value of theta such that the initial conditions in Eq. (14) satisfy the corresponding matching conditions between the in and out maps. Hence, theta must satisfy the nonlinear constraint below
  (2/sqrt(3 - C_J)) ( cos(theta) sin^2(theta - psi) - cos(psi) cos(theta - psi) ) + sin^3(theta - psi) = 0                                  (17)
such that cos(theta - psi) > 0, i.e., the spacecraft is moving away from circle C. Finally, the total time required to come back close to the initial condition defined in Eq. (14) is
  T = 2 pi q + O(mu^alpha)                                                                                                                (18)".
Mapping to the sources (READ + COMPUTED): Eq. 15 = Barrabes & Gomez 2002 Eq. 27 (and 2003 Eq. 10); Eq. 16 = 2003 Eq. 46 (planar form of 2002 Eq. 26 with phi_0 = 0); Eq. 17 = 2003 Eq. 45, since 1/sqrt(c) = 2/sqrt(3 - C_J) with c = (3 - C_J)/4 (COMPUTED identity);
Eq. 18 = 2002/2003 Eq. 4 up to the sphere-transit correction (2003 digest section 6). Eq. 14 = 2002/2003 Eq. 16/5 with varphi = 0 (planar). The 2008 paper never states alpha; it states only "alpha in (1/3, 1/2)" (p.9, "with an error of the order mu^(1 - alpha), for some alpha in (1/3, 1/2)"); the only printed alpha anywhere in the three papers is 0.4 (Barrabes & Gomez examples).
Corrected orbits (p.9), verbatim: "Since these approximate initial conditions are only valid for very small values of mu and do not yield perfectly periodic orbits, the initial conditions had to be differentially corrected using a grid of C_J values. Only a few initial conditions per p-q pair at mu = 10^-6 yield periodic orbits.
Differentially correcting such orbits presented several challenges. First, the small value of the mass parameter requires significant accuracy in the integration and very small step in the continuation parameter." (2010 adds: grid of 500 values of C_J, periodicity error below 1e-10, double precision, RKF 7-8 with local error below 1e-13, CPU 2.85 h (1-2), 5.66 h (2-1), 8.03 h (3-2), 23.49 h (7-3) for all 500, 2.6 GHz Opteron, MATLAB [2010 p.1628].)
Designation (p.9): "The designation of the cycler trajectories obtained via p-q resonant orbits is composed of 2 digits (the p-q values) followed by a letter (from 'a' to 'z'). For given p-q resonant orbits, there are only a discrete set of resonant orbits that will be periodic and, thus, the letter 'a' would represent the cycler trajectory within a p-q family that has the lowest Jacobi constant value."

### Table 2 (p.10), the seed table at mu = 1e-6, transcribed digit by digit
Caption (verbatim): "Cycler designation, p-q values, Jacobi constant C_J, period and initial conditions, application and stability index for several cycler trajectories". The caption does not state the mass ratio. COMPUTED (section 6.1): the printed states reproduce the printed C_J to 1e-16 at mu = 1e-6 and fail by 0.7 to 9 at mu = 0.01215, so the table is at mu = 1e-6 (consistent with p.9 text and the Figure 4 caption "mu = 1e-6").
Columns: designation; p-q; C_J (first line) and T (second line); X_i = (x, y) (initial position, x and y); V_i = (u, v) (initial velocity, x and y); k. First line of each row is C_J, x, u; second line is T, y, v. y = 0 for every row, so the initial point is a crossing of the section y = 0.
| Desig. | p-q | C_J | T | x_i | y_i | u_i | v_i | k |
|---|---|---|---|---|---|---|---|---|
| 12a | 1-2 | -0.4048949508278787 | 12.5729004699819580 | -0.9997842236429277 | 0.0000000000000000 | -1.0475686168407430 | 1.5221048236500363 | 994.8214 |
| 21a | 2-1 | 0.3044238301466371 | 6.2807379868905375 | -0.9994542367695188 | 0.0000000000000000 | 0.0000000000159674 | 1.6429377286480178 | 2.0214 |
| 23a | 2-3 | -1.4624706218555543 | 18.8645742008117736 | -0.9997040086087932 | 0.0000000000000000 | -0.0000000000029493 | 2.1140593042383320 | 1.8405 |
| 23b | 2-3 | -0.3856270265962789 | 18.8497995179817757 | -0.9953809710199844 | 0.0000000000000000 | -0.9554688344071062 | 1.5726409617865593 | -0.7753 |
| 32a | 3-2 | -0.3403221450450835 | 12.5363376944500722 | -0.9999035023356472 | 0.0000000000000000 | 0.0000000000131068 | 1.8333742370247279 | 6.3818 |
| 32b | 3-2 | 2.0635340336761394 | 12.5660196280208911 | -1.0188148478462549 | 0.0000000000000000 | -0.7183851038146349 | 0.6492612209948345 | 1.5681 |
| 52a | 5-2 | 1.0461882704974470 | 12.5651492405106922 | -0.9994423797258251 | 0.0000000000000000 | 0.0000000000459148 | 1.3990717541095201 | 2.0366 |
| 54a | 5-4 | -0.5902501452788234 | 25.1304852528305673 | -0.9988849982363450 | 0.0000000000000000 | -0.2709574260065351 | 1.8758004354491455 | -1.3192 |
| 54b | 5-4 | -0.6598717597930506 | 25.1321450585584110 | -0.9905419593393083 | 0.0000000000000000 | -0.1198961389475177 | 1.9094434196183308 | 1.9943 |
| 73a | 7-3 | 0.8957501590757784 | 18.8492803402344329 | -0.9954265899784440 | 0.0000000000000000 | -0.2486030886355384 | 1.4293154529931373 | 1.8799 |
Digit note: the PDF image of 21a x_i is hard to read at the 11th decimal; the `.txt` text layer gives -0.9994542367695188, which I use (my first image reading had a different digit there). All other digits agree between image and text layer for the cells I compared.
Footnotes: x_i, y_i initial position; u_i, v_i initial velocity; k as Eq. 8. The "application (communication/navigation or transportation)" column promised in the caption does not appear in the table.
Figure 4 (p.11): the nine orbits plotted (12a, 21a, 23a, 23b, 32a, 32b, 52a, 54a, 73a; the printed panels are (a) 12a, (b) 21a, (c) 23a, (d) 23b, (e) 32a, (f) 32b, (g) 52a, (h) 54a, (i) 73a; 54b is in the table but is not a panel), "Periodic orbits computed for different p-q resonances with mu = 1e-6".
Table 2 has ten rows and Figure 4 nine panels (my count of the page).

## 5. Continuation in mass (pp.9-12), as printed in 2008, versus 2010
2008 text (p.9-10), verbatim: "Then, to look for resonant periodic orbits in the Earth-Moon system, continuation in the mass parameter mu was performed at fixed T. While this continuation is the most natural to preserve the resonant character of these orbits, it leads in most cases to drive the periodic orbits to the singularity at the center of the Moon, showing that the initial resonance relation on these particular families is in fact not preserved. However, such resonance relation may be present in families with different phasing and presenting more distant fly-bys.
For example, the cycler shown on Fig. 4(f) cannot be continued to the mass parameter of the value of the Earth-Moon system, even though it is very similar to the 2:3 resonance presented in Fig. 2. These cyclers are in fact different only by their longitude of perigee." (Fig. 4(f) is 32b: C_J = 2.0635, perigee/apogee shifted.)
Then (p.10): "Another example is presented in Figs. 5 and 6, where the orbit represented in blue has been generated with the above patched conic approximation while the red orbit correspond to its continuation at fixed period as a function of T [sic; of mu]. The longitude of the perigee is shifted by 90 deg and the 'cycler' looses actually its property of fly-by the backside of the Moon."
Order in 2008: "As with the 2-body approximation, one can also continue these families in C_J so as to increase the periselene distance. Then for a family member sufficiently far from the Moon, we can continue the family at fixed C_J until the mass ratio matches the value for the Moon. The resulting periodic orbit can then be continued at fixed mu by varying C_J and investigate the various higher order resonance relations present. As opposed to the pure classification of second species solutions, not all the allowable resonance relation indicated by the two-body approximation are allowed for cyclers in the Earth-Moon systems."
2008 order: (1) C_J continuation at mu = 1e-6 to raise the periselene; (2) mu continuation at fixed C_J to mu_M; (3) C_J continuation at mu_M to find resonances.
2010 order ([2010] p.1629): "To continue a (p, q) resonant orbit from small mu to mu_M, the direct strategy is to fix the period T, to preserve the resonance, while varying mu. However, in most cases, such continuation leads to a periodic orbit that impacts the moon (especially for resonant orbits with high p/q), showing that resonant cyclers with the seed (p, q) do not exist uniformly in mu in the family. Alternatively, mu continuation at fixed C_J does not preserve the resonance relation, but allows the continuation of cyclers from small mu to mu_M in many cases, but for cycler trajectories that fly by the back of the moon, such continuation leads to an orbit that impacts the moon. Instead a three-step continuation strategy was necessary. When the mu continuation begins producing orbits leading to lunar impact, we switch to C_J continuation with mu fixed so as to increase the periselene distance; then the mu continuation with fixed C_J is resumed until mu = mu_M."
"The continuation strategies just discussed produce periodic orbits with mu = mu_M, but because the strategies do not preserve the (p, q) resonance of the seed, there is a final step. For each periodic orbit, which has a particular C_J or the period T, a segment of a characteristic curve of the periodic orbit family (family segment for short) is generated by C_J continuation, in both directions, with fixed mu = mu_M, until the family has a natural termination (i.e., one or more of the following quantities grows without limit: dimension of the orbit, C_J, and T). ... If the family segment intersects the p/q of the seed, indicated by a thin horizontal line in Fig. 3, then there is a resonant orbit for mu = mu_M with the (p, q) of the seed. For example, in Fig. 3d, there are four family segments for the seed 7-3 resonance, but only three intersections with the 7/3 line; thus, only three of the continuations led to orbits with the 7-3 resonance." [2010 p.1629; wording checked against the text layer; the paper also says the periodic-orbit family in (mu, C_J, T) is "a two-parameter family", so a one-parameter continuation must choose a path in the (mu, C_J) or (mu, T) space.]
UNCLEAR (2010): in the first step of the three-step strategy the printed sentence does not say which of T or C_J is held fixed during the first mu continuation; the sentences before it name the two options (fixed T: hits the Moon; fixed C_J: does not preserve the resonance and hits the Moon for back-of-the-Moon cyclers), and the three-step strategy "begins" with whichever mu continuation produced the lunar impact. INFERRED: first leg at fixed T (the resonance-preserving choice), the resumed leg at fixed C_J ("the mu continuation with fixed C_J is resumed").
Figure 5-6 (p.12): argument of perigee shift with mu continuation at fixed T; Figure 6 gives x_i, y_i, u_i, v_i, T, k against mu from 0 to 0.012 (k about 3.4 to 3.9 on the plot axes; the orbit's C_J is not stated; qualitative).
Cost, [2010] p.1629-1630: "the average computation time for the cyclers with p = 2 and q = 1 was 35-40 min"; 2010 ends (p.1630): "the approach of continuing from the small-mu perturbation approximation has a couple of benefits. One is that it is more systematic and thorough in finding resonant orbits for mu = mu_M. The other is to compute resonant tight cyclers, meaning cyclers that have both small perigee and periselene distances. We have shown that a 7-3 resonant orbit for small mu can be continued to a 7-3 resonant orbit for the Earth-Moon mass ratio."

## 6. Checks from printed numbers (COMPUTED, Python double precision)
### 6.1 Mass ratio of Table 2 (instrument check)
Recomputing C_J = 2 Omega - (u^2 + v^2) from each printed state:
- at mu = 1e-6: difference from the printed C_J is below 6e-17 for all ten rows (the printed digits are reproduced to the last printed digit);
- at mu = 0.01215: differences 1.98, 2.04, 2.00, 3.18, 1.95, 0.74, 2.04, 2.15, 8.98, 3.16 (rows in table order), nowhere near zero.
So Table 2 is at mu = 1e-6 (Moon at x = mu - 1 = -0.999999). The Moon is at distance r_2 from the state: 2.15e-4, 5.45e-4, 2.95e-4, 4.62e-3, 9.55e-5, 1.88e-2, 5.57e-4, 1.11e-3, 9.46e-3, 4.57e-3 (12a, 21a, 23a, 23b, 32a, 32b, 52a, 54a, 54b, 73a): each state is a y = 0 crossing close to the Moon, in the second-species regime.
### 6.2 Semi-major axis and period of the printed states
Earth-centred inertial speed (u, v + x - mu) and distance |x - mu| give a = 1/(2/r - V^2):
| row | a computed | (q/p)^(2/3), B&G labelling | (p/q)^(2/3), 2008 Eq. 10 labelling | T - 2 pi q | C_J - C_J1 (Eq. 15 lower end) | r_2 |
|---|---|---|---|---|---|---|
| 12a | 1.58677 | 1.58740 | 0.62996 | +6.530e-3 | +1.306 | 2.148e-4 |
| 21a | 0.63011 | 0.62996 | 1.58740 | -2.447e-3 | +1.700e-3 | 5.448e-4 |
| 23a | 1.31786 | 1.31037 | 0.76314 | +1.502e-2 | -1.332e-3 (outside) | 2.950e-4 |
| 23b | 1.31039 | 1.31037 | 0.76314 | +2.436e-4 | +1.076 | 4.618e-3 |
| 32a | 0.76598 | 0.76314 | 1.31037 | -3.003e-2 | +1.019e-2 | 9.550e-5 |
| 32b | 0.76312 | 0.76314 | 1.31037 | -3.510e-4 | +2.414 | 1.882e-2 |
| 52a | 0.54306 | 0.54288 | 1.84202 | -1.221e-3 | -8.828e-4 (outside) | 5.566e-4 |
| 54a | 0.86219 | 0.86177 | 1.16040 | -2.256e-3 | +8.195e-2 | 1.114e-3 |
| 54b | 0.86182 | 0.86177 | 1.16040 | -5.962e-4 | +1.233e-2 | 9.457e-3 |
| 73a | 0.56846 | 0.56844 | 1.75921 | -2.756e-4 | +0.1179 | 4.572e-3 |
Reading: the Table 2 designation p-q means p spacecraft revolutions about the Earth per q Moon revolutions, a = (q/p)^(2/3), T close to 2 pi q: the Barrabes & Gomez labelling, not the 2008 Eq. 9/10 one. The deviation of a from the Keplerian value grows as r_2 shrinks (23a, 32a), as expected for a point near the Moon. Interval (Eq. 15) values at mu = 1e-6: 1-2: (-1.711013, 2.970934); 2-1: (0.302724, 2.872078); 2-3: (-1.461139, 2.987424); 3-2: (-0.350508, 2.971249); 5-2: (1.047071, 2.636960); 5-4: (-0.672200, 2.992994); 7-3: (0.777805, 2.740616).
23a (C_J = -1.46247) and 52a (C_J = 1.04619) lie just OUTSIDE the first-order interval by 1.3e-3 and 8.8e-4: allowed, since the seed is only O(mu^(1-alpha)) accurate (1e-6^0.6 = 2.5e-4 times a constant).
### 6.3 The four psi = 90 degree rows
21a, 23a, 32a, 52a have u_i of order 1e-11 (perpendicular crossings: psi = 90 degrees), and sin(psi) = 0.999798, 1.000166, 0.998737, 1.000090 from Eq. 16, i.e. they sit at the lower end C_J1 of Eq. 15 where sin psi = 1 exactly (the C_J endpoint, excluded by Barrabes & Gomez 2003 p.149 because cos^2 phi_0 sin^2 psi_0 = 1 there). The other six rows (12a, 23b, 32b, 54a, 54b, 73a) have u_i not equal to zero.
Direction from Eq. 16 for the six: psi branches 55.3206 / 124.6794 (12a, table 124.5372), 58.8306 / 121.1694 (23b, table 121.2811), 40.1071 / 139.8929 (32b, table 137.8934), 81.7788 / 98.2212 (54a, table 98.2195), 86.8206 / 93.1794 (54b, table 93.5930), 80.7476 / 99.2524 (73a, table 99.8668).
The table value of psi lies on the second branch (psi > 90 degrees) in all six rows, agreeing with the branch value to 0.1 to 2 degrees.
INFERENCE for P2: a uniform 500-point grid strictly inside the interval will rarely land a seed on the symmetric (perpendicular) crossings that sit within about 1e-2 of the lower end; the endpoint neighbourhood must be scanned, and the corrector must tolerate sin(psi) slightly above 1 (clip it to 1 and take theta from Eq. 17 at psi = 90 degrees, or move C_J inward).
### 6.4 Table 1 and Eq. 11
See section 2. Table 1 digits check; Eq. 11's 2.8284 is the epsilon = 0 value.

## 7. Section V (homoclinic type cyclers, pp.11-17): content and Table 3
Not needed for P2; recorded for completeness. Class: connections of the unstable L1 Lyapunov orbit, "energy refers to the Hamiltonian, that is 0.5 C_J". Method (p.13): continuation of the system
H(x) - h = 0; g_1(x) = 0; phi_T(x) - x = 0; ||v^u||^2 - 1 = 0; D phi_T(x) v^u - Lambda^u v^u = 0; ||v^s||^2 - 1 = 0; D phi_T(x) v^s - Lambda^s v^s = 0; g_2(phi_T^u(psi^u(theta^u, xi_0))) = 0; g_2(phi_T^s(psi^s(theta^s, xi_0))) = 0; phi_T^u(psi^u(theta^u, xi_0)) - phi_T^s(psi^s(theta^s, xi_0)) = 0
(unknowns h, T, x, Lambda^u, v^u, Lambda^s, v^s, theta^u, T^u, theta^s, T^s; xi_0 small, e.g. 1e-6; method of Barrabes, Mondelo & Olle, "in preparation" 2008, ref 26). Fig. 7-8: L1 Lyapunov orbit of energy -1.5921, four connections He_1..He_4 (He_1, He_3 self-symmetric about y = 0; He_2, He_4 mirror images); Fig. 9: He_1 family at energies -1.5653, -1.5511, -1.5305, -1.5005.
Quote (p.14): "This makes these connections interesting as cycler trajectories. An spacecraft can be placed in the Lyapunov orbit, where it will be regularly flying by the Moon. When desired, the homoclinic connection can be taken, and this will give five Earth flyby opportunities before going back to the Lyapunov p.o." and "the Lyapunov p.o. needs station-keeping since it is unstable."
Printed text: Lyapunov period 29.1640 days; connection flight time periselene 1 to 19 is 113.6319 days; periselene 1 distance 6325 km (Moon two-body: aposelene of an ellipse with a = 3320 km, e = 0.905); Earth perigees (labels 8 and 12) elliptic a = 201832 km, e = 0.66; LEO rendezvous from circular radius 67808 km needs 703 m/s.
### 2008 Table 3 (p.18), "Flight times and orbital elements of the pericenters and apocenters of the homoclinic connection of energy -1.4502" (caption; the table header says energy -1.450162), transcribed from image and text layer (agree)
Units: T_flight days from periselene 1; r, a in km; v in km/s; omega in degrees counterclockwise from +x. Moon-centred rows (Moon symbol) or Earth-centred rows (Earth symbol).
p.o. of energy -1.450162 (Moon-relative):
| lbl | T_flight | r | v | a | e | omega |
|---|---|---|---|---|---|---|
| 1 | 0.000 | 6325.459 | 0.271 | 3320.262 | 0.90511 | -0.000 |
| 2 | 8.983 | 259919.372 | 0.952 | -5649.305 | 2.11119 | 70.184 |
| 3 | 14.582 | 138300.201 | 1.466 | -2359.010 | 59.62637 | 0.000 |
| 4 | 20.181 | 259919.372 | 0.952 | -5649.305 | 2.11119 | -70.184 |
connection of energy -1.450162:
| lbl | body | T_flight | r | v | a | e | omega |
|---|---|---|---|---|---|---|---|
| 1 | Moon | 0.000 | 6573.556 | 0.249 | 3429.393 | 0.91708 | -0.760 |
| 2 | Moon | 8.565 | 245836.046 | 0.952 | -5658.111 | 1.00021 | 70.381 |
| 3 | Moon | 13.730 | 142540.921 | 1.476 | -2323.020 | 62.31635 | 3.139 |
| 4 | Moon | 19.935 | 300136.002 | 0.940 | -5763.960 | 6.04401 | -67.513 |
| 5 | Moon | 29.939 | 49872.612 | 0.820 | -10300.291 | 1.07234 | -54.055 |
| 6 | Earth | 34.250 | 70971.632 | 3.063 | 215364.185 | 0.67046 | -62.458 |
| 7 | Earth | 39.973 | 357925.556 | 0.609 | 214728.164 | 0.66695 | 41.848 |
| 8 | Earth | 45.666 | 67808.075 | 3.128 | 201832.084 | 0.66404 | 147.477 |
| 9 | Earth | 51.236 | 353798.282 | 0.602 | 210749.263 | 0.67892 | -106.097 |
| 10 | Earth | 56.816 | 68883.803 | 3.122 | 218603.529 | 0.68489 | -0.000 |
| 11 | Earth | 62.396 | 353798.284 | 0.602 | 210749.266 | 0.67892 | 106.097 |
| 12 | Earth | 67.965 | 67808.081 | 3.128 | 201832.090 | 0.66404 | -147.477 |
| 13 | Earth | 73.659 | 357925.559 | 0.609 | 214728.168 | 0.66695 | -41.848 |
| 14 | Earth | 79.382 | 70971.637 | 3.063 | 215364.190 | 0.67046 | 62.458 |
| 15 | Moon | 83.693 | 49872.615 | 0.820 | -10300.290 | 1.07234 | 54.055 |
| 16 | Moon | 93.697 | 300136.005 | 0.940 | -5763.960 | 6.04401 | 67.513 |
| 17 | Moon | 99.902 | 142540.921 | 1.476 | -2323.020 | 62.31635 | -3.139 |
| 18 | Moon | 105.067 | 245836.045 | 0.952 | -5658.111 | 1.00021 | -70.381 |
| 19 | Moon | 113.632 | 6573.556 | 0.249 | 3429.393 | 0.91708 | 0.760 |
NOTE: this "Table 3" of the 2008 paper is the homoclinic table; the 2010 paper's Table 3 is the nine-plus-seven-row Earth-Moon cycler table (the one the catalogue rows come from). Do not confuse them.

## 8. What fails and why (as stated, plus consequences)
- Two-body seed corrected at the physical mass fails for a below about 0.6 (apogee near the Moon) and gives no control of the longitude of perigee (pp.7, 19).
- mu continuation at fixed T drives the orbit onto the Moon (pp.9-10) because the family does not exist uniformly in mu; 2010 adds that this is "especially for resonant orbits with high p/q".
- The p = q = 1 case has no second-species seed (2003 digest, Henon 1997).
- Small mu needs high-accuracy integration and very small continuation steps (p.9).
- Only a few C_J in the seed interval give genuinely periodic orbits at mu = 1e-6 (p.9), because of the close flyby; the corrector must be run on a grid of C_J values, not continued along a one-parameter family at fixed mu = 1e-6.
- The seed's resonance is not preserved by the continuation (2010 p.1629), so final resonance membership is decided at mu_M by a C_J walk and an intersection with the line p/q.
- The 2008 paper contains no equation giving the corrector and no stopping rule for the continuation; Table 2 gives no Moon passage data (periselene, perigee) and the "application" column is absent.

## 9. Recipe for the P2 build (task #899 / P2 in `2026-10-04-897-technique-synthesis-papers-to-problems.md`)
Goal: reproduce the catalogued Casoliva et al. 2010 rows (ids `casoliva-1-2c-em-resonant-po-2010`, `casoliva-1-2d-...`, `casoliva-1-2e-...`, `casoliva-2-1a-...`, `casoliva-2-1b-...`, `casoliva-3-2c-...`, `casoliva-7-3a-em-cycler-2010`, `casoliva-7-3b-...`, `casoliva-7-3c-...`) from seeds at mu = 1e-6, continued to mu_M = 0.0121529529 ([2010] p.1629), without using the printed Table 3 states in any start file.
Labelling for the whole build: designation p-q means p spacecraft revolutions per q Moon revolutions; seed Kepler a = (q/p)^(2/3); T about 2 pi q. (Code-level check: with (p, q) = (1, 2) the C_J interval must be (-1.711013183, 2.970934233).)
Step 1, interval and grid. For each (p, q) coprime with p/q <= 2 sqrt 2: C_J1, C_J2 from Eq. 15; a grid of 500 values of C_J strictly inside it ([2010] p.1628); plus a refined scan within 1e-2 of C_J1 (section 6.3) and within the same distance of C_J2.
Step 2, two seeds per C_J. psi = asin(s) and pi - asin(s), s = (2 - C_J + (p/q)^(2/3))/(2 sqrt(3 - C_J)) (Eq. 16); clip s to [-1, 1] only inside the endpoint tolerance. For p = q = 1 stop (no seed).
Step 3, theta. Solve Eq. 17 for theta in [0, 2 pi): two roots differing by pi; keep the one with cos(theta - psi) > 0 (2003: exactly one). Closed form for Eq. 17 not given; root-find on a 3600-point scan with bisection (COMPUTED, works: section 6.3 values).
Step 4, state. Choose alpha in (1/3, 1/2) (the only printed example is 0.4); x = mu - 1 + mu^alpha cos(theta), y = mu^alpha sin(theta); speed from the EXACT Jacobi relation v^2 = 2 Omega(x, y) - C_J at mu = 1e-6 (not from sqrt(3 - C_J)); u = v cos(psi), v_y = v sin(psi). Guess period T = 2 pi q (Eq. 18; the O(mu^alpha) correction is not given as a coefficient, do not invent one).
Step 5, correction. Differentially correct the periodic-orbit equations (the unknowns of the section II.C system with mu fixed = 1e-6: x_0 on a Poincare section, T) at fixed C_J until the closure error ||X(T) - X(0)|| is below 1e-10 ([2010] p.1628), with an adaptive integrator at tolerance below 1e-13 (RKF 7-8 in the paper), regularised about the Moon if the project's `core/cr3bp_regularized.py` is used. "Only a few initial conditions per p-q pair ... yield periodic orbits" (p.9): record the success rate; expect most grid points not to close.
Step 6, minimal period and phase checks. Closed orbit must have minimal period near 2 pi q (not a multiple), a passage of the Moon below about 1e-1 (second-species regime), and the Earth-centred a within a few percent of (q/p)^(2/3) (Table 2: within 6e-3).
Step 7, three-step continuation to mu_M (2010): leg 1, continue in mu at fixed T (INFERRED fixed T first) until the periselene distance falls below the trigger (lunar impact threatens; the 2010 text gives no number: choose a trigger of about 1e-2 lunar distances and record it, REGISTER AS A DECISION); leg 2, continue in C_J at fixed mu to raise the periselene; leg 3, resume mu continuation at fixed C_J to mu = mu_M. Use pseudo-arclength (`search/mu_continuation.py`) with a tangent-cosine guard; step size must be very small at small mu.
Step 8, family walk at mu_M. From each orbit reached, continue in C_J in both directions at mu = mu_M until natural termination (dimension, C_J or T unbounded); record the family segment (the 2010 Fig. 3 plots a quantity labelled "Resonance relation p/q" against C_J, with the seed's p/q as a thin horizontal line; how that ratio is computed from an orbit is not stated, INFERRED: the ratio of the Moon's period to the Earth-centred Keplerian period of the orbit, whose natural form is 2 pi divided by the sidereal period; choose and record a definition, and require that it equals p/q = 1/2, 2, 3/2, 7/3 on the printed Table 3 rows); find intersections with the seed's p/q. Orbits at the intersections are the candidate rows; compute k and compare.
Step 9, controls. Start control, 2008 Table 2: build the seed from its C_J and (p, q) (the ten printed C_J at mu = 1e-6), correct to closure, and compare the corrected x_i, u_i, v_i, T, k at the y = 0 crossing and the Jacobi constant. Because the corrected orbits and the table are both converged solutions of the same system, they must agree to the tolerance (not to O(mu^(1-alpha))). A symmetric orbit is found only if a seed converges to it; compare modulo the time-reversal and mirror symmetries and the choice of crossing (the table lists one y = 0 crossing; orbits have several).
End control, 2010 Table 3 at mu_M (printed C_J, period, k for the catalogued rows; these values are in `search/earth_moon_resonant_families.py::TABLE3_ROWS`, transcribed there, and I compared the nine rows' C_J with the 2010 page image: 1-2c 1.5691874798, 1-2d 2.5803060666, 1-2e 2.7629814961, 2-1a 0.4887353098, 2-1b 1.1964188553, 3-2c 0.7089330385, 7-3a 1.0215696153, 7-3b = 7-3c 1.0687623900, all period 2 pi q or 18.8495559215 = 6 pi for the 7-3 rows (6 pi = 18.849555921538759, COMPUTED)). Compare C_J, T, k and the stability class of each reached orbit with the row, and recount segments for 7-3: "four family segments ... only three intersections" ([2010] Fig. 3d). All nine admitted C_J lie inside the seed intervals above (COMPUTED).
Do NOT map designations between papers: 2008 letters (12a, 21a, 23a/b, 32a/b, 52a, 54a/b, 73a) are ordered by C_J at mu = 1e-6; 2010 letters (1-2a to 1-2e etc.) are ordered by C_J at mu_M within the continued set. 12a is not 1-2a; 73a at mu = 1e-6 has C_J = 0.8958 while 7-3a at mu_M has C_J = 1.0216.
Honest scope of the small-mass control: 2008 Table 2 holds one seed each for 1-2 (12a), 2-1 (21a) and 7-3 (73a) and two for 3-2 (32a, 32b), 2-3 (two), 5-2, 5-4 (two). The 2010 Fig. 3 shows 5, 4, 4 and 4 family segments for 1-2, 2-1, 3-2 and 7-3 respectively (letters a to e, a to d, a to d, a to d), so the 2010 scan used more seeds than 2008 lists. Therefore "at least seven of the nine clean rows" (P2 pre-registered pass) needs the FULL 500-point scan of step 1, not the ten 2008 seeds. The ten seeds alone can at most test the start-side machinery.
Note on the 2010 Fig. 3 segment counts: the letters in the legends are 12a to 12e (five), 21a to 21d, 32a to 32d, 73a to 73d, read from the page image.

## 10. How it applies to the project's problems
- It is the published one-moon generator of cyclers by seed plus continuation in mass (the project never built it; `#780` used only the two-body route the authors report fails below a = 0.6).
- The seed side is only a starting guess: the controls are about the corrected orbits and the continuation, not about the seed formulas' exact numbers.
- Out of reach of this method as published: more than one flyby per period, the Ross & Roberts-Tsoukkas families above C_J = 3 (the seed needs C_J < 3), and two moons (a collision arc about one moon only).

## 11. Proposed corpus-index row
| casoliva-mondelo-villac-mease-barrabes-olle-2008-families-cycler-trajectories-earth-moon-AIAA-2008-6434.pdf | 2026-10-04-digest-casoliva-2008-aiaa-families-cycler-trajectories-seeds.md | AIAA 2008-6434: Table 2 = ten corrected second-species seeds at mu = 1e-6 (12a to 73a; reproduce the printed C_J to 1e-16 at mu = 1e-6), Eq. 14-18 seed formulas, homoclinic He1 table; p-q labelling in Eq. 9/Table 1 is opposite to Table 2 and Eq. 15-16 (Barrabes & Gomez labelling); recipe for the P2 build; text layer (10021 words); full digest with checks. | DIGESTED (#920) |
