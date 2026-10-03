# Digest - Kumar (2025), "Analysis of Unstable Resonant Orbits for Saturn Tour Design: Between Titan and Rhea" (IAC-25-C1.9.6)

**Digested:** 2026-10-03 (text-layer PDF, 15 pages; read in full by text extraction, and every
page carrying a figure or table was also read as a rendered image). **Citation:** Bhanu Kumar
(Department of Mathematics, University of Michigan), paper IAC-25-C1.9.6, 76th International
Astronautical Congress (IAC), Sydney, Australia, 29 Sep - 3 Oct 2025, 15 pages. Copyright
International Astronautical Federation. The CrossRef DOI 10.52202/083087-0076 is carried in the
file name only; the paper's own text does not print a DOI for itself. Filed in the private paper
corpus as
`kumar-2025-analysis-unstable-resonant-orbits-saturn-tour-design-titan-rhea-IAC-25-C1.9.6-doi-10.52202-083087-0076.pdf`
(md5 154852a05593afd3a6eb275640e5cebc). All page numbers below are the printed "Page N of 15".

**Relation to other Kumar items already in the corpus.** The paper builds on Kumar, Anderson and
de la Llave (secondary resonance overlap, AAS 2023, its [13]; whiskered tori, CMDA 2022, its [20]),
Kumar's multi-shooting parameterization paper (its [12], 2025) and an in-preparation
multiple-shooting paper (its [24]). Those have their own digests; this note covers only what this
paper prints.

## What it is
A conference paper in two parts. Part 1 (Sections 4 and 5) studies, in two separate planar
circular restricted three-body problems (Saturn-Titan, Saturn-Rhea), unstable resonant periodic
orbits, their stable and unstable manifolds on Poincare sections, and zero-delta-v heteroclinic
connections between resonances at stated Jacobi constants, for the "pump-down" region between
Titan and Rhea. Part 2 (Section 6) takes one family, Titan 3:2 low/medium-eccentricity, and
continues it into the concentric circular restricted four-body problem (CCR4BP: Saturn, Titan as
base moon, Rhea as perturber) to look for invariant tori and Rhea-induced secondary resonances.
Abstract: "we study resonances in the region between Titan and Rhea. After first studying
single-moon pump-down tour segments using restricted 3-body models for each moon separately, we
then study the effect of Rhea on some Titan resonant orbits using a restricted 4-body model."

## (1) Models and constants exactly as printed
- **PCRTBP** (Section 2.1, p.2): Hamiltonian H0 = (px^2+py^2)/2 + px*y - py*x - (1-mu)/r1 - mu/r2,
  synodic frame centred on the m1-m2 barycentre, units with m1-m2 distance 1, period 2*pi,
  G(m1+m2) = 1. Jacobi constant C = -2H. "C is also approximately (within order O(mu)) equal to
  the well-known Tisserand parameter T, which for higher energies (C < 3) is related [15] to the
  traditional patched-conic m2-relative flyby v_inf as C ~ T = 3 - v_inf^2" (in normalised units).
- **Mass parameters printed:** Saturn-Titan mu = 2.36639e-4; Saturn-Rhea mu = 4.05841e-6 (p.2).
  In the CCR4BP, "for m1, m2, and m3 being Saturn, Titan, and Rhea, mu3 = 4.05746e-6 while mu
  remains 2.36639e-4" (p.2). mu3 is defined in the text as the m3 mass over (m1+m2), the
  typeset fraction is garbled in extraction but the printed numbers are as quoted. The Rhea
  value 4.05841e-6 (Saturn-Rhea two-body) and 4.05746e-6 (as the perturber in the Titan
  frame) differ in the fourth significant figure; the paper does not comment.
- **CCR4BP** (Section 2.2, p.2-3): three masses m1 >> m2, m3; m2 and m3 move on coplanar
  concentric circles about m1 with radii r12 and r13; "m2 has no effect on the motion of m3 nor
  vice versa". Angular velocities by Kepler's third law, Omega_i = sqrt(G(m1+m_i) r1i^-3). The
  justification printed is "those moons' near-zero inclinations (~0.3 deg) and eccentricities
  (~0.03 and 0.001)" (Titan, Rhea). Normalisation: G(m1+m2), r12 and Omega2 all 1 (Titan
  units). Rhea's phase angle theta3(t) = (Omega3 - 1) t + theta3,0; m3 position
  (-mu + r13 cos(theta3), r13 sin(theta3)). Equations of motion Eq. [4] in
  position-momentum form (x, y, px, py) (derivation attributed to Blazevski and Ocampo [16]);
  time-periodic Hamiltonian Eq. [5].
- **Base moon / perturber:** in the four-body part Titan is the base moon (m2) and Rhea the
  perturber (m3). The roles are not swapped anywhere. The reverse (Titan perturbing Rhea) is
  listed as future work (p.14).
- **Orbit radius of the perturber:** "Rhea... has an orbital radius of about 0.4315 normalized
  units" (Section 6.1, p.12; Titan units). Forcing quantities: "Tp = 2.48376 (in Titan-normalized
  time units) being the Titan-relative synodic period of Rhea" (p.12), with Tp = 2*pi/|Omega3 - 1|
  (p.3). The paper does not print an explicit period ratio of Titan to Rhea; from the printed
  r13 = 0.4315 and Kepler's third law one could compute it, but the paper does not.
- **Rhea collision threshold used:** closest approach to the centre of Rhea above 0.00145
  normalised units (Saturn-Rhea units, p.9). The Titan collision threshold in normalised units
  is not stated.
- **Units:** dimensionless synodic units throughout. Time-of-flight in Figures 4, 6, 9, 10 is
  labelled "TU" (not defined in the text).

## (2) Three-body content (all PCRTBP, Sections 4 and 5)
**Saturn-Titan (Section 4).** Resonances studied: 1:2 (exterior), 1:1, 3:2 (interior); "an m:n
resonant orbit is called exterior here if it has semi-major axis larger than that of Titan
(m < n)". Each of 1:2 and 3:2 has two unstable families, a low/medium-e and a high-e family
(Figures 1, 2). Bifurcation diagrams (Fig. 2, symmetric x-intercept vs C), axis ranges read by
eye: 1:2 panel C about 2.85 to 3.15, x-intercept about -1.6 to -2.3; 3:2 panel C about 2.9 to 3.06,
x-intercept about -0.45 to -0.75. The 1:1 family is the planar L1 Lyapunov family (Fig. 3); "the
Jacobi constant decreases monotonically as the orbit size increases"; a member at C ~ 2.968 is
shown in the inertial frame in Fig. 4 (figure title reads C = 2.9678 by eye). The high-e
families are stable at very high e and low C, become unstable as e decreases, and have a fold
bifurcation where C takes a maximum (p.6). No periodic-orbit initial conditions, periods or
multipliers are printed for Titan.
- Manifolds: periapse section (l = 0) in synodic Delaunay coordinates, L = sqrt((1-mu)a),
  G = L sqrt(1-e^2), l = mean anomaly, g = argument of periapse in the rotating frame (p.7).
  Fig. 5 caption: "Titan 1:2 (red/blue) and 3:2 (orange/green) unstable and stable manifolds
  respectively", C = 2.957, globalized 3 revolutions each. The body text then reads: "Intersections between the orange 1:2
  unstable and blue 3:2 stable manifolds are clearly visible, indicating the presence of
  heteroclinic transfers from 1:2 to 3:2 Titan MMRs, corresponding to a decrease in semi-major
  axis" (p.8). The colours named in the text (orange 1:2 unstable, blue 3:2 stable) do not match the
  caption colour assignment; this note does not resolve which is meant. Axis ranges read by eye: g 0 to 2*pi, L about 0.8 to 1.6.
- Fig. 6 (p.8): one explicit 1:2 (high-e, exterior) to 3:2 (high-e, interior) heteroclinic
  trajectory at C = 2.957 passing through a 1:1 Lyapunov-like segment. "While we did not compute
  its manifolds, it can be deduced that this 1:2 to 3:2 transfer is in fact a heteroclinic formed
  by chaining 1:2 to 1:1 and 1:1 to 3:2 heteroclinic connections." The figure title carries an
  elapsed time of 35.499997 TU (read by eye).
- The paper also gives a patched-conic reading of these trajectories (flyby bends velocity back
  outbound), qualitatively only.

**Saturn-Rhea (Section 5).** Only exterior resonances: 1:2, 2:3, 5:7, 3:4 (p.8). Two families
per resonance (low/medium-e and high-e), Fig. 7. No Rhea-surface-colliding orbits.
**Table 1 (p.9), transcribed and checked against the rendered page:**

| m | n | Family | Min C | Max C |
|---|---|--------|-------|-------|
| 1 | 2 | low/mid-e | 2.9718 | 3.1497 |
| 1 | 2 | high-e | 2.9662 | 2.9698 |
| 2 | 3 | low/mid-e | 2.9879 | 3.0524 |
| 2 | 3 | high-e | 2.9810 | 2.9865 |
| 5 | 7 | low/mid-e | 2.9910 | 3.0364 |
| 5 | 7 | high-e | 2.9802 | 2.9899 |
| 3 | 4 | low/mid-e | 2.9931 | 3.0266 |
| 3 | 4 | high-e | 2.9861 | 2.9921 |

Text around it: "for each Rhea MMR, there are ranges of C values (e.g., (2.9865, 2.9879) for Rhea 2:3)
where no appropriate orbit exists"; all mid-e Rhea orbits found lie entirely outside Rhea's orbit
circle; "the 3:4 high-e orbits become stable for C < 2.9861 (and the mid-e ones pass through
Rhea)" (p.11).

**Rhea manifolds and heteroclinics (apoapse section, (L, g) coordinates, Fig. 8, p.10).**
Fig. 8 plots 2:3 (red/blue), 5:7 (magenta/cyan) and 3:4 (orange/green) unstable/stable
manifolds at five Jacobi constants. The text gives these globalization levels (Poincare
mappings, "revs"):

| C | globalization | outcome as printed |
|---|---------------|--------------------|
| 3.00 | over 25 | manifolds "well-separated"; homoclinic tangles visible for 3:4 and 5:7; no heteroclinics; same for all C > 3.00 |
| 2.99 | 15 | stronger transport to/from 3:4; 2:3 isolated from 3:4; no 5:7 orbit exists at this C (Table 1) |
| 2.989 | 10 | 5:7 reappears; 3:4 and 5:7 manifolds come near 2:3 but do not intersect; 2:3 to 5:7 connections need order-20 resonances, over 20 revs |
| 2.988 | 5 | "clear intersections of 2:3 unstable and 5:7 stable, and 5:7 unstable and 3:4 stable manifolds"; implies 2:3 to 3:4 via 5:7; 2:3 orbit is high-energy mid-e family near the impact threshold 2.9879, 5:7 and 3:4 orbits are high-e |
| 2.9863 | 5 | direct intersections of 2:3 and 3:4 manifolds |

Also: heteroclinic intersections "do not exist for C = 2.9885" (studied the same way).
Computed trajectories (Section 3.1.3 method, Eq. [8]):
- C = 2.988, 2:3 to 3:4 via 5:7: passes through Rhea's surface (Fig. 9 top row), "this transfer is
  not possible in reality".
- C = 2.988, mid-e 2:3 to high-e 3:4 via a 7:10 segment ("identifiable by its 7 loops
  corresponding to periapses"): no Rhea impact; "requiring 3 extra revs total" relative to the
  5:7 route (Fig. 9 bottom row; figure titles read by eye give elapsed time 113.097336 TU for both).
- C = 2.9863, 2:3 to 3:4 direct: collides with Rhea (Fig. 10 top); five revolutions
  shadowing a 5:7 high-e orbit: no collision (Fig. 10 bottom; titles read by eye 50.265482 TU and
  about 94.24778 TU).
- Summary printed (p.11): "the key Jacobi constant values for transfers from 2:3 to 3:4 orbits
  are in the vicinity of C = 2.988 for mid-e 2:3 to high-e 3:4 transfers via the 7:10 MMR, and
  between 2.9861 < C < 2.9865 for transfers via 5:7." Only the range 2.9861 <= C <= 2.9865
  remains for the 5:7 route. A fine grid of C was not examined. 3:4 manifolds are said to show
  strong transport to lower semi-major axis, so 3:4 to 1:1 heteroclinics "likely" exist (not
  computed).

## (3) Four-body content (Section 6, pp.11-14), in full
**What is computed.** Only one family: Titan 3:2 low/medium-e PCRTBP periodic orbits, continued
in the perturber's mass mu3 from 0 (PCRTBP) to "realistic mu3 for Rhea" in the Saturn-Titan-Rhea
CCR4BP, using the stroboscopic map (Rhea's synodic period Tp = 2.48376 Titan time units). "Unless
otherwise specified, in this section all plots will display objects invariant under F rather
than under the flow." Orbit-family plots use adapted coordinates (l + g, -3L + 2G), in which "for
mu3 = 0 trajectories will lie largely on horizontal lines".
**Crossing the perturber's orbit.** "None of the 3:2 Titan resonant PCRTBP periodic orbits
displayed in Figure 1, nor the one shadowed in Figure 6, intersect that of Rhea... However, some
of the 3:2 orbits at higher and moderate eccentricities do approach Rhea's orbit fairly closely
at periapse" (p.12). So no Rhea-orbit-crossing orbit is treated; the paper states that whether
"PCRTBP-like" quasiperiodic analogues exist at moderate or high e "is not immediately obvious".
**Rotation number and tori.** Eq. [9] F_mu3(K(theta)) = K(theta + omega), rotation number
omega = 2*pi*Tp/T. Torus algorithm: quasi-Newton method of Kumar, Anderson and de la Llave [20],
solving for K together with stable/unstable/central directions and Floquet matrix. During mass
continuation "the distance r13 must be varied to keep the third body's frequency Omega3 and thus
omega constant". The paper warns that if omega/(2*pi) = Tp/T is rational the orbit is at a
secondary resonance and periodic orbits (not tori) should be computed, and that "if neighboring
secondary resonance regions overlap and destroy the tori between them, this continuation will
fail before reaching the desired final mu3 value."
**Secondary resonances found.** For the 3:2 low/mid-e orbits "the PCRTBP, Tp/T ranges from
0.1898475 to 0.1952735". A Farey sequence with q < 50 lists ratios in that range "including
4/21, 5/26, 6/31, 7/36, 8/41, and 9/47". Fig. 11 plots PCRTBP period T against Jacobi constant
(C axis about 2.94 to 3.06, T axis about 12.7 to 13.1, read by eye) with red dashed lines at the
secondary-resonance periods, labelled by the inverse ratio T/Tp: 21/4, 47/9, 26/5, 31/6, 36/7,
41/8 (top to bottom). The T values themselves are not tabulated. (Arithmetic from the printed
Tp: 21/4 x 2.48376 is about 13.04, 41/8 x 2.48376 is about 12.73; these are this note's
arithmetic, not paper values.) The curve in Fig. 11 is non-monotonic in C (read by eye): it
falls from about 13.1 near C = 2.95, has a local minimum near C = 2.97 and a local maximum
near C = 3.01, falls again, and rises steeply near C = 3.06. Secondary resonance period is
p*T = q*Tp: "if Tp/T = p/q, then the equivalent CCR4BP orbits' periods will be pT", a q-iteration
periodic orbit of F (p.4).
**Outcome.** "We successfully numerically continued non secondary resonant PCRTBP 3:2 Titan
orbits from a wide range of Jacobi constants as quasiperiodic orbits into the CCR4BP." The F-
invariant tori are in Fig. 12 (Cartesian, title carries mu3 = 4.057e-06 and mu = 2.366e-04) and
as horizontal curves in Fig. 13. "The extent and shape of most tori is by and large similar to
that of the corresponding PCRTBP periodic orbit family, covering a wide range of eccentricities
from near-circular to near-collision with Titan." "However, many mid-e orbits having periods
near the aforementioned secondary resonances failed to continue as tori." "Secondary resonant
long periodic orbits and separatrices for such orbits were computed and plotted in the
higher-energy part of the family" (Fig. 13; the text says "bottom of Figure 13" although Fig. 13
is a single panel). Fig. 13 shows tori, secondary-resonant periodic-orbit points (x and o
markers) and "separatrices (short red/blue curves)". "Such secondary resonance regions and
separatrices are bracketed and separated by 3:2 tori that persist into this CCR4BP... the
displayed secondary resonances do not overlap - which would have destroyed the tori between them -
and are not the dominant driver of the dynamics inside this 3:2 family despite its close passes to
Rhea's orbit. Outside of a few secondary resonance regions, the dynamics are quite 'PCRTBP-like'
still, consisting of quasiperiodic orbits. Given Rhea's small size, this is not altogether
unexpected." Conclusion (p.14): "While this family was mostly found not to experience large
changes due to effects from Rhea, some of the highest-energy orbits in the system did show
non-negligible effects despite Rhea's small mass." Introduction: "finding small but
non-negligible Rhea-induced secondary resonant phenomena." No table of values, no mu3 at which
any particular torus fails, and no list of which of the six ratios gained a second pair of weak
multipliers (and hence separatrices) is printed.
**Manifolds and connections in the four-body model.** Separatrices are the only four-body
manifold objects: "1D submanifolds (of the orbit's 2D stable & unstable manifolds) associated
with these weak directions", computed by a parameterization method for F (Eq. [11]) and shown as
short curves in Fig. 13. The paper does NOT compute the stable or unstable manifolds of any
resonant invariant torus, and does NOT compute any heteroclinic or homoclinic connection in the
CCR4BP. All heteroclinic content (Fig. 5, 6, 8, 9, 10) is PCRTBP only. The sentence that frames
the whole paper's four-body scope is the last paragraph of the Introduction: "we start with some
Titan 3:2 resonant PCRTBP periodic orbits from part 1 of the study and attempt to transition
them to the CCR4BP, finding small but non-negligible Rhea-induced secondary resonant
phenomena." The Titan 1:2 and 1:1 families, and the Titan high-e 3:2 family, are not continued
to the four-body model; this is listed as future work ("including 1:1 and 3:2 high-e", p.14).

## (4) Methods as stated
- **Stroboscopic map** (Section 2.2.1, p.3): F is the time-2*pi/|Omega3 - 1| map of the CCR4BP
  flow on the extended phase space (x, y, px, py, theta3); theta3 advances by exactly +/-2*pi so it
  can be fixed (theta3,f = 0), giving a 4D map. "The stroboscopic map preserves all the dynamical
  properties of the CCR4BP, including invariant sets and manifolds". Valid also for mu3 = 0.
- **PCRTBP periodic orbits** (3.1.1): start from the Kepler problem orbit of semi-major axis
  (n/m)^(2/3), argument of periapse and true anomaly both 0 (exterior) or both pi (interior);
  continue from mu = 0 to the target mu with the mirror theorem (Roy and Ovenden [21]); the same
  symmetric-orbit continuation, detailed in Kumar and Moreno [22], then follows the family using
  x or ydot at the perpendicular x-axis crossing as parameter.
- **Manifolds** (3.1.2): parameterization method [11] as developed by Kumar [12]; functional
  equation Phi_tau(k)(W(k, s)) = W(k+1 mod m, lambda*s), W as m Taylor series in s with W0 the
  section crossing points and W1 the scaled monodromy eigenvectors, lambda the m-th root of the
  multiplier; local curves integrated to the section, then globalized by iterating the Poincare
  map. Section: apoapse (true anomaly pi) for exterior MMRs, periapse (true anomaly 0) for
  interior MMRs.
- **Heteroclinics** (3.1.3): solve W^u_1p(k1, s1) = W^s_2p(k2, s2), Eq. [8], by discretised search
  then bisection. "If a heteroclinic intersection is found after propagating the local unstable
  and stable manifolds of two orbits by M map iterations each, the resulting heteroclinic transfer
  will require approximately 2M orbital revolutions" (3.1.4, TOF proxy).
- **CCR4BP** (3.2): tori by quasi-Newton on F(K(theta)) = K(theta + omega) from [20];
  secondary-resonant periodic orbits by a new multiple-shooting method adapted from the torus
  algorithm, F(X(k)) = X(k+p mod q), "solves for the periodic orbit simultaneously with its
  various Floquet directions and matrix", details deferred to [24] (in preparation, 2024);
  continuation in mu3 by that method; separatrices by a parameterization method analogous to the
  PCRTBP case (Eq. [11]). Chirikov's overlap criterion [18] is cited for the role of resonance
  overlap in transport.
- **Tolerances, integrator, step sizes, Newton convergence criteria, number of shooting
  segments:** not stated.

## (5) Reproduction targets actually printed
The paper prints no initial conditions, no state vectors, no periods or Floquet multipliers for
any named orbit, no rotation numbers other than the Tp/T range, and no mass-continuation
table. What is printed:
1. Constants: mu(Saturn-Titan) = 2.36639e-4; mu(Saturn-Rhea) = 4.05841e-6; mu3(Rhea in Titan
   frame) = 4.05746e-6; Rhea orbital radius 0.4315 (Titan units); Tp = 2.48376 (Titan units);
   Rhea collision floor 0.00145 (Saturn-Rhea units). (pp.2, 9, 12)
2. Table 1: Jacobi-constant extent of non-colliding unstable Rhea orbits per family (p.9).
3. Jacobi constants at which Rhea connections exist or fail: none for C >= 3.00 up to 25 revs;
   C = 2.99 isolated 2:3; C = 2.989 no 2:3 to 5:7 within 10 revs; C = 2.988 (5 revs) 2:3 to 5:7 to
   3:4 intersections, with a Rhea-colliding 5:7 route and a non-colliding 7:10 route (3 extra
   revs); C = 2.9885 no intersections; C = 2.9863 (5 revs) direct 2:3 to 3:4 (collides) and a
   non-colliding five-revolution 5:7-shadowing route; window 2.9861 <= C <= 2.9865 for the 5:7
   route; 3:4 high-e stable below C = 2.9861; no 2:3 non-colliding orbit in (2.9865, 2.9879). (pp.9-11)
4. Titan manifold demonstration at C = 2.957 (3 revs) with a 1:2 to 3:2 heteroclinic via 1:1 (pp.7-8).
5. Titan 3:2 low/mid-e: Tp/T range [0.1898475, 0.1952735]; Farey ratios Tp/T with q < 50 including
   4/21, 5/26, 6/31, 7/36, 8/41, 9/47; the corresponding secondary-resonant period is p*T = q*Tp
   (p.13).
6. Figure-only values (read by eye; not printed as numbers): Fig. 2 bifurcation curves; Fig. 11
   T(C) curve; Figs. 12-13 tori; elapsed times in figure titles (Section 2 and 3 notes above).
7. Qualitative result usable as a pass/fail: Titan 3:2 low/mid-e tori persist at the real mu3 for
   most Jacobi constants, fail near the listed secondary resonances, and the secondary-resonance
   regions do not overlap.

## (6) What remains to be done, and other moons and planets
Printed (p.14): "complete the study of Rhea PCRTBP resonance-to-resonance heteroclinics, by
studying the pump-down from 1:2 to 2:3, and the pump-down from 3:4 to 1:1. Second is the matter
of finishing the characterization of Rhea's effects on Titan resonant orbits, including 1:1 and
3:2 high-e. And third is to study the effect of Titan on Rhea resonant orbits - effects that may
be much more significant than those of Rhea on Titan's orbits, given Titan's much more massive
size. Heteroclinics from Titan resonant to Rhea resonant orbits would be a natural end goal of
this study as well. Finally, similar investigations also remain to be carried out for the other
moons in the Saturnian system relevant to Enceladus mission design - Dione, Tethys, and Enceladus
itself." Intended use: tour design toward Enceladus (decadal survey and Takubo et al [8],
Brandenburg et al [9]).
**Uranus.** One mention in the body: "...as the orbits used in those studies are not of the very
high energies where traditional patched-conic methods work best (e.g. the Uranian tour of [10])"
(p.1). Reference [10] is Landau, Davis and Karimi, "Trajectory options for a uranus orbiter and
probe", AAS/AIAA Astrodynamics Specialist Conference, 2023. No Uranian moon, no Uranian
computation. No other planet is treated; Jupiter appears only in references [1]-[6] and [13]
(Jupiter-Ganymede 4:3 plus Europa case study).

## (7) Text search results (whole paper, lines joined, ligatures normalised, case-insensitive)
cycler 0; homoclinic 2 (one in body, one in the title of reference [4]); heteroclinic 45;
"secondary resonan" 43; torus 3; tori (whole word) 23; stroboscopic 9; Hyperion 0; Dione 1
(the future-work sentence above); Enceladus 7; Uranus 1 (title of reference [10]) plus
Uranian 1 (the body sentence above); Oberon 0; Titania 0; Umbriel 0.
- homoclinic, body: "While homoclinic tangle phenomena are visible for the 3:4 and 5:7 orbit
  manifolds, no heteroclinics occur even after such a long TOF." (p.9, Rhea PCRTBP, C = 3.00).
  The other hit is the title of Anderson and Lo [4], "Flyby design using heteroclinic and
  homoclinic connections of unstable resonant orbits".
- Hyperion: no occurrence.
- Uranus/Uranian: quoted in (6).

## (8) Relevant bibliography (full, as printed)
[1] R. L. Anderson and M. W. Lo, "Role of invariant manifolds in low-thrust trajectory design,"
J. Guidance, Control, and Dynamics, 32(6), 1921-1930, Nov. 2009.
[2] R. L. Anderson and M. W. Lo, "Dynamical systems analysis of planetary flybys and approach:
Planar europa orbiter," JGCD 33(6), 1899-1912, Nov. 2010.
[3] R. L. Anderson and M. W. Lo, "A dynamical systems analysis of planetary flybys and approach:
Ballistic case," J. Astronautical Sciences 58(2), 167-194, Apr. 2011.
[4] R. L. Anderson and M. W. Lo, "Flyby design using heteroclinic and homoclinic connections of
unstable resonant orbits," Spaceflight Mechanics, Proc. 21st AAS/AIAA Space Flight Mechanics
Meeting, New Orleans, Feb. 13-17, 2011, Advances in the Astronautical Sciences vol. 140, Univelt,
2011, pp. 321-340.
[5] R. L. Anderson and M. W. Lo, "Spatial approaches to moons from resonance relative to
invariant manifolds," Acta Astronautica 105, 355-372, 2014.
[6] R. L. Anderson, "Approaching moons from resonance via invariant manifolds," JGCD 38(6),
1097-1109, Jun. 2015.
[7] National Academies, Origins, Worlds, and Life: A Decadal Strategy for Planetary Science and
Astrobiology 2023-2032, National Academies Press, 2022, DOI 10.17226/26522.
[8] Y. Takubo, D. Landau, B. Anderson, "Automated tour design in the saturnian system," Celestial
Mechanics and Dynamical Astronomy 136(1), 8, 2024, DOI 10.1007/s10569-023-10179-8.
[9] W. Brandenburg, R. P. Russell, M. Shaw, "Pathfinding for pump-down flyby tours at Saturn
using directed graphs," AAS/AIAA Astrodynamics Specialist Conference, Broomfield CO, Aug. 2024.
[10] D. Landau, A. Davis, R. Karimi, "Trajectory options for a uranus orbiter and probe,"
AAS/AIAA Astrodynamics Specialist Conference, Big Sky MT, Aug. 2023.
[11] A. Haro, M. Canadell, J. Figueras, A. Luque, J. Mondelo, The Parameterization Method for
Invariant Manifolds, Springer, 2016, vol. 195.
[12] B. Kumar, Multi-shooting parameterization methods for invariant manifolds and heteroclinics
of 2 DOF Hamiltonian Poincare maps, with applications to celestial resonant dynamics, 2025.
[13] B. Kumar, R. L. Anderson, R. de la Llave, "4th body-induced secondary resonance overlapping
inside unstable resonant orbit families: A Jupiter-Ganymede 4:3 + Europa case study," AAS/AIAA
Astrodynamics Specialist Conference, Aug. 2023.
[14] A. Celletti, Stability and Chaos in Celestial Mechanics, Springer, 2010,
DOI 10.1007/978-3-540-85146-2.
[15] S. Campagnola, R. P. Russell, "Endgame problem part 2: Multibody technique and the
tisserand-poincare graph," JGCD 33(2), 476-486, 2010, DOI 10.2514/1.44290.
[16] D. Blazevski, C. Ocampo, "Periodic orbits in the concentric circular restricted four-body
problem and their invariant manifolds," Physica D 241(13), 1158-1167, 2012,
DOI 10.1016/j.physd.2012.03.008.
[17] A. Morbidelli, Modern celestial mechanics: aspects of solar system dynamics, Taylor and
Francis, 2002.
[18] B. V. Chirikov, "Resonance processes in magnetic traps," Soviet J. Atomic Energy 6(6),
464-470, 1960, DOI 10.1007/BF01483352.
[19] B. Kumar, R. L. Anderson, R. de la Llave, "High-order resonant orbit manifold expansions for
mission design in the planar circular restricted 3-body problem," Commun. Nonlinear Sci. Numer.
Simul. 97, 105691, 2021, DOI 10.1016/j.cnsns.2021.105691.
[20] B. Kumar, R. L. Anderson, R. de la Llave, "Rapid and accurate methods for computing whiskered
tori and their manifolds in periodically perturbed planar circular restricted 3-body problems,"
Celestial Mechanics and Dynamical Astronomy 134(1), 3, 2022, DOI 10.1007/s10569-021-10057-1.
[21] A. E. Roy, M. W. Ovenden, "On the occurrence of commensurable mean motions in the solar
system: The mirror theorem," MNRAS 115(3), 296-309, Jun. 1955, DOI 10.1093/mnras/115.3.296.
[22] B. Kumar, A. Moreno, "Networks of periodic orbits in the Earth-Moon system through a
regularized and symplectic lens," AAS/AIAA Astrodynamics Specialist Conference, Aug. 2025.
[23] J. S. Parker, R. L. Anderson, Low-Energy Lunar Trajectory Design, JPL Deep Space
Communications and Navigation Series, Wiley, Jun. 2014, vol. 12.
[24] B. Kumar, "A new fast multiple-shooting method for computing periodic orbits in symplectic
maps leveraging simultaneous floquet vector computation to avoid large linear systems," in
preparation, 2024.

## What it does NOT contain
- No four-body heteroclinic, homoclinic or other connection, and no stable or unstable manifold of
  any invariant torus. The only four-body manifold objects are the separatrices of
  secondary-resonant periodic orbits (Fig. 13), shown as short curves, with no connection computed.
- No torus for Titan resonances other than 3:2 low/medium-e; no torus for any Rhea resonance; no
  Rhea orbit perturbed by Titan.
- No initial conditions, periods, Floquet multipliers, rotation numbers per orbit, or tabulated
  continuation results for any periodic orbit or torus; no mu3 value at which a torus is lost.
- No Hyperion, Dione, Enceladus, Tethys computation; no Uranus system computation; no cyclers.
- No tolerances, integrator, or convergence statistics.
- No computation of tours or delta-v budgets; the "tour design" framing is motivation only.

## What a reproduction attempt could use
- Constants and model definition: pp.2-3 (mu = 2.36639e-4, mu3 = 4.05746e-6, Eq. [4]-[5], Omega
  from Kepler's third law, Rhea radius 0.4315 and Tp = 2.48376 on p.12).
- Titan 3:2 low/mid-e secondary-resonance ratios 4/21, 5/26, 6/31, 7/36, 8/41, 9/47 and the Tp/T range
  0.1898475 to 0.1952735 (p.13) as a check that an independent code reproduces the same
  resonance list, then the pass/fail statement that tori persist between them at real mu3 (pp.13-14).
- Rhea Table 1 (p.9): the C extent of each non-colliding family is directly checkable by an
  independent PCRTBP periodic-orbit code (needs the Rhea collision floor 0.00145).
- Rhea heteroclinic Jacobi-constant windows (C = 2.988 via 7:10 and 5:7; 2.9861 to 2.9865 via 5:7;
  none at 2.9885; none for C >= 3.00 within 25 revs), pp.9-11.
- Titan C = 2.957 1:2 to 3:2 heteroclinic (p.8).
- Not available from this paper: any four-body connection to reproduce. The four-body
  reproducible content is limited to tori persistence and secondary-resonance structure of one
  Titan family perturbed by Rhea.
