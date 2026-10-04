# Digest: Neelakantan & Ramanan (2022), "Two-impulse transfer to multi-revolution halo orbits in the Earth-Moon elliptic restricted three body problem framework"

J. Astrophys. Astron. 43:50 (2022), DOI 10.1007/s12036-022-09830-x, 16 pages (article pp1-15, references pp15-16).
Rithwik Neelakantan and R. V. Ramanan, Department of Aerospace Engineering, IIST, Thiruvananthapuram, India. MS received
29 October 2021, accepted 18 February 2022, published online 12 August 2022. Footnote (p1): "Presented in 2020 AAS/AIAA
Astrodynamics Specialist Conference; AAS 20-607."
Filed in the private paper corpus as
neelakantan-ramanan-2022-two-impulse-transfer-multi-revolution-halo-orbits-earth-moon-ER3BP-jaa-43-50-doi-10.1007-s12036-022-09830-x.pdf

Digested 2026-10-04 from all 16 pages. Each statement is marked READ (seen on the page, with page, section, equation or
table) or INFERRED (our reading or arithmetic on printed entries). Page numbers are the journal's printed "Page n of 16".
Tables 2 to 9 are transcribed digit by digit; where a digit count in a long string of zeros is hard to see it is flagged.

## 0. What the paper is

READ (abstract, p1): a direct, single-segment, two-impulse transfer from a circular Earth parking orbit (EPO) to a
multi-revolution (MR) halo orbit about Earth-Moon L1 in the Earth-Moon elliptic restricted three-body problem (ERTBP), with no
manifolds and no bridge segment. "The location of insertion into the MR halo orbit and the components of the insertion
velocity are treated as unknowns and obtained using differential evolution". "The optimal solutions indicate that there
exist trajectories with lower cost and for significantly lower time of flight than those reported in the literature for
similar problems." The method is benchmarked first at e = 0 (planar-symmetric CRTBP halo, Az = 15,000 km) against Rausch
(2005), then applied to the MR orbit M5N2 and four further MR orbits.

## 1. The model as printed (section 2, pp3-5)

READ. Primaries move on Keplerian ellipses of semi-major axis a and eccentricity e about their barycentre; frame is "a
non-uniformly rotating and pulsating frame", origin at the barycentre, "instantaneously normalized by the distance between
the primaries" (p4, Figure 1, drawn with m1 on the left of the origin and m2 on the right at positive x). Primary separation
(eq. 1, p4): R(v) = a(1 - e^2)/(1 + e cos v), v the true anomaly. Independent variable change (eq. 2, p4):
d/dt = (dv/dt) d/dv, dv/dt = (1 + e cos v)^2/(1 - e^2)^(3/2).

Equations of motion (eqs. 3-4, p4), dots are d/dv:

    x'' - 2 y' = omega_x,   y'' + 2 x' = omega_y,   z'' = omega_z
    omega(x,y,z,v) = (1 + e cos v)^(-1) * Omega_bar(x,y,z)
    Omega_bar = (1/2)(x^2 + y^2) + (1-mu)/r1 + mu/r2 + (1/2) mu (1-mu) - (1/2) e cos(v) z^2

with r1, r2 the distances to the larger and smaller primaries. Values used (p5): "The values of mass ratio and eccentricity
for the Earth-Moon system employed in this study are mu = 0.0122 and e = 0.0554, respectively."

INFERRED: this is the same system as `core/er3bp.py` and the Shu-Lin 2025 eq. (1) (the z-equation z'' = -(e cos v z + grav)/(1 + e cos v)
equals Z'' + Z = dOmega/dZ with the full (1/2)Z^2 term). The paper does not print the primary positions as coordinates; from
the L1 halo initial conditions (x0 = 0.8229... for the e = 0 halo, below, and the L1 lying near x = 0.837 in a frame with the
Moon at 1 - mu) the frame is the Szebehely one used by the project (larger primary at x = -mu, smaller at 1 - mu). The
mass ratio 0.0122 differs from the canonical Earth-Moon 0.0121505843 used elsewhere in the project; for any reproduction the
printed 0.0122 and 0.0554 must be used literally.

## 2. What is computed (sections 3-5)

### Periodicity in the ERTBP (section 3, pp2, 5)

READ (p2): "The halo orbits in the ERTBP framework are periodic only after multiple revolutions around the Lagrangian point
and have periods which are rationally commensurable with that of the primaries (Peng & Xu 2015a). They are known as
multi-revolution (MR) orbits". The initial point on the periodic orbit "must be chosen such that the primaries are at an
apse" (p2). Periodicity condition in the planar case: Moulton et al. (1920), "two perpendicular crossings with the syzygy axis".

READ (p5): "the target MR orbits are denoted by the notation MaNb, where a and b are integers denoting the values of M and N,
respectively. The third body makes 'a' revolutions around the Lagrangian point in these orbits, whereas the primaries
complete 'b' revolutions around their barycentre. For example, in the orbit M5N2, the spacecraft makes five revolutions
around the Lagrangian point while the primaries complete two revolutions around their barycentre." Commensurability (eq. 5,
p5, after Peng and Xu 2015b): T_E = 2 N pi = M T_C, with T_E the period of the MR orbit and T_C the period of the CRTBP halo.
READ (p5): the MR orbit "cannot be characterized uniquely by the Az amplitude. The average out-of-plane amplitude is used to
represent the orbits (Rithwik & Ramanan 2019)."

Generation of an MR orbit (p5): the state at t = 0 is [x0, 0, z0, 0, ydot0, 0] and at half period [x_T/2, 0, z_T/2, 0, ydot_T/2, 0]
(orthogonal crossing of the x-z plane, same symmetry as the CRTBP, "despite being a non-autonomous system"). Unknowns
[x0, z0, ydot0] are found by a differential-evolution (DE) single-segment search (Rithwik & Ramanan 2019), with search bounds
"fixed in the neighborhood of the optimal design parameters that are obtained under the CRTBP framework". Two MR halo orbits
used here are therefore designed by DE in the earlier papers (Rithwik & Ramanan 2019, 2021), not by continuation in e. READ
(p2): that earlier work "designed MR orbits in a single-level approach based on differential evolution, avoiding the need for
numerical continuations on eccentricity and mass ratio ... generated both Lyapunov and halo solutions for the same MR orbit
and captured multiple solutions for each of the halo/Lyapunov MR orbits".

The paper describes Peng and Xu 2015a as follows (p2, READ): "Peng & Xu (2015a) constructed MR halo orbits in the Earth-Moon
system using the halo orbit conditions in CRTBP as the initial guess. They divided the MR halo orbit into multiple segments
and employed numerical continuation on eccentricity to generate the design." Nothing from that paper's tables is reproduced.

### Transfer construction (section 4, pp5-7)

READ. Backward design: choose a location on the MR orbit (true anomaly v of the Moon around the Earth at insertion) and a
halo-orbit-insertion (HOI) impulse [dxdot, dydot, dzdot]; add it to the MR-orbit state; integrate the equations (3) backward
in time with a Runge-Kutta-Fehlberg 7/8 integrator (relative and absolute tolerance 1e-12) until the first close approach to
Earth (CAA); compute the trans-halo impulse dV_EPO as the difference of the geocentric velocity at the CAA from the speed of
the chosen circular EPO. Four unknowns [v, dxdot, dydot, dzdot], all in the DE population (p6). Objective (eq. 6, p6):

    OBJ = W_h * |CAA_achieved - CAA_desired| / (E-M distance) + W_v * |dV_HOI|,   dV_HOI = sqrt(dxdot^2 + dydot^2 + dzdot^2)   (eq. 7)

with average Earth-Moon distance 384,400.0 km. DE parameters (p7): np = 40, F = 0.5, CR = 0.8, W_h = 10, W_v = 0.5, initial
step h = 0.01, convergence tolerance eps = 1.0E-5, search bounds for velocity components [-3000 m/s, 3000 m/s]; FORTRAN95
code, Linux, Intel Core i5 2.5 GHz, 8 GB. Random numbers from the GFORTRAN RAND generator.
INFERRED: because the trajectory is integrated backward from the target orbit, the Earth-departure time (phase of the
primaries at departure) is not a free variable; it is set by the insertion true anomaly v.

## 3. Every table, transcribed

### Table 1 (p7): comparison of transfers to a CRTBP halo orbit (benchmark, e = 0)

| Parameter | Direct transfer reported in Rausch (2005) | Direct transfer using the proposed technique |
| --- | --- | --- |
| Az of target halo orbit | 15,000 km | 15,000 km |
| CAA/size of circular EPO | 200 km | 200 km |
| Total dV required | 3.700 km/s | 3.6355 km/s |
| Flight duration | 4.19 days | 4.08367 days |

Benchmark halo orbit (p7, text): L1, Az = 15,000 km, northern family, initial condition
[0.8229570002125, 0.0, 0.04241855682667, 0.0, 0.152044631998602, 0.0] (the last digits of the first and last entries are small
in the scan; treat as READ to about 12 digits), period 12.05848 days. Search bounds for HOI velocity components [-3000, 3000]
m/s. Result in the text (p7): HOI impulse 586.48 m/s, flight duration 4.08367 days, total dV 3.6355 km/s from a 200 km
circular EPO. Other quoted references (p7): Parker and Born (2008) manifold transfer to a similar halo, 3.7289 km/s for
8.2 days; Mingotti et al. (2011) halo with Az = 8000 km, 200 km EPO, 3.659 km/s for 67.3 days.
Internal inconsistency in the paper: the same paragraph ends "The proposed technique produces a transfer that needs a lower
total velocity impulse (3.5289 km/s) for a shorter flight duration (4.088956 days)", which differs from Table 1 and the
preceding text (3.6355 km/s, 4.08367 days). We do not know which is intended; the benchmark closing figure in Table 1 and
the first statement agree with each other.

### M5N2 primary case (section 5.2, pp7-9)

READ. Target MR halo orbit M5N2, northern, average Az = 47,924 km (Rithwik and Ramanan 2019), initial conditions
[0.851666641652152, 0.0, 0.183285539178136, 0.0, 0.25828972225268, 0.0], period 4 pi (p8; also Table 8). CAA altitude 185 km,
EPO 185 km circular. Optimal two-impulse transfer (p8): total dV = 3.2802 km/s, of which 0.4882 km/s for MR halo insertion,
flight duration 4.92426 days, insertion at Moon true anomaly 359.82 degrees ("the apogee of the third revolution by the
spacecraft"); 107 m/s lower than the lowest value (3.388 km/s) of Peng and Xu (2015b).
CRTBP comparison orbit (p8): halo with Az = 47,924 km (the same average amplitude) with initial condition
[0.845272317414636, 0.0, 0.170480760011887, 0.0, 0.265160667222108, 0.0], period 11.459994 days; optimal transfer total dV =
3.4270 km/s (of which 0.4678 km/s for halo insertion), insertion at the apogee of the orbit, flight duration 4.79290 days.
The ERTBP result is lower by "about 147 m/s", attributed to "the inclusion of effect of eccentricity".

### Table 2 (p9): DE performance with different weights (M5N2, seed 5055)

| (W_h, W_v) | CAA obtained (km) | dV_HOI (km/s) | Iterations for convergence | Computational time (s) |
| --- | --- | --- | --- | --- |
| (10, 0.1) | 184.967 | 0.48824 | 406 | 1645 |
| (10, 0.2) | 185.062 | 0.48823 | 333 | 1387 |
| (10, 0.5) | 185.019 | 0.48823 | 187 | 905 |
| (10, 1.0) | 185.157 | 0.48822 | 172 | 702 |
| (10, 2.0) | 186.373 | 0.48822 | 198 | 946 |

### Table 3 (p9): different search bounds for velocity components

| Search bounds (m/s) | HOI location (deg) | dV_HOI (km/s) | Iterations | Computational time (s) |
| --- | --- | --- | --- | --- |
| [-3000, 3000] | 359.82 | 0.48822 | 187 | 905 |
| [-2000, 2000] | 362.14 | 0.48826 | 171 | 767 |
| [-1000, 1000] | 362.17 | 0.48821 | 355 | 702 |

### Table 4 (p9): different seeds for the random number generator

| Seed | HOI location (deg) | dV_HOI (km/s) | Iterations | Computational time (s) |
| --- | --- | --- | --- | --- |
| -5055 | 359.82 | 0.48822 | 187 | 905 |
| -258410 | 362.20 | 0.49322 | 189 | 954 |
| -8545523 | 362.17 | 0.49124 | 179 | 925 |

READ (p9) the paper's reading: "the DE-based algorithm converges to nearly the same solution irrespective of the search
bounds for velocity components and seeds". INFERRED, for context: the three seed rows differ by about 5 m/s (0.48822 to 0.49322
km/s) so "nearly the same" is at the few-m/s level, and the weight table shows the HOI cost agrees to 2e-5 km/s. The HOI
constraint is "the third revolution of the spacecraft around the apogee" (p9). Note that 359.82 and 362.1x degrees differ by about
2.3 degrees in the true-anomaly variable.

### Section 5.2.2 to 5.2.4 results quoted in text (pp10-11), with figures (not tabulated)

READ. Transfers to different insertion locations on M5N2, flight duration constrained below 10 days (Figs 5, 6): dV_EPO between
2.758 and 2.877 km/s, dV_HOI between 0.488 and 1.1 km/s, total dV between 3.280 and 3.747 km/s, flight durations between
3.92 and 5.75 days; five peaks in the curves matching the five revolutions. Peng and Xu (2015b) manifold transfers to M5N2:
total dV between 3.388 and 3.934 km/s, flight durations 56 to 71 days. Direct comparison at the insertion location
v = 126.648 degrees (from Fig. 10 of Peng and Xu 2015b): their flight duration 57.42011 days and total dV 4.6121 km/s; this
paper's transfer at that location 3.8215 km/s, "less by about 787 m/s". (Arithmetic on the printed pair gives 790.6 m/s; a small
rounding or typo-level difference.) Fixed flight durations 4 and 5 days, locations at 30 degree intervals (Fig. 7): difference
between minimum and maximum HOI impulse 445.4 m/s (4 days) and 640.4 m/s (5 days); maximum difference in total dV about
192 m/s. Different flight durations (Fig. 8, 5 to 90 days): maximum total dV about 4.90 km/s at 15 days, cyclic trend.

### Table 5 (p12): characteristics of the transfers in Figs 9 and 10 (M5N2, 185 km CAA)

| Orbit | Flight duration (days) | dV_TOTAL (km/s) | No. of close passes to Earth | Difference between each close pass (days) | CAA of each close pass (km) |
| --- | --- | --- | --- | --- | --- |
| a | 6.0 | 3.6566 | 1 | - | 185 |
| b | 15.0 | 4.900 | 2 | 8.41 | 185, 3982 |
| c | 36.9 | 3.955 | 5 | 8.27, 8.74, 7.67, 8.17 | 185, 616, 764, 5862, 1465 |
| d | 54.5 | 3.4183 | 7 | 8.23, 7.29, 7.93, 9.72, 8.28, 8.48 | 185, 955, 2550, 1008, 3852, 11732, 14242 |

READ (p11): the trajectories "are in the neighborhood of Earth for most of the flight duration", close passes repeat "at an
interval of 7-10 days", the spacecraft reaches "up to a maximum distance of 60 Earth radii".

### Table 6 (p13): optimal transfer to different revolutions of M5N2 (CAA 185 km)

| Revolution number | HOI point (true anomaly, deg) | dV_HOI (km/s) | dV_EPO (km/s) | dV_TOTAL (km/s) | Flight duration (days) |
| --- | --- | --- | --- | --- | --- |
| 1 | 71.960 | 0.50566 | 2.8304 | 3.3360 | 5.1476253 |
| 2 | 217.342 | 0.46744 | 2.8273 | 3.2947 | 4.8984034 |
| 3 | 359.820 | 0.48828 | 2.7917 | 3.2802 | 4.9241521 |
| 4 | 506.894 | 0.46785 | 2.8206 | 3.2885 | 5.1687698 |
| 5 | 648.561 | 0.45535 | 2.8293 | 3.2847 | 4.9714350 |

READ (p13): minimum HOI impulse occurs at the apogee of each revolution (Earth distance, Fig. 11 labels
approximately 3.81E05, 3.787E05, 3.823E05, 3.785E05, 3.81E05 km at the five insertion points) and at the closest approach to L1
(Fig. 12 labels approximately 6.022E04, 5.908E04, 6.08E04, 5.901E04, 6.205E04 km; these label digits are small in the scan).
INFERRED: row 3 of Table 6 equals the section 5.2.1 result (3.2802 km/s, 359.82 degrees; the flight duration 4.9241521 versus
4.92426 in the text differ in the fourth decimal place only, a text rounding).

### Table 7 (p14): M5N2 transfers for different CAAs

| CAA (km) | HOI point (true anomaly, deg) | dV_HOI (km/s) | dV_EPO (km/s) | dV_TOTAL (km/s) | Flight duration (days) | Velocity at CAA (km/s) |
| --- | --- | --- | --- | --- | --- | --- |
| 185 | 359.82 | 0.48822 | 2.7917 | 3.2802 | 4.924152 | 10.7373 |
| 200 | 73.947 | 0.45451 | 2.8171 | 3.2716 | 4.976615 | 10.7225 |
| 500 | 75.620 | 0.44524 | 2.5840 | 3.0292 | 5.094445 | 10.4894 |
| 1000 | 505.705 | 0.43795 | 2.2062 | 2.6441 | 4.937150 | 10.1113 |
| 2000 | 506.808 | 0.42438 | 1.5753 | 1.9997 | 5.040203 | 9.4808 |
| 3000 | 506.798 | 0.41456 | 1.0431 | 1.4576 | 5.031087 | 8.9484 |
| 4000 | 649.421 | 0.40764 | 0.5895 | 0.9972 | 5.055669 | 8.4945 |

Note: row 1 here prints dV_HOI 0.48822 and Table 6 row 3 prints 0.48828 for the same orbit and location; both are printed.

### Table 8 (p14): initial conditions and period of target MR orbits (mu = 0.0122, e = 0.0554 assumed)

| MR orbit and class | x0 | z0 | ydot0 | Period |
| --- | --- | --- | --- | --- |
| M2N1 Lyapunov | 0.804125156956177 | 0.000000000000000029 | 0.31182413982453 | 2 pi |
| M3N1 halo | 0.875404052867664 | 0.201620653468241 | 0.21551063837955 | 2 pi |
| M4N2 Lyapunov | 0.804504659626012 | 0.00000000000000020 | 0.31826866733409 | 4 pi |
| M4N2 halo | 0.895820897947402 | 0.194415672884192 | 0.34702042423622 | 4 pi |
| M5N2 halo | 0.851666641652152 | 0.183285539178136 | 0.25828972225268 | 4 pi |

The two Lyapunov z0 entries are numerical zeros of the planar orbits (the count of leading zeros, about 1e-17 and 2e-16, is
the least certain part of the scan). The M5N2 row matches the M5N2 initial conditions in section 5.2 digit for digit.
Initial-state phase is not printed: INFERRED that the Moon is at an apse at t = 0 (p2 requirement), and between perigee
(v = 0) and apogee (v = pi) the paper gives no choice; both should be tried.

### Table 9 (p15): transfer designs for different MR orbits (CAA 185 km, 2 to 10 days flight)

| MR orbit and class | HOI point (true anomaly, deg) | dV_HOI (km/s) | dV_EPO (km/s) | dV_TOTAL (km/s) | Flight duration (days) | Velocity at CAA (km/s) |
| --- | --- | --- | --- | --- | --- | --- |
| M2N1 Lyapunov | 91.056 | 0.44931 | 2.7858 | 3.2351 | 4.3708185 | 10.6910 |
| M3N1 halo | 66.025 | 0.39580 | 2.8377 | 3.2335 | 5.0100358 | 10.7422 |
| M4N2 Lyapunov | 88.504 | 0.42921 | 2.7840 | 3.2132 | 4.2512113 | 10.6895 |
| M4N2 halo | 237.636 | 0.60700 | 2.7953 | 3.4023 | 5.2342646 | 10.7465 |
| M5N2 halo | 359.820 | 0.48828 | 2.7917 | 3.2802 | 4.9241521 | 10.7373 |

READ (p14): "the MR halo orbit M3N1 is preferable over the MR Lyapunov orbit M2N1 because the total velocity impulse required for
transfer is less. Also, among orbits with N = 2, the MR Lyapunov orbit M4N2 has the least dV_TOTAL for transfer."

## 4. Method note: what the velocity numbers mean, and a caution

READ: all dV values are in km/s, flight durations in days. The paper prints no Earth gravitational parameter, Earth radius,
or time/length/velocity scale factors. INFERRED caution for reproduction: subtracting the printed dV_EPO from the printed
"velocity at CAA" gives an implied EPO speed of about 7.9456 km/s for the 185 km M5N2 rows and about 7.905 km/s for the other
rows of Tables 7 and 9 that we checked (the two M4N2 and M3N1 cases lie between those). A circular speed at 185 km or 200 km
altitude is lower than either. We have not resolved this; the paper does not say how the EPO speed is defined (frame, radius
or scale constants). dV_HOI and flight duration, which depend only on the dynamics and the objective, are the more trustworthy
reproduction targets; dV_EPO and dV_TOTAL need the unprinted constants.

## 5. Conclusions as printed (p15)

READ: "without the use of manifold theory, the transfers from Earth to halo orbits could be generated"; "The minimum total
velocity impulse required for transfer to the MR halo orbit M5N2 is 3.2802 km/s for a flight duration of 4.92426 days, which
are lower than the values reported in the literature (3.388 km/s and 65 days, respectively)". (The 65 days here, and 56 to 71
days in section 5.2.2, are the manifold-based values of Peng and Xu 2015b; the paper prints both forms.) Variation of the HOI
impulse with CAA only about 80 m/s from 4000 km down to 185 km. "Any location on MR halo orbit can be reached with a flight
duration between 4 and 6 days."

## 6. What this means for the project

INFERRED throughout.

Positive-control candidates (ordered by how well-defined the recipe is):
1. Table 8 initial conditions as a periodicity control for `core/er3bp.py`. Recipe: mu = 0.0122 and e = 0.0554 literal;
   start at the printed [x0, 0, z0, 0, ydot0, 0] with the primaries at an apse (try true anomaly 0 and pi), integrate for the
   printed period (2 pi for M2N1 and M3N1; 4 pi for the M4N2 and M5N2 rows) in true-anomaly time, and check state closure and
   the x-z plane crossing at half period (ydot, x-dot, z-dot pattern of p5). These are five 15-digit, same-model targets with
   a stated period and an integer resonance. They are quoted with DE-solver accuracy that the paper does not state; closure
   to the order of the printed digits cannot be assumed, only to the orbit's own design tolerance. The M5N2 halo is the most
   valuable: Peng and Xu 2015b is the comparison source for the M5N2 transfers (INFERRED that the same Earth-Moon M5N2 orbit
   is the one treated by Peng and Xu 2015a and 2015b; this paper does not say so).
2. The e = 0 benchmark halo (p7): mu = 0.0122, ICs [0.8229570002125, 0, 0.04241855682667, 0, 0.152044631998602, 0], period 12.05848
   days (the paper gives the period in days only). A plain CR3BP half-period differential
   correction check; independent of e.
3. The CRTBP comparison orbit for M5N2 (p8): [0.845272317414636, 0, 0.170480760011887, 0, 0.265160667222108, 0], period
   11.459994 days.
4. Transfer numbers (Tables 5, 6, 7, 9; 3.2802 km/s, 4.92426 days) require implementing the unprinted EPO and unit constants
   (section 4) and a DE search whose seeds give HOI costs differing by a few m/s (Table 4); weak as a golden test, usable as
   an order-of-magnitude check. Do not use dV_TOTAL as an exact target.
A caution in the project's own terms: items 1 to 3 are DE-designed orbits from the authors' earlier papers; the physical-orbit
check is whether they close in the project's own propagator, which is exactly the independent test wanted.

Overlap with open tasks: #435 and #437 (e > 0 continuation and fold-aware continuation of multi-revolution families). Direct
overlap through the targets: Table 8 gives five Earth-Moon MR orbits (M2N1 Lyapunov, M3N1 halo, M4N2 Lyapunov, M4N2 halo,
M5N2 halo) at one (mu, e) point. These are same-model landing points a fold-aware continuation from the CRTBP orbit family
should be able to reach or explain; if the continuator reaches the M5N2 halo at e = 0.0554 and the printed IC closes, that is
independent evidence that the continuation lands on a branch that exists in the real Earth-Moon problem. The paper also says
(p2) that a DE single-segment search "captured multiple solutions for each of the halo/Lyapunov MR orbits", a statement of multiplicity of solutions; no folds or bifurcation values are printed. It contains no
stability indices, no Floquet multipliers, no continuation in e.

Peng and Xu 2015 gap: this paper does NOT fill it. It cites Peng and Xu 2015a only in the references and in the descriptive
sentences quoted in section 2 above (it describes the method, not the tables), and it quotes numbers from Peng and Xu 2015b
(ASR 55:1015) and its Fig. 10: transfer costs 3.388 to 3.934 km/s, durations 56 to 71 days, location v = 126.648 degrees with
57.42011 days and 4.6121 km/s. Those numbers are from 2015b, not from the corrupt 2015a paper (the stability paper).
It does show, indirectly, that 2015a and 2015b are Earth-Moon papers using the M5N2 label, with M = revolutions about L1 and
N = primary revolutions (p5).

Catalogue classes: nothing here is a cycler, quasi-cycler or flyby-bearing resonant periodic orbit. The MR halo orbits are
resonant periodic orbits about L1 (commensurable with the primaries' period) but have no flybys of a third body; the transfer
trajectories (Fig. 9, Table 5) make repeated Earth close passes at 185 km to several thousand km, but they are single transfers,
not recurring cycles. Out of catalogue scope; relevant to the ERTBP infrastructure (periodic-orbit targets and the M:N
resonance bookkeeping) only.

## 7. References cited by the paper on multi-revolution halo orbits and elliptic-problem continuation (full as printed, pp15-16)

The reference list prints journal, volume and page only (no titles) for most items; titles below are filled in from the text
of the paper where it names the work, otherwise left blank.
- Peng, H., Xu, S.: Celest. Mech. Dyn. Astron. 123, 279 (2015a), https://doi.org/10.1007/s10569-015-9635-2 [Stability of two groups of multi-revolution elliptic halo orbits in the ERTBP; the corrupt-copy paper].
- Peng, H., Xu, S.: Adv. Space Res. 55, 1015 (2015b), https://doi.org/10.1016/J.ASR.2014.11.013 [direct transfers to MR halo orbits via stable manifolds, Earth-Moon ERTBP; source of the 3.388 to 3.934 km/s, 56 to 71 day comparison numbers].
- Peng, H., Xu, S.: Astrophys. Space Sci. 357, 1 (2015c) [low-energy transfers to MR halo orbits using a patched Sun-Earth/Moon CRTBP and Earth-Moon ERTBP model; page number as printed].
- Rithwik, N., Ramanan, R. V.: Design of halo orbits in the framework of Elliptic Restricted Three Body Problem using differential evolution. 11th IAA Symposium on the Future of Space Exploration, International Academy of Astronautics (2019) (source of the M5N2 orbit and its Az = 47,924 km).
- Rithwik, N., Ramanan, R. V.: J. Astrophys. Astron. 42 (2021), https://doi.org/10.1007/s12036-020-09651-w (single-segment DE design of MR orbits).
- Nath, P., Ramanan, R. V.: Adv. Space Res. 57, 202 (2016), https://doi.org/10.1016/j.asr.2015.10.033 (the Sun-Earth CRTBP direct two-impulse method this paper adapts).
- Ferrari, F., Lavagna, M.: Nonlinear Dynamics 93, 453 (2018), https://doi.org/10.1007/s11071-018-4203-4 (resonant ERTBP periodic orbits categorised by number of revolutions).
- Baresi, N., Olikara, Z. P., Scheeres, D. J.: J. Astronaut. Sci. 65, 157 (2018) (invariant tori and quasi-periodic orbits in Earth-Moon, QPT theory).
- Pernicka, H. J.: The numerical determination of nominal libration point trajectories and development of a station-keeping strategy, Ph.D. thesis, Purdue University (1990) (quasi-periodic orbits in Sun-Earth ERTBP, two-step differential correction).
- Broucke, R.: AIAA J. 7, 1003 (1969), https://doi.org/10.2514/3.5267 (stability of planar ERTBP periodic orbits).
- Moulton, F. R., Buchanan, D., Buck, T., Griffin, F. L., Longley, W. R., MacMillan, W. D.: Periodic Orbits. Carnegie Inst., Washington (1920).
- Lei, H., Xu, B., Sun, Y.: Adv. Space Res. 51, 917 (2013).
- Rausch, R. R.: Earth to halo orbit transfer trajectories. MS thesis, Purdue University (2005).
- Parker, J. S., Born, G. H.: J. Astronaut. Sci. 56, 441 (2008). Mingotti, G., Topputo, F., Bernelli-Zazzera, F.: J. Guid. Control Dyn. 34, 1644 (2011). Topputo, F.: Celest. Mech. Dyn. Astron. 117, 279 (2013). Cipriano, A. M., Dei Tos, D. A., Topputo, F.: Front. Astron. Space Sci. 5, 29 (2018).
- Storn, R., Price, K.: J. Global Optim. 11, 341 (1997) (differential evolution).
- Richardson, D. L.: J. Guid. Control Dyn. 3, 543 (1980), https://doi.org/10.2514/3.56033. Howell, K. C., Pernicka, H. J.: Celest. Mech. 41, 107 (1987). Paffenroth, R., Doedel, E., Dichmann, D.: AAS/AIAA Astrodynamics Specialist Conference (2001), AAS 01-303 (pseudo-arclength continuation of libration-point orbits).

Of these, for the project the targets are Peng and Xu 2015b and 2015c (corpus status not checked here; check CORPUS_INDEX before filing) and
Rithwik and Ramanan 2019 and 2021 (the source of the Table 8 orbits and of the DE design method).

## 8. Slips noted, stated factually

- Section 5.1 text (p7) gives 3.5289 km/s and 4.088956 days where Table 1 and the preceding sentences give 3.6355 km/s and
  4.08367 days; the two pairs are not reconciled in the paper.
- p10 "787 m/s" versus the difference of the printed totals (4.6121 minus 3.8215 km/s, 790.6 m/s); p9 "147 m/s" and p8 "107 m/s"
  agree with the printed totals.
- Table 7 row 1 and Table 6 row 3 give dV_HOI of 0.48822 and 0.48828 km/s for the same case.
- The flight duration 4.92426 days (text, p8 and p15) versus 4.9241521 or 4.924152 (Tables 6, 7, 9).
Each is of the size expected from separate runs or rounding and does not affect the paper's conclusions.
