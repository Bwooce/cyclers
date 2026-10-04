# Digest: Leiva & Briozzo 2005, fast periodic transfer orbits in the Sun-Earth-Moon Quasi-Bicircular Problem

Date: 2026-10-04. Task: `#884` literature follow-up (see `docs/notes/2026-10-04-884-literature-check.md`) and the `#892` positive-control search.

Paper: A. M. Leiva and C. B. Briozzo, "Fast periodic transfer orbits in the Sun-Earth-Moon Quasi-Bicircular
Problem", Celestial Mechanics and Dynamical Astronomy 91:357-372 (2005), DOI 10.1007/s10569-004-7818-3,
16 pages. Filed in the private paper corpus as
leiva-briozzo-2005-fast-periodic-transfer-orbits-sun-earth-moon-quasi-bicircular-problem-cmda-91-357-doi-10.1007-s10569-004-7818-3.pdf.

Marking: READ = stated or printed in the paper (journal page given). INFERRED = my deduction, not printed.
COMPUTED = a short check I ran for this note (stated with its result). Page numbers are the journal's
(357-372). The paper has no numbered tables; every printed number is in the running text and is
transcribed below digit by digit from the page images. No digit was unclear.

## 0. Summary

- READ (abstract, p357): "Starting from the identification and classification of a family of fast periodic transfer orbits in the Earth-Moon planar circular Restricted Three Body Problem (RTBP), and using analytic continuation techniques, we find two unstable periodic orbits in the Sun-Earth-Moon Quasi-Bicircular Problem (QBCP). The orbits found perform periodic Earth-Moon transfers with a period of approximately 29.5 days."
- READ: one RTBP family (a single symmetric family, found by a grid search on one section), one member of it with period exactly one Sun synodic period T_sun = 6.7911939, continued in epsilon to the QBCP at two Sun epochs (t0 = 0 and t0 = T_sun/2, chosen by a symmetry argument). The result is two QBCP periodic orbits of period T_sun, both unstable (|s1| = 5.496, 5.531), mild by the standards of this literature.
- READ: planar only, QBCP only. The authors say these orbits do "not approach the primaries as close as the orbits in our previous work" (p359).
- INFERRED (strong, from the printed h range): the orbits are very high-energy for the project's purposes (Hamiltonian h about +0.0105, Jacobi constant C = -2h about -0.021 for the RTBP member; the project's cycler families sit at C of about 3.0 to 3.2). They are not members of the project's C32, C31, C21, C11 or C33 families and are not in Table 1 of the 2008 follow-up. See sections 4 and 6.
- READ: this paper is the first of a pair. The 2008 paper (digest: `docs/notes/2026-10-04-digest-leiva-briozzo-2008-rtbp-to-qbcp-periodic-transfer-orbits.md`) generalises it from one hand-picked orbit and a symmetry-chosen epoch to 34 low-energy candidate orbits with the epoch chosen by a first-order necessary condition. Same Hamiltonian family, same units, same frame.

## 1. The model (Sect. 2, p359-363)

READ unless marked.

- Units (p360): "Adimensional units are used, with mE = 1 - mu and mM = mu, distance between the primaries unity, and orbital period 2 pi (giving unit angular frequency). For the Earth-Moon system this gives mu ~ 0.012150, time units of ~104 h (one sidereal month/2 pi), length units of ~384,400 km, and velocity units of ~1024 m/s. The values for all Solar System parameters are taken from the Jet Propulsion Laboratory (JPL) ephemeris." mu is printed to four significant figures only, so closure tests are mu-sensitive at the 1e-5 level of mu (see section 5).
- Two frames, both synodic and rotating (p360): for the analysis of the RTBP family, Earth at (x_E, y_E) = (-mu, 0), Moon at (x_M, y_M) = (1 - mu, 0); for the continuation to the QBCP, Earth at (mu, 0), Moon at (-1 + mu, 0). So the QBCP stage is the frame rotated by pi from the RTBP stage, the same as the 2008 paper. Section 4.3 (p368) repeats: "we have now passed from the synodic coordinate system with xE = -mu and xM = 1 - mu of Sections 4.1 and 4.2, to the synodic coordinate system with xE = mu and xM = -1 + mu."
- RTBP Hamiltonian (Eq. 1, p360): H_RTBP = (1/2)(px^2 + py^2) + y px - x py - (1 - mu)/r1 - mu/r2, with px = xdot - y, py = ydot + x, r1^2 = (x - xE)^2 + y^2, r2^2 = (x - xM)^2 + y^2; "H_RTBP = h, the Jacobi integral". So h is the Hamiltonian value (negative at low energy) and C = -2h. COMPUTED: using this formula the printed localisation orbit (x0 = 1.107569, ydot0 = -1.644251, mu = 0.0121505816, Earth at -mu) gives h = -0.245295 against the printed -0.245294, so the sign and normalisation of h are confirmed.
- Quasi-bicircular primaries (Sect. 2.2.1, p361-362): Jacobi coordinates r (Earth to Moon) and R (Earth-Moon barycentre to Sun), complex z and Z, Eqs. 2-3, Kepler relations Eq. 4, Fourier solution Eqs. 5-6 in powers of epsilon = a_i/a_e ~ 1/389. Constants in RTBP units (p362): "G = 1, mE + mM = 1, mS = 328900.54, alpha = mu, beta = 1 - mu, a_i = 1, a_e = 388.81114, n_i = 1, and n = 1 - n_e = 0.9251959855". The b_j and c_j are "in (Andreu, 1998)".
- Auxiliary functions (Eq. 7, p362): alpha_k(t) = alpha_k0 + sum_{j>=1} alpha_kj cos(j n t) for k = 1, 3, 4, 6, 7 and alpha_k(t) = sum_{j>=1} alpha_kj sin(j n t) for k = 2, 5, 8.
- QBCP Hamiltonian (Eq. 8, p363): H_QBCP = (1/2) alpha_1 (px^2 + py^2) + alpha_2 (x px + y py) + alpha_3 (y px - x py) + alpha_4 x + alpha_5 y - alpha_6 ((1 - mu)/r1 + mu/r2 + mS/rS), with px = (xdot - alpha_2 x - alpha_3 y)/alpha_1, py = (ydot - alpha_2 y + alpha_3 x)/alpha_1, r1^2 = (x - mu)^2 + y^2, r2^2 = (x + 1 - mu)^2 + y^2, rS^2 = (x - alpha_7)^2 + (y - alpha_8)^2. Note: as printed here r1 is the Earth distance (Earth at +mu) and r2 the Moon distance (Moon at -1 + mu); this differs from Eq. 1 (where r1 is also the Earth and r2 the Moon, but in the other frame) only through the frame. The alpha_6 factor multiplies the whole potential, the same as the 2008 paper, the 2018 Jorba-Cusco et al. paper and the corrected `core/qbcp.py` (INFERRED from the printed form; the module itself was not checked here).
- Coefficients: "can be found in Table 1.5 of (Andreu, 1998; see 1.5.1.)" (p363). Not printed in this paper; the number of Fourier terms used is not stated.
- Sun sense and phase (p363, Fig. 2b p362): "Note that in the Earth-Moon synodic system the motion of the Sun is retrograde, and that at t = 0 the primaries are collinear on the x axis in the sequence Moon-Earth-Sun (larger masses to the right)." Fig. 2 caption: "In (b) the Sun (mS) is shown at t = 0, when crossing the x axis." So in the paper frame the Sun is on the +x axis at t = 0 and moves clockwise.
- Periodicity (p363): "The Hamiltonian for the infinitesimal mass is T_sun-periodic with T_sun = 6.7911939 (i.e. one solar month in RTBP time units), thus in the QBCP periodic orbits in the synodic system must be k T_sun-periodic with k integer." COMPUTED: 2 pi / 0.925195985520347 = 6.79119387, so T_sun = 2 pi / n to the printed 8 digits (the 7th decimal rounds up).
- Continuation family (Eq. 9, p364): H = H_RTBP + epsilon (H_QBCP - H_RTBP), 0 <= epsilon <= 1. The same as Eq. 5 of the 2008 paper. As there, epsilon scales the whole difference of Hamiltonians, not the Sun's mass alone.
- Slip (INFERRED, harmless): on p362 the text prints "alpha = mu, beta = 1 - mu", whereas p361 defines alpha = mE/(mE + mM) and beta = mM/(mE + mM), which in these units are 1 - mu and mu. The printed assignment looks swapped, probably a typesetting or transcription slip. It has no consequence for the results, since the QBCP stage uses Andreu's coefficients and not alpha and beta directly.

## 2. The three-body (RTBP) family (Sect. 3-4.2, p363-367)

### 2.1 How the family was found (Sect. 4.1, p365)

READ. Section Sigma = {y = 0, xdot = 0, ydot < 0} in the frame with Earth at -mu ("we looked for orbits crossing the x axis perpendicularly"). Grid of initial conditions (x0, y0, xdot0, ydot0) = (1.1 + m Dx, 0, 0, -n Dy'), Dx = Dy' = 1e-3, m and n positive integers, "so looking for retrograde orbits passing behind the Moon (but not too close to it)". Integrator: variable-step Bulirsch-Stoer, relative precision 1e-14. Candidates: those returning to Sigma within 1e-3, with return time tau in (5, 8) because the search was for T_sun-periodic orbits (k = 1). Newton-Raphson refinement to a fixed point within relative error 1e-10. Symmetry: any orbit found crossing perpendicularly is symmetric (Szebehely 1967).

The selected orbit (p365): "x0 = 1.107569 and ydot0 = -1.644251. This PO has a Jacobi constant h = -0.245294 and a period tau = 6.372441, fairly close to T_sun." This is the seed, not the T_sun-periodic member.

### 2.2 Reconstruction of the family (Sect. 4.2, p365-367)

READ. Continuation in h in both directions, steps Dh = +/-1e-5, y0 = xdot0 = 0 kept (symmetric), x0 taken from the previous orbit, ydot0 from the new h by Eq. 1, Newton-Raphson refinement to relative precision 1e-10 in the return to Sigma, until "the natural termination of the family, which in this case begins and ends at a collision".

- Range (p366): "The family lies in the range hb <= h <= he, with hb = -0.770396 and he = 0.244529." (INFERRED: in C = -2h, this is C from 1.5408 down to -0.4890.)
- Shape (Fig. 3, p366; panels at h = -0.769996, -0.496536, -0.206536, 0.077513): "The family begins at hb with a collision with mE (Figure 3a). As h grows, this collision develops into two loops around mE, whose closest distance to mE increments progressively (Figure 3b, c, d). At the same time, the closest distance to mM diminishes, until at h = he the loop passing behind mM degenerates into a collision."
- Classification (Fig. 4, p366-367): h(x_Sigma) has two branches; "the leftmost branch starts at h = hb with x = xE with infinite slope, and grows monotonously with x until h = he for x ~ 0.8 with vanishing slope. The rightmost branch starts at h = hb for x = 1.23 with infinite slope, and decreases monotonously with increasing x until h = he for x -> xM. Thus the family is open." (Fig. 4 is in the Earth-at-minus-mu frame.)
- Period-in-family (Fig. 5, p367): "a monotonously growing function of h. So, there is only one member of this family with period T_sun which can be continued onto the QBCP." Read from the plot (INFERRED, plot reading only): T* rises from about 6.27 at hb to about 7.85 at he; the T_sun point is at h slightly above zero.
- Stability (p366): "the stability changes at h1 = -0.034536 and h2 = 0.244514. For hb <= h < h1 the fixed points are elliptic (stable orbits); for h1 < h < h2 we have reflection hyperbolic points (unstable POs); and for h2 < h <= he the fixed points are again elliptic." The T_sun member lies between h1 and h2, so it is unstable: "It is clear from this figure that the RTBP PO from which we will start the analytic continuation is unstable, but a priori this does not imply the stability or instability of the QBCP orbit to be obtained." (p367)

### 2.3 The T_sun-periodic member: printed initial condition (Sect. 4.3, p368)

READ. In the frame with Earth at +mu, Moon at -1 + mu (the continuation frame), the RTBP T_sun-periodic orbit is given by (p368): x0 = -1.02379270, y0 = 0, xdot0 = 0, ydot0 = 1.91110553. Its period is T_sun = 6.7911939 (printed to 8 digits). It is the x-axis crossing with ydot > 0 and x < x_M (that is, behind the Moon, beyond it as seen from the Earth).

Its Jacobi constant is not printed. COMPUTED from the printed state with Eq. 1 (mu = 0.0121505816): h = +0.010463, C = -0.020926. This lies between h1 = -0.034536 and h2 = 0.244514 as required. COMPUTED: integrating the printed state in the plain RTBP for exactly 2 pi / n (rtol 1e-13, mu = 0.0121505816) returns to within 4.7e-5 of the initial state (mostly in xdot); the miss is 6.2e-5 with mu = 0.0121505482 and 3.2e-4 with mu = 0.01215. So the printed RTBP state is a good, not a high-precision, T_sun-periodic orbit; the paper does not claim a closure figure for it. COMPUTED: along this RTBP orbit the minimum distance to the Earth's centre is 131,710 km and to the lunar surface 12,079 km (lunar radius 1737.4 km, 384,400 km per length unit).

## 3. The two periodic orbits in the QBCP (Sect. 4.3 and 5, p367-370)

### 3.1 Choice of epoch (Sect. 4.3, p367-368)

READ. The initial point is the x-axis crossing with ydot > 0 and x < x_M. To preserve the RTBP symmetry (x, y, xdot, ydot) -> (x, -y, -xdot, ydot) when epsilon > 0, "taking into account the chosen position of the initial point on the orbit, this amounts to asking invariance under the symmetry (x, y, t) -> (x, -y, -t), and the QBCP Hamiltonian will be invariant under this symmetry only if t0 = 0 or t0 = T_sun/2, that is, if at the initial time all three primaries lie on the x axis." Hence two continuations, "taking t0 = 0 in one case and t0 = T_sun/2 in the other."

Method: epsilon increased in steps D-epsilon = 3.5e-5; at each step the previous initial state was integrated from t0 to t0 + T_sun (Bulirsch-Stoer, relative error 1e-14) and Newton-Raphson refined until the return was within 1e-8; polynomial extrapolation from ten successive values allowed jumps in epsilon up to 3e-3; at epsilon = 1 a final refinement to relative error below 1e-11 (Sect. 4.3, p368). The start epoch was held fixed along the continuation.

### 3.2 The two orbits: everything printed (p368-370)

READ. Frame: Earth at (mu, 0), Moon at (-1 + mu, 0) (Sun on the +x axis at QBCP time zero, sequence Moon-Earth-Sun). Variables are positions and time derivatives (x, y, xdot, ydot) of the synodic coordinates; INFERRED that the velocities are time derivatives and not canonical momenta, because Eq. 8 defines the momenta separately from the velocities and the paper writes xdot, ydot throughout (the momentum py = (ydot + alpha_3 x)/alpha_1 at y = 0 differs from ydot at the percent level, so the reading matters; the coordinator has shown that the 2008 Table 2 orbits close under the time-derivative reading).

| | Orbit 1 | Orbit 2 |
|---|---|---|
| start time t0 | 0 | T_sun/2 (= 3.39559695 with T_sun = 6.7911939) |
| x0 | -1.01950751115 | -1.01940558303 |
| y0 | 0 | 0 |
| xdot0 | 0 | 0 |
| ydot0 | 1.97782573253 | 1.97615899365 |
| period | T_sun = 6.7911939 | T_sun = 6.7911939 |
| closure | relative error below 1e-11 | below 1e-11 |
| minimum distance to the Moon's surface | 10,431 km | 10,400 km |
| speed at that point | 2023 m/s | 2024 m/s |
| minimum distance to the Earth's centre | 134,588 km | 134,381 km |
| speed at that point | 2443 m/s | 2446 m/s |
| stability parameters | abs(s1) = 5.496, abs(s2) = 2.069 | abs(s1) = 5.531, abs(s2) = 2.070 |

Printed text (p368): "The two POs found are retrograde, and they are T_sun-periodic with a relative error < 1e-11. This means the infinitesimal mass returns to within 4e-3 m of its initial position, and within 1e-8 m/s of its initial velocity." (Check: 1e-11 times 384,400 km is 3.8e-3 m; 1e-11 times 1024 m/s is 1e-8 m/s. Consistent.) Printed (p369): "Both orbits are unstable, the absolute values for the stability parameters being |s1| = 5.496, |s2| = 2.069 for t0 = 0 and |s1| = 5.531, |s2| = 2.070 for t0 = T_sun/2."

Checks (COMPUTED): the Moon sits at x = -1 + mu in the paper frame, so the lunar-surface distance of the t0 = 0 state is |x0 + 1 - mu| x 384,400 km - 1737.4 km = 10,431.9 km (mu = 0.0121505816), matching the printed 10,431 km to about 1 km (0.5 km of this is the printed rounding of the lunar radius, which is not stated in the paper). For orbit 2 the same formula gives 10,392.8 km against the printed 10,400 km; the 7 km difference is INFERRED to be the Earth-Moon distance scale factor of the pulsating QBCP frame at t = T_sun/2 (the homothety is not unity there), which the formula omits. The check does not discriminate between mu = 0.01215, 0.0121505 and 0.0121505816 (they differ by 0.03 km). It also confirms that the minimum lunar distance occurs at the start point, an x-axis crossing behind the Moon.

Stability (p369-370). Each orbit's largest stability parameter is about 5.5, so the largest eigenvalue is about 5.3 (INFERRED: eta + 1/eta = 5.496 gives eta = 5.30). The authors tested this numerically with a 4-ball of 832 initial conditions of radius 1e-7 around the fixed point at t0 = 0: "even after 6 periods the distance to the fixed point had grown at most to ~6e-2". Over a longer time (Fig. 7): "15 T_sun is approximately the destabilization time for both orbits ... starting the orbit with an error <= 1e-11 as in this case, we can expect that after a number of periods n ~ 11/log(5.3) ~ 15 the errors grow to order unity." Both orbits behave alike to 15 T_sun and then differ; "after destabilizing, both orbits make some passages very close to Earth."

### 3.3 The figures

- Fig. 3 (p366): four panels of the RTBP family in the Earth-at-minus-mu frame, axes -1.5 to 1.5. described from the page image (INFERRED, shapes only): (a) h = -0.769996: two large lobes side by side; (b) h = -0.496536: two large loops crossing with a small loop at the crossing; (c) h = -0.206536: a three-loop figure; (d) h = 0.077513: a large outer loop enclosing two smaller inner loops (the family member closest to the T_sun one, which has h about 0.0105).
- Fig. 6 (p369): the t0 = 0 orbit over one period in QBCP synodic coordinates (Moon left, Earth right, axes -2 to 2 for the orbit, the Sun on the top and right axes in AU). Solid curve: the orbit, a large outer loop with two smaller inner loops around the Earth; circles mark the orbit at times labelled 0 to 6, squares mark the Sun's position at the same times (a dotted circle of radius 1 AU). The Sun starts on the +x axis (label 0) and the labels 1 to 6 run clockwise (lower right, bottom, lower left, upper left, top, upper right), consistent with a retrograde Sun; the angular spacing of the labels (about 51 degrees) suggests the labels are T_sun/7 apart (INFERRED from the picture; the caption says only "at the indicated times"). The orbit's label 0 is at the far left, at x of about -1.02, behind the Moon.
- Fig. 7 (p370): the two orbits (top: t0 = 0; bottom: t0 = T_sun/2) over 0 to 15 T_sun (left; they look like the Fig. 6 curve, slightly thickened) and 15 to 20 T_sun (right; erratic loops, several close to Earth).
- Fig. 4 and Fig. 5 as in section 2.2.

### 3.4 Are they cycler-type orbits?

INFERRED from the printed numbers and Fig. 6; the paper does not use the word cycler. In the paper's own sense they are "transfer" orbits: "we refer to an orbit as a transfer one if it encircles both primaries and passes between them, even if it does not pass particularly close to either" (p359). Facts:

- They encircle both primaries: yes, the outer loop passes behind the Moon at x of about -1.02 (the Moon is at -0.988) and the Earth is inside the loops.
- Lunar encounters per period: one close pass per T_sun, at the symmetric start point t0 (10,431 km and 10,400 km above the surface, speed about 2.02 km/s in the rotating frame). The rest of the loop stays far from the Moon (INFERRED from Fig. 6). Orbit 2 starts at the x-axis crossing half a period later but the geometry is the same, so also one pass.
- Earth encounters: none close. Minimum Earth distance 134,588 km and 134,381 km, about 0.35 length units, about 21 Earth radii; the orbit never approaches the Earth the way a low-energy Earth-Moon cycler does.
- Energy: the RTBP member has C of about -0.02, far above the L1, L2 and L3 energies (C of 3.19, 3.17 and 3.01); the project's cycler families are at C of about 3.0 to 3.2 (see the 2008 digest, section 8). Such orbits sit well outside the low-energy neck regime.

So they are periodic orbits with one flyby of the Moon per synodic month and no Earth flyby; I would not call them cyclers in the project's sense (a near-ballistic repeating itinerary with flybys of both bodies), though they are periodic Earth-Moon "transfers" in the authors' sense. The authors themselves say (p359) their orbit "does not approach the primaries as close as the orbits in our previous work".

## 4. Relation to the 2008 paper, the 2006 atlas and the Ross and Roberts-Tsoukkas "(3,2)" family

READ: the paper identifies the family only as "a family of unstable periodic orbits in the Earth-Moon RTBP which, though much faster than the orbits in (Leiva and Briozzo, 2005) and having none of the disadvantages mentioned above, does not approach the primaries as close as the orbits in our previous work" (p359). The "previous work" is Leiva & Briozzo, "Control of chaos and fast periodic transfer orbits in the Earth-Moon CR3BP", Acta Astronautica (in press in 2005). The 2006 atlas of families (numbered families such as 357 and 037, which the 2008 paper cites as "Leiva and Briozzo 2006b") post-dates this paper, so this paper contains no family number. Hence "which family of the 2008 Table 1 or the 2006 atlas" is not stated in the paper.

INFERRED:

- Not in the 2008 Table 1. The 2008 candidate list is restricted to "low energies (h <= -1.58617) and periods shorter than 6 months", every row has h between -1.59413 and -1.52791, and the T_sun-periodic orbit here has h of about +0.0105. The whole 2005 family has h from -0.770 to +0.2445, so no member of it meets the 2008 energy criterion. No 2008 row has tau/T_sun = 1. The 2008 paper (as digested) does not discuss the 2005 orbits at all, apart from referring to the 2005 paper for the alpha_kj.
- Not the Ross and Roberts-Tsoukkas "(3,2)" family. The project's C32 family is at C of about 3.15 to 3.18 and period about 2.4 to 2.8 T_sun (2008 digest, section 8); the 2005 orbits are at C of about -0.02 and period exactly 1 T_sun, with Earth distances above 134,000 km. Neither the energy nor the period nor the lunar-pass count matches. The 2005 paper therefore cannot be the source of the statement that the (3,2) member "persists even under solar perturbation". Combined with the 2008 finding (C32 appears there only as periodic arcs, Table 3), no printed Leiva and Briozzo result shows a C32 true periodic orbit in the QBCP. Which paper Ross and Roberts-Tsoukkas had in mind still cannot be settled from these two papers; the 2006 atlas paper (Leiva and Briozzo 2006b, not yet held as far as I know) and the Acta Astronautica 2005 paper are the remaining candidates.
- Relation to the 2008 method. 2005: one orbit, epoch chosen by the reversing symmetry (t0 = 0 or T_sun/2), orbit found by a coarse grid on a section; 2008: many orbits, epoch chosen by the first-order necessary condition (their Eq. 31) with four phases spaced T_sun/4 apart. The two 2005 epochs (0 and T_sun/2) are two of the four 2008 times for q = 1 and a symmetric orbit (INFERRED, since for a symmetric orbit the symmetric phases satisfy tan(2 phi) = 0 in the 2008 formula, giving phi = 0 and pi/2 modulo pi, that is times 0, T_sun/4, T_sun/2, 3 T_sun/4; the 2005 paper uses the two that preserve the symmetry t -> -t).

## 5. Positive controls for the project's corrected modules

All conversions to the project frame (Earth at -mu, Moon at 1 - mu, as in `src/cyclerfinder/core/qbcp.py` and `core/bcr4bp.py`): (x, y, xdot, ydot)_project = -(x, y, xdot, ydot)_paper. The Sun's sense of motion is clockwise in both. The paper's Sun on the +x axis at QBCP time zero becomes the Sun on the -x axis (sequence Sun-Earth-Moon) at QBCP time zero in the project frame. Per the 2008 note, `core/qbcp.py` evaluates the Fourier series at the printed time (the clock t is Andreu's, with the Sun on the paper-frame +x axis at t = 0), so the start time carries the phase; no separate phase input is needed.

Control Q1 (QBCP, t0 = 0). Start state in the project frame: x = +1.01950751115, y = 0, xdot = 0, ydot = -1.97782573253 at t = 0. Integrate one period T = 2 pi / 0.925195985520347 = 6.79119387 (use this and not the printed 6.7911939, which differs by 3e-9 and would add about 3e-9 times the speed to the miss). Expected: return to the initial state; printed closure relative error below 1e-11, instability multiplier 5.3 per period, so rounding of the printed 12 digits alone gives about 5e-11. INFERRED: with the corrected alpha tables the residual should be dominated by coefficient truncation and by mu (printed only as 0.012150), probably between 1e-9 and 1e-6; a miss of 1e-3 or more would indicate a model error of the kind found in `#891` and `#892`. A wrong Sun sense or a wrong epoch (INFERRED) should miss by of order the solar perturbation over a period, 1e-3 to 1e-1. Secondary checks: minimum distance to the lunar surface 10,431 km at the start (printed 10,431; the instant is a perilune); minimum Earth distance 134,588 km; abs(s1) = 5.496, abs(s2) = 2.069 from the monodromy matrix of the four-dimensional map over T (s = eigenvalue + 1/eigenvalue).

Control Q2 (QBCP, t0 = T_sun/2). Start state in the project frame: x = +1.01940558303, y = 0, xdot = 0, ydot = -1.97615899365 at t = T/2 = 3.395596936 (T as above, half period). Integrate one period T. Expected as for Q1; abs(s1) = 5.531, abs(s2) = 2.070; lunar-surface distance near 10,400 km (the pulsating scale factor shifts the plain formula by about 7 km); Earth distance 134,381 km.

Both Q1 and Q2 test the epoch and velocity reading that the 2008 closure result also used (state at the printed time, velocities as time derivatives), and they do so with an unambiguous epoch: the paper states "t0 = 0" and "t0 = T_sun/2" explicitly, and the symmetric state has y = xdot = 0. They are the best-conditioned printed QBCP control found so far: multiplier 5.3 against 155 or more for the 2008 Table 2 orbits, and 12 printed digits against 9. The pair also gives a symmetry test: the reversing map (x, y, xdot, ydot, t) -> (x, -y, -xdot, ydot, -t) with t0 = 0 should map the orbit to itself.

Control B1 (RTBP, Sun off; tests the CR3BP stage only). State in the project frame: x = +1.02379270, y = 0, xdot = 0, ydot = -1.91110553, mu = 0.0121505816, period T. COMPUTED: closure 4.7e-5 (see section 2.3). The printed state is a good approximation of the T_sun-periodic RTBP orbit and nothing more; do not use it to certify a solver at the 1e-9 level. The RTBP orbit found by the project's own corrector from this state (period near T) is the right comparison, with C = -2 x 0.010463 = -0.020926 as the identification check.

Control B2 (the BCR4BP, `core/bcr4bp.py`). Not a precision control: INFERRED from the 2008 digest (Eq. 6 there, |H_QBCP - H_BCP| of order 0.03 |H_BCP|) and from the module docstring (O(epsilon^2) deviation) that the incoherent BCR4BP closes these orbits only to the size of the QBCP-versus-BCP difference, not to 1e-8. Use: Sun phase theta_sun0 = pi (project frame, Sun on -x at t = 0) for Q1; for Q2 the Sun angle at the start is pi - n T/2 = 0, that is theta_sun0 = 0. Compare the QBCP residual with the BCR4BP residual: a QBCP miss of 1e-8 against a BCR4BP miss of 1e-3 to 1e-2 is expected, while a BCR4BP miss of 1 or more would indicate a Sun-sense or phase error.

Not reproducible from what is printed: the alpha_kj table and number of terms (Andreu 1998, Table 1.5), mu beyond four figures, the lunar radius used for the "distance to the Lunar surface", the Earth-Moon distance scale used for the km values.

## 6. What this means for `#884`

Already in this paper (READ):

- The idea of continuing an Earth-Moon RTBP periodic orbit, whose period is commensurate with the Sun's synodic period, to the QBCP by a homotopy in epsilon, with a start epoch chosen by symmetry. The method is the 2008 method in its simplest case. For the project this is prior art for the outline of `#884` (period commensurate with T_sun, epoch from symmetry, continuation to the physical Sun), but only for a different, high-energy family and with k = 1.
- Retrograde, symmetric, T_sun-periodic QBCP orbits exist for that family at both symmetric epochs (t0 = 0 and T_sun/2), unstable with a largest multiplier of about 5.3, lifetime of about 15 periods under a 1e-11 error (p370).
- The same-epoch pair behaving alike for 15 periods and differently after destabilisation (p370): evidence that the two symmetric-phase equivalents of one RTBP orbit are distinct orbits that are close in phase space, which is the project's "two equivalents at the symmetric phases" observation in a high-energy case (INFERRED parallel).

Not in this paper:

- Any cycler-class family of the project (C32, C31, C21, C11, C33 or a catalogued Ross, Braik or Casoliva orbit), or any orbit at C near 3. Note for novelty: this paper adds no prior art against the project's low-energy `#884` results beyond what the 2008 paper already contains.
- Anything about the Ross and Roberts-Tsoukkas (3,2) statement. The 2005 orbits are not (3,2); the statement remains unsupported by the printed Leiva and Briozzo papers.
- A Melnikov-type or first-order phase condition (introduced in 2008), resonances with q >= 2, periodic arcs, folds or continuation in the Sun's mass, three dimensions, BCR4BP, an eccentric Moon.
- The alpha_kj tables.

Use for the corrected-model effort (`#891`, `#892`): this paper's two orbits supply a third independent printed QBCP control (after the 2018 POL1 and POL2 and the 2008 Table 2 orbits), with stated epochs, a modest instability multiplier (about 5.3) and 12 printed digits. The expected result is closure far better than the 2008 rows allow. Run them in `core/qbcp.py` with the printed time as the clock. They exercise the Sun's sense, the coefficient tables and the epoch convention together at a lunar distance of 12,000 km, so they probe the model where the Sun's tidal field matters (large orbit, loop enclosing Earth and Moon, period equal to the Sun's synodic period), which the L1 and L2 libration-point controls do not.
