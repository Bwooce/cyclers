# Digest: Leiva & Briozzo 2008, extension of fast periodic transfer orbits from the Earth-Moon RTBP to the Sun-Earth-Moon QBCP

Date: 2026-10-04. Task: `#884` literature follow-up (see `docs/notes/2026-10-04-884-literature-check.md`).

Paper: A. M. Leiva and C. B. Briozzo, "Extension of fast periodic transfer orbits from the Earth-Moon
RTBP to the Sun-Earth-Moon Quasi-Bicircular Problem", Celestial Mechanics and Dynamical Astronomy
101:225-245 (2008), DOI 10.1007/s10569-008-9134-9, 21 pages. Filed in the private paper corpus as
leiva-briozzo-2008-extension-fast-periodic-transfer-orbits-earth-moon-rtbp-to-sun-earth-moon-qbcp-cmda-101-225-doi-10.1007-s10569-008-9134-9.pdf.

Marking: READ = stated or printed in the paper (journal page given). INFERRED = my deduction, not printed.
Page numbers are the journal's (225-245). Tables were transcribed digit by digit from the page images and
cross-checked against the PDF text layer; the two agreed everywhere except where flagged. Minus signs are
written as "-". Nothing in any table is computed by me.

## 0. Summary

- READ (abstract, p225): "Starting from 80 families of low-energy fast periodic transfer orbits in the Earth-Moon planar circular Restricted Three Body Problem (RTBP), we obtain, by analytical continuation 11 periodic orbits and 25 periodic arcs with similar properties in the Sun-Earth-Moon Quasi-Bicircular Problem (QBCP)."
- READ: the selection condition (Sect. 3, p231-233) is a first-order necessary condition; it exists only when the RTBP period is (p/q) times the Sun's synodic period with q = 1 or 2. It returns the Sun phase at which to start the continuation (Eq. 30/31).
- READ: planar only, QBCP only (no 3D, no eccentric Moon, no BCR4BP), continuation in epsilon of H = H_RTBP + eps (H_QBCP - H_RTBP) (Eq. 5), 34 candidate orbits, 6 continued to eps = 1 as orbits (11 phase-specific orbits), 25 more as "periodic arcs" (not true periodic orbits, see Sect. 5 below).
- READ: every orbit and arc found is unstable. The paper never uses the word "cycler" and never names the Ross/Braik family labels.
- INFERRED (strong, from Jacobi constants): the project's C32 family at the 5/2 resonance is the paper's families 180A_1 and 180A_2; the project's C31 family at 5/2 is the paper's family 357. In this paper these appear only as periodic arcs (Table 3), not as periodic orbits (Table 2). See section 7.

## 1. The model (Sect. 2, p227-229)

READ unless marked.

- Units (p228): "Nondimensional units are used, with mE = 1 - mu and mM = mu, distance between the primaries ai = 1, and orbital period 2 pi (giving a mean motion ni = 1). For the Earth-Moon system this gives mu ~ 0.0121505, time units of ~104 h (one sidereal month/2 pi), length units of ~384400 km, and velocity units of ~1024 m/s." Constants for the Solar System taken from JPL ephemerides.
- Frame (p228, Fig. 1): synodic, origin at the Earth-Moon barycentre, planar. Earth at (x_E, y_E) = (mu, 0), Moon at (x_M, y_M) = (-1 + mu, 0). So the Moon is on the left, the Earth to its right. The Sun (Fig. 1b, "shown at t = 0, when crossing the x axis") moves retrograde.
- RTBP Hamiltonian (Eq. 1, p228): H_RTBP = (1/2)(px^2 + py^2) + y px - x py - mE/rE - mM/rM, with px = xdot - y, py = ydot + x. Jacobi integral C = -2h = -2 H_RTBP. So h is the Hamiltonian value (energy-like, negative), C = -2h.
- BCP Hamiltonian (Eq. 2, p228): H_BCP = H_RTBP - mS/rS + (mS/a_e^2)(x cos(theta) + y sin(theta)), theta = Omega t + phi. Constants (p228): mS = 328900.54, a_e = 388.81114, Omega = -0.925195985520347, phi the initial Sun phase. The solar potential is arranged so a particle at the barycentre is in free fall (zero solar acceleration at x = y = 0).
- QBCP (Sect. 2.3, p229): Andreu's (1998, 2002, 2003) Fourier solution. Auxiliary functions (Eq. 3): alpha_k(t) = alpha_k0 + sum_{j>=1} alpha_kj cos(j n t) for k = 1, 3, 4, 6, 7; alpha_k(t) = sum_{j>=1} alpha_kj sin(j n t) for k = 2, 5, 8. Hamiltonian (Eq. 4):
  H_QBCP = (1/2) alpha_1 (px^2 + py^2) + alpha_2 (x px + y py) + alpha_3 (y px - x py) + alpha_4 x + alpha_5 y - alpha_6 ( mE/rE + mM/rM + mS/rS ),
  with px = (xdot - alpha_2 x - alpha_3 y)/alpha_1, py = (ydot - alpha_2 y + alpha_3 x)/alpha_1, rE^2 = (x - mu)^2 + y^2, rM^2 = (x + 1 - mu)^2 + y^2, rS^2 = (x - alpha_7)^2 + (y - alpha_8)^2 (subscripts 7 and 8 as read from the page image; consistent with Eq. 3, where alpha_7 is a cosine and alpha_8 a sine series), and n = -Omega = 0.925195985520347 "(the Sun's mean motion in the BCP)".
  This is the same Hamiltonian as `src/cyclerfinder/core/qbcp.py` (alpha_6 on the whole potential) -- INFERRED from the printed form, not from the project's module being checked here.
- Sun sense and phase (p229): "In the Earth-Moon synodic system the motion of the Sun is retrograde, and at t = 0 the primaries are collinear on the x axis in the sequence Moon-Earth-Sun (larger masses to the right)." The Hamiltonian is T_sun-periodic with T_sun = 6.7911939 "(i.e. one synodic month in RTBP time units)".
- Coefficient tables: NOT printed. The paper defers to Andreu (1998) and Leiva & Briozzo (2005) for the alpha_kj; the number of Fourier terms used is not stated. The Sun's mass, a_e and n are given only for the BCP; for the QBCP the paper says "Assuming the same units for masses, distance, and time introduced in Sect. 2.1".
- mu: printed only as ~0.0121505 (6 digits). INFERRED: the printed L1 abscissa -0.836915310 (Table 1 note) is reproduced by mu = 0.0121505482 (my root-find of the collinear equation; mu = 0.0121505 gives 0.83691555, mu = 0.0121505816 gives 0.83691515). So the tables were generated with mu near 0.01215055, not the project's 0.0121505816. The difference is 3.4e-8 and matters for closure tests at the 1e-7 level.
- Continuation family (Eq. 5, p230): H = H_RTBP + eps (H_QBCP - H_RTBP), 0 <= eps <= 1. For choosing the phase only, they replace it with Eq. 7, H = H_RTBP + eps (H_BCP - H_RTBP), justified by Eq. 6: |H_QBCP - H_BCP| <~ 0.03 |H_BCP| (stated "by direct numerical computation or by a careful analysis"). Note this eps scales the whole difference of Hamiltonians (primaries' motion and the Sun), unlike the project's continuation in the Sun's mass.

## 2. The selection condition (Sect. 3, p229-233)

READ. Derivation steps:

1. (p231) Hill-type truncation of the solar potential (Eq. 12): U_sun(r, t) = (1/2)(mS/a_e^3)[ r^2 - 3 (s_hat . r)^2 + O(a_e^-1) ], valid because r/a_e <~ 1/389. Eq. 13: H = H_RTBP + eps U_sun, T_sun-periodic in t.
2. A tau-periodic solution must have tau and T_sun commensurate: "say p T_sun = q tau with p and q co-prime. The minimal common period of U_sun and x(t) is then T* = p T_sun = q tau". Necessary condition (Eq. 14, p231): the work done by the Sun over one period T* vanishes, eps * integral_0^T* grad U_sun . rdot dt = 0.
3. Poincare-expansion in eps (Eqs. 21-24, p232), Fourier expansion of the RTBP orbit z(t) = x + i y (Eq. 25). The O(eps) integrals vanish identically unless (n + m) q - 2 p = 0, i.e. m = -n + 2p/q (Eq. 28), "which has integer solutions for m only when q = 1 or q = 2, since p and q are co-prime".
4. Result (Eq. 29): Im[ e^{-2 i phi} sum_n z_0n z_{0,-n+2p/q} ] = 0, equivalently (Eq. 30, p233)
   tan(2 phi) = Im[ sum_n z_0n z_{0,-n+2p/q} ] / Re[ sum_n z_0n z_{0,-n+2p/q} ].
   Quote (p233): "Note that this is a condition on the solar phase phi at t = 0, and provides the required phase for continuation as an explicit function of the Fourier coefficients z_0n of the unperturbed solution, which is the PO in the RTBP." And: "When q != 1, 2 no condition is obtained at first order in the eps expansion. ... This course will not be pursued in the present work".
5. Numerical form (Eq. 31, p235): discrete Fourier transform of N = 2^m samples (m = 10, 11, ... until nine significant figures stop changing; m = 12 sufficed in most cases), tan(2 phi) = Im[ sum_{n=0}^{N-1} z_0n z_{0,(-n+2p/q) mod N} ] / Re[ same ]. Because tan(2 phi + k pi) = tan(2 phi) there are four phases in [0, 2 pi): phi_k = phi_1 + k pi/2, k = 2, 3, 4 (p235: "four phase values in the interval [0, 2 pi), that we denote phi_i, i = 1, ..., 4, and which are related by phi_k = phi_1 + k pi/2, k = 2, 3, 4").
6. Conversion to a time (Eq. 32, p235): t_i = (-phi_i / Omega) mod T_sun, i = 1..4, "the corresponding times t_i at which the Sun crossed the positive x axis". "So instead of starting integrations at time t = 0 with solar phase phi_i, we started each integration at time t_i with solar phase phi_i = 0. ... we will designate each continuated PO by its number in Table 1 followed by t1, t2, etc."

Answers to the specific questions:

- Is it a Melnikov-type or symmetry argument? INFERRED: it is a first-order averaging (Poincare-Melnikov-type) condition, the vanishing of the first-order work integral of the quadrupole (tidal) term; the paper calls it only "necessary (though not sufficient)" (p230) and "only necessary and approximate" (p242). It is not a symmetry argument. The four admissible phases come from the quadrupole term varying as 2 phi (period pi in phi), hence the pi/2 spacing of the four t_i in time (T_sun/4 = 1.69779845 TU). This is the same structure as the project's #884 "Melnikov amplitude" argument; the project's zeros at {0, pi/2}, {pi/2, pi}, {pi/3, 2 pi/3} are the same family of conditions in a different phase convention.
- Role of q: the condition exists only for q = 1, 2, so resonances 8/3 (q = 3) and 5/6 of the project's survey are outside the paper.
- Phase is "kept constant along the continuation" (p235-236); the displacement of x0 is left unconstrained (p236).

## 3. Table 1 (p234): candidate RTBP orbits, 34 rows

READ. Caption: "Identification number, initial conditions h, y, and ydot, and period tau, of the POs in the RTBP selected to attempt continuation to the QBCP." Note under the table: "The initial conditions are given in synodic coordinates at the crossing with xdot < 0 of a Poincare surface Sigma at x = -0.836915310".

So: section Sigma = {x = -0.836915310} (the L1 abscissa in their frame, Moon at -0.98785, Earth at +0.01215), crossing with xdot < 0; x is fixed, y and ydot are tabulated, xdot is fixed by h through Eq. 1 (INFERRED: they print only (h, y, ydot) and the sign of xdot). h is the RTBP Hamiltonian value of Eq. 1, C = -2h (INFERRED consistency check: 2 x 1.59005198 = 3.18010396 and 2 x 1.57583831 = 3.15167662, matching the project's C32 values 3.18010 and 3.15168 to the printed digits). mu is not printed to better than 0.0121505 (see section 1). The identification numbers are the family numbers of the 2006 atlas (Sect. 4, p233: "The POs are identified by the number of the family in (Leiva and Briozzo 2006b) to which they belong"; "a trailing number has been added" to distinguish several in one family). tau/T is tau divided by T_sun = 6.7911939.

The 80 families are those with "both low energies (h <= -1.58617) and periods shorter than 6 months (tau <= 40 in synodic units)" (p233); "note that in this plot some families have been continuated to higher values of h" (p233) -- hence the rows with h > -1.58617 (357, 053_1, 077_2, 084_1, 180A_2, ...). Candidate selection is the 34 intersections of Fig. 2 with horizontal lines tau = (p/q) T_sun, q = 1, 2.

| Number | h | y | ydot | tau/T_sun |
|---|---|---|---|---|
| 053_1 | -1.56587846 | 0.0746237099 | 0.0382328315 | 5/2 |
| 053_2 | -1.58753537 | 0.0348725952 | -0.0149583173 | 5/2 |
| 077_1 | -1.58840213 | -0.0041344629 | 0.0619882192 | 5/2 |
| 077_2 | -1.57183324 | -0.0446025645 | 0.125546118 | 5/2 |
| 084_1 | -1.57183324 | 0.0485238196 | 0.0516615071 | 5/2 |
| 084_2 | -1.58840213 | 0.0357836953 | -0.0150779150 | 5/2 |
| 180A_1 | -1.59005198 | 0.00520342002 | 0.0479318298 | 5/2 |
| 180A_2 | -1.57583831 | -0.0283283340 | 0.117872065 | 5/2 |
| 357 | -1.52791268 | 0.129037155 | 0.115396082 | 5/2 |
| 146A | -1.59170073 | 0.00684686737 | 0.0389091502 | 3 |
| 018A | -1.58890359 | -0.0348199799 | 0.0299219506 | 7/2 |
| 144A | -1.58963300 | -0.0414498069 | -0.00372442160 | 7/2 |
| 147 | -1.59196333 | 0.00534538143 | 0.0491647641 | 7/2 |
| 157A_1 | -1.59096752 | 0.0111493986 | 0.0260247945 | 7/2 |
| 157A_2 | -1.58960218 | 0.0207084632 | 0.00687076626 | 7/2 |
| 187A | -1.59288010 | -0.0159741427 | -0.0215813955 | 7/2 |
| 187B_1 | -1.59003721 | 0.0220349307 | -0.0202753029 | 7/2 |
| 187B_2 | -1.59003721 | -0.0359084333 | -0.0350860129 | 7/2 |
| 188A_1 | -1.59064811 | 0.0209945181 | -0.0103532292 | 7/2 |
| 188A_2 | -1.59316603 | 0.00892212868 | 0.00616109122 | 7/2 |
| 194 | -1.59203408 | -0.0208044831 | -0.0327233527 | 7/2 |
| 244 | -1.59203408 | 0.000196827658 | 0.0549913768 | 7/2 |
| 305 | -1.59029796 | -0.0347974421 | 0.0281011828 | 7/2 |
| 013 | -1.58740571 | -0.0399746624 | -0.0441622383 | 4 |
| 020 | -1.58710323 | -0.0369266173 | -0.00906635816 | 4 |
| 021 | -1.58710323 | -0.0446655750 | -0.0470435962 | 4 |
| 171_1 | -1.58930426 | 0.0275869044 | -0.0153042817 | 4 |
| 171_2 | -1.59262227 | 0.00719782489 | 0.0385102976 | 4 |
| 081A | -1.58708764 | -0.0494577087 | -0.0323811597 | 9/2 |
| 172 | -1.58857380 | 0.0325459456 | -0.0264123993 | 9/2 |
| 032B_1 | -1.58703219 | -0.0280775876 | 0.00719954961 | 5 |
| 032B_2 | -1.58703219 | -0.0191844708 | -0.00487293371 | 5 |
| 287 | -1.59413574 | 0.00344758201 | -0.000982857301 | 5 |
| 238 | -1.59331122 | -0.00774529571 | -0.0196143210 | 11/2 |

Counts by tau/T_sun: 5/2: 9 rows; 3: 1; 7/2: 13; 4: 5; 9/2: 2; 5: 3; 11/2: 1 (total 34).

## 4. Table 2 (p238): initial conditions of the periodic orbits in the QBCP, 11 rows

READ. Caption: "Initial conditions in synodic coordinates for the 3, 4, and 5 T_sun-periodic orbits in the QBCP obtained by analytical continuation of 3, 4, and 5 T_sun-periodic orbits in the RTBP." Columns: Number, t_i, x, xdot (first line of each row) and Period, y, ydot (second line). The tabulated variables are therefore positions (x, y) and velocities (xdot, ydot) in the QBCP synodic frame (Moon on the left, Earth to the right), not canonical momenta. The momenta follow from the printed relation p_x = (xdot - alpha_2 x - alpha_3 y)/alpha_1 etc. (p229).

| Number | Period | t_i | x | xdot | y | ydot |
|---|---|---|---|---|---|---|
| 146A_t3 | 3T | 3.32657957 | -0.833881068 | -0.0658016572 | -0.00176277248 | 0.0366759999 |
| 146A_t4 | 3T | 6.72217651 | -0.833912081 | -0.0658968551 | -0.00207184538 | 0.0369412534 |
| 013_t3 | 4T | 1.92708674 | -0.841058432 | -0.0710601802 | -0.0415648661 | -0.0231934953 |
| 013_t4 | 4T | 5.32268367 | -0.841255581 | -0.0709901229 | -0.0417347404 | -0.0224675395 |
| 020_t1 | 4T | 1.37929365 | -0.840861762 | -0.0890884586 | -0.0359687313 | 0.0101908036 |
| 020_t2 | 4T | 4.77489058 | -0.841778236 | -0.0889071284 | -0.0358040150 | 0.0138070651 |
| 171_2_t3 | 4T | 1.73158585 | -0.837975852 | -0.0330903804 | 0.0107192583 | 0.0366443351 |
| 171_2_t4 | 4T | 5.12718279 | -0.838292090 | -0.0322885877 | 0.0111104984 | 0.0373817317 |
| 032B_1_t4 | 5T | 6.35357835 | -0.824049581 | -0.112839541 | -0.0276354646 | -0.0392670699 |
| 053d_2_t3 | 5T | 3.28799799 | -0.833040875 | -0.103358636 | 0.0365770176 | -0.0240383643 |
| 053d_2_t4 | 5T | 6.68359492 | -0.832339057 | -0.102739753 | 0.0372962615 | -0.0265393587 |

"053d" is printed with a superscript d (the 5T orbit 053_2 continued to a true orbit, as opposed to the (5/2)T arc 053_2 of Table 3).

Epoch of these states (READ + INFERRED). READ (p235): the RTBP initial condition of Table 1 was integrated "up to t = t_i, so obtaining the initial conditions at solar phase phi = 0, and the integration time was reset to zero for consistency with the QBCP model"; each QBCP integration started "at time t_i with solar phase phi_i = 0". p238 calls (x, y, velocities, initial time) the "initial conditions". INFERRED reading (A): the tabulated state is the QBCP state at the instant the Sun is on the positive x axis (QBCP time zero in Andreu's convention, sequence Moon-Earth-Sun), and t_i is the RTBP-orbit time elapsed since the Table 1 section crossing (the x values -0.82 to -0.85 are not on the Table 1 section x = -0.8369, consistent with a state propagated to a different time). The paper does not say this in one sentence. Alternative (B): state at QBCP time t_i. A closure test would tell which is right; (A) should be tried first.

Check I made on the printed t_i: within every orbit the two tabulated epochs differ by 3.3956 = T_sun/2 to the printed digits (for example 146A: 6.72217651 - 3.32657957 = 3.39559694; T_sun/2 = 3.39559695). So of the four phases phi_k only a T_sun/2-separated pair appears for each orbit; the other pair (T_sun/4 away) is not tabulated, presumably because it failed.

Precision: initial conditions are stated to fractional error below 1e-9 (p235: "resulting in all cases in a fractional error <1e-9 for the initial conditions of the continuated orbit and the components of M").

## 5. Tables 3 and 4 (p240): periodic arcs, 13 + 12 = 25 rows

What a "periodic arc" is (READ, p239): "In these cases we reduced the integration interval to [0, tau], obtaining -- continuation now succeeded -- periodic arcs in the QBCP instead of periodic orbits. These arcs are periodic in the sense that after a time tau equal to the period of the RTBP PO, their coordinates and velocities return to their initial values, but being q = 2 the solar phase at t = tau is pi instead of zero." Consequently the continued object is a fixed point of the map over tau = (p/2) T_sun with the Sun then at phase pi; it is not a periodic orbit of the QBCP (which needs p T_sun). For arcs "the integration of the evolution equations for the fundamental solution matrix gives only m(tau), instead of the monodromy matrix M = m(p T_sun)"; the s_i of Table 5 for arcs come from m(tau) and "still provide useful information on the arc stability". All arcs have q = 2 (periods 5/2 and 7/2 T_sun). Initial conditions are "given with fractional errors below 1e-9" (p239). Same column conventions as Table 2 (first line: Number, t_i, x, xdot; second line: y, ydot).

Table 3 (caption: "(5/2) T_sun-periodic arcs in the QBCP obtained by analytical continuation of (5/2) T_sun-periodic orbits in the RTBP"), 13 rows:

| Number | t_i | x | xdot | y | ydot |
|---|---|---|---|---|---|
| 053_1_t3 | 2.94372263 | -0.823291283 | -0.204967189 | 0.0729037633 | 0.0157937233 |
| 053_1_t4 | 6.33931957 | -0.822940134 | -0.205554219 | 0.0728395206 | 0.0160284110 |
| 053_2_t3 | 3.28799799 | -0.832638313 | -0.102970156 | 0.0369785750 | -0.0257859571 |
| 053_2_t4 | 6.68359492 | -0.832737737 | -0.103123506 | 0.0368880646 | -0.0248040354 |
| 077_1_t1 | 1.13012471 | -0.807338411 | -0.141817858 | -0.00903646943 | -0.00511122361 |
| 077_1_t2 | 4.82781242 | -0.839269824 | -0.0853997094 | 0.000509872176 | 0.0631047327 |
| 084_2_t3 | 3.34190762 | -0.837184965 | -0.0936535576 | 0.0340469056 | -0.0151707186 |
| 084_2_t4 | 6.73750455 | -0.837104209 | -0.0935333914 | 0.0341170567 | -0.0146415891 |
| 180A_1_t1 | 1.51986327 | -0.838273181 | -0.0700359357 | 0.0103226185 | 0.0420907706 |
| 180A_1_t2 | 4.91546020 | -0.837710299 | -0.0705053126 | 0.00978465371 | 0.0414271551 |
| 180A_2_t1 | 1.17849132 | -0.845704702 | -0.126198028 | -0.0202900217 | 0.146995322 |
| 357_t3 | 2.82792318 | -0.816461641 | -0.296940242 | 0.129114141 | 0.0853142743 |
| 357_t4 | 6.22352012 | -0.815847321 | -0.298659169 | 0.128445946 | 0.0850832278 |

Note: 053_2_t3 and 053_2_t4 here share t_i (3.28799799, 6.68359492) with the 053d_2 orbits of Table 2, but the states differ (e.g. x = -0.832638313 vs -0.833040875): the arcs are the (5/2)T versions and the 053d rows are the 5T orbits (p239: "We present them twice in these Tables").

Table 4 (caption: "(7/2) T_sun-periodic arcs ... "), 12 rows:

| Number | t_i | x | xdot | y | ydot |
|---|---|---|---|---|---|
| 018A_t1 | 1.36400740 | -0.838283891 | -0.0626291135 | -0.0327661089 | 0.0319939392 |
| 018A_t2 | 4.75960434 | -0.837560836 | -0.0634742122 | -0.0331814869 | 0.0289125536 |
| 144A_t3 | 1.89115633 | -0.837604853 | -0.0458798745 | -0.0426290387 | -0.00335455397 |
| 144A_t4 | 5.28675326 | -0.837485334 | -0.0461119947 | -0.0426245640 | -0.00379001677 |
| 147_t1 | 1.65992167 | -0.838588024 | -0.0358277014 | 0.00979166411 | 0.0511566376 |
| 147_t2 | 5.05551861 | -0.838334830 | -0.0364152988 | 0.00928545481 | 0.0511342220 |
| 187A_t1 | 0.731217088 | -0.837880215 | -0.0371102305 | -0.0371102305 | -0.0156566365 |
| 187A_t2 | 4.12681402 | -0.838137701 | -0.0374160755 | -0.0175746564 | -0.0143312619 |
| 188A_1_t3 | 3.38816739 | -0.839117139 | -0.0908157085 | 0.0174832264 | -0.00186802548 |
| 188A_1_t4 | 6.78376433 | -0.838350309 | -0.0899598374 | 0.0184230777 | -0.00458069786 |
| 188A_2_t1 | 0.157104258 | -0.834216258 | -0.0478467375 | 0.00272691115 | 0.00189692582 |
| 188A_2_t2 | 3.55270119 | -0.834027198 | -0.0473041023 | 0.00229340183 | 0.00144173781 |

FLAG, possible printing error: in Table 4 row 187A_t1 the printed xdot (-0.0371102305) and the printed y (-0.0371102305) are identical to all ten digits, which is not credible for two independent variables. Neither the page image nor the text layer distinguishes them. Do not use 187A_t1 as a control without testing both readings. (The row 187A_t2 has y = -0.0175746564, xdot = -0.0374160755.)

## 6. Table 5 (p241): minimal distances and stability parameters, 36 rows

READ. Caption: "Minimal distances d_E to the Earth and d_M to the Moon, and stability parameters s_i, for each periodic orbit and arc in the QBCP obtained by analytical continuation of periodic orbits in the RTBP". s_i = lambda_i + 1/lambda_i, k = 1, 2 (p230 and p239) where lambda are eigenvalues of M (orbits) or m(tau) (arcs). The paper does not say whether d_E and d_M are to the centre of the body or the surface; the d_M of 725 km for 032B_1_t4 is below the lunar radius (1737 km), consistent with the statement (p244) "some being Moon colliders" if centre distances. INFERRED.

| Number | d_E (km) | d_M (km) | abs(s1) | abs(s2) |
|---|---|---|---|---|
| 053_1_t3 | 111191 | 18702 | 5.9 | 1957.8 |
| 053_1_t4 | 111038 | 18728 | 5.9 | 1959.1 |
| 053_2_t3 | 121705 | 4593 | 62.6 | 6.1 |
| 053_2_t4 | 121558 | 4610 | 61.8 | 6.1 |
| 077_1_t1 | 118266 | 3110 | 220.8 | 0.3 |
| 077_1_t2 | 118325 | 3150 | 216.4 | 0.3 |
| 084_2_t3 | 119167 | 4683 | 29.5 | 33.1 |
| 084_2_t4 | 119052 | 4705 | 32.7 | 29.9 |
| 180A_1_t1 | 118458 | 7371 | 49.9 | 12.4 |
| 180A_1_t2 | 118471 | 7365 | 49.5 | 12.4 |
| 180A_2_t1 | 106956 | 23966 | 5153.0 | 5.2 |
| 357_t3 | 96931 | 14713 | 5365.7 | 5.2 |
| 357_t4 | 96755 | 14716 | 5377.8 | 5.2 |
| 146A_t3 | 121237 | 7311 | 336.6 | 29.1 |
| 146A_t4 | 121361 | 7346 | 341.1 | 29.3 |
| 018A_t1 | 133112 | 12440 | 371.2 | 34.3 |
| 018A_t2 | 133120 | 12448 | 368.8 | 34.4 |
| 144A_t3 | 137294 | 10432 | 8423.6 | 6.3 |
| 144A_t4 | 137295 | 10442 | 8344.1 | 19.2 |
| 147_t1 | 118788 | 3701 | 796.2 | 591.8 |
| 147_t2 | 118786 | 3696 | 859.6 | 603.3 |
| 187A_t1 | 126512 | 6669 | 4061.2 | 17.4 |
| 187A_t2 | 126668 | 6684 | 3809.2 | 19.8 |
| 188A_1_t3 | 122518 | 5908 | 25.6 | 7.9 |
| 188A_1_t4 | 122431 | 5816 | 29.1 | 8.0 |
| 188A_2_t1 | 124181 | 6473 | 94.7 | 9.2 |
| 188A_2_t2 | 124332 | 6519 | 97.4 | 9.1 |
| 013_t3 | 137125 | 2729 | 155.8 | 35.9 |
| 013_t4 | 137126 | 2733 | 155.8 | 35.9 |
| 020_t1 | 134724 | 3253 | 766.4 | 8.3 |
| 020_t2 | 134709 | 3250 | 767.7 | 8.3 |
| 171_2_t3 | 119714 | 4199 | 1330.9 | 28.6 |
| 171_2_t4 | 119714 | 4197 | 1339.6 | 29.7 |
| 032B_1_t4 | 129841 | 725 | 278.1 | 1.1 |
| 053d_2_t3 | 121518 | 4542 | 4031.9 | 32.4 |
| 053d_2_t4 | 121518 | 4542 | 3923.3 | 35.4 |

The 11 orbit rows are 146A_t3, 146A_t4, 013_t3, 013_t4, 020_t1, 020_t2, 171_2_t3, 171_2_t4, 032B_1_t4, 053d_2_t3, 053d_2_t4; the other 25 rows are arcs (READ from the table layout and Sect. 4.3).

Reading of the |s| values (READ p244): "All the periodic orbits and arcs found in the QBCP are unstable, but their stability parameters |s_i| vary greatly as shown in Table 5." Smallest max(|s1|,|s2|) among arcs: 188A_1_t3 (25.6, 7.9); among true orbits: 013_t3 and 013_t4 (155.8, 35.9) and 032B_1_t4 (278.1, 1.1). Statement (p244): d_E "approximately in the range 96000-138000 km"; d_M "most ... pass fairly close to the Moon"; capture around the Moon achievable "by applying modest braking impulses at the point of closest approach", minimal impulses from 10 m/s (013_t3, 013_t4) "to some 50 m/s in the worst cases".

## 7. Which orbits continued, which failed, and why (Sect. 4.2-4.3, p236-239)

READ.

- Totals: of the 34 candidates "Only six ... could be continued up to eps = 1 by the procedure of Sect. 4.1, and only for some of the four possible values of the solar phase given by Eq. 31, giving a total of 11 POs" (p236). The six RTBP orbits are 146A, 013, 020, 171_2, 032B_1, 053_2 (the last at 5T, as an orbit, for two phases). Then "proceeding as described above 25 periodic orbits failing to continuate as such could be continuated to periodic arcs" (p239; the paper's "periodic orbits" here means the RTBP orbits whose orbit continuation failed).
- Failure modes (p236-239), with examples from Fig. 3 (curves eta_i(eps) that do not reach eps = 1: 287_t1, 081A_t1, 305_t3):
  1. Failing to start: Newton-Raphson did not converge at the first nonzero eps even after reducing the step below 1e-4. Stated reading: "Eq. 31 gives only a necessary condition ... this kind of continuation failure is in fact a strong indication that it does not [satisfy the sufficient conditions]".
  2. Reaching a bifurcation or intersection: the slope d eta_i / d eps became infinite at some eps < 1 (example 287_t1). "No attempt to continue past these points has been made."
  3. Reaching a collision with the Moon at some eps < 1 (081A_t1, 305_t3). Regularisation would allow going past it; not attempted.
  4. Failing to find the iterated fixed point: instability grew with eps until accumulated numerical error pushed the iterate outside the Newton-Raphson convergence region. Remedy would be quadruple precision; not attempted. Section 4.3: "In general (p/q)T_sun-periodic orbits with q = 2 proved difficult to continue ... these orbits have p >= 5, thus T* >= 5 T_sun. The exponential divergence along the unstable manifold ... after such a long integration time"; "periodic orbits with q = 1 and p > 5 showed to be very difficult to continue." The remedy used: shorten the integration to [0, tau], giving arcs.
- Numerical method (p235-236): start eps step 1e-4; at each eps integrate 0 to p T_sun, build the Poincare map and monodromy matrix, Newton-Raphson to a fractional error < 6e-8; polynomial extrapolation lets steps grow to 1e-2; at eps = 1 a final Newton-Raphson to fractional error < 1e-9; Bulirsch-Stoer integrator, fractional accuracy 1e-14 per step. The phase is held fixed; x0 free (not constrained orthogonal to the orbit).
- "The only (p/2)T_sun-periodic orbit that we were able to extend to periodic orbits in the QBCP was the RTBP PO 053_2, for two values of the initial solar phase" (p239).
- No orbit in Table 2 has period 5/2, 7/2, 9/2 or 11/2 T_sun; the half-integer-period members survive only as arcs, if at all (no 9/2 or 11/2 row appears in Tables 2-4 at all: families 081A, 172, 238 all failed).

## 8. Correspondence to the project's cycler families (item 7)

INFERRED from Jacobi constants, with the periselene as a secondary check; the paper itself never names a cycler family. The 2006 atlas digest (`docs/notes/2026-07-28-744-broucke-leiva-barrabes-earth-moon-lineage-digest.md`) prints only families 357 and 037, so it does not map family numbers beyond 357. The project's C values are from `docs/notes/2026-10-04-884-sun-forced-em-cyclers.md` (family table, 5/2 rows); the paper's from C = -2h of Table 1.

| Project family (period) | Project C, periselene | Paper row | C = -2h (from Table 1) | Paper d_M, QBCP arc |
|---|---|---|---|---|
| C32 (Ross-RT 32 = Braik-Ross C32), 5/2 | 3.18010, 8,288 km | 180A_1 (h = -1.59005198) | 3.18010396 | 180A_1_t1: 7371 km; _t2: 7365 km |
| C32, 5/2 | 3.15168, 23,304 km | 180A_2 (h = -1.57583831) | 3.15167662 | 180A_2_t1: 23966 km |
| C31 (Ross-RT 31), 5/2 | 3.05583, 14,294 km | 357 (h = -1.52791268) | 3.05582536 | 357_t3: 14713 km; _t4: 14716 km |

Agreement of the Jacobi constants to the digits the project note prints (6 significant figures) at both C32 members and at C31, together with the matching 5/2 resonance and the periselene ordering (a small periselene for the higher-C member), makes the identification very likely but it is still an inference. It would be removed by one run of the project's CR3BP corrector from the Table 1 section data (see section 10).

Where this leaves the Ross and Roberts-Tsoukkas statement (quoted in `docs/notes/2026-10-04-884-literature-check.md`): the paper's two C32 members (180A_1, 180A_2) are NOT in Table 2 (they did not continue as periodic orbits); they appear only in Table 3 as periodic arcs (180A_1_t1, 180A_1_t2, 180A_2_t1), i.e. fixed points of the (5/2)T map with the Sun at phase pi at return. So what this paper actually shows for C32 is: a continuation in eps from the CR3BP member to the QBCP at eps = 1 exists for two Sun epochs of 180A_1 and one of 180A_2, as arcs, all unstable (|s1| = 49.9, 49.5, 5153.0). The phrase "persists even under solar perturbation" as a periodic orbit is not supported by the C32 rows of this paper. Which orbit Ross and Roberts-Tsoukkas meant cannot be decided from what is printed here: the only true QBCP periodic orbits are at 3T, 4T and 5T (Table 2), which are not the C32 period (the project's C32 family spans 2.42 to 2.83 T_sun, INFERRED from "T 16.46-19.25 TU" in the project note divided by 6.7911939). The 2005 companion paper (CMDA 91:357-372) is a candidate for the source and is not held.

Other orbits encircling both primaries with close lunar passes: every row in Tables 2-5 encircles both (the paper's loose sense of "transfer", p227: "we refer to an orbit as a transfer one if it encircles both primaries and passes between them, even if it does not pass particularly close to either"). Close lunar passes: d_M from 725 km (032B_1_t4) to 23966 km (180A_2_t1); orbits below 3000 km are 032B_1_t4 (725), 013_t3 (2729), 013_t4 (2733). I cannot map the remaining families (146A, 013, 020, 171_2, 032B_1, 053_2 and the arcs) to named cycler families from what is printed; the matches above rest on Jacobi constants, which are only available for the 5/2 rows.

## 9. What this means for #884

Already in this paper (READ):

- The method of the project's survey in outline: select Earth-Moon periodic orbit families by commensurability of their period with the Sun's synodic period; choose Sun phases by a first-order necessary condition; continue in a homotopy parameter to the full Sun; report the QBCP orbits that result and their stability. Planar only.
- The selection rule: tau = (p/q) T_sun with q = 1 or 2, four Sun phases pi/2 apart (as T_sun/4 apart in time) from Eq. 31.
- The finding that the continued orbits are unstable (all of them) and that continuation fails for several reasons (fold, collision, loss of convergence from instability).
- Specific members: 180A_1 and 180A_2 (C32 at 5/2) and 357 (C31 at 5/2), as arcs; no C21, C11 or C33 member appears to be involved (no 3/1 row, none at 3/2).

Not in this paper (READ or absence):

- Any cycler-family labelling, any use of the word cycler, any Ross-RT or Braik-Ross family.
- True periodic orbits for the 5/2 members of C32 and C31 (arcs only).
- Resonances with q >= 3 (8/3, 5/6 etc. of the project's survey): excluded by Eq. 28.
- Three-dimensional orbits, eccentric lunar orbit, BCR4BP or any ephemeris model.
- Continuation in the Sun's mass (this paper's eps scales H_QBCP - H_RTBP, which also switches on the primaries' non-circular motion and the QBCP's own frequencies) and a physical-Sun-mass fold analysis.
- Periselene, stability and forced-orbit census at all commensurate members of all catalogued families, and the Melnikov-amplitude analysis of the phases.
- The QBCP coefficient tables and the number of terms (printed nowhere in this paper).

For novelty: the claim "Sun-forced equivalents of Earth-Moon cycler-class orbits exist at the physical Sun" is not novel in outline (prior art: this paper, for q <= 2, planar, arcs for the half-integer members), and the project should not describe its 5/2 results for C32 and C31 as new in kind. What stands apart from this paper (if the project's recomputed results survive `#891`/`#892`) is: true periodic orbits (not arcs) at the 5/2 members, the 8/3 members, the C21/C11/C33 families, the 3D families, and the BCR4BP.

## 10. Positive controls

All conversions below to the project frame (Earth at -mu, Moon at 1 - mu; `src/cyclerfinder/core/qbcp.py` reflects x and y by pi relative to Andreu): (x, y, xdot, ydot)_project = -(x, y, xdot, ydot)_paper. The Sun sense is unchanged by the rotation (retrograde, clockwise in both). The paper's Sun at t = 0 on +x (sequence Moon-Earth-Sun) becomes the Sun on -x (sequence Sun-Earth-Moon, left to right) in the project frame.

Control A (three-body, no Sun; tests the CR3BP stage and the Jacobi sign): Table 1 rows 180A_1 and 180A_2 (C32 at 5/2) and 357 (C31 at 5/2).
- Recipe: mu = 0.0121505482 (the value implied by the printed L1 abscissa; also try the project's 0.0121505816 and report the difference). Section x = -0.836915310 (project frame +0.836915310), crossing with xdot < 0 in the paper frame (xdot > 0 in the project frame). Set y and ydot from Table 1; solve xdot from h through Eq. 1 with the stated sign. Period tau = (5/2) T_sun = (5/2) x 6.7911939 = 16.97798475 (with the printed T_sun; T_sun is printed to 8 digits, so a period tolerance near 1e-7 is warranted).
- Expected: C = -2h = 3.18010396 (180A_1), 3.15167662 (180A_2), 3.05582536 (357); closure of the periodic orbit at the printed 9-10 digit precision amplified by the monodromy multiplier (project note: 80 for the C = 3.18010 member, 4.3e3 for C = 3.15168). Independent of the project's own code: the Table 1 numbers are from the paper. Caveat: the project's C values were found by the project's walk, so agreement to 6 digits is an identification check, not a golden for the project's solver; the golden is the paper's (y, ydot, h) triple.

Control B (QBCP, true periodic orbit, smallest instability): Table 2 rows 013_t3 and 013_t4 (4T period, |s1| = 155.8, |s2| = 35.9, d_M 2729 and 2733 km, d_E about 137125 km).
- Recipe (reading A of the epoch): QBCP clock zero with the Sun on the +x axis in the paper frame (the project's Sun at its epoch zero convention, Moon-Earth-Sun in the paper frame). Initial state 013_t3: x = -0.841058432, y = -0.0415648661, xdot = -0.0710601802, ydot = -0.0231934953. Convert velocity to momenta with the printed p_x, p_y relations if the project integrates in canonical variables. Integrate T* = 4 T_sun = 27.1647756 and compare the final state with the initial one. Repeat for 013_t4 (x = -0.841255581, y = -0.0417347404, xdot = -0.0709901229, ydot = -0.0224675395). If reading A fails, repeat with the Sun phase at the start equal to Omega t_i (reading B), t_i = 1.92708674 or 5.32268367.
- Expected closure: the printed state has about 9 decimals; with an instability multiplier of about 155 for one period the residual from rounding alone is of order 5e-10 x 155 ~ 1e-7 (INFERRED), plus any model difference from the coefficient truncation. A residual of about 1e-7 supports the model; a residual of 1e-3 to 1e-2 means a model difference (or the wrong epoch reading). The paper's own closure is 1e-9 fractional.
- Other orbit candidates with printed numbers: 146A_t3 (3T, |s| = 336.6 and 29.1), 020_t1 (4T, 766.4 and 8.3), 171_2_t3 (4T, 1330.9 and 28.6); 032B_1_t4 (5T, 278.1 and 1.1) has d_M = 725 km and is a lunar collision candidate, so it is a poor control.

Control C (QBCP, periodic arc, the C32/C31 members relevant to the project): Table 3 rows 180A_1_t1 (|s1| = 49.9, |s2| = 12.4), 180A_1_t2, 180A_2_t1 (|s1| = 5153.0), 357_t3, 357_t4 (|s1| about 5366). The arc closes after tau = (5/2) T_sun = 16.97798475 with the Sun then at phase pi; the test is x(tau) = x(0) in the four phase-space variables, not a QBCP period. The most benign arc is 188A_1_t3 (Table 4, 7/2 T_sun; |s| = 25.6, 7.9) but it is not a C32 or C31 member. For C32, 180A_1_t1 is the most relevant (|s1| = 49.9, residual rounding amplification about 5e-10 x 50 = 3e-8, INFERRED).

Ranking (best first for reproducibility and relevance): (1) 013_t3 / 013_t4 for the QBCP model check; (2) 180A_1_t1 for the C32 5/2 link; (3) Control A rows for the CR3BP and the C32 identification. Do not use 187A_t1 (printing flag).

Not reproducible from what is printed: the QBCP alpha_kj table and number of terms (the project's tables must be the same as Andreu's 1998 / 2018 versions for the test to be meaningful, which `qbcp.py` states they are), and any member not in Tables 2-4.

## Addendum 2026-10-04 (`#896`): corrections found when the printed numbers were made tests

Test: `tests/core/test_leiva_briozzo_2008_tables.py` (commits `ebbb8d47`, `f0e757e1`). Every digit
of Tables 1 to 5 agrees with the page images (pp234, 238, 240, 241).

- Section 5, row 187A_t1 (the FLAG above): the misprint is in y, not xdot. Keeping the printed xdot
  and solving for y alone gives y = -0.0175001 and a return to 1.4e-5; keeping y and solving for
  xdot leaves 8.6e-2. The row is omitted from the tests.
- Section 4, epoch: the digest recommends reading (A) first. Reading (B), the state at QBCP clock
  time t_i, is the one that closes, for all 24 usable arcs (as for Table 2 in `test_qbcp.py`).
- Section 6, Table 5: d_E and d_M are centre distances in plain units of 384,400 km. All 35 d_E agree;
  surface and pulsating-scale readings fail. Seven d_M are printed larger than the computed minimum
  (1.4 to 7.5 km; 032B_1_t4 prints 725 km against 483.4 km, inside the Moon): strict expected
  failures; probably a minimum sampled at integration steps (INFERRED).
- Section 1, mu: confirmed. The section abscissa is L1 at mu = 0.0121505482 (1.5e-10); Table 1
  closes 12 to 1100 times better there than at the project default.
- Section 10, Control B expects about 1e-7 for 013_t3; measured 5e-6 at the paper's mu.
