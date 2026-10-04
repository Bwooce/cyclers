# Digest: Danby 1965, "Matrizant of Keplerian Motion"

Evidence tags: READ (journal page) is what the printed page says, read on 150 and 300 dpi page images (the matrices M and N on a 300 dpi crop); COMPUTED is my own numerical check of 2026-10-05 (scratch only, no test files); INFERRED is my reading across sources.

Reference: J. M. A. Danby (Yale University), "Matrizant of Keplerian Motion", AIAA Journal 3(4):769-770, April 1965 (Technical Note, received 29 October 1964; work supported in part by the Office of Naval Research and the Air Force Office of Scientific Research), DOI 10.2514/3.2976. Filed in the private paper corpus as `danby-1965-matrizant-keplerian-motion-aiaa-j-3-769-doi-10.2514-3.2976.pdf` (2 PDF pages: journal page 769 carries the end of an unrelated plasma-accelerator note, with Figure 3 and equations (11) to (15) of that note, which are not part of this paper; the Danby note starts at the bottom right of p.769 and occupies p.770). Neighbouring digest: `docs/notes/2026-10-04-digest-deprit-deprit-bartholome-1968-kepler-matrizants.md`.

## 1. What it is

A two-page note giving two closed-form expressions for the 6 x 6 state transition matrix ("matrizant", also "error matrix" or "guidance matrix") of three-dimensional Keplerian motion, valid, per the note, for "any type of orbit" and usable in the integration of perturbed orbits (reference to Danby 1964, AIAA J. 2:13). It replaces the longer formulas of Danby's earlier paper "The matrizant of Keplerian motion" (AIAA J. 2:16, 1964; "I was rash enough to comment that there was no general formula applicable to every type of orbit, except over time spans short enough for series expansions": Goodyear's formulas, which use Herrick's universal variables and are related to Bower's, were then sent to him; reference [2] is Goodyear and Guseman, IBM/RTCC Report 12-008, 1964, and [3] is Herrick, Astrodynamical Report 7, 1960; neither is held). The note gives no numerical example, no table, no figure and no derivation beyond a sketch ("eschewing the algebra, which is straightforward").

## 2. Formulation (READ p.770)

Frame and variables:
- 3D, Cartesian, inertial-time derivative, no universal variables, no Kepler equation: the matrices depend only on the current Cartesian position and velocity in the orbit frame, the elapsed time t (measured "from an arbitrary epoch"), the angular momentum per unit mass h, the gravitational parameter mu, and (for the inverse) the semimajor axis a.
- Orbital reference system (the one used in Danby's 1964 paper): X toward perihelion, Y in the orbit plane 90 degrees ahead, Z along the angular momentum vector. A prime denotes d/dt; X'' = -mu X / R^3, Y'' = -mu Y / R^3, R the distance from the focus.
- Relation (1): dW = M(t) dG = M(t) M^-1(t0) dW0. dW is the six-vector of position and velocity departures at t (Cartesian components, ordered X, Y, Z, X', Y', Z', positions then velocities), dW0 the same at t0, dG a vector of six independent geometrical-element changes.
- The elements (Eckert-Brouwer, as modified by Danby): the first set is dG^T = [dl0 + dr, dp, dq, e dr, da/a, de], with l0 the mean anomaly at an arbitrary epoch, dp, dq, dr infinitesimal rotations about the reference axes, a the semimajor axis, e the eccentricity (the paper writes "dr" for the rotation about the Z axis, and "delta r" is distinct from the distance R); the note observes that the combinations (HX + KX') and (HY + KY') of Eckert and Brouwer can be written a(YX'/h - 1) and -a X X'/h, with h the angular momentum per unit mass. The matrix M in (2) is for the modified set dG1^T = [(dl0 + eta dr)/n, a de/h, dp, da/2a, a e dr/h, dq], with n the mean motion and eta = (1 - e^2)^(1/2). The element sets are only an intermediate device: the final product M M0^-1 involves no elements except a (and h, mu).

## 3. The matrizant written out

### 3.1 Form 1: the orbital-frame matrix M, eq. (2) (READ p.770, 300 dpi crop)

Rows and columns ordered (X, Y, Z, X', Y', Z') for the rows and the six modified elements for the columns:

    row 1 (X)   :  X'    ,  Y X' - h    ,  0  ,  2X - 3 t X'   ,  Y Y'            ,  0
    row 2 (Y)   :  Y'    ,  -X X'       ,  0  ,  2Y - 3 t Y'   ,  -Y X' - 2h      ,  0
    row 3 (Z)   :  0     ,  0           ,  Y  ,  0             ,  0               ,  -X
    row 4 (X')  :  X''   ,  Y' X' + Y X'' ,  0  ,  -X' - 3 t X'' ,  Y'^2 + Y Y''   ,  0
    row 5 (Y')  :  Y''   ,  -X'^2 - X X'' ,  0  ,  -Y' - 3 t Y'' ,  -X' Y' - Y X'' ,  0
    row 6 (Z')  :  0     ,  0           ,  Y' ,  0             ,  0               ,  -X'

Inverse (eq. 2 and 2a): M^-1 = A Phi M^T Phi^T, with A = diag(a/mu, a/(mu h), 1/h, a/mu, a/(mu h), 1/h) (the paper prints "A = [a/mu  a/mu h  1/h  a/mu  a/mu h  1/h] I", read as a diagonal matrix of those six entries) and Phi = [[0, -I], [I, 0]] (I the 3 x 3 identity, so Phi is 6 x 6). The matrizant is R(t; t0) = M(t) M^-1(t0).

### 3.2 Form 2: the radial/transverse/normal matrix N, eq. (3), (4) (READ p.770, 300 dpi crop)

Displacements resolved along the radius vector (dR), the forward transverse direction in the orbit plane (dT), and the normal (dZ); dR' and dT' (and dZ') are the VELOCITY displacements RESOLVED along the R and T directions, "and not the errors in the radial and transverse components of velocity". (3): (dR, dT, dZ, dR', dT', dZ')^T = N N0^-1 (dR0, dT0, dZ0, dR0', dT0', dZ0')^T, with

    row 1 (dR)  :  R'     , -hX/R           , 0  ,  2R - 3 t R'     , -hY/R           , 0
    row 2 (dT)  :  h/R    , -R X' + h Y/R   , 0  ,  -3 h t / R      , -R Y' - h X/R   , 0
    row 3 (dZ)  :  0      , 0               , Y  ,  0               , 0               , -X
    row 4 (dR') :  -mu/R^2 , h X'/R         , 0  ,  -R' + 3 mu t/R^2 , h Y'/R          , 0
    row 5 (dT') :  0      , -R' X' - R X''  , 0  ,  -h/R            , -R' Y' - R Y''  , 0
    row 6 (dZ') :  0      , 0               , Y' ,  0               , 0               , -X'

R' = dR/dt. The inverse is printed as N^-1 = A Phi N^T Phi^T "as before" (eq. 4); the same diagonal A of (2a) is used (COMPUTED: this is correct, section 5). The note claims that in the second set of formulas X, Y and their derivatives may be components "referred to any pair of rectangular axes fixed in the plane of the orbit", because in the product N N0^-1 they occur only in rotationally invariant combinations (X X0 + Y Y0, X Y0' - Y X0'), and that these formulas therefore "can lead to no difficulties when the eccentricity is very small" (the perihelion direction, which defines the axes of the first form, is ill defined as e tends to 0).

### 3.3 Where the generality comes from, and its one limit

- Valid for elliptic and hyperbolic orbits, as a closed form in the state and the elapsed time with no transcendental functions and no Kepler solver (the state at t must be known, from any propagator).
- A contains the semimajor axis a, so the exactly parabolic case (a infinite) is not covered by the formulas as printed; the note's claim "any type of orbit" is thus true for a finite a of either sign. COMPUTED: the matrix degrades smoothly as a grows (section 5, item C).
- The vector products R . R and R x R0' structure of the N N0^-1 product, also used in the note's own remark, are the same invariants as in Deprit and Deprit-Bartholome.

## 4. Comparison with Deprit & Deprit-Bartholome 1968 and with `shepperd_stm`

| | Danby 1965 | Deprit & Deprit-Bartholome 1968 | `core/kepler_stm.shepperd_stm` |
|---|---|---|---|
| Dimension | 3D (6 x 6) | planar only (4 x 4); 3D announced as future work | 3D (6 x 6) |
| Frame | orbital axes (perihelion, in-plane normal, angular momentum) for form 1; any in-plane axes for form 2 (radial/transverse/normal components) | any inertial Cartesian frame (Tables I, II), plus the orbital and intrinsic frames by constant rotations | inertial Cartesian |
| Variables | Cartesian state, t, h, mu, a (a only in the inverse) | Cartesian state, t, mu, and the integrals H, G, P, Q | universal variables (Stumpff functions) |
| Method | element-derivative matrix M (Eckert-Brouwer) times the inverse of its value at t0 | Jacobi last multiplier: adjoint solutions from integrals plus one quadrature; B = A^-1 given explicitly | closed-form universal-variable STM with Stumpff c-functions |
| Derivation shown | none ("the algebra is straightforward") | full | per the module docstring, standard Battin/Shepperd form |
| Singular cases | a infinite (parabolic); form 1 at e = 0 (no perihelion) | G = 0 (rectilinear); H = 0 (parabolic); circular excluded in the derivation | none claimed (multi-revolution, backward, parabolic handled) |
| Printed slips found | none | three (b44, two intrinsic-frame entries) | not applicable |

Both papers make the matrizant a function of the current state and the elapsed time with no solver of Kepler's equation, but in different variables: Danby's M uses explicit Cartesian derivatives of the orbital elements; Deprit's A is built from the integrals. The Danby note is the shorter route to a 3D matrix; Deprit's gives the inverse B in closed form with no a. Deprit's digest says the 3D matrizant was future work; Danby 1964/1965 is the 3D matrizant that predates that remark. Goodyear's universal-time matrizant (which Danby cites) is the ancestor of the universal-variable family to which `shepperd_stm` belongs (INFERRED from the module docstring and from the note; Goodyear 1965 is not held).

## 5. Numerical checks (COMPUTED, scratch only)

Units mu = 1; states given in an inertial frame, rotated into the orbital axes (perihelion from the eccentricity vector, Z along h) for form 1, with the product R = M(t) M^-1(t0) rotated back to inertial axes; the reference is a central finite difference (step 1e-6) of `cyclerfinder.core.kepler.propagate` and the 6 x 6 matrix of `cyclerfinder.core.kepler_stm.shepperd_stm`.

A. Form 1 (eq. 2, 2a), six states (three genuinely 3D with inclined velocity, one retrograde planar, one planar hyperbola of e = 1.24, one near-circular e = 0.02 with a small out-of-plane velocity, plus an e = 0.53 and e = 0.56 ellipse and a hyperbola e = 1.51), dt = 0.7, 2.5 and 9.0 (several revolutions), and two time origins t0 = 0 and 3: relative residual against the finite difference between 6.3e-11 and 4.2e-10 (the size of the finite-difference noise); the two time origins give identical results to the digits printed; against `shepperd_stm` the absolute difference is at most 5.0e-14. So the transcription of M, A, Phi is right and the printed matrix needed no correction.

B. Form 2 (eq. 3, 4): the same cases and the same dt, with the in-plane axes at the perihelion and again rotated by 0.8 rad about Z (the claimed axis independence), using the diagonal A of (2a) in N^-1 = A Phi N^T Phi^T and also the numerical inverse of N0, the components (dR, dT, dZ) and (dR', dT', dZ') rotated to the inertial frame: both give the same residuals as A (6e-11 to 4e-10 relative). So the axis-independence claim holds, N^-1 = A Phi N^T Phi^T holds with the same A, and the printed N has no slip.

C. Near-parabolic hyperbolas (inclined orbit, energy chosen so that e - 1 = 4e-3 to 4e-7 in the orbit, a from -2.5e2 to -2.5e6): against `shepperd_stm` the absolute difference grows with |a|: 5.2e-12, 3.6e-11, 5.0e-9 respectively. This is the loss of significance expected from the a factors of A (INFERRED cause: the factors a/mu in A multiply cancelling terms), so the formulas are usable but degrade as 1/|1 - e|; the exactly parabolic case cannot be evaluated (a infinite). Observation, not investigated: for e - 1 about 1e-8 the finite-difference propagation through `core/kepler.py` raised `KeplerConvergenceError` (chi = -4.4e5); whether that is the finite-difference offset flipping the sign of the energy at the elliptic/hyperbolic boundary or a convergence issue in the universal-variable Newton iteration was not examined.

Not tested: the first element set with dG (the intermediate), the Eckert-Brouwer combinations (HX + KX') = a(YX'/h - 1) and (HY + KY') = -aXX'/h (the note states them, not checked), exactly circular orbits in form 1 (undefined perihelion by construction; form 2 with arbitrary axes was tested only to e = 0.02).

## 6. Printed numbers usable as sourced tests

No numbers: the note has no example. Testable structure: R = M M0^-1 equals the Jacobian of Kepler propagation (form 1) and, rotated to radial/transverse/normal components, N N0^-1 (form 2) equals the same; M^-1 = A Phi M^T Phi^T (that is, A Phi M^T Phi^T M = I; as in the check A), and the axis-independence of form 2. These are the same checks as for Deprit and Deprit-Bartholome and add a 3D case.

## 7. Techniques applicable to the project's problems

**#928 (analytic Kepler STM as the P = 0 control for the KS propagator's transition matrix).** The project already has `core/kepler_stm.py` (Shepperd, 3D); Deprit's closed form is an independent planar second analytic route, and Danby's two forms give an independent 3D one, with 3D controls for a Levi-Civita/KS code that is planar-only at first (the planar block of the 6 x 6 is the control, as in the Deprit digest). Useful points: (i) Form 2 (radial/transverse/normal components) uses only R, R', h, mu, t and the orbit-plane Cartesian state, so it is the easiest to compare with a regularised propagator's output in polar-like variables; (ii) because it is a function of the state at the two times and the elapsed time t only, the physical time at the end of a KS integration (read off from the time element) is the one input to supply; (iii) Danby's M is a product of an element-derivative matrix and its inverse, which is also the structure of a regularised-element STM, so the KS STM in the P = 0 limit factorises the same way (INFERRED, not tested). Constraint: neither paper supplies the radial-fall (G = 0) control that the #928 plan lists, since both exclude rectilinear motion (Danby's h appears in denominators of A; Deprit excludes G = 0 explicitly).

**#899 (the Perko eq. 50 constant at C = -1, V1 = 2).** Perko's a_ij are entries of the matrizant composed with a rotation (Deprit digest, section 5); Danby's form 2 is the natural way to read them in Perko's frame: the (i, j) axes are fixed in the orbit plane (Danby shows the axes may be any in-plane pair), the components dR, dT, dR', dT' are position and velocity components resolved along fixed in-plane directions (the velocity components are resolved displacements, matching the digest's reading of Perko's "j-velocity at t1"), and N N0^-1 evaluated between the perpendicular crossing at t = 0 and the tangency at t1 gives the four needed entries once the axes are rotated to put i along V1. Cross-check target: the entries computed in the Deprit digest section 7 for the retrograde unit circle (a21 = -2, a24 = 2, a41 = -3.712389, a44 = 2.712389 under the stated conventions). The two closed forms must agree (the 3D 6 x 6 of `shepperd_stm` already agrees with Deprit's to 4e-14). This does not settle Perko's frame conventions (the open follow-up 1 of the Deprit digest). One caveat: for the circular type-3 orbit the perihelion axis of form 1 is undefined, so use form 2 or `shepperd_stm`.

## 8. Source honesty and limits

- No derivation, example or numerical value is printed; every number above except the matrices is my computation.
- The formulas were checked numerically and need no correction; this is the only reliable claim about their validity here, and a re-derivation was not attempted.
- Not held (referenced): Danby 1964 (AIAA J. 2:16 and 2:13), Goodyear and Guseman 1964 (IBM/RTCC 12-008), Herrick 1960, Eckert and Brouwer 1937, Brouwer and Clemence 1961 (reference 6: the book is a standard text; whether the corpus holds it was not checked).
- The scan carries two unrelated items: the first page ends a plasma-accelerator note (Hall accelerator, equations (11) to (15), Fig. 3, references 1 to 7 of that note) that has no bearing on this project.

## 9. Follow-ups

1. Cross-check Deprit's section 7 entries against Danby's form 2 for the retrograde unit circle (planar block of the 6 x 6; same numbers expected).
2. If an analytic Kepler STM is wanted in `src/` for inclined or near-circular orbits, `shepperd_stm` already covers it; Danby's form 2 is a candidate only as an independent test fixture (no singularity at e = 0).
3. Investigate whether `core/kepler.propagate` fails near the parabolic boundary (e - 1 about 1e-8, chi = -4.4e5, `KeplerConvergenceError`) or whether this is only the finite-difference offset crossing the boundary; relevant to #928's near-collision work and to any near-parabolic arc.
4. Candidate acquisitions if a derivation is wanted: Danby 1964 (AIAA J. 2:16), Goodyear and Guseman 1964 (IBM/RTCC 12-008), Goodyear 1965 (Astron. J. 70:189). Check CORPUS_INDEX first; none is claimed held.
