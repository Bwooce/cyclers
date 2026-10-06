# Digest: Lacomba and Llibre 1988, "Transversal Ejection-Collision Orbits for the Restricted Problem and the Hill's Problem with Applications" (#960 batch 32)

Ernesto A. Lacomba (Universidad Autonoma Metropolitana, Iztapalapa, Mexico) and Jaume Llibre (Universitat Autonoma de Barcelona),
Journal of Differential Equations 74:69-85 (1988), received 29 January 1987, revised 29 September 1987.
DOI 10.1016/0022-0396(88)90019-8. 17 pp. The title matches the page image.
- Source scan: `84c5547a-lacomba1988.pdf`, md5 b71fe88d4a0f103673ec689976c9f0a0. It has a good publisher text layer.
- Proposed corpus filename:
  `lacomba-llibre-1988-transversal-ejection-collision-orbits-restricted-hill-problem-jde-74-69-doi-10.1016-0022-0396(88)90019-8.pdf`
- How I read it: text layer for the structure; every theorem, the key equations, the numbers of Part B and the integral
  evaluations on pp.79-80 were read on 110 to 300 dpi page images (pp.1, 4, 7, 11, 12, 16, 17). Checks are in
  `check_lacomba.py` with output `check_lacomba.out`.

## 0. Verdict

**This is the proof that the planar circular restricted problem and Hill's problem have transversal ejection-collision (EC)
orbits, and the source of the "no extendable regular integral" corollary.** It is the paper behind several lines of the held
digests (Llibre 1982, Rodriguez del Rio thesis, which cite it as "no C^1-extendable regular integrals").
- Part A (restricted problem): an analytic proof, Melnikov-like, for **mu small and C large**.
- Part B (Hill's problem): **numerical only**, at the single level C = 5. The paper says so (p.70: "we can only prove
  numerically"; Assumption B.2 is "proved numerically later on").
- Proposed catalogue implication: none. It is a theory paper. It gives one numerical EC orbit of Hill's problem (sec. 3)
  that is a usable positive control for an EC solver. I reproduced it.

## 1. Setting (pp.69-72)

- Planar circular restricted problem, rotating frame of frequency 1. Larger primary `m1 = 1 - mu` at the origin, smaller
  `m2 = mu` at `e2 = (-1, 0)`, `mu` in [0, 1/2]. Hamiltonian (0.1):
  `H = |p|^2/2 + q2 p1 - q1 p2 - 1/|q| + mu(1/|q| - 1/|q - e2| - p2)`. The `- mu p2` is correct for rotation about the barycentre.
- Hill's problem (0.2): `H = |p|^2/2 + q2 p1 - q1 p2 - 1/|q| - q1^2 + q2^2/2`. Jacobi constant `C = -2H`. In these units
  L1, L2 are at `3^(-1/3)` from the origin (p.82), so lengths are in Hill units `(mu/3)^(1/3)` times the primary separation.
- Regularisation: polar coordinates, then McGehee variables `v = r^(1/2) rdot`, `u = r^(3/2) thetadot`, time `dt = r^(3/2) dxi`.
  The collision manifold at `r = 0` is the torus `(u^2 + v^2)/2 = 1 - mu` (restricted) or `= 1` (Hill). The flow on it is the
  Kepler flow: two circles of equilibria `u = 0, v = +-[2(1 - mu)]^(1/2)` (Fig. 1). The usual McGehee variables fail for the
  restricted problem because the potential is not homogeneous. The method of Llibre 1982 is used instead.
- Regions: C above the L2 level so that the Hill region has a bounded component with only the larger primary (restricted),
  or `C > 3^(4/3) = 4.326749` (Hill; my arithmetic).

## 2. Results, as printed

- **Proposition A.1.** For mu in [0, 1/2] and C > C_2(mu), the extended level component is a closed solid torus whose
  boundary is the collision manifold.
- **Proposition A.2 and A.3.** At mu = 0 all ejection orbits are collision orbits: the circle `r = 2/C, v = 0, u = -(2/C)^(3/2)`
  (zero angular momentum, apocentre at `2/C`). For small mu > 0 the two curves `gamma^u_C(mu)`, `gamma^s_C(mu)` stay circles and,
  by the symmetry S, meet at `theta = 0` and `theta = pi`: at least two EC orbits near `(2/C, 0, 0, -(2/C)^(3/2))` and
  `(2/C, 0, pi, -(2/C)^(3/2))`.
- **Theorem A.4.** The restricted problem has transversal EC orbits for mu small enough and C large enough. The proof
  computes the first variation along the mu = 0 heteroclinic orbit (an explicit cosh and tanh solution, p.76). It reduces
  to `I = 24 C^(-5/2) (26/45 - 49 pi/192) + O(C^(-7/2))`, which is nonzero.
  - I evaluated the four integrals on p.79 and the combination: `pi/32`, `1/5`, `(3 pi - 7)/9`, `3 pi/64` all agree to 10
    digits, and the combination is -0.2239828474 = 26/45 - 49 pi/192. The sign is negative, so I is nonzero.
  - "By analyticity, Theorem A.4 holds for almost all mu in (0, 1/2] and C > C_2(mu)" (p.80).
  - Remark: **Pinol** (thesis, 1987) shows the integrand `cosh s (sinh s - cosh s) f(s; 0)` has one sign. That would remove the
    "C large" condition (three-component Hill region suffices). It is a private thesis citation; not proved here.
- **Theorem A.5.** For mu > 0 small the restricted problem has no extendable regular integral (C^1, constant on the
  collision manifold, regular on levels). Proof: a transversal EC orbit would make `W^s + W^u` three-dimensional inside a
  codimension-1 level set of the integral, which is a contradiction. mu = 0 has one (angular momentum).
- **Theorem B.5 and B.6.** Hill's problem has transversal EC orbits (numerical, C = 5) and no extendable regular integrals.

## 3. Numbers of Part B (pp.83-85) and my check

- At C = 5 the EC orbit starts on the q1 axis with `p0 = (r0, v, theta, u) = (0.429043..., 0, 0, -0.302848...)` (p.84), by the symmetry
  also `p_pi` with theta = pi. `r0` is the first-return root of a Poincare map on the annulus `v = 0` (Fig. 6, found from the computed return map,
  6 digits). The transversality derivative of the curve `r = p(theta)` at theta = 0 is **0.026...** (p.85).
- **Check 1, algebra.** From the Jacobi relation (B.2) `u^2 + v^2 - 2 = -C r + 3 r^3 cos^2(theta)` at v = 0, theta = 0, r = 0.429043,
  C = 5 I get `u0 = -0.302848` (matches the printed digits; residual 1e-16).
- **Check 2, integration.** I integrated equations (B.1) with DOP853, rtol 1e-13, from `p0` and from `r0 +- 1e-4`.
  From `p0` the orbit reaches `r = 1e-4` at `xi = 6.993` with `v = -1.4140` (limit `-sqrt 2 = -1.4142`), `u = -2.1e-4`, `theta = -0.348`.
  The two neighbours end with `u = -4.1e-2` and `+4.0e-2`. So `u` changes sign across `r0` and `r0 = 0.429043` is the collision root,
  confirming that the printed `p0` is an EC orbit (apocentre 0.429 Hill units).
- The value 0.026 was **not** reproduced. The paper does not define the scaling of `p(theta)` closely enough
  for me to recompute it without guessing.
- Printed-constant note: the intro says `C = -2H` with H of (0.1), but the L4/L5 level is stated as 3 for every mu and the upper bound
  4.25 (p.72). With `C = -2H` the L4 level is `3 - mu(1 - mu)` (my algebra: 3, 2.91, 2.84 ... 2.75 for mu = 0 ... 0.5).
  The quoted 3 and 4.25 are Szebehely's convention (`+ mu(1 - mu)`). Be careful if you compare level values across sources.
  Hill's problem and the large-C results are not affected. The definition of `I_C` on p.72 also writes "= C" where H is meant.

## 4. What it gives the project (#899)

- **Existence and transversality of EC orbits** at large C (Theorem A.4), the base for Chenciner-Llibre 1988 (next digest),
  which removes the small-mu condition and adds the torus consequence.
- **A reproducible positive control.** Hill's problem, C = 5, `r0 = 0.429043`, apocentre on the q1 axis, collides after xi = 6.99
  in McGehee time. Use it to test any Hill-problem EC solver or the project's regularised integrator. Judge by the same
  criterion as the paper: `u -> 0`, `v -> -sqrt 2` as `r -> 0`.
- It bounds scope: Part A is perturbative (mu small, C large). No value of `mu_0(C)` or `C_0` is given.

## 5. Citation mining

- Llibre 1982 (Celest. Mech. 28:83-105, restricted problem, small mu): HELD (`llibre-1982-restricted-three-body-problem-mass-parameter-small-...`).
- Devaney 1981 "Singularities in classical mechanical systems" (in Katok, Birkhauser, pp.221-333): not held. The held Devaney file
  is a different paper (CMP 80:465). Wanted-list row 36.
- Llibre and Simo 1980 (Math. Ann. 248:153-184, oscillatory solutions): not held, not on the wanted list. Llibre 1984 (JDE 54:221,
  extendable integrals of the n-body problem): not held, not listed.
- McGehee 1974 (Invent. Math. 27:191): not held. Moser 1973 "Stable and random motions" (book): not held.
- Pinol 1987 thesis (UAB): not found; listed with Llibre-Pinol collision orbits, wanted-list row 27.
- Siegel 1936 (algebraic integrals): not held. Gravalos 1941 thesis: not held. Szebehely 1967 "Theory of Orbits": HELD.
- Chenciner and Llibre 1988 (cited as "to appear"): in this batch, wanted-list row 36 (shared with this paper).

*Filed as `cyclers_pdf/papers/lacomba-llibre-1988-transversal-ejection-collision-orbits-restricted-hill-problem-jde-74-69-doi-10.1016-0022-0396-88-90019-8.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
