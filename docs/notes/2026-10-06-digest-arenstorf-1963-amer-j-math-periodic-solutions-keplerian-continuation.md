# Digest: Arenstorf 1963, "Periodic Solutions of the Restricted Three Body Problem Representing Analytic Continuations of Keplerian Elliptic Motions" (#960 batch 25)

R. F. Arenstorf (NASA Marshall Space Flight Center and University of Alabama), American Journal of
Mathematics 85(1):27-35 (January 1963), doi 10.2307/2373181. Received 13 September 1962. The wanted
list also records it as NASA TN D-1859.
- Filed as `cyclers_pdf/papers/arenstorf-1963-periodic-solutions-restricted-three-body-analytic-continuations-keplerian-elliptic-motions-amer-j-math-85-27-doi-10.2307-2373181.pdf`.
  JSTOR copy, 10 pages (a JSTOR cover plus 9 pages), text layer, md5 2c531bdc1db85ad905a0910e00c35803.
  Supplied by the owner.
- I read the whole paper from the text layer. The equations were checked for structure only: the layer
  garbles symbols, so read any formula on the page image before quoting it.
- This is the full existence proof behind the held AIAA J 1:238 note
  (`arenstorf-1963-existence-periodic-solutions-passing-near-both-masses-...`). The owner's re-sent
  `Arenstorf63b.pdf` is another download of that held note and was NOT filed.

## 0. Verdict: what the theorem says (citable support for "theorem-generic" claims)

**Setting.** The planar circular restricted problem in the rotating frame, masses 1 - mu and mu,
0 <= mu <= 1 (eq. 1).

**Generating orbits.** At mu = 0, a Kepler ellipse about the origin (eq. 2):
- semi-major axis a, eccentricity 0 < eps < 1;
- started at apocentre, z(0) = a(1 + eps);
- sidereal period T0 = 2 pi a^(3/2) commensurable with 2 pi: a = (m/k)^(2/3), with k and m coprime,
  m > 0, and k > 0 for direct and k < 0 for retrograde motion.

In the rotating frame this is periodic with synodic period T* = 2 pi m, closing after |k - m|
revolutions.

**Theorem** (stated on pp.27-28 and proved on pp.29-32). Fix a = (m/k)^(2/3).
- **Exceptional values.** At most finitely many eps in (0, 1) are exceptional (eq. 23):
  - (i) eps = sqrt(1 - a^-3), for a > 1 and direct motion (c* > 0);
  - (ii) eps such that the generating orbit **collides with the second mass**: x*(t) = 1 for some t
    in 0 <= t <= m pi, the half period, which suffices by symmetry (eq. 23, read on the p.33 image).
    The paper proves that only finitely many eps qualify, and none when 2a <= 1 (p.34).
- **Existence.** For every closed eps-interval I free of exceptional values, there is mu* > 0 such
  that for every 0 <= mu < mu* there is a family of **symmetric** periodic solutions (symmetric about
  the x1-axis, the line of the primaries). The family depends analytically on eps in I and on mu, and
  reduces to the Kepler ellipse at mu = 0.
- **Properties.** The synodic period T(mu, eps) and the Jacobi constant are holomorphic in eps and mu,
  and depend on eps.
- **Second family.** Starting at pericentre (eps -> -eps) gives a second family. It is distinct as a
  set of curves when k - m is odd (p.34).

**What it covers:**
- These are Poincaré's **periodic solutions of the second kind**: continuations of elliptic (eps > 0)
  Kepler orbits. The paper says earlier "proofs" (Poincaré, Schwarzschild, Charlier) were shown invalid
  by Staeckel and Wintner. The method is Poincaré continuation with Birkhoff's **symmetric** periodicity
  condition (eqs. 4-5, 13) in place of the general one.
- The non-degeneracy determinant, eq. (21), read on the p.32 image, is
  D* = 3 eps (-1)^m (1 - eps(-1)^k)^2 m pi / ((c* - a^-1) c*^3), which is nonzero for 0 < eps < 1 and
  c* != a^-1. Exceptional value (i) is exactly a c* = 1 (equivalently 2h* + 3c* = 0). There the
  generating Jacobi constant j* = -(1/2) w^(2/3) - sign(c*) w^(-1/3) sqrt(1 - eps^2) (with
  w = |a^(-3/2)|) is at its maximum, direct case only (p.33).
- **Isoenergetic continuation is not obtained:** for fixed Jacobi constant the relevant determinants
  vanish (pp.32-33). The paper cites Wintner: for small eps, isoperiodic solutions do not exist. The
  continuation is in mu with T and h adjusting.
- First-kind solutions (eps = 0) also follow with the same method, except at a few ratios.

**What it does NOT cover:**
- **Orbits that pass THROUGH or arbitrarily close to the second mass.** Collision orbits are excluded,
  and mu* depends on the eps-interval, shrinking as I approaches a collision value. The orbits "passing
  repeatedly near both masses" (the Arenstorf / Apollo-type figure-eight family) are mentioned as
  numerical motivation only (p.29). "The calculations indicate their existence for values of mu at least
  as large as" Earth-Moon. **That claim is numerical, not part of the theorem.** For near-collision
  (second-species) orbits, use Breakwell-Perko, Perko 1974 and Hitzl 1977 (held), not this paper.
- Large mu: the theorem is for small mu only, with no explicit bound on mu*.
- Non-symmetric periodic orbits, the spatial problem, and the elliptic restricted problem.

## 1. Use in the open routes

- `#944` X2, `#946` X3, `#948` R4: cite Arenstorf 1963 for the claim that "every non-colliding,
  commensurable rotating Kepler ellipse continues to a symmetric periodic orbit family for small mu
  (except finitely many eccentricities)".
  - It is a small-mu, non-collision, symmetric-orbit existence result.
  - It does NOT guarantee existence at a specific planetary mu (e.g. Earth-Moon 0.01215, or Sun-planet
    values) unless that mu lies below the unquantified mu*.
  - It does NOT cover close-approach orbits.
  - A "theorem-generic" claim at a given mu needs a numerical continuation as well.

## 2. Citation mining (refs 1-10)

- [1] Birkhoff 1915 (Rend. Circ. Mat. Palermo 39:265-334): not held.
- [2] Charlier, Mechanik des Himmels: not held.
- [3] Koopman 1927 (Trans. AMS 29:287-331): not held.
- [4] Moser 1953 (Math. Ann. 126:325-335): not held.
- [5] Poincaré, Méthodes nouvelles: not held.
- [6] Schwarzschild; [7] Siegel; [8] Staeckel; [9]-[10] Wintner: not held.

All are classical existence theory. One low-priority row is added to the wanted list (Birkhoff 1915,
Moser 1953 and Koopman 1927 as the closest to symmetric-continuation and many-revolution periodic
orbits). The rest are textbooks or historical.
