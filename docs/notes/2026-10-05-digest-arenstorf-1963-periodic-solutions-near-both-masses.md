# Digest: Arenstorf 1963, "Existence of Periodic Solutions Passing Near Both Masses of the Restricted Three-Body Problem" (#960)

R. F. Arenstorf (NASA Marshall Space Flight Center, Computation Division), "Existence of Periodic Solutions
Passing Near Both Masses of the Restricted Three-Body Problem", AIAA Journal 1(1):238-240 (January 1963),
Technical Notes and Comments, DOI 10.2514/3.1516.
- Crossref-confirmed 2026-10-05: 1(1), pp. 238-240.
- Presented at the ARS 17th Annual Meeting, Los Angeles, 13-18 Nov 1962.
- Filed as `cyclers_pdf/papers/arenstorf-1963-existence-periodic-solutions-passing-near-both-masses-restricted-three-body-aiaa-j-1-238-doi-10.2514-3.1516.pdf`.
  3 pages, text layer, md5 66f3c6970ba7e72300bfb984223827d6.
- Page 238 opens with the tail of an unrelated MHD note, and page 240 continues with an unrelated note
  (Porter & Hatfield). Those parts were ignored.
- Pages 239-240 were read as page images.
- The note has no reference list of its own. The references printed on p.238 belong to the MHD article.

Evidence tags: READ, COMPUTED, INFERRED.

## 0. Verdict

This is an announcement with a proof outline. The full proof is "to appear elsewhere"; it appeared as Amer.
J. Math. 85:27-35 (1963), doi 10.2307/2373181, which is not held. The result concerns periodic orbits of the
second KIND (near rotating Kepler ellipses), not second species.

- For any integers m, k with a^(3/2) = m/k, and for all e in a closed interval that avoids finitely many
  collision values, there is, for fixed sufficiently small mu > 0, a one-parameter family of
  synodically-closed symmetric solutions near the rotating ellipse.
- Choosing a and e by Eq. (9) makes them pass near BOTH Earth and Moon.
- Earth-Moon mu is asserted, not proved. It is illustrated numerically at mu = 1/82.

The note also proposes, in 1963, an Earth-Moon "ferry vehicle perpetually on such a path". That is an
Earth-Moon cycler concept and should be recorded as prior art.

## 1. Statement (READ p.238-239)

- Abstract: "There exist in the restricted three-body problem with small mass ratio one-parametric families
  of synodically closed solution curves, which are near rotating Keplerian ellipses with arbitrary rational
  sidereal frequencies and appropriate positive eccentricities. By suitable selection of the parameter
  values, these periodic solutions can be made to come close to both attracting bodies." (p.238)
- Generating solution x*(t) at mu = 0: a rotating ellipse with T0 = 2 pi |a|^(3/2), and a^(3/2) = m/k with m
  natural and k a nonzero integer, positive for direct and negative for retrograde. The synodic period is
  T* = 2 pi m = |k| T0. "it closes after k - m revolutions around the origin" (p.238).
- Restrictions on e (p.239): "For given a = (m/k)^(2/3), there are at most finitely many e in 0 < e < 1 with
  e = (1 - a^-3)^(1/2) or x*(t) = 1 at least once in 0 <= t <= T*. In the latter case, P collides with M.
  Such e have to be omitted."
- Result (p.239): "for every closed e interval I not containing such exceptional values and for fixed
  sufficiently small mu > 0, there exists a family of periodic solutions ... depending continuously upon the
  parameter e in I, which transfers into the corresponding family of solutions x*(t) for mu -> 0".
- Proof route (p.239):
  - Symmetric periodicity: two perpendicular crossings of the real axis, at t = 0 and t = T/2.
  - The implicit system [4] is solved for T and eta through a non-zero functional determinant
    D* = 3 m pi eta* e[(-1)^k - e]/(c* - a^-1) c*^2 (Eq. 6).
  - Variables H (the Jacobi integral), F, U, V and c are introduced (Eqs. 7-8).
  - eta -> 0 needs a modified condition and "can occur only when 2a <= 1".

## 2. Near both masses, and Earth-Moon mu (READ p.239-240)

- Eq. (9): choose a = (m/k)^(2/3) and e with "delta/2 < a(1 - e) < 2 delta, 1 < a(1 + e) < 1 + delta", "for
  instance, with suitable small delta > 0".
- "this can be achieved for x*(t) and thus also for the corresponding x(t), which is reasonably near x*(t)
  still for mu = 1/82 ~ the value of mu for the case E = Earth, M = moon". This is an assertion; the proof is
  for sufficiently small mu.
- "m = 2, k = +5 gives for a^-1 - 1 < e < 1 a highly interesting family of solutions passing close to E and
  M, with no exceptional e values in this range" (p.240).
- When k - m is odd, solutions passing near E and M can also be started at the perigee point.
- "m = 1, k = +2 yields promising flight paths for solar probes returning to Earth after a little less than
  1 yr" (p.240).
- Figs. 1-4, "calculated for the case mu = 1/82", are pictures only, with no numbers:
  - Fig. 1: rotating-frame closed path, m = 1, k = 2.
  - Fig. 2: rotating-frame closed path, m = 2, k = 5, with five numbered loops.
  - Fig. 3: inertial view of the Fig. 1 orbit, showing "capture and rejection of P by M".
  - Fig. 4: Earth-centred view of the Fig. 2 orbit.
- Prior-art passage (p.240): "a radiation protected heavy Earth-moon ferry vehicle perpetually on such a
  path, to be supplied or boarded by passengers after rendezvous with much smaller crafts near E or M ...
  Thus tourist trips to the moon and back become practical and more economical."

## 3. Gate answers for `#948` R4

- **Theorem class:** second-kind, symmetric, planar, circular orbits for small mu (with mu_1 unspecified).
  They are near Kepler ellipses whose collision values of e are excluded. The Moon passage distance stays
  of order delta, which is fixed. It does not tend to 0 with mu, so these are not second-species orbits.
- **Earth-Moon mu:** not covered by the theorem. The note asserts it and illustrates it at mu = 1/82.
- **Which orbits can pass near both masses (COMPUTED):**
  - General bound: any ellipse with a small perigee r_p that reaches the Moon distance r = 1 needs
    a >= (1 + r_p)/2 > 1/2. This holds for both of Arenstorf's constructions: the apogee-at-Moon
    construction of Eq. (9), and the p.240 variant for k - m odd that comes close to M 'at a later time
    (not at apogee)'.
  - The Eq. (9) instance:
  - Both inequalities together need 1 < a(1 + e) and a(1 - e) < 2 delta. Adding them, 2a lies between
    1 + delta/2 and 1 + 3 delta, so a is just above 1/2.
  - Hence m/k is just above 2^(-3/2) = 0.35355.
  - The paper's example m/k = 2/5 = 0.4 gives a = 0.5429.
- **Relation to Genova & Aldrin 2015:** they wrote "such a cycler was not shown to exist in the restricted
  three-body problem: Arenstorf shows a 3:1 resonance orbit but without the required close Earth passes".
  - A 3:1 orbit (three spacecraft revolutions per lunar month) has m/k = 1/3, a = 0.48075 and 2a = 0.9615,
    which is less than 1 (COMPUTED).
  - At mu = 0 such an ellipse cannot reach the Moon's orbit, so Eq. (9) can never be met. Arenstorf's
    theorem cannot supply a 3:1 orbit passing near both bodies.
  - Convention checked: Genova & Aldrin's 3:1 is three spacecraft revolutions per lunar month. Their
    lunar encounter is every 26 days and their Earth perigees are every 7-10 days. Their own text has the
    CR3BP apogee 'below lunar distance' (`docs/notes/2026-06-10-genova-aldrin-2015-mining.md` sec. 2), so
    m/k = 1/3 is right.
  - This is consistent with Genova & Aldrin's statement.
  - The Arenstorf 3:1 figure they cite is not in this note. It is probably in the Proc. XIV IAC 1963 paper
    or NASA TN D-1859 (INFERRED; neither is held).
  - The same a > 1/2 bound appears in Hitzl & Henon's eq. 47 for second-species C_ij families: a > 1/2
    requires j < 2 sqrt(2) i.
- **Consequence for R4:** an Earth-grazing 3:1 member at Earth-Moon mu would have to come from second-species
  (near-collision) dynamics, where the Moon's gravity raises the apogee. That is Bruno's e = 1 arcs and
  Gomez-Olle double collisions, not Arenstorf's second-kind continuation. R4's framing is unchanged.

## 4. Relation to the held corpus

- It is cited by Perko 1974 (ref. [2]) and Hitzl 1977 (via ref. 2, the Amer. J. Math. version).
- It is the origin of the "Arenstorf orbits" that Hitzl 1977 Figs. 1-4 reproduce at mu = 1/82.30, with the
  Davidson and Causey computations.
- Broucke 1968 (held) is the numerical Earth-Moon family source of the same era.

## 5. Positive controls

There are no numbers. The figures are qualitative: m/k = 1/2 and 2/5 at mu = 1/82.

## 6. Citation mining (policy step 4)

The note itself has no reference list. Cited works by implication, not held:
1. Arenstorf, R. F. (1963), "Periodic Solutions of the Restricted Three Body Problem Representing Analytic
   Continuations of Keplerian Elliptic Motions", Amer. J. Math. 85:27-35, doi 10.2307/2373181 (CONFIRMED).
   The full proof. It is also NASA TN D-1859 (per Hitzl 1977 ref. 2).
2. Arenstorf, R. F. (1963), "Periodic trajectories passing near both masses of the restricted three-body
   problem", Proc. XIV IAF Congress (Paris 1963), Springer. No DOI found (UNCONFIRMED). Most likely home of
   the 3:1 figure that Genova & Aldrin cite.
