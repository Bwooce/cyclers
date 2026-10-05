# Digest: Perko 1974, "Periodic Orbits in the Restricted Three-Body Problem: Existence and Asymptotic Approximation" (#960)

L. M. Perko (Northern Arizona University), "Periodic Orbits in the Restricted Three-Body Problem: Existence
and Asymptotic Approximation", SIAM Journal on Applied Mathematics **27(1)**:200-237, July 1974, DOI
10.1137/0127016. Received 12 Dec 1972, revised 13 Jul 1973.

The issue is 1, not 2:
- The PDF title page reads "Vol. 27, No. 1, July 1974".
- Crossref (2026-10-05) gives volume 27, issue 1, pp. 200-237, published 1974-07.
- The `#938` note's "27(2)" was a slip, corrected in sec. 5 item 3 by `#960`.

Filed as `cyclers_pdf/papers/perko-1974-periodic-orbits-restricted-three-body-existence-asymptotic-approximation-siam-j-appl-math-27-200-doi-10.1137-0127016.pdf`.
- 38 pages, text layer (scan with OCR), md5 76debe07c1d115322c20c6c230d5806f.
- Journal page = PDF page + 199.
- I read all pages from the text layer, which is heavily garbled in the formulas. Theorems 1-4 (pp.221,
  224, 230, 232), Proposition 4 (p.225), the Appendix opening (p.233) and the references (p.237) were read
  as page images.

Evidence tags: READ = printed (journal page). INFERRED = my reading.

## 0. Verdict

This is the second-species existence theorem that every held Perko paper builds on. It covers symmetric
planar circular-problem orbits that pass near the smaller primary once (type A) or twice (type B) per
period, for all sufficiently small mu > 0. The upper limit mu_1 is not specified and depends on the
generating ellipse.

It supports the label "theorem-generic" only for that class and only in the limit of small mu. It does
not, by itself, cover:
- a fixed physical mass ratio. No explicit mu range is given, so the Earth-Moon value is not covered.
- asymmetric orbits.
- chains with more than the type-A or type-B number of passages.
- Earth-collision (e = 1) arcs.
- the elliptic problem.

## 1. Setting (READ pp.200-205)

- Planar circular restricted problem in Earth-centred inertial coordinates, Eq. (4). mu = mu*/(1 + mu*),
  where mu* is the Moon-to-Earth mass ratio. Unit Moon radius and angular velocity.
- The proof uses the boundary-layer (matched-asymptotic) approximation of Breakwell & Perko 1965 and Perko
  1964, with the error estimates of Perko 1964 and 1967 (summarised in the Appendix, pp.232-237).
- Periodicity condition: two perpendicular crossings of the Earth-Moon line (Birkhoff 1915; Arenstorf
  1963). **Only symmetric orbits are treated** (p.204).
- Generating ellipses (sec. 2):
  - **Type B, ratio m/k:** a0 = (m/k)^(2/3) and x0(t1) = xm(t1). The ellipse meets the Moon again after
    T0 = 2 pi m, with m Moon orbits and k particle orbits in T0 (p.201).
  - **Type A, ratio m/k:** x0(t1) = xm(t1), plus a timing condition F_mk(a0, e0) = 0 (Eq. 3) or its
    variant. These are solvable "for a continuum" of (a0, e0) when dF/de is not 0, for example if
    a0(1 - e0) < 1 < a0(1 + e0) for direct orbits (pp.202-203).
  - For type A, the apogee (k even) or perigee (k odd) lies on the Earth-Moon line (p.203).
  - Arenstorf's Figs. 1-7 of his 1963 IAC paper come from type A. His Figs. 8-9 come from type B with
    ratio 2/4 ("ratio 1/2, order 2") (p.204).

## 2. The theorems (READ, page images)

**Theorem 1 (p.221), type A.** Let m and k be positive integers and (a0, e0) in E0 define a type-A
generating ellipse with apogee or perigee time t0 and x0(t1) = xm(t1). Let a0(1 - e0) < 1 < a0(1 + e0) and
a0(1 - e0^2) not equal to 1.

Then:
- For all but possibly a finite number of such (a0, e0) in any compact subset of E0, with a21 a24 not 0,
- and given k0 > 0 and eps > 0, there is mu_1 > 0 such that for mu in (0, mu_1)
- there exists a one-parameter family of solutions of (4), periodic in rotating coordinates.
- The initial conditions are x(t0) = x0(t0) + dx0 and xdot(t0) = xdot0(t0) + dxdot0, with either:
  - dxdot0 = (0, dv), |dv| <= mu k0, and dx0 given to first order by [mu(K0 - K1) - a24 dv]/a21; or
  - dx0 = (dr, 0) and dxdot0 given by [mu(K0 - K1) - a21 dr]/a24.
- The period in rotating coordinates is T = 2[t*(...) - t0] = T0 + 2 dt1 + O(mu^(2-eps)), where
  |V1| dt1 = a11 dr + a14 dv + (mu/V1^2)|ln mu e1| + mu K2.

Remarks 3-4 (pp.221-222): a0 or e0 may be used as the family parameter instead. As mu -> 0 the orbit tends
to two pieces of the generating ellipse joined at a corner at t1, which is a SINGLE corner at the Moon in
rotating coordinates.

**Theorem 2 (p.224).** Under Theorem 1's hypotheses, for mu in (0, mu_1), denumerably many of these
periodic solutions are periodic in both rotating and inertial coordinates. If (t* - t0)/pi = p/q, the
inertial period is qT = 2 p pi. Example (p.225): the direct type A ratio 1/2 with t* - t0 = 2 pi/3 has
period 3T = 4 pi (Fig. 5).

**Proposition 4 (p.225), type B.** If a0(1 - e0) < 1 < a0(1 + e0), the relative velocities at the two Moon
meetings are equal, V1 = V2, with |V1| > A0 = |sqrt(a0(1 - e0^2)) - 1|. The initial perpendicular crossing
is set at the perilune (Eq. 9).

**Theorem 3 (p.230), type B.** Let a0(1 - e0) < 1 < a0(1 + e0) and a0(1 - e0^2) not equal to 1.

Then:
- For all but possibly finitely many (a0, e0) in any compact subset of F0, with a23 a24 not 0,
- there is mu_1 such that for mu in (0, mu_1) a one-parameter family of solutions exists, periodic in
  rotating coordinates.
- It is parametrised by dV1 or d alpha1 near the perilune: r_p = |Delta1|(e1 - 1)^(1/2)/(e1 + 1)^(1/2),
  with Delta1 = +-mu |tan alpha1|/V1^2 and e1 = sec alpha1.
- T = 2 T0 + 2 dt2 + O(mu^(2-eps)) (p.231).

Remark 7: as mu -> 0 the orbit tends to two intersecting type-B ellipses, with TWO corners at the Moon in
rotating coordinates, at t1 and t2.

**Theorem 4 (p.232).** Under Theorem 3's hypotheses, denumerably many members are also inertially
periodic. The inertial period is of order 1/(mu |ln mu|).

## 3. What the theorems cover and exclude, for `#944` X2, `#946` X3 and `#948` R4

Covered: symmetric, planar, circular-problem second-species orbits whose limit is:
- a single type-A Kepler arc closed at one Moon collision per period, or
- a type-B pair with two collisions,

with crossing ellipses (perigee inside and apogee outside the Moon orbit), for mu small enough.

| Item | In or out of Perko 1974 | Where it is covered instead |
|---|---|---|
| Explicit mass range | Out. mu_1 is existential and depends on (a0, e0, k0, eps). | No numerical mu is proved, so "theorem-generic at the Earth-Moon or Titan mass" needs separate evidence such as continuation from small mu. |
| Exceptional (a0, e0) | Finitely many per compact set; those with a21 a24 = 0 or a23 a24 = 0. | Tangential encounters, a0(1 - e0^2) = 1, are excluded explicitly. |
| Asymmetric orbits | Out. | Bolotin, and Font, Nunes & Simo. Relevant to X2: an asymmetric S-arc chain is NOT "Perko theorem-generic". |
| Chains with three or more distinct Moon passages per period, or several different ellipses | Out. | Gomez & Olle; Bolotin; Perko 1976, 1977, 1981 for O(mu) near-Moon passages. |
| Elliptic restricted problem | Out (circular only). | X3 must cite Bolotin 2005, not Perko. |
| Earth-collision or Earth-grazing arcs | Out. The proof needs a0(1 - e0) > 0 with perigee bounded away from Earth (Lemma 1 proof, p.207), and e0 < 1 strictly. | R4's Bruno e = 1 arcs and Gomez-Olle double-collision orbits are outside Perko 1974. |
| Spatial problem | Out (planar). | |

INFERRED reading for policy: Perko 1974 makes symmetric single- and double-passage second-species members
"theorem-generic as mu -> 0". A specific member at a physical mass is a continuation of a theorem-guaranteed
family. It is not itself a theorem consequence.

## 4. Relation to the held Perko and companion papers

- Perko 1967, SIAM J. Appl. Math. 15:738 (ref. [10]), held and digested: the error-estimation method used
  in the proof.
- Breakwell & Perko 1974, Celest. Mech. 9:437 (ref. [6], "to appear"), held: second-order matching. The
  O(mu^2 ln^2 mu) terms are named in Remark 1 (p.209).
- Perko 1976, 1976b, 1977, 1981, 1981b, held:
  - the later O(mu) near-Moon-passage families (1976, 1977, 1981).
  - the first/second-species bifurcation theory (1981b). Its p.201 bridge statement explicitly combines with
    "Perko 1974" to give a complete theory for families with only first- and second-species orbits.
- Arenstorf 1963, AIAA J 1:238 (ref. [2]), filed by `#960`: the second-KIND (near-ellipse) existence
  counterpart.
- Arenstorf 1963, Amer. J. Math. 85:27 (ref. [1]), not held: the full proof.
- Hitzl & Henon 1977 call these limits "generating orbits by Perko" (Hitzl 1977 p.1410).

## 5. Positive controls

The paper gives no numerical orbit and no table. The only printed geometry is figure-level:
- Figs. 1, 3, 4, 5 sketch type B 1:1, type A 1:2 and the inertially periodic example with period 4 pi.
- The Theorem 2 example (ratio 1/2, t* - t0 = 2 pi/3, inertial period 4 pi) is an analytic check of the
  inertial-period formula qT = 2 p pi. It is not a numeric control.

## 6. Citation mining (policy step 4)

The references [1]-[12] (p.237) were checked against `CORPUS_INDEX.md`, the digests and the filenames.

Held:
- [2] Arenstorf, AIAA J 1963 (filed by `#960`).
- [6] Breakwell & Perko 1974.
- [10] Perko 1967.

Not held, in priority order:
1. [1] Arenstorf, R. F. (1963), Amer. J. Math. 85:27-35, doi 10.2307/2373181 (CONFIRMED). The second-kind
   existence proof that both Perko and Hitzl cite.
2. [5] Breakwell, J. V. & Perko, L. M. (1965), "Matched asymptotic expansions, patched conics and the
   computation of interplanetary trajectories", Proc. XVI IAC; also AIAA Paper 65-689, doi
   10.2514/6.1965-689 (CONFIRMED for the AIAA version). First-order matching, Eqs. (A.2)-(A.13).
3. [9] Perko, L. M. (1964), "Asymptotic matching in the restricted three-body problem", PhD, Stanford. No DOI
   (thesis). The source of the error estimates.
4. [3] Arenstorf, R. F. (1963), "Periodic trajectories passing near both masses of the restricted
   three-body problem", Proc. XIV IAC, Springer. No DOI found (UNCONFIRMED). Figs. 1-9 are referenced in
   sec. 2.
5. [4] Birkhoff, G. D. (1915), "The restricted problem of three bodies", Rend. Circ. Mat. Palermo
   39:265-334, doi 10.1007/BF03015982 (CONFIRMED). The symmetric periodicity condition.
6. Textbooks, low priority: [7] Coddington & Levinson 1955; [8] Graves 1956; [11] Poincare, Methodes
   Nouvelles vol. 3; [12] Whittaker, Analytical Dynamics.
