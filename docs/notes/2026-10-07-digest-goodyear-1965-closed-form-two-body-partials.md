# Digest: Goodyear 1965, "Completely General Closed-Form Solution for Coordinates and Partial Derivatives of the Two-Body Problem" (#960 batch 30)

W. H. Goodyear (IBM Federal Systems Division, Houston), The Astronomical Journal 70(3):189-192 (April 1965),
ADS 1965AJ.....70..189G. Received 20 August 1964, revised 18 January 1965. Work for NASA Manned Spacecraft Center
(contract NAS 9-996).
- Given file: supplied file `1965AJ70189G.pdf`, 4 pp., md5 366dcade5586ac17fb863ea736edeecd. The text layer is poor OCR
  (formulas garbled). Page 1 starts with the end of the previous paper's references (a three-body paper; "Danby 1962",
  "Liapounov 1889", "Lovett 1911", "Su-Shu Huang 1960"). Those lines are not Goodyear's. The title follows.
- Filed as `cyclers_pdf/papers/goodyear-1965-completely-general-closed-form-solution-coordinates-partial-derivatives-two-body-problem-aj-70-189-ads-1965AJ-70-189G.pdf`.
- How I read it: all four pages as 170 dpi images; the matrices of Eq. 16, its velocity companion on p. 191, and Eq. 19 on 300 dpi
  crops (170 dpi for the rest). One smudge after "where psi vanishes at the reference time t0" (p. 190, right column) hides about one line
  of text; the sentence is complete without it. Nothing I quote depends on it.
- Wanted-list row 42 (Danby 1964; Goodyear 1965; Goodyear and Guseman 1964). Companion of the held Danby 1965 digest
  (`2026-10-05-digest-danby-1965-matrizant-keplerian-motion.md`, line 101 lists this paper as a candidate acquisition).

## 0. Verdict

**This is the root of the "universal-variable closed-form Kepler STM" that the project uses through Shepperd 1985.
It prints all 36 partials of the two-body state with respect to the initial state, valid for every conic, with no
element singularities. It gives no numerical example.**
- It is a four-page letter-style paper with formulas only: the universal Kepler solution (Eqs. 3-8), and closed forms for
  the 6x6 state transition matrix (Eqs. 11-16) and for the partials with respect to mu (Eqs. 18-19).
- Shepperd 1985 (held, `shepperd-1985-universal-keplerian-state-transition-matrix-celest-mech-35.pdf`) states in his first
  section that he improves on Goodyear. The Shepperd digest (line 67) already carries "the four 3x3 sub-blocks (Goodyear form)".
- Catalogue: no row depends on it. **Proposal only:** none. Cite it in the doc-string of
  `src/cyclerfinder/core/kepler_stm.py` as the origin of the formulation, if the owner wants the attribution (the module
  names Shepperd 1985 and Battin and does not mention Goodyear: grep of `src` for "goodyear" gives no hit).

## 1. Universal Kepler solution (READ p. 189, image)

- Constants (Eq. 3): r0 = |r0|, sigma0 = r0 . v0, alpha = v0.v0 - 2 mu / r0. Note the paper's alpha is **twice the
  specific energy**, i.e. alpha = -mu/a; it is not 1/a.
- Universal functions of psi (Eq. 4): s_n = sum_j alpha^j psi^(2j+n) / (2j+n)!, n = 0..3. Eq. 12 adds s4, s5 the same way.
- Kepler equation (Eq. 5): t = t0 + r0 s1 + sigma0 s2 + mu s3. Radius (Eq. 6): r = r0 s0 + sigma0 s1 + mu s2.
- Lagrange coefficients (Eq. 7): f = 1 - mu s2 / r0, g = (t - t0) - mu s3, fdot = -mu s1 / (r r0), gdot = 1 - mu s2 / r.
  State (Eq. 8): [r v] = [r0 v0][f fdot; g gdot].
- psi is Stumpff's regularising variable (Eq. 10: dt/dpsi = r), or in the elliptic case psi = (E - E0)/(-alpha)^(1/2)
  (Eq. 9; the page shows a square root of -alpha, so alpha < 0 for ellipses).
- Credits: Stumpff 1947 and 1959; Herrick 1960; Herget 1948; Goodyear says his form is "both valid for all values of mu
  and algebraically simpler". mu may be negative (repulsive force) or zero.

## 2. State transition matrix (READ pp. 190-191, image)

With U = s2 (t - t0) + mu (psi s4 - 3 s5) (Eq. 13), a = -mu r / r^3 (Eq. 14), a0 = -mu r0 / r0^3 (Eq. 15), and
the 3x3 blocks (columns of r then v are the 3x2 matrix P = [r v]; Q = [r0; v0]^T rows are 2x3):
- dr/dr0 = f I + U v a0^T + P [ -(fdot s1 + (f-1)/r0)/r0 , -fdot s2 ; (f-1) s1 / r0 , (f-1) s2 ] Q
- dr/dv0 = g I - U v v0^T + P [ -fdot s2 , -(gdot-1) s2 ; (f-1) s2 , g s2 ] Q
- dv/dr0 = fdot I + U a a0^T + P [ -fdot (s0/(r r0) + 1/r^2 + 1/r0^2) , -(fdot s1 + (gdot-1)/r)/r ;
  (fdot s1 + (f-1)/r0)/r0 , fdot s2 ] Q
- dv/dv0 = gdot I - U a v0^T + P [ -(fdot s1 + (gdot-1)/r)/r , -(gdot-1) s1 / r ; fdot s2 , (gdot-1) s2 ] Q
- Eq. 19: dr/dmu = P [ -s2/r0 ; U/r0 - s3 ], dv/dmu = [r v a] [ -s1/(r r0) ; s2/r0 ; U/r0 - s3 ].
- Eq. 17 lists the inverse-matrix symmetry identities (partials with t and t0 swapped). Goodyear says Eq. 16 was built to
  satisfy them. Eq. 11 gives the variational equations they solve (second derivative = (mu/r^3)(3 r r^T / r^2 - I) times the partial).
- Dots on f and g are small on the scan. I fixed which entries carry a dot by the finite-difference test below, and
  the signs by the 300 dpi crops (minus signs sit before the long fractions in the velocity blocks).

## 3. Check (COMPUTED)

Script `check_goodyear_stm.py` (outputs `check_goodyear_stm.out`, `check_goodyear_signs.out`): universal propagation
with Eqs. 3-8 (Newton on Eq. 5), central finite differences of the 6x6 map, units mu = 1.
- Elliptic (t = 2.3), hyperbolic (t = 1.7) and a long elliptic arc (t = 9.0, more than one revolution).
- The structure "diagonal + U v a0^T + P C Q" matches finite differences to 1e-10 in the first two cases
  (1e-6 in the third; I did not separate finite-difference error from series truncation) for all four blocks, and Eq. 19 matches to 5e-10.
- Every entry of C matches the printed form once the leading minus signs in the velocity blocks, and the dots
  listed in section 2, are used. Without those minus signs four entries differ in sign. Eq. 19 with the opposite sign
  fails by 8 (elliptic).
- So the printed formulas are right as read, and the formulation reaches elliptic, hyperbolic and multi-revolution arcs.
  This is a self-consistency check against finite differences, not against a published number: Goodyear prints none.

## 4. Does the project have an equivalent?

- Yes. `src/cyclerfinder/core/kepler_stm.py` implements the same universal-variable STM (`shepperd_stm`, quadrant form
  of Ellison 2018 Eq. 17, secular term C, Stumpff c4 and c5), validated by central differences against
  `core/kepler.propagate`. It is the Shepperd 1985 form, which the Shepperd digest lists as an improvement on Goodyear's. `core/_stumpff.py` has c2 and c3; c4, c5 are in `kepler_stm.py`.
- Mapping of variables (my derivation from Eq. 6 against the project's r = r0 U0 + (r0.v0/sqrt(mu)) U1 + U2): psi = chi / sqrt(mu)
  and Goodyear's alpha psi^2 = -z. Goodyear's alpha is -mu/a; `kepler_stm.py` uses alpha = 1/a. So the s_n are the
  project's U_n up to powers of sqrt(mu) and a sign in z. They are not drop-in compatible.
- Goodyear's s_n (Eq. 4) are power series. The project evaluates Stumpff functions with closed forms and a small-z series guard.
- No test fixture needs this paper. Danby 1965 gives the same matrix in orbital elements and Cartesian state (digest). A third
  independent form is therefore available for cross-checks of `shepperd_stm`.

## 5. Citation-mining

- Not held: Stumpff 1947 (Astron. Nachr. 275:108), 1959 (Himmelsmechanik, ch. V), 1962 (NASA TN D-1415); Herrick 1960 (ASTIA
  AD 250 757); Sconzo 1963 (Mem. Soc. Astron. Ital. 34(2), the 36-term Jacobian); Herget 1948; Sundman 1912 (Acta Math. 36:127).
  None is on the wanted list except through row 42 (Goodyear and Guseman 1964; Danby 1964).
- HELD: Shepperd 1985 (above). Danby 1965 is HELD (`danby-1965-matrizant-keplerian-motion-aiaa-j-3-769...`).
- New candidate that matters: Stumpff 1962 NASA TN D-1415 (NTRS, free) for the original universal-variable ephemeris
  formulas. Low priority.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
