# Digest: Hitzl 1977, "Generating Orbits for Stable Close Encounter Periodic Solutions of the Restricted Problem" (#960)

D. L. Hitzl (Lockheed Palo Alto Research Laboratory), "Generating Orbits for Stable Close Encounter Periodic
Solutions of the Restricted Problem", AIAA Journal 15(10):1410-1418 (October 1977), DOI 10.2514/3.60808.
Crossref-confirmed 2026-10-05: 15(10), pp. 1410-1418.
- It was presented as AIAA Paper 76-840 at the AIAA/AAS Astrodynamics Conference, San Diego, 18-20 Aug 1976.
  Submitted 7 Sept 1976, revised 17 June 1977 (p.1410).
- Filed as `cyclers_pdf/papers/hitzl-1977-generating-orbits-stable-close-encounter-periodic-solutions-restricted-problem-aiaa-j-15-1410-doi-10.2514-3.60808.pdf`.
  9 pages, text layer (scan with OCR), md5 e828b8052c2f81a8b4ce85299a60f624.
- I read all pages from the text layer. Table 1 (p.1417) was read as a page image.

Companions (held, digested 2026-10-04):
- `docs/notes/2026-10-04-digest-hitzl-henon-1977-critical-generating-orbits.md`. This is "paper I", Celest.
  Mech. 15:421, the full critical-orbit table.
- `docs/notes/2026-10-04-digest-hitzl-henon-1977b-stability-second-species-orbits.md`. This is "paper II",
  Acta Astronautica 4:1019, the mu = 0 stability theory and the p.1039 plans.

Evidence tags: READ, COMPUTED, INFERRED as in the companion digests.

## 0. Verdict

This is an expository AIAA-audience summary of the Hitzl-Henon programme. It contains no new mu > 0
stability computation.

What it contains:
- The new result stated in the abstract: "for mu = 0, second species orbits are critical if, and only if,
  the Jacobi Constant C has an extremum" (p.1410). The text adds that the multiplicative factors introduce no
  extraneous solutions (p.1416-1417).
- The Henon-Guyot 1970 family-m example, which is the only mu > 0 stability data shown. It is reused from
  Ref. 11, pp.354 and 367-368.
- A six-row excerpt of the paper I critical-orbit table.

The `#938` R4 claim that the stability side of the programme was "announced, never delivered" therefore
STANDS for this paper. Its closing position is the same as paper II's: the critical generating orbits are
"the correct generating solutions from which to start an analytical (and numerical) search for stable close
encounter periodic solutions of the restricted problem for mu > 0" (p.1410). No such search is reported.

## 1. Content (READ)

- Context (pp.1410-1411):
  - Close-encounter periodic orbits "in general ... are highly unstable" (refs. 3, 4).
  - An "early theoretical stability analysis [Abraham 1967] was a failure".
  - Matched asymptotic expansions to O(mu^2) (Breakwell-Perko 1974) were extended "to obtain a stability
    analysis valid through O(mu)". But "very few stable close encounter orbits were known".
- Figs. 1-4 show two Earth-Moon (mu = 1/82.30) periodic orbits of Arenstorf, Davidson and Causey, in the
  rotating and inertial frames:
  - m/l = 1/2, n = 1. In the second-species limit this is the direct orbit of family C12 at tau/pi = 1,
    eta/pi = 2 (Fig. 1 caption).
  - m/l = 1/2, n = 2, near the critical generating orbit C24(1) at tau/pi = 1.99362, eta/pi = 3.98751
    (Fig. 4 caption).
  - These are figures only, with no states.
- Commensurability for Arenstorf's second-kind orbits: |l| T0 = m TM = 2 pi m, so a = |m/l|^(2/3), Eqs.
  (2)-(4) (p.1411).
- Henon-Guyot "single exceptional orbit" (pp.1411-1412):
  - The family m of retrograde simple-periodic orbits about both primaries has critical orbits m1 and m2.
  - It is stable between m1 and m2 and outside m3 for 0 < mu < mu0 = 0.327 (Fig. 6 caption, from Ref. 11,
    p.354).
  - Fig. 7 shows "the very narrow domain of stability". As mu -> 0, m1 and m2 coalesce to the critical
    generating orbit A0(-1), "a local maximum for the Jacobi constant C" (p.1412).
  - Fig. 8 traces five orbits of m1 at mu*.
- Stability index (pp.1412-1414):
  - k = AD + BC from displaced half-orbits. Numerical cautions: double precision near 1e-14, central
    differences, AD - BC = 1.
  - Or k = (Tr - 2)/2 from the monodromy matrix. Stability for |k| < 1.
  - The 3:1 (k = -1/2) and 4:1 (k = 0) resonances need finer analysis (p.1414).
- Critical generating orbits (pp.1414-1415):
  - Switch functions sigma0, sigma1, sigma2 (Eq. 15) and the timing condition F0 = 0 (Eq. 22, from Henon
    1968).
  - Solutions exist for all tau > 0.51500... and eta > 0 (p.1414).
  - C = 2 sigma1 sqrt(a(1 - e^2)) + 1/a (Eq. 26).
  - Extremal C along F0 = 0 by a Lagrange multiplier gives G* = 0 (Eqs. 28-30).
  - "the infinite set of extremal orbits possessing jumps in the stability index k is now accessible
    analytically" (p.1415).
  - Fig. 9 reproduces Henon's characteristics with "the approximate location of 62 critical orbits".
- Direct stability analysis (p.1416): a recurrence across one moon passage (Eqs. 31-38) gives the necessary
  condition Eq. (40). The reduction (Eq. 41) is "both necessary and sufficient for these second species
  limiting orbits to be critical". "Full details are available elsewhere [16]", which is paper II.
- "the direct orbits C_ij(1) all occur at maxima of C whereas all of the retrograde orbits C_ij(2) occur at
  minima of C" (p.1417). "The deflection angle delta is zero, in fact, only for orbits with integer values
  for tau/pi and eta/pi" (p.1417).

## 2. Table 1 (READ p.1417, page image), "Critical second species orbits for mu = 0"

| Name | Fig. | n* | tau/pi | eta/pi | s0 s1 s2 | a | e | x0 | x1 | C | T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C12(1) | 10 | 3 | 0.98725 | 1.97522 | + + - | 0.63289 | 0.58182 | -1.00112 | 0.26466 | 2.87412 | 6.20310 |
| C23(1) | 11 | 2 | 1.99451 | 2.99244 | - + + | 0.76343 | 0.30997 | -0.52679 | 1.00007 | 2.97130 | 12.53191 |
| C24(2) | 12 | 6 | 2.12158 | 3.77226 | - - - | 0.69628 | 0.57791 | 1.09867 | -0.29389 | 0.07423 | 13.33030 |
| C25(1) | 13 | 4 | 1.99366 | 4.97863 | - + + | 0.54459 | 0.83813 | -0.08815 | 1.00103 | 2.64132 | 12.52652 |
| C38(1) | 14 | 7 | 2.99629 | 7.98179 | + + - | 0.52112 | 0.92046 | -1.00078 | 0.04145 | 2.48322 | 18.82627 |
| C38(2) | 15 | 11 | 3.01234 | 7.94148 | + - - | 0.52636 | 0.91527 | -1.00812 | 0.04460 | 1.31530 | 18.92711 |

Cross-check against paper I's Table I as transcribed in the 2026-10-04 digest (COMPUTED comparison):
- C12(1) agrees in every printed field.
- Last-digit differences:

| Row | Field | This paper | Paper I (digest) |
|---|---|---|---|
| C25(1) | C | 2.64132 | 2.64131 |
| C38(1) | a | 0.52112 | 0.52111 |
| C38(1) | e | 0.92046 | 0.92048 |
| C38(1) | x1 | 0.04145 | 0.04144 |
| C38(1) | C | 2.48322 | 2.48318 |
| C38(2) | e | 0.91527 | 0.91526 |
| C38(2) | C | 1.31530 | 1.31528 |

- These differences are below 5e-5, which is a rounding or recomputation level. Paper I is the authority;
  the digest transcription should be re-checked against paper I's page image before either is used as a
  5-digit golden.
- C23(1) and C24(2) were not in the digest excerpt and were not compared.

## 3. Gate answer: `#948` R4, "announced, never delivered"

- Paper II (p.1039) promised two papers. One would give numerical k = +-1 boundaries at mu > 0, "especially
  ... mu* = 1/82.30 ... Earth-Moon ... examined in detail". The other would give second-order analytic
  boundaries.
- This 1977 AIAA paper is not either one:
  - It predates or parallels paper II: it cites paper II as "to appear" (ref. 16) and paper I as "to appear"
    (ref. 14).
  - Its only mu > 0 stability content is the reuse of Henon-Guyot 1970.
  - Its title promises "Generating Orbits for Stable ... Solutions", and that is what it delivers: mu = 0
    generating (critical) orbits, not stable orbits at a positive mass.
- **The R4 claim stands for the held corpus plus this paper.**
- Follow-up: forward citations of paper II after 1977, for example Hitzl's later Celest. Mech. or AIAA
  papers. This was not searched here.

## 4. Relation to `#906`/`#937` (the demanded-turn gate)

- The paper's direct stability analysis (p.1416) uses the moon-centred hyperbola's v_inf, impact parameter
  b, e_H and the deflection delta. Its statement that delta = 0 only at integer (tau/pi, eta/pi) (p.1417) is
  the second-species limit of a demanded turn.
- This is consistent with `#937`'s finding that the turn must be measured from body-relative velocities
  (inertial, Moon-centred), not rotating-frame velocities. The paper's p.1416 definitions are all
  Moon-centred inertial hyperbola quantities (INFERRED; the paper does not discuss frames).
- It supplies no turn-angle numbers beyond that. Nothing here changes the `#937` ruling.

## 5. Positive controls

- Table 1 above at mu = 0. Use paper I's 5-decimal values as the authority. These are the same orbits.
- The Henon-Guyot family-m criticality boundary for 0 < mu < 0.327 is figure-only here (Fig. 6/7). The
  source table is Henon & Guyot 1970, not held.

## 6. Citation mining (policy step 4)

The references 1-16 (p.1418) were checked against `CORPUS_INDEX.md`, the digests and the filenames.

Held:
- 1: Perko 1974 (filed by `#960`).
- 3: Broucke 1968, TR 32-1168.
- 6: Breakwell & Perko 1974.
- 10: Szebehely 1967.
- 12: Henon 1968.
- 14, 16: Hitzl & Henon 1977a and 1977b.
- 2: the AIAA J note by Arenstorf 1963 is filed by `#960`; the Amer. J. Math. paper is not held (below).

Not held, in priority order:
1. Henon, M. & Guyot, M. (1970), "Stability of Periodic Orbits in the Restricted Problem", in Giacaglia
   (ed.), Periodic Orbits, Stability and Resonances, pp.349-374, doi 10.1007/978-94-010-3323-7_33. Already
   `#938` sec. 5 item 12, where it is CONFIRMED. It is the only mu > 0 stability source used here.
2. Arenstorf, R. F. (1963), "Periodic Solutions of the Restricted Three Body Problem Representing Analytic
   Continuations of Keplerian Elliptic Motions", Amer. J. Math. 85:27-35, doi 10.2307/2373181
   (CONFIRMED). The full existence proof behind the 1963 AIAA J note.
3. Henon, M. (1965), "Exploration numerique du probleme restreint II. Masses egales, stabilite des orbites
   periodiques", Annales d'Astrophysique 28:992-1007. No DOI found (UNCONFIRMED). The Crossref search
   returned only part III (Bull. Astron. 1966, doi 10.3406/bastr.1966.14442).
4. Abraham, R. (1967), Foundations of Mechanics, Benjamin, pp.219-230. Book; the "failed" stability
   analysis.
5. Low priority (numerical methods and historical background): Fehlberg 1963 (Liege symposium) and NASA
   TR R-248 (1966); Darwin 1911, Scientific Papers IV; Poincare, Methodes Nouvelles (NASA TT F-452).
