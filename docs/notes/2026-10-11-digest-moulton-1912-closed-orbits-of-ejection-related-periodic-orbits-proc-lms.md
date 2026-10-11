# Digest: Moulton 1912, "Closed Orbits of Ejection and Related Periodic Orbits. Painleve's Theorem" (Proc. London Math. Soc. s2-11:367-397)

**Digested:** 2026-10-11, owner upload, read in full. The PDF is a 400 dpi JBIG2 scan with a hidden OCR layer from the
source. The layer garbles most equations, so every equation used here was read on page renders (110 to 260 dpi crops)
and the layer served as a digit check only.
**Citation:** Forest Ray Moulton (Research Associate of the Carnegie Institution of Washington, Chicago), Proceedings of
the London Mathematical Society, series 2, vol. 11, pp. 367-397. Received and read 13 June 1912 (the running heads read
"1912.]" and "[June 13,"). The year is 1912 and the pages are 367-397, 31 pages, as the wanted list's Crossref note says. The
file name carries 1913; the paper itself says 1912. DOI 10.1112/plms/s2-11.1.367 is not printed on
the scan; it is the one recorded as confirmed on the wanted list. The page heading prints the title with "Painleve's Theorem" added.
**Filed as** `cyclers_pdf/papers/moulton-1912-closed-orbits-of-ejection-related-periodic-orbits-proc-lms-s2-11-367-doi-10.1112-plms-s2-11.1.367.pdf` (31 pages, scan with OCR layer, md5 ef24c394f72c11783b719b00c4e97239).
The check script and its output are filed beside it as `cyclers_pdf/papers/moulton-1912-closed-orbits-of-ejection-related-periodic-orbits-proc-lms-s2-11-367-doi-10.1112-plms-s2-11.1.367-check_moulton.py` and `-check_moulton.out`.
Value transcription: `data/sources/moulton-1912-values.yaml`.

## Verdict

An analytic existence paper with no tables. It proves, for a small mass ratio, that the restricted problem has closed
orbits of ejection (ejection from the big body and collision with it, in the opposite direction, along the line of the
two bodies), and that each is the limit of two families of periodic orbits that pass near the big body. It also extends
the periodic orbits to ones that recede and return many times, with period 2 pi p/q. Its numbers are series
coefficients, closed-form conditions and two comparisons with older numerical orbits (Burrau and Darwin). No orbit in it is
at a finite mass ratio with printed initial conditions. It is not a golden source. It is the earliest source in the corpus for
the arcs that Bruno later calls the K-families.

## What the paper does

- **Two-body ejection (secs. 2 to 3).** Radial motion under the big body, series in tau = (t - t_j)^(1/3). The series
  (9) has alternating signs. It has a three-branch point at ejection, so the analytic continuation through collision
  turns complex and the closed orbit of ejection is not itself periodic.
- **Closed orbit condition (sec. 4, eq. 16).** In rotating axes the two-body radial orbit returns to the big body in the
  opposite direction if its ejection-to-collision time is a multiple 2 j pi of the frame period. That fixes the greatest
  distance: xi0^(3/2) = 2 j sqrt(2(1 - mu)). Fig. 1 draws j = 1, 2, 3.
- **Three bodies (secs. 5 to 8).** Levi-Civita's regularising variables remove the singularity at ejection. The solution is
  a power series in the initial constant a, the time parameter delta and mu. The two conditions that the orbit cross the
  x axis perpendicularly at half the period have a nonzero Jacobian (eqs. 41 to 44), so by the implicit function theorem a
  closed orbit of ejection exists for |mu| small, for every integer j, for ejection in either direction.
- **Painleve's theorem (sec. 9).** For small mu there are two orbits of ejection through any chosen point with any chosen
  speed above the fixed-axes speed of that point. Moulton gives this as his form of the theorem Painleve conjectured and
  Levi-Civita proved near a body.
- **Large mu (sec. 10).** Continuation argument, not a proof. The real orbit persists for mu from 0 to 1 unless it folds
  to pass near the other body, which Moulton calls extremely improbable but not certain. He therefore says closed orbits
  of ejection exist from both bodies. For mu near 1 the loops around the small body become circles of radius j^(2/3).
- **Related periodic orbits (sec. 11).** For mu = 0, orbits that cross the x axis perpendicularly near the big body with
  period 2 j pi are periodic and symmetric. As the nearest approach goes to zero they tend to the closed orbit of
  ejection. The direct-motion ones have loops (Fig. 3, periods 2 pi and 4 pi). The same Jacobian argument gives the
  periodic orbits for small mu, as power series in mu. Two families, direct and retrograde near the body, share one limit.
- **Many near approaches (sec. 12).** Orbits whose period in fixed axes is P = 2 pi p/q, so that q spacecraft
  revolutions take p frame rotations, and whose rotating-frame period is 2 pi p. The same existence argument applies.
  The simplest is p = 1, q = 2 (Fig. 4). If the near-approach distance a goes to zero, the limit is Darwin's orbit.
  Moulton states that nothing shows the second apsidal distance b vanishes at the same time, and doubts that two
  apsidal distances vanish together because of the lack of symmetry.

## Values transcribed and our checks

All in `data/sources/moulton-1912-values.yaml`. No tables exist.
- **Series (9).** Coefficients 3/7, 23/63, 1894/4851 and 3293/7007 (with signs and powers of a/c) reproduce exactly (ours:
  rational-arithmetic solution of eq. 7). Three witnesses.
- **0.18 pi (eq. 38).** Reproduces (0.1817 pi, ours). The OCR layer misreads it as 0.187 pi.
- **Darwin's orbit.** C-bar = 1.716 from eq. (30) with mu = 1/11, xi0 = 2^(1/3) reproduces (1.7158, ours). The paper's
  own value for Darwin's constant is 20/11 = 1.818 in Moulton's units.
- **Burrau's orbit does not reproduce.** Moulton prints C = 2.38 against Burrau's -2 C_B = 2.2528, for mu = 0.5 and
  xi0^(3/2) = 2 sqrt(2(1 - mu)) = 2, with the undisturbed greatest distance 2. With eq. (45) as printed I get 1.75, 1.88 or
  1.65 for the three readings of xi0, never 2.38. The two printed statements about xi0 also differ (xi0^(3/2) = 2 and xi0 = 2).
  The digits are as printed on a 260 dpi crop. This is a probable slip or a missing factor in the comparison, not an
  erratum claim, and it does not touch the existence proof.

## Relevance to this project

- **Bruno 1994 and the second-species generating families.** Bruno's arcs of Chapter IX that collide with P1, the
  rectilinear c = 0, e = 1 K-families (KS_k, K_k^j, ...; Appendix Tables 5 to 8), are Moulton's closed orbits of ejection at
  mu = 0 in modern form: one radial two-body arc from the big body, returning in the opposite direction after j frame
  rotations. Moulton gives the closed-form greatest distance (eq. 16) and the Jacobi constant (eq. 29) for them, and the
  existence proof at small mu. Bruno's second-species orbits are different: their generating arcs collide with the small
  body P2 (Henon's S families, Bruno's T families). Moulton's proof covers ejection from the big body only. Ejection from
  the small body is covered by his continuation argument in sec. 10, not by a proof, so it is the weaker half and it is the
  half that second species needs. So Moulton is the origin of the K-family arcs and a precursor of the second-species
  idea, not a source of second-species generating data. Krasinskii 1963, named in the same wanted row, is the later
  K-family source.
- **Collision-orbit background rows (Olle 1989, Llibre and Pinol, Martinez-Seara et al. 2021).** Their ejection-collision
  orbits are the multi-collision generalisation of this paper's single ejection-collision orbit. Moulton's own open
  question in sec. 12, whether both apsidal distances of the p = 1, q = 2 orbit can vanish together (a double ejection-
  collision orbit), is the n-ejection-collision problem those papers study. Moulton also gives the implicit-function route
  (nonzero Jacobian from the Levi-Civita regularised flow) that later existence results follow.
- **Resonant or repeated-encounter objects.** Yes in a limited sense, with the project's p:q label defined as spacecraft
  revolutions to lunar revolutions. (a) The closed ejection orbits and the sec. 11 periodic orbits have fixed-axes period j
  frame rotations, so 1:j, one near pass of the big body per period. (b) The sec. 12 orbits have q near passes of the
  big body per period 2 pi p, a q:p resonance, so Moulton's p = 1, q = 2 is the project's 2:1, with two near passes
  per period. Darwin's limit orbit has period pi, which is 2:1 as well. In Earth-Moon terms the big body is the Earth, so
  these are near-Earth passes, and for mu small no lunar pass is claimed. No orbit here has a lunar encounter, so none can
  collide with a project Earth-Moon family (A, B, D or the symmetric ones) as a literal row. They are the mu to 0 limit
  of the Casoliva class-1 orbits with a very small perigee, in the same way Bruno's K-families are.
- **No catalogue change.** No numeric orbit, no collision, and no candidate-novel row is touched.

## Figures (described)

- **Fig. 1 (p. 373).** Three closed orbits of ejection in rotating axes for small mu, j = 1, 2, 3. Solid for ejection toward
  +x, dotted for -x, each symmetric about the x axis and crossing it perpendicularly at the mid-point. Loops grow with j and the
  curves are drawn at different scales. The big body is marked 1 - mu and the small body mu on the x axis.
- **Fig. 2 (p. 391).** The mu = 1 - epsilon case, j = 1 and j = 2: a loop around the near-tiny big body and a large
  near-circular loop around the other body.
- **Fig. 3 (p. 392).** Periodic orbits for periods 2 pi (j = 1) and 4 pi (j = 2) passing near the big body, with a small loop at
  the near pass for the direct-motion class and a large heart-shaped outer curve.
- **Fig. 4 (p. 396).** The p = 1, q = 2 orbit for retrograde motion near the big body: two large lobes, above and below the
  x axis, joined at a small loop near the big body with the near-approach points a and b marked on the axis.

## Citation mining

Held: Szebehely 1967 (cites Moulton, per its digest); Burrau 1913 is held but is not the memoir cited here.
Not held and not wanted before today (one new low-priority row, Darwin 1909 is already in the remaining-items row):
- Levi-Civita, T. (1906), "Sur la resolution qualitative du probleme restreint des trois corps", Acta Math. 30:305-327
  (the paper cites 306-327 and 305-327 in two footnotes). The regularisation Moulton depends on, and the proof of
  Painleve's conjecture near a body.
- Painleve, P., Lecons sur la theorie analytique des equations differentielles (Stockholm lectures), p. 583. The conjecture.
- Burrau, C. (1894), "Recherches numeriques concernant des solutions periodiques d'un cas special du probleme des trois
  corps", Astron. Nachr. 135 (no. 3280) and 136 (no. 3251). The numerical orbit compared in sec. 8.
- Darwin, G. H. (1909), "On certain families of periodic orbits", Mon. Not. R. Astron. Soc. 70:134. The orbit with period pi.
