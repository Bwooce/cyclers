# Digest: Szebehely 1967, "Burrau's Problem of Three Bodies" (PNAS) (#960 wanted-list row 34)

V. Szebehely (Yale University Observatory), Proceedings of the National Academy of Sciences 58(1):60-65 (July 1967),
doi 10.1073/pnas.58.1.60, PMC335596. Communicated by G. M. Clemence, May 12, 1967. Astronomy section.
- **File: our own assembly, not a publisher PDF.** The PDF routes (PNAS, Europe PMC render, PMC "pdf/pnas00677-0079.pdf")
  are bot-gated (fetch-log). The PMC article page shows the six scanned pages as PNG images on the NCBI CDN
  (`cdn.ncbi.nlm.nih.gov/pmc/blobs/0428/335596/.../pnas00677-0079.png` to `-0084.png`, 816 x 1261 px, 4-bit grey, about
  100 dpi). We joined them with `img2pdf` (md5 3d62d142eda09f2d1a2ab5b2b6abbe75, 603,785 bytes, 6 pp.) and made an
  `ocrmypdf --force-ocr -O 3 -l eng --image-dpi 100` copy (md5 806a8592e42bd83adec541dba851aa2c, 754,991 bytes). The md5s
  are of our assembly. PNG md5s: 0079 68803bfac5101c03b46bbb48f947c02d, 0080 8093f8c8b11b3c7c1b47535016ff53d1, 0081
  5cb167fdcdcbfe50066106b8ca6b82e9, 0082 e17ae064484b7edddc2eca41e384407f, 0083 e70c4dd561a6df47f206b956f78992c5, 0084
  32ba3262e88ceb6d0a37c3ad63dfb5bb.
- Dark-pixel fractions at 100 dpi, assembly/OCR: p.1 0.0881/0.0881, p.3 0.0941/0.0941, p.6 0.0386/0.0386. Page 1 image
  shows the title and author.
- Proposed filename: `szebehely-1967-burraus-problem-three-bodies-pnas-58-60-doi-10.1073-pnas.58.1.60.pdf` (the OCR copy).
  If the owner prefers a publisher PDF, PNAS or Europe PMC will give one in a browser.
- PDF page n = journal page 59 + n. **How it was read:** all six pages on the page images; the text layer only to find things.
- Checks: `work/extra_checks.py/.out` (unit of time).

## 0. Verdict

**A short announcement, with no numerical table: the Pythagorean (Burrau) problem has an elliptic-hyperbolic asymptotic
solution, and the near-collision at t about 15.8 suggests a periodic orbit.** It is "Paper 1" of the Szebehely-Peters
pair (both held and digested). It adds history, a symmetry proof and a physical scaling; it adds no data.
- **It is not a witness for the 1967b x_2 slip.** The 1967b digest asked whether this paper repeats the periodic-orbit
  initial conditions (to confirm -0.0129612186 against the printed -0.0129612126). It does not: the paper has no
  coordinates, and it predates the periodic orbit. Proposal: close that open item as "PNAS checked; no table".
- **Citation correction (PROPOSAL):** both held AJ papers cite it as "Proc. Natl. Acad. Sci. U. S. 57, July" (read on the
  image of the 1967b reference list, p.4 of the scan; the 1967a list has the same entry in its text layer). The published
  paper is vol. 58, no. 1 (July 1967), pp. 60-65. "57, July" was an in-press reference: the July 1967 issue opened volume
  58. Same paper, not a different one: the 1967b note calls it "Paper 1" on Burrau's problem, and this paper is that.
  (The 1967a digest attributes "restricted-problem experience" with regularisation to it; this paper has no
  regularisation content, so that attribution belongs to Szebehely 1967, AJ 72:370, the other citation in the same
  sentence; INFERRED.)
- Nothing for the catalogue. `#931(a)` keeps the 1967b periodic orbit as its control.

## 1. Content

- **Sec. 1-2 (p.60).** Meissel (1893) proposed the problem; Burrau (1913) integrated it. Definition: masses m1 = 3, m2 = 4,
  m3 = 5 at rest at the vertices of a 3-4-5 right triangle, each mass opposite the side of the same length. Angular
  momentum L = 0. L = 0 is necessary, not sufficient, for a triple collision (Wintner pp.253, 265); so triple collision
  is possible in principle, and existence is settled by numerical experiments at several accuracies with binary-collision
  regularisation.
- **Sec. 3 (pp.60-61).** Sundman's and Poincaré's convergent series are "of questionable value" for describing motion;
  computer solutions suffer from end-figure errors. Definition adopted: "The solution is a qualitative and quantitative
  description of the general and detailed behavior of a dynamical system for the time -inf to +inf" (italic in the
  original), with Birkhoff's emphasis on qualitative behaviour.
- **Sec. 4 (pp.61-62).** History of periodic and asymptotic solutions of the general problem: Lagrange; Wintner,
  Whittaker, Lovett; E. Strömgren 1900-1922 (his first and second papers describe an asymptotic orbit related to this
  solution); Merman's term "elliptic-hyperbolic"; Zumkley (1941) and García (1966); Bohlin (1923) inconclusive. Recent
  quasi-periodic and periodic results: Arnold, Moser, Jefferys and Moser, Arenstorf. "The existence of periodic orbits in
  the general problem of three bodies, when no restrictions are placed on the masses and distances, seems to be still
  limited essentially to orbits of librational type" (italic).
- **Sec. 5 (p.62).** The solution: after several close approaches m2 and m3 form a close binary; the binary's centre of
  mass and m1 depart in opposite directions on hyperbolic relative motion. By symmetry x(t) = x(-t). Between -t0 and t0
  the largest mass m3 has alternating close approaches with the two smaller ones; none between the smaller masses. One
  m3-m2 approach is almost a collision with nearly rectilinear orbits while m1 is nearly at rest; this suggests a periodic
  orbit with zero initial velocities near (not at) the Pythagorean configuration, with a collision of m2 and m3 at T/2 and
  "T/2 < t0". Meissel's proposition "deserves admiration".
- **Sec. 5 end (p.63).** Details of the asymptotic solution, the integration method and the close approaches are
  deferred to a paper with C. F. Peters "which will be submitted" to the Astronomical Journal (the held 1967a paper).
- **Sec. 6 (p.63).** G. M. Clemence and "Dr. L. H. Stretensky of the U.S.S.R. Academy of Sciences" (as printed) asked
  about the symmetry. Theorem: "If the solution of the gravitational
  problem of n bodies with initially zero velocities exists, then this solution is symmetric with respect to the initial
  configuration." Sketch: Taylor series in t about t = 0; with xdot(0) = 0 every odd derivative vanishes (third: (df/dx)
  xdot = 0, and so on); continuation through binary collisions is assumed.
- **Sec. 7 (p.64).** GM T^2/L^3 = 1 in the numerical work. With M = M_sun/4 and L = 25 AU the triangle is 75, 100, 125 AU
  with masses 3/4, 1 and 5/4 M_sun, and the unit of time sqrt(L^3/GM) "≅ 40 years". The binary forms at t0 = 2400 years;
  the system tends to periodic motion with T = 1320 years, "although at T/2 the close approach mentioned in this paper will
  in all probability be catastrophic".
- Footnote 2 (p.64): integrations to t = 23 by M. Standish (Yale, 1966); to t = 23 in regularised variables by R. Spinelli
  with M. Lecar and Szebehely (NASA Institute for Space Studies, 1966); to t = 16 by L. Stanek under E. Stiefel (ETH
  Zürich, 1966).

## 2. Cross-checks (ours)

- Unit of time (`extra_checks.out`, Gaussian constant): sqrt(25^3/0.25) / k / 365.25 = 39.789 yr. Printed "≅ 40": agrees.
- t0 = 2400 yr = 60.32 units. The 1967a digest has the final binary formation and escape at t = 59.2 to 60: agrees.
- T = 1320 yr = 33.17 units. The 1967b periodic orbit has T = 31.8229622453 units = 1266 yr. The PNAS figure is an
  estimate made before the periodic orbit was computed, 4 percent long. "T/2 < t0" holds (15.9 < 60).
- The masses, triangle and zero velocities agree with the held Burrau 1913 digest and the 1967a digest. The quoted
  near-collision (m2-m3, m1 nearly at rest) is the t10 = 15.830 approach of 1967a (r23 about 4e-4).
- Burrau's other venue (footnote 1): "Publikationer og mindre Meddelelser fra Københavns Observatorium, No. 13 (1913)",
  the same paper as A.N. 195:113 (held).

## 3. Citation mining

Checked with `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md`:

| Reference (footnote) | Held? | Wanted list |
|---|---|---|
| 1 Burrau 1913, A.N. 195:113 (also Copenhagen Publ. No. 13) | HELD (`burrau-1913-...`) | n/a |
| 3 Wintner 1941, Analytical Foundations | NOT held | no row (textbook) |
| 4 Sundman 1912, Acta Math. 36:105 | NOT held (index "sundman" hits are about the Sundman transformation) | no row; low priority |
| 5 Poincaré 1884, Acta Math. 4:201 | NOT held | no row |
| 6 Birkhoff 1920, Science 51:51 | NOT held | no row |
| 7 Arnold 1963, Usp. Mat. Nauk 18:13 | NOT held | no row |
| 8 Moser 1966, SIAM Rev. 8:145 | NOT held | no row |
| 9 Jefferys and Moser 1966, AJ 71:568 | NOT held (Jefferys 1965 and 1974 are held) | no row. Quasi-periodic three-body orbits; low priority |
| 10 Arenstorf 1966, J. Reine Angew. Math. 221:113 | NOT held (Arenstorf 1963 x2 held) | no row; low priority |
| 11-12 Whittaker 1900, Encyklopädie 1905-1923 | NOT held | no row |
| 13 Lovett 1911, Quart. J. Math. 42:252 | NOT held | no row |
| 14-19 E. Strömgren 1900, 1909 (A.N. 182:181), 1919 (MN 80:12), 1921, 1922, 1923 | NOT held (Strömgren 1933 Bull. Astron. 9:87 held) | no row. 1909 and 1919 are the asymptotic-orbit precursors; low priority |
| 20 Merman 1958, Bull. Inst. Theor. Astron. 6:687 | NOT held | no row |
| 21-22 Zumkley 1941, A.N. 272:66; García 1966 | NOT held | no row |
| 23 Bohlin 1923 | NOT held | no row |
| 24 Szebehely 1966, AJ 71:968 (power series, restricted problem) | NOT held | no row; low priority |

- Wanted-list row 34 can be removed (paper in hand). The not-found note of fetch2-sonnet (captcha) is resolved through
  the PMC page images.

*Filed as `cyclers_pdf/papers/szebehely-1967-burraus-problem-three-bodies-pnas-58-60-doi-10.1073-pnas.58.1.60.pdf`.*
