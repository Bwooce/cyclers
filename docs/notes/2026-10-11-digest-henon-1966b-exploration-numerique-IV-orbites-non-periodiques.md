# Digest: Hénon 1966b, "Exploration numérique du problème restreint. IV. Masses égales, orbites non périodiques (fin)" (#960)

M. Hénon (Institut d'Astrophysique, Paris), Bulletin Astronomique (3e série) 1, fasc. 2, pp. 49-66 (1966),
DOI 10.3406/bastr.1966.14447. Manuscript received 15 February 1966. In French, with English, German and Russian
abstracts. PDF page n = journal page n + 46 (two Persée cover pages).
- Filed as `cyclers_pdf/papers/henon-1966b-exploration-numerique-probleme-restreint-IV-masses-egales-orbites-non-periodiques-bull-astron-1-49-doi-10.3406-bastr.1966.14447.pdf`.
  This is the Persée scan as uploaded by the owner (md5 5f085fc52203785d6232e786d7969946, 20 pp.).
- OCR copy: `...-french.pdf` (`ocrmypdf --force-ocr -l fra+eng`).
- English companion: `...-en-translation.pdf` (+ `.tex`, 20 pp.): an UNOFFICIAL full translation (private research
  copy) from the page images. All 15 equations are typeset and all 12 figure captions translated. The figures are
  page renders. Misprints are marked [sic].
- Checks filed beside the PDF: `...-check_henon1966b.py/.out` (our integration) and `...-witness.py/.out`.
- Printed numbers: `data/sources/henon-1966b-printed-numbers.yaml`.
- **Label clash.** In this paper's own bibliography, "Hénon 1966b" is C.R. Acad. Sci. 262:312 (force-free magnetic
  fields). In this project, "Hénon 1966b" means paper IV.
- Wanted list: the paper IV item of row 51; with this upload the row is complete.

## 0. Verdict

**Paper IV completes the mu = 1/2 map. It covers the outer orbits round BOTH bodies (classes l and m) and gives the
overall classification of all solutions (Fig. 12).**
- **Class l** (direct orbits of large radius round both bodies, C > 0) has a region of invariant curves. It is small
  compared with the accessible region, and the rest escapes.
  - The curves end at l2 (C = 3.3026), where l becomes unstable and meets Kristiansson's asymmetric class K.
  - They survive past l1 (3.428, a = -1), where l is unstable in the strict sense but still surrounded by closed
    curves.
- **Class m** (retrograde orbits round both bodies, C < 0) has a large region of curves. It shrinks regularly and
  disappears at m1 (C = -0.6164), the extremum of C on m.
- **The right-hand outer region (x > 0, far out) is all escape.** Eqs. (1)-(10) explain it with a two-body model in
  rotating axes: escape iff A - C/2 >= 0.
- **Overall picture (Fig. 12):**
  - Each of the six principal classes (f, g, h, i, l, m) is surrounded by a region of curves that ends where its
    orbits become unstable.
  - A semi-ergodic region exists round M1 and M2 for C > 3.4568.
  - Everything else is escape.
  - Regions round the other (mostly unstable) classes are very small.
- **Theory (pp.64-65).** The KAM bounds of Arnold (M < 10^-333) and Moser (M < 10^-48) are far below mu = 1/2, yet the
  curves exist. Theory does not explain this, nor where the region of curves ends.
  - Hénon's empirical rule: the region often ends at an unstable periodic orbit.
- **Erratum (p.66).** It confirms the paper III Fig. 12/13 caption interchange that the paper III digest marked.
- **Catalogue impact: none.** mu = 1/2. Section 5 gives what bears on our orbits round two bodies.

## 1. Content, section by section

- **Introduction (p.50).** Paper III covered the central region (f, g, h, i). This paper covers the rest of the
  section plane.
- **C positive, neighbourhood of l (pp.50-58, Figs. 1-7).**
  - At C = 4.5 the left outer region holds the stable l orbit (x near -4.5) and closed curves round it. Elsewhere the
    third body escapes.
  - The right outer region (x > 0) is all escape.
  - Eqs. (1)-(10): the two primaries are replaced by a unit mass at the origin. Then C = 2A - 2E (5), and escape is
    A - C/2 >= 0 (6). On the section this is (7). For x > 0 it always holds (9). For x < 0 it reduces to the dashed
    escape limit (10), xdot^2 >= 2/|x| - C^2/(4x^2).
  - At C = 4 (Fig. 2) the curves shrink, and two semi-ergodic orbits fill areas inside the escape limit. Hénon
    expects them to escape eventually.
  - Figs. 4-7 follow the region down to C = 3.30. Fig. 3 plots its x-limits against C.
  - **Origin symmetry (mu = 1/2 only).** l is symmetric about the origin. So every n-periodic class that meets l comes
    with its mirror, and the pair acts like one 2n-periodic class (Moser's beta vanishes). The a = 0 and a = -1/2 cases
    are therefore not special. Exception: at C = 3.348 (a = -1/2) the 3-periodic class eta_3 is its own mirror, so l is
    unstable there (the paper III case).
  - C = 3.7: 8 islands (two mirror chains of 4, asymmetric 4-periodic nuclei). C = 3.5: 5 exterior islands (the
    narrowing in Fig. 3).
  - **l1, C = 3.428.** a reaches -1, and l meets the double-periodic pair epsilon_2 and epsilon'_2. This is Moser's
    4-periodic exceptional case. The numerics (Fig. 5) show instability, yet closed curves of finite size surround l.
    That is "stability in a wider sense", foreseen by Moser 1958, Fig. 2b. The nearby asymmetric double-periodic
    classes zeta_2 and zeta'_2 hold them.
  - **l2, C = 3.3026.** Type 4, so the curves stretch vertically (Fig. 7). l meets Kristiansson's asymmetric class K,
    whose stable orbits give two small regions of curves (Fig. 7, C = 3.30).
- **C negative, neighbourhood of m (pp.58-61, Figs. 8-11).**
  - m is stable from C = -infinity to -0.6164.
  - At C = -4 (Fig. 8) the whole plane is accessible, and the picture resembles Fig. 1, since (10) is even in C.
  - As C rises the region contracts regularly and the curves lengthen vertically (Figs. 10-11). Fig. 9 plots its
    x-limits.
  - a passes 0 and -1/2 twice, and by the origin symmetry nothing special happens.
  - At m1 (C extremum) the region vanishes, bounded by the unstable m branch.
- **Other classes and overall picture (pp.61-63, Fig. 12).**
  - A complete classification "appears absolutely impossible".
  - Only the six principal classes have large regions of curves.
  - Fig. 12 merges III Figs. 10 and 19 with IV Figs. 3 and 9 into one (C, x) map.
- **Applications: double stars (p.63).**
  - Matter ejected at a collinear Lagrange point is always in the semi-ergodic or the escape region. Since
    semi-ergodic orbits pass close to a star from time to time, the matter either hits a star or escapes. No stable
    orbit round both stars is possible for it.
  - The habitability conclusion of paper II is confirmed: near each stable principal class there is a wide
    quasi-periodic region.
  - Earth-Moon and Sun-Jupiter mass ratios would be worth doing but much longer, because small mu is near-degenerate.
- **Theory (pp.64-65).**
  - Stable periodic orbits are the "nuclei" of the ordered motion. Unstable ones often fix the edge of a region of
    curves.
  - The same behaviour appears in Hénon-Heiles, in accelerators, in satellite oscillations and in force-free fields.
  - The Arnold and Moser theorems (eq. 11) need tiny perturbations: eqs. (12)-(14) give M < 10^-333 from Arnold 1963a,
    and (15) gives M < 10^-48 from Moser 1962. The curves exist at perturbations of order one, and the theory does
    not explain this.
  - Nothing predicts where the curves end.
  - The restricted problem's complexity, against the simple two-body problem, makes a closed-form general solution
    hopeless.

## 2. Printed numbers and checks

There are no tables. The printed values (24 entries) are in `data/sources/henon-1966b-printed-numbers.yaml`.
- **Witnesses:**
  - A: tesseract, fra, 300 dpi.
  - B: macOS Vision, 300 dpi.
  - O: the ocrmypdf text layer.
  - C: our reading of 300 dpi crops, which decides.
- Every value was read by C at 300 dpi. The one exception is the Fig. 1 panel label "C = 4.5", read at 150 dpi; its
  caption was read at 300 dpi.
- The machine counts differ in places, from dropped, merged or double-counted lines (listed in the YAML header). We did
  not compare every machine reading digit by digit.
- **Independent check (ours; same model as the paper III check).** Symmetric orbits are found by single shooting, with
  a = (trace of the 4x4 monodromy - 2)/2.

| Printed (paper IV) | Claim | Ours | Hénon 1973, six digits (held) |
|---|---|---|---|
| 3.428... | l1, a reaches -1 | min a = -1.0000000 at C = 3.427719, x0 = -2.131850 | l1 3.427719, -2.131850 |
| 3.348... | a = -1/2 on l | C = 3.348278, x0 = -1.980702 | not tabulated |
| 3.3026... | l2, a = +1 (type 4) | C = 3.302593, x0 = -1.907962 | l2 3.302594, -1.907964 |
| -0.6164... | m1, extremum in C | fold between -0.6165 (two roots, a = 0.900 and 1.105) and -0.6163 (none) | m1 -0.616363, -0.617017 |
| 10^-333 (eq. 14) | Arnold bound with n = 2 | e^4 328^-4 = 4.717e-9; to the 40th power = 10^-333.05 | - |

- l at C = 4.5 is at x0 = -4.539 (Fig. 1 shows it near -4.6). m at C = -4 is at x0 = -4.421 (Fig. 8 shows it near
  -4.4).
- Eqs. (5), (7) and (10) were checked by substitution. Eq. (8) is the weaker condition with xdot^2 dropped.
- Our l1 and l2 agree with Hénon 1973 to six digits, which makes them positive controls. All of paper IV's quoted
  values agree with Hénon 1973 to the printed digits.

## 3. Figures (described, not digitised)

- **Fig. 1:** the C = 4.5 left outer region, with curves round l, the accessible limit (dash-dot) and the escape limit
  (dashed).
- **Fig. 2:** C = 4, with curves round l and two semi-ergodic orbits (dots and crosses).
- **Fig. 3:** region of curves against C for l (C 3.3-4.0, x -3.2 to -1.8), with the characteristics of epsilon_2 and
  eta_3 and the points l1, l2.
- **Figs. 4-7:** upper half sections for C from 3.9 down to 3.30. The 8 and 5 island chains are at 3.7 and 3.5. Fig. 5
  shows zeta_2/zeta'_2 and the epsilon_2 points. Figs. 6-7 show the eta_3 points and the K island.
- **Fig. 8:** C = -4, with curves round m and the escape limit.
- **Fig. 9:** region of curves against C for m, with m1.
- **Figs. 10-11:** C = -3 to -0.68.
- **Fig. 12:** the whole (C, x) classification. Regions of curves (vertical hatching), non-accessible (horizontal),
  semi-ergodic (dotted) and escape (white), with all six characteristics, the critical orbits f3, g1, h3, i1, l2, m1,
  and L1, L3.

## 4. Misprints in the original (marked in the translation)

- **p.64:** "Arnold (1936b)" should be 1963b.
- **p.64:** "Zlatousov" in the text against "Zlatoustov" in the bibliography. The bibliography entry has no year; the
  text says 1964.
- **p.65:** "ramarquable".
- **p.52:** "l'étoile" (the star) where the third body is meant. This is probably a slip from the double-star framing.
- **p.66 erratum (Hénon's own):** paper III's p.72 and p.73 captions are interchanged.

## 5. Relevance to this project

- **What IV adds to III.**
  - III covered motion near one primary (f, g, h, i).
  - IV covers motion round BOTH primaries: l (direct) and m (retrograde), and the far-field escape geometry.
  - The orbits that encircle the pair are the mu = 1/2 analogue of "circumbinary" or outer orbits, not of cyclers.
  - The only orbits in III and IV that visit both bodies in turn are the semi-ergodic ones (III Fig. 12, and the two
    in IV Fig. 2, which probably escape).
  - No periodic orbit that shuttles between M1 and M2 appears in either paper's figures. The paper names the
    double-periodic and 3-periodic classes epsilon_2, zeta_2 and eta_3 near l, but only as bounds and islands of the
    l region.
- **A rule worth reusing in our stability screens.** "Unstable in the strict sense" (|a| > 1 at l1, or l beyond l2
  near K) need not mean "leaves the neighbourhood". Closed curves of finite size can surround an unstable periodic
  orbit when nearby stable classes exist (Moser 1958, Fig. 2b; III Fig. 21).
  - So a monodromy-only verdict on a cycler candidate's long-term confinement is conservative.
  - The converse, a stable index, does not guarantee a large confined region either: the l region is small at C = 4.5.
- **The escape criterion A - C/2 >= 0 (eqs. 5-10).** It is a clean far-field test for whether an orbit far from both
  bodies is bound. It applies to our two-body-like outer arcs as a sanity check. It is a mu-independent two-body
  approximation, so it is valid only when the third body is far from both.
- **Mass-ratio scope.** Hénon says the work should be redone for Earth-Moon and Sun-Jupiter but is much harder at
  small mu. That is the gap the later Hill-case papers V-VI (held) and his 1997/2001 generating-family books (held)
  address. Nothing here feeds #974 directly.
- **Catalogue: no impact, no novelty claim.**

## 6. Citation mining (reference list, pp.65-66)

| Reference | In CORPUS_INDEX? | On the wanted list? | Proposal |
|---|---|---|---|
| Hénon I, II, III; Hénon-Heiles 1964 | held | - | none |
| Kristiansson 1933, Publ. Copenhague Obs. 89 (asymmetric class K) | no | no | low: source of class K (also cited by paper II); added to the bundled low-priority row |
| Moser 1955, 1958; Arnold 1963b; Siegel 1956 | no | proposed in the paper III digest | as there |
| Arnold 1961 Soviet Math. Dokl. 2:247; Arnold 1963a Russ. Math. Surveys 18:9; Kolmogorov 1954 Dokl. 98:527; Moser 1962 Göttingen Nachr. | no | no | low: KAM theory background, no numbers we use |
| Hayli 1965 Ann. Ap. 28:49; Laslett-Symon 1956; Symon-Sessler 1956; Zlatoustov et al. 1964; Hénon 1966 C.R. 262:312 | no | no | none (other-field analogues) |
| Kuiper 1941 ApJ 93:133; Kopal 1956 Ann. Ap. 19:298; Gould 1959 AJ 64:136 | no | no | none (binary-star mass transfer) |
| Birkhoff 1927 (Dynamical Systems); Poincaré 1892 | Birkhoff 1915 held; 1927 book not held | no | none (textbooks) |

Proposed wanted-list rows: no new row beyond one bundled low-priority row: Moser 1955 and 1958, Contopoulos 1965,
Szebehely 1963 (from paper III) and Kristiansson 1933 (from this paper).
