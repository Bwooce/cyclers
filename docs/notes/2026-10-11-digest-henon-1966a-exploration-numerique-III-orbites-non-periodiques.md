# Digest: Hénon 1966a, "Exploration numérique du problème restreint. III. Masses égales, orbites non périodiques" (#960)

M. Hénon (Institut d'Astrophysique de Paris), Bulletin Astronomique (3e série) 1, fasc. 1, pp. 57-79 (1966),
DOI 10.3406/bastr.1966.14442. Manuscript received 29 November 1965. In French, with English, German and Russian
abstracts. The Persée header gives pp. 57-79 (not 57-80). PDF page n = journal page n + 54 (two Persée cover pages).
- Filed as `cyclers_pdf/papers/henon-1966a-exploration-numerique-probleme-restreint-III-masses-egales-orbites-non-periodiques-bull-astron-1-57-doi-10.3406-bastr.1966.14442.pdf`.
  This is the Persée scan as uploaded by the owner (md5 d0922998ce5c30787f350ad413ccfe0b, 25 pp., text layer on the
  cover pages only).
- OCR copy: `...-french.pdf` (`ocrmypdf --force-ocr -l fra+eng`, the flags used for papers I and II).
- English companion: `...-en-translation.pdf` (+ `.tex`, 28 pp.): an UNOFFICIAL full translation (private research
  copy) from the page images. All 10 equations are typeset. All 21 figure captions are translated, and the figures are
  given as page renders. Misprints are marked [sic].
- Checks filed beside the PDF: `...-check_henon1966a.py/.out` (our independent integration) and
  `...-witness.py/.out` (the three-witness counts).
- Printed numbers: `data/sources/henon-1966a-printed-numbers.yaml`.
- Wanted list: this is paper III of row 51 (Hénon III, IV and Matukuma). Paper IV (1966b, Bull. Astron. 1(2):49,
  DOI 10.3406/bastr.1966.14447) is still wanted.

## 0. Verdict

**Paper III of the Hénon mu = 1/2 series: the first surface-of-section map of ALL orbits for equal masses.** Papers I
and II gave the periodic orbits and their stability. This paper shows what fills the rest of phase space.
- **Three orbit types.** Quasi-periodic orbits lie on invariant curves of the section map. Semi-ergodic orbits fill an
  area of the section but not all of the accessible region (Hénon's word for chaotic orbits). Escape orbits leave to
  infinity once C < 3.4568... (the L2/L3 value).
- **Stable periodic orbits organise everything.** Each region of invariant curves surrounds a stable periodic orbit of
  papers I and II (f, g, h, i). It ends where that orbit goes unstable or where a resonant class bounds it.
- **The retrograde class h keeps a bounded quasi-periodic region far below the escape threshold.** The region exists
  down to h1 (C = 1.2175), and again from h2 (1.0594) to h3 (0.9556). The escape threshold is 3.4568. Hénon's
  conclusion (p.74): "it is not possible to conclude that a satellite is unstable by examining the Jacobi constant
  C alone."
- **The direct class i loses its region much earlier,** at i1 (C = 3.7388). There is a small second region between
  i14 (3.5664) and i15 (3.4473).
- **The paper has no numeric tables.** The brief expected tables of orbit classes and sections. There are none: only
  the numbers in the text and the figure labels. They are transcribed with a three-witness check in the YAML.
- **Catalogue impact: none.** mu = 1/2 is far from any Solar System pair. This is method and classification
  background (see sec. 5).

## 1. Content, section by section

- **Introduction (p.58).** Szebehely (1963) noted that earlier numerical work computed particular orbits. Hénon wants the
  overall picture. The method is the section y = 0, ydot > 0 in the plane (x, xdot) at fixed C (eq. 1). One figure shows
  many orbits at one C.
- **Accessible region (pp.58-60, eqs. 2-3, Fig. 1).**
  - Hénon's Jacobi constant has no mu(1-mu) term, so C_Hénon = C_std - 0.25. M1 is at x = -1/2 and M2 at x = +1/2.
  - There are 4 separate regions for C > 4, three for 3.4568... < C < 4, and one for C < 3.4568...
  - Our integration gives C(L1) = 4 exactly and C(L2) = C(L3) = 3.4567962 in this convention.
- **Central region (pp.59-62, Figs. 2-3).**
  - At C = 4.5 the stable orbits f, g, h, i are isolated fixed points, each ringed by closed curves.
  - h and f are retrograde about M1 and M2. i and g are direct about M1 and M2.
  - Fig. 3 is a quasi-periodic rosette about M1 (a precessing retrograde ellipse) inside an annulus that contains the
    periodic orbit h.
  - The section has two symmetries:
    - xdot -> -xdot, from the time-reversal symmetry (paper II Fig. 1);
    - near left-right symmetry, which is exact under the origin symmetry for mu = 1/2.
- **Rotation angle (pp.61-63, eqs. 4-8, Fig. 4).**
  - The points on a curve advance by a rotation angle alpha. When alpha/2pi = p/q the orbit is periodic with q
    revolutions. When it is irrational the orbit is quasi-periodic.
  - Near a stable symmetric periodic orbit, the linear map (eq. 4) with ad - bc = 1 (eq. 5) and a = d (eq. 6) gives
    cos alpha0 = a (eq. 7). The curves are ellipses with axes along x and xdot, with axis ratio sqrt(-c/b) (eq. 8).
  - Critical types 1, 2, 3 and 5 (c = 0) make the ellipses flat horizontally. Types 4 and 6 (b = 0) make them flat
    vertically.
  - Near an unstable orbit the curves are hyperbolae.
  - Fig. 4: alpha/2pi about h at C = 4.5 is between 0.06 and 0.08, with its minimum alpha0 at h.
- **Interior islands (p.64, Fig. 5).**
  - Fig. 2 looks like a full family of curves, which suggests a second isolating integral. That integral is close to the
    angular momentum about the near body (Contopoulos 1965 series).
  - Near every rational alpha/2pi = p/q, a curve breaks into a chain of q islands. Fig. 5 shows 9 islands about i at
    C = 4.5. There are islands within islands without end (Arnold 1963).
- **Semi-ergodic region (pp.64-70, Figs. 6-12).**
  - Below C of about 4.5, one orbit's points fill an area (Fig. 6, C = 4). Hénon calls these orbits semi-ergodic.
  - There is a 6-island exterior chain about a periodic orbit of six revolutions about M1 (Fig. 8 is the ribbon orbit).
  - The transition from curves to chaos is abrupt (Fig. 7, an orbit of transitional type).
  - The region of curves about i shrinks (Fig. 9) and vanishes at i1, C = 3.7388. Fig. 10 maps it against C. On the
    left it is bounded by the unstable branch of class i.
  - At C = 3.5 (Fig. 11) the semi-ergodic region is large. One orbit (Fig. 12) wanders near M1 and near M2 alternately,
    since the two regions join for C < 4.
- **Escape region (pp.72-74, Figs. 13-15).**
  - At C = 3 (Fig. 13) the region about h is shrinking, with a 7-island chain.
  - Below 3.4568 every orbit outside the curves escapes (Figs. 14-15, C = 2.85, to x = -30 and beyond). So the
    semi-ergodic region becomes the escape region.
  - The curves are not affected. Orbits on them "totally ignore the possibility of escape".
- **The commensurabilities 3 and 4 (pp.74-79, Figs. 16-21, eqs. 9-10).** This section follows the family h as C falls.
  - **a = 0 (alpha0 = 2pi/4) at C = 1.8669.**
    - This is Moser's fourth-order case (eq. 10): the orbit is stable if |p| < 1.
    - The curves stay closed, so h stays stable.
    - The 4-periodic class alpha_4 (Bartlett's "class alpha") bifurcates here. Its stable part gives the 4-island
      chain at C = 2 (exterior) and C = 1.9 (interior). The chain is then absorbed into the curves about h.
  - **a = -1/2 (alpha0 = 2pi/3) at C = 1.5183.**
    - This is Moser's third-order case (eq. 9, three-branched curves, Fig. 20a): h is unstable at this one C.
    - The region of curves shrinks to a point.
    - The 3-periodic class i_3 crosses here. Its unstable branch bounds the region, and its stable branch gives a
      3-island chain at C = 1.55 and 1.5183.
  - **h1 at C = 1.2175 (type 6).** h goes unstable here and meets the double-periodic class n'_2. The ellipses stretch
    vertically. The unstable n'_2 orbits bound the region.
  - **h2 to h3 (C = 1.0594 to 0.9556).** h is stable again.
    - The double-periodic class b_2 is stable from b_2 7 (1.135) to b_2 8 (1.0594). It gives 2 islands (Fig. 21).
    - At C = 1.07, curves surround h even though h is unstable, because the stable b_2 orbits are close by.
    - Past h3 (an extremum of C, type 1) the region of curves is gone for good.
  - The next article will treat the neighbourhoods of classes l and m. That is paper IV, 1966b, not held.

## 2. Printed numbers and checks

There are no tables. All 28 printed values are in `data/sources/henon-1966a-printed-numbers.yaml`, with page and
witness counts.
- **Witnesses:**
  - A: tesseract, fra, 300 dpi.
  - B: macOS Vision, 300 dpi.
  - O: the ocrmypdf text layer.
  - C: our reading of the 300 dpi crops, which decides.
- **Result:** no witness misread a digit. Where a witness dropped a line, C read the value on the crop. The dropped
  lines were the Russian abstract, the p.74 figure labels, and three p.78 values in B.
- **Independent check (ours; script beside the PDF).** This is a single-shooting symmetric family-h orbit at fixed C,
  with the 4x4 monodromy and a = (trace - 2)/2, rtol 1e-12.

| Printed (paper III) | What it claims | Ours | Hénon 1973, six digits (held) |
|---|---|---|---|
| 3.4568... | C at L2/L3 | 3.4567962 | - |
| 1.8669... | a = 0 on h | a = 0 at C = 1.866733 | not tabulated |
| 1.5183... | a = -1/2 on h | a = -0.5 at C = 1.518282, x0 = -0.95497 | not tabulated (paper II: x0 = -0.955, a = -0.5) |
| 1.2175... | h1, a = -1 | a = -1 at C = 1.216891, x0 = -1.0887 | h1 1.216891, x0 -1.088744 |
| 1.0594... | h2, a = -1 | a = -0.99998, x0 = -1.21037 | h2 1.059408, x0 -1.210359 |
| 0.9556... | h3, a = +1 (C extremum) | two roots, a = 1.084 and 0.918 (fold just below) | h3 0.955552 |

- **Our h1 equals Hénon 1973 to six digits.** It is a positive control for the method. The 1965-66 value 1.2175 is
  6e-4 high, as the Hénon 1973 digest already notes.
- **Our a = 0 point is at C = 1.86673,** 1.7e-4 below the printed 1.8669. No higher-precision published value is held,
  so we report this and do not call it an erratum.
- **Hénon 1973 against the paper-II values quoted here:**
  - i1, i14, h2, h3 and b_2 8 agree to the printed digits.
  - i15 (3.4473 against 3.445853) and b_2 7 (1.135 against 1.133574) differ by 1.4e-3.
  - h1 differs by 6e-4.

## 3. Figures (described, not digitised)

- **Fig. 1:** accessible region in (x, xdot) for the three C ranges.
- **Fig. 2:** the C = 4.5 section with f, g, h, i.
- **Fig. 3:** a quasi-periodic rosette about M1. The plot is printed rotated.
- **Fig. 4:** alpha/2pi against x about h.
- **Fig. 5:** the 9 interior islands.
- **Fig. 6:** the C = 4 section with the semi-ergodic sea and the 6 exterior islands.
- **Fig. 7:** an orbit of transitional type, at x 0.5 to 0.9, about g.
- **Fig. 8:** the 6-island ribbon orbit. The plot is printed rotated.
- **Fig. 9:** curves about i at C = 3.9, 3.85, 3.8 and 3.75.
- **Fig. 10:** region of curves against C for class i.
- **Fig. 11:** the C = 3.5 section with one numbered semi-ergodic orbit.
- **Fig. 12 (printed under the Fig. 13 caption, p.73):** that orbit in (x, y), visiting M1 and M2.
- **Fig. 13 (printed under the Fig. 12 caption, p.72):** the C = 3 section with the 7-island chain.
- **Figs. 14-15:** the escape orbit at C = 2.85.
- **Figs. 16-18, 21:** curves about h at C from 2.5 down to 0.97.
- **Fig. 19:** region of curves against C for class h, with the characteristics h, b_2, i_3 and alpha_4 and the
  critical orbits h1, h2, h3.
- **Fig. 20:** curve shapes at the 2pi/3 and 2pi/4 resonances.
- Any number read off Figs. 10 or 19 is approximate. None is used as data.

## 4. Misprints in the original (marked in the translation)

- **The Fig. 12 and Fig. 13 captions are interchanged (pp.72-73).** The text calls the C = 3 section Fig. 13 and the
  M1/M2 orbit Fig. 12.
- **p.75:** "au voisinage du point n" should be "point h".
- **p.57 résumé:** "les orbites demi-ergodiques" should be "semi-ergodiques". The English and German abstracts say
  semi-ergodic.
- **Fig. 16:** the x labels (1.1 ... 0.6) lack minus signs. The region is x = -1.1 to -0.6.
- **Persée cover pages:** the transcribed abstracts drop words ("s'ouvre les orbites", "occupée des orbites"). These
  are Persée slips, not Hénon's. The translation follows journal p.57.

## 5. Relevance to this project

- **#974 (Hénon 2003 asymmetric Hill families, scaled to moons): no direct input.**
  - #974 needs Hill-limit (mu -> 0) families. This paper is mu = 1/2.
  - Its Hill-case counterpart is paper VI (Hénon 1970, held: docs/notes/2026-10-07-digest-henon-1970-hill-case-non-periodic-orbits.md).
    Paper VI uses the same three-way classification (quasi-periodic, ergodic or semi-ergodic, escape) and the same
    section method.
  - What paper III adds is the method lineage and the rule: a region of invariant curves about a stable periodic orbit
    survives far past the zero-velocity opening, and only the orbit's own stability index bounds it. #974's stability
    claims about Hg and the asymmetric families can cite this as the origin of that rule.
- **The mu = 1/2 background cited by Hénon 1969 (H5).** H5 cites 1965a/b and 1966a/b for the equal-mass picture. With
  III now held, only IV (1966b, the neighbourhoods of l and m) is missing from that background.
- **Resonant and repeated-encounter orbits:**
  - Every island chain is a resonant periodic orbit about ONE primary. The chains are: 9 islands about i at C = 4.5;
    6 islands about M1 at C = 4; 7 islands at C = 3; alpha_4 (4-periodic) at C = 1.8669; i_3 (3-periodic) at
    C = 1.5183; b_2 (2-periodic) at C of about 1.07-1.135.
  - The paper ties each chain to alpha0/2pi = m/n on the parent family. This is the island-chain view of the n-periodic
    bifurcations of paper II.
  - The only orbit that visits both primaries is the semi-ergodic orbit of Figs. 11-12 (C = 3.5). It is chaotic
    transit, not a periodic or cycler-like orbit, and the paper finds no periodic M1-M2 cycler structure.
  - So there is no catalogue candidate and no novelty claim here.
- **A reusable caution for our own screens.** "Escape is energetically allowed" (C below the L2 value) does not mean
  "the orbit escapes". Hénon shows long-lived bounded motion down to C = 0.96 about a primary. This matches the
  retrograde-satellite results of H5/H6 and Hénon 1973.

## 6. Citation mining (reference list, p.79)

| Reference | In CORPUS_INDEX? | On the wanted list? | Proposal |
|---|---|---|---|
| Arnold 1963, Uspekhi Mat. Nauk 18:91 (Russian Math. Surveys 18:85) | no (only Arnold-diffusion titles) | no | low: KAM background; no numbers we use |
| Bartlett 1964, Publ. Copenhague Obs. 179 | no | yes, row 19 | none (already wanted); it is also the source of Hénon's "class alpha" (alpha_4) |
| Contopoulos 1965, ApJ 142:802 | no | no | low: third-integral series; background only |
| Hénon 1965a (I), 1965b (II); Hénon and Heiles 1964 | held | - | none |
| Moser 1955, Nachr. Akad. Wiss. Göttingen IIa 6:87 | no (Moser 1953 held) | no | low: the 2pi/3 and 2pi/4 exceptional cases (eqs. 9-10) |
| Moser 1958, AJ 63:439 | no | no | low-medium: English statement of the same resonance cases; free on ADS |
| Siegel 1956, Vorlesungen über Himmelsmechanik | no | no | none (textbook; later Siegel-Moser editions are standard) |
| Szebehely 1963, Space Sci. Rev. 2:219 | no | no | low: survey of numerical work up to 1963 |

Proposed wanted-list rows (not added by this digest; DOIs not looked up): Moser 1958 AJ 63:439, and Moser 1955
Göttingen Nachr. together in one low-priority row; Contopoulos 1965 ApJ 142:802 and Szebehely 1963 SSR 2:219 as
background, lowest priority.
