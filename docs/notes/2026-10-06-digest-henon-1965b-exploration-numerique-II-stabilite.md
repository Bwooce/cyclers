# Digest: Hénon 1965b, "Exploration numérique du problème restreint. II. Masses égales, stabilité des orbites périodiques" (#960 batch 21)

M. Hénon (Institut d'Astrophysique, Paris), Annales d'Astrophysique 28(6):992-1007 (1965), ADS
1965AnAp...28..992H. In French, with English and Russian abstracts.
- Filed as `cyclers_pdf/papers/henon-1965b-exploration-numerique-probleme-restreint-II-masses-egales-stabilite-orbites-periodiques-ann-astrophys-28-992-ads-1965AnAp-28-992H-french.pdf`.
  - This is an `ocrmypdf --force-ocr -l fra+eng` copy of the image-only ADS scan. The source scan is
    md5 81e76dbdfbe797bde1e612043c3155b9, 16 pp.
- The English companion is `...-en-digest.pdf` (+ `.tex`, 11 pp.). It is a detailed English digest, NOT
  a full translation. The translator subagent declined a sentence-by-sentence translation of the
  copyrighted article (see sec. 3).
  - It reproduces the English abstract verbatim, all 38 numbered equations and all 10 tables.
  - Its check log: 336 table cells checked against 300 dpi crops, 0 mismatches.
  - Figures are described with their numeric labels, not embedded.
- I read the class h and class i tables (p.1002) on the page image myself, and they agree cell for cell
  with the digest PDF. I read the class f/h text (pp.1001-1002) from the OCR.
- Was on the wanted list (rank 18 before batch 21); removed in this batch.

## 0. Verdict

**This is the published numeric substitute for the absent Bruno & Varin Table 8 (family h at
mu = 0.5).** KIAM preprint 64/2005 says its Table 8 "corresponds" to the left table on p.1002, with values
differing by about one unit in the fourth digit. Conversions: Hénon's x0 = x1(0) - 0.5 and
C_H = C - 0.25, and his stability index is a = Tr/2.

**Class h critical orbits, mu = 1/2** (p.1002, read on the page image; decimal commas converted):

| Orbit | C | x0 | Type | Bruno-Varin k |
|---|---|---|---|---|
| h1 = n'_2 5 | 1.2175 | -1.0875 | 6 | 1 |
| h2 = b_2 8 | 1.0594 | -1.2103 | 5 | 2 |
| h3 | 0.9556 | -1.4410 | 1 | 3 |
| h4 | 1.938 | -2.246 | 1 | 6 |
| h5 | 1.934 | -2.259 | 6 | 7 |
| h_2 6 | 0.183 | -2.563 | 6 | 9 |
| h_2 7 | 0.178 | -2.600 | 1 | 10 |
| h_2 8 | 1.536 | -3.387 | 1 | 13 |
| h_2 9 | 1.535 | -3.392 | 6 | 14 |
| h_3 10 | -0.054 | -3.640 | 6 | 16 |
| h_3 11 | -0.055 | -3.653 | 1 | 17 |

The k column is the Bruno-Varin correspondence from KIAM 64/2005 sec. 7.

- **The class h orbits are retrograde satellites of the first body.** They are stable down to
  C = 1.2175 (h1). Beyond that they are only slightly unstable, then stable again from h2 to h3. "Stable
  out to distances comparable with the separation of the two bodies."
- C = 1.2175 is far below C = 4 (where the regions around the two bodies connect) and C = 3.4568 (where
  they connect to infinity). "The Jacobi integral alone is quite insufficient to study satellite
  stability."
- Beyond h3 the types follow a periodic pattern 1-1-6-6, probably to infinity.
- h1 described twice is n'_2 5. h2 described twice is b_2 8. The critical-orbit type depends on the class
  in which the orbit is counted (h2 is type 5 in h and type 1 in b_2).

**Class i (direct satellites):** i1 C = 3.7388, x0 = -0.1813, type 1 (the stability limit); i2 3.859 /
-0.065 / 1; i3 3.858 / -0.060 / 6; i14 3.5664 / -0.2099 / 1; i15 = gamma_2 1, 3.4473 / -0.3228 / 5.
Retrograde satellites are much more stable than direct ones (Jackson 1913).

## 1. Content (READ, from the OCR and the English digest)

- **Method:**
  - Stability of the symmetric simple-periodic orbits of paper I (Hénon 1965a), from the Poincaré map T
    on y = 0 at fixed C.
  - Stability index a = (trace of the monodromy matrix)/2; the orbit is stable when |a| < 1.
  - Critical orbits (a = +-1) are classified into types 1-6, and the intersections with double-periodic
    classes are identified.
- Classes a, b, c (around the Lagrange points) are almost always strongly unstable.
- f and h are the retrograde satellites of each body. g and i are the direct satellites.
- New classes are found (j, k, l, m, n and others), with tables.
- **Application:** habitability of planets in double-star systems. The approximately circular orbits
  around one star, or around both, are stable over wide regions.
- **Typos in the original, as found by the digest agent:**
  - eq. (12) has a spurious factor of 2;
  - the Fig. 16 caption says "13" where l2 is meant;
  - running-head "PROLBEME".

## 2. Relevance

- `#944`/`#956` symmetric-family controls at mu = 1/2. These are the only published family h numbers at
  mu = 0.5 that the corpus holds.
- The type 1-6 critical-orbit taxonomy is the one used by Hénon & Guyot 1970 (held) and Bruno-Varin.

## 3. The translation request

- The owner asked for an English translation PDF.
- The translator subagent produced a complete-data English digest instead. Its reason: a full rendering
  would reproduce the whole copyrighted article.
- The same happened for Varin KIAM 16/2008.
- This has been reported to the lead for the owner's decision.

## 4. Citation mining

- Darwin 1911, Jackson 1913, Bartlett 1964 (on the wanted list) and Hénon & Heiles 1964 (held). Paper I
  is Hénon 1965a (on the wanted list).
- No new rows are added.
