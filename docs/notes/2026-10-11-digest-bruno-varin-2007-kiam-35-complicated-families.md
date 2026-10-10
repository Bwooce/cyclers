# Digest: Bruno & Varin 2007, "Complicated families of periodic solutions of the restricted problem" (KIAM Preprint 35/2007) (#960, #971 R14)

A. D. Bruno and V. P. Varin (Keldysh Institute of Applied Mathematics RAS), "Сложные семейства периодических
решений ограниченной задачи" ("Slozhnye semeistva periodicheskikh reshenii ogranichennoi zadachi"; library
English title "Complicated Families of Periodic Solutions of the Restricted Problem"), Preprint of KIAM No. 35,
Moscow, 2007. In Russian. UDK 521.1+531.314. Supported by RFBR grant 05-01-00050.
- **Identifiers:** no DOI (Crossref query 2026-10-11 found none; KIAM preprint DOIs start much later). Library id
  2007-35 (pid 2857). ISSN 2071-2898 (print), 2071-2901 (online), on the cover.
- **Library page:** https://library.keldysh.ru/preprint.asp?lg=e&id=2007-35 (HTTP 200). Authors "Bruno A. D., Varin
  V. P."; language Russian; **18 pages**; research direction "Mathematical problems and theory of numerical methods";
  CC BY 4.0 badge (link to creativecommons.org/licenses/by/4.0). The PDF cover also says "Статья доступна по лицензии
  Creative Commons Attribution 4.0 International".
- **File fetched:** https://keldysh.ru/papers/2007/source/prep2007_35.pdf, HTTP 200, application/pdf, 316 619 bytes,
  a real PDF 1.6 (Adobe Acrobat Pro 2017), 19 PDF pages = library cover + preprint pp.1-18.
  **md5 ef50600a0c6b4015417f2640d9f2ce65.**
- **Text layer:** unusable. The body is set in legacy Type 3 Cyrillic fonts (pdftotext gives "slovnye semejstwa
  periodi~eskih re{enij"). I made a Russian OCR sidecar with `ocrmypdf --force-ocr -l rus+eng --sidecar`
  (`prep2007_35-ocr-rus.txt`, 50 405 bytes). The OCR PDF copy (`prep2007_35-ocr.pdf`) is a by-product. File the
  original PDF, as was done for 34/2007.
- **Proposed corpus filenames:**
  - `cyclers_pdf/papers/bruno-varin-2007c-complicated-families-periodic-solutions-restricted-problem-kiam-preprint-35-2007-russian.pdf`
    (the original PDF; the stems `bruno-varin-2007-` and `bruno-varin-2007b-` are taken);
  - the same stem with `-ocr-rus.txt` for the sidecar;
  - the same stem with `-tables.yaml` for the transcription (or `data/sources/bruno-varin-2007c-kiam-35-tables.yaml`),
    and `-compare_ssr2008.py/.out` for the comparison script.
- **How I read it:** the whole OCR sidecar; the page images (200-300 dpi) of secs. 6.2-6.5 (pp.5-7), Tables 2-3
  (pp.17-18) and the reference list (pp.15-16). A contact sheet of all 19 pages at 50 dpi shows **no figures**: every
  page is text or a table. I compared it with the held JAMM 71:933 text layer (secs. 6-7), checked the four JAMM differences in sec. 0 on JAMM page images (170 dpi, pp.949, 951, 953), and compared JAMM Table 2
  (p.946) and Table 3 (p.952) on page images.

## 0. Verdict

**This is the Russian preprint of JAMM 71:933 secs. 6-7, with the full 51-item reference list. The JAMM digest's
guess ("likely but unverified") is now confirmed.**
- It continues KIAM 34/2007. Its sections are numbered 6 and 7, it opens at "6. Семейство i", and its
  bibliography is JAMM's list item for item (51 entries, same order).
- **Content:** sec. 6, family i (6.1 generating family, pieces K1-K20; 6.2-6.3 the closed families i1-i4; 6.4 the
  Hypothesis and Table 3; 6.5 the external annulus); sec. 7, horseshoes and tadpoles (eqs. 7.1-7.16 and the closing
  Remark). Figs. 8-16 are cited but **not printed**. JAMM prints them.

**For #971 R14 it adds no new numbers.**
- The cascade numbers are the same four pairs (mu'_k, mu''_k), k = 1-4, at 4-6 digits. Varin 2008 (KIAM 16/2008)
  and Varin & Bruno 2009 SSR 43(3):253 Table 2 remain the numeric source, at 9 digits
  (`data/sources/varin-bruno-2009-closed-families-tables.yaml`).
- There is no construction, proof or mechanism beyond JAMM sec. 6. The cascade is stated as a **Hypothesis** (sec.
  6.4), with "empirical asymptotics" k^(-8/3) in Table 3. So the R14 wording "predicted by Varin & Bruno's empirical
  asymptotics" stands.
- Sec. 6.5 attributes the Sun-Neptune (mu = 5.178e-5) computation to its ref. [9], Voyatzis & Kotoulas 2005 PSS
  (held). The last of those families "is closed ... the value of mu is not small and already sufficient for
  self-bifurcation (as for family i)". This is the R14 attribution already recorded in the Voyatzis-Kotoulas digest.

**What it does add (all read on the page images):**
1. **The JAMM Table 3 misprint is the journal's, not the authors'.** The preprint prints mu'_4 = 9.543·10^-5
   (p.18). JAMM prints 9.543·10^-4.
2. **JAMM sec. 6.5 has a translation error.** The preprint says family i covers the near-circular direct orbits
   "с a < 1, т.е. из внутреннего кольца" ("with a < 1, i.e. from the inner annulus", p.6). The JAMM English says
   "with a > 1, that is, from an internal annulus" (p.953, page image). a < 1 is right: the sentence contrasts the inner
   annulus with the external one.
3. **JAMM sec. 6.2 drops two mu values.** The preprint shows i1 in Fig. 11a at mu = 5e-3, mu_M = 1.2155e-2, 2.3e-2
   and 3.5e-2. JAMM p.951 (page image) lists only mu_M and 2.3e-2 for the closed curves.
4. **JAMM sec. 6.1 garbles a citation.** The preprint reads "([11], § 2; [20], п. 10.2.7; [47], п. 1.1)"; JAMM p.949 (page image) reads
   "(Ref. 11, § 20, Section 10.2.7 and Ref. 47, Section 1.1)".
5. **Reference-list questions from the JAMM digest (sec. 1.3, item 6), now checked on the preprint image:**
   - [31] Kogan: the preprint prints **1988**, T. 26, N. 6, pp.813-818. So JAMM's "1998" is a journal slip.
   - [44] Hitzl & Hénon: the preprint prints Acta Astronaut. **1977, V. 4, N. 9-10, pp.1019-1039**, which matches the
     held file. So JAMM's "15(4):421-52" is a journal copy error.
   - [37] Hénon VI: the preprint also prints "Astron. and Astrophys. 1969. V. 1. N. 1. P. 24-36". The error is in
     the authors' own list. The held file is A&A 9:24 (1970).
6. **Table 2 last digits differ from JAMM in 33 cells (sec. 2).** Neither print is wrong. The preprint mostly
   truncates the 5-digit values to 4 decimals. JAMM rounds about half of those cells instead. Use SSR 2008 Table 2
   for 5-digit values.
7. **Page count.** The library and the PDF give **18 pp.** The wanted row (from KIAM 51/2007 ref. [8]) says 27 pp.
   Both are recorded here; I do not choose. A possible source of "27": KIAM 36/2006 is 27 pp. (35/2007 ref. [47])
   and is cited next to it. This is a guess.

**Proposals (lead decides):**
- File the PDF, the OCR sidecar and the YAML under the names above. Add a CORPUS_INDEX row.
- **Wanted list row 16:** remove (received). Note "18 pp." against the "27 pp." in the row.
- **JAMM digest:** change "likely but unverified" to "confirmed"; add the four JAMM slips above (items 1-4) and the
  three reference-list answers (item 5).
- **R14:** no change to the plan or the likelihood. Cite this preprint with JAMM only for priority of wording.
- **Hénon VI [37]:** the corpus already holds A&A 9:24 (1970). Cite that form, not the Bruno-Varin form.
- **Full English translation made** (34/2007 conventions, tectonic, 13 pp.):
  `bruno-varin-2007c-complicated-families-periodic-solutions-restricted-problem-kiam-preprint-35-2007-russian-en-translation.tex/.pdf`.
  Sec. 6-7 text, references and Tables 2-3 (tables generated from the YAML). Equations were read on 130 dpi page
  images. Translator notes carry the JAMM slips (items 1-4), the reference-list errors, and two preprint points:
  - K18 says orbit 36 "coincides with orbit 20". By Table 2 (rows 22 and 36 share C, a~, e~) and the pattern of K4
    and K10 it should be orbit 22. The JAMM text layer prints "identical to orbit 22". Probable preprint misprint.
  - "rho_2 = z_2^2 + z_3^2" (p.8) is as printed; the definitions give 2 rho_2. JAMM text layer prints the same.
  Propose filing it beside the PDF.

## 1. Content

**Sec. 6.1, generating family i (mu = 0), pp.3-5.** Family i starts as direct circular orbits of infinitely small
radius about P1. It is made of pieces K1-K20 (Id, E±_{(p+1)/p}, C_{k,k+1} + l B1, C_{k,k+1}), grouped in cycles: K1-K4,
K5-K10, K11-K18, and K19-K20 begin the fourth. Each cycle starts and ends on the same orbit of a piece Id_p. Fig. 8
(not printed) shows the characteristics in (a~, e~); Fig. 9 (not printed) shows the zig-zag counts. The text matches
JAMM pp.949-950 piece by piece, except the citation in item 4 above.

**Sec. 6.2-6.3, the cascade, pp.5-6 (page image).**
- i1 is born at mu = mu'_1 ≈ 4.1313·10^-3, where two parts of the right characteristic meet (Fig. 10 d). It shrinks
  to a single orbit at mu''_1 ≈ 3.66863·10^-2.
- i2: mu'_2 ≈ 6.61705·10^-4 to mu''_2 ≈ 5.27272·10^-3. i3: mu'_3 ≈ 2.15292·10^-4 to mu''_3 ≈ 1.88241·10^-3.
  i4: mu'_4 ≈ 9.54305·10^-5 to mu''_4 ≈ 8.86552·10^-4.
- Another bifurcation happens for mu in (mu'_1, 5·10^-3) between the unclosed characteristics (p.6).

**Sec. 6.4.** The Hypothesis (two sequences mu'_k < mu''_k decreasing to 0; closed families i_k exist only on
[mu'_k, mu''_k]); Table 3. At mu_J, i2 and i3 were computed in [10, 11] (Colombo & Franklin 1968; Bruno KIAM
66/1993); related closed locally multiple families in [12] (KIAM 67/1993). At mu_M ≈ 1.2155092·10^-2 the closed family
i1 exists and was "not indicated" by Broucke 1968 [15].

**Sec. 6.5.** Direct near-circular orbits of the outer annulus lie on a countable set of families. At mu = 5.178·10^-5
families containing pieces Id_p = {p/(p-1) > N > (p+1)/p}, p = -2, ..., -7, were computed [9]. Their characteristics
stay near the generating ones, except the last, which is closed.

**Sec. 7, horseshoes and tadpoles, pp.7-13.** Same as JAMM sec. 7: local Poincaré coordinates near the circle of
fixed points, the reduced system (7.8)-(7.9), level curves of H (Fig. 12), the period function T~(v0) (eq. 7.14,
Fig. 14), the bifurcation condition T~(v0) = n sqrt(mu) (eq. 7.16, page image; with eq. 7.15 n = T~(v0)/sqrt(mu)), and the global structure of FH and FT families
(Figs. 15-16). The tadpole limit T~ → 2/(3 sqrt 3) = 0.385 is printed on p.10. The JAMM digest already checked it.
The closing Remark (p.12-13) says these families have no generating solutions in the main limit problem.

## 2. Tables (transcribed; `bruno-varin-2007c-kiam-35-tables.yaml`)

- **Table 2** (39 critical orbits of generating family i, mu = 0; columns k, T/(2π), C, Tr, a~(0), e~(0), a~(T/2),
  e~(T/2)). **Table 3** (mu'_k, mu''_k, their ratios, k^(-8/3), k = 1-4). The sec. 6.2-6.3 text values are listed beside
  Table 3.
- **Witnesses:** the page image decides. The ocrmypdf sidecar contains every numeric string of 3 or more characters
  in Table 2 (172 of 172; `compare_ssr2008.out`). It drops or garbles only symbols: the ∓ on row 18_1, the Tr "2" of
  rows 22 and 28, the infinity signs and the subscripts on k.
- **Conventions** (yaml `conventions:`): mu = 0, so the (1-mu) period factor and the +mu Jacobi shift of the SSR 2009
  family-h tables do not arise. The header prints T/(2π) with T the full period. C = -2H = 1/a + 2 sqrt(a) on the
  circular rows.

**Comparison with JAMM Table 2 and SSR 2008 Table 2** (`compare_ssr2008.py/.out`; SSR values from
`data/sources/bruno-varin-2008-ssr-tables.yaml`):
- JAMM differs from the preprint in 33 cells, each by one unit in the last digit. In 32 of them the preprint value
  is the SSR 5-digit value truncated, and the JAMM value is the same number rounded. Examples: row 9 C 3.0925 / 3.0926
  (SSR 3.092586); a~ on rows 21, 22, 38, 39_2: 0.7990, 0.8254, 0.8617 / 0.7991, 0.8255, 0.8618 (SSR 0.79906, 0.82548,
  0.86177). The 33rd is row 14_2 T~: 1.973 / 1.974. SSR prints 1.97277 (a probable misprint per the SSR yaml note;
  KIAM 36 and 51 print 1.97370).
- The rule is not uniform in either print. Over the 224 numeric cells, the preprint truncates 62, rounds 20, and has
  131 where both agree. JAMM truncates 30 and rounds 51. So "the preprint truncates, JAMM rounds" is a tendency, not
  a convention to rely on.
- 11 preprint cells match neither rule against SSR. JAMM prints the same value in each of them. They are the cells
  where SSR 2008 changed a value from KIAM 36/51: row 1 C (an SSR misprint per the SSR yaml), rows 25_2 T~, a~(T/2),
  e~(T/2) (KIAM 36 prints 3.08850, 0.75483, 1.92141), row 18_1 a~(T/2), e~(T/2) (23.872 against SSR 23.93659), row 26_2
  e~(0), and rows 32_3, 34_3 a~(T/2), e~(T/2). The preprint follows KIAM 36/51 in these. They are SSR-side questions,
  not preprint misprints. Row 18_1 is not explained by any note I found.
- **Table 3** agrees with JAMM except row 4 mu'_k (sec. 0, item 1). Its ratios check against the printed mantissas
  (JAMM digest checks.out).

## 3. Citation mining (51 references)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i` over CORPUS_INDEX.md (both run 2026-10-11; they agree). The list is JAMM's, so the JAMM digest's sec. 3
applies. Changes since that digest:
- **Now held:** Perko 1982 I and II [39, 40] and Perko 1983 [41] (`perko-1982a-...`, `perko-1982b-...`,
  `perko-1983-...`); Voyatzis & Kotoulas 2005 PSS [9]; Colombo, Franklin & Munford 1968 [10]; Hénon 1970 VI [37]
  (as A&A 9:24).
- **Still wanted, already on rows:** [7] Kotoulas & Voyatzis 2004, [8] Voyatzis, Kotoulas & Hadjidemetriou 2005,
  [16] Papadakis & Goudas 2006, [6] Franklin 1970: row 52. [31]-[35] Kogan, Lidov & Vashkov'yak, Benest 1976, and
  [4], [48]-[51] Llibre & Ollé, Barrabés & Mikkola, Schanzle, Taylor, Henrard: row 62. [11]-[14] Bruno KIAM 66, 67,
  68/1993 and 93/1996: row 51.
- **Not held, not on a row:** [3] Bruno 2006, "Периодические решения системы Гамильтона", Kosmich. Issled.
  44(3):258-271 (Russian pages as printed here). English form: Cosmic Research 44(3):245-257, doi
  10.1134/S0010952506030099 (Crossref, CONFIRMED 2026-10-11). It has the same title as Bruno's 1972 "Researches I"
  (see `digest-kiam18-1972-kiam93-1996-metadata.md`). Proposed as an addition to row 51, low priority. Also [26],
  [27], [28], [38] (power geometry, desingularisation): the JAMM digest's Tier D proposal stands.
- Page counts printed here for row 51: KIAM 66/1993 27 pp., 67/1993 29 pp., 68/1993 23 pp., 93/1996 32 pp.

*Filed as `cyclers_pdf/papers/bruno-varin-2007c-complicated-families-periodic-solutions-restricted-problem-kiam-preprint-35-2007-russian.pdf`. Check scripts, outputs, witness tables and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/bruno-varin-2007c-kiam-35-tables.yaml`.*
