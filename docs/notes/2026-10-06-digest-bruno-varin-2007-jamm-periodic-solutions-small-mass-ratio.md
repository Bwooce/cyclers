# Digest: Bruno & Varin 2007, "Periodic solutions of the restricted three-body problem for a small mass ratio" (JAMM 71:933-960) (#960 batch 29)

A. D. Bruno and V. P. Varin (Moscow), J. Appl. Math. Mech. 71(6):933-960 (2007; (c) 2008), doi
10.1016/j.jappmathmech.2007.12.012. English translation (by "E.L.S.") of Prikl. Mat. Mekh. 71(6):1034-1066 (2007).
Received 25 January 2007.
- Supplied by the owner as `5b3729d1-bruno2007.pdf`: 28 pp., Elsevier text layer, md5 93227424a1521aca7b91f3f4623fe369.
  Journal page = PDF page + 932.
- Filed as `cyclers_pdf/papers/bruno-varin-2007-periodic-solutions-restricted-three-body-problem-small-mass-ratio-jamm-71-933-doi-10.1016-j.jappmathmech.2007.12.012.pdf`.
- **How I read it:**
  - I read the whole text layer.
  - I read these on 200 dpi page images: p.939 (Fig. 2, Table 1), p.946 (Table 2), p.947, p.948 (Figs. 6-7 and
    the sec. 5.1 text), p.952 (Fig. 11, Table 3) and p.957 (eq. 7.16).
  - I compared it section by section with the held KIAM preprint 34/2007: the Russian OCR `.txt` and our English
    `.tex` translation.
  - Arithmetic checks: `cyclers_pdf/papers/<pdf stem>-checks.py`, output in `cyclers_pdf/papers/<pdf stem>-checks.out` (filed beside the PDF).
- Wanted-list row 18. In batch 29 the row was cut down to the follow-up preprint only (sec. 1.4). (Wanted-list row numbers in this digest are the batch-28 numbering; the list was renumbered in batch 29.)

## 0. Verdict

**This IS the journal form of KIAM preprint 34/2007, with two extra sections and the reference list.**
- Secs. 1-5 of the two are the same text, sentence for sentence. The JAMM English is a different translation
  from ours and loses a few sentences (sec. 1.3).
- The preprint stops at sec. 5.2. JAMM goes on with:
  - sec. 6: family i, with Table 2, Table 3 and Figs. 8-11;
  - sec. 7: horseshoe and tadpole orbits, with Figs. 12-16;
  - a 51-item reference list.
- **Row 18's gap is closed.** JAMM gives the reference list and all the figures that 34/2007 cites but does not
  print (Figs. 1-7).
- The preprint's citation numbers [1]-[46] map one to one onto JAMM refs 1-46. JAMM adds 47-51, which are cited only
  in secs. 6-7.

**Family h at finite mu: no new numbers.**
- The only family-h data at mu > 0 is **Fig. 7**, which shows the characteristics in (a~, e~) at mu = 0.1, 0.3 and
  0.5. It is a picture with no tables. It is the figure that was missing from 34/2007.
- It could be digitised as a rough shape check. It does not replace the absent Tables 4-8 of preprints 67 and 64
  (row 17 stays).
- Table 1 (generating family h, mu = 0) is the same as the preprint's Table 1, cell by cell.
- The only finite-mu numbers in the paper are the **family-i bifurcation values mu'_k and mu''_k** (sec. 6.2-6.3,
  Table 3). They agree with the held Varin 2008 (KIAM 16/2008) Table 2.
- **Table 3 has a misprint.** It prints mu'_4 = 9.543·10^-4. It should be 9.543·10^-5 (sec. 2.3).

**What it gives the project:**
- A citable, English, peer-reviewed source for:
  - the four mu -> 0 limit problems;
  - the generating families h and i (Tables 1-2);
  - the closed-family cascade i_k (Table 3).
- **The Earth-Moon closed family i_1** (p.952, Fig. 11a). At mu_M = 1.2155092e-2, a closed family i_1 exists
  "which was not indicated when calculating the family i for mu_M" by Broucke 1968 (held). So Broucke's TR 32-1168 is
  incomplete for family i at Earth-Moon mass ratio. Our use of Broucke as a completeness control should note this.
- A horseshoe and tadpole perturbation theory (sec. 7). It gives a closed-form condition for where tadpole families
  bifurcate (eq. 7.16).

**Catalogue implication (PROPOSAL only):** none directly. No catalogue row cites these families as cycler
sources. One proposal: if any `#944` X2 or `#956` R9 work uses Broucke 1968 as a family-i completeness reference at
Earth-Moon mu, add the i_1 caveat.

## 1. Q1: identity with KIAM preprint 34/2007

### 1.1 Section map

| Part | KIAM 34/2007 (23 pp.) | JAMM 71:933 (28 pp.) |
|---|---|---|
| Abstract | family h only ("changes little") | adds family i (infinitely many self-bifurcations, closed subfamilies) and the horseshoe/tadpole theory |
| 1.1 Statement | yes; has an extra sentence ("the aim of this article is to give a survey ...") | same, without that sentence |
| 1.2 Contents | announces secs. 5, 6, 7 and says 6.2-6.4 and 7 are new | same |
| 1.3 Properties | yes | same |
| 2.1-2.3 Limit problems | yes; Fig. 1-2 cited, not printed | same; **Fig. 1 (p.936) and Fig. 2 (p.939) printed** |
| 3 Main limit problem | yes; Fig. 3-4 cited | same; **Fig. 3 (p.940), Fig. 4 (p.944) printed** |
| 4 Generating families | yes; Fig. 5 cited | same; **Fig. 5 (p.945) printed**; Table 2 is placed on p.946 but belongs to sec. 6 |
| 5.1 Generating family h | Table 1 (p.22); Fig. 6 cited | Table 1 (p.939); **Fig. 6 (p.948) printed** |
| 5.2 h for mu -> 1/2 | Fig. 7 cited; ends "continued in the preprint 'Complex families ...'; references there" | **Fig. 7 (p.948) printed**; no "continued" line |
| 6 Family i (6.1-6.5) | absent | pp.949-953: Table 2 (generating family i, 39 critical orbits), Figs. 8-11, the hypothesis 6.4.1, Table 3, and sec. 6.5 (external annulus at mu = 5.178e-5) |
| 7 Horseshoes and tadpoles | absent | pp.953-959: eqs. 7.1-7.16, Figs. 12-16, and a closing Remark |
| References | absent | 51 items, pp.959-960 |
| Acknowledgement | RFBR 05-01-00050 | same, plus a dedication to Euler's 300th anniversary |

**Citation numbers align.** Some pairs where the preprint's bracket number sits in the same sentence as the
matching JAMM entry:
- "[14], §4 ... IR+J" = JAMM 14, Bruno KIAM 93/1996.
- "[45], §1; [22], §6" = JAMM 45 (Bruno-Varin KIAM 67/2005) and 22 (Bruno-Varin 2006 CMDA).
- "mu = 0.3, 0.4, 0.5 [46]" = JAMM 46, KIAM 64/2005.
- "Hill's problem [29]" = JAMM 29, Hill 1878.
- "Hénon [30, 37]" = JAMM 30 and 37, Hénon 1969 V and Hill VI.
- "Euler [1]" = JAMM 1, Euler 1772.

So the JAMM list is 34/2007's list. It is not merely a similar one.

**The "Complex families" follow-up preprint.** The 34/2007 tail says the continuation is in the preprint "Complex
families of periodic solutions of the restricted problem". JAMM secs. 6-7 are very likely that continuation in
journal form: 34/2007 sec. 1.2 already announces secs. 6 and 7. This is **likely but unverified**, because the
follow-up preprint was not seen. Proposal: drop that preprint's priority on the wanted list to "attribution only".

### 1.2 Table 1: cell-by-cell comparison

- I read JAMM Table 1 (p.939) on the page image. All 16 rows and 8 columns agree with the 34/2007 Table 1
  transcription in our `.tex` (which was itself read on the preprint's page image).
- Only typography differs. JAMM row 6, e~(T/2), prints "-1,370" with a decimal comma. It is -1.370 in the preprint.

### 1.3 Text differences (JAMM defects; prefer the preprint for secs. 1-5)

1. **Sec. 5.1, p.948: a dropped sentence changes the meaning.**
   - JAMM: "From orbit 6, the family h continues as the family A1 up to orbit 11. Here, the Jacobi constant C reaches
     a minimum in orbit 9."
   - The Russian preprint (OCR lines 928-930) says three things: C has its minimum on orbit 7; orbit 9 is a collision
     with P2 (the zero-mass body); C has its maximum on orbit 10.
   - **Table 1 supports the preprint:**
     - C(7) = -1.785103 is below C(6) = -1.711013 and C(8) = -1.491531.
     - Row 9 has a~(T/2) undefined and e~(T/2) = -+inf, like the P2-collision rows 2 and 16.
     - C(10) = 2.929162 is above C(11) = 2.929161.
   - I cannot tell whether the loss is in the PMM Russian or only in the translation. The PMM original is not held.
2. Sec. 5.1 also drops "Many details of this family, including orbits, were given ([45], §1; [22], §6)".
3. Sec. 2.3 drops "[Family c] is obtained from the family a by the mapping (2.9)".
4. Sec. 1.1 says "the body P1 executes a Kepler motion with respect to the body P1". The preprint has P2 relative to
   P1.
5. Sec. 4.4:
   - The heading reads "[Refs. 14, 20; Ch. 1]". The preprint has "[14], [20], Ch. 10".
   - Family b runs "from e~ = -1 to e~ = -1". The preprint has "from e~ = 1 to e~ = -1".
   - "E1/2(t)" should be E1/2(1).
6. Reference-list slips. These are checked against the held files or against the wanted list. They are not from
   memory.
   - [44] Hitzl & Hénon 1977b is printed "Acta Astronaut 1977;15(4):421-52". That copies [43]'s volume and pages. The
     held file is Acta Astronaut. 4:1019 (`hitzl-henon-1977b-...-acta-astronaut-4-1019-...`).
   - [37] Hénon, "Numerical exploration ... VI. Hill's case: nonperiodic orbits", is printed "A&A 1969;1(1):24-36".
     Wanted row 33 has A&A 9:24-36 (1970). Row 33's citation was "not checked". ADS should settle it before anyone
     cites either form.
   - [31] Kogan, "Distant satellite orbits ...", Kosmich. Issled., is printed "1998;26(6):813-8". Volume 26 of
     Kosmicheskie Issledovaniya is 1988, so 1998 is probably a typo. This is unverified; Kogan is not held.
7. Errata kept from the preprint (our `.tex` translator notes): (2.6) and (2.7) print the index i where the variables
   carry j; (3.4) prints m(2a)^-1/2. JAMM is garbled at these equations in the text layer, so I did not image-check
   whether JAMM corrected them.

## 2. Q2: family h, or any family, at finite mu

- **Family h, mu > 0:**
  - The only numbers are the mu values themselves: mu_J = 0.00095388; 0.1 and 0.2 [45]; 0.3, 0.4 and 0.5 [46];
    ~0.012 [15].
  - Fig. 7 (p.948, image checked) plots the characteristics at mu = 0.1 (solid), 0.3 (dashed) and 0.5 (dash-dot) for
    a~ in [-6, 3] and e~ in [-2, 2].
  - The text says there are "no new singularities" and "self-bifurcations do not occur". This matches the 64/2005
    conclusion recorded in the KIAM digest.
  - There are no initial conditions, no Jacobi constants and no periods at mu > 0. **Held-digest status is unchanged:**
    the 67/64 tables at mu = 0.1-0.5 are still not recovered (wanted row 17).
- **Family h, mu = 0 (Table 1, p.939, image checked):**
  - The same 16 critical orbits as 34/2007 Table 1 and as CMDA 2006 Tables 1-3 (held).
  - Our existing digests already quote it in full (`2026-10-06-digest-bruno-varin-kiam-preprints-2005-2010.md` sec. 3).
    It is not repeated here.
- **Family i at mu = 0 (Table 2, p.946, image checked):**
  - 39 critical orbits, k = 1 to 39_2, with T/(2 pi), C, Tr, a~(0), e~(0), a~(T/2), e~(T/2).
  - The subscript on k is the number of arc-solutions in the orbit.
  - Examples: k = 1: T~ = 1/2, C = 3.4668, a~(0) = -0.4807, e~ = 1. k = 18_1: T~ = 2, C = 1.4845, a~(0) = -0.6736,
    e~ = -+2, a~(T/2) = 23.872. k = 39_2: T~ = 4, C = 2.7447, a~(0) = -0.8618, e~(0) = 0.4786.
- **Family i, finite mu (secs. 6.2-6.3, p.950-952; Table 3, p.952, image checked):**

| k | mu'_k (text) | mu''_k (text) | Table 3 mu'_k | Table 3 mu''_k | Varin 2008 Table 2 (held) |
|---|---|---|---|---|---|
| 1 | ≈ 4.1313e-3 | ≈ 3.66863e-2 | 4.131e-3 | 3.669e-2 | 4.13129887e-3 / 3.66863029e-2 |
| 2 | ≈ 6.61705e-4 | ≈ 5.27272e-3 | 6.617e-4 | 5.273e-3 | 6.61705554e-4 / 5.27272358e-3 |
| 3 | ≈ 2.15292e-4 | ≈ 1.88241e-3 | 2.153e-4 | 1.882e-3 | 2.15292269e-4 / 1.88241384e-3 |
| 4 | ≈ 9.54305e-5 | ≈ 8.86552e-4 | **9.543e-4 (misprint)** | 8.866e-4 | 9.54304953e-5 / 8.86552296e-4 |

- **The Table 3 misprint at k = 4.** As printed, mu'_4 = 9.543e-4 would be larger than mu''_4 = 8.866e-4. That breaks
  the paper's own condition mu'_k < mu''_k. It would also give mu'_4/mu'_1 = 0.231, not the printed 0.023. The p.952
  text and Varin 2008 both give 9.54305e-5. That value gives 0.023, as printed (`cyclers_pdf/papers/<pdf stem>-checks.out`).
- The other ratios and k^(-8/3) (0.157, 0.053, 0.025) check, as printed.
- The text truncates rather than rounds: mu'_2 is 6.61705e-4 against Varin's 6.61705554e-4.
- These numbers come from the same authors' work, so the match is consistency, not independent confirmation.
- Other finite-mu statements:
  - The families i_2 and i_3 were computed at mu_J by Colombo & Franklin 1968 and by Bruno 1993 [10, 11].
  - The closed family i_1 exists at mu_M (Fig. 11a).
  - In sec. 6.5, families containing Id_p pieces for p = -2 to -7 were computed at mu = 5.178e-5. Only the shapes are
    described; there are no numbers.

### 2.1 Arithmetic checks (`cyclers_pdf/papers/<pdf stem>-checks.py`)

- **Table 1, t = 0.** At mu = 0, x1 = a~(2-|e~|), y2|y2| = e~/x1 and C = -y2^2 + 2 x1 y2 + 2/|x1|.
  - The printed C is reproduced from (a~(0), e~(0)) in all 16 rows. The largest difference is 1.6e-5 (row 9). That
    is within the 5-decimal rounding of a~ and e~.
  - The x1(0) computed this way matches the CMDA 2006 Table 1 x1(0) to ±1e-5 in every row.
  - Closed forms: 2^(2/3) = 1.58740, 4^(2/3) = 2.51984, 2^(-2/3) = 0.629961 and 4^(-2/3) = 0.396850, as printed.
- **Table 1, T/2 columns: not verified for rows 7, 8, 14 and 15.**
  - Rows 3, 4, 6, 10, 11 and 13 give C at T/2 within 1e-3 (the 3-decimal rounding).
  - Rows 7, 8 and 15 have |e~(T/2)| > 2, which is outside the strip (3.14).
  - Row 14, read as (3.12), gives C = -2.541 against the printed -2.141393.
  - So these cells are not plain (3.12) coordinates. Also, CMDA Table 3 gives a~*(T/2), e~*(T/2) for the same rows
    with different values: for example, row 14 is -0.19450 / 28.10268.
  - This is an open convention question. It is the same in both prints and does not affect the t = 0 data.
- **Table 2, circular rows on Id (e~ = 1).** Kepler gives T~ = 1/(N-1) with N = a^-1.5, and C = 1/a + 2 sqrt(a).
  - Rows 1, 2, 9, 10, 21, 22, 37 and 38 reproduce T~ = 1/2 to 4 within ±0.001 and C within ±1e-4.
- **Table 2, all other non-collision rows.** C from (a~(0), e~(0)) agrees within 2e-4, with one exception.
  - Row 26_2 differs by 0.003. Its e~(0) = -1.9985 sits next to the P1-collision line |e~| = 2, where C is very
    sensitive to e~. Matching the printed C needs e~ ≈ -1.9986, one or two units in the last printed digit. I do not
    count this as an error.
- **Sec. 7 check.** As v0 -> 1, T~(v0) -> 2/(3 sqrt 3) = 0.3849, and the paper prints ≈ 0.385.
  - This equals the classical small-amplitude long-period L4 libration period, 2 pi / sqrt(27 mu / 4), in time scaled
    by sqrt(mu). That is an independent check.
  - So the tadpole bifurcation condition T~ = n sqrt(mu) (7.16) needs n > 2/(3 sqrt(3 mu)), as printed on p.957.

## 3. Citation mining (51 references)

The checks use `ls cyclers_pdf/papers | grep -i` plus a grep of `docs/notes/CORPUS_INDEX.md`. Wanted-list rows are
numbered as in `2026-10-05-960-wanted-papers.md` today.

**HELD:**
- [2] Bruno 1994 de Gruyter (`bruno-1994-restricted-3-body-problem-...`). It was filed this batch, but **wanted row 19
  still lists it**, so row 19 is stale.
- [10]-[14]: none held (see below).
- [15] Broucke 1968 (`broucke-1968-...`).
- [17] Hénon 1965a (`henon-1965a-...`).
- [18] Hénon 1968 (`henon-1968-...`).
- [19] Bruno-Varin KIAM 10/2005 (`bruno-varin-2005a-...`, TeX source and translation).
- [20] Hénon 1997 (`henon-1997-...`).
- [21] Hénon 2001 (`henon-2001-...`).
- [22] Bruno-Varin 2006 (`bruno-varin-2006-...`).
- [23] Szebehely 1967 (`szebehely-1967-...`).
- [25] Bruno 1981 (`bruno-1981-...`).
- [30] Hénon 1969 V (`henon-1969-...`).
- [42] Perko 1981 (`perko-1981-...`).
- [43] Hitzl & Hénon 1977 (`hitzl-henon-1977-...`).
- [44] Hitzl & Hénon 1977b (`hitzl-henon-1977b-...`).
- [45] and [46] Bruno-Varin KIAM 67 and 64/2005: held as TeX source only (`bruno-varin-2005b-...`, `-2005c-...`). The
  tables are absent; row 17 stays.

**Already on the wanted list:**
- [37] Hénon VI: row 33. Note the year and volume conflict in sec. 1.3.
- [10] Colombo & Franklin 1968 and [7] Kotoulas & Voyatzis 2004: row 62.
- [11] Bruno KIAM 66/1993, [14] KIAM 93/1996 and [47] Bruno-Varin KIAM 36/2006: row 61.
- [1]-[46] as a list: row 18. JAMM itself closes row 18.

**Not held, new candidates (none is on the wanted list):**
- **Perko 1983**, "Periodic solutions of the restricted problem that are analytic continuations of periodic solutions
  of Hill's problem for small mu > 0", Celest. Mech. 30:115-132 [41]. This is the Hill-to-restricted continuation
  theorem behind family-f/g persistence (DRO-type, `#945` R2). **The best new candidate; proposed as its own row.**
- Perko 1982, "Families of symmetric periodic solutions of Hill's problem I, II", Amer. J. Math. 104:321-352 and
  353-397 [39, 40]. Proposed as one row with Perko 1983.
- Bruno KIAM 67/1993 and 68/1993 [12, 13] (double and multiple periodic solutions, Sun-Jupiter). Proposed addition to
  row 61.
- Bruno & Petrovich, KIAM 53/2006, "Desingularization of the restricted three-body problem" [28]. Also Bruno 2006,
  Cosmic Research 44(3):245-257 [3]; Bruno 2000, Power Geometry (Elsevier) [26]; Bruno 1999, J. Math. Sci.
  95(5):2483-2512 [27]; Bruno KIAM 91/1978 [24] (the content is held as Bruno 1981); Bruno 1994 in Seimenis (ed.),
  Hamiltonian Mechanics, Plenum, pp.43-49 [38]. Low priority; proposed as one Tier D row.
- Quasi-satellite and distant-retrograde background [31]-[35]:
  - Kogan 1988(?), Kosmich. Issled. 26(6):813-818;
  - Lidov & Vashkov'yak 1990 (Kazan), 1993 (Kosmich. Issled. 31(2):75-99) and 1994 (Pis'ma Astron. Zh. 20(3):229-240);
  - Benest 1976, Celest. Mech. 13:203-215. The one CORPUS_INDEX hit for "benest" is in an unrelated row, so Benest is
    not held.
  - Relevant to `#945` R2 DRO controls. Proposed as one row; Benest 1976 is the most accessible.
- Horseshoe and co-orbital [4, 5, 48-51]:
  - Llibre & Ollé 2001, A&A 378:1087-1099;
  - Llibre & Ollé 2004, in "New Advances in Celestial Mechanics and Hamiltonian Systems", Kluwer, pp.137-152;
  - Barrabés & Mikkola 2005, A&A 432:1115-1129;
  - Schanzle 1967, AJ 72:149-157;
  - Taylor 1981, A&A 103:288-294;
  - Henrard 2002, CMDA 83:291-302.
  - None is held. Proposed as one Tier D row (co-orbital family structure). It is not a cycler source.
- [8], [9] Voyatzis, Kotoulas & Hadjidemetriou 2005, CMDA 91:191-202, and Voyatzis & Kotoulas 2005, Planet. Space Sci.
  53:1189-1199. These are **not on row 62**, which has Kotoulas & Voyatzis 2004 only. They were named in the CMDA 2006
  digest's not-held list but never added. Proposed addition to row 62.
- [16] Papadakis & Goudas 2006, Astrophys. Space Sci. 305:99-124 (the mu = 0.4 survey). Proposed addition to row 62.
- [6] Franklin & Colombo 1970, Icarus 12:338-347. Proposed addition to row 62.
- Not proposed: [1] Euler 1772; [29] Hill 1878 (classical; background only); [36] Wintner 1941 (textbook).
