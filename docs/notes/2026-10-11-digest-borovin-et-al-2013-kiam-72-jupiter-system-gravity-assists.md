# Digest: Borovin, Golubev, Grushevskii, Koryanov & Tuchin 2013, "Flights in the Jupiter System Using Gravity Assists near Galilean Satellites" (KIAM Preprint 72/2013) (#960, #943)

G. K. Borovin, Yu. F. Golubev, A. V. Grushevskii, V. V. Koryanov, A. G. Tuchin (KIAM RAS), "Полёты в системе Юпитера с
использованием гравитационных манёвров около галилеевых спутников", Preprint of KIAM No. 72, Moscow, 2013, 32 pp.
In Russian. UDK 629.78. Supported by RAS Presidium Programme 22.
- **Identifiers:** no DOI (Crossref query 2026-10-11 found none). Library id 2013-72 (pid 3652).
- **Library page:** https://library.keldysh.ru/preprint.asp?lg=e&id=2013-72 (HTTP 200). English title "Flights in
  the Jupiter System Using Gravity Assists near Galilean Satellites"; keywords gravity assist, Jovian system,
  Tisserand graph, radiation dose; 32 pages; "Theoretical and applied problems of mechanics". The page carries the
  CC BY 4.0 badge. **The PDF cover has no licence line** (unlike the 2007 covers). Record both.
- **File fetched:** https://keldysh.ru/papers/2013/prep2013_72.pdf, HTTP 200, application/pdf, 1 914 160 bytes, a
  real PDF 1.6 (Word 2010 / Acrobat 9.5.5), 33 PDF pages = library cover + 32 pp. **md5 f097c337598b8972384e030496388a71.**
  Clean Unicode text layer (no OCR needed).
- **Proposed filename:** `cyclers_pdf/papers/borovin-golubev-grushevskii-koryanov-tuchin-2013-flights-jupiter-system-gravity-assists-galilean-satellites-kiam-preprint-72-2013-russian.pdf`.
- **How I read it:** the whole text layer; Fig. 5 (p.18) and the cover on page images.

## 0. Verdict

**The preprint form of Golubev, Grushevskii, Koryanov & Tuchin 2014, JCSSI 53:445 (held, digested). #943 collision
check: none.**
- Same mission frame (Laplace-P Ganymede lander), same phases (JOI, period-reduction "debut", v_inf reduction,
  endgame, GOI), same Tisserand-Poincaré and "reflection / re-reflection" (billiard) language, and the same reading
  of the JUICE event table: **G5-C11-G12, G12-C24-G25, G25-C29-G30** as v_inf-reducing links (p.16).
- **No repeating G-C sequence.** The G-C-G links are one-shot: a cross ("косой") flyby K12 of the second moon puts
  the craft back on the first moon with a lower v_inf. The chain grammar is printed as
  R11 ... R1n K12 K21 R1(n+1) ... (eq. 22) and R11 ... R1n K12 R21 ... R2k K21 R1(n+1) ... (eq. 23).
  Neither repeats.
- **No v_inf values for G or C, no periods for a G-C leg.** The debut is a Ganymede-resonant sequence of periods
  6→5→5→4→3→2.5→2 Ganymede periods (p.11; "approximately 7.155 days" per Ganymede period, p.4; both read on page images). Fig. 5 (isoinfine bands on Tisserand
  graphs for Ganymede and Callisto) has colour bands at "a fixed step" in v_inf but **no v_inf labels** (page image).
  So nothing can be compared with gc-1 (37.57 d, v_inf G/C 2.40/1.81 km/s) or gc-2 (3.62/3.04).
- **Collision with gc-1 / gc-2: none.** "Related" only at the concept level (two-moon G-C flybys), as already
  recorded for JCSSI 2014 and ISSFD 2017. Consistent with both digests.

**What it gives the project:** nothing new beyond the JCSSI paper. Two small items:
- Table 1 (p.7; values from the clean text layer, not image-checked): dV_max and chi_mod = dV_max / V_pl for planets, dwarf planets and moons (e.g. Io 1.809, Europa
  1.433, Ganymede 1.949, Callisto 1.725 km/s; Titan 1.867; the Moon 1.680). This is the clean text-layer version of
  the table that the JCSSI digest calls "scrambled". Background only; not a golden source.
- The radiation-dose argument (pp.21-23) for keeping a Laplace-P tour to Ganymede and Callisto only (150-250 krad
  behind 10 mm Al). Context for any #943 G-C mission framing.

**Proposals:** file it under the name above with a CORPUS_INDEX row ("preprint form of JCSSI 53:445; no #943
collision"). Remove Borovin et al. 2013 from wanted row 50.

## 1. Citation mining (39 references)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i` over CORPUS_INDEX.md (both run 2026-10-11; they agree).
- Held (file checked): Minovitch 1963 TR 32-464 [2]; Minovitch 1972 JSR 9(10) [3]; Strange & Longuski 2002 [31];
  Strange, Russell & Buffington AAS 07-277 [22]; Szebehely 1967 (English original of [26]).
- Endgame I and II [20, 21] are held only in their AAS 09-224 / 09-227 conference forms
  (`campagnola-russell-2009-endgame-partA-...`, `-partB-...`), not as the JGCD 33(2) papers cited here.
- Not held: Uphoff, Roberts & Friedman 1976 JSR 13(6) [29] (Galileo-era tour design; row 39 is the Jovian background
  row; not checked whether it names this paper). Low priority.
- On wanted row 50 already: Woolley 2010 thesis [17]; Boutonnet & Schoenmaekers AAS 12-207 [4]; Campagnola et al.
  AAS 12-185 [24]; Labunsky, Papkov & Sukhanov 1998 [30]. Campagnola, Skerritt & Russell AAS 11-425 [23] has its own
  Tier B row.
- New, low priority, not proposed: KIAM internal reports [15, 16]; Ovchinnikov, Trofimov & Shirobokov, KIAM preprint
  9/2012 [18]; Myatt et al. 2004 Ariadna report [19]; Golubev, Grushevskii & Khairullin KIAM preprint 78/1993 and Kosmich. Issled. 34(2) 1996 [34, 35];
  Grushevskii KIAM preprint 76/2007 [36]; Podzolko & Getselev 2013 radiation [38]. None bears on a cycler.

*Filed as `cyclers_pdf/papers/borovin-golubev-grushevskii-koryanov-tuchin-2013-flights-jupiter-system-gravity-assists-galilean-satellites-kiam-preprint-72-2013-russian.pdf`.*
