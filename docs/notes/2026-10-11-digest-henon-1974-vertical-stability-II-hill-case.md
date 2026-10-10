# Digest: Hénon 1974, "Vertical Stability of Periodic Orbits in the Restricted Problem. II. Hill's Case" (#960 batch 42, #1017)

M. Hénon (Observatoire de Nice), Astronomy & Astrophysics 30:317-321 (1974), ADS 1974A&A....30..317H.
Received 26 October 1973. No DOI (Crossref has no record; checked 2026-10-11).
- Proposed filename: `henon-1974-vertical-stability-periodic-orbits-restricted-problem-II-hill-case-aa-30-317-ads-1974AA-30-317H.pdf`.
  - The file to file is an `ocrmypdf --force-ocr -l eng` copy, md5 1449147d7e1291f8a2dc0ecd8a76e7d2.
  - The source is an image-only ADS JBIG2 scan, 5 pp., md5 8822c83329f1e01b50d2a576b9861c15 (owner upload).
  - Page images survived the OCR. Dark-pixel fractions at 100 dpi, original against OCR copy: p.1 0.062/0.068,
    p.3 0.054/0.059, p.5 0.026/0.028. The copy is slightly darker, with no blank pages.
- **How it was read:** all five pages were read on page images. Equations and tables were read on 400 dpi crops.
  The OCR text layer was used only to find things. It garbles the g′₂ names ("£53v", "253") and drops signs.
- **Companion `...-tables.yaml`** (`henon-1974-hill-vertical-stability.yaml`): eq. (7) and Tables 1-6, 52 rows.
  - The format follows `data/sources/henon-2003-hill-families-tables.yaml`.
  - Three witnesses: the image, tesseract and macOS Vision. The image decides.
  - Comparison file: `witness-comparison.tsv`, 138 numeric cells.
- Checks: `recompute_av.py/.out`, `crosscheck.py/.out`, `ejection_check.py/.out`, `gprime2_av.py/.out` and `hill.py`.
  All values in these files are ours.

## 0. Verdict

**This is the vertical-stability companion to H5 (Hénon 1969) for the simple Hill families a, c, f, g and g′.**
- It prints the vertical index a_v for these families:
  - all their vertical critical orbits (Table 1);
  - exact a_v at 20 Γ values (Tables 2, 3, 5, 6);
  - the L2 start value (eq. 7);
  - the asymptotic laws for Γ → ±∞.
- We reproduced every printed a_v to the printed digits (sec. 3): eq. (7), all 20 exact values and all 10 critical orbits.
  Three cells in the asymptotic-formula columns do not reproduce. They are printed slips (sec. 1).

**For #974 it is partial input only.**
- #974 needs the Hénon 2003 families Ha-Hf, Hg and g3, plus the asymmetric second-species parents (R15).
- This paper gives no a_v for any of those. It is the "printed a_v of the Hill families" that #1017 and #994 asked
  for, but only for the parent families.
- **Overlap with Hénon 2003: no shared family. One shared orbit.**
  - 1974 Table 1 g′2 (Γ 4.271428, ξ0 0.073664) is the maximum-Γ row of 2003 Table VI, family Ha.
  - Hénon 2003 calls that row "critical orbit of g′ described twice".
  - The 2003 Hc intersection with g′ described three times (Γ 4.279538, a = −1/2) lies on the g′ segment that this
    paper finds three-dimensionally stable (g′2v 4.280603 to g′2 4.271428). This paper does not print it.
  - Our a_v for g′2 is −0.944739. So the Ha maximum (g′2 described twice) has a_v = 2a_v² − 1 = 0.785065, and it is
    vertically stable. This value is ours, not printed.
- **Vertical critical orbits found** (Table 1; Type from Paper I eq. 22: A_v or D_v means a_v = −1, B_v or C_v means
  a_v = +1):

| Orbit | ξ0 | ξ1 | Γ | η̇0 η̇1 | Type | a_v |
|---|---|---|---|---|---|---|
| a1v | 0.581265 | 0.774657 | 4.005313 | + − | C_v | +1 |
| a2v | 0.212661 | 0.949136 | 1.228063 | + − | B_v | +1 |
| a3v | 0.092988 | 1.218163 | −0.029389 | + − | A_v | −1 |
| g1v | 0.310842 | −0.310842 | 3.057455 | + − | A_v−D_v | −1 (tangent) |
| g2v | 0.164715 | −0.164715 | 1.383094 | + − | C_v | +1 |
| g′1v | 0.086462 | −0.570854 | 4.285184 | + − | A_v | −1 |
| g′2v | 0.082216 | −0.576960 | 4.280603 | + − | D_v | −1 |
| g′₂3v | 0.004399 | −0.464697 | 3.390160 | − − | C_v | +1 |
| g′₂4v | 0.285816 | −0.114197 | 0.477157 | − − | B_v | +1 |
| g′₂5v | 0.408116 | −0.059949 | −0.219529 | − − | D_v | −1 |

Family c has the same critical orbits, mirrored through M2 (p.318). Family f has none.

- **Stability results (text and figures, pp.317-319):**
  - **Family a (and c):** vertically stable from L2 (a_v = 0.976516) to a1v, unstable to a2v, stable to a3v, then
    unstable to the end (a_v → −∞). The family is never three-dimensionally stable, because it is horizontally
    unstable throughout (H5).
  - **Family f (retrograde, the DRO family):** "all orbits of family f are vertically stable, and also
    three-dimensionally stable" (p.318). a_v → +1 at both ends. The minimum is below 0.1 near Γ ≈ −0.3, read from
    Fig. 3, so it is approximate.
  - **Family g (direct):** vertically stable until it touches −1 at g1v. It is unstable after g2v and stays so.
    It is three-dimensionally stable only "between the origin and g1".
  - **Family g′:** vertically stable at first, then unstable on g′1v-g′2v (very short) and on g′₂3v-g′₂4v, and
    unstable for good after g′₂5v (a_v → −∞). It is three-dimensionally stable only "between g1 and g′1v and
    between g′2v and g′2".
- **Use:**
  - Gates on any 3D claim about lifts of the simple Hill families (`#994`, R15 controls). Lifts of these parents
    are known-class under `#438`/`#444`.
  - The 10 vertical critical orbits are the bifurcation points to the 3D Hill families.
  - f (DRO) is 3D-stable everywhere in Hill's limit.

## 1. Content

**Equations** (sec. 2, p.317, read on the page image).
- Hill's limit of Paper I's equations, by eq. (1) x = 1 − μ + μ^{1/3}ξ, y = μ^{1/3}η, z = μ^{1/3}ζ and eq. (2).
- Eq. (3): ξ̈ = 2η̇ + Ω′_ξ, η̈ = −2ξ̇ + Ω′_η, ζ̈ = Ω′_ζ.
- Eq. (4): Ω′ = 3/2 ξ² − 1/2 ζ² + 1/ρ.
- Eq. (5): Γ = 2Ω′ − (ξ̇² + η̇² + ζ̇²).
- Eq. (6), the vertical variational equation: d²Δζ/dt² = −[1 + (ξ² + η²)^{−3/2}] Δζ.
  - **Note:** the "1 +" is essential. The task brief's form ζ̈ = −ζ/r³ gives 0.514955 at L2 instead of the printed
    0.976516, so we rejected it.
- Planar Γ and units are those of H5 and Hénon 2003. M2 is at the origin and M1 is on the −ξ side.

**Index.** "we define as in I a vertical stability index a_v; an orbit will be vertically stable if |a_v| < 1" (p.317).
- Paper I (A&A 28:415, pp.416-417, read on its page images) defines it as follows.
  - Eq. (9) is the full-revolution map [[a_v, b_v], [c_v, d_v]].
  - Eq. (21) gives the half-orbit coefficients A_v, B_v, C_v, D_v, from crossing 0 to the next symmetric crossing.
  - Eq. (16): a_v = d_v = A_v D_v + B_v C_v.
  - Eq. (22): a_v − 1 = 2 B_v C_v, a_v + 1 = 2 A_v D_v.
- So a_v is half the trace of the 2×2 vertical monodromy (cos φ in the 1974 eq. 14).
- In the 1974 eq. (18), a_v ≈ ½ Δζ̇1/Δζ̇0 for large |a_v|.

**Asymptotics** (sec. 4).
- Families f and g as Γ → +∞: eqs. (8)-(14) give a_v = cos φ, φ = 2π(ω + 1)/ω, ω = ±ρ0^{−3/2} − 1, Γ = 1/ρ0 ± 2ρ0^{1/2}.
  - The paper says this is "identical with the value found for the horizontal stability index a in V".
  - Tables 2 and 3 compare it with the exact a_v.
- Family f as Γ → −∞: a_v → +1 (eqs. 15-16).
- Families a, c and g′ as Γ → −∞: a_v = Q4 (−Γ)^{3/2} (eq. 20), with Q4 = −8 sin³t1 cos t1 / (1 + 3 sin²t1) (eq. 21).
  - t1 comes from H5 Table 7. n = 1 is a and c, n = 2 is g′.
  - Table 4 gives Q4 for n = 1 to 10. Tables 5 and 6 compare it with the exact values.
  - The paper says "the second term of the expansion would be a constant". Our check: exact minus formula is about
    −3.4 to −3.9 for a and about −5.5 to −5.8 for g′.

**Printed oddities** (kept as printed in the yaml).
- Table 2 at Γ = 1 and Γ = 4: the formula values 0.1798 and 0.8627 do not follow from the exact root of eq. (10),
  which gives 0.1801 and 0.8629.
  - They are copied from H5 Table 11, where ρ0 was solved about 1e-4 off (H5 transcription, oddity 2).
- Table 4, n = 3: eq. (21) gives Q4 = −0.24070548, which rounds to −0.240705. The printed value is −0.240706, one
  unit in the last place.
- The g′₂5v ξ1 cell has a small stray dot printed before its minus sign. It is a print mark.
- Table 1 g′ sign convention differs from H5.
  - Simple-periodic part: ξ0 is H5's second ξ column and ξ1 is minus its first.
  - g′₂ rows: these start at ξ0 > 0 with η̇0 < 0. They are the images, through M2, of H5's start.
- H5 Table 5 prints the g′ ejection orbit at Γ 3.68448. This paper prints 3.684489.
  - Our extrapolation of the family to ξ0 = 0 gives 3.684488 ± 2e-6, with far crossing −0.509084.
  - So the 1974 value stands. H5's 3.68448 looks like 3.684489 truncated rather than rounded.

## 2. Witnesses (Table 1-6 transcription)

- 138 numeric cells. Names, signs, Types and infinity signs are not counted.
- **A, tesseract** (400 dpi, psm 6, digit whitelist):
  - 113 agree;
  - 19 lose only the minus;
  - 5 drop a decimal point;
  - 1 misreads (Table 6, Γ = −36: "7190" for −77.90).
- **A2, tesseract without the whitelist:** 132 agree. It differs on the same 6 cells, with every sign right.
- **B, macOS Vision** (VNRecognizeTextRequest, .accurate, no language correction; Swift CLI built with the
  MacOSX26.5 SDK, because the default SDK does not match swiftc 6.3.3):
  - 119 agree;
  - 17 single-character cells (the Γ and n columns) give no token;
  - 2 merge a neighbouring mark.
- No cell has A and B agreeing with each other against the image. Every disagreement was settled on the 400 dpi crop.

## 3. Cross-checks

**(a) Against Hénon 2003** (`data/sources/henon-2003-hill-families-tables.yaml`):
- No family is in both papers.
- One Γ coincides in the whole set: 1974 g′2 against 2003 Ha. Γ 4.271428 and ξ0 0.073664 are identical in both
  papers.
- 2003 lists this orbit with ξ0 = ξ1 = 0.073664, because N = 2 (g′ described twice). 1974 lists ξ1 = −0.587694,
  the other crossing of the single loop.
- We integrated the 1974 state. It closes with ξ1* = −0.5876939 and T/2 = 1.151034. H5 gives 1.15103.

**(b) Against Part I** (equal masses; digest 2026-10-06 and the `-tables.txt`):
- Same index. Part I defines a_v and this paper uses it "as in I". The stable range |a_v| < 1 is the same.
- Same Table 1 layout: Name, x0, x1, C, ẏ0, ẏ1, Type becomes Name, ξ0, ξ1, Γ, η̇0, η̇1, Type.
- Same naming rule: the subscript is half the number of crossings, then a number and "v".
- Same Type letters, A_v-D_v by the vanishing half-orbit coefficient.
- The Lagrange-point start value plays the same role: Part I eq. (23) at L3 gives a_v = cos 2πs_z/s_i = 0.936196,
  and this paper's eq. (7) at L2 gives 0.976516.
- Differences: Γ is Hill's reduced Jacobi constant, not C. Part I has no asymptotic tables.

**(c) a_v recomputed** (`recompute_av.out`).
- Method: integrate eq. (3) planar with eq. (6) for the 2×2 vertical matrix, using DOP853 with rtol 1e-13 and
  atol 1e-14.
- Formula: a_v = A_v D_v + B_v C_v at the half orbit (crossing N), checked against half the trace at crossing 2N.
  The two agree to 1e-7 relative or better (det M − 1 ≤ 1e-9).
- **Eq. (7) L2:** cos(2T) with T = 2π/2.071594 gives **0.976516**, as printed. Small Lyapunov orbits at
  Γ = 4.3267, 4.326 and 4.32 give 0.976520, 0.976569 and 0.976993, which tend to it.
- **Tables 2, 3, 5, 6, all 20 exact a_v reproduced to the printed digits.** We root-found the family member at the
  printed Γ, starting from H5's ξ. For example:
  - f, Γ = 1: 0.390749 (printed 0.3907);
  - g, Γ = 5: 0.538845 (0.5388);
  - a, Γ = −4: −8.214774 (−8.215);
  - a, Γ = −100: −544.2367 (−544.2);
  - g′, Γ = −36: −77.902509 (−77.90).
  - The g′ rows had to start at the crossing near M2. At the far crossing, ξ̇ at the near-collision crossing is too
    steep in ξ0 to bracket a root.
- **Table 1, all 10 vertical critical orbits reproduced:**
  - From the printed state, a_v is within 2e-5 of ±1.
  - At the printed Γ, the vanishing coefficient matches the printed Type in 10 of 10 rows.
  - For the 9 transversal crossings, the Γ* where a_v = ±1 matches the printed Γ to all 6 decimals, and so do its ξ0
    and ξ1: 27 of 27 cells.
  - g1v is a tangency: the minimum of a_v + 1 is 6e-15 at Γ 3.057455, with A_v = D_v = 1.3e-8. This fits its printed
    type A_v−D_v.
- **Formula columns** (`crosscheck.out`):
  - Table 3: 4 of 4 reproduce from eq. (14).
  - Table 2: 4 of 6 reproduce. The other two are the H5 copies noted in sec. 1.
  - Table 4: 9 of 10 reproduce. The exception is n = 3, noted in sec. 1.
  - Tables 5-6: 10 of 10 reproduce from Q4·(−Γ)^{3/2}.
  - The Table 1 entries that are also in H5 Tables 2, 4, 5 and 10 round to H5's five digits. The only exception is
    the ejection Γ (sec. 1).

## 4. Citation mining

References (p.321), checked with `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md`:

| Reference | Held? |
|---|---|
| Benest, D. 1971, A&A 13:157 (retrograde satellites, Sun-Jupiter) | NOT held (no `benest` file; not in CORPUS_INDEX). Low priority: planar stability only. |
| Benest, D. 1974, A&A "to be published" (retrograde satellites at other mass ratios; probably A&A 32:39, not verified) | NOT held. Low priority. The wanted list has only Benest 1976 (row 63). |
| Hénon, M. 1969, A&A 1:223 (H5, "V") | HELD (`henon-1969-...-aa-1-223-ads-1969AA-1-223H.pdf` + `-tables.txt`) |
| Hénon, M. 1970, A&A 9:24 | HELD (`henon-1970-...-aa-9-24-ads-1970AA-9-24H.pdf`) |
| Hénon, M. 1973, A&A 28:415 (Part I) | HELD (`henon-1973-...-aa-28-415-ads-1973AA-28-415H.pdf` + `-tables.txt`) |
| Szebehely, V. 1967, Theory of Orbits | HELD (`szebehely-1967-theory-of-orbits-...-book.pdf`) |

- No new must-have lead. The Benest pair adds nothing to vertical stability.
- The real gap for #974 and #994 remains: no published a_v exists for the Hénon 2003 families or the asymmetric
  second-species parents. A 3D claim on those must compute a_v itself. `hill.py` here is a tested starting point:
  eq. (6), the Paper I index, and positive controls against the 31 values reproduced above.

*Filed as `cyclers_pdf/papers/henon-1974-vertical-stability-periodic-orbits-restricted-problem-II-hill-case-aa-30-317-ads-1974AA-30-317H.pdf` (force-OCR copy). Table transcription: `data/sources/henon-1974-hill-vertical-stability.yaml`. Witness comparison, recomputation scripts and outputs are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<file name>`.*
