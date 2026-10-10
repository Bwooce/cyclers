# Digest: Huang 1962, "Preliminary Study of Orbits of Interest for Moon Probes" (#960 batch 43, wanted row 23)

Su-Shu Huang (Goddard Space Flight Center, NASA), The Astronomical Journal 67(5):304-310 (June 1962),
doi 10.1086/108730, ADS 1962AJ.....67..304H. Received January 20, 1962.
- DOI checked on Crossref (2026-10-11): title "Preliminary study of orbits of interest for moon probes", AJ 67, p.304,
  author Su-Shu Huang, issued 1962-06.
- Proposed filename: `cyclers_pdf/papers/huang-1962-preliminary-study-orbits-interest-moon-probes-aj-67-304-doi-10.1086-108730-ads-1962AJ-67-304H.pdf`.
  - The file to file is the `ocrmypdf --force-ocr -l eng` copy, md5 04f53dd34001f5db62dcefd9d5e57f4c, 7.3 MB.
  - Source: ADS scan, 7 pp., md5 fd14a29869dc4329dab764a3ec0da68b, 0.56 MB (200 dpi colour JPX plus a 600 dpi JBIG2
    mask per page). PDF page n = journal page 303 + n.
  - Render check at 100 dpi (`render_check.out`). Dark-pixel fractions, original against OCR copy: p.1 0.0546/0.0606,
    p.2 0.0373/0.0415, p.5 0.0549/0.0610. Same pixel size. The copy is about 11 per cent darker, with no blank page.
    A 400 dpi crop of Table I from the OCR copy is fully legible.
  - The OCR copy is 13 times the size of the scan (force-OCR rasterises the page). PROPOSAL: re-make it with `-O 3`
    added at filing, as was done for Causey 1964.
- **How it was read:** every page was read on page images.
  - Pages 304, 305, 308 and 310 (all text, Table I, eqs. 1-12, the references) were read on 400 dpi crops.
  - Pages 306, 307 and 309 are figures only (Figs. 3 and 4). They were viewed whole at about 100 dpi.
  - The OCR text layer was used only to find things. It reads "n/m = 2/3" as "n/m=2" and the ellipses in Table I as
    letters.
- Companion `huang-1962-tables.yaml`: Table I (46 values) and the initial conditions of eqs. (11) and (12).
- Scripts (all values in them are ours):
  - `integrate_huang.py/.out`: eqs. (11), (12) and the Huang & Wade 1963 rows.
  - `witness_compare.py/.out`: witnesses, plus the eq. (3) check.
  - `compare_lineage.py/.out` and `extra_checks.py/.out`: family matches and distances.
- Part II is Huang & Wade 1963 (`digest-huang-wade-1963.md`). Read the two together: Part II tabulates the two
  families that this paper finds.

## 0. Verdict

**This is the first printed member of two Earth-Moon periodic-orbit families of the 2/3 type.** In both, the probe
makes three revolutions about the Earth while the Moon makes two. Both orbits enclose the Earth and the Moon, and both
pass the far side of the Moon once per period.
- **Retrograde orbit** (eq. 11): x = −0.39215, y = 0, dx/dt = 0, dy/dt = 2.35164.
  - The paper integrates it to t = 240 ("38.20 sidereal months") and sees 19 lunar encounters, "a closed one for
    all practical purposes".
  - It finds the orbit stable to ±0.10 per cent changes in the injection speed (Fig. 3).
- **Direct orbit** (eq. 12): dy/dt = −1.61025. It leaves the resonance after the fourth lunar passage. The paper
  calls it unstable (Fig. 4).
- **Our check** (DOP853, μ = 0.01215, `integrate_huang.out`):
  - Eq. (11) returns to its start after T = 12.686626 TU (55.17 d). The velocity error is 5e-5. The far-side
    perpendicular crossing is at x = 1.16135. This is a pass at the paper's precision.
  - Eq. (12) is a five-digit rounding of a strongly unstable orbit. Its period difference p2 − p1 is +0.016.
  - Part II prints both orbits to eight digits. Both then close.
- **What it gives the project:**
  - A dated (June 1962) published source for both families.
  - The retrograde family is the project's `#997` F4 2/3-retrograde family. We matched it member by member (Part II
    digest, sec. 4).
  - The direct family is Causey 1964's "ratio 2/3 order 1". We joined them by continuation (Part II digest, sec. 4).
  - It does NOT give an Earth-grazing orbit. The perigee is 0.35 L, about 134,000 km from the Earth's centre. So it
    is not prior art for the object of `#948` R4 (sec. 4).
- **PROPOSALS only:**
  1. File the paper under the name above. Remove wanted-list row 23 together with Part II.
  2. Cite Huang 1962 as the earliest printed member of the `#997` F4 2/3-r family and of Causey's 2/3 order 1
     family. See the Part II digest for the family matches and the `#997` "direct 2/3 did not close" correction.
  3. Table I has one printed slip: m = 5, n = 2 prints 0.5401, and eq. (3) gives 0.5407. Keep it as printed in the
     yaml.

## 1. Content

**Starting ellipses** (sec. II, pp.304-305; read on the image).
- The Moon's orbit has semimajor axis 1 and period P0. The vehicle's Earth-only ellipse has semimajor axis a and
  period P.
- "P/P0 = n/m" (1), with m and n integers.
- The encounter with the Moon is at apogee: "(1+e)a = α" (2).
- "a = (1−μ)^(1/3) (n/m)^(2/3)" (3), with μ = 0.01215.
- Table I lists a for m = 1..12 and n = 1..6. It has 46 values for relatively prime (m, n) and "..." elsewhere.
- Limits:
  - α < 2a (4), so that the orbit is an ellipse;
  - α > a (5);
  - a < ½(α + γ) (6), where γ is the largest acceptable perigee.
- With α in 1.08-1.20 and γ = 0.5, only a in 0.54-0.85 is left. The paper takes n/m = 2/3 (a = 0.7600), "the case
  n/m = ½ having been studied by Message and by Newton".

**Model** (sec. III, p.308).
- Planar restricted three-body problem in the rotating frame, with the origin at the barycentre.
- Unit length is the Earth-Moon separation and total mass is 1. The Earth has mass 1 − μ and the Moon has mass μ.
  "The period of the moon around the earth is 2π."
- Eqs. (7)-(8):
  - ẍ − 2ẏ = x − (1−μ)(x−x1)/r1³ − μ(x−x2)/r2³
  - ÿ + 2ẋ = y − (1−μ)y/r1³ − μy/r2³
- Perigee start: "x = −a(1−e) − μ, y = 0" (9). So the Earth is at x = −μ and the start is on the far side of the
  Earth from the Moon.
- Start speed: "dx/dt = 0, dy/dt = ±[(1+e)/(a(1−e))]^(1/2) + a(1−e)" (10).
  - The minus sign gives direct orbits and the plus sign gives retrograde orbits.
  - The "+a(1−e)" term is the frame rotation.
- With α = 1.14 and a = 0.7600: e = 0.5 and a(1−e) = 0.38 (ours), so x = −0.39215, which is the printed value.
  - Eq. (10) gives dy/dt = 2.3668 (retrograde) and −1.6068 (direct) (ours).
  - Successive approximation moved these to the printed 2.35164 and −1.61025.
- Integration: Runge-Kutta with Δt = 0.01 on the IBM 7090, by Clarence Wade, Jr. It took about ½ hour to reach
  t = 240 in double precision.

**Results** (pp.308-310).
- Retrograde (eq. 11):
  - The final α is "about 1.16".
  - It has 19 encounters in t = 0-240, and "the orbit repeats itself every time after each encounter" (Fig. 1).
- Direct (eq. 12):
  - After three passages the orbit "deviates greatly after the fourth passage" before t = 60, and the final α is
    "about 1.07" (Fig. 2).
  - The text gives the reason: the line of apsides turns at each encounter. Also, the real lunar eccentricity would
    break the timing, so the direct orbit "should not be used".
- Encounter interval: "slightly more than two sidereal months" for retrograde orbits and "less than two sidereal
  months" for direct orbits (p.309).
  - Our periods agree: 12.69 TU > 4π for retrograde and 11.64 TU < 4π for direct.
- Stability (sec. IV, p.310):
  - Retrograde: dy/dt changed by ±0.05, ±0.10 and ±0.15 per cent, integrated to t = 60, Fig. 3. All keep the
    regularity except −0.15 per cent, whose third encounter is "on the wrong side of the moon".
  - Direct: changed by ±0.005, ±0.010 and ±0.015 per cent, Fig. 4. Each makes "only two encounters", so the orbit is
    unstable.
  - A slightly larger injection speed is "less damaging" than a slightly smaller one.
- Figs. 3-4 mark encounter times only. They carry no table, so nothing was transcribed from them.

## 2. Witnesses (Table I)

- 46 numeric cells, 26 ellipsis cells and 12 row labels.
- The image decides. The second witness is tesseract 5, psm 6, on a 400 dpi crop (`witness-tesseract/w62_T1.txt`).
- Tesseract agrees on all 46 values and all 12 labels. No normalisation was needed. It reads the ellipses as "tee",
  "see" and so on, and these cells are not counted.
- Cell file: `witness-comparison-huang1962.tsv`.
- **Eq. (3) check** (ours): 45 of 46 printed values equal (1−μ)^(1/3)(n/m)^(2/3) rounded to 4 decimals.
  - m = 5, n = 2 prints 0.5401. The formula gives 0.540676.
  - Tesseract also reads 0.5401, and the 400 dpi crop confirms it. So this is a printed slip, not a misreading.
- Every "..." is at gcd(m, n) > 1, and every value is at gcd(m, n) = 1.

## 3. Cross-checks

- **Frame and μ.** Eqs. (11) and (12), integrated with the Earth at −μ and the Moon at 1−μ and μ = 0.01215, give the
  behaviour the text describes.
  - Eq. (11) closes.
  - Eq. (12) has its first far-side pass at x = 1.07624 with ẋ ≈ 6e-5, and drifts afterwards.
  - Part II's eight-digit rows reproduce every printed period. This confirms the frame (Part II digest, sec. 3).
- **The 38.20 sidereal months** equal 240/(2π) = 38.197 (ours).
- **Part II consistency.** Part II's Table II last pair (x0 −0.39215, ẏ0 2.3516410/2.3516411) and Table I first pair
  (−1.6102480/−1.6102481) round to this paper's 2.35164 and −1.61025.

## 4. #948 R4 collision check

- `#948` R4 (`docs/notes/2026-10-05-938-fable-corpus-novel-paths-review.md` sec. 3 R4) is about **Earth-grazing**
  second-species cyclers: Bruno's e = 1 arcs and the Gómez-Ollé double-collision orbits, continued to the Earth-Moon
  mass. This needs regularisation at both primaries.
- This paper has two orbits that pass both bodies (ours; L = 384,400 km; C without the μ(1−μ) term):

| orbit | x0 | ẏ0 | C | T (TU) | closest Earth | closest Moon |
|---|---|---|---|---|---|---|
| eq. (11), retrograde | −0.39215 | 2.35164 | −0.159610 | 12.6866 | 0.3546 (136,300 km) | 0.1573 (60,460 km) |
| eq. (12), direct | −0.39215 | −1.61025 | 2.777696 | 11.6430 (not closed) | 0.3488 (134,060 km) | 0.0884 (33,980 km) |

- **Neither is Earth-grazing.** The perigees are about 134,000 km from the Earth's centre. So **this paper is not prior
  art for the R4 object**. The R4 "why unpublished" argument is untouched.
- It is prior art for the two 2/3 families that the project has already met (`#997`, Causey 1964). See the Part II
  digest sec. 4 for the orbit-level matches.
- **Not the same orbits as Newton 1959 or Schwaniger 1963.**
  - Newton's orbits are type 1/2, with a period of about 2π, and the paper says so itself (p.305).
  - Schwaniger's orbit has a period of 26 d (6.0 TU).
  - Huang's orbits are type 2/3, with periods of 11.6-12.8 TU (51-56 d).

## 5. Citation mining (p.310, read on the image)

Held status was checked with `ls cyclers_pdf/papers | grep -i <name>` and a grep of `CORPUS_INDEX.md`.

| Reference | Held? | Wanted list |
|---|---|---|
| Egorov, V. A. 1958, *The Russian Literature of Satellites* (International Physical Index Inc., New York), p.115 | no. No file, no index line. Named only in the Genova-Aldrin digest's history line | not listed. PROPOSAL: add as low priority. The title suggests a translated collection; the item is probably Egorov's Earth-Moon trajectory work (not checked). |
| Message, P. J. 1959, Astron. J. 64, 226 | no (no file, no index line) | not listed. Crossref: "Some periodic orbits in the restricted problem of three bodies and their stabilities", doi 10.1086/107925. PROPOSAL: add it. It is the n/m = ½ Earth-Moon study that this paper names, and Part II uses it for stability. |
| Moulton, F. R. 1914, *An Introduction to Celestial Mechanics* (Macmillan), Chap. 8 | no | not needed (textbook) |
| Newton, R. R. 1959, Smithsonian Contrib. Astrophys. 3, 69 | **yes**: `newton-1959-periodic-orbits-planetoid-passing-close-two-gravitating-masses-...-ads-1959SCoA-3-69N.pdf` (index line 259) | removed in batch 30 |

- Wanted-list numbering: the 2026-10-05 wanted list has these papers as row 23. The 2026-10-11 DOI-lookup note numbers
  the same entry 24.

*Filed as `cyclers_pdf/papers/huang-1962-preliminary-study-orbits-interest-moon-probes-aj-67-304-doi-10.1086-108730-ads-1962AJ-67-304H.pdf`. Check scripts, outputs, witness tables and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/huang-1962-tables.yaml`.*
