# Digest: Wagner 1968, "The Restricted Problem of Three Bodies (III)", NASA Technical Report No. 4 (NASA CR-95832, NTRS 19680020622) (#960, owner-supplied)

Curtis A. Wagner, Department of Physics, University of Illinois, Urbana. NASA Technical Report No. 4, NASA Grant NsG 280-63, dated 1 June 1968 (cover page). NTRS record
19680020622, titled "The restricted problem of three bodies /3/ Final report", also numbered TR-4 and NASA-CR-95832, grant NSG-280-63, publication date 1 June 1968 (NTRS API,
checked 2026-10-11). Sponsor: NASA, under grant NsG 280-63 (the grant ended 31 August 1967 after two extensions). The last two pages are the final report of that grant.

**Series.** This is part III of the Bartlett-Wagner series on the restricted problem. Part I is Bartlett 1964 (also NASA Technical Report No. 1); part II is Bartlett and Wagner 1965
(also NASA Technical Report No. 2). The report's own reference list cites both with venues (sec. 4). It was found as a report; no journal form is cited for part III.

- **Filed as** `cyclers_pdf/papers/wagner-1968-restricted-problem-three-bodies-III-symmetric-solutions-asymptotic-limiting-orbits-ntrs-19680020622.pdf`.
  Numbers: `data/sources/wagner-1968-numbers.yaml`.
  - Filed copy: `ocrmypdf --redo-ocr -l eng -O1`, md5 03d6971426be365891d095cbb7d8b004, 4,321,069 bytes, 24 pages, of the received scan (md5 67aa8032872efc31697b079587e91cc9, 4,384,626
    bytes), whose text layer is letter-spaced and poor. A plain run refused ("page 2 already has text").
- **How it was read:** the OCR text and every text page as an image at 90 dpi (the page images are sharp), the references and final-report pages, and Figs. 2, 5 and 7 as images; the other
  five figure pages were seen at thumbnail level only. Twenty-five numeric values were compared between the OCR text and the images; all agree.

## 0. Verdict

**A short numerical study (9 text pages and 8 figures) of how Strömgren's periodic-asymptotic limiting orbits at L4 and L5 change with the mass ratio, and of one new hybrid class.**
It has no table of orbits (the relevant table is in Bartlett 1964, part I, which is not held) and no orbit with printed initial conditions. It adds mass-ratio thresholds at which
classes of symmetric periodic orbits close up, merge or change termination. It is classification background, not data for the project's frames.

## 1. Content

- **Setting (pp.1-4):** Thiele coordinates (E, F); the mass-ratio parameter gamma = (m1 - m2) / (m1 + m2); the energy parameter T = K - gamma^2 where K is the Jacobi constant in the paper's scaling. At L4,
  K = 11 + gamma^2, so T = 11 for every gamma, and all periodic-asymptotic limiting orbits lie on the T = 11 line. L4 is at E = -pi/2, F = sinh^-1(sqrt 3) = 1.316958. Spiral motion
  at L4 exists for |gamma| < gamma_0 = sqrt(23/27) = 0.922958 (the Routh critical mass ratio, mu = 0.0385209), closed ellipses beyond. The slope of the outgoing spiral at E = -pi/2 is given
  in closed form (p.3), which is how the limiting orbits are launched from L4.
- **Scan of the limiting orbits (pp.4-7, Figs. 1-3):** the 14 limiting orbits of Strömgren that are symmetric about the xi axis, tabulated in Bartlett 1964 at gamma = 0, were traced as functions of
  gamma by a continuation on F_i. Findings (values in the numbers file):
  - Pairs of originally similar orbits coalesce smoothly: I and II at gamma about -0.059 (closing of classes (k) and (mu), which begin and end only on I and II); III and IV at gamma about -0.24
    (closing of class (beta), which becomes an infinite nested set of small closed classes shrinking to points, the outermost vanishing near gamma about -0.83); XI and XIV at gamma about -0.582.
  - Orbit VIII reaches a minimum gamma and continues as a conjugate orbit VIII*; they merge at gamma about -0.70, after which class (g) persists with different termination and, through a second branch, (f) symmetry.
  - Orbits V and VI persist over almost the whole range, down to gamma = -0.85, with no partner. The calculations stop near |gamma| = 0.85 for accuracy.
- **eta-axis symmetric orbits (pp.7-10, Figs. 4-8):** Strömgren's second eta-axis limiting orbit (eta_0 = 1.7008 at gamma = 0, encircling L1 and spiralling into L4) has a compact class traced by the
  "coiling" property (Fig. 5, extreme orbits at T = 11.00 and T = 9.486659). A **hybrid class** is found for the combination of Strömgren's first and second eta-axis orbits (first has eta_0 about 1.72025):
  Figs. 6-7 follow nine members, T from 10.7753 down to 10.57482, up to a maximum of about 12.5, and back, with cusps becoming loops that encircle L5. The class begins and ends on the same hybrid
  limiting orbit. The author suggests that the number of possible multi-membered hybrid limiting orbits is almost limitless and that they could generate "by far the most abundant and important" classes.
  He states, but does not check, that the simple and hybrid classes exist over wide ranges of gamma.
- **Final report (PDF pp.23-24):** a one-page summary of the grant's four technical reports (sec. 4).

## 2. Figures (described)

- **Fig. 1 (PDF p.15):** F_i against gamma for 11 (the text says 10 of 14) limiting orbits, a tightly nested band widening as gamma tends to plus or minus gamma_0.
- **Fig. 2 (PDF p.16):** final position E_f against gamma for orbits III, IV, VIII and VIII*; III and IV join at the bottom of their curve (E about 0.75, gamma about -0.24, as the text says) and VIII and VIII* join at the foot
  of theirs (E about -2.27, gamma about -0.70). **Fig. 3:** F_f against gamma for V, VI, XI, XIV.
- **Fig. 4:** Strömgren's second limiting orbit and two extreme orbits of its class. **Fig. 5 (PDF p.19):** the T against F profile of that class, a single curve ending in a small coil.
- **Figs. 6-7 (PDF pp.20-21):** the nine orbits of the hybrid class in the (E, F) plane around L4 and L5; Fig. 7 (read as an image) shows the lower orbits 6 to 9 growing loops around L5. **Fig. 8:** T against F
  profiles of the hybrid class.
- A text inconsistency: the Fig. 1 caption says 11 orbits and the text says 10 of the 14. The figure was not counted here.

## 3. Relevance to this project

- **Symmetric-solution classes against the held Strömgren 1933 and Hénon 1965 papers.**
  - The classes named here, (k), (mu), (beta), (g), (f), (alpha-delta), are the Copenhagen classes. The held Strömgren 1933 translation (digest 2026-10-06) is the source of the class lettering (a to n,
    with k, g, f, m named there) and the held Hénon 1965 I and II (digests 2026-10-06) merge Strömgren and Bartlett 1964 and 1965 into 21 classes at mu = 1/2. Hénon 1965a notes that all the limiting periodic-asymptotic
    orbits have C = 11/4.
  - **Consistency check (ours):** Wagner's K = 11 at L4 for gamma = 0 is 4 C with C = 11/4, because C at L4 is 3 - mu (1 - mu) = (11 + gamma^2)/4. So Wagner's T = K - gamma^2 = 11 line is the C = 11/4 line of
    Hénon at gamma = 0. The same orbits, in another normalisation.
  - What this report adds to Hénon 1965 is the dependence on mass ratio: Hénon fixes mu = 1/2 (gamma = 0) and Bartlett and Wagner 1965 extend to unequal masses (cited, not held). Wagner 1968 gives the
    thresholds at which classes close (gamma about -0.059, -0.24, -0.83, -0.582, -0.70).
- **Which mass ratio:** gamma from -0.85 to +0.85, i.e. mu = (1 - gamma) / 2 from 0.075 to 0.925, with gamma = 0 equal masses and the spiral regime ending at the Routh value mu = 0.0385. By symmetry gamma and
  -gamma swap the two masses, so the classes are sampled for mu from 0.075 to 0.5. Two consequences:
  - Earth-Moon (mu = 0.0122), Sun-Jupiter and the satellite systems of the outer planets have mu below the Routh value, so they lie in the closed-ellipse regime, where this paper's spiral limiting orbits do not exist.
  - A system with mass ratio above 0.0385 is in the range. Pluto-Charon (m_Charon / (m_Pluto + m_Charon) about 0.109, from standard masses, ours) is, so these classes are the right comparison for the #998 Pluto-Charon lane;
    the paper does not mention it.
- **Resonant or repeated-encounter orbits?** No. The orbits are asymptotic to the triangular points and the classes lie around L1, L4 and L5 and the masses; none is a mean-motion resonant orbit and none
  is described as repeating encounters with a body. Class (g) is said to begin at one of the masses (a small orbit around the mass, not a close-encounter cycler). No cycler content.
- **Collision with the project's rows:** none. No orbit is tabulated with initial conditions, so nothing can be checked against a row or a family beyond the thresholds above. Part I (Bartlett 1964) has the
  table (14 limiting orbits, initial and final conditions, class lists) and would be the useful source for a numerical check.
- **Catalogue implication:** none.

## 4. Citation mining (against CORPUS_INDEX.md and the wanted list, 2026-10-11)

| Ref. | Item | Status |
|---|---|---|
| 1 | Strömgren, "Fortgesetzte Untersuchungen über asymptotische Bahnen im Probleme Restreint", Mat. Fys. Medd. Dan. Vid. Selsk. 10(7) (1930) | neither; **proposed row** (low to medium): the original of the periodic-asymptotic limiting orbits; the held Strömgren 1933 (Bull. Astron. 9:87) is the later summary |
| 2 | Bartlett, "The Restricted Problem of Three Bodies", Mat. Fys. Skr. Dan. Vid. Selsk. 2(7) (1964); also NASA Technical Report No. 1, NsG 280-63 (1964) | **wanted (row 19, part I)**; stays wanted; the NASA technical-report form may be on NTRS (not searched) |
| 3 | Bartlett and Wagner, "The Restricted Problem of Three Bodies (II)", Mat. Fys. Skr. Dan. Vid. Selsk. 3(1) (1965); also NASA Technical Report No. 2 | neither; **proposed row** (medium): the unequal-mass extension that this report builds on. Row 19 names only part I |
| final report | Bartlett, motion under a periodic cubic force, Mat. Fys. Medd. Dan. Vid. Selsk. 36(11) (1968); NASA Technical Report No. 3 (with Wagner) | neither; low; not celestial |

- **Wanted-list handling:** row 19 is Bartlett 1964, part I. Part III is recorded as received here; part I stays wanted; part II is proposed as an added item of the same row. The wanted list is not edited until
  the public-repo window opens.
- Held neighbours (CORPUS_INDEX): Strömgren 1933 (French and English translation), Hénon 1965a and 1965b (French and translations).

*Filed as `cyclers_pdf/papers/wagner-1968-restricted-problem-three-bodies-III-symmetric-solutions-asymptotic-limiting-orbits-ntrs-19680020622.pdf`. Numbers: `data/sources/wagner-1968-numbers.yaml`.*
