# Digest: Ross (c. 2021), "A cycler-quartet between Venus and Sol/Terra L1" (#960, #942)

D. R. Ross (independent, USA), "A cycler-quartet between Venus and Sol/Terra L1", unrefereed manuscript,
academia.edu item 45489864. Undated; the text mentions "early March 2021" (p.4). No DOI.
- Filed as `cyclers_pdf/papers/ross-2021-cycler-quartet-venus-sol-terra-l1-academia-45489864-manuscript.pdf`.
  11 pages, born-digital text layer, md5 16e1f79e878c0593d23a0d54a8680ec5. Supplied by the owner.
- I read all pages from the text layer, and checked the result block (p.6) on the page image.

## 0. Verdict

**No collision with `#942` ev-C.** This closes the "possible class collision" in
`2026-10-06-942-evC-prior-art-search.md`; Ross is the nearest unrefereed class relative.
- **"2L4" is a POWERED cycler.** Abstract: "It is a powered cycler that needs adjustment at aphelion".
- **It never encounters Earth.** Its aphelion is 0.979 AU, inside Earth's perihelion of 0.98329 AU
  (p.5). It targets the Sun-Earth L1 point (STL1) instead.
- **Its Venus flyby is infeasible as stated.** The required turn is 108.9 deg at v_inf 3,864 m/s, with
  closest approach 4,983 km from Venus's centre, inside the planet (radius 6,050 km). The author
  says: "A course-correction cannot be done naively" (sec. 5, p.7).
- ev-C, by contrast, is ballistic: Venus hosts the returns, there is a real Earth encounter (with almost
  no Earth turn), and v_inf is V 13.17 / E 9.07 km/s.
- The shared 2-synodic period (1,167.84 d) and the matching 123-d Earth-to-Venus leg are a same-period
  coincidence (lead's reading, confirmed). The v_inf, the Earth encounter and the ballistic status all
  differ.

## 1. Content (READ)

- **Framing (secs. 1-2):** adapts McConaghy, Longuski & Byrnes's nPr method to Venus as the home planet,
  with Earth's Hill sphere as the "outer planet". Assumptions are circular, coplanar and conic, with
  synodic period S = 13/5 Venus years. Notes the 8-yr Earth / 13-yr Venus near-commensurability and the
  Hohmann "5(1.0)10" Venus-Earth cycler.
- **Method check (sec. 3):** Gooding's Lambert solver, in C# with single-precision vectors. It reproduces
  McConaghy's Aldrin 1L1: aphelion 2.229 AU, v_inf at Earth 6,537 m/s, at "Mars" 9,730 m/s, 146.8 d,
  turn 83.7 deg, against McConaghy's 2.23, 6,540, 9,750, 146 and 84.
- **2L4 (sec. 4, p.5-6; p.6 block on the page image):**
  - left branch of 2P4, semimajor axis 0.84 (Venus AU);
  - aphelion 0.979 AU;
  - v_inf at Venus 3,864 m/s; v_inf at STL1 1,950 m/s;
  - "Shortest transfer time 161 days" (the abstract says 159 d);
  - required turn at Venus 108.9 deg.
  - A time-of-flight "fudge" shortens the transfer angle to 1.25 rad, against 1.256637, to absorb the
    8:13 mismatch.
  - 3S6 (right branch, v_inf 2,441 / 809 m/s) reaches only 0.87 AU. 1L2 has v_inf 16,530 / 11,426 m/s
    and perihelion 0.4278 AU.
- **Correction (sec. 5):** r_p = (1/sin delta - 1) GM/v_inf^2 with sin delta = 0.813646, so r_p =
  4,983 km. (The "delta" in that formula is half the 108.9-deg turn.) The author proposes aphelion
  maneuvers on the model of Byrnes, Longuski & Aldrin 1993 (about 230 m/s).
- **Quartet (secs. 6-8):** the Venus-to-Earth departure (mid-1999 example) and an Earth(STL1)-to-Venus
  twin of 123 d (282-d period minus 159). Downtime about 1,008 d. Two synods = 1,167.84 d. Four vehicles,
  against ten for the Hohmann cycler.

**My checks:**
- r_p = (1/0.813646 - 1) x 324,859 / 3.864^2 = 0.229 x 21,757 = 4,983 km. Matches.
- 159 + 1,008 = 1,167 d, about 2 x 583.9 d.

## 2. Gate relevance

- **`#942` R1(b), ev-C: no collision** (sec. 0). It is cited as an unrefereed Venus-hosted 2-synodic
  concept that is powered, has no Earth encounter, and requires a sub-surface Venus flyby as stated.
- **ev-A, ev-B:** no relation (those are two-working-body).

## 3. Positive controls

- The Aldrin 1L1 reproduction (sec. 3) restates McConaghy et al. 2002 values, which are held. Not
  independent.

## 4. Citation mining

- Hollister & Menning, AIAA 69-931 (1969), "Interplanetary Orbits for Multiple Swingby Missions",
  doi 10.2514/6.1969-931. This is the conference version of the held H&M 1970 JSR paper ("presented as
  Paper 69-931", p.1193), so the content is held.
- Hollister, W. M. (1969), "Castles in Space", Astronautica Acta 14(2):311-316. Not held; no DOI found
  by Crossref. Added to the wanted list (Tier D, history).
- Morrison, O. (2018), "Use of Manifolds in the Insertion of Ballistic Cycler Trajectories", MS thesis,
  Cal Poly, doi 10.15368/theses.2018.80. Free. Added to the wanted list (Tier D). Received in batch 29: the thesis uses a southern
  halo about Sun-Earth L2, not L1 (see `2026-10-06-digest-morrison-2018-manifolds-insertion-ballistic-cycler-ms-thesis.md`).
- David, H. (2007), "The Case for Venus" (web page). Not acquired.
