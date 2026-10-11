# Digest: Kotoulas and Voyatzis 2004, "Comparative study of the 2:3 and 3:4 resonant motion with Neptune: an application of symplectic mappings and low frequency analysis" (CMDA 88:343-363) (#960, #975 R14)

T. Kotoulas and G. Voyatzis (Department of Physics, University of Thessaloniki), Celestial Mechanics and Dynamical Astronomy
88:343-363 (2004), DOI 10.1023/B:CELE.0000023391.85690.31. Received 22 May 2003, revised 5 September 2003, accepted 30 October 2003.
Wanted-list row 48 (Kotoulas-Voyatzis is the first item; Bray and Goudas 1967 and Voyatzis, Kotoulas and Hadjidemetriou 2005 stay wanted).

- **Filed as** `cyclers_pdf/papers/kotoulas-voyatzis-2004-comparative-study-2-3-and-3-4-resonant-motion-neptune-symplectic-mappings-cmda-88-343-doi-10.1023-B-CELE.0000023391.85690.31.pdf`.
  Table transcription: `data/sources/kotoulas-voyatzis-2004-tables.yaml`.
  - The file is the publisher PDF as received, md5 0743699683d9296e5f577efa9f314560 (709,285 bytes), 21 pages, A4, text layer
    present. No OCR was needed. PDF page n = journal page 342 + n.
- **How it was read:** the whole text layer (all 21 pages), then the page renders at 70 dpi for every figure page, and Table I
  (the only table) at 300 dpi, every cell by eye. Three witnesses on Table I: the text layer, the 300 dpi image read by eye,
  and tesseract on the crop. All 58 printed values agree. A script also compared the 29 rows of the yaml with the text layer: no
  difference.
- **Read first for #975:** `voyatzis-kotoulas-2005-...-r14-verdict.md` (private papers directory) and
  `docs/notes/2026-10-07-digest-voyatzis-kotoulas-2005-exterior-resonances-neptune.md`.

## 0. Verdict

**A Kuiper-belt stability study of the 2:3 and 3:4 exterior resonances with Neptune, done with two 6-dimensional symplectic
mappings and a low-frequency spectral indicator. It has no periodic-orbit table, no interior family and no closed family,
and it does not touch the R14 target.**
- The models are mappings M23 and M34 of the three-dimensional elliptic restricted problem (Sun, Neptune, a massless body)
  at mu = 5.178e-5, built from an averaged Hamiltonian with correction terms (Sec. 2, Eqs. 1-2, Appendix A).
- It then classifies about 700 trajectories per resonance on an (a, e) grid at several inclinations as regular or chaotic and
  draws gray-scale "LFP maps" (Figs. 7-11). The conclusion is that the 3:4 phase space is more chaotic than the 2:3 one, which fits
  the lack of 3:4 objects.
- The only numbers are Table I (29 mapping coefficients for each resonance), a handful of initial conditions of mapping orbits,
  grid and run parameters, and a few stability limits. All are in the tables file.

## 1. Content

- **Sec. 1:** introduction. It cites Hadjifotinou and Hadjidemetriou 2002 and Kotoulas and Hadjidemetriou 2002 for the periodic
  orbits of the planar circular problem at 2:3 and 3:4 and their stability, with bifurcation points of the same type for both
  (families then exist in the planar elliptic problem and the 3D circular one). Those orbit results are cited, not reproduced.
- **Sec. 2, mappings:** the generating function is F = sigma_n S_(n+1) + sigma_z,n Sz_(n+1) + nu_n N_(n+1) + T H-bar(...), with T the
  synodic period of the resonant periodic orbit (6 pi for 2:3 and 8 pi for 3:4). H-bar = H0 + mu H1 + mu e H2 + mu z^2 H3 +
  mu e z^2 H4 + Hcor, z = sin(i/2). Corrections make the periodic orbits of the original problem coincide with the mapping fixed
  points in position and stability, and move the Kozai resonance to e about 0.25 (2:3) and 0.20 (3:4).
  - Setting e = 0 and i = 0 gives a 2-dimensional map that is the Poincare section of the planar circular problem. Resonant
    periodic orbits are the fixed points at sigma = 0 and sigma = pi. The stable one (sigma = pi) lies at S = S0 and constant
    eccentricity e0.
  - The stable fixed points form one-parameter families F23 and F34 in e, at a_R about 39.405 AU (2:3) and 36.415 AU (3:4).
    In the planar circular problem both families extend to high e except for collision orbits. In the mappings F23 stays stable to
    e = 0.63 and F34 to e = 0.46. These are the only family-level numbers. **No initial conditions (x0, ydot0), periods or Jacobi
    constants are printed for F23 or F34.**
- **Sec. 3, orbit types:** five orbit types T0 to T4 for 2:3 at i = 5 deg (regular, weakly chaotic, strongly chaotic, sticky,
  stable chaos), with the (a, e) of each (tables file). Escape condition |a - a_R| > 1.5 AU. Grids of 22 x 31 (2:3) and 23 x 30 (3:4)
  points, integrated to 4 Gy (Fig. 4).
- **Sec. 4, low-frequency analysis:** the spectrum of the action S_n of 131,072 iterations; b = b8 is the log of the power in the
  eight lowest frequencies over the total. Eq. 4 gives the perturbative form of the spectrum with small denominators.
- **Secs. 5-6, LFP maps:** for each of the planar circular, planar elliptic, 3D circular and 3D elliptic cases. The 2:3 region is
  stable for e up to 0.63 at i below 3 deg, with Kozai stabilisation to i about 20 deg and Pluto at the edge (Fig. 8c). The 3:4 region
  is stable for e up to 0.45, is destroyed above i of 15 deg, and is more affected by Neptune's eccentricity and by inclination.
- **Sec. 7, conclusions:** secondary resonances destroy low-e stability, the high-order resonances f_sigma shrink the stable area, the Kozai
  resonance protects, and the 3:4 dynamics is "significantly more chaotic" than the 2:3.

## 2. Figures (described)

- **Fig. 1:** four panel pairs, the mapping sections above and the planar circular Poincare sections below, at e0 = 0.1 and 0.4 for each
  resonance. Closed invariant curves around the stable fixed point at sigma = pi, a separatrix from the unstable point at sigma = 0.
- **Figs. 2-3:** time histories of a, e, i and sigma for T0 to T4 (to 4 Gy).
- **Fig. 4:** the (a, e) plane for 2:3 and 3:4 with the families F23 and F34, the separatrix, high-order and secondary resonance curves, the
  Kozai borders, and marks for escape time (squares under 100 My, triangles to 1 Gy, stars to 4 Gy).
- **Figs. 5-6:** power spectra. **Figs. 7-10:** LFP maps on (a, e). **Fig. 11:** LFP maps on (e, i).

## 3. Relevance to this project

**What the 2004 paper adds for #975 (route R14, the Varin closed-family cascade i_5-i_8 at the satellite mass ratios): nothing
that collides with it, and no numbers it can use.**
- **Resonant families at the Neptune mass ratio:** two only, 2:3 and 3:4 (exterior, a > 1), each a family of stable symmetric fixed
  points in e at mu = 5.178e-5. They are named and bounded (stable to e = 0.63 and 0.46 in the mappings) but not tabulated. There
  are no periodic-orbit initial conditions, no periods beyond T = 6 pi and 8 pi, no Jacobi constants and no stability indices.
  Hadjifotinou and Hadjidemetriou 2002 and Kotoulas and Hadjidemetriou 2002 are the sources of such numbers and are not held.
- **Stability results:** only for mapping orbits and only as regular or chaotic maps over (a, e, i). No Floquet or stability index of a
  periodic orbit appears. The mapping coefficients of Table I are model constants, not orbit data.
- **Interior i_5-i_8 or any planar periodic orbit with numbers:** no. Both resonances are exterior (a_R = 39.4 and 36.4 AU). The word
  "interior" does not occur as an orbit class. No orbit is a periodic orbit with printed state, and no value of a Jacobi constant is printed.
- **Collision with the R14 target:** none, agreeing with the 2005 verdict. The target is interior closed families branching from
  family i at mu = 2.5e-5 to 2.4e-4. The 2004 paper has no closed family and no interior family, and it is at one mu, 5.178e-5.
  The 2005 paper's only closed family (exterior, p = -7, figure only) is not in the 2004 paper either.
- **Link to the 2005 paper (inference):** the 2:3 and 3:4 families of 2004 are the first-order exterior resonant families that the 2005
  paper follows far past the mapping range (I and II of 2/3 and 3/4, in its curves 1 to 3 of Fig. 1, which the 2005 digest places on
  open curves, not on the closed p = -7 one). The 2004 mapping stops at e = 0.63 and 0.46. The families of 2005 go on to Neptune
  and Sun collisions. This is a reading of the two papers, not a statement in either.
- **Possible small use, not needed for #975:** the printed a_R values and the 0.63 and 0.46 stability limits give a loose shape check
  for a planar circular solver at mu = 5.178e-5 on the 2:3 and 3:4 stable families, in the same spirit as the 2005 positive control.
  The limits are mapping limits, which the paper says are "absolute limitations in the validity of the mappings", not limits of the
  true families. a_R is about 39.405 and 36.415 AU against Kepler 39.403 and 36.427 AU (our check). Not a digit control.
- **Catalogue implication:** none. No row is a cycler, and the paper is exterior Kuiper-belt dynamics.

## 4. Citation mining (against CORPUS_INDEX.md and the wanted list, 2026-10-11)

- **Held neighbour:** Voyatzis and Kotoulas 2005 (PSS 53:1189), with its R14 verdict.
- **Added as a wanted row (row 61):** Hadjifotinou and Hadjidemetriou 2002, CMDA 84:135-157 ("A symplectic mapping model for the study of
  2:3 resonant trans-Neptunian motion") and Kotoulas and Hadjidemetriou 2002, Earth, Moon, and Planets 91:63-93 ("Resonant periodic
  orbits of trans-Neptunian objects"). These are the two papers this one cites for the planar circular periodic orbits and their
  stability at 2:3 and 3:4, so they are the likely source of Neptune-mu periodic-orbit numbers for #975. Low to medium priority.
- **Not added (Kuiper-belt stability background, no periodic-orbit tables expected):** Malhotra 1996 (AJ 111:504), Morbidelli 1997 and
  1999, Morbidelli, Thomas and Moons 1995, Nesvorny and Roig 2000 and 2001 (Icarus 148 and 150), Gallardo and Ferraz-Mello 1998,
  Melita and Brunini 2000, Yu and Tremaine 1999, Duncan, Levison and Budd 1995, Michtchenko and Ferraz-Mello 1995 and 2001.
- **Not added (method):** Hadjidemetriou 1998 (symplectic maps, book chapter), Voyatzis and Ichtiaroglou 1992 (spectral analysis),
  Laskar 1993, Froeschle et al. 1997, Sidlichovsky 1991, Press et al. 1992, Murray and Dermott 1999.
- **Cited as "submitted":** Kotoulas and Hadjidemetriou, "Dynamical evolution of trans-Neptunian objects near the 3/4 resonance", Earth,
  Moon, and Planets (submitted 2003). Not searched; the journal reference is not printed beyond "submitted".

*Filed as `cyclers_pdf/papers/kotoulas-voyatzis-2004-comparative-study-2-3-and-3-4-resonant-motion-neptune-symplectic-mappings-cmda-88-343-doi-10.1023-B-CELE.0000023391.85690.31.pdf`. Table transcription: `data/sources/kotoulas-voyatzis-2004-tables.yaml`.*
