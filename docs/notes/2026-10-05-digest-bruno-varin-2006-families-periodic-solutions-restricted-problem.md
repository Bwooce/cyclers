# Digest: Bruno & Varin 2006, "On families of periodic solutions of the restricted three-body problem" (#960)

A. D. Bruno and V. P. Varin (Keldysh Institute of Applied Mathematics, Moscow), "On families of periodic
solutions of the restricted three-body problem", Celestial Mechanics and Dynamical Astronomy 95:27-54
(2006), doi 10.1007/s10569-006-9021-1.
- Received 17 Nov 2005; published online 15 Aug 2006.
- Filed as `cyclers_pdf/papers/bruno-varin-2006-families-periodic-solutions-restricted-three-body-problem-cmda-95-27-doi-10.1007-s10569-006-9021-1.pdf`.
  28 pages, text layer, md5 7a73f9b6e3b0545918078edbbe328560. Journal page = PDF page + 26.
- I read the text layer in full. Tables 1-2 (p.40) and Table 4 (p.46) were read as page images.

## 0. Verdict

This is a programme-opening paper.
- It sets out methods for computing the symmetric periodic families of the planar circular restricted
  problem for all mu in [0, 1/2].
- It announces a cycle covering nine families, a, b, c, f, g, h, i, l and m, at
  mu = 0, mu_J, 0.1, 0.2, 0.3, 0.4 and 0.5 (Eq. (1.4), p.29).
- It delivers data for ONE family, h (retrograde about the heavier primary P1), at mu = 0 and at
  mu = 0.00095 (Sun-Jupiter). Larger mu is only summarised, and the detail is in two KIAM preprints
  (Bruno & Varin 2005b, c).
- So it does NOT itself tabulate the symmetric families for all mu in [0, 1/2]. The `#938` note's sec. 5
  item 11 description ("symmetric families for all mu in [0, 1/2]") overstates what this paper contains.
- The paper does not mention Leiva & Briozzo; the atlas is later work.

## 1. Content (READ)

- Model and coordinates (pp.27-31):
  - Planar circular restricted problem as a Hamiltonian, Eqs. (1.1)-(1.2), in synodic coordinates centred
    on P1 (heavier) with P2 at (1, 0).
  - Symmetric solutions cross the symmetry plane x2 = y1 = 0 twice.
  - Four coordinate systems are used for the characteristics: two global and two local "astronomical"
    systems (a, e about P1 or P2) (Sect. 2).
  - Plane and vertical traces Tr and Trv are the stability indices, computed after Henon & Guyot 1970
    Formula (18).
- Under mu -> 1 - mu the families map as a, b, c, f, g, h, i, l, m -> c, b, a, h, i, f, g, l, m (p.30),
  the same map as Henon & Guyot Eq. (12).
- Generating families (Sect. 3): the limit mu -> 0. Families of first species (circular Ir), second kind
  (elliptic EN) and second species (collision arcs A0, A1, ...), from Bruno 1990/1994 and Henon 1997.
- Computation (Sects. 4-5):
  - Fixed-step RK5 with shooting and Newton continuation.
  - Thiele-Burrau and Levi-Civita regularisation near collisions.
  - Fourier/FFT representation for difficult parts, and optional arbitrary precision.
  - Jacobi drift over a half period below 1e-10 (p.33).
- **Family h at mu = 0 (Sect. 6):** 16 critical orbits (collision, Tr = +-2, or Trv = +-2) in Tables 1-3
  (pp.40-41). Table 1 columns: k, x1(0), y2(0), v1(T/2), v2(T/2), normalised period T~ = T/2pi.
  Examples (READ p.40, image):
  - k = 2: x1(0) = -2.30499, y2(0) = -0.43384, v1 = -0.56616, v2 = 0, T~ = 1.50000; Table 2 C = 2.679465.
  - k = 7: x1(0) = -2.57909, y2(0) = 0.45608, T~ = 2.19007; C = -1.785103.
- **Family h at mu = 0.00095 (Sect. 7, p.46 ff.):** 35 critical orbits in Tables 4-6.
  - The section heading says mu = 0.00095388, but "all computations for the family h were done for
    mu = mu_J = 0.00095 for technical reasons" (p.46).
  - Table 4 columns: k, x1(0), y2(0), x1(T/2), y2(T/2), T~.
  - Examples (READ p.46, image):
    - k = 1: -0.93273, 1.03212, 0.92705, -1.04392, T~ 0.47360.
    - k = 12: -3.17449, 0.00113, 0.00000, infinity, T~ 2.00034.
    - k = 35: -5.85779, -0.17069, 1.00000, infinity, T~ 5.50206.
  - Table 7 maps the mu = 0 critical orbits to the mu_J ones.
- **Evolution to mu = 1/2 (Sect. 8, pp.51-52):** a summary only.
  - h has no self-bifurcations as mu grows.
  - It intersects families a and e at the k = 1 and 2 orbits.
  - For mu > 0.3, complete linear stability coincides with plane stability.
- **Possible erratum, offered with respect (p.29):** "For mu = mu_M = 0.1215585, corresponding to the Earth
  (P1) - Moon (P2) case, by Broucke (1968)". The Earth-Moon value is 0.01215585; one zero is missing.

## 2. Relation to the held corpus

- **Bruno 1981 ("On periodic flybys of the Moon", Celest. Mech. 24:255), held:** cited here. That paper
  gives the mu = 0 Earth-Moon flyby arcs. This 2006 paper is the general-family programme that builds on
  Bruno 1990/1994 (book, not held; doi 10.1515/9783110901733).
- **Henon & Guyot 1970, filed in this batch:** the stability formula and the critical-orbit sub-families
  for h. Section 8 cites H&G pp.369-371 for the Tr = +-2 sub-families. The H&G family-h tables are the
  mu > 0 critical-orbit control for h.
- **Hitzl & Henon 1977a, b; Henon 1997, 2001; Perko 1981:** cited as the generating-family and
  second-species sources.
- **The Leiva-Briozzo atlas cross-check:** this paper gives family h at mu = 0 and 0.00095 only. For the
  atlas cross-check at other mu, the KIAM preprints 2005b (mu = 0.1, 0.2) and 2005c (mu = 0.3, 0.4, 0.5) are
  the data sources. Neither is held.

## 3. Positive controls

- Tables 1-3 (mu = 0) and Tables 4-6 (mu = 0.00095): 16 and 35 critical orbits of family h, with x1(0),
  y2(0), period, C and traces, to 5-6 decimals.
- Coordinates are P1-centred synodic (P2 at (1, 0)). The H&G tables are barycentric. Convert before
  comparing.
- Use mu = 0.00095 exactly for Tables 4-6, not 0.00095388.

## 4. Citation mining (references pp.53-54)

Held:
- Breakwell & Perko 1974.
- Broucke 1968.
- Bruno 1981.
- Henon & Guyot 1970 (filed by `#960`).
- Henon 1968.
- Henon 1997 (Generating Families, LNP m52).
- Hitzl & Henon 1977a, b.
- Perko 1981.
- Szebehely 1967.
- Henon 2001 (LNP m65), held as `henon-2001-generating-families-...-doi-10.1007-3-540-44712-1.pdf`.
- Bruno 1978a, b (Celest. Mech. 18), held under the transliteration "Brjuno" (`brjuno-1978-...` and
  `brjuno-1978b-...`). Correction 2026-10-05: the first version of this digest listed both as not held.

Not held, in priority order:
1. Bruno, A. D. & Varin, V. P. (2005b), "Family h of periodic solutions of the restricted problem for small
   mu", KIAM Preprint 67; and (2005c) "... for big mu", KIAM Preprint 64. Russian. The family-h data for
   mu = 0.1-0.5.
2. Bruno, A. D. (1994), "The Restricted 3-Body Problem: Plane Periodic Orbits", de Gruyter, doi
   10.1515/9783110901733 (CONFIRMED in batch 1).
3. Bartlett, J. H. (1964), "The restricted problem of three bodies (1)", Kong. Dan. Vidensk. Selsk.
   Mat.-Fys. Skr. The mu = 1/2 families.
4. Bray & Goudas 1967 (3D doubly-symmetric orbits); Voyatzis & Kotoulas 2005 and Voyatzis, Kotoulas &
   Hadjidemetriou 2005 (Neptune exterior resonances); Kotoulas & Voyatzis 2004; Henon 1965a, b; Henon
   1973 (vertical stability, A&A 28); Bruno 1972, 1993, 1996 KIAM preprints; Bruno 1998/2000 (Power
   Geometry); Abalakin et al. 1971.
