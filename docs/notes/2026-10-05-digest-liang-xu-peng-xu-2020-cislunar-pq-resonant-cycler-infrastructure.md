# Digest: Liang, Xu, Peng & Xu 2020, "A cislunar in-orbit infrastructure based on p:q resonant cycler orbits" (#960)

Yuying Liang, Ming Xu (Beihang University), Kun Peng (CAST) and Shijie Xu (Beihang), "A cislunar in-orbit
infrastructure based on p:q resonant cycler orbits", Acta Astronautica 170:539-551 (2020), DOI
10.1016/j.actaastro.2020.02.029.
- Crossref-confirmed 2026-10-05: 170, pp. 539-551. Received 24 Jan 2019, accepted 16 Feb 2020.
- Filed as `cyclers_pdf/papers/liang-xu-peng-xu-2020-cislunar-in-orbit-infrastructure-pq-resonant-cycler-orbits-acta-astro-170-539-doi-10.1016-j.actaastro.2020.02.029.pdf`.
  13 pages, text layer, md5 df0a02c0d38d2d17e7e19a5cabfed7c7.
- I read all pages from the text layer.

**Name collision:** the first author is Yuying Liang (Beihang). The held `liang-2024-callisto-ganymede-
europa-triple-cyclers-JGCD.pdf` and the 2026 ice-giant review are by Guoliang Liang (NUAA). They are
different people and different groups. Do not merge their citation keys.

Evidence tags: READ, COMPUTED (my check, 2026-10-05, scratch script, not committed), INFERRED.

## 0. Verdict

This is a mission-concept paper, a cislunar "railway system" of cargo tugs parked on a p:q resonant orbit.
- The orbit work is:
  - p:q resonant orbits from a two-body guess, Eqs. (11)-(12).
  - Differential correction in the planar CR3BP.
  - Refinement in the bicircular model (BCM) by minimum-norm multiple shooting, where they are only
    quasi-periodic.
- Only ONE numeric member is given: a 2:1 resonant orbit with a two-number initial state.
- The family is Casoliva et al. 2010's Class 1 resonant family (2,1). The paper cites Casoliva for the
  definition.
- Its closest lunar approach, about 81,000 km, is outside the lunar SOI. So under the project's schema it is
  a `resonant_po`, not a lunar-encounter cycler.
- It does not duplicate any catalogue row.

## 1. Content (READ)

- Models (pp.541):
  - PCR3BP with the Earth at (-mu, 0) and the Moon at (1 - mu, 0).
  - BCM with the Sun's mass 328900.55, angular velocity 0.925196 and distance 388.811143 (normalised),
    with the Sun "assumed to rotate clockwise", Eq. (1).
- Multiple shooting (pp.541-542): N patch points. The underdetermined Newton step is solved as a
  minimum-norm problem, ||dQ||^2 with a Lagrange multiplier, using a block-tridiagonal recursion, Eqs.
  (5)-(10), to a tolerance of 1e-10.
- Resonance (p.542): following Casoliva et al. [19], "q T_M = p T_s"; "the spacecraft traverses an
  (inertial) elliptical orbit p times, while the Moon completes q revolutions", with a_s = (q/p)^(2/3)
  (Eq. 11).
- Geometric family bounds, Eq. (12): (1 + R_M)(p/q)^(2/3) - 1 < e < 1 - R_E (p/q)^(2/3).
  - Safety margins: 2000 km above the Earth and 500 km above the Moon.
  - The apogee is on the far side of the Moon.
  - "a family of planar resonant cycler orbits can be defined by only three parameters: the resonance
    relation, the eccentricity e, and the argument of perigee omega" (p.542).
- **The numeric member (p.542):** "Its initial position is 80927.9 km away from the Moon, and the initial
  velocity is 1.102 km/s perpendicular to the Earth-Moon line."
  - It was chosen "near the red one" of Fig. 2, the 2:1 member that exactly meets the Earth margin.
  - No frame, sign or epoch is printed for this state. Fig. 3(a) shows it in the PCR3BP.
  - Fig. 3(b) shows it refined in the BCM over 10 periods. It is "no longer strictly periodic in the BCM
    but still preserves the safety margin", with ||F|| < 1e-10 after 9 iterations (p.543).
- Mission numbers, which are not cycler invariants:
  - Table 1, modified-Lambert transfer costs: EPO -> SRCO 1660.7, SRCO -> MPO 2369.4, MPO -> SRCO 1392.1,
    SRCO -> EPO 2434.3 m/s.
  - Table 2, end-to-end costs: Earth -> Moon 4066.9 m/s over 18.12 d; Moon -> Earth 4296.7 m/s over
    17.55 d.
  - Launch and rendezvous windows: Figs. 7-17.

## 2. My check of the printed state (COMPUTED)

Method:
- Planar CR3BP with the project registry mu = 0.01215058439, 384,400 km and 375,190.26 s.
- Start at x0 = 1 - mu + 80927.9/384400 (far side), y0 = 0, xdot0 = 0, with a rotating-frame ydot0.

Results:
- With ydot0 = -1.102 km/s (-1.0756 nondim), the trajectory returns to y = 0 at t = 3.108, x = -1.2575,
  with xdot = 0.032. That is near, but not at, a perpendicular crossing. The printed state does not close as
  a symmetric periodic orbit to the printed precision in this model.
- The inertial-velocity reading and the +1.102 sign do not give a 2:1 return at all.
- Correcting ydot0 alone for a perpendicular half-period crossing gives:

| Quantity | Value |
|---|---|
| ydot0 | -1.04414 nondim = -1.0698 km/s |
| Half period | 3.0949, so T = 6.1898 nondim, about 26.9 d |
| Jacobi constant C | 2.0934 |
| Perigee radius | 13,872 km (altitude 7,494 km) |
| Apogee radius | 480,145 km |
| Minimum Moon distance | 80,928 km, at the start point |

- The corrected member's perigee is far above the paper's 2000 km "red" margin. So either the printed
  velocity uses a different frame or convention, or the chosen member is not as close to "the red one" as
  the text implies.
- INFERRED: the 3 percent velocity mismatch is too large for rounding. This is a convention gap, not an
  erratum claim.
- Use the printed two numbers only as a loose shape target. Use my corrected member only as a COMPUTED
  stand-in, not as a golden.

## 3. Gate answers: literal-collision checks

- **Families and model:** only the (2,1) Casoliva Class 1 family in this paper, with one numeric member
  (above).
  - Eq. (12) bounds the general p:q family, and Fig. 2 shows two 2:1 members at the margins, as figures.
  - Model: planar CR3BP, with quasi-periodic BCM versions.
  - No other p:q value is computed.
- **Catalogue duplication (read only):** no row duplicates it.
  - There is no `liang-2020` row.
  - The nearest rows are `casoliva-2-1a-em-resonant-po-2010` (C = 0.4887, periselene 90,471 km, perigee
    72,142 km) and `casoliva-2-1b-em-resonant-po-2010` (C = 1.1964, periselene 92,590 km, perigee 6,790 km,
    apogee 485,181 km). Both are `resonant_po` and both are the same (2,1) family at different energies.
  - The Liang member (C about 2.09 by my check) would be a third point on that family, with the same
    `resonant_po` verdict, because 80,928 km is above the 66,183 km SOI.
- **`#947` R3 (Zhou et al. 2025 Table 6 "cycler-like" fixed points): no literal collision.**
  - Zhou's n = 2 row is x0 = 0.693, with xdot0 = -0.0212 (not perpendicular), on the Earth side.
  - Liang's member is a symmetric far-side 2:1 orbit.
  - Add the Liang 2:1 signature (C about 2.09, far-side perpendicular crossing at x about 1.198) to the R3
    match list alongside the Casoliva rows.
- **`#948` R4 (Earth-grazing second-species): no collision.**
  - The Liang orbits are high-energy, Keplerian-like and far from the Moon (no near-collision).
  - The Earth-margin end of the family (2000 km altitude, Fig. 2 red) is an Earth-grazing resonant orbit of
    the second KIND, not second species.
- **`#956` R9 (exterior-realm cyclers through L2): no literal collision.**
  - The Liang member's apogee is beyond the Moon, but C about 2.09 is far below C(L2), about 3.17. The Hill
    region is wide open, so this is not an L2-tube transit, which is the R9 class of Ross &
    Roberts-Tsoukkas.
  - It is still published prior art for "p:q orbits with apogee beyond the Moon" at high energy. Cite it if
    an R9 result is presented at C near 2.
- **Prior-art lineage the paper records:** Earth-Moon cyclers "first explained by Broucke in 1968 [4]". It
  does not mention Arenstorf 1963's ferry concept (see the Arenstorf digest).

## 4. Relation to the held Liang papers

- Guoliang Liang et al. 2024 (CGE triple cyclers) and the 2026 ice-giant review are a different first
  author (see the name collision above). Neither cites this paper (not checked in full). There is no
  method overlap beyond multiple shooting.
- Yuying Liang's own earlier paper, Liang, Xu & Xu 2017 (Acta Astro 133:282, cislunar polygonal-like
  periodic orbits, ref. [18]), is not held (below).

## 5. Positive controls

- Printed: the 2:1 member. Its distance of 80,927.9 km from the Moon is the most reliable number. The
  velocity of 1.102 km/s is frame-unclear (sec. 2).
- Tables 1-2 transfer costs, which are mission numbers and not invariants.

## 6. Citation mining (policy step 4)

The references [1]-[30] (pp.550-551) were checked against `CORPUS_INDEX.md`, the digests and the
filenames.

Held:
- [4] Broucke 1968.
- [9] Russell & Strange 2009.
- [13] Byrnes, Longuski & Aldrin 1993.
- [15] McConaghy, Longuski & Byrnes 2004.
- [16] Rogers et al. 2015.
- [19] Casoliva et al. 2010.
- [20] Szebehely 1967.

Not held, in priority order. The DOI was checked through the Crossref API on title, authors, volume and
pages.
1. Liang, Y., Xu, M. & Xu, S. (2017), "The cislunar polygonal-like periodic orbit: Construction, transition
   and its application", Acta Astronautica 133:282-301, doi 10.1016/j.actaastro.2017.01.028 (CONFIRMED).
   Earth-centred resonant periodic orbits and their transition conditions; R3 and R9 collision context.
2. Anderson, R. L., Campagnola, S. & Lantoine, G. (2016), "Broad search for unstable resonant orbits in the
   planar circular restricted three-body problem", CMDA 124:177-199, doi 10.1007/s10569-015-9659-7
   (CONFIRMED; Crossref issue year 2015). A resonant-orbit enumeration method.
3. Gomez, G. & Mondelo, J. M. (2001), "The dynamics around the collinear equilibrium points of the RTBP",
   Physica D 157:283-321, doi 10.1016/S0167-2789(01)00312-8 (CONFIRMED). The multiple-shooting method source.
4. Pelle, S. et al. (2019), "Earth-Mars cyclers for a sustainable human exploration of Mars", Acta
   Astronautica 154:286-294, doi 10.1016/j.actaastro.2018.04.034 (CONFIRMED).
5. Lower priority (Earth-Moon mission background, free-return trajectories, libration-point constellations
   and reports; DOIs not checked):
   - Farquhar & Dunham 1981, JGC 4:192.
   - Bao, Li & Baoyin 2018, ASR 61:97.
   - Luo, Yin & Han 2013, JGCD 36:263.
   - Xu, Liu & Xu 2013, MPE.
   - Burns et al. 2013, ASR 52:306.
   - Zimmer 2013, Acta Astro 90:119.
   - Lo, PAMM 2007.
   - Aldrin 1985 SAIC presentation (no DOI).
   - Niehoff & Friedlander 1985 study (no DOI).
   - Xu, Tan & Xu 2014.
   - Biesbroek & Janin 2000.
   - McVay et al. 2016.
   - Gill 2019 NASA report.
   - Koon et al. 2006 / Wang et al. book (duplicate refs. [21] and [29]).
   - Kuehn 2015.
   - Matsumoto 2006; Wang & Liu 2016; Anselmo 2000.

## 7. Novelty-gate anchor: deliberately NOT added

- I tried a `KNOWN_CORPUS` anchor for this paper (body_set {Moon}, authors "Liang", "Xu", "Peng") and
  removed it.
- With it, `tests/search/test_literature_check.py::test_new_corpus_entries_flagged_published` failed. The
  matcher strong-linked the fixture's Guoliang Liang 2024 JGCD hit (doi 10.2514/1.G008387) to this Yuying
  Liang anchor through the shared surname. That is the name collision from the top of this note, now
  reproduced in code.
- The existing Casoliva 2010 anchor already covers the (2,1) Class 1 family.
- If an anchor is wanted later, the matcher first needs author disambiguation, for example by first name or
  DOI.
