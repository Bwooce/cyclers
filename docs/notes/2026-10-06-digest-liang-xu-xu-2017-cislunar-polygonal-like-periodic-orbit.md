# Digest: Liang, Xu & Xu 2017, "The cislunar polygonal-like periodic orbit: Construction, transition and its application" (#960 batch 23; #947 R3)

Yuying Liang, Ming Xu and Shijie Xu (Beihang University), Acta Astronautica 133:282-301 (2017), doi
10.1016/j.actaastro.2017.01.028.
- Filed as `cyclers_pdf/papers/liang-xu-xu-2017-cislunar-polygonal-like-periodic-orbit-construction-transition-application-acta-astro-133-282-doi-10.1016-j.actaastro.2017.01.028.pdf`.
  20 pages, text layer, md5 00b64d09bf4179ceb3bfa19d5a7685fa. Supplied by the owner.
- I read secs. 1, 3.2 (result), 4-6, Tables 1-2 and the references from the text layer. Sec. 2 (the
  Hori-Lie machinery) was skimmed. The numbers below are from the text layer; they are simple decimals,
  but read the page before using any as a golden value.

## 0. Verdict

The direct predecessor of the held Liang, Xu, Peng & Xu 2020 (Acta Astro 170:539, p:q resonant cycler
infrastructure).
- It defines **polygonal-like periodic orbits (PLPOs)**: m:n Earth-Moon resonant orbits, in the planar
  CR3BP rotating frame, whose apses form a polygon.
  - They are constructed analytically with lunar gravity treated as a perturbation (Hori-Lie
    averaging).
  - The periodicity condition is in three modified Kepler elements (a, e, omega).
- **Representatives (sec. 4, Fig. 5), in Earth-Moon units:**
  - **5:2 PLPO:** a = 0.5620, e = 0.4660, omega = 1.3634, C = 3.0996.
  - **7:3 PLPO:** a = 0.5604, e = 0.3430, omega = 4.1493, C = 3.1858.
- **These are NOT lunar-flyby cyclers.**
  - Apoapsis a(1 + e) is about 0.82 (5:2) and 0.75 (7:3), inside the Hill-region boundary. They stay
    in the Earth's region.
  - A small rotation in omega (0.0674 rad for 5:2, 0.0034 rad for 7:3) breaks the "frozen tori" and
    gives cislunar transit orbits.
  - The "parking apron" concept uses a PLPO as a dormant orbit, and a degenerated transit orbit as the
    working orbit to the Moon.
- **#947 R3 (Zhou 2025 Table 6 "cycler-like" fixed points):** add the 5:2 and 7:3 PLPOs to the match
  list as m:n resonant Earth-region orbits. They are not Moon-encountering, so a Zhou row that does
  encounter the Moon cannot be a PLPO. **No collision expected; check (a, e, C) when R3 runs.**

## 1. Content (READ)

- **Method (secs. 2-3):**
  - Kepler-element initial conditions, transformed to rotating Cartesian coordinates (sec. 3.1).
  - Lunar gravity is averaged with the Hori-Lie perturbation method.
  - Periodicity conditions, eqs. (23)-(24), in (a, e, omega).
- **Transition (sec. 4):** Poincaré sections Sigma1 (fixed a0, e0, theta0) and Sigma2 (apsides,
  r_e . r_e-dot = 0), each propagated for 10 yr.
  - 5:2: frozen tori for omega in [1.365, 1.395]; the torus breaks at 1.408; critical non-transit at
    1.417; transit above 1.425, and below 1.265. The critical apoapsis line is a(1 + e) = 0.83.
  - 7:3: frozen tori for omega in [4.1500, 4.1513]; critical omega 4.1520; transit-set boundaries
    4.1475 / 4.1520; a(1 + e) = 0.78.
  - So the 7:3 orbit is more sensitive than the 5:2. Sec. 4.3 gives a five-critical-value description
    (omega_1-omega_5).
- **Parking apron (sec. 5):**
  - The dormant orbit is the 5:2 PLPO (a = 0.5620, e = 0.4660, omega = 1.3634, theta = 0). The working
    orbit is the transit orbit at omega = 1.0461, with 6 loops before entering the Moon's region.
  - Table 1 gives the dV and TOF from the 7 intersection points: 47.1-610.1 m/s and 7.4-77.6 d.
  - Table 2 gives Hohmann transfers from a 200 km LEO to 5 apoapsides: dV 3.25-3.35 km/s,
    TOF about 3.88-3.90 d.
  - The paper contrasts this with the "Cislunar Bus System" of Xu & Liang (IAC-13), which used Lambert
    transfers costing about 3 km/s to de-orbit.

## 2. Relation to the held Liang papers

- **Liang, Xu, Peng & Xu 2020** (held; digest `2026-10-05-digest-liang-xu-peng-xu-2020-...`): extends
  from PLPOs (Earth-region only, lunar gravity as a perturbation) to p:q resonant orbits that do
  encounter the Moon (cyclers). The 2020 digest already cites this 2017 paper as its predecessor.
- **Liang 2024 JGCD** (Callisto-Ganymede-Europa triple cyclers) and **Liang et al. 2026** (ice-giant
  review): same group, different systems. This 2017 paper adds no Jovian content.

## 3. Citation mining (refs 1-32)

Held:
- Barrabés & Gómez 2002 and 2003 [10, 11].
- Antoniadou & Voyatzis [13, 14] (the 2013 and 2014 items).
- Vaquero 2013 PhD (covering [17, 22, 24] in part).
- Lo & Parker 2004 [18].
- Casoliva 2008 AIAA [31].
- Font-Nunes-Simó 2002 [9].
- Szebehely 1967 [28].

Not held:
- **[32] Xu, M. & Liang, Y., "A Cislunar In-Orbit Infrastructure Using Cycler Trajectories in the
  Earth and Moon System", IAC-13-A3.P.57 (2013).** Earth-Moon cyclers ("Cislunar Bus System"). Added
  to the wanted list for the `#947` match list.
- [4] Carrico et al., lunar-resonant trajectory design; [5] McComas et al. 2011 (IBEX long-term stable
  orbit); [15] Dichmann et al., 3:1 resonance dynamics. All three are IBEX/TESS P/2 and 3:1 lunar
  resonance; one low-priority row.
- [12] Barrabés & Gómez 2004, "A note on second species solutions generated from p-q resonant
  orbits", CMDA 88(3):229-244: not held (the held Barrabés-Gómez files are the 2002 and 2003 papers).
  Added to the same low-priority row.
- [26] Belbruno & Marsden 1997 (AJ 113, resonance hopping in comets); [27] Yang & Baoyin (Kleopatra):
  background, not listed.
