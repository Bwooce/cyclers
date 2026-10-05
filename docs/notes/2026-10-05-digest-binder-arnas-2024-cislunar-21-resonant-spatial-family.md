# Digest: Binder & Arnas 2024, "Reliable and Repeatable Transit Through Cislunar Space Using the 2:1 Resonant Spatial Orbit Family" (#960)

A. Binder and D. Arnas (Purdue University).
- Journal version: J. Guidance, Control, and Dynamics 47(9):1973-1979 (2024), doi 10.2514/1.G007800
  (Crossref-confirmed by the lead).
- The filed copy is arXiv 2304.13584v2, 34 pages, text layer, md5 65112e2fec9d24cd54f40bf7cc8de4b6, as
  `cyclers_pdf/papers/binder-arnas-2024-reliable-repeatable-transit-cislunar-21-resonant-spatial-family-arxiv-2304.13584v2-jgcd-47-9-1973-doi-10.2514-1.G007800.pdf`.
- I read the text layer in full. Page cites are arXiv pages. The 7-page JGCD version may be condensed; I
  did not check it.

## 0. Verdict

- The paper publishes a SPATIAL Earth-Moon 2:1 resonant family, R2:1-S. It bifurcates by period doubling
  from both the planar prograde 2:1 family (R2:1-P) and the planar retrograde 2:1 family (R2:1-R), and it
  links them continuously in three dimensions (Sect. IV.B, p.12-13).
- Members pass close to both the Earth and the Moon: an Earth pass about every 14 d and a Moon pass about
  every 28 d. The synodic period is about twice the lunar period, because of the period doubling.
- Some members of the R2:1-S2 half are linearly stable ("marginally stable", Fig. 11).
- Exactly ONE numeric member is printed, Table 2 (p.26):

| x0 (km) | z0 (km) | ydot0 (km/s) | C | P_syn (d) | stability index nu |
|---|---|---|---|---|---|
| 368966 | 15789 | -0.82796 | 2.78959 | 54.47 | 1418.34 |

- It is measured at a crossing of surface A (the y = 0 plane between L1 and the Moon, z > 0), with
  xdot = zdot = 0 at the crossing.
- Model constants (Table 1, p.3): mu1 = 398600.4415, mu2 = 4902.8005821478 km^3/s^2, l* = 384400 km,
  t* = 375190.26 s, mu = 0.012150585.
- Everything else is figure-only: hodographs Figs. 9-11, and the Broucke stability diagram, Fig. 7.

## 1. Gate answers

- **`#957` R10 (spatial second-species seeds at a physical mass ratio):**
  - Binder & Arnas is published prior art for spatial 2:1 Earth-Moon resonant orbits that pass near both
    primaries.
  - Any R10 spatial 2:1 Earth-Moon member must be checked against R2:1-S. Use the Table 2 state as an
    anchor point, and the family's period-doubling origin from the planar 2:1 families.
  - R2:1-S members are high-energy (C about 2.79 for the printed member). They are not second-species
    (near-collision) orbits as constructed, but the family does pass the Moon at a range of distances
    (Fig. 10b).
- **`#947` R3 (Zhou 2025 Table 6 "cycler-like" planar fixed points):**
  - The planar parents R2:1-P and R2:1-R (Arenstorf 1963, Broucke 1968) are the planar families involved.
  - The paper prints no planar initial conditions. It is a family-identity check for R3, not a numeric one.
- **Catalogue (read only):**
  - The row `em-cycler-21-3d-spatial-2026` is a project-computed spatial (2,1) Earth-Moon orbit. It has
    C = 3.0258, period 18.17 nondimensional, and comes from a vertical-critical bifurcation of the
    low-energy Ross-RT (2,1) planar family.
  - Its data-gap note says "No published IC for THIS Earth-Moon spatial member (Antoniadou & Libert 2019 is
    mu=0.001 ...)".
  - Binder & Arnas is a published Earth-Moon spatial 2:1 family. It is probably a DIFFERENT family: high
    energy, period-doubling origin, C about 2.79 against 3.03 (INFERRED from the printed C and origin).
  - It should still be cited in that row's class attribution as Earth-Moon prior art. Proposal only; no
    catalogue write.

## 2. Positive control

- Table 2 state above (p.26). Propagate in the CR3BP with the Table 1 constants. A closure after the
  54.47 d synodic period, and nu of about 1418, are the checks.
- Convert x0 from km to nondimensional units with l* = 384400 km. The text implies the origin is at the
  barycentre: surface A lies between L1 and the Moon, at about 0.836-0.988 LU. 368966/384400 = 0.95985
  (COMPUTED), which is inside that range, so the reading is consistent.

## 3. Citation mining (references [1]-[34])

Held:
- [8] Arenstorf 1963 AIAA J (filed by `#960`).
- [9] Broucke 1968.
- [17] Szebehely 1967.
- [27] Koon et al. book.
- [34] Byrnes et al. 1993.

Not held, in priority order. DOIs were not checked unless noted.
1. [14] Gupta, M. (2020), "Finding Order in Chaos: Resonant Orbits and Poincare Sections", MS thesis,
   Purdue. The planar 2:1 families in detail.
2. [12] Vaquero & Howell (2014), Acta Astronautica 94(1):302-317, doi 10.1016/j.actaastro.2013.05.006
   (printed in the reference list); [13] Vaquero (2013), PhD thesis, Purdue.
3. [15], [16] Gupta, Howell & Frueh (2021, 2022), resonant orbits for cislunar surveillance.
4. [10], [11] de Almeida Prado & Broucke (1995, 1996), JGCD 18(3):593 and 19(4):929, dois
   10.2514/3.21428 and 10.2514/3.21720 (printed in the reference list).
5. [24] Campbell (1999), PhD thesis, Purdue, bifurcations in the CR3BP; [25] Grebow (2006), MS thesis,
   Purdue; [33] Haapala & Howell, higher-dimensional Poincare maps; [32] Gomez et al. 2004, connecting
   orbits and invariant manifolds; [22] Broucke 1969 AIAA J (held:
   `broucke-1969b-stability-...-aiaa-j-7-1003`); [23] Howard & MacKay 1987.
6. Tethers and background, low priority: [29]-[31] tether papers; [1]-[7] programme documents.
