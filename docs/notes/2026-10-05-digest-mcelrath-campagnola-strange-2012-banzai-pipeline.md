# Digest: McElrath, Campagnola & Strange 2012, "Riding the Banzai Pipeline at Jupiter" (#960)

T. P. McElrath (JPL), S. Campagnola (ISAS/JAXA) & N. J. Strange (JPL), "Riding the Banzai Pipeline at
Jupiter: Balancing Low Delta-V and Low Radiation to Reach Europa", AIAA 2012-4809, AIAA/AAS Astrodynamics
Specialist Conference, 2012, doi 10.2514/6.2012-4809 (wanted-list CONFIRMED).
- Filed as `cyclers_pdf/papers/mcelrath-campagnola-strange-2012-banzai-pipeline-jupiter-low-dv-low-radiation-europa-aiaa-2012-4809-doi-10.2514-6.2012-4809.pdf`.
  13 pages. The upload (md5 09bb6a151ccabddca4345cfe7285f063) has a text layer except on pp.9-11, which
  are page images. I OCR'd those three pages (`ocrmypdf --skip-text`); the filed PDF carries that layer.
  Former wanted-list row 4.
- I read pp.1-8 and 12-13 from the text layer, and pp.6 and 9-11 on the page images. Table 1 (p.6)
  agrees between the text layer and the image.

## 0. Verdict

The "Banzai pipeline" is a ONE-WAY descent from Callisto to Europa orbit: C-C backflip, C-to-G inclined
Hohmann, G-G backflip, G-to-E inclined Hohmann, then Europa orbit insertion ("strictly a one-way descent",
p.4 footnote). Nothing repeats. So there is no literal collision for either `#943` X1 pair:
- **Ganymede-Callisto:** one G4-C5 transfer and one C-to-G Hohmann leg in the pipeline. No repeating G-C
  path.
- **Ganymede-Europa:** one G-to-E Hohmann leg. No repeating G-E path.

It is still relevant prior art for X1. It publishes (p.2-3) that Ganymede is the only Galilean moon that
can rotate a v_inf between the Callisto-Ganymede and Ganymede-Europa Hohmann transfers in a
flyby-backflip-flyby, and it gives the G-C and E-G-C phasing cycles.

## 1. Content (READ)

- **Geometry (Fig. 1, pp.2-3), patched conic with circular coplanar moon orbits:**
  - Inclined Hohmann transfers between moon pairs; their periapsis and apoapsis velocities trace arcs out
    of the orbit plane.
  - A Europa v_inf of 1.711 km/s at just over 22 deg to the orbit plane is the smallest that meets all
    the conditions.
  - The Ganymede backflip uses a single flyby at 100 km minimum altitude. The G-C end is "at the backflip
    limit": "If Ganymede was somewhat smaller ... this sequence would lose its attractiveness."
  - A backflip (Uphoff, Roberts & Friedman 1976) is a 180-deg transfer with an inclination change.
- **Flight order (p.4):** C-C-G-G-E, then EOI.
- **Phasing (pp.5-6, Fig. 5):**
  - E and G are near 2:1, with a -5.2 deg inertial offset per repeat. The E-G alignment direction turns
    a full circle in about 487 d (438 d Sun-relative).
  - G and C are near 7:3. At the E-G alignment, the G-C phase nearly repeats every 16 weeks (average
    drift 3.15 deg per cycle); 7- and 9-week repeats match to within +/-22 deg.
  - So a pipeline arrival is available every 7-9 weeks if a 22-deg G-C misalignment is acceptable, and
    about every 7.5 months within minor variations (p.7).
  - Check: Ganymede 7 x 7.155 d = 50.09 d; Callisto 3 x 16.689 d = 50.07 d. This matches the "offset
    pairs every 50 days" of G-out/C transfers (p.11).
- **Table 1 (p.6), full ephemeris, each case from Callisto to Europa:**

| Final G flyby | Europa arrival | Plane | v_inf at E, km/s | EOI to 100 km, km/s | G-C offset, deg |
|---|---|---|---|---|---|
| 18 Apr 2027 | S | down | 1.703 | 1.211 | -18.1 |
| 20 Jun 2027 | S | up | 1.665 | 1.186 | +5.1 |
| 11 Oct 2027 | N | down | 1.701 | 1.209 | +8.8 |
| 29 Nov 2027 | N | down | 1.720 | 1.222 | -11.9 |
| 1 Feb 2028 | N | down | 1.695 | 1.206 | +13.4 |
| 21 Mar 2028 | S | up | 1.654 | 1.179 | -7.7 |

  Offsets up to about 18 deg do not change the cost much.
- **Full tour (pp.8-12, Figs. 7-11):**
  - 20 months, with Ganymede flybys G0-G4 (805 m/s JOI, 110 m/s PJR).
  - G4-C5 transfer at about 18 Rj perijove; 30 krad to C5.
  - Callisto sequence C5-C14: 3:1, 2:1, 3:2, near-1:1, with ΔV leveraging. 139 m/s actual against 94 m/s
    phase-free.
  - Pipeline: C14 to G15 to G16 to EOI. G15 and G16 are 3.5 d apart.
  - Total dose 89 krad. Europa v_inf 1.71 km/s. 223 m/s from PJR to EOI, of which 100 m/s is the cost of
    real phasing.
- **Extensions (p.12):**
  - Not to Io: the Europa turn is too small.
  - The Uranian satellites are dynamically similar, but there is no radiation motive there.

## 2. Gate relevance

- **`#943` X1 Ganymede-Callisto:** no collision. Two useful facts for a G-C cycler check (INFERRED):
  - the 7:3 G-C commensurability gives a 50-d near-repeat of the G-C phase;
  - Ganymede's turn is the binding limit at the G-C Hohmann (Fig. 1).
- **`#943` X1 Ganymede-Europa:** no collision. The E-G alignment regresses through a full circle in about
  487 d. Any repeating G-E path must close against this drift.
- **Cite with:** Campagnola et al. 2019 (GCGC, held), Anderson et al. 2021 (Ganymede-Europa endgame,
  held), and Buffington, Campagnola & Petropoulos 2012 (filed today).

## 3. Positive controls

- Table 1: six full-ephemeris pipelines with Europa v_inf 1.654-1.720 km/s. The patched-conic ideal is
  1.711 km/s.
- Fig. 1 construction: E v_inf 1.711 km/s at just over 22 deg, with a 100-km Ganymede backflip limit.
  This is reproducible in a circular coplanar patched-conic model.

## 4. Citation mining (references 1-11, p.13)

Held:
- none directly.
  - Ref 7 (Campagnola, Russell & Skerritt, "Flybys in the PCR3BP", CMDA, to be published) is not held.
  - Ref 3's Petropoulos-Kloster-Landau AAS 09-354 is not held.

Not held:
1. Uphoff, C., Roberts, P. H. & Friedman, L. D. (1976), "Orbit Design Concepts for Jupiter Orbiter
   Missions", JSR 13(6):348-355, doi 10.2514/3.57096. The backflip origin. Already in the wanted list (Clipper and Jovian background row).
2. Campagnola, S., Russell, R. P. & Skerritt, P. (2012), "Flybys in the planar, circular, restricted,
   three-body problem", CMDA 113:343-368. DOI not checked. Multi-body flyby model (X1 context).
3. Heaton, A. F., Strange, N. J., Longuski, J. M. & Bonfiglio, E. P. (2002), "Automated Design of the
   Europa Orbiter Tour", JSR 39(1):17-22. DOI not checked.
4. Strange, N. J. & Longuski, J. M. (2002), "Graphical Method for Gravity-Assist Trajectory Design",
   JSR 39(1):9-16. DOI not checked. The Tisserand-graph method.
5. Johannesen & D'Amario (1999), AAS 99-360; Petropoulos, Kloster & Landau (2009), AAS 09-354;
   Grebow, Petropoulos & Finlayson (2011), AAS 11-427; Petropoulos, Longuski & Bonfiglio (2000), JSR
   37(6):776-783. Europa tour background.
6. Garrett et al. (2003), JPL 03-006 (radiation model); Tisserand (1896); Labunsky et al. (1998).
   Background.
