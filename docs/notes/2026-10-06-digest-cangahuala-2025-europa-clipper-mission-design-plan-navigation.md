# Digest: Cangahuala et al. 2025, "Europa Clipper Mission Design, Mission Plan, and Navigation" (#960, #943)

L. A. Cangahuala, S. Campagnola, B. K. Bradley, D. R. Boone, B. B. Buffington, J. M. Ludwinski, S. Nandi
(JPL) & C. J. Scott (APL), "Europa Clipper Mission Design, Mission Plan, and Navigation", Space Science
Reviews 221:22 (2025), doi 10.1007/s11214-025-01140-2. Open access (CC BY 4.0).
- Filed as `cyclers_pdf/papers/cangahuala-campagnola-et-al-2025-europa-clipper-mission-design-plan-navigation-space-sci-rev-221-22-doi-10.1007-s11214-025-01140-2.pdf`.
  53 pages, text layer, md5 0484d448254f5fc5125838c1a0a9f1e9.
- Scope: I read sec. 2.4-2.5 (pp.10-16, Jupiter arrival and the 21F31_V6 tour). I read Table 2 (p.13,
  the targeted-flyby list) on the page image at 300 dpi; its text layer is unusable. I searched the full
  text for cycler and Ganymede-Callisto content. The mission plan and navigation sections (3-5) were
  skimmed.

## 0. Verdict

The flown reference tour 21F31_V6 has **no cycler**. The paper names cyclers only as one option for the
hemisphere switch (sec. 2.5.4, p.14): "petal rotation using Europa flybys ..., Ganymede flybys, or Callisto
flybys; different types of cyclers; a Europa pi-transfer; or a 'switch-flip.' For 21F31_V6, a petal
rotation with Callisto flybys is used". It cites Campagnola et al. 2019 (held), which is where the
Ganymede-Callisto (GCGC) cycler option is described.
- **`#943` cell gc:** no literal collision. The switch subphase (Table 2) is labelled "Ganymede-Callisto
  Petal Rotation". It is mainly Callisto non-resonant transfers with two Ganymede flybys, and it does not
  repeat.

## 1. Content (READ)

- **Arrival (sec. 2.4, pp.10-11):**
  - Ganymede flyby at 200 km, then about 12 h later JOI: a 900 m/s deterministic finite burn of 6.1 h
    at a periapsis range of 11.27 RJ, into a 203-day capture orbit.
  - Periapsis raise maneuver of 46-53 m/s, targeting G01 at 8.3 km/s.
- **Tour (sec. 2.5, pp.11-16):** 21F31_V6 has 53 Europa flybys over 4.3 yr, plus 7 Ganymede and 9
  Callisto flybys. Maximum inclination 7.5 deg; TID 2.97 Mrad. The phases are:
  - pump-down: G01-G03, E01-E02, G04-G05;
  - Europa Campaign 1 (26 Europa flybys, anti-Jovian, crank-over-the-top sequences);
  - the transition (the G-C petal rotation);
  - Europa Campaign 2 (leading-hemisphere petals of 5:1+ and 4:1-, then COT-3 to COT-5);
  - Ganymede impact (G08) for disposal.
  - The final pre-launch tour (21F31_V7) is "almost identical"; details are in Campagnola et al.
    2025a, b (JAS; not held).
- **Table 2, transition "Ganymede-Callisto Petal Rotation" (p.13, page image):**

| Enc. | In/Out | Date (TDB) | Alt, km | v_inf, km/s | Orbit res. | Days to next |
|---|---|---|---|---|---|---|
| 34C01 | O | 29-Sep-2032 15:52:37 | 1332.99 | 4.81 | 1.73:1.46 | 28.81 |
| 36G06 | I | 28-Oct-2032 11:13:26 | 1653.9 | 5.35 | 1.31:0.83 | 9.37 |
| 37C02 | I | 06-Nov-2032 20:06:00 | 203.2 | 3.59 | 1:1 | 16.68 |
| 38C03 | I | 23-Nov-2032 12:28:04 | 4641.32 | 3.58 | 1.55:1.55 | 25.87 |
| 39C04 | O | 19-Dec-2032 09:24:05 | 603.79 | 3.59 | 1.27:1.27 | 21.15 |
| 41C05 | I | 09-Jan-2033 12:59:43 | 1460.55 | 3.58 | 2.59:2.59 | 43.17 |
| 43C06 | O | 21-Feb-2033 17:05:07 | 10055.64 | 3.53 | 1:1 | 16.68 |
| 44C07 | O | 10-Mar-2033 09:30:36 | 3064.5 | 3.56 | 1.49:1.81 | 24.79 |
| 46G07 | O | 04-Apr-2033 04:26:39 | 969.26 | 4.3 | 3.5:1.18 | 25.06 |
| 47C08 | O | 29-Apr-2033 05:59:00 | 627.54 | 5.12 | 1:1 | 16.68 |
| 48C09 | O | 15-May-2033 22:24:18 | 3678.53 | 5.1 | 0.71:0.57 | 11.74 |

  - The numbers were read at 300 dpi. Re-read them on the page before using them as a control.
  - The column labels are as printed ("V-inf.", "Orbit Res.", "Enc. Dur." in days).

## 2. Gate relevance (`#943` cell gc)

- The transition is a Callisto petal rotation bracketed by Ganymede flybys G06 and G07. It is one-shot,
  not a cycler: the Callisto v_inf changes from 4.81 to 3.5x to 5.1 km/s, and Ganymede appears only twice
  in 7.5 months.
- **Against gc-1** (2.397/1.807) and **gc-2** (3.617/3.039): the Callisto v_inf of 3.53-3.59 km/s in
  the middle stretch is 0.5 km/s above gc-2's and 1.7 km/s above gc-1's. The Ganymede v_inf (5.35,
  4.30) is well above both. No collision.
- The "different types of cyclers" option points to Campagnola et al. 2019 (GCGC, held), which sec. 6.11
  of the `#942`/`#943` note has already excluded.

## 3. Positive controls

- Table 2 (the full 21F31_V6 flyby list) after a page re-read. It is not needed for gc.

## 4. Citation mining

Not mined in full. New items for the wanted list:
- Campagnola et al. 2025a, b (JAS; the 21F31 reference tour; "Design of the 21F31 Reference Tour" is
  already in Tier A).
- Scott et al. 2017 and 2025 (capture and pump-down; the 2025 JAS paper is doi 10.1007/s40295-025-00534-2,
  noted in `2026-10-06-943-gc-prior-art-search.md`).
- Campagnola et al. 2023 (eclipse-driven crank direction).
