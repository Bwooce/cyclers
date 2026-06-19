# #345 digest — Voyager (Kohlhase-Penzo 1977) + Mariner-10 (Giberson-Cunningham 1975)

**Date:** 2026-06-19 AET
**Task:** #345 / #387 classic-mission mga_tour admission candidates
**Verdict (both):** NEGATIVE for V0 mga_tour catalogue admission — no per-encounter V∞ published. Filed as supporting references; V∞ routes to the dedicated navigation literature or NAIF SPK direct evaluation.

Read every page of both (directly, in the main session — subagent quota
partially exhausted). Both confirm the now-firm #345 pattern: mission-
description / project-management papers publish encounter dates +
closest-approach geometry + maneuver ΔV budgets + science, but NOT the
hyperbolic-excess-velocity (V∞) tuples the §14 V0 mga_tour standard
requires (set by Heaton-Longuski 2003 U00-01, confirmed-negative by
Wolf-Smith 1995 Cassini #355, Dunne-Burgess 1978 Mariner SP-424 #382,
Bourke 1971 Voyager #382, Lam 2008 Juno #382, Bellerose 2018 Cassini #382).

## Kohlhase & Penzo 1977 — "Voyager Mission Description"

* **Citation:** Space Science Reviews 21(2):77-101 (1977). C.E. Kohlhase,
  P.A. Penzo (JPL/Caltech). Cached: `cyclers_pdf` commit `26e66d0`.
* **What it is:** the canonical Voyager (MJS77) flight-design summary —
  the flown two-spacecraft Jupiter-gravity-assist-to-Saturn design with
  the Voyager 2 Saturn-gravity-assist continuation option to Uranus.
* **Trajectory data published (Table IV p.90, Figs 6-11):**
  - JST: launch 1977-09-01, Jupiter 1979-03-05 (4.9 R_J), Saturn 1980-11-13
  - JSX(→Titan): launch 1977-08-20, Jupiter 1979-07-09 (10.0 R_J), Saturn 1981-08-27
  - JSX(→Uranus): same launch, Saturn 1981-08-27, Uranus 1986-01-30
  - Per-satellite closest-approach radii + encounter times (Io/Europa/
    Ganymede/Callisto; Titan/Rhea/Tethys/Enceladus/Dione/Mimas)
  - Post-Jupiter ΔV: 60 m/s (JSX) / 100 m/s (JST)
  - Occultations, comm distances, one-way light times
* **V∞: ABSENT.** Table IV has no V∞ column; closest-approach is in
  planet radii, not asymptotic speed. Figs 7-10 give km closest-approach
  distances + satellite periods only.
* **The real V∞ source (from this paper's own reference list, p.101):**
  - **McKinley, E.L. & Van Allen, R.E. 1976**, "Mariner Jupiter-Saturn
    1977 Navigation Strategy," *J. Spacecraft & Rockets* 13:494-501 —
    the navigation paper most likely to carry encounter V∞.
  - Jacobson, Synnott, Dixon 1976 (nav accuracy analysis)
  - Bourke, Miles, Penzo, Van Dillen, Wallace 1972 (AIAA 72-943, MJS
    preliminary mission design)
  - Penzo 1974 (AIAA 74-780, science objectives + mission design)

## Giberson & Cunningham 1975 — "Mariner 10 mission to Venus and Mercury"

* **Citation:** Acta Astronautica 2(7-8):715-743 (1975), DOI
  10.1016/0094-5765(75)90012-0. W.E. Giberson (Project Manager),
  N.W. Cunningham (Program Manager), JPL. Cached: `cyclers_pdf` commit
  `28d8531`.
* **What it is:** the contemporaneous JPL project-manager mission paper
  for MVM'73 / Mariner 10 — the first dual-planet gravity-assist mission
  (Venus → Mercury, with the 2:1 resonant Mercury re-encounters).
* **Trajectory data published (text + Figs 1-2, 8-9):**
  - Launch 1973-11-03, Atlas/Centaur, launch weight 499 kg
  - Venus flyby 1974-02-05, 5000 km altitude (aimed 4800 km)
  - Mercury I 1974-03-29/30, 1000 km altitude (aimed 960 km)
  - Mercury re-encounter 176 days after the first
  - TCM-1 7.8 m/s (L+10 d); TCM-2 1.37 m/s (1974-01-21); total
    propulsion capability 119 m/s
* **V∞: ABSENT.** Dates + flyby altitudes + TCM ΔVs + extensive Venus/
  Mercury science, but no per-encounter V∞.
* **The real V∞ source (from this paper, p.716 + refs):**
  - **Bourke & Beerer 1970** — the baseline MVM'73 preliminary mission
    design (cited p.716 as the design baseline; this vindicates the
    earlier #387 Bourke-Beerer pointer that had been flagged unconfident)
  - Minovitch 1963 (the foundational gravity-assist analysis)
  - Sturms 1965 / 1966a / 1966b; Eckman 1969 (trajectory studies)

## Net #345 status after this digest

Both Voyager and Mariner-10 join the #345 confirmed-negative set: the
mission-overview literature is insufficient for V0 admission. The
catalogue cannot admit Voyager or Mariner-10 mga_tour rows until either:
1. A dedicated navigation paper with a per-encounter V∞ table is acquired
   (McKinley-Van Allen 1976 JSR for Voyager; Bourke-Beerer 1970 for
   Mariner-10), OR
2. The **NAIF SPK direct-evaluation path** is built — derive V∞ at each
   encounter from the mission's archived SPK kernels (same fallback noted
   for Cassini #361). This is a methodology build, not an acquisition;
   it would unblock Voyager 1/2, Mariner-10, Galileo, Cassini, Pioneer,
   Juno all at once and is probably the highest-leverage way to clear the
   #345 backlog wholesale.

The two papers ARE valuable as catalogue corroborating-source / KNOWN_CORPUS
context for the missions; only the V0 numeric-admission bar is unmet.
