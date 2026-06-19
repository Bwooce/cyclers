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

## McKinley & Van Allen 1976 — "Mariner Jupiter/Saturn 1977 Navigation Strategy" (added 2026-06-19)

* **Citation:** J. Spacecraft & Rockets 13(8):494-501 (1976), DOI
  10.2514/3.57113. E.L. McKinley, R.E. Van Allen (JPL). Cached:
  `cyclers_pdf` commit `aca2a36`. The DEDICATED Voyager (MJS77)
  navigation paper — the Tier-1 V∞ source flagged in the section above.
* **What it is:** trajectory-correction-maneuver (TCM) strategy +
  propellant costs + delivery accuracies + planetary-quarantine analysis.
* **Data published:** Table 1 (delivery-accuracy reqs: Jupiter 1500/600 km,
  Saturn 4000/600 km), Table 2 (TCM execution errors), Tables 3-5 (TCM
  ΔV: Earth-Jupiter leg 53.8 m/s, Jupiter-Saturn 39.6 m/s; 155 m/s total
  nav allocation), Tables 6-7 (PQ contamination + Titan retarget), Fig 7
  (mission ΔV99 vs post-Jupiter deterministic ΔV), Fig 8 (propellant vs
  launch date). Saturn B-plane sensitivity 50,500 km/m/s.
* **V∞: ABSENT — decisive finding.** Even the dedicated navigation-strategy
  paper (the best-case venue for V∞) does NOT tabulate per-encounter
  hyperbolic excess velocity. Its V∞ values live in internal JPL documents
  it cites and we cannot obtain (Wallace 1974 EM 392-117; Curkendall 1974
  PD618-115; Gates 1963 TR 32-504). **This confirms the Voyager V∞ is not
  in the public literature in clean form — the SPK-direct path is the ONLY
  reliable route.**
* Verdict: NEGATIVE for V0 admission, but strongly motivates #390.

## Rinker, Jacobson & Wood 1976 — off-scope for #345

`cyclers_pdf` commit `aca2a36`. "Statistical Analysis of Trim Maneuvers in
Low-Thrust Interplanetary Navigation," JSR (1976). Solar-electric-propulsion
navigation trim-maneuver methodology for comet/asteroid missions. A
low-thrust-navigation reference for #309 / #359, NOT a #345 mission-tour
admission paper. No catalogue impact.

## Net #345 status after this digest

Voyager (×2 papers), Mariner-10 (×2 papers) all join the #345
confirmed-negative set: NEITHER the mission-overview NOR the dedicated
navigation literature publishes per-encounter V∞ in a V0-admissible
tuple. The catalogue cannot admit these mga_tour rows until:

**The NAIF SPK direct-evaluation path (#390) is built** — derive V∞ at
each encounter from the mission's archived SPK kernels using the flyby
epochs we already have. This is now the ONLY remaining route (the
acquisition route is exhausted: even the best nav paper lacks V∞, and the
true sources are unobtainable internal JPL docs). It is a methodology
build, not an acquisition; it unblocks Voyager 1/2, Mariner-10, Galileo,
Cassini, Pioneer, Juno all at once and is the highest-leverage way to
clear the #345 backlog wholesale. The #361/#384/#387 acquisition tasks
for V∞ are effectively SUPERSEDED by #390.

The two papers ARE valuable as catalogue corroborating-source / KNOWN_CORPUS
context for the missions; only the V0 numeric-admission bar is unmet.
