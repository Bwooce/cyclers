# Digest — D'Amario & Byrnes 1983 "Interplanetary Trajectory Design for the Galileo Mission" (AIAA-83-0099)

**Date**: 2026-06-19 AET (Agent dispatched by parent #384; corpus-document-policy digest)
**Verdict (TL;DR)**: **NOT a cycler; NOT catalogue-admissible as a cycler/mga_tour row.** This is the **interplanetary-leg companion** to the already-digested Diehl-Kaplan-Penzo AIAA-83-0101 satellite-tour paper (`2026-06-17-digest-diehl-1983.md`). It documents the **pre-Challenger 1986-launch DIRECT Earth→Jupiter trajectory** (broken-plane, single deep-space maneuver), which was NEVER FLOWN — the flown mission is the 1989 VEEGA (`2026-06-17-digest-damario-1992-galileo.md`, catalogue row `damario-1992-galileo-veega`). It publishes Jupiter-arrival V_inf (VHP) but it is a direct E→J two-leg trajectory with **no planetary gravity-assist encounters** (only a deep-space broken-plane maneuver), so there are no per-encounter heliocentric V_inf magnitudes of the kind the mga_tour rows require. Value to the corpus is **provenance / errata corroboration**, not new catalogue data.

---

## 1. Header

- **Title** (verbatim, cover & p.1): "Interplanetary Trajectory Design for the Galileo Mission"
- **Paper number**: AIAA-83-0099 (cover; running header A83-17908)
- **Authors** (verbatim, p.1 byline + footnotes): **Louis A. D'Amario** (Member Technical Staff, Member AIAA), **Dennis V. Byrnes** (Consultant, Member AIAA)
- **Affiliation**: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California 91109
- **Venue**: AIAA 21st Aerospace Sciences Meeting, January 10-13, 1983, Reno, Nevada (cover page)
- **Copyright**: American Institute of Aeronautics and Astronautics, Inc., 1983 (p.1 footer)
- **Source**: AIAA 83-0099
- **Page count**: 10 pages (cover + 9 body), text layer present

## 2. What the paper actually is

A JPL mission-design analysis of the **1986 direct Earth→Jupiter** Galileo trajectory, after the mission was "reprogrammed to use a direct Earth-Jupiter trajectory with a May 1986 launch date and with arrival at Jupiter occurring in mid-1988" (Abstract, p.1). The reprogramming was driven by the wide-body Centaur-F upper stage reinstatement in July 1982, replacing the previous two-stage IUS baseline that needed a 1985 ΔVEGA (Earth-gravity-assist) trajectory launched May 1985, arriving late 1989 (Introduction, p.1; Figs 1-2, p.1).

**Three trajectory concepts appear in the paper:**
1. **1985 ΔVEGA** (Fig 1, p.1): prior baseline — Earth launch May 1985, Earth flyby July 1987, Jupiter arrival Nov 1989. (Superseded; this is the EGA concept, NOT the flown VEEGA.)
2. **1986 direct broken-plane** (Fig 2, p.1): the new baseline — Earth launch May 1986, Jupiter arrival Aug 1988, single deep-space broken-plane maneuver Jan 1987. **This is the paper's subject.**
3. The flown trajectory (1989 VEEGA) does NOT appear here — it postdates this paper (Challenger forced the IUS substitution and the 1989 VEEGA; documented in D'Amario-Bright-Wolf 1992).

The "broken-plane" concept (Figs 2, 4): a near-180° Earth→Jupiter Type I transfer with a single mid-course maneuver (~8 months after launch, over 2/3 of transfer angle) that supplies the heliocentric inclination needed to reach out-of-ecliptic Jupiter, rather than paying it all at launch via high C3 (Section II.B, p.2).

## 3. Extracted data (sourced, with page citations)

### Launch energy / departure (1986 direct trajectory)
- **Best Type I ballistic C3 ≈ 84 km²/s²**; best Type II ballistic C3 ≈ 81 km²/s² (p.2). Centaur capability = **80 km²/s²** for the injected spacecraft mass of 2550 kg (p.2, p.9). ⇒ ballistic transfers INFEASIBLE in 1986 without a broken-plane maneuver.
- Minimum C3 for near-Hohmann ecliptic transfer would be **~78 km²/s²** with ~2-yr transfer (p.2) — but Jupiter is out of ecliptic, raising the requirement.

### Reference broken-plane trajectory — Table 1 (p.3), launch 5/21/86, arrival 8/27/88
| Parameter | Ballistic | Broken-plane |
|---|---|---|
| C3 | **118 km²/s²** | **80 km²/s²** |
| DLA (Earth equator) | -48° | -22° |
| DLA (ecliptic) | -42° | -12° |
| Pre-maneuver inclination (ecliptic) | 11° | 2.9° |
| Post-maneuver inclination (ecliptic) | — | 2.4° |
| Broken-plane maneuver | — | **231 m/s** |
| **VHP at Jupiter** | **6.1 km/s** | **5.8 km/s** |

VHP = hyperbolic approach velocity at Jupiter = **Jupiter-arrival V_inf**. This is the only per-encounter V_inf the paper publishes (and Jupiter is the only encounter; the broken-plane maneuver is a deep-space ΔV, not a flyby).

### Other sourced values
- Broken-plane maneuver magnitude across the 10-day launch period (May 21-31): **150-230 m/s** for the Aug 27 arrival (p.7, Fig 9/10).
- Jupiter VHP across launch/arrival space: ~5.7-6.1 km/s, lowest for later arrival dates (Fig 6, p.4).
- Example trajectory (launch 5/21/86, arrival 8/27/88, C3 118): heliocentric inclination 11°, DLA -48° (equator) / -42° (ecliptic), VHP 6.1 km/s (p.2, ballistic).
- Optimal broken-plane: leaves Earth with ecliptic inclination 2.9°, DLA -22° (equator)/-12° (ecliptic); maneuver ~1/3 of flight time after launch (over 2/3 of transfer angle); changes ecliptic inclination by 0.5° and node by 18° (p.2-3).
- Spacecraft mass breakdown (Table 4, p.5): Orbiter dry 1138 kg, Probe 335 kg, RPM propellant 932 kg, adapter 145 kg, **total injected 2550 kg**.
- Propellant-margin example (Table 5, p.6, launch 5/25/86 arrival 8/27/88): JOI ΔV **668 m/s**, BPM 167 m/s, PJR 351 m/s, Tour ΔV 205 m/s, final orbiter mass 1171 kg, margin 33 kg.
- Post-Io perijove radius 4 R_J; JOI 90 min after probe entry; initial orbit period 200 days; PJR establishes first in-orbit perijove at 11.5 R_J (p.4-5).

### Targets of opportunity (Section IV, Table 7, p.8-9)
Preliminary asteroid/comet flyby search on the 1986 direct trajectory (~4000 objects screened):
| Name | Type | Flyby date | Jupiter arrival | Propellant penalty |
|---|---|---|---|---|
| 6647 P-L71 | Palomar-Leiden asteroid | Aug 27 1986 | Aug 27 1988 | 26 kg |
| Britta | Main belt asteroid | Oct 27 1986 | Oct 26 1988 | 21 kg |
| Yi Xing | Main belt asteroid | Nov 19 1986 | Oct 26 1988 | 34 kg |
| 4698 P-L71 | Palomar-Leiden asteroid | Nov 21 1986 | Sep 16 1988 | 15 kg |
| Bus | Comet | Feb 24 1988 | Oct 23 1988 | 37 kg |

No V_inf published for these candidate flybys — only propellant-margin penalties (p.8-9). Preliminary/non-exhaustive (p.8).

## 4. Per-encounter V_inf? — only Jupiter VHP
**The only encounter on this trajectory is Jupiter arrival.** VHP (V_inf) at Jupiter = **5.8 km/s** (broken-plane) / **6.1 km/s** (ballistic), Table 1 p.3. There are **no planetary gravity-assist flybys** on the 1986 direct trajectory (the EGA was in the *superseded* 1985 ΔVEGA baseline, Fig 1, and the GAs only return with the post-Challenger 1989 VEEGA). The deep-space broken-plane maneuver is a propulsive ΔV, not a flyby V_inf. So the per-encounter V_inf data density that admitted D'Amario 1992 as an mga_tour row is **absent** here.

## 5. References (p.9) — corpus cross-links
1. Mitchell, R.T., "Project Galileo Mission Design," IAF-80-G-291, 1980.
2. Casani, J.R., "Galileo: Mission to Jupiter," IAF-81-204, 1981.
3. O'Neil, W.J., "The Galileo ΔVEGA Mission to Jupiter," IAF-82-192, 1982.
4. O'Neil, W.J., "Galileo Mission Overview," AIAA Paper 83-0096, Jan 1983.
5. **D'Amario, L.A. and Byrnes, D.V., "Interplanetary Trajectory Optimization with Application to Galileo," J. Guidance, Control and Dynamics, Vol. 5, No. 5, Sep-Oct 1982, p.465.** *(the optimizer this paper relies on — sibling to the 1981 MOSES paper digested separately)*
6. Gill, P.E. and Murray, W., Numerical Methods for Constrained Optimization, Academic Press, 1974, pp.29-66. *(the NPL Newton method also cited by the 1981 MOSES paper)*

The Diehl-Kaplan-Penzo AIAA-83-0101 satellite-tour paper cites THIS paper (AIAA-83-0099) as its interplanetary-leg companion (per `2026-06-17-digest-diehl-1983.md`).

## 6. Catalogue / KNOWN_CORPUS relevance

### 6.1 Catalogue admission — DO NOT ADMIT as a cycler/mga_tour row
- Not a cycler (single direct E→J transfer, not repeating).
- Not an mga_tour: no chain of gravity-assist flybys; the only encounter is Jupiter arrival; the broken-plane is a deep-space ΔV. The per-encounter V_inf table that admitted D'Amario 1992 (Venus/Earth-1/Earth-2 GAs) does not exist here.
- The trajectory was **never flown** (pre-Challenger 1986 direct concept). The catalogue already carries the FLOWN mission as `damario-1992-galileo-veega`. Admitting this concept paper would duplicate-by-design-epoch a superseded variant with thinner data.

### 6.2 KNOWN_CORPUS / errata corroboration (the real value)
- This is the **third independent Galileo-attribution corroborator** alongside Diehl-Kaplan-Penzo 1983 (AIAA-83-0101) and D'Amario-Bright-Wolf 1992. It confirms the **D'Amario & Byrnes** author pairing for the Galileo interplanetary leg, AIAA-83-0099, 1983 Reno.
- Reinforces the `2026-06-17-digest-damario-1992-galileo.md` §5 finding that the KNOWN_CORPUS "Diehl-Belbruno-Roberts Galileo VEEGA design (1986)" anchor is mis-attributed: the 1983 interplanetary paper is **D'Amario-Byrnes**, the 1983 tour paper is **Diehl-Kaplan-Penzo**, and the flown VEEGA is **D'Amario-Bright-Wolf 1992**. "Belbruno" appears in none of these three. Recommend the parent fold this corroboration into the KNOWN_CORPUS anchor cleanup already proposed in the D'Amario 1992 digest §5.
- Documents the design-epoch sequence cleanly: **1985 ΔVEGA (EGA concept) → 1986 direct broken-plane (this paper) → 1989 VEEGA (flown, D'Amario 1992).** Useful for any future "Galileo mission design lineage" note.

### 6.3 Optimizer-lane note
Section II.A (p.2) states the trajectories were generated with the **PLATO (PLanetary Trajectory Optimization) program**, a constrained-optimization scheme using the NPL Newton algorithm of Gill & Murray (ref 6) — the same algorithmic family as the D'Amario-Byrnes-Stanford 1981 MOSES optimizer (digested separately). Confirms the Newton + analytic-derivative + bounds-constrained lineage is the JPL house method behind both Galileo design and the multiple-flyby optimizer relevant to the project's #380/#347 corrector lane.
