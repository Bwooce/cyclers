# Bellerose-Roth-Wagner 2018 Cassini reconstruction deep-read digest (#382 paper 4)

Per #382 task brief 2026-06-17 AET (#361 Cassini admission unblock candidate). All 10 pages read; this digest is sourced strictly from the paper itself.

## 1. Header (verbatim)

- **Title:** "The Cassini Mission: Reconstructing Thirteen Years of the Most Complex Gravity-Assist Trajectory Flown to Date"
- **Authors:** Julie Bellerose (Navigation Engineer, Outer Planet Navigation Group, Mission Design and Navigation Section), Duane Roth (Delivery Manager, Mission Design and Navigation Section), Sean Wagner (Navigation Engineer, Flight Path Control Group, Mission Design and Navigation Section). All Jet Propulsion Laboratory / Caltech, Pasadena, US, 91109.
- **Venue:** SpaceOps Conferences, 15th International Conference on Space Operations, 28 May - 1 June 2018, Marseille, France.
- **DOI:** `10.2514/6.2018-2646`. Available via arc.aiaa.org (the paper is labeled "Downloaded by UNIVERSITY OF MINNESOTA on May 26, 2018").
- **Copyright:** California Institute of Technology, published by AIAA.
- **Page count:** 10 (pp. 1-10).

## 2. What the paper actually is

A **post-mission-reconstruction methodology paper**, NOT a per-flyby trajectory-publishing paper. The paper describes the **uniform-reconstruction effort** undertaken by the Cassini Navigation team (2017-2018) to deliver a single, consistent set of trajectory reconstructions using one uniform Saturn-system model + one uniform satellite-ephemeris set, replacing the **~100+ different models** used across 172 navigation deliveries during the 13-year mission.

The four substantive sections (§I-IV):

1. **§I Introduction** (pp. 1-2) — Cassini mission timeline; prime + Equinox + Solstice mission counts; SOI 2004-07-01; impact 2017-09-15; mission-profile Figure 1 with year-by-year flyby counts.
2. **§II Uniform Reconstruction Structure** (pp. 3-7) — Navigation background (arc structure, OD + Maneuver teams), Saturn pole and Titan position uncertainty evolution over time (Fig 3), inputs preparation (Table 1: 5000 DSN passes, 80000 ionosphere entries, 10000 radiometric data edits, 2243 optical images, 2253 small-thruster events, 80 encounter thrusting events, 102 spacecraft attitude files), uniform-reconstruction tool design (`auto_recon.py`, Fig 4).
3. **§III Preliminary Results** (pp. 8-9) — Trajectory differences between uniform reconstruction (using sat389 model) and operational reconstructions: average position differences <7 km / velocity 150 mm/s through T16 (~2006), then <4.5 km / 70 mm/s for the rest of the tour.
4. **§IV Conclusion** (p. 10) — Cassini uniform reconstruction to be published on NAIF in Summer 2018.

Fifteen references (Refs 1-15, p. 10):
- Cassini Navigation Plan + Extended Navigation Plan + Solstice Navigation Plan (Refs 1-3, all JPL Technical Report D-11621 series — JPL-internal but commonly cited).
- **Cassini NAIF website https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/** (Ref. 4) — this is the actual trajectory delivery, NOT in this paper.
- Buffington 2010 IAC "Proposed End-of-Mission" (Ref. 5).
- Antreasian 2008 AIAA "Cassini OD Results January 2006 - End of Prime" (Ref. 6).
- Pelletier 2012 SpaceOps "Cassini OD Results July 2008 - December 2011" (Ref. 7).
- Bellerose 2016 AAS-16-142 "Cassini Navigation: Road to Consistent sub-Kilometer Accuracy" (Ref. 8).
- Vaquero 2014 AIAA-2014-4348 "Cassini Maneuver Experience for 4th Year of Solstice Mission" (Ref. 9).
- Vaquero 2018 AAS-18-151 "Flying Cassini Through the Grand Finale Orbits" (Ref. 10).
- Wagner 2016 AAS-16-305 "Cassini Maneuver Performance" (Ref. 11).
- Cassini MONTE + ODP tools (Refs 12-13).
- Criddle 2017 AAS-17-625 "Optical Navigation during Cassini's Solstice Mission" (Ref. 14).
- Bellerose 2018 AAS-18-152 "Cassini OD Operations through Final Titan Flybys and Mission Grand Finale" (Ref. 15).

## 3. Per-encounter data extracted

### 3.1 Body sequence and tour-level facts (sourced)

| Item | Value | Source |
|---|---|---|
| Launch year | 1997 | abstract p. 1 |
| SOI date+time | **2004-07-01 14:00 ET** | p. 4 |
| Equinox mission start | 2008-09 | p. 1 |
| Saturn equinox | 2009-08-11 | p. 1 |
| Solstice mission span | 2010-2017 | p. 1 |
| Cassini equatorial-phase return | 2012-06 | p. 1 |
| Final Titan T125 (grazed F-ring outer edge) | (date not given) | p. 1 |
| Final Titan T126 (Grand Finale gate) | (date not given) | p. 1 |
| Grand Finale F-ring/D-ring period | 6.4 days × 22 times | p. 1 |
| **Saturn atmospheric impact (EOM)** | **2017-09-15 11:54 ET** | p. 4 |
| Total Cassini orbits at Saturn | "almost 300" | abstract p. 1 |
| Total maneuvers executed | **360 maneuvers** out of 492 designed | p. 5 |
| Total tracking passes (DSN) | ~5000 | Table 1 p. 5 |
| Total nav arcs | **172 delivered** (157 in uniform reconstruction) | p. 4 |
| Total optical images for navigation | 2243 post-SOI + 609 pre-SOI not processed = 2852 total | Table 1 p. 5 |
| Total small-thruster events | 2253 (in 80 encounters where thrusting was used) | Table 1 p. 5 |

### 3.2 Per-mission-year flyby counts (Fig 1 p. 2) — sourced

| Mission Phase | Year of Tour | Years | Orbits | Titan flybys | Enceladus | Other icy <10,000 km |
|---|---|---|---|---|---|---|
| Prime | 1 | 2004-2005 | 11 | (count from Fig 1 dots) | 0 | Phoebe |
| Prime | 2 | 2005-2006 | 15 | (dots) | (dots) | Hyperion, Dione, Tethys, Rhea, Telesto |
| Prime | 3 | 2006-2007 | 22 | (dots) | 0 | Rhea |
| Prime | 4 | 2007-2008 | 27 | (dots) | 0 | Iapetus, Epimetheus |
| Equinox | 5 | 2008-2009 | 39 | (dots) | (many dots) | Mimas, Rhea, Helene, Dione, G ring arc |
| Equinox | 6 | 2009-2010 | 21 | (dots) | (dots) | Rhea, Helene, Dione |
| Solstice | 7 | 2010-2011 | 16 | (dots) | (dots) | (none under 10,000 km) |
| Solstice | 8 | 2011-2012 | 19 | (dots) | (dots) | Dione, Tethys, Methone, Telesto |
| Solstice | 9 | 2012-2013 | 25 | (dots) | (dots) | Rhea |
| Solstice | 10 | 2013-2014 | 12 | (dots) | 0 | Dione, Tethys |
| Solstice | 11 | 2014-2015 | 12 | (dots) | 0 | Dione, Tethys |
| Solstice | 12 | 2015-2016 | 20 | (dots) | (dots) | Dione, Epimetheus, G arc |
| Solstice | 13 | 2016-2017 | **56** (Grand Finale) | (dots, none in Grand Finale) | 0 | (none) |

**Mission-totals per §I:**
- **Prime (2004-2008):** 45 Titan + 4 Enceladus + 9 other icy (Huygens on Titan T1)
- **Equinox (2008-2010):** 26 Titan + 12 Enceladus + 12 other icy
- **Solstice (2010-2017):** 46 Titan + 12 Enceladus + 12 other icy
- **Tour total:** **117 Titan + 28 Enceladus + 33 other icy = 178 satellite flybys**

(Compare to Wolf-Smith 1995 sample tour: 33 Titan + 1 Enceladus + 5 other = 38 — the actual flown mission is 4-5x larger than the 1995 sample tour.)

### 3.3 V_∞ per encounter — **NOT PUBLISHED**

**Critical finding:** the paper does NOT publish V_∞ at any flyby. Fig 1 shows a per-mission-year science profile with **dot-counted flyby tallies** but no numerical per-flyby table at all. The paper's substance is the **uniform-reconstruction process** (auto_recon.py tool, Fig 4 software diagram), NOT a per-flyby trajectory table.

The published numerical trajectory data is:
- **Position-difference statistics** between uniform reconstruction (sat389) and operational reconstructions: <7 km / 150 mm/s through T16 (~2006), <4.5 km / 70 mm/s thereafter (§III p. 8).
- **Saturn pole right ascension over time:** 40.57 → 40.586 deg (Fig 3 p. 4).
- **Titan position uncertainty over time:** 2.7 km in 2005 → <5 km by T_c (end 2004) → ~0 by 2007 (Fig 3 p. 4 — Titan uncertainty fell from 200 km pre-mission to <5 km after first Titan encounters).
- **Estimated/consider parameter uncertainties** (Table 2 p. 7): Cassini state <5 km, <20 cm/s; ΔV small burns 0.25-1.2 mm/s; OTM(ME) 0.02% proportional + 3.5 mm/s fixed; OTM(RCS) 0.4% proportional + 0.5 mm/s fixed; Saturn 0.2 km; Saturn J2-J14 <0.01%; Saturn pole (RA, DEC) <0.01%, <0.0001%; CD Saturn 100%; satellites <0.1 km Titan, <km icys (Hyperion/Iapetus/Phoebe 2-5x higher).

### 3.4 The actual flyby tables ARE in NAIF SPK kernels, NOT in this paper

The paper explicitly directs readers to **NAIF kernels for the actual trajectory data:**

> "Cassini's trajectory reconstructions are already publicly available on the NAIF website [4]; deliveries were usually made every few months." (p. 2)
> "Now that the mission is done, ... the paper reports on the uniform reconstruction of the entire Cassini orbital mission" (abstract p. 1)

Ref. 4: `https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/`

**The V_∞ data is recoverable from the NAIF SPK kernels by evaluating Cassini's Saturn-relative state at each flyby epoch and computing satellite-relative V_∞ from the satellite SPK at the same epoch.** This is **derived, not sourced** in the strict cycler-catalogue sense — though SPICE-kernel-derived values are higher fidelity than any tabulated V_∞ would be.

## 4. Catalogue admission verdict

### 4.1 V0 admissibility test

Per the Heaton-Longuski 2003 U00-01 admission precedent and the Wolf-Smith 1995 Cassini negative-verdict precedent (`docs/notes/2026-06-17-digest-wolf-smith-1995-cassini.md`) and the three prior negatives in this digest sweep (Bourke 1971, Dunne-Burgess 1978, Lam 2008), the V0 evidence standard for an `mga_tour` row is:

**Per-flyby epoch + V_∞ tuple + body sequence are minimum.**

Bellerose-Roth-Wagner 2018 publishes:

- Per-flyby epoch: **NO** at the paper level. The paper publishes only **tour-level epochs** (SOI 2004-07-01 14:00 ET, EOM 2017-09-15 11:54 ET) and dot-counted flyby-count tallies per mission-year. Per-flyby epoch is in the NAIF SPK kernels (Ref. 4), NOT in this paper.
- Body sequence: PARTIAL. The bodies are visible in Fig 1 (Titan, Enceladus, Phoebe, Hyperion, Dione, Tethys, Rhea, Telesto, Iapetus, Epimetheus, Mimas, Helene, G ring arc, Methone) — but the **temporal sequence is not encoded numerically** in the paper.
- V_∞ tuple: **NO**. Not published anywhere in this paper.

**Verdict: the paper is TWO+ COLUMNS SHORT of the V0 `mga_tour` minimum** (worse than Wolf-Smith 1995 Cassini, which at least published Tour-18-5 per-flyby epochs to the second).

### 4.2 #361 unblock test

The brief explicitly asks: **does Bellerose 2018 publish per-flyby V_∞ tables (which would close #361)?**

**Answer: NO.** Bellerose 2018 is a **methodology paper about the uniform-reconstruction software tool**, not a tabulated per-flyby publication. The actual trajectory data (V_∞, B-plane, closest-approach altitude per flyby) is on the **NAIF SPK kernels** referenced as Ref. 4, NOT in any tabulated form in this paper.

**#361 stays acquisition-gated. Bellerose 2018 does NOT unblock it.** The acquisitions needed to actually unblock #361 are listed in §7 below — at minimum:
- Pelletier 2012 SpaceOps "Cassini OD Results July 2008 - December 2011" (Ref. 7) — per-flyby OD residuals from the Equinox+early-Solstice missions.
- Bellerose 2016 AAS-16-142 (Ref. 8) — covers the Solstice mission consistency.
- Bellerose 2018 AAS-18-152 (Ref. 15) — covers the final Titan flybys + Grand Finale.
- The actual NAIF SPK kernels themselves — but those are derived V_∞ post-eval.

The Strange-Russell-Buffington 2007 / Yam-Davis-Longuski 2009 / Valerino 2014 trio currently in our literature anchor are still post-launch sources that **may or may not** publish per-flyby V_∞ tables — they are designated for the Cassini end-of-mission ΔV analysis and Saturn-impact analysis, not for tabulated flyby V_∞ for the prime+Equinox+Solstice main tour.

### 4.3 Recommendation

**RECOMMEND: do NOT admit Cassini as a `mga_tour` row from Bellerose 2018 alone.** V_∞ data is not published; per-flyby epochs are not tabulated; the body sequence is only dot-counted in Figure 1.

Instead:

1. **File Bellerose 2018 as a KNOWN_CORPUS reconstruction-methodology anchor** for the Cassini mission. It is the canonical reference for the **uniform-reconstruction effort** and the **177-flyby tour-total statistics** (45 Titan + 4 Enceladus + 9 other in Prime + 26 + 12 + 12 in Equinox + 46 + 12 + 12 in Solstice), and a citable source for tour-level epochs (SOI 2004-07-01 14:00 ET; EOM 2017-09-15 11:54 ET).

2. **(Action — acquisition gate, HIGH PRIORITY for #361)** **Pelletier-Antreasian-Ardalan 2012 SpaceOps "Cassini OD Results (July 2008 - December 2011)"** (Ref. 7) — this is a series of OD-results papers covering the prime + Equinox + Solstice missions:
   - Antreasian 2008 AIAA-2008-6747 "Cassini OD Results: January 2006 - End of Prime" (Ref. 6) — Prime mission OD.
   - Pelletier 2012 SpaceOps "Cassini OD July 2008 - December 2011" (Ref. 7) — Equinox/Solstice.
   - Bellerose 2018 AAS-18-152 (Ref. 15) — Final Titan + Grand Finale.
   These are the per-flyby OD-residuals papers that may publish V_∞ tables.

3. **(Action — alternative path)** **Direct evaluation from NAIF SPK kernels.** The Cassini SPK is publicly available at `https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/` and contains all the Cassini state vectors throughout the mission. We can evaluate Cassini's Saturn-relative state at any flyby epoch and compute V_∞ at each satellite. **This is derived-not-sourced for V_∞** — but it is the highest-fidelity source available (no truncation, no interpolation, the actual flown trajectory). For #361 / V0 admission, this would violate `golden-tests sourced-only` — but the **derived V_∞ values can be USED in catalogue cells with `vinf_source: derived` + `vinf_fidelity: real-ephemeris` per the schema v4 fidelity columns**, and the row would be admitted at V0 with the body sequence + epoch tuple from a sourced paper (Bellerose 2018 + Wolf-Smith 1995 for the 1995 sample tour body set).

4. **If a V_∞ table surfaces** from any of Refs 6-7-15, admit a `cassini-flown-tour-2004-2017` row at V0 with the structural template below.

### 4.4 Structural template (for future use, NOT for writeback now)

If/when V_∞ data surfaces (most likely from Antreasian 2008 or Pelletier 2012, or from NAIF-SPK-derived values):

```yaml
- id: cassini-flown-tour-2004-2017
  name: "Cassini Saturn satellite tour, flown (Bellerose-Roth-Wagner 2018 reconstruction; 178 flybys / 360 maneuvers / 13 years)"
  source: literature
  trajectory_regime: ballistic   # 178 ballistic-conic flyby arcs interrupted by 360 OTMs; impulsive ΔV
  model_assumption: real-ephemeris   # NAIF SPK / sat389 uniform reconstruction per Bellerose 2018
  cycler_class: multi-arc
  orbit_class: mga_tour
  epoch_locked: true
  n_returns: 1
  validity_window:
    start: "2004-07-01T14:00:00Z"   # SOI per Bellerose 2018 p. 4
    end:   "2017-09-15T11:54:00Z"   # Saturn atmospheric impact per Bellerose 2018 p. 4
  launch_epoch: "1997-10-15T00:00:00Z"   # Cassini-Huygens launch (per NASA fact sheets; not in Bellerose 2018)
  validation_level: V0
  source_ephemeris: "sat389 Saturn-system model (Bellerose 2018) + DE-series solar-system ephemeris"
  bodies: ["Saturn", "Titan", "Enceladus", "Rhea", "Dione", "Iapetus", "Tethys", "Mimas", "Phoebe", "Hyperion", "Telesto", "Helene", "Epimetheus", "Methone", "Pallene"]
  sequence_canonical: <from NAIF SPK trajectory evaluation>
  vinf_kms_at_encounters: <NEEDS-MORE-DATA: not published in Bellerose 2018; derivable from NAIF SPK>
  flyby_totals:
    titan: 117   # 45 Prime + 26 Equinox + 46 Solstice
    enceladus: 28   # 4 + 12 + 12
    other_icy: 33   # 9 + 12 + 12
    grand_total: 178
  maneuvers_executed: 360
  maneuvers_designed: 492
  orbits_at_saturn: 294   # "almost 300" per abstract; precise count not given
  navigation_arcs: 172   # 157 used in uniform reconstruction
  first_published:
    authors: ["Wolf, A. A.", "Smith, J. C."]
    year: 1995
    title: "Design of the Cassini Tour Trajectory in the Saturnian System"
    venue: "Control Engineering Practice, Vol. 3, No. 11, pp. 1611-1619"
    doi: "10.1016/0967-0661(95)00172-7"
    notes: "Pre-launch sample tour (Tour 18-5); the flown tour is 4-5x larger."
  corroborating_sources:
    - authors: ["Bellerose, J.", "Roth, D.", "Wagner, S."]
      year: 2018
      title: "The Cassini Mission: Reconstructing Thirteen Years of the Most Complex Gravity-Assist Trajectory Flown to Date"
      venue: "SpaceOps 2018, Marseille, France"
      doi: "10.2514/6.2018-2646"
    - authors: ["Antreasian, P. G.", "Ardalan, S. M.", "Bordi, J. J."]  # et al.
      year: 2008
      title: "Cassini Orbit Determination Results: January 2006 - End of Prime Mission"
      venue: "AIAA/AAS Astrodynamics Specialist Conference"
      paper_id: "AIAA 2008-6747"
    - authors: ["Pelletier, F. J.", "Antreasian, P.", "Ardalan, S. M."]  # et al.
      year: 2012
      title: "Cassini Orbit Determination Results (July 2008 - December 2011)"
      venue: "SpaceOps 2012, Vol 1, p. 167"
    - authors: ["Strange, N. J.", "Russell, R.", "Buffington, B. B."]
      year: 2007
      title: "Mapping the V-infinity globe"
      venue: "AAS 07-277"
    - authors: ["Yam, C. H.", "Davis, D. C.", "Longuski, J. M.", "Howell, K. C.", "Buffington, B."]
      year: 2009
      title: "Saturn Impact Trajectories for Cassini End-of-Mission"
      venue: "Journal of Spacecraft and Rockets"
      doi: "10.2514/1.38760"
    - authors: ["Valerino, P."]
      year: 2014
      title: "Updating the Reference Trajectory for the Cassini Solstice Mission"
      venue: "SpaceOps 2014"
      doi: "10.2514/6.2014-1880"
```

This template combines Wolf-Smith 1995 (pre-launch design, `first_published`) with Bellerose 2018 + the OD-results trio (post-flight reconstruction, `corroborating_sources`). Writeback only when at least one of the V_∞-tabulated sources is in hand.

## 5. KNOWN_CORPUS impact

The existing Cassini anchor (per `docs/notes/2026-06-17-digest-wolf-smith-1995-cassini.md` §5 / `docs/notes/2026-06-17-349-cassini-anchor-topology-label.md`) at `src/cyclerfinder/search/literature_check.py` is:

- Name: `"Cassini-Huygens Saturn-Titan satellite tour design"`
- body_set: `{Titan, Enceladus, Rhea, Dione, Iapetus}` (per the most-recent version per Wolf-Smith 1995 update)
- topology_label: `{pump-tour, mga-tour}`
- Citation: includes Strange / Yam / Valerino / Wolf+Smith.

**RECOMMENDED amendments from Bellerose 2018:**

1. **Add Bellerose 2018 as the AUTHORITATIVE post-flight reconstruction reference.** Suggested rewording of the citation field to add: `"Bellerose, Roth & Wagner, 'The Cassini Mission: Reconstructing Thirteen Years of the Most Complex Gravity-Assist Trajectory Flown to Date,' SpaceOps 2018 DOI 10.2514/6.2018-2646 -- post-flight uniform reconstruction: 178 flybys (117 Titan + 28 Enceladus + 33 other icy), 360 maneuvers, ~300 orbits, SOI 2004-07-01 14:00 ET, atmospheric impact 2017-09-15 11:54 ET; published on NAIF SPK kernels."`

2. **Extend body_set to add all the bodies named in Bellerose 2018 Fig 1:** Phoebe, Hyperion, Telesto, Epimetheus, Mimas, Helene, Methone, Pallene. Suggested: `frozenset({"Titan", "Enceladus", "Rhea", "Dione", "Iapetus", "Tethys", "Mimas", "Phoebe", "Hyperion", "Telesto", "Helene", "Epimetheus", "Methone"})`. (Pallene is not in Fig 1's "Other Icy Satellites (under 10,000 km)" list and may not have had a sub-10000-km flyby; conservative.)

3. **Add Bellerose 2018 + Pelletier 2012 + Antreasian 2008 to authors tuple:** `("Wolf", "Smith", "Bellerose", "Roth", "Wagner", "Antreasian", "Pelletier", "Strange", "Russell", "Buffington", "Yam", "Davis", "Longuski", "Valerino")`.

These are RECOMMENDATIONS — not edits in this session.

## 6. Errata

Versus the task brief and the paper text:

- **Brief said "this is a POST-MISSION reconstruction paper. Wolf-Smith 1995 published the PRE-LAUNCH base tour but no V_∞ table (#355). Bellerose 2018 documents the actual flown tour across ~300 orbits + 13 years."** Confirmed. Bellerose 2018 IS a post-mission reconstruction paper.
- **Brief said "CRITICAL: does it publish per-flyby V_∞ tables (which would close #361)? If yes, #361 unblocks and Cassini becomes the 4th catalogue mga_tour row. If no (e.g. paper is primarily about uniform-reconstruction methodology), REPORT honestly and #361 stays acquisition-gated."** **Confirmed: NO V_∞ tables. The paper is about uniform-reconstruction methodology. #361 stays acquisition-gated.**
- **Brief said "1997 launch."** Paper confirms abstract p. 1 "Cassini launched in 1997"; launch date is widely known as 1997-10-15 (NASA fact sheets) but is NOT in this paper.
- **Mission counts in §I:**
  - "The prime mission included 45 flybys of Titan, 4 of Enceladus, and 9 of other icy satellites." (p. 1)
  - "The Equinox mission added 26 Titan flybys, and 12 more Enceladus and icy satellite flybys." (p. 1)
  - "The Solstice mission, from 2010 to 2017, added 46 Titan flybys, 12 Enceladus, and 12 other icy satellite flybys." (p. 1)
  - **Implied total:** 45+26+46 = 117 Titan; 4+12+12 = 28 Enceladus; 9+12+12 = 33 other icy = **178 satellite flybys total**. This is consistent with widely-cited Cassini mission statistics.
- **Saturn pole right ascension precision:** Fig 3 shows 40.57-40.586 deg with sub-mdeg precision over 11 years — useful as a Saturn-system gravitational anchor.
- **Wolf-Smith 1995 sample tour scale:** Wolf-Smith 1995 had 38 targeted satellite flybys; Bellerose 2018 confirms the flown tour had 178 = **4.7x larger** than the 1995 sample. This confirms that admitting Wolf-Smith 1995 as a catalogue row would NOT actually represent the flown Cassini trajectory (it represents only a fragment).
- **Cross-paper:** the brief mentions "#346 Davis pattern" as the structural template — but the actual template I have followed is the Wolf-Smith 1995 Cassini negative-verdict template (`docs/notes/2026-06-17-digest-wolf-smith-1995-cassini.md`) plus the Heaton-Longuski 2003 admission YAML structure from `data/catalogue.yaml`. Both match Davis et al.

## 7. Action items

For the parent / #382 / #361 / #345:

1. **DO NOT writeback a `mga_tour` row for Cassini from Bellerose 2018 alone.** No V_∞ tables; no per-flyby epoch tuples; the body sequence is only dot-counted in Figure 1.

2. **(Action — #361 STAYS ACQUISITION-GATED)** The acquisition needed to actually unblock #361 is at minimum one of:
   - **Antreasian-Ardalan-Bordi 2008 AIAA-2008-6747 "Cassini OD Results: January 2006 - End of Prime Mission"** (Bellerose 2018 Ref. 6).
   - **Pelletier-Antreasian-Ardalan 2012 SpaceOps "Cassini OD Results (July 2008 - December 2011)"** (Bellerose 2018 Ref. 7).
   - **Bellerose 2018 AAS-18-152 "Cassini OD Operations through Final Titan Flybys and Mission Grand Finale (February 2016 - September 2017)"** (Bellerose 2018 Ref. 15).
   - **Buffington 2010 IAC-10 "Proposed End-of-Mission for the Cassini Spacecraft: Inner D Ring Ballistic Saturn Impact"** (Bellerose 2018 Ref. 5).
   These are highly specific JPL/AAS papers; AIAA ARC and AAS digital library should have them.

3. **(Action — alternative path)** **Open a task to evaluate per-flyby V_∞ from the NAIF SPK kernels** (`https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/`). This would produce derived V_∞ values at higher fidelity than any tabulated paper source, with `vinf_source: derived` + `vinf_fidelity: real-ephemeris`. Combined with the sourced epoch tuples (from any of the above OD-results papers), this could ratify a Cassini catalogue row at V0 with derived V_∞ acceptable under the schema's `vinf_fidelity` column. **Requires #361 prerequisites to specify whether derived-from-SPK V_∞ is acceptable for V0** — if not, the schema/spec needs amendment.

4. **(Action — KNOWN_CORPUS amendment)** Open a small follow-up task to amend the existing Cassini literature anchor (`src/cyclerfinder/search/literature_check.py`) to:
   - Add Bellerose 2018 as a CITATION on equal footing with Strange / Yam / Valerino.
   - Extend body_set to add Phoebe, Hyperion, Telesto, Epimetheus, Mimas, Helene, Methone (and confirm with Wolf-Smith 1995 + Bellerose 2018 cross-check).
   - Extend authors tuple to add Bellerose, Roth, Wagner.
   - Keep topology_label `{pump-tour, mga-tour}` (Wolf-Smith 1995 + Bellerose 2018 confirm).
   - Update the body_set frozenset and authors tuple to include the full Bellerose 2018 list.

5. **(Action — Cassini tour-total fact)** Bellerose 2018's 178 satellite flybys / 360 maneuvers / ~300 orbits / 13 years are useful as a reference benchmark for the catalogue's expected "scale of a flown mga_tour". For comparison: Voyager 2 (JSUN, ~9 satellite flybys), Galileo (~30 flybys), Cassini (~178 flybys) — Cassini is the largest by an order of magnitude.

6. **(Action — Wolf-Smith 1995 sample-tour scale documented)** The 4.7x ratio between Wolf-Smith 1995 sample (38 flybys) and Bellerose 2018 flown (178 flybys) is a real difference of substance. If Wolf-Smith 1995 is ever admitted as a catalogue row at V0, the row name should clearly indicate "pre-launch sample tour (Tour 18-5)" not "Cassini tour" to avoid conflation with the flown trajectory.

End of digest paper 4 of 4.
