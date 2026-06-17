# Wolf-Smith 1995 Cassini Tour deep-read digest (Agent E)

Per Agent E task brief 2026-06-17 AET (#349 / classic-mission `mga_tour` admissions). All 9 pages read; this digest is sourced strictly from the paper itself.

## 1. Header (verbatim)

- **Title:** "Design of the Cassini Tour Trajectory in the Saturnian System"
- **Authors:** A.A. Wolf and J.C. Smith, Members of Technical Staff, Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA, USA.
- **Venue:** *Control Engineering Practice*, Vol. 3, No. 11, pp. 1611-1619 (1995).
- **Publisher:** Pergamon / Elsevier Science Ltd.
- **Article ID:** 0967-0661(95)00172-7. DOI form: `10.1016/0967-0661(95)00172-7`.
- **Submission:** Received February 1995; in final form September 1995.
- **Page count:** 9 (pp. 1611-1619).

## 2. What the paper actually is

A **pre-launch base-tour design study** for the Cassini orbiter's four-year nominal satellite tour at Saturn. Cassini did not launch until October 1997 and SOI did not occur until 2004-07-01; this paper is the 1995 pre-launch reference design. The paper is published in a Cassini "special section" of *Control Engineering Practice* (referenced explicitly in §1: "The scientific objectives of the tour are discussed in the introduction to this special section").

The paper has two distinguishable contributions:

1. **Methodology transfer from Galileo.** §2 lays out the tour-design concepts (gravity-assist vector geometry, orbit pumping vs cranking, line-of-apsides rotation rules, transfer-orbit inclination constraints at 180/360 deg transfer angles, occultation geometry). Galileo's tour-design methodology (Wolf-Byrnes 1993; D'Amario-Bright-Wolf 1992 — already in our acquisitions corpus) is named as the foundation; Cassini extended it for Saturn-specific constraints (Titan-only same-body pumping because non-Titan satellites are too light, atmospheric drag floor at Titan, ring-plane crossing geometry).

2. **A sample tour ("Tour 18-5" is not named in the text — see §5 below — but Table 2 fully publishes it).** §5 walks through the sample tour qualitatively then Table 2 gives the encounter list. §4 describes the three-stage design pipeline (interactive initial design → optimizer that absorbs third-body / oblateness perturbations and drives breakpoint ΔVs to zero → precision integrator producing a flyable nominal estimate).

The level of published detail is much higher than the Wolf-Byrnes 1993 Galileo equivalent abstract suggests. Table 2 is a **per-encounter table with full date+time, altitude, latitude, west-longitude, post-flyby inclination, and a narrative-purpose comment for 54 encounters**, which is more than I expected for a 9-page paper. There is no per-encounter V∞ column.

## 3. Body sequence and per-encounter data

### 3.1 Tour-level facts (sourced)

- **Saturn Orbit Insertion (SOI) date:** 2004-07-01 (§3 p. 1615, "The orbiter arrives at Saturn on 1 July, 2004"). The paper does not give the SOI hour/minute.
- **Tour end date:** 2008-07-01 (§5 p. 1619, "The tour ends on 1 July, 2008, 4 years after insertion into orbit about Saturn, during orbit 63"). Confirmed by Table 2 last entry Titan 54 at 080504 18:12:14, with subsequent occultations and a closing aimpoint targeting an extended-mission Titan flyby.
- **Arrival inclination to Saturn equator:** ~17 deg (§3 p. 1615).
- **Phoebe approach flyby:** Mentioned in §3 ("permits a flyby of Phoebe on approach to Saturn") but Phoebe does NOT appear in Table 2 — Phoebe is pre-SOI / pre-orbit-1.
- **Huygens probe release:** Described qualitatively in §1 (separation maneuver after Titan impact targeting, probe enters atmosphere, relays data through orbiter overhead). The paper says "The probe mission is described in detail in Section 4" — **but §4 of this paper is actually titled "Methods and Software"** and contains no Huygens release details. **This is a paper-internal cross-reference error: the Huygens detail is presumably in a sibling paper in the same Cassini special section** (the §1 description of "Section 4" is for the special section, not this paper). The Huygens release happens around the first Titan encounter (Titan 1, 2004-11-27 per Table 2); the paper does not publish a precise release date.
- **Orbit count:** 63 orbits over the 4-year nominal tour (§5 p. 1616).
- **Targeted-flyby count:** 38 (§5 p. 1616: "This tour contains 38 close satellite flybys and 63 orbits. Of these flybys, 33 are of Titan and 5 of other satellites. One targeted flyby each of Enceladus, Tethys, Dione, Rhea, and Iapetus occurs"). This sits inside the 35-40 range advertised in the abstract.
- **Maximum inclination reached:** 76 deg (§5 p. 1618, end-of-tour high-inclination sequence Titan 54).
- **Minimum Titan flyby altitude (design constraint):** 950 km (§3 p. 1615, sample-tour assumption). Several Table 2 encounters hit this floor exactly: Titan 35, 36, 37, 38, 43, 46.
- **Ring-plane crossing constraint:** ≥ 2.7 R_S (§3 p. 1615).
- **Minimum inter-targeted-flyby gap:** 16 days for a few orbits, 19-20 days for most (§3 p. 1615 — operational rather than dynamical).
- **Occultation totals:** 32 Earth-by-Saturn occultations (5 near-equatorial + 27 high-inclination) and 10 Earth-by-Titan occultations (§5 p. 1619).

### 3.2 ΔV budget

The paper publishes **only one ΔV number**: ~500 m/s for the entire tour (§1 p. 1611, "the total ΔV available from the orbiter's thrusters is about 500 m/s for the entire tour containing 35-40 encounters"). There is **no per-phase ΔV breakdown table**. The optimizer (§4) is described as driving most breakpoint maneuvers to zero by varying flyby times and aimpoints, with the deterministic ΔV resulting from this process; no number is given for the sample tour's actual deterministic ΔV total.

### 3.3 Body sequence (Table 2 p. 1617)

Sample tour, 54 encounters (38 targeted + 16 nontargeted, by my count), all dates `yymmdd hhmmss`:

| # | Encounter | Date/time (UTC) | Alt km | Lat deg | W.Lon deg | Post-flyby inc deg | Purpose |
|---|---|---|---|---|---|---|---|
| 1 | Titan 1 | 04-11-27 15:05:29 | 1500 | 61 | 105.1 | 10.5 | Reduce period, inclination |
| 2 | Titan 2 | 05-02-15 08:28:23 | 1250 | 61.9 | 87.9 | 2.1 | (same) |
| 3 | Titan 3 | 05-04-04 04:27:47 | 2397 | 17.1 | 73.2 | 0.5 | (same) |
| 4 | Rhea 4 | 05-05-03 05:54:08 | 999 | -73.1 | 270.5 | 0.56 | Rotate orbit ccw |
| 5 | Dione 4N | 05-05-03 12:07:12 | 13377 | -11 | 290.6 | 0.5 | Imaging (nontargeted) |
| 6 | Titan 5 | 05-06-02 18:28:38 | 4408 | 1.4 | 286.9 | 0.4 | Rotate orbit ccw |
| 7 | Titan 6 | 05-07-09 16:48:05 | 4161 | -0.3 | 71.6 | 4 | (same) |
| 8 | Dione 6N | 05-07-11 09:06:46 | 83464 | 0.6 | 342.8 | 0.4 | Imaging (nontargeted) |
| 9 | Dione 7 | 05-08-07 10:15:29 | 1005 | -12.8 | 287 | 0.4 | Rotate orbit ccw |
| 10 | Titan 8 | 05-09-07 07:53:50 | 4957 | -0.4 | 287.6 | 0.3 | (same) |
| 11 | Titan 9 | 05-10-14 04:05:52 | 4752 | 1.7 | 70.9 | 0.3 | (same) |
| 12 | Titan 11 | 05-12-12 19:17:17 | 1852 | 76.8 | 275.8 | 10.1 | Increase inc for Iapetus flyby |
| 13 | Titan 12 | 06-01-13 16:48:03 | 2132 | 49.7 | 275.7 | 15.5 | (same) |
| 14 | Iapetus 13 | 06-02-18 10:51:49 | 931 | -21.2 | 204.3 | 15.4 | Iapetus imaging |
| 15 | Titan 13 | 06-03-02 09:13:41 | 1511 | 11.6 | 109.8 | 20.2 | Occultations of Saturn, rings |
| 16 | Titan 16 | 06-05-21 02:45:17 | 1125 | -67.6 | 314.2 | 8.5 | Reduce inc, rotate cw |
| 17 | Titan 17 | 06-06-22 00:00:41 | 1740 | -52.4 | 106 | 0.4 | Rotate clockwise |
| 18 | Titan 18 | 06-07-11 18:53:33 | 2229 | 0.9 | 250.5 | 0.4 | (same) |
| 19 | Rhea 19N | 06-08-21 04:53:23 | 68841 | 1.3 | 297.4 | 0.4 | Imaging (nontargeted) |
| 20 | Titan 19 | 06-08-23 13:28:37 | 1958 | 0 | 107.9 | 0.4 | Rotate clockwise |
| 21 | Titan 20 | 06-09-12 08:16:55 | 11924 | 0.7 | 66.5 | 0.4 | (same) |
| 22 | Enceladus 20N | 06-09-14 03:43:48 | 5421 | 10.4 | 20.5 | 0.4 | Imaging (nontargeted) |
| 23 | Rhea 20N | 06-09-14 15:45:10 | 33970 | -0.9 | 93.5 | 0.4 | Imaging (nontargeted) |
| 24 | Tethys 21 | 06-10-04 16:18:43 | 648 | 81.6 | 257.4 | 0.4 | Rotate clockwise |
| 25 | Rhea 22N | 06-10-25 10:26:42 | 58934 | -0.9 | 281.3 | 0.4 | Imaging (nontargeted) |
| 26 | Titan 22 | 06-10-26 18:48:09 | 11308 | -0.1 | 293.2 | 0.4 | Rotate clockwise |
| 27 | Titan 23 | 06-11-15 13:30:03 | 2224 | 0 | 251.1 | 0.4 | (same) |
| 28 | Titan 24 | 06-12-28 07:31:12 | 1028 | -1.1 | 108.3 | 0.1 | Target to Enceladus |
| 29 | Enceladus 25 | 07-01-16 18:45:29 | 605 | -8.3 | 215.1 | 0.2 | Enceladus imaging |
| 30 | Enceladus 27N | 07-02-28 09:38:45 | 87264 | 0.2 | 267.3 | 0.2 | Imaging (nontargeted) |
| 31 | Dione 27N | 07-02-28 14:24:07 | 76634 | 1 | 276.8 | 0.2 | Imaging (nontargeted) |
| 32 | Titan 27 | 07-03-02 03:51:01 | 16784 | 4.4 | 292.2 | 0.4 | Rotate clockwise |
| 33 | Titan 28 | 07-03-21 22:01:39 | 2189 | -0.4 | 251.4 | 0.4 | (same) |
| 34 | Dione 29N | 07-05-01 16:40:26 | 24268 | -0.5 | 347.4 | 0.4 | Imaging (nontargeted) |
| 35 | Titan 29 | 07-05-03 15:46:59 | 2027 | -0.9 | 106.5 | 0.3 | Rotate clockwise |
| 36 | Titan 30 | 07-05-23 10:22:59 | 1000 | 85.3 | 328.7 | 14.3 | Occultations of Saturn, rings |
| 37 | Mimas 31N | 07-06-18 05:46:19 | 97383 | -28 | 208 | 14.4 | Imaging (nontargeted) |
| 38 | Titan 32 | 07-07-10 06:29:13 | 1050 | -59 | 241.3 | 2.3 | Reduce inc. |
| 39 | Titan 33 | 07-08-11 03:41:47 | 2098 | -11.4 | 68.7 | 0.3 | Position node for hi-inc seq |
| 40 | Enceladus 33N | 07-08-13 02:09:14 | 27267 | 2.3 | 251.6 | 0.3 | Imaging (nontargeted) |
| 41 | Enceladus 34N | 07-09-02 05:25:49 | 53700 | -0.6 | 96.6 | 0.3 | Imaging (nontargeted) |
| 42 | Titan 35 | 07-09-24 14:36:50 | 950 | -52.9 | 120.3 | 14.2 | Hi-inclination sequence |
| 43 | Titan 36 | 07-10-10 13:31:45 | 950 | -76 | 174.8 | 29.7 | (hi-inc; Saturn-ring occs + min-alt Titan + Titan occs) |
| 44 | Titan 37 | 07-10-26 12:07:07 | 950 | -65.5 | 176.3 | 40.9 | (hi-inc) |
| 45 | Dione 38N | 07-11-09 02:32:59 | 99682 | 0.8 | 305.9 | 40.9 | Imaging (nontargeted) |
| 46 | Titan 38 | 07-11-11 10:35:36 | 950 | -26 | 140.8 | 51.5 | (hi-inc) |
| 47 | Titan 43 | 08-01-14 04:48:01 | 950 | -16.1 | 145.6 | 58.7 | (hi-inc) |
| 48 | Dione 45N | 08-02-01 18:06:59 | 91316 | 5.4 | 36.3 | 58.7 | Imaging (nontargeted) |
| 49 | Tethys 46N | 08-02-12 09:48:43 | 35957 | -22.3 | 104 | 58.7 | Imaging (nontargeted) |
| 50 | Titan 46 | 08-02-15 01:49:24 | 950 | -17.5 | 157.9 | 64.7 | Hi-inclination sequence |
| 51 | Tethys 49N | 08-03-11 17:59:16 | 33630 | -12.8 | 178.4 | 64.7 | Imaging (nontargeted) |
| 52 | Titan 51 | 08-04-02 21:29:11 | 1022 | 29.6 | 133 | 70.3 | Hi-inclination sequence |
| 53 | Titan 53 | 08-04-18 19:49:51 | 950 | -19.9 | 177.2 | 72.1 | (hi-inc) |
| 54 | Titan 54 | 08-05-04 18:12:14 | 958 | 32.7 | 150.1 | 76 | (hi-inc; aimpoint targeted to extended-mission Titan) |

**Targeted-flyby count check:** 33 Titan + Rhea 4 + Dione 7 + Iapetus 13 + Tethys 21 + Enceladus 25 = 38 targeted. Consistent with §5.

**Nontargeted count:** 16 in Table 2 (every row with `N` suffix). The paper does not give the total nontargeted count.

**Bodies visited (Table 2 sequence_canonical equivalent):** Titan (33×), Rhea (1 targeted + 4 nontargeted), Dione (1 targeted + 5 nontargeted), Iapetus (1 targeted), Tethys (1 targeted + 2 nontargeted), Enceladus (1 targeted + 4 nontargeted), Mimas (1 nontargeted). **Phoebe is mentioned for pre-SOI approach but does not appear in the published tour table.**

### 3.4 What is NOT published

**V∞ at any encounter is not tabulated.** Table 2 publishes (altitude, latitude, west-longitude, post-flyby inclination) — the satellite-relative aimpoint geometry — but no V∞ magnitudes anywhere. The paper's only V∞ statement is qualitative (§1 p. 1611, "A single targeted flyby can change the orbiter's Saturn-relative velocity by hundreds of m/s"). Equation (1) p. 1612 (the bending equation) is given in terms of `r_p, V∞, μ` symbolically without numerical V∞ values for any specific flyby.

**Per-phase ΔV is not tabulated.** Only the ~500 m/s tour total appears.

**Tour ID is NOT given.** The paper calls it "a sample tour" or "the sample tour" throughout (§1 p. 1611, §5 first sentence). There is no internal designation like "T18-5" or analogous to Heaton-Longuski's "U00-01". The closest the paper comes is calling it "the sample tour presented here" (§6 p. 1619). Any tour-ID label we apply (e.g. for a catalogue row id) must be a project-internal designation, NOT a paper-internal one.

**Time of flight per leg is not separately tabulated**, but is trivially derivable from Table 2 date/time deltas between consecutive same-tour encounters.

## 4. Catalogue admission verdict

### 4.1 V0 admissibility test

The V0 evidence standard per spec §16.7.12 / §14 (per the brief and the Heaton-Longuski admission precedent) for an `mga_tour` row is: **per-flyby epoch + V∞ tuple + body sequence are minimum.**

Wolf-Smith 1995 publishes:

- Per-flyby epoch: YES (Table 2, second column, full date+time to seconds).
- Body sequence: YES (Table 2, first column, 54 encounters).
- V∞ tuple: **NO**. Table 2 has no V∞ column. The paper does not publish V∞ at any individual flyby.

Verdict: **the paper as published is one column short of the V0 mga_tour minimum.** It is per-flyby-rich on aimpoint geometry (altitude / lat / W.lon / post-flyby inclination) but per-flyby-poor on V∞.

### 4.2 Three options

**Option A — Admit at V0 with V∞ derived from real-eph SPICE.** The encounter epochs are published to the second; Saturn-system SPICE kernels (SAT441 / NAIF) cover the 2004-2008 era. We have the Saturn moon SPK and SPICE infrastructure used in #344 (Saturn-Titan-Rhea Phase 2). For each Table 2 encounter, compute the orbiter heliocentric/Saturnocentric state by lying it on the published aimpoint at the published epoch (closest-approach altitude + lat/lon + the bending equation closed by the previous-leg energy), then propagate to compute V∞ at the satellite. **This is derived-not-sourced for V∞** and therefore violates the [golden-tests sourced-only feedback](MEMORY.md). Not acceptable at V0.

**Option B — Admit at V0 with V∞ from a sibling paper.** The Cassini special section likely has a sibling paper that publishes V∞ values for the same sample tour (the paper's §1 cross-reference to "Section 4" hints at this). The post-launch trio in `papers/` (Strange-Russell-Buffington 2007, Yam-Davis-Longuski-Howell-Buffington 2009, Valerino 2014) describes the **post-launch flown tour**, not the 1995 sample tour, so the V∞ values would not necessarily match Wolf-Smith 1995 Table 2. The 1995 pre-launch sample tour was **superseded** by the actual flown tour years before SOI; per-flyby V∞ values for this specific pre-launch sample tour may exist only in Cassini-project internal docs, not in the open literature. **Verdict: NEEDS-MORE-DATA. The required sibling V∞ source has not been located.**

**Option C — Admit Wolf-Smith 1995 at V0 with the body sequence + epoch tuple ONLY, declaring V∞ as "not published by Wolf-Smith 1995; published only as ~500 m/s tour total".** This is a relaxation of the v4.7 mga_tour V0 standard (the Heaton-Longuski precedent had Table-3 V∞ at heliocentric encounters and Table-5 V∞ at every Uranian-tour event). **Verdict: would set a precedent that mga_tour rows can omit V∞ entirely, which weakens the schema.** Not recommended without #345 explicit ratification of a lower V0 bar for "epoch-and-aimpoint-only" tour designs.

### 4.3 Recommendation

**RECOMMEND: do NOT admit Wolf-Smith 1995 as a V0 `mga_tour` row in this session.** The paper is missing the V∞ column the Heaton-Longuski precedent established as the V0 minimum. Instead:

1. **File Wolf-Smith 1995 as a KNOWN_CORPUS literature anchor** for the Cassini Titan-pump tour design (§5 of this digest). The paper is the canonical pre-launch reference for Cassini-tour topology and should be cited; it does NOT have to be a catalogue row to play that role.
2. **Open a follow-up acquisition task** for the sibling-paper Cassini special section + the Cassini Mission Plan (JPL D-5564 series; the 1990s/2000s flight-design internal-but-citable docs Strange-Russell-Buffington cites) to find a published V∞ table for either the 1995 sample tour or the flown tour.
3. **If a V∞ table surfaces**, admit a `cassini-wolf-smith-1995-sample-tour` row OR a `cassini-flown-tour-2004-2008` row at V0 with the structural template below.

### 4.4 Structural template (for future use, NOT for writeback now)

If/when V∞ data surfaces, the row construction would look like:

```yaml
- id: cassini-wolf-smith-1995-sample-tour
  name: "Cassini sample satellite tour (Wolf-Smith 1995, pre-launch reference)"
  source: literature
  trajectory_regime: ballistic
  model_assumption: analytic-ephemeris  # patched-conic + third-body perturbations per Wolf-Smith §4 optimizer
  cycler_class: multi-arc  # 63 orbits + 38 targeted satellite flybys = 38+ distinct heliocentric/Saturnocentric arcs
  orbit_class: mga_tour
  epoch_locked: true
  n_returns: 1  # single SOI / single tour completion, no cycler return
  validity_window:
    start: "2004-07-01T00:00:00Z"  # SOI per §3 p. 1615
    end:   "2008-07-01T00:00:00Z"  # tour end per §5 p. 1619
  launch_epoch: "2004-07-01T00:00:00Z"  # tour start (SOI); interplanetary launch was 1997 but pre-SOI is outside the tour
  validation_level: V0
  source_ephemeris: "Wolf-Smith 1995 §4 optimizer + precision integrator (Saturn J2 + Sun third-body)"
  bodies: ["Saturn", "Titan", "Rhea", "Dione", "Iapetus", "Tethys", "Enceladus", "Mimas"]
  sequence_canonical: <encoded from Table 2 54-row sequence>
  vinf_kms_at_encounters: <NEEDS-MORE-DATA: not published>
  first_published:
    authors: ["Wolf, A. A.", "Smith, J. C."]
    year: 1995
    title: "Design of the Cassini Tour Trajectory in the Saturnian System"
    venue: "Control Engineering Practice, Vol. 3, No. 11, pp. 1611-1619"
    doi: "10.1016/0967-0661(95)00172-7"
```

This template should be written **only** once the V∞ gap is closed. Do not commit a placeholder row.

## 5. KNOWN_CORPUS impact

The existing Cassini anchor in `src/cyclerfinder/search/literature_check.py` (lines 1010-1045) is:

- Name: `"Cassini-Huygens Saturn-Titan satellite tour design"`
- body_set: `{Titan, Enceladus, Rhea, Dione, Iapetus}` — **does NOT include Tethys or Mimas**, which Wolf-Smith 1995 Table 2 confirms are also in the sample tour (Tethys 21 targeted + Tethys 46N/49N nontargeted; Mimas 31N nontargeted).
- topology_label: `{pump-tour, mga-tour}` — correct per Wolf-Smith 1995.
- Citation: Strange-Russell-Buffington 2007 AAS-07-277 + Yam-Davis-Longuski 2009 JSR + Valerino 2014 SpaceOps — all **post-launch** sources.
- DOI: Yam-Davis-Longuski 2009.

**RECOMMEND two amendments:**

1. **Add Wolf-Smith 1995 as the LEADING pre-launch citation**, with the Strange/Yam/Valerino trio retained as post-launch refinement references. Wolf-Smith 1995 is the canonical pre-launch base design and should anchor the citation. Suggested rewording of the `citation` field: `"Wolf & Smith, 'Design of the Cassini Tour Trajectory in the Saturnian System,' Control Engineering Practice 3(11):1611-1619 (1995) -- pre-launch sample tour: 63 orbits, 38 targeted flybys, 33 of Titan + 5 of (Enceladus, Tethys, Dione, Rhea, Iapetus), tour ID not given by paper; Strange, Russell & Buffington, 'Mapping the V-infinity globe' (AAS 07-277, 2007) -- same-body Titan-pump method used in Cassini extended mission; Yam, Davis, Longuski, Howell & Buffington, 'Saturn Impact Trajectories for Cassini End-of-Mission,' JSR DOI 10.2514/1.38760 (2009); Valerino, 'Updating the Reference Trajectory for the Cassini Solstice Mission,' SpaceOps 2014 DOI 10.2514/6.2014-1880."` DOI field: prefer the Wolf-Smith DOI `10.1016/0967-0661(95)00172-7`.

2. **Extend body_set to include Tethys and Mimas** (Wolf-Smith Table 2 confirms both). Suggested: `frozenset({"Titan", "Enceladus", "Rhea", "Dione", "Iapetus", "Tethys", "Mimas"})`. Phoebe is also mentioned by Wolf-Smith (pre-SOI approach flyby, §3) — include it if the anchor is meant to cover the full Cassini mission rather than just the tour proper.

3. **Add Wolf-Smith 1995 to authors tuple**: `("Wolf", "Smith", "Strange", "Russell", "Buffington", "Yam", "Davis", "Longuski")`.

These are RECOMMENDATIONS — not edits in this session.

## 6. Errata

Versus the Agent E task brief and our internal notes:

- **Brief said "Cassini tour visited Titan, Enceladus, Iapetus, Rhea, Dione, Hyperion, Tethys, Phoebe."** Wolf-Smith 1995 Table 2 includes **Mimas** (orbit 31N) but **does NOT include Hyperion** in the sample tour. Hyperion appears nowhere in the paper. Phoebe is mentioned in §3 (pre-SOI approach) but does NOT appear in Table 2. So the sample-tour bodies are {Titan, Rhea, Dione, Iapetus, Tethys, Enceladus, Mimas}, not the brief's set. (The actual flown Cassini mission did include a Hyperion flyby and a Phoebe flyby, but the 1995 sample tour as published does not.)
- **Brief said "Titan ~45+ flybys for pump."** Wolf-Smith 1995 sample tour has **33 Titan flybys (targeted)**, per §5 p. 1616 explicit count. The 45+ number must be from the flown tour (with extended mission), not the 1995 sample.
- **Brief said abstract describes "Huygens probe release into Titan atmosphere on approach."** The abstract itself says "a Saturn orbiter and a Titan probe" without "on approach"; the §1 narrative describes release after SOI and approach to the first Titan flyby, not on the interplanetary approach to Saturn. Minor framing nuance.
- **Brief said paper is 9 pages** — confirmed (1611-1619).
- **Paper-internal cross-reference error:** §1 last paragraph says "The probe mission is described in detail in Section 4" but §4 of this paper is "Methods and Software", not the probe mission. This is almost certainly a stale cross-reference from an editing pass — likely meant to refer to a sibling paper in the same Cassini special section.

## 7. Action items

For the parent / #349 / #345:

1. **DO NOT writeback an `mga_tour` row for Wolf-Smith 1995 in this session.** V∞ is missing; the paper falls one column short of the v4.7 V0 mga_tour minimum standard.

2. **(Action — KNOWN_CORPUS amendment)** Open a small follow-up task to amend the existing `Cassini-Huygens Saturn-Titan satellite tour design` anchor (`src/cyclerfinder/search/literature_check.py` line 1011) to:
   - lead with Wolf-Smith 1995 (DOI `10.1016/0967-0661(95)00172-7`)
   - extend body_set to add Tethys + Mimas
   - extend authors tuple to add Wolf + Smith
   - keep Strange/Yam/Valerino as post-launch refinements
   - keep topology_label `{pump-tour, mga-tour}` (Wolf-Smith confirms).

3. **(Action — acquisition gate)** Acquisition task: locate a Cassini-tour V∞ table. Candidates:
   - The 1995 *Control Engineering Practice* Cassini special section sibling papers (the §1 cross-reference suggests one exists).
   - Cassini Mission Plan (JPL internal but commonly cited; D-5564 / D-11431 series).
   - Smith-Buffington 1996 AAS papers on the Cassini tour (likely successors to Wolf-Smith 1995 with V∞ tables).
   - **NOT** Strange-Russell-Buffington 2007 / Yam 2009 / Valerino 2014: these describe the post-launch flown tour and the end-of-mission refinements, not the 1995 sample tour or its V∞ structure.

4. **(Action — Galileo parity)** Agent D's parallel D'Amario-Bright-Wolf 1992 Galileo digest is on the matched paper. The Galileo Wolf-Byrnes 1993 AAS-93-567 is referenced by Wolf-Smith 1995 §1 as the foundational tour-design paper — if our acquisitions list does not have AAS-93-567, open an acquisition task for it. (Cross-check with Agent D's digest before duplicating.)

5. **(Action — body_set audit)** The Tethys + Mimas omission in the current KNOWN_CORPUS anchor body_set was a real gap (#349 added the topology_label discriminator but the body_set itself was incomplete). Sweep other moon-tour anchors for similar omissions when the underlying source is consulted.

End of Agent E digest.
