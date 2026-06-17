# Dunne-Burgess 1978 Mariner-10 deep-read digest (#382 paper 2)

Per #382 task brief 2026-06-17 AET (#345 / #361 admission candidates for `mga_tour`). Textbook scope discipline applied to a 168-page NASA Gov't-Printing-Office popular-science publication: chapters with mission-trajectory content read fully (1, 2, 5, 6, 7, 8); chapters on spacecraft engineering (3, 4) sampled for ΔV / propulsion data; chapter 9 + appendices sampled for additional epoch references.

## 1. Header (verbatim)

- **Title:** "The Voyage of Mariner 10: Mission to Venus and Mercury"
- **Authors:** James A. Dunne, Eric Burgess
- **Publisher:** National Aeronautics and Space Administration, Scientific and Technical Information Office, Washington, D.C., 1978
- **Document ID:** NASA SP-424. Prepared by Jet Propulsion Laboratory, California Institute of Technology.
- **Source format:** This PDF is a NASA-History-Office HTML-to-PDF conversion (note "PDF makers note: Book converted from webpage (HTML) to pdf by PDF-creator, and supplied as is" on pp. 2, 14, 30). The original 1978 paper book has more page-level data; the HTML conversion preserves the text but loses some tabular formatting if any was present in the original.
- **Page count (PDF):** 168 pages (cover + 6 frontmatter + 9 chapters + 4 appendices + index).

## 2. What the paper actually is

A **NASA-History-Office popular-science book** documenting the Mariner-10 mission (launch 1973-11-03; Venus flyby 1974-02-05; three Mercury encounters 1974-03-29, 1974-09-21, 1975-03-16; end-of-mission 1975-03-24). The book is written for an interested public audience — it has substantial astronomy / planetary-science narrative content (Chapter 1 is about ancient Greek astronomy; Chapter 9 is about planetary science synthesis). Chapters 3 and 4 are spacecraft engineering for general readers (instrument descriptions, thermal control, the "sail vs V-tilt solar panel" tradeoff).

**The trajectory content (Chapters 2, 5, 7, 8) is qualitative and event-narrative, not numerical-table-driven.** Specifically: there is **no V_∞ table at any flyby**. The book publishes:

- Launch date and launch-window length (Chapter 4 p. 29).
- Closest-approach altitudes for each flyby (qualitative + a few sourced numbers).
- Closest-approach dates and times in PDT.
- Trajectory-correction-maneuver (TCM) ΔV magnitudes (Chapter 5, 7, 8 — small, ~1-18 m/s each).
- Trip-totals (~500 days, "more than a billion kilometers").

**The book does NOT publish:**

- V_∞ at Earth departure.
- V_∞ at Venus flyby.
- V_∞ at any of the three Mercury flybys.
- Heliocentric ToF per arc.
- C3 at launch.
- Bend-angle / aim-point / B-plane parameters.

The Chapter 2 "Suggestions for Further Reading" + Appendix C (Spacecraft and Science Teams) + Appendix D (Award Recipients) **do not list the JPL Mariner-10 Mission Plan or any flight-design AIAA paper** that would contain the per-flyby V_∞ data. This is consistent with the book's popular-science scope.

The chapters in detail:

1. **Chapter 1 "Earth's Sister and the Twilight Planet"** (pp. 1-9 / PDF pp. 13-22) — Ancient astronomy of Mercury and Venus. No trajectory data.

2. **Chapter 2 "Mariner Venus-Mercury Mission"** (pp. 10-17 / PDF pp. 23-32) — Mission concept history; Mariner-spacecraft lineage; launch vehicle Atlas SLV-3D/Centaur D-1A; spacecraft mass 533.6 kg (1175 lb). Mentions "the 1973 opportunity offered one of the lowest launch energies to swing by Venus and subsequently encounter Mercury" (Chapter 4 p. 29) — qualitative C3 only.

3. **Chapter 3 "Mariner's Payload"** (pp. 18-27 / PDF pp. 32-43) — Science instruments (TV imaging, IR radiometer, plasma science, magnetometer, charged-particles, EUV, radio occultation). Spacecraft engineering, no trajectory.

4. **Chapter 4 "Spacecraft, Scientists, and Schedules"** (pp. 28-43 / PDF pp. 44-59) — Launch window opens "November 2, 1973 (November 3, on the East Coast), 1.5 hours" (p. 29). Project formally initiated December 1969; 4-year development. Spacecraft trim ΔV capability: "originally 56 m/sec, increased to 122 m/sec (401 ft/sec)" by larger Pioneer-10/11 propellant tank (p. 30). Nitrogen gas tank size 2.45 → 3.62 kg.

5. **Chapter 5 "Venus Bound — Success and Near Failure"** (pp. 44-58 / PDF pp. 60-74) — Mission narrative from launch through Venus encounter. **Trajectory-correction maneuver details:**
   - TCM-1: shortly after launch (Nov 1973), 17.41 Hz Doppler shift = "about 1.3 m/sec (4 ft/sec)" velocity change (p. 56 / Ch.5 page 12).
   - Mariner 10 was "640,000 km (about 400,000 mi) from Venus and approaching the planet at a speed of over 29,600 km/hr (18,400 mi/hr)" on February 4, 1974 (p. 58 / Ch.5 page 13).
   - **This 29,600 km/hr ≈ 8.22 km/s is the Venus-relative approach speed at 640,000 km — NOT the asymptotic V_∞ at Venus, but a finite-distance approximation that captures most of V_∞ because 640,000 km is well outside Venus's sphere of influence (~616,000 km). This is the closest the book comes to a Venus V_∞.**
   - TCM-2: January 21, 1974, 24-min after rolls; 3.8-sec rocket burn; "within 27 km (17 mi) of the aim point" (p. 57). ΔV magnitude not explicitly given in m/s in this paragraph (the 1.3 m/s value reported earlier was TCM-1 result).
   - Venus closest approach: **5,794 km (3,600 mi)** at **10:01 a.m. PDT February 5, 1974** (p. 62 / Ch.6 page 1). Approach was from the dark side; passed over the sunlit side; Venus gravity bent the path from the original 16,000 km offset to the 5,794 km flyby (Chapter 5 p. 57).
   - Bias / desired flyby point: per Fig. 5-10 the "biased injection aim point" was 55,000 km (3 hrs late if uncorrected); TCM-1 brought it to 1,380 km off and 2 min early; TCM-2 brought it to the final desired 5,784 km / 27-km error.

6. **Chapter 6 "Best Seen in Black Light"** (pp. 59-69 / PDF pp. 75-85) — Venus science results. No new trajectory numbers beyond Chapter 5's 5,794 km closest approach. Mariner-10 was acquired by Goldstone 64-m at 45 million km from Earth at 9:21 a.m. PDT Feb 5; closest approach at 10:01 a.m. PDT; Earth occultation 10:07 a.m. PDT (= 6 min after closest approach).

7. **Chapter 7 "Mercury, Moonlike and Earthlike"** (pp. 70-87 / PDF pp. 86-103) — Cruise to Mercury and **first Mercury encounter (Mercury I)**.
   - **TCM-3:** "On March 16, [1974], at 04:54 a.m. PDT, the propulsion system was ignited and burned for 51 sec to change the velocity of Mariner by **17.8 m/sec (59 ft/sec)** directly away from the Sun. This would change the Mercury flyby from the sunlit to the dark side of the planet" (p. 72).
   - **Mercury I closest approach:** "The flyby was expected to be 200 km (124 mi) closer to Mercury than planned. ... All ... requirements of the science experiments at Mercury, no additional maneuvers were planned" (p. 72-73). The actual closest approach was at **1:46 p.m. PDT on March 29, 1974, with the flight path targeted to pass behind the planet on the night side** (p. 75). Mercury I closest-approach altitude is documented in the Appendix A approach figures (Fig. A-3 caption mentions the 20,700 km distance "a half-hour before Mariner made its first close flyby of Mercury, March 1974" / p. 112).
   - **Note:** the explicit Mercury I closest-approach altitude is NOT printed in Chapter 7 (the book describes only that the flyby was "200 km closer to Mercury than planned" and the planned altitude is not given in the text). Other sources (NASA fact sheets) give Mercury I altitude as 703 km (or 700 km after the trim), but the Dunne-Burgess text does not source this number directly.

8. **Chapter 8 "Return to the Innermost Planet"** (pp. 88-99 / PDF pp. 104-115) — **Mercury II and Mercury III encounters.**
   - **TCM-4** (May 9, 1974, p. 89) and **TCM-5** (July 1974, p. 89-90) — described qualitatively; ΔV not explicitly given in m/s here.
   - **Mercury II closest approach:** Fig. 8-4 caption: aim point AFTER TCM-5 was at 50,000 km (Mercury II distant flyby for sunlit-side TV coverage). **"Point of closest approach during the second encounter occurred at 1:59 p.m. PDT on September 21, 1974"** (p. 91). The 50,000-km altitude is explicit in Fig. 8-4.
   - The Mercury II orbit period was tuned so the spacecraft "would [re-encounter Mercury] approximately 17 min later than the time desired" after TCM-3 — the implication is Mercury's orbital period of 88 days × 2 ≈ 176 days = ~the actual Mercury I → Mercury II gap of 176 days (March 29 → September 21, 1974 = 176 d). **This is the famous 2:1 Mariner-10 resonance** (spacecraft heliocentric period = 2 × Mercury's heliocentric period — but this resonance relationship is not stated in the book in those terms; the book just narrates the events).
   - **Mercury III:** "the closest planetary flyby yet accomplished" (p. 96). **"On March 16, 1975"** Mercury III happened (Fig. 8-10 caption sources March 16, 1975; one picture taken from 67,000 km, another from 65,000 km, and finally from 19,000 km "34 min after the spacecraft swept past Mercury"). The exact closest-approach altitude is not given numerically in the text. Magnetometer table on p. 98 gives bow shock at 3:31 PDT, magnetopause 3:39, maximum field 3:49, magnetopause 3:54/3:56, bow shock 3:58/3:59 — so Mercury III closest approach was at ~03:49 PDT on March 16, 1975.
   - **End-of-mission:** "The end came on March 24, 1975, when the final depletion of the nitrogen supply was signalled by the onset of an unprogrammed pitch turn" (p. 98).

9. **Chapter 9 "A Clearer Perspective"** (pp. 100-104 / PDF pp. 117-122) — Planetary-science synthesis. No trajectory data.

10. **Appendices A, B, C, D** (PDF pp. 122-159) — Mercury mosaics, image processing, science teams, award recipients. No trajectory data; some incidental altitudes appear in figure captions but they are observation distances at picture-taking time, not closest-approach altitudes (e.g., Fig. A-1 photomosaic taken "13-min period when Mariner was 200000 km from Mercury", Fig. A-2 mosaic "taken from a distance of 210,000 km" 6 hours after the flyby).

## 3. Per-encounter data extracted

### 3.1 Body sequence (canonical)

E → V → M → M → M, with Earth launch 1973-11-03 and three Mercury encounters in the famous 2:1 spacecraft:Mercury heliocentric resonance.

### 3.2 Epoch tuple (sourced)

| # | Encounter | Date | Time (PDT) | Closest-approach altitude (sourced) | Source page |
|---|---|---|---|---|---|
| 0 | Earth launch | 1973-11-03 | (1.5-hr launch window; East Coast date 1973-11-03; West Coast 1973-11-02) | n/a (launch) | p. 29 / Ch.4 |
| 1 | Venus flyby | 1974-02-05 | 10:01 a.m. | **5,794 km (3,600 mi)** | p. 62 / Ch.6 |
| 2 | Mercury I | 1974-03-29 | 1:46 p.m. | NOT EXPLICITLY GIVEN (text says "200 km closer than planned"; planned altitude not stated) | p. 73, 75 |
| 3 | Mercury II | 1974-09-21 | 1:59 p.m. | **~50,000 km** (Fig. 8-4 aim point post-TCM-5) | p. 91 / Ch.8 |
| 4 | Mercury III | 1975-03-16 | ~03:49 a.m. (magnetometer max field; closest approach) | NOT GIVEN ("closest planetary flyby yet accomplished") | p. 96-98 / Ch.8 |
| 5 | End of mission | 1975-03-24 | n/a | n/a | p. 98 |

### 3.3 ΔV tuple (sourced) — these are trim maneuvers, NOT mission-design deterministic ΔV

| # | Maneuver | Date | ΔV (m/s) | Notes |
|---|---|---|---|---|
| TCM-1 | trim post-launch | 1973-11-13 (approx; "shortly after launch") | 1.3 m/s (4 ft/s) | (technically the measured Doppler-shift result of TCM-1) — p. 56 / Ch.5 |
| TCM-2 | pre-Venus | 1974-01-21 | NOT EXPLICITLY GIVEN (3.8-s burn; 27-km flyby-point error remaining) | p. 57 / Ch.5 |
| TCM-3 | post-Venus, pre-Mercury I (Sun-line) | 1974-03-16 04:54 PDT | **17.8 m/s (59 ft/s)** | retrograde from Sun, 51-s burn; moved Mercury I from sunlit to dark side | p. 72 / Ch.7 |
| TCM-4 | post-Mercury I | 1974-05-09 | NOT EXPLICITLY GIVEN | p. 89 / Ch.8 |
| TCM-5 | pre-Mercury II | 1974-07 (mid-July) | NOT EXPLICITLY GIVEN (50,000-km aim point shift per Fig. 8-4) | p. 89-90 / Ch.8 |
| TCMs-6,7,8 | pre-Mercury III | 1974-10 → 1975-03 | NOT EXPLICITLY GIVEN ("three trajectory correction maneuvers were successfully completed during this period") | p. 96 / Ch.8 |
| Solar sailing | continuous | 1974-04 → 1975-03 | (no propellant; attitude only) — quantitative gas savings: "to some 25% of normal cruise usage" | p. 95 / Ch.8 |

### 3.4 V_∞ — NOT PUBLISHED

**Critical finding:** the book does NOT publish V_∞ at any encounter. The Venus-approach speed of **29,600 km/hr (8.22 km/s) at 640,000 km** (p. 58) is the closest the book comes to a V_∞ value, and it is a Venus-relative speed at finite distance, not the asymptote.

Numeric trajectory quantities that ARE published:
- Total mission distance: "more than a billion kilometers" (Foreword p. iv).
- Mission duration to Mercury I: ~500 days (Foreword p. iv: "a little over 500 days").
- Earth → Venus distance traveled: "236 million kilometers" by Feb 5, 1974 (Foreword p. iv).
- Mercury orbital period implication: Mercury I → Mercury II gap = 176 d ≈ 2 × Mercury's 87.97-day period. Mercury I → Mercury III gap = 352 d ≈ 4 × Mercury's period. (The book does NOT name this as a 2:1 resonance; this is the well-known Mariner-10 design but the book describes it only as "approximately 17 min later than the time desired" after TCM-3.)

## 4. Catalogue admission verdict

### 4.1 V0 admissibility test

Per the Heaton-Longuski 2003 U00-01 admission precedent and the Wolf-Smith 1995 Cassini negative-verdict precedent and the Bourke 1971 sibling-paper-1 verdict in this same digest sweep, the V0 evidence standard for an `mga_tour` row is:

**Per-flyby epoch + V_∞ tuple + body sequence are minimum.**

Dunne-Burgess 1978 publishes:

- Per-flyby epoch: YES for Venus (1974-02-05 10:01 PDT), Mercury I (1974-03-29 1:46 PDT), Mercury II (1974-09-21 1:59 PDT). PARTIAL for Mercury III (1975-03-16 ~03:49 PDT, inferred from magnetometer bow-shock timing).
- Body sequence: YES (E → V → M → M → M).
- V_∞ tuple: **NO**. No V_∞ at any encounter. The 8.22 km/s Venus-approach speed at 640,000 km is the closest data point, but it is finite-distance not asymptotic.

**Verdict: the book is ONE COLUMN SHORT of the V0 `mga_tour` minimum** (same shortfall as Wolf-Smith 1995 Cassini and Bourke 1971 — third negative verdict in this digest sweep).

### 4.2 Three layered problems

1. **V_∞ absent.** Same blocker as Wolf-Smith 1995 and Bourke 1971. V_∞ could in principle be back-derived from SPICE / DE440 reconstruction at the published epochs (Mercury and Venus SPICE kernels cover 1973-1975 trivially), but that is **derived, not sourced** — violates the [`golden-tests sourced-only`](MEMORY.md) rule.

2. **Mercury I and Mercury III closest-approach altitudes not published.** Mercury I is described only as "200 km closer than planned" without the planned altitude. Mercury III is described as "the closest planetary flyby yet accomplished" without a number. Other widely-cited sources (NASA fact sheets, the Mariner-10 mission summary, and Dunne's own 1974 *Mariner 10 Mercury Encounter* Science article) give Mercury I = 703 km, Mercury III = 327 km, but those numbers are not in THIS book.

3. **TCM ΔVs incomplete.** Only TCM-1 (1.3 m/s) and TCM-3 (17.8 m/s) are given numerically. TCM-2, TCM-4 through TCM-8 are described qualitatively without a m/s number. Per-leg deterministic ΔV totals are not given.

### 4.3 Recommendation

**RECOMMEND: do NOT admit Mariner-10 from Dunne-Burgess 1978 as a `mga_tour` row in the catalogue from this paper alone.** V_∞ data is missing (V0 minimum fails), and Mercury I/III altitudes are missing.

Instead:

1. **File Dunne-Burgess 1978 as a KNOWN_CORPUS literature anchor** for the Mariner-10 mission (popular-science archetype; canonical NASA-History-Office reference for the mission narrative and epoch tuple). This is a high-quality CITABLE PRIMARY SOURCE for the epoch tuple (Venus 1974-02-05, Mercury I 1974-03-29, Mercury II 1974-09-21, Mercury III 1975-03-16) even though the V_∞ data is missing.

2. **Open follow-up acquisition tasks** for the Mariner-10 flight-design papers:
   - **Bourke-Beerer 1970 "Mariner Venus/Mercury 1973 Preliminary Mission Design"** AIAA-70-1049 (referenced by Bourke 1971 as Ref. 7 — see `/home/bruce/dev/cyclers/docs/notes/2026-06-17-digest-bourke-1971.md` §1). This is the pre-launch JPL Mariner-10 trajectory design and likely has the C3 + V_∞ table.
   - **Dunne 1974** "Mariner 10 Mercury Encounter" *Science* 185:141-142.
   - **JPL Mariner Venus/Mercury 1973 Mission Plan** (D-series JPL internal but commonly cited).
   - **Stuart-Roy 1976 (or similar)** mission-reconstruction papers on the actual flown Mariner-10 trajectory.

3. **If a V_∞ table surfaces** (most likely from Bourke-Beerer 1970 AIAA-70-1049 — which is in the parallel #382 wave of acquisitions), admit a `mariner-10-vmm-1973` row at V0 with the structural template below.

### 4.4 Structural template (for future use, NOT for writeback now)

If/when V_∞ data surfaces (most likely from Bourke-Beerer 1970 or Dunne 1974):

```yaml
- id: mariner-10-vmm-1973-grand-tour
  name: "Mariner 10 Earth-Venus-Mercury-Mercury-Mercury tour (1973-1975 flown)"
  source: literature
  trajectory_regime: ballistic   # heliocentric Earth-Venus arc + Venus-Mercury arc + 2:1 Mercury-resonant orbit (no propulsion after TCM-3 in nominal mission)
  model_assumption: analytic-ephemeris   # patched-conic for design; DE-series ephemeris for reconstruction
  cycler_class: multi-arc   # 4 distinct heliocentric arcs (E-V, V-M1, M1-M2, M2-M3) plus end-of-mission
  orbit_class: mga_tour   # schema v4.7; Venus + 3-Mercury sequence is the canonical first-use of gravity assist for a flown mission
  epoch_locked: true   # the 1973-11-03 launch window was selected for the favorable Earth-Venus-Mercury geometry; not repeatable
  n_returns: 1   # single Earth-launch + 3 Mercury encounters + end-of-mission; not a cycler
  validity_window:
    start: "1973-11-03T00:00:00Z"   # launch
    end:   "1975-03-24T00:00:00Z"   # end-of-mission
  launch_epoch: "1973-11-03T00:00:00Z"
  validation_level: V0
  source_ephemeris: <to be sourced from Bourke-Beerer 1970 or Dunne 1974>
  bodies: ["E", "V", "Mercury"]
  sequence_canonical: "E-V-Mercury-Mercury-Mercury"
  vinf_kms_at_encounters: <NEEDS-MORE-DATA: not published in Dunne-Burgess 1978>
  encounter_epochs:
    - body: "V"
      epoch: "1974-02-05T18:01:00Z"   # 10:01 PDT = 18:01 UTC per Dunne-Burgess p.62
      closest_approach_km: 5794
      source_quote: "Mariner 10 made its closest approach of 5794 km (3600 mi) at 10:01 a.m. PDT" (Dunne-Burgess 1978 p. 62)
    - body: "Mercury"  # Mercury I
      epoch: "1974-03-29T20:46:00Z"   # 1:46 PDT
      closest_approach_km: <NOT IN BOOK — Dunne 1974 Science gives ~703 km but verify>
      source_quote: "the closest approach of about 5790 km (3600 mi)" — but THIS IS A TYPO in Dunne-Burgess p.63 (confused with Venus); Mercury I was 703 km per NASA fact sheets
    - body: "Mercury"  # Mercury II
      epoch: "1974-09-21T20:59:00Z"   # 1:59 PDT
      closest_approach_km: 50000   # Fig. 8-4 aim point post-TCM-5
      source_quote: "Point of closest approach during the second encounter occurred at 1:59 p.m. PDT on September 21, 1974" (p. 91)
    - body: "Mercury"  # Mercury III
      epoch: "1975-03-16T10:49:00Z"   # ~03:49 PDT inferred from magnetometer max-field timing
      closest_approach_km: <NOT IN BOOK — "closest planetary flyby yet accomplished", NASA fact sheets give 327 km>
      source_quote: "On March 16, 1975 ... a multiple impact feature" (Fig. 8-10) + magnetometer max field 3:49 PDT (p. 98)
  first_published:
    authors: ["Bourke, R. D.", "Beerer, J. G."]
    year: 1970
    title: "Mariner Venus/Mercury 1973 Preliminary Mission Design"
    venue: "AAS/AIAA Astrodynamics Conference, Santa Barbara, CA, August 19-21, 1970"
    paper_id: "AIAA Paper 70-1049"
  corroborating_sources:
    - authors: ["Dunne, J. A.", "Burgess, E."]
      year: 1978
      title: "The Voyage of Mariner 10: Mission to Venus and Mercury"
      venue: "NASA SP-424, NASA History Office"
    - authors: ["Dunne, J. A."]
      year: 1974
      title: "Mariner 10 Mercury Encounter"
      venue: "Science 185(4146):141-142"
```

Dunne-Burgess 1978 would then be a `corroborating_source` for the Mariner-10 row but NOT the `first_published` row (Bourke-Beerer 1970 is the actual flight design); per the same logic as Bourke 1971, the popular-science book sits below the AIAA flight-design paper in the citation hierarchy.

## 5. KNOWN_CORPUS impact

Search the existing literature_check.py for Mariner-10 / Mercury / Venus-Mercury anchors. If none, adding a Mariner-10 anchor is reasonable:

- Name: `"Mariner-10 Earth-Venus-Mercury tour (Bourke-Beerer / Dunne-Burgess 1978)"`
- body_set: `frozenset({"E", "V", "Mercury"})`
- topology_label: `{"mga-tour"}` — Mariner-10 is the canonical Venus-Mercury gravity-assist mission and the **first-ever flown gravity-assist**.
- Citation: `"Bourke & Beerer, 'Mariner Venus/Mercury 1973 Preliminary Mission Design,' AIAA Paper 70-1049 (1970) -- pre-launch JPL flight design; Dunne & Burgess, 'The Voyage of Mariner 10: Mission to Venus and Mercury,' NASA SP-424 (1978) -- canonical mission narrative + flown epoch tuple (Venus 1974-02-05, Mercury I 1974-03-29, Mercury II 1974-09-21, Mercury III 1975-03-16, end-of-mission 1975-03-24)."` DOI: Bourke-Beerer has no public DOI; cite by AIAA Paper number.
- Authors: `("Bourke", "Beerer", "Dunne", "Burgess")`.

This anchor would prevent a future false-novelty claim for "first gravity-assist mission" or "Venus-Mercury tour" searches.

## 6. Errata

Versus the task brief and the book text:

- **Brief said "extract per-encounter V_∞ + ToF + body sequence + epoch for the Mariner-10 V-V-M tour."** The book does NOT publish V_∞. The Foreword identifies the body sequence as "V-V-M" but the trajectory is actually **V-M-M-M** (one Venus + three Mercury). Mariner-10 was a Venus-Mercury(-Mercury-Mercury) tour, not a Venus-Venus-Mercury tour. The "V-V-M" string in the brief is a typo (Venus-Mercury-Mercury is at least the first two encounters past Venus; the third Mercury was the extended-mission third look).
- **Brief said "If the V_∞ data isn't in this book (it's a Gov't popular-science publication), REPORT that and flag the need for the actual flight design paper (Dunne 1972 'Mariner 10 Spacecraft' or similar)."** REPORTED. The flight-design paper is **Bourke-Beerer 1970 AIAA Paper 70-1049** (the brief's "Dunne 1972 Mariner 10 Spacecraft" reference does not match a standard JPL flight-design paper title; the correct citation appears to be Bourke-Beerer 1970, which is referenced by paper 1 of this digest set — Bourke 1971 §IX Ref. 7).
- **Page-43 oddity:** Chapter 5 page 13 in the PDF says "the cameras took the first picture of the planet ... at 9:50 a.m. PDT this picture was displayed on the monitor screens. The photo showed the lighted cusp of Venus at the north pole (Fig. 6-1) just 12 min before Mariner 10 made its closest approach of about 5790 km (3600 mi) above the surface of the planet" — wait, this is **Venus** at 5790 km, not Mercury. The book is self-consistent here; Venus closest approach was 5790 / 5794 km (paper rounds inconsistently: 5790 km on p. 63 and 5794 km on p. 62). **Both numbers refer to the same Venus flyby.** Mercury I altitude is NOT 5790 km despite the apparent visual similarity to the Venus value in some other contexts.
- **Foreword p. iv mistakes:** "On February 5, 1974, after traveling 236 million kilometers, Mariner 10 skimmed past Venus within 12 kilometers of the preplanned aim point." — "12 km" is the aim-point error magnitude (TCM-2 trim quality) ROBOT not the closest-approach altitude; an inexpert reader could conflate these. The closest-approach altitude was 5,794 km per Chapter 6 p. 62.
- **Cross-paper:** Bourke 1971 (paper 1 of this digest set) Ref. 7 cites Bourke-Beerer 1970 AIAA-70-1049 as "Mariner Venus/Mercury 1973 Preliminary Mission Design" — that is the flight-design paper to acquire if Mariner-10 catalogue admission is wanted.

## 7. Action items

For the parent / #382 / #345:

1. **DO NOT writeback an `mga_tour` row for Mariner-10 from Dunne-Burgess 1978 alone.** V_∞ data is missing (V0 minimum fails); Mercury I and Mercury III closest-approach altitudes are also missing from this book.

2. **(Action — KNOWN_CORPUS amendment)** Open a small follow-up task to add a `Mariner-10 Earth-Venus-Mercury tour` literature anchor to `src/cyclerfinder/search/literature_check.py` covering E/V/Mercury with topology_label `{mga-tour}` and the citation per §5 above. This prevents future false-novelty claims for "Venus-Mercury gravity assist" / "first gravity assist mission" searches.

3. **(Action — acquisition gate)** **HIGH PRIORITY** acquisition: **Bourke-Beerer 1970 AIAA Paper 70-1049 "Mariner Venus/Mercury 1973 Preliminary Mission Design"** — this is the most likely source for per-flyby V_∞ for the Mariner-10 tour. AIAA papers from 1970 are usually available via the AIAA Aerospace Research Central archive (arc.aiaa.org). DOI form: probably `10.2514/6.1970-1049`.

4. **(Action — secondary acquisition)** Dunne 1974 "Mariner 10 Mercury Encounter" *Science* 185(4146):141-142 — likely has Mercury I closest-approach altitude in publishable form.

5. **(Action — historical anchor)** Dunne-Burgess 1978 IS a valuable epoch-tuple source for Mariner-10 (Venus and Mercury flyby dates to the minute, all sourced). When the V_∞ source surfaces, this book becomes a `corroborating_source`.

6. **(Action — typesetting / popular-science scope)** This is a useful negative datum: a NASA Gov't-Printing-Office popular-science book DOES NOT typically publish V_∞ tables. For future digest sweeps, when a paper is identified as "NASA SP-NNN" or "NASA History Office", expect popular-science scope and plan acquisition of the underlying AIAA / JPL flight-design papers in parallel.

End of digest paper 2 of 4.
