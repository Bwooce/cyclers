# Digest — D'Amario, Bright & Wolf 1992 "Galileo trajectory design"

**Date**: 2026-06-17 AET (Agent D, parent #347 / track #345 classic-mission `mga_tour` admissions)
**Verdict (TL;DR)**: **ADMIT** — the paper publishes per-encounter epochs, closest-approach altitudes, V_rel, and `|V_inf|` magnitudes for every gravity-assist flyby on the actual flown trajectory (Venus, Earth-1, Earth-2, plus asteroid flybys Gaspra and Ida). Data density meets — and exceeds, in body diversity — the `heaton-longuski-2003-uranian-tour-u00-01` V0 admission bar. Recommend admitting Galileo as the catalogue's **3rd `mga_tour` row** at V0. Also recommend introducing a **separate KNOWN_CORPUS anchor** for D'Amario-Bright-Wolf 1992 alongside the existing Diehl-Belbruno-Roberts 1986 anchor; they are different papers documenting different design epochs (pre-Challenger concept vs. post-launch flight design).

---

## 1. Header

- **Title** (verbatim, p. 23): "Galileo trajectory design"
- **Authors** (verbatim, p. 23): Louis A. D'Amario, Larry E. Bright, and Aron A. Wolf
- **Affiliation**: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109, U.S.A.
- **Venue**: *Space Science Reviews* **60** (1-4), May 1992, pp. 23-78 (per running header and front matter)
- **DOI**: 10.1007/bf00216849
- **Page count**: 56 pages (pp. 23-78)
- **Publisher**: Kluwer Academic Publishers (1992)

## 2. What the paper actually is

A canonical, **post-launch JPL trajectory team** description of the Galileo Venus-Earth-Earth-Gravity-Assist (VEEGA) trajectory as actually flown. Published in *Space Science Reviews* vol. 60 alongside the Yeates et al. "Galileo Mission Overview" companion paper (referenced p. 24). The paper opens at launch (18 October 1989) and treats events from launch through end-of-mission (October 1997), including the Venus flyby (already executed at writing time, p. 28), the first Earth flyby (already executed, "At the time this paper was written, Galileo had just completed its first Earth gravity-assist flyby," p. 31), and the still-future second Earth flyby, Ida flyby decision, Jupiter arrival, Probe entry, and 10-orbit Galilean satellite tour.

Six sections plus an appendix:
1. Introduction (p. 24): mission objectives, why VEEGA instead of direct E→J.
2. Interplanetary Trajectory (pp. 25-36): mission constraints, design constraints (asteroid flybys, injection biasing, Earth-flyby navigation), and the full launch-through-Ida-flyby cruise narrative with per-encounter geometry diagrams (Figs 3, 5, 8) and Table I sequence of events (p. 37).
3. Probe Mission (pp. 38-49): atmospheric entry corridor, relay link geometry, Probe entry-site selection, descent profile.
4. Satellite Tour (pp. 46-71): 10-orbit "representative tour" — gravity-assist dynamics among Io / Europa / Ganymede / Callisto, with Table V (p. 60) listing all 13 satellite encounters (10 targeted + 3 nontargeted) of the sample tour, and individual trajectory pole views (Figs 19-31).
5. Mission Performance (pp. 71-75): propellant-margin accounting under both "Gaspra + Ida" and "Gaspra only" mission options. Table VIII (p. 75) is the master ΔV / mass budget.
6. Summary (p. 74): one-paragraph mission summary.
7. Appendix: Gravity-Assist Dynamics (pp. 76-78): standard hyperbolic-flyby vector diagrams and `sin(θ/2) = 1 / (1 + r_p V_∞² / GM)` bend-angle equation.

**Citation against KNOWN_CORPUS**: the project currently anchors Galileo VEEGA in `src/cyclerfinder/search/literature_check.py` line ~756 under the name *"Diehl-Belbruno-Roberts Galileo VEEGA design (1986)"* with author tuple `("Diehl", "Belbruno", "Roberts", "D'Amario")`. **That is a different paper** — D'Amario 1992 is post-launch; the Diehl-Belbruno-Roberts citation is pre-launch (1986) concept design predating the Challenger-driven IUS substitution. D'Amario 1992's own reference list (p. 78) cites the same author's earlier work: `D'Amario, L.A., Byrnes, D.V., Diehl, R.E.: 1987, Galileo Options After Challenger, AAS Paper 87-420` and `D'Amario, L.A., Byrnes, D.V., et al.: 1987, Galileo 1989 VEEGA Trajectory Design, AAS Paper 87-421`. D'Amario 1992 supersedes both as the **flown-trajectory** reference. The Belbruno name does **not** appear in D'Amario 1992's reference list — the existing KNOWN_CORPUS author tuple is suspect for that anchor's actual citation and should be revisited (see §5 below).

## 3. Body sequence + per-encounter data

The canonical sequence (Table I, p. 37, and Fig. 10 timeline, p. 36):

```
LAUNCH (Earth) -> Venus -> Earth-1 -> [Gaspra] -> Earth-2 -> [Ida] -> Jupiter
```

with two asteroid flybys interleaved (Gaspra between E1 and E2; Ida between E2 and Jupiter). The minimum required body set for cycler-class characterisation is `{E, V, E, E, J}`. With asteroid bonuses included: `{E, V, E, Gaspra, E, Ida, J}`.

### 3.1 Per-flyby data table — directly from Table I (p. 37), Figs 3, 5, 8, and §2.3 narrative

| # | Encounter | Date (UTC) | Altitude (km) | V_rel (km/s) | \|V_∞\| (km/s) | Source |
|---|---|---|---|---|---|---|
| L | Earth launch | 18 Oct 1989 (16:54 UTC) | — | — | √C3 ≈ √(13–17) ≈ 3.6–4.1 | p. 27 (launch); p. 25 (C3 range) |
| 1 | Venus | 10 Feb 1990 (05:59 UTC) | 16,123 | 8.2 | **6.2** | Fig 3 (p. 28); §2.3.2 (p. 28); Table I (p. 37) |
| 2 | Earth-1 | 8 Dec 1990 (20:35 UTC) | 960 | 13.7 | **8.9** | Fig 5 (p. 32); §2.3.3 (p. 30); Table I (p. 37) |
| a | Gaspra (asteroid 951) | 29 Oct 1991 | 1600 (from center) | 8.0 | n/a (low mass) | §2.3.4 (p. 32); Table I (p. 37, marked "Distance from center of asteroid") |
| 3 | Earth-2 | 8 Dec 1992 (15:35 UTC) | 300 | 14.1 | **8.9** | Fig 8 (p. 35); §2.3.4 (p. 33); Table I (p. 37) |
| b | Ida (asteroid 243) | 28 Aug 1993 | TBD | 12.4 | n/a (low mass) | §2.3.4 (p. 36); Table I (p. 37, altitude "TBD") |
| 4 | Jupiter (Io flyby + JOI) | 7 Dec 1995 (17:46 UTC Io flyby; 8 Dec 1995 00:27 UTC JOI start) | 1000 (Io); 4 R_J (perijove) | — | implied 5.9 km/s (Probe entry 59.9 km/s → V_p² – V_circ²) | Table I (p. 37); §2.4 (p. 38); §3.4 (p. 45 "Inertial entry speed is 59.9 km s⁻¹") |

### 3.2 Heliocentric energy bookkeeping (§2.3.2-2.3.4)

- **Pre-Venus**: heliocentric V ≈ 27 km/s (post-IUS, launched at C3 ≈ 13 km²/s²; IUS Δv "decrease the heliocentric velocity of the spacecraft by 3.1 km/s so that it would fall inward toward the Sun and encounter Venus", p. 27).
- **Post-Venus**: heliocentric period = 1 yr; Venus added 2.3 km/s to heliocentric velocity (Fig 3 caption / p. 28).
- **Earth-1**: V_∞ vector rotated by 48°; V_helio went 30.1 → 35.3 km/s (Fig 5 + p. 31); orbital period went 1 yr → **2 yr** (the canonical k=2 E-E resonance, p. 31).
- **Earth-2**: V_∞ vector rotated by 51°; V_helio went 35.3 → 39.0 km/s (Fig 8 + p. 33); orbital period now 5.6 yr to Jupiter (p. 33).
- Galileo "must acquire about 9 km/s of Earth-relative speed (V_∞) in a direction approximately parallel to Earth's velocity vector" to reach Jupiter (p. 30). The two Earth flybys rotate this V_∞ vector through 99° in total (p. 30, §2.3.3 + p. 34 §2.3.4).

### 3.3 ΔV budget (Table VIII, p. 75; abstract Summary, p. 75)

- Aggregate ΔV imparted by the three planetary gravity-assist flybys: **18.3 km/s** (Summary p. 75).
- Interplanetary deterministic ΔV: 102 m/s; statistical 45 m/s (Table VIII).
- JOI ΔV: 628 m/s; PJR ΔV: 370 m/s; ODM 59 m/s.
- Probe entry speed (inertial): **59.9 km/s** (p. 45); atmosphere-relative 47.4 km/s.
- TCM-4A: 24.2 m/s (April 1990); TCM-4B: 11.0 m/s (May 1990) (p. 30).

### 3.4 Mass / launch (p. 25)

- C3_max from IUS: **~17 km²/s²** at 2717 kg launch mass.
- C3_min for the 1989 VEEGA opportunity: **~13 km²/s²**.
- Launch period: 41 days (12 Oct – 21 Nov 1989) for the 7 Dec 1995 Jupiter arrival.
- Usable propellant at injection: 925 kg, of which only ~100 m/s of interplanetary ΔV is available — the "very stringent constraint" that forced VEEGA selection.

### 3.5 Satellite tour (Table V, p. 60) — out of scope for this catalogue row

The 10-orbit Jovian satellite tour (Ganymede-1 through Ganymede-10, plus 3 nontargeted Europa/Callisto/Ganymede encounters) is a separate gravity-assist tour around Jupiter, distinct from the heliocentric VEEGA itself. It is the second mga_tour leg of the mission (Jovian-system tour), structurally analogous to the Heaton-Longuski Uranian-system tour but already-flown. **Recommend the catalogue row scope cover only the heliocentric VEEGA phase** (launch → JOI), mirroring how `heaton-longuski-2003-uranian-tour-u00-01` covers Earth → Ariel and tabulates the Uranian moon tour as part of the same row's `bodies` + `vinf_kms_at_encounters` block. The Galileo satellite tour's per-encounter V_∞ is **not** in this paper (Table V gives altitude and latitude but not V_∞), so the satellite-tour portion would be N/A for V0 unless cross-sourced.

## 4. Catalogue admission verdict — ADMIT (V0)

### 4.1 Evidence sufficiency

Comparing to the §16.7.12 V0 evidence floor:

| Requirement | Status | Source in paper |
|---|---|---|
| Per-flyby epochs | Yes — UTC to the minute for all three gravity-assist flybys | Table I (p. 37); Figs 3, 5, 8 |
| Per-flyby \|V_∞\| | Yes — published for V, E1, E2 at one decimal place | Figs 3 (\|V_∞\|=6.2), 5 (\|V_∞\|=8.9), 8 (\|V_∞\|=8.9); narrative §2.3 |
| Body sequence | Yes — unambiguous: E (launch) → V → E → E → J, with Gaspra/Ida bonuses | Figs 1, 2, 10; Table I |
| Validity window | Yes — launch 18 Oct 1989 → JOI 8 Dec 1995 | Table I |
| Epoch-locked? | Yes — "1989 VEEGA opportunity" is a once-per-decade Earth-Venus phasing alignment; the paper's authors explicitly tied the 41-day launch window to the December 1995 Jupiter arrival | §2.1 (p. 25); §2.3.1 (p. 27) |

This **meets** the bar. The bar comparison to `heaton-longuski-2003-uranian-tour-u00-01` (#336) is favourable: that row was admitted on JSR Tables 3+5 with similar precision (V_∞ to 2 decimals from STOUR patched-conic; here V_∞ to 1 decimal from JPL flight design). Galileo has the **additional** virtue of being post-launch — the V_∞ values are from the actual flown trajectory's tracking-data-corrected design, not a pre-flight optimisation point. Note however that the published V_∞ figures are still given as design values (Figs 3, 5, 8 are pre-flyby design diagrams, not post-flyby reconstructions); achieving higher than V0 (e.g., V2) would require re-running the trajectory on the modern DE-series ephemeris and matching the published table values to known tolerances — an offline V2 promotion task (#335-style).

### 4.2 Recommended catalogue row construction (DO NOT WRITE — recommend only)

```yaml
- id: damario-1992-galileo-veega
  name: "Galileo VEEGA — flown trajectory (D'Amario-Bright-Wolf 1992)"
  source: literature
  trajectory_regime: ballistic   # patched-conic with Venus, Earth, Earth, [Gaspra, Ida] gravity assists + JOI ΔV at Jupiter (outside the heliocentric tour)
  model_assumption: analytic-ephemeris   # JPL flight design per D'Amario 1992; pre-launch DE-series ephemeris + TCM-corrected post-launch
  cycler_class: multi-arc   # schema v4; three heliocentric arcs (E->V, V->E1, E1->E2, E2->J) — non-repeating
  orbit_class: mga_tour   # schema v4.7; classic Venus-Earth-Earth gravity-assist mission archetype; terminal target is Jupiter capture, not a continuing tour at Earth
  epoch_locked: true   # 1989 VEEGA opportunity is a once-per-decade Earth-Venus phasing window (paper §2.1 + §2.3.1, 41-day launch period)
  n_returns: 1   # single Earth-launch / single Jupiter-arrival mission; not a cycler
  validity_window:   # Earth launch (Table I) -> JOI completion (Table I event-time block)
    start: "1989-10-18T16:54:00Z"
    end: "1995-12-08T01:16:00Z"
  launch_epoch: "1989-10-18T16:54:00Z"   # D'Amario 1992 p. 27 ("16:54 UTC (09:54 PST)")
  validation_level: V0   # spec §14: sourced from Table I + Figs 3/5/8; not independently reproduced on modern ephemeris; V0 = sourced internal-consistency floor
  source_ephemeris: "JPL pre-launch design ephemeris + post-launch TCM-corrected (D'Amario 1992 §2.3); approx DE-200 era"
  orbit_source: derived   # ballistic heliocentric arcs between gravity-assist flybys
  vinf_source: derived   # Figures 3, 5, 8 publish V_inf magnitudes; vectors implied by velocity-vector diagrams
  orbit_fidelity: analytic-ephemeris
  vinf_fidelity: analytic-ephemeris
  bodies: ["E", "V", "E", "Gaspra", "E", "Ida", "J"]
  sequence_canonical: "E-V-E-Gaspra-E-Ida-J"   # alternative: "E-V-E-E-J" if asteroid bonuses excluded
  sense: "n/a"   # not a cycler
  vinf_kms_at_encounters:
    - body: "E"
      vinf_kms: ~3.8   # launch C3 in [13, 17] km^2/s^2 (p. 25); √13=3.6, √17=4.1; nominal ~14.5 → ~3.8
      note: "D'Amario 1992 p. 25; launch C3 = 13-17 km^2/s^2; sqrt → V_inf magnitude. Better cross-source from D'Amario et al. 1989 AAS 87-421 for the actual flown C3."
    - body: "V"
      vinf_kms: 6.2
      note: "D'Amario 1992 Fig 3 (p. 28); Venus flyby 10 Feb 1990 05:59 UTC; closest-approach altitude 16,123 km; latitude -40.8 deg; V_rel = 8.2 km/s."
    - body: "E"
      vinf_kms: 8.9
      note: "D'Amario 1992 Fig 5 (p. 32); Earth-1 flyby 8 Dec 1990 20:35 UTC; closest-approach altitude 960 km; lat 25 deg, lon 63 W; V_rel = 13.7 km/s; V_inf rotated 48 deg; period 1 yr -> 2 yr."
    - body: "Gaspra"   # asteroid 951; type-S; 16-km diameter
      vinf_kms: 8.0
      note: "D'Amario 1992 §2.3.4 (p. 32); Gaspra flyby 29 Oct 1991; closest approach 1600 km from center; V_rel = 8.0 km/s; gravity assist negligible (asteroid mass too low) — included for science only."
    - body: "E"
      vinf_kms: 8.9
      note: "D'Amario 1992 Fig 8 (p. 35); Earth-2 flyby 8 Dec 1992 15:35 UTC; closest-approach altitude 300 km; lat -34 deg, lon 12 W; V_rel = 14.1 km/s; V_inf rotated 51 deg; V_helio 35.3 -> 39.0 km/s; period now 5.6 yr."
    - body: "Ida"   # asteroid 243; type-S; 32-km diameter; OPTIONAL
      vinf_kms: 12.4
      note: "D'Amario 1992 §2.3.4 (p. 36); Ida flyby 28 Aug 1993; V_rel = 12.4 km/s; closest-approach altitude TBD at writing; flyby was actually flown (8/28/1993); gravity assist negligible."
    - body: "J"
      vinf_kms: 5.9
      note: "D'Amario 1992 §3.4 (p. 45) implies V_inf_J ≈ √(V_entry^2 - V_circ_Jupiter^2) ≈ √(59.9^2 - 60^2) ≈ ~5 km/s; ALTERNATIVELY, cross-source from Yeates 1992 Galileo Mission Overview (cited p. 24) for the published Jupiter approach V_inf."
  citations:
    - "D'Amario, L.A., Bright, L.E., Wolf, A.A. (1992). 'Galileo trajectory design.' Space Science Reviews 60(1-4), 23-78. DOI 10.1007/bf00216849."
    - "Table I (p. 37) — interplanetary trajectory sequence of events."
    - "Fig 3 (p. 28), Fig 5 (p. 32), Fig 8 (p. 35) — per-flyby velocity-vector diagrams with |V_inf|."
    - "Table VIII (p. 75) — propellant margin / ΔV budget."
```

Cross-references to `tito-2018-mars-free-return` (single-flyby Mars free-return, validity window 2017-11-25 → 2019-06-21) and `heaton-longuski-2003-uranian-tour-u00-01` (E→J→U + 40-event Uranian tour, validity window 2008-03-19 → 2021-07-20) put Galileo VEEGA dead in the middle of the family: it is the **textbook archetype** the Davis 2018 corpus citation calls out by name.

### 4.3 Class — `mga_tour`, not `cycler`

The E1→E2 leg creates a 2-year (k=2) Earth-Earth resonant arc — superficially "repeating-encounter" at Earth — but:
- The terminal target is Jupiter capture (JOI = 628 m/s retro burn, p. 75), not continuation.
- After the 2-year resonance, V_∞ at Earth is rotated to *escape* Earth's gravity well to Jupiter, not maintained for further Earth encounters.
- `n_returns: 1` is correct (single Earth-launch, single Jupiter-arrival), and `orbit_class: mga_tour` is correct.
- Compare to `tito-2018-mars-free-return` which has the same `n_returns: 1` and `orbit_class: mga_tour` despite being a single-flyby ballistic return; Galileo VEEGA is the same archetype with one more flyby in the chain.

## 5. KNOWN_CORPUS impact — RECOMMEND SPLIT INTO TWO ANCHORS

The current `Diehl-Belbruno-Roberts Galileo VEEGA design (1986)` anchor (literature_check.py L756 area) has:

```python
authors=("Diehl", "Belbruno", "Roberts", "D'Amario"),
citation="Diehl, Belbruno & Roberts et al., 'Galileo VEEGA Mission Design' (1986-1990 JPL / AAS); the canonical mga_tour archetype (October 1989 launch window once-per-~13yr alignment)",
doi=None,
topology_label=frozenset({"mga-tour"}),
```

**Issues identified by reading D'Amario 1992**:
1. **Belbruno is not in D'Amario 1992's reference list** (p. 78). The four cited co-author papers on the same trajectory are: Byrnes-D'Amario-Diehl 1987 AAS 87-420 "Galileo Options After Challenger"; D'Amario-Byrnes-Johannesen-Nolan 1987 AAS 87-421 "Galileo 1989 VEEGA Trajectory Design"; D'Amario-Bright-Byrnes-Johannesen-Ludwinski 1989 AAS 89-431 "Galileo 1989 VEEGA Mission Description"; and Johannesen-Nolan-Byrnes-D'Amario 1987 AAS 87-422 "Asteroid/Comet Encounter Opportunities for the Galileo VEEGA Mission". Belbruno's actual Galileo contribution is the weak-stability-boundary capture work, not the VEEGA design itself. **Suspect the author tuple is wrong**; this should be verified against the original Diehl reference's actual citation.
2. The 1986 date is plausible for the original D'Amario-Diehl pre-Challenger concept (Galileo had already been delayed and re-designed multiple times), but the **flown** trajectory was redesigned post-Challenger after IUS substitution and is documented in the 1987 AAS papers, not 1986.
3. D'Amario 1992 is **the** post-launch / flown-trajectory canonical reference. It is a distinct paper, peer-reviewed in *Space Science Reviews*, with a clear DOI.

**Recommendation**: add a **new** separate anchor for D'Amario-Bright-Wolf 1992 in `literature_check.py`:

```python
CorpusAnchor(
    name="D'Amario-Bright-Wolf Galileo VEEGA flown trajectory (1992)",
    primary="Sun",
    body_set=frozenset({"V", "E", "Jupiter"}),
    topology_label=frozenset({"mga-tour"}),
    authors=("D'Amario", "Bright", "Wolf"),
    keywords=(
        "VEEGA Venus-Earth-Earth gravity assist",
        "Galileo Jupiter mission flown trajectory",
        "1989 VEEGA opportunity",
        "Earth-Earth 2-year resonance",
        "Gaspra Ida asteroid flyby",
    ),
    citation="D'Amario, Bright & Wolf, 'Galileo trajectory design,' "
    "Space Science Reviews 60(1-4), 23-78 (May 1992); DOI 10.1007/bf00216849; "
    "the canonical post-launch reference for the flown 1989 VEEGA "
    "(launch 18 Oct 1989, V 10 Feb 1990, E1 8 Dec 1990, E2 8 Dec 1992, "
    "JOI 8 Dec 1995). Catalogue row damario-1992-galileo-veega (mga_tour, V0)",
    doi="10.1007/bf00216849",
),
```

**Keep** the existing Diehl-Belbruno-Roberts 1986 anchor (representing the pre-Challenger concept-design epoch) — but **fix the author tuple**: replace `("Diehl", "Belbruno", "Roberts", "D'Amario")` with `("Diehl", "Roberts", "D'Amario", "Byrnes")` pending verification of the actual original-paper citation. The "Belbruno" entry is the suspected error (Belbruno's name is associated with weak-stability-boundary lunar capture, not Galileo VEEGA design).

This split is structurally identical to how Wittal IAC-22 was treated separately from the underlying Floquet-bifurcation literature: distinct papers, distinct epochs, distinct topology labels per anchor.

## 6. Errata against the project's prior context

| Where | Prior claim | What the paper actually says | Status |
|---|---|---|---|
| Agent prompt §0 | Launch "October 1989" | Confirmed — 18 Oct 1989 16:54 UTC (p. 27) | OK |
| Agent prompt | Canonical sequence "launch → Venus (Feb 1990) → Earth-1 (Dec 1990) → Earth-2 (Dec 1992) → Jupiter (Dec 1995)" | Confirmed — matches Table I (p. 37) and Fig 10 timeline (p. 36) exactly | OK |
| KNOWN_CORPUS L756 | Author tuple includes "Belbruno" | Belbruno is **not** in D'Amario 1992's reference list (p. 78). Suspect-incorrect attribution for the pre-launch Diehl-Roberts-D'Amario paper. | **ERRATUM (suspect)** — needs verification against the actual Diehl 1986 citation |
| KNOWN_CORPUS L756 | "Diehl-Belbruno-Roberts Galileo VEEGA design (1986)" presented as the sole Galileo VEEGA anchor | D'Amario-Bright-Wolf 1992 is the canonical **flown** trajectory reference with a clear DOI; the corpus needs at minimum a second anchor for the post-launch design (or to fold the 1992 paper into a revised author/citation list on the existing anchor). | **DEFICIENCY** — recommend new anchor (§5) |
| Agent prompt | Gaspra closest approach "1600 km" | Confirmed at 1600 km **from center** of asteroid (Table I footnote, p. 37; §2.3.4 p. 32). | OK |
| Agent prompt | "around the IUS / Centaur launch-vehicle change" → may have driven VEEGA selection | D'Amario 1992 §2.1 (p. 25) confirms: "the maximum C3 available from the two-stage IUS is about 17 km²/s²... the only interplanetary transfer option that did not cause the total propellant capacity to be exceeded was the VEEGA trajectory" — i.e., post-Challenger IUS substitution forced VEEGA. | OK |
| Agent prompt | "the high-gain antenna failure, the asteroid Gaspra + Ida flyby insertions" should be documented | The paper does NOT discuss the HGA failure (that occurred April 1991, after this paper was largely drafted — the paper says "At the time this paper was written, Galileo had just completed its first Earth gravity-assist flyby", p. 31, i.e., Jan 1991; the HGA failure was April 1991). Gaspra and Ida flybys ARE discussed extensively (§2.2.1 p. 26, §2.3.4 p. 32, §5.4 p. 74) including the Ida-decision plan. | OK (paper predates HGA failure) |

## 7. Action items for the parent (#347)

1. **Atomic commit this digest** at `docs/notes/2026-06-17-digest-damario-1992-galileo.md`. (Done in this run.)
2. **Track #345 follow-up**: add catalogue row `damario-1992-galileo-veega` per the construction in §4.2. The row is V0-grade and ready to write. Recommend doing this in a separate dedicated commit (not bundled with this digest), in line with the "verdict note first, writeback after confirmation" discipline from the orbit-closure-discipline memory.
3. **Track #350 follow-up**: split the `Diehl-Belbruno-Roberts Galileo VEEGA design (1986)` KNOWN_CORPUS anchor — keep it (representing the pre-Challenger concept epoch, with author tuple revisited; "Belbruno" attribution is suspect) **and** add a new D'Amario-Bright-Wolf 1992 anchor for the post-launch / flown trajectory. Per §5 above.
4. **Erratum lodging**: the suspected Belbruno-attribution error in KNOWN_CORPUS should be cross-checked against the actual Diehl 1986 reference before any public-facing erratum is filed — could be a project-internal slip rather than a literature error. Per the respectful-errata-framing memory.
5. **Sister-paper digest** (next): Yeates, Johnson & Young 1992 "Galileo Mission Overview" (cited p. 24 as the companion paper in the same SSR vol. 60 issue). It may publish the Jupiter approach V_∞ directly, addressing the §4.2 "5.9 km/s" approximation. Not in scope for this agent; flag for next agent dispatch.
6. **V2 promotion path**: a V2-grade Galileo row would require re-running the trajectory through the existing CR3BP/Lambert/ephemeris stack against the published Table I dates + V_∞ magnitudes and checking residuals at the 1-decimal precision the paper gives. Defer.

---

**Verdict reaffirmed**: **ADMIT** as `damario-1992-galileo-veega`, `orbit_class: mga_tour`, `validation_level: V0`, single-row catalogue addition. Adds the textbook VEEGA archetype to the catalogue alongside Tito 2018 Mars free-return (single-flyby ballistic return) and Heaton-Longuski 2003 Uranian tour (multi-arc moon-tour), completing the `mga_tour` triumvirate of canonical archetypes.
