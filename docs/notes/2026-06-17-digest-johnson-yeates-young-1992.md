# Digest — Johnson, Yeates & Young 1992 "Space Science Reviews Volume on Galileo Mission Overview"

**Date**: 2026-06-17 AET (digest of `#358` — sister paper to D'Amario-Bright-Wolf 1992; parent #347 / track #345 classic-mission `mga_tour` admissions)
**Verdict (TL;DR)**: **HONEST NEGATIVE — does NOT materially extend D'Amario 1992's trajectory data.** This is a high-level science-objectives overview, not a trajectory design paper. Per-flyby V_∞ values, JOI orbital elements, and Jupiter approach V_∞ are **NOT published here**. The paper's mission-events table (Table I, p. 4) is a coarse year/month milestone list. The §4.2 "5.9 km/s" Jupiter approach V_∞ approximation in the #354 D'Amario digest is therefore **NOT confirmable from this companion paper**; cross-sourcing remains pending. **Two minor cross-checks** are positive: (1) Venus flyby date "February 9, 1990" (p. 4) is **consistent** with D'Amario's "10 Feb 1990 05:59 UTC" once timezone is accounted for (Feb 9 21:59 PST = Feb 10 05:59 UTC); (2) Probe entry velocity "48 km/s relative to the atmosphere" (p. 17) is **consistent** with D'Amario's 47.4 km/s atmosphere-relative (one-decimal vs zero-decimal precision). The current #356 Galileo catalogue admission verdict stands as-is.

---

## 1. Header

- **Title** (verbatim, p. 3): "SPACE SCIENCE REVIEWS VOLUME ON GALILEO MISSION OVERVIEW"
- **Authors** (verbatim, p. 3): T. V. Johnson, C. M. Yeates (Jet Propulsion Laboratory, Caltech), and R. Young (NASA Ames Research Center)
- **Venue**: *Space Science Reviews* **60**, 3–21, 1992 (page footer p. 3)
- **DOI**: not on the article face; the SSR vol. 60 DOI root is shared with D'Amario 1992
- **Page count**: 19 pages (pp. 3–21)
- **Publisher**: Kluwer Academic Publishers (1992); printed in Belgium

## 2. What the paper actually is

A **science-objectives overview**. Per the abstract (p. 3, verbatim): "The Galileo Mission is an extremely complex undertaking. This paper provides a brief historical overview, a discussion of broad scientific objectives, and a description of the spacecraft and trajectory characteristics."

Five sections plus front matter and acknowledgements:
1. Introduction (pp. 3–5): mission history (Pioneer/Voyager predecessor lineage, J. Van Allen working group recommendation, Jupiter Orbiter Probe rename, Challenger-driven IUS substitution forcing VEEGA), launch milestone narrative.
2. Objectives (pp. 5–7): why Jupiter (composition, atmosphere, satellite system, magnetosphere).
3. Scientific Objectives (pp. 7–11): Table II science objectives list; per-discipline narrative (Atmosphere, Satellites, Magnetosphere).
4. Galileo Spacecraft (pp. 11–19): Orbiter design (Telecommunications, Propulsion, Command and Data System, AACS, Probe Deceleration / Descent / Communications / Power / Command / Mass / Relay Radio Hardware).
5. Mission Design (pp. 19–20): VEEGA concept narrative + "10 cushion billiard shot" / "100 m/s" tour ΔV claim.

**It is the SISTER paper to D'Amario-Bright-Wolf 1992** ("Galileo trajectory design," same SSR vol. 60 issue, pp. 23–78, DOI 10.1007/bf00216849). The two papers were edited together by C. Russell (Acknowledgements, p. 21). The trajectory paper is the technical companion; this paper is the science companion.

**Citation against KNOWN_CORPUS**: Like D'Amario 1992, this paper does NOT cite Belbruno. The Bibliography (p. 21) is **intentionally sparse** — verbatim: "Due to the extremely broad nature of the science and engineering background discussed in this paper, a detailed reference list would be far too massive to be included easily." Lists 12 references total. Of relevance: **Yeates et al. 1985 "Galileo: Exploration of Jupiter's System" NASA SP-479** and **Landano-Jones 1983 AIAA 83-0097** and **O'Neil 1990 AIAA 90-3854**. None of these is by Belbruno; none is by Diehl-Roberts in the 1986 timeframe. **Cannot confirm or refute** Agent D's hypothesis about the "Diehl-Belbruno-Roberts 1986" KNOWN_CORPUS attribution from this paper alone.

## 3. What's in the paper — trajectory-data inventory

### 3.1 Table I (p. 4) — Major mission events

The full table reads (verbatim transcription):

| Year | Months | Event |
|---|---|---|
| 1989 | October | Launch |
| 1990 | February | Venus encounter |
| 1990 | December | Earth encounter |
| 1991 | October | Gaspra encounter |
| 1992 | December | Earth encounter |
| 1993 | August | Ida encounter |
| 1995 | August | Probe release |
| 1995 | December | Io encounter / Probe relay / Jupiter orbit insertion |
| 1996 | July | First Ganymede encounter |
| 1996 | October | Second Ganymede encounter |
| 1997 | October | Completion of 10 targeted encounter tour and tail petal |

**No epoch precision beyond year/month. No altitudes. No V_rel. No V_∞.** D'Amario 1992's Table I (p. 37) provides all of these to UTC-minute and one-decimal precision.

### 3.2 Body-text trajectory references

The paper's body text references the trajectory at narrative resolution only:

- p. 3–4 (Introduction): "the Galileo spacecraft and IUS, was successfully launched from the Kennedy Space Center... heading initially toward a flyby of Venus. The Venus flyby was accomplished on **February 9, 1990** and the first Earth flyby on **December 8, 1990**." — note the Feb 9 vs D'Amario's Feb 10 — see §5 below.
- p. 17 (Probe Deceleration Module): "Entry velocity relative to the atmosphere is **48 km s⁻¹**, far higher than any atmospheric entry attempted to date." — vs D'Amario p. 45 atmosphere-relative 47.4 km/s.
- p. 17: "heat shield... subjected to a hot atmospheric shock layer (**14000 K**)... heat transfer at the nose of the vehicle at peak heating exceeds **42 kW cm⁻²**... approximate mass of the forebody heat shield is **145 kg**, of which about 90 kg is expected to be lost by ablation during entry."
- p. 19 (Mass): "total mass of the Probe is about **331 kg**, that of the Probe adapters is **6.8 kg**, and that of the RRH equipments is **23.2 kg**."
- p. 19–20 (Mission Design): "Once captured by Jupiter's gravity, the Orbiter would repeat its **initial 200-day orbit** if nothing were done... For Galileo to be utilized more effectively during its limited lifetime, the orbital period must be shortened and the spacecraft targeted to make very close flybys of the Galilean satellites."
- p. 20: tour ΔV claim — verbatim: "**the entire 'satellite tour' can be flown so that rockets need supply only about 100 m s⁻¹ of velocity change** — 60 times less than what would be needed without the satellites' help!" Consistent with D'Amario Table VIII per-encounter ΔTCM allocations (#354 §3.3).

### 3.3 Figure 1 (p. 4) — VEEGA trajectory diagram

A schematic of the heliocentric VEEGA path with the same body sequence as Fig. 10 of D'Amario 1992 (Launch 10/18/89, Venus 2/10/90, E1 12/08/90, Gaspra 10/29/91, E2 12/08/92, Ida 8/28/93, Probe Release ~7/12/95, Jupiter Arrival 12/7/95, end 11/7/97). The diagram caption gives time ticks of "S/C = 30 days, EARTH = 30 days, VENUS = 30 days, JUPITER = 100 days." **No quantitative trajectory data on the figure itself.**

### 3.4 Probe + Spacecraft mass / power (Table III, p. 15; §4 narrative)

- Orbiter mass: **2380 kg** (includes 103-kg science payload), HGA 4.8 m, peak data rate 134 kbps (p. 11–12 and Fig. 2 caption p. 13).
- Probe mass: **335 kg** (Fig. 2 p. 13 caption; includes 28-kg avionics); Descent module 118 kg.
- Probe descent battery: Li-SO₂; **21 A·hr** capacity (about 730 W·hr), required 16.3 A·hr (p. 18). "Battery lifetime 55–60 min" (Fig. 2 p. 13).
- Probe diameter: 125 cm; height 86 cm (Fig. 2 caption p. 13).
- Probe receiver acquisition probability: **0.995** within 50 s; false-lock probability less than **5 × 10⁻⁵** (p. 19).

These are spacecraft-engineering figures, not trajectory parameters. No new V0-grade evidence.

## 4. Does this paper extend D'Amario 1992's per-encounter data?

**No.** Direct comparison:

| Quantity | D'Amario 1992 | Johnson-Yeates-Young 1992 | Extends? |
|---|---|---|---|
| Per-flyby epoch (UTC) | yes — to the minute | year/month only (Table I) | No, less precise |
| Per-flyby \|V_∞\| | yes (Figs 3/5/8) | not published | No |
| Per-flyby altitude | yes (Table I) | not published | No |
| Per-flyby V_rel | yes (Figs 3/5/8) | not published | No |
| Per-flyby V_∞ rotation angle | yes (Figs 5/8: 48°, 51°) | not published | No |
| Pre-/post-flyby V_helio | yes (Figs) | not published | No |
| Probe inertial entry speed | yes (p. 45: **59.9 km/s**) | not published — only atmosphere-relative 48 km/s | No |
| Probe atmosphere-relative entry | yes (p. 45: **47.4 km/s**) | yes — "48 km s⁻¹" (p. 17), 1-digit | Consistent, not extending |
| Jupiter approach \|V_∞\| | implied via 59.9 km/s minus circular-Jupiter speed (~5–6 km/s) | not published | No |
| JOI orbital elements | yes — JOI ΔV 628 m/s, perijove 4 R_J, period 200-day | partial: "initial 200-day orbit" (p. 20) | Confirms 200-day period; no other elements |
| JOI ΔV | yes — 628 m/s (Table VIII) | not published | No |
| C3_launch | yes — 13–17 km²/s² | not published | No |
| Tour ΔV total | yes (Table VIII per-stage) | yes — "100 m/s" total (p. 20), confirms order of magnitude | Consistent, not extending |
| Ida flyby altitude | TBD at writing in D'Amario | not published | Neither extends |
| Gaspra closest approach | 1600 km from center | "1600 km" in Note Added in Proof (p. 21) | Consistent, not extending |

**Conclusion**: this paper publishes **no new quantitative trajectory data** that D'Amario 1992 doesn't already have, with the single exception of the **Note Added in Proof** (p. 21, verbatim): "On October 29, 1991 Galileo completed one of its first scientific milestones when it flew by the asteroid Gaspra at a range of **1600 km**. Preliminary data, including four images returned in early November 1991, indicate that this first ever spacecraft encounter with an asteroid was highly successful." — but D'Amario also gives 1600 km, so this is corroboration not extension.

## 5. Cross-check: Venus flyby date discrepancy resolved

Johnson-Yeates-Young (p. 4): "The Venus flyby was accomplished on **February 9, 1990**."
D'Amario 1992 (p. 28, Fig 3): "Venus, 10 Feb 1990 05:59 UTC."

These are the **same event**: 10 Feb 1990 05:59 UTC = 9 Feb 1990 21:59 PST. JPL routinely reports times in PST (local-pacific) in mission-overview prose while engineering trajectory documentation uses UTC. Both citations are correct; the Johnson paper drops the time portion entirely. **No erratum.**

## 6. Belbruno-attribution check

Agent D's #354 verdict (§5) flagged that the current `Diehl-Belbruno-Roberts Galileo VEEGA design (1986)` anchor in `KNOWN_CORPUS` (literature_check.py ~L756) has `authors=("Diehl", "Belbruno", "Roberts", "D'Amario")` and that "Belbruno" is suspect — D'Amario 1992's reference list (p. 78) does not cite Belbruno.

**This paper's reference list (p. 21) is intentionally sparse** (verbatim p. 21: "a detailed reference list would be far too massive"). It contains 12 references, none by Belbruno, none by Diehl-Roberts in the 1986 timeframe. **The Belbruno-attribution question is therefore NEITHER confirmed NOR refuted by this paper.** Agent D's hypothesis remains pending an actual lookup of the original 1986 Diehl-Roberts-D'Amario paper(s) cited in the project's KNOWN_CORPUS — the work the project's anchor purportedly references.

## 7. Catalogue impact — none

Does this paper enable a V0→V1 promotion of the `damario-1992-galileo-veega` row (per Spec §14, V1 = like-for-like same-model reproduction of published values)? **No.** V1 requires reproducing per-encounter V_∞ on a sourced ephemeris model — this paper publishes no such values. The reproducibility evidence remains D'Amario 1992 alone, and V0 is the appropriate ceiling for that row.

Does this paper provide additional KNOWN_CORPUS keywords beyond what Agent D recommended for a new D'Amario-Bright-Wolf 1992 anchor? Minor additions only: "Jupiter Orbiter Probe (JOP)" — the pre-1977 mission name (p. 3); "STS-34 Atlantis" (p. 3); "Project Galileo" (p. 3). None of these is a trajectory keyword. **No KNOWN_CORPUS changes are mandated by this paper.**

Does this paper introduce a separate catalogue admission (e.g., the Probe entry as a distinct entry?)? **No.** The Probe entry is a single-event atmospheric-entry capsule, not a multi-body trajectory; out of scope for the catalogue's four classes (cycler / quasi_cycler / precursor_mga / mga_tour).

## 8. Errata against the project's prior context

| Where | Prior claim | What the paper actually says | Status |
|---|---|---|---|
| Agent prompt (this digest) | "may publish the Jupiter approach V_∞ directly, addressing the §4.2 5.9 km/s approximation" | Paper does NOT publish Jupiter approach V_∞. The §4.2 5.9 km/s approximation remains the best D'Amario-derived value; needs Yeates et al. 1985 SP-479 or a JPL post-mission report for a sourced direct citation. | **Honest negative** — flagged in this digest |
| #354 verdict §4.2 (5.9 km/s Jupiter approach V_∞) | "ALTERNATIVELY, cross-source from Yeates 1992 Galileo Mission Overview" | The Yeates 1992 paper (this one) does not publish Jupiter approach V_∞. Cross-source elsewhere. | **Action item revised** — see §9 below |
| #354 verdict §6 errata table | Said "Note paper does not discuss HGA failure (drafted Jan 1991)" | This paper confirms — drafted similarly; Note Added in Proof p. 21 only adds the Gaspra-flyby completion (Oct 1991), still does not mention HGA failure (April 1991). Suggests both papers were typeset for the 1992 SSR vol. 60 issue in mid-1991 with a 4-month proof window. | **OK** — consistent |
| #354 verdict §2 | D'Amario 1992 sister paper "may publish the Jupiter approach V_∞" | Does not. | **Erratum** to the #354 hypothesis — revised below |
| Project precedent (D'Amario 1992 catalogue row §4.2) | `vinf_kms_at_encounters` for body J entry: "5.9 km/s — D'Amario 1992 §3.4 (p. 45) implies V_inf_J ≈ √(V_entry² - V_circ_Jupiter²); ALTERNATIVELY cross-source from Yeates 1992" | **The Yeates 1992 cross-source falls through.** Update the row's note to: "5.9 km/s (derived); D'Amario 1992 §3.4 implies V_inf_J ≈ √(59.9² - 60²) ≈ 5–6 km/s. Yeates et al. 1992 (sister SSR paper) does NOT publish Jupiter approach V_∞ — no cross-source available in the SSR vol. 60 issue. For a sourced Jupiter approach V_∞ cross-check, consult the JPL Galileo mission post-flight report (1996+) or D'Amario et al. 1989 AAS 89-431." | **Erratum to D'Amario row construction** — see §9 |

## 9. Action items for the parent (#347)

1. **Atomic commit this digest** at `docs/notes/2026-06-17-digest-johnson-yeates-young-1992.md`. (Done in this run.)
2. **#356 / #358 closure**: confirm to the parent that **Johnson-Yeates-Young 1992 does NOT materially extend D'Amario 1992**. The Galileo `mga_tour` catalogue admission verdict #354 stands as written. No V0→V1 promotion enabled by this digest.
3. **D'Amario row construction revision** (when #356 writeback happens): the `vinf_kms_at_encounters` Jupiter body J entry should be revised away from the Yeates-1992 cross-source aspiration. Use the wording recommended in §8 above. The "5.9 km/s" value remains an order-of-magnitude derivation, not a sourced value — flag it for V2-promotion attention if any.
4. **Belbruno-attribution check remains open**. Neither D'Amario 1992 nor Johnson-Yeates-Young 1992 cites Belbruno for the VEEGA design. The "Diehl-Belbruno-Roberts 1986" attribution in KNOWN_CORPUS still wants verification against the actual 1986 paper. Track as a separate to-do (probably handled in #364 KNOWN_CORPUS errata wave already in flight).
5. **Probe entry as Probe-mission row?** This paper documents the Probe atmospheric-entry mission alongside the Orbiter — atmosphere-relative entry 48 km/s, inertial 59.9 km/s (cross-source from D'Amario), 14000 K shock layer, 42 kW/cm² peak heat flux, 145 kg ablating heat shield. **Out of scope for this catalogue** (atmospheric entry isn't a trajectory class). Note as a citation candidate if the catalogue's scope ever extends to atmospheric-entry probes.
6. **Sister-paper digest follow-up**: Yeates et al. 1985 NASA SP-479 ("Galileo: Exploration of Jupiter's System") is the 1985 pre-launch reference cited p. 21. May (or may not) publish the Jupiter approach V_∞. If acquisitions can source SP-479 in PDF, queue a follow-on digest. **Recommendation**: not high priority; the catalogue admission verdict stands at V0 with current sourcing.

---

**Verdict reaffirmed**: **HONEST NEGATIVE**. Johnson-Yeates-Young 1992 is a science-objectives overview that does NOT publish per-flyby trajectory data, JOI orbital elements, or Jupiter approach V_∞ beyond what D'Amario 1992 already gives. The Galileo VEEGA `mga_tour` catalogue admission at V0 stands; the only revised action item is correcting the D'Amario row's `vinf_kms_at_encounters` Jupiter-entry note to drop the Yeates 1992 cross-source aspiration. No KNOWN_CORPUS changes. No catalogue rows added.
