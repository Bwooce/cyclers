# Mining Note: Wittal, Miaule & Asher (2022), IAC-22-C1.6.6

**Source:** Matthew M. Wittal, Stéphane Miaule, Benjamin W. Asher, "Earth-Moon Cycler
Mission Design for Lunar Logistics," 73rd International Astronautical Congress, Paris,
2022, paper IAC-22-C1.6.6 (NTRS 20220013595). 7 pp. (6 pp. body + references).

**Mined:** 2026-06-11 from the full PDF. This is the primary source of catalogue row
`wittal-2022-em-cycler-family`, which until now was ingested from the NTRS landing-page
abstract only.

**Verdict: USEFUL (catalogue enrichment only — the CR3BP backfill does NOT unblock).**
The full text contains **zero initial conditions** of any kind: no state vectors, no
rotating-frame ICs, no Jacobi constants, no nondimensional periods, no stability
indices/eigenvalues, no epochs. The orbits exist only as figures. The NO_SOURCED_IC
blocker on this row is therefore **confirmed by the full text** — it is a publication
gap in the primary source, not an access gap. What the paper does give: family taxonomy
and stability domains, encounter cadences, and two small numeric tables (launch masses,
ΔV requirements) — enough to enrich the citation-only row, not to reproduce any member.

---

## 1. IC-completeness verdict (mining lens a) — NEGATIVE, precisely

- **Numbers in the paper:** exactly two tables — Table 1 (launch-vehicle delivered
  masses, p. 3) and Table 2 (ΔV requirements, p. 5) — plus scattered scalars in the
  text (perilune radii, period 27 ± 3 d, inclinations, burn sizes). Trajectories appear
  only as rendered figures (Figs. 1–5) and per-maneuver bar/line charts (Figs. 6–9,
  values readable only approximately from axes).
- **No reproducible state anywhere:** no Cartesian or elemental state vector, no epoch,
  no Jacobi constant, no nd-period, no monodromy/stability data, in any frame.
- **Dynamical model: never stated.** The paper nowhere names CR3BP, an ephemeris model,
  a force model, or a simulation tool. Internal evidence points to an ephemeris-class
  numerical simulation rather than CR3BP: the 2-petal period is quoted as
  "approximately 27 ± 3 days" (Sec. 5, p. 4) — a CR3BP periodic orbit would have an
  exact period; maintenance ΔV accumulates per-flyby across 5 cycles with discrete MCC
  events (Figs. 6 and 8, pp. 5–6); Earth-flyby inclination and periapsis vary
  cycle-to-cycle (Figs. 7 and 9, pp. 5–6); Fig. 4 (p. 4) shows the same trajectories in
  the "Earth-Moon Rotating frame" and the "Earth inertial frame" as tool screenshots.
  But strictly: **model unstated** — record it that way.
- **Catalogue flag:** our row carries `model_assumption: cr3bp`. The full text does NOT
  support that tag (model unstated; evidence leans ephemeris-class). Recommend
  revisiting the tag (e.g., to a model-unstated annotation) at the next schema pass.

**Backfill consequence:** `wittal-2022-em-cycler-family` stays NO_SOURCED_IC, now
*permanent with respect to this source*. The only remaining source-axis lead is the
companion paper, ref [9] (see Sec. 5).

## 2. Family geometry, taxonomy, stability (mining lens b)

This is the same near-polar family set our catalogue row describes (see Sec. 4). The
families are "petal" cyclers in lunar resonance, explored at lunar-relative inclination
near 90° (abstract, p. 1; Conclusions, p. 5).

- **2-petal** (Sec. 2.1, p. 2): orbits twice per lunar sidereal period (27.3 d); the
  most stable family; 3-petal cyclers tend to collapse into this form and stay there
  "more-or-less indefinitely". Period of the simulated 2-petal cycler: **≈ 27 ± 3 days**
  (Sec. 5, p. 4). Lower perilune ⇒ higher Earth-relative inclination (Sec. 2.1, p. 2).
- **3-petal** (Sec. 2.2, p. 2; Sec. 5, p. 5): period ideal and NRHO/south-pole access
  good, but **unstable** — collapses readily into a 2-petal (or, less frequently,
  5-petal) form; periapsis consistently on the order of **200,000 km** with relatively
  low (< 35°) inclination; expensive to maintain in any configuration; not analyzed in
  detail by the authors (some collapses recoverable with low-thrust, Sec. 5, p. 5).
- **5-petal** (Sec. 2.3, p. 2; Fig. 3, p. 3): lunar flyby once every **53 days**, access
  to the *same* pole once every **163 days**; family examined over perilune
  **7,500–15,000 km** (Fig. 3 caption, p. 3); greatest accessibility flexibility but
  highly variable inclination (Sec. 5, p. 5).
- **4-petal:** "no stable four-petal cyclers were found that could meet the period and
  inclination requirement" (Sec. 2, p. 2). **> 5 petals:** considered infeasible due to
  the long periods (Sec. 2, p. 2). Explicit negatives — record both.
- **Stability domain** (Sec. 5, pp. 4–5): the 2-petal cycler is "relatively stable and
  periodic so long as the perilune passage was greater than 7500 km and less than
  ~15,000 km in the ±ẑ direction"; cyclers starting outside this domain had excessive
  maintenance cost or destabilized within a few cycles. Perilune **< 7,500 km** leads to
  instability — pushed into the 2-petal formation or a non-moon-resonant periodic orbit
  (Sec. 2.3, p. 2). At **20,000 km** perilune, lunar gravity was insufficient to
  maintain high lunar inclination (Sec. 5, pp. 4–5).
- **Encounter cadence / polar access:** perilune passes alternate between north and
  south lunar poles, so a given pole (and hence Gateway/south-pole departure geometry)
  is reachable **every other cycle** (Sec. 4, p. 4).
- **Specific Earth-launch geometry:** 2-petal with 20,000 km perilune has Earth-relative
  inclination ≈ **28.5°** (KSC latitude; apogee-scenario launch, Sec. 3, p. 3). The
  5-petal high-inclination variant with 10,000 km perilune has perigee altitudes as low
  as **600 km** at Earth inclination **80.5°**, allowing direct perigee injection but
  requiring a VSFB polar launch (Sec. 3, p. 3).
- **Prograde + retrograde pair:** Fig. 4 (p. 4) shows a specific 5-petal
  prograde/retrograde cycler pair with similar periapses providing NRHO-apolune access;
  caption notes the trajectories are un-optimized with no mid-course corrections "yet
  demonstrate a high degree of repeatability", and that the eccentricity direction
  points outside the Moon's SOI at these velocities.

## 3. ΔV and operations budgets (mining lens b)

**Table 2, p. 5 — ΔV requirements (m/s) for the three primary configurations:**

| Event | 2-petal (20,000 km) | 5-petal low-inc (25,000 km) | 5-petal high-inc (10,000 km) |
|---|---|---|---|
| Cycler rendezvous | 349.7 | 285.0 | 0 |
| LLO transfer | 1,254.1 | 1,291.0 | 1,481.5 |
| NRHO transfer | 446.3 | 529.7 | 742.2 |

(The high-inc 5-petal's zero rendezvous ΔV is the direct perigee-injection scenario,
Sec. 3, p. 3. The 2-petal apogee scenario costs 349.7 m/s plus a 5.5-day transit from
launch to cycler rendezvous, Sec. 3, p. 3.)

- **Maintenance / station-keeping:** lifetime ΔV over 5 cycles is **< 50 m/s for all
  5-petal family members**, "about half the cost of two-petal cyclers when considering
  that each cycle of [a] five-petal cycler is nearly twice as long as that of a
  two-petal" (Sec. 5, p. 5, citing Fig. 8; 2-petal totals in Fig. 6, p. 5, same < 50 m/s
  axis scale). Course-correction burns occur during perilune passage and during perigee
  passage between the second and third petal, assumed impulsive (Sec. 2.3, p. 2).
  Figs. 6–9 give the per-maneuver breakdowns (MCC1…MCC5, Flyby1…Flyby4) and the
  per-flyby Earth inclination/periapsis variability, but only graphically — no table.
- **Gateway rendezvous** (Sec. 4.1, p. 4): three burns — lunar flyby burn, NRHO
  insertion burn, phase-correction burn. The cycler's geometry yields an ascending node
  ~180° from the desired NRHO; a post-lunar-flyby orbit over two revolutions about the
  Moon (~12 days) rotates the node before insertion into the **9:2 lunar synodic
  resonance NRHO** (Gateway). NRHO-transfer costs in Table 2 above.
- **Lunar landing** (Sec. 4.2, p. 4): simplified three-burn solution — lunar flyby burn
  targeting **100 km perilune at 90° inclination**, circularization lowering apolune to
  100 km, then the final landing burn (landing itself excluded from the ΔV requirement).
- **Launch capability** (Table 1, p. 3, kg delivered to the three cyclers above, in the
  same column order): Falcon 9 RTLS 1925/1925/1440; Vulcan VC0 2260/2260/1695; Falcon 9
  ASDS 3520/3520/2640; Vulcan VC2 6160/6160/4620; Falcon Heavy (recover)
  7030/7030/N/A(5270); New Glenn 7095/7095/N/A(5320); Vulcan VC4 8825/8825/6615; Vulcan
  VC6 11,180/11,180/8385; Falcon Heavy (expend) 15,545/15,545/N/A(11,655). The N/A
  values: the high-inclination orbit needs a VSFB launch for which Falcon Heavy / New
  Glenn infrastructure does not exist (footnote 1, Table 1, p. 3). Apogee injection from
  KSC gives ~33% more upmass than near-polar perigee injection from VSFB (Sec. 3, p. 3).

## 4. Dedup vs catalogue row (mining lens c)

- **Same source, same family set.** The catalogue row `wittal-2022-em-cycler-family`
  cites NTRS 20220013595, which is this exact paper (IAC-22-C1.6.6). This mining
  upgrades the row's provenance from abstract-only to full-text; it does not introduce a
  new or extended family relative to the row's intent. No dedup action needed.
- The row's note "The full IAC PDF was not accessible to this ingest" is now obsolete
  and the period note can be sharpened: the **2-petal** member has period ≈ 27 ± 3 d
  (Sec. 5, p. 4); the 5-petal cadence is 53 d per flyby (Sec. 2.3, p. 2). Inclination
  ~90° lunar-relative is confirmed as the family's defining parameter (abstract;
  Conclusions, p. 5), with the stability domain perilune 7,500–15,000 km.
- **Adjacent rows:** Genova & Aldrin 2015 is cited as prior work (ref [6]) — Wittal's
  3-petal family is geometrically reminiscent of the Genova-Aldrin 3-petal cycler but is
  treated here as an unstable transitional form; keep the rows separate.

## 5. Actions / leads

1. **Backfill verdict (the headline):** the CR3BP backfill for
   `wittal-2022-em-cycler-family` stays blocked — NO_SOURCED_IC confirmed against the
   full text. Additionally flag the row's `model_assumption: cr3bp` as unsupported (the
   paper never states its model; internal evidence leans ephemeris-class simulation).
2. **Catalogue enrichment (v4.2 checklist pass):** period note (2-petal ≈ 27 ± 3 d),
   stability domain (perilune 7.5–15 Mm), maintenance ΔV (< 50 m/s per 5 cycles class),
   encounter cadence (5-petal: flyby/53 d, same pole/163 d), and the Table 2 ΔV block
   are all now full-text-sourced and quotable. `fleet_size`: still not stated by the
   source — keep null. V∞-at-encounters: inapplicable (these are bound Earth-Moon
   orbits, not hyperbolic-encounter cyclers) — record as n/a rather than not-extracted.
3. **Acquisition lead (the one remaining source-axis hope for ICs):** ref [9] — Wittal,
   M. M., Smith, J. D., and M., C. A., "Mission Design Considerations for Robotic Lunar
   and Gateway Payload Return," AAS/AIAA Astrodynamics Specialist Conference, 2021. The
   companion paper underlying the accessibility quantification (cited in Sec. 2, p. 2);
   it may state the simulation environment and possibly orbit states.
4. **Negatives recorded:** no stable 4-petal meeting period+inclination requirements;
   > 5 petals infeasible (both Sec. 2, p. 2); 3-petal unstable in all configurations
   tried (Sec. 5, p. 5); no epochs, no states, no Jacobi constants, no model statement
   anywhere in the paper.
