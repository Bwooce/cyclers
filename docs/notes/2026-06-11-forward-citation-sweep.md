# Forward-citation sweep — newer works citing our load-bearing holdings (2026-06-11)

**Task:** forward-trace the key papers behind the catalogue (2015–2026, emphasis post-2020) and
triage which citing works are worth acquiring.
**Method:** per-seed web searches (Google-Scholar-style + site-targeted), abstract pages followed
where reachable. The Semantic Scholar citations API was unavailable during the session (fetch
tooling outage), so coverage relies on search-engine recall; a follow-up API pass would tighten
completeness, particularly for Russell & Ocampo (the most-cited seed).
**Access legend:** `free` = arXiv / NTRS / university PDF / open repository; `gated` = AIAA ARC,
Springer, ScienceDirect (human acquisition needed).
**Triage lens (strict):** HIT only if it provides (a) new cycler members with published numbers,
(b) real-ephemeris reproduction / maintenance-ΔV of known cyclers, (c) a new construction method,
or (d) Earth-Moon / moon-tour cycler extension. Passing citations and numbers-free mission
concepts are marked marginal or omitted.

---

## Seed 1 — Russell & Ocampo 2004/2005 (Earth-Mars cycler classification, JGCD + PhD)

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Fornari, E. & Pontani, M., "A global exploration method to identify families of cycling Earth–Mars trajectories," *Aerotecnica Missili & Spazio*, DOI 10.1007/s42496-020-00050-6 | 2020, AM&S (Springer) | Independent global-search METHOD for Earth-Mars cycler families — directly comparable to (and a cross-check of) the Russell-Ocampo enumeration | **HIT** (method + possible new family numbers) | gated (Springer; AM&S sometimes posts PDFs) |
| Pelle, S., Gargioli, E., Berga, M., Pisacreta, J., Viola, N., Dalla Sega, A., Pagone, M., "Earth-Mars cyclers for a sustainable human exploration of Mars," *Acta Astronautica* 154:286–294 | 2019, Acta Astro | Cycler-class trade-off + baseline selection for a crewed architecture; cites the classification work; may contain per-class trajectory parameters | marginal-HIT (architecture, but with a class trade table) | gated (ScienceDirect) |
| "Optimization of Earth-Mars transfer trajectories with launch constraints," *Astrophys. Space Sci.*, DOI 10.1007/s10509-021-04025-2 | 2021 | Transfer optimization, cites cycler survey in passing | NOT a hit | gated |
| Pontani, M. & Conway, B.A., "Optimal Trajectories for Hyperbolic Rendezvous with Earth–Mars Cycling Spacecraft," *JGCD* 41(2) | 2017/2018, JGCD | Taxi-to-cycler hyperbolic rendezvous optimization (ops layer, not new members); companion "Optimal low-thrust hyperbolic rendezvous for Earth-Mars missions," *Acta Astro* (2018, pii S0094576517316454) | marginal (ops/GNC-adjacent, but quantitative) | gated (AIAA / ScienceDirect) |

## Seed 2 — McConaghy, Longuski et al. 2002–2006 (S1L1 / two-synodic cyclers)

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Wilde, M., Patel, B., Kish, B., "Design Considerations for an Earth–Mars Cycler Spacecraft Using the S1L1 Cycler," *J. Spacecraft & Rockets* 59(3), DOI 10.2514/1.A35160 | 2021 (online) / 2022, JSR | Journal follow-on of the FIT S1L1 vehicle work; uses the nominal 154-day S1L1 leg + 1032-day contingency; may restate S1L1 itinerary numbers in a citable journal form | marginal-HIT (S1L1 cross-reference; mostly vehicle sizing) | gated (AIAA) |
| Wilde, M. et al., "Parametric Design of a Crew Transfer Vehicle for Earth–Mars Cyclers," *JSR*, DOI 10.2514/1.A34637 | ~2020, JSR | Taxi-vehicle sizing for cycler rendezvous; cites the V∞ tables | NOT a hit (vehicle design, no new trajectories) | gated (AIAA) |
| "Development of Practical Earth-Mars Cycler Trajectories," ERAU Discovery Day (Prescott), presentation 34 | 2023, ERAU student forum | Student re-derivation of practical cyclers | marginal (likely no archival numbers) | free (commons.erau.edu) |

## Seed 3 — Genova & Aldrin 2015 (Earth-Moon 3-petal cycler) + Earth-Moon cycler lane

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Ross, S.D. & Roberts-Tsoukkas, M., "Stable, Low-Energy Prograde Earth-Moon Cycler Orbits," AAS 25-621, AAS/AIAA Astrodynamics Specialist Conf. | 2025, AAS | **New class of fully ballistic, STABLE, prograde Earth-Moon cyclers** (first stable ballistic E-M cyclers); systematic family construction with specified Earth/Moon orbit counts; resonant + near-chaotic regimes | **HIT (top)** | free (https://ross.aoe.vt.edu/papers/ross-roberts-tsoukkas-2025-AAS-25-621.pdf) |
| Roberts-Tsoukkas, M. & Ross, S.D., "Stable Prograde Earth-Moon Multi-Orbiter Cyclers via Three-Body Dynamics" (journal version of the above) | 2026 | Extended journal treatment of the same family — likely the better golden-test source | **HIT (top)** | free (https://vsgc.odu.edu/wp-content/uploads/2026/04/Roberts-Tsoukkas_Michael_Cycler-Journal-Paper.pdf) |
| Wittal, M.M., Miaule, S., Asher, B.W., "Earth-Moon Cycler Mission Design for Lunar Logistics," IAC-22-C1.6.6 | 2022, IAC | Families of Earth-Moon cyclers with near-polar lunar inclination; logistics application; cites Genova-Aldrin | **HIT** (Tier-2 members) | free (NTRS 20220013595) |
| "A cislunar in-orbit infrastructure based on p:q resonant cycler orbits," *Acta Astronautica* 170:539–551 | 2020, Acta Astro | p:q resonant E-M cycler orbits built by differential correction in PCR3BP, refined by multiple shooting in the bicircular model; Lambert links to LEO/LLO parking | **HIT** (Tier-2 members + multiple-shooting method) | gated (ScienceDirect, pii S0094576520300916) |
| "Preliminary investigation and proposal of periodic orbits and their utilization for logistics in the cislunar regime," (pii S0265964624000262) | 2024 | 14 sample cislunar periodic orbits w/ stability indices for logistics | marginal (survey-grade) | gated |
| "Integrated Orbital Design Method for Manned Lunar Exploration with Relaxed Temporal Constraints," *JSR*, DOI 10.2514/1.A36062 | ~2024/2025, JSR | Multistage lunar trajectory splicing; cites circumlunar cycler work in passing | NOT a hit | gated (AIAA; pdf link exposed on ARC) |
| Gupta, M., "Navigating Chaos: Resonant Orbits for Sustaining Cislunar Operations," Purdue PhD (Howell group) | 2024 | Resonant-orbit catalogue for cislunar ops — background for the E-M resonant cycler lane | marginal-HIT (method/back-fill, not cyclers per se) | free (engineering.purdue.edu Howell publications page) |
| "Cislunar Resonant Transport and Heteroclinic Pathways: From 3:1 to 2:1 to L1," arXiv:2509.12675 | 2025, arXiv | Resonant transport mechanics underpinning E-M cycler family transitions | marginal-HIT | free (arXiv) |
| "Mass-Optimal Low-Thrust Forced Periodic Trajectories in the Earth-Moon CR3BP," arXiv:2502.05140 | 2025, arXiv | Forced (low-thrust) periodic orbits in E-M CR3BP — method for powered-cycler analogues | marginal-HIT | free (arXiv) |

## Seed 4 — Byrnes, Longuski, Aldrin 1993 (Aldrin cycler)

Forward citations are dominated by the same works as Seeds 1–2 plus encyclopedic coverage
(Wikipedia/Grokipedia/Marspedia). The only NEW quantitative threads found: the Pony Express line
(Seed 9-adjacent, below) which uses the Aldrin-class cycler as the relay backbone, and Howe 2025
(Seed 10, already held). No additional acquisition-worthy citer unique to this seed.

## Seed 5 — Landau & Longuski (human Mars transportation / cycler V∞ tables)

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Pelle et al. 2019 (see Seed 1) | 2019 | Architecture trade citing the V∞ tables | marginal-HIT | gated |
| Wilde et al. 2021/2020 (see Seed 2) | 2020–2022 | Vehicle design consuming the V∞ tables | marginal | gated |
| Donahue & Duggan, Boeing IAC-22 Mars flyby (already held) | 2022 | — | held | — |

No new trajectory-numbers citer found beyond what other seeds surfaced. The Landau-Longuski
citation graph mostly fans out into architecture/ISRU papers (not hits under the lens).

## Seed 6 — Jones, Hernandez, Jesick 2017 (VEM triple cyclers, AAS 17-577)

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| (Yang, H. group, NUAA), "Callisto–Ganymede–Europa Triple Cyclers," *JGCD*, DOI 10.2514/1.G008387 | 2024, JGCD | **New CGE triple-cycler members** — direct methodological descendant of the 2017 VEM/IEG triple-cycler papers; Jovian-moon cycler numbers (moon-tour lane) | **HIT** | free PDF circulating on ResearchGate (pub. 383986230); gated at AIAA |
| Hernandez, S., Jones, D.R., Jesick, M., "One Class of Io-Europa-Ganymede Triple Cyclers," AAS 17-462 (Adv. Astronaut. Sci. 162, pp. 973–984) | 2017, AAS | Companion to AAS 17-577 — NOT a forward citation but un-held and feeds the same moon-tour cycler lane (1:2:4 Laplace resonance, exact repeating cyclers) | **HIT** (backfill) | semi-free (Semantic Scholar hosts PDF) |

## Seed 7 — Hughes, Edelman & Longuski 2014 (fast Mars free returns, AIAA 2014-4109)

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Hughes, K.M. (et al.), "Fast Free Returns to Mars and Venus with Applications to Inspiration Mars," *JSR*, DOI 10.2514/1.A33293 | ~2015–2017, JSR | Archival journal version of the held conference paper — extended tables of free-return candidates (incl. EVME 2021 / 2023 opportunities) usable as golden free-return targets | **HIT** (validation-grade tables) | gated (AIAA) |

## Seed 8 — Şaloğlu & Taheri 2023/2025 (iso-impulse)

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Şaloğlu, K., "On the Connection of Optimal Impulsive and Low-Thrust Trajectories," PhD dissertation, Auburn Univ. (etd 10415/10334) | 2026 | Consolidates the iso-impulse line + impulsive→low-thrust homotopy; likely contains the fullest statement of the ΔV-allocation machinery and extra worked examples | **HIT** (method; supersedes both held papers as a single source) | free (Auburn etd) |
| "Acceleration-Based Switching Surfaces for Impulsive Trajectory Design Between Cislunar Libration Point Orbits," *J. Astronaut. Sci.*, DOI 10.1007/s40295-024-00432-z | 2024, JAS | Same group; impulsive design between cislunar LPOs — tangent to Tier-2 | marginal | gated (Springer) |

## Seed 9 — Rogers et al. 2012 (cycler establishment / V∞-leveraging) + Pony Express line

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Rogers et al., "Establishing cycler trajectories between Earth and Mars," *Acta Astronautica* 112:114ff | 2015, Acta Astro | Archival journal version of the held AIAA 2012-4746 — establishment ΔV via VILM or low thrust; likely refined/extended numbers | **HIT** (validation: establishment-ΔV) | gated (ScienceDirect, pii S0094576515000831) |
| Pascarella, A. et al., "Low-thrust trajectory optimization for the solar system pony express," *Acta Astronautica* 203:280ff, DOI 10.1016/j.actaastro.2022.11.046 (approx.) | 2023, Acta Astro | Journal version of held AAS 22-015: SEP insertion into an E-M cycler (36 kg prop, 500 kg s/c) + 2 kg for 8 maintenance flybys over 6 yr — real-ephemeris maintenance-ΔV of a known cycler | **HIT** (maintenance-ΔV validation; held AAS version may suffice — diff before buying) | gated (journal); AAS version free at ai.jpl.nasa.gov |
| "Low-thrust roundtrip trajectories to Mars with one-synodic-period repeat time," *Acta Astronautica* (pii S0094576515000107) | 2015, Acta Astro | Low-thrust one-synodic cycler-like roundtrips | marginal-HIT (powered-cycler lane) | gated |
| Miguel, N., Colombo, C., Dalla Vedova, F., "Systematic Construction of Solar-Sail-Based Stopover Cyclers," *J. Astronaut. Sci.* 70:6, DOI 10.1007/s40295-023-00372-0 | 2023, JAS | Systematic METHOD for sail-based stopover cyclers between two bodies (demonstrated Earth–Main Belt) — extends the stopover/cycler construction toolbox | **HIT** (method; low-thrust/sail cycler lane) | gated (Springer); ProQuest openview page exists |

## Seed 10 — Howe 2025 (ICES-2025-555 escalator/cycler, held)

Too recent for academic forward citations. Found only popular coverage (Space Settlement
Progress blog, Nov 2025: "Novel design of a Mars Cycler"). Re-check in 12 months.

## Bonus (method lane, surfaced during sweep)

| Work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Ozaki, N. et al., "Asteroid Flyby Cycler Trajectory Design Using Deep Neural Networks," *JGCD* 45(8), DOI 10.2514/1.G006487 | 2022, JGCD | NN surrogate replacing inner trajectory optimization in a global cycler search — ML-seeded periodic multi-flyby construction (genome roadmap) | **HIT** (method) | free (arXiv:2111.11858) |
| "Global Optimality in Multi-Flyby Asteroid Trajectory Optimization: Theory and Application Techniques," *JGCD*, DOI 10.2514/1.G009335 | ~2025, JGCD | Global-optimality theory for multi-flyby problems | marginal-HIT (method) | gated (check arXiv) |

---

## Ranked acquisition shortlist (by expected catalogue / validation payoff)

1. **Ross & Roberts-Tsoukkas, AAS 25-621 + 2026 journal version** (Earth-Moon stable prograde
   cyclers) — new ballistic, *stable* E-M cycler families with a systematic construction method;
   feeds Tier-2 catalogue directly with publishable member numbers. Both PDFs free. **Fetchable.**
2. **Callisto–Ganymede–Europa Triple Cyclers, JGCD 2024 (10.2514/1.G008387)** — new triple-cycler
   members with numbers; extends the held VEM triple-cycler line to the Jovian moon-tour lane.
   Free PDF on ResearchGate. **Fetchable.**
3. **Wittal, Miaule & Asher, IAC-22-C1.6.6** (Earth-Moon cycler lunar logistics, near-polar
   families) — Tier-2 members + ops context; NTRS. **Fetchable.**
4. **Şaloğlu 2026 Auburn PhD dissertation** — single consolidated source for the iso-impulse →
   low-thrust connection; richer worked examples than the two held papers. Free etd. **Fetchable.**
5. **Rogers et al., Acta Astronautica 112 (2015)** — archival establishment-ΔV (VILM + low-thrust)
   numbers for Earth-Mars cyclers; upgrades the held conference version for validation use.
   Gated (human acquisition).
6. **Fornari & Pontani, AM&S 2020 (10.1007/s42496-020-00050-6)** — independent global-search
   method for Earth-Mars cycler families; cross-check potential against the Russell-Ocampo
   enumeration. Gated (Springer).
7. **Miguel, Colombo & Dalla Vedova, JAS 70:6 (2023)** — systematic stopover-cycler construction
   (sail); method feed for the powered/low-thrust cycler lane. Gated.
8. **Hughes et al., JSR (10.2514/1.A33293)** — archival fast free-return tables (Inspiration-Mars
   era opportunities); golden free-return targets. Gated.
9. **Ozaki et al., JGCD 2022** — ML-surrogate global cycler search (genome roadmap). arXiv free.
   **Fetchable.**
10. **"A cislunar in-orbit infrastructure based on p:q resonant cycler orbits," Acta Astro 170
    (2020)** — p:q resonant E-M cyclers w/ multiple-shooting bicircular refinement; Tier-2
    members + method. Gated.
11. **Wilde, Patel & Kish, JSR 2021 (10.2514/1.A35160)** — S1L1 journal cross-reference (low
    priority: trajectory content likely duplicates held sources). Gated.
12. **Pascarella et al., Acta Astro 203 (2023)** — only if it adds numbers beyond the held
    AAS 22-015 (diff abstract/tables first). Gated.

### Sweep counts
- Seeds traced: **10** (+1 bonus method lane)
- Web searches run: **18**
- Citing/related works triaged: **~30**
- Strict HITs: **13** (incl. 2 backfill/bonus)
- Marginal: **~10**; non-hits noted and dropped: ~7
- Already-held items re-encountered and deduplicated: Howe 2025, AAS 22-015, AIAA 2012-4746,
  AIAA 2014-4109, Genova-Aldrin 2015 (both), McConaghy 2004, Russell 2004, Patel 2019,
  Saloglu 2023/2025, Jones 2017, Agrawal 2022, Donahue 2022.

### Caveats
- Search-engine recall only; the Semantic Scholar citations API was unreachable during the
  session. A follow-up API pass on DOIs 10.2514/1.1011, 10.2514/1.15215, 10.2514/1.A35091,
  10.2514/6.2014-4109 would close the completeness gap for 2024–2026 citers.
- Ranks 1–4 and 9 are immediately fetchable (free PDFs); ranks 5–8 and 10–12 need human-gated
  acquisition (AIAA ARC / Springer / ScienceDirect).
