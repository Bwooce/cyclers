# #943 cell gc: prior-art search for two-working-body Ganymede-Callisto cyclers (2026-10-06)

Requested by the lead after twobody-gen-opus's cell-gc result (commit `e3cab2c7`; candidates gc-1 and
gc-2 in `docs/notes/2026-10-05-942-943-two-working-body-generator.md` sec. 6.11). Search by
corpus-file-opus.

The candidates, both k = 3 (37.57 d, 3 G-C synodic periods), both moons bending:
- gc-1: v_inf G / C = 2.397 / 1.807 km/s; G-G 1-rev low, G->C->G on one conic, Callisto 1:1.
- gc-2: v_inf G / C = 3.617 / 3.039 km/s; G-G 1-rev low, G->C, C-C 1-rev high, C->G; Callisto turns
  2 x 6.87 deg at 9,800 km (ratio 0.259).

## 1. Method

- **Citation graph (OpenAlex):**
  - all 33 works citing Russell & Strange 2009 (JGCD 32(1), doi 10.2514/1.36610);
  - all 10 works citing AAS 07-118 (OpenAlex W1503601394).
  - Crossref's cited-by data is not public, so OpenAlex replaced it.
- **OpenAlex topic searches:** "ganymede callisto cycler", "double cycler jupiter moons", "two-moon
  cycler", "jovian moon cycler", "callisto ganymede resonant transfer tour".
- **Web searches:** G-C cycler, "Ganymede-Callisto cycler", JUICE and Europa Clipper tour cyclers,
  Lynam, Grushevskii, Lantukh, Hernandez/Jones/Jesick follow-ups, "massive target" cycler, 2025-2026
  papers.
- **arXiv API:** cycler with Ganymede, Callisto, Jovian or Jupiter. One hit, Braik & Ross 2026
  (arXiv 2605.31543), which is Earth-Moon only.
- **NTRS API and site search:** no cycler hits. **JPL TRS:** AAS 07-118 only.
- **Open-access full texts searched:** Cangahuala et al. 2025 (Space Sci. Rev.) and Boutonnet et al.
  2024 (Space Sci. Rev.).
- **Held papers re-read for G-C content:** Hernandez-Jones-Jesick AAS 17-608, Lynam & Longuski 2011,
  Liang et al. 2024, Russell & Strange 2007/2009, Yang et al. 2023 (review), Russell 2012 (survey).
- **Held status:** every hit was run through `scripts/check_wanted_vs_corpus.py`.

## 2. Verdict

**No published two-working-body Ganymede-Callisto cycler matching gc-1 or gc-2 was found.** The only
published two-working-body G-C cycler is still Campagnola et al. 2019 GCGC (v_inf 3.5/4.5, a different
class; already excluded in sec. 6.11).

There is one flag that needs the owner's explicit comparison:
- **Russell & Strange GanCal#1 and GanCal#5 (held)** have the SAME period as gc-1/gc-2: 37.6 d,
  3 synodic periods (AAS 07-118 Table 3; R-S 2009). In them Ganymede is the flyby body and Callisto is
  massless.
- Their v_inf are G/C 3.18/3.26 (#1) and 3.24/3.34 (#5) km/s.
- **gc-2 is within 0.5 km/s of both at both moons:**
  - against #1: 0.44 at G, 0.22 at C;
  - against #5: 0.38 at G, 0.30 at C.
- The leg structures differ:
  - GanCal#5 = g(1.504, 541.5 deg, L) G(3.747, 628.9 deg, U), with no Callisto return.
  - gc-2 has a C-C 1-rev-high return, and Callisto turns only 6.87 deg twice (ratio 0.259).
- So gc-2 may be a two-working-body relative of the R-S GanCal family rather than an unrelated object
  (INFERRED). Sec. 6.11 records GanCal#5 as LITERAL (recovered in-run) but does not state the gc-2 to
  GanCal distance; it should.
- gc-1 (2.40/1.81) is more than 0.78 km/s away from both GanCal members at both moons.

## 3. Findings, by relevance

| # | Work | Held? | G-C content | Collides with gc-1 / gc-2? |
|---|---|---|---|---|
| 1 | Russell & Strange 2007 (AAS 07-118) / 2009 (JGCD 32(1), 10.2514/1.36610) | held | GanCal#1, #5: one-working-body (Callisto massless), 37.6 d | no literal collision; gc-2 near (sec. 2) |
| 2 | Campagnola et al. 2019 (JGCD 42(12), 10.2514/1.G004309) | held | GCGC two-working-body cycler, v_inf 3.5/4.5, ~90 deg rotation per cycle | no (different class; sec. 6.11) |
| 3 | Liang, Yang, Li, Bai & Qin 2024 (JGCD, 10.2514/1.G008387) | held | CGCEC switched "double cycler": the C-G half is ~50 d (4 synodic), multi-rev legs, both moons flown by; v_inf C 5.67, G 6.99 (Table 3) | no (different period, v_inf ~2x higher, three-moon sequence) |
| 4 | Boutonnet, Langevin & Erd 2024, "Designing the JUICE Trajectory", Space Sci. Rev. 220:67, 10.1007/s11214-024-01093-y (open access) | not held | "Callisto-Ganymede-Callisto round trips" (p.33, Figs. 26-27): one-shot v_inf-reduction endgame, C 5.22 -> 3.5 -> 2 km/s, G to ~1.5 km/s | no (non-periodic, pumping); the nearest one-shot precursor in v_inf to gc-1 |
| 5 | Cangahuala et al. 2025, "Europa Clipper Mission Design, Mission Plan, and Navigation", Space Sci. Rev. 221:22, 10.1007/s11214-025-01140-2 (open access) | not held | p.14: the sub-Jovian switch "can be accomplished ... [by] different types of cyclers"; the baseline 21F31_V6 uses a Callisto petal rotation | no |
| 6 | Campagnola, Buffington, Anderson, Pellegrini, Restrepo & Lam 2025, "Europa Clipper Mission Design: Design of the 21F31 Reference Tour", JAS 72, 10.1007/s40295-025-00527-1 (closed) | not held | the reference-tour paper; the likeliest home of the G-C cycler details promised in Campagnola 2019 (abstract does not mention cyclers) | unknown; acquire |
| 7 | Takubo, Campagnola, Pellegrini & Anderson 2026, JGCD, 10.2514/1.G009868 (also AIAA SciTech 2026-1262) | not held | contingency escapes from the Clipper multi-moon tour | unlikely; acquire to confirm |
| 8 | Lantukh & Russell 2012, "Automated Inclusion of n-pi Transfers in Gravity-Assist Flyby Tour Design", AIAA 2012-4749, 10.2514/6.2012-4749 | not held | n-pi tour pathfinding with Jupiter examples; cites R-S 2009 | unknown; acquire |
| 9 | Golubev, Grushevskii, Koryanov & Tuchin 2014, J. Comput. Syst. Sci. Int. 53:445-463, 10.1134/S1064230714030083 | not held | "crossed" G->C->G gravity assists for capture and pump-down (one-shot) | no (per abstracts); background |
| 10 | Lynam 2012, PhD dissertation, Purdue (AAI3545317) | not held | Laplace-resonant (Io-Europa-Ganymede) triple cyclers and capture sequences; no Callisto cyclers (abstract) | no |
| 11 | Scott, Ellison, Bokelmann & Ozimek 2025, "Europa Clipper Mission Design: Pump Down Trajectory Design", JAS, 10.1007/s40295-025-00534-2 (closed) | not held | pump-down optimisation (the abstract names no moons) | unlikely |
| 12 | Hernandez, Jones & Jesick 2017 (AAS 17-608); Lynam & Longuski 2011 (Acta 69) and 2010 (AIAA 2010-8256) | held (2010 as the 2011 journal form) | Io-Europa-Ganymede only; Lynam notes Callisto "rarely" encountered | no |
| 13 | Kumar, Anderson & de la Llave 2023 (Acta 211:76) | held | Ganymede-Europa resonant tori (CR4BP) | no (G-E cell context) |
| 14 | Russell 2012 survey (JGCD 35(3)); Yang et al. 2023 review; Anderson et al. 2014/2018 petal rotation | held | reviews and petal rotation | no |

No relevant hits in the other OpenAlex citers of R-S 2009, which are Earth-Moon, asteroid, Enceladus,
Titan and methods papers.

## 4. Recommendations

- **Owner adjudication of gc-2:** compare it against GanCal#1/#5 explicitly (same k = 3 period; v_inf
  within 0.5 km/s at both moons). Question: is gc-2 the two-working-body continuation of a GanCal member?
- **Acquire,** in priority order:
  1. Campagnola et al. 2025 JAS 21F31.
  2. Lantukh & Russell 2012 AIAA 2012-4749.
  3. Takubo et al. 2026 JGCD.
  4. Golubev et al. 2014.
  5. Lynam 2012 dissertation.
- **File:** the two open-access Space Science Reviews papers (Boutonnet 2024; Cangahuala 2025) were
  downloaded and can be filed and digested now.

## 5. Follow-up 2026-10-06: the two open-access papers filed and digested

- Boutonnet, Langevin & Erd 2024 (JUICE): `2026-10-06-digest-boutonnet-langevin-erd-2024-designing-juice-trajectory.md`.
  - Single C-G-C round trip: G v_inf about 3.3 km/s; C arrival 1.9-2.3 and outbound 1.8-2.4 km/s; a
    50-day Ganymede-to-Callisto window cycle.
  - One-shot; no collision. Nearest one-shot relative by v_inf (Callisto overlaps gc-1; Ganymede within
    0.3 km/s of gc-2).
- Cangahuala et al. 2025 (Clipper): `2026-10-06-digest-cangahuala-2025-europa-clipper-mission-design-plan-navigation.md`.
  - The flown 21F31_V6 switch is a "Ganymede-Callisto Petal Rotation" (Table 2, page image): C01, G06,
    C02-C07, G07, C08, C09; Callisto v_inf 3.53-3.59 mid-sequence, Ganymede 5.35 and 4.30.
  - Non-resonant and non-repeating, so not a cycler; no collision with gc-1 or gc-2.
- The gc-2 / GanCal#1/#5 flag (sec. 2) was sent to twobody-gen-opus at the lead's request, with the
  pre-registered Callisto-mass-homotopy question.
