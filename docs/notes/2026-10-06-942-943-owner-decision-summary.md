# #942 / #943 two-working-body cyclers: decision summary for the owner (2026-10-06)

Author: twobody-gen2-opus. Source of truth: `docs/notes/2026-10-05-942-943-two-working-body-generator.md`
(sections cited as "6.x"). This note asks one decision per candidate under the `#875` novelty policy
(`docs/notes/2026-09-07-875-novelty-policy-decision.md`): candidate-novel with attribution (to whom), or not.
No catalogue rows have been written.

## 1. How the candidates were made and checked

- Generator: the Hollister & Menning 1970 date-residual corrector, generalised to two working bodies.
  - Ideal circular-coplanar models: R-S 2007/2009 Table 2 constants for the Jovian and Venus-Mars
    cells; Venus 0.61520 yr for Earth-Venus; Mars 1.875 yr for Earth-Mars.
  - Return catalogue: full-rev n:m, half-rev, and generic 1-rev returns.
  - Exact zeros (< 1e-8 km/s).
  - The #888/#937 demanded-turn gate at the registry floors, with minimax free directions.
  - An independent DOP853 re-fly. Sec. 6.1.
- Gauntlet per cell:
  - literal collisions with R-S, H&M, Rall, Campagnola 2019, Jones 2017, Pisarevsky 2008 and the
    catalogue;
  - the offline literature_check;
  - SOI self-consistency.
- Recall controls passed in-run or targeted:
  - VenMar#45 (vm, vm2, vmn);
  - Hollister 1H/2H (ev);
  - GanCal#5 (gc);
  - GanEur#43 and EurGan#131 (ge);
  - R-O 2.5.1.+0 (em, 6.43);
  - GanCal#1 and GanEur#5 (targeted, sec. 5);
  - GanEur#316 (targeted, 6.28).
- Real-ephemeris rung (d), patched conic, 5 epochs 2030-2056 unless stated.
  - Heliocentric: Standish J2000 mean elements, plus DE440 where noted. Controls: VenMar#45 (C1), the
    H&M endpoint orbits 1/2/11/12/13/15, D1 (Hollister 1H, direct route).
  - Jovian: NAIF jup365, 10 cycles. Controls at R-S's own epochs:
    - GanEur#316@2019, half-rev path, gate pass (6.38, 6.50);
    - GanCal#1@2013, full-rev path, closes, gate indeterminate 0.9917. The owner's post-hoc ruling (6.45)
      accepted it as a path control.
- Fixes made on the way that bear on the results:
  - kepler_step (6.3);
  - the solver GM (6.23);
  - half-rev root selection (6.34; the ev set unchanged at 31/31);
  - fixed-leg extents (6.30);
  - the long-chain stall (6.33);
  - chain-length continuation (6.35), the anchored shoot (6.48).
- Restart-based "no other closure" statements made before 6.35 are WEAK.

## 2. Candidates

All V_inf in km/s. "Ratio" = demanded / available turn at the registry floor (gate pass < 1).
Ideal-model dates are the representative zero's Lambert-leg starts (model epoch: all bodies at angle 0
at t = 0).

### 2.1 gc-1 (Ganymede-Callisto, both bend) — real-ephemeris PASS

- Structure: k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s.
  - Period 37.570 d (3 G-C synodic periods). Dates 1.178102, 11.756379, 34.957696 d.
  - Legs: G-G 1-rev low return; G->C and C->G on one conic; Callisto 1:1 full-rev.
- V_inf: G 2.397, C 1.807.
- Turns: Ganymede 2 x 29.35 deg (ratio 0.646, 2,437 km); Callisto 2 x 40.15 deg (ratio 0.737, 1,801 km).
- Distance from Jupiter: 888,745-2,294,675 km (the r_max includes the full-rev leg, 6.30).
- Real ephemeris: jup365, 10 cycles (376 d), PASS at 5/5 epochs, worst 0.755-0.761; DOP853 re-fly
  <= 0.19 km (6.26). The full-rev path caveat was LIFTED by the owner's C4 ruling (6.45).
- Prior art (web, OpenAlex, arXiv and NTRS search, `2026-10-06-943-gc-prior-art-search.md`; plus
  21F31 ISSFD 2024, Golubev 2014 and 2017, Lam 2018, Buffington 2014, Minovitch 1972, Cangahuala 2025,
  Boutonnet 2024): no collision.
  - The only published two-working-body G-C cycler is Campagnola et al. 2019 GCGC (V_inf 3.5/4.5,
    alternating, high-fidelity), a different class.
  - Nearest one-shot relative: the JUICE C-G-C round trip (Boutonnet 2024), whose Callisto V_inf overlaps.
- Caveats: patched conic only (the Jovian n-body lane has no positive control, #968). Ideal margin 26 %.
- **Decision asked:** candidate-novel?
  - Proposed attribution: method, Hollister & Menning 1970 (date-residual corrector); ideal model and
    return taxonomy, Russell & Strange 2007/2009 (who named the massive target as future work, AAS
    07-118 p.18); class context, Campagnola et al. 2019 (the first two-working-body G-C cycler).

### 2.2 gc-2 (Ganymede-Callisto, both bend) — real-ephemeris PASS

- Structure: k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/0s.
  - Period 37.570 d. Dates 0.838679, 11.684556, 14.289304, 35.803635 d.
  - Legs: G-G 1-rev low; G->C and C->G of 2.605 d each on one conic; C-C 1-rev high.
- V_inf: G 3.617, C 3.039.
- Turns: Ganymede 2 x 19.72 deg (0.788, 1,024 km); Callisto 2 x 6.87 deg (0.259, 9,800 km).
- Distance: 791,455-2,337,392 km.
- Real ephemeris: jup365, 10 cycles, PASS at 5/5 epochs, worst 0.813-0.823; re-fly <= 2.4e-5 km
  (Lambert-only, 6.14).
- Prior art: as gc-1, no collision. Near relative, the R-S GanCal family (6.20-6.21):
  - gc-2 = GanCal#5's skeleton (shared G-G 1l leg and period), with GanCal#5's 24.25-d G->C leg
    replaced by a G->C leg plus a C-C 1-rev return: one extra Callisto encounter.
  - GanCal#5's leg passes Callisto at 4.5 SOI (170,800 km) near that time, so the inserted encounter
    has no counterpart on it.
  - gc-2 needs at least 0.214 of Callisto's GM.
  - V_inf distance: 0.38/0.30 to GanCal#5, 0.44/0.22 to GanCal#1.
  - In the patched conic no Callisto-mass path joins them.
- Caveats: patched conic only (#968); a continuous-gravity check of the GanCal relation is an open owner
  option.
- **Decision asked:** candidate-novel?
  - Proposed attribution: as gc-1, plus R-S GanCal#5 as the one-working-body skeleton it extends.
  - Or a "GanCal-family relative"? (This is the closest call of the set.)

### 2.3 ev-C (Earth-Venus, Earth almost massless; R-S architecture at Venus) — real-ephemeris PASS

- Structure: k2|LE>V/0s|LV>V/1h|LV>E/0s.
  - Period 1167.889 d (2 E-V synodic periods). Dates 70.459681, 193.288309, 1017.239797 d.
  - Legs: E->V 122.83 d; V-V 1-rev high, 823.95 d; V->E 221.11 d.
- V_inf: E 9.075, V 13.166.
- Turns: Venus 2 x 19.65 deg (0.746, 3,055 km); Earth turn 0 in the ideal model (0.06-0.11 on real
  elements).
- Distance: 0.509-1.651 AU.
- Real ephemeris: PASS at 5/5 on Standish (worst 0.869-0.914, re-fly <= 0.093 km) and on DE440 (the same
  ratios, re-fly <= 0.116 km), after the GM fix (6.23).
- Prior art (`2026-10-06-942-evC-prior-art-search.md`):
  - No refereed Venus-hosted Earth-Venus cycler.
  - AAS 07-118 ran no Earth-Venus set (the R-S 2009 claim is a citation slip).
  - Nearest unrefereed class relative: D. Ross "2L4" (powered, never reaches Earth, aphelion 0.979 AU,
    targets Sun-Earth L1, Venus flyby below the surface; no collision).
  - VanderVeen 1969 affirms the 8-yr E-V repeat. corpus-file-opus's own arithmetic gives a -2.4 deg slip
    per 8 yr, which the 16-yr chain absorbs. Not a published objection.
- Minovitch 1963 (JPL TR 32-464, 31 Oct 1963; key pages read on the images by corpus-file-opus, full
  digest in progress): NO collision.
  - Every return trajectory is one-shot.
  - Table 14 has ~1-yr Earth-Venus-Earth free returns with no Venus-Venus leg; no trajectory returns
    to Venus.
  - Concept priority, recorded: p.15 "example 5" imagines an E-V-M-E-V-M-E free fall that "repeats
    the same flight" (an idea only; whether it is possible is stated as unknown, and none is computed).
    The p.50-51 "space bus" chain (Table 23, E-V-M-E-M-E-V-E, 1970-75) does not repeat.
  - Cite Minovitch 1963 as the earliest statement of the repeating free-fall (cycler) idea.
- **Decision asked:** candidate-novel under `#875` (ii)?
  - The known architecture (R-S one-working-node) at a never-treated pair.
  - Proposed attribution: Russell & Strange 2007/2009, explaining its re-application to Earth-Venus
    with Venus as the host.

### 2.4 ev-A (Earth-Venus, both bend) — Standish PASS; DE440 near-ballistic

- Structure: k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s.
  - Period 1167.889 d. Dates 583.944387, 984.221555, 1576.257792 d.
  - Legs: E->V; Venus 1:1 full-rev; V-V 1-rev high return; V->E.
- V_inf: E 4.893, V 10.364.
- Turns: Venus 6.57/21.60/15.03 deg (worst 0.574); Earth 31.76 deg (0.347).
- Distance: 0.512-1.243 AU.
- Real ephemeris:
  - Standish direct route: PASS at 5/5 (worst 0.579-0.580; re-fly <= 1.9e-3 km).
  - DE440: no gate-passing ballistic member found (WEAK, old restarts). Closure from the minimax
    solution costs 1.56-2.53 m/s mid-course per 16-yr chain (6.25).
- Prior art: outside Hollister's 3.2-yr FR/SY itineraries. Menning 1968 estimates at least 1024
  Earth-Venus FR/SY orbits and names half-rev and order variations. ev-A carries a generic (non-FR,
  non-SY) Venus return, so it is not one of those. Web prior-art search DONE 2026-10-07
  (`2026-10-07-942-evAB-prior-art-search.md`): no collision; nearest is the Hollister-Menning class.
- Minovitch 1963 (JPL TR 32-464): NO collision. All its returns are one-shot (Table 14 E-V-E free
  returns; Tables 16-17 one-way E-V-M-E chains); the repeating free-fall idea (p.15) is stated but not
  computed (see ev-C).
- **Decision asked:** candidate-novel?
  - Proposed attribution: Hollister 1969 and Hollister & Menning 1970 (method and the E-V
    two-working-body class).
  - Or "a member of the Hollister-Menning class"? (Web search done: no collision.)

### 2.5 ev-B (Earth-Venus, both bend) — Standish 2/5

- Structure: k3|RE/1:1|LE>V/0s|LV>V/1l|RV/3:2|LV>E/0s.
  - Period 1751.833 d. Dates 547.596902, 597.693177, 1778.796557 d.
- V_inf: E 8.012, V 10.932.
- Turns: Venus 23.49/14.23/9.26 deg (worst 0.673); Earth 2 x 21.81 deg (0.375).
- Distance: 0.610-1.686 AU.
- Real ephemeris:
  - Standish direct: PASS at 2/5 epochs (0.933, 0.757; re-fly <= 1.6e-2 km); gate fails at 2;
    no convergence at 1.
  - DE440: Venus legs as ev-A. The Earth 1:1 figures are an upper bound (lunar reflex in DE440; 6.25).
- Prior art: as ev-A (search DONE, no collision). Near: the VESTA web concept's "1752-d Earth-Venus
  cycle" (no sources) has ev-B's period (3 synodic periods). But its E->V leg is 109 d, against ev-B's
  50.1 d: a period coincidence only.
- Minovitch 1963 (JPL TR 32-464): NO collision. All its returns are one-shot (Table 14 E-V-E free
  returns; Tables 16-17 one-way E-V-M-E chains); the repeating free-fall idea (p.15) is stated but not
  computed (see ev-C).
- **Decision asked:** as ev-A, with the weaker real-ephemeris standing.

## 3. Recorded, not proposed as finds

| Item | Ideal model | Real ephemeris | Status |
|---|---|---|---|
| ge-1 k1\|LG>E/0s\|RE/1:1\|LE>G/0s, G 3.734 / E 8.184, worst 0.550 | gate pass | rung (d) NEGATIVE (conditional): closes at 5/5 epochs, all gate-failing (2.57-6.60) | ideal-model member only |
| ge-2 k3\|LG>E/1h\|RE/3:2\|LE>G/0s, G 1.371 / E 1.620, worst 0.612 | gate pass | NEGATIVE (conditional): one 10-cycle closure, 7.93; fails from k = 2 | ideal-model member only |
| ge-3 k3\|LG>G/1h\|LG>E/1h\|RE/1:1\|LE>G/1l, G 3.881 / E 8.413, worst 0.883 | gate pass | NEGATIVE (conditional): 3/5 epochs, 2.88-3.32 | ideal-model member only |
| ge-4/5/6 | near-one-body (turn < 1 deg, ratio < 0.05 at one moon) | not laddered | R-S-class relatives |
| em-1, em-2 (k3, E 5.333 / M 4.713, Mars generic return 2 x 38.4 deg), em-3 (E 4.684 / M 4.539) | worst 0.939-0.983 | rung (d) not passed: converged landings fail at Mars (2.9-24) | ideal-model members only |
| em-4, em-5 (Mars 1:1, M 9.8-10.0) | worst 0.988-0.996 | not laddered | Rall-adjacent, not Rall members |
| vm2-1 k3\|LV>V/2l\|LV>M/0s\|LM>V/0s, V 6.086 / M 4.849, worst 0.981 | passes the registry floor, fails Rall/H&M's 1.1-radius rule | single-cycle failure; the 7-cycle negative is path-limited | ideal-model curiosity (6.10); in Rall's attempted class (thesis sec. 4.4); weak related prior art: the one-shot E-V-M-E chains of Minovitch 1963 (Venus V_inf 5.37-6.11 km/s; no M->V or V->V leg) |
| vm2n-2 k4\|RV/2:1 x 2\|LV>M/0s\|LM>V/0s, V 6.221 / M 4.864, worst 0.421 | vm2-1 neighbour with margin | rung (d) not passed (1/5 converged, Mars 5.91) | ideal-model member only |
| vm-1, the em one-body rows, the ge one-body rows, vmn | | | R-S / R-O class members |

- The "negatives" are conditional on the patched-conic model, the epochs and the methods used. They are
  not proofs of absence.
- The blend's lambda = 1 gate passes for ge were an artefact: the full-rev legs miss the moon by
  748-3,192 km (6.49).
- Prior art for em: Rall 1969 and Rall & Hollister 1971 (no member), Pisarevsky 2008 (outside its
  covered diagrams; Fig. 14 not digitised), the catalogue's 273 E-M rows (no V_inf match).

## 4. Open items

- Minovitch JPL TR 32-464 (1963): full digest in progress (corpus-file-opus); the key-page verdict above is no collision.
- Web prior-art search: done for all (`2026-10-07-942-943-ge-em-prior-art-search.md` for ge and em: no
  collision). To acquire before any em row is proposed: Fornari & Pontani 2020 (AM&S,
  10.1007/s42496-020-00050-6), a global search for cycling E-M families, not held. Pending for
  ev-A/ev-B: Hollister 1963 Sc.D. and Crocco 1956 (corpus-file-opus queue).
- The Jovian n-body lane (#968): no positive control, so gc-1 and gc-2 have no n-body check.
- gc-2: a continuous-gravity check of the GanCal#5 relation (owner option).
- Pisarevsky 2008 Fig. 14 (class I.1 points) is not digitised.
- vm2n-2 and em: a ramp-continuation real-ephemeris attempt (not run; other landings untested).
- The owner's 6.45 amendment and the 6.18/6.24 conventions (DE440 for full-rev rows is descriptive
  only) apply to the readings above.

## 5. Decisions asked of the owner (#875)

| Candidate | Real-ephemeris standing | Literal collision | Proposed `our_status` | Attribution if candidate-novel |
|---|---|---|---|---|
| gc-1 | PASS 5/5 (jup365, 10 cycles) | none | candidate-novel | Hollister & Menning 1970 (method); Russell & Strange 2007/2009 (model, future-work statement); Campagnola et al. 2019 (class) |
| gc-2 | PASS 5/5 (jup365, 10 cycles) | none; near the R-S GanCal family | candidate-novel, or GanCal-family relative | as gc-1, plus R-S GanCal#5 (skeleton) |
| ev-C | PASS 5/5 (Standish and DE440) | none (D. Ross 2L4 nearest, unrefereed; Minovitch 1963 none) | candidate-novel under `#875` (ii) | Russell & Strange 2007/2009, re-applied to Earth-Venus with Venus hosting |
| ev-A | Standish 5/5; DE440 near-ballistic | none (Minovitch 1963 none; web search done) | candidate-novel, or Hollister-Menning class member | Hollister 1969; Hollister & Menning 1970 |
| ev-B | Standish 2/5 | none (Minovitch 1963 none; web search done) | as ev-A | as ev-A |
| ge-1..3, em-1..5, vm2-1, vm2n-2 | not passed / negative (conditional) | none | not candidate-novel (ideal-model members only) | — |
