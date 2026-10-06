# #942 R1(b) candidate ev-C: prior-art search (2026-10-06)

Requested by the lead. The candidate is from `docs/notes/2026-10-05-942-943-two-working-body-generator.md`
secs. 6.7 and 6.14; gauntlet `data/942_cell_ev_gauntlet.json`. Search by corpus-file-opus, with the same
method as `2026-10-06-943-gc-prior-art-search.md`.

**ev-C:**
- k = 2, period 1167.89 d (2 Earth-Venus synodic periods).
- Key `k2|LE>V/0s|LV>V/1h|LV>E/0s`. Legs from the gauntlet x_days:
  - E->V zero-rev, 122.83 d;
  - V->V one-rev high return, 823.95 d;
  - V->E zero-rev, 221.11 d.
- v_inf E / V = 9.07 / 13.17 km/s. Venus turns 19.65 deg twice (ratio 0.746).
- Earth turn ratio 0.06-0.11 on the real-element chain: Earth is almost massless. This is the R-S
  architecture at Venus, with Venus hosting the returns.
- Rung (d): ballistic over 5 cycles (16 yr) on Standish mean elements.

## 1. Method

- **OpenAlex citer lists:**
  - Hollister & Menning 1970 (10.2514/3.30134): 15 citers.
  - Hollister 1969 (10.2514/3.29664): 45 citers, about 17 of them AIAA book chapters.
  - Jones, Hernandez & Jesick 2017 (AAS 17-577, OpenAlex W2953427882): 1 citer.
  - R-S 2007/2009: reused from the gc search.
- **Topic queries** (OpenAlex and web): "Earth-Venus cycler", "Venus cycler trajectory", "Earth Venus
  periodic orbit spacecraft", "Venus free return trajectory resonant", "Venus gravity assist cycler
  Earth", "VISIT orbit Venus", "Earth Venus Mars cycler", Venus-hosted generic/free-return cyclers, and
  Byrnes/McConaghy/Longuski Venus variants.
- **Full or abstract reads:** Wikipedia "Cycler" (Earth-Venus content cites only Hollister 1969);
  venautics.space VESTA "2.4 Cyclers" (habitat concept, a 1752-d Earth-Venus cycle, no sources); arXiv
  2006.04900 (Venus human-science white paper, no cyclers).
- **Held status:** every hit was run through `scripts/check_wanted_vs_corpus.py`.

## 2. Verdict

**No refereed publication of an Earth-Venus cycler with Venus hosting the returns was found.** The
refereed Earth-Venus cycler record is still Hollister 1969 and Hollister & Menning 1970 (Earth and Venus
both working). R-S 2007/2009 ran no Earth-Venus set, and the Purdue Venus work (Hughes et al.) is one-shot
free returns.

**One possible class collision (unrefereed), to acquire and adjudicate:**
- D. Ross, "A cycler-quartet between Venus and Sol/Terra L1", manuscript on academia.edu (item 45489864;
  undated; not peer-reviewed as far as visible). I could not read it: academia.edu returns 403 to the
  fetch tool.
- From the search snippets of its abstract:
  - "a new cycler between Venus and Earth ... denoted '2L4', which returns to its planet after two
    synodic periods";
  - it "passes Venus and may use its well; it falls short of Earth proper but approaches its first
    Lagrangian halo" (Sun-Earth L1);
  - "From Venus, the Earthbound trip time is 159 days, then 1008 days downtime before return to Venus";
  - "the Earth-to-Venus cycler taking 123 days, launching ten days after Venus-to-Earth";
  - a quartet of four cyclers (outbound and inbound twins, each a synodic period apart).
- Comparison with ev-C, INFERRED from the snippets only:
  - same class: Venus hosts, Earth is not a working body, k = 2. The 2L4 period is 159 + 1008 = 1167 d,
    against 1167.89 d for ev-C;
  - the 2L4 Earth-to-Venus leg (123 d) is within 0.2 d of ev-C's E->V leg (122.83 d);
  - the Venus-to-Earth legs differ: 159 d for 2L4 against 221.11 d for ev-C;
  - 2L4 stops at Sun-Earth L1 rather than Earth, and no v_inf is quoted in the snippets.
  - So 2L4 may be ev-C's family, or a neighbour with a different V->E branch. It cannot be judged without
    the manuscript.

## 3. Findings

| # | Work | Held? | Earth-Venus content | Collides with ev-C? |
|---|---|---|---|---|
| 1 | Hollister 1969, JSR 6(4):366 (10.2514/3.29664); Hollister & Menning 1970, JSR 7(10):1193 (10.2514/3.30134); Menning 1968 thesis | held | the Earth-Venus cycler families with both bodies working (the in-run control) | no: ev-C is Earth-massless, outside the H&M itineraries (sec. 6.7) |
| 2 | Ross, D., "A cycler-quartet between Venus and Sol/Terra L1", academia.edu 45489864 | not held | "2L4": a Venus-hosted, 2-synodic Venus-"Earth" (L1) cycler; legs 159 d and 123 d | POSSIBLE class collision (sec. 2); acquire |
| 3 | Russell & Strange 2007 (AAS 07-118) / 2009 (JGCD 32(1)) | held | the one-body architecture; VenMar#45 (Venus hosts, Mars massless); no Earth-Venus set | no; the R1(b) gate is OPEN against them |
| 4 | Jones, Hernandez & Jesick 2017 (AAS 17-577) | held | VEM triple cyclers at low v_inf (seeds below 5 km/s at Earth and Mars) | no (v_inf far from 9.07/13.17) |
| 5 | Hughes, Edelman, Saikia & Longuski 2015, "Fast Free Returns to Mars and Venus with Applications to Inspiration Mars", JSR 52(6):1712-1735 (10.2514/1.A33293) | not held (the 2014 AIAA 2014-4109 version is held) | one-shot Earth-Venus-Earth free returns for a human Venus flyby | no (not periodic; Earth hosts) |
| 6 | Hughes 2016, "Gravity-Assist Trajectories to Venus, Mars, and the Ice Giants", PhD dissertation, Purdue (AAI10248372) | not held | the same Venus free returns, extended | unlikely; acquire for completeness |
| 7 | Gillespie & Ross 1967; Sohn 1964; VanderVeen 1969; Hollister & Prussing 1965 | held | one-shot Venus swing-bys and E-V-M-V-E triples | no |
| 8 | Sanchez Net et al. 2022 (JSR); Pisarevsky 2008; Rall 1969/1971 | held | Earth-Mars | no |

The other OpenAlex citers of H&M 1970 and Hollister 1969 are Earth-Mars cycler, Jovian or textbook
items; none has an Earth-Venus cycler.

## 4. Recommendations

- **Acquire** the Ross "cycler-quartet" manuscript (academia.edu 45489864; needs a logged-in browser
  download by the owner). Adjudicate 2L4 against ev-C on: the period, the E->V and V->E leg times, the
  v_inf at Venus and at Earth or L1, and whether Earth (or L1) is approached ballistically with zero turn.
- Until then, record ev-C as "candidate, pending owner adjudication, NOT novel", with the 2L4 manuscript
  as a possible class precedent (unrefereed).
- **Acquire** Hughes et al. 2015 JSR and Hughes 2016 PhD for completeness (low collision risk).
