# #942 R1(b) candidates ev-A and ev-B: prior-art search (2026-10-07)

Requested by the lead. The method follows `2026-10-06-943-gc-prior-art-search.md` and
`2026-10-06-942-evC-prior-art-search.md`. Search by twobody-gen2-opus. The candidates are in
`docs/notes/2026-10-06-942-943-owner-decision-summary.md` sec. 2.4-2.5.

- ev-A: k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s, period 1167.89 d (2 E-V synodic periods). V_inf E 4.893 /
  V 10.364. Both planets bend: Venus 6.6/21.6/15.0 deg, Earth 31.8 deg.
- ev-B: k3|RE/1:1|LE>V/0s|LV>V/1l|RV/3:2|LV>E/0s, period 1751.83 d (3 synodic periods). V_inf E 8.012 /
  V 10.932. Venus 23.5/14.2/9.3 deg, Earth 2 x 21.8 deg; E->V leg 50.1 d.
- Both carry a GENERIC (non-full-rev, non-symmetric) Venus return. That is what puts them outside the
  Hollister-Menning FR/SY itineraries (results note sec. 6.7, 6.12).

## 1. Method

- **OpenAlex citer lists (re-run 2026-10-07):**
  - Hollister & Menning 1970 (10.2514/3.30134): 15 citers.
  - Hollister 1969 (10.2514/3.29664): 45 citers; 27 have indexed titles in the API page.
  - Pisarevsky, Kogan & Guelman 2008 (10.2514/1.30046): 2 citers.
  - Titles screened for Venus, cycler, periodic, swing-by, free return; abstracts read where relevant.
  - Menning 1968 (MIT S.M. thesis, "Freefall Periodic Orbits Connecting Earth and Venus") has no
    OpenAlex record, so its citers come only through the two papers above.
- **OpenAlex title/abstract searches (14 queries):** "venus cycler", "earth venus cycler trajectory",
  "earth venus periodic orbit spacecraft", "venus free return periodic", "venus gravity assist periodic
  orbit", "two planet periodic trajectory", "planetary cycler venus", "venus earth cyclic trajectory",
  "venus resonant return trajectory", "earth venus repeating trajectory", "venus swingby periodic orbit",
  "cycler orbit venus habitat", "earth-venus shuttle orbit", "interplanetary cycler venus mars earth".
- **NTRS API (6 queries):** "venus cycler", "earth venus periodic orbit", "periodic swingby orbit venus",
  "venus earth cycling spacecraft", "venus free return repeated", "earth venus shuttle cycling orbit".
- **arXiv API:** ti:cycler (3 hits); abs:cycler AND abs:venus (0); abs:cycler AND earth AND mars (0).
- **Web searches (6):** "Earth-Venus cycler", Venus cycler gravity assist, Hollister-Menning follow-ups,
  "Venus cycler" repeated encounters, E-V periodic orbits at 2-3 synodic periods, Menning 1968.
- Google Scholar is not reachable by the tools; its role is covered by the OpenAlex citer graph.
- Raw outputs: scratch only (not committed).

## 2. Verdict

**No publication of an Earth-Venus cycler with a generic Venus return was found, and no collision with
ev-A or ev-B.** The published Earth-Venus cycler record is still:
- Hollister 1969;
- Hollister & Menning 1970 (Menning 1968): both planets working, with full-revolution (FR) and
  symmetric (SY) returns only. Menning estimates at least 1024 such orbits and names half-rev and order
  variations without computing them. ev-A and ev-B are neither (generic Venus returns).
- Jones, Hernandez & Jesick 2017: VEM triple cyclers (three planets, low V_inf).

## 3. Findings

| # | Work | Earth-Venus content | Collides with ev-A / ev-B? |
|---|---|---|---|
| 1 | Hollister 1969, JSR 6(4):366; Hollister & Menning 1970, JSR 7(10):1193; Menning 1968 thesis (held) | the E-V two-working-body cycler class, FR/SY returns | no: ev-A/B use generic Venus returns, outside the FR/SY itineraries; nearest published class |
| 2 | Minovitch 1963, JPL TR 32-464 (held; results note 6.52) | one-shot E-V-E free returns; the "repeats the same flight" idea, not computed | no |
| 3 | Jones, Hernandez & Jesick 2017 (held) | VEM triple cyclers | no (three planets; V_inf far) |
| 4 | Pisarevsky, Kogan & Guelman 2008 and IAC-06 (held) | the two-planet periodic method; no E-V numbers | no |
| 5 | Hughes et al. 2014/2015 (held); Hughes 2016 PhD (not held) | one-shot E-V-E free returns (Inspiration Mars) | no |
| 6 | D. Ross 2021 "cycler-quartet" (held) | powered Venus-to-L1 "2L4", k = 2 | no (ev-C's nearest relative; not E-V two-working-body) |
| 7 | VESTA "2.4 Cyclers" (venautics.space, web concept, no sources) | an E-V cycle of 1752 d (3 synodic periods), 5 revolutions, E->V 109 d; no V_inf, no flyby data | no: ev-B has the same period but an E->V leg of 50.1 d, not 109 d; a period coincidence only |
| 8 | Sanchez Net et al. 2020/2022 (Pony Express); Rogers 2014 PhD; McConaghy 2004 PhD; Russell 2004/2005 | Earth-Mars cyclers only | no |
| 9 | selenianboondocks "EVMVE-2034" (blog) | a one-shot E-V-M-V-E circumnavigation | no |
| 10 | Sun-Venus CR3BP part 2 (2024, Arch. Appl. Mech.) | periodic orbits around Venus in the CR3BP | no (other model and object) |

The other citers of H&M 1970 and Hollister 1969 are Earth-Mars or Jovian cycler papers or methods
papers. None has an Earth-Venus cycler.

## 4. Pending (corpus-file-opus queue)

- Hollister 1963 Sc.D. thesis (Venus-swingby Mars round trips).
- Crocco 1956 (E-M-V-E).
Their verdicts will be added here when corpus-file-opus reports.

## 5. Recommendation

ev-A and ev-B have no literal collision. Under `#875` the choice for the owner is "candidate-novel,
attributing Hollister 1969 and Hollister & Menning 1970 (method and class; Minovitch 1963 for the
idea)" or "a member of the Hollister-Menning class". The generic Venus return is the distinguishing
feature, and Menning's own list of variations (half-rev, order) does not include it.
