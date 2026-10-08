# #1025: the literature gate (v3-A1) on the 42 clean #973 members

Status: DONE. Pre-registration secs. 1-2 (commit 05a18d39, before any run); results sec. 3; prior art by
family secs. 4-5. No catalogue writes; nothing in this note is called novel.

Members: the clean two-working-body members of `docs/notes/2026-10-07-973-two-working-body-k-extension.md`.
- gc k = 4: 10 (sec. 10).
- gc k = 5: 6 (sec. 11).
- gc k = 6: 4 (sec. 12).
- ev k = 4: 11 (sec. 13).
- ev k = 5: 11 (sec. 14).

Gate: `literature_check` v3-A1 (`docs/notes/2026-10-08-1045-literature-gate-v3-preregistration.md`
secs. 1-8), offline corpus backend.

## 1. Method

Driver: `scripts/litcheck_1025_973_members.py`. It imports `scripts/litcheck_942_943_scope.py` and uses
its `sig_of`, `from_gauntlet` and `run` unchanged, so each signature is built exactly as for gc-1 and
ev-A:
- primary: Jupiter (gc) or Sun (ev);
- sequence: the time-ordered encounters of one cycle;
- period: k;
- V_inf: per encounter;
- topology: {"repeated-moon"};
- working_bodies: "two" when each body turns >= 0.05 deg somewhere;
- return types: from the cycle key and the leg flight times (FR, FR-n:m, HR, SY, GEN).

The gate runs once per member. In the same run it also re-checks the gate state, with expected values
from the v3-A1 record (`data/942_943_litcheck_scope.json`):
- GanCal#5 and Hollister 1H: "published";
- gc-1, gc-2, ev-C: "inconclusive";
- ev-A, ev-B: "not-found".

The labels were computed first, without the gate (`--labels-only`), so the expectations below can be
specific:
- All 42 are "two".
- No gc member is strictly alternating.
- Every ev member carries at least one FR-n:m return with n:m other than 1:1 (2:1, 3:2 or 2:3), and
  some also a GEN return.

## 2. Expected statuses (fixed before running)

- **gc, 20 members: "inconclusive", naming the Russell & Strange Ganymede-Callisto anchor through F7.**
  - That anchor's body set contains {Ganymede, Callisto}, and it is excluded on working bodies ("one"
    against "two"). This is the gc-1/gc-2 result.
  - The Campagnola 2019 GCGC anchor (two working bodies, `alternating_scope`) is excluded because no gc
    member alternates.
  - **One named risk, gc6-0.** Its time-ordered sequence is (Ganymede, Callisto, Ganymede): two
    distinct Ganymede flybys and one Callisto flyby per cycle, i.e. the cyclic order G G C. The gate's
    #972 F3 rule reads a final encounter equal to the first as a catalogue "closing repeat" and drops
    it. That would make (G, C) look alternating, and the GCGC anchor might then not exclude it.
    - So gc6-0 may come out "published" (GCGC) or with a GCGC-led "inconclusive".
    - If it does, the pre-registered result stands as recorded, and is reported as a
      signature-convention artefact.
    - A diagnostic re-run with the same cycle started at the Callisto encounter, (Callisto, Ganymede,
      Ganymede), is then reported beside it, labelled as a diagnostic. It does not replace the result.
- **ev, 22 members: "not-found".**
  - The H&M anchor (two working bodies, return types {FR, SY, HR}) is excluded on return types. Every
    ev member has an FR-n:m return outside that set, and some also GEN.
  - F7 does not fire: the working labels agree ("two"), so working bodies are not among the exclusion
    reasons. This is the ev-A/ev-B result.
  - No ev member has an HR return.
- **Checks:** as listed in sec. 1, all as expected.

Stop rule: if a check is not as expected, or a member's status differs from the above (other than the
named gc6-0 risk), the run's output is kept, nothing is re-run, and the lead is told.

## 3. Results (one run, 2026-10-08 19:51 AEDT)

Data: `data/1025_litcheck_973_members.json` (each member's status, confidence, citation, notes, and the
excluded anchors with their reasons) and `data/1025_litcheck_973_members.log`.

- **Checks: all as expected.**
  - GanCal#5 and Hollister 1H: published.
  - gc-1, gc-2, ev-C: inconclusive.
  - ev-A, ev-B: not-found.
- **gc, 20 members: 19 "inconclusive" naming Russell & Strange 2009 (Ganymede-Callisto) through F7,
  as pre-registered. gc6-0 is "published" (0.85), citing the Campagnola 2019 GCGC anchor: the named
  risk of sec. 2.**
  - Diagnostic for gc6-0 (`data/1025_litcheck_gc6-0_diagnostic.json`; NOT the result): the same
    signature with the cycle started at the Callisto encounter, (Callisto, Ganymede, Ganymede), gives
    "inconclusive" naming R-S Ganymede-Callisto through F7, like the other 19.
  - Reading: gc6-0 meets Ganymede twice in a row per cycle (G G C), with V_inf G 1.673 / C 1.257
    km/s. The GCGC cycler strictly alternates G C G C at 3.5 / 4.5 km/s.
  - The gate's #972 F3 rule drops a final encounter equal to the first as a catalogue "closing
    repeat". On the generator's time-ordered list of distinct flybys, that turns G G C into an
    alternating (G, C).
  - The "published" result is a signature-convention artefact, not a collision. It is recorded as
    run; the lead decides whether the gate or the signature builder changes (papercut filed).
- **ev, 22 members: all "not-found", as pre-registered.** The H&M anchor is excluded on return types in
  every case (each member has an FR-2:1, FR-3:2 or FR-2:3 Earth or Venus return); F7 does not fire.

| Family | Members | Status | Named anchor |
|---|---|---|---|
| gc k = 4 | 10 | inconclusive (all) | R-S 2009 Ganymede-Callisto (F7: working bodies) |
| gc k = 5 | 6 | inconclusive (all) | the same |
| gc k = 6 | 4 | 3 inconclusive; gc6-0 published (GCGC; convention artefact, diagnostic inconclusive) | R-S 2009 G-C; Campagnola 2019 GCGC for gc6-0 as run |
| ev k = 4 | 11 | not-found (all) | none; H&M excluded on return types |
| ev k = 5 | 11 | not-found (all) | the same |

None of this makes any member novel:
- "inconclusive" goes to a human (the owner);
- "not-found" is necessary, not sufficient (`#875` policy);
- the prior-art reading by family follows.

## 4. Prior art, gc family (20 members, k = 4, 5, 6)

Sources already searched or read:
- `2026-10-06-943-gc-prior-art-search.md`:
  - sec. 2: the only published two-working-body G-C cycler is Campagnola 2019 GCGC (alternating,
    3.5 / 4.5 km/s);
  - sec. 9: Buffington 2014, Clipper 13F7-A21, a one-way G-C pump-down;
  - sec. 10: Niehoff 1971, Io-commensurate capture ellipses, no moon-to-moon legs;
  - sec. 11: six Lynam capture papers, one-shot, Ganymede V_inf >= 4.6 and Callisto >= 7.1 km/s.
  - None collides with gc-1 or gc-2. None of these sources has a closed periodic two-working-body G-C
    cycle at 4-6 synodic periods either.
- `2026-10-07-971-fable-corpus-review-2.md` R11 names the k = 4 (50.09-d) relatives:
  - Boutonnet, Langevin & Erd 2024: the JUICE 50-day cycle, one-shot;
  - Lynam 2015: capture windows 50.07-50.13 d apart;
  - Lam et al. 2018 17F13: Callisto 3:4 and 3:5 legs of 50.0 d;
  - Liang et al. 2024: the C-G-C half of a three-moon CGCEC cycle, at 5.67 / 6.99 km/s, open, not
    closed on itself.
  - The #973 note sec. 3.2 reproduces Liang's half as an open segment in the gc cell. No gc member is
    within 0.3 km/s of Liang's per-moon V_inf (#973 gauntlet).
- Russell & Strange 2009 (Table 1: Ganymede -> Callisto, Callisto massless) is the treated system that
  makes every gc member "inconclusive". Their GanCal rows are k = 3. Our digest does not record the
  largest period they searched (#973 note sec. 4).
- The gc members are periods 4-6 of the same two-working-body architecture as gc-1. Per the #971 note
  (R11, "Policy"), the owner ruled gc-1's class candidate-novel; a k >= 4 member would be another member
  of that class. That is the owner's reading to confirm; this note does not extend it.

Continuous gravity (`docs/notes/2026-10-08-1025-sigma-batch.md` sec. 10, commit 2272c1f8; the R-S
ideal model with both moons massive):
- EXISTS at full mass (15): gc4-1, 2, 6, 8, 10, 25, 28, 38, 41; gc5-3, 5, 6; gc6-0, 1, 2.
- FOLDS (1): gc6-3 (sigma_f 0.748).
- IMPACT (2): gc5-0 (sigma 0.438), gc5-4 (sigma 0.908).
- NUMERICAL STOP (2): gc4-19, gc5-7 (untested, not negative).
- Real ephemeris (jup365, continuous gravity): untested for every member. Per the sigma note sec. 10.4,
  the lane cannot yet pose full-rev legs on jup365 (#1044/#1046).

## 5. Prior art, ev family (22 members, k = 4, 5)

Sources: `2026-10-07-942-evAB-prior-art-search.md` (the citer graphs of Hollister 1969 and H&M 1970,
14 OpenAlex, 6 NTRS and the arXiv queries; secs. 2-3).
- The published Earth-Venus cycler record is Hollister 1969 and H&M 1970 / Menning 1968. Both planets
  work there, with FR (1:1) and SY returns only, plus Menning's named but uncomputed half-rev and
  order variations.
- Every #973 ev member has an n:m full-rev return with n:m other than 1:1 (2:1, 3:2, 2:3), which is
  outside those itineraries. Some also have a generic Venus return (the ev-A/ev-B distinguishing
  feature).
- Structurally, no k = 4 or k = 5 one-visit structure equals an H&M orbit or a sub-period of one (#973
  note sec. 6.5). The nearest H&M row by V_inf is 0.39 km/s (orbit 7, from ev5-0/1) and otherwise
  0.69 km/s or more.
- k = 5 is the 8-yr (8:13) repeat. VanderVeen 1969 names the 8-yr geometric repeat; Ross c.2021's
  "5(1.0)10" is a powered, uncomputed Venus-Earth concept. Both are lineage for the period, not the
  orbit.
- None of the other works in the evAB note's table (Minovitch 1963, Jones 2017 VEM, Pisarevsky 2008,
  Hughes, VESTA, the Earth-Mars cyclers) is an Earth-Venus two-working-body cycler (evAB note sec. 3).
  The VESTA concept's 1,752-d period is k = 3, not 4 or 5. No new search was run for the k = 4/5
  members; the evAB search's queries were not restricted by period.
- ev5-10 is a sun-grazer (r_min 0.058 AU, 12.5 solar radii). It is listed, not credible as a design.
- Continuous gravity and real ephemeris: untested for all 22. #1007's ev positive control (Hollister
  1H, Standish ramp) failed at 5/5 epochs (commit 3c30d7d7), so the heliocentric chain lane has no
  validated ev route yet.
