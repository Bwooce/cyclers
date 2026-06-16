# Phase 4 of #289 — precursor MGA matcher for Aldrin / S1L1 (#302)

Phase 4 of task #289 (Track-A epoch-locked trajectory substrate) closes the
matcher's discovery probe: given a steady-state cycler row in
`data/catalogue.yaml`, what one-shot Earth-launched MGA chains could insert
a spacecraft INTO the cycler at its first encounter?

Published cycler papers (Aldrin 1985, Russell-Ocampo 2004, McConaghy 2002 /
2006, Braik-Ross 2026, Roberts-Tsoukkas-Ross 2026, Spreen 2020) document the
*cycler* invariants in detail but only rarely the *precursor* insertion
trajectory.  Phase 4 was scoped on the hypothesis that the insertion region
is a low-corpus-coverage region of the literature where genuine
literature-fresh candidates are plausible.  This note records the honest
verdict.

## Parts

| Part | Deliverable | Result |
|---|---|---|
| A | `cyclerfinder.search.precursor_matcher` module + `tests/search/test_precursor_matcher.py` | 14 / 14 tests pass; 68 / 68 across touched surfaces |
| B | `scripts/run_302_aldrin_precursor.py` + `data/precursor_302_aldrin.jsonl` | 394 / 400 candidates survive closure; best closure 0.45 km/s; **0 literature-fresh per KNOWN_CORPUS** |
| C | `scripts/run_302_s1l1_precursor.py` + `data/precursor_302_s1l1.jsonl` | 394 / 400 candidates; best closure 2.07 km/s; **0 literature-fresh** |
| D | This doc | Verdict + Phase-5 path |

## The matcher in one paragraph

`find_cycler_precursors(cycler_id, catalogue, ephemeris, ...)` reads the
target row's first encounter body + seed V_inf
(`vinf_kms_at_encounters[0]`), enumerates Earth-launched chains via the
Phase-2 Tisserand-Poincaré multi-shell BFS filtered to terminate at the
target body with V_inf bin close to the seed, validates each chain via the
Phase-1 closure driver with Phase-3 TOF optimisation, and runs the Phase-2
`literature_check` against `KNOWN_CORPUS` on each survivor.  The returned
`PrecursorMatch` list is ranked by
`quality_score = vinf_match_residual + closure_residual + 2 * flyby_continuity`
(matching the `optimise_chain_tofs` loss weighting).

The matcher is READ-ONLY on the Phase 1-3 modules — `epoch_aware_genome.py`
and `tisserand_mga_window.py` are unchanged.  The matcher is also
READ-ONLY on `catalogue.yaml` — no writeback, even for literature-fresh
candidates.  Per the discipline preamble in
`cyclerfinder.search.literature_check`, "not-found" is
necessary-not-sufficient for novelty; a fresh candidate is a CANDIDATE for
the V0-V5 gauntlet (when extended to epoch-locked classes; future task),
not a novelty claim.

## Aldrin precursor scan

Target row: `aldrin-classic-em-k1-outbound`, Earth V_inf seed 6.5 km/s
(Russell 2004 dissertation Table 3.4 cycler 1.0.1.-1).  Launch window
2030-01-01 to 2034-12-31 (~2.4 Earth-Mars synodic periods).  Intermediate
bodies: Venus + Earth.  V_inf grid (4, 5, 6, 7, 8) km/s.  TOF box (80, 500)
days per leg.  Max chain length 3 legs.  Multi-shell BFS on; pump-envelope
factor 1.0; a-range (0.3, 2.5) AU.

### Numerical verdict

| Statistic | Value (km/s) |
|---|---|
| Candidates surfaced by BFS | 400 (geometric proposals matching `\|V_inf - 6.5\| <= 0.8`) |
| Survivors of closure + TOF opt + wide gates (100 km/s) | 394 |
| Best closure residual | 0.45 |
| Best V_inf match residual | 0.21 |
| Closure residual median | 12.51 |
| Closure residual p75 | 21.99 |
| Flyby continuity median | 17.48 |
| Flyby continuity p75 | 27.69 |
| Closure < 1 km/s count | 2 |
| Closure < 1 AND flyby < 1 km/s count | **0** |

The closure-residual distribution shows that ~ half of all geometric BFS
proposals achieve closure within 12.5 km/s at the optimised launch epoch +
per-leg TOFs.  Only 2 candidates close inside 1 km/s — the published-grade
threshold — and ZERO satisfy both publication-grade closure AND
publication-grade ballistic continuity simultaneously.

### Chain topology

The surviving 394 partition into two sequences:

  * `E -> V -> E` : 304 candidates (the Venus-flyby precursor archetype)
  * `E -> E -> V -> E` : 90 candidates (an extra Earth resonant return)

The E-V-E pattern is the canonical Venus-leveraging precursor used by
Galileo VEEGA, Cassini VVEJGA, MESSENGER, BepiColombo, and the published
VEM-cycler insertion literature (Jones-Hernandez-Jesick AAS 17-577, Hughes
et al. 2014, Genova-Aldrin 2015).

### Literature-check verdict

**0 / 394 candidates are literature-fresh per KNOWN_CORPUS.**

Every surviving candidate's structural fingerprint
(`primary=Sun, sequence=(E,V,E,…)`, V_inf regime 4-8 km/s) was flagged by
the offline corpus matcher as covered by the published VEM cycler family
(citation: Jones, Hernandez & Jesick, AAS 17-577 "VEM triple cyclers");
several anchors also reach: Hughes et al. 2014 VEM extensions, Strange-Russell
2007 Tisserand pump-tour graph, and the Aldrin / Russell-Ocampo / McConaghy
cycler corpus itself.  The matcher's `literature_check` returns `published`
with confidence 0.95 for all 394 candidates.

### Best survivor

  * sequence: `('E', 'V', 'E')`
  * launch_epoch_utc: `2034-05-25T…Z`
  * leg_tofs_days: ~173, ~190
  * Earth V_inf at terminal: 7.30 km/s (vs seed 6.5; residual 0.80 km/s)
  * closure_residual: 0.45 km/s
  * flyby_continuity (Venus flyby ΔV): 6.58 km/s — NOT ballistic
  * literature: `published` (matched to VEM cycler corpus)

The 6.58 km/s flyby continuity ΔV at Venus means this candidate is not a
true ballistic precursor — it requires a powered Venus flyby (or an
appropriately-placed DSM on the leg into Venus) to actually close.

## S1L1 precursor scan

Target row: `s1l1-2syn-em-cpom`, Earth V_inf seed 5.65 km/s (per spec.md
§9, with PROVENANCE CAVEAT on the catalogue data_gaps).  Same launch window
and search box as Aldrin; V_inf grid (3, 4, 5, 6, 7) km/s.

### Numerical verdict

| Statistic | Value (km/s) |
|---|---|
| Candidates surfaced by BFS | 400 |
| Survivors of closure + TOF opt | 394 |
| Best closure residual | 1.79 (with TOF opt) — note: V_inf match 1.81 at the best closure |
| Best V_inf match residual | 0.06 (different candidate) |
| Closure residual median | 13.51 |
| Closure residual p75 | 22.99 |
| Flyby continuity median | 17.59 |
| Flyby continuity p75 | 27.69 |
| Closure < 1 km/s count | 0 |
| Closure < 1 AND flyby < 1 km/s count | **0** |

S1L1's distribution is structurally similar to Aldrin's but uniformly
slightly worse (closure median 13.5 vs 12.5, best closure 1.79 vs 0.45 km/s).
The same E-V-E (304) / E-E-V-E (90) topology partition.

### Literature-check verdict

**0 / 394 candidates are literature-fresh per KNOWN_CORPUS.**  Same
mechanism as Aldrin — every survivor's E-V-E structural footprint matches
the VEM cycler family.

## The honest scientific verdict

**Zero literature-fresh precursor candidates were surfaced by this Phase 4
discovery probe.**  Both Aldrin and S1L1 precursor searches yielded only
candidates that the offline literature corpus correctly identifies as
covered by the published VEM-cycler precursor literature — primarily
Jones-Hernandez-Jesick 2017 AAS 17-577 ("VEM triple cyclers"),
Hughes-Edelman-Longuski 2014 VEM extensions, and the Strange-Russell 2007
Tisserand pump-tour graph.

This is a CLEAN NEGATIVE for the Phase 4 scope.  It is not a failure of
the matcher (which surfaced 800 candidates, closes 788 of them to a 100
km/s wide gate, and correctly invokes the literature check on every
survivor); it is the substantive answer to the question Phase 4 was scoped
to ask.

### Search-coverage caveat

The literature-fresh count of zero is conditional on:

  * **chain topology**: chains of length 2-3 legs starting at Earth and
    using {V, E} as intermediate bodies.  Mars-as-intermediate (M as a
    flyby on the way to the cycler), low-thrust insertion, or
    asteroid-leveraging tours are out-of-scope for this probe.
  * **launch window**: 2030-2034 (~2.4 Earth-Mars synodic periods).  The
    surviving 394 candidates concentrate at 2034-05-25 to 2034-05-29 —
    a single ~5-day window where Earth-Venus phasing closes.  An
    expanded window (say 2030-2050, ~10 synodic periods) would surface
    additional launch opportunities but each would be variations on the
    same E-V-E topology against the same VEM corpus.
  * **V_inf grid**: 3-8 km/s in 1 km/s bins.  The bin width sets the
    geometric tolerance the BFS uses to bridge between hetero flybys;
    finer bins would catch sub-bin candidates but at substantial BFS-cost
    quadratic in bin count.
  * **closure model**: ballistic Lambert per leg, no automated DSM
    placement, no multi-rev Lambert (n_revs=0).  A multi-rev Lambert
    would surface 1-rev and 2-rev E-V transfers (typical Mariner /
    Venus-Express geometry); automated DSMs would close the Venus
    flyby's 6.58 km/s continuity deficit on the best Aldrin candidate.
  * **insertion target**: the cycler's FIRST encounter (the natural
    "boarding" point).  A precursor that joins the cycler at a *later*
    encounter (post-Mars-flyby reentry into the steady state) is a
    different shape that would need a separate matcher.
  * **literature corpus**: `KNOWN_CORPUS` rev `568d8a4` — the curated 30
    anchor entries.  Live WebSearch would expand the corpus and could
    flag *additional* rediscoveries but is unlikely to surface a
    not-found candidate the curated corpus missed (E-V-E precursors are
    well-published).

The "0 literature-fresh" result is the honest answer to the question Phase
4 was scoped to ask within these bounds.  A genuinely novel precursor
trajectory would require either a topologically different chain (e.g.
asteroid-leveraging, low-thrust, three-body capture) or an unexplored
insertion target.

## Phase 5 path (recommendation)

The matcher's closure quality is gated by two structural deficiencies the
current substrate cannot bridge:

  1. **Automated DSM placement.**  The best Aldrin candidate's Venus flyby
     leaves 6.58 km/s of continuity deficit — far above the ballistic
     gate.  A DSM placed on the leg INTO Venus (per Vasile-Conway 2006
     MGA-DSM transcription) can close that deficit at a fraction of the
     ΔV cost of a powered flyby.  Phase 3's `DSMSpec` extension supports
     a hand-placed DSM but the placement itself is a Phase-5 task.
  2. **Multi-rev Lambert.**  The BFS proposes n_revs=0 chains;
     Earth-Venus transfers historically span n_revs=1 or n_revs=2 (the
     "fast" vs "slow" branches that mission designers swap between
     depending on the launch window).  A multi-rev sweep at the Lambert
     level would expand the surviving population without expanding the
     V_inf or epoch grid.

Phase 5 (task #303 when scoped) should:

  * Automate DSM placement per Vasile-Conway 2006 §3.2 transcription
    (free DSM position + magnitude as continuous decision variables;
    NLOpt or scipy.optimize.differential_evolution as the outer loop).
  * Add multi-rev Lambert to the closure driver — the `lambert(...)` API
    already supports `max_revs > 0`; the Phase-1 closure currently uses
    `max_revs=0`.
  * Extend the V_inf-leveraging-leg substrate so the closure can target
    a SPECIFIC cycler-cadence phase at the terminal Earth encounter
    (today the matcher reports `epoch_alignment_score=0` because no
    target phase window is supplied by the catalogue rows).
  * Re-run the Aldrin + S1L1 scans through the upgraded substrate.
    Expected outcome: closure quality drops to publication grade
    (≤1 km/s), flyby continuity to ≤0.1 km/s, and the literature-check
    correctly flags the upgraded candidates against the published
    DSM-equipped VEM literature (still likely 0 / N novel).

If after Phase 5 the count is still zero literature-fresh, the structural
conclusion is that **precursor MGA insertion into the known Earth-Mars
cycler families is fully published** and the literature-fresh discovery
region must be sought elsewhere — for instance asteroid-leveraging
precursors (the #226 FBS substrate), low-thrust insertion (the #168
low-thrust cycler maintenance substrate extended to insertion), or
non-Earth-Mars-cycler insertion (the Russell-Ross stable cycler families
or the Hernandez-Jones-Jesick Io-Europa-Ganymede triple cyclers).

## Discipline checklist

  * **NO catalogue writeback** ✓ — both JSONLs are discovery-probe output;
    catalogue.yaml unchanged
  * **NO novelty claim** ✓ — the doc explicitly reports "0 literature-fresh"
    and frames each candidate as covered by published literature
  * **literature_check ran on every survivor** ✓ — 788 / 788 (Aldrin 394 +
    S1L1 394) returned a structured verdict; per
    `feedback_literature_novelty_check_baseline` the not-found gate is
    necessary-not-sufficient (and was never reached anyway)
  * **independent cross-check** ✓ — inherited from `close_epoch_locked`;
    `independent_check_residual_kms` recorded on every JSONL row (best
    Aldrin: 3.7e-12 — round-off floor, as expected for the Lambert /
    universal-Kepler agreement on the SAME conic)
  * **READ-ONLY on Phase 1-3 modules** ✓ — `epoch_aware_genome.py` and
    `tisserand_mga_window.py` are unchanged on this branch
  * **READ-ONLY on catalogue.yaml** ✓ — only `by_id[cycler_id]` accessed
  * **pre-commit tests + ruff + mypy passed** ✓
  * **incremental progress reports** ✓ — every 50 candidates the matcher
    emitted a stdout line with running survivor count and ETA
  * **respectful errata framing** — N/A, no defect claims raised
  * **work on main, never branch** ✓
