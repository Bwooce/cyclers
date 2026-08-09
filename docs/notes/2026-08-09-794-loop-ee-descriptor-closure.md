# #794 — loop-ee closure via primary-source Russell leg descriptors: full sweep, 14/14 rows

**Date:** 2026-08-09
**Task:** `#794` (registered by `#793`, dispatched 2026-08-09) — close the remaining
`loop-ee`/`loop-ee-N` `#54` `data_gaps` on the 14 catalogue rows carrying
`free_return_arcs[]` descriptors, with the descriptor-to-segment mapping resolved from
PRIMARY sources per row (not the flagged-approximate `descriptor_to_phsi` structural map).

## Outcome: all 14 rows closed (23 loop segments), every value cross-checked two ways

The blocker `#793` flagged — "no published crosswalk between the per-arc descriptor and
the materialized segments" — dissolved on actually reading the two primary sources:

1. **Russell 2004 dissertation, §4.8 (pp. 125–127)**: "The transit times and Mars v∞ are
   calculated using the designated transit leg, as indicated by an uppercase descriptor
   letter... All legs reported are Earth-Earth free-returns." The Table 4.9–4.13 rows print
   the complete per-leg decomposition **in itinerary order** (Leg 1, Leg 2, ...). So for
   every `russell-ch4-*` row: the uppercase leg IS `out-em`+`ret-me` (split by the Mars
   encounter), and the lowercase legs, in cyclic order after it, ARE the `loop-ee-N`
   segments. Nothing structural needs to be guessed.
2. **McConaghy, Russell & Longuski 2005** (JSR 42(4), DOI 10.2514/1.8123 — the
   nomenclature paper the thesis's own p. 126 cites as the origin of the descriptor
   syntax, "Ref. 25"): gives the complete, closed-form semantics of every descriptor
   argument — g(t_f, θ, ε) Lagrange-formulation Lambert (their Eqs (1)–(12)),
   f(M:N, ϕ, λ) resonant full-rev with a = a_E(M/N)^(2/3) and the v_out direction from
   (ϕ, λ) (Eqs (14)–(16)), h(t_f, N, ε, i′) backflip (Eqs (17)–(19)). For the two
   `russell-ocampo-*` rows, its **Fig. 4 explicitly labels Cycler 2.5.1.+0's four legs**
   (Leg 1 g designated w/ Mars encounter; Leg 2 f(1:1,74.9,−144.1); Leg 3 h(1/2,0,U,15.1);
   Leg 4 f(1:1,74.9,35.9)) and **Table 2 prints Cycler 4.3.1.−5's full label** — the
   per-row arc-option mapping the dispatch asked for is *printed in the source*.

`#793`'s framing was thus slightly off in a useful way: the descriptor list is NOT a list
of "candidate realizations" for these rows — it is the itinerary, leg by leg. (What IS
approximate remains `cycler_assembly.py::descriptor_to_phsi`'s map onto Russell's p.h.s.i
integers — untouched by this task and still flagged there.)

## Method (per arc type)

All solves in canonical units (r_E = 1 AU, μ_sun = 1.32712440018e11 km³/s²,
v_E = 29.7847 km/s, year = Earth's period 365.2569 d). Scratch scripts:
`solve_794.py`, `propagate_ocampo.py`, `lambert_crosscheck_794.py` (session scratchpad;
results fully recorded here and in the row notes).

- **g(t_f, θ, ε)** — generic return: solve the paper's Eq (1)
  √μ t_f = a^{3/2}[2πN + α − β − sinα + sinβ] for a on the ε branch
  (U = upper α = 2π−α₀; L/Ls/Ll = lower; Ls/Ll = short-/long-period root when the lower
  branch has two), N = floor(θ/360°); v_out from Eqs (9)–(11); (a, e) from (r, v_out).
  No 0°/180° degeneracy: θ is printed, so the transfer geometry is never ambiguous —
  this is exactly what sidesteps `#793`'s near-resonant Lambert wall.
- **f(M:N, ϕ, λ)** — full-rev resonant: a = (M/N)^(2/3) exactly (M = years = ToF,
  N = spacecraft revs — see the descriptor.py bug below); v_out = √(μ(2 − 1/a)) directed
  by (ϕ, λ) per Eq (16); λ = 0/±180 ⇒ in-ecliptic (e only), otherwise the loop is
  genuinely inclined (e AND i; i recorded in segment notes — the segment schema has no
  inclination field).
- **h(t_f, N, ε, i′)** — half-rev backflip: a from the 180° Lambert (t_f = 0.5 yr ⇒ the
  circular a = 1 boundary solution), v_out per Eqs (17)–(19): a = 1, e = 0, i = i′.

## Cross-checks (all passed)

1. **Independent published V∞ match, 23/23 arcs**: every derived loop's emergent Earth
   V-infinity matches the row's separately-printed v∞E column to ≤ 0.02 km/s (print
   precision). The solve consumes only (ToF, angle) or (M:N, ϕ, λ); v∞E is never an
   input for g arcs and only enters f/h arcs through the printed ϕ/i′ — so this is a
   genuine consistency check on the mapping AND the solver. It also uniquely selects the
   correct root in both double-root (Ls/Ll) cases, agreeing with the printed subscript.
2. **Second solver, 13/13 g-arcs**: `cyclerfinder.core.lambert.lambert(max_revs=8)` on
   each g-arc's geometry reproduces (a, e) to 4 decimals and supplies the project-
   vocabulary (n_revs, branch) labels stored in the segments.
3. **Propagation of Russell's own printed reproduction vectors (ocampo rows)**:
   Kepler-propagating Table 3.5 / Table 3.7's printed initial-v∞ + per-flyby Δv vectors
   through the full cycle reproduces every loop arc at print precision:
   - 2.5.1+0: 366 d arc → a=1.0003, e=0.2118, i=9.005° (derived 1/0.2107/8.99);
     182 d arc → a=0.9996, e=0.0010, i=15.077° (derived 1/0/15.081);
     365 d arc → a=1.0018, e=0.2115, i=8.975° (derived 1/0.2107/8.99).
     This *independently confirms the leg ORDER* (f, h, f) printed in Fig. 4.
   - 4.3.1−5: the whole 548 d loop emerges as ONE conic, a=0.9998, e=0.0006, i=5.954°,
     re-passing Earth's position at the 1-yr mark to 0.006 AU with zero turn — confirming
     that its f(1:1,84.039,∓90) and h(0.5,0,U,±5.961) legs are the SAME circular inclined
     orbit (λ=±90 at a=1 ⇒ e=0, i=5.961°) continuing through an unpowered Earth passage.
     This resolves why Russell lists only 2 flybys and Table 3.7 no encounter at ~2948 d.
   - Bonus: the propagation also emerges the designated g-arc conics (2.5.1+0:
     a=1.5651/e=0.4010, independently matching `#596`'s AR+V∞ inversion values
     1.5633/0.4001 on `out-em`; 4.3.1−5: a=1.2524/e=0.2049, aphelion 1.509 AU ↔ printed
     AR 0.99×1.52).

## What was written back (data/catalogue.yaml)

Per-row loop segments now carry tof_days / n_revs / branch / a_au / e (+ inclination in
the note where nonzero), with DERIVE provenance and the arc-option mapping cited:

| row | loop segment(s) | source leg | a (AU) | e | i (deg) | emergent V∞ (printed) |
|---|---|---|---|---|---|---|
| russell-ch4-4.991gG2 | loop-ee 533.7 d, n1 low | g(1.4612,526.02,Ll) | 1.0512 | 0.1704 | 0 | 4.995 (4.99) |
| mcconaghy-2006-em-k2 | loop-ee (same, model caveat) | g(1.4612,526.02,Ll) | 1.0512 | 0.1704 | 0 | 4.995 (row anchor 4.7 — see gap) |
| russell-ch4-8.049gGf2 | loop-ee-1 365.3 d n1 resonant; loop-ee-2 546.1 d n1 low | f(1:1,74.468,−180); g(1.4951,538.24,Ll) | 1.0000; 1.0837 | 0.2678; 0.2708 | 0; 0 | 8.049/8.050 (8.05) |
| russell-ch4-9.353Gg2 | loop-ee 930.3 d n2 low | g(2.5469,916.9,L) | 1.0555 | 0.3093 | 0 | 9.352 (9.35) |
| russell-ch4-3.64gGg3 | loop-ee-1 365.3 d n1 resonant; loop-ee-2 907.5 d n2 low | f(1:1,82.995,−180); g(2.4845,894.42,Ll) | 1.0000; 1.0213 | 0.1220; 0.1227 | 0; 0 | 3.639/3.640 (3.64) |
| russell-ch4-3.78Gg3 | loop-ee 1279.1 d n3 low (RESTRUCTURED 2→1) | g(3.5018,1260.65,L) | 1.0157 | 0.1270 | 0 | 3.783 (3.78) |
| russell-ch4-5.30gGf3 | loop-ee-1 1095.8 d n2 resonant; loop-ee-2 535.0 d n1 low | f(3:2,82.487,118.851); g(1.4646,527.25,Ll) | 1.3104; 1.0545 | 0.2447; 0.1807 | 6.589; 0 | 5.301/5.303 (5.30) |
| russell-ch4-9.94Gg3 | loop-ee 1718.1 d n5 high (RESTRUCTURED 2→1) | g(4.7037,2053.31,Ls) | 0.9015 | 0.3494 | 0 | 9.943 (9.94) |
| russell-ch4-3.66gfF3 | loop-ee-1 878.9 d n2 high; loop-ee-2 365.3 d n1 resonant | g(2.4062,866.21,Ls); f(1:1,82.955,87.388) | 0.9790; 1.0000 | 0.1249; 0.0056 | 0; 7.038 | 3.665/3.660 (3.66) |
| russell-ch4-5.30ggF3 | loop-ee-1 535.0 d n1 low; loop-ee-2 709.2 d n0 single | g(1.4646,527.25,Ll); g(1.9416,338.97,U) | 1.0545; 1.5828 | 0.1807; 0.3711 | 0; 0 | 5.303/5.301 (5.30) |
| russell-ch4-5.75ggF3 | loop-ee-1 693.5 d n0 single; loop-ee-2 915.8 d n2 low | g(1.8987,323.54,U); g(2.5074,902.67,L) | 1.5790; 1.0339 | 0.3753; 0.1926 | 0; 0 | 5.751/5.750 (5.75) |
| russell-ch4-6.44Gg3 | loop-ee 762.3 d n3 low (RESTRUCTURED 2→1) | g(2.087,1111.33,L) | 0.7574 | 0.3462 | 0 | 6.442 (6.44) |
| russell-ocampo-2.5.1+0 | loop-ee-1/2/3 (366/182/365 d) | f(1:1,74.919,∓144.069); h(0.5,0,U,±15.081); f(1:1,74.919,±35.931) | 1.0; 1.0; 1.0 | 0.2107; 0; 0.2107 | 8.99; 15.081; 8.99 | 7.817 (7.8) |
| russell-ocampo-4.3.1-5 | loop-ee 548 d n1 resonant (single continuous conic) | f(1:1,84.039,∓90) + h(0.5,0,U,±5.961) | 1.0000 | 0.0000 | 5.961 | 3.097 (3.1) |

Branch vocabulary: Lambert-solvable g arcs use the project's own solver labels
(`single`/`low`/`high`, matched by running `core.lambert` on the same geometry); the
non-Lambert resonant/backflip arcs use `"resonant"`/`"half-rev"` (segments schema is
permissive; each note defines the label).

### Structural corrections (sourced from the printed leg columns)

- **Over-materialized** (template's "3-synodic ⇒ 2 loops" heuristic was wrong):
  `3.78Gg3`, `9.94Gg3`, `6.44Gg3` each print only TWO legs (one designated + ONE long
  multi-rev loop). Collapsed to a single `loop-ee`; empty `loop-ee-2` slots removed;
  `sequence_canonical` corrected E-E-E-M-M → E-E-M-M; maneuver boundary ids updated.
- **Under-materialized**: `8.049gGf2` prints THREE legs (gGf) ⇒ TWO loops. `loop-ee` →
  `loop-ee-1` (f) + `loop-ee-2` (g); `sequence_canonical` E-E-M-M → E-E-E-M-M; an interior
  Earth-flyby maneuver added.
- **Id misnomer**: `russell-ch4-3.64gGg3`'s printed third leg is f(1:1,82.995,−180.000),
  so Russell's own shorthand rule gives **3.64gGf**, not gGg. Id left unchanged (stable
  identifier; renames ripple into cyclers.space etc.) — noted in the row, follow-up task
  registered (see below).
- **mcconaghy-2006-em-k2**: loop closed with the Russell circular-coplanar realization
  (its own `orbit_source` and `free_return_arcs` ARE Russell Table 4.9 row 1), with an
  explicit model caveat: the emergent 4.995 km/s is Russell's 4.99 anchor, NOT this row's
  McConaghy-2006 ephemeris-flavored 4.7; the gap is narrowed (kind: unknown), not deleted.

### data_gaps bookkeeping

- The 12 `russell-ch4-*` + `mcconaghy` loop-specific `#54` gap entries: removed (closed)
  or narrowed (mcconaghy model caveat). The `out-em`/`ret-me` `#54-backfill` gaps remain
  (different problem: they are halves of the designated arc split by the Mars encounter).
- The two ocampo rows' `segments[*].a_au` / `segments[*].n_revs` gaps: narrowed to
  ret-me-only, with the propagation-emergent whole-g-arc conics recorded in the gap notes
  as a head start for the Mars-phase split treatment.

## Code fix: `search/descriptor.py` M:N convention was reversed (real bug)

Both primary sources agree: in f(M:N, ...), **M = Earth years (= ToF), N = spacecraft
revs** (thesis p.126: "the number following the colon ... represent[s] the number of
revolutions by the spacecraft"; paper: "M also equals the transfer time of flight in
years... N, is the number of spacecraft revolutions", a = a_E(M/N)^{2/3}).
`descriptor.py::arc_to_leg_topology` read revs from M and `arc_tof_seed_days` read years
from N — both reversed. Harmless for 1:1 arcs; wrong for 3:2 / 2:1. Numerically decisive
cross-check: only M=years makes each row's printed legs sum to n × 2.1354-yr synodic
(e.g. 5.75ggF3: 1.8987 + 2.5074 + 2 = 6.406 = 3S with F(2:1)=2 yr; the reversed reading
gives 5.41). Fixed + tests updated (`test_descriptor_arctype.py`, `test_descriptor_tof.py`).

**Blast-radius check (per the bugfix-invalidates-past-searches rule):** the parser feeds
`seed_ladder` Rung 1 (seed only, refined by correctors) — rows affected are exactly the
non-1:1-resonance descriptor rows (5.30gGf3, 3.66gfF3, 5.30ggF3, 5.75ggF3). Any past
negative search verdict that consumed Rung-1 descriptor seeds for those four rows was
seeded with a wrong ToF/rev-count for the f/F leg; registered as follow-up `#813` rather
than re-run here (seeds feed correctors, so the practical impact needs checking, not
assuming, in either direction).

## Also observed (registered, not executed here)

- `#814`: `russell-ch4-3.77Gh3` (and other h-loop ch4 rows like 8.165Gfh-f2, plus the
  Table 4.10/4.11 siblings not yet carrying `free_return_arcs[]`) are closable by exactly
  this machinery; several rows' `free_return_arcs` were simply never ingested.
  Also covers the 3.64gGg3 → 3.64gGf naming-slip disposition.
- `#815`: the ocampo ret-me/out-em designated-arc split: Table 3.5/3.7 printed vectors
  give the whole-arc conics at print precision (recorded in the gap notes); a per-SEGMENT
  writeback needs a defined Mars-phase split convention. This is the honest scope
  boundary of the dispatch's "separate, optional" `ret-me` item — NOT built this pass.
- Sign-pairing nuance: the thesis prints 4.3.1−5's h-arc as ±5.961 paired with f's ∓90.0
  (the combination that makes the zero-turn continuation exact, confirmed by propagation);
  the JSR paper's Table 2 appears to print ∓5.961 (same-sign pairing), which would demand
  a third, unlisted ~12° flyby. The thesis/catalogue version is the physically consistent
  one; the paper also notes "all permutations of the plus and minus are valid" for the
  hemisphere symmetry. Recorded here for the errata trail; no public claim warranted.

## Verification

- `uv run pytest tests/data tests/search -q` — full ratchet run 2026-08-09 ~11:12–12:15 AET
  under heavy load (self-hosted CI job + two other agents' concurrent suites): ONE failure,
  `tests/search/test_crnbp_real_ephemeris_consistency.py::test_jupiter_spice_moon_state_fn_europa_sma_sane`
  — a Jupiter/Europa SPICE-state sanity test with zero overlap with this task's changed
  files; re-run in isolation (`-n0`, single test AND whole file) passes cleanly, i.e. the
  documented contention-flake class `#809` just mitigated, not a regression. All XFAIL/XPASS
  entries in the run are the pre-existing documented cross-platform set (#584/#631/... class).
- `uv run pytest tests/search/test_descriptor_arctype.py test_descriptor_tof.py
  test_descriptor_parse.py test_descriptor_catalogue.py tests/data/test_arc_descriptors.py -q`
  — 44 passed.
- `uv run ruff check .` — clean; `ruff format --check .` — 1148 files already formatted.
- `uv run mypy src tests` — Success, no issues in 842 source files.
- YAML parse + 14-row loop-field completeness sweep — green (23/23 segments filled;
  maneuver boundary ids validated against segment ids on all restructured rows).
