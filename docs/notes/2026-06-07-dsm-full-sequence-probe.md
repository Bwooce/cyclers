# Full multi-arc sequence probe for 6.44Gg3 — task #153 (follow-up to #150)

**Date:** 2026-06-07
**Code:** `tests/search/test_dsm_leg.py` (additive; `src/cyclerfinder/search/dsm_leg.py`
NOT edited — the probe is pure test-side composition over the existing surface).
**Context read (in order):** `docs/notes/2026-06-07-dsm-leg-genome.md` (the #150
clean E→M→E negative + its follow-up: drive the FULL multi-leg sequence),
`src/cyclerfinder/search/dsm_leg.py`, `src/cyclerfinder/search/mbh.py` (read-only),
`data/catalogue.yaml` `russell-ch4-6.44Gg3`, `docs/notes/2026-06-07-mbh-wrapper.md`
(the sourced anchors + failed single-ellipse / minimal-chain attempts),
`docs/notes/multi-arc-classification.md` (descriptor semantics: first descriptor
number = TOF in YEARS).

## What this probe is

#150 ran the **minimal** `E→M→E` two-leg chain (the control) and floored at
9.40 km/s, predicting the missing piece was the **full multi-arc topology** — the
two generic Earth-Earth free-return arcs the row's descriptor actually encodes,
each as its own DSM leg. This task drives that full sequence.

## The sequence the descriptor encodes (derivation)

`russell-ch4-6.44Gg3` `free_return_arcs[]`:

| arc | raw_descriptor | tof_years |
|---|---|---|
| g | `g(2.087,1111.33,L)` | 2.087 |
| G | `G(4.3191,1194.88,L)` | 4.3191 |

The two arc ToFs **sum to 2.087 + 4.3191 = 6.406 yr = the row's `period.years`
(6.41)**, so the two generic Earth-Earth free returns ARE the whole cycle (the
262-day E→M and M→E transits are the Mars-bracketing *pieces* interior to the
arcs, not additional legs). Unrolling each generic E-E arc into its
Mars-bracketing pieces gives the period-matching sequence

```
E -M-> E -M-> E      (sequence E-M-E-M-E, 4 legs)
  \____g____/  \____G____/
   2.087 yr      4.3191 yr
```

| leg | bodies | ToF (d) | source of the ToF |
|---|---|---|---|
| 1 | E→M | 262.0 | `invariants.transit_times_days[0]` / `segments[out-em]` |
| 2 | M→E | g_arc − 262 = 500.4 | g arc 2.087 yr (`free_return_arcs[0]`) minus the outbound transit |
| 3 | E→M | 262.0 | `invariants.transit_times_days[0]` |
| 4 | M→E | G_arc − 262 = 1315.5 | G arc 4.3191 yr (`free_return_arcs[1]`) minus the outbound transit |

(The alternate reading — two loop arcs as *native* E→E legs plus the two transits,
sequence `E-E-E-M-E` — was also run, Config 3 below; it sums to 7.84 yr, NOT the
sourced 6.41, and floors even harder. The period-matching `E-M-E-M-E` is the
defensible reading.)

## Seeding (every quantity sourced)

| seeded quantity | value | source |
|---|---|---|
| departure V∞ magnitude | 6.44 km/s | `vinf_kms_at_encounters[0]` (Earth anchor; direction-free) |
| α₀ (departure azimuth) | +π/2 | tangential prograde at t=0 on the circular backend (refined by MBH; only sets the half-plane) |
| β₀ (elevation) | 0 | coplanar seed |
| leg ToFs | (262, 500.4, 262, 1315.5) d | descriptor arc decomposition above |
| η per leg | 0.5 | Takao default DSM fraction (free to move) |
| t0 | 0 | epoch-free circular family |

SOURCED anchors (the **EXPECTED** side, never imposed): V∞ E = 6.44, M = 3.74 km/s.
The emerged V∞ is **EVIDENCE**. Constraint-vs-evidence separation is absolute —
nothing our code computes is ever an EXPECTED value.

## Configurations run (MBH: cauchy, rng_seed=6, ≤120 hops, stall 60)

All driven through the existing `make_dsm_chain_step` MBH adapter. Per-gene
absolute hop scales: t0 ±5 d, V∞ ±0.5 km/s, α ±0.2, β ±0.1, ToF ±20 d, η ±0.1.

### Config 1 — E-M-E-M-E, ToF bounded ±30%, η free — verbatim

```
feasible=False   total_dV = 29.94 km/s
hops attempted/accepted = 61/0   (stopped on stall)
emerged V_inf_in  {1(E→M):12.24, 2(M→E):15.89, 3(E→M):6.69, 4(M→E):7.55}   sourced E=6.44, M=3.74
per-leg dV_DSM = (8.66, 8.13, 11.60, 1.54) km/s
eta per leg    = (0.448, 0.491, 0.545, 0.502)
tof days/leg   = (262.0, 500.4, 262.0, 1315.5)
wall = 22.7 s
```

### Config 2 — E-M-E-M-E, ToF FROZEN, only V∞/α/β/η move — verbatim

```
feasible=False   total_dV = 29.92 km/s
hops attempted/accepted = 61/0
emerged V_inf_in  {1:12.22, 2:15.88, 3:6.72, 4:7.57}
per-leg dV_DSM = (8.64, 8.15, 11.60, 1.52) km/s
eta per leg    = (0.447, 0.491, 0.546, 0.502)
wall = 29.2 s
```

Config 2 floors **identically** to Config 1 — the floor is NOT a ToF-search
failure; it is structural.

### Config 3 — E-E-E-M-E (loop arcs as native E→E legs), ToF ±30% — verbatim

```
feasible=False   total_dV = 49.18 km/s
hops attempted/accepted = 61/0
emerged V_inf_in  {1:13.92, 2:6.39, 3:11.27, 4:7.35}
per-leg dV_DSM = (7.53, 25.87, 8.23, 7.55) km/s
wall = 19.5 s
```

(For reference, the #150 minimal E→M→E control floors at 9.40 km/s.)

## Verdict — a NEGATIVE with a PRECISE mechanism (not a basin-selection miss)

Every full-sequence configuration floors at **30–49 km/s** with **0 of 61 hops
accepted** and emerged encounter V∞ (6–16 km/s) far from the sourced 6.44 / 3.74.
Crucially, this is **worse** than the #150 minimal chain (9.40 km/s), and the
mechanism is now identified exactly:

**`dsm_leg`'s back-arc Lambert is SINGLE-revolution (`max_revs=0`), but the
generic loop arcs are MULTI-revolution.** A 1.27–1.54 AU heliocentric ellipse has
a period of ~1.4–1.9 yr, so the g arc (2.087 yr) is ~1.5 revolutions and the G arc
(4.3191 yr) is ~3 revolutions. Single-rev Lambert over a >1-period ToF is forced
onto the degenerate near-radial high-energy branch. Direct comparison on the
G-arc M→E return (1315.5 d), circular backend — verbatim:

```
single-rev (max_revs=0):           |v1| = 28.18 km/s   <- the branch dsm_leg is locked to
multi-rev  (max_revs=3) best:       |v1| = 15.51 km/s
```

The single-rev branch is ~13 km/s worse than the best multi-rev branch. The DSM
impulse on the long legs (legs 2 and 4) is dominated by this Lambert degeneracy,
not by any real maneuver the cycler requires.

### Correction to the #150 follow-up framing

#150 concluded the missing piece "is not the interior-impulse mechanic but the
**multi-leg sequence** … a sequence/seed change, not new mechanics." This probe
shows that is **only half right**. The evaluator does string arbitrary legs
correctly, and the sequence/seed are now sourced honestly — but each leg's
back-arc Lambert is single-rev, so the *defining feature* of this multi-arc
cycler (its multi-revolution loop arcs) **cannot be represented by the current
primitive**. The real frontier blocker is mechanical after all: **`dsm_leg` needs
a multi-rev Lambert branch (`max_revs > 0`) and a branch-selection coordinate**
before the full topology is representable. The single-rev primitive is correct for
sub-period transfer legs (Gates 1–3 prove that) and is the right tool for the E→M
/ M→E transits; it is the wrong tool for the resonant loop arcs.

This is the **honest bound**: the full multi-arc sequence floors at 30 km/s under
the single-rev primitive, with the floor traced to a specific, fixable mechanical
limitation rather than to basin selection. No close-and-match; **no catalogue
writeback** (the row stays multi-arc with its sourced anchors).

## Tests added (behind `@pytest.mark.slow`, qualitative assertions only)

`tests/search/test_dsm_leg.py` — the existing 7 stay green; 2 added (9 total):

- `test_dsm_644gg3_full_sequence_probe` — drives the E-M-E-M-E full sequence via
  MBH; asserts only that the search RUNS and returns a finite, fully-audited
  result (the scientific verdict is reported here, not gated). Prints the
  emerged-vs-sourced evidence.
- `test_dsm_leg_single_rev_floors_on_multirev_loop_arc` — CONSTRUCTED: the G-arc
  M→E single-rev Lambert |v1| is strictly worse (>1 km/s) than the best multi-rev
  branch, and more multi-rev branches exist — the precise mechanism behind the
  negative. The reference is the multi-rev Lambert family, NOT a DSM-evaluator
  output (golden-rule separation).

Verbatim: `uv run pytest tests/search/test_dsm_leg.py -m 'slow or not slow'` →
**9 passed in ~39 s**. ruff check + ruff format --check + mypy on both files: clean.

## Follow-up (the now-precise next step)

**Add a multi-rev Lambert branch to `dsm_leg`** (`max_revs` parameter + a
branch-selection genome coordinate, e.g. an integer revs count per leg or a
discrete branch index). Re-run this exact full-sequence probe with the loop-arc
legs free to take 1–3 revolutions; that is the first configuration in which the
6.44Gg3 multi-arc topology is even *representable*. Until then the single-ellipse
(#137) and single-rev one-DSM (#150/#153) genomes both provably lack the
sourced-anchor basin for this row — a consistent, twice-confirmed bound.

The s1l1-2syn-em-cpom two-arc probe (5.65 / 3.05; arcs g 1.4612 yr + G 2.8096 yr)
was left for the optional slot but is the **same** structural situation: those arcs
(1.46 / 2.81 yr) are also multi-revolution, so it would floor for the same reason
and is better deferred until the multi-rev primitive exists.
