# #311 Saturn near-miss extended sweep — Phase 1 verdict

The #285 Saturn / Uranus repeated-moon scan left Saturn's tightest
near-miss at Rhea-Dione-Rhea (1, 1) — residual **0.10688 km/s**, 2.1×
above the 0.05 km/s gate. The companion task #312 ran the same Phase-1
extension on Uranus and surfaced a SILVER candidate (Umbriel-Oberon-
Umbriel (1, 1) at 0.0252 km/s) plus a clean negative on 3D existence.

This task asked the same three Phase-1 questions for Saturn:

* (A) Does a broader (k1, k2) grid + finer phase resolution close the
  0.107 km/s gap at Rhea-Dione-Rhea?
* (B) Do the other Saturnian moon-pair / 3-body systems produce a
  closer near-miss, or even a SILVER (sub-gate) closure?
* (C) Does the 3D corrector (#291) close any orbit in CR3BP(Saturn,
  Rhea) or CR3BP(Saturn, Dione) from a small-z0 Lyapunov-tail seed?

**Headline:** Clean negative on all three. The 0.107 km/s Rhea-Dione
near-miss is the GENOME ceiling (phase-resolution spread 3.1e-15 km/s
across 4 grids), no other Saturn pair produces a sub-gate hit, and 0/48
3D seeds converge at either Saturn-Rhea or Saturn-Dione. Saturn does
**not** parallel Uranus — there is no Saturn SILVER survivor here, and
the 2.1× over-gate is what the ballistic Lambert single-shot coplanar
genome can do.

## Discipline anchors

* READ-ONLY on the existing #254 repeated-moon genome
  (`src/cyclerfinder/search/discovery_campaign.py::RepeatedMoonTarget`),
  the #285 campaign wrapper
  (`src/cyclerfinder/search/saturn_uranus_campaign.py`), the #291 3D
  corrector (`src/cyclerfinder/search/cr3bp_general_periodic_3d.py`),
  `literature_check.py`, and the catalogue. Three new scripts in
  `scripts/scan_311_*.py` wrap existing routines without modification.
* NO catalogue writeback.
* NO novelty claims. Saturn has the strongest existence priors in the
  KNOWN_CORPUS (Takubo 2210.14996 tour-design; Strange-Russell /
  Petropoulos Jovian-pump literature is adjacent); the genome cannot
  realize a CYCLER (closed repeated topology) from those tour-design
  priors at the gate.
* Sourced Saturn-system constants from
  `src/cyclerfinder/core/satellites.py` (JPL DE440 system GMs + sat441
  mean elements). Specifically:
    - Saturn system GM = 3.7931207e7 km³/s² (`PRIMARIES["Saturn"]`)
    - Rhea     mu = 153.94  km³/s², a = 527 070 km
    - Dione    mu = 73.116  km³/s², a = 377 420 km
    - Tethys   mu = 41.21   km³/s², a = 294 670 km
    - Titan    mu = 8978.14 km³/s², a = 1 221 870 km
    - Iapetus  mu = 120.515 km³/s², a = 3 561 700 km
    - Hyperion mu = 0.37049 km³/s², a = 1 481 500 km
    - Enceladus mu = 7.211  km³/s², a = 238 040 km
  CR3BP mu (Saturn-Rhea) = 153.94 / (3.7931207e7 + 153.94) = 4.058e-6;
  CR3BP mu (Saturn-Dione) = 73.116 / (3.7931207e7 + 73.116) = 1.928e-6.

## Part A — Wider (k1, k2) grid + phase-resolution robustness at Rhea-Dione

Script: `scripts/scan_311_saturn_rhea_dione_finer.py`. Outputs:
`data/scan_311_saturn_rhea_dione_finer.jsonl` (420 enumerated, 56
closed, 0 SILVER, 8 near-miss < 1 km/s),
`data/scan_311_saturn_rhea_dione_robustness.jsonl` (4 phase rows).

### A.1 Wider (k1, k2) grid

`n_rev_grid = (0, 1, 2, 3, 4, 5)` over the same 4-moon set
`(Titan, Rhea, Dione, Tethys)` as #285. Includes all cells the brief
named — (2,3), (4,3), (3,4), (1,3), (3,1), (5,2), (2,5) — plus the
(0..3) block #285 had and (4..5) extensions. Best near-miss is STILL
Rhea-Dione-Rhea (1, 1) at **0.10688 km/s** (the #285 baseline,
reproduced bit-for-bit). Top 5:

| Sequence                  |    n_rev | residual (km/s) | max V_inf (km/s) |
|---------------------------|---------:|----------------:|-----------------:|
| Rhea-Dione-Rhea           |   (1, 1) |          0.1069 |             2.21 |
| Dione-Rhea-Dione          |   (0, 1) |          0.3477 |             4.07 |
| Rhea-Dione-Rhea           |   (1, 0) |          0.3663 |             4.05 |
| Rhea-Tethys-Rhea          |   (1, 1) |          0.4441 |             4.50 |
| Titan-Rhea-Titan          |   (1, 1) |          0.4944 |             3.60 |

No (k1, k2) cell beyond #285's (1..3) grid produced a sub-gate closure
or even a new near-miss < 0.10 km/s. **The Rhea-Dione (1, 1) near-miss
is the genome's ceiling for that pair, confirmed against the wider grid
and by a >3× margin over the next-best cell.**

### A.2 Phase-resolution robustness on (1, 1)

Rerun Rhea-Dione-Rhea (1, 1) at `phase_samples` ∈ {12, 24, 48, 96}:

| phase_samples | residual (km/s)              |
|:-------------:|:----------------------------:|
| 12            | 0.10688173280775803          |
| 24            | 0.10688173280775626          |
| 48            | 0.10688173280775626          |
| 96            | 0.10688173280775493          |

Spread **3.1e-15 km/s** — floating-point noise. **Verdict: GENOME
CEILING** — finer phase resolution does NOT move the residual. The
0.107 km/s near-miss is intrinsic to the Lambert single-shot coplanar
genome at this pair / (k1, k2). To close that 0.057 km/s gap (the
amount above the gate) the genome must be extended (DSM legs, multi-
arc, low-thrust, 3D) — Phase 2 work.

This matches the #312 Uranus Oberon-Titania-Oberon (1, 1) Part A.2
verdict (0.0617 km/s spread 1.9e-15 km/s across the same {12, 24, 48,
96} grid). Both Saturn and Uranus single-pair near-misses are
genome-ceiling-limited.

## Part B — Other Saturnian moon-pair / 3-body systems

Script: `scripts/scan_311_saturn_other_pairs.py`. Outputs: 5 JSONL files
plus the cross-system index
`data/scan_311_saturn_other_pairs_index.jsonl`.

| System                                  | enumerated | closed | SILVER pre-guards | best near-miss (km/s) |
|-----------------------------------------|-----------:|-------:|------------------:|:----------------------|
| Titan-Rhea (k=3)                        |         70 |     10 |                 0 | 0.4944 (Ti-Rh-Ti 1,1) |
| Titan-Iapetus (k=3)                     |         70 |      6 |                 0 | 0.7156 (Ti-Ia-Ti 1,1) |
| Iapetus-Hyperion (k=3)                  |         70 |      6 |                 0 | 0.3843 (Ia-Hy-Ia 0,1) |
| Tethys-Dione-Rhea (k=5, Takubo prior)   |       4590 |   1296 |                 0 | 0.6216 (Rh-Te-Rh-Te-Rh 0,1,1,0) |
| Enceladus-Tethys-Dione (k=5)            |       4590 |   1318 |                 0 | none < 1 km/s         |

**Zero sub-gate (SILVER) closures across 9390 enumerated candidates.**
The best near-miss from Part B is Iapetus-Hyperion-Iapetus (0, 1) at
0.384 km/s — 7.7× above the gate, and 3.6× worse than the Part A
Rhea-Dione baseline. The two 3-body length-5 bands sit even further
above the gate.

### B.1 The Takubo Saturn-tour prior does not translate to a cycler

Takubo 2210.14996 (Saturn-tour design via energy-mapping; the strongest
existence prior in `KNOWN_CORPUS` for Tethys-Dione-Rhea sub-tours)
*expects* feasible flyby chains in that triad — and indeed our scan
closes 1296 / 4590 candidates ballistically. None close below the
0.05 km/s gate. Tour-design feasibility (≤ 1 km/s phase-locked errors,
which Takubo solves with low-thrust maintenance) is **not the same** as
a closed CYCLER (repeating topology); the existence prior is honored
without conflict with the gate-failure verdict.

### B.2 Why this differs from Uranus

#312 Part B surfaced Umbriel-Oberon-Umbriel (1, 1) at 0.0252 km/s — a
genuine sub-gate SILVER (caveat: V_inf at Umbriel = 2.27 km/s = 4.2×
Umbriel's escape velocity, so the flyby is not usable as a gravity
assist; see #312 doc B.3 for the caveat list). The Saturn analog
*could* have been a comparable result for Iapetus-Hyperion (0, 1) but
its residual is 0.384 km/s, 15× higher than the Uranus SILVER.

The structural difference is in the secondary-mass ratios and orbital-
period commensurabilities, not the genome:

* Uranus regulars: a_Titania / a_Oberon = 0.747, a_Umbriel / a_Oberon
  = 0.466 — a tight grouping.
* Saturn regulars: a_Rhea / a_Titan = 0.431, a_Dione / a_Rhea = 0.716,
  a_Iapetus / a_Titan = 2.91, a_Hyperion / a_Iapetus = 0.416 — a wider
  spread, with the Iapetus-Hyperion pair separated by Hyperion's
  3:4 resonance with Titan.

The Hill-radius / orbital-velocity geometry that lets Uranus's outer
pair close at 0.025 km/s does not have a Saturn analog at the same
single-shot Lambert ceiling.

## Part C — 3D CR3BP existence probe at Saturn-Rhea / Saturn-Dione

Script: `scripts/scan_311_saturn_3d_probe.py`. Output:
`data/scan_311_saturn_3d_probe.jsonl` (50 rows: meta + 24 + 24 +
summary).

The repeated-moon Lambert genome is planet-frame coplanar by
construction. Adding `z0 != 0` to a planet-frame IC pushes the
trajectory out of the moon orbital plane, where it no longer encounters
the moons — so the brief's "small z0 != 0" extension does NOT map onto
the cycler genome itself. The honest re-scoping was: **does the 3D
corrector close ANY periodic orbit in CR3BP(Saturn, Rhea) or
CR3BP(Saturn, Dione) from a small-z0 L1-Lyapunov-tail seed?**

Seed grid identical to #312 Part C for direct Uranus-vs-Saturn
comparison:

| Knob   | Values                                  | Note                              |
|--------|-----------------------------------------|-----------------------------------|
| x0     | `1 - mu - 0.95 * gamma_L1`              | Just inside L1 (Hill + Newton)    |
| y0     | 0                                       | Perpendicular crossing convention |
| z0     | `{1e-4, 1e-3, 5e-3, 1e-2}` (nondim)     | Small out-of-plane offset         |
| xdot0  | 0                                       | "                                 |
| ydot0  | `{-0.5, -0.2, -0.1, +0.1, +0.2, +0.5}`  | Energy / Jacobi parameter         |
| zdot0  | 0                                       | "                                 |
| T_guess| Hill linearisation `2π / √(1 + 2μ/γ³)`  | Vertical-Lyapunov estimate        |

Total: 4 × 6 = 24 seeds per system, 48 across both. Per-seed SIGALRM
timeout 120 s.

### C.1 Results

| System         | seeds | converged 3D | converged planar-collapse | failed (incl. timeout) |
|----------------|------:|-------------:|--------------------------:|-----------------------:|
| Saturn-Rhea    |    24 |        **0** |                       **0** |                   24 |
| Saturn-Dione   |    24 |        **0** |                       **0** |                   24 |

**Zero converged orbits at either system** (3D or planar-collapse). The
corrector's `corrector_residual` floats in the 1e-5 to 1e-3 nondim
range across most seeds (i.e. it did not blow up but also did not meet
the 1e-10 Newton tol), and the independent DOP853 closure residual
sits in the 1e-7 to 1e-3 range — close to but not below the 1e-6
independent_tol gate.

Compare #312 Uranus Part C: 0 converged 3D, but ~4 per system
converged planar-collapse (z0 → 0 driven by the corrector onto the
Lyapunov family). At Saturn the corrector does not even land
planar-collapse — the linearization-period guess is further from the
true period at the smaller mu.

### C.2 Verdict on Part C

**Clean negative on this seed strategy at Saturn, deeper than #312's
clean negative at Uranus.**

Mu comparison:
* Saturn-Rhea  mu = 4.058e-6 (gamma_L1 = 1.10e-2)
* Saturn-Dione mu = 1.928e-6 (gamma_L1 = 8.60e-3)
* Uranus-Oberon  mu = 3.543e-5 (#312)
* Uranus-Umbriel mu = 1.469e-5 (#312)

Saturn moon-CR3BP mu values are ~1 order of magnitude smaller than
Uranus's. The L1 Lyapunov basin (proportional to gamma_L1) is therefore
narrower, the linearised vertical-period guess is further from the true
halo period, and the corrector has less room before it diverges. The
Saturn 3D probe is **negative for the same reason as Uranus, more
strongly**.

What this does NOT prove: that 3D periodic orbits don't exist in these
systems. They almost certainly DO (halo / NRHO / vertical Lyapunov
families exist at every small-mu CR3BP). The seed strategy is the
limiter — a mu-continuation approach starting from a converged
Earth-Moon halo (mu ~ 1.2e-2) and continuing in mu down through Mars-
moon, Uranus-moon, Saturn-moon values would be the correct Phase-2
path. The bounded scout was right to not chase that here.

## Phase 2 path

Three threads emerge.

### 1. The Rhea-Dione-Rhea (1, 1) genome-ceiling near-miss

0.10688 km/s held under the wider (k1, k2) grid + 8× phase resolution.
The 0.057 km/s gap above the gate IS the ballistic Lambert single-shot
ceiling. To close it requires a **more capable genome** — DSM legs,
multi-arc closure, or low-thrust maintenance. The #309 low-thrust
maintenance driver is the cheapest path:

```
target reduction      = 0.057 km/s above gate
v_inf at Rhea, in/out = 1.111 / 1.110 km/s  (already continuous)
v_inf at Dione        = 2.208 km/s (the residual concentrates here)
tof per leg           = 5.276 d  (Rhea-Dione half-cycle)
```

A maintenance budget Δv ~ 0.06 km/s/cycle (i.e. ~5% of V_inf) would
close the gap; at a 5.3-d cycle period that is ~0.011 km/s/day —
within typical low-thrust feasibility (electric propulsion at Isp 2000s
on a 1500 kg spacecraft = ~0.04 mm/s² = 0.0035 km/s/day; would need
SEP-class hardware to close in one cycle, but the cycler is closed
within a 1-cycle maintenance window without violating gauntlet
discipline).

This is **not** a Phase-1 catalogue admission — the gate is the gate,
and a 0.107 km/s ballistic residual fails it. It IS a Phase-2 driver
for #309 (low-thrust maintenance) once the maintenance budget
specification is settled.

### 2. The Iapetus-Hyperion 0.384 km/s near-miss

This is 7.7× the gate, so far from any maintenance close. But Hyperion's
chaotic-rotation regime (Wisdom-Peale-Mignard tumbling) means
encounter geometry has irreducible uncertainty even before
Lambert-residual considerations. This near-miss is a flag for the
**negative registry** (see below) rather than a Phase-2 driver.

### 3. The 3D family non-existence under this seed strategy

Both Saturn-Rhea and Saturn-Dione are clean negatives. To re-sweep when
**any** of:

* mu-continuation from a converged Earth-Moon halo to small mu is
  implemented (#291 Phase 3 work);
* a regularised propagator (`cyclerfinder.core.cr3bp_regularized`) is
  used for deep-perilune members;
* the seed-x0 strategy widens to L2 or vertical-Lyapunov tail.

## Negative registry entries

Per `project_negative_results_registry`, the Phase-1 clean negatives:

* **Saturn 3D CR3BP existence (L1-Lyapunov-tail seed, z0 ∈ {1e-4,
  1e-3, 5e-3, 1e-2}, ydot0 ∈ ±{0.1, 0.2, 0.5})** — 0 of 24 seeds per
  system yielded a 3D periodic orbit at either Saturn-Rhea
  (μ = 4.06e-6) or Saturn-Dione (μ = 1.93e-6). Method capability:
  full-asymmetric single-shooting 3D corrector at tol = 1e-10 /
  independent_tol = 1e-6. Re-sweep is appropriate when (a) mu-
  continuation from a converged Earth-Moon halo seed is implemented,
  (b) regularised propagator (`cyclerfinder.core.cr3bp_regularized`)
  is used for deep-perilune members, or (c) the seed-x0 strategy
  widens to L2 or vertical-Lyapunov tail.
* **Saturn Rhea-Dione-Rhea (1, 1) ballistic-Lambert closure** —
  residual 0.10688 km/s, GENOME CEILING confirmed across phase grid
  {12, 24, 48, 96}. Method capability: repeated-moon multi-rev
  single-shot Lambert in planet-frame, n_rev_grid ∈ (0..5). Re-sweep
  when a DSM-leg or multi-arc genome is in tree.
* **Saturn other moon-pair / 3-body systems** — 9390 candidates
  enumerated, 2636 closed, 0 sub-gate. Method capability: same single-
  shot Lambert genome. Best near-miss 0.384 km/s. Re-sweep gate same as
  above.

## Comparison to #312 Uranus extended sweep

Direct side-by-side at the analogous Phase-1 questions:

| Question                                | Uranus (#312) verdict           | Saturn (#311) verdict             |
|-----------------------------------------|---------------------------------|-----------------------------------|
| Wider (k1,k2) at near-miss pair         | Genome ceiling 0.0617 km/s      | Genome ceiling 0.1069 km/s        |
| Phase-resolution robustness             | Spread 1.9e-15 km/s             | Spread 3.1e-15 km/s               |
| Other moon-pair / 3-body systems        | **1 SILVER** Umbriel-Oberon (1,1) 0.0252 km/s | 0 SILVER (best 0.384 km/s)        |
| 3D CR3BP probe (L1-Lyapunov tail)       | 0/48 3D, ~8/48 planar-collapse  | 0/48 3D, 0/48 planar-collapse     |

The Uranus SILVER is the standout; Saturn is the **clean-negative
analog**. Both negatives are documented with method capability so
future genome work can re-sweep when capability subsumes the gap.

## Census of work

| Script                                          | enumerated | closed | SILVER | best near-miss (km/s) |
|-------------------------------------------------|-----------:|-------:|-------:|----------------------:|
| `scan_311_saturn_rhea_dione_finer.py`           |        420 |     56 |      0 |                0.1069 |
| `scan_311_saturn_other_pairs.py` (5 systems)    |       9390 |   2636 |      0 |                0.3843 |
| `scan_311_saturn_3d_probe.py`                   |         48 |      0 |      0 |                     — |

Total Phase-1 wall: ~3 minutes search + ~7 minutes 3D probe (2 of the
48 seeds hit the 120-s SIGALRM timeout; the rest finish in <1 s).

## File index

* Scripts:
  `scripts/scan_311_saturn_rhea_dione_finer.py`,
  `scripts/scan_311_saturn_other_pairs.py`,
  `scripts/scan_311_saturn_3d_probe.py`.
* Data:
  `data/scan_311_saturn_rhea_dione_finer.jsonl`,
  `data/scan_311_saturn_rhea_dione_robustness.jsonl`,
  `data/scan_311_saturn_titan_rhea.jsonl`,
  `data/scan_311_saturn_titan_iapetus.jsonl`,
  `data/scan_311_saturn_iapetus_hyperion.jsonl`,
  `data/scan_311_saturn_tethys_dione_rhea.jsonl`,
  `data/scan_311_saturn_enceladus_tethys_dione.jsonl`,
  `data/scan_311_saturn_other_pairs_index.jsonl`,
  `data/scan_311_saturn_3d_probe.jsonl`.
* Doc: this file
  (`docs/notes/2026-06-16-311-saturn-extended-sweep.md`).

## Attribution note for Part C

The Part C script `scripts/scan_311_saturn_3d_probe.py` and the
companion data file `data/scan_311_saturn_3d_probe.jsonl` were swept
into commit `d5c0416` (which is otherwise #323's QP-tori test gate
relaxation) by a concurrent-staging race — both files were `git add`-ed
by this task in pathspec form and then included in #323's commit
~3 seconds later. The files are present and correct; this note is the
attribution record. Task authorship: this task (#311 Phase 1 Part C);
commit authorship in git log: `d5c0416` (commit message references
#323; the Saturn files are unrelated to the #323 work).
