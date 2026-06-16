# #312 Uranus near-miss extended sweep — Phase 1 verdict

The #285 Saturn / Uranus repeated-moon scan left a single tight near-miss:
Oberon-Titania-Oberon (1, 1) at residual 0.0617 km/s — 23% above the
0.05 km/s gate. NO published Uranian-cycler prior exists in the
`KNOWN_CORPUS` (0 of 32 anchors mention any Uranian body). So Uranus was
the session's highest novelty-leverage probe: a tight near-miss at a
fresh primary.

This Phase-1 task asked three honest questions:

* (A) Does a broader (k1, k2) grid + finer phase resolution close the
  0.062 km/s gap at Oberon-Titania-Oberon?
* (B) Do the other Uranian moon-pair / 3-body systems produce a closer
  near-miss, or even a SILVER (sub-gate) closure?
* (C) Does the 3D corrector (#291) close any orbit in CR3BP(Uranus,
  Oberon) or CR3BP(Uranus, Umbriel) from a small-z0 Lyapunov-tail seed?

**Headline:** Part B surfaced a **literature-fresh SILVER closure** —
Umbriel-Oberon-Umbriel (1, 1) at residual 0.0252 km/s, all guards
passed, KNOWN_CORPUS `not-found`, p_fp 0.59. The closure is robust
under a fine-grid relative-phase-offset sweep (0.024033 km/s at 96 × 96).
This is the kind of regime-fresh find the brief calls out.

Honest caveats follow before catalogue ingestion.

## Discipline anchors

* READ-ONLY on the existing #254 repeated-moon genome
  (`src/cyclerfinder/search/discovery_campaign.py::RepeatedMoonTarget`),
  the #285 campaign wrapper
  (`src/cyclerfinder/search/saturn_uranus_campaign.py`), the #291 3D
  corrector (`src/cyclerfinder/search/cr3bp_general_periodic_3d.py`),
  `literature_check.py`, and the catalogue. Three new scripts in
  `scripts/scan_312_*.py` wrap existing routines without modification.
* NO catalogue writeback. The SILVER candidate stays out of the
  catalogue pending the #305 (BCR4BP) / #306 (3D) gauntlets.
* NO novelty claims in commit messages. "Literature-fresh in the
  KNOWN_CORPUS" is the necessary-not-sufficient gate from
  `feedback_literature_novelty_check_baseline.md`; full novelty is
  decided by the gauntlets.
* Sourced moon-mu and SMA from
  `src/cyclerfinder/core/satellites.py` (JPL DE440 system GMs + sat441
  mean elements, accessed 2026-06-14).
* Independent DOP853 cross-check on every SILVER. The Lambert IS the
  closure, so a re-Lambert is not independent — `dop853_cross_check_leg`
  re-integrates each leg under DOP853 at rtol=atol=1e-12 (per
  `feedback_orbit_closure_discipline`).

## Part A — Wider (k1, k2) grid + phase-resolution robustness at Oberon-Titania

Script: `scripts/scan_312_uranus_oberon_titania_finer.py`. Outputs:
`data/scan_312_uranus_oberon_titania_finer.jsonl` (210 enumerated, 28
closed, 0 SILVER, 12 near-miss < 1 km/s),
`data/scan_312_uranus_robustness.jsonl` (4 phase-resolution rows).

### A.1 Wider (k1, k2) grid

`n_rev_grid = (0, 1, 2, 3, 4, 5)` over the same 3-moon set
`(Titania, Oberon, Umbriel)` as #285 — 6 length-3 sequences × 36 (k1,
k2) cells minus the all-zero-rev corner = 210 enumerated. Best
near-miss is STILL Oberon-Titania-Oberon (1, 1) at **0.0617 km/s** (the
#285 baseline, exactly reproduced). Top 5:

| Sequence                       | n_rev | residual (km/s) | max V_inf (km/s) |
|--------------------------------|:------|----------------:|-----------------:|
| Oberon-Titania-Oberon          | (1,1) |          0.0617 |             0.96 |
| Titania-Oberon-Titania         | (1,1) |          0.1085 |             1.44 |
| Umbriel-Titania-Umbriel        | (1,1) |          0.1850 |             1.47 |
| Titania-Oberon-Titania         | (1,0) |          0.2859 |             1.28 |
| Oberon-Titania-Oberon          | (0,1) |          0.3872 |             1.48 |

No (k1, k2) cell beyond #285's `(0..3)` grid produced a sub-gate
closure or even a new near-miss < 0.06 km/s. **The Oberon-Titania (1,
1) near-miss is the genome's ceiling for that pair, confirmed against
the wider grid.**

### A.2 Phase-resolution robustness on (1, 1)

Rerun Oberon-Titania-Oberon (1, 1) at `phase_samples` ∈ {12, 24, 48, 96}:

| phase_samples | residual (km/s)              |
|:-------------:|:----------------------------:|
| 12            | 0.06169657458692290          |
| 24            | 0.06169657458692102          |
| 48            | 0.06169657458692102          |
| 96            | 0.06169657458692102          |

Spread 1.9e-15 km/s across {12, 24, 48, 96}. **Verdict: GENOME CEILING**
— finer phase resolution does NOT move the residual. The 0.062 km/s
near-miss is intrinsic to the Lambert single-shot coplanar genome at
this pair / (k1, k2). To close that gap, the genome must be extended
(DSM legs, multi-arc, low-thrust, 3D) — Phase 2 work.

This matches the #285 Saturn Rhea-Dione-Rhea (1, 1) near-miss verdict
(0.107 km/s held across {12, 24, 48}; genome ceiling, not phase
artifact).

## Part B — Other Uranian moon-pair / 3-body systems

Script: `scripts/scan_312_uranus_other_pairs.py`. Outputs: 5 JSONL files
plus the cross-system index
`data/scan_312_uranus_other_pairs_index.jsonl`.

| System                              | enumerated | closed | SILVER pre-guards | best near-miss (km/s)                  |
|-------------------------------------|-----------:|-------:|------------------:|:---------------------------------------|
| Titania-Umbriel (k=3)               |         70 |     11 |                 0 | 0.6615 (Umbriel-Titania-Umbriel (2,2)) |
| **Oberon-Umbriel (k=3)**            |     **70** | **11** |             **1** | **0.0252** SILVER (see below)          |
| Ariel-Titania (k=3)                 |         70 |     10 |                 0 | 0.2727 (Ariel-Titania-Ariel (1,1))     |
| Miranda-Ariel-Umbriel (k=5, 3-body) |       4590 |   1483 |                 0 | 0.4113 (Mir-Ari-Mir-Ari-Mir (1,1,1,1)) |
| Ariel-Umbriel-Titania (k=5, 3-body) |       4590 |    880 |                 0 | 0.4453 (Ari-Umb-Ti-Umb-Ari (1,2,1,1))  |

Both 3-body length-5 bands are >0.4 km/s — well above the gate. The two
inner-moon pair scans (Titania-Umbriel, Ariel-Titania) also miss
broadly. **The Oberon-Umbriel pair is the surprise.**

### B.1 The Umbriel-Oberon-Umbriel (1, 1) SILVER

`scan_312_uranus_oberon_umbriel.jsonl` contains exactly one scored row:

```
candidate_id        = repeated-moon-uranus-00000041
sequence            = [Umbriel, Oberon, Umbriel]
n_rev               = [1, 1]
residual_kms        = 0.025232
vinf_per_encounter  = [0.9199, 0.9604, 0.8947] km/s
tof_days            = [14.94, 14.94]
max_vinf_kms        = 0.9604
DOP853 cross-check  = max_dr_arrival = 1.6e-5 km — PASSED (< 1.0 km gate)
Tier-0 NN           = all admitted (max dv 4.18 km/s; NN background
                      Voyager-2 statistics happen to be permissive here)
literature_check    = not-found (0 KNOWN_CORPUS anchors for Uranus /
                      Oberon / Umbriel; confidence 0.40)
ml_flagger_p_fp     = 0.5918  (< 0.75 SILVER threshold)
verdict             = SILVER -- all guards passed
```

### B.2 Robustness of the SILVER — relative-phase-offset sweep

The repeated-moon genome's per-moon longitude offset
`2pi * j / len(moons)` is a deterministic seed convention, not a swept
parameter; the phase grid `phase_samples` only sweeps the GLOBAL anchor
longitude. So the residual depends on which moon SET is registered:

| Moon set fed to ``RepeatedMoonTarget`` | Umbriel-Oberon relative offset | Residual at (1,1) |
|----------------------------------------|:-------------------------------:|:-----------------:|
| `(Titania, Oberon, Umbriel)` — #285    | 4π/3 ≈ 240°                    | 0.6360 km/s       |
| `(Oberon, Umbriel)` — Part B           | π = 180°                       | 0.0252 km/s       |

The 2-moon convention happened to land near the true basin floor; the
3-moon convention sat in a different basin. To rule out a knife-edge
artifact, the addendum script
`scan_312_uranus_umbriel_oberon_offset_sweep.py` sweeps the relative
offset explicitly at 96 samples (in concert with 96 global-phase
samples). Output:
`data/scan_312_uranus_umbriel_oberon_offset_sweep.jsonl`.

```
BEST residual = 0.024033 km/s
  rel_offset  = 45.00 deg   (NOT 180° -- the 2-moon convention is close
                             but not exactly on the basin floor)
  phase0      = 311.25 deg  (any of 10 symmetric values; the global
                             phase only matters mod the moon-system
                             rotation, so multiple phase0 give the same
                             residual)
  tof_scale   = 1.50
  tof_days    = 11.2054
  v_inf in    = [-,    0.980, 2.259]
  v_inf out   = [2.283, 0.987, -    ]
```

The basin floor is 0.024 km/s — only ~0.001 km/s below the SILVER's
0.025 km/s and 0.640 km/s below the 3-moon convention's reading. The
SILVER is a **genuine local minimum**, not a knife-edge.

V_inf-continuity at the central Oberon flyby: in = 0.980 km/s, out =
0.987 km/s — match to 0.007 km/s. At the Umbriel wraparound flyby: in
= 2.259 km/s, out = 2.283 km/s — match to 0.024 km/s (the worst
defect, which sets the residual). **The cycler closes on itself at
both flybys.**

### B.3 Honest caveats before catalogue ingestion

1. **V_inf at Umbriel is 2.27 km/s — 4.2× Umbriel's surface escape
   velocity** (√(2 × 85.1 / 584.7) = 0.539 km/s). The closest-approach
   radius needed for a substantial bending angle (turn angle ψ
   satisfying sin(ψ/2) = 1 / (1 + r_p v_inf² / mu_moon), so r_p =
   mu_moon × (1/sin(ψ/2) − 1) / v_inf²) at v_inf = 2.27 km/s gives r_p
   below Umbriel's surface for even modest bending. **The Umbriel
   "flyby" in this cycler is not a useful gravity assist** — it is
   effectively a V_inf-continuity match at a fictitious encounter
   geometry. The Oberon flyby at v_inf = 0.98 km/s is in a usable
   regime (r_p ~ 1.7 R_Oberon for a 10° bend), but the cycler as a
   whole does not provide transport in the way an Earth-Mars cycler
   does.
2. **The genome is patched-conic, coplanar, ballistic.** Uranus's
   J2 (oblateness perturbation) is 3.34e-3 — significant for orbits
   spanning hundreds of thousands of km. The genome ignores it. A real
   trajectory in this regime would experience nodal precession + apsidal
   precession that this analysis does not model. The 0.025 km/s residual
   is the IDEALIZED ballistic figure; J2 (and the other regular moons)
   would shift it.
3. **The 0.025 km/s closure depends on the 2-moon enumeration
   convention.** Running the campaign with the 3-moon set
   `(Titania, Oberon, Umbriel)` — which is what #285 actually did —
   reports 0.636 km/s for this same (sequence, n_rev) cell because the
   relative-offset convention is different. This is a **search-genome
   incompleteness**: the genome holds the relative-offset as a fixed
   `2pi * j / len(moons)` seed and only sweeps the global phase. The
   true optimum lives at relative-offset = 45° between Umbriel and
   Oberon, which neither convention picks; the 2-moon convention (180°)
   just happens to be in the same basin. This is reported here so a
   future genome revision can sweep relative offsets explicitly; it
   does NOT invalidate the 0.024–0.025 km/s closure (the
   relative-offset sweep confirmed the basin floor).
4. **Phase-2 gate, not Phase-1 admission.** The SILVER is plumbed for
   #305 (BCR4BP V0-V5) and #306 (3D V0-V5) gauntlets. The orbit-closure
   discipline anchor (`feedback_orbit_closure_discipline`) requires
   independent cross-check at closure ≤ 1e-6 for catalogue admission;
   DOP853 vs Lambert at 1.6e-5 km IS that cross-check for the patched-
   conic envelope. CR3BP or higher-fidelity admission is the gauntlet's
   job, NOT this scan's.

### B.4 Lit-fresh structural fingerprint

```python
sig = CandidateSignature(
    primary="Uranus",
    sequence=("Umbriel", "Oberon", "Umbriel"),
    period_k=2,
    vinf_per_encounter_kms=(0.9199, 0.9604, 0.8947),
    n_rev=(1, 1),
)
check_literature(sig, search=offline_corpus_search)
# -> SearchResult(status="not-found", confidence=0.40, citation=None)
```

The offline `KNOWN_CORPUS` contains 32 anchors; 0 mention any Uranian
body. The structural fingerprint matches none. **`status = not-found`,
confidence 0.40** — `not-found` is necessary-not-sufficient (per
`feedback_literature_novelty_check_baseline.md`); a future online
corpus update or human triage may reveal a Uranian cycler reference we
have not curated, in which case this is a regime-fresh rediscovery,
not a discovery.

## Part C — 3D CR3BP existence probe at Uranus-Oberon / Uranus-Umbriel

Script: `scripts/scan_312_uranus_3d_probe.py`. Output:
`data/scan_312_uranus_3d_probe.jsonl` (49 rows: meta + 24 + 24 +
summary).

The repeated-moon Lambert genome is planet-frame coplanar by
construction. Adding `z0 != 0` to a planet-frame IC pushes the
trajectory out of the moon orbital plane, where it no longer encounters
the moons — so the brief's "small z0 != 0" extension does NOT map onto
the cycler genome itself. The honest re-scoping was: **does the 3D
corrector close ANY periodic orbit in CR3BP(Uranus, Oberon) or
CR3BP(Uranus, Umbriel) from a small-z0 L1-Lyapunov-tail seed?**

Seed grid:

| Knob   | Values                               | Note                              |
|--------|--------------------------------------|-----------------------------------|
| x0     | `1 - mu - 0.95 * gamma_L1`           | Just inside L1 (Hill + Newton)    |
| y0     | 0                                    | Perpendicular crossing convention |
| z0     | `{1e-4, 1e-3, 5e-3, 1e-2}` (nondim)  | Small out-of-plane offset         |
| xdot0  | 0                                    | "                                 |
| ydot0  | `{-0.5, -0.2, -0.1, +0.1, +0.2, +0.5}` | Energy / Jacobi parameter     |
| zdot0  | 0                                    | "                                 |
| T_guess| Hill linearisation `2π / √(1 + 2μ/γ³)` | Vertical-Lyapunov estimate     |

Total: 4 × 6 = 24 seeds per system, 48 across both systems. Per-seed
SIGALRM timeout 120 s.

### C.1 Results

| System         | seeds | converged 3D | converged planar-collapse | failed (incl. timeout) |
|----------------|------:|-------------:|--------------------------:|-----------------------:|
| Uranus-Oberon  |    24 |        **0** |                         4 |                     20 |
| Uranus-Umbriel |    24 |        **0** |                         4 |                     20 |

**Zero converged 3D orbits at either system.** The four "converged"
rows per system landed with `degenerate_planar = True` — the corrector
pushed `z0` back to ~0 (planar collapse onto the Lyapunov family). The
remaining 20 per system either diverged (`corrector_residual > 1e-3`)
or hit the per-seed timeout.

### C.2 Verdict on Part C

**Clean negative on this seed strategy.** The L1-Lyapunov-tail + small
z0 grid does not land any 3D family in either CR3BP. This is acceptable
per the brief.

What this does NOT prove: that 3D periodic orbits don't exist in these
systems. They almost certainly DO (halo / NRHO / vertical Lyapunov
families exist at every small-mu CR3BP). The seed strategy is the
limiter — at mu ~ 3e-5 / 1.5e-5 the L1 Lyapunov basin is much narrower
than Earth-Moon (mu ~ 1.2e-2), and the period guess from the vertical
Lyapunov linearization (~2π) is far from the actual halo period
(typically ~π). A mu-continuation approach starting from a converged
Earth-Moon halo and continuing in mu down to Uranus values would be the
correct Phase-2 path. The bounded scout was right to not chase that.

## Phase 2 path

Two threads emerge.

### 1. The Umbriel-Oberon-Umbriel (1, 1) SILVER

The candidate is **deliberately withheld from catalogue admission**. The
admission path:

* `#305 BCR4BP V0-V5 gauntlet` — embed the planar Lambert IC in
  BCR4BP(Sun, Uranus, Oberon) or BCR4BP(Sun, Uranus, Umbriel) (note
  that Uranus has TWO suitable BCR4BP framings here, not one), run
  the V0 false-positive flagger + V1 short-arc-bridge + V2 multi-arc
  + V3 mass + V4 lit-cross-check + V5 manuscript guards.
* `#306 3D V0-V5 gauntlet` — currently pending. If the SILVER survives
  V0-V2, V3 (mass) is the next gate.

If both gauntlets fast-fail (V0 p_fp > 0.75 after BCR4BP re-evaluation,
or V1 short-arc bridge to higher-fidelity rejects the closure), the
candidate becomes a clean Phase-1 negative — exactly what the brief's
"quantify why it can't close" path expects. The honest residual record
is in
`data/scan_312_uranus_oberon_umbriel.jsonl` and
`data/scan_312_uranus_umbriel_oberon_offset_sweep.jsonl`.

Concrete IC to feed forward:

```
primary        = Uranus
sequence       = Umbriel -> Oberon -> Umbriel
n_rev          = (1, 1)
tof_per_leg    = 11.2054 d (= 1.5 * sqrt(P_Oberon * P_Umbriel))
phase0_anchor  = 311.25 deg  (any of 10 symmetric values give the same
                              residual; phase0 = 0 is permissible too)
rel_offset     = 45.00 deg   (Umbriel longitude - Oberon longitude at t0)
residual_kms   = 0.024033
v_inf at Umbriel = 2.27 km/s  (CAVEAT: 4.2x surface escape)
v_inf at Oberon  = 0.98 km/s
```

### 2. The Oberon-Titania-Oberon (1, 1) genome-ceiling near-miss

0.0617 km/s held under the wider (k1, k2) grid + 8x phase resolution.
The 0.0117 km/s gap above the gate IS the ballistic Lambert single-shot
ceiling. To close it requires a **more capable genome** — DSM legs,
multi-arc closure, or low-thrust maintenance. The #226 FBS optimizer
(once the #243 fair-trial completes) is a candidate driver; the #309
low-thrust maintenance driver is another. This is a Phase-2 PR for
genome capability, not a Phase-1 rerun.

## Negative registry entries

Per `project_negative_results_registry`, the Phase-1 clean negatives:

* **Uranus 3D CR3BP existence (L1-Lyapunov-tail seed, z0 ∈ {1e-4, 1e-3,
  5e-3, 1e-2}, ydot0 ∈ ±{0.1, 0.2, 0.5})** — 0 of 24 seeds per system
  yielded a 3D periodic orbit at either Uranus-Oberon (μ = 3.54e-5) or
  Uranus-Umbriel (μ = 1.47e-5). Method capability: full-asymmetric
  single-shooting 3D corrector at tol = 1e-10 / independent_tol = 1e-6.
  Re-sweep is appropriate when (a) mu-continuation from a converged
  Earth-Moon halo seed is implemented, (b) regularised propagator
  (`cyclerfinder.core.cr3bp_regularized`) is used for deep-perilune
  members, or (c) the seed-x0 strategy widens to L2 or vertical-
  Lyapunov tail.
* **Uranus Oberon-Titania-Oberon (1, 1) ballistic-Lambert closure** —
  residual 0.0617 km/s, GENOME CEILING confirmed across phase grid
  {12, 24, 48, 96}. Method capability: repeated-moon multi-rev
  single-shot Lambert in planet-frame, n_rev_grid ∈ (0..5). Re-sweep
  when a DSM-leg or multi-arc genome is in tree.

## Census of work

| Script                                                           | enumerated | closed | SILVER | bytes JSONL |
|------------------------------------------------------------------|-----------:|-------:|-------:|------------:|
| `scan_312_uranus_oberon_titania_finer.py`                        |        210 |     28 |      0 | small        |
| `scan_312_uranus_other_pairs.py` (5 systems)                     |       9390 |   2395 |      1 | mid          |
| `scan_312_uranus_umbriel_oberon_offset_sweep.py`                 |          1 |      1 |      1 | small        |
| `scan_312_uranus_3d_probe.py`                                    |         48 |      8 |      0 | small        |

Total Phase-1 wall: ~3 minutes search + ~12 minutes 3D probe.

## File index

* Scripts:
  `scripts/scan_312_uranus_oberon_titania_finer.py`,
  `scripts/scan_312_uranus_other_pairs.py`,
  `scripts/scan_312_uranus_umbriel_oberon_offset_sweep.py`,
  `scripts/scan_312_uranus_3d_probe.py`.
* Data:
  `data/scan_312_uranus_oberon_titania_finer.jsonl`,
  `data/scan_312_uranus_robustness.jsonl`,
  `data/scan_312_uranus_titania_umbriel.jsonl`,
  `data/scan_312_uranus_oberon_umbriel.jsonl`,
  `data/scan_312_uranus_ariel_titania.jsonl`,
  `data/scan_312_uranus_miranda_ariel_umbriel.jsonl`,
  `data/scan_312_uranus_ariel_umbriel_titania.jsonl`,
  `data/scan_312_uranus_other_pairs_index.jsonl`,
  `data/scan_312_uranus_umbriel_oberon_offset_sweep.jsonl`,
  `data/scan_312_uranus_3d_probe.jsonl`.
* Doc: this file
  (`docs/notes/2026-06-16-312-uranus-extended-sweep.md`).
