# #341 Neptune Proteus-Triton extended sweep — Phase 1 verdict

The #320 Vector B repeated-moon scan left a single tight near-miss at
Neptune: Proteus-Triton-Proteus (1, 1) at residual 0.0584 km/s — 16%
above the 0.05 km/s gate. NO published Neptunian-cycler prior exists in
the `KNOWN_CORPUS` (0 anchors mention any Neptunian body). This was the
session's strongest lit-fresh near-miss signal — and the analog of the
#312 → #339 Uranus workflow that produced the catalogue's first
computed quasi_cycler row (`umbriel-oberon-1-1-uranian-quasi-cycler`).

This Phase-1 task asked three honest questions, structurally identical
to #312:

* (A) Does a broader (k1, k2) grid + finer phase x rel_offset
  resolution close the 0.0584 km/s gap at Proteus-Triton-Proteus?
* (B) Do the other Neptune-system pairs (Triton-Nereid, Nereid-Triton,
  Proteus-Triton-Nereid) produce a closer near-miss?
* (C) Does the 3D corrector (#291) close any periodic orbit in
  CR3BP(Neptune, Triton) or CR3BP(Neptune, Proteus) from a small-z0
  Lyapunov-tail seed?

**Headline:** clean negative — no SILVER closure surfaces. The Part A
near-miss residual is a **genome ceiling** above the basin's Nyquist
resolution; Part B is **non-executable under sourced discipline**
(Nereid GM not determined); Part C **fails to produce any 3D variant**
(72 seeds, 0 non-planar converged). The 0.0584 km/s residual additionally
**fails the physical-sanity gate** even if it were sub-gate, because
Proteus's mass cannot deliver a useful flyby bend.

Frame: Neptune extended sweep (#341) — 70 + 4 + 72 = 146 cells / seeds
evaluated, best closure 0.05839 km/s, all 6 near-miss cells lit-fresh
class (KNOWN_CORPUS anchor overlap = 0), zero SILVER, zero physical-
gate passes.

## Discipline anchors

* READ-ONLY on the existing #254 repeated-moon genome
  (`src/cyclerfinder/search/discovery_campaign.py::RepeatedMoonTarget`),
  the #320 Vector B sweep
  (`scripts/scan_320_epoch_aware_moon_systems.py`), the #291 3D corrector
  (`src/cyclerfinder/search/cr3bp_general_periodic_3d.py`),
  `literature_check.py`, and the catalogue. Three new scripts in
  `scripts/scan_341_neptune_*.py` wrap existing routines without
  modification.
* NO catalogue writeback. The 0.0584 km/s near-miss stays out of the
  catalogue (not SILVER under any test).
* NO novelty claims in commit messages. Lit-fresh anchor overlap = 0
  is the necessary-not-sufficient gate from
  `feedback_literature_novelty_check_baseline.md`.
* Sourced moon-mu and SMA from
  `src/cyclerfinder/core/satellites.py` (JPL DE440 system GMs + sat441
  mean elements + sat-phys_par, accessed 2026-06-14). Neptune system GM
  6.836527100580e6 km³/s² (JPL DE440 gm_de440). Triton mu 1428.49546 /
  a 354800 km. Proteus mu 2.58342 / a 117600 km.
* `feedback_golden_tests_sourced_only`: Nereid GM is JPL SSD = 0.0
  (mass not determined), so we do not fabricate a value. Part B's
  Nereid-involving sweeps cannot run.

## Part A — Wider (k1, k2) grid + phase x rel_offset robustness

Script: `scripts/scan_341_neptune_proteus_triton_finer.py`. Outputs:
`data/scan_341_neptune_proteus_triton_finer.jsonl` (70 cells enumerated,
6 near-miss < 1 km/s, 0 SILVER),
`data/scan_341_neptune_robustness.jsonl` (4 phase-resolution rows).

The script directly ports `_sweep_one_cycle` from
`scripts/scan_320_epoch_aware_moon_systems.py` so the residuals are
identically comparable to the #320 baseline. Using
`RepeatedMoonTarget.close` would have been wrong: it freezes the
relative phase offset at `2π / len(consts)` (the deterministic
moon-index offset) and only sweeps the global phase, which misses the
240/255 deg minimum where #320's 0.0584 km/s lives.

### A.1 Wider (k1, k2) grid

`n_rev_grid = (0, 1, 2, 3, 4, 5)` over Neptune (Triton, Proteus) — 2
length-3 cycles × 35 (k1, k2) cells (excluding the trivial (0, 0)) = 70
cells. 6 near-miss < 1 km/s, all centered at `tof_scale=2.0`. Best 5:

| Sequence                  | n_rev | residual (km/s) | phase0° | rel_off° | max bend° | phys gate |
|---------------------------|:------|----------------:|--------:|---------:|----------:|:---------:|
| Proteus-Triton-Proteus    | (1,1) |          0.0584 |     240 |      255 |      38.0 |   FAIL    |
| Proteus-Triton-Proteus    | (1,0) |          0.1725 |     240 |       60 |       3.9 |   FAIL    |
| Triton-Proteus-Triton     | (1,1) |          0.2251 |     285 |      105 |      10.4 |   FAIL    |
| Triton-Proteus-Triton     | (1,0) |          0.5172 |     345 |      240 |      12.1 |   FAIL    |
| Triton-Proteus-Triton     | (0,1) |          0.5848 |     255 |      330 |      12.4 |   FAIL    |

No (k1, k2) cell beyond #320's `(0..3)` grid produced a sub-gate
closure or even a new near-miss < 0.058. **Wider grid does not close
the 0.058 → 0.05 km/s gap.** This matches the Part A.1 outcome at
Uranus Oberon-Titania (#312).

### A.2 Phase x rel_offset robustness on (1, 1)

Re-run Proteus-Triton-Proteus (1, 1) at `n_phase = n_offset ∈ {12, 24,
48, 96}` (the #320 Vector B basin convention):

| ps  | residual (km/s)              | phase0° | rel_off° |
|:---:|:----------------------------:|--------:|---------:|
| 12  | 0.5824284648253311           |     120 |      240 |
| 24  | 0.058387519970644064         |     240 |      255 |
| 48  | 0.05838751997063496          |     127 |      255 |
| 96  | 0.05838751997063496          |     127 |      255 |

The ps=12 row is below the basin's Nyquist (30° grid spacing misses
the 240/255 deg minimum entirely; the next-best grid point produces
0.582 km/s). ps=24, 48, 96 land at 0.05838751997 km/s flat to 1e-14.
**Verdict above the resolution floor: GENOME CEILING** — finer phase
× rel_offset grid does NOT move the residual; the 0.058 km/s near-miss
is intrinsic to the Lambert single-shot coplanar genome at this pair /
(k1, k2). The phase-resolution-PARTIAL verdict-row label in the JSONL
reflects the spread including ps=12, but the resolved-state verdict is
the same as #312 / #285: ceiling.

Note: the residual at ps=24 vs ps≥48 differs in the 14th significant
digit (0.058387519970644 vs 0.058387519970635) — well below any
physical interpretation threshold; both round to the same 5-decimal
near-miss.

### Critical complementary finding: physical-sanity gate FAILS at the best

The 0.0584 km/s near-miss row has `physical_gate_passed = False`. The
Triton flyby itself delivers a healthy 38.0° bend at V_inf ≈ 1.43 km/s
— but the Proteus flybys at V_inf ≈ 1.93 / 1.92 km/s only bend the
trajectory 0.24° / 0.26°, far below the 5° minimum-useful-bend floor.
Proteus's GM is 2.58 km³/s² vs Triton's 1428 km³/s² — a ratio of 553×.
Even if a basin search closed below 0.05 km/s, the resulting
trajectory would not be a physically useful flyby tour: Proteus is too
small to be the anchor body of a Triton-cycler. The brief flagged this
("the GM ratio Triton/Proteus is ~550 so a Proteus flyby contributes
negligible bending"); the verdict-row data confirms it quantitatively.

## Part B — Other Neptune-system pair / 3-body sweeps

Script: `scripts/scan_341_neptune_other_pairs.py`. Output:
`data/scan_341_neptune_other_pairs.jsonl` (1 verdict row).

Three of four briefed configurations (Triton-Nereid, Nereid-Triton,
Proteus-Triton-Nereid) require Nereid, which is INTENTIONALLY OMITTED
from `src/cyclerfinder/core/satellites.py` because JPL SSD lists its GM
as 0.0 (mass not determined; see `satellites.py:175-176`). Per
`feedback_golden_tests_sourced_only` and the task discipline ("sourced
golden for Neptune system mu values"), the Nereid GM is not fabricated.
These sweeps are recorded as a registered empty region (per
`project_negative_results_registry`).

The fourth configuration (Triton-Proteus) is the reverse of Part A's
moon set; the Part A `_sweep_one_cycle` enumerator covers BOTH closed
length-3 cycles on this pair — Proteus-Triton-Proteus AND
Triton-Proteus-Triton — so it is fully captured by
`data/scan_341_neptune_proteus_triton_finer.jsonl` (best
Triton-Proteus-Triton residual 0.2251 km/s at (1,1)).

**Reactivation condition for the empty region:** re-sweep when Nereid
GM becomes a sourced quantity (JPL SSD updated above 0.0, or an
independent peer-reviewed determination). Additionally, Nereid's
e=0.75 makes the circular-coplanar moon-orbit assumption used by the
current genome unsuitable even with a sourced GM; extending the genome
to eccentric moons is its own task.

## Part C — 3D CR3BP existence probe at Neptune-Triton / Neptune-Proteus

Script: `scripts/scan_341_neptune_3d_probe.py`. Output:
`data/scan_341_neptune_3d_probe.jsonl` (72 seed rows, 1 verdict).

72 seeds (6 z0 ∈ {1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1} × 6 ydot0 ∈
{-0.5, -0.2, -0.1, 0.1, 0.2, 0.5}) at the L1-Lyapunov tail
(x0 = x_L1 − 0.95 γ) in CR3BP(Neptune, Triton) and CR3BP(Neptune,
Proteus). Period guess from the Hill-linearized vertical-Lyapunov
formula. Full asymmetric 3D corrector with `FREE_VARS_FULL_ASYMMETRIC`.

| System          | n_seeds | converged 3D | converged planar | failed | best 3D ind_res |
|-----------------|--------:|-------------:|-----------------:|-------:|----------------:|
| Neptune-Triton  |      36 |        **0** |               23 |     13 |             N/A |
| Neptune-Proteus |      36 |        **0** |                1 |     35 |             N/A |

**Verdict: NO 3D periodic orbit closes** from the L1-Lyapunov-tail seed
convention. The corrector collapses every converged orbit to the planar
manifold (z0 → 0). This matches the #312 Uranus Part C outcome
(`scan_312_uranus_3d_probe.jsonl`), where the same seed convention also
produced 0 / 48 non-planar convergences at CR3BP(Uranus, Oberon) and
CR3BP(Uranus, Umbriel).

The Triton system has uniform converged-planar Jacobi ≈ 2.99979 across
the 23 collapses (well-defined planar Lyapunov family at the seed
energy); the Proteus system at much smaller mu has mostly failures (35
/ 36) because the L1 Hill seed is closer to the primary and steeper.

**Critical caveat:** this is a NARROW question (do 3D orbits exist
from this specific seed?). A negative does NOT prove that 3D periodic
orbits do not exist in CR3BP(Neptune, Triton). To rule out the family
broadly would require continuation from a known planar Lyapunov member
into the vertical-bifurcation tangent, which is Phase 2 territory.

**Answer to the brief's "does inclining the orbit close the residual?"
question: NO** — at the L1-Lyapunov-tail seed convention used here, no
3D variant is produced, so none can be tested against the 0.058 km/s
gap.

## Phase 2 path

The original analog (#312 Uranus) found its SILVER in **Part B** —
on a DIFFERENT moon pair (Oberon-Umbriel rather than the original
near-miss Oberon-Titania). At Neptune, Part B is **structurally
non-executable** within sourced discipline: we have only two registered
Neptune moons, the briefed alternative pair (any Nereid pairing) is
data-blocked. The 3-body length-5 bands that #312 Part B also covered
require a third moon, which doesn't exist for Neptune in the registry.

Two remaining levers, both Phase 2:

1. **Low-thrust maintenance extension** (#309 machinery). The 0.0584
   km/s − 0.05 km/s = 0.008 km/s gap is plausibly closable with sub-mm/s
   continuous-thrust maintenance. BUT — the physical-sanity gate
   already fails on the 0.0584 km/s configuration (Proteus too weak to
   bend, 0.24° vs 5° floor). Low-thrust maintenance cannot fix an
   un-physical flyby tour; it can only close residual on a tour where
   each flyby is physically genuine. **Recommend: skip the low-thrust
   probe at this configuration; it doesn't address the binding
   constraint.**

2. **Acquire Nereid GM** (data-acquisition Phase 2 lever, parallel to
   the literature-acquisitions work). If a sourced Nereid GM emerges
   (peer-reviewed determination from gravitational signature or stellar
   occultation), the Part B Triton-Nereid sweep becomes available. The
   high-eccentricity (e=0.75) makes the current circular-coplanar
   genome unsuitable, so this also implies a genome extension. Two-step
   dependency.

Without one of these, Neptune is a **clean Phase-1 negative** for
ballistic-coplanar quasi-cycler discovery. The Triton retrograde
geometry remains structurally interesting for inclination-trade work
but isn't reachable by the current genome.

## Pointers

* Scripts:
  * `scripts/scan_341_neptune_proteus_triton_finer.py`
  * `scripts/scan_341_neptune_other_pairs.py`
  * `scripts/scan_341_neptune_3d_probe.py`
* Data:
  * `data/scan_341_neptune_proteus_triton_finer.jsonl`
  * `data/scan_341_neptune_robustness.jsonl`
  * `data/scan_341_neptune_other_pairs.jsonl`
  * `data/scan_341_neptune_3d_probe.jsonl`
* Substrate (READ-ONLY): `src/cyclerfinder/core/satellites.py`,
  `src/cyclerfinder/search/discovery_campaign.py`,
  `src/cyclerfinder/search/cr3bp_general_periodic_3d.py`,
  `src/cyclerfinder/search/physical_sanity.py`,
  `src/cyclerfinder/search/literature_check.py`.
* Reference baseline: `scripts/scan_320_epoch_aware_moon_systems.py`
  + `data/scan_320_epoch_aware_neptune.jsonl` (commit 0e6f3f2).
* Reference Uranus analog: `docs/notes/2026-06-16-312-uranus-extended-sweep.md`,
  `docs/notes/2026-06-16-327-umbriel-silver-verification.md`,
  `docs/notes/2026-06-17-339-silver-quasi-cycler-admission.md`.
