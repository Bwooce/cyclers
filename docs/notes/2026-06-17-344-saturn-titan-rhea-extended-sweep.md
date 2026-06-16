# #344 Saturn Titan-Rhea extended sweep — Phase 1 verdict

The #320 Vector B repeated-moon scan flagged ONE SILVER hit at Saturn:
Titan-Rhea-Titan (1, 1) at residual 0.0316 km/s — below the 0.05 km/s
gate AND physical-sanity PASS (Titan flybys bend the trajectory
~50 deg at V_inf ≈ 1.72 km/s, Rhea flyby bends 6.8 deg at V_inf ≈ 1.68
km/s, all above the 5 deg minimum-useful floor). But the KNOWN_CORPUS
anchor-overlap count was already 1 at scan time (Davis-Phillips-McCarthy
Saturn tulip-shaped orbits, Acta Astronautica 143:16-28, 2018), so the
candidate was NOT lit-fresh class even before #334's corpus update.

This Phase-1 task asked four honest questions, structurally analogous
to #312 → #339 (Uranus Umbriel-Oberon SILVER that produced the
catalogue's first computed quasi_cycler row) and #341 (Neptune
Proteus-Triton, which became a clean negative):

* (A) Does a broader (k1, k2) grid + finer phase × rel_offset
  resolution close DEEPER than the 0.0316 km/s baseline?
* (B) Do the other Saturn-system pairs (Titan-Iapetus, Titan-Dione,
  Titan-Tethys, Rhea-Dione, 3-body Titan-Rhea-Dione) produce a
  closer SILVER, ideally lit-fresh?
* (C) Does the 3D corrector (#291) close any periodic orbit in
  CR3BP(Saturn, Titan) or CR3BP(Saturn, Rhea) from a small-z0
  Lyapunov-tail seed, opening an inclination-trade lever?
* (D) Post-#334 lit-recheck on the strongest hit (whichever it is):
  is it still structural lit-fresh, or is it rediscovery, or is it
  inconclusive (matched on body-set but not on (k1,k2) topology or
  V_inf tuple)?

**Headline:**

* Part A.1 (wider grid): NO new (k1,k2) cell closes below the SILVER
  Titan-Rhea-Titan (1,1) baseline at the 24×24 phase grid.
* **Part A.2 (phase robustness): the residual DROPS to 0.0102 km/s at
  ps=96** — the 0.0316 km/s SILVER was a phase-resolution artifact, not
  the genome ceiling. The cell closes ~3× deeper than the #320 SILVER.
* Part B (other Saturn pairs + 3-body): **0 NEW SILVER**. 5 pair/3-body
  sweeps (280 length-3 cells + 4608 length-5 cells = 4888 cells), best
  per pair (km/s): Titan-Tethys 0.0757 (lit-fresh, physical-gate FAIL),
  Rhea-Dione 0.107 (anchored), Titan-Dione 0.111 (anchored),
  Titan-Iapetus 0.231 (anchored), Titan-Rhea-Dione 3-body 0.223
  (anchored). None reaches the 0.05 km/s gate.
* Part C (3D existence probe): **0 / 72 seeds produce a 3D periodic
  orbit.** Saturn-Titan 0 / 36 (24 planar collapse, 12 fail);
  Saturn-Rhea 0 / 36 (0 planar, 36 fail — small μ basin too narrow).
  Same outcome as #312 Uranus and #341 Neptune Part C; the
  L1-Lyapunov-tail seed convention does not produce non-planar
  convergence in the Saturnian inner satellite systems.
* Part D (post-#334 lit-recheck): the 0.0102 km/s Titan-Rhea-Titan
  (1, 1) closure is at the SAME cell as #320 (phase≈97.5°,
  rel_off≈292.5°, tof_scale=2.0). KNOWN_CORPUS anchor overlap at the
  post-#334 corpus state is 2: Davis-Phillips-McCarthy 2018 (body_set
  ⊇ {Titan, Rhea}) and Cassini-Huygens Saturn-Titan satellite tour
  (body_set ⊇ {Titan, Rhea}). Both are body-set-only matches, not
  (k1,k2)-topology or V_inf-tuple matches; **but body-set match is the
  discipline's necessary-not-sufficient gate, so the candidate is NOT
  lit-fresh class** under `feedback_literature_novelty_check_baseline`.

Frame: Saturn extended sweep (#344) — 70 (A.1) + 4 (A.2) + 4888 (B) +
72 (C) = **5034 cells / seeds evaluated**, best closure 0.0102 km/s,
lit-anchor overlap = 2, **zero new admissible SILVER**.

## Discipline anchors

* READ-ONLY on the existing #254 repeated-moon genome
  (`src/cyclerfinder/search/discovery_campaign.py::RepeatedMoonTarget`),
  the #320 Vector B sweep
  (`scripts/scan_320_epoch_aware_moon_systems.py`), the #291 3D corrector
  (`src/cyclerfinder/search/cr3bp_general_periodic_3d.py`),
  `literature_check.py`, `physical_sanity.py`, and the catalogue. Three
  new scripts in `scripts/scan_344_saturn_*.py` wrap existing routines
  without modification.
* NO catalogue writeback. The 0.0316 → 0.0102 km/s closure does NOT
  promote to catalogue: it is anchored on body-set (the necessary-not-
  sufficient gate of `feedback_literature_novelty_check_baseline`), so
  even though it passes the SILVER residual + physical-sanity gates, it
  fails the lit-fresh prerequisite for a #339-class admission.
* NO novelty claims in commit messages. Anchor overlap = 2 (Davis-
  Phillips-McCarthy 2018, Cassini-Huygens tour) — at minimum the
  candidate's body set is the published Saturn-moon repeating-tour
  archetype.
* Sourced moon-mu and SMA from `src/cyclerfinder/core/satellites.py`
  (JPL DE440 system GMs + sat441 mean elements + JPL SSD phys_par,
  accessed 2026-06-14). Saturn system GM 3.7931207e7 km³/s² (NASA
  Saturn fact sheet / JPL DE440 gm_de440). Titan GM 8978.14 / a
  1221870 km. Rhea GM 153.94 / a 527070 km. Tethys GM 41.21 / a
  294670 km. Dione GM 73.116 / a 377420 km. Iapetus GM 120.51511 / a
  3561700 km (JPL SSD phys_par SAT441).

## Part A — Wider (k1, k2) grid + phase × rel_offset robustness

Script: `scripts/scan_344_saturn_titan_rhea_finer.py`. Outputs:
`data/scan_344_saturn_titan_rhea_finer.jsonl` (70 cells enumerated, 6
near-miss < 1 km/s, 1 SILVER),
`data/scan_344_saturn_robustness.jsonl` (4 phase-resolution rows).

The script directly ports `_sweep_one_cycle` from
`scripts/scan_320_epoch_aware_moon_systems.py` so residuals are
identically comparable to the #320 baseline. Smoke-tested at the SILVER
cell: 0.03161954212289819 km/s reproduced bit-for-bit at the published
phase0=90°, rel_off=285°, tof_scale=2.0, V_infs [1.7528, 1.6769,
1.7212] km/s.

### A.1 Wider (k1, k2) grid

`n_rev_grid = (0, 1, 2, 3, 4, 5)` over Saturn (Titan, Rhea) — 2
length-3 cycles × 35 (k1, k2) cells (excluding the trivial (0, 0)) = 70
cells. 6 near-miss < 1 km/s, all centered at tof_scale=2.0. 1 SILVER
(Titan-Rhea-Titan (1,1) at 0.0316 km/s, unchanged from #320). Best 5:

| Sequence              | n_rev | residual (km/s) | phase0° | rel_off° | max bend° | phys gate |  anchors  |
|-----------------------|:------|----------------:|--------:|---------:|----------:|:---------:|:---------:|
| Titan-Rhea-Titan      | (1,1) |          0.0316 |    90.0 |    285.0 |      50.5 |    PASS   |     2     |
| Rhea-Titan-Rhea       | (1,1) |          0.5582 |    30.0 |     75.0 |      19.8 |    FAIL   |     2     |
| Titan-Rhea-Titan      | (0,1) |          0.5848 |   225.0 |    330.0 |      21.0 |    FAIL   |     2     |
| Titan-Rhea-Titan      | (1,0) |          0.7158 |    15.0 |     75.0 |       6.7 |    FAIL   |     2     |
| Rhea-Titan-Rhea       | (1,0) |          0.8356 |   225.0 |     30.0 |      11.8 |    FAIL   |     2     |

No (k1, k2) cell beyond #320's (0..3) grid produced a sub-SILVER
closure at the 24×24 phase grid. **Wider grid alone does not deepen
beyond the published 0.0316 km/s SILVER** — same outcome shape as
the #341 Neptune A.1 / the #312 Uranus A.1.

### A.2 Phase × rel_offset robustness on (1, 1)

Re-run Titan-Rhea-Titan (1, 1) at `n_phase = n_offset ∈ {12, 24, 48,
96}` (the #320 Vector B basin convention):

| ps  | residual (km/s)              | phase0° | rel_off° | tof_scale |
|:---:|:----------------------------:|--------:|---------:|:---------:|
| 12  | 0.0542909979365096           |   240.0 |    300.0 |    2.0    |
| 24  | 0.03161954212289819          |    90.0 |    285.0 |    2.0    |
| 48  | 0.011163625043707937         |    97.5 |    292.5 |    2.0    |
| 96  | 0.010188096573990224         |   273.7 |    288.8 |    2.0    |

**Verdict: PHASE_RESOLUTION_CLOSED_GAP** — spread 0.0441 km/s across
{12, 24, 48, 96}, but the ps=48 and ps=96 rows both close at residual
below 0.0316 km/s (the published SILVER) and well below the gate. The
0.0316 km/s residual at ps=24 was a basin-sampling artifact, NOT the
genome ceiling.

This differs fundamentally from the #341 Neptune outcome
(GENOME_CEILING — spread < 1e-13 across {48, 96}). At Saturn
Titan-Rhea, **the basin is narrower than the 24×24 phase grid
resolves**, and finer sampling closes ~3× deeper. The basin floor at
ps=96 (≈0.0102 km/s) is still positive — the cell does not close
exactly ballistically in Lambert single-shot — but it is comfortably
inside the SILVER gate.

The best-cell phase shifts from (90°, 285°) at ps=24 to (97.5°, 292.5°)
at ps=48 to (273.7°, 288.8°) at ps=96: ps=48 and ps=96 are close in
rel_off but the global phase0 differs by ≈ π — the basin has two
near-equal-depth minima at 90° offset by a half-cycle of the global
phase, which is the expected symmetry of the Titan-Rhea-Titan cycle.
The 0.0102 km/s ≈ 0.011 km/s drift between ps=48 and ps=96 indicates
the basin is still not fully resolved at ps=96 (a finer grid would
likely close another ~10-20% deeper but not change the order of
magnitude).

### V_inf continuity across resolutions

| ps  | V_inf [Titan, Rhea, Titan] (km/s)     | tof_days (per leg) |
|:---:|---------------------------------------|-------------------|
| 12  | [1.7198, 1.6876, 1.7741]              | 16.977            |
| 24  | [1.7528, 1.6769, 1.7212]              | 16.977            |
| 48  | [1.7269, 1.6383, 1.7381]              | 16.977            |
| 96  | [1.7375, 1.6463, 1.7273]              | 16.977            |

The V_inf tuple is stable across resolutions (Titan flybys 1.72-1.78
km/s, Rhea flyby 1.64-1.69 km/s); the residual reduction is from the
phase-sampling resolving the basin floor, not from a different physical
configuration. Same TOF (16.977 days/leg, tof_scale=2.0) at every
resolution — both legs are at the geometric-mean of the two moons'
periods, doubled, which is also the synodic-period basin used by #320.

Discipline note: a 0.0102 km/s closure at the same physical
configuration as the #320 SILVER (same cell, same V_infs to first
order, same flyby bends to first order) is still a *finer measurement
of the same basin*, not a new candidate. The lit-recheck verdict
applies identically.

Discipline note: a 0.0102 km/s closure at the same physical
configuration as the #320 SILVER (same cell, same V_infs to first
order, same flyby bends to first order) is still a *finer measurement
of the same basin*, not a new candidate. The lit-recheck verdict
applies identically.

## Part B — Other Saturn-system pair / 3-body sweeps

Script: `scripts/scan_344_saturn_other_pairs.py`. Outputs per pair:

| Pair                    | Cells | Best residual (km/s) | Best cycle             | n_rev | V_inf tuple (km/s)                        | phys gate | anchors | lit-fresh |
|-------------------------|------:|---------------------:|------------------------|-------|-------------------------------------------|:---------:|:-------:|:---------:|
| Titan-Iapetus           |    70 |               0.2306 | Titan-Iapetus-Titan    | (1,1) | [5.05, 1.79, 5.28]                        |   FAIL    |    1    |    NO     |
| Titan-Dione             |    70 |               0.1108 | Dione-Titan-Dione      | (1,1) | [8.45, 3.44, 8.34]                        |   FAIL    |    2    |    NO     |
| Titan-Tethys            |    70 |           **0.0757** | Titan-Tethys-Titan     | (1,1) | [2.50, 6.68, 2.51]                        |   FAIL    |  **0**  |  **YES**  |
| Rhea-Dione              |    70 |               0.1069 | Rhea-Dione-Rhea        | (1,1) | [1.11, 2.21, 1.11]                        |   FAIL    |    2    |    NO     |
| Titan-Rhea-Dione 3-body |  4608 |               0.2229 | Dione-Titan-Dione-Titan-Dione | (1,1,1,1) | [8.45, 3.44, 8.34, 3.48, 8.23] |   FAIL    |    2    |    NO     |

Notes:

* **Titan-Tethys is the ONLY lit-fresh-class pair**: Davis-Phillips-
  McCarthy body_set is {Titan, Enceladus, Rhea, Dione}, Cassini-Huygens
  is {Titan, Enceladus, Rhea, Dione, Iapetus}; neither includes Tethys.
  But the best residual 0.0757 km/s is **above the 0.05 km/s gate**,
  AND physical-gate FAILS (Titan flybys bend 30.2° / 30.0° = PASS,
  Tethys flyby at V_inf 6.68 km/s bends only 0.17°, far below the 5°
  floor — Tethys's GM 41.2 km³/s² is too low to bend a 6.7 km/s
  trajectory). Even if a finer phase-resolution sweep closed the
  0.0076 km/s gap to 0.05, the Tethys flyby would not deliver a
  physically useful tour. Lit-fresh near-miss but not promotable.
* **Titan-Iapetus**: NEW coverage (Iapetus is NOT in #320's set).
  Iapetus is in Cassini-Huygens body_set but not Davis-Phillips-
  McCarthy. Best 0.231 km/s is far above gate; the Iapetus flyby V_inf
  ≈ 1.79 km/s is low (Iapetus is a slow-moving outer regular satellite
  with a ~79-day period; the geometry is hostile to a Titan-Iapetus
  ballistic cycler). Anchored on Cassini tour (Cassini did 1 Iapetus
  flyby in 2007).
* **Rhea-Dione**: best 0.1069 km/s reproduces #285's 0.10688 km/s
  baseline to 5 significant digits — finer methodology (#320 Vector B
  basin convention vs #285's coarser 12-sample grid) does NOT close
  the residual further. Same cell (Rhea-Dione-Rhea (1,1), V_inf 1.11 /
  2.21 / 1.11 km/s). Anchored.
* **Titan-Dione** (revisit of #320 at (0..5)): best 0.111 km/s at
  Dione-Titan-Dione (1,1). The Titan flybys are at V_inf 8.45 / 8.34
  km/s (genuinely high-energy Saturn-Titan encounters in the Cassini
  band); Dione flyby at 3.44 km/s is intermediate. Anchored.
* **3-body Titan-Rhea-Dione** (length-5): best 0.223 km/s at
  Dione-Titan-Dione-Titan-Dione (1,1,1,1) — a 4-leg cycle that goes
  Dione → Titan → Dione → Titan → Dione. This is essentially the
  Titan-Dione pair pumped twice; Rhea is not used in the best cycle.
  The 18 length-5 closed cycles on {Titan, Rhea, Dione} include
  combinations like Titan-Rhea-Dione-Rhea-Titan and
  Rhea-Titan-Dione-Titan-Rhea, but the deepest closure surfaces on the
  pair-pumped Titan-Dione cycle, not on a true 3-body topology.
  Anchored.

**Part B verdict: ZERO new SILVER, ZERO admissible lit-fresh candidate
in the Saturn-other-pair territory.** The only lit-fresh pair
(Titan-Tethys) clears the body-set gate but fails the residual gate by
50% AND fails the physical-sanity gate (Tethys can't bend).

## Part C — 3D CR3BP existence probe at Saturn-Titan / Saturn-Rhea

Script: `scripts/scan_344_saturn_3d_probe.py`. Output:
`data/scan_344_saturn_3d_probe.jsonl` (72 seed rows + 1 verdict row).

72 seeds (6 z0 ∈ {1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1} × 6 ydot0 ∈
{-0.5, -0.2, -0.1, 0.1, 0.2, 0.5}) at the L1-Lyapunov-tail
(x0 = x_L1 − 0.95 γ) in CR3BP(Saturn, Titan) and CR3BP(Saturn, Rhea).
Full asymmetric 3D corrector with `FREE_VARS_FULL_ASYMMETRIC`.

| System         | n_seeds | converged 3D | converged planar | failed | best 3D ind_res |
|----------------|--------:|-------------:|-----------------:|-------:|----------------:|
| Saturn-Titan   |      36 |        **0** |               24 |     12 |             N/A |
| Saturn-Rhea    |      36 |        **0** |                0 |     36 |             N/A |

**Verdict: NO 3D periodic orbit closes** from the L1-Lyapunov-tail seed
convention at either Saturn-Titan or Saturn-Rhea. The corrector
collapses every converged orbit on Saturn-Titan to the planar manifold
(z0 → 0); at Saturn-Rhea every seed fails outright.

This matches the same pattern as #312 Uranus and #341 Neptune Part C
(both 0 / 48 and 0 / 72 non-planar convergences respectively). The
Saturn-Titan basin has uniform Jacobi ≈ 2.99976 across the 24 collapsed
seeds (well-defined planar Lyapunov family at the seed energy). Saturn-
Rhea's mu (4.06e-6) puts L1 only γ ≈ 0.0114 units from Rhea, and the
seed convention's `x0 = x_L1 − 0.95γ` is too close to the secondary for
the corrector to stabilise — every seed yields a non-converged result.

**Critical caveat (same as #341):** this is a NARROW question (do 3D
orbits exist from this specific Lyapunov-tail seed?). The negative does
NOT prove that 3D periodic orbits don't exist in CR3BP(Saturn, Titan).
Brinckerhoff & Howell AAS 09-129 documents Saturn-Titan halo and
vertical-Lyapunov families; their seed convention is bifurcation-
continuation from a known planar Lyapunov member into the vertical-
bifurcation tangent, which is Phase 2 territory. The Phase 1 question
"does inclining the orbit close the 0.0102 km/s residual?" is **NO**
at the L1-Lyapunov-tail seed convention — no 3D variant is produced
at all, so none can be tested against the basin floor.

## Part D — Post-#334 lit-recheck on the Titan-Rhea-Titan (1,1) closure

The 0.0102 km/s closure at the same cell as the #320 SILVER is the
strongest Saturn signal in this sweep. The post-#334 KNOWN_CORPUS state
returns **2 anchors** matching its structural signature
(primary=Saturn, body_set={Titan, Rhea}):

1. **Davis-Phillips-McCarthy "Tulip-shaped orbits in the Saturn system,"
   Acta Astronautica 143:16-28 (2018), DOI 10.1016/j.actaastro.2017.11.011.**
   body_set = {Titan, Enceladus, Rhea, Dione}. Keywords: "Saturn
   tulip-shaped orbit / Titan period-multiplying orbit / Saturnian
   periodic orbit Np petal." This is a CR3BP-class family (petal-count
   periodic orbits in the planet-moon CR3BP), NOT a moon-flyby
   ballistic cycler. Match is BODY-SET ONLY — Davis 2018 does not
   document Lambert-stitched Titan-Rhea-Titan (k1=k2=1) cycles. The
   #320 doc (line 130-132) and the #285 doc (line 105-107) both record
   the match as a structural-fingerprint hit, not confirmed
   rediscovery.

2. **Cassini-Huygens Saturn-Titan satellite tour design**
   (Strange, Goodson, Yam, Buffington — JGCD / AAS 2010-2017; Goodson
   et al., JGCD 2008). body_set = {Titan, Enceladus, Rhea, Dione,
   Iapetus}. The published archetype is a Titan-pumped MGA tour
   (~127 Titan flybys, 9 Enceladus, 4 Rhea, 2 Dione, 1 Iapetus over
   the Cassini mission); Rhea was a one-off tour target from a Titan
   flyby, NOT a repeating Titan-Rhea-Titan ballistic cycler.
   Match is BODY-SET ONLY (necessary not sufficient).

**Lit-recheck depth analysis (analog of the #339 Umbriel-Oberon
discussion):**

| Match depth                                                    | Verdict        |
|----------------------------------------------------------------|----------------|
| (primary, body_set) — every Saturn-Titan-Rhea candidate hits   | 2 anchors      |
| (k1, k2) topology — (1, 1) two-leg pump cycler                 | NOT documented |
| V_inf tuple at IC (Titan 1.72-1.75 km/s, Rhea 1.68 km/s)       | NOT documented |
| Period (≈33.95 days), Jacobi, family continuation              | NOT documented |
| Independent cross-check at closure ≤ 1e-6 nondim               | NOT attempted  |

The #339 Umbriel-Oberon SILVER also matched only on body_set
(Earth-Moon CR3BP corpus anchored {Moon} for the Tisserand-family
adjacent literature), but its structural fingerprint at (k1, k2) and
V_inf was NOT documented in the matched anchor — so under the
discipline of `feedback_literature_novelty_check_baseline.md` the
candidate was admitted as lit-fresh class (necessary-not-sufficient
satisfied), then went through V0-V5 gauntlet and ratched into
catalogue at #340.

**The Titan-Rhea-Titan case is materially different from #339 in TWO
respects:**

1. **Anchor overlap COUNT is 2, not 0.** #339's `_candidate_anchors`
   query returned zero matches; the Saturn case returns two. The
   discipline (line 16 of `feedback_literature_novelty_check_baseline`)
   says "no cycler called novel until search/literature_check.py clears
   it vs the PUBLISHED record" — at #339 the matcher returned 0; at
   #344 it returns 2. The *prerequisite* for a novelty claim is unmet
   regardless of match depth.

2. **The matched anchors include the Cassini-Huygens published mission
   tour at the exact body set.** Cassini-Huygens IS a published Saturn
   Titan-and-Rhea-included MGA tour with extensive design heritage
   (Strange et al., Goodson et al., JGCD 2008-2017). Even if the (1,1)
   topology specifically isn't a published Cassini tour element, the
   structural fingerprint *is* the Cassini archetype's body set.
   `_candidate_anchors`'s `seq_set <= anchor.body_set` test (line 954
   of literature_check.py) is exactly the structural-fingerprint test
   the discipline defines.

**Verdict: the Titan-Rhea-Titan (1, 1) closure is NOT lit-fresh class
under the project discipline.** It is a body-set anchored candidate
under TWO independent published Saturn-system literature streams
(Davis-Phillips-McCarthy 2018 tulip family + Cassini-Huygens tour
design). The 0.0102 km/s closure is a genuine basin-floor measurement
deeper than the #320 SILVER, but it does not promote to a #339-class
catalogue admission.

This is the same disposition as the #341 erratum on Neptune
Proteus-Triton-Proteus (1,1): once the matcher fires, the candidate is
no longer in the lit-fresh class regardless of residual depth, and the
documentation must reflect the post-corpus state honestly.

### Negative-results registry entry

Per `project_negative_results_registry`, the empty/anchored region for
this method is:

* **Region:** Saturn moon-flyby ballistic quasi-cyclers, body set ⊆
  {Titan, Enceladus, Rhea, Dione, Iapetus} (the Cassini tour archetype
  body set). 
* **Method:** Lambert single-shot coplanar 2-leg repeated-moon genome
  (#254 RepeatedMoonTarget), 24×24 → 96×96 basin grid, n_rev ∈ (0..5),
  tof_scale ∈ {0.5, 1.0, 1.5, 2.0}, physical-sanity gate (min useful
  bend 5°), KNOWN_CORPUS gate at post-#334 corpus state (sha db54476).
* **Outcome:** Best closure 0.0102 km/s (Titan-Rhea-Titan (1,1) at
  ps=96). SILVER residual + physical-sanity PASS, but body-set
  anchored to Davis-Phillips-McCarthy 2018 + Cassini-Huygens tour
  design. NOT lit-fresh class. NO catalogue admission.
* **Reactivation conditions:** (a) acquire and read Davis-Phillips-
  McCarthy 2018 to confirm whether Titan-Rhea-Titan (1,1) is
  specifically documented or only the Saturn-Titan tulip family; if
  Davis 2018 only covers the tulip family (a different orbit class),
  the anchor is conservative and could be tightened to body_set =
  {Titan} only; this would unblock the Titan-Rhea-Titan signature for
  lit-fresh class. (b) An equivalent acquisition of the Strange et al.
  2010-2017 Cassini Titan-pump tour design papers to confirm whether
  Titan-Rhea-Titan (1,1) is in the published Cassini tour design space
  or only Titan-Titan resonance hops. (c) An independent peer-reviewed
  paper specifically titled / abstracted as Titan-Rhea-Titan repeating
  ballistic cycler in the published literature.

## Phase 2 path

* **The 0.0102 km/s closure does NOT promote to V0-V5 gauntlet** —
  lit-fresh prerequisite unmet. Per discipline, V2/V3/V4 are
  rediscovery-checking guards; the lit-anchor gate is the gate that
  rules a candidate IN to the gauntlet in the first place.

* **Targeted acquisition** of Davis-Phillips-McCarthy 2018 (Acta
  Astronautica 143:16-28) is the lever that could either confirm
  rediscovery (close out the region) or refine the anchor (re-open
  lit-fresh class). The paper costs negligibly to acquire
  (Elsevier ScienceDirect article, library-accessible). This is the
  Phase 2 follow-up: it is a literature-acquisition task, not a
  computational task.

* The candidate's V_inf tuple at Titan (1.72-1.75 km/s) is in the
  Cassini Titan-flyby V_inf operating band (~5.8 km/s typical
  Cassini-Titan encounter, with some prime-mission flybys at lower
  V_inf for spectroscopy targeting). The (1.72 km/s) is much LOWER
  than typical Cassini operations, and looks more like a low-energy
  Tisserand-pump-tour-class regime — but this is speculation pending
  the source check.

* Saturn-system other-pair sweep (Part B) confirms the result: zero
  new SILVER, and the only lit-fresh-class pair (Titan-Tethys) fails
  the residual + physical-sanity gates. The Saturn-other-pair
  territory has no admissible candidates that improve on the anchored
  Titan-Rhea-Titan baseline.

* The 3D existence probe (Part C) returns the same negative as #312
  Uranus and #341 Neptune at the L1-Lyapunov-tail seed convention:
  no non-planar convergence at Saturn-Titan (24 planar collapses, 12
  failures), and total failure at Saturn-Rhea (mu too small). The
  inclined-orbit lever does not open at this seed convention; Saturn-
  Titan halo families need bifurcation-continuation from a known
  planar Lyapunov member (Phase 2).

## Pointers

* Scripts:
  * `scripts/scan_344_saturn_titan_rhea_finer.py`
  * `scripts/scan_344_saturn_other_pairs.py`
  * `scripts/scan_344_saturn_3d_probe.py`
* Data:
  * `data/scan_344_saturn_titan_rhea_finer.jsonl`
  * `data/scan_344_saturn_robustness.jsonl`
  * `data/scan_344_saturn_titan_iapetus.jsonl`
  * `data/scan_344_saturn_titan_dione.jsonl`
  * `data/scan_344_saturn_titan_tethys.jsonl`
  * `data/scan_344_saturn_rhea_dione.jsonl`
  * `data/scan_344_saturn_titan_rhea_dione_3body.jsonl`
  * `data/scan_344_saturn_other_pairs_index.jsonl`
  * `data/scan_344_saturn_3d_probe.jsonl`
* Substrate (READ-ONLY): `src/cyclerfinder/core/satellites.py`,
  `src/cyclerfinder/search/discovery_campaign.py`,
  `src/cyclerfinder/search/cr3bp_general_periodic_3d.py`,
  `src/cyclerfinder/search/physical_sanity.py`,
  `src/cyclerfinder/search/literature_check.py`.
* Reference baseline: `scripts/scan_320_epoch_aware_moon_systems.py`
  + `data/scan_320_epoch_aware_saturn.jsonl` (commit 0e6f3f2).
* Reference analogs: `docs/notes/2026-06-17-341-neptune-extended-sweep.md`
  (Neptune Proteus-Triton clean-negative analog),
  `docs/notes/2026-06-17-339-silver-quasi-cycler-admission.md` (Uranus
  Umbriel-Oberon admitted-SILVER analog — body-set-only match at
  ANCHOR_COUNT = 0 admitted, vs. Saturn here at ANCHOR_COUNT = 2 not
  admitted).
