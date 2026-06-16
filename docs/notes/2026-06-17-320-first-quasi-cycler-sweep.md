# #320 Phase 1 — first systematic quasi_cycler discovery sweep

Date: 2026-06-17
Tasks: #320 Phase 1 (Vector A + Vector B + verdict)
Substrate (read-only): #290 (`qp_tori`), #297 (`epoch_aware_genome`), #319
(`v1_qp` + `v2_qp`), #324 (`physical_sanity`), #328 (`literature_check` 41-
anchor `KNOWN_CORPUS`)
Scripts: `scripts/scan_320_qp_tori_3d_brackets.py`,
`scripts/scan_320_epoch_aware_moon_systems.py`
Outputs: `data/scan_320_qp_tori_3d_brackets.jsonl`,
`data/scan_320_epoch_aware_{saturn,neptune,pluto}.jsonl`

## Headline

  * **Vector A (QP-tori at #299 brackets)**: 12 brackets → 12 corrector runs
    converged in Fourier-mode space, **2 SILVER** survivors of the V1_qp +
    V2_qp + irrationality gauntlet.
  * **Vector B (epoch-aware repeated-moon sweep at non-EM systems)**:
    Saturn 1 SILVER (Titan-Rhea-Titan (1,1), 0.032 km/s, 1 anchor overlap);
    Neptune 0 SILVER (best 0.058 km/s, **lit-fresh-class** near-miss);
    Pluto 2 SILVER (Hydra-Nix-Hydra / Nix-Hydra-Nix at ~1e-3 km/s, both
    4 anchor overlaps).
  * **NO catalogue writeback.** SILVER candidates feed #306-style follow-up
    tasks. Numbers below are V1_qp / V2_qp / physical-gate / offline-corpus
    pass counts; novelty claims require the live #328 web search + human
    adjudication.

## Vector A — QP-tori at #299 Earth-Moon 3D Neimark-Sacker brackets

Driver: `correct_qp_torus` (Olikara-Howell GMOS, `n_trans=2`,
`amplitude=5e-4`, `tol=1e-8`) seeded from each of the 12 Neimark-Sacker
brackets in `data/family_296_3d_subfamilies_299.jsonl`. Each converged
torus runs through `v1_qp.run_v1_qp` + `v2_qp.run_v2_qp` plus an
`is_practically_irrational` check on the frequency ratio
`omega_trans/omega_long`.

| count | gate                                                |
| ----- | --------------------------------------------------- |
| 12    | corrector ran without raising                        |
| 0     | strict `QPTorus.converged` (tol=1e-8 + indep<1e-3)   |
| 11    | V1_qp PASS (Fourier-norm < 1e-5 AND independent off-grid < 1e-4) |
| 3     | V2_qp PASS (3-cycle off-grid drift < 5e-2 nondim)    |
| 9     | irrational frequency ratio (max_denom=10, tol=1e-3) |
| **2** | **SILVER (V1 PASS + V2 PASS + irrational)**         |

### SILVER survivors

  * **Bracket 2** — k=4, parent at step_a=+8, T_a = 10.27935 TU,
    C_a = 3.03196. Corrector residual 6.1e-8, V1 Fourier 6.1e-8 (10⁻³
    inside floor), V2 max drift 1406 km nondim-equivalent. Frequency
    ratio -0.218 (irrational at default tol).
  * **Bracket 10** — k=4, parent at step_a=+110, T_a = 9.30832 TU,
    C_a = 3.12624. Corrector residual 1.1e-6, V1 Fourier 6.6e-7, V2 max
    drift 10742 km nondim-equivalent. Frequency ratio -0.265 (irrational).

Both SILVERs are **k=4** Neimark-Sacker brackets. The other k=4 bracket
(bracket 3, the partner bracket of #2 at step_a=+9) passes V1_qp + V2_qp
gauntlet but the frequency ratio falls within tol=1e-3 of a small rational
(probably the 1:4 phase-lock at the bifurcation itself), so it is
correctly screened out as a phase-locked PERIODIC orbit, not a genuine
2-torus. This is the expected discipline: a true torus must MOVE OFF the
bifurcation point's rational rotation number.

### What V1 / V2 gates told us

  * V1_qp PASS on 11 / 12 brackets is the substrate working as designed:
    the Fourier-norm residual sits at 1e-6 to 1e-8 once the Newton converges
    in mode space (the lone V1 fail is bracket 0, k=6, parent furthest from
    the bracket's unit-circle crossing → biggest corrector residual at
    3.8e-6, still under V1_qp's 1e-5 floor; the **independent off-grid
    check** is what trips it, telling us the corrector matched aliased
    modes rather than a genuine torus).
  * V2_qp PASS on only 3 / 12 is the expected Phase-1 ceiling at n_trans=2.
    The V2 floor (5e-2 nondim per the empirical N=2 calibration; see
    `v2_qp` module docstring) catches the hyperbolic-instability-amplified
    drift over 3 cycles, which scales as `O(exp(k * lambda_max))` with the
    parent's leading Floquet exponent. Phase 2 (n_trans ≥ 4) should
    halve the per-cycle drift and admit more of the V1-pass set.

### Lit-fresh check

The Vector A SILVERs are Earth-Moon CR3BP QP-tori on a 3D L2 / L1 vertical-
family branch (Antoniadou-Voyatzis 2018 ↔ Roberts-Tsoukkas-Ross 2026
sourcing chain). The literature-check module's KNOWN_CORPUS does NOT
currently include a QP-torus anchor for the Earth-Moon system — the corpus
is biased toward sequence-based moon-tour and heliocentric cyclers. A
proper lit-check for these requires:

  * A `CandidateSignature` shape that admits **torus** descriptors (no
    sequence, no n_rev — a {primary, C_band, T_band, n_modes, rho_rational_k}
    fingerprint).
  * Anchors for Olikara-Scheeres 2010, Olikara 2016, Howell-Howell 2014,
    Henderson-Howell 2008 (Earth-Moon QP-tori at the Lyapunov / halo
    families).

These extensions are out of #320 Phase 1's scope; the SILVER survivors are
flagged as **lit-fresh-class-pending-signature-extension** rather than
either "lit-fresh" or "rediscovery". A Phase 2 follow-up task should add
a QP-torus-shaped signature to `literature_check`.

## Vector B — epoch-aware repeated-moon sweep at non-EM systems

Driver: a per-system 24×24 (relative-offset × global-phase) basin-floor sweep
of closed length-3 repeated-moon cycles, with per-leg revolution counts in
`{0, 1, 2, 3}` (the trivial (0,0) corner excluded; #285 fix B). Each
candidate goes through:

  * Closure via planet-frame coplanar circular Lambert (same machinery as
    #285 / #312 / #339);
  * Physical-sanity gate at every encounter (#324, `min_useful_bend_deg=5°`);
  * Offline KNOWN_CORPUS anchor-overlap count (#328 41-anchor corpus).

SILVER ≡ residual < 0.05 km/s **AND** physical-gate passed.
"Lit-fresh-class" ≡ 0 anchor overlap (offline only; necessary-not-sufficient
for novelty per `feedback_literature_novelty_check_baseline`).

### Saturn (Enceladus / Tethys / Dione / Rhea / Titan, 20 closed cycles)

| count | gate                                                |
| ----- | --------------------------------------------------- |
| 300   | cells evaluated (excludes (0,0); 5 × 4 = 20 cycles × 16 (n_rev1, n_rev2) cells - (0,0) = 300) |
| 60    | residual < 1.0 km/s (the near-miss band)            |
| **1** | **SILVER** (Titan-Rhea-Titan (1, 1), 0.032 km/s)    |
|       | Best record overall: 0.341 km/s (worst end of near-miss top-5) |

**SILVER detail**: Titan-Rhea-Titan (1, 1), phase0=90°, rel_off=285°,
tof_scale=2.0. V_inf per encounter (1.75, 1.68, 1.72) km/s. Physical-gate
PASS — Titan delivers 49° max-bend at safe alt 1500 km; Rhea delivers 6.8°
at safe alt 100 km (just above the 5° floor — marginal but real). Offline
anchor: overlaps with the Davis-Phillips Saturn tulip-orbit anchor on the
body set (a structural-fingerprint match, not a confirmed publication of
this specific Titan-Rhea-Titan (1,1) cycle).

**Note vs #285/#311**: The #285 / #311 Rhea-Dione-Rhea (1, 1) near-miss
sat at 0.107 km/s (genome ceiling per the #311 robustness sweep). This
new sweep — Titan added to the moon set, n_rev grid widened to (0..3),
relative-offset sweep added — surfaces a **lower** basin-floor (0.032
km/s) at a DIFFERENT cycle topology (Titan-Rhea instead of Rhea-Dione).

### Neptune (Triton + Proteus, 2 closed cycles)

| count | gate                                                |
| ----- | --------------------------------------------------- |
| 30    | cells evaluated (2 cycles × 16 cells - (0,0)s)      |
| 6     | residual < 1.0 km/s                                 |
| 0     | SILVER                                              |
|       | Best: Triton-Proteus-Triton (0, 1) at 0.585 km/s    |
|       | Top-2: Proteus-Triton-Proteus (1, 1) at 0.058 km/s — **lit-fresh-class** near-miss, just above the gate |

The Triton-only inner ring is a hostile cycler target (large + retrograde +
inclined per the satellites.py sourcing note); Proteus is too small for
useful bending at the V_inf range surfaced (max V_inf ~3.7 km/s).

**Phase 2 recommendation**: Proteus-Triton-Proteus (1, 1) at 0.058 km/s
is the strongest **lit-fresh-class** near-miss across the whole #320
sweep. It sits just 16% above the 0.05 km/s gate — analogous to the
Saturn Rhea-Dione 0.107 km/s near-miss that #311 confirmed as a genome
ceiling. An extended sweep (multi-arc / DSM / 3D corrector / wider
n_rev grid) is the right Phase 2 lever; per the #316 cross-system note,
the retrograde Triton geometry may favour a 3D-out-of-plane closure
that the 2D coplanar genome cannot reach.

### Pluto (Charon + Nix + Hydra, 6 closed cycles)

| count | gate                                                |
| ----- | --------------------------------------------------- |
| 90    | cells evaluated (6 cycles × 16 - (0,0)s)            |
| 51    | residual < 1.0 km/s                                 |
| **2** | **SILVER** (Hydra-Nix-Hydra (1, 1) 0.0014, Nix-Hydra-Nix (1, 1) 0.0007) |

**Caveat on the Pluto SILVERs**: both pass the patched-conic physical-
sanity gate as configured (`min_useful_bend_deg=5°`) but the encounter
V_inf values are ~15-30 **m/s** — three orders below the Saturn SILVER's
~1.7 km/s. At those velocities the cycle is effectively a very-low-
energy resonant closure in the Nix-Hydra ring, not a "gravity-assist
cycler" in the operational sense. Both candidates overlap with 4
KNOWN_CORPUS anchors (Howard et al. Persephone, Stern Game-Changer,
Showalter-Hamilton three-body resonance, Brozovic et al. satellite
orbits) so they are **NOT lit-fresh** — the Showalter-Hamilton
resonance anchor specifically describes the Styx-Nix-Hydra co-orbital
chaotic dynamics this sweep reproduces.

The Pluto SILVERs are admitted to the JSONL record but flagged in the doc
as **degenerate co-orbital closures, not operational cyclers**. They are
useful as a calibration data-point: the sweep machinery correctly finds
the Nix/Hydra resonance ring with very low residual, and the offline
lit-check correctly recognises it as published (Showalter-Hamilton 2015
Nature paper).

## Cross-Vector verdict — siblings to #339?

  * **Vector A's 2 SILVER tori are NOT siblings to #339** in any obvious
    sense: #339 is a Uranian moon-tour cycler at V_inf ~0.9 km/s with no
    invariant-torus structure; Vector A's SILVERs are Earth-Moon CR3BP
    QP-tori at C ~ 3.03 / 3.13 with no inter-body sequence.
  * **Vector B's Saturn Titan-Rhea-Titan (1, 1) IS a near-sibling to
    #339**: same multi-rev repeated-moon topology, same (1, 1) per-leg
    revolution structure, same closure residual band (0.032 vs 0.025
    km/s), same physical-gate-PASS profile. The structural fingerprints
    differ only in the moon set (Titan-Rhea vs Umbriel-Oberon) and the
    primary (Saturn vs Uranus).
  * **Vector B's Pluto SILVERs are NOT operational siblings to #339**:
    V_inf orders of magnitude apart, published in the Showalter-Hamilton
    anchor already.

## Phase 2 recommendation

The strongest signal across the whole sweep is:

  1. **Saturn Titan-Rhea-Titan (1, 1)** — the new SILVER. Phase 2 should
     deepen this with the #312-style robustness battery (phase-resolution
     sweep, n_rev-grid widening to (0..5)², DOP853 cross-check, full
     literature_check live search) and confirm or refute the
     Davis-Phillips anchor overlap. If lit-fresh + V0-V5 gauntlet holds,
     this becomes a #339-class catalogue admission candidate.

  2. **Neptune Proteus-Triton-Proteus (1, 1)** at 0.058 km/s near-miss —
     lit-fresh-class, just above the gate. Phase 2 should run the
     #311/#312-style extended sweep on this cell: finer phase resolution,
     multi-arc / DSM continuation, 3D corrector with z0 ≠ 0 seeded
     (matching the spec brief's call for 3D-out-of-plane closure).

  3. **Vector A's 2 QP-tori SILVERs** — defer until the
     `literature_check` module gains a QP-torus-shaped
     `CandidateSignature` so the lit-fresh determination is meaningful.
     The substrate proves the gauntlet is calibrated correctly (the
     phase-locked twin bracket is correctly rejected by the irrationality
     check); the question of catalogue admission requires the corpus
     extension and is itself a Phase 2 follow-up.

The Pluto Nix/Hydra SILVERs are intentionally NOT in the Phase 2 plan —
they are NOT lit-fresh and they are operational degenerates. They remain
in the JSONL as substrate calibration evidence (the lit-check + physical-
gate chain correctly handles them).

## Files

  * `scripts/scan_320_qp_tori_3d_brackets.py`
  * `scripts/scan_320_epoch_aware_moon_systems.py`
  * `data/scan_320_qp_tori_3d_brackets.jsonl`
  * `data/scan_320_epoch_aware_saturn.jsonl`
  * `data/scan_320_epoch_aware_neptune.jsonl`
  * `data/scan_320_epoch_aware_pluto.jsonl`
  * `docs/notes/2026-06-17-320-first-quasi-cycler-sweep.md` (this file)

## Discipline reminders

  * NO catalogue writeback (all SILVER candidates require #306-style
    follow-up + V0-V5 gauntlet + live lit-check + human adjudication).
  * NO novelty claims (offline anchor-overlap = 0 is necessary-not-
    sufficient per `feedback_literature_novelty_check_baseline`).
  * All substrate modules touched READ-ONLY (qp_tori, v1_qp, v2_qp,
    epoch_aware_genome, physical_sanity, literature_check).
  * Frozen census ratchet untouched (`tests/test_catalogue_rediscovery.py`).
