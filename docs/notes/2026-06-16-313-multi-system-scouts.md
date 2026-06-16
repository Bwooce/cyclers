# #313 Phase 1: Multi-system tulip + BCR4BP scouts (Mars moons / Sun-Jupiter moons)

Cheap parallel scout using existing genomes at fresh mu values. No new
corrector code, no catalogue writeback. Two scripts, four systems, three
honest verdicts.

## Discipline anchors

- READ-ONLY on `src/cyclerfinder/genome/*.py` and `src/cyclerfinder/core/bcr4bp.py`.
  Both scout scripts wrap existing routines (`find_tulip_at_system`,
  `multi_shoot_periodic`, `correct_bcr4bp_periodic`,
  `continue_bcr4bp_family_in_musun`, `correct_symmetric_free_period`).
- NO catalogue writeback. JSONL deliverables only.
- NO novelty claims. Tulip families at Mars moons and L1 Lyapunov families at
  Jupiter moons are mu-scalings of published Earth-Moon (Koblick 2023) / EM L1
  Lyapunov families (Braik-Ross 2026 / generic Conley-McGehee). Lit-fresh
  hits would be candidates for a future BCR4BP / QP-tori gauntlet, not
  discoveries.
- CR3BP-limit anchor (BCR4BP @ mu_sun=0 → CR3BP) IS the structural-correctness
  test for the Sun-Jupiter-moon BCR4BPSystem constants. Both scouts pass it.

## Part A — Mars-Phobos / Mars-Deimos tulip + multi-rev

Script: `scripts/scan_313_mars_moons.py`. Output:
`data/scan_313_mars_phobos.jsonl` (6 rows), `data/scan_313_mars_deimos.jsonl`
(6 rows).

### Mars-Phobos (mu ≈ 1.65e-8, L = 9375 km)

| Attempt | Outcome | Notes |
|---|---|---|
| `tulip_direct` (Tier A/B at Np=2)        | `seed_no_converge` | Replicates #281 |
| `multishoot_n=2` (Koblick Np=2 paper seed)| `RuntimeError: step size below spacing` | DOP853 chokes on Phobos perilune graze (mu~1e-8, perilune ~ secondary radius) |
| `multishoot_n=3`                         | `RuntimeError: step size below spacing` | same |
| `multishoot_n=4`                         | `RuntimeError: step size below spacing` | same |
| `tulip_paperseed_np3` (Np=3 paper seed)  | `seed_no_converge`     | |
| `tulip_paperseed_np4` (Np=4 paper seed)  | `direct_seed_match` — **but flagged spurious below** | |

### Mars-Deimos (mu ≈ 2.25e-9, L = 23457 km)

| Attempt | Outcome | Notes |
|---|---|---|
| `tulip_direct`                           | `seed_no_converge` | |
| `multishoot_n=2`                         | **converged** — flagged spurious below | |
| `multishoot_n=3`                         | `RuntimeError: step size below spacing` | |
| `multishoot_n=4`                         | `RuntimeError: step size below spacing` | |
| `tulip_paperseed_np3`                    | `seed_no_converge` | |
| `tulip_paperseed_np4`                    | `direct_seed_match` — flagged spurious below | |

### Honest verdict — Mars-moon CR3BP tulip basin

**The three "converged" rows are spurious matches of the Tier A petal_count
gate to planar orbits**, not genuine multi-petal halo-tulip geometry:

- Mars-Phobos `tulip_paperseed_np4`: x0=1.023731 (the paper-row fixed
  crossing), z0=-2e-12 (collapsed to planar), ydot0=-0.0472, T_TU=6.2779
  (~ 2π), C=2.9994, petal_count=2. The corrector flattened the seed onto
  the xy-plane at Mars-mu — find_tulip_at_system's Tier A gates on
  petal_count, which fires "2" on the planar 2:1 Lyapunov-like orbit. NOT a
  halo-shaped tulip (z0 ≈ 0).
- Mars-Deimos `multishoot_n=2`: x0=0.9994, z0=0.0011 (negligible), ydot0=0.0013,
  T_TU=2.421, petal_count=1. Multi-shoot landed an interior orbit near
  Deimos with one perilune passage — also not a tulip.
- Mars-Deimos `tulip_paperseed_np4`: same z0≈0 collapse as Phobos.

Closure residuals were strong (<1e-11), confirming these ARE valid CR3BP
periodic orbits — just not the tulip family the Tier A gate thinks they are.

**Discovery verdict: the Earth-Moon Koblick tulip basin is NOT directly
reachable at Mars-mu without continuation in mu first.** The corrector's
basin of attraction at mu~1e-8 / 1e-9 collapses the perpendicular-crossing
seed to planar orbits, and multi-shooting at n_segments≥3 fails
catastrophically (DOP853 step-size collapse near Phobos perilune — at
mu~1e-8 the secondary's Hill sphere is so small that nearly any planar
incursion grazes the body). A mu-continuation from Earth-Moon down to Mars
through Pluto-Charon (mu~0.108) and Jupiter-Io (mu~5e-5) would be the
disciplined Phase 2 path; this scout's clean negative justifies that
investment if the discovery question is "do Mars-moon tulip orbits exist?".

The hostile-mu result also confirms the #281 single-shot finding: the tulip
genome's basin of attraction is a stronger function of mu than the
Conley-McGehee L1 Lyapunov family (which DOES converge at all these mu —
Part B uses it).

## Part B — Sun-Jupiter-Europa / Sun-Jupiter-Io BCR4BP

Script: `scripts/scan_313_sun_jupiter_moons.py`. Outputs:
`data/scan_313_sun_jupiter_europa.jsonl` (54 rows: constants, seed, anchor,
50 members, summary), `data/scan_313_sun_jupiter_io.jsonl` (same shape).

### BCR4BP constants (sourced + computed)

|  | Sun-Jupiter-Europa | Sun-Jupiter-Io |
|---|---|---|
| mu              | 2.528e-5     | 4.704e-5     |
| L (km)          | 671100       | 421800       |
| TU (s, days)    | 48844, 0.565 | 24338, 0.282 |
| a_sun (in LU)   | 1159.80      | 1845.28      |
| mu_sun          | 1047.57      | 1047.57      |
| omega_sun       | 0.99918      | 0.99959      |
| Sun synodic TU  | 6.288        | 6.286        |

Sources: GM_Sun (JPL DE440 / IAU 2015), GM_Jupiter_sys + moon GMs + moon SMAs
(JPL SSD gm_de440 / phys_par / mean elements), Jupiter SMA (JPL DE440
heliocentric), AU (IAU 2012 exact). All values trace to constants already in
the repo registry (`core/satellites.py`, `core/constants.py`).

### CR3BP-limit anchor (sourced-golden discipline)

| Moon | corr_res | indep_closure | converged |
|---|---|---|---|
| Europa | 2.5e-14 | 4.3e-13 | ✓ |
| Io     | 7.9e-14 | 2.0e-12 | ✓ |

The CR3BP L1 Lyapunov re-converges through the BCR4BP corrector at mu_sun=0
to machine precision — the BCR4BPSystem constants are CORRECT for these
systems. (Per the brief: "the CR3BP-limit anchor MUST pass to numerical
precision, same as #292 Phase 1. If it doesn't, the BCR4BPSystem constants
are wrong.")

### mu_sun continuation (0 → 1047.57)

| Moon | n_members | x0 final | vy final | T (TU) final | Drift in IC | Indep closure final | Stability |
|---|---|---|---|---|---|---|---|
| Europa | 50/50 | 0.976811 | 0.023953 | 3.0853 | < 1e-9 | 2.0e-7 | hyperbolic_pair throughout |
| Io     | 50/50 | 0.972088 | 0.024085 | 3.0403 | < 1e-9 | 6.1e-8 | hyperbolic_pair throughout |

### Honest verdict — Sun-Jupiter-moon L1 substitute

**The L1 Lyapunov substitute family survives the full Sun perturbation
continuation at both moons.** Quantitatively, the survival is much stronger
than at Earth-Moon (the #303 EM L1 Lyapunov continuation moved x0 ~1e-4 over
the same mu_sun range, here Δx0 < 1e-9). The reason is geometric: a_sun_LU
~1100-1850 vs Andreu's 388.8 for EM. The Sun's direct acceleration per test
particle scales as mu_sun / a_sun² (and the indirect term as mu_sun / a_sun³),
so at Jupiter's distance the Sun perturbation is ~(388.8/1160)² ≈ 11x weaker
than at Earth despite mu_sun being only ~3x smaller in mass-ratio units (Sun
~1047 jovian masses vs ~328900 EM-masses).

This is a **clean qualitative negative for "Sun perturbation matters at
Jupiter for L1-substitute orbits"** — at corrector precision, Sun-Jupiter-
moon L1 Lyapunov geometry is indistinguishable from the unperturbed
CR3BP L1 Lyapunov geometry. The same conclusion does NOT necessarily hold for
larger-amplitude libration members (this scout used C = C_L1 - 5e-4, a
small-amplitude Lyapunov); a follow-on study at C = C_L1 - 0.05 (deeper
into the L1 substitute family) might reveal the perturbation re-asserting.

## Surprises

- **Mars-Deimos `multishoot_n=2`** converged where Mars-Phobos n=2 didn't.
  At mu~2.25e-9 (smaller than Phobos at 1.65e-8) the perilune graze is even
  more violent, but the integrator landed an interior orbit AT the secondary
  rather than a near-secondary tulip. Expected easier-to-fail; converged
  instead — but to a planar 1-petal orbit, not a tulip. Spurious convergence,
  not a discovery.
- **Mars-Phobos `tulip_paperseed_np4` Tier A "direct_seed_match"**: the
  corrector flattened z0 from 0.138427 (Koblick Np=4 paper seed) to ~0 at
  mu=1.65e-8. The `petal_count` gate fired "2" on the resulting planar
  orbit and returned `success=True`. **The `find_tulip_at_system` Tier A
  gate has a known false-positive mode at very small mu** — the planar
  L4/L5 / 2:1 resonance basin can absorb the perpendicular-crossing seed
  while still satisfying petal_count==Np_target.
- **Sun-Jupiter-Europa and -Io continuation drift is below 1e-9 in
  (x0, vy)** across 50 mu_sun steps from 0 to 1047.57. The Sun perturbation
  at Jupiter is so geometrically weak that L1 Lyapunov orbits are
  effectively unperturbed.
- **Easiest expected, failed**: Mars-Phobos multi-shoot at n=2 (mu only one
  order of magnitude below Pluto-Charon's 1e-1 where things worked) — RHS
  step-size collapse before any progress, every n_segments.
- **Hardest expected, succeeded**: Sun-Jupiter-Europa BCR4BP CR3BP-limit
  anchor at completely fresh constants. The BCR4BP corrector hit
  corrector_residual=2.5e-14 on the first try. The genome is robust under
  parameter swap.

## Literature-fresh count

Zero literature-fresh hits across all four systems:

- Mars-Phobos / Mars-Deimos: no genuine tulip orbits converged (the
  "successes" are spurious Tier A matches to planar orbits, not tulip family
  members).
- Sun-Jupiter-Europa / -Io: the L1 Lyapunov family is generic
  Conley-McGehee, published exhaustively; survival under Sun perturbation
  in BCR4BP is qualitatively expected (Simo-Jorba-Gomez methodology), even
  if the quantitative drift size here is a new datapoint.

The Phase 1 scout deliverable is **two clean negatives + structural
correctness confirmation of the BCR4BPSystem at fresh mu**, NOT a
discovery.

## Phase 2 path (if pursued)

If a future task wants to push these scouts:

1. **Mars-moon tulip via mu-continuation**: continue the tulip family from
   Earth-Moon (mu=1.215e-2) down to Pluto-Charon (1.08e-1)... wait — that's
   UP. Pluto-Charon is the largest mu in the registry. The continuation needs
   to go EM (1.2e-2) → Jupiter-Galileans (2.5e-5 to 5.7e-5) → Mars-Phobos
   (1.65e-8). Several orders of magnitude in mu; needs intermediate
   anchor systems and probably switching to a SUNDMAN-REGULARIZED multi-shoot
   STM (the current `multi_shoot_periodic` uses regularized states + plain
   STM, which is what fails at mu~1e-8).
2. **Sun-Jupiter-moon L1 deeper into the family**: probe at C = C_L1 - 0.05
   (large-amplitude Lyapunov) and the L1 halo family (3D out-of-plane) to
   check whether the qualitative survival result holds where Sun-perturbation
   geometry is potentially first-order. This would also let `correct_bcr4bp_periodic`
   exercise the `FREE_VARS_HALO` mask (#304) at fresh-mu.

Neither is required for the Phase 1 deliverable — both clean negatives
stand.

## Files

- `scripts/scan_313_mars_moons.py` — Part A driver
- `scripts/scan_313_sun_jupiter_moons.py` — Part B driver
- `data/scan_313_mars_phobos.jsonl` (6 rows)
- `data/scan_313_mars_deimos.jsonl` (6 rows)
- `data/scan_313_sun_jupiter_europa.jsonl` (54 rows)
- `data/scan_313_sun_jupiter_io.jsonl` (54 rows)
- This note

## Test gate

`uv run pytest tests/genome/test_tulip.py tests/genome/test_bcr4bp_genome.py
tests/genome/test_bcr4bp_continuation.py -x --timeout=180` — 22 passed
(green pre- and post-scout).
