# Seed-run result: Cuevas del Valle et al. 2026 Scenario II target → EM L1 southern halo (#190)

**Date:** 2026-06-13. Closes the two adoption actions of
`docs/notes/2026-06-10-cuevas-del-valle-2026-cr3bp-mpc-mining.md` (Sec. 5 items 1–2),
as amended by the L1 relabel in
`docs/notes/2026-06-11-cuevas-del-valle-2023-floquet-mining.md` (Sec. 3.3).

## 1. Jacobi-gap ΔV lower bound — ADOPTED

`cyclerfinder.core.cr3bp.jacobi_gap_dv_min(speed, delta_c)` implements the
Jacobi-gap minimum-ΔV technique of Cuevas del Valle, Urrutxua & Solano-López 2026
(CEAS EuroGNC 2026, CEAS-GNC-2026-012, Sec. 7.2): an impulse at fixed position moves
C only through −v², so ΔC = v₀² − v_f² and the cheapest impulse is the tangential
burn, ΔV_min = |√(v₀² − ΔC) − v₀|. Tight for a single impulse at the given speed;
raises when ΔC > v₀² (zero-velocity ceiling). Tests (hand-derived algebra,
consistency with `jacobi_constant`, sampled-direction lower-bound property, edge
cases) in `tests/core/test_cr3bp.py`. The paper's own 0.2734 m/s Scenario I figure
is NOT reproducible (relative states, untabulated μ) and is not used as a golden.

## 2. Seed run — CONVERGED (clean positive, L1 relabel confirmed)

Seed: the Scenario II target state (paper Sec. 7.3 p. 18, full precision)

```text
[0.824024728136525, 0, -0.054501847320725, 0, 0.164671964079122, 0]
```

run through `correct_periodic` with the physical Earth-Moon μ = 0.0121505844
(`cr3bp_system("Earth", "Moon")`) and period guess 2.7549 nd (the paper's
approximate CHASER period — the only period printed).

| Quantity | Result |
|---|---|
| Converged | yes; closure residual 7.7e-13 |
| Seed correction \|ds\| | 2.3e-3 nd |
| Period T | 2.760206 nd (+0.19% vs the chaser's ~2.7549) |
| Jacobi C | 3.151692 (ours; none published) |
| Radau crosscheck | closed, ΔC drift 1.1e-14 |
| x range | [0.82396, 0.86775] — straddles x_L1 = 0.836915, never near Moon/L2 |
| z range | [−0.05295, +0.04359] — dominant excursion southern |
| Monodromy λ_max | 1.577e3 → ν = ½(λ + 1/λ) ≈ 788.6 (strongly unstable) |

Verdicts:

- **The 2023 note's L1 relabel is dynamically confirmed:** the converged orbit
  encircles L1 (x_L1 = 0.836915) and stays Earth-side of the Moon throughout — an
  Earth-Moon **L1 southern halo**, not L2 as the 2026 paper's text says.
- **Provenance: sourced-seed only.** μ, C and the target's period are not published;
  the converged (state, T, C) tuple above is OUR corrector's output and must not be
  treated as a golden. The published ~2.7549 nd chaser period (adjacent northern
  member) is the only sourced comparison and agrees to 0.19%.
- Pinned by `test_cuevas2026_em_l1_southern_halo_seed_converges` in
  `tests/search/test_cr3bp_periodic.py` (qualitative topology/stability assertions +
  the 1% sourced-period check).
- **No catalogue writeback** (unchanged from the mining notes): no published
  (μ, state, T, C) tuple; nothing meets v4.2 standards.
