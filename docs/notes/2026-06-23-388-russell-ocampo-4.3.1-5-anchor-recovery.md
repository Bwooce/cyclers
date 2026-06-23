# #388 — russell-ocampo-4.3.1-5 multi-arc closure trial: ANCHOR RECOVERED (the wall is energy-selective)

**2026-06-23.** The one #365 row that was descriptor-available-but-not-ingested
(`negative_results.yaml` resweep gate). Its per-arc descriptor is published in
McConaghy-Russell-Longuski 2005 Table 2 (in-corpus, digested
`2026-06-17-digest-mcconaghy-2005.md`) but `free_return_arcs` was `null`. Ingested it
and ran the multi-arc DSM closure lane — the last genuinely-un-run direct-lane trial.

## Descriptor ingestion (sourced, validated)

McConaghy-2005 Table 2 formal label: `4 g(7-1/14, 5-1/14 rev, L) f(1:1, 84.039°, ∓90.0°)
h(0.5, 0, U, ±5.961°)` — n=4 synodic, three legs (generic / full-rev 1:1 resonant /
half-rev). Ingested as `free_return_arcs`:

| arc | type | tof_years | resonance | raw_descriptor |
|---|---|---|---|---|
| 1 | generic | 7.071429 (=7-1/14) | — | `g(7-1/14, 5-1/14 rev, L)` |
| 2 | full-rev | null (M:N-determined) | '1:1' | `f(1:1, 84.039, ∓90.0)` |
| 3 | half-rev | 0.5 | — | `h(0.5, 0, U, ±5.961)` |

**Repeat-constraint check (Eq. 28, McConaghy-2005): Σt_f = 7.0714 + 1 + 0.5 = 8.5714 =
4 × (2-1/7) = 4S — EXACT.** This confirms the mixed-number reading ("7-1/14" = 7+1/14,
not 7−1/14) and the leg t_f assignments. (Caught a YAML sexagesimal trap on the way:
unquoted `1:1` parses to 61 in YAML 1.1 — resonance is force-quoted `'1:1'`.)

## Closure trial result

`close_multiarc_row(row, Ephemeris("astropy"), n_starts=25, gradient="lambert",
tol_kms=0.1)`, sequence E-E-M-M:

- **best residual = 0.4852 km/s** — a TRUE FLOOR (identical at max_nfev 400 and 900;
  not an optimiser plateau). NOT converged to the 0.1 km/s ballistic gate.
- **total ΔV ≈ 0.995 km/s.**
- **Emerged V∞ recovers the sourced anchor:** Earth out ≈ 3.04 (sourced Russell 2004
  Table 3.4: **3.1**), Mars ≈ 2.24–2.58 (sourced **2.5**). Within ~0.05–0.08 km/s.

## Why this matters — the #388 wall is energy-selective, not universal

Every prior direct-lane row collapsed OFF-anchor: S1L1 (anchor 4.7/5.0) and the V3
regressions all emerged at ~half their anchor V∞ (the low-energy basin). **4.3.1-5 does
NOT collapse — it recovers its published anchor.** The distinguishing variable is energy:
4.3.1-5 is "one of the lowest-V∞ E-M cyclers known" (Russell §3.8, near-Hohmann,
AR=0.99), so its genuine anchor *is* in the low-energy basin the lane gravitates toward.
The wall is therefore **family/energy-selective**: high-V∞ cyclers drift off-anchor;
low-energy near-Hohmann ones recover. This is a new, sharper characterization than
"universal seed-basin wall".

## Promotion status — NOT a clean V0→V1 (yet)

The trial recovers the anchor but does **not** ballistically close (residual 0.485 >
0.1 gate; total ΔV ~0.995 km/s). This is consistent with Russell's own classification
of 4.3.1-5 as **near-ballistic — "requires a small powered nudge"** (AR=0.99, aphelion
just short of Mars's circular orbit). So the ~0.485 km/s floor is plausibly the genuine
near-ballistic DSM, not pure non-convergence. Either way it is not a ballistic V0→V1.

**Open next step (B3):** Russell's note says it is "genuinely ballistic crossing Mars's
orbit at perihelion" in the real-eph (eccentric-Mars) model — the lane's structured-epoch
search may not have hit a Mars-perihelion-favorable epoch. A bounded epoch-targeted
re-run (Mars at perihelion) tests whether it closes below the gate. If it does → genuine
reproduction (mandatory independent cross-check before any promotion, per orbit-closure
discipline). If it floors at ~0.485 again → record as characterized near-ballistic
anchor-recovery (no promotion), and the residual is the real DSM.

## Disposition

- Descriptor ingested into `catalogue.yaml` (clean sourced data win regardless of closure).
- `negative_results.yaml` 4.3.1-5 entry updated: now RUN; outcome = anchor-recovered
  near-ballistic (0.485 km/s floor), distinct from the off-anchor collapse; promotion
  pending the Mars-perihelion epoch trial + cross-check.
