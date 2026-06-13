# Repeated-moon multi-rev genome — build + reproduce-before-search gate (#254)

Build record for the Track-A repeated-moon multi-revolution cycler genome of
`docs/notes/2026-06-14-repeated-moon-multirev-genome-design.md`. Module:
`src/cyclerfinder/search/moon_cycler_genome.py`; tests:
`tests/search/test_moon_cycler_genome.py`.

Scope of this build: design **steps 1-4** (registry/graph, genome
representation, periodicity corrector, and the mandatory reproduce-before-search
gate). The open-ended **search (step 5)** was deliberately NOT run — that long
combinatorial compute belongs in the #253 discovery-campaign daemon, not a
one-shot agent.

## What was built

1. **Moon-system registry + Tisserand/V_inf graph** (`MoonSystem`,
   `jupiter_system`, `moon_vinf_to_tisserand`/`moon_tisserand_to_vinf`,
   `moon_linkable`, `vinf_graph_edges`). Reuses the body-agnostic Tisserand
   machinery (`search/tisserand.py`, #75) centre-aware via `mu = PRIMARIES[
   "Jupiter"]` and the satellite registry (`core/satellites.py`, #76/#117). No
   new physics — a thin moon-scale wrapper. Galilean periods derive from the
   registry SMA + Jupiter GM (Io 1.770 d, Europa 3.552 d, Ganymede 7.155 d,
   Callisto 16.691 d — match JPL SSD to <1e-2 d).

2. **Genome representation + decision vector** (`LegGene`, `EncounterGene`,
   `MoonCyclerGenome`). Moon sequence `[m_1..m_k]` (repeated each cycle),
   per-leg `(p:q resonance, n_rev)`, per-encounter flyby B-plane angle, overall
   epoch + perijove scale. `to_vector`/`from_vector` are a lossless round-trip;
   `is_valid` is the cheap structural filter. The three published Liang members
   (`liang_member_genome`) encode validly: the
   Callisto-Ganymede-Callisto-Europa-Callisto sequence (5 flybys / 4 legs), all
   `n_rev = 1` one-rev Lambert arcs, with the member's perijove scale.

3. **Repeated-sequence periodicity corrector** (`PeriodicityResidual`,
   `liang_periodicity_residual`). ONE canonical residual: the worst V_inf-
   magnitude continuity defect across the cycle's flybys (the ballistic-cycler
   periodicity condition), in km/s, evaluated at the Eq. 16-anchored cumulative
   flyby epochs — the same scalar the daemon would minimise for a new
   candidate. Routes through the trusted same-model reconstruction
   `cge_scaffold.reproduce_member` (#222) so the corrector inherits the audited
   planet-centric Lambert legs rather than duplicating them.

4. **Reproduce-before-search gate** (`reproduce_before_search_gate`,
   `GateResult`). See result below.

## REPRODUCE-BEFORE-SEARCH GATE RESULT — **PASS**

The genome reproduces all three numerically-printed Liang CGE members. EXPECTED
side = Liang et al. 2024 published Tables 3/5/7 (sourced); ACTUAL side = the
genome's same-model reconstruction. Tolerance = the SOURCED print-precision
floor `cge_scaffold.vinf_print_tolerance_kms` (dominated by the 4-decimal
Table 1 mean-motion quantization), never a hand-tuned number.

| catalogue id | member | pass | worst V_inf residual (km/s) | sourced V_inf tol (km/s) | worst periodicity residual (km/s) |
|---|---|---|---|---|---|
| liang-2024-cgcec-111-highperijove | A | True | 1.516e-02 | 1.372e-01 | 8.094e-03 |
| liang-2024-cgcec-110-highperijove | B | True | 1.372e-02 | 1.377e-01 | 7.809e-03 |
| liang-2024-cgcec-111-lowperijove  | C | True | 4.819e-02 | 1.382e-01 | 9.034e-03 |

Every worst per-flyby V_inf residual is 2.8x-9x INSIDE its sourced tolerance,
and every periodicity (V_inf-continuity) residual is ~1e-2 km/s — also inside
tolerance. The members are published ballistic (residual defect Delta-v below
1e-8 m/s, paper p. 13); our reconstruction is input-precision-limited, so the
achievable residual is the print floor (~1e-2 km/s), not 1e-11 — exactly as the
#222 reproduction established.

**Member D** (`liang-2024-cgcec-ephemeris-2033`, V0) is intentionally NOT a gate
golden: its per-flyby data are published as FIGURES ONLY (no printed per-flyby
V_inf/ToF), so it is structurally unreproducible from sourced numbers. The gate
covers exactly the three numerically-printed V1 members.

## Verdict

The genome's validation gate is **settled PASS** — the genome and corrector
faithfully recover Liang's published CGE members from sourced data in a
topology (repeated-moon, multi-rev, Callisto 3x/cycle) the prior zero-rev
single-encounter genome could not represent. This clears the design's
reproduce-before-search precondition: the #253 discovery-campaign daemon is now
warranted to run the open-ended enumerate-and-close search (design step 5) over
the Jovian moon-seq x resonance box, with every SILVER survivor fed to the
unchanged V0-V5 gauntlet.

## Discipline notes

- NO catalogue writeback. The gate validates the genome only; a future closed
  search candidate is SILVER, not a validated cycler.
- Goldens are Liang's PUBLISHED values (Tables 3/5/7), routed through the
  audited #222 same-model reconstruction; the residual is our model's output,
  asserted against the sourced print tolerance — not circular.
- The search (step 5) and daemon-host (step 6) are explicitly out of scope for
  this build (anti-hang: long compute -> #253 daemon).
