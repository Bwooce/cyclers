# Moon-tour Tier-1 complete (task #76)

**Date:** 2026-06-08
**Plan:** `docs/superpowers/plans/2026-06-06-moontour-tier1-patched-conic.md`
**Design:** `docs/superpowers/specs/2026-06-06-moontour-planetcentric-design.md`

Tier-1 of the moon-tour milestone is shipped. It makes the catalogue's
patched-conic **moon-system** cycler rows computable on the *same* dynamical
model the heliocentric catalogue already uses (Kepler conics + impulsive
gravity-assist V∞ rematch), by adding a **central body = a planet, flyby bodies =
its moons** axis — without disturbing a single heliocentric byte.

## What Tier-1 delivered

| Layer | Module | What landed |
|---|---|---|
| Registry | `core/satellites.py` | `SatelliteData` + `SATELLITES`/`PRIMARIES`; Galilean four + Saturnian midsize + Titan; JPL-SSD-sourced; mean motion derived at import via Kepler III. Full-moon-name code scheme (reconciled against the catalogue's reserved names, no 2-letter aliases). |
| Ephemeris | `core/ephemeris.py` | `_CentredCircularBackend` + `Ephemeris(center=...)` — moon states about the primary, km-scaled. Heliocentric backends byte-identical. |
| Corrector | `search/correct.py` | Two Sun-couplings lifted: `mu_central` plumbed into Lambert; `_max_bend_deg` resolves moon codes via `SATELLITES`. Sun-default keeps the heliocentric solver byte-identical. |
| Tisserand | `search/tisserand.py` | `_a_p_km` + `mu=` resolve a moon; `T = 3 − v∞²` round-trips about Jupiter; `linkable` prunes a Jovicentric pair. |
| VILM | `search/vilm.py` | n:m_K± taxonomy; Eq.(9) V̄∞-efficiency root; Eq.(13) ΔV-min quadrature (no-GA + GA-routed); Europa 3-VILM endgame scalar; admissible ΔV-floor for search pruning. |
| Matching | catalogue/matcher | `(model_assumption, primary)` pool pre-filter — a Jovicentric V∞ never compares to a heliocentric one. |
| Gauntlet | `verify/fidelity.py` + tests | Axis-B persistence about a primary; Axis-A VILM-vs-Lambert agreement; Axis-D wrong-central-μ falsification guard. |
| Catalogue | `data/catalogue.yaml` | Two Jovian rows re-tagged `non-keplerian` → `multi-arc` (gauntlet routes them to invariants, not CR3BP); the Saturnian row keeps `non-keplerian` with an honest Titan-Tier-1 / midsize-Tier-2 split note. |

## Golden discipline

Every EXPECTED side traces to a published source, never to a value our own code
computed:

- **Registry construction** golden = published Endgame Part-1 Table 3 (ã_M /
  Ṽ_M), sourced independently of the JPL-SSD values that built the registry.
- **VILM ΔV / V̄∞** goldens = published Part-1 Tables 1–3 + the worked Europa A6
  scalar (154 m/s / 46 d). The two flagged suspect Part-2 Table 1 cells are
  EXCLUDED (never goldens).
- The **corrector's V∞ is NON-GOLDEN** — the Russell-Strange / Hernandez rows are
  family-seed null-numeric records with no sourced Jovicentric V∞ multiset.
  Tier-1 makes those rows *computable* and *VILM-feasibility-gated*, not
  *V∞-rediscovery-gated*. This is the honest boundary.

## Deviations recorded (honest-risk)

1. **Phase-3 I-E-G closure is bend-INFEASIBLE in the coplanar-circular model.**
   The chain closes about Jupiter (V∞-continuity + periodicity residual → 0,
   converged — proving the corrector is centre-correct), but at the closed
   geometry V∞ lands ~9–12 km/s, forcing 100–150° required turns while the small
   moons bend only 2–5°. Recorded as a strict xfail; tol was NOT loosened and no
   seed was fabricated. Bend-feasible Jovian tours need the VILM layer +
   Laplace-resonance phasing.
2. **VILM admissible ΔV-floor = escape+capture, not the no-GA quadrature.** The
   plan's "no-GA quadrature is a valid lower bound" is backwards — a gravity
   assist *reduces* ΔV, so the no-GA value is not ≤ the with-GA value. The
   implemented floor is the escape+capture energy bound, which is genuinely
   admissible (≤ any real tour-leg ΔV).
3. **Axis-A reframed to VILM-vs-Lambert** on a Jovicentric Hohmann V∞ (two
   independent code paths) rather than VILM-vs-corrector, because the corrector's
   closed family is the higher-V∞ non-bend-feasible one (deviation 1).

## Census

The catalogue is **268 rows** including the moon rows (3 Earth-primary, 2
Jupiter, 1 Saturn). The `MULTI_ENCOUNTER_SEQUENCE` ratchet moved **192 → 223**
this session — but from the Russell 2004 Table 3.4 (+16) and Rall 1970 MIT TE-34
(+15) ingests, **not** from the moon-tour re-tag (which re-classed three existing
rows within their buckets and added no rows).

## Verification gates (all green 2026-06-08)

- `uv run pytest -m "not slow"` — exit 0 (full suite).
- Moon-tour fast tests — 80 passed.
- Moon-tour slow tests — I-E-G closure PASS; Axis-A agreement PASS; Axis-D
  falsification PASS; bend-feasibility strict-xfail (honest-risk).
- `ruff check` + `ruff format --check` clean; `mypy --strict src tests` clean on
  all moon-tour files (one residual tree error is a concurrent free_return WIP
  file, out of this lane).

## What Tier-2 holds (OUT of scope here)

- **CR3BP modelling** — Earth-Moon Arenstorf/Genova/Wittal rows + the Saturnian
  midsize-moon (Mimas/Enceladus/Tethys) members, currently citation-only. Needs
  a CR3BP propagator, Jacobi constant, and the `orbit_elements.cr3bp{}` backfill.
- **The T > 3 ballistic-transfer region** (Endgame Part-2) — needs CR3BP
  reachability (Tier-1.5/Tier-2).
- **JPL three-body catalog pull** — accepted for the later milestone, not this
  one.
- **Bend-feasible Jovian moon-tour seeds** — feeding the VILM layer's
  Laplace-resonance phasing into the corrector so the I-E-G chain closes
  bend-feasibly (flips deviation 1's xfail to a strict pass).
- **Individual Titan cyclers** — split out of the family-seed Saturnian row as
  sourced members are ingested (each takes `cycler_class: multi-arc` /
  `model_assumption: circular-coplanar`).
