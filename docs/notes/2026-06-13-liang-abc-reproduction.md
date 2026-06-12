# Liang et al. 2024 CGE Triple Cyclers — Members A/B/C Idealized-Model Reproduction (#222)

**Source:** Liang, Yang, Li, Bai & Qin, "Callisto-Ganymede-Europa Triple Cyclers,"
*Journal of Guidance, Control, and Dynamics* (Engineering Note), 2024,
DOI 10.2514/1.G008387 — Tables 1–7.
**Mining note:** `docs/notes/2026-06-11-liang-2024-cge-triple-cyclers-mining.md`
(incl. the 2026-06-12 transcription rescan that re-verified every number used here).
**Code:** `src/cyclerfinder/search/cge_scaffold.py`;
golden tests `tests/search/test_liang_cge_reproduction.py` (commit 70a1cc5).

**VERDICT: all three idealized members (A, B, C) REPRODUCE in the paper's own
model.** The printed (phases, leg ToFs, V∞) tuples are self-consistent in the
circular-coplanar two-body patched-conic model to within the precision the
printed inputs permit. Worst |V∞ − printed| over every inbound/outbound side of
every flyby: **A 1.52e-2, B 1.38e-2, C 4.82e-2 km/s** — each 2.8–10× inside the
per-flyby print-quantization tolerance, and 20–200× below the ~1 km/s
wrong-branch/wrong-anchor separation. V0 → V1 promotion evidence assembled for
all three rows (writeback HELD, see end).

## 1. Reproduction route (same-model golden convention)

THEIR model exactly: circular, coplanar moon orbits about Jupiter at the
Table 1 mean motions (Europa 1.7693, Ganymede 0.8782, Callisto 0.3765 rad/day);
conic Jupiter-two-body arcs; instantaneous zero-radius flybys; ballistic. No
real ephemeris anywhere. Orbit radii DERIVED from μ_J + mean motion (the paper
never prints μ_J or the radii); μ_J = the registry Jupiter GM 1.26686534e8
km³/s² (`cyclerfinder.core.satellites.PRIMARIES`). Derived r_Europa − 10000 km
= **660,993.5 km** vs the paper's printed "about 660988 km" (p. 14) — 5.5 km
off, inside both the "about" and the μ_J ambiguity (r ∝ μ^(1/3): ±45 km over
the plausible GM spread).

Route: place the moons by mean motion + Table 2/4/6 phases at the flyby
epochs; solve the Jupiter-frame Lambert problem between consecutive flyby
positions at the printed per-leg ToFs (revolution counts 0–4 enumerated, both
multi-rev branches, prograde); select per leg the solution closest to the
printed V∞ pair (identification only — residuals are then asserted, so a wrong
selection cannot pass); check V∞-magnitude continuity at every flyby and every
printed V∞ on both sides.

### Geometry conventions established (and numerically validated)

1. **Phases at perijove departure.** The Tables 2/4/6 phases are the moons'
   polar angles at t = 0, the perijove departure of the initial conic (period
   T = T_Callisto, perijove on the perifocal +x axis). Since T = T_C, the conic
   semi-major axis equals Callisto's orbit radius, so it crosses Callisto's
   circle exactly at eccentric anomaly E = ±π/2; structure flag 1 (Fig. 2)
   selects E = +π/2, giving t_c0 = (π/2 − e)/n_C. Self-check: Callisto's
   printed phase puts it AT the crossing at t_c0 to **6.5e-6 rad (A/B)** /
   **5.3e-5 rad (C)** — i.e. the printed phases themselves round-trip.
2. **KEY FINDING — the table origin is t_c1 = T + t_c0, not t_c0.** Per the
   Eq. 16 epoch law with n_cycle = 1, the first tabulated Callisto flyby
   (transit time 0 in Tables 3/5/7) occurs ONE FULL SPACECRAFT REVOLUTION
   after perijove departure. The encounter state is identical at t_c0 and
   t_c0 + T (spacecraft and Callisto are both T-periodic), but Ganymede's and
   Europa's downstream phases are not: anchoring at t_c0 instead fails the
   reproduction by ~1 km/s at every later flyby (tested both ways; the
   Eq. 16 guess epochs also rule the t_c0 reading out analytically:
   t_g1 − t_c1 = 4·S_C,G + t_g0 − (T + t_c0) ≈ 31.91 d reproduces the printed
   31.8973 d only with the t_c1 anchor, using the corrected S_C,E erratum
   reading for the Europa term).
3. **All legs are 1-rev Lambert solutions** (n_revs = 1, branch pattern
   high/low/high/low), matching the paper's "all transfer arcs are forced to
   more than one revolution" (p. 9). Selection margins (cost gap to the
   runner-up Lambert solution): **1.00–3.43 km/s** — identification is
   unambiguous by ≥ 7× the largest tolerance.

In-print defects honoured, not imported (`data/errata.yaml`): Eq. 16's
S_G,E → S_C,E subscript (`liang-2024-eq16-synodic-subscript`; only the t_c1
term of Eq. 16 is load-bearing here, but the corrected reading is what makes
the epoch-law cross-check in (2) close); Eq. 13's ν₂²+ν₂² → ν₂²+ν₃²
(`liang-2024-eq13-denominator`; Eq. 13 unused — no defect-Δv optimisation in
the reproduction route).

## 2. Tolerance rationale (derived before asserting)

Per-flyby tolerance: `tol(t, moon) = 2·v_moon·(5e-5 rad/day)·t + 1e-3 km/s`.

* **Dominant term — Table 1 print quantization:** mean motions are printed to
  4 decimals (half-ULP 5e-5 rad/day), so a moon's polar angle at epoch t
  carries up to 5e-5·t rad of irreducible uncertainty; that rotates both the
  Lambert endpoint and the moon velocity, each contributing ≤ v·δθ to the V∞
  difference. Over the flyby epochs (19–119 d) this gives 0.017–0.138 km/s.
  Sensitivity scan (Member C, the worst member): perturbing a SINGLE mean
  motion by ±5e-5 swings the worst residual between **1.6e-2 and 1.5e-1 km/s**
  (n_E +5e-5 → 1.60e-2; n_E −5e-5 → 1.01e-1; n_C +5e-5 → 1.46e-1; n_G ±5e-5 →
  ~5.1e-2) — bracketing the observed 4.82e-2. The misfit is fully explained by
  input print precision.
* **μ_Jupiter (unprinted, mining note §9.4):** velocities scale ~μ^(1/3) at
  fixed mean motions/ToFs. Scan: μ_J ∈ {1.266e8, 1.267e8, 1.26712764e8}
  moves the worst Member C residual only between 4.74e-2 and 5.09e-2 km/s
  (~±2e-3) — an order of magnitude below the mean-motion term, confirming the
  registry-GM choice is immaterial. Folded into the 1e-3 km/s floor together
  with the Table 2/4/6 phase and Tables 3/5/7 ToF/V∞ print quantizations.
* The mining-note §5 guidance "validate V∞ at ~1e-3 km/s" assumed μ_J was the
  only gap; it did not account for the accumulated mean-motion quantization,
  which dominates by ~50×. The 1e-3 level IS asserted where no time
  accumulation applies: the initial-conic flyby-0 V∞ (pure (a, e, r_p) conic
  geometry, no Lambert): residuals **1.6e-5 (A/B)**, **4.1e-4 (C)** km/s.
* Tolerances were not widened at any point to admit a member; a member failing
  them would have been recorded as a documented negative.

## 3. Per-member residual tables

Epochs in days from perijove departure (t = 0). `in`/`out` = reconstructed
|V∞| inbound/outbound; `printed` = Tables 3/5/7; `cont` = |in − out|
(ballistic ⇒ 0); `tol` = per-flyby tolerance (§2). "alt~" = our reproduction
of the paper's no-Jupiter defect-altitude FICTION (p. 16 caveat; informational
only, registry moon GM/radius, never a catalogue anchor).

### Member A — 1-1-1, high perijove (r_p = 660,993.5 km, e = 0.648881, t_c0 = 2.44865 d)

| # | Moon | t (d) | printed | in | out | \|in−p\| | \|out−p\| | cont | tol | alt~ (printed) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | Callisto | 19.1371 | 5.6730 | 5.6730 | 5.6705 | 1.6e-5 | 2.5e-3 | 2.5e-3 | 1.7e-2 | 1.05e6 (1.90e6) |
| 1 | Ganymede | 51.0344 | 6.9919 | 6.9870 | 6.9816 | 4.9e-3 | 1.0e-2 | 5.3e-3 | 5.7e-2 | 6.9e6 (9.8e5) |
| 2 | Callisto | 69.2041 | 5.6698 | 5.6643 | 5.6724 | 5.5e-3 | 2.6e-3 | 8.1e-3 | 5.8e-2 | 32,578 (33,839) |
| 3 | Europa | 99.1384 | 4.6685 | 4.6823 | 4.6837 | 1.4e-2 | 1.5e-2 | 1.4e-3 | 1.4e-1 | 6,403 (6,241) |
| 4 | Callisto | 119.1131 | 5.8721 | 5.8746 | — | 2.5e-3 | — | — | 9.9e-2 | — (19,765) |

Cycle ToF 99.9760 d (printed Fig. 5d range 99.86–100.14). Legs: revs/branch =
1/high, 1/low, 1/high, 1/low; selection margins 1.14, 2.90, 2.75, 3.16 km/s.
**Worst residual 1.52e-2 km/s (flyby 3 out). REPRODUCES.**

### Member B — 1-1-0, high perijove (same conic as A; diverges at the Europa branch)

| # | Moon | t (d) | printed | in | out | \|in−p\| | \|out−p\| | cont | tol | alt~ (printed) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | Callisto | 19.1371 | 5.6730 | 5.6730 | 5.6705 | 1.6e-5 | 2.5e-3 | 2.5e-3 | 1.7e-2 | 1.05e6 (1.90e6) |
| 1 | Ganymede | 51.0344 | 6.9919 | 6.9870 | 6.9816 | 4.9e-3 | 1.0e-2 | 5.3e-3 | 5.7e-2 | 6.9e6 (9.8e5) |
| 2 | Callisto | 69.2041 | 5.6698 | 5.6643 | 5.6721 | 5.5e-3 | 2.3e-3 | 7.8e-3 | 5.8e-2 | 57,217 (60,877) |
| 3 | Europa | 99.4891 | 4.4853 | 4.4978 | 4.4990 | 1.3e-2 | 1.4e-2 | 1.2e-3 | 1.4e-1 | 11,215 (10,825) |
| 4 | Callisto | 119.1725 | 5.7914 | 5.7935 | — | 2.1e-3 | — | — | 9.9e-2 | — (36,258) |

Cycle ToF 100.0354 d. First two flybys bit-identical to Member A (the paper's
own observation for Table 5). Legs 1/high, 1/low, 1/high, 1/low; margins 1.14,
2.90, 3.05, 3.43 km/s. **Worst residual 1.38e-2 km/s (flyby 3 out). REPRODUCES.**

### Member C — 1-1-1, LOW perijove (r_p = 330,496.8 km, e = 0.824440, t_c0 = 1.98235 d)

| # | Moon | t (d) | printed | in | out | \|in−p\| | \|out−p\| | cont | tol | alt~ (printed) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | Callisto | 18.6708 | 7.6433 | 7.6429 | 7.6394 | 4.1e-4 | 3.9e-3 | 3.5e-3 | 1.6e-2 | 5.5e5 (1.33e6) |
| 1 | Ganymede | 50.9250 | 10.4922 | 10.4853 | 10.4763 | 6.9e-3 | 1.6e-2 | 9.0e-3 | 5.6e-2 | 5.4e5 (6.3e5) |
| 2 | Callisto | 68.7377 | 7.6409 | 7.6316 | 7.6276 | 9.3e-3 | 1.3e-2 | 4.0e-3 | 5.7e-2 | 28,730 (27,438) |
| 3 | Europa | 99.8644 | 12.0213 | 11.9805 | 11.9731 | 4.1e-2 | 4.8e-2 | 7.4e-3 | 1.4e-1 | 3,383 (3,636) |
| 4 | Callisto | 118.8006 | 7.7838 | 7.7677 | — | 1.6e-2 | — | — | 9.9e-2 | — (12,516) |

Cycle ToF 100.1298 d. Legs 1/high, 1/low, 1/high, 1/low; margins 1.00, 1.79,
1.21, 2.64 km/s. **Worst residual 4.82e-2 km/s (flyby 3 out). REPRODUCES**
(the doubled-V∞ high-energy contrast member also carries the largest
quantization sensitivity — consistent with §2's perturbation bracketing).

### Altitude-fiction cross-check (informational)

For the tight flybys (#2, #3; turn angles O(0.5–1 rad)) our reproduction of
the paper's no-Jupiter altitude convention lands within ~10% of the printed
columns across all three members (ratios 0.93–1.05) — corroborating the
turn-angle structure. The near-zero-turn flybys (#0, #1) diverge as expected
(the fiction is hyper-sensitive as δ → 0; 1.05e6 vs 1.90e6 km etc.); they are
deliberately not compared. None of these values is ingested anywhere
(paper's p. 16 caveat; mining note §9.3).

## 4. What the reproduction does and does not establish

* DOES: the printed (Table 2/4/6 phases, Table 3/5/7 leg ToFs) tuples, run
  through the paper's own idealized model with the paper's own construction
  conventions, return the printed per-flyby V∞ on BOTH sides of every flyby
  and V∞-magnitude continuity at print precision — the full first cycle of
  each member is internally consistent and ballistic-compatible. All binding
  printed constraints are in the residual (every table V∞, both flyby sides,
  continuity, the flyby-0 conic anchor, the phase convention, the cycle ToF
  envelope); no subset-closure.
* DOES NOT: re-demonstrate the 10-cycle (~1000 d) continuation (the paper
  prints per-flyby data for the first cycle only; later cycles exist as
  figure traces) nor the ephemeris Member D (no numeric per-flyby anchors —
  V0 ceiling stands, mining note §4.4). The rows' "≥ 1000 d demonstrated, not
  indefinite" caveat is unchanged.
* Model-discipline notes: same-model golden (their idealized model, not real
  ephemeris); altitude fiction reproduced in their convention, never "fixed";
  both in-print equation defects honoured from the errata ledger.

## 5. PROPOSED WRITEBACK (HELD FOR REVIEW) — do not apply without user sign-off

V0 → **V1** for all three idealized rows (criterion: reproduced in-repo from
printed inputs in the source's own model; this is the mining note §5
"candidate for a reproduction script + V1 on ingest" milestone — V2 would
additionally require matching Table 3/5/7 at full print precision, which the
4-dp Table 1 mean motions make information-theoretically unreachable, §2).

For each of `liang-2024-cgcec-111-highperijove`, `liang-2024-cgcec-110-highperijove`,
`liang-2024-cgcec-111-lowperijove`:

1. `validation_level: V0` → `validation_level: V1` with inline comment:
   `# V1: reproduced in-repo in the source's own idealized model (#222, commit 70a1cc5); see docs/notes/2026-06-13-liang-abc-reproduction.md`
2. Append to `notes` (per member, worst figures from §3):
   * A: `REPRODUCED 2026-06-13 (#222): same-model reconstruction (cge_scaffold.py) matches Table 3 per-flyby V_inf on both flyby sides to worst 1.52e-2 km/s and V_inf continuity to 8.1e-3 km/s, inside the Table 1 print-quantization tolerance (0.017-0.14 km/s per flyby); all 4 legs 1-rev multi-rev Lambert as published; golden tests tests/search/test_liang_cge_reproduction.py.`
   * B: same text with `Table 5`, worst `1.38e-2`, continuity `7.8e-3`.
   * C: same text with `Table 7`, worst `4.82e-2`, continuity `9.0e-3`.
3. No other field changes: anchors, data_gaps (mass_ratio n/a; μ_J unprinted),
   tof_days_bounds, source_ephemeris: null all remain correct as written.
4. NOT proposed: any change to `liang-2024-cgcec-ephemeris-2033` (Member D
   stays V0 — no numeric per-flyby anchors exist to reproduce against).

## 6. Follow-ups surfaced (no action taken)

* The reproduction-led sweep gate of mining note §6.2 now has its seed
  validated: Member A's reconstruction is a working same-model scaffold for
  perturbation campaigns (other 0/1 structures, perijoves, switched-pair
  sequences, Saturnian analogue) once multi-rev/repeated-moon genes land in
  the genome.
* The t_c1 = T + t_c0 anchoring subtlety (§1.2) is exactly the kind of
  unprinted construction convention worth checking first in any future
  Lynam/Longuski or Hernandez triple-cycler reproduction (mining note §8
  acquisitions).
