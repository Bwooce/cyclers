# `#822`: Vaquero 2:1 <-> 3:1 Earth-Moon free transfer — computed and verified at all 13 grid points

**Task:** `#822`, registered 2026-08-10 (found during `#811`), dispatched 2026-08-11. Compute an
actual free transfer — an unstable-to-unstable, same-Jacobi-constant heteroclinic invariant-manifold
connection — between Vaquero (2013)'s own 2:1 and 3:1 Earth-Moon "Periodic Cycler" families
(Purdue Ph.D. dissertation, Sec. 4.4.7, pp. 169-172; corpus md5 `fdcbf871322b87cd1dd3448059cb2596`)
in their printed Jacobi-constant overlap band `C ∈ [2.54, 2.66]`.

**Verdict: COMPUTED AND VERIFIED — a genuine, independently-Radau-cross-checked, zero-Δv
heteroclinic connection Wu(2:1) → Ws(3:1) was converged and passed the full verification battery
at ALL 13 shared 0.01-grid Jacobi constants of the overlap band, plus the reverse direction
Wu(3:1) → Ws(2:1) demonstrated directly at C=2.60.** This REPRODUCES (supplies the missing
numbers for) Vaquero's own qualitative pp. 171-172 assertion that a free transfer between the
two families exists in exactly this band. **Nothing here is novel and nothing is claimed
novel** — the existence claim is Vaquero's own, and quantitative 3:1→2:1 heteroclinics in the
Earth-Moon PCR3BP are independently published at OTHER energies (`C_J ∈ [3.00, 3.15]`) by
Kumar, Rawat, Rosengren & Ross (*Adv. Space Res.* **77**(3):3815, 2026, DOI
`10.1016/j.asr.2025.12.005`, arXiv:2509.12675 — already in the corpus, digested 2026-06-20).
The mandatory `literature_check.py` gate was run (Sec. 6): verdict `published` (confidence
0.95), anchored on the Casoliva 2010 / Vaquero lineage — exactly consistent with this task's
reproduction framing.

Module: `src/cyclerfinder/search/vaquero_em_cycler_connections.py`
(tests: `tests/search/test_vaquero_em_cycler_connections.py`, 12 tests; driver:
`scripts/screen_822_vaquero_em_free_transfer.py`; full per-C record incl. every evidence field:
`data/found/822_vaquero_em_free_transfer/results.json`).

---

## 1. What Vaquero actually asserts, and what was therefore reproducible

Vaquero prints NO transfer trajectory, NO intersection state, and NO Δv for this connection
anywhere in Sec. 4.4.7 — only the prose existence claim (pp. 171-172, including the
"2.66 < C < 2.54" typesetting slip `#787`'s digest documents) plus its precondition: both
families unstable at shared C. `#799` already established the precondition quantitatively (2:1
`|λ| ∈ [3.15, 5.73]`, 3:1 `|λ| ∈ [11.26, 12.78]` across the band). With no published digits,
success is a genuine connection with this project's own full self-consistency battery — the same
honest framing as `#766`/`#767`/`#781` (targets without published digit-grade states), NOT a
number-matching exercise.

## 2. Method (why not `#783`'s lane, and what was reused)

`#783`'s `earth_moon_resonant_connections.py` implements Barrabés-Mondelo-Ollé full-state
continuation shooting, and its clean negative was on a categorically different target (Casoliva's
own hardest He1 HOMOCLINIC point, cold-started — the "conditioning wall"). Method fit here was
treated as open, per the dispatch. The lane that fits is this project's own Poincaré-section
Newton machinery (`genome/heteroclinic_cycle.correct_connection`, the `#754`/`#759`/`#767`/`#781`
chain): a 2-D Newton on the manifold phases `(tau_u, tau_s)` with the section gap `(Δx, Δxdot)`
on `{y=0}` as residual — reused unchanged. In the event it converges cleanly everywhere in the
band.

New module pieces (all promoted to first-class, tested code):

- **Node re-derivation** (`build_vaquero_overlap_node`): each family member at the target C is
  re-converged at call time by the same fixed-Jacobi corrector `#799` used (`tol=1e-12`,
  `half_crossings=3` topology asserted), from `#799`-archived guesses vendored as
  DERIVE-provenance constants; the recomputed Floquet magnitude is checked against the archived
  `|λ|` to 1e-6 relative (staleness guard). Both nodes land at exactly the target C — the
  equal-energy precondition holds to machine precision.
- **Whole-crossing-set manifold curves** (`manifold_section_crossings`): `#781`'s scratchpad
  close-approach method, promoted: one propagation per phase-grid point collects EVERY `{y=0}`
  crossing (all `x`/`ydot` signs, all k) out to `n_periods` — O(n_tau) integrations per (node,
  branch) instead of an intractable `(branch, k_u, k_s)` product scan (the connections live at
  crossing indices k~27-54, i.e. ~4.5-9 orbit periods per leg; a brute-force k-grid of
  grid-scanned `correct_connection` calls would take days).
- **The ydot-sign gate** (`find_connection_seeds` + `verify_connection`): an `(x, xdot)` match on
  `{y=0}` at fixed C determines only `|ydot|`, NOT its sign — the 2-D section residual could in
  principle "converge" onto two crossings with opposite `ydot`. Candidate seeds are pre-filtered
  to same-`ydot`-sign pairs, and every verified connection is additionally gated on the FULL
  planar 4-state gap (which contains `Δydot` explicitly) plus an explicit `ydot`-sign match.
- **Seed diversity** (`find_free_transfer`): the distance-ranked top of the seed list is often
  `n_tau`-adjacent near-duplicates of ONE close approach; without a diversity filter a fragile
  geometry burns the whole refinement budget on the same failure (observed directly at C=2.61:
  1174 seeds, 8 refined, 0 converged — all eight the same crossing neighborhood; with the
  0.05-radius diversity skip, the same C verifies). This turned 5 of the first-pass "negatives"
  into verified connections.
- **Verification battery** (`verify_connection`, all hard gates): full-4-state gap ≤ 1e-4;
  ydot signs match; ghost distance to BOTH orbits' own (all-sign-union) section points >
  `GHOST_GUARD_DELTA=1e-3`; independent-Radau (`rtol=atol=1e-11`) re-derivation of both legs'
  crossings, gap and disagreement vs. DOP853 ≤ 1e-4; backward re-approach ≤ 1e-5; forward
  re-approach ≤ 0.5; per-leg Jacobi drift ≤ 1e-9; seed energy-surface offset ≤
  `max(1e-9, 100·ε²)` (below). Gate ceilings are calibrated generously against the observed
  evidence floor — the e-6..e-5 floor is integrator-tolerance-limited (local error ~1e-11
  amplified by the legs' own Floquet growth), the same arithmetic as `#783`'s conditioning
  wall, benign here because legs are re-derived, not shot through.

The 2:1 family's saddle is genuinely NEGATIVE real across the whole band (λ ∈ [-5.73, -3.15];
`#781`'s 4:5-saddle situation — branch labels are related by the period map's own sign flip);
the 3:1's is positive real. No special-case code needed; recorded for correct branch reading.

### Two mid-task method findings worth recording

1. **Manifold-offset ε as an evidence-quality control, not a free knob.** At the sibling-module
   default `ε=0.5e-5`, high-C legs span ~7-9 periods of amplification and the re-derivation
   evidence floor rises to ~1-2e-4 — OVER the (unchanged) gates, so converged Newton
   intersections at C=2.61/2.63/2.65 initially failed verification honestly. Raising ε to 1e-4
   shortens the legs (k~27-37 instead of 42-54), pulling the floor back to e-7..e-5 and letting
   the SAME gates pass. The gates were never loosened; ε is recorded per row. (At the C=2.54
   band edge the same move also fixes the horizon problem: the 2:1 saddle is weakest there,
   `|λ|=3.15`, and an ε=0.5e-5 offset cannot reach O(1) separation within 9 periods.)
2. **Jacobi-gate metric semantics** (`feedback_verify_metric_semantics_before_ranking` applies
   to gates too): the original gate compared the crossing's Jacobi to the NODE's, which
   conflates true integration drift (~1e-12..3.5e-11 everywhere, gated at 1e-9) with the
   linear manifold seed's own O(ε²) energy-surface offset. The saddle eigenvector is
   energy-orthogonal to FIRST order (energy invariance of the linearized period map:
   `∇C·Mv = ∇C·v`, so `(λ-1)∇C·v = 0`), leaving a quadratic offset — measured ~6e-12 at
   ε=0.5e-5 and ~2.6e-9 at ε=1e-4, the clean (20x)² = 400x scaling. `verify_connection` now
   measures drift against the SEED's own conserved value and gates the construction offset
   separately at `max(1e-9, 100·ε²)`.

## 3. The primary transfer (C = 2.60, mid-band)

`Wu(2:1) → Ws(3:1)`, `branch_u=+1, branch_s=+1`, `k_u=45, k_s=30`, `ε=0.5e-5`:

```
node 2:1: x0=1.0666538160256032  ydot0=-0.8231049641748784  T=5.841315589133051   C=2.6 (exact)
node 3:1: x0=0.8970164526223031  ydot0=-0.8032821103215304  T=6.276485216101793   C=2.6 (exact)
converged: tau_u=3.6477847...  tau_s=5.8460807...   Newton residual 8.0e-11
crossing (y=0): x=+0.865972  xdot=-0.479921  ydot=-0.607359  (ydot<0 BOTH legs)
transit: t_u=+43.45 nd (188.7 d)   t_s=-30.84 nd (133.9 d)
Δv = 0 (heteroclinic — position AND velocity continuous at the match)
```

Evidence: full 4-state gap `5.7e-6`; Radau independent gap `4.1e-6` (worst vs-DOP853 `4.9e-6`);
ghost margins `0.520 / 0.481` (520x/481x the guard); backward re-approach `5.5e-8`; forward
`2.0e-2`; per-leg drift `2.1e-12 / 1.7e-11`; seed offsets `6.5e-12 / 4.7e-12`. Two further
independent connections converged and fully verified at this C during exploration (`k=(42,32)`
at `x=-0.2244, xdot=-1.6068`; `k=(46,30)` at `x=-0.3640, xdot=+1.4740`) — the intersection
structure is rich, not an isolated point. The primary crossing sits near-but-genuinely-off the
3:1 orbit's own neighborhood — the expected "intersection near the target orbit" geometry (cf.
Anderson & Lo 2011's own Table 3 "near the 3:4 orbit").

**Transit-time caveat**: manifold transit clocks are ε-convention-bound (trajectories are
asymptotic to both orbits; times are measured from/to the ε-scale seed points), so the ~150-370 d
totals below are transits of the ε-parametrized representatives, not unique physical transfer
times. Vaquero prints no ToF to compare against.

## 4. Band sweep — all 13 shared grid points VERIFIED, plus the reverse direction

Driver settings per row (recorded in `results.json`): `n_tau=36` (2.62: 54), `n_periods=9`
(2.62: 10), `max_refine=12` (2.62: 14), `tol=1e-9`, ε per the strategy in Sec. 2. All 13 rows
pass every gate of the battery; residuals 8.0e-11..8.7e-10; ydot signs match on every row;
per-leg drift ≤ 3.5e-11 everywhere.

| C | ε | branches | (k_u,k_s) | residual | crossing (x, xdot) | 4-state gap | Radau gap | ghost (from, to) | t_u + t_s (d) |
|---|---|---|---|---|---|---|---|---|---|
| 2.54 | 1e-4 | (-1,+1) | (40,27) | 3.3e-10 | (+0.233225, -1.785596) | 4.5e-6 | 1.8e-6 | 0.243, 1.800 | 168.7+120.0 |
| 2.55 | 5e-6 | (-1,+1) | (54,31) | 2.6e-10 | (-1.074291, +0.239636) | 1.0e-5 | 1.0e-5 | 0.258, 1.108 | 228.1+140.4 |
| 2.56 | 5e-6 | (-1,+1) | (49,32) | 6.8e-10 | (-0.328920, +1.620448) | 3.2e-5 | 1.9e-5 | 0.110, 1.655 | 208.2+144.4 |
| 2.57 | 5e-6 | (-1,+1) | (48,30) | 2.3e-10 | (+0.866609, -0.489803) | 2.6e-6 | 1.8e-6 | 0.534, 0.491 | 203.6+134.4 |
| 2.58 | 5e-6 | (+1,+1) | (46,28) | 2.8e-10 | (+0.919769, -0.324518) | 1.9e-6 | 1.9e-6 | 0.360, 0.325 | 197.6+126.0 |
| 2.59 | 5e-6 | (-1,+1) | (42,31) | 2.1e-10 | (-1.058394, +0.264720) | 1.2e-5 | 1.1e-5 | 0.283, 1.103 | 176.3+139.7 |
| 2.60 | 5e-6 | (+1,+1) | (45,30) | 8.0e-11 | (+0.865972, -0.479921) | 5.7e-6 | 4.1e-6 | 0.520, 0.481 | 188.7+133.9 |
| 2.61 | 1e-4 | (-1,+1) | (27,28) | 1.1e-10 | (+0.399809, +1.255917) | 5.4e-7 | 5.3e-7 | 0.110, 1.314 | 115.6+120.2 |
| 2.62 | 1e-4 | (+1,+1) | (31,34) | 2.5e-10 | (-0.281971, +0.889276) | 8.0e-6 | 7.6e-6 | 0.482, 0.938 | 131.9+185.8 |
| 2.63 | 1e-4 | (+1,+1) | (32,30) | 5.7e-10 | (-0.320287, -0.014537) | 3.2e-5 | 3.1e-5 | 0.827, 0.338 | 150.4+150.6 |
| 2.64 | 5e-6 | (-1,+1) | (43,38) | 8.7e-10 | (+0.504232, +1.003143) | 4.3e-5 | 3.6e-5 | 0.288, 1.076 | 214.1+170.7 |
| 2.65 | 1e-4 | (+1,-1) | (44,44) | 7.7e-10 | (+0.869008, -0.447578) | 1.7e-7 | 1.9e-7 | 0.479, 0.448 | 222.4+182.7 |
| 2.66 | 1e-4 | (+1,+1) | (37,32) | 4.5e-10 | (+0.869164, -0.575466) | 6.5e-6 | 3.5e-6 | 0.599, 0.576 | 180.5+150.6 |

Reverse direction `Wu(3:1) → Ws(2:1)` at C=2.60 (ε=1e-4): VERIFIED — residual `8.1e-10`,
crossing `(-0.286956, -1.354329)`, `k=(29,32)`, branches `(-1,+1)`, 4-state gap `1.1e-5`,
Radau gap `4.9e-6`, ghost `1.387 / 0.069`, transit 124.2+140.8 d. The CR3BP time-reversal
symmetry `(x, y, xdot, ydot, t) → (x, -y, -xdot, ydot, -t)` maps `Wu(A)∩Ws(B)` connections into
`Wu(B)∩Ws(A)` ones, so the reverse's existence is guaranteed given the forward hit; it is
demonstrated directly anyway (never merely asserted).

Honest search-history note (recorded per `orbit-closure discipline`, not hidden): the first
sweep pass (no diversity filter, ε=0.5e-5 everywhere, `max_refine=8`) verified only 5 of 13
points. The three fixes above (seed diversity, per-region ε, the Jacobi-gate metric correction)
each flipped specific honest negatives to verified connections — every fix is a diagnosed
mechanism, not a tolerance loosened to force a pass; the gates themselves never moved. The
intermediate negatives and their diagnoses are preserved in this note and the git history.

## 5. Relation to Kumar-Rawat-Rosengren-Ross 2026 (the closest published quantitative work)

Their published 3:1→2:1 heteroclinics sit at `C_J ∈ {3.00, 3.05, 3.10, 3.15}` — OUTSIDE
Vaquero's [2.54, 2.66] band, with orbit ICs anchored at those higher energies — so `#822`'s
in-band numbers do not duplicate their tables; they corroborate the same transport mechanism at
neighboring energies. Grounded against the paper's own digest
(`docs/notes/2026-06-20-digest-kumar-2025.md`) per `feedback_ground_citations_against_content`;
their Table-6 orbit ICs already reproduce to ~1e-13 in this project (`resonance_network.py`,
`#598`). Their Table-5 printed intersection states are a true digit-grade reproduction target
for this machinery — registered as `#827`.

## 6. Literature gate (mandatory floor — run, not skipped, despite reproduction framing)

Signature: `CandidateSignature(primary="Earth", sequence=("Moon",), resonances=("2:1", "3:1"),
topology_label={"resonant", "repeated-moon"})`. `build_queries` generated 36 queries; run
against the real WebSearch tool this task (live queries including "Moon cycler resonance 2:1
3:1" and "two classes cycler trajectories Earth-Moon system"), with the retrieved results
replayed into `check_literature`. Verdict: **`published`, confidence 0.95**, citation Casoliva
et al. 2010 (DOI `10.2514/1.46856` — the curated corpus anchor for this exact Earth-Moon
resonant-cycler lineage; Vaquero 2013 is the dissertation continuation of it), "NOT
novelty-claimable" — the correct verdict for a reproduction. The live search also independently
surfaced the Kumar et al. 2026 paper (Sec. 5).

## 7. Follow-ups registered

- `#827`: digit-grade reproduction of Kumar-Rawat-Rosengren-Ross 2026's own printed Table-5
  3:1→2:1 heteroclinic intersection states (`C_J ∈ {3.00, 3.05, 3.10, 3.15}`) with `#822`'s
  machinery — unlike Vaquero, they print exact intersection states.
- `#828`: adjudicate whether `#822`'s verified connections upgrade any `#811` catalogue rows'
  validation level (the `#822` registration's own "natural upgrade path toward V2/V3-class
  connection evidence"), and if so perform the writeback under full catalogue-edit ratchet
  discipline.

## 8. Verification run

- `tests/search/test_vaquero_em_cycler_connections.py`: 12 tests (sourced-range derivation,
  node re-derivation + saddle signs, the primary hit's convergence + full battery with real
  margins incl. the ydot-sign and seed-offset gates, seed plumbing, honest-rejection paths) —
  all pass.
- Full `uv run pytest tests/data tests/search tests/scripts -q`: exit 0, zero FAILED/ERROR
  (pre-existing documented XPASS cross-platform notes only).
- `uv run ruff check .`, `uv run ruff format --check .`, full `uv run mypy src tests`: all clean.
- Full sweep record: `data/found/822_vaquero_em_free_transfer/results.json` (13 forward rows +
  reverse demo, every row re-recorded under the FINAL code version for uniform evidence shape).
