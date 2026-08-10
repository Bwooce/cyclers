# #818 — reusable passive-node self-consistency gate + re-run over `#816`'s 785 passive roots (2026-08-10)

**Verdict: gate built, tested, and run. Joint result CONFIRMS `#817`'s operative claim — NO
passive-node root survives once the gate joins the battery (0/785), and the gate alone
auto-rejects every root that any existing gate could not (all 6 turn-feasible passive roots,
including both entries of `#817`'s hand-adjudicated object at 33%). But `#817`'s literal
universal hypothesis ("ALL 785 fail the same way") is REFUTED under the gate's declared
fraction-of-budget criterion: 698/785 pass it, because their working-budget denominator is a
physically impossible 98°–178° turn. Under the complementary ABSOLUTE reading, all 785/785 do
fail. Both numbers reported honestly below.**

Pure post-processing of stored numbers; no physics re-run, no catalogue write, no stamp.

## 1. The gate

`cyclerfinder.search.physical_sanity.passive_node_is_self_consistent(body, vinf_kms,
working_turn_budget_deg, *, max_parasitic_fraction=0.02) -> PassiveNodeVerdict`

- **Deflection law reused verbatim**: `cyclerfinder.core.flyby.max_bend` (BMW §6.4,
  `sin(δ/2) = 1/(1 + r_p·V∞²/μ)`) — the SAME implementation whose safe-altitude values `#817`
  verified reproduce the stored `max_bend_deg_per_encounter` bit-for-bit. The gate evaluates
  it at the body's **Laplace SOI** radius `a·(μ/μ_p)^(2/5)` instead of the safe periapsis.
  A parity test (`test_817_deflection_law_is_core_flyby_max_bend_verbatim`) pins the gate to
  that implementation at rel 1e-15.
- **Why the SOI boundary**: any pass that counts as an encounter in patched-conic semantics
  has periapsis inside the SOI, and deflection is monotonically decreasing in `r_p`, so
  `δ(r_SOI)` is the *minimum* deflection any real encounter imparts. Rejecting on the lower
  bound makes the rejection robust. The verdict record also carries the Hill-radius value
  (`a·(μ/(3μ_p))^(1/3)`, the more lenient bound, `#817`'s 24% row) informationally.
- **Which radius the project uses** (task question): both existed. `data/transfer_network.py::
  r_soi_km` implements the true Laplace SOI (design §5); `search/tour_self_consistency.py::
  soi_km` (`#480`'s geometric rule) implements the **Hill** formula while its docstring said
  "Laplace SOI" — a docstring mislabel, corrected in this task (behavior unchanged; for `#480`'s
  containment purpose the lenient Hill bound is fine and arguably right). The new gate uses the
  Laplace SOI as verdict-bearing, matching `#817`'s own decisive 33% row and `transfer_network`'s
  precedent, and reports Hill alongside. `laplace_soi_km` / `hill_radius_km` handle both
  planets (about the Sun) and registry satellites (about their primary).
- **Reject rule**: inconsistent unless
  `δ_parasitic(SOI) ≤ max_parasitic_fraction · working_turn_budget`, where the budget is the
  required turn at the trajectory's working (non-passive) node(s). Non-positive budget →
  automatic reject (fraction ∞).
- **Passivity criterion**: `PASSIVE_NODE_TURN_MAX_DEG = 0.05°`, the same value as
  `scan_816_unequal_tof_discrete_roots.py`'s `TURN_TRIVIAL_DEG`, asserted equal in a test so
  gate and classifier can never silently disagree on what "passive" means.

## 2. Threshold: 2% of the working turn budget — and why

`DEFAULT_MAX_PARASITIC_TURN_FRACTION = 0.02`. A judgment threshold (NOT a sourced physical
constant), grounded in the project's own precedents rather than invented fresh:

- `#817`'s adjudication itself used **"under 2% of the working turn"** as its negligibility
  line (the 160,000-km-standoff row of its table), and still treated 6.5% as disqualifying-
  side. This gate adopts the adjudication's own line.
- Same spirit as `#324`'s 5° useful-bend floor: a parasitic turn at ~2% of the working budget
  is at the level of targeting/TCM noise; above it the closure's arithmetic is materially wrong.
- Calibration against the two anchor cases: the Russell-Strange-2009 genuinely-negligible
  passive-target architecture (Enceladus at ~4 km/s vs a tens-of-degrees Titan working turn)
  sits below 0.5%; `#817`'s Oberon case is 33%. Two orders of magnitude of separation — the
  verdict on both anchors is insensitive to the exact value over roughly 0.5%–10%.
- Sensitivity on `#816`'s data (see §4): the joint outcome (0 survivors) is IDENTICAL at 2%,
  5%, or 10%, and identical under the lenient Hill-radius bound. The only roots whose
  *per-gate* verdict flips between 2% and 5% (the 2.85% group, §4) are already dead on the
  `#324` bend gate.

## 3. Controls (tests/search/test_physical_sanity.py, `test_817_*` / `test_818_*`)

- **Positive control** (`#817`'s object, stored numbers): Oberon passive at V∞ 1.31114 km/s
  vs the 4.2327790° Ariel working turn → Laplace SOI 9,678 km, parasitic **1.3968°** =
  **33.0%** → REJECT; Hill 13,288 km, 1.0207° = 24% — every figure matching `#817`'s table.
- **Negative control** (R-S 2009-style, `#817` §2's own calibration case): Enceladus passive
  at 4 km/s vs a 30° Titan budget → parasitic ~0.11° ≈ 0.37% → ADMIT. The gate does not
  over-reject the legitimate passive-science-target architecture.
- Plus: verbatim-law parity (rel 1e-15), Hill parity with `tour_self_consistency.soi_km`,
  Laplace < Hill ordering for all registry moons, planet lookup (Earth SOI ≈ 9.24e5 km
  textbook value), zero-budget reject, unknown-body `KeyError`, input validation, and the
  borderline 2.85% `#816` object (§4) rejected at the 2% floor.

## 4. Re-run over `#816`'s 785 stored passive-node roots (pure post-processing)

`scripts/postprocess_818_passive_node_gate.py` → `data/found/818_passive_node_gate/
gate_results.json`. Inputs are the stored per-root `(passive body, V∞, required turns)`;
for `anchor_passive` roots the larger of the two stored anchor V∞ entries is used (higher V∞
→ smaller parasitic deflection → most lenient bound, so every rejection is robust).

| population | n | gate verdict |
|---|---|---|
| all passive-node roots | 785 | 87 REJECT / 698 ADMIT at 2% |
| turn-feasible passive roots (working node CAN deliver its turn) | 6 | **6/6 REJECT** |
| passive roots passing ALL `#816` physical gates | 2 | **2/2 REJECT** (both at 33.0%) |
| **joint survivors (physical gates AND this gate)** | **0** | — |

Cross-tab `(turn_feasible, passes_physical_gates, self_consistent)`:
`(F,F,T)` 698 · `(F,F,F)` 81 · `(T,F,F)` 4 · `(T,T,F)` 2.

The 6 turn-feasible rejects are **two physical objects** (each found in multiple
direction/mirror conventions):

1. Ariel-Oberon `q=13 n_rev=(2,2)` — `#817`'s adjudicated object. Passive Oberon at
   1.3111 km/s → parasitic 1.3968° = **33.00%** of the 4.2328° Ariel budget. Rejected at any
   threshold ≤ 33%; the gate reproduces `#817`'s verdict automatically.
2. Ariel-Oberon `q=13 n_rev=(3,3)` (4 entries). Passive **Ariel** at 4.6215 km/s → parasitic
   0.2024° = **2.85%** of the 7.0981° Oberon budget. Rejected at 2% with a thin 1.42×
   margin (would be admitted at 5%); independently already dead on the `#324` bend gate
   (achievable bend at the passive Ariel node is 0.71° < 5° floor), so the joint outcome does
   not hinge on the threshold. Reported honestly as the nearest-to-the-line case.

**`#817`'s hypothesis, adjudicated by the numbers:**

- *Literal/universal form* ("that every passive-node root in `#816`'s box fails the same
  way") — **REFUTED under the fraction criterion**: 698/785 (89%) pass the parasitic gate,
  min fraction 0.064%. Not because their parasitic deflections are small in Oberon-object
  terms (range 0.063°–6.58° at SOI; V∞ range 0.68–8.26 km/s, wider than the note's "1–3 km/s")
  but because their working budgets are 98°–178° — an order of magnitude beyond what any
  Uranian moon can deliver, which is exactly why every one of those 698 was already rejected
  by `#816`'s required-turn wall. A big denominator bought by an impossibility is not
  admissibility; it just means this gate is not the one that kills them.
- *Operative form* (the reason `#818` exists: "would have auto-classified `#816`'s object out
  with no adjudication needed") — **CONFIRMED**: the gate rejects every passive root the rest
  of the battery cannot, joint survivors 0/785, and `#816`'s
  `n_passive_node_physical_gate_passers: 2` anomaly is now closed automatically.
- *Absolute reading* (reported, not verdict-bearing): all **785/785** have parasitic
  SOI-boundary deflection above the 0.05° passivity threshold itself — i.e. every stored
  passive root models the node as turning less than the body must physically impart. In that
  self-inconsistency sense `#817`'s "all fail the same way" is true; the fraction-of-budget
  gate is simply the stricter, budget-aware formulation of it.

No exception requiring a new adjudication was found: both turn-feasible objects fail the
gate, and everything else was already dead on pre-existing gates. No follow-up task
registered; no stamp change needed (`#816`'s `n_survivors: 0` stands, now with the passive
branch closed by gate rather than by hand).

## 5. Reusability + discipline

- The gate is genome-agnostic: any future search producing a near-zero-required-turn node
  (planet or registry moon) can call `passive_node_is_self_consistent` with the node's body,
  V∞, and the trajectory's working turn budget. Not wired into existing search modules (per
  the registration's "small" scope); `#819`, if ever dispatched, should call it from day one.
- `catalogue.yaml` / `empty_regions.jsonl` untouched. `data/found/818_passive_node_gate/`
  holds the reproducible per-root results.
- Docstring correction landed in `search/tour_self_consistency.py::soi_km` (Hill formula was
  labeled "Laplace SOI"; behavior unchanged, consumers unaffected — `nbody/jovian.py` derives
  its own Hill value inline, `data/transfer_network.py::r_soi_km` was already correct).
