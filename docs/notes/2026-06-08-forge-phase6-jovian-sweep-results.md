# Forge Phase 6 — first novelty campaign results: Jovian VILM moon-system sweep (task #172)

**Date run:** 2026-06-09
**Verdict:** EMPTY (the most-likely, expected outcome per the design note §5) — a
rigorous, bounded, method-versioned negative. **This is a complete, valid result.**
**Orchestrator:** `scripts/forge_phase6_moon_run.py`
**Empty-region record:** `data/empty_regions.jsonl` (region_id
`jovian-IEG-vilm-2026-06-09`)
**Review queue:** *no entries* (zero SILVER survivors — correct).

---

## The family swept + why

The planet-centric **Galilean Jovian VILM space** — the Io-Europa-Ganymede(-Io)
chain about Jupiter (#76 substrate, never before swept forward). Chosen per the
design note §5 as the single family with the best novelty odds + lowest cost: the
substrate (`SATELLITES`/`PRIMARIES`, the centred ephemeris, the centre-agnostic
corrector, centre-aware Tisserand, the VILM ΔV-floor, the `(model_assumption,
primary)` dedup bucket) all shipped with #76, and its catalogue bucket is
null-numeric so closures read `novel`. The open question the campaign tested:
**does VILM gating + Laplace-resonance phasing surface a bend-feasible Jovian tour
the #76 no-leveraging closure could not?**

## Search extent (bounded)

| field | value |
|---|---|
| topologies | I-E-G-I at `period_k` ∈ {1, 2}, all-single Lambert legs |
| epochs | 64 per topology (span 8 d, centred on the #76 converging seed) |
| points_total | 128 |
| ephem | `circular`, center=`Jupiter`, `mu_central=PRIMARIES[Jupiter]` |
| prune gates | `vilm_dv_floor<=budget` (50 km/s), `linkable(Jovicentric)`, `max_bend_deg feasibility` |
| budget | 50 km/s (generous — above the VILM floor, so the prune does not mask the result) |
| wall | 6.7 s, 8 cores |

## Outcome

- **12 closed** I-E-G-I chains (the corrector converges about Jupiter — the #76
  centre-correctness result, reproduced through the full novelty pipeline).
- **12 / 12 `novel`** against the null-numeric Jovian catalogue bucket (as §3a
  predicted — there is nothing populated in the `(circular-coplanar, Jupiter)`
  bucket to match).
- **0 / 12 bend-feasible.** Every closure lands at V∞ 8.3–26.8 km/s, far above the
  ~6 km/s feasible band the small moons can turn at. Best max-V∞ = 26.8 km/s;
  **gap to the V∞ floor = 20.8 km/s.**
- **12 / 12 REJECTED** (auto-falsified in `evaluate_closure` because
  bend-infeasible — the firewall: a flyby that cannot deliver the required turn is
  an impossible bend, so it never reaches SILVER).
- **0 SILVER, 0 queued, 0 panel-killed.**

This is the **#76 honest-risk finding generalised**: the no-leveraging
single-ellipse genome closes a *higher-V∞* family that is structurally
bend-infeasible; VILM gating did not surface a bend-feasible Laplace-resonant tour
within this genome's reach. The classic empty-set lesson (Hughes/Jones: a thorough
sweep that finds nothing is real science) — captured here as a first-class,
machine-readable, re-sweepable record rather than buried in prose.

## n-body cross-check (Task 5.1)

**No-op.** Zero SILVER survivors → the n-body ARTIFACT rung was not exercised.
Recorded for honesty: nothing reached the rung.

## The firewalls, verified live

- **Dedup firewall:** all 12 closures read `novel` (sparse bucket), yet NONE was
  promoted past REJECTED — bend-infeasibility caps them below SILVER. A re-derived
  closure is not a discovery; GOLD is structurally unreachable by the loop.
- **Re-sweep gate (capability-subsumption, design §6b):** a second run with the
  *same* single-ellipse method SKIPs (`should_sweep → False`, no duplicate
  record). The recorded empty is method-versioned — a later **multi-arc / n-body /
  low-thrust** method (which the partial order ranks strictly more capable) will
  re-sweep this exact region, while a weaker/equal method is correctly skipped.
  The #163-reopens-#137 lesson is encoded as a gate: the negative never
  permanently forecloses a future capability.
- **No catalogue writeback.** The empty-region log and the review queue are both
  non-catalogue (`is_catalogue_source() → False`).
- **No tolerance loosening.** `tol_kms`, the budget, and the bend cap were left at
  their honest values; the empty set was NOT manufactured by relaxing a gate. The
  search_extent validation rejects unbounded negatives, so an over-pruning
  artefact masquerading as empty would have been refused.

## The sharpened hypothesis

The single-ellipse no-leveraging genome is the *only* family that closes the
Galilean I-E-G chain, and it closes structurally bend-infeasible (V∞ ≫ the
turn-feasible band). A bend-feasible Jovian resonant cycler — if one exists —
requires a strictly-more-capable genome: V∞-leveraging (VILM) maneuvers as
*solved DOF* (not merely a prune gate), and/or a multi-arc / low-thrust /
broken-plane representation. Those methods, when built, will re-sweep this region
automatically via the §6b gate. Until then, **the Galilean I-E-G VILM space is
empty of bend-feasible novelty as far as the single-ellipse method can reach** —
a conditional, bounded, reproducible negative, exactly as "empty" is always
conditional on the method (§6a).

## Honest yield, stated plainly

A discovery was **not** expected on the first sweep, and we did not get one. The
value delivered: (a) the first *forward* novelty sweep of a freshly-built space,
(b) a rigorous bounded negative now on disk as a re-sweepable record, and (c)
end-to-end proof the Phase 6 wiring (VILM prune → centred scan → bridge →
signature → match → Axis-A → gauntlet → empty-region emit + re-sweep gate) works
on a non-heliocentric centre. A clean EMPTY first sweep is the success criterion.
