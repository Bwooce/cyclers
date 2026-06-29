# #480 — Level-3 high-fidelity maintenance-ΔV optimizer (SCOPE, not yet executed)

> **APPROACH A EXECUTED 2026-06-30** — see
> `docs/notes/2026-06-30-480-eggie-maintenance-verdict.md`. The cycle-by-cycle maintenance
> lane (`chain_cycles`, validated on Liang #223) was run on a feasibly-discovered real-eph
> EGGIE member (`scripts/eggie_maintenance_480.py`): **ballistic for exactly 2 cycles, then
> large impulses (~170-760 m/s/cycle), cumulative ~3.3 km/s over 10** — robust across
> retarget budgets. Quantitatively reproduces the paper's "ballistic for two cycles, after
> which large impulses are required" claim (a NOVEL curve at level-2, not a printed-number
> reproduction). The positive control (the EIGE ~30 m/s/10 paper figure) hits the 1-rev
> B-plane wall (`...-eige-realeph-maintenance-verdict.md`) so the method's positive control
> is Liang #223. Approach A is COMPLETE at level-2; Approach C (level-3 B-plane NLP) remains
> the only un-done lever and stays "weeks, last resort". The original scope follows.

**Status:** scoping doc for a fresh, optional effort. The #480 science is already settled
(see the 7 dated `docs/notes/2026-06-*-480-eggie-*` notes): construction bug fixed + guarded,
the IEG ballistic-cycler CLASS reproduced unguided in real ephemeris (~0.5 km/s below
Table-4), exact Table-4 member not reproduced. This scopes the ONE remaining non-redundant
direction — quantifying the per-cycle MAINTENANCE ΔV — and is a new project, not a
continuation of the corrector lane (which is exhausted).

## The gap this would close

The paper (Hernandez-Jones-Jesick 2017, AAS 17-608) reports the real-ephemeris tour is
ballistic for ~2 cycles then needs **maintenance ΔV** (the EIGE example grows to **~30 m/s
over 10 cycles**, p.10-11). We cannot currently quote a maintenance number: our
multiple-shooting corrector (FD and analytic-STM, real-eph, correct ballistic seed)
**plateaus at ~0.1-0.4 km/s continuity and does not converge to a periodic orbit**
(`2026-06-30-480-eggie-level3-nbody.md` + its correction; `cc4f241`). 

## The reframe (the crux — why the current approach can't give the number)

A "maintenance ΔV" is a **station-keeping budget**, NOT a single converged periodic orbit.
The paper does NOT close one periodic orbit; it propagates the near-ballistic real-eph tour
forward cycle by cycle in the high-fidelity model and measures the per-cycle impulse needed
to re-acquire the cycler. Our closed-periodic multiple-shooting is the WRONG framing for this
quantity — which is exactly why it plateaus and yields no number. The fix is to measure
maintenance the way the paper defines it, not to force a periodic-orbit closure.

## Approaches

### A — Cycle-by-cycle maintenance (RECOMMENDED)
Forward-propagate the real-eph ballistic EGGIE seed one cycle at a time in high fidelity,
re-targeting each next cycle to the cycler geometry and charging the powered-flyby ΔV;
report maintenance ΔV vs cycle number (the curve to compare to the paper).
- **Reuse (already built + validated):** `nbody/jovian.py::shoot_cycle` / `chain_cycles` /
  `optimize_cycle` — the Liang CGCEC lane that ALREADY does per-cycle chained correction with
  free epoch offsets and inherited inbound V∞ (validated on Liang Member D, #223). This is the
  right machinery; it was built for exactly this chained-maintenance pattern.
- **Seed:** the real-eph ballistic EGGIE member (V∞≈9.0/6.6/7.7, the discovered family member;
  `search/eggie_ballistic.py` + the periapsis-seed construction).
- **Output:** maintenance ΔV per cycle over ~10 cycles; compare to paper EIGE (~30 m/s/10).
- **Effort:** medium (days) — mostly wiring the EGGIE seed into `chain_cycles` (which targets
  CGCEC today) + a per-cycle ΔV accounting + the comparison. Reuses validated infra.
- **Risk:** moderate — `chain_cycles` is CGCEC-shaped; generalising the sequence + the
  re-target objective to EGGIE is the work. The per-cycle re-target may itself hit the same
  flyby-feasibility limits; report honestly if so.

### B — Full-closure with more DOF (cheap secondary lever)
Extend `jovian_shoot`'s free vector from states-only to **[states + node epochs + free ToFs]**
(the heliocentric `nbody/shooter.py::shoot` already has this layout; `jovian_shoot` does not)
+ the analytic STM (`cc4f241`). Test whether epoch/ToF freedom breaks the ~0.1-0.4 km/s
plateau where states-only could not.
- **Effort:** small (extend the free-var packing + residual epoch handling).
- **Risk:** epoch freedom did NOT help in the ideal model (Stage 2b,
  `2026-06-29-480-eggie-stage2-nbody-verdict.md`); may not in real-eph either. Cheap to try
  before committing to A; if it converges, it gives a closed-orbit ΔV directly.

### C — Faithful SNOPT-style high-fidelity NLP (heavy, last resort)
Replicate the paper's full optimization: minimize total ΔV over all flyby epochs + tangential
maneuvers with continuity + B-plane + altitude (25-70,000 km) constraints. We have an NLP
backend (FBS analytic gradients `core/fbs_match_point.py` + scipy/SLSQP; paper used SNOPT).
- **Effort:** large (weeks). **Only if A and B both fail to match the paper.**

## Recommendation
Do **B first** (cheap: it may simply break the plateau and give a closed-orbit number), then
**A** (the definitional match for "maintenance ΔV" + reuses validated chained-cycle infra).
Reserve **C** for a genuine miss against the paper.

## Mandatory disciplines (carry over)
- **Positive control first:** validate the maintenance method reproduces a KNOWN maintenance
  figure (the Liang Member D chained ΔV, or the paper's EIGE ~30 m/s/10) BEFORE trusting an
  EGGIE number ([[feedback_verify_gauntlet_with_positive_control]]).
- **Per-encounter self-consistency** on any reconstructed cycle
  ([[feedback_constructed_tour_per_encounter_self_consistency]]).
- **Foreground-bounded compute**, chunked, runlog-flushed; no detached subagent jobs; `trf`
  not `lm` ([[feedback_long_agents_commit_incrementally]]).
- **Sourced comparison** to the paper's printed maintenance figures; **no exact-Table-4
  reproduction claim** (the member sits ~0.5 km/s below Table-4); **no catalogue self-admission**
  (a reproduced published tour is human-admitted V4-ceiling at most).

## Payoff vs cost (honest)
Payoff: closes the #480 quantitative gap (a maintenance-ΔV curve vs the paper) and a
generally-useful Jovian moon-tour station-keeping capability. Cost: A is days, B is small, C is
weeks. The #480 CLASS reproduction + bug fix already stand; this is refinement with **no
catalogue impact**. Worth doing only if the maintenance-ΔV number is itself a goal — otherwise
#480 is complete as-is.
