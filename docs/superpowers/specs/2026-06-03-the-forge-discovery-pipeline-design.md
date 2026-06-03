# The Forge — Discovery + Intensive Cross‑Validation Pipeline (Design Spec)

**Status:** design (2026-06-03). Foundation = the data‑validation‑hardening plan (Axis C/D); pipeline phases = the companion Forge plan.

## Purpose
Move the project from *reproducing catalogued cyclers* to *discovering new ones with a defensible, multi‑axis confidence tier on every result*. The deliverable is **trust**, not raw discoveries: a candidate is only promoted to the degree that **independent methods, fidelities, and sources agree**, and every gate is **falsifiable**.

## Governing principles (distilled from the 2026‑06 work)
1. **Independence is the validation currency.** Agreement only counts when the things agreeing are independent (different solver / different model fidelity / different published source / different reader).
2. **Construction‑first, not free optimisation.** Resonance‑anchored seeds (`construct_resonant_cycler`) avoid the off‑family basins that defeated S1L1; the optimiser refines, it does not blindly search.
3. **Never compare across fidelity.** Tag and enforce circular‑coplanar vs analytic‑ephemeris vs real‑DE440 (the S1L1 5.65‑vs‑4.99 bug class).
4. **Gates must have teeth.** Deliberately‑wrong input must fail (falsification), or the gate proves nothing.
5. **Novel ≠ sourced.** A discovery with no published counterpart can never be cross‑checked against a source; it is held at "machine‑confirmed, unsourced" until a human (or GMAT) signs off.

## Stages
1. **Enumerate** — `feasible_cells` over (bodies × resonance k × sequence × revs/branch), 3D‑Tisserand pruned, on a deepening frontier with **loop‑until‑dry** so coverage is systematic.
2. **Seed & solve** — construction‑first seed → idealised refine → real‑eph (`optimise_cell_ephemeris`) refine.
3. **Cross‑validation gauntlet** — the intensive core (below).
4. **Classify & score** — validation tier + confidence; novelty match vs catalogue (`canonical_signature`/`match`).
5. **Ledger + ratchets** — every candidate recorded with provenance + re‑runnable derivation trail; frozen exclusion *and* validation‑tier census ratchets.
6. **Human gate** — novel‑confirmed candidates surface for human / GMAT sign‑off before being called a discovery.

## The cross‑validation gauntlet (4 independent axes)
A candidate earns a tier only where the axes that *can* apply agree.

| Axis | Independent checks | Failure mode caught |
|---|---|---|
| **A · Code path** | in‑house Lambert vs lamberthub izzo vs gooding; resonance‑construction vs free‑optimiser vs Kepler re‑propagation | solver / code bugs |
| **B · Fidelity ladder** | coplanar → analytic‑ephemeris → real‑DE440; V∞ persists or varies *predictably* across rungs | model artifacts, cross‑fidelity confusion |
| **C · Source** | ≥2 *independent* papers' anchors when it matches a published cycler | data / transcription error |
| **D · Adversarial** | second‑agent independent re‑derivation; falsification perturbation; digit‑for‑digit re‑extraction of scanned sources | single‑reader error, no‑op gates |

## Validation tiers (the promotion ladder)
- **GOLD (rediscovery, cross‑validated):** matches a published cycler; Axes A+B+C+D all clear (independent sources agree). The strongest claim.
- **SILVER (novel, machine‑confirmed):** no published counterpart; Axes A+B+D clear (method‑independence + fidelity‑persistence + adversarial reproducibility). Explicitly *unsourced* → human gate.
- **BRONZE (consistency‑checked):** closes, but the available checks share a source/fidelity (catches transcription only).
- **REJECTED / UNVALIDATED:** fails a gate, or no independent check exists.

## Rediscovery vs novel (the critical fork)
- **Rediscovery** leans on Axis C (the published anchor) — the gold standard.
- **Novel** has no Axis C, so trust rests entirely on A+B+D; it is never silently promoted past SILVER. This honesty is the whole point: the pipeline must not manufacture confidence it can't justify.

## Component map (≈70% exists)
- Stages 1–2: `feasible_cells`, `construct_resonant_cycler`, `optimise_cell_idealized` / `optimise_cell_ephemeris`, 3D ephemeris + 3D Tisserand — **built**.
- Axis A: lamberthub V1 crosscheck + Kepler re‑propagation — **built**; add construction‑vs‑optimiser agreement.
- Axis B: three fidelities exist — **add** the persistence gate.
- Axes C/D + provenance/fidelity tags — **the data‑validation‑hardening plan**.
- Stage 5–6: ledger + `canonical_signature`/`match` — **built**; add tier ratchet + human gate.

## Orchestration
A workflow that fans candidates out and runs each through the gauntlet in parallel (multi‑agent adversarial‑verify), with **loop‑until‑dry** on the novelty frontier and **per‑finding adversarial panels** (N independent verifiers; majority‑refute kills it). This institutionalises the two saves of this session: the orbit‑15 source‑typo catch (checksum/re‑extraction) and the verify‑agent's missed xfail‑flip (independent re‑verification).

## Non‑goals
- High‑fidelity GMAT integration (V4) is a stretch sign‑off step, not in‑pipeline.
- Low‑thrust cyclers (deferred #37) and planet‑centric moon‑tours (#76) are separate scopes.
- Automated "discovery" *claims* without the human gate — explicitly forbidden.
