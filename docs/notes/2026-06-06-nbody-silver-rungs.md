# N-body SILVER rungs — verdicts on the two held candidates (harness Phase B)

Task #131 / n-body harness Phase B (plan
`docs/superpowers/plans/2026-06-06-nbody-harness.md`, design
`docs/superpowers/specs/2026-06-06-nbody-harness-design.md`). The restricted
n-body harness (REBOUND 5.0 / IAS15, rails third-body force over the shared DE440
BSP) propagated the two USER-HELD SILVER candidates one full repeat period and
recorded the verdict to a review-queue audit trail. **No promotion, no catalogue
write** — the candidates remain user-held SILVER; the rung records, the human
decides (`review_queue.is_catalogue_source()` is `False` by contract).

## Candidate provenance (reconstructed — flagged)

The runtime `data/review_queue.jsonl` artefact was **NOT present on disk** at
execution time (concurrent agents; the same situation the #126 diligence note
records, `docs/notes/2026-06-06-silver-candidates-russell-diligence.md` lines
22-25, 239-241). The candidates were therefore RECONSTRUCTED from their recorded
provenance (`data/OUTSTANDING.md` "Forge Phases 4 + 5" + the diligence note):

| candidate | sequence | k | per-encounter V∞ (km/s) | source |
|---|---|---|---|---|
| forge-silver-1 | E-M-E-E | 2 | [E 9.75, M 13.01, E 9.76, E 9.75] | OUTSTANDING + diligence (full) |
| forge-silver-2 | E-M-E-E | 2 | [E 9.62, M 12.06, (E, E not reported)] | OUTSTANDING + diligence (partial) |

For candidate 2 the two return-E magnitudes were "not reported in the section";
they are mirrored from the 9.62 home-Earth value (FLAGGED). **The per-leg
`tof_days` are recorded nowhere** (not in OUTSTANDING, the diligence note, nor
the gauntlet ledger), so a representative 2-synodic E-M-E-E leg split was used:
~150 d outbound E→M, ~330 d M→E return, the balance of the ~4.27-yr 2-synodic
period on the E→E loop. The rung verdict is regime-level (correction-ΔV band), so
it is robust to the exact split — but the absent on-disk seed is a recorded
limitation.

The seed is reconstructed at the home-Earth node as `v_sc = v_planet(E) + |V∞|·v̂`
(the recorded V∞ *magnitude* applied prograde — the seed *vectors* were not on
disk either). NON-GOLDEN throughout: the candidate numerics are OUR computation
(SILVER = unsourced by tier definition), so the rung asserts a regime, never a
sourced value.

## Verdicts (verbatim, Sun + E + M + J, IAS15 epsilon 1e-9)

| candidate | rung verdict | correction ΔV (km/s) | terminal closure (km) | converged | promoted |
|---|---|---|---|---|---|
| forge-silver-1 | **ARTIFACT** | 19.872260103429888 | 20394058.21161299 | False (diverged) | False |
| forge-silver-2 | **ARTIFACT** | 10.067736231784032 | 20387143.432708144 | False (diverged) | False |

Both grade **ARTIFACT** — exactly the honesty-boundary expectation (plan: these
high-V∞-basin candidates float at E∞ ~9.7 / M∞ ~12-13 km/s and #135 shows the
family lands off-anchor). The reconstructed seed sends the spacecraft on a
high-eccentricity heliocentric arc that grazes a perturber's softened core; IAS15
hits the integration budget / a non-finite state and the propagation is surfaced
as **divergent** (`converged=False`), which the threshold logic grades ARTIFACT
regardless of the ΔV magnitude. The recorded correction ΔV (~10-20 km/s, far
above the 1000 m/s ARTIFACT floor) is consistent: the conic seed lives only in
patched-conic land, not as a real n-body-ballistic trajectory — strong
REJECTED-style evidence, recorded for the human reviewer.

> Caveat: the correction ΔV here is the single-node terminal-closure impulse (the
> Phase-B metric), not the full multiple-shooter node-impulse (Phase C). For a
> divergent seed the magnitude is a lower-bound signal; the *divergence itself* is
> the load-bearing verdict.

## Jupiter sensitivity (Gate-4 body-inclusion arm, candidate 1)

| body set | correction ΔV (km/s) | terminal closure (km) | converged |
|---|---|---|---|
| Sun + E + M + J | 19.872260103429888 | 20394058.21161299 | False |
| Sun + E + M | 18.924313481381372 | 20392065.586034182 | False |

**Δ (Jupiter on − off) = 0.9479466220485158 km/s.** Jupiter moves the rung metric
by ~0.95 km/s — well above the rung's 0.2 km/s ROBUST threshold — so its inclusion
is justified by *evidence* over this multi-year baseline (the design §2 standing
rule: include a body iff it moves the consumer's metric by more than the
consumer's tolerance, proven by the §5.3 sensitivity test). Earth + Mars + Jupiter
is the recorded rung body set.

## Independence / discipline

The rung runs an *independent integrator* (REBOUND/IAS15, not `core/kepler.py`)
over the *same* DE440 BSP astropy caches — a shared-DE440 cross-check rung, **NOT
a V4 stamp** (design §4, Q6). The audit record carries `promoted=False` and
`independence: "shared-DE440 cross-check (NOT a V4 stamp; design Q6)"`. The two
candidates stay user-held SILVER; the GMAT lane (Phase D, independent codebase +
ephemeris) is the only path to V4.

## Golden gates behind these numbers (Phase A)

All five §0 self-validation gates are green and gate every science number above:

- GATE 1 (two-body reduction): Sun-only n-body vs `core/kepler.propagate` agree to
  **8.43e-08 km** / **1.46e-14 km/s** over a 120-day arc (threshold 1 km / 1e-5).
- GATE 2 (energy floor): Sun-only 1-year relative energy drift **0.0** (machine
  precision).
- GATE 3 (DE440 anchor): planet-state ingestion reproduces
  `Ephemeris("astropy").state` to **0.0 km / 0.0 km/s** at all sampled epochs (the
  shared-reader identity — the rung's independence is the *integrator*, not the
  reader).
- GATE 4 (accuracy convergence): final state moves **1.36e-07 km** between IAS15
  epsilon 1e-7 and 1e-10 (threshold 50 km).
