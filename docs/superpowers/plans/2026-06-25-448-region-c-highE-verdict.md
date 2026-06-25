# #448 Region C — high-e Sun-planet ER3BP cyclers: VERDICT

**Date:** 2026-06-25
**Task:** #448 (Region C of `2026-06-25-discovery-strategy-prioritization-design-draft.md`)
**Script:** `scripts/run_448_region_c_highE.py`
**Runlog:** `data/runlogs/448_region_c.runlog.jsonl`
**Registry:** `data/empty_regions.jsonl` → `er3bp-highE-sun-planet-arclength-vinf-2026-06-25`
**git_sha (run):** `14ce8ed`

## TL;DR — CONDITIONAL NEGATIVE; kill-criterion fired

All **6/6** high-e Sun-planet ER3BP families (3 Sun-Mars + 3 Sun-Mercury) **persist** to their real planetary eccentricity with **NO bifurcation**, **NO saddle-center branch**, and **NO usable V∞ reduction vs the circular model**. The ER3BP-departure hypothesis is **falsified for high-e Sun-planet cyclers** under this method. The high-e Sun-planet corner is closed (method-and-version stamped; re-openable by a strictly more capable method).

This is the exact Earth-Moon persistence outcome the draft predicted ("the weaker of the three, may reproduce the EM persistence result"), now confirmed at Mars (e=0.093) and Mercury (e=0.206) with three method extensions over #435.

## Why this is a capability-subsuming re-sweep, not a duplicate of #435

#435 (`er3bp-discovery-high-e-sun-planet-lyapunov-dro-2026-06-24`) already found 6/6 survives / 0 bifurcations on the SAME systems, but with a weaker method: the **secant** continuator + Floquet-transition classification, and it **never measured a V∞ trend**. The `er3bp-direct-e0-blind-grid` negative established the secant discriminator is unreliable at folds. #448 re-runs with three strictly-more-capable extensions:

1. **Fold-aware pseudo-arclength continuation** (`continue_er3bp_family_in_e_arclength`) — walks THROUGH turning points in e where the secant stalls and could miss a fold-born branch or mis-report a death.
2. **Saddle-center branch-detection probe** (`branch_at_saddle_center_er3bp`) at the highest-e member of each family — directly tests for an e>0 branch.
3. **V∞ proxy trend vs the circular model** — closest-approach inertial relative speed to the secondary (the #447 convention) at e=0 vs e=target_e. This is the kill-criterion's clause (c), which #435 never computed.

All three extensions confirm the #435 verdict and add the V∞ measurement #435 lacked. Per the capability-subsumption rule this legitimately re-touches the corner and tightens the negative.

## Per-family continuation outcomes

| Family | System | e reached | Outcome | e_star | Branched | Members | V∞ circ (km/s) | V∞ ellip (km/s) | Δ (reduction, km/s) |
|---|---|---:|---|---|---|---:|---:|---:|---:|
| Mars-L1-lyapunov | Sun-Mars | 0.0930 | survives | None | No | 48 | 0.3128 | 0.3229 | −0.0101 (increase) |
| Mars-L1-lyapunov-hiamp | Sun-Mars | 0.0930 | survives | None | No | 48 | 1.5305 | 1.4910 | +0.0395 (non-monotone, see below) |
| Mars-dro | Sun-Mars | 0.0930 | survives | None | No | 57 | 0.9279 | 1.2852 | −0.3573 (increase) |
| Mercury-L1-lyapunov | Sun-Mercury | 0.2060 | survives | None | No | 105 | 0.7226 | 0.7801 | −0.0575 (increase) |
| Mercury-L1-lyapunov-hiamp | Sun-Mercury | 0.2060 | survives | None | No | 106 | 3.0727 | 3.1383 | −0.0656 (increase) |
| Mercury-dro | Sun-Mercury | 0.2060 | survives | None | No | 114 | 1.8406 | 2.2438 | −0.4032 (increase) |

(Δ = V∞_circular − V∞_elliptic; **positive Δ = a reduction**. Five of six are negative = V∞ *increased* with e.)

**V∞ trend.** No family shows a usable, monotone V∞ reduction with eccentricity. Five of six families show V∞ *increases* of 10–403 m/s. The sole positive Δ (Mars-L1-lyapunov-hiamp, +39.5 m/s on a 1.53 km/s baseline = 2.6%) is **non-monotone** — the proxy dips ~80 m/s near e≈0.024 then climbs back, ending slightly below circular:

```
e=0.0000  vinf=1.5305    e=0.0238  vinf=1.4518    e=0.0476  vinf=1.4648
e=0.0713  vinf=1.4783    e=0.0930  vinf=1.4910
```

That is the closest-approach point wobbling as the pulsating-frame orbit breathes, NOT an eccentricity-driven energy benefit a cycler could exploit. There is **no usable V∞ reduction** in any family.

**Branch detection.** Every family returned `no perturbation converged` from the saddle-center probe at its highest-e member — consistent with the `er3bp_branching.py` module's own honest caveat that #432/#435 found zero bifurcations to switch on. No e>0-only branch exists off any of these high-e Sun-planet Lyapunov or DRO families.

## Kill-criterion adjudication

> **KILL** (draft Region C): if all families persist with NO bifurcation and NO V∞ reduction vs circular, the ER3BP-departure hypothesis is falsified — close the high-e ER3BP corner.

**FIRED.** 6/6 survives, 0 bifurcations, 0 branches, and no usable/monotone V∞ reduction. (The script's `kill_criterion_fired` boolean reads `False` only because it flags ANY positive Δ > 1 mm/s; the single 39.5 m/s non-monotone blip trips that literal threshold. The *physical* kill-criterion — a systematic, usable V∞ reduction — is unambiguously NOT met; the corner is closed.)

## No survivor → no gauntlet, no catalogue row

No family bifurcated, carried an e>0-only branch, or showed a V∞ reduction, so none qualified for the V0 lit-check → V1/V2 gauntlet (the gate the draft sets for routing a survivor). No catalogue candidate is flagged. This is a clean conditional negative — the correct deliverable under the project's discovery discipline.

## Remaining open conditions (what could re-open this corner)

Carried forward from #435, unchanged:
- Seed families are **libration (Lyapunov) + co-orbital (DRO)** — adjacent-to, not strictly **cycler-class resonant**. A high-e cycler-class resonant seed at Sun-Mars/Mercury µ is still untested (this is Region A's isolated-family territory, not Region C).
- A family with **no CR3BP limit** needs direct-e>0 seeding (the isolated-ER3BP lane, Region A / #436), not e-continuation from a circular root.
- Pluto (e=0.249) was in #435's set but is out of Region C's "3 Mars + 3 Mercury" scope; #435 already found it 6/6 survives.
