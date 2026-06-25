# #444 — C21 stable 3D candidate: gauntlet PASS, novelty CLOSED (known class)

**Date:** 2026-06-25. Disposition of the #438 headline candidate — the stable
out-of-plane 3D extension of the Braik-Ross **C21 (2,1)** Earth-Moon cycler.

**Candidate (exact):** `braik_ross_system()` (mu = 0.012150584270572),
state0 = `[0.7440212218499672, 0, -0.2057098355650995, 0, 0.35368280201143637, 0]`,
period T = 18.167169790651315 TU, Jacobi C = 3.025791614769996, winding
(k1,k2,k_z) = (2,0,10), Floquet tag **stable**.

## What was run (all gates the candidate genuinely PASSES)

| Gate | Result | Margin |
|---|---|---|
| **Floquet stability** (independent re-derivation) | all 6 |ev| = 1.0000 | n/a — on the unit circle |
| **V1** same-model (`run_v1_3d`, full-asymmetric 6D closure + independent Radau) | **PASS**, n_iter=1, `degenerate_planar=False` | indep. closure 3.2e-12 km/s vs 1e-3 floor (9 orders) |
| **V2** long-span (`run_v2_3d`, 3- and 6-cycle bounded drift) | **PASS** | 6-cycle drift 1.7 mm vs 50 000 km floor (7 orders); LINEAR growth |
| **JPL** full-EM-family novelty sweep (halo/vertical/axial + 7 resonant branches) | **no match** | closest 13/23 resonant overlap period but planar + C≤2.97 |

V1's n_iter=1 confirms the IC was already a genuine full 6D periodic orbit (not
merely symmetric at the perpendicular crossing). V2's *linear* (non-exponential)
drift is the dynamical fingerprint of the all-|ev|=1 Floquet verdict. These are
real, correct results.

## The binding gate — publication novelty — returns KNOWN CLASS

A focused literature check (3 WebSearch queries) is **decisive**:

> Spatial (3D) resonant periodic orbits **bifurcate from the vertical critical
> orbits of the planar resonant family**, and **stable** spatial members exist
> for the **2/1** MMR (and 3/2, 5/2, 3/1, 4/1, 5/1), for prograde and retrograde
> motion across broad inclination ranges.

Established by **Antoniadou & Voyatzis (2013, CMDA, DOI 10.1007/s10569-012-9457-4;
2014)**, **Voyatzis, Antoniadou & Tsiganis (2014)**, and **Antoniadou & Libert
(2019, MNRAS 483(3):2923, DOI 10.1093/mnras/sty3195)** — the last is *already a
`known_corpus_3d` anchor*. Supporting: arXiv 1811.09442, 1211.0964, 1805.00288.

The candidate is the **Earth-Moon-µ realization** of this published class — a
correctly-computed stable spatial 2:1-resonant periodic orbit — **not a
discovery**. Its JPL absence is a catalogue-**coverage** gap (the JPL DB does not
enumerate every spatial-resonant family), not undiscovered structure.

## Verdict

**NOT NOVELTY-CLAIMABLE. No catalogue writeback.** The same-model gauntlet,
stability, and JPL non-match all hold but are necessary-not-sufficient; the
publication-novelty gate — the binding one for a discovery — returns
published-class. This **corrects the over-optimistic "strongest novel candidate"
framing** in the #438 verdict and the earlier #444 registry entries.

Per the discipline (hard-scrutinise "it closed!"; a clean novelty negative is a
success): the thread is closed honestly.

## What stands / is reusable

- The V1+V2 3D gauntlet was exercised on a real candidate end-to-end — concrete
  progress toward **#306** (3D V0-V5 gauntlet). The REBOUND IAS15 + `frames.py`
  rotating↔inertial infrastructure is mapped and ready for a periodic-orbit V3
  (genuinely independent integrator architecture) when #306 needs it.
- The candidate is a good **validation exemplar**: its IC + T + C + Floquet
  signature can be cross-checked against Antoniadou & Libert 2019 — but any
  golden's EXPECTED side must be **sourced from the paper**, never our computed
  IC (golden-sourced-only discipline).

## Registry

Four method-versioned entries record the full thread (scrutiny → JPL fullsweep →
V1V2-gauntlet → LITNOVELTY), the last superseding the framing of the first three.

## Disposition

- **#444: closed** — honest negative on novelty; gauntlet engineering retained.
- The genuine novel-cycler frontier is the **isolated-family / from-scratch**
  search (#441 Phase 2, #442 Phase 3 Earth-Moon-cycler seeding), not 3D
  extensions of known planar roots.
