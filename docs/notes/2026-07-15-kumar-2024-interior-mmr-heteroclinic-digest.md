# Kumar-Rawat-Rosengren-Ross 2024 — Interior MMR + heteroclinic connections, Earth-Moon (IAC-24-C1.9.5)

**Date filed/digested:** 2026-07-15 (light digest); **full mining pass:** 2026-07-15 (#597/#15).
**Trigger:** same `ross.aoe.vt.edu/papers/` review as the AAS-25-569 companion. Confirmed genuinely
new via title-page read; full 15-page body read page-by-page for this pass.

**Citation:**
> Kumar, B., Rawat, A., Rosengren, A. J., and Ross, S. D. (2024). "Investigation of Interior Mean
> Motion Resonances and Heteroclinic Connections in the Earth-Moon System." 75th International
> Astronautical Congress (IAC), Milan, Italy, 14-18 October 2024. IAC-24-C1.9.5.
> No DOI found (IAC conference paper, not indexed in CrossRef as of this filing).

Registered as `CorpusAnchor` `kumar-2024-interior-mmr-iac` in
`src/cyclerfinder/search/literature_check.py` (provenance `verified-against-source`).

## Model and method

Planar CR3BP, Earth-Moon, `mu = 1.2150584270572e-2`. Hamiltonian form in synodic Delaunay
coordinates `(L, G, l, g)` with `L = sqrt((1-mu) a)`, `G = L sqrt(1-e^2)`, `g` = argument of periapse
w.r.t. the rotating-frame +x axis. At `mu=0` these are action-angle (Kepler-in-rotating-frame), so
perturbation-theory intuition (pendulum separatrices) carries over qualitatively to `mu>0`.

**Poincaré section**: perigee (Earth-relative osculating true/mean anomaly = 0), NOT the more common
fixed-`x`/fixed-`y` sections — chosen because mean anomaly is (almost always) monotonically
increasing, so the section rarely develops tangencies with the flow (fixed-x/y sections do, causing
manifold-plot discontinuities in prior literature). Perigee-crossing detection uses a dot-product
sign-change test (position-velocity orthogonality) rather than a direct anomaly threshold, with an
`l~0 or 2*pi` confirmation guard against false positives near the Moon (where `l` can transiently
decrease).

**Periodic-orbit continuation**: seed from the Earth-Kepler two-body problem at
`a=(n/m)^(2/3)` (m:n MMR), continue in `mu` via perpendicular x-axis crossings (symmetric periodic
orbits only) to `mu=1.2150584270572e-2`, then continue through the family using the x-intercept (or
Jacobi constant `C` at fold points) as the continuation parameter.

**Manifolds**: parameterization method (Haro et al. 2016) extended to the m-iteration (period-`m`
Poincaré map) case — an `m:n` resonant orbit crosses the perigee section `m` times per period, so its
points are period-`m`, not period-1, fixed points of the map.

## Key quantitative/topological findings (4:1, 3:1, 2:1 interior)

- **4:1**: stable+unstable, prograde+retrograde ALL belong to a single continuous family (joined by
  fold bifurcations, where Floquet multipliers pass through 1). Apoapsis is **always below the
  Moon's orbit** (`x < 1-mu`) at every eccentricity tested, even the highest-e members. Exists over a
  wide `C` range (checked down to `C=2.85`, where periapsis already grazes Earth's surface).
- **3:1**: stable and unstable orbits belong to **two separate families** (not one, unlike 4:1). Also
  never reaches the Moon's orbit at apoapsis. Exists over a comparably wide `C` range — the lowest-`C`
  3:1 unstable orbit found is at `C=3.45`ish territory extending down through the whole 3.00-3.45+
  band tested.
- **2:1**: qualitatively different and richer. A single continuous family links non-resonant
  near-circular orbits -> stable prograde 2:1 -> (through the Earth singularity) unstable retrograde
  -> stable retrograde -> ... -> non-resonant again (matches Broucke 1968 "Family BD", not previously
  identified as 2:1-resonant). Has members with **apoapsis above the Moon's orbit** as well as below.
  Crucially, **2:1 UNSTABLE prograde orbits only exist for `C <~ 3.15-3.16`** (higher eccentricity/
  lower energy than 3:1 or 4:1) — there is a large `C` gap `(3.16, 3.45]`-ish where 3:1 (and 4:1)
  unstable orbits exist but 2:1 does not.

## Heteroclinic / resonance-overlap results (the paper's main payload)

- For `C` in the gap before 2:1 unstable orbits appear (single-resonance regime), 3:1 manifolds
  closely resemble an ideal-pendulum separatrix — consistent with standard Hamiltonian perturbation
  theory (no overlap).
- **As soon as 2:1 unstable orbits appear (`C <= 3.15`), heteroclinic connections between 3:1 and 2:1
  appear IMMEDIATELY** — no low-energy "coexist without overlap" window of the kind seen in typical
  low-`mu` systems (e.g. Jupiter-Europa). The paper attributes this to the Earth-Moon system's large
  `mu ~ 0.0122` breaking the small-`mu` perturbation-theory assumption. Confirmed present for
  `C = 3.15, 3.10, 3.05, 3.00` (all four explicitly checked).
- For `C <= 3.10`, some 2:1-manifold points on the perigee section **escape toward `L = sqrt((1-mu)a)
  = 1` (the Moon's semimajor axis)** — i.e. natural (zero-cost) transport from a 2:1 unstable orbit to
  the Moon vicinity, without needing the 3:1 intermediary.
- **4:1 has NO heteroclinic connection to 3:1 at ANY tested energy** (`C = 3.15, 3.00, 2.85` all
  checked; the `C=2.85` orbit already has periapsis below Earth's surface, i.e. as extreme as the
  family gets). The 4:1 stable/unstable manifolds nearly coincide (like an ideal-pendulum
  stable=unstable degenerate separatrix) rather than crossing.
- A **rotational invariant circle (RIC)** — a persisting (deformed) 1D KAM torus from the Keplerian
  limit — is identified in the `(g, L)` Poincaré section at `C=3.10`, sitting at `L~0.78` between the
  4:1 and 3:1 resonance islands (`L=0.625` and `L=0.68` respectively). This RIC is a literal
  topological transport barrier: it explains WHY 4:1 never connects to 3:1 (similar barrier suspected,
  less clearly, at `C=3.05`; unclear whether one persists at `C<=3.00`).
- **Conclusion (paper's own words): "we have identified the 3:1 mean motion resonance as truly being
  a 'gateway to the Moon' for lower-energy mission design"** — reachable from 2:1 (heteroclinic) and
  from 2:1's own Moon-ward manifold escape, while 4:1 is dynamically isolated by the RIC barrier.
- Companion/follow-on work (their ref [22], `= 2026-07-15-rawat-2024-resonance-widths-chaotic-zones-digest.md`
  in this corpus) is explicitly cited as extending this with a semi-analytical (perturbed-Kepler)
  comparison, and is reported there to UNDERESTIMATE the 2:1/3:1 resonance widths relative to this
  paper's full-CR3BP computation.

## Scope limits (stated by the authors)

Planar only. Real cislunar missions in these resonances (TESS 2:1, IBEX 3:1) actually orbit
**inclined**, not planar — the authors flag the spatial CRTBP extension as the most significant future
direction, noting it needs new tools (Poincaré map remains a 4D symplectic map even after reducing to
a fixed-energy submanifold in the spatial problem, vs. 2D in the planar case).

## Relevance to this project

- Directly on-topic for #314 (heteroclinic-cycle genome, currently validated only against the
  Sun-Jupiter-Oterma system)/#411 (cross-system SE<->EM closure)/#496/#503. None of this project's
  existing heteroclinic-search machinery (`genome/heteroclinic_cycle.py`,
  `genome/cross_system_cycle.py`) currently targets Earth-Moon interior-MMR-to-interior-MMR chains —
  this paper (plus its two AAS companions) is the concrete literature template for that if it is ever
  built.
- **Directly informs `src/cyclerfinder/search/resonance_network.py`'s (#267 Track-B tier 3) own
  documented data gap.** That module is built against Kumar-Rawat-Rosengren-Ross's 2026 ASR journal
  paper ("Cislunar Resonant Transport and Heteroclinic Pathways: From 3:1 to 2:1 to L1," the anchor
  already registered as `kumar-2025-arxiv-2509.12675`) and its docstring states the paper PDF is
  **not held in the local mirror**, so the exact common Jacobi constant, unstable-member periods, and
  metric form used by that journal paper are NOT independently sourceable — its reproduce-before-trust
  gate is `xfail`-marked for this reason. This IAC paper is an earlier conference-stage iteration of
  the SAME author group on the SAME 3:1/2:1/4:1 chain (companion to AAS-24-368 below), and does supply
  concrete, sourced numbers: `C<=3.15` for the 3:1<->2:1 heteroclinic window, the specific RIC-barrier
  explanation for why 4:1 is excluded from the chain (matching `resonance_network.py`'s implicit "3:1
  -> 2:1 -> L1, not 4:1" scope), and the perigee-Poincaré-map/parameterization-method methodology the
  module's own docstring describes. **This does not close `resonance_network.py`'s reproduce-before-
  trust gap** (it is still a different, earlier paper than the exact one cited there, and does not
  give that paper's specific table values) — but it is independent corroborating evidence for the
  module's core physical claims and a candidate partial substitute if the exact arXiv:2509.12675 PDF
  remains unobtainable. Flagged as a follow-on opportunity, not actioned in this pass (out of scope
  for a digest-only task).

**Status: full mining pass complete.** No catalogue impact, no writeback (this is a
methods/results reference, not a specific cycler with citable elements). `CorpusAnchor` registered.
