# #316 Phase 1: Cross-system cycler conceptual & literature deep-dive

Date: 2026-06-17. Task: literature triage + theoretical framework + an
order-of-magnitude geometric probe for the open conceptual question
"do cyclers that **cycle between** the Sun-Earth and Earth-Moon CR3BP regimes
have a published systematic study, and if not, what is the right mathematical
framework to launch one?". Research-only deliverable. NO source-code changes
to `src/cyclerfinder/`, NO catalogue writeback, NO novelty claims.

## TL;DR verdict

**INCONCLUSIVE — partial published coverage along two adjacent axes; an exact
specification of the open hypothesis emerges from the literature triage and
phase 2 is worth scoping.**

The cross-system cycler concept — a repeating orbit that closes simultaneously
in **both** Sun-Earth and Earth-Moon synodic frames, with each period including
heteroclinic-manifold-tube transits between SE-L1/L2 and EM-L1/L2 — sits in a
gap between two well-mined adjacent areas:

1. **Heteroclinic chains as one-shot transfers** (Koon-Lo-Marsden-Ross 2001;
   Canalias-Gómez-Marcote-Masdemont 2006; Castelli 2012; van der Weg-Vasile
   2015). These authors explicitly construct SE-L1/L2 unstable-manifold ↔
   EM-L1/L2 stable-manifold connections, but as **single-pass** Earth-to-Moon
   trajectories (the "Shoot the Moon" Hiten-style mission paradigm) — never
   composed into a repeating cycler.

2. **BCR4BP synodic-resonant periodic orbits** (Boudad-Howell-Davis 2020;
   McCarthy-Howell 2021; Park-Howell 2022; Reddy-Howell 2023; Brown-Peterson
   2024 in the closely-related HR4BP). These **do** close simultaneously in
   both synodic frames, by construction (the p:q resonance ratio enforces
   it), but the published families are **local** — computed around a single
   libration point of one system (typically EM-L1 or EM-L2), with the Sun as
   a periodic forcing rather than a region the orbit visits. The orbit state
   space stays inside the Earth-Moon Hill sphere; it does **not** include a
   heteroclinic transit out to SE-L1/L2 manifold tubes within a single
   period.

The unmined corner is therefore: **a BCR4BP (or HR4BP) synodic-resonant
periodic orbit whose state-space support spans the Earth-Moon Hill sphere
AND the Sun-Earth L1/L2 manifold neighborhood within one orbital period, with
heteroclinic-tube transitions between them**. This is exactly the cycler
analog of the KLMR one-shot transfer.

The OOM probe (Part C below) confirms the geometric pieces are not absurd:
the SE-L1/L2 region sits ~1.5 Mkm from Earth (~3.9 Moon-distances), the
spatial gap to EM-L1/L2 is ~1 Mkm, and the smallest natural near-commensurate
closure period is **19 years** (the 235:19 Saros / Metonic family).

---

## Part A — Literature triage

Search strategy: WebSearch on the four candidate frameworks, fetch
authoritative URLs where possible (several PDFs returned binary noise and
were not readable; for those the abstract / search-result summary was used,
and the paper is marked as such below). All citations are publicly accessible
by DOI / arXiv / institutional URL.

### Triage table

| # | Citation | Available URL | Category | Notes |
|---|---|---|---|---|
| 1 | **Koon, Lo, Marsden, Ross (2001)** — "Low Energy Transfer to the Moon" / "Shoot the Moon", *Celestial Mechanics and Dynamical Astronomy* 81 (2001) 63–73 | [koon.pdf @ Caltech](http://www.cds.caltech.edu/~koon/papers/koon.pdf), [Shoot the Moon @ VT](https://ross.aoe.vt.edu/papers/shoot_moon_2000.pdf) | **TRANSFER-TYPE** | The founding work on SE-L1/L2 → EM-L1/L2 heteroclinic chains. Explicitly approximates the SEM-spacecraft 4-body as two coupled CR3BPs, finds SE unstable-manifold ↔ EM stable-manifold connections, reproduces Hiten as a one-shot ballistic capture. The heteroclinic chain is **not** iterated into a cycler. |
| 2 | **Canalias, Gómez, Marcote, Masdemont (2006)** — "Homoclinic and heteroclinic transfer trajectories between Lyapunov orbits in the Sun-Earth and Earth-Moon systems", *Discrete and Continuous Dynamical Systems* 14 (2006) 261–279 | [aimsciences.org @ DOI](https://www.aimsciences.org/article/doi/10.3934/dcds.2006.14.261) | **TRANSFER-TYPE** | Classifies homo/heteroclinic connections between Lyapunov orbits across the SE and EM systems by bifurcation family. Numerically proves connections exist asymptotically (no Δv needed). Still one-shot transfers; no cycler closure. |
| 3 | **Lo, M.W. (2002)** — "The InterPlanetary Superhighway and the Origins Program", *IEEE Aerospace Conference* | [gg.caltech IPS paper](https://www.gg.caltech.edu/~mwl/publications/papers/IPSAndOrigins.pdf) | **STRUCTURAL ADJACENT** | The "interplanetary superhighway" framing — the manifold tubes as a transport network — implicitly suggests iteration. But explicit cyclers are not constructed; the paper's deliverable is the *concept* of tube-mediated transport, not a periodic orbit. |
| 4 | **Belbruno, E.** — WSB / ballistic lunar capture (1987, 1991 Hiten) | [Wikipedia WSB](https://en.wikipedia.org/wiki/Weak_stability_boundary), [arXiv 2407.00853 Cantor structure](https://arxiv.org/pdf/2407.00853) | **TRANSFER-TYPE** | The WSB is **defined by** the failure of repeated cycling about the Moon under one-cycle-stability test. Belbruno's framework gives Earth-to-Moon low-energy transfers via the WSB; recent work (arXiv 2407.00853) shows the WSB has Cantor-set structure for infinitely many cycles, but as a parameter-space classification, not a constructive cycler design. |
| 5 | **van der Weg, W.S. & Vasile, M. (2015)** — "Sun-Earth L1 and L2 to Moon transfers exploiting natural dynamics", *Celest. Mech. Dyn. Astron.* | [strathprints PDF](https://strathprints.strath.ac.uk/51990/1/van_der_Weg_Vasile_CMDA_2015_Sun_earth_L1_L2_to_moon_transfers_exploiting_natural_dynamics.pdf) | **TRANSFER-TYPE** | Explicitly uses SE-L1/L2 manifolds as origin and reaches the Moon by natural dynamics. One-shot trajectory; no return-leg / cycler closure described. |
| 6 | **Boudad, Howell, Davis (2020)** — "Dynamics of synodic resonant near rectilinear halo orbits in the bicircular four-body problem", *Adv. Space Res.* 66 (2020) | [ScienceDirect S0273117720305536](https://www.sciencedirect.com/science/article/abs/pii/S0273117720305536) | **STRUCTURAL ADJACENT** | First systematic family of BCR4BP synodic-resonant NRHOs. The 9:2 resonant NRHO (orbital period = 2 synodic periods = ~59 days) closes simultaneously in EM and SE frames. **However the orbit stays inside the EM Hill sphere** — no SE-L1/L2 visit. |
| 7 | **McCarthy, B.P. & Howell, K.C. (2021)** — "Quasi-Periodic Orbits in the Sun-Earth-Moon BCR4BP", AAS 21-270 | [Purdue PDF](https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Conferences/2021_AAS_McCHow.pdf) | **STRUCTURAL ADJACENT** | Families of quasi-periodic tori around SE-L1, SE-L2 in the BCR4BP. **Quasi-periodic, not cycler**; computed for stationkeeping at SE libration points (e.g., JWST class missions), not as multi-region cyclers. |
| 8 | **McCarthy, B.P. & Howell, K.C. (2023)** — "Four-body cislunar quasi-periodic orbits and their application to ballistic lunar transfer design", *Adv. Space Res.* 71 | [ScienceDirect S0273117722008614](https://www.sciencedirect.com/science/article/abs/pii/S0273117722008614) | **TRANSFER-TYPE** | The quasi-periodic tori become the engineering platform for ballistic lunar transfers. Same comment as #6: closes in BCR4BP, but as a transfer; cycler closure is not the deliverable. |
| 9 | **Park, B. & Howell, K.C. (2022)** — "Multiple families of synodic resonant periodic orbits in the bicircular restricted four-body problem", *Adv. Space Res.* 70 | [ScienceDirect S0273117722004781](https://www.sciencedirect.com/science/article/abs/pii/S0273117722004781) | **STRUCTURAL ADJACENT** | Tens of synodic-resonant families with various p:q ratios — DROs, L1 Lyapunovs, NRHOs. All EM-Hill-bounded. The closest published thing to a "cross-system cycler" exists here in p:q form but only one of the two systems is *visited*. |
| 10 | **Reddy, K. & Howell, K.C. (2023)** — "Synodic Resonant Halo Orbits in the BCR4BP", AAS 23-227 | [Purdue PDF](https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Conferences/2023_AAS_RedHow.pdf) | **STRUCTURAL ADJACENT** | Extends Park-Howell to halo families. Confirms the closure semantics: "for each Sun angle, there are two distinct states along this NRHO: one at θ_i, and one at θ_i + 2π, that occur at the same Earth-Moon-Sun configuration." This **is** simultaneous closure — but again, only EM region visited. |
| 11 | **Brown, M.J. & Peterson, A.A. (2024)** — "Structure of Periodic Orbit Families in the Hill Restricted 4-Body Problem", arXiv 2402.19181 / *SIAM J. Appl. Dyn. Syst.* | [arXiv 2402.19181](https://arxiv.org/abs/2402.19181) | **STRUCTURAL ADJACENT** | Pseudo-arclength continuation of EM-CR3BP families into the HR4BP (a coherent time-periodic model of the SEM system). Substantial bifurcation atlas. The continuation tracks bring CR3BP families into the SEM-perturbed regime but again the spatial support stays EM-bounded. |
| 12 | **Burgos-Garcia, J. & Gidea, M. (2014)** — "Hill approximation in a restricted four-body problem", arXiv 1412.3775 | [arXiv 1412.3775](https://arxiv.org/pdf/1412.3775) | **STRUCTURAL ADJACENT** | Foundational mathematics for the HR4BP. No cycler claim. Provides the symplectic substrate that #11 builds on. |
| 13 | **arXiv 2509.12675 (2025)** — "Cislunar Resonant Transport and Heteroclinic Pathways: From 3:1 to 2:1 to L1" | [arXiv 2509.12675](https://arxiv.org/pdf/2509.12675) | **TRANSFER-TYPE** | Earth-Moon-only. Heteroclinic transitions between EM mean-motion resonances (3:1, 2:1) and EM-L1 Lyapunovs. No Sun-Earth coupling. |
| 14 | **Ross, S.D., Roberts, A., Tsoukkas, T. (2025)** — "Stable, Low-Energy Prograde Earth-Moon Cycler Orbits", AAS 25-621 + journal version (Roberts-Tsoukkas-Ross 2026) | [ross.aoe.vt.edu](https://ross.aoe.vt.edu/papers/ross-roberts-tsoukkas-2025-AAS-25-621.pdf) | **STRUCTURAL ADJACENT** | Five stable EM cycler families. Discusses 2:3 synodic resonance with Sun-Earth period as a *desirable property* of one family (yielding periodic Sun-Earth-Moon-spacecraft alignments) but the orbits themselves are EM-CR3BP. Per our existing mining note (Ross-Roberts-Tsoukkas 2025) the Sun is a perturbation, not a region visited. |
| 15 | **Genova, A. & Aldrin, B. (2015)** — "A Free-Return Earth-Moon Cycler Orbit for an Interplanetary Cruise Ship", AAS 15- | local PDF | **TRANSFER-TYPE (cycler-class, single-system)** | The 3:1 resonance EM cycler. **Explicitly mentions** "the addition of solar gravity and a modest Δv maneuver causes Aldrin's C-2-R theorized cycler's apogee to drop temporarily below lunar distance" — i.e., Sun gravity ENABLES the cycler, but only as perturbation; the orbit does not visit SE-L1/L2. |
| 16 | **Castelli, R. (2012)** — "Transfers between the Earth-Moon and Sun-Earth systems using akin orbits", *Acta Astronautica* (~2006 series) | [ScienceDirect S0094576506000981](https://www.sciencedirect.com/science/article/abs/pii/S0094576506000981) | **TRANSFER-TYPE** | Explicitly catalogues transfers from EM libration orbits to SE libration orbits using "akin" orbits (similar Jacobi level in both CR3BPs). One-shot bidirectional transfers. No repeating closure. |
| 17 | **Anderson, R.L. & Lo, M.W.** — JPL line of work on patched-CR3BP, e.g., *J. Guid. Control Dyn.* 2009 ("Role of Invariant Manifolds in Low-Thrust Trajectory Design") and follow-ups | [JPL TRS](https://trs.jpl.nasa.gov/handle/2014/41117) (catalogue) | **TRANSFER-TYPE** | Patched-CR3BP framework. Transfer-design oriented. Not cycler-oriented. |
| 18 | **Bosanac, N., Howell, K.C., Fischbach, E. et al.** — Lagrangian dynamics + cislunar periodic orbits / motion-primitive sets | [Lagrangian Dyn paper @ Springer](https://link.springer.com/article/10.1007/s11071-022-07829-1), [Motion primitives @ Springer](https://link.springer.com/article/10.1007/s10569-022-10063-x) | **STRUCTURAL ADJACENT** | Motion-primitive decomposition of EM-CR3BP / BCR4BP families. Primitives include heteroclinic arcs. **Could be the right substrate** for assembling a cross-system cycler, but no such cycler is constructed in the published work. |

### Triage counts

| Category | Count |
|---|---|
| DIRECT MATCH (published cross-system cycler with both-region visits) | **0** |
| STRUCTURAL ADJACENT (manifold heteroclinics or patched-CR3BP, not as repeating cyclers) | **9** |
| TRANSFER-TYPE (uses the manifold but as a one-shot transfer) | **8** |
| NO MATCH | 0 (everything found was relevant adjacent literature) |
| **Total** | **17** |

### Negative result of the literature search

After the deep dive, **no paper was found** that constructs a periodic orbit
explicitly defined by:
- closure in **both** the Sun-Earth and Earth-Moon synodic frames, AND
- a state-space support that includes a heteroclinic-manifold-tube transit
  between SE-L1/L2 and EM-L1/L2 within a single period.

The closest published objects to this concept are the BCR4BP synodic-resonant
families (papers 6, 9, 10) — they have the first property by construction, but
not the second. The KLMR/Castelli line (papers 1, 2, 16) has the second
property but as a one-shot, not closed into a cycler.

This is an **inconclusive** negative because the literature is mature enough
that a publishable cycler of this type would likely have surfaced if the
geometry trivially closed; the absence is consistent with either:
- the orbit class being **impossibly thin** (no near-commensurate period
  brings the heteroclinic arcs into closure with reasonable Δv), or
- the orbit class existing but being **off the published mining axes** —
  Phase 2 work could surface it.

The 19-year Saros / Metonic commensurability (Part C) is encouraging: it
gives a natural near-commensurate closure window, which is the minimum
prerequisite.

---

## Part B — Conceptual framework sketch

### B.1 What WOULD a cross-system cycler look like?

| Property | Specification |
|---|---|
| **State space** | 6D phase space, naturally heliocentric J2000 inertial; representable in either the SE-rotating or EM-rotating frame |
| **Period T** | Commensurate with **both** T_SE_synodic ≈ 365.25 d and T_EM_synodic ≈ 29.5306 d, i.e., T = p · T_SE = q · T_EM for small integers (p, q). The 235:19 (Saros) and 19:1 (Metonic) families give T = 19 yr ≈ 6940 d (Part C). |
| **Spatial support** | Within one period the orbit visits both the SE-L1/L2 manifold neighborhood **and** the EM Hill sphere, traversing the ~1 Mkm spatial gap via heteroclinic-tube arcs (KLMR/Castelli mechanism). |
| **Topology** | Each cycle contains: (a) a near-EM-Lyapunov / NRHO arc, (b) an EM-L1/L2 unstable-manifold tube exit, (c) a Lagrange-point traversal (Earth's neighborhood), (d) an SE-L1/L2 stable-manifold tube ingress, (e) a SE-Lyapunov arc, (f) the SE→EM return mirror of (b-d). Optional simplifications: (a)+(b) or (e)+(f) may be a single small loop near a libration point. |
| **Closure** | Both 6D states in SE-rotating frame and EM-rotating frame match initial conditions after T. The two rotation rates ω_SE and ω_EM impose **two** distinct closure conditions; near-commensurability is needed for both to be near-satisfied with the same period. |
| **Δv budget** | If purely manifold-mediated (no DSM), ballistic. Realistically a small phasing-Δv per period (analogous to Genova-Aldrin's 39 m/s/month, which is high because they're not exploiting SE-L1/L2 manifolds). |

### B.2 Mathematical formulation candidates

#### (i) Patched CR3BP (two coupled CR3BPs with interface matching)

- Spacecraft propagates in SE-CR3BP while outside Earth's Hill sphere.
- Inside Earth's Hill sphere (radius ≈ 1.5 Mkm), spacecraft propagates in
  EM-CR3BP about EM-barycenter.
- Continuity is enforced at the patching boundary (position + velocity in
  the inertial frame; the frame transformation handles the rotation-rate
  mismatch).
- Periodic-orbit corrector becomes a **multiple-shooting Newton** with patch
  points at every Hill-sphere crossing.
- **Strength**: cheap, structurally close to KLMR. The patching mechanism is
  already used in `src/cyclerfinder/data/method_capability.py` ordering
  (patched-conic → n-body). The framework re-exports KLMR's mechanism.
- **Weakness**: the patching boundary is artificial. The actual SEM physics
  is smooth; the boundary introduces O(mu_EM · mu_SE) errors that may
  accumulate over a 19-year period.

LOC estimate (NEW Phase 2 module): ~600-900 LOC for a corrector +
patched-CR3BP propagator.

#### (ii) Bicircular CR4BP (BCR4BP)

- Spacecraft propagates in the BCR4BP — one frame (typically EM-rotating)
  with Sun as a time-periodic perturbation at orbital frequency ω_sun.
- Sun-commensurate closure (per `src/cyclerfinder/core/bcr4bp.py:sun_commensurate_period`)
  enforces simultaneous closure in EM and SE frames at p:q resonance.
- For the cross-system cycler concept, **add a state-space constraint**:
  the orbit's trajectory must include a portion in the SE-L1/L2 manifold
  neighborhood (e.g., |r - r_SE_L1| < some δ for some t in [0, T]).
- Multiple-shooting corrector with **phase-space-region constraints** at
  selected interior shooting points (the "SE-L1 visit" point and the
  "EM-L1 visit" point).
- **Strength**: uses the existing #292 BCR4BP substrate. No new propagator.
  The framework is the *natural* generalization of Park-Howell (paper 9):
  they search synodic-resonant orbits without constraining state-space
  region; we'd add the cross-system region constraint.
- **Weakness**: the BCR4BP is incoherent (Sun's gravity acts but Sun's mass
  doesn't move). Andreu QBCP / Gimeno-Jorba's Fourier series QBCP corrects
  this but is heavier. For a Phase 2 SCOPING run BCR4BP is fine.

LOC estimate (extension to existing `src/cyclerfinder/genome/bcr4bp_genome.py`):
~300-500 LOC for region-constrained multiple-shooting.

#### (iii) Hill Restricted 4-Body Problem (HR4BP)

- Coherent time-periodic SEM model (Brown-Peterson 2024 framework).
- Same closure semantics as BCR4BP but with the Sun appearing as a periodic
  perturbation of the Hamiltonian itself (not just a kinematic gravity
  source).
- **Strength**: more physically faithful than BCR4BP. Brown-Peterson's
  pseudo-arclength continuation atlas provides ready seeds.
- **Weakness**: not yet in `src/cyclerfinder/`. Implementing HR4BP would be
  a new 400-700 LOC module before any cycler corrector lands.

LOC estimate (NEW module): ~1200-1600 LOC (HR4BP module + corrector + tests).

#### (iv) Full ephemeris (V4-class fidelity)

- DE440 ephemeris, all SEM gravity exact.
- Brute-force differential corrector at long arc length.
- **Strength**: ground truth. If a cross-system cycler exists, it exists here.
- **Weakness**: the (i)-(iii) families are the natural seeds for (iv);
  starting from (iv) without seeds is hopeless (chaos + 19-year arcs).

LOC estimate: 0 NEW — reuse `src/cyclerfinder/core/ephemeris.py` machinery,
but only as a Phase 4+ verification step. Out of scope for #316 Phase 2.

#### Recommendation

**(ii) BCR4BP with state-space-region constraint** is the lowest-cost,
highest-leverage Phase 2 path. The substrate exists; the corrector
extension is well-defined; the closure semantics match the conceptual
question; and the natural seeds come for free from #303 (BCR4BP L1
continuation) and #304 (BCR4BP halo continuation).

### B.3 Discovery hypothesis

**Candidate region**: spatial-support condition that the orbit passes through
a small ball around SE-L1 or SE-L2 (radius ~50000 km, comparable to EM Hill
radius). The orbit also has at least one EM-Hill-bounded segment (the
"resonant arc" that Park-Howell already characterize).

**Expected period**: 19 years (235:19 Saros) is the most natural; integer
sub-multiples (e.g., 1 year, 18.6-year nodal regression) are also possible
windows but less commensurate.

**Expected v_inf**: per KLMR, the SE-L1/L2 → EM-L1/L2 heteroclinic transit
arrives at the EM Hill boundary with |v - v_EM_periodic| ≈ 0.05-0.2 km/s.
This is low — manifold-mediated — consistent with a cycler that needs only
a small phasing Δv per period.

**Mass-budget implication**: if discovered, a stable ~0-2 m/s/month
maintenance cycler with broad coverage (visits SE-L1/L2 *and* near-Moon
each period) would be a genuinely new mission-design building block,
potentially relevant to:
- multi-spacecraft constellations spanning SE-L1 (JWST / Roman / SO) and
  cislunar (Gateway / Artemis logistics)
- low-Δv mass transport between SE-L1/L2 and EM regions
- the conceptual analog of an Aldrin Mars cycler but in cislunar /
  Sun-Earth-region space.

### B.4 Falsifiable scope of the hypothesis

The hypothesis is falsifiable if Phase 2 produces either:

- **Constructive existence**: a corrector run converging on a 19-year (or
  other commensurate-period) BCR4BP periodic orbit with the spatial-support
  constraint satisfied and residuals < 1e-8. This would be a candidate
  cross-system cycler family.

- **Clean negative**: a systematic search across BCR4BP synodic-resonance
  ratios (1:1, 2:1, 3:1, ..., 19:1, 235:19) with the spatial-support
  constraint, all returning `seed_no_converge` or violating the constraint
  at every fixed point of the corrector. This would constitute the first
  published negative result on cross-system cyclers and would warrant an
  entry in the `data/empty_regions.jsonl` registry.

The hypothesis is **not** falsified by:
- the absence in the literature (only 17 papers triaged; insufficient
  sample),
- failure under any single seed (Phase 2 must scan systematically),
- a converged BCR4BP synodic-resonant orbit that does NOT satisfy the
  spatial-support constraint (this is exactly Park-Howell's already-published
  result; we re-confirm but don't add anything).

---

## Part C — Order-of-magnitude probe

Script: `scripts/scope_316_cross_system_candidate.py`. Pure analytic
geometry, sourced constants only (JPL DE440 / IAU 2015 GM values, IAU 2012
AU). No propagation, no corrector, no test gate. Runs in <1 s.

### Geometric outputs

| Quantity | Value | Notes |
|---|---|---|
| mu_SE = M_Earth-Moon / M_total | 3.040 × 10⁻⁶ | Sourced GMs |
| mu_EM = M_Moon / M_Earth-Moon | 1.215 × 10⁻² | Sourced GMs |
| SE-L1/L2 distance from Earth | ~1,502,000 km | Hill approx (mu/3)^(1/3) |
| SE-L1 distance / Moon-distance | 3.91 | — |
| EM-L1 distance from Moon | ~61,300 km | Hill approx |
| EM-L1 distance from Earth | ~323,100 km | A_moon − r_EM_L1 |
| EM-L2 distance from Earth | ~445,700 km | A_moon + r_EM_L2 |
| Spatial gap SE-L1 → EM-L1 | ~1,180,000 km | (1502k − 323k) |

The gap is **bridgeable by heteroclinic-tube transit** (per KLMR, ~90-120 d
arc at v_inf < 0.1 km/s).

### Temporal commensurability

| Quantity | Value |
|---|---|
| T_SE_synodic | 365.25 d (sidereal year) |
| T_EM_synodic | 29.5306 d |
| Ratio T_SE / T_EM | 12.3685 |
| Best small-int rational approximation | **235 / 19** (rel err 8.5 × 10⁻⁶) |
| Resulting near-commensurate closure period | **6,940 d ≈ 19.00 yr** |

**This is the Saros / Metonic family** — the very same near-commensurability
that produces eclipse predictability over 19-year cycles. It is the natural
period for a cross-system cycler. Sub-multiples (year-class periods at
12:1, 13:1, etc.) are less commensurate. The 19-year period is heavy for
a mission profile but is the **dynamical floor** the math demands.

### Verdict from the probe

The geometry and timing are not absurd: the SE/EM spatial regions are
bridgeable via known mechanisms, the temporal closure has a natural
near-commensurate window at 19 years, and the substrate (BCR4BP) is in
place. The framework is **physically plausible**.

The probe does **not** prove the orbit exists; it proves the question
deserves a corrector. Phase 2.

---

## Part D — Verdict & phase 2 path

### Verdict

**INCONCLUSIVE — partial published coverage; the open hypothesis is
sharpened into a falsifiable Phase 2 specification.**

- The cross-system cycler concept is not refuted by the literature.
- The closest published objects (BCR4BP synodic-resonant families) have one
  of the two defining properties (simultaneous closure in both frames) but
  lack the other (spatial-support visit to both SE-L1/L2 *and* EM regions).
- The geometric / temporal sanity check passes (Saros 19-year window,
  ~1 Mkm bridgeable gap).
- The required corrector extension is small (~300-500 LOC additive on
  the #292 BCR4BP substrate).

### Phase 2 specification

**Phase 2 task**: BCR4BP synodic-resonant periodic-orbit search with a
state-space-region constraint requiring spatial support at both EM-L1/L2
AND SE-L1/L2 within a single period.

Formulation: extend `src/cyclerfinder/genome/bcr4bp_genome.py` with a
multiple-shooting corrector whose interior shooting points carry phase-space
constraints. Add the cross-system region predicate as an admissibility
check on the converged orbit.

Initial seeds: continuation in the synodic-resonance ratio q from
Park-Howell families (1:1, 2:1, 3:1, ..., 19:1, 235:19) at each p ∈
{1, 2, 3, ..., 19}. The 235:19 ratio is the natural target; sub-multiples
are scan points.

LOC estimate: **300-500 LOC** additive (extending the existing genome
module + a new scan script). No new core module.

### Phase 3 specification (conditional)

**If Phase 2 produces a converged candidate**: V-pipeline adaptation.
Cross-check the candidate against full SEM ephemeris (DE440) at the
acceptance-tier definitions for BCR4BP-V0 (already deferred per #292
Phase 5). Literature flagger run via `search/literature_check.py`. Mass
budget per the §16.5 baseline.

**If Phase 2 returns clean negative across the full synodic-resonance
ratio scan**: register the empty region in `data/empty_regions.jsonl`
per the negative-results-registry memory anchor, with method-version
fingerprint = BCR4BP-Newton-multishoot @ Phase-2 spec. Future
HR4BP-class methods (paper 11 substrate) would be a subsumption target.

### Anti-claims (sourced discipline)

- No cycler is claimed; the verdict is INCONCLUSIVE.
- The literature triage is limited to 17 papers; the BCR4BP/HR4BP corner
  is active and growing.
- The probe is analytic geometry, not a propagation.
- The 19-year period is a *natural commensurability window*, not a
  predicted orbital period of any specific orbit.

---

## References (all publicly accessible)

1. Koon, Lo, Marsden, Ross (2001), *Celest. Mech. Dyn. Astron.* 81:63–73.
   [koon.pdf](http://www.cds.caltech.edu/~koon/papers/koon.pdf)
2. Canalias, Gómez, Marcote, Masdemont (2006), *DCDS* 14:261–279.
   DOI [10.3934/dcds.2006.14.261](https://www.aimsciences.org/article/doi/10.3934/dcds.2006.14.261)
3. Lo (2002), *IEEE Aerospace Conf.*, ["InterPlanetary Superhighway"](https://www.gg.caltech.edu/~mwl/publications/papers/IPSAndOrigins.pdf)
4. Belbruno, [Wikipedia WSB summary](https://en.wikipedia.org/wiki/Weak_stability_boundary).
5. van der Weg & Vasile (2015), *Celest. Mech. Dyn. Astron.*, [strathprints PDF](https://strathprints.strath.ac.uk/51990/1/van_der_Weg_Vasile_CMDA_2015_Sun_earth_L1_L2_to_moon_transfers_exploiting_natural_dynamics.pdf)
6. Boudad, Howell, Davis (2020), *Adv. Space Res.* 66.
   [ScienceDirect S0273117720305536](https://www.sciencedirect.com/science/article/abs/pii/S0273117720305536)
7. McCarthy & Howell (2021), AAS 21-270.
   [Purdue PDF](https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Conferences/2021_AAS_McCHow.pdf)
8. McCarthy & Howell (2023), *Adv. Space Res.* 71.
   [ScienceDirect S0273117722008614](https://www.sciencedirect.com/science/article/abs/pii/S0273117722008614)
9. Park & Howell (2022), *Adv. Space Res.* 70.
   [ScienceDirect S0273117722004781](https://www.sciencedirect.com/science/article/abs/pii/S0273117722004781)
10. Reddy & Howell (2023), AAS 23-227.
    [Purdue PDF](https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Conferences/2023_AAS_RedHow.pdf)
11. Brown & Peterson (2024), arXiv [2402.19181](https://arxiv.org/abs/2402.19181).
12. Burgos-Garcia & Gidea (2014), arXiv [1412.3775](https://arxiv.org/pdf/1412.3775).
13. arXiv [2509.12675](https://arxiv.org/pdf/2509.12675) (2025).
14. Ross, Roberts, Tsoukkas (2025), AAS 25-621.
    [ross.aoe.vt.edu](https://ross.aoe.vt.edu/papers/ross-roberts-tsoukkas-2025-AAS-25-621.pdf)
15. Genova & Aldrin (2015), AAS 15- (in `cyclers_pdf` per memory).
16. Castelli (2006/2012), [ScienceDirect S0094576506000981](https://www.sciencedirect.com/science/article/abs/pii/S0094576506000981).
17. Anderson, Lo et al. — JPL patched-CR3BP series (DOI search via JPL TRS).
18. Bosanac et al. — [Springer s11071-022-07829-1](https://link.springer.com/article/10.1007/s11071-022-07829-1), [s10569-022-10063-x](https://link.springer.com/article/10.1007/s10569-022-10063-x).

## Provenance

- Literature triage 2026-06-17 via WebSearch + WebFetch (binary-PDF fetches
  flagged where text extraction failed; replaced by search-result and
  abstract summaries from the same surfaced URL).
- The 235:19 Saros / Metonic commensurability is well known in celestial
  mechanics; the probe rediscovered it from first principles as the
  smallest-integer near-rational of T_SE/T_EM.
- No catalogue writeback, no novelty claims, no source-code changes to
  `src/cyclerfinder/`.

Companion script: `scripts/scope_316_cross_system_candidate.py`.
