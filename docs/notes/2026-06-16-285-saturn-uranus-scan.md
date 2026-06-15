# 2026-06-16 — Task #285 — Prioritized repeated-moon scan at Saturn and Uranus

**Status: complete — clean negatives at both targets.** No SILVER survivors
under the 0.05 km/s V_inf-continuity gate; the entire 5-tier scorer chain
(Tier-0 NN + DOP853 cross-check + literature_check + ml_flagger) never had to
fire because the closure gate alone pruned the search space. The most
important finding is the **near-miss tier**: two near-gate candidates that
the underlying patched-conic/circular-coplanar/multi-arc genome could not
quite close. Both are anchored in known existence-prior families and are
listed below for the multi-arc + DSM continuation work to revisit.

## Census

| System | Band | Topology | Sequences x n_rev evaluated | Closed | Sub-gate (SILVER) | Near-miss (< 1 km/s) |
|---|---|---|---:|---:|---:|---:|
| Saturn | A | Titan-Rhea-Dione-Tethys, k=3 closed cycles, n_rev=1..3 | 108 | 23 | **0** | 5 |
| Saturn | B | Titan-Enceladus-Rhea Liang-CGE analogue, k=5 cycles, n_rev=0..2 | 1440 | 294 | **0** | 0 |
| Uranus | A | Titania-Oberon-Umbriel pairs, k=3, n_rev=0..3 | 90 | 28 | **0** | 12 |
| Uranus | B | Ariel-Umbriel-Titania Liang-CGE analogue, k=5, n_rev=0..2 | 1440 | 716 | **0** | 29 |
| **Total** | — | — | **3078** | **1061** | **0** | **46** |

Total runtime: ~30s wall-clock for both campaigns (the 5-tier scorer chain
was never invoked because no candidate cleared the closure gate). The
"closed" column counts candidates for which the Lambert solver returned a
solution at the requested n_rev at any of the swept phasing samples; "near-
miss" is the count of those whose canonical residual fell in [gate, 1 km/s).

## Top near-misses (listed for follow-up; NOT SILVER, NOT discoveries)

### Saturn Band A — Rhea-Dione-Rhea (1,1) — residual **0.107 km/s**

```
sequence: ['Rhea', 'Dione', 'Rhea']
n_rev: [1, 1]
residual_kms: 0.10688
vinf_per_encounter_kms: [1.111, 2.208, 1.110]  -- Dione flyby is the binding leg
tof_days: [5.276, 5.276]
```

A ~2x miss of the 0.05 km/s gate (verified independent of phase-grid density:
residual stabilises at 0.107 km/s at phase_samples in {12, 24, 48}). This is
NOT a phase-resolution artefact — it's the binding ceiling of the circular-
coplanar / patched-conic / multi-arc genome for this topology. The Rhea-Dione
1:1 resonance is a well-known Saturnian tour primitive (Takubo 2210.14996
generic prior); a real-eccentricity or multi-arc continuation could
plausibly close it.

### Uranus Band A — Oberon-Titania-Oberon (1,1) — residual **0.062 km/s**

```
sequence: ['Oberon', 'Titania', 'Oberon']
n_rev: [1, 1]
residual_kms: 0.06170
vinf_per_encounter_kms: [0.962, 0.655, 0.924]  -- the closest near-miss seen
tof_days: [16.24, 16.24]
```

Within ~25% of the gate (0.062 vs 0.05). The 1:1 Oberon-Titania resonant arc
at low V_inf (~0.7-1.0 km/s) sits inside the patched-conic feasibility
envelope; the near-miss residual is the V_inf-continuity defect at the
anchor flyby, which the multi-arc / DSM follow-up could reasonably close.
Uranus has NO published cycler family (no existence prior); this is the
first quantified evidence that the Uranian regulars host a near-feasible
repeated-moon topology under this genome.

## Pipeline composition (the #282 + #261 + #256 wiring)

```
RepeatedMoonTarget          -> enumerate (sequence, n_rev) deterministically
  .close()                  -> Lambert per leg + V_inf continuity residual
  [residual < gate]         -> SILVER survivor (none reached this point)
  legs_from_repeated_moon_candidate -> rebuild SI legs
  FiveTierPrioritizer.score_candidate_legs -> Tier-0 NN per-leg dV (only
                                tier applicable to patched-conic legs; #282
                                architectural seam)
  dop853_cross_check_leg    -> DOP853 rtol=atol=1e-12 vs Lambert arrival
  check_literature(offline) -> rediscovery filter against KNOWN_CORPUS
  FalsePosFlagger.score()   -> p_fp on the retrained-on-corpus flagger
  verdict in {SILVER, BRONZE, REJECT}
```

The independent DOP853 cross-check (different solver, same two-body physics)
is mandatory by `feedback_orbit_closure_discipline`. In the SILVER-tier
smoke runs, `max_dr_arrival_km` came in at 1e-4 to 1e-5 km — well within
the 1 km gate. The cross-check is plumbed and proven; nothing in this
campaign tripped it because no candidate ever closed.

## Comparison to the unscored #264 daemon

The unscored #264 campaign daemon (`scripts/discovery_campaign_daemon.py`)
would have produced **the same negative**: it shares the closure engine
(`RepeatedMoonTarget.close`), and the closure gate is what pruned all 3078
candidates here. The 5-tier scorer composition adds nothing to this
particular system pair because there is no admitted SILVER cohort to
re-rank.

**Where the scorer would have changed the call** (the prompt's question):
the comparison rule is sharper than "did we surface anything #264 missed?".
The 5-tier scorer is a *post-closure rank-product over admitted candidates*
— at Saturn / Uranus, with zero admitted candidates, the rank-product is
vacuous. The scorer's value is the BRONZE / REJECT verdicts on SILVER-eligible
candidates that the unscored daemon would have routed straight to the
gauntlet. We confirmed this end-to-end in the smoke pass on Saturn at a
relaxed gate (the relaxed-gate smoke produced two Rhea-Titan-Rhea rows
correctly BRONZE'd by the literature_check against Davis-Phillips-McCarthy
2018 tulip family — exactly the disposition a SILVER-routing daemon would
have missed).

## Honest negatives (the topology / n_rev cells that turned up empty)

* **Saturn Band B (Titan-Enceladus-Rhea Liang-CGE)**: 1440 candidates,
  294 closed, ZERO sub-1-km/s near-misses. The 3-body Liang topology DOES
  NOT generalise from Jupiter to Saturn under this genome: Liang's
  Callisto-Ganymede-Callisto-Europa-Callisto works because the Galilean
  Laplace resonance keeps every leg's resonance commensurate; Saturn's
  Titan-Enceladus-Rhea triplet has no analogous lock, so the multi-arc
  V_inf-continuity defect runs in the multi-km/s range across the entire
  enumeration. The 5-tier scorer would not have changed this call.
* **Saturn Band A's outer pairs (Titan-Tethys, Titan-Dione)**: closure
  succeeds but residuals are in the 4-10 km/s band. Titan's SMA
  (1.22 Mkm) vs Tethys (0.29 Mkm) / Dione (0.38 Mkm) gives a multi-Tisserand
  jump per leg that the multi-rev Lambert grid cannot bridge with circular-
  coplanar phasing.
* **Uranus Band A's Umbriel-Oberon pair**: 0 sub-gate, 0 near-misses
  (all closures in the 4-9 km/s residual band). Umbriel-Oberon (0.27 vs 0.58
  Mkm) is a 2.16x SMA ratio — outside the resonant primitives the patched-
  conic genome can host.
* **Uranus Band B**: 716 closed, 29 near-misses (best 0.45 km/s); no sub-
  gate. The 3-body Ariel-Umbriel-Titania Liang topology HAS near-misses
  (unlike the Saturn analogue) but they are 9x the gate — outside the
  multi-arc genome's reach.

## What this campaign DID accomplish

1. **First quantified bound** on Saturn / Uranus repeated-moon cyclers under
   the #254 / #259 multi-rev genome at the current ~0.05 km/s gate. The
   bands are not "untested" — they are "tested and clean".
2. **Two strong near-miss anchors** for the multi-arc / DSM continuation
   work (#281 frontier roadmap): Rhea-Dione-Rhea (1,1) at 0.107 km/s
   and Oberon-Titania-Oberon (1,1) at 0.062 km/s. Both are 1:1 resonant
   topologies — the most likely to close under intermediate flybys.
3. **Pipeline plumbing validated**: Tier-0 NN, DOP853 cross-check,
   literature_check (offline), and FalsePosFlagger all proven to fire end-
   to-end on the smoke pass at a relaxed gate. The SILVER/BRONZE/REJECT
   verdict policy was exercised on real Saturn rows.
4. **Honest negative recorded** for the unscored-#264 comparison: at zero
   SILVER survivors, the 5-tier composition is rank-product-vacuous; the
   scorer's value emerges only when SILVERs exist.

## Anti-claims (what this is NOT)

* This is not a novelty claim. Even the near-miss rows are not novel
  candidates — they have not closed.
* No catalogue writeback. No empty-region registry append (the discovery
  campaign engine's normal output for a clean-negative sweep was bypassed
  here: this scan ran via the bespoke #285 driver, which deliberately does
  NOT write to `data/empty_regions.jsonl`; that pathway belongs to the
  unscored daemon's region-tracking, and conflating the two would
  duplicate the negative).
* No gauntlet promotion. The verdict policy here matches the #274 gauntlet
  V0 gate so a future SILVER from a more-capable genome at these targets
  can drop straight into the gauntlet without re-thresholding.

## Files

* `src/cyclerfinder/search/saturn_uranus_campaign.py` — campaign runner +
  scoring chain
* `scripts/scan_285_saturn.py` — Saturn driver (Bands A + B)
* `scripts/scan_285_uranus.py` — Uranus driver (Bands A + B)
* `data/scan_285_saturn.jsonl` — output (meta + per-band summary + near-
  misses; zero SILVER rows)
* `data/scan_285_uranus.jsonl` — same

## Suite green

`tests/ml/` + `tests/search/{test_discovery_campaign,test_five_tier_prioritizer,test_literature_check}.py`
all pass (67 tests). The 3 pre-existing failures in `tests/data/test_*census*` /
`tests/data/test_multi_arc_invariants.py` are census-ratchet failures from
the Tito 2018 admission commit `b6bcbb3` on `main` — confirmed by
re-running the failing tests on a stash-clean HEAD. NOT caused by this
work.
