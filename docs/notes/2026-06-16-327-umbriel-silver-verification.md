# 2026-06-16 — #327 Umbriel SILVER verification at gate-passing V_inf

## What was open

The #324 physical-sanity gate (5 deg max-bend floor) was added AFTER the
#312 Umbriel-Oberon-Umbriel SILVER row was minted. The #312 doc's
caveat-1 had cited V_inf = 2.27 km/s at Umbriel (4.2x escape, 2.70 deg
max bend) as the headline pathology — but **2.27 km/s is the
offset-sweep refinement record** at rel_offset = 45 deg, NOT the stored
SILVER row. The actual SILVER row
(`repeated-moon-uranus-00000041` in
`data/scan_312_uranus_oberon_umbriel.jsonl`) sits at V_inf = 0.92 / 0.96
/ 0.89 km/s with max bend 14.7 deg at Umbriel — it **passes the #324
gate cleanly**.

The #324 doc itself (§ Results) records that the SILVER passes the gate
at its coarse-grid V_inf and that the 10 offset-sweep records (residual
0.024 km/s at V_inf 2.28 km/s) fail. That part is correct. What was
*unverified* until now: **is the SILVER's 0.025 km/s closure real at
gate-passing V_inf, or is it only the unphysical refinement that gets
to low residual?**

## Approach

`scripts/verify_327_umbriel_silver.py` wraps the production primitives
(no modification of `discovery_campaign.py`, `saturn_uranus_campaign.py`,
or `physical_sanity.py`) and answers seven questions on the SILVER's
stored IC. Output: `data/silver_327_verified.jsonl` (10 rows).

## Findings

### Closure under each convention (production `RepeatedMoonTarget.close`)

| Convention                     | Moons (sorted seed offsets) | Residual (km/s) | Max V_inf (km/s) | Gate-passing? |
|--------------------------------|---|----------------:|-----------------:|:-------------:|
| 2-moon                         | `(Oberon@0, Umbriel@180 deg)`        |    **0.025232** |       **0.9604** |     yes       |
| 3-moon                         | `(Oberon@0, Titania@120, Umbriel@240)` |       0.635986 |           2.1108 |      no       |

Both reproduce #312's stored values byte-perfect. The convention split
is the deterministic-offset seed; same physics, different `theta0` per
moon.

### Relative-offset basin floor (confirmation sweep, 24x24)

| Region                                                    | Residual (km/s) | rel_off (deg) | Max V_inf (km/s) | Gate-passing? |
|-----------------------------------------------------------|----------------:|--------------:|-----------------:|:-------------:|
| Overall basin floor (NO V_inf constraint)                 |    **0.024033** |          45.0 |       **2.2827** |      no       |
| Best gate-passing sample (max V_inf <= 1.0 km/s)          |    **0.025232** |         180.0 |       **0.9604** |     yes       |

The overall basin floor at 0.024 km/s reproduces #312's 96x96 finding
even at 1/16th the grid density (24x24). Critically, **the
gate-passing region has its OWN local floor at 0.025232 km/s** — and
that floor is the SILVER row itself. The 0.024 km/s "deeper" floor only
exists at unphysical V_inf. Within the physically usable envelope, the
SILVER IS the optimum.

### n_rev sweep at the 2-moon convention

Only `(n_rev_1, n_rev_2) = (1, 1)` closes below the 0.05 km/s gate.
Top 5 by residual:

| n_rev | Residual (km/s) | Max V_inf (km/s) |
|:-----:|----------------:|-----------------:|
| (1, 1)  |        **0.025232** |       0.9604 |
| (0, 1)  |        1.064904 |       4.9004 |
| (1, 0)  |        1.179092 |       1.5640 |

Adding revolutions does NOT close the gap further — they hurt. The
gate-passing residual is intrinsic to the (1, 1) genome cell.

### DOP853 independent cross-check (rtol = atol = 1e-12)

| Leg                  | dr_arrival (km) | dv_arrival (km/s) | Passed (<= 1 km) |
|----------------------|----------------:|------------------:|:----------------:|
| Umbriel -> Oberon    |       3.902e-6  |          2.723e-11 |       yes        |
| Oberon -> Umbriel    |       1.596e-5  |          2.296e-10 |       yes        |

`max dr_arrival = 1.596e-5 km = 2.736e-11 nondim` (LU = Oberon SMA =
583511 km). This is **5 orders of magnitude below** the 1e-6 nondim
orbit-closure discipline gate (`feedback_orbit_closure_discipline`).
The closure is independently confirmed in the patched-conic envelope.

### Physical-sanity gate at stored V_inf (#324 wrapper)

| Encounter   | V_inf (km/s) | Max bend (deg) | Useful? |
|-------------|-------------:|---------------:|:-------:|
| Umbriel #1  |        0.9199 |          14.71 |   yes   |
| Oberon      |        0.9604 |          23.70 |   yes   |
| Umbriel #2  |        0.8947 |          15.45 |   yes   |

All three encounters clear the 5 deg floor; the central Oberon flyby
clears it 4.7x over. Gate verdict: **PASSED**.

### Literature check (offline_corpus_search, 35+ KNOWN_CORPUS anchors)

Signature:
`primary=Uranus, sequence=(Umbriel, Oberon, Umbriel), period_k=2,
vinf=(0.9199, 0.9604, 0.8947) km/s, n_rev=(1, 1)`

`check_literature`: `status = not-found`, `confidence = 0.40`,
`citation = None`. No KNOWN_CORPUS anchor names any Uranian body.
This is the same offline result as #312 — necessary-not-sufficient
for novelty, per `feedback_literature_novelty_check_baseline`.

### ML flagger (logistic regression over hand-crafted features)

`p_fp = 0.591774` — IDENTICAL to the #312 stored row (the flagger
features only depend on the stored IC, not the convention used to
verify it). Production `P_FP_SILVER_MAX = 0.75` from `#274`; `0.59 <
0.75` is a SILVER pass, not a BRONZE downgrade.

## Verdict

**`REAL_CANDIDATE_AWAITING_306`** — every guard passes at the
gate-passing IC:

* Closure: 0.025232 km/s (gate 0.05 km/s) — PASS
* Independent DOP853 cross-check: 2.7e-11 nondim (gate 1e-6) — PASS
* Physical-sanity (#324): all 3 encounters useful — PASS
* Basin-floor gate-passing region floor: 0.025232 km/s (same as SILVER) — PASS
* Literature: not-found — PASS (necessary, not sufficient)
* ML flagger: p_fp = 0.5918 < 0.75 — PASS

The candidate is **not a coarse-grid artifact** at the gate-passing IC.
The 0.024 km/s offset-sweep refinement is a TRUE physical artifact of
its own (it pulls V_inf at Umbriel to 2.28 km/s where the bend is 2.7
deg), but the SILVER row stands on its own at gate-passing V_inf with
0.025 km/s residual, 5-orders-of-magnitude-below-required DOP853
agreement, and useful bends at every encounter.

## Path forward

* **NO catalogue admission.** Per the discipline anchors, the SILVER is
  plumbed for the #306 3D V0-V5 gauntlet (Uranus moons; currently
  planar-EM-only). Admission must wait for #306 to gauntlet the planar
  Lambert IC against full 3D propagation in CR3BP(Uranus, Oberon) and
  CR3BP(Uranus, Umbriel).
* **#312 caveats still apply.** The doc's caveats 2 (J2 ignored) and 3
  (2-moon convention dependence) are unaffected by this verification.
  J2 effects would shift the residual; the gauntlet's higher-fidelity
  step is where that gets resolved.
* **NO novelty claims** until #306 closes AND a wider literature pass
  (Anderson 2025 Uranus-cycler design-space pull) returns clean.

## Anti-claim ledger

What this verification does NOT establish:

* That the candidate is novel — `not-found` in the offline 35-anchor
  corpus is necessary-not-sufficient. The Uranus-moon-cycler literature
  is small but non-empty; a full online sweep is the V4 gauntlet's job.
* That the candidate exists as a true periodic orbit — the patched-conic
  ballistic closure is a topological signal, not an admission. CR3BP /
  J2-augmented closure at the gauntlet's tolerance is the next gate.
* That the deep basin floor (0.024 km/s @ V_inf 2.28) is invalid as
  GENOMES — it could still surface a real cycler in a multi-flyby genome
  where the unphysical encounter is replaced by a chemical / leveraging
  arc (e.g. the #226 FBS optimizer's territory). That is not in scope
  here.

## File index

* Script: `scripts/verify_327_umbriel_silver.py`
* Output: `data/silver_327_verified.jsonl`
* Source SILVER row: `data/scan_312_uranus_oberon_umbriel.jsonl`
  (`repeated-moon-uranus-00000041`)
* Doc: this file (`docs/notes/2026-06-16-327-umbriel-silver-verification.md`)
* Prior context: `docs/notes/2026-06-16-312-uranus-extended-sweep.md`,
  `docs/notes/2026-06-16-324-physical-sanity-gate.md`
