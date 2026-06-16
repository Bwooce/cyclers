# 2026-06-16 — Physical-sanity flyby gate (#324)

## Motivation: the #312 Umbriel-Oberon-Umbriel manual catch

#312's Phase-1 Part B surfaced one **SILVER**-verdict candidate that passed
every automated guard — closure + cross-check + NN-tier + lit-corpus check +
ML-flagger — and that the agent's manual physical-sanity review correctly
rejected. The reasoning (copied here from
`docs/notes/2026-06-16-312-uranus-extended-sweep.md` §B.3, caveat 1):

> V_inf at Umbriel is 2.27 km/s — 4.2× Umbriel's surface escape velocity
> (√(2 × 85.1 / 584.7) = 0.539 km/s). The closest-approach radius needed
> for a substantial bending angle … at v_inf = 2.27 km/s gives r_p below
> Umbriel's surface for even modest bending. The Umbriel "flyby" in this
> cycler is not a useful gravity assist — it is effectively a V_inf-
> continuity match at a fictitious encounter geometry.

The 2.27 km/s figure is the encounter V_inf magnitude at the *offset-sweep*
best-residual record (`data/scan_312_uranus_umbriel_oberon_offset_sweep.jsonl`,
`kind: best_overall`, residual 0.024 km/s at rel_offset 45°). The SILVER-row
itself (the `candidate_id: repeated-moon-uranus-00000041` entry in
`data/scan_312_uranus_oberon_umbriel.jsonl`) sits at the *coarser* phase grid
with V_inf ≈ 0.92 km/s at Umbriel and a 0.025 km/s residual; both rows belong
to the same family, but at the refined optimum the V_inf rises and the geometry
becomes unphysical. The genome ceiling and the geometric ceiling are entangled.

The gap that needed closing: the **automated** pipeline had no equivalent of the
agent's "is this flyby physically useful?" question — so a candidate could be
admitted (or refined into a tighter-residual neighbour) at a V_inf where the
patched-conic bend is essentially zero, and survive every other guard.

## The physics

At a hyperbolic flyby with minimum-safe periapsis `r_p = r_body + alt_min`,
the maximum ballistic deflection angle is (BMW §6.4; also documented in
`src/cyclerfinder/core/flyby.py:40-70`)

```
sin(delta_max / 2) = mu_body / (mu_body + r_p * V_inf^2)
```

A "useful" flyby must rotate the V_inf asymptote by some finite amount. The
threshold adopted here is **5° at the body's safe-altitude default**. This is
a judgment call — the standing rule
(`feedback_golden_tests_sourced_only`) forbids us from labelling a
non-sourced threshold as a golden, so the threshold is a **knob** rather than
a constant:

* Real cycler / tour flybys engineered in the open literature sit comfortably
  above any sensible floor: Galileo's Earth flybys at V_inf ≈ 6.232 km/s
  reach ≈ 74.6° max bend at the 300 km safe altitude (D'Amario, Bright &
  Wolf, Space Sci Rev 60:23, 1992); Cassini's Venus flybys at V_inf ≈ 7 km/s
  reach ≈ 61.4° (Peralta & Flanagan, AAS 95-117, 1995); the classic Aldrin
  Mars flyby at V_inf ≈ 5.5 km/s reaches ≈ 32.2° (Russell & Ocampo 2005;
  McConaghy 2002). 5° is **two orders of magnitude** below any of these.
* Below ~5° the flyby contributes deflection at the level of trajectory-
  correction-manoeuvre (TCM) noise; calling it a gravity assist is
  effectively sleight of hand. The Umbriel-at-2.27 km/s case gives 2.70°.
* 5° is **stricter than zero** (catches the Umbriel pathology) but **looser
  than every operationally interesting flyby**, so a gate-PASS does NOT
  certify usefulness — it only certifies not-pathological. Candidates that
  pass remain subject to lit-check, ML, and the full V0–V5 gauntlet.

## Implementation

`src/cyclerfinder/search/physical_sanity.py` — a thin wrapper around
`cyclerfinder.core.flyby.max_bend` (the patched-conic formula above). No new
physics, no modification of `core/flyby.py`. The module exposes:

* `flyby_is_useful(body, vinf_kms, *, min_safe_altitude_km=None,
  min_useful_bend_deg=5.0) -> FlybyPhysicalVerdict` — per-encounter scalar
  verdict.
* `candidate_passes_physical_gate(sequence, vinf_kms_per_encounter, *,
  min_useful_bend_deg=5.0, per_body_min_safe_altitude_km=None) ->
  (bool, list[FlybyPhysicalVerdict])` — sequence-level check; a multi-leg
  tour is admitted iff every encounter passes.

Body lookup tries `PLANETS` then `SATELLITES`; unknown body raises `KeyError`
(the gate must NEVER silently admit an unknown body, per the same discipline
that the constants module enforces).

Tests in `tests/search/test_physical_sanity.py` (26 cases; all pass — 18
direct + 8 parametrized):

* Motivating Umbriel case at V_inf = 2.27 km/s → reject (max bend 2.70°)
* Actual #312 SILVER row's V_inf = 0.92 km/s at Umbriel → admit (14.71°)
* Galileo Earth flyby (V_inf 6.232 km/s) → admit (74.57°)
* Aldrin Mars flyby (V_inf 5.5 km/s) → admit (32.16°)
* Cassini Venus flyby (V_inf 7.0 km/s) → admit (61.42°)
* Sequence-level worst-case Umbriel-Oberon-Umbriel @ (2.27, 0.98, 2.27) →
  reject (Umbriel#1 fails first)
* Sequence-level actual SILVER @ (0.92, 0.96, 0.89) → admit
* Threshold + altitude knob behaviour
* Parametrized consistency check vs `core/flyby.py::max_bend` for 8 (body,
  V_inf) pairs spanning planets + moons.

## Re-run sweep (per `feedback_bugfix_invalidates_past_searches`)

The standing rule mandates that every correctness fix sweep over past
results: buggy / missing-guard solvers are false-negative generators (and
false-positive admitters). The #324 gate is bug-fix-equivalent — it would
have rejected candidates that earlier reports admitted.

`scripts/rerun_324_physical_gate.py` runs the gate over every recent
session JSONL with per-encounter V_inf:

| File | Notes |
|---|---|
| `scan_285_saturn.jsonl`, `scan_285_uranus.jsonl` | repeated-moon scans |
| `scan_312_uranus_*.jsonl` (10 files) | per-pair JSONLs + offset_sweep + 3D probe + robustness |
| `scan_313_mars_phobos.jsonl`, `scan_313_mars_deimos.jsonl` | tulip-orbit ICs (no encounter V_inf — skipped) |
| `scan_313_sun_jupiter_europa.jsonl`, `scan_313_sun_jupiter_io.jsonl` | CR3BP ICs — skipped |
| `scan_309_low_thrust_em.jsonl`, `scan_309_low_thrust_vem.jsonl` | low-thrust EM/VEM |
| `scan_298_galileo_veega.jsonl` | Galileo VEEGA pre-screen |
| `precursor_302_aldrin.jsonl`, `precursor_302_s1l1.jsonl` | precursor-MGA |

Output: `data/rerun_324_physical_gate.jsonl` (registry).

**Results:**

```
rows_scanned          = 1013
rows_with_vinfs       = 811
rows_skipped_no_vinfs = 169   # CR3BP tulip ICs, etc.
passed_gate           = 781
failed_gate           =  30
SILVER passed gate    =   1   # the #312 row at coarse-grid V_inf=0.92 km/s
SILVER failed gate    =   0
```

* **The single SILVER survivor passes the gate** at its actual coarse-grid
  V_inf = 0.92 km/s at Umbriel (max bend 14.7°). The gate's verdict is
  consistent with the agent's caveat: at the coarse-grid magnitudes the
  geometry is fine; the candidate fails for OTHER reasons (genome ceiling /
  unmodelled J2 — see #312 doc §B.3 caveats 2-3). The gate is a floor, not
  a ceiling; admit ≠ admissible-as-a-cycler.
* **The 10 offset-sweep rows that DO trigger the gate** are exactly the
  family that surfaced V_inf ≈ 2.28 km/s at Umbriel (max bend 2.67°) under
  the finer 96×96 grid — the precise records that the agent's manual review
  flagged. Their residual is 0.024 km/s (TIGHTER than the SILVER), so under
  a future "smaller residual = better" sorter they would have surfaced as
  the family representative; the gate now catches them automatically.
* The remaining 20 failures cluster in scan_298 (Galileo VEEGA Phase-2
  shell-seed rows with non-converged 23-44 km/s V_inf at Earth — expected,
  per the file's own notes) and the precursor_302 search-history rows
  (Aldrin / S1L1 closures with V_inf ≥ 24 km/s at Earth). None were ever
  carried as SILVER. The gate's contribution this time is the
  **forward-looking guard against future false-positives**, not a
  retrospective catch.
* The gate has no sibling SILVER false positives in the present session
  output.

## Phase 2 recommendation: integrate into the discovery daemon

`scripts/discovery_campaign_daemon.py` should call
`candidate_passes_physical_gate` between the closure / cross-check tier and
the lit-check tier. Concretely: after a candidate clears closure +
cross-check + NN, BEFORE it consumes a lit-check query (which is the most
expensive automated step), run the physical-sanity check; if it fails,
short-circuit the candidate with a structured rejection reason
(`PHYSICALLY_UNUSABLE_FLYBY`) and skip lit/ML entirely. This both prevents
the #312-shaped trap AND saves lit-corpus API budget on candidates that
shouldn't proceed.

The gate's threshold (`min_useful_bend_deg`) should be exposed on the daemon
CLI so an operator can run a tighter sweep (e.g. 20° to discard not-just-
pathological-but-uninteresting flybys) without re-deploying.

## Discipline

* NO catalogue writeback. Gate-pass is not a SILVER promotion; it is a
  floor that other guards still run on top of.
* NO novelty claims. Per `feedback_literature_novelty_check_baseline`,
  novelty is decided downstream of literature_check.py; the gate is upstream.
* `core/flyby.py` not modified — wrapped. The threshold lives in
  `physical_sanity.py` so the cost-only primitives stay clean of judgment
  constants.
* Existing JSONLs not modified — registry is a new file.
* Sourced anchors in the tests are mission-shaped V_inf magnitudes (Galileo,
  Cassini, Aldrin), not sourced *thresholds*. The 5° floor is openly a
  judgment call, per `feedback_respectful_errata_framing`.
