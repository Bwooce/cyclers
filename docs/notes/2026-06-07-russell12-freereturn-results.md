# Russell-12 like-for-like — #137 free-return genome results (2026-06-07)

Recovered run (the #137 agent was killed before reporting; its fix commits
37b9bbd/8fc2ffb survived — this note records the re-run of its committed
campaign). Circular model, like-for-like, `--genome free-return`,
128 epochs × 14 workers, ~5 s wall.

## Acceptance gate
`tests/search/test_russell12_likeforlike_probe.py` — **4 passed** (10.5 s):
the radial-crossing genome makes the sourced geometry representable
(truth-residual ≈ 0 where the old genome measured 3.2–37.5 km/s).

## Per-row outcomes (V∞ DERIVED, never imposed — non-circular evidence)
| row | outcome |
|---|---|
| mcconaghy-2006-em-k2 | CLOSE-AND-MATCH |
| russell-ch4-4.991gG2 | CLOSE-AND-MATCH |
| russell-ch4-8.049gGf2 | CLOSE-AND-MATCH |
| russell-ch4-9.353Gg2 | CLOSE-OFF-ANCHOR |
| russell-ch4-3.64gGg3 | CLOSE-AND-MATCH |
| russell-ch4-3.78Gg3 | CLOSE-AND-MATCH |
| russell-ch4-5.30gGf3 | CLOSE-AND-MATCH |
| russell-ch4-9.94Gg3 | CLOSE-AND-MATCH |
| russell-ch4-3.66gfF3 | CLOSE-MATCH-SYMMETRIC-ONLY |
| russell-ch4-5.30ggF3 | CLOSE-MATCH-SYMMETRIC-ONLY |
| russell-ch4-5.75ggF3 | CLOSE-MATCH-SYMMETRIC-ONLY (vinf Δ 0.04/0.01) |
| russell-ch4-6.44Gg3 | NO-CLOSE (res 0.2666 — near the 0.1 floor) |

**Counts (original 256-point scan): 7 CLOSE-AND-MATCH, 3 SYMMETRIC-ONLY,
1 OFF-ANCHOR, 1 NO-CLOSE** (was 0/9-11/1-3 on the lambert genome, #125/#135).

**Counts (Part 3 dense phase scan, the current default): 8 CLOSE-AND-MATCH,
3 SYMMETRIC-ONLY, 1 OFF-ANCHOR, 0 NO-CLOSE** — `9.353Gg2` (OFF-ANCHOR → MATCH:
phase-grid resolution) and `6.44Gg3` (NO-CLOSE → OFF-ANCHOR: genuinely different
family). See Part 3 below.

## Part 1 — §14 V1 writeback (#137 follow-up, 2026-06-07)

USER-approved principle: like-for-like reproduction qualifies for §14 V1 when the
§14 mechanics pass. Applied the LITERAL §14 V1 mechanics to each matched row's
closed free-return geometry, like-for-like on the CIRCULAR ephemeris (a
circular-coplanar reproduction of a circular-coplanar source):

* path (a) — every leg re-solved with `lamberthub` izzo2015 + gooding1990,
  agreement < `V1_TOLERANCE_MPS` (1e-3 m/s);
* path (c) — Kepler forward re-propagation residual < `KEPLER_REPROP_TOL_KM` (1 km);

both reused verbatim from `cyclerfinder.verify.agreement.crosscheck_code_paths`
(§14 V1 is exactly paths a+c; path b is an extra Forge witness, not part of §14).
Module: `src/cyclerfinder/search/free_return_v1.py`; evidence test:
`tests/search/test_free_return_v1_mechanics.py` (10 passed).

**The decisive honesty gate.** A single free-return ellipse forms a CLOSED,
V_inf-continuous E->M->E cycler only when its descending Earth crossing coincides
with where Earth actually is on the circular orbit. The return-Earth phase mismatch
splits the rows cleanly:

| row | outcome | return-Earth Δphase | Mars V_inf continuity | §14 V1 |
|---|---|---|---|---|
| russell-ch4-5.30gGf3 | CLOSE-AND-MATCH | -0.2° | 0.01 km/s | **PASS → V1** |
| russell-ch4-9.94Gg3 | CLOSE-AND-MATCH | -0.5° | 0.04 km/s | **PASS → V1** |
| russell-ch4-5.75ggF3 | SYMMETRIC-ONLY | -6.1° | 0.18 km/s | **PASS → V1** |
| mcconaghy-2006-em-k2 | CLOSE-AND-MATCH | 175.2° | 24.35 km/s | refused |
| russell-ch4-4.991gG2 | CLOSE-AND-MATCH | 176.0° | 24.58 km/s | refused |
| russell-ch4-3.64gGg3 | CLOSE-AND-MATCH | 179.1° | 22.79 km/s | refused |
| russell-ch4-3.78Gg3 | CLOSE-AND-MATCH | 179.2° | 23.52 km/s | refused |
| russell-ch4-3.66gfF3 | SYMMETRIC-ONLY | 175.9° | 21.93 km/s | refused |
| russell-ch4-5.30ggF3 | SYMMETRIC-ONLY | 170.5° | 23.68 km/s | refused |
| russell-ch4-8.049gGf2 | CLOSE-AND-MATCH | n/a | Lambert-singular | refused |

For the three PASS rows the reconstructed E->M->E arc is a genuine closed cycler:
lamberthub agreement ~1e-8 m/s (≪ 1e-3), Kepler reprop ~3-9e-4 km (≪ 1 km), Mars
V_inf continuous to ≤ 0.18 km/s. For the six 175-180° rows the single free-return
ellipse is genuinely MULTI-ARC — its return needs the eliminated phasing loops, so
forcing a Lambert return leg lands a ~24 km/s Mars V_inf discontinuity (a broken,
non-physical trajectory whose lamberthub/Kepler self-consistency would pass
VACUOUSLY). `VINF_CONTINUITY_TOL_KMS = 0.5` (the campaign match tolerance, not
tightened) rejects exactly that. `8.049gGf2` (deep aphelion, 93-d transit) is
single-rev-Lambert-singular on reconstruction. None of the six promote.

The evidence chain per PASS row: SOURCED aphelion + transit seed → free-return
closure → DERIVED V_inf matches sourced (already shown) → §14 V1 mechanics pass on
a closed, V_inf-continuous reconstruction. **Scope (honest):** this is a
circular-coplanar reproduction of a circular-coplanar source — like-for-like; it is
NOT a real-ephemeris (V3) result.

Writeback (same commit): `_LEVEL_EVIDENCE` extended (3 new entries),
`_LEVEL_BY_ID` extended, `scripts/backfill_validation_level.py` applied
(V0→V1 on the 3 rows), the v4.5 census ratchet
(`tests/data/test_schema_v45_fields.py`) and `data/README.md` updated.

**Validation census (the headline): 6 rows are V1** — the Aldrin pair
(outbound + inbound, #125) and four Russell free-return rows: the three above
(#137 Part 1) plus `russell-ch4-9.353Gg2` (#137 Part 3, promoted by the dense
phase scan; see Part 3). Everything else is V0 / untagged. No row claims V2+.

## Caveats
- The runner's `--probe-at-truth` section still prints the OLD lambert-genome
  diagnostics (script wiring: the probe path predates `--genome`); the
  authoritative truth-residual≈0 pin is the committed acceptance test.
  (Part 3 fixes this — the probe now honours `--genome`.)
- Implication: the Jones elimination chain's negatives were computed on the
  old genome — re-validation on the free-return genome is now warranted
  (#133), as is the near-miss survey. (Part 2 assesses the extension.)

## Part 2 — Jones VEM extension assessment (STOP/report; NOT modest)

ASSESSED whether the free-return (radial-crossing) genome extends to the Jones
VEM topologies (`jones-2017-vem-emevve-outbound` E-M-E-V-V-E /
`...-meevem-inbound` M-E-E-V-E-M). **Verdict: the extension is NOT modest — a
half-genome was deliberately NOT improvised (per the brief's STOP rule).**

### Why the genome's primitive does not carry over

The free-return genome's representable primitive is **one heliocentric ellipse
that radially crosses TWO DISTINCT body radii** (inner, outer); the per-body V_inf
and the leg ToFs EMERGE from the single shape DOF `(a, e)` riding the Mars-V_inf
ridge. The entire E->M->E free return is that ONE ellipse (two halves of a
symmetric arc). Three structural facts break this for the Jones VEM rows:

1. **Multiple distinct transfer ellipses coupled through a flyby.** EMEVVE/MEEVEM
   each contain BOTH an Earth-Mars arc AND an Earth-Venus arc — two ellipses with
   independent `(a, e)`, coupled at the Venus (and Earth) flybys by V_inf-magnitude
   continuity. The genome's single-ridge premise (everything emerges from one
   ellipse's `(a, e)`) has no DOF for a second ellipse coupled through a flyby.
2. **Same-body resonant legs.** EMEVVE has a `V->V` leg; MEEVEM has an `E->E` leg.
   These are full-revolution resonant loops returning to the SAME body radius —
   categorically NOT a two-radius crossing, so the radial-crossing primitive
   cannot represent them at all. (They are exactly the E-E loops the OLD Lambert
   genome modelled — a different primitive.)
3. **Venus is a powered/bending intermediary, and the gap is its flyby V_inf
   MAGNITUDE.** The #110/#120/#122 surveys already established that reaching this
   family is a leg-|V_inf|-magnitude / leg-topology problem, refuted for
   3D-inclination, vector-residual, and B-plane targeting. A genome that derives
   V_inf from one ellipse cannot independently set the Venus-flyby bend.

A control probe confirms the E-M *sub-arc* is NOT the blocker: the free-return E-M
primitive reaches a Mars-V_inf floor of ~2.81 km/s (sweeping `(a, e)`), comfortably
spanning the Jones sourced Mars V_inf (2.42-3.12) — i.e. the genome is a genuinely
lower-V_inf representation than the Lambert-chain corrector (the #137 point). The
blocker is purely the **multi-ellipse-coupled-through-a-bend-feasible-Venus-flyby**
structure, which the radial-crossing primitive does not provide.

### Design questions (handed to OUTSTANDING / M-ED, not improvised here)

A principled VEM free-return-style genome would need to answer:

1. **Multi-ellipse coupling.** How are the E-M and E-V ellipses' shapes co-solved
   so the Venus and Earth flybys are V_inf-continuous AND bend-feasible — i.e.
   the Venus bend is a FREE input (a flyby-propagation shooter: choose the bend at
   the sourced r_p, propagate V_inf-out), not a Lambert/crossing OUTPUT? (This is
   exactly the deferred Phase-3 n-body shooting architecture.)
2. **Same-body resonant legs.** What primitive represents the `V->V` / `E->E`
   resonant loops (resonance-ratio-parameterised full-rev ellipses) and how does
   it compose with the radial-crossing transfer arcs in one corrector?
3. **Phase closure across three bodies.** The E-M free return's Term A is a single
   relative-phase constraint; a VEM chain needs simultaneous V/E/M phase closure
   over the 12.8-yr repeat period — what is the residual vector and its rank?

These are real-design items (multi-arc-per-leg / flyby shooter), squarely the
existing M-ED front-runner — NOT a modest additive extension of `free_return.py`.
No Jones hunt was run on a free-return genome because no defensible VEM genome
exists yet; running the OLD Lambert/Lambert-chain corrector again would only
re-confirm #110/#120/#122. The Jones headline xfail therefore stays **xfail** —
tolerance NOT loosened, xfail NOT flipped (no genuine criterion was met).

## Part 3 — the two stragglers + the probe fix

### Probe wiring fix
`--probe-at-truth` now HONOURS `--genome` (the results-note caveat). New
`probe_at_truth_free_return` seeds the SOURCED ellipse `(a, e)` (from aphelion +
transit), scans t0 for the best phase, runs `free_return_correct`, and reports the
end-to-end residual — making the truth-residual≈0 visible for the free-return
genome (the lambert probe path is byte-unchanged for `--genome lambert`).

### `russell-ch4-9.353Gg2` (was CLOSE-OFF-ANCHOR) — PHASE-GRID RESOLUTION
Diagnosis: not wrong topology. The seed from sourced aphelion (2.21) + transit
(85 d) already derives E∞ 9.45 / M∞ 10.55 (sourced 9.35 / 10.52) at tof_em = 85.0 d
— truth IS representable. But this deep-aphelion, high-e (e≈0.43) row has a NARROW
t0 residual basin: a 256-point phase grid steps over it (seed-res-at-best-phase
0.54 → corrector drifts to E∞ 9.95, Δ0.60 > 0.5 → off-anchor). Denser scans:
1024 → 0.16 (E∞ Δ0.04); **4096 → 0.0048 (E∞ Δ0.09, M∞ Δ0.02) = CLOSE-AND-MATCH.**
Fix: a dense phase floor (`FR_PHASE_EPOCHS_FLOOR = 4096`) on the (cheap,
Lambert-free) free-return t0 scan in `run_row_free_return` + the probe. The row is
now CLOSE-AND-MATCH by default, AND it clears §14 V1 mechanics (closed,
V_inf-continuous arc) → promoted to **V1** (4th free-return V1, same writeback
machinery as Part 1).

### `russell-ch4-6.44Gg3` (was NO-CLOSE, res 0.2666) — GENUINELY DIFFERENT FAMILY
Diagnosis: the NO-CLOSE was ALSO a phase-grid artifact (it closes cleanly at 4096:
res 0.0017), but it closes **OFF-ANCHOR** — and that off-anchor is real:
- The seed from aphelion (1.54) + transit (262 d) gives a=1.27, e=0.213,
  tof_em=226 d → derived E∞ 3.01 / M∞ 3.06, FAR from sourced E∞ 6.44 / M∞ 3.74.
- A free-return E-M ellipse CAN match sourced (6.44, 3.74) — but only at
  a=1.225, e=0.259, **tof_em = 166 d**, NOT the 226-262 d the row's aphelion/transit
  imply. The row's sourced aphelion + transit and its sourced V∞ describe DIFFERENT
  free-return ellipses.
- The 262-d transit (vs the ~85-175 d of every matching row) + aphelion 1.54
  (barely reaching Mars, apo≈Mars sma) is a slow near-tangent transfer — a
  different arc than a fast free-return. turn_ratio 0.95, E∞ > M∞.

Verdict: a descriptor/topology mismatch, NOT a phase or seeding deficiency. The
free-return *family* exists for this V∞ (166 d ellipse), but the catalogue row's
aphelion + transit point at a different, low-V∞ arc. No promotion; stays
CLOSE-OFF-ANCHOR. (Candidate follow-up: re-derive the seed for this row from its
V∞ rather than aphelion+transit — an M-ED descriptor-interpretation item, not a
genome defect.)

### Updated counts (default `--genome free-return --model circular`)
**8 CLOSE-AND-MATCH** (9.353Gg2 joined), 3 SYMMETRIC-ONLY, 1 CLOSE-OFF-ANCHOR
(6.44Gg3). 0 NO-CLOSE. **Validation census: 6 rows V1** — Aldrin pair + four
free-return rows (5.30gGf3, 9.94Gg3, 5.75ggF3, 9.353Gg2).
