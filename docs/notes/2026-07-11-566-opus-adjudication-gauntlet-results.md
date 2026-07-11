# #566 — Opus adjudication of the 5-representative V2→V4-strict gauntlet results

**Date**: 2026-07-11
**Task**: #566 adjudication (P0, judgment-only — the follow-up pass the #566 result paragraph
itself asked for; no new computation).
**Model**: Opus (trust-bearing verdict + catalogue-writeback recommendation, per this project's
model-tiering policy — the raw pass/fail table from #566 must not stand as the final word without
a judgment pass, matching the #561–#565 discipline for this thread).
**Inputs read directly**: `data/gauntlet_566_five_representatives.jsonl` (all 63 records, per
stage per candidate), `src/cyclerfinder/data/validation/v4_uranus_strict.py` (full module),
`docs/notes/2026-07-11-565-fable-corrected-adjudication.md`, and the #559 / #338 / #330–#332 /
#535 / #557 OUTSTANDING.md entries for the negative-baseline comparison.

---

## 0. Verdict in one paragraph

The uniform 5/5 `PASS_AS_QUASI_CYCLER` is **genuine** — it is not a gauntlet malfunction, and it
is real evidence that all six non-Miranda Uranian moon-pair directions produce short-arc
quasi-cyclic structures that survive J2 + real-ephemeris perturbation over their 3–25 day legs,
exactly as #312 does. **But the gauntlet's power to *reject* a real exact-closure member of this
class is weak-to-zero and now demonstrably so**, so a PASS should be weighted as *confirms class
membership*, **not** as independent strong existence evidence for any individual candidate. The
single decisive gap is epoch coverage: the 5 were tested at **one** epoch (2000-06-21) — which is
precisely #312's own known-favorable anchor — and #559 has already shown that #312's V4-strict
pass/fail *does* vary across epochs (≈86–90% daily PASS), with the FAILs being **confirmed
numerical artifacts** (Lambert branch-flip + DOP853 planet-crossing), not real dynamics.
**Recommendation: do NOT write any of the 5 to `data/catalogue.yaml` yet.** Gate writeback on
#567 (below): apply the already-diagnosed #560 robustness fixes, fix the audit-field bug, then run
the per-representative epoch-sensitivity scan so the family's real-ephemeris evidence matches the
provenance #312 actually earned.

---

## 1. Data verification (I did not trust the #566 summary)

Read every record in `data/gauntlet_566_five_representatives.jsonl`. The per-stage numbers hold:

| candidate (jsonl line) | V2 drift (km) | V2 closure (km/s) | V3 agree (km) | V4 agree (km) | V4-strict agree-vs-V3 (km) | chain |
|---|---|---|---|---|---|---|
| Titania-Oberon-Titania (57) | 7.7e5–8.7e5 | ~1e-14 | 5e-8 | 1.7e4–2.0e4 | 2.3e3–3.1e3 | FAIL-V2 / PASS V3+V4+V4s |
| Ariel-Umbriel-Ariel (2) | 3.8e5 | ~1e-14 | 1.6e-8 | 2.2e3 | 3.7e3 | FAIL-V2 / PASS V3+V4+V4s |
| Ariel-Titania-Ariel (12) | 3.8e5 | ~1e-14 | 1.2e-8 | 2.3e3 | 4.6e3–5.2e3 | FAIL-V2 / PASS V3+V4+V4s |
| Ariel-Oberon-Ariel (18) | 3.1e5–3.8e5 | ~1e-14 | 6–7e-8 | 1.9e3 | 2.8e3–3.4e3 | FAIL-V2 / PASS V3+V4+V4s |
| Umbriel-Titania-Umbriel (26) | 2.9e5–5.3e5 | ~1e-14 | 1.1e-8 | 4.1e3–6.9e3 | 2.5e3–6.0e3 | FAIL-V2 / PASS V3+V4+V4s |

Each is run at n_cycles ∈ {3, 5, 10}; all three sub-runs agree per candidate. The V2 FAIL +
machine-precision closure = the `FAIL_QUASI_BOUNDED` admission pattern #330 established (exact
closure by construction, drift overflow = *defining property of the class*, not a rejection). The
5/5 uniformity is exactly as reported. **Confirmed.**

## 2. The audit-field bug — confirmed cosmetic for the verdict, but a provenance defect at writeback

Verified in `v4_uranus_strict.py` directly. Lines **575–578** hardcode
`_moon_state_spice("Umbriel", et_launch)` and `_moon_state_spice("Oberon", et_launch)` to compute
the four recorded `eccentricity_used_* / inclination_used_*` fields — **regardless of the
candidate's actual sequence**. The jsonl confirms it empirically: every one of the 5 candidates
records the identical `e_u=0.0041, e_o=0.0016`, including Ariel-Titania-Ariel (which involves
neither Umbriel nor Oberon).

The physics path is genuinely independent of these fields: `_cycle_v4_strict` (lines 375–441)
SPICE-samples each Lambert endpoint via `_moon_state_spice(moon, ...)` over the **actual**
`sequence`, and the perturber set is the full 5-moon `URANIAN_PERTURBER_MOONS` default — neither
reads the four audit scalars. `passes_v4_strict` (lines 673–679) is a function of
`drift_agreement_vs_v3`, leg convergence, and `bounded_drift_survives` only. **So the pass/fail
verdict is unaffected — the characterization "physics generic, audit fields cosmetic" is true.**

Where it is *not* merely cosmetic: **at writeback**. #312's provenance was frozen into
`data/silver_327_v*_verdicts.jsonl` and pinned by `tests/verify/test_silver_327_v4_strict_passes.py`.
If these 5 verdicts are frozen the same way, four of five rows would bake a **wrong** eccentricity/
inclination into their permanent provenance record — a recorded-value defect of exactly the kind
`[[feedback_digest_not_adoption]]` warns about (a sourced/recorded number silently diverging from
what the row actually is). **Non-blocking for *this* adjudication; MUST be fixed before any
writeback that freezes these verdicts.** (Fixing it is in scope for #560/#567, not this pass.)

## 3. The core question: genuine broad validation, or a non-discriminating gauntlet?

**Both readings are partly correct, and the negative baseline resolves the tension.** Decompose
what each stage actually tests on this candidate class:

- **V2 closure ~1e-14 is tautological.** The #563 enumeration *constructs* exact symmetric V∞-
  magnitude closures; the V2 closure residual re-measuring ~0 is guaranteed, not evidence.
- **V2 drift overflow is the class definition.** Every quasi_cycler in this family drifts; FAILing
  the 50,000 km floor is what *makes* them quasi_cyclers rather than strict cyclers. Uniform by
  construction.
- **V3 agreement ~1e-8 km is a numerical self-consistency check**, not independent physics — it
  confirms REBOUND IAS15 and the analytic Kepler propagator agree on the *same idealized model*.
- **V4 / V4-strict agreement (2–7e3 km vs a 50,000 km floor)** confirms that J2 + real-ephemeris
  perturbations don't materially bend the **short** Lambert arcs over the window. Over 3–25 day
  legs deep in Uranus's gravity well, third-body/J2 perturbations *are* physically small — so
  passing is the expected outcome for **every** member, not a discriminating filter.

So option (b) from the task framing is substantially true: **the chain has little power to
separate "real" from "artifact" closures within this specific population** — everything that is an
exact short-arc symmetric closure at Uranian moon-pair energies passes. But option (a) is *also*
true and is the more important read: the reason everything passes is **benign** — the short-arc
physics genuinely is benign for all six directions, so a PASS is a faithful confirmation of
quasi_cycler class membership, identical in kind to what #312 itself earned.

**The negative baseline confirms this is the right synthesis, not hand-waving.** What does a REAL
negative look like in this stack?

- **#535 Earth Hill-sphere corridor**: a ~15,000 km knife-edge that **collapsed entirely** under
  ER3BP eccentricity (e=0.0167 vs idealized e=0). A genuine physical rejection — the structure did
  not survive a more faithful model. None of the 5 do anything like this; their V4-strict-vs-V3
  agreement is 2–6e3 km, three orders under the wall #535 hit.
- **#559 daily-epoch scan of #312 itself**: ≈86–90% daily PASS, and the ≈10–14% FAILs were
  **diagnostically confirmed to be numerical artifacts** — a Lambert branch-selection flip with no
  continuity tracking (`_cycle_v4_strict` ~line 418; terminal miss jumps 23,448→1,340 km between
  adjacent *hours*) and a DOP853 stiff-death on non-physical near-parabolic planet-crossing arcs
  (perijove *inside* Uranus, no collision guard). **Every observed V4-strict rejection of a real
  #312-class closure to date traces to a numerical artifact, not to real dynamics.**

That is the crux: the gauntlet's only demonstrated "rejections" in this exact population are
artifacts. So a clean-epoch PASS is trustworthy as class-membership confirmation, **and** the
chain has essentially zero demonstrated power to physically reject a real member — precisely
because there is nothing physical to reject over these short arcs. The 5/5 is genuine *and* the
gauntlet is non-discriminating, for the same benign reason.

## 4. Is single-epoch V4-strict (2000-06-21 only) sufficient confidence? No.

This is the one place I diverge from a "just write them all in" reading, and #559 makes the case
concretely rather than on general caution:

1. **The tested epoch is #312's own favorable anchor.** 2000-06-21 (DOY 172) is the epoch #338
   selected *because* #312 passes there, and #559 re-confirmed it PASSes in both 2000 and 2030
   (~14k km drift). Testing the 5 new candidates at exactly this epoch is testing at the single
   point already known to be favorable for the sibling. It is the epoch most likely to pass, for
   all of them, for reasons that may be partly phase-configuration luck rather than intrinsic
   robustness.
2. **#312 is known to FAIL V4-strict at other epochs.** #338 found 2000-01-15 FAILs (91,000 km
   drift); #559 found ≈10–14% of daily epochs FAIL. Even granting those specific FAILs are
   artifacts, the point stands: a *single* epoch tells us nothing about whether a given candidate
   sits in a wide tolerant band or on an artifact-prone phase, and the 5 candidates have **no**
   epoch-sensitivity data at all.
3. **A naive per-candidate daily scan run *today* would be untrustworthy.** The two artifact
   mechanisms #559 diagnosed live in `v4_uranus_strict.py` and are **candidate-agnostic** — they
   would contaminate the 5 candidates' epoch sweeps exactly as they contaminated #312's, producing
   the same ~10–14% artifact-driven FAIL spikes and telling us nothing real. The #560 fixes
   (Lambert-branch continuity + perijove/collision guard) are the prerequisite for a *meaningful*
   epoch scan. This is the `[[feedback_isolated_sweep_flips_suspect_artifact]]` lesson applied
   forward: don't run the scan until the known artifact generators are fixed, or you manufacture
   the same false negatives.

So: a single untested-for-sensitivity epoch, at the one epoch known to be favorable, under code
with two confirmed artifact generators, is **not** sufficient confidence to grant these 5 the V4
validation_level #312 legitimately holds. #312 earned V4 via the #338 annual sweep + the #559
daily characterization; the 5 have neither.

## 5. Catalogue-writeback recommendation

**Hold. Do not write any of the 5 to `data/catalogue.yaml` in the current evidence state, and do
not upgrade #312's row on the strength of #566 alone.** The writeback these results *point
toward* is real and worth doing — but doing it now, at V4 parity with #312 on one favorable epoch
and under the audit-field/artifact defects, would over-claim relative to #312's actual provenance
and would freeze a wrong eccentricity/inclination field into four of the rows.

**When the gate (#567) clears, the writeback should follow #564 §5's reframing**: admit #312 and
the validated representatives as **"the first-documented members of a 30-member Uranian
symmetric-closure quasi_cycler family"** (enumerated-30 idealized-exact kept distinct from
validated-N real-ephemeris), N = the number that clear the *fixed-code, epoch-robust* gate — not a
blanket "all 30 validated." The Titania-Oberon representative (line-57 / row 23) must cite
Canales-Howell-Fantino (arXiv:2110.03683, 2308.10029) and Kumar 2025 (**arXiv:2509.03655 §6.2** —
the corrected ID from #565 §5) as adjacent-but-not-colliding prior art; the other four carry only
the #328 Uranian baseline. Each row's prose should state honestly that the gauntlet confirms
*class membership* (short-arc structure survives real ephemeris at the validated epoch band), not
independent per-candidate existence — the exact closure is a construction property (#563), and the
family shares one dynamical mechanism.

## 6. Recommended next task

**#567** — Trustworthy epoch-robustness gate for the #566 representatives + #312, as the
catalogue-writeback prerequisite. Scope: (1) apply the already-fully-diagnosed #560 robustness
fixes to `v4_uranus_strict.py` (Lambert branch-selection continuity tracking; perijove/collision
guard that distinguishes a physical FAIL from a DOP853 stiff-death on a planet-crossing arc); (2)
fix the hardcoded Umbriel/Oberon audit-field sampling (§2) so recorded e/i track the actual
sequence — cheap, and prevents freezing wrong provenance at writeback; (3) under the fixed code,
run the #338-style annual + #559-style daily epoch-sensitivity scan on all 5 representatives and
re-confirm #312, reporting each candidate's true PASS band width (wide-tolerant vs narrow) rather
than a single point; (4) THEN the writeback decision, at a validation_level justified by the
band-width evidence, per §5. Full detail in the OUTSTANDING.md #567 entry.

---

## 7. Summary

1. **5/5 PASS is genuine, not a malfunction** — all six pair-directions produce short-arc
   quasi-cyclic structures that survive real-ephemeris perturbation over their legs, same as #312.
2. **The gauntlet is non-discriminating for this class, for a benign reason** — the short-arc
   physics is genuinely benign for all members, so a PASS confirms *class membership*, not
   independent per-candidate existence (the exact closure is a #563 construction property). The
   negative baseline (#535 real collapse; #559 artifact-only #312 FAILs) confirms this synthesis.
3. **Audit-field bug: cosmetic for the verdict, provenance defect at writeback** — must be fixed
   before freezing any of these verdicts into a catalogue row.
4. **Single-epoch V4-strict at #312's favorable anchor is insufficient** — no epoch-sensitivity
   data exists for the 5, #312 is known to vary across epochs, and a scan run under today's code
   would reproduce #559's confirmed artifacts.
5. **Recommendation: HOLD writeback; gate on #567** (apply #560 fixes + audit-field fix + per-
   representative epoch scan), then admit #312 + the epoch-robust survivors as the first-documented
   members of the 30-member symmetric-closure family per #564 §5's reframing.
