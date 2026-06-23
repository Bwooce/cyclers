# M7 Phase 1 results — per-row real-ephemeris maintenance-ΔV (horizon TCM) is built

**2026-06-23 (#423).** Phase 1 of the M7 build (plan: `2026-06-23-m7-implementation-plan.md`)
is complete: `real_closure.RealClosureResult.horizon_tcm_mps` (hardcoded `0.0` since M6b)
is now a real measurement under the opt-in `compute_tcm` flag, golden-validated on the
real S1L1 cycler.

## What was built

- **`nbody/maintenance_shoot.target_leg`** — n-body fixed-time position-targeting Newton
  on the propagator's co-integrated STM (`d r_f/d v_0 = STM[0:3,3:6]`). Solves for the
  departure velocity whose (perturbed) propagation hits the next encounter planet — the
  n-body analogue of a fixed-time Lambert / the flyby B-plane targeting as position
  targeting. Golden: Sun-only ≡ two-body `lambert()` to mm/s; Mars-perturbed sub-km.
- **`nbody/maintenance_shoot.continuous_maintenance_chain`** — walks a cycler's encounter
  sequence as ONE continuous trajectory, targets each leg, and bills each interior flyby's
  maintenance ΔV via `core.flyby.flyby_dv_for` (the part a free ballistic flyby cannot
  bridge). Honest divergence → `inf` (row stays V0, never forced).
- **`real_closure.verify_real_closure(compute_tcm=True)`** + `_compute_horizon_tcm` — builds
  the N-cycle node sequence (incl. the inter-cycle home flyby), seeds per-leg from the
  constructed cycler's rev-correct `leg.v_depart`, runs the chain, populates
  `horizon_tcm_mps`/`per_cycle_tcm_mps`, which feed `v3_class_split_verdict` (#424).

## Three design decisions (with evidence)

1. **Leg-by-leg position-targeting, not global B-plane multiple-shooting.** Lighter, reuses
   the propagator STM + `flyby_dv` + `real_closure`'s Cycler builder; does not touch
   `shooter.py`. The infra audit confirmed soundness.
2. **Per-leg perturbers = system bodies − the leg's endpoint flyby bodies.** Naive continuous
   integration of an endpoint body from its own node centre diverges (the spacecraft starts
   in the softened core) = the patched-conic handoff artifact, not fuel. Endpoint exclusion
   removes it while keeping every genuine non-endpoint third-body perturbation — so M7 works
   for moon tours / Venus-flyby chains, not just Earth-Mars (where it reduces to Sun-cruise).
   Validated: naive diverges; default converges; a non-endpoint Jupiter perturber shifts the
   TCM (real signal).
3. **Per-leg seeds from the rev-correct sourced/constructed departure velocity.** Decisive:
   seeding a multi-rev resonant cycler from a single-rev Lambert lands the Newton on a
   high-energy single-rev transfer — S1L1 read **530,000 m/s** and the powered rows diverged.
   With sourced-v∞ seeds S1L1 reads **40.1 m/s**, reproducing the corrected #169 proxy.

## Golden / validation

- S1L1 full 22-node App-C sequence, sourced-seeded, 7-cycle: **~40 m/s**, all 21 legs
  converge — matches the corrected Sun-only proxy (~40.2 m/s, post-#198) to abs 5 m/s.
  M7 position-targeting and the proxy's sourced-direction accounting measure the same
  ballistic physics on a real resonant cycler. Under the 120 m/s V3 budget → S1L1
  essentially_ballistic (unchanged).
- Wiring test: `_compute_horizon_tcm` on a real-eph E-M-E cycler returns a finite 2-cycle
  TCM (incl. the home flyby) and drives `v3_class_split_verdict` programmatically — the
  manual #175 convention is retired into code.

## Scope / limits

- **Heliocentric lane only** (μ_Sun, Sun-central propagator) — exactly the `real_closure`
  cycler lane (the ~215 V0 Earth-Mars SnLm census + heliocentric multi-planet cyclers).
  Planet-central moon-tour M7 (Saturn/Uranus central body) is a separate generalisation;
  those rows validate via their own V3/V4 lanes (`v2_moontour`/`v3_3d`/`v4_uranus`), not
  touched here.
- **Phases 2-3 remain** (not started): run M7 on the banded + V1+ rows, then a catalogue-wide
  detached/checkpointed batch (days-scale) = the V0→V3 mass-promotion / validation-ceiling
  lift. Expect a residue of unshootable (off-anchor / high-V∞) rows that stay V0 honestly.

## Phase 2 coverage scan (2026-06-23) — M7 is DATA-gated, not compute-gated

Ran `scripts/m7_phase2_coverage.py` (`verify_real_closure(compute_tcm=True, n_cycles=2)`)
over all 318 cycler-type rows. Result: **measured (finite TCM) = 8, diverged = 5,
skip (can't construct / no window / error) = 305.**

The decisive finding: **305/318 rows cannot even construct a real-eph cycler** — their
full per-leg node sequences are not tabulated (incomplete leg data, the "M7's
catalogue-completion concern" already noted in `real_closure.EXPECTED_SKIPS`). And the
8 that *do* return a finite number are **non-trustworthy**: 11k–126k m/s, all
`v3-real-closure-fail` (e.g. `rall-1970-m4-1` = 126,417 m/s). Those are
`real_closure`'s AUTO-constructed cyclers run WITHOUT sourced per-leg v∞ seeds — the
same single-rev / wrong-basin pathology that gave S1L1 530,000 m/s before seeding.

So a meaningful catalogue-wide M7 measurement is gated on TWO data inputs, not compute:
1. **Full-cycle leg sequences** (305 rows lack them) — the catalogue-completion task.
2. **Sourced per-leg v∞ seeds** so the targeting stays in the rev-correct basin — only
   the App-C-block rows (S1L1, #188, #192) carry these today.

The only trustworthy computed M7 numbers come from the sourced App-C blocks:
- **S1L1**: 40 m/s/7cyc (ballistic; the validated golden).
- **#188 / #192**: honestly diverge (powered; not ballistically maintainable — rest on
  their published 420 / 1678 m/s budget).

Note S1L1's computed continuous TCM (40 m/s) lands in the **low_maintenance** band
(10–300), while its *sourced* band (#417) is **essentially_ballistic** (<10). Not a
contradiction — sourced = best-window ΔV, M7 = continuous-from-one-seed TCM (different
bases, exactly the #424 rationale). This is a sourced-vs-measured band MISMATCH to
surface for review, not silently overwrite.

**Phase 3 (catalogue-wide V0→V3) is therefore blocked on leg-data + seed reconstruction,
not on the shoot.** It should be sequenced after (or merged with) the full-cycle
leg-completion effort; running the shoot blindly now produces garbage on auto-constructed
rows. Recommend NOT mass-promoting; instead drive M7 from sourced node+seed data per row
as that data is reconstructed.

## Option C reproduction attempt (2026-06-23) — STOP & REASSESS

Directive: "scope and execute option C; if we don't reproduce, stop and reassess."
Attempted reconstruction of the Russell-2006 strictly-ballistic census rows (the V0
ceiling) → M7. **Verified blocker, three independent confirmations:**

1. **No per-member data in the catalogue.** e.g. `russell-2006-117-5.225Ggg3`:
   `invariants` all null, `transit_times_days` [], `free_return_arcs` [],
   `trajectory.segments` []. Sole sourced content = the descriptor `5.225Ggg3` + the
   claim "one of 9 cyclers < 1 m/s". Exactly the validation-ceiling publication gap:
   the per-member reproducible state was never printed.
2. **Current reconstruction tooling returns None** for all 9 strictly-ballistic rows
   (`seed_dsm_chain_from_descriptor` → None; only S1L1 / mcconaghy-2006-em-k2 returns a
   chain — and it already reproduces, 40 m/s).
3. **Reconstructing the specific member from the bare descriptor hits the #388
   family-selection wall** (documented, `2026-06-08-multiarc-basin-selection-results.md`):
   the reconstruction relaxes to off-anchor families (emerged V∞ ~9–16 km/s vs published
   ~5.2) — the published cycler is NOT uniquely recovered. The descriptor underdetermines
   the family.

**Conclusion:** M7 (the measurement) is sound and validated, but we **cannot reproduce
the per-member real-eph state of the Russell census** with current inputs/tools — it is
gated by (a) a publication gap (data never printed) and (b) the #388 family-selection
wall on descriptor-only reconstruction. Reproducible today: S1L1 (40 m/s ballistic) +
the 2 App-C powered rows (diverge → published budget). The other ~215 stay V0 honestly.

**Reassessment — the tractable route is DATA, the hard route is ALGORITHM:**
- **(A) Data path (the S1L1 recipe, tractable):** Russell 2004 dissertation Appendix C
  (pp.201–245, readable, 43 dense pages) holds per-leg "DATA NECESSARY TO REPRODUCE"
  blocks, parent-indexed — but only #188/#192 are transcribed. If the strictly-ballistic
  parents (#54/#117/#177/…) have App-C blocks (gated by the Table5↔App-C parent-number
  bridge, the #170 caveat), transcribing them makes those rows reconstructible + M7-
  measurable exactly as S1L1 was. Bounded corpus transcription, not research.
- **(B) Algorithm path (#388 frontier, open):** a global/family-targeted reconstruction
  that recovers the specific published member from the descriptor. Multi-week; #388
  characterized this to the bottom as hard.
- **(C) Accept the ceiling:** census rows are V0 by an irreducible publication gap;
  reproduce only what has sourced data.

Discovery stays gated (the "reproduce reliably first" bar is NOT met). Recommend (A).

## S1L1 published-number reproduction attempt (2026-06-23) — epoch is NOT the gap

Tried to close S1L1's M7 (40 m/s) vs published essentially_ballistic (<10 m/s) by
launch-epoch/phasing optimization. **Verified dead end:** an epoch sweep of ±780 d
(one E-M synodic) in 60 d steps converges ONLY at the App-C design epoch (40.1 m/s);
EVERY other shift diverges. A cycler's phasing is locked to its synodic alignment, so
shifting the launch date breaks the fixed-ToF resonant structure and no ballistic chain
exists. The App-C epoch is already optimal — epoch is not the 40-vs-<10 gap.

**So the gap is the optimization METHOD, not phasing.** M7 does **flyby-only**
maintenance at **exact** node positions (the most constrained model) → it is a
CONSERVATIVE UPPER BOUND on the true maintenance. Russell's published <10 m/s comes
from his **joint trajectory optimizer** (free mid-leg DSMs, corridor relaxation, jointly
minimized). M7's 40 is consistent with (an upper bound on) <10 — it confirms
ballistic-CLASS / V3-pass (≪120 m/s budget) but does NOT reproduce the published
sub-tier. Reproducing the published *number* requires re-implementing the joint
optimizer (free DSMs) — the #388 full-shooter frontier (compute-heavy, previously
walled), or at least extending the M7 chain with a per-leg mid-leg-DSM degree of freedom
(`search/dsm_leg.py` exists) and joint minimization.

**Net:** no multi-leg cycler is reproduced to its published maintenance tier. M7 gives a
sound conservative bound (class, not tier). Closing to the published number is an
optimizer build, not a tweak — STOP and reassess per directive.

## Gate status

Per `project_dvband_validation_coupling_gate`: the gate to resume discovery+validation
campaigns was "land #424 + ≥ M7 Phase 1". Both are now done — the gate is **cleared**.
Phases 2-3 are the ceiling-raise, not the gate.
