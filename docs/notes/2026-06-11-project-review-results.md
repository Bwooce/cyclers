# Full-project review — results and dispositions (2026-06-11)

**Scope:** the first full-project correctness review: all of
`src/cyclerfinder/` (core math, search correctors, data+verify layer,
performance) plus an independent numerical math-verification pass.
**Method:** 4 scoped review agents (core math / search correctors /
data+verify / performance) + a numerical math-verification agent running
8 independent probe suites — finite-difference vs STM, ∇C·f = 0 (Jacobi
gradient orthogonality), Lambert vs independently-coded Kepler, flyby
vector identities, Tisserand invariance, frame conventions, and
rotating↔inertial round-trips. Every flagged finding was adjudicated
against ground truth (re-runs, source re-transcription, independent
derivations) before any disposition.

## Headline verdict: the equations are right

Every independently-derived formula agreed with the implementation at or
near machine precision. **All confirmed defects live in solver plumbing
and conventions** — bookkeeping, sign/rotation conventions, residual
scaling, guard wiring — and every one of them lived in a stratum that no
published golden test exercised. The golden-test discipline protects
exactly the strata it covers; the defects pooled in the gaps between
goldens (refinement-loop internals, diagnostic fields, backend rotation
conventions, derivative algebra on branches the goldens never hit).

## Findings and dispositions

| Finding | Confidence | Disposition | Commit(s) |
|---|---|---|---|
| `joint_epoch_tof_close` epoch double-shift (refinement re-applied a stale grid-relative offset; no-improvement iterations re-centred up to 2× off the best point) | Confirmed (deterministic regression fails on pre-fix code) | FIXED; #181 6-row closure re-run: **6/6 UNCHANGED to every recorded digit**, V1 writeback STANDS | `58674f6` fix, `d7f0c87` adjudication |
| Vacuous lon diagnostics in the same closer (`residual_lon_deg ≡ 0`, `sc_lon_deg` duplicated `mars_lon_deg` — dressed `on_family.lon_ok` up as an independent gate) | Confirmed | FIXED — diagnostics-honesty defect only; the binding terms for this path were always the two v∞ bands | `58674f6`, `d7f0c87` |
| Data-layer guard bypasses: nested `validation.level` routed around the over-claim validator AND `apply_v0_v1/v2/v3` never checked the evidence registry (both ends open); registry drift between hand-maintained dicts undetected; duplicate-id silent last-write-wins; non-atomic catalogue writes | Confirmed | FIXED (validator + writeback helpers both check `_LEVEL_EVIDENCE`; preflight drift gate; duplicate-id `ValueError`; sibling-tmp + `os.replace` atomic writes). 448 data tests green | `091783a` |
| R_x(−i) mirror in `_InclinedCircularBackend` + the ramped-elements `_tilt` (Standish ascending-node convention requires R_x(+i); orbit normals sat 2×inc off DE440 — Venus 6.789°, Mars 3.699°) | Confirmed (doubly: independent n̂ formula + DE440 h = r×v anchor) | FIXED; post-fix normals ≤0.0004° from DE440. **Blast radius adjudicated:** the #120 3D-inclination negative STANDS (its decisive DE440 control was mirror-free); Tisserand is mirror-invariant; continuation results valid (the λ=1 step is no longer a plane flip) | `278ff1a` (search), `1d6ad1b` (core) |
| Lambert dT/dz spurious √C factor + Illinois root-solve residual overflow on multi-rev high branches (log-compression rescues the dropped branches) | Confirmed | FIXED; full suite green | `f6a0460` |
| `crosscheck_leg` shared Lambert endpoints with the primary path (consistency dressed as independence) | Confirmed | FIXED — independent endpoint re-query default-ON, poisoned-input fault test proves the teeth, gate-classification convention recorded in `validate.py` | `ba55b2e` |

## Zero retractions

Every adjudication **confirmed** the existing results — nothing was
retracted:

- **Joint-closer re-run (#195):** all 6 #181 rows reproduce the recorded
  (v∞, ToF) tuples to every recorded digit; the double-shift corrupted
  which grid points the refinement *evaluated*, never the returned
  tuple's internal consistency, and the real rows stayed in-basin. The
  approved V1 writeback stands
  (`docs/notes/2026-06-10-tof-fix-closure-results.md`, adjudication
  section).
- **Mirror blast radius:** the #120 "3D-inclination refuted" negative
  was re-examined because the mirrored backend fed it; its DE440 control
  (full 3D + eccentric, mirror-free) independently bounds the result, so
  the negative stands.
- **McConaghy Table 7.1 tail (#193/#94):** the rows-18–24 anomaly is a
  **SOURCE print defect** — Table 7.1 prints Table 7.5's dates with an
  orphaned V∞/closest-approach tail incompatible with them. Our
  transcription is character-exact (rendered page + OCR layer agree);
  the DE440-emerged tail matches Table 7.5's printed values (≤0.05 km/s,
  CA to 104–511 km on shared dates). Reproduction verdict upgraded
  PARTIAL-CONFIRMED → **CONFIRMED** (all 23 legs vs printed goldens);
  Table 7.1 rows 20–24 flagged DEFECTIVE for golden use, Table 7.5
  recorded as the sourced golden for that segment. **#94 closed.**
  Commits `53c0a92` (per-leg reproduction), `a8c0928` (adjudication).

## The false-consensus doctrine (3rd incident)

The review's probe suites confirmed a **63 s UTC/TDB epoch-conversion
offset shared between the primary path and its "independent"
cross-check** — the third false-consensus incident, after the #180
shared-ToF bug (three "independent" methods inheriting one upstream
defect) and the #197 crosscheck endpoint sharing. Agreement between
gates is only worth what they do NOT share. Doctrine, now operational:

1. **Gate-classification tiers** — every `_LEVEL_EVIDENCE` gate is
   declared CONSISTENCY (same inputs, different algorithm) or
   INDEPENDENCE (independently re-derived inputs); every promotion
   requires at least one true independence gate (convention recorded in
   `src/cyclerfinder/data/validate.py`'s module docstring).
2. **"shared with primary path:" declarations** — each cited gate names
   what it shares with the construction under test; if the line cannot
   be written, the gate is not independent.
3. **Fault injection** — poison an input and prove the gate fails
   (the #197 poisoned-input test is the template).
4. **Positive controls** — each gate demonstrated to catch a planted
   defect of the class it claims to catch.
5. **Per-interface external anchors** — every convention boundary
   (frames, time scales, node conventions) pinned to an external source
   (the DE440 orbit-normal anchor is the template).

## Open items (tracked)

- **#198** — the 63 s UTC/TDB epoch-conversion offset (probe-confirmed;
  fix in flight).
- **#212** — Earth-Moon μ double-count in `cr3bp_system()` (−1.2%;
  found independently by the Ross & Roberts-Tsoukkas mining pass,
  `docs/notes/2026-06-11-ross-roberts-tsoukkas-2025-mining.md` §7; fix
  in flight; the 2026-06-10 Earth-Moon backfill needs a re-run after.
  Saturnian/Jovian pairs expected unaffected — system-vs-planet GM
  relative error ≤ ~2e-4 — pending quantification).
- **#206** — Lambert blast-radius re-run (which historical results
  exercised the fixed dT/dz / Illinois branches).
- **#201** — performance batch (review's perf findings, none
  correctness-relevant).
- **#202** — fault-injection harness (doctrine item 3, generalised).
