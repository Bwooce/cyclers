# #433 Jupiter Galilean quasi-cycler sweep — Plan

> REQUIRED SUB-SKILL on execution: superpowers:subagent-driven-development. Design approved 2026-06-24 (parallel with #435).

**Goal:** Extend the repeated-moon "moontour" quasi-cycler sweep (the path that produced the Umbriel catalogue row) to the **Jupiter Galilean + inner-moon system** — untouched by #320 — at fine phase resolution and wider n_rev, the two levers that already paid off on Saturn (#344: a 0.032→0.0102 km/s phase-alias correction). Gate any sub-0.05 km/s survivor through the V2-moontour/QP gauntlet + ML flagger + literature_check.

**Architecture:** New `scripts/scan_433_jupiter_galilean.py` mirroring the committed `scripts/scan_320_epoch_aware_{saturn,neptune,pluto}.py` epoch-aware repeated-moon sweep, configured for the Jupiter Galilean system (Io/Europa/Ganymede/Callisto, + Amalthea inner) with a finer phase grid (48×48) and wider per-leg n_rev (0–5). Reuses the existing sweep + closure machinery; the only new content is the Jupiter config + grid. Survivors gated through the existing moontour/QP validation + ML + lit-check. Report-only until a survivor passes the gauntlet (then the Umbriel admission template applies — a SEPARATE follow-on, not this task).

**Confirmed context (from scout):**
- Galilean moons registered in `src/cyclerfinder/core/satellites.py`: Io, Europa, Ganymede, Callisto, Amalthea.
- #320 sweep template: `scripts/scan_320_epoch_aware_saturn.py` (20 cycles, 24×24 phase, n_rev∈{0,1,2,3}); outputs `data/scan_320_epoch_aware_*.jsonl`. READ this script — it is the structural template.
- Gauntlet entry points: `run_v1_qp` (`src/cyclerfinder/data/validation/v1_qp.py:289`), `run_v2_qp` (`v2_qp.py:221`); moontour V2 variant at `v2_moontour.py`; ML flagger `FalsePosFlagger().score()` (`src/cyclerfinder/ml/falsepos_flagger.py:96`); `check_literature` (`src/cyclerfinder/search/literature_check.py`).
- Risk: Galilean Laplace resonance (Io:Europa:Ganymede 1:2:4) → many phase-locked periodics; the rotation-number irrationality filter (`is_practically_irrational`, `qp_tori.py:756`, max_denominator=10 tol=1e-3) and the moontour resonance handling must screen these.

**Conventions:** main; `uv run` ruff+mypy; no Co-Authored-By; pathspec; subagents finish through commit, no self-spawned reviewers.

---

## Task 1: Jupiter Galilean sweep config + smoke

**Files:** Create `scripts/scan_433_jupiter_galilean.py`

- [ ] READ `scripts/scan_320_epoch_aware_saturn.py` end-to-end; copy its structure (the epoch-aware repeated-moon enumerator, the phase grid, the Lambert closure, the per-cell JSONL record, the offline literature corpus, `_print_progress`). Change: central body = Jupiter; moon set = [Io, Europa, Ganymede, Callisto] (+ Amalthea if it registers cleanly); enumerate closed length-3 moon cycles; **phase grid 48×48**; **per-leg n_rev ∈ (0,1,2,3,4,5)** excluding trivial all-zero. Output `data/scan_433_jupiter_galilean.jsonl`. Print per-cycle best closure + dv_band + physical-gate (max-bend) + literature status. Keep the resonance/irrationality screening that scan_320 uses.
- [ ] Smoke: run with a TINY grid (one moon cycle, 8×8 phase, n_rev (0,1,2)) to confirm it enumerates Galilean cycles, closes Lamberts, and writes records without crash. Confirm the Galilean moon states resolve (Jupiter GM + SAT441/satellites.py elements).
- [ ] ruff; commit `search/#433: Jupiter Galilean repeated-moon quasi-cycler sweep config`.

## Task 2: deliverable sweep run (controller)

- [ ] Launch `scripts/scan_433_jupiter_galilean.py` detached (full 48×48 grid, n_rev 0–5, all Galilean length-3 cycles) + harness-tracked waiter. Expect minutes-to-tens-of-minutes (O(cycles × 48² × n_rev combos) Lambert solves; cheap per cell).
- [ ] Harvest → rank cells by closure residual; identify any sub-0.05 km/s closures that pass the physical max-bend gate (>5°) and are literature-fresh.

## Task 3: gauntlet survivors + verdict + registry (controller)

- [ ] For each sub-0.05 km/s physical-gate-passing survivor: run the moontour/QP gauntlet (`run_v2_qp`/`v2_moontour`), ML flagger (`FalsePosFlagger().score()` < 0.5 else human-flag), `check_literature`. Do NOT write to catalogue here — catalogue admission of any genuine survivor is a separate follow-on (the Umbriel #339/#340 template).
- [ ] Verdict note `docs/superpowers/plans/2026-06-24-433-jupiter-galilean-verdict.md`: best closures per cycle, how many cleared the gates, lit-fresh status, resonance-lock screening counts. Honest framing: a clean "no sub-gate physical lit-fresh survivor" is a registry-grade negative for the Jupiter Galilean repeated-moon quasi-cycler region.
- [ ] Method-versioned `data/empty_regions.jsonl` entry (Jupiter Galilean was unswept by #320; this fills that gap). `uv run pytest tests/data -q -k "empty_region or registry"` before commit.
- [ ] If a genuine catalogue-able quasi_cycler survives all gates → file a follow-on admission task (do not admit inline).

## Self-review
- Reuses the proven moontour sweep + gauntlet machinery; new content = Jupiter config + finer grid (the two levers that worked on Saturn). ✓
- Resonance-lock screening retained (Galilean Laplace risk). ✓
- Report-only; survivor admission is a gated follow-on, not inline (closure discipline). ✓
- Golden/lit discipline: lit-fresh is necessary-not-sufficient; ML advisory; human + gauntlet govern. ✓
