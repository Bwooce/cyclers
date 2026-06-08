# GMAT V4 with B-plane flyby targeting — generator → parser → manual protocol → writeback gate (#171)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:executing-plans` or
> `superpowers:subagent-driven-development`. Checkbox steps; strict TDD (write
> failing test → run **red** → minimal impl → run **green** → commit). Work on
> `main` — **do NOT branch** (project rule). uv-managed venv (no pip). Lint/type
> gate before **every** commit:
> `uv run ruff check .` · `uv run ruff format --check .` · `uv run mypy src tests`.
> Fast suite: `uv run pytest -m "not slow"`. **None of these tests are `slow`** —
> they are all string-templating / pure-geometry / fixture-parse, runnable without
> GMAT installed.
>
> This plan is the task-level expansion of the **approved** design
> `docs/notes/2026-06-08-gmat-v4-design.md` — read it first (§0 the GMAT-not-installed
> finding, §3 the B-plane kernel, §4 the V4 gate). It **extends, does not duplicate**,
> the GMAT V4 lane already designed in
> `docs/superpowers/plans/2026-06-06-nbody-harness.md` Phase D (Tasks D.0–D.3:
> `scripts/gmat_v4_aldrin.py`, `scripts/parse_gmat_report.py`, the 5% Aldrin band,
> Q7 manual out-of-CI). Where this plan and the nbody Phase D could diverge, the
> nbody-harness design + its Approval win; this plan only **adds** the B-plane
> flyby-targeting layer (#169) and **generalises** the single-row lane to S1L1.
>
> **THE GATING DEPENDENCY (binding, top of mind):** GMAT is **NOT installed in this
> environment** (`which gmat` → exit 1; design §0). Every task below is buildable
> and unit-testable **WITHOUT GMAT**. The GMAT *run* and the V4 *writeback* are
> **manual human steps** documented in the run-book (Phase 4) and explicitly OUT of
> CI (Q7). **No task in this plan installs GMAT, runs GMAT, or writes a V4 row.**
>
> **GOLDEN/HONESTY (binding).** The only SOURCED anchors are (a) each row's published
> v∞ nodes and (b) the Jones AAS 17-577 continuity tolerance (1e-3 km / 1e-6 km/s,
> §2.5). Our reference ΔVs — Aldrin **2.9138 km/s**, S1L1 **62 m/s / 7 cycles** —
> are OUR computed values; GMAT is their **independent external check**. They are the
> *reference under check*, NEVER an EXPECTED-from-source assertion. The B-plane
> geometry test is a **self-consistency round-trip** (Jones tabulates no worked
> example), not a golden. **No catalogue writeback in this plan.**

---

## Goal

Build the buildable half of the canonical high-fidelity V4 lane for **confirmed**
cyclers (Aldrin #134 first, then S1L1 `russell-ch4-4.991gG2`): a generator that
hands GMAT a confirmed row's seed + flyby sequence + **per-flyby Jones B-plane
targets at flyby periapse**, a parser whose V4 predicate is **NLP convergence AND
maintenance ΔV within the declared band**, and a **manual run protocol** with the
GMAT-install prerequisite flagged. The B-plane targeting layer is the piece #169
proved is missing (the naive patched-conic→continuous handoff diverges; ~122 km/s is
an artifact, not fuel; real maintenance needs B-plane-targeted flyby station-keeping).

Outcome of executing this plan: the generator + parser + B-plane kernel + run-book,
all green without GMAT. The V4 **verdict** and any `_LEVEL_EVIDENCE` writeback are a
downstream human step after a real GMAT run — explicitly NOT in this plan.

---

## VERIFY-FIRST (live state to re-read before any edit)

The design cites these on 2026-06-08; re-read before editing in case of drift:

1. **GMAT availability** — re-run `which gmat GMAT GMAT_console` and the `find`
   probe (design §0). If GMAT has appeared, the run-book's "manual" framing still
   stands (Q7: GMAT stays out of CI regardless), but note it in the run-book.
2. **Phase D prior art** — `docs/superpowers/plans/2026-06-06-nbody-harness.md`
   lines 1234–1354 (Phase D Tasks D.0–D.3). If `scripts/gmat_v4_aldrin.py` /
   `scripts/parse_gmat_report.py` already exist (from an nbody-harness Phase-D run),
   **extend them in place** rather than recreating — this plan adds the B-plane
   block and the S1L1 generalisation to whatever Phase D landed. As of authoring,
   `ls scripts/*gmat* scripts/*parse*` → none exist yet.
3. **Jones B-plane kernel** — `docs/notes/2026-06-07-jones-aas17-577-method-deepdive.md`
   §2.2 Eqs.1-5 (Ŝ,T̂,R̂, δ, r_p, θ_B). `core/flyby.py` already has turn-angle /
   `max_bend` / `r_p` geometry — re-read `flyby.py` and **reuse** its (μ, rp_min)
   registry resolution (the way `flyby_dv_for` does) rather than re-literal-ing
   constants.
4. **S1L1 seed** — `src/cyclerfinder/search/s1l1_corrected.py`: `APPC_LEGS`,
   `APPC_EPOCH_DAYS`, `continuous_chain(start_leg_no=2)` (leg-2 Earth departure
   2026-12-15, v∞=(−2.278, 5.322, 0.574) km/s). The per-Mars-flyby (v∞⁻, v∞⁺) pairs
   are the consecutive App-C nodes. **Read node vectors from `APPC_LEGS`**, never a
   self-computed attribute.
5. **Aldrin reference** — 2.9138 km/s per-synodic maintenance ΔV
   (`docs/notes/2026-06-07-aldrin-continuation-v3-evidence.md`,
   `maintain.py`/`bvp.py` powered-periodic lock). OUR value, reference-under-check.
6. **V4 writeback gate** — `src/cyclerfinder/data/validate.py` `_LEVEL_EVIDENCE`
   (`:480`) + `validate_validation_level` over-claim guard (`:661`); the existing
   S1L1 V3 entry (`:615`, ending "canonical V4 (GMAT) remains OPEN"). The writeback
   shape a V4 pass would take — **not written by this plan**.

---

## Phase 0 — B-plane targeting kernel (LOW; pure geometry, no GMAT)

> The new piece Phase D lacks. Jones AAS 17-577 Eqs.1-5 as a pure function from
> `(v∞⁻, v∞⁺, μ_body, rp_min)` to the GMAT B-plane goal `(B·R, B·T, r_p, θ_B)`.
> Self-consistency tested (no worked Jones example → not a golden). New module
> `scripts/gmat_bplane.py` (script-side helper, mirrors the out-of-CI script lane;
> NOT a `src/` production module — it serves the manual GMAT generator only).

- [ ] **Task 0.1** — `test_bplane_frame_orthonormal` (RED→GREEN). CONSTRUCT a
  `v∞⁻` and assert `bplane_frame(vinf_minus)` returns `(Ŝ,T̂,R̂)` that are
  unit-length and mutually orthogonal (`Ŝ = v̂∞⁻`, `T̂ = (Ŝ×k̂)/‖·‖`, `R̂ = Ŝ×T̂`,
  `k̂=(0,0,1)` — Jones Eq.4). Reference = the orthonormality identity, not an
  evaluator output. Label: mechanics.
- [ ] **Task 0.2** — `test_bplane_target_roundtrip` (RED→GREEN). CONSTRUCT a
  feasible `(v∞⁻, v∞⁺)` pair (equal magnitude, turn δ strictly inside
  `max_bend(μ_mars, rp_min, |v∞|)`). Compute `bplane_target(...)` → `(B·R, B·T,
  r_p, θ_B)` via Jones Eqs.1,2,5 and `|B| = r_p·√(1+2μ/(r_p v∞⁻²))`. Assert the
  recovered B-vector reproduces the intended turn (the projection of `v̂∞⁺` onto the
  B-plane lies along `−B̂` to < 1e-9; Jones' "the flyby bends v∞⁺ such that its
  B-plane projection is along −B"). Self-consistency, NOT a golden. Label: mechanics.
- [ ] **Task 0.3** — `test_bplane_charges_infeasible_turn` (RED→GREEN). CONSTRUCT a
  turn BEYOND `max_bend`; assert `bplane_target` flags it infeasible (a subsurface /
  out-of-window `r_p`, surfaced for the generator to reject or to add a TCM). Reuse
  `core/flyby.py`'s feasibility, do not re-derive. Label: mechanics.
- [ ] **Task 0.4** — Lint/type gate clean. **Commit (code+test):**
  `scripts/gmat_bplane: Jones B-plane flyby-target kernel (Ŝ,T̂,R̂,θ_B,B·R,B·T) (GMAT V4 #171)`.

## Phase 1 — Generator: Aldrin powered-periodic script (LOW–MED; no GMAT)

> Extends nbody Phase D Task D.0. Beeson Algorithm-1 two-section skeleton (init +
> mission sequence with one flyby Target). Single Mars flyby; recovers 2.9138 km/s.
> New/extended `scripts/gmat_v4_aldrin.py`; test `tests/scripts/test_gmat_aldrin_gen.py`.

- [ ] **Task 1.1** — `test_aldrin_script_force_model_and_epoch` (RED→GREEN). Call
  `generate_aldrin_script(out, epoch_iso="2030-01-01T00:00:00")`; assert the text
  has `ForceModel`, `Sun`, `Earth`, `Mars`, `Jupiter`, and the epoch (mirrors nbody
  D.0). Label: scripts.
- [ ] **Task 1.2** — `test_aldrin_script_has_flyby_target_block` (RED→GREEN). Assert
  the script contains a `Target`/`Achieve` block with `BdotR` and `BdotT` achieve
  goals (from Phase 0's `bplane_target` for the Aldrin Mars flyby) and the BC
  propagated to `Mars.Periapsis` (Beeson: BC at flyby periapse, gravity in EOM).
  Assert the varied TCM components are present (the maintenance impulse). The numeric
  B-plane goals come from Phase 0 applied to the Aldrin (v∞⁻, v∞⁺); the row's v∞
  nodes are the sourced input. Label: scripts.
- [ ] **Task 1.3** — `test_aldrin_script_provides_initial_guess` (RED→GREEN). Assert
  the "Provide an Initial Guess" slot (Beeson Algorithm-1) is filled with the seed
  departure state + epoch (not empty). Label: scripts.
- [ ] **Task 1.4** — Lint/type gate clean. **Commit (code+test):**
  `scripts/gmat_v4_aldrin: Aldrin powered-periodic generator + B-plane flyby Target (GMAT V4 #171)`.

## Phase 2 — Generator: S1L1 flyby-station-keep chain (MED; no GMAT)

> Same templating, loop unrolled over the 7 App-C Mars flybys (Beeson: "completely
> expanded in script form"). `mode="flyby-station-keep"`. New
> `scripts/gmat_v4_s1l1.py` (or a `row`/`mode` parameter on a shared generator —
> implementer's call, but keep Aldrin's tests green); test
> `tests/scripts/test_gmat_s1l1_gen.py`.

- [ ] **Task 2.1** — `test_s1l1_script_seeds_appc_leg2` (RED→GREEN). Assert the
  generated script's initial guess carries the App-C leg-2 Earth departure
  (2026-12-15, v∞=(−2.278, 5.322, 0.574) km/s) read from `s1l1_corrected.APPC_LEGS`,
  NOT a hardcoded literal in the test (import the constant). Label: scripts.
- [ ] **Task 2.2** — `test_s1l1_script_has_per_mars_flyby_targets` (RED→GREEN).
  Assert the script contains **one B-plane `Target`/`Achieve` block per Mars flyby**
  (7 for the 7 App-C Mars encounters), each with `BdotR`/`BdotT` goals from Phase 0
  applied to that node's (v∞⁻, v∞⁺), BC at `Mars.Periapsis`. Count the blocks ==
  number of Mars nodes in `APPC_LEGS`. Label: scripts.
- [ ] **Task 2.3** — `test_s1l1_generalises_from_row_descriptor` (RED→GREEN). Assert
  the generator takes a row descriptor (id, sequence, nodes, mode) so Aldrin and
  S1L1 flow through the same entry point with different `mode` — the generalisation
  contract (design §2). Label: scripts.
- [ ] **Task 2.4** — Lint/type gate clean. **Commit (code+test):**
  `scripts/gmat_v4_s1l1: S1L1 per-Mars-flyby B-plane station-keep generator (GMAT V4 #171)`.

## Phase 3 — Parser + two-part V4 predicate (LOW; fixture strings, no GMAT)

> Extends nbody Phase D Task D.1. The predicate is **NLP convergence (Jones 1e-3 km
> / 1e-6 km/s, sourced) AND maintenance ΔV within ±5% of our reference (where a
> reference exists)**. `scripts/parse_gmat_report.py`; test
> `tests/scripts/test_gmat_parse.py`.

- [ ] **Task 3.1** — `test_parse_maintenance_dv_sums_tcm` (RED→GREEN). Against a
  FIXTURE report string with multiple `Maneuver.TotalDV = ...` lines, assert
  `parse_maintenance_dv` returns their per-cycle sum. Label: scripts.
- [ ] **Task 3.2** — `test_parse_convergence` (RED→GREEN). Against fixture strings,
  assert `parse_convergence` returns True iff every Target/Optimize block reported
  NLP convergence (and False on a "did not converge" line). Beeson's primary
  acceptance = convergence. Label: scripts.
- [ ] **Task 3.3** — `test_v4_pass_two_part_predicate` (RED→GREEN). Assert
  `v4_pass(2.95, ref_dv=2.9138, converged=True)` is True (within 5%, converged);
  `v4_pass(2.95, ref_dv=2.9138, converged=False)` is False (not converged);
  `v4_pass(3.20, ref_dv=2.9138, converged=True)` is False (outside the
  2.768–3.060 band). Docstring states 2.9138 is OUR value under external check, the
  ±5% band declared up front never back-fit. Label: scripts.
- [ ] **Task 3.4** — `test_v4_pass_convergence_only_when_no_reference` (RED→GREEN).
  For the S1L1 Mars-perturbed continuous arm (no prior reference ΔV, design §4/§6),
  assert `v4_pass(produced_dv, ref_dv=None, converged=True)` is True (convergence-
  only) and records `produced_dv` as the figure; `converged=False` → False. You
  cannot band a number you do not yet have. Label: scripts.
- [ ] **Task 3.5** — Lint/type gate clean. **Commit (code+test):**
  `scripts/parse_gmat_report: convergence + ΔV-band V4 predicate, convergence-only fallback (GMAT V4 #171)`.

## Phase 4 — Manual run protocol + writeback gate doc (LOW writing; no GMAT)

> Run-book documenting the MANUAL out-of-CI procedure (Q7). New
> `docs/notes/2026-06-08-gmat-v4-runbook.md`. **No CI test invokes GMAT.**

- [ ] **Task 4.1** — Write the run-book. Sections:
  1. **PREREQUISITE — install GMAT** (pin a version, e.g. R2022a; record the install
     path; this is the gating dependency — GMAT is NOT in this environment).
  2. **Generate** the script: `uv run python scripts/gmat_v4_aldrin.py` (Aldrin) or
     `scripts/gmat_v4_s1l1.py` (S1L1) → a `.script`.
  3. **Run GMAT** batch on the script (feasibility pass first, then optional
     optimality; homotopy Sun+planets → +moon if the corrector stalls — Beeson
     §3.4 / Jones §2.5). Manual human step.
  4. **Parse** the report: `uv run python scripts/parse_gmat_report.py <report>` →
     converged?, maintenance ΔV.
  5. **Apply the V4 predicate** (Phase 3): converged AND within band (Aldrin
     2.768–3.060 km/s; S1L1 58.9–65.1 m/s patched-conic; convergence-only for the
     Mars-perturbed continuous arm).
  6. **Record** the one-line V4 result. State the golden caveat (2.9138 / 62 m/s are
     OUR values; GMAT is the external check) and the Jones continuity tolerance
     (1e-3 km / 1e-6 km/s) as the sourced convergence bar.
  7. **Writeback gate (human, separate step):** on V4-PASS, add `(id, "V4")` to
     `_LEVEL_EVIDENCE` (`validate.py:480`) with an evidence string citing GMAT
     version + script + report + parsed ΔV + pass result (shape: the existing V3
     S1L1 entry at `:615`); this replaces that entry's "canonical V4 (GMAT) remains
     OPEN" line. **NOT done by this plan — flagged as the downstream human step.**
  8. State explicitly: **this is NOT in CI; it is a once-per-confirmed-row stamp**
     (Q7); generator/parser/B-plane kernel are unit-tested without GMAT so the logic
     is covered when GMAT is absent.
- [ ] **Task 4.2** — **Commit (doc):**
  `docs: GMAT V4 B-plane run-book — manual Aldrin + S1L1 protocol, writeback gate (GMAT V4 #171)`.

## Phase 5 — Full gate

- [ ] `uv run pytest tests/scripts/ -q` green (string/geometry/fixture only, no GMAT).
- [ ] `uv run pytest -m "not slow"` green (no regression).
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests` clean.
- [ ] Confirm **no CI test invokes GMAT** (Q7). Commit if cleanup needed:
  `scripts: GMAT V4 lane lint/type gate; confirm GMAT stays out of CI (GMAT V4 #171)`.

---

## Go/No-Go gates (decisive)

1. **GO Phase 0→1 iff** the B-plane kernel round-trips (Task 0.2: recovered B-vector
   reproduces the intended turn to < 1e-9). NO-GO → the targeting math is wrong;
   stop, do not generate scripts with bad goals.
2. **GO Phase 1→2 iff** the Aldrin generator emits a force model + a flyby
   `Target`/`Achieve` block + an initial guess (Phase 1 tests). NO-GO → fix the
   single-flyby template before unrolling the chain.
3. **GO Phase 2→3 iff** the S1L1 generator emits one B-plane Target per App-C Mars
   flyby seeded from `APPC_LEGS` (Tasks 2.1–2.2). NO-GO → the generalisation
   contract is broken.
4. **GO Phase 3→4 iff** the parser's V4 predicate is two-part (convergence AND band)
   with a convergence-only fallback (Tasks 3.3–3.4). NO-GO → the gate would
   mis-classify a non-converged run as V4.
5. **TERMINAL GATE (manual, OUT of this plan):** a human installs GMAT, runs the
   generated script, parses the report, and — only on V4-PASS — writes the
   `_LEVEL_EVIDENCE` row. This plan ships everything *up to* that gate.

## The sourced acceptance criterion (one line)

**V4-PASS iff** GMAT (independent codebase + ephemeris, spec §14 V4) reports **NLP
convergence** of every flyby B-plane Target to Jones' sourced continuity tolerance
(**1.0e-3 km position / 1.0e-6 km/s velocity**, AAS 17-577 §2.5) **AND** the summed
per-flyby maintenance ΔV reproduces our reference within the **self-declared ±5%
band** (Aldrin 2.768–3.060 km/s about 2.9138; S1L1 58.9–65.1 m/s about the
patched-conic 62 m/s; **convergence-only** for the Mars-perturbed continuous arm,
whose TCM GMAT produces for the first time).

## Out of scope (explicit)

- **Installing or running GMAT** — manual human step (Q7); GMAT is not in this
  environment.
- **Any catalogue / `_LEVEL_EVIDENCE` writeback** — downstream human step after a
  real GMAT run (Phase 4 Task 4.1 §7 documents the shape only).
- **Any CI test that invokes GMAT** — OUT (Q7).
- **Re-deriving the patched-conic 62 m/s or the Aldrin 2.9138 km/s** — those exist
  (#169, #134); this plan hands them to GMAT as references under check, it does not
  recompute them.
- **The nbody substrate / shooter / SILVER rungs** — owned by
  `2026-06-06-nbody-harness.md` Phases A–C; this plan is the GMAT B-plane extension
  of its Phase D only.
- **Any `src/` production edit** — the lane is `scripts/` + `tests/scripts/` + docs;
  the only `src/` interaction is read-only reuse (`core/flyby.py`,
  `s1l1_corrected.APPC_LEGS`) and the documented (not executed) `validate.py`
  writeback gate.
