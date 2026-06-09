# Self-seeding longitude-rendezvous construction (task #173)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:executing-plans` or
> `superpowers:subagent-driven-development`. Checkbox steps; strict TDD (write
> failing test → run **red** → minimal impl → run **green** → commit). Work on
> `main` — **do NOT branch** (project rule). uv-managed venv (no pip). Lint/type
> gate before **every** commit:
> `uv run ruff check .` · `uv run ruff format --check .` · `uv run mypy src tests`.
> Fast suite: `uv run pytest -m "not slow"`; n-body / scan probes are `slow`.
>
> This plan is the task-level expansion of the **approved** design
> `docs/notes/2026-06-08-self-seeding-construction-design.md` — read it first,
> especially §1 (the epoch/longitude search), §1.3 (the full residual), §2 (the
> on-family guard + honest risk), §5 (the S1L1 validation gate). Where this plan
> and the design diverge, the design wins, EXCEPT where a design claim no longer
> holds against live code — then the implementer flags it and overrides with the
> verified state.
>
> **GOLDEN/HONESTY (binding).** The EXPECTED side of the S1L1 gate is the row's
> OWN App-C-confirmed geometry (Russell App-C #83 / #166 / #167): epoch 2026-12-15,
> Mars longitude 201.0° on 2027-06-13, Mars v∞ 5.248 (real-eph breathing 3.2–8.0),
> 7 encounters in the 3-SOI band (≈0.0116 AU). The self-seed's found epoch / v∞ /
> miss are EVIDENCE, never imposed. The independent n-body band is the SAME
> #165/#167 band, **never loosened**. **No catalogue writeback anywhere in this
> plan.** A clean FAIL or OFF-FAMILY is a SUCCESS (design §5, §7).

---

## Goal

Replace the App-C printed seed with a SEARCH: given a row's descriptor + the
corrected g/G topology, FIND the Earth-departure epoch and v∞-vector where real DE440
Mars sits at the encounter longitude (the longitude-rendezvous constraint App-C
supplied), then drive the existing #167 corrected-topology tail + independent n-body
confirm. **Prove the method on S1L1 FIRST (known-answer gate, design §5) — that is the
cheapest-decisive go/no-go. Only on PASS, apply to ONE unsourced member.**

Outcome is decisive either way:
- S1L1 self-seed PASS → method validated → one unsourced attempt;
- S1L1 self-seed FAIL → method killed cheaply, failure mode reported → STOP;
- (post-gate) unsourced row CONFIRMED / PARTIAL / OFF-FAMILY — each first-class.

All new code is **additive** in a NEW module
`src/cyclerfinder/search/self_seed.py` (+ its tests). It **imports read-only** from
`search/s1l1_corrected.py`, `search/free_return_chain.py`,
`search/continuation_chain.py`, `core/lambert.py`, `core/ephemeris.py`,
`core/kepler.py`, `core/constants.py`. **No edit to any existing src/ or tests/ file**
(except adding new test files). **No `core/` edit. No catalogue writeback.** The
concurrent agent (#176) owns `scripts/gmat`, `verify/bplane`, `nbody/` — **do NOT
touch those.**

---

## VERIFY-FIRST (live state to re-read before any edit)

The design cites these against live code on 2026-06-08; re-read before editing:

1. `search/s1l1_corrected.py` — `APPC_EPOCH_DAYS` (9325.742435), `APPC_LEGS`,
   `APPC_MARS_TRANSIT` (the answer key: leg 3 → (179.8, 5.248)), `build_seeded_arcs`,
   `reconstruct_mars_encounters`, `g_arc_clearances`, `continuous_chain`. These are
   the #167 tail the self-seed feeds, AND the App-C answer key the S1L1 gate checks
   against. Confirm signatures unchanged.
2. `search/free_return_chain.py` + `search/continuation_chain.py` — the descriptor →
   two-arc (a,e,n_rev) coplanar→DE440 geometry (Stage A). Re-read the entry point and
   what it returns (the closed elements `arc1 a/e/n_rev`, `arc2 a/e/n_rev`,
   encounter epoch `t0`). Confirm the continuation's `t0` is the lever the longitude
   solve will move.
3. `core/lambert.py` — `lambert(...)` at line 511 and `LambertSolution` at line 92.
   Confirm the call signature (r1, r2, tof, mu, multi-rev / branch args) for the
   Stage-B refinement onto true DE440 Mars position. If multi-rev is supported, note
   how `n_rev`/branch is passed.
4. `core/ephemeris.py` — `Ephemeris("de440")` / `.state(body, t_sec)` returning
   (r_km, v_km/s) ecliptic; how `s1l1_corrected` constructs its ephemeris in the
   tests. Confirm `state("M", t)` and `state("E", t)`.
5. The S1L1 synodic period for the scan window: Earth–Mars synodic ≈ 779.9 d
   (~2.135 yr). Confirm any synodic constant in `core/constants.py`; if absent,
   compute from `PLANETS["E"]`/`PLANETS["M"]` periods (do NOT hardcode a literal
   without sourcing it in a comment).
6. `tests/nbody/test_s1l1_corrected_nbody.py` — the independent REBOUND/IAS15
   confirm pattern + the 3-SOI band constant. The S1L1-self-seed n-body confirm
   REUSES this pattern (new test, same band). **Do not edit the nbody dir or its
   files — read the pattern and replicate in a NEW test under `tests/search/` that
   imports the confirm helper, or gate the confirm via the existing `continuous_chain`
   Sun-only path which lives in `s1l1_corrected.py`.** (#176 owns `nbody/`.)

---

## Phase 1 — The longitude residual + synodic scan (NEW, mechanics RED first)

> Build the one genuinely-new primitive (design §1.2): `residual_lon(t_depart)` and
> the synodic-period bracket-and-refine scan. Mechanics gates first — each "expected"
> value is defined BY CONSTRUCTION (a known geometry), never by the scan's own output.
> Tests in NEW `tests/search/test_self_seed.py`.

- [ ] **Task 1.1** — `residual_lon` primitive. In NEW `search/self_seed.py`, add a
  function that, given a G-arc geometry (a, e, n_rev, the descriptor ToF_G) and a
  candidate `t_depart`, propagates the arc (Kepler-Sun from real DE440 Earth at
  `t_depart`, using the #167 recipe `v_sc = v_planet + v∞`) to `t_depart + ToF_G`,
  and returns `lon_sc_encounter − lon_Mars_DE440(t_depart + ToF_G)` (degrees, wrapped
  to (−180,180]). Pure; reuses `_lon_deg`-style helper (define locally — do NOT import
  a private from s1l1_corrected; re-implement the 3-line `arctan2` longitude).
- [ ] **Task 1.2 (RED→GREEN, mechanics)** —
  `test_residual_lon_zero_at_appc_seed_epoch`. CONSTRUCT the G-leg-2 geometry from the
  **App-C answer key** (this is the known-answer fixture, sourced): the App-C
  departure state at 2026-12-15. Assert `residual_lon` at the App-C `t_depart` is
  `< 0.5°` (the App-C seed IS a true longitude rendezvous — #166 measured 4e-7°; allow
  margin for the Kepler-vs-published-arc reconstruction). This proves the residual is
  correct against a sourced known answer. Label: mechanics.
- [ ] **Task 1.3 (RED→GREEN, mechanics)** —
  `test_residual_lon_large_off_phase`. Assert `residual_lon` evaluated a half-synodic
  period away from the App-C epoch is large (> 60°) — the residual actually
  discriminates phase (guards against a degenerate always-zero bug). Label: mechanics.
- [ ] **Task 1.4** — `synodic_longitude_scan`. Add a function that scans `t_depart`
  across one synodic period (window centred on a supplied guess, e.g. the
  continuation `t0`), coarse step (default 10 d), brackets every sign change of
  `residual_lon`, and Brent/bisection-refines each bracket to a root. Returns the list
  of candidate epochs (ALL of them — design §2.2: enumerate, don't optimise).
- [ ] **Task 1.5 (RED→GREEN, mechanics)** —
  `test_scan_surfaces_appc_epoch_for_s1l1`. With the S1L1 G-leg-2 geometry, assert the
  scan returns at least one root within ±5 d of the App-C epoch (2026-12-15). This is
  the core of the validation gate, isolated as a mechanics test. Label: mechanics
  (may be `slow` if the scan is heavy — mark accordingly).
- [ ] **Task 1.6** — Lint/type gate clean. **Commit:**
  `search: longitude-rendezvous residual + synodic scan (self-seed #173)`.

## Phase 2 — Assemble the on-family gate + the full self-seed driver (NEW)

> Wire the scan to the descriptor→geometry (Stage A) and the full multi-term
> on-family residual (design §1.3), producing a `SelfSeedResult` per candidate epoch.

- [ ] **Task 2.1** — `self_seed_g_leg`. Given a row descriptor (the two-arc closed
  elements from `continuation_chain` — Stage A, reused) + a DE440 `Ephemeris`, run the
  synodic scan, and for each candidate epoch optionally refine with a `core/lambert`
  solve onto **true DE440 Mars position** at the arrival epoch (Stage B option 2,
  design §1.2). Return a `SelfSeedResult` dataclass holding, per candidate:
  `t_depart`, `vinf_vec`, `residual_lon`, emerged `vinf_E`/`vinf_M`, `tof_G`,
  `tof_g`, `g_closest_mars_au`, `flyby_bend_deg`/`max_bend_deg`. EVIDENCE fields; the
  descriptor anchors are EXPECTED, supplied separately for the gate.
- [ ] **Task 2.2** — `on_family(result, anchors, *, band)` predicate. Returns True iff
  ALL of design §1.3 hold at the SAME epoch: `|residual_lon| < ~0.5°`, `vinf_E`/`vinf_M`
  within a band admitting real-eph breathing (default ±1.0 of the anchor — wide on
  purpose; design §2.3 inverse-#164 trap), both ToFs in band, g-arc sub-Mars
  (`g_closest_mars_au` ≫ encounter band; aphelion < ~1.4 AU), flyby bend ≤ max. Return
  a structured verdict (which terms pass/fail) — NOT a bare bool — so a partial close
  is diagnosable, not mistaken for a result.
- [ ] **Task 2.3 (RED→GREEN, mechanics)** —
  `test_on_family_true_for_appc_s1l1_seed`. Build a `SelfSeedResult` from the App-C
  S1L1 seed (the answer key) and assert `on_family(...)` returns all-pass against the
  real-eph anchors. CONSTRUCTED from sourced data. Label: mechanics.
- [ ] **Task 2.4 (RED→GREEN, mechanics)** —
  `test_on_family_false_for_off_phase_seed`. A seed a half-synodic away fails the
  `residual_lon` term (and likely v∞). Asserts the gate rejects the off-family basin —
  the #165/2026-06-04 failure mode, now caught by the explicit longitude term. Label:
  mechanics.
- [ ] **Task 2.5** — Lint/type/regression gate clean; `uv run pytest
  tests/search/test_self_seed.py -m "slow or not slow"` green. **Commit:**
  `search: self-seed driver + on-family gate (#173)`.

## Phase 3 — THE VALIDATION GATE: S1L1 self-seed from scratch (slow, decisive)

> The known-answer go/no-go (design §5). Run the FULL self-seed on S1L1 WITHOUT ever
> reading the App-C block, then check it recovers the App-C-confirmed geometry and
> drives the #167 tail to the same 7-in-band confirmation.

- [ ] **Task 3.1 (slow)** — `test_s1l1_self_seed_recovers_appc_geometry`. Drive
  `self_seed_g_leg` on S1L1 using ONLY the descriptor-derived two-arc geometry (Stage
  A via `continuation_chain`) — **do NOT pass any APPC_* constant into the search**
  (the search must be App-C-blind; APPC_* is used ONLY on the assert/EXPECTED side).
  Assert the found seed recovers, against the App-C answer key:
  - departure epoch within tolerance of 2026-12-15 (target ≤ a few days; **report the
    exact delta**);
  - Mars-flyby longitude rendezvous Δlon < ~0.5° (cf. App-C 201.0°);
  - emerged Mars v∞ in the real-eph band (cf. 5.248, breathing 3.2–8.0).
- [ ] **Task 3.2 (slow)** — `test_s1l1_self_seed_drives_167_tail_in_band`. Feed the
  FOUND seed (not App-C) into the #167 tail — reconstruct the G-leg encounter and the
  g-arc clearance via the reused `s1l1_corrected` primitives (or replicate the recipe
  in `self_seed.py` to keep the dependency one-way) — and the independent Sun-only
  confirm (`continuous_chain` / the REBOUND pattern, SAME 3-SOI band, NOT loosened).
  Assert the encounter is in-band at the real-eph v∞. (If reaching all 7 cycles
  requires the App-C per-leg nodes, scope this gate to the first G encounter — the one
  the self-seed actually found — and state that scope explicitly in the results note.)
- [ ] **Task 3.3** — Apply the **binary gate** (design §5) and write the verdict to
  `docs/notes/2026-06-08-self-seeding-s1l1-gate-results.md`:
  - **PASS:** scan surfaced the App-C epoch (≤ tol), on-family residual small there,
    found seed drives the tail in-band → method validated; proceed to Phase 4. Record
    the exact epoch delta, Δlon, emerged v∞, miss.
  - **FAIL:** scan missed the App-C epoch / landed a different basin / longitude
    residual never zeroes at an anchor-consistent epoch → method NOT validated;
    record WHICH failure mode (the App-C key diagnoses it), STOP. **Do NOT run
    Phase 4.** A clean FAIL is a successful task outcome (design §5, §7).
- [ ] **Task 3.4** — Lint/type gate clean. **Commit (code+test):**
  `search: S1L1 self-seed validation gate (#173)`. **Commit (results note):**
  `docs: S1L1 self-seed gate verdict (#173)`.

## Phase 4 — Prove-on-ONE unsourced member (slow; ONLY if Phase 3 PASSES)

> Gated on Phase 3 PASS. Apply to exactly ONE unsourced row — never a batch
> (design §6). Honest three-way verdict, each first-class.

- [ ] **Task 4.1 (slow)** — Pick ONE unsourced member and record WHY: a russell-ch4
  row from the 7 NOT-REACHABLE (#170) — prefer a near-ballistic coplanar one (best V3
  odds), e.g. `russell-ch4-6.44Gg3` (the continuation companion, already has
  `free_return_arcs[]`); OR one ocampo member with a clean descriptor. Confirm its
  descriptor + anchors from `data/catalogue.yaml` (read-only). State the row in the
  test name.
- [ ] **Task 4.2 (slow)** — `test_self_seed_one_unsourced_member`. Run the full
  self-seed + on-family gate + #167-tail reconstruction + independent Sun-only confirm
  (SAME 3-SOI band) + the continuous-TCM budget (`continuous_chain` pattern, V3 budget
  120 m/s — cf. #169/#170). Assert only that the search RUNS and returns a finite,
  fully-audited result; the scientific verdict is REPORTED, not gated.
- [ ] **Task 4.3** — Apply the three-way gate (design §6) and write
  `docs/notes/2026-06-08-self-seeding-one-member-results.md`:
  - **CONFIRMED** — on-family seed found, encounter(s) in-band at real-eph v∞, TCM ≤
    120 m/s ⇒ first unsourced row reachable without a published seed; recommend V3,
    **held to main-session review** (no writeback).
  - **PARTIAL** — in-band but TCM over budget (the #170 powered pattern) ⇒ PARTIAL,
    no writeback.
  - **OFF-FAMILY / EMPTY-SET** — no epoch hosts the rendezvous at the anchor v∞ ⇒
    clean quantified negative; record (feeds the negative-results registry,
    method-versioned).
- [ ] **Task 4.4** — Lint/type gate clean. **Commit (code+test):**
  `search: self-seed one unsourced member probe (#173)`. **Commit (results note):**
  `docs: self-seed one-member verdict (#173)`.

---

## Decisive gate (one-line summary)

Does the App-C-blind synodic longitude scan surface S1L1's App-C departure epoch
(2026-12-15, ≤ a few days) with the full on-family residual small there, and does the
found seed drive the #167 tail to encounters in the SAME 3-SOI band at the real-eph
v∞? **Yes → method validated, attempt ONE unsourced member. No → method killed
cheaply, failure mode reported, STOP.** No path returns "inconclusive": the App-C
answer key makes every S1L1 outcome diagnostic, and a clean FAIL/OFF-FAMILY/EMPTY-SET
is a first-class success.

## Out of scope (explicit)

- Any batch over the 7 russell-ch4 / 194 ocampo members — Phase 4 is ONE member only.
- Any catalogue writeback or V-level change (S1L1 stays its #167-recommended level;
  this task tests a METHOD).
- Any edit to existing src/ or tests/ files; any `core/` edit; anything under
  `scripts/gmat`, `verify/bplane`, `nbody/` (concurrent agent #176).
- Phase 4 entirely, unless Phase 3 fires PASS.
