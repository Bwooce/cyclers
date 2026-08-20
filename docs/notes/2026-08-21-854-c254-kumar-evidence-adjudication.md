# `#854`: does `#827`'s C=2.54 Kumar Table-5 reproduction warrant its own catalogue annotation, and does it narrow `#840`?

**Task:** `#854`, registered 2026-08-21 (found during `#827`'s advisor review, not dispatched).
Adjudicate whether `#827`'s digit-grade reproduction of Kumar-Rawat-Rosengren-Ross (2026) Table 5's
C=2.54 `Wu(3:1) -> Ws(2:1)` intersection state — which touches the catalogued
`vaquero-31-c254-em-cycler-2013` row — warrants its own comment-only `ADDED EVIDENCE` annotation
(parallel to `#828`'s two), and whether it narrows `#840`'s registered round-trip scope to just the
C=2.66 edge.

---

## Verdict (read this first)

**Annotate: yes, comment-only, no level change.** `#827`'s C=2.54 result is a genuinely different
evidence kind from `#828`'s existing annotation on the same row (a digit-grade match to a state a
peer-reviewed paper itself printed, at the paper's own mass ratio — vs `#828`'s self-consistency-only
Newton closure at this project's own mass ratio) and is not redundant with it. Per `#828` Sec. 2's own
ruling this is still not a tier-promotion question; the annotation is comment-only.

**`#840`'s scope: does NOT narrow to just C=2.66 — `#827`'s note overclaimed this.** `#840` is a
DIRECTION question (forward vs reverse) *inside a fixed model* (this row's own registry mass ratio);
`#827`/`#854` is an EVIDENCE-SOURCE question (self-consistency vs published-digit match) computed at
a *different* mass ratio (Kumar's own printed value, `Δmu ≈ 1.24e-10` from this project's registry
value). Those are different axes. `#840`'s C=2.54 leg — a reverse-direction `find_free_transfer` run
at THIS row's own registry mu, verified against THIS row's own recorded state to the tight tolerance
`#828` achieved at matched mu (`~1e-15`) — is still fully open. What changes is its *character*:
existence at C=2.54 is now independently pre-confirmed at a neighbouring mu (on top of the
time-reversal argument `#840` already cited), so the C=2.54 leg of `#840` is now a confirmation run,
not a blind search. The C=2.66 leg is untouched by any of this and remains a full open item, same as
before.

---

## 1. Independent re-verification of `#827`'s C=2.54 result (not taken on say-so)

Per this project's standing discipline, `#827`'s own finding-task framing was not trusted directly.
Two independent checks were run this task, both from scratch, neither reading a value out of
`results.json` and treating it as ground truth:

**(a) Seed-level reconstruction** (mirrors the existing `test_connection_reconverges_and_matches_
print_at_2_86_type1` / `..._3_15_type2` pattern in `tests/search/test_kumar_em_resonant_
heteroclinics.py`, extended here to C=2.54, which had no such test): rebuilt both Table-6 nodes at
C=2.54 via `build_kumar_node` (which itself re-derives against the print and rejects disagreement),
seeded `vcc.refine_connection` with *only* the recorded phase indices (`branch_u=-1, branch_s=-1,
k_u=19, k_s=14, tau_u=3.7202948069488033, tau_s=5.98472243980346` — never the converged crossing
state), and ran `vcc.verify_connection` + `keh._match_against_print` fresh.

Result: `conn.converged=True`, `residual=3.363e-11`; `evidence.passed=True`,
`ydot_signs_match=True`, `full_state_gap=4.651e-07`; `match_distance=7.610531887910121e-07`,
`runner_up_distance=1.307`. All bit-identical to `results.json`'s recorded row — expected, not
notable on its own (deterministic Newton from an identical seed on the same machine reproduces
exactly); the actual independence here is that the seed indices were typed in from the recorded
record, not read back from a live connection object, and the node re-derivation and match are both
freshly computed.

**(b) End-to-end reconstruction from Kumar's printed state alone** (the stronger check, closing the
one gap (a) leaves — that the seed indices `(19, 14, -1, -1)` came from an honest search rather than
being back-fit): ran `keh.reproduce_table5_intersection(system, 2.54)` directly, which seeds *only*
from `kumar_table5_state6(2.54)` (Kumar's own printed Table-5 row, legitimate published input) and
performs its own candidate search — no `results.json` read anywhere in this path.

    2026-08-21T08:10:53+10:00 start; 50.0 s elapsed (results.json recorded 123.1 s on the prior run's
    machine — a speed difference, not a discrepancy)
    matched=True, match_distance=7.610531887910121e-07, runner_up_distance=1.307
    n_candidates=11, n_refined=1, n_converged=1
    connection: k_u=19, k_s=14, branch_u=-1, branch_s=-1, tau_u=3.7202948069488033,
    tau_s=5.98472243980346

Landed on the identical seed indices and match distance as `results.json`, from a completely
independent search starting only from the printed digits. This closes the non-circularity question:
the C=2.54 row is a real, reproducible result, not an artifact of the checkpoint file.

## 2. Node identity: confirmed the same physical orbit, at the expected mu-shift scale

`node31 = build_kumar_node(system, 3, 2.54)` (Kumar's own mu, `1.2150584270572e-2`) vs.
`vaquero-31-c254-em-cycler-2013`'s recorded `state_nd`/`period_nd`/`jacobi_constant` (this project's
own DE440-registry mu, `0.01215058439469525`):

| quantity | Δ |
|---|---|
| `x0` | 3.398e-10 |
| `ydot0` | 1.061e-09 |
| `period_nd` | 3.345e-10 |
| `jacobi_constant` | 5.329e-15 |

These sit at the `Δmu ≈ 1.24e-10` scale the module's own docstring documents (`KUMAR_MU`'s comment:
"differ by ~1.24e-10 absolute, which is above the digit-grade comparison floor") — **not** the
`~1e-15` scale `#828` Sec. 3 got comparing `#822`'s re-converged node to this same row *at matched
mu*. This is the correct signature for "the same orbit, computed at a slightly different mu," and it
is the evidence this row and Kumar's Table-6 3:1 node at C=2.54 are the same physical object, not a
coincidentally nearby one.

## 3. Why this is a different evidence kind from `#828`'s existing annotation (not redundant)

`#828`'s existing block on this row records `#822`'s Wu(2:1)->Ws(3:1) *forward* connection at C=2.54,
computed entirely in this project's own model (registry mu), verified only for internal
self-consistency (Newton residual, ghost-guard margins, an independent-Radau cross-check — all
checks the connection is a real solution of *this project's own* equations of motion, not that it
matches any external source).

`#827`'s C=2.54 result is the Wu(3:1)->Ws(2:1) *reverse* connection, computed at Kumar et al.'s own
printed mu, matched digit-grade to a state a peer-reviewed paper (DOI 10.1016/j.asr.2025.12.005)
itself published in its own Table 5. That is evidence a *second, independent source* — not just a
second run of this project's own code — puts a heteroclinic connection at this exact orbit and this
exact Jacobi constant. Neither the direction, the mu, nor the evidence standard (self-consistency vs.
external-source digit match) overlaps `#828`'s block. It is not redundant, and it is not a tier
question (per `#828` Sec. 2, a cross-object connection is not what spec §14's ladder measures at any
rung) — a comment-only annotation is the correct, established response
(`ADDED EVIDENCE (..., no level change)`, the same convention `#828`/`#834` used).

## 4. `#840`'s scope, examined directly (not inherited from `#827`'s framing)

`#840`'s registration is explicit that its deliverable is a reverse-direction run using `#822`'s own
`find_free_transfer` machinery, "so each of the two in-band catalogued rows … has a demonstrated
ROUND TRIP **at its own C**." Read plainly, "at its own C" means at the row's own model — this
project's own registry mu, the value literally recorded in the row's `mass_ratio` field — because
that is the model `#822`'s forward connection and this row's own V1 evidence are both computed in.
A round-trip transport claim recorded on a row has to hold in that row's own model to be a round-trip
claim about *that row*.

`#827`/`#854`'s reverse-direction connection at C=2.54 does not hold in that model — it holds in a
neighbouring one, Kumar's own printed mu, offset by `Δmu ≈ 1.24e-10`. Sec. 2 shows the resulting node
is the same orbit to `~1e-9..1e-10`, not to the `~1e-15` `#828` achieved at matched mu — the shift
is real and it is exactly the scale the project's own digit-grade discipline (`KUMAR_MU`'s docstring)
already treats as significant, for a different purpose (a published-digit match) here now extended
to a second purpose (a same-model transport claim). `#840`'s literal deliverable — a reverse-direction
`find_free_transfer` run at this row's own registry mu, verified against this row's own recorded
state to `#828`'s own tight tolerance — is therefore still fully open at C=2.54, exactly as it is at
C=2.66.

**What genuinely changes:** `#840`'s own registration already argued existence is guaranteed by
CR3BP time-reversal symmetry given the forward hit. `#827`/`#854` adds an actual independent
computation — at a neighbouring mu, not the row's own — that a reverse connection does exist near
this exact node and Jacobi constant, with a wide margin (match_distance 7.6e-07 vs. `KUMAR_MATCH_TOL`
1e-4, runner-up separation four orders of magnitude larger). That does not satisfy `#840`, but it
does mean `#840`'s C=2.54 leg is now a **confirmation run** — a near-certain, low-risk re-derivation
at the row's own mu — rather than an open search with only a symmetry argument behind it. The C=2.66
leg is untouched by any of this: `#827`'s Table-5 rows never include C=2.66 (`#827`'s own table:
`{2.54, 2.70, 2.86, 3.00, 3.05, 3.10, 3.15}`), so no analogous pre-confirmation exists there.

**Conclusion on the overclaim question:** `#827`'s note states this "directly overlaps `#840`'s
registered scope … for the C=2.54 edge specifically" and asks `#854` to adjudicate "whether it
narrows `#840`'s remaining scope to just the C=2.66 edge." The overlap claim is fine (it does touch
the same row and the same edge); the narrowing claim, read as "the C=2.54 half of `#840` is now
covered," is **not correct** and is corrected here — `#840` needs its own dispatch at BOTH C=2.54 and
C=2.66, unchanged in scope, with C=2.54's registration text updated to record the pre-confirmation.

## 5. Catalogue edit

Comment-only `ADDED EVIDENCE (#854 …)` block appended to `vaquero-31-c254-em-cycler-2013`'s
`validation_level` comment, immediately following `#828`'s existing block on the same field. The new
block's opening two sentences explicitly reconcile with `#828`'s existing "no round-trip transport
claim is available at THIS row's own C" caveat, so a reader does not see the two blocks as
contradicting each other (`#828`'s caveat is scoped to this row's own model; `#827`'s reverse
connection is at a neighbouring model).

`yaml.safe_load` of the catalogue before and after this edit was compared programmatically and
found **identical** (parsed data unchanged — a pure comment edit).

## 6. Verification run

- Independent seed-level reconstruction and end-to-end `reproduce_table5_intersection(2.54)` run:
  Sec. 1, this task, scratchpad script + direct `uv run python -c` invocation (not committed —
  ephemeral verification, not a new artifact; the module and its existing tests already cover this
  code path, this task only exercised the untested C=2.54 point).
- Node-identity comparison: Sec. 2, this task.
- `python3 -c "import yaml; assert yaml.safe_load(open('<pre>')) == yaml.safe_load(open('<post>'))"`
  equivalent comment-only check: pass.
- Full `uv run pytest tests/data tests/search -q` ratchet (never a subset, per
  `[[feedback_catalogue_edits_run_all_ratchets]]`): see commit log for exit status.
- `uv run ruff check .`, `uv run ruff format --check .`, full `uv run mypy src tests`: see commit
  log for exit status.
