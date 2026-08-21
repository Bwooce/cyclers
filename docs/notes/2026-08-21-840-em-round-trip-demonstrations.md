# `#840`: Earth-Moon 2:1<->3:1 free-transfer round trip at the two catalogued band edges

**Task:** `#840`, registered 2026-08-12 (found during `#828`), dispatched 2026-08-21 (corrected
2026-08-21 by `#854`, scope unchanged). `#822` demonstrated the REVERSE-direction heteroclinic
connection Wu(3:1)->Ws(2:1) at C=2.60 only -- a Jacobi constant where neither family node is a
catalogued row. `#828` recorded the FORWARD direction Wu(2:1)->Ws(3:1) touching two catalogued
rows at their own band edges: `vaquero-31-c254-em-cycler-2013` (C=2.54) and
`vaquero-21-c266-em-cycler-2013` (C=2.66), both self-consistency-only evidence at this project's
own registry mu (`0.01215058439469525`). This task runs the REVERSE direction at those SAME two
Jacobi constants, at the same registry mu, so each row ends up with a demonstrated round trip
(both directions, same C, same mu) in this project's own model.

**Verdict: BOTH C values converged and verified.** C=2.54's reverse connection lands almost
exactly on the CR3BP time-reversal-symmetry prediction of the forward entry's own crossing
(~6e-9) -- genuinely the mirror-image manifold intersection. C=2.66's reverse connection required
a diagnosed epsilon adjustment and converges on a DIFFERENT, but equally genuine and fully
verified, manifold intersection (this family pair's intersection structure is rich, `#822`'s own
note already documents this). Both rows now carry a demonstrated round trip in this project's own
model. Per `#828` Sec. 2 (still binding, re-affirmed by `#854`), this is evidence completeness,
**not** a validation-tier promotion or `orbit_class` change -- `data/catalogue.yaml` is untouched.

Module: `src/cyclerfinder/search/vaquero_em_cycler_connections.py` (new: `find_free_transfer_reverse`,
`find_free_transfer_reverse_near_prediction`, `predict_reverse_crossing`; refactored
`_find_free_transfer_between` shared core). Driver: `scripts/screen_840_em_reverse_round_trip.py`.
Full per-C record: `data/found/840_em_reverse_round_trip/results.json`. Registry writeback:
`data/manifold_connections.yaml` (two new entries + `reverse_of`/`round_trip_note` updates on the
two existing forward entries + the Kumar entry's `round_trip_note`). Tests:
`tests/search/test_vaquero_em_cycler_connections.py` (4 new tests), `tests/data/test_manifold_connections.py`
(3 new/updated tests).

---

## 1. Method

`#822` Sec. 4's own C=2.60 reverse demo is the direct template: the SAME
`manifold_section_crossings` / `find_connection_seeds` / `refine_connection` / `verify_connection`
pipeline, with `node31` (3:1) as the unstable origin and `node21` (2:1) as the stable destination
(roles swapped from the forward direction). This is now first-class module code
(`find_free_transfer_reverse`), factored via a shared `_find_free_transfer_between` core so the
reverse direction reuses the forward direction's own logic verbatim rather than duplicating it.

**A new technique, promoted from this task's own diagnosis**: both orbits are x-axis symmetric (built
by `correct_symmetric_fixed_jacobi`), so the CR3BP time-reversal symmetry
`(x, y, xdot, ydot, t) -> (x, -y, -xdot, ydot, -t)` predicts where the REVERSE connection's own
`{y=0}` crossing sits directly from the ALREADY-VERIFIED forward connection's own crossing
(`predict_reverse_crossing`), without any search. `find_free_transfer_reverse_near_prediction`
ranks candidate seeds by proximity to this prediction (rather than native `(x, xdot)`
close-approach distance) and sweeps `epsilon` upward when the closest candidate converges but
fails `verify_connection`'s evidence-quality gates -- the same `#822`-documented tradeoff
("Manifold-offset ε as an evidence-quality control"), never loosening the gates themselves.

## 2. What was tried first, and the diagnosed pathology at C=2.66

A first blind `find_free_transfer_reverse` scan at the sibling-module default `epsilon=1e-4`
converged 6/11 refined candidates at C=2.66, but **every one failed `verify_connection`'s
forward-reapproach gate** (`forward_distance > 0.5`), including the seed closest to the
symmetry-predicted crossing itself (`forward_distance=1.58`). Every OTHER gate passed with good
margin on that same candidate (full-state gap 6.8e-6, Radau gap 3.7e-6, backward re-approach
1.3e-6, ghost margins 0.58/0.60) -- confirming this is a genuine, zero-Δv heteroclinic geometry
hitting an evidence-quality floor at this epsilon, not a spurious Newton false-positive
(`tests/search/test_vaquero_em_cycler_connections.py::test_reverse_c266_at_default_epsilon_fails_forward_gate_not_other_gates`
pins this diagnostic finding as a permanent regression check).

Raising epsilon shortens each leg (fewer Floquet-amplification periods, the same mechanism `#822`
already used at several forward-direction C points). C=2.54 needed no adjustment at all --
epsilon=1e-4 converges immediately on the very first/closest-to-predicted candidate. C=2.66
needed epsilon=2e-4, found via `find_free_transfer_reverse_near_prediction`'s own systematic
sweep (trying up to 5 predicted-ranked candidates per epsilon).

## 3. Results

### C=2.54 -- genuinely the mirror image of the forward connection

`Wu(3:1) -> Ws(2:1)`, `branch_u=-1, branch_s=+1`, `k_u=26, k_s=38`, `ε=1e-4`:

```
node 3:1 (row vaquero-31-c254-em-cycler-2013): x0=0.9013301668020125  ydot0=-0.8462249954358775  T=6.269604424886022   C=2.54
node 2:1 (uncatalogued):                       x0=1.0905363960533268  ydot0=-0.8231863180949408  T=5.941227735609639   C=2.54
converged: tau_u=4.258037700729926  tau_s=0.316633984139843   Newton residual 5.003e-10
crossing (y=0): x=+0.233225082  xdot=+1.785596427  ydot=+1.552406 (ydot>0 BOTH legs)
transit: t_u=+25.90 nd  t_s=-37.06 nd
```

Evidence: full 4-state gap `3.37e-06`; Radau independent gap `1.35e-06`; ghost margins
`1.800 / 0.243` (1800x/243x the guard); backward re-approach `4.40e-08`; forward `1.77e-03`;
per-leg drift `2.07e-11 / 2.15e-12`; seed offsets `7.83e-09 / 1.25e-08`. **This crossing matches
the CR3BP time-reversal-symmetry prediction of the FORWARD entry's own crossing
(`(0.23322507758762506, -1.785596430399187)` -> predicted reverse `(0.23322507758762506,
+1.785596430399187)`) to `~5.7e-9`** -- genuinely the mirror-image manifold intersection of the
same forward connection, not a coincidentally nearby different one.

### C=2.66 -- a genuine, different intersection point

`Wu(3:1) -> Ws(2:1)`, `branch_u=-1, branch_s=+1`, `k_u=29, k_s=33`, `ε=2e-4`:

```
node 3:1 (uncatalogued):                       x0=0.8919375112041409  ydot0=-0.7577709445440712  T=6.283952207405823   C=2.66
node 2:1 (row vaquero-21-c266-em-cycler-2013): x0=1.0338302047346954  ydot0=-0.9089334377051435  T=5.662843584779122   C=2.66
converged: tau_u=4.122557122466014  tau_s=2.595563242258196   Newton residual 1.379e-10
crossing (y=0): x=+0.959028989  xdot=+0.327193834  ydot=-1.014989 (ydot<0 BOTH legs)
transit: t_u=+33.61 nd  t_s=-38.29 nd
```

Evidence: full 4-state gap `2.86e-06`; Radau independent gap `1.77e-06`; ghost margins
`0.334 / 0.336` (334x/336x the guard); backward re-approach `2.08e-09`; forward `5.19e-02`
(well under the `0.5` ceiling); per-leg drift `1.33e-11 / 2.21e-12`; seed offsets
`3.77e-08 / 1.78e-08`. **This crossing is `~0.26` away from the symmetry prediction of the
forward entry's own crossing** (`(0.8691637553536629, -0.5754661265104016)` -> predicted reverse
`(0.8691637553536629, +0.5754661265104016)`) -- a genuinely different intersection point, not the
mirror image, confirming `#822`'s own note that this family pair's manifold intersection
structure is "rich", not a single point. Still a fully independent, verified connection at the
same C and mu, which is all `#840`'s round-trip claim requires.

## 4. Node identity, independently re-checked this task

Per `[[feedback_orbit_closure_discipline]]`, node identity was independently re-verified this
task (fresh `build_vaquero_overlap_node` calls, not reused from `#828`'s own numbers):

| row | Δ`state_nd` (x0, ẏ0) | Δ`period_nd` | Δ`jacobi_constant` | ΔBarden `nu` (rel) |
|---|---|---|---|---|
| `vaquero-31-c254-em-cycler-2013` | 0.0, 3.109e-15 | 1.066e-13 | 5.329e-15 | 1.400e-08 |
| `vaquero-21-c266-em-cycler-2013` | 0.0, 4.441e-15 | 5.151e-14 | 7.994e-15 | 4.603e-13 |

These match `#828`'s own reported deltas exactly (3.1e-15/1.1e-13/5.3e-15/1.4e-08 and
4.4e-15/5.2e-14/8.0e-15/4.6e-13 respectively), confirming the SAME catalogued orbits, re-derived
from a code path built by a later task.

## 5. Registry writeback (`data/manifold_connections.yaml`)

Two new entries added: `em-vaquero-hetero-wu31c254-ws21c254-2026` and
`em-vaquero-hetero-wu31c266-ws21c266-2026`, each with `reverse_of` pointing at its forward
counterpart. The two EXISTING forward entries' `reverse_of` fields were set to point back
(round trip is symmetric -- both sides cross-reference), and their `round_trip_note` fields
updated from "NOT demonstrated at this C" to the new demonstrated status. The Kumar entry
(`em-kumar-hetero-wu31c254-ws21c254-2026`, a DIFFERENT model -- Kumar's own mu, not this
project's registry mu) keeps `reverse_of: null` (no same-model forward counterpart exists in
that model) but its `round_trip_note` was corrected to no longer describe `#840` as an open
follow-up.

Schema validation: `uv run check-jsonschema --schemafile data/manifold_connection.schema.json
data/manifold_connections.yaml` -- clean. `uv run python -m cyclerfinder.data.validate_connections
check` (schema + semantic + referential layers) -- clean, zero violations.

**Sanitizer check**: every free-text field written or edited this task (both new entries'
`identity_evidence`/`derivation`/`evidence_class`/`round_trip_note`, and the three edited
`round_trip_note` fields on pre-existing entries) was extracted and run through
`cyclers.space/src/lib/catalogue.ts`'s `sanitizeCatalogueText` standalone (the same technique a
prior session used to catch 4 real grammar breaks in this registry). No `#NNN`-as-grammatical-
subject pattern was used anywhere in this task's own new text (task references are either dropped
entirely or referenced via entry ids, e.g. "see em-vaquero-hetero-wu31c254-ws21c254-2026", which
contain no `#` and are untouched by the sanitizer's regexes) -- every field this task wrote passes
through byte-identical (module trailing-newline trimming aside).

## 6. Verification run

- `tests/search/test_vaquero_em_cycler_connections.py`: 16 tests (12 existing + 4 new: symmetry
  prediction plumbing, C=2.54 reverse convergence + symmetry match, C=2.66 reverse convergence at
  the diagnosed epsilon + confirmed-different crossing, and the C=2.66 default-epsilon
  forward-gate-failure regression pin) -- all pass, none `@pytest.mark.slow`.
- `tests/data/test_manifold_connections.py`: 20 tests (updated entry-count assertion to 5 ids,
  plus 2 new: `#840`-artifact transcription guard, forward/reverse `reverse_of` symmetric
  cross-reference check) -- all pass.
- `uv run ruff check .` / `uv run ruff format --check .`: clean.
- Full `uv run mypy src tests scripts/screen_840_em_reverse_round_trip.py`: clean (858 source
  files, zero issues).
- Full `data/found/822_vaquero_em_free_transfer/results.json` was read-only this task (never
  modified) -- the forward connections it records are untouched.

## 7. OUTSTANDING.md

`#840` bullet added: ✓ DONE, both C values converged and verified, honest per-point results as
above (C=2.54 the literal mirror image, C=2.66 a genuinely different intersection), round-trip
status of both rows now demonstrated, no catalogue.yaml or validation-tier change (per `#828`
Sec. 2, re-affirmed here).
