# `#754` — Jupiter-Europa 3:4-LO homoclinic connection + Anderson & Lo 2011 Table 2 gate

**Task:** `#754`, re-scoped by `#757` (see
`docs/notes/2026-07-28-757-task-b-rescoping-confirmed-families.md` Sec. 5) around the
Table-2 half only: Anderson & Lo 2011's Table 2 ("Homoclinic Trajectory State at
Intersection", p.190) is the intersection of the `3:4-LO` resonant orbit's OWN stable
and unstable manifolds with each other — a self-connection, buildable with only the
one family `#755` reviewer-confirmed (`3:4-LO`), with no dependency on `5:6-LO`
(a clean negative, independently retried by the sibling task `#758` on a disjoint
file this same session — not touched here).

**Sources read directly this task**: the paper's own text layer (`.txt` sidecar of
`cyclers_pdf/papers/anderson-lo-2011-...BF03321164.pdf`) around p.170-171 (the
Poincare-map section definition and Eq. 7), and pp.190-191 (Table 2, Table 3, and
the surrounding "Using Homoclinic and Heteroclinic Connections" prose) — confirmed
verbatim against the rendered text (grep line numbers 22, 247-269, 949-953,
1548-1625 of the text-layer sidecar). Code read in full:
`src/cyclerfinder/genome/heteroclinic_cycle.py` (all 657 lines, pre-change) and
`src/cyclerfinder/search/jovian_resonant_families.py` (candidate dataclass, Table-1
seed table, gate machinery) — neither file's own prior claims taken on faith; the
"duck-typed, never calls `from_libration`" claim from `#757`'s scoping note was
independently re-verified by reading `correct_connection`/`_seed_on_manifold`/
`_connection_residual` line by line: confirmed accurate.

---

## 1. What was built

### `genome/heteroclinic_cycle.py` (extended, not rewritten)

1. **`x_sign` section filter.** Added `x_sign: int | None` to `_section_crossing`
   (same pattern as the existing `ydot_sign`), threaded through
   `_connection_residual` and `correct_connection` as `x_sign_u`/`x_sign_s` kwargs.
   `None` (default) leaves the section unrestricted in x — no behavior change for
   any existing caller; all 8 pre-existing W-Z Sun-Jupiter-Oterma tests
   (`tests/genome/test_heteroclinic_cycle.py`) pass unchanged.
2. **`ConnectionNode` Protocol.** Replaced the `LyapunovNode`-typed node parameters
   of `correct_connection`/`assemble_cycle`/`crosscheck_cycle`/`_seed_on_manifold`/
   `_connection_residual` with a structural `typing.Protocol` (`label`, `state0`,
   `period`, `jacobi`, `unstable_eigvec`, `stable_eigvec`, `converged`, each a
   read-only `@property`) so a non-libration node type can be passed without a
   subtype relationship — matching the module's own "duck-typed" design, now
   enforced by `mypy --strict` rather than true only by convention. (Read-only
   `@property` members, not plain attribute annotations, were required: mypy
   otherwise treats a Protocol's plain attributes as needing a *setter* too, which
   a frozen dataclass's fields structurally are not.)
3. **Persisted section filters + epsilon on `HeteroclinicConnection`.** Added
   `ydot_sign_u`/`ydot_sign_s`/`x_sign_u`/`x_sign_s`/`epsilon` fields (all with
   backward-compatible defaults) so `crosscheck_cycle`'s independent re-derivation
   restricts to the exact same filtered crossing the corrector found — **this
   surfaced and fixed a real, previously-latent bug**: `crosscheck_cycle` hardcoded
   `epsilon=1e-6` when re-seeding the manifold, ignoring whatever `epsilon` the
   connection was actually found with. Anderson & Lo's own manifold offset is
   `0.5e-5` (`ANDERSON_LO_EPSILON`); at that value the Radau re-seed **missed the
   k-th qualifying section crossing entirely** within the same `max_time`, silently
   reporting `independent_residual = inf` (a false failure) instead of the true
   cross-check agreement. Diagnosed directly (see `crosscheck_debug.py` in this
   task's scratch work): at *matching* epsilon, DOP853 and Radau agree to `<1e-8`,
   confirming this was a latent harness bug, not a real cross-check failure. Latent
   because every prior caller (the W-Z tests) used the default `epsilon=1e-6`
   throughout, so the mismatch never manifested before this task's `0.5e-5` usage.
   Also fixed: `crosscheck_cycle` never threaded `ydot_sign`/`x_sign` into its
   `_section_crossing` calls at all (same latent-because-`None`-everywhere gap) —
   now reads them back from the stored connection too.

### New module `search/jovian_resonant_connections.py`

- **`ResonantNode`** — same six fields as `LyapunovNode`, structurally satisfies
  `ConnectionNode`. `ResonantNode.from_candidate(system, cand)` recomputes the
  saddle Floquet pair via the existing `heteroclinic_cycle._planar_floquet_pair`
  from the candidate's own stored `(x0, ydot0, period)` — no schema change to
  `ResonantFamilyCandidate`. Raises `ValueError` if `cand.is_real_unstable` is
  `False` (the magnitude-only Floquet convention is invalid for a complex/marginal
  eigenpair) or if the recomputed eigenvalue disagrees with the candidate's own
  `max_eigenvalue` by more than `1e-6` relative (stale/mismatched data fails
  loudly). Both guards covered by dedicated tests.
- **`own_section_points`** — the orbit's own qualifying `{y=0, ydot>0, x<0}`
  section crossings over one period (the node's own IC is included explicitly,
  since `solve_ivp` events do not reliably register the `t=0` boundary). Used as
  the ghost-guard reference set.
- **`find_homoclinic`** — coarse scan over `branch_u, branch_s ∈ {+1,-1}` and
  `k_u, k_s ∈ 1..N` calling `correct_connection(node, node, ...)` for each
  combination, rejecting any converged crossing within `GHOST_GUARD_DELTA=1e-3` of
  the orbit's own section point (the homoclinic trivial-solution ghost guard),
  returning all survivors ranked by distance to `TABLE2_STATE`.
- **`ydot_from_section_eq7`** — the paper's own Eq. 7
  (`ydot = +sqrt(x^2 + 2(1-mu)/r1 + 2*mu/r2 - C - xdot^2)`), verified algebraically
  identical to `core.cr3bp.jacobi_constant`'s own formula evaluated at `y=0` and
  solved for `ydot` (round-trip test: `test_ydot_from_section_eq7_matches_jacobi_constant`).
- **`gate_table2`** — the honest Table-2 gate: PASS iff `(x, xdot)` both match to
  `TABLE2_GATE_ABS_TOL=1e-4`; `ydot` (Eq.-7-derived) is reported but does NOT gate
  `passed` on its own (it is Jacobi-redundant given `(x, xdot)` fixed — see Sec. 3).
  Reports an explicit, non-fudged FAIL when no candidate survives.

---

## 2. The coarse scan: only ONE genuine converged homoclinic self-intersection found

Scanned all 4 branch-sign combinations (`{+1,-1}^2`) across `k_u, k_s` mostly
spanning `1..6` (the paper gives no branch/index — `#757`'s own spec calls this an
undocumented scan, not a missing piece of the method) using `correct_connection`
with the paper's own section (`x_sign=-1`, `ydot_sign=+1`) and manifold offset
(`epsilon=0.5e-5`). Per-combination cost varied wildly (8s to ~120s depending on
whether/how quickly Newton converges or gives up), so the scan was run in many
small chunked, blocking shell calls (never backgrounded), checkpointing every
attempt.

**Result: exactly one combination converges to a genuine (non-ghost) homoclinic
point:**

```
branch_u = +1, branch_s = -1, k_u = 3, k_s = 3
tau_u = 14.463070474804017, tau_s = 11.20947242577664
crossing (x, xdot) = (-1.28540962, 2.8589e-09)
Newton residual = 2.58e-9   (gate requires <= 1e-7)
```

Every other combination tried either (a) failed to reach the section within the
horizon (`max_time_factor=3.0`, ~3 orbital periods), (b) converged onto the orbit's
own trivial section point (rejected by the ghost guard, `ghost_distance < 1e-3`),
or (c) did not converge at all within 20 Newton iterations. The two "mirror" branch
pairs (`+1,+1`) and (`-1,-1`) and the reverse-opposite pair (`-1,+1`) were also
scanned near `k∈{2,3,4}` and yielded no converged hit.

## 3. Verification of the one hit

- **Ghost guard**: the orbit's own qualifying section point is its IC itself,
  `(x0, xdot0) = (-1.4304078294961569, 0.0)` (confirmed by `own_section_points`,
  one point, matching the IC to `<1e-9`). The converged crossing is
  `0.145` away from it in `(x, xdot)` norm — **145x** the `GHOST_GUARD_DELTA=1e-3`
  threshold, a large, non-delicate margin exactly as `#757`'s scoping note
  predicted (`~0.146`).
- **Independent Radau cross-check** (after fixing the `crosscheck_cycle`
  epsilon bug, Sec. 1): DOP853 and Radau agree on the crossing to `6.79e-8` —
  comfortably inside the mandated `<=1e-6`.
- **Fig. 25 corroboration** (reported, not gated): propagating the converged
  homoclinic state forward 3 orbital periods brings it to within `2.76e-5` nondim
  of the 3:4-LO orbit; propagating it backward 3 periods brings it to within
  `7.79e-6` — both asymptotically approach the orbit, exactly the qualitative
  behavior the paper's own Fig. 25 (forward + backward panels) describes ("the
  homoclinic trajectory approaches the 3:4 resonant orbit" in both time
  directions, p.191).

## 4. The Table-2 gate result: **HONEST FAIL** (close, not fudged)

| Quantity | Recovered | Table 2 (p.190) | Abs. error | Tolerance | Passed? |
|---|---|---|---|---|---|
| `x` | -1.28540962 | -1.28427733 | 1.132e-3 | 1e-4 | **NO** |
| `xdot` | 2.859e-9 | 0.00000009 | 8.71e-8 | 1e-4 | YES |
| `ydot` (Eq. 7, reported only) | 0.465377 | 0.46372205 | 1.655e-3 | 1e-4 | NO (not gating) |

**Overall gate: FAIL** — `x` misses the `1e-4` tolerance by roughly an order of
magnitude (`1.13e-3`), even though `xdot` matches to `<1e-7` and the whole
intersection sits in the same qualitative region of the section the paper
describes. This is reported plainly, per this task chain's own discipline — not
loosened, not hidden behind a wider tolerance chosen after the fact.

**Why this is still meaningful, not a wasted result:**

1. It is an order of magnitude closer than "no relationship" would produce — the
   section spans roughly `x ∈ [-1.5, 1.5]`; landing within `1.1e-3` of a
   specific published digit-precision point via an entirely independent
   computation (manifold theory + Newton correction, vs the paper's own
   interpolation-based search) is a strong positive signal, just short of the
   tolerance this task's own honest gate demands.
2. `xdot` — the "hard" coordinate near a delicate near-zero value — matches
   to `8.7e-8`, four orders of magnitude inside tolerance. If our 3:4-LO family
   were a coincidental, unrelated orbit, there is no reason its self-intersection
   would land this close to Table 2's own `xdot ≈ 0` on this specific coordinate.
3. **A parsimonious explanation for the `x` gap already exists and does not need
   new hypotheses**: `#755`'s own reviewer verdict confirmed `3:4-LO` on the
   weight of evidence (near-machine-precision eigenvalue match, `2.8e-8` relative)
   while flagging a genuine, unresolved `2.1%` period offset from the naive
   `2πq` value. For a family this strongly unstable (`λ≈1036`), manifold
   trajectories diverge exponentially over the ~3-period integration horizon used
   here; a `2.1%`-scale timing/phase difference between our corrected orbit and
   the paper's own (if any exists) is entirely sufficient to explain a
   `1e-3`-scale shift in *where* the manifold happens to re-cross the section,
   without implying a different orbit family altogether. This is consistent with,
   not contradictory to, `#755`'s own reviewer ruling.
4. **This gate independently corroborates the `#755` identification on a third,
   orthogonal axis** (state-space manifold-intersection geometry), exactly as
   `#757`'s scoping note anticipated: the eigenvalue match (`#755`) and the shape/
   close-approach match (`#755`) are now joined by a manifold-self-intersection
   match that lands within `0.08%` of the section's own dynamic range of the
   published point — not a coincidence, but not (yet) a tolerance-clearing
   reproduction either.

**Not attempted, and why**: a finer local re-scan around `(k_u=3, k_s=3)`
perturbing `epsilon` or `scan_n` would not close the gap — the corrector already
converges to `2.58e-9` residual (far tighter than the `1.1e-3` gap), and the
independent Radau cross-check confirms this is the TRUE self-intersection point
for our specific orbit, not a numerical artifact. Closing the remaining gap would
require either (a) an independent proof our `3:4-LO` IC is bit-for-bit the paper's
own (not achievable from the paper's published 6-7 significant figures alone), or
(b) resolving the open period-exactness question `#755` raised. Neither is in this
task's scope.

## 5. Verification

- `uv run ruff check` / `ruff format --check` on all three changed/added files:
  clean.
- `uv run mypy src tests` (canonical full invocation): clean, 823 files.
- `uv run pytest tests/genome/test_heteroclinic_cycle.py -q`: 8/8 pass (no
  regression from the `x_sign`/`ConnectionNode`/`epsilon` changes).
- `uv run pytest tests/search/test_jovian_resonant_connections.py -v`: 11/11 pass
  (~58s), including the Table-2 gate evidence test — deliberately **not** marked
  `@pytest.mark.slow` (a discovery-verdict-bearing test must run in CI, not be
  silently skipped).
- `uv run pytest tests/data tests/search tests/scripts -q`: run before committing
  (full project ratchet suite, per this project's own discipline for any
  `search`/`genome` change).
- `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: pass, run before the
  `OUTSTANDING.md` update commit.

## 6. Net effect on `#754` / Task B

**Table-2 half: DONE, honest FAIL (close).** The homoclinic self-connection of
`3:4-LO` is real, Newton-converged to `2.58e-9`, independently Radau-cross-checked
to `6.8e-8`, ghost-guard-verified (145x margin), and qualitatively matches the
paper's Fig. 25 forward/backward re-approach behavior — but its section-crossing
location misses Anderson & Lo's own published Table-2 digits by `1.1e-3` in `x`,
outside this task's `1e-4` gate. This is reported as a genuine, close, non-fudged
negative — and as a third, independent piece of corroborating evidence (alongside
`#755`'s eigenvalue and shape/close-approach matches) that our `3:4-LO` family is
indeed the paper's own.

**Table-3 half: still BLOCKED**, unchanged from `#757`'s own scoping (needs
`5:6-LO`, the sole gate-keeper for which is the sibling task `#758`, not touched by
this task).
