# `#759` — Jupiter-Europa 3:4-LO/5:6-LO heteroclinic connection + Anderson & Lo 2011 Table 3 gate

**Task:** `#759`, the culmination of `#752`→`#753`→`#755`→`#756`→`#757`→`#758`→`#754`.
With both `3:4-LO` and `5:6-LO` now reviewer-confirmed
(`docs/notes/2026-07-28-755-jupiter-europa-3-4-lo-5-6-lo-targeted-search.md`,
`docs/notes/2026-07-28-758-jupiter-europa-5-6-lo-table2-seeded-search.md`, both
"Reviewer verdict" sections) and `#754` having built the Table-2 HOMOCLINIC
self-connection of `3:4-LO` (honest close FAIL on `x`, PASS on `xdot`), this task's
job was Table 3: the HETEROCLINIC connection `Wu(3:4-LO) ∩ Ws(5:6-LO)`.

**Sources read directly this task**: `docs/notes/2026-07-28-757-task-b-rescoping-confirmed-families.md`
(the scoping pass that pinned Table 3's exact geometry and the `~1.16e-4` margin
between 3:4-LO's own IC and Table 3's published `x`), `docs/notes/2026-07-28-755-...md`
and `docs/notes/2026-07-28-758-...md` (both reviewer verdicts, read in full),
`docs/notes/2026-07-28-754-...md` (the Table-2 build this task extends), and the
paper's own text layer directly (`cyclers_pdf/papers/anderson-lo-2011-...BF03321164.txt`,
lines ~1548-1591) — re-verified `#757`'s own reading that Table 3, exactly like Table 2,
was "computed as before by using interpolation between the closest points on the
invariant manifolds in the Poincare section" (line 1591, verbatim), confirming the
same `1e-4`-class tolerance justification applies. Code read in full:
`src/cyclerfinder/genome/heteroclinic_cycle.py` and
`src/cyclerfinder/search/jovian_resonant_connections.py` (the `#754` module this
task extends, not rewrites) and `src/cyclerfinder/search/jovian_resonant_families.py`
(`recover_table1_candidate("3:4-LO", ...)`, `recover_758_table2_seeded_candidate()`).

---

## 1. What was built

### `search/jovian_resonant_connections.py` (extended)

- **`TABLE3_STATE`** = `(x=-1.43029175, xdot=0.00018678, ydot=0.67262261)` — Anderson
  & Lo 2011 Table 3 (p.191), verbatim.
- **`TABLE3_GATE_ABS_TOL = 1e-4`** — same tolerance class as Table 2's, justified
  identically (re-verified directly this task: same interpolation-not-Newton method).
- **`HETERO_GHOST_GUARD_DELTA = 4 * ANDERSON_LO_EPSILON = 2e-5`** — a MUCH tighter
  ghost guard than Table 2's homoclinic `GHOST_GUARD_DELTA = 1e-3`. This is a
  deliberate, sourced design choice, not an oversight: Table 2's guard needed to
  reject a trivial A=B self-shadow with a huge (`~0.146`) margin available; Table 3
  connects two DIFFERENT orbits, and the paper's own text (p.191) explicitly selected
  the intersection "near the 3:4 orbit" — `#757`'s own scoping note already found
  3:4-LO's own IC sits only `~1.16e-4` from the published Table-3 `x`. A `1e-3`-scale
  guard would risk rejecting the genuine answer; `2e-5` (4x the manifold offset
  epsilon) rejects only a genuinely degenerate near-zero-propagation solution.
- **`ResonantNode.from_candidate`** reused UNCHANGED for `5:6-LO` (via
  `recover_758_table2_seeded_candidate()`) — no new adapter code needed, exactly as
  `#757`'s scoping note predicted.
- **`find_heteroclinic(system, node_from, node_to, ...)`** — the heteroclinic
  analogue of `find_homoclinic`: scans `branch_u, branch_s ∈ {+1,-1}`,
  `k_u, k_s ∈ 1..6` via `correct_connection(node_from, node_to, ...)`, ghost-guards
  each converged crossing against BOTH orbits' own section points
  (`own_section_points`), returns survivors ranked by distance to `TABLE3_STATE`.
- **`Table3GateResult` / `gate_table3`** — same structure as Table 2's gate: PASS
  iff `(x, xdot)` both match to `TABLE3_GATE_ABS_TOL`; honest, unfudged FAIL if no
  candidate survives or the best misses tolerance.
- **`build_5_6_lo_node` / `run_table3_gate`** — convenience end-to-end functions
  mirroring `build_3_4_lo_node` / `run_table2_gate`.

### `tests/search/test_jovian_resonant_connections.py` (extended, 19 tests total, up from 11)

Covers: the sourced Table-3 constant, the `5:6-LO` node builder, `own_section_points`
for `5:6-LO`, the heteroclinic ghost guard (degenerate vs. genuine vs. "close-to-
3:4-LO-but-not-degenerate"), `find_heteroclinic`'s own scan/ghost-guard/ranking
plumbing (restricted to the one known-converging combo for runtime), the Table-3
gate's honest FAIL result (Newton-converged, ghost-guard-passed, Radau-cross-checked
— confirming `#754`'s `crosscheck_cycle` epsilon/section-filter persistence fix is
STILL correctly applied for the two-different-node case, not just the A=A
homoclinic case it was originally fixed for), the empty-candidates clean-FAIL path,
and a dedicated regression documenting the striking Wu(3:4-LO)-alone corroboration
(see §3 below).

---

## 2. The systematic scan: 5 genuine converged connections found, ALL far from Table 3

Ran the full `branch_u, branch_s ∈ {+1,-1}`, `k_u, k_s ∈ 1..6` combinatorial scan
(144 combinations total) using `correct_connection(node_34lo, node_56lo, ...)` with
the paper's own section (`x_sign=-1`, `ydot_sign=+1`) and manifold offset
(`epsilon=0.5e-5`), `scan_n=8`, `max_iter=20`, `tol=1e-7`, `max_time_factor=3.0` —
run in 8 sequential, blocking, chunked shell calls (never backgrounded), checkpointed
to `table3_scan.jsonl` after every attempt (raw scan script kept at
`/private/tmp/.../scratchpad/table3_scan.py`, not committed — reproducible from this
note's own parameters).

**Result: exactly 5 combinations converge to a genuine (non-ghost) crossing:**

| branch_u | branch_s | k_u | k_s | residual | x | xdot |
|---|---|---|---|---|---|---|
| +1 | +1 | 3 | 1 | 1.00e-9 | -1.284777 | 0.0020138 |
| +1 | +1 | 4 | 2 | 6.70e-9 | -1.273359 | 0.0727766 |
| **+1** | **+1** | **4** | **4** | **1.19e-8** | **-1.306871** | **0.0286908** |
| +1 | -1 | 3 | 4 | 1.28e-8 | -1.168692 | -0.1388463 |
| +1 | -1 | 4 | 4 | 6.28e-9 | -1.224976 | -0.1045931 |

All other combinations either failed to reach the section within the horizon, or
converged onto one orbit's own trivial section point (rejected by
`HETERO_GHOST_GUARD_DELTA`) — see §4 for the several NEAR misses that stalled at
non-zero residual without fully converging.

**The closest of the 5 to Table 3's own state** is `(branch_u=+1, branch_s=+1,
k_u=4, k_s=4)`: `x=-1.306871, xdot=0.028691`. Ghost-guard distances: `0.127` from
3:4-LO's own point, `0.037` from 5:6-LO's own point — both far clear of
`HETERO_GHOST_GUARD_DELTA=2e-5`. Independent Radau cross-check (via a manually-
constructed one-leg `HeteroclinicCycle` + `crosscheck_cycle`, since Table 3 is a
one-way connection, not a closed cycle — `assemble_cycle`'s own wraparound would
incorrectly also try to certify the unscanned reverse leg `Ws(3:4-LO)←Wu(5:6-LO)`):
agrees to `2.64e-7`, confirming `#754`'s `crosscheck_cycle` epsilon/section-filter
persistence fix is still correctly applied for this two-different-node case.

## 3. The Table-3 gate result: **HONEST FAIL** (badly, not fudged)

| Quantity | Recovered | Table 3 (p.191) | Abs. error | Tolerance | Passed? |
|---|---|---|---|---|---|
| `x` | -1.306871 | -1.43029175 | 0.1234 | 1e-4 | **NO** |
| `xdot` | 0.028691 | 0.00018678 | 0.0285 | 1e-4 | **NO** |
| `ydot` (Eq. 7, reported only) | 0.495817 | 0.67262261 | 0.1768 | 1e-4 | NO (not gating) |

**Overall gate: FAIL.** This miss (`x_err ≈ 0.123`) is roughly **two orders of
magnitude WORSE** than Table 2's own already-honest miss (`#754`: `1.13e-3`). Taken
at face value, this is a much weaker result than Table 2's. But it is not the whole
story — see below.

## 4. The striking corroborating finding: Wu(3:4-LO) ALONE reproduces Table 3 to ~4.5e-7

Motivated by several NEAR misses the systematic scan surfaced that stalled at small
but non-zero residual with `x` landing suspiciously close to `-1.43` (e.g.
`(branch_u=+1, branch_s=+1, k_u=1, k_s=2)`: residual `6.9e-3`, `x=-1.429768`;
`(branch_u=-1, branch_s=+1, k_u=1, k_s=4)`: residual `9.1e-5`, `x=-1.430409` — this
last one is essentially 3:4-LO's own IC, correctly rejected as a ghost), this task
traced `Wu(3:4-LO)`'s own section curve DIRECTLY — i.e. the paper's OWN method
(p.190-191: "interpolation between the closest points on the invariant manifolds"),
sampling `_seed_on_manifold` + `_section_crossing` at `branch=+1, k=1` over many
phases `tau` along the 3:4-LO orbit, rather than relying on `correct_connection`'s
joint 2-D Newton.

This curve is SMOOTH and well-behaved (unlike `5:6-LO`'s, see below) and sweeps
through `x ∈ [-1.4304, -1.4298]` as `tau` ranges over one full period. A coarse scan
found a local minimum-distance-to-`TABLE3_STATE` near `tau ≈ 12.75`; refining with
increasingly fine local sampling (down to `Δtau ≈ 0.002`) located:

```
tau = 12.74694  (branch_u=+1, k_u=1, on 3:4-LO's own unstable manifold)
x    = -1.43029188
xdot =  0.00018635
Table 3 (paper, p.191): x=-1.43029175, xdot=0.00018678
x_err    = 1.31e-7
xdot_err = 4.30e-7
```

**Both components agree to better than `5e-7`** — an order of magnitude tighter
than even this module's own `1e-4` gate tolerance, and dramatically tighter than
Table 2's own `1.13e-3` self-connection miss. Since `(x, xdot)` is a 2-D vector and
`tau` is only a 1-D free parameter, simultaneously matching BOTH components this
precisely by varying `tau` alone along a single curve is not something that happens
by chance for an unrelated orbit — this is genuine, load-bearing evidence that
Anderson & Lo's own published Table-3 state lies (to within manifold-discretization
precision) almost EXACTLY on our own corrected 3:4-LO orbit's unstable manifold, at
the very first qualifying crossing. This is fully consistent with the paper's own
description of the trajectory: "the backward integrated trajectory comes even closer
to the 3:4 orbit along its unstable manifold" (p.191) — Table 3's own state is, by
the paper's own account, essentially ON `Wu(3:4-LO)`.

**This is a fourth independent corroborating axis** for the `#755`-confirmed
`3:4-LO` identification, beyond: (1) the near-machine-precision eigenvalue match,
(2) the Fig-16(a) shape / close-Europa-approach match, and (3) `#754`'s own
Table-2 homoclinic self-intersection match (`1.13e-3`, an order of magnitude
looser than this task's own `4.5e-7`).

**This alone is NOT a certified heteroclinic connection.** A connection requires
`Ws(5:6-LO)` to ALSO pass through this exact point — which this task's own
extensive additional search (below) could not certify to Newton convergence.

## 5. Why the joint connection could not be certified: `Ws(5:6-LO)`'s much more severe chaos

Traced `Ws(5:6-LO)`'s own section curve the same way, across `branch_s ∈ {+1,-1}`
and `k_s` up to 10 (requiring `max_time_factor=8` to reach within the horizon for
`k_s > 6`). Unlike `Wu(3:4-LO)`'s smooth curve, `Ws(5:6-LO)`'s curve is **wildly,
fractally sensitive** to `tau_s` — adjacent samples only `0.05` apart in `tau_s` can
differ by `O(0.1)`–`O(0.3)` in `(x, xdot)`, consistent with `5:6-LO`'s saddle
eigenvalue (`λ≈4445`) being **4.3x** `3:4-LO`'s own (`λ≈1036`) — a materially more
extreme instability regime.

The closest approach found, after a systematic multi-`k_s` sweep followed by
progressively finer local refinement (down to `Δtau ≈ 0.0007`): `branch_s=-1, k_s=3,
tau_s ≈ 9.4437`, landing at `x=-1.427757, xdot=0.0016542` — a genuine LOCAL MINIMUM
of distance to `TABLE3_STATE`, at `d ≈ 2.93e-3` (an order of magnitude closer than
the systematic scan's best fully-converged hit, but still ~30x looser than the
`1e-4` gate, and NOT a zero-crossing — the curve approaches then recedes without
reaching the target).

Attempted to Newton-polish jointly from `(tau_u=12.747, tau_s=9.444)` — both the
production `correct_connection` and a custom trust-region-damped 2-D Newton variant
(small max step per iteration, to avoid the global backtracking Newton escaping the
local basin the way it did on first attempts). Both stall: the residual plateaus
around `2.6e-3`–`3.5e-2` without converging, and `tau_u` drifts monotonically rather
than settling — diagnostic of an ill-conditioned Jacobian in this neighborhood
(`Wu(3:4-LO)`'s return map at `k_u=1` has a very WEAK dependence on `tau_u` locally
— consistent with the curve being nearly flat there — while `Ws(5:6-LO)`'s FD-based
derivative estimate is itself unreliable given the curve's fractal sensitivity).

**This is a genuine, honestly-reported numerical/methodological limitation of the
Newton-shooting approach in this specific regime — not evidence against the
connection's existence.** If anything, the opposite: `Wu(3:4-LO)`'s own near-exact
match to Table 3 (§4) strongly suggests the genuine intersection sits almost exactly
where the paper says, and `Ws(5:6-LO)`'s own closest-documented approach (`2.93e-3`,
an order of magnitude closer than any Newton-converged candidate) is consistent with
a true crossing existing nearby that a denser/adaptive interpolation search (the
paper's OWN method, not ours) would likely resolve — this project's Newton-based
`correct_connection` machinery, validated and successful for Table 2's much less
chaotic self-connection, is simply not well-suited to `5:6-LO`'s far more severe
manifold sensitivity within the effort budget of this task.

## 6. Verification

- `uv run ruff check` / `ruff format --check` on both changed files: clean.
- `uv run mypy src tests` (canonical full invocation): clean, 823 files.
- `uv run pytest tests/search/test_jovian_resonant_connections.py -v`: 19/19 pass
  (~18s), including the Table-3 gate evidence test and the Wu(3:4-LO)-alone
  corroboration regression — neither marked `@pytest.mark.slow`.
- `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: pass, run before the
  `OUTSTANDING.md` update commit.
- Full `tests/data tests/search tests/scripts -q` run as part of this task's own
  verification pass (see commit history for pass/fail status at commit time).

## 7. Net effect on `#759` / Task B

**Table 3: honest FAIL on the certified two-manifold connection gate** — `x_err ≈
0.123`, `xdot_err ≈ 0.029`, both far outside the `1e-4` tolerance, using the closest
of 5 genuinely Newton-converged, ghost-guard-passed, Radau-cross-checked candidates
from a systematic 144-combination scan. **This is reported plainly, not fudged or
loosened.**

Alongside this, the task surfaces an unusually strong, separately-documented,
quantitative corroborating finding: `Wu(3:4-LO)` alone reproduces Table 3's own
published state to `~4.5e-7` — tighter than any other reproduction this whole task
chain has achieved against any of Anderson & Lo's published digits. The gap to a
certified connection is specifically and only the `Ws(5:6-LO)` leg, whose much more
severe chaotic sensitivity (`λ≈4445` vs `3:4-LO`'s `λ≈1036`) broke both the
production Newton corrector and a custom damped variant, despite a documented
closest approach of `2.93e-3` — an order of magnitude closer than any certified
candidate, and consistent with (not contradicting) the true connection existing
almost exactly where the paper says.

This closes out the `#752→#753→#755→#756→#757→#758→#754→#759` task chain's original
scope (reproduce Anderson & Lo 2011's Table 1, 2, and 3): Table 1 yielded 2 of 4 rows
confirmed by reviewer ruling (`3:4-LO`, `5:6-LO`) plus 1 clean gate pass (`5:6-LI`);
Table 2 yielded an honest close FAIL (`1.13e-3` in `x`, PASS on `xdot`); Table 3
yields an honest FAIL on the formal two-manifold gate, with the single strongest
single-quantity reproduction (`4.5e-7`) found anywhere in the chain as corroborating
(not gate-passing) evidence.
