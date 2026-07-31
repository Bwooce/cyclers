# `#767` — Saturn-Titan 3:4 resonant-orbit homoclinic self-connection

**Task:** `#767`, the Saturn-Titan Task-B analog, mirroring `#754`'s Jupiter-Europa Table-2
homoclinic self-connection build and `#766`'s own self-consistency-gated build at an energy
with no published state to gate against — but here for the newly-confirmed (`#765`) Saturn-Titan
3:4 resonant orbit at Vaquero 2013's own sourced Jacobi constant `C = 3.010000`.

**Critical honesty framing (read this before the numbers, exactly as `#766`'s own note states):**
Vaquero 2013's own thesis describes a homoclinic connection of the 3:4 resonant orbit (Sec. 4.3.1,
Fig. 4.9) as a **figure only** — no state table, unlike Anderson & Lo 2011's own digit-grade
Table 2 for the analogous Jupiter-Europa case. There is nothing published to reproduce here.
Every number below is therefore **self-consistency** evidence — Newton residual, ghost-guard
margin, independent Radau cross-check, forward/backward re-approach — never a reproduction of a
published number. A qualitative visual comparison against Fig. 4.9 is offered at the end as soft
corroboration only, not a gate.

---

## Verdict (read this first)

**A genuine, non-trivial homoclinic self-intersection EXISTS at Vaquero's own `C = 3.010000` and
was found — in fact FOUR independent ones, via three different `(branch_u, branch_s)` sign
combinations, including a genuine reflection-symmetric mirror pair.** All four Newton-converged
to residuals `<1e-8` (three below `1e-9`), independently Radau-cross-checked to `<2e-7`, and land
`219x`–`307x` past the ghost-guard threshold (`GHOST_GUARD_DELTA = 1e-3`) — a real, non-delicate
margin. Forward/backward re-approach is honest and asymmetric (see Sec. 5): the backward leg is
tight (`<5e-8` for all four), the forward leg is looser (`0.011`–`0.080`, still far below the O(1–2)
trajectory scale) — explained directly by this orbit's own strong instability (`|lambda|~2129.8`,
the strongest of any orbit this task chain has built a self-connection for).

---

## 1. The orbit itself (from `saturn_titan_resonant_families.py`, `#765`)

Read directly from `recover_table41_candidate("3:4")`'s own return value this task (not from
`#765`'s note by hand):

```
x0       = 1.0301662783998498
ydot0    = 0.085643311781831
period   = 26.140797240249157   (period/2pi = 4.1604370971486295)
jacobi   = 3.01 exactly
lambda   = 2129.807723811082    (real saddle, is_real_unstable=True)
```

`#765`'s own gate: `eigenvalue_rel_err = 1.069e-6`, `period_rel_err = 1.727e-4`, both well inside
`TABLE41_EIGENVALUE_GATE_REL_TOL`/`TABLE41_PERIOD_GATE_REL_TOL` (`1e-3`) — `passed = True`,
re-confirmed fresh this task via `saturn_titan_resonant_connections.build_34_node()`.

`ResonantNode.from_candidate` (reused UNCHANGED from `jovian_resonant_connections.py`, confirmed
fully system-agnostic per `#764`'s own finding) recomputes the saddle Floquet pair from this
stored `(x0, ydot0, period)` with no schema change needed.

## 2. Why this task built a NEW sibling module, not an extension of `jovian_resonant_connections.py`

Read that module in full this task. Its own truly system-agnostic machinery
(`ResonantNode`, `own_section_points`, `correct_connection` itself via `heteroclinic_cycle`,
`_full_state_crossing`, `_ghost_distance`, `HomoclinicReapproachResult`) is reused directly,
unchanged. But the module's own TOP-LEVEL constants and convenience builders (`TABLE2_STATE`,
`ANDERSON_LO_C_FLYBY`, `SECTION_X_SIGN=-1`, `jupiter_europa_system`, `build_3_4_lo_node`, etc.)
are Jupiter-Europa-specific — importing them for a Saturn-Titan orbit would be a false-cognate
reuse, not genuine sharing. New sibling module:
`src/cyclerfinder/search/saturn_titan_resonant_connections.py`, exactly mirroring how `#765`'s
own `saturn_titan_resonant_families.py` is a thin sibling of `jovian_resonant_families.py`
rather than an edit to it.

## 3. The section convention had to be independently derived — NOT copied from Anderson & Lo

Anderson & Lo's own one-sided section (`{y=0, x<0, ydot>0}`) is specific to THEIR 3:4-LO orbit's
own geometry (IC at negative `x`). This task's own 3:4 orbit sits at **positive** `x`
(`x0=1.0301663`, confirmed directly this task) — reusing `x_sign=-1` verbatim would exclude the
orbit's own IC from `own_section_points`'s reference set entirely, **silently breaking the ghost
guard** (an empty reference set makes `_ghost_distance` return `inf` unconditionally, so nothing
would ever be rejected as a trivial self-shadow — a real, caught-before-shipping bug, not
hypothetical: an early attempt at exactly this using `x_sign=+1` only found essentially no
crossings past `k=1` within many orbital periods, since this orbit visits the single-sided
`{x>0, ydot>0}` quadrant only once per period).

Direct inspection of this orbit's own `{y=0}` crossings over one period (this task) found 4
crossings, at exactly TWO **perpendicular** (`xdot=0`) points — `x0=+1.0301663` (the IC) and
`x=-1.3666368` (the corrector's own half-period return target,
`_HALF_CROSSINGS['3:4']=2` in `saturn_titan_resonant_families.py`) — both `ydot>0` — plus a
non-perpendicular mirror pair at `x=+1.1064570` (`ydot<0`). This module's own section is
`{y=0, ydot>0}`, **x UNRESTRICTED** — the natural generalization keeping BOTH of the orbit's own
perpendicular reference points (`own_section_points` unions `jrc.own_section_points` at
`x_sign=+1` and `x_sign=-1`). `correct_connection`'s own `x_sign_u`/`x_sign_s` kwargs already
support `None` (unrestricted) natively — confirmed by reading `heteroclinic_cycle.py` in full
this task — no change to any existing file was needed.

Under this convention, a direct trace of the unstable manifold (`tau=0`, `branch=+1`) found
dozens of genuine, well-separated crossings within a 10-period horizon (vs. essentially none
under the single-sided restriction) — this section choice is what made the scan below tractable.

## 4. The scan: FOUR genuine converged hits found across THREE branch-sign combinations

Scanned `branch_u, branch_s in {+1,-1}` (all 4 combinations) with `k_u, k_s` mostly spanning
`1..6` (`#754`'s original range, per the dispatch note's own instruction to start there),
`max_time_factor=3.0` (this orbit's eigenvalue, `|lambda|=2129.8`, sits closer to `#754`'s own
`C_flyby` regime, `|lambda|=1036`, than `#766`'s weaker `|lambda|=54.6` case — and indeed
`max_time_factor=3.0` sufficed here, no widening needed). Light triage settings first
(`scan_n=6, max_iter=10, tol=1e-6`, `#766`'s own two-stage approach), refined afterward
(`scan_n=12, tol=1e-9, max_iter=60, fd_step=1e-7` — `fd_step=1e-7` empirically found tighter than
the module's own default `1e-6`, see Sec. 6) for every combination that showed promise.

| branch_u | branch_s | k_u | k_s | residual | crossing (x, xdot) | ghost_distance (×guard) |
|---|---|---|---|---|---|---|
| **+1** | **+1** | **5** | **5** | `7.18e-10` | `(-1.673201, -7.7e-10)` | `0.3066` (**307×**) |
| −1 | −1 | 4 | 5 | `9.50e-10` | `(0.844031, -0.114813)` | `0.2187` (219×) |
| −1 | −1 | 5 | 4 | `7.44e-10` | `(0.844031, +0.114813)` | `0.2187` (219×) |
| +1 | −1 | 4 | 4 | `1.67e-10` | `(-1.126473, 0.132115)` | `0.2741` (274×) |

The `(-1,-1,4,5)`/`(-1,-1,5,4)` pair are genuine reflection-symmetric mirrors of each other
(equal `x`, `xdot` sign-flipped, matching to `<1e-6`) — two independent intersections related by
the CR3BP's own time-reversal symmetry, not the same point found twice (see
`test_mirror_pair_reflection_symmetry`). The `(+1,+1,5,5)` hit sits ON the symmetry axis
(`xdot≈0` to `7.7e-10`) — structurally the SAME point-type as Anderson & Lo's own published
Table 2 state and `#766`'s own primary hit — chosen as this task's **PRIMARY** result for that
reason, and because it has the tightest evidence bundle of the four (see Sec. 5/6). The
`(+1,-1,4,4)` hit is a fourth, independently-found, off-axis intersection via yet another
branch-sign combination.

Branch pair `(-1,+1)` was ALSO scanned across `k=1..6` and yielded no genuine hit (every
combination either converged onto/near the orbit's own trivial section points, or failed to
converge at all under light-triage settings) — an honest negative for that one quadrant, not
suppressed.

**This is a richer scan result than either `#754` (exactly one hit) or `#766` (three hits,
one primary + one mirror pair)** — three of the four `branch`-sign combinations tried here
produced genuine intersections, versus one for `#754` and effectively 1.5 (two combos out of
four scanned) for `#766`.

## 5. Verification of all four hits

**Ghost guard**: the orbit's own two qualifying section points (this section's own convention,
Sec. 3) are `(1.0301663, 0)` (the IC) and `(-1.3666368, 0)` (the half-period point). All four
found crossings sit `0.2187`–`0.3066` away — `219×`–`307×` the `GHOST_GUARD_DELTA=1e-3`
threshold, a real margin by a wide berth, comparable to or exceeding `#754`'s own `145×` and
well past `#766`'s own `37×`–`67×`.

**Independent Radau cross-check** (`assemble_cycle`/`crosscheck_cycle`, `rtol=atol=1e-11`):
all four agree between DOP853 and Radau to `1.92e-8`–`4.67e-7` — comfortably inside the
mandated `<=1e-6` for every hit.

| hit | residual | ghost_distance | Radau independent_residual |
|---|---|---|---|
| PRIMARY (+1,+1,5,5) | `7.18e-10` | `0.3066` | `1.92e-8` |
| MIRROR_A (−1,−1,4,5) | `9.50e-10` | `0.2187` | `1.98e-7` |
| MIRROR_B (−1,−1,5,4) | `7.44e-10` | `0.2187` | `2.00e-7` |
| TERTIARY (+1,−1,4,4) | `1.67e-10` | `0.2741` | `4.67e-7` |

## 6. Forward/backward re-approach: honest, asymmetric, and explained

This orbit's own unstable eigenvalue (`|lambda|~2129.8`) is the **strongest instability** of any
orbit this task chain has built a homoclinic self-connection for (`#754`'s `C_flyby`: `1036`;
`#766`'s `C=3.0041`: `54.6`). Direct empirical testing this task found `jrc.homoclinic_reapproach_check`'s
own hardcoded `rtol=atol=1e-12` crossing-detection integrator (not threaded through from its own
public kwargs) left the recomputed full states at the unstable- and stable-leg crossings
(`y_u`, `y_s`) differing by `~1e-6` for this orbit — an order of magnitude looser than `#766`'s
own `<1e-8` at that same tolerance, because this orbit's own crossing search integrates over a
comparable elapsed time under MUCH stronger exponential growth. This task's own
`saturn_titan_resonant_connections.homoclinic_reapproach_check` threads a tighter
`rtol=1e-13, atol=1e-14` through EVERY internal integration (the manifold seeding, the crossing
search, AND the final re-propagation) — recovering roughly an order of magnitude improvement
(e.g. the PRIMARY hit's own `forward_distance` improved from `0.42` at the untightened pipeline
to `0.011`) — reported honestly, not tightened further by fiat.

| hit | `t_u` (periods) | `t_s` (periods) | backward_distance | forward_distance |
|---|---|---|---|---|
| PRIMARY (+1,+1,5,5) | `2.715` | `−2.715` | `4.13e-8` | `0.0112` |
| MIRROR_A (−1,−1,4,5) | `1.875` | `−2.714` | `9.27e-9` | `0.0366` |
| MIRROR_B (−1,−1,5,4) | `2.714` | `−1.875` | `1.05e-8` | `0.0367` |
| TERTIARY (+1,−1,4,4) | `2.174` | `−2.362` | `7.40e-10` | `0.0797` |

**Why the asymmetry is real, not a bug**: `backward_distance` propagates the SAME leg's own
recovered crossing state back to its OWN originating seed — close to an exact numerical
time-reversal of the identical integration, so it is consistently tight (`<5e-8` for all four
hits, tighter than `#766`'s own `5.6e-8` primary result in three of four cases). `forward_distance`
propagates that SAME state (derived from the unstable leg alone) FORWARD to independently
reproduce the STABLE leg's own separately-computed seed — a genuine cross-leg consistency test
that inherits the corrector's own finite full-state residual floor (`~1e-6`–`1e-9`, tighter than
the raw `(x,xdot)` Newton residual alone since ydot enters too), amplified by
`exp((ln|lambda|/period) * elapsed_time)`. For the PRIMARY hit this amplification factor is
`exp(0.293 * 71.0) ≈ 7×10^7` — an order-of-magnitude looser starting residual than `#766`'s own
case, amplified by a COMPARABLE exponential factor, naturally produces a materially looser
(but still `<1.2` orders of magnitude below the O(1–2) trajectory scale) `forward_distance`.
This is the honest, explained consequence of this orbit's own much stronger instability — not
fudged, not hidden. All four `forward_distance` values remain far below the trajectory's own
O(1–2) nondimensional scale (`1–8%` of it) — a genuine, non-coincidental re-approach, just not
machine-tight.

## 7. Qualitative Fig. 4.9 corroboration (soft evidence only, not a gate)

Vaquero 2013 Fig. 4.9 (p.112) shows the 3:4 orbit's own homoclinic manifold trajectory as a
looping, multi-lobed excursion that departs the orbit near one side and re-approaches from the
other, consistent qualitatively with a transversal self-intersection existing at moderate
`x`-range excursions similar in scale to this task's own found crossings (`x` in roughly
`[-1.7, 0.85]`, vs. the orbit's own IC/half-period points at `x = 1.03`/`-1.37`). This is offered
as soft, qualitative corroboration ONLY — Fig. 4.9 has no digit-grade state to compare against,
per the critical honesty framing above, so this is not and cannot be a quantitative gate.

## 8. Code delivered

* `src/cyclerfinder/search/saturn_titan_resonant_connections.py` (new sibling module):
  `own_section_points` (x-unrestricted union), `build_34_node`, `HomoclinicCandidate`,
  `find_homoclinic` (residual-ranked, no published target), `homoclinic_reapproach_check`
  (tighter end-to-end integration precision than `jrc`'s own version — see Sec. 6). No changes
  to any existing file — `jovian_resonant_connections.py`, `heteroclinic_cycle.py`, and
  `saturn_titan_resonant_families.py` are all reused unchanged, confirming `#764`'s own finding
  that this whole pipeline is system-agnostic.
* `tests/search/test_saturn_titan_resonant_connections.py`: 23 tests covering the section
  convention derivation, all four known hits' convergence/ghost-margin/Radau-cross-check/
  forward-backward-re-approach, the mirror-pair symmetry, `find_homoclinic`'s own scan/ranking
  behavior, and an honest-empty-scan regression. None marked `@pytest.mark.slow` (a
  discovery-verdict-bearing result must run in CI).

## 9. Verification

* `uv run ruff check` / `ruff format --check` on both new files: clean.
* `uv run mypy src tests` (canonical full invocation): clean, 827 source files.
* `uv run pytest tests/search/test_saturn_titan_resonant_connections.py -v`: 23/23 pass (~47s).
* `uv run pytest tests/data/test_outstanding_structure.py tests/data/test_outstanding_header_body_consistency.py -q`:
  run before the `OUTSTANDING.md` update commit (see commit history for pass/fail status
  recorded at commit time).

## 10. Net effect on `#767`

**DONE — a genuine, honestly-evidenced positive result.** Four independent, ghost-guard-passed,
Newton-converged (`<1e-8` residual, three `<1e-9`), independently Radau-cross-checked (`<5e-7`)
homoclinic self-intersections of the Saturn-Titan 3:4 resonant orbit exist at Vaquero's own
`C=3.010000`, including a genuine reflection-symmetric mirror pair. Forward/backward re-approach
is honestly reported as asymmetric (backward tight, forward looser but still well below the
trajectory's own O(1–2) scale), explained directly by this orbit's own strong instability
(the strongest of any self-connection this task chain has built). This is self-consistency
evidence only — no published Table-2-style state exists to gate against (Fig. 4.9 is a figure
only) — framed explicitly as such throughout, per this task chain's own honesty discipline.
Unblocks `#768` (the 3:4↔6:5 resonant-chain reproduction), which needed this connection built
first.
