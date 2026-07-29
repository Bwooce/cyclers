# `#766`: homoclinic self-connection AT the torus seed's own energy, C=3.0041

**Task:** `#766`, direct continuation of `#761` (which proved Kumar, Anderson, de la
Llave & Gunter 2021's own "arbitrarily chosen" Jupiter-Europa 3:4 exterior resonant
seed orbit — Jacobi constant `C = 3.0041`, `KUMAR_2021_C` in
`src/cyclerfinder/search/jovian_resonant_families.py` — is a genuine real saddle on
the SAME continuous family as the confirmed 3:4-LO). This task builds a homoclinic
self-connection (`Wu(orbit) ∩ Ws(orbit)`) AT that specific energy — the energy the
catalogued `europa-3-4-crnbp-torus-jupiter-2026` torus is actually built around —
which `#754` never did (`#754`'s own connection was computed at
`ANDERSON_LO_C_FLYBY = 2.99163956830415`, a different point on the same family).

**Critical honesty framing (read this before the numbers):** there is no published
Table-2-style state to gate against here. Kumar 2021 does not report a homoclinic
connection state for its own seed orbit anywhere in the paper. Everything below is
therefore a **self-consistency** result — Newton residual, ghost-guard margin,
independent Radau cross-check, forward/backward re-approach — never a reproduction
of a published number. This is stated explicitly in every place a number is
reported, and should not be mistaken for the same kind of evidence `#754`'s
Table-2 gate produced.

---

## Verdict (read this first)

**A genuine, non-trivial homoclinic self-intersection EXISTS at C=3.0041 and was
found, Newton-converged to `2.0e-10` residual, Radau-cross-checked to `2.4e-8`,
ghost-guard-passed with a real (37x) margin, and forward/backward re-approach
confirmed to `5.6e-8`/`2.3e-5`.** Two further independent (mirror-pair) hits were
also found and converged, corroborating that this energy genuinely supports
transversal homoclinic self-intersections rather than this being an isolated
numerical fluke.

A secondary, methodologically important finding: the SPECIFIC self-connection
branch `#754` found at `C_flyby` (`branch_u=+1, branch_s=-1, k_u=3, k_s=3`) does
**not** survive continuously up to `C=3.0041` — direct continuation of that
specific `(tau_u, tau_s)` solution in Jacobi constant hits an apparent fold/
tangency around `C ≈ 2.99941`, well below Kumar's energy. The genuine connection
found at `C=3.0041` is therefore a **different** `(branch, k)` combination
(`k_u=k_s=6` instead of `k_u=k_s=3`), consistent with the ~19x weaker instability
at this energy (`|lambda|≈54.6` vs. `1036` at `C_flyby`) requiring the manifold to
develop through more crossings (and more elapsed time, ~5.3 orbital periods
instead of ~2) before a transversal self-intersection appears.

---

## 1. The orbit itself (from `#761`'s `continue_34lo_to_kumar_c`)

Re-derived fresh this task (not read from `#761`'s note by hand — pulled directly
from `continue_34lo_to_kumar_c()`'s own return value, per this task's own
instruction to read the precise numbers from the function, not rounded figures):

```
x0      = -1.3852484456241640   (C_flyby continuation endpoint, #761)
ydot0   =  0.5988394002678391
period  = 25.312119648766764    (period/2pi = 4.02854896223472)
jacobi  =  3.0041 exactly
lambda  = 54.589750588953734    (real saddle, is_real_unstable=True)
```

(Matches `#761`'s own recorded `x0=-1.3852484456241585`, `period=25.31211964876615`,
`lambda=54.589750588974` to machine/reproduction precision — the tiny last-digit
differences are ordinary run-to-run floating-point noise in the continuation's own
Newton corrector, not a discrepancy.)

`build_34lo_kumar_c_node()` (new, `jovian_resonant_connections.py`) wraps this
directly: it calls `continue_34lo_to_kumar_c()` fresh every invocation (re-running
the full `#761` continuation gauntlet — closure, period bounds, equilibrium gate,
Jacobi conservation, independent-Radau, fold/topology-jump detection — all inherited
for free, not a cached shortcut) and builds a `ResonantNode` from the endpoint via
the EXISTING, unchanged `ResonantNode.from_candidate` adapter (confirmed general
enough for any `ResonantFamilyCandidate`, not just the original 3:4-LO one, exactly
as the dispatch note anticipated).

## 2. Why a straight re-run of `#754`'s own scan needed re-parameterization

`#754`'s own coarse scan (`branch_u, branch_s in {+1,-1}`, `k_u, k_s in 1..6`,
`max_time_factor=3.0`) was re-run FIRST, unchanged, at `C=3.0041`. Every combination
either failed to converge or converged onto (or very near) the orbit's own trivial
section point — the ghost-guard ceiling. Direct diagnosis: at `max_time_factor=3.0`
(`max_time ≈ 76` nondim, ~3 periods), any combination with `k_u` or `k_s >= 4` fails
to reach the section AT ALL within the horizon (`residual=inf`, the leg never gets
there) — the ~19x weaker instability at this energy means the manifold needs
several MORE periods to separate from the orbit by a macroscopic distance than the
`C_flyby` family did. A direct probe confirmed this: seeding the unstable manifold
at `tau=0` and searching for the k-th qualifying crossing with a generous
`max_time = 15 * period` reached crossings up to `k=8` in well under a second each
— the crossings exist, they are just further out in TIME than `#754`'s own
3-period horizon allows.

## 3. The continuation experiment: `#754`'s specific branch folds before reaching C=3.0041

Before re-scanning, this task tried the more principled approach first: continue
`#754`'s OWN converged `(tau_u, tau_s)` solution (`branch_u=+1, branch_s=-1,
k_u=3, k_s=3`, `tau_u=14.463070474804017, tau_s=11.20947242577664`,
residual `2.58e-9`) in Jacobi constant, using `#761`'s own family-continuation
members as the underlying orbit at each step and re-Newton-correcting `(tau_u,
tau_s)` from the previous step's converged values (a direct, un-scanned Newton
re-solve from the previous phase converges in well under 1 second per step, versus
tens of seconds for a blind coarse re-scan — confirmed: re-solving the exact
`C_flyby` point from its own known `(tau_u, tau_s)` seed with no grid scan
reproduces the `#754` result to `8.4e-10` residual in 0.6s).

**This continuation succeeded smoothly from `C=2.99164` up to `C≈2.99941`**
(monotone `x0` drift, smoothly decaying `|lambda|` from `1006` down to `~190`,
adaptive step size growing back up to the full `5e-4` after the first few
cautious steps) **and then hit a wall**: step-halving down to `dC ~ 1.5e-7` still
failed to converge past `C = 2.9994100780240007`. A local residual-norm grid probe
around a nearby high-k candidate confirmed the qualitative signature of a
fold/tangency directly: gridding `(tau_u, tau_s)` in a `9x9` neighborhood around a
near-converged `(branch_u=+1, branch_s=+1, k_u=4, k_s=4)` point showed the `dx`
and `dxdot` residual components' sign-change boundaries becoming IDENTICAL over a
large sub-region (rows where both components flip sign together, rather than along
independent, transversally-crossing curves) — the residual vector has a local
MINIMUM there, not a zero: the two manifold curves are tangent, not transversal,
at that neighborhood.

**Conclusion: `#754`'s specific `(k_u=3, k_s=3)` self-connection branch is not the
right one to look for at `C=3.0041` — it disappears (folds) around `C≈2.9994`,
about half the remaining distance to Kumar's energy.** A genuinely different
`(branch, k)` combination was needed, motivating the wider coarse re-scan below.

## 4. The successful coarse re-scan

Scanned `branch_u, branch_s in {+1,-1}`, `k_u, k_s in 4..8` (the low-`k` range
1..3 was already established as ghost-adjacent from the initial re-run of `#754`'s
own scan) at `max_time_factor=8.0` (empirically sufficient for `k` up to 8 to reach
the section), light triage settings (`scan_n=6`, `max_iter=10`) first to identify
promising combinations quickly, refined with tighter settings
(`tol=1e-9`, `max_iter=40`) afterward. Three genuine, ghost-guard-passed,
Newton-converged hits were found (of the ~40 combinations tried before stopping —
this was a targeted, not exhaustive, scan; the full `4 x 5 x 5 = 100`-combination
grid for `k in 4..8` was not completely covered, since three independent positive
hits were already in hand):

| branch_u | branch_s | k_u | k_s | residual | crossing (x, xdot) | ghost_distance |
|---|---|---|---|---|---|---|
| +1 | +1 | 6 | 6 | `1.97e-10` | `(-1.42207150, -2.02e-10)` | `0.0368` |
| +1 | +1 | 5 | 6 | `3.33e-10` | `(-1.40661497, +0.06351614)` | `0.0670` |
| +1 | +1 | 6 | 5 | `6.80e-10` | `(-1.40661497, -0.06351614)` | `0.0670` |

The `(5,6)`/`(6,5)` pair are mirror images of each other under the CR3BP's
time-reversal reflection symmetry (`x, xdot` equal; `xdot`-sign of the crossing
flipped) — genuinely two intersections related by symmetry, not the same one found
twice. The `(6,6)` hit sits ON the symmetry axis (`xdot ≈ 0` to `2e-10`) —
structurally the SAME type of point as `#754`'s own Table-2 self-connection
(`Anderson & Lo`'s published Table 2 state also has `xdot ≈ 0`, `9e-8`) — chosen as
this task's PRIMARY reported result for that reason.

All ghost distances are far clear of `GHOST_GUARD_DELTA = 1e-3` (37x-67x margin) —
a real margin, smaller than `#754`'s own 145x at `C_flyby` (expected: the weaker
instability at this energy means less "unfolding" of the manifold per crossing,
so the genuine intersection sits somewhat closer to the orbit's own trivial point
in absolute terms than the strongly-unstable `C_flyby` case does), but not remotely
razor-thin.

## 5. Verification of the primary hit `(branch_u=+1, branch_s=+1, k_u=6, k_s=6)`

```
tau_u = 10.72913431175392,  tau_s = 14.582984371438714
crossing (x, xdot) = (-1.4220714951697728, -2.016580027963677e-10)
Newton residual = 1.9723783076056544e-10   (an order of magnitude tighter than
                                             #754's own 2.58e-9 at C_flyby)
ghost_distance = 0.03682304954560878   (37x the 1e-3 guard threshold)
```

* **Ghost guard**: the orbit's own qualifying section point is its IC itself,
  `(x0, xdot0) = (-1.3852484456241640, 0.0)` — the found crossing sits `0.0368`
  away, far outside `GHOST_GUARD_DELTA=1e-3`.
* **Independent Radau cross-check** (via `assemble_cycle`/`crosscheck_cycle`,
  `rtol=atol=1e-11`): DOP853 and Radau agree to `2.42e-8` — comfortably inside the
  `<=1e-6` mandate, and in fact tighter than `#754`'s own `6.79e-8` cross-check at
  `C_flyby`.
* **Forward/backward re-approach** (new module function
  `homoclinic_reapproach_check`): the found intersection state's FULL 6-vector was
  re-derived directly (via a new `_full_state_crossing` helper returning the full
  state and elapsed time, not just the `(x, xdot)` section-plane point
  `heteroclinic_cycle._section_crossing` reports), then:
  - propagated BACKWARD by exactly the unstable leg's own elapsed transit time
    (`t_u = 134.269...` nondim, ≈ 5.31 orbital periods) and compared to the
    ORIGINAL epsilon-scale unstable-manifold seed at `tau_u` —
    **`backward_distance = 5.638e-8`**;
  - propagated FORWARD by `|t_s| = 134.269...` and compared to the stable-manifold
    seed at `tau_s` — **`forward_distance = 2.252e-5`**.

  Both distances are far smaller than the O(1) trajectory scale (positions of
  order 1-2 nondim units), confirming the found state genuinely retraces BOTH the
  unstable leg backward and the stable leg forward — the actual definition of a
  homoclinic point, not a coincidental close pass. They are not machine-precision
  (unlike the Newton residual itself) because re-propagating a ~5.3-period,
  `|lambda|≈54.6`-per-period hyperbolic excursion amplifies ordinary integration
  roundoff by a factor of order `54.6^5.3 ≈ 10^9`; `5.6e-8` and `2.3e-5` after that
  amplification are excellent agreement, not a red flag. (An earlier hand-copied
  manual check using truncated printed digits for the intermediate state, rather
  than the full-precision value the production code carries, gave a much cruder
  `~7e-3` agreement — a lesson in why this check needed to be **production code**,
  re-deriving the state at full float precision, not a scratch script retyping
  console output; retained here as an honest account of the debugging path.)

## 6. Corroboration: the mirror pair also converges

The `(k_u=5, k_s=6)` and `(k_u=6, k_s=5)` combinations were independently
Newton-converged (residuals `3.33e-10` and `6.80e-10`), ghost-guard-passed
(`0.0670` margin each, 67x the threshold), and reproduce each other's `(x, xdot)`
to `<1e-6` under the expected reflection symmetry (`x` equal, `xdot` sign-flipped).
Forward/backward re-approach on these two: `backward_distance` `2.45e-7` /
`4.98e-6` (both tight), `forward_distance` `2.09e-2` / `2.19e-2` (looser than the
primary hit's `2.3e-5`, but still three orders of magnitude below the O(1)
trajectory scale — reported honestly as a looser but still meaningfully small
re-approach, not claimed as tight as the primary hit's).

Three independent, differently-indexed `(branch, k)` combinations all converging
to genuine, ghost-guard-passed points is strong evidence that C=3.0041 genuinely
supports transversal homoclinic self-intersections — not an isolated numerical
fluke at one specific index choice.

## 7. Code delivered

* `src/cyclerfinder/search/jovian_resonant_connections.py`:
  - `build_34lo_kumar_c_node()` — wraps `#761`'s `continue_34lo_to_kumar_c()` into
    a `ResonantNode` (reusing `ResonantNode.from_candidate` unchanged), mirroring
    `build_3_4_lo_node`/`build_5_6_lo_node`'s own pattern.
  - `find_homoclinic()` extended with optional `target: NDArray | None = None` and
    `rank_by_residual: bool = False` keyword arguments. Every pre-existing caller
    is unaffected (`target=None, rank_by_residual=False` reproduces the exact old
    behaviour, ranking by distance to `TABLE2_STATE`); a caller scanning an energy
    with no published state to target (this task) passes
    `rank_by_residual=True`, ranking by Newton-residual tightness instead and
    reporting `dist_to_table2=nan` on every returned candidate (never a distance
    to an irrelevant target).
  - `_full_state_crossing()` (private) — like `heteroclinic_cycle._section_crossing`
    but returns `(elapsed_time, full_6vector)` at the k-th qualifying crossing,
    needed for the re-approach check below (which needs the full state and the
    time already elapsed reaching it, not just where it crossed the section).
  - `HomoclinicReapproachResult` / `homoclinic_reapproach_check()` — the
    forward/backward re-approach self-consistency check (Sec. 5), built as
    reusable module code (not a one-off script): re-derives the intersection's
    full state from the candidate's own stored `tau_u/tau_s/k_u/k_s/branch_u/
    branch_s`, propagates backward/forward by the legs' own elapsed transit times,
    and reports the distance to each leg's original manifold seed.
* `tests/search/test_jovian_resonant_connections.py`: a new `node_kumar_c`
  module-scoped fixture and nine new tests covering: the continued node's own
  numbers (matching `#761`'s recorded values), the orbit's own section point, the
  primary hit's convergence + real (>10x) ghost margin, the independent Radau
  cross-check, the forward/backward re-approach distances, the mirror-pair
  corroboration, `find_homoclinic`'s new `rank_by_residual` behaviour, and a
  regression test confirming the extension does NOT change any pre-existing
  caller's default behaviour. None are marked `@pytest.mark.slow` (this is a
  discovery-verdict-bearing result; per this chain's own discipline it must run
  in CI, not be silently skipped).

## 8. Verification

* `uv run ruff check` / `ruff format --check` on both changed files: clean.
* `uv run mypy src tests` (canonical full invocation): clean, 825 source files.
* `uv run pytest tests/search/test_jovian_resonant_connections.py -q`: 27/27 pass
  (~95s wall).
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: run before the
  `OUTSTANDING.md` update commit (see commit message for result).

## 9. What this means for the torus's transport-utility question

This is the concrete piece `#761`'s own Sec. 3 flagged as the remaining gap: a
genuine homoclinic self-connection now exists AT the torus's own cited seed
energy, independently corroborated (Radau, ghost-guard, forward/backward
re-approach, a mirror-pair second finding) — not merely at a different point on
the same family (`#754`'s own `C_flyby` result). **This does NOT, by itself,
certify that this specific manifold intersection is USED by the catalogued
`europa-3-4-crnbp-torus-jupiter-2026` torus's own construction** — that torus's
numerical seed is still the project's own two-body-limit proxy (`C=2.9040` at
project `mu`), not Kumar's `C=3.0041` orbit directly (the same honest gap `#750`
and `#761` both already flagged, unchanged by this task). Whether/how to connect
this finding to the catalogue row is, per the dispatch note, explicitly a
SEPARATE, later decision — `catalogue.yaml` is not touched by this task.
