# `#775` — genuine continuation attempt on the Saturn-Titan resonant-chain closure

**Task:** `#775`, the third attempt at closing Vaquero 2013's own Saturn-Titan "resonant chain"
periodicity correction (Fig. 4.10) — a homoclinic self-connection of the 3:4 resonant orbit,
corrected into an exactly periodic new orbit. Per the dispatch note's own explicit instruction,
this task's job was NOT to retry cold-start Newton/shooting with yet another seed (`#773` already
tried that twice, honestly), but to attempt genuine CONTINUATION — a sequence of small,
individually-easy, warm-started steps — starting from `#767`'s own already-converged, well-verified
homoclinic self-connection.

**Verdict up front**: **honest, sharply-evidenced continued NEGATIVE.** Two genuine continuation
avenues were built and run in good faith; both fail decisively, and the failure is MORE sharply
characterized (not merely "eventually stalls" but "cannot make a single accepted step at ANY
tested scale, across nearly 10 orders of magnitude") than either of `#773`'s own findings. My own
honest technical assessment (requested explicitly by the dispatch note) is that this specific
closure is very likely **intractable with this project's current Newton/shooting-based numerical
toolkit, in any formulation** — see the final section, and the accompanying recommendation on
`#774` below.

---

## Sources read in full before starting

`docs/notes/2026-08-01-773-resonant-chain-periodicity-closure.md` (the immediately preceding
attempt — read in full, including its own permanent `t_cross`-drift branch guard, which this task
reuses unmodified), `docs/notes/2026-07-31-768-saturn-titan-resonant-chain.md`, and
`docs/notes/2026-07-31-767-saturn-titan-homoclinic-connection.md` (the CONFIRMED, genuinely
converged starting point this task continues from), plus `#753`'s and `#761`'s own continuation
methodology (`continue_family`/`continue_34lo_to_kumar_c`, `src/cyclerfinder/search/
cr3bp_continuation.py`, `src/cyclerfinder/search/jovian_resonant_families.py`) as the precedent for
"small adaptive steps, warm-started, gated by fold/topology-jump detection." All of `#767`'s own
CONFIRMED homoclinic-connection code
(`src/cyclerfinder/search/saturn_titan_resonant_connections.py`: `attempt_chain_closure`,
`attempt_chain_closure_multiple_shooting`, the `t_cross` guard, `resonant_chain_target_point`,
`rank_by_proximity_to_65`) was read and reused directly, not reimplemented.

---

## The starting point: `#767`'s own near-6:5 candidate, re-derived and independently verified

Re-ran `#768`'s own closest-to-6:5 homoclinic self-connection directly this task
(`branch_u=branch_s=-1`, `k_u=4`, `k_s=5`), reproducing its exact numbers independently:

```
converged=True, residual=9.6277e-10
crossing_xv = (x=0.91407251, xdot=-0.09173657)
dist_to_65 = 0.09404305044147829
t_u = 39.80475100319532, t_s = -70.69486131241366
t_target = t_u + |t_s| = 110.49961231560899
backward_distance = 1.2656e-09, forward_distance = 0.0053082
```

This exactly reproduces `#768`'s own reported `(t_u=39.805, t_s=-70.695, t_target=110.4996,
backward=1.27e-9, forward=0.0053)` — confirming this is a stable, reproducible, genuine,
well-converged homoclinic self-connection, the correct starting point per the dispatch note's own
instruction ("pick whichever of the 4 is closest to 6:5's own fixed point" — this is it, `dist=0.094`,
the closest of `#767`'s own four hits).

Its OWN periodicity-map residual (the quantity `attempt_chain_closure`'s corrector tries to drive
to zero, evaluated at the fixed `crossing_index=12` nearest `t_target`, `n_events_seed=13`,
`t_cross=106.3328`) is **`R0 = (-2.16251, -0.06076)`, `‖R0‖ = 2.1634`** — large, confirming (as
`#773` already found) that a genuine homoclinic self-connection is NOT automatically close to a
solution of the DIFFERENT periodicity-fixed-point problem.

---

## Avenue 1: artificial homotopy in the periodicity-map's own residual target

### Rationale

`#773`'s own final recommendation was continuation from an already-converged solution. The only
genuinely already-converged, well-verified nearby solution available is `#767`'s own homoclinic
connection above — but it solves a DIFFERENT nonlinear system (manifold matching, `Wu(tau_u) =
Ws(tau_s)`) than the periodicity-map fixed point `attempt_chain_closure` needs (`x_ret(x,xdot) = x`,
`xdot_ret(x,xdot) = xdot` at the fixed crossing index). This task built a **Newton/artificial-parameter
homotopy** bridging the two: define

```
target_residual(s) = (1 - s) * R0,   s in [0, 1]
```

At `s=0` the seed trivially satisfies "residual == target_residual(0) == R0" (no correction
needed — this is exactly where the walk starts, at the already-converged homoclinic point). At
`s=1`, `target_residual = 0`, the TRUE periodicity condition. Walking `s` from 0 to 1 in small,
adaptively-sized steps, each step warm-started from the PREVIOUS step's own converged `(x, xdot)`
and Newton-corrected onto the new intermediate target (reusing `_chain_map_step`'s own STM-based
Jacobian and the SAME `#773` `max_t_cross_drift` branch-drift guard), is a standard, legitimate
numerical-continuation technique (an "artificial parameter"/Newton-homotopy method) for exactly
this situation: a good starting point exists, but the target is far away in a badly-conditioned
direction. This is a genuinely different algorithm from `#773`'s own single-shot damped Newton (see
the module docstring's own detailed rationale) — it takes a SEQUENCE of small, individually
Newton-converged steps, not one large damped step toward the full target.

Shipped as `continue_chain_closure_homotopy` in
`src/cyclerfinder/search/saturn_titan_resonant_connections.py`, with its own
`ChainHomotopyStep`/`ChainHomotopyResult`/`ChainHomotopyStopReason` dataclasses, reusing
`_chain_map_step` directly (not reimplemented) and the `#773` branch-drift guard at every
micro-step.

### Direct diagnostic: the map has NO usable step-size window near this seed

Before running the full walk, a direct single-step diagnostic (varying only the target-shrinkage
parameter `s`, taking exactly ONE full Newton step from the seed toward `target_residual(s)`, no
backtracking) shows the map's local linear model breaks down at every practically achievable scale:

| `s_try` | step norm `‖Δ(x,xdot)‖` | `t_cross` (seed: `106.3328`) | actual `‖mod_residual‖` | linear-predicted (`s·‖R0‖`) | ratio |
|---|---|---|---|---|---|
| `0.05` | `3.47e-2` | — | **map step failed** (bad radicand / index) | `0.108` | — |
| `0.01` | — | — | **map step failed** | `0.0216` | — |
| `0.001` | `6.93e-4` | `113.16` | `0.4066` | `0.00216` | **188x** |
| `1e-4` | `6.93e-5` | `102.62` | `0.2825` | `2.16e-4` | **1308x** |
| `1e-5` | `6.93e-6` | `106.41` | `1.523e-2` | `2.16e-5` | **705x** |
| `1e-6` | `6.93e-7` | `106.334` | `1.744e-4` | `2.16e-6` | **81x** |
| `1e-7` | `6.93e-8` | `106.3328` | `9.545e-7` | `2.16e-7` | **4.4x** |
| `1e-8` | `6.93e-9` | `106.3328` | `5.198e-7` | `2.16e-8` | **24x** (worse, not better) |

This directly, quantitatively confirms and SHARPENS `#773`'s own Finding 2 (map sensitive at the
"8th significant digit"): there is no scale at which a single linearized Newton step tracks the
actual nonlinear residual — steps large enough to matter (`>~1e-5`) are already 80x-1300x off the
local linear prediction (well outside ANY plausible trust region), while steps small enough to
approach the linear regime (`<~1e-7`) are already comparable to or below this problem's own
`rtol=atol=1e-13/1e-14` integration-noise floor (note the `s=1e-8` row is WORSE than `s=1e-7`, not
better — a clean signature of numerical noise dominating, not genuine curvature). **There is no
window between "big enough to matter" and "small enough to be linear" for this seed's own
periodicity map** — a materially sharper diagnosis than `#773`'s own framing.

### Full homotopy walk: literally zero steps accepted, across 10 orders of magnitude

Ran the shipped `continue_chain_closure_homotopy` from this exact seed (`t_target=110.4996`,
default `ds_init=0.05`, `ds_min=1e-9`, `max_inner_iter=8`, `max_backtrack=8`, the default
`0.5 * node.period` branch-drift cap):

```
converged=False
stop_reason=STEP_FLOOR
s_reached=0.0
steps=[seed only]        (26 outer-step attempts, ds shrunk 0.05 -> ~7.5e-10, 290 total map calls)
wall time ~92s
```

**Not one single step — at any tested `ds`, from `0.05` down to `7.5e-10`, spanning nearly 10
orders of magnitude — was ever accepted.** This is a materially sharper negative than `#773`'s own
reported "genuine on-branch stall, residual `~0.16`-`1.9`" (which DID accept several early steps
before stalling) or the multiple-shooting "decelerating crawl" (which made continuous, if
decelerating, progress over hundreds of iterations). Here, the very FIRST micro-step already fails,
and shrinking the step 9 further orders of magnitude never rescues it.

### Positive control and sanity check (confirms this is a real finding, not a bug)

Two checks confirm the homotopy MACHINERY itself is correct, so the zero-progress result above is a
genuine finding about this specific problem, not an implementation bug:

1. **Positive control** — seeding at `node`'s own IC with `t_target = node.period` (`node` is
   itself exactly periodic, so `R0 ~ 4.1e-10` already): the walk sails through `s=0..1` and reports
   `converged=True, stop_reason=REACHED_TARGET, s_reached=1.0` (21 accepted steps).
2. **Sanity check on a genuinely easier, deliberately non-physical toy case** — a seed perturbed
   `1e-3` off `node`'s own IC, targeting the SAME short horizon (`node.period`, i.e. far less
   compounded growth than the real `~4.2`-period chain target): the walk DOES make real, if slow,
   incremental progress (`s` reached `~0.106` after 8 outer steps, each step genuinely accepted,
   not degenerate) — confirming the harness correctly makes progress when the underlying
   conditioning allows it, and that the near-6:5 seed's own zero-progress result is a genuine
   property of THAT specific `~4.2`-period, `|lambda|~2129.8` compounded-instability regime, not an
   artifact of the code.

---

## Avenue 2: is there an easier, less-unstable nearby Jacobi constant to bootstrap from?

`#773`'s own honest technique assessment suggested, as an alternative: "start from an
ALREADY-converged, much less unstable nearby periodic solution (if one exists at a lower Jacobi
constant...) and continue in whatever parameter increases the instability up to this regime." This
requires first checking whether such an easier regime actually EXISTS nearby, before investing in
re-deriving a homoclinic connection + chain closure attempt there. Reused `#753`'s own
`cr3bp_continuation.continue_family` directly (no new code needed — fully system-agnostic, per
`#764`'s own confirmed finding) to continue node's OWN 3:4 family in Jacobi constant `C`, in both
directions from `C=3.010000`:

| Direction | Result | `\|lambda\|` trend |
|---|---|---|
| **increasing `C`** | `stop_reason=TOPOLOGY_JUMP` after 5 members, last clean member at **`C=3.014000`** (`x0=1.037655`, `period=27.0824`, `lambda=1615.105`) | `2129.8 -> 1615.1` (dropping) as `C` rises `3.010 -> 3.014`, then the branch's own topology jumps |
| **decreasing `C`** | `stop_reason=MAX_STEPS`, 21 members reaching `C=2.990` cleanly | `2129.8` drops to a shallow MINIMUM `~1300.2` near `C≈2.992`, then rises again to `1311.3` by `C=2.990` |

**Finding: no substantially easier regime is reachable within this family's own
topologically-connected range.** The best available reduction is `|lambda|` `2129.8 -> ~1300`
(≈39% smaller, at `C≈2.992`) — NOT an order of magnitude. Re-estimating the compounded growth over
a comparable `~4.2`-period loop at that easier point (`exp(ln(1300)/25.72 * 4.2*25.72) ≈ 1.2e13`)
is only ~10x more tractable than the `1.2e14` at `C=3.01` — nowhere near enough to change the
qualitative outcome (the homotopy diagnostic above shows the usable step window vanishes across
~10 orders of magnitude at `C=3.01`; a mere 10x improvement in compounded growth would not open a
usable window). There is no nearby "easy" Jacobi constant to bootstrap a converged chain from
within this family's own connected branch.

**A suggestive, NOT proven, structural coincidence worth flagging explicitly**: node's OWN 3:4
family's topology-jump point in the increasing-`C` direction lands EXACTLY at `C=3.014000` — the
same value Vaquero's own thesis names as the chain family's suspected termination bound
(`"...ends for a value of Jacobi constant C < 3.01400"`, Fig. 4.12). This task did NOT continue the
actual chain family itself (never having a converged member to continue from), so this is only a
coincidence in node's OWN, related-but-distinct family — offered as a genuinely interesting lead
for whoever next investigates `#774`'s own termination-claim question, not as evidence resolving
it. It IS consistent with the hypothesis that Vaquero's own chain family is continuable only in ONE
direction from its own `C=3.01` construction point (upward, toward `C=3.014`, where some genuine
structural change occurs), which would also explain why `C=3.01` — nominally "below" the `3.014`
threshold — is nonetheless where the constituent 3:4/6:5 orbits and the homoclinic-selection
construction are anchored (Table 4.1's own stated value): it may sit very close to where the
family's own branch begins, not deep inside a broad existence interval.

---

## Honest technique assessment (per the dispatch note's own explicit request)

Both avenues were tried in good faith, reusing this project's own existing corrector/continuation
machinery throughout (never reimplemented), and both are clean, well-instrumented, quantified
negatives — sharper and more decisive than `#773`'s own findings, not merely a repeat of them:

1. The artificial-homotopy continuation — the closest analogue to "continue from an
   already-converged solution" that this problem's own structure actually supports — cannot accept
   a single step at ANY tested scale from `0.05` down to `7.5e-10`. The direct single-step
   diagnostic shows why: there is no step size at which this specific seed's own periodicity map is
   simultaneously (a) large enough to matter and (b) well-approximated by its own local Jacobian.
   Steps above `~1e-5` already produce residual changes `80x`-`1300x` worse than the linear
   prediction; steps below `~1e-7` are swamped by this problem's own double-precision integration
   noise floor. This is a genuinely different, sharper diagnosis than "eventually stalls after real
   progress" — it is "cannot begin at all, at any achievable resolution."
2. The Jacobi-constant sensitivity survey shows no substantially easier nearby regime exists to
   bootstrap from within the family's own topologically-connected range (`|lambda|` only varies
   `2129.8` down to `~1300`, an order of magnitude short of what would matter).

**My own honest final assessment**: this specific resonant-chain closure — correcting a homoclinic
self-connection at this particular `~4.2`-period / `|lambda|~2129.8` compounded-instability
combination into an exactly periodic orbit — appears **genuinely intractable with this project's
current Newton/shooting-based numerical toolkit**, in every formulation tried across `#773` and
`#775` combined: single-shooting (stalls), multiple-shooting (decelerates, never accelerates),
and now artificial-homotopy continuation (cannot take a single step at any scale). The common
thread across all three failures is the SAME underlying fact, now quantified three independent
ways: this system's compounded exponential growth over the ~4.2-period loop
(`exp(0.293*110.5)≈1.2e14`) leaves no numerically-resolvable window between "a step big enough to
matter" and "a step small enough to trust the local linear model" — the map's own fractal structure
at this scale is simply below what double-precision arithmetic (and this task's own best-available
starting point) can resolve.

This does NOT mean the chain orbit does not exist — Vaquero's own thesis demonstrates it does (Fig.
4.10, continued in Fig. 4.11). It means finding it with THIS toolkit's methods, from THIS starting
point, is not achievable within any reasonable further effort along these lines. Per `#773`'s own
flagged-but-unexplored escape hatch: Vaquero cites Lo & Parker's own "iterative multi-patchpoint
refinement" methodology (ref. [68] in her own thesis) for this exact class of problem — closer in
spirit to `#773`'s own multiple-shooting attempt (which came closest of the three techniques tried,
showing real if decelerating progress) than to either single-shooting formulation, but apparently
requiring a materially different refinement strategy (adaptive patch-point placement/insertion,
not merely more uniform segments — `#773`'s own `n_segments=8->16` test found finer UNIFORM
segmentation alone did not qualitatively help) that neither `#773` nor this task attempted. This
would be a genuinely separate, substantial undertaking (a new corrector architecture, not a
parameter tweak to what exists), not a natural next increment of the work done across
`#767`-`#775`.

---

## Recommendation on `#774` (a technical judgment for the coordinating session to weigh, not a decision made here)

`#774` (the continuation-in-`C` campaign to confirm/refute Vaquero's own `C < 3.01400` termination
claim) is gated on a converged chain orbit at `C=3.01` to continue from. Across `#773` and this
task combined, THREE independent, well-instrumented technique classes (single-shooting,
multiple-shooting, artificial-homotopy continuation) have all failed to produce one, each failing
in a genuinely different, decisive way, with the underlying obstruction (no resolvable
numerical step-size window at this compounded-instability scale) now characterized quantitatively
rather than merely observed as a stall. My own honest technical judgment is that `#774` should be
considered **closed-out/abandoned as currently scoped**, rather than merely "blocked" — not because
the underlying question (does the chain family really terminate near `C=3.014`?) is uninteresting
(the structural coincidence found in Avenue 2 above is, if anything, a reason to think it's a real
phenomenon), but because the gating prerequisite (a converged `C=3.01` chain orbit from this
toolkit) is very unlikely to materialize through further effort along the same lines. If this
question is revisited, it should be scoped as a genuinely new undertaking — either (a) building
Lo & Parker's own multi-patchpoint iterative refinement technique from scratch (a new corrector
architecture, not a further attempt with existing tools), or (b) testing the Avenue-2 structural
coincidence directly by attempting to continue node's OWN (non-chain) 3:4 family's topology-jump
point at `C=3.014` further, to see whether IT is itself the signature of whatever Vaquero's own
Fig. 4.12 curve is describing — a smaller, cheaper, and independently informative investigation
that does not require ever closing the chain orbit itself. Both are explicitly NOT attempted here,
registered only as candidate directions for a future, separately-scoped task, per this task's own
dispatch note's own escape valve.

---

## Code delivered

* `src/cyclerfinder/search/saturn_titan_resonant_connections.py` (extended, not replaced):
  `ChainHomotopyStopReason`, `ChainHomotopyStep`, `ChainHomotopyResult`, and
  `continue_chain_closure_homotopy` — the artificial-parameter homotopy continuation (Avenue 1
  above), reusing `_chain_map_step`'s own STM-based Jacobian and the `#773` `max_t_cross_drift`
  branch-drift guard at every micro-step. No changes to `attempt_chain_closure`,
  `attempt_chain_closure_multiple_shooting`, or any other existing function — purely additive.
  Avenue 2 (the Jacobi-constant sensitivity survey) needed no new code: it directly reuses
  `cr3bp_continuation.continue_family` and `saturn_titan_resonant_families.recover_table41_candidate`
  as-is, confirming `#764`'s own system-agnostic finding once again.
* `tests/search/test_saturn_titan_resonant_connections.py` (extended, 35 tests total, up from 33):
  `test_continue_chain_closure_homotopy_positive_control_trivial_seed` (confirms the machinery is
  correct: a trivially-already-periodic seed converges immediately) and
  `test_continue_chain_closure_homotopy_near65_seed_makes_zero_progress` (the headline finding,
  regression-guarded with a bounded `ds_min=1e-3` to keep CI runtime reasonable — asserts
  `stop_reason=STEP_FLOOR`, `s_reached=0.0`, exactly one step (the seed) ever accepted). Neither
  marked `@pytest.mark.slow`.

## Verification

* `uv run ruff check` / `ruff format --check` on both changed files: clean.
* `uv run mypy --strict` on both changed files: clean.
* `uv run mypy src tests` (canonical full invocation): clean, 827 source files.
* `uv run pytest tests/search/test_saturn_titan_resonant_connections.py -q`: 35/35 pass.
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: both pass (run again before the
  `OUTSTANDING.md` commit — see commit history).

## Net effect on `#775`

**DONE — honest, sharply-evidenced continued NEGATIVE.** Genuine continuation (an artificial
residual-homotopy, warm-started, gated, reusing `#773`'s own branch-drift guard) was built and run
from `#767`'s own already-converged near-6:5 homoclinic candidate, per this task's own explicit
mandate — and it fails more decisively than either of `#773`'s own cold-start attempts: zero
accepted steps at any tested scale across ~10 orders of magnitude, with a direct diagnostic
explaining exactly why (no resolvable step-size window between "matters" and "linear" for this
seed's own periodicity map). A secondary Jacobi-constant survey found no easier nearby regime to
bootstrap from instead. My own honest assessment: this specific closure is very likely intractable
with this project's current Newton/shooting toolkit; `#774` should likely be considered
closed-out/abandoned as currently scoped rather than merely blocked (recommendation only, for the
coordinating session to weigh) — see the dedicated section above for the full reasoning and two
candidate directions if this thread is ever revisited.
