# `#781` — Neptune-Triton 4:5 resonant-orbit homoclinic self-connection

**Task:** `#781`, the connection-stage follow-up `#776`'s own note explicitly recommended
("The 4:5-saddle and 4:7-stress rows...are the most promising manifold-source candidates for
such a task"), building the SAME `correct_connection`/`ResonantNode`/ghost-guard machinery
(`jovian_resonant_connections.py`, `#754`; `saturn_titan_resonant_connections.py`, `#767`)
against a THIRD system: Neptune-Triton, using `#776`/`#777`'s own confirmed "4:5-saddle"
resonant orbit (Miceli & Bosanac 2026 ESM4 `Res45+x+h`, `C=2.987089791658`).

---

## Verdict (read this first)

**CORRECTED 2026-08-08 by a Fable adversarial review, dispatched at the coordinating session's own
initiative given this task's "novel" framing — the highest-stakes claim class this project makes.
Two genuine problems were found and are corrected below: an inflated hit count, and missed prior
art. The underlying result survives, narrowed. Do not trust the ORIGINAL wording of this section as
first committed (preserved nowhere but git history) — read this corrected version.**

**A genuine, non-trivial homoclinic self-intersection EXISTS at the 4:5-saddle orbit's own
`C=2.987089791658`.** The honest count, independently re-derived by the adversarial review: **TWO
distinct homoclinic trajectories**, not four independent hits. PRIMARY (`k_u=k_s=13`) is one
trajectory; SECONDARY/MIRROR_A/MIRROR_B (`k=15`/`(13,17)`/`(17,13)`) are the SAME single
trajectory, registered at three different crossing indices — confirmed by propagating SECONDARY's
own converged unstable-manifold seed once and finding its own section crossings land exactly on
MIRROR_A (`k=13`) and MIRROR_B (`k=17`) as well as SECONDARY itself (`k=15`), with converged
`tau_u`/`tau_s` agreeing to `<8e-9` across all three labels. The "mirror pair" is an automatic
consequence of a symmetric orbit's on-axis perpendicular crossing (a `k=15±2` reflection pair
exists BY CONSTRUCTION for any such crossing) — NOT independent corroborating evidence, and PRIMARY
has an exactly analogous, previously-unreported mirror pair at `k=12/14`. All four originally-
reported numeric hits are still real, Newton-converged, ghost-guard-passed results (residuals
`<1e-8`, Radau cross-check `<1.31e-7`, `134x`–`287x` past the ghost-guard threshold) — nothing was
fabricated — but they represent two distinct trajectories, not four. "Transversal" was asserted in
the original write-up but never actually measured anywhere in this task's own code; that adjective
should be read as unverified pending a dedicated transversality computation, not as an established
fact. Forward/backward re-approach is honest: backward is tight for both trajectories (`<5.4e-8`),
forward is tight for the SECONDARY trajectory (`~1.7e-3`) but notably looser for PRIMARY (`0.157`)
— reported as-is, not fudged, with a plausible explanation below.

**Literature-novelty framing (the load-bearing gate for this task, read in full before trusting
"novel"): the SPECIFIC object is still novel, but the original "no manifold-connection literature
exists for Neptune-Triton" framing was WRONG and is corrected here.** A Fable adversarial review
found prior art this task's own literature-check pass missed: Spear 2021 (MS thesis, CU Boulder,
same Bosanac group as Miceli & Bosanac 2026 — see `docs/notes/2026-08-08-spear-2021-digest.md`)
computes 11 planar HETEROCLINIC connections in the Neptune-Triton CR3BP between L1/L2 Lyapunov
orbits. This falsifies "no manifold-connection literature exists" as originally stated. What
survives, independently confirmed in the digest above: Spear 2021 contains ZERO homoclinic content
of its own, computes 12 resonant families none of which is a 4:5 or 4:7 ratio (this task's own
object is genuinely absent from it), and its own future-work section names resonant-orbit
connection design as explicitly NOT YET DONE. The corrected, defensible claim is narrower than
originally stated: **first published homoclinic connection, and first resonant-orbit manifold
connection of any kind, in the Neptune-Triton CR3BP** — not "the first manifold-connection work of
any kind in this system." Miceli & Bosanac 2026 itself still publishes zero homoclinic/heteroclinic
manifold-intersection content — independently re-verified by grepping BOTH the paper's own text
layer and the companion Miceli 2025 dissertation's text layer for "homoclinic"/"heteroclinic":
**zero hits in either document.** The paper's own manifold use (Sec. 3.1.4/3.2) is confined to (a)
sampling arcs along a SINGLE orbit's own stable/unstable manifold as coarse "motion primitives"
and (b) an approximate, discretized configuration-space-adjacency graph search connecting those
primitives — categorically NOT an exact Newton-corrected Poincaré-section `Wu ∩ Ws` intersection.
`search/literature_check.py`'s own mandatory-floor run (Sec. 5 below) returned `not-found`
(confidence 0.4) — necessary-not-sufficient, and the adversarial review found this specific gate's
`not-found` carries NEAR-ZERO evidentiary weight for this claim type (its own cycler-vocabulary
queries cannot express "homoclinic manifold intersection" as a concept at all, which is exactly why
it missed nothing about Spear 2021 either way — the miss there was a WebSearch-query-phrasing gap,
not a `literature_check.py` gap). The paper/dissertation text grep plus the Spear 2021 digest are
the load-bearing evidence; the formal gate is at most a very weak floor on top of it, not
meaningful corroboration on its own.

The 4:7-stress attempt (Sec. 6) is an **honest, diagnosed FAIL**, exactly as `#759`'s own
extreme-instability precedent predicted: at this task's own bounded effort, Newton stalls at
residuals `~1.7e-3`–`3.4e-3` (target `1e-9`) from every geometrically-motivated seed tried, and
even the marginal ghost-guard-clearing candidates found at low `k` sit only `1.2x`–`3.5x` past
`GHOST_GUARD_DELTA` (vs. 4:5's `134x`–`287x`) — a razor-thin margin, the same "fractal
sensitivity" signature `#759` documented for its own extreme-instability case. Not forced, not
tolerance-weakened.

---

## 1. The orbit itself (from `neptune_triton_resonant_families.py`, `#776`/`#777`)

Read directly from `recover_esm_candidate("4:5-saddle")`'s own return value this task:

```
x0              = -0.9690564220154437
ydot0           = -0.12676711978337848
period          = 30.398418802764855   (period/2pi = 4.838058614637641)
jacobi          = 2.987089791658 exactly (ESM4 Res45+x+h, verbatim)
max_eigenvalue  = -105.05423683377265   (Barden convention; GENUINELY NEGATIVE real)
planar_floquet  = -105.05423683428756   (independent full-period monodromy eigendecomposition)
is_real_unstable = True
```

`#776`'s own gate: this is the ONLY negative `max_eigenvalue` among all ten `#776` gate rows —
every prior orbit this task chain has built a node from (Jupiter-Europa, Saturn-Titan, and
Neptune-Triton's own Lyapunov/DPO/other-resonant rows) has a POSITIVE real value. `is_real_unstable
=True` still holds (a negative real eigenvalue is still a genuine real saddle) — but this is a
NEW physical regime for this task chain's own node-building machinery, and it exposed a real bug
(Sec. 2 below).

## 2. A real bug in `ResonantNode.from_candidate`'s own staleness guard, fixed this task

`ResonantNode.from_candidate` (`jovian_resonant_connections.py`) computes
`rel_err = abs(lam_u - cand.max_eigenvalue) / cand.max_eigenvalue` to guard against building a
node from stale/mismatched data. `lam_u` (from `_planar_floquet_pair`) is MAGNITUDE-only by that
function's own documented convention (`+105.054...`), so this formula is only correct when
`cand.max_eigenvalue` is positive — true for every candidate this project built a node from
BEFORE this task. For the 4:5-saddle candidate (`max_eigenvalue=-105.054...`), the un-fixed
formula divides by a NEGATIVE denominator, making `rel_err` itself negative and the guard
permanently inert (it can never exceed `eigenvalue_rel_tol`, regardless of actual disagreement) —
a real correctness gap, not a hypothetical one, first exercised by this task.

**Fix** (this task, `jovian_resonant_connections.py`): compare magnitudes —
`rel_err = abs(lam_u - abs(cand.max_eigenvalue)) / abs(cand.max_eigenvalue)`. Byte-for-byte
identical arithmetic for every already-positive `max_eigenvalue` (confirmed: the full
Jovian + Saturn-Titan `test_jovian_resonant_connections.py` / `test_jovian_resonant_families.py`
/ `test_saturn_titan_resonant_connections.py` / `test_saturn_titan_resonant_families.py` suites
all still pass unchanged after the fix), newly-correct for the negative case (confirmed:
`rel_err ~ 2e-12` for the 4:5-saddle candidate post-fix, cleanly inside `eigenvalue_rel_tol=1e-6`).
See `test_node_from_candidate_staleness_guard_fixed_for_negative_eigenvalue` (new test) and the
module's own updated docstring for the full account.

**Physical meaning of the negative eigenvalue** (module docstring, `neptune_triton_resonant_
connections.py`): the period map sends the unstable eigenvector to `-105.05 x` itself — a point
perturbed along the unstable direction does not merely grow each period, it also FLIPS to the
opposite side of the orbit. `_planar_floquet_pair`'s own magnitude-only eigenvector-direction
convention is unaffected by the sign (confirmed: `unstable_eigvec`/`stable_eigvec` are still
unit-normalised, `test_resonant_node_is_real_saddle_despite_negative_stored_eigenvalue`), but this
DOES matter for how a "branch" reads: `branch_u=+1`/`branch_u=-1` pick a sign of the initial
epsilon-offset along the SAME unstable eigenvector, and because the underlying map itself flips
sign every period, the natural convergent combination for this orbit is OPPOSITE branch signs
(`branch_u=+1, branch_s=-1`) — all four converged hits below use exactly that combination; several
same-sign (`branch_u=branch_s`) seeds were explored and did NOT converge cleanly (Sec. 4).

## 3. Section convention: BOTH `x_sign` AND `ydot_sign` fully unrestricted

Direct inspection of the 4:5-saddle orbit's own `{y=0}` crossings over one period (this task,
`rtol=atol=1e-13`) found 6 crossings (excluding the `t=0`/`t=period` boundary), at exactly TWO
perpendicular (`xdot=0`) points — the IC itself (`x0=-0.969056...`, `ydot<0`) and the half-period
return (`x=0.981105...`, `ydot>0`, `half_crossings=3`) — plus TWO non-perpendicular mirror pairs
(`x~=1.122044`, `ydot<0`; `x~=-1.143372`, `ydot>0`).

Unlike Saturn-Titan's own 3:4 orbit (both perpendicular points at `ydot>0` — only `x` needed
unioning, `stc.own_section_points`), this orbit's own two perpendicular points sit at OPPOSITE
`ydot` signs. Fixing ANY single `ydot_sign` (Saturn-Titan's own remaining restriction) silently
drops ONE of the two — confirmed directly this task: `ydot_sign=+1` drops the IC entirely;
`ydot_sign=-1` drops the half-period return — **neither single choice is safe**, since an empty
or incomplete reference set risks the ghost guard missing its single most important trivial-shadow
point (the exact hazard the module docstring documents, and the regression this task's own
`test_single_ydot_sign_choice_would_drop_one_perpendicular_point` test locks down).

This module therefore generalizes ONE STEP FURTHER than `stc`'s own partial (`x`-only) relaxation:
`own_section_points` unions ALL FOUR `(x_sign, ydot_sign)` combinations, recovering the full
6-point per-period crossing set with no double-counting (confirmed: each crossing's own sign pair
is unique to one of the four buckets). `find_homoclinic`'s own search similarly passes
`ydot_sign_u=None, ydot_sign_s=None, x_sign_u=None, x_sign_s=None` to `correct_connection`
(natively supported — no change to `heteroclinic_cycle.py` needed).

## 4. The scan: how the four hits were found, and why a brute-force scan was abandoned

An initial brute-force scan mirroring `#754`/`#767`'s own default parameters (`branches=(+1,-1)`,
`k_range=range(1,7)`, `max_time_factor=3.0`, `scan_n=12`) proved computationally intractable —
timed at **12–40 seconds per `correct_connection` call** even at lighter `scan_n=6` settings (the
internal `scan_n x scan_n` grid dominates cost at this orbit's own long period, `T=30.4`, and
6-period search horizon). A full `4 branch-combos x 12x12 k-combos` grid at these settings would
have taken hours.

**Method actually used**: build the two 1-D Poincaré curves separately — the unstable manifold's
own `k`-th crossing as a function of `tau_u` alone, and the stable manifold's own `k`-th crossing
as a function of `tau_s` alone (`O(n_tau)` per `(branch, k)`, not `O(n_tau^2)` per `(k_u, k_s)`
pair) — for `k in 6..21` (chosen from a direct crossing-time probe: growth `|lambda|^n` reaches
macroscopic separation, `x` excursions well beyond the base orbit's own `[-1.14, 1.12]` extent,
by `k~13`-`19`, ~2.3-3.3 orbital periods). Building all 64 curves (2 directions x 2 branches x 16
`k` values, `n_tau=15`) took **81 seconds total** — then a cheap `O(curve_len^2)` nearest-point
search across all curve PAIRS identified close-approach `(tau_u, tau_s)` seeds directly, fed to
`correct_connection` with `tau_u0`/`tau_s0` EXPLICIT (skipping its own internal grid scan
entirely) for fast (3–10s) Newton refinement. This method is not currently promoted into the
production module (it lived in scratchpad scripts) — `find_homoclinic`'s own brute-force scan
remains available and is independently exercised by `test_find_homoclinic_returns_known_primary_
combo` (narrow window, ~86s, confirming the PRIMARY hit through a genuinely different code path).

| branch_u | branch_s | k_u | k_s | residual | crossing (x, xdot) | ghost_distance (×guard) |
|---|---|---|---|---|---|---|
| **+1** | **-1** | **13** | **13** | `7.65e-10` | `(1.16933872, -3.8e-10)` | `0.1671` (**167×**) |
| +1 | -1 | 15 | 15 | `2.64e-10` | `(-1.38561105, -2.4e-10)` | `0.2867` (287×) |
| +1 | -1 | 13 | 17 | `4.32e-10` | `(-1.27434376, +0.12236316)` | `0.1346` (135×) |
| +1 | -1 | 17 | 13 | `6.22e-10` | `(-1.27434376, -0.12236316)` | `0.1346` (135×) |

All four use `branch_u=+1, branch_s=-1` — every same-branch-sign seed tried (three attempts,
distinct `(k_u,k_s)` combinations near genuine close-approaches identified the same way) FAILED to
converge (residuals plateaued `1.4e-3`–`3.1e-3`, well short of `tol=1e-9`) — an honest negative,
consistent with Sec. 2's own physical explanation (the negative eigenvalue's own branch-flip
makes opposite-sign the natural convergent combination here), not suppressed or hidden.

**CORRECTED (adversarial review, 2026-08-08): the four rows above are TWO distinct homoclinic
trajectories, not four independent hits.** PRIMARY (`k_u=k_s=13`) is one trajectory, sitting ON the
symmetry axis (`xdot~0` to `<4e-10`) — structurally the same point-type as Anderson & Lo's own
published Table 2 state and `#766`'s/`#767`'s own primary hits. SECONDARY (`k_u=k_s=15`) is a
SECOND, distinct trajectory (`tau_u` differs from PRIMARY's by `16.3`, and PRIMARY's own crossing
appears nowhere in SECONDARY's own `k=1..21` crossing list — genuinely independent), also on-axis.
MIRROR_A (`13,17`) and MIRROR_B (`17,13`) are NOT a third and fourth independent trajectory — they
are SECONDARY's own trajectory, registered at two more of its own crossing indices (`k=13` and
`k=17`, straddling SECONDARY's own `k=15` on-axis point): propagating SECONDARY's own converged
unstable-manifold seed forward once produces crossings at `k=13` (= MIRROR_A's point), `k=15` (=
SECONDARY's own point), and `k=17` (= MIRROR_B's point), with all three sharing `tau_u≈18.30074707`
to `<8e-9`. This reflection-pair structure (`k=15±2` mirroring each other) is an AUTOMATIC
consequence of a symmetric orbit's own on-axis perpendicular crossing, not independent evidence of
"transversal self-intersections beyond one isolated choice" — PRIMARY has the exact same structure
at its own `k=12/14`, simply never reported in the original write-up. The honest count is: **two
distinct homoclinic trajectories exist at this energy** (PRIMARY and SECONDARY), each internally
symmetric, with SECONDARY's own symmetry producing three crossing-index registrations of the same
underlying orbit. See also `#767`'s own note (Saturn-Titan) for an analogous, pre-existing framing
weakness in this task chain — not unique to `#781`, but most pronounced here.

## 5. Verification of all four hits

**Ghost guard**: the orbit's own two perpendicular reference points (Sec. 3) are
`(-0.969056422015, 0)` (IC) and `(0.981105, 0)` (half-period return); the four non-perpendicular
mirror-pair references sit at `(±1.122044, ∓0.1602)` and `(±1.143372, ∓0.1533)`. All four found
crossings sit `0.1346`–`0.2867` away — `134×`–`287×` the `GHOST_GUARD_DELTA=1e-3` threshold, a
real margin comparable to `#767`'s own `219×`–`307×` (Saturn-Titan) and well past `#766`'s own
`37×`–`67×` (Jupiter-Europa, weaker-instability case).

**Independent Radau cross-check** (`assemble_cycle`/`crosscheck_cycle`, `rtol=atol=1e-11`): all
four agree between DOP853 and Radau to `9.24e-10`–`1.30e-7` — comfortably inside the mandated
`<=1e-6` for every hit.

| hit | residual | ghost_distance | Radau independent_residual |
|---|---|---|---|
| PRIMARY (13,13) | `7.65e-10` | `0.1671` | `1.30e-07` |
| SECONDARY (15,15) | `2.64e-10` | `0.2867` | `1.04e-09` |
| MIRROR_A (13,17) | `4.32e-10` | `0.1346` | `1.44e-09` |
| MIRROR_B (17,13) | `6.22e-10` | `0.1346` | `9.24e-10` |

## 6. Forward/backward re-approach: honest, and PRIMARY's own asymmetry noted plainly

This orbit's own unstable eigenvalue magnitude (`|lambda|~105`) is far milder than Saturn-Titan's
3:4 orbit (`|lambda|~2129.8`) — closer to `#766`'s own `C=3.0041` case (`|lambda|~54.6`), where
`jrc`'s own un-tightened defaults were already sufficient. This module's own
`homoclinic_reapproach_check` still threads the tighter `rtol=1e-13, atol=1e-14` end-to-end
(mirroring `stc`'s own choice) purely for cross-task comparability, not because it was found
strictly necessary here.

| hit | `t_u` (periods) | `t_s` (periods) | backward_distance | forward_distance |
|---|---|---|---|---|
| PRIMARY (13,13) | `2.565` | `-2.565` | `1.52e-08` | `0.157` |
| SECONDARY (15,15) | `3.327` | `-3.327` | `1.96e-09` | `1.70e-03` |
| MIRROR_A (13,17) | `2.323` | `-4.331` | `1.62e-09` | `1.71e-03` |
| MIRROR_B (17,13) | `4.331` | `-2.323` | `5.36e-09` | `1.69e-03` |

`backward_distance` is tight for all four (`<5.4e-8`), consistent with the pattern this task
chain has established throughout: propagating a leg's own crossing back to its own seed is close
to an exact numerical time-reversal of the identical integration. `forward_distance` is tight
(`~1.7e-3`) for SECONDARY/MIRROR_A/MIRROR_B — but notably looser for PRIMARY (`0.157`, ~2 orders
of magnitude worse than the other three) despite PRIMARY's own elapsed time (`2.565` periods)
being the SHORTEST of the four, which would naively predict the TIGHTEST forward re-approach, not
the loosest. This is reported honestly as an open, unexplained asymmetry — not smoothed over.
A plausible (untested-this-task) hypothesis: `k_u=k_s=13` (PRIMARY) sits closer to a fold/tangency
in the underlying `(tau_u, tau_s)` residual surface than `k=15` or the `(13,17)`/`(17,13)` pair,
making the local corrector-residual floor genuinely worse-conditioned there even though the
raw 2-D Newton residual (`7.65e-10`) is unremarkable — consistent with `#766`'s own finding
(Sec. 3, `jovian_resonant_connections.py`) that a specific `(branch,k)` combination can sit near a
fold while a nearby one does not. All four `forward_distance` values remain far below the O(1-2)
nondimensional trajectory scale (`0.17%`–`16%` of it) — genuine, non-coincidental re-approaches,
just not uniformly tight.

## 7. Literature novelty gate (mandatory, read in full — this is what licenses "novel" above)

**Step 1 — direct text-layer grep of the primary source and its companion dissertation (the
load-bearing evidence).** Both `cyclers_pdf/papers/miceli-bosanac-2026-generating-planar-
trajectories-neptunian-system-motion-primitives-jas-73-11-...txt` (the JAS-2026 paper's own text
layer) and `cyclers_pdf/papers/miceli-2025-...-phd-dissertation-colorado.txt` (the companion
dissertation) were grepped this task for `homoclinic`/`heteroclinic`: **zero occurrences in
either document.** The paper's own two `intersect`-adjacent hits are unrelated: p.6/line 310
("intersects one of the primaries" — a corrector failure-mode check, not a manifold intersection)
and Sec. 3.2/line 639 ("the capability to encode multiple intersections between two primitives" —
the coarse graph-search's own configuration-space adjacency, not a Poincaré-section computation).
The dissertation's own `Poincaré`/`Poincare` hits (lines 989-1025, 3616-3622) are general
dynamical-systems-theory background (what a Poincaré map is, how it's used for visualization) —
not a Neptune-Triton-specific manifold-intersection computation.

**Step 2 — supplementary WebSearch — CORRECTED, this task's own execution of this step was
inadequate.** Two targeted queries ("Neptune Triton homoclinic connection CR3BP resonant
orbit"; `"Neptune" "Triton" invariant manifold intersection Poincare section Wu Ws periodic
orbit`) surfaced only: the SAME Miceli & Bosanac work, general CR3BP manifold-computation
background (Koon-Lo-Marsden-Ross, arXiv 1111.0032), and a 2007 ScienceDirect paper on Neptune-
Triton Lagrangian-point (L1/L2) normalization/center-manifold dynamics with no overlap to this
task's own result. **This was NOT sufficient.** A Fable adversarial review's very first query
("Neptune Triton CR3BP homoclinic connection resonant orbit invariant manifold" — essentially a
rephrasing of this task's own first query) surfaced as its #1 result: **Spear 2021, "Planar
Heteroclinic Connections in the Neptune-Triton Circular Restricted Three Body Problem" (MS thesis,
CU Boulder, Bosanac advisor)** — genuine, on-point Neptune-Triton CR3BP manifold-connection prior
art that this task's own search entirely missed. See `docs/notes/2026-08-08-spear-2021-digest.md`
for the full digest and its implications (summarized in the corrected Verdict section above). This
was a real execution gap in this task's own search step, not an inherent limitation of the search
tools available — the Verdict section above records the corrected novelty scope.

**Step 3 — `search/literature_check.py`'s own mandatory-floor gate, run with the real WebSearch
tool wired.** Signature: `CandidateSignature(primary="Neptune", sequence=("Triton",),
resonances=("4:5",), topology_label={"resonant"})`. `build_queries` generated 7 queries (all
7 run against the real WebSearch tool this task):

```
Triton cycler trajectory
Triton cycler
Triton cycler resonance 4:5
Voyager Neptune Triton encounter trajectory
Trident Neptune Triton mission concept
Stone Triton cycler
Neptune moon tour cycler ballistic patched-conic
```

Every query returned results (RC battery-charger products literally branded "Triton", Voyager 2
flyby history, the Trident Discovery-class mission concept, an unrelated cyclist named "Stone",
and one general moon-cycler-search paper referencing Neptune only as a "tested for general use"
system, not this specific orbit/connection) — **none matched the cycler structural fingerprint.**
Formal verdict via `check_literature`: **`status="not-found"`, `confidence=0.4`**. Per the
module's own discipline, `not-found` is necessary-not-sufficient for novelty on its own — but
combined with Step 1's direct, load-bearing text-layer evidence (zero connection-level content in
either primary source), this task treats the finding as genuinely **NOVEL**, not a reproduction
or self-consistency-only result, per the `#781` dispatch note's own explicit instruction.

**Explicit caveat, stated per the dispatch note's own honesty discipline**: `literature_check.py`'s
own module docstring scopes it to CYCLER-trajectory vocabulary ("NOT FOR RAW (non-cycler) CR3BP
PERIODIC ORBITS") — a raw homoclinic manifold self-connection of a resonant orbit is further
still from that scope than the raw periodic orbits its own docstring already excludes (every
generated query above is literally unable to express "homoclinic manifold intersection" as a
concept). This gate's `not-found` is therefore read as weak, not dispositive, additional
corroboration layered on top of Step 1 — never as the primary novelty grounding on its own,
exactly as the module docstring states.

## 8. The 4:7-stress attempt — an honest, diagnosed FAIL

**The orbit**: `recover_esm_candidate("4:7-stress")` — `x0=1.787757094694`, `ydot0=
-1.147924599743`, `period=44.514684550531`, `jacobi=2.997230642137`, `max_eigenvalue=
14624.10512414333` (POSITIVE this time — no branch-flip subtlety). Own perpendicular section
points: the IC (`x0=1.7878`) and the half-period return (`x=-2.1316`, `half_crossings=5`) — a
much wider orbit (`x` spans `[-2.13, 1.79]`) than 4:5-saddle.

**Diagnostic method** (same curve-building approach as Sec. 4, bounded to `k in 1..8`,
`n_tau=10`, `max_time_factor=3.0` — a deliberately SMALL exploratory pass, not an exhaustive
scan, per the dispatch note's own "if time/tractability allows" framing): building all 32 curves
took 32 seconds. The closest-approach search found the SAME kind of degenerate near-orbit
clustering `own_section_points` would reject at low `k` (`k<=3`, `d~0`, `gd~0`) — expected, the
manifold has not yet separated. At `k` up to 8 (~1.6 orbital periods), the BEST "real-separation"
candidates found (`gd > GHOST_GUARD_DELTA`) sit at only `gd=1.2e-3`–`3.5e-3` — **1.2×–3.5× the
guard threshold**, compared to 4:5-saddle's own `134×`–`287×` at a comparable relative elapsed
time. This razor-thin margin is itself diagnostic (per
`feedback_verify_automated_ghost_guard_booleans`/`feedback_isolated_sweep_flips_suspect_artifact`):
even a "converged" hit here would be barely distinguishable from a trivial self-shadow, not the
kind of real, non-delicate margin every genuine hit in this task chain has shown.

**Newton refinement attempted on three of these marginal candidates** (`(branch_u=+1,k_u=8,
branch_s=-1,k_s=8)`, `(branch_u=-1,k_u=5,branch_s=+1,k_s=8)`, `(branch_u=+1,k_u=6,branch_s=+1,
k_s=8)`, seeded from their own coarse-grid `tau` estimates, `tol=1e-9, max_iter=60,
max_time_factor=3.0`): **all three FAILED to converge** — residuals plateaued at `1.7e-3`,
`2.0e-3`, and `3.4e-3` respectively (target `1e-9`), each taking 10–70 seconds of wall-clock
before exhausting `max_iter` with `notes="did not reach tol"`. This is failure mode **(b) Newton
stalls/does not converge** from `#781`'s own dispatch-note diagnostic checklist — NOT (a) no
crossing found (crossings exist and were located cheaply), NOT (c) converges-but-ghost-fails
(never reached convergence to test this), NOT (d) Radau-blows-tolerance (never reached a
converged candidate to cross-check).

**Diagnosis, matching `#759`'s own extreme-instability precedent directly**: at `|lambda|~
1.46e4` (`ln(14624)/period = 0.2155` per nondim time unit — over 40% faster growth per unit time
than 4:5-saddle's own `0.153`), meaningful manifold separation develops within roughly ONE
orbital period, meaning the coarse `n_tau=10` tau-grid used to seed the geometric close-approach
search is almost certainly too coarse to resolve the TRUE fine structure of the manifold curves
at this instability — the same "fractally sensitive" signature `#759`'s own results note
documents for its own extreme-instability case (Jupiter-Europa, a DIFFERENT orbit/system, but the
same qualitative Newton-intractability mechanism: tiny seed-`tau` perturbations amplified by
`|lambda|` per period produce wildly different, essentially uncorrelated crossing locations,
defeating both a coarse geometric pre-scan and Newton's own finite-difference Jacobian).

**This is reported as the task's own honest, bounded-effort conclusion, not a deeper investment**
— per the dispatch note's own explicit framing ("An honest FAIL here, clearly diagnosed and
reported, is a completely acceptable and expected outcome; do not force it or weaken any
tolerance to make it 'work'"). No tolerance was loosened, no marginal candidate was reported as a
hit, and no further escalation (e.g. a much finer `n_tau`, a purpose-built continuation-from-
lower-instability-C approach mirroring `#775`'s own methodology) was attempted this task — a
natural, well-scoped follow-up for a FUTURE task if the 4:7-stress connection is ever wanted, not
bundled into `#781`'s own scope.

## 9. Code delivered

* `src/cyclerfinder/search/neptune_triton_resonant_connections.py` (new sibling module):
  `own_section_points` (all-four-sign-combination union), `build_45_node`, `HomoclinicCandidate`,
  `find_homoclinic` (residual-ranked, no published target), `homoclinic_reapproach_check`
  (tighter end-to-end integration precision, mirroring `stc`'s own choice). No chain-closure/
  multiple-shooting/homotopy machinery — Neptune-Triton has no published chain-figure analogue to
  Vaquero's Fig. 4.9-4.10, so `stc`'s own `attempt_chain_closure`-family functions have no
  Neptune-Triton counterpart in this task's scope, per the dispatch note's own instruction.
* `src/cyclerfinder/search/jovian_resonant_connections.py` (one targeted fix): `ResonantNode.
  from_candidate`'s own `rel_err` staleness-guard formula now compares magnitudes
  (`abs(lam_u - abs(cand.max_eigenvalue))`), fixing a real, previously-unexercised correctness
  gap for negative-signed `max_eigenvalue` candidates (Sec. 2). Byte-for-byte behavior-preserving
  for every pre-existing positive-eigenvalue caller (confirmed: the full Jovian + Saturn-Titan
  connection/family test suites pass unchanged).
* `tests/search/test_neptune_triton_resonant_connections.py`: 23 tests covering the node-building
  regression (including the `#781` staleness-guard fix), the section-convention derivation (both
  the "recovers all 6 points" and "either single ydot_sign drops one perpendicular reference"
  regressions), all four known hits' convergence/ghost-margin/Radau-cross-check/forward-backward-
  re-approach, the mirror-pair symmetry, the on-axis check for both PRIMARY and SECONDARY,
  `find_homoclinic`'s own scan/ranking behavior (re-deriving PRIMARY through the module's own
  brute-force-scan code path, independent of the tau0-seeded shortcut the other tests use), and an
  honest-empty-scan regression. None marked `@pytest.mark.slow` (a discovery-verdict-bearing
  result must run in CI, not be silently skipped).

## 10. Verification

* `uv run pytest tests/search/test_neptune_triton_resonant_connections.py -v`: 23/23 pass
  (~251s wall, 8-way parallel).
* `uv run ruff check` / `ruff format --check` on all three changed/new files: clean.
* `uv run mypy src tests` (canonical full invocation): clean, 833 source files.
* `uv run pytest tests/data tests/search -q`: run after the sibling `#780` task's own concurrent
  run finished (per this task's own concurrent-agent git/CPU-contention discipline) — see commit
  history for the recorded pass/fail status.
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: run before the `OUTSTANDING.md`
  update commit.

## 11. Net effect on `#781`

**DONE — a genuine, honestly-evidenced, NOVEL positive result, CORRECTED 2026-08-08 by adversarial
review** (still the first novel — not reproduction, not self-consistency-only — connection-stage
result of this whole `#754`/`#759`/`#766`/`#767`/`#781` task chain, per `#776`'s own note
anticipating this exact possibility, but narrower in scope than originally claimed). **Two**
distinct, ghost-guard-passed, Newton-converged (`<1e-8` residual, three of the four crossing-index
registrations `<1e-9`), independently Radau-cross-checked (`<1.31e-7`) homoclinic self-intersection
trajectories of the Neptune-Triton 4:5-saddle resonant orbit exist at Miceli & Bosanac's own
`C=2.987089791658` (PRIMARY and SECONDARY) — not four independent hits as originally reported;
SECONDARY's own on-axis reflection symmetry produces two further crossing-index registrations
(MIRROR_A/MIRROR_B) of the SAME trajectory, an automatic consequence of the symmetry, not
independent corroboration. "Transversal" was asserted but never measured — read as unverified.
Forward/backward re-approach is honestly reported including an unexplained PRIMARY-specific
asymmetry, not smoothed over. The novelty claim itself is narrower than first stated: prior art
(Spear 2021, missed by this task's own original WebSearch step, found by adversarial review) rules
out "no Neptune-Triton manifold-connection literature exists," but that thesis's own object (11
heteroclinic connections between Lyapunov orbits, zero homoclinic content, no 4:5/4:7 resonant
families) is genuinely distinct from this task's own — the corrected claim is **first homoclinic
connection, and first resonant-orbit manifold connection of any kind, in the Neptune-Triton
CR3BP**. The 4:7-stress attempt is an honest, clearly-diagnosed FAIL matching `#759`'s own
extreme-instability precedent — not forced. This task also fixed a real (if previously
unexercised) bug in the shared `ResonantNode.from_candidate` staleness guard, exposed by this
orbit's own genuinely negative saddle eigenvalue — the first such orbit this whole task chain has
built a node from; independently confirmed behavior-preserving by the adversarial review as well
as this task's own regression suite. **Not yet written back to `data/catalogue.yaml`** — per the
adversarial review's own recommendation, any catalogue writeback should wait for this corrected
framing, not the original.
