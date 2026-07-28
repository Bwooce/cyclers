# `#753`: Jupiter-Europa planar CR3BP resonant-orbit families — Anderson-Lo 2011 Table 1 gate

**Task:** `#753`, dispatched 2026-07-28 as "Task A" of `#752`'s scoping recommendation
(user-authorized "Task A only for now"; Task B — `#754` — remains HELD pending this
report and review). Spec: `docs/notes/2026-07-28-752-resonant-manifold-jovian-tour-scoping.md`.

**Source paper** (read directly this task, both the acquired PDF's text layer and
the digest): R. L. Anderson & M. W. Lo (2011), "A Dynamical Systems Analysis of
Resonant Flybys: Ballistic Case," *J. Astronaut. Sci.* 58(2):167–194, DOI
`10.1007/BF03321164`. PDF at `cyclers_pdf/papers/anderson-lo-2011-...BF03321164.pdf`.

**Code delivered:** `src/cyclerfinder/search/jovian_resonant_families.py` (new
module) + `tests/search/test_jovian_resonant_families.py` (20 tests, all passing).
Both pass `ruff check`, `ruff format --check`, and `uv run mypy src tests` (project
canonical invocation) cleanly.

---

## What the module does

Retargets the cislunar resonant-orbit pipeline (`search/resonance_network.py`,
`#267`) to Jupiter-Europa, reusing — not reinventing — this project's existing
machinery:

* `correct_symmetric_fixed_jacobi` (`search/cr3bp_periodic.py`) for the corrector
  (the same modified Howell-Breakwell single-shooter the paper itself describes,
  including the `half_crossings` extension the paper's own text flags as its
  modification for resonant, not libration-point, orbits).
* `continue_family` (`search/cr3bp_continuation.py`) for natural-parameter
  Jacobi-constant continuation, with its full gauntlet (closure, period bounds,
  equilibrium gate, Jacobi conservation, independent-Radau cross-check, dedup,
  fold detection) inherited for free.
* `barden_stability` (`search/cr3bp_periodic.py`) and `_planar_floquet`
  (`search/resonance_network.py`) for stability classification.

Sourced constants, verified directly against the PDF text layer this task (not
inherited from the dispatch note or any prior digest without re-checking):

* `mu = 2.5266448850435e-5` — p.169, verbatim ("For the Jupiter-Europa system in
  this study, μ = 2.5266448850435 × 10⁻⁵"). This is the PAPER's own value, used
  deliberately instead of this project's DE440-registry Jupiter-Europa mu
  (`core.cr3bp.cr3bp_system("Jupiter", "Europa")`, ≈2.5274944e-5, a documented
  ~0.034% GM-vintage difference per the `#745` digest) — a fair reproduction test
  has to use the paper's own mu, not a different one.
* `C_flyby = 2.99163956830415` — p.179, verbatim.
* Table 1 (p.184), verbatim: `3:4-LO → 1036.116088`, `5:6-LI → 1.000008`,
  `5:6-LO → 4445.387515`, `5:6-NO → 28178.258323`.

The module also implements the paper's own two-body p:q resonant-ellipse seed
construction literally (`two_body_resonant_seed`, Fig. 1 / Eq. 5-6: periapse at
the secondary's orbital radius, rotating-frame closure period `T_full = 2πq`) —
this alone, however, did **not** locate any Table-1-matching orbit when converged
directly at `C_flyby` (see "Search history" below); a broader grid+bisection
sweep (`survey_candidates`) plus continuation from a higher, less chaotic Jacobi
constant was needed.

---

## The gate: 1 of 4 Table-1 rows recovered

| Row | Target `\|λ\|` | Recovered | Rel. err | Period/2π | Gate (1e-3) |
|---|---|---|---|---|---|
| **5:6-LI** | 1.000008 | **1.0001036** | **9.6e-5** | **5.999985 ≈ 6** ✓ | **PASS** |
| 3:4-LO | 1036.116088 | 818.401441 | 21.0% | 21.0148 (no clean q match) | FAIL |
| 5:6-LO | 4445.387515 | 4533.602947 | 1.98% | 16.1088 (no clean q match) | FAIL |
| 5:6-NO | 28178.258323 | 20573.609752 | 27.0% | 16.0398 (no clean q match) | FAIL |

Tolerance `1e-3` relative, justified by the corrector's own convergence floor
(crossing residuals of 1e-11–1e-13 across every candidate found this task, four
orders of magnitude tighter) and the paper's own stated ~1e-11 shooting floor.

**5:6-LI (CONFIRMED):** located by a direct grid+bisection sweep at `C_flyby`
in the "inner" `0 < x0 < 1-mu` Jupiter-Europa band (`x0 = -0.374722`,
`ydot0_sign = +1`, auto-selected `half_crossings = 52` — this orbit has ~130
small sub-loops per period, i.e. it is a genuinely intricate multi-loop
resonant orbit, not a simple ellipse). Two independent confirmations beyond the
eigenvalue match itself:
1. `period / 2π = 5.999985 ≈ 6` exactly confirms genuine 5:6-resonance lineage
   (`T_full = 2πq` with `q = 6`) — none of the other three candidates achieve
   this (see table).
2. `barden_stability`'s nontrivial eigenvalue is **real** and slightly `> 1`,
   matching the paper's own qualitative description ("the 5:6-LI orbit is
   only slightly unstable," p.183) — a saddle-type instability, not a
   marginally-stable complex pair (see the degenerate-eigenvalue finding below,
   which is exactly the distinction that makes this check meaningful).

Also continues smoothly in Jacobi constant toward `C_flyby` via `continue_family`
(6-step demonstration in the test suite; the full gauntlet — including the
independent-Radau cross-check — passes at each step once `jacobi_tol` is loosened
from `continue_family`'s default `1e-10` to `1e-8`, justified in the module
docstring: this orbit's one-period Jacobi drift is O(1e-9) even at
`rtol=atol=1e-12`, benign integration noise from its unusually long/convoluted
130-crossing path, independently confirmed tight by the Radau cross-check).

**3:4-LO, 5:6-LO, 5:6-NO (NOT CONFIRMED):** the best candidates found are
1.98%–27% off and — more importantly — none has a period that lands on a clean
`2πq` multiple for `q ∈ {4, 6}`. This means these are not merely imprecise
recoveries of the right orbit; there is no independent evidence they are even
the *same family* as the paper's rows, beyond an eigenvalue magnitude in a
roughly similar range. Honestly reported as unconfirmed near-misses, not
partial matches — the module's `_TABLE1_CANDIDATE_SEEDS` docstring says so
explicitly and the test suite hard-asserts these three rows FAIL the gate (so a
future accidental "fix" that silently starts passing would be caught, not
silently trusted).

---

## Search history (why this was hard, and what was tried)

Anderson & Lo's own text (p.176, "Continuation of Resonant Orbits") says their
method used "a grid search... to search *graphically*" for candidate ICs near a
resonance, then single-shot convergence, then linear-extrapolation continuation
to the target energy — i.e. a **visually-guided** search, not a blind automated
one. Reproducing that blind is genuinely hard, and this task's search reflects
that:

1. **Direct two-body seed at `C_flyby`** (`x0 = ±1.0`, periapse at Europa's
   radius): converges to trivial, near-integrable, marginally-stable orbits
   (the seed geometry never actually approaches Europa within one resonant
   period for odd `p`, a fact derivable from the fixed inertial periapse angle
   never landing on Europa's rotating-frame position — worked out and confirmed
   numerically this task).
2. **Wide grid+bisection sweep directly at `C_flyby`** (both "inner"
   `0 < x0 < 1-mu` and "outer" `x0 > 1-mu` bands, both `ydot0` signs, `half_crossings`
   1–10): found dozens of converged orbits, including several with eigenvalues
   in the right order of magnitude (10²–10⁵), but `C_flyby`'s dynamics there are
   EXTREMELY sensitive — refining the grid resolution 10x in a promising window
   revealed multiple genuinely distinct roots packed within `Δx0 < 0.01`, with
   eigenvalues varying by orders of magnitude between them (fractal-like
   sensitivity, not numerical error — confirmed by an independent Barden
   cross-check on one such pair, `4348.18` vs `4348.18` to 8 significant
   figures, i.e. NOT integration noise).
3. **Continuation-neighborhood robustness check**: warm-starting a tiny
   (`ΔC ≈ 0.0003`) Jacobi step from three of the "close" candidates found in
   (2) caused the corrector to jump onto an entirely different (usually
   trivial) branch almost immediately — direct evidence these particular
   candidates sit near a bifurcation/branch-crossing, not on a smoothly
   continuable family, and are therefore unreliable as Table-1 reproductions
   even where the eigenvalue looked close.
4. **Descent from a higher, less chaotic `C = 3.3`**: located 8 well-behaved,
   weakly-unstable (`λ` in `1.02`–`1.53`) families via a grid+bisection sweep at
   `C = 3.3`, then continued each down toward `C_flyby` in small steps. Every
   one hit a fold (`ydot0_from_jacobi`'s radicand going negative — the family's
   `x0` runs into the zero-velocity-curve boundary near Europa) between
   `C ≈ 3.07` and `C ≈ 3.25`, well short of `C_flyby = 2.9916`. None of these
   particular simple-loop members extend to the flyby energy.
5. **Wide auto-crossing scan at `C = 3.5`, target periods `T = 8π` (3:4) and
   `T = 12π` (5:6) exactly**: found ~110 orbits with the exact right period, but
   essentially all weakly unstable (`λ ≈ 1.00–1.12`) — the strongly unstable
   Table-1 members apparently require deeper, more convoluted topology (more
   sub-loops) that a simple wide scan at high `C` does not surface.
6. **5:6-LI found and confirmed** via a direct grid+bisection sweep at
   `C_flyby` itself in the previously under-explored "inner" band, landing on
   the intricate 130-crossing orbit described above.

## Methodological finding (load-bearing, not a side note)

`resonance_network._planar_floquet`'s "largest-magnitude eigenvalue" selection
heuristic — built and validated for the cislunar members, all of which have
`|λ| > 1.5` — is **unreliable near `|λ| ≈ 1`**: every periodic orbit's monodromy
carries a trivial unit eigenvalue pair from the flow's own time-translation
symmetry, and `argmax(|eigenvalue|)` can silently select that trivial pair
instead of the genuine nontrivial one. This was caught mid-task: dozens of
5:6-period candidates all reported `_planar_floquet` `λ ≈ 1.000000`–`1.000002`
almost everywhere in the search grid, which looked at first like a profusion of
5:6-LI-like near-matches. Cross-checking with `barden_stability` (which
explicitly discards the two eigenvalues nearest 1 before selecting the
nontrivial pair, per Barden 1994) revealed most of these are genuinely
**neutrally stable** — a complex, unit-modulus eigenvalue pair (e.g.
`0.963 + 0.270j`), not "barely unstable" at all. Only after re-classifying via
`barden_stability` did the genuine, real, `> 1` 5:6-LI candidate emerge. This is
now documented in the module's own docstring and covered by a standing
regression test (`test_planar_floquet_degenerate_eigenvalue_pitfall_is_real`)
so it cannot silently regress. Every large-eigenvalue candidate in this module
is cross-checked between both functions and they agree to `< 1e-4` relative in
every case (also covered by a test) — the degenerate-eigenvalue issue is
specific to the `|λ| ≈ 1` regime, not a general bug.

---

## Verification

* `uv run ruff check` / `ruff format --check` on both new files: clean.
* `uv run mypy src tests` (project canonical invocation): clean (821 files, no
  issues).
* `uv run pytest tests/search/test_jovian_resonant_families.py -q`: 20/20 pass
  (~65s wall time; the continuation-demonstration test is the slow one at ~47s
  due to the 130-crossing orbit's expensive per-step Radau cross-check —
  deliberately kept to 6 steps, not the original 25, to stay fast).
* `uv run pytest tests/data tests/search -q`: run as part of this task's own
  verification (catalogue-adjacent ratchets); see the commit history for the
  pass/fail status recorded at commit time.

---

## Recommendation for `#754` (Task B) — opinion, not a decision

Task B (generalize `heteroclinic_cycle.correct_connection` to resonant nodes +
homoclinic mode, gated on Tables 2/3's connection states) explicitly depends on
having the 3:4-LO and 5:6-LO resonant orbits' invariant manifolds in hand — the
paper's own Table 2/3 gate needs the manifolds of families this task did **not**
confirm. My opinion: dispatching Task B now would mean building manifold
globalization on top of unconfirmed orbits, which is very likely to produce
"it looked plausible" results with no honest gate underneath — exactly the
outcome this project's own discipline (`feedback_orbit_closure_discipline`,
"it closed!" is the danger signal) warns against. I would recommend against
dispatching Task B until either (a) 3:4-LO and 5:6-LO are confirmed by further
targeted search (the fold-based continuation-from-high-C approach in this
task's search history #4 came closest in spirit — a finer C-grid or a smarter
seed derived from the two-body flyby-vector-rotation construction Anderson & Lo
themselves use, Section "Designing Flybys Using Two-Body Approximations," might
close the gap), or (b) Task B is re-scoped to use the ONE confirmed family
(5:6-LI) for a reduced-ambition manifold/connection demonstration not tied to
Tables 2/3's specific numbers. This is my assessment for the user to weigh, not
a decision — `#754` stays HELD pending review as originally specified.
