# `#757` — Re-scoping Task B (`#754`) around the confirmed families: what Tables 2/3 actually require

**Task:** `#757` (research/scoping only — no code, no catalogue changes). Question: with
`5:6-LO` now a clean negative after three search passes (`#753`/`#755`/`#756`) and only
`3:4-LO` (reviewer-confirmed, `#755`) + `5:6-LI` (clean gate pass, `#753`) in hand, which of
Anderson & Lo 2011's Table-2 (homoclinic) and Table-3 (heteroclinic) connection gates are
actually reachable, and what should `#754` become?

**Sources read directly this pass** (not just digests): Anderson & Lo 2011 (JAS 58:167) —
full text layer of
`cyclers_pdf/papers/anderson-lo-2011-dynamical-systems-resonant-flybys-ballistic-case-jas-58-167-doi-10.1007-BF03321164.txt`
around the Poincaré-map methodology (journal pp. 170-171), family-selection prose
(pp. 183-184), Table 1 + footnote 4 (p. 184), and the "Using Homoclinic and Heteroclinic
Connections" section with Tables 2/3 (pp. 190-191), **confirmed against the rendered PDF
pages 24-25** (journal pp. 190-191) via direct page reads. Code read in full:
`src/cyclerfinder/genome/heteroclinic_cycle.py` (all 657 lines),
`src/cyclerfinder/search/jovian_resonant_families.py` (module constants, candidate
dataclass, seed table, gate machinery), plus the `#756` sweep checkpoint
`data/found/756_jupiter_europa_5_6_lo_relaxed_period/candidates.jsonl`.

---

## 1. What the paper's Tables 2 and 3 actually are

### The decisive footnote (p. 184, footnote 4 — verified in text layer, line 1161)

> "The 3:4-LO and 5:6-LO resonant orbits will be referred to as the 3:4 and 5:6 resonant
> orbits hereafter."

Every subsequent "3:4 orbit" / "5:6 orbit" in the connection sections means **3:4-LO** and
**5:6-LO** specifically. The selection prose immediately before Table 1 (pp. 183-184) is
explicit about how the flavors were chosen:

> "the 5:6-LI orbit is only slightly unstable, so it was **removed from consideration**. The
> eigenvalues of both the 5:6-LO and 5:6-NO orbits indicate that they are very unstable. A
> comparison of these two types of orbits with the flyby trajectory and a preliminary
> analysis using Poincare sections indicated that the **5:6-LO orbit was most relevant** to
> this problem."

So: **5:6-LI — the family we confirmed — is the one flavor the paper itself explicitly
discarded** before computing any manifolds or connections. 5:6-NO was also discarded (in
favor of 5:6-LO). No connection in the paper involves 5:6-LI or 5:6-NO.

### Table 2 (p. 190): "Homoclinic Trajectory State at Intersection"

| x | y | ẋ | ẏ |
|---|---|---|---|
| −1.28427733 | 0.0 | 0.00000009 | 0.46372205 |

Verified against the rendered PDF page: the four column headers are exactly *x, y, ẋ, ẏ*
(nondimensional planar rotating-frame CR3BP state), one data row, C = C_flyby =
2.99163956830415. The surrounding prose (pp. 190-191):

- It is the intersection of the **stable and unstable manifolds of the 3:4 orbit
  (= 3:4-LO) with itself** — "the intersections of the stable and unstable manifolds of a
  given resonant orbit are points that will approach that orbit when they are integrated
  into the future as well as the past."
- "In both cases, the homoclinic trajectory approaches the **3:4 resonant orbit**"
  (Fig. 25 forward + backward panels). **`#745`'s digest claim is independently CONFIRMED:
  Table 2 involves 3:4-LO alone. No 5:6 family of any flavor is needed.**
- The intersection point itself happens to sit "almost exactly at the location of the 5:6
  orbit with a difference in x position from the 5:6 orbit of approximately 8.0 × 10⁻⁵"
  (p. 184, re-stated at p. 190 as "intersect almost at the 5:6 orbit"). This is a
  geometric coincidence-of-interest at this energy, **not** a dependency — the 5:6-LO
  orbit plays no role in computing the Table-2 connection. (It does hand us a free gift;
  see §4.)
- The intersection was "computed using **interpolation** with the closest points on the
  invariant manifolds" — NOT a Newton-corrected true intersection. The published digits
  are therefore interpolation-limited (manifold discretization + their ε ≈ 0.5×10⁻⁵
  offset), which bounds how tight a reproduction gate can honestly be (see §5 gate spec).

### Table 3 (p. 191): "Heteroclinic Trajectory State at Intersection"

| x | y | ẋ | ẏ |
|---|---|---|---|
| −1.43029175 | 0.0 | 0.00018678 | 0.67262261 |

Same format, verified against the rendered page. The prose is unambiguous:

> "an intersection of the **unstable manifold of the 3:4 orbit** and the **stable manifold
> of the 5:6 orbit** near the 3:4 orbit was selected... the forward integrated trajectory
> approaches the 5:6 orbit along its stable manifold, and the backward integrated
> trajectory comes even closer to the 3:4 orbit along its unstable manifold."

With footnote 4, this is exactly: **Wu(3:4-LO) ∩ Ws(5:6-LO)**, on the section, near the
3:4-LO orbit's own section point. **Table 3 requires 5:6-LO — the one family that is a
clean negative after three dedicated search passes. It is NOT reachable with the confirmed
families.** 5:6-LI is not a legitimate stand-in (the paper explicitly removed it), and see
§3 for why it is also *practically* unusable as a manifold node.

### What a "connection state" numerically consists of (for gating)

- **Section**: {y = 0} along the **negative x-axis opposite Europa**, **one-sided ẏ > 0**
  (pp. 170-171: "the surface of section is defined by y = 0 along the negative x-axis
  opposite Europa... a one sided Poincare map with y > 0 was used" — the "y > 0" there is
  ẏ > 0, per Eq. (7) context: on the section x and ẋ are the free coordinates and ẏ is
  recovered from the Jacobi constant with the + sign). Both table rows have ẏ > 0 and
  x < 0 — consistent.
- A connection state is therefore **2 independent numbers (x, ẋ)** at fixed
  C = 2.99163956830415, with y ≡ 0 and ẏ derived via the paper's Eq. (7):
  ẏ = +sqrt(x² + 2(1−µ)/r₁ + 2µ/r₂ − C − ẋ²) (planar). The published ẏ column is a
  redundancy/consistency check, not a third degree of freedom.
- Manifold computation parameters (for faithfulness): eigenvector offsets taken in **both
  ± directions**, magnitude ≈ **0.5 × 10⁻⁵** (~3.36 km) (p. 181-182); trajectories
  generated with RKF7(8); intersections found by interpolating the two manifolds' section
  curves; the section "actually intersects the unstable orbit" (p. 181) — i.e. the
  resonant orbit itself crosses the section, unlike the typical libration-orbit setup.

## 2. Code-side: what the `#754` generalization actually requires

### `genome/heteroclinic_cycle.correct_connection` (read in full)

The solver is **duck-typed** on its node objects: the body only touches `.state0`,
`.period`, `.jacobi`, `.label`, `.unstable_eigvec`, `.stable_eigvec` — it never calls
`LyapunovNode.from_libration`. Concretely:

- **"Resonant-member node type"** = a `ResonantNode` dataclass (same six fields) with a
  `from_candidate()` constructor that takes a `ResonantFamilyCandidate` and computes the
  saddle eigenpair via the module's existing `_planar_floquet_pair(system, state0, period)`
  (one STM propagation). Near-zero change to the solver itself. The
  `_planar_floquet_pair` sign caveat (magnitudes only, positive-real saddle assumed) is
  satisfied: 3:4-LO's Barden eigenvalue is real positive (1036.116117,
  `is_real_unstable=True`, `#755`) — assert this at node-construction time.
- **"One-sided {y=0} section"** already exists: `_section_crossing` takes
  `ydot_sign` (threaded through as `ydot_sign_u`/`ydot_sign_s` kwargs of
  `correct_connection`). The paper's section = `ydot_sign=+1` on both legs. What is
  genuinely missing is the **negative-x restriction** (the paper's section is the negative
  x-axis only): add an optional `x_sign: int | None` filter to `_section_crossing`
  (same pattern as `ydot_sign`, ~5 lines, threaded through `_connection_residual` /
  `correct_connection` / `crosscheck_cycle`).
- **"Homoclinic A=B mode"** is *structurally* already allowed — `assemble_cycle`'s own
  docstring says "a single node is a degenerate homoclinic cycle O->O", and
  `correct_connection(system, node, node)` passes the Jacobi equality trivially. What is
  genuinely missing is the **trivial-solution ghost guard**: for A=B, a "converged
  connection" can be the orbit shadowing itself (seed still ε-near the orbit crossing the
  section at the orbit's own section point). Required guard: reject any converged crossing
  whose (x, ẋ) lies within δ of the *orbit's own* section-crossing set (computable by
  propagating the orbit one period and collecting its qualifying section points), plus the
  existing minimum-time floor. For Table 2 this guard has huge margin: the genuine
  homoclinic point (x ≈ −1.2843) is ~0.146 away in x from 3:4-LO's own section point
  (x ≈ −1.4304). Port the `#701`/`#702` seed-anchored-reference discipline as documented
  in `ccr4bp_heteroclinic_search` (discipline only, not architecture — `#752`'s "wrong
  donor" verdict re-confirmed on this read).
- `epsilon` default is 1e-6; the paper used 0.5e-5 — pass explicitly (same order; fine).
- The independent-integrator cross-check (`crosscheck_cycle`, Radau vs DOP853) already
  exists and needs only the same `x_sign` threading.

### `search/jovian_resonant_families.py` (read: constants, dataclass, seeds, gate)

- Confirmed candidates are stored as **hardcoded seed tuples**
  `_TABLE1_CANDIDATE_SEEDS = {"5:6-LI": (-0.374722, 1.0, 52, 37.6990), "3:4-LO":
  (-1.430408, 1.0, 6, 25.6725), ...}` re-converged on demand by
  `recover_table1_candidate(label)` → `ResonantFamilyCandidate`.
- `ResonantFamilyCandidate` stores `x0, ydot0, period, jacobi, ...` and eigen**values**
  (`barden_eigenvalue`, `max_eigenvalue`, `planar_floquet_eigenvalue`) but **NOT
  eigenvectors**. No schema change needed: `state0 = (x0, 0, 0, 0, ydot0, 0)` + `period`
  fully determine the monodromy, so `ResonantNode.from_candidate` just recomputes the
  saddle pair with `_planar_floquet_pair` — deterministic, one STM propagation.

## 3. Verdicts

### Table 2 — **BUILDABLE NOW** with only the confirmed 3:4-LO

All ingredients exist; the deltas are the three small items above (node adapter, `x_sign`
filter, homoclinic ghost guard) plus a driver + gate. Two bonuses:

1. **The Table-2 gate doubles as an independent, state-space corroboration of the `#755`
   reviewer ruling** that our 3:4-LO *is* the paper's 3:4-LO — an axis (manifold
   intersection geometry) fully orthogonal to the eigenvalue/period evidence the ruling
   weighed. If our family's manifolds self-intersect at the published (x, ẋ) to ~1e-4,
   the identification is corroborated a third way; if they don't, that is genuine,
   reportable counter-evidence. Either outcome is valuable.
2. Already-in-hand consistency check: our 3:4-LO's own section IC x0 = −1.430408 sits
   1.16×10⁻⁴ from Table 3's x = −1.43029175 — squarely matching the paper's own statement
   (p. 184) that the 5:6 manifolds fold/intersect "near the 3:4 resonant orbit... the
   difference in x is approximately 2 × 10⁻⁴."

### Table 3 — **NOT reachable with confirmed families**; do NOT substitute 5:6-LI

Table 3 is Wu(3:4-LO) ∩ Ws(**5:6-LO**). 5:6-LO is a clean negative (`#753`/`#755`/`#756`).
Substituting 5:6-LI would not be a reproduction of anything the paper published — and it is
also **practically unworkable**: 5:6-LI's max eigenvalue is 1.000008, so the manifold
e-folding time is ~T/ln λ ≈ 37.7/8×10⁻⁶ ≈ 4.7×10⁶ nondimensional time units
(~10⁵ orbit periods) — there is effectively no computable manifold to intersect. This is
precisely why the paper removed it from consideration. **Recommendation: no exploratory
3:4-LO↔5:6-LI demonstration either** — it would be neither a reproduction nor a feasible
computation. (An exploratory homoclinic of 5:6-LI is ruled out for the same reason.)

## 4. NEW FINDING — the paper indirectly publishes a digit-grade 5:6-LO seed

Combining Table 2 with the p. 184 prose gives coordinates no prior search pass had:

- The Table-2 intersection (x = −1.28427733, ẋ = 0.00000009) lies "almost exactly at the
  location of the 5:6 orbit" on the section, "difference in x position... approximately
  8.0 × 10⁻⁵". Therefore **the 5:6-LO orbit's own section point is pinned at
  x ≈ −1.28427733 ± ~8×10⁻⁵, ẏ ≈ +0.4637, with small |ẋ|** (the intersection's ẋ is
  9×10⁻⁸; the orbit's exact ẋ there is unstated but visually coincident in the paper's
  (x, ẋ) section plots).
- All Table-1 families were generated by the paper's own *symmetric* corrector, so
  5:6-LO is x-axis-symmetric; if this section point is (near) one of its perpendicular
  crossings, `correct_symmetric_fixed_jacobi` seeded at x0 ≈ −1.2843, ydot0_sign=+1
  reaches it directly; if the crossing is non-perpendicular, the existing general
  asymmetric fixed-C corrector (`cr3bp_general_periodic`) can be seeded with the full
  4-state (−1.2843, 0, ~0, +0.4637).
- **This region is genuinely unexplored at digit precision**: the `#756` checkpoint
  (`candidates.jsonl`, all 159 converged candidates) contains NOTHING in
  x0 ∈ (−1.35, −1.20) except neutral λ≈1.0000000 island orbits at x0 ≈ −1.335, −1.245,
  −1.202. Every serious 5:6-LO hunt worked the x0 ≈ −1.42…−1.43 hotspot (by analogy with
  3:4-LO) or coarse wide grids. This is exactly the "genuinely different seed strategy"
  `#756`'s own closing opinion called for — and it is sourced, not guessed.

## 5. Concrete task specs

### Task `#754` (re-scoped) — "3:4-LO homoclinic connection + Anderson-Lo Table-2 gate"

Spec-complete, Sonnet-tier (`[[feedback_subagent_model_tiering]]` — spec-complete TDD with
a deterministic sourced gate). No new dynamics; planar CR3BP only.

1. **`genome/heteroclinic_cycle.py` (minimal edits)**: add optional `x_sign: int | None`
   to `_section_crossing` (same pattern as `ydot_sign`: skip crossings whose
   `sign(x) != x_sign`), thread through `_connection_residual`, `correct_connection`, and
   `crosscheck_cycle`. No behavior change when `None` (default) — existing W-Z tests must
   stay green untouched.
2. **New module `search/jovian_resonant_connections.py`**:
   - `ResonantNode` dataclass (`label, state0, period, jacobi, unstable_eigvec,
     stable_eigvec, converged` — same field names as `LyapunovNode` so
     `correct_connection` accepts it unchanged) with
     `from_candidate(system, cand: ResonantFamilyCandidate)` using
     `heteroclinic_cycle._planar_floquet_pair`; raise if `not cand.is_real_unstable`
     (sign-stripping caveat) or if `|λ|` disagrees with `cand.max_eigenvalue` by >1e-6
     relative.
   - `own_section_points(system, node, *, ydot_sign=+1, x_sign=-1)`: the orbit's own
     qualifying section crossings over one period (for the ghost guard and reporting).
   - `find_homoclinic(system, node, *, epsilon=0.5e-5, ...)`: coarse scan over
     `branch_u, branch_s ∈ {+1,−1}²` and `k_u, k_s ∈ 1..6` (paper gives no k/branch;
     document the scan), each via `correct_connection(node, node,
     ydot_sign_u=+1, ydot_sign_s=+1, x_sign=-1, epsilon=0.5e-5)` with a generous
     `max_time_factor` (3:4-LO period ≈ 25.67; the manifolds traverse to x ≈ −1.28);
     **ghost guard**: reject any converged crossing within δ = 1e-3 (nondim, in (x, ẋ)
     norm) of `own_section_points`, and reject legs whose qualifying-crossing time is
     below the existing t-floor. Return ALL surviving converged connections, ranked by
     distance to the Table-2 state.
3. **Gate (sourced, no self-reference)** — constants verbatim in the module with page
   cites: `TABLE2_STATE = (x=−1.28427733, xdot=0.00000009, ydot=0.46372205)` (p. 190).
   PASS iff a surviving connection has: Newton residual ≤ 1e-7; |x − x_T2| ≤ 1e-4 AND
   |ẋ − ẋ_T2| ≤ 1e-4 (tolerance justified by the paper's interpolation-not-Newton
   intersection + their ε = 0.5e-5 manifold offset — document this; do NOT chase 1e-8
   digits the source itself cannot support); ẏ recovered from Eq. (7) (+ branch) matches
   ẏ_T2 to the propagated tolerance; independent Radau re-derivation of both legs at the
   converged (τ_u, τ_s, k, branch) agrees to ≤ 1e-6 (existing `crosscheck_cycle`
   pattern). **Corroboration (reported, not gated)**: forward AND backward integrations
   from the intersection state both re-approach the 3:4-LO orbit (min orbit-distance
   metric over a few periods each way — the paper's Fig. 25 behavior); crossing x lands
   ~8×10⁻⁵-order from Table 2 per the p. 184 coincidence statement. Honest FAIL
   reporting per `#753`'s style — a miss is reported with the best-found state, never
   loosened to pass.
4. **Verification**: gate evidence test NOT `@pytest.mark.slow`
   (`[[feedback_delegation_fresh_agent_not_fork]]`); run `uv run pytest tests/data
   tests/search tests/scripts -q`, full `mypy src tests`, `ruff check` + `format --check`.

### New task (proposed, `#758`) — "5:6-LO recovery from Table-2's implied section point"

Small, cheap, high information value; independent of (and parallelizable with) the
homoclinic build — different files (`jovian_resonant_families.py` + note).

1. Fine symmetric-corrector sweep: x0 ∈ −1.28427733 ± 2×10⁻⁴ (dense, e.g. 2000 points,
   emphasizing ±8×10⁻⁵), `ydot0_sign=+1`, `half_crossings ∈ 1..12`, period_guess near
   2π·6 ≈ 37.70 but **no period filtering during search** (`#756` discipline; period is
   after-the-fact corroboration only, per the `#755` reviewer precedent).
2. If (1) misses: general asymmetric fixed-C corrector (`cr3bp_general_periodic`) seeded
   at the full 4-state (−1.28427733, 0, 0.0, +ẏ(Eq. 7)) and small perturbations, in case
   the pinned crossing is non-perpendicular.
3. Gate: Table-1 5:6-LO eigenvalue 4445.387515 at the module's existing 1e-3 relative
   tolerance (Barden-authoritative); corroboration: period_over_2pi vs 6,
   `europa_closest_approach` (expect a genuine close-flyby signature like 3:4-LO's
   ~1,641 km, NOT `#756`'s 12,000-28,000 km misses).
4. **If it passes, Table 3 (original full Task B scope) unlocks immediately** — the
   heteroclinic Wu(3:4-LO)∩Ws(5:6-LO) gate becomes a straight application of the same
   `find_connection` machinery (`correct_connection(node_34LO, node_56LO, ...)` with the
   Table-3 state (−1.43029175, 0.0, 0.00018678, 0.67262261) as the gate). If it fails,
   Table 3 stays honestly unreachable and is reported as such — no 5:6-LI substitution.

## 6. Overall recommendation

Dispatch **two tasks** (parallelizable, disjoint files):

1. **`#754` re-scoped** to the Table-2 homoclinic build above — genuinely well-founded
   *today* with only confirmed 3:4-LO, carries a sourced digit-grade gate, AND
   independently stress-tests the `#755` reviewer ruling.
2. **`#758`** — the Table-2-derived 5:6-LO seed retry. This is the first 5:6-LO strategy
   with *sourced digit-grade coordinates* rather than analogy/grid guesswork, and it is
   the sole gate-keeper for Table 3.

Hold the Table-3 heteroclinic gate itself until (2) resolves; if (2) is negative, `#754`'s
heteroclinic half remains blocked on 5:6-LO and should stay held — explicitly do **not**
build a 3:4-LO↔5:6-LI "demonstration" (non-reproductive AND computationally infeasible,
λ = 1.000008 ⇒ manifold e-folding ~4.7×10⁶ time units).
