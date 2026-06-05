# M-3D — Full-3D / inclination-lift design (Venus i=3.39°, Mars i=1.85°)

**Status:** design draft (brainstorming, 2026-06-05). No code, no task plan yet.
**Contract:** `docs/superpowers/plans/2026-06-02-multirev-3d-vem-ephemeris-roadmap.md:69-84`
(M-3D — "Ready for task-level planning: NO — needs design first").
**Scope decision carried in:** Geometry = **Full 3D inclination**; idealized
circular model *may stay 2D as a fast pre-filter* (roadmap scope table,
`...roadmap.md:18`).
**Non-goals (locked, see §7):** planet-centric / CR3BP (task #76, T-P graph);
real-ephemeris blind discovery / TCM budget (M-ED); multi-rev Lambert solver
math (M-L).

---

## 0. Headline finding (the brief's premise is half-built already)

The task brief expected the coplanar assumption to be spread across the
ephemeris, Tisserand, optimiser and construct layers. **Most of that 3-D
substrate already exists and is sourced.** Auditing read-only on 2026-06-05:

- The **real-ephemeris path is already fully 3-D.** `_AstropyBackend.state`
  returns true heliocentric DE440 vectors with a non-zero z-component, rotated
  ICRS→J2000-ecliptic (`core/ephemeris.py:210-237`); the rotation comment that
  it is "bundled" is the only inaccuracy and is already corrected in the roadmap
  note. There is **no coplanar shortcut in the DE440 data layer.**
- The **Lambert solver is already 3-D-native** — it takes 3-D `r1, r2` and
  returns 3-D `v1, v2` with no in-plane projection (`core/lambert.py:1-7`,
  signature `core/lambert.py:91-115`; `construct_cycler` feeds it whatever
  `ephem.state` returns, `search/construct.py:124-156`). Construction carries
  z end-to-end; no 2-D flattening anywhere in `construct.py`.
- An **inclined *circular* analytic backend already exists**:
  `_InclinedCircularBackend` rotates the in-plane circular state by
  `R_z(+lan) @ R_x(-inc)` (`core/ephemeris.py:114-169`). It is wired into
  `_CircularBackend.state` behind an **exact** `inc_deg != 0.0` test
  (`core/ephemeris.py:106-111`) so the live default stays byte-identical.
- **Sourced J2000 inclinations / nodes are already in the table** as record
  fields (`PlanetData.inc_deg/lan_deg/ecc`, `core/constants.py:118-146`) — Venus
  `inc=3.39467605°, lan=76.67984255°`; Mars `inc=1.84969142°, lan=49.55953891°`
  (Standish & Williams Table 1, quoted at `core/constants.py:181-186,209-212`).
  **They are deliberately set to `0.0` in the live `PLANETS` dict** so the
  circular backend stays coplanar; the real values are exercised only via
  injected `PlanetData` into `_InclinedCircularBackend`.
- A **sourced 3-D Tisserand predicate already exists**: `linkable_3d`
  (`search/tisserand.py:471-680`), built on the `cos(i)` Tisserand term
  (`search/tisserand.py:6-8`), with a coplanar short-circuit and a
  monotonicity guarantee (coplanar-True ⇒ 3-D-True). `tisserand_feasible` already
  threads it behind a `use_3d` flag from the presence of an `ephem`
  (`search/sequence.py:300-363`).

So M-3D is **not** a green-field "lift the coplanar assumption" build. It is a
**wiring, activation, and frame-correctness** milestone: turn on the inclined
elements where it is physically correct, fix the **one** place the coplanar
assumption is still load-bearing (the closure frame), and decide the
fast-pre-filter policy. The genuine unsolved physics is narrow and named in §4.

---

## 1. Where the coplanar assumption actually lives (audited, file:line)

Two kinds of "coplanar": (A) **data** that is flattened to z=0, and (B) **frame
algebra** that assumes the orbit plane *is* the ecliptic. Sorted by whether the
lift is real work or a switch.

### A — Data flattening (mostly already liftable)

1. **Live `PLANETS` inc/lan = 0.0** (`core/constants.py:181-186` Venus,
   `:209-212` Mars). The *only* reason the circular backend is coplanar. Sourced
   real values sit in the same comment. **Lift = a policy decision** (§3), not
   new math.
2. **`_CircularBackend.state` exact-zero gate** (`core/ephemeris.py:106-111`).
   Coplanar path is byte-identical *by construction* when `inc_deg == 0.0`; flip
   happens automatically once (1) carries non-zero inc. No code change needed in
   the backend itself — it already delegates to the inclined backend.
3. **`vinf_to_tisserand` / `tisserand_to_vinf`** fix `i=0` (`cos i = 1`)
   (`search/tisserand.py:88-110`). Scalar V∞↔T conversion at a *single* body —
   this is genuinely i-independent at fixed V∞ (the V∞↔T identity has no i term;
   i enters only the (a,e,i) decomposition). **Correct as-is; not a coplanar
   bug.** See §5.
4. **`linkable` / `vinf_contour` / `_contour_roots_in_u`** are coplanar by
   contract (`search/tisserand.py:1-30, 88-95, 223, 331`). **Superseded, not
   fixed:** `linkable_3d` already exists alongside (§5). The coplanar `linkable`
   stays frozen for byte-stable goldens.

### B — Frame algebra (the real crux — see §2)

5. **`closure_residual` uses the uniform `to_rotating`** about ecliptic **+z**
   at a single scalar `omega` (`model/cycler.py:235-252`, transform
   `core/frames.py:90-146`). This is exact **only** for circular-coplanar (the
   docstring says so, `model/cycler.py:233`). With inclined/eccentric states the
   spacecraft has real z-motion that this frame does not co-rotate — closure
   "breathes" out of plane. **This is the one place the coplanar assumption is
   still load-bearing in idealized mode.**
6. **`synodic_omega`** returns a scalar mean-motion about +z
   (`core/frames.py:203-230`). Plane-agnostic in magnitude but assumes the spin
   axis is the ecliptic normal.
7. **The dynamic frame `to_rotating_dynamic`** reads `θ(t)=atan2(r_y,r_x)` and
   rotates about **+z**, preserving z but **not** tilting to the orbit normal
   (`core/frames.py:300-389`). It is z-aware (keeps the z-component) but
   **plane-fixed** (frame normal = ecliptic +z, not body[0]'s orbit normal).
   For real DE440 states this is what M-ED's drift metric runs in
   (`verify/propagate.py:289-310`, `verify/real_closure.py:85,775`). Whether
   that is "good enough" or needs an inclined-frame variant is an **open
   question** (§8).

### Plane-agnostic already (no lift needed — confirmed)

- **`core/flyby.py`** — `max_bend`, `bend_angle`, `flyby_dv`, `dv_from_turn_deficit`
  all operate on full 3-D V∞ vectors via `np.dot`/`np.linalg.norm`
  (`core/flyby.py:117-141, 182-231`). The bend is the unsigned 3-D angle between
  two vectors; the powered-flyby surrogate is already plane-free. **The
  out-of-plane steering DOF the spec wants (§12.1) is geometrically present in
  the V∞ vectors — what is missing is decomposition/reporting, not the cost
  model.** See §6.
- **`construct_cycler`** (`search/construct.py:124-204`) — carries 3-D states
  through Lambert into `Encounter.r`, `vinf_in/out`, `Leg.v_depart/arrive`. No
  flattening.
- **`optimise_cell_idealized`** (`search/optimize.py:960-1136`) — optimises
  *interior epochs only*; geometry is whatever `construct`/`ephem` produce. The
  objective is `closure_residual + Σ flyby_dv` — so it inherits the coplanar
  assumption **transitively through `closure_residual` (B-5)**, not directly.
  No 2-D shortcut in the optimiser itself.

**Condensed audit (the one-screen version):**

| # | Location | Coplanar? | Lift cost |
|---|---|---|---|
| 1 | `constants.py:181,209` live inc/lan=0 | data flattened | policy switch (§3) |
| 2 | `ephemeris.py:106-111` zero-gate | auto-lifts via (1) | none |
| 3 | `tisserand.py:88-110` V∞↔T | i-free at fixed V∞ | correct as-is |
| 4 | `tisserand.py` `linkable`/contour | by contract | superseded by `linkable_3d` |
| 5 | **`cycler.py:235-252` closure frame** | **load-bearing** | **real work (§2)** |
| 6 | `frames.py:203-230` `synodic_omega` | +z spin axis | follows from §2 |
| 7 | `frames.py:300-389` dynamic frame | z-aware, plane-fixed | open Q (§8) |

---

## 2. The crux: the closure frame in 3-D

Idealized closure (`model/cycler.py:closure_residual`) asks: after one period,
does the spacecraft's departure velocity, viewed in a synodic rotating frame,
return to its initial value? Today the frame spins about the **ecliptic +z** at
a constant scalar `omega` (`to_rotating`). For a coplanar circular system the
orbit plane *is* the ecliptic, so this is exact. For inclined orbits the
spacecraft trajectory has genuine z-excursion that a +z-spinning frame does not
absorb, so a geometrically-perfect inclined cycler shows a spurious non-zero
`closure_residual`.

Three valid formulations of the inclined closure frame (this is the decision the
design exists to make):

- **F1 — Keep ecliptic +z, tolerate the z-residual.** Change nothing in
  `frames.py`; re-baseline `closure_residual` thresholds and treat the
  out-of-plane component as "geometric breathing" exactly as spec §12(c)
  prescribes for the *ephemeris* frame (`docs/spec.md:233`). Cheapest; honest;
  but conflates idealized-mode (which is supposed to close *exactly*) with
  ephemeris-mode tolerance.
- **F2 — Spin about the system invariable-plane normal.** Define the frame
  z-axis as the (mass-or-angular-momentum-weighted) normal of the participating
  bodies' orbit planes, spin about *that*. A single rotation `R` maps ecliptic →
  invariable plane; apply it before/after the existing `to_rotating`. Exact for a
  *common* inclined plane; approximate when V and M planes differ (they do:
  3.39° vs 1.85°, different nodes).
- **F3 — Per-body orbit-normal frame, generalised `omega` to a vector.** Replace
  the scalar `omega ẑ` with an angular-velocity **vector** `ω⃗` along the anchor
  body's instantaneous orbit normal; the Coriolis term `v − ω⃗ × r` becomes a
  full 3-D cross product (today it hardcodes the z-only form,
  `core/frames.py:132-136`). Most physically faithful; matches what the dynamic
  frame should arguably also do (§8). Most code change.

**Recommendation: F3 for the transform primitive, with F1 as the honest interim
acceptance criterion.** Generalise `to_rotating`/`from_rotating` to a vector
`ω⃗` (a strict superset — `ω⃗ = ω ẑ` reproduces today bit-for-bit, preserving
every M3 golden), and anchor `ω⃗` on the home body's orbit normal computed from
its `(r,v)` (`(r×v)/|r|²`, which `synodic_omega_dynamic` already computes the
z-component of, `core/frames.py:286-297`). This is the smallest change that makes
inclined idealized closure *mean* something, and it unifies the uniform and
dynamic frames under one vector-ω primitive. Until the F3 frame is validated,
report the inclined `closure_residual` against an F1-style tolerance so we never
claim exact closure we cannot prove.

---

## 3. The inclined-idealized middle rung

The fidelity ladder (`verify/fidelity.py`, `data/provenance.py:95-103`) already
names **three** tiers: `circular-coplanar` → `analytic-ephemeris` → `real-de440`.
The middle rung `analytic-ephemeris` is currently **unavailable** —
`solve_at_fidelity` raises `FidelityRungUnavailableError` for it
(`verify/fidelity.py:62-68,157`, "no in-house backend"). The inclined backend
(`_InclinedCircularBackend`) is the natural inhabitant of a rung between flat
circular and full DE440: **circular but inclined** (real i, Ω, mean sma;
eccentricity optionally on via `PlanetData.ecc`).

Two ways to slot it in:

- **R-A — A new fourth tier `circular-inclined`** between `circular-coplanar`
  and `analytic-ephemeris`. Cleanest semantically (the inclined-circular model is
  *not* a mean-element analytic ephemeris — it has no eccentric anomaly solve),
  but it widens the `Fidelity` literal, which is wired into provenance,
  cross-fidelity persistence, and the catalogue schema. Larger blast radius.
- **R-B — Make `analytic-ephemeris` resolvable via the inclined backend**
  (inc+ecc on). Reuses the existing tier; honest *if* we also turn on `ecc`
  (otherwise "analytic-ephemeris" with zero eccentricity is mislabelled). The
  inclined backend is circular-only today; adding eccentricity is a Kepler solve
  on top of the existing rotation (the machinery exists in `core/kepler.py`).

**Recommendation: R-A (`circular-inclined` as its own rung), inc-only first,
ecc as a follow-on.** Rationale: it keeps each rung a *single* honest physical
assumption (the whole point of the fidelity ladder per `verify/fidelity.py:1-8`,
the S1L1 5.65-vs-4.99 confusion class), and the inclined-circular backend is
*already built and sourced* — we are labelling a capability we have, not faking
one. Eccentricity is a separable later lift that promotes `circular-inclined` →
`analytic-ephemeris` cleanly. **Interaction with the Phase-1 fidelity ladder
(being built now):** M-3D should *register* the rung and its `solve_at_fidelity`
dispatch, but **not** redefine the ladder's persistence semantics — that is the
Forge's. M-3D supplies a rung; the Forge decides what "persists across rungs"
means. Keep the change additive.

**Trade-off vs jumping straight to DE440 (M-ED's territory):** the inclined
rung is a *pre-filter and a diagnostic*, not a discovery engine. It answers "does
this geometry survive real inclination at all?" cheaply (closed-form, no Lambert
multi-start, no launch-epoch search) before M-ED pays for a DE440 phase-match.
That is exactly the roadmap's "idealized circular model may stay 2D as a fast
pre-filter" decision (`...roadmap.md:18,78`) — except the *inclined-circular*
rung is a strictly better pre-filter than the flat one for any row whose
feasibility is inclination-sensitive (chiefly VEM, where V's 3.39° dominates).

---

## 4. The genuine unsolved physics (narrow, named)

After §1's audit, the *new* physics M-3D must actually get right is small:

1. **The inclined closure frame** (§2) — choosing F1/F2/F3 and validating that
   `ω⃗ = ω ẑ` reproduces every M3 golden bit-for-bit.
2. **Plane-change cost attribution at flybys** (§6) — decomposing the (already
   3-D) V∞ bend into in-plane + out-of-plane so the optimiser can *exploit* the
   free plane-change the spec promises, and so we can report it.
3. **Node/argument phasing for inclined construction seeds.** The inclined
   backend places each body on its ascending node at `t_sec=0`
   (`core/ephemeris.py:117-123`). When two bodies have different `lan`, their
   `t=0` phases are no longer the simple `θ=0` the flat model assumes — the
   construction seed (`search/optimize.py` `_free_return_seed`) may need the node
   offset folded in, or the optimiser's epoch search must cover it. Likely
   absorbed by the existing multi-start + DE pass, but must be checked, not
   assumed.

Everything else is wiring.

---

## 5. Tisserand in 3-D — what changes, honest limits

The 3-D predicate is **already implemented and sourced** (`linkable_3d`,
`search/tisserand.py:471-680`; Strange & Longuski 2002, cited
`search/tisserand.py:500`). The mechanism: at fixed V∞, each body's Tisserand
equation `T_p = a_p/a + 2 cos(i) √((a/a_p)(1−e²))` analytically fixes
`cos(i_sc)` at every `(a,e)` (`search/tisserand.py:452-468`); a 2-D `(a,e)` scan
asks whether the two bodies agree on a single reachable `cos(i_sc)` within
tolerance, gated on physical orbit-crossing and `i_sc ≤ i_sc_max`
(`search/tisserand.py:471-499`). **What "changes" going to 3-D is the
intersection geometry:** coplanar `linkable` intersects two 1-D contours in the
`(a,e)` plane; `linkable_3d` intersects two 2-D surfaces in `(a,e,i)` projected
back to a 1-D agreement locus — strictly *more* pairs link (inclination opens
options), and the monotonicity guarantee (coplanar-True ⇒ 3-D-True,
`search/tisserand.py:544-559`) makes it a conservative superset.

**Honest limits (must be stated, not hidden):**

- The V∞↔T scalar identity (`vinf_to_tisserand`) is **genuinely i-independent**
  (`T = 3 − V∞² a_p/μ`). It is *not* a coplanar bug; i enters only the (a,e,i)
  decomposition, which `linkable_3d` handles. Do not "fix" `vinf_to_tisserand`.
- `linkable_3d` tests *energetic linkability with an inclination budget*, not
  *node compatibility*. Two bodies can share a reachable `i_sc` magnitude while
  their orbit planes' **lines of nodes** make the transfer phasing-infeasible.
  The predicate is a necessary, not sufficient, screen — same character as
  coplanar `linkable`. This limit must be documented at the gate, not papered
  over.
- `i_sc_max_deg` default 30° (`search/tisserand.py:476`) is a *modelling
  choice*, not sourced. For V/E/M cyclers (i ≤ 3.4°) it is generous; flag it as
  a tunable, not a physical constant.
- M-3D's Tisserand work is therefore **mostly wiring + honest-limits
  documentation**, not new solver math. The lift is to (a) make
  `tisserand_feasible` consult `linkable_3d` for the inclined rung consistently,
  and (b) pin the behaviour as reviewed diffs (§5 test strategy in §6 below).

**Explicit scope fence vs the T-P graph (adjacent, distinct):** the
Tisserand–Poincaré graph (`docs/notes/2026-06-05-endgame-tisserand-mining.md`,
Campagnola & Russell 2010) generalises Tisserand to **CR3BP, planet-centric**
moon systems (T can exceed 3; ballistic transfers linked-conics calls
impossible). That is **task #76 / Forge Pt-2 territory and explicitly NOT
M-3D.** M-3D stays heliocentric linked-conic; it adds *inclination* to the
existing Tisserand graph, not *CR3BP reachability*. Keep the two from bleeding
together: the 3-D term here is `cos(i)` in a linked-conic T; the T-P graph's
extension is a different feasibility region entirely.

---

## 6. Plane-change steering at Venus, and the test strategy

### Plane-change steering (spec §12.1)

The spec claims VEM closure lives in the cheap out-of-plane steering a 3-D
gravity assist provides "for free" — each flyby has two steering DOF (turn angle
+ out-of-plane node) and "pays most of the plane-change for free"
(`docs/spec.md:252`). The good news from §1: **`core/flyby.py` already operates
on full 3-D V∞ vectors** — `bend_angle` is the unsigned 3-D angle, `flyby_dv`
the 3-D magnitude+bend-deficit surrogate. The out-of-plane bend is *already
costed correctly*; the flyby model does not need a coplanar lift.

What is **missing** is *attribution and exploitation*:

- **Decomposition (reporting):** split the V∞-in→V∞-out rotation into the
  in-plane component (changes a,e) and the out-of-plane component (changes i, the
  node). Both are already inside `bend_angle`'s single number; the optimiser and
  the catalogue should *see* the split so a reviewer can confirm Venus is doing
  the plane-change work the spec predicts. New helper, e.g. `bend_decompose(v_in,
  v_out, orbit_normal) -> (delta_inplane, delta_outofplane)`; pure function,
  alongside `flyby.py`'s existing primitives.
- **Exploitation (optimiser):** the idealized optimiser today varies only
  interior epochs (`search/optimize.py:1048`). To *use* the free plane-change it
  must be allowed to vary the flyby b-plane orientation (the out-of-plane node
  DOF) at the inclined-rung body — a new continuous DOF, gated to the inclined
  rung. This is the one optimiser-surface change M-3D introduces, and it is
  what could make a coplanar-infeasible VEM row close. **Decision needed**: add
  the b-plane DOF in M-3D, or defer the optimiser change to M-ED and have M-3D
  ship only the decomposition/diagnostic? (Open Q, §8.)

### Test strategy (sourced anchors are almost absent — be honest)

The pinned coplanar baselines become **reviewed diffs**, exactly as the M8-Core
tests anticipate: `tests/search/test_sequence_multibody.py:1-7` *explicitly*
pins the coplanar `tisserand_feasible` baseline "so the M-3D inclination lift is
a reviewed diff, not a silent behaviour change." M-3D's first test obligation is
to flip those with documented physical reasons, not to invent new green tests.

Sourced 3-D anchors we can use:

- **The inclinations/nodes themselves** are sourced (Standish & Williams,
  `core/constants.py:181-186,209-212`) — a 3-D-state assertion that the inclined
  backend reproduces V at i=3.39°, M at i=1.85° at the node is a legitimate
  *sourced* golden (the input is the source).
- **DE440 z-components** are ground truth for the astropy backend (already
  asserted indirectly). The inclined-circular backend can be cross-checked
  against DE440 *near the node* for order-of-magnitude agreement (a
  consistency check, NOT a golden — DE440 is eccentric/perturbed; only the sign
  and ~i-scale of the z-excursion should match).
- **Physics invariants, again** (per project memory `golden_tests_sourced_only`):
  the honest backbone of M-3D's tests. (a) `linkable_3d` monotonicity
  (coplanar-True ⇒ 3-D-True) — already a gate (`search/tisserand.py:544-559`).
  (b) `to_rotating` with `ω⃗ = ω ẑ` reproduces the scalar form bit-for-bit
  (regression invariant). (c) Round-trip identity of any new inclined frame
  transform (`from(to(x)) == x`), mirroring the dynamic-frame gate
  (`core/frames.py:399-407`). (d) Tisserand `cos(i)` term reduces to the
  coplanar value at i=0.
- **What we do NOT have:** published 3-D V∞/ΔV closure numbers for any V/E/M
  cycler at the inclined-circular rung. There is **no sourced "inclined cycler
  closes at X km/s" anchor.** M-3D must not manufacture one; closure values it
  computes are *our* numbers and can only be used as regression pins (reviewed
  diffs), never as golden EXPECTED values (project memory
  `golden_tests_sourced_only`).

---

## 7. Recommended approach (the 2–3 alternatives + pick)

**Approach 1 — "Activate + frame-fix" (RECOMMENDED).** Turn on the existing
inclined substrate behind an explicit opt-in (do **not** mutate live `PLANETS`);
fix the closure frame with the F3 vector-`ω⃗` generalisation; register the
`circular-inclined` fidelity rung; consult `linkable_3d` consistently; add the
flyby bend-decomposition *diagnostic*; **defer the optimiser b-plane DOF to a
gated follow-on**. Smallest correct change; every existing golden stays
bit-for-bit because the lift is opt-in and the frame generalisation is a strict
superset. Ships the honest rung and the diagnostic; does not over-claim VEM
closure.

**Approach 2 — "Full inclined optimiser".** Approach 1 plus the b-plane
out-of-plane DOF in `optimise_cell_idealized`, attempting actual inclined VEM
closure in M-3D. Higher value (could close a coplanar-infeasible row) but
higher risk: it overlaps M-ED's optimiser work, and a failed close would be hard
to attribute (frame? DOF? phasing?). Better as a *measured* follow-on once the
rung + frame are validated.

**Approach 3 — "Mutate the model to inclined globally".** Set live `PLANETS`
inc/lan to the sourced values and let everything go 3-D. Rejected: breaks every
coplanar golden silently, violates the byte-stability discipline the codebase
has carefully preserved (`core/constants.py:181-186` comment, `ephemeris.py:99`),
and conflates the rungs the fidelity ladder exists to separate.

**Pick: Approach 1**, with Approach 2's b-plane DOF scoped as an explicit
follow-on decision (§8 Q4).

### Phased sketch (no task plan yet — order only)

1. **Frame primitive** — generalise `to_rotating`/`from_rotating` to vector
   `ω⃗`; prove `ω⃗ = ω ẑ` reproduces M3 goldens bit-for-bit; add round-trip
   gate. (Crux; foundational.)
2. **Opt-in inclined ephemeris** — a sanctioned way to get inclined states
   (inclined backend with sourced `PlanetData`, **not** by mutating `PLANETS`);
   3-D-state sourced assertions.
3. **Closure in 3-D** — `closure_residual` uses the anchor body's orbit-normal
   `ω⃗`; re-baseline thresholds with documented F1 interim tolerance.
4. **Fidelity rung** — register `circular-inclined`; wire `solve_at_fidelity`
   dispatch additively; do not touch persistence semantics.
5. **Tisserand consistency** — `tisserand_feasible` consults `linkable_3d` for
   the inclined rung; flip the M8-Core pinned baselines as reviewed diffs.
6. **Flyby bend-decomposition diagnostic** — pure helper + reporting; no
   optimiser change.
7. **(Gated follow-on, pending Q4)** — optimiser b-plane DOF.

### Explicit non-goals

- Planet-centric / CR3BP / moon tours / T-P graph — **task #76**, see §5 fence.
- Real-ephemeris blind discovery, `optimise_cell_ephemeris`, TCM budget —
  **M-ED**.
- Multi-rev Lambert solver math — **M-L** (M-3D consumes whatever Lambert
  returns; it is already 3-D-native).
- Eccentricity in the analytic rung — separable follow-on (promotes
  `circular-inclined` → `analytic-ephemeris`); not required for the i-lift.
- Mutating live `PLANETS` to inclined — see Approach 3 rejection.

---

## 8. Open questions for the user

1. **Closure frame (§2):** F1 (tolerate z-residual, re-baseline), F2 (invariable
   plane), or F3 (vector-`ω⃗` orbit-normal)? Recommendation is F3 primitive +
   F1 interim acceptance — confirm, or pick F1/F2 if you want minimal change.
2. **Fidelity rung (§3):** add a new `circular-inclined` tier (R-A, recommended)
   or repurpose the existing `analytic-ephemeris` rung (R-B)? R-A widens the
   `Fidelity` literal and touches provenance/schema — acceptable?
3. **Eccentricity:** inclination-only for M-3D (recommended), or fold
   eccentricity in at the same time (turns the rung into a true
   `analytic-ephemeris`)?
4. **Optimiser b-plane DOF (§6):** ship the out-of-plane steering DOF in M-3D
   (Approach 2, attempt real inclined VEM closure), or M-3D ships only the
   bend-decomposition *diagnostic* and the closure attempt waits for M-ED
   (Approach 1, recommended)?
5. **Dynamic frame (§1 item 7):** does the real-ephemeris drift metric
   (`to_rotating_dynamic`, used by M-ED) need an inclined (orbit-normal) variant
   under M-3D, or is the z-aware-but-plane-fixed dynamic frame left for M-ED to
   own with its own tolerance? (Scope-boundary call between M-3D and M-ED.)
6. **`i_sc_max_deg` (§5):** keep the 30° default as the inclined-rung tunable, or
   set it tighter (≈5°) for the V/E/M regime to make the screen sharper?

---

## Approval (2026-06-05)

User-approved with all recommendations accepted: (Q1) F3 vector-ω⃗
orbit-normal closure-frame primitive with F1 interim acceptance; (Q2) R-A —
new `circular-inclined` fidelity rung (Fidelity literal/schema widening
accepted); (Q3) inclination-only for M-3D, eccentricity as a separable
follow-on; (Q4) Approach 1 — bend-decomposition diagnostic only, inclined
closure attempts stay with M-ED; (Q5) the dynamic-frame inclined variant is
M-ED's to own; (Q6) keep the 30° default as the rung tunable. M-3D is NOT on
M-ED's critical path (per the M-ED design §6.1 finding); plan written when
scheduled.
