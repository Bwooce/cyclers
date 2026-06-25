# Design draft — #305 BCR4BP V0–V5 validation gauntlet adaptation

Status: **DESIGN-DRAFT (scoping, not implementation).** User reviews before any build.
Date: 2026-06-25
Author scope: maps the V0–V5 gauntlet onto BCR4BP periodic orbits, reusing the existing
BCR4BP genome (#292/#303/#304/#334) and the QP-tori gauntlet (#319) as the adaptation
template. No code is written by this draft; no catalogue writeback is proposed.

---

## 0. TL;DR

- BCR4BP candidates are **libration-point families** (planar L1 Lyapunov #303, 3D halos
  #304, POL1 substitute), **not flyby cyclers**. The gauntlet here validates *periodic
  orbits in a Sun-perturbed Earth–Moon model*, not patched-conic itineraries.
- **BUILT today:** only the in-corrector dual-closure check (DOP853 corrector residual +
  independent Radau re-propagation) — this is effectively a **V0 + a partial V1**, living
  inside `correct_bcr4bp_periodic`. There is **no `data/validation/v*_bcr4bp.py` lane**, no
  driver scripts, no frozen tests, no `_LEVEL_EVIDENCE` registration. **V0..V5 are all
  MISSING as standalone gauntlet tiers.**
- **Recommended first build: V1-BCR4BP** (`v1_bcr4bp.py`), mirroring `v1_qp.py` — it
  promotes the corrector's existing independent Radau check into a standalone frozen
  verdict with a sourced floor. Smallest, highest-confidence, unblocks everything above it.
- **Key BCR4BP-specific risk (the load-bearing one):** the model has **two angular
  frequencies** — the Earth–Moon synodic rate (= 1, nondim) and the Sun synodic rate
  `ω_S ≈ 0.925196 rad/TU`. Their ratio is irrational, so a BCR4BP orbit is **strictly
  periodic only when its period is Sun-commensurate** (`T = 2π·n/ω_S`, integer `n`); a
  generic family member is **quasi-periodic** with a residual `sun_phase_drift = |ω_S·T −
  2π·n|`. This directly governs how V2 (multi-lap periodicity) must be defined — see §5.

---

## 1. What exists (BUILT vs MISSING), with sources

### 1a. The BCR4BP genome (read; READ-ONLY for this gauntlet)

| Component | File | Role |
|---|---|---|
| Core model (EOM, STM, propagator, constants) | `src/cyclerfinder/core/bcr4bp.py` | `bcr4bp_eom`, `bcr4bp_stm_eom`, `propagate_bcr4bp`, `BCR4BPSystem`, `andreu_default()`, `sun_commensurate_period(omega_sun, n)` |
| Periodic-orbit corrector | `src/cyclerfinder/genome/bcr4bp_genome.py` | `correct_bcr4bp_periodic(...) -> BCR4BPPeriodicOrbit` |
| μ_sun continuation | `src/cyclerfinder/genome/bcr4bp_continuation.py` | `continue_bcr4bp_family_in_musun(...) -> BCR4BPFamily` |
| Multi-system registry | `src/cyclerfinder/genome/bcr4bp_systems.py` | `REGISTRY`, `build_bcr4bp_system(...)`, `derive_from_sources(...)` |

Model facts that matter to the gauntlet (from `core/bcr4bp.py`, sourced Rosales-Jorba 2023
Table 3 / Gimeno-Jorba 2018):
- Frame: Earth–Moon synodic rotating, nondim. Synodic rate normalized to 1.
- Sun on a planar circle at `θ_S = θ_S0 + ω_S·t`, `ω_S = 0.925195985520347 rad/TU`,
  `T_S = 2π/ω_S ≈ 6.7912 TU ≈ 29.5 d` (lunar synodic month).
- Implements the **incoherent** standard BCR4BP (Simó/Jorba/Gómez four-constant model),
  **not** Andreu's coherent QBCP. Reduces *exactly* to CR3BP at `mu_sun = 0`.
- Non-autonomous: `t0` (equivalently `θ_S0`) is a real parameter — closure is **epoch
  (Sun-phase) dependent**.

The corrector's existing convergence predicate (this is the seed of V0/V1):
```
converged = (corrector_residual < tol)            # DOP853 single-shooting Newton
            AND (independent_closure_residual < independent_tol)   # Radau re-propagation
```
plus it already computes `sun_phase_drift = |ω_S·T − 2π·n|`. Result dataclass
`BCR4BPPeriodicOrbit` carries: `state_initial, period_nondim, sun_commensurate_n,
sun_phase_drift, converged, corrector_residual, independent_closure_residual, n_iter,
system, free_vars, residual_indices, is_half_period_residual, notes`.

Families that exist as driver scripts (candidate sources for the gauntlet):
- `scripts/run_303_bcr4bp_l1_continuation.py` → `data/bcr4bp_l1_family_303.jsonl` (planar L1 Lyapunov)
- `scripts/run_304_bcr4bp_halo_continuation.py` → `data/bcr4bp_halo_family_304.jsonl` (L1 southern halo)
- POL1 closes as a corrector seed in `tests/genome/test_bcr4bp_genome.py` (to a *nearby*
  incoherent-BCR4BP orbit, not the published QBCP numbers). **POL2 is NOT implemented**
  (a declared #292 negative — did not converge with the symmetric pattern).

### 1b. The V0–V5 spec (canonical) — `docs/spec.md` §14

Verbatim/near-verbatim tier definitions (lines 388–511). These are written for
**patched-conic / CR3BP cyclers**, so each must be *reinterpreted* for a BCR4BP periodic
orbit (the same move #319 made for QP-tori):

- **V0 — Internal consistency:** hard constraints met; V∞ preserved across each flyby;
  closure residual ≤ tol (idealized).
- **V1 — Solver cross-check:** every leg re-solved with an independent solver, agreement
  < 1e-3 m/s; full trajectory re-propagated with an integrator *other than* the one that
  built it; positions met < tol.
- **V2-ballistic — Multi-lap periodicity:** ≥3 continuous laps; bounded drift in the
  defining rotating frame (tolerant of geometric breathing), evaluated in the row's
  defining model.
- **V3 — Ephemeris realisation / independent integrator:** confirmed on an integrator
  independent of the finding solver; bounded over 3–5 laps.
- **V4 — High-fidelity external:** independent codebase + ephemeris (NASA GMAT, or
  Tudat/pykep n-body) reproduces trajectory within tol.
- **V5 — Novelty + expert review:** canonical signature misses catalogue + literature;
  human review.

Trust gate (§14, verbatim): "only V3+ candidates are 'credible'; only V5 + catalogue/
literature miss may be called a discovery." **#391 Hill pre-screen caveat:** CR3BP spectral
stability does NOT imply real-eph survival — a CR3BP-stable orbit can clear V1–V3 and still
fail a Sun-on model by 4–5 orders. *This is the entire motivation for a Sun-perturbed
intermediate rung; BCR4BP is structurally the model that closes that gap below DE440.*

### 1c. The QP-tori gauntlet (#319) — the adaptation TEMPLATE

Files under `src/cyclerfinder/data/validation/`: `v1_qp.py`, `v2_qp.py` (QP lane reaches
only V1+V2 today; V3+ deferred). The reusable pattern:

- **Standalone `run_vN_*()` per tier** → frozen `VNVerdict*` dataclass with a `passes_vN`
  headline bool. No central dispatcher; **tiers chain by argument passing** (`run_v3`
  takes `v2_verdict=`, etc.). Each tier asserts `len(prior.per_cycle) >= n_cycles`.
- Driver script per task → writes `data/*_verdict.jsonl` → **frozen pytest** asserts the
  JSONL → the JSONL string is registered in `validate.py::_LEVEL_EVIDENCE[(id, level)]`,
  which is the *only* promotion gate (a census over-claim guard, line ~502+). No gauntlet
  module ever writes the catalogue.
- Sourced-or-labelled floor constants as module-level `Final`.
- **The QP reinterpretation diff** (the precedent we copy): strict `||X(T)−X(0)||`
  closure does not exist on a torus, so V1 became **invariance modulo rotation** (Fourier
  residual + an independent **off-grid sample** check with a *different RNG seed* —
  statistical independence replacing integrator independence); V2 became "does the state
  still lie on the torus after k laps" with an empirically-calibrated nondim floor.

This QP precedent is the direct structural analogue for BCR4BP, because **a BCR4BP orbit
is also quasi-periodic in general** (two frequencies) — see §5/§7.

### 1d. The V4 real-eph infrastructure

- `v4_uranus_strict.py` (#335) is the only fully-built real-ephemeris V4 lane: **Python +
  `spiceypy` → scipy DOP853**, kernels furnished from the GMAT R2022a install
  (`~/GMAT/R2022a/data/.../*.bsp|.tls|.tpc`). spiceypy binds the same JPL SPICE C library
  GMAT compiles against. It is **moon-tour / planet-frame**, not heliocentric.
- GMAT R2022a confirmed present at `~/GMAT/R2022a`; headless `env -u DISPLAY ./GmatConsole
  --run <script>`. GMAT loads **DE405**; the in-house independent-n-body lane uses DE440.
- **There is no Sun–Earth–Moon 4-body real-eph V4 module yet.** This is the largest gap.

### 1e. BUILT-vs-MISSING map (the headline)

| Tier | BCR4BP meaning (this draft) | Status | Where the partial lives |
|---|---|---|---|
| **V0** | Internal consistency: corrector residual ≤ tol, state finite, r_p ≥ r_p_min vs Earth/Moon, `sun_commensurate_n` integer + `sun_phase_drift` recorded | **PARTIAL** (logic exists inside corrector; not a standalone tier, no floors gated, no constraint checks) | `correct_bcr4bp_periodic` |
| **V1** | Same-model closure under an **independent integrator** (Radau vs DOP853), at the masked residual indices, < floor | **PARTIAL** (Radau check exists in-corrector; not promoted to a frozen verdict, no sourced floor, no km/s gate) | `correct_bcr4bp_periodic.independent_closure_residual` |
| **V2** | Long-span **bounded drift** in the coherent 4-body model over ≥3 Sun-commensurate laps (quasi-periodic — see §5) | **MISSING** | — |
| **V3** | Independent integrator over the long span (REBOUND IAS15 vs scipy), same BCR4BP model | **MISSING** | — |
| **V4** | Real-ephemeris **Sun–Earth–Moon** DE440/SPICE realisation | **MISSING** (no heliocentric/SEM 4-body lane exists; `v4_uranus_strict` is planet-frame) | — |
| **V5** | Novelty + literature miss + human review | **MISSING** (and note #303/#304/#334 found no new species — V5 has nothing to promote yet) | — |

---

## 2. Scope decision (read before architecting)

**The existing structural findings strongly shape what this gauntlet is *for*.**
#303/#304/#326/#334 are a chain of clean negatives: BCR4BP families continue *smoothly*
from CR3BP with **no new species, no bifurcations**, and the Sun perturbation is material
only at Sun–Earth–Moon (`μ_sun/a_sun³ ≳ 1e-4`); it is machine-precision-negligible at
Sun–Jupiter and most triples (k ≈ 2.89 scaling, #334). So:

- The gauntlet's near-term value is **NOT** "validate a novel BCR4BP cycler for discovery"
  (none exists). Its value is twofold:
  1. **A deterministic, paper-reproducible regression rung** ("did the Sun term stay wired
     in correctly?") that chaotic V3/V4 cannot bisect — exactly the Andreu-digest verdict.
  2. **The missing intermediate fidelity rung for the #391 Hill problem**: a CR3BP-stable
     cycler that the Hill pre-screen flags MARGINAL can be *run through the SEM BCR4BP V2*
     to see whether the Sun actually unbinds it — cheaper and more diagnosable than DE440.
- **Recommendation: build V1→V2→V3 against BCR4BP periodic orbits first; treat V4 (SEM
  4-body real-eph) as a separately-scoped sub-project; do not build V5 until a candidate
  with novelty potential exists.** This matches the QP lane's own choice (it stopped at V2).

---

## 3. Mapping V0..V5 to BCR4BP-periodic meaning

| Tier | Generic spec | BCR4BP reinterpretation |
|---|---|---|
| V0 | hard constraints + idealized closure | corrector converged at `tol`; finite state; `r_p` to Earth and Moon ≥ floors; `sun_commensurate_n` is a positive integer and `sun_phase_drift` recorded (the periodicity is *conditional* on commensurability — V0 must surface this, not hide it) |
| V1 | re-solve legs + re-propagate with a different integrator | **same-model closure under an independent integrator**: re-propagate the converged IC under Radau (≠ the DOP853 the corrector used) over the masked period; closure residual < floor, both nondim and km-converted |
| V2 | ≥3 continuous laps, bounded drift in defining model | ≥3 **Sun-commensurate** laps (`k·T`, `k=1..n_cycles`) propagated in the coherent BCR4BP at fixed `θ_S0`; bounded position drift `||X(kT) − X(0)||` ≤ floor — **but only meaningful when the orbit is Sun-commensurate**; for non-commensurate members V2 measures quasi-periodic recurrence, see §5 |
| V3 | independent integrator, long span | re-run the V2 long span under a *different integrator family* (REBOUND IAS15 / scipy LSODA vs DOP853); per-lap terminal positions agree with V2 within floor → bounded-drift signature is integrator-independent |
| V4 | independent codebase + real ephemeris | **Sun–Earth–Moon DE440/SPICE**: seed the BCR4BP IC into a real-eph SEM integration (spiceypy Sun+Earth+Moon point masses, true lunar eccentricity/inclination, epoch from `θ_S0`); bounded over the span; cross-checked vs an independent codebase (GMAT script emitter, DE405, for human confirmation) |
| V5 | novelty + literature + human | canonical signature misses catalogue + BCR4BP literature anchors (`search/literature_check.py` already holds #334 anchors); human review |

---

## 4. Architecture per level (module/function signatures + assertion + floor)

New package files mirror the QP lane: `src/cyclerfinder/data/validation/v{0,1,2,3,4}_bcr4bp.py`,
exported in `__init__.py`. All take a converged `BCR4BPPeriodicOrbit` (re-verify, never
re-solve). All return a frozen `V*VerdictBCR4BP` with a `passes_vN` headline bool. All are
READ-ONLY on the genome. None writes the catalogue.

### V0 — `v0_bcr4bp.py`
```python
def run_v0_bcr4bp(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    *,
    corrector_floor: float = V0_BCR4BP_CORRECTOR_FLOOR,       # 1e-10 (matches corrector tol)
    rp_floor_earth_km: float = ...,                           # sourced per-body flyby floor
    rp_floor_moon_km: float = ...,                            # flyby_altitude_references.yaml
    phase_drift_floor: float = V0_BCR4BP_PHASE_DRIFT_FLOOR,   # see §5; LABELLED convention
    notes: str = "",
) -> V0VerdictBCR4BP
```
- **Asserts:** `orbit.converged`; `corrector_residual <= corrector_floor`; state vector all
  finite; min periapsis radius to Earth and to Moon over one period ≥ sourced floors;
  `sun_commensurate_n` is a positive integer; `sun_phase_drift <= phase_drift_floor` **or**
  the verdict explicitly tags the orbit `quasi_periodic=True` (does not silently pass).
- **Floor:** `V0_BCR4BP_CORRECTOR_FLOOR = 1e-10` (sourced from corrector `tol`); periapsis
  floors from `flyby_altitude_references.yaml` (per the digest≠adoption discipline); phase
  drift floor is a **labelled convention** (§5), not a physical constant.

### V1 — `v1_bcr4bp.py` (RECOMMENDED FIRST BUILD)
```python
def run_v1_bcr4bp(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    *,
    independent_floor_nondim: float = V1_BCR4BP_FLOOR_NONDIM,  # 1e-6 (matches v1_3d)
    v1_floor_kms: float = V1_FLOOR_KMS,                        # 1e-3 km/s, reuse spec floor
    rtol: float = 1e-12, atol: float = 1e-12,
    notes: str = "",
) -> V1VerdictBCR4BP
```
- **Asserts (mirrors `v1_3d.run_v1_3d`):** re-propagate `orbit.state_initial` under **Radau**
  (independent of the corrector's DOP853) over the period (full or half per
  `is_half_period_residual`); masked closure residual ≤ `independent_floor_nondim` AND its
  km/s conversion ≤ `v1_floor_kms`. Distinct from the corrector's own check by *recomputing
  freshly* and recording in the verdict (the QP-lane "fresh confirmation" discipline).
- **Floor:** `V1_BCR4BP_FLOOR_NONDIM = 1e-6` (matches periodic V1, Olikara-Howell band);
  `V1_FLOOR_KMS = 1e-3` reused verbatim from spec §14 / `v1_3d`.
- **Why first:** the Radau machinery already exists inside the corrector; this tier *promotes
  and freezes* it with a sourced floor and a standalone verdict. Lowest risk, unblocks V2+.

### V2 — `v2_bcr4bp.py`
```python
def run_v2_bcr4bp(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    *,
    n_cycles: int = V2_BCR4BP_N_CYCLES_MIN,                   # 3
    drift_floor_kms: float = V2_BCR4BP_DRIFT_FLOOR_KMS,       # LABELLED, see note
    require_commensurate: bool = True,                        # see §5 / §7 risk
    rtol: float = 1e-12, atol: float = 1e-12,
    notes: str = "",
) -> V2VerdictBCR4BP
```
- **Asserts:** with `θ_S0` fixed at the orbit's epoch, propagate `k·T` for `k=1..n_cycles`
  in the coherent BCR4BP; max `||X(kT) − X(0)||` (km) ≤ `drift_floor_kms`. If
  `require_commensurate` and the orbit is non-commensurate (`sun_phase_drift` above the V0
  convention), **FAIL with reason `non_commensurate_no_strict_period`** rather than
  reporting a meaningless drift — OR (config) fall back to the QP-style recurrence metric
  (§5). Records per-cycle drift list (so V3 can chain on `len(per_cycle) >= n_cycles`).
- **Floor:** `V2_BCR4BP_N_CYCLES_MIN = 3` (spec); `V2_BCR4BP_DRIFT_FLOOR_KMS` is a
  **judgment-call labelled constant** like the QP V2 floor — start from `v2_3d`'s
  `50_000 km` and recalibrate empirically against #303/#304 family members at N≥3.

### V3 — `v3_bcr4bp.py`
```python
def run_v3_bcr4bp(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    *,
    v2_verdict: V2VerdictBCR4BP,
    n_cycles: int = V3_BCR4BP_N_CYCLES_MIN,                   # 3
    agreement_floor_kms: float = V3_BCR4BP_AGREEMENT_FLOOR_KMS,  # ~100 km, mirror v3_3d
    notes: str = "",
) -> V3VerdictBCR4BP
```
- **Asserts:** re-propagate the V2 long span with an **independent integrator family**
  (REBOUND IAS15 with the Sun-acceleration as an `additional_force`, or scipy LSODA);
  per-cycle terminal positions agree with `v2_verdict.per_cycle` within `agreement_floor_kms`
  → the bounded-drift signature is integrator-independent. Chains on
  `len(v2_verdict.per_cycle) >= n_cycles`.
- **RISK (see MEMORY `reference_rebound_variation_custom_force_gotcha`):** REBOUND
  variational particles do NOT auto-differentiate a Python `additional_forces` callback.
  V3 here only needs *state* re-propagation (no STM), so the gotcha is avoided **as long as
  the Sun term is applied as a force, not relied on through variational eqs**. Gate with a
  PERTURBED parity test (Sun-only-correct passes even with the bug — must test *with* the
  perturber active). Floor: `V3_BCR4BP_AGREEMENT_FLOOR_KMS ≈ 100 km` (mirror `v3_3d`).

### V4 — `v4_sem_realeph.py` (separately-scoped; largest risk)
```python
def run_v4_sem_realeph(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    launch_epoch_utc: str,                                    # maps θ_S0 -> real epoch
    *,
    v3_verdict: V3VerdictBCR4BP,
    n_cycles: int = V4_BCR4BP_N_CYCLES_MIN,                   # 3
    agreement_floor_kms: float = V4_BCR4BP_AGREEMENT_FLOOR_KMS,
    spice_kernel_paths=None,                                  # default -> ~/GMAT/R2022a/data
    notes: str = "",
) -> V4SEMVerdict
```
- **Asserts:** convert the BCR4BP nondim IC + `θ_S0` to a J2000 inertial state at
  `launch_epoch_utc` (the epoch chosen so the real Sun–Earth–Moon geometry matches the
  bicircular `θ_S` — this is the hard part, see §7); integrate under **spiceypy** Sun +
  Earth + Moon point masses (real DE440 lunar eccentricity/inclination) over the span;
  bounded within `agreement_floor_kms` vs the V3 BCR4BP span. Optional GMAT script emitter
  (`scripts/gmat_v4_sem_generate.py`, DE405) for independent-codebase human confirmation.
- **This is where the #391 Hill caveat bites:** a BCR4BP orbit may itself fail V4 by orders
  once true lunar eccentricity + solar radiation-free real geometry is applied. **A V4 fail
  is an expected, informative outcome, not a build defect.** Floor mirrors
  `V4_AGREEMENT_FLOOR_KMS = 50_000 km`; recalibrate.
- **Build dependency:** needs a new frame-conversion helper (nondim synodic ↔ J2000
  inertial at epoch) and a SEM epoch-matching routine. Reuse `v4_uranus_strict.py`'s
  spiceypy furnish/cleanup pattern (`_spice_furnsh_all`, `spkezr`, `kclear` in try/finally).

### V5 — not built now
Defer. No BCR4BP candidate with novelty potential exists (#303/#304/#334 negatives). When
one does: reuse `search/literature_check.py` (already holds #334 BCR4BP anchors) + human
review, register in `_LEVEL_EVIDENCE`.

### Orchestration & promotion (identical to QP lane)
- No central dispatcher. Driver scripts per task: `scripts/run_305_v1_bcr4bp_gauntlet.py`
  etc., reading `data/bcr4bp_*_family_*.jsonl` → calling `run_vN_bcr4bp` → writing
  `data/bcr4bp_v{N}_verdict.jsonl`.
- Frozen pytest under `tests/data/test_v{N}_bcr4bp.py` asserts the JSONL.
- Promotion only via `validate.py::_LEVEL_EVIDENCE[(id, "V{N}")]` string registration; no
  gauntlet module writes the catalogue. (And per MEMORY `catalogue_edits_run_all_ratchets`:
  any catalogue row change must run `uv run pytest tests/data tests/search -q`.)

---

## 5. The two-frequency / quasi-periodicity problem (the central risk) — detailed

This is the BCR4BP-specific issue that most shapes the design, and the reason the QP lane
(#319) is the right template rather than the strict-periodic `v1_3d`/`v2_3d` lane.

**The physics.** BCR4BP carries two angular frequencies:
- Earth–Moon synodic rate, normalized to **1** (nondim).
- Sun synodic rate **`ω_S ≈ 0.925196 rad/TU`** (`T_S ≈ 6.7912 TU`).

`ω_S` is irrational. A BCR4BP trajectory closes strictly (`X(T) = X(0)` at the same
Sun-phase) **only if `ω_S·T = 2π·n` for integer `n`** — i.e. the period is
*Sun-commensurate*. The corrector already exposes this as `sun_commensurate_n` and the
residual `sun_phase_drift = |ω_S·T − 2π·n|`. The #303 L1 Lyapunov family (`T ≈ 2.95 TU`)
and #304 halos are **not** generically Sun-commensurate; only the POL substitutes are built
to close at `T = T_S` (`n = 1`).

**Why this breaks a naive V2/V3.** "≥3 continuous laps, bounded drift" assumes a strict
period to lap against. If the orbit is not Sun-commensurate, after one nominal period the
Sun is at a *different phase*, so `X(T) ≠ X(0)` even for a perfect quasi-periodic orbit —
the "drift" you measure is the Sun-phase mismatch, **not** instability. Reporting that as a
V2 drift is a false negative.

**Three options for V2/V3 (recommend Option A as default, B as the general path):**

- **Option A — Commensurate-only V2 (recommended first):** restrict V2/V3 to
  Sun-commensurate orbits (`sun_phase_drift` below the V0 convention). Lap on `k·T`,
  measure true `||X(kT) − X(0)||`. Clean and correct for POL1 and any commensurate family
  member. Non-commensurate members FAIL V2 with reason `non_commensurate_no_strict_period`
  (honest — they have no strict multi-lap period). This is the smallest correct V2.

- **Option B — QP-recurrence V2 (general):** for non-commensurate orbits, adopt the #319
  metric directly: the orbit lives on a 2-torus; measure whether the state, after `k` laps,
  still lies on the torus invariant under the Sun-phase rotation — i.e. compare `X(k·T_S)`
  (lap on the *Sun* period, not the orbit period) and check bounded recurrence. This needs
  a torus representation the BCR4BP genome does not yet emit (it emits single periodic
  orbits, not invariant tori), so Option B is a larger build — possibly a *merge* of #305
  with the QP-tori genome, out of scope for the first pass.

- **Option C — Stroboscopic map V2:** sample the trajectory at the Sun period `T_S`
  (stroboscopic section) and test boundedness of the resulting map orbit. Middle-ground
  build; correct for both commensurate and quasi-periodic cases; recommended as the V2/V3
  upgrade after Option A lands.

**V0/V1 are unaffected** — they validate a single closed period (the corrector already
handles the commensurate-`n` bookkeeping), so the two-frequency issue does not block the
recommended first build. It bites at **V2 onward**, which is why V2's signature above carries
`require_commensurate` and an explicit non-commensurate FAIL path.

**The phase-drift floor is a labelled convention, not a sourced physical constant.** Pick a
threshold (e.g. `sun_phase_drift < 1e-6 rad` ⇒ "commensurate") and label it as a convention
in the module, per the digest≠adoption discipline — do not present it as sourced.

---

## 6. Bite-sized build sequence (TDD, V1 first)

Each step is a red→green→refactor unit with a frozen test, mirroring the QP lane. No
catalogue writeback at any step. Run `uv run ruff check . && ruff format --check .` and
`uv run pytest tests/data -q` before each commit.

1. **V1-BCR4BP (`v1_bcr4bp.py`)** — frozen `V1VerdictBCR4BP`; test: POL1 (commensurate,
   converged) PASSES; a deliberately perturbed IC FAILS the Radau floor. *Highest
   confidence — reuses existing Radau check.*
2. **V0-BCR4BP (`v0_bcr4bp.py`)** — constraints + commensurability surfacing; test: POL1
   passes; an IC with sub-floor Earth periapsis FAILS; a non-integer-`n` orbit is tagged
   `quasi_periodic`. (V0 after V1 is fine — V1 is the load-bearing first cut.)
3. **V2-BCR4BP Option A (`v2_bcr4bp.py`)** — commensurate-only multi-lap drift; test: POL1
   (`n=1`, lap on `T_S`) bounded over 3 laps PASSES; a non-commensurate #303 member FAILS
   with the explicit reason.
4. **V3-BCR4BP (`v3_bcr4bp.py`)** — independent-integrator agreement; test: V2-pass POL1
   agrees within floor under LSODA/IAS15; **include the PERTURBED parity test** for the
   REBOUND custom-force gotcha.
5. **V2/V3 Option C upgrade (stroboscopic)** — optional, after 1–4 land; extends V2/V3 to
   quasi-periodic members. Separate task.
6. **V4-SEM real-eph (`v4_sem_realeph.py`)** — separately scoped; first build the
   nondim-synodic↔J2000 epoch-matching helper with its own golden test (a known
   transformation), *then* the integration. Largest risk; do not start before 1–4 land.
7. **Driver scripts + `_LEVEL_EVIDENCE` registration** — only once a tier is frozen and the
   user approves promotion.

---

## 7. Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| **Two-frequency quasi-periodicity** breaks naive V2/V3 multi-lap (false negatives) | **HIGH** | V0/V1 unaffected (single period). V2 starts commensurate-only (Option A) with explicit non-commensurate FAIL; upgrade to stroboscopic (Option C) later. Phase-drift floor labelled as convention. (§5) |
| V4 SEM 4-body real-eph: **epoch-matching** nondim synodic `θ_S` ↔ real DE440 Sun–Earth–Moon geometry | **HIGH** | Build + golden-test the frame/epoch helper *before* the integrator; reuse `v4_uranus_strict` spiceypy pattern; accept that a BCR4BP orbit may legitimately FAIL V4 (#391 Hill caveat — informative, not a defect). Largest separately-scoped piece. |
| REBOUND **variational/custom-force gotcha** at V3 | MED | V3 needs only state re-prop (no STM) → gotcha avoided if Sun applied as a force; gate with a PERTURBED parity test (MEMORY `reference_rebound_variation_custom_force_gotcha`). |
| **No novel candidate to validate** — #303/#304/#334 are negatives | MED (scope) | Reframe gauntlet value as (a) deterministic regression rung, (b) #391 Hill intermediate fidelity rung; stop at V3 until a candidate exists; do not build V5. (§2) |
| Floor constants masquerading as sourced | LOW | V2 drift + phase-drift floors are explicitly **labelled judgment-call conventions**; V1 floors reuse the sourced spec values (1e-3 km/s, 1e-6 nondim). (digest≠adoption) |
| Incoherent BCR4BP vs Andreu coherent QBCP model gap | LOW | The genome is the *incoherent* model by design; POL1 closes only to a *nearby* orbit. Golden tests target the same model the corrector uses (never the published QBCP numbers) — same-model golden discipline. |

---

## 8. Recommended decision for the user

1. **Approve V1-BCR4BP as the first build** (`v1_bcr4bp.py`, mirroring `v1_qp.py`/`v1_3d.py`)
   — smallest, highest-confidence, promotes the existing Radau check into a frozen tier.
2. **Adopt V2 Option A (commensurate-only)** for the first V2, with the explicit
   non-commensurate FAIL path; defer the stroboscopic (Option C) general path.
3. **Treat V4 SEM real-eph as a separate sub-project**, gated on the epoch-matching helper.
4. **Do not build V5**; reframe the gauntlet's near-term purpose as a regression + #391-Hill
   intermediate rung, not discovery validation (no novel BCR4BP candidate exists).

The single most important BCR4BP-specific risk to internalize: **a BCR4BP orbit is strictly
periodic only when Sun-commensurate; in general it is quasi-periodic with two incommensurate
frequencies, so "≥3 laps, bounded drift" (V2+) must be defined against Sun-commensurability
or a stroboscopic/torus recurrence metric — not a naive strict-period lap.**
