# #306 — Completing the 3D V0–V5 validation gauntlet (DESIGN DRAFT)

**Date:** 2026-06-25. **Status:** DRAFT — scoping only, for review before any
implementation. No code, no commits.

**Motivation (concrete):** the #444 C21 candidate
(`braik_ross_system()`, state0 = `[0.7440212218499672, 0, -0.2057098355650995,
0, 0.35368280201143637, 0]`, T = 18.167169790651315 TU, C = 3.025791614769996,
Floquet **stable**) was just driven through V1 + V2 in the 3D gauntlet end-to-end
and PASSED both. V3 was *noted* as needing a periodic-orbit IAS15 cross-check —
but `v3_3d.py`'s headline `run_v3_3d` is **moontour-specific** (Lambert-leg
IAS15) and does NOT apply to a CR3BP periodic orbit. That is the nearest gap and
the first build target here.

This draft (a) maps every gauntlet level V0..V5 to its meaning for a **3D CR3BP
periodic orbit**, (b) states BUILT vs MISSING per level, (c) gives buildable
architecture + signatures for the MISSING levels, (d) gives a TDD build
sequence, and (e) records risks and open questions.

---

## 0. Scope clarification — what kind of object the 3D gauntlet validates

The 2D/heliocentric gauntlet (`verify/`, `data/validation/v2_moontour.py`,
`v3_3d.run_v3_3d`, `v4_uranus.py`, `real_closure.py`) validates **patched-conic
Lambert-leg cyclers / moontours**: objects defined by a body sequence + leg ToFs
+ V∞ continuity. Its V3/V4 are *Lambert-leg* re-propagations.

The 3D gauntlet (`v1_3d.py`, `v2_3d.py`) validates a **full-3D CR3BP periodic
orbit**: a 6-vector IC `(x,y,z,ẋ,ẏ,ż)` + period `T`, closed under
`correct_general_periodic_3d`. There are **no Lambert legs**. The two ladders
share level *names* (V0..V5) and *floors* (spec §14) but NOT machinery. The C21
case is the second kind. **This is the key reason `v3_3d.run_v3_3d` does not
apply** — confirmed by reading the module: it takes `sequence`, `leg_tofs_days`,
`vinf_tuple_kms`, a `V2MoontourVerdict`, and re-solves Lambert legs.

There is, however, an in-file stub already pointing the right way:
`v3_3d.run_v3_periodic_regression` (lines 605–644) re-propagates a periodic IC
under `cr3bp.propagate` at rtol=atol=1e-13 and compares closure to a stored
value. **But its own docstring admits it is INTRA-model** — same scipy DOP853
propagator, "asserts the same propagator gives the same closure when asked
again" — and explicitly flags "A future REBOUND-IAS15 CR3BP backend
(force-callback + Coriolis) would be a stronger independent check; this
implementation is the wireable baseline." That future is exactly V3-3D below.

---

## 1. BUILT-vs-MISSING map, per level (3D CR3BP periodic orbit)

| Level | 3D-CR3BP-periodic meaning | Status | Where |
|---|---|---|---|
| **V0** | Internal consistency: corrector residual ≤ tol (idealized); Jacobi constant `C` constant along the arc; `degenerate_planar=False` if a genuinely 3D claim; Floquet multipliers on/off the unit circle consistent with the stability tag. | **BUILT (implicitly)** — but not packaged as a `run_v0_3d`. Every quantity exists (`corrector_residual`, `degenerate_planar` on `Periodic3DOrbit`; Jacobi via `cr3bp`; Floquet via the monodromy STM). | scattered: `v1_3d.run_v1_3d` precondition + `search/cr3bp_general_periodic_3d` |
| **V1** | Same-model closure: re-close IC under the 3D corrector AND independently re-propagate under **Radau** at rtol=atol=1e-12; both < spec 1 m/s floor. | **BUILT + EXERCISED** (C21 PASS, indep. closure 3.2e-12 km/s). | `data/validation/v1_3d.py` `run_v1_3d` |
| **V2** | Long-span bounded drift: ≥3 consecutive periods in same CR3BP model; max cumulative position drift ≤ 50,000 km floor. | **BUILT + EXERCISED** (C21 PASS, 6-cycle drift 1.7 mm, linear growth). | `data/validation/v2_3d.py` `run_v2_3d` |
| **V3** | **Independent integrator ARCHITECTURE**: re-propagate the IC in the CR3BP rotating frame under **REBOUND IAS15** (Gauss-Radau, architecturally distinct from scipy DOP853/Radau) and assert the same one-period 6D closure + the same V2 bounded-drift signature. | **MISSING for periodic orbits.** `run_v3_3d` is moontour-only; `run_v3_periodic_regression` is intra-model (same DOP853), explicitly a placeholder. | gap — new module `v3_3d_periodic.py` |
| **V4** | Real-ephemeris 3D: realise the rotating-frame IC as an inertial Earth-Moon-barycentric state at a real epoch, propagate under **DE440 + SPICE** real Earth/Moon (+Sun third body) with an independent integrator (REBOUND IAS15 / GMAT), assert bounded drift over ≥3 cycles in the dynamic rotating frame. The §14 V4 "spectral CR3BP stability does NOT imply V4 survival" note applies directly — this is where an all-|ev|=1 CR3BP orbit can still escape. | **MISSING for EM periodic orbits.** `v4_uranus.py` is a planet-frame Kepler+J2+third-body *moontour* gate (no SPICE Uranian kernels → documented scipy fallback); `real_closure.py` is heliocentric Lambert-cycler DE440. Neither realises a **rotating-frame CR3BP IC** into the EM inertial frame. `frames.py` transforms are **Sun-anchored heliocentric**, NOT Earth-Moon barycentric. | gap — new module + a new EM rotating↔inertial transform |
| **V5** | Novelty + expert review: canonical signature misses catalogue AND literature; human expert review; ideally independent group reproduction. | **PARTIALLY BUILT (process, not 3D-wired).** The C21 case already ran the *novelty* half by hand (JPL fullsweep + 3 WebSearch lit queries → KNOWN CLASS). No `run_v5_3d` packaging; the lit/JPL check is manual. V5 is intrinsically human — the deliverable is a *checklist + evidence bundle*, not an automatic pass. | gap — light: a V5 evidence-bundle assembler, not a solver |

**One-line summary:** V0 (latent, unpackaged), V1, V2 are BUILT for 3D periodic
orbits. **V3 (independent-integrator IAS15-in-CR3BP), V4 (real-eph EM 3D), and
V5 (novelty bundle) are MISSING.** V3-3D is the nearest, concretely motivated by
C21.

---

## 2. Architecture for the MISSING levels

### 2.1 V3-3D — REBOUND IAS15 in the CR3BP rotating frame (NEAREST GAP)

**Goal.** A genuinely independent integrator *architecture* (not the same scipy
solver re-run) reproduces (a) the one-period 6D closure V1 asserts and (b) the
≥3-cycle bounded-drift signature V2 asserts. Agreement ⇒ the periodicity is a
real property of the model, not a scipy-DOP853/Radau artefact.

**Why IAS15-in-CR3BP is non-trivial (the load-bearing risk).** REBOUND
integrates inertial N-body. The CR3BP rotating-frame EOM has Coriolis +
centrifugal terms that are **velocity-dependent** and live in an
`additional_forces` Python callback. Per the project memory *REBOUND variation +
custom force gotcha* and `nbody/propagator.py` lines 294–306, REBOUND's native
variational particles do **NOT** differentiate a Python `additional_forces`
callback. For V3-3D we do **not** need the STM (closure + drift are state-only),
so the variational hazard is sidestepped — but the velocity-dependent force
callback itself is the real work. Two candidate realisations:

* **(A) Rotating-frame force callback (preferred, truly independent).** Drive
  REBOUND with a single massless particle and an `additional_forces` callback
  implementing the full CR3BP rotating acceleration
  `a = -∇Ω(x,y,z) - 2 ω × v` (two point-mass primaries fixed at `(-µ,0,0)` and
  `(1-µ,0,0)`, ω = ẑ in nondim units). IAS15's Gauss-Radau stepper is the
  independent architecture vs scipy DOP853/Radau. This is the strongest V3.
  Risk: Coriolis is velocity-dependent → must confirm IAS15 handles a
  velocity-dependent callback stably (it does, but tolerance/step behaviour
  needs a parity gate, §4).
* **(B) Inertial two-primary N-body (fallback).** Realise the rotating IC into a
  *nondim inertial* frame (the CR3BP's own synodic→inertial rotation at ω=1,
  trivial closed form), add Earth + Moon as two massive bodies on circular
  orbits, integrate inertial IAS15 (no custom force — REBOUND native), rotate
  back each period. Cleaner numerically (no Python callback) but introduces the
  circular-orbit realisation of the primaries, which is *exactly* the CR3BP
  assumption, so it is still same-model. Acceptable as the baseline if (A)'s
  callback proves fragile.

Both share the verdict math with V1/V2. Recommend building **(A)** with **(B)**
as the documented fallback (mirrors the honest IAS15/LSODA fallback pattern
already in `v3_3d._ias15_propagate_planet_frame`).

**New module:** `src/cyclerfinder/data/validation/v3_3d_periodic.py`

```python
V3_PERIODIC_CLOSURE_FLOOR_NONDIM: Final[float]   # = 1e-7 (matches the existing
    # run_v3_periodic_regression placeholder; ~38 m EM — tighter than V1's 1e-6
    # but it is an AGREEMENT bar, not the spec floor)
V3_PERIODIC_DRIFT_AGREEMENT_FLOOR_KMS: Final[float]  # = 1.0 km — V3(IAS15) vs
    # V2(DOP853) per-cycle drift must agree to << the 50,000 km V2 floor

@dataclass(frozen=True)
class V3PeriodicVerdict3D:
    candidate_id: str
    integrator: str                       # "REBOUND IAS15 (rotating callback)" | "...inertial" | fallback label
    closure_residual_nondim_ias15: float  # ||X(T)-X(0)|| under IAS15
    closure_residual_kms_ias15: float
    per_cycle_drift_kms_ias15: tuple[float, ...]
    per_cycle_drift_kms_dop853: tuple[float, ...]   # from a V2Verdict3D, sliced
    drift_agreement_kms: float            # max_k |ias15[k] - dop853[k]|
    closure_floor_nondim: float
    drift_agreement_floor_kms: float
    n_cycles_propagated: int
    passes_v3: bool                       # closure < floor AND drift_agreement < floor AND n>=3
    degenerate_planar: bool               # carried from the candidate
    notes: str = ""

def run_v3_3d_periodic(
    candidate_id: str,
    state0: NDArray[np.float64],
    period_nondim: float,
    system: cr3bp.CR3BPSystem,
    *,
    v2_verdict: V2Verdict3D,              # supplies the DOP853 drift series to agree with
    n_cycles: int = 3,
    ias15_epsilon: float = 1e-12,
    closure_floor_nondim: float = V3_PERIODIC_CLOSURE_FLOOR_NONDIM,
    drift_agreement_floor_kms: float = V3_PERIODIC_DRIFT_AGREEMENT_FLOOR_KMS,
    prefer_rotating_callback: bool = True,  # (A); False → (B) inertial
    notes: str = "",
) -> V3PeriodicVerdict3D: ...
```

**Internal helper** (mirrors `_ias15_propagate_planet_frame`):

```python
def _ias15_propagate_cr3bp_rotating(
    state0_nondim: NDArray[np.float64],
    t_nondim: float,
    mu: float,
    *,
    epsilon: float = 1e-12,
    mode: str = "rotating_callback",   # | "inertial_twobody"
) -> tuple[NDArray[np.float64], str]:
    """Propagate a 6D CR3BP rotating state for t_nondim under REBOUND IAS15.
    Returns (state_f_nondim, integrator_label). Falls back to scipy LSODA
    (multistep BDF — distinct family from DOP853) with an honest label if
    rebound is unimportable, exactly as the moontour V3 does."""
```

**What V3-3D asserts + floor:** (i) IAS15 one-period closure < 1e-7 nondim
(~38 m EM) AND (ii) IAS15-vs-DOP853 per-cycle drift agreement < 1.0 km over ≥3
cycles. PASS ⇒ the V1/V2 verdict is integrator-architecture-independent.

### 2.2 V4-3D — real-ephemeris Earth-Moon 3D

**Goal.** Realise the rotating-frame CR3BP IC as an inertial Earth-Moon
**barycentric** state at a real epoch, propagate under real DE440 Earth+Moon
geometry (+Sun as third body — the §14 V4 note's "Sun switched on" killer),
independent integrator, and assert **bounded** drift over ≥3 cycles in the
*dynamic* rotating frame. This is the gate that can *fail* an all-|ev|=1 CR3BP
orbit (spec §14 note, lines 405–420 of spec.md).

**The missing infrastructure piece — an EM rotating↔inertial transform.**
`core/frames.py` only has **Sun-anchored heliocentric** transforms
(`to_rotating`/`from_rotating`, `to_rotating_dynamic` anchored on a heliocentric
body). There is **no Earth-Moon barycentric** equivalent. V4-3D needs:

```python
# new in core/frames.py (or a sibling core/frames_em.py)
def cr3bp_to_em_inertial(
    state_nondim: NDArray[np.float64],   # rotating CR3BP (x,y,z,ẋ,ẏ,ż), nondim
    epoch_et: float,                     # SPICE ET seconds
    system: cr3bp.CR3BPSystem,           # supplies l_km, t_s, mu for de-nondim + barycentre
    *,
    ephem,                               # SPICE / DE440 Earth+Moon state provider
) -> tuple[NDArray, NDArray]:            # (r_km, v_km_s) EM-barycentric inertial (J2000)
    """Pulsating-rotating CR3BP frame → EM-barycentric inertial at a real epoch.
    Uses the instantaneous Earth→Moon vector to define the rotating x-axis,
    the EM orbital angular velocity for ω (NON-uniform — pulsating frame, §11(c)),
    and l_km(t) = instantaneous |Earth-Moon| for the length scale. Inverse:
    em_inertial_to_cr3bp."""
```

This is the EM analogue of the heliocentric `to_rotating_dynamic`/`from_…` pair
already in `frames.py`, and follows the §11(c) "dynamic ephemeris frame +
tolerant verification" amendment (non-uniform ω, bounded-drift test). It is the
single genuinely new piece of *core* infra V4-3D needs.

**New module:** `src/cyclerfinder/data/validation/v4_3d_earthmoon.py`

```python
@dataclass(frozen=True)
class V4Verdict3DEarthMoon:
    candidate_id: str
    epoch_utc: str
    integrator: str                  # "REBOUND IAS15 DE440 (Earth+Moon+Sun)" | "GMAT R2022a" | fallback
    n_cycles_propagated: int
    per_cycle_drift_kms: tuple[float, ...]   # in the dynamic EM rotating frame
    max_drift_kms: float
    drift_floor_kms: float           # spec §14 V2-real = 200,000 km (real-eph breathing budget)
    hill_escape: bool                # did it leave the EM Hill sphere? (the §14-note failure mode)
    sun_third_body: bool             # Sun included (the killer term)
    passes_v4: bool
    notes: str = ""

def run_v4_3d_earthmoon(
    candidate_id, state0, period_nondim, system, *,
    epoch_utc: str = "2030-01-01T00:00:00",   # or scan a few epochs like real_closure
    n_cycles: int = 3,
    backend: str = "rebound",                  # | "gmat"
    include_sun: bool = True,
    drift_floor_kms: float = 200_000.0,
    notes: str = "",
) -> V4Verdict3DEarthMoon: ...
```

**Compose existing infra:**
* `verify/spice_kernels.py` (`ensure_leapseconds_kernel`, kernel furnish) +
  `verify/mission_spk.py`'s spiceypy `furnsh`/`spkezr`/`str2et` pattern — for
  Earth + Moon DE440 states. DE440 + LSK are already cached (per memory).
  Earth/Moon are core DE440 bodies (unlike the missing Uranian satellite
  kernels that forced `v4_uranus.py` into fallback) — so **V4-3D can run a real
  GMAT/SPICE path, not just a fallback.**
* `nbody/propagator.py` `RestrictedNBody("rebound")` IAS15 with a real-ephemeris
  Earth+Moon+Sun force model (the `nbody/forces.py` rails pattern) — the
  independent integrator.
* GMAT R2022a at `~/GMAT/R2022a` via `env -u DISPLAY ./GmatConsole --run`
  (memory: GMAT install) + `scripts/gmat_v4_generate.py` / `parse_gmat_report.py`
  as the template for the highest-fidelity V4 cross-check (Earth+Moon are in
  GMAT's bundled DE405/421/424 — no missing-kernel problem).

**What V4-3D asserts + floor:** bounded drift ≤ 200,000 km (spec §14 V2-real
floor, absorbing real-eph breathing) over ≥3 cycles in the *dynamic* EM rotating
frame, with the Sun included, and NO Hill-sphere escape. A clean *negative*
(C21 escapes once the Sun is on) is a valid, expected outcome per the §14 note —
it would mirror the heliocentric Hill-stability story and is a real result, not
a failure of the gate.

### 2.3 V5-3D — novelty + expert-review evidence bundle

V5 is intrinsically human (spec §14: "human expert review; ideally independent
reproduction by a separate group"). The build is **not a solver** — it is an
*evidence-bundle assembler* that packages what the C21 case did by hand:

```python
def assemble_v5_bundle_3d(
    candidate_id, state0, period_nondim, system, *,
    jpl_sweep_result,           # search/.../jpl_periodic_orbits novelty sweep
    literature_check_result,    # search/literature_check.py verdict (mandatory gate, memory)
    v1, v2, v3, v4,             # the four lower verdicts
) -> V5Bundle3D: ...           # → a markdown/JSON dossier for human sign-off
```

It enforces the discipline already exercised on C21: JPL fullsweep + the
mandatory `search/literature_check.py` published-record check (memory:
*literature-novelty check is baseline*), and refuses to stamp "novel" — it
*assembles the case* for a human. Golden/exemplar cross-checks must be
**sourced from the paper** (memory: golden-sourced-only), never our computed IC.

---

## 3. Bite-sized TDD build sequence

**Phase A — V3-3D (do first; concretely motivated by C21).**
1. **RED:** test `_ias15_propagate_cr3bp_rotating` against a *known* CR3BP arc —
   propagate the C21 IC one period, assert IAS15 final state matches the
   `cr3bp.propagate` (DOP853) final state to < 1e-7 nondim. (Parity gate for the
   force callback — this is the load-bearing correctness check.)
2. **GREEN:** implement helper, mode (A) rotating callback; LSODA fallback.
3. **RED:** test `run_v3_3d_periodic` on C21 → `passes_v3=True`, closure < 1e-7,
   drift agreement vs the existing C21 `V2Verdict3D` < 1 km.
4. **GREEN:** implement `run_v3_3d_periodic`, `V3PeriodicVerdict3D`.
5. **RED/GREEN:** negative test — a deliberately *non-periodic* IC must FAIL V3
   (closure or drift-agreement blows the floor). Guards against a vacuous pass.
6. Wire `__all__` + `data/validation/__init__.py`; run
   `uv run ruff check . && ruff format --check .` and
   `uv run pytest tests/data -q` (memory: ruff pre-commit; run-all-ratchets).

**Phase B — EM rotating↔inertial transform (core infra for V4).**
7. **RED:** round-trip identity test `em_inertial_to_cr3bp(cr3bp_to_em_inertial(x))
   ≈ x` to ~1e-10 (mirrors the existing `from_rotating_omega_vec` round-trip
   test), plus a parity test that at a *circular* synthetic EM ephemeris the
   transform reduces to the uniform CR3BP synodic frame.
8. **GREEN:** implement the EM transform pair in `core/frames.py` (or
   `core/frames_em.py`).

**Phase C — V4-3D.**
9. **RED:** test that `run_v4_3d_earthmoon` propagates the C21 IC ≥3 cycles on
   DE440 with the Sun on and reports a finite drift + a definite
   `hill_escape` boolean. (Assert the *machinery*, not a specific pass — the
   physics verdict is whatever the math says.)
10. **GREEN:** implement REBOUND DE440 backend (Earth+Moon+Sun) + dynamic-frame
    drift.
11. **(stretch) GMAT cross-check** via `scripts/gmat_v4_generate.py` template.

**Phase D — V5-3D bundle.**
12. `assemble_v5_bundle_3d` + a test that it *refuses* to stamp novel when the
    literature check returns a hit (C21 regression: KNOWN CLASS).

---

## 4. Risks + open questions

* **(R1, load-bearing) Rotating→inertial / force-callback fidelity for IAS15.**
  V3-3D mode (A) puts the velocity-dependent Coriolis term in a Python
  `additional_forces` callback. Memory (*REBOUND variation + custom force
  gotcha*) warns this is exactly where REBOUND's native machinery is blind — but
  that gotcha is about *variational particles*, which V3-3D does **not** use
  (state-only). The residual risk is IAS15 step/tolerance behaviour under a
  velocity-dependent callback. **Mitigation:** the Phase-A step-1 parity gate
  (IAS15 vs DOP853 one-period) is mandatory before any verdict; if (A) is
  fragile, fall back to mode (B) inertial two-primary (no callback) and document
  it honestly, exactly as the moontour V3 documents its LSODA fallback.

* **(R2) `frames.py` has no Earth-Moon transform.** Confirmed: all existing
  transforms are Sun-anchored heliocentric. V4-3D needs a *new* EM-barycentric
  pulsating-rotating transform. This is real new core work, not a wrapper. It
  must be round-trip + circular-parity tested (Phase B) before V4 trusts it.

* **(R3) V4 is where an all-|ev|=1 CR3BP orbit can legitimately die.** The §14
  note (spec.md 405–420) is explicit: spectral CR3BP stability does NOT imply V4
  survival; the Sun tide is ~30% of Earth's gravity at the lunar distance. **A
  C21 V4 FAIL is a plausible, valid, publishable-as-negative outcome**, not a
  bug. The cheap predictor (Hill-amplitude classifier in the note) could be run
  first to set expectations.

* **(R4) C21 is the natural first test candidate but is KNOWN-CLASS** (#444:
  Antoniadou & Libert 2019). That is fine for *exercising the gates* (it is a
  correct stable spatial 2:1 EM orbit) — but it must NOT be re-stamped novel;
  the V5 bundle's literature gate is the guardrail. For a *golden* test the
  EXPECTED IC/period/Jacobi must come from the paper, never our computed value
  (memory: golden-sourced-only).

* **(R5) Catalogue discipline.** No level here admits to `catalogue.yaml`. V3/V4
  PASS are necessary-not-sufficient; admission needs V4 + V5 + the literature
  miss. Any future catalogue row edit triggers the full ratchet suite (memory:
  catalogue-edits-run-all-ratchets) — out of scope for this gauntlet build, but
  noted so it is not forgotten if a genuinely novel candidate ever clears V5.

* **(Q1) Epoch handling for V4.** Single fixed epoch vs a small scan (as
  `real_closure._resolve_real_t_start` does)? A periodic CR3BP orbit is
  epoch-agnostic in idealised terms, but real-eph realisation is not — the EM
  geometry/Sun phase at injection matters. Recommend a small epoch scan
  (a few values across a year) and report the best, mirroring the heliocentric
  phase-match convention.

* **(Q2) Is a `run_v0_3d` worth packaging?** V0 quantities all exist but are
  scattered. Low effort to package (corrector residual + Jacobi-constancy +
  `degenerate_planar` + Floquet-tag consistency) into one `V0Verdict3D` for a
  clean ladder. Optional, low priority vs V3/V4.

---

## 5. Recommended first build

**V3-3D (`run_v3_3d_periodic` + `_ias15_propagate_cr3bp_rotating`).** Nearest
gap, concretely motivated by the C21 case that already PASSED V1+V2 this session,
and it composes cleanly on infra that already exists (the REBOUND IAS15 setup in
`v3_3d.py`/`nbody/propagator.py`, the V2 drift series to agree with). **Key
risk: R1** — the velocity-dependent Coriolis force callback under IAS15; gate it
with the mandatory IAS15-vs-DOP853 one-period parity test (Phase-A step 1) and
keep the inertial two-primary mode (B) as the documented fallback.
