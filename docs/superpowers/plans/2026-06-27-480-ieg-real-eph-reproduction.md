# #480 M1 — IEG Real-Ephemeris Reproduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce the Hernandez-Jones-Jesick 2017 (AAS 17-608) Io-Europa-Ganymede triple cycler in **real Galilean-moon ephemeris** and validate it against the paper's published invariants, ending with a V4-tier known-reproduction reported for human admission.

**Architecture:** Extend the proven V4-Uranus validation pipeline to Jupiter, but feed it a **real-ephemeris Galilean-moon position provider** (`jup365.bsp`) so the spacecraft is corrected against true Io/Europa/Ganymede positions — not idealized circular orbits (the #473 wall). The existing `nbody/shooter.py` already accepts an `Ephemeris` + `bodies` and does multiple-shooting with an analytic STM Jacobian; the new work is (a) a Galilean `Ephemeris`, (b) a Jupiter `v4_jupiter` gauntlet mirroring `v4_uranus.py`, (c) a seed adapter from the #473 IEG geometry, and (d) a sourced invariant golden.

**Tech Stack:** Python 3.11, uv, spiceypy (jup365.bsp + de440 + naif0012.tls), scipy (DOP853, least_squares LM), the existing `cyclerfinder.nbody` shooter + `cyclerfinder.data.validation` V3/V4 framework, numba-accelerated kernels (#475).

**Spec:** `docs/superpowers/specs/2026-06-26-480-resonance-lock-moontour-generator-design.md`
**Grounded golden source:** `docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md`

**Scope (YAGNI):** IEG only; ballistic (no powered DSM); no discovery enumerator; no cross-system continuation; **no catalogue self-admission** (V4-ceiling known-reproduction → report for human admission).

---

## File Structure

| File | Responsibility | Create/Modify |
|---|---|---|
| `src/cyclerfinder/nbody/galilean_ephem.py` | Real-eph `Ephemeris` serving Io/Europa/Ganymede + Jupiter states from `jup365.bsp` | Create |
| `src/cyclerfinder/verify/spice_kernels.py` | add `ensure_jup365_kernel()` (path/fetch helper) | Modify |
| `src/cyclerfinder/data/validation/v4_jupiter.py` | Jupiter J2 + Galilean-third-body gauntlet (analog of `v4_uranus.py`), using the real-eph ephem | Create |
| `src/cyclerfinder/search/ieg_seed.py` | #473 IEG geometry → `ShootingSeed` | Create |
| `tests/nbody/test_galilean_ephem.py` | ephem parity vs jup365 / NAIF spkezr | Create |
| `tests/data/test_v4_jupiter.py` | v4_jupiter end-to-end structural test (pattern: `test_v4_uranus.py`) | Create |
| `tests/verify/test_ieg_reproduction_golden.py` | invariant golden vs the Hernandez-17-608 digest, citation-anchored | Create |

The resonance-locked corrector is **not a new module** — it is the existing `nbody/shooter.py::shoot(...)` driven with the Galilean ephem + the IEG seed. Unit "3" in the spec is realized by a thin driver inside `v4_jupiter.py` / the golden test, not new corrector code.

---

## Task 0: Confirm the `Ephemeris` interface (read-only spike)

**Files:** none (investigation). Output: a 6-line note appended to the bottom of this plan under "## Task 0 findings".

- [ ] **Step 1: Locate and read the `Ephemeris` class** the shooter consumes.

Run: `grep -rn "class Ephemeris\|def __init__\|def state\|def position\|def body_state\|astropy" src/cyclerfinder/nbody/ | grep -i ephem`
Then read the file it lives in (likely `src/cyclerfinder/nbody/ephem.py` or similar).

- [ ] **Step 2: Record verbatim** in "## Task 0 findings" at the bottom of this plan: the `Ephemeris` constructor signature, the method that returns a body's state at an epoch (name + signature + return type/units + frame), and how `bodies` strings map to ephemeris lookups. This is the contract `galilean_ephem.py` must satisfy or subclass.

- [ ] **Step 3: Confirm jup365 availability.** Run: `ls -la ~/dev/references/kernels/jup365.bsp` and `uv run python -c "import spiceypy as s; s.furnsh('/home/bruce/dev/references/kernels/jup365.bsp'); print(s.spkezr('502', 0.0, 'J2000', 'NONE', '599')[0])"` (502=Europa, 599=Jupiter barycenter). Expected: a 6-vector (km, km/s). Record the NAIF IDs used (Io 501, Europa 502, Ganymede 503, Jupiter 599 / 5 for barycenter).

---

## Task 1: jup365 kernel helper

**Files:**
- Modify: `src/cyclerfinder/verify/spice_kernels.py`
- Test: `tests/nbody/test_galilean_ephem.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/nbody/test_galilean_ephem.py
from pathlib import Path
from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel

def test_jup365_kernel_path_exists():
    p = ensure_jup365_kernel()
    assert Path(p).is_file()
    assert p.endswith("jup365.bsp")
```

- [ ] **Step 2: Run it to verify it fails** — `uv run pytest tests/nbody/test_galilean_ephem.py::test_jup365_kernel_path_exists -v` → FAIL (no `ensure_jup365_kernel`).

- [ ] **Step 3: Implement** `ensure_jup365_kernel()` in `spice_kernels.py` mirroring `ensure_leapseconds_kernel`: return the local path `~/dev/references/kernels/jup365.bsp` if present (expanduser), else raise a clear `RuntimeError` naming the NAIF download URL (`https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp`). Do NOT auto-download a 1GB+ kernel; the floor is "present-or-clear-error".

- [ ] **Step 4: Run** the test → PASS.

- [ ] **Step 5: Commit** — `git add src/cyclerfinder/verify/spice_kernels.py tests/nbody/test_galilean_ephem.py && git commit -m "verify: #480 jup365 Galilean-satellite kernel helper"`

---

## Task 2: Galilean real-ephemeris provider

**Files:**
- Create: `src/cyclerfinder/nbody/galilean_ephem.py`
- Test: `tests/nbody/test_galilean_ephem.py`

The `Ephemeris` contract from Task 0 dictates the exact subclass/shape. The module furnishes `jup365.bsp` + de440 + naif0012.tls and serves Jupiter-centered J2000 states for Io/Europa/Ganymede/Jupiter via `spkezr`.

- [ ] **Step 1: Write the failing parity test** — assert the provider's Europa state at a fixed epoch matches a direct `spiceypy.spkezr("502", et, "J2000", "NONE", "599")` to < 1 km / < 1e-6 km/s, for 3 epochs.

```python
def test_galilean_ephem_matches_spkezr():
    import spiceypy as s
    from cyclerfinder.nbody.galilean_ephem import GalileanEphemeris, ensure_furnished
    ensure_furnished()
    eph = GalileanEphemeris()
    for et in (0.0, 7.05 * 86400.0, 28.22 * 86400.0):
        ours = eph.body_state("Europa", et)            # contract from Task 0
        ref, _ = s.spkezr("502", et, "J2000", "NONE", "599")
        assert max(abs(a - b) for a, b in zip(ours[:3], ref[:3])) < 1.0    # km
        assert max(abs(a - b) for a, b in zip(ours[3:], ref[3:])) < 1e-6   # km/s
```

- [ ] **Step 2: Run → FAIL** (module missing).

- [ ] **Step 3: Implement** `galilean_ephem.py`: `ensure_furnished()` furnishes jup365 (Task 1) + de440 (`astropy_de440_bsp_path`) + leapseconds (`ensure_leapseconds_kernel`) idempotently; `GalileanEphemeris` implements the Task-0 `Ephemeris` contract, mapping {"Io":"501","Europa":"502","Ganymede":"503","Jupiter":"599"} and returning Jupiter-centered J2000 states. Match whatever method name/units/frame Task 0 recorded (adapt km↔canonical if the contract requires).

- [ ] **Step 4: Run → PASS.** Also add `test_galilean_ephem_io_period_sanity`: propagate Io one period (1.769 d) and assert it returns near its start position (< 5000 km), confirming real-eph orbital motion (not a frozen point).

- [ ] **Step 5: Commit** — `git add src/cyclerfinder/nbody/galilean_ephem.py tests/nbody/test_galilean_ephem.py && git commit -m "nbody: #480 Galilean real-ephemeris provider (jup365 J2000 states)"`

---

## Task 3: IEG seed adapter

**Files:**
- Create: `src/cyclerfinder/search/ieg_seed.py`
- Test: `tests/data/test_v4_jupiter.py` (shared test file; this task adds the seed tests)

Convert the #473 IEG geometry (sequence `("Io","Europa","Ganymede","Io")`, resonance-locked leg ToFs at the EGGIE-like configuration, V∞ from the digest) into a `ShootingSeed` (`nbody/shooter.py:92`). The EGGIE target invariants (digest Table 4): leg ToFs (days) Ganymede 1.59, Ganymede 8.60, Io 7.34, Europa 10.69; V∞ Europa 9.12, Ganymede 7.07, Ganymede 7.07, Io 8.38; T_syn = 7.05 d; total ToF 28.22 d.

- [ ] **Step 1: Write the failing test** — `ieg_eggie_seed()` returns a `ShootingSeed` with `sequence == ("Europa","Ganymede","Ganymede","Io","Europa")`, `len(node_states)==5`, `tofs == [1.59, 8.60, 7.34, 10.69]` (within 0.01), `period_days == pytest.approx(28.22, abs=0.05)`, and finite 6-vectors per node (built from circular Galilean states at the encounter epochs as the warm start; the corrector refines to real-eph).

- [ ] **Step 2: Run → FAIL.**

- [ ] **Step 3: Implement** `ieg_seed.py`: build per-node Cartesian states from `SATELLITES[moon].sma_km` circular orbits at the resonance-locked encounter phases (warm start only — exact positions come from the corrector), populate `vinf_in`/`vinf_out` from the digest V∞ tuple, set `slack_leg` to the longest leg, `period_days=28.22`. Constants sourced from the digest (cite it in a module docstring); the seed is a *guess*, so circular-orbit construction is acceptable here (NOT a sourced golden).

- [ ] **Step 4: Run → PASS.**

- [ ] **Step 5: Commit** — `git add src/cyclerfinder/search/ieg_seed.py tests/data/test_v4_jupiter.py && git commit -m "search: #480 IEG EGGIE seed adapter (#473 geometry -> ShootingSeed)"`

---

## Task 4: Real-eph corrector convergence (shooter + Galilean ephem)

**Files:**
- Test: `tests/data/test_v4_jupiter.py`
- (No new src — drives `nbody/shooter.py::shoot` with `GalileanEphemeris`.)

- [ ] **Step 1: Write the test (load-bearing, math decides — do NOT assert PASS/FAIL of reproduction yet, per `feedback_orbit_closure_discipline`)**

```python
def test_ieg_seed_corrects_against_real_galilean_ephemeris():
    from cyclerfinder.nbody.shooter import shoot
    from cyclerfinder.nbody.galilean_ephem import GalileanEphemeris, ensure_furnished
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed
    ensure_furnished()
    res = shoot(ieg_eggie_seed(), ephem=GalileanEphemeris(),
                bodies=("Io","Europa","Ganymede"), jacobian="stm",
                max_nfev=50, max_wall_sec=300.0)
    # Structural assertions only:
    assert res.sequence == ("Europa","Ganymede","Ganymede","Io","Europa")
    assert len(res.corrected_states) == 5
    import math
    assert math.isfinite(res.defect_norm)
    assert res.defect_norm <= res.seed_defect_norm        # the corrector did not diverge
```

- [ ] **Step 2: Run** — `uv run pytest tests/data/test_v4_jupiter.py::test_ieg_seed_corrects_against_real_galilean_ephemeris -v --timeout=600`. Expected on first run: it may FAIL if the seed is too far for the shooter — that is a real finding, not a bug. If it fails to converge, record the `defect_norm` vs `seed_defect_norm` and the failure mode; the homotopy fallback (spec Unit 3 escalation: ramp from circular→real-eph) is the next rung, added as Task 4b only if needed.

- [ ] **Step 3: If non-convergent, implement the homotopy rung** (Task 4b): a helper in `ieg_seed.py` that continues the seed from circular-moon ephem to `GalileanEphemeris` over a perturbation parameter λ∈[0,1], re-seeding `shoot` at each step. Add `test_ieg_homotopy_reaches_real_ephemeris`. (Skip this task entirely if Step 2 converges.)

- [ ] **Step 4: Run → PASS** (structural).

- [ ] **Step 5: Commit** — `git add tests/data/test_v4_jupiter.py [src/cyclerfinder/search/ieg_seed.py] && git commit -m "data: #480 IEG corrects against real Galilean ephemeris (structural)"`

---

## Task 5: v4_jupiter gauntlet

**Files:**
- Create: `src/cyclerfinder/data/validation/v4_jupiter.py`
- Test: `tests/data/test_v4_jupiter.py`

Mirror `v4_uranus.py` (`run_v4_uranus -> V4UranusVerdict`, the per-cycle verdict, the drift/agreement floors) with Jupiter constants and the Galilean ephem. Constants (sourced, cite in docstring): `JUPITER_J2 = 1.4736e-2` (Jacobson JUP310/JUP365 — confirm the exact value + source in the kernel comments during impl; the v4_uranus pattern cites Jacobson), `JUPITER_R_EQ_KM = 71492.0`, `JUPITER_MU = PRIMARIES["Jupiter"]`, `JOVIAN_PERTURBER_MOONS = ("Io","Europa","Ganymede","Callisto")`, reuse `V4_AGREEMENT_FLOOR_KMS = 50_000.0`, `V4_N_CYCLES_MIN = 3`.

- [ ] **Step 1: Write the failing structural test** (clone `tests/data/test_v4_uranus.py::test_silver_v4_runs_end_to_end_and_produces_verdict` shape): run V2→V3 to get a `V3Verdict3D` for the IEG tour, call `run_v4_jupiter(...)`, assert `isinstance(v, V4JupiterVerdict)`, `len(v.per_cycle)==n_cycles`, cycle-0 drift==0.0, V3 series carried + matched. Do NOT assert pass/fail.

- [ ] **Step 2: Run → FAIL** (`run_v4_jupiter` missing).

- [ ] **Step 3: Implement** `v4_jupiter.py` by adapting `v4_uranus.py`: same `run_v4_*` signature + frozen `V4JupiterVerdict`/`V4CycleVerdictJupiter` dataclasses, swap Uranus constants→Jupiter, swap `URANIAN_PERTURBER_MOONS`→`JOVIAN_PERTURBER_MOONS`, and source moon positions from `GalileanEphemeris` (the real-eph difference vs v4_uranus's circular fallback). Keep the DOP853 + J2 + third-body force model.

- [ ] **Step 4: Run → PASS** (structural). Add `test_v4_jupiter_constants_sourced` asserting the constant values match their docstring sources.

- [ ] **Step 5: Commit** — `git add src/cyclerfinder/data/validation/v4_jupiter.py tests/data/test_v4_jupiter.py && git commit -m "data/validation: #480 v4_jupiter gauntlet (Jupiter J2 + Galilean real-eph third-body)"`

---

## Task 6: Invariant golden — reproduction vs Hernandez 17-608

**Files:**
- Test: `tests/verify/test_ieg_reproduction_golden.py`

The success gate. Compares the converged + V3/V4-validated IEG tour's invariants to the **sourced** digest values, anchored on the verified citation. EXPECTED side traces to the paper (NOT code-computed) per `feedback_golden_tests_sourced_only`.

- [ ] **Step 1: Write the golden test**

```python
# tests/verify/test_ieg_reproduction_golden.py
import math, pytest
from cyclerfinder.search.literature_check import anchor_for_key

# SOURCED from docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md
# (Hernandez-Jones-Jesick 2017, AAS 17-608, Table 4 — EGGIE 4-synodic ballistic cycler)
EGGIE_TSYN_DAYS = 7.05
EGGIE_TOTAL_TOF_DAYS = 28.22
EGGIE_GANYMEDE_VINF_KMS = 7.07          # both Ganymede flybys equal
EGGIE_TOTAL_DV_MS = 0.70                 # near-ballistic
LAPLACE_RATIO = (1, 2, 4)                # Ganymede:Europa:Io

def test_ieg_citation_is_decision_grade():
    a = anchor_for_key("hernandez-2017-ieg-608")     # raises if unresolved/unverified (#484/#486)
    assert a.system == "jovian"
    assert {"Io","Europa","Ganymede"} <= set(a.body_set)

def test_eggie_reproduction_matches_published_invariants():
    from cyclerfinder.search.ieg_seed import ieg_eggie_seed
    from cyclerfinder.nbody.shooter import shoot
    from cyclerfinder.nbody.galilean_ephem import GalileanEphemeris, ensure_furnished
    ensure_furnished()
    res = shoot(ieg_eggie_seed(), ephem=GalileanEphemeris(),
                bodies=("Io","Europa","Ganymede"), jacobian="stm",
                max_nfev=80, max_wall_sec=600.0)
    assert res.converged, f"IEG EGGIE did not converge in real ephemeris: defect={res.defect_norm}"
    # Reproduction: V∞ at the Ganymede encounters within tolerance of the published 7.07 km/s
    gany_vinf = [v for moon, v in zip(res.sequence, res.vinf_per_encounter_kms) if moon == "Ganymede"]
    for v in gany_vinf:
        assert v == pytest.approx(EGGIE_GANYMEDE_VINF_KMS, abs=0.5)
    # Near-ballistic: correction ΔV small (paper 0.70 m/s; real-eph maintenance grows, allow a band)
    assert res.correction_dv_kms * 1000.0 < 50.0     # m/s; EIGE-class real-eph maintenance ceiling
```

- [ ] **Step 2: Run → likely FAIL initially** (convergence + tolerance reality). This is the real scientific test: if it converges and matches the published V∞ + stays near-ballistic, the reproduction succeeds. If not, record the gap honestly (it may need the homotopy rung from Task 4b, or it characterizes a real-eph divergence — both legitimate outcomes; do NOT loosen tolerances to force a pass).

- [ ] **Step 3: Make it pass legitimately** — only via the homotopy/seed-refinement rungs already specced, never by relaxing the sourced tolerances below physical defensibility (V∞ abs=0.5 km/s, maintenance < ~50 m/s over the validated horizon are the documented bands; if reality exceeds them, that is the reported result).

- [ ] **Step 4: Run → PASS (or documented characterized result).**

- [ ] **Step 5: Commit** — `git add tests/verify/test_ieg_reproduction_golden.py && git commit -m "verify: #480 IEG reproduction golden vs Hernandez 2017 (AAS 17-608) invariants"`

---

## Task 7: Verdict note + handoff (NO catalogue self-admission)

**Files:**
- Create: `docs/notes/2026-06-27-480-ieg-reproduction-verdict.md`

- [ ] **Step 1: Write the verdict note** — record: did the IEG tour reproduce in real Galilean ephemeris? the converged V∞/ToF/altitude vs the digest tables; the V3/V4 verdicts; the reproduction tolerance outcome; and the honest framing (V4-ceiling known-reproduction — a published tour cannot reach V5). State explicitly that **no catalogue.yaml row was self-admitted**; the reproduction is reported here for separate human admission.
- [ ] **Step 2: Run the full ratchet** — `uv run pytest tests/data tests/search tests/nbody tests/verify -q 2>&1 | tee /tmp/480_ratchet.log | tail -5`. Expected: green; the new tests added, no existing result changed; no census/tier change (no catalogue edit).
- [ ] **Step 3: Commit** — `git add docs/notes/2026-06-27-480-ieg-reproduction-verdict.md && git commit -m "docs: #480 M1 IEG real-eph reproduction verdict"`

---

## Self-Review

**Spec coverage:** Unit 1 (Jovian eph) → Tasks 1-2 (+ the real-eph provider, the A-fork addition). Unit 2 (seed adapter) → Task 3. Unit 3 (resonance-locked corrector) → Task 4 (realized as shoot()+ephem, not new code). Unit 4 (v4_jupiter gauntlet) → Task 5. Unit 5 (invariant golden) → Task 6. Staging/scope (no admission) → Task 7. The exact-match stretch (M1-stretch) is deliberately out of this plan (data-gated). ✓

**Placeholder scan:** Task 0 is an explicit investigation whose findings parameterize Tasks 1-2 (the `Ephemeris` contract is the one interface not fully mapped pre-plan — Task 0 closes it before any code). All other tasks carry real signatures from the grounding exploration. The `JUPITER_J2` value is flagged for confirmation against the kernel comments during Task 5 impl (not invented).

**Type consistency:** `ShootingSeed`/`ShootResult` (shooter), `V4JupiterVerdict` (mirrors `V4UranusVerdict`), `GalileanEphemeris.body_state` (pending Task-0 contract — the one name to reconcile), `anchor_for_key` (literature_check) — used consistently across tasks.

**Risk honesty:** Tasks 4 and 6 are explicitly "math decides" — the plan does NOT assume reproduction succeeds; a characterized real-eph divergence is a legitimate, reportable M1 outcome (and would itself validate #473's "model is the wall" with the *correct* model). Tolerances are sourced and must not be loosened to force a pass.

## Task 0 findings
_(to be filled by the implementer in Task 0 before Task 1)_
