# Family-Targeted Reproduction (Resonance-Anchored Construction) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Work directly on `main` — do NOT branch (project rule).

**Goal:** Reproduce S1L1 (and, at scale, the ~200 SnLm catalogue rows) by **constructing** each cycler from its *sourced* heliocentric orbit and differential-correcting to the real ephemeris — instead of free-optimising legs to hit V∞ targets (which lands off-family).

**Why this approach (and how S1L1 was originally found):** The 2002/2004 discovery method (McConaghy/Longuski/Byrnes; Russell 2004; Hollister & Menning 1970) was **not** free optimisation. It was resonance-first construction: pick the spacecraft's resonant heliocentric orbit, place the Earth/Mars encounters at the orbit's crossings under V∞-matching (ballistic flybys rotate but preserve |V∞|), solve in the circular-coplanar model, then **differentially correct** the coplanar seed onto the inclined/elliptic real ephemeris (Newton / steepest-descent). Our free optimiser keeps landing in off-family basins precisely because it discards this structure. We adopt their method: S1L1's V∞ is a *property* of its sourced orbit (a=1.30 AU, e=0.257), not something to search for.

**Architecture:** A new `construct_resonant_cycler` builds the heliocentric ellipse from sourced `(a, e)`, finds the true anomalies where `r == r_planet`, derives encounter epochs + leg ToFs from Kepler's equation, and computes per-encounter V∞ analytically. A golden gate cross-checks the *computed* V∞ against the *independently sourced* anchors (5.65/3.05 for S1L1 — two separately published quantities, so this is a legitimate, non-circular golden test). That coplanar seed then initialises the existing real-eph optimiser (`optimise_cell_ephemeris`, which already accepts `tof_seed_days`) for in-family DE440 closure.

**Tech Stack:** Python 3.11, numpy-only production path; `core/kepler.py` (has `propagate`), `core/constants.py` (`PLANETS` with `sma_au`, `ecc`, `mean_motion_deg_day`), `core/ephemeris.py`, `search/optimize.py::optimise_cell_ephemeris`, `search/phase_match.py` (multi-rev resolver). pytest + uv + ruff + mypy.

**Provenance contract:** `(a, e)` are sourced (Rogers 2012 Table 1 / spec §9); V∞ 5.65/3.05 are *independently* sourced (spec §9). Constructing from `(a,e)` and checking V∞ against the independent anchors is a cross-validation, not circular. Encounter epochs/ToFs the construction derives are COMPUTED and labelled so.

---

## Background (read first)

- `src/cyclerfinder/core/kepler.py:65` — `propagate(...)` (Kepler). No elements→state helper exists yet; Task 1 adds one.
- `src/cyclerfinder/core/constants.py` — `PLANETS[code]` carries `sma_au`, `ecc`, `mean_motion_deg_day`, `radius_eq_km`; `MU_SUN_KM3_S2`, `AU_KM`, `SECONDS_PER_DAY`.
- `src/cyclerfinder/search/optimize.py` — `optimise_cell_ephemeris(cell, ephem, *, vinf_cap, priority_date_iso, vinf_targets_kms, tof_seed_days, ...)` already takes a ToF seed; this plan supplies a *family-correct* one.
- `scripts/diagnose_s1l1_realeph.py` — current evidence the free search lands off-family (V_E≈26.5).
- Sourced S1L1: a=1.30 AU, e=0.257, peri=0.97, apo=1.64, V∞_E=5.65, V∞_M=3.05, E→M ~154 d, 2-synodic period.

---

## Task 1: `coe_to_rv` — orbital elements → heliocentric state

**Files:** Create helper in `src/cyclerfinder/core/kepler.py`; Test `tests/core/test_coe_to_rv.py`.

- [ ] **Step 1: Failing test** — a circular orbit (e=0) at a=1 AU gives |r|=AU, |v|=sqrt(mu/r), r·v=0; an e>0 orbit at perihelion (ν=0) gives r=a(1-e), and speed = vis-viva.

```python
import numpy as np
from math import sqrt
from cyclerfinder.core.kepler import coe_to_rv
from cyclerfinder.core.constants import MU_SUN_KM3_S2, AU_KM

def test_coe_to_rv_circular():
    r, v = coe_to_rv(a_km=AU_KM, e=0.0, true_anom_rad=0.0, mu=MU_SUN_KM3_S2)
    assert abs(np.linalg.norm(r) - AU_KM) < 1.0
    assert abs(np.linalg.norm(v) - sqrt(MU_SUN_KM3_S2 / AU_KM)) < 1e-6
    assert abs(float(np.dot(r, v))) < 1e-3

def test_coe_to_rv_perihelion():
    a = 1.30 * AU_KM
    r, v = coe_to_rv(a_km=a, e=0.257, true_anom_rad=0.0, mu=MU_SUN_KM3_S2)
    assert abs(np.linalg.norm(r) - a * (1 - 0.257)) < 1.0
```

- [ ] **Step 2:** Run → FAIL (no `coe_to_rv`).
- [ ] **Step 3:** Implement planar `coe_to_rv(a_km, e, true_anom_rad, mu, *, arg_peri_rad=0.0)` via the standard perifocal→inertial formulae (planar: inclination 0). Return `(r[3], v[3])` numpy.
- [ ] **Step 4:** Run → PASS.
- [ ] **Step 5:** Commit `core/kepler: add coe_to_rv (orbital elements -> heliocentric state)`.

## Task 2: `construct_resonant_cycler` — sourced orbit → encounter schedule + V∞

**Files:** Create `src/cyclerfinder/search/resonant_construct.py`; Test `tests/search/test_resonant_construct.py`.

- [ ] **Step 1: Failing test (the S1L1 golden cross-check)** — construct from sourced (a,e); assert the *computed* V∞ at the Earth (r=1 AU) and Mars (r=1.524 AU) crossings match the *independently sourced* anchors within tolerance.

```python
from cyclerfinder.search.resonant_construct import construct_resonant_cycler

def test_s1l1_vinf_from_sourced_orbit_matches_independent_anchors():
    # a/e sourced (Rogers 2012 Table 1); V_inf 5.65/3.05 INDEPENDENTLY sourced (spec §9).
    res = construct_resonant_cycler(a_au=1.30, e=0.257, bodies=("E", "M"))
    assert abs(res.vinf_kms["E"] - 5.65) < 0.3
    assert abs(res.vinf_kms["M"] - 3.05) < 0.3
```

- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Implement: build the ellipse (a,e); for each body solve `r(ν)=a(1-e²)/(1+e cos ν)=r_body` for the crossing true anomalies; `coe_to_rv` at each crossing; V∞ = |orbit_velocity − planet_circular_velocity| (coplanar, planet on circular orbit at `sma_au`). Return a dataclass with `vinf_kms` per body, the crossing true anomalies, and leg ToFs (via Kepler's time-of-flight between crossings). **Compute only; never assert against the (a,e) you were given.**
- [ ] **Step 4:** Run → PASS if the orbit's geometry is self-consistent with the published V∞ (it should be — both are properties of the real S1L1). If it does NOT match, that is a *finding* (the published a/e and V∞ are from different model fidelities); record it honestly and keep the test xfail rather than forcing it.
- [ ] **Step 5:** Commit `search: resonant-construction of cycler V_inf from sourced orbital elements`.

> **This task alone may be the S1L1 idealised-reproduction win:** it validates S1L1's published V∞ from its published orbit, with no free optimisation — exactly the McConaghy/Russell construction.

## Task 3: Family-correct seed → real-eph closure

**Files:** Modify `src/cyclerfinder/search/optimize.py` (a thin `seed_from_resonant_construction` helper or reuse `tof_seed_days`); Test `tests/search/test_s1l1_real_rediscovery.py` (the existing xfail gate).

- [ ] **Step 1:** Feed the Task-2 leg ToFs + a priority epoch into `optimise_cell_ephemeris` (E-M-E-E cell, multi-rev L1) as `tof_seed_days`, with `vinf_targets_kms={"E":5.65,"M":3.05}`. The seed now starts *in-family*.
- [ ] **Step 2:** Run the S1L1 real-eph gate. If V∞ now lands within 0.4 km/s of 5.65/3.05 → **flip the gate from xfail to a passing assert** (sourced-anchor close, the headline win). If it still lands off-family, keep xfail and record the achieved V∞ + the gap (do NOT loosen tolerance) — the construction seed narrowed it but DE440 differential correction needs more (note the next lever: a Newton differential-corrector on epochs preserving V∞-matching, per Russell 2004).
- [ ] **Step 3:** Commit (message reflects the actual outcome — pass or narrowed-xfail).

## Task 4: Scale to the SnLm catalogue rows

**Files:** Test `tests/data/test_snlm_resonant_reproduction.py`; possibly `data/discover.py` integration.

- [ ] **Step 1:** Parametrise `construct_resonant_cycler` over the catalogue rows that carry sourced `(a, e)` AND sourced per-body V∞ (the constructible SnLm subset). Assert each row's computed V∞ matches its sourced anchors within tolerance — a batch golden cross-check turning anchored-but-inert rows into validated reproductions.
- [ ] **Step 2:** Record how many rows pass (V0-construction reproduced) vs need real-eph correction. Update `data/OUTSTANDING.md` with the new reproduced count.
- [ ] **Step 3:** Commit.

## Task 5: Full-suite + lint + type gate

- [ ] `uv run pytest -q` (+ `-m slow` for the real-eph gates); `uv run ruff check . && uv run ruff format --check . && uv run mypy src`. Commit any fixups.

---

## Self-Review

**Spec coverage:** Task 1 (elements→state), Task 2 (resonance-anchored construction + the S1L1 golden cross-check — likely the idealised win), Task 3 (in-family real-eph closure — the S1L1 gate flip if it lands), Task 4 (catalogue-scale reproduction), Task 5 (gates).

**Provenance:** (a,e) and V∞ are independently sourced, so constructing from one and checking the other is non-circular. Derived epochs/ToFs are COMPUTED-labelled. No tolerance loosening; xfail stays honest if a gate doesn't genuinely pass.

**Risk:** Task 2's cross-check could *fail* if the published a/e and V∞ come from different model fidelities (circular-coplanar vs real). That is itself a valuable finding (record, don't force). Task 3 may narrow but not close S1L1 on DE440 — the documented next lever is a Newton differential-corrector (Russell 2004 method), a follow-on milestone.

---

## Findings from execution (2026-06-03)

**Tasks 1–2 DONE & pushed** (`coe_to_rv`, `construct_resonant_cycler`; commit `35dff13`). The resonance‑anchored construction reproduces, from each cycler's *sourced* orbit, with no optimisation:
- **S1L1** (a=1.30, e=0.257) → V∞ 4.90 / 4.98 km/s, E→M leg 152.6 d — matches Russell 2004 coplanar (4.99/5.10) & McConaghy 2006 (4.7/5.0) and the sourced ~154 d. ✅
- **Aldrin** (out/in) → V∞ 6.58 / 9.75 km/s — matches sourced 6.5 / 9.7. ✅

**Task 3 (spec 5.65/3.05) re‑scoped:** these are a *higher‑fidelity* figure, not S1L1's coplanar signature. The Mars 3.05 specifically requires **eccentric Mars** (coplanar geometry forces ~5.0). So closing to 5.65/3.05 is a real‑ephemeris constrained search over Mars's orbital phase, seeded by the construction (which already lands E→M at ~154 d) — not a coplanar task.

**Task 4 (catalogue scale) is DATA‑LIMITED:** only **3 of 233** rows carry both a cycler‑level (a,e) *and* sourced V∞ (see `scripts/batch_resonant_reproduction.py`). The ~200 SnLm/Russell rows have V∞ but no cycler‑level (a,e); the peri/apo rows (VISIT, Case 1, U0L1, Hollister) lack V∞. **Unlock:** populate per‑cycler (a,e) for the ballistic Russell rows from the Russell 2004 dissertation tables (already in `docs/refs/`), then the batch script validates them. That is the next data task before Task 4 yields at scale.
