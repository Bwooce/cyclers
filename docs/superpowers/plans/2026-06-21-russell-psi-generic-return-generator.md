# Russell ψ Generic-Return Generator + Global Cycler Search — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build Russell's (2004) circular-coplanar generic-return generator, his turn-angle-min-max cycler assembly, and the global `p.h.s.i` cycler search, validated against his Tables 2.2/2.3/3.4, then wire it to assemble the catalogue's SnLm rows.

**Architecture:** Three new pure modules under `src/cyclerfinder/search/` — `generic_return.py` (the numerical multi-rev-Lambert return generator), `cycler_assembly.py` (full/half-rev geometry + turn-angle min-max), `cycler_search.py` (Eq 3.1 + Fig 3.9 global loop) — plus a `RussellModel` constants object and a closer/batch driver. Reuses `core/lambert.lambert(max_revs=…)` and `core/kepler.propagate`.

**Tech Stack:** Python 3, numpy, scipy.optimize, pytest, uv. `from __future__ import annotations`, full type hints, ruff+mypy clean (pre-commit enforces).

**Source of truth:** `docs/superpowers/specs/2026-06-21-russell-psi-generic-return-generator-design.md` — **all equations are verbatim in its Appendix A** (A.1 assembly/Eq 3.1; A.2 generation; A.3 Eqs 2.12–2.25). Cite "spec App A.x" rather than re-deriving.

**Model constants (LOAD-BEARING):** Earth period 1.0 yr; **Mars 1.875 yr**; synodic τ = 1/(1/1−1/1.875) = 15/7 yr (compute, don't hardcode 2.14); canonical μ_sun=1 AU³/TU², 1 TU = 58.1324409 d, 1 AU = 149597871 km, μ_earth = 3.00348960e-6. Goldens only match in this model.

---

## File structure

- Create `src/cyclerfinder/search/generic_return.py` — `RussellModel`, `GenericReturn`, ψ geometry, grid generator, sub-family binning, target-|v∞| interpolate+refine.
- Create `src/cyclerfinder/search/cycler_assembly.py` — full/half-rev return geometry (Eqs 2.12–2.25), `f_count` (Table 3.2), `omega_minimax` (Eqs 3.2–3.8), `group_half_years` (§3.6), `descriptor_to_phsi`, `assemble_cycler`.
- Create `src/cyclerfinder/search/cycler_search.py` — `cycler_tof` (Eq 3.1), `Cycler`, AR/TR, `search_cyclers` (Fig 3.9).
- Create `scripts/russell_cycler_batch.py` — assemble + report the catalogue SnLm rows (no writeback).
- Tests: `tests/search/test_generic_return.py`, `tests/search/test_cycler_assembly.py`, `tests/search/test_cycler_search.py`.

---

## PHASE A — Generic-return generator

### Task A1: RussellModel constants

**Files:** Create `src/cyclerfinder/search/generic_return.py`; Test `tests/search/test_generic_return.py`

- [ ] **Step 1: Failing test**
```python
from cyclerfinder.search.generic_return import RussellModel
import math
def test_russell_model_constants():
    m = RussellModel()
    assert m.tu_days == 58.1324409
    assert m.au_km == 149597871.0
    assert m.mu_sun == 1.0
    assert m.period_yr("E") == 1.0
    assert m.period_yr("M") == 1.875
    # synodic = 1/(1/1 - 1/1.875) = 15/7
    assert abs(m.synodic_yr("E", "M") - 15.0/7.0) < 1e-12
    # circular speed of Earth in canonical AU/TU: sqrt(mu/a), a_E from period via Kepler III
    aE = (1.0**2) ** (1/3)  # P^2 = a^3 in canonical (P in yr? no—use TU). See impl note.
    assert m.body_circular_speed("E") > 0.0
```
- [ ] **Step 2: Run → FAIL** (`ModuleNotFoundError`). `uv run pytest tests/search/test_generic_return.py -v`
- [ ] **Step 3: Implement `RussellModel`** — a frozen dataclass holding the canonical units + per-body period (E=1.0, M=1.875 yr). Methods: `period_yr(body)`; `synodic_yr(a,b)=1/abs(1/period_yr(a)-1/period_yr(b))`; semi-major axis from period via Kepler III in canonical units (year↔TU conversion using `tu_days`, `YEAR_DAYS=365.25`); `body_circular_speed(body)=sqrt(mu_sun/a_canonical)`; `body_state(body, theta)` returning coplanar circular `(r_vec, v_vec)` in canonical AU, AU/TU at mean-anomaly angle θ. Keep Earth/Mars only for now (extensible map).
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `git add -A && git commit -m "search/#388: RussellModel circular-coplanar constants (Mars 1.875 yr, canonical units)"`

### Task A2: ψ geometry + GenericReturn

**Files:** `generic_return.py`; test same file
- [ ] **Step 1: Failing test** — ψ is the in-plane angle of the v∞ vector referenced to v_B, positive toward r_B (spec App A.2):
```python
import numpy as np
from cyclerfinder.search.generic_return import RussellModel, psi_of_vinf_vec, GenericReturn
def test_psi_reference_geometry():
    m = RussellModel()
    r_B, v_B = m.body_state("E", 0.0)
    # a v∞ aligned with +v_B has psi 0; rotating toward r_B increases psi
    vinf_along_v = v_B / np.linalg.norm(v_B) * 0.1
    assert abs(psi_of_vinf_vec(vinf_along_v, r_B, v_B)) < 1e-9
    # 90 deg toward r_B:
    rhat = r_B / np.linalg.norm(r_B)
    vinf_along_r = rhat * 0.1
    assert abs(abs(psi_of_vinf_vec(vinf_along_r, r_B, v_B)) - np.pi/2) < 1e-6
def test_generic_return_dataclass():
    g = GenericReturn(psi_deg=114.0, tof_body_periods=1.25, a_au=0.804, n_revs=1, branch="slow", vinf=0.5)
    assert g.branch == "slow"
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement** `psi_of_vinf_vec(vinf_vec, r_B, v_B)` — build the in-plane basis (v̂_B, and the in-plane perpendicular toward r̂_B), return `atan2(component toward r_B, component along v_B)`. `GenericReturn` frozen dataclass with the fields in the test.
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: psi referencing-angle geometry + GenericReturn`

### Task A3: Multi-rev Lambert grid generator

**Files:** `generic_return.py`; test same file
- [ ] **Step 1: Failing test** — coarse grid produces binnable solutions of multiple (N, branch):
```python
from cyclerfinder.search.generic_return import RussellModel, generate_generic_returns
def test_generate_returns_coarse_grid():
    m = RussellModel()
    rs = generate_generic_returns(m, "E", max_tof_body_periods=6.0, dtheta_deg=2.0, max_revs_cap=4)
    assert len(rs) > 50
    assert {g.n_revs for g in rs} & {1, 2}            # multiple rev counts present
    assert {g.branch for g in rs} == {"fast", "slow"} # both branches present
    assert all(g.vinf > 0 for g in rs)
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement** `generate_generic_returns(model, body, *, max_tof_body_periods=6.0, dtheta_deg=0.5, refine_dtheta_deg=1/24, max_revs_cap=15)` per spec App A.2: step ToF over the grid (ToF↔transfer-angle via the body's mean motion); at each ToF solve `core/lambert.lambert(r_B(t0), r_B(t0+ToF), tof, mu=model.mu_sun, prograde=True, max_revs=max_revs_cap)`; for each returned solution record departure v∞ = `v1 − v_B(t0)`, its `|v∞|`, ψ via `psi_of_vinf_vec`, `a` from vis-viva, `n_revs=sol.n_revs`, `branch` from `sol.branch` (map "low"/"high"→"slow"/"fast" per spec App A.2: '+'=slow,'−'=fast). Posigrade only. Skip degenerate/`LambertError`.
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: multi-rev Lambert generic-return grid generator`

### Task A4: Sub-family binning + interpolate + refine

**Files:** `generic_return.py`; test same file
- [ ] **Step 1: Failing test**
```python
from cyclerfinder.search.generic_return import RussellModel, generate_generic_returns, bin_sub_families, returns_at_vinf
def test_bin_and_query_at_vinf():
    m = RussellModel()
    rs = generate_generic_returns(m, "E", dtheta_deg=2.0, max_revs_cap=4)
    bins = bin_sub_families(rs)
    assert all(isinstance(k, tuple) and len(k) == 2 for k in bins)  # (n_revs, branch)
    got = returns_at_vinf(m, "E", 0.5, dtheta_deg=2.0, max_revs_cap=4)
    assert got  # at least one sub-family brackets |v∞|=0.5
    assert all(abs(g.vinf - 0.5) < 1e-3 for g in got)  # refined to target |v∞|
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement** `bin_sub_families(returns)` → `dict[(n_revs,branch), sorted-by-vinf list]`; `returns_at_vinf(model, body, vinf, ...)` — per sub-family find the bracketing pair around `vinf`, linearly interpolate (ψ, ToF, a), then 1-D refine: fix |v∞|, Newton on (ToF, ψ) — propagate the departure state (built from ψ, |v∞|) for ToF via `core/kepler.propagate`, residual = body-position miss at arrival; iterate to `< 1e-9` (machine precision per spec App A.2). Return converged `GenericReturn`s (vinf exactly the target).
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: sub-family binning + target-vinf interpolate + 1-D refine`

### Task A5: GOLDEN — reproduce Russell Tables 2.2 & 2.3

**Files:** Test `tests/search/test_generic_return.py`
- [ ] **Step 1: Write the golden test** (EXPECTED = dissertation, spec App A.2 / mining note §2.2). Table 2.2 is at |v∞| = ½ Earth circular speed; Table 2.3 at |v∞| = 0.1838 AU/TU. Assert representative sourced rows match (ψ°, ToF body-periods, a AU, signed N) within tol (ψ ±1°, ToF ±0.05, a ±0.03):
```python
import pytest
from cyclerfinder.search.generic_return import RussellModel, returns_at_vinf
def _find(rs, n_revs, branch):
    cand = [g for g in rs if g.n_revs == n_revs and g.branch == branch]
    return cand
@pytest.mark.slow
def test_golden_table_2_2_N1_returns():
    m = RussellModel()
    vinf = 0.5 * m.body_circular_speed("E")   # Table 2.2: |v∞| = half body circular velocity
    rs = returns_at_vinf(m, "E", vinf, dtheta_deg=0.5, max_revs_cap=6)
    # Sourced Table 2.2 rows (psi_deg, tof_body_periods, a_au, N): e.g.
    #   (114.0, 1.250, 0.804, +1), (-139.5, 0.899, 0.540, +1), (-73.19, 3.336, 3.114, -1)
    psis = sorted(round(g.psi_deg, 1) for g in rs)
    assert any(abs(p - 114.0) <= 1.0 for p in psis)
    assert any(abs(p - (-139.5)) <= 1.0 for p in psis)
@pytest.mark.slow
def test_golden_table_2_3_vinf_5p5():
    m = RussellModel()
    rs = returns_at_vinf(m, "E", 0.1838, dtheta_deg=0.5, max_revs_cap=15)
    # Table 2.3 spot rows: #14 (-86.96, 1.466, 1.086, -1); #34 (104.2, 1.348, 0.921, -1)
    assert any(abs(g.psi_deg - (-86.96)) <= 1.0 and abs(g.a_au - 1.086) <= 0.05 for g in rs)
```
- [ ] **Step 2: Run → likely FAIL first** (units/sign/branch-mapping mismatches are expected; this golden is the acceptance gate). `uv run pytest tests/search/test_generic_return.py -m slow -v`
- [ ] **Step 3: Fix the generator** (NOT the test): reconcile units (AU/TU vs km/s), ψ sign convention, fast/slow→branch mapping, ToF-in-body-periods, until the sourced rows match within tol. Do NOT loosen the tolerances past (ψ ±1°, ToF ±0.05, a ±0.03) — if a row won't match, that is a real finding; stop and report rather than widen.
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: GOLDEN — reproduce Russell Tables 2.2 & 2.3 generic returns`

---

## PHASE B — Re-initiating returns + turn-angle min-max assembly

### Task B1: Full-rev / half-rev return geometry (Eqs 2.12–2.25)

**Files:** Create `src/cyclerfinder/search/cycler_assembly.py`; Test `tests/search/test_cycler_assembly.py`
- [ ] **Step 1: Failing test** — the full-rev circle z-height (Eq 2.17) and half-rev Battin components (2.18/2.19) computed from the model:
```python
import numpy as np
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.cycler_assembly import full_rev_circle_z, half_rev_components
def test_full_rev_z_and_half_rev_components():
    m = RussellModel()
    vB = m.body_circular_speed("E"); vinf = 0.1838
    # Eq 2.17: z_F = (v_F^2 - vinf^2 - vB^2)/(2 vB); pick a feasible v_F (Eq 2.13, N=1,M=1 -> v_F=vB)
    vF = vB
    z = full_rev_circle_z(vF, vinf, vB)
    assert abs(z - (vF**2 - vinf**2 - vB**2)/(2*vB)) < 1e-12
    vhr, vht = half_rev_components(m, "E", a=1.0)   # Eqs 2.18/2.19 (r1=r2=aE)
    assert vht > 0
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement per spec App A.3** (verbatim eqs): `full_rev_feasible_vF(model, body, n, M)` (Eq 2.13), `full_rev_circle_z(vF, vinf, vB)` (Eq 2.17), `half_rev_components(model, body, a)` → `(v_Hr, v_Hθ)` (Eqs 2.18/2.19 with r₁=r₂=body radius), `half_rev_intersection(model, body, a, vinf)` → the two post-flyby v∞ tips (Eqs 2.23–2.25 + K). All canonical units.
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: full/half-rev return geometry (Russell Eqs 2.12-2.25)`

### Task B2: Table 3.2 f_count + turn-angle min-max + grouping

**Files:** `cycler_assembly.py`; test same file
- [ ] **Step 1: Failing test** (spec App A.1):
```python
import numpy as np
from cyclerfinder.search.cycler_assembly import f_count, omega_minimax, group_half_years
from cyclerfinder.search.generic_return import RussellModel
def test_f_count_table_3_2():
    assert [f_count(h) for h in range(0, 9)] == [1, 2, 2, 2, 3, 4, 4, 4, 5]
def test_omega_minimax_branches():
    m = RussellModel()
    vB = m.body_circular_speed("E"); vinf = 0.1838
    vinf_in = np.array([vinf, 0.0, 0.0])      # some incoming v∞ vector (canonical)
    # f_j = 1 path: omega_c = pi - 2|phi_GR|
    w1 = omega_minimax(m, "E", vinf, vinf_in, f_j=1)
    assert 0 < w1 < np.pi
    w3 = omega_minimax(m, "E", vinf, vinf_in, f_j=3)
    assert 0 < w3 < np.pi
def test_group_half_years_sums_to_h():
    hs, omega_max = group_half_years(m=RussellModel(), body="E", vinf=0.1838,
                                     vinf_in=np.array([0.1838,0,0]), h=10, s=3)
    assert sum(hs) == 10 and len(hs) == 3
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement per spec App A.1**: `f_count(h_j)` (even: h/2+1; odd: 2·int(h/4+1)); `omega_minimax(model, body, vinf, vinf_in_vec, f_j)` (φ_FR Eq3.3, φ_GR Eq3.4, decision rule: f=1→ω_c Eq3.8; f=2→acos(sinφ_GR sinφ_FR); f>2→ω_MIN Eq3.2 if ≥ω_a Eq3.5 else ω_b via iterating Eq3.6 then Eq3.7); `group_half_years(...)` (§3.6: equal split if ω_c≥ω_minimax else pile on first; return tuple + ω_MAX).
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: Table 3.2 f_count + turn-angle min-max + half-year grouping`

### Task B3: GOLDEN — Cycler 4.9.2.-1 assembly

**Files:** test `tests/search/test_cycler_assembly.py`
- [ ] **Step 1: Golden test** (spec App A.1 worked example): for the 4.9.2.-1 generic return, h split `{9,0}`; first-return flyby turn angles ≈ `[83°,45°,45°,45°,45°,83°]`, second `≈24°` (tol ±3°). Build the generic return via `returns_at_vinf` at that cycler's |v∞|, group h=9 across s=2, compute per-flyby turn angles, assert the sequence matches.
- [ ] **Step 2: Run → FAIL likely** (`uv run pytest tests/search/test_cycler_assembly.py -m slow -v`)
- [ ] **Step 3: Fix assembly** (NOT the test) until the turn angles match within ±3°. No widening past ±3° — a mismatch is a real finding to report.
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: GOLDEN — Cycler 4.9.2.-1 turn-angle assembly`

---

## PHASE C — Global p.h.s.i search

### Task C1: cycler_tof (Eq 3.1) + Cycler + AR/TR

**Files:** Create `src/cyclerfinder/search/cycler_search.py`; Test `tests/search/test_cycler_search.py`
- [ ] **Step 1: Failing test**
```python
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.cycler_search import cycler_tof, aphelion_ratio, turn_ratio
def test_cycler_tof_eq_3_1():
    m = RussellModel()
    tau = m.synodic_yr("E", "M")
    # Eq 3.1: TOF = (tau*p - h/2)/s, years
    assert abs(cycler_tof(m, p=4, h=3, s=1) - (tau*4 - 3/2)/1) < 1e-9
def test_ratios():
    m = RussellModel()
    assert aphelion_ratio(1.64) == 1.64/1.52          # AR = max aphelion / 1.52 AU
    assert turn_ratio(max_allowable=1.0, omega_max=0.5) == 2.0
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement** `cycler_tof(model,p,h,s)` (Eq 3.1, years), `aphelion_ratio(aph_au)=aph_au/1.52`, `turn_ratio(max_allowable, omega_max)=max_allowable/omega_max` (max allowable = 200 km Earth-flyby bend, computed from `core/flyby` bend at the leg's v∞), `Cycler` dataclass (`p,h,s,i`, generic_return, turn_angles, ar, tr, vinf_e, vinf_m, sequence, epochs).
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: Eq 3.1 cycler TOF + AR/TR gates + Cycler dataclass`

### Task C2: search_cyclers (Fig 3.9 loop)

**Files:** `cycler_search.py`; test same file
- [ ] **Step 1: Failing test** — a small search returns ≥1 AR/TR-passing cycler:
```python
from cyclerfinder.search.cycler_search import search_cyclers
from cyclerfinder.search.generic_return import RussellModel
def test_search_small():
    cs = search_cyclers(RussellModel(), p_max=2, ar_min=0.9, tr_min=0.85, dtheta_deg=1.0)
    assert isinstance(cs, list)
    assert all(c.ar > 0.9 and c.tr > 0.85 for c in cs)
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement** the Fig 3.9 loop (spec App A.1): `DO p=1..p_max; h=1..5·p_max; s=1..3·p_max`: `tof=cycler_tof`; if `tof>0`, `generate_generic_returns`/`returns_at_vinf`; `DO i` over available (n_revs, branch); `group_half_years`; per-group `omega_minimax` → ω_MAX; AR from generic return aphelion, TR from ω_MAX; record `Cycler` iff `tr>tr_min and ar>ar_min`. Parameterize `dtheta_deg` for test speed.
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: global p.h.s.i cycler search (Fig 3.9 loop)`

### Task C3: GOLDEN — Russell Table 3.4 (44 cyclers)

**Files:** test `tests/search/test_cycler_search.py`
- [ ] **Step 1: Golden test** (`@pytest.mark.slow`): run `search_cyclers(p_max=6, ar_min=0.9, tr_min=0.85, dtheta_deg=0.5)`; assert it recovers the sourced anchor rows — Cycler **4.3.1.-5** (v∞_E≈3.10, v∞_M≈2.53, AR≈0.992) and the **Aldrin 2.1.1.+2** — within tol (v∞ ±0.15 km/s, AR ±0.02), and total count is in the neighbourhood of 44 (assert `>= 30` to allow grid-resolution variance; document the gap if < 44).
- [ ] **Step 2: Run → FAIL likely**
- [ ] **Step 3: Fix search** (NOT the test) until the anchor cyclers appear with matching v∞/AR. If the full 44 aren't recovered at the chosen resolution, `log` what's missing (no silent cap) and report — do not loosen the anchor tolerances.
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: GOLDEN — recover Russell Table 3.4 cyclers (4.3.1.-5, Aldrin)`

---

## PHASE D — Closer wiring to catalogue rows

### Task D1: descriptor_to_phsi

**Files:** `cycler_assembly.py`; test `tests/search/test_cycler_assembly.py`
- [ ] **Step 1: Failing test**
```python
from cyclerfinder.search.cycler_assembly import descriptor_to_phsi
from cyclerfinder.data.catalog import load_catalog
def test_descriptor_to_phsi_maps_known_row():
    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    spec = descriptor_to_phsi(row)
    assert spec is not None
    assert spec.s >= 1 and spec.vinf_e_kms > 0
def test_descriptor_to_phsi_none_for_ocampo():
    cat = load_catalog()
    ocampo = next(e for e in cat.entries if e.id.startswith("russell-ocampo"))
    assert descriptor_to_phsi(ocampo.raw) is None
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement** `descriptor_to_phsi(row)` → `PhsiSpec | None`: parse `free_return_arcs` (g/G generic returns → s and the rev count i; f/h → h half/full-rev returns) + `invariants.transit_times_days` + the sourced `vinf_kms_at_encounters` (E/M anchors); return None when no per-arc descriptor (ocampo rows). Dataclass `PhsiSpec(p,h,s,i,vinf_e_kms,vinf_m_kms,...)`.
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `search/#388: descriptor -> p.h.s.i mapping for catalogue rows`

### Task D2: assemble_cycler + batch report

**Files:** `cycler_assembly.py`; Create `scripts/russell_cycler_batch.py`
- [ ] **Step 1: Failing test** — `assemble_cycler` builds a Cycler from a PhsiSpec at the sourced |v∞|:
```python
from cyclerfinder.search.cycler_assembly import descriptor_to_phsi, assemble_cycler
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.data.catalog import load_catalog
def test_assemble_known_row():
    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    spec = descriptor_to_phsi(row)
    cyc = assemble_cycler(RussellModel(), spec)
    assert cyc.vinf_e > 0 and cyc.vinf_m > 0
```
- [ ] **Step 2: Run → FAIL**
- [ ] **Step 3: Implement** `assemble_cycler(model, phsi)` (build the generic return at the sourced |v∞| via `returns_at_vinf`, assemble re-initiating flybys via Phase B); write `scripts/russell_cycler_batch.py` iterating descriptor-bearing rows → assemble → report emerged v∞ vs sourced anchor / AR / TR / turn angles + a runlog `data/runs/russell-cycler-<ts>.jsonl`. NO catalogue writeback; print proposed promotions HELD.
- [ ] **Step 4: Run → PASS**, then run the batch: `uv run python scripts/russell_cycler_batch.py`
- [ ] **Step 5: Commit** `search/#388: assemble_cycler + catalogue-row batch (held, no writeback)`

### Task D3: Results note

**Files:** Create `docs/notes/2026-06-21-russell-cycler-assembly-results.md`
- [ ] **Step 1** Write the note: per-row outcome from the batch (assembled v∞ vs sourced anchor, AR, TR, ballistic?), which catalogue rows the generator-built assembly reproduces within the 0.5 km/s gate, and the held promotion verdict (no writeback). State plainly if rows still don't assemble — that is the honest finding.
- [ ] **Step 2: Commit** `git add docs/notes/2026-06-21-russell-cycler-assembly-results.md data/runs/ && git commit -m "research/#388: Russell cycler-assembly batch results (held)"`

---

## Final verification
- [ ] `uv run pytest tests/search/test_generic_return.py tests/search/test_cycler_assembly.py tests/search/test_cycler_search.py -q` (incl. `-m slow`) — all PASS.
- [ ] `uv run ruff check . && uv run ruff format --check .` — clean.
- [ ] No `data/catalogue.yaml` / `validate.py` edit (held-for-review honesty gate).
- [ ] Report golden status (Tables 2.2/2.3/3.4 + 4.9.2.-1) and the catalogue-row assembly verdict.
