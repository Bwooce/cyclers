# Russell ψ Generic-Return Generator + Global Cycler Search + Closer Wiring — Design Spec

**Date:** 2026-06-21
**Status:** Design approved (brainstorming), pending implementation plan.
**Tracks:** #388 / #307. **Deferred sibling:** #414 (out-of-plane 3D v∞ sphere).
**Source:** Russell, R. P., "Global Search and Optimization of Free-Fall Cycler
Trajectories," Ph.D. dissertation, UT Austin, 2004. Method mined in
`docs/notes/2026-06-07-russell-2004-dissertation-method-mining.md` plus two 2026-06-21
deep extractions (Ch.3 §3.4–3.8 assembly; Eq 3.1 + §2.7 generation) recorded in this
spec's Appendix A.

## Why

The DSM descriptor-seed lane cannot close the Russell/McConaghy SnLm cyclers even on
the circular-coplanar model (#388 spikes, `project_dsm_closure_modeljump_blocker`):
its same-body legs are plain Lambert arcs with no free-return |v∞|-preservation, so the
corrector collapses to a degenerate zero-v∞ basin. The v∞-preserving generic free
return *exists* and is near-ballistic — but it must be built **by Russell's ψ
construction**, not discovered by local least-squares. This spec builds Russell's
actual method: the generic-return generator, his cycler assembly (turn-angle min-max
flyby placement), the global `p.h.s.i` search, and wiring it to close the catalogue
rows. Validated against Russell's own sourced tables (2.2, 2.3, 3.4).

## Scope

IN: coplanar (ecliptic-plane) generic-return generation; (|N|, fast/slow) sub-family
binning + interpolation + 1-D refine; full-rev/half-rev re-initiating-return geometry;
turn-angle min-max flyby assembly; global `p.h.s.i` enumeration with AR/TR gates;
golden validation vs Tables 2.2/2.3 (returns) and 3.4 (44 cyclers); wiring to assemble
+ close the descriptor-bearing catalogue rows.

OUT: out-of-plane / 3D v∞ sphere (#414); powered-flyby continuation to real DE440
ephemeris (the existing dsm/correct lanes own real-eph; this spec is the
circular-coplanar generator Russell's parents live in); any catalogue writeback.

## Model constants (LOAD-BEARING — the goldens only match in Russell's model)

Russell's simplified model, NOT our DE440 defaults:
- Earth period 1.0 yr; **Mars period 1.875 yr** (chosen so geometry repeats every
  15 yr; NOT the real 1.881). Synodic period **τ = 1/(1/1 − 1/1.875) = 15/7 ≈
  2.142857 yr** (derive unrounded from the periods, never hardcode 2.14).
- Canonical units: μ_sun = 1 AU³/TU²; 1 TU = 58.1324409 days; 1 AU = 149597871 km;
  μ_earth = 3.00348960e-6 AU³/TU². Circular-coplanar zero-sphere-of-influence patched
  conic; Earth and Mars on circular coplanar orbits.

These live in a single `RussellModel` constants object so the generator and the golden
tests share one source of truth. The existing `Ephemeris("circular")` backend is close
but uses real periods; this spec uses Russell's 1.875-yr Mars — so the generator takes
an explicit model object rather than reusing the default ephemeris (golden fidelity).

## Architecture — five phases, each independently testable

### Phase A — Generic-return generator (`search/generic_return.py`)

The numerical core (Russell §2.7.3–2.7.5; Appendix A.2). Generic returns are NOT
analytic — they are the multi-rev Lambert grid + binning + interpolation + refine.

```python
@dataclass(frozen=True)
class GenericReturn:
    psi_deg: float           # in-plane v∞-sphere angle, ref to v_B, + toward r_B
    tof_body_periods: float  # ToF in flyby-body periods (Russell's table unit)
    a_au: float              # transfer semi-major axis
    n_revs: int              # |N| ≥ 1 (0 = single-rev)
    branch: str              # "fast" | "slow"  (Russell's −/+ sign)
    vinf: float              # |v∞| (canonical AU/TU)

def psi_of_vinf_vec(vinf_vec, r_B, v_B) -> float        # ψ from a v∞ vector (geometry)
def generate_generic_returns(model, body, *, max_tof_body_periods=6.0,
                             dtheta_deg=0.5, refine_dtheta_deg=1/24,
                             max_revs_cap=15) -> list[GenericReturn]
```
- Step the transfer over a grid (Russell: ½° transfer-angle interval, refined to 1/24°
  near intersections, ~180k solves over 6 body periods); at each, solve the
  multi-rev Lambert `r_B(t0) → r_B(t0+ToF)` (reuse `core/lambert.lambert(max_revs=…)`),
  posigrade only; record each (|N|, fast/slow) solution's departure-side
  (|v∞|, ψ, a, N). `dtheta_deg`/`max_revs_cap` are parameters (coarse default for fast
  tests; Russell's fine values for fidelity). `N_MAX` is operational: the largest rev
  count whose multi-rev minimum ToF still fits the 6-body-period cap (=15 in Russell's
  worked example) — not a closed form.

```python
def bin_sub_families(returns) -> dict[tuple[int,str], list[GenericReturn]]
def returns_at_vinf(model, body, vinf, *, ...) -> list[GenericReturn]
```
- `bin_sub_families`: group by `(n_revs, branch)` — the "same sub-family" Russell's
  interpolation caveat requires (interpolating across sub-families is invalid).
- `returns_at_vinf`: per sub-family, find the two grid points bracketing the target
  |v∞|, **linearly interpolate (ψ, ToF, a)**, then **1-D refine** — fix |v∞|, Newton on
  (ToF, ψ) via Kepler propagation until the arc re-encounters the body to tolerance
  (machine precision). Returns the converged v∞-preserving free returns BY
  CONSTRUCTION.

**Golden A:** reproduce **Table 2.2** (N=1, |v∞| = ½ body circular velocity → 11
points) and **Table 2.3** (|v∞| = 0.1838 AU/TU = 5.5 km/s → 40 points): the emerged
(ψ, ToF, a, N, fast/slow) match the tabulated cells within tolerance. EXPECTED traces
to the dissertation; never self-computed.

### Phase B — Re-initiating returns + turn-angle min-max assembly (`search/cycler_assembly.py`)

Russell §2.5–2.6 (full/half-rev return geometry) + §3.4–3.6 (turn-angle min-max);
Appendix A.1.

- `half_rev_returns` / `full_rev_returns(model, body, vinf, n_half_years)`: the analytic
  re-initiating returns (Russell Eqs 2.12–2.25, transcribed verbatim in Appendix A.3).
  Full-rev: re-encounter after integer years, identical velocity diagram, free
  longitude on the full-rev circle. Half-rev: after odd half-years, v_e flips → v∞ must
  be re-pointed.
- `f_count(h_j) -> int`: Table 3.2 (Appendix A.1) — even `h_j`: `f_j = h_j/2 + 1`; odd
  `h_j`: `f_j = 2·INT(h_j/4 + 1)`; with the tabulated inter-flyby year-spacings.
- `omega_minimax(model, vinf, vinf_in_vec, f_j) -> float`: the turn-angle min-max
  (Eqs 3.2–3.8 + decision rule, Appendix A.1): `f_j=1 → ω_c = π − 2|φ_GR|`;
  `f_j=2 → acos(sin φ_GR sin φ_FR)`; `f_j>2 → ω_MIN if ω_MIN≥ω_a else ω_b` (iterate Eq
  3.6 for λ). `φ_FR = −asin(v∞/(2 v_e))` (3.3), `φ_GR = π/2 − acos(v̂∞·v̂_e)` (3.4).
- `group_half_years(h, s, ...) -> tuple[int,...]`: §3.6 — try `h_j = INT(h/s)` equal
  split; if `ω_c ≥ ω_minimax` keep equal (remainder on the last), else pile all on the
  first (`h_1=h`, rest 0). Returns `ω_MAX` too.

**Golden B:** reproduce Cycler **4.9.2.-1**'s assembly (Appendix A.1 worked example):
`{h₁=9, h₂=0}`, first-return flyby turns `[83°, 45°, 45°, 45°, 45°, 83°]`, second `24°`
(matches Table 3.4 row). And spot-check `f_count` against Table 3.2 rows.

### Phase C — Global p.h.s.i search (`search/cycler_search.py`)

Russell §3.7 / Fig 3.9 (Appendix A.1).

- `cycler_tof(model, p, h, s) -> float`: **Eq 3.1** `TOF = (τ·p − h/2)/s` (years). A
  cycler is feasible iff `TOF > 0` AND the generic-return Lambert at that ToF has
  `N_MAX > 0`.
- `search_cyclers(model, *, p_max, ar_min=0.9, tr_min=0.85) -> list[Cycler]`: the
  Fig 3.9 loop — `DO p=1..p_max; DO h=1..5·p_max; DO s=1..3·p_max`: compute TOF; if
  feasible, generate the `2·N_MAX+1` generic returns; `DO i=−N_MAX..N_MAX`: group `h`
  across `s` (Phase B), place + optimize the re-initiating flybys (Phase B), compute
  per-flyby turn angles, **AR** (max aphelion / 1.52 AU) and **TR** (max allowable turn
  @ 200 km Earth flyby / ω_MAX); record `p.h.s.i` + properties iff `TR>tr_min AND
  AR>ar_min`.
- `@dataclass Cycler`: `p,h,s,i`, the generic return, per-flyby turn angles, AR, TR,
  v∞ at Earth/Mars, leg sequence + epochs.

**Golden C:** rediscover **Table 3.4** (44 cyclers at AR_MIN=0.9, TR_MIN=0.85): the
search returns the same `p.h.s.i` set with matching turn angles / v∞ / AR. Spot-anchor
Cycler 4.3.1.-5 (v∞_E=3.10, v∞_M=2.53, AR=0.992) and the Aldrin cycler 2.1.1.+2.

### Phase D — Closer wiring + catalogue-row closure (`search/cycler_assembly.py` + batch)

- `descriptor_to_phsi(row) -> PhsiSpec | None`: map a catalogue row's `g/G/f/h`
  arc descriptor + `invariants.transit_times_days` onto the `p.h.s.i` structure (s
  generic returns from the g/G arcs; h from the half/full-rev `f/h` arcs / transit
  gaps; i from the rev count). The descriptor IS the per-arc decomposition, so this is
  a structural map, not a search. Returns None if the row's descriptor can't be mapped
  (recorded, not crashed).
- `assemble_cycler(model, phsi) -> Cycler`: build the specific row via Phase B/C
  assembly at the row's sourced |v∞| anchors.
- Re-run a closure batch (extend/replace `scripts/dsm_closure_batch.py` or a sibling):
  for each descriptor-bearing row, assemble via this lane, report converged / emerged
  v∞ vs sourced anchor / AR / TR / turn angles. Success = the assembly reproduces the
  row's sourced v∞ within the campaign gate (0.5 km/s) with TR>1 (ballistic). NO
  writeback; promotions held for session review (orbit-closure-discipline).

## Data flow

model (Russell constants) → `generate_generic_returns` (grid+bin) →
`returns_at_vinf` (interp+refine) → [Phase B re-initiating returns + turn-angle
min-max assembly] → [Phase C p.h.s.i search OR Phase D descriptor→p.h.s.i map] →
Cycler (AR/TR-gated) → golden checks / held closure report.

## Error handling

- Sub-family with no bracketing point at target |v∞| → yields nothing (not an error).
- Degenerate Lambert at grid edges / `N_MAX=0` ToF → skipped, logged.
- `descriptor_to_phsi` returns None for un-mappable rows (recorded out-of-scope).
- Infeasible turn-angle (TR≤tr_min) → recorded as a non-ballistic / failed assembly,
  not a crash and not a tolerance-loosen.

## Testing (TDD)

- Phase A: ψ geometry unit test; `generate_generic_returns` produces binnable
  sub-families; **golden Tables 2.2 & 2.3** (coarse grid for CI speed + a `@slow`
  fine-grid run).
- Phase B: `f_count` vs Table 3.2; `omega_minimax` decision-rule branches;
  **golden 4.9.2.-1 turn angles**.
- Phase C: `cycler_tof` = Eq 3.1; **golden Table 3.4** (`@slow`; spot rows inline);
  AR/TR gate logic.
- Phase D: `descriptor_to_phsi` maps `4.991gG2` correctly; assembly reproduces a
  sourced row's v∞; un-mappable row → None; the negative (no ballistic assembly) is
  recorded, not forced.

## Honesty gates (orbit-closure-discipline — non-negotiable)

1. Golden EXPECTED = Russell's tabulated cells (2.2/2.3/3.4), never self-computed.
2. |v∞| anchors emerge from the construction; never imposed.
3. No tolerance/AR/TR loosening to force a match.
4. No catalogue writeback; any promotion held for session review.
5. Russell's model constants (Mars 1.875 yr, canonical units) used for all goldens —
   a match in the wrong model is not a match.

## Out of scope / deferred

- 3D / out-of-plane v∞ sphere → #414.
- Real-DE440 continuation of the assembled cyclers → existing dsm/correct lanes; a
  later rung, not this spec.
- (Resolved 2026-06-21) The full/half-rev geometry Eqs 2.12–2.25 are now transcribed
  verbatim in Appendix A.3 — no remaining un-sourced equations; the spec is gap-free.

## References
- `docs/notes/2026-06-07-russell-2004-dissertation-method-mining.md`
- This spec Appendix A (the 2026-06-21 Ch.3 + Eq 3.1/§2.7 extractions).
- `src/cyclerfinder/core/lambert.py` (multi-rev Lambert), `core/ephemeris.py`
  (`circular` backend), `search/dsm_descriptor_seed.py` (the lane this supersedes for
  resonant legs).
- Memory: `project_dsm_closure_modeljump_blocker`, `feedback_published_rounded_values_are_display`,
  `feedback_golden_tests_sourced_only`, `feedback_orbit_closure_discipline`,
  `project_s1l1_nomenclature`.

---

## Appendix A — Sourced method (verbatim extractions, 2026-06-21)

### A.1 Cycler assembly (Russell Ch.3 §3.4–3.8; print pp. 68–85)

**Eq 3.1 (TOF, years):** `TOF = (τ·p − h/2) / s`, τ = E-M synodic period (15/7 yr in
the 1.875-yr-Mars model). Constrains the cycle to p synodic periods.

**`p.h.s.i`:** p = synodic periods; h = half-years for full/half-rev re-initiating
returns; s = identical generic returns; i = signed rev count (sign = upper/lower
Lambert solution curve = slow/fast).

**Table 3.2 (`h_j → f_j`):** even `h_j`: `f_j = h_j/2 + 1`, all inter-flyby spacings
1 yr. Odd `h_j`: `f_j = 2·INT(h_j/4 + 1)`; spacings `t_k−t_{k−1}=1` for `k=2..f/2`,
`=½·MOD(h_j,4)` for `k=f/2+1`, `=1` for `k=f/2+2..f`. (h_j=0→f=1; 1→2; 2→2; 3→2;
4→3; 5→4; 6→4; 7→4; 8→5.)

**Turn-angle min-max (Eqs 3.2–3.8):**
- φ_FR = −asin(v∞/(2 v_e))  (3.3)
- φ_GR = π/2 − acos(v̂∞₁₋ · v̂_e)  (3.4)
- ω_MIN = acos(cos φ_FR cos φ_GR + sin φ_FR sin φ_GR)  (3.2)
- ω_a = acos(cos²φ_FR cos λ_a + sin²φ_FR), λ_a = π/(f_j−2)  (3.5)
- ω_b: solve Eq 3.6 iteratively for λ (λ_b=(π−2λ)/(f_j−2)), then
  ω_b = acos(cos φ_GR cos φ_FR cos λ + sin φ_GR sin φ_FR)  (3.7)
- ω_c = π − 2|φ_GR|  (3.8)
- Decision: f_j=1 → ω_c; f_j=2 → acos(sin φ_GR sin φ_FR); f_j>2 → ω_MIN if
  ω_MIN≥ω_a else ω_b.

**§3.6 grouping:** `Σ h_j = h`; try `h_j=INT(h/s)`; if `ω_c ≥ ω_minimax` equal-split
(remainder on last), else `h_1=h`, rest 0. ω_MAX = ω_minimax-1.

**Fig 3.9 loop + gates:** `h_MAX=5·p_MAX`, `s_MAX=3·p_MAX`; AR = max aphelion / 1.52
AU; TR = max-allowable-turn (200 km Earth flyby) / ω_MAX; record iff TR>TR_MIN AND
AR>AR_MIN. Table 3.4 used AR_MIN=0.9, TR_MIN=0.85.

**Worked example 4.9.2.-1:** `{h₁=9, h₂=0}`; first return flybys `[83°,45°,45°,45°,45°,83°]`,
second `24°`.

### A.2 Generic-return generation (Russell §2.7.3–2.7.5; print pp. 41–53)

ψ = in-ecliptic-plane angular coordinate of the v∞ solution on the constant-|v∞|
sphere, referenced to v_B, positive toward r_B. Generic returns are numerical: solve
multi-rev Lambert at ½° transfer-angle intervals (1/24° near intersections) up to 6
body periods (~180k solutions), posigrade only; bin by (|N|, fast/slow); sort each
sub-family by |v∞|; for a target |v∞| interpolate (ψ,ToF,a) between same-sub-family
neighbours; refine with a 1-D solver fixing |v∞| and iterating (ToF, ψ) in Kepler
integration to machine precision. N_MAX = largest rev count whose min ToF fits the
6-period cap (=15 in the worked example). Tables 2.2 (N=1) and 2.3 (|v∞|=0.1838 AU/TU)
are the goldens.

### A.3 Full-rev / half-rev return geometry (Russell §2.7.1–2.7.2; print pp. 32–39)

Verbatim (transcribed 2026-06-21). **Frame:** origin at the tip of v_B, z-axis along
v_B, y-axis opposed to the body's angular-momentum vector; the half-rev primed frame
has x′ along r_B. Symbols: v_B = body heliocentric speed; v∞ = hyperbolic excess; v_F =
full-rev post-flyby speed; a (a_F) = transfer semi-major axis; a_B = body semi-major
axis; r (r₁,r₂) = body distance(s) from primary; N/M = s/c-vs-body revolution counts;
γ = body flight-path angle; μ = gravitational parameter.

**Full-rev returns:**
- (2.12) vis-viva: `v_F = sqrt(2μ/r − μ/a_F)`.
- (2.13) feasible post-flyby speed: `v_F = sqrt(2μ/r − μ·(N/M)^(2/3)/a_B)` (from a_F =
  a_B·(N/M)^(2/3)). Feasible (N,M,r,μ,a_B) = those making v_F real; N = 1..N_max.
- (2.14) intersection feasibility (logical): `|v_B − v∞| ≤ v_F ≤ v_B + v∞`.
- (2.15) v∞ sphere: `x² + y² + z² = v∞²`.
- (2.16) full-rev sphere: `x² + y² + (z + v_B)² = v_F²`.
- (2.17) full-rev circle z-height: `z_F = (v_F² − v∞² − v_B²) / (2 v_B)`. The full-rev
  circle is the slice of (2.15) at z = z_F.

**Half-rev returns** (outbound terminal velocity, Battin):
- (2.18) radial: `v_Hr² = μ·[2/(r₁+r₂) − 1/a]` (negative root = fast, positive = slow).
- (2.19) transverse: `v_Hθ² = 2μ·r₂/(r₁² + r₁ r₂)` (independent of a / TOF).
- (2.20) half-rev circle (primed): `(y′)² + [z′ + v_B cos γ]² = v_Hθ²`.
- (2.21) half-rev circle plane: `x′ = v_Hr − v_B sin γ`.
- (2.22) v∞ sphere (primed): `(x′)² + (y′)² + (z′)² = v∞²`.
- (2.23) `x = [v_Hr − v_B sin γ] cos γ + K sin γ / (2 v_B cos γ)`
- (2.24) `y = ± sqrt{ v∞² − [v_Hr − v_B sin γ]² − K²/(2 v_B cos γ)² }` (± = the two
  intersection points; real-y required for existence)
- (2.25) `z = [v_Hr − v_B sin γ] sin γ − K/(2 v_B)`
- (unnumbered, defines K): `K = 2 v_Hr v_B sin γ + v_B²·[1 − 2 sin²γ] − v_Hθ² − v_Hr² + v∞²`

**Half-rev a constraint:** `a ≥ a_min` (= (r₁+r₂)/2 in the symmetric r₁=r₂ case);
below a_min, Eq 2.18's v_Hr goes imaginary. Half-rev returns exist only at the discrete
a where the ToF curve meets the required ToF.

(Numbering note: in the source 2.13 is the lone feasible-v_F formula, 2.14 the
feasibility *constraint*, 2.15 the v∞ sphere — the full-rev *circle* is its z_F slice
via 2.15+2.16+2.17. Physics matches §2.7.5's references exactly.)
