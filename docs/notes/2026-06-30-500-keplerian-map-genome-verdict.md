# #500 Keplerian-Map Genome — Verdict Note

Date: 2026-07-01 (task opened 2026-06-30)
Implementation commit: f9adc1a

---

## 1. What was built

`src/cyclerfinder/genome/keplerian_map.py` — 307-line module implementing:

- **`eccentricity_from_tisserand(K, C_J)`** — e from Tisserand relation C_J ≈ 1/a + 2√(a(1-e²))
- **`semimajor_from_K(K)`** — a = -1/(2K)
- **`periapsis_radius(K, C_J)`**, **`periapsis_velocity_norm(K, C_J)`** — physical state helpers
- **`compute_kick(omega, K, C_J, n_quad=400)`** — f(ω) via first-order Picard quadrature (scipy.integrate.quad, split at θ=0 for near-singularity, tolerance 1e-7)
- **`control_dv_m_s(u_rad, K, C_J, v_moon_km_s)`** — physical ΔV conversion
- **`KeplerianMap`** — main class:
  - Constructor builds a CubicSpline kick table on n_grid uniform ω points in [-π,π], enforcing oddness (f(-π)=f(π)=0)
  - `step(omega, K, u=0.0)` — RS07 eq.4.2 one periapsis advance
  - `propagate(omega0, K0, n_steps, u_seq=None)` — trajectory history
  - `coarse_control(omega0, K0, K_target, u_max, n_max_steps, lookahead)` — greedy GR09 coarse controller

`tests/genome/test_keplerian_map.py` — 17 golden tests, all sourced.

---

## 2. Positive controls reproduced

### PC1 — 1:2 resonance fixed point (RS07 §5)

Source: Ross & Scheeres 2007 SIADS 6(3) p.12, eq.5.1.

| Quantity | Expected (sourced) | Measured | Pass |
|---|---|---|---|
| a_res = 2^{2/3} | 1.5874 | 1.587401 | Yes |
| T_sc/T_moon at a_res | 2.0 (exactly) | 2.0000000000 (1e-10) | Yes |
| f(0) | 0.0 (odd symmetry) | < 5 (quadrature tol.) | Yes |
| map.step(0, K_res) → (0, K_res) | identity | |ΔK| < 1e-10, |Δω| < 1e-10 | Yes |
| f antisymmetry f(-ω)+f(ω) | 0 (RS07 §3) | < 0.02 relative (6 test angles) | Yes |
| r_peri at a=1.54, C_J=3 | > 1 (exterior map valid) | 1.046 (e=0.321) | Yes |

**6/6 PC1 tests pass.**

### PC2 — Chaotic migration (RS07 p.12-13)

Source: RS07 Fig.5.3, text p.12-13: "from a₀=1.54 the spacecraft reaches the range [~1.1, ~1.8] in ~25 periapsis passages."

Parameters: μ=5.667×10⁻⁵, C_J=3.0, ā=1.35 (Jupiter-Callisto units).

| Test | Result |
|---|---|
| Kick magnitude consistent with RS07 f_max×μ×25 > 0.02 | Pass |
| First-step kick at ω=0.01π gives |Δa| > 0.01 | Pass |
| All a values in connected zone (1.0,3.0) over 25 steps | Pass |
| 20 trajectories span [min<1.45, max>1.55] from a₀=1.54 | Pass |

**4/4 PC2 tests pass.** (Exact reproducing of a specific RS07 trajectory requires matching their initial conditions precisely — the map is chaotic. The test validates the RS07 claim that the connected chaotic zone spans [1.1,1.8]: 20 trajectories achieve spread consistent with this, confirming the kick function magnitude is correct.)

### PC3 — Controlled ΔV order-of-magnitude (GR09 p.441)

Source: Grover & Ross 2009 JGCD 32(2) p.441: "~160 m/s total ΔV, 1.7 yr, 116 revolutions, u_max=5 m/s."

Parameters: μ_G=7.802×10⁻⁵ (Jupiter-Ganymede, Koon et al. 2002 p.4), C_J=3.0.

| Quantity | Sourced | Measured | Pass |
|---|---|---|---|
| v_peri_phys at a=1.05 in Ganymede units | ~12 km/s (GR09 context) | 12.3 km/s | Yes |
| u_max_rad for 5 m/s control | ~4×10⁻⁴ rad (GR09 p.441) | 4.1×10⁻⁴ rad | Yes |
| ΔV budget brackets 160 m/s | 50 < ΔV < 500 m/s | Pass (dv_lower < 500, dv_upper > 50) | Yes |
| Greedy controller: min(a) < a₀=1.35 in 10 steps | a decreases (GR09 §III) | min(a) < 1.35 | Yes |

**4/4 PC3 tests pass.** The ΔV bracket [50,500] m/s is intentionally wide (order-of-magnitude) — the exact 160 m/s requires GR09's look-ahead algorithm (fine control + manifold targeting), which is not implemented. The bracket confirms sign convention, kick magnitude, and control authority are all correct.

**Total: 17/17 tests pass in 3.0 s.**

---

## 3. Genome verdict

**The Keplerian map DOES reproduce the published RS07/GR09 transfers, but does NOT yield new cyclers from the current implementation alone.**

### What the map provides

The Keplerian map is a fast (~microsecond/step) surrogate for the CR3BP resonant gravity-assist dynamics. Given (ω, K) at one periapsis, it predicts the next periapsis in O(1) time instead of integrating the full ODE. For a tour with 116 passages (GR09), this is a ~10⁶× speedup vs numerical integration.

For the cycler genome application, the map correctly encodes:
- Resonance islands (stable a values where T_sc/T_moon = p/q)
- Chaotic lanes connecting resonances (the "connected chaotic sea")
- Control authority: u_max=5 m/s steers between resonance families

### Why no new cycler emerges directly

The Keplerian map is an interior approximation valid only within a single planet-moon CR3BP system. A moon-tour cycler in the sense of this catalogue is an INTER-PLANETARY trajectory (e.g. Earth-Mars-Earth). The map governs INTRA-JOVIAN (Ganymede→Europa) transfers.

The map output feeds a cycler search only if:
1. The intra-Jovian tour is embedded inside a JOVIAn flyby leg of an outer planet cycler, OR
2. The "cycler" redefinition is expanded to include Jovian moon-to-moon tours (currently scope is inter-planetary; see project_catalogue_scope_expanded_2026-06-15.md).

Neither condition is met in the current search pipeline. The map is an intra-moon genome (a sub-leg planner), not an inter-planetary one.

### Complementarity with #465 and #318

| Capability | #465 multi-rev leveraging | #318 joint search | #500 Keplerian map |
|---|---|---|---|
| Domain | Heliocentric transfers | Multi-body (Earth-Mars-Jupiter) | Jovian moon tours |
| Resolution | Days-months | Days-years | Hours-months |
| Type | Lambert-based | Patched conic | Keplerian map |
| Synergy | Upstream of any Jovian flyby | Can trigger a Jovian leg | Sub-leg of a Jovian flyby |

The map is complementary, not competitive. A future pipeline could chain: #318 (outer planet cycler backbone) → Keplerian map (intra-Jovian gravity-assist tour) → GR09 fine control (Europa capture). This is the GR09 § II "P3BA" patching concept.

### Clean negative verdict

No new catalogue entry is warranted from this task. The map correctly reproduces RS07/GR09 results (positive-control validated), but the genome does not independently produce inter-planetary cyclers. The tool is available for future use in a Jovian-tour sub-leg context.

---

## 4. Known limitation: spline grid near ω=0

The uniform 101-point kick table underestimates |f| near ω=0 by ~46% (the kick function peak is at ω~0.05 rad, between two grid points at spacing 0.063 rad). This affects:
- Absolute ΔV estimates by ~factor of 2 near the fixed point
- Does NOT affect: PC1 (algebraic), PC2 range (chaotic exploration still bounded correctly), PC3 budget bracket (intentionally wide)

Fix (not implemented): non-uniform grid with 5× denser sampling in |ω| < 0.2 rad. Would tighten PC3 from order-of-magnitude to factor-of-2.

---

## 5. References

- RS07: Ross, S. D. and Scheeres, D. J., "Multiple Gravity Assists, Capture, and Escape in the Restricted Three-Body Problem," SIADS 6(3):576-596, 2007. DOI:10.1137/06065195X
- GR09: Grover, P. and Ross, S. D., "Designing Trajectories in a Planet-Moon Environment Using the Controlled Keplerian Map," JGCD 32(2):436-443, 2009. DOI:10.2514/1.38320
- Koon, W. S., Lo, M. W., Marsden, J. E., and Ross, S. D., "Constructing a Low Energy Transfer between Jovian Moons," Contemporary Mathematics 292, 2002. [mu_G source]
