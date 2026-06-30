# Gómez–Koon–Lo–Marsden–Masdemont–Ross 2004 — Spatial CR3BP manifolds (digest)

Digested 2026-06-30 (Task #499 heteroclinic/transport pass).

**Source (cite exactly):**
Gómez, G., Koon, W. S., Lo, M. W., Marsden, J. E., Masdemont, J. and Ross, S. D.,
"Connecting orbits and invariant manifolds in the spatial restricted three-body
problem," *Nonlinearity* **17** (2004) 1571–1606.
DOI: 10.1088/0951-7715/17/5/002. AAS preprint: AAS 01-301.

Full text-layer PDF read. 36 pages; all sections including Appendix (15th-order
normal-form computation). Appendix gives complete NHIM/manifold algorithm detail.

---

## 1. One-line contribution

Extends the planar CR3BP manifold-tube/heteroclinic-connection framework (Koon
et al. 2000 *Chaos*) to the **spatial (3D) case**: constructs 3D invariant manifold
tubes (topologically S³×ℝ), proves their role as transit/non-transit separatrices,
builds a 3D Petit Grand Tour (Ganymede→Europa with high-inclination Europa capture),
and finds heteroclinic connections between L1 and L2 libration-point orbits (Lissajous,
quasi-halo) in the Sun–Earth system with zero ΔV.

---

## 2. Model and method

**Model:** Spatial circular restricted three-body problem (CR3BP). Mass ratio µ is
the only parameter. Equations of motion in rotating frame (Eq.1); Jacobi integral C.
The planar system (z = ż = 0) is an invariant subsystem.

**Energy surface:** 5D manifold (six phase-space coords minus one integral). Manifold
tubes are S³×ℝ inside the 5D energy surface. Poincaré section at x = 1−µ (hyperplane
through the moon) gives a 3D intersection S³ in ℝ⁴.

**Key structural result:** Linear flow near any collinear libration point is
saddle×centre×centre (eigenvalues ±λ, ±iν, ±iω). The NHIM (normally hyperbolic
invariant manifold) is topologically S³ for each energy level h; its stable and
unstable manifolds are S³×ℝ tubes. Transit orbits lie *inside* the tubes; non-transit
orbits lie *outside*.

**Algorithm for spatial itinerary construction (§5):**
1. Compute 15th-order normal-form expansion near L₁ and L₂ (Birkhoff/Lie-series, see
   Appendix). Initial conditions for manifold tube computed in normal-form coordinates
   (q₁,p₁,q₂,p₂,q₃,p₃) then transformed to CR3BP coordinates.
2. Fix q₁ = p₁ = ±ε (small) to select transit orbits; scan circles in (q₃,p₃) plane
   (parametrized by r_v = √(q₃²+p₃²)); solve for (q₂,p₂) on the energy surface.
   About 10⁶ points used per computation.
3. Integrate to the Poincaré section x = 1−µ. The 3D object C₁^{±m,j} (intersection
   of manifold tube with the section) has topology S³ in ℝ⁴ = (y,ẏ,z,ż).
4. **The spatial reduction trick:** for a desired z-amplitude value (z',ż'), the
   intersection of C₁^{u,2} and C₁^{s,1} at that z-slice produces closed curves γ^j_{z'ż'}
   in the (y,ẏ) plane (analogues of the planar Poincaré-cut circles). Intersection of
   their interiors = transit orbits with the prescribed itinerary and z-amplitude.
5. For the Petit Grand Tour, match the Jupiter–Ganymede L₁ unstable manifold tube with
   the Jupiter–Europa L₂ stable manifold tube at a common Poincaré section; a ΔV
   bridges the energy gap between the two three-body systems.

**Heteroclinic L1↔L2 connections (§8):** Fix Jacobi constant C; scan L₁ unstable
manifold forward and L₂ stable manifold backward to the section x = 1−µ, ẋ > 0.
I₁⁻ (one crossing) is empty; I₂⁻ (two crossings) is not empty → connections exist
for ≥2 crossings. Connections found between Lissajous–Lissajous, Lissajous–quasi-halo,
and halo–halo pairs.

---

## 3. Sourced goldens (exact numbers from the paper)

**Mass parameters:**
- Jupiter–Ganymede: µ_G = 7.802 × 10⁻⁵ (p.1577, text)
- Jupiter–Europa: µ_E = 2.523 × 10⁻⁵ (p.1577, text)
- Sun–Earth: µ = 3.040 423 398 444 176 × 10⁻⁶ (p.1592, §7)

**Orbital eccentricities (Jovian moons, p.1577):**
- Ganymede: 0.0006
- Europa: 0.0101
- Orbital plane inclination between Ganymede and Europa: within 0.3°

**3D Petit Grand Tour ΔV budget (§6, p.1591–1592):**
- Patch-point manoeuvre (Ganymede→Europa transfer): ΔV = **1214 m/s**
- Final insertion into high-inclination circular Europa orbit: ΔV = **446 m/s**
- **Total ΔV = 1660 m/s**
- Comparison: Hohmann (patched two-body) Ganymede→Europa: ΔV = **2822 m/s**
- Savings: 1660/2822 = **43% of Hohmann value** (57% savings)
- Transfer flight time: **~25 days** (Ganymede→Europa leg)
- Final Europa orbit: **48.6° inclination**, **100 km altitude**
- Ganymede close approach altitude: **100 km**
- Trajectory: begins on Jovicentric orbit beyond Ganymede; one loop around Ganymede;
  then Europa capture

**Patch-point z-range (§6, Fig.12, p.1591):**
- (z', ż') in range **(0.0160 ± 0.0008, ±0.0008)** [normalized CR3BP units]
- Approximately **1000 km in z position** and **20 m/s in ż velocity**

**Jacobi constants for Sun–Earth L1/L2 Poincaré maps (§7, Fig.13–14, p.1592–1593):**
- C = 3.000 85
- C = 3.000 826 459 043 28
- C = 3.000 802 915 133 64  (heteroclinic connections found at this level; halo
  z-amplitude 0.2 normalized units per notation of Gómez–Masdemont–Simó 1998)
- C = 3.000 785 158 376 34

**Example (X;M,I) transit orbit in Jupiter–Europa system (§5, p.1588–1589):**
- Jacobi constant: C = **3.0028**
- z-slice for initial construction: **(z', ż') = (0.0035, 0)** (Fig.10)

**Normal-form order:** N = 15 (Appendix, p.1398)

**Linearized frequencies at collinear libration points (Appendix, p.1602):**
λ² = (c₂ − 2 + √(9c₂² − 8c₂))/2 (saddle exponent)
ν² = (2 − c₂ + √(9c₂² − 8c₂))/2 (in-plane centre frequency)
ω² = c₂ (out-of-plane centre frequency)
where c₂ = µ|k−1+µ|⁻³ + (1−µ)|k+µ|⁻³.

---

## 4. Reuse verdicts

### #291 — 3D periodic orbit families (vertical Lyapunov, halo, quasi-halo)

**DIRECT MATCH.** §7 (pp.1592–1596) maps exactly the structure needed: how
the Poincaré section on z=0 reveals the full zoo of 3D libration orbits (central
Lissajous → quasi-halo → halo) at each Jacobi constant. Four specific C values
with qualitatively different pictures (Fig.13/14). The normal-form algorithm (Appendix)
is the production tool for seeding 3D orbit families. The 15th-order Lie-series
implementation (Fortran, §App, p.1398) is cited to Jorba–Masdemont 1999 (*Physica D*
132 189–213) which is the canonical Barcelona-school reference.

**Import:** (a) The c₂-based frequency formulae (§App eq.A1) seed initial conditions
for 3D families from the mass ratio alone. (b) The Poincaré-section topology (fixed
point = near-vertical Lyapunov; halo = two fixed points closer to boundary; Lissajous =
invariant curves around vertical Lyapunov fixed point) is the diagnostic to classify
which family a numerically continued orbit belongs to. (c) The specific C values give
concrete test cases for a corrector validation.

### #306 — 3D manifold tubes and transit orbit geometry

**DIRECT MATCH.** The paper is the canonical 3D-tube reference. The manifold
topology (S³×ℝ inside 5D energy surface) and the transit/non-transit separatrix role
are both proved and numerically demonstrated. The z-slicing algorithm (§5) is exactly
the computational method needed to find spatial transit orbits at prescribed
out-of-plane amplitude. The tube-intersection criterion (int(γ¹_{z'ż'}) ∩ int(γ²_{z'ż'}))
is the 3D analogue of the planar lemon-shaped intersection.

**Import:** The z-slicing trick collapses the hard 4D intersection (S³ ∩ S³) to a
1D family of 2D intersections (closed curves in (y,ẏ) plane at each z-slice). This
is the practical algorithm for our #306 spatial manifold-search code.

### #314 — Cross-system heteroclinic connections (E-M, J-G-E)

**PARTIAL MATCH.** §6 gives the patched-three-body Ganymede→Europa construction in
full detail (§6, Fig.12, ΔV = 1214 m/s at patch). The method — seek intersection of
Gan W⁺ᵤ(S¹) and Eur W⁺ₛ(S²) at a common Poincaré section, then apply a ΔV to bridge
the energy mismatch — is the template for any cross-system heteroclinic in a 3D
patched-three-body model. The Sun–Earth L1→L2 zero-ΔV connections (§8) demonstrate
*intra-system* heteroclinics that require no bridge ΔV.

**Import for #314 (E-M cross-system):** The cross-system case always requires a ΔV
at the patch point (the two subsystems have different Jacobi constants). The Sun–Earth
intra-system case shows zero-ΔV is possible only within a single three-body system.
For an Earth–Moon cycler context: a low-energy Earth-to-Moon transfer via Sun–Earth
L1/L2 manifolds interacting with Earth–Moon L1/L2 manifolds needs the cross-system
patching with ΔV (referenced in §Future Work p.1214 as the "Belbruno–Miller" case).

### #411 — 3D spatial extension of planar heteroclinic search (#494)

**DIRECT MATCH.** This paper is *the* 3D lift of the planar work. The algorithm in
§5 (z-slicing of 4D Poincaré cuts) is the exact generalization to implement for #411.
The key insight: instead of looking for S¹ ∩ S¹ curves in the planar case, we look
for S³ ∩ S³ in 4D — but reduce it to a 1-parameter family of S¹ ∩ S¹ by slicing at
(z',ż'). The practical implementation (scatter ≈10⁶ points on the normal-form tubes,
integrate to section, then post-process the (z,ż) and (y,ẏ) projections) is described
in §5. Heteroclinic examples between Lissajous, halo, quasi-halo orbits shown in §8
(Figs.17–19).

### #308 — asteroid transport / quasi-Hilda leveraging

**NOT DIRECTLY APPLICABLE.** The paper does not address transport rates or
set-oriented methods. However, §1 and the general manifold-tube framework are cited
background for the Dellnitz 2005 paper (see that digest). The Petit Grand Tour concept
(sequential moon capture) is distinct from asteroid transport.

---

## 5. Key citations (from paper's reference list)

- [26] = Koon–Lo–Marsden–Ross 1999 AAS 99-451 (Genesis/heteroclinic, digested in this
  same batch)
- [27] = Koon et al. 2000 *Chaos* 10, 427–469 (planar CR3BP heteroclinic, canonical)
- [23] = Jorba–Masdemont 1999 *Physica D* 132, 189–213 (centre-manifold algorithm,
  the primary software reference)
- [30] = Koon et al. 2002 *Contemp. Math.* 292, 129–145 (planar Petit Grand Tour)
- [13] = Gómez–Koon–Lo–Marsden–Masdemont–Ross 2001 AAS 01-301 (3D spatial preview)

---

## 6. Status

**Digested** (full paper read, all sourced numbers extracted). Not yet cross-referenced
against code. Do NOT add to catalogue (no cycler trajectories). CORPUS_INDEX update
deferred to the calling task's commit step.
