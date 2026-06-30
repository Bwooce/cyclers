# Digest: Naik, Lekien, Ross (2017) — Lobe Dynamics and Rate of Escape

**Full citation:** Naik S., Lekien F., Ross S.D., "Computational Method for Phase Space Transport
with Applications to Lobe Dynamics and Rate of Escape," *Regular and Chaotic Dynamics*, Vol. 22,
No. 3, 2017, pp. 272–297. DOI 10.1134/S1560354717030078

**Status: DIGESTED**

---

## 1. Contribution

Presents a complete numerical pipeline for computing **lobe dynamics** in 2D area-preserving maps:
(1) intersection point classification (PIPs, SIPs) via orientation function; (2) lobe boundary
identification via equivalence classes and signed adjacency; (3) lobe area via Green's theorem on
piecewise-linear curves; (4) nontransverse intersection handling via the interior winding-number
integral; (5) an adaptive curve densifier for near-tangent manifolds. Released as the open-source
**Lober** software (C + MATLAB). Applied to two test cases: Oscillating Vortex Pair (OVP) fluid
flow and ship capsize (2-DOF conservative system).

---

## 2. Algorithm Summary

### Intersection Classification

- **PIP/SIP** (primary/secondary intersection point): determined by orientation
  ρ(pi) = sgn(sin θ) ‖t1(pi)‖ · ‖t2(pi)‖  where θ = angle between tangent vectors
- Equivalence classes by signed adjacency γ(pi, pj) partition the intersection points
- Lobe boundaries = arcs between consecutive PIPs on alternating curves

### Lobe Area (Green's Theorem)

[Ai] = (1/2) Σk (yk · xk+1 − xk · yk+1)  [Eq. 2.6]

where the sum is over the piecewise-linear boundary of lobe Ai.

### Error Estimate

δ[A1 \ A2] = (1/4) |Q2 − Q3 − Q1|

where Q1, Q2, Q3 are partial sub-area sums from the three-corner decomposition.

### Nontransverse Case (Interior Function)

For near-tangent or nontransverse intersections:

Ji(x0, y0) = Im ∫_{Ci} (dx + i dy) / (x − x0 + i(y − y0))

Activated via `lober -light` option (see below).

### Transport Rate Formula

a_{2,1}(n) = μ(L_{2,1}(1)) − Σ_{m=1}^{n} μ(L_{2,1}(1) ∩ f^{m-1}(L_{1,2}(1)))  [Eq. 5.7]

T_{1,2}(n) = [B(q) \ (B(q) ∩ f^{-n}(B(q)))]  [Eq. 5.8 — boundary method]

---

## 3. Lober Software

### Repositories

- **Curve densifier (densifier tool):** https://github.com/Shibabrat/curve_densifier  (C + MATLAB)
- **Lober (main lobe solver):** https://github.com/shibabrat/lober

Note: two distinct repos; densifier is called from within Lober via -DENS option.

### Command Syntax

```
# Transverse intersections (standard):
lober <c1> <c2> <rslt> [ -DENS <nPass> <nDens> ]

# Nontransverse / boundary method:
lober -light <c1> <c2> <rslt> [ -DENS <nPass> <nDens> ]
```

`<c1>`, `<c2>`: input curve files; `<rslt>`: output file to create.

### Input Format (Tecplot ASCII)

```
VARIABLES="x""y"
ZONE T="the curve C1"
0.2    0.4
0.23   0.45
...
```

### Output Files

- `<rslt>`: one line, 4 numbers: [area_inside_lobes, area_outside_lobes, rel_err1, rel_err2]
- 6 additional Tecplot ASCII files (N×2 arrays):
  - `c10.dat`, `c20.dat` — intersection points
  - `c11.dat`, `c22.dat` — boundaries of C1∩C2 and C2∩C1
  - `c12.dat`, `c21.dat` — boundaries of C1\C2 and C2\C1

---

## 4. Sourced Goldens

### OVP Test Case

| Quantity | Value |
|---------|-------|
| Perturbation amplitude | ε = 0.1 |
| Circulation strength γ (near-orthogonal) | 0.5 |
| Circulation strength γ (near-tangent, needs -light) | 1.81 |
| PIP at q (first intersection) | (0.0, 2.065) |
| Lobe area sign change | γ ≥ 1.5 (geometry inverts) |

Validation: Lober output agrees with Rom-Kedar, Leonard, Wiggins (J. Fluid Mech. 214, 1990)
Fig. 9 brute-force calculations.

### Ship Capsize Test Case (Section 5.2)

System (rescaled Lagrangian):
- V(x,y) = (1/2)x² + y² − x²y
- Equations: ẍ = −x + 2xy, ÿ = −R²y + (1/2)R²x²
- Parameter R = ωθ/ωφ = **1.6** (pitch/roll frequency ratio)
- Poincaré SOS: U1 = {(y, vy) | x=0, vx>0}

Critical energy and saddle equilibrium points:
- Critical energy: **Ee = 0.25**
- Saddle equilibrium points: **(±1, 0.5, 0, 0)**

Escape rate computation (Table 5):
- Energy above critical: e = 0.28, Δe = 0.00307

| Iterate | Exit via left (entered left) | Exit via left (entered right) |
|---------|-----------------------------|-----------------------------|
| 1       | 0%                          | 0%                          |
| 2       | 0%                          | 11.5%                       |
| 3       | 2.93%                       | 0.016%                      |
| 4       | 1.87%                       | 1.441%                      |

| Iterate | Exit via right (entered left) | Exit via right (entered right) |
|---------|-------------------------------|-------------------------------|
| 1       | 0%                            | 0%                            |
| 2       | 0%                            | 0%                            |
| 3       | 11.2%                         | 2.90%                         |
| 4       | 0.0246%                       | 0.278%                        |

SOS energy values shown in Fig. 17: e = 0.22, 0.25, 0.28.

---

## 5. Reuse Assessment

### Transport / lobe dynamics (primary use)
**USE DIRECTLY — this IS the tool.** Lober is the reference software implementation for computing:
- Lobe areas from manifold curves (stable/unstable manifolds of hyperbolic fixed points)
- Transport rates a_{i,j}(n) across partial barriers
- Boundary method for multilobe / self-intersecting turnstiles

Both the OVP pip location q=(0.0, 2.065) and the Table 5 ship capsize percentages are
**validation goldens** for any lobe dynamics implementation. Before trusting any lobe area
computation, verify it reproduces Table 5 at e=0.28, R=1.6.

### #267 resonance_network.py
**INDIRECT (important).** Transport rates across MMR partial barriers (BIPs from Rawat 2026)
require exactly this pipeline: BIP → manifold curves → Lober → lobe area → transport rate.
The lobe area per Jacobi-C level quantifies how efficiently trajectories can leave one resonance
island and enter an adjacent one. This is the outer-loop rate scoring that #267 would use.

### #292 BCR4BP / #293 ER3BP
**INDIRECT.** In both time-periodic systems, the stroboscopic Poincaré map is a 2D area-preserving
map; the transit/non-transit manifold curves from Fitzgerald 2022 are exactly the curve inputs that
Lober needs. The pipeline: Fitzgerald 2022 generates manifold curves → Lober computes lobe areas
and transport rates for the BCP or ER3BP system.

---

## 6. Companion Papers Cited

- **Hiraiwa, Bando, Nisoli, Sato (2024):** "Designing Robust Trajectories by Lobe Dynamics in
  Low-Dimensional Hamiltonian Systems," *Phys. Rev. Research* 6, L022046.
  [Ref 12 in Rawat 2026 = Hiraiwa already in KNOWN_CORPUS; confirms the lobe-design connection]
- **Naik & Ross (2017):** "Geometry of Escaping Dynamics in Nonlinear Ship Motion,"
  *Commun. Nonlinear Sci. Numer. Simul.* 47, 48–70. [Companion paper with full ship capsize analysis]
- **Ross & Scheeres (2007):** Multiple gravity assists and capture/escape in the CR3BP,
  *SIAM J. Appl. Dyn. Syst.* 6(3), 576–596. [Ref 10 in Rawat = lobe dynamics in CR3BP context]

---

## 7. Mined vs Digested

**DIGESTED.** No catalogue records affected (purely methodological, no cyclers). CORPUS_INDEX
entry needed for both this paper and the companion Naik-Ross 2017.
