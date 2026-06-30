# Digest: Onozaki, Yoshimura, Ross (2017) — Tube Dynamics in the 4-Body System

**Full citation:** Onozaki K., Yoshimura H., Ross S.D., "Tube dynamics and low energy Earth-Moon
transfers in the 4-body system," *Advances in Space Research*, Vol. 60, No. 10, 2017, pp. 2117–2132.
DOI 10.1016/j.asr.2017.07.024

**Status: DIGESTED**

---

## 1. Contribution

Extends tube-manifold / invariant-manifold transit-orbit design from the CR3BP to the bicircular
restricted 4-body problem (BCR4BP: Sun–Earth–Moon–S/C). Key innovation: decomposes the 4-body
system into two *coupled* perturbed 3-body systems (Moon-perturbed Sun–Earth–S/C for departure;
Sun-perturbed Earth–Moon–S/C for arrival), computes FTLE fields to approximate LCS (stable/unstable
manifold analogs), and patches departure and arrival trajectories with a mid-course maneuver ΔVP.
Result: a 100-day LEO→LLO transfer totalling ΔVT = 3.880 km/s — saves 0.099 km/s vs Hohmann,
costs 0.043 km/s more than WSB but reduces flight time by 60 days.

---

## 2. Model

**Bicircular 4-body:**
- Moon-perturbed system (departure): Sun–Earth–S/C in S–BEM (Sun–Barycenter–Earth–Moon) rotating
  frame, perturbed by Moon
- Sun-perturbed system (arrival): Earth–Moon–S/C in E–M rotating frame, perturbed by Sun
- Coupling: patch points matched at Poincaré section U = {(x,y,vx,vy,hM) | x > 1−lS, y=0, vy>0}
  at the same angle: hSE_M = hEM_M

**FTLE grid:** 1000 × 1000 points, integration time T = 7 (nondimensional Earth-Moon time units).
Ridge detection uses second-derivative estimation with truncation error control (Appendix A);
tolerance TOL applied to max|εij|; central difference step Δw refined until TOL satisfied.

---

## 3. Sourced Goldens

### Physical Constants

| Quantity | Value |
|----------|-------|
| lS (Sun mass ratio)  | 3.02319×10^-6 |
| lM (Moon mass ratio) | 1.21536×10^-2 |
| aS (Earth–Sun distance) | ≈ 1.49598×10^8 km |
| aM (Earth–Moon distance) | ≈ 3.84400×10^5 km |
| ωS | ≈ 1.99640×10^-7 rad/s |
| ωM | ≈ 2.66498×10^-6 rad/s |
| mS | ≈ 1.99976×10^30 kg |
| mE | ≈ 5.97219×10^24 kg |
| mM | ≈ 7.34767×10^22 kg |

### Parking Orbits

| Orbit | Altitude | Speed |
|-------|---------|-------|
| LEO   | 167 km  | 7.80713 km/s |
| LLO   | 100 km  | 1.63346 km/s |

### Departure Energy Bounds (Moon-perturbed system, S–BEM frame)

| Bound | E^SE value | ΔVE (km/s) |
|-------|-----------|-----------|
| Minimum energy (LCS enters transit) | E^SE_Dmin = 1.50043 | ΔVE = 3.189 |
| Maximum energy used in design | E^SE_Dmax = 1.50027 | ΔVE = 3.202 |

Moon-perturbed LCS IC (repelling): wl = (0.999368, 0, −0.005, −0.0958139) in (x, y, vx, vy)

### Arrival Energy Bounds (Sun-perturbed system, E–M frame)

| Bound | E^EM value | ΔVM (km/s) |
|-------|-----------|-----------|
| Minimum (lower transit bound) | E^EM_Amin = 851.528 | ΔVM = 0.634 |
| Maximum | E^EM_Amax = 851.493 | ΔVM = 0.650 |

LLO energy (Sun-perturbed frame): E^EM_LLO = 852.703
Sun-perturbed attracting LCS IC: wl = (0.981683, 0, −0.1, −1.92901) in (x, y, vx, vy)

### Specific Transfer Design (from Section 5.3)

**Departure trajectory initial point** (S–BEM frame, nondim):
(xD, yD, vxD, vyD, hMD) = (0.999922, 0, 0, 0.368963, 0)   [E^SE_t0 = 1.50027, ΔVE = 3.202 km/s]

**Poincaré section crossing angles:**
- Departure family: hM = 0 rad
- Arrival family: hM = 2.55 rad

**Patch point — departure side:**
(xSE, ySE, vxSE, vySE, hSE_M) = (1.00507, 0, 0.0113680, 0.0175449, 4.96074)

**Patch point — arrival side:**
(xEM, yEM, vxEM, vyEM, hEM_M) = (1.00507, 0, 0.0121205, 0.0184927, 4.96074)

**Midcourse maneuver:** ΔVP = 0.036 km/s   [at x = 1.00507]

**Arrival trajectory endpoint** (S–BEM frame):
(xA, yA, vxA, vyA, hMA) = (0.997900, 0.00140873, 0.0249961, 0.0372041, 2.55)
[E^EM_t0 = 851.511, ΔVM = 0.642 km/s]

### Transfer Comparison (Table 1)

| Transfer | ΔVE (km/s) | ΔVM (km/s) | ΔVP (km/s) | ΔVT (km/s) | T (days) |
|---------|-----------|-----------|-----------|-----------|---------|
| Hohmann | 3.141 | 0.838 | — | 3.979 | 5 |
| WSB (Belbruno–Miller 1993) | 3.161 | 0.648 | 0.029 | 3.838 | 160 |
| Proposed 4-body | 3.202 | 0.642 | 0.036 | **3.880** | **100** |

Savings vs Hohmann: 0.099 km/s; extra vs WSB: 0.043 km/s; time saving vs WSB: 60 days.

---

## 4. Reuse Assessment

### #292 BCR4BP
**USE DIRECTLY.** This paper is the primary reference implementation for the BCR4BP departure–arrival
patching strategy. The FTLE grid specification (1000×1000, T=7 nondim) is directly reusable as a
parameter baseline. All mass/distance constants provide the calibration layer for any BCR4BP
implementation. The departure/arrival energy bounds and LCS ICs are the first numerical targets for
any #292 validation.

Key reuse items:
- Coupling strategy: hSE_M = hEM_M at the Poincaré section U (same x-coordinate, y=0, vy>0)
- FTLE truncation-error adaptive grid (Appendix A) for reliable ridge detection
- Table 1 as a three-way comparison baseline (Hohmann, WSB, 4-body)

### #293 ER3BP
Indirect. The BCR4BP decomposition is philosophically similar to splitting the ER3BP into
a periodic perturbation; the LCS detection methodology ports across.

### #267 resonance_network.py
Not directly relevant (paper does not discuss MMR widths or resonant orbits).

### Transport / lobe dynamics
Indirect. The tube manifolds used here are the objects that lobe dynamics quantifies — but this
paper computes them via FTLE/LCS rather than algebraically. The FTLE approach complements the
algebraic lobe approach (Naik 2017) for the non-autonomous 4-body problem.

---

## 5. Mined vs Digested

**DIGESTED.** No catalogue records affected. CORPUS_INDEX entry needed.
