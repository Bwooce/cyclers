# Digest: Rawat, Kumar, Rosengren, Ross (2026) — Cislunar MMR Widths

**Full citation:** Rawat A., Kumar B., Rosengren A.J., Ross S.D., "Cislunar Mean-Motion Resonances:
Definitions, Widths, and Comparisons," *Journal of Guidance, Control, and Dynamics*, Vol. 49, No. 4,
Apr 2026, pp. 1064–1082. DOI 10.2514/1.G009336

**Status: DIGESTED**

---

## 1. Contribution

Defines and quantifies the 2:1 and 3:1 Earth-Moon MMR zones in the PCR3BP using a Poincaré map
at perigee. Two distinct zone types are provided: (a) *stable* resonance width (outermost librational
quasi-periodic torus) and (b) *unstable/chaotic resonance zone* (region enclosed by manifolds of the
unstable resonant periodic orbit, bounded by BIPs). The widths are expressed in geocentric orbital
elements (a, e) and validated against TESS (stable 2:1), IBEX (stable 3:1), and Spektr-R (chaotic
3:1) using JPL Horizons high-fidelity ephemeris propagation.

Companion paper: Kumar–Rawat–Rosengren–Ross 2025, "Cislunar Resonant Transport and Heteroclinic
Pathways: From 3:1 to 2:1 to L1," *Advances in Space Research*, DOI 10.1016/j.asr.2025.12.005
(Ref [42] herein, directly feeds #267).

---

## 2. Model

- **Problem:** PCR3BP, Earth-Moon system, coplanar
- **Mass ratio:** μ = 1.2150584270571545×10^-2
- **Moon mean orbit:** am = 383397.7725 km, em = 0.055545526
- **Poincaré section:** ΣC = {X ∈ MC | l = 0} (perigee, mean anomaly = 0)
- **Coordinates:** cylindrical (ϖ, a) ∈ S¹ × I on ΣC; equivalent to synodic Delaunay variables (g, G)
- **BIPs:** first intersection of stable/unstable manifolds of the unstable resonant periodic orbit;
  parameterise the partial barrier; define chaotic zone boundary
- **Jacobi C scan:** C = 1.60 to 3.54, ΔC = 0.02
- **Tisserand approximation (inclined):** C(a,e,i) ≈ 1/a + 2√(a(1−e²)) cos i
- **Perigee detection function** (avoids false positives near apogee for large μ):
  h(l) = cos(l) + (1/4) sin(l) − 1  [Eq. A1]; zero with positive slope only at l=0

---

## 3. Sourced Goldens

### Resonance Family Extents

| Resonance | Stable family (no Earth collision) | C range (stable) |
|-----------|-------------------------------------|------------------|
| 3:1       | C = 2.44 to 3.47                    | 2.10–3.47 total (Earth collision C ≤ 2.44) |
| 2:1       | C = 1.94 to ?                       | starts C = 1.60 (Earth collision C ≤ 1.94) |

### Resonance Widths in Semi-Major Axis (averaged across all Jacobi C)

| Resonance | Stable width (km) | Unstable/chaotic zone width (km) |
|-----------|-------------------|----------------------------------|
| 2:1       | **28,328**        | **83,982**                       |
| 3:1       | **20,281**        | **34,619**                       |

Source: Section VII, p. 1075, lines 616–620: "The stable resonance zones of 2:1 and 3:1 span
approximately 28,328 and 20,281 km in the semi-major axis, respectively, when averaged across all
Jacobi constants. The regions of influence of 2:1 and 3:1 unstable resonances respectively span
approximately 83,982 and 34,619 km."

### Spacecraft Jacobi Constants (from 2-year high-fidelity propagation)

| Spacecraft | Cmean   | imean (deg) | MMR      | Regime   |
|-----------|---------|-------------|----------|----------|
| TESS      | 2.69889 | 28.8        | 2:1      | stable   |
| IBEX      | 3.10281 | 12.83       | 3:1      | stable   |
| Spektr-R  | 2.57688 | 54.32       | 3:1      | chaotic  |

Source: Section VIII / Fig. 15 caption, p. 1079, line 951-952.

### Poincaré–Cartan Integral Invariants
The oriented sum of areas enclosed by librating islands over canonical-conjugate plane pairs
{(g,G), (h,H), (l,L)} is invariant through sequential perigee mappings; l=0 section forces the
(l,L) projection to zero so the full area appears in (g,G)+(h,H). Numerical confirmation for
TESS and IBEX shown (Appendix A.3).

---

## 4. Reuse Assessment

### #267 resonance_network.py / Kumar 2025 heteroclinic scorer
**USE DIRECTLY.** The two stable width values (28,328 km / 20,281 km) and chaotic zone
widths (83,982 km / 34,619 km) are the resonance scoring inputs that resonance_network.py
and the Kumar 2025 heteroclinic connection paper require. Specifically:
- Stable widths → libration island half-widths for torus enclosure tests
- Chaotic zone widths → extent of manifold-transport corridors connecting 3:1 to 2:1 to L1
- BIP construction methodology → template for defining partial barriers in a wider width atlas
- Spacecraft validation (TESS/IBEX/Spektr-R) → positive control for verifying the width code

**To-do for #267:** ingest 28,328 km (2:1 stable) and 20,281 km (3:1 stable) as width constants;
add chaotic zone extents 83,982 / 34,619 km as corridor bounds.

### #292 BCR4BP
Indirect. The PCR3BP widths are the zeroth-order baseline; BCR4BP adds a Sun perturbation that
will shift and modulate them. Useful as the unperturbed reference.

### #293 ER3BP
Indirect, same as #292.

### Transport / lobe dynamics
**USE DIRECTLY.** The BIP construction is exactly the partial-barrier parameterization that Naik
2017 (Lober) needs as input. The BIP at the first manifold intersection is the turnstile-defining
point q in Naik Eq. (5.5). Rawat 2026 provides the manifold geometry; Naik 2017 provides the
lobe area computation — they compose.

---

## 5. Mined vs Digested

**DIGESTED.** Numerical goldens above are from the paper. No catalogue entries implicated (resonance
widths are inputs to scoring, not cycler records). No search to re-run. CORPUS_INDEX entry needed.
