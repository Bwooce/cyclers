# Digest: Fitzgerald, Ross (2022) — Transit Orbits in Periodically-Perturbed RTBP

**Full citation:** Fitzgerald J., Ross S.D., "Geometry of transit orbits in the periodically-perturbed
restricted three-body problem," *Advances in Space Research*, Vol. 70, 2022, pp. 144–156.
DOI 10.1016/j.asr.2022.04.029

**Status: DIGESTED**

---

## 1. Contribution

Extends the Conley–McGehee equilibrium-region transit-orbit theory (originally for autonomous CR3BP)
to **time-periodic** systems — specifically the Bicircular Problem (BCP) and the Elliptic Restricted
Three-Body Problem (ER3BP). The Lagrange equilibrium points L1/L2 no longer exist in these systems,
so the paper defines and computes their **Lagrange periodic orbit** dynamical replacements via
differential correction. A stroboscopic (period-T) Poincaré map replaces the continuous flow;
the monodromy matrix K gives the transit/non-transit classification in a symplectic eigenbasis.
Result: an effective Hamiltonian H̃2 separating phase space into transit and non-transit orbits,
directly generalising the autonomous Lyapunov orbit saddle geometry.

---

## 2. Model

### Bicircular Problem (BCP)
Sun–Earth–Moon–S/C where the Moon perturbs the Sun–Earth CR3BP with circular Moon orbit.
Nondimensionalised by Earth–Moon distance (l) and angular velocity (ωm0).

**BCP parameters (nondimensional):**
- l = 0.01215 (Sun mass parameter in BCP context, from text)
- l0 = 328900.54 (Earth–Moon distance in km, characteristic length)
- a0 = 388.81114 (Sun–Moon distance nondim)
- ωm0 = 0.925195985520347 (Moon angular velocity nondim)

### Elliptic Restricted Three-Body Problem (ER3BP)
Earth–Moon–S/C with true eccentricity of Moon's orbit.

**ER3BP parameter:**
- Earth–Moon eccentricity: e = 0.0549006

### Stroboscopic map and monodromy
Period-T Poincaré map; monodromy matrix K has eigenvalues {r, 1/r, e^{iw}, e^{-iw}}.
In the symplectic eigenbasis (q1, p1, q2, p2):

H̃2 = k̃ q1 p1 + (m̃/2)(q2² + p2²)

where k̃ = (1/T) ln(r) and m̃ = (1/T) w.  [Eq. B.4: r = e^{k̃T}, w = m̃T]

Transit ↔ {q1 p1 < 0, H̃2 > 0} through the equilibrium region;
non-transit ↔ {q1 p1 > 0 or H̃2 ≤ 0}.

---

## 3. Sourced Goldens

### BCP L1 Lagrange Periodic Orbit (phase h=0)

Initial condition [x, y, px, py] in rotating-frame canonical coordinates:

```
[0.837595408485656, 0, 0, 0.827678389393936]
```

Source: Appendix A (p. 154, lines 729–733). Note: authors thank Jorba and Rosales (Barcelona)
for providing this IC.

### ER3BP L1 Lagrange Periodic Orbit (phase h=0)

Initial condition [x, y, px, py]:

```
[0.792718947200736, 0, 0.000001145970495, 0.886145419995798]
```

Source: Section 6 / Table in paper. Obtained via continuation from BCP L1 IC as e increases
from 0 to 0.0549006 (Appendix C, Fig. C.1).

### BCP L1 Monodromy Eigenvalues

| Quantity | Value |
|---------|-------|
| r (saddle eigenvalue) | 4.2874×10^8 |
| w (center eigenvalue, rad/period) | 3.0273 |
| k̃ = (1/T)ln(r) | (1/T)·ln(4.2874×10^8) |
| m̃ = (1/T)·w | (1/T)·3.0273 |

### ER3BP L1 Monodromy Eigenvalues

| Quantity | Value |
|---------|-------|
| r | 8.3659×10^7 |
| w | 1.9863 |

### Transit Classification Thresholds

| System | H̃2 threshold | c threshold |
|--------|--------------|-------------|
| BCP    | 10^-6        | 10^-4       |
| ER3BP  | 10^-8        | 4×10^-5     |

Transit orbit: |H̃2| < threshold AND q1 p1 < 0 AND |q1 p1| > c.
Non-transit: |H̃2| < threshold AND q1 p1 > 0.

### Effective Hamiltonian System Matrix (Appendix B)

A = diag-block: [[k̃, 0, 0, 0], [0, −k̃, 0, 0], [0, 0, 0, m̃], [0, 0, −m̃, 0]]

Solution: x(t) = e^{At} x(0) with diagonal/rotation blocks:
- q1(t) = e^{k̃t} q1(0), p1(t) = e^{-k̃t} p1(0)
- (q2, p2): rotation at frequency m̃

---

## 4. Reuse Assessment

### #293 ER3BP
**USE DIRECTLY — primary reference.** The paper provides exactly what #293 needs:
- ER3BP L1 Lagrange periodic orbit IC (14 significant figures) for initialising the corrector
- Monodromy eigenvalues (r=8.3659×10^7, w=1.9863) for building the transit gate H̃2
- Continuation methodology from BCP → ER3BP (Appendix C) as the construction path
- Transit/non-transit classification threshold (H̃2=10^-8, c=4×10^-5) as a validation golden

**To-do for #293:** implement the BCP→ER3BP continuation starting from the BCP L1 IC above,
parameterised by eccentricity from e=0 to e=0.0549006; verify against the ER3BP IC above.

### #292 BCR4BP
**USE DIRECTLY — BCP layer.** The BCP L1 periodic orbit IC and monodromy eigenvalues are the
BCR4BP analog. Note that Onozaki 2017 (digested separately) handles the full BCR4BP; Fitzgerald
2022 provides the L1 equilibrium-region geometry that Onozaki uses implicitly.

### #267 resonance_network.py
Indirect. The transit gate H̃2 defines which trajectories can pass through L1 into the resonance
network, so it is an upstream gate for #267 — but no MMR width data is provided.

### Transport / lobe dynamics
**USE DIRECTLY.** The monodromy eigenbasis (q1, p1) saddle plane is exactly where the turnstile
lobes live. H̃2 = const level sets define the Conley–McGehee boxes whose boundaries are the
lobe manifolds that Naik 2017 (Lober) computes areas for. The two tools compose: Fitzgerald 2022
provides the equilibrium-region structure; Naik 2017 provides the lobe area computation.

---

## 5. Cross-references

- Continuation IC sourced from Jorba, Jorba-Cuscó, Rosales 2020 (Celestial Mech. Dyn. Astron.
  132(2)) — if Barcelona group IC ever needs updating, consult that paper directly.
- Onozaki 2017 (digested) is cited here as [Onozaki 2017]; their BCR4BP FTLE approach is the
  trajectory-design complement to this paper's equilibrium-region geometry.
- Naik–Lekien–Ross 2017 (digested) is the lobe-area computation tool; applies to the manifold
  curves whose structure Fitzgerald 2022 characterises.

---

## 6. Mined vs Digested

**DIGESTED.** No catalogue records affected. CORPUS_INDEX entry needed.
