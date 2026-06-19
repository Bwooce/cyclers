# Digest: Wilczak & Zgliczyński (2003, 2005) - Heteroclinic Connections in the PCR3BP

## Source Documents
1. **Part I:** Wilczak & Zgliczyński (2003). *Heteroclinic Connections between Periodic Orbits in Planar Restricted Circular Three Body Problem - Part I*. Comm. Math. Phys. (arXiv:math/0201278).
2. **Part II:** Wilczak & Zgliczyński (2005). *Heteroclinic Connections between Periodic Orbits in Planar Restricted Circular Three Body Problem - Part II*. (arXiv:math/0401146, DOI 10.1007/s00220-005-1374-x).

## System Parameters
- **Mass Ratio ($\mu$):** 0.0009537 (Sun-Jupiter system)
- **Jacobi Constant ($C$):** 3.03 (Oterma comet energy level)
- **Poincaré Section ($\Theta$):** $\{y = 0\}$, parameterized by $(x, \dot{x})$.

## Key Results & Contributions
This two-part work provides a rigorous computer-assisted topological proof of the existence of homoclinic and heteroclinic connections between the Lyapunov orbits $L_1^*$ and $L_2^*$ in the PCR3BP. These trajectories represent resonance transitions that the comet Oterma makes between the interior (Sun) and exterior regions of the Jovian system.

### Part I: The Heteroclinic Cycle and 4-Symbol Dynamics
The authors established the foundational topological method based on covering relations of "h-sets" (parallelograms with defined stable and unstable edges). 
By computing rigorous bounds on Poincaré maps along numerically approximated trajectories, they proved:
- **Lyapunov Orbits:** The exact existence of $L_1^*$ at $x \approx 0.9208$ and $L_2^*$ at $x \approx 1.0819$, and their strict hyperbolicity.
- **Heteroclinic Connections:** The existence of transversal heteroclinic connections from $L_1^*$ to $L_2^*$ (and vice-versa by symmetry).
- **Homoclinic Connections (4-Symbol):** The existence of a homoclinic orbit to $L_1^*$ in the interior region (corresponding to a 3:2 resonance) and a homoclinic orbit to $L_2^*$ in the exterior region (corresponding to a 1:2 resonance). 

These form a 4-symbol dynamical system where the comet can transition between $L_1^*$, $L_2^*$, the interior region $S$, and the exterior region $X$.

### Part II: The 6-Symbol Dynamics and Higher-Order Resonances
Part II expands the symbolic dynamics by finding new, geometrically distinct homoclinic orbits:
- **5:3 Interior Homoclinic:** A new homoclinic connection to $L_1^*$ representing a 5:3 resonance with Jupiter.
- **2:3 Exterior Homoclinic:** A new homoclinic connection to $L_2^*$ representing a 2:3 resonance.

With these additional connections, the topological framework was extended to a 6-symbol dynamical system $\{L_1, L_2, X, E, I, S\}$. This mathematically proves that a comet (like Oterma) can undergo an infinite number of random transitions between the 1:1, 3:2, 5:3, 1:2, and 2:3 resonances. The authors also introduced a method to detect an infinite number of symmetric periodic and homoclinic orbits by exploiting the time-reversal symmetry of the PCR3BP.

## Extracted Golden Dataset
The system constants, fixed points, stable/unstable eigenvectors, and the $(x, \dot{x})$ crossing coordinates for the heteroclinic and homoclinic orbits (both 4-symbol and 6-symbol) have been transcribed directly into `data/golden/wz_oterma_heteroclinic.yaml` for validation testing of the `#314` heteroclinic-cycle corrector.
