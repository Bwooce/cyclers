# Digest: Canales, Howell, Fantino, Gilliam 2023

**Paper:** "Transfers between moons with escape and capture patterns via Lyapunov exponent maps" (arXiv:2308.10029)
**Authors:** David Canales, Kathleen C. Howell, Elena Fantino, Annika J. Gilliam
**Date:** Aug 19, 2023

## Key Values for FTLE Testing (Ganymede-Europa Setup)

### System Mass Ratios (CR3BP)
- **Jupiter-Europa (J-E):** $\mu = 2.528 \times 10^{-5}$
- **Jupiter-Ganymede (J-G):** $\mu = 7.804 \times 10^{-5}$

### Departure Map (Jupiter-Ganymede)
- **Poincaré Section:** $x = 0.965$ (in J-G rotating frame, near $L_1$)
- **Jacobi Constant ($J_d$):** $3.00754$
- **Grid Bounds:** $y \in [-0.006, 0.015]$, $\dot{y} \in [-0.01, 0.02]$
- **Horizon ($t_d$):** $-10$ normalized J-G time units (backward integration)
- **Reference Figure:** Fig. 4(b) - Departure FTLE map

### Arrival Map (Jupiter-Europa)
- **Poincaré Section:** $x = 1.028$ (in J-E rotating frame, near $L_2$)
- **Jacobi Constant ($J_a$):** $3.00240$
- **Grid Bounds:** $y \in [-0.018, 0.05]$, $\dot{y} \in [-0.04, 0.01]$
- **Horizon ($t_a$):** $+10$ normalized J-E time units (forward integration)
- **Reference Figure:** Fig. 4(d) - Arrival FTLE map

## Methodology
The paper combines Finite-Time Lyapunov Exponent (FTLE) maps with the Moon-to-Moon Analytical Transfer (MMAT) method. FTLE maps are constructed on Poincaré sections at the CR3BP gateways ($L_1$/$L_2$). The maximum stretching direction determines the FTLE scalar field. Isolines of specific behaviors (e.g., capture, transit) are identified from the FTLE structures. MMAT provides the geometric intersection conditions for conics outside the moons' SoI. Matching arrival constraints (from MMAT) with the departure FTLE map yields feasible single-impulse moon-to-moon transfers with specific endgame characteristics.
