# Digest: Ross–BozorgMagham–Naik–Virgin 2018, PRE 98, 052214

**File:** `ross-bozorgmagham-naik-virgin-2018-phase-space-conduits-transition-PRE-98.pdf`
**Citation:** S. D. Ross, A. E. BozorgMagham, S. Naik, L. N. Virgin, "Experimental validation of
phase space conduits of transition between potential wells," *Phys. Rev. E* **98**, 052214 (2018).
DOI: 10.1103/PhysRevE.98.052214.
**Digested:** 2026-06-30.

---

## 1. Model and Method

**N DOF Hamiltonian system with a rank-1 saddle.** Near the saddle, the Hamiltonian in normal
(eigenvector) coordinates is:

    H₂(q₁, p₁, …, qN, pN) = λ q₁ p₁ + Σ_{k=2}^{N} (ωk/2)(qk² + pk²)     … (1)

where λ is the real eigenvalue of the saddle (reactive coordinate), ωk are the centre frequencies
(bath coordinates). N=2 is the main case.

**NHIM (normally hyperbolic invariant manifold):** In the equilibrium region, the (2N−3)-sphere
M^{2N-3}_h = {q₁=p₁=0, Σ_k (ωk/2)(qk²+pk²) = h} is the NHIM; for N=2 it is an unstable
periodic orbit (Lyapunov orbit). Acts as the "anchor" for the separatrices.

**Separatrices (manifold tubes):** Stable manifolds W^s±(M^{2N-3}) have topology (2N-2)-tube
(S^{2N-3} × ℝ); for N=2 a cylinder R¹×S¹. Transit trajectories lie INSIDE the tube; nontransit
bounce OUTSIDE.

**Experimental system (Sec. III/IV):** Ball rolling on a machined 4-well surface in 2 DOF. Height
function:

    H(x,y) = α(x²+y²) − β(√(x²+γ) + √(y²+γ)) − ξxy + H₀

Parameters (α, β, γ, ξ, H₀) = (0.07, 1.017, 15.103, 0.00656, 12.065) in cm units.
Surface accurate to 0.003 mm. Camera: Prosilica GC640 at 50 Hz, 0.16 cm pixel resolution.
120 trajectories of ~10 s; ~4000 Poincaré section intersections recorded; 400 transition events.

---

## 2. Sourced Goldens

### Transition area scaling law (Sec. III, p. 052214-4)

For fixed excess energy ΔE > 0 above the saddle critical energy Ee:

    A_trans = T_po * ΔE     (to leading order in ΔE)         … from text between eqs (11)–(13)

where T_po = 2π/ω is the **period of the unstable periodic orbit** (Lyapunov orbit) in the
bottleneck at energy E, and ω is the imaginary part of the complex-conjugate eigenvalue pair of
the linearisation about the saddle.

**Full transition fraction (eq. 13):**

    p_trans = A_trans / A_E
            = (T_po / A₀) * ΔE * [1 − (τ/A₀) * ΔE + O(ΔE²)]

where:

    A₀ = 2 ∫_{r_min}^{r_max} (14/5) [Ee − gH(r)] √(1 + 4H_r²(r)) dr      … eq. (11)

    τ  = ∫_{r_min}^{r_max} (14/5) (1 + 4H_r²(r)) / [Ee − gH(r)] dr         … eq. (12)

(for the specific rolling-ball system; the leading-order formula A_trans = T_po * ΔE is general for
any N=2 Hamiltonian with rank-1 saddle — ref [37] MacKay 1990).

**Experimental growth rate (Sec. V / Fig. 6b):** Linear fit slope for small ΔE:

    (1.0 ± 0.23) × 10⁻³ (s/cm)²

**Predicted growth rate from theory:** T_po/A₀ ≈ 0.87 × 10⁻³ (s/cm)² (stated in text, p. 052214-4).

**Agreement:** Experimental regions of transition account for > 99% of observed transition trajectories.
Agreement to within 1% (abstract).

**Lyapunov time of system (Sec. IV):** ≈ 0.4 s (used as mixing time before recording intersections).

**Damping ratio (Sec. III):** ζ ≈ 0.025 (small; justifies short-time conservative dynamics for
ΔE/E << 1 at E > 1000 cm²/s²).

### The law stated compactly
For any N=2 Hamiltonian system with a rank-1 saddle, the **tube cross-section area** on a
Poincaré section is, to leading order in ΔC = C₁ − C (where C₁ = C_saddle):

    A_tube(ΔC) = T_po(C_saddle) * ΔC     [units: phase-space area]

This is the MacKay 1990 flux formula (ref [37]); Ross et al. experimentally confirm it.

---

## 3. Key Physical Insights

1. **Tube = separatrix.** The (S¹)-curve on a 2D Poincaré section IS the transition boundary:
   all points inside transit, all outside bounce. There is no grey zone (in the conservative limit).

2. **Linear growth in ΔC.** For small excess energy above the bottleneck UPO, the tube's
   cross-sectional area (and hence the fraction of initial conditions that transit) grows **linearly**
   in ΔC. Nonlinear corrections are O(ΔC²).

3. **Robustness to dissipation.** The tube persists and bounds transit to within 1% even for
   ζ ≈ 0.025 damping, over relevant timescales. KAM tori and periodic orbits are fragile to
   damping; the tube is not.

4. **Higher DOF extension.** For N=3 (4D Poincaré section), the NHIM has topology S³; tube
   cross-section is S³; projects to two transverse planes each showing a disk — standard method
   from Gabern et al. 2005 (ref [28]).

---

## 4. Reuse Assessment

### #494 (k₁,k₂) cycler construction — seed-window sizing
**This is the direct reference for the ΔC tube-cross-section growth law cited as ref [28]
in the Ross-RT 2026 paper.** Confirms:

    A_tube ∝ T_po(C_J) * (C_J − C)     for C just below C_J

where C_J is the Jacobi constant at L_J (the libration point energy) and T_po is the UPO period.

Practical implication: when sizing the seed window for the (k₁,k₂) construction, the manifold
tube's Poincaré-section cross-section has LINEAR radius scaling ~ √(T_po * ΔC). A doubling of
ΔC means √2 larger tube radius. The window should scale accordingly to guarantee intersection
coverage. The leading-order formula A_trans = T_po * ΔC gives the AREA; radius in each
transverse direction scales as √(T_po * ΔC / π) for a circular cross-section.

### #314 `heteroclinic_cycle.py` (Lyapunov manifold corrector)
**Supporting reference.** Confirms the tube-crossing detection logic: a Poincaré section cut of
the manifold tube appears as a closed curve; points INSIDE the curve transit through the bottleneck.
The area inside this closed curve = A_trans = T_po * ΔC to leading order. Use this as a sanity
check on the tube cross-section size extracted by #314's manifold integration.

### #411 cross-system Newton stall at |R|=0.59 rad
**Indirect.** Ross-BozorgMagham 2018 confirms that the tube IS a genuine codimension-1 separatrix
in 2 DOF (N=2), so the two manifolds MUST intersect on the Poincaré section if both are at the
same energy. The |R|=0.59 rad stall in #411 is not a symptom of missing intersection (the manifolds
do cross); it is a corrector basin / branch-tracking problem. This paper gives no algorithmic fix for
the Newton stall. See Braik-Ross 2025 for the multi-crossing strategy.

### #496 two-phase corrector port
**Supporting reference.** Confirms that the tube cross-section on a Poincaré section is a CLOSED
CURVE (S¹ topology) for N=2, so "tube vs tube" intersection search on the Poincaré section is
well-posed. The closest-pair initial guess strategy (from Braik-Ross 2025) exploits exactly this
topology: each manifold's tube appears as one closed curve, and their intersection is a point (or
set of points) on the section.

---

## 5. Corpus Status Recommendation
**`digested`** — cited as ref [28] in Ross-RT 2026 OUTSTANDING thread and in the #494 task
description. Paper is a physical-experiment validation of tube dynamics; no catalogue rows are
mined from it. Status: `digested`.
