# Grover & Ross 2009 — Controlled Keplerian Map (JGCD 32-2) — digest

Digested 2026-06-30. Focus: controlled Keplerian map for moon-tour gravity-assist
trajectory design in a planet-moon environment; primary input for **#500 (SPECULATIVE,
KEYSTONE)** — the map as a moon-tour cycler genome.

**Source (cite exactly):**
Grover, P., & Ross, S. D., "Designing Trajectories in a Planet–Moon Environment Using
the Controlled Keplerian Map," *Journal of Guidance, Control, and Dynamics*, Vol. 32,
No. 2, March–April 2009, pp. 436–443. DOI: 10.2514/1.38320. Virginia Polytechnic
Institute and State University. © 2008 AIAA.

> 8 pp., clean digital typeset. Builds on Ross & Scheeres 2007 SIAM J. Appl. Dyn.
> Sys. 6(3) 576–596 (ref [9]) as the foundational map derivation. Uses Patched
> Three-Body Approximation (P3BA) introduced in Ross-Koon-Lo-Marsden AAS 03-143.
> Numerical example: Jupiter–Europa–Ganymede four-body system.

---

## 0. HEADLINE VERDICT

**The paper IS the engineering implementation specification for #500.** It provides:
(a) the explicit map equations with control input, (b) the kick-function integral,
(c) a coarse/fine two-phase control algorithm, (d) a patched three-body switching
procedure, and (e) a worked Ganymede→Europa tour with sourced ΔV and time numbers.
No additional derivation is required to implement the #500 genome — this paper
plus Ross-Scheeres 2007 (foundational theory) are sufficient.

---

## 1. METHOD / MODEL

### 1.1 The Keplerian Map (uncontrolled)

A 2D symplectic twist map on the cylinder A = S¹ × ℝ, updating periapse angle ω
and Keplerian energy K at each periapse passage (PCR3BP):

    ω_{n+1} = ω_n + 2{-2K_n}^{-3/2}   [mean motion × period]
    K_{n+1} = K_n + μ · f(ω_n; C_J, K̄)

where:
- K = Keplerian energy = -1/(2a), so a = -1/(2K)
- ω = angle of periapse in rotating frame
- μ = mass ratio of perturbing moon (small parameter)
- C_J = Jacobi constant (≈ twice CR3BP Hamiltonian with sign flip)
- K̄ = average value of K over the map's valid domain (approximation parameter)
- f(ω) = kick function (odd in ω), computed by quadrature (see §1.3)

The map is valid when the spacecraft periapse distance exceeds the perturber's
circular orbit radius. An **apoapse map** exists for the interior realm (spacecraft
apoapse < perturber's orbit). Both are called "Keplerian maps" in context.

Phase space (a vs ω) has a "swiss cheese" structure: connected chaotic sea with
embedded stable resonant islands (Fig. 1b). The biggest energy kicks occur for a
very narrow range of ω near π (periapse slightly behind perturber → positive kick;
periapse slightly ahead → negative kick). These regions are labelled A+ and A−.

### 1.2 The Controlled Keplerian Map

Controlled version adds scalar input u_n ∈ U = [-u_max, u_max]:

    F: (ω_n, K_n; u_n) → (ω_{n+1}, K_{n+1})

    ω_{n+1} = ω_n + 2{-2K_n}^{-3/2}
    K_{n+1} = K_n + μ · f(ω_n + u_n; C_J, K̄)

u_n is applied as an instantaneous ΔV to the periapse angle (approximated as
equivalent to a velocity increment at the kick moment). The β parameter
(= ε·C_J·K̄) governs the size of regions A+ and A−; it is chosen as a design
parameter controlling aggressiveness of the control.

### 1.3 Kick Function (Appendix)

The kick function is:

    f(ω) = -(1/√p) ∫_{-π}^{π} (r³/r₂²) sin(ω - θ) · t dθ
           - ∫_{0}^{π} sin(ω + 2π - θ) · t dθ

where:
    r = p/(1 + e·cos θ)   [orbit radius; p = a(1-e²) = -1/(2K)·(1-e²)]
    r₂ = √(1 + r² - 2r·cos(ω - t))  [distance to perturber at time t]

The Jacobi constant invariance yields f = f(ω; C_J, K̄) — once C_J and K̄ are
fixed, f is a function of ω alone and can be precomputed and stored.

### 1.4 Patched Three-Body Approximation (P3BA)

Multi-moon transfers use two separate Keplerian maps:
- **Apoapse map** in the Jupiter–M1 (outer moon) system: valid when a > a_M1
- **Periapse map** in the Jupiter–M2 (inner moon) system: valid when a < a_M2

**Switching region**: defined in *phase space* (not configuration space) at the
intersection of the two systems' constant-C_J energy contours in the a–e plane.
A switching orbit transitions the dominant perturber from M1 to M2; no rocket
burn is required for this switch.

**Search procedure** (5 steps in §III):
1. Choose starting point (a₀, e₀) on the JGS contour, outside the switching region.
2. Narrow ω_g (Ganymede periapse angle) to the A+ region near ω_g ≈ π.
3. Use predictor-corrector sensitivity to find ω_e (Europa periapse angle) giving
   a further Europa kick at the next periapse.
4. Iterate forward (periapse map in JES) and backward (apoapse map in JGS) to
   find switching-region entry/exit pairs P1f and P1b.
5. Verify in the full four-body equations and cycle through nearby (a,e) values.

### 1.5 Control Algorithm

**Coarse control** (large a-decrements, Ganymede→switching region):
- At each map iterate, look ahead n_max steps with u=0.
- If a future iterate falls in A−, apply a control sequence (coarse discretization)
  to steer it toward ω_opt (maximum Δa decrease per iteration).
- β parameter controls the size of A+ and A− windows.

**Fine control** (switching region→Europa capture):
- Forward-backward targeting: aim for the interior of the first cut of Europa's
  stable invariant manifold of the L2 periodic orbit on the Poincaré section.

---

## 2. SOURCED GOLDENS (verbatim with page references)

All numbers from the JGCD paper unless noted.

| Quantity | Value | Source |
|---|---|---|
| Mass ratio (Jupiter–Callisto example) | μ = 5.667×10⁻⁵ | p. 437, Fig. 1b caption |
| Jacobi constant (example) | C_J = 2.995 | p. 437, Fig. 1b caption |
| Mean semimajor axis (example) | ā = 1.35 (Jupiter–Callisto units) | p. 437, Fig. 1b caption |
| Control discretization step | 1 m/s (ΔV equivalent) | p. 441 |
| Maximum control input | u_max = ±5 m/s (ΔV equivalent) | p. 441 |
| Lookahead horizon | n_max = 5 iterates | p. 441 |
| Ganymede periapse angle search range | 0.90π < ω_g < 0.99π | p. 439 |
| Sample trajectory ΔV (Ganymede→Europa) | **160 m/s** | p. 441 |
| Sample trajectory transfer time | **1.7 years** | p. 441 |
| Number of Jupiter orbits (periapse/apoapse passages) | **116 revolutions** | p. 441 |
| Time-of-flight ratio (controlled vs zero-fuel) | **~10%** (controlled takes ~10% of uncontrolled time) | p. 441, also Abstract |
| Uncontrolled trajectory time (same system) | implied ≈ 17 years (from 10% ratio) | computed from above two |

**Note on 10% ratio**: "The time taken for this mission is less than 10% of that
taken for the optimal (zero) fuel trajectory for the same four-body system" (p. 441).
The uncontrolled zero-fuel trajectory is extremely long (years) because the spacecraft
gets stuck in resonances.

**Periapse angle convention**: kicks are significant for ω near π (behind perturber
is positive; ahead of perturber is negative in the periapse-map).

---

## 3. REUSE VERDICTS

### #500 (SPECULATIVE, KEYSTONE) — Controlled Keplerian map as moon-tour cycler genome

**STATUS: PRIMARY IMPLEMENTATION SOURCE. MINED.**

This paper delivers everything needed to build the #500 genome:

1. **Map equations** (Eq. 1 and 2): explicit, implementable; no additional derivation
   needed.
2. **Kick function** (Appendix): integral formula with all variables defined; can be
   precomputed as a lookup table in ω for fixed (C_J, K̄).
3. **Periapse vs apoapse map selection**: periapse map for exterior realm (a > a_moon),
   apoapse map for interior realm (a < a_moon). Both derived by the same method.
4. **Control parameterization**: scalar u_n per step; coarse (energy-pumping) + fine
   (manifold-targeting) two-phase strategy; β controls aggressiveness.
5. **Switching region / patching**: phase-space-defined switching between moon systems;
   the 5-step procedure in §III is the genome's inter-moon segment planner.
6. **Worked example**: Ganymede→Europa in Jupiter system, 160 m/s, 1.7 yr, 116 orbits.

For the cycler application, the map would be run as a *cyclic* iterative genome:
a sequence of (moon, periapse angle, control input) specifying each gravity-assist
pass. The resonance hopping in the chaotic sea is the natural tour structure;
the control input steers between resonance islands.

**Key upstream reference to also acquire/check**: Ross & Scheeres 2007 SIAM J. Appl.
Dyn. Sys. 6(3):576–596 ("Multiple Gravity Assists, Capture, and Escape in the
Restricted Three-Body Problem") — this is where the kick function f(ω) is *derived*
(cited here as [9]). It is the foundational theory behind the map; the 2009 paper is
the engineering application.

**Limitation noted in the paper**: maps are 2D (planar CR3BP). Extension to 3D yields
4D maps in (ω, Ω, K, i) — noted as future work (§V). For the moon-tour cycler genome,
3D may be required for non-coplanar moon systems.

### #267 (resonance_network / resonant transport)

**STATUS: RELEVANT — but this paper is the APPLICATION layer; the foundational
resonance tube theory is in Koon-Lo-Marsden-Ross 2000 (Chaos 10:427) and 2001
(CMDA 81:27–38, digested separately).**

The Keplerian map directly encodes the resonance hopping mechanism: stable resonant
islands in the phase space correspond to p:q mean-motion resonances with the moon;
the chaotic lanes between islands are the paths a spacecraft uses to drift through
resonance space. The "lobes" mediating movement between resonances (ref [6]) are
explicitly mentioned. For the resonance_network: the map gives a rapid propagator for
exploring which resonances are reachable from a given initial (ω, K) at a given energy.

### #494 (k1,k2 binary cyclers / circumbinary)

**STATUS: NOT RELEVANT.** This paper treats a single planet + one (or two) moons;
there is no binary primary. The controlled Keplerian map applies to the planet-moon
system where the spacecraft orbits the planet and is perturbed by moons. It does not
address circumbinary orbits or the binary-primary paradigm needed for (k1,k2) families.

---

## 4. ACQUISITION NOTES

- Ross & Scheeres 2007 SIAM (ref [9]) — foundational map derivation — **acquire/check
  if not yet in corpus**: "Multiple Gravity Assists, Capture, and Escape in the
  Restricted Three-Body Problem."
- Koon et al. 2002 Contemporary Mathematics (ref [14]) — "Constructing a Low Energy
  Transfer Between Jovian Moons" — closely related implementation; may have additional
  numbers.
- The paper's Figs 5/8 contain trajectory histories (semimajor axis vs time) but these
  are continuous curves, not tables; the only sourced scalar numbers are those in §2
  above.

---

## 5. STATUS

**DIGESTED.** Not catalogue-eligible (methodology paper; no cycler trajectory to ingest).
Maps to #500 implementation directly.
