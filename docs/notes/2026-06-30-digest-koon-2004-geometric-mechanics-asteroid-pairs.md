# Koon, Marsden, Ross, Lo & Scheeres 2004 — Geometric Mechanics and Dynamics of Asteroid Pairs (NYAS 1017) — digest

Digested 2026-06-30. Focus: full two-body problem (F2BP) for binary asteroids —
geometric mechanics reduction, restricted F2BP phase space structure (realms, tubes,
lobes), and transport rates. Relevance check for **#494 (k1,k2 binary cyclers /
circumbinary)**.

**Source (cite exactly):**
Koon, W. S., Marsden, J. E., Ross, S. D., Lo, M., & Scheeres, D. J., "Geometric
Mechanics and the Dynamics of Asteroid Pairs," *Annals of the New York Academy of
Sciences*, Vol. 1017, 2004, pp. 11–38. DOI: 10.1196/annals.1311.002. Caltech,
JPL, Univ. Michigan. © 2004 New York Academy of Sciences.

> 28 pp. (incl. extensive references), clean digital typeset. Proceedings paper from
> a 2003/2004 NYAS meeting. Self-described as a "preliminary results" paper and
> programmatic overview; the main contribution is the systematic geometric mechanics
> framework and the RF2BP (restricted full two-body problem) results, not a definitive
> data paper. Many statements are prospective ("we intend to...").

---

## 0. HEADLINE VERDICT

**This paper is about the dynamics of the binary pair itself (the "full two-body
problem"), NOT about a spacecraft orbit around a binary primary.** The
Restricted F2BP (RF2BP) is the limit in which one binary member is massless —
analogous to the CR3BP — and its phase space (exterior realm, interior realm,
saddle-point tubes) is qualitatively similar to the CR3BP. However, the paper does
not develop circumbinary spacecraft trajectory families, does not parameterize
(k1,k2) resonance indices, and does not address the free-flyer cycler problem.

**For #494**: indirect background only. The SE(3) reduction and energy-momentum
method establish the theoretical language; the RF2BP Jacobi integral and 4-saddle
topology provide a qualitative analogue to the CR3BP. No direct formulas for
circumbinary (k1,k2) cycler design emerge.

---

## 1. METHOD / MODEL

### 1.1 Full Two-Body Problem (F2BP)

Two extended rigid bodies interacting gravitationally. Configuration space:
Q = SE(3) × SE(3). Symmetry group: SE(3) (diagonal left action — overall translations
and rotations). Reduction gives shape space Q/G ≅ SE(3), coordinatized by relative
attitude A = A₂ᵀA₁ and relative position R = A₂ᵀ(r₁ - r₂).

Equations of motion (simplified, one body spherical; Eq. 1):

    r″ + 2ω × r′ + ω′ × r + ω × (ω × r) = ∂U/∂r
    I·ω′ + ω × I·ω = -μ · r × ∂U/∂r

where ω = rotational velocity (body-fixed frame), r = relative position, I = specific
inertia tensor, U = gravitational potential. Conserved quantities: angular momentum K,
energy E (Eq. 2).

Free parameter: μ = M₁/(M₁+M₂) (same notation as CR3BP mass ratio).

### 1.2 Restricted Full Two-Body Problem (RF2BP)

Limit μ → 0 (massless test body in rotating distended body's gravity field).
Equations (normalized, Eq. 3):

    ẍ - 2ẏ = -∂U/∂x
    ÿ + 2ẋ = -∂U/∂y

where:
    U(x,y) = -1/√(x²+y²) - ½(x²+y²) + U₂₂
    U₂₂ = -3C₂₂(x² - y²) / (x²+y²)^(5/2)

**Jacobi integral** exists (Eq. 4): E = ½(ẋ² + ẏ²) + U(x,y) — exactly analogous
to the CR3BP Jacobi constant.

**Equilibrium points**: 4 of them, symmetrically placed along the x and y axes, each
at radius R ~ 1 + O(C₂₂). Those on the x-axis (along long axis) are saddle points;
those on y-axis are stable (center × center). Qualitatively the same saddle structure
as L1/L2 in the CR3BP.

### 1.3 Phase Space: Realms and Transport

**Energy threshold E_S**: energy of the symmetric saddle points (x-axis equilibria).
- For E > E_S: bottleneck appears, enabling movement between interior realm (R < 1,
  near asteroid) and exterior realm (R > 1, away from asteroid).
- For E ≤ E_S: no such movement possible.

**Tube dynamics** (between realms): 2D stable/unstable manifold tubes of the periodic
orbits around the saddle points partition the energy surface. Particles inside a tube
transit between realms; particles outside do not. The tube structure (Fig. 3A) is
directly analogous to the L1/L2 tubes in the CR3BP.

**Lobe dynamics** (within a realm): transport between regions within a realm is
mediated by turnstile lobes enclosed by segments of stable/unstable manifolds of
hyperbolic fixed points. Flux per iterate and cumulative transport computed from lobe
areas (Wiggins 1992 formalism).

**Mixed phase space**: stable KAM tori + chaotic sea (Fig. 5). "Lobes of ejection"
(finger-like structures in exterior realm phase space, Fig. 5A) determine escape rate.

---

## 2. SOURCED GOLDENS (verbatim with page references)

| Quantity | Value | Source |
|---|---|---|
| Fraction of near-Earth asteroids that are binaries | ~20% (Margot et al. 2002) | p. 14 |
| Ellipticity parameter C₂₂ range | 0 to 0.05 for physical systems | p. 23, §RF2BP |
| Example C₂₂ (figures) | C₂₂ = 0.05 | p. 24, Fig. 3A caption |
| Equilibrium point radius | R ~ 1 + O(C₂₂) | p. 23–24 |
| Energy threshold for realm crossing | E > E_S (saddle-point energy along x-axis) | p. 24 |
| Example system energy | "slightly negative" (total energy of sphere + tri-axial ellipsoid of equal mass) → Hill stable (cannot mutually escape) | p. 15, Fig. 2 caption |
| 4 equilibrium points | symmetrically on x-axis (2 saddles) and y-axis (2 centers) | p. 23–24 |

**Key qualitative result with quantitative content**: Fig. 5B shows that the argument
of periapse with respect to the rotating asteroid has fewer surviving particles in the
fourth quadrant (270°–360°), citing agreement with Scheeres et al. 1996 numerical
results for asteroid 4769 Castalia. This is a geometrically-derived selectivity result
for ejecta direction.

**No new trajectory numbers**: this is a framework paper. The only numbers are
dimensional parameters for illustration.

---

## 3. REUSE VERDICTS

### #494 (k1,k2 binary cyclers / circumbinary)

**STATUS: BACKGROUND THEORY ONLY — NOT A DIRECT IMPLEMENTATION SOURCE.**

What the paper provides that is potentially relevant:
- **SE(3) symmetry reduction**: systematic language for describing dynamics in a
  body-fixed frame of the binary. If the binary primary is the "extended body" in
  an F2BP, the reduced equations give the correct frame transformation.
- **RF2BP Jacobi integral**: if a spacecraft is treated as the massless particle (μ→0)
  in the gravity field of the rotating distended binary, the RF2BP applies. The Jacobi
  integral E = ½(ẋ²+ẏ²) + U(x,y) would constrain allowable regions of motion,
  including exterior realm access (circumbinary regime).
- **4-saddle topology**: qualitative analogue of L1/L2 for a binary extended body.
  Saddle-point periodic orbits and their manifolds would bound circumbinary orbits
  the same way L1/L2 bound Jupiter-family orbits.
- **Energy-momentum method**: for determining stability of relative equilibria —
  relevant to which circumbinary resonance families are stable vs unstable.

What the paper does NOT provide:
- No parameterization of (k1,k2) resonance indices for spacecraft circumbinary orbits.
- No family of periodic/quasi-periodic circumbinary trajectories with computed periods
  or stability eigenvalues.
- No mass ratio (μ) grid search results.
- The RF2BP treatments are primarily about ejecta and collision dynamics (Fig. 5),
  not about exterior-realm spacecraft tour design.

**Assessment**: The Ross-RT 2026 (k1,k2) work is almost certainly using a different
but related paradigm — likely the P-type (circumbinary) restricted 3-body problem with
the binary pair as the primary doublet, and (k1,k2) as the resonance ratio between
the spacecraft period and the binary orbital period. The F2BP paper does not address
this directly. However, the Jacobi integral structure and the 4-saddle topology in the
RF2BP are the correct background theory to understand WHY certain circumbinary resonant
families exist and what their stability-boundary energies are.

**Bottom line**: cite this paper as theoretical background for #494 (the RF2BP Jacobi
integral and saddle-point topology as the binary-primary analogue of the CR3BP L1/L2
structure), but it is not the source of the (k1,k2) parameterization or cycler design
formulas.

### #500 (Controlled Keplerian map as moon-tour cycler genome)

**STATUS: NOT DIRECTLY RELEVANT.** The paper's tube dynamics (between realms of the
binary's gravity field) are conceptually related to the L1/L2 tubes used in the MMO
framework (Grover-Ross 2009), but the F2BP is a *different* physical model (extended
body gravity vs. point-mass moon). The lobe/tube methods described are the same
computational tools used in the Keplerian map's supporting theory.

### #267 (resonance_network / resonant transport)

**STATUS: TANGENTIALLY RELEVANT — lobe/tube dynamics framework only.** The paper
applies the same Koon-Lo-Marsden-Ross 2000/2001 tube dynamics to the binary asteroid
F2BP setting. The generic machinery (tube partitions, lobe transport, Poincaré section
analysis in Delaunay variables) is common. No new resonance-network results for the
Sun-Jupiter system. The "pinch points" in Fig. 2C of the F2BP phase space are an
interesting structural feature (resonances between rigid-body rotation and orbital
motion) but this is a different system from the interplanetary resonance network.

---

## 4. NOTES ON KEY RELATIONSHIPS

- **Scheeres 2002 (Celest. Mech. Dyn. Astron. 83:155–169)**: the Hill and Lagrange
  stability conditions for the F2BP — these give the boundary between Hill-stable
  (bound, no mutual escape) and Hill-unstable (possible escape) binary states. For #494,
  these stability thresholds are the energy bounds for circumbinary orbit families.
- **Koon et al. 2000 (Chaos 10:427)** and **2001 (CMDA 81:27–38)**: the same tube
  dynamics framework applied to comets (digested separately); those papers provide
  the foundational resonance-transport theory.
- **Jaffé, Ross et al. 2002 (Phys. Rev. Lett. 89:011101)**: "Statistical theory of
  asteroid escape rates" (ref [35]) — uses the same methods for actual escape-rate
  calculations; provides quantitative transport numbers not in this paper.
- **Scheeres 2004 NYAS 1017**: companion paper in the same volume on stability of
  relative equilibria in the full two-body problem (ref [3]).

---

## 5. STATUS

**DIGESTED.** Not catalogue-eligible (binary asteroid dynamics / framework paper;
no repeating spacecraft cycler trajectory). Background theory for #494; not a direct
implementation source.
