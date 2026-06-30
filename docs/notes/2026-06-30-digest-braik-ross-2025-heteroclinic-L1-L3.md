# Digest: Braik–Ross 2025, AAS 25-716

**File:** `braik-ross-2025-heteroclinic-transfer-L1-L3-earth-moon-AAS-25-716.pdf`
**Citation:** A. Braik, S. D. Ross, "Heteroclinic Transfer Between L1 and L3 in Earth-Moon System,"
AAS 25-716, AAS/AIAA Space Flight Mechanics Meeting (2025).
**Digested:** 2026-06-30.
**Note:** Conference paper (AAS), no arXiv source. Companion preprint is
`braik-ross-2026-orbital-networks-three-body-problem-arxiv-2605.31543.pdf` (also in corpus).

---

## 1. Contribution

First demonstration of **heteroclinic connections between L₁ and L₃ Lyapunov orbits in the
Earth-Moon CR3BP**. Prior work focused exclusively on L₁ ↔ L₂. Two distinct transfer families
(H⁺ and H⁻) are identified and characterised over a range of Jacobi constants.

Key novelty:
- L₃ (opposite side of Earth from Moon) is quasi-stable, has lower unstable eigenvalues than L₁/L₂;
  some L₃ periodic orbits can remain stable for up to 250 days without control.
- L₃ stationkeeping cost cited: ΔV = 0.02 m/s over 60 days (ref [13], Conti & Circi 2025).

---

## 2. Model

Planar CR3BP, Earth-Moon system. Mass ratio μ = m₂/(m₁+m₂) (Moon/total; standard EM value,
not stated explicitly but implied by the critical Jacobi constants). Jacobi constant definition:
C(x,y,ẋ,ẏ) = −2Ū(x,y) − (ẋ²+ẏ²) where Ū = −½(x²+y²) − (1−μ)/r₁ − μ/r₂.
Note: their definition OMITS the μ(1−μ) constant (differs from conventions that put C(L₄/L₅)=3).

---

## 3. Sourced Goldens

### Critical Jacobi Constants — Table 1, p. 3
| Lagrange point | C_i |
|---|---|
| L₁ | **3.1883** |
| L₂ | **3.1722** |
| L₃ | **3.0122** |
| L₄ | **2.9880** |
| L₅ | **2.9880** |

These are sourced (tabulated, not figure-derived). Convention: dropping the μ(1−μ) term.

### Heteroclinic connection existence window (Sec. "Selection of Jacobi Constant", p. 4)
- **Necessary condition:** C < C₃ = 3.0122 (L₃ neck open).
- **Sufficient (manifold intersection) condition:**

      C < C_max ≈ 3.0015                     … eq. (6)

  Above this value, no intersection exists on the Poincaré section. This bound was found by
  numerical exploration; not analytically derived.

- **Practical search range:** C = 2.89 to C = 2.99 (Sec. "Results", p. 13). Upper limit chosen
  to avoid near-tangency at C_max; lower limit because larger orbits demand excessive seeding.

### Lyapunov orbit corrector parameters (Sec. "Lyapunov Orbit Generation", p. 5)
- Initial guess: from JPL Three-Body Periodic Orbits Database (ref [18]).
- Convergence threshold: **ε_L = 10⁻¹²** (|ẋ₁| at half-period crossing).
- Energy constraint enforcement: ẏ₀ = √(−2Ū(x₀,0) − C_target), positive root.
- STM propagated simultaneously: Φ̇ = Df(X(t)) Φ, Φ(0) = I₄ₓ₄.

### Manifold seeding parameters (Sec. "Invariant Manifold Computation", p. 7)
- Perturbation amplitude: **ε_M = 10⁻⁶ nondim ≈ 0.384 km**.
- Phase sampling: N points evenly spaced τ_k = (k−1)/(N−1), k=1…N; N chosen based on orbit
  amplitude.
- Both "+" and "−" branches seeded.

### Integration tolerance (Sec. "Heteroclinic Connection Detection", p. 10)
- Runge-Kutta absolute and relative tolerance: **10⁻¹⁴** for heteroclinic trajectory construction.

### Poincaré section (Sec. "Poincaré Sections", p. 4)
- Section: **Σ_C = {(y, ẏ) | x = −μ, C = constant}**
  (i.e. at the Earth location in nondimensional rotating coordinates).
- Direction: all crossings recorded (multi-crossing strategy — see below).

### Orbital periods (Sec. "Lyapunov Orbit Generation", p. 6 / Fig. 3b)
- **L₃ family period:** essentially constant at **≈ 27 days** across the C range studied.
  (figure-derived from Fig. 3b; the text states "approximately 27 days").
- **L₁ family period:** decreases steadily with increasing C (smaller orbits).

### Transfer time results (Sec. "Transfer Time Analysis", p. 15 / Fig. 11)
- Search range: C ≈ 2.89 to 2.99.
- **H⁻ family consistently faster than H⁺ by ≈ 4–5 days** at any given C (text, p. 16):
  "an H⁺ transfer typically takes 4.5 days longer than the corresponding H⁻ transfer."
- Discrete jumps in TOF curve: each jump ≈ **one L₃ orbit period ≈ 27−28 days**
  (text, p. 15: "the size of the jump in each curve is approximately equal to one L₃ orbit period").
- Jump locations: H⁻ curve jumps near C ≈ 2.98; H⁺ curve near C ≈ 2.97 (text, p. 15).
- Time unit: 1 nondim time = **4.345 days** (sidereal period 2π nondim = 27.3 days).

### Manifold intersection geometry (Sec. "Results and Discussion", p. 13)
- All trajectories (both families, all C) cross x ≈ −μ near:

      (x, y) ≈ (−μ, −0.93)    [nondimensional rotating frame]     (text, p. 13)

  This "narrow band" is independent of Jacobi constant — the connection corridor.

### L₃ quasi-stability (Introduction, p. 2)
- Asteroid capture ΔV ≈ 20 m/s using L₃ stable manifold (ref [17], Jorba & Nicolás 2021).

---

## 4. Key Algorithm: Multi-Crossing Poincaré Map (Sec. "Poincaré Section Recording", p. 7–9)

**Problem with single-crossing:** Standard approach (record only FIRST crossing with prescribed ẋ
sign) misses manifold segments that execute "local loops" near x = −μ in the L₁ ↔ L₃ geometry.
These loops cause multiple rapid section piercings; a first-crossing filter leaves gaps in the (y, ẏ)
map, destroying the ability to detect intersections.

**Solution — two-event strategy:**
1. **Multi-crossing event:** Record EVERY crossing of x = −μ, regardless of ẋ sign.
2. **Region-limited integration:** Terminate each trajectory once it exits a predefined x-bounds
   neighbourhood of x = −μ. Prevents contamination from distant returns.

This "multi-crossing Poincaré map" fills the gaps and reveals the complete manifold structure.

**Closest-pair intersection detection (Sec. "Heteroclinic Connection Detection", p. 10):**
1. Propagate dense W₁ˢ⁻ (L₁ stable, "−" branch) and W₃ᵘ⁺ (L₃ unstable, "+" branch).
2. For each stable-manifold point p_s, compute Euclidean distance in (y, ẏ) to all unstable-
   manifold points p_u WHERE sign(ẋ_s) = sign(ẋ_u) (same direction of section crossing).
3. Select the pair with MINIMUM distance as the heteroclinic point candidate.
4. At each C, EXACTLY TWO closest pairs exist — one at ẏ > 0 and one at ẏ < 0 — giving two
   distinct transfer families H⁺ and H⁻.

**Refinement noted as future work:** fit smooth curves to the Poincaré cuts and use root-finding
(differential correction) to converge to exact intersection. Closest-pair used here with negligible
separation sufficient for trajectory patching.

---

## 5. Transfer Family Geometry

**H⁻ family (ẏ < 0 crossing pair):**
- Spacecraft departs L₃, executes an extra loop on the L₃ / far-Earth side BEFORE crossing x=−μ.
- More direct path once past x=−μ toward L₁.
- Faster: lower average semimajor axis (closer to Earth on average).

**H⁺ family (ẏ > 0 crossing pair):**
- Spacecraft crosses x=−μ more directly from L₃ side, then executes a wide loop on the L₁ / near-
  Moon side before settling into L₁.
- For lower C values, the H⁺ loop swings near the **L₅ vicinity** before turning back.
- ~4.5 days slower than H⁻ at same energy.

**Time-reversal symmetry** (eq. 16: (x,y,ẋ,ẏ,t) → (x,−y,−ẋ,ẏ,−t)) gives H₁→₃±
as the mirrors of H₃→₁±. Four connections total at each C.

---

## 6. Reuse Assessment

### #411 cross-system Newton stall at |R|=0.59 rad
**High-value analog. This paper addresses the closest structural analog to #411's problem.**

Both scenarios involve two libration points in a single PCR3BP system (L₁ and L₃ here; EM and SE
orbits in #411 after being bridged to a common inertial section). The key finding:

1. **The single-crossing Poincaré map fails in the L₁↔L₃ geometry** — exactly the same failure
   mode could be causing the Newton stall in #411 if the initial manifold scan misses the true
   intersection due to loop-induced crossing gaps.

2. **Fix:** multi-crossing strategy (record all crossings in a bounded neighbourhood) and
   closest-pair search for the initial guess. This is cheaper and more robust than a cold Newton.

3. **The (y, ẏ) plane closest-pair approach** gives an initial guess within the same basin as the
   true intersection; THEN Newton (or differential correction) converges from there.

**Direct application to #411 stall:** The current #411 corrector starts Newton from a phase-space
guess at the patch section with |R|=0.59 rad residual. Braik-Ross 2025's lesson is: populate the
Poincaré section DENSELY (multi-crossing), then use the closest pair as the Newton seed — this
gives a much better basin entry. The correction algorithm (2×2 Newton over c_em, c_se) can be
retained; only the SEED needs improvement.

### #496 two-phase corrector port
**Direct blueprint.** Braik-Ross 2025 implements exactly the two-phase corrector:
- Phase 1: multi-crossing Poincaré map → closest-pair initial guess.
- Phase 2: differential correction from that initial guess to refine to machine precision.

The "improved intersection precision" is explicitly called out as future work (Conclusion item 2),
confirming the current paper is Phase 1 only and that Phase 2 (differential correction) is the
natural next step. The #496 port should implement this two-phase structure.

Specific implementation detail: sign(ẋ_s) = sign(ẋ_u) filter on the closest-pair search prevents
false pairings between manifold arcs that cross the section in opposite directions.

### #494 (k₁,k₂) cycler construction — tube cross-section seeding
**Indirect but useful.** The intersection corridor at (x,y) ≈ (−μ, −0.93) for ALL C in 2.89–2.99
illustrates that the geometric location of the tube-tube intersection shifts very little with energy —
the corridor is essentially energy-independent. For #494 seed sizing: the window placement does
not need to track the Jacobi constant; a fixed geometric patch covers the relevant range. The
AREA of the intersection (tube cross-section) does scale with ΔC per the Ross-BozorgMagham 2018
law, but the LOCATION is stable.

### #314 `heteroclinic_cycle.py`
**Confirms methodology for single-system connections.** The L₁↔L₃ geometry uses the same
Floquet-seeding + Poincaré-section approach as #314's L₁↔L₂ implementation. The key difference
is L₃'s quasi-stability (lower unstable eigenvalues) and the need for the multi-crossing filter. The
manifold seeding parameters (ε_M = 10⁻⁶, N proportional to orbit amplitude) and corrector
tolerance (ε_L = 10⁻¹²) are directly reusable as reference values.

---

## 7. Corpus Status Recommendation
**`digested`** — no catalogue rows mined (this is an Earth-Moon CR3BP study, not a cycler paper).
Status: `digested`. Companion arXiv paper (braik-ross-2026-orbital-networks) already in corpus
and should be digested separately.
