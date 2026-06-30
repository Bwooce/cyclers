# Koon, Lo, Marsden, Ross 2002 — Constructing a Low Energy Transfer Between Jovian Moons (Contemp. Math. 292)

Digested 2026-06-30. Petit Grand Tour paper — the Ganymede-to-Europa manifold-tube
transfer construction. Filed in relation to #318 (joint moon-tour) and #500 (genome).

**Source (cite exactly, no file path):**
Koon, W.S., Lo, M.W., Marsden, J.E., Ross, S.D., "Constructing a Low Energy Transfer
between Jovian Moons," in Celestial Mechanics: Dedicated to Donald Saari for his 60th
Birthday, Contemporary Mathematics, Vol. 292, American Mathematical Society, Providence,
RI, 2002, pp. 129–145.
(Conference presentation: Celestial Mechanics, Evanston, IL, December 1999.)

> Text layer present. Full paper (14 pp. + references) read in its entirety. Key
> Poincaré section figures (Figs. 3.2, 2.4) are figure-only with key coordinates called
> out in captions; numerical values read from captions and surrounding text.

---

## 1. Contribution in three lines

Introduces the **"Petit Grand Tour"** (precursor to the MMO): starting beyond Ganymede,
a spacecraft is **temporarily captured by Ganymede** (one loop, ballistic), transfers to
**ballistic capture at Europa**, using a single maneuver at the manifold intersection
point. The construction uses the **coupled PCR3BP** (two nested Jupiter-moon-spacecraft
3-body problems) and finds intersections between invariant manifold tubes on a Poincaré
section, yielding a transfer ΔV of 1208 m/s vs. 2822 m/s for a Hohmann (57% savings).

---

## 2. Model and method

**Primary model:** Coupled Planar Circular Restricted 3-Body Problem (PCR3BP).
- Jupiter-Ganymede-spacecraft (µ_G = 7.802 × 10⁻⁵) for Ganymede leg
- Jupiter-Europa-spacecraft (µ_E = 2.523 × 10⁻⁵) for Europa leg
- Both coupled: 4-body bicircular model for verification; spacecraft in common orbital
  plane; Ganymede and Europa on coplanar circular orbits about Jupiter

**Mass parameters** (explicitly stated, eq. adjacent to Fig.2.1):

    µ_G = 7.802 × 10⁻⁵  (Jupiter-Ganymede)
    µ_E = 2.523 × 10⁻⁵  (Jupiter-Europa)

**Construction algorithm (§3):**
1. Compute unstable manifold of Ganymede L₁ Lyapunov orbit (interior region, amplitude
   0.35 in L-point-to-moon units, CG = 3.0576)
2. Compute stable manifold of Europa L₂ Lyapunov orbit (exterior region, amplitude 0.26
   in L-point-to-moon units, CE = 3.0028)
3. Find their intersection on the Poincaré section U4 in the Jupiter-Europa rotating
   frame (θ_E = θ_G = π, i.e. all 4 bodies on the x-axis at t=0)
4. Apply ΔV at the intersection to bridge the Jacobi-constant gap between the two
   manifold tubes
5. Integrate backward to Ganymede and forward to Europa in the full 4-body bicircular
   model

**Poincaré sections used** (§3, Fig.2.3):
- U1 = (y=0, x<0) in interior region
- U2 = (x=1-µ, y<0) in moon region
- U3 = (x=1-µ, y>0) in moon region
- U4 = (y=0, x<-1) in exterior region

---

## 3. Sourced goldens (exact numbers)

All numbers are from the numerical example in §3 unless noted.

### 3.1 Moon physical parameters

| Quantity | Value | Source |
|---|---|---|
| µ_G (Jupiter-Ganymede mass parameter) | **7.802 × 10⁻⁵** | p.4 (eq. label) |
| µ_E (Jupiter-Europa mass parameter) | **2.523 × 10⁻⁵** | p.4 (eq. label) |
| Ganymede orbital eccentricity | **0.0006** | p.4 text |
| Europa orbital eccentricity | **0.0101** | p.4 text |
| Inclination of orbital planes w.r.t. each other | **< 0.3°** | p.4 text |

### 3.2 Lyapunov orbit parameters (chosen for tube intersection)

| Quantity | Value | Source |
|---|---|---|
| Europa L₂ Jacobi constant (Jupiter-Europa frame) | **C_E = 3.0028** | p.10 |
| Europa L₂ Lyapunov orbit amplitude | **0.26** (in units of L₂-to-Europa distance) | p.10 |
| Ganymede L₁ Jacobi constant (Jupiter-Ganymede frame) | **C_G = 3.0576** | p.10 |
| Ganymede L₁ Lyapunov orbit amplitude | **0.35** (in units of L₁-to-Ganymede distance) | p.10 |

Note: these are chosen to produce a known tube intersection in position space (from
numerical experiments). "These particular Jacobi constants ... are chosen because they
are known, from numerical experiments, to lead to a dynamical channel intersection." (p.10)

### 3.3 Poincaré section intersection point and ΔV

| Quantity | Value | Source |
|---|---|---|
| Intersection point (x, ẋ) in Europa rotating frame | **(x ≈ -1.22, ẋ ≈ -0.005)** | p.11, Fig.3.2a caption |
| Jacobi-integral estimate of ΔV (y-direction) | **~1209 m/s** | p.11 text |
| Actual 4-body bicircular ΔV at patch point | **1208 m/s** | p.11 text, Fig.3.3 |
| Hohmann transfer Ganymede→Europa (baseline) | **2822 m/s** | p.11 text |
| Ratio (3-body / Hohmann) | **42.9%** | p.11 text |
| Transfer flight time (patch-point to Europa) | **~25 days** | p.12 text |

The agreement between the Jacobi estimate (1209) and the 4-body result (1208) is
< 0.1% — validates the coupled 3-body approximation for this system.

### 3.4 Additional maneuvers (total budget)

| Maneuver | ΔV | Purpose | Source |
|---|---|---|---|
| Near Ganymede (to add one loop) | **157 m/s** | Ganymede L₁ orbit too large; adds 1 revolution | p.12 text |
| At Europa closest approach | **87 m/s** | Near-circular capture orbit | p.12 text |
| **Total mission ΔV** | **1452 m/s** | Main ΔV + Ganymede loop + Europa circularize | p.12 text |

1452 m/s is still substantially less than Hohmann (2822 m/s).

### 3.5 Symbolic dynamics theorem (§2, p.6)

Existence theorem: for any admissible bi-infinite sequence of symbols {I, M, X} (Interior,
Moon, Exterior), there exists an orbit near the homoclinic/heteroclinic chain whose
itinerary matches. Numerical construction: truncate to finite itineraries by successive
Poincaré map intersections.

Example itinerary constructed: **(M, X; M, I, M)** in Jupiter-Ganymede system.

### 3.6 Earth-Moon "Shoot the Moon" comparison (Fig.4.1, §4)

Transfer from Earth orbit to ballistic lunar capture (same concept, Sun-Earth system):
- ΔV at Sun-Earth L₂ manifold intersection: **34 m/s** (Fig.4.1c label)
- (Belbruno-Miller 1991 Hiten mission used this class of trajectory)

This number is a reference for the concept applied to the E-M system, not a cycler row.

---

## 4. Key construction details for our codebase

**Phase choice for Ganymede relative to Europa:** At t=0, all four bodies (Jupiter,
Europa, Ganymede, spacecraft) aligned on x-axis in Europa rotating frame (θG' = 0).
This sets the relative phase of the moons at the patch point. General tours require
searching over the moon-moon phase angle as an additional parameter.

**Manifold tube cross-section shape:** On the Poincaré section, the tube cross-section
appears as a closed curve in (x, ẋ)-space. In (x, ẏ)-space, the Ganymede and Europa
manifolds lie on roughly constant lines (constant Jacobi integral in each system). The
gap in ẏ between the manifolds directly gives the ΔV to bridge them (eq. 2.2 applied).

**Jacobi constant gap:** The Ganymede L₁ channel has minimum Jupiter-Europa Jacobi value
3.0576 (significantly > Europa L₂ channel at 3.0028). The ΔV to close this gap:
    ΔV² ≈ Δ(C_J) × (scale factor from eq.2.2)
This estimate gave 1209 m/s vs. the actual 1208 m/s — extremely accurate.

**Itinerary notation:** The semicolon separates past from future; e.g. (X; M, I) = came
from exterior, in moon region, going to interior. Transit through a moon region is always
bounded by two Lagrange point crossings.

---

## 5. Reuse verdict vs. open threads

### #318 (joint moon-tour search — joint_cell.py, joint_sobol.py) — **HIGH**

This paper contains the **most algorithmic detail** of the three Ross-group papers for
implementing a manifold-tube moon-to-moon transfer leg. Specifically:

**Direct reuse:**
- The Poincaré section intersection algorithm (§3) is exactly what a "joint cell" leg
  between two moons computes: find intersection of Ganymede-L₁ unstable tube with
  Europa-L₂ stable tube at a chosen section angle θ
- The Jacobi-constant gap formula → ΔV estimate is a **fast analytical filter** for
  which moon pairs can possibly produce a transfer: |C_G - C_E| → ΔV ≥ threshold
- The bicircular 4-body refinement (starting from the 3-body intersection as initial
  guess) maps exactly to our multiple-shooting corrector
- The moon-phase angle θG' = θE - θG at t=0 is the key free parameter to search over;
  the paper shows the construction for θG' = 0

**For joint_sobol.py:** The Sobol search over initial conditions maps to searching over
(Lyapunov orbit amplitude, moon phase angle) for each tube-intersection; the two
parameters that determine whether an intersection exists are exactly C_G and C_E (which
fix the tube sizes).

**Key limitation:** The paper constructs only one tour example (Ganymede→Europa). The
generalization to arbitrary moon pairs requires:
- Computing Lyapunov orbits and their manifolds for arbitrary L₁/L₂ Jacobi constants
- Searching the 2-parameter family (C_G, C_E) for intersections
- Handling the moon-phase angle as an additional search dimension

These are all well-posed computations building on this paper's method.

### #465 (multi-rev leveraging) — **LOW-MODERATE**

This paper is primarily about the single-pass manifold tube intersection, not multi-pass
resonant GA. The Ganymede loop addition (157 m/s) is an impulsive fix, not a multi-rev
leveraging design. The resonant GA mechanism is in the MMO 2003 paper and Keplerian
map 2007 paper. However, the Poincaré section method here is compatible with multi-rev:
if multiple loops around Ganymede are desired, one selects from tube trajectories that
execute more revolutions before exiting.

### #500 (Keplerian-map genome) — **SUPPORTING**

This paper provides the **capture endpoint** that the Keplerian-map genome must reach.
Specifically:
- The Europa L₂ stable manifold tube is the capture condition that a Keplerian-map
  trajectory migrating from large a to near Europa must enter
- The Jacobi constant C_E = 3.0028 and the Europa L₂ tube define the specific
  (ω, a) exit region in the Keplerian map phase space for Europa capture
- For a multi-moon Keplerian map genome: the inter-moon ΔV (the 1208 m/s gap) represents
  the energy cost of **switching** from one moon's Keplerian map (at Jacobi constant C_G)
  to the next moon's map (at C_E). This is the "switching cost" in the genome.

### #494 (binary/circumbinary) — **NOT RELEVANT**

Single central body (Jupiter) throughout. No binary content.

---

## 6. Corpus-index status

**Verdict: `mined`** — full paper read; all sourced numbers extracted, construction
algorithm fully described. No catalogue rows (jovian moon tour, not an Earth-Mars cycler
or class-in-scope object). Cross-reference: cited as Koon et al. [2002] in Ross-Scheeres
2007 digest (reference [15] there).

---

## 7. Relationship between the three Ross-group papers

The three papers form a coherent suite:

| Paper | Role | Algorithmic depth |
|---|---|---|
| Koon-Lo-Marsden-Ross 2002 (this paper) | Manifold tube intersection, single-ΔV inter-moon transfer | HIGH — complete construction algorithm |
| Ross-Koon-Lo-Marsden 2003 (MMO) | Extension to 3-moon sequence with resonant GA | MEDIUM — building blocks described, example shown |
| Ross-Scheeres 2007 (Keplerian map) | Analytical map for resonant GA dynamics | HIGH — fully derived, map equations extractable |

For our project: Koon 2002 gives the **capture/escape construction** (one leg); MMO 2003
gives the **multi-leg tour architecture** (sequence of legs); Keplerian map 2007 gives
the **fast surrogate** for the resonant GA inter-moon transfer sub-legs.

---

## 8. Not tabulated / figures-only

- Fig.2.4 (tube intersection in U3 section): (x, ẋ) position of intersection labeled
  approximately in caption; exact (x, ẋ, ẏ) at intersection not in text
- Fig.3.2b (x, ẏ) Poincaré section: manifolds lie on roughly constant ẏ lines;
  the ẏ gap visually corresponds to 1209 m/s but not tabulated as a (ẏ_G, ẏ_E) pair
- Fig.4.1 (Shoot the Moon): Sun-Earth L₂ manifold cross-sections labeled qualitatively;
  ΔV = 34 m/s is stated in the caption to panel (c)
