# Ross & Roberts-Tsoukkas 2025 (AAS 25-621) — stable ballistic prograde Earth-Moon cyclers (mining)

**Date:** 2026-06-11
**Source:** S. D. Ross and M. Roberts-Tsoukkas, "Stable, Low-Energy Prograde
Earth-Moon Cycler Orbits," AAS 25-621, AAS/AIAA Astrodynamics Specialist
Conference, 2025 (16 pp., LaTeX preprint dated 2025-07-17; publicly posted at
https://ross.aoe.vt.edu/papers/ross-roberts-tsoukkas-2025-AAS-25-621.pdf).
All table/page/equation numbers below are from that AAS manuscript. A 2026
journal version exists ("Stable Prograde Earth-Moon Multi-Orbiter Cyclers via
Three-Body Dynamics," Roberts-Tsoukkas & Ross) — see §10; re-verify numbering
against it before embedding citations in code comments.
**Writeback: NONE.** Extraction + proposed rows only; catalogue edits are
review-gated. No code changes in this pass (one code bug *found*, flagged §7).

**Verdict: USEFUL (top of the sweep, as ranked).** This is the first
publication of STABLE, fully ballistic, prograde Earth-Moon cyclers — exactly
the class our CR3BP lane lacks (our own Saturnian discoveries were unstable
Lyapunovs, error ×~2000/period). The paper delivers (1) five (k1,k2)-cycler
families with 15–16-digit Jacobi constants and periods for the stable
representative of each, (2) a complete, implementable construction method
(tube-overlap seeding → fixed-Jacobi symmetric differential correction →
half-period monodromy stability → pseudo arc-length continuation), and (3) six
high-precision critical-tangency constants usable as manifold-tube goldens.
The honest gap: **no full initial state vector is printed for any member** —
but the symmetric-orbit structure makes every IC recoverable from the
published (μ, C, T) by a one-dimensional solve (§5). Mining also surfaced a
**real μ bug in our `cr3bp_system("Earth","Moon")`** (−1.2%, masked by a test
tolerance) — see §7.

---

## 1. Model (lens 2 — binding)

Everything in the paper is **pure planar CR3BP (PCR3BP)** in the standard
barycentric rotating frame. No bicircular/ephemeris refinement anywhere in
this AAS version (solar/eccentricity perturbations are explicitly future work,
Conclusion p. 14–15).

- Mass ratio (p. 3): **μ = 1.2150584270572 × 10⁻²** (Earth at (−μ, 0), Moon at
  (1−μ, 0)). This equals GM_Moon / GM_Earth-system for the JPL DE440 values
  (4902.800118 / 403503.2356254802 = 0.01215058439…, matches to 1e-8 relative).
- Length unit (p. 3): a_m = **384,400 km** (Earth-Moon distance = 1).
- Time unit (p. 3): sidereal month T_m = **27.321661 days = 2π TU**, i.e.
  1 TU = 27.321661/(2π) = **4.348377401631 days** (≈ 375,699.8 s). Verified:
  Table 4's TU and day columns agree under exactly this conversion.
- Jacobi convention (Eq. 3, p. 3): C = −2Ū − (ẋ²+ẏ²) with
  Ū = −½(x²+y²) − (1−μ)/r1 − μ/r2, i.e.
  **C = (x²+y²) + 2(1−μ)/r1 + 2μ/r2 − v²**, explicitly *omitting* the
  μ(1−μ) additive constant some authors use to normalize C(L4) = 3.
  **This is exactly our `jacobi_constant` convention** in
  `src/cyclerfinder/core/cr3bp.py` — Jacobi values transfer with no offset.
- Equilibrium values (p. 13 text): **C1 = C(L1) = 3.188341105401253**,
  **C2 = C(L2) = 3.172160450399808** (paper's μ). Hill-region cases (Fig. 3,
  pp. 3–4): Case 2 (C1 > C > C2) opens the L1 neck only; Case 3 (C2 > C > C3)
  also opens L2 to the exterior realm.
- Moon radius for collision termination (p. 10): **R_m = 1,740 km**.

Orbit class context: the orbits are periodic in the rotating frame,
alternately encircling Earth and Moon (temporary capture each cycle),
**prograde**, and classified per Broucke (1968, TR 32-1168) symmetric-orbit
classes 2, 3, or 5 depending on crossings (p. 8).

## 2. (k1,k2)-cycler classification (Def. 1, p. 5)

Four Poincaré sections on the x-axis (p. 5, Fig. 4):

- U1⁻ = {(x, ẋ) | y = 0, x < −μ, ẏ < 0} (Earth realm, descending)
- U1⁺ = {(x, ẋ) | y = 0, x > −μ, ẏ > 0} (Earth realm, ascending)
- U2⁻ = {(x, ẋ) | y = 0, x < 1−μ, ẏ < 0} (Moon realm)
- U2⁺ = {(x, ẋ) | y = 0, x > 1−μ, ẏ > 0} (Moon realm)

**Definition 1 (p. 5):** a (k1,k2)-cycler is a periodic solution crossing U1⁻
exactly k1 times and U2⁺ exactly k2 times per cycle (k1, k2 ∈ ℤ⁺). The sign
conventions select prograde motion. For symmetric cyclers, perpendicular
crossings may fall on either signed section of each pair (p. 8).

## 3. Construction method (lens 3 — implementable end-to-end)

1. **Tube seeding (pp. 6–8).** For C < C1 the L1 Lyapunov orbit's
   stable/unstable manifold tubes are co-dimension-one separatrices; their
   cuts on U1 (Earth branch) / U2 (Moon branch) bound all realm-transit
   trajectories (Koon et al. 2000 Global Orbit Structure Theorem, their
   ref [13] — the paper's framing: stable cyclers live *inside* the tube
   region but beyond the prior unstable constructions). As C decreases below
   C1 the cut areas grow ∝ ΔC = C1 − C (their ref [16]). At a critical
   C^{u1}_{k1} the k1-th stable- and unstable-manifold cuts first become
   tangent on the symmetry line {ẋ = 0}, at location x^{u1}_{k1} — the first
   possible emergence of symmetric (k1,·)-cyclers (Fig. 6; tangency for
   k1 = 1 at x ≈ −0.768). Same construction Moon-side gives (C^{u2}_{k2},
   x^{u2}_{k2}). Symmetry-reduced candidate sets: S^{u1}_{k1} (Eq. 6),
   S^{u2}_{k2} (Eq. 7) — intersections of the cut interiors with {ẋ = 0}.
2. **Upper bound (Eq. 8, p. 8):** C_(k1,k2) ≡ min{C^{u1}_{k1}, C^{u2}_{k2}}.
   Below it, flow S^{u1}_{k1} forward to U2 and intersect with S^{u2}_{k2}
   (Γ^{u2} = P(S^{u1}_{k1}) ∩ S^{u2}_{k2}); pull back to U1 for IC guesses.
3. **Differential correction at fixed Jacobi constant (pp. 9–10).** Symmetric
   IC (x0, 0, 0, ẏ0) with ẏ0 = ẏ0(x0, C) from C(x0, 0, 0, ẏ0) = C (Eq. 9);
   time-reversal symmetry s_x: (x, y, ẋ, ẏ, t) → (x, −y, −ẋ, ẏ, −t) (Eq. 10).
   Integrate to the first crossing near S^{u2}_{k2} at t1 and require a
   perpendicular crossing y(t1) = 0, ẋ(t1) = 0 (Eq. 11) — a 1-unknown
   Newton-Raphson in x0. Update (Eq. 12, as printed; **re-verify subscripts
   against the journal before implementing** — render is small):

   ```
   δx0 = ẋ1 · ( ẏ1 / (ẍ1 Φ21) ) ·
         [ 1 − (Φ24/Φ21)·(1/ẏ0ᵍ(n))·∂Ū/∂x(0)
             − (ẏ1/ẍ1)·(1/Φ21)·( Φ31 − Φ34·(1/ẏ0ᵍ(n))·∂Ū/∂x(0) ) ]⁻¹
   ```

   with ẍ1 from the EOM (1) at the crossing and ∂Ū/∂x(0) at the initial
   point; then re-enforce C: ẏ0ᵍ(n+1) = ±√(−2Ū(x0ᵍ(n+1), 0) − C).
   Tolerance ε = 10⁻⁸ nd (~10⁻⁵ m/s); convergence typically n ≤ 5
   iterations; period **T = 2 t1** (p. 9).
4. **Stability via half-period monodromy (Eqs. 13–15, p. 10).**
   M = Φ(T, T/2) Φ(T/2, 0) with Φ(T, T/2) = G Φ(T/2, 0)⁻¹ G,
   G = diag(1, −1, −1, 1) (Barden 1994, their ref [18]) — no full-period
   integration needed. Nontrivial eigenpair (λ, 1/λ); **stability parameter
   ν ≡ ½(λ + 1/λ)** (Eq. 15); linearly stable iff |ν| < 1 (complex-conjugate
   unit-circle pair); ν = 0 ⇔ λ = ±i defines the stable-subfamily midpoint.
5. **Family continuation (Eq. 16, p. 10):** pseudo arc-length continuation
   (their refs [23, 24], incl. Zhong & Ross 2021) traces the one-parameter
   family F(k1,k2) = {(x(s), C(s))} while monitoring ν; endpoints s_min/s_max
   may terminate in lunar-surface collision (r_p < R_m, R_m = 1,740 km).

## 4. Published numbers — complete transcription (lens 1)

### Table 1 (p. 8) — Earth-side critical Jacobi constants and tangency x (U1⁻ crossings)

| k1 | C^{u1}_{k1} | x^{u1}_{k1} |
|----|-------------------|-----------------|
| 1 | 3.151763728314920 | −0.767856324800 |
| 2 | 3.129751730201047 | 0.723754610150 |
| 3 | 3.188341092440989 | −0.332153924455 |

Text note (p. 8): C^{u2}_{k2} increases monotonically with k2; the Earth-side
C^{u1}_{k1} are *not* monotonic. Note C^{u1}_3 sits only 1.3e-8 below C1.

### Table 2 (p. 8) — Moon-side critical values (U2⁺ crossings)

| k2 | C^{u2}_{k2} | x^{u2}_{k2} |
|----|-----------------|--------------|
| 1 | 3.1833333078762 | 1.0016252150 |
| 2 | 3.1840565764573 | 0.8611415325 |
| 3 | 3.1845534633380 | 1.0110341410 |

### Table 3 (p. 11) — Earth-Moon prograde (k1,k2)-cycler periodic orbits

C_(k1,k2) = theoretical upper bound (Eq. 8); C^max = Jacobi at which the
family actually first emerges (C^max < C_(k1,k2) < C1); C^stable = the
ν = 0 midpoint of the largest stable subfamily; T^stable its period;
Δp_m = perilune-distance width of the largest stable subfamily (km).

| k1 | k2 | C_(k1,k2) | C^max_(k1,k2) | C^stable_(k1,k2) | T^stable (TU) | Δp_m (km) |
|----|----|------------------|-------------------|-------------------|--------------------|---------|
| 1 | 1 | 3.1517637283600 | 3.151175879916394 | 3.151175879508174 | 10.29206921007976 | 0.13 |
| 2 | 1 | 3.1297495000000 | 3.129389531092325 | 3.129389531088256 | 19.44043166795154 | 4.23 |
| 3 | 1 | 3.1833333078762 | 3.161796247265416 | 3.161784147013429 | 14.78849241668140 | 253.70 |
| 3 | 2 | 3.1840565764573 | 3.182762785398336 | 3.182762663084288 | 17.90058010350006 | 42.08 |
| 3 | 3 | 3.1845534633380 | 3.183379082936385 | 3.177224018696528 | 18.14546057589189 | 2041.34 |

Period conversions (1 TU = 4.348377401631 d): 44.7538 / 84.5343 / 64.3059 /
77.8385 / 78.9033 days respectively (rows 1–3 match Table 4's day column
exactly).

### Table 4 (p. 13) — (k1,1) summary with deltas (Eq. 17 definitions)

ΔC_(k1,1) = C1 − C_(k1,1); ΔC^max = C_(k1,1) − C^max; ΔC^stable = C^max − C^stable.

| k1 | ΔC_(k1,1) | ΔC^max_(k1,1) | ΔC^stable_(k1,1) | T^stable (TU) | T^stable (days) | Δp_m (km) |
|----|--------------|--------------|--------------|-----------|------------|--------|
| 1 | 3.657738e-02 | 5.878484e-04 | 4.082197e-10 | 10.292069 | 44.753800 | 0.13 |
| 2 | 5.859161e-02 | 3.599689e-04 | 4.069189e-12 | 19.440432 | 84.534335 | 4.23 |
| 3 | 1.272710e-02 | 1.381776e-02 | 1.210025e-05 | 14.788492 | 64.305944 | 253.70 |

Energy regime (p. 13 text): all (k1,1) examples satisfy C_(k1,1) < C2 =
3.172160450399808 but C > C3 — **Case 3** (exterior escape energetically
open). Our derived observation (not stated in the paper): the (3,2) and (3,3)
stable members have C^stable > C2 (3.18276…, 3.17722…) — **Case 2**, exterior
escape energetically closed for those.

### Per-family text/figure values

- **(1,1)** (p. 11, Fig. 1): stable subfamily extremely narrow, perilune width
  ≈ 0.1 km (Δp_m = 0.13 km). Midpoint period T^stable = **44.7538 days**,
  within 1% of the 2:3 resonance with the Earth-Moon **synodic** month
  29.530588 d (2×44.7538 = 89.508 vs 3×29.530588 = 88.592; we compute 1.03%) —
  near-repeating Sun-Earth-Moon-spacecraft geometry. Fig. 1 label: 45 days.
- **(2,1)** (p. 11, Fig. 7): stable subfamily slightly wider (4.23 km).
  T^stable = 84.534 d ≈ **3.09 × 2π TU** — not an integer multiple of the
  sidereal month, so the orbit does not close in the Earth-centered inertial
  frame (Fig. 7b). Fig. 1 label: 84 days.
- **(3,1)** (p. 12, Fig. 8): family parameterized by perilune altitude, 0 km
  to > 15,000 km. **Stable window ≈ 750–1,000 km perilune altitude**
  (Δp_m = 253.70 km; nd perilune radius ≈ 0.00648–0.00713). T^stable =
  64.3059 d (Fig. 1 label: 64 days). Fig. 8(b) U2 Poincaré insets give
  **basin-of-stability widths** (largest surrounding resonant torus, in km):
  23, 195 (within window), and 155, 146, 92, 227 along the family; the
  bottom-left inset shows a *linearly unstable* member still surrounded by a
  stable torus (global stability without linear stability). Entry/exit of the
  stable window marked by a **period-doubling bifurcation** (ν = −1 side) and
  a **tangent bifurcation** (ν = +1 side); further interior bifurcations at
  classical 1/4 and 1/2 branching ratios (Greene's criterion, their ref [25]).
  Basin width is reported only for the (3,1) family (p. 11).
- **(3,2)** (p. 13, Fig. 9): **two distinct stability windows**, the larger
  ≈ 40 km wide in perilune altitude (Δp_m = 42.08 km;
  C^max = 3.182762785 per Fig. 9a label). Stable orbits' period close to the
  **2:5 synodic resonance** (5×29.530588/2 = 73.826 d); the resonance is
  satisfied *exactly* by two nearby **unstable** members with stability
  parameters of order 10–100. Fig. 9b representative: **74 days** — note this
  is the near-resonant window member, NOT the Table-3 ν = 0 member (whose
  T^stable converts to 77.84 d). A member of this family was previously
  identified by Leiva & Briozzo (their refs [10, 11]; p. 13 text spells
  "Leiva and Brizzolara") and shown to **persist under solar perturbation**
  (Sun-Earth-Moon quasi-bicircular model) — the only prior sighting, as a
  single unstable orbit.
- **(3,3)** (p. 14, Fig. 10): **five distinct stability windows**; the largest
  spans perilune altitudes ≈ **4,200–6,200 km** (Δp_m = 2041.34 km),
  corresponding to perigee altitudes ≈ **112,400–113,500 km** (~3× GEO
  distance). Family terminates at both ends in lunar-surface collision. Two
  smaller windows near C^max. Fig. 10 axis labels: C^max = 3.183379083;
  Fig. 10(b) gives the family's initial condition curve in (x0, C): x0
  (= perigee crossing on U1⁻) runs ≈ **−0.330 to −0.315 nd**. Generic
  observation (p. 14): *every* (k1,k2) family examined has a stable region
  near C^max (generic Hamiltonian bifurcation structure, their ref [26],
  Golubitsky et al.).

### Internal cross-table inconsistencies (flag for journal re-check)

We verified Tables 1–4 against each other (arithmetic, 2026-06-11):

1. Table 3's C_(2,1) = 3.1297495000000 (trailing zeros suspicious) disagrees
   with Eq. 8 + Table 1: min{3.129751730201047, 3.1833333078762} =
   3.129751730201047 (Δ = 2.23e-6). Table 4's ΔC_(2,1) is consistent with the
   *Table 3* value.
2. Table 3's C_(3,1) = 3.1833333078762 (= C^{u2}_1, the Eq. 8 min) is
   inconsistent with Table 4: C1 − ΔC_(3,1) = 3.188341105401253 − 1.272710e-2
   = 3.175614005, and ΔC^max_(3,1) = 1.381776e-02 only reconciles with
   *that* value (3.175614005 − 3.161796247265416 = 1.3817758e-2). So Table 4
   implies C_(3,1) = 3.1756140, Table 3 prints 3.1833333 (Δ = 7.7e-3).
3. Table 3's C_(1,1) = 3.1517637283600 vs Table 1's C^{u1}_1 =
   3.151763728314920 — differ at the 11th digit (benign rounding).

Rows (1,1)/(2,1) of Tables 3 and 4 are mutually consistent; all ΔC^stable
values check exactly. Resolve against the 2026 journal version before any of
the C_(k1,k2) *bound* columns is used as a golden; the **C^max, C^stable,
T^stable columns are internally consistent** and are the catalogue-grade
numbers.

## 5. Completeness of (μ, state0, period, Jacobi) tuples

Census: **5 families, 1 stable representative each** (the ν = 0 midpoint),
plus 6 critical-tangency constants and 2 figure-level family curves.

| Item | μ | C | T | state0 | Status |
|---|---|---|---|---|---|
| (1,1), (2,1), (3,1), (3,2), (3,3) stable members | exact (printed) | 15–16 digits | 16 digits (TU) + days | **not printed** | complete up to a 1-D solve |
| Tangency constants (Tables 1–2) | exact | 13–16 digits | n/a | x on section printed | complete as section-geometry goldens |
| (3,3) family curve (Fig. 10b) | exact | curve | n/a | x0 ∈ [−0.330, −0.315] (figure-read) | seed bracket only |

The missing state0 is *structurally* recoverable: every catalogued member is a
symmetric cycler with IC of the form **(x0, 0, 0, ẏ0)** where ẏ0 follows from
the published C (Eq. 9). Fixing C = C^stable, the only unknown is x0, and the
published T^stable provides the independent acceptance check (T = 2t1 must
match to the printed 16 digits). This is far better than the usual V0
publication gap: the tuple closes with a one-dimensional root-solve in the
same model. (Golden discipline: the EXPECTED values are the published
(μ, C^stable, T^stable); any x0 we recover is a *derived* quantity stored with
provenance, never itself a golden.)

## 6. Relationship to known orbits — dedup/supersession (lens 4)

- **Arenstorf cyclers** (their ref [1]; our `arenstorf-em-figure8-1963`):
  cited as prior art; the new class differs (prograde, temporary capture at
  both primaries, stable subfamilies). **No overlap, no supersession.** Note
  Arenstorf demos use the test-problem μ = 0.012277471; this paper uses the
  physical 0.012150584270572 — same-model golden rule applies per row.
- **Aldrin cyclic concepts** (ref [2]) and **Uphoff & Crouch 1993** (ref [3]):
  background only.
- **Genova & Aldrin 2015** (their ref [6]; our `genova-aldrin-2015` row,
  reclassified bicircular 2026-06-10): the (3,1) cycler's rotating-frame
  geometry (3 Earth-realm petals + lunar loop, Fig. 1) superficially resembles
  the Genova/Aldrin 3-petal cycler, but **they are not the same orbit**: the
  3-petal Genova orbit does not exist in the pure CR3BP (needs solar gravity
  + station-keeping), whereas the (3,1) family is a genuine ballistic PCR3BP
  periodic orbit with a stable window. The (3,1) row would be the first
  *pure-CR3BP* 3-petal-geometry entry; keep both, cross-reference, do not
  merge.
- **Leiva & Briozzo 2006/2008** (refs [10, 11]): identified a single
  *unstable* orbit resembling the (3,2) class and showed persistence under
  solar perturbation; this paper supplies the continuous family + stable
  windows around it. If we ever ingest Leiva-Briozzo, its member should be
  parented to the (3,2) family.
- **Wittal/Miaule/Asher IAC-22** (their ref [9]; our `wittal-2022` family-seed
  row): cited as prior Earth-Moon cycler family work; different construction;
  no number overlap.
- **Casoliva et al. 2008/2010** (refs [4, 5], "two classes of cycler
  trajectories in the Earth-Moon system"): nearest prior family-level work;
  not in our catalogue; acquisition candidate if we want the unstable
  complement.
- Theoretical footing: Koon-Lo-Marsden-Ross 2000 (ref [13]) heteroclinic/
  Global Orbit Structure Theorem; this paper's stable cyclers sit inside the
  L1/L2 tube-enclosed region but beyond the theorem's unstable-orbit
  neighborhood — a genuinely new regime, not a rediscovery.
- Companion AAS paper: Braik & Ross, AAS 25-716, L1–L3 heteroclinic transfer
  (ref [17]) — source of the differential-correction loop details.

## 7. Side-finding: registry μ bug in `cr3bp_system("Earth","Moon")` (action: fix + retest)

While checking the paper's μ against ours:
`PRIMARIES["Earth"]` in `src/cyclerfinder/core/satellites.py` is the JPL
**Earth-system GM** (Earth+Moon, 4.0350323562548019e5 km³/s²), but
`cr3bp_system()` computes μ = GM2/(GM1+GM2) — **double-counting the Moon**:

- produced: μ = 4902.800118/(403503.2356+4902.800118) = **0.0120047** (−1.20%)
- correct: μ = 4902.800118/403503.2356254802 = **0.01215058439**, which
  matches this paper's 0.012150584270572 to 1.0e-8 relative.

`tests/core/test_cr3bp.py::test_earth_moon_mu_physical` asserts
`mu ≈ 0.01215 (abs=2e-4)` and the error is 1.46e-4 — **the tolerance masks the
bug**. Affected consumers: `scripts/cr3bp_backfill.py` (Earth-Moon backfill of
2026-06-10, which even prints "~0.01200 (JPL GMs)" without noticing),
and any Earth-Moon use of `cr3bp_moontour_run.py`. The Arenstorf golden is
unaffected (explicit μ = 0.012277471). The same double-count pattern exists
for Jupiter/Saturn systems but the relative error there is ≤ ~2e-4 (system GM
vs planet GM). Also note `t_s` inherits the error (≈ 0.6% on the Earth-Moon
time unit: 372,931 s vs ≈ 375,700 s). **Fix is one line plus a tightened test
tolerance; Earth-Moon backfill results of 2026-06-10 need a re-run after the
fix.** Review-gated like everything else in this note.

## 8. Proposed catalogue rows (NO writeback — review-gated)

Five Tier-2 Earth-Moon rows, one per family, anchored on the stable (ν = 0)
representative. Common fields: `model_assumption: cr3bp`;
`mass_ratio: 1.2150584270572e-2` (printed, p. 3); frame: Earth-Moon rotating,
barycentric, planar; `center`: Earth-Moon barycenter; `source_ephemeris`: n/a
(pure CR3BP, no ephemeris dependency); length unit 384,400 km; time unit
sidereal month/2π = 4.348377401631 d (p. 3); Jacobi convention = ours;
`tof_days_bounds`: [period, period]. Source cite for every number:
Ross & Roberts-Tsoukkas 2025, AAS 25-621, Table 3 (p. 11) + per-family text.

| proposed id | family | C^stable | period (TU / days) | encounter cadence | stable-window data |
|---|---|---|---|---|---|
| ross-rt-em-cycler-11-2025 | (1,1) | 3.151175879508174 | 10.29206921007976 / 44.753800 | 1 perigee + 1 perilune pass per period | Δp_m 0.13 km; ~2:3 synodic resonance (p. 11) |
| ross-rt-em-cycler-21-2025 | (2,1) | 3.129389531088256 | 19.44043166795154 / 84.534335 | 2 perigee + 1 perilune | Δp_m 4.23 km; 3.09 sidereal months (p. 11) |
| ross-rt-em-cycler-31-2025 | (3,1) | 3.161784147013429 | 14.78849241668140 / 64.305944 | 3 perigee + 1 perilune | window 750–1,000 km perilune alt; Δp_m 253.70 km; torus widths 23–227 km (Fig. 8) |
| ross-rt-em-cycler-32-2025 | (3,2) | 3.182762663084288 | 17.90058010350006 / 77.838478* | 3 perigee + 2 perilune | Δp_m 42.08 km; 2 windows; near 2:5 synodic (p. 13) |
| ross-rt-em-cycler-33-2025 | (3,3) | 3.177224018696528 | 18.14546057589189 / 78.903311* | 3 perigee + 3 perilune | window 4,200–6,200 km perilune alt / 112,400–113,500 km perigee alt; Δp_m 2041.34 km; 5 windows (p. 14) |

(* day values are our conversion via the printed time unit; Table 4 prints
days only for the (k1,1) rows — mark the (3,2)/(3,3) day figures as derived.)

Row-level extras worth carrying: `stability: linearly stable (|ν|<1, ν=0 at
listed member; Eq. 15 convention ν = ½(λ+1/λ))`; family C^max values
(Table 3); the two bifurcation types bounding the (3,1) window (Fig. 8b);
Case 2/Case 3 energy-regime flag (C2 = 3.172160450399808). Validation level
at ingest: the (μ, C, T) tuple is sourced and complete; state vector pending
our 1-D recovery ⇒ start at the family-seed level (like `wittal-2022`) and
promote once the corrector closes on (C^stable, T^stable) — that closure is a
*same-model* reproduction, eligible to raise validation per the ladder rules.

Also worth a non-orbit registry entry (or `data_gaps` note): the six
critical-tangency constants (Tables 1–2) as method-level goldens.

## 9. CR3BP-lane seeds and goldens

System construction: **do not use `cr3bp_system("Earth","Moon")` until §7 is
fixed**; build `CR3BPSystem(mu=1.2150584270572e-2, l_km=384400.0,
t_s=375699.8)` explicitly (paper's μ verbatim — same-model golden target).

Seeds (all ICs of the form (x0, 0, 0, ẏ0(x0, C)); T_guess from Table 3):

1. **(1,1):** C = 3.151175879508174, x0 bracket around −0.7679 (tangency
   x^{u1}_1 = −0.767856324800; C^stable is only 5.88e-4 below C^{u1}_1, so the
   perigee crossing stays near the tangency point). T_guess = 10.292069 TU.
2. **(2,1):** C = 3.129389531088256, x0 bracket around 0.7238
   (x^{u1}_2 = 0.723754610150 — note positive x, U1⁺-side perpendicular
   crossing). T_guess = 19.440432 TU.
3. **(3,1):** C = 3.161784147013429, x0 near −0.332 (x^{u1}_3 =
   −0.332153924455). Secondary check: perilune altitude must land in
   750–1,000 km (r_p ≈ 0.00648–0.00713 nd). T_guess = 14.788492 TU.
4. **(3,2):** C = 3.182762663084288, x0 near −0.332. T_guess = 17.900580 TU.
5. **(3,3):** C = 3.177224018696528, x0 ∈ [−0.330, −0.315] (Fig. 10b);
   perilune 4,200–6,200 km. T_guess = 18.145461 TU.

Corrector route: our `correct_periodic` (7-unknown full-state shooting) can
consume these directly, but the clean implementation is the paper's
**fixed-Jacobi symmetric corrector** — 1 unknown (x0), constraints
y(t1) = ẋ(t1) = 0 at the half period, C enforced algebraically (§3 step 3).
That variant also buys the Barden half-period monodromy
(M = G Φ(T/2)⁻¹ G Φ(T/2)) for stability at half the integration cost, and ν
for stable-window scans — the capability our lane lacks (everything we have
found so far is unstable; a |ν| < 1 verdict has never been exercised).

Acceptance (sourced-golden rule, EXPECTED = published values only):

- period: corrected T vs Table 3 T^stable (16 digits printed; expect
  agreement limited by our integrator tolerance, target ≤ 1e-10 nd);
- Jacobi: re-evaluated C(state0) vs C^stable (consistency, since C is
  enforced);
- stability: |ν| < 1 at the recovered member, and ν ≈ 0 (midpoint property,
  Table 3 definition p. 11);
- topology: k1 U1⁻ crossings + k2 U2⁺ crossings per period (Def. 1) — the
  orbit-closure-discipline "verify topology vs source first" step;
- windows: perilune-altitude window membership for (3,1)/(3,3).

Method-level goldens (independent of any orbit row): C^{u1}_1 =
3.151763728314920 with tangency at x = −0.767856324800 (Table 1) — a
15-digit target for an L1-Lyapunov manifold-tube cut + tangency-detection
routine, which our lane's Lyapunov family machinery is one step away from.
Avoid the Table-3 C_(k1,k2) *bound* column as a golden until the §4
inconsistencies are resolved against the journal.

## 10. What the 2026 journal version should add (lens 5)

Stated future work (Conclusion, pp. 14–15): **non-symmetric cyclers; 3D
families; perturbation models (lunar eccentricity, solar gravity); extension
to the exterior realm** (cyclers reaching toward the Earth SOI edge; asteroid
capture relevance, their ref [29]); spacecraft deployment/transfer planning;
and the open question of basin-width ≈ stable-subfamily-width ("a
relationship we intend to investigate further", p. 14). The journal title in
our sweep ("Stable Prograde Earth-Moon Multi-Orbiter Cyclers via Three-Body
Dynamics") suggests multi-spacecraft/constellation framing. Expect: possibly
more (k1,k2) combinations ("we explored several, but certainly not all",
p. 14), hopefully printed initial conditions, resolution of the §4 table
inconsistencies, and persistence results under solar/eccentricity
perturbation (which would gate any bicircular/ephemeris-fidelity claims —
cross-fidelity rule: until then, every number here is **PCR3BP-only**).

## 11. Reference leads worth pulling later

- [4]/[5] Casoliva, Mease, Mondelo, Villac, Barrabes, Olle (2008 AAS / 2010
  JGCD 33(5) 1623–1640): "families of cycler trajectories in the Earth-Moon
  system" — the unstable-cycler family predecessor; candidate acquisition.
- [10]/[11] Leiva & Briozzo (Acta Astronautica 58(8) 2006 379–386; CMDA 101
  2008 225–245): the (3,2)-resembling unstable member + quasi-bicircular
  persistence — the bridge to higher-fidelity validation of the (3,2) row.
- [12] Broucke, NASA TR 32-1168 (1968): symmetric-orbit classes 2/3/5; the
  classical Earth-Moon periodic-orbit census these families thread through.
- [17] Braik & Ross, AAS 25-716 (2025): the differential-correction loop
  source; companion paper, likely has its own golden-grade constants.
- [18] Barden, Purdue M.S. thesis (1994): half-period monodromy identity —
  cite when implementing Eq. 13–14.
- [24] Zhong & Ross, Applied Mathematical Modelling 97 (2021) 81–95:
  differential correction + arc-length continuation how-to (implementation
  reference for the continuation step).

---

## TRANSCRIPTION RESCAN (2026-06-12)

Character-by-character re-read of the source PDF (Ross & Roberts-Tsoukkas
2025, AAS 25-621) against this note and every repo artifact that inherited its
numbers (the 5 `ross-rt-em-cycler-*-2025` catalogue rows,
`tests/search/test_cr3bp_ross_families.py`,
`tests/core/test_cr3bp.py::test_earth_moon_mu_physical`, and the
`_LEVEL_EVIDENCE` strings in `cyclerfinder/data/validate.py`). Each
load-bearing table was read TWICE in independent passes; all arithmetic
cross-checks re-run exactly (Decimal).

**VERDICT: the golden set is CLEAN.** Every 15–16-digit value below is a
three-way MATCH (PDF = note = repo):

- **μ = 1.2150584270572e-2** (p. 3, read twice) — matches §1, `ROSS_MU`, all
  five rows' `mass_ratio`, and the `test_earth_moon_mu_physical` expected.
- **Table 3 (p. 11), C^stable and T^stable, all 5 families, all digits** —
  match §4/§8, `_FAMILIES`, the rows' `jacobi_constant`/`period_nd`, and the
  `_LEVEL_EVIDENCE` strings. THE critical set: no digit slips anywhere.
- **Table 3 C^max, all 5** — match §4 and the rows' `jacobi_max`.
- **Table 3 Δp_m** (0.13 / 4.23 / 253.70 / 42.08 / 2041.34 km) — match; the
  (3,1) 750–1,000 km window (p. 12 + Fig. 8 caption), the (3,3)
  4,200–6,200 km perilune / 112,400–113,500 km perigee windows (p. 14), and
  the Fig. 8(b) basin widths 23/195/155/146/92/227 km all confirmed.
- **Tables 1–2 (p. 8), all 6 tangency constants + x-locations** — match §4
  and `ROSS_TABLE1_CU1`.
- Supporting constants confirmed: C1 = 3.188341105401253 and
  C2 = 3.172160450399808 (p. 13); R_m = 1,740 km (p. 10); ν = ½(λ+1/λ),
  G = diag(1,−1,−1,1) (p. 10); a_m = 384,400 km, T_m = 27.321661 d (p. 3);
  synodic 29.530588 d, 3.09×2π (p. 11); Fig. 9(a) label 3.182762785 and
  Fig. 9(b) "74 days" (p. 13); Fig. 10 label 3.183379083 and x0 axis
  −0.33…−0.315 (p. 14); Table 4 (p. 13) all 18 entries as transcribed in §4.

**The two flagged C-bound conflicts are CONFIRMED SOURCE-INTERNAL (PDF
genuinely self-inconsistent; not our misreading — both sides re-read twice):**

1. C_(2,1): Table 3 prints 3.1297495000000; Eq. 8 + Table 1 give
   3.129751730201047. Exact Δ = **2.2302e-6** (this note's §4 had rounded it to
   2.24e-6; correct 3-sig-fig value is 2.23e-6 — prose-only slip, also echoed
   in the (2,1) catalogue row's data_gap note; BOTH FIXED 2026-06-12 in the
   rescan-cleanup pass; the pinning test's
   2.235e-6 ± 1e-7 assertion brackets the true value and passes).
   Table 4's ΔC_(2,1) = 5.859161e-2 sides with the Table 3 value
   (C1 − 3.1297495 = 5.8591605e-2; C1 − eq8 = 5.858938e-2). As §4 said.
2. C_(3,1): Table 3 prints 3.1833333078762 (= C^{u2}_1, the correct Eq. 8
   min); Table 4's ΔC_(3,1) = 1.272710e-2 and ΔC^max_(3,1) = 1.381776e-2
   both reconcile only with C_(3,1) = 3.175614005. Exact Δ = **7.7193e-3**
   (§4's 7.7e-3 ✓; pinning test's 7.72e-3 ✓ exact).

Both stay data_gap/`kind: conflict` for the 2026 journal; neither touches the
C^stable/C^max/T^stable golden columns, which re-verify as internally
consistent (all Table 4 ΔC^stable/ΔC^max/ΔC identities re-checked exactly:
4.08220e-10, 4.069e-12, 1.2100252e-5, 5.8784844e-4, 3.5996891e-4,
3.6577377e-2 vs the printed values — all agree to printed precision).

Two new micro-findings (no action, recorded for completeness):

- **Table 4 day column is the rounded-TU-column conversion**: e.g.
  10.292069 TU × 4.348377401631 d = 44.753800 d (printed), whereas the
  full-precision T^stable × TU = 44.753801 d; likewise 84.534335 vs
  84.534334 and 64.305944 vs 64.305946. Last-printed-digit only. §4's "match
  Table 4's day column exactly" holds at the 4-dp precision quoted there. The
  (k1,1) rows' `tof_days_bounds` use Table 4's printed (sourced) values; the
  derived (3,2)/(3,3) conversions 77.838478 / 78.903311 d re-verify exactly.
- §4 item 3 (C_(1,1) Table 3 vs Table 1): exact Δ = 4.508e-11, at the 11th
  *decimal place*; strictly not a rounding of the Table 1 value (which would
  round to …315, not …3600) but benign in magnitude — characterization
  unchanged.

**Bottom line: zero NOTE-ERRORs and zero WIRING-ERRORs in any golden or
catalogue value; the only note defect was the 2.24e-6-vs-2.23e-6 prose
rounding above (fixed 2026-06-12). The 5 V1 rows and the golden tests stand
as written.**
