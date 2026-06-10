# Şaloğlu/Taheri/Landau 2023 — iso-impulse multiplicity & ΔV-allocation (mining)

**Date:** 2026-06-10
**Source:** K. Şaloğlu, E. Taheri, D. Landau, "Existence of Infinitely Many
Optimal Iso-Impulse Trajectories in Two-Body Dynamics," *Journal of Guidance,
Control, and Dynamics*, Vol. 46, No. 10, 2023. doi:10.2514/1.G007409. Read from
the AAS 23-307 conference version (33rd AAS/AIAA Space Flight Mechanics Meeting,
Austin TX, Jan 2023; 38 pp.); **page numbers below are AAS-manuscript pages;
equation/table/figure numbers are from that version and may be renumbered in the
JGCD print** — re-verify eq. numbers against the journal copy before citing them
in code comments.
**Writeback: NONE.** Extraction + adoption assessment only; no code or catalogue
edits in this pass.

**Verdict: USEFUL** on three of the five mining lenses, with one important
negative. (1) The ΔV-allocation machinery is fully extractable and trivially
implementable (a quadratic per phasing orbit). (2) The impulse-count bounds are
real but apply to *rendezvous/transfer maneuvers between orbits with loiter
freedom* — they do **not** directly bound our epoch-anchored, flyby-bounded
cycler legs (the honest negative), though they do supply a cheap analytic
feasibility test (Eq. 16) for whether a multi-rev leg admits iso-ΔV impulse
splits. (3) The phasing-orbit primer-vector reinterpretation documents — and
fixes — the exact singular failure mode our `verify/primer.py` two-point
construction has on near-integer-revolution coasts. (4) Three fully published
boundary-value problems (Earth–Dionysus, GTO–GEO, LEO–GEO) with input states
*and* output ΔVs are clean candidate goldens for the multi-impulse lane.

---

## 1. Core construct: impulse anchor positions + the fundamental two-impulse arc

The empirical observation driving the whole paper (Sec. II.A, Fig. 1, pp. 4–5):
in long-time-horizon multi-revolution minimum-ΔV impulsive solutions (the
Earth–Dionysus benchmark of Taheri & Junkins, "How many impulses redux," their
ref [43]), **all intermediate impulses occur at only two distinct spatial
positions** — coined *impulse anchor positions* — and all impulses at a given
anchor point in the **same direction** (the primer direction there). In Fig. 1
(six-impulse, N_rev = 5 solution): five impulses totalling 7.521545 km/s at one
anchor, one impulse of 2.38588 km/s at the Dionysus-orbit intersection.

Every multi-impulse minimum-ΔV solution decomposes into four coast phases
(Conclusions, p. 34): (1) coast on the initial orbit (late departure), (2) a
to-be-determined number of *phasing orbits* through the anchor position, (3) a
**phase-free two-impulse arc** (the fundamental solution), (4) coast on the
target orbit (early arrival).

The fundamental arc is found by the parametric problem (Eq. 1, p. 5):

```
minimize_{θ_i, θ_f, t_pf}   ‖Δv_i‖ + ‖Δv_f‖
```

over departure/arrival true anomalies on the two terminal orbits and the
transfer time, evaluated with a zero-revolution Lambert solver (their refs
[50] Gooding 1990, [51] Prussing 2000). This deliberately **circumvents
switching-surface generation and TPBVP continuation** — they report
orders-of-magnitude speedup over the minimum-fuel/minimum-integral-acceleration
TPBVP routes (Contributions list, p. 3).

## 2. The ΔV-allocation machinery (Sec. II.B–D, pp. 6–11) — implementable

Given the fundamental arc, multi-impulse solutions are *generated* (not
searched for) by splitting the anchor impulse Δv (the 7.521545 km/s one) into
same-direction sub-impulses, each injecting onto a phasing orbit that returns
to the anchor position. All recipes below are algebraic.

- **Timing constraint, one phasing orbit** (Eq. 2):
  `M·T_E + N·T_p + N_D·T_D + N_pf·T_pf = TOF`, where M, N, N_D, N_pf are
  integer revolution counts on the Earth orbit, phasing orbit, Dionysus orbit,
  and the phase-free-arc orbit; TOF is the *surplus* time (total mission time
  minus the three common arcs), not the mission TOF.
- Phasing-orbit period/energy: `T_p = 2π√(a_p³/μ)` (Eq. 3), vis-viva at the
  anchor radius r (Eq. 4).
- **Split parameter:** Δv₁ = α·Δv, α ∈ (0,1), applied along the fundamental
  impulse direction (= primer direction, ‖p‖=1 there). Velocity matching
  (Eq. 5): `‖v_E + Δv₁‖ = v_p`; period as a function of α (Eqs. 6–7):

  ```
  T_p(α) = 2πμ √( −1 / [ 8 (v_p²/2 − μ/r)³ ] ),   v_p² = ‖v_E + αΔv‖²
  ```

- **The allocation quadratic** (Eqs. 8–9): given a required T_p,

  ```
  (ΔvᵀΔv)·α² + 2(v_EᵀΔv)·α + v_Eᵀv_E + 2(μ²π²/(2T_p²))^{1/3} − 2μ/r = 0
  ```

  Sum of roots is negative (a > 0, b > 0 here), so **the positive root is the
  unique feasible α** (p. 8–9). Workflow: pick integers (M, N, …), solve Eq. 2
  for T_p, reject if T_p < period of the initial orbit or negative, else solve
  Eq. 9 for α and accept if α ∈ (0,1).
- **Four impulses / two phasing orbits** (Eqs. 10–12): one constraint, two
  unknowns (α₁, α₂) ⇒ a one-parameter **infinite family**; discretize α₁, get
  `T₂(α₂) = (TOF − M·T_E − N₁T₁(α₁) − N_D T_D − N_pf T_pf)/N₂` (Eq. 12), solve
  Eq. 9 for α₂. Ordering constraint (Eq. 11): `0 < α₁ < α₂ < 1`.
- **General n_p phasing orbits** (Eqs. 13–15):
  `Σ_{k=1}^{n_p} N_k T_k + M T_E + N_D T_D + N_pf T_pf = TOF` (Eq. 13);
  `0 < α₁ < … < α_{n_p} < 1, T_i > 0` (Eq. 14); last α from the quadratic
  (Eq. 15).
- **Existence test for a given (M, N_k, N_D, N_pf) combination** (Eq. 16):

  ```
  T(α=0) ≤ (TOF − N_D T_D − N_pf T_pf − M T_E) / Σ N_k ≤ T(α=1)
  ```

  where T(α=0) is the initial-orbit period and T(α=1) the phase-free-arc
  orbit period. Violated ⇒ no feasible family for that integer combination.
- **Direction restriction** (Remark, p. 20): phasing-orbit periods must lie in
  (T(α=0), T(α=1)); a retrograde split (shortening the period below the
  initial orbit's) always costs extra ΔV. The split never changes total ΔV —
  **iso-ΔV by construction** — only mission time and impulse magnitudes.

### Solution-space envelopes (Sec. III.D, Eqs. 17–23, pp. 24–27)

Feasible (T₁, …, T_{n_p}) regions are polygons whose corners are analytic:
left-most corners at T₁,min = T_E; `T_{n_p},max` capped at T(α=1);
bifurcation point where all periods equal at `T = TOF/n_p` (general corner
formulas Eq. 21, cap corner Eq. 22, adjacent-corner Eq. 23). Five-impulse
Earth–Dionysus example: T₂,max = T₃,min = (2490.48 − 365.25)/2 = 1062.62 d
(Eq. 17), T₃,max = 1161.47 d (Eq. 18), T₂,min = 963.77 d (Eq. 19), corner-4 at
664.51 d (Eq. 20), corner-5 at TOF/3 = 830.16 d. Useful as closed-form
*bounds on phasing-orbit periods* if we ever enumerate these families.

## 3. Impulse-count bounds (Sec. III.E, pp. 27–30) — lens 2, the careful part

The paper provides analytic lower (required) and upper (allowable) impulse
counts for three maneuver classes in inverse-square gravity:

1. **Fixed-terminal-time rendezvous** (time-dependent boundary states):
   lower bound 3 for Earth–Dionysus (no phase-*constrained* two-impulse
   solution exists at the optimal ΔV = 9.907425 km/s; p. 12), upper bound
   from the lower limb of Eq. 16: `n_p < TOF/T(α=0)` ⇒ for Earth–Dionysus
   2490.48/365.25 = 6.82 ⇒ n_p ≤ 6 ⇒ **3 ≤ impulses ≤ 8** (p. 29). Applicable
   to "any large, but finite-time-horizon maneuver."
2. **Free-terminal-time rendezvous:** TOF can be inflated by integer multiples
   of the target period (β·T_D), so the upper count is unbounded; for any
   *given* finite TOF the maximum is always computable from Eq. 16.
3. **Phase-free transfer** (match elements except true anomaly): lower bound
   **2** (the fundamental arc); upper bound unbounded with unbounded time
   (add κ·T_E, κ ∈ ℕ; p. 29). Total ΔV identical across all classes' families.

### Does this bound our per-leg DSM count? Mostly NO — record the negative.

- The bounds presuppose the **anchor structure**: terminal coast freedom
  (late-departure/early-arrival on the *bodies' orbits*) plus phasing orbits
  re-encountering a fixed spatial anchor. Our cycler legs are
  **epoch-anchored, flyby-bounded point-to-point arcs**: the spacecraft cannot
  loiter on Earth's orbit (it arrives/departs hyperbolically with a v∞ fixed
  by the adjacent legs), and the terminal impulses of the fundamental arc are
  absorbed by flybys, not burned. The paper's maneuver classes (rendezvous /
  orbit-to-orbit transfer) simply do not contain our boundary conditions, so
  **the 3..8-style bounds cannot replace the one-DSM-per-leg fiat**.
- What *does* transfer: the structural theorem that extra impulses beyond the
  fundamental solution are **iso-ΔV redistributions, never ΔV improvements**,
  *whenever* the optimum has the anchor structure (surplus TOF spent coasting).
  For a short (< 1 rev) cycler leg, Eq. 16 fails for every n_p ≥ 1 (no phasing
  orbit fits), i.e. the iso-ΔV split mechanism is unavailable and a single
  interior impulse is not leaving anchor-type savings on the table. For
  **long multi-rev legs — exactly our S1L1-type resonant intervals** — Eq. 16
  is a one-line feasibility test for whether the leg admits same-direction
  impulse splits at the encounter position without ΔV penalty. Note the S/L
  Earth-to-Earth resonant intervals of the S1L1 nomenclature *are* this
  paper's phasing orbits (N revs returning to the Earth-encounter anchor);
  the difference is our "impulses" there are flybys with bounded turn, not
  free same-direction burns, so the allocation freedom is narrower.
- Caveat the other way: for **fixed-time** legs whose optimum does *not* have
  the anchor structure, the true minimum-ΔV impulse count is governed by
  Taheri & Junkins (their ref [43]) switching-surface analysis, not by this
  paper. The one-DSM fiat remains a transcription choice justified (or not)
  by primer diagnostics per leg — which is what `verify/primer.py` +
  `diagnose_impulse_schedule` already do.
- **Seeding (standing problem b): USEFUL.** The generation procedure is a
  deterministic enumeration of multi-impulse basins from one cheap two-impulse
  solution: integer grid over (M, N_k, N_D, N_pf) filtered by Eq. 16, then
  α from a quadratic. No MBH needed to *find* these families; MBH only needed
  when the anchor assumptions break. This is a credible seeding recipe for
  multi-rev legs before handing to the corrector.

## 4. Primer vector on phasing orbits (Sec. III.F, pp. 31–34) — lens 3

- Standard machinery restated: p = −λ_v (minimum-principle sign), ‖p(t_i)‖ = 1
  with p along Δv at impulses (Eq. 24); ṗ = λ_r (Eq. 25); interior impulses
  need ṗ(t_i) = ‖p‖-rate zero, `ṗᵀp = 0` (Eq. 26); coast-arc mapping by the
  6×6 STM partitioned [Φ₁ Φ₂; Φ₃ Φ₄] (Eq. 27); initial primer rate recovered
  from endpoint primer values via Glandorf's inversion (Eq. 28):
  `λ_r(t0) = Φ₂⁻¹Φ₁ λ_v(t0) − Φ₂⁻¹ λ_v(t_f)` (their ref [54], Glandorf 1969).
  This Eq. 28 construction is **exactly our `primer_on_coast`** solve
  `ṗ(0) = Φ_rv⁻¹ (p₁ − Φ_rr p₀)`.
- **The documented failure mode (p. 32):** over one full revolution of a
  phasing orbit the STM is the identity ⇒ Φ₂ singular ⇒ Eq. 28 has no
  solution. Our `primer_on_coast` will hit the same `np.linalg.solve`
  singularity (or severe ill-conditioning) on near-integer-rev coasts — i.e.
  precisely the multi-rev resonant legs we care about. **Their fix:** use
  continuity of p and ṗ at the phasing-orbit entry/exit (where p is known
  from the impulse directions) and propagate the primer *through* the phasing
  orbit with Eq. 27 directly, instead of re-solving the two-point problem on
  the singular arc.
- **New interpretation — "potential impulse opportunities":** on a phasing
  orbit with k revolutions, ‖p(t)‖ touches 1 (with ṗ = 0) once per revolution
  at the anchor passage (Fig. 22a/b: 3-rev case, three unit-touches inside the
  shaded phasing interval; TU = 5022750.126364 s for that plot, solution of
  Fig. 6a). All touches satisfy Lawden's necessary conditions, but an impulse
  is applied only at entry/exit of the phasing sequence. Geometric picture
  (Fig. 23): at each opportunity p sits on the unit sphere and ṗ lies in the
  tangent plane; all impulses at one anchor share a single supporting plane.
- **Consequence for our diagnostics:** in `diagnose_impulse_schedule`,
  `max‖p‖ = 1` attained in the *interior* of a multi-rev coast must NOT be
  read as "IMPROVABLE — add impulse"; it can be a degenerate iso-ΔV impulse
  opportunity fully consistent with optimality (adding the impulse with the
  primer direction changes nothing; adding it anywhere else costs ΔV). Our
  current `tol`-guarded `> 1 + tol` verdict is technically safe (touching 1 is
  not flagged), but the caveat string and the near-integer-rev Φ_rv
  conditioning deserve an explicit guard + this citation.

## 5. Published numbers — candidate goldens (lens 4)

All three problems are two-body, fully specified by Cartesian boundary states +
TOF, with published ΔV decompositions. These satisfy the sourced-golden rule
(published values, not values our code computed).

### 5.1 Earth–Dionysus benchmark (heliocentric; Table 1, p. 11; Sec. III.A pp. 12–13)

| Quantity | Value |
|---|---|
| Initial position (km) | [−3637871.082, 147099798.784, −2261.441] |
| Initial velocity (km/s) | [−30.265, −0.849, 5.053e−5] |
| Final position (km) | [−302452014.884, 316097179.632, 82872290.075] |
| Final velocity (km/s) | [−4.533, −13.110, 0.656] |
| Time of flight | 3534 days |

Fundamental phase-free two-impulse solution (Sec. III.A, p. 12): t_pf = 348.46 d,
Δv₁ = 7.521545 km/s (anchor impulse, at Earth-orbit position), Δv₂ = 2.38588
km/s (Dionysus-orbit intersection), **total ΔV = 9.907425 km/s** (matches the
four extremals of their ref [43] to seven digits); impulse true anomalies
θ_i = 179.27°, θ_f = 149.20°. Common coasts: Earth orbit 193.24 d, Dionysus
orbit 501.81 d ⇒ surplus TOF = 3534 − (193.24 + 348.46 + 501.81 = 1043.51) =
**2490.48 d**. Periods: T_E = 365.25 d, T_D = 1191.88 d, T_pf = 1161.47 d.
A phase-*constrained* two-impulse Lambert between the same states always
exceeds 9.907425 km/s (p. 12) — itself a testable inequality.

Spot-check solutions (Sec. III.A, Figs. 6–11):
- 3-impulse, M=0, N₁=3: phasing period 830.16 d (Fig. 6a).
- 3-impulse, M=5, N₁=1: phasing period 664.23 d (Fig. 6b).
- 3-impulse, M=0, N₁=3, N_pf=1: TOF→1329.01 d, phasing period 443.00 d (Fig. 6d).
- 26 three-impulse families total (6/9/5/6 with 4/5/6/7 revolutions, p. 15);
  Table 2 (p. 19) family census: 3→26, 4→42, 5→37, 6→21, 7→7, 8→1 (134
  families, one anchor); two-anchor analysis collapses to **5 solutions**
  (pp. 18–19).
- Shorter-mission variant (pp. 19–20): subtract T_D ⇒ 2342.12 d; also drop the
  final coast ⇒ 1840.31 d. **No minimum-ΔV solution exists with the last
  impulse inside a 1840.31-d window** — the time-feasibility statement
  (matches Table 10, N_rev = 4 of their ref [43]). Four-impulse short
  solution: N₁=1, N₂=2 (Fig. 12), with T₂ vs T₁ range in Fig. 21.

### 5.2 GTO→GEO (geocentric; Table 3, p. 21)

| Quantity | Value |
|---|---|
| Initial position (km) | [6721.95652173912, 0, 0] |
| Initial velocity (km/s) | [0, 10.0384360619658, 1.23256496402036] |
| Final position (km) | [−42165, 0, 0] |
| Final velocity (km/s) | [0, −3.07462812005026, 0] |
| Time of flight | 2 d (Case 1), 12.011 d (Case 2) |

Tangential orbits ⇒ single-impulse base solution **ΔV = 1.4873346 km/s** at GTO
apogee (p. 20). Case 1: three impulses (N₁=N₂=1, Fig. 13; 1 DU = 6378 km).
Case 2: 17 impulses (16 phasing orbits, one rev each, Fig. 14); same total ΔV;
first two impulses ≈ 0.2 km/s, magnitudes decaying along the sequence (Fig. 15
— note 17-impulse bars plotted ×5). Demonstrates the operational point: stretch
TOF until max individual impulse falls below the propulsion limit (p. 22).

### 5.3 LEO→GEO (geocentric, non-intersecting, 28° plane change; Table 4, p. 23)

| Quantity | Value |
|---|---|
| Initial position (km) | [−6677.99994822088, −0.734261017649844, −0.390413508385778] |
| Initial velocity (km/s) | [0.000962087495564221, −6.82150753761816, −3.62705989590144] |
| Final position (km) | [42163.9999975184, −0.457455672469879, 0] |
| Final velocity (km/s) | [3.33583804464932e−5, 3.07466457999987, 0] |
| Time of flight | 1 d |

Two-impulse base: phase-free arc 0.220 d, **total ΔV = 4.2206849 km/s**, no
terminal coasts (boundary anomalies coincide with the impulse anomalies), so
phasing surplus = 0.780 d (p. 23). The minimum-ΔV split distributes the plane
change across both impulses (p. 23, citing Curtis Chp. 6). Four-impulse
solution (N₁=N₂=1, Fig. 17): largest impulse 1.944 km/s (first phasing
injection, p. 24).

**Golden assessment:** Earth–Dionysus is the best catalogue-grade golden —
heliocentric, full Cartesian boundary states, published total + per-impulse ΔV,
plus *derived* checkable quantities (t_pf, coast durations, T_pf, the 2490.48-d
surplus, the 830.16/664.23/443.00-d phasing periods, the n_p ≤ 6 bound). The
two geocentric cases are unit-level goldens for the allocation machinery (just
change μ); they exercise plane change (LEO–GEO) and the tangential
single-impulse degenerate base (GTO–GEO). Caveats: μ values are not printed in
the paper (standard GM_Sun/GM_Earth assumed); reproducing the fundamental-arc
numbers requires our own Eq.-1 parametric optimizer (the published t_pf/ΔV/θ
values are the targets, ~5–7 significant digits); figure-read values
(e.g. 0.2 km/s first impulses) are approximate — anchor tests on the tabled
and texted numbers only.

## 6. Honest applicability limits (lens 5)

- **Dynamics:** inverse-square two-body only. (CR3BP extension of the
  switching-surface side is the authors' AAS 22-838, their ref [49] — separate
  source.) No perturbations, no mass/Isp modeling (pure ΔV).
- **Geometry:** *not* restricted to coplanar or circular terminals — arbitrary
  3D orbits (Earth–Dionysus is inclined/eccentric; LEO–GEO has 28° plane
  change). Genuinely general within two-body.
- **Boundary conditions:** rendezvous (fixed/free time) and phase-free
  orbit-to-orbit transfer. **Not** epoch-fixed point-to-point with hyperbolic
  flyby boundaries — our cycler-leg class is outside all three maneuver
  classes (Sec. III.E list, p. 29). Terminal-coast freedom (late
  departure/early arrival) is load-bearing throughout.
- **Iso-ΔV only:** the machinery never lowers ΔV below the fundamental
  two-impulse (or single-impulse) value; it redistributes impulses and trades
  mission time. It is a multiplicity/realizability tool (impulse-magnitude
  caps, time-feasibility), not an optimizer.
- **Anchor-position premise is empirical-then-verified:** anchors are
  discovered from the fundamental arc's impulse locations and verified by PVT;
  the paper notes these solutions are *not* directly generable from the
  optimal-switching-surface construct (Contributions item 3, p. 3).
- **NOT-USEFUL for us:** (a) the impulse-count bounds as a per-leg DSM-count
  replacement (wrong boundary-condition class — see §3); (b) the two-anchor
  enumeration detail (pp. 18–19) and the envelope-corner taxonomy (Eqs. 17–23)
  beyond period bounds — catalogued here for completeness, no current
  consumer; (c) nothing here helps the FD-Jacobian bottleneck or basin
  selection for *non*-anchor-structured legs.

## 7. Adoption assessment for cyclerfinder

1. **Primer diagnostics hardening (cheapest, highest value):** add to
   `verify/primer.py` (a) a conditioning guard on the Φ_rv solve in
   `primer_on_coast` for near-integer-rev coasts, with the continuity-based
   fallback (propagate p, ṗ from the adjacent impulse via Eq. 27 instead of
   the two-point Eq. 28 inversion); (b) a "potential impulse opportunity"
   clause in the multi-rev caveat: interior ‖p‖ → 1 touches with ṗ ≈ 0 on
   near-resonant coasts are iso-ΔV degeneracies, not improvability signals
   (Sec. III.F).
2. **Iso-ΔV split feasibility probe:** Eq. 16 as a one-line test on long
   multi-rev legs (S/L resonant intervals) for whether the leg's anchor
   impulse can be split across phasing revolutions without ΔV change —
   relevant to impulse-magnitude-capped variants and to explaining iso-ΔV
   plateaus the MBH lane finds.
3. **Deterministic multi-impulse seeding:** for legs that *do* have anchor
   structure (encounter-anchored multi-rev phasing), enumerate
   (M, N_k, N_D, N_pf) under Eq. 16 + the α-quadratic (Eq. 9/15) to seed the
   corrector in known-distinct basins, instead of relying on MBH stochastic
   hops. Directly addresses standing problem (b) for that leg class.
4. **Goldens:** stage the three §5 problems as multi-impulse-lane validation
   targets (Earth–Dionysus fundamental arc + one 3-impulse family check;
   LEO–GEO two-impulse base; GTO–GEO single-impulse base + 3-impulse split).
   Per the backfill checklist: center = Sun (Dionysus case) / Earth (others),
   tof_days_bounds from the tables, source_ephemeris = n/a (Cartesian states
   published directly — no ephemeris dependency, which is what makes these
   unusually clean).
5. **Keep the fiat, justified:** one-DSM-per-leg stands for short legs (Eq. 16
   infeasible ⇒ no iso-ΔV splits exist; any genuine improvement claim must
   come from primer diagnostics, unchanged).

## 8. Reference leads worth pulling later

- [43] Taheri & Junkins, "How many impulses redux," *J. Astronautical
  Sciences* 67(2), 2020, 257–334 — the switching-surface machinery and the
  Earth–Dionysus extremal tables (Table 10 cited for the N_rev = 4 case);
  the proper source for fixed-time impulse-count questions on
  non-anchor-structured problems.
- [49] Saloglu & Taheri, "Acceleration-based switching surfaces for impulsive
  trajectory design in restricted three-body dynamics," AAS 22-838, 2022 —
  CR3BP extension; relevant to the CR3BP moon-tour lane.
- [41] Landau, "Efficient maneuver placement for automated trajectory design,"
  JGCD 41(7), 2018 — PVT-driven maneuver placement for many-rev problems
  (same third author); candidate for the DSM-placement seeding story.
- [51] Prussing, "A class of optimal two-impulse rendezvous using
  multiple-revolution Lambert solutions," *J. Astronautical Sciences* 48(2),
  2000 — multi-rev Lambert branch handling for the fundamental-arc solver.
- [54] Glandorf, "Lagrange multipliers and the state transition matrix for
  coasting arcs," *AIAA Journal* 7(2), 1969 — the Eq.-28 primer-rate
  inversion our `primer_on_coast` implements; the citation to attach to the
  conditioning guard.
