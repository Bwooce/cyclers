# Mining Note: Şaloğlu & Taheri (2025) — Iso-Impulse 3D Classification & Feasibility

**Source:** K. Şaloğlu, E. Taheri, "Classification and Feasibility Assessment of Infinitely
Many Iso-Impulse Three-Dimensional Trajectories," *The Journal of the Astronautical
Sciences* (2025), DOI: 10.1007/s40295-025-00528-0; preprint arXiv:2501.01583 (3 Jan 2025,
40 pp. — page/equation numbers below refer to the arXiv version).

**Mined:** 2026-06-10. This is the 3D follow-up to Saloglu, Taheri & Landau, JGCD Vol. 46
No. 10 (2023) [their Ref. 52], which is mined separately; this note covers only what is
NEW here.

**TRANSCRIPTION RESCAN (2026-06-12):** all load-bearing values re-read against the
arXiv PDF (two passes per load-bearing table) — **MATCH throughout**, and the two
flagged source defects are confirmed in print. Verified: Table 1 p. 18 (initial
[7000, 0.02, 60°, 0°, 0°, θ_i]; target [105000, 0.3, 12°, 0°, 0°, θ_f]); two-impulse
base p. 18 verbatim ("Δv_total is **3.9618011** km/s with of Δv₁ = **2.8246140**
km/s and Δv₂ = **1.1371871** km/s"); three-impulse base p. 19 (**3.8641159** =
2.9390 + 0.6815 + 0.2436 km/s, particleswarm); coast 6.4738 d, T_f = 3.9191 d,
0.73 ≤ Σ N_{k,1} ≤ 58.09 (p. 20). Earth–Mars Sec. VI.A: 5.5865 / 5.5873 km/s,
313.2430 / 801.4034 d, Lambert 6.047 km/s, 793 + 686.9658 = 1479.9658 d /
759.6963 d — all as noted. **SOURCE-DEFECT-CONFIRMED (both):** p. 30 prints
"becomes 701.2694 days and 1405.8939 days" in prose and then
"793 − 720.2694 = 72.7306 < T₀" in the feasibility check on the same page
(arithmetic supports 720.2694); Δv_total printed **5.6109** km/s on p. 30 (Fig. 21c
sentence) vs **5.6108** km/s on p. 31 — genuinely inconsistent in print, exactly as
recorded in §7. The #207 goldens wired from this paper in
`tests/verify/test_dv_bracket.py` (3.9618011 and its breakdown) are
transcription-clean.

**Verdict: USEFUL (secondary, two narrow takeaways).** The headline machinery — generating
infinitely many equal-ΔV solutions by inserting phasing orbits at impulse anchor positions
(APs) — is structurally inapplicable to our epoch-anchored, time-fixed cycler legs (their
own Earth–Mars example demonstrates the failure mode explicitly, p. 30). What survives for
us: (1) the phase-free two-/three-impulse **base solution as a rigorous per-leg ΔV lower
bound** ("certificate"), cheap to compute and bracketing our MBH/DSM results from below
with Lambert from above; (2) the **three-impulse base-solution parameterization (Eq. 7)**
as a candidate seed generator for the queued broken-plane genome (#3). One self-contained
geocentric golden extracted (Table 1 + ΔV breakdowns); the Earth–Mars numbers are
catalogued but depend on Ref. 49 for the boundary orbits.

---

## 1. What the 3D extension adds over the 2023 paper

The 2023 JGCD paper [52]: two-impulse, time-free, phase-free base solutions; ΔV-allocation
at a **single** AP (the one with the largest impulse); planar examples. This paper adds
(contributions list, pp. 3–4):

1. **Three-impulse base solutions** (Sec. III). The base optimization becomes

   > minimize over (θ_i, θ_f, t_pf,1, t_pf,2, r_mid, φ, λ):
   > ΔV_total = Σ_{i=1..3} ||Δv_i||   (Eq. 7, p. 8)

   where θ_i, θ_f are true anomalies of the (phase-free) departure/arrival points on the
   terminal orbits, t_pf,1/t_pf,2 are the two connecting-arc coast times, and the
   **midcourse impulse position is parameterized in spherical coordinates** (r_mid,
   azimuth φ, elevation λ). Solved globally with MATLAB `particleswarm` (p. 19).
   Available TOF generalizes Eq. 2 to Eq. 8 (p. 8):
   TOF = (t_f − t_0) − (t_c1 + t_pf,1 + t_pf,2 + t_c2).

2. **Simultaneous ΔV-allocation at 2–3 APs** (Sec. IV): the single-AP time constraint
   (Eq. 3, p. 6) generalizes to

   > Σ_{i=1}^{n_i} Σ_{k=1}^{n_p,i} N_{k,i} T_{k,i}(α_{k,i}) + N_0 T_0 + N_f T_f
   > + N_pf,1 T_pf,1 + N_pf,2 T_pf,2 = TOF   (Eq. 10, p. 15)

   with per-AP impulse-ratio orderings 0 < α_{1,i} < … < α_{n_p,i} < 1 (Eq. 11). The
   allocation stays fully analytic: discretize α's at all but one slot, get periods from
   vis-viva (Eq. 5, p. 7: T_p = 2πμ·√(−1/(8((v_p²/2) − μ/r)³)) with
   v_p = ||v_0 + α_k Δv||), close the last α by the quadratic
   v_0ᵀv_0 + 2α v_0ᵀΔv + α²ΔvᵀΔv + 2(μ²π²/(2T²))^{1/3} − 2μ/r = 0 (Eq. 6, p. 8).
   Impulse **direction** is never changed (fixed by primer vector theory); only magnitude
   is split (Sec. II.B, p. 7).

3. **Selection strategy between two- and three-impulse base solutions + ΔV-optimality
   certificate** (Sec. III.A, Figs. 4–5; Sec. VI) — see §2 below.

4. **Four-layer classification** of all iso-impulse solutions (Fig. 6, p. 14) — see §4.

### Time-feasibility relations (the implementable core, Sec. IV.B)

Time available for intermediate phasing orbits only:

> TOF_p = TOF − N_f T_f − N_pf,1 T_pf,1 − N_pf,2 T_pf,2 − N_0 T_0   (Eq. 12, p. 16)

Two-sided feasibility with 2–3 APs (Eq. 13, p. 16):

> T_0 ≤ [TOF_p − Σ_{i=2,3} Σ_k N_{k,i} T_{k,i}(α_{k,i}=0)] / Σ_k N_{k,1} ≤ T_{k,1}(α_{k,1}=1)

Eqs. 14–15 (p. 16) extract upper/lower bounds on the total revolution counts Σ N_{k,i} per
AP (setting the other APs' counts to 1). Caveat (pp. 16–17): if the AP's total impulse
*decreases* orbital energy, T_{k,i}(α=1) is the *minimum* period and the bound flips.
Solution-envelope corner points are closed-form (Appendix Eqs. A1–A9, pp. 34–36; also
verifiable as a small LP, p. 18 — design vector = all phasing periods, equality = Eq. 10,
inequalities T_{k,i} ≤ T_{k+1,i}).

---

## 2. Base-solution selection and the ΔV-optimality certificate

**Circle-to-circle (i, β) map (Fig. 4, p. 10; β = r0/rf):** a reproduction (with added
data points) of Marec's impulse-count partition for coaxial inclined circular orbits.
Calibration points: A(i=0°, β=0.084), B(37.389°, 0.149), C(40.841°, 0.1719),
D(40.841°, 0.2473), E(60.185°, 1.0), F(0°, 1.0); planar Hohmann/bi-elliptic separatrix
11.94 ≤ rf/r0 ≤ 15.58 (p. 9, citing Curtis Fig. 6.8). Practical reading (bullets,
p. 10): three-impulse base solutions win only for large inclination change
(i ≳ 40°) or extreme radius ratios; for β > 0.2473 and i < 40.841° either can win and
TOF breaks ties.

**Selection flowcharts (Fig. 5, p. 13):**
- *Time-free* (5a): pick min ΔV_total first; tie-break on shorter TOF.
- *Time-fixed rendezvous* (5b): time-feasibility filter **first** (TOF ≥ T_0, i.e. at
  least one phasing orbit of the smallest period must fit, p. 11), then ΔV-optimality.
  Outcomes: feasible base → theoretical minimum-ΔV is recoverable; "relax the time
  constraint" (extend mission by κT_f, Eq. 9, p. 11 — adding whole target-orbit periods
  preserves phasing) → recoverable at the cost of mission time; "shorter mission time" →
  **no base solution**, fall back to Lambert/NLP, and the base ΔV_total is then a strict
  lower bound on what any time-fixed solution can achieve (pp. 12, 31–32).

This is the **certificate**: for a general 3D time-fixed rendezvous, base-solution
ΔV_total (min of Eq. 1 and Eq. 7 solutions) lower-bounds the achievable ΔV; the two-impulse
Lambert solution at the fixed TOF upper-bounds it (p. 32). The certificate tells you
*whether* the theoretical minimum is attainable at your mission time, before running any
time-fixed optimizer.

---

## 3. ΔV-allocation under per-impulse caps (Δv ≤ Δv_max)

Mechanism (Secs. II.B, IV.A): split the impulse at an AP into n_p smaller collinear
impulses, each injecting into a phasing orbit; the spacecraft spends integer revolutions
N_{k,i} on each. ΔV_total is exactly preserved; the price is **always paid in flight
time** (multiples of phasing-orbit periods). With 2–3 APs all large impulses can be capped
simultaneously (p. 4: "one can reduce all individual impulse magnitudes ... such that
individual impulses satisfy the Δv ≤ Δv_max constraint"; conclusion p. 33).

**Relevance to our maintenance/DSM impulse splitting: NOT applicable in-leg.** A cycler
leg is phase-fixed at both ends with TOF set by the encounter epochs; inserting phasing
revolutions breaks the encounter schedule. The only place the mechanism fits is cycler
*establishment/insertion* (spreading a large injection burn over phasing orbits before the
first Earth departure), which is standard practice we already get without this formalism.

---

## 4. Four-layer classification — enumeration tool or theory?

Layers (Fig. 6, p. 14): 1) base solution (two- or three-impulse); 2) feasible solution
space = integer points (Σ N_{k,1}, Σ N_{k,2}[, Σ N_{k,3}]) inside a polygon/polytope from
Eq. 13 (Figs. 11, 14, 17); 3) solution families = partitions of each revolution total
across distinct phasing orbits; 4) solution envelopes = polygonal period ranges per
phasing orbit (Fig. 7; Appendix).

**Assessment for MBH seeding: theory only, for us.** The enumeration enumerates *phasing
orbit* structures at APs — degrees of freedom our legs do not possess. Our family indexing
(Lambert arc revolution counts, DSM placement) is disjoint from this hierarchy. The
polytope counting is elegant (1817 feasible solutions for the geocentric 2-AP case, p. 24;
4747 for 3-AP, p. 26) but counts things we cannot use. Negative recorded.

---

## 5. Published inputs/outputs (candidate goldens)

All dynamics are pure two-body. DU = Earth equatorial radius for geocentric cases (p. 8).

### 5.1 Geocentric three-impulse problem (Sec. V.A — best self-contained golden)

Inputs (Table 1, p. 18), elements [a (km), e, i, ω, Ω, θ]:
- Initial: [7000, 0.02, 60°, 0°, 0°, θ_i free]
- Target: [105000, 0.3, 12°, 0°, 0°, θ_f free]  (co-axial, 48° inclination change)

Outputs:
- Two-impulse base (p. 18): ΔV_total = **3.9618011 km/s** (Δv1 = 2.8246140,
  Δv2 = 1.1371871 km/s); departs initial-orbit perigee, arrives target apogee; plane
  change concentrated at the second (apogee) impulse. Primer magnitude ≤ 1 throughout
  (Fig. 8b) — extremal.
- Three-impulse base (Eq. 7 via particleswarm, p. 19): ΔV_total = **3.8641159 km/s**
  (Δv1 = 2.9390, Δv2 = 0.6815, Δv3 = 0.2436 km/s); coast total 6.4738 days
  (p. 20); primer touches 1 at all three impulses (Fig. 9b). Three-impulse beats
  two-impulse by ~98 m/s, consistent with Fig. 4 at large Δi.
- Allocation feasibility with κ = 1 added target period (T_f = 3.9191 days): TOF_p =
  3.9191 d, bounds 0.73 ≤ Σ N_{k,1} ≤ 58.09 (p. 20).
- 2-AP case (p. 24): T_0 = 0.0675 d, T_{k}(α=0) = 5.3616 d; 6.3096 ≤ Σ N_{k,2} ≤ 7.2970;
  Σ N_{k,1} ≤ 501.4617; 1817 feasible integer solutions (Fig. 14).
- 3-AP case (p. 26): 5.5786 ≤ Σ N_{k,1} ≤ 443.3749, 5.5786 ≤ Σ N_{k,2} ≤ 6.5660,
  7.2639 ≤ Σ N_{k,3} ≤ 8.6147; 4747 feasible solutions (Fig. 17).
- Envelope corners (p. 28–30, Fig. 20): T_{1,1,max} = 3.1978 d; T_{2,1} ∈ [0.0675,
  4.2413] d, T_{3,1}, T_{4,1} ∈ [0.0675, 5.3616] d; second-AP periods bounded by
  T_{k,2}(α=1) = 7.5606 d; T_{1,3} ∈ [3.9191, 7.0492] d.

Golden quality: inputs fully published, outputs published to 7 figures, model unambiguous
(2BP). To use as an independent golden we must implement the Eq. 1/Eq. 7 phase-free
optimizers (global search; they used particleswarm) — EXPECTED side traces to the
publication, satisfying the sourced-only rule. Note μ_Earth value is not printed; treat
the last 2–3 digits as solver-dependent and validate at ~1e-4 km/s level.

### 5.2 Earth-to-Mars (Sec. VI.A, pp. 30–32, Fig. 21) — partially self-contained

Problem taken from Taheri & Junkins, "How Many Impulses Redux," JAS 67(2), 2020 [their
Ref. 49] — **boundary orbit elements are not reprinted here**; a golden needs Ref. 49's
data (coplanar-looking in Fig. 21a/b, but Fig. 21d shows a 3D trajectory with Z ~0.06 AU).
Published numbers:
- Two-impulse base: ΔV_total = **5.5865 km/s**, flight time 313.2430 days.
- Three-impulse base: ΔV_total = **5.5873 km/s**, flight time 801.4034 days (two-impulse
  wins on both criteria).
- Time-fixed rendezvous, mission time 793 days: with required terminal coasts the
  two-impulse base needs a total time printed as 701.2694 days in one sentence and used
  as 720.2694 in the feasibility check "793 − 720.2694 = 72.7306 < T_0 = 365.25" (p. 30)
  — **internal inconsistency in the paper** (the arithmetic supports 720.2694); either
  way TOF ≥ T_0 fails, so neither base solution is time-feasible at 793 d.
- Time-fixed solution (hybrid PVT-based NLP, initial guess from thrust-acceleration
  continuation): ΔV_total = **5.6109 km/s** (printed 5.6108 on p. 31 — 0.1 m/s
  inconsistency).
- Two-impulse Lambert at 793 d: **6.047 km/s** (from Ref. 49). Bracket demonstrated:
  5.5865 ≤ 5.6109 ≤ 6.047 km/s.
- Time-relaxation example: adding one Mars period, 793 + 686.9658 = 1479.9658 d mission
  time gives available TOF 759.6963 d, recovering minimum ΔV (p. 31).

### 5.3 Circle-to-circle, 45° inclination change, β = 0.5 (Sec. VI.B, p. 32)

- Two-impulse: 1.7036 km/s over 1.3905 days; three-impulse: 1.6853 km/s over 4.2767 days
  (three-impulse ΔV-optimal, as Fig. 4 predicts at i = 45°).
- Time-fixed thresholds (initial orbit period 1.5137 d): Δt > 5.7904 d → three-impulse
  base feasible (min ΔV recoverable); Δt ∈ (2.9042, 5.7904) d → only two-impulse
  feasible, ΔV-optimality sacrificed; Δt < 2.9042 d → no base solution, solve Lambert.
- **Not golden-grade:** absolute radii are not printed (only β = 0.5 and km/s outputs;
  the implied DU is recoverable only by reverse-engineering, which violates the
  sourced-only rule).
- Remark (p. 33): mirror-image solutions about the y–z plane are iso-ΔV; pick the one
  with the lowest terminal coast times for time-fixed use.

### 5.4 Earth-to-Dionysus (two-AP allocation, Sec. V.B.1, pp. 21–23)

Base solution data are from the 2023 paper [52] (θ_i = 179.27°, θ_f = 149.20°,
t_c1 = 193.24 d, t_pf = 348.46 d, t_c2 = 501.81 d, ΔV_total = 9.907425 km/s, mission time
3534 d, TOF = 2490.48 d — restated at Eq. 2, p. 6; goldens belong to the 2023 mining
note). New here (2-AP allocation): with T_0 = 365.25 d and T(α_{k,1}=1) = 1161.47 d,
Σ N_{k,2} ≤ 1.83, Σ N_{k,1} ≤ 3.64 (at Σ N_{k,2} = 1), lower bounds ≥ 1.14 each; only
(2,1) and (3,1) support infinitely many solutions (Fig. 11), versus 134 families for the
single-AP version in [52] (p. 22). Envelope corners (Fig. 19, p. 28): at T_{1,1} = T_0,
T_{2,1} ∈ [365.25, 481.88] d, T_{3,1} ∈ [466.67, 568.11] d; T_{1,1,max} = 443 d;
T_{1,2} ∈ [1161.47, 1191.88] d with an extra corner at T_{1,1} = 432.87 d.

---

## 6. Honest applicability limits for our lane + #3 verdict

Our legs are epoch-anchored at **both** ends (planet states from real ephemerides at
fixed encounter epochs), TOF fixed by the cycler period, phases fixed. Against that:

1. **Iso-impulse generation (the paper's headline): NOT applicable to leg design.** The
   entire infinite-family construction monetizes slack TOF into phasing revolutions.
   Cycler legs have no such slack — and the paper's own Earth–Mars 793-day case lands in
   exactly the "no base solution / solve Lambert" branch (p. 30). Expect every one of our
   legs to fail the TOF ≥ T_0 certificate by construction.
2. **The certificate as a lower bound: APPLICABLE and cheap.** Solving Eq. 1 and Eq. 7
   phase-free between the actual inclined departure/arrival orbits (elements frozen at
   leg epochs) gives a rigorous per-leg ΔV floor independent of our DSM machinery. Useful
   as a sanity band for MBH outputs: base ΔV ≤ our DSM solution ≤ ballistic Lambert.
   Caveat: with ephemeris orbits the "orbit" changes between epochs, so the bound is
   exact only in the frozen-element approximation — fine as a diagnostic, not a
   constraint.
3. **#3 broken-plane genome verdict: PARTIALLY USEFUL.**
   - The Eq. 7 three-impulse base parameterization (θ_i, θ_f, t_pf,1, t_pf,2, r_mid,
     φ, λ) is exactly a phase-free broken-plane/bi-elliptic search and can seed the
     time-fixed broken-plane MBH with a good midcourse-impulse location (their
     line-of-nodes initialization tip, p. 12, citing Vinh et al.: nodal points are good
     initial guesses though not optimal).
   - Fig. 4 answers "is a third impulse ever ΔV-justified?" For Earth–Mars geometry
     (β ≈ 0.66, plane changes ~1.85–3.4°) the answer is firmly **two-impulse-optimal
     time-free** — i.e., the broken-plane DSM's value in our lane comes from the
     *time-fixed phasing constraint*, not from base-solution ΔV structure. That is a
     useful negative: it predicts broken-plane gains vanish for legs whose epochs happen
     to sit near the phase-free optimum, and concentrates where phasing forces an
     off-optimum Lambert. (Fig. 4 is strictly coaxial circle-to-circle — heuristic only
     for elliptic, non-coaxial real orbits.)
   - The out-of-plane impulse *placement* question of #3 is not answered here beyond the
     parameterization; no analytic placement condition is given (contrast Iorfida 2016,
     which gives the out-of-plane primer geometry).
4. **Per-impulse caps for maintenance ΔV: NOT applicable** in-leg (see §3); the split is
   paid in phasing revolutions we cannot afford.
5. **Four-layer classification: theory only** for us (see §4).

## 7. Negatives registry entries

- Iso-impulse phasing-orbit families: structurally excluded for time-fixed, phase-fixed
  cycler legs (TOF ≥ T_0 unsatisfiable); do not revisit unless we ever design insertion/
  establishment phases.
- Four-layer polytope enumeration: counts phasing-orbit structures, disjoint from our
  family indexing; not an MBH seeding tool.
- Circle-to-circle Sec. VI.B numbers: not golden-grade (absolute radii unpublished).
- Paper contains two internal numeric inconsistencies (701.2694 vs 720.2694 d, p. 30;
  5.6109 vs 5.6108 km/s, pp. 30–31) — quote with care if ever used as goldens.
  **UPDATE 2026-06-13 (#232): both fixed in the VoR by section removal — see §8.**

---

## 8. Version-of-record diff (2026-06-13, #232)

The open-access J. Astronaut. Sci. version of record (DOI 10.1007/s40295-025-00528-0;
J. Astronaut. Sci. (2025) 72:54; 38 journal pages, PDF created 2025-10-17) was diffed
against the arXiv:2501.01583 preprint these notes were mined from. The two flagged
internal inconsistencies (errata-ledger ids `saloglu-taheri-2025-p30-total-time` and
`saloglu-taheri-2025-dv-total-mismatch`) **do not survive — both are fixed in the VoR**,
and the mechanism is a section revision rather than a digit correction:

- **The entire arXiv Sec. VI.A Earth–Mars rendezvous worked example was removed in
  production.** It is replaced by a new CAPSTONE-inspired **orbit-raising example**
  (VoR §5.4 "An Orbit-Raising Example", journal pp. 29–32, Table 5 + Figs. 19–20).
  The VoR's worked examples are now: geocentric three-impulse (§5.1), geocentric
  two-/three-AP (§5.2), Earth–Dionysus (§5.2.1), and the orbit-raising example (§5.4).
  The VoR Conclusion (§6) lists exactly these — no Earth–Mars example.
- A full-text search of the VoR for every Earth–Mars number from §5.2 above
  (5.5865, 5.5873, 313.2430, 801.4034, 6.047, 701.2694, 720.2694, 72.7306, 5.6109,
  5.6108, 759.6963, 1479.9658, 686.9658) and for the mission time "793" returns
  **zero hits**. Both defective cells are simply absent from the version of record.
- This is the routine, no-fault kind of preprint→publication revision; the underlying
  work is excellent and the defects were typesetting-class in the first place.

**New numeric content in the VoR not in the preprint mining (the replacement example,
recorded for completeness — not wired):** orbit-raising problem, Table 5 elements
initial `[15265.9970, 0.6451, 0.7553, 0.0689, 5.5234, 0]` and target
`[47312.5193, 0.8854, 0.8225, 0.0689, 5.5234, 0]` (a in km, angles in rad; CAPSTONE
JPL Horizons states); Δv_total = 1.0885 km/s; Δv_max = 0.5443 km/s (Δα = 0.5); 6-day
mission; n_p,1 = 4 phasing orbits; feasibility 15.23 < ΣN_{k,1} < 27.62; chosen
(N_{1,1},N_{2,1},N_{3,1},N_{4,1}) = (6,4,4,3); Fig. 19a periods {0.2730, 0.3413,
0.3822, 0.4892} d with impulses {0.2394, 0.1959, 0.0870, 0.1661, 0.3999} km/s; Fig. 19b
periods {0.2193, 0.3, 0.3323, 0.8064} d with impulses {0.0109, 0.3157, 0.2068, 0.3077,
0.2475} km/s. This is the same Δv-allocation machinery already assessed as
NOT-applicable-in-leg (§3, §6.4); the CAPSTONE framing is an establishment/insertion
apogee-raising case, exactly the one niche §3 flagged — still outside our epoch-anchored,
phase-fixed cycler legs. No new takeaway for our lane; nothing to wire.

**Page re-cite (preprint→VoR pagination) for the values these notes DO use:** the
self-contained geocentric goldens (§5.1 here) survive verbatim in the VoR. Two-impulse
base Δv_total = 3.9618011 km/s (Δv1 = 2.8246140, Δv2 = 1.1371871) and three-impulse
base Δv_total = 3.8641159 km/s (Δv1 = 2.9390, Δv2 = 0.6815, Δv3 = 0.2436) are printed
on **VoR journal p. 19** (Figs. 7–8 text), matching the arXiv p. 18–19 transcription
exactly. The Earth–Dionysus base restated from Ref. [40] (Δv1 = 7.521545,
Δv2 = 2.38588 km/s, t_pf = 348.46 d, θ_i = 179.27°, θ_f = 149.20°) is on **VoR p. 21**
(§5.2.1; note the VoR cites this as Ref. [40], renumbered from the arXiv's Ref. [52]).
The arXiv-only Earth–Mars numbers in §5.2 above have no VoR pagination because the
example was removed. The `tests/verify/test_dv_bracket.py` goldens (3.9618011 and its
breakdown) are therefore VoR-confirmed.
