# Mining Note: Liang, Yang, Li, Bai & Qin (2024) — Callisto-Ganymede-Europa Triple Cyclers

**Source:** G. Liang, H. Yang, S. Li (NUAA), X. Bai (Rutgers), L. Qin (NUAA),
"Callisto-Ganymede-Europa Triple Cyclers," *Journal of Guidance, Control, and Dynamics*
(Engineering Note), DOI: 10.2514/1.G008387. Page numbers below refer to the 22-page
author draft ("Draft, Journal of Guidance, Control, and Dynamics"). A portion was
presented as Paper IAC-23,C1,9,9,x76777 at the 74th IAC, Baku, 2–6 October 2023
(Acknowledgments, p. 20).

**Mined:** 2026-06-11. Deep mine, with the Forge Phase 6 Jovian empty-region registry
(`jovian-IEG-vilm-2026-06-09`, `jovian-perm-vilm-2026-06-09`) as the explicit lens.

**Verdict: USEFUL (primary — first sourced numeric anchors for the Jovian moon-cycler
bucket).** The paper publishes **four concrete ballistic CGE triple-cycler members**
(three in an idealized circular-coplanar two-body model, one in a SPICE-ephemeris
patched-conic model), each sustained ≥ 10 cycles (~1000 days), with per-flyby V∞, ToF,
and initial moon phases tabulated. The idealized members are fully reproducible from
printed data (up to a μ_Jupiter choice). This is exactly the "new sourced data" trigger
for the §6b re-sweep gate — but the members live in a **topology and capability class
our swept genome does not cover** (repeated-moon CGCEC sequences with multi-revolution
legs), so the existing zero-rev EMPTY verdicts remain valid *as conditioned*; see §6.

---

## 1. The near-resonance idea (Sec. II.A, pp. 2–3)

Lynam & Longuski's triple cyclers [9] (their spelling "Lynan") sit on the exact 1:2:4
Laplace resonance of Io–Europa–Ganymede, using the gravity assist to precess the
argument of periapsis at the constant Laplace-configuration repeat rate. Callisto is
not in the Laplace resonance, so that machinery is unavailable for any Callisto-bearing
triple. Liang et al. instead exploit a **near-resonance** induced by the Laplace
resonance:

- Synodic periods: S_C,G = **12.5232 d** (Callisto–Ganymede), S_G,E = **7.0509 d**
  (Ganymede–Europa); ratio S_C,G/S_G,E = 1.7761 ≈ 7/4 (Eq. 1, p. 3).
- Mismatch per resonance period: M_C,G,E = 4·S_C,G − 7·S_G,E = 50.0928 − 49.3563 =
  **0.7365 d** (Eq. 2, p. 3) → the three-moon angular configuration approximately
  repeats every **~50 days** (the "near synodic period" T_S).
- T_S is near-commensurate with all three orbit periods: T_S/T_C ≈ 3, T_S/T_G ≈ 7,
  T_S/T_E ≈ 14 (Eq. 3, p. 3) — every ~50 d the moons return near their inertial
  positions, not just their relative geometry.

Why C-G-E and not Io: Io's radiation environment is hostile to long-duration spacecraft
(p. 2, citing [13] Campagnola et al., Europa Clipper); Europa's orbit is higher; and the
Laplace resonance makes S_C,G ≈ S_G,E-cluster commensurabilities possible (p. 2,
citing [12] Lynam/Kloster/Longuski for the ≈16:3 Callisto–Ganymede vs Ganymede–Io
synodic ratio, restated p. 3).

Because the resonance is imperfect, phase differences accumulate ~100° in ~150 days
(p. 8), which is what kills the naive "visit all three per revolution" design and what
makes the previous Laplace-resonant methods [9, 10] **inapplicable** here — the relative
phases never reset exactly (p. 2).

## 2. Models — fidelity tags (Sec. II.B–C, pp. 3–6)

All dynamics are **two-body Keplerian about Jupiter** in both models; flybys are
impulsive patched-conic velocity rotations in both. Per-result tags:

1. **Idealized model** (results §4.1–4.3 below): circular, planar moon orbits; mean
   motions from Table 1, p. 4 — Europa **1.7693**, Ganymede **0.8782**, Callisto
   **0.3765 rad/day**; orbit radii derived from μ_Jupiter + mean motion (μ_J value not
   printed — the one reproduction gap); moon initial phases free parameters.
   Tag: `circular-coplanar, patched-conic, two-body(Jupiter), ballistic, multi-rev`.
2. **High-precision ephemeris model** (result §4.4): moon states from SPICE kernel
   **JUP365.bsp**, ecliptic J2000.0 Jupiter-centered inertial (p. 4) — non-circular,
   non-coplanar moons; spacecraft legs still Keplerian about Jupiter with impulsive
   flybys (the paper calls it the "patched-conic ephemeris model", pp. 2, 18).
   Tag: `ephemeris(JUP365), patched-conic, two-body(Jupiter), ballistic, multi-rev,
   3D (moon orbits only)`.

Gravity-assist model (Sec. II.C, Eqs. 5–13, pp. 4–6, after [15] Chen/Baoyin/Li 2014 and
[16] Li/Baoyin/Jiang): minimum-Δv flyby-defect formulation in a P-ξηζ frame aligned
with the inbound V∞; turn bound sin(δ/2) = μ/((R+h)V∞² + μ) (Eq. 7, p. 5); optimal
defect angles ψ* = arctan2(ν₂,ν₃), δ* = min(δ_max, γ) (Eqs. 12–13, p. 6). Minimum flyby
altitude constraint **50 km** (p. 11). A solution is declared ballistic when the
optimizer's residual Δv < **1e-8 m/s** (p. 13).

## 3. Construction method (Sec. III, pp. 6–13) — vs our moontour genome

### 3.1 Initial conic guess (III.A)

Three conditions bound the spacecraft orbit period (pp. 7): (1) T_S is an integer
multiple of the spacecraft period; (2) the starting moon's period is an integer ratio
to the spacecraft period; (3) period large enough that apojove clears Callisto's orbit.
Given T, the semi-major axis follows from a = (T√μ_J/2π)^(2/3) (Eq. 14, p. 7);
eccentricity is replaced by the perijove radius, bounded r_p ∈ [R_Jup, r_Europa]
(Eq. 15, p. 7). In all published examples the spacecraft period equals **Callisto's
period** (T = T_C ≈ 16.69 d; 3T ≈ T_S).

Discrete structure encoding: in the first revolution the spacecraft notionally cuts
each moon's circle twice; a 0/1 flag per moon (flyby above/below the x-axis of the
perifocal frame, spacecraft departing perijove on the +x axis) defines the initial
phase structure, e.g. **1-1-1** or **1-1-0** for (Callisto, Ganymede, Europa)
(Fig. 2, p. 8). Choosing the structure fixes the moons' initial phases.

### 3.2 Switched-double-cycler strategy (III.B) — the key trick

A triple cycler is decomposed into **two double cyclers that alternate**: one C-G
double cycler completed in one ~50 d cycle period, then a switch to a C-E double
cycler for the next ~50 d, giving a triple-cycler repeat period of **2 cycle periods
≈ 100 d** (p. 9, Fig. 3). Rationale (pp. 8–9): when only a *pair* configuration
(Europa–Callisto or Ganymede–Callisto) must recur, the accumulated phase error is far
smaller than for the triple configuration, and switching between the two double cyclers
is achieved mainly by shifting encounter epochs (i.e. changing the semi-major axis via
the flyby) rather than rotating the apse line. **All transfer arcs are forced to more
than one revolution** (multi-rev Lambert legs); the multi-rev coast is explicitly "a
process to reshape the phase of the spacecraft" (p. 9) — the gravity assist mainly
modifies semi-major axis, not argument-of-periapsis precession as in the Laplace-
resonant designs (pp. 18–19).

For the worked sequence **CGCEC** (one cycle = Callisto→Ganymede→Callisto→Europa→
Callisto, ~100 d; Fig. 4, p. 10), the approximate flyby epochs per cycle n_cycle are
(Eq. 16, p. 11):

> t_c1 = (6n−5)T + t_c0; t_g1 = (8n−4)S_C,G + t_g0; t_c2 = (2n+2)T + t_c0;
> t_e1 = 22n·S_G,E + t_e0; t_c3 = (6n+1)T + t_c0

with t_e0, t_g0, t_c0 the ToFs from perijove to each moon in the initial conic guess
and T = T_Callisto. **Probable typo:** the t_e1 term is printed with S_G,E, but Fig. 4
labels the same interval "11·2 S_C,E" and the text (p. 10) says the 50-d structure
repeat equals "4 C-G synodic periods or 11 C-E synodic periods"; 22·S_C,E ≈ 99.2 d fits
the epoch ordering t_c2 < t_e1 < t_c3, while 22·S_G,E ≈ 155 d does not. Read S_G,E in
Eq. 16 as S_C,E (Callisto–Europa synodic, ≈ 4.511 d from Table 1 mean motions).

### 3.3 Optimization (III.C)

Cycle-by-cycle **local** optimization, chained: variables = flyby epochs t_GA1…t_GAn
and return epoch t_f (Eq. 17, p. 11; the starting epoch/V∞ are inherited from the
previous cycle and not optimized); objective = Σ‖Δv_f‖ of the flyby defects (Eq. 18);
legs solved as (multi-rev) Lambert problems between epochs; solver = IMODE (improved
multi-operator differential evolution [17]). Initial epoch guesses from Eq. 16,
searched within several days; when Europa's accumulated phase error makes the guess
stale (total Δv rising to ~200 m/s after several hundred days), the Europa epoch is
shifted by ±integer synodic periods and the variable range expanded until ballistic
(p. 12). Every local solution found was ballistic, so no global problem is ever solved
(p. 13).

Ephemeris-model additions: the start epoch is found by a traversal search minimizing
Δθ(t) = |θ^i_CG − θ^r_CG(t)| + |θ^i_GE − θ^r_GE(t)| (Eq. 19, p. 12) — i.e. find the
real-ephemeris date whose relative phase structure best matches the idealized member;
and the optimization variables are switched from epochs to per-leg **ToFs** with bounds
±2–6 days around the idealized values (Eq. 20, p. 13), because epoch variables fail to
converge after several cycles in the ephemeris model.

### 3.4 Capability comparison vs our Phase-6 genome

| Capability | Our swept genome (`jovian-*-vilm-2026-06-09`) | Liang et al. |
|---|---|---|
| Leg type | single-ellipse free-return, zero-rev legs | multi-rev Lambert legs (forced > 1 rev) |
| Sequence class | simple 3-moon loops (I-E-G, E-G-C, …), no repeats | repeated-moon CGCEC (Callisto thrice per cycle) |
| Flyby role | bend-feasibility gate on a fixed ellipse | SMA-changing GA, per-flyby defect minimized |
| Phasing | epoch grid, period_k 1–3 | near-resonance epoch law (Eq. 16) + DE per cycle |
| Moon model | circular coplanar | circular coplanar AND SPICE JUP365 |
| Verdict | EMPTY above 6 km/s V∞-floor at 26.8 km/s best | ballistic members at max-V∞ 6.99 km/s (high perijove) |

Their method does not merely tighten our search — it **changes topology class**
(repeated-moon sequences) and **adds a leg capability we lack** (multi-rev legs used as
phase-reshaping reservoirs). It strictly subsumes nothing of ours and ours nothing of
theirs; the two searches are disjoint in genome space. See §6.

## 4. Published members — complete transcription

All four members share sequence **CGCEC** per ~100 d cycle, propagated 10 cycles
(~1000 d), ballistic (residual Δv below the 1e-8 m/s threshold; ephemeris member
max 1.0383e-7 m/s per cycle, p. 18 — "can be neglected"). Per 10-cycle solution:
Europa ×10, Ganymede ×10, Callisto ×21 encounters (+1 each counting the initial conic
guess) (p. 15).

**Caveat on the "Flyby Altitude" columns** (paper's own, p. 16): the altitudes in
Tables 3, 5, 7 are computed from the *required defect Δv* on the moon-centered
hyperbola ignoring Jupiter's gravity — they are a reflection of how little turning is
needed, not physical periapses (e.g. 1.9e6 km vastly exceeds Callisto's sphere of
influence; that "flyby" needs essentially zero turn). In the ephemeris solution all
flybys occur above 100 km real altitude (p. 18).

### 4.1 Member A — idealized, 1-1-1 structure, high perijove (Fig. 5, Tables 2–3, pp. 14–15)

Perijove r_p = r_Eu − 10000 km ≈ **660,988 km** (p. 14). Initial spacecraft orbit:
period = T_Callisto, departs perijove at t = 0 (perifocal +x). Initial moon phases
(Table 2, p. 14): Europa **0.063748 rad**, Ganymede **0.70579 rad**, Callisto
**1.3550 rad**.

First-cycle flyby summary (Table 3, p. 15) — ToF is per-leg, days:

| # | Moon | ToF (d) | V∞ (km/s) | "Altitude" (km) |
|---|---|---|---|---|
| 0 | Callisto | 0 | 5.6730 | 1,900,851 |
| 1 | Ganymede | 31.8973 | 6.9919 | 978,172 |
| 2 | Callisto | 18.1697 | 5.6698 | 33,839 |
| 3 | Europa | 29.9343 | 4.6685 | 6,241 |
| 4 | Callisto | 19.9747 | 5.8721 | 19,765 |

Cycle ToF over 10 cycles: 99.86–100.14 d (Fig. 5d); per-pair ToF drift over 10 cycles:
C-G 31.9→33.0 d, G-C 18.2→17.0 d (≈1 d range); C-E and E-C vary ~4 d (Fig. 5c, p. 15) —
the C-G ToF is systematically more regular than C-E (p. 15).

### 4.2 Member B — idealized, 1-1-0 structure, high perijove (Fig. 6, Tables 4–5, pp. 15–17)

Same perijove as A. Initial phases (Table 4, p. 15): Europa **−0.48843 rad**, Ganymede
**0.70579 rad**, Callisto **1.3550 rad** (only Europa's flag/phase differs from A).

First-cycle flyby summary (Table 5, p. 17):

| # | Moon | ToF (d) | V∞ (km/s) | "Altitude" (km) |
|---|---|---|---|---|
| 0 | Callisto | 0 | 5.6730 | 1,900,851 |
| 1 | Ganymede | 31.8973 | 6.9919 | 978,172 |
| 2 | Callisto | 18.1697 | 5.6698 | 60,877 |
| 3 | Europa | 30.2850 | 4.4853 | 10,825 |
| 4 | Callisto | 19.6834 | 5.7914 | 36,258 |

(First two rows identical to Member A — the trajectories diverge at the Europa branch.)
Observed in both A and B: perijove drops, eccentricity grows, apse line rotates
clockwise over the 10 cycles, interpreted as easing the Europa encounter (p. 17).

### 4.3 Member C — idealized, 1-1-1 structure, low perijove (Fig. 7, Tables 6–7, pp. 17–18)

Perijove r_p = (r_Eu − 10000)/2 ≈ **330,494 km**. Initial phases (Table 6, p. 18):
Europa **0.94281 rad**, Ganymede **1.3883 rad**, Callisto **1.7936 rad**.

First-cycle flyby summary (Table 7, p. 18):

| # | Moon | ToF (d) | V∞ (km/s) | "Altitude" (km) |
|---|---|---|---|---|
| 0 | Callisto | 0 | 7.6433 | 1,329,021 |
| 1 | Ganymede | 32.2542 | 10.4922 | 629,226 |
| 2 | Callisto | 17.8127 | 7.6409 | 27,438 |
| 3 | Europa | 31.1267 | 12.0213 | 3,636 |
| 4 | Callisto | 18.9362 | 7.7838 | 12,516 |

Lower perijove ≈ doubled V∞ — confirms the high-perijove family is the low-energy one.
Across all three idealized members: ToF(G→C)+ToF(C→G) ≈ 50 d and likewise for the C/E
pair; cycle-ToF spread ≤ 0.35 d; all cycle ToFs ≈ 100 d (p. 18).

### 4.4 Member D — ephemeris model (Fig. 8, pp. 19–20)

Ballistic CGCEC, 10 cycles, **SPICE JUP365.bsp** patched-conic. Departs Callisto
**25 September 2033 18:04:43**, returns to Callisto **22 June 2036 01:44:39**
(~1000 d). Max required Δv for one cycle **1.0383e-7 m/s** (negligible → ballistic);
all flybys above **100 km** real altitude (p. 18). Cycle ToF 99.4–100.5 d; per-pair ToF
behaviour matches the idealized model (C-G less variable than C-E) (Figs. 8d–e, p. 20).
**No per-flyby table is printed for this member** — flyby epochs/V∞ exist only as
figure traces; data completeness is epoch+sequence+kernel only.

## 5. Proposed catalogue rows (NO writeback — sketches for review)

All four as members under a new family seed `liang-2024-cge-triple-family`
(Jovian sibling row alongside `hernandez-2017-jovian-ieg-triple-family`, which already
cites this paper as a corroborating source — promote that corroboration to its own
family + member rows). Common fields: center: Jupiter; sequence_canonical:
"Callisto-Ganymede-Callisto-Europa-Callisto (per ~100 d cycle)"; class: triple cycler
(switched-double-cycler construction); ballistic; demonstrated duration ≥ 1000 d
(10 cycles, not indefinite — near-resonance, paper's own caveat p. 20).

1. `liang-2024-cgcec-111-highperijove` — Member A. Model: circular-coplanar two-body
   patched-conic. Inputs: Table 1 mean motions, Table 2 phases, r_p = r_Eu − 10000 km,
   T_sc = T_Callisto, perijove departure. Anchors: Table 3 (5 flybys: ToF, V∞).
   tof_days_bounds per cycle ≈ [99.86, 100.14]. Reproducible up to μ_J choice →
   candidate for a reproduction script + V1 on ingest, V2 if our reproduction matches
   Table 3.
2. `liang-2024-cgcec-110-highperijove` — Member B. As above with Table 4 phases,
   Table 5 anchors.
3. `liang-2024-cgcec-111-lowperijove` — Member C. As above with Table 6 phases,
   r_p = (r_Eu − 10000)/2, Table 7 anchors. Useful as the high-V∞ contrast member.
4. `liang-2024-cgcec-ephemeris-2033` — Member D. Model: ephemeris(JUP365)
   patched-conic; source_ephemeris: JUP365.bsp. Anchors: departure/arrival epochs only;
   per-flyby data unpublished (figures only) → V0/V1 ceiling until/unless the authors'
   data is obtained. data_gaps: per-flyby epochs, V∞, altitudes.

v4.2 backfill checklist: center = Jupiter (all four); tof_days_bounds = per-cycle
~[99.4, 100.5] (D) / [99.86, 100.15] (A–C); source_ephemeris = null (A–C, idealized) /
JUP365.bsp (D).

**Golden candidates** (sourced-only rule satisfied): Tables 2+3, 4+5, 6+7 as
EXPECTED sides for an idealized-model reproduction (same-model target: circular-
coplanar two-body about Jupiter, mean motions exactly Table 1). Gap: μ_J unprinted —
treat trailing digits as solver/constant-dependent; validate ToF at ~0.01 d and V∞ at
~1e-3 km/s level. The ephemeris member is NOT golden-grade (no numeric per-flyby
anchors).

## 6. Jovian empty-region impact (registry: `jovian-IEG-vilm-2026-06-09`, `jovian-perm-vilm-2026-06-09`)

**Does this paper contradict our EMPTY verdicts? No.** The registry negatives are
conditional on: simple non-repeating 3-moon loops, per_leg_revs = [0,0,0], single-
ellipse free-return genome, period_k 1–3. Liang et al.'s members are (a) a different
topology class — repeated-moon CGCEC with Callisto visited 3× per cycle; (b) built on
legs of **more than one revolution** (explicitly forced, p. 9); (c) closed by per-cycle
epoch optimization, not a repeating ellipse. None of their members is a point in our
swept region. In fact the paper *corroborates* our negative's interpretation: the naive
visit-all-three-per-revolution design fails on phase drift (~100° per 150 d, p. 8), and
our zero-rev closures landing at V∞ 8.3–26.8 km/s is consistent with the feasible
low-V∞ family living elsewhere (multi-rev, ~100 d cycles, max V∞ ≈ 6.99 km/s).

**Does it REOPEN the region? Yes — via the "new sourced data" arm of the §6b gate, and
it defines the capability bar for the method arm.** Specifically:

1. **Sourced anchors now exist** in the previously anchor-free (circular-coplanar,
   Jupiter) bucket: Members A–C are fully-specified, table-anchored, same-model
   (circular-coplanar patched-conic — our sweep's own model). The registry's
   `source_anchors: "none populated"` statement is now stale for the Jovian bucket
   and should be amended to cite this paper when rows land.
2. **Justified re-sweep:** only after the genome gains (i) multi-rev Lambert legs with
   per-leg rev counts as genome genes, and (ii) repeated-moon sequence topologies
   (CGCEC-class). A re-sweep with the current single-ellipse genome would re-confirm
   the same negative and is NOT justified. The natural first campaign is a
   *reproduction-led* sweep: seed at Member A's phases/perijove, verify ballistic
   closure, then perturb (other 0/1 structures, other perijoves, other switched-pair
   sequences — the paper only publishes CGCEC; e.g. CECGC, Ganymede-started variants,
   and the Saturnian analogue are open).
3. **Saturnian transfer:** the same near-resonance construction may apply to the
   Saturn registry entries (`saturnian-titan-vilm-2026-06-09`,
   `saturnian-titan-endgame-vilm-2026-06-10`) if a comparable integer synodic-ratio
   pair exists among Enceladus/Dione/Rhea/Titan — the paper itself notes Russell &
   Strange [8] found Titan–Enceladus double cyclers. Worth a synodic-ratio scan before
   any Saturn re-sweep.
4. **Topology lesson for the negative-results registry** (Forge #172 design): record
   `per_leg_revs` and "repeated-moon allowed?" as first-class region coordinates —
   this paper is the concrete case where the empty region and the populated region
   differ *only* in those coordinates.

## 7. Relationship to the Laplace resonance — summary

- Lynam & Longuski [9] / Hernandez et al. [10]: exact 1:2:4 Io-Europa-Ganymede Laplace
  resonance; phases reset exactly each Laplace period; GA holds (a, e) and precesses ω
  by a constant angle. Breaks down in ephemeris models because real precession is
  non-constant and the tight per-leg ToF constraints can't absorb the targeting error
  (pp. 18, citing [10]).
- Liang et al.: 7:4 *near*-resonance of (C-G, G-E) synodic periods, ~50 d quasi-repeat,
  0.7365 d/period mismatch absorbed by switching between two double cyclers and by
  multi-rev legs that re-phase via semi-major-axis changes. Because the GA adjusts
  *a* rather than ω, the construction survives the ephemeris model's irregular
  precession — their headline qualitative claim (pp. 18–19), demonstrated by Member D.
- Trade: Laplace-resonant triples include Io (radiation); CGE avoids Io entirely and
  reaches Callisto, which the Laplace machinery cannot (p. 2).

## 8. References worth acquiring (precise citations from the paper's list, pp. 20–22)

- **[9] Lynam, A. E., and Longuski, J. M.**, "Laplace-resonant triple-cyclers for
  missions to Jupiter," *Acta Astronautica* Vol. 69, No. 3–4, 2011, pp. 158–167,
  DOI 10.1016/j.actaastro.2011.03.011. (Text misspells "Lynan".) **HIGH priority** —
  the three Laplace-resonant IEG triple-cycler members with different periods; would
  populate the `hernandez-2017-jovian-ieg-triple-family` sibling line with per-member
  data.
- **[10] Hernandez, S., Jones, D. R., and Jesick, M.**, "One Class of Io-Europa-
  Ganymede Triple Cyclers," AAS/AIAA Astrodynamics Specialist Conference, Vol. 162,
  Univelt, Stevenson WA, 2017, pp. 973–984. **HIGH** — our existing family-seed row is
  V0 precisely for lack of this paper's tables; also contains the initial-guess
  energy-constraint strategy Liang builds on (Sec. III.A "similar to the method in
  reference [10]", p. 6) and the ephemeris breakdown analysis Liang cites (p. 18).
- **[8] Russell, R. P., and Strange, N. J.**, "Cycler Trajectories in Planetary Moon
  Systems," *JGCD* Vol. 32, No. 1, 2009, pp. 143–157, DOI 10.2514/1.36610. (Ref list
  typo "Sttange".) **HIGH** — already a family seed (`russell-strange-2009-jovian-
  multimoon-family`); per-member double-cycler tables would feed both Jovian and
  Saturnian buckets.
- **[12] Lynam, A. E., Kloster, K. W., and Longuski, J. M.**, "Multiple-satellite-aided
  capture trajectories at Jupiter using the Laplace resonance," *Celestial Mechanics
  and Dynamical Astronomy* Vol. 109, 2011, pp. 59–84, DOI 10.1007/s10569-010-9307-1.
  MEDIUM — source of the 16:3 C-G/G-Io synodic ratio underpinning Eq. 1.
- **[6] Jones, D. R., Hernandez, S., and Jesick, M.**, "Low excess speed triple cyclers
  of Venus, Earth, and Mars," AAS/AIAA Astrodynamics Specialist Conference, Vol. 162,
  2017, pp. 3387–3398. Already mined (AAS 17-577). Note: Liang's in-text "Drew et
  al. [6]" refers to this — D. R. (Drew) Jones cited by first name; not a separate
  "Drew et al." paper. Do not create a phantom acquisition target.
- **[7] Ozaki, N., Yanagida, K., Chikazawa, T., Pushparaj, N., Takeishi, N., and
  Hyodo, R.**, "Asteroid Flyby Cycler Trajectory Design Using Deep Neural Networks,"
  *JGCD* Vol. 45, No. 8, 2022, pp. 1496–1511, DOI 10.2514/1.G006487. (In-text "Naoya
  et al." — again a first-name cite, Naoya Ozaki.) MEDIUM — asteroid-cycler bucket
  opener if/when we add one.
- **[11] Yang, H., Hu, J., Bai, X., and Li, S.**, "Review of Trajectory Design and
  Optimization for Jovian System Exploration," *Space: Science & Technology* Vol. 3,
  2023, 0036, DOI 10.34133/space.0036. LOW (review; same group) — useful citation map.
- **[17] Sallam, K. M., Elsayed, S. M., Chakrabortty, R. K., and Ryan, M. J.**,
  "Improved Multi-operator Differential Evolution Algorithm for Solving Unconstrained
  Problems," IEEE CEC 2020, pp. 1–8, DOI 10.1109/cec48606.2020.9185577. LOW — solver
  only; any global DE/MBH of ours is interchangeable.

## 9. Errata / quote-with-care list

1. Eq. 16 (p. 11): t_e1 coefficient printed as S_G,E; Fig. 4 and the p. 10 text imply
   S_C,E (Callisto–Europa synodic ≈ 4.511 d). Use S_C,E.
2. In-text author names "Lynan" [9], "Drew et al." [6], "Naoya et al." [7], and ref-list
   "Sttange" [8] are misspellings/first-name cites — map to Lynam, Jones, Ozaki,
   Strange.
3. "Flyby Altitude" columns in Tables 3/5/7 are Δv-equivalent fictions (paper's own
   caveat, p. 16) — never ingest them as physical flyby altitudes; the 50 km constraint
   applies to the optimizer, the >100 km statement only to the ephemeris member.
4. μ_Jupiter is never printed; idealized-model reproductions should treat trailing
   digits as constant-dependent.
5. The cyclers are **not** guaranteed indefinite repetition (near-resonance + chained
   local optimization); the published claim is ≥ 1000 days / 10 cycles in both models
   (Conclusion, p. 19–20). Catalogue rows must carry that bound, not "infinite".
