# Russell 2004 dissertation — methodology mining (ψ, generic returns, global search, ephemeris solve)

Mined 2026-06-07 (Task #142). Methodology chapters only (Ch.2 generic-return /
ψ machinery; Ch.3 turn-angle optimization + global search; Ch.5 ephemeris
continuation / multiple-shooting / powered-SOI). This is the ORIGINAL of the
free-return / corrector machinery we reverse-engineered
(`src/cyclerfinder/search/free_return.py`, `search/correct.py`); his actual solve
differs from ours instructively.

**Source (cite exactly, no file path):**
Russell, R. P., "Global Search and Optimization of Free-Fall Cycler Trajectories
[Earth-Mars cyclers]," Ph.D. dissertation, The University of Texas at Austin,
2004. (Advisor C. Ocampo.) 268 pp.

> Clean digital typeset; equations and all tables read unambiguously. Vision
> read of front matter + Ch.2.7 (pp.41-55), Ch.3.5-3.8 (pp.71-79), Ch.5.2-5.4
> (pp.144-157).

---

## 1. The algorithm in 3 lines

CIRCULAR-COPLANAR GLOBAL SEARCH: enumerate cyclers as integer tuples `p.h.s.i`
(p synodic periods, h half-years of full/half-rev returns, s identical generic
returns, i = the multi-rev Lambert solution index); for each, solve the
multiple-revolution Lambert "generic return" on a constant-|v∞| **sphere
referenced by angle ψ** to get free-return geometry, then **analytically minimize
the maximum flyby turn angle** (`ω_MINIMAX`) over the re-initiating flybys; keep
cyclers with aphelion-ratio `AR>AR_MIN` and turn-feasibility `TR>TR_MIN`.
EPHEMERIS SOLVE: take each circular-coplanar parent and **continue** (homotopy)
to the real ephemeris via a **multiple-shooting** root-find — `5n` unknowns
(per-leg planet-referenced v∞ vector + t0 + tf), constraints = position-match +
time-continuity + flyby-Δv=0, driven by **SNOPT** (elastic mode). Mapping
verdict: **our `correct.py` IS Russell's Ch.5 multiple-shooting corrector** (same
per-leg independence, same v∞-continuity+position+time constraints); our
free-return genome is a *symmetric-ellipse special case* of his ψ generic-return
construction.

---

## 2. The algorithm stated precisely

### 2.1 The ψ referencing angle + generic-return solve (Ch.2.7, pp.41-53)
ψ = "the angular coordinate ... [of the v∞ solution's] location in the ecliptic
plane on the v∞ sphere **referenced to v_B**, with positive being aligned with
r_B" (p.41). (v_B = body velocity, r_B = body position.)
- **Three free-return types**: half-rev (odd-nπ), full-rev (even-nπ), and
  **generic** (non-integer-π, the general case). Full-rev solutions are analytic
  (lie on concentric spheres); half-rev lie on constant-diameter circles ("a
  tube"); **generic returns lie on a series of non-uniform arcs in the ecliptic
  plane** and are "entirely numerical in origin" (Ch.2.8, p.54).
- **Generic-return generation** (p.52-53): solve the multiple-revolution Lambert
  problem at discrete transfer-angle intervals (1/2°, refined to 1/24°) up to a
  max ToF (6 body periods); bin solutions by N (revs) and fast/slow; **~180,000
  solutions over 6 years** at a ~half-hour interval. For a target |v∞|,
  interpolate ψ, ToF, a between the two bracketing points that **lie on the same
  sub-family** (critical caveat — interpolating across sub-families is invalid).
  Refine each interpolated seed with a 1-D solver that fixes |v∞| and iterates on
  **(ToF, ψ)** in a Keplerian integration to machine precision.
- **The v∞ sphere**: intersections of the half-rev tubes / full-rev spheres /
  generic arcs with the constant-|v∞| sphere are the feasible post-flyby v∞+
  directions. Realistic |v∞| ≈ 0.1-0.25 × body velocity (the example uses
  |v∞|=0.1838 AU/TU ≈ **5.5 km/s** — Fig.2.22 / Table 2.3).

### 2.2 Worked generic-return tables (Tables 2.2, 2.3) — tabulated anchors
- **Table 2.2 (p.42), N=1 generic returns** (|v∞| = half the body's circular
  velocity, exaggerated): 11 intersection points, columns ψ(deg) | ToF(body
  periods) | a(AU) | N(rev signed): −139.5/0.899/0.540/+1; −96.33/1.571/1.245/−1;
  −81.04/2.406/2.146/−1; −73.19/3.336/3.114/−1; −68.37/4.299/4.099/−1;
  −65.06/5.276/5.090/−1; 74.29/5.655/2.941/+1; 78.20/4.621/2.433/+1;
  83.74/3.567/1.920/+1; 92.73/2.469/1.392/+1; 114.0/1.250/0.804/+1.
- **Table 2.3 (p.51), generic returns for |v∞|=0.1838 AU/TU (=5.5 km/s)**: 40
  solutions (ψ | ToF | a | N). Representative rows: #1 −158.6/5.944/0.668/+8;
  #14 −86.96/1.466/1.086/−1; #34 104.2/1.348/0.921/−1; #40 161.6/2.048/0.663/−3.
  (Full 40-row set transcribable if a golden is wanted — see §4.)

### 2.3 Turn-angle min-max optimization (Ch.3.4-3.5, pp.71-77; Eqs.3.2-3.8)
Given a generic return needing `f_j` flybys to re-initiate (Table 3.2 maps
`h_j → f_j`), the flybys are spaced on the full-rev circle to **minimize the
maximum turn angle** `ω_MINIMAX-j`. Closed-form (latitude/longitude on the v∞
sphere, z-axis = v_e):
- `ω_MIN = acos(cos φ_FR cos φ_GR + sin φ_FR sin φ_GR)` (Eq.3.2)
- `φ_FR = −asin(v∞/(2 v_e))` (Eq.3.3); `φ_GR = π/2 − acos(v∞1−·v_e/(v∞ v_e))` (3.4)
- `ω_a = acos(cos²φ_FR cos λ_a + sin²φ_FR)`, `λ_a = π/(f_j−2)` (Eq.3.5)
- `ω_b` from iterating Eq.3.6 for λ; `ω_c = π − 2|φ_GR|` (Eq.3.8, f_j=1 case).
- **Decision rule** (pp.75-76): f_j=1 → ω_c; f_j=2 → acos(sinφ_GR sinφ_FR);
  f_j>2 → if ω_MIN≥ω_a use ω_MIN, else compute ω_b.
- **Grouping h half-years across s identical generic returns** (Ch.3.6, p.77):
  distribute h_j to minimize max(ω_MINIMAX-j); `Σh_j = h`; for h_j≥1 more
  half-years lower ω, for h_j=0 it depends on geometry.

### 2.4 The global-search algorithm (Ch.3.7, Fig.3.9, pp.78-79)
```
Choose AR_MIN, TR_MIN, p_MAX; h_MAX=5·p_MAX, s_MAX=3·p_MAX
DO p=1..p_MAX
  DO h=1..h_MAX
    DO s=1..s_MAX
      Calculate TOF (function of p,h,s)
      IF TOF>0:
        multi-rev Lambert → 2·N_MAX+1 generic free-return solutions (i=−N_MAX..N_MAX)
        IF N_MAX>0:
          DO i=−N_MAX..N_MAX
            determine h_j for each s_j
            optimize intermediate flybys (turn-angle min-max) for each group
            compute velocities before/after each flyby + turn angles
            Calculate TR and AR
            IF TR>TR_MIN AND AR>AR_MIN: record cycler p.h.s.i + properties
```
- `AR` = aphelion ratio = max ecliptic-plane cycler aphelion / 1.52 AU (Mars
  radius). **AR>1 ⇒ the cycler can intercept Mars without a powered maneuver**
  (p.79). Used AR_MIN=0.9, TR_MIN=0.85 (Table 3.4 caption).
- **Non-iterative** broad scan (it just enumerates and discards infeasible
  p,h,s,i). Results: **24 ballistic cyclers (2-4 synodic periods), 92 ballistic
  (5-6 synodic), hundreds of near-ballistic** (Ch.3.1 abstract).

### 2.5 Ephemeris solve — continuation + multiple shooting (Ch.5.4, pp.146-156)
THE real-model transition, in order:
- **Continuation/homotopy** (5.4.1, p.146): "solve a sequence of sub-problems
  where the fidelity ... is increased ... The solution to each sub-problem becomes
  the initial guess for each successive sub-problem." Circular-coplanar solution
  is sub-problem #1. Shortcoming stated: "the problem may become infeasible at
  one or more of the intermediate steps." Cannot jump straight to the accurate
  ephemeris — "the gap ... prohibits an immediate jump" (p.147).
- **Why NOT the analytic/integer-structured transfers** (5.4.2-5.4.3, pp.147-148)
  — **DIRECTLY RELEVANT to our free-return-genome pivot**: the analytic full/half/
  generic transfers each need ≥2 integer inputs (revs, fast/slow) that "fix the
  structure of each leg". As the model moves toward the real case, "changing the
  value of the integer inputs ... is a difficult if not insurmountable obstacle
  when dealing with a gradient-based optimizer." So Russell uses the **non-analytic
  approach**: "numerically searching for the spacecraft velocity vector that leads
  to an intercept of the desired planet ... eliminates the integer programming
  problem, and the solution structure is free to morph as necessary." He notes
  analytic fixed-structure "works very well for cyclers with favorable turning
  angle requirements, such as the **S1L1 cycler**" but is too rigid in general.
- **Multiple shooting** (5.4.4, pp.149-151): each leg is made **completely
  independent** (avoids the erratic sensitivity of chaining a multi-decade
  trajectory). **Per-leg 7 params**: initial planet, final planet, t0, tf, and the
  **3 components of the initial planet-referenced spacecraft velocity**. Planets
  fixed; so **`5n` unknowns** for n legs (the 2 planet IDs are fixed). Block-
  diagonal sparse Jacobian → SNOPT.
- **Gradient method / SNOPT** (5.4.3): SQP; tried VF13 (Harwell) but it "often
  fails when the constraints become infeasible"; **SNOPT** chosen for its
  "non-linear elastic mode" that minimizes weighted constraint violations even
  when infeasible — "the continuation method can therefore continue 'walking'
  toward the accurate ephemeris model even if favorable solutions disappear along
  the path." Rejected GA/Monte-Carlo as too slow over all cases/dates.
- **Constraints (Table 5.3, p.154)** for `5n` params: `r_f(i) = r_fplan(i)` (3n
  nonlinear equality — position match), `t_f(i-1) = t_0(i)` (n-1 linear equality
  — time continuity, solved by SNOPT each major iter), `Δv(i)=0` (n-1 nonlinear
  equality). Total 5n−2 constraints, 5n unknowns.
- **Powered-SOI flyby maneuver** (5.4.5, Eqs.5.1-5.5, pp.152-155): a ballistic
  flyby fixes |v∞| (Eq.5.1) and bounds the turn `ω_req ≤ ω_avail` (Eq.5.2), with
  `ω_req = acos(v∞-·v∞+/(|v∞-||v∞+|))` (Eq.5.3) and
  `ω_avail = 2 asin(μ_plan/(μ_plan + r_p,min v∞,small²))` (Eq.5.4 — **our
  `core/flyby.max_bend`**). Instead of imposing both as constraints, replace with
  ONE conditional Δv (Eq.5.5):
  `Δv = |v∞+ − v∞-|` if `ω_req ≤ ω_avail`, else
  `Δv = sqrt(v∞+² + v∞-² − 2 v∞+ v∞- cos(ω_req − ω_avail))` —
  i.e. the gravity assist supplies as much turn as it can; the powered Δv covers
  the magnitude mismatch (Case A) or the leftover turn deficit (Case B). Maneuver
  at the SOI (zero-radius), before the flyby if v∞-<v∞+ else after; powered burn
  *inside* the hyperbola is excluded as too risky. **This is exactly the right way
  to surface "how powered is this flyby" as a continuous Δv rather than a hard
  pass/fail bend gate — a direct upgrade for our `bend_feasible` post-hoc check.**

### 2.6 Period / cycle bookkeeping (Ch.5.3, Tables 5.1-5.2, pp.143-145)
- Cyclers propagated for **7 cycles** by default (only row of Table 5.1 with all
  bold/clean residuals); for tighter optimization: 15 cycles (32 yr) for
  1-synodic, 11 (47 yr) for 2-synodic, 5 or 10 (32/64 yr) for 3-synodic, 7 or 11
  (60/94 yr) for 4-synodic.
- **Table 5.2 — integer multiples of circular-coplanar cycler repeat times**
  (yr): the 1/2/3/4-synodic-period base repeat times are **2.14 / 4.27 / 6.41 /
  8.54 yr** (× cycle count). E.g. 3-synodic × 5 cycles = 32.03 yr ≈ 32 (best
  resonance). **Confirms our project's S/L period basis: the E-M synodic is 2.14
  yr; "p-synodic period" cycler base repeat = p × 2.135 yr.** (Cross-checks the
  Jones-note finding that VEM "synodic" is a *different* 6.4-yr beat — Russell's
  Earth-Mars synodic is 2.135 yr.)

---

## 3. Maps to our X / does not map

| Russell construct | Our code / concept | Verdict |
|---|---|---|
| Multiple-shooting ephemeris solve: per-leg independent, `5n` unknowns (v∞ vec + t0 + tf), constraints position-match + time-continuity + Δv | `search/correct.py` (free vars `[t0, leg ToFs]`, residual = flyby V∞-continuity + periodicity) | **MAPS — our corrector is his Ch.5 method, simplified.** He carries the full 3-vector v∞ per leg + position match; we carry magnitudes + ToFs. His Δv=0 constraint = our magnitude-continuity residual. |
| ψ generic-return on constant-|v∞| sphere (free-return geometry emerges from Lambert+ψ) | `search/free_return.py` radial-crossing `(a,e)` genome (v∞, ToF, ν emerge) | **PARTIAL MAP — ours is the symmetric special case.** His generic return is the general (asymmetric, multi-rev, fast/slow) free-return; our radial-crossing ellipse is the symmetric E-M-E free-return (his half-rev/symmetric subset). His ψ is the descriptor we hold; his *solve* (interpolate on the v∞ sphere over 180k Lambert solutions) is NOT what we do. |
| **Non-analytic velocity-vector search (reject integer-structured transfers because structure can't morph under a gradient optimizer)** | our pivot from free-Lambert-genome → radial-crossing `(a,e)` genome (`free_return.py` docstring) | **MAPS — and VALIDATES our pivot.** §5.4.2 is the same lesson: a fixed-integer-structure transfer is brittle on the real ephemeris; let the shape morph. He even names S1L1 as the favorable-turn case where analytic fixed-structure DOES work — consistent with our S1L1 being a special, well-conditioned case. |
| Powered-SOI conditional Δv (Eq.5.5): GA supplies max turn, Δv covers the rest | `core/flyby.max_bend` + `flyby_bend_slacks` (binary feasible/not) | **PARTIAL — clear upgrade.** We compute the same `ω_avail` (Eq.5.4) but gate binary; Russell turns the deficit into a continuous Δv. Adopting Eq.5.5 would let `correct.py` report a Δv cost instead of pass/fail, and matches how our publication layer wants "how ballistic". |
| Continuation circular-coplanar → ephemeris, re-seed each sub-problem | `search/seed_ladder.py` (partial) | **MAPS** — same homotopy as Rall; Russell adds SNOPT elastic-mode to survive intermediate infeasibility (a robustness trick we lack — `least_squares` has no elastic mode). |
| SNOPT elastic mode (keep walking when constraints infeasible) | `scipy.least_squares` (no elastic mode) | **DOES NOT MAP — robustness gap.** Our solver gives up where SNOPT would keep minimizing violation. |
| Global search `p.h.s.i` enumeration + AR/TR gates (Fig.3.9) | `search/scan.py`, `search/sequence.py` (sequence enumeration) | **PARTIAL** — we enumerate sequences; his AR (aphelion≥Mars) + TR (turn feasibility) gates + the h/s half-year bookkeeping are a more structured circular-coplanar pre-filter we could adopt. |
| Turn-angle min-max optimization (Eqs.3.2-3.8) over re-initiating flybys | (absent — we check each bend independently) | **DOES NOT MAP** — same gap as flagged in the Rall note; this is the closed-form version of "minimize the worst flyby bend across the series". |
| AR = aphelion/1.52 AU ≥ 1 ⇒ ballistic Mars intercept | (implicit; we check Mars-reach in `free_return._residuals` term B) | **MAPS conceptually** — his AR is a clean scalar gate equivalent to our "orbit reaches Mars" margin. |

---

## 4. Candidate test anchors (tabulated → golden-eligible)

EXPECTED side traces to the dissertation (published) and to closed-form Lambert
geometry — usable as candidate reference anchors (circular-coplanar; confirm
before promoting to hard goldens per project golden-discipline).

**Anchor U1 — Table 2.3 (p.51): 40 generic returns at |v∞|=5.5 km/s
(0.1838 AU/TU).** Each row (ψ deg, ToF in body periods, a in AU, N signed) is a
self-consistent Lambert free-return — a clean golden for a "ψ generic-return
solver" if we build one. (Spot rows in §2.2.) Note this is **exactly the
|v∞|=5.5 km/s** that the project's S1L1 coplanar 154-d Lambert work references
(memory `project_s1l1_realeph_closure_blocker`: "coplanar 154d Lambert reaches
Vinf_E~5.55"). Table 2.3 is the parent solution set for that regime.

**Anchor U2 — Table 2.2 (p.42): 11 N=1 generic returns** (exaggerated |v∞|);
closed-form-traceable, smaller/simpler than U2.

**Anchor U3 — Table 5.2 (p.144): circular-coplanar base repeat times** 2.14 /
4.27 / 6.41 / 8.54 yr for 1/2/3/4-synodic cyclers (× cycle count). Pure
arithmetic on the E-M synodic (2.135 yr) — usable as a period-basis sanity
golden for the catalogue's S/L period fields.

**Anchor U4 — Table 3.4 + 3.5-3.11 (pp.83-92, NOT read this pass):** the
2-3-4-5-6-synodic ballistic cycler lists (AR_MIN=0.9, TR_MIN=0.85) with per-cycler
properties — the catalogue-relevant member tables. Flagged for a follow-up mining
pass if catalogue rows tracing to Russell-2004 need source values. (Appendix C,
pp.201-245, gives the full ephemeris-model cycler trajectories.)

---

## 5. Single most implementable finding (this paper)

**The powered-SOI conditional Δv (Eq.5.5) replacing the binary bend gate.** Our
`flyby_bend_slacks` / `bend_feasible` already compute Russell's `ω_avail`
(Eq.5.4 = our `max_bend`); Eq.5.5 turns the *deficit* (`ω_req − ω_avail`, or the
magnitude mismatch when feasible) into a scalar Δv via the law of cosines. Wiring
this into `correct.py` would (a) give the publication/review layer a continuous
"how-powered-is-this-cycler" Δv instead of pass/fail, (b) match Russell's exact
definition so any cross-check is non-circular, and (c) is ~10 lines reusing
existing geometry. Runner-up: his §5.4.2 rationale is the canonical written
justification for our radial-crossing-genome pivot — cite it in the
`free_return.py` docstring's "why not fixed-structure Lambert" note.

---

## 6. v4.2 backfill checks

- **center**: all heliocentric / Sun-centred; planet-referenced v∞ is per-leg
  local. No center ambiguity.
- **tof_days_bounds**: Tables 2.2/2.3 give ToF in **body periods** (×365.25 d for
  Earth-referenced returns; these are Earth-Earth return legs). The generic-return
  search caps ToF at **6 body periods**. Real transfer-leg ToF bounds for specific
  cyclers are in Tables 3.x / Appendix C (not read). If a Russell row is added,
  source its tof_days_bounds from those.
- **source_ephemeris**: Russell's ephemeris solve uses **mean orbital elements
  referenced to J2000** (Table 5.4 "Mean Elements at J2000") for the accurate
  model — NOT a numerically-integrated DE ephemeris. The circular-coplanar parents
  are `model_assumption: circular-coplanar`; the optimized cyclers are
  `analytic-ephemeris (mean-element J2000)`. Any catalogue row tracing to
  Russell-2004 must carry `source_ephemeris: mean-element J2000 (Russell 2004)`,
  not "DE4xx".

---

## 7. Honest "not extractable" list

- The 5n-parameter ephemeris cycler member lists (Tables 3.4-3.11, 5.5, 5.6) and
  Appendix C trajectories were NOT transcribed this pass (method focus). They are
  the catalogue-member source if Russell rows are added — needs a follow-up.
- Per-cycler orbital elements (a,e,i) appear in Appendix C tables only (not read).
- The analytic-gradient partial expressions (Ch.5.4.6, Eqs.5.4.6.1-3) are derived
  in pp.156-164 (read only the §5.4.6 intro, p.156-157); the full STM/partial
  formulas not transcribed.
- Fig.2.18-2.21 (ψ-vs-v∞ and ToF-vs-v∞ for N=0..15) are plots; only the Table
  2.2/2.3 interpolated intersection points are tabulated.
