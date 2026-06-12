# Load-bearing numbers — READ + REPRODUCE verification pass (#203)

**Date:** 2026-06-12
**Scope:** every value about to feed a test EXPECTED side (golden-wiring #187/#191/#207)
plus spot-checks of values already wired. REPRODUCED with the repo's own core
(`core/lambert`, `core/kepler.propagate`, `core/_stumpff`, `core/cr3bp`,
`search/s1l1_corrected`) or, where the repo has no matching function, with an
independent throwaway probe (closed-form / series / RK4 / scipy DOP853 in `/tmp`).
**Writeback: NONE** — verification note only; no production code or catalogue edits.
Golden discipline preserved: every EXPECTED value below is the *published* number;
my arithmetic is used only to confirm or to FLAG a suspect printed value, never to
mint an expected value.

---

## HEADLINE

**2 NEW SOURCE DEFECTS found (both BLOCK the corresponding golden until handled),
plus 1 transcription/derivation hazard, plus re-confirmation of 1 previously
flagged defect. Everything else reproduces clean.**

1. **SOURCE DEFECT — Vallado 1991 TR-91-6 Stumpff table, z = −39.47842 row (p.E-12).**
   Printed C = 5.83559577, S = 0.97444596. The repo's `_stumpff`, an independent
   power series, and the closed cosh/sinh form ALL give **C = 6.75677528,
   S = 1.05406777**. No nearby z reproduces both printed values. The other four
   Stumpff rows (z = 0, 0.57483, 39.47842, 50.0) reproduce to ≤4e-9.
   **BLOCKS:** the `_stumpff` golden wiring (§4.1 of the Vallado note) — drop the
   z=−39.47842 row or replace its expected with 6.75677528 / 1.05406777 (sourced
   instead from BMW closed form). The four good rows are safe.

2. **SOURCE DEFECT (or trailing-digit error) — Vallado 1991 TR-91-6 Kepler case
   D.3-3 hyperbolic output (p.E-14).** Inputs Ro=(0.3,1,0), Vo=(3,0,0), Dt=5,
   mu=1. Printed R=(13.9623306, **−0.1172043**, 0). Repo `propagate` (single-shot)
   AND repo stepped (1000 substeps) AND an independent high-precision RK4 all give
   R_y = **−0.11822049** (R_x = 13.96228122 agrees to ~5e-5 rel). The printed
   elements a=−0.1411563, e=8.0473056 reproduce exactly from the inputs, so the
   geometry is right; only the printed propagated R_y is wrong (Δ ≈ 1.0e-3). On an
   e=8 hyperbola propagated to r=14 DU this exceeds the report's own stated ~1e-5
   cross-machine tolerance. **BLOCKS:** wiring D.3-3 as a tight Kepler golden — use
   it only as a coarse (>2e-4) sanity check, or skip. The four other Kepler cases
   reproduce clean.

3. **DERIVATION HAZARD (not a value error) — Shakouri 2019 case-1 Remark 2.7
   t_f = 2315 s (Table 2).** The ΔV goldens J_c=1.5210, J_m=0.9878 reproduce
   EXACTLY from Eq. 16; but my closed-form perigee→apogee reconstruction gives
   t_f = 5214.5 s (the transfer-ellipse half-period), NOT the printed 2315 s. The
   ΔV values are radius-only and confirmed; t_f depends on the actual swept arc
   (departure/arrival anomalies θ₁₂=270°, θ₃₄=30°) which my reconstruction does not
   capture. **CONSEQUENCE:** wire J_c/J_m as goldens (confirmed), but do NOT wire
   t_f from a perigee→apogee assumption — the golden-wiring task must reproduce the
   full arc geometry first, or treat t_f as out of reproducible scope.

4. **RE-CONFIRMED DEFECT — Vallado 1991 TR-91-6 Interplanetary Mars row (p.E-36).**
   Printed dist 227.8e6 km with vh1=4.1745, tof=311.804 d. Hohmann from 1 AU to
   227.8e6 km gives vh1=2.94, tof=258.7 d; the printed outputs match a distance of
   277.8e6 km exactly (227.8↔277.8 transposition). The existing DO-NOT-USE flag in
   the mining note is correct and stands. Venus row checks out (tof 146.1 vs printed
   145.983). **Already flagged — no new action.**

---

## Vallado 1991 TR-91-6 Appendix E  (PRIORITY 1; gates #187/#191)

Canonical constants used: 1 TU = 806.81106492270 s, μ⊕ = 398600.5, mu=1 in DU/TU.
Tolerance target: ~1e-5 relative (author's stated cross-machine floor, p. iv).

| value | cite | reproduced with | achieved vs expected | verdict |
|---|---|---|---|---|
| Stumpff z=0 | p.E-12 | `_stumpff` | C d=0, S d=−3e-9 | PASS |
| Stumpff z=0.57483 | p.E-12 | `_stumpff` | C d=−1e-9, S d=−3e-9 | PASS |
| Stumpff z=39.47842 | p.E-12 | `_stumpff` | C d=5e-16, S d=4e-9 | PASS |
| Stumpff z=50.0 | p.E-12 | `_stumpff` | C d=2e-9, S d=−3e-9 | PASS |
| **Stumpff z=−39.47842** | p.E-12 | `_stumpff` + series + cosh/sinh | C 6.7568 vs printed 5.8356; S 1.0541 vs 0.9744 | **SOURCE-DEFECT** |
| Lambert D.4-2 short | p.E-20 | `lambert` mu=1 | worst |dv|=9.4e-7 | PASS |
| Lambert D.4-6 short | p.E-22 | `lambert` | worst |dv|=1.0e-7 | PASS |
| Lambert 5.11b short | p.E-24 | `lambert` | worst |dv|=1.2e-7 | PASS |
| Lambert D.4-1 long | p.E-20 | `lambert` | worst |dv|=9.3e-8 | PASS |
| Lambert 5.11e long | p.E-25 | `lambert` | worst |dv|=1.0e-7 | PASS |
| Lambert A423 35TU | p.E-27 | `lambert` | worst |dv|=4.7e-8 | PASS |
| Kepler BMW p210 | p.E-13 | `propagate` | dr rel 7e-7, dv 1.4e-6 | PASS |
| **Kepler D.3-3 hyp** | p.E-14 | `propagate`+RK4 | R_y −0.11822 vs printed −0.1172043 (Δ1e-3) | **SOURCE-DEFECT** |
| Kepler D.3-4 back-multirev | p.E-15 | `propagate` | dr rel 2e-8, dv 4e-8 | PASS |
| Kepler Kaplan p307 hyp | p.E-16 | `propagate` | dr rel 2e-7, dv 2e-7 | PASS |
| Kepler p225 #4.18 radial | p.E-17 | `propagate` | dr rel 2e-5, dv 7e-6 | PASS (loose; radial escape) |
| Interplanetary Mars row | p.E-36 | closed-form Hohmann | printed match 277.8e6 not 227.8e6 | SOURCE-DEFECT (re-confirmed, already flagged) |

Lambert note: short/long way maps to `prograde=True/False` per case (the solver's
short/long selection is geometry-driven); each case matched one flag to <1e-6.

## Shakouri 2019 / Iorfida 2016  (PRIORITY 2; μ pinned first)

| value | cite | reproduced with | achieved vs expected | verdict |
|---|---|---|---|---|
| Shakouri case-1 Remark 2.7 J_c | Table 2 p.14 | closed-form Eq.16 (μ=398600.x) | 1.5210 vs 1.5210 | PASS (μ-insensitive over 398600–398600.5) |
| Shakouri case-1 Remark 2.7 J_m | Table 2 p.14 | closed-form | 0.9878 vs 0.9878 | PASS |
| Shakouri case-1 Remark 2.7 t_f | Table 2 p.14 | half-period of transfer ellipse | 5214.5 s vs printed 2315 s | DERIVATION-HAZARD (see headline 3) |
| Iorfida Table 5 r_0 (e=0.3,a=1,ν=30°) | Table 5 p.1067 | perifocal conic | [0.626,0.361] vs [0.63,0.36] | PASS; confirms caption typo a=0.1 (gives 10× too small) |
| Iorfida Table 5 r_f (ν=100°) | Table 5 p.1067 | perifocal conic | [−0.167,0.945] vs [−0.17,0.95] | PASS |

μ-pin (Shakouri): the ΔV values are insensitive to μ across 398600–398600.5 (both
arc speeds scale √μ and the differences are dominated by the radius geometry), so μ
is not load-bearing for J_c/J_m here — Curtis μ⊕ confirmed adequate. The Iorfida
Table 5 m/D/x optimality values are `a`-independent (Eq.23) and require the Eqs.24–43
line-ellipse machinery the repo does NOT implement; only the perifocal endpoint
positions were reproducible in-repo. The full Table 5 must be wired by the
golden task with its own discriminant code (CANT-REPRODUCE in-repo today).

## Saloglu 2023 / 2025  (PRIORITY 3; #207 already wired — independent confirm)

| value | cite | reproduced with | achieved vs expected | verdict |
|---|---|---|---|---|
| Saloglu 2025 §5.1 two-impulse base ΔV | Table 1 p.18 | bi-tangential closed-form, Curtis μ | 3.96594 vs published 3.9618011 (gap +4.14 m/s) | PASS — confirms #207 (gap < 10 m/s is the optimiser's true-anomaly freedom; μ-pin valid: solar μ misses by orders of magnitude) |
| Saloglu 2023 Earth-Dionysus N₁=3 phasing period | Fig.6a | Eq.2/3 timing, surplus/3 | 830.16 d vs 830.16 d | PASS |
| Saloglu 2023 rev-count bound n_p≤6 | §III.E | floor(2490.48/365.25)=6 | 6 vs 6 | PASS (arithmetic; matches #207 test) |

The wired #207 `test_dv_bracket.py` assertions (`total >= 3.9618011`,
`total − 3.9618011 < 0.010`, solar-μ rejected) are all independently reproduced.
No drift.

## Spot-checks of already-wired goldens

| value | cite | reproduced with | achieved vs expected | verdict |
|---|---|---|---|---|
| Ross-Roberts-Tsoukkas 2025 μ_EM = 1.2150584270572e-2 | AAS 25-621 p.3 | `cr3bp_system` + registry GMs | 0.01215058439 vs published, rel 1.0e-8 | PASS (registry GM_Moon/GM_system; #212b) |
| Russell App-C Mars per-leg v∞ (5.248…8.046) | App-C char. table | `|APPC_LEGS v∞ vector|` vs `APPC_MARS_TRANSIT` scalar | match to 4 dp on all 7 nodes (flyby conserves |v∞|) | PASS — two transcriptions mutually consistent |
| Russell App-C mean Mars v∞ = 5.475 | App-C avg row | mean of 7 published | 5.4749 vs 5.475 | PASS |
| Arenstorf x0/vy0/T | Hairer "Solving ODEs I" B5 | scipy DOP853 rtol 1e-12, one period | closure 1.4e-9 | PASS — golden matches canonical B5 to full double precision |

**Note on attribution:** the task brief called the s1l1_corrected per-leg v∞ a
"McConaghy App-C" value; the in-repo source (`search/s1l1_corrected.py`,
`test_s1l1_corrected.py`) attributes it entirely to **Russell 2004 Appendix C #83**,
not McConaghy. The values themselves are confirmed; only the brief's label was off.

## Items the repo cannot reproduce (hand to golden-wiring task)

- Iorfida Table 5 m/D/x₁/x₂ (needs Eqs.24–43 line-ellipse discriminant; not in core).
- Shakouri t_f (needs the full swept-arc geometry, not perigee→apogee half-period).
- Saloglu phase-free fundamental-arc t_pf/θ values (need the Eq.1/Eq.7 optimiser).
  Only the algebraic by-products (phasing period, rev-count bound, base ΔV via the
  named bi-tangential construction) were reproducible — and they all passed.

---

## PDF ADJUDICATION (2026-06-12) — the two flagged Vallado rows re-read from source

The #203 pass above flagged two rows by *reproduction* but had only checked them
against our mining note. Both have now been adjudicated by re-reading the actual
report PDF (Vallado 1991, USAFA TR-91-6, Appendix E) character-by-character, twice
independently per page, the same way the McConaghy Table 7.1 row was settled.

**Headline: 0 of the 2 are our transcription errors — BOTH are genuine
SOURCE-PRINT DEFECTS. Neither row is rescued; both stay dropped.**

1. **Stumpff z = −39.47842 (p. E-12, FindCandS table) — SOURCE-DEFECT.** The PDF
   prints exactly `−39.47842  5.83559577  0.97444596`. Our note transcribed it
   faithfully; the z value is also printed exactly (not mis-printed). Independent
   reproduction gives C = 6.75677528, S = 1.05406777, and no nearby z reproduces
   both printed values. The published cell is wrong. **NOT a transcription error →
   row NOT rescued; stays dropped** from the #187 `_stumpff` golden set. The four
   other Stumpff rows re-read identically and remain usable.

2. **Kepler D.3-3 hyperbolic R_y (p. E-14, KEPLER case 4 / BMW D.3,3) —
   SOURCE-DEFECT.** The PDF prints the full input AND output state exactly as our
   note records, including R_y = −0.1172043; no single digit/sign in the printed
   R_y or any input explains the gap (inputs reproduce a, e, |R| exactly).
   Independent `propagate` + RK4 give R_y = −0.11822049 (Δ ≈ 1.0e-3, above the
   report's own ~1e-5 floor on this e=8 hyperbola). **NOT a transcription error →
   row NOT rescued; stays dropped** as a tight #187 golden (coarse >2e-4 sanity
   use only). The other Kepler cases reproduce clean.

Mining note `docs/notes/2026-06-10-vallado-1991-tr916-mining.md` annotated in place
(dated CORRECTION/ADJUDICATION blocks at §4.1 and §3 case 4).
