# Mining Note: Iorfida, Palmer & Roberts (2016), JGCD — Perpendicular-Thrust Primer Geometry

**Source:** E. Iorfida, P. L. Palmer, M. Roberts, "Geometric Approach to the Perpendicular
Thrust Case for Trajectory Optimization," *Journal of Guidance, Control, and Dynamics*,
Vol. 39, No. 5, May 2016, pp. 1059–1068. DOI: 10.2514/1.G001525.

**Mined:** 2026-06-10 from full PDF (10 pp.).

**Verdict: USEFUL (secondary).** Not a low-thrust paper despite the title — it is *impulsive*
primer-vector theory (Lawden) for thrust perpendicular to the orbital plane. It feeds the
broken-plane half of queued sub-project #3 and gives a closed-form cross-check for our
STM-based primer diagnostics (`verify/primer.py`). It contributes nothing to the
continuous-low-thrust half of #3.

---

## 1. What problem it solves exactly

Lawden's primer vector `p(t)` obeys `p̈ = G(r) p` (Eq. 1, p. 1060) and normally requires
integrating a 6×6 STM along the transfer (the procedure our `verify/primer.py`
`primer_on_coast` implements). The paper's contribution chain:

1. **Polar decoupling (Sec. III, pp. 1060–1061).** The gravity gradient `G` is
   eigen-decomposed in the local `{ê_r, ê_θ, ê_h}` basis: eigenvalues `+2μ/r³` (radial) and
   `−μ/r³` (twice, transverse + out-of-plane) — Eqs. 6–10. Writing
   `p = p_r ê_r + p_θ ê_θ + p_h ê_h` (Eq. 12), the primer equation splits (Eq. 15) so the
   **out-of-plane component `p_h` fully decouples** from the in-plane pair.

2. **Analytic out-of-plane solution (Sec. IV, p. 1061).** `p̈_h = −(μ/r³) p_h` (Eq. 18) has
   the same form as the 2BP equation of motion itself, so on a conic transfer orbit

   > `p_h(ν) = A·x(ν) + B·y(ν)`  (Eq. 19)

   where `(x, y)` are the perifocal position components of the transfer orbit (Eq. 17). The
   coefficients follow algebraically from the boundary values `p_h0, p_hf` (Eqs. 20–21), with
   singular geometry only when `ν_f − ν_0 = kπ` (Eq. 22). **No STM, no integration.** `A, B`
   scale as `1/a` (Eq. 23), so the optimality structure is **independent of the transfer
   semimajor axis** — it depends only on eccentricity and the two true anomalies.

3. **Perpendicular-thrust case (Sec. V, pp. 1062–1066).** If both impulses are purely
   out-of-plane (`p_r0 = p_θ0 = p_rf = p_θf = 0`, so `p_h0 = ±1, p_hf = ±1`, Eq. 25), the
   whole Lawden optimality question reduces to plane geometry: `p_h`-isolines are parallel
   straight lines `y = m x + D p_h` over the transfer ellipse in the perifocal frame
   (Eq. 24). Line slope/intercept closed forms: Eq. 27 (`m_L, D_L`, same-sign boundary
   conditions `p_h0 = p_hf`), Eq. 28 (`m_I, D_I`, opposite-sign). Tangency of the
   `p_h = ∓1` line to the transfer ellipse is the optimal/non-optimal boundary; the
   line–ellipse intersection discriminant is

   > `Δ = −D² + (1−e²)m² + 2eDm + (1−e²)`  (Eq. 36; general form Eq. 35)

   with the transfer optimal iff `Δ ≤ 0` plus a direction condition (clockwise needs
   `x_0 > x_f`, anticlockwise `x_0 < x_f`; Fig. 9). The opposite-sign case uses two
   discriminants `Δ±` (Eqs. 40–41) that define *existence* (Eq. 42), with optimality decided
   by the root ordering of Eq. 43 via Table 3 + Algorithms 1–2 (p. 1066). Step-by-step
   procedures: Table 2 (p. 1065, same-sign) and Table 4 (p. 1066, opposite-sign).

4. **Profile structure theorem (Sec. V.D, p. 1067).** The `p_h` profile has at most three
   regions / three turning points; a sub-optimal trajectory has exactly two minima and one
   maximum (cf. Fig. 6). So the number of candidate intermediate out-of-plane impulses is
   bounded a priori by geometry.

### Connection to primer-vector optimality conditions

This *is* primer-vector theory — it is a specialisation, not an alternative. Lawden's four
necessary conditions (Sec. II, p. 1060) are taken as given; the paper's novelty is that for
the out-of-plane channel the boundary-value problem Eq. 3–5 (the part our code solves by STM
integration) collapses to the two-coefficient algebra of Eq. 21. The in-plane primer pair
(`p_r, p_θ`) remains coupled (Eq. 15) and gets **no** analytic solution here.

---

## 2. Applicability to the cycler programme — honest assessment

**Sub-project #3 (low-thrust / broken-plane genome):**

- **Low-thrust half: NOT applicable.** "Perpendicular thrust" means an impulsive ΔV normal
  to the orbital plane, not continuous thrust. There is no finite-burn, no thrust-arc, no
  control-law content anywhere in the paper.
- **Broken-plane half: directly applicable as a screen.** A broken-plane Earth–Mars leg is
  exactly the regime where the dominant cost is out-of-plane. Given a candidate transfer
  conic and two node true anomalies, Eq. 36 (plus the Table 2/4 procedure) answers in closed
  form whether a two-impulse perpendicular plane-change pair is Lawden-optimal or whether the
  genome should allocate an intermediate out-of-plane DSM — without any STM integration, and
  independent of `a`. This is a cheap per-genome-evaluation feasibility/optimality gate, the
  kind of pruning the broken-plane genome will need inside a population loop. Caveat: it
  covers only the *pure* out-of-plane case (`p_r = p_θ = 0` at both ends); real broken-plane
  legs have mixed in-plane/out-of-plane impulses, for which only the decoupled `p_h = Ax+By`
  diagnostic (item below) applies, not the full optimality tables.

**Existing primer diagnostics (`verify/primer.py`, `verify/primer_refine.py`):**

- `p_h(ν) = A x + B y` (Eqs. 19–21) is an **independent analytic cross-check** of the
  out-of-plane component of our numerically integrated `primer_on_coast` on any conic coast
  arc: project our STM-propagated primer onto `ê_h` and compare against the closed form.
  Different method, same model — good test material.
- The Sec. V.D structure theorem (≤3 turning points for `p_h`) bounds where
  `primer_refine._primer_peak` can legitimately find out-of-plane-driven peaks.

**Maintenance schedules:** marginal. Cycler maintenance impulses are rarely pure
out-of-plane (cf. the mixed 13/14/8 m/s split in the Genova–Aldrin note), so the full
perpendicular-case machinery seldom applies verbatim; the decoupled analytic `p_h` is the
reusable piece.

---

## 3. Published numeric examples (golden material)

### 3.1 Sec. V.C worked example — RECOMMENDED golden (Table 5, p. 1067)

Inputs (body text, p. 1067): transfer ellipse `e = 0.3`, `a = 1.0` DU, `ν_0 = 30°`,
`ν_f = 100°`, giving `r_0 = [0.63, 0.36, 0]ᵀ` DU, `r_f = [−0.17, 0.95, 0]ᵀ` DU.

Published outputs (Table 5, p. 1067):

| Case | `m` | `D` | `x_1` | `x_2` |
|------|-----|-----|-------|-------|
| `p_h0 = p_hf` (Eq. 27 + Eq. 39 roots) | −0.74 | 0.82 | 0.08 | −1.29 |
| `p_h0 ≠ p_hf` (Eq. 28 + Eq. 43 roots) | 2.85 | −1.42 | 0.21 | −0.79 |

Published verdicts: same-sign case is non/sub-optimal in **both** directions (point lies in
the `Δ > 0` region of the `D–m` map, Fig. 14a); opposite-sign case is **anticlockwise
optimal, clockwise sub-optimal** (Table 3 + Algorithms 1–2, p. 1067).

Mining verification (arithmetic re-derivation during this mining pass, not project code):
all eight Table 5 values reproduce to the printed 2 decimals from Eqs. 27/28 and the
line–ellipse root formulas (Eqs. 34/39/43), with `Δ(m_L, D_L) = +0.364 > 0` confirming the
published non-optimal verdict, and `Δ_+ = 3.85, Δ_− = 8.70 > 0` confirming existence for the
opposite-sign case. The example is internally consistent and suitable as an independent
golden for any implementation of Eqs. 24–43.

**Caption typo warning:** the Table 5 caption prints `a = 0.1 DU`; the body text says
`a = 1.0` DU and the printed `r_0, r_f` only reproduce with `a = 1.0`. (The `m, D, x`
results are `a`-independent anyway, per Eq. 23/Sec. IV.)

### 3.2 Fig. 2 / Table 1 example (p. 1062) — secondary, lower precision

Transfer: `e = 0.5`, `a = 1.0` DU, zero inclination, `ν_0 = 0°`, `ν_f = 90°` (canonical
units `μ = 1 DU³/TU²`). Impulses `ΔV_0 = [0, 0, 0.60]` DU/TU and `ΔV_f = [0, 0, −0.45]`
DU/TU (stated as 35% of the transfer-orbit speed at each endpoint). Departure/arrival
elements (Table 1): departure `e = 0.68, a = 1.58 DU, i = 19.29°, Ω = 180°, ω = 180°`;
arrival `e = 0.56, a = 1.26 DU, i = 21.37°, Ω = 270°, ω = 105.91°`. Useful as a
construction/scenario fixture (the maneuver mainly rotates RAAN), but the 2-significant-figure
ΔVs make it weaker golden material than Table 5.

---

## 4. Implementation notes if/when #3 broken-plane work picks this up

- Required state per evaluation: transfer eccentricity `e` and perifocal endpoint positions
  only. Procedure: Table 2 (p. 1065) for `p_h0 = p_hf`; Table 4 (p. 1066) for
  `p_h0 ≠ p_hf`; both are ~10 lines of arithmetic.
- Singular geometries to guard: `ν_f − ν_0 = kπ` (Eq. 22 denominator) and `D_{L,I} = 0`
  (lines through origin / coincident-line degeneracy discussed p. 1063).
- Scope limits: two-body, conic transfer arcs, elliptic (`0 ≤ e < 1`; parabolic limit
  discussed, open orbits excluded — p. 1062); both terminal impulses must be purely
  out-of-plane for the optimality tables to be valid.
- Suggested first artifact: a `verify/`-side test comparing `primer_on_coast`'s `ê_h`
  projection against `p_h = Ax + By` (Eqs. 19–21) on a conic coast, plus a golden test
  pinned to Table 5 (Sec. 3.1 above).
