# Vasile & Campagnola — DFET method body mining (transcription, MGA global search)

Mined 2026-06-07 (Task #142). The DFET (Direct Finite Elements in Time) method
body (pp.3-21), deliberately skipped in prior scans (font-broken). Complements
the data-table notes `2026-06-05-vasile-tables-retranscription.md` and
`2026-06-05-vasile-hiraiwa-scan.md`. Informs our low-thrust v2 transcription and
the global-search / phasing machinery.

**Source (cite exactly, no file path):**
Vasile, M. & Campagnola, S., "Design of Low-Thrust Multi-Gravity Assist
Trajectories to Europa," Journal of the British Interplanetary Society (JBIS)
(arXiv:1105.1823). 27 pp.

> Font-broken PDF, but the Read-tool vision render of pp.3-20 was clean; all
> equations and Tables 1-6 read unambiguously off the rendered pages.

---

## 1. The algorithm in 3 lines

THREE-STAGE pipeline: (BS1) solve the **phasing problem** — given a GA sequence,
find launch/flyby dates that minimize Δv, via a closed-form `t_INITIAL(k)` (Eq.19)
+ branch-and-prune over integer rev-counts `k` (merit Eq.21); (BS3) refine the
promising sequences as **optimal impulsive transfers with one DSM per leg**
(Lambert subarcs, `6N−6` unknowns, min Σ|Δv|, Eq.50) — impulsive optimum used as
the low-thrust first guess; (DFET/DITAN) transcribe the low-thrust legs with
**Direct Finite Elements in Time** — weak-form weighted-residual collocation
(Eqs.52-67) → sparse NLP solved by SQP. The Jovian tour uses a separate
**Synchronous-Orbit-Tour** integer program (BS2, Eqs.22-41). Mapping verdict:
**DFET is an ALTERNATIVE low-thrust transcription to our Sims-Flanagan (Ozimek/
Yam) one — does NOT map to our SF leg model; the BS1 phasing solve and BS3
one-DSM-per-leg DOES map to our free-return phasing + the Forge's impulsive
seeding.**

---

## 2. The method stated precisely

### 2.1 Dynamic model (Ch.2.1, Eqs.1-6)
- `ṙ = v`; `v̇ = ∇U(r) + ∇U_B(r) + u/m` (Eq.1) — point mass under primary +
  third-body disturbing potential + thrust control `u`.
- `U(r) = μ/|r|` (Eq.2); third-body `U_B(r) = μ_B(1/|d| − ⟨d,ρ⟩/ρ³)` (Eq.3).
- State `x = [r_x,r_y,r_z,v_x,v_y,v_z,m]^T`, control `u = [u_x,u_y,u_z]^T` (Eq.4)
  — **same 7-state as our SF leg and as Ozimek.**
- Thrust bound `T_min ≤ |u| ≤ T_max` (Eq.5); `T_min = 1e-4·T_max` (avoids Hessian
  singularity in min-mass problems — cf. Ozimek's mass-leak λ, same purpose).
- `ṁ = −u/(I_sp g_0)` (Eq.6).

### 2.2 Swing-by + power models (Eqs.7-15, Table 1)
- Linked-conic, zero-SOI: `r_i = r_0 = r_P` (Eq.7); `|ṽ⁻| = |ṽ⁺|` (Eq.8);
  `⟨ṽ⁺,ṽ⁻⟩ = ṽ² cos β` (Eq.9); `β = 2 asin(μ̃/(ṽ² r̃_p + μ̃))` (Eq.10) — **the
  bend formula = our `max_bend`, = Russell Eq.5.4, = Jones Eq.2.**
- SEP power: `F_max = η_e P_in F_sp` (Eq.11); `P_in = min(P_in*, P_max)` (Eq.15);
  array power `P_eff = η_S (P_1AU/R_S²)[1 − C_T(T_S − T_0)] cos α_ss` (Eq.13);
  panel temp `T_S = (S_0 α_s cos α_ss/(R_S² σ κ ε_s))^0.25` (Eq.14).
- **Table 1 power-system constants**: η_e=0.9, η_S=0.9, C_T=3e-4 K⁻¹, T_0=290 K,
  κ=1.8, ε_s=1.0, α_s=0.8, T_max=423 K, P_SS=300 W.

### 2.3 BS1 — the phasing problem (Ch.3.1.1, Eqs.16-21)
"Given a sequence of swing-bys we want to find the set of dates ... that minimises
the Δv ... the **phasing problem**." Coplanar-circular reduced model:
- Per-leg transfer time = odd multiple of Hohmann time:
  `ΔT = (2m+1)ΔT_H = (2m+1)π sqrt((a_A+a_B)³/(8μ_s))` (Eq.16), `m ∈ ℕ`.
- Planet angular positions `θ(t) = θ_0 + ωt` (Eq.17); phasing constraint
  `θ_B,FINAL − θ_A,INITIAL = (2k+1)π` (Eq.18) → closed form
  **`t_INITIAL(k) = ((2k+1)π − ω_B ΔT − Δθ_0)/(ω_B − ω_A)`** (Eq.19).
- Launch date per phase: `t_Launch(k) = t_INITIAL − Σ_{j<i} ΔT^(j)` (Eq.20).
- Merit `F(k^(1)..k^(N)) = Σ (t_Launch^(i) − t_Launch^(1))²` (Eq.21), minimized by
  **fast branch-and-prune over the N integer rev-counts `k`**. Coded as **BS1**.
> **This is the analytic phasing solve that our free-return corrector approximates
> numerically** — Eq.19 gives the date(s) closing the phase in closed form for the
> circular-coplanar model. Directly comparable to `free_return._residuals` Term A
> (planet-span = transfer-angle), but Russell/Vasile solve it for the date, we
> drive a residual.

### 2.4 BS2 — Synchronous Orbit Tour integer program (Ch.3.1.2, Eqs.22-41)
The moon-tour leg (also maps to resonant cyclers): a sequence of resonant orbits
`ρ_i = n_i/m_i = T_i/T_M`, `a_i = ρ_i^{2/3} a_M` (Eq.22). Per-GA outgoing speed
constraints (Eqs.23-24); deviation angle `β_0 ≤ β_i ≤ β_MAX` (Eq.25), β_MAX from
Eq.10. Feasible resonance band `ρ_min,i ≤ ρ_i ≤ ρ_MAX,i` (Eqs.29-31). Total tour
time `T = T_M Σ n_i` (Eq.34). **Optimal SOT = min Σ n_i** s.t. resonance/feasibility
(Eqs.35-40) — an integer program; pruned by GA efficiency `η_i = β_i/β_MAX`
(Eq.41). Coded as **BS2**.

### 2.5 BS3 — optimal impulsive transfer w/ one DSM per leg (Ch.3.1.3, Eqs.48-51)
Each phase split into two Lambert subarcs joined by a **deep-space maneuver
(DSM)**. Unknowns `X = [t_P1..t_PN, t_DSM,1..t_DSM,N-1, r_DSM,1..r_DSM,N-1,
r̃_p,1..r̃_p,N-2]^T` → **`6N−6` components** (Eq.49). Minimize
`F_BS3 = Σ Δv_i = Σ |v⁺_DSM,i − v⁻_DSM,i|` (Eq.50) s.t. flyby Eqs.8-9 +
departure/arrival velocity caps (Eq.51). Key assumption: "**optimal impulsive
transfers can be regarded as a limit case of a minimum mass, low-thrust transfer
with no limit on the thrust level**" → the min-Δv impulsive solution is the
low-thrust first guess (cites Okutsu for the same). **This is exactly our Forge's
ballistic/impulsive-first seeding philosophy, and the SF-leg zero-thrust limit.**

### 2.6 DFET transcription (Ch.3.2, Eqs.52-67) — THE namesake method
Decompose the trajectory into M **phases**, each with time domain `D^j`, dynamic
vars x, control vars u, parameters p. Per phase:
- objective `J^j = φ^j(x_0^b, x_f^b, t_f, p) + ∫ L^j(x,u,p) dt` (Eq.52);
- dynamics `ẋ − F^j(x,u,p,t) = 0` (Eq.53); path constraints `G^j ≥ 0` (Eq.54);
  boundary `ψ^j ≥ 0` (Eq.55); **inter-phase link constraints** (Eq.56).
- Each phase's time domain split into **N finite time elements**
  `D^j = ∪ D^j_i(t_{i-1}, t_i)`; within each element states/controls expanded in
  **polynomials of order (p−1)**: `[x; u] = Σ_{s=1}^p f_s(t)[x_s; u_s]` (Eq.57),
  `f_s ∈ P^{p-1}(D^j_i)` (Eq.58). A finite element = sub-domain `D^j_i` + parameter
  subset `[x_s, u_s, p]`.
- **Weak-form weighted-residual transcription** of the dynamics (Eq.59):
  `∫_{t_i}^{t_{i+1}} {ẇ^T x + w^T F^j} dt − w_{i+1}^T x_{i+1}^b + w_i^T x_i^b = 0`,
  i=1..N−1, with **generalized weight (test) functions**
  `w = Σ_{s=1}^{p+1} g_s(t) w_s` (Eq.60), `g_s ∈ P^p(D^j_i)` (Eq.61).
- Find `x_s ∈ ℝ^{p·m}`, `u_s ∈ ℝ^{p·n}`, p, boundary `[x_f^b, x_0^b]` satisfying
  Eq.59 + `G^j ≥ 0` (Eq.62) + `ψ^j = 0` (Eq.63).
- Each integral replaced by a **q-point Gauss quadrature sum, q = p** (Eq.64);
  state continuity at matching points `x_i^b = x_{i+1}^b` (Eq.64); algebraic
  constraints collocated at Gauss nodes (Eq.65).
- Assemble all phases → **NLP** `min J(y)` s.t. `c(y) ≥ 0`, `b_l ≤ y ≤ b_u`
  (Eqs.66-67), `y = [x_s, u_s, x_0^b, x_f^b, t_0, t_f, p]^T`. Implemented in
  **DITAN** (Direct Interplanetary Trajectory ANalysis); uses **JPL DE405**
  ephemeris; sparse SQP.

### 2.7 Design procedure (Ch.4, p.14) + system budget (Table 2)
BS1 → pool of candidates (min phasing-constraint violation) → BS3 (C_3 cap, no
arrival-velocity limit) → DITAN (DFET) per leg → BS2 SOT → add SOT to transfer →
DITAN whole-trajectory optimize. **Table 2 spacecraft**: 1500 kg total, max
C_3=3.16 km²/s², 7×150 mN ion thrusters, 40 kW EOL, 133 m² array, 187 kg
propellant, 416 kg solar array, 773 kg total SEP.

---

## 3. Maps to our X / does not map

| Vasile-Campagnola element | Our code / concept | Verdict |
|---|---|---|
| DFET weak-form weighted-residual transcription (Eqs.59-65), polynomial elements + Gauss collocation | `core/sims_flanagan.py` (midpoint-impulse SF transcription) | **DOES NOT MAP — a DIFFERENT transcription.** We use Sims-Flanagan (impulses + Kepler coast); DFET uses polynomial finite elements + weighted residuals. Both end in a sparse SQP NLP, but the leg discretization is fundamentally different. DFET is the alternative we did NOT pick; note it if SF accuracy ever proves limiting. |
| 7-state `[r,v,m]`, thrust bound, `T_min=1e-4 T_max` anti-singularity | our SF leg 7-state; Ozimek's mass-leak λ | **MAPS (state) / PARALLEL (anti-singularity).** Same state vector; the `T_min` floor plays Ozimek's λ role. |
| BS1 phasing: closed-form `t_INITIAL(k)` (Eq.19), branch-prune over integer revs | `search/free_return.py` Term A (planet-span = transfer-angle), `search/phase_match.py` | **MAPS — analytic version of our phasing residual.** Eq.19 solves the circular-coplanar phasing in closed form; we drive it numerically. Could seed our corrector with Eq.19 dates. |
| BS3 one-DSM-per-leg impulsive optimum as low-thrust first guess (`6N−6` vars, min Σ|Δv|) | Forge ballistic/impulsive seeding; SF zero-thrust limit | **MAPS — same "impulsive optimum seeds low-thrust" philosophy.** Our SF leg's all-zero schedule = ballistic; the Forge seeds impulsively. BS3's explicit min-Δv-impulsive→low-thrust limit is the documented justification. |
| Swing-by bend β = 2 asin(μ/(ṽ²r_p+μ)) (Eq.10) | `core/flyby.max_bend` | **MAPS exactly** (4th independent source for this formula: also Russell 5.4, Jones Eq.2, Ozimek Eq.17). |
| BS2 Synchronous Orbit Tour integer program (resonance ρ=n/m, min Σn_i) | `search/resonance.py`, `search/tisserand.py` | **PARTIAL MAP — relevant to #76 Tier 2 / resonant cyclers.** The SOT is a resonant-flyby-chain integer program on the Tisserand plane (Fig.12); our resonance/tisserand modules cover the same space for the cyclers context. |
| DITAN whole-NLP, JPL DE405, SQP | our verify stack; no NLP optimizer | **DOES NOT MAP** — we have no DFET/DITAN-equivalent whole-trajectory low-thrust optimizer. |

---

## 4. Candidate test anchors (tabulated)

EXPECTED side traces to the published paper (DITAN-optimized, DE405). Per
golden-discipline these are candidate reference anchors (font-broken source →
"confirm vs clean JBIS before promoting to hard golden"; the vision read was
clean but it's still a single raster source — same caveat as the prior
Vasile-table notes).

- **Tables 3 & 4 (pp.16-17): EVMEJ first-guess vs DITAN-optimized encounter
  dates** (MJD2000). Solution 1: Earth dep 3692→3719, Venus 4130→4112, Mars
  4328→4284, Earth 4558→4667, Jupiter arr 5625→5962. Solution 2: Earth dep
  1356→1343, Venus 1794→1751, Mars 1990→1942, Earth 2218→2190, Jupiter arr
  3232→3128. (FGS = BS1 mean-elements; optimized = DITAN DE405.)
- **Table 5 (p.20): Ganymede SOT** — apocenter (FGS→opt km) + resonance n:m:
  1.90e7→1.80e7 / 75:1→69:1; 4.97e6→5.92e6 / 10:1→7:1; 3.13e6→3.19e6 / 5:1→5:1;
  2.22e6→2.22e6 / 3:1→3:1.
- **Table 6 (p.20): Europa SOT** — pericenter + resonance: 6.40e5→6.36e5 / 3:1;
  6.25e5→6.25e5 / 5:2; 6.02e5→6.05e5 / 2:1; 5.72e5→5.77e5 / 8:5.
- **Table 1 (power constants), Table 2 (mass budget)** — system-design anchors.
- (Tables 7 & 8, the full TOF + GA-characteristics, already transcribed in
  `2026-06-05-vasile-tables-retranscription.md`.)

> Note: this is an Earth→Jupiter/Europa MGA mission, NOT an Earth-Mars cycler.
> All anchors are low-thrust-method / resonant-tour anchors, not cycler-catalogue
> rows.

---

## 5. Single most implementable finding (this paper)

**BS1's closed-form phasing solution `t_INITIAL(k)` (Eq.19) as an analytic seed
generator for our phasing/corrector.** For the circular-coplanar model our
free-return corrector currently *searches* for the epoch that closes the
planet-span constraint; Eq.19 gives it in closed form (per integer rev-count k),
and Eq.21's branch-and-prune merit ranks the candidate phasings. Wiring Eq.19
into `search/phase_match.py` / `free_return.py` as the initial t0 (instead of a
guessed seed) would make the corrector start on-phase and converge faster — and
gives a non-circular cross-check against our numeric phasing. Runner-up: BS3's
explicit "min-Δv impulsive = low-thrust limit" is the documented rationale for
the Forge's impulsive-first seeding.

---

## 6. v4.2 backfill checks

- **center**: cruise modeled in J2000 mean-ecliptic Sun-centred; capture/tour in
  J2000 mean-equatorial Jupiter-centred (Fig.1). Multi-center, but the cyclers-
  relevant cruise leg is heliocentric — no catalogue center ambiguity for any
  Earth-Mars-relevant extraction.
- **tof_days_bounds**: BS1 caps per-leg ToF at odd multiples of Hohmann time
  (Eq.16); the SOT total `T = T_M Σ n_i`. Specific leg ToFs are in Table 7 (prior
  note). No new bounds for catalogue.
- **source_ephemeris**: **DITAN uses JPL DE405**; BS1 uses **3-D analytical mean
  orbital elements** (the FGS↔optimized discrepancy in Tables 3-6 is exactly the
  mean-element-vs-DE405 gap — the paper says so, p.15). Any anchor promoted must
  carry `source_ephemeris: DE405 (DITAN-optimized)` for the optimized columns,
  `mean-element analytic (BS1)` for the FGS columns.

---

## 7. Honest "not extractable" list

- No Earth-Mars cycler trajectories (this is an Earth→Europa MGA mission). The
  EVMEJ / EVMEJ sequences include Earth-Mars legs but as part of a Jupiter-bound
  MGA, not a repeating cycler.
- DFET polynomial basis functions `f_s`, `g_s` are given structurally (P^{p-1},
  P^p spaces) but the specific orders p used per phase are not tabulated.
- Thrust/coast time histories are plotted (Figs.6-7 right panels) not tabulated.
- The font-broken source means all numeric anchors carry the standard
  "confirm-vs-clean-copy" caveat even though this vision read was unambiguous.
