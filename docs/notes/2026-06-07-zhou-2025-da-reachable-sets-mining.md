# Zhou, Armellin, Qiao & Li 2025 — single-impulse reachable set via DA polynomials (method mining)

Mined 2026-06-07 (full mine; triage in
`docs/notes/2026-06-07-marginal-papers-triage.md` §3 graded this FULL-MINE
WORTHY). Companion to the primer-vector diagnostic
(`docs/notes/2026-06-07-primer-vector-diagnostic.md`, #144) and the n-body
harness (`docs/notes/2026-06-06-nbody-silver-rungs.md`). Primer answers WHERE to
burn; a reachable set answers WHAT a bounded burn buys.

**Source (cite exactly, no file path):**
X. Zhou, R. Armellin, D. Qiao, X. Li, "Single-Impulse Reachable Set in Arbitrary
Dynamics Using Polynomials," arXiv:2502.11280v1 [astro-ph.IM], 16 Feb 2025.
Beijing Institute of Technology / University of Auckland. 34 pp.

> Clean digital typeset; equations, figures, and all four tables read
> unambiguously. Vision read of the full paper (pp.1-34), including references.

---

## 1. The method in 5 lines

Given a fixed-epoch velocity impulse of **arbitrary direction** and bounded
magnitude `‖Δv‖ ≤ Δv_max` under **arbitrary** ODE dynamics `dx/dt = f(x)`
(explicitly the CRTBP, not just two-body), build by Differential Algebra (DA) a
high-order Taylor map of the final state `x ≈ T_x(α,β)` in the two impulse
angles (elevation `α`, azimuth `β`), the impulse fixed at `‖Δv‖ = Δv_max` since
the RS boundary lives on the max-impulse sphere (Eq.5). Automatic Domain
Splitting (ADS) cuts the `(α,β)` domain into sub-domains so each Taylor poly
stays accurate. Project the 3-D RS onto an auxiliary plane orthogonal to the
nominal velocity (Eqs.8-20) to get a 2-variable curve `(x_p, z_p)(α,β)`; the RS
**boundary** is the *envelope* of that family of curves, found by root-solving
the envelope equation (Eq.34) reduced to a one-variable polynomial along each
ray (Eqs.43-45, `scipy.optimize.fsolve`). A high-order **local polynomial
approximation** expanded from a few anchor points (Eqs.46-51) replaces most
root-solves, cutting envelope-solve CPU time by >84%.

---

## 2. The method stated precisely (equations + pages)

### 2.1 DA + ADS preliminaries (Sec. II, pp.4-5)
- DA = computer algebra on truncated Taylor polynomials; functions carry their
  derivatives, so a single DA integration yields the *n*-th-order Taylor
  expansion of the flow about an expansion point (Fig.1, p.4). Compose, invert,
  solve nonlinear systems, differentiate/integrate all in the poly algebra.
- **ADS** (Fig.2, p.5): assumes Taylor coefficients decay exponentially with
  order; predicts the (n+1) coefficient size by exponential fit to estimate
  truncation error; **splits** the domain into two sub-domains whenever the
  estimated error exceeds threshold `ε` during propagation. Output = a
  collection of localized Taylor polys, each over an auto-determined subset of
  the initial domain. This is the validity-domain / domain-splitting handling.

### 2.2 RS model in state space (Sec. III.A, pp.5-9)
- Dynamics general nonlinear ODE `dx/dt = f(x)`, `x=[r;v]∈R⁶` (Eq.1-2, p.5).
- Impulse `Δv₀` in spherical coords (Eq.4); on the boundary `‖Δv₀‖ = Δv_max`
  (Eq.5, p.6) → `Δv₀ ≈ T_Δv₀(α,β)` (Eq.6), `α∈[−π/2,π/2]`, `β∈[−π,π]`.
- Final state poly `x ≈ T_x(α,β)` (Eq.7). **Integrator: 8th-order Runge-Kutta
  (RK78), rel tol 1e-12, abs tol 1e-12** (p.6).
- **Auxiliary plane** orthogonal to nominal velocity, transform matrix
  `T(x̄)` built from angular momentum `h̄=r̄×v̄`, velocity `v̄`, `h̄×v̄`
  (Eq.8, Fig.3 p.7). Project both maneuvered and nominal position; plane
  constraint `δy_p = y_p − ȳ_p = 0` (Eq.11). RS on plane parameterized by
  `(x_p, z_p)`.
- Two ways to get the map `(α,β) → (x_p,z_p)`:
  - **Partial map inversion** (Sec.III.A.1, pp.8-9): add `δt_f` as a DA
    variable, integrate the τ-normalized ODE `dx/dτ = f(x)·(t_f−t₀)` once
    (Eq.15); augment + invert the map (Eqs.16-17) to solve `δt_f(α,β)` enforcing
    `δy_p=0` (Eq.18), giving `x_p≈T_{x_p}(α,β)`, `z_p≈T_{z_p}(α,β)` (Eqs.19-20).
    **One integration only — fast, but the inversion can hit a numerical
    singularity in highly nonlinear regions (e.g. NRHO perilune).**
  - **Newton's iteration** (Sec.III.A.2, pp.9-10): avoids inversion; iterates
    the final epoch `t_f^{(k+1)} = t_f^{(k)} − (1/A)·δy_p` (Eqs.21-24) until
    `|δy_p| ≤ η`. Robust but needs multiple integrations (slower). Practical
    recipe (p.10): try inversion first; substitute inverted map back into direct
    map (should be identity) — if any 2nd+-order coeff > 1e-3, fall back to
    Newton.

### 2.3 RS model in observation space (Sec. III.B, p.11)
Same machinery, but the target variables are the line-of-sight azimuth `γ` and
elevation `δ` from a single observer (Eqs.25-26), giving polys `γ≈T_γ(α,β)`,
`δ≈T_δ(α,β)` (Eqs.27-28). Everything downstream (envelope, root-find, local
poly) is identical with `(x_p,z_p)` replaced by `(γ,δ)`.

### 2.4 Envelope-boundary extraction (Sec. IV, pp.11-17) — the root-finding step
- ADS splits `D={α∈[−π/2,π/2], β∈[−π,π]}` into `N` sub-domains `D_i` (Eqs.29-31);
  each gives sub-polys `x_{p,i}≈T(δα,δβ)`, `z_{p,i}≈T(δα,δβ)` over the unit box
  `δα,δβ∈[−1,1]` (Eqs.32-33).
- **Envelope equation** (envelope theory, Eq.34, p.12):
  `g(δα,δβ) = ∂x_p/∂δα · ∂z_p/∂δβ − ∂x_p/∂δβ · ∂z_p/∂δα = 0`.
  The four partials are themselves Taylor polys (Eqs.36-39), so `g≈T_g(δα,δβ)`
  is a 2-variable poly (Eq.40); its roots form the characteristic curve whose
  image is the sub-envelope.
- **Reduction to a 1-variable poly** (Eqs.42-45, Fig.4 p.14): take initial
  guesses on the box boundary `δα,δβ=±1`; for each, write the point in polar
  `(θ,r)` with `θ=atan2(δα,δβ)` (Eq.43), `r=√(δα²+δβ²)` (Eq.44); **fix θ,
  solve for r** along that ray: `T_g*(r) = T_g(r sinθ, r cosθ) = 0` — a
  one-variable polynomial solved with **`scipy.optimize.fsolve`**, `r*∈[0, r̄]`
  (`r̄=√2` on a corner ray, `=1` on an axis ray). If no feasible root, the guess
  itself is a characteristic point (boundary-limited).
- Caveat acknowledged (p.15): valid only if a **unique** root per ray; ADS keeps
  each sub-envelope simple and they observed no multiplicity, but cannot prove
  uniqueness. Mitigations: finer ADS threshold, or pick root nearest the guess.

### 2.5 The >84% speedup — local polynomial approximation (Sec. IV.B, pp.15-17)
Root-solving cost grows linearly with the number of guesses. Instead: solve the
envelope eq numerically only at a few **anchor points** `(θ_k, r_k*)`; add
variations `δθ_k, δr_k*` as DA variables (Eqs.46-49) and **invert** the augmented
map enforcing `g=0` to get a **local poly** `δr_k* ≈ T(δθ_k)` (Eq.50). For any
nearby azimuth, the root is then *evaluated analytically*: `θ≈θ_k+δθ_k`,
`r*≈r_k*+T(δθ_k)` (Eq.51, Figs.5-6). Sub-envelopes from all sub-domains are
unioned and interior points deleted with the **alpha-shape** method (a convex-
hull generalization) to get the final envelope. Overall procedure = **Table 1
pseudocode (p.18, 11 steps)**; inputs `Δv_max, f(x), t₀, t_f, ε, n, N_p, N_a`.

---

## 3. Library / tooling verdict

| component | tool | provenance (p.) | license note |
|---|---|---|---|
| DA core | **DACE** (Differential Algebra Core Engine), C++ | github.com/dacelib/dace (p.4) | check before adoption; DACE is open-source but verify terms |
| Python wrapper | **DACEyPy** | github.com/giovannipurpura/daceypy (p.4) | the actual implementation language used here |
| 1-var root solve | `scipy.optimize.fsolve` | p.14 | already a project dep (scipy) |
| interior-point deletion | **alphashape** (PyPI) | pypi.org/project/alphashape (p.17) | new dep, MIT-style typical |

**Integrator-coupling verdict (the load-bearing tooling question for us):** the
method needs the dynamics `f(x)` *as a DA-evaluable right-hand side* so the flow
can be Taylor-expanded by integrating in the DA algebra (Eq.15). It does **not**
mandate their specific RK78 — any integrator that steps DA numbers works — but it
**does require the force model to be expressed in DACE's DA arithmetic, not as a
black-box numeric propagator.** This is the key adoption blocker: our
**REBOUND/IAS15 rails propagator cannot host the expansion** as-is, because
REBOUND integrates ordinary floats, not DA polynomials. To use this method we
would have to re-express our heliocentric n-body force model (Sun + planets on
DE440 rails) inside DACEyPy and integrate it there — REBOUND becomes the
independent cross-check, not the host. That is real work, but the force model is
simple (point-mass third bodies on prescribed rails), so it is tractable; the
hard part is feeding the DE440 rails ephemeris into a DA-time integration (the
body positions are functions of the DA time variable when `t_f` is itself a DA
variable, Eq.13/15).

---

## 4. Benchmark results (golden-eligible as method cross-checks)

CRTBP, Earth-Moon `μ = 0.0121505839` (Eq.52), two NRHOs near EM L2. User params
(Table 4, p.19): `Δv_max = 10 m/s`, poly order `n = 6`, `N_p = 51` guesses/bound,
`N_a = 6` anchors/bound. Laptop: 2.5 GHz i5-12500H, 16 GB.

**Linearly stable NRHO (Tables 2,5; Figs.7-20).** `t_f = 1T ≈ 2.667 TU`
(~9.84 d). ADS → 29 sub-domains at `ε=1e-5`. Local-poly vs exact agreement
~1e-12 absolute.

| metric | value | page |
|---|---|---|
| CPU all-points (exact envelope solve, 29 sub-domains) | **12.3184 s** | p.21 |
| CPU anchors + local-poly approximation | **1.9497 s** | p.21 |
| **CPU reduction** | **84.17%** | p.21 |
| max relative error `P = d_max²/S_RS` over 100 epochs | **0.0658%** | p.25 |
| average relative error over 100 epochs | **0.0032%** | p.25 |
| fraction of epochs with `P < 0.01%` | **91%** | p.25 |
| single-epoch best `P` (Sec.V.B example) | **7.4350e-8** | p.23 |

Error-index `P = d_max²/S_RS × 100%` (Eq.55, p.22), `S_RS` = RS area, `d_max` =
max distance of an outside-MC point from the envelope. Table 5 (p.26) lists 12
epochs; e.g. `t_4`: `S_RS=1.7466e-4`, `d_max=7.5522e-6`, `P=3.2656e-5`%.
Threshold sensitivity (short arc `t_f=0.1T`, Table 6 p.29): tightening `ε`
1e-5→1e-6 grows sub-domains 8→18 and cuts relative error **76.14%** (6.8445e-4 →
1.6328e-4).

**9:2 NRHO (Table 3; Figs.23-24, pp.29-30).** `t_f = 0.5T` reaching perilune —
highly nonlinear. **Partial map inversion fails / overestimates the RS** at
perilune (the singularity, Fig.24 — its red dashed envelope balloons far outside
the MC cloud); **Newton's iteration gives the correct tight envelope.** This is
the documented inversion-singularity failure mode and the reason both methods are
provided.

**Observation space (Figs.25-30, pp.30-32).** Observer on 9:2 NRHO views target
on stable NRHO at `t_f=1T`, `ε=1e-3`, 13 sub-domains; relative error
**4.8443e-7** (`S_RS=0.0986 rad²`, `d_max=2.1857e-4 rad`).

These NRHO/CRTBP numbers are **NOT directly golden for us** (different dynamics,
different orbits) but the *method invariants* are reusable acceptance gates if we
re-implement: (a) local-poly vs exact-solve envelope agreement ~1e-12; (b) MC
cloud fully contained by the envelope with `P` well under 0.1%; (c) the
inversion-vs-Newton split must reproduce the perilune overestimate failure.

---

## 5. Maps-to-our-X verdict

**Concept fit: strong. "TCM reachability for an Aldrin cycler leg."** Inputs we
already have: a rails / real-ephemeris node state (heliocentric `[r;v]` at an
encounter or mid-leg epoch) and an n-body force model. Pick a TCM epoch on a
leg, set `Δv_max` to the maintenance/TCM budget (our Aldrin maintenance is
~2.9 km/s total but per-TCM corrections are far smaller), propagate to the next
node epoch `t_f`, and the method returns the **reachable encounter-plane
footprint** — the 2-D set of arrival positions (on a plane through the target)
attainable for that bounded single burn. The encounter-success test becomes
geometric: **does the target node (Mars B-plane point / next-node position) lie
inside the reachable envelope?** If yes, the leg is recoverable with ≤Δv_max from
that epoch; the envelope size vs budget gives a **tolerance band** directly.

**How it strengthens existing chains:**
- **V2 / encounter-success tolerance bands:** today encounter success is a
  closure threshold (km at the node). A reachable set converts that into a
  *budgeted* statement — "the node is inside the Δv_max=X m/s single-burn
  reachable set from epoch t" — which is a stronger, falsifiable evidence
  artefact than a bare miss-distance, and is exactly the "given ≤X m/s at time
  t, what is reachable?" framing the triage wanted.
- **#148 recoverable-Δv:** the primer diagnostic flags coast 0 (E→M) as
  IMPROVABLE (`max|p|≈1.122`) but *cannot quantify* recoverable Δv. A reachable
  set is the complementary quantifier: the smallest `Δv_max` whose envelope just
  contains the target node is a lower bound on the single-burn correction cost
  from that epoch — a concrete recoverable-Δv number where primer gives only a
  yes/no.

**Scope estimate: MED.** Not LOW because the force model must be re-expressed in
DACEyPy DA arithmetic (REBOUND cannot host it) and DE440 rails must be fed into a
DA-time integration; not HIGH because the algorithm is fully specified
(Table 1, 11 steps), the deps are small (DACEyPy + scipy + alphashape), and our
n-body force law is simple point-mass-on-rails. A first cut: two-body heliocentric
DA RS (validate against an analytic/MC golden), then add third-body rails terms,
then wire to a real node state. The envelope/root-find/local-poly layer is pure
geometry on the polys and ports verbatim.

---

## 6. Honest limits (verify-on-use)

- **Single impulse, fixed epoch only.** The paper solves *only* scenario 1 of
  Xue et al.'s three (fixed epoch, arbitrary direction). Multi-burn or
  free-epoch TCM windows need extension (compose maps across burns, or add the
  epoch as another DA variable — the latter is partly shown via `δt_f` but only
  to hit the auxiliary plane, not as a free maneuver epoch). Our maintenance
  schedule is multi-impulse, so this is a genuine gap, not a wrapper away.
- **Expansion validity over long coasts.** All examples are *short* arcs
  (`t_f ≤ 1T ≈ 9.84 d` for the stable NRHO; `0.5T` for 9:2). ADS controls the
  truncation error but at a cost: sub-domain count climbs with arc length
  (Fig.15, 8→~29 over one period) and *explodes* near strong nonlinearity
  (perilune forces Newton's iteration). **Our cycler legs are 100-600 days**
  (E→M ~132 d, M→E ~599 d in the Aldrin note) — far longer and multi-rev-scale.
  This is the *same long-arc fragility* the primer-vector survey flags
  (`...primer-vector-diagnostic.md` lines 110-118, Guzman 2002). Expect heavy
  ADS splitting / possible non-convergence over a full 599-day leg; the method
  is most credible applied to a *short final segment* near the target node
  (the last few days of approach), which is also where a TCM actually lives.
- **Multi-rev / uniqueness caveat.** The 1-variable ray root-solve assumes a
  unique root per azimuth (p.15); the authors observed none but cannot prove
  it. Long multi-rev arcs are exactly where multiplicity could appear —
  unverified for our regime.
- **Inversion singularity is real and silent-ish.** Partial map inversion
  *overestimates* the RS at high nonlinearity (Fig.24) rather than erroring
  loudly; the 1e-3 identity-check (p.10) is the guard. Any re-implementation
  must ship that check, else a perilune-like region yields a plausible-looking
  but wrong (too large) reachable set — a false "recoverable" verdict.
- **Reachable set is a boundary approximation, not exact** (Table 1 note, p.18):
  accuracy is set by `ε` and order `n`; tightening both costs CPU. Fine for an
  evidence band, not for a hard feasibility proof at the envelope edge.
