# `#761`: is Kumar 2021's torus seed on the confirmed 3:4-LO family? — continuation tractability

**Task:** `#761`, re-scoped 2026-07-29 (dispatched: "get a fable subagent to do the
investigation into 761") as a mini-investigation with a real numerical-experiment
component: determine whether a genuine, traceable continuous family of
Jupiter-Europa 3:4 exterior ("Outer"/LO) resonant periodic orbits connects

* **Kumar, Anderson, de la Llave & Gunter 2021** (AAS 21-651, arXiv:2109.14815)'s
  own "arbitrarily chosen" seed point — Jacobi constant `C = 3.0041`,
  ~22,052 km Europa closest approach, **no IC table published anywhere in the
  paper** (`#750`'s direct-read finding) — the point the catalogued
  `europa-3-4-crnbp-torus-jupiter-2026` torus's cited lineage points at; and
* this project's own **confirmed 3:4-LO** (Anderson & Lo 2011 Table 1
  reproduction, `#753`→`#755` chain, reviewer-confirmed) —
  `C = 2.99163956830415` (`ANDERSON_LO_C_FLYBY`), `x0 = -1.4304078294961569`,
  `|λ| = 1036.116117`, ~1,641 km Europa closest approach.

`#750` left this question **genuinely indeterminate** (its own brief attempt hit
the `#753`-documented "fractal branch-jump" fragility). This task exploited what
`#750` did not have pinned down: BOTH exact endpoints, in the SAME dynamical
system (Kumar 2021 uses the bit-identical `μ = 2.5266448850435e-5` as Anderson &
Lo 2011 — `#750` §2 derived this from Kumar's own Table 1 `Gm` values), so a
continuation in Jacobi constant between them is a true same-system experiment.
No catalogue changes.

---

## Verdict (read this first)

**SAME CONTINUOUS FAMILY — continuation succeeds cleanly, bidirectionally, with
an independent published-number corroboration at the far endpoint.**

The confirmed 3:4-LO continues smoothly from `C = 2.99163956830415` up to
`C = 3.0041` in 24 natural-parameter steps (`ΔC = 5×10⁻⁴`), every member passing
`continue_family`'s full gauntlet (closure / period bounds / equilibrium gate /
Jacobi conservation / independent-Radau cross-check / dedup / fold and
topology-jump detection). No fold, no branch jump, no eigenvalue discontinuity
anywhere on the path. The member converged at exactly `C = 3.0041` has a Europa
closest approach of **22,035.8 km vs. Kumar 2021's own published 22,052 km**
(p.8, the `μ₃ = 0` PCRTBP value — i.e. exactly this model) — **0.073 % relative**
— independent, paper-sourced numeric evidence that the continuation lands on
Kumar's own specific seed orbit, not merely some same-energy orbit.

The previously-feared "fractal" hazard did appear, but in a now-characterized,
benign form: it is a **first-step predictor/basin-width artifact at the
weakly-unstable end** (cold starts near `C ≈ 3.004` with `ΔC ≥ 2.5×10⁻⁴` jump to
a different root on the very first step, cleanly caught by the topology-jump
gate), resolved by `ΔC = 10⁻⁴` — **not** a break in the family.

---

## 1. The numbers

### Forward continuation (`C_flyby → 3.0041`, `ΔC = 5×10⁻⁴`, ~8 s wall)

24 steps, stop reason `jacobi_bound` (i.e. it ran out of allowed range, not out
of family), 25 gauntlet-passing members. Smooth monotone evolution throughout:

| C | x0 | period | \|λ\| |
|---|---|---|---|
| 2.9916395683 | -1.4304078295 | 25.672529 | 1036.12 |
| 2.9941395683 | -1.4217559589 | 25.511335 | 515.65 |
| 2.9966395683 | -1.4130683007 | 25.415095 | 310.78 |
| 2.9991395683 | -1.4041774120 | 25.374836 | 199.13 |
| 3.0016395683 | -1.3948582362 | 25.348685 | 109.01 |
| 3.0036395683 | -1.3870752433 | 25.319539 | 62.09 |

Per-step corrector residuals 5×10⁻¹⁴ – 1×10⁻¹¹; `x0` strictly monotone; period
varies by < 1.4 % over the whole walk; `|λ|` decays smoothly 1036 → 62 with no
jump anywhere (each adjacent ratio ≤ ~1.16⁻¹).

### Endpoint member at exactly `C = 3.0041` (secant-extrapolated + Newton-corrected)

```
x0      = -1.3852484456241585
ydot0   =  0.598839400267831        (ydot0_sign = +1, half_crossings = 6)
period  = 25.31211964876615         (period/2pi = 4.028549 -- cf. 4.085910 at C_flyby)
C       = 3.0041 exactly; crossing residual 4.97e-14
lambda  = 54.589750588974 (real saddle; Barden vs _planar_floquet agree to 2e-10)
Radau cross-check: closure + Jacobi conservation pass, dJ = 4.0e-15
Europa closest approach = 0.0328356 nondim = 22,035.8 km
  (sampling-converged: 22035.95 / 22035.77 / 22035.75 km at n = 6000/20000/80000)
  vs Kumar 2021 p.8's own stated 22,052 km (mu3=0) -> 0.073% relative, 16 km
```

### Reverse continuation (hysteresis check, `3.0041 → C_flyby`)

* `ΔC = 5×10⁻⁴` and `2.5×10⁻⁴`: **fail on the very first step** — with no secant
  history, the zeroth-order predictor (`x0_pred = seed x0`) lands the corrector
  on a *different* converged root (`x0 ≈ -1.5798`, `T ≈ 18.85` vs. the family's
  ~25.32) and `continue_family`'s topology-jump gate correctly stops the branch.
  This is the direction-dependent, cold-start face of the `#750`/`#753`
  "fractal" behavior — a basin-width issue, not family structure.
* `ΔC = 10⁻⁴`: **clean** — 125 members walk all the way down (~40 s), and the
  final direct correction at `C_flyby` re-lands on
  `x0 = -1.4304078294961544` vs. the confirmed `-1.4304078294961569` —
  **agreement to 2.4×10⁻¹⁵ (machine precision)**, `|λ| = 1036.1161170` (3×10⁻⁹
  relative from the confirmed 1036.116116696). No hysteresis: forward and
  reverse trace the same curve (midpoint spot-check at `C = 2.9976395683`
  agrees with the forward member to 10 digits).

### Interior-point check (extension past Kumar's C)

From the `C = 3.0041` member, the family continues smoothly upward at
`ΔC = 10⁻⁴` to at least `C = 3.0141` (`|λ|` decaying 54.6 → 5.5, `x0`
→ -1.3407): Kumar's "arbitrarily chosen" point is an ordinary **interior**
member of the family, not an endpoint or near-fold point — consistent with the
authors' own "arbitrarily chosen" framing.

---

## 2. Why this succeeded where `#750`'s brief attempt failed

`#750`'s attempt continued from **the catalogued torus's own numerical seed**
(`x0 = +1.61`, `C = 2.9040`, at this project's DE440-derived `μ`) — a *third*
object (the two-body-limit proxy, ~409,500 km Europa approach) at a *different*
mass ratio, and lost convergence around `C ≈ 2.926`. This task instead continued
**between the two paper-anchored endpoints at the papers' own shared `μ`**
(`ANDERSON_LO_MU`, which `#750` §2 showed is bit-identical in both papers) —
a shorter, better-conditioned path (`ΔC = 0.0125` total) along a strongly
unstable but numerically well-behaved stretch of the family, with the confirmed
3:4-LO's exact IC as the anchor. The fractal-sensitivity hazard is real but
turned out to be confined (on this path) to cold-start first steps near the
weakly-unstable end, where the corrector's Newton basin is narrow relative to
the per-step `x0` drift (~`3.9 × ΔC`).

---

## 3. What this means for the torus's transport-utility question (`#761`'s parent framing)

1. **The torus's cited lineage is now verified family-consistent**: Kumar
   2021's seed orbit — the object the catalogued
   `europa-3-4-crnbp-torus-jupiter-2026` row's `seed_lineage` cites — IS a
   member of the same continuous family as the confirmed Anderson-Lo 3:4-LO.
   The `#750` "plausible, unconfirmed common ancestry" hypothesis is now
   **confirmed** at the family level (the specific-orbit level was already
   answered NO by `#750` — different energies, different members).
2. **`#754`'s homoclinic machinery is directly relevant in principle**: the
   family member at the torus seed's own energy (`C = 3.0041`) is a genuine
   real saddle (`|λ| ≈ 54.6`), so the `Wu ∩ Ws` homoclinic-connection approach
   `#754` built for 3:4-LO applies at that energy *in kind*. **But nothing is
   finished by this result**: `#754`'s existing connection was computed at
   `C_flyby = 2.99164`, not at `C = 3.0041` — a connection at the torus seed's
   own energy would need to be built fresh (and the ~19× weaker instability
   means slower manifold departure/arrival; tractable-looking, not guaranteed).
   This task deliberately built nothing — that is a separate decision.
3. **Honest remaining gap (unchanged from `#750`)**: the catalogued torus's own
   *numerical* seed is still the project's two-body-limit proxy
   (`C = 2.9040` at project `μ`, ~409,500 km approach), NOT Kumar's `C = 3.0041`
   orbit. This task connects the two *papers'* points; whether the catalogued
   proxy seed's own family connects onto this branch (across both a `ΔC ≈ 0.1`
   span and a 0.054 % `μ` difference) remains open — `#750`'s attempt at that
   harder question broke near `C ≈ 2.926` and nothing here re-answers it. Any
   connection-building should therefore target the paper-anchored `C = 3.0041`
   member (now exactly known: `x0 = -1.3852484456241585`), not the proxy seed.

---

## 4. Code delivered

* `src/cyclerfinder/search/jovian_resonant_families.py`: sourced constants
  `KUMAR_2021_C = 3.0041` and `KUMAR_2021_CLOSEST_APPROACH_KM = 22052.0`
  (both verbatim from Kumar 2021 p.8), new function
  `continue_34lo_to_kumar_c()` (the full forward continuation + exact-C
  endpoint correction, with the cold-start basin caveat documented in its
  docstring), and a `#761` module-docstring update.
* `tests/search/test_jovian_resonant_families.py`: two new tests —
  sourced-constant check, and the full finding as a standing regression
  (continuation reaches the bound with ≥ 20 monotone members; endpoint is a
  real saddle at exactly `C = 3.0041` with `|λ| ≈ 54.59`,
  `x0 ≈ -1.3852484456`; closest approach reproduces the published 22,052 km
  to < 0.2 %). Golden expected side is Kumar's published number, never our own
  computation. ~10 s wall — not marked slow.
* Exploratory probes were ad hoc scratchpad scripts (not committed); everything
  evidentiary above is reproducible from `continue_34lo_to_kumar_c()` and this
  note's recorded parameters (reverse walk: seed the `C = 3.0041` member,
  `direction = -1`, `d_jacobi = 1e-4`, `half_crossings = 6`, `ydot0_sign = +1`).

## Verification

* `uv run ruff check` + `ruff format --check` on both changed files: clean.
* `uv run mypy src tests` (project canonical, strict config): clean, 823 files.
* `uv run pytest tests/search/test_jovian_resonant_families.py -q`: 45/45 pass
  (~55 s).
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: run before the
  commit touching `data/OUTSTANDING.md` (see commit message).
