# #557 admission criterion: Sun-Jupiter quasi-Hilda transient-capture `quasi_cycler` candidates

Settled in writing BEFORE any sweep code is built, per the #339-style-criterion-trap
discipline and exactly as #535 did for Earth
(`docs/notes/2026-07-03-535-quasi-cycler-transient-drift-admission-criterion.md`). This note
is the single source of truth the #557 search code must implement exactly. It records the
re-derivation of #535's Earth criterion to Jupiter's timescale, per the plan
(`docs/superpowers/plans/2026-07-11-557-jupiter-quasi-hilda-transient-capture-plan.md`) §1/§3
and the catalogue-scope correction
(`docs/notes/2026-06-16-catalogue-scope-taxonomy.md`, "system-relative window" subsection,
2026-07-11).

## 0. The one thing that changes vs #535, and why

#535's Earth criterion was written in **absolute years** (1-yr floor, 10-15-yr window) only
because the Sun-Earth system-period happens to be 1 year. The criterion is really
**dimensionless, in units of one rotating-frame period of the CR3BP system under study**:
floor = 1 system-period, window = 10-15 system-periods, count 3-15 returns, geometry ratio
3.0. This is the corrected reading now recorded in the taxonomy doc. For Sun-Jupiter the
system-period is Jupiter's heliocentric period.

**Practical consequence for the code:** everything is encoded in **periods**, where one
period = `2*pi` nondimensional CR3BP time units (one full revolution of the rotating frame).
The year-to-nondim conversion is NEVER done through a hardcoded `2*pi = 1 yr` constant (that
is the exact latent bug the plan §5 / Fable review flagged — true only for Sun-Earth). Years
appear ONLY in human-readable reporting, via the explicitly-verified Jupiter period below.

## 1. Verified physical constants (checked in-repo, not recalled)

All computed from the project's own sourced constants (`core/constants.py`) and
`cr3bp.cr3bp_system("Sun","Jupiter")`, this session:

| Quantity | Value | Source / check |
|---|---|---|
| Sun-Jupiter mass ratio mu | 9.5388115e-4 | `cr3bp_system("Sun","Jupiter").mu` (DE440 system GM, Park et al. 2021); matches Koon 2001 p.29 (9.537e-4) to 3 sig figs |
| Jupiter Hill radius R_hill | 0.0682534 nondim = 0.3551 AU | `(mu/3)**(1/3)`; matches the #527 empty-region entry's own 0.0682534 |
| Jupiter orbital period T_J | **11.868 yr** | 360 / mean_motion_deg_day from SMA 5.202887 AU (Standish & Williams), / 365.25 d/yr. Agrees with the published sidereal 11.862 yr to 0.05%. (The plan's hint 11.8618 and Koon's "~12 yr" are both within rounding; **use 11.868 for reporting**.) |
| Jupiter eccentricity e | **0.04838624** | `PLANETS["J"].ecc` (Standish & Williams Table 1 J2000 mean, in-repo sourced). The plan/task's 0.0489 is a slightly-higher NASA-fact-sheet mean; Koon 2001 uses 0.0483. Adopt the in-repo sourced 0.04838624; the 0.0005 spread cannot change a collapse/survive verdict. |
| Collinear libration Jacobi | **C_L1 = 3.038761, C_L2 = 3.037489** | computed this session (brentq on dOmega/dx at the sourced mu) |

**Critical energy fact:** temporary Jovian capture (the trajectory entering the Jupiter
region through the L1/L2 necks, Koon 2001) is energetically possible only when the necks are
OPEN, i.e. **C < C_L1 = 3.0388**. The #527/#530/#531 Hilda work all sat at C = 3.14 (necks
CLOSED, 3.48x outside the Hill sphere by construction) — that is the wrong energy band for
this phenomenon and is why those entries found zero Hill encounters. The #557 scan targets
**C in [3.00, 3.038]**, strictly below C_L1.

## 2. Return-separation floor = 1 Jupiter period (plan §1)

`min_separation` = **1 Jupiter period = 2*pi nondim = 11.868 yr**. This is the strict
analogue of #535's Earth 1-yr floor (which was Earth's own 1-period noise filter). Its job is
to merge intra-episode Hill-sphere grazing wiggles into one return, not to encode any assumed
recurrence period.

Rejected alternatives (per plan §1): the Hilda 3:2 libration period (~250-300 yr) presupposes
the object is IN the resonance and librating — exactly the structure a transient-capture
screen must not assume (quasi-Hildas transition BETWEEN resonances, Koon 2001), and it would
over-merge; the comet's own 7.91-yr heliocentric period is resonance-specific and varies
across the population. T_J is population-invariant and cleanly defined. Use T_J.

Informational cross-check (NOT the floor): documented real captures recur on multi-decade
spacing (Gehrels 3 captured ~1970-1973, re-capture predicted decades later; the Hilda-Jupiter
synodic period is 23.7 yr). The 11.868-yr floor sits just below real return spacing — correct:
it filters sub-orbital grazing without merging genuine distinct captures.

## 3. Admission window = 10-15 Jupiter periods (plan §3, Option A — RESOLVED)

`window` = **[10, 15] Jupiter periods = [119, 178] yr** (10*2pi to 15*2pi nondim). This is
Option A of the plan §3 crux, **resolved before this build** per the #557 dispatch and the
taxonomy-doc correction: the "10-15 yr" catalogue figure was Earth-calibrated, is enforced
nowhere in code/schema (`validity_window` is an unconstrained ISO date pair; `validate.py`
checks only `epoch_locked` + finite `n_returns >= 1`), and the catalogue's one existing
`quasi_cycler` row (#312) already runs an 83-yr window. So generalizing to a system-relative
window is a documentation fix, not a schema migration, and it is made.

**Which `validity_window` semantics (taxonomy-doc consequence #1):** this criterion uses the
window as an **ENCOUNTER-ACCUMULATION span** (all 3-15 distinct returns must fall inside the
sliding window itself), the same semantics #535 used — NOT #312's launch-epoch-span semantics.
State this explicitly so a future reader does not compare across incompatible definitions.

A candidate is ADMISSIBLE iff some window of length in [10, 15] periods, sliding across the
full propagated span, contains between 3 and 15 distinct returns (per §2) inclusive. Report
the window bounds and the actual return epochs used; never silently cherry-pick.

## 4. Propagation horizon

**55 Jupiter periods** (= 55*2pi nondim = ~653 yr). Chosen as ~4x the window midpoint
(12.5 periods), matching #535's own horizon/window ratio (50 yr / 12.5 yr = 4x) so a 15-period
window has ~40 periods of sliding room. In nondimensional revolutions this is the SAME
~10-60-revolution regime as #535's Earth run — i.e. the same per-point CR3BP cost, NOT a
~12x multiplier — because 1 Jupiter period = 2*pi nondim, exactly like 1 Earth year. Sample
density is set per REVOLUTION (500 samples/period), never per calendar year.

## 5. Scale-invariant parameters (transfer as-is, plan §3)

- **Encounter distance** = `r_hill = (mu/3)**(1/3) = 0.0682534` nondim. Defined relative to
  the Hill sphere itself; the detector's `dist < r_hill` test transfers with no code change.
  The 35x-larger absolute Jupiter Hill sphere is fully absorbed by passing the right `r_hill`.
- **Return count** `n_returns in [3, 15]` — dimensionless, unchanged.
- **Bounded-geometry ratio** = 3.0 (loosest return's closest approach <= 3x the window's
  tightest). Dimensionless, unchanged; keep #535's "tighten if the real data clusters more
  tightly, and record why" provisional stance.

## 6. Positive control (plan §2a / §3) — the sourced anchor

**82P/Gehrels 3**, chosen because it is one of Koon 2001's two worked temporary-capture
examples, is low-inclination (i = 1.13 deg, planar-screen-suitable), and its elements are
sourced from a citable primary source.

**Sourced osculating heliocentric elements — JPL Small-Body Database (SBDB) API**
(`ssd-api.jpl.nasa.gov/sbdb.api?sstr=Gehrels+3`, retrieved 2026-07-11; solution orbit_id 19,
JPL, DE431, epoch JD 2452484.5, 143 obs, data arc 1975-2020):

    a = 4.13541208620034 AU      e = 0.1232542077837754
    q = 3.625705145656267 AU     i = 1.126554815725817 deg
    Tisserand param w.r.t Jupiter (t_jup) = 3.027   <-- SBDB-reported, phase-independent

**Conversion to the planar Sun-Jupiter rotating-frame IC** (proven #535/#523 vis-viva-at-
perihelion path): at perihelion the radial velocity is zero (so `xdot0 = 0`, matching #523's
own construction), the heliocentric speed is `v_p = sqrt((1-mu)*(2/q_nd - 1/a_nd))` in nondim
units (distance unit = Jupiter SMA), the comet is placed on the +x axis at barycentric
`x0 = -mu + q_nd`, and the rotating-frame tangential velocity is `ydot0 = v_p - x0` (inertial
minus frame rotation omega x r). With `q_nd = q/5.202887 = 0.69686`, `a_nd = 0.79483`:

    x0 = 0.69591   xdot0 = 0.0   ydot0 = 0.57308   ->   Jacobi C = 3.02943

**Independent cross-check the conversion is correct:** the Jacobi constant computed from the
constructed IC (3.02943) matches the independently-sourced SBDB Tisserand parameter (3.027) to
0.002 — the Tisserand parameter is the CR3BP quasi-invariant that approximates the Jacobi
constant, so this is a genuine independent validation of the elements->IC pipeline, not a
circular self-check. And 3.02943 < C_L1 = 3.0388: the neck is open, so capture is
energetically allowed — consistent with Gehrels 3's documented temporary Jovian capture.

**Positive-control verdict (verified this session before any scan):** the constructed anchor
IC, propagated in the planar Sun-Jupiter CR3BP, ENTERS Jupiter's Hill sphere, reaching a
closest approach of 0.0081 nondim = **0.12 R_hill** (deep temporary capture) — reproducing the
documented capture behavior of the real object from sourced elements. The pipeline (elements ->
IC -> propagation -> detector) is validated. As with #535's RH120 anchor, the anchor motivates
the search but is not necessarily itself admissible (a single deep capture is not >=3 distinct
returns); the broad scan looks for recurrent-capture ICs nearby.

## 7. ER3BP sensitivity gate (plan §5 Risk 1) — run EARLY

Before trusting any broad-scan admissible hit, the anchor (and any hit) is re-propagated
through the project's ER3BP core (`core/er3bp.py`) at Jupiter's real eccentricity
(e = 0.04838624), with an e=0 positive control first (must reproduce the CR3BP result). #535's
Earth corridor TOTALLY COLLAPSED under this exact check (e=0.0167). Koon 2001 (p.29) states
Jupiter's eccentricity "plays little role during the fast resonance transition" because
capture is a fast (<1 Jupiter period) tube-mediated transit, unlike Earth's slow multi-year
horseshoe — so Jupiter MAY be less fragile, but this is a hypothesis to TEST, not assume
(Jupiter's e is ~2.9x Earth's, so the perturbation is larger and the robustness argument must
actually win empirically). If the anchor's capture does not survive real e, that is a decisive
result: stop the broad scan and report it.

## 8. What this criterion does NOT decide (unchanged from #535 §6)

- No encounter-velocity / flyby-quality (`dv_band`) requirement — this is a discovery screen,
  not a V0-V5 gauntlet pass.
- Planar restriction: first-pass planar CR3BP (Gehrels 3 i=1.13 deg is low, and Koon 2001's
  model is planar), but a planar null does not rule out inclined-capture structure.
- Any admissible hit needs an independent Fable second-opinion pass before any
  catalogue-adjacent writeback (the #557 dispatch's explicit "second pair of eyes" convention
  for discovery verdicts), and then the standard catalogue admission pipeline.
- A clean null is registered in `data/empty_regions.jsonl` ONLY under this Option-A,
  system-period-relative criterion, with the criterion version tagged explicitly so it is
  never confused with a fixed-years or Option-B structural-empty result.
