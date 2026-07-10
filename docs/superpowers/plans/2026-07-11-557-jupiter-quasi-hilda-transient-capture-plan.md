# #557 Sun-Jupiter Quasi-Hilda Transient-Capture `quasi_cycler` Screen — Design Plan (review-first)

> **Status: PLANNING ONLY.** This is a design document for review, mirroring how #538's
> plan (`docs/superpowers/plans/2026-07-09-538-qbcp-cross-system-cycler-correction.md`) was
> written and reviewed before any implementation code was dispatched. **No sweep script has
> been written, nothing has been run, and no `data/catalogue.yaml` / `data/empty_regions.jsonl`
> edits have been made.** Implementation model and go/no-go are for the user to decide AFTER
> reviewing this document, per the explicit instruction in the #557 OUTSTANDING entry. The
> criterion re-derivation below is the #339-style-criterion-trap discipline applied to Jupiter:
> settle the admission rule in writing first, exactly as #535 did for Earth.

**Goal:** Decide whether — and if so, exactly how — to point #535's validated Hill-sphere-return
`quasi_cycler` screen (`src/cyclerfinder/search/hill_sphere_return_detector.py`) at the
Sun-Jupiter quasi-Hilda population, which (unlike Earth's one-shot RH120) contains REAL comets
with documented, REPEATING temporary satellite captures (Oterma, P/Gehrels 3,
111P/Helin-Roman-Crockett).

**Headline recommendation (read this first):** Before any sweep is worth running, one decision
dominates all the others and must be made by the user: **the catalogue `quasi_cycler`
validity-window (10–15 yr) is Earth-timescale-calibrated and is structurally too short for
Jupiter's natural return spacing (~12–24 yr).** Keeping the literal 10–15 yr window makes the
Jupiter screen a *near-certain structural empty* — a null that says nothing about the dynamics,
only about a unit mismatch. The physically honest fix is to rescale the window (and the
separation floor) by the system timescale ratio T_J/T_E ≈ 11.86, giving a ~120–180 yr window and
a ~12 yr floor — but that is a **class-definition / scope decision** (it departs from the
catalogue's fixed mission-relevance window), analogous to the #320/#535 `quasi_cycler` scoping
discussions, not a coding choice. Everything below is conditional on resolving this.

---

## 0. What was checked before writing (grounding)

- **Detector read in full** (`src/cyclerfinder/search/hill_sphere_return_detector.py`): the
  `find_returns` / `find_admission_windows` core is pure numpy, operates on a plain
  `(t, position-relative-to-body)` time series, is explicitly documented "no CR3BP-specific code",
  and takes `r_hill`, `min_separation`, `window_lo/hi`, `n_returns_lo/hi`, `geometry_factor` all as
  **parameters**. It is already fully system-agnostic. Nothing Earth-specific is baked in.
- **#535 sweep script read in full** (`scripts/run_535_earth_transient_quasi_cycler_search.py`):
  the Earth specifics live ONLY in the constants block (mu via `cr3bp.cr3bp_system("Sun","Earth")`,
  `r_hill = (mu/3)**(1/3)`, `MIN_SEPARATION_YEARS=1.0`, `WINDOW_LO/HI=10/15`, `X0_GRID`, `C_VALUES`,
  the `2*pi rad = 1 yr` time-unit conversions) and in the CR3BP propagation. The scan/detect/report
  scaffold is reusable.
- **#535 ER3BP sensitivity check read in full** (`scripts/run_535_er3bp_sensitivity_check.py`):
  the exact structure needed for the "run real-eccentricity check EARLY" mitigation is already
  built and just needs re-parameterizing for Jupiter (`e=0.0489`).
- **#535 admission-criterion note read in full**
  (`docs/notes/2026-07-03-535-quasi-cycler-transient-drift-admission-criterion.md`): §3 explicitly
  flags the 1-yr floor as Earth-orbital-period-motivated and NOT transferable to Jupiter without
  re-derivation; §5 flags the 3× geometry factor as a provisional new numeric choice; §4 ties the
  window/count to the catalogue schema.
- **Empty-region registry verified directly** (`data/empty_regions.jsonl`, not trusting the #557
  summary): the three Sun-Jupiter Hilda-adjacent entries all use DIFFERENT methods —
  `sun-jupiter-hilda-32-mmr-dahotm` (`poincare-section-enumeration`, periodic-orbit family),
  `coorbital-resonance-unstable-manifold-encounters` (`invariant-manifold` from certified periodic
  orbits), `hilda-c3.14-homoclinic-connection-hill-encounter` (`heteroclinic`, #314 machinery). NONE
  carries the `hill-sphere-return-detection` capability tag #535's method uses, and — decisively —
  all three live at C=3.14, ENERGETICALLY BELOW the L1/L2 neck (see §3). A broad Hill-sphere-return
  sweep across the neck-open quasi-Hilda energy band is genuinely unscreened territory. Confirmed
  by reading the entries, not the summary.
- **Physical constants computed inline** (not recalled): Sun-Jupiter μ=9.5388e-4;
  R_hill=(μ/3)^(1/3)=0.06825 nondim = 0.355 AU = 53.1 Mkm (35× Earth's 1.5 Mkm in absolute terms,
  6.82× in nondim rotating-frame units); C_L1=3.0388, C_L2=3.0375; Jupiter period 11.862 yr; Hilda
  (3:2) heliocentric period 7.91 yr; Hilda–Jupiter synodic period 23.7 yr; Jupiter e=0.0489.
- **Corpus digest read** (`docs/notes/2026-06-30-digest-koon-2001-resonance-capture-jupiter-comets.md`):
  Koon-Lo-Marsden-Ross 2001 (CMDA 81:27–38) is the canonical source for the Jupiter-comet
  temporary-capture / 3:2↔2:3 resonance-transition mechanism, worked for Oterma and Gehrels 3, in
  the SAME planar Sun-Jupiter CR3BP this screen would use. States μ=9.537e-4 (p.29), Jupiter e=0.0483
  "plays little role during the fast resonance transition" (p.29) — directly relevant to the
  eccentricity-robustness risk in §5.

---

## 1. Return-separation timescale for Jupiter (criterion §3 re-derivation)

**Recommendation: the return-separation FLOOR = Jupiter's orbital period ≈ 11.86 yr.**

Reasoning, not analogy-by-default:

- #535's Earth 1-yr floor was chosen as a **numerical-noise filter** — "a gap shorter than one
  primary-secondary revolution reflects a grazing sub-cycle, not a genuine new approach" — and was
  *deliberately decoupled* from the actual (then-unknown, ultimately single-shot) return timescale.
  The floor's job is to merge intra-episode wiggles, not to encode the expected recurrence period.
  The natural unit of "one revolution of the rotating frame" is the primary-secondary orbital
  period. At Jupiter that is **11.86 yr**. This is the strict analogue of the Earth choice and
  carries the identical justification.
- **Rejected alternative — the Hilda resonance libration period (~250–300 yr):** this presupposes
  the object is *in* the 3:2 resonance and librating, which is exactly the structure the screen must
  not assume (quasi-Hildas transition *between* resonances via temporary capture; Koon 2001). Using
  the libration period as a floor would bake in the answer and is far longer than any plausible
  return spacing, so it would over-merge. Reject.
- **Rejected alternative — a small fraction of T_J (e.g. the 7.91 yr Hilda heliocentric period):**
  the comet's own heliocentric period is a candidate noise scale, but it is resonance-specific
  (7.91 yr only for exact 3:2) and varies across the quasi-Hilda population as the object hops
  resonances. T_J is population-invariant and cleanly defined. Prefer T_J.
- **Cross-check against the genuine physical recurrence (informational, NOT the floor):** the
  Hilda–Jupiter synodic period is 23.7 yr, and the documented real-object capture recurrences are
  multi-decade (Gehrels 3 captured ~1970–1973, next capture predicted ~2060s–2090s; 111P similar).
  So genuine distinct returns are spaced ~12–24+ yr apart. **The floor (11.86 yr) sits just below
  the real return spacing — correct: it filters sub-orbital grazing without merging genuine distinct
  captures.** This same fact is what forces the window-rescaling decision in §3/§0.

**One caveat to record in the criterion note if this proceeds:** unlike Earth (where 1 yr ≪ the
10–15 yr window, so the floor was dynamically irrelevant to admissibility), at Jupiter the floor
(11.86 yr) is *comparable to the entire un-rescaled catalogue window*. That is the quantitative
statement of the structural tension in §0.

---

## 2. Seeding strategy

**Recommendation: option (c) — literal real-object anchor first, then a broad neck-open-energy
quasi-Hilda scan around it** — mirroring #535's own RH120-first-then-nearby-`(x0,C)`-scan structure.
This is the project's established discipline (positive control before broad search) and it is what
distinguishes a trustworthy null from an instrumentation artifact.

**2a. Positive-control anchor — real object, sourced.** The target object must (i) have a documented
temporary Jupiter capture and (ii) be low-inclination enough for a first-pass PLANAR CR3BP screen.

- **Primary choice: P/Gehrels 3** (i ≈ 1.1°, low — planar-screen-suitable; documented temporary
  capture ~1970–1973, re-capture predicted; it is one of Koon 2001's two worked examples). **Source
  the osculating heliocentric elements (a, e, i, Ω, ω, M, epoch) from the JPL Small-Body Database
  (`ssd.jpl.nasa.gov/sbdb.cgi`, a citable primary source), NOT from memory.** Convert to a planar
  Sun-Jupiter rotating-frame IC at perihelion via vis-viva (same procedure #535 used to re-derive
  the RH120 IC from `a=0.998625, e=0.019833`, which validated to 4 sig figs — that conversion path
  is proven in-repo).
- **Fallback if elements cannot be sourced offline at build time:** use **Comet Oterma**, whose
  qualitative 2:3↔3:2 transition behavior is *published in-corpus* (Koon 2001, digest already read)
  and therefore constitutes a citable, reproducible positive control even without fresh SBDB access;
  OR construct the anchor directly in CR3BP energy space at C just below C_L1 (see §3) near the 3:2
  semimajor axis, treating the real-object match as a later refinement. **I could not fetch live
  SBDB elements within this planning pass** (no network retrieval performed for this document); the
  build task should do so as its first step, and if it cannot, fall back to Oterma-documented
  behavior rather than inventing elements.
- **Do NOT reuse the #527/#530 Hilda periodic-orbit seed (C=3.14, x0≈0.7615) as the anchor.** It is
  a certified *periodic* orbit at C > C_L1 (necks closed) that #530 already showed stays 3.48× outside
  the Hill sphere — by construction it produces ZERO Hill encounters and would be a broken positive
  control. It is the wrong dynamical object (this is the whole #535/#523/#527 lesson: the periodic
  family core is not where the close approaches live).

**2b. Broad scan.** Once the anchor reproduces a documented capture episode, scan `(x0, C)` across the
neck-open quasi-Hilda band (§3) exactly as #535 scanned around RH120. Expect narrow/filamentary
admissible structure (the #535 Earth finding was a ~15,000 km filament; chaotic-transport phase
space is generically filamentary). Use a **coarse-first** grid (see §5) before any fine pass.

**Why not (a) or (b) alone:** (a) literal-seed-only repeats the #535 RH120 mistake (a single real
object is one point and may itself be single-shot or currently between captures — it motivates but
rarely *is* the admissible candidate); (b) blind-scan-only forfeits the positive control that makes
a null trustworthy (the #535 discipline, and this project's "verify a filter with a positive control
before trusting an all-negative" rule). (c) is #535's own proven structure.

---

## 3. Adapting the admission thresholds to Jupiter's scale

The detector's five parameters split cleanly into **scale-invariant** (transfer as-is) and
**timescale-dependent** (must be rescaled):

| Parameter | Earth (#535) | Jupiter | Verdict |
|---|---|---|---|
| Encounter distance = `r_hill` | 0.01 nondim | 0.06825 nondim (0.355 AU) | **Scale-INVARIANT** — defined as `(μ/3)^(1/3)`, i.e. relative to the Hill sphere itself. The `find_returns` distance test is `dist < r_hill`; feed Jupiter's `r_hill` and it transfers with no code change. The 35×-larger absolute Hill sphere is fully absorbed by passing the right `r_hill`. |
| Geometry factor (loosest/tightest ≤ 3×) | 3.0 | 3.0 | **Scale-INVARIANT** — a dimensionless RATIO of closest-approach distances within a window. Transfers as-is. (Keep #535's "tighten if the real data clusters more tightly, and record why" provisional stance.) |
| Return count `n_returns ∈ [3,15]` | 3–15 | 3–15 | **Dimensionless (a count) — transfers as-is**, BUT only meaningful jointly with the window (see below). |
| Separation floor `min_separation` | 1.0 yr | **11.86 yr** | **Timescale-DEPENDENT — rescale to Jupiter's orbital period** (§1). |
| Admission window `[window_lo, window_hi]` | 10–15 yr | **decision required (see §0)** | **Timescale-DEPENDENT and the crux.** |

**The window decision (the crux, restated concretely).** With floor = 11.86 yr and the requirement
of ≥3 *distinct* returns each separated by ≥11.86 yr, the minimum time to accumulate 3 distinct
returns is ≥ ~24 yr — already exceeding the entire un-rescaled 10–15 yr window. **The literal
catalogue window is structurally unsatisfiable at Jupiter.** Two honest options for the user:

- **Option A (recommended if the screen runs at all): rescale the window by T_J/T_E = 11.86** →
  `window ≈ [119, 178] yr`, floor 11.86 yr, count 3–15 unchanged. This preserves the *dimensionless
  shape* of #535's criterion (window ≈ 10–15 primary-secondary orbital periods; floor ≈ 1 period;
  count 3–15) and is the physically self-consistent transfer. **Cost: it departs from the catalogue's
  fixed 10–15 yr mission-relevance `validity_window`, so any admissible Jupiter candidate would NOT
  be catalogue-admissible as a `quasi_cycler` under the current schema without a class-scope decision**
  (the ~150 yr "validity window" is not a spacecraft-mission-planning window in the sense the schema
  intends). This is a `quasi_cycler`-class-definition question for the user, exactly parallel to the
  #320/#535 scoping discussions — flag it, do not silently pick it.
- **Option B: keep the literal 10–15 yr window** → the screen is a near-certain *structural* empty.
  This is NOT worth running as a discovery search (a guaranteed null teaches nothing), and if run it
  must NOT be registered as a dynamical empty region (it would be a units artifact, not a negative
  result about Jupiter's phase space). Only defensible use: as an explicit demonstration of the
  timescale mismatch — cheaper to state in this document than to compute.

**Recommendation: resolve to Option A (rescaled window) as a scope decision BEFORE building, or
decline the search.** Do not build under Option B.

---

## 4. Build vs. reuse

**Reuse verbatim (zero changes):**
- `src/cyclerfinder/search/hill_sphere_return_detector.py` — `find_returns`,
  `find_admission_windows`, `Return`, `AdmissionWindow`. Fully system-agnostic; every scale enters
  through parameters. This is the bulk of the intellectual content and it is done.
- The CR3BP EOM / STM core (`core/cr3bp.py`), the ER3BP core (`core/er3bp.py`), the
  `preflight_search` gate + `MethodCapability` plumbing, and the empty-region registry round-trip.

**Reuse with re-parameterization only (new script `scripts/run_557_jupiter_quasi_hilda_search.py`,
structurally a clone of `run_535_...py`):**
- `system = cr3bp.cr3bp_system("Sun", "Jupiter")` instead of `("Sun","Earth")` → μ, and
  `r_hill=(μ/3)**(1/3)` recomputes automatically (formula already generic).
- Constants block: `MIN_SEPARATION_YEARS = 11.86`; `WINDOW_LO/HI_YEARS = 119/178` (Option A);
  `N_RETURNS_LO/HI = 3/15` (unchanged); `GEOMETRY_FACTOR = 3.0` (unchanged).
- Propagation horizon: #535 used 50 yr (≈ 3–4× its window). At Jupiter, ≈ 3–4× a ~150 yr window →
  **~400–600 yr propagation per point** — but per §5 (corrected 2026-07-11), this is the SAME
  ~10-60 nondimensional revolutions as #535's Earth horizon, NOT a cost multiplier, because 1 Jupiter
  system-period = 11.86 yr, not 1 yr. **Time-unit conversions do NOT carry over unchanged** — #535's
  literal `2*pi rad = 1 yr` constant is Sun-Earth-specific (1 nondim time unit = 1/(2*pi) of the
  SUN-JUPITER period here, not the Sun-Earth period); every place #535 converts years to nondim time
  or vice versa must go through Jupiter's own period (`T_JUPITER_YEARS ≈ 11.86`), not the literal
  `2*pi = 1 yr` constant #535 hardcoded. Getting this wrong reproduces the ~10x cost overrun §5
  originally (incorrectly) estimated as unavoidable. Sample density (`t_eval` spacing) must likewise
  be set per REVOLUTION, not per calendar year, or the output arrays inflate ~12x for no benefit.
- `X0_GRID` / `C_VALUES`: **target the neck-open energy band `C ∈ ~[3.00, 3.038]`** (below
  C_L1=3.0388 so the L1/L2 necks are open and Jupiter-region capture is energetically possible), NOT
  the C=3.14 family. `x0` spanning the 3:2 interior resonance semimajor-axis region down through the
  Jupiter-capture bubble. Exact grid set from the anchor's converged (x0, C) once §2a is done.

**Genuinely new (small, well-scoped):**
- Real-object element sourcing + heliocentric-elements→rotating-frame-IC conversion for the anchor
  (§2a). The conversion math is the proven #535 RH120 path; the new part is fetching/citing SBDB
  elements. ~1 short function + a sourced-value comment.
- An early ER3BP sensitivity check for Jupiter — re-parameterize
  `scripts/run_535_er3bp_sensitivity_check.py` with `e=0.0489` and the Jupiter constants (§5). This
  is a clone, not new machinery.
- Nothing in the detector, corrector, or criterion CORE needs building. The criterion note itself
  needs a new companion note (`docs/notes/2026-07-11-557-jupiter-quasi-hilda-admission-criterion.md`)
  recording the §1/§3 re-derivations as the source-of-truth the sweep must match — the same
  "settle-in-writing-first" artifact #535 produced.

**Net: ~90% reuse.** The build is a re-parameterized clone + a small sourced-seed function + a
criterion note. The scientific risk is entirely in the §0/§3 window decision and the §5 fragility
question, not in engineering.

---

## 5. Cost / risk estimate and what must differ from #535

**Compute cost — CORRECTED 2026-07-11 (a Fable review of this plan found the estimate below was
wrong by ~10x, and traced it to a real latent bug, not just pessimism).** CR3BP integration cost
scales with **nondimensional rotating-frame revolutions**, not calendar years. #535's Earth screen
integrated ~50 revolutions (50 yr horizon, 1 rev = 1 yr for Sun-Earth). A Jupiter screen under the
corrected dimensionless window (§1: 10-15 system-periods, ~119-178 yr at Jupiter's 11.86-yr period)
integrates the SAME ~10-60 revolutions — i.e. **the same per-point cost as #535's Earth run**, not
8-12x more. A #535-sized grid should cost about what #535's did (~1 s/point); the abandoned
61x51=3111-point 2D scan analog is **~1-2 CPU-hours**, not the 1-3 hours this section originally
estimated from a years-based (not revolutions-based) accounting.

**The original 8-12x estimate is not just pessimistic — it is the exact symptom of a real bug to
avoid when building this.** #535's own script (`scripts/run_535_earth_transient_quasi_cycler_search.py`)
encodes `t_max = YEARS * 2*pi`, i.e. it hardcodes "1 revolution = 1 year" — true only because
Sun-Earth's period happens to be 1 year. **A naive clone that keeps this constant and simply sets
`YEARS=500` for Jupiter would actually integrate ~500 REVOLUTIONS (≈5,900 Jupiter-years), reproducing
this section's original 8-12x-too-slow estimate as a real, not estimated, overrun.** The Jupiter
build must convert horizon-in-years to horizon-in-revolutions via Jupiter's own orbital period
(`t_max = (YEARS / T_JUPITER_YEARS) * 2*pi`) and set `t_eval` sample density per REVOLUTION, not per
calendar year (arrays otherwise inflate ~12x too). This is a build-time correctness check, not a
performance nicety — get it wrong and the screen silently costs 10x more AND samples the wrong
physical horizon.

**Grid-density compounding — also checked, not expected.** Corridor width in nondimensional phase
space is set by mu, not by window length. Sun-Jupiter's mu is ~318x Earth's, its Hill radius ~0.068
nondim vs Earth's ~0.010, and the capture mechanism itself is Koon 2001's FAST tube-mediated transit
rather than Earth's slow knife-edge horseshoe — all of which argue Jupiter's admissible corridors
should be WIDER, not narrower, permitting a COARSER grid, not a finer/compounding one. No
multiplicative blowup from grid density is expected on top of the (corrected, revolution-based)
per-point cost above.

A `preflight_search` timing pilot (per the #521 discipline, mandatory — the AST ratchet
`tests/scripts/test_scripts_call_preflight.py` enforces the call) must still measure real s/point on
this horizon before committing to grid size — and, as a side benefit, running it will immediately
catch the unit bug above if it's present (it would report ~10x the expected s/point).

**Risk 1 — repeating #535's fragility (idealized-model knife-edge collapses under real dynamics).**
This is the highest-value lesson. #535 found a real admissible Earth corridor that *totally collapsed*
under ER3BP e=0.0167. Mitigation, and a genuine reason for cautious optimism specific to Jupiter:
- **Run the real-eccentricity ER3BP check EARLY — on the ANCHOR, before the broad scan, not as an
  afterthought.** #535 ran its sensitivity check last, after all the width-characterization work; here
  it should gate the broad scan. If the documented-real-object anchor's capture episodes do not
  survive Jupiter's e=0.0489, do not invest in the broad scan.
- **Jupiter's case may be genuinely less fragile than Earth's, for a sourced reason.** Koon 2001
  (p.29) states Jupiter's eccentricity "plays little role during the fast resonance transition"
  because the capture is a *fast* (< 1 Jupiter period) tube-mediated transit, whereas Earth's RH120
  corridor was a *slow* multi-year horseshoe libration acutely sensitive to the annual perturbation.
  Different dynamical mechanism → different fragility profile. This is a hypothesis to TEST early, not
  an assumption — but it means the #535 collapse does not automatically predict a Jupiter collapse.
- Note Jupiter e=0.0489 is ~2.9× Earth's, so the ER3BP perturbation magnitude is larger; the
  fast-transition robustness argument must actually win against that, empirically.

**Risk 2 — repeating #535's instrumentation-gap abandonment (the wide scan died to unbuffered stdout
on a nohup process, never reached a conclusion).** Mitigations, all from the project's own memory:
- **`python3 -u` / explicit `flush=True` from the very first launch** — non-negotiable given #535's
  exact failure here.
- **Incremental checkpoint runlog** (append + flush per grid row), per-row timestamps
  (`date -Iseconds`-style), running admissible-count and ETA — never a black box that only prints at
  the end. The ~1–3 hr runtime makes this essential, not optional.
- **Coarse-first grid**, then refine only around any admissible hit — do not launch the full fine 2D
  grid blind (that is what #535's abandoned scan did).
- **Run synchronously / foreground-with-checkpoint, not detached-and-polled** — the project has
  repeatedly lost detached agents to quota walls; a checkpointed foreground run is recoverable.

**Risk 3 — the structural-window trap (§0/§3).** If the build proceeds under Option B (literal
window), the result is a foregone structural empty and any "0 admissible" MUST NOT be registered in
`data/empty_regions.jsonl` as a dynamical negative (it would be a units artifact — a false negative
by construction, exactly the "verify with a positive control before trusting an all-negative" trap).
Mitigation: **make the Option-A-vs-B window decision explicit and get user sign-off before running.**

**Risk 4 — planar restriction.** The first-pass screen is planar CR3BP; real quasi-Hildas have
nonzero inclination (Gehrels 3 ~1.1°, 111P ~4.2°). A planar screen is defensible for a first pass on
low-i objects (and matches Koon 2001's planar model), but a planar null does not rule out
inclined-capture structure. Record this scope limit; do not over-claim from a planar result.

**What would make this worth doing (net):** the phenomenon is REAL and sourced (documented repeating
temporary captures, a published mechanism), which is a materially stronger starting position than
Earth's single-shot RH120. The dominant risks are a *scope decision* (§0/§3, free to resolve on
paper) and *fragility* (§5 Risk 1, cheaply testable early via the anchor + ER3BP check). The
engineering is ~90% reuse. This is a favorable risk profile FOR A SCREEN — provided the window
decision is made first and the ER3BP anchor check gates the broad scan.

---

## 6. Recommended sequence (if the user approves proceeding — NOT auto-fired)

1. **User decision on §0/§3 window** (Option A rescaled ~120–180 yr, or decline). Blocking.
2. Write the companion criterion note recording the §1/§3 re-derivations (source-of-truth first).
3. Source P/Gehrels 3 elements from JPL SBDB (or fall back to Oterma); convert to the rotating-frame
   anchor IC; reproduce a documented capture episode as the positive control.
4. **ER3BP e=0.0489 sensitivity check on the anchor** (early gate). If it collapses, stop.
5. `preflight_search` timing pilot on the ~400–600 yr horizon; size the coarse grid accordingly.
6. Coarse `(x0, C)` scan over the neck-open band with full `-u`/flush/checkpoint instrumentation;
   refine around hits.
7. Any admissible hit → ER3BP re-check on that exact IC before any catalogue/gauntlet consideration;
   any clean null → register in `empty_regions.jsonl` ONLY if run under the physically-consistent
   (Option A) criterion, never under Option B.

**Implementation model and go/no-go: deferred to the user, per the #557 instruction. This document
recommends nothing about which model builds it.**
