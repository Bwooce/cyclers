# #589: Fable design review of the proposed Gurfil-Kasdin z²-maximization search (GATING, pre-code)

**Date:** 2026-07-14
**Reviewer:** Fable (independent design read, per #581/#583's gating precedent).
**Scope reviewed:** the `#589` entry in `data/OUTSTANDING.md` (3-stage plan: positive-control
reproduction of Gurfil & Kasdin 2002, CMAME 191:2141-2158, DOI 10.1016/S0045-7825(01)00380-2,
via a new z²-maximizing fitness on the existing niching-GA + `characterize()` stack; then a
multi-seed cluster-everything widening; then the #588 dedup/literature/adjudication pipeline).
**Evidence base:** the source paper read IN FULL from
`cyclers_pdf/papers/gurfil-kasdin-2002-out-of-ecliptic-trajectories-deterministic-crowding-...pdf`
(text-extracted, all 18 pages), `scripts/run_581_gurfil_reproduction.py`,
`scripts/run_583_widened_bounded_drift_search.py`,
`src/cyclerfinder/data/validation/er3bp_drift_classifier.py`,
`src/cyclerfinder/core/er3bp_geocentric.py` (constants/signatures),
`docs/notes/2026-07-14-588-gurfil-companion-papers-digest.md`,
`docs/notes/2026-07-14-590-followup-results.md`, plus light verified arithmetic (below).

## Verdict: NO-GO on the 3-stage plan as scoped; SCOPE-DOWN to a half-day analytic note

Stages 2-3 (widened niching harvest + dedup/literature/adjudication) should not be built:
the novelty target is ill-posed by the paper's own construction, the proposed infrastructure
reuse is a category error for this trajectory class, and two factual premises in the #589
entry are wrong against the paper's own text. Stage 1 (the positive control) is, on its own,
a weak gate that mostly re-derives conic arithmetic; it should only be built if a cheap
analytic pre-check (recommended below, ~half a day, no GA) leaves a residual worth chasing.

---

## 1. Two factual corrections to the #589 entry (checked against the paper itself)

### 1a. "a niching GA it never ran" is false — deterministic crowding IS the paper's method

The paper's title is "Characterization and design of out-of-ecliptic trajectories using
**deterministic crowding genetic algorithms**." Section 3 describes the Mahfoud
deterministic-crowding tournament explicitly (their ref [8] — the same niching method
`search/niching_ga.py` implements), run at population 150 x 300 generations (Table 1) on
each of 8 characterization sets + 2 design runs. What the paper genuinely did NOT do:
multi-seed replication, cluster-everything analysis of its final populations (it reports
only the single best-fitness individual per run, Fig. 2 aside), and any literature check.
The residual gap is real but much narrower than the entry claims — and Section 2 below
argues the un-run harvest would not have found additional structure anyway.

### 1b. "search box out to 1-5 AU from Earth" is a misreading of Table 2

Table 2's `(x0)lim` column is the **initial-condition** box, and it is tiny and near-Earth:
per-component limits `[5.7735e-4 AU, 0.3]` interleaved as (position, velocity) per axis.
Verified numerically against the paper's own constraint columns:

- `5.7735e-4 = 1e-3/sqrt(3)` exactly, i.e. `(r0)max = 150,000 km ~ 1e-3 AU` split across
  3 components — the IC position box is **|x0|,|y0|,|z0| <= ~86,400 km**, with the
  hypersphere constraints `6378 km <= r0 <= 150,000 km` enforced by rejection (their
  Eqs. 17-20).
- `0.3 = 15.4 km/s / (29.78 km/s * sqrt(3))` (their normalized velocity unit is Earth's
  mean orbital speed); sets 5-8 halve it to 0.15 (7.7 km/s).

The 1-5 AU figures are `rmax` — a **path inequality constraint** `r(x0,t) <= rmax` enforced
over the whole 5-year window (their Eq. 17), not a search-box dimension. The design cases
narrow the IC box further still: the low-energy run departs a **200-km LEO parking orbit**
(`(r0)max = 6578 km`), the high-energy run departs GEO (42,378 km). So "widen via multi-seed
harvest across the same 1-5 AU search box the paper explored" as written would search a
domain the paper never searched, while missing the actual structure-carrying axis (the
`rmax` constraint sweep — see 2b).

Also model-level: the paper uses the **circular** CR3BP (Earth-centered rotating frame,
`mu = muE/muS`, their Eqs. 7-9) — autonomous, no `theta0`. The existing `SUN_EARTH_ER3BP`
carries `e = 0.0167` and a `theta0` launch phase. Reproduction must run at `e = 0` with no
theta0 gene, or attribution gets muddied for no benefit.

## 2. Q1: is "maximize z²" a coherent niching objective? Mostly no — and the paper's own Fig. 2 says so

The objective is not literally unbounded: `z <= r <= rmax` caps it, so the `rmax` path
constraint (death-penalty rejection) is the real shaping device. But it IS a pure
max-excursion objective, and such objectives generically peak on the **boundary of the IC
box**: more departure energy and a larger allowed excursion always help. The paper's own
convergence analysis (Fig. 2, discussed p. 2147) documents exactly this: the final
population sits at box corners — upper bound on `x0`, lower bounds on `y0`, `x'0` — and the
ONLY surviving "niche structure" is the trivial `+-(z0, z'0)` mirror pair guaranteed by
their symmetry Eq. (14). `y'0` alone converged interior. That is the same
boundary-convergence signature #582/#585 diagnosed as a search-box artifact for a different
objective, here appearing as the **published, expected behavior** of this one.

Consequences:

- A multi-seed niching harvest over this objective would collect mirror copies and
  `rmax`-constraint-active corner variants, not distinct dynamical basins. Deterministic
  crowding preserves diversity only where the landscape is genuinely multimodal; a
  monotone-toward-the-corner landscape leaves it nothing to preserve (compare #583's
  deep-Hill collapse finding: crowding cannot protect lesser niches against a smoothly
  reachable dominant basin).
- The paper's real structure axis is the **`rmax` sweep**: 5 AU -> Type I, <= 2 AU ->
  Type II, ~3 AU -> Type III (their summary, p. 2155). That is a one-parameter constraint
  continuation, not a multimodal search. The coherent "extension the paper never ran" is a
  finer `rmax x (v0)max` grid with frequency-ratio characterization — not a niching harvest.
  Whether even that is worth running is Section 4's question (no).

## 3. Positive control (stage 1) is a weak gate: the paper prints no ICs, and conic arithmetic reproduces its headline numbers

Unlike the family-census paper (Table 3 ICs to 14 digits, Table 4 features — what made
#581's 11/14 reproduction a real control), THIS paper publishes **no optimal initial
conditions at all**. The reproducible targets are two scalars plus qualitative labels:
z_max = 0.223 AU (rmax = 2 AU, v0 <= 12.7 km/s from 200-km LEO, Type II) and z_max = 0.374 AU
(rmax = 3 AU, v0 <= 10.7 km/s from GEO, Type III), heliocentric ranges `0.984 <= q <= 1.124`
and `0.958 <= q <= 2.274`, and the frequency-ratio type.

Checked arithmetic (two-body, no CR3BP needed): departures at these speeds are immediately
hyperbolic geocentric (v_esc = 11.0 km/s at 6578 km, 4.3 km/s at 42,378 km), giving
v_inf = 6.33 and 9.78 km/s. Dumping v_inf into a pure heliocentric plane change
(`sin(i/2) = v_inf/(2*29.78)`) yields i = 12.2 deg -> z_max ~ 0.211 AU (vs published 0.223,
**95%**) and i = 18.9 deg -> 0.324 AU (vs 0.374, **87%**, the remainder consistent with the
eccentricity their own `q <= 2.274` reports — larger heliocentric radius at high ecliptic
latitude). So a GA that walks to the IC-box corner and converts v_inf to inclination
reproduces the paper. A "successful" stage 1 therefore validates almost nothing beyond the
plumbing; a failed stage 1 would almost certainly indicate an infra bug, not a science
finding. As a gate it has little discriminating power in either direction — the mirror of
[[feedback_verify_gauntlet_with_positive_control]]'s "judge by the RIGHT criterion" rule:
there is no sharp criterion available to judge by.

## 4. Q2: "novel" is ill-posed here — three ways

1. **No 4th type exists by construction.** Type I/II/III are the three orderings of one
   frequency ratio: `omega_z < omega_inplane`, `>`, `~=` (their Eqs. 22-24). The taxonomy is
   exhaustive over qualitative orderings; the only unlisted outcome is non-quasi-periodic
   (chaotic) motion, which a z²-max objective does not select for and which would not be
   "a new type" in the paper's sense anyway.
2. **The regime is near-Keplerian, so "new trajectories" are optimization increments, not
   structures.** Per Section 3's arithmetic, ~90-95% of the achievable z is conic
   plane-change arithmetic; the CR3BP correction is a perturbation on a 5-year drift arc.
   There is no family concept, no reproducible-IC census, nothing for the #256-style
   flagger or the `our_status` machinery to attach to. #588's adjudication chain was long
   and only partially conclusive **with** a well-defined family target; here the chain
   would run with no defined success predicate at all.
3. **No catalogue home.** A non-returning, z-maximizing observatory drift trajectory is
   none of the four admitted classes (cycler / quasi_cycler / precursor_mga / mga_tour —
   [[project_catalogue_scope_expanded_2026-06-15]]). Even a "successful" widened harvest
   yields rows `catalogue.yaml` cannot take. The end-state deliverable is a note, which the
   analytic pre-check below produces for ~1% of the compute.

## 5. Q3: the #590 finding and the infrastructure mismatch — it's not a fitness swap

#590 established that "bounded" is horizon-relative even for the PUBLISHED bounded families
(Family J's own Table 3 IC escapes at ~30 yr; Mikkola et al. 2006 predict generic eventual
escape for inclined QS/DRO motion). For #589 this cuts deeper: the target class is not even
nominally bounded — the paper's trajectories **drift away from Earth by design** (their own
text: "the trajectory drifts away from Earth ... a Type II trajectory results", p. 2150),
with a 5-year `t_f` that is an arbitrary mission lifetime, not a dynamical property.
Concretely, every load-bearing piece of the proposed reuse is wrong for this class:

- `escape_radius = 0.5 AU` is a **terminal death event** hard-wired into all three layers:
  `gurfil_kasdin_fitness` (default arg, fitness 0 on escape), `characterize()` (hard-coded
  `0.25 = 0.5²` event), and `classify_bounded_drift` (`ESCAPE_RADIUS_NORM_DEFAULT`). Every
  single paper trajectory reaches 1-5 AU from Earth and would be killed/scored divergent
  on sight. The existing stack does not merely mis-rank this class — it cannot represent it.
- `characterize()`'s outputs (1yr/5yr rmin/rmax annulus, `practically_stable <= 0.1 AU`,
  DRO/DPO/ERO/DEO typing by velocity sign) and the `TABLE34` 14-family match criteria
  (`IC_PROXIMITY_TOL`/`FEATURE_FACTOR_TOL`) are all bounded-geocentric-annulus concepts,
  meaningless at 2-5 AU drift scales.
- The #586 cluster-everything harvest gates clusters on `classify_bounded_drift` verdicts —
  for a class where "divergent" is the intended behavior, that gate inverts into a
  reject-everything filter (or, if bypassed, filters nothing).
- What the paper's class actually needs is a **new characterization layer**: FFT/peak
  estimation of `omega_z` vs `omega_x, omega_y` (Type I/II/III label), z_max, endurance
  above a height threshold, heliocentric `q` range, and the `r(t) <= rmax` path-constraint
  handling inside the fitness — plus a `e = 0` circular-CR3BP configuration with no theta0.
  That is a #581-stage-2-sized capability build, not "a new fitness function on existing
  infrastructure." The #589 entry's own "not yet designed" caveat is doing a lot of work.
- Cost note (CPU-constrained box): each fitness eval is a 5-year integration spanning
  6578-km perigees to 5 AU (stiff, 4 orders of magnitude in r) — substantially heavier than
  the current 1e6-1e7 km annulus evals — times pop 150 x 300 gens x sets x seeds.

## 6. Q4/Q5: recommendation

**Do not build stages 2-3.** The strategic picture from #588/#590 — enormous adjudication
effort on a WELL-defined family space ending at "narrowed, not resolved" — argues against
opening a strictly worse-defined space. This is the shinier-but-shakier target the review
brief suspected.

**Do this instead (half a day, no GA, CPU-negligible):**

1. Commit an analytic characterization note: the two-body v_inf/plane-change/eccentricity
   accounting of Section 3 (expanded to all 8 characterization sets if desired), the
   taxonomy-exhaustiveness argument of Section 4, and the Table 2 encoding derivation of
   Section 1b. This IS the reproduction, at the fidelity the paper's published data
   supports (scalars + types, no ICs printed).
2. Correct the #589 entry's two factual premises (done alongside this review).
3. File the regime in the negative-results registry sense: empty of catalogue-relevant
   structure **by construction** (non-returning class, exhaustive taxonomy, near-Keplerian
   dynamics), re-open condition = a reframed objective that ties z-excursion to a
   RETURNING/bounded class.

**If a reproduction is still wanted afterwards**, gate it on the analytic note leaving a
genuinely unexplained residual (it currently leaves ~5-13% on two scalars, comfortably
attributable to eccentricity + CR3BP perturbation), and spec it as its own small
circular-CR3BP build with the Section 5 characterization layer — judged by "matches the
paper's scalars/types," never by the bounded-drift stack.

**Salvageable reframing worth a separate task (flagged, not scoped here):** "maximize
out-of-ecliptic excursion SUBJECT TO remaining bounded/returning" — e.g. a max-|z| /
max-inclination figure-of-merit sweep across the already-validated bounded DRO/QS families
from #581-#586. That keeps the observatory-utility idea, has a defined success predicate,
reuses the existing infrastructure honestly (bounded class, existing classifier, existing
family census as anchors), and produces catalogue-relevant annotations. That is a different
task than #589 as written.
