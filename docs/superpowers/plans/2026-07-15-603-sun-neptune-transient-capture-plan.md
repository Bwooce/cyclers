# #603 Sun-Neptune Transient-Capture `quasi_cycler` Screen — Design Plan (review-first)

> **Status: PLANNING ONLY.** This is a design document for review, mirroring how #557's
> plan (`docs/superpowers/plans/2026-07-11-557-jupiter-quasi-hilda-transient-capture-plan.md`)
> was written and reviewed before any implementation code was dispatched. **No sweep script has
> been written, nothing has been run, and no `data/catalogue.yaml` / `data/empty_regions.jsonl`
> edits have been made.** Implementation model and go/no-go are for the user to decide AFTER
> reviewing this document, per the explicit #603 OUTSTANDING entry ("plan first, user decision
> on scope after"). Everything numeric below was computed in-repo or verified by live web search
> this session, not recalled from training data — every claim below states its verification
> method.

**Goal:** Decide whether — and if so, exactly how — to point #535/#557's validated Hill-sphere-
return `quasi_cycler` screen (`src/cyclerfinder/search/hill_sphere_return_detector.py`) at
Sun-Neptune, extending the machinery proven twice (Earth, Jupiter) to a third system.

**Headline recommendation (read this first): DECLINE to build a full discovery sweep under the
naive extended recipe. If any Neptune result is wanted at all, the only defensible option is a
small, cheap, explicitly-caveated idealized-CR3BP-only characterization with NO positive control
and NO catalogue path — not a real discovery search.** Two independent findings, either one of
which alone would be a serious yellow flag, stack here:

1. **The window-rescaling problem is not "bigger," it is qualitatively different.** Neptune's own
   orbital period is **164.89 yr** (computed below), so the *same* dimensionless logic #557 used
   for Jupiter (floor = 1 system period, window = 10-15 system periods) gives a floor of **~165
   yr** and a window of **~1649-2473 YEARS**. That is not merely a bigger number than Jupiter's
   already-flagged 119-178 yr departure from the catalogue's 10-15 yr mission-relevance
   `quasi_cycler` definition — it is ~13.9x further still, deep into territory (a sliding window
   longer than most of recorded human history) where "admissible under this window" no longer
   plausibly means anything a spacecraft mission, or arguably any bounded scientific claim about a
   *repeating cycler*, could use. This is flagged honestly below (§2) as a harder class-definition
   problem than Jupiter's, not solved by more careful arithmetic.
2. **No documented real anchor object exists.** Unlike Jupiter (Koon et al. 2001's two worked
   examples, Oterma and Gehrels 3, both real comets with *documented, repeating, Hill-sphere-
   entering* temporary Jovian captures), a live literature search (§3 below) found **no published
   case of any real minor body entering Neptune's actual Hill sphere** in a temporary-capture
   episode analogous to a quasi-Hilda comet. The closest real, documented Neptune-associated
   transient objects — (309239) 2007 RW10 and the "four temporary Neptune co-orbitals" — are
   **quasi-satellite / 1:1-resonance co-orbitals, a structurally different mechanism that by
   definition stays OUTSIDE the Hill sphere**, and Horner & Evans (2006) explicitly document that
   Neptune (like Uranus) has "great difficulty" capturing Centaurs into even this weaker 1:1
   coupling, compared to Jupiter/Saturn. There is no fallback anchor here the way Oterma was a
   fallback for Gehrels 3 at Jupiter — the fallback tier itself is empty.

---

## 0. What was checked before writing (grounding)

- **#557's plan and criterion note read in full**
  (`docs/superpowers/plans/2026-07-11-557-jupiter-quasi-hilda-transient-capture-plan.md`,
  `docs/notes/2026-07-11-557-jupiter-quasi-hilda-admission-criterion.md`) — this document follows
  their section structure and terminology directly; every place Neptune's numbers depart from a
  straight rescale of Jupiter's is flagged explicitly rather than silently substituted.
- **#535's original criterion note read in full**
  (`docs/notes/2026-07-03-535-quasi-cycler-transient-drift-admission-criterion.md`) — confirms the
  1-yr Earth floor was always understood as "one primary-secondary orbital period," explicitly
  flagged in its own §6 as needing re-derivation for any other system, and that the floor's job is
  merging intra-episode grazing, never encoding an assumed recurrence period.
- **Detector re-read** (`src/cyclerfinder/search/hill_sphere_return_detector.py`) — confirmed
  (again) it is pure numpy, takes `r_hill`, `min_separation`, `window_lo/hi`, `n_returns_lo/hi`,
  `geometry_factor` all as parameters, nothing system-specific hardcoded. No code changes needed
  for Neptune specifically; the same caveat #557 flagged about time-unit conversions in the
  *calling script* (not the detector) applies identically here (§5).
- **#557's actual run artifacts read directly, not summarized from memory**
  (`scripts/run_557_jupiter_quasi_hilda_search.py`, `scripts/run_557_er3bp_sensitivity_check.py`
  both exist; the `sun-jupiter-quasi-hilda-transient-capture-quasi-cycler` entry in
  `data/empty_regions.jsonl` was read in full). Confirmed #557 actually ran (query in the task
  prompt was correct): 264-point idealized CR3BP scan found 16/264 admissible points (6%,
  genuinely richer than Earth's single corridor), but the ER3BP real-eccentricity gate
  (e=0.04838624) killed 15/16 outright and characterized the 1 apparent survivor as a knife-edge
  coincidence (collapses under any of e/x0/C perturbed by as little as 0.005) — a clean,
  well-characterized negative, exactly as the task brief described. Cost was **169.5 s wall for the
  full 264-point sweep** — this actual, cheap number is the load-bearing fact for the cost estimate
  in §6 below (Neptune should cost about the same, not more, per the same revolutions-not-years
  argument #557's own §5 corrected).
- **`data/empty_regions.jsonl` searched directly for `neptune`** (not trusting the task prompt's
  characterization): four entries found —
  `repeated-moon-neptune-sweep` (#254/#253, patched-conic Lambert repeated-moon-flyby search),
  `uranus-neptune-regular-moon-endgame-vilm-2026-06-23` (#465, leveraging-VILM moon-tour endgame),
  `neptune-proteus-triton-proteus-multirev-leveraging-2026-06-26` (#465, powered multi-rev
  leveraging releg), and `neptune-triton-proteus-symmetric-closure-599-2026-07-15` (#599, direct
  symmetric-closure enumeration). **Confirmed directly by reading each entry's `method_capability`
  and `centre` fields: every one of these is Neptune-CENTRIC** (spacecraft/moon-tour trajectories
  around Neptune, using Triton/Proteus as flyby bodies, patched-conic or Lambert-releg genomes).
  **None carries the `hill-sphere-return-detection` capability tag** #535/#557's method uses, and
  none is a *heliocentric* small-body screen — the question here (does a Sun-orbiting minor body
  get transiently captured by Neptune's gravity) is genuinely unswept by any existing entry. This
  directly confirms the task brief's own claim rather than merely repeating it.
- **`core/constants.py` read for Neptune's sourced elements** — `_NEPTUNE_SMA_AU = 30.06992276`
  (Standish & Williams), `PLANETS["N"].ecc = 0.00859048` (J2000 mean, Standish & Williams Table 1,
  the same in-repo sourced table Jupiter's e=0.04838624 came from).
- **All Neptune-specific physical constants computed in-repo this session** (not recalled — see
  §1/§4 for the exact values and how each was derived/cross-checked).
- **Live web search performed for a real Neptune temporary-capture anchor** (§3) — this is the one
  piece of #557's grounding that could not simply be re-run in-repo; searched directly rather than
  relying on training-data recall of an "exactly right" obscure astrodynamics fact, per this
  project's own standing discipline.

---

## 1. Verified physical constants (computed in-repo, not recalled)

| Quantity | Value | How computed / cross-check |
|---|---|---|
| Sun-Neptune mass ratio μ | **5.151118e-5** | `cr3bp.cr3bp_system("Sun","Neptune").mu` (DE440 system GM, Park et al. 2021, same source class as Jupiter's) |
| Neptune orbital period T_N | **164.8947 yr** | `360 / PLANETS["N"].mean_motion_deg_day / 365.25`, from the in-repo sourced SMA (30.06992276 AU). Cross-checked independently via plain Kepler's third law, `T ≈ a^1.5` (Sun-dominated, Neptune's mass negligible): `30.06992276^1.5 ≈ 164.89 yr` — matches to 4 sig figs. Also consistent with the widely-published ~164.8 yr sidereal period. |
| Neptune Hill radius R_hill | **0.0257984 nondim = 0.775757 AU** | `(mu/3)**(1/3)`, same formula as Earth/Jupiter — scale-invariant, transfers with no code change |
| Neptune eccentricity e | **0.00859048** | `PLANETS["N"].ecc`, in-repo sourced (Standish & Williams Table 1 J2000 mean) — **~5.6x smaller** than Jupiter's 0.04838624 |
| Collinear libration Jacobi constants | **C_L1 = 3.0058187, C_L2 = 3.0057500** | Computed this session: `brentq` root of `dΩ/dx` at the sourced μ, then `cr3bp.jacobi_constant()` (the project's own function) at that root — verified against Jupiter's known values (3.038761/3.037489) using the identical procedure before trusting the Neptune numbers, so the method itself is checked, not just the Neptune output. |
| T_Jupiter/T_Earth (for comparison) | 11.868 | From #557's own note |
| **T_Neptune/T_Earth** | **164.8947** | Directly — Earth's own period is 1 yr by definition |
| **T_Neptune/T_Jupiter** | **13.90** | `164.8947 / 11.868` — Neptune's system-period-relative window is ~13.9x further from the catalogue's literal 10-15 yr figure than Jupiter's already was |

---

## 2. The window/floor re-derivation for Neptune (the crux)

**Applying #557's Option-A logic literally: floor = 1 Neptune period ≈ 165 yr; window =
10-15 Neptune periods ≈ 1649-2473 yr.**

This is *not* a new derivation choice — it is the literal, unmodified application of the same
dimensionless statement #557 settled on ("floor = 1 system period, window = 10-15 system periods,
count 3-15 returns, geometry ratio 3.0 — a statement about ROTATING-FRAME REVOLUTIONS, inherently
system-agnostic because CR3BP nondimensional time IS defined in units of the system's own period").
In that sense the *method of derivation* transfers perfectly and there is nothing new to invent
here — which is exactly why the result is worth taking seriously rather than dismissing as an
arithmetic accident: **the criterion was always dimensionless; Jupiter's 119-178 yr number only
"looked like" a modest, defensible departure from the catalogue's 10-15 yr figure because
Jupiter's own period (11.9 yr) happens to be of the same order as a human career.** Neptune's
period (164.9 yr) is not, and the same formula pushed through honestly gives 1649-2473 yr.

**Did I look for a different, better-motivated floor/window specific to Neptune, rather than just
scaling Jupiter's approach up? Yes — and I did not find one.** Candidates considered and rejected:

- **A local-resonance libration period** (the Jupiter-Hilda 3:2 libration period, ~250-300 yr, was
  rejected by #557 for presupposing the object is *in* a specific resonance). The analogous
  Neptune exterior resonances relevant to scattered-disk/Centaur dynamics (e.g. 2:1, 5:2 with
  Neptune) have libration/superperiod timescales that are *also* set by multiples of Neptune's own
  ~165 yr period (few hundred to few thousand years) — so this does not give a shorter, more
  tractable number; if anything it is comparable-to-longer than the naive system-period floor.
  Rejected for the same reason #557 rejected it: it presupposes the resonant structure the screen
  is trying to discover.
- **The real, documented recurrence spacing of an anchor object** (Jupiter used Gehrels 3's
  ~1970-1973 capture / multi-decade recapture prediction as an *informational* cross-check, not the
  floor itself, but it was reassuring corroboration that the floor sat just below real return
  spacing). **No such cross-check is available for Neptune** — §3 found no documented Hill-sphere-
  entering real capture episode to check the floor against at all. This is not a failure of this
  plan's derivation; it is a genuine absence of grounding data that #557 had and this system does
  not.
- **A co-orbital/quasi-satellite libration timescale** (2007 RW10's ~7,500-12,500 yr libration
  period, §3): if anything this argues the *real* Neptune-coupled transient dynamics run on
  timescales *longer* than the naive 165-yr-period floor, not shorter — another data point against
  finding a smaller, more tractable number by looking harder, not for one.

**Conclusion: 1649-2473 yr is the best-motivated number available, and it is a genuinely different
kind of problem than Jupiter's 119-178 yr, not just a bigger one.** At Jupiter, the departure from
the catalogue's stated mission-relevance window was defensible as "the same *shape* of criterion,
scaled to the system" (#557's own framing) while still describing something a human civilization
could plausibly observe recur (a comet's capture cycle within centuries). At Neptune, a "recurring"
structure needing ~2000+ years between windows describes something on the order of the length of
recorded civilization itself — it is fair to ask whether a candidate "admissible" at this window
is testing the same scientific question `quasi_cycler` was meant to capture at all, or is really
answering a different question (long-timescale resonant/chaotic transport structure, which the
project already has other methods for — e.g. the #562/#563 symmetric-closure and resonance-hopping
machinery) wearing the `quasi_cycler` label. This is flagged as an open class-definition question
for the user, not resolved here, exactly as #557 flagged (and left open) its own smaller version of
this question.

**A second, genuinely new (not present at Jupiter) rescaling requirement, found while checking the
grid design, not assumed:** the neck-open Jacobi-constant band (`C < C_L1`, where capture through
the L1/L2 necks is energetically possible) scales with μ. Jupiter's #557 scan covered `C ∈
[3.000, 3.038761]`, a band **0.0388 wide**; Neptune's analogous band is `C ∈ [~3.000, 3.0058187]`,
only **~0.0058 wide — about 6.7x narrower** (consistent with the libration-point Jacobi-constant
offset scaling as ~μ^(2/3): `(mu_Jupiter/mu_Neptune)^(2/3) ≈ 6.98`, close to the directly-measured
ratio). **Practical consequence: reusing #557's Jacobi-constant grid step (0.005) at Neptune would
leave room for only ~1 usable grid point below C_L1 before the neck is already closed** — the
energy-axis grid, not just the time axis, needs re-derivation, roughly ~7x finer
(e.g. step ≈ 0.0007-0.001) to resolve a comparably-sampled neck-open band. This is a genuinely new
finding this plan surfaces that #557 did not have to confront (Jupiter's larger μ happened to give
it a comfortably wide neck-open band).

---

## 3. Does a real anchor object exist? (live web search performed, not recalled)

**Searched this session** (queries and sources below); **conclusion: no.** No published case of a
real minor body entering Neptune's actual Hill sphere in a temporary-capture episode analogous to
Koon et al. (2001)'s Jupiter-family-comet mechanism was found.

**What was found instead, and why none of it substitutes:**

- **(309239) 2007 RW10** — de la Fuente Marcos & de la Fuente Marcos (2012), *A&A* 545, L9
  ("(309239) 2007 RW10: a large temporary quasi-satellite of Neptune," arXiv:1209.1577): a real,
  ~250 km object, currently in a **temporary quasi-satellite state** (has been for ~12,500 yr, will
  remain for another ~12,500 yr, previously an L5 Trojan, will become an L4 Trojan next), with
  mean-longitude libration period ~7,500 yr. **This is a 1:1 mean-motion-resonance (co-orbital)
  phenomenon.** By construction, a quasi-satellite orbits the Sun in resonance with the planet
  while staying well OUTSIDE the planet's Hill sphere — it is defined by never being gravitationally
  bound to the planet, which is the exact opposite of the `< r_hill` encounter this detector looks
  for. Propagating 2007 RW10's real sourced elements would almost certainly show **zero** Hill-
  sphere encounters — the same "wrong dynamical object, guaranteed-zero-encounters-by-construction"
  problem #557 flagged for the C=3.14 certified-periodic Hilda seed (§2b of the Jupiter plan). It
  cannot serve as a positive control for THIS detector.
- **The "four temporary Neptune co-orbitals"** — de la Fuente Marcos & de la Fuente Marcos,
  "Four temporary Neptune co-orbitals: (148975) 2001 XA255, (310071) 2010 KR59, (316179) 2010 EN65,
  and 2012 GX17" (arXiv:1210.3466): same structural issue — documented temporary 1:1-resonance
  co-orbitals (horseshoe/quasi-satellite transitions), not Hill-sphere-entering captures.
- **Horner & Evans (2006)**, "The Capture of Centaurs as Trojans" (arXiv:astro-ph/0511791;
  *MNRAS Letters* 367, L20) — the directly relevant negative finding: **"Uranus and Neptune seem to
  have great difficulty capturing Centaurs into the 1:1 resonance, while Jupiter captures some, and
  Saturn the most (~80%)."** Even the *weaker*, Hill-sphere-avoiding co-orbital-capture phenomenon
  is documented as rare at Neptune specifically, relative to Jupiter/Saturn. This is independent
  literature evidence — not an assumption — that whatever underlying dynamical coupling drives
  temporary-capture phenomena generally is documented as weaker at Neptune, not merely
  under-studied.
- **General Centaur-giant-planet close-encounter literature** (2060 Chiron's dynamical history,
  Centaur "severity of close encounter" scales, Nesvorný/Vokrouhlický irregular-satellite-capture
  papers) — these document real close encounters and even permanent captures (the irregular
  moons), but as **formation-epoch** (Nice-model-era, ~4 Gyr ago) one-off binary-disruption capture
  events, not present-day *repeating* transient captures of the quasi-Hilda kind. No specific,
  named real object with a documented **present-day, Hill-sphere-entering, repeating** temporary
  Neptune capture was found in any of these.
- **Tisserand-parameter-with-respect-to-Neptune literature** (used to classify Centaur injection
  pathways, e.g. arXiv:2511.03021) discusses population-level dynamics, not a specific documented
  repeating-capture object either.

**Honest bottom line for §3: this is a materially weaker starting position than Jupiter's, and
weaker than even Earth's single-shot RH120 case (which was at least a genuine, if single-shot,
real Hill-sphere-entering minimoon capture).** Neptune has no known real object that has ever been
documented entering its Hill sphere in a capture episode of the relevant kind, repeating or not.

---

## 4. Seeding strategy (given §3's finding)

**No real-object anchor is available.** The #535/#557 "positive control first" structure (2a: real
anchor; 2b: broad scan around it) cannot be followed as written, because step 2a has no candidate
to source elements for. Two honest paths:

- **(Recommended if anything proceeds) Idealized-model-only scan, explicitly flagged as having NO
  positive control.** Scan `(x0, C)` directly across the neck-open band (§2, with the finer C-grid
  §2 identifies) using the same detector, same rescaled window/floor, same ER3BP survival gate —
  but report any result (admissible or not) as **unvalidated by any known real object**, exactly
  the caveat #557's own document said it would need to attach if its Gehrels-3-or-Oterma fallback
  had also failed to materialize (it did not, at Jupiter — but the contingency language exists in
  that plan precisely for a case like this one). Any "0 admissible" result under this path is
  weaker evidence than #557's clean, control-validated negative and **must not** be registered in
  `data/empty_regions.jsonl` with the same confidence tier as #557's entry; if registered at all,
  it should be explicitly labeled as method-unvalidated (no positive control available).
- **(Not recommended, but recorded for completeness) Construct a synthetic "closest known
  candidate" anchor from 2007 RW10 or one of the four temporary co-orbitals anyway**, purely to
  exercise the pipeline mechanically (elements → IC → propagation). This would almost certainly
  return zero Hill-sphere encounters (§3), which tells the builder nothing about Neptune's actual
  capture dynamics — it would only validate that the elements→IC conversion code runs, not that the
  screen's premise is sound. Not worth building for that alone.

**Neither path produces a genuine positive control.** This is the second independent reason (after
§2's window problem) the headline recommendation is to decline, not merely to proceed cautiously.

---

## 5. Threshold-adaptation table (scale-invariant vs. timescale-dependent)

| Parameter | Earth (#535) | Jupiter (#557) | Neptune (this plan) | Verdict |
|---|---|---|---|---|
| Encounter distance = `r_hill` | 0.01 nondim | 0.0682534 nondim | **0.0257984 nondim = 0.775757 AU** | **Scale-invariant** — `(mu/3)^(1/3)`, transfers with no code change |
| Geometry factor | 3.0 | 3.0 | 3.0 | **Scale-invariant** (dimensionless ratio) |
| Return count `n_returns` | 3-15 | 3-15 | 3-15 | **Dimensionless**, transfers as-is |
| Separation floor `min_separation` | 1.0 yr | 11.868 yr | **164.8947 yr** | **Timescale-dependent** — rescale to Neptune's own orbital period (§2) |
| Admission window | 10-15 yr | 119-178 yr | **1649-2473 yr** | **Timescale-dependent and the crux** (§2) — an order of magnitude further from the catalogue's stated 10-15 yr figure than Jupiter's departure |
| Propagation horizon | 50 yr (~50 rev) | 653 yr (55 rev) | **~9069 yr (55 rev)** | Same **revolution count**, not a cost multiplier — see §6 |
| Neck-open Jacobi-constant band width (`C_L1 - 3.00`) | n/a (no libration-point neck concept used) | 0.0388 | **0.0058** (~6.7x narrower) | **Genuinely new for this system** — grid step must be rescaled, not just window/floor (§2) |
| Real positive-control anchor | RH120 (single-shot minimoon) | Gehrels 3 (repeating, documented) | **None found** (§3) | **Missing entirely** — the one row with no transfer path at all |

---

## 6. Build vs. reuse

**Reuse verbatim (zero changes), same as #557's own assessment:**
- `src/cyclerfinder/search/hill_sphere_return_detector.py` — fully system-agnostic, nothing to
  change.
- `core/cr3bp.py`, `core/er3bp.py`, `preflight_search`/`MethodCapability` plumbing, the
  empty-region registry round-trip.

**Reuse with re-parameterization only, IF built (structurally a clone of
`scripts/run_557_jupiter_quasi_hilda_search.py`):**
- `system = cr3bp.cr3bp_system("Sun", "Neptune")` → μ and `r_hill` recompute automatically.
- Constants: `MIN_SEPARATION_YEARS = 164.8947`; `WINDOW_LO/HI_YEARS = 1648.9/2473.4`;
  `N_RETURNS_LO/HI = 3/15`; `GEOMETRY_FACTOR = 3.0` — all unchanged in *kind*, just recomputed.
- Propagation horizon `≈ 55 periods ≈ 9069 yr` (mirroring #557's own 4x-window-midpoint choice) —
  **must go through Neptune's own period, never a hardcoded `2*pi = 1 yr` constant** (the exact
  latent bug #557's Fable review caught; the fix that script already has must be preserved, not
  reintroduced as a fresh mistake for a third system).
- `X0_GRID`/`C_VALUES`: target `C ∈ [~3.000, 3.0058187]`, **with a step roughly 7x finer than
  #557's 0.005** (§2/§5) to resolve the much narrower neck-open band — this is new, not a copy.

**Genuinely new, and this is where Neptune differs from #557's "~90% reuse" verdict:**
- No anchor-sourcing function is buildable at all (§3/§4) — the ~1-function piece #557 needed
  (SBDB fetch + vis-viva conversion) has nothing to attach to here. This is a *subtraction* from
  #557's build list, not an addition, but it removes the pipeline's only validation step.
- The C-grid-density rescaling (§2) is a real, non-trivial design question #557 never had to
  answer, because Jupiter's larger μ gave it a comfortably wide neck-open band by luck, not by any
  general property of the method.
- A criterion note, if this proceeds, would need its own explicit "NO POSITIVE CONTROL AVAILABLE"
  section — a structurally different document from #535's and #557's, which both opened with a
  validated anchor.

**Net: reuse fraction is high on the mechanical/code side (~85-90%, similar to #557), but the
*evidentiary* structure is fundamentally thinner — the one thing that cannot be reused from #535/
#557 is the part that made either of their results trustworthy (a validated positive control).**

---

## 7. Cost/risk estimate

**Compute cost — apply #557's own hard-won lesson, stated explicitly for Neptune.** CR3BP
integration cost scales with **nondimensional rotating-frame revolutions**, not calendar years.
#557's actual measured cost for a 264-point, 55-period-horizon scan was **169.5 s wall** (read
directly from its `empty_regions.jsonl` entry, not estimated). A Neptune scan at the same
revolution count (55 periods = ~9069 yr, vs Jupiter's 55 periods = ~653 yr) should cost **about the
same, not ~14x more** — the ~9069-year number sounds enormous but is the *same* 55 rotating-frame
revolutions as Jupiter's run. **A naive clone that keeps any hardcoded `2*pi = 1 yr` conversion and
sets `YEARS=9069` directly would integrate ~9069 REAL revolutions (over a million Neptune-years) —
reproducing the exact ~10x-and-worse unit bug #557's own Fable review caught and fixed once
already.** This is the single most important engineering trap to avoid, restated for the third
time now across Earth/Jupiter/Neptune, and worth flagging precisely because "surely nobody would
make the same mistake twice" is exactly the kind of complacency that causes a mistake to be made a
third time.

**Risk 1 — no positive control means any result (positive or negative) is unvalidated.** This is
the dominant risk, already covered in §3/§4. Any admissible idealized hit would need to be treated
as *purely* a mathematical curiosity of the model (no evidence any real object exhibits it); any
"0 admissible" result would be a weaker, unvalidated negative, not a clean one like #557's.

**Risk 2 — the eccentricity-fragility question is genuinely untested, and the "smaller e might
survive" intuition in the #603 OUTSTANDING entry does not hold up under the project's own prior
evidence.** The #603 OUTSTANDING motivation notes Neptune's e=0.00859 is ~5.7x smaller than
Jupiter's e=0.04838624, suggesting the ER3BP collapse mode "may not fire." **Checked against the
project's own track record, not assumed: Earth's own eccentricity (e=0.0167) is itself ~2.9x
smaller than Jupiter's, yet Earth's corridor still totally collapsed under its own real e (#535).**
Two prior data points spanning e=0.0167 to e=0.0484 both collapsed; there is no established
monotonic "smaller e survives better" pattern in this project's own evidence to extrapolate from.
Neptune's much-smaller e is a plausible reason for optimism but is **not** independently supported
by the two existing results — it would need its own empirical test (were there anything worth
testing it on).

**Risk 3 — grid-density compounding on the energy axis, not just cost.** §2/§5's finding that the
neck-open C-band is ~7x narrower than Jupiter's is a genuine, previously-unconsidered complication:
a grid that isn't rescaled correctly on this axis would not just be "coarser than ideal" but could
plausibly miss the entire physically-relevant band (or, if the band is missed by using #557's
literal step, silently scan a region where the necks are already closed, C > C_L1, near-certainly
returning a guaranteed-empty structural result masquerading as a dynamical one — precisely #557's
own §3 Option-B trap, in the energy dimension instead of the time dimension).

**Risk 4 — the class-definition question (§2) is not really a "risk" in the technical sense; it is
the real reason to decline.** Even a perfectly-executed, correctly-costed, correctly-gridded sweep
under the rescaled window answers a question ("does an idealized point exist that's admissible over
a 1649-2473-yr window") whose relevance to the catalogue's `quasi_cycler` class is itself in serious
doubt at this timescale, independent of whether the sweep finds anything.

---

## 8. Recommended sequence (if the user overrides the headline recommendation and wants a minimal
run for completeness)

1. **User decision on whether the ~1649-2473 yr window question is even worth answering** (§2) —
   blocking, and this decision is harder than #557's because there is no obviously-acceptable
   "Option A" here the way there was at Jupiter.
2. If proceeding: write a criterion note explicitly headed "NO POSITIVE CONTROL AVAILABLE" (§3/§4),
   record the §2/§5 C-grid rescaling derivation before any code is written.
3. `preflight_search` timing pilot on the ~9069-yr/55-period horizon, sized against #557's actual
   169.5 s/264-point baseline — this pilot will also immediately catch the unit bug in §7 if
   present (it would report a wildly larger s/point than the Jupiter baseline).
4. Coarse `(x0, C)` scan over the rescaled neck-open band with the finer C-grid (§2/§5), full
   `-u`/flush/checkpoint instrumentation (per this project's own #535 instrumentation-loss lesson).
5. **ER3BP e=0.00859048 sensitivity check on any idealized hit** — with the §7 Risk 2 caveat that
   there is no principled reason yet to expect it to fare better than Earth's or Jupiter's did.
6. Any result (hit or clean) gets registered, if at all, with an explicit
   `positive_control: none available` field distinguishing it from #535's and #557's
   control-validated entries — never silently presented at the same evidence tier.

**Implementation model and go/no-go: deferred to the user, per the #603 instruction. This document
recommends DECLINING outright as the primary path, with the above as the explicitly-weaker fallback
if the user wants a documented characterization anyway.**
