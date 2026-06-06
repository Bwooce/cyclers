# Jones-family corrector variant — DESIGN DRAFT

**Status:** design draft (no code), 2026-06-06. Task #121. Sole write target is
this file. Sibling/parent designs: `2026-06-05-m-ed-realeph-multiarc-discovery-design.md`
(the N-arc corrector this builds on) and the M-3D inclined-backend work
(`core/flyby.py::bend_decompose`, `Ephemeris(model="inclined-circular")`).

**Question.** Why does `ballistic_correct` (`src/cyclerfinder/search/correct.py`)
reach the S1L1 / Sanchez E-M families but **not** the Jones VEM triple-cycler
family, and what corrector variant closes the gap?

**Evidence chain.** `data/OUTSTANDING.md` (#110 density-refuted, #120
3D-inclination-refuted), `docs/notes/2026-06-05-jones-aas17-577-vem-mining.md`
(what Jones actually computed), `src/cyclerfinder/search/correct.py` +
`search/scan.py` + `search/seed_ladder.py` (the current corrector), and
`tests/test_vem_rediscovery.py` (the headline gate).

---

## 0. LEAD FINDING — the headline gate solves ONE repeat period but compares against TWO (a test/data-contract bug, prior to any research)

The IMPORTANT first analysis step asked: does the headline gate compress the
itinerary? It does — **but not in the way the task framing guessed**, and the
mismatch is sharper and more actionable than a "5-leg approximation of a 10-leg
itinerary" story. Two distinct facts, both verified directly against the data:

### 0.1 The period level is NOT compressed — `sequence_canonical` IS one full repeat period

The Jones rows carry a 6-encounter `sequence_canonical` (`E-M-E-V-V-E` outbound,
`M-E-E-V-E-M` inbound) and a 10-segment / 11-encounter `trajectory.segments`
list. The task hypothesised the canonical sequence might be a 5-leg compression
of a 10-leg cycle. **It is not.** Measured (`data/catalogue.yaml`,
`jones-2017-vem-emevve-outbound`):

```
per-segment tof_days: [309, 842, 1668, 1384, 468 | 259, 922, 498, 1939, 1053]
sum of first 5 segments (one E-M-E-V-V-E cycle) = 4671 d = 12.79 yr
sum of all 10 segments                          = 9342 d = 25.58 yr
```

The sourced repeat period is **T = 12.8 yr** (Jones p.9 "the repeat period T is
12.8 years"; catalogue `period.years = 12.8`, `k = 2`). One `sequence_canonical`
cycle (5 legs) already sums to 12.79 yr — i.e. it **is** the full 12.8-yr repeat
period. The 11-encounter table is **two repeat periods** (25.6 yr), exactly as
the table caption states ("over two repeat periods", mining note §1.2). So the
Cell the corrector solves (5 legs, period pinned to 12.8 yr via `period_k=2` on
the VEM beat) is the **correct single-repeat-period model**. Hypothesis (b)'s
"are we solving half the itinerary?" premise is **false at the period level**.

> Caveat worth a plan-time check, not a blocker: the corrector pins
> `period_sec` to `beat(VEM) * k = 6.406 * 2 = 12.81 yr` and then frees `[t0, 4
> leg ToFs]` with the 5th (longest) leg as slack. The sourced single-cycle legs
> are `[309, 842, 1668, 1384, 468]` — a Venus-resonant 1668-d arc dominates, and
> that is the slack leg by `argmax`. The seed is sourced and the period is
> sourced; the geometry is representable. No compression bug here.

### 0.2 The V∞ COMPARISON is mis-contracted — 6 computed values vs 11 sourced values, `strict=True`

This is the real defect. The headline assertion
(`test_vem_rediscovery.py:337-343`) is:

```python
got = sorted(max(|vinf_in|, |vinf_out|) for e in result.best_cycler.encounters)  # 6 values (one cycle)
expected = _sourced_vinf_multiset(entry_id)                                       # 11 values (TWO cycles)
for g, x in zip(got, expected, strict=True):  # ValueError before any V∞ compare
```

A one-repeat-period solve yields **6** per-encounter V∞; the sourced multiset is
the **11**-encounter, two-repeat-period span. `zip(..., strict=True)` raises
`ValueError` on the length mismatch **before a single V∞ value is compared**.
Even a corrector that reproduced the family perfectly would fail this assertion.

The xfail reason text already half-acknowledges this ("the strict-multiset
compare also surfaces a length gap"), but it is buried under the physics
narrative and framed as a *secondary* symptom. It is not secondary: it is a
**test contract bug** sitting on top of the physics finding, and it must be
fixed before the physics result can even be read off this gate.

**Verdict on (b): it is a bug-class finding at the comparison layer, NOT a
research direction at the period layer.** The fix is small and belongs in the
plan as Phase 0:

- Compare **per-repeat-period multisets**: take the sourced 11-encounter list,
  fold it to one cycle (the first 6, or the canonical de-duplicated cycle), and
  compare against the corrector's 6. OR
- Have the corrector emit two cycles (solve the chain, then propagate a second
  cycle) and compare all 11. The first is simpler and matches the single-cycle
  Cell already built.

This does **not** make the gate pass — the physics gap below is real — but it
removes a false failure mode so the gate measures what it claims to.

### 0.3 Why this matters for the design

Once 0.2 is fixed, the gate measures the genuine physics gap (#110/#120): the
closed families floor at 11–18 km/s vs sourced 2.4–7.0, with **zero
bend-feasible**. Everything below is about closing *that* gap. But the design
must lead with 0.2 because a reader running the gate today sees a `ValueError`,
not a 14-km/s V∞ gap, and would mis-diagnose the problem.

---

## 1. Why the corrector reaches S1L1/Sanchez but not Jones — the mechanism

Reading `correct.py` against the Jones method (mining note §2), the gap is **not**
scan density (#110), **not** inclination (#120), and **not** itinerary length
(0.1). It is the **residual definition**. Three nested reasons, in order of
importance:

### 1.1 The residual is |V∞|-magnitude continuity ONLY — it does not match the V∞ *vector*

`_residual_vector` (`correct.py:127-148`) drives, per intermediate encounter,
`|V∞_in| − |V∞_out|` to zero, plus one magnitude closure term. **A flyby
rotates V∞ at constant magnitude**, so magnitude-continuity is *necessary* for
ballistic closure — but it is far from *sufficient*. It says nothing about
whether the required rotation (the bend from `V∞_in` to `V∞_out`) is physically
achievable at that body, nor about the *direction* the next leg needs.

Bend feasibility is checked **post-hoc** (`_bend_feasible`, `correct.py:214-229`)
and is explicitly **never in the residual** (module docstring line 6: "bend
feasibility checked post-hoc, never in the residual"). Consequence: the
least-squares root-find converges to *any* magnitude-continuous chain — and the
nearest such basin from a near-Hohmann seed is a **high-V∞, high-bend, powered**
family where the required turns exceed the gravity-assist maximum. That is
exactly the #110/#120 signature: hundreds of "closed" families, **zero
bend-feasible**, V∞ floored at 11–18 km/s.

**This is the central architectural finding.** The Jones method (mining §2,
pp.5-7) does the opposite: it matches the **full v∞ vector** to within
100–200 m/s (`Δv∞^max`), with the bend realised by a **B-plane-targeted powered
flyby** (turning angle δ = ∠(v∞⁻, v∞⁺), Eq.1; periapsis radius solved from the
required turn, Eq.2; B-plane unit vectors Ŝ/T̂/R̂ and angle θ_B, Eqs.4-5). The
feasibility constraint (altitude window, Eq.7) is *inside* the search, not a
post-hoc filter. Our corrector has the feasibility test (`flyby.py::max_bend`,
`is_ballistic_feasible`, `bend_decompose` all exist) but does not let it *steer*
the solve.

### 1.2 Why S1L1/Sanchez tolerate this but VEM does not

S1L1 (E-M-E-E) and the Sanchez E-M / EEM families are **two-body** (Earth+Mars),
near-resonant, with small required bends and a **wide ballistic-feasible
plateau** (the parent design §0.2 quotes `maintain.py:564` "flat ΔV plateau").
On such a plateau, magnitude-continuity and bend-feasibility nearly coincide:
the nearest magnitude-closed basin *is* bend-feasible, so the cheap residual
happens to land the right family. The S1L1 prototype's own near-anchor hits
(`correct_s1l1_twoarc.py::main`, "near anchors" tag) work for this reason.

VEM is **three-body** with Venus legs. Venus encounters demand large,
specifically-directed turns to bridge the inner-system geometry; the
magnitude-feasible set and the bend-feasible set **diverge**, and the
magnitude-only residual walks into the wrong (powered, high-V∞) basin every time.
The plateau that saved S1L1 does not exist for VEM. This is why #110's dense scan
and #120's inclination both moved the floor by noise (−0.1 to −0.37 km/s): they
re-sampled and re-projected the **same magnitude-only basin structure**. Neither
touched the residual that selects the basin.

### 1.3 Single-ellipse-per-leg is real but secondary

"Single-ellipse-per-leg" (one Lambert conic per leg, `_vinf_nodes`
`correct.py:115-116`) is accurate but is **not** the primary exclusion — the
Jones legs *are* single Lambert arcs between consecutive encounters (mining §2:
"Lambert's problem is solved to determine legs connecting consecutive
encounters"). Jones does **not** chain sub-arcs within a leg; the encounters
*are* the patch points. So a leg = a conic is the right model. What Jones adds is
not more arcs per leg but **vector v∞-matching + B-plane targeting + in-search
feasibility** at the encounters. This reframes hypothesis (a) substantially
(see §3a).

---

## 2. The verdicts on (a)–(d), in brief

| # | Direction | Verdict |
|---|---|---|
| (a) | Multi-arc-per-leg (sub-arcs joined at intermediate unpowered flybys / DSMs) | **REJECT as framed.** Jones legs are single Lambert conics; encounters already are the patch points. The catalogue's 10 segments = two repeat periods, not 10 distinct legs in one cycle (§0.1). There is no missing intra-leg arc to add. *Reframed* version (treat each encounter as a powered-flyby patch with vector matching) collapses into (c). |
| (b) | Full 11-encounter-sequence solving | **BUG, not research (lead finding).** Period level is not compressed (§0.1); the 6-vs-11 mismatch is a `strict=True` zip contract bug in the comparison (§0.2). Fix the compare to per-repeat-period; do not "solve 11 legs". |
| (c) | Eccentric/inclined intermediate flybys with **B-plane targeting + full v∞-vector residual** | **RECOMMEND (primary).** Directly addresses the §1.1 mechanism: replace magnitude-only continuity with vector v∞-matching and bring bend-feasibility into the solve via B-plane targeting — exactly the Jones method. Inclination alone is refuted (#120); the *vector residual* is the lever, and it subsumes (a)-reframed. |
| (d) | N-body shooting refinement (SNOPT-style) seeded by the patched-conic solution | **DEFER (Jones's final step, gated on (c)).** Jones differentially corrects the <200 m/s v∞ mismatch to *ballistic in full n-body* (mining §2, SNOPT). This is the last 200 m/s, and it presupposes a patched-conic solution already in the right family. Useless before (c) lands; aligns with the GMAT stretch goal, out of v1 scope. |

---

## 3. Detailed analysis per direction (trade-offs + evidence)

### (a) Multi-arc-per-leg — REJECT as framed, fold into (c)

- **For:** the project's own S1L1 closure is "multi-arc" (the E-E resonant
  intervals patched by flybys); project memory
  `project_s1l1_realeph_closure_blocker` records S1L1 as genuinely multi-arc.
- **Against (decisive):** the multi-arc-ness is *between* encounters that already
  exist as Cell nodes, not *within* a leg. Jones (mining §2) solves one Lambert
  per consecutive-encounter pair — no deep-space patch points, no intra-leg
  sub-arcs. The catalogue's 10 segments are two 5-leg cycles (§0.1), not a
  10-leg single cycle. Adding sub-arcs would invent structure the source does
  not have and would risk a golden-discipline violation (manufacturing arc
  breakpoints not in Jones).
- **Net:** there is no "missing arc" to add. The thing (a) is reaching for —
  modelling each encounter as a real powered-flyby patch where the bend is
  solved, not assumed — is exactly (c). Reframe and merge.

### (b) Full encounter-sequence solving — the §0 lead finding

- Already settled in §0. The Cell→chain construction does **not** drop
  encounters at the period level; it correctly models one 12.8-yr repeat period.
  The only defect is the 6-vs-11 comparison contract (§0.2). Phase-0 fix.
- One genuine sub-question for the plan: the corrector emits **6** per-encounter
  V∞ for a 6-encounter cycle, but the cycle's two ends (`b0`, `bn`) are the
  *same physical Earth* at the repeat boundary, carrying only one leg each
  (`_per_encounter_vinf`, `correct.py:203-211`). So the de-duplicated cycle has
  **5** distinct flyby encounters, not 6. The per-repeat-period multiset fold in
  §0.2 must account for this boundary identification (don't double-count the
  wrap Earth). This is a comparison-layer detail, still not research.

### (c) B-plane targeting + full v∞-vector residual — RECOMMENDED

- **For (decisive):**
  - It attacks the actual mechanism (§1.1): the magnitude-only residual selects
    the wrong basin. A vector residual `(v∞_out_required − v∞_out_achievable)`
    that *includes the B-plane-targeted bend* selects the ballistic-feasible
    family directly, the way Jones's <200 m/s v∞-match does.
  - The feasibility machinery already exists: `flyby.py::max_bend`,
    `is_ballistic_feasible`, `bend_angle`, `bend_decompose` (the M-3D landed
    work). Today they are post-hoc; (c) moves them *inside* the residual.
  - It subsumes (a)-reframed and is the necessary precursor to (d).
- **Against / risks:**
  - A vector residual over a chain has more ways to be locally infeasible; the
    least-squares may need the B-plane angle θ_B as an explicit free variable per
    flyby (Jones Eq.5), enlarging the free-var vector from `[t0, ToFs]` to
    `[t0, ToFs, {θ_B,i}]`. More DOF = more local minima; the seeding ladder
    (`seed_ladder.py`) becomes load-bearing.
  - Inclination is **refuted as a standalone lever** (#120): (c) must be sold as
    "vector matching + B-plane", not "add inclination". The inclined backend is a
    *fidelity* knob applied *after* the vector residual selects the family, not
    the thing that selects it. State this as an explicit non-goal (§6).
  - Real-eccentricity intermediate flybys: the DE440 control in #120 (full
    3D+eccentric, floor 18.16, zero bend-feasible) shows eccentricity *alone*
    does not help either — consistent with §1.2 (the residual, not the model,
    picks the basin). Eccentric Mars/Venus states come for free from DE440; they
    are not the lever.
- **Gate:** the sourced Jones multisets (Tables 2/3, per-repeat-period fold) within
  `VEM_VINF_TOL_KMS = 0.5`, **with bend-feasible True** (the post-hoc check must
  now pass, since the residual targets the feasible family). The McConaghy
  4.7/5.0 S1L1 cross-check (pending #94) remains the regression anchor for the
  E-M family the corrector already reaches.

### (d) N-body shooting refinement — DEFER

- **For:** it is literally Jones's final step (mining §2: "<200 m/s
  differentially corrected to ballistic in full n-body via SNOPT"). It is what
  makes the patched-conic v∞-match *truly* ballistic.
- **Against:** it operates on the last 100–200 m/s and **presupposes a
  patched-conic solution already in the Jones family** — which only (c)
  produces. Running (d) on today's magnitude-only output would refine a
  *wrong-family* 18 km/s chain to high precision: precise nonsense.
- **Scope:** maps to the GMAT / full-n-body stretch goal (spec §2). Out of v1.
  Sketch it as Phase 3, explicitly gated on (c) converging first.

---

## 4. Recommendation — staged combination, (c)-led

**Recommend (c) as the corrector variant, with (b)'s comparison fix as a
mandatory Phase 0 and (d) as a deferred Phase 3.** Phased sketch:

- **Phase 0 — fix the gate contract (small, do first).** Make
  `test_jones_vem_ballistic_rediscovers_sourced_multiset` compare
  per-repeat-period multisets (fold the sourced 11 → one cycle accounting for the
  wrap-Earth identity, §3b), not `zip(6, 11, strict=True)`. Keep xfail; this only
  removes the false `ValueError` so the gate measures physics. No tolerance
  change, no xfail flip (per the STOP/report discipline already in the file).

- **Phase 1 — vector v∞-residual.** Replace `_residual_vector`'s magnitude terms
  with **full v∞-vector** continuity: at each intermediate encounter, residual =
  `v∞_out_achievable_via_feasible_bend − v∞_out_required_by_next_leg` (3 comps
  per encounter). Use the existing `bend_angle` / `max_bend` to clamp the
  achievable bend; residual penalises the *unachievable remainder*. This alone
  may already pull the solve toward the feasible basin without new free vars —
  cheapest test of the §1.1 hypothesis. Gate: does the EMEVVE floor drop below
  10 km/s and produce *any* bend-feasible closed family? (Compare directly
  against #110's 17.86 / 0-bend-feasible baseline — same grid, new residual.)

- **Phase 2 — B-plane DOF (if Phase 1 underdetermines).** Add per-flyby B-plane
  angle θ_B as a free variable (Jones Eqs.4-5); extend the seeding ladder to seed
  θ_B from the sourced periapsis-altitude column (Tables 2/3 give r_p directly —
  a sourced seed, not a fitted one). Gate: rediscover the sourced Jones multiset
  within 0.5 km/s, bend-feasible True. **This is the xfail-flip gate.**

- **Phase 3 — n-body refinement (DEFER, stretch).** Seed a full-n-body shooting
  refiner (GMAT or a homegrown multiple-shooting TPBVP) from the Phase-2
  patched-conic solution; drive the residual <200 m/s → ballistic. Gated on
  Phase 2 converging. Out of v1.

**What gates the whole effort:** the two sourced Jones multisets (AAS 17-577
Tables 2/3, folded to one repeat period), with bend-feasible required; the
McConaghy 4.7/5.0 S1L1 cross-check (#94) as the don't-break-what-works
regression; and golden discipline throughout (EXPECTED = published Jones values
only; the sourced periapsis altitudes seed θ_B but are never the EXPECTED side).

---

## 5. Honest feasibility risks

1. **Phase 1 may not be enough.** If the magnitude→vector swap still
   underdetermines the bend direction, the wrong basin may persist until B-plane
   DOF (Phase 2). Phase 1 is a cheap discriminator, not a guaranteed fix.
2. **More DOF → more local minima.** A per-flyby θ_B vector makes the
   least-squares landscape rougher; convergence may hinge on the sourced-altitude
   seed being good. If the seed is poor, we are back to scanning — but now over a
   higher-dimensional space (cost risk).
3. **The 0.5 km/s tolerance vs model fidelity.** Jones's own values are
   *n-body-ballistic*; our patched-conic vector match tops out at the 100–200
   m/s Jones himself needed SNOPT to remove. Hitting 0.5 km/s on the patched
   conic may be feasible (it is 2-4× the residual Jones tolerated pre-SNOPT), but
   if not, the honest outcome is xfail-stays + "needs Phase 3", **not** loosening
   the tolerance.
4. **Golden-discipline trap.** Seeding θ_B from sourced r_p is legitimate
   (sourced seed). Deriving r_p from our own Lambert and then "matching" it would
   be circular. The plan must keep the EXPECTED side purely Jones-published.
5. **No new science guaranteed.** (c) is the well-motivated front-runner, but
   #110/#120 are a reminder that confidently-predicted levers have floored at
   noise. The Phase-1 gate is deliberately a fast falsifier.

---

## 6. Explicit non-goals

- **Not** "add inclination" — refuted standalone (#120). Inclination/eccentricity
  are fidelity knobs *after* family selection, not the selector.
- **Not** more scan density — refuted (#110); 2816 points/row already saturates.
- **Not** intra-leg sub-arc chaining (a-as-framed) — Jones legs are single Lambert
  conics; no missing arcs exist (§0.1, §3a).
- **Not** flipping the xfail or loosening `VEM_VINF_TOL_KMS` outside the Phase-2
  gate — the STOP/report discipline in `test_vem_rediscovery.py` stands.
- **Not** full n-body / SNOPT / GMAT in v1 — Phase 3, deferred, stretch-goal.
- **Not** ingesting Table 4 as a cycler — it is a family-mixed architecture table
  (mining §1.4), not a single-cycler signature.

---

## 7. Open questions (verbatim, for the owner / plan author)

1. **Phase-0 fold semantics:** when folding the sourced 11-encounter multiset to
   one repeat period for the comparison, do we (i) take the first 6 encounters,
   (ii) de-duplicate to the 5 distinct flyby encounters accounting for the
   wrap-Earth identity, or (iii) have the corrector emit two cycles and compare
   all 11? (§0.2 / §3b — I recommend (ii), but it changes what "the multiset" means
   and should be an explicit owner decision before the gate is rewritten.)

2. **Does Phase 1 (vector residual, no new DOF) alone move the EMEVVE floor below
   10 km/s and produce any bend-feasible family**, or is the per-flyby θ_B DOF
   (Phase 2) strictly required? (This is the single highest-value experiment; it
   decides whether the variant is a residual swap or a free-var-space expansion.)

3. **Is the patched-conic vector match physically capable of 0.5 km/s on the
   Jones members, or is the 100–200 m/s residual Jones removed with SNOPT an
   irreducible floor** that forces Phase 3 before the xfail can ever flip?

4. **Seed provenance for θ_B:** the sourced periapsis altitudes (Tables 2/3)
   pin r_p, which constrains the bend magnitude but not the B-plane angle θ_B
   directly. Is there enough in AAS 17-577 to seed θ_B without circularity, or
   does θ_B have to be scanned (and if scanned, does that reintroduce the
   density problem #110 in a new dimension)?

5. **Cross-check coupling:** the McConaghy 4.7/5.0 S1L1 anchor (#94) is still
   pending. Should the (c) variant be gated on *also* not regressing the S1L1 /
   Sanchez E-M families it currently reaches, and is the magnitude→vector residual
   swap guaranteed to preserve those (the §1.2 plateau argument says yes, but it
   is untested)?

6. **Does the inclined-circular backend interact with the vector residual** in a
   way the coplanar #120 run could not show — i.e. is there a regime where vector
   matching + inclination *together* help even though inclination alone did not?
   (Low prior, but cheap to check once Phase 1 exists.)

---

## Approval (2026-06-06)

User-approved with all recommendations accepted: (lead) the headline gate's
per-repeat-period comparison bug (6 computed vs 11 sourced V∞, zip strict
ValueError before any comparison) is fixed FIRST as Phase 0; (a) multi-arc-
per-leg rejected as framed; (b) confirmed bug-class, not research; (c)
ADOPTED — full v∞-vector residual + B-plane targeting inside the solve,
staged: Phase 1 vector residual as the cheap falsifier (does any bend-
feasible family appear below 10 km/s?), Phase 2 per-flyby θ_B free variables
seeded from sourced periapsis altitudes as the xfail-flip gate; (d) n-body
refinement deferred to Phase 3 (GMAT-stretch adjacency). §7's six open
questions resolve per the staged falsifier results rather than upfront.
Execution queued behind the Forge Phase 4/5 lane (it holds correct.py live
in the editable install).
