# 3D-inclined "releg" genome — re-open the INCLINED moon-tour dead regions (DESIGN DRAFT, #451)

**Date:** 2026-06-26
**Status:** DESIGN-DRAFT — for user review. No production code written by this doc.
**Issue:** #451. Ranked **capability-lever #5** by the discovery-strategy
prioritization draft
(`2026-06-25-discovery-strategy-prioritization-design-draft.md`, §6 lever 5):
*"3D-inclined releg genome — re-opens inclined moon-tours (Amalthea is
inclined)."* It is the **out-of-plane sibling of #449** (the low-thrust / DSM
releg genome).

This draft answers ONE question: **what is the smallest, well-bounded capability
that makes a moon-tour releg INCLINATION-AWARE (full 3D position/velocity
matching across an out-of-plane flyby) — so that tours whose moons live on
DIFFERENT orbital planes (Amalthea at the Jovian inner edge; the inclined
ice-giant satellites) become representable at all — and how does it plug into
#449's `Releg` seam, the #291/#287 3D genome, the gauntlet, and the lit-check
without rebuilding the optimiser?**

---

## 0. TL;DR

- **Goal:** make the moon-tour leg geometry **3D**. Today both the discovery-side
  leg solver (`search/discovery_campaign.py::_moon_state`, **line 248**) and the
  validation-side mirror (`data/validation/v2_moontour.py`, via the same
  `_moon_state` import, **line 78**) place every moon on a **circular-COPLANAR**
  orbit with a **hardcoded `z = 0.0`** (`pos = [sma·cosθ, sma·sinθ, 0.0]`,
  `vel = [-v·sinθ, v·cosθ, 0.0]`). A moon-tour cycle therefore *cannot represent*
  an encounter between two moons on different planes — the out-of-plane component
  is structurally absent. The 3D-inclined releg generalises `_moon_state` to a
  full inclined Keplerian element set (`a, e≈0, i, Ω`) and generalises the flyby
  closure from a pure **V∞-magnitude** continuity check to the **3D V∞-VECTOR
  turn-cone** check the geometry actually requires (`max_bend` / `bend_angle` are
  already vector-ready in `core/flyby.py`).
- **Re-opens (exact registry re-open keys):** the moon-tour negatives whose
  `interpretation`/re-open phrase names **"3D/inclined relegs"** — see §1.2.
  Headline key:
  `jupiter-galilean-amalthea-repeated-moon-quasicycler-2026-06-24`
  (verbatim re-open: *"3D/inclined relegs (Amalthea is inclined)"*), plus the
  Jovian/Saturnian `*-vilm` negatives and the `uranus-neptune-regular-moon-…`
  and `repeated-moon-*-sweep` entries that list "3D/inclined relegs" alongside
  low-thrust.
- **Reuse over rebuild — almost everything exists.** The Lambert solver
  (`core/lambert.py`) is **already fully 3D** — it takes `r1, r2` *vectors* and
  solves the transfer plane from them (docstring line 65 even discusses the
  out-of-plane / multi-rev branch). The flyby turn geometry
  (`core/flyby.py::max_bend`, `bend_angle`, `dv_from_turn_deficit`,
  `dv_powered_flyby_periapsis`) is **already V∞-vector-based** and works in 3D
  unchanged. The 3D periodic-orbit machinery (#287/#291,
  `search/cr3bp_general_periodic_3d.py`, `cr3bp_3d_family_tracer.py`) and the
  out-of-plane fingerprint (`genome/known_corpus_3d.py`,
  `spatial_novelty_prefilter.py`) exist. **#449 builds the `Releg` protocol +
  the swap-seam injection.** The new work in #451 is **(a) an inclined moon
  ephemeris** (`a,e,i,Ω` instead of `z=0`), **(b) a 3D-vector flyby-continuity
  residual** (turn-cone feasibility, not just |ΔV∞|), and **(c) one new `Releg`
  backend / flag that switches the leg to the inclined model** — NOT a new
  optimiser, corrector, or Lambert.
- **The exact coplanar seam to generalise (load-bearing):**
  `search/discovery_campaign.py::_moon_state` **line 248** (`z = 0.0` in both
  `pos` and `vel`) — the SINGLE source of the coplanar assumption, consumed by
  BOTH the discovery driver (`_close_one_phasing`, line 487) AND the validation
  mirror (`v2_moontour.py` imports it, line 78, used at line 267). Generalising
  this one function (behind an inclined-aware variant) makes both sites 3D in
  lockstep. The second seam is the **flyby continuity model**
  (`_close_one_phasing` lines 503-528) which today compares only V∞ *magnitudes*
  and must additionally verify the **turn-angle between the in/out V∞ VECTORS is
  inside the achievable bend cone** (the out-of-plane bend is the new binding
  constraint).
- **Composition with #449 (orthogonal axes).** #449 and #451 are **two
  orthogonal axes of the same `Releg` seam**: #449 is *powered-vs-ballistic*
  (does the leg spend a budgeted ΔV?), #451 is *inclined-vs-coplanar* (does the
  leg's geometry carry an out-of-plane component?). A releg may be ANY of the
  four combinations: `{ballistic|powered} × {coplanar|inclined}`. #451's
  `InclinedReleg` composes with #449's `DsmReleg`/`LowThrustReleg` by sharing the
  same `RelegResult` contract and the same two injection sites — an inclined leg
  that also needs to spend ΔV is just `DsmReleg` running on inclined moon states
  with the 3D turn-cone residual. **#451 does NOT duplicate #449's seam; it
  extends the ephemeris + residual that the seam consumes.**
- **Validation / golden target (sourced, non-circular).** The decisive
  inclination-aware analytic reference is **Strange, Russell & Buffington 2007,
  "Mapping the V∞ Globe" (AAS 07-277)** — the V∞ *globe* is precisely the
  out-of-plane V∞-sphere geometry (a flyby rotates the V∞ vector over a sphere;
  pumping changes |V∞|, *cranking* changes its declination/inclination). It is in
  `KNOWN_CORPUS` and cross-referenced in
  `2026-06-17-349-cassini-anchor-topology-label.md` (§"Strange-Russell"). The
  golden asserts the **crank/pump invariants** on the sphere (a ballistic flyby
  preserves |V∞| and rotates the vector within the achievable cone; the
  inclination change per flyby is bounded by `max_bend`). A second, *numeric*
  golden is **Heaton & Longuski 2003 (Galileo-style Uranian satellite tour,
  JSR)** — OCR'd, with published V∞ (≤7.5 km/s arrival, 0.5 km/s Tisserand
  increments) AND an explicit inclination-reduction tour (14 deg reduced via
  Titania flybys; insertion ΔV 0.92 km/s to Ariel) — the inclined-tour ΔV the
  3D releg must bracket. **HONESTY:** neither source is a *state-level*
  per-leg-V∞-vector golden (Cassini Wolf-Smith 1995 has the inclination column
  but **no V∞ column at all**, confirmed §1.3), so the golden is an *invariant*
  golden (crank/pump sphere geometry + bracketed tour ΔV), not a trajectory
  replay — see §6 and the acquisition flag §8.
- **#449 dependency finding: #451 PREFERABLY composes onto #449's merged `Releg`
  seam, but does NOT strictly need it.** The two share the same two injection
  sites; if #449 lands first, #451 adds one backend + the inclined ephemeris +
  the 3D residual against a stable seam (lowest friction). If #451 proceeds
  first or in parallel, it must itself introduce the minimal seam (the same
  dependency-injected `lambert`/leg callable that `_close_one_phasing` ALREADY
  accepts, line 463) — so #451 is *buildable against the design* even with #449
  unmerged. Recommendation: **sequence #451 after #449's seam merges**, to avoid
  two agents racing the same two call sites. See §5.
- **CRITICAL HONESTY — buildability:** **MEDIUM-HIGH** that the *capability*
  (inclined ephemeris + 3D turn-cone residual + composed releg) is buildable from
  currently-digested material — every primitive (3D Lambert, vector flyby cone,
  3D PO machinery, #449 seam) exists. **MEDIUM** on the *golden*: the
  inclination-aware references are digested at the **invariant** level
  (crank/pump sphere, bracketed tour ΔV) but NOT at the **per-leg state** level.
  The FIRST plan task is therefore an **acquisition/transcription** task to firm
  the golden (transcribe the Strange-Russell V∞-globe crank/pump relations to a
  golden YAML, and the Heaton-Longuski Uranian per-flyby V∞/inclination table)
  before any inclined-releg code is asserted against it. See §8.

---

## 1. Which dead regions this re-opens — and why a PLANAR releg cannot represent them

### 1.1 Why a coplanar releg structurally cannot close an inclined tour

The current moon-tour leg geometry places EVERY moon on the **same plane**:
`_moon_state` (discovery_campaign.py:248) returns
`pos = [a·cosθ, a·sinθ, 0.0]`, `vel = [-v·sinθ, v·cosθ, 0.0]` — the third
component is a literal `0.0`. The validation mirror imports the SAME function
(v2_moontour.py:78). Consequences:

1. **The out-of-plane defect is invisible, not small.** When two moons orbit on
   different planes (mutual inclination `Δi`), the *real* leg must absorb an
   out-of-plane velocity component proportional to `V·sin(Δi)`. The coplanar
   model sets this to zero by construction, so it either (a) silently reports a
   feasible tour that does not exist in 3D, or (b) — because the moons are forced
   coplanar — never poses the inclined leg at all. Either way the inclined tour
   is **not representable**: it is not that the residual is large, it is that the
   relevant degree of freedom is absent.
2. **The flyby continuity check is magnitude-only.** `_close_one_phasing`
   (lines 503-528) compares only `|V∞_in|` vs `|V∞_out|` (a scalar). A real
   flyby can rotate the V∞ *vector* only within an achievable **bend cone**
   (`max_bend(μ, rp_min, V∞)`); an out-of-plane tour demands a turn the cone may
   not provide. The coplanar model cannot even pose this question because both
   V∞ vectors lie in `z = 0`. The 3D releg replaces "magnitudes match" with
   "the in→out V∞-VECTOR turn is inside the cone (or the deficit is paid)".

This is the **dual** of #449's argument. #449 said: the ballistic leg has *no
knob to spend ΔV* (a budget problem). #451 says: the coplanar leg has *no
out-of-plane degree of freedom* (a representation problem). They are orthogonal:
#449 widens the *energy* axis, #451 widens the *geometry* axis.

### 1.2 Exact re-open keys (from `data/empty_regions.jsonl`)

| `region_id` (verbatim) | Why coplanar/ballistic fails | Re-open phrase (verbatim) |
|---|---|---|
| `jupiter-galilean-amalthea-repeated-moon-quasicycler-2026-06-24` | Amalthea is **inclined** (~0.38 deg to Jupiter's equator, but the Galilean moons are near-equatorial — the mutual-plane mismatch is the unrepresented DOF); fresh Amalthea cycles infeasible | *"Re-sweep warranted only if a NEW capability subsumes this: **3D/inclined relegs (Amalthea is inclined)**, low-thrust relegs …"* |
| `uranus-neptune-regular-moon-endgame-vilm-2026-06-23` | disjoint Tisserand contours (STRUCTURAL); Uranian system is tilted ~98 deg, moons are in the tilted equatorial plane | strategy-draft D3: *"Re-open only with a low-thrust or **3D/inclined releg** genome (Region F)"* (structural caveat below) |
| `jovian-IEG-vilm-2026-06-09` | high-V∞ basin | strategy-draft §1.3 row 1 re-open: *"**3D/inclined relegs**, low-thrust"* |
| `jovian-perm-vilm-2026-06-09` | high-V∞ basin | same |
| `saturnian-titan-vilm-2026-06-09` | high-V∞ basin | same |
| `saturnian-titan-endgame-vilm-2026-06-10` | high-V∞ basin | same |
| `repeated-moon-jupiter-sweep` | empty as far as the coplanar ballistic genome reached | re-open: *"capability-subsuming genome"* (an inclined genome strictly subsumes the coplanar one — coplanar is the `i=0` limit) |
| `repeated-moon-saturn-sweep` | empty | same |
| `repeated-moon-uranus-sweep` | empty | same |
| `repeated-moon-neptune-sweep` | empty | same |
| `repeated-moon-mars-sweep` | empty | same (Phobos i≈1 deg, Deimos i≈1.8 deg to Mars equator) |

**The headline target is Amalthea.** It is the ONE registry entry whose re-open
phrase names *inclined relegs FIRST* and explicitly ties the re-open to
inclination (*"Amalthea is inclined"*). The Galilean+Amalthea system is the
primary campaign target; the others are re-tested with the subsuming inclined
method either way (capability-subsumption rule: the coplanar genome is the `i=0`
special case of the inclined genome, so the inclined genome strictly subsumes it
and may re-sweep).

**Honest structural caveat (Uranus/Neptune).** Same as #449's caveat: if the
Tisserand contours are disjoint at *every* V∞, an inclined leg does not bridge
them either — inclination is orthogonal to the radial-spacing problem. #451
re-tests these as a *stronger* negative (now coplanar AND inclined both fail),
and re-stamps the registry. Inclination helps where the **planes** differ
(Amalthea, the tilted ice-giant systems' mutual inclinations), NOT where the
**radii/contours** are disjoint.

### 1.3 The golden-source reality check (why this is honest, not optimistic)

The natural state-level golden would be a published inclined moon-tour with
per-flyby V∞ vectors. **It does not exist in the digested corpus:**

- **Cassini (Wolf-Smith 1995)** publishes a per-encounter table with
  **post-flyby inclination** for 54 encounters (max 76 deg) — the richest
  inclination data we have — but **NO V∞ column anywhere** (digest
  `2026-06-17-digest-wolf-smith-1995-cassini.md` line 117: *"V∞ at any encounter
  is not tabulated"*). It is a *cranking* exemplar without the numbers needed to
  validate a leg's ΔV.
- **Davis 2018 (Saturnian ocean worlds)** is Saturn-Titan polar **tulip** orbits
  + Enceladus NRHO halos — **CR3BP periodic-orbit families**, not an inclined
  inter-moon *tour* with per-leg V∞ (digest `2026-06-17-346-davis-2018-deep-read`).
- **Heaton-Longuski 2003 (Uranian)** has published V∞ (≤7.5 km/s arrival, 0.5
  km/s Tisserand increments) and an inclination-reduction narrative (14 deg via
  Titania flybys, 0.92 km/s to Ariel) — usable as a **bracketed tour-ΔV** golden,
  not a state-level replay.
- **Strange-Russell-Buffington 2007 ("Mapping the V∞ Globe")** is the analytic
  *invariant* reference: it defines the V∞-sphere and the crank (inclination) /
  pump (energy) decomposition — exactly the out-of-plane geometry #451 adds.
  Usable as an **invariant** golden (the relations the 3D residual must satisfy).

So the golden is built at two honest levels — **invariant** (Strange-Russell
crank/pump sphere relations) + **bracketed tour ΔV** (Heaton-Longuski Uranian) —
and the FIRST plan task transcribes them. There is no fabricated state-level
golden.

---

## 2. What a "3D-inclined releg" is, precisely

Reuse #449's `Releg` contract **unchanged** (the seam is shared):

```
releg(r_a, v_a_moon, r_b, v_b_moon, tof_s, mu, n_rev) -> RelegResult
    outputs: vinf_out, vinf_in, dv_kms, feasible
```

#451 changes **what produces `r_a, v_a_moon, r_b, v_b_moon`** and **what
`feasible` means**:

- **Inclined ephemeris (the new input).** Replace the coplanar `_moon_state`
  with `_moon_state_3d(theta0, n, t, a, mu, i, raan)` returning the full inclined
  circular state: rotate the in-plane `(x, y, 0)` state by the rotation
  `R_z(Ω)·R_x(i)` so the moon sits on its real plane. For `i = 0, Ω = 0` this is
  bit-for-bit the current `_moon_state` (regression-preserving, the coplanar
  special case). The inclination/RAAN come from the satellites registry (an
  acquisition-flagged field — §8, the registry currently stores only
  `sma_km`/mean-motion, no `i`/`Ω`).
- **3D V∞-vector turn-cone residual (the new feasibility).** The leg's Lambert
  (already 3D) yields V∞ **vectors** at each moon. The flyby continuity check
  becomes: at each interior flyby, (1) `||V∞_in| − |V∞_out||` ≤ magnitude gate
  (as today), AND (2) the **turn angle** `δ_req = angle(V∞_in, V∞_out)` ≤
  `max_bend(μ_planet, rp_min, |V∞|)` (the achievable cone), else the deficit
  `dv_from_turn_deficit(...)` is the powered cost (composing with #449). The
  out-of-plane component lives entirely in (2): a coplanar tour has `δ_req` in
  the `z = 0` plane; an inclined tour's `δ_req` has an out-of-plane component the
  cone must cover.

**The `InclinedReleg` backend** is then just: run the (3D) Lambert on inclined
moon states, compute the V∞ vectors, evaluate the turn-cone residual.
`dv_kms = 0` if the turn is inside the cone (ballistic-inclined); `dv_kms =
dv_from_turn_deficit(...)` (or the Oberth-credited `dv_powered_flyby_periapsis`)
if it exceeds the cone — at which point it is *also* powered, i.e. the
`{powered × inclined}` quadrant that composes #449 × #451.

**Scoring change.** Identical structure to #449: minimise Σ `dv_kms` over the
cycle subject to feasibility, but now feasibility includes the **3D turn-cone**
constraint at every flyby. A coplanar-only tour is the `i = 0` slice; an inclined
tour is a hit if its per-cycle ΔV (the sum of out-of-plane bend deficits, if any)
lands inside the powered dv-band — OR zero, if the inclined flyby geometry closes
ballistically.

---

## 3. Architecture, new files, and the plug-in seam

```
   tour skeleton (sequence, per-leg tof, n_rev, moon ephemerides + i, Ω)
   [from repeated-moon-*-sweep / Amalthea negative — EXISTING enumerations]
                              │
                              ▼
        ┌───────────────────────────────────────────────────────────┐
        │  core/satellites.py  (EXTEND — add i_deg, raan_deg fields) │
        │   inclination + node from JPL SSD orbital tables (§8 acq)  │
        └───────────────────────────┬───────────────────────────────┘
                                    │ inclined element set
                                    ▼
        ┌───────────────────────────────────────────────────────────┐
        │  search/moon_ephemeris_3d.py  (NEW — small)                │
        │   _moon_state_3d(theta0, n, t, a, mu, i, raan)             │
        │     = R_z(Ω)·R_x(i) · [coplanar state]                    │
        │   i=0,Ω=0 ⇒ EXACTLY discovery_campaign._moon_state (regr.) │
        └───────────────────────────┬───────────────────────────────┘
                                    │ inclined r,v per moon
                                    ▼
        ┌───────────────────────────────────────────────────────────┐
        │  search/releg_solver.py  (#449 — EXTEND with InclinedReleg)│
        │   InclinedReleg backend:                                   │
        │     - 3D Lambert (core.lambert, already 3D) on inclined r,v│
        │     - V∞ VECTORS at each moon                              │
        │     - 3D turn-cone residual (core.flyby.max_bend/bend_angle│
        │       /dv_from_turn_deficit — already vector-ready)        │
        │   composes with #449 DsmReleg (powered × inclined)         │
        └───────────────────────────┬───────────────────────────────┘
                                    │ per-leg RelegResult (3D-aware)
                                    ▼
        ┌───────────────────────────────────────────────────────────┐
        │  search/releg_moontour.py  (#449 driver — EXTEND)          │
        │   - inclined-aware flyby-continuity (turn-cone, not |ΔV∞|) │
        │   - sum out-of-plane bend-deficit ΔV per cycle             │
        └───────────────────────────┬───────────────────────────────┘
                                    │
                                    ▼
   EXISTING  verify/dv_band_acceptance.py (powered window)
             search/literature_check.py + genome/known_corpus_3d.py (V0 3D)
             genome/spatial_novelty_prefilter.py (out-of-plane fingerprint)
             data/validation/v2_moontour.py (V2 — inclined-ephemeris-aware)
                                    │
                                    ▼
        SILVER holding (novel inclined tour) | reproduction | empty-restamp
```

**New / changed files:**

| File | Responsibility | Reuses |
|---|---|---|
| `src/cyclerfinder/core/satellites.py` (EXTEND) | Add `inclination_deg` + `raan_deg` (or `node_deg`) fields to `SatelliteData`, sourced from JPL SSD satellite orbital tables. Default `0.0` keeps every existing row coplanar (regression-preserving). | — |
| `src/cyclerfinder/search/moon_ephemeris_3d.py` (NEW, small) | `_moon_state_3d` inclined circular ephemeris; the `i=0,Ω=0` regression equals `discovery_campaign._moon_state`. | `discovery_campaign._moon_state` (regression anchor) |
| `src/cyclerfinder/search/releg_solver.py` (#449 file — EXTEND) | Add `InclinedReleg`; reuse `RelegResult`. Composes with `DsmReleg` for the powered×inclined quadrant. | #449 `Releg` protocol, `core.lambert`, `core.flyby` |
| `src/cyclerfinder/search/releg_moontour.py` (#449 file — EXTEND) | Inclined-aware flyby-continuity (3D turn-cone); sum out-of-plane bend deficit. | #449 driver, `core.flyby` |
| `data/golden/strange_russell_vinf_globe.yaml` (NEW) | Sourced golden: crank/pump invariants (AAS 07-277) + Heaton-Longuski Uranian bracketed tour ΔV/inclination (JSR 2003). EXPECTED side traces only to the papers (`feedback_golden_tests_sourced_only`). | — |

**The plug-in seam (shared with #449, generalised).** The same two call sites —
`discovery_campaign.py:487`(`_close_one_phasing`, which ALREADY accepts an
injected `lambert`/leg callable, line 463) and `v2_moontour.py:267` — consume the
`Releg` backend. #451 routes the inclined ephemeris into them by passing
`_moon_state_3d` (with the registry `i, Ω`) instead of `_moon_state`, and the
`InclinedReleg` backend instead of `BallisticReleg`. **No optimiser, corrector,
or Lambert is written.** The coplanar path is preserved exactly as the
`i=0,Ω=0,BallisticReleg` configuration (regression test).

**Catalogue-scope boundary.** Identical to #449: a 3D-inclined moon-tour that
closes (ballistically or powered) is a `quasi_cycler`/`mga_tour` SILVER candidate
only after V0 lit-check (now the **3D** lit-check `known_corpus_3d` +
`spatial_novelty_prefilter` out-of-plane fingerprint) + the V2/V3 gauntlet + the
mandatory literature-novelty gate. Amalthea's Callisto-Ganymede-Callisto-class
known structures stay V0-known admissions, not novelty claims.

---

## 4. Data flow

1. **Input:** an inclined tour skeleton `(primary, sequence, per-leg tof, n_rev,
   moon ephemerides INCLUDING i, Ω)` — from the existing Amalthea / repeated-moon
   enumerations, now carrying inclination.
2. **VILM-floor prefilter (cheapest first, unchanged from #449).** Skip legs
   whose analytic leveraging floor already exceeds the band. (Inclination does
   not change the *radial* contour test; it adds the *plane* test next.)
3. **Inclined-ephemeris placement.** Place each moon with `_moon_state_3d` on its
   real `(i, Ω)` plane.
4. **Per-leg 3D Lambert + V∞-vector turn-cone solve.** The (already-3D) Lambert
   gives V∞ vectors; the `InclinedReleg` checks the turn-cone at each flyby. If
   inside the cone → ballistic-inclined (`dv=0`); if the out-of-plane bend
   exceeds the cone → pay `dv_from_turn_deficit` (composes #449).
5. **Cycle close + wrap** (same definition as today, now 3D-vector continuity at
   every flyby incl. the anchor wrap). Sum delivered ΔV.
6. **dv-band gate / novelty pipeline:** 3D V0 lit-check
   (`known_corpus_3d` + `spatial_novelty_prefilter`) → signature dedup →
   `literature_check` → V2 moontour (inclined-aware) → SILVER holding.
7. **Output:** classified ledger; an empty band **re-stamps** the corresponding
   `empty_regions.jsonl` entry with the subsuming inclined method + version
   (capability-subsumption record), never silent deletion.

---

## 5. #449 dependency finding (explicit)

**Finding: #451 SHARES #449's two injection sites and PREFERABLY composes onto
#449's merged `Releg` seam, but is buildable against the design if #449 is
unmerged. Recommendation: sequence #451 AFTER #449's seam lands.**

Evidence and reasoning:

- **The seam is the same two call sites** (`discovery_campaign.py:487`,
  `v2_moontour.py:267`). If two agents edit these concurrently they race.
  Sequencing #451 after #449 eliminates the race and lets #451 reuse the
  `RelegResult` contract, the `releg_solver.py` module, and the
  `releg_moontour.py` driver verbatim — #451 then adds exactly one backend
  (`InclinedReleg`) + the inclined ephemeris + the 3D residual.
- **#451 is buildable against the design alone.** `_close_one_phasing` ALREADY
  takes a dependency-injected `lambert` callable (line 463), so #451 could
  introduce the minimal seam itself if needed. The `Releg` protocol is a thin
  contract; #451 can define `InclinedReleg` against the *documented* contract
  before #449's file exists. So #451 is **not hard-blocked** by #449.
- **The two are ORTHOGONAL capability axes**, both attached to the same seam:
  #449 = powered-vs-ballistic (energy), #451 = inclined-vs-coplanar (geometry).
  The `{powered × inclined}` quadrant is their composition: `DsmReleg` running on
  `_moon_state_3d` states with the 3D turn-cone residual. Neither subsumes the
  other; together they span the full `{ballistic|powered} × {coplanar|inclined}`
  releg space.

**Conclusion:** **soft dependency.** Build #451 on top of #449's merged seam
(lowest friction, no call-site race); fall back to building #451's own minimal
seam against the design only if #449 slips. Either way #451 adds the
inclined-ephemeris + 3D-residual + one backend, never a new optimiser.

---

## 6. Validation strategy (sourced/golden proof the inclined releg works)

**Invariant golden (decisive, analytic — Strange-Russell-Buffington AAS 07-277).**
The V∞-globe relations the 3D residual must satisfy:

- **Pump/crank decomposition.** A ballistic flyby preserves `|V∞|` and rotates
  the V∞ vector over the sphere; the achievable rotation per flyby is bounded by
  `max_bend(μ_planet, rp_min, |V∞|)`. The golden asserts: for a fixed `|V∞|`, the
  `InclinedReleg` reports a flyby **feasible** iff the required V∞-vector turn is
  ≤ `max_bend`, and reports the **inclination change** per flyby consistent with
  the globe's crank relation. EXPECTED traces to the paper's globe geometry, not
  to a number our code computed.
- **Coplanar-limit regression.** With `i=0, Ω=0` for all moons, `_moon_state_3d`
  ≡ `_moon_state` and `InclinedReleg`'s magnitude residual ≡ today's
  `_close_one_phasing` residual **bit-for-bit** on a fixed skeleton (the
  capability is a strict superset; the `i=0` slice must reproduce the existing
  coplanar genome exactly).

**Bracketed tour-ΔV golden (numeric — Heaton-Longuski 2003 Uranian).** Run the
inclined releg on a Uranian (Titania/Oberon/Ariel) skeleton; the per-leg V∞ must
fall under the published ≤7.5 km/s arrival ceiling and the inclination-reduction
ΔV must bracket the published values (14 deg reduction; 0.92 km/s Ariel
insertion). This is a *bracket*, not a state replay (the paper does not publish
per-leg V∞ vectors). EXPECTED traces to the OCR'd Heaton-Longuski tables.

**Structural-emptiness golden (the honesty test).** Run the inclined driver on a
Uranian Ariel→Umbriel leg whose Tisserand contours are disjoint; it must STILL
report **unbridgeable** — inclination does not bridge disjoint radial contours.
Reproduces `uranus-neptune-regular-moon-endgame-vilm-2026-06-23` as an
*inclined* re-test that *also* finds it empty. Proves #451 does not fabricate a
bridge from inclination where the radial physics forbids one.

**Positive-control golden.** A Jovian skeleton that includes Amalthea (the
headline target): the inclined releg must POSE the out-of-plane leg the coplanar
genome could not, and either close it within the cone (ballistic-inclined) or
report the out-of-plane bend deficit as the ΔV (powered-inclined). The
*capability proof* is that the inclined leg is **representable and scored at
all** — which the coplanar genome structurally could not do.

**Why golden-clean.** Every EXPECTED value traces to Strange-Russell's globe
relations / Heaton-Longuski's printed tables — never a number the releg solver
computed (`feedback_golden_tests_sourced_only`). The genome "works" iff it (a)
satisfies the crank/pump sphere invariants, (b) reproduces the coplanar genome
exactly in the `i=0` limit, (c) brackets the Heaton-Longuski Uranian tour ΔV, and
(d) honestly reports disjoint-contour legs unbridgeable even when inclined.

**Honesty on the golden's reach.** Passing it proves the *3D-inclined releg
capability* and re-opens the inclined moon-tour regions (Amalthea first) by making
their out-of-plane legs representable. It does NOT, by itself, yield a *novel*
catalogue cycler — whether a feasible inclined Amalthea cycle lands inside the
band AND clears the 3D lit-novelty gate is the open empirical question a separate
campaign answers.

---

## 7. Risks + kill-criteria

| Risk | Mitigation / kill-criterion |
|---|---|
| **Inclined ephemeris drifts the coplanar regression** (the `i=0` slice no longer reproduces `_moon_state`). | Regression test asserts `_moon_state_3d(i=0,Ω=0)` ≡ `_moon_state` bit-for-bit on a fixed `(θ0, n, t, a, μ)` grid. Any drift fails the build. |
| **Inclination adds a DOF but no real moon has enough mutual inclination to matter** (Galileans are near-equatorial; Amalthea ~0.4 deg). | This is a real possibility and a legitimate NEGATIVE. KILL: if across the Jovian+Amalthea and Mars-moon skeletons the out-of-plane bend deficit is always below the magnitude gate (inclination contributes < 0.05 km/s), the inclined DOF is *physically inert* for these systems — re-stamp the registry with the inclined method finding "inclination representable but immaterial here", and reserve #451 for the genuinely tilted systems (ice-giant mutual inclinations). A clean "representable but immaterial" is a success under the discipline. |
| **Uranus/Neptune stays empty even inclined** (disjoint contours, expected). | Legitimate STRONGER negative — re-stamp `uranus-neptune-regular-moon-endgame-vilm` with the inclined method+version. Inclination is orthogonal to radial-contour disjointness; the VILM-floor prefilter makes this cheap. |
| **The 3D turn-cone residual "closes" a tour that is not a cycler** (no repeat). | V2-moontour gate (≥3 cycles, drift-bounded) + closed-cycle wrap continuity stay mandatory; one inclined leg-close is necessary-not-sufficient. |
| **No state-level golden → an invariant golden could pass a subtly-wrong residual.** | Triangulate: (1) coplanar-limit bit-for-bit regression, (2) crank/pump sphere invariant, (3) Heaton-Longuski bracket, (4) disjoint-contour negative. A wrong residual must violate at least one. PLUS the FIRST plan task transcribes the goldens before code (§8). |
| **The two swap sites drift out of sync** (discovery vs validation ephemeris differ). | The `_moon_state_3d` function + the `Releg` protocol are the single sources of truth; both sites import them. Regression asserts both produce identical inclined states on a fixed skeleton. |
| **Scope creep into a novelty claim from a known inclined tour** (Cassini cranking is published). | Hard boundary: the 3D lit-novelty gate (`known_corpus_3d` + `spatial_novelty_prefilter` + §16.5) stays mandatory. Cassini-class cranking is a V0-known admission, not novelty. |

---

## 8. Prerequisite corpus / environment acquisitions (honest flags)

1. **REQUIRED — golden transcription (FIRST plan task).** Unlike #449 (whose
   VILM golden already passes in-repo), #451's inclination-aware golden is NOT
   yet a YAML. Two cheap transcriptions, both from already-OCR'd/digested sources:
   - **Strange-Russell-Buffington 2007 "Mapping the V∞ Globe" (AAS 07-277)** —
     transcribe the crank/pump sphere relations (the invariant the 3D residual
     must satisfy) into `data/golden/strange_russell_vinf_globe.yaml`. The paper
     is in `KNOWN_CORPUS`; the relations are summarised in
     `2026-06-17-349-cassini-anchor-topology-label.md`. **Confirm the per-relation
     numbers from the paper body before asserting** (the cassini note is a
     pointer, not a full transcription — flag if the AAS 07-277 PDF body needs a
     deeper read for the exact crank-angle relation).
   - **Heaton-Longuski 2003 (Uranian, JSR)** — transcribe the published per-flyby
     V∞ / inclination table (the .txt is OCR'd at
     `cyclers_pdf/papers/heaton-longuski-…2.3981.txt`; V∞ ≤7.5 km/s, 0.5 km/s
     Tisserand increments, 14 deg inclination reduction, 0.92 km/s Ariel
     insertion are visible) into the same golden YAML as the bracketed tour-ΔV
     target.
2. **REQUIRED — registry inclination fields.** `core/satellites.py` stores
   `sma_km` + derived mean motion but **no `i`/`Ω`**. Add `inclination_deg` +
   `node_deg` from the JPL SSD satellite orbital-element tables (the same SSD
   source already cited for `sma_km`, line 93-99 — non-circular, sourced). Default
   `0.0` preserves every existing coplanar row. Small, mechanical, sourced.
3. **No NEW PDF acquisition is strictly required** — both golden sources are
   already in-corpus (Strange-Russell in KNOWN_CORPUS; Heaton-Longuski OCR'd).
   The work is transcription + sourcing, not acquisition. BUT: if the deeper read
   of AAS 07-277 reveals the crank-angle relation is not numerically reproducible
   from the digested material, flag a **targeted deep-read of the AAS 07-277 PDF
   body** as the one possible new corpus task (it is in `cyclers_pdf`, just not
   fully digested at the relation level).
4. **No environment/install gate** — pure-Python on the existing scipy stack
   (3D Lambert, flyby cone, rotations are all numpy).

**Net:** the inclined-releg *capability* is buildable from currently-digested
material (every primitive exists); the *golden* needs a one-task transcription of
two in-corpus sources to the invariant + bracketed-ΔV level. The FIRST plan task
is that transcription — no inclined code is asserted before the golden exists.

---

## 9. Confidence

- **Capability buildable from digested material:** **MEDIUM-HIGH** — 3D Lambert,
  vector flyby cone, 3D PO machinery, and the #449 seam all exist; the new work
  is an inclined ephemeris (one rotation), a 3D turn-cone residual (vector angle
  vs `max_bend`), and one composed backend. No new optimiser.
- **Re-opens real dead regions with a sourced validation target:** **MEDIUM-HIGH**
  — the registry names "3D/inclined relegs" as the explicit re-open key (Amalthea
  first); the golden sources are in-corpus but need transcription to the
  invariant level (the honesty discount).
- **Composition with #449:** **HIGH** — orthogonal axes on a shared seam; the
  `{powered × inclined}` quadrant is a clean composition.
- **#449 dependency:** **soft** — preferably sequence #451 after #449's seam
  merges (avoids the two-call-site race); buildable against the design otherwise.
- **Yields a NOVEL catalogue cycler:** **LOW-MEDIUM and honestly so** — the
  capability makes inclined Amalthea legs representable for the first time, but
  the near-equatorial Galilean planes mean the inclined DOF may be physically
  immaterial there (an explicit kill-criterion, §7). The genuine inclination
  payoff is in the strongly-tilted ice-giant mutual-inclination systems; whether
  any inclined tour there lands inside the band AND clears the 3D lit-novelty
  gate is the open bet. The capability win is likely; the novel-cycler payoff is
  speculative.

---

## 10. Open questions for the user

1. **Build order vs #449:** sequence #451 after #449's `Releg` seam merges
   (recommended, avoids the call-site race), or build #451's own minimal seam in
   parallel against the design? Recommendation: **after #449**.
2. **Golden depth:** is the invariant + bracketed-ΔV golden (Strange-Russell
   crank/pump sphere + Heaton-Longuski Uranian bracket) sufficient for the
   capability gate, or is a targeted deep-read of AAS 07-277 warranted to firm
   the crank-angle relation to a state-level number first?
3. **Headline target scope:** Amalthea-Galilean first (the named re-open key), or
   lead with the strongly-tilted ice-giant mutual-inclination systems where the
   inclined DOF is more likely *material* (per the §7 "immaterial here"
   kill-criterion)?
4. **Catalogue policy for an inclined tour:** confirm a novel inclined moon-tour
   that closes is a `quasi_cycler`/`mga_tour` SILVER candidate following the same
   V0-known vs SILVER-novel split as #449 (Callisto-Ganymede-Callisto exemplar).
5. **Registry inclination source:** confirm JPL SSD satellite orbital-element
   tables (the existing `sma_km` source) for `i`/`Ω`, defaulting `0.0` to
   preserve all coplanar rows.
