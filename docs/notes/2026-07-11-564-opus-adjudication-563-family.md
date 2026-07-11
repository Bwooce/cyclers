# #564 — Opus adjudication of the #563 30-closure Uranian symmetric-cycler family

**Date**: 2026-07-11
**Task**: #564 (P0, judgment-only — Stage-3/writeback gate for the #558→#563 thread)
**Model**: Opus (trust-bearing discovery-verdict judgment, per #561's precedent for this
territory). Fable second-opinion pass to follow per #561's two-model discipline.
**Scope**: adjudicate the full 30-closure #563 symmetric-closure survivor table
(`data/enumerate_563_symmetric_closures.jsonl`, read directly, all 60 raw passes → 30
deduped physical loops). Produce a prioritized gauntlet plan, a #312 reframing
recommendation, and a go/no-go on the 2 residual asymmetric near-closures.
**Discipline**: research/adjudication only — NO V1-V4-strict gauntlet run, NO
`data/catalogue.yaml` writeback, NO asymmetric-search build. Recommendations only.

---

## 1. What the 30 closures actually are

#563 directly constructed and gate-checked every symmetric-closure candidate
(`rel_offset ∈ {0°, 180°}`, matched leg revolution counts, exactly-commensurate
`tof = n·T_syn/2`) across the 6 non-Miranda Uranian moon pairs, both anchor-flyby
directions. 60 raw gate-passing rows dedupe (direction-mirror pairs) to **30 unique
physical symmetric closures**. All 30 close to machine precision (worst residual
1.6e-13 km/s), with DOP853 two-body cross-checks passing identically. This is the
"perpendicular-crossing" symmetric periodic-orbit signature — **exact by construction,
not numerical luck** — so residual/exactness is now a dead axis for triage (a change
from #561/#562, where residual was the primary discriminant).

The 30, organized by pair (V∞ = the two distinct per-encounter magnitudes km/s;
`max bend` = the largest single-encounter bend deg = the binding flyby-feasibility
constraint; `†` = contains #312 or its known sibling; `NEW` = one of #563's 11
coverage-gap recoveries; `LIT` = Titania-Oberon/Oberon-Titania direction, mandatory
Canales/Kumar citation if it proceeds):

| # | pair | n | n_rev | tof (d) | V∞ (km/s) | max bend° | flag |
|---|---|---|---|---|---|---|---|
| 1 | Ariel-Umbriel | 1 | (0,0) | 3.216 | 0.979 / 1.300 | 13.98 | shortest tof of all 30 |
| 2 | Ariel-Umbriel | 2 | (1,1) | 6.432 | 1.123 / 1.305 | 10.93 | |
| 3 | Ariel-Umbriel | 3 | (2,2) | 9.648 | 0.884 / 1.492 | 16.71 | |
| 4 | Ariel-Umbriel | 3 | (2,2) | 9.648 | 1.547 / 1.310 | 8.31 | NEW |
| 5 | Ariel-Titania | 3 | (0,0) | 5.321 | 1.231 / 1.719 | 9.62 | |
| 6 | Ariel-Titania | 5 | (1,1) | 8.868 | 1.106 / 1.031 | 23.41 | |
| 7 | Ariel-Oberon | 5 | (0,0) | 7.752 | 1.521 / 1.829 | 8.07 | NEW |
| 8 | Ariel-Oberon | 9 | (1,1) | 13.953 | 1.430 / 1.407 | 13.02 | |
| 9 | Ariel-Oberon | 10 | (1,1) | 15.504 | 1.670 / 1.615 | 10.15 | NEW |
| 10 | Umbriel-Titania | 1 | (0,0) | 3.954 | 1.230 / 1.006 | 24.36 | |
| 11 | Umbriel-Titania | 3 | (1,1) | 11.863 | 1.173 / 1.766 | 10.19 | NEW |
| 12 | Umbriel-Titania | 4 | (0,0) | 15.818 | 1.617 / 2.342 | 5.59 | NEW — lowest bend of all 30 |
| 13 | Umbriel-Titania | 4 | (1,1) | 15.818 | 1.633 / 1.706 | 9.76 | NEW |
| 14 | Umbriel-Titania | 4 | (1,1) | 15.818 | 1.314 / 2.139 | 8.27 | NEW |
| 15 | Umbriel-Titania | 4 | (2,2) | 15.818 | 0.771 / 0.978 | 25.48 | SUPERSEDE |
| 16 | Umbriel-Oberon | 2 | (0,0) | 5.987 | 1.691 / 1.281 | 15.35 | |
| 17 | Umbriel-Oberon | 3 | (0,0) | 8.980 | 1.080 / 1.520 | 11.83 | |
| 18 | Umbriel-Oberon | 5 | (0,0) | 14.967 | 1.494 / 2.124 | 6.50 | NEW |
| 19 | Umbriel-Oberon | 5 | (1,1) | 14.967 | 0.926 / 0.965 | 24.69 | **† = #312** |
| 20 | Umbriel-Oberon | 6 | (1,1) | 17.960 | 1.224 / 1.995 | 9.42 | |
| 21 | Umbriel-Oberon | 7 | (1,1) | 20.954 | 1.393 / 1.660 | 9.65 | NEW |
| 22 | Umbriel-Oberon | 7 | (2,2) | 20.954 | 0.930 / 0.672 | 42.11 | † #312 n=7 sibling |
| 23 | Titania-Oberon | 1 | (0,0) | 12.316 | 1.799 / 1.608 | 10.23 | NEW, LIT |
| 24 | Titania-Oberon | 1 | (0,0) r180 | 12.316 | 2.162 / 1.968 | 7.03 | LIT |
| 25 | Titania-Oberon | 1 | (1,1) | 12.316 | 2.169 / 1.255 | 15.91 | LIT |
| 26 | Titania-Oberon | 2 | (0,0) | 24.632 | 0.918 / 1.524 | 28.13 | LIT |
| 27 | Titania-Oberon | 2 | (1,1) | 24.632 | 0.827 / 0.959 | 32.92 | NEW, LIT |
| 28 | Titania-Oberon | 2 | (1,1) r180 | 24.632 | 1.092 / 1.507 | 21.30 | LIT |
| 29 | Titania-Oberon | 2 | (2,2) | 24.632 | 1.217 / 0.616 | 47.21 | LIT |
| 30 | Titania-Oberon | 2 | (2,2) r180 | 24.632 | 0.657 / 0.335 | 87.61 | LIT |

Family spans all 6 non-Miranda pairs and all four moons (Ariel, Umbriel, Titania,
Oberon) as both anchor and flyby body. #312 (row 19) is one member; its known n=7
sibling (row 22) is another.

---

## 2. Triage axes now that exactness is uniform

With residual dead, the discriminating axes are:

- **max bend angle** (binding flyby-feasibility + real-ephemeris robustness margin).
  Low-V∞ members buy their commensurability with deep flybys → high bend. Rows 29-30
  (bend 47°, 88° at V∞ 0.62/0.34) demand near-surface/impossible periapses and are the
  **most fragile under real perturbations** and least mission-interesting. Rows with
  max bend ≲ 10° (4, 5, 7, 12, 18, 20, 21, 24) have the cleanest geometric margin.
- **V∞ magnitude**: the family clusters ~0.9-2.3 km/s. Moderate V∞ (~1.0-1.7) is the
  sweet spot — low enough to be ballistic-comfortable, high enough to keep bend modest.
- **tof**: short-tof members accumulate less real-ephemeris drift over the V4-strict
  window and give the fairest V4 shot. This matters concretely: #312's own V4-strict
  failures already cluster at the URA111 kernel edge (1900-2099), and drift scales with
  integration span. Prefer the shortest-tof clean member per pair.
- **literature risk**: the 8 Titania-Oberon closures (rows 23-30) sit in the largest,
  most literature-exposed basin (§4).

**Exactness in the idealized (patched-conic + #324-bend + DOP853 two-body) model is NOT
the same as V4-strict bounded-drift under DE440/URA111.** That gap is exactly why we
cannot declare all 30 "validated" from the enumeration alone, and why a representative
real-ephemeris gauntlet subset still has real evidentiary value.

---

## 3. Gauntlet recommendation — 5 candidates, one best-representative per remaining pair

**The claim the gauntlet must support** is "#312 is the first-documented member of a
symmetric-closure family spanning the non-Miranda Uranian moon pairs." That claim is
substantiated by **at least one V4-strict-validated real-ephemeris member per pair**,
not by validating every algebraic sibling within a pair. #312 (Umbriel-Oberon) is
already fully validated (V2 10-cycle, REBOUND, V4-scipy/strict, 85/85 interior epochs).
So the marginal evidence needed is **one validated representative in each of the other
5 pairs**. That is the natural, non-arbitrary cutoff: it is the minimum set that
upgrades the catalogue claim from "one validated member + 29 idealized-exact closures"
to "family demonstrated across all 6 pairs."

Selection rule per pair: **shortest-tof, lowest-order, flyby-feasible member** (best V4
odds + simplest to defend as "the representative"), with the Titania-Oberon pick chosen
for best geometric quality since it also has to carry the literature test.

| priority | # | pair | n | n_rev | tof (d) | V∞ | max bend° | why this one |
|---|---|---|---|---|---|---|---|
| P1 | 24 | Titania-Oberon | 1 | (0,0) r180 | 12.316 | 2.162/1.968 | 7.03 | LIT-risk pair — must be tested to discharge Canales/Kumar; cleanest bend of the 8 T-O closures |
| P1 | 1 | Ariel-Umbriel | 1 | (0,0) | 3.216 | 0.979/1.300 | 13.98 | shortest tof of all 30 → best V4 shot; lowest order |
| P1 | 5 | Ariel-Titania | 3 | (0,0) | 5.321 | 1.231/1.719 | 9.62 | only geometrically-clean A-T member (n=5 is 23.4° bend); short tof |
| P1 | 7 | Ariel-Oberon | 5 | (0,0) | 7.752 | 1.521/1.829 | 8.07 | lowest order + lowest bend for A-O; NEW recovery |
| P1 | 10 | Umbriel-Titania | 1 | (0,0) | 3.954 | 1.230/1.006 | 24.36 | shortest tof for U-T (bend 24.4° ≈ #312's own 24.7°, known-feasible) |

**Recommended primary gauntlet set: 5 candidates (rows 24, 1, 5, 7, 10).** Together with
the already-validated #312 this gives one real-ephemeris-validated member in each of the
6 non-Miranda pairs — the complete, defensible basis for the "family across the pairs"
catalogue claim.

**Tier-2 (optional, only if the P1 set all pass and extra budget exists)** — to showcase
that individual pairs carry multiple validated members, not just intra-pair diversity on
paper:
- Row 12 (Umbriel-Titania n=4 (0,0), max bend **5.59° — the cleanest closure in the
  whole family**) as a "best-geometry showcase," and/or
- Row 22 (the #312 n=7 sibling) which is already essentially known and cheap to confirm.

**Explicitly NOT recommended for the gauntlet at any tier:** the extreme-bend, low-V∞
Titania-Oberon rows 29 (bend 47°) and 30 (bend 88°, V∞ 0.335) and Umbriel-Oberon row 22's
partner-class row-47 analog (bend 42°). They are algebraically exact but demand
near-surface flybys, are the least robust to real perturbations, and add nothing to the
family-existence claim that a clean sibling in the same pair does not already carry.

**Cutoff justification (not an arbitrary top-N):** the 30 members are exact-by-
construction from a single symmetric-closure principle; within a given pair they differ
only in `(n, n_rev, rel_offset)` and are all guaranteed-exact in the idealized model.
Real-ephemeris validation of a second, third, … member of an *already-represented* pair
is evidentially redundant for the census claim (it re-tests the same idealized→real gap
in the same dynamical neighborhood). The information-bearing quantity is **one validated
member per previously-unrepresented pair** — hence 5, not 30, and not an arbitrary "top-8
by V∞."

---

## 4. Literature risk — Titania-Oberon direction (grounded, not assumed)

Verified against the corpus (`docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md`,
`docs/notes/2026-06-20-digest-canales-howell-2023.md`,
`docs/notes/2026-06-20-digest-kumar-2025.md`, `CORPUS_INDEX.md`):

- **Canales-Howell-Fantino(-Gilliam) MMAT / FTLE work** (arXiv:2110.03683 MMAT;
  arXiv:2308.10029 FTLE-maps): moon-to-moon transfers via **halo/libration-orbit
  invariant-manifold patching** (2BP-CR3BP MMAT), applied to a Titania-Oberon halo-to-halo
  transfer. Structurally adjacent (Uranian moon-to-moon transfer) but a **different
  topology** — a one-shot manifold transfer between libration-point orbits, not a periodic
  ballistic repeated-symmetric-flyby cycler.
- **Kumar 2025** (arXiv:2509.12675 and the Uranus-Oberon multi-shooting MMR work surfaced
  in #328): **single-moon** resonant orbits (Uranus-Oberon 3:4/4:5/5:6 exterior, 4:3/5:4/6:5
  interior MMRs) with heteroclinic transitions. Again structurally adjacent (Uranian
  resonant transport) but **single-moon-resonance based**, not a two-moon symmetric-flyby
  cycler.

**Verdict:** #562's literature-triage stands mechanically — **neither is a direct topology
match**, but any Titania-Oberon candidate that proceeds to the gauntlet (i.e. P1 row 24)
**must cite both** Canales-Howell-Fantino 2021 and Kumar 2025 as adjacent prior art, per
#562's own verdict. The other 5 pairs carry only the general Uranian-system baseline
(#328) and no direction-specific obligation.

---

## 5. #312 reframing — YES, reframe, but keep two counts distinct

**Recommendation: reframe #312's catalogue prose** to place it as the first-documented
representative of the now-enumerated family — **but the prose must not conflate "enumerated
symmetric closures" with "validated quasi-cyclers."** These are two different, both-true
counts:

- **30** = symmetric-closure family members, **exact in the idealized model**
  (patched-conic + #324 bend + DOP853 two-body), exhaustively enumerated for the symmetric
  class within #558's tested tof range (`tof_scale ≤ 3.0`). This is a *mathematical* census.
- **1 (→ up to 6)** = members **gauntlet-validated on real ephemeris**: #312 today; the P1
  subset (§3) would bring it to 6 if they pass.

Writing "a family of 30 validated Uranian quasi-cyclers" would **overclaim** — 29 of the
30 have not been shown to hold bounded drift under DE440/URA111, and the idealized→real
gap is real (#312's own V4 failures at kernel edges prove it is not automatic). The
honest reframing, for the eventual (separate) writeback step, is roughly:

> "#312 is the first-documented representative of a family of **30 symmetric periodic
> closures** (rel_offset ∈ {0°, 180°}, exactly-commensurate `tof = n·T_syn/2`) spanning
> the six non-Miranda Uranian moon pairs (Ariel/Umbriel/Titania/Oberon), enumerated
> exhaustively for the **symmetric-closure class within the searched TOF range** (#563).
> The enumeration is exact in the idealized (patched-conic) model; #312 remains the only
> member validated on real ephemeris through the full V1-V4-strict gauntlet.
> **Qualifier:** this is the symmetric-class census only — it is NOT proven exhaustive over
> asymmetric closures (two genuinely asymmetric near-closures are known; see #564)."

This reframing is **safe to make once at least the P1 T-O and one non-U-O member have been
gauntlet-validated** (so the "family" language rests on ≥2-3 validated pairs, not one). If
the writeback is wanted sooner, the prose above is already defensible as written because it
explicitly separates the enumerated-30 from the validated-1. **The actual edit is out of
scope for this pass** (per #564) — this is the recommended text + the guardrail, to be
applied in the deliberate writeback step.

---

## 6. Asymmetric near-closures — NO-GO / defer; symmetric census closes the thread

The two residual asymmetric near-closures #563 correctly declined to enumerate:
- Oberon-Titania n=3 n_rev=(0,0), tof 36.95 d, rel_offset 114.15°, residual **7.9e-3 km/s**
  (also sits just outside the strict `tof_scale ≤ 3.0` bound);
- Oberon-Titania n=2 n_rev=(0,1), tof 24.63 d, rel_offset 268.19°, residual **3.5e-2 km/s**.

**Recommendation: NO-GO on a dedicated adaptive/basin-width-aware asymmetric search;
the symmetric census is sufficient to close out the #558→#563 thread.** Reasoning:

1. **The original question is already answered.** "Is #312 isolated or a family?" —
   definitively a family (30 symmetric members). The asymmetric question is a *new,
   lower-value* research direction, not a loose end in the original thread.
2. **Both candidates are in the single most literature-exposed basin** (Titania-Oberon),
   where the novelty margin against Canales/Kumar is thinnest. Finding an asymmetric T-O
   closure risks re-deriving a structure already implicit in the MMAT/MMR literature — poor
   novelty ROI.
3. **The method would be fundamentally more expensive with no exhaustiveness guarantee.**
   The symmetric class is provably finite/enumerable (that is what made #563 clean and
   cheap). The asymmetric class is a continuous 2D `(rel_offset × tof)` basin-width-aware
   adaptive search with no finite enumeration and no proof of completeness — it can only
   ever produce "found some / found none so far," never a census. That is a categorically
   weaker deliverable than the symmetric result already in hand.
4. **residual 3.5e-2 km/s is near the 0.05 gate floor** — the second candidate may be a
   genuine near-miss rather than a true closure; chasing it is speculative even before the
   novelty question.

**If it is ever greenlit** (a standalone future task, not this thread), scope it as:
Titania-Oberon-only, adaptive 2D `(rel_offset, tof)` refinement **seeded at the two known
near-closures**, **gated behind first resolving the Canales/Kumar novelty question** (no
point localizing an asymmetric T-O closure that turns out to be a known manifold/MMR
structure), and registered as method-versioned in the negative-results registry if it
comes up empty. Rough size: a genuine multi-day adaptive search + a literature-adjudication
gate — materially larger than #563, and **P2/P3 at best** given items 1-4. Recommendation
stands at **defer**.

---

## 7. Summary of recommendations

1. **Gauntlet subset = 5** (rows 24, 1, 5, 7, 10 in §3): one shortest-tof/cleanest
   representative per each non-Umbriel-Oberon pair; with #312 this validates one member in
   every non-Miranda pair. T-O row 24 is mandatory (discharges the Canales/Kumar test).
   Tier-2 optional (rows 12, 22). Do NOT gauntlet all 30 — 24 of them are exact-by-
   construction siblings within already-represented pairs (redundant for the census claim).
2. **Reframe #312** as first-documented member of a **30-member symmetric-closure family**
   (§5 draft prose), keeping the enumerated-30 (idealized-exact) count strictly distinct
   from the validated-1(→6) (real-ephemeris) count, with the explicit symmetric-class +
   searched-tof-range qualifier. Apply in the separate writeback step, ideally after ≥2-3
   pairs are gauntlet-validated. Edit NOT made here (out of scope).
3. **Asymmetric near-closures: defer (NO-GO).** The symmetric census closes the thread;
   the asymmetric case is a lower-value, more-expensive, non-exhaustive, high-lit-risk
   standalone future task (P2/P3), gated behind the Canales/Kumar novelty question if ever
   pursued.

Two-model discipline: this Opus verdict should get the #561-precedent **Fable
second-opinion** before any Stage-3 dispatch or writeback.
