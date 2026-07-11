# #565 — Fable-corrected re-adjudication of the #563 30-closure Uranian symmetric-cycler family

**Date**: 2026-07-11
**Task**: #565 (P0, judgment-only — corrects #564 §2/§3 before any Stage-3 dispatch / writeback)
**Model**: Opus (correcting its own #564 output against Fable's adversarial second opinion, per
the #561 two-model discipline).
**Supersedes**: the geometric-quality ranking in **#564 §2** and the gauntlet-representative
picks / exclusion reasoning in **#564 §3**, plus two citation/contradiction nits. The original
`docs/notes/2026-07-11-564-opus-adjudication-563-family.md` stays as the historical record;
**this note is the corrected reference going forward.**
**Leaves standing from #564** (Fable re-verified these; NOT reopened here): §1 data verification
(all 30 rows match `data/enumerate_563_symmetric_closures.jsonl`, 14-CONFIRM/5-SUPERSEDE/11-NEW
against #562, #312 identity/status); the "**5 not 30**" cutoff count and its logic (one validated
representative per previously-unrepresented pair); §5 #312-reframing prose; §6 asymmetric-closures
deferral (NO-GO); both asymmetric near-closure numbers.

---

## 0. The error being corrected (one sentence)

#564 read `max_bend_deg_per_encounter` as a **demanded** turn angle (so it ranked *high* bend as
"fragile/near-impossible" and *low* bend as "clean/comfortable"). The field is the opposite: it is
the **maximum achievable ballistic deflection** at the safe-altitude periapsis — a pure function of
(body, V∞) — gated as a **floor** (`flyby_is_useful()` passes iff `max_bend_deg >= 5.0°`, the #324
gate; `src/cyclerfinder/search/physical_sanity.py:80,100-103`). **Higher = more capability margin
above the reject floor; lower = closer to the floor (more marginal).** Nothing in the pipeline ever
computes a *required* turn angle — `residual_at_point` in the #558/#563 scripts matches V∞ magnitudes
only. So every #564 ranking built on the bend axis is inverted, and the correct margin statistic is
the **minimum per-encounter achievable bend** across a loop's encounters (the encounter closest to
the floor is the binding one — this is the same "min per-encounter" statistic #563's own NEW-closures
characterization used).

---

## 1. Corrected geometric-quality axis (supersedes #564 §2)

Recomputed directly from `data/enumerate_563_symmetric_closures.jsonl` (60 raw passes → 30 deduped
physical loops; `min = min(max_bend_deg_per_encounter)`). The full 30 rows keep #564 §1's row
numbering; only the quality reading changes. **`min bend°` = minimum per-encounter achievable
deflection = capability margin above the 5.0° floor. Bigger is more robust, not less.**

Ranked most-MARGINAL (near the floor) → most-ROBUST:

| min bend° | row | pair | n | n_rev | tof (d) | reading |
|---|---|---|---|---|---|---|
| 5.13 | 16 | Umbriel-Oberon | 2 | (0,0) | 5.987 | most marginal of all 30 |
| 5.21 | 9 | Ariel-Oberon | 10 | (1,1) | 15.504 | |
| 5.39 | 12 | Umbriel-Titania | 4 | (0,0) | 15.818 | **3rd-most marginal — NOT the "cleanest" (see §3)** |
| 5.49 | 13 | Umbriel-Titania | 4 | (1,1) | 15.818 | |
| 6.03 | 4 | Ariel-Umbriel | 3 | (2,2) | 9.648 | |
| 6.09 | 18 | Umbriel-Oberon | 5 | (0,0) | 14.967 | |
| 6.22 | 7 | Ariel-Oberon | 5 | (0,0) | 7.752 | A-O rep (tof/order, §2) |
| 6.23 | 25 | Titania-Oberon | 1 | (1,1) | 12.316 | lowest-margin T-O |
| 6.27 | 24 | Titania-Oberon | 1 | (0,0) r180 | 12.316 | **#564's T-O pick — 2nd-worst T-O margin (see §2)** |
| 6.40 | 14 | Umbriel-Titania | 4 | (1,1) | 15.818 | |
| 6.51 | 3 | Ariel-Umbriel | 3 | (2,2) | 9.648 | |
| 6.85 | 20 | Umbriel-Oberon | 6 | (1,1) | 17.960 | |
| 6.99 | 8 | Ariel-Oberon | 9 | (1,1) | 13.953 | **larger margin than A-O rep row 7** |
| 7.41 | 21 | Umbriel-Oberon | 7 | (1,1) | 20.954 | |
| 8.37 | 2 | Ariel-Umbriel | 2 | (1,1) | 6.432 | |
| 8.43 | 1 | Ariel-Umbriel | 1 | (0,0) | 3.216 | A-U rep (shortest tof of all 30, §2) |
| 8.85 | 23 | Titania-Oberon | 1 | (0,0) | 12.316 | **corrected T-O rep — best margin of the short-tof T-O trio (§2)** |
| 9.16 | 11 | Umbriel-Titania | 3 | (1,1) | 11.863 | |
| 9.25 | 5 | Ariel-Titania | 3 | (0,0) | 5.321 | A-T rep (tof/order, §2) |
| 9.35 | 10 | Umbriel-Titania | 1 | (0,0) | 3.954 | U-T rep (shortest tof, §2) |
| 11.24 | 6 | Ariel-Titania | 5 | (1,1) | 8.868 | **larger margin than A-T rep row 5** |
| 11.28 | 26 | Titania-Oberon | 2 | (0,0) | 24.632 | |
| 11.33 | 17 | Umbriel-Oberon | 3 | (0,0) | 8.980 | |
| 11.50 | 28 | Titania-Oberon | 2 | (1,1) r180 | 24.632 | |
| 15.43 | 22 | Umbriel-Oberon | 7 | (2,2) | 20.954 | **#312 n=7 sibling — robust, ≈#312's own margin (see §4)** |
| 15.55 | 19 | Umbriel-Oberon | 5 | (1,1) | 14.967 | **#312 itself — sits at the ROBUST end** |
| 17.77 | 29 | Titania-Oberon | 2 | (2,2) | 24.632 | **#564 excluded as "fragile" — actually more robust than #312** |
| 21.20 | 15 | Umbriel-Titania | 4 | (2,2) | 15.818 | |
| 24.94 | 27 | Titania-Oberon | 2 | (1,1) | 24.632 | **highest-margin practical member; ≈#312 V∞ profile (§2 alt)** |
| 45.29 | 30 | Titania-Oberon | 2 | (2,2) r180 | 24.632 | **#564 excluded as "impossible" — most robust of all 30** |

**The single most important consequence**: **#312 (row 19) itself sits near the robust end
(min 15.55°)**, and the two rows #564 excluded as "most fragile / near-surface / impossible"
(rows 29 at 17.77° and 30 at 45.29°) have **larger** capability margins than #312 does. The
inverted axis exactly reversed the family's robustness ordering.

Which axes are *actually* discriminating now that exactness is uniform (corrected §2):
- **min per-encounter achievable bend** — capability margin above the 5° reject floor. Marginal
  members (rows 16, 9, 12, 13 at ~5.1–5.5°) sit just above the floor and have the *least* real-
  ephemeris robustness cushion; robust members (rows 30, 27, 29, 15, 19, 22 at ~15–45°) have the
  most. This is a *floor-clearance* margin, **not** a "demanded turn" — a low value is not "hard to
  fly," it just has less headroom before the #324 gate would reject it.
- **V∞ magnitude** — clusters ~0.9–2.3 km/s; moderate (~1.0–1.7) is the ballistic sweet spot.
  (Very low V∞, e.g. row 30 at 0.335 km/s, is *why* its achievable bend is huge — slow encounters
  bend easily — but such weak encounters are their own mission-utility question, distinct from the
  falsely-attributed "impossible periapsis.")
- **tof** — shortest-tof members accumulate least real-ephemeris drift over the V4-strict window;
  #312's own V4 failures cluster at the URA111 kernel edge and scale with integration span. This
  axis is unchanged by the correction and remains the primary V4-odds driver.
- **literature risk** — the 8 Titania-Oberon rows (23–30) sit in the most literature-exposed basin
  (§5); unchanged.

---

## 2. Corrected gauntlet representatives (supersedes #564 §3 table)

The cutoff **logic and count are unchanged** (Fable confirmed): one V4-strict-validated real-
ephemeris representative per previously-unrepresented pair, giving **5 candidates** which — with
the already-validated #312 (Umbriel-Oberon) — cover all six non-Miranda pairs. Only the
**Titania-Oberon pick** changes, and two of the other four **stated justifications** are corrected
(the picks themselves stand — see §2.1).

| priority | row | pair | n | n_rev | tof (d) | V∞ (km/s) | min bend° | why this one (corrected) |
|---|---|---|---|---|---|---|---|---|
| P1 | **23** | Titania-Oberon | 1 | (0,0) | 12.316 | 1.608 / 1.799 | **8.85** | **corrected T-O pick**: best margin of the three shortest-tof T-O rows (8.85 vs 24's 6.27, 25's 6.23), most moderate V∞ of that trio, NEW recovery. Carries the mandatory Canales/Kumar test (§5). |
| P1 | 1 | Ariel-Umbriel | 1 | (0,0) | 3.216 | 0.979 / 1.300 | 8.43 | shortest tof of all 30 → best V4 shot; lowest order. Margin 8.43° comfortable. Unchanged. |
| P1 | 5 | Ariel-Titania | 3 | (0,0) | 5.321 | 1.231 / 1.719 | 9.25 | shortest tof + lowest order for A-T (n=3 vs n=5). Margin 9.25° comfortable. *(Corrected: NOT "the only geometrically-clean A-T member" — row 6 actually has the larger margin, 11.24°; the pick rests on tof/order, not bend.)* |
| P1 | 7 | Ariel-Oberon | 5 | (0,0) | 7.752 | 1.521 / 1.829 | 6.22 | lowest order + shortest tof for A-O; NEW recovery. *(Corrected: NOT "lowest bend for A-O" — row 8 has the larger margin, 6.99°; the pick rests on lowest order + shortest tof. Margin 6.22° clears the floor comfortably.)* |
| P1 | 10 | Umbriel-Titania | 1 | (0,0) | 3.954 | 1.230 / 1.006 | 9.35 | shortest tof for U-T. Margin 9.35° comfortable. *(Corrected: #564's "bend 24.4° ≈ #312's" parenthetical was the max-bend value and is moot under the correct min-bend axis; the pick rests on shortest tof.)* |

**Recommended primary gauntlet set: 5 candidates (rows 23, 1, 5, 7, 10).** With #312 this validates
one real-ephemeris member in every non-Miranda pair — the complete, defensible basis for the
"family across the pairs" catalogue claim. (Count and logic unchanged from #564; only row 24 → 23.)

### Why row 23 over row 24 (and over the high-margin alternative row 27)

The T-O representative must (a) be a strong V4 candidate and (b) carry the literature test. Among
the **eight** T-O closures:

- The three **shortest-tof** T-O rows (all n=1, tof 12.316 d — half the n=2 rows' 24.632 d) are
  23 (min 8.85°, V∞ 1.61/1.80), 24 (min 6.27°, V∞ 1.97/2.16), 25 (min 6.23°, V∞ 1.26/2.17).
  **Row 23 dominates both 24 and 25 at the same tof**: larger margin *and* more moderate V∞.
  #564 picked row 24 believing it had the "cleanest bend"; under the correct axis row 24 is the
  **2nd-worst** T-O margin (only row 25 is lower). Row 23 is the correct short-tof pick.
- The **highest-margin** T-O member is row 27 (min 24.94°, V∞ 0.827/0.959) — essentially **#312's
  own validated V∞/bend profile transplanted to Titania-Oberon** (#312 is V∞ 0.926/0.965, min
  15.55°), so it is the closest dynamical analog to the already-validated member and the most
  robust-looking T-O candidate. Its cost is **2× the tof** (24.632 d), which by #564's own
  (unchallenged) tof→drift logic is the weaker V4 bet.

**Verdict**: **row 23 primary** — it keeps the short-tof / V4-drift advantage while being the
best-margin, most-moderate-V∞ member of that shortest-tof tier, and it is a NEW coverage recovery.
**Row 27 is the named high-margin alternative**: if the T-O slot is ever re-run for maximum
geometric robustness / closest-to-#312 profile and the 2× tof can be absorbed, row 27 is the pick.
The margin/tof trade is explicit — row 23 buys tof (V4 odds) at 8.85° margin; row 27 buys margin
(24.94°, #312-like) at double the tof.

### 2.1 The other four picks (rows 1/5/7/10) STAND — with two reasoning corrections

The task premise was that rows 1/5/7/10 "were justified on tof/order, not bend." Verified against
#564 §3's actual stated reasoning: **the picks all stand** (tof/order dominates and independently
confirms each), but **two of the four stated justifications did contain inverted bend sub-claims**,
now corrected in the table above:
- **Row 1** — clean: "shortest tof of all 30; lowest order." No bend claim. Stands as written.
- **Row 5** — the "only geometrically-clean A-T member (n=5 is 23.4° bend)" sub-claim is inverted
  (row 6 has the *larger* margin, 11.24° vs row 5's 9.25°). **Pick still stands** on shorter tof
  (5.321 vs 8.868 d) + lower order (n=3 vs n=5); reasoning restated to tof/order.
- **Row 7** — the "lowest bend for A-O" sub-claim is inverted (row 8 has the *larger* margin,
  6.99° vs row 7's 6.22°). **Pick still stands** on lowest order (n=5) + shortest tof (7.752 d);
  reasoning restated to tof/order.
- **Row 10** — the "bend 24.4°" parenthetical was the max-bend (not min) value and is moot under
  the corrected axis. **Pick stands** on shortest U-T tof (3.954 d); min bend 9.35° is comfortable.

---

## 3. Row 12 "best-geometry showcase" — WITHDRAWN

#564 §3 Tier-2 called **row 12** (Umbriel-Titania n=4 (0,0)) "the cleanest closure in the whole
family (max bend 5.59°)" and proposed it as a best-geometry showcase. **This is withdrawn.** Under
the correct axis row 12 is the **3rd-most-marginal** of all 30 (min 5.39°, barely above the 5.0°
floor) — the *opposite* of a best-geometry showcase; its low bend value means it sits *closest* to
the reject floor, not that it has the "cleanest" geometry.

If a robustness/best-geometry showcase is ever genuinely wanted, the correct picks are the
**highest-min-bend** members — row 30 (45.29°), row 27 (24.94°), or row 29 (17.77°) — **but all
three are in already-represented pairs** (Titania-Oberon), so per the unchanged "one representative
per unrepresented pair" census logic (§2) they are same-pair-redundant and add nothing to the
family-existence claim. Net: **no separate best-geometry showcase is recommended**; if intra-pair
robustness diversity is ever demonstrated, do it with a genuinely high-margin member, never row 12.

---

## 4. Rows 29/30 exclusion — RESTATED (redundancy, not "impossible periapsis")

#564 §3 excluded rows 29 (min 17.77°) and 30 (min 45.29°) "at any tier" as demanding
"near-surface/impossible periapses… the most fragile under real perturbations." **This exclusion
reason is false and withdrawn.** No required periapsis is computed anywhere; `max_bend_deg` is the
*achievable* deflection at the *safe* altitude, and rows 29/30 have the **two largest capability
margins in the entire family** — both larger than #312's own (15.55°). They are the *most* robust
members by this axis, not the least.

**Their valid exclusion reason is same-pair redundancy**: rows 29 and 30 are additional
Titania-Oberon closures, and the T-O pair is already carried by the selected representative (row 23,
§2). Validating a second/third member of an already-represented pair is evidentially redundant for
the census claim (it re-tests the same idealized→real gap in the same dynamical neighborhood). That
— **not** any feasibility or periapsis-altitude concern — is why they are outside the primary-5 set.

(The one legitimate *physical* caveat, distinct from the withdrawn one: row 30's V∞ is very low,
0.335/0.657 km/s — weak, slow encounters — which is a mission-utility question, not an "impossible
periapsis." It is *why* its achievable bend is 45°, and it is unrelated to the fragility claim #564
made.)

---

## 5. Literature risk / citation fix — Kumar arXiv ID corrected to 2509.03655

The substantive literature verdict (#564 §4) **stands unchanged and is confirmed**: the two closest
adjacencies to the Titania-Oberon direction — Canales-Howell-Fantino MMAT/FTLE (arXiv:2110.03683,
2308.10029; halo/manifold one-shot moon-to-moon transfer) and Kumar 2025 (single-moon Uranus-Oberon
MMR heteroclinics) — are both **structurally adjacent but neither a direct topology match**, so the
T-O representative (now row 23) must **cite both** as adjacent prior art but faces no novelty
collision. The other five pairs carry only the general #328 Uranian baseline.

**Citation ID fix**: #564 §4 cited Kumar 2025 as **arXiv:2509.12675**, which is a *different* paper
(Kumar-Rawat-Rosengren-Ross, Earth-Moon cislunar resonant transport — zero Uranian content). The
correct Uranus-Oberon MMR paper is **arXiv:2509.03655**, "Multi-shooting parameterization methods
for invariant manifolds and heteroclinics of 2-DOF Hamiltonian Poincaré maps," **§6.2** (MMR overlap
in the Uranus-Oberon PCRTBP), per `docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md:64,86`.
This is the same-author/same-year citation-collision trap flagged in
[[feedback_ground_citations_against_content]]. Any T-O writeback must use **arXiv:2509.03655 §6.2**.

---

## 6. Row-22 contradiction — RESOLVED (excluded on same-pair redundancy; NOT on bend)

#564 §3 stated two contradictory things about **row 22** (Umbriel-Oberon n=7 (2,2), the #312 n=7
sibling): (a) Tier-2 "already essentially known and cheap to confirm" (recommend), and (b) "NOT
recommended at any tier… Umbriel-Oberon row 22's partner-class **row-47** analog (bend 42°)…
demand near-surface flybys." The "row 47" does not exist (only 30 rows); statement (b) was in fact
describing **row 22 itself** (its *max* bend is 42.11°) under the inverted axis, and colliding with
statement (a).

**Resolution — one consistent statement:** Row 22's **min** bend is **15.43°** — a large, robust
margin, essentially equal to #312's own (15.55°); it is emphatically **not** a near-surface /
fragile flyby. But it is the **#312 n=7 sibling in the already-represented Umbriel-Oberon pair**, so
it is **excluded from the primary-5 on same-pair redundancy** (U-O is carried by #312), exactly like
rows 29/30 for T-O. The "row 47 / bend 42° / near-surface" language is an artifact of the inverted
axis and is dropped entirely. Row 22 remains the **natural first Tier-2 add** *if* intra-pair
multiplicity is ever demonstrated (it is robust, cheap, and the closest sibling to #312) — but it is
not part of the census-minimum set.

---

## 7. Corrected summary

1. **Bend axis inverted in #564; corrected here.** `max_bend_deg` is achievable-deflection margin
   above the 5° floor (higher = more robust), and the binding statistic is **min per-encounter
   bend**. #312 (row 19, 15.55°) sits at the *robust* end; the marginal members are rows 16/9/12/13
   (~5.1–5.5°, just above the floor).
2. **Gauntlet set unchanged in count/logic (5, not 30), one pick corrected:** **rows 23, 1, 5, 7,
   10** (was 24, 1, 5, 7, 10). **Titania-Oberon: row 24 → row 23** (best margin of the shortest-tof
   T-O trio, moderate V∞, NEW; row 24 was actually the 2nd-worst T-O margin). **Row 27** is the
   named high-margin / #312-profile alternative at 2× tof. Rows 1/5/7/10 **picks stand**; rows 5 and
   7 had inverted bend sub-claims in their justifications, now restated on tof/order grounds.
3. **Row 12 "best-geometry showcase" WITHDRAWN** — it is the 3rd-most-marginal member, not the
   cleanest.
4. **Rows 29/30 exclusion RESTATED** — same-pair redundancy with the T-O representative, **not**
   "impossible periapsis" (they are the two *most robust* members in the family, both above #312).
5. **Kumar citation fixed** to **arXiv:2509.03655 §6.2** (was the wrong 2509.12675). Substantive
   adjacency verdict unchanged.
6. **Row-22 contradiction RESOLVED** — excluded from the primary-5 on same-pair redundancy (U-O
   carried by #312); robust (min 15.43°), not fragile; the nonexistent "row 47 / bend 42°"
   near-surface language dropped.

**Unchanged and still standing from #564** (Fable-confirmed): §1 data verification, the 5-not-30
cutoff count/logic, the §5 #312 reframing prose (enumerated-30 idealized-exact kept distinct from
validated-1→6 real-ephemeris), and the §6 asymmetric-closures deferral (NO-GO). This note corrects
only the bend-axis ranking, the T-O pick, the Kumar ID, and the row-22 contradiction.
