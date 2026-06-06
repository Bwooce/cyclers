# Literature diligence: Forge SILVER novelty candidates vs Russell 2004

Task #126. Diligence on the Forge Phase 4 run's first two machine-confirmed,
bend-feasible, panel-survived SILVER candidates (held pending human review;
explicitly NON-catalogue). Question: does the literature — above all Russell's
exhaustive circular-coplanar enumeration — already document this family?

Primary source: **Russell 2004 dissertation (UT Austin, hdl:2152/1253)**.

## The two candidates (from `data/OUTSTANDING.md` "Forge Phases 4 + 5" section + the Forge plan completion notes)

Both are closed ballistic **2-synodic (k=2)** Earth-Mars cyclers on **real DE440
ephemeris**, topology **E-M-E-E** with **multi-rev Earth-Earth loop legs**, found
by `discover_novel` (construction-first, scan_parallel) and routed to SILVER:

| # | revs descriptor | Earth V∞ | Mars V∞ | other-encounter V∞ |
|---|---|---|---|---|
| 1 | (0,0,1) branches (s,s,low) | 9.75 | 13.01 | E 9.76, E 9.75 km/s |
| 2 | (0,0,2) | 9.62 | 12.06 | (rest not reported in the section) |

Per-encounter V∞ vector candidate 1: **[E 9.75, M 13.01, E 9.76, E 9.75] km/s**.
Both `match=novel` against the catalogue (supersession-aware, R1 delta 3). The
review_queue.jsonl artefact was not present on disk at diligence time (concurrent
agents are editing src/tests); candidate numerics are taken from OUTSTANDING.md
and the plan completion notes, which agree.

Diligence band per the task: Earth V∞ **9.4–10.0 km/s** and/or Mars V∞
**11.5–13.5 km/s**.

---

## 1. Russell rows in (or near) the V∞ band

### Russell's leg-descriptor taxonomy (needed for the topology comparison)

All Russell free-return legs are Earth-Earth returns. Three string types, each
beginning with a letter:

- **g / G** — *generic* return (G = the designated transit leg used to compute
  transit time and Mars V∞; lowercase = not the transit leg).
- **h / H** — *half-rev* return (the leading number is revs, e.g. `h(0.5,…)`).
- **f / F** — *full-rev* return; the `M:N` after the letter is the resonant
  revolution count, e.g. `f(1:1,…)`, `f(2:1,…)`, `f(3:2,…)`.

Cycler shorthand: `v∞E xx…x` where each x is the first letter of each leg, e.g.
`4.99gG`, `8.05gGf`. The `p.h.s.i` descriptor in Chapter 3 (e.g. `2.5.1.+0`) has
the leading integer **p = number of synodic periods** (so 2-synodic = `2.*`).

Model: **circular-coplanar Earth-Mars** throughout Chapters 3–4 (Chapter 5 walks
the same parent families into accurate ephemeris).

### Table 4.9 — Ballistic two-synodic period cyclers (doc p.127; PDF p.142) — THE direct comparison set

Russell's dedicated 2-synodic ballistic table. **Exactly five** 2-synodic
ballistic cyclers exist in his solution space; four are tabulated (the fifth is a
±-permutation variant). All rows, verbatim:

| V∞E | V∞M | tout/tin (d) | aphel (AU) | TR | legs |
|---|---|---|---|---|---|
| 4.99ᵃ | 5.10 | 150/150 | 1.64 | 2.65 | g(1.4612,526.02,Ll) G(2.8096,651.46,U) — the "S1L1" McConaghy cycler |
| 8.05 | 10.02 | 93/93 | 2.19 | 1.45 | g(1.4951,538.24,Ll) G(1.7757,279.24,U) f(1:1,74.468,-180.000) — "8.05gGf" |
| 8.16ᵇ | 10.06 | 92/92 | 2.20 | 1.08 | G(1.7708,277.48,U) f(1:1,74.244,143.198) h(0.5,0,U,±15.756) f(1:1,74.244,-36.802) |
| **9.35ᵃ** | **10.52** | 85/85 | 2.21 | 1.70 | G(1.7238,260.58,U) g(2.5469,916.9,L) — "9.35Gg" |

(ᵃ documented in Ref. 15 = McConaghy/Longuski/Byrnes; ᵇ documented in Chapter 3.)

**Highest Earth V∞ in the entire 2-synodic ballistic table = 9.35 km/s (Mars
10.52).** This is the single closest Russell 2-synodic row to our candidates — a
**near-miss, not a match**: Earth 9.35 vs our 9.62/9.75 (Δ 0.27–0.40), and Mars
10.52 vs our 12.06/13.01 (Δ 1.5–2.5 km/s — Mars V∞ is well outside the band on
the Russell side). No Russell 2-synodic row reaches Earth ≥ 9.4 or Mars ≥ 11.5.

Topology note: row 2 (`8.05gGf`) and row 3 (`8.165Gfh-f`) ARE the same E-M-E-E /
multi-arc-with-full-rev-loop topology class as our candidates (generic transit +
1:1 full-rev E-E loop). So Russell's taxonomy **does admit** our topology — but
his enumeration places it at 8.0–8.2 km/s Earth V∞, not 9.6–9.8, and with 1:1
(single-rev) loops, never the 2-rev loop of candidate 2.

### Table 3.4 — Two/three/four-synodic ballistic or near-ballistic (doc p.83; PDF p.98)

The 2-synodic (`2.*`) rows here: `2.1.1.+2` (E 4.1 / M 2.0), `2.3.1.+1`
(E 5.4 / M 5.3), `2.5.1.+0` (E 7.8 / M 9.9). None in band. The only Table 3.4
rows touching the band at all are higher-period: `4.7.1.-2` (E 6.6 / **M 11.4**,
just under), `4.14.1.+0` (E 14.1 / **M 12.7**), `4.10.1.-3` (E 10.2 / M 3.6) —
all 4-synodic, none a topology or V∞ match.

### Table 4.13 — Near-ballistic three-synodic (doc p.134; PDF p.149)

Verdict-relevant: Table 4.14 confirms **zero** near-ballistic 2-synodic cyclers
exist in Russell's space. The band-touching rows here are 3-synodic:
`7.51 / 12.70` and `7.51 / 10.00`. Not k=2; not a match.

### Tables 4.10–4.11 — Ballistic three-synodic, generic legs (doc pp.128–129)

3-synodic, not k=2, but recorded since they straddle the band. Rows in band:
`6.59 / 11.35`, `8.72 / 13.07`, `9.60 / 11.15`, `9.92 / 10.75`, `9.94 / 10.76`.
The `9.60/11.15` and `8.72/13.07` rows are the closest 3-synodic V∞ neighbours
to our candidates, but they are k=3 (wrong period class) — not a match.

### Chapter 3 Tables 3.9–3.11 — Five/six-synodic (doc pp.90–92)

5- and 6-synodic (k=5,6); recorded for completeness only. Band rows include
`5.10.2.+0` (E 8.6 / M 13.0), `5.11.1.-2` (E 9.5 / M 12.5), `5.12.1.+1`
(E 11.5 / M 13.2), `6.14.1.+1` (E 8.6 / M 13.1), `6.20.1.+0` (E 8.6 / M 13.0),
`6.21.1.-1` (E 12.7 / M 12.9). Wrong period class; not matches.

### Russell's exhaustiveness claim (Table 6.1, doc p.188; PDF p.203) — the decisive context

Russell's Chapter 4 "global claim" is explicitly bounded to a solution space:
**"Having v∞ at Earth less than 10 km/s"**, **"Total repeat times of 3 or less
synodic periods"**, transits on generic or full-rev legs, up to four generic
returns, all in the **circular-coplanar model**. Table 4.14 reports the complete
2-synodic census: **4 generic-ballistic + 1 full-rev = 5 ballistic, 0 full-rev
beyond that, 0 near-ballistic** for 2-synodic.

Our candidates (Earth V∞ 9.62/9.75 < 10; k=2 ≤ 3) sit **inside** Russell's
claimed-complete bounds. They are therefore NOT excluded by his stated scope —
yet his exhaustive 2-synodic enumeration does not contain them. The reconciling
difference is the **dynamical model**: Russell is circular-coplanar; our
candidates close on real **DE440** (eccentric, inclined) ephemeris. A real-eph
multi-arc chain can carry V∞ outside the circular-coplanar fixed points — exactly
the cross-fidelity boundary the Forge refuses to paper over (which is why the
single-ellipse resonance construction path was correctly demoted to *unavailable*
for these chains, per the plan notes).

---

## 2. Topology comparison

Our candidates: **E-M-E-E** with multi-rev Earth-Earth loop legs (1 rev for
candidate 1, 2 revs for candidate 2). In Russell's taxonomy these are a generic
transit leg (`G`) plus **full-rev** Earth-Earth loop legs (`f(1:1,…)` for 1 rev;
a 2-rev loop would be `f(2:2,…)` or a 2:1 resonant `f(2:1,…)`).

- His taxonomy **admits the topology**: Table 4.9 rows 2 and 3 (`8.05gGf`,
  `8.165Gfh-f`) are exactly generic-transit + 1:1 full-rev-loop 2-synodic
  cyclers, catalogued here as `russell-ch4-8.049gGf2` and `russell-ch4-8.165Gfh-f2`.
- But the matching descriptor family in his tables tops out at **Earth 8.16 /
  Mars 10.06 km/s for the f(1:1) loop case**, and Table 4.14 shows **zero**
  2-synodic full-rev cyclers beyond what's in Table 4.9. A **2-rev** (candidate 2)
  Earth-Earth loop has **no** 2-synodic representative anywhere in his tables.

So: topology class present in Russell; the specific high-V∞ (and 2-rev-loop)
members are not. Consistent with the model-boundary explanation above.

---

## 3. Gapped catalogue rows that could "match by gap"

Checked `data/catalogue.yaml` (237 entries). Every Russell Table 4.9 2-synodic
cycler is already individually catalogued with V∞ populated:
`russell-ch4-4.991gG2` (4.99/5.10), `russell-ch4-8.049gGf2` (8.05/10.02),
`russell-ch4-8.165Gfh-f2` (8.16/10.06), `russell-ch4-9.353Gg2` (9.35/10.52),
plus the Table 3.4 `russell-ocampo-*` rows. None of these has a V∞ gap, and the
highest (9.35/10.52) is still a near-miss, not a match.

k=2 Earth-Mars rows with **null/gapped V∞** (the "failed only because of a gap"
candidates):

| row id | k | seq | V∞ status | could it be one of ours? |
|---|---|---|---|---|
| `mcconaghy-2005-em-case1` | 2 | E-E-M-M | [null, null] | **Unlikely.** It is the S2L1 "Case 1" ballistic 2-synodic cycler; Rogers 2012 Table 3 pins the 4:3(2)- establishment flyby V∞ at **3.378 km/s** and the family sits with Cases 2/3 at ~2.0–5.3 km/s. Low-energy SnLm family, ~6–7 km/s below our candidates. Gap is on the steady-state value + return-leg ToF, not on the energy regime. Not a plausible host. |
| `russell-2004-ch4-broadclass-family` | null | (seed) | [] | **No.** Family-seed for the Chapter-4 tables; its members ARE the populated Table 4.9 rows above (max 9.35/10.52). Carries no hidden high-V∞ member. |
| `mcconaghy-2005-em-snlm-broadclass-family` | null | E-M | [null, null] | **Unlikely.** SnLm taxonomy seed; per-member V∞ spans ~2.5 (Case 2) to ~12.6 km/s (U0L1 *establishment*, an up-escalator single-leg branch, not a 2-synodic E-M-E-E multi-loop). The 12.6 figure is an establishment-leg flyby V∞, wrong topology and not a steady-state 2-synodic cycler. No 2-synodic E-M-E-E member at ~9.7/13. |

Other gapped EM rows (`aldrin-4-3-2-establishment` k=2 Earth 5.509; the niehoff/
arenstorf/genova/wittal Earth-Moon rows; `mcconaghy-2005-em-u0l1` k=1) are out of
period class, body set, or energy regime. **No gapped catalogue row is a
plausible match for either candidate.**

---

## 4. VERDICT per candidate

### Candidate 1 — (0,0,1) (s,s,low), V∞ [E 9.75, M 13.01, E 9.76, E 9.75] km/s

**NOT-FOUND (novelty case strengthened).**

- Closest Russell row: Table 4.9 `9.35Gg` (E 9.35 / M 10.52) — near-miss on
  Earth (Δ 0.40), clear miss on Mars (Δ 2.5; Russell's Mars V∞ is below the
  11.5–13.5 band entirely). Topology class (generic transit + 1:1 full-rev loop)
  IS present in Russell at rows 2–3, but at 8.0–8.2 km/s Earth V∞, not 9.75.
- Russell's complete 2-synodic census (Table 4.14: 5 ballistic, 0 near) contains
  no member at this energy. The candidate falls inside Russell's stated global
  bounds (E V∞ < 10, k ≤ 3) yet is absent — explained by Russell being
  circular-coplanar vs our DE440. Not a rediscovery.
- No gapped catalogue row hosts it.

### Candidate 2 — (0,0,2), V∞ [E 9.62, M 12.06, …] km/s

**NOT-FOUND (novelty case strengthened).**

- Same near-miss neighbour (Table 4.9 `9.35Gg`, E 9.35 / M 10.52). Additionally,
  candidate 2's **2-rev** Earth-Earth loop has *no* 2-synodic representative
  anywhere in Russell's tables (Table 4.14: zero 2-synodic full-rev cyclers
  beyond the single f(1:1) member already in Table 4.9). Mars 12.06 is outside
  Russell's 2-synodic range.
- No gapped catalogue row hosts it.

Both verdicts are conservative: reported as NOT-FOUND-in-the-literature with the
near-miss numbers explicit. This corroborates the Forge's SILVER hold —
machine-confirmed, unsourced, model-distinct from the circular-coplanar canon —
and does NOT by itself assert a publishable novel discovery (an independent
*source* would be needed for GOLD, and none exists). The most likely human
resolution remains the one recorded in OUTSTANDING.md: a different / higher-energy
(possibly real-eccentricity-enabled) family than the sourced near-ballistic
Sanchez/S1L1 cyclers.

---

## Coverage statement (honest)

**Checked (text-layer extraction, `-layout`, verified verbatim):**
- Russell Table of Contents / List of Tables (front matter) — table index + captions.
- Table 3.4 (doc p.83 / PDF p.98) — full, all 44 rows scanned for the band.
- Table 4.9 (doc p.127 / PDF p.142) — full, the dedicated 2-synodic table.
- Tables 4.10–4.11 (doc pp.128–129 / PDF pp.143–144) — full, band rows extracted.
- Table 4.13 (doc p.134 / PDF p.149) — near-ballistic 3-synodic.
- Table 4.14 (doc p.134 / PDF p.149) — summary census (the exhaustiveness numbers).
- Tables 3.9–3.11 (doc pp.90–92 / PDF pp.105–107) — 5/6-synodic, band-filtered scan.
- Table 6.1 (doc p.188 / PDF p.203) + §6.2/§6.3 — the global-claim bounds + model.
- Chapter 5 ephemeris discussion + Table 5.5/5.6 captions (PDF pp.192–195) —
  confirmed Chapter 5 re-optimizes the SAME Chapter-4 parent families into
  accurate ephemeris (∆v / averaged V∞), introduces no new high-V∞ 2-synodic family.
- Leg-descriptor taxonomy text (PDF pp.141–142) for g/G, h/H, f/F semantics.
- `data/catalogue.yaml` — all 237 entries parsed; all Earth-Mars k=2 and all
  gapped Earth-Mars rows inspected.

**NOT checked / limitations:**
- Did not read the *full* row-by-row contents of Tables 4.10–4.11 (141 rows) or
  3.9–3.11 beyond a column-5/6 band filter — only band-touching rows were pulled
  verbatim. (These are 3/5/6-synodic, off our k=2 class, so low risk.)
- Did not exhaustively read Chapter 5's per-cycle ephemeris solution tables
  (Table 5.5/5.6 bodies) row-by-row — confirmed only that they re-optimize
  Chapter-4 parents, not that no single averaged-V∞ entry lands near 9.7/13.
- Ref. 25 (the external cycler-nomenclature source Russell cites for the
  descriptor strings) was not consulted; descriptor semantics taken from
  Russell's own in-text summary.
- The live `data/review_queue.jsonl` SILVER artefact was absent at diligence time
  (concurrent agents editing); candidate numerics rely on OUTSTANDING.md + the
  Forge plan completion notes (which are mutually consistent).
- Text extraction was clean throughout; the PDF vision (pdftoppm) fallback was
  not required.
