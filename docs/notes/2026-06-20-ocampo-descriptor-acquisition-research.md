# Ocampo n.m.k descriptor acquisition — are the per-arc free-return descriptors recoverable?

**Date:** 2026-06-20 AET. **Type:** read-only acquisition research. **Writeback:** NONE.
**Question:** Can the per-arc free-return descriptors (the g/G/f/h arc geometry) — or
per-member orbital state — for the ~200 `russell-ocampo-*` `n.m.k` rows be RECOVERED
from the published Russell/McConaghy corpus, or is the publication gap genuine?

## Verdict

**NO — genuine publication gap for the bulk.** The per-arc free-return descriptor is
**recoverable for at most ~7 specific cyclers** (the ones individually worked-out in a
detail table), of which **only 2 map to catalogue `russell-ocampo-*` rows**. For the
remaining ~190+ `n.m.k` rows the per-arc geometry is **not printed in any source we
hold, nor in any source the published corpus points to** — it lives only in Russell's
unpublished electronic file ("contact the author"). So:

- **Descriptor-acquisition cannot un-gate the bulk.** The data does not exist in print.
- The realistic descriptor-recoverable set is **~7 cyclers total**, **2 of which are
  catalogue `russell-ocampo-*` rows** (`2.5.1+0`, `4.3.1-5`) — and both already have a
  promotion path that does NOT need the descriptor (see below). So the *net new* rows
  un-gated by descriptor acquisition is effectively **zero**.

A separate, larger un-gating lever exists but is **not descriptor acquisition**: the
single-ellipse `closer_sweep_v1` substrate (#365) already promoted 7 ocampo rows to V1
by treating the high-V∞ short-transit members as Lambert-equivalent radial-crossing
ellipses. That is a *solver* path, explicitly out of scope for this "is the data
recoverable" question, and is bounded by topology (only single-ellipse-reachable
members close), not by missing descriptors.

## What a row needs to be descriptor-seedable (the contrast)

A descriptor-bearing row carries `free_return_arcs` with `arc_type` / `resonance` /
`tof_years` / `raw_descriptor`. Concrete example — `mcconaghy-2006-em-k2` (the in-scope
control), verbatim from `data/catalogue.yaml`:

```
free_return_arcs:
  - arc_type: generic
    tof_years: 1.4612
    raw_descriptor: "g(1.4612,526.02,Ll)"
  - arc_type: generic
    tof_years: 2.8096
    raw_descriptor: "G(2.8096,651.46,U)"
```

That `g(t_f, θ, ε)` form is exactly the McConaghy-Russell-Longuski 2005 per-leg formal
descriptor (g/f/h leg grammar, §digest below). A multi-arc closer can seed each arc from
`(t_f, θ, ε)`. **The 200 `russell-ocampo-*` rows carry `free_return_arcs: absent/[]`**
— confirmed by a YAML scan: `0 of 200` ocampo rows have any `free_return_arcs`, whereas
`11 of 14` `russell-ch4-*` rows do. The ch4 rows carry descriptors because Russell's
Chapter-4/Appendix-C uses the g/G/f/h shorthand per parent cycler; the Chapter-3 `n.m.k`
member tables do not.

## Evidence — what each source actually prints

Checked against the digests/transcriptions (all sourced from the private corpus
`cyclers_pdf/papers`, vision-read, page-cited) and the catalogue:

| Source (corpus file / table) | What it prints for `n.m.k` cyclers | Per-arc descriptor? |
|---|---|---|
| **Russell 2004 dissertation Tables 3.4 / 3.9 / 3.10 / 3.11** (the 201-row source of the `russell-ocampo-*` rows; `russell-2004-t34` / `t39_311`) | Per row: Aphelion Ratio, Turn Ratio, E→M transit (days), V∞ at Earth, V∞ at Mars, per-flyby geocentric turn angles. **Cycle-level summary only.** | **NO.** No (a,e,ω,ν), no per-arc (t_f, θ, ε). `2026-06-07-russell-2004-member-tables-transcription.md` transcribes all four tables verbatim — none has a per-arc column. |
| **Russell 2004 dissertation Tables 3.5–3.8** | Per-encounter Δv vectors (Δvx,Δvy,Δvz km/s) + times + initial r_mars, for **4 specific cyclers**: `2.5.1+0`, `3.1.2+1`, `4.3.1-5`, `4.5.2-2`. Enough to simulate one cycle. | Partial: per-encounter *state*, not the g/f/h descriptor — but reproducible-grade for those 4 only. |
| **Russell-Ocampo 2003 (AAS-03-145 / JGCD 27(3))** Tables 4 / A1 / A2 | Same columns as dissertation T3.4 (AR/TR/transit/V∞/turns). | **NO** per-arc geometry. |
| **Russell-Ocampo 2003 Tables 5–8** | Same per-encounter Δv state as dissertation T3.5–3.8, same 4 cyclers (`2-5-1-3`,`3-1-2-11`,`4-3-1-20`,`4-5-2-12`). | Partial state for the same 4 only. |
| **McConaghy-Russell-Longuski 2005 (JSR 42(4)) Table 2** | The formal per-leg g/f/h descriptor labels — the EXACT descriptor format. But **only 7 "well-known" cyclers**: Aldrin, VISIT-1, VISIT-2, ballistic S1L1, Byrnes case-3, **Cycler 2.5.1.+0**, **Cycler 4.3.1.−5**. | **YES — but for 7 cyclers only.** (digest `2026-06-17-digest-mcconaghy-2005.md` §3 transcribes all 7.) |
| **McConaghy-Longuski-Byrnes 2004 (JSR 41(4))** Tables 1-3 / 4 / 5 / 6 | Tables 1-3 = graphical drawings (no numbers). Table 4 = 21 promising cyclers (AR/TR/V∞/transit). Table 6 = S1L1 DE405 22-encounter itinerary (one cycler). | **NO** broad per-arc descriptors; only the single S1L1 itinerary is reproducible-grade. |
| **McConaghy 2006 (JSR 43(2))** | Two-synodic S2/S1L1 focus; maps entirely to the one `mcconaghy-2006-em-k2` row. | One cycler. |

**Bottom line of the table:** the per-arc descriptor (g/f/h) is printed in exactly ONE
place — McConaghy 2005 Table 2 — for exactly **7 named cyclers**. Everywhere else the
`n.m.k` family is published as a cycle-level summary (AR/TR/V∞/transit/turn-angles).
This is the descriptor-gating the task describes, and it is **real**.

## Which rows are descriptor-recoverable, and what would be ingested

Of the 7 McConaghy-2005-Table-2 cyclers, only 2 map to catalogue `russell-ocampo-*` IDs:

1. **`russell-ocampo-2.5.1+0`** ← `2g(1-11/14, 11/14 rev, U) f(1:1, 74.919°, ∓144.069°)
   h(0.5, 0, U, ±15.081°) f(1:1, 74.919°, ±35.931°)` (4 arcs: g,f,h,f). No free vars.
2. **`russell-ocampo-4.3.1-5`** ← `4g(7-1/14, 5-1/14 rev, L) f(1:1, 84.039°, ∓90.0°)
   h(0.5, 0, U, ±5.961°)` (3 arcs: g,f,h). No free vars.

These two could have `free_return_arcs` ingested from McConaghy 2005 Table 2 (the
descriptor format matches the schema exactly). **However, both already have a non-
descriptor promotion record:** `russell-ocampo-2.5.1+0` is **already V1** (#365, single-
ellipse closer; it is one of the 7 currently-V1 ocampo rows), and `russell-ocampo-4.3.1-5`
is a documented #365 NO-CLOSE multi-arc negative. So ingesting the descriptor would
*enable a multi-arc closure attempt* on `4.3.1-5` (genuine new capability for that one
row) but would not by itself promote anything — it would feed the multi-arc closer lane,
whose closure is a separate solver question, not an acquisition question.

The other 5 Table-2 descriptors (Aldrin = already V3; VISIT-1/2 = have a FREE parameter
ϕ, so NOT fully reproducible without a Niehoff source; Byrnes case-3 = has free λ) are
either already covered or not cleanly catalogue-`russell-ocampo` rows.

The 4 Δv-state cyclers (T3.5–3.8 / RO T5–8) overlap heavily with the above
(`2.5.1+0`, `3.1.2+1`→`russell-ocampo-3.1.2+1`, `4.3.1-5`, `4.5.2-2`→`russell-ocampo-4.5.2-2`)
and were all run by #365: 1 PASS (`2.5.1+0`), 3 honest multi-arc negatives. They are
state-reproducible but their state is already mined; no descriptor acquisition needed/possible.

## The acquisition target that WOULD broaden recovery (and why it probably won't)

The only published pointer to broader per-member geometry is:

- **McConaghy, Yam, Landau & Longuski, "Two-Synodic-Period Earth-Mars Cyclers with
  Intermediate Earth Encounter," AAS 03-509** (Aug 2003) — ref 14 in McConaghy 2005,
  cited as the S1L1 primary source. **NOT in the corpus** (`cyclers_pdf/papers` scan:
  not found). Also **AAS 03-508** (Russell-Ocampo, ref 13) — not in corpus.
  - Acquisition value: these are 2-synodic-focused; even if acquired they would at most
    add descriptor detail for the S1L1/S2 family already covered by `mcconaghy-2006-em-k2`
    and `russell-ch4-4.991gG2` (V3). They are **unlikely to print the per-arc geometry
    for the broad ~190-member `n.m.k` set** — that set is a 2502-cycler enumeration whose
    per-member geometry Russell explicitly keeps in an unpublished electronic file.
- **Russell 2004 dissertation Appendix C** explicitly states (p.201, per
  `…member-tables-transcription.md` §App.C): an electronic file of all 4263 ephemeris
  solutions is available **"contact the author"** — i.e. the canonical machine-readable
  per-member source is NOT in any PDF. This is the real location of the missing data.

So the highest-leverage acquisition is not a paper at all but **R. P. Russell's archival
electronic file** (author contact), which is outside the published-corpus mining lane.

## Honest bottom line

- **Realistically un-gateable by descriptor ACQUISITION: ~0 net new rows.** The per-arc
  descriptor is published for only 7 cyclers (McConaghy 2005 Table 2); only 2 are
  `russell-ocampo-*` rows; one of those is already V1 and the other already has its
  Δv-state mined and a #365 multi-arc negative on record. Acquiring the
  `4.3.1-5` descriptor would enable a *multi-arc closure attempt* on that single row
  (capability gain), not an automatic promotion.
- **The ~190 `n.m.k` rows remain V0 by a genuine publication gap.** Every printed source
  gives them as a cycle-level summary (AR/TR/V∞/transit/turn angles). No table or appendix
  in the held corpus, or in the corpus's cited siblings, prints their per-arc free-return
  geometry or per-member orbital state. The data exists only in Russell's unpublished
  electronic file.
- **Specific acquisition required to broaden recovery:** R. P. Russell's electronic
  4263-solution file (dissertation Appendix C, "contact the author"). Secondary,
  lower-confidence: AAS 03-509 and AAS 03-508 (neither in corpus; both 2-synodic-scoped,
  so unlikely to cover the broad member set).
- This **confirms** the `project_validation_ceiling` memory: the ocampo V0 wall is a
  PUBLICATION gap (n.m.k summary-only), not a solver or laziness gap. Past the ceiling is
  new-input-gated — and the specific new input is an author-held data file, not a paper.

## Sources read for this note
- `docs/notes/2026-06-07-russell-2004-member-tables-transcription.md` (Tables 3.4–3.11 + 3.5–3.8 + App.C verbatim)
- `docs/notes/2026-06-17-digest-russell-ocampo-2003.md` (Tables 4/5-8/A1/A2)
- `docs/notes/2026-06-17-digest-mcconaghy-2005.md` (Table 2 — the only per-arc descriptor table; §3 transcription)
- `docs/notes/2026-06-17-digest-mcconaghy-2004.md` (Tables 4/6)
- `docs/notes/2026-06-08-self-seeding-triage-results.md` (#177: 204 OFF-FAMILY-NO-DESCRIPTOR — the gating mechanism)
- `docs/notes/2026-06-17-365-phase-d-results.md` (#365 promotion wave; the single-ellipse closer path)
- `docs/notes/2026-06-07-russell-2004-dissertation-method-mining.md`; `docs/notes/CORPUS_INDEX.md`
- `data/catalogue.yaml` (`mcconaghy-2006-em-k2` descriptor control vs `russell-ocampo-*` rows; YAML scan of free_return_arcs + validation levels)
- Corpus listing `cyclers_pdf/papers/` (confirmed AAS 03-509 / 03-508 absent)
