# Schema-v3 `legs[]` → `trajectory.segments` Backfill Inventory & Worklist

Task #66 — the sweep/inventory step (spec §16.6.5 steps 3–4) of the **lazy,
sweep-driven** backfill. This is the actionable worklist that ranks the
remaining legacy entries by migration-readiness. It is NOT a mandate to
bulk-rewrite: per `data/README.md` ("Backfill") a big-bang rewrite "would write
mostly nulls and obscure which nulls are genuine gaps."

Regenerate the underlying numbers with:

```sh
uv run --with pyyaml python scripts/_backfill_sweep.py
```

(`scripts/_backfill_sweep.py` is a throwaway analysis script; the authoritative
machine query is `cyclerfinder.data.catalog.find_data_gaps()` / the
`python -m cyclerfinder.data.catalog gaps` CLI.)

## Headline counts (as of 2026-06-01, post task #66 second pass)

| Class                       | Count | Meaning                                                              |
| --------------------------- | ----: | ------------------------------------------------------------------- |
| **migrated**                |    20 | has `trajectory.segments` (see "Migrated this pass" below)          |
| **legacy** (`legs[]` only)  |   186 | still on the flat top-level `legs[]`                                |
| **neither** (citation-only) |    13 | family-seed / citation-only; no `legs[]`, no segments → **skip**    |
| **total**                   |   219 |                                                                     |

### Legacy breakdown by migration-readiness

| Rank             | Count | Definition                                                                                          | Action                                                   |
| ---------------- | ----: | --------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| **full**         |     0 | legs cover the encounter topology, every leg has `tof_days`, a single orbit-level `(a,e)` applies   | migrate losslessly (all such entries already done)       |
| **topology-gap** |   185 | every present leg has `tof_days`, but `len(legs) < n_encounters` from `sequence_canonical`          | needs s1l1-style `data_gaps[]` markers, NOT a naive copy |
| **partial**      |     0 | some legs lack `tof_days`                                                                            | —                                                        |
| **sparse**       |     1 | no legs, or no per-leg `tof_days` (`aldrin-4-3-2-establishment`)                                     | low priority; source-gated                               |

**Key finding (confirms the README's lazy-backfill rationale):** after the two
lossless Aldrin migrations below, **zero** legacy entries are "full." All 185
remaining legacy entries are *topology-gap* — their `legs[]` records only the
source-tabulated transit(s) (typically the single E→M leg, or an E→M / M→E
symmetric pair), while `sequence_canonical` (e.g. `"E-E-M-M"`, 4 encounters)
implies an intermediate Earth-loop + return that the literature tables do not
break out. A mechanical rewrite would invent null return/loop segments,
defeating the whole point of `data_gaps[]`. Each of these needs the same
deliberate, source-gated treatment as the `s1l1-2syn-em-cpom` exemplar.

## Migrated this pass (task #66)

### First pass — lossless Aldrin pair

| id                              | why it qualified                                                                                                                                      |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `aldrin-classic-em-k1-outbound` | 1-synodic `E-M` cycle = exactly 2 encounters; both legs (146 d out / 634 d return) source-attested; single Aldrin ellipse `a=1.60, e=0.393` on both arcs. Lossless — no `data_gaps[]`. |
| `aldrin-classic-em-k1-inbound`  | Time-mirror of the outbound; same lossless conditions (`M-E`, 146 d / 634 d, same ellipse).                                                            |

### Second pass — attested-segments + `data_gaps[]` (the s1l1 exemplar shape)

All 17 below follow the `s1l1-2syn-em-cpom` pattern: attested transit leg(s)
populated with real values, every missing return / intermediate-loop segment
and missing per-segment `(a,e)` registered as an explicit `data_gaps[]`
known-unknown. **None are lossless.**

| id                          | shape migrated                                                                                                  |
| --------------------------- | --------------------------------------------------------------------------------------------------------------- |
| `mcconaghy-2006-em-k2`      | `E-E-M-M`; 2 attested 153-d legs (out/ret, abstract) + 1 `loop-ee` (derive) + per-segment `(a,e)` gaps          |
| `mcconaghy-2005-em-case1`   | `E-E-M-M`; 1 attested E→M leg (365 d, Rogers Table 4) + `ret-me` (unknown) + `loop-ee` (derive); orbit `(a,e)` 1.22/0.238 on attested arc |
| `mcconaghy-2005-em-u0l1`    | `E-M` (L=0, **no loop**); 1 attested E→M leg (97 d, Rogers Table 4) + `ret-me` (unknown); orbit `(a,e)` 2.05/0.563 on attested arc |
| `russell-ch4-*` (14 rows)   | 2 attested legs (`t_out`/`t_in`, Russell 2004 Table 4.x) + 1 (`k=2`) or 2 (`k=3`) `loop-ee*` (derive) + per-segment `(a,e)` gaps; orbit `(a,e)` null |

All 20 migrated entries are verified **signature-preserving**: every
`signature_hash` is bitwise-identical before/after (the V∞ multiset and
orbit-level `(a,e)` — the only signature inputs — are untouched, and the loader
dedupes leg `(a,e)` to the same multiset regardless of leg count), and the
aggregate sha1 over the sorted list of all non-null entry signature_hashes is
unchanged (`25d8a711…`). The loader's `_segments_as_legs` reads the new
`trajectory.segments` and produces the identical canonical signature it
previously derived from `legs[]`.

## Ranked worklist — remaining legacy entries (185 topology-gap + 1 sparse)

All remaining migrations are **source-gated** and must follow the
`s1l1-2syn-em-cpom` pattern: record the attested segment(s) with real values,
and register every missing return / intermediate-loop segment (and any missing
per-segment `(a,e)`) as an explicit `data_gaps[]` known-unknown. Do **not**
write null segments without a matching `data_gaps[]` entry.

Priority order (highest-value first):

1. ~~**`mcconaghy-2006-em-k2`**~~ — **DONE** (task #66 second pass). Migrated as
   2 attested 153-d segments + `loop-ee` (derive) + per-segment `(a,e)` gaps,
   the exemplar shape. The leg decomposition still depends on AIAA 2002-4420 /
   JSR 2006 for the loop arc and per-segment elements — those are the registered
   `data_gaps[]`.
2. ~~**`russell-ch4-*` family (14 entries)**~~ — **DONE** (task #66 second pass).
   Each migrated as 2 attested `t_out`/`t_in` segments + intermediate-loop
   segment(s) (1 for `k=2` `E-E-M-M`, 2 for `k=3` `E-E-E-M-M`) marked `derive`,
   plus per-segment `(a,e)` `unknown` gaps. orbit-level `(a,e)` left null.
3. ~~**`mcconaghy-2005-em-case1`, `mcconaghy-2005-em-u0l1`**~~ — **DONE** (task
   #66 second pass). case1 (`E-E-M-M`): 1 attested 365-d E→M segment + `ret-me`
   (unknown) + `loop-ee` (derive). u0l1 (`E-M`, L=0): 1 attested 97-d E→M segment
   + `ret-me` (unknown), **no loop segment** (zero intermediate loops). Both
   inherit the Rogers Table 1 cycler ellipse on the attested arc(s).
4. **`russell-ocampo-*` family (~185 entries, the bulk)** — single E→M transit
   leg from Russell 2004 Tables 3.4 / 4.x; `E-E-M-M` (or higher-k) topology with
   `orbit_elements.a_au/e = null` (only Aphelion Ratio tabulated). These are the
   lowest-yield: migrating each writes 1 real segment + 3+ `data_gaps[]`
   markers. **Still DEFERRED** until the multi-rev Lambert solver (task #54) can
   *derive* the missing loop arcs, at which point they become `kind: "derive"`
   gaps that get filled programmatically rather than from literature. These are
   the bulk of the remaining 185 legacy entries.
5. **`aldrin-4-3-2-establishment`** (sparse; `analytic-ephemeris`) — 1 leg with
   no `tof_days`. **Source-blocked** on Rogers 2012 Tables 3/4 extraction (also
   blocks `v_infinity_leveraging_dv_kms`). Lowest priority; not migrated.

## "Done" criterion

Per spec §16.6.5: the backfill is complete when **no entry retains a top-level
`legs[]`**. Tracked continuously by `find_data_gaps()`, which emits a synthetic
`{path: "trajectory.segments", kind: "derive"}` gap for every un-migrated
legacy row. Today that is 186 entries (down from 203 after the task #66 second
pass; the bulk `russell-ocampo-*` family remains, deferred to task #54).
