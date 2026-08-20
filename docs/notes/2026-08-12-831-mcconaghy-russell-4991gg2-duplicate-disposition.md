# `#831` — `mcconaghy-2006-em-k2` / `russell-ch4-4.991gG2`: disposition

**Date:** 2026-08-12
**Registered by:** `#826` (`docs/notes/2026-08-11-826-russell12-closure-adjudication.md` §5)
**Verdict:** **KEEP BOTH ROWS.** No merge, no row deleted, no census count changed.
The two rows are one physical object carried as **two source realizations**, and that is a
lawful catalogue state — not a duplication defect.

---

## 1. The question

`#826` established that `mcconaghy-2006-em-k2` and `russell-ch4-4.991gG2` share, byte for
byte, `orbit_source: russell-2004-t49_413`, the `free_return_arcs` pair
`g(1.4612,526.02,Ll)` + `G(2.8096,651.46,U)`, `aphelion_au: 1.64`, `turn_ratio: 2.65`,
`period` k=2 / 4.27 yr and `sequence_canonical: "E-E-M-M"`. The McConaghy row's own
`orbit_elements.note` records Russell tagging Table 4.9 row 1 as parent cycler 4.991gG2 (#83)
and stating "Also known as the 'S1L1' cycler".

They are the same physical cycler. That was never in doubt. The question `#831` had to answer
is narrower: **does the catalogue's own deduplication rule require them to be one row?**

## 2. The spec answer: the merge rule does not fire

This is the load-bearing finding, and it is what makes keep-both *spec-conformant* rather
than merely convenient.

Spec §16.2 defines the canonical signature as `(normalised body sequence, period k, sorted
multiset of (body, V∞) pairs, sorted multiset of leg (a,e))`, with **V∞ binned to 0.05 km/s**
and (a,e) binned to 0.01 AU / 0.01. §16.3 then branches on it:

- **exact** (hash match) → log a re-derivation, *do not create a new entry*;
- **probable** → `probable-match-NEEDS-HUMAN`, and "the fuzzy stage is deliberately
  conservative — **it never auto-merges, only flags**";
- **novel** → the gauntlet.

Run the pair through it. Sequence, `k`, and the leg (a,e) multiset all match. The V∞ multiset
does **not**: 4.7 vs 4.99 is ~6 bins apart at 0.05 km/s, and 5.0 vs 5.10 is 2 bins apart. The
hashes differ, so **the "exact" branch never fires**. The pair is at most a
`probable-match-NEEDS-HUMAN` — and `#831` *is* that human confirmation step. Its answer is
recorded here.

Note also that `transit_times_days` (153 vs 150), the pair's other visible difference, is not
a signature input at all, so it neither joins nor splits them.

**Correcting a wishful claim in the rows themselves.** Both rows previously asserted that "M7
canonical-signature matching tolerances should be wide enough (~0.5 km/s on V_inf, ~10 d on
ToF) to collapse these two entries under one signature." That is not what the spec does: the
binning is 0.05 km/s, an order of magnitude tighter, and ToF is not in the signature. The
claim has been corrected in place on both rows rather than left to mislead a future reader
into thinking a merge is already mandated.

## 3. The affirmative case for keeping both

The spec argument shows a merge is not *required*. Three further grounds show it would be
actively harmful.

**(a) A merge would launder V0 evidence into a V3 row.** `#826` held
`mcconaghy-2006-em-k2` at **V0** on a sharp ground: `build_genome` reads `free_return_arcs` —
which are *Russell's* — and never touches the 4.7/5.0 fields, so the probe on the McConaghy
row **is** the `4.991gG2` probe (identical emerged V∞ 5.008 / 5.107 to three decimals, same
measured TR 2.658). It supplies zero independent evidence for this row, and the row's own
cited anchor is unreproduced: `vinf_source: mcconaghy-2006` gives 4.7 km/s at Earth against an
emerged 5.008 — a ~6.5% miss on a two-significant-figure published value.

A merge has only two possible outcomes, and both destroy that finding. Either the surviving
row absorbs the McConaghy anchor at V3 — validating by absorption exactly what `#826` refused
to launder — or the anchor is dropped, and with it the record that it is *unreproduced*.
**The V0/V3 split is information.** It says: the Russell circular-coplanar realization is
V3-validated; the McConaghy ephemeris-flavoured realization is not, and nobody has closed it.
One row cannot carry two validation levels.

**(b) A merge subtracts and adds nothing.** Every argument for merging is already satisfied
without one:

| merge would supply | already true today |
|---|---|
| the McConaghy citation on the Russell row | it is already in `russell-ch4-4.991gG2.corroborating_sources` |
| correct attribution under §16.4's earliest-`priority_date` rule | already on the Russell row: `2002-08-05` (McConaghy/Longuski/Byrnes, AIAA 2002-4420 — Russell's Ref. 15, the true first publication), earlier than the McConaghy row's `2003-08-01` |
| an explicit statement that they are one object | both rows already cross-reference each other by id in `notes` |

Meanwhile a merge would **delete** payload that exists only on the McConaghy row:
`dv_band: essentially_ballistic` + `dv_band_source: mcconaghy-2006` (SOURCED from the JSR
paper's ~10 m/s / ~30 yr statement), the 20-node `flyby_altitudes_km` block
(`flyby_altitudes_source: computed-m7`), its own `data_gaps` entries keyed to
DOI 10.2514/1.15215, and the 4.7 / 5.0 / 153-d realization itself. Information flows out of
the catalogue, none flows in.

**(c) There are genuinely ≥3 numerical realizations of this object, not 2.** Russell's own
Table 5.5 gives a *third* set for #83 4.991gG2 — a 7-cycle ephemeris optimisation with avg
V∞ E = 5.37, avg V∞ M = 5.48 km/s, avg E-M transit 165 d. The catalogue carries the two that
are separately sourced and separately anchored. A schema that cannot hold one object under
several published realizations would be the wrong schema for this literature; the row is the
right unit for *a sourced realization*, and `orbit_source` / `vinf_source` are exactly the
fields that say which.

## 4. Why not the reverse merge (McConaghy row as survivor)

Considered and rejected on every axis: the Russell row holds the earlier and more correct
`first_published` (2002 AIAA 2002-4420 vs the 2006 JSR paper), the earlier `priority_date`,
the higher validation level (V3 vs V0), the derived `loop-ee` segment with its
`source_quotes`, a populated `delta_v_kms: 0.0`, and the internally consistent anchor set —
its circular-coplanar `model_assumption` matches its circular-coplanar 4.99/5.10 values,
which is precisely what the McConaghy row cannot say of its own (see §6).

## 5. Census implications

**None.** No row is added or removed by this disposition; the frozen 240-id `MULTI_ARC_ALLOWLIST` in
`tests/data/test_cycler_class_census.py` is unchanged; `tests/data/test_validation_tier_census.py`
is unchanged. This is worth stating explicitly, because it is the *reason* `#826` deferred the
question rather than settling it inline: a merge would have rippled through every frozen
census ratchet plus **~184 references to `mcconaghy-2006-em-k2` across ~75 files**, including
live source (`src/cyclerfinder/verify/real_closure.py`), `data/negative_results.yaml`, and a
dedicated regression suite (`tests/verify/test_365_mcconaghy_2006_em_k2_v1.py`,
`scripts/run_365_mcconaghy_2006_em_k2_v1.py`). That blast radius is not the *reason* for the
verdict — §2 and §3 are — but it does confirm the cost/benefit is not close.

The one honest cost of keep-both: any naive "how many distinct physical cyclers do we hold"
count over rows is inflated by one for this pair. That is a property of counting rows rather
than signatures, and it is now documented at the rows themselves rather than left implicit.

## 6. Follow-up registered

**`#844`** — `mcconaghy-2006-em-k2` declares `vinf_fidelity: circular-coplanar` and
`model_assumption: circular-coplanar`, while its own `data_gaps` note states the 4.7/5.0
values are ephemeris-flavoured and that the circular-coplanar realization "cannot honestly be
claimed to close to 4.7 km/s (they close to 4.99, the sibling's anchor)". The row's declared
fidelity therefore contradicts its own sourced content. This is a more actionable defect than
the duplication ever was, and it is deliberately **not** folded in here: it is a fidelity-field
correction with its own evidence question (which anchor does McConaghy 2006 actually publish,
and under what model), not a disposition question.

Also noted, low priority: neither row carries `our_status`; per §16.4 both are
`known-reproduction`.

## 7. What changed

1. Both rows' `notes` corrected: the "M7 tolerances should be wide enough to collapse them"
   claim replaced with the §16.2 binning fact and the `#831` verdict, dated and referenced so
   the question is not re-litigated.
2. Three ratchet tests added to `tests/data/test_multi_arc_invariants.py` pinning the
   disposition mechanically — the shared physics must stay identical, the source anchors must
   stay distinct (including an explicit assertion that they do not collide under the §16.2
   0.05 km/s binning), and the V0/V3 split must not be laundered.
3. No catalogue row added, removed, or renumbered.
