# Self-seeding multi-rev extension + unsourced-row triage — RESULTS (task #177)

**Date:** 2026-06-09
**Code:** `src/cyclerfinder/search/self_seeding.py` (multi-rev Stage-A extension),
`tests/search/test_self_seeding.py` (mechanics + slow rescue test),
`scripts/triage_self_seeding.py` (cheap triage driver),
`scripts/validate_self_seeding_reachable.py` (reachable-subset full-tail validator).
**Runlogs:** `data/runs/self-seeding-triage-20260609T0205Z.jsonl` (212 rows),
`data/runs/self-seeding-reachable-20260609T0210Z.jsonl` (6 reachable rows).
**Prior:** `docs/notes/2026-06-08-self-seeding-results.md` (the S1L1 PASS),
`docs/notes/2026-06-08-self-seeding-one-member-results.md` (6.44Gg3 OFF-FAMILY #173).
**Writeback:** NONE. Recommendations only, held to main-session review.

---

## 0. What was built (the multi-rev Stage-A extension, #177 build 1)

#173's Stage A took only the SHORT-way (inbound) Mars radial crossing of the converged
G arc — one transit ToF. 6.44Gg3 was OFF-FAMILY partly because that single short-way
transit (~131 d) is half the row's long-transit signature (262 d). The extension
(`g_arc_branches`) enumerates ALL Mars-crossing branches of the SAME converged `(a, e)`
shape — the **short-way** (before aphelion), the **long-way** (after aphelion, longer
ToF), and **k full-revolution** branches (`k` in `[0, 2]`, reusing `core/lambert`'s
multi-rev solver for the Stage-B refinement) — so a long-transit row is matched against
ALL its branch shapes before being declared OFF-FAMILY. `triage_transit_match` gates
REACHABLE vs OFF-FAMILY on whether ANY branch transit lands within `tol_days` of the
row's tabulated transit.

Mechanics gate (green): 6.44Gg3's long-way branch (~292 d) reproduces its real-eph
262-d signature the short-way 131-d branch missed (`test_multirev_branch_rescues_
long_transit_6_44gg3`). The base short branch stays `g_arc_branches[0]` (byte-compatible
with #173's single `g_arc_shape`).

The gate tolerance is FIXED at **30 d** across all rows. It is CALIBRATED to the
validated S1L1 known-answer PASS (coplanar short-way 169 d vs tabulated 150 d, Δ +19 d)
and NOT widened per row to inflate REACHABLE (brief honesty rule).

---

## 1. Triage of the 212 unsourced rows (#177 build 2, the bulk)

Scope: the 12 non-V3 `russell-ch4` rows + all 200 `russell-ocampo` members (the
catalogue grew from the brief's 7+194 to 12+200). Cheap per-row, NO n-body — just the
coplanar-branch-vs-tabulated-transit comparison.

| outcome | count | meaning |
|---|---|---|
| **REACHABLE** | **6** | a coplanar G-arc branch transit lands within 30 d of the tabulated transit |
| OFF-FAMILY-NO-CLOSE | 2 | the descriptor shape does not reach Mars (`5.30gGf3`, `5.75ggF3`) |
| OFF-FAMILY-NO-DESCRIPTOR | 204 | no 2-arc g/G free-return descriptor to derive the coplanar G-arc shape from |
| **total** | **212** | |

The 6 REACHABLE (all `russell-ch4` rows that carry a 2-arc g/G descriptor):

| row | tabulated transit (d) | best branch | branch ToF (d) | Δ (d) |
|---|---|---|---|---|
| `russell-ch4-9.353Gg2` | 85 | short | 77 | −8 |
| `russell-ch4-9.94Gg3` | 82 | short | 91 | +9 |
| `russell-ch4-5.30ggF3` | 143 | short | 129 | −14 |
| `russell-ch4-3.78Gg3` | 171 | short | 188 | +17 |
| `russell-ch4-3.64gGg3` | 175 | short | 194 | +19 |
| `russell-ch4-6.44Gg3` | 262 | **long** | 292 | +30 |

6.44Gg3 is the multi-rev rescue: its short branch (134 d) is OFF-FAMILY vs 262 d, but
the **long** branch (292 d) lands within the gate. Without the #177 extension it would
have stayed OFF-FAMILY at the triage stage.

### Why 204 are OFF-FAMILY-NO-DESCRIPTOR (a real, important bound)

The `russell-ocampo` rows (Russell 2004 §3 Table 3.4) carry only an aphelion (often a
sub-Mars `aphelion_ratio < 1`), a simple-model transit, and v_inf anchors — NOT the g/G
free-return arc ToFs the coplanar G-arc shape must be derived from. The aphelion-ratio<1
members are sub-Mars: a single coplanar ellipse with `aphelion = a(1+e)` literally does
not reach Mars (they encounter Mars only near Mars's perihelion, an eccentric-Mars
geometry the coplanar shape cannot represent). They are recorded honestly as
inapplicable-gate negatives — NOT silently dropped, NOT counted REACHABLE. Reaching them
needs a descriptor-recovery / eccentric-Mars step out of this task's scope. This bounds
the catalogue's reach under the self-seed-multirev-v1 method and feeds the registry.

---

## 2. Validation of the REACHABLE subset (#177 build 3) — the decisive tail

For each of the 6 REACHABLE rows: full self-seed (best branch) → corrected-topology
longitude rendezvous → INDEPENDENT REBOUND/IAS15 Sun-only confirm (the #167 arbiter,
3-SOI band ≈ 0.0116 AU, NEVER loosened). The decisive scientific term is the EMERGED
v_inf vs the row's SOURCED anchor (the term that failed for 6.44Gg3 in #173).

| row | branch | v_inf E emerged / anchor | v_inf M emerged / anchor | n-body miss (AU) | n-body v_inf M | verdict |
|---|---|---|---|---|---|---|
| `9.353Gg2` | short | 12.45 / 9.35 | 21.01 / 10.52 | 3.6e-13 | 21.01 | OFF-FAMILY-AT-ANCHOR-VINF |
| `9.94Gg3` | short | 9.94 / 9.94 | 17.56 / 10.76 | 3.8e-13 | 17.56 | OFF-FAMILY-AT-ANCHOR-VINF |
| `5.30ggF3` | short | 6.33 / 5.30 | 11.90 / 5.44 | 7.7e-13 | 11.90 | OFF-FAMILY-AT-ANCHOR-VINF |
| `3.78Gg3` | short | 4.55 / 3.78 | 7.34 / 4.63 | 1.8e-12 | 7.34 | OFF-FAMILY-AT-ANCHOR-VINF |
| `3.64gGg3` | short | 3.96 / 3.64 | 6.26 / 4.59 | 1.3e-12 | 6.26 | OFF-FAMILY-AT-ANCHOR-VINF |
| `6.44Gg3` | **long** | 10.45 / 6.44 | **7.83** / 3.74 | 1.8e-12 | 7.83 | OFF-FAMILY-AT-ANCHOR-VINF |

### V3-candidate list: **EMPTY (0 of 6).**

Every REACHABLE row achieves the longitude rendezvous (n-body Mars miss at machine
precision, deep inside the unloosened 3-SOI band) but NONE reaches its sourced v_inf
anchors. The longitude/transit gate is necessary but not sufficient: the emerged v_inf
at the rendezvous epoch is off the family anchors in every case. This is the honest
minority-or-none yield the brief anticipated — a clean "mostly OFF-FAMILY, 0 reachable
to V3" is success, NOT a shortfall. No tolerance was loosened to manufacture a PASS.

PER-ROW, no batch-trust: each row was confirmed independently on the same-model
arbiter; one row's result never transferred to another.

---

## 3. The multi-rev extension's effect — did it rescue any rows?

**It improved 6.44Gg3 materially, but did not rescue it to V3-candidate.** The #173
short branch put the rendezvous at v_inf_M ≈ **10.9** km/s (3× the 3.74 anchor, above
the real-eph ~8 km/s ceiling). The #177 long branch (292 d, matching the 262-d
signature) drops it to v_inf_M ≈ **7.83** km/s — now UNDER the real-eph breathing
ceiling, a real step toward the family. But it is still off the 3.74 anchor, and the
departure v_inf rose to 10.45 (off the 6.44 anchor). So the long branch finds a
longer, lower-energy Mars arrival (the right qualitative family direction) without
closing the v_inf anchors. The slow test `test_multirev_long_branch_improves_6_44gg3_
mars_vinf` pins both halves: long-branch v_inf_M < short-branch v_inf_M (rescue
direction) AND long-branch v_inf_M still outside the 3.74 anchor band (honest residual).

No row was rescued from OFF-FAMILY to V3-candidate by the extension. Its value is (a)
admitting long-transit rows to the validation stage that the single-branch method would
have rejected at triage, and (b) quantifying — on the arbiter — how close the long
branch comes (6.44Gg3: within the real-eph ceiling, off the anchor).

---

## 4. Honest yield

- **REACHABLE: 6 / 212** (≈ 3%); the rest OFF-FAMILY (2 no-close, 204 no-descriptor).
- **V3-candidates: 0 / 6 reachable.** All six are OFF-FAMILY-AT-ANCHOR-VINF.
- The catalogue's reach under self-seed-multirev-v1 is bounded: only the 6 ch4 rows with
  a 2-arc g/G descriptor are even gateable, and none of those closes the v_inf anchors.

This matches the brief's stated honest expectation (a minority reachable; 6.44Gg3
failed) and the orbit-closure discipline (a clean negative is success; ALL binding
constraints — here the v_inf anchors, not just longitude/transit — must close).

---

## 5. Registry hand-off (#172)

The OFF-FAMILY negatives in both runlogs are method-versioned
(`self-seed-multirev-v1/coplanar-branch-transit-gate` for the triage,
`self-seed-multirev-v1/full-tail` for the validation). **#172's empty-region registry
(`data/empty_regions.jsonl`) should ingest these runlogs** as method-versioned empty
results — per the negative-results-registry memory, "empty" is conditional on the
method, and a later method that recovers ocampo descriptors or models eccentric-Mars
encounters would subsume this sweep.

---

## 6. Honesty ledger

- EXPECTED side = each row's SOURCED v_inf anchors + tabulated transit (Russell 2004
  Tables 3.4 / 4.9 / 4.13); emerged epoch / v_inf / miss are EVIDENCE.
- The 30-d transit gate is calibrated to the S1L1 known PASS and held FIXED; the
  1.5 km/s v_inf band is the #173 breathing band; the 3-SOI miss band is the #165/#167
  constant. NONE loosened to inflate REACHABLE or manufacture a V3-candidate.
- The independent REBOUND/IAS15 Sun-only integrator is the arbiter for every geometric
  arrival; the scientific verdict turns on the SOURCED v_inf anchor, reported not gated.
- No catalogue writeback. Nothing CONFIRMED on the search's say-so.
