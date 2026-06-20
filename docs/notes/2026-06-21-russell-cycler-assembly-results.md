# #388 — Russell ψ generic-return generator + global cycler search: results

Date: 2026-06-21. Outcome of building Russell's (2004) generic-return generator,
turn-angle cycler assembly, and global `p.h.s.i` search (spec
`2026-06-21-russell-psi-generic-return-generator-design.md`, plan
`2026-06-21-russell-psi-generic-return-generator.md`), then applying it to the
catalogue's descriptor-bearing rows. **No catalogue writeback; promotions held.**

## What was built (all golden-validated)

| phase | module | golden | status |
|---|---|---|---|
| A | `search/generic_return.py` | Russell Tables 2.2 & 2.3 (generic returns) | PASS (`cb262f1`) |
| B | `search/cycler_assembly.py` | Cycler 4.9.2.-1 grouping/flyby counts | PASS (`20c56e2`) |
| C | `search/cycler_search.py` | Russell Table 3.4 — recovers 4.3.1.-5 (3.098 km/s, AR 0.9926) + Aldrin 2.1.1.+2 | PASS (`bee2ed1`) |
| D | descriptor→p.h.s.i + batch | catalogue-row assembly | this note (`b20e598`, `6d20fd3`) |

The generator reproduces Russell's circular-coplanar generic returns and **rediscovers
his published 44-cycler global-search table** from first principles (multi-rev Lambert
grid → (N, fast/slow) sub-family binning → interpolate + 1-D Kepler refine →
turn-angle-min-max flyby assembly → AR/TR-gated `p.h.s.i` enumeration). This is the
faithful Russell method the failed DSM-Lambert lane lacked.

### Two substantive interpretation findings (sourced)
1. **Russell's "a (AU)" column in Tables 2.2/2.3 is the transfer PERIOD** (= a_sma^1.5
   in body periods), not the semi-major axis — the ^1.5 relation holds across 7
   independent anchor rows with ψ and |v∞| matching exactly. `GenericReturn.a_au`
   carries this column value; true SMA = `a_au^(2/3)`.
2. **Generic-return slow/fast branch** is NOT the Lambert energy branch; it is the
   phase of the final partial revolution (`frac = (tof/P_transfer) mod 1`; `frac>0.5`
   = Russell slow (N>0), else fast (N<0)). This partitions all anchor rows correctly.

## Catalogue-row batch (Russell's IDEALIZED model — circular-coplanar, Mars 1.875 yr)

12 descriptor-bearing rows; all 12 assemble a generic-return cycler at the row's
sourced Earth |v∞|. Note: the |v∞|-match is **by construction** (the assembler queries
at the sourced |v∞|), so it is NOT an independent validation. The emergent quantities
are AR (does the orbit reach Mars) and TR (are the re-initiating flybys ballistic):

| row | AR | TR | ballistic + reaches Mars? |
|---|---|---|---|
| **mcconaghy-2006-em-k2** (V0) | **1.238** | **1.694** | **yes** |
| russell-ch4-3.64gGg3 | 0.669 | 1.263 | no (AR<1) |
| russell-ch4-3.78Gg3 | 0.672 | 1.231 | no (AR<1) |
| russell-ch4-4.991gG2 | 0.768 | 0.530 | no |
| russell-ch4-8.049gGf2 | 0.660 | 0.702 | no |
| russell-ch4-9.353Gg2 | 0.673 | 0.872 | no |
| russell-ch4-5.30gGf3 | 0.780 | 0.501 | no |
| russell-ch4-9.94Gg3 | 0.682 | 0.664 | no |
| (4 of 12 have TR>1) | | | |

## Honest caveats — why this is NOT a promotion

1. **Idealized model only.** This is Russell's circular-coplanar patched-conic with
   Mars at 1.875 yr — NOT real DE440. A V0→V1 catalogue promotion requires real-eph
   closure (the existing dsm/correct lanes), which the spec explicitly scoped OUT.
2. **The assembly is approximate, not the catalogue cycler's exact geometry.** Multiple
   generic returns exist at a given Earth |v∞|; the descriptor→p.h.s.i map (no published
   crosswalk exists) selects one heuristically by (rev count, branch). The mismatch is
   visible: `russell-ch4-4.991gG2` assembles at AR 0.768 (aphelion ~1.17 AU) while the
   catalogue tabulates aphelion 1.64 AU (AR ~1.08) — the assembler picked a different
   generic return at that |v∞| than the catalogued cycler. So per-row AR/TR here are
   indicative, not the catalogued orbit's properties.
3. `mcconaghy-2006-em-k2` assembling as AR 1.238 / TR 1.694 (ballistic, reaches Mars in
   the idealized model) is **encouraging but not sufficient** — it is the idealized
   parent, and (per `project_s1l1_nomenclature`) the catalogue sequence may itself be
   mis-modeled. It stays V0; no writeback.

## Status

The Russell generic-return generator + global cycler search is real, tested, committed,
and **golden-validated against Russell's Tables 2.2/2.3/3.4** — a faithful capability
the project previously lacked. Its application to promote the catalogue SnLm rows is
**idealized-model-only and approximate**; real-ephemeris closure remains the gating
step for any V0→V1 promotion (future #388 work, or feeding these idealized parents as
seeds into the real-eph continuation). Out-of-plane (3D v∞ sphere) is deferred to #414.

## References
- spec / plan: `docs/superpowers/specs|plans/2026-06-21-russell-psi-generic-return-generator*.md`
- runlog: `data/runs/russell-cycler-20260620T160656Z.jsonl`
- Memory: `project_dsm_closure_modeljump_blocker`, `feedback_golden_tests_sourced_only`,
  `project_s1l1_nomenclature`, `feedback_orbit_closure_discipline`.
