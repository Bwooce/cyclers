# #388 option 3 — §14-V1-for-f/h-legs: well-posedness + data verdict (NOT BUILT)

**2026-06-23.** Investigated whether a §14 V1 mechanic could be built to promote the
f/h-leg census cyclers (russell-ocampo-4.3.1-5, -2.5.1+0, …) whose descriptors were
ingested this session. **Verdict: not well-posed AND data-blocked — not built.** A hollow
V1 gate would violate golden discipline, so it is deliberately not implemented.

## Why §14 V1 does not extend to f/h legs

The spec §14 V1 gate is a per-leg *reconstruction-faithfulness* check: re-solve each leg
with two independent Lambert solvers (lamberthub izzo2015 + gooding1990) and require
agreement < V1_TOLERANCE_MPS (1e-3 m/s), plus a Kepler re-propagation. It earned the
*generic-leg* russell-ch4 rows their V1.

1. **Resonant (`f`) legs have no independent second solver.** A full-rev M:N leg is
   analytically determined: a = a_E(M/N)^(2/3) (McConaghy-2005 Eq 14), |v_out| =
   √(μ(2/r_E − 1/a)) (Eq 15). Our `search/resonant_construct.construct_resonant_cycler`
   computes |v_out| by exactly this vis-viva — i.e. the *same* formula. There is no
   izzo-vs-gooding analog: a "V1 agreement" check would compare the formula to itself
   (circular). Half-rev (`h`) legs are the singular Lambert case (θ=π) — Lambert is
   undefined, so izzo/gooding can't be applied either.

2. **The genuinely-independent cross-check is cross-SOURCE**, per the McConaghy-2005
   digest: the McConaghy-2005 Table 2 `g/f/h` labels and the Russell-Ocampo 2003
   Tables 5-8 per-encounter ΔV *vectors* are two independent published representations of
   the same family. Reproducing BOTH would be a defensible V1. But those ΔV vectors are
   NOT ingested (the rows' source_quotes carry only V∞ scalars + leg ToFs).

## Data-level block (decisive, independent of the semantic one)

Even the construction path is starved: `construct_resonant_cycler` requires per-segment
`(a,e)`, but both rows have `orbit_elements.a_au/e = null` and an explicit `data_gaps`
entry: *"per-segment (a,e) of every arc not tabulated by Russell 2004 — only the
aggregate aphelion ratio is published"* (#54-backfill). Russell 2004 publishes only the
aggregate aphelion ratio + per-encounter V∞ scalars + cumulative-epoch ToFs — never the
per-arc conic elements or per-encounter ΔV vectors.

## Conclusion — the same ceiling, reconfirmed from a new angle

Both this session's parallel tracks land on the same bedrock:
- **Discovery probe (option 1):** Uranus/Neptune leveraging-VILM endgame is a structural
  empty (contours disjoint) — registry-recorded.
- **V1-build (option 3):** the f/h-leg census rows are **publication-gap-limited, not
  capability-limited.** No buildable mechanic promotes them, because the per-arc geometry
  (per-segment a,e / per-encounter ΔV vectors) is simply not in the published record.

This is the [[project_validation_ceiling]] finding reconfirmed: past the descriptor-bearing
generic-leg rows, promotion is **new-input-gated** (a source that publishes per-arc
geometry — e.g. a future Russell-school detailed-table release), not iteration-gated. A
hollow V1 gate is NOT built; the rows stay V0 with their descriptors ingested (data
hygiene) and this verdict on record.
