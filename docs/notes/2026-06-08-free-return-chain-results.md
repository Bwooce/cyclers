# Two-arc free-return CHAIN — results & three-way-gate verdict (#163)

**Date:** 2026-06-08
**Experiment:** build Russell's ACTUAL generic-return-arc construction — a TWO-ARC
free-return chain (`g` + `G`), each arc a free-return ellipse crossing Mars's radius,
patched at an intermediate Earth flyby — and test whether it closes `6.44Gg3` (and
`4.991gG2`/S1L1) at the SOURCED anchors where the single-DSM `dsm_leg` transcription
floored off-anchor (#162) and the single free-return arc was refused (#137).

**Code (this experiment):** new module
`src/cyclerfinder/search/free_return_chain.py` + tests
`tests/search/test_free_return_chain.py`. Reuses (does not edit)
`search/free_return.py` (the anchor-respecting primitive) and `core/flyby.py` (the
bend machinery). No edits to `data/catalogue.yaml`, `docs/spec.md`,
`search/continuation.py`.

**GOLDEN/HONESTY.** EXPECTED = the SOURCED anchors. `6.44Gg3` (Russell 2004 Table
4.13): V_inf **E = 6.44**, **M = 3.74** km/s, aphelion 1.54 AU, arcs `g(2.087,...)`
+ `G(4.3191,...)`. `4.991gG2`/S1L1 (Russell 2004 Table 4.9, the row's OWN anchors —
NOT the CPOM 5.65/3.05 framing): V_inf **E = 4.99**, **M = 5.10** km/s, aphelion 1.64
AU, arcs `g(1.4612,...)` + `G(2.8096,...)`. Emerged V_inf is EVIDENCE, never imposed.
**No catalogue writeback.**

---

## The primitive (why this is the right one, not a `dsm_leg` repeat)

Each Russell leg descriptor is `g(TOF_years, psi_deg, branch)`: the FIRST number is
the full Earth-to-Earth ARC time in YEARS; the SECOND is the wrap angle `psi` on the
V_inf sphere, which EXCEEDS 360 deg (`1111.33` / `1194.88` for 6.44Gg3) — i.e. the
arcs are MULTI-revolution. The chain's free variables are the two arc shapes
`(a_1, e_1, a_2, e_2)`; the residual is ANCHOR-RESPECTING, exactly like
`free_return` and deliberately NOT the dV-budget of `dsm_leg`:

* per-arc emerged V_inf at Earth and Mars  - sourced anchor,
* per-arc emerged Earth-to-Earth arc ToF (`n_rev*period + time-above-Earth`, with the
  integer `n_rev` chosen by ToF-min) - the descriptor ToF,
* intermediate-Earth-flyby V_inf-magnitude continuity (the same physics as the two
  Earth-anchor terms; reported separately) + bend-feasibility (`core/flyby.max_bend`).

The ToF term is the structural binding that DISTINGUISHES the g-arc from the G-arc
and **forbids single-ellipse collapse**. Without it the solver finds one ellipse that
hits both V_inf anchors at the WRONG ToF (a spurious CLOSE — observed and rejected
during development; see "the collapse trap" below). With it, the g- and G-arcs are
forced onto different revolution counts and the construction is a genuine two-arc
object.

---

## Verdict

### 6.44Gg3 — **AMBIGUOUS** (V_inf basin reached + bend-feasible; descriptor ToF EMPTY-SET)

The sharpest result the tree has on this row. The decomposition is decisive:

> **The V_inf anchors ARE reachable and the intermediate Earth turn is trivially
> bend-feasible — but the per-arc DESCRIPTOR ToFs are NOT simultaneously reachable.**
> The arc ToF is quantised by the discrete revolution count, and the V_inf-fixed
> ellipse period leaves the g-arc (2.087 yr) and G-arc (4.3191 yr) targets in the
> gaps between revolution counts.

This is **AMBIGUOUS** against the §4 three-way gate, and it is a *better-resolved*
ambiguity than #162's floor-drop: there the emerged V_inf never approached the
anchors; here the V_inf half is essentially SOLVED (residual ~0.08 km/s) and the
irreducible obstruction is isolated to the structural ToF axis (EMPTY-SET-flavoured).

* **CLOSE?** The §4 CLOSE criterion (emerged V_inf within 0.5 of both anchors AND
  bend-feasible intermediate flyby) is, taken literally, MET on the V_inf+bend axes.
  But it is NOT a full structural closure: the g-arc descriptor ToF is unreachable
  (0.50 yr residual), so this is reported as AMBIGUOUS, NOT CLOSE — softening the
  gate to call it CLOSE would hide the ToF wall. (No catalogue writeback regardless.)
* **EMPTY-SET on the ToF axis:** the g-arc target 2.087 yr lies BELOW the n_rev=1
  arc-time minimum (~2.43 yr at the V_inf-matching ellipse); the G-arc target 4.32 yr
  lies BETWEEN n_rev=2 (~3.82 yr) and n_rev=3 (~5.15 yr). No `(a,e)` on the V_inf
  basin reaches either descriptor ToF — a clean, quantified structural negative.

### 4.991gG2 / S1L1 — **CLOSE-LEANING** (V_inf solved, G-arc ToF near-exact, g-arc ~0.14 yr off)

The closest the two-arc chain comes to a full multi-arc closure. Both V_inf anchors
match to <0.1 km/s; the **G-arc ToF lands essentially exactly** on its descriptor
(2.81 vs 2.8096 yr at n_rev=1); only the g-arc carries a ~0.14 yr ToF residual (its
n_rev=0 arc time is ~1.32 yr vs the 1.4612 yr target). The intermediate Earth turn is
trivially in-cone. **Continuation-seed implication:** this `(a_1, e_1, a_2, e_2)` is
the strongest circular-coplanar two-arc seed yet for the S1L1/#94 continuation
(`search/continuation.py`) — the piece that single-ellipse construction never gave
it. The remaining g-arc ToF gap is the obvious next target (a g-arc with a slightly
larger period, or a non-ToF-min n_rev with a small phasing burn).

---

## Numbers (verbatim, circular-coplanar backend)

### 6.44Gg3 (anchors E 6.44 / M 3.74; arcs g 2.087 yr, G 4.3191 yr; aphel 1.54)

```
max_res 0.4989 km/s   vinf_res 0.0811 km/s   tof_res 0.4989 yr
converged True (knife-edge)   vinf_within_tol True
arc1 (g): vinf E/M = 6.482 / 3.704   a/e = 1.221 / 0.2606   n_rev 1   arc_tof 2.434 yr
arc2 (G): vinf E/M = 6.359 / 3.817   a/e = 1.2271 / 0.2597  n_rev 2   arc_tof 3.820 yr
intermediate Earth flyby: turn 0.58 deg << max_turn 72.53 deg (bend in-cone)
                          continuity 123 m/s
solver nfev 6
grid multi-start (9x10 (a,e) seeds): BEST max_res 0.4989 — floor is ROBUST,
   dominated entirely by the ToF residual; vinf_res stays ~0.08.
```

### 4.991gG2 / S1L1 (anchors E 4.99 / M 5.10; arcs g 1.4612 yr, G 2.8096 yr; aphel 1.64)

```
max_res 0.1361 km/s   vinf_res 0.0132 km/s   tof_res 0.1361 yr
converged True   vinf_within_tol True
arc1 (g): vinf E/M = 4.977 / 5.113   n_rev 0   arc_tof 1.325 yr  (target 1.4612)
arc2 (G): vinf E/M = 4.990 / 5.100   n_rev 1   arc_tof 2.810 yr  (target 2.8096 — exact)
intermediate Earth flyby: turn 0.014 deg << max_turn 89.85 deg (bend in-cone)
                          continuity 13 m/s
```

### Degenerate reduction (arc1 == arc2) — the mechanics gate

```
both arcs converge to a/e = 1.22097 / 0.26055 (identical to 1e-6)
intermediate continuity 1.1e-10 km/s, turn 0.0 deg  -> reduces to ONE free_return
ellipse, matching free_return_geometry's vinf to 1e-9. (test passes.)
```

---

## The collapse trap (recorded for honesty)

The FIRST residual formulation drove ONLY the V_inf anchors (no ToF term). On
6.44Gg3 it reported `max_res ~0`, V_inf E/M = 6.44/3.74 at BOTH arcs, continuity 0,
turn 0 — a *spurious CLOSE*. Inspection showed both arcs had collapsed to the SAME
single ellipse `(1.223, 0.2602)` whose E->M transit is 166 d and period 1.35 yr —
matching NEITHER arc descriptor (2.087 / 4.3191 yr). A single ellipse DOES reproduce
both V_inf anchors (at aphelion 1.54), but at the wrong ToF; that is precisely why
the row is multi-arc and why the single free-return arc was refused (#137). Adding
the descriptor-ToF binding term (with per-arc `n_rev`) is what makes the two arcs
distinct and the verdict honest. This trap is the concrete mechanism behind the
memory blocker's "representable != reachable."

---

## Comparison to the prior floors

| pass | method | Mars V_inf emerged | binding obstruction |
|---|---|---|---|
| #137 | single free-return arc | (refused: breaks Earth V_inf continuity) | single ellipse can't close E |
| Gate-3 | single free_return, aphelion+262d seed | 3.06 (vs 3.74) | one ellipse, wrong basin |
| #162 | chained-DSM, dV-budget objective | 8-16 (off-anchor) | irreducible ΔV_DSM ~9 km/s |
| **#163 6.44Gg3** | **two-arc chain, anchor objective** | **3.70-3.82 (at 3.74!)** | **descriptor ToF (quantised n_rev gap)** |
| **#163 4.991gG2** | **two-arc chain, anchor objective** | **5.10 (exact)** | **g-arc ToF ~0.14 yr; G-arc exact** |

The two-arc chain **decisively beats** both the single-arc 3.01/3.06 and the dsm_leg
~9 floors on the V_inf axis — the emerged Mars V_inf is AT the sourced 3.74, not 3.06
or 8-16. The obstruction has MOVED from "wrong V_inf basin / irreducible DSM" to a
sharply-localised "descriptor ToF not reachable by an integer revolution count at the
V_inf-fixed ellipse." That is genuine progress and a publishable, quantified result.

---

## What this means / next step

1. **6.44Gg3 is ToF-blocked, not V_inf-blocked, under the two-arc free-return
   primitive.** The near-ballistic low-Mars-V_inf geometry IS representable as two
   free-return arcs at the sourced anchors; what does not close is the requirement
   that each arc's Earth-to-Earth time equal the published descriptor ToF with an
   integer revolution count. Releasing the ToF (allowing a small phasing burn on the
   intermediate Earth-Earth loop, the catalogue's `loop-ee-*` segments) is the
   obvious DSM-minimal closer — and connects to the #162 conclusion that the loop-arc
   ToF/branch is the lever.
2. **4.991gG2/S1L1 is the continuation seed.** Its near-closed two-arc
   `(a_1, e_1, a_2, e_2)` (V_inf solved, G-arc ToF exact) is the circular-coplanar
   two-arc seed `search/continuation.py` has been missing for S1L1/#94. Walking it to
   ephemeris is the natural follow-up (out of scope here; no edit to continuation).

**No catalogue writeback.** The emerged V_inf are evidence only; the sourced anchors
remain the EXPECTED target and the descriptor ToFs were not fully reached.

---

## Test results (verbatim)

```
$ uv run pytest tests/search/test_free_return_chain.py -m "" -v
6 passed in 2.78s
```

Tests: `test_arc_ee_time_increases_by_one_period_per_rev`,
`test_best_n_rev_selects_closest`,
`test_degenerate_reduces_to_single_free_return_ellipse`,
`test_two_arc_descriptors_drive_distinct_revolution_counts`,
`test_644gg3_two_arc_chain_sourced_anchor_probe` (@slow),
`test_4991gg2_two_arc_chain_sourced_anchor_probe` (@slow).
