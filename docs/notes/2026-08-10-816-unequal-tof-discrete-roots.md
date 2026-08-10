# #816 — unequal-leg-time discrete asymmetric-closure enumeration at Uranus (2026-08-10)

**Verdict: CLEAN NEGATIVE for the question as posed — the unequal-tof formulation IS discrete
(1645 isolated roots enumerated, isolation restored exactly as the closed form predicts), 436 of
them are genuine dual-flyby asymmetric candidates, and ALL 436 fail the physical
required-turn-vs-achievable-bend wall, 0 survivors.** One (and only one) object in the entire
census passes every physical gate, but one of its two nodes requires *exactly zero* turn — a
one-node-working, passive-crossing architecture outside this genome's two-sided bend-usefulness
semantics; it is reported and flagged for adjudication (NOT claimed as a discovery), see §5.

Reproducible artifact: `scripts/scan_816_unequal_tof_discrete_roots.py` (chunked
`--directions`/`--merge` CLI; ~18 min pool-parallel / ~1.5 h serial; ruff clean). Full root list:
`data/found/816_unequal_tof_asymmetric_roots/roots.json` (+ per-chunk partials alongside).

## 1. Formulation (from the `#792` scoping closed form, one constraint relaxed)

`#792`'s scoping pass (`docs/notes/2026-08-10-792-scoping-vs-680.md` §2-§3,
`scripts/check_792_manifold_closed_form.py`) proved the equal-leg-time system is a 1-D continuum
that true periodicity collapses onto the symmetric `#563`/`#569` family. §4 identified the ONE
formulation in the lineage that is both genuinely asymmetric-capable and structurally
discrete-capable: independent leg times. This task ran it:

- Unknowns `(beta, tof0, tof1)`; equalities: E-match + h-match (== the two |v∞|-magnitude
  residuals, by the Tisserand-affine invertibility argument, which is leg-time-independent) +
  the pattern-repeat resonance `(n_a−n_b)·(tof0+tof1) ≡ 0 (mod 360°)`.
- The repeat condition is used to *eliminate* `tof1` exactly: `T = tof0 + tof1 = q·T_syn`,
  integer resonance order `q`. What remains is a square 2-D root-solve
  `F(beta, tof0) = [r_mid, r_periodic] = 0` per (pair, direction, q, n_rev0, n_rev1) — a
  finite, closed-form-informed enumeration (deflated Newton from grid local minima,
  `search/deflated_newton.py`), NOT a new adaptive grid search.
- Residual definition, `_leg_options` Lambert-branch convention, `#324` bend gate, and DOP853
  cross-check are reused from `scan_558_uranus_all_pairs_offset_sweep.py` verbatim; the only
  generalization anywhere is `tof0 ≠ tof1`.

Box: circular-coplanar patched-conic; moons Ariel/Umbriel/Titania/Oberon (the `#563` census
scope; Miranda excluded there on mass-deficiency), all 12 ordered (anchor, flyby) directions;
per-leg n_rev ∈ [0,3] (`#558` spec box, lineage min-departure-|v∞| branch selection); per-leg
durations ∈ [0.15, 3.6]·√(P_a·P_b) (superset of `#680`'s extended box); all q with T inside the
box (3-13 per direction, 82 (direction, q) solves total).

## 2. Positive controls (both pass; run before the sweep, asserted in-code)

- **PC1 (faithfulness)**: with `tof0 == tof1` the unequal-tof residual equals
  `scan_558.residual_at_point`'s residual at the catalogued `#312` point to 0.0 (bitwise).
- **PC2 (root recovery + isolation)**: the pipeline recovers the catalogued Ariel-Umbriel
  (0,0) symmetric golden (`data/enumerate_563_symmetric_closures.jsonl`,
  tof=3.216088179066208 d) as a root at beta≡0, tof0=tof1=T/2, residual 3e-11, with
  **cond(J) = 3.3e2** — versus 1e8-1e12 for the same point under `#680`'s free-(beta,tof)
  formulation. Fixing `T = q·T_syn` transversally cuts the along-manifold null direction and
  restores isolation, exactly as the closed form predicts.

## 3. Census result (all 12 directions, all q, all n_rev combos)

1645 validated isolated roots (residual ≤ 1e-9 km/s re-evaluated, max cond(J) = 9.9e4 —
everything isolated; max Tisserand-congruence defect |ΔE|,|Δh| ≤ 5e-5, i.e. every root's return
leg is a congruent copy of the outbound conic, confirming the structure prediction):

| class | count | meaning |
|---|---|---|
| symmetric_equal_tof | 424 | tof0=tof1, beta≡0/180 — the already-catalogued `#563`/`#569` family, recovered on the diagonal (expected) |
| asymmetric_candidate | 436 | tof0≠tof1, BOTH nodes require genuine turn — the objects this task exists to find |
| flyby_passive / anchor_passive | 395 / 390 | one node requires ~zero turn — one-node-working trivial-branch generalizations (see §5) |
| trivial_ballistic_resonant | 0 | fully-ballistic resonant Keplerian orbits (none converged inside the box) |

**All 436 genuine dual-flyby asymmetric candidates fail the physical gates — every single one
on required-turn infeasibility** (required turns 26°-179° per node vs achievable bends 4.6°-88°
at these moons' minimum safe periapses; most also by large margins). This is the same physical
wall the `#792` scoping measured along the equal-tof manifold (137°-156° required), now
confirmed at the discrete unequal-tof roots. The `#324` achievable-bend floor itself passes for
many candidates (the moons can bend 5°+); what kills them is that the *required* turn at an
asymmetric closure is an order of magnitude beyond what Ariel/Umbriel/Titania/Oberon can
deliver. 0 survivors. This is the registration's own predicted, acceptable outcome.

The `#565`/`#680` "necessary-not-sufficient" gap (V∞-magnitude matching without a required-turn
check) is CLOSED in this formulation: required turns at both nodes are computed per root
(moon-local frame, cycle-rotation-corrected at the anchor) and gated against the achievable
bend at that node's own V∞.

## 4. Why this settles the lineage's asymmetric question

- Equal-tof: closed-form continuum, collapses onto the symmetric family under true periodicity
  (`#680` + `#792` scoping — no isolated asymmetric closures EXIST).
- Unequal-tof (this task): discrete and well-posed, isolation confirmed numerically — and every
  genuinely-asymmetric isolated root is physically unreachable at these moons (bend wall).
- Together: within the circular-coplanar patched-conic anchor-flyby-anchor genome at Uranus, no
  physically-achievable genuinely-asymmetric true-periodic closure exists in any leg-time
  formulation. The asymmetric-closure question in this lineage is now closed at both the
  existence level (equal-tof) and the achievability level (unequal-tof).

## 5. The one all-physical-gates object: passive-node, flagged for adjudication (→ `#817`)

One physical trajectory (found twice, once from each direction convention — same object):

```
Ariel-Oberon  q=13 n_rev=(2,2) beta=13.575407° tof0=18.783412 d tof1=21.526054 d  (Oberon passive)
Oberon-Ariel  q=13 n_rev=(2,2) beta= 7.211811° tof0=21.526054 d tof1=18.783412 d  (same, mirrored)
residual 7.5e-11 km/s, cond(J)=145, DOP853 max |dr| = 4.9e-5 km,
V∞ = (1.311, 1.327, 1.311) km/s, required turns: Ariel 4.2327790° ≤ 8.04° achievable,
Oberon 0.000000000° (exactly zero to 9 decimals)
```

Structure (matches the closed form's complementary-arc branch exactly): both legs are
complementary arcs of ONE conic (a ≈ 414,500 km; arc durations sum to one orbital period;
T = 40.309 d = 5 full spacecraft periods = 13·T_syn(Ariel-Oberon)). The Ariel flyby performs a
4.23° radial-velocity-flip that maps the orbit onto its mirror image; the Oberon encounter
needs NO turn — the trajectory passes through Oberon's position ballistically. Its propulsive
skeleton is therefore a **single-moon (Ariel-only) repeated-flyby cycler whose coast arc is
phased to also encounter Oberon every cycle passively** — the Russell & Strange (2009)
passive-science-target cycler architecture (their Titan-Enceladus census), at Uranus. That
architecture is explicitly OUTSIDE this genome's two-sided bend-usefulness semantics (see the
`#571` stamp's own caveat recording exactly this structural distinction), so it is NOT a
survivor of this census and NOT a claimed discovery.

- `literature_check.py` (live WebSearch injected, 5-query trail): **not-found**
  (necessary-not-sufficient; nearest frame is Russell & Strange's moon-cycler architecture at
  Jupiter/Saturn — no Uranus pair, no Ariel-Oberon cycler in the searchable record).
- Registered as `#817` for Opus+Fable adjudication: is a one-node-working
  Ariel-cycler-with-passive-Oberon-rendezvous worth admitting under the expanded catalogue
  scope (quasi_cycler / precursor semantics), or is it subsumed by the repeated-moon census +
  passive-crossing being a free phasing condition? NOT my call; parked with full numbers here
  and in `roots.json`.

## 6. Discipline

- No catalogue write; no novelty claim anywhere.
- Gates untouched: `#324` floor 5.0°, DOP853 < 1 km, residual gate as lineage. The only
  addition is the required-turn feasibility check — a STRICTER, pre-declared physical gate.
- The `anchor_passive` class was added mid-run when the flagged object surfaced (the initial
  classifier only had the flyby-side passive class); partials are retrofit-reclassified at
  merge by a pure function of stored per-root numbers — no physics re-run, and the flagged
  object is reported, not suppressed.
- Empty-region stamped: `uranus-unequal-tof-asymmetric-discrete-roots-2026-08-10` in
  `data/empty_regions.jsonl` (method-versioned, capability-tagged, conditional on the
  circular-coplanar patched-conic model and the box in §1).
