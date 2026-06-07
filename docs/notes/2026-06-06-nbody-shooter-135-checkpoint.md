# N-body Jones shooter — C.0 #135-verdict checkpoint

Task #136 Phase C, Task C.0. Read the landed #135 like-for-like verdict
(`docs/notes/2026-06-06-russell12-likeforlike.md`) and decided the shooter seeding
strategy from it BEFORE writing any solver code, as the plan mandates.

## The #135 verdict (re-read at run time)

The like-for-like diagnostic ran coplanar-vs-coplanar on known-solvable
Russell/McConaghy instances and found **0 CLOSE-AND-MATCH**: the corrector closes
geometrically but lands **off-anchor** (our V∞ ≈ 9–28 km/s vs sourced 3–10 km/s)
even with the model mismatch removed. The verdict is **seeding/basin, not solver
deficiency**.

This is corroborated, for the specific Jones VEM rows, by the chain of refutations
recorded in `data/OUTSTANDING.md`:

- #110: dense scan (2816 pts/row) — floor ~17.9 / 18.5 km/s, 0 bend-feasible.
- #120: 3D inclination — moved the floor < 0.4 km/s; REFUTED.
- #122: vector residual (in-solver bend hinge) — collapsed the powered basin to
  ~0 closures; the one bend-feasible family sits at ~21 km/s; REFUTED.
- #137 free-return: the E-M **sub-arc** reaches the Jones Mars-V∞ floor (~2.81),
  but the multi-ellipse-through-Venus-flyby structure is the missing piece —
  exactly the flyby-propagation shooter, not a conic-genome addition.

## Decision (binding for C.1–C.5)

The verdict has NOT changed from the plan-authoring read. The shooter is therefore
seeded from the **#133 near-miss conic survey** (`near_miss_survey`): the
lowest-V∞ conic chains within a relaxed near-miss tolerance (0.5 km/s — Jones's
own stage-1 <=200 m/s analogue, AAS 17-577 method deep-dive §2.2), ranked by
max-V∞ (closest to the Jones basin first), NEVER from a blind equispaced/coplanar
scan. The `cyclerfinder.nbody.shooter` module docstring states this mandate and
cross-references #135; `tests/nbody/test_c0_135_checkpoint.py` guards it.

## #142 method deep-dive incorporated

The #142 mining note `docs/notes/2026-06-07-jones-aas17-577-method-deepdive.md`
landed during this run and was read. It confirms the Phase C blueprint
(broad-search near-Hohmann seeding -> B-plane powered-flyby feasibility -> SNOPT
n-body correction) and supplies the published tolerances now used as the shooter's
provenance: seed window <=50 d from Hohmann + seed v∞ <= 5 km/s; flyby altitude
window 100–100,000 km; interior-flyby Δv∞ tolerance 100–200 m/s; **SNOPT
continuity 1.0E-3 km / 1.0E-6 km/s** (adopted verbatim as the shooter's
ballistic-continuity acceptance floor, `shooter._POS_CONTINUITY_KM` /
`_VEL_CONTINUITY_KMS`). The flyby-feasibility kernel (their Eqs.1-3 powered-flyby
periapsis r_p) is implemented as `_flyby_periapsis_hinge_km`.
