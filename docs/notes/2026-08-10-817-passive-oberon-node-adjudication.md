# #817 — adjudication of `#816`'s passive-Oberon-node object (2026-08-10)

**Verdict: (b) NOT ADMISSIBLE — but on a THIRD ground, not either of the two the task
enumerated.** The object is neither subsumed by prior work nor excluded by any class
definition. It is disqualified because **it is not a physically realizable trajectory under
the model's own semantics**: the exact-zero required turn at Oberon and the closure's own
constraint that the spacecraft be *at* Oberon cannot both hold. There is nothing to admit.

Adjudicated at HEAD `49e4bb0b`. Adjudication only: no compute, no catalogue write, no
`empty_regions.jsonl` stamp, no novelty claim.

## The object (restated from `data/found/816_unequal_tof_asymmetric_roots/roots.json`)

```
Ariel-Oberon  q=13  n_rev=(2,2)  beta=13.575407°  tof0=18.783412 d  tof1=21.526054 d
T = 40.309467 d = 13·T_syn(Ariel,Oberon)      residual 7.5e-11 km/s   cond(J)=145
V_inf: Ariel 1.32705 / Oberon 1.31114 km/s    DOP853 cross-check 4.9e-5 km
required turns: Ariel 4.2327790° (achievable 8.0403°)   Oberon 0.000000000°
```

Derived here from the stored per-root invariants (`E = -6.988348 km²/s²`,
`h = 1305030.775 km²/s`, `mu_Uranus = 5794556.4 km³/s²` from `SATELLITES['Oberon']`):

- one conic, `a = 414,587 km`, `e = 0.53951`, `r_p = 190,915 km`, `r_a = 638,259 km`,
  `T_sc = 8.0645 d`. Periapsis sits essentially **on Ariel's orbit** (Ariel `sma = 190,929
  km`, 14 km out); apoapsis clears Oberon's orbit (`sma = 583,511 km`).
- `energy_leg0 == energy_leg1` and `h_leg0 == h_leg1` to ~1e-12 → **both legs are arcs of the
  same conic**, so the Oberon junction is velocity-continuous *by construction*. The 0° is
  structural, not a numerical coincidence.
- `T/T_Ariel = 15.9935`, `T/T_Oberon = 2.9935` (difference exactly `q = 13`). Neither is an
  integer: the object closes **up to a rotation** (≈ −2.34°/cycle), which is exactly what
  `scan_816_unequal_tof_discrete_roots.py:352-355` encodes by de-rotating `vin_a2` by
  `theta = n_a·t_total` before measuring the anchor turn.

## 1. Class fit — the schema WOULD accommodate it; the search genome never has

Answering sub-question 1 affirmatively rather than as a rejection:

- The five classes (`docs/notes/2026-06-16-catalogue-scope-taxonomy.md` §"The four classes",
  `docs/spec.md`, `data/README.md`) are defined on **periodicity / epoch-lock / n_returns**.
  Not one of them says anything about bend-usefulness at every node. Nothing in the schema
  excludes a one-working-node architecture.
- Closes-up-to-rotation + epoch-locked + finite returns ⇒ `quasi_cycler` is the only
  candidate class, and it is the *same* semantics as the object's own sibling row
  `ariel-oberon-1-1-uranian-quasi-cycler-2026` (same pair, same primary, same
  anchor-flyby-anchor genome, `orbit_class: quasi_cycler`, V4-windowed).
- The two-sided bend rule is a **search-genome gate** (`#324`
  `DEFAULT_MIN_USEFUL_BEND_DEG`), not a class definition. The `#571` empty-region stamps say
  so in their own words: the passive-target architecture is *"conditional on this project's
  own gate POLICY, not a universal physical claim"* and is *"a STRUCTURALLY DIFFERENT,
  uncovered case for this genome"*. **Uncovered ≠ out of scope.** The project has never
  searched this architecture; it has never ruled it inadmissible either.

So the correct reading is: a genuine Russell-Strange-style passive-target cycler at Uranus
*would* be admissible as `quasi_cycler`. This particular object is not one — see §2.

## 2. Russell & Strange 2009 is the wrong frame FOR THIS OBJECT (decisive)

R-S 2009 (`docs/notes/2026-06-30-digest-russell-strange-2009-planetary-moon-cyclers.md`) is
self-consistent because the **passive target is tiny and the V_inf is high**: Titan is the
flyby body ("the largest moon in the Saturn system by almost 2 orders of magnitude"),
Enceladus/Dione/Mimas/Rhea/Tethys are targets. Enceladus (`mu = 7.2`) at ~4 km/s deflects a
100-km-altitude pass by ~0.15° — genuinely negligible.

This object **inverts every term of that**:

| | working node | "passive" node |
|---|---|---|
| body | Ariel | Oberon |
| `mu` (km³/s²) | 83.5 | **205.3 (2.46×)** |
| V_inf (km/s) | 1.32705 | **1.31114 (lower)** |
| achievable bend at safe altitude | 8.0403° | **14.7426° (largest in the trajectory)** |

Oberon is the **most** massive body of the pair, met at the **lowest** V_inf, and therefore
carries the **greatest** bend authority of any node. It is the least defensible passive
target in the pair.

**Quantitatively.** Using the same patched-conic deflection the project's own gate uses,
`delta = 2·asin(1 / (1 + r_p·V_inf²/mu))` — which reproduces the stored
`max_bend_deg_per_encounter` values **bit-for-bit** (14.742638012954531° at Oberon,
8.040293699481882° at Ariel), so this is the repo's own physics, not an outside model:

| Oberon pass distance `r_p` | unmodelled deflection | as % of the 4.2328° Ariel turn budget |
|---|---|---|
| 811 km (safe altitude) | 14.7426° | 348% |
| 5,000 km | 2.6734° | 63% |
| **9,678 km (Laplace SOI)** | **1.3968°** | **33%** |
| **13,288 km (Hill radius)** | **1.0207°** | **24%** |
| 50,000 km | 0.2730° | 6.5% |
| 160,000 km | 0.0855° | 2.0% |

The genome places the spacecraft **at Oberon by construction**, not by one equation among
three. `_states` (script line 176-183) returns Oberon's state at `t = tof0` and both legs are
Lambert arcs `r0 → r1 → r2` through it; the two-component residual
(`residual_vec_unequal`, line 186-203) is `[opt0.vinf_in - opt1.vinf_out,
opt0.vinf_out - opt1.vinf_in]` — a pair of **|V_inf|-magnitude** matches, *not* a radius
residual (the name `r_mid` is `residual`_mid, not `radius`_mid; grounded against the
producing function per `[[feedback_verify_metric_semantics_before_ranking]]`). So position
coincidence at Oberon is structural to the genome and holds at every root, while the required
turn there is **exactly 0°** — which in patched-conic terms means periapsis at infinity.
These are contradictory. The fork is total:

- **Treat it as a real encounter** (anywhere inside Oberon's SOI, ≤ 9,678 km): Oberon injects
  ≥ 1.40° of deflection the closure never modelled — ≥ 33% of the entire turn budget Ariel
  supplies. The closure is invalid as computed.
- **Treat the 0° literally**: the spacecraft must stay ≳ 160,000 km from Oberon to keep the
  parasitic deflection under 2% of the working turn — ~12 Hill radii, ~17 SOI radii, ~27% of
  Oberon's own orbital radius. At that separation there is no encounter at all, and the
  object is a bare single-moon Ariel free-return with no second body involved.

This is exactly the class of self-consistency failure `[[feedback_constructed_tour_per_encounter_self_consistency]]`
pins from `#480`: *every encounter must be within the SOI of its node; an analytic match at
the node is not sufficient.* Here the node geometry and the node dynamics contradict each
other.

*(Scope: verified for this pair only. The generalization — that every passive-node root in
`#816`'s box fails the same way, all four moons being `mu ∈ [83, 235]` at 1-3 km/s — is
plausible but **not checked here**; it is the hypothesis motivating `#818`.)*

## 3. Then what IS it? — both branches lead to "no row" (sub-question 2)

**Not subsumed, strictly.** `repeated-moon-uranus-sweep` is the `#254` **two-working-node**
A-B-A genome with no bend gate at all (residual gate 0.05 km/s, 9019 closed / 0 routed); it
would never have flagged a passive-node object. No single-moon free-return Uranian sweep
exists anywhere in `data/empty_regions.jsonl` (the closest, `uranus-ariel-umbriel-ariel-`/
`uranus-titania-oberon-titania-multirev-leveraging-2026-06-26`, are different pairs and a
different endgame genome). So "already known via a prior stamp" is **not** the reason.

The reason is the fork:

- **Strip the Oberon coincidence** → a bare single-moon Ariel resonant free-return
  (periapsis on Ariel's orbit, 5 spacecraft revs, apse-line-mirror flyby). `#816`'s own
  module docstring **pre-declared** this before the run: *"the zero-apsidal-rotation branch
  is a trivial resonant Keplerian orbit … classified out as `trivial_ballistic_resonant`,
  not a discovery"*, and `classify_root`'s docstring names `anchor_passive`/`flyby_passive`
  as *"the two one-node-working **generalizations of the registration's pre-declared trivial
  branch**"*. That the exclusion was declared **before** the run, not rationalized after, is
  what makes it binding.
- **Keep the Oberon coincidence** → physically impossible, per §2.

Either way: no catalogue row.

Supporting, not load-bearing: the passive branch is not rare in this box —
`flyby_passive` (395) + `anchor_passive` (390) = **785 of 1645 roots (48%)**. The object is
distinguished only by being the one member (found twice, once per direction convention)
whose single working node's required turn is achievable. And per `#577`'s precedent, a
known-class member found by our own enumeration at a new body pair is not novelty-bearing.

## 4. Honest confidence (sub-question 3)

**High confidence this is not a meaningful physical discovery; high confidence it is a
formulation artifact.** But *not* for the reason the task's option-(c) framing offers. The
near-zero turn is **not** "a coincidental crossing any dense enumeration would eventually
produce" — it is exact (0.000000000°), structural (one conic, `E` and `h` equal to 1e-12),
and populates a whole 785-root branch. It is a *degeneracy of the formulation*, which is a
sharper and more damning diagnosis than coincidence.

`literature_check.py`'s NOT-FOUND remains necessary-not-sufficient and is **not relied on
here** in either direction: the verdict rests on physical self-consistency, which is
upstream of novelty. Nothing here says the record contains this object; it says there is no
realizable object to look for.

Confidence is lowest on one point, stated plainly: if a future genome carries a **b-plane
degree of freedom** at the science-target node, a properly-modelled Ariel-cycler-with-Oberon-
science-pass at Uranus is a legitimate, currently-unsearched object. That is registered as
`#819` and is **not** implied or recommended by this verdict.

## 5. Relationship to `#816`

This **strengthens** `#816`'s call, it does not correct it. `#816` declined to claim the
object on *policy* grounds ("outside this genome's two-sided bend-usefulness semantics") and
parked it with full numbers — exactly right. This adjudication upgrades that to a *physical*
disqualification with reproducible arithmetic.

`uranus-unequal-tof-asymmetric-discrete-roots-2026-08-10` needs **no amendment**: the object
was never stamped as a survivor there (`n_survivors: 0`, reported separately as
`n_passive_node_physical_gate_passers: 2`). Any addendum recording this verdict is the
coordinating session's call; this task touched neither `catalogue.yaml` nor
`empty_regions.jsonl`.

## 6. Follow-ons registered (not dispatched)

- **`#818`** — a reusable passive-node self-consistency gate: for any root with a
  near-zero required turn at a node, compute the parasitic deflection at that body's own
  SOI/Hill boundary and reject unless it is below a declared fraction of the working turn
  budget. Generalizes the `#480` SOI rule to the *dynamical* side. Would have auto-classified
  this object out with no adjudication needed. Small: one function + test + a `#816` re-run
  of the 785 stored passive roots (pure post-processing of stored numbers).
- **`#819`** — speculative, explicitly NOT implied by this verdict: a genuine R-S-style
  Uranian cycler with a *properly modelled* second-moon science pass (finite altitude,
  deflection carried as a working constraint, extra b-plane DOF). A campaign, not a check.
