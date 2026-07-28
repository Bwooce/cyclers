# `#755`: Targeted search for Anderson-Lo 2011 Table-1 rows 3:4-LO and 5:6-LO

**Task:** `#755`, a direct continuation of `#753` (see
`docs/notes/2026-07-28-753-jupiter-europa-resonant-families-table1-gate.md` for
the full prior search history, sourced constants, and the module's
methodological findings — this note does not repeat that background). `#753`
confirmed 1 of Anderson & Lo 2011's Table 1's four rows (5:6-LI); this task's
job was the two remaining in scope: 3:4-LO (target `|λ|=1036.116088`) and
5:6-LO (target `|λ|=4445.387515`). 5:6-NO was explicitly out of scope.

**Code:** `src/cyclerfinder/search/jovian_resonant_families.py` (extended, not
rewritten) + `tests/search/test_jovian_resonant_families.py` (25 tests, all
passing, up from 20). Read the source paper directly this task (PDF pages
13-14 and 17-18, i.e. printed pp.179-180 and 183-184) to see Figs. 9-11 and 16
as images, not just OCR text — this was essential (see below).

---

## Result summary

| Row | Target `\|λ\|` | Recovered | Rel. err (eigenvalue) | Period/2π | q | Period rel. err | Gate |
|---|---|---|---|---|---|---|---|
| 3:4-LO | 1036.116088 | **1036.116117** | **2.8e-8** | 4.08591 | 4 | **2.1%** | **FAIL** (period) |
| 5:6-LO | 4445.387515 | 982457 (closest found near-integer period) | 220x too big | 6.25006 | 6 | 4.2% | FAIL (both) |

Both rows are honestly reported as **NOT CONFIRMED** under the task's own
mandated dual criterion (eigenvalue match AND period on a clean `2πq`
multiple). But the two negatives are qualitatively very different, and that
difference is the main finding of this task.

---

## 3:4-LO: a striking, well-evidenced near-miss — not a confirmation

Extending `#753`'s own strategy 1 recommendation ("a finer C-grid... a
different simple-loop starting family... needs real numerical exploration"),
a much finer `(x0, half_crossings)` grid at `C_flyby` directly, in the same
`x0 ≈ -1.42` to `-1.45` region `#753` had already flagged as a "fractal
sensitivity" hotspot (huge, wildly x0-sensitive eigenvalues packed within
`Δx0 < 0.01`), located:

```
x0 = -1.4304078294961569, ydot0_sign = +1, half_crossings ∈ {5, 6, 7}
period = 25.672528919046933  (period/2π = 4.085910)
crossing_residual = 2.07e-13
max_eigenvalue (Barden, real) = 1036.116116695996
planar_floquet cross-check    = 1036.116095551334   (agree to <1e-7 relative)
target (paper, Table 1)       = 1036.116088
relative error                = 2.77e-8
```

This is essentially an **exact** reproduction of the paper's own published
digits — several orders of magnitude tighter than the `1e-3` gate, and using
BOTH independent eigenvalue extractions in the module (Barden vs
`_planar_floquet`) agreeing to <1e-7 relative. It was found via a genuinely
different, much finer region of `(x0, half_crossings)` space than `#753`
sampled (`n_grid` up to 2500-4000 in `Δx0 ≈ 0.01`-`0.09`-wide windows, vs
`#753`'s coarser sweeps), confirming `#753`'s own hypothesis that its
"fractal sensitivity" hotspot region harbored the real family — it just
needed finer resolution to hit the exact point.

**Independent corroboration beyond the eigenvalue** (not part of the formal
digit gate, but strong supporting evidence this is genuinely the paper's own
3:4-LO orbit, not a coincidental fractal neighbor):

- **Spatial envelope matches Fig. 16(a) closely.** The trajectory spans
  `x ∈ [-1.430, 1.258]`, `y ∈ [-1.384, 1.384]` — the paper's own plotted
  "flower" orbit at `C_flyby` (Fig. 16a, axes `-1.5` to `1.5`) is visually
  almost exactly this size and shape (a near-circular loop around Jupiter
  with several small cusps/loops, one right next to Europa).
- **Genuine close Europa approach.** Minimum distance to Europa along the
  orbit is `0.00247` nondim (`≈ 1657 km`, `≈ 97 km` surface altitude) — matching
  the paper's own qualitative mechanism for why these orbits are so unstable
  ("close flyby of Europa" causing large local Lyapunov exponents, p.177-178).

**What does NOT confirm:** `period/2π = 4.08591`, a real, tightly-converged
(`crossing_residual = 2e-13`, not tolerance noise) `~2.1%` offset from the
naive `q=4` value. This fails even a generous `1%` period-tolerance (this
module's own `TABLE1_PERIOD_REL_TOL = 1e-2`).

**Why this offset might be expected (a hypothesis, not a resolved fact):**
`T_full = 2πq` is the closure condition for the paper's own **two-body**
seed construction (Eq. 5-6) — it is the period a `p:q` two-body resonant
ellipse *must* have to close in the rotating frame. But a genuinely periodic
**three-body** orbit only needs `X(T) = X(0)` for *some* `T`; nothing in the
CR3BP forces that `T` to equal `2πq` exactly once the orbit is no longer
close to the two-body integrable limit. The paper's own text says as much
directly (p.171, right after introducing the `p:q` notation): *"In the
two-body problem, resonances are precisely related by integers, but here
[i.e. for the actual CRTBP orbits/Poincaré-map islands] ... **a precise
relationship does not exist**, and islands in Poincaré sections representing
quasi-periodic orbits are used to determine resonance."* `5:6-LI` — the ONE
family `#753` confirmed — is explicitly described by the paper as "only
slightly unstable" (`λ ≈ 1.000008`, barely off the marginal-stability
boundary), i.e. very close to the two-body/quasi-periodic limit, which is
exactly the regime where `T ≈ 2πq` would be expected to hold tightly (and it
does, to `2.5e-6` relative). `3:4-LO`, by contrast, is **extremely** unstable
(`λ ≈ 1036`, three orders of magnitude further from neutral stability) — a
family that has moved far from the two-body limit by the time it reaches
`C_flyby`, for which there is no a priori reason its period should still
track the naive two-body value closely.

If this hypothesis is right, the "period lands on a clean `2πq` multiple"
criterion (validated only on the one family, `5:6-LI`, where it happens to
hold almost exactly) may not be the correct universal confirmatory test for
the *strongly* unstable rows. This task does **not** unilaterally decide
that question — `TABLE1_PERIOD_REL_TOL`/`GateRow.passed` in the module keep
the strict dual criterion exactly as specified, so `3:4-LO` is reported as
**NOT CONFIRMED**, full stop. But the eigenvalue match is so precise, and the
qualitative distinction from `#753`'s original near-misses so large (those
were `2%`-`27%` off on eigenvalue itself, with periods bearing NO relation to
any plausible `q` — `21.01`, `16.11`, `16.04`), that this row deserves
explicit flagging as a much stronger, better-evidenced candidate that a human
reviewer should weigh, not lumped in with an ordinary miss.

---

## 5:6-LO: unchanged, genuinely NOT confirmed

Applied the identical strategy (wide `x0 ∈ [-1.6,-0.85]` scans with
`half_crossings` from 1 to 25, both `ydot0` signs, then repeated fine-grid
refinement — up to `n_grid=4000` in `Δx0 ≈ 0.02`-wide windows — around the
same `x0 ≈ -1.42` "fractal hotspot" that yielded `3:4-LO`'s match) to 5:6-LO.
Unlike 3:4-LO, this did **not** turn up anything close:

- Many candidates in the hotspot region have huge eigenvalues (`10²`-`10⁶`,
  wildly x0-sensitive, confirming the same fractal structure `#753` already
  documented) — but **none** land within an order of magnitude of the target
  `4445.387515` while simultaneously having `period/2π` anywhere near `6`.
- The single closest simultaneous approach found: `x0 = -1.424533`,
  `λ = 982457` (`223×` too big), `period/2π = 6.25006` (`4.2%` off `q=6`).
- Candidates with `period/2π` landing exactly on integers (e.g. `6.00017` at
  `x0=-0.869280`) all have small, weakly-unstable eigenvalues (`~1.9`), the
  same "trivial near-2-body-limit" pattern `#753` already found and rejected
  for this row.

This is a genuinely different outcome from 3:4-LO — not a symmetric partner
finding. `5:6-LO`'s hardcoded seed is left unchanged from `#753` (the best
candidate remains the pre-existing one, still not confirmed).

---

## Strategy 2 (implemented, did not itself find a match): flyby-vector-rotation seed

Implemented Anderson & Lo's own two-body flyby-VECTOR-ROTATION construction
(Section "Designing Flybys Using the Two-Body Approximations", pp.172-174,
Fig. 2 — read directly from the PDF text this task, distinct from and more
sophisticated than the plain resonant-ellipse seed `#753`'s module already
had): `two_body_flyby_rotation_seed` models a spacecraft on a `p:q` two-body
resonant ellipse encountering the secondary and having its hyperbolic excess
velocity (`V∞`) ROTATED by the turn angle a hyperbolic flyby at a chosen
periapsis radius implies (standard patched-conic geometry,
`δ = 2·arcsin(1/(1 + r_p·V∞²/μ))`), producing a new two-body orbit in a
different resonance — literally how the paper's own text describes
constructing its cycling 3:4↔5:6 flyby trajectory. `flyby_rotation_symmetric_seed`
propagates the (generally non-symmetric) post-flyby state forward to its
next perpendicular x-axis crossing, producing an `(x0, jacobi)` seed usable by
`converge_candidate`/`survey_candidates`.

Verified geometrically correct (`V∞` magnitude is preserved under the
rotation to `<1e-9` relative; turn angle correctly decreases with periapsis
radius). Hit the SAME numerical hazard the paper's own text documents for its
"crudest method" (p.174: close approaches to the singularity at Europa kept
the differential corrector from converging) — the exact periapsis point sits
only `μ ≈ 2.5e-5` from the CRTBP's actual singularity at `1-μ`, causing
DOP853 to grind for tens of seconds per near-singular passage; fixed with the
same tactic the paper itself used ("slightly modifying the patchpoints just
prior to Europa approach", p.175) via a `safety_margin` back-off parameter.

Did not itself locate a Table-1 match in the time available this task
(genuinely reproducing the paper's own full multi-patchpoint iterative
refinement procedure, p.175-176, is a substantially larger undertaking than
fit in this task's budget). Kept as a documented, tested, reusable seed
strategy — a real second tool for a future task, not dead code.

---

## Verification

- `uv run ruff check` / `ruff format --check` on both changed files: clean.
- `uv run mypy src tests` (project canonical invocation, per
  `feedback_verify_must_include_full_mypy`): clean, 821 files.
- `uv run pytest tests/search/test_jovian_resonant_families.py -q`: 25/25
  pass (~75s wall time).
- `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: pass (run
  before committing the `OUTSTANDING.md` update).

---

## Reviewer verdict (coordinating session, 2026-07-28)

Independently re-read PDF p.171 directly (not just this task's quote) to check the Eq. 6
citation in full context before ruling: confirmed accurate and unambiguous — *"In the
two-body problem, resonances are precisely related by integers, but here \[i.e. the actual
CRTBP\] ... a precise relationship does not exist, and islands in Poincaré sections
representing quasi-periodic orbits are used to determine resonance."* The paper itself,
not a hypothesis invented after the fact, defines CR3BP resonance as approximate
(`p·n_p ≈ q·n_q`).

**Ruling: 3:4-LO is CONFIRMED**, on the weight of evidence — not by loosening
`TABLE1_PERIOD_REL_TOL` retroactively to force a pass (that would be exactly the kind of
tolerance-fudging this task was explicitly told not to do; the module's own strict
dual-criterion gate correctly and honestly still reports `FAIL`, and stays as-is). The
basis for the ruling:
1. Eigenvalue match to `2.8e-8` relative — tighter than the already-confirmed `5:6-LI`'s
   own `9.6e-5` by four orders of magnitude, at the corrector's own convergence floor.
2. Independent corroboration: trajectory shape matches the paper's own Fig. 16(a), and the
   close-Europa-approach mechanism/altitude (~97 km) matches the paper's own qualitative
   description of this family's instability source.
3. The period offset (2.1%) is not an isolated excuse — it fits a coherent, physically
   principled pattern rather than special-pleading for this one case: `5:6-LI` (weakly
   unstable, `λ≈1.000008`, close to the two-body limit) tracks `2πq` almost exactly
   (`2.5e-6` relative), while `3:4-LO` (extremely unstable, `λ≈1036`, three orders of
   magnitude further from marginal stability) is exactly the regime the paper's own Eq. 6
   says should NOT be expected to hold a precise integer period ratio. This is the opposite
   quality of evidence from `#753`'s original rejected candidates, which were wildly off on
   BOTH eigenvalue (2%-27%) AND period (`21.01`/`16.11`/`16.04` vs `4`/`6`/`6`) simultaneously
   — those were correctly rejected as a different family entirely; this one is not that case.

**Operational implication for future `5:6-LO` search**: if strongly unstable families
genuinely drift from the naive `2πq` period (as this ruling concludes for `3:4-LO`), then
`#755`'s own `5:6-LO` search — which explicitly favored candidates with `period/2π` near an
integer — may have been filtering out the very candidate it was looking for. A follow-up
`5:6-LO` search should relax the period-proximity search criterion (not just the final gate)
alongside the eigenvalue-magnitude search, consistent with this ruling.

**Net effect on `#754`**: 2 of the 2 families Task B needs are no longer both unconfirmed —
`3:4-LO` now stands confirmed, `5:6-LO` remains a genuine, unchanged gap. `#754` still needs
`5:6-LO` (or a re-scope around only `3:4-LO`+`5:6-LI`) before its own Table 2/3 gate is
fully in reach.

---

## Opinion on `#754` (Task B) — not a decision

`#753`'s own recommendation was to hold Task B until 3:4-LO/5:6-LO confirm or
a re-scope happens. This task's result is a genuine, well-evidenced PARTIAL
step forward on exactly one of those two rows (3:4-LO's eigenvalue, to
extreme precision) but an unresolved question about whether the
period-exactness criterion is the right bar for that row, and a continued
clean negative on the other (5:6-LO). My opinion: this does not yet clear
`#753`'s own bar for dispatching Task B as originally scoped (which needs
BOTH 3:4-LO and 5:6-LO's manifolds for the Table 2/3 gate) — 5:6-LO remains
entirely unconfirmed, and 3:4-LO's confirmation status hinges on a genuine,
unresolved methodological question (is period-exactness even the right test
for a strongly unstable family?) that deserves a human decision, not a
unilateral call by this module. If a reviewer decides 3:4-LO's near-machine-
precision eigenvalue match, corroborated by its Fig-16(a)-matching shape and
close-Europa-approach mechanism, is sufficient despite the period offset,
that would put one of the two Task-B-blocking rows in hand — still short of
both. This is my assessment for the user to weigh, not a decision.
