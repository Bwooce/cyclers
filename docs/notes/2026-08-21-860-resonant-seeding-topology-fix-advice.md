# #860 — Resonant Atlas seeding-topology blocker: diagnosis + path-forward advice

**Date**: 2026-08-21
**Task number**: `#860`, verified next-free at write time (grep of `data/OUTSTANDING.md` and
`docs/notes/` finds `#859` as the highest allocated number; no `#860`+ anywhere).
**Scope**: advice only. No code written or run, no file modified except this note.
**Sources read in full**: `src/cyclerfinder/search/resonant_atlas_stage_a.py` (the `#859`
harness), `jovian_resonant_families.py` (module docstring incl. the `#753`/`#755`/`#756`/`#758`
update history; `two_body_resonant_seed`; `two_body_flyby_rotation_seed` +
`flyby_rotation_symmetric_seed`; `converge_candidate`/`_classify`), `cr3bp_continuation.py`
(full), `cr3bp_jacobi_arclength.py` and `mu_continuation.py` and `deflated_newton.py` and
`pseudo_arclength.py` and `deflated_variational_periodic_orbit.py` (docstrings/contracts),
`earth_moon_resonant_families.py` (header), the `#776`/`#777`/`#859`/`#849` bullets in
`data/OUTSTANDING.md`, `docs/notes/2026-08-21-859-resonant-atlas-pilot-harness.md` (full), the
`#728` Anderson & Kumar AAS 24-288 digest, and the flyby-seed tests in
`tests/search/test_jovian_resonant_families.py`.

## 1. The question

`#859`'s smoke test (independently re-verified by the coordinating session) found that the
`two_body_resonant_seed(p, q, x0_sign=-1)` → `correct_symmetric_fixed_jacobi` →
`continue_family` pipeline converges cleanly at every tested Uranus-Oberon published ratio
(4:5, 5:6, 4:3, 5:4, 6:5) but lands on a **stable, |λ| ≈ 1.0 near-unit-circle family every
time**, never on the published unstable saddle family the campaign exists to survey. `#776`
documented the same failure independently on Neptune-Triton. Since the whole point of
`#789`/`#859` is systems with **no published IC table to seed from**, this is a campaign-level
blocker: fix the seeding, redirect, or shelve?

## 2. Recommendation up front

**Fix it — do not shelve, and do not run Stage A as currently seeded.** The mechanism is now
well understood, and the two best fixes are (i) a ~15-line *conjugate-apse* ("encounter-phase")
seed variant that targets the missing branch directly from the numerically safe side, and
(ii) replacing `continue_family`'s natural-parameter walk with the already-built, already-proven
pseudo-arclength fold-turning walk (`cr3bp_jacobi_arclength.py`, `#249`) so the reliably-found
stable branch becomes an on-ramp to the *whole* family curve, saddle segments included. Both are
cheaply validatable against the Oberon positive control (published per-family Jacobi ranges,
AAS 24-288) in minutes-to-hours of compute before any novel-system CPU is spent; total effort
~2-4 days. If the Oberon gate still fails after both (plus the mu-continuation cross-check in
Sec. 4d), *then* shelve Stage A — a sweep whose misses are uninterpretable is worse than no
sweep.

## 3. Mechanism diagnosis

Three interacting causes, all visible in the code and the project's own prior task record. This
is **not** a bug in the corrector/continuation wiring (the `#859` note's fold-detection checks
confirm the wiring independently) and **not** simple bad luck — it is systematic.

**(a) The seed's phase choice structurally selects the encounter-avoiding branch.**
`two_body_resonant_seed` places the ellipse's r = 1 apse (periapsis for exterior ratios a > 1;
for interior ratios a < 1 the stored `e` goes negative and r = 1 is actually the *apoapsis* —
worth an explicit comment someday) on the x-axis. `x0_sign` picks which side: `+1` puts that
apse at the secondary's own longitude at t = 0 (the conjunction/close-encounter configuration),
`-1` puts it at opposition. The harness surveys `x0_sign = -1` **only**, for a real documented
reason (DOP853 step-collapse at the singularity for `+1` seeds). But the instability of the
target families **is** the repeated close flyby: Anderson & Lo attribute it explicitly
(pp. 177-178), and this project's own `#758` evidence is decisive — the genuine Jupiter-Europa
5:6-LO makes a 668 km Europa approach and 3:4-LO a 1641 km approach, while *every* wrong
candidate from `#756`'s wide sweep stays 12,000-28,000 km away with |λ| nowhere near target.
A resonant orbit whose geometry avoids the secondary is precisely the near-integrable,
elliptic-island-center kind: Barden returns a complex unit-modulus pair, |λ| ≈ 1.0 dead flat —
exactly what the smoke test measured (1.0000000000000053, ...). So the pipeline is aimed at the
stable branch *by construction*: the only phase it ever seeds is the one that continues into the
encounter-avoiding family.

**(b) Basin asymmetry finishes the job.** `correct_symmetric_fixed_jacobi` is a scalar Newton
solve in x0 at fixed (C, crossing index, ẏ0 sign). Multiple genuinely distinct branches coexist
under one p:q label at nearby (x0, C) — Anderson & Lo's own Table 1 has four "5:6"/"3:4"
branches spanning |λ| = 1.000008 to 2.8e4; `#776` found two topologically distinct "4:3"
branches at nearby C; `#504` found a stable-but-wrong-topology root at Pluto-Charon (the
motivating case written into `deflated_newton.py`'s docstring). The stable near-circular roots
have wide smooth basins; the saddle roots' basins are narrow and "fractal" (`#753`'s documented
x0 ≈ -1.42 "fractal sensitivity" hotspot; `#758`: the bracket/sign-flip scan misses the true
5:6-LO root even at coarse grid resolution *centered on it*, and its direct-Newton basin is only
~±2e-4 wide in x0). A generic seed therefore converges to the stable root with overwhelming
probability even when a saddle root exists nearby.

**(c) Possible C mismatch (unverified, cheap to check).** The harness corrects at the two-body
seed's *own* Jacobi constant. Published saddle families occupy narrow C ranges (AAS 24-288
prints them per family, e.g. 3:4 ∈ [2.9916, 3.0261] at µ = 3.54326e-5). If the seed's natural C
falls outside the target family's existence range for some ratios, *no* x0 basin at that C
contains the saddle and convergence to some other family is guaranteed. I did not run the
numbers (no-code scope); flagged as a contributing factor to check in Step 1 below.

Bottom line: **basin-of-attraction and structural-phase, together** — the `-1` seed is both
geometrically constructed at the stable branch's own configuration and sitting in its wide
basin, while the branch the campaign wants is only reachable from the excluded `+1` phase (or
by walking the family curve to it).

## 4. Fix candidates assessed

### (a) Multi-start / basin diversity around the two-body seed — **weak, mostly disproven**

The project already ran this experiment at scale on a *known-answer* target: `#756` swept
x0 ∈ [-1.9, 1.7], both ẏ0 signs, half-crossings 2-8 (159 converged candidates,
`data/found/756_...`) hunting Jupiter-Europa 5:6-LO — and did **not** find it. It was found
only by `#758`'s *paper-sourced* x-position (Table 2's homoclinic intersection, ±1e-4 window).
Narrow fractal basins defeat grid multi-start exactly for the branch class we want; small
perturbations of the `-1` seed stay inside the stable root's wide basin. So: the target family
generally does *not* surface "somewhere in a wider two-body-seeded basin search" at practical
grid density. The one multi-start variant with a mechanism behind it is surveying the *other
phase* — which is candidate (c1) below, done safely.

### (b) The flyby-rotation seed — **right physics, unproven in the field; second tier**

`two_body_flyby_rotation_seed` + `flyby_rotation_symmetric_seed` encode exactly the missing
ingredient (a genuine close-encounter V∞ rotation — it even starts internally from the
`x0_sign=+1` ellipse and backs off by `safety_margin`). But its test coverage is
geometry-consistency only (V∞ magnitude preserved, turn-angle monotonicity, input validation,
singularity backoff — no test ever converges it onto a known family), and its one field trial
(`#755`, aimed at *known* Table-1 rows) "did not itself locate a Table-1 match in the time
available." Also: its output is non-symmetric and needs a slow propagation to a perpendicular
crossing (documented near-1:1 hazard). Verdict: keep as a secondary seed source; if revived,
first re-validate it on Jupiter-Europa where `#758` now provides the answer key that `#755`
lacked. Do not make it the primary fix.

### (c) Strategies already proven elsewhere in this codebase — **the primary fix lives here**

**(c1) Conjugate-apse ("encounter-phase") seed — new, tiny, mechanism-directed.** The
numerically-safe equivalent of the excluded `+1` seed: take the *same physical ellipse* as
`x0_sign=+1` (r = 1 apse aligned with the secondary) but seed the corrector from its **other
apse**, which sits at x0 = -(2a - 1) on the far side, well away from the singularity, and is a
perpendicular symmetric crossing. ~15 lines next to `two_body_resonant_seed` (grep confirms no
apoapsis-phase construction exists anywhere in `search/`). This targets the encounter-adjacent
branch directly while keeping the harness's safe-seeding convention. Risk: the corrected orbit
at the seed's natural C may still fall into a stable basin — unknown until tried, and the
Oberon control answers it in minutes. While here, also try correcting at the *published* Oberon
per-family C values instead of the seed's natural C (kills mechanism (c) as a confound).

**(c2) Pseudo-arclength fold-turning from the stable branch — the strongest candidate, already
built and proven.** `cr3bp_jacobi_arclength.py` (`#249`) documents in its own docstring *this
exact failure mode and its fix*: "the 1-DOF perpendicular-x-crossing symmetric corrector lands
only on the stable branch; natural-parameter Jacobi continuation diverges at the fold" — and it
recovered the previously-unrecoverable unstable Earth-Moon cycler members (C11a, C21) by
turning the saddle-center fold in (x0, C). `#776`'s 4:3 `FOLD_REVERSAL` family-mixing finding
is direct evidence the same stable/saddle fold structure exists in resonant families too. Under
this fix the naive two-body seed's *reliability at finding the stable branch becomes an asset*:
it is a dependable on-ramp onto the family curve; the arclength walk then traverses the curve in
both directions through folds, classifying Barden |λ| at every member (the module already
computes `nu`/`abs_lambda` per member via the same gauntlet discipline). This *also* fixes the
`#859` note's independent Sec. 4 defect — the current `d_jacobi=5e-4 × 9 steps` default samples
a C-span of 4.5e-3, nowhere near a family's existence range — because an arclength walk
terminated by its own natural boundaries *is* the "spans the existence range" survey `#859`'s
registration asked for. Two fixes for the price of one. Caveat: fold-connectivity between the
recovered stable branch and the target saddle branch is not guaranteed per ratio (disconnected
components and branch points, rather than folds, both occur in PCR3BP family trees) — which is
exactly what the Oberon gate measures. Integration cost: the module is written in the
(k1,k2)-cycler idiom but uses the same `half_crossings`/`ydot0_sign`/`correct_symmetric_fixed_jacobi`
substrate; adapting `stage_a_worker` to call it is glue work, not new numerics.

**(c3) Deflated Newton in x0 at fixed C — the backstop for disconnected branches.**
`deflated_newton.py` (`#524`) was written *for this*: its own docstring names `#504`'s
stable-but-wrong-topology root as the motivating case. Scalar deflation on
`correct_symmetric_fixed_jacobi`'s residual enumerates coexisting roots at a given C without
re-finding the stable one; swept over a modest C grid (published ranges where available) it
catches saddle families that no fold connects to the stable branch. Cheap per solve; narrow
basins still require the deflated Newton to pass nearby, so pair it with a moderate x0 grid.
Note `#648`'s `deflated_variational_periodic_orbit.py` (the `#606` seedless spectral corrector
turned distinct-family enumerator, built + positive-controlled) exists as a heavier third-tier
cross-check if (c1)-(c3) all disappoint — it is the only genuinely seed-free enumerator in the
codebase.

For the record, candidate 2's premise checks out historically: **every** prior
family-confirmation success was published-anchor-guided — `#753`/`#758` Jupiter-Europa used
Anderson & Lo's own C_flyby, Table 1 eigenvalue targets, and Table 2's homoclinic x-position;
`#765` Saturn-Titan used Vaquero's Table 4.1; `#776`/`#777` vendored the Miceli-Bosanac ESM
12-decimal ICs; `#780` Earth-Moon vendored Casoliva Table 3. None bootstrapped table-free.

### (d) mu-continuation from a solved system — **sound, proven machinery, best as cross-check**

`mu_continuation.py` (pseudo-arclength in (x0, C, mu), fold-turning, Radau cross-checked,
topology-jump guarded) already exists and was used to self-discover the Roberts-Tsoukkas
binary-regime cyclers up to mu = 0.5 — a far larger mu excursion than any needed here.
Physically the idea is sound for these targets: Poincaré second-kind resonant families persist
for all small mu, and the relevant spans are modest (Triton 2.09e-4 → Oberon 3.54e-5 is ~6x;
Triton → Titania 3.92e-5 similar; Europa 2.53e-5 → Ganymede ~7.8e-5 ~3x; Rhea ~4e-6 is a ~50x
descent from Triton — larger risk of en-route bifurcation, guarded by the existing
period-continuity/topology-jump gates). Its unique strength: the source member's topology is
*table-verified* (Neptune-Triton ESM rows, Miceli-Bosanac's ~20 resonant labels 1:2…4:7; the
Jupiter-Europa `#758` seeds), so what arrives at the target mu carries verified identity — the
one thing (c1)/(c2) can't fully guarantee. Its limits: per-(p,q) coverage only where a solved
system has that label (nowhere near all 43 coprime pairs), one member per walk (still needs a
C-continuation at target mu afterwards), and eigenvalues must be re-measured at target mu (fine
— that's Stage A's job). Verdict: not the primary engine for a 43-ratio × 4-system grid, but
the ideal *independent verification* of whatever (c1)/(c2) recover on Oberon, and the fallback
primary if fold-connectivity fails broadly.

## 5. Concrete next-step scope (Oberon-gated, in order)

1. **(~0.5 day build, minutes of compute)** Add the conjugate-apse seed (c1) to
   `jovian_resonant_families.py`; on Uranus-Oberon, converge + Barden-classify it at all six
   published ratios, at both the seed's natural C and 2-3 C values inside each family's
   published range (AAS 24-288 p.7-8). Also record whether the naive seed's natural C even lies
   inside the published ranges (kills or confirms mechanism 3c).
2. **(~1-2 days glue + hours of compute)** Wire `cr3bp_jacobi_arclength`-style fold-turning
   into a Stage A' worker (c2): from every converging seed (stable is fine), walk the family
   curve both directions to its natural boundaries, classify |λ| per member, report in-band
   segments. This supersedes the `d_jacobi`/`n_c_steps` recalibration question in the `#859`
   note — do not spend a decision on that knob separately.
3. **(the gate)** Success criterion on Oberon before ANY novel-system dispatch: recover, for
   ≥4 of the 6 published families, an unstable segment whose C-range overlaps the paper's
   printed range AND whose topology matches the label (period ≈ 2πq within the few-% tolerance
   the `#755` reviewer ruling established, winding count = p, close-approach signature
   present). Cross-check ≥1 family by mu-continuing the same-ratio Neptune-Triton
   table-verified saddle down to Oberon's mu (d) and confirming both methods land on the same
   family. If a ratio resists, apply (c3) deflated-Newton at the published C before declaring
   it missing.
4. **Decision rule**: gate passes → run Stage A' on the pilot's 3 novel systems with the fixed
   seeding, misses now honestly stampable (method-conditional, per the negative-results
   registry discipline). Gate fails ≥3 of 6 families after (c1)+(c2)+(c3)+(d) → **shelve Stage
   A** and say so in OUTSTANDING.md: at that point the tractable path runs through `#790`'s
   corrector/alphabet build, not more seeding tweaks, and unfixed-seeding sweeps would only
   manufacture unusable negatives.

Total: ~2-4 days effort, trivial compute, before committing the multi-day sweep. No catalogue
writeback anywhere in this scope.

## 6. Interaction with `#849`

None. `#849` (closed 2026-08-21) fixed a positional-array-reading bug in
`dsm_descriptor_seed._descriptor_params` / `search/descriptor.py` — catalogue-row
`free_return_arcs` posing in the DSM closure lane. It shares no code with
`two_body_resonant_seed`, `correct_symmetric_fixed_jacobi`, `cr3bp_continuation`, or anything
this note recommends touching. Different failure class (data plumbing vs. branch topology),
different lane, already done. This recommendation neither depends on nor conflicts with it.

## 7. What I am not certain about

- **Fold-connectivity** (the load-bearing assumption of c2) is supported by `#249` and `#776`
  precedents but is conjectural per ratio/system; some saddle families may be disconnected or
  joined via branch points an arclength walk passes straight through. The Oberon gate exists to
  measure exactly this.
- The **phase-selects-branch analysis** in Sec. 3a is geometric reasoning grounded in the code
  and the `#758` close-approach evidence, not numerically verified here (no code was run, per
  scope). The conjugate-apse seed could still land stable at some ratios.
- Whether AAS 24-288 prints **eigenvalue-vs-C data precise enough** for a quantitative gate
  (the digest indicates C-ranges and characterization "as a function of Jacobi constant",
  likely figures) — the C-range overlap + topology + close-approach gate above deliberately
  avoids depending on digitized eigenvalues; digitization is an optional tightening.
- The naive seed's natural C vs. published family ranges (mechanism 3c) is unchecked; Step 1
  checks it for free.
- Effort estimates are reading-based (module contracts, `#859`'s measured per-cell costs), not
  prototyped.
