# #630: `perimoon_passage.py` audit against existing bend-gate negatives and admitted (k1,k2) binary-cycler rows

Task #630 (dispatched 2026-07-18). Audit only — no search, no corrector build, no
`data/catalogue.yaml` writeback. Checks whether `src/cyclerfinder/search/perimoon_passage.py`
(built for #627's Titan-continuation pilot, validated once against the admitted Pluto-Charon
(3,2) cycler, and never applied anywhere else) changes or refines the characterization of two
existing bodies of work where an encounter-quality judgment was already made by a different
method.

## What `perimoon_passage.py` actually computes

`find_perimoon_passage(system, state0, period, secondary_radius_km, ...)` propagates a **full
period of a planar CR3BP symmetric periodic orbit** (DOP853, dense output), coarse-samples +
locally refines to the global minimum distance to the secondary (fixed at `(1-mu, 0)` in the
rotating frame), and returns the periapsis distance/altitude and the rotating-frame relative
speed there. Its required inputs are a `cyclerfinder.core.cr3bp.CR3BPSystem` (a two-primary
restricted-three-body system with a single mass ratio `mu`), a periodic `state0`, and a `period`.
This is a **restricted-three-body, single-secondary, full-period-propagation** tool — not a
patched-conic/Lambert flyby evaluator.

Its own validation case (`tests/search/test_627_titan_pilot.py::test_627_perimoon_passage_on_admitted_pluto_charon_cycler`)
runs it against the admitted Pluto-Charon (3,2) cycler (`ross-rt-pc-cycler-32-2026`), asserting a
small positive altitude (not below Charon's surface, not a distant non-encounter) and refinement
stability. The test uses a slightly different Jacobi constant (C=3.579222016200) than the row's
final catalogue-admitted value (C=3.57951501972907, an ~2.9e-4 difference — likely an earlier
`nu=0` window estimate from the #494/#627 development process, superseded by #505's refined
midpoint). This doesn't affect anything below since Target 2 re-runs the check at the *exact*
catalogue IC, not the test fixture's.

## Target 1 — bend-gate-limited negatives (#607, #609, #571 via #625's certification)

Read `#625`'s bullet and the `data/empty_regions.jsonl` entries for `saturn-titan-{mimas,enceladus,
tethys,dione}-analytically-empty-571`, `mars-phobos-deimos-symmetric-closure-609-2026-07-16`, and
`smallbody-multimoon-symmetric-closure-mass-limited-607-2026-07-16`.

**These are architecturally incompatible with `perimoon_passage.py` as written.** All of #607/#609/
#571's candidates come from the `#563`/`#600` circular-coplanar point-mass-primary + patched-conic
Lambert construction (see `search/physical_sanity.py`, `core/flyby.py::max_bend`) — a genuinely
different genome from the restricted-three-body CR3BP periodic orbits `perimoon_passage.py` was
built for. There is no `CR3BPSystem`/periodic `state0`/`period` for a 2-or-3-moon Lambert-arc tour
candidate; running `find_perimoon_passage` against them literally requires an entirely different
tool (a patched-conic hyperbola periapsis evaluator), which does not exist in this codebase and is
out of the scope of "a small adapter."

More importantly: applying `perimoon_passage.py`'s *underlying question* — is the passage a
genuine encounter (not a near-collision, not a meaninglessly distant pass) — to how these entries'
own bend gate is actually computed shows the question is **already structurally foreclosed**.
`core/flyby.py::max_bend(mu_planet, rp_min, vinf)` evaluates the *maximum achievable* deflection at
`rp_min = radius_eq_km + safe_alt_km` — i.e., at the **closest safety-permitted periapsis**, not at
whatever periapsis a specific constructed trajectory happens to use. Two direct consequences,
confirmed against the actual `empty_regions.jsonl` numbers:

* **Not a collision, by construction**: `rp_min` uses a sourced `safe_alt_km` floor (e.g. Mimas/
  Enceladus/Tethys/Dione and Phobos/Deimos radii + safe altitudes from `core/satellites.py`), so
  the assumed periapsis is never below the surface.
* **Not a distant/degenerate pass, by construction**: `rp_min` *is* the closest distance the gate
  ever considers — the single most-favorable-to-feasibility geometry. #571's own Saturn rows
  document exactly this: e.g. Titan-Mimas's `bend_ceiling_vinf_kms=0.429` (the V∞ at which 5° bend
  is achievable at `rp_min`) vs. `min_achievable_vinf_kms=4.5425` (the real Hohmann-floor V∞, 10.6x
  higher) → `bend_at_min_achievable_vinf_deg=0.0466`. The bend gate already assumes the tightest
  possible pass and still fails because the real achievable V∞ is too high — never because the
  construction let the passage happen far away or below the surface.

So for Target 1: **the answer is "bend was the sole issue," not "additionally a degenerate
passage."** `perimoon_passage.py`'s distinctive concern (verifying an *actual propagated*
trajectory's periapsis isn't a collision or a vacuous distant pass) doesn't surface anything new
here, because the gate these negatives were built on already assumes the single closest,
non-collisional periapsis by construction — there is no alternate, worse-or-better periapsis for a
propagated trajectory to reveal. No adaptation of `perimoon_passage.py` was built (none would be
mechanical — it would be a new tool, not an adapter), and no `empty_regions.jsonl` entry needs a
follow-up characterization change.

## Target 2 — #494/#549 admitted (k1,k2) binary-cycler family members

`#494` admitted 5 rows (`e7bca1b`, 2026-06-30, closing `#315`/`#252`/`#255`): 4 abstract-mu Table-I
representatives from Ross & Roberts-Tsoukkas 2026 (`ross-rt-mu001-cycler-11-2026` mu=0.001 (1,1),
`ross-rt-mu01-cycler-32-2026` mu=0.1 (3,2), `ross-rt-mu03-cycler-31-2026` mu=0.3 (3,1),
`ross-rt-mu05-cycler-11-2026` mu=0.5 (1,1)) plus the one real physical system,
`ross-rt-pc-cycler-32-2026` (Pluto-Charon, mu=0.10876473603280369, V2).

**`perimoon_passage.py` runs against this genome natively — no adaptation needed.** Each row's
`orbit_elements.cr3bp` block is exactly a `(mu, state_nd, period_nd)` CR3BP periodic-orbit
specification, i.e. precisely `find_perimoon_passage`'s expected input (`system, state0, period`).

**The 4 abstract-mu rows correctly never claim a physical "useful encounter" in the first place** —
their `orbit_elements` blocks explicitly read "Keplerian elements inapplicable" / "no real bodies,"
`lunit_km`/`tunit_s` are nominal 1.0 placeholders, and `vinf_kms_at_encounters` is null with a note
that the conserved quantity is the Jacobi constant, not a patched-conic V∞. There is no real
secondary radius to check a periapsis altitude against, so there is nothing for
`perimoon_passage.py`'s altitude/collision check to add — this is correctly out of scope, not a
gap. Admission for these rows rests on `winding_topology`'s `reaches_secondary`/`(k1,k2)` check
(a purely topological crossing-count test that needs only `mu`+`state0`+`period`, no real body
scale) plus corrector convergence + Barden stability + independent-Radau crosscheck
(`tests/search/test_ross_rt_2026_mu_family.py`) — never a real-encounter-geometry claim, so nothing
was "assumed by construction" that needed checking.

Running `find_perimoon_passage` against all 4 abstract rows anyway (nondimensional units only,
`secondary_radius_km=0`, informational — no catalogue claim depends on this) confirms the
topological admission is backed by an actual close geometric approach in every case, not just a
crossing-count artifact:

| row | mu | topology | reaches_secondary | r2 (nd, P1-P2 units) |
|---|---|---|---|---|
| mu001-cycler-11 (1,1) | 0.001 | (1,1) | True | 0.00177 |
| mu01-cycler-32 (3,2) | 0.1 | (3,2) | True | 0.02659 |
| mu03-cycler-31 (3,1) | 0.3 | (3,1) | True | 0.12120 |
| mu05-cycler-11, fundamental T1=T/3 | 0.5 | (1,1) | True | 0.01969 |

(The mu=0.5 row's *published* period is the 3rd iterate; run at the full published T the topology
correctly resolves to (3,3), matching the catalogue's own note — checked directly, not assumed.)

**The one real-body row, Pluto-Charon, was independently checked — both previously (approximately,
by #627's own validation test) and now (exactly, at the catalogue's admitted IC, in this audit).**
Running `find_perimoon_passage` at the *exact* catalogue-admitted state
(`x0=-0.693198287043369, C=3.57951501972907, T=11.8334625170346`, Jacobi self-consistency checked
to 2.7e-15) against the real Charon radius (606.0 km, `core/satellites.py`) gives:

```
topology (3,2), reaches_secondary=True
r2_km = 1086.76 km   (altitude 480.76 km above Charon's 606 km mean radius, NOT below surface)
speed_rel_kms = 0.382 km/s
```

This is a genuine, real-scale close passage — not a collision, not a distant non-encounter. So for
Target 2: **"useful encounter" was not merely assumed by construction for the one row where the
concept applies; it holds up under `perimoon_passage.py`'s own quantitative check, run here at the
row's exact admitted IC** (not an approximation). For the other 4 rows, no physical claim was ever
made that needed checking, and the geometric topology admission is backed by genuine close nd-scale
approaches in every case.

## Bottom line

Both targets: **nothing changes.** No adaptation of `perimoon_passage.py` was required for Target
2 (native fit); Target 1 is architecturally out of the module's scope (different genome — patched-
conic Lambert tours, not CR3BP periodic orbits) and, more substantively, the existing bend gate's
own construction (evaluate at the closest safety-permitted periapsis) already structurally
forecloses the "degenerate passage" question `perimoon_passage.py` would otherwise raise. Existing
characterizations for `#607`/`#609`/`#571`/`#625` and for `#494`/`#549`'s admitted rows all hold.
No catalogue or `empty_regions.jsonl` edits made. No new reusable code was written (this was
read-only analysis plus ad hoc verification snippets run inline, not committed as a module), so no
new tests are needed and the ruff/pytest ratchets were not re-run (nothing in `src/`/`tests/` was
touched — confirmed via `git status`).
