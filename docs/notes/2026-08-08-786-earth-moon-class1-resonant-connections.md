# `#786` — Earth-Moon Class 1 (Casoliva p:q-resonant orbit) homoclinic self-connections

**Task:** `#786`, split from `#785`'s own second avenue (per `#783`'s own note that it does NOT
depend on `#783`'s connection-reproduction result), applying this project's own established
Poincaré-section Newton-shooting connection machinery (`jovian_resonant_connections.py`-style —
the same method already used for Jupiter-Europa `#754`, Saturn-Titan `#767`, and Neptune-Triton
`#781`) to Earth-Moon's own Class 1 resonant orbits (`#780`'s own Table 3, Casoliva et al. 2010).
This project's own resonant-connection machinery had never been pointed at Earth-Moon at all
before this task — genuinely unexplored territory for this specific method.

---

## Verdict (read this first)

**A filename collision, found before writing any code.** The dispatch note asked for a module
named `earth_moon_resonant_connections.py` — but `#783` had already committed a module under that
exact name (its own Class 2/Barrabés-Mondelo-Ollé continuation-of-homoclinic-connections attempt,
a categorically different algorithm). The dispatch note's own author did not know about `#783`'s
filename choice when writing `#786`'s instructions. Rather than silently overwrite `#783`'s
committed work, this task uses `src/cyclerfinder/search/earth_moon_class1_resonant_connections.py`
instead — the same disambiguation an `advisor()` review call caught and recommended before any
code was written.

**CLEAN NEGATIVE on the connection search itself.** `build_node` works cleanly for both targeted
orbits (`7-3b`, `7-3c`) — eigenvalue cross-checks against `#780`'s own `table3_gate_report`
confirm to `<3e-10` relative, comfortably inside this task chain's own ~50-2000 Newton-tractable
`|lambda|` band. But an extensive, multi-strategy search for a genuine homoclinic self-intersection
(`Wu(orbit) ∩ Ws(orbit)`) never converged for either orbit. Every near-converged Newton attempt
plateaus at, or very near, one of the orbit's own 20 natural `{y=0}` section points — the
diagnostic signature of an unusually dense crossing structure (this orbit crosses `{y=0}` 20 times
per period, vs. 2-6 for every prior sibling target) starving the correction of a genuine
transversal-intersection basin within the region searched (`k<=34`, ~1.5-3.5 periods per leg).
This is reported as an honest, well-diagnosed negative — not forced, not fudged, no tolerance
loosened.

**Literature-novelty finding (the load-bearing gate result of this task, independent of the
connection-search outcome).** Vaquero's 2013 Purdue dissertation — already in this project's own
corpus, `#765`'s own primary source, but only its Saturn-Titan chapter (Sec. 4.3) was ever
digested — has an entire, never-previously-read Earth-Moon chapter (Sec. 4.4, "Resonance
Transition in the Earth-Moon System", pp.133-172) that explicitly computes GENUINE homoclinic
self-connections of planar Earth-Moon p:q resonant orbits: **1:2 and 2:3, at `C=2.8284`** (Sec.
4.4.1, Figs. 4.22-4.24), via her own family-continuation orbit-generation lineage — a categorically
different method from both Casoliva's Class 1 (elliptical/second-species) and Class 2
(L1-Lyapunov). **This means the underlying PHENOMENON — a genuine, Newton-tractable homoclinic
self-connection of an Earth-Moon p:q resonant orbit — is established prior art, not novel.**
This task's own specific target (Casoliva's `7-3b`/`7-3c`, `C~1.0687`) was never computed by
either Vaquero or Casoliva — had a genuine hit been found here, it would have been reported as
corroboration-and-extension of Vaquero's own demonstrated phenomenon to a new specific resonance,
never as a novel discovery, per the dispatch's own explicit instruction ("report the result as
reproduction/corroboration, not novel discovery — that's still a legitimate, useful result"). Since
the search itself is a clean negative, this question is moot for THIS task's own object — but the
corpus-completeness gap this finding exposed is real and independently useful, registered as `#787`.

---

## 1. Target selection

`#780`'s `table3_gate_report` (all 16 Casoliva 2010 Table 3 rows) was re-derived directly this task
to identify candidates. Filtering for (a) `exists_in_em_system=True` (physically real, does not fly
through the Earth), (b) `satisfies_resonance=True` (genuinely satisfies its own labelled p:q
resonance), and (c) a REAL in-plane saddle (`|k_signed| > 2`, `k_signed == k_par`):

| Row | k_signed | lambda (approx) | exists_in_em_system | satisfies_resonance | Chosen? |
|---|---|---|---|---|---|
| `1-2d` | 4.8573 | ~4.64 | True | True | No — below this chain's own tractable band |
| `2-1b` | 2.0374 | ~1.19 | True | True | No — essentially marginal instability |
| `2-1c`/`2-1d` | 103-119 | ~103-119 | **False** | (mixed) | No — unphysical, flies through Earth |
| `7-3b` | 57.356 | ~57.34 | True | True | **PRIMARY** |
| `7-3c` | 57.043 | ~57.02 | True | True | **SECONDARY** |

`7-3b`/`7-3c` sit squarely inside the ~50-2000 `|lambda|` band this whole task chain has repeatedly
found Newton-tractable (`#766`'s `C=3.0041`: `|lambda|~54.6`; `#781`'s 4:5-saddle: `|lambda|~105`).
They share the same Jacobi constant and period (`C~1.0687`, `T~18.8496`) but are confirmed DISTINCT
periodic orbits (different `(x0, xdot0, ydot0)` triples, different section-crossing sequences,
`#780`'s own printed IC differs). `1-2d`/`2-1b` are registered as a follow-on (see Sec. 6), not
attempted this task.

## 2. A new geometric wrinkle: no perpendicular crossing at all

Every prior sibling module's target orbit was built via `correct_symmetric_fixed_jacobi` — a
genuine perpendicular `{y=0, xdot=0}` crossing IS the orbit's own IC by construction. Casoliva's
own Class 1 orbits are recovered via `#780`'s `recover_table3_row`, a GENERAL full-state Newton
corrector (`correct_periodic`, 6 free state components + period, no perpendicularity constraint)
seeded from Casoliva's own printed generic (non-perpendicular) Poincaré-crossing IC — her own text:
"we have not used this property [perpendicular crossings] in this paper." Direct inspection this
task (`rtol=atol=1e-13`) confirms `7-3b` crosses `{y=0}` **20 times per period**, at NO
perpendicular point whatsoever — the closest approach is `|xdot| ~= 0.2587` (row 3 of the crossing
table, `t~=2.987`). This is a materially richer and qualitatively different crossing structure than
any prior sibling orbit (Jovian: 2/period, both perpendicular; Saturn-Titan: 4/period, 2
perpendicular; Neptune-Triton: 6/period, 2 perpendicular + 2 non-perpendicular mirror pairs).

Two consequences, both handled explicitly in the new module (not silently assumed away):

1. **The corrector's own converged IC does not sit exactly on `{y=0}`** (`correct_periodic`'s
   min-norm 7-unknown/6-equation Newton step has no `y0=0` constraint — the converged `7-3b` IC
   lands at `y0 = -1.60e-05`). A new `_snap_to_y0` helper phase-shifts the converged IC onto the
   section exactly (`dt = -y/vy`, iterated to convergence) — re-verified this task: periodicity
   closure at the snapped IC is `3.6e-11`/`1.5e-09` for `7-3b`/`7-3c`, far tighter than `#780`'s
   own `tol=1e-8` gate, confirming this is a phase re-parametrization of the SAME solution.
2. **No single-sign-restricted section, nor even Neptune-Triton's own 4-combo union, covers a
   fundamentally different crossing density** — the same 4-combo union (`#781`'s own fully-
   unrestricted convention) is reused directly here, and confirmed this task to recover the full
   20-point crossing set with no double-counting.

`ResonantNode` is NOT built via `ResonantNode.from_candidate` (that classmethod hardcodes
`state0=[x0,0,0,0,ydot0,0]`, shaped to a perpendicular-crossing candidate). `build_node` in the new
module constructs the node directly from the snapped full 6-state, reusing
`_planar_floquet_pair` for the eigenvector derivation (the same routine `from_candidate` itself
uses internally) — no change to shared machinery needed.

## 3. Manifold offset: raised from the reused default, and why

Every prior sibling module reused `jrc.ANDERSON_LO_EPSILON` (`0.5e-5`) verbatim. This module uses
`EPSILON=1e-4` (20x larger) instead, a deliberate and documented deviation. A direct diagnostic
(propagating a single unstable-manifold seed forward, tracking `(k, ghost-distance)` pairs) found:
at `epsilon=0.5e-5`, the manifold does not clear `GHOST_GUARD_DELTA` until `k~13-16`; at
`epsilon=1e-4`, it clears comfortably by `k~10` (ghost distance `5.7e-3`, six times the guard). The
smaller default combined with this orbit's minimum inter-crossing time gap (`~0.0686` nondim time
units) would force an uncomfortably tight margin against `correct_connection`'s own hardcoded
`max_step=horizon/500` crossing-detection resolution at the short horizons a smaller epsilon would
require. A direct resolution check (counting `{y=0}` events at `horizon/500` vs. `horizon/5000`
over 2- and 3-period horizons) found IDENTICAL counts both resolutions (41 events/2 periods, 61
events/3 periods) — confirming `correct_connection`'s own default resolution is safe at the
`max_time_factor<=3-4` this task actually used. (A further `horizon/20000` check found one fewer
event over the 2-period horizon — an unresolved minor discrepancy, most plausibly a near-tangential
root-finding artifact at that much finer step count, not evidence against the resolution this
task's own scans actually rely on.) `GHOST_GUARD_DELTA` itself is UNCHANGED (`1e-3`, identical to
every sibling module) — only the manifold offset that reaches it sooner was raised.

## 4. Search methodology and the honest negative

Four independent strategies were used, in increasing order of sophistication, as each prior one
came up empty:

**(a) Blind diagonal `(k_u=k_s, branch_u, branch_s)` grid+Newton scan.** `k=10..20` (11 values,
except one batch ran `10..14`), all 4 branch-sign combinations, `scan_n=6-8` internal grid,
`tol=1e-7`, `epsilon=1e-4`. 40 total `correct_connection` calls. **Zero converged.** Every
near-converged residual (as low as `8.2e-6`) landed at `(x,xdot)` values matching one of the
orbit's own 20 natural section-point crossings to 5-6 significant figures — a ghost, not a
homoclinic intersection.

**(b) Direct manifold-tube pre-scan (a cheaper, independent cross-check).** Rather than blind
Newton, this samples the unstable manifold (24 `tau` values x 2 branches) forward for 3 periods and
the stable manifold backward for 3 periods, collects ALL `{y=0}` crossings with their `(tau, branch,
k, x, xdot)` label (no Newton — pure propagation, ~22s total vs. tens of minutes for (a)), filters
out points within `5x GHOST_GUARD_DELTA` of the orbit's own 20 points, and ranks the remaining
`(unstable, stable)` pairs by raw Euclidean `(x,xdot)` distance. For `7-3b`, the closest such pair
sits `3.29e-4` apart (genuinely close, not epsilon-scale); several more sit `3.7e-4` to `5.2e-4`
apart. For `7-3c`, the closest is `2.36e-4`.

**(c) Targeted Newton refinement seeded exactly at (b)'s own best candidates.** Seeding
`correct_connection` with `tau_u0`/`tau_s0` from the closest non-ghost pairs (no internal grid
scan needed): for `7-3b`, the search plateaus after 1-25 iterations at residuals `1.7e-5` to
`1.85e-4`, never reaching `tol`. Measuring the ACTUAL matched crossing this seed resolves to
(`correct_connection`'s own `_section_crossing` k-indexing, not this task's own independent
tube-scan k-count) found it landed MUCH closer to one of the orbit's own 20 points than the raw
propagated seed did — `5.37e-5` ghost distance for the `7-3b` best candidate (well inside the
guard), `2.39e-3` for `7-3c`'s (just outside the `1e-3` guard, inside the pre-scan's own `5e-3`
filter margin). This is reported exactly as measured: the original working hypothesis (a genuine
non-ghost fold/tangency) is CORRECTED here — the more likely explanation is that this task's own
from-scratch tube-scan crossing count and `_section_crossing`'s own internal count are not
perfectly aligned at these high `k` values (a discretization/counting mismatch, not a dynamical
fold), and the true object `correct_connection` searches near is closer to a ghost than the
pre-scan's own raw-point filter indicated.

**(d) `fd_step` sweep (ruling out finite-difference noise).** Sweeping `fd_step` over `1e-6` to
`1e-4` (two orders of magnitude) at the `7-3b` best candidate reproduces the IDENTICAL
non-improvement (`n_iter<=2`, same residual to high relative precision every time) — the plateau is
not a finite-difference artifact.

**Both `7-3b` and `7-3c` show the same qualitative behaviour** across `k=10-34` and both targeted
and blind seeding strategies: this orbit's own unusually dense (20-crossing/period) section
geometry appears to make the trivial ghost solution the dominant nearby attractor for
`correct_connection`'s own local 2-D Newton search across the region explored, starving it of a
genuine transversal-intersection basin. This is reported as a clean, well-diagnosed negative on the
SEARCH conducted (not a proof that no genuine self-intersection exists anywhere for this orbit —
per this project's own discipline, "no X found" is conditional on the search formulation, never
proof of absence).

## 5. Literature novelty gate

**Step 1 — Casoliva's own text (load-bearing for the METHOD question).** Confirmed directly
(module docstring, `#780`'s own digest): Casoliva's Class 1 method (elliptical/second-species
differential correction) never touches manifolds or homoclinic connections at all — Class 2's own
Barrabés-Mondelo-Ollé continuation is the ONLY homoclinic-connection content either Casoliva paper
contains, and that is scoped entirely to the L1-Lyapunov He1 family (`#783`'s own already-negative
target).

**Step 2 — corpus-index grep + digest review (the load-bearing finding of this task).**
`docs/notes/CORPUS_INDEX.md` was grepped for `vaquero|anderson.*lo|koon.*lo.*marsden|resonant.*
manifold|homoclinic`, surfacing Lo & Parker 2004 (already in corpus) and Parker-Davis-Born 2010
(already in corpus) as the closest candidates. Direct text-layer reading of Lo & Parker 2004's own
Sec. III.A ("Heteroclinic Transfers") confirmed it computes ONLY heteroclinic transfers BETWEEN
DIFFERENT orbit families (its own Family F, "resonant lunar flyby orbits", the direct
classification-scheme predecessor of Casoliva's Class 1, only appears in a heteroclinic A-F1/F1-A
transfer to a DIFFERENT family, never a homoclinic self-connection of F1 itself). Parker-Davis-Born
2010 similarly chains heteroclinic connections between different orbits' patch points.

A subsequent WebSearch pass ("resonant orbit Earth Moon invariant manifold homoclinic
self-intersection Poincaré section") surfaced Vaquero & Howell 2014 (*Acta Astronaut.* 94:302-317)
— already flagged in this project's own acquisition backlog (item 85, paywalled, unacquired) as
closely related to Vaquero's own 2013 dissertation. Checking the dissertation directly (already
fully acquired in this project's corpus) found the load-bearing content: **Sec. 4.4.1 ("Planar
Natural Transfers: Resonant Homoclinic and Heteroclinic Connections"), pp.134-138**, explicitly
states: *"If the intersection corresponds to a stable and an unstable manifold associated with the
same periodic orbit, then a homoclinic connection can be calculated"* and computes exactly this for
BOTH a `1:2` and a `2:3` planar Earth-Moon resonant orbit at `C=2.8284` (Figs. 4.22-4.24), via her
own family-continuation orbit-generation lineage (Fig. 3.17), plus a heteroclinic connection
between them. Sec. 4.4.7 ("Earth-Moon Periodic Cyclers") separately catalogues STABLE `2:1`/`3:1`
resonant cyclers (a family-continuation study, not a homoclinic-connection one) — neither section
touches `7:3` or Casoliva's own elliptical/second-species generating-family method.

This is a genuine corpus-completeness gap: Vaquero 2013's own Ch. 4 has TWO independent system case
studies (Saturn-Titan Sec. 4.3, already digested for `#765`; Earth-Moon Sec. 4.4, never digested
until this task needed it). Registered as `#787` (see Sec. 6).

**Step 3 — `search/literature_check.py`'s own mandatory-floor gate**, run with real WebSearch
results fed in (8 queries per `MAX_QUERIES=8`, signature
`CandidateSignature(primary="Earth", sequence=("Moon",), resonances=("7:3",),
topology_label={"resonant"})`). Verdict: **`status="published"`, `confidence=0.75`**, citing "A
FREE-RETURN EARTH-MOON CYCLER ORBIT FOR AN INTERPLANETARY CRUISE SHIP" (an Uphoff-Crouch-style
free-return cycler — a completely different object class, patched-conic, no CR3BP resonant-orbit
manifold content whatsoever). This is a KNOWN FALSE POSITIVE, consistent with the module's own
documented scope limitation ("NOT FOR RAW (non-cycler) CR3BP PERIODIC ORBITS... over-matches on
bare body-name tokens") — every prior sibling task in this chain has hit an analogous limitation of
this tool for raw resonant-orbit-manifold objects. Reported here per the mandatory-floor
instruction, but NOT treated as informative evidence either way; the load-bearing evidence is Step
2 above.

**Verdict for THIS task's own object**: moot, since the search itself is a clean negative — no
"novel" or "corroboration" claim is made for a connection that was not found. Had a genuine hit
been found, Step 2's finding would have required framing it as corroboration-and-extension of
Vaquero 2013's own established phenomenon, never as novel discovery.

## 6. Registered follow-up

**`#787`** — a dedicated digest of Vaquero 2013's Earth-Moon chapter (Sec. 4.4 in full, not just
the 4.4.1 finding this task needed), cross-checked against Casoliva's own Table 3/`#780` and this
task's own `7-3b`/`7-3c` work for any direct numeric overlap, plus a citation-mining pass per the
standing corpus-document-policy discipline. Not dispatched.

`1-2d`/`2-1b` (the two weaker-eigenvalue Table 3 rows, `lambda~4.6`/`1.2`) remain untried — outside
this task chain's own demonstrated ~50-2000 Newton-tractable band, a genuine attempt would need a
much longer scan horizon than this task's own budget allowed. Not registered as a separate task
number (folded into this note as a documented gap, since it is a narrow within-scope extension of
this task rather than a new avenue).

## 7. Verification

1. `uv run pytest tests/search/test_earth_moon_class1_resonant_connections.py -v` — 14 passed
   (~3.6 min wall time).
2. `uv run pytest tests/data tests/search -q` — full ratchet suite. 2 pre-existing failures
   (`tests/search/test_eggie_ballistic.py::test_gate_b_table4_vinf_reached_but_subsurface`,
   `tests/search/test_504_pluto_charon_kk_sweep.py::test_504_sweep_33`), both topically unrelated
   to this task (Europa/Pluto-Charon), in files this task never touched, and independently observed
   this session as already under active investigation by a concurrent sibling task. Numerous
   `XPASS`/`XFAIL` entries pre-date this task (documented cross-platform DOP853/BLAS
   non-bit-reproducibility class, `#584`/`#631`/`#632`/`#635`/`#731`/`#782`/`#784` precedent).
3. `uv run ruff check .` / `uv run ruff format --check .` — clean on both new files.
4. `uv run mypy src tests` — 837 source files, clean.
5. `tests/data/test_outstanding_structure.py`/`test_outstanding_header_body_consistency.py` —
   both pass after the `#786`/`#787` `data/OUTSTANDING.md` edits.
