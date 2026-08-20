# `#839`: targeted Wu(3:1) <-> Ws(2:1) Earth-Moon heteroclinic search AT C=3.13 — GENUINE CONNECTION FOUND

**Task:** `#839`, registered 2026-08-12 (found during `#828`, gated on `#827`). `#827` digit-grade
reproduced Kumar, Rawat, Rosengren & Ross (2026)'s own SEVEN printed Table-5 `Wu(3:1) -> Ws(2:1)`
Earth-Moon CR3BP heteroclinic states (`C_J in {2.54, 2.70, 2.86, 3.00, 3.05, 3.10, 3.15}`), but
`vaquero-31-c313-em-resonant-po-2013` sits at **C=3.13**, which is not one of those seven printed
rows (`keh.kumar_table5_state6(3.13)` raises `KeyError`, `#827`'s own dedicated test). This task is
a **fresh search**, not a reproduction: is there a genuine `Wu(3:1) <-> Ws(2:1)` connection at
THIS row's own Jacobi constant, seeded from this project's own converged nodes, with no printed
digit target to match against?

**Why it matters:** `vaquero-31-c313-em-resonant-po-2013` is classed `resonant_po` (not `cycler`)
specifically because `#453`/`#811` found "no demonstrated transport utility" for it. A genuine
transport connection touching this row's own orbit is evidence relevant to that `orbit_class`
question. This task is the **compute half only** (per the `#822`/`#828` and `#827`/`#854` split
precedent) — no catalogue writeback here. **New follow-up `#855` registered** for the
`orbit_class` adjudication.

**Module:** `src/cyclerfinder/search/vaquero_c313_targeted_search.py` (see its own docstring for
the full method account: mu identity, node construction, the search strategy, and the epsilon
calibration).
**Data:** `data/found/839_c313_targeted_search/results.json` (both nodes' full provenance, both
independent search runs' full verification batteries, and the leg-minimum-radius diagnostic).
**Evidence tests:** `tests/search/test_vaquero_c313_targeted_search.py` (11 tests, NOT slow,
~50s) — from-scratch node re-derivation and two independent from-scratch connection
reconstructions, seeded only by recorded phase indices.

---

## Result: a genuine, independently-verified connection exists at C=3.13

**mu:** the catalogue row's own `mass_ratio` is THIS PROJECT'S registry Earth-Moon mu
(`0.01215058439469525`), NOT Kumar's own printed mu (`1.2150584270572e-2`, differing at the
~1.24e-10 absolute level) — every node and manifold here is built at the project's own registry
mu, matching the catalogue row exactly.

**Node31** (the 3:1 orbit) **is the catalogue row itself**: re-converged from its own recorded
`state_nd`/`period_nd` with the same fixed-Jacobi symmetric corrector `#799`/`#811` used
(`half_crossings=3`), landing back on the catalogue's own digits to ~1e-13 (well under the
`NODE31_IC_ABS_TOL=1e-8` gate).

**Node21** (the 2:1 orbit) has no catalogued or Vaquero-sourced counterpart at C=3.13 — Vaquero's
own 2:1 family (`VAQUERO_C_RANGE_21`) tops out at C=2.66. Kumar's own Table 6 prints 2:1 rows at
C=3.10 and C=3.15, a wider continuation of the SAME family; node21 was built by re-converging
those two prints at THIS project's own mu (an external bracket, not trusted raw) and
step-continuing 0.001-Jacobi-at-a-time from the C=3.10 anchor to C=3.13 (**a 0.01 step was tried
first and measured to jump onto a different branch** — period jumped 6.8 -> 11.8 -> 6.1 -> 3.5 nd
across three 0.005-scale steps; the module docstring records this). The final C=3.13 member lands
within 0.0037 of a linear interpolation of the two own-mu bracket endpoints (well under the 0.02
continuity-check tolerance).

**Node21 identity cross-check (advisor-caught, resolved same session):** the C=3.10 bracket
anchor was checked directly against `#827`'s own `build_kumar_node(kumar_system(), 2, 3.10)` at
Kumar's own mu — period `6.798978144793446` and `|lambda|=333.04100533262795` there, vs
`6.7989781409144125` / `333.0410072133632` here (agreement to ~1e-8 relative despite the
~1.24e-10 absolute mu difference). **node21 IS the same physical Kumar 2:1 family member `#827`
already reproduced**, not a wrong branch — this was the one construction step with no direct
catalogue anchor, and it is now anchored.

### The search

`#822`'s unchanged Poincare-section Newton connection machinery
(`vaquero_em_cycler_connections`), direction `Wu(3:1) -> Ws(2:1)` matching Kumar's own Table-5
convention, `epsilon=KUMAR_EPSILON=1e-4` (reused verbatim from `#827` for the same reason: the
2:1 node's saddle here is extreme, `|lambda|` in the hundreds, which amplifies the sibling default
epsilon past the forward-reapproach gate). UNSEEDED — no printed state exists at this C to select
which crossing to converge, so every same-ydot-sign close approach between the two manifolds'
whole crossing sets is a candidate, closest first.

Run **twice**, at two independent phase-grid resolutions:

| n_tau | n_seeds | n_refined | n_converged | k_u | k_s | branch_u | branch_s | full_state_gap | radau_gap | passed |
|---|---|---|---|---|---|---|---|---|---|---|
| 48 | 738 | 1 | 1 | 17 | 14 | -1 | +1 | 1.653e-08 | 1.597e-08 | **true** |
| 64 | 1078 | 1 | 1 | 20 | 13 | -1 | +1 | 7.674e-08 | 5.181e-08 | **true** |

Both runs' FIRST refined candidate (closest seed) both converged and passed the full, unmodified
`#822` verification battery — full planar-4-state gap, ydot-sign hard gate, ghost guard,
independent-Radau re-derivation, forward/backward re-approach, Jacobi drift on both legs — with
excellent margins (worst gate: forward-reapproach at n_tau=48, 8.0e-5 vs a 0.5 ceiling; every
other gate 3-5 orders of magnitude under its ceiling). Full per-gate values are in
`results.json`.

**Honest framing of the two runs:** they land on DIFFERENT specific manifold crossings (different
`k_u`/`k_s`, different physical crossing locations) — the two runs are not claimed to be the same
object, and neither run enumerated the full intersection set. What is established is **existence**,
corroborated independently twice: two different phase-grid constructions each find *a* connection
that both converges and passes the same unmodified battery. This is stronger than a single hit,
but weaker than "the" connection at this C.

### Physical character (data for `#855`, not a gate here)

Both legs were densely re-sampled (DOP853, rtol=atol=1e-13, 20000 points) to find the minimum
geocentric and selenocentric radius along the full transfer:

| n_tau | transfer time (nd) | min geocentric (km) | min selenocentric (km) |
|---|---|---|---|
| 48 | 45.78 | 57935.3 | **46247.6** |
| 64 | 46.98 | 57941.7 | **46168.2** |

Both transfers pass through a selenocentric radius of ~46,200 km — **well inside** the lunar SOI
(66,182.9 km) and even inside the Hill radius (61,524.1 km), unlike the catalogue row's own
periodic orbit (periselene 66,995.2 km, 1.2% OUTSIDE the SOI, the reason it is currently
`resonant_po`). The heteroclinic TRANSFER trajectory itself dips much closer to the Moon than the
row's own orbit does. This is reported as raw physical data for `#855`'s adjudication — it is not
itself a claim about `orbit_class` (SOI entry along a manifold leg is not the same question as
SOI entry of the periodic orbit itself, and `#855` will need to work out which one the schema's
"transport utility" criterion actually asks about).

## `#855` — new follow-up registered

Adjudicate whether this genuine, independently-verified `Wu(3:1) <-> Ws(2:1)` connection —
touching `vaquero-31-c313-em-resonant-po-2013` at its OWN Jacobi constant via its OWN orbit as
node31 — constitutes "demonstrated transport utility" under schema v4.9/`#453`'s criterion, and
if so whether `orbit_class` should move `resonant_po -> cycler` (and how that interacts with
`#811`'s own SOI-marginal boundary call on the periodic orbit itself, given the transfer leg's
much deeper lunar approach found above). No catalogue writeback performed in this task.

## Literature-novelty

Not re-run live here (no catalogue writeback) — `#822`'s own live-WebSearch mandatory-floor run
against this exact paper's connection concept already returned `published`, confidence 0.95. A
genuine hit at C=3.13 is a new point in an already-published family, not a new physical claim;
nothing here is claimed novel.

## Verification

`tests/search/test_vaquero_c313_targeted_search.py` (11 tests, NOT slow, ~50s) —
`uv run pytest tests/search/test_kumar_em_resonant_heteroclinics.py
tests/search/test_vaquero_c313_targeted_search.py -q` also run together, clean. `ruff check .` /
`ruff format --check .` clean on the touched files. `uv run mypy src tests` clean except one
pre-existing, unrelated `pypdfium2` stub gap in `src/cyclerfinder/verify/ocr.py` (untouched by
this task, present before it). No `data/catalogue.yaml` change in this task, so the full
`tests/data tests/search -q` ratchet was not required by the dispatch.
