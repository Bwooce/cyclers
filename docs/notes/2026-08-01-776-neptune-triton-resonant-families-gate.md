# `#776`: Neptune-Triton planar CR3BP resonant/Lyapunov/DPO/LPO family confirmation gate

**Task:** `#776`, the second concrete per-system task of `#760`'s new-system discovery campaign,
spec-complete in `docs/notes/2026-08-01-771-neptune-triton-scoping.md` §6 (recommended by
`#771`'s own scoping pass). Family confirmation ONLY, mirroring `#765`'s own Task-A-only scope
-- no connection/heteroclinic/chain work, no retrograde families, no catalogue writeback.

**Source data** (acquired, filed, digested this task; see
`docs/notes/2026-08-01-776-miceli-bosanac-2026-neptune-triton-digest.md` for the full digest and
citation-mining pass): Miceli, G.E. & Bosanac, N. (2026), "Generating Planar Trajectories for
Neptunian System Exploration Using Motion Primitives," *J. Astronaut. Sci.* 73:11, DOI
`10.1007/s40295-025-00545-z` (open access, CC-BY 4.0), plus its own four Springer Electronic
Supplementary Material files (Online Resources 1-4; ESM2-4 are machine-readable text files of
12-decimal nondimensional CR3BP periodic-orbit states/periods/Jacobi constants) and the companion
Miceli 2025 CU Boulder PhD dissertation (Table 2.1's equilibrium-point positive control). All six
files re-downloaded and md5-verified this task to match `#771`'s own scoping-note record exactly
(JAS-2026 PDF `576ce77e832861c1011efd08749b654e`; dissertation PDF
`c3ed8d2ae0824f1e8bdd59420128f7ff`; ESM2/ESM3/ESM4 line counts 93/297/752). Filed at
`cyclers_pdf/papers/` (see `docs/notes/CORPUS_INDEX.md`'s own `#776` entries).

**Code delivered:** `src/cyclerfinder/search/neptune_triton_resonant_families.py` (new module, a
thin sibling of `saturn_titan_resonant_families.py`, reusing its corrector/classification/
continuation machinery directly -- no reimplementation) + `tests/search/test_neptune_triton_resonant_families.py`
(47 tests, all passing). Both pass `ruff check`, `ruff format --check`, and the project's
canonical `uv run mypy src tests` (829 files, clean).

---

## Unlike Vaquero 2013: this source's data is already nondimensional at 12 decimals

Vaquero's own Table 4.1 (the Saturn-Titan module's source) printed dimensional `(x [km], ydot
[km/s], T [days])` at 6 significant figures, requiring this project's own registry-derived `l*`/
`t*` and a self-validation step. Miceli & Bosanac's own ESM files print `(x0, y0, z0, xdot0,
ydot0, zdot0, Integration Time, Jacobi Constant)` already nondimensional to 12 decimals, under a
mass ratio stated to 16 digits in every file header. There is therefore NO `l*`/`t*` sensitivity
anywhere in this module's primary gate -- `MICELI_L_KM`/`MICELI_T_S` (the paper's own stated
`l* = 354,760 km` / `t* ~= 8.081353e4 s`, independently re-grepped from the JAS-2026 paper's own
body text, not just the dissertation) enter ONLY for human-readable period-in-days reporting.

## The ten vendored gate rows -- a genuinely clean 10/10 sweep

| Label | Source | Periodicity resid. | x0/ydot0/T rel. err | Barden vs. planar_floquet | is_real_unstable | max &#124;λ&#124; | Passed |
|---|---|---|---|---|---|---|---|
| 1:7 | ESM3 `Res17+x+h` | 9.0e-11 | 5.3e-13 / 5.1e-13 / 1.8e-13 | 2.9e-07 | False | 1.0 | **PASS** |
| 3:2-start | ESM4 `Res32-x+h` | 1.1e-11 | 0 / 2.4e-12 / 7.8e-13 | 6.1e-07 | False | 1.0 | **PASS** |
| L2-lyapunov-target | ESM4 `L2LyapunovTarget` | 1.9e-11 | 0 / 2.4e-12 / 8.9e-12 | 8.0e-10 | True | 332.6 | **PASS** |
| L1-lyapunov | ESM2 `L1Lyapunov` | 1.7e-09 | 6.6e-13 / 5.8e-11 / 7.1e-13 | 8.4e-11 | True | 2024.7 | **PASS** |
| L2-lyapunov | ESM2 `L2Lyapunov` | 8.3e-10 | 4.7e-13 / 3.9e-11 / 5.3e-13 | 2.7e-12 | True | 1741.7 | **PASS** |
| DPO | ESM4 `DPO` | 4.8e-11 | 0 / 1.9e-11 / 1.7e-11 | 3.2e-11 | True | 3.03 | **PASS** |
| LPO | ESM4 `LPO` | 1.2e-10 | 0 / 3.7e-11 / 2.0e-11 | 5.2e-14 | False | 1.0 | **PASS** |
| 4:5-saddle | ESM4 `Res45+x+h` | 2.4e-09 | 5.7e-13 / 3.0e-12 / 3.2e-13 | 4.9e-12 | True | -105.05 | **PASS** |
| 4:7-stress | ESM4 `Res47-x+h` | 9.8e-08 | 0 / 3.7e-13 / 8.1e-14 | 6.5e-08 | True | 14624.1 | **PASS** |
| 4:3-saddle | ESM4 `Res43+x+h` | 1.2e-10 | 5.0e-13 / 1.3e-12 / 6.7e-14 | 1.7e-10 | True | 20.4 | **PASS** |

**All ten rows pass every one of the three gate criteria**: (a) periodicity self-consistency
(`PERIODICITY_GATE_TOL = 1e-6` nondim; worst observed 9.8e-8, the strongly-unstable 4:7-stress
row, still more than an order of magnitude inside), (b) reproduction (`REPRODUCTION_GATE_REL_TOL
= 1e-6`; worst observed 5.8e-11, four-plus orders of magnitude inside), (c) internal
Barden-vs-`_planar_floquet` eigenvalue cross-check (`CROSSCHECK_GATE_REL_TOL = 1e-5`; worst
observed 6.5e-7, more than an order of magnitude inside). This is a genuinely clean sweep, unlike
the Saturn-Titan module's own two honest partial fails (6:5's eigenvalue near-miss, L2's
period-transcription anomaly) -- explained entirely by this data's own already-nondimensional
12-decimal precision, with no unit-conversion/rounding loss like Vaquero's dimensional km/s
table. Per the dispatch note's own framing, a partial result would have been an equally
fine outcome; this is simply what the evidence shows.

**Precision caveat, stated honestly**: eigenvalues are labelled `is_real_unstable`/`max_eigenvalue`
and cross-checked ONLY internally (Barden vs. an independent full-period monodromy
eigendecomposition) -- there is NO published eigenvalue table for this system (unlike Vaquero's
own Table 4.1), so nothing here is called a "reproduction" on that axis, only a self-consistency
check, per the dispatch note's own explicit instruction to be precise about that distinction.

Three rows (1:7, 3:2-start, LPO) are genuinely near-unit-circle (`is_real_unstable=False`,
`max_eigenvalue=1.0` to <1e-6) -- NOT saddles. This is irrelevant to `passed` (the gate never
gates on saddleness, per `#771`'s own scoping note) and matches `#771`'s own eigenvalue survey,
which characterized these labels' own general ranges without claiming every row is unstable.

## Textual cross-check: the "1:7" row IS the paper's own named target orbit

The JAS-2026 paper's own body text (Sec. 4.1) states "a 1:7 resonant orbit with C_J = 1.8 and a
period of 41.14 days." ESM3's own `Res17+x+h` row states `C_J = 1.806962818639` (rounds to 1.8)
and `Integration Time = 43.981049667607` nondim, which converts via this module's own `l*`/`t*`
to **41.137 days** -- matching "41.14 days" to <0.01% relative. This is independent,
paper-text-level corroboration that the vendored row is genuinely the paper's own named Scenario-1
target orbit, not merely a same-labeled coincidence.

## Equilibrium-point positive control (dissertation Table 2.1)

| Point | Target x (Table 2.1) | Recovered x | Rel. err |
|---|---|---|---|
| L1 | 0.959217 | 0.9592169511... | 5.1e-08 |
| L2 | 1.041493 | 1.0414934478... | 4.3e-07 |
| L3 | -1.000087 | -1.0000870626... | 6.3e-08 |
| L4 | 0.499791 | 0.4997910497... | 9.9e-08 |
| L5 | 0.499791 | 0.4997910497... | 9.9e-08 |

All five reproduce to <5e-7 relative, via `lagrange_collinear_x` (L1/L2, reused directly from
`search/reachable_representatives.py`), a new small L3 root-find (the same `dUbar/dx = 0`
equation, `search/cr3bp_periodic._ubar_grad_x_at_axis`, a different bracket), and the closed-form
`x = 0.5 - mu` for L4/L5.

## Item (d): `two_body_resonant_seed` lineage + continuation onto two multi-member families

**Two-body seed lineage -- honest, well-characterized negative.** `two_body_resonant_seed(4, 3,
x0_sign=-1)` and `two_body_resonant_seed(2, 3, x0_sign=-1)`, each converged directly at its own
natural Jacobi constant, both converge CLEANLY (residuals ~1e-14/1e-15) but to the WRONG
resonance topology -- the "4:3" seed lands on a `period/2pi ~= 4.0` orbit, and the "2:3" seed
lands on a `period/2pi ~= 7.0` orbit, neither matching its own `p:q` label. This is a third,
independent confirmation of the same qualitative finding Anderson & Lo (Jupiter-Europa) and
Vaquero (Saturn-Titan) each document for the analogous naive attempt in their own systems.
`x0_sign=+1` was NOT attempted as part of any gate: it places `x0` exactly at the secondary's own
orbital radius, the DOP853 step-collapse hazard already documented in
`jovian_resonant_families.survey_candidates`'s own docstring (both `(4,3,+1)` and `(2,3,+1)`
confirmed to hang past a 20s wall-clock budget in scratchpad testing this task -- not run as part
of any committed code path).

**Continuation onto two multi-member families -- two clean confirmations, one honest negative.**
The ESM4 "Res23-x+h" family (3 printed members) and "Res43+x+h" family (4 printed members) are
the multi-member families `#771`'s own dispatch note calls out by name.

* **2:3 family** (`continue_23_family`): seeded from the lowest-C printed member, walking in
  `d_jacobi=2e-4` steps toward the highest-C member. Reaches `StopReason.JACOBI_BOUND` cleanly --
  229 gauntlet-passing members, 0 rejections -- landing on the 2nd/3rd printed members to
  6.0e-5/5.7e-7 and 5.6e-4/3.6e-5 relative (x0/period) respectively. A genuine, tight,
  independent seed-lineage confirmation of the whole family.
* **4:3 family, saddle half** (`continue_43_saddle_family`): the printed "Res43+x+h" rows split
  across TWO topologically distinct branches at nearby C -- an honest finding, not a data error
  (see below). The genuine-saddle half (`half_crossings=2`, `|lambda| ~= 16-20`) continues
  CLEANLY: `StopReason.JACOBI_BOUND`, 51 gauntlet-passing members, 0 rejections, landing on the
  printed endpoint to 1.0e-4 relative x0 / 9.4e-6 relative period. A second, independent
  multi-member family confirmation.
* **4:3 family, near-unit half** (`continue_43_near_unit_family_fold_reversal`, documented
  negative): the OTHER half of the printed "Res43+x+h" rows (`half_crossings=1`, near-unit-circle,
  NOT saddles) does NOT cleanly continue to its own higher-C printed member -- the natural-
  parameter walk hits `StopReason.FOLD_REVERSAL` after 95 members, stopping at C=3.0316, well
  short of the target C=3.0408. This is reported honestly as a genuine, well-characterized
  negative: the "4:3" label spans (at least) two topologically distinct branches at nearby C, not
  one smooth curve -- exactly the kind of finding this project's own standing practice treats as
  valuable, not a gate failure to paper over. (The genuine-saddle half above, which IS one smooth
  branch, gives the required "at least two multi-member families" confirmation on its own,
  together with the clean 2:3 family.)

---

## Verification

* `uv run ruff check` / `ruff format --check` on both new files: clean.
* `uv run mypy src tests` (project canonical invocation): clean, 829 files.
* `uv run pytest tests/search/test_neptune_triton_resonant_families.py -q`: 47/47 pass.
* `uv run pytest tests/data/test_outstanding_structure.py tests/data/test_outstanding_header_body_consistency.py -q`:
  run before committing the `OUTSTANDING.md` update (see commit history for pass/fail status
  recorded at commit time).

---

## Explicitly out of scope (per the dispatch note)

No manifold/homoclinic/heteroclinic/chain work of any kind, no retrograde families (no
digit-grade data for them -- `#771` Sec. 1), no catalogue writeback (nothing here is
catalogue-eligible until confirmed AND novelty-checked via `search/literature_check.py`, a
separate future task), and no vendoring of the remaining ~40 dataset rows (this task vendored 10
primary gate rows + 6 more continuation-family-only rows, ~16 total, deliberately narrow per the
dispatch note's own "prove narrow first" instruction).

## Registered follow-up (not dispatched)

`#777` -- registered, NOT dispatched: vendor the remaining ~34 rows of the Miceli & Bosanac
2026 ESM dataset not already covered by `#776` (the full ~50-row dataset spans ~20 resonant
labels 1:2...4:7 plus additional DPO/manifold-arc primitives this task's own narrow 10+6-row
subset does not touch), extending the family-confirmation gate to near-complete dataset coverage.
An explicit, deliberate follow-up per the user's own "option 1" choice (prove the narrow ~10-row
version first) -- not run this task.

## Recommendation for a Task-B analog -- opinion, not a decision

Neither `#776`'s own dispatch note nor this task attempted any connection-stage work. Unlike the
Saturn-Titan thread (`#768`/`#773`/`#775`), which found chain-closure work genuinely hard even
with a good IC source, this system's own literature contains NO published homoclinic/heteroclinic
connection states at all (`#771` Sec. 2) -- so a future connection-stage task here would be
genuinely NOVEL-territory work (subject to `search/literature_check.py`'s own novelty gate),
distinct from the Saturn-Titan case's reproduction-shaped Fig. 4.9-4.12 target. The 4:5-saddle and
4:7-stress rows (both genuine, tightly-confirmed real saddles, `|lambda| ~= 105` and `~= 1.46e4`
respectively) are the most promising manifold-source candidates for such a task, per `#771`'s own
eigenvalue survey. This is my assessment for the user to weigh, not a decision -- a
connection-stage task is a separate dispatch, as is `#777` above.
