# Digest: Miceli & Bosanac (2026), "Generating Planar Trajectories for Neptunian System Exploration Using Motion Primitives"

**Task:** `#776`, acquisition + digest step (per `[[feedback_corpus_document_policy]]` /
`[[feedback_per_paper_digest_todo]]`) preceding the actual gate-reproduction work — see
`docs/notes/2026-08-01-776-neptune-triton-resonant-families-gate.md` for the code/gate results.
Full citation, filing paths, and md5s are in `docs/notes/CORPUS_INDEX.md`'s own entry; not
repeated here. This continues the `#760` new-system discovery campaign
(Jupiter-Europa `#752`→`#759`, Saturn-Titan `#764`→`#765`→`#767`→`#773`→`#775`, Neptune-Triton
scoped by `#771` and dispatched as `#776`).

**Scope of this digest:** full-page read of the JAS-2026 paper itself (a paper, not a thesis —
full corpus-policy scope, not chapter-summary) plus its own Online Resources 2-4 (the
machine-readable ESM data files, the load-bearing anchor for this task's own gate); the
companion Miceli 2025 PhD dissertation is treated at chapter-summary scope (per corpus policy
§2's thesis rule), reading only Ch.2 (CR3BP/characteristic-quantity background, cross-checked
against this project's own `core/cr3bp.py`) and Table 2.1 (the equilibrium-point positive
control), plus a light TOC/reference-list spot-check.

## The paper's own content

Builds on Smith & Bosanac's own "motion primitive" trajectory-design method (short, pre-computed
arcs — periodic-orbit segments, manifold arcs — chained via graph search into a full transfer)
and applies it, for the first time, to the Neptune-Triton system (previously only cislunar/
Jupiter-Europa applications existed in this lineage). Two worked design scenarios: Scenario 1
(Sec. 4.1, high-energy arrival from an interplanetary trajectory to a 1:7 resonant target orbit,
`C_J = 1.8`, period 41.14 days) and Scenario 2 (Sec. 4.2, low-energy transfer from a 3:2 resonant
orbit to a Triton-region Low Prograde Orbit). Sec. 3 walks through the method itself using a
smaller L1/L2 Lyapunov example. Table 2.1-equivalent equilibrium-point content (this paper omits
the numeric table itself, deferring to the dissertation) and the full CR3BP formulation are
standard (`mu = 2U* - v^2` Jacobi convention, this project's own convention exactly — independently
verified this task, see below).

**What makes this paper different from Vaquero 2013 (the Saturn-Titan module's own source):** it
ships FOUR supplementary data files as formal Springer Electronic Supplementary Material —
genuinely a first for this project's own literature (Vaquero's own Table 4.1 was a printed,
6-significant-figure dimensional table; this paper's own ESM files are a machine-readable,
12-decimal, ALREADY-NONDIMENSIONAL text dump). This is independently confirmed a materially
stronger IC/period/Jacobi-constant anchor than any prior system this project has gated against —
at the cost of publishing no per-orbit eigenvalue table (unlike Vaquero's own Table 4.1), so the
stability-classification axis of this task's own gate is internally-cross-checked (Barden vs
`_planar_floquet`), never "reproduced" against a published target.

## Independent re-verification this task (not inherited from `#771`'s own scoping note)

* **md5s**: re-downloaded all six files from the exact URLs `#771`'s own note records (Springer
  `link.springer.com` for the main PDF, `static-content.springer.com` ESM URLs for the four
  supplementary files, `colorado.edu/faculty/bosanac` for the dissertation). Every md5 matches
  `#771`'s own scoping-note md5 EXACTLY: JAS-2026 PDF `576ce77e832861c1011efd08749b654e` (40 pp);
  dissertation PDF `c3ed8d2ae0824f1e8bdd59420128f7ff`. Line counts of ESM2/ESM3/ESM4 (93/297/752)
  also match exactly.
* **Mass ratio**: every one of the three ESM text files' own header states, verbatim,
  `Mass ratio: 2.089503183689124e-04`; both the JAS-2026 paper's own body text (p.5) and the
  dissertation's own body text independently state the rounded display `mu ~= 0.00020895`,
  consistent.
* **Characteristic quantities**: JAS-2026 p.5 body text (NOT just the dissertation, a stronger
  source than `#771`'s own note implied) states, verbatim, `l* = 354, 760 km` and
  `t* ~= 8.081353 x 10^4 s` — used in `search/neptune_triton_resonant_families.py` for
  period-in-days reporting only.
* **Table 2.1 (dissertation, p.~59)**: re-grepped directly, verbatim: L1 `x=0.959217`,
  L2 `x=1.041493`, L3 `x=-1.000087`, L4/L5 `x=0.499791` (y omitted from the gate, +-sqrt(3)/2).
  This module's own `lagrange_collinear_x` (L1/L2, reused directly from
  `search/reachable_representatives.py`) plus a new tiny L3 root-find (same `dUbar/dx=0`
  equation, `search/cr3bp_periodic._ubar_grad_x_at_axis`, different bracket) and the closed-form
  `x = 0.5 - mu` for L4/L5 reproduce all five to <5e-7 relative — see the gate module.
* **1:7 target orbit cross-check**: the paper's own body text (Sec. 4.1) states "a 1:7 resonant
  orbit with C_J = 1.8 and a period of 41.14 days" — the ESM3 file's own `Res17+x+h` row states
  `C_J = 1.806962818639` (rounds to 1.8) and `Integration Time = 43.981049667607` nondim, which
  converts via this module's own `l*`/`t*` to 41.137 days (matches "41.14 days" to <0.01%
  relative) — an independent, paper-text-level corroboration that ESM3's own `Res17+x+h` row IS
  the paper's own named Scenario-1 target orbit (not merely a same-labeled coincidence).
* **µ registry delta**: this project's own DE440-registry Neptune-Triton mu
  (`core.satellites` registry, `GM_Triton / (GM_Neptune-system + GM_Triton) = 1428.49546 /
  6837955.59604 = 2.08907e-4`) differs from the paper's own value by ~2.08e-4 relative — the same
  GM-vintage-delta class the Jovian/Saturn-Titan modules each document for their own systems;
  this module uses the paper's own value, per the standing per-system pattern.

## Citation-mining pass (full reference list, pp.42-43, 62 entries)

Per `[[feedback_corpus_document_policy]]`'s mandatory step 3 — full findings and the resulting
backlog registration are in `docs/notes/2026-07-27-730-acquisition-backlog-master-list.md`'s
2026-08-01 update (not duplicated here). Summary: reference [36] is Vaquero 2013 — already
corpus (`#765`), confirming (as `#771`'s own note already flagged) the whole resonant-orbit
lineage traces back to the same method family. Reference [62] is the AIAA-2024 precursor,
already confirmed paywalled/superseded by `#771`. Reference [27] (Gillespie, Miceli & Bosanac
2025 AAS, the cislunar sibling application) was already flagged not-found-open by `#771`. TWO
genuinely new candidates surfaced and registered (not acquired, out of `#776`'s own narrow
scope): Restrepo & Russell 2018 (a multi-system planar periodic-orbit database, DOI
`10.1007/s10569-018-9844-6`, medium priority) and Marley 2010 (a NASA JPL Neptune-Triton-KBO
mission-architecture technical report, no DOI, low priority/background only). Everything else in
the 62-entry list is either general-methods background (graph search, clustering/HDBSCAN,
reinforcement learning trajectory design — a different methodological branch this project does
not use) or already-corpus astrodynamics foundations (Szebehely 1967, Park/Folkner/Williams/
Boggs DE440/441, Murray & Dermott).

## Verdict

**ADOPT (already adopted this task).** The three ESM data files are the primary validation
anchor for `#776`'s own gate module (`search/neptune_triton_resonant_families.py`) — not
background reading; the dissertation contributes the Table 2.1 equilibrium-point positive
control. All ten vendored gate rows pass every one of the three gate criteria (periodicity
self-consistency, reproduction, internal Barden-vs-planar_floquet cross-check) — see the
results note for the full row-by-row report. The paper's own graph-search/motion-primitive
method itself (Secs. 3-4, the transfer-DESIGN machinery) is NOT adopted anywhere in
`cyclerfinder` (a structurally different capability — chaining pre-computed arcs via A*/k-shortest-
path graph search — from this project's own resonant-family/manifold/connection pipeline); noted
for completeness, not a gap, since `#776`'s own scope is family confirmation only.
