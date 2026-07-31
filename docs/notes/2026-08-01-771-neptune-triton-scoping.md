# Scoping `#771`: Neptune-Triton deeper pass — a genuine digit-grade anchor EXISTS, and `#764`'s "no IC table" finding is OVERTURNED

**Task:** `#771` (research/scoping only — no repo code, no catalogue changes), the second
per-system candidate pass of `#760`'s new-system discovery campaign, mirroring `#764`'s own
scoping format. `#764`'s first look at Neptune-Triton (from the AIAA-2024 full text) found real
recent literature but "no IC, period, or eigenvalue tables ... only one 6-digit `C_J`" and
ranked the system second, behind Saturn-Titan, with a "weak gate" caveat. This pass re-verified
that finding independently instead of inheriting it — and it no longer holds.

**Headline: the journal version of the same work (Miceli & Bosanac, *J. Astronaut. Sci.* 73:11,
2026, open access, CC-BY 4.0) ships FOUR supplementary files (Online Resources 1-4), three of
which are machine-readable text files containing 12-decimal nondimensional initial states,
integration times (verified this pass to be full orbital periods), and Jacobi constants for
~50 planar Neptune-Triton CR3BP periodic orbits across ~20 resonant-orbit labels plus
L1/L2 Lyapunov, DPO, and LPO families — under a mass ratio printed to 16 digits
(`µ = 2.089503183689124e-04`).** This is a BETTER-precision IC anchor than Vaquero 2013's
Table 4.1 (6 significant figures, `µ ≈` 5 sig figs) on every axis except one: no per-orbit
eigenvalues are published. Verdict below: **GO**, with a `#765`-shaped family-confirmation-only
first task (`#776`).

---

## 1. Sources acquired and read first-hand this pass

All downloads live in this session's scratchpad; **nothing is filed into `cyclers_pdf/` yet** —
filing + `CORPUS_INDEX.md` registration + per-paper digest todos are Step 0 of the recommended
first task, exactly as `#764` handled Vaquero (downloaded/verified at scoping, filed by `#765`).

* **Miceli, G.E. & Bosanac, N., "Generating Planar Trajectories for Neptunian System
  Exploration Using Motion Primitives," *The Journal of the Astronautical Sciences* 73:11
  (2026), DOI `10.1007/s40295-025-00545-z`** — open access (CC-BY 4.0), downloaded from
  Springer; md5 of the copy read this pass: `576ce77e832861c1011efd08749b654e` (40 pp).
  Read via full text extraction; every load-bearing claim below grepped/read from the text
  layer directly.
* **Online Resources 1-4** (Springer ESM, downloaded this pass from
  `https://static-content.springer.com/esm/art%3A10.1007%2Fs40295-025-00545-z/MediaObjects/40295_2025_545_MOESM{1..4}_ESM.{pdf,txt}`):
  * **ESM1** (PDF, 2.4 MB): flowcharts + governing-parameter tables for the motion-primitive
    method itself — no orbit data.
  * **ESM2** (txt, 93 lines): the Sec. 3 walkthrough example — 1 L1 Lyapunov + 1 L2 Lyapunov
    orbit + 75 manifold-arc primitives at `C_J = 3.013995763057`.
  * **ESM3** (txt, 297 lines): Scenario 1 (high-energy arrival → **1:7 resonant target orbit at
    `C_J = 1.806962818639`**, the JAS version's target — the AIAA-2024 version's 3:4 target at
    `C_J = 1.75598` was superseded) — 21 resonant periodic-orbit rows spanning 1:2 … 4:7 at
    `C_J ≈ 1.56–2.46`, + 245 4:5-resonant-manifold arcs + the sourced interplanetary-arrival
    state (`C_J = 0.963140720704`, matching `#764`'s quoted `0.963141`).
  * **ESM4** (txt, 752 lines): Scenario 2 (low-energy 3:2 resonant → Triton LPO) — 32 resonant
    rows at `C_J ≈ 2.93–3.67`, 8 DPOs, 1 LPO, 1 L2 Lyapunov target, + 673 manifold arcs.
  * File-header format: `Primitive  x0 y0 z0 xdot0 ydot0 zdot0 IntegrationTime JacobiConstant`,
    12 decimals, with `System: Neptune-Triton CR3BP` and `Mass ratio: 2.089503183689124e-04`
    stated in each header.
* **Miceli, G.E. (2025), "Data-Driven Spacecraft Trajectory Design for Planetary System
  Exploration," PhD dissertation, University of Colorado Boulder (advisor N. Bosanac)** —
  freely hosted on the Bosanac group site
  (`https://www.colorado.edu/faculty/bosanac/sites/default/files/2025-08/Miceli_PhDDissertation.pdf`);
  md5 of the copy read this pass: `c3ed8d2ae0824f1e8bdd59420128f7ff`. Adds a THIRD
  (medium-energy, 4:3-resonant-start) design scenario the JAS paper drops, and Table 2.1's
  equilibrium-point locations (`L1 x=0.959217`, `L2 x=1.041493`, `L3 x=-1.000087`,
  `L4/L5 x=0.499791`) — a cheap positive-control gate row. No per-orbit eigenvalue tables
  found in the dissertation either (checked its full tables list).
* **AIAA-2024 (`10.2514/6.2024-1280`), the version `#764` originally grepped:** paywalled at
  ARC (HTTP 403 this pass, no open copy found on the group site). NOT required — the JAS
  version + ESM data supersede it. One content delta worth recording: AIAA-2024 described 3:1
  and 4:1 **retrograde** resonant families; the JAS ESM files contain **prograde (`+h`) rows
  only**, so retrograde families have no digit-grade anchor.
* Related-work chain for the eventual digest pass: AAS 24-161, AAS 25-138 (conference
  precursors, not found open-access); Gillespie/Miceli/Bosanac AAS 2025 (the primitive-
  generation method, cislunar); Smith & Bosanac (Earth-Moon primitives). Notably, the paper's
  own resonant-orbit reference [36] is **literally Vaquero 2013** — the same thesis already in
  corpus (`#765`), confirming the whole lineage is one method family.

## 2. First-hand verification of the supplementary data (the load-bearing step)

`#764`'s Neptune-Triton verdict rested on the printed paper alone. This pass verified the ESM
data itself, from the raw files, with a throwaway scratchpad script (standard planar CR3BP,
DOP853 at `rtol=atol=1e-12`, using the ESM header's own 16-digit `µ` — no repo code touched):

* **Jacobi-constant convention check** (6 rows spanning all three data files, resonant +
  Lyapunov + DPO): recomputing `C_J = 2U* − v²` from each row's own state reproduces the row's
  stated Jacobi constant to **`1.3e-14`–`3.1e-12`**. Their convention is exactly this repo's.
* **Periodicity check** (same 6 rows): propagating each row's IC for its stated Integration
  Time returns to the initial state to **`1.4e-12`–`1.9e-9`** (nondim). So "Integration Time"
  for periodic-orbit rows IS the full period, and these are genuine periodic orbits under this
  repo's own dynamics — a real cross-implementation reproduction, not a formatting guess. (The
  paper states their multiple-shooting corrector converges to `1e-10`; consistent.)
* **Eigenvalue survey** (all 53 resonant + LPO periodic-orbit rows, 4x4 planar STM): the
  resonant rows are a mix of near-unit-circle members (consistent with AIAA-2024's "in-plane
  stability indices ... close to 2" language `#764` quoted) and **genuine real saddles**:
  `1:2−x…1:6−x` at `|λ| ≈ 1.2–1.6`, 4:3 members at `|λ| ≈ 16–25`, 2:5/3:5/3:7 members at
  `|λ| ≈ 58–80`, the Scenario-1 4:5 manifold-source family at `|λ| ≈ 105`, and one 4:7 row at
  `|λ| ≈ 1.5e4`. The Lyapunov/DPO rows are strong saddles (`|λ| ≈ 333–3667`). Every magnitude
  is inside or below the range the existing Newton machinery has already handled well
  (`#753`/`#765`: `~1e3–2e3`; trouble historically began near `~4.4e3` only for *chain* work,
  `#759`/`#768`) — family-stage tractability risk: LOW.
* Naming semantics (from the paper, Sec. 3.1.1): `p:q` = p spacecraft revs about Neptune per
  ~q Triton revs (same shape as this project's Anderson-Lo convention; checked: `Res17`'s
  period `43.981 ≈ 7·2π`); `±x` = periapse/apoapse of the defining state on the `+x̂`/`−x̂`
  side of Neptune; `±h` = orbital angular momentum along `±ẑ` (prograde/retrograde).

**Net: the "no digit-grade gate" objection is gone.** What remains genuinely unpublished:
per-orbit eigenvalues/stability indices (qualitative prose only), any homoclinic/heteroclinic
connection states (they build manifold-arc primitives and graph-search transfers, never
connection orbits — connection-stage work in this system would be NOVEL, subject to
`search/literature_check.py`), and retrograde-family data.

## 3. Catalogue and corpus state (re-checked this pass)

* `data/catalogue.yaml`: zero Neptune-Triton rows; the only "Neptune" hits are Voyager 2's
  Grand-Tour `mga_tour` row. No seed/µ duplication risk. (Confirms `#764`.)
* `cyclers_pdf/papers/` + `CORPUS_INDEX.md`: no Miceli/Bosanac/Triton entry; the only Neptune
  documents are the Voyager-2 Science encounter papers (`#429`). Both papers above are new to
  the corpus.

## 4. Ranking check: does the deeper look change `#764`'s ordering?

Only in Neptune-Triton's favor. Saturn-Enceladus and Pluto-Charon gained nothing this pass
(nothing new searched for or found on them — `#772` stays parked on `#764`'s reasoning).
Neptune-Triton now holds a digit-grade anchor **stronger on the IC axis than the one that made
Saturn-Titan tractable** (12-decimal machine-readable states + exact 16-digit `µ`, vs Table
4.1's 6 significant figures + "`µ ≈`" 5 sig figs — the exact imprecision behind `#765`'s
honest 6:5 eigenvalue near-miss and `#769`). Its one weaker axis: no published eigenvalue
targets, so the eigenvalue side of a gate is internal-cross-check-corroborated rather than
source-gated. The `#764` retrograde-frame subtlety stays resolved (isolated planar CR3BP is
orientation-agnostic; the fidelity caveat about Neptune J2/inclined real-ephemeris contexts is
a writeback-tier caveat, not a pipeline blocker).

## 5. Pipeline repurposability (per `#764` §5, re-confirmed against this system)

Nothing about Neptune-Triton needs structural changes. `µ = 2.0895e-4` sits between
Jupiter-Europa (`2.528e-5`) and Saturn-Titan (`2.3658e-4`), both proven. The generic stack
(`two_body_resonant_seed`, `converge_candidate`/`survey_candidates`,
`correct_symmetric_fixed_jacobi`, Barden + `_planar_floquet`, `basin_robustness_scan`,
`ResonantNode`/`correct_connection`/`find_homoclinic`) is system-agnostic — `#765`/`#767`
already demonstrated the thin-sibling pattern end to end. A Neptune-Triton module is
constants + sourced rows + gate, exactly like `saturn_titan_resonant_families.py`. Two
system-specific notes for the implementing task:
* The ESM data is **already nondimensional**, so unlike `#765` there is NO `l*`/`t*`
  sensitivity anywhere in the primary gate — dimensional quantities (km, days; the paper's own
  `l* = 354,760 km`, `m* ≈ 1.02457e26 kg`) enter only for human-readable reporting. The
  registry-vs-paper `µ` delta (registry `1428.49546/6.836527100580e6 = 2.08951e-4` per `#764`
  vs the paper's exact `2.089503183689124e-04`) should be measured and documented, and the
  SOURCE's value used, per the standing per-system pattern.
* Many sourced rows are NOT saddles (near-unit eigenvalues). `ResonantNode.from_candidate`
  correctly refuses those — irrelevant for the family-confirmation task (its gate is
  IC/period/`C_J`, not saddleness), but the later connection-stage task must pick its nodes
  from the saddle subset (Sec. 2's survey identifies them).

## 6. Recommendation: **GO.** Concrete first task (`#776`, dispatchable as written)

**Task `#776` — "Neptune-Triton planar CR3BP periodic-orbit families vs the Miceli-Bosanac
JAS-2026 supplementary-data gate"** (family confirmation ONLY, mirroring `#765`'s Task-A-only
scope; Sonnet-tier per `[[feedback_subagent_model_tiering]]`):

0. **Acquire + file the anchors**: the JAS-2026 PDF (CC-BY 4.0; md5 above), all four ESM files
   (the data files ARE the anchor — file them alongside the PDF in `cyclers_pdf/papers/`), and
   the Miceli 2025 dissertation (md5 above) → `CORPUS_INDEX.md` + per-paper digest todos per
   `[[feedback_per_paper_digest_todo]]`. Record AIAA-2024 (`10.2514/6.2024-1280`) as the
   paywalled, superseded conference precursor and AAS 24-161/25-138 as not-found-open. Re-verify
   the md5s and every gate number below against the files directly, per
   `[[feedback_ground_citations_against_content]]` — this note is one hop removed.
1. **New module `search/neptune_triton_resonant_families.py`** (thin sibling of
   `saturn_titan_resonant_families.py`): `neptune_triton_system()` at the SOURCE's own
   `µ = 2.089503183689124e-04` (ESM headers), `l* = 354,760 km` / paper `m*` for
   dimensional reporting only; measure + document the registry-vintage `µ` delta. Vendor the
   selected gate rows as sourced constants (provenance: file + row label + `C_J`).
2. **Gate rows** (~10, spanning all three data files and both object classes; each verbatim
   from the ESM files): the three scenario-defining orbits — 1:7 target
   (`C_J = 1.806962818639`), 3:2 start (`C_J = 3.028835529717`), L2 Lyapunov target
   (`C_J = 3.003706759619`) — plus ESM2's L1/L2 Lyapunov pair, one DPO, the Triton LPO, and
   three saddle-class resonant members incl. the Scenario-1 4:5 manifold source
   (`C_J = 2.987089791658`, `|λ|≈105`) and the strongest-instability 4:7 row
   (`C_J = 2.997230642137`, `|λ|≈1.5e4`, the deliberate stress row).
3. **Gate (dual-criterion, honest per-row report, `GateRow`-style)**: (a) data
   self-consistency — propagate each row's IC for its stated period, gate the return distance
   (scoping pass saw `1e-9`–`1e-12`; tolerance to be justified against the paper's own stated
   `1e-10` corrector tolerance, not assumed); (b) reproduction — re-converge each row with OUR
   corrector at the row's own `C_J` and gate `x0`/`ẏ0`/period against the printed 12-decimal
   values; (c) independent Barden vs `_planar_floquet` cross-check on every row (standing
   discipline — computed eigenvalues have NO published target here and must be labelled as
   internally-cross-checked, never "reproduced"); (d) seed-lineage check on at least two
   resonant families: this project's own `two_body_resonant_seed(p,q)` → converge → continue
   in `C_J` to the row's own value → confirm it lands on the printed member (the paper builds
   its seeds with literally this construction, so this closes the loop family-to-family; the
   multi-member families in ESM4 — e.g. four 4:3 rows, three 2:3 rows — give continuation
   targets at several `C_J` values each). Positive control: L1/L2/L3 locations vs dissertation
   Table 2.1 (6-decimal x-values above).
4. **Explicitly out of scope for `#776`**: manifolds, homoclinic/heteroclinic connections and
   chains (novel-territory work in this system — a separate later task, and per the
   `#768`/`#773`/`#775` lesson, chain-closure work is genuinely hard even with a good source
   and must never be bundled into the family task); retrograde families (no digit-grade data);
   any catalogue writeback (nothing here is catalogue-eligible until confirmed AND
   novelty-checked via `search/literature_check.py`).

**Honest counterweight:** `#776` is reproduction-shaped by design, and its anchor is
supplementary *data* rather than a peer-typeset table — mitigated by the CC-BY journal of
record explicitly describing these files ("initial state, integration time, and Jacobi
constant of all primitives," Data Availability statement) and by this pass's own independent
`1e-12`-level convention/periodicity verification. The eigenvalue axis stays source-ungated
(unlike Vaquero) — the gate's discriminating power is the 12-decimal IC/period/`C_J`
reproduction plus family-lineage continuation, which is materially stronger than the
"structural-grade only" assessment `#764` gave this system. Discovery upside after the gate:
connection/chain states in a system with NO published connection literature found (novelty
gate should bite less than Saturn-Titan's), and first-ever Neptune-Triton catalogue rows.
