# Discovery-capability paper sweep — triage summary (2026-06-13)

**Task:** #250. A literature sweep + digest hunting papers that bring NEW
discovery capability — a novel cycler-search **genome**, a discovery **method**,
or a new cycler **family/topology** our current genome can't represent — NOT
just another cycler to catalogue as a validation anchor. Framing per
`docs/notes/2026-06-13-discovery-program-spec.md`.
**Discipline:** cite by publication only; NO catalogue writeback (digests are
notes + plans). Read-before-judging — no title-dismissals.

## Network limitation hit (material)

`WebFetch` was **DENIED** in this environment for the entire run; `WebSearch`
worked normally. Consequence: I could read every candidate's **title, authors,
venue, abstract, and search-engine summaries**, but could **not transcribe PDF
tables / equations / initial-condition tuples** for any paper. Every digest is
abstract-and-metadata grade; the one ADOPT-track item (Roberts-Tsoukkas 2026
journal) is flagged for a **user PDF fetch (#116)** to transcribe its new
numeric tuples. This is a tooling limit, not a network outage — I did not
retry-loop or background-and-poll (per the anti-hang rule).

## Candidate list + verdicts

| # | Paper | Class | New capability? | Verdict | Access |
|---|---|---|---|---|---|
| 1 | **Roberts-Tsoukkas & Ross 2026**, "Stable Prograde Earth-Moon Multi-Orbiter Cyclers via Three-Body Dynamics" (journal ext. of AAS 25-621) | ballistic cycler family | **Coverage/topology extension** — cycler families at HIGHER mass parameter incl. binary-star range; "universal stable subfamily" claim | **HIT — ADOPT (gated on user PDF fetch)** | free (VSGC PDF) |
| 2 | **Fu, Wu, Gong & Shi 2026**, "Data Mining-Based Cislunar Escape-Family Analysis," arXiv:2601.11881 | escape families (CR3BP + BCR4BP) | clustering-based family extraction + cross-fidelity persistence METHOD (Track-B-adjacent) — but output is escapes, not cyclers | **BACKGROUND (method watch)** | free (arXiv) |
| 3 | **Merrill, Kulik, Bryan & Savransky 2025**, "Mass-Optimal Low-Thrust Forced Periodic Trajectories in the Earth-Moon CR3BP," arXiv:2502.05140 | low-thrust forced-periodic | low-thrust reachable-set ordering; forced-periodic = out of cycler scope; primitives already held | **BACKGROUND (Track-A ref)** | free (arXiv) |
| 4 | **Tagliaferri, Blazquez, Acciarini & Izzo 2024** (ESA-ACT), "Global Optimization for Trajectory Design via Invariant Manifolds in the EM CR3BP," arXiv:2405.18916 | halo↔vert-Lyapunov transfers | MBH+SQP manifold transfers — every ingredient already in toolkit | **BACKGROUND (no gap filled)** | free (arXiv) |
| 5 | **Takubo, Landau & Anderson 2023/24**, "Automated Tour Design in the Saturnian System," arXiv:2210.14996 (CMDA 136:8) | moon-tour transfer | VILM-database + DP tour assembly; endgame primitive already held (#179); not a cycler | **BACKGROUND (lane backfill)** | free (arXiv) |

### Lighter items surfaced by forward-citation expansion (no dedicated digest)

| Paper | Why noted | Verdict |
|---|---|---|
| **Braik & Ross 2026**, "Orbital Networks in the Three-Body Problem," arXiv:2605.31543 | reachable-set family-accessibility network; (3,2)-cycler is the dominant hub | **ALREADY MINED + ADOPTED** (`2026-06-13-braik-ross-2026-orbital-networks-mining.md`, #230) — not re-mined |
| **"Identifying Fixed Points … High-Order Transfer Map," arXiv:2509.12671** | DA/HOTM exhaustive EM PO search; hybrid DRO–Lyapunov family | **ALREADY MINED** (`2026-06-13-high-order-transfer-map-2509.12671-mining.md`) |
| **arXiv:2512.11064**, "An Open Benchmark of One Million High-Fidelity Cislunar Trajectories" (LLNL/SSAPy, Dec 2025) | 1e6 propagated cislunar trajectories, high-fid (Earth/Moon/solar gravity + SRP), 6-yr spans, public (LLNL Green Data Oasis / UDL) | **BACKGROUND — DATA ASSET.** Not discovery capability, but a ready-made corpus for the **ML false-positive flagger** and the #240 surrogate (the spec's "corpus is CR3BP-corrector data, not flyby-cycler" gap — this is high-fidelity cislunar, complementary). Note for the ML-flagger task. |
| **Li/Izzo et al.**, "Learning Optimal Control & Dynamical Structure … Diffusion Models," arXiv:2410.02976 | ML-generative global trajectory-search priors | **BACKGROUND** (downstream of Track B, #240-adjacent; same tier as Ozaki/generative-PO) |
| **arXiv:2605.23770**, "Reachability for Low-Thrust Trajectories via Maximum Initial Mass" (2026) | low-thrust reachability sibling of #3 | **BACKGROUND** (low-thrust track ref) |
| **arXiv:2507.19928**, "Periodic orbit tracking in cislunar space: finite-horizon" | control/tracking, not discovery | **SKIP** (out of scope — station-keeping/control) |

## Ranked discovery-relevance order

1. **Roberts-Tsoukkas & Ross 2026 journal (multi-orbiter)** — the ONLY genuine
   discovery-capability item: a coverage/topology extension into a mass-parameter
   regime the catalogue has never held, found by a method we already validated
   and adopted (Ross corrector). Immediately reproducible once tuples are
   transcribed. **Highest leverage by a wide margin.**
2. Fu 2601.11881 (data-mining escape families) — the only *method* with Track-B
   (prioritizer/clustering) transfer value, but applied to escapes; a technique
   to lift later, not results to use now.
3. arXiv:2512.11064 (1M cislunar benchmark) — not capability, but the most
   useful *data asset* found, for the ML false-positive flagger / surrogate.
4–6. Merrill / Tagliaferri / Takubo — all BACKGROUND: each restates a capability
   already in the earned toolkit (low-thrust primitives, MBH manifold transfers,
   VILM endgame), or produces an out-of-scope output (forced-periodic / transfers
   / tours, not ballistic cyclers).

## Headline finding

**The discovery frontier is genome/coverage-gated, exactly as the spec predicts
— and the literature confirms there is no off-the-shelf new ballistic-cycler
GENOME to adopt.** No 2024–2026 paper introduces a repeated-moon multi-rev,
broken-plane, or inclined ballistic-cycler search genome we lack (explicitly
searched; none found). The single highest-leverage discovery-capability item is
**Roberts-Tsoukkas & Ross 2026's higher-mass-parameter multi-orbiter cyclers**:
it does not give us a new genome, but it hands us **known-non-empty new
catalogue territory (cycler families at μ values we've never represented)
reachable by a corrector we already trust** — the cleanest available path to a
novel-but-sourced catalogue row, pending the user PDF fetch.

## Actions / follow-ons

- **USER FETCH (free, #116):** Roberts-Tsoukkas & Ross 2026 journal PDF
  (vsgc.odu.edu) — transcribe higher-μ (μ, C^stable, T^stable) tuples + any
  printed ICs under two-pass rescan discipline → reproduce with adopted Ross
  corrector (same-model golden) → proposed new V1 rows. (#1 priority.)
- **ML-flagger task input:** record arXiv:2512.11064 (1M cislunar trajectories,
  LLNL Green Data Oasis / UDL) as a candidate high-fidelity corpus.
- **Track-B watch:** Gong's BUAA group (data-mining cislunar) and the Fu
  clustering pipeline — revisit if/when the Track-C daemon needs an automated
  family-clustering + cross-fidelity-persistence layer.
- **No re-mine needed** for Braik-Ross (2605.31543) or the high-order transfer
  map (2509.12671) — both already have dedicated notes.
- **Clean negatives (do not redo without a new method):** no new ballistic
  broken-plane / inclined / repeated-moon cycler genome 2024–2026; no new
  Earth-Mars cycler family with sourced rows; no external Liang-CGE citers
  beyond the group's own follow-on.
