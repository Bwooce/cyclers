# Cycler Orbit Finder — Project Definition & High-Level Design

**Audience:** an implementing agent/engineer (e.g. Claude Code) starting from an empty repo.
**Goal of this doc:** enough definition and structure to scaffold the project and begin building in priority order, with concrete validation anchors so correctness is checkable at every step.

> This document is the **canonical specification** as authored by the project owner. It is preserved verbatim. Implementation-level decisions made during planning (toolchain, scope of the first slice, etc.) live in [overview.md](overview.md) and per-phase plans under `phases/`.

---

## 1. Purpose

Build a tool that **systematically finds and verifies planetary cycler trajectories** — heliocentric orbits that repeatedly re-encounter two or more planets on a fixed schedule, maintained by gravity assists with little or no propellant. Primary targets: Earth–Mars cyclers (well-documented, used as validation) and the under-explored **Venus–Earth–Mars (VEM)** triple-cycler space (where genuinely new orbits may exist).

The tool must do what blind optimisation cannot: use **structured, Tisserand-guided sequence enumeration** with **V∞-matched leg construction**, then optimise and verify. A prior prototype confirmed that throwing a global optimiser at a hand-picked flyby sequence does **not** find ballistic cyclers — that failure mode must be designed out.

---

## 2. Goals and non-goals

**Goals**

- Enumerate feasible flyby sequences for a chosen body set and period bound, ranked by energetic feasibility.
- Construct patched-conic cycler candidates (circular-coplanar first, then real ephemeris).
- Quantify each candidate: total maintenance ΔV, per-encounter V∞, flyby feasibility (bend vs max), radial span, period/resonance error, and a "taxi" hyperbolic-rendezvous cost.
- Optimise promising seeds to (near-)ballistic closure.
- **Verify** closure independently (multi-lap propagation; Lambert cross-check; optional GMAT round-trip).
- **Check novelty** against a catalogue of published cyclers.
- Visualise: Tisserand graphs, pork-chops, trajectory loops, multi-lap periodicity.

**Non-goals (v1)**

- Full n-body / low-thrust optimal control (patched-conic + impulsive flybys only; low-thrust is a stretch goal).
- Spacecraft/subsystem design (that lives in the mission concept docs).
- Crewed-cycler taxi mission design (compute the taxi V∞ as a metric only).

---

## 3. Background the implementer needs (concise)

- **Synodic period**: time between successive same-geometry alignments of two planets. Earth–Mars ≈ 2.135 yr; Earth–Venus ≈ 1.599 yr.
- **Cycler period** is an integer multiple *k* of a synodic period. A VEM cycler must satisfy multiple synodic commensurabilities; the natural lowest beat is **≈ 6.4 yr (3× E–M ≈ 4× E–V)**, but longer commensurabilities (e.g. 12.8 yr, 32 yr) also support closure — Jones, Hernandez, Jesick (2017) report VEM triple cyclers at the 6.4-yr and 12.8-yr periods (see §9.1).
- **V∞** (hyperbolic excess speed) at a flyby = |spacecraft heliocentric velocity − planet velocity|. A ballistic flyby preserves |V∞| and rotates its direction by up to a max bend δ_max, where `sin(δ_max/2) = 1 / (1 + r_p·V∞²/μ_planet)` (r_p = min safe flyby radius, μ_planet = planet GM).
- **Closure**: a cycler repeats in the synodic rotating frame. After total period T, the spacecraft's outgoing state must match its initial state rotated by the frame's advance (Earth's mean motion × T). Residual = the maintenance ΔV.
- **Key physical fact**: Earth and Venus are strong steerers (~60° bend at 7 km/s); **Mars is weak (~24°)** — so Venus assists must do the heavy shaping in VEM cyclers.

---

## 4. Architecture

```
cyclerfinder/
├── core/
│   ├── constants.py      # μ_sun, AU, planet μ/radii/a, safe-altitude defaults
│   ├── ephemeris.py      # unified planet-state interface: circular-coplanar | astropy(JPL)
│   ├── lambert.py        # universal-variable solver (single + multi-rev) + lamberthub cross-check
│   ├── flyby.py          # max bend, ballistic feasibility, powered-flyby ΔV
│   ├── kepler.py         # universal-variable propagator (for sampling/plots/verification)
│   └── frames.py         # synodic rotating-frame transforms
├── search/
│   ├── tisserand.py      # V∞-graph contours; pairwise linkable-region detection
│   ├── resonance.py      # synodic periods; k-synodic candidates; multi-body beat finder
│   ├── sequence.py       # enumerate body orderings + rev counts, filtered by Tisserand feasibility
│   ├── construct.py      # V∞-matched patched-conic build of a sequence at given times
│   └── optimize.py       # global (pygmo/scipy DE) + local (NM/SLSQP) over seed parameters
├── model/
│   ├── cycler.py         # Cycler / Leg / Encounter dataclasses; closure residual
│   └── score.py          # metrics: ΔV, V∞, radial span, period error, taxi cost
├── verify/
│   ├── propagate.py      # multi-lap propagation; periodicity / re-encounter miss
│   ├── crosscheck.py     # Lambert vs lamberthub (izzo, gooding)
│   └── gmat_bridge.py    # (stretch) export script + parse results for GMAT verification
├── data/
│   └── catalog.py        # known cyclers (Aldrin, McConaghy 2-synodic, Russell–Ocampo set) + novelty check
├── viz/
│   └── plots.py          # tisserand, porkchop, trajectory loop, multi-lap rotating frame
├── cli.py                # orchestration / commands
├── configs/*.yaml        # run configs
└── tests/                # pytest; validation anchors below
```

**Reference prototypes (proven, port these):** the user has working scripts — `lambert.py` (universal-variable solver, validated to 0.000 m/s vs lamberthub izzo & gooding), and `scan.py`–`scan6.py`/`scan_vem.py` (Tisserand graph, pork-chops, flyby feasibility, closure optimisation, multi-lap propagation). Port their validated logic into the module structure above rather than rewriting from scratch.

---

## 5. Algorithm pipeline

1. **Resonance** (`resonance.py`): compute synodic periods for the body set; list integer-k periods and multi-body beat periods (the VEM 6.4-yr case). These set the candidate total periods T.
2. **Tisserand enumeration** (`tisserand.py` + `sequence.py`): build V∞ contours per body; find energy regions where consecutive bodies are linkable at a common V∞; enumerate flyby sequences (body order + per-leg revolution count) that stay within those regions and within a V∞ cap. **This replaces blind structure-guessing** — only energetically plausible sequences proceed.
3. **Seed construction** (`construct.py`): for each sequence at the resonance period, build legs with **Lambert**, placing planets from `ephemeris.py`. Prefer constructions that satisfy flyby V∞-magnitude continuity *by design* (search encounter times so |V∞_in| ≈ |V∞_out| at each node).
4. **Score & filter** (`score.py`): rank seeds by total maintenance ΔV, max V∞, radial span (must reach all target bodies), period error, taxi cost. Keep top-N.
5. **Optimise** (`optimize.py`): global DE then local polish over each surviving seed's free parameters (encounter times, phases), minimising closure ΔV with a low-V∞ regulariser and bend-feasibility penalties. Period either locked to k·synodic (strict cycler) or lightly freed (then reported).
6. **Refine on ephemeris** (`ephemeris.py` astropy/pykep backend): re-solve the best circular-coplanar solutions on real planet states; report the change.
7. **Verify** (`verify/`): multi-lap propagation in the rotating frame (laps must overlap → periodic); Lambert cross-check; optional GMAT round-trip.
8. **Novelty** (`catalog.py`): compare closing solutions to the known-cycler catalogue (period, sequence, V∞ signature) and flag matches vs candidates-for-new.
9. **Report** (`viz/`, `io`): figures + a JSON/CSV record per candidate.

---

## 6. Key interfaces (sketch)

```python
# ephemeris.py
class Ephemeris:
    def __init__(self, model: str = "circular"):  # "circular" | "astropy"
    def state(self, body: str, t_sec: float) -> tuple[np.ndarray, np.ndarray]:  # r[km], v[km/s]

# lambert.py
def lambert(r1, r2, tof, mu=MU_SUN, prograde=True, revs=0, branch="low") -> tuple[v1, v2]
def lambert_crosscheck(r1, r2, tof, **kw) -> dict   # {mine, izzo, gooding, max_diff_mps}

# flyby.py
def max_bend(mu_planet, rp_min, vinf) -> float           # rad
def flyby_dv(vin_vec, vout_vec, mu_planet, rp_min) -> float  # km/s (0 if ballistic & feasible)

# model/cycler.py
@dataclass
class Encounter: body:str; t:float; r:np.ndarray; v_planet:np.ndarray; vinf_in:np.ndarray; vinf_out:np.ndarray
@dataclass
class Cycler:
    bodies:list[str]; period:float; encounters:list[Encounter]; legs:list
    def maintenance_dv(self) -> float
    def closure_residual(self) -> float
    def radial_span(self) -> tuple[float,float]
    def max_vinf(self) -> float

# search/optimize.py
def find_cyclers(bodies, k_synodic, vinf_cap, n_keep=20, ephem="circular") -> list[Cycler]
```

---

## 7. Tech stack & environment

- **Python 3.11** recommended (lets `pykep`/`pygmo` install; the project's earlier sandbox was 3.12 where pykep had no wheel). Provide a 3.12 fallback path using `astropy` + `scipy` only.
- Core deps: `numpy`, `scipy`, `astropy` (JPL ephemerides), `lamberthub` (independent Lambert cross-check), `matplotlib`, `pyyaml`, `pytest`.
- Optional: `pykep`/`pygmo` (global opt + MGA), `poliastro` (cross-check).
- External (stretch): **NASA GMAT** for gold-standard verification via `gmat_bridge.py` (generate a GMAT script, run headless, parse output).
- Config-driven runs (YAML); results to `./out/<run-id>/` as JSON + figures.

---

## 8. Milestones (build in this order; each ends with passing tests)

- **M0 — scaffold:** repo, packaging, CI, pytest harness, constants.
- **M1 — core mechanics:** `ephemeris` (circular), `lambert` (single-rev) + `kepler`. *Gate: Lambert matches lamberthub izzo & gooding < 1e-3 m/s on 3 test legs.*
- **M2 — flyby + maps:** `flyby`, `tisserand`, `resonance`. *Gate: synodic E–M = 2.135 yr, E–V = 1.599 yr; VEM beat ≈ 6.406 yr; Mars max bend ≈ 24° at 7 km/s.*
- **M3 — model + construction:** `cycler`, `frames`, `construct`; reproduce the **Aldrin cycler** (a ≈ 1.66 AU, e ≈ 0.41, P = 1 synodic) and a **2-synodic E–M cycler**. *Gate: closure residual computed correctly; Aldrin elements reproduced.*
- **M4 — enumeration + scoring:** `sequence`, `score`; Tisserand-filtered sequence lists. *Gate: enumerator rejects energetically infeasible sequences.*
- **M5 — optimisation:** `optimize`; rediscover a published **2-synodic E–M ballistic cycler from scratch**. *Gate: matches published low-V∞ values (~5.65 km/s Earth, ~3.05 km/s Mars) within tolerance.*
- **M6 — ephemeris + verify:** astropy backend, multi-lap `propagate`, `crosscheck`. *Gate: best E–M cycler verified periodic over ≥3 laps.*
- **M7 — catalogue/novelty:** `catalog`; flag known vs candidate-new. *Gate: correctly tags the rediscovered E–M cyclers as known.*
- **M8 — VEM campaign + UX:** run the enumerator on `[Venus, Earth, Mars]` at the 6.4-yr beat; `cli`, `viz`, reporting, docs.
- **Stretch — GMAT bridge & low-thrust** continuation.

---

## 9. Validation anchors (hard numbers to test against)

- Lambert vs lamberthub (izzo, gooding): max Δ|v| < 1e-3 m/s.
- Earth–Mars synodic 2.135 yr; Earth–Venus 1.599 yr; VEM beat 3×E–M ≈ 4×E–V ≈ 6.40 yr.
- Aldrin cycler (Byrnes/Longuski/Aldrin 1993, as tabulated in Rogers et al. 2012): a ≈ 1.60 AU, e ≈ 0.393, perihelion ≈ 0.97 AU, aphelion ≈ 2.23 AU, E→M leg ≈ 146 d. (Earlier spec drafts cited a=1.659 / e=0.41 — see §9.1 for the reconciliation.)
- Gravity-assist max bend at V∞ = 7 km/s, computed with `r_p = R_eq + 300 km` per `constants.PLANETS[code].safe_alt_km`: Mars ≈ 22.1°, Earth ≈ 66.6°, Venus ≈ 61.4°. (Earlier spec drafts gave Earth/Venus ≈ 60–63° and Mars ≈ 24°; these used inconsistent per-body altitudes — see §9.1.)
- Published 2-synodic E–M cycler encounter speeds ≈ 5.65 km/s (Earth), 3.05 km/s (Mars) — the optimiser should be able to reach this class, not a high-V∞ degenerate one.
- Degenerate-solution guard: reject/penalise solutions whose "ballistic" closure relies on V∞ > ~11 km/s (the known cheat mode).

### 9.1 Numerical reconciliation

The original §9 anchor numbers (Aldrin elements and bend angles) were
revisited during M2/M3 implementation. Three discrepancies surfaced:

1. **Aldrin a/e.** The original a=1.659, e=0.41 is a resonance-construction
   choice (forcing P_cycler = T_syn yields a=1.6582 AU). It is internally
   inconsistent with the also-quoted 146-day E→M leg, which corresponds to
   a different ellipse. The literature consensus (Rogers 2012 Table 1,
   citing Byrnes/Longuski/Aldrin 1993): a=1.60 AU, e=0.393, perihelion
   ≈ 0.97 AU, aphelion ≈ 2.23 AU, E→M ≈ 146 d. M3's patched-conic
   constructor reproduces these to ±0.002 on every element.

2. **Bend anchors.** The original (Mars ≈ 24°, Earth/Venus ≈ 60–63° at
   V∞=7) cannot be satisfied by any single r_p_min choice: Mars 24°
   implies altitude −60 km (subsurface); Earth 60–63° implies 1056–1757 km.
   Russell 2004 used 200 km Earth altitude (computed 67.1°). The
   project adopts an unambiguous `R_eq + 300 km` per
   `constants.PLANETS[code].safe_alt_km`, yielding Mars 22.1°,
   Earth 66.6°, Venus 61.4°. Tests in M2 assert these values.

3. **Attribution.** The VEM triple-cycler paper cited in §16.4 is
   AAS 17-577 by Jones, Hernandez, Jesick (JPL), not Longuski et al.

Full investigation: `docs/errata-investigation.md`.

---

## 10. Risks & known-hard parts

- **Search landscape is brutal.** Mitigate with Tisserand-guided enumeration + V∞-matched construction; never rely on blind DE over a guessed sequence (proven to fail at ~26 km/s).
- **Multi-rev Lambert** branch selection — use lamberthub for multi-rev; test against known cases.
- **Degenerate high-V∞ "closures"** — guard explicitly (see anchors).
- **VEM strict closure** is genuinely open research — frame M8 as "search + report best candidates," not "guaranteed novel cycler." Any closing VEM solution must pass GMAT verification and catalogue check before any novelty claim.
- **pykep on Python 3.12** — pin 3.11 or use the astropy-only path.
- **Closure frame correctness** — unit-test the rotating-frame transform; a wrong frame silently fakes/breaks closure.

---

## 11. Definition of done (v1)

A CLI run `cyclerfinder find --bodies E,M --k 2 --vinf-cap 7` returns ranked, verified Earth–Mars cyclers with figures and a JSON record, correctly tagging known ones; and `--bodies V,E,M --period beat` runs the VEM enumeration end-to-end, returning ranked candidates with honest feasibility/closure metrics and novelty flags — no candidate asserted as new without GMAT verification + catalogue check.

---

## 12. Review-driven refinements (v1.1 — supersedes noted earlier text)

A technical review raised four valid issues at the idealized→real transition. All are accepted and folded in here.

**(a) The periodicity paradox — two optimisation modes (amends §5 step 6, §8 M6).** Circular-coplanar gives exact geometric closure; real eccentricity/inclination breaks it, so a cycler that closes perfectly idealized will drift over decades. `optimize.py` and `model/cycler.py` therefore support two targets:

- *Idealized mode:* strict periodic closure over one period T (used to discover geometries).
- *Ephemeris mode:* minimise the **summed trajectory-correction-manoeuvre (TCM) ΔV over a finite horizon of 3–5 laps (~20–30 yr)**, not perfect closure. This is the realistic maintenance metric and is what gets reported for any flight candidate.

**(b) Multi-rev Lambert branching (amends §6, §10).** For N ≥ 1 revolutions there are two branches (high/low energy) up to N_max, i.e. several discrete valid trajectories for the same dates. The Lambert solver returns **all** valid solutions for a TOF; `n_revs` and `branch` become discrete nodes in the search tree, evaluated by the sequence ranker. Revised interfaces:

```python
def lambert(r1, r2, tof, mu=MU_SUN, prograde=True, max_revs=2) -> list[LambertSolution]
# LambertSolution: {n_revs:int, branch:str, v1:np.ndarray, v2:np.ndarray}

@dataclass
class LegCandidate:
    departure_body: str; arrival_body: str
    n_revs: int; branch: str            # "low" | "high"
    vinf_in: np.ndarray; vinf_out: np.ndarray
```

**(c) Dynamic ephemeris frame + tolerant verification (amends §4 `frames.py`, §9, §10).** Define the verification frame as a **non-uniform rotating frame anchored to instantaneous Sun–Earth (or Earth–Mars) positions**, not a constant angular velocity. A perfectly optimised real-ephemeris trajectory will still show small "geometric breathing" in this frame; periodicity tests must set a non-zero tolerance and check **bounded** lap-to-lap drift, not exact overlap.

**(d) Hard constraints, not soft regularisers (amends §5 step 5).** Replace the low-V∞ weight with a **constrained** objective: minimise total maintenance ΔV (including powered-flyby components) subject to **hard inequalities** V∞ ≤ V∞_cap and r_p ≥ r_p_min at every flyby (barrier/penalty or an NLP solver such as SLSQP). This structurally removes the high-V∞ "free closure" cheat rather than discouraging it.

**Long-term verification interface (amends §4 `verify/propagate.py`):**

```python
def verify_long_term_stability(cycler: Cycler, n_laps: int, ephem: Ephemeris) -> dict:
    """Propagate continuously through n_laps on real ephemeris.
    Returns {'stable': bool, 'total_tcm_dv': float, 'per_lap_dv': list[float],
             'max_drift_km': float}."""
```

### 12.1 The idealized→ephemeris bridge: planetary phase alignment (new module `search/phase_match.py`)

This answers the central transition question. The idealized solution does **not** give calendar dates — it gives a *relative-phase signature*: the inter-planet angular configuration required at each encounter, plus the leg times. The bridge is a four-step phase-matching procedure:

1. **Extract the phase signature** from the idealized cycler: the required heliocentric longitude separations between the bodies at departure and at each encounter (e.g. "Venus +X° from Earth at departure, Mars +Y° at first Mars encounter"), and the leg durations.
2. **Search real dates for matching geometry.** Scan the astropy/JPL ephemeris across a multi-decade launch range for epochs that minimise the mismatch to that signature. Because the three synodic frequencies are incommensurate and orbits are eccentric/inclined, good windows recur **near the beat period (~6.4 yr) but each is slightly different in quality** — so this is itself a small optimisation over launch epoch, returning a ranked list of candidate real launch dates.
3. **Seed and re-optimise in ephemeris mode.** Use (matched epoch + idealized leg times) as the initial guess for the ephemeris-mode optimiser (§12a), minimising finite-horizon TCM ΔV. Real 3D states are used directly; eccentricity and the ~1.85°/3.4° inclinations of Mars/Venus are absorbed into small TCMs and into the **flyby b-plane** (each 3D gravity assist has two steering DOF — turn angle and out-of-plane node — which pays most of the plane-change for free).
4. **Over-generate and filter.** The idealized search is a *generator* of candidate geometries; some will have no good real window in the target era. Phase-match **all** surviving idealized seeds and keep those whose best real window yields acceptable horizon TCM. This filter is what turns abstract cyclers into ones with actual, schedulable launch dates — and the matched epochs are themselves a mission-planning deliverable (the real launch-window calendar).

```python
# search/phase_match.py
def phase_signature(cycler: Cycler) -> dict          # required inter-body angles + leg times
def find_real_windows(signature, ephem, date_range, n=10) -> list[Epoch]  # ranked launch epochs
def to_ephemeris_seed(cycler, epoch) -> dict          # initial guess for ephemeris-mode optimise
```

**Milestone impact:** split **M6** into *M6a — idealized closure verification* and *M6b — phase-match + ephemeris-mode TCM minimisation over 3–5 laps*, with a gate that the reported flight candidate's horizon TCM is bounded and its frame-drift stays within tolerance.

---

## 13. Search-space decomposition, coverage & resumability

The search space is unbounded, so the design must make an unbounded run *productive and provably-covering* rather than aimless. The decomposition is by **discrete structure, not by time**.

### 13.1 The unit cell

The atomic unit of the search is a **structural cell**:

```
cell = (body_set, flyby_sequence, period_k, per_leg_revs[], per_leg_branch[])
```

— e.g. `({V,E,M}, "E-V-M-E-M-E", k=3, revs=[0,0,1,0,1], branch=[low,low,low,low,low])`. Everything else (encounter epochs, phases) is *continuous* and lives **inside** a cell as a bounded sub-problem. The discrete structure is the combinatorial part; the timing is a low-dimensional bounded optimisation per cell.

Time windows do **not** segment discovery — a given cell's cycler recurs every era, so windowing discovery would re-find the same families repeatedly. Time enters only at realisation (§12.1 phase-matching), as a bounded launch-epoch sub-search per idealized solution.

### 13.2 Iterative-deepening frontier (this is what makes "forever" safe)

Cells are enumerated by increasing complexity under three caps:

- `L_max` — max encounters in the sequence,
- `k_max` — max period in synodic multiples,
- `N_max` — max revolutions per leg.

Under fixed caps the cell set is **finite and fully enumerable**. The search runs as iterative deepening: complete the frontier at the current caps, then raise a cap and continue. At any moment you can state precisely what has been *exhaustively structurally covered* and what has not. "Run forever" = monotonically expanding this frontier — every cell is eventually visited, and progress is auditable, not luck.

### 13.3 Tisserand pruning (shrinks the effective space by orders of magnitude)

Before any continuous search, each cell is tested against the V∞ graph: if consecutive bodies in the sequence cannot be linked at a common V∞ within `V∞_cap`, the cell is **energetically impossible** and is discarded without optimisation. The vast majority of cells die here, so compute is spent only on viable structures.

### 13.4 Making the inner (timing) search near-deterministic

Within a surviving cell, do **not** rely on a single blind optimiser run (the proven failure mode). Instead:

- Use free-return / resonance construction so most timing parameters are fixed by the cell's structure, leaving a **low-dimensional root-find / bounded grid**.
- Cover the remaining continuous DOF with a **fixed, reproducible multi-start grid** plus local polish, so coverage within a cell is systematic, not stochastic.

This yields *practical* completeness: every structure attempted, every viable structure densely sampled. (Honest limit: this is not a formal proof of exhaustiveness — a continuous global optimum can still hide in a narrow basin. It is the standard, and as strong as this problem class allows.)

### 13.5 Deduplication & novelty (what makes a long run accumulate value)

Every closing solution is reduced to a **canonical signature**:

```
signature = (normalised body sequence, period in synodic units,
             ordered tuple of encounter V∞ magnitudes,
             leg (a,e) set)   — all rounded to tolerance
```

Two solutions with matching signatures are the same cycler. A persistent **catalogue** (seeded with Aldrin, the Russell–Ocampo set, McConaghy's 2-synodic, etc.) is checked on every hit: known → logged as a re-find; not present → flagged **candidate-novel** and queued for GMAT verification + literature confirmation. This is why running indefinitely is productive: it grows the catalogue and surfaces only genuinely new structures.

### 13.6 Work queue, checkpointing, parallelism

The cell enumeration is a **work queue**; each cell is an independent job → embarrassingly parallel across cores/machines. A persistent **ledger of completed cells** (with their results/signatures) makes the run **resumable and non-redundant** — interrupting and restarting never repeats work, and the frontier state is always known.

### 13.7 Prioritisation (get novel results sooner)

Order the frontier to front-load the under-explored, higher-novelty regions rather than re-confirming well-mapped ones:

1. VEM and other multi-body cells before pure Earth–Mars (the latter are largely catalogued).
2. Longer periods (`k ≥ 4`) that human-mission studies truncated for transit-time reasons but a robotic cycler tolerates.
3. Low-thrust-assisted and modest-inclination variants (stretch).

So: a continuously-running finder is well-posed — **structure-tiled, Tisserand-pruned, frontier-deepened, signature-deduplicated, checkpointed and parallel** — and will surface novel cyclers as it goes, with a clear and honest account of exactly how much of the space has been covered at any time.

### 13.8 Execution recipe (cell IDs, ledger, deepening loop)

**Cell ID** — a deterministic, sortable string so every cell is addressable and dedupable:
`{bodyset}|{sequence}|k{K}|r{revs}|b{branches}` e.g. `VEM|E-V-M-E-M-E|k3|r00101|blllll`.

**Ledger** — an append-only store (SQLite or JSONL) that is the single source of truth and makes runs resumable/parallel:
`cell_id, status[pending|pruned|searched|solved|failed], n_solutions, best_dv, signatures[], validation_level, t_done, host`. On start, load the ledger and skip/claim cells atomically; interrupting never repeats work.

**Deepening loop:**

```
for L in 3..L_max:
  for k in 1..k_max:
    for cell in enumerate_cells(bodyset, length=L, k=k, N_max):
      if ledger.has(cell): continue
      if not tisserand_feasible(cell, vinf_cap): ledger.mark(cell,"pruned"); continue
      sols = inner_search(cell)                 # free-return seed + fixed grid + local polish
      for s in sols:
          level = validate(s)                    # §14 gauntlet (auto stages)
          ledger.record(cell, s, signature(s), level)
    checkpoint()
raise caps; repeat                                # monotonic frontier expansion
```

Workers pull cells from the queue; the ledger coordinates. Order the frontier per §13.7 (VEM and long-period cells first).

---

## 14. Candidate validation pipeline (the gauntlet, V0–V5)

Every emitted candidate carries a **validation level** = the highest gate it has passed. Gates run cheapest-first; a candidate is only as trustworthy as its level. V0–V3 are automatic; V4 is batched for promising candidates; V5 is human.

| Level | Gate | Check | Tooling |
|-------|------|-------|---------|
| **V0** | Internal consistency | hard constraints met (V∞ ≤ cap, r_p ≥ r_p_min, bend ≤ max); V∞ magnitude preserved across each flyby; closure residual ≤ tol (idealized) | in-house |
| **V1** | Solver cross-check | every leg re-solved with **lamberthub izzo + gooding**, agreement < 1e-3 m/s; full trajectory re-propagated with the **Kepler** propagator (not the Lambert that built it), planet positions met < tol | lamberthub, kepler.py |
| **V2** | Multi-lap periodicity | ≥3 continuous laps; **bounded** drift in the dynamic rotating frame (tolerant of geometric breathing) | propagate.py |
| **V3** | Ephemeris realisation | phase-matched to a real launch window; ephemeris-mode horizon TCM over 3–5 laps (~20–30 yr) bounded and within ΔV budget | astropy backend |
| **V4** | High-fidelity external | **independent codebase + ephemeris** (NASA GMAT, or Tudat/pykep n-body) reproduces trajectory and maintenance ΔV within tol | GMAT bridge |
| **V5** | Novelty + expert review | canonical signature misses catalogue **and** literature; human expert review; ideally independent reproduction by a separate group | catalog.py + human |

**Trust gating:** only **V3+** candidates are "credible"; only **V5 + catalogue/literature miss** may be called a *discovery* or submitted for publication. The auto-pipeline tags V0–V3 on every hit; V3 candidates that miss the catalogue are queued for V4 (GMAT); V4 passers go to human V5. This is the spit-out → trust ladder, and it is what protects against publishing a re-derived or numerically-faked "novel" cycler.

---

## 15. Dissemination — open catalogue & a `cyclers.space` library

**Recommended, with discipline.** There is no single canonical public database of cycler trajectories — they're scattered across papers — so an open, reproducible library is a genuine community contribution and a credibility/visibility asset that pairs naturally with open-sourcing the finder (the catalogue *is* the tool's output).

Design it conservatively:

- **Generated, not hand-curated:** the site is a static build from the validated ledger (data → static pages), open-source repo, cheap to host (e.g. static hosting) and low-maintenance.
- **Honesty front-and-centre:** every entry shows its **validation level (V0–V5)**, orbital elements, V∞ signature, figures, and **references / prior art**. Publish only **V3+**. Tag each entry explicitly as *known (cite source)*, *verified candidate*, or *published-elsewhere* — **never present a re-derived known cycler as a discovery** (recall the 2-synodic Earth–Mars cycler we closed is well-known).
- **Start minimal, grow later:** v1 = reproducible dataset + static catalogue pages; later add interactive Tisserand exploration and a 3D orbit viewer.
- **Practical:** register the domain (check availability first), host as a static site from the repo, and treat the catalogue schema as the same record the finder emits — so discovery → validation → publication is one continuous, auditable data flow.

The payoff: the finder, the validation gauntlet, and the public library become a single pipeline — candidates flow in, earn a validation level, and only sufficiently-verified, prior-art-checked entries surface publicly, with full reproducibility.

---

## 16. Catalogue schema, identity matching & attribution

One record type flows through everything — finder output, the search ledger, the validation gauntlet, and the public site are all the *same object*. Below: the schema, the canonical signature that gives a cycler a stable identity, the matching that collapses accidental re-derivations, and the attribution model that keeps credit correct.

**File locations.** The canonical catalogue is at [`data/seed_cyclers.yaml`](../data/seed_cyclers.yaml) — sole source of truth. See [`data/README.md`](../data/README.md) for conventions, attribution rules, the `primary:` schema extension for non-heliocentric entries, and the regenerable cross-reference table command (`scripts/render-catalogue.py`). See [`data/OUTSTANDING.md`](../data/OUTSTANDING.md) for the long-form research-questions / open-source-access log.

### 16.1 Shared catalogue record

```jsonc
{
  "id": "cyc-emem-k2-0007",              // stable human slug
  "signature": "sha1:…",                 // canonical identity hash (see 16.2)
  "signature_fields": {                  // the basis of the hash (human-readable)
    "bodies": ["E","M"],
    "sequence_canonical": "E-M-E-M",     // cyclic-canonical flyby order (16.2)
    "sense": "outbound",                 // direction / escalator branch
    "period": {"pair":"E-M","k":2,"years":4.27},
    "vinf_multiset_kms": [["E",5.65],["E",5.65],["M",3.05],["M",3.05]], // sorted, binned
    "leg_elements_AU": [{"a":1.66,"e":0.41}]            // sorted multiset, binned
  },
  "trajectory": {
    "model": "ephemeris",                // circular | ephemeris
    "encounters": [{"body":"E","t_rel_days":0,"vinf_in":null,"vinf_out":[…],
                    "bend_deg":46.5,"bend_max_deg":61.9,"rp_km":7000}, …],
    "legs": [{"from":"E","to":"M","tof_days":146,"n_revs":0,"branch":"low",
              "a_AU":1.66,"e":0.41}, …],
    "radial_span_AU": [0.97, 2.34]
  },
  "metrics": {
    "maintenance_dv_mps_idealized": 0,
    "horizon_tcm_dv_mps": 120, "horizon_laps": 5, "horizon_years": 21,
    "max_vinf_kms": 6.2
  },
  "validation": { "level": "V4",
    "gates": {"V0":{"pass":true}, "V1":{"max_diff_mps":0.0},
              "V2":{"max_drift_km":1.2e4}, "V3":{"horizon_tcm_mps":120},
              "V4":{"tool":"GMAT","date":"…"}} },
  "reproducibility": { "finder_version":"…","config_hash":"…","seed":7,
    "run_id":"…","cell_id":"EM|E-M-E-M|k2|r0000|bllll","launch_epoch":"2032-11-02" },
  "discovery": {
    "source": "literature",              // literature | this-project | both
    "first_published": {"authors":["McConaghy","Landau","Longuski"],"year":2006,
                        "title":"Notable Two-Synodic-Period Earth–Mars Cycler",
                        "doi":"…","venue":"JSR"},
    "priority_date": "2006-01-01",       // earliest established date → owns attribution
    "our_status": "known-reproduction",  // known-reproduction | candidate-novel | verified-novel
    "rediscoveries": [{"run_id":"…","cell_id":"…","date":"…"}]
  },
  "publication": { "visibility":"public",   // public | internal
    "tag":"published-elsewhere",            // known | published-elsewhere | candidate | verified-novel
    "references":[{"authors":[…],"year":2006,"doi":"…"}],
    "assets":["fig_loop.png","fig_porkchop.png"] }
}
```

### 16.2 Canonical signature — a cycler's stable identity

A cycler's physical identity is invariant to absolute epoch/phase, to where on the loop you choose to "start," and to small numerical noise. The signature is built to be invariant to exactly those, and nothing else:

- **Sequence** → reduce the cyclic flyby string to its **lexicographically-minimal rotation** (a loop has no privileged start). Keep a separate `sense` field (outbound/inbound escalator) so direction variants are distinguishable but groupable; optionally also fold under reversal and record the sense separately.
- **Period** → integer `k` (exact) plus years.
- **V∞** → a **sorted multiset** of `(body, V∞)` pairs (rotation-invariant), each V∞ **binned to 0.05 km/s**.
- **Leg geometry** → sorted multiset of `(a, e)`, binned to 0.01 AU / 0.01.

`signature = sha1(canonical_json(signature_fields))`. The binning absorbs numerical noise; the canonicalisation absorbs the loop's rotational/reflective symmetry. Two records with the same hash are the same cycler.

### 16.3 Matching: collapsing accidental re-derivations

```python
def match(candidate, catalog):
    sig = canonical_signature(candidate)
    if sig.hash in catalog.by_hash:                 # exact (within binning)
        return ("known", catalog.by_hash[sig.hash])
    pool = catalog.filter(bodies=sig.bodies, k=sig.period.k)   # coarse prefilter
    near = [(e, signature_distance(sig, e)) for e in pool]
    near = [(e,d) for e,d in near if d < TAU_NEAR]  # weighted L1 over period,
    if near:                                         #   V∞ multiset, (a,e) multiset
        return ("probable-match-NEEDS-HUMAN", min(near, key=lambda x:x[1])[0])
    return ("novel", None)
```

Three outcomes: **exact** → log as a re-derivation (append to `rediscoveries`, inherit the existing record's attribution, do **not** create a new entry); **probable** → flag for human confirmation (this is the safeguard against a literature orbit described with slightly different rounding slipping through as "new"); **novel** → proceed up the validation gauntlet (§14). The fuzzy stage is deliberately conservative — it never auto-merges, only flags.

### 16.4 Attribution & literature ingest

The catalogue is seeded with, and continuously ingests, **published** cyclers — each with full citation — so finder hits match against prior art from day one.

- **Ingest:** a curated seed file plus ongoing literature review adds each published cycler with its parameters (as published, possibly approximate) and `first_published` citation. A signature is computed for every literature entry so it participates in matching.
- **Priority rule:** attribution goes to the **earliest `priority_date`**. If a `first_published` exists, the entry is `known-reproduction` and credit is to those authors **regardless of our finding it independently** — we never claim a cycler that has prior publication.
- **Retroactive correction (the "we missed it" case):** if a finder result was tagged `candidate-novel` and a later-ingested literature entry matches its signature, the pipeline **auto-downgrades** it to `known-reproduction`, attaches the citation, and re-tags the public entry to `published-elsewhere`. Symmetrically, if *we* verify-novel (V5) and publish, our entry becomes the `first_published` source for it.
- **Hard rule for any novelty claim:** publish as `verified-novel` only if **V5** *and* no exact match *and* no unresolved probable-match *and* a documented literature review returned nothing — all four.

**Seed catalogue (with citations):** Aldrin / Byrnes–Longuski–Aldrin (1993); Russell & Ocampo's 24 ballistic cyclers + nomenclature (2004–05); McConaghy et al. *Notable Two-Synodic-Period Earth–Mars Cycler* (2006); Niehoff VISIT cyclers; Jones, Hernandez, Jesick *Low Excess Speed Triple Cyclers of Venus, Earth, and Mars* (AAS 17-577, 2017). These anchor the matcher so the very first runs correctly tag known families (e.g. our re-derived 2-synodic Earth–Mars cycler resolves to the McConaghy entry, credited accordingly).

### 16.5 Discovery workflow — `source: this-project` entries

The §16.4 ingest path covers cyclers we *find in the literature*. The
mirror path is cyclers *we discover with the finder*. Both produce
records in the same catalogue with the same schema; the difference is
the `source:` field plus a `discovery_run:` block describing how the
finder found it.

**Lifecycle of a finder-discovered cycler:**

1. **Search emits a candidate.** The deepening loop (§13.8) finds a
   closing trajectory in some cell. The auto-validation gauntlet (§14)
   runs V0 (internal consistency) → V1 (lamberthub + Kepler cross-check)
   → V2 (multi-lap bounded drift) → V3 (ephemeris-mode TCM under
   horizon ΔV budget).
2. **Entry is written** to `data/seed_cyclers.yaml` with:
   - `source: "this-project"`
   - `our_status: "candidate-novel"`
   - `validation.level: "V3"` and the per-gate results
   - `discovery_run: {finder_version, run_id, cell_id, discovery_date}`
   - `first_published: null` (not yet published)
   - `priority_date: <discovery_date>` (revised to publication date if/when published)
   - Catalogue matcher run; exact match downgrades to
     `known-reproduction` and inherits attribution from §16.4 — entry
     never gets the `candidate-novel` tag in that case.
3. **V4 high-fidelity check** (independent codebase + ephemeris;
   GMAT / Tudat / pykep N-body) batched manually for promising
   candidates. On pass: `validation.level: "V4"`.
4. **V5 expert review** (human astrodynamicist; documented literature
   search returning nothing; ideally independent reproduction by a
   separate group). On pass: `our_status: "verified-novel"`.
5. **Publication.** When the cycler is published in a peer-reviewed
   venue with the project listed as the source:
   - `first_published: {authors: [...], year, title, doi, venue}`
   - `priority_date: <publication date>` (locks in priority date for
     the matcher; subsequent literature finds bow to this date)
   - Our entry becomes the literature-citable `first_published` record
     for any future re-derivations of the same cycler.

**Authorship of project discoveries.** Per the project's commit and
attribution rules (no AI authorship), `first_published.authors` for
discovered cyclers lists the **human contributors** who made
substantive engineering or scientific input — typically the project
owner plus any collaborators who designed the search, reviewed V5
output, or contributed to the discovering algorithm. AI-assisted prose
drafting or boilerplate code generation does not confer authorship.
The project's own self-citation lives in `CITATION.cff` at the repo
root.

**Why split `discovery_run:` from `reproducibility:`.** The
§16.1 `reproducibility:` block carries the finder configuration needed
to *re-run* the search and reproduce the bit-identical record. The
new `discovery_run:` block carries the *provenance* of the original
discovery — the run that first emitted this candidate, distinct from
the canonical-current `reproducibility:` block which a later, cleaner
re-run may overwrite. Without this split, re-running the finder would
silently overwrite the discovery's provenance.

**Workflow worked example:** see [`data/README.md`](../data/README.md)
"Recording a discovery" subsection for a step-by-step walkthrough
including an example entry.
