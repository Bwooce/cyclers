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
- **M3 — model + construction:** `cycler`, `frames`, `construct`; reproduce the **Aldrin cycler** (a ≈ 1.60 AU, e ≈ 0.393, P = 1 synodic — reconciled literature anchor; see §9 / §9.1) and a **2-synodic E–M cycler**. *Gate: closure residual computed correctly; Aldrin elements reproduced.*
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

### 12.2 Representation framework — why we don't standardize on one model

Cyclers are *concepts* (cyclic re-encounter patterns), not *orbits*. A single cycler has multiple legitimate representations, and the catalogue carries each as a distinct artifact rather than picking one canonical form.

**The three representations:**

| Representation | Where it lives | What it captures |
|---|---|---|
| **Idealized form** | `data/catalogue.yaml` `orbit_elements` (circular-coplanar; `model_assumption: circular-coplanar`) | The cycler's *signature shape* — `(a, e, perihelion, aphelion)` from the simplified model literature publishes in. This is the cycler's *identity*. |
| **Real-ephemeris instances** | `cyclers.space/src/data/windows.json` (regenerated weekly by `phase_match.find_real_windows` against JPL DE440) | Per launch window: actual departure date + real V∞ at injection. Many instances per idealized form. |
| **Mission-context analytic** | A handful of catalogue entries with `model_assumption: analytic-ephemeris` (e.g. Rogers 2012 establishment variants — Aldrin 4:3(2), 3:2(1)) | Mid-fidelity — eccentricity effects retained, not full N-body. Used when the source paper itself worked at this fidelity. |

A fourth representation — `model_assumption: cr3bp` — covers entries whose source uses Circular Restricted 3-Body dynamics (Arenstorf, future Saturnian manifolds). CR3BP entries have an incommensurable signature (Jacobi constant, not V∞) and M7 matching must not cross-compare with patched-conic entries.

**Why no standardization:**

- **Convert all to circular-coplanar** loses analytic-ephemeris eccentricity info and can't express CR3BP at all.
- **Convert all to real-ephemeris** requires per-entry M6b optimisation × multiple launch epochs; the result isn't an "orbit" but a *family of orbits*, one per real launch window. The idealized form remains useful precisely because it's *epoch-agnostic*.
- **Each representation answers a different question.** "What cycler is this?" → idealized. "When can I fly it?" → real-ephemeris instance. "What does Rogers 2012 say about its establishment cost?" → analytic-ephemeris row.

**The `model_assumption` field (§16.1 schema v2, 2026-06-01) is the structural fix.** It tells consumers which model each entry came from. Consumers:

- **M7 matching** (signature collision detection): partition by `model_assumption`; only cross-compare within partition.
- **M5 optimisation** (rediscovery from cell): consumes `circular-coplanar` only.
- **M6b real-ephemeris optimisation** (TCM budget per window): consumes `circular-coplanar` as the *seed* and produces a `real-ephemeris-instance` view stored separately.
- **Site / public catalogue** (cyclers.space): renders idealized fields from `orbit_elements` and links to the real-ephemeris launch-window table.

**Derived views, not conversion.** When a downstream consumer needs all entries in one representation (e.g. "every entry needs a 2027 launch date"), the answer is to generate a *derived view* from the catalogue, not rewrite the catalogue. `windows.json` is the prototype derived view. Future derived views (TCM budgets, c₃ tables, validation-level rankings) follow the same pattern.

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
| **V2-ballistic** | Multi-lap periodicity (ballistic) | ≥3 continuous laps; **bounded** drift in the dynamic rotating frame (tolerant of geometric breathing), evaluated **in the row's defining model** (for a circular-coplanar row that is the idealized propagation; the like-for-like model scope is recorded in the evidence, same convention as the V1 scoping) | propagate.py |
| **V2-powered** | Multi-cycle maintenance periodicity (powered) | ≥3 consecutive cycles where **(a)** every planned encounter is achieved within the documented encounter tolerance **with the documented per-cycle maintenance applied**, AND **(b)** intra-cycle drift versus the cycle's planned trajectory stays **bounded** (reset at each maneuver) | maintain.py + propagate.py |
| **V3** | Ephemeris realisation | phase-matched to a real launch window; ephemeris-mode horizon TCM over 3–5 laps (~20–30 yr) bounded and within ΔV budget | astropy backend |
| **V4** | High-fidelity external | **independent codebase + ephemeris** (NASA GMAT, or Tudat/pykep n-body) reproduces trajectory and maintenance ΔV within tol | GMAT bridge |
| **V5** | Novelty + expert review | canonical signature misses catalogue **and** literature; human expert review; ideally independent reproduction by a separate group | catalog.py + human |

**Trust gating:** only **V3+** candidates are "credible"; only **V5 + catalogue/literature miss** may be called a *discovery* or submitted for publication. The auto-pipeline tags V0–V3 on every hit; V3 candidates that miss the catalogue are queued for V4 (GMAT); V4 passers go to human V5. This is the spit-out → trust ladder, and it is what protects against publishing a re-derived or numerically-faked "novel" cycler.

> **Note — the V2 class-split and the "≥3 laps" convention.** V2 is split into
> **V2-ballistic** and **V2-powered** because a single drift metric cannot judge
> both regimes honestly. A *ballistic* cycler is meant to be geometrically
> periodic, so the rotating-frame-repeat drift over ≥3 laps is the right
> instrument. A *powered* cycler is **retargeted every cycle by design** — its
> maintenance maneuver shapes velocity, not where the planets are — so the
> cross-cycle rotating-frame-repeat metric is structurally unsatisfiable (it
> measures whether the spacecraft returned to the *same place relative to the
> incommensurately-breathing planets*, which a per-cycle-retargeted trajectory
> never does). The (a)+(b) reformulation closes the "limp into the encounters"
> loophole the other way: a powered cycler earns V2-powered only if every cycle
> actually *achieves its planned encounters with the maintenance applied* **and**
> the intra-cycle trajectory tracks its plan (drift reset at each maneuver), so a
> trajectory that quietly diverges between maneuvers and is yanked back at the
> last moment cannot pass.
>
> **Why ≥3 laps (a stated convention, not a magic number).** ≥3 laps = two
> intervals — the minimum that distinguishes *secular accumulation* (drift that
> grows lap-over-lap) from *bounded periodic breathing* (drift that oscillates
> within a band). The full geometric-modulation horizon (~7 laps / ~15 yr for an
> E–M cycler, where the slow eccentricity/nodal modulation completes) is **V3's**
> burden, not V2's. V3's phase-matched horizon-TCM gate (3–5 laps, ~20–30 yr)
> naturally extends V2-powered's per-cycle maintenance accounting to that horizon
> — V2-powered's per-cycle (a)+(b) is the per-cycle unit V3 sums and budgets.
>
> **Amendment (2026-06-07, the (c)+(a+b) class-split).** This split was adopted
> in response to task #134's evidence: the powered Aldrin E–M cycler, propagated
> against the *original* single V2 gate, drifts **~4.14e8 km over 3 laps (≈2072×
> the 200,000 km real tolerance)** — not a marginal miss but a structural one,
> because the gate measured the wrong thing for a retargeted cycler (see
> `tests/verify/test_aldrin_v2_v3_campaign.py`). Under the amended V2-powered
> gate the same physics **passes** (per-cycle encounter V∞-continuity ≤1e-6 km/s,
> intra-cycle Kepler-reprop residual ≤0.002 km, in-family maintenance ΔV
> 2.76–2.91 km/s/cycle over 3 consecutive cycles).

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

**File locations.** The canonical catalogue is at [`data/catalogue.yaml`](../data/catalogue.yaml) — sole source of truth. See [`data/README.md`](../data/README.md) for conventions, attribution rules, the `primary:` schema extension for non-heliocentric entries, and the regenerable cross-reference table command (`scripts/render-catalogue.py`). See [`data/OUTSTANDING.md`](../data/OUTSTANDING.md) for the long-form research-questions / open-source-access log.

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

  // 2026-06-01 schema v2: six additive optional field categories.
  // All optional; see data/README.md "Schema v2" for full defaults +
  // backfill rules. Consumers that omit any new field treat the entry
  // as if model_assumption="circular-coplanar" and every other v2 field
  // is null. None of these participate in §16.2 canonical signatures.
  "model_assumption": "circular-coplanar",  // FIDELITY of the numbers: circular-coplanar | analytic-ephemeris | cr3bp; default circular-coplanar
  // schema v4 (2026-06-03): cycler_class is the orbit's STRUCTURAL kind — see §16.7.
  // single-ellipse: one repeating Kepler ellipse (S1L1, Aldrin); orbit_elements (a,e) is authoritative.
  // multi-arc: different ellipse per leg (Russell generic-return E-E-M-M); no single (a,e) — per-arc lives
  //   in trajectory.segments[], cycle-level identity lives in invariants{}. Top-level a/e MUST be null.
  // non-keplerian: rotating-frame/CR3BP periodic orbit (Arenstorf); Kepler elements inapplicable, left null.
  "cycler_class": "single-ellipse",         // single-ellipse | multi-arc | non-keplerian; default single-ellipse
  "delta_v_kms": 0.0,                       // per-cycle maintenance ΔV; 0 for strict ballistic; null for near-ballistic / undetermined
  "v_infinity_leveraging_dv_kms": null,     // establishment ΔV (Rogers 2012 4:3, 3:2 variants); null when not applicable
  "fleet_size": null,                       // integer vehicle count required for the cadence; null where not stated
  "flyby_mechanics": [                      // per-encounter geometry, parallel to vinf_kms_at_encounters; null/missing when not extracted
    {"body":"E","turning_angle_deg":60,"min_altitude_km":200,"rp_km":6578},
    {"body":"M","turning_angle_deg":24,"min_altitude_km":200,"rp_km":3596}
  ],

  "trajectory": {
    "model": "ephemeris",                // circular | ephemeris
    "encounters": [{"body":"E","t_rel_days":0,"vinf_in":null,"vinf_out":[…],
                    "bend_deg":46.5,"bend_max_deg":61.9,"rp_km":7000}, …],
    "legs": [{"from":"E","to":"M","tof_days":146,"n_revs":0,"branch":"low",
              "a_AU":1.66,"e":0.41}, …],
    "radial_span_AU": [0.97, 2.34]
  },

  // orbit_elements gains six v2 children (all optional, default null).
  // The *_km pair is parallel to perihelion_au / aphelion_au for use
  // with non-Sun primaries (Earth-Moon, Jovian, Saturnian). The
  // orientation fields enable full 3D state for M6+ ephemeris-mode work.
  "orbit_elements": {
    // (existing) "perihelion_au", "aphelion_au", "a_au", "e",
    //            "inclination_deg" ...
    // schema v4 (2026-06-03): frame/units-tagged so elements scale to non-Sun
    // centres. reference_frame names the frame (cf. JPL SBDB's `equinox`);
    // the *_au fields are the heliocentric-inertial default, the *_km fields
    // the planet-centric form, selected by reference_frame. See §16.7.
    "reference_frame": "heliocentric-inertial", // heliocentric-inertial | planetcentric-inertial | rotating-synodic; default heliocentric-inertial
    "center": "Sun",            // body the elements are referenced to; default Sun (= primary for non-helio)
    "periapse_km": null,         // parallel to perihelion_au for non-heliocentric (planetcentric-inertial) entries
    "apoapse_km": null,          // parallel to aphelion_au for non-heliocentric entries
    "raan_deg": null,            // Right Ascension of Ascending Node, Ω
    "arg_periapsis_deg": null,   // Argument of Periapsis, ω
    "true_anomaly_deg": null,    // True anomaly at epoch, ν
    "epoch_iso8601": null        // ISO-8601 epoch for the anomaly; null = generic / circular-coplanar baseline
  },
  // schema v4: present ONLY for cycler_class=multi-arc. The cycle-level
  // identity descriptors a source actually publishes when no single (a,e)
  // exists (Russell). First-class + testable, not buried in notes. See §16.7.
  "invariants": {
    "aphelion_ratio": null,      // outbound-arc aphelion / inbound-arc aphelion (Russell)
    "transit_times_days": null,  // [outbound, inbound] Earth↔Mars transit durations
    "turn_ratio": null           // flyby turn-angle ratio, where published
  },
  // schema v4: present ONLY for cycler_class=non-keplerian (CR3BP / rotating-frame
  // periodic orbits, e.g. Arenstorf). Mirrors JPL's three-body periodic-orbit
  // catalog: identity is (jacobi, period, stability) + a state vector in the
  // rotating synodic frame normalised by lunit/tunit. Kepler orbit_elements stay
  // null here (structurally inapplicable). See §16.7.
  "cr3bp": {
    "jacobi_constant": null,     // dimensionless conserved energy-like quantity (identity)
    "period_nd": null,           // period in normalised time units (tunit)
    "stability_index": null,     // ≤1 stable; >1 unstable
    "mass_ratio": null,          // μ = m2/(m1+m2) of the primary pair
    "libration_point": null,     // 1–5, where applicable; null for non-libration families
    "family": null,              // halo | lyapunov | dro | figure-eight | ...
    "state_nd": null,            // [x,y,z,vx,vy,vz] in rotating synodic frame, normalised
    "lunit_km": null, "tunit_s": null  // normalisation constants for de-normalising
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

The six v2 fields added 2026-06-01 (`model_assumption`, `delta_v_kms`, `v_infinity_leveraging_dv_kms`, `fleet_size`, `flyby_mechanics`, and the orbit_elements 3D-orientation extensions) do **not** participate in the canonical signature. They are descriptive metadata, not identity carriers — including them would silently invalidate matches against pre-v2 records and against literature entries that omit them. M7 *may* additionally pre-filter the matcher pool by `model_assumption` (so a `cr3bp` finder hit doesn't match circular-coplanar literature), but this is a *pool filter*, not a signature input.

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
2. **Entry is written** to `data/catalogue.yaml` with:
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

### 16.6 Storage format — rationale & comparison to standard orbit catalogues

This section records *why* the catalogue is shaped the way it is, how it
relates to the established orbit-storage standards, and the three schema
elements added 2026-06-01 to close the gaps those standards exposed: the
OCM-aligned `trajectory{}` block, the `family{}` linkage, and the
`data_gaps[]` known-unknown marker.

#### 16.6.1 What the established standards store

| System | Representation | Epoch | Multi-segment + maneuvers | Fit for us |
|---|---|---|---|---|
| **CCSDS OCM** (Orbit Comprehensive Message, CCSDS 502.0-B-3 §6) | State *or* elements; logical blocks `TRAJ` → `PHYS` → `COV` → `MAN` → `PERT` → `OD` → `USER` (Table 6-1 order); `CENTER_NAME` / `TRAJ_REF_FRAME` are per-TRAJ-block, `TIME_SYSTEM` / `EPOCH_TZERO` are message-global | **Required** (`EPOCH_TZERO`) | **Yes** — this is its core strength | Trajectory sub-tree: excellent. Whole record: poor (no attribution/novelty/signature; would dump 80% into `USER`). |
| **CCSDS OEM** (Ephemeris Message) | Sampled state vectors + interpolation metadata | Required | Segments yes, maneuvers no | Over-specified; we don't sample states. |
| **NORAD TLE / CCSDS OMM** | Mean elements bound to the SGP4 theory | Required | No | Theory-bound; meaningless without its propagator. |
| **JPL SBDB / Horizons** | Osculating Keplerian + covariance; SPK Chebyshev kernels | Required | No | Single-body, epoch-anchored. |
| **MPC** | Osculating elements (a,e,i,Ω,ω,M,H) at epoch | Required | No | Single-body, epoch-anchored. |
| **JPL CR3BP periodic-orbit catalog** | Dimensionless rotating-frame state; **indexed by Jacobi constant; organized into families** | **None** (epoch-free) | Family continuation, no maneuvers | **Closest analog** — validates epoch-free, invariant-indexed, family-organized storage. |
| **AstDyS proper elements** (Milani–Knežević) | Time-averaged quasi-invariant (a, e, sin i) | Suppressed (averaged out) | No | Precedent for deliberately epoch-free storage. |

Two observations drive the design:

1. **Every *operational* standard requires an epoch.** OCM, TLE/OMM,
   SBDB, MPC all answer "where is this real object." Our circular-coplanar
   cycler *families* have no epoch — they are geometric/dynamical
   solutions, not tracked objects — so they legitimately fit none of
   these as a top-level container.
2. **The *scientific* catalogues that store epoch-free family solutions
   validate our approach.** The CR3BP periodic-orbit catalog indexes by a
   conserved invariant (Jacobi constant) and organizes by family; we
   index by the §16.2 canonical signature (V∞ multiset + ToF + period_k)
   and — as of 2026-06-01 — organize by an explicit `family{}` field.
   For the patched-conic / V∞ regime specifically there is *no*
   established standardized schema (the literature stores these as paper
   tables — Russell's, McConaghy's taxonomy), so this catalogue fills a
   real gap rather than reinventing one.

#### 16.6.2 Decision: borrow OCM's trajectory model, not its envelope

We **do not** migrate the catalogue to OCM. We **do** reshape the
trajectory sub-tree to mirror OCM's `TRAJ` + `MAN` decomposition, kept
inside our own YAML envelope, and treat **canonical OCM (KVN/XML) as an
*export target* for the epoch-anchored subset** (`model_assumption:
analytic-ephemeris` entries), where real interoperability lives. There is
no standardized YAML serialization of OCM (CCSDS defines KVN and XML; JSON
is an in-progress Nav-WG effort), so a YAML projection is *our* mapping —
it buys structural discipline and readability, not free interoperability.
Idealized circular-coplanar families cannot emit a meaningful OCM, so the
exporter applies only to the concrete subset.

The §16.1 record sketch already gestured at a `trajectory{}` block; the
2026-06-01 revision formalizes it and reconciles it with the on-disk
`data/catalogue.yaml` projection:

```jsonc
"trajectory": {
  // OCM metadata header
  "center": "Sun",            // OCM CENTER_NAME — the gravitational primary (== top-level primary)
  "ref_frame": "ECLIPJ2000",  // OCM REF_FRAME
  "time_system": "TDB",       // OCM TIME_SYSTEM; null for idealized circular-coplanar
  "epoch_tzero": null,        // OCM EPOCH_TZERO; null = epoch-free family (see 16.6.1)

  // === OCM TRAJ === ordered conic arcs, leg-by-leg
  "segments": [
    {"id":"out-em","from":"E","to":"M","traj_type":"keplerian-arc",
     "tof_days":154,"n_revs":0,"branch":"single","a_au":1.30,"e":0.257},
    {"id":"ret-me","from":"M","to":"E","traj_type":"keplerian-arc",
     "tof_days":null,"n_revs":null,"branch":null},   // known-unknown — see data_gaps
    {"id":"loop-ee","from":"E","to":"E","traj_type":"keplerian-arc",
     "tof_days":null,"n_revs":null,"branch":null}    // the "L1" intermediate Earth loop
  ],

  // === OCM MAN === per-encounter flyby / ΔV mechanics
  // (absorbs the v2 top-level delta_v_kms + flyby_mechanics placeholders)
  "maneuvers": [
    {"at_segment_boundary":["out-em","ret-me"],"body":"M","type":"flyby-ballistic",
     "dv_kms":0.0,"turning_angle_deg":null,"periapsis_alt_km":null}
  ]
}
```

- `traj_type` mirrors OCM's `TRAJ_TYPE`: `keplerian-arc` (conic, the
  ballistic default) or `cartesian-state` (boundary states, for
  ephemeris-mode entries).
- `branch` ∈ `single | low | high` selects the multi-revolution Lambert
  branch (`n_revs > 0` arcs have two time-of-flight branches). This is the
  field the multi-rev solver work (task #54) populates for the S1L1
  Earth-loop segment.
- `maneuvers[].type` ∈ `flyby-ballistic | flyby-powered | launch |
  arrival`. Ballistic flybys carry `dv_kms: 0.0`; the V∞-continuity
  condition at each flyby is the cycle-closure constraint. A
  `flyby-powered` maneuver is one whose geometrically *required* turn
  (`turning_angle_deg`) exceeds the *achievable* turn
  (`max_turning_angle_deg`) at its `periapsis_alt_km`; the shortfall is
  paid as a periapsis ΔV. The classic Aldrin (1L1) is powered at the
  **Earth (geocentric)** flyby: ~84° required vs ~72° achievable at a
  200 km Earth flyby (McConaghy 2002 Table 4 / dissertation).
- `maneuvers[].turning_angle_deg` is the *required* geocentric/planetocentric
  turn (sourced from the orbit geometry); `maneuvers[].max_turning_angle_deg`
  (additive, schema v3) is the *achievable* ballistic turn at the flyby.
- `maneuvers[].periapsis_alt_km` is the flyby periapsis altitude above the
  body's equatorial radius and is **configurable per orbit and per body**.
  Where a source specifies it (e.g. Aldrin's 200 km Earth flyby), record the
  sourced value; otherwise the mechanics code falls back to the conservative
  per-body default `constants.PLANETS[code].safe_alt_km` (300 km). Smaller
  altitude → tighter flyby → larger achievable turn, so this value directly
  sets `max_turning_angle_deg` and thus the powered/ballistic classification.
- **Signature integrity:** `segments` does **not** participate in the
  §16.2 hash beyond what `legs[]` already contributed (the per-leg `(a,e)`
  multiset, deduped). `vinf_kms_at_encounters[]`, `period{}`, and
  `orbit_elements{}` remain the signature carriers and are untouched.
  `maneuvers[]` is pure descriptive metadata (per §16.2's v2 rule).

#### 16.6.3 `family{}` — explicit family linkage (CR3BP-catalog pattern)

The CR3BP catalog's one clear advantage over our prior schema is
first-class **family** organization. Previously the relationship between,
e.g., the McConaghy "notable" 2-synodic variant and the S1L1 cycler lived
in prose `notes:`. The additive optional `family{}` field makes it
queryable:

```jsonc
"family": {
  "id": "s1l1-em",                 // stable family slug
  "name": "S1L1 Earth-Mars 2-synodic",
  "nomenclature": "Russell-McConaghy SnLm",
  "continuation_param": {"name":"k_synodic","value":2}  // along-family parameter (cf. Jacobi constant)
}
```

`family{}` is **not** a signature input (members are distinguished by
their own signatures); it is a grouping/navigation aid for the matcher's
human-review stage and the public site. A missing `family{}` is treated
as "ungrouped."

#### 16.6.4 `data_gaps[]` — "we don't know it yet" ≠ "not applicable"

A bare `null` is overloaded: it conflates **structurally not applicable**
(e.g. `inclination_deg` for a coplanar model — there is no value to know)
with **known-unknown** (e.g. the S1L1 return-leg ToF — a value exists in
the literature, we simply have not extracted it). These are
*semantically different* and must be distinguishable for honest
provenance and for sweep-driven lazy backfill.

`data_gaps[]` resolves this **without** reinterpreting the ~220 existing
nulls. It does not flip the meaning of a bare `null`; it adds an explicit
register on top of it. Rule:

- **A `data_gaps[]` entry = known-unknown (tracked).** A value is
  expected to exist; we have not filled it yet. This is the explicit,
  machine-readable, sweepable marker the backfill uses, and it carries
  the TODO context (`note`, `source_hint`, `todo_ref`).
- **A bare `null` with no `data_gaps[]` entry keeps its legacy meaning**
  (`data/README.md`: "not in source / not yet derived"), with genuinely
  not-applicable cases additionally flagged in the entry's `notes:` (the
  existing CR3BP / coplanar pattern — e.g. CR3BP entries have no
  Keplerian RAAN/ω/ν). Promoting such a null to an explicit `data_gaps[]`
  entry is precisely how a vague "missing" becomes a precise, actionable
  known-unknown.

```jsonc
"data_gaps": [
  {
    "path": "trajectory.segments[ret-me].tof_days",  // dotted/keyed path into this record
    "kind": "unknown",                               // unknown | uncertain | derive
    "note": "return M->E leg ToF not tabulated in the abstract",
    "source_hint": "McConaghy/Longuski/Byrnes 2002, AIAA 2002-4420, Table 2",
    "todo_ref": "#54-backfill"                        // task / issue reference
  }
]
```

- `kind`: `unknown` (no value on hand), `uncertain` (have a value but it
  is provisional / single-source), or `derive` (computable from other
  fields once a dependency lands — e.g. a leg `(a,e)` derivable once the
  multi-rev solver closes the arc).
- **Sweepable.** `data_gaps[]` is the single machine-readable source of
  known-unknowns. `cyclerfinder.data.catalog.find_data_gaps(catalog)`
  (and its CLI) enumerate every gap across the catalogue, so the gap
  inventory is a query, not a manual audit. This is the mechanism behind
  "do sweeps to identify gaps and lazy-fill."
- **Lazy fill.** When a value arrives: populate the field, update its
  `source_quotes:` per the §16.4 / README provenance rule, and remove the
  matching `data_gaps[]` entry. The no-fabrication rule (README rule 2)
  is unchanged — a known-unknown stays `null` + `data_gaps` until a real
  source backs the value.
- This does **not** weaken rule 2: previously a gap and a not-applicable
  null were indistinguishable, which understated how much is genuinely
  *missing*. `data_gaps[]` surfaces that debt explicitly.

#### 16.6.5 Backfill plan (legs[] → trajectory{})

The ~220 existing entries carry the flat top-level `legs[]` (and the v2
`flyby_mechanics` / `delta_v_kms` placeholders). Migration to
`trajectory{}` is **lazy and sweep-driven**, not a big-bang rewrite
(most entries have only partial leg data; a bulk rewrite would write
mostly `null`s and obscure which nulls are real gaps):

1. **Loader compatibility (task #64).** `catalog.py` reads
   `trajectory.segments` when present and falls back to the legacy
   `legs[]` otherwise, producing an identical canonical signature either
   way. Both forms are valid on disk during the transition.
2. **Exemplar first (task #65).** Migrate `s1l1-2syn-em-cpom` fully:
   populate the known 154-d E→M segment, add the return + Earth-loop
   segments with `data_gaps[]` markers for their unknown
   `tof_days`/`n_revs`/`branch`. This is the worked example for the rest.
3. **Sweep to inventory.** `find_data_gaps` + a per-entry "has the entry
   been migrated to `trajectory{}`?" check produce the backfill worklist.
   An entry still on `legs[]` is itself a (structural) migration gap.
4. **Lazy fill, source-gated.** Migrate remaining entries opportunistically
   — when an entry is touched for any other reason, or when a source is
   read that supplies the missing legs. Each migration: move `legs[]` →
   `trajectory.segments`, add `maneuvers[]` from `flyby_mechanics`, mark
   genuine unknowns in `data_gaps[]`, update provenance.
5. **Completion check.** The migration is "done" when no entry retains a
   top-level `legs[]` and the only `data_gaps[]` remaining are true
   known-unknowns awaiting source access (tracked, not blocking).

Backfill execution is deferred work tracked under task #65; only the S1L1
exemplar is populated immediately.

### 16.7 Scaling to multi-arc, planet-centric & n-body cyclers (schema v4)

Schema v1–v3 quietly assumed every cycler is **one repeating heliocentric
Kepler ellipse**. That is true for S1L1 and Aldrin, but false for two whole
classes the catalogue already contains as citation-only rows, and false for
everything the multi-body (M8) and planet-centric (#76) milestones will add.
Schema v4 makes the record *structurally honest* about what kind of orbit it
holds, and removes the heliocentric/2-body biases that block n-body and
moon-tour entries. **All v4 fields are additive, optional, and default to the
v3 behaviour** — a record that omits them reads exactly as before, so the
loader, site, and tests keep working until each entry is backfilled. None of
the v4 fields participate in the §16.2 canonical signature (same rule as v2).

#### 16.7.1 `cycler_class` — the structural kind

| `cycler_class` | What it is | Where identity lives | Example |
|---|---|---|---|
| `single-ellipse` *(default)* | one Kepler ellipse the vehicle repeats; flybys re-match V∞ | top-level `orbit_elements{a,e}` | S1L1, Aldrin |
| `multi-arc` | a *different* ellipse per leg (generic-return); no single `(a,e)` exists | per-arc `(a,e)` in `trajectory.segments[]` + cycle-level `invariants{}` | Russell E-E-M-M |
| `non-keplerian` | rotating-frame / CR3BP periodic orbit; not a Kepler ellipse at all | `orbit_elements.cr3bp{}` | Arenstorf E-Moon |

The class tells every consumer how to dispatch: the resonance constructor
only applies to `single-ellipse`; the validator picks which anchors are
checkable (§16.7.4); the site renders the right columns instead of a row of
em-dashes.

**Physical invariant (a loader/validation gate, no source needed):** a
`multi-arc` or `non-keplerian` record MUST NOT carry a non-null top-level
`orbit_elements.a_au`/`e`. A single semimajor axis on a multi-arc cycler is a
data error by construction — the kind of mistake that produced the original
"null but misleading" Russell rows.

#### 16.7.2 Frame/units-tagged elements (planet-centric scaling)

`orbit_elements` gains `reference_frame` + `center`, following JPL SBDB, which
names the frame at the orbit level (its `equinox`/`epoch`) and lets each
quantity carry its own units. Heliocentric entries keep `a_au`/`perihelion_au`
(frame `heliocentric-inertial`, the default); planet-centric entries use the
existing `periapse_km`/`apoapse_km` pair under frame `planetcentric-inertial`
with `center` set to the primary (Earth, Mars, Jupiter…). This is what lets an
Earth–Moon or Mars–Phobos–Deimos cycler store real elements instead of nulls.
The body registry (`core/constants.py`, currently V/E/M only) is extended with
moon `μ`/radius/ephemeris under task **#76**; the *schema* is ready now, the
*compute* lands with that milestone.

#### 16.7.3 `period.basis` — n-body beat periods

`period{pair,k,years}` is a 2-body synodic relation. An n-body cycler (V-E-M
and beyond) repeats on the **beat** of several synodic pairs, with no single
`pair`. v4 adds an optional `period.basis`: a list of `{pair,k}` whose periods
beat to `period.years`; the legacy `pair`/`k` remain valid as the 2-body
special case (and as a derived convenience when `basis` has one entry). This
aligns with the M8 `period_basis` design and makes the schema n-body-ready
immediately.

```yaml
period:
  years: 6.40
  basis:                 # n-body: beat of multiple synodic pairs
    - {pair: E-M, k: 3}
    - {pair: V-E, k: 5}
# 2-body entries keep the flat form: {pair: E-M, k: 2, years: 4.27}
```

#### 16.7.4 `invariants{}` and the `cr3bp{}` block — source-matched anchors

The point of capturing these is **testable expected outputs that match what
the source actually published** (the golden-discipline rule: anchors trace to
a source, never to our own compute):

- **`multi-arc` → `invariants{}`** carries the descriptors Russell tabulates
  when no `(a,e)` exists — `aphelion_ratio`, `transit_times_days`, `turn_ratio`
  — promoted from prose `notes` to first-class fields the validator can assert.
- **`non-keplerian` → `orbit_elements.cr3bp{}`** mirrors the **JPL three-body
  periodic-orbit catalog** exactly: identity is `(jacobi_constant, period_nd,
  stability_index)` plus a `state_nd` vector in the rotating synodic frame with
  `mass_ratio`, `libration_point`/`family`, and `lunit_km`/`tunit_s` to
  de-normalise. This replaces the Arenstorf row's "all-null + honest note" with
  the field's *real* published identity.

#### 16.7.5 Validation dispatch by class

The tiered gauntlet (data-validation-hardening plan / the Forge) dispatches on
`cycler_class` so it never applies the wrong check:

| class | reproduce via | expected-output anchors checked |
|---|---|---|
| `single-ellipse` | `construct_resonant_cycler` from `(a,e)` | V∞ multiset, `(a,e)`, period |
| `multi-arc` | multi-leg solver (no single ellipse) | V∞ multiset, period, `invariants{}` |
| `non-keplerian` | CR3BP differential corrector (M8+) | `jacobi_constant`, `period_nd`, `stability_index` |

This closes the gap that let multi-arc rows be "validated" only on V∞ while
their published invariants went unchecked, and stops the constructor from being
mis-applied to orbits that have no single `(a,e)`.

#### 16.7.6 Prior art consulted

The v4 shape deliberately tracks established orbit catalogues rather than
inventing one: **CCSDS OCM** (multi-`segment` trajectory model — already
borrowed in §16.6.2) for the per-arc decomposition; **JPL SBDB** (frame named
at orbit level, per-element units, per-element `sigma` uncertainty) for the
frame/units-tagged elements; and the **JPL three-body periodic-orbit catalog**
(`jacobi`/`period`/`stability` + synodic state vector + `lunit`/`tunit`) for
the `cr3bp{}` block. Adopting `sigma`-style published uncertainty on golden
anchors is a noted future extension for the provenance work, not part of v4.

Backfill is lazy and source-gated (as in §16.6.5): tag existing rows with their
`cycler_class` first (a mechanical sweep — single-ellipse is the default,
Russell rows → multi-arc, the 6 non-heliocentric rows → non-keplerian), then
populate `invariants{}`/`cr3bp{}` opportunistically as sources are read.

#### 16.7.7 `free_return_arcs[]` — arc-type descriptors (schema v4.1)

**Storage design: option (a) — separate `free_return_arcs[]` list.**

Russell's *arc* decomposition is genuinely distinct from the catalogue's
*encounter-segment* decomposition (§16.6.2):

- **`trajectory.segments[]`** (OCM model): encounter legs — one segment per
  planetary encounter boundary (E→M transit, M→E return, E→E intermediate loop).
- **`free_return_arcs[]`** (Russell model): Earth-to-Earth free-return arcs —
  each arc is a complete Earth-departure → Earth-arrival trajectory that happens
  to pass near Mars. A single Russell arc spans what the catalogue models as two
  or more encounter segments.

These are *different decompositions of the same physical orbit* and must not
be conflated. Design (b) — adding arc fields directly onto `trajectory.segments[]`
— was rejected because there is no 1-to-1 correspondence between Russell arcs
and encounter segments, and forcing one would fabricate structure not in the source.

**Descriptor field meanings** (sourced from Russell 2004 dissertation p.126–127):

| Field | Source quote / meaning |
|---|---|
| `arc_type` | "Each string type begins with a letter indicating either half-rev, full-rev, or generic return." g/G → `generic`, h/H → `half-rev`, f/F → `full-rev` |
| `tof_years` | "In all three cases, the first number in the parenthesis is the time of flight **in years** for that Earth-Earth leg." Applies to g/h arcs only. |
| `resonance` | "The number following the colon in the full-rev strings represent the number of revolutions by the spacecraft, thus the full-rev return is an M:N resonant orbit." (e.g. `"3:2"`, `"2:1"`, `"1:1"`). Only on f/F arcs; `null` for g/h. For full-rev arcs the TOF is determined by the M:N resonance condition, so `tof_years` is `null`. |
| `raw_descriptor` | Verbatim Russell token, e.g. `"g(1.4612,526.02,Ll)"` or `"F(3:2,82.487,180.000)"`. Uppercase letter = designated transit leg (transit times and Mars V∞ are computed from this leg). |

**Known puzzle resolved:** the second parameter in g/h descriptors (e.g. 526.02
in `g(1.4612,526.02,Ll)`) is **ψ in degrees** (the referencing angle on the V∞
sphere per Russell §2.7.3), NOT time of flight in days. The TOF is the *first*
parameter in years (1.4612 yr ≈ 533 d ≠ 526.02 d).

**Schema v4.1 field:**

```yaml
free_return_arcs:   # null when no descriptor is available (data gap, not error)
  - arc_type: generic          # generic | half-rev | full-rev
    resonance: null            # M:N string for full-rev; null for g/h
    tof_years: 1.4612          # TOF in years for g/h arcs; null for f/F
    raw_descriptor: "g(1.4612,526.02,Ll)"   # verbatim Russell token
```

**Backfill coverage (2026-06-03):** 12 entries with explicit descriptors from
Russell 2004 Tables 4.9–4.13 (pp.127–134); 3 entries with incomplete
descriptors gapped (`russell-ch4-8.165Gfh-f2`, `russell-ch4-3.77Gh3`,
`russell-ch4-5.66Gfh3`); all `russell-ocampo-*` entries gapped (Russell Ch3
tables carry AR/TR summary, not leg descriptors).

#### 16.7.8 Data-architecture design rationale (prior art, decisions recorded)

Four bodies of prior art shaped the catalogue's data architecture. The
decisions below were each made earlier in the project (§16.6, §16.7.6,
§16.7.7); this subsection records them against their sources in one place.

> **Provenance note (updated 2026-06-05):** CCSDS 502.0-B-3 and the
> Campagnola–Russell Endgame pair (AAS 09-224/09-227 preprints of the JGCD
> papers) have now been read in full — see
> `docs/notes/2026-06-05-ccsds-odm-502-mining.md` and
> `docs/notes/2026-06-05-endgame-tisserand-mining.md`. Acton 1996 remains
> cited from its bibliographic record only (paywalled). No golden-anchor
> numbers derive from these sources; they inform *structure* only.

1. **Seeds, not tracks** — *contra* SPICE SPK/DAF (Acton, C.H., "Ancillary
   data services of NASA's Navigation and Ancillary Information Facility,"
   *Planetary and Space Science* 44(1):65–70, 1996, DOI
   `10.1016/0032-0633(95)00107-7`). SPICE stores *sampled trajectories*
   (Chebyshev-fitted state histories in DAF segments) because its job is
   replaying a specific mission's as-flown path. Our job is the opposite:
   the catalogue stores the *generative seed* (elements, invariants, arc
   descriptors, encounter sequence) and the code regenerates the trajectory
   on demand against any ephemeris. Storing SPK-style sampled states would
   freeze each row to one ephemeris realisation and balloon the repo.
   SPK-kernel *references* (for cross-checking against published mission
   kernels) were a v4.2 candidate; their sourced essence — the ephemeris
   model the source's numbers were computed against — is now built as the
   `source_ephemeris` field (see §16.7.9).

2. **OCM over OEM** — CCSDS Orbit Data Messages (CCSDS 502.0-B-3, Blue
   Book, April 2023) defines four message types: OPM (single state), OMM
   (theory-bound mean elements), OEM (sampled state vectors + HERMITE /
   LINEAR / LAGRANGE interpolation metadata), and — new in issue 3 — the
   OCM (§6 of the same standard; an earlier draft of this subsection
   mis-cited it as "CCSDS 504.0-B", which is actually the Attitude Data
   Messages standard) with multi-segment trajectories and maneuver blocks.
   OEM was considered and rejected for the same reason as SPK: it is a
   *track* format (sampled states + interpolation), and we don't sample
   states. OCM's `TRAJ`/`MAN` logical-block decomposition is the model the
   `trajectory{}` block borrows (§16.6.2); the full comparison table is in
   §16.6. A projection from a fully-populated row to canonical OCM KVN/XML
   remains possible by construction (see `data/README.md`). Reading the
   standard itself (2026-06-05) validated two v4.2 choices: OCM §6.2.5.4(d)
   explicitly allows `CENTER_NAME` to change per TRAJ block (our
   `segments[].center`), and OCM metadata records `CELESTIAL_SOURCE` (e.g.
   `JPL_DE_FILES`) — the standard's twin of our `source_ephemeris`. One
   export caveat: `TIME_SYSTEM`/`EPOCH_TZERO` are message-global, so an OCM
   projection must hold one time system across segments even when centers
   differ.

3. **Descriptor as genome** — Russell 2004 §2.7/Ch.4 descriptor strings.
   Encoding a trajectory's *structure* as a short string over a small
   alphabet makes families enumerable and searchable independent of any
   numeric realisation. `free_return_arcs[].raw_descriptor` (§16.7.7)
   adopts this directly: the verbatim Russell token is first-class data,
   with `arc_type`/`resonance`/`tof_years` as its parsed projection.
   *(A secondary reference originally cited here as "Campagnola, Skelton &
   Lantoine 2014, Global Search for Gravity-Assist Trajectories" could not
   be verified in any database — no such paper exists; likely a mix-up of
   Anderson, Campagnola & Lantoine 2015, "Broad search for unstable
   resonant orbits in the planar circular restricted three-body problem,"
   Cel. Mech. Dyn. Astron. 124(2):177–199, DOI `10.1007/s10569-015-9659-7`.
   Dropped; Russell 2004 — held and verified — carries this concept alone.)*

4. **Tisserand graph** — Campagnola, S. & Russell, R.P., "The Endgame
   Problem Part 1: V∞-Leveraging Technique and the Leveraging Graph" (DOI
   `10.2514/1.44258`) and "Part 2: Multibody Technique and the
   Tisserand–Poincaré Graph" (DOI `10.2514/1.44290`), *Journal of Guidance,
   Control, and Dynamics* 33(2), 2010. Representing flyby reachability as
   a graph over (body, V∞) nodes is what the Tisserand module implements;
   running shortest-path search (Dijkstra/A*) over that graph to *propose*
   candidate sequences is a recorded Forge-pipeline enhancement candidate,
   not yet built.

> **SPICE-as-validation-instrument (added 2026-06-06):** rejecting SPICE *as
> storage* (point 1 above) does **not** reject SPICE *as an independent
> measuring instrument*. The catalogue still stores seeds, not tracks — that is
> intact. Separately, `tests/verify/test_ephemeris_crosscheck.py` (optional
> `validation` extra; skips cleanly without spiceypy) points NAIF's CSPICE
> reader (spiceypy) at the *same* DE440 kernel astropy already caches and
> compares `Ephemeris("astropy").state(body, t)` against `spkgeo` over a
> 2025–2045 grid for all eight planets. The risk it closes is the
> single-point-of-failure where every internal cross-check (Axis A's
> lamberthub/Kepler paths, the corrector, the gauntlet) consumes the *same*
> astropy states — a frame-convention, time-scale (TDB/UTC/TT), or
> jplephem-reader bug would shift them all consistently and pass every existing
> gate. An independent reader over identical data catches exactly that class.
> Confirmed agreement is sub-km / sub-mm/s (same data, different reader =
> numerical precision); a larger systematic would be the bug, reported, not
> tolerance-widened. This is consistency validation of the *generator's
> ephemeris layer*, not stored trajectory data — seeds-not-tracks is untouched.

#### 16.7.9 Schema v4.2 — segment center, TOF bounds, source ephemeris

The v4.2 sub-rev adds three **additive, optional, non-signature** fields. None
is a signature field, none changes the canonical signature, and no existing row
is required to carry any of them. Backfill is lazy and source-gated (as in
§16.7.6): the structure ships now; values are populated only as sources are read
and a published number is in hand.

1. **`trajectory.segments[].center`** (optional string) — the body the
   segment's conic/arc is centred on. Absent ⇒ `"Sun"` (heliocentric, the
   status quo). It is a free string, consistent with the top-level
   `orbit_elements.center` (which is likewise a free string defaulting to
   `"Sun"`); deliberately **no enum**, so planet-centric (moon-tour) cycler
   segments can name their primary. This is schema-readiness for tracked
   issue #76 (planet-centric / moon-tour segments).

2. **`trajectory.segments[].tof_days_bounds`** (optional `[min, max]`, days) —
   a *published* time-of-flight range, for sources that state a range rather
   than a point value. Exactly two numbers, both `> 0`, with `min <= max`.

   **Non-containment decision.** `tof_days_bounds` is deliberately **not**
   required to contain the segment's `tof_days` (when both are present). The
   motivating example: the Aldrin outbound segment carries `tof_days = 146`
   (the circular-coplanar idealization) while Rogers et al. 2012 Table 4
   publishes a STOUR real-ephemeris range of **161–172 d** for the same leg.
   These are different *model framings* of one physical leg, and both are
   sourced. Forcing the point value to lie inside the published range would
   reject valid, sourced data, so the validator checks shape and ordering only.

3. **`source_ephemeris`** (top-level, optional string) — the ephemeris model
   the source paper states its published numbers were computed against (e.g.
   `"DE405"`, `"DE430"`, `"STOUR ephemeris"`). When present it must be a
   non-empty string. This is the reduced, sourced essence of the earlier
   "SPK-kernel references" v4.2 candidate (§16.7.8 item 1): rather than store
   sampled kernels, we record *which* ephemeris a row's anchor numbers trace
   to, which is what cross-checking actually needs.

**Schema v4.2 fields:**

```yaml
source_ephemeris: "DE430"     # ephemeris the source's numbers were computed against
trajectory:
  segments:
    - id: out
      from: E
      to: M
      tof_days: 146           # circular-coplanar idealization (point value)
      tof_days_bounds: [161, 172]   # Rogers et al. 2012 STOUR range; NOT required to contain tof_days
      center: Sun             # absent => Sun; name the primary for planet-centric segments
```

**Validation split.** The JSON Schema (`data/catalogue.schema.json`) enforces
the structural shape it can express — `tof_days_bounds` is exactly two numbers
each with `exclusiveMinimum: 0`, `center` is a string. The Python semantic gate
(`validate_schema_invariants`) enforces what JSON Schema cannot: `min <= max`,
non-empty `source_ephemeris`, and (by *omission*) the non-containment rule
above. The permissive `additionalProperties: true` posture is preserved
throughout. No rows are backfilled in this rev (structure + validation only).

#### 16.7.10 Schema v4.3 evaluation (2026-06-05)

CCSDS 502.0-B-3 mining (`docs/notes/2026-06-05-ccsds-odm-502-mining.md`)
surfaced five v4.3 candidates. Each is gated on the YAGNI doctrine of §16.7.6 /
§16.7.9: a field ships only with both a demonstrable in-repo consumer and
sourced data to backfill at least one row today. Four are deferred; one is
adopted.

1. **Row supersession links — `supersedes` / `superseded_by` — ADOPT.** CCSDS
   carries `PREVIOUS_MESSAGE_ID` / `NEXT_MESSAGE_ID` (Table 6-2/6-3) to chain
   message revisions. The catalogue has a live need: the `vem-emeeve-3syn` row
   was premise-invalidated on 2026-06-05 (Jones AAS 17-577 found no feasible
   one-synodic EMEEVE cycler; the attested truth lives in the two-synodic
   members `jones-2017-vem-emevve-outbound` and `jones-2017-vem-meevem-inbound`).
   That supersession was expressed only in prose `notes` and `data_gaps[]`,
   where nothing can verify the named siblings still exist. The deciding
   consumer is **referential-integrity validation**: a first-class link lets
   `validate_schema_invariants` check (new Rule 6) that every target resolves to
   an existing row id and is not a self-reference — a machine-checkable guarantee
   prose cannot give, and a capability genuinely distinct from the grouping role
   of `family{}`. Both gates pass: the consumer is the cross-row validator built
   in this rev, and the sourced backfill is the `vem-emeeve-3syn` row, whose
   `superseded_by` now names its two attested replacements. Shipped as an
   additive, optional, non-signature array pair; the JSON Schema enforces
   array-of-non-empty-strings, the Python gate enforces resolution and
   no-self-link.

2. **`TRAJ_BASIS=SIMULATED` marker — DEFER (redundant).** CCSDS `TRAJ_BASIS`
   distinguishes as-flown (DETERMINED/TELEMETRY) from design-study (SIMULATED)
   trajectories. The catalogue already carries this distinction twice: every row
   is `source: literature` (none is as-flown telemetry) and `model_assumption`
   names the dynamical framing (`circular-coplanar` / `analytic-ephemeris`). A
   `traj_basis` field would encode no information not already derivable from
   those two, with no consumer asking for it. Trigger to revisit: ingestion of
   an *as-flown* mission kernel alongside a design-study row, where the two must
   be told apart at the field level.

3. **Per-element sigma uncertainty — DEFER (no sourced data).** JPL SBDB
   publishes per-element 1-σ; adopting σ on golden anchors was already noted as a
   future extension (§16.7.6). The gate is sourced data, and none exists: a sweep
   of the held mining notes finds only formula symbols (the `±` in the Endgame
   leveraging equations) and tour-accuracy percentages (`±5%`/`±10%` Tisserand-vs-
   CR3BP agreement), never a published 1-σ on a catalogue anchor value. No data ⇒
   defer. Trigger: a source that tabulates an explicit uncertainty on a number we
   store as a golden anchor.

4. **`CCSDS_OCM_VERS`-style version semantics — DEFER (already satisfied).**
   CCSDS uses a `major.minor` version string whose minor bumps for corrections
   and major for breaking changes. `catalogue.schema.json` already carries
   `"version": "<major>.<minor>"` as a string and the project applies exactly
   that bump discipline (v4 → v4.1 → v4.2 → v4.3). Nothing further is needed; no
   consumer reads structured version components.

5. **OCM exporter time-system caveat in `data/README.md` — ADOPT (doc only).**
   The message-global `TIME_SYSTEM` / `EPOCH_TZERO` constraint is recorded in
   §16.7.8 item 2 but was absent from the `data/README.md` OCM-export paragraph.
   Added there as a two-line caveat so the exporter contract is stated where the
   export is described. No schema or code change.

**Net for v4.3:** one schema field pair (`supersedes` / `superseded_by`) with a
cross-row referential-integrity rule and one sourced backfill row; one doc
caveat. Candidates 2–4 deferred with explicit trigger conditions. The schema
version is bumped 4.2 → 4.3; the canonical signature is unchanged (the links are
non-signature provenance metadata).

#### 16.7.11 Schema v4.4 — per-field provenance tags (2026-06-05)

Forge phase 0 task 3 back-fills the per-field provenance vocabulary defined by
`src/cyclerfinder/data/provenance.py` (the `SOURCE_REGISTRY` source keys and the
`Fidelity` tiers) onto the catalogue. Four additive, optional, non-signature
top-level row keys are added — `orbit_source`, `vinf_source`, `orbit_fidelity`,
`vinf_fidelity` — plus an optional declared `validation_tier`.

The deciding consumer (the YAGNI gate of §16.7.6/§16.7.9) is the already-shipped
forward-compatible `validate_provenance_tags`: until this rev the function was a
no-op because no row carried tags. The tags make the cross-source vs same-source
distinction, and the cross-fidelity refusal (the S1L1 5.65-vs-4.99 bug class),
machine-checkable rather than prose-only. Both gates pass: the consumer is live,
and the sourced backfill is mechanical (every tag is derived from metadata
already on the row — `model_assumption` → fidelity; an explicit Russell table
reference in the field's own `note:`/`source_quotes:` text → source key, with
author-token and `first_published.doi` fallbacks — by
`scripts/backfill_provenance_tags.py`, whose docstring records the rules). No
external source is consulted and no new physics value is introduced; this is
provenance metadata only.

Coverage: 224 of 237 rows are tagged; the 13 rows whose own citations name no
`SOURCE_REGISTRY` paper (Jones 2017, Sanchez Net 2022, Arenstorf, Genova,
Wittal, Hernandez, the two Russell-Strange moon-tour families, and two
broad-class family stubs) are left **untagged** — an absent tag is the explicit
"unknown" marker, and such a row classifies `UNVALIDATED` rather than carrying an
invented source key.

The JSON Schema constrains `*_source` to the `SOURCE_REGISTRY` enum and
`*_fidelity` to the `Fidelity` enum (enforced structure, hence the 4.3 → 4.4
bump); the Python semantic gate adds the rule it cannot express — a declared
`validation_tier`, when present, must equal the tier `classify_validation`
computes from the row's sources + fidelities, so a row cannot over-claim
`cross_validated`. `validation_tier` is normally left absent: the tier is a
*computed* classification, frozen as a monotone census ratchet in
`tests/data/test_validation_tier_census.py` (the validation-strength sibling of
the `cycler_class` census), not stored data — storing it would duplicate derived
state and risk drift. The canonical signature is unchanged (the tags are
non-signature provenance metadata). The live distribution at this rev is
`{cross_validated: 5, consistency_checked: 218, unvalidated: 14}`; the five
cross-validated rows each pair two genuinely different citations at one fidelity
(e.g. the Aldrin classic cycler: orbit from Rogers 2012 Table 1, V∞ from Russell
2004 Table 3.4).

#### 16.7.12 Schema v4.5 — validation_level (the §14 gauntlet level) (2026-06-06)

Forge phase 3 / task #104 adds one additive, optional, non-signature top-level
row key: `validation_level` (enum `V0`..`V5`), the §14 gauntlet level — the
highest gate a row has **mechanically** passed. The level is back-filled by
`scripts/backfill_validation_level.py` (idempotent, `--dry-run`), which applies
§14 to RECORDED in-repo test evidence, never aspirationally.

At this rev exactly one row earns above the floor: `aldrin-classic-em-k1-outbound`
is **V1** because its real-DE440 cycler clears §14 V1 — every leg re-solved with
`lamberthub` izzo2015 + gooding1990 agrees to < `V1_TOLERANCE_MPS`, AND the Kepler
forward re-propagation residual passes (the "re-propagated with the Kepler
propagator, planet positions met < tol" half of V1). This is demonstrated with
teeth by `tests/verify/test_agreement_lamberthub.py`
(`test_report_includes_lamberthub_path`,
`test_real_eph_paths_a_and_c_pass_b_flags_model_mismatch`). The rows the gauntlet
machinery has exercised at the internal-consistency floor are stamped **V0**: the
Aldrin INBOUND twin (no test builds/cross-checks it on real ephemeris) and the
other M6b regression rows (all `EXPECTED_SKIPS` — incomplete leg data or wrong
topology, so they do not pass real-closure). Every other row is left untagged —
an absent `validation_level` is the explicit V0 floor.

**V2 (2026-06-07, the §14 V2 class-split amendment).** Exactly one row reaches
V2: `aldrin-classic-em-k1-outbound` is **V2-powered**. Under the original single
V2 gate the powered Aldrin drifted ~4.14e8 km / 3 laps (≈2072× tolerance) — a
structural fail, because the cross-cycle rotating-frame-repeat metric is the
wrong instrument for a per-cycle-retargeted cycler. The amended V2-powered gate
(≥3 consecutive cycles, every planned encounter achieved with the per-cycle
maintenance applied **and** bounded intra-cycle drift vs the planned trajectory)
**passes**: over 3 consecutive in-family cycles the Mars-flyby V∞ continuity is
≤1e-6 km/s (clause a, ≤ the 0.5 km/s encounter tolerance), the intra-cycle
Kepler forward-reprop residual is ≤0.002 km (clause b, ≤ the 1 km bound), with a
strictly-positive in-family maintenance ΔV (≈2.76–2.91 km/s/cycle). Evidence:
`tests/verify/test_aldrin_v2_powered.py`. The Aldrin **inbound** twin stays V1:
its real-window optimiser lands on a ballistic ΔV≈0 neighbour rather than the
in-family powered solve (the recorded off-family resolver issue), so the
"per-cycle maintenance applied" half of V2-powered is not demonstrated for it.
The four #137 free-return rows (`russell-ch4-{5.30gGf3,9.94Gg3,5.75ggF3,
9.353Gg2}`) are **not** promoted to **V2-ballistic**: they are `cycler_class:
multi-arc` — the V1 evidence closes a single E→M→E free-return *ellipse slice*
(arc span ≈0.3–0.4 of the catalogue 4.27/6.41 yr period), but the full cycler
includes Earth-to-Earth resonant phasing intervals (e.g. the 3:2 full-rev
return) the single ellipse does not represent, so no continuous ≥3-lap trajectory
exists to propagate. They stay V1.

The JSON Schema constrains the value to `V0`..`V5` (hence the 4.4 → 4.5 bump); the
Python semantic gate `validate_validation_level` adds the over-claim rule it
cannot express — a row may declare a level above `V0` only when it appears in the
in-repo mechanical-evidence registry (`_LEVEL_EVIDENCE` in
`src/cyclerfinder/data/validate.py`). No row silently promotes itself off the
floor without a sourced, in-repo evidence pointer (golden discipline — when in
doubt, V0). The canonical signature is unchanged (the level is derived metadata).
This level mirrors the ledger's §13.8 `validation_level` field and the §13.6
spit-out → trust ladder: catalogue and ledger carry the same record the finder
emits.
