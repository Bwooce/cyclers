# #679 — Second fresh discovery-strategy pass (2026-07-22), downstream of #661

Analysis-only (no code, no catalogue writes, no dispatches), mirroring `#605`/`#645`/`#661`'s
format. This pass exists because `#661`'s entire 5-item shortlist is exhausted (`#662`
formulation-blocked, `#664`/`#665`/`#666`/`#667` all executed) and the concurrent W-Z
proof-machinery arc (`#636`, `#668`-`#678`) is pure capability-building against the
already-published Oterma benchmark — zero new catalogue rows even at full success. The job:
surface something that WOULD produce new candidate rows, without re-proposing anything from
`#661`'s list even in modified form.

Inputs actually read (not summarized from memory): `#661`/`#662`/`#663`/`#664`/`#665`/`#666`/
`#667`'s full bullets; `#605`/`#623`/`#645`'s full bullets + the `#645` note; the `#563`/`#564`/
`#565`/`#566`/`#569` symmetric-closure-family arc (including `#564` §3's asymmetric-closure
DEFER and `#565`'s Fable confirmation of it); `#607`'s clean-negative diagnosis; the CURRENT
STATE dashboard; a module-level capability map of `src/cyclerfinder/{core,search,nbody,verify,
data}` (129 search modules; read-only survey agent, spot-verified directly for every
load-bearing claim below); all 92 `data/empty_regions.jsonl` entries; catalogue census by
`orbit_class`/`model_assumption`/system; `docs/notes/CORPUS_INDEX.md` spot-checks against every
externally-surfaced paper below; `scripts/gmat_v4_uranus_generate.py`'s actual force model;
`search/lobe_overlap_scorer.py` and `search/cislunar_bct_search.py` module docstrings; plus
seven time-boxed web searches (BCR4BP/HR4BP periodic-orbit families 2024-2026; Earth-Moon
cycler literature 2023-2026; ER3BP families; low-thrust cyclers; LCS/superhighway transport;
lobe dynamics; ballistic capture at Mars; Phobos/MMX QSO families; recent arXiv cycler work).

## Survey headline: the obvious "new" leads are already in the building

The live web survey's four most promising-looking leads all turned out to be already absorbed —
worth recording so the next pass doesn't re-surface them:

- **Ross/Roberts-Tsoukkas 2025/2026 stable prograde E-M cyclers** (AAS 25-621, arXiv
  2606.29189): in corpus, mined, Table-I goldens committed
  (`data/golden/ross_rt_2026_cycler_families.yaml`), μ-continuation / Floquet-branch / 3D-lift
  campaigns all executed and registry-stamped (`#389`/`#392`/`#393`, `#434`/`#438`/`#444` —
  the 3D-lift novelty axis is formally CLOSED).
- **Braik-Ross 2026 "Orbital Networks in the Three-Body Problem"** (arXiv 2605.31543,
  reachable-set network over EM periodic-orbit families): in corpus, mined + KNOWN_CORPUS; the
  catalogue's EM cycler nodes derive from this lineage, and `#650` built this project's own
  transfer network.
- **Hiraiwa-Bando lobe-dynamics transfer design** (PRR 6:L022046 2024; Acta Astronautica 248,
  arXiv 2602.17444): in corpus, KNOWN_CORPUS, and `search/lobe_overlap_scorer.py` (#278)
  already implements the 2026 paper's flux-weighted-graph framework.
- **BCR4BP/HR4BP synodic-resonant continuation** (Brown-Peterson-Henry-Scheeres, SIADS 2024,
  arXiv 2402.19181): `core/bcr4bp.py` (#292) + `core/qbcp.py` + a μ_sun-continuation reach
  spike (#412) already exist, and the V4 DE440 real-ephemeris gauntlet already answers
  solar-perturbation survivability per-row more strongly than any bicircular model would
  (`#389`: far-amplitude EM cycler branch destroyed by solar tide, all 100 epochs; low-amplitude
  branches characterized in `#392`/`#393`).

The genuinely available items below are therefore mostly *internal* frontiers this project's own
recent results opened, not imports.

## Ranked shortlist

### 1. Asymmetric-closure census at Uranus — reopen the `#564` §3 deferral on `#663`'s new evidence

**What.** Enumerate the genuinely-asymmetric (free `rel_offset`) exact closures of the
`#558`/`#563`-lineage direct-construction system, starting Ariel-Umbriel and Titania-Oberon,
via multi-start reduced-formulation Newton/LM + deflation (`search/deflated_newton.py`) over
the (β, tof) box per (pair, direction, n_rev), with `#663`'s interval machinery certifying
empty sub-boxes for a bounded-exhaustiveness statement.

**Why genuinely new — and why it needs an explicit user GO.** The symmetric census (`#563`, 30
members, provably exhaustive *for the symmetric class only*) is the ONLY method lineage that
ever produced this project's one confirmed-novel finding (`#312`). `#562`'s own data showed
asymmetric near-closures (Oberon-Titania n=3, rel=114.15°, residual 7.9e-3; n=2, rel=268.19°,
3.5e-2) that the symmetric enumeration by construction cannot contain. `#564` §3 deferred the
asymmetric search (NO-GO: expensive 2D adaptive search, no exhaustiveness guarantee, highest
lit-risk basin, gated on Canales/Kumar novelty) and `#565` (Fable) confirmed the deferral —
so proposing this is a deliberate, evidence-based REOPENING of a standing dual-adjudicated
NO-GO, and the user must approve that reversal explicitly. The evidence that changed
(2026-07-19, `#663`, all post-dating the deferral): (a) an EXACT non-symmetric root exists in
**Ariel-Umbriel** at β≈74.3° (50-dps Newton, residual → ~1e-29) — asymmetric members are real
and NOT confined to the high-lit-risk Titania-Oberon basin; (b) the grid-gap risk was
*realized* (a genuine exact closure found between symmetric grid points near `#600`'s
near-miss, independently DOP853-confirmed); (c) the cost/rigor objections are materially
reduced: `#663` built the well-conditioned reduced formulations (cond(J)~3e3-5e3 away from the
symmetric degeneracy), the interval-safe Lambert/Kepler machinery, and the LM+DOP853
verification chain; (d) the deferral's stated novelty-gating precondition (Canales/Kumar) was
discharged by `#566`/`#569` (family written to catalogue 2026-07-11).

**Positive controls (in-repo, no acquisition needed).** Re-find `#663`'s β≈74.3° Ariel-Umbriel
root blind from multi-start; recover the committed symmetric goldens as the β∈{0°,180°}
degenerate cases; hit `#562`'s two recorded T-O near-closures as first seeded targets.

**Novelty ceiling / honest risk.** HIGH conditional on gates: any gate-passing asymmetric
closure is a new family member of the `#312` type in a symmetry class nobody (including us) has
enumerated. Realistic risk: both of `#663`'s exact finds are gate-hostile (the `#600`-adjacent
closure fails the `#324` bend gate at ~0.83°; the β≈74.3° root's gate status is unchecked), so
the realistic outcome is a census where most members fail physical gates — but the symmetric
census yielded 30 passers, so the prior is not negligible. Exhaustiveness is honest-bounded,
not proven: `#663` §4d documents the interval certificate's near-root width/depth failure mode.

**Cost.** 1-2 weeks. Opus for the formulation/well-posedness layer; Sonnet mechanical sweep
behind the existing deterministic gates; standard Opus+Fable adjudication for any survivor.

### 2. Interplanetary WSB repeating-capture quasi-cycler at Mars

**What.** Extend `#378`'s cislunar ballistic-capture-transfer chain search (`core/wsb.py`,
`genome/bct_transfer.py`, `search/cislunar_bct_search.py`) from Sun-Earth-Moon to
Sun-Earth-Mars: search for a repeating capture↔escape chain whose return leg re-acquires the
MARS weak-stability-boundary set each cycle — a `quasi_cycler`-class object whose Mars
"encounters" are temporary captures (weeks-months dwell) rather than hyperbolic flybys.

**Why genuinely new.** The repeating-capture object class has only ever been searched cislunar
(`cislunar-bct-wsb-quasicycler-2026-06-26`, clean negative). The Mars edition is a different
dynamical regime (the capture sets live in Sun-Mars space, arrival modulated by the Earth-Mars
synodic cycle) and has real mission pull (long Mars dwell with no insertion burn). No published
"ballistic-capture cycler" surfaced in the survey — the capture literature (Topputo-Belbruno
2015; the 2023 "ballistic capture corridors at Mars" co-orbital time-varying-manifold work) is
one-shot transfers only.

**Positive control (named, acquisition needed).** Reproduce Topputo & Belbruno, "Earth-Mars
transfers with ballistic capture," CeMDA 121:329 (2015) — NOT currently in corpus (only the
Belbruno 2004 textbook is); acquire + digest + index per corpus policy before building.

**Novelty ceiling / honest risk.** HIGH if anything closes; risk also high — the cislunar
analog was negative, Belbruno 2004 Thm 3.58 (capture on W is chaotic) cuts against clean
periodicity, and real Sun-Earth-Mars geometry makes any chain epoch-locked at best. Expected
outcome is an honest, registry-stamped negative that genuinely extends the anti-catalogue to a
new region×method cell.

**Cost.** 2-3 weeks including acquisition/digest and the un-genericized parts of the BCT
machinery (the W-set predicate is currently cislunar-specific).

### 3. Quasi-periodic "cycler corridor" census around the stable prograde EM cyclers

**What.** For each linearly-STABLE member of the Braik-Ross/Ross-RT EM cycler families (C21:
107/201 stable; C32: 164/201 stable at z0=0.24, per the `#438` registry entries; plus the
stable planar goldens), compute the surrounding KAM/quasi-periodic torus corridor with the
existing torus machinery (`search/variational_qp_torus.py`, `qp_torus_fixed_jacobi_
continuation.py`) and measure its extent — the measure-positive volume of naturally cycling,
station-keeping-free trajectories around each cycler.

**Why genuinely new (for this project and mostly for the literature).** This is exactly
`#444`'s own named redirect (b): "cycler-USABILITY in regions where the dynamical family is
known but its transport utility is not characterized." The RT/Braik-Ross papers (2025-2026)
establish the periodic orbits and their stability but publish no corridor-width/torus
characterization; stability ⇒ tori exist, but nobody has measured how much cycling volume
surrounds each member. Output: corridor extents per member + possible `quasi_cycler` rows.

**Positive control.** Reproduce a published EM CR3BP quasi-periodic torus family result
(Olikara-Scheeres 2012 GMOS EM families — the standard benchmark; in-repo `#612` torus
validation controls as the mechanical cross-check).

**Novelty ceiling.** LOW-MODERATE, stated honestly: this characterizes known objects rather
than discovering a species; but it is cheap, safe, directly mission-relevant, and the most
likely of the four to produce catalogue-adjacent output (new quasi_cycler rows or corridor
fields) this month.

**Cost.** ~1 week; heaviest reuse of any item here.

### 4. Periapse Poincaré-map cartography for repeating temporary-capture itineraries (Saturn-Titan first)

**What.** Build the one classical seedless discovery map this codebase genuinely lacks
(survey-verified: no periapsis/apsis Poincaré-map machinery anywhere in `search/`): periapse
maps in the planet-moon CR3BP (Davis-Howell lineage), used to find capture-lobe →
escape-lobe → re-capture repeating itineraries geometrically — the same object class as item 2
but found by exact lobe geometry rather than a W-set predicate, and at a moon (Titan) rather
than a planet.

**Why new & honest overlap note.** Method-shape absent from the inventory; but its discovery
content overlaps `#664`'s set-oriented pipeline (the statistical sibling of the same transport
question). Periapse maps give per-trajectory exact structure with no Monte-Carlo noise; `#664`
gives measure/residence statistics. Ranked last for exactly that overlap.

**Positive control.** Reproduce a published periapse-map structure: Davis & Howell 2012 (JGCD
35(1)) periapse-map transit/capture regions, or Villac & Scheeres 2003 (Hill problem escape
lobes) — either is a concrete, figure-level reproduction target. Acquisition needed (neither in
corpus).

**Novelty ceiling.** MODERATE-LOW. **Cost.** 1-2 weeks.

## Standing follow-on deliberately NOT claimed by this pass

`#664`'s own bullet leaves its natural next step explicitly on the table: applying the
validated GAIO pipeline to one of this project's own systems at equal rigor (new
section/energy/region-indicator glue + calibration). That is arguably the highest-expected-value
single dispatch currently available, but it belongs to `#661`-item-2's lineage — proposing it
here would be re-proposing an executed `#661` item in modified form, which this pass is
explicitly barred from. Flagged for the coordinator's dispatch decision on its own merits, with
its own open schema question (how a metastable SET gets a catalogue row) still unresolved.

## Considered and explicitly rejected (with reasons — do not re-surface without new evidence)

- **BCR4BP/HR4BP synodic-resonant continuation of the EM cyclers** (Brown et al. 2024
  pipeline): V4 DE440 gauntlet already answers solar-perturbation survivability per-row more
  strongly (`#389`/`#392`/`#393`); BCR4BP core + continuation spike exist (#292/#412); coherent-
  model periodic counterparts add elegance, not rows.
- **Lobe-dynamics turnstile design**: already built (#278, `lobe_overlap_scorer.py` implements
  Hiraiwa et al. 2026); paper KNOWN_CORPUS.
- **Reachable-set orbital-network mining**: Braik-Ross 2026 already mined; `#650` built the
  network; `#434`/`#438`/`#444` mined the nodes and formally closed the 3D-lift axis.
- **Small-body J2/C22 re-sweep of `#607`'s negatives** (the `#665`-SRP pattern-match — primary
  oblateness at Sylvia/Kleopatra is order 0.02-0.18, a genuinely dominant term the model
  omits): rejected after verifying the metric semantics ([[feedback_verify_metric_semantics_
  before_ranking]] discipline): `#607`'s 0/97,664 is 100% bend-gate-bound (0.1-8° achievable vs
  the ≥5° floor, a pure function of the MOONS' tiny GMs) — primary-field fidelity cannot flip a
  flyby-bend wall. The negative is model-robust, not model-conditional.
- **Planetary-J2 axis for the Uranian family**: already carried — `gmat_v4_uranus_generate.py`
  propagates degree-2 Uranus gravity (J2=3.34343e-3, Jacobson 2014); the validated members are
  J2-validated at the tier that matters.
- **Koopman spectral methods**: re-affirm `#661`'s rejection; `#664` delivered the
  operator-theoretic discovery value.
- **Diffusion/transformer generative seeding**: `#614` clean negative at n=54k; `#623`'s
  explicit rejection stands; no new evidence.
- **FLI/LCS "superhighway" cartography** (Todorović-Wu-Rosengren 2020): a triage/seeding layer
  that feeds the same corrector pipelines already walled by W1-W3; `ftle_scorer.py` partially
  covers it; no new object class.
- **Aerogravity assist / tethers / sail-class-β cyclers**: control/design-dependent, not
  discoverable invariant objects (`#519` closed, `#645` re-affirmed); a sail-class-β re-sweep
  is `#665`'s own registry-documented future extension (user-gated), not a new idea this pass
  may claim.
- **New-body-pair sweeps** (Jupiter-Saturn cyclers, Earth-Venus extension, Haumea, Pluto small
  moons, Janus-Epimetheus): `#645`'s W7 rejection stands; the small-GM systems are structurally
  bend-gate-dead per `#607`/`#667`'s lesson.
- **Phobos/MMX QSO census**: rich published MMX-era QSO-family literature (Hill problem with
  ellipsoidal Phobos, out-of-plane bifurcations, tori) makes this reproduction-only — the exact
  profile `#666` already demonstrated; the useful residue is a threshold decision (below).
- **Retrograde 1:-1 co-orbitals** (Morais-Namouni lineage, Ka'epaoka'awela): same
  weak-encounter-utility + "known object, relabeled" profile `#661` ranked last and `#666`
  confirmed; also risks reading as `#666`-in-modified-form.
- **ER3BP isolated-family push**: already executed (`#440`-`#442` gap analysis, isolated-seed
  modules, 5 registry entries 2026-06-24/25).

## User decision points (flagged, not assumed)

1. **GO/NO-GO on reopening the `#564` §3 asymmetric-closure deferral** (item 1). This reverses
   a standing dual-adjudicated NO-GO; the case rests on `#663`'s post-deferral evidence. Not
   dispatchable without an explicit user decision.
2. **Quasi-cycler admissibility for capture-chain objects** (item 2): minimum captures per
   cycle, dv_band ceiling for deterministic correction per cycle, and the periodicity tolerance
   that separates "quasi-cycler" from "sequence of transfers." Schema-shaped, analogous to
   `#664`'s still-open SET-row question.
3. **Torus-corridor schema** (item 3): new `quasi_cycler` rows vs. corridor-width fields on the
   parent cycler rows.
4. **`#667`'s 0.3-Hill-fraction "close" ceiling**: at small-Hill bodies it excludes the
   beyond-Hill QSO/DRO regime by construction (Phobos: the ceiling sits inside the body).
   Revisit or document as intended scope?
5. **Dispatch priority of the un-claimed `#664` own-system application** vs. this pass's items.

## Recommended dispatch order

**1 → 3 → 2 → 4**, with item 1 contingent on decision point 1 (if declined, promote 3 to
immediate dispatch — it is the cheapest and most likely to touch the catalogue). Item 4 only
if the `#664` own-system follow-on is not chosen instead (they compete for the same
transport-discovery niche).
