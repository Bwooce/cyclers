# #686 — Third fresh discovery-strategy pass: N≥4-body discovery (2026-07-22)

Analysis-only (no code, no catalogue writes, no dispatches), mirroring `#661`/`#679`'s format,
but — per the task's own mandate — leading with an honest TRACTABILITY verdict before any
shortlist. The question, in the user's own framing: the project has thoroughly mined 3-body
(CR3BP/ER3BP) discovery; is a genuine N≥4-body discovery search tractable at all within this
project's realistic time/compute budget, and if so, how do we "attack that impossible problem"?
Explicitly NOT in scope (barred by the dispatch): `#679`'s rejected BCR4BP/HR4BP
robustness-continuation item, and the `#533`-`#620` SE-L2↔EM-L2 torus-connection approach.

Inputs actually read (not summarized from memory): `#686`'s and `#679`'s full bullets; the full
`#533`/`#537`/`#538`/`#544` arc bullets including the `#537` DO-NOT-CERTIFY stamp; `#611`/`#612`
(seedless spectral wall-crossings); `#617`/`#618`/`#619`/`#620` in full (the arc's terminal
diagnoses); `#607`/`#609` (hierarchical/small-body clean negatives); `#522`'s full bullet
(linking-number build + the `#534`/`#548` shelve stamps); `src/cyclerfinder/core/bcr4bp.py`
(all 430 lines) and `core/qbcp.py` (structure + coefficient tables); the `#500` Keplerian-map
genome verdict note in full; the in-corpus Kumar–Anderson–de la Llave SIADS digest in full;
`search/literature_check.py`'s Jovian anchors (Liang CGE, Hernandez IEG, Russell-Strange
Ganymede-Europa/Ganymede-Io); the catalogue's `hernandez-2017-jovian-ieg-triple-family` row and
primary-body census (381 rows: 32 Earth-primary planet-moon + heliocentric lineage, 21 Saturn,
18 Jupiter, 1 Pluto, 4 P1-generic); `data/empty_regions.jsonl`'s Jovian entries
(`jovian-IEG-vilm-2026-06-09`, `jovian-perm-vilm`, `#433` Galilean quasi-cycler,
`#318`/`#501` real-ephemeris Sobol EGE/GCG/EGCE/CGCEC shoots); plus eight time-boxed web
searches and four abstract-level fetches (arXiv 2109.14815, arXiv 2309.06073, Bhanu Kumar's
publication page, Kumar-related 2024-2026 follow-up work).

## 1. What "N≥4-body" can even mean here — the taxonomy, ruled in/out

(a) **Coherent restricted 4-body models of the Sun-Earth-Moon system (BCR4BP/QBCP/HR4BP).**
Fully built in-repo (`core/bcr4bp.py` #292, `core/qbcp.py` #533/#592) and fully exercised: the
`#533`-`#620` arc took this machinery to its genuine capability ceiling (see §2). Synodic-
resonant periodic families in these models (Brown-Peterson-Henry-Scheeres SIADS 2024) are
μ_sun-continuations of CR3BP orbits — perturbation-check territory, rejected by `#679` and
barred here. RULED OUT (exhausted or barred).

(b) **Literal general N-body periodic solutions (choreographies, figure-eight, free-fall
orbits).** Real mathematics, zero spacecraft relevance: every known choreography-class solution
requires comparable (usually equal) masses and has no solar-system realization for a restricted
test particle to exploit. No positive control with mission meaning exists. RULED OUT.

(c) **Equilateral/Trojan restricted 4-body problems** (Sun-Jupiter-Trojan + spacecraft;
Lucy-adjacent). Genuinely 4-body-native families exist in the mathematical literature, but
Trojan GMs are ~1e-9-1e-7 of Jupiter's — any repeated-encounter object is bend-gate-dead by
the same structural mass wall `#607`/`#609`/`#433`-Amalthea established as model-robust
(`max_bend` is a pure function of (μ, r_p, V∞); no formulation fidelity can raise it). And
the non-encounter structures (tadpole/horseshoe librations) are parking, not cycling. RULED OUT.

(d) **Hierarchical compositions** (`#609` cycler-of-cyclers; `#607` multi-moon small bodies).
Both closed as clean negatives on mass grounds, and the composition idea remains blocked on
missing partners: the catalogue census (re-checked live) has NO heliocentric cycler to any
giant planet (heliocentric rows are Earth-Mars/Venus-Earth lineage only), and Mars has no
gate-passing moon-cycler (`#609` step 1). Revival condition, named for the record: an
Earth-Jupiter or Earth-Saturn heliocentric cycler entering the catalogue. RULED OUT for now
(blocked on data, not method).

(e) **Periodically-forced 3-body systems where the 4th body is a second MOON of the same
primary — the concentric circular restricted 4-body problem (CCR4BP).** Jupiter + Europa +
Ganymede + spacecraft, with the moons taken at their 2:1 Laplace commensurability so the model
is TIME-PERIODIC and the whole stroboscopic-map machinery applies. This is the one branch of
the taxonomy that is (i) genuinely 4-body-native (see §3), (ii) NOT touched by any prior task
in this project (grep: "CCR4BP"/"concentric" appears nowhere in `src/` and only incidentally in
one 2004-dissertation mining note), and (iii) carried by a published, reproducible 2021-2023
positive-control chain (Kumar-Anderson-de la Llave). RULED IN — the single tractable lane. §3.

(f) **Quasi-periodically-forced / full-ephemeris discovery.** Two incommensurate forcing
frequencies (real moon periods, Sun + Moon simultaneously, DE440) kill the stroboscopic
reduction, push invariant tori up a dimension, and have NO published discovery-scale success
anywhere in the literature surveyed. This is where "impossible problem" genuinely starts.
Discovery stays in reduced models; full ephemeris remains what it already is here — the V4
per-row VALIDATION lane. RULED OUT for discovery, by design not by defeat.

## 2. Tractability verdict (first, per the mandate)

**The general N≥4-body discovery problem remains intractable for this project — but the reason
is sharper than "chaos" or "dimensionality," and the project's own `#533`-`#620` wall does NOT
generalize to everything 4-body. There is exactly one bounded, genuinely 4-body-native lane
that clears every diagnosed obstruction: the Jupiter-Europa-Ganymede CCR4BP.**

The honest decomposition:

**(a) Dimensionality.** Adding a 4th body generically breaks the Jacobi integral and adds the
forcing phase as a dimension. That is survivable exactly when the forcing is TIME-PERIODIC:
the stroboscopic map reduces the extended 5D phase space (planar case) back to a 4D map, tori
become 1D invariant circles, their whiskers 2D — and two 2D manifolds in a 4D map space
intersect GENERICALLY in points. Heteroclinic connections are therefore findable objects, not
measure-zero accidents (this is Kumar et al.'s own dimension-count, and it is why their search
works). With incommensurate forcing (case (f) above) this reduction dies and the problem
genuinely earns the word "impossible" at this project's scale. So: tractability is not about
N; it is about whether the model reduces to a periodically-forced PCRTBP.

**(b) Chaos / conditioning — what the `#533`-`#620` wall actually was.** The arc's terminal
diagnoses are quantitative and specific, not generic pessimism: EM-L2's monodromy spectral
radius ~1540 (`#612`), one-stroboscopic-period amplification ~2-3e4 in the QBCP (`#619`),
achievable torus invariance residual ~9.5e-4 (`#618`) — hence the unstable-manifold DIRECTION
is unknowable to tens of degrees by post-hoc one-period STM extraction (`#619`'s decisive
finding), and the integration-free arc-collocation alternative drowns in ghost minima
(`#620`). Root drivers: (i) the Sun is an enormous perturbation of cislunar space (μ_S =
328,900.5 in Earth+Moon=1 units — an O(1) structural forcing, not a small parameter); (ii) the
target region (EM-L1/L2) is violently unstable; (iii) the SE↔EM geometry spans a 388:1 scale
gap in a single coherent frame; and (iv) — decisively for tractability — the target object
(an SE-L2↔EM-L2 cycler) has NEVER been published by anyone: there was no positive control at
the end of that chain, only at its intermediate stages.

**(c) Does the wall generalize?** Apply those four drivers as a test to any proposed 4-body
lane. For the JEG CCR4BP: (i) forcing is WEAK — the perturbing moon's mass ratio is 2.5e-5
(Europa) / 7.8e-5 (Ganymede), a genuinely perturbative parameter nine orders below the Sun's
cislunar forcing, so tori persist robustly and continuation from the PCRTBP is honest; (ii)
the anchor objects are mean-motion-resonant orbits in the chaotic sea, whose per-period
instability is mild compared to EM-L2's ~1540 (and this is measurable cheaply on day one —
if it were EM-L2-like, Kumar et al.'s published computations would not exist); (iii) no scale
gap — Europa and Ganymede orbit the same primary at a 1.6:1 radius ratio in one frame; (iv)
published positive controls exist at EVERY stage: tori (AAS 21-651 / arXiv 2109.14815),
whiskers + connections in the periodically-perturbed setting (SIADS 24(1):219-258, 2025 —
already in corpus, digested), and an actual Ganymede-4:3 → Europa-3:4 CCR4BP transfer (Acta
Astronautica 211:76-87, 2023). Additionally, the specific `#619` failure MECHANISM (post-hoc
STM eigenvector extraction from an imperfect torus) has a published structural fix in exactly
this literature: the parameterization method solves torus + stable/unstable BUNDLES
SIMULTANEOUSLY in one invariance equation (SIADS paper Eq. 4.4, per the in-repo digest), so
the whisker directions come out of the solve with their own residual instead of through one
violently-amplified propagation. The wall does not transfer on any of its four axes.

**Verdict, stated plainly:** a genuine N≥4-body discovery search is NOT tractable in general —
not in coherent Sun-Earth-Moon models (exhausted, `#533`-`#620`), not in general N-body
(physically unrealized), not under incommensurate forcing (no reduction, no controls) — but it
IS tractable in exactly one configuration, the Laplace-commensurate Jupiter-Europa-Ganymede
CCR4BP, which is also (not coincidentally) the only solar-system configuration combining a real
moon-moon commensurability with perturbing masses large enough to matter. The attack on the
"impossible problem" is to refuse the general problem and take the one reduction the solar
system actually offers.

## 3. The one shortlist item

### 1. Laplace-locked repeating resonant tour (ballistic quasi-cycler) in the Jupiter-Europa-Ganymede CCR4BP

**What.** Search for a closed heteroclinic CYCLE among whiskered resonant tori of the planar
CCR4BP — e.g. Jupiter-Ganymede 4:3 torus → (unstable manifold ∩ stable manifold) →
Jupiter-Europa 3:4 torus → back to the Ganymede 4:3 torus — i.e. a repeating, low-energy,
ballistic (or small-deterministic-dv) itinerary that revisits both moons' resonant
neighbourhoods indefinitely, phase-locked to the 2:1 forcing. The published record stops at
ONE-WAY transfers (Acta Astronautica 2023); the closed cycle — the cycler-class object — is
unclaimed (checked live against Bhanu Kumar's own publication list through 2025: one-way
transfers, the AAS 23-397 secondary-resonance study, an Oberon PCRTBP survey, cislunar MMR
work; no cycle, no repeating tour).

**Why genuinely 4-body-native (the task's core bar).** Three independent grounds. (1) The
anchor objects are not perturbed 3-body orbits in any useful search sense: in the CCR4BP the
PCRTBP resonant periodic orbits cease to exist as periodic orbits at all and become 2D
whiskered tori in a 5D extended phase space with no Jacobi constant — the connection search
lives in a phase space with different dimension counts than any CR3BP search this project has
run. (2) Stronger: AAS 23-397 (Kumar-Anderson-de la Llave 2023, abstract fetched live) shows
the 4th body generates SECONDARY resonances (11/34, 12/37, 23/71, 25/77 lockings between the
4:3 orbits' internal frequency and Europa's forcing) whose islands and overlap "cause a
complete structural change of the higher-energy unstable 4:3 orbits" — new invariant objects
with NO 3-body analog whatsoever, born of the 4th body, and their overlap is itself a
4-body-native transport mechanism. (3) The cycle's closure condition is a phasing
commensurability between TWO moons' geometry — a constraint that is meaningless in any 3-body
model. Nothing here could be found by perturbing an existing 3-body cycler, because no single
CR3BP contains an object that alternates which secondary governs it.

**Why this is not re-proposing already-swept ground (checked against the registry, not
assumed).** This project's Jovian coverage is extensive but entirely in OTHER regimes:
patched-conic flyby/VILM sweeps (`jovian-IEG-vilm`, `jovian-perm-vilm` — empty above a 6 km/s
V∞ floor with ~20 km/s gaps), the `#433` Galilean quasi-cycler sweep (best closure
Callisto-Ganymede-Callisto 8.3 m/s but literature-KNOWN; all fresh closures Amalthea-dead),
and the `#318`/`#501` real-ephemeris Lambert-seeded Sobol shoots (EGE/GCG/EGCE/CGCEC — all
"compute-bounded empty," 0 converged from top-k shots). Every one of those `empty_regions`
verdicts is method-conditional by the registry's own doctrine, and none had manifold-guided
low-energy capability: the CCR4BP resonant-torus lane is precisely the "new method whose
capability subsumes the prior sweep's" that the negative-results-registry doctrine names as
the legitimate re-entry condition. The catalogue's existing Jovian rows (Hernandez IEG,
Russell-Strange Ganymede-Europa/Ganymede-Io, Liang CGE — all anchored in
`literature_check.py`) are patched-conic flyby cyclers, a different species in a different
energy regime.

**Why it dodges the `#533`-`#620` wall.** The four-axis test of §2(c): weak forcing
(μ ~ 2.5e-5-7.8e-5 vs the Sun's O(1) cislunar forcing), mildly-unstable anchor objects (vs
EM-L2's ~1540-2e4), no cross-scale frame gap (1.6:1 vs 388:1), and published positive controls
at every stage (vs a never-published end target). Plus the structural fix for `#619`'s
specific mechanism: whiskers via the simultaneous bundle solve, not post-hoc STM extraction.
Honest caveat: mild instability is an expectation from the published record, not yet a
measured in-repo number — Stage B's first gate (below) measures it before anything expensive
is built on it.

**Named positive controls (staged).**
- PC0 (in-repo, free): the `#500` Keplerian-map goldens — 17 committed tests reproducing
  Ross-Scheeres 2007 and Grover-Ross 2009 (Jupiter-Callisto/Ganymede resonant-kick dynamics).
- PC1 (in-corpus, digested): reproduce the SIADS 24(1):219-258 whiskered-torus + connection
  results in the periodically-perturbed Jupiter-Europa PERTBP (the paper's own validation
  system; PDF filed at `kumar-anderson-delallave-2025-...-arxiv-2109.14814.pdf`, digest notes
  its validation tables are in the later pages, "read on demand" — that demand is now).
- PC2 (acquisition needed): reproduce the Acta Astronautica 211:76-87 (2023) Ganymede-4:3 →
  Europa-3:4 CCR4BP transfer. NOT in corpus. Open-access routes exist: arXiv 2109.14815 (the
  AAS 21-651 predecessor, 20pp), NTRS 20230005667, and the Georgia Tech repository copy;
  the AAS 23-397 secondary-resonance paper (arXiv 2309.06073) should be acquired in the same
  pass. File + digest + CORPUS_INDEX per corpus policy before building against them.

**Staged plan and cost (honest).**
- **Stage A — composed Keplerian-map screen (~2-4 days, no user GO needed).** Compose the
  existing positive-controlled `genome/keplerian_map.py` into an alternating two-map system
  (Europa map ⊗ Ganymede map — the "P3BA patching" the `#500` verdict itself named as the
  map's future use) and search for periodic itineraries of the composed map: cheap existence
  evidence + seed geometry (which resonance pairs, which phasings) before any expensive build.
  Honest caveat: each map ignores the other moon, so this is screen-grade heuristics only —
  a negative here would NOT be registry-grade; a positive is a seed, not a result. Note
  `#604`'s prior negative on Keplerian-map chaining was about INTERPLANETARY arrival V∞
  compatibility (Jupiter-arrival ~order of magnitude above the map's regime) and does not
  touch this intra-Jovian use, which is the map's native validated regime.
- **Stage B — the real build (~3-5 weeks total, USER GO REQUIRED — largest single build since
  the Track-A frontier).** (1) Planar CCR4BP EOM+STM module patterned on `core/bcr4bp.py`
  (exact 2:1 commensurability → time-periodic; structural tests: reduces exactly to
  Jupiter-Europa PCRTBP at μ_Gan→0 and Jupiter-Ganymede PCRTBP at μ_Eur→0, mirroring
  bcr4bp's μ_S→0 CR3BP test). (2) Adapt `#617`'s pseudospectral torus corrector — the
  problem shape is IDENTICAL (non-autonomous periodically-forced, θ1 locked to the forcing
  clock, θ2 free; 4D planar state instead of 6D) — this is the legitimate "EOM swap" case,
  unlike the CR3BP→QBCP jump `#611` documented as structural. (3) Whiskers: FIRST measure the
  one-period amplification of the target resonant tori and run `#619`'s own
  perturbation-robustness diagnostic on cheap post-hoc STM extraction (check, don't guess);
  build the Kumar bundle-solve only if that gate fails — this is the main cost fork
  (±1-2 weeks). (4) Manifold globalization + CPU mesh-intersection search in the 4D
  stroboscopic space (Kumar et al. used GPUs for a 5-7x speedup; CPU + days-long detached
  runs is fine per standing practice). (5) Cycle closure: mirror a found connection through
  the model's symmetry, solve the phasing commensurability, close the loop with the existing
  multiarc-BVP discipline; mandatory independent-integrator cross-check, Fable adversarial
  review, and V4 real-ephemeris continuation attempt before any catalogue write ("it closed!"
  is the danger signal).

**Novelty ceiling / honest risks.** MODERATE-HIGH, honestly bounded. The body pair has
published cyclers in a different regime (patched-conic flyby — Russell-Strange 2009 has a
literal Ganymede-Europa cycler anchor in our own `literature_check.py`), so a found cycle is a
new SPECIES on a published body pair, not a new body pair — analogous to how the prograde EM
cyclers coexist with classical ones. The unclaimed-cycle window is real but at genuine scoop
risk: Kumar/Anderson own this machinery and system and could publish a cycle at any
conference cycle. Risks: (i) the cycle may need small dv at patch points → a `quasi_cycler`
row with a dv_band, not a pure cycler — acceptable, that class exists; (ii) the CCR4BP is an
idealized commensurate model (true Europa:Ganymede ratio 2.0149, not 2:1) — the found object
is model-native, epoch-locked to an idealization, and the V4 real-ephemeris continuation may
honestly fail (`#501`'s real-eph EGE empties show Jovian real-eph closure is hard); even
then, a validated CCR4BP cycle + a characterized real-eph failure is a publishable,
registry-grade result on a genuinely new region×method cell; (iii) planar first pass only;
(iv) a clean negative (connections exist one-way, no phasing-compatible cycle) is a fully
acceptable final outcome and extends the anti-catalogue with a genuinely new method stamp.

**No second item is offered.** Every other lane surveyed fails the tractability test or an
exclusion, per §1 and §4 — a 1-item shortlist is the honest size.

## 4. Considered and explicitly rejected (with reasons — do not re-surface without new evidence)

- **BCR4BP/HR4BP synodic-resonant continuation** — barred by the dispatch; `#679`'s rejection
  stands (V4 DE440 gauntlet answers solar survivability per-row more strongly).
- **Any revisit of the SE-L2↔EM-L2 cross-system connection** (any technique) — `#619`'s
  conditioning bound and `#620`'s ghost-minima obstruction are quantitative and unaddressed by
  anything in this survey; the end target remains publication-free (no positive control).
  Note the parameterization-method bundle solve WOULD address `#619`'s extraction mechanism,
  but not the ~2e4 amplification-vs-residual arithmetic underlying it, and `#620`'s
  independent obstruction stands; the burden of proof for reopening is not met.
- **General N-body choreographies / free-fall periodic orbits** — no solar-system realization,
  no restricted-problem relevance, no mission-meaningful positive control. Pure mathematics.
- **Trojan / equilateral R4BP cyclers** — structurally bend-gate-dead (Trojan GMs), same
  model-robust mass wall as `#607`/`#609`/`#433`-Amalthea; libration structures aren't
  catalogue classes.
- **`#609` cycler-of-cyclers revival** — still blocked on missing partners at both ends
  (no giant-planet heliocentric cycler row; no Mars moon-cycler). Revival condition named in
  §1(d).
- **EM L4/L5 BCR4BP dynamical substitutes (Kordylewski-adjacent stability islands)** —
  genuinely 4-body-native (the Sun destabilizes the CR3BP triangular points and moves the
  stable structure), but parking, not transport; no repeating-encounter object class.
- **Saturn CCR4BP editions (Enceladus-Dione 2:1, Titan-Hyperion 4:3, Mimas-Tethys)** — the
  same machinery applies in principle (real commensurabilities exist), but the perturbing
  masses are 1-2+ orders smaller (Enceladus GM 7.2, Dione 73, Hyperion ~0.4 km³/s² vs
  Ganymede 9888, Europa 3203): kick amplitudes and resonance widths shrink with μ, and the
  Saturnian VILM registry entries already show the flyby-regime emptiness. JEG is the unique
  system with both a commensurability AND mass. Named as the natural SECOND target only if
  the JEG lane validates, not before.
- **Uranian CCR4BP (home-turf system)** — no present-day MMRs among the major Uranian moons →
  no commensurate forcing → no time-periodic reduction; the quasi-periodic-forcing
  generalization is over the tractability boundary (§1(f)). Kumar-Anderson's 2024 Oberon
  survey is PCRTBP-level for the same reason.
- **Direct seedless discovery in full-ephemeris N-body** — no invariant structure to anchor a
  search; find-in-model → V4-validate remains the correct architecture (validation-ceiling
  doctrine).
- **Low-thrust multi-moon tours** (2026 Saturnian low-thrust tour literature surfaced in the
  survey) — control-dependent design, not discoverable invariant objects (`#519`/`#645`
  doctrine).
- **Lunar-flyby-augmented Earth-Mars cycler establishment (Sun-Earth-Moon-Mars)** — design
  optimization on existing rows (precursor_mga territory, Rogers 2012 lineage), not a new
  object class.
- **Keplerian-map interplanetary chaining** — `#604`'s V∞-compatibility negative stands;
  Stage A's intra-Jovian composition is the map's native regime and is not that idea.

## 5. User decision points (flagged, not assumed)

1. **GO/NO-GO on Stage B** (~3-5 wk, the largest single build proposed since the Track-A
   frontier; cost is a column not a verdict per standing policy, but this is explicitly a
   portfolio call — the W-Z arc precedent shows multi-task corridors can absorb weeks with
   zero rows). Stage A + PC1 can run first without commitment and would sharpen the estimate.
2. **Schema**: `model_assumption` currently enumerates `circular-coplanar |
   analytic-ephemeris | cr3bp` — a CCR4BP object needs a new value (and the epoch-locking
   semantics of "locked to an idealized 2:1 commensurability" need a `validity_window`
   convention). Analogous to `#664`'s still-open SET-row question.
3. **Acquisition approval**: Acta Astronautica 211:76-87 author-manuscript routes + arXiv
   2109.14815 + arXiv 2309.06073 → corpus policy filing (OCR/digest/index) before PC2.
4. **Catalogue-row policy for model-native objects**: if a CCR4BP cycle validates in-model but
   fails V4 real-ephemeris continuation, does it get a row (new model_assumption value) or
   registry-only treatment? Decide before the campaign, not after ("hold writeback till
   confirmed" discipline).

## 6. Recommended dispatch order

Stage A (cheap screen; Opus for the composed-map design judgment, Sonnet for the sweep) and
the PC1 reproduction (Opus — numerical-methods tier, same class as `#611`/`#612`/`#617`) can
both start now, in either order or parallel. Stage B only on decision point 1's GO, informed
by Stage A/PC1. If the user declines Stage B, the honest position recorded here stands on its
own: the N≥4-body frontier has exactly one tractable lane, it is now named, scoped, and
positive-control-mapped, and everything else is exhausted, unphysical, or over the
quasi-periodic-forcing boundary.
