# #645 — Creative discovery-strategy pass (2026-07-18), the #605 sequel

Analysis-only (no code, no catalogue writes, no dispatches), mirroring `#605`/`#623`'s format.
Inputs actually read (not summarized from memory): `#605`'s and `#623`'s full bullets +
`docs/notes/2026-07-17-623-strategic-review.md`; the dedicated bullets for `#606`-`#620`
(seedless-corrector arc + SE↔EM negatives), `#624`/`#628`/`#642`/`#643` (generative-model arc,
cross-μ falsification and its Opus-adjudicated root cause), `#627`/`#629`/`#633`/`#638`/`#639`/
`#640` (RRT real-moon-μ arc), `#641` (Sun-Jupiter census), `#607`/`#609`/`#610`/`#625`
(mass-limited negatives + interval certification), `#635`-`#637` parking lot, `#644` (running);
`src/cyclerfinder/search/deflated_newton.py` and `variational_periodic_orbit.py` module headers;
`#570`/`#543`'s cycler-network scoping; CORPUS_INDEX spot-checks; plus four time-boxed web
searches (covariant-Lyapunov-vector computation; deflation-for-periodic-orbits precedent; recent
cycler literature; conditioned generative orbit models).

## The wall inventory this pass generates against (nothing below may silently re-hit one)

- W1 **Family-selection/basin wall** (`#605`'s diagnosis). Genuinely attacked by the
  `#606`-`#618` seedless spectral correctors — but `#606`'s own honest caveat stands: the
  seedless method has its OWN selection bias (slides onto the wrong, already-known family).
- W2 **Manifold-conditioning wall** (`#619`): at EM-L2 the one-period stroboscopic map amplifies
  ~2e4×, so a one-shot one-period STM eigendecomposition from a torus point known to ~1e-3
  cannot recover the UNSTABLE direction (tens of degrees of error); stable direction is robust.
  Named by `#619` as what blocks SE↔EM closure; predicts the same for any strongly-unstable
  libration point.
- W3 **Ghost minima** (`#620`/`#626`): integration-free collocation of open ARCS has spurious
  machine-zero-residual minima with O(1) physical loop defect; can VERIFY a near-real seed,
  cannot SEARCH; even `#619`'s real 166,000 km near-miss seed lands on ghosts.
- W4 **No μ-conditioning in the generative model** (`#642`, Opus-adjudicated): cross-μ transfer
  FALSIFIED — μ-mismatched Earth-Moon-geometry seeds collapse to L4/L5 equilibria under
  `correct_periodic` at foreign μ; structural, not small-N. In-distribution lift (~13-27×) real.
- W5 **Mass-limited physical-bend gates** (`#599`/`#607`/`#609`, several interval-certified via
  `#610`/`#625`): small moons cannot bend at the system's natural velocity scale. Physics.
- W6 **Fragmented near-C_L1 bifurcation structure at small μ** (`#627`/`#629`/`#633`): 0/16,375
  on-topology-and-stable at Saturn-Titan; the ν≈1.042 near-miss snaps topology under
  continuation; `#638` closed not-warranted on the strength of it.
- W7 **"Sweep another body" marginal value is low** (`#633`, `#641`): explicitly not re-proposed.
- W8 **Validation ceiling / new-input-gated** ([[project_validation_ceiling]]): further V-tier
  progress and much of discovery is gated on new external inputs, not more iteration.

---

## Ranked shortlist

### 1. Recover the EM-L2 unstable direction with segment-anchored discrete-QR (covariant-Lyapunov-vector) extraction, then re-run `#619`'s connection corrector
**What.** `#619` extracted manifold directions from a SINGLE one-period STM eigendecomposition
(forward, and a backward variant flagged "for study only"). The standard cure for exactly its
failure mechanism — "the STM is effectively linearized along a wildly-diverged trajectory" — has
existed since Benettin (1980) and is the basis of every modern Lyapunov-vector computation
(Ginelli et al. 2007; Dieci–Van Vleck discrete-QR): split the period into 10-20 short segments,
and because `#618`'s torus is an explicit functional representation, RE-PROJECT the base
trajectory onto the torus surface at every segment boundary, so the linearization stays along a
near-torus arc everywhere. Compose the segment STMs with QR re-orthonormalization; the dominant
subspace of the product gives the unstable direction at the departure point. Per-segment
amplification is chosen ~3-10× (not 2e4×), and base-point error is reset to O(torus residual)
each segment instead of compounding, so the direction estimate error should scale like
O(residual) — ~0.05° — rather than O(residual × 2e4) — the observed tens of degrees. Validate on
SE-L2 against trusted GMOS manifold data (the `#619` positive-control protocol, reusable
verbatim), re-run `#619`'s own perturbation-robustness test, then feed the new directions into
the EXISTING 12-unknown/18-residual `#538` corrector and re-run.
**Walls check.** This attacks W2 head-on rather than around it (what `#620` tried); W3 does not
apply (shooting, not arc collocation); W1 applies to the connection corrector as before, but the
`#619` floor (norm 0.855, 166,000 km) was reached WITH directions now known to be wrong by tens
of degrees — the seed quality changes materially. `#619` considered one-shot backward extraction
(a structurally different torus approximation, |dot|~0.46 even on clean SE-L2) — that is NOT
this; segment-anchored forward QR was never tried and is not covered by `#619`'s impossibility
argument, which is specifically about one-period one-shot linearization.
**Honest risks.** (a) The connection may genuinely not exist — three independent methods failed,
but each at a method-level obstruction (direction error / ghosts / ghosts-from-real-seed), none
excluding existence; (b) the torus's own 9.5e-4 surface error may dominate in a way segmenting
cannot fix (empirical question — `#619`'s robustness test answers it in a day); (c) even a
recovered direction only restores `#619`'s method to a fair trial, not a guaranteed closure.
**Cost.** ~2-4 days (extractor + SE-L2 positive control + robustness quantification + corrector
re-run). **Confidence.** MEDIUM the direction becomes recoverable to ≪1°; LOW-MEDIUM (~15-25%)
that closure follows. Highest information-per-cost on the program's most consequential live
negative: either the first genuine SE↔EM cislunar cycler (after mandatory Fable adjudication),
or the `#538`-`#626` negative upgraded from "conditioning-limited" to "fails even with correct
directions" — a materially stronger, more citable final state. Model: Opus.

### 2. Deflation in the seedless spectral corrector's coefficient space — an anti-basin family enumerator from two never-combined in-house capabilities
**What.** `deflated_newton.py` (`#524`) implements Farrell-Birkisson-Funke deflation as a
generic primitive but has only ever been aimed at shooting/scalar residuals — i.e., at
basin-restricted methods (its own docstring lists it among them). `#606`'s seedless spectral
corrector removed the seed-basin problem but has its OWN documented selection bias (its cold
start slid onto the vertical-Lyapunov family instead of the halo). Combining them — deflate
already-found solutions out of the SPECTRAL least-squares residual and re-solve from the same
cold start — turns the seedless corrector into a systematic enumerator of distinct families at
fixed (μ, C): exactly the missing answer to W1's residue. Design points (real but tractable):
the deflation distance must be gauge-invariant (min over time-shift phase — cheap via FFT
cross-correlation of coefficient vectors — and mod the discrete symmetries), and every deflated
find must pass the module's existing independent-Radau cross-check before being counted
(deflation finds more minima; the cross-check keeps them honest). Web check found no prior
astrodynamics-standard application of deflation to CR3BP family enumeration — fresh in the
field, not just in this codebase.
**Positive control.** Earth-Moon at a fixed C where the family census is known (JPL-verifiable):
one cold start + repeated deflation should re-enumerate the known distinct families.
**First real target.** The `#633` Titan-corridor near-miss: the ν≈1.042 (3,3)-cluster snapped to
different topologies under 1D continuation (a path-tracking failure mode); a spectral+deflation
census at that (μ, C) region has no path to snap and would settle whether `#633`'s
stability-negative hides thin-basin members — a method-independent check of a fresh headline
negative, in the spirit of [[feedback_bugfix_invalidates_past_searches]] without needing a bug.
**Walls check.** No manifolds (W2 n/a); W3's ghost pathology was specific to open arcs — closed
periodic loops in `#606`/`#611` validated cleanly against independent integration, and the Radau
gate stays mandatory; W6 is the thing being TESTED, not assumed away (honest expectation: `#633`'s
171 on-topology members were all wildly unstable — more members found by a better enumerator are
not more STABLE members; LOW odds the negative flips).
**Cost.** ~2-3 days. **Confidence.** MEDIUM-HIGH it works mechanically; LOW-MEDIUM it changes any
standing conclusion. Its durable value is the reusable capability: every future census/anchor
hunt gets a basin-bias-free enumeration mode. Model: Opus for design, Sonnet behind it.

### 3. Ingest the JPL SSD Three-Body Periodic Orbits catalog as a new external input: known-family gate + multi-μ goldens
**What.** `#641` discovered two things at once: (a) JPL's public periodic-orbit database
(`ssd-api.jpl.nasa.gov/doc/periodic_orbits.html` — family-labeled, multi-system, the modern
Doedel-tradition catalog) is exactly the authoritative reference this project keeps needing, and
(b) `search/literature_check.py` is structurally unfit for raw (non-cycler) periodic-orbit
candidates — it matched 5 physically-distinct Sun-Jupiter families to the same single cycler
paper (a documented, real tooling gap). Ingest the catalog into a local registry and build a
deterministic known-orbit matcher (nearest catalogued member in (μ, C, T, symmetry/family class,
amplitude) with an explicit distance threshold), plus promote selected rows to sourced golden
anchors at multiple μ (satisfies [[feedback_golden_tests_sourced_only]] — external published
values, not our own code's output).
**Why now.** `#644` (running) and any future census will hit the same adjudication gap `#641`
hit; this is the cheapest genuinely NEW external input available (W8 says new-input-gated), and
it is the load-bearing prerequisite for any cross-μ work regardless of idea 4's outcome.
**Walls check.** Mechanical; no dynamical walls. It discovers nothing by itself — stated plainly.
**Cost.** ~1-2 days (API fetch + registry + matcher + tests). **Confidence.** HIGH on feasibility
and utility as tooling. Model: Sonnet.

### 4. The μ-conditioning question (explicitly asked): architectural conditioning is data-gated and mostly moot; the one live cheap test is physics-based coordinate normalization
**Verdict on the architectural gap.** You cannot learn a μ-conditional from a corpus with zero
variance in μ — the `#608` corpus is single-μ (Earth-Moon), so "add μ-conditioning to the model"
is not an architecture decision, it is a DATA acquisition problem. The only realistic multi-μ
labeled corpus is idea 3's JPL catalog — but once that catalog is in hand, nearest-neighbor +
the existing continuation machinery from catalogued members at the nearest μ almost certainly
beats any generative model as a cross-μ seeder (the model's only edge is density between
catalogued families, which `#614` showed it does not track tightly). Honest conclusion:
**cross-μ generative modeling is probably a permanently dead end for this program** — not
because a conditioned model couldn't be trained, but because its niche is occupied by a simpler,
stronger classical tool the moment the data needed to train it exists.
**The one cheap live test.** A COORDINATE fix needs no new data: re-express the genome in
μ-adapted coordinates — position relative to the secondary/nearest libration point scaled by the
Hill radius (μ/3)^(1/3), energy as ρ=(C−3)/(C_L1(μ)−3) (this project's OWN `#629` design read
found ρ nearly μ-invariant across the sourced RRT anchors — direct in-house evidence for this
normalization), period in local characteristic time — and decode at the target μ by inverting
with the target's scales. This targets W4's verified mechanism directly (seeds land in the RIGHT
geometry relative to the target μ's libration points instead of Earth-Moon absolute
coordinates, which is what fed the L4/L5 equilibrium collapse). Pilot = re-run `#624`'s exact
protocol (with `#642`'s equilibrium filter) at μ=0.001 + Sun-Earth with the transform inserted.
**Honest limits.** Hill scaling is only valid in the secondary's dynamical neighborhood; the
corpus mixture is unknown, and global/resonant/L4-L5-adjacent families will not transform
sensibly (may need a cheap geometric classifier to transform only the near-secondary subset).
`#628`'s cluster-reweighting negative does not cover this (that was output reweighting, not a
coordinate change).
**Cost.** Hours-to-a-day. **Confidence.** LOW-MEDIUM (~25-40%) that a real transferred lift
appears; the pilot is decisive either way and, if negative, closes the entire cross-μ question
permanently (with idea 3 as the standing replacement). Model: Sonnet.

### 5. Inter-cycler transfer-compatibility network over the 361-row catalogue — a genuinely new object class from the project's unique asset
**What.** The M5/M6-tier work `#570`'s scoping explicitly deferred ("a genuine inter-cycler
transfer cost... is new search work, a separate task, not something to fake into this schema").
Nobody — here or in the literature — has computed which of THIS catalogue's 361 rows are
mutually transfer-compatible: for every pair sharing an encounter body, gate on V∞
magnitude/Tisserand compatibility (the same machinery class as `#604`'s compatibility check) and
phase-window feasibility, producing a graph over the catalogue (edges = feasible low-ΔV
cycler-to-cycler taxi transfers, hubs, connected components). Output: populated
`cycler_networks.yaml` derived rows (needs `#570`'s still-open schema or a simplified first
slice) + a citable network-structure artifact. A cheap multi-cycler ITINERARY that no single
published cycler provides would be a genuinely novel object class.
**Walls check.** None of W1-W6 apply (uses catalogued encounter data + closed-form gates). The
real risk is different and named: phase alignment between independent cyclers is measure-zero
without powered patching, so most edges will carry real ΔV, and the honest outcome may be "no
cheap edges exist" — itself a citable, registry-worthy negative about cycler-network feasibility
(this program's stated main product).
**Cost.** ~1 week including the `#570` schema slice. **Confidence.** MEDIUM-HIGH a real artifact
results; LOW a headline novel itinerary emerges. Model: Sonnet behind a Fable design read on the
between-cycler cost definition (the one real judgment call).

---

## Checked and NOT recommended (with reasons — do not re-propose without new evidence)

- **Rectangular/over-collocation ghost suppression** for the arc corrector: `#626` established it
  only converts false convergence into honest non-convergence at ~1e6 km defects; no seed in the
  true basin exists. Revisit ONLY if idea 1 produces a genuinely better near-miss.
- **Any new body/system sweep with existing methods**: W7; `#633`/`#641` just measured this
  pattern's marginal value. Includes Venus-inclusive/Crocco-style rosters unless a new METHOD
  motivates them.
- **Bigger/deeper generative models or heuristic family-tag conditioning** on the existing
  corpus: `#614` clean negative; premise unchanged.
- **W-Z computer-assisted-proof build** (`#636`): stays parked per its own bullet (user-decision
  gated). Noted honestly: the SE↔EM triple-method negative is the closest thing to a forcing
  function yet, and if idea 1 restores correct directions and closure STILL fails, that negative
  becomes both stronger and more certification-worthy — the natural trigger to revisit `#636`.
- **Low-thrust/solar-sail cycler classes as catalogue objects**: `#519` closed (Sims-Flanagan
  infeasible); corpus already covers the low-thrust cycler literature (Rauwolf/Chen-Longuski/
  Pony Express); powered/sail members are design-dependent controls, not discoverable invariant
  objects — a poor fit to the V-gauntlet's validation philosophy. Not a discovery direction.
- **p:q Earth-Moon resonant-cycler census** (Acta Astronautica 2020 infrastructure paper found
  in this pass's web check): mission-design USE of known E-M cycler classes, territory the
  catalogue/corpus already covers; at most a minor corpus-completeness check, not a direction.
- **QBCP refit for other systems / bigger grids / #520 revival**: all rejected with reasons in
  `#623`; nothing has changed.

## Recommended dispatch order
1. Idea 1 (segment-QR unstable-direction recovery + `#619` re-run) — Opus, ~2-4 days.
2. Idea 3 (JPL SSD catalog ingestion + known-family gate) — Sonnet, ~1-2 days; unblocks `#644`
   follow-ups and any future census adjudication.
3. Idea 2 (deflation × spectral enumerator) — Opus design + Sonnet build, ~2-3 days.
4. Idea 4 (Hill/ρ-normalized cross-μ pilot) — Sonnet, hours; decisive close-out of the cross-μ
   question either way.
5. Idea 5 (catalogue transfer-compatibility network) — Fable design read first, then Sonnet.

Task numbers deliberately NOT allocated here (next-unused stays `#646`); the coordinating
session registers/dispatches per [[project_task_numbering_convention]].

Sources (web checks this pass): Ginelli et al. covariant-Lyapunov-vector computation
(arxiv.org/abs/1105.5228; link.springer.com/article/10.1007/s00332-012-9126-5), Farrell-
Birkisson-Funke deflation (SIAM J. Sci. Comput. 37(4), 2015 — no astrodynamics family-enumeration
application found), JPL SSD periodic-orbit API (ssd-api.jpl.nasa.gov/doc/periodic_orbits.html),
Litteri et al. CMDA 138:25 (link.springer.com/article/10.1007/s10569-026-10299-x), p:q resonant
cycler infrastructure (sciencedirect.com/science/article/abs/pii/S0094576520300916).
