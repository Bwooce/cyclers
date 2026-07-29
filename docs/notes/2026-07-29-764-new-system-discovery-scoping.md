# Scoping `#764`: which system gets the resonant-manifold pipeline next — anchor search + recommendation

**Task:** `#764` (research/scoping only — no code, no catalogue changes, no runs), the first
concrete step of `#760`'s new-system discovery campaign after the user go-ahead ("go on 760").
For each candidate system — Saturn-Titan/Enceladus, Neptune-Triton, Pluto-Charon — determine:
(a) whether a published Anderson-&-Lo-2011-style validation anchor exists (digit-grade
eigenvalues / Jacobi constants / ICs / connection states the project's own corrector could gate
against); (b) whether the catalogue already carries relevant rows; (c) whether
`search/jovian_resonant_families.py` + `search/jovian_resonant_connections.py` are directly
repurposable (mu + seed geometry only) or need structural changes; then recommend ONE system and
a spec-complete first task, mirroring `#752`'s own format.

**Sources checked this pass:** the full `#752→#759` task-chain notes + both jovian modules read
in full (all 1307 + 721 lines); `docs/notes/CORPUS_INDEX.md` + `cyclers_pdf/papers/` filename
grep for Titan/Enceladus/Triton/Pluto/Charon; `data/catalogue.yaml` grep for all six body names;
`core/satellites.py` registry GMs; the `#730` acquisition-backlog master list; and a real web
literature search per system. One decisive primary source was downloaded and text-verified
first-hand this pass (Vaquero 2013, see below) — every load-bearing number quoted from it was
read from the extracted text layer, not from an abstract.

---

## Headline verdict: **Saturn-Titan first.** It is the only candidate with a genuine
digit-grade, Anderson-Lo-Table-1-style anchor — and the anchor is free to acquire.

| System | mu | Digit-grade family anchor? | Connection anchor? | Pipeline fit | Verdict |
|---|---|---|---|---|---|
| Saturn-Titan | ≈2.3658e-4 (thesis' own value) | **YES** — Vaquero 2013 Table 4.1: ICs + periods + unstable eigenvalues for 3:4, 6:5 resonant + L1/L2 Lyapunov orbits at C=3.010000 | Structural only (homoclinic + 3:4↔6:5 chain published as figures + a falsifiable family-termination claim, no state tables) | Direct repurpose (constants + seeds) | **GO — first target** |
| Neptune-Triton | 0.00020895 (paper's own value) | NO — families computed in print, but no IC/eigenvalue tables (one 6-digit C for one target orbit) | NO | Direct repurpose (retrograde is a non-issue in an isolated planar CR3BP) | Second choice; weak gate |
| Saturn-Enceladus | 1.9002485658670e-7 (in-corpus sourced) | NO for resonant families (published anchors are halo/DPO/DRO bifurcation structure) | NO | Corrector fine, but close-flyby resonant-instability regime is razor-thin (Hill radius ≈948 km vs 252 km moon radius) and real dynamics need Saturn J2 | Deprioritize |
| Pluto-Charon | 0.10876473603280369 (catalogue `#505`) | NO for unstable resonant families (Ross-RT 2026 Table I is digit-grade but for STABLE cyclers — already reproduced in-repo `#494`/`#505`) | NO | **Structural changes needed** (two-body seed construction breaks at mu≈0.109) | Do NOT attempt first |

---

## 1. Saturn-Titan — a real anchor, found and verified first-hand

**Anchor: Vaquero Escribano, T.M. (2013), "Spacecraft Transfer Trajectory Design Exploiting
Resonant Orbits in Multi-Body Environments," PhD dissertation, Purdue University (advisor
K.C. Howell), August 2013.** Freely downloadable:
`https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2013_Vaquero.pdf`
(35 MB, born-digital text layer, no OCR needed; md5 of the copy downloaded and read this pass:
`fdcbf871322b87cd1dd3448059cb2596`). Journal companion (the citable venue, same content class):
Vaquero & Howell, "Transfer Design Exploiting Resonant Orbits and Manifolds in the Saturn-Titan
System," *J. Spacecraft & Rockets* 50(5):1069-1085 (2013), DOI `10.2514/1.A32412`. The `#730`
backlog already lists the thesis in its §9 low-priority tail as "Saturn-Titan-Hyperion
resonant-manifold precedent" — this pass confirms it is much more than background: it is the
per-system gate source.

**Table 4.1 ("Initial State, Period, and Unstable Eigenvalue for Selected Periodic Orbits,
C = 3.010000"), read verbatim from the thesis text layer (§4.3.1, pp.~104-110):**

| Orbit | x (km) | ẏ (km/s) | T (days) | λu |
|---|---|---|---|---|
| 3:4 resonant | 1.25869e6 | 0.477301 | 66.3312 | **2,129.81** |
| 6:5 resonant | 1.14214e6 | 0.545759 | 71.2638 | **191.641** |
| L1 Lyapunov | 1.15897e6 | 0.447315 | 8.2829 | 1,004.72 |
| L2 Lyapunov | 1.25231e6 | 0.549329 | 79.7260 | 892.850 |

This is precisely the Anderson-Lo-Table-1 shape (in fact better on one axis: Anderson & Lo
published NO ICs for their families — Vaquero prints the seed x and ẏ directly, removing the
single largest source of `#753`/`#755`/`#756`'s multi-task search pain). The thesis states
µ ≈ 2.3658e-4 for Saturn-Titan (§4.4 text; consistent with this project's registry:
8978.14/3.7931207e7 = 2.3670e-4 — a GM-vintage delta of ~0.05%, same class as the 0.034%
Jupiter-Europa delta the Jovian module already documents and handles by using the SOURCE's own
value). The eigenvalues are dimensionless — a unit-ambiguity-free primary gate quantity; the
dimensional x/ẏ/T secondary gates need the thesis's own characteristic quantities
(l\*, t\* — the first task must extract them from the thesis's appendix/tables or derive them
and bound the vintage sensitivity explicitly).

**Connection-stage anchors (structural, not digit-grade — flagged honestly):** §4.3.1 also
publishes (as figures + prose, no state tables): a homoclinic connection of the 3:4 resonant
orbit (Fig. 4.9), a **"periodic resonant chain" that cycles indefinitely between the 3:4 and
6:5 resonances** (Fig. 4.10) — i.e. exactly the Anderson-Lo-style resonance-cycling periodic
orbit, in Saturn-Titan — a whole FAMILY of such chains continued in C (Fig. 4.11), and a
falsifiable published termination claim: "it is suspected that this family of periodic resonant
chains ends for a value of Jacobi constant C < 3.01400" (Fig. 4.12: the manifold gap at
C=3.014000). Also published: maneuver-free resonant→Lyapunov transfer scenarios (multiple
shooting), and a Hyperion application (Hyperion's 3:4 orbit reconverged in the CR3BP at
C=3.00937 is linearly STABLE; access from an unstable 3:5 resonant orbit at ~mm/s ΔV).

**`#759` tractability carry-forward (favorable):** the anchor families' instabilities are
λu=191.6 (6:5, comfortably tractable), 892.9/1004.7 (Lyapunovs, right at Jupiter-Europa
3:4-LO's λ≈1036 which the Newton machinery handled well), and 2,129.8 (3:4) — 2x the level that
worked, half the λ≈4445 that broke Newton shooting in `#759`. Family-stage risk: LOW.
Connection-stage risk: MODERATE (3:4's manifold is the required unstable leg for the chain);
the thesis's own published knobs (30 km manifold offset, 100,000 fixed points along the orbit,
50 TU integration) give a sourced fallback densification recipe if Newton struggles, exactly
the "paper's own denser interpolation method" mitigation `#759` recommended.

**Novelty ledger (per `feedback_literature_novelty_check_baseline`):** the Saturn-Titan
3:4↔6:5 resonance-cycling chain, homoclinic 3:4 connection, resonant↔Lyapunov transfers, and
Hyperion access are all PUBLISHED (Vaquero 2013; also Gawlik/Marsden-lineage "Titan Trajectory
Design Using Invariant Manifolds and Resonant Gravity Assists," AAS Spaceflight Mechanics 2010,
Caltech — additional Saturn-Titan manifold prior art). So the first task is deliberately
reproduction-shaped — that is what "validation anchor" means, mirroring how `#753` worked.
The discovery upside afterward: Newton-certified connection STATES (the thesis publishes
map-interpolation figures, no digits), other resonance pairs at other energies, and
catalogue-eligible cycler-class rows (the published chain family itself is catalogue-eligible
as a literature row — zero Saturn-Titan CR3BP rows exist today; the only Saturnian rows are
the Russell-Strange 2009 patched-conic family seed).

**Catalogue state:** `russell-strange-2009` Saturnian moon-cycler family seed (patched-conic,
Titan-Enceladus focus) is the ONLY Saturn row family — no CR3BP Saturn-Titan rows, no overlap
with this pipeline's object class.

## 2. Neptune-Triton — real recent literature, but no digit-grade gate

**Best sources found:** Miceli, Bosanac, Stuart & Alibay, "Motion Primitive Approach to
Spacecraft Trajectory Design in the Neptune-Triton System," AIAA SciTech 2024
(DOI `10.2514/6.2024-1280`; full text obtained and grepped this pass), and the journal version
Miceli & Bosanac, *J. Astronaut. Sci.* 73:11 (2026), DOI `10.1007/s40295-025-00545-z` (open
access). They state **µ = 0.00020895** (matches registry: 1428.49546/6.836527100580e6 =
2.0895e-4), compute planar resonant families 1:2, 1:3, 1:4, 1:5, 2:3, 3:5, 4:5, 3:4 prograde
and 3:1, 4:1 retrograde **using literally this pipeline's own two-body resonant-ellipse seed
construction**, and globalize manifolds off monodromy eigenvectors. But: **no IC, period, or
eigenvalue tables.** The only digit-grade numbers are the target 3:4 orbit's C_J = 1.75598 and
the arrival state's C_J,0 = 0.963141 (high-energy, open-Hill regime); their Table 1 lists
family-member Jacobi constants at 2-3 decimals only, and stability is described qualitatively
("in-plane stability indices ... close to 2" — i.e. mildly unstable, Newton-friendly if ever
pursued). No homoclinic/heteroclinic connections computed → connection-stage discovery there
would be novel, but the family-stage gate would be an order-of-magnitude/structural check, not
a digit gate.

**Retrograde subtlety (resolved — non-issue at this scoping level):** Triton's orbit is
retrograde relative to Neptune's spin, but an ISOLATED planar CR3BP is orientation-agnostic —
the frame simply co-rotates with Triton, and the published Neptune-Triton CR3BP work above does
exactly that. `#599`'s retrograde machinery matters for multi-moon (CCR4BP: Proteus orbits
opposite-sense to Triton) and inclined/real-ephemeris contexts, which is where the planar-CR3BP
fidelity caveat for this system genuinely bites — a writeback-tier caveat, not a pipeline
blocker. Catalogue state: zero Neptune-Triton rows (only Voyager-2's Grand-Tour mga_tour row
mentions Neptune).

## 3. Saturn-Enceladus — wrong physical regime for THIS pipeline, and no resonant anchor

In-corpus sources give a sourced µ = 1.9002485658670e-7 (Frauenfelder-Koh-Moreno 2023, SIADS
22:3284, digested `#744`) and Saturn-Enceladus periodic-orbit content — but it is all
libration-point/Halo/DPO/DRO bifurcation structure (also Moreno et al. 2024 bifurcation graphs,
digested `#728`, incl. the Appendix-A 29 km/14 km-altitude halo family), NOT p:q resonant
families. Fresh search found the same picture (Enceladus science-orbit papers add Saturn-J2 +
Enceladus-oblateness terms — a model extension this pipeline does not have); the
endgame/moon-tour literature reaches Enceladus via Titan leveraging and halos, not via unstable
Enceladus-resonant orbits. Physically this is not an accident: Enceladus's Hill radius is
≈ a(µ/3)^(1/3) ≈ 948 km against a 252 km moon radius, so the close-flyby mechanism that makes
these resonant families strongly unstable (Anderson & Lo p.177-178; `#755`/`#758`'s own
corroboration axis) operates only in a razor-thin shell above the surface. A resonant-manifold
campaign here would have no anchor AND a physically marginal object class. Deprioritize; the
Enceladus-relevant unlock in the corpus (halo bifurcation graphs) belongs to a different
pipeline.

## 4. Pluto-Charon — no anchor for unstable resonant families, and real structural changes

**No published unstable-resonant-orbit/manifold-connection source found** (searched: Zotos
orbit-classification and capture papers, PSS 2018/2019; Jbara 2025 arXiv:2510.13479 ZVC/chaos
study — no periodic-orbit tables; Giuliatti Winter et al. "Sailboat island" is a STABLE
first-kind family). The digit-grade Pluto-Charon-adjacent source the project already has —
Ross & Roberts-Tsoukkas 2026 (arXiv:2606.29189) Table I, reproduced in-repo `#494` and
instantiated at the true Pluto-Charon µ = 0.10876473603280369 as catalogue row
`ross-rt-pc-cycler-32-2026` (V2, `#505`) — is for STABLE prograde cyclers: no saddle
eigenstructure, no manifolds, wrong object class for this pipeline (and
`ResonantNode.from_candidate` correctly refuses non-real-saddle candidates).

**Structural changes genuinely required at µ ≈ 0.109** (the dispatch bullet's suspicion is
confirmed): `two_body_resonant_seed`'s construction (barycentric GM=1 two-body ellipse,
periapse at r=1 = "the secondary's radius") assumes µ << 1 — at µ=0.109 the secondary sits at
x = 1-µ = 0.891, the two-body-about-the-barycenter approximation is poor everywhere near the
pair, and the p:q resonance labeling itself blurs (no near-integrable Kepler limit). A credible
route exists — seed at small µ and continue in µ upward (the multiple-shooting +
mass-parameter-continuation "system translation" technique is published in Vaquero 2013 §4.5,
and `#494`/`#505` already demonstrated the corrector + C-sweep working at this µ for a stable
family) — but that is a new seed-strategy build with NO gate at the far end. Highest risk of
the three; do not attempt first.

## 5. What is reusable vs system-specific (from the full module reads)

**Already generic (takes a `system`/`mu` argument; zero changes needed):**
`two_body_resonant_seed` (pure p:q geometry), `two_body_flyby_rotation_seed` (µ is a
parameter), `converge_candidate`, `survey_candidates`, `basin_robustness_scan`, `_classify`
(Barden-authoritative + `_planar_floquet` cross-check), `europa_closest_approach` (misnamed —
it is generically "distance to the secondary at 1-µ"), the whole
`correct_symmetric_fixed_jacobi`/continuation/Barden core stack, and on the connection side
`ResonantNode.from_candidate`, `correct_connection`, `find_homoclinic`/`find_heteroclinic`,
`own_section_points`, `ydot_from_section_eq7` (section signs and ghost-guard radii are
parameters/defaults, not hardwired).

**Jupiter-Europa-specific (the thin layer a new system re-instantiates):** the sourced
constants (`ANDERSON_LO_MU`, `ANDERSON_LO_C_FLYBY`, `TABLE1_TARGETS`, Table-2/3 states),
`_EUROPA_SMA_KM` + `jupiter_europa_system()`, the empirically-located seed tables, the gate
tolerances (each justified against Anderson-Lo's own print precision — a NEW system needs its
own tolerance justification against ITS source's precision), and the paper-specific section
convention (`SECTION_X_SIGN=-1` is Anderson-Lo's own negative-x section; Vaquero's maps use
both x<0 and x>0 regions, so the Saturn-Titan module must pick and document its own per-target
section). Conclusion: a Saturn-Titan module is a THIN sibling of `jovian_resonant_families.py`
— constants + seeds + gate — exactly the relationship that module has to
`resonance_network.py`.

## 6. Recommendation: **Saturn-Titan.** Concrete first task (`#765`, dispatchable as written)

**Task `#765` — "Saturn-Titan 3:4/6:5 unstable resonant families + L1/L2 Lyapunov + Vaquero
Table-4.1 gate"** (spec-complete; Sonnet-tier per `[[feedback_subagent_model_tiering]]`, same
tier `#753` ran at):

0. **Acquire + register the anchor first**: download the Vaquero 2013 thesis PDF (URL + md5
   above) into `cyclers_pdf/papers/`, register in `CORPUS_INDEX.md` (born-digital text layer,
   no OCR), file its own digest todo per `[[feedback_per_paper_digest_todo]]`; record the JSR
   companion (DOI `10.2514/1.A32412`) as the citable venue. Re-verify every Table-4.1 number
   and the µ value against the PDF directly (this note's extraction is one hop removed — the
   task must not inherit them unverified, per `[[feedback_ground_citations_against_content]]`).
1. **New module `search/saturn_titan_resonant_families.py`** mirroring the Jovian module's
   structure: `saturn_titan_system()` at the THESIS's own µ (2.3658e-4 as printed; use more
   digits if the thesis appendix provides them), with l\*/t\* taken from the thesis's own
   characteristic quantities if printed (else registry Titan sma 1,221,870 km + GM sum, with
   the vintage delta measured and documented exactly as the Jovian module documents its 0.034%
   µ delta).
2. **Seeds**: the existing `two_body_resonant_seed(3,4)`/`(6,5)` sweep at C=3.010000 via the
   existing `survey_candidates`/`converge_candidate` — PLUS direct Newton seeds from Table
   4.1's own printed (x, ẏ) nondimensionalized (the `#758` lesson: paper-sourced seed windows
   beat blind grid scans; here the source prints the seeds, so lead with them, and
   `basin_robustness_scan` around them).
3. **Gate (sourced, dual-criterion, honest per-row report)**: primary gate on the dimensionless
   unstable eigenvalues λu ∈ {3:4: 2129.81, 6:5: 191.641} at a stated tolerance justified by
   the thesis's own ~6-significant-digit printing plus a MEASURED µ/l\*/t\* sensitivity bound
   (expect ~1e-4..1e-3 relative; state it, don't assume it); secondary criteria on T in days
   {66.3312, 71.2638} and the dimensional (x, ẏ) match (unit-conversion-dependent, so
   secondary); and the L1/L2 Lyapunov rows {1004.72, 892.850} as cheap positive controls run
   through the SAME classification path (they catch system-constant errors independently of
   the resonant search). A miss on any row is reported as a miss, `GateRow`-style — never
   silently dropped.
4. **Explicitly out of scope for `#765`** (its own Task-B analog, registered only after the
   families confirm): the connection stage — homoclinic Wu(3:4)∩Ws(3:4) and the published
   3:4↔6:5 periodic-resonant-chain family. Its gates are STRUCTURAL only (existence, Fig-4.10
   geometry, and the falsifiable published termination claim "family ends for C < 3.01400") —
   a weaker gate class than Anderson-Lo's Tables 2/3, flagged as such up front. λu≈2130 on the
   3:4 manifold leg is the known risk point (`#759` carry-forward); the thesis's own manifold
   parameters (30 km offset, 100k fixed points, 50 TU) are the sourced densification fallback.

**Honest counterweight:** like `#753`, `#765` is reproduction-shaped by design; the discovery
plays (Newton-certified connection states, other resonance pairs/energies, catalogue rows for
the published chain family and any new ones, all through `search/literature_check.py`) come
after the gate passes, and the Saturn-Titan prior-art surface (Vaquero-Howell, Gawlik/Marsden
lineage) is bigger than Jupiter-Europa's was — expect the novelty gate to bite harder here.
That is the correct trade against Neptune-Triton's near-absent gate and Pluto-Charon's
no-gate-plus-structural-work: validation first, discovery second, exactly the pattern that made
the Jovian chain honest.
