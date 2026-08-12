# `#836`: adjudication of `#810`'s Pluto-Charon (5,1) stable member — **ADMITTED at V1**

**Task:** `#836`, registered 2026-08-11 during `#810`. Opus adjudication + catalogue-writeback
decision for the stable Pluto-Charon (5,1) prograde CR3BP cycler member found by `#810`'s
fixed-hc C-sweep. Scope per the `#659`/`#817` adjudication precedent: (a) independent
reproduction beyond `#810`'s own in-task checks, (b) the narrow Charon clearance as an explicit
adjudication axis, (c) novelty framing under a live literature gate, (d) tier + writeback or a
reasoned refusal.

**VERDICT: ADMITTED** to `data/catalogue.yaml` as `pc-cycler-51-2026`, `validation_level: V1`,
`our_status: known-class-member`, **NOT claimed novel**.

The clearance axis, which was the one thing that could have flipped this, is resolved
decisively in the member's favour — and for a structural reason, not a lucky sampling. Three
corrections to the surrounding record are recorded below; none of them touch the physics.

---

## (a) Independent reproduction — every number reproduces

Rebuilt from scratch against the library primitives in a standalone script that deliberately
does **not** import `scripts/run_810_pc_51_fixed_hc_sweep.py`, so the reconstruction shares no
control flow with the run under adjudication. `#810`'s own recorded bracket was used only as a
brentq bracket and as linear continuation guesses; the root itself was re-found.

**Mandatory positive control first.** `#504`'s `sweep_32_positive_control()`, UNMODIFIED,
re-run in this session: `stable_found=True`, C=3.5795150197296888, x0=-0.6931982870433999,
T=11.833462517008149 TU, nu=-1.2e-07, `topology_ok`, `crosscheck_ok`. Against the committed
`ross-rt-pc-cycler-32-2026` row: dC=6.2e-13, dx0=3.1e-14, dT=2.6e-11. The anchor holds.

**The (5,1) member, independently re-derived:**

| quantity | `#836` (this pass) | `#810` (recorded) | agreement |
|---|---|---|---|
| C (brentq nu=0 root) | **3.167935964707404** | 3.167935964707279 | 1.2e-13 |
| x0 | **-0.7058054139667359** | -0.7058054139668293 | 9.3e-14 |
| ydot0 | -0.6722667009874522 | -0.6722667009872901 | 1.6e-13 |
| T (TU) | **24.305715846918527** | 24.305715846921732 | 3.2e-12 |
| Barden nu | 5.29e-11 | 5.98e-10 | both ~0 |
| winding | (5,1), (w1,w2)=(5.000000,1.000000), prograde, reaches_secondary | same | exact |
| Radau crosscheck | PASS, dJ=7.945e-13 | PASS, dJ=7.8e-13 | same |
| stable window nu=+1 | C=3.167773861779837 | 3.167773862 | ~1e-10 |
| stable window nu=-1 | C=3.168099843905549 | 3.168099844 | ~1e-10 |
| window width | 3.259821e-04 | 3.26e-04 | same |

Corrector crossing residual 1.6e-12. Four **perturbed restarts** (x0 guess offset by ±1e-3 and
±3e-3, period guess off by 0.2%) all reconverge to the same x0 within 3.1e-14 — the member is a
genuine isolated fixed point of the corrector, not an echo of the starting value.

`#656`'s seed independently re-verified with `hc=None`: reconverges at residual 5.4e-13 with
winding exactly (5.0000, 1.0000), prograde. `#810`'s correction of the `#656` bullet's prose
"hc=4" (that was the LOST (4,0) branch's index; the seed's own is hc=5) is confirmed.

**An algorithm-independent stability check `#810` did not do.** Barden's index comes from a
half-period STM and a G-matrix factorization. I also computed the **full-period monodromy
directly** by variational propagation and took its eigenvalues:

```
 3.022819172514786e-07 + 0.9999999999993328j    |lam|=0.999999999999378  arg=+1.570796024513
 3.022819172514786e-07 - 0.9999999999993328j    |lam|=0.999999999999378  arg=-1.570796024513
 0.999999995277145     + 9.720985107100692e-05j |lam|=1.000000000002023  arg=+0.000097209851
 0.999999995277145     - 9.720985107100692e-05j |lam|=1.000000000002023  arg=-0.000097209851
 det(M4) = 1.0000000000028044   (symplectic => 1)
```

All four eigenvalues sit on the unit circle (linearly stable), and the **nontrivial pair sits at
arg = ±pi/2 to seven digits**, so nu = cos(pi/2) = 0 *structurally* — an independent confirmation
of the nu=0 midpoint by a different algorithm. It agrees with the Barden index to 3e-7, which is
full-period integration noise (the trivial pair's own |lam| is off unity by 2e-12).

**Self-correction, recorded per house style.** My first version of this check selected the
"nontrivial" eigenvalue by maximising `| |lam| - 1 |`. That is the correct selector for an
*unstable* orbit and exactly the wrong one for a stable one — with every eigenvalue on the unit
circle it returned the trivial pair and printed a meaningless `nu=1.000000`. Caught because a
maximally-stable orbit reporting nu=1 is self-contradictory; fixed by selecting on `|arg|`. (A
second, trivial slip: my first run used `SweepResult.c_stable`/`x0_stable`, which do not exist —
the fields are `jacobi_mid`/`x0_mid`/`period_mid`/`nu_mid`.)

## (b) The narrow Charon clearance — the adjudication axis, and it is not a sampling artifact

`min_body_clearance_km` takes a plain `min()` over `solve_ivp` output points with
`max_step = period/4000`. That is a **sampled minimum**, i.e. only an upper bound on the true
closest approach, and at Charon-encounter speeds the spacing between samples corresponds to
hundreds of km of arc — enough, in principle, to hide the whole 139 km margin. So the recorded
745.5 km could not be taken at face value.

**Refined it three ways, and it holds exactly.** Brent local minimization on dense output,
across two integrators (DOP853, Radau) and rtol/atol from 1e-11 to 1e-13, with grids from 4001
to 200001 points:

| body | grid min (km) | Brent-refined min (km) | t* | local minima |
|---|---|---|---|---|
| Charon | 745.5244 | **745.5244** | 0.500000 T | 5 |
| Pluto | 5128.1485 | **5128.1485** | 0.304350 T | 8 |

Total spread of the Charon figure across all refinements: **1.9e-9 km**.

The reason it is exact is structural, and worth stating plainly: **the closest approach to
Charon occurs at t = T/2**, the orbit's second perpendicular x-axis crossing
(|y| = 1.6e-11, xdot = 0 identically) — and Charon sits *on the x-axis* at (1-mu, 0). So the
closest-approach distance is not an integrated minimum at all, it is the single coordinate
difference

```
|x(T/2) - (1-mu)| = 0.0380369589206897 nd = 745.5243948455181 km
```

By the symmetry of the orbit this is also a true periapsis with respect to Charon (the velocity
at T/2 is purely along +y, perpendicular to the Charon radius vector, which lies along x). The
sampling concern is fully retired — not because the sampling got lucky, but because the quantity
is a symmetry-point coordinate.

**Against the 606.0 km radius (Nimmo 2017, the `#660` gate's sourced value): altitude 139.52 km,
1.230 Charon radii from centre. PASSES.**

**The margin cannot be bought back with stability.** Across the *entire* |nu|<1 window the
clearance barely moves:

| C | nu | closest approach (km) | altitude (km) |
|---|---|---|---|
| 3.167773861779837 (nu=+1 edge) | +0.999992 | 746.492 | 140.492 |
| 3.16785 | +0.528966 | 746.038 | 140.038 |
| **3.167935964707404 (member)** | **0.000000** | **745.524** | **139.524** |
| 3.168 | -0.392056 | 745.141 | 139.141 |
| 3.168099843905549 (nu=-1 edge) | -1.000000 | 744.543 | 138.543 |

So there is no "safer" member of the stable window to prefer — the whole window flies at
138.5-140.5 km. The trade does not exist, which is itself a finding.

**Weighing "passes the gate" against "admissible".** These are not automatically the same
question, and the `#659` Antiope precedent is the reason to ask. But the two cases are
categorically different regimes, not two points on one scale: `#659` REJECTED candidates that
passed **30-38 km BELOW** the secondary's surface — those trajectories do not exist, they
intersect the body. This member passes **139.5 km ABOVE** the surface. A trajectory that clears
by 1.23 body radii is a real trajectory; the question it raises is *operational margin*, not
physical existence. The honest resolution is therefore **admit and record the thinness
prominently**, which is what the row does — in `validation_level`'s own evidence text, in the
`_LEVEL_EVIDENCE` entry, and in a dedicated NARROW CHARON CLEARANCE paragraph in `notes`.

The caveats that genuinely attach (recorded on the row, none of them blockers): this is an ideal
**point-mass** CR3BP figure; Charon is not spherical; the CR3BP assumes a circular secondary
orbit where the true PC eccentricity is ~0.0022; and there is no ephemeris or solar
perturbation. A real mission would care about all four at a 139 km margin.

Also verified, since k2=1 is what makes a single period bound the whole encounter geometry:
**exactly one** Charon close approach per period. The five local minima of the Charon distance
over one period are at 10739 / 14319 / **745.5** / 14319 / 10739 km — one genuine encounter,
four far passes, consistent with the (5,1) winding.

**External physical cross-check.** C(L1) = 3.62101825776712 (computed) against Jbara 2025
(arXiv:2510.13479), which independently publishes C_L1 ~ 3.6210 at mu~0.109 — MATCH. And
C = 3.16794 < C(L1), so the L1 gateway is **open** at this energy, which is a precondition for
any trajectory that passes between Pluto and Charon. This is the closest thing to a true
INDEPENDENCE gate available in this lane (see the tier discussion).

## (c) Novelty framing — **not claimed novel**, and that is the right call

**The literature gate was re-run live by this task on 2026-08-12** (real web search, not
inherited from `#810`). Signature: `primary=Pluto`, `sequence=(Charon,)`, `resonances=("5:1",)`,
`topology_label={repeated-moon}`; the module's own `build_queries` trail was used. **Status:
not-found** — necessary-not-sufficient per `[[feedback_literature_novelty_check_baseline]]`.

Prior art grounded **against the actual artifacts** per
`[[feedback_ground_citations_against_content]]`, not inherited:

- **Ross & Roberts-Tsoukkas 2026, arXiv:2606.29189** (v2, 10 Jul 2026; HTML *and* PDF fetched).
  Computes families **(1,1), (3,3), (3,2), (3,1)** only, at **mu in {0.001, 0.012150584270572,
  0.1, 0.3, 0.5}** (Sun-Jupiter, Earth-Moon, three abstract mu). The strings "Pluto" and
  "Charon" appear nowhere; **no k1=5 family appears anywhere**. `#810`'s characterization is
  accurate. Its conjecture, verbatim: *"We conjecture that saddle-center birth is universal
  among cycler families, suggesting that stable cyclers are a generic feature of three-body
  dynamics"*, alongside *"every cycler family contains a subfamily that is linearly stable to
  both planar and out-of-plane perturbations."*
- **Jbara 2025, arXiv:2510.13479** — surfaced live and DISTINGUISHED. It *is* Pluto-Charon
  CR3BP at mu~0.109, but it studies zero-velocity structures, Lagrange-point instability and
  tadpole/horseshoe chaos: no cyclers, no periodic-orbit families or IC tables, no winding
  numbers, nothing near C=3.168, and the word "cycler" does not appear. Already a known corpus
  anchor in this repo (it is `#494`/`#504`'s C_L1 cross-check source) — verified via
  `docs/notes/`, per `[[feedback_corpus_check_index_not_filenames]]`; no new corpus item.
- **Antoniadou & Libert 2018, arXiv:1805.00288 (CMDA)** — surfaced live and DISTINGUISHED. Its
  "5/1" is a **mean-motion resonance** of a massless body about a star in a star+giant-planet
  system (mu~1e-3), not a (k1,k2) cycler winding; it does not study orbits that alternately
  encounter both primaries, never uses the word "cycler", and never approaches mu~0.109. The
  numeral collision with our "(5,1)" is exactly the concept-collision trap that memory warns
  about.
- **JPL SSD Three-Body Periodic Orbits catalog** does not index Pluto-Charon among its 7
  systems, so `jpl_family_check` cannot adjudicate this row either way.

**The framing decision.** `our_status: known-class-member`, NOT `candidate-novel`. The reasoning,
stated so it can be challenged:

The published class is **generically indexed by winding number** — Ross & Roberts-Tsoukkas
define prograde symmetric cyclers as a (k1,k2)-classified family built from L1 Lyapunov
manifold-tube intersections, and then *conjecture that a linearly stable subfamily exists in
every such family*. Computing the (5,1) entry at mu=0.1087 therefore fills in a slot in a
published, explicitly open-ended classification. It is a **confirming instance of the authors'
own conjecture at a family+mu they never computed** — which is genuinely interesting and worth
cataloguing, but it is not a new dynamical species and calling it "novel" would over-claim.
This is the same logic the `braik-ross-c21-3d-corridor-01-2026` row already applies to itself
("KAM theory GUARANTEES some quasi-periodic torus family exists around any linearly-stable
periodic orbit; `#682`/`#684` measure how BIG it is, not whether it exists"), and it matches the
house style of `#822`'s own careful non-claim. `known-class-member` is also
**not novelty-claimable by construction** (spec §685), which is the correct permanent status
here.

What *is* fairly claimable, and is what the row says: **the first documented member of the (5,1)
prograde cycler family at any mass ratio in this project's records, and the first Pluto-Charon
cycler beyond (3,2)**.

## (d) Tier — V1, assigned against the written criteria

**The `#836` bullet's premise is wrong and this matters for the reasoning.** The bullet states
"the (3,2) precedent `ross-rt-pc-cycler-32-2026` is V1". **It is V2.** Its V1 base is `#494`
(corrector closure + Barden + independent Radau + topology), and `#505` separately upgraded it
to **V2-ballistic** with a 100-period inertial REBOUND/IAS15 run. (The V1 rows are the
*abstract-mu* siblings, e.g. `ross-rt-mu01-cycler-32-2026`.) So V1 here is **not** "the same as
the (3,2) precedent" — it is "the (3,2) precedent's V1 **base**, which is exactly what I
verified and no more".

Against spec §14 as written:

- **V0** (internal consistency, closure residual ≤ tol): met — 1.6e-12.
- **V1** (independent re-solve + re-propagation by a propagator that did not build the
  trajectory): met in this lane's established form — fixed-Jacobi corrector closes, winding
  topology confirmed, Barden |nu|<1, and an **independent Radau** re-propagation passes
  (closure <1e-6, dJ=7.9e-13), plus the algorithm-independent full-period monodromy above.
  This is the identical evidence structure `#494` used for the (3,2) row's V1 base.
- **V2-ballistic** (≥3 continuous laps with bounded rotating-frame drift, evaluated in the row's
  defining model): **NOT met — never attempted for this member.** `#505` did that work for
  (3,2); nothing equivalent exists here. Registered as `#846`.

**V1.** Assigning V2 by analogy to the sibling row would have been precisely the over-claim the
`_LEVEL_EVIDENCE` mechanism exists to prevent.

**The independence declaration, written honestly.** `validate.py`'s own module docstring
(lines 47-52) requires each new entry to name what the cited gate *shares* with the construction
under test, and warns that "every promotion requires at least one TRUE independence gate". The
honest statement, which is in the entry verbatim: the Radau crosscheck and the full-period
monodromy **both share mu, the CR3BP equations of motion, and `state0`** with the construction —
they differ only in integrator/algorithm, so by the file's own taxonomy they are **CONSISTENCY
gates, not INDEPENDENCE gates**. The CR3BP lane has no upstream ephemeris to re-derive, so the
independence available is confined to (i) the SOURCED inputs — mu and lunit_km from
`satellites.py`'s DE440 GM/a values, the clearance radii from Nimmo 2017 — and (ii) the one
external-source physical check, C(L1) vs Jbara 2025. That is the same structural limit the
(3,2) row's V1 base sits under, so the tier is precedent-consistent; the limit is recorded
rather than papered over.

---

## Corrections to the surrounding record

Stated plainly per the task's own instruction not to paper over anything to keep the discovery
narrative clean. None of these affect the member's physics or the verdict.

1. **The `#836` bullet's "the (3,2) precedent is V1" is wrong — it is V2** (V1 base `#494`,
   V2-ballistic upgrade `#505`). Detailed above. The V1 assignment survives, for a better
   reason.
2. **`#810`'s census framing understates what is unknown.** See the next section — this is the
   most substantive finding of the adjudication.
3. **The paper's title is wrong in the catalogue, across ~6 rows.** The title printed on page 1
   of the arXiv PDF is **"Stable Ballistic Prograde Cyclers in the Three-Body Problem"**; the
   arXiv *listing metadata* renders it as "Stable Families of Ballistic Prograde Cyclers in the
   Restricted Three-Body Problem"; the catalogue rows carry a third, paraphrased string,
   "Families of Stable Prograde Cycler Orbits in the Circular Restricted Three-Body Problem",
   which matches neither. `#810`'s note inherits the paraphrase. Those rows also cite **v1**
   while the live paper is **v2** (10 Jul 2026), which is what I grounded against. The new row
   uses the PDF-printed title and records the variants; sibling rows are **not** touched here
   (concurrent sessions are editing `catalogue.yaml`) — `#845` is registered to reconcile all of
   them at once.

## The census is NOT "1 of 15" — it is "1 confirmed, 10 certified-empty, 3 unresolved"

`#810`'s note asserts: *"Other 8 higher topologies' certified-empty negatives are untouched
(their grids found no seed at all; nothing downstream to lose)."* I audited that claim against
the `#504`/`#549`/`#656`/`#807` records, classifying every one of the 15 (k2<=k1<=5) topologies
by its **negative mechanism**: (A) no seed found at all — not vulnerable to the `hc=None`
auto-redetection branch-loss gap; (B) a seed or intermediate orbit *was* found under `hc=None`
and the branch was then lost or no window was found — vulnerable to exactly the gap that hid the
(5,1) member; (C) positive.

**As narrowly stated, `#810`'s sentence is accurate**: it is scoped to `#656`'s nine higher-k1
topologies minus (5,1), and all eight of those genuinely had `_grid_seed_search` return nothing.

**But the census prose built on it is not.** The narrow scoping silently excludes `#504`'s three
*lower*-k1 negatives, which used the same vulnerable `hc=None` machinery:

- **(3,1)** — the strongest concern. Its own record shows a **converged, stable orbit found and
  then discarded for wrong topology** (`reaches_secondary=False`, near-primary). That is the
  *structurally identical signature* to the (5,1) case that has just turned out to be a real
  false negative. Never re-swept with fixed hc.
- **(3,3)** — `#807` **directly measured** the mu-continuation leaving the (3,3) branch at its
  very first step under `hc=None`, landing on an unrelated family. `#807` only fixed the
  *reporting* of that ("stable/wrong-topology" → "clean negative"), not the underlying
  continuation. The true (3,3) family at PC mu has never been explored past the branch-loss
  point.
- **(1,1)** — weakest of the three (no confirmed wrong-topology capture on record, just "no
  stable window" from an `hc=None` C-sweep of unknown branch fidelity), but the identical code
  path with no independent evidence ruling out branch loss.

Final classification: **bucket A (genuinely certified-empty, no seed) = 10** —
(2,1),(2,2),(4,1),(4,2),(4,3),(4,4),(5,2),(5,3),(5,4),(5,5); **bucket B (unresolved, same bug
class) = 3** — (1,1),(3,1),(3,3); **bucket C (positive) = 2** — (3,2) and now (5,1).

So the honest census is **"2 stable members confirmed of 15; 10 certified-empty; 3 unresolved
under the same method gap that just flipped (5,1)"** — *not* a clean "1 of 15 beyond (3,2)",
which would imply the other fourteen are settled. Per
`[[feedback_bugfix_invalidates_past_searches]]` a limited search path is a false-negative
generator, and the (5,1) result is direct evidence that this particular one generated at least
one. The `#504`/`#549`/`#656` census prose is rewritten accordingly, and **`#844` is registered
to re-sweep (1,1)/(3,1)/(3,3) with fixed hc**, prioritising (3,1) and (3,3) whose branch-loss
evidence is already on record.

`#810`'s other structural claims check out: PC (3,2) is indeed no longer structurally unique at
this mu, the `pluto-charon-kk-45-cycler-sweep-2026-07-19` empty-region stamp is correctly left
unedited (append-only, and its (5,1) row already said UNSETTLED-not-certified-empty), and no new
empty-region stamp is warranted (a positive is not an empty region).

## Follow-ups registered

- **`#844`** — fixed-hc re-sweep of PC (1,1), (3,1), (3,3): the three bucket-B negatives whose
  `hc=None` branch-loss exposure is the same one that hid the (5,1) member. Highest-value item
  this adjudication produced.
- **`#845`** — reconcile the Ross & Roberts-Tsoukkas 2026 citation across all ~6 rows citing
  arXiv:2606.29189: title string (three variants in play) and v1-vs-v2 version drift.
- **`#846`** — optional V2-ballistic upgrade campaign for `pc-cycler-51-2026` by `#505`'s
  100-period REBOUND/IAS15 method, the only thing standing between this row and V2.

## Verification

`uv run pytest tests/data tests/search -q` (full ratchet, mandatory after the `catalogue.yaml`
edit), `tests/scripts` (the `validate.py` edit + the `preflight_search` AST ratchet),
`ruff check .` / `ruff format --check .`, full `uv run mypy src tests` — see the `#836` commits.
