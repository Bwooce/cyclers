# #714 — CCR5BP/CRNBP (N≥5-body) discovery-strategy pass (2026-07-27)

Analysis-only (no code, no catalogue writes, no dispatches), mirroring `#686`'s format and —
per the dispatch's own mandate — leading with an honest TRACTABILITY verdict before any
shortlist. The question: with the CCR4BP arc (`#689`-`#708`) complete and one genuine novel
discovery banked (`umbriel-1-2-torus-homoclinic-uranus-2026`), is a genuine N=5 CRNBP
discovery build worth attempting now, grounded on the three papers acquired in
`#710`/`#711`/`#712`?

Inputs actually read (not summarized from memory): the three digests in full
(`docs/notes/2026-07-26-digest-negri-prado-2022-crnbp.md`,
`docs/notes/2026-07-26-digest-gilliam-bettinger-2024-crnbp-jovian.md`,
`docs/notes/2026-07-26-710-digest-aryan-fitzgerald-2024-jovian-pccfbp.md`); the `#686`
strategy pass in full (`docs/notes/2026-07-22-686-nbody-discovery-strategy-pass.md`);
`src/cyclerfinder/core/ccr4bp.py` in full (all 441 lines);
`search/ccr4bp_heteroclinic_search.py` and `search/variational_ccr4bp_torus.py` (docstrings +
structure); the `#693`-`#696`, `#694`, `#701`-`#703` OUTSTANDING.md bullets for what the
Jovian/Uranian CCR4BP campaigns actually concluded; plus two live numeric checks against the
in-repo JPL SSD registry (`core.satellites`) computed during this pass (reported inline below
— these are this pass's own arithmetic, re-runnable in one line, not sourced claims).

## 1. The three tractability sub-questions, answered concretely

### 1(a) The EOM extension: genuinely small, with one real transcription hazard

Judged against the actual `ccr4bp.py` source, not the digests' word. The N=4 module's
structure is: base CR3BP RHS + per-perturber direct term `-mu_gan*(r-r_gan)/|r-r_gan|^3` +
indirect term `-mu_gan*r_gan/a_gan^3`, with the STM's A-matrix adding only the direct term's
Hessian block (`_ganymede_second_deriv_block`) because the indirect term is independent of
spacecraft position. Extending to N=5 per Negri & Prado Eq. 11 requires:

1. **A second `(mu_j, R_j, psi_j)` perturber triple** — mechanical; `#695`/`#696`/`#701`
   already demonstrated the parameterization pattern generalizes across systems.
2. **The pairwise inner-sum coupling term** (the genuinely new N≥5 physics — mutual
   correction for each perturber's pull on the frame origin via the other perturber). The
   load-bearing code fact, verified from the digests' transcription of Eq. 11: this term
   depends ONLY on `(R_j, psi_j, R_k, psi_k)` — it is INDEPENDENT of the spacecraft state.
   It is therefore a known time-periodic forcing added to the accelerations and contributes
   **exactly zero to the variational Jacobian**: `ccr4bp_stm_eom`'s A-matrix needs only a
   second copy of the existing direct-term Hessian block for the new moon, nothing else. The
   STM extension is as mechanical as the RHS extension.
3. **NOT a naive superposition**: the Gilliam digest's own warning stands — one cannot just
   call `_ganymede_acceleration` twice; the cross term must be transcribed from Eq. 11
   verbatim.

**One real hazard found during this pass, previously unflagged:** the two digests transcribe
the coupling term with **opposite signs**. Negri & Prado digest §2:
`- sum_j mu_j [ direct_j  -  sum_{k!=j} mu_k (...)/rkj^3 ]`; Gilliam digest §1:
`- sum_j mu_j [ direct_j  +  sum_{k!=j} mu_k (...)/rkj^3 ]` — with the numerator
`(Rj*cos(psi_j) - Rk*cos(psi_k))` written identically in both. Both digests claim
term-for-term identity with the other, so at least one transcription has a sign slip. This
is exactly the class of error the structural-reduction tests would catch late and painfully;
it must be resolved by reading Eq. 11 and thesis Eqs. 25-27 directly from the two source PDFs
BEFORE implementation, and then adjudicated by the digit-grade checks in §1(b). With that
resolution step included, the extension is a bounded, well-specified, ~2-4-day core build —
the digests' "small change" claim survives direct inspection of the code it would change.

### 1(b) The positive-control question: what tier actually exists (be precise)

The task's suspicion about the Gilliam tables is **correct**. Precisely tiered:

- **Tier 0 (structural, in-house, free):** `mu_4 -> 0` must reduce the N=5 module
  byte-exactly to `ccr4bp.py`'s N=4 (whose own `mu_gan -> 0` -> CR3BP chain is already
  ratcheted). Also adjudicates the §1(a) sign question: with one perturber the coupling term
  must reduce to the existing direct+indirect form exactly.
- **Tier 1 (digit-grade, literature-sourced, ALGEBRAIC only):** Gilliam's thesis Tables 5-6
  "Lagrange Box" numbers (Jupiter-Europa E1 dx/dy = 57.86/44.07 km, Uranus-Oberon
  77.99/128.07 km, etc., six systems). These are **static instantaneous force-balance loci**
  — CRNBP-perturbed analogues of the CR3BP libration points, found by Newton-Raphson on the
  acceleration field. They are NOT periodic orbits, NOT connections, NOT even dynamical
  substitutes in the orbit sense — exactly "library-point existence" tier, as the dispatch
  suspected. What they DO validate, digit-grade: the N≥5 force field AND (via the
  Newton-Raphson solve) the analytic Jacobian — i.e., precisely the two things §1(a) builds.
  Genuinely valuable as an EOM-implementation control; worthless as a corrector/search
  control.
- **Tier 2 (periodic/quasi-periodic orbit): NOTHING exists in the literature.** Gilliam's
  own explicit admission ("not currently known if any periodic trajectories exist in the
  CRNBP" outside resonance-locked cases) stands; Negri & Prado's Jupiter-Ganymede
  vertical-Lyapunov epsilon-homotopy continuation is figure-only, no IC table — fails
  `[[feedback_golden_tests_sourced_only]]` as a golden target.
- **Tier 3 (connection-level): nothing at all.**

**But the existence risk is materially narrower than Gilliam's quote suggests — for one
specific configuration.** Gilliam's open question is about the GENERIC CRNBP, where two
perturbers' synodic rates are incommensurate, the forcing is quasi-periodic, and exact
periodic orbits generically cannot exist (this is `#686` §1(f)'s "over the boundary" regime
— it killed the general N≥4 lane then and it still does). The one escape is the same one
`#686` found at N=4: a configuration where the forcing is genuinely TIME-PERIODIC, so the
target objects are 2D invariant tori of a periodically-forced planar system — the EXACT
mathematical setting `#690`-`#694` already solved, with KAM-persistence logic intact and the
perturbing masses still perturbatively small. §1(c) shows the solar system offers exactly
one such configuration, and it is better than expected.

Given periodicity, the honest positive-control position is: **no literature-grade dynamical
control exists or can exist yet (the frontier is genuinely unclaimed); the control must be
SELF-GENERATED by continuation from this project's own validated N=4 limits.** Concretely:
`mu_Io: 0 -> physical` continuation of `#690`'s converged Jupiter-Europa 3:4 CCR4BP torus,
where the `mu_Io = 0` endpoint IS the known-good (machine-precision reproduction of `#690`'s
torus, itself validated against Kumar 2021's object class), and every continuation step is
seeded from a converged neighbor. This is the same pattern `#690`'s own mass continuation
already used, and the same epsilon-homotopy method Negri & Prado themselves demonstrate.
It satisfies `[[feedback_verify_gauntlet_with_positive_control]]` in the only way available
— and the discipline record for any N=5 task must say plainly that the control is
continuation-anchored, not literature-anchored. A solver run NOT reachable by continuation
from a validated limit would have no control at all and must not be trusted on a 0/N.

### 1(c) Which real systems qualify: exactly one — and the Laplace lock is stronger than the papers used

For N=5 time-periodicity the requirement is NOT merely "three moons with pairwise
near-resonances": the TWO perturbers' synodic rates relative to the base secondary must be
commensurate. A pairwise two-moon resonance only ever buys the N=4 CCR4BP (already built).
Commensurability of two synodic rates against a third body requires a genuine THREE-moon
resonance chain. Census:

- **Jupiter, Io-Europa-Ganymede (Laplace resonance): the only qualifying trio in the solar
  system — and it qualifies exactly.** The key fact, verified numerically this pass: the
  Laplace relation `n_Io - 3*n_E + 2*n_G = 0` algebraically forces the two perturber synodic
  rates to an exact small-integer ratio in EVERY base frame of the trio: Europa base
  `omega_Io/omega_Gan = -2` exactly; Ganymede base ratio `3`; Io base ratio `3/2`. From
  observed sidereal periods (Io 1.769137786 d, Europa 3.551181041 d, Ganymede 7.15455296 d)
  the Europa-frame ratio computes to **-2.000000001** (Laplace residual 5.6e-10 per n_E) —
  the physical system enforces the commensurability to observational precision. So the N=5
  Jovian CRNBP in the Jupiter-Europa frame is time-periodic at PHYSICAL frequencies, with
  fundamental period = Ganymede's synodic period — **the SAME single forcing clock
  `theta1` that `#690`'s corrector already uses**, with Io's phase an exact slaved function
  `theta_io(t) = theta_io0 - 2*(theta_gan(t) - theta_gan0)`. No 3-torus, no second free
  clock, no structural corrector jump. Notably, Negri & Prado's own worked example idealized
  the PAIRWISE periods to 1:2:4 (a ~7e-3 relative frequency distortion); the synodic-rate
  route needs only a projection onto the Laplace constraint (~2e-4 relative on the
  registry-SMA-derived rates, whose deviation is registry rounding, not physics; the
  observed-period route is 5.6e-10) — ~35x milder than the idealization the framework paper
  itself used. Two bonuses: (i) registry-computed `mu_Io = 4.70e-5` (Jupiter+Europa
  normalizer) sits between `#695`'s and `#690`'s validated forcing strengths — squarely in
  the proven perturbative regime; (ii) the Laplace libration (`lambda_Io - 3*lambda_E +
  2*lambda_G` librating about 180 deg with small amplitude, ~0.06 deg literature value, to
  be re-confirmed at build time) PINS the relative initial phase `theta_io0` given
  `theta_gan0` — the model has ONE free epoch phase, not two, a physically-enforced
  reduction of the search space.
- **Callisto is NOT in the Laplace chain** (period ratio to Ganymede ~2.33, non-resonant):
  adding it (N=6) breaks periodicity — excluded, and with it the Gilliam CRNBP6/CRNBP10
  cases as discovery targets (they are sensitivity studies, correctly so).
- **Saturn:** pairwise MMRs only (Mimas-Tethys 2:1, Enceladus-Dione 2:1, Titan-Hyperion
  4:3); no three-moon chain; and `#693` already quantified Saturnian `mu_pert` at 40-8,000x
  below the JEG reference. Excluded twice over.
- **Uranus:** no present-day MMRs among the major moons (`#686`/`#693`) — an N=5 Uranian
  model has irreducibly quasi-periodic forcing. Excluded; the Uranian CCR4BP discovery
  (`#701`) has NO N=5 extension lane, which is worth stating explicitly since Uranus is this
  project's home turf.
- **Neptune, Pluto:** Triton retrograde/inclined + small moons mass-dead (`#693`, `#599`);
  Pluto's near 1:3:4:5:6 chain with Charon involves perturber masses ~1e-9 (Gilliam's own
  Pluto case is a model-necessity demo, not a structure claim). Excluded.

So the answer to the dispatch's question is: the Galilean Laplace trio is not merely the
obvious first candidate — it is the ONLY realistic one, and (the genuinely new observation
of this pass) its qualification is exact-to-observation rather than approximate. Whether the
"Laplace relation => exact synodic commensurability => strictly periodic N=5 CRNBP at
physical frequencies" observation is itself unremarked in the literature is a lit-check
question, not a claim made here — Negri & Prado's choice to idealize pairwise suggests they
did not use it, but the underlying celestial mechanics is classical.

## 2. Tractability verdict (first, per the mandate)

**TRACTABLE WITH CAVEATS — as a single-lane, staged, mostly-reuse extension; NOT tractable
as a general CRNBP discovery program.** Stated against `#686` §2(c)'s own four-axis wall
test: (i) forcing weak — `mu_Io = 4.70e-5`, inside the regime `#690`-`#701` already proved
convergent, and the coupling term is O(mu_j*mu_k) ~ 4e-9, a tiny correction to a small
forcing; (ii) anchor objects are the SAME tori `#690`/`#691` already converged and
whisker-validated (segmented-CLV worst-case 0.0083 deg); (iii) no scale gap — Io at 0.629
and Ganymede at 1.595 Europa-SMA in one frame; (iv) positive controls — **the one degraded
axis**: digit-grade at EOM tier only (Gilliam tables, static), nothing dynamical in the
literature, control must be self-generated by continuation from the validated N=4 limit.
3.5 of 4 axes clear. The `#533`-`#620` wall does not transfer here any more than it did to
the CCR4BP; the QBCP arc's specific lesson (no positive control at the END target ⇒ weeks
of undiagnosable convergence walls) is mitigated — not eliminated — by the continuation
architecture, because every step of a continuation has the previous step as its control, so
a wall is immediately localized in `mu_Io` rather than ambiguous between physics and
tooling.

The caveats, plainly:
1. **No literature golden test exists at torus or connection tier, and none can be
   manufactured by more searching** — the frontier is open (which is also the opportunity).
   Every N=5 result stands on continuation provenance + the existing ghost-guard
   discipline + independent-integrator checks. This must be stated in any writeback.
2. **The digest sign discrepancy (§1(a)) must be resolved from the source PDFs before any
   code is written.** Cheap, but mandatory sequencing.
3. **Exact periodicity requires projecting the registry-derived synodic rates onto the
   Laplace constraint** (~2e-4 relative adjustment) — a model idealization, though one the
   physical system justifies to 5.6e-10; the same class of one-value change `ccr4bp.py`'s
   own docstring already reserves for the exact-2:1 option, and far milder than the pairwise
   idealization in the framework paper itself.
4. **The discovery upside is conditional**: the headline questions are "does the first-known
   N=5 CRNBP invariant torus exist at physical `mu_Io`" and "does `#694`'s
   ghost-guard-verified JEG homoclinic connection survive genuine 5-body forcing." A YES is
   a strong, likely-publishable first (subject to the mandatory
   `search/literature_check.py` gate before the word "novel" or "first" is used); a
   characterized breakdown is a registry-grade negative with a named mechanism. Both are
   acceptable outcomes; neither is a new cycler row by itself.
5. **Everything outside the Laplace-locked Jovian lane is NOT tractable** and should be
   declined without new evidence: no other qualifying trio exists (§1(c)), and unlocking
   the frequencies (real-rate quasi-periodic forcing) crosses `#686` §1(f)'s boundary into
   the regime where Gilliam's existence question genuinely bites and no reduction exists.

## 3. Shortlist (ranked; scope/effort/reuse stated per item)

### 1. CRNBP core: EOM + STM extension with digit-grade validation (GATE for everything else)

**What.** Resolve the Eq. 11 sign question from both source PDFs (Negri & Prado Eq. 11 vs
Gilliam thesis Eqs. 25-27) and record the reconciliation in the digest notes; then a
`core/crnbp.py` (or a perturber-list generalization of `CCR4BPSystem` — implementer's
choice, but `ccr4bp.py` itself stays untouched per the `#689` discipline) implementing the
general-N EOM + STM with the coupling term; structural ratchet tests (`mu_4 -> 0` byte-exact
to `ccr4bp.py`; chain down to CR3BP); a small Newton-Raphson equilibrium solver reproducing
Gilliam's Tables 5-6 Jupiter-Europa E1/E2 numbers digit-grade (this validates BOTH the force
field and the analytic Jacobian — the thesis's Uranus-Oberon row is a free cross-system
second check); a Laplace-projection constructor
(`jupiter_europa_io_ganymede_default()`) documenting the ~2e-4 rate projection with the
observed-period 5.6e-10 justification. **Effort:** ~2-4 days. **Reuse:** the whole `#689`
test/discipline pattern. **Risk:** none research-grade; spec-complete —
Sonnet-tier per `[[feedback_subagent_model_tiering]]`, with the sign-reconciliation read
done first and checked by the coordinating session.

### 2. First N=5 invariant torus: `mu_Io` continuation of the validated JEG torus (GATED on 1)

**What.** Swap the forcing evaluation in `search/variational_ccr4bp_torus.py`
(`_ganymede_on_theta1` -> a two-moon `_moons_on_theta1` with Io's phase slaved to `theta1`
by the exact -2 ratio — `omega1` and the whole corrector structure unchanged, a smaller
delta than the `#617`->`#690` EOM swap was) and continue `#690`'s converged Europa 3:4 torus
in `mu_Io: 0 -> 4.70e-5`. Positive control: the `mu_Io = 0` endpoint must reproduce `#690`'s
torus to machine precision (Tier-0), every step seeded from its converged neighbor. Outcome:
either the first known N=5 CRNBP invariant torus (then run the mandatory literature-novelty
check before any "first" language) or a characterized breakdown curve in `mu_Io` (registry-
grade negative with mechanism). Sweep `theta_io0` over the Laplace-pinned value ± the
libration amplitude as a cheap robustness check. **Effort:** ~1-2 weeks. **Reuse:** `#690`
corrector verbatim minus the forcing function; `#689` constructor pattern. **Risk:** the
genuine open-physics question, deliberately taken in its cheapest falsifiable form.

### 3. Connection survival: does `#694`'s JEG homoclinic connection survive Io? (GATED on 2)

**What.** Re-run the `#691` whisker -> `#694` globalize -> heteroclinic search -> ghost-guard
chain on the surviving N=5 torus, asking specifically whether the ghost-guard-verified
Europa-3:4 homoclinic connection (residual ~1e-14, Radau/DOP853 <1 km, mesh-refinement-
stable) persists, deforms, or dies under Io's forcing — the N=5 analogue of the project's
banked N=4 result, and the first connection-level object in any CRNBP if it survives.
Requires threading the system-parameterization through the `#694` modules (they import
`core.ccr4bp` directly) — and while in there, fix the hardcoded `_L_KM`/`_v_unit_km_s`
Europa constants properly (the `#694`/`#696`-documented unit bug currently worked around
per-driver) rather than adding a third workaround. Optional extension, only if the main
question resolves cleanly: the Io-base-frame N=5 model (ratio 3/2) to see whether `#695`'s
honest Io-Europa near-miss changes character under Ganymede's added forcing. **Effort:**
~1-2 weeks. **Reuse:** all four `#691`/`#694` modules, refactored not rebuilt. **Risk:**
moderate; every verdict passes the existing ghost-guard + independent-integrator + mesh-
refinement discipline, which `#701`/`#702` already battle-tested.

No fourth item. The Saturn/Uranus/Neptune/Pluto N=5 lanes and the N=6 (+Callisto) lane fail
§1(c) structurally; a 3-angle (quasi-periodically-forced) torus corrector would be a
structural build against the `#611`-documented jump with zero controls at any tier — both
declined without new evidence, mirroring `#686` §4's discipline.

## 4. Considered and explicitly rejected (do not re-surface without new evidence)

- **Generic CRNBP discovery in any non-Laplace system** (incl. all-Uranian N=5): forcing
  irreducibly quasi-periodic; Gilliam's periodic-existence open question applies in full;
  no reduction, no control tier above algebraic. This is `#686` §1(f) restated with sharper
  evidence.
- **N=6+ Jovian (add Callisto), CRNBP10-style:** breaks periodicity (Callisto non-resonant);
  Gilliam's own use of these cases is divergence/sensitivity study, not structure.
- **3-angle torus corrector for quasi-periodic forcing:** structural jump (`#611`),
  dimension-raised tori, no published success at discovery scale, no control. The genuinely
  "impossible" part of the problem; refusing it is the strategy.
- **Reproducing Negri & Prado's worked examples as positive controls:** all three are
  figure-only (their §3; digest §4) — fails the sourced-golden discipline. The
  epsilon-homotopy METHOD is adopted (shortlist item 2); the figures are not targets.
- **Reproducing Gilliam Ch. VI Poincaré maps:** qualitative map-structure comparisons, no
  digit-grade content; adds nothing over the Tables 5-6 equilibrium check.
- **Treating the Gilliam equilibrium tables as a dynamical control:** they are static
  force-balance loci (§1(b)); using them to "validate" a torus corrector would be exactly
  the tier confusion the dispatch warned about.

## 5. User decision points (flagged, not assumed)

1. **GO/NO-GO on the staged lane** (items 1->2->3, ~3-5 weeks total if all gates pass —
   smaller than `#686` Stage B was, because the pipeline is reused rather than built).
   Item 1 is cheap and stands alone as pure model-capability + corpus-grounding work even
   under a NO on 2/3.
2. **Novelty-claim handling:** "first N=5 CRNBP torus/connection" claims are strong; the
   `literature_check.py` gate is necessary-not-sufficient as always, and the
   Laplace-commensurability observation itself needs a targeted lit-check before being
   called new (respectful-errata/novelty framing discipline).
3. **Control provenance disclosure:** any N=5 result's writeup must state its positive
   control is continuation-anchored (self-generated from the validated N=4 limit), not
   literature-anchored — a weaker tier than `#690`'s Kumar-anchored control, disclosed, not
   hidden.
4. **Schema:** an N=5 model-native object needs a `model_assumption` value beyond `#707`'s
   CCR4BP work (and the Laplace-projection idealization needs the same
   `validity_window`-style disclosure the exact-2:1 option would have needed).
5. **Gilliam-table mismatch policy:** Tables 5-6 are thesis content (public, but not the
   peer-reviewed CMDA subset). If our reproduction disagrees, investigate both sides before
   assuming our bug — the published-rounded-values and errata disciplines both apply.

## 6. Recommended dispatch order

Item 1 first, alone (it contains the mandatory sign-reconciliation read and the only
digit-grade external check; everything downstream inherits its correctness). Item 2 on its
completion, as a single owned long-run task with checkpointed continuation state (long-runs
doctrine; no subagent backgrounding). Item 3 only on item 2's torus surviving to physical
`mu_Io`, dispatched with the `#694`-module refactor scoped in. If item 2 finds a breakdown
instead, item 3 is cancelled and the breakdown curve itself is written up + registry-stamped
— the honest position being that the solar system's one N=5-periodic configuration was
tested at the right tier and answered, which is what a discovery program is for.
