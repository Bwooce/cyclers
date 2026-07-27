# #724 — Final confirmation pass: N=5 CRNBP torus narrow-novelty claim (2026-07-27)

**Task:** independently skeptical re-check of `#722`'s refined verdict (a narrow
residual novelty claim survives for `#720`/`#723`'s phase-corrected N=5 torus),
NOT a restatement of it. Re-read both TCP source papers in full hunting for
anything `#722` under-weighted; re-ran the full continuation pipeline in the
foreground; hunted for counter-evidence and follow-on literature; sanity-checked
the Kumar et al. 2021 grounding. All work synchronous, 2026-07-27 AET.
Reproduction script: `scripts/verify_724_rerun_continuation.py` (verbatim the
run quoted in §2).

## VERDICT: CONFIRMED — with one substantive correction to the claim language

The narrow novelty claim SURVIVES this independent re-check. Neither TCP paper
computes, sweeps through, or anticipates-in-executed-form a mean-motion-resonant
torus substitute; no follow-on work closing the gap was found; the delivered
object is exactly what `#722` said it is, at the corrected phase, reproduced
here bit-for-bit.

**But `#722`'s (and this task's own dispatch's) description of the orbit as
"interior" is WRONG and must not enter any catalogue writeback.** The Kumar
et al. 2021 Jupiter-Europa 3:4 resonant orbit is EXTERIOR to Europa:

- `resonance_semimajor(3, 4) = (4/3)^(2/3) = 1.2114` Europa-SMA units > 1
  (spacecraft makes 3 revolutions per 4 Europa revolutions, so it is the
  SLOWER, outer body). Verified directly from
  `src/cyclerfinder/genome/composed_moon_map.py` and the actual seed orbit
  (`x0 = 1.6102` at apoapse-type crossing, SMA 1.2114).
- Kumar et al. 2021's own abstract (arXiv:2109.14815, p.1): "we compute tori
  corresponding to **exterior Jupiter-Europa** and **interior Jupiter-Ganymede**
  PCRTBP resonant periodic orbits" — the authors' own terminology, checked
  against the source PDF this pass.
- The orbit lies BETWEEN the two moons (1.2114 Europa units = 0.759 Ganymede
  units, i.e. near Ganymede's interior 3:2 — the m:n vs m:2n Laplace
  correspondence Kumar 2021's own §"Relations between..." explains). Calling it
  "interior" full stop is factually false; "interior" appears to have bled in
  from the UNRELATED Kumar/Rawat Earth-Moon interior-MMR papers elsewhere in
  corpus.

## 1. Independent re-read of both TCP papers — the absence is textual, not inferred

Both papers re-read IN FULL from the PDFs (not from `#722`'s digest):
`cyclers_pdf/papers/baresi-owen-scheeres-2023-...-AAS-23-201.pdf` (16 pp) and
`owen-baresi-scheeres-2024-...-issfd2024.pdf` (7 pp). Complete enumeration of
every computed object in each:

**Paper 1 (AAS 23-201):**
- Figs. 3-5: L1/L2 equilibrium dynamical substitutes (periodic orbits), all
  three Jupiter-moon frames.
- Figs. 6-7: manifolds of those equilibrium substitutes (no intersections found
  in JCI — their stated motivation to move to torus substitutes).
- Fig. 8a-d: quasi-periodic torus substitutes of the **30th/50th L2 Jupiter-
  Europa, 50th L1 Jupiter-Ganymede, 30th L1 Jupiter-Io Lyapunov planar family
  members** — the caption and body text ("different members of the L1 and L2
  Lyapunov planar families") identify every torus as Lyapunov-family.
- Figs. 9-11: Floquet directions + manifolds of those same Lyapunov tori.
- No appendix, no supplementary tables (Table 1 = system constants only).
- The words "resonance/resonant" occur ONLY for the Laplace resonance of the
  model construction itself — never for a spacecraft mean-motion-resonant
  orbit/torus.

**Paper 2 (ISSFD 2024):**
- §V (its own words): "The QPOs at Ganymede and Europa were quasi-periodic
  continuations of the **70th and 75th members of their Lyapunov families**,
  respectively." Planar only (§IV: "we will only be considering ... the planar
  tri-circular restricted five-body problem").
- "Resonance" in its Results = resonant TRANSFER ARCS (apoapsis-reduction,
  TOFs at integer multiples of Ganymede's 7.16-d period, Ross et al. 2004
  lineage) — a property of the connecting trajectories, not a resonant torus.
- **The absence is confirmed by the paper's own future-work text, not merely
  inferred from silence** (§VI): "There are many periodic orbits beyond
  Lyapunov, include [sic] those with vertical components, which may have
  corresponding dynamical substitutes in the tri-circular restricted three-body
  problem." I.e. the authors themselves state that orbits beyond the Lyapunov
  family were NOT computed and are future work. The extension direction is
  anticipated; the computation is not executed.

**Structural point `#722` noted that this pass sharpens:** TCP's periodicity
(Eq. 3, `f(t̃+2πk, X) = f(t̃, X)`, period `2πk`) EXISTS ONLY BECAUSE Table 1
locks all synodic rates to exact rationals (J-E frame: Io = +1, Ganymede = -0.5
exactly). Physical, non-rationalized rates are not a variant TCP could adopt
without abandoning its own `2πk` periodicity and fixed rotation number
`ρ = T/k`. Our formulation keeps Ganymede's PHYSICAL registry synodic rate
(-0.5035527; ephemeris-period value -0.5036474, see §2) and gets periodic
forcing from a single clock with Io exactly Laplace-slaved
(`ω_Io = -2·ω_Gan`) — same object class, but at the physical fundamental
frequency (~0.73% off TCP's idealized one, ~40× their own Table-1 parameter
precision). This is a genuine, structural (not cosmetic) rate-model difference.

**Phase is NOT a differentiator and must not be claimed as one:** TCP Table 1's
Jupiter-Europa column (`φ_Io0 = π, φ_Gan0 = 0`) is exactly the physical
libration-center phase our post-`#723` default now uses. They had it right
first; `#723` brought us into agreement with them (and with Sinclair 1975 /
Gilliam). "At the physical Laplace phase" may appear descriptively but never as
a novelty axis.

## 2. Reproduction of the delivered object (foreground, this pass)

Re-ran the exact `#720`/`#723` pipeline (same solver settings as
`tests/search/test_variational_crnbp_torus.py::test_continuation_reaches_physical_mu_io_at_n1_2`):
Europa 3:4 symmetric resonant CR3BP orbit (perp residual 7.2e-13) → `#690`
CCR4BP torus, n1=2/n2=20 (residual_rms 1.2210263312e-04) → mu_io=0 CRNBP seed at
`theta_io0 = π` (1.2209943365e-04) → 8-step mu_Io continuation to the physical
4.7043400305e-05:

```
step 1: mu_io=5.880e-06  residual_rms=1.2208996410e-04  rot=0.496471201
...monotone (one -9.5e-10 re-solve dip at step 1, as #723 documented)...
step 8: mu_io=4.704340e-05  residual_rms=1.2343143649e-04  rot=0.496468269
```

- Final residual **1.2343143649e-04** — matches `#723`'s reported value to
  4e-11 relative. Closure residual 2.321e-3 (inside the 5e-3 gate; same
  n2=20 truncation floor as the `#690` baseline — `#721` §2's off-grid caveat
  stands: quote ~2.1e-3 off-grid invariance alongside the on-grid 1.23e-4 in
  any writeback).
- Seed provenance confirmed in code: `_resonant_symmetric_orbit(mu, 3, 4)` →
  `vt.discover_ccr4bp_torus_from_resonant_orbit` (`#690`'s corrector) →
  `vc.discover_crnbp_torus_from_ccr4bp_seed` — the Kumar-class object, not a
  Lyapunov orbit, at every stage.
- Rates confirmed in code AND against an independent ephemeris check (IAU/JPL
  sidereal periods 1.769137786 / 3.551181041 / 7.15455296 d):
  `omega_gan(code) = -0.5035526857` vs ephemeris-period value `-0.5036473892`
  (1.9e-4 relative, registry-SMA rounding, as `crnbp.py`'s docstring documents);
  both are ~7.3e-3 away from TCP's idealized -0.5. Ephemeris Io/Gan synodic
  ratio = -2.000000001 (the Laplace lock the code projects exactly).
- Phase confirmed: `theta_io0 = π`, `phi_L = theta_io0 + 2·theta_gan0 = 180°`
  exactly, in the delivered final system.

## 3. Counter-evidence hunt (what would overturn the claim) — nothing found

Actively looked for ways the gap could already be closed:

1. **Could a TCP Lyapunov continuation sweep through the resonant
   configuration as a side effect?** No. The Lyapunov planar family members
   computed (up to the 75th) are libration-region objects near x ≈ L1/L2 of
   each moon (their Figs. 8-10 stay in the moon's vicinity); the Europa 3:4
   resonant orbit is a Jupiter-circling orbit at a ≈ 1.211 Europa units,
   BETWEEN the moons' orbits — a topologically and geographically different
   region. A planar Lyapunov member does not pass through or approximate it at
   any family index shown, and neither paper computes any bifurcated branch.
2. **Follow-on TCP literature.** Web-searched (multiple phrasings, 2025-2026
   windows): no journal version of either TCP paper, no Owen thesis output,
   and no third TCP paper computing beyond-Lyapunov substitutes has appeared.
   Owen's other published work (knot-theory heteroclinic screening, in corpus)
   is CR3BP-based, not TCP. Paper 2's future work (spatial case,
   beyond-Lyapunov) remains, as far as found, unexecuted in print.
3. **Kumar lineage / N=5 anticipation.** Kumar et al. 2021 itself mentions Io
   and the Laplace resonance ONLY to explain the m:n ↔ m:2n Europa/Ganymede
   correspondence; its stated future work is heteroclinic connections in the
   N=4 CCR4BP — no N=5/multi-perturber extension attempted or announced. That
   future work was executed as Kumar, Anderson & de la Llave, "Transfers
   between Jupiter-Ganymede and Jupiter-Europa resonant tori in a concentric
   circular restricted 4-body model," **Acta Astronautica 211 (2023)** — still
   N=4 by its own title/abstract (GPU manifold near-intersections,
   Ganymede 4:3 → Europa 3:4). **This paper is NOT in corpus** — flagged as an
   acquisition candidate (medium-high: it is the resonant-torus analogue of the
   TCP transfer paper and the direct model for any future `#714`-item-3-style
   N=5 connection work). Its existence does not touch the N=5 torus claim.
4. **Gilliam §2.2.5.1** ("literature results produced for this model can be
   replicated in the... CRNBP with ease") is a model-equivalence remark, not a
   computed resonant torus; no other thesis chapter computes one.

## 4. The precise, defensible claim (proposed language)

**Headline (narrow, defensible):**

> First computed quasi-periodic invariant-torus substitute of a mean-motion-
> resonant periodic orbit in a Laplace-locked Jupiter-Io-Europa-Ganymede
> restricted five-body model: the planar Jupiter-Europa 3:4 resonant orbit of
> Kumar et al. 2021 (exterior to Europa, a ≈ 1.211 Europa SMA), continued to
> the physical Io mass with Ganymede at its physical (non-rate-idealized)
> synodic rate and Io exactly Laplace-slaved (ω_Io = -2·ω_Gan), at the
> physical libration-center phase (Φ_L = 180°), via 2D pseudospectral Fourier
> collocation (n1=2, n2=20; on-grid residual RMS 1.23e-4, off-grid invariance
> ~2.1e-3).

Mandatory qualifiers wherever the claim is made:
- The MODEL class is Baresi, Owen & Scheeres's (AAS 23-201, 2023; ISSFD 2024)
  — cite both; the Tri-Circular/Laplace-locked idea is theirs, and they
  computed the first N=5 tori (L1/L2 Lyapunov-family substitutes), Floquet
  stability, manifolds, and Europa↔Ganymede transfers in it.
- Distinct from TCP on exactly TWO evidence-backed axes: (i) orbit family
  substituted (mean-motion-resonant, Kumar-class, vs L1/L2 Lyapunov planar);
  (ii) rate model (physical Ganymede synodic rate vs both rates idealized to
  exact rationals — TCP's own periodicity construction requires the
  idealization). Method (pseudospectral PDE vs stroboscopic map+DFT) is a
  difference but NOT independently novelty-bearing. Phase is NOT a
  differentiator (TCP already used 180°).
- Say "exterior Jupiter-Europa 3:4" or plain "Jupiter-Europa 3:4"; NEVER
  "interior" (see verdict).
- No stability, manifold, or transfer computation exists for OUR object yet
  (TCP has all three for theirs) — do not imply otherwise.
- Never: "first N=5 CRNBP torus", "first torus in the Laplace-locked model",
  or any framing `#722` §5's not-to-claim list already bars.
- Novelty is conditional on the literature searched to date (`#721`'s sweep +
  `#722`'s full read + this pass's follow-on search); the mandatory
  `search/literature_check.py` gate still applies before any catalogue
  writeback (NOT run here — writeback is explicitly a separate, future task).

## 5. Items for the coordinating session

1. **Correction to propagate:** `#722` §4/§5's "interior" wording (and the
   `#724` dispatch example language) is factually wrong — exterior orbit, per
   Kumar 2021's own abstract and our own seed geometry. Any digest amendment
   should fix `2026-07-27-722-...-digest.md` §4/§5 (owned by the coordinating
   session, not edited here).
2. **Acquisition candidate:** Kumar, Anderson & de la Llave, Acta Astronautica
   211 (2023) (S0094576523002813) — CCR4BP resonant-tori transfers; not in
   corpus; direct positive-control material for any future N=5
   connection/manifold task.
3. `#721`'s recommended TCP Lyapunov-substitute positive control (reproduce
   their Fig. 8 torus with our corrector at their Table-1 idealized rates)
   remains unexecuted and would further strengthen any writeback.
4. Verification artifacts this pass: `scripts/verify_724_rerun_continuation.py`
   (the §2 reproduction, foreground, ~30 s).
