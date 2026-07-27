# #721 — Adversarial verification of #720's N=5 CRNBP torus (2026-07-27)

**Task:** independent, skeptical re-verification of `#717`'s N=5 CRNBP zero-coupling
claim and `#720`'s "smooth clean mu_Io continuation, no disqualifying literature"
torus result, explicitly modeled on the `#701`→`#702` precedent (do not re-run the
same code and nod along; re-derive, cross-check with a different method, and try to
BREAK the novelty claim). Verification code:
`scripts/verify_721_crnbp_adversarial_checks.py` (all runs foreground, 2026-07-27,
outputs quoted below verbatim).

## VERDICT

**Split, and the split matters:**

1. **The COMPUTATION is CONFIRMED.** `#717`'s zero-coupling claim is independently
   re-derived from raw Newton and holds to machine precision, and `#720`'s torus is a
   genuinely converged invariant torus of its stated model — verified with an
   independent integrator pair (Radau vs DOP853) under the full coupling-included
   EOM, at exactly the same fidelity as the already-validated `#690` baseline.
2. **The RESULT AS FRAMED is NOT CONFIRMED**, on two independent grounds:
   - **The novelty claim is REFUTED by published prior art.** Baresi, Owen &
     Scheeres (AAS 23-257, 2023) and Owen, Baresi & Scheeres (ISSFD 2024) already
     computed **two-dimensional quasi-periodic invariant tori, their Floquet
     stability, and their stable/unstable manifolds in exactly this Laplace-locked
     time-periodic Jupiter-Io-Europa-Ganymede restricted five-body model** (their
     "Tri-Circular Problem", TCP) — including in the **Jupiter-Europa frame** —
     plus Europa↔Ganymede manifold-to-manifold transfer design. Gilliam's own
     thesis (this project's `#717` positive-control source) describes this work
     explicitly in §2.2.5.1 and states "The CRNBP in the case of five bodies is
     essentially equal to this dynamical model."
   - **The delivered torus sits at a Laplace-forbidden perturber phase.** The repo
     default `theta_io0 = theta_gan0 = 0` fixes the Laplace argument
     `theta_io + 2*theta_gan ≡ 0°`; the physical system librates about **180°**
     with ~0.064° amplitude. The all-aligned configuration `#720` computed is one
     the real Jovian trio can never occupy — the exact wrong-phase-anchoring bug
     class (`#701`/`#702`'s `ref_vec`) this task was dispatched to hunt.
     Correctable: the continuation re-run at the physical phase (`theta_io0 = pi`)
     converges equally well (see §3).

**Gate consequence:** the V0-V5 vetting/writeback chain contemplated on a CONFIRMED
verdict is NOT authorized on the present framing. A corrected, much narrower framing
(see §5) could be re-adjudicated later.

## 1. Independent re-derivation of the zero-coupling claim — CONFIRMED

Method (`eom` subcommand; no reuse of `#717`'s derivation, tests, or code paths):
build the exact relative acceleration of the massless particle w.r.t. the M1-M2
barycenter directly from Newton's law, as
`a_rel = a_sc − (M1·a1 + M2·a2)/(M1+M2)`. On this route the extra-body-to-extra-body
coupling **never appears at all** — each primary's acceleration is linear in each
extra body's mass. Compared at 200 random 3D configurations against (a) Negri &
Prado's Eq. 8 bracket form WITH the inner-sum coupling and (b) the superposition
form (coupling dropped, k=1,2 indirect kept):

```
raw-Newton vs Eq.8(+coupling) vs superposition, 2 extras: max rel deviation = 2.009e-16
raw-Newton vs Eq.8(+coupling) vs superposition, 3 extras: max rel deviation = 3.114e-16
```

All three agree to machine precision, for two AND three extra bodies. The deeper
reason (sharper than `#717`'s antisymmetry argument, and confirming it): Negri &
Prado's coupling term is an artifact of their Eq. 4-6 bookkeeping route
(`s_CM12 = −Σ Mj·sj/(M1+M2)` via total-momentum conservation), which re-expresses
the primaries' accelerations through the extras' accelerations. The physical
relative acceleration of a massless spacecraft contains **no** mu_j·mu_k product
terms, period; Eq. 8 is exactly equal to it, so the coupling sum is identically a
zero-sum re-arrangement. `#717`'s claim, its Eq. 7-9 transcription (verified again
directly against the source PDF text this pass), and `#720`'s omission of the term
from the corrector residual are all correct.

Also quantified (not assumed): `crnbp_eom`'s inherited idealisations vs a faithful
Eq. 9/10 transcription (barycentric `−mu2` perturber-position shift + full
`mu1`/`mu2` k=1,2 indirect terms):

```
crnbp_eom vs faithful Eq.9/10: max |delta a| = 3.024e-08  (Io forcing scale ~4.6e-04)
```

~4,000× below the torus residual floor; the documented "inherited, not new"
idealisation stance in `crnbp.py`'s docstring is accurate.

## 2. Independent convergence cross-check — CONFIRMED (in-model)

Rebuilt `#720`'s full pipeline from scratch (3:4 resonant orbit → `#690` CCR4BP
torus at n1=2/n2=20 → mu_io=0 CRNBP seed → 8-step mu_Io continuation). Reproduces
`#720`'s numbers: seed `residual_rms = 1.2210e-4`, final (physical
mu_Io = 4.70434e-5) `1.2392e-4`, smooth and monotonic, rotation number stable
(0.4964720 → 0.4964683).

Then the genuinely independent check (`flow` subcommand, the `ghost_guard`
pattern): sample random torus angles, propagate the analytic torus state through
the FULL coupling-included `crnbp_eom` with **both DOP853 and Radau** at
rtol=atol=1e-12, over horizons up to a **full forcing period** (50× longer than
the corrector's own 0.02-period closure check), and measure the gap to the torus's
own predicted point `u(θ1+ω1t, θ2+ω2t)`:

| torus | FD invariance RMS (off-grid) | drift @0.05P (med) | @0.25P | @1.00P | Radau−DOP853 @1P |
|---|---|---|---|---|---|
| mu_io=0 seed (=`#690` baseline) | 2.141e-3 | 1.756e-3 | 5.614e-3 | 2.947e-2 | 3.0e-12 |
| physical mu_Io, θ_io0=0 (`#720`) | 2.146e-3 | 1.749e-3 | 5.860e-3 | 3.333e-2 | 3.3e-12 |
| physical mu_Io, θ_io0=π (physical phase) | 2.146e-3 | 1.756e-3 | 5.811e-3 | 2.737e-2 | 3.0e-12 |

Readings:
- The physical-mu_Io torus is statistically **indistinguishable from the validated
  `#690` baseline** at every horizon — Io's forcing is genuinely absorbed by the
  corrector; this is not a numerical ghost. Radau/DOP853 agree to ~4e-12 (nondim),
  so none of the drift is integrator artifact.
- **Honest fidelity caveat (pre-existing, inherited from `#690`, not an Io or N=5
  effect):** the invariance residual at random OFF-collocation-grid points is
  ~2.15e-3 RMS — ~17× the reported on-grid `residual_rms` — and free-flight drift
  from the Fourier surface reaches ~3e-2 nondim (~20,000 km at 1 DU = 671,100 km)
  over one full forcing period. That is the truncation floor of the n2=20
  representation of this eccentric 3:4 orbit, identical in the mu_io=0 baseline.
  Any future writeback should quote the off-grid figure, not only the on-grid one.

## 3. The Laplace-phase finding (the `#702`-class bug of this arc)

`lambda_Io − 3·lambda_Eu + 2·lambda_Gan` librates about **180°** (amplitude
~0.064°; Sinclair 1975, Paita et al. 2018). In Europa-synodic phases
(`psi_j = lambda_j − lambda_Eu`) this is `theta_io(t) + 2·theta_gan(t) ≡ 180°`,
which the projected lock `omega_io = −2·omega_gan` conserves exactly
(`d/dt = 0`, verified). Three independent corroborations that 180° (not 0°) is
physical:
1. Gilliam thesis §3.x: the resonance "prohibits [Io, Europa, Ganymede] from
   aligning on the same side of Jupiter"; "Europa and Ganymede at the 0° position
   only when Io is at the 180° position."
2. Baresi/Owen/Scheeres AAS 2023, Table 1, Jupiter-Europa column: `ϕ10 = π`
   (Io), `ϕ30 = 0` (Ganymede) — exactly `(π, 0)`.
3. The classical libration-center literature cited above.

The repo default (`jupiter_europa_io_ganymede_default`, inherited by every `#720`
test and the delivered torus) uses `(0, 0)` — Laplace argument 0°, the antipode of
the libration center, a configuration the real trio never occupies. Because
`theta_io0 + 2·theta_gan0` is invariant under time shift AND under the corrector's
θ1-origin freedom, this is a genuinely different periodic forcing, not a
relabeling. Measured effect: the two phase variants' converged coefficient sets
differ by 2.0e-3 (L2), the same order as Io's entire effect on the torus
(1.6–2.3e-3 vs the seed) — i.e., **the phase matters at exactly the order of the
new physics being claimed**.

`#717` documented the free choice as safe for its equilibrium control (true —
Gilliam's own Ch. IV/V "test case 1" is the theoretical all-aligned case, which is
what the Tables 5-6 reproduction needed). `#714` §1(c) and shortlist item 2
explicitly required pinning/sweeping `theta_io0` for the torus task; `#720`
skipped that step. Mitigation, verified this pass: the continuation re-run at
`theta_io0 = π` converges identically well (final `residual_rms = 1.2343e-4`,
smooth, monotone, same rotation number to 1e-8) — so the correction is cheap, but
the delivered object is not the physical-configuration one.

## 4. Adversarial literature clearance — FAILED (novelty refuted)

What was searched (fresh, deliberately trying to disprove novelty): Baresi/Owen/
Scheeres author-trail (the corpus already held Owen-Baresi 2024's knot-theory
paper — same authors, one hop away); "tri-circular" + Jupiter-Io-Europa-Ganymede;
"exploiting the Laplace resonance" trajectories; Kumar/Anderson/de la Llave 5-body
extensions; "five-body invariant tori Jupiter Laplace"; GTOC Galilean-tour
material (GTOC6 is patched-conic flyby mapping — different object class); plus a
line-by-line read of Gilliam's related-work §2.2.5.1 and of the Negri & Prado
introduction.

**Disqualifying prior art found (primary):**
- **Baresi, N., Owen, D., Scheeres, D. J., "Exploiting the Laplace Resonance for
  Designing Trajectories in the Jupiter-Io-Europa-Ganymede System," AAS 23-257,
  AAS/AIAA Astrodynamics Specialist Conference, Big Sky MT, Aug 2023** (open
  access via Surrey Open Research). Introduces the **Tri-Circular Problem (TCP)**:
  Jupiter-centered synodic frame with the other two Laplace moons on concentric
  circles, time-periodic with period `2πk` (k=4/2/1 for Jupiter-Io/-Europa/
  -Ganymede frames; their Table 1 for the Jupiter-Europa frame: `mu_Io=4.7043e-5`,
  `a_Io=0.6287`, `a_Gan=1.5949` — numerically our system to ~1e-3, with rates
  idealized to exact −0.5/1 where we keep the physical 0.5036 magnitude).
  Computes: periodic dynamical substitutes of L1/L2; **2D quasi-periodic invariant
  tori as substitutes of CRTBP L1/L2 Lyapunov planar orbits in the Jupiter-Europa,
  Jupiter-Ganymede, and Jupiter-Io frames** (stroboscopic invariant-curve/DFT
  collocation, GMOS lineage — a different method from `#720`'s PDE pseudospectral,
  same object class); Floquet stability of the tori; stable/unstable manifolds;
  Europa↔Ganymede transfer identification.
- **Owen, D., Baresi, N., Scheeres, D. J., "Transfer Trajectory Design in the
  Jupiter-Io-Europa-Ganymede Tri-circular Problem," 29th ISSFD, Darmstadt, 2024.**
  Extends the above to low-energy Europa↔Ganymede transfers via manifold
  intersections (ΔV ~500 m/s, TOF ~30 d) — substantially pre-empting `#714`
  shortlist item 3 (connection survival under N=5 forcing) as well.
- **Gilliam thesis §2.2.5.1** (the `#717` control source itself): describes the
  TCR5BP, cites both papers ([27],[28]), and states the N=5 CRNBP "is essentially
  equal to this dynamical model," so "literature results produced for this model
  can be replicated in the more general case of the CRNBP with ease."

**Framing sanity-check (dispatch item 4):** Gilliam's oft-quoted admission reads,
in full: "it is not currently known if any periodic trajectories exist in the
CRNBP **beyond those few cases leveraging resonance**" (§4.5). The carve-out IS
the Laplace-locked case — the one configuration `#714`/`#720` targeted, and the
one where the TCP papers already operate. `#714`'s premises "no literature-grade
dynamical control exists or can exist yet for any N=5 CRNBP torus" and "the
frontier is genuinely unclaimed" are therefore **factually wrong**, and its
"whether the Laplace-commensurability observation is itself unremarked in the
literature is a lit-check question" resolves to: it is the founding idea, and the
title, of the AAS 2023 paper. `#720`'s live WebSearch missed all of this; the
Gilliam digest (`2026-07-26-digest-gilliam-bettinger-2024-crnbp-jovian.md`) never
mentioned §2.2.5.1/TCR5BP — the disqualifying reference sat inside the project's
own source PDF.

**What may remain narrowly novel (NOT adjudicated here, do not claim without a
fresh gate):** the TCP papers compute Lyapunov-family torus substitutes; `#720`'s
object is a substitute of an **interior 3:4 Europa-resonant (Kumar-class) orbit**,
which does not appear in their figures; the PDE-pseudospectral (rather than
stroboscopic-collocation) route in this model may also be methodologically fresh.
That is an incremental family extension inside a published framework — "first
3:4-resonant-family torus substitute in the TCP/N=5 CRNBP, computed at the
physical (non-idealized) Ganymede synodic rate and physical Laplace phase" — not
"first N=5 CRNBP invariant torus." A positive-control opportunity falls out for
free: reproduce a TCP Lyapunov-substitute torus (their Fig. 8) with our corrector,
which would give the N=5 lane the literature-anchored dynamical control `#714`
believed could not exist.

## 5. Code-review notes (adversarial pass over `crnbp.py` / `variational_crnbp_torus.py`)

- Coupling-term sign, per-j isolation test design, `_perturber_second_deriv_block`,
  Coriolis signs, `theta1 → t` clock mapping (incl. the `omega_gan < 0`
  convention), Io's slaved second-harmonic phase (`+2·theta1` at
  `period_multiple=1`), gauge rows, and the residual/Jacobian block structure all
  check out against both the source PDFs and independent FD/flow evidence.
  `#720`'s n1=1→n1=2 representation-capacity finding is real (Io's forcing is pure
  second-harmonic in θ1) and correctly diagnosed.
- The one real defect found is the **epoch-phase anchoring** of §3 (model-level,
  not algebraic). Everything downstream of `jupiter_europa_io_ganymede_default()`
  inherits it.
- Minor reporting hazard, worth a docstring note someday: `residual_rms` is an
  on-collocation-grid quantity; at these n2 values the off-grid residual is ~17×
  larger (§2). Not a bug — standard collocation behavior — but "residual 1.2e-4"
  undersells the truncation floor if quoted alone.

## 6. Recommended follow-ups (for the coordinating session; none started here)

1. Acquire both TCP papers into `cyclers_pdf` + digest + `CORPUS_INDEX.md`
   (AAS 23-257 PDF is openly downloadable from Surrey Open Research; ISSFD 2024
   likewise). Amend the Gilliam digest with §2.2.5.1.
2. Fix `jupiter_europa_io_ganymede_default()` to pin
   `theta_io0 = pi − 2·theta_gan0` (document the 0.064° libration amplitude as the
   validity window), keeping the all-aligned option available for the Gilliam
   Tables 5-6 control (which deliberately uses the theoretical aligned case).
   Re-stamp `#720`'s tests/torus at the physical phase — this pass already
   verified the continuation converges there.
3. If the narrow claim of §4 is worth pursuing: run the TCP Lyapunov-substitute
   positive control first, then re-adjudicate novelty with the TCP papers cited.
4. Process lesson (memory-worthy): the literature gate must include the
   related-work sections of the sources already in hand, and an author-trail
   sweep of corpus authors — the disqualifying citation was one hop from two
   different in-corpus documents and was missed by both `#711`'s digest and
   `#720`'s live search.
