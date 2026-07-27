# Digest: Baresi, Owen & Scheeres — the Jupiter-Io-Europa-Ganymede "Tri-Circular
# Problem" (TCP), both papers (`#722`)

**Task:** `#722`, dispatched off `#721`'s adversarial verification of `#720`'s
N=5 CRNBP torus (`docs/notes/2026-07-27-721-n5-crnbp-torus-adversarial-
verification.md`), which found that Baresi/Owen/Scheeres already computed 2D
quasi-periodic invariant tori, Floquet stability, manifolds, and Europa<->
Ganymede transfers in exactly this project's own N=5 Laplace-locked Jupiter-
Io-Europa-Ganymede model. This note reads both papers in full and delivers the
residual-novelty verdict `#721` deferred.

## 0. Identity, acquisition, and one correction to `#721`'s own citation

**Paper 1:**
Baresi, N., Owen, D., Scheeres, D. J., "Exploiting the Laplace Resonance for
Designing Trajectories in the Jupiter-Io-Europa-Ganymede System," AAS/AIAA
Astrodynamics Specialist Conference, Big Sky, MT, USA, Aug 13-17, 2023.
**Paper number printed on the PDF itself is AAS 23-201**, not AAS 23-257 as
`#721`'s report and this task's own dispatch text stated — corrected here
against the primary source (the running header on p.1 of the acquired PDF
reads "AAS 23-201"; the Surrey Open Research landing page's own metadata
citation independently confirms this). No DOI. Acquired open-access from
Surrey Open Research (`openresearch.surrey.ac.uk`, record `99815166602346`),
Author's Accepted Manuscript / "textimage" version, text-layer PDF, 16 pp.
Filed as `cyclers_pdf/papers/baresi-owen-scheeres-2023-exploiting-laplace-
resonance-tri-circular-problem-AAS-23-201.pdf` (private corpus commit
`ad5b0dd`).

**Paper 2:**
Owen, D., Baresi, N., Scheeres, D. J., "Transfer Trajectory Design in the
Jupiter-Io-Europa-Ganymede Tri-circular Problem," 29th International
Symposium on Space Flight Dynamics (ISSFD), Darmstadt, Germany, Apr 22-26,
2024. No DOI. Acquired open-access from Surrey Open Research
(`99877064602346`), Author's Accepted Manuscript, text-layer PDF, 7 pp.
**Note on the acquisition route:** the issfd.org proceedings page's own link
for this title (`ISSFD2024_19-5.pdf`) is mislabeled/stale — it actually serves
Owen & Baresi's separate knot-theory paper (DOI `10.1007/s42064-024-0201-0`,
already in this corpus, `owen-baresi-2024-knot-theory-...pdf`). The correct
paper was retrieved from Surrey's own repository instead; verified by reading
the acquired PDF's own title page. Filed as `cyclers_pdf/papers/owen-baresi-
scheeres-2024-transfer-trajectory-design-tri-circular-problem-issfd2024.pdf`
(private corpus commit `ad5b0dd`).

Both are Author's Accepted Manuscripts (not publisher typeset), clean
text-layer PDFs, `pdftotext -layout` extraction verified clean, read in full
(all 16 + 7 pages).

## 1. The TCP model formulation vs this project's `core/crnbp.py`

**Structurally the SAME restricted N=5 problem, same physical system, but a
DIFFERENT (and importantly, non-general) equations-of-motion presentation.**

TCP's Eq. (1) (identical in both papers, reproduced verbatim in the ISSFD
paper's Eq. 1):

```
r'' + 2*omega x v + omega x (omega x r) = -muJ/r^3 * r
    - sum_{i=1..3} mu_i * [ r_i/r_i^3 + r_iJ/r_iJ^3 ]
```

i.e. Jupiter's direct pull plus, for EACH of the three Galilean moons
independently, its own direct + indirect (barycentric-correction) term — a
straight PER-MOON SUPERPOSITION. **There is no pairwise "inner sum" coupling
term between extra bodies anywhere in TCP's Eq. (1)** — contrast this
project's own `core/crnbp.py`, which implements Negri & Prado's fully general
N-body form INCLUDING the extra-body-to-extra-body coupling term (`#717`
Step 0), then separately PROVES (algebraically, numerically, and by
independent re-derivation) that this coupling term sums to exactly zero for
the total acceleration for any N. **TCP's Eq. (1) is therefore not a
different/simpler physics model — it is exactly this project's own
`_perturbers_acceleration` proven-reduced form**, arrived at independently by
a different author group without ever deriving (or needing) the general
coupling term. This is a genuine, useful cross-validation: two independent
derivations (Negri & Prado's general form reduced by proof here vs.
Baresi/Owen/Scheeres's own from-scratch Eq. 1) land on the identical reduced
equation.

**Two concrete formulation differences, both bearing on the residual-novelty
question below:**

1. **Rate idealization.** TCP Table 1 (reproduced in both papers, identical
   numbers) gives, for the Jupiter-Europa frame (their reference pair = this
   project's own base CR3BP pair): `mu_Io = 4.7043e-5`, `a_Io = 0.6287`,
   `n_Io(tilde) = 1` (exact), `phi_Io0 = pi`; `mu_Gan = 7.8047e-5`,
   `a_Gan = 1.5949`, `n_Gan(tilde) = -0.5` (exact), `phi_Gan0 = 0`. TCP LOCKS
   BOTH synodic rates to exact rational values (Io at +1, Ganymede at -0.5 —
   i.e. an idealized EXACT 4:2:1 commensurability). This project's
   `jupiter_europa_io_ganymede_default()` instead keeps Ganymede's rate at its
   PHYSICAL, JPL-registry-derived value (`omega_gan ~ -0.5036`, ~0.7% off the
   idealized -0.5) and Laplace-projects only Io's rate onto EXACTLY
   `-2 * that physical omega_gan` (not onto the idealized -0.5). The masses
   and semi-major axes agree with TCP's Table 1 to ~1e-3 (registry vs. their
   sourced JPL constants); the rates diverge at the ~2e-4 to ~7e-3 level
   depending on which idealization is compared. `#721`'s own §4 already
   flagged this as a candidate residual-novelty axis ("computed at the
   physical (non-idealized) Ganymede synodic rate") — confirmed directly
   against the source here.
2. **Phase.** TCP's `phi_Io0 = pi`, `phi_Gan0 = 0` (Jupiter-Europa frame) is
   EXACTLY the physical libration-center phase this project's `#723`
   correction now defaults to (`theta_io0 = pi - 2*theta_gan0`, giving
   `theta_io0 = pi` at `theta_gan0 = 0`). This is now a THIRD independent
   corroboration (after Sinclair 1975 and Gilliam's own thesis prose) that
   180 deg, not 0 deg, is correct — directly from the source PDF, not a
   secondhand citation.

## 2. Torus computation method, Floquet stability, manifolds (Paper 1, AAS 23-201)

**Method:** two-stage numerical continuation, NOT this project's 2D
pseudospectral PDE-collocation approach:
- **Equilibrium substitutes (periodic orbits replacing L1/L2):** 1D collocation
  (Legendre-polynomial, `N` subintervals x `m`-degree, Eqs. 4-17) grown from
  the CR3BP L1/L2 state as a degenerate ("flat") initial guess, converged by
  Newton's method on the collocation residual + periodicity constraint. Linear
  stability from the monodromy matrix built as a byproduct of the collocation
  Jacobian (Eq. 16-17).
- **Quasi-periodic torus substitutes of CR3BP PERIODIC ORBITS (their genuinely
  novel-to-them contribution):** homotopy continuation in a coupling
  parameter `epsilon` (Eq. 19: `h = g + epsilon*(f - g)`, `g` = base
  Jupiter-single-moon CR3BP vector field, `f` = full TCP field) from a
  CR3BP periodic orbit (`epsilon=0`) to a full 2D quasi-periodic invariant
  torus of the TCP (`epsilon=1`), using the STROBOSCOPIC-MAP + DFT
  "GMOS"-lineage method (Olikara 2016 PhD thesis; Baresi, Olikara & Scheeres
  2018 JAS; the SAME Olikara-Scheeres/Baresi method lineage this project's own
  `variational_qbcp_torus`/`variational_ccr4bp_torus`/`variational_crnbp_torus`
  reuses, but Baresi/Owen/Scheeres apply it via multiple-shooting +
  stroboscopic invariant-circle matching + DFT rotation operator `R_{-rho}`,
  NOT the 2D tensor-product Fourier pseudospectral collocation `#617`/`#690`/
  `#720` use). Quasi-periodicity constraint (their Eq. 20/Owen's Eq. 6):
  `R_{-rho}{phi_T(X_i)} - {X_i} = 0`. Rotation number is FIXED (not solved
  for) at `rho = T/k`, `k` = least common multiple of the synodic rates —
  the torus's poloidal frequency is a KNOWN CONSTANT of the resonance, not a
  free continuation output. This is a structurally different but
  mathematically equivalent object class to this project's own tori.
- **Floquet stability:** Floquet matrix `B = R_{-rho} * Phi(T+t0, t0)`
  (Eq. 23), 6N eigenvalues on 6 concentric circles in the complex plane (the
  standard symplectic-pair structure for Hamiltonian quasi-periodic tori).
  **Result: EVERY torus/equilibrium-substitute the paper computed is
  hyperbolic** — at least one Floquet circle has radius > 1, hence a paired
  circle < 1 by the symplectic reciprocal-pair property. No linearly stable
  case was found among any of the L1/L2 dynamical substitutes or Lyapunov
  quasi-periodic tori tested.
- **Orbit families computed:** L1/L2 dynamical substitutes (periodic-orbit
  analogues of the CR3BP equilibria) for all three Jupiter-moon pairs
  (Figs. 3-5); quasi-periodic torus substitutes of the **30th, 50th, 70th
  Lyapunov PLANAR orbit family members** in the Jupiter-Europa (L1 and L2),
  Jupiter-Ganymede (L1), and Jupiter-Io (L1) systems (Fig. 8a-d) — i.e. ALL
  computed orbit families in this paper are members of the L1/L2 LYAPUNOV
  PLANAR family (libration-point-centered), indexed by a family-member number
  (larger member = larger amplitude/period along the family), NOT interior
  mean-motion-resonant orbits.
- **No digit-grade energy/Jacobi-constant table.** The TCP is non-autonomous
  (time-periodic), so it has NO conserved Jacobi integral in the classical
  CR3BP sense — same structural loss this project's own CCR4BP/CRNBP already
  documents. No torus-amplitude numeric table is given (Figs. 8/9/10/11 are
  plots only); Table 1 is the only digit-grade content (system constants, not
  torus output).
- **Manifold/transfer scoping:** exterior unstable manifolds from the L1/L2
  EQUILIBRIUM substitutes do NOT intersect in the Jupiter-Centered Inertial
  frame (Fig. 7) — the periodic perturbation spreads them apart over time.
  Manifolds from the quasi-periodic TORUS substitutes DO overlap in
  configuration space (Fig. 11: 50th Jupiter-Ganymede unstable manifold vs.
  50th Jupiter-Europa stable manifold) — this motivates, but does not itself
  execute, the transfer design that becomes the ISSFD 2024 follow-up.
  Explicit future-work statement: "searching through this database of
  stable/unstable manifold trajectories in order to identify minimum-fuel
  transfer opportunities... between different Galilean moon systems."

## 3. Europa<->Ganymede transfer results (Paper 2, ISSFD 2024)

**Method:** patches the UNSTABLE manifold of a Ganymede quasi-periodic torus
(a QPO — "quasi-periodic orbit" in their terminology) to the STABLE manifold
of a Europa QPO, WITHIN the same TCP dynamical system (unlike Koon et al.'s
classic multi-moon-orbiter approach, which patches manifolds across two
DIFFERENT autonomous CR3BP systems). Because the TCP is time-periodic,
patching requires matching not just position/velocity but ALSO the moons'
relative orientation at the patch time — states are extracted from each
manifold only at four discrete synodic phases (`t~=0, pi/2, pi, 3pi/2`), then
a nearest-neighbor search (position match within 1000 km in the
Jupiter-Centered Inertial frame) finds candidate patch-point pairs, refined by
nonlinear nonlinear-optimization differential correction (continuity
constraint threshold `1e-6`) to a continuous impulsive-DV transfer.

**Orbit families used (their own explicit statement):** "QPOs at Ganymede and
Europa were quasi-periodic continuations of the **70th and 75th members of
their Lyapunov families**, respectively" — again the L1/L2 LYAPUNOV PLANAR
family, member-indexed, NOT an interior mean-motion-resonant family.

**Numeric results (digit-grade, both example transfers reported in full):**
- Transfer A: Ganymede L1 QPO -> Europa L2 QPO, **DV = 839 m/s, TOF = 28.6
  days** (Fig. 7).
- Transfer B: Ganymede L1 QPO -> Europa L2 QPO, **DV = 618 m/s, TOF = 50.1
  days** (Fig. 6).
- Broader set (Fig. 5, DV-vs-TOF scatter): DV as low as ~500 m/s (abstract),
  clustering at integer multiples of Ganymede's orbital period (~7.16 days,
  matched to within an hour) — attributed to resonant apoapsis-reduction,
  echoing Ross et al. 2004's finding in the patched-3BP approach.
- Baseline comparison: two-body Ganymede->Europa Hohmann transfer ~2,800 m/s.
  Koon et al. (1999/2002) patched-3BP manifold transfers: ~1,452 m/s (TOF <1
  month) down to ~470 m/s (TOF >=1 year) but ending in LOW-EUROPA ORBIT
  (requiring a large injection burn, 96% of total DV in their case), not a
  libration/quasi-periodic orbit.
- **Scope: planar (2D) only.** Explicit stated future work: "extension of
  this method into the spatial case... periodic orbits beyond Lyapunov,
  include those with vertical components."

## 4. Direct comparison against `#720`'s own delivered torus

`#720`'s object (per `#721`'s independently-reproduced numbers): a 2D
pseudospectral torus (`variational_crnbp_torus.py`, `n1=2`/`n2=20` Fourier
representation) built by `mu_Io`-continuation from `#690`'s already-validated
Jupiter-Europa **3:4 mean-motion-resonant** CCR4BP torus — itself a
continuation of **Kumar, Anderson, de la Llave & Gunter 2021's** (AAS 21-651,
arXiv:2109.14815, already in this corpus) UNSTABLE RESONANT PERIODIC ORBIT of
the base Jupiter-Europa CR3BP (an interior orbit near Europa's 3:4
mean-motion resonance with Ganymede's period, NOT a libration-point family
member), converged at physical `mu_Io = 4.70434e-5`, `residual_rms = 1.2392e-4`.

**Point-by-point:**

| Axis | Baresi/Owen/Scheeres (TCP, both papers) | `#720` |
|---|---|---|
| Model | TCP Eq. (1), proven-identical to this project's coupling-cancelled `crnbp.py` reduction | `core/crnbp.py` (full Negri-Prado form, coupling proven zero) |
| Rate idealization | Both moons' rates locked to exact rational values (Io=+1, Gan=-0.5) | Ganymede kept at PHYSICAL registry rate; only Io Laplace-projected onto it |
| Phase | `(phi_Io0, phi_Gan0) = (pi, 0)` — physical | `theta_io0 = pi` (per `#723`) — SAME physical phase, independently corroborated |
| Orbit family substituted | L1/L2 **Lyapunov planar** family (member-indexed: 30th/50th/70th/75th) | **Kumar et al. 2021's interior 3:4 Europa-Ganymede mean-motion-resonant orbit** |
| Torus method | Stroboscopic map + DFT + multiple-shooting collocation (Olikara/Baresi lineage) | 2D tensor-product Fourier pseudospectral PDE collocation (`#617`/`#690` lineage) |
| Stability/manifolds | Computed (Floquet, all hyperbolic) + manifold transfers (DV/TOF table above) | Not computed (`#720` explicitly out of scope: "No whisker/manifold globalization... later, gated work") |

**No overlap found.** Neither TCP paper computes, mentions, or shows a figure
resembling an interior mean-motion-resonant orbit family (3:4, 3:2, 7:5, or
any other Kumar-class resonance) anywhere in either paper — every single
family-member index cited in both papers (30th/50th/70th/75th) is explicitly
identified as belonging to "their Lyapunov families," the standard
libration-point planar family, which is geometrically and topologically
distinct from an interior resonant orbit near a moon's own path. A
Lyapunov-family member growing in amplitude does not pass through or
approximate the Kumar-class resonant orbit's geometry (the Kumar orbit stays
close to Europa's own path at a 3:4 period ratio with Ganymede; a large
Lyapunov-planar member instead grows outward from the L1/L2 saddle region).
This is a genuine family-class difference, not a labeling artifact.

## 5. Verdict

**`#720`'s specific torus is NOT fully subsumed by Baresi/Owen/Scheeres —
genuine, narrow residual novelty survives, on the SAME grounds `#721`'s §4
already flagged and this reading now confirms directly against both source
papers:**

1. **Orbit family:** `#720` substitutes the Kumar et al. 2021 INTERIOR 3:4
   Europa-Ganymede resonant orbit; neither TCP paper computes any interior
   resonant family — both are confined to the L1/L2 Lyapunov planar family.
   This is the load-bearing difference; it is unambiguous from a full read of
   both papers (every orbit-family reference in both is explicitly "Lyapunov").
2. **Rate model:** `#720`/`crnbp.py` keeps Ganymede's PHYSICAL (non-idealized)
   synodic rate and Laplace-projects Io onto it; TCP idealizes BOTH rates to
   exact rational values. A secondary, smaller difference (~1e-3 level) but a
   genuine one.
3. **Method:** 2D pseudospectral PDE collocation (`#720`) vs. stroboscopic-map
   + DFT + shooting (TCP) — a different, independently-useful numerical route
   to the same object class, not itself novelty-bearing on its own but
   corroborating that this is not a reproduction.

**What IS subsumed, unambiguously, and should NOT be re-claimed:**
- The GENERAL claim "no literature-grade dynamical control exists for any N=5
  CRNBP torus" (`#714`'s premise) — refuted, exactly as `#721` found.
- The Laplace-locked-resonance framing/model itself ("Tri-Circular Problem")
  as a novel idea — it is the TCP papers' own founding contribution, not
  available to claim.
- Any claim of the FIRST quasi-periodic torus, Floquet stability computation,
  or manifold/transfer design in this N=5 model in general — all three exist
  in Baresi/Owen/Scheeres, predating `#720` by roughly three years (2023) and
  one-plus years (2024) respectively.
- A "first N=5 CRNBP torus" framing of ANY kind, generic or resonant-family.

**Best-evidence framing for any future writeback** (not adjudicated or gated
here, per this task's scope — `#721`'s §4 caveat still applies): "first
quasi-periodic torus substitute of an interior mean-motion-resonant
(Kumar-class, Europa 3:4) periodic orbit in the N=5 Jupiter-Io-Europa-
Ganymede Laplace-resonant (Tri-Circular) model, computed via 2D pseudospectral
collocation at the physical (non-rate-idealized) Ganymede synodic rate and the
physical Laplace phase." This is a narrow, incremental extension inside an
already-published framework, not a new frontier — consistent with `#721`'s
own framing and this project's `[[project_novel_findings_status]]` memory
(novel hits are rare; census, not new species, is the expected outcome).

## 6. Recommended follow-ups (not started here, out of `#722`'s scope)

1. Any future writeback of `#720`'s torus must cite BOTH TCP papers and use
   the narrow framing in §5, not the "first N=5 CRNBP torus" framing `#714`
   used.
2. A genuinely strong positive control now exists and is unexploited: this
   project's own pseudospectral corrector could reproduce a TCP Lyapunov-
   substitute torus (e.g. the 50th Jupiter-Europa L2 member, Fig. 8a/9/10) at
   the SAME idealized rates/phase Table 1 specifies, giving a literature-
   digit-grade positive control this project's N=5 lane currently lacks
   entirely (per `#721`'s §4 and the Gilliam digest's own §5 "weaker
   discovery-strategy grounding" finding).
3. Given the manifold/transfer machinery now exists in the literature for the
   SAME N=5 model, `#714`'s shortlist item 3 (connection survival under N=5
   forcing) is now a scoped, positive-control-backed follow-on rather than an
   open question — but building it is a separate, larger task, not part of
   `#722`.

## 7. Citation-mining pass (added scope, dispatched by the coordinating session)

Re-read the introduction/background sections of the two OTHER N=5-adjacent
corpus papers specifically for THEIR OWN citation trails, applying the same
lesson this task is itself a response to (a paper's own related-work section
is a separate search surface from generic keyword search, and self-coined
terminology can hide directly-relevant prior art).

### 7a. Negri & Prado 2022 ("Circular Restricted n-Body Problem," JGCD,
DOI 10.2514/1.G006430, arXiv:2307.10881)

Re-read the paper's Introduction/background citations (not the full paper
again). Full reference-relevant list of everything topically overlapping this
project's N-body/CRNBP/CCR4BP/Jovian-moon search domain, cross-checked against
`CORPUS_INDEX.md`:

- **Iuliano, A., & Gomes, V. M., "The circular restricted (N+1)-body problem
  formulation," Astrophysics and Space Science (2019).** The paper this work
  explicitly supersedes/corrects ("we found and corrected inaccuracies in
  Iuliano's formulation"). **NOT in corpus.** Directly relevant — it is the
  ORIGINAL circular-restricted-(N+1)-body EOM this project's `crnbp.py` line
  descends from (via Negri & Prado's correction), and its own erroneous
  formulation would be a useful negative/error-pattern control given this
  project already had to resolve its own sign-transcription ambiguity on the
  Negri & Prado formula (`#717`'s Step 0). **Flagged as an acquisition
  candidate**, not pursued here.
- **Negri, R. B., & Prado, A. F. B. A., "A study about a Bicircular Model
  Applied to Analyze Perturbations in Trajectories," (2020) / their own
  earlier BCR4BP-correction paper.** This is Negri & Prado's own PRIOR work
  that the current N-body paper generalizes (self-citation, not third-party).
  Not independently checked against corpus by name here — likely
  already-covered ground (this project's `ccr4bp.py`/`bcr4bp.py` lineage
  already traces through the Negri-Prado BCR4BP correction per the existing
  `2026-07-26-digest-negri-prado-2022-crnbp.md`). Low priority.
- **Koon, Lo, Marsden, Ross** (multiple citations, low-energy transfer /
  invariant-manifold lineage) — ALREADY heavily represented in corpus
  (multiple 1999-2006 KLMR papers digested, see CORPUS_INDEX's "Outer-planet"
  section). No gap.
- **Scheeres, D. J.** (general multi-body dynamics citations) — already
  represented via Kumar/Baresi/Olikara-Scheeres lineage papers in corpus. No
  gap identified.

**Net: one clear gap — Iuliano & Gomes 2019 (Astrophys. Space Sci.), the
predecessor (N+1)-body formulation Negri & Prado's own paper corrects.**
Flagged for a future acquisition task; not acquired here (out of `#722`'s
scope per the coordinator's explicit instruction).

### 7b. Aryan & Fitzgerald 2024 ("Four Body Invariant Structures and Chaos
Analysis for Jovian Multi-Moon Ballistic Transfers," AAS 24-103)

Re-read the paper's introduction/background section. Its citation trail is
the PCCFBP (periodically-perturbed CFBP) lineage, largely already anchored:

- **Koon, Lo, Marsden, Ross** low-energy transfer / multi-moon-orbiter papers
  — already in corpus (multiple entries, see CORPUS_INDEX "Outer-planet"
  section, e.g. `ross-koon-lo-marsden-2003-multi-moon-orbiter-AAS-03-143.pdf`,
  `koon-lo-marsden-ross-2002-low-energy-transfer-jovian-moons-contemp-
  math-292.pdf`). No gap.
- **Anderson, R. L., & Lo, M. W.** invariant-manifold/resonant-transition work
  (Jovian-system tour design lineage) — the Kumar-Anderson-de la Llave-Gunter
  2021/2023 CCR4BP papers already in corpus are the direct descendants of
  this lineage (Anderson is a co-author). No new standalone Anderson & Lo
  paper was surfaced as a MISSING distinct citation in Aryan & Fitzgerald's
  own reference list — their citations to this lineage route through papers
  already anchored.
- **Campagnola, S., Buffington, B. B., et al.**, Europa Clipper tour-design
  methods (Tisserand-Poincare graph, endgame problem) — same lineage already
  noted as cited by Baresi/Owen/Scheeres too (§8 below); NOT currently in
  corpus as a standalone paper (only referenced, not filed). Topically
  adjacent (moon-tour graph search, not invariant-torus/CRNBP dynamics) —
  lower priority than the CRNBP-specific gaps, but flagged.
- No self-coined/hidden-terminology risk found in this paper analogous to
  TCP's "Tri-Circular Problem" — its own model name (PCCFBP) is a standard
  extension of the already-well-covered CCR4BP naming convention, and its own
  reference list does not reveal any single overlooked paper as
  consequential as the TCP papers were for Gilliam's thesis.

**Net: no new high-priority CRNBP/N-body gap found in Aryan & Fitzgerald's own
citations** beyond the Campagnola et al. tour-design lineage already flagged
as lower-priority/orthogonal (moon-TOUR combinatorial search, not
torus/invariant-manifold dynamics).

### 7c. Cross-check: did the TCP papers themselves reveal any further gap?

Both TCP papers' own reference lists (§ REFERENCES above, reproduced in full
in the extraction) were checked against `CORPUS_INDEX.md`. Findings:
- **Olikara, Z. P., "Computation of quasi-periodic tori and heteroclinic
  connections in astrodynamics using collocation techniques," PhD thesis,
  U. Colorado Boulder, 2016.** The foundational method thesis both TCP papers
  cite repeatedly (their torus/manifold algorithm's primary source). **NOT in
  corpus.** This is the single most consequential gap found in this whole
  citation-mining pass — it is the METHODS thesis underlying the exact
  stroboscopic-map+DFT torus algorithm both TCP papers (and, one level
  removed, Baresi's later real-ephemeris work) use, a different route to the
  same object class this project's own pseudospectral corrector targets.
  Flagged as a HIGH-PRIORITY acquisition candidate for a future task (likely
  freely available via ProQuest/CU Boulder's institutional repository, not
  checked here — out of `#722`'s scope).
- **Baresi, N., "Spacecraft formation flight on quasi-periodic invariant
  tori," PhD thesis, U. Colorado Boulder, 2017.** Companion methods thesis,
  same lineage. **NOT in corpus.** Same flag, secondary priority to Olikara's.
- **Baresi, N., Olikara, Z. P., & Scheeres, D. J., "Fully numerical methods
  for continuing families of quasi-periodic invariant tori in astrodynamics,"
  JAS 65(2), 2018.** The specific continuation-method paper cited for the
  phase-condition/parametrizing-equation machinery (their Eqs. 7/21). **NOT
  in corpus.** Flagged, same lineage/priority as the two theses above.
- Campagnola & Russell 2010 (Tisserand-Poincare endgame) and Boutonnet &
  Schoenmaekers 2016 (JUICE tour) — same moon-TOUR lineage flagged in §7b,
  orthogonal to CRNBP dynamics, lower priority.

**Summary of all flagged acquisition candidates (NOT pursued — surfaced only,
per explicit scope instruction):**

| Candidate | Why relevant | Priority |
|---|---|---|
| Olikara 2016 PhD thesis (CU Boulder) | Foundational stroboscopic-map torus/heteroclinic method both TCP papers build on | High |
| Baresi, Olikara & Scheeres 2018 (JAS 65:2) | Specific continuation-method paper for the TCP's phase-condition machinery | High |
| Baresi 2017 PhD thesis (CU Boulder) | Companion QP-torus formation-flying thesis, same lineage | Medium |
| Iuliano & Gomes 2019 (Astrophys. Space Sci.) | Predecessor (N+1)-body EOM Negri & Prado's own paper corrects | Medium |
| Campagnola & Russell 2010 / Boutonnet & Schoenmaekers 2016 | Jovian moon-tour combinatorial design lineage, cited by both TCP papers and Aryan & Fitzgerald | Low (orthogonal to torus/CRNBP dynamics) |
