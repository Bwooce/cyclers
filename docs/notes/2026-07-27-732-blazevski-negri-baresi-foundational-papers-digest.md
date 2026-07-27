# #732 — Three user-supplied foundational papers: Blazevski & Ocampo 2012,
# Negri & Prado 2020, Baresi/Olikara/Scheeres 2018

**Task:** `#732`, the top three items (ranked #1, #4, #5 in `#730`'s §2
cluster) of the `#730` citation-mining acquisition backlog, PDFs supplied
directly by the user (page-1-verified by the coordinating session before
dispatch). All three are grounding-chain papers for this project's own
`core/ccr4bp.py` / `core/bcr4bp.py` / `core/crnbp.py` and
`search/variational_qp_torus.py` / `search/variational_crnbp_torus.py`
lineage. Each carries a specific, checkable cross-check question against this
project's own code — answered in full below, not deferred.

## 0. Acquisition, filing, text-layer check

All three filed in the private `cyclers_pdf` corpus (`/Users/bruce/dev/cyclers_pdf/papers/`,
separate repo, never committed to the public `cyclers` repo):

| Paper | Filename | Pages | Text layer |
|---|---|---|---|
| Blazevski & Ocampo 2012 | `blazevski-ocampo-2012-periodic-orbits-ccr4bp-invariant-manifolds-physica-d-241-1158-doi-10.1016-j.physd.2012.03.008.pdf` | 10 | Native, clean (`pdffonts` shows embedded Type 1 fonts; `pdftotext -layout` extraction verified readable including all equations) |
| Negri & Prado 2020 | `negri-prado-2020-generalizing-bicircular-restricted-four-body-problem-jgcd-43-6-1173-doi-10.2514-1.G004848.pdf` | 7 | Native, clean (embedded Type 1C fonts; extraction verified, though the PDF's own two-column layout interleaves oddly under `pdftotext -layout` in a few spots — cross-checked visually where load-bearing) |
| Baresi, Olikara & Scheeres 2018 | `baresi-olikara-scheeres-2018-fully-numerical-methods-continuing-qp-tori-jas-65-157-doi-10.1007-s40295-017-0124-6.pdf` | 26 | Native, clean (embedded Type 1C fonts; extraction verified, all four methods' equations legible) |

No OCR needed for any of the three — all are digitally-typeset originals with
embedded, extractable text layers. All three read in full from the extracted
text (cross-checked against the PDF directly for the load-bearing equations).

---

## 1. Blazevski & Ocampo 2012 — "Periodic orbits in the concentric circular
restricted four-body problem and their invariant manifolds," *Physica D* 241,
1158-1167

### 1.1 What the paper actually does

Computes **periodic orbits** (not tori) in the planar CCR4BP for the
Jupiter-Europa-Ganymede-spacecraft system, in TWO senses:

1. **§5.1 "Stable retrograde orbits":** orbits found by experimentation +
   Newton's method that circulate around a primary (m2 or m3) in a rotating
   frame — easy to find, not derived from any special method.
2. **§3/§5.2 "Unstable orbits near the collinear libration points":** the
   paper's real contribution — an algorithm to construct **unstable
   Lyapunov-LIKE PERIODIC orbits** near the L1/L2 collinear libration points
   of m2 (Europa) and m3 (Ganymede), exploiting the Europa-Ganymede
   near-2:1 commensurability (`ω2 ≈ 2ω3`) to get a good initial guess via
   repeated perpendicular-crossing bisection (Eqs. 9-12, using
   `TL2,m2,lin(m3=0) ≈ TL2,m2,lin(m3≠0)`, i.e. the CRTBP-linearized libration
   period survives approximately into the full 4-body system). §4 computes
   the (linearly unstable) monodromy eigenvalues/eigenvectors and their
   stable/unstable manifolds; §5.3's headline finding is that these
   manifolds do **not** intersect (in contrast to Koon-Marsden-Ross-Lo 1999
   patched-system manifolds, which do and give the low-energy transfer that
   paper is famous for) — a genuine, explicitly-flagged negative result, not
   a construction failure.

**No torus of any kind — quasi-periodic or otherwise — is computed anywhere
in this paper.** Every object is a strictly periodic orbit (period `n·Tp`,
`Tp = 2π/|ω3-ω2|` the model's forcing period) or its 1D stable/unstable
manifold tube.

### 1.2 Cross-check (a): EOM fidelity against `core/ccr4bp.py` — CONFIRM the
concept, REFUTE literal equation identity

The paper's own EOM (inertial frame, Eq. 3; m1-m2 rotating frame, Eq. 5) are
**NOT algebraically identical** to `ccr4bp.py`'s equations, on two specific,
checkable points:

1. **Frame origin.** Blazevski & Ocampo explicitly hold **m1 (Jupiter) fixed
   at the origin** in both the inertial frame ("in an m1-centered inertial
   reference frame m1 is fixed") and the m1-m2 rotating frame ("centered at
   m1"). This is stated as an approximation ("the assumed motion of bodies
   m1, m2 and m3 does not satisfy Newton's equations of motion, but is merely
   an approximation," §2, directly after Eq. 2). `core/ccr4bp.py`, by
   contrast, is built on the standard **barycentric** CR3BP convention
   (`cr3bp.cr3bp_eom`): Jupiter at `(-mu, 0, 0)`, Europa at `(1-mu, 0, 0)`,
   origin = the Jupiter-Europa barycentre — both primaries orbit their
   mutual barycentre exactly (the classical CR3BP two-body reduction is
   EXACT, not an approximation). Fixing m1 at the origin instead of at its
   own barycentric position is a strictly cruder model: it drops the O(mu)
   reflex motion of the primary that the standard CR3BP retains exactly.
2. **Mass normalization.** Blazevski & Ocampo normalize `Gm1 = 1` (Eq. 2) and
   use raw mass ratios `m2, m3` (i.e. `GM_i/GM_1`), not the CR3BP's
   `(1-mu, mu)` barycentric split. `ccr4bp.py`'s indirect Ganymede term
   (`-mu_gan * r_gan/a_gan^3`) is specifically constructed, per its own
   docstring, to compensate for the **Jupiter-Europa BARYCENTRE's**
   acceleration toward Ganymede — a barycentric-frame indirect term with no
   analogue in Blazevski & Ocampo's m1-fixed formulation, where by
   construction there is no barycentre acceleration to compensate for (m1 is
   simply pinned).

**Verdict:** `ccr4bp.py` is **not a literal implementation of Blazevski &
Ocampo's own 2012 equations.** It shares their genuinely founding CONCEPTUAL
contribution — the "concentric circular restricted four-body problem" idea
itself (extra bodies on independent, concentric, coplanar, non-mutually-coupled
circular orbits about the dominant primary, producing a non-autonomous but
time-periodic rotating-frame system) — this paper coined the term and the
idea, and `ccr4bp.py`'s own docstring's framing ("concentric" perturber,
non-autonomous but time-periodic) is exactly this. But the actual EQUATIONS
`ccr4bp.py` implements descend algebraically from the LATER,
barycentric-CR3BP-consistent Simo/Gomez/Jorba/Masdemont-style BCR4BP
formulation (via the Kumar/Anderson/de la Llave and Negri/Prado lineage
already in this project's corpus — see §2 below), not from Blazevski &
Ocampo's own m1-fixed, non-mu-normalized equations. So the answer to "is
`ccr4bp.py` a faithful reduction/generalization of Blazevski & Ocampo's
original formulation" is: **faithful to the founding IDEA, not to the
literal EQUATIONS** — a nuanced but fully evidence-grounded distinction, not
a simple yes/no.

(Minor aside, not pursued further: the paper's own Eq. 15 numeric parameter
table pairs `m2 = 7.79850e-05` with the Europa-radius entry `r12 = 0.627009`
and `m3 = 2.52805e-05` with the Ganymede-radius entry `r13 = 1` — but Fig. 4's
own captions confirm m2 = Europa, m3 = Ganymede, and the REAL JPL mass ratios
are `mu_Europa ≈ 2.53e-5`, `mu_Ganymede ≈ 7.80e-5` — i.e. the printed m2/m3
values appear numerically SWAPPED relative to the real moons they're
captioned against. This could be a genuine 1-line erratum in the 2012 paper,
or an OCR/print-column misread on this pass's part (the table is printed as
two independent unlabelled columns, not row-paired cells) — flagged
respectfully per `[[feedback_respectful_errata_framing]]`, not asserted, and
not load-bearing for any conclusion above.)

### 1.3 Cross-check (b): Laplace-resonance-orbit relevance to `#724`'s N=5
torus novelty claim — the gap SURVIVES, and this paper is actually the
EARLIEST known articulation of the whole N=5 idea

Two separate findings, both strengthening (not weakening) `#724`'s verdict:

1. **Scope mismatch, same pattern already found for the later TCP papers.**
   This paper's Laplace-resonance-based construction method (§3/§5.2) builds
   **unstable PERIODIC orbits near the L1/L2 COLLINEAR LIBRATION POINTS** of
   Europa and Ganymede — not quasi-periodic tori of any kind, and not
   interior mean-motion-resonant orbits away from the libration points. This
   is exactly the same family-class distinction `#722`'s digest already
   established for Baresi/Owen/Scheeres's 2023/2024 TCP papers (libration-point
   Lyapunov-family objects, never a Kumar-class interior/exterior
   mean-motion-resonant orbit): the SAME distinction, now confirmed to hold
   all the way back to the 2012 origin paper of this entire lineage. Nothing
   in Blazevski & Ocampo's own constructed object set resembles, anticipates,
   or bears on `#720`/`#724`'s specific delivered object (a torus substitute
   of the Kumar et al. 2021 exterior Jupiter-Europa 3:4 resonant PERIODIC
   orbit).
2. **This paper's own final Conclusions paragraph is the historical origin
   of the "N=5 Laplace-locked" idea itself**, predating Kumar et al. 2021 and
   Baresi/Owen/Scheeres 2023 by roughly a decade: "we exploited the
   Europa-Ganymede resonance to obtain the periodic orbits in this paper. As
   is well-known, this resonance is part of an even more interesting
   resonance, namely the Laplace resonance involving Io-Europa-Ganymede...
   If one assumes, for instance, that ωI = 2ωE then the five body
   Jupiter-Io-Europa-Ganymede-spacecraft system will still be periodic in
   rotating frames and their may be periodic orbits, though if we use the
   exact value of ωI then the equations will be quasi-periodic in any
   rotating frame and hence there are no periodic orbits. **It would be
   interesting to see if their are periodic orbits in the five-body system
   under the assumption that ωI = 2ωE.**" This is stated purely as
   speculative future work, never executed in this paper (no N=5
   computation appears anywhere in the body) — and even hypothetically, it
   speculates about PERIODIC orbits (not quasi-periodic tori), still
   confined to the same libration-point-adjacent object class the paper
   already established at N=4.

**Net effect on `#724`:** this paper neither computes nor anticipates, in any
executed or hypothesized form, `#720`/`#724`'s Kumar-class resonant-orbit
torus substitute. If anything it PUSHES BACK the search horizon for the "has
anyone ever done this" question to 2012 with the same negative result,
strengthening rather than threatening the narrow novelty claim `#722`/`#724`
already defended. No change to `#724`'s verdict language is needed; this is
an additional, deeper-history data point supporting it, worth citing in any
future writeback as the earliest on-point negative-control precedent.

---

## 2. Negri & Prado 2020 — "Generalizing the Bicircular Restricted Four-Body
Problem," *JGCD* 43(6), 1173-1179

### 2.1 What the paper actually does

Presents THREE variants of the BCR4BP, differing only in how the third body
`M3` (e.g. the Sun in Earth-Moon-Sun) is assumed to move relative to the
`M1-M2` primary pair:

- **§II.A "Binary case"** (the "usual"/classical BCR4BP): `M3` is assumed to
  orbit the `M1-M2` BARYCENTRE (`CM12`) on a circle. This is the standard
  Simo/Gomez/Jorba/Masdemont formulation. Final EOM, Eq. 7 (canonical units,
  `μ = M2/(M1+M2)`, `μ3 = M3/(M1+M2)`, `R3` = distance CM12-to-M3):

  ```
  x'' - 2y' - x = -(1-μ)(x+μ)/r1³ - μ(x-1+μ)/r2³ - μ3(x - R3 cosψ)/r3³ - μ3 cosψ/R3²
  y'' + 2x' - y = -(1-μ)y/r1³ - μy/r2³ - μ3(y - R3 sinψ)/r3³ - μ3 sinψ/R3²
  ```

  The indirect term is `-μ3·(R3 cosψ, R3 sinψ)/R3³ = -μ3·r_M3/R3³` — i.e.
  **exactly** the classical direct+indirect Sun/third-body term.
- **§II.B "Nonbinary case":** for systems where `M3` is NOT `>>` the M1-M2
  separation (Sun-Jupiter-Earth, Sun-Jupiter-Saturn, triple asteroids), the
  paper instead assumes `M3` orbits `M1` alone (physically more reasonable
  when `M3` is close), but this "comes at the expense of neglecting indirect
  effects of the smaller primary" — the resulting indirect term
  (`-μ3 cosψ/R3²`, Eq. 12a) has NO `(1-μ)` factor, because it assumes
  `sCM12 ∝ M1·p3/|p3|³` only.
- **§II.C "General case":** restores the dropped `M2` correction via an
  approximate two-body-plus-perturbation treatment (`M3` orbits `M1` for the
  circular-orbit IMPOSITION, but the acceleration of `CM12` toward `M3`
  includes BOTH `M1`'s and `M2`'s pull, Eq. 13-14), giving Eq. 15's indirect
  term `-μ3(1-μ)cosψ/R3² + μμ3·R3cosψ/(R3²+1-2R3cosψ)^{3/2}` — an O(μ)
  correction on top of the nonbinary case's indirect term. **Eq. 15 reduces
  EXACTLY to Eq. 7 (the binary case) when `R3 >> 1`** — the paper says this
  explicitly.

**On "corrects Huang 1960":** the `#730` master list's one-line
characterization ("corrects Huang 1960's BCR4BP indirect term") is **not
directly supported by this paper's own text** — Huang 1960 (ref [9]) and
Cronin et al. 1964 (ref [10]) are cited only as prior "closer full
derivations of the BCR4BP" that "none... go through a detailed derivation";
the paper never claims Huang's own formula is wrong, only that the classical
BCR4BP (which Huang and Cronin also used) is a poor APPROXIMATION for
nonbinary systems where `R3` is not `>>1`. Flagging this per
`[[feedback_ground_citations_against_content]]` — the "corrects Huang"
framing appears to be an inherited imprecision from an earlier note, not a
claim this paper itself makes. (The Iuliano correction described in the
already-corpus Negri & Prado 2022 CRNBP digest is a separate, later,
textually-confirmed claim — about a DIFFERENT paper, Iuliano 2016/2019 — not
this one.)

### 2.2 Cross-check: does `core/bcr4bp.py` use the correct term? YES —
already the paper's own validated "binary case," and the physically correct
choice for Sun-Earth-Moon

`bcr4bp.py`'s `_sun_acceleration`:

```python
ax = -system.mu_sun * dx / d3 - system.mu_sun * sx / a_sun3
```

with `d3 = |r - r_sun|³`, `a_sun3 = a_sun_nondim³` — this is **term-for-term
identical** to Negri & Prado's Eq. 7a "binary case" (`-μ3(x-R3cosψ)/r3³ -
μ3cosψ/R3²`, noting `μ3 cosψ/R3² ≡ μ3·(R3cosψ)/R3³ = μ3·sx/a_sun³`). So
`bcr4bp.py` implements the classical **binary case**, not the nonbinary or
general case.

Is that the right choice? **Yes, unambiguously, by this paper's own stated
criterion.** `bcr4bp.py`'s system is Sun-Earth-Moon with `a_sun_nondim =
388.81` EM-distance units — an extreme `R3 >> 1` regime. The paper's own
§IV text states plainly that for such systems "the binary case is still
better than the nonbinary" and that the general/binary difference "vanishes
as R3 becomes large" (Eq. 15 "reduce[s] to Eqs. (7) ... when applied to a
binary system of primaries (R3 >> 1)"). Quantitatively: the general case's
extra correction relative to the binary case is an `O(μ_EM)` relative
adjustment (`μ_EM ≈ 0.01215`) applied to an indirect term that is already
`O(μ_sun/a_sun³) ≈ O(3.29e5/388.8³) ≈ O(5.6e-6)` of the leading gravity term
— i.e. an absolute correction of order `1e-7`, several orders below the
`O(eps²)` coherent-QBCP correction (the alpha_i Fourier tables) that
`bcr4bp.py`'s own docstring already documents as its dominant, currently-unmodeled
fidelity gap. **`bcr4bp.py` is not using an outdated or incorrect term — it
already matches this paper's own validated best-choice formulation for
exactly the regime it targets.** (Note: the "nonbinary"/"general" case
distinction this paper actually introduces would matter for a hypothetical
future Sun-Jupiter-Europa-style BCR4BP, or for `ccr4bp.py`'s Ganymede term at
a much larger mass ratio than Europa's tiny `mu ≈ 2.5e-5` — not for the
dormant Sun-Earth-Moon module actually checked here.)

### 2.3 Citation-mining pass (background/related-work section)

Re-read the Introduction in full. Every citation topically overlapping this
project's search domain, cross-checked against `CORPUS_INDEX.md`:

- **Huang, S.-S., "Very Restricted Four-Body Problem," NASA TN D-501 (1960)**
  — cited as one of the "closer full derivations" of the BCR4BP. **Not in
  corpus.** A 1960 NASA Technical Note; likely freely available via NASA
  NTRS (not checked this pass). Low-medium priority — background/historical
  value only, since (per §2.1 above) this paper does not actually claim
  Huang's formula is erroneous, only that the general BCR4BP class (Huang's
  included) is a poor fit for nonbinary systems.
- **Cronin, J., Richards, P. B. & Russell, L. H., "Some Periodic Solutions of
  a Four-Body Problem," *Icarus* 3(5-6), 423-428 (1964),** DOI
  `10.1016/0019-1035(64)90003-X` — the other "closer full derivation" cited
  alongside Huang. **Not in corpus.** Same low-medium priority/rationale as
  Huang above.
- **Gabern, F. & Jorba, A., "A Restricted Four-Body Model for the Dynamics
  Near the Lagrangian Points of the Sun-Jupiter System," *DCDS-B* 1(2),
  143-182 (2001)** — cited [11] as a prior nonbinary-regime BCR4BP
  application (Sun-Jupiter-Earth). **Not in corpus.** Topically adjacent
  (nonbinary BCR4BP dynamics near libration points) but orthogonal to this
  project's own CCR4BP/CRNBP Jovian-moon focus (this is Sun-Jupiter-scale,
  not Jupiter-moon-scale). Low priority.
- **Negri, R. B., Sukhanov, A. & de Almeida Prado, A. F. B., "Lunar Gravity
  Assists Using Patched-Conics Approximation, Three and Four Body Problems,"
  *Adv. Space Res.* 64(1), 42-63 (2019)** — self-citation, own prior work.
  Not independently checked; low priority (cislunar patched-conics, already
  well-covered ground per this project's own Earth-Moon lineage).
- **Assadian, N. & Pourtakdoust, S. H., "On the Quasi-Equilibria of the
  BiElliptic Four-Body Problem with Non-Coplanar Motion of Primaries," *Acta
  Astronaut.* 66(1-2), 45-58 (2010)** — cited [14] as the source of the
  triple-integral-of-motion comparison method (Eqs. 16-20) this paper reuses.
  **Not in corpus.** Method-adjacent (a different four-body variant,
  non-coplanar/bielliptic) but not directly on this project's CCR4BP/CRNBP
  model-definition critical path. Low priority.
- All other references ([1-8], [12-13], [15-16]) are either already-covered
  low-thrust/manifold-transfer application papers orthogonal to the
  model-definition question, or (ref [9] Negri-Prado's own earlier BiCircular
  paper) already noted as covered ground by the existing
  `2026-07-26-digest-negri-prado-2022-crnbp.md`'s own citation-mining pass.

**Net: two clear historical-interest gaps (Huang 1960, Cronin et al. 1964),
neither urgent** — flagged for the `#730` backlog, not pursued further here.

---

## 3. Baresi, Olikara & Scheeres 2018 — "Fully Numerical Methods for
Continuing Families of Quasi-Periodic Invariant Tori in Astrodynamics," *J.
Astronaut. Sci.* 65, 157-182

### 3.1 What the paper actually does

A head-to-head numerical comparison of **four** methods for computing 2D
quasi-periodic invariant tori of an AUTONOMOUS Hamiltonian system (`ẋ =
f(x)`, Eq. 1), tested first in the Earth ECEF co-rotating frame (accuracy
test against the analytic two-body solution) and then on a family of
Earth-Moon PCRTBP distant retrograde orbits (DROs):

1. **PDE(CD)** — solves the torus invariance PDE (`ω0·∂u/∂θ0 + ω1·∂u/∂θ1 =
   f(u)`, Eq. 2) over a full 2D `(θ0,θ1)` grid, approximating the angle
   derivatives with 2nd-order **central differences** (Schilder, Osinga &
   Vogt 2005's original method).
2. **PDE(DFT)** — the same invariance PDE, but the angle derivatives are
   computed via the **Discrete Fourier Transform** of the gridded torus
   values (i.e. spectral/pseudospectral differentiation) instead of finite
   differences.
3. **KKG** — a two-point boundary value problem (TPBVP) for an **invariant
   curve of the Poincaré map** at a surface of section (Kolemen, Kasdin &
   Gurfil 2012), a reduced-dimension (1D-curve, not 2D-grid) shooting method.
4. **GMOS** — a TPBVP for an **invariant curve of the stroboscopic map**
   (Gómez & Mondelo 2001; Olikara & Scheeres 2012), also a reduced-dimension
   shooting method, propagating the invariant-circle points forward one
   period and matching the rotated result via DFT (Eq. 20-23).

**Findings:** in the ECEF accuracy test, GMOS and PDE(DFT) both
significantly outperform PDE(CD) and KKG (Fig. 9, machine-precision-level
errors vs. `~1e-4` for the other two). In the more practical Earth-Moon DRO
torus-family test, **GMOS is found to be both MORE ACCURATE and FASTER than
PDE(DFT)** (Table 1: e.g. `N=25`, GMOS 72.8 s total vs. PDE(DFT) 67.7 s total
but PDE(DFT) spends 49.7 s of that in Newton's-method linear algebra alone,
vs. GMOS's 0.05 s — and PDE(DFT)'s Newton-solve cost scales far worse with
grid refinement, `N0=N1=51` → 3325 s vs. GMOS(MS) 647 s). The paper's own
stated conclusion: **"the GMOS algorithm... [is] our preferred choice for
future investigations of practical astrodynamics problems"** — cited reasons:
better accuracy/speed in their own tests, Floquet/stability information as a
free byproduct of the collocation Jacobian, and a fundamentally
lower-dimensional (curve, not full grid) representation that scales better
to higher-dimensional tori.

**Important scope caveat, explicit in the paper itself:** both of the
paper's own TEST PROBLEMS use a **STABLE** parent periodic orbit (monodromy
eigenvalues on the unit circle — "assume that M admits one pair of complex
conjugate eigenvalues with unitary modulus," the ECEF planar circular orbit
and the Earth-Moon DRO are both well-known STABLE families). The paper does
not test, or even discuss, a torus continuation problem with a violently
UNSTABLE (hyperbolic monodromy) parent orbit.

### 3.2 Cross-check (c): which of the four methods does this project's own
torus corrector resemble, and do the paper's findings validate or challenge
that choice?

**This project's `variational_qp_torus.py`/`variational_ccr4bp_torus.py`/
`variational_crnbp_torus.py` is unambiguously a PDE(DFT)-class method**, not
a shooting/invariant-curve method:

- It solves the SAME invariance PDE (`omega1*du/dtheta1 + omega2*du/dtheta2 =
  f(u)`, verbatim the module's own docstring Eq. (*), identical in form to
  this paper's Eq. 2) on a full 2D `(theta1, theta2)` grid.
- Angle derivatives are computed **spectrally** via a real tensor-product
  Fourier-series representation (`u_c = sum C[c,a,b]*phi_a(theta1)*phi_b(theta2)`,
  analytic term-by-term differentiation) — mathematically the SAME operation
  as this paper's DFT-based derivative approximation (a change of basis
  between gridpoint values and Fourier coefficients is exactly what a DFT
  is; this project's real cos/sin basis and this paper's complex DFT
  coefficients are equivalent representations of the identical
  pseudospectral-differentiation idea).
- Like PDE(DFT) (and unlike GMOS/KKG), this project's corrector never
  integrates the true nonlinear flow inside the search — only in bootstrap
  seed generation and the independent closure check, exactly mirroring this
  paper's own PDE(DFT) vs. GMOS distinction (GMOS's residual IS a forward
  propagation; PDE's is not).

**So: this project's own choice is the paper's LESS-preferred of the two
"good" methods** (PDE(DFT), not GMOS) — a genuine, evidence-backed finding
that partially CHALLENGES the choice, but with an important, checkable
mitigating fact:

**This project switched from GMOS to PDE specifically to escape a pathology
this paper's own test suite never exercises.** `variational_qp_torus.py`'s
own docstring (§"Motivation and the wall this crosses") documents that the
EM L1 quasi-halo torus's parent periodic orbit has a **monodromy spectral
radius of ~1540** (violently unstable/hyperbolic) — the GMOS stroboscopic-map
residual's one-period forward propagation (and every finite-difference
column of its Jacobian) is amplified ~1540x by this instability, causing
catastrophic ill-conditioning that degrades monotonically with torus
amplitude and fails to converge above amplitude ~0.015-0.02. Both of
Baresi/Olikara/Scheeres's own test cases — the ECEF planar circular orbit and
the Earth-Moon DRO — are explicitly **STABLE** parent orbits (their own
Eq. 1 setup assumes "one pair of complex conjugate eigenvalues with unitary
modulus," i.e. an elliptic, not hyperbolic, center); their comparison never
probes a parent orbit anywhere near this project's ~1540x-amplification
regime. **The paper's own preference for GMOS is therefore conditioned on
stable, well-conditioned parent orbits — exactly where GMOS has no
propagation-amplification problem to begin with — and does not directly
speak to the specific unstable-parent-orbit wall this project built the PDE
corrector to cross.**

**Net assessment:**
- **Partially validates:** for STABLE tori (where this project's OWN existing
  GMOS-lineage corrector, `genome/qp_tori.py`, already converges cleanly —
  documented as working below amplitude ~0.01), this paper's finding that
  GMOS is faster AND more accurate than a PDE(DFT)-style approach is
  credible and likely transfers; this project has never benchmarked its own
  PDE corrector against its own GMOS corrector on the same low-amplitude,
  stable torus to confirm this directly — a concrete, cheap follow-up this
  paper's own Table 1/Fig. 9 methodology suggests and this project currently
  lacks.
- **Does not challenge, and is not contradicted by,** the PDE corrector's
  specific reason for existing: crossing a documented shooting-fragility
  wall tied to a violently unstable parent orbit, a regime this paper's own
  comparison never tests. Nothing in this paper suggests GMOS would succeed
  where `#611`/`#612` documented it failing (the ill-conditioning is a
  property of the stroboscopic-map RESIDUAL's forward propagation, which
  this paper's own GMOS description confirms is intrinsic to the method, not
  an implementation artifact).
- **A concrete, better-founded design going forward** (not previously
  articulated this precisely): a HYBRID policy — prefer GMOS
  (`genome/qp_tori.py`) whenever the parent orbit's monodromy spectral
  radius is small enough for it to converge cleanly (this paper's own
  finding: faster, more accurate, gives stability for free), and fall back
  to the PDE/pseudospectral corrector only above whatever spectral-radius
  threshold GMOS's amplification makes impractical — rather than treating
  the two as strictly ordered alternatives. This project's own PDE
  correctors do not currently compute Floquet/stability information at all
  (explicitly out of scope per `#690`/`#720`'s own discipline notes); GMOS
  gives this for free per this paper's finding — a second concrete reason to
  prefer GMOS wherever it converges.
- **Secondary structural note:** this paper's own analysis is scoped to
  AUTONOMOUS systems only (its own Eq. 1, `ẋ=f(x)`, no explicit time
  dependence); it explicitly flags that the same four methods extend to
  non-autonomous systems only by pointing to a companion reference (Baresi &
  Scheeres 2016, IAC). This project's own CCR4BP/CRNBP torus correctors
  already operate on genuinely non-autonomous (time-periodic) systems
  directly — a generalization beyond this 2018 paper's own literal scope,
  consistent with (not contradicting) its comparison.

### 3.3 Citation-mining pass (background/related-work section)

Re-read the Introduction and Methodologies opening in full. Topically
relevant citations, cross-checked against `CORPUS_INDEX.md`:

- **Schilder, F., Osinga, H. M. & Vogt, W., "Continuation of quasi-periodic
  invariant tori," *SIAM J. Appl. Dyn. Syst.* 4(3), 459-488 (2005)** — the
  originating PDE(CD) method paper. **Not in corpus.** Already flagged in
  `#730`'s own reference list (this project's `variational_qp_torus.py`
  docstring already cites it by name) — no new gap, just confirming it is
  still unacquired.
- **Gómez, G. & Mondelo, J. M., "The dynamics around the collinear
  equilibrium points of the RTBP," *Physica D* 157(4), 283-321 (2001)** —
  the original stroboscopic-map GMOS precursor method. **Not in corpus.**
  Foundational to the GMOS lineage this project's own `genome/qp_tori.py`
  already implements a version of. Medium priority.
- **Olikara, Z. P. & Scheeres, D. J., "Numerical method for computing
  quasi-periodic orbits and their stability in the restricted three-body
  problem," *Adv. Astronaut. Sci.* 145, 911-930 (2012)** — the direct GMOS
  algorithm paper this project's `genome/qp_tori.py` is built on. **Not in
  corpus** (only referenced/summarized secondhand so far). Medium-high
  priority — this is the actual algorithm source, not just a lineage
  citation.
- **Kolemen, E., Kasdin, N. J. & Gurfil, P., "Multiple Poincaré sections
  method for finding the quasiperiodic orbits of the restricted three body
  problem," *CMDA* 112(1), 47-74 (2012)** — the KKG method source. **Not in
  corpus.** Lower priority (this paper's own finding is that KKG is the
  weakest of the four methods; useful only as a documented negative
  comparator, not a method this project would adopt).
- **Olikara, Z. P., "Computation of Quasi-periodic Tori and Heteroclinic
  Connections in Astrodynamics Using Collocation Techniques," PhD thesis,
  U. Colorado Boulder (2016)** — cited [22] as the source of "a more
  efficient version of the GMOS algorithm based on collocation." **Already
  independently flagged as the single highest-priority gap in `#730`'s §2
  item 2** (also independently surfaced by `#722`'s own citation-mining
  pass, §7c) — this is a THIRD independent citation trail converging on the
  same thesis. No new action; reinforces the existing HIGH-priority flag.
- **Jorba, A. & Villanueva, J., "On the persistence of lower dimensional
  invariant tori under quasi-periodic perturbations," *J. Nonlinear Sci.*
  7(5), 427-473 (1997)** — cited for the two-parameter-family degeneracy
  argument motivating the extra phase/parametrizing constraints (this
  project's own gauge-fixing rows mirror this exactly). **Not in corpus.**
  Medium priority — direct theoretical grounding for this project's own
  gauge-fixing scheme (`variational_qp_torus.py`'s phase/amplitude anchor
  rows).
- **Haro, A. & de la Llave, R., "A parameterization method for the
  computation of invariant tori and their whiskers in quasi-periodic maps:
  explorations and mechanisms for the breakdown of hyperbolicity," *SIAM J.
  Appl. Dyn. Syst.* 6(1), 142-207 (2007)** — cited [12] as one of the
  excluded analytical/semi-analytical methods ("limited region of
  convergence"). Already covered ground (the Haro/de la Llave parameterization-method
  lineage is heavily represented in corpus via the Kumar-lineage papers);
  no new gap.
- **Jorba, A., "Numerical computation of the normal behaviour of invariant
  curves of n-dimensional maps," *Nonlinearity* 14, 943-976 (2001)** —
  already cited by name in this project's own `variational_qp_torus.py`
  references list. **Not in corpus** (already flagged, `#730`'s own
  reference set does not separately list it but it is the same class of
  gap as Schilder et al. above). Low-medium priority.

**Net: reinforces the existing HIGH-priority Olikara 2016 thesis flag (now a
3rd independent hit), and surfaces one genuinely new, medium-high-priority
gap not previously flagged: Olikara & Scheeres 2012 (AAS 145, "Numerical
method for computing quasi-periodic orbits and their stability in the
restricted three-body problem")** — the actual GMOS algorithm source paper
this project's own `genome/qp_tori.py` builds on, currently held only
secondhand via later citing papers.

---

## 4. Summary for the coordinating session

All three papers filed + digested + citation-mined per
`[[feedback_corpus_document_policy]]`. Direct answers to the three mandatory
cross-check questions:

1. **Blazevski & Ocampo 2012 EOM fidelity:** `ccr4bp.py` faithfully carries
   forward this paper's founding CONCEPT (concentric, coplanar, circular,
   non-mutually-coupled extra perturbers around the dominant primary) but
   is NOT a literal implementation of its equations — those equations fix
   the primary at the origin with a non-barycentric mass normalization,
   while `ccr4bp.py` is built on the standard barycentric CR3BP. **Laplace-
   resonance relevance to `#724`:** none — the paper's own Laplace-resonance
   method builds unstable PERIODIC orbits near LIBRATION POINTS only (no
   torus, no interior resonant orbit), the same distinction already
   established for the later TCP papers, now confirmed all the way back to
   this 2012 origin paper; its own final-paragraph N=5 speculation is
   unexecuted and still confined to periodic (not quasi-periodic) orbits.
   `#724`'s novelty claim is unaffected and, if anything, reinforced by a
   deeper negative-control history.
2. **`bcr4bp.py`'s indirect term:** CORRECT, not outdated. It matches Negri &
   Prado 2020's own "binary case" (the classical Simo/Gomez/Jorba BCR4BP
   term) exactly, and this paper's own analysis confirms the binary case is
   the right choice for `bcr4bp.py`'s actual system (Sun-Earth-Moon,
   `R3=388.8 >> 1`) — the paper's "general case" correction targets a
   different regime (nonbinary systems, small `R3`) that does not apply
   here. Also: the master list's "corrects Huang 1960" framing is not
   directly supported by this paper's own text (Huang is cited as a prior
   derivation, not shown erroneous) — flagged as an inherited imprecision.
3. **Baresi/Olikara/Scheeres 2018 method comparison:** this project's own
   torus corrector (`variational_qp_torus.py`/`variational_crnbp_torus.py`)
   is a **PDE(DFT)**-class method (2D grid, spectral/Fourier
   differentiation, no forward integration in the search) — the paper's
   OWN second-choice method, not its preferred GMOS. The paper's finding
   (GMOS faster+more accurate, gives stability for free) is credible for
   STABLE parent orbits but was never tested against the violently unstable
   (~1540x monodromy amplification) parent-orbit regime this project's PDE
   corrector was specifically built to handle — so the finding partially
   validates a GMOS-preference for the STABLE regime (an untested but
   promising follow-up: benchmark this project's own PDE corrector against
   its own GMOS corrector on a shared low-amplitude torus) without
   invalidating the PDE choice for the unstable-parent-orbit wall it exists
   to cross.

Citation-mining across all three surfaced two genuinely new gaps beyond
what `#730`/`#722` already flagged: **Olikara & Scheeres 2012** (AAS 145,
the direct GMOS algorithm source, medium-high priority) and, at lower
priority, **Huang 1960** / **Cronin, Richards & Russell 1964** (historical
BCR4BP derivations) and **Gabern & Jorba 2001** (Sun-Jupiter nonbinary
BCR4BP application). None acquired here, per this task's scope (digest +
flag, not acquire).
