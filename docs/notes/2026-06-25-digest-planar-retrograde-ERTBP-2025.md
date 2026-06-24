# Digest: Martínez-Cacho, Gil Calvo, Bombardelli & Baresi 2025 — Planar retrograde periodic orbits (QSOs/DROs) in the ER3BP

Single-paper digest. Read all 36 pages of the PDF on 2026-06-25 AET.
Commissioned to test the #435 DRO-into-e>0 negative ("survive, no bifurcation,
no novel family") and the #432/#436 isolated-family question (challenged by
Antoniadou & Libert 2018). READ-and-EXTRACT only — no commits, no shared-file
edits; the controller consolidates CORPUS_INDEX.md.

## 1. Header (confirmed from title page p. 430)

- **Title (verbatim)**: *Planar retrograde periodic orbits in the elliptic
  restricted three-body problem*
- **Authors (verbatim)**: Alicia Martínez-Cacho ᵃ, Daniel Gil Calvo ᵃ,
  Claudio Bombardelli ᵃ, Nicola Baresi ᵇ
  - ᵃ School of Aerospace Engineering, Universidad Politécnica de Madrid (UPM),
    Madrid, 28040, Spain
  - ᵇ Surrey Space Centre, University of Surrey, Guildford, GU2 7XH, United
    Kingdom
  - Corresponding author: A. Martínez-Cacho (alicia.martinez.cacho@upm.es)
- **Venue**: *Acta Astronautica* 229 (2025) 430–465
- **DOI**: 10.1016/j.actaastro.2025.01.006 (PII S0094-5765(25)00008-6)
- **Received / revised / accepted**: 16 July 2024 / 18 December 2024 /
  4 January 2025; available online 12 January 2025
- **Licence**: open access, CC BY 4.0 (© 2025 the authors, Elsevier on behalf
  of IAA)
- **Length**: 36 pages (pp. 430–465)
- **Keywords (verbatim)**: elliptic restricted three body problem; resonant
  orbits; curvilinear coordinates; quasi-satellite orbits; circular Hill three
  body problem

## 2. What the paper is

A comprehensive recent (2025) study of **planar RETROGRADE periodic orbits in
the ER3BP** — explicitly the family class our #435 DRO seeds probe. The paper's
own framing (Abstract, §1): planar retrograde periodic QSOs "beyond the limit of
validity of Hill's approximation are analyzed in detail starting from Hénon's
*f* family and including symmetric and asymmetric solutions up to a multiplicity
of degree seven." The terminology is important and must be mapped to ours:

- **QSO = quasi-satellite orbit = DRO = distant retrograde orbit** — the paper
  states this identity in the very first sentence (§1 p. 430): "Quasi-satellite
  orbits (QSOs), also known as distant retrograde orbits (DROs)..." So this
  paper IS the recent comprehensive ER3BP study of the family our #435 campaign
  continued. The retrograde QSO family is **Hénon's *f* family** (the *f* family
  is "the same as one already studied by Jackson", a Jordan-curve [non-self-
  intersecting] orbit; §1 p. 430).
- **Two top-level classes** (§4, §5):
  - **steady QSO** = non-self-intersecting (Jordan-curve) guiding-centre orbit;
    the "DRO" proper. (§4 p. 434)
  - **swing QSO** = self-intersecting LVLH trajectory; the guiding centre
    itself moves (multi-revolution / "swing" QSOs). (§5 p. 435)

The model is built in a **pulsating curvilinear (Nechville-type) coordinate
formulation** (§2.2–2.3), continued with a **modified Lara–Peláez 2002
predictor–corrector** (§3). The continuation parameter is the mass ratio μ
(Hill → CR3BP) and then the eccentricity e (CR3BP → ER3BP) (§3 p. 433).

## 3. Model, systems, eccentricity ranges

### 3.1 Model (§2)
Full ER3BP in pulsating curvilinear coordinates (ρ, θ, z), true anomaly ν as
independent variable; length unit `LU = a·p/γ` with `p = 1−e²`,
`γ = 1+e cos ν` (§2.2 p. 431). Setting e=0 in Eq. (10) recovers the CR3BP in
curvilinear coords (Eq. 11); Hill's approximation is Eq. (12) (§2.4 p. 432).
Energy parameter `Γ = 1 − C_J/C_J^{L1}` (Hill case) and `Γ_μ` (Eq. p. 440),
so Γ=0 at L1 energy, larger for higher energy. **This is a same-model match to
our CR3BP/ER3BP curvilinear stack** — the Jacobi-equivalent Γ is directly
comparable.

### 3.2 Systems (Table 6, p. 440) — verbatim μ, e, a, T
| System | μ | e | a (km) | T (days) |
|--------|---|---|--------|----------|
| Mars–Phobos (M-Ph) | 1.661025566426e-08 | 0.0151 | 9.3772e3 | 0.31908773969512 |
| **Sun–Mars (S-M)** | 3.227154604242 29e-07 | **0.0934** | 2.279392e8 | 686.9715762940817 |
| Jupiter–Europa (J-Eu) | 2.528582285697873e-05 | 0.009 | 6.709e5 | 3.55043931445928 |
| Saturn–Titan (S-Ti) | 2.366992384915317e-04 | 0.0288 | 1.221865e6 | 15.946220549988796 |

Note: **the high-e system studied here is Sun–Mars at e=0.0934** — exactly one
of our #435 high-e Sun-planet targets (we also did Mercury 0.206, Mars 0.093,
Pluto 0.249). Hill's approximation is also computed (μ-independent) as the
generating model. They do NOT study Mercury or Pluto, so the very-high-e end
(e≳0.2) of our #435 sweep is NOT covered here — only e≤0.0934.

### 3.3 Bounds on the search (§4.2 p. 434, §5.3 p. 436, Eq. 43)
- Steady QSOs: resonance ratio `N21:N32` with `(N32)_max = 10`,
  `(N21/N32)_min = 3/10` → 22 resonant steady QSOs (the 1:1 upper bound is the
  steady-QSO ceiling).
- Swing QSOs: period-multiplying multiplicity `m = N32/Ns` up to `m=7`,
  `(N21)_max = 11`, `(Ns)_max = 6` (Eq. 43 p. 436).

## 4. Do any retrograde families exist ONLY at e>0 (isolated / no circular limit)? — DIRECT #435 / #432 / #436 test

**This is the decisive finding for our purposes, and it CORROBORATES the #435
negative.**

The paper's entire construction is a **bottom-up continuation chain**:

> Hill's approximation (μ=0 limit)  →  CR3BP (continue in μ)  →  ER3BP
> (continue in e from 0 to the system value).  (§3 p. 433; §5.5 p. 439; §6 p. 440)

Every family in the paper is **born in the circular problem** (Hill / CR3BP) and
then carried into e>0. Two explicit statements:

1. Steady QSOs (§4.3 p. 435): "the steady QSOs correspond to planar orbits...
   the out-of-plane motion can be omitted... an additional continuation stage
   is required until the desired resonant period has been reached" — i.e. they
   are continued *from* the two-body / CR3BP resonant orbit into ER3BP. No
   isolated e>0 steady QSO is reported.

2. Swing QSOs (§5.4 p. 437, p. 437 col. 2): "**all the families investigated in
   this work directly bifurcate from the family of steady QSOs in Hill's
   approximation through period-multiplying bifurcations.**" And again (§7
   Conclusions p. 450): the swing families "bifurcate from the family of steady
   QSOs" (period-multiplying in Hill / μ-dependent CR3BP; tangent in the general
   CR3BP for asymmetric ones).

So **no retrograde QSO family in this paper exists only at e>0 with no circular
limit**. Every symmetric and asymmetric retrograde family traces back to the
circular (Hill/CR3BP) steady-QSO (DRO) trunk via a bifurcation that already
exists at e=0. This is the SAME structural picture our #435 sweep reported:
DRO seeds *survive* the e>0 continuation, and the families that branch off do so
from circular-limit bifurcations — **there is no isolated elliptic-only
retrograde family hiding at high e in these systems.** The Antoniadou & Libert
2018 "isolated high-e elliptic family with no circular limit" phenomenon is **NOT
reproduced** for the planar retrograde QSO class in this paper. (Caveat in §6
below: Antoniadou & Libert's isolated families are a different object —
eccentric *mean-motion-resonant* families of the planetary/general ER3BP, not
co-orbital 1:1 retrograde QSOs — so the two results are not in contradiction;
this paper simply confirms the retrograde-QSO class behaves the way #435 found.)

The paper notes (and our #435 should record) that the *existence* of specific
resonant triplets is μ- and e-dependent: some triplets cease to exist as μ
grows (§6.2.7 p. 450 — e.g. the 1:4:1 NA QSO exists only for μ ≤ 4.9e-9, the
2:7:1 NA only for μ ≤ 2.4e-11) and some cannot be closed in ER3BP. But these are
families *terminating* as parameters increase, not families *appearing* with no
circular limit. The direction is loss-of-existence, never spontaneous-elliptic-
birth.

## 5. Bifurcation structure, folds, turning points in e — DIRECT #436 / #437 test

This is the second high-value finding, and it gives #437 (fold-aware follow-on)
a **concrete, sourced, same-model fold to target.**

### 5.1 The fold bifurcation in eccentricity (§6.2.1 p. 441)
> "for the **Sun–Mars** system, a **fold bifurcation was identified at
> e = 0.0324**, marking a maximum in eccentricity for the branch ... This
> explains why the orbit could not be continued up to Sun–Mars eccentricity
> (e = 0.0934). Notably, if this branch is further continued past the fold
> bifurcation back to the CR3BP (e = 0), it connects to a swing QSO discussed in
> the next subsection."

The orbit is the **2:3:1 IDS (interior-doubly-symmetric) swing QSO** in Sun–Mars
(§6.2 p. 440: "the 2:3:1 IDS swing QSO cannot be closed in the ER3BP for the
Sun–Mars system. Remarkably, this orbit connects to the 2:3 steady QSO through a
fold bifurcation"). So the **e-continuation of a retrograde swing QSO turns
around at e=0.0324 and folds back** — a continuation that naively pushes toward
e=0.0934 fails precisely because of a fold/turning-point in e. This is exactly
the **continuation-fragility our #436 negative described**, and it is the
**fold-aware target #437 was set up to chase**: the fold is real, named, located
(e=0.0324, Sun–Mars), and connects two distinct retrograde families (2:3 steady
QSO ↔ a swing QSO) across the turning point.

### 5.2 Bifurcation cascade (§5.4 p. 437–439, Figs 2–6)
- All symmetric swing families bifurcate from steady QSOs via **period-
  multiplying** bifurcations (Hill & μ-dependent CR3BP). Asymmetric families
  bifurcate via **tangent** bifurcations from the symmetric swing families in the
  μ-dependent CR3BP (in Hill they bifurcate directly from steady QSOs) (§5.4
  p. 439 footnote 9; §6.2.7 p. 450).
- m=3 case (Fig. 2 p. 438) is special: **two bifurcation points giving rise to
  two branches sharing one common branch** — EDS orbits on the inner common
  branch, IDS on the two outer branches. m=5, m=7(1) give doubly-symmetric
  orbits from two steady-QSO bifurcation points each; m=4,6,7(2) give both
  symmetric and asymmetric families.
- Multiplicity / tangent-bifurcation behaviour is **μ-dependent**: "As μ
  increases, the tangent bifurcating point diverges further from the period
  multiplying one, limiting the range of periods for which the family exists"
  (§5.5 p. 439). This is the mechanism behind the μ-dependent family
  termination in §6.2.7 (e.g. 1:4:1 NA only for μ≤4.9e-9). **Same fragility
  signature as #436.**
- ER3BP closure failures: many swing QSOs "could not be closed" in ER3BP
  (Tables 7–9 use "c" = collision, blank = not-found / numerical challenge).
  §5.5 p. 439: closure can fail "due to numerical challenges ... or the fact
  that a specific resonance ceases to exist when going from CR3BP to ER3BP ...
  one specific case has been identified as a fold bifurcation." So the fold at
  e=0.0324 is the *worked example* of a generic e-continuation failure mode.

### 5.3 General trend in e and μ (§7 Conclusions p. 450)
> "resonant asymmetric QSOs include resonance triplets that exist only in the
> Hill's approximation model and for four different primaries mass ratio μ (three
> of them) or only for two or three different values of μ (two triplets)."
> "Resonant symmetric QSOs, likewise, may also cease to exist for higher μ."
> "some resonant symmetric QSOs will disappear but also opens the door for the
> emergence of new feasible triplets. A detailed and extensive analysis is
> recommended." — i.e. the authors explicitly leave open whether *new* triplets
> emerge at higher μ/e, but within the studied range every family has a circular
> limit.

## 6. Stability; tabulated ICs / periods / Jacobi-equivalent (same-model golden)

### 6.1 Stability
Stability is tracked via monodromy-matrix eigenvalues (§3 Eq. 20 p. 433); the
geometric classification (IDS/EDS/single-symmetric LWSS/RWSS/LNSS/RNSS,
asymmetric WA/NA) is derived from symmetry + crossing structure (§5.2 p. 436).
Steady QSOs (the *f*-family DRO) are the stable trunk ("planar QSOs are known to
be particularly stable", §1 p. 430); swing QSOs born at period-multiplying
bifurcations inherit/lose stability across the bifurcation. The paper does NOT
publish a per-orbit stability-index table; stability is shown structurally via
the Γ-vs-period bifurcation diagrams (Figs 2–6) and the "Bif" markers.

### 6.2 Tabulated ICs — Appendix B (pp. 454–456+) — SAME-MODEL GOLDEN MATERIAL
Appendix B prints **machine-precision (15-digit) initial conditions** for
steady QSOs in both CR3BP and ER3BP, per system. Initial state form
`s0 = (ρ0, 0, 0, θ'0, 0)` at periapsis (ν=0), x-axis crossing. Examples (verbatim,
Appendix B.1):

- **B.1.1 Hill's approximation** (steady QSOs, `ξ0, η0`): e.g.
  - 1:2 → ξ0 = −0.79704264475925, η0 = 2.19615537549956
  - 2:3 → ξ0 = −1.06588835216964, η0 = 2.52237853252648
  - 9:10 → ξ0 = −1.82614498818691, η0 = 3.79938920930938
  (22 resonances 1:2 … 9:10 tabulated.)
- **B.1.3 Sun–Mars** steady QSOs (CR3BP `ρ0, θ'0` | ER3BP `ρ0, θ'0`), e.g.:
  - 1:2 → CR3BP (−0.00546327973042, 0.01515612806123) |
          ER3BP (−0.00606928578541, 0.01529217473349)
  - 2:3 → CR3BP (−0.00730522383292, 0.01744618027522) | **ER3BP ( — , — )**
    [blank: the 2:3 steady QSO could NOT be closed in ER3BP for Sun–Mars — this
    is the orbit tied to the e=0.0324 fold, §6.2.1]
  - 9:10 → CR3BP (−0.01251796699906, 0.02645723131910) |
           ER3BP (−0.01374285495579, 0.02768366969028)
- **B.1.2 Mars–Phobos, B.1.4 Jupiter–Europa, B.1.5 Saturn–Titan** likewise give
  full CR3BP+ER3BP `ρ0, θ'0` columns for all 22 resonances (1:2 … 9:10).
- **B.2 Swing QSOs** — Hill `ξ0, η0` for the symmetric resonances (1:5:1 IDS,
  2:3:1 IDS, 3:6:2 IDS, … 11:18:6 IDS, etc.), e.g.:
  - 1:5:1 IDS → ξ0 = −0.32382381506633, η0 = 2.24501073437046
  - 2:3:1 IDS → ξ0 = −0.97716364840151, η0 = 2.44850205523003

These are directly ingestible as a **same-model golden** for our curvilinear
ER3BP corrector: pick (system, resonance), seed `(ρ0, 0, 0, θ'0, 0)` at ν=0, and
our integrator should reproduce a periodic orbit of period `T = N32·T_steady`
(Eq. 32, 38). The Sun–Mars 2:3 blank ER3BP entry is the **fold-bifurcation
golden** — our continuation should also fail to close 2:3 IDS at e=0.0934 and
should detect the turning point at e≈0.0324 if #437's fold detection works.
(The Appendix A coordinate transform Eqs A.1–A.7, pp. 453–454, converts these
curvilinear ICs to the LVLH/Cartesian L-frame for cross-checking.)

## 7. For our use: CORROBORATE or CHALLENGE #435?

**CORROBORATES the #435 DRO negative, decisively, with a same-model recent
publication — and supplies the fold #437 should target.**

1. **#435 DRO-into-e>0 ("survive, no bifurcation, no novel family")**:
   CORROBORATED. The paper's retrograde QSO (=DRO) families all have a circular
   limit; DRO/steady-QSO seeds survive e-continuation; the only new families are
   period-multiplying / tangent bifurcations whose roots already exist at e=0.
   No isolated elliptic-only retrograde family is found. Our "survive, no novel
   family" verdict matches the published structure. **One refinement #435 should
   absorb**: "no bifurcation" is too strong as a blanket statement — the steady
   QSO trunk *does* throw period-multiplying bifurcations into swing-QSO families
   (Figs 2–6), but those are circular-limit bifurcations, not e>0-only ones, and
   they are self-intersecting *swing* orbits, not new *steady* DRO families. So
   #435's claim is correct for the steady/Jordan-curve DRO class; the swing
   bifurcations are a known adjacent structure, also circular-rooted. Note also:
   the paper's e-coverage tops out at Sun–Mars e=0.0934; it does NOT reach our
   Mercury (0.206) or Pluto (0.249) — those remain ours alone, but the trend
   (loss-of-existence, never elliptic-birth) extrapolates against finding an
   isolated family there.

2. **#432 / #436 "isolated high-e elliptic family" challenge (Antoniadou &
   Libert 2018)**: this paper does NOT find Antoniadou–Libert-style isolated
   elliptic-only families *in the retrograde co-orbital QSO class*. Important
   scope distinction to record: Antoniadou & Libert's isolated families are
   eccentric mean-motion-resonant families of the general/planetary ER3BP
   (resonances like 2:1, 3:1 in the eccentricity direction), a DIFFERENT object
   from the 1:1 co-orbital retrograde QSO. So this paper neither reproduces nor
   refutes them — it confirms the *retrograde-QSO* class is circular-rooted,
   which is the class #435 actually probed. The #432/#436 isolated-family hunt,
   if continued, should target eccentric MMR families (2:1, 3:1, …) in the
   general ER3BP, NOT the retrograde QSO continuation — this paper shows the QSO
   route is circular-rooted and will keep returning the #435 negative.

3. **#436 continuation-fragility & #437 fold-aware follow-on**: STRONGLY
   SUPPORTED and given a concrete target. The **e=0.0324 fold bifurcation in the
   Sun–Mars 2:3:1 IDS swing QSO** (§6.2.1 p. 441) is exactly the turning-point-
   in-e that #436 hit blindly. #437's fold-aware continuation should (a)
   reproduce this fold as a regression golden, (b) detect the turning point and
   pseudo-arclength around it (the paper used "a standard pseudo-arclength
   continuation method [39]" for the Hill families, §5.4 p. 437), and (c) note
   that past the fold the branch reconnects to the 2:3 steady QSO at e=0 — i.e.
   the fold links two named retrograde families, it does not birth a new one.

### Net recommendation
- Keep the #435 DRO negative; upgrade its provenance from "our sweep" to "our
  sweep + corroborated by Martínez-Cacho et al. 2025 (Acta Astro 229:430), the
  recent comprehensive same-model ER3BP retrograde-QSO study; no isolated
  elliptic-only retrograde family exists in the studied range e≤0.0934".
- Do NOT chase isolated families via the retrograde-QSO route; redirect any
  Antoniadou–Libert-style hunt to eccentric MMR families in the general ER3BP.
- Feed #437 the e=0.0324 Sun–Mars 2:3:1 IDS fold as a sourced same-model golden;
  use Appendix B ICs (steady-QSO `ρ0, θ'0` per system, CR3BP & ER3BP) as the
  corrector regression set, with the Sun–Mars 2:3 blank ER3BP entry as the
  expected fold-failure case.
- Corpus-anchor candidate for `literature_check.py` (recommend, don't implement):
  Martínez-Cacho et al. 2025, retrograde QSO / DRO families, systems
  {Mars–Phobos, Sun–Mars, Jupiter–Europa, Saturn–Titan} + Hill, resonant
  triplets `N21:N32:Ns` up to m=7, planar ER3BP. A discovery candidate that is a
  planar retrograde QSO in one of these systems at a tabulated resonance →
  status="published".

## 8. Key references the paper builds on (for our corpus cross-check)
- Lara & Peláez 2002 — the predictor–corrector continuation algorithm (§3,
  ref [38]).
- Hénon 2003 (and the three-article sequel) — the *f*/*g*/*c* family
  nomenclature and the g3 family of period-multiplied QSOs (§1, refs [16–18]).
- Antoniadou et al. (cited [25], §1 p. 431) — "addressed bifurcation aspects of
  symmetric and asymmetric resonant periodic orbits in the ER3BP" — this is the
  Antoniadou-school reference; worth pulling to settle the #432/#436 isolated-
  family scope precisely.
- Lidov & Vashkov'yak [26,27]; Kogan [28]; Lara [29]; Baresi et al. [30] (the
  "Linear Model" LM, Eq. 33–34) — QSO analytic / averaging heritage.
- Pushparaj et al. [5]; Chen et al. [22] — symmetric/asymmetric QSOs around
  Phobos.
- Broucke & Boggs [23], Broucke [24] — general (non-restricted) elliptic
  three-body planar POs incl. a few QSOs.

End of digest.
