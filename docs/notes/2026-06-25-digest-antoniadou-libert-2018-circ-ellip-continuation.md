# Digest: Antoniadou & Libert 2018 — circular→elliptic continuation of resonant POs (ER3BP)

Single-paper digest. Read 26/26 pages of the PDF on 2026-06-25 AET.
Acquired to give the ER3BP discovery campaigns (#432 / #435 / #436) a
PUBLISHED reference for the question: *do e>0-only elliptic families exist
with no circular limit?*

## 1. Header — AUTHORSHIP CORRECTION

- **Title (verbatim)**: *Origin and continuation of 3/2, 5/2, 3/1, 4/1 and 5/1
  resonant periodic orbits in the circular and elliptic restricted three-body
  problem*
- **Authors (verbatim, title page)**: **Kyriaki I. Antoniadou and Anne-Sophie
  Libert**, NaXys, Department of Mathematics, University of Namur, 8 Rempart de
  la Vierge, 5000 Namur, Belgium (kyriaki.antoniadou@unamur.be)
- **THE FILENAME GUESS "voyatzis-antoniadou" IS WRONG.** Voyatzis (George
  Voyatzis, Thessaloniki) is NOT an author of this paper. He is a co-author on
  several *cited* prior works (Hadjidemetriou-Voyatzis 2000; Voyatzis-Kotoulas
  2005; Antoniadou-Voyatzis 2013/2016; Voyatzis et al. 2009/2018). The correct
  short cite is **Antoniadou & Libert (2018)**. Recommend renaming the corpus
  PDF / index entry accordingly.
- **Venue**: "Accepted for publication in Celestial Mechanics and Dynamical
  Astronomy" (stated on title page).
- **DOI**: 10.1007/s10569-018-9834-8 (printed as the authenticated-version link
  on the title page)
- **arXiv**: 1805.00288 (v3, 22 May 2018)
- **Length**: 26 pages text (front matter through conclusions) + references.

## 2. What the paper actually is

A systematic study of **symmetric resonant periodic orbits (POs)** of an
interior-MMR planetary configuration — star `P0`, inner massless body `P1`,
outer planet `P2` (`m2 = 0.001`, `mu = m2/(m0+m2) = 0.001`) — in BOTH the
**circular** restricted three-body problem (CRTBP) and the **elliptic** RTBP
(ERTBP). The five resonances studied are **3/2, 5/2, 3/1, 4/1, 5/1** (all
*interior* MMRs, `p > p+q` indexing with the giant at `a2 = 1`). The paper:
1. computes the **circular family** of POs and locates where elliptic resonant
   families bifurcate off it (Sec. 3);
2. continues those families into the ERTBP, parametrised by the eccentricities
   `e1` (inner body) and `e2` (the giant `P2`), to high `e` (up to ~1) (Sec. 5);
3. classifies linear stability + maps regular vs chaotic phase-space domains
   with DS-maps (de-trended Fast Lyapunov Indicator, DFLI) (Sec. 4);
4. reports the **isolated** elliptic families (not bifurcating from circular).

The motivation (Sec. 1) is exoplanetary: identifying stable phase-space
neighbourhoods that could host terrestrial companions in single-giant systems;
it is a follow-on to Antoniadou & Libert (2018, the 93-system DS-map paper).

## 3. THE CRUX — do e>0-only families exist with no circular limit?

**YES.** The paper reports two distinct origins for ERTBP resonant families,
and one of them is genuinely elliptic-only:

### 3.1 Families that DO continue from the circular problem (the majority)

The dominant mechanism. From the circular family of POs (Fig. 2, the curve of
`(mean-motion ratio, e1≈0)` circular orbits), at every rational `(p+q)/p` where
`T = 2π(a1^{-3/2} − 1)^{-1}`, **two branches of resonant elliptic families
bifurcate** — one stable, one unstable — by the **Poincaré–Birkhoff fixed-point
theorem** (Sec. 3, p. 6, verbatim):

> "At rational values of the mean-motion ratio, namely `(p+q)/p`, where the
> periodic orbits are of period `T = 2π(a1^{-3/2} − 1)^{-1}` (when `a2 = 1.0`),
> there bifurcate two branches of resonant families of symmetric elliptic
> periodic orbits in the CRTBP; one being stable (depicted by blue colour) and
> one unstable (depicted by red), due to the Poincaré-Birkhoff fixed point
> theorem. One corresponds to the location of `P1` at pericentre (denoted
> hereinafter by `I`) and one to its location at the apocentre (denoted `II`)."

These elliptic families **emanate from the circular family** when `e1 > 0` and
are continued to the ERTBP (`e2` increasing from 0). This is the standard
circular→elliptic bifurcation.

### 3.2 Families that DO NOT continue from the circular problem — ISOLATED families

The directly relevant finding. The paper repeatedly computes **isolated
families** — explicitly defined as **"not continued from bifurcation points"**
(Abstract) / **"not linked with bifurcation points"** (Sec. 3 end) — that exist
in the ERTBP and have **no circular-family origin**. Verbatim, Abstract:

> "Moreover, new isolated (not continued from bifurcation points) families are
> computed in the elliptic restricted problem. The majority of the new families
> mainly consist of stable circular periodic orbits at high eccentricities."

[NB: "stable circular periodic orbits at high eccentricities" is the Abstract's
own phrasing; "circular" here is loose — Sec. 5 makes clear these are the
**high-`e` stable elliptic** members of each isolated family.]

Sec. 3 (p. 7), verbatim, stating prior precedent and the new contribution:

> "It has been shown that there also exist families that do not bifurcate from
> periodic orbits and are isolated. Antoniadou and Voyatzis (2013) computed
> spatial isolated symmetric families in 2/1 MMR in 3D-GTBP resulting by
> foldings of the spatial families as the mass of the inner body increased from
> zero. Voyatzis and Hadjidemetriou (2005); Voyatzis et al. (2009); Antoniadou
> and Voyatzis (2016) presented planar asymmetric isolated families in 2/1 MMR
> in GTBP resulting by collision bifurcations as the planetary mass ratio
> varied. In this paper, we present some new isolated families of symmetric
> periodic orbits in the ERTBP."

**An isolated family is reported for EVERY one of the five MMRs studied**, in
each case at high eccentricity in a specified symmetric configuration:
- **3/2** (Sec. 5.1, p. 9): "an isolated family which consists of highly
  eccentric stable periodic orbits for both `P1` and `P2` and belongs to the
  configuration `(0,π)`."
- **5/2** (Sec. 5.2, p. 13): "an isolated family which consists of highly
  eccentric stable periodic orbits for both the planet and the small body and
  belongs to the configuration `(0,π)`." (also a broad-stable `(π,0)` family
  with `e1 > 0.96`).
- **3/1** (Sec. 5.3, p. 17): "in the configuration `(π,0)`, we computed a new
  isolated family which possesses only stable periodic orbits where both the
  small body and the planet are highly eccentric."
- **4/1** (Sec. 5.4, p. 20): "in the configuration `(π,0)`, we computed a new
  isolated family which possesses only stable periodic orbits where both the
  small body and the planet are highly eccentric."
- **5/1** (Sec. 5.5, p. 23): "in the configuration `(π,0)`, we computed a new
  isolated family which possesses only stable periodic orbits where both the
  small body and the planet are highly eccentric."

Conclusions (Sec. 6, p. 26) confirm: "We also computed some isolated families,
in particular at high eccentricities, for each MMR."

**VERDICT for our purposes: the literature DOES support the existence of
genuine high-`e` elliptic families with no circular limit** — they are the
"isolated" families, present in all five resonances, consistently in a single
symmetric configuration `(0,π)` or `(π,0)`, consistently stable, consistently
at high `e` for both bodies. This is exactly the "e>0-only family" object our
#432/#435/#436 campaigns hunted for and returned NEGATIVE on.

### 3.3 IMPORTANT scope caveat — this is the SUN-PLANET (RTBP) frame, not Earth–Moon

The paper's massless body is the **inner terrestrial planet**; the perturber is
the **outer giant planet** with `mu = m2/(m0+m2) = 0.001` (≈ Jupiter/Sun). The
"eccentricity" that is continued is the **giant planet's `e2`** (ERTBP primary
eccentricity) and the inner body's `e1`. So the isolated families live in the
**Sun–giant-planet** ER3BP at `mu ≈ 1e-3` — directly comparable to our **#435
high-`e` Sun-planet** campaign, and to **#432** (continue known cyclers into
e>0) if those cyclers are Sun-planet. It is NOT the Earth–Moon `mu ≈ 0.0121`
regime; the isolated-family mechanism is mu-dependent (see §6 below). Our #411
EM-libration work (commit 3f51739: "EM-libration seed stays EM-bounded") is a
different mu and frame and is not addressed here.

## 4. The circular→elliptic BIFURCATION mechanism (for the #432/#435 Floquet monitor and #436 branch-switcher)

This is the precise structure our Floquet monitor should be watching for.

- **Where elliptic families branch off** (Sec. 3, p. 6, Fig. 2): along the
  **circular family** (the 1-parameter curve of circular POs, `e1 ≈ e2 ≈ 0`,
  swept by mean-motion ratio). At each resonant `(p+q)/p`:
  - **First-order MMRs (`q = 1`; e.g. 2/1, 3/2)**: the circular family **breaks
    and a GAP forms** in the perturbed (`mu ≠ 0`) case. The two resonant
    elliptic branches (`I` at pericentre, `II` at apocentre) emanate from the
    circular family across that gap (Fig. 2; Sec. 3 p. 6). For 3/2 the gap is
    at `x = a1 = (3/2)^{-2/3} ≈ 0.763143` (Sec. 5.1).
  - **Second-order MMRs (`q = 2`; e.g. 3/1, 5/3)**: the circular family does
    NOT break — instead it develops a **small UNSTABLE segment** (red in Fig. 2,
    where `e1 ≈ 0`). At the **ends of that unstable segment** there are
    bifurcation points to **doubly-symmetric** ERTBP families (this is the
    paper's **Scheme II**, period `T = 2 T0`). (Sec. 3 p. 6; Sec. 5.3 for 3/1.)
  - **Higher-order (`q > 2`; 4/1 `q=3`, 5/1 `q=4`)**: Scheme-II bifurcation
    points generate orbits described `q` times, `T = q T0` (Sec. 3 p. 7).

- **Two distinct continuation Schemes** (Sec. 3, pp. 6–7):
  - **Scheme I — from the CRTBP**: take a *resonant* PO of the CRTBP (a member
    of branch `I`/`II`) and continue it mono-parametrically into the ERTBP,
    increasing `e2` from 0, when its period equals `T = k T0` or `T = k T0/m`
    (`m` = multiplicity). Bifurcation points labelled `B^{(p+q)/p}_{F,#}`.
    *This is the workhorse — most ERTBP families in the paper come from Scheme I.*
  - **Scheme II — from the circular family directly**: at the ends of the
    `q = 2` unstable circular segments, doubly-symmetric POs emanate **directly
    into the ERTBP** with `T = 2 T0`. The NUMBER of POs generated (2 or 4)
    depends on the perturbation magnitude (the masses).

- **The bifurcation diagnostic the paper uses** (Figs. 4b/7b/10b/12b/14b):
  bifurcation points are located where the **CRTBP family period `T` crosses
  the resonant commensurate value** (`T0 = 4π` for 3/2; `4π/3` for 5/2;
  `2π = 2 T0` for 3/1; `2π/3` for 4/1; `2π = 4 T0` for 5/1). This is a
  **period-commensurability** detector, not (directly) a Floquet/monodromy
  eigenvalue crossing. *This matters for our monitor design* — see §5.

## 5. Continuation METHOD vs our #436 stack (the #437 fold-aware question)

The paper does NOT give algorithmic pseudocode for the corrector, but the
method is the standard **Hadjidemetriou-school predictor–corrector
continuation** that the Thessaloniki/Namur group uses throughout the cited
corpus. What we can extract and contrast:

- **State / unknowns**: ERTBP symmetric POs are pinned by the perpendicular-
  crossing conditions `x'(T)=x'(0)`, `x(T)=x(0)`, `ẏ(T)=ẏ(0)` with
  `ẋ(0)=ẏ... =0` at the section (Eqs. 3, 7, p. 4–5). The free element `x'(0)`
  encodes `e2` (the giant's eccentricity): `x'(0) = a2(1−e2)` at pericentre or
  `a2(1+e2)` at apocentre (Eq. 7, p. 5). So **continuation is in `e2`** (or in
  `e1` along a CRTBP family) — a *physical* parameter, with the period locked to
  `T = k T0` by the resonance. This is mono-parametric natural-parameter
  continuation.
- **Stability bookkeeping**: linear stability via the **monodromy-matrix
  eigenvalues in reciprocal pairs**; stability indices `b1 = λ1+λ2`,
  `b2 = λ3+λ4`; stable iff `|b_i| < 2` (Sec. 4, p. 7). They explicitly track
  the 6 instability types of Broucke (1969)/Marchal (1990). **This is exactly
  the Floquet content our monitor needs** — the paper's blue/red colouring of
  every family IS a continuously-tracked stability index along the continuation.
- **Fold/turning-point handling**: the figures (e.g. Fig. 7b 5/2 branch `II`
  showing an S-fold in `(T−T0, e1)`; Fig. 12b 4/1 showing folded `II`) show the
  families **turning back in the natural parameter** — i.e. the continuation
  must be **fold-aware** (pseudo-arclength), because pure natural-parameter
  `e2`-stepping would stall at these folds. The paper does not name
  pseudo-arclength but the diagrams are unambiguous: their continuator survives
  folds (multiple `e` values per period, period turning points).

**Comparison to #436's secant + 2-var symmetric corrector (which flipped
classification on step count):**
- The paper's method tracks the **full reciprocal-pair monodromy spectrum**,
  not a single 2-variable residual. The stability index `b_i` is computed from
  the monodromy matrix at every continuation step — a flip would be visible as a
  `b_i` crossing ±2, not an artifact of step count.
- The **period is LOCKED to the resonance** (`T = k T0`) rather than free. Our
  #436 corrector's instability under step-count change is consistent with the
  classic symptom of a natural-parameter continuator hitting a **fold it cannot
  see** (it overshoots/undershoots the turning point depending on step size).
- **Implication for #437 (fold-aware pseudo-arclength)**: the paper's families
  CONTAIN the exact folds (5/2 `II`, 4/1 `II`, the `e1=0` configuration changes
  where a family passes through `(e1,e2)=(0,0)` and flips configuration, e.g.
  3/2 `(0,π)→(π,π)` at `e1=0`, Sec. 5.1) that defeat a naive secant continuator.
  A robust reproduction REQUIRES: (a) pseudo-arclength parametrisation so the
  step survives `de/ds → ∞` turning points; (b) continuous monodromy-eigenvalue
  tracking (not a 2-var residual sign) for stability; (c) period locked to
  `k T0` as a constraint, with `e2` free. This is a direct argument FOR the #437
  fold-aware build — the published families are unreachable without it.

## 6. Resonances, e-ranges, prograde/retrograde, stability tabulated

- **System**: interior MMRs only, `mu = 0.001`, `a2 = 1`. Inner massless body
  = terrestrial planet (subscript 1), outer giant = perturber (subscript 2).
  All POs are **symmetric** (the paper excludes asymmetric POs here; asymmetric
  2/1 families are in the companion Antoniadou-Libert 93-system paper).
- **Resonances**: 3/2 (`q=1`), 5/2 (`q=3`), 3/1 (`q=2`), 4/1 (`q=3`), 5/1
  (`q=4`). **5/1 is studied for the first time in any RTBP** (previously only in
  the general TBP by Michtchenko et al. 2006); Sec. 1, p. 2, and Sec. 6.
- **Four symmetric configurations** per resonance, by the libration centres of
  `(θ_i, θ_j)` ∈ {`(0,0)`, `(0,π)`, `(π,0)`, `(π,π)`}, with apsidal difference
  `Δϖ = 0` or `π` (Sec. 2, p. 4; sign convention: negative `e1`/`e2` plotted
  when the first/second bracket argument librates about `π`). Resonant-angle
  pairs: `(θ1,θ2)` for odd `q` (3/2, 5/2, 4/1); `(θ3,θ1)` for 3/1 (`q=2`,
  Eq. 5); `(θ4,θ1)` for 5/1 (`q=4`, Eq. 6).
- **Eccentricity range**: families continued over **`e1, e2 ∈ [0,1]`** — the
  whole eccentricity range, explicitly to HIGH `e` (the paper's headline
  contribution: extending previously-known families, considered totally
  unstable, to high `e` and finding stable segments; and the new high-`e`
  isolated families). DS-maps computed on `(e1,e2)` and `(ϖ2,e2)` planes
  (Figs. 5,6,8,9,11,13,15) and a representative `e1=0.2` `(a1/a2,e2)` map
  (Fig. 3).
- **Prograde / retrograde**: all interior MMRs here are **prograde** (the
  massless body orbits the star in the same sense as the giant; `θ_i` are the
  standard prograde resonant angles). The paper does NOT treat retrograde
  resonances — it cites Morais/retrograde and Kotoulas-Voyatzis exterior-
  resonance work but its own families are prograde interior MMRs. *So this paper
  does NOT cover the retrograde-resonance leg of our search space.*
- **Stability headline** (Sec. 6, p. 26): "The majority of the new families we
  found in each MMR has highly eccentric stable periodic orbits." Previously
  totally-unstable families were continued and found to have **stable segments
  at high `e`** (e.g. 3/1 family `Ic` stable for `0.75 < e1 < 0.98`, Sec. 5.3;
  5/2 broad-stable `(π,0)` for `e1 > 0.96`, Sec. 5.2; 4/1 `IC` stable `e1>0.74`,
  Sec. 5.4). Structurally this means: **a family that is unstable in its low-`e`
  segment can re-stabilise at high `e`** as eigenvalues re-enter the unit circle
  — the unstable label is energy-/eccentricity-local, not a property of the
  whole family.
- **Per-MMR bifurcation-point inventory** (the `B` points justified in the `b`
  panels): 3/2: `B^{3/2}_{I,1}`, `B^{3/2}_{IIS,1}` (latter newly reported).
  5/2: `B^{5/2}_{I,1}`, `B^{5/2}_{I,2}`. 3/1: `B^{3/1}_{I,1}`, `B^{3/1}_{II,1}`.
  4/1: `B^{4/1}_{I,1}`, `B^{4/1}_{I,2}` (stable family), `B^{4/1}_{II,1}`
  (unstable; newly reported). 5/1: `B^{5/1}_{II,1}` (only one, on the unstable
  family). New, previously-unreported bifurcation points are claimed for 3/2,
  4/1, 5/1.

## 7. FOR OUR USE — re-adjudicating the #432/#435/#436 negatives

### 7.1 Does this paper CONFIRM or CHALLENGE our negatives?

**It CHALLENGES them — partially and importantly.** Our campaigns concluded
NEGATIVE on "a novel e>0-only family with no circular/CR3BP limit." This paper
demonstrates that **such families DO exist in the published literature** — the
**isolated families**, present in all five of 3/2, 5/2, 3/1, 4/1, 5/1, at high
`e` for both bodies, predominantly STABLE, in a fixed symmetric configuration
(`(0,π)` for the first-order/odd cases, `(π,0)` for the others). They are
*by construction* not reachable by continuation from the circular family — they
are not connected to it. So the *existence* premise our negatives were testing
is **affirmed by the literature**, not refuted.

But two qualifications keep our negatives partly intact:
1. **Method-reachability, not non-existence.** Our #432 (continue cyclers into
   e>0) and #436 (branch-switch at the circular→elliptic bifurcation) are
   **continuation-from-circular** methods. The isolated families are
   **definitionally unreachable by those methods** — they have no bifurcation
   link to continue along. So a #432/#436 NEGATIVE is the *expected* result for
   an isolated family and does NOT prove the family absent; it proves the method
   blind to it. The negatives are **correct verdicts about the method's reach,
   not about the dynamics.** (Cf. memory `project_negative_results_registry`:
   "empty" is always conditional on the method.) This is a textbook
   capability-subsumption flag: a new method (direct high-`e` seeding) can
   subsume the prior continuation sweep.
2. **Frame/mu match needed.** The isolated families are Sun–giant ER3BP at
   `mu ≈ 1e-3`. #435 (high-`e` Sun-planet seeds) is the matching campaign and
   its negative is the one most directly challenged. If #435 used
   **direct seeding** (not continuation) and still missed them, the gap is in
   **where** it seeded: the paper pins the isolated families to a SPECIFIC
   high-`e`, SPECIFIC symmetric configuration `(0,π)`/`(π,0)`, with the apsidal
   difference at `0`/`π` and BOTH bodies eccentric. A seed sweep that did not
   target that exact corner (high `e1` AND high `e2`, anti-aligned/aligned
   apsides, the right resonant-angle libration centre) would miss them even
   though they exist. **Recommend re-running #435/#436 as a TARGETED seed at the
   paper's reported configuration**, not a blind continuation.

### 7.2 Concrete re-examination actions

1. **Reframe the #432/#435/#436 negatives in the registry** from "no e>0-only
   family exists" to "no e>0-only family reachable by circular→elliptic
   continuation, AND none found by the prior (untargeted) high-`e` seed sweep."
   The literature (Antoniadou & Libert 2018) shows isolated high-`e` families
   exist; our methods were not constructed to reach them. (See memory
   `feedback_bugfix_invalidates_past_searches` analogue: a method blind-spot,
   not a bug, but the same "re-examine prior negatives" trigger.)
2. **Build the targeted direct-seed test (subsumes #435/#436):** seed at high
   `e1` AND high `e2` (both bodies eccentric, e.g. `e1, e2 ≳ 0.5–0.9`), in the
   symmetric configuration `(0,π)` for 3/2 & 5/2 / `(π,0)` for 3/1, 4/1, 5/1,
   with `Δϖ = π`/`0` per the configuration, and run the corrector with the
   period LOCKED to the resonant `k T0`. This is the only construction that can
   land on an isolated family — it is the paper's own recipe.
3. **#437 (fold-aware pseudo-arclength) is JUSTIFIED and necessary** — but note
   it solves the *connected*-family reachability (continuing branch `I`/`II`
   through folds), NOT the isolated families. The isolated families need
   **direct seeding**, not better continuation. So #437 and a new
   "targeted-isolated-seed" task are **complementary, both needed**:
   - #437 fixes the #436 step-count instability on the *connected* elliptic
     families (the folds in 5/2 `II`, 4/1 `II`, and the `e1=0` configuration
     flips are the concrete failure cases).
   - a new direct-seed task is the only way to reproduce the *isolated*
     families.
4. **Golden-validation target available.** This paper is a clean **same-model
   golden source** for an ER3BP continuator: `mu = 0.001`, `a2 = 1`, the five
   MMRs, the bifurcation-point period values (`T0 = 4π`, `4π/3`, `2π`, `2π/3`,
   `2π`), the specific `x = a1` gap locations (3/2 at `0.763143`, 5/2 at
   `0.542884`, 4/1 at `0.39685`), and the DS-map `a1` values (Fig. 3 caption:
   3/2 `a1=0.7631`, 5/2 `0.5428`, 3/1 `0.4807`, 4/1 `0.3968`, 5/1 `0.3419`).
   These are PUBLISHED expected values (satisfies memory
   `feedback_golden_tests_sourced_only`). Recommend a golden test that
   reproduces a CRTBP→ERTBP bifurcation point at one of these
   `(MMR, a1, T0)` triples before trusting any new ER3BP continuator.
5. **Scope honesty.** This paper does NOT address: retrograde resonances; the
   Earth–Moon mu regime (#411); 3D/spatial families (planar only); asymmetric
   POs (symmetric only). Our negatives in those sub-spaces are NOT touched by
   this paper and remain as they were.

### 7.3 One-line adjudication

The literature (Antoniadou & Libert 2018, CMDA, DOI 10.1007/s10569-018-9834-8)
confirms **genuine high-`e` elliptic-only ("isolated") resonant families exist
with no circular limit** in the Sun–giant ER3BP (`mu≈1e-3`) for all of
3/2, 5/2, 3/1, 4/1, 5/1 — therefore our #432/#435/#436 negatives are
**method-blind-spot negatives, not non-existence proofs**, and should be
re-examined with a TARGETED high-`e1`+high-`e2`, fixed-symmetric-configuration
DIRECT SEED (the isolated families are unreachable by continuation by
construction), with #437 fold-aware pseudo-arclength fixing the *connected*-
family reach as a separate, complementary track.

## 8. Reference appendix — key citations in Antoniadou & Libert 2018

- Antoniadou, K.I. & Libert, A.-S. (2018) — the companion 93-single-giant-system
  DS-map paper; this work is its dynamical-systems backbone (cited throughout;
  2/1 families overplotted in Fig. 16).
- Hadjidemetriou, J.D. (1992, 1993a, 1993b) — prior 3/1, 4/1 resonant-family
  computations the paper extends to high `e` (Secs. 5.3–5.4).
- Ferraz-Mello et al. (1992, 1993, 2006) — averaged-Hamiltonian /
  apsidal-corotation-resonance (ACR) treatments of 3/2, 5/2, 3/1, 4/1.
- Michtchenko et al. (2006) — the ONLY prior study of 5/1 (general TBP);
  this paper is first to do 5/1 in the RTBP.
- Antoniadou & Voyatzis (2013, 2016); Voyatzis & Hadjidemetriou (2005);
  Voyatzis et al. (2009) — prior ISOLATED-family results (spatial 2/1; planar
  asymmetric 2/1 by collision bifurcation) that motivate the isolated-family
  search here.
- Hénon (1973, 1997) — vertical stability and the fundamental symmetry
  `Σ:(t,x,y)→(−t,x,−y)` (Eq. 2).
- Voyatzis (2008) — the DFLI chaos indicator used for the DS-maps (Eq. 9).
- Broucke (1969); Marchal (1990) — the 6 instability types for `|b_i| > 2`.

End of digest.
