# Digest: Anderson & Kumar 2024 (AAS 24-288)

**Paper:** "A Survey of Oberon Mean Motion Resonant Unstable Orbit Properties and Connections for
Uranian Tours." 2024 AAS/AIAA Astrodynamics Specialist Conference, AAS 24-288, 19 pages.
**Authors:** Bhanu Kumar (Postdoctoral Researcher, Institute for Mathematics, Heidelberg University),
Rodney L. Anderson (Technologist, Jet Propulsion Laboratory, Caltech). [Author order on the title
page is Kumar first, Anderson second — matches the paper's own byline, not the task prompt's
"Anderson, Kumar" ordering.]
**PDF source:** direct S3-hosted conference-proceeding link,
`https://s3.amazonaws.com/amz.xcdsystem.com/A464D031-C624-C138-7D0E208E29BC4EDD_abstract_File24217/FinalPaperUpload_288_0930052440.pdf`.
**Filed:** `anderson-kumar-2024-oberon-mmr-unstable-orbit-survey-aas-24-288.pdf` (private `cyclers_pdf`
repo).
**Acquired/digested:** 2026-07-27 (`#728`).
**Text layer:** native LaTeX-produced text-layer PDF confirmed via `pdffonts` (embedded Type-1
Computer Modern + NimbusRoman fonts, plus non-embedded Helvetica for a table/label) and
`pdftotext -layout` (clean extraction, 887 lines). **No OCR needed.**

**Why in corpus:** flagged by `#598`'s citation-mining pass (`data/OUTSTANDING.md` line ~15280-15281)
as "directly relevant to the project's Uranian-system work (`#312`/`#569` family)." This task (`#728`)
executes that acquisition + the mandated cross-check against `#312`/`#569` AND the more recent
`#701`-`#708` CCR4BP Umbriel-Titania arc (not yet existing when `#598` wrote that flag).

## 1. What the paper actually does

Studies unstable mean-motion-resonant (MMR) periodic orbits of Uranus's outermost large moon,
Oberon, in two models:

1. **Uranus-Oberon planar circular restricted 3-body problem (PCRTBP)** — computes unstable resonant
   periodic-orbit families for the 3:4, 4:5, 5:6 exterior and 4:3, 5:4, 6:5 interior Oberon
   resonances, their stable/unstable manifolds (parameterization method), and heteroclinic
   connections between adjacent-resonance families, characterized as a function of Jacobi constant
   C ∈ [3.00, 3.01] (§"Uranus-Oberon PCRTBP Resonant Orbit Analysis," pp.6-12).
2. **Uranus-Oberon-Titania concentric circular restricted 4-body problem (CCR4BP)** — the same
   Blazevski & Ocampo 2012 [18] CCR4BP model class already used by the Kumar-Anderson-de la Llave
   Jovian-system papers in this corpus (`kumar-anderson-delallave-*`) — continues a subset of the
   PCRTBP resonant orbits to the physical Titania mass (µ3 = 3.91677e-5) as either quasi-periodic
   tori (irrational period ratio with Titania's synodic period) or secondary-resonant periodic orbits
   (rational ratio), studying how Titania's forcing perturbs/destroys the Oberon resonant-orbit
   structure (§"Oberon Resonant Orbits in the Uranus-Oberon-Titania CCR4BP," pp.12-18).

**Key PCRTBP results (sourced, §"Unstable Oberon Resonant Periodic Orbits," p.7-8):**
- Computed Jacobi-constant ranges for each resonant family: 3:4 ∈ [2.9916, 3.0261], 4:5 ∈
  [2.9914, 3.0157], 5:6 ∈ [2.9921, 3.0104], 4:3 ∈ [2.9836, 3.0279], 5:4 ∈ [2.9902, 3.0165], 6:5 ∈
  [2.9949, 3.0109].
- L1/L2 libration-point Jacobi constants at this µ = 3.54326e-5: C(L1) = 3.00454, C(L2) = 3.00450.

**Key heteroclinic-transfer results (§"Resonant Orbit Stable/Unstable Manifolds...," pp.8-12,
PCRTBP-only, Titania NOT included):**
- Exterior: no heteroclinics for C > 3.01 even at 15-rev manifold globalization. First 4:5↔5:6
  heteroclinic at C = 3.0072 (10-rev globalization). First direct 3:4↔5:6 heteroclinic at C = 3.0059
  (10 revs); drops to a 6-rev / 12-rev-TOF connection by C = 3.0052; a 4-rev / 8-rev-TOF connection
  requires C = 3.0039.
- Interior: no heteroclinics for C > 3.0080. First 6:5↔5:4 heteroclinic at C = 3.0065 (persists down
  to 7-rev/14-rev-TOF globalization). First 5:4↔4:3 heteroclinic at C = 3.0028 (7 revs/14-rev TOF);
  the same C also yields a DIRECT 6:5↔4:3 heteroclinic (needs ≥10-rev globalization / 20-rev TOF).
  This 6:5↔4:3 direct connection shortens to 6-rev/12-rev TOF by C = 3.0011, and 5-rev/10-rev TOF by
  C = 3.0010.
- Rough TOF-per-rev bounds: Oberon period Tob = 13.46 days; exterior transfers bounded by [Tob,
  4/3·Tob] per rev, interior by [0.75·Tob, Tob] per rev.

**Key CCR4BP (Titania-forced) results (§"Oberon Resonant Orbits in the Uranus-Oberon-Titania
CCR4BP," pp.12-18) — the section most relevant to this project's own CCR4BP work:**
- **Exterior 3:4 and 5:6 families:** mostly survive continuation to physical Titania mass as
  quasiperiodic tori with only minor perturbation. Exception: a strong 3:10 secondary resonance
  (near L = 1.036) inside the 5:6 family PREVENTS torus continuation in a moderately wide region
  around it — the only clearly-identified secondary-resonance disruption on the exterior side.
- **Interior 4:3 family:** the PCRTBP orbit family physically INTERSECTS Titania's orbit (unlike the
  6:5 family, which comes close but does not cross it) — an infinite-perturbation singularity for any
  µ3 > 0, so NO quasiperiodic torus continuation is possible at all for most of this family; only
  secondary-resonant periodic orbits (and their own manifolds) can exist there, and their computation
  is explicitly left as "a work in progress" (p.16, not completed in this paper).
- **Interior 6:5 family — the dramatic finding:** most of the family fails to continue as tori due to
  overlapping secondary resonances with Titania's synodic period. Six low-order secondary resonances
  (listed p.17, order = p+q: 25:69, 21:58, 17:47, 30:83, 13:36, 22:61, 9:25 — orders 94, 79, 64, 113,
  49, 83, 34 respectively) are shown to CONSECUTIVELY OVERLAP (Fig. 12) between the 21:58 orbit
  (Jacobi C = 3.007714) and the 9:25 orbit (C = 3.005723) — meaning NO PCRTBP 6:5 orbit with Jacobi
  constant between those two endpoints can persist as a quasiperiodic torus under Titania's real
  forcing; the family "undergoes a complete structural change" there, replaced by secondary-resonant
  orbits. Below C = 3.007714 it is "all but certain" no individual 6:5 orbit survives as a torus.
  This directly implies (explicitly stated, p.18) that the MOST USEFUL 6:5 orbits for heteroclinic
  transfers (per the PCRTBP finding that 6:5↔5:4 heteroclinics only start at C = 3.0065, already
  BELOW the 3.007714 threshold) fall exactly in this secondary-resonance-dominated regime — a genuine
  obstacle for any future 6:5-based Uranian tour design using this resonance.
- **Explicit scope limit, stated in Conclusion (p.19):** "In this study, we do not yet look at the
  heteroclinics between mean motion resonances in the CCR4BP" — the heteroclinic-connection results
  above are ALL PCRTBP-only (Titania excluded); the CCR4BP section studies only orbit-family/torus
  STRUCTURE and secondary-resonance overlap, not actual manifold-based heteroclinic or homoclinic
  connections in the 4-body model. The paper's own Conclusion explicitly lists "the same [heteroclinic
  characterization] should also be done in the CCR4BP" as unfinished future work, contingent on "first
  requir[ing] an understanding of the resonant orbit families' structure in the 4-body model" — i.e.
  this paper is the prerequisite structural study, not the connection-finding study itself.
- **Explicit future-work list (Conclusion, p.19):** (1) finish characterizing Titania's effect on
  interior Oberon resonances (including the still-incomplete 4:3 secondary-resonant-orbit
  computation); (2) **"there are three other large moons of Uranus - Titania, Umbriel, and Ariel - for
  which this study should be repeated"**; (3) compute actual CCR4BP heteroclinics (not just PCRTBP
  ones); (4) study interaction with libration-point orbits.

## 2. Cross-check — Oberon MMR survey vs. this project's `#312`/`#569` and `#701`-`#708` Uranian results

This is the single most important cross-check in this task; read carefully, grounded directly in
`data/catalogue.yaml` rows and the `#701`-`#708` docs/notes, not inference from the abstract.

### 2a. vs. `#312`/`#569` (Umbriel-Oberon-Umbriel patched-conic quasi-cycler, `umbriel-oberon-1-1-uranian-quasi-cycler-2026`)

**No overlap — different object class entirely, nothing to cross-check numerically.** `#312`/`#569`'s
catalogue row is a **multi-arc, patched-conic gravity-assist tour**: two Kepler/Lambert legs
(Umbriel→Oberon, Oberon→Umbriel, each T_leg = 14.940560615336594 days per the row's
`transit_times_days`), matched via CR3BP V∞-continuity at each moon flyby, discovered by this
project's own repeated-moon-encounter symmetric-closure genome search (`#254`/`#285`/`#312`/`#563`
lineage) — explicitly documented in the row's own notes as "NOT a single CR3BP rotating-frame
periodic orbit... a multi-arc tour across two different Uranus-moon pairs." This Anderson & Kumar
paper, by contrast, studies **single-moon (Oberon-only) resonant periodic orbits and their invariant
manifolds** within one fixed Uranus-Oberon(-Titania) system — a genuinely different mathematical
object (a continuous PO family + its manifolds, vs. a discrete two-flyby Lambert-arc tour). Neither
paper's specific numbers (Jacobi constants, TOF, resonance ratios) map onto the other's; there is
nothing here to numerically compare against `#312`/`#569`. The methods are also unrelated (patched
Keplerian-arc symmetric closure vs. CR3BP continuation + manifold globalization).

### 2b. vs. `#701`-`#708` (Umbriel 1:2-exterior resonant torus + homoclinic connection, Titania-forced CCR4BP; `umbriel-1-2-torus-homoclinic-uranus-2026`)

**This IS the directly comparable result — same model class, same perturbing moon, different target
moon and different specific resonance. No contradiction; genuinely complementary, and in one narrow
respect the project's own result is methodologically AHEAD of what this published paper achieves.**

Grounding facts pulled directly from the catalogue row and its `ccr4bp_provenance` block:
- `model_assumption: ccr4bp`, `base_resonance: "spacecraft:Umbriel=1:2-exterior"` — i.e. an unstable
  spacecraft periodic orbit resonant 1:2 with Umbriel (exterior, larger semimajor axis), in a
  Uranus-Umbriel base restricted 3-body system, continued to a quasi-periodic torus under Titania's
  physical forcing (`mu_gan`/`a_gan`/`omega_gan` fields — literally the Kumar-lineage CCR4BP
  Ganymede-perturber naming convention reused for Titania, confirming this project's own
  `core/ccr4bp_umbriel_titania.py` module is built on the SAME Blazevski-Ocampo/Kumar-Anderson CCR4BP
  formulation this Oberon paper (and its cited Jovian sibling papers already in corpus) uses).
- The row's finding is a **torus + its own homoclinic manifold connection** (unstable manifold of the
  torus reconnecting to its stable manifold, `orbit_class: torus_homoclinic`), verified with an
  independent Radau-vs-DOP853 integrator cross-check (`integrator_delta_km` ~1e-7 km) and a
  real-ephemeris (URA111/DE440) recurrence check across 10 epochs spanning 2000-2083, all 10
  "comparable" (`n_epochs_comparable: 10`, `fraction_epochs_comparable: 1.0`).

**Point-by-point comparison:**

1. **Target moon differs:** this paper studies Oberon (outermost large Uranian moon) resonant with
   Uranus, with Titania as the perturber. `#701`'s object studies UMBRIEL (a different, more interior
   moon) resonant with Uranus, ALSO with Titania as the perturber. This paper never computes, sweeps,
   or mentions any Umbriel resonant orbit anywhere in its text — confirmed by a full-text read (the
   only Umbriel mention anywhere in the paper is the one-line future-work item quoted above, "this
   study should be repeated" for Umbriel and Ariel). So there is no direct numerical overlap to check
   agreement/disagreement against.
2. **Resonance ratio differs:** this paper's Oberon study covers 3:4, 4:5, 5:6 (exterior) and 4:3,
   5:4, 6:5 (interior) — all near-1:1 (low-order, m/n close to 1) resonances. `#701`'s Umbriel torus is
   a 1:2-exterior resonance — a substantially higher-order commensurability than anything studied in
   this paper. No direct numeric match is expected or found.
3. **Object type differs in a way that matters — this project's result is methodologically ahead of
   what this paper delivers for ITS OWN target moon:** this paper explicitly states (Conclusion, p.19,
   quoted in full above) that it does **not** compute actual CCR4BP heteroclinic/homoclinic
   connections for Oberon — only PCRTBP-only heteroclinics (no Titania) and CCR4BP torus/orbit-family
   STRUCTURE (no connections). Computing genuine CCR4BP-native connections is explicitly listed as
   unfinished future work, contingent on first understanding the 4-body orbit-family structure (which
   is what this paper itself delivers, for Oberon). `#701`'s result — a converged CCR4BP torus WITH
   its own real, independently-verified homoclinic manifold connection, for Umbriel — is exactly the
   next step this paper explicitly says has not yet been done for Oberon. This does not mean `#701`'s
   result is "more novel than a published paper" in a strong sense (different moon, different
   resonance, and this project's own novelty gate for `#701` was already run independently per
   `[[project_novel_findings_status]]` and the `#701`/`#702` DONE bullets in `data/OUTSTANDING.md`) —
   but it is a genuine, citable observation: this project has independently reached, for a sibling
   Uranian moon, a class of CCR4BP result (an actual manifold connection, not just torus structure)
   that this specific 2024 AAS paper's own authors flag as not yet achieved for their own target moon.
4. **No contradiction is possible or found:** since neither the specific moon, resonance, nor
   connection type overlaps, there is nothing in this paper that could contradict `#701`'s finding.
   The two results are complementary members of the same emerging Uranian CCR4BP research program
   (independently, using the same authors' own model formulation).
5. **A genuinely new methodological observation worth flagging for a future task (inferential, not
   proven — stated honestly as such):** this paper's own secondary-resonance theory (used repeatedly
   for Oberon 6:5, §"6:5 Oberon resonance," p.16-18) is that CCR4BP tori become unreliable or fail to
   exist near LOW-ORDER (small p+q) rational values of the stroboscopic-map rotation number relative
   to the perturber's synodic period — order 34-94 secondary resonances were shown to consecutively
   overlap and destroy Oberon 6:5 tori across a real, measurable C-range. The `#701` Umbriel torus's
   own `torus_rho_strob` field is **5.995568015306847** — per this project's own
   `search/variational_ccr4bp_torus.py` docstring, `rho_strob = omega2 * period` is explicitly defined
   as "the advance of theta2 per [perturber] synodic period (the stroboscopic-map rotation number)" —
   i.e. the SAME quantity Kumar's own secondary-resonance framework is built on. 5.995568 differs from
   the exact integer 6 (a 6:1 secondary resonance, order p+q = 7 — a MUCH LOWER, and per Kumar's own
   "resonance strength decreases exponentially with order" principle, much STRONGER-effect order than
   any of the 34-94-order resonances that were shown to actually destroy Oberon 6:5 tori) by only
   ~0.074%. This project's own `torus_closure_residual` (1.427e-4) and the independent integrator
   cross-check (~1e-7 km) both indicate a genuinely converged, well-resolved torus — so this is NOT
   raising doubt that `#701`'s object is real. But it IS a well-grounded, citable observation that the
   `#701` torus sits unusually close to what this paper's own framework would flag as a
   low-order/strong secondary resonance, and this paper's method (computing the secondary-resonance
   PERIODIC ORBIT and its separatrix width at exactly rho=6, per its own §"Secondary Resonant Periodic
   Orbits" methodology) could directly test whether the `#701` torus is safely outside that
   resonance's libration zone or (less likely, given the clean closure residual) uncomfortably close to
   its edge. Flagging this as a genuinely useful, citation-grounded follow-up idea for a future task —
   NOT executed here (out of this task's scope), and NOT a claim that `#701` is wrong.

**Direct answer for the dispatching session:** this paper does not overlap with, extend, or
contradict either `#312`/`#569` (different object class, no numeric comparison possible) or
`#701`-`#708` (different moon and resonance ratio within the same Uranus-Titania-forced CCR4BP model
class — complementary, not overlapping; no contradiction is possible or found). The one substantive
connection is that `#701`'s achieved CCR4BP homoclinic-connection result is precisely the class of
result this paper's own authors explicitly flag as NOT YET DONE for their own target moon (Oberon),
and this paper's own secondary-resonance methodology suggests a concrete, well-motivated robustness
check worth running against `#701`'s near-integer rotation number in a future task.

## 3. Mandatory citation-mining pass

All 25 references read; cross-checked against `docs/notes/CORPUS_INDEX.md` (2026-07-27) and the
`cyclers_pdf/papers/` filename listing directly.

**Already in corpus** (no action needed — confirmed by direct grep, not assumed):
- [15] Kumar, Anderson, de la Llave 2023, AAS 23-397 — `kumar-anderson-delallave-2023-secondary-resonance-overlap-ganymede-4-3-ccr4bp-AAS-23-397-arxiv-2309.06073.pdf`.
- [16] Kumar, Anderson, de la Llave 2022, CMDA 134(1):3 — likely the same paper already flagged as a
  "possible-distinct-candidate" by the sibling `#727` digest (title matches "Rapid and accurate
  methods for computing whiskered tori..."); not re-flagged here, already tracked.
- [21] Kumar, Anderson, de la Llave 2021, CNSNS 97:105691 — `kumar-anderson-delallave-2021-highorder-resonant-manifold-expansions-cnsns-arxiv-2109.14800.pdf`.
- [22] Parker & Anderson 2014, *Low-Energy Lunar Trajectory Design* — already flagged (low priority) by
  the `#727` digest's own citation-mining pass; not re-flagged.
- [24] Haro, Canadell, Figueras, Luque, Mondelo 2016 — already flagged (medium priority) by `#727`'s
  citation-mining pass; not re-flagged.
- [18] Blazevski & Ocampo 2012 — already flagged HIGH priority by `#727`'s citation-mining pass (the
  CCR4BP model-definition paper); this paper independently confirms the same gap (this project's own
  `#701` CCR4BP code lineage traces to this same un-corpused source paper via TWO independent Kumar
  papers now) — reconfirmed HIGH priority, not double-counted as a new item.
- [14] Heaton & Longuski 2003, "Feasibility of a Galileo-Style Tour of the Uranian Satellites," JSR
  40(4):591-596 — **confirmed already in corpus** (`heaton-longuski-2003-feasibility-galileo-style-tour-uranian-satellites-jsr-doi-10.2514-2.3981.pdf`,
  status `mined-by-catalogue + KNOWN_CORPUS`, OCR'd) via direct `CORPUS_INDEX.md` grep — checked
  before flagging per policy, genuinely NOT a gap.
- [20] Chirikov 1960 — foundational resonance-overlap theory, already well-represented per the sibling
  `#727` digest's own citation note; not flagged as a new acquisition need.

**Genuinely new candidates, flagged NOT acquired** (per task instructions):
1. **Anderson, R.L. & Lo, M.W. 2009**, "Role of Invariant Manifolds in Low-Thrust Trajectory Design,"
   JGCD 32(6):1921-1930 [ref 1]. Medium priority.
2. **Anderson, R.L. & Lo, M.W. 2010**, "Dynamical Systems Analysis of Planetary Flybys and Approach:
   Planar Europa Orbiter," JGCD 33(6):1899-1912 [ref 2]. Medium-high priority — already independently
   flagged by the sibling `#727` digest's own citation-mining pass (same reference, ref [8] there);
   reconfirmed, not double-counted.
3. **Anderson, R.L. & Lo, M.W. 2011**, "A Dynamical Systems Analysis of Planetary Flybys and Approach:
   Ballistic Case," JAS 58:167-194 [ref 3]. Medium-high priority — also independently flagged by `#727`
   (ref [1] there); reconfirmed.
4. **Anderson, R.L. & Lo, M.W. 2011**, "Flyby Design using Heteroclinic and Homoclinic Connections of
   Unstable Resonant Orbits," Adv. Astronaut. Sci. 140:321-340 [ref 4]. Medium-high priority — also
   independently flagged by `#727` (ref [2] there); reconfirmed. Directly relevant title (heteroclinic
   AND homoclinic connections of unstable resonant orbits — exactly `#701`'s object class).
5. **Anderson, Campagnola, Koh, McElrath & Woollands 2021**, "Endgame Design for Europa Lander:
   Ganymede to Europa Approach," JAS 68(1):96-119 [ref 5]. **HIGH priority** — already independently
   flagged HIGH by `#727`'s citation-mining pass; reconfirmed by a second, independent paper (now
   flagged by two of the two most recent Kumar-lineage digests).
6. **Anderson, R.L., Campagnola, S. & Buffington, B.B. 2018**, "Analysis of Petal Rotation Trajectory
   Characteristics," JGCD 41(4):827-840 [ref 6]. Low-medium priority — new (not previously flagged),
   tour-design-technique paper.
7. **Anderson, R.L. & Lo, M.W. 2014**, "Spatial Approaches to Moons from Resonance Relative to
   Invariant Manifolds," Acta Astronautica 105:355-372 [ref 7]. Medium priority — new, directly
   relevant title (spatial/3D moon-approach-via-resonance, a natural extension of this project's own
   currently-planar Uranian CCR4BP work).
8. **Anderson, R.L. 2015**, "Approaching Moons from Resonance via Invariant Manifolds," JGCD
   38(6):1097-1109 [ref 8]. Medium priority — already independently flagged by `#727` (ref [3] there);
   reconfirmed by a second paper.
9. **Strange, N.J., Landau, D.F. & Longuski, J.M. 2014**, "Design of Initial Inclination Reduction
   Sequence for Uranian Gravity-Assist Tours," Adv. Astronaut. Sci. 150:1469-1485 [ref 12]. **HIGH
   priority** — new, directly Uranian-tour-design, and this paper's own intro cites it as the standard
   approach for the inclination-reduction phase preceding an Oberon-first tour (this project's own
   Uranian moon-tour scoping threads, `#552`/`#571`-`#579`, are inclination-relevant).
10. **Landau, D., Davis, A. & Karimi, R. 2023**, "Trajectory Options for a Uranus Orbiter and Probe,"
    AAS 23-460 [ref 13]. **HIGH priority** — new, a very recent (2023) Uranus mission-design study
    directly relevant to this project's own Uranian-system focus.
11. **Howell, K.C., Davis, D.C. & Haapala, A.F. 2012**, "Application of Periapse Maps for the Design
    of Trajectories Near the Smaller Primary in Multi-Body Regimes," Math. Probl. Eng. 2012:351759
    [ref 23]. Medium priority — new, a periapse-map design methodology potentially reusable for this
    project's own moon-encounter-design work.
12. **Celletti, A. 2010**, *Stability and Chaos in Celestial Mechanics*, Springer [ref 17]. Low-medium
    priority — new, a foundational textbook on the resonance/stability theory this paper's secondary-
    resonance analysis is built on.
13. **Morbidelli, A. 2002**, *Modern Celestial Mechanics: Aspects of Solar System Dynamics*, Taylor &
    Francis [ref 19]. Low-medium priority — new, another foundational resonance-theory textbook (cited
    for the "resonance order/strength decreases exponentially" principle used throughout this paper's
    §"6:5 Oberon resonance" analysis).

None of the above were acquired this pass (per task scope: flag only). Two items (Strange-Landau-
Longuski 2014, Landau-Davis-Karimi 2023) are new HIGH-priority, Uranian-specific gaps not previously
flagged by any prior digest — the strongest new leads from this pass. Blazevski & Ocampo 2012 remains
the single most-corroborated gap across the whole Kumar-lineage citation-mining effort (now flagged by
three independent digests: `#727` and both papers in this `#728` pass reference it, though only the
Oberon paper cites it directly here — the Kumar-Moreno paper does not cite it at all, being pure CR3BP
not CCR4BP).

## Summary answers (for the dispatching session)

- **Oberon-vs-`#312`/`#569`/`#701`-`#708` cross-check:** no overlap with `#312`/`#569` (different
  object class — patched-conic multi-arc tour vs. single-moon resonant-orbit family, nothing to
  compare numerically). Complementary (not overlapping, not contradicting) with `#701`-`#708`: same
  CCR4BP model class and perturbing moon (Titania), different target moon (Oberon here vs. Umbriel in
  `#701`) and different resonance ratio (near-1:1 here vs. 1:2-exterior in `#701`); `#701`'s achieved
  CCR4BP homoclinic connection is explicitly the kind of result this paper's own authors say is NOT
  YET done for their own target moon. A concrete, citation-grounded (not proven) robustness-check idea
  is flagged: `#701`'s torus rotation number (5.9956) sits close to a low-order (order-7) secondary
  resonance per this paper's own framework, worth a future targeted check, not treated as a problem
  here.
- **Periodic-orbit-network-vs-`#570` cross-check:** no overlap. `#570`'s schema is a fleet/downlink
  concept for sets of already-catalogued cyclers; Kumar & Moreno's "network" is a bifurcation graph
  linking periodic-orbit families within one CR3BP system. Different objects sharing only a name; no
  schema correction warranted; no specific orbit in the paper maps to an existing catalogue row.
