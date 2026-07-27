# Digest — Casoliva et al., Earth-Moon cycler family papers (#725)

**Papers:**
1. Casoliva, J., Mondelo, J. M., Villac, B. F., Mease, K. D., Barrabés, E. & Ollé, M.,
   *"Two Classes of Cycler Trajectories in the Earth-Moon System,"* Journal of Guidance, Control,
   and Dynamics 33(5), Sept-Oct 2010, pp. 1623-1640. DOI `10.2514/1.46856`.
2. Same authors, *"Families of Cycler Trajectories in the Earth-Moon System,"* AIAA 2008-6434,
   AIAA/AAS Astrodynamics Specialist Conference, Honolulu HI, 18-21 Aug 2008 (the JGCD paper's
   conference precursor).

**Filed:** `/Users/bruce/dev/cyclers_pdf/papers/casoliva-mondelo-villac-mease-barrabes-olle-2010-two-classes-cycler-trajectories-earth-moon-jgcd-33-5-doi-10.2514-1.46856.pdf`
(+ `.txt`) and `...-2008-families-cycler-trajectories-earth-moon-AIAA-2008-6434.pdf` (+ `.txt`), commit
`3a0dfe5` in the private `cyclers_pdf` repo. **OCR status: both text-layer** (native LaTeX PDFs;
`pdftotext -layout` yields 355,842 chars/18 pp for the 2010 JGCD and 110,676 chars/20 pp for the
2008 AIAA — both far above the 10-char/page floor; no OCR needed).

**Context:** flagged as an acquisition candidate in `docs/notes/2026-06-11-ross-roberts-tsoukkas-2025-mining.md`
§7/§11 (paywalled; user directly supplied both PDFs 2026-07-27, task #713 → #725). Hypothesis under
test: this is the **unstable complement** to the "stable" Ross-Roberts-Tsoukkas Earth-Moon
`(k1,k2)`-cycler lineage already deeply mined (`docs/notes/2026-06-30-digest-ross-roberts-tsoukkas-2026-stable-ballistic-cyclers.md`).

## 1. Dynamical model

Both papers use the **planar circular restricted three-body problem (PCR3BP/PCRTBP)**, confirming
the abstract. Precisely:

- **2008 (AIAA 2008-6434):** strictly planar 4-state EOM, `X = [x, y, u, v]ᵀ` (Eq. 1). Moon-on-the-left
  convention, Earth at `(µ,0)`, Moon at `(µ−1,0)`, **µ = 0.01215** (rounded). Jacobi constant
  `C_J = 2Ω − ‖V‖²` (Eq. 3), `C = −2H` (Eq. 5). Stability index `k = λ+1/λ` (Eq. 8), stable iff `|k|≤2`.
  Continuation system: fix `(H(x)−h, g(x), φ_T(x)−x)=0` over unknowns `(µ,h,T,x0)`, one-parameter
  continuation on µ, h, or T (§II.C).
- **2010 (JGCD):** upgrades to the **full 6-state spatial EOM** `X=[x,y,z,u,v,w]ᵀ` (Eq. 1) *specifically
  so out-of-plane stability can be assessed*, even though every computed orbit in the paper is
  planar (`z=w=0`). Same Jacobi constant form (Eq. 2), `C_J = −2h` (Eq. 3; OCR of the actual paper
  drops the minus sign — confirmed by sign-consistency: their homoclinic-cycler energies `h≈−1.45`
  to `−1.59` map to `C_J≈2.89–3.18`, matching the RRT-lineage `C≈3.0–3.7` band, not a negative
  range). Stability splits into **in-plane `k∥=κ+1/κ`** and **out-of-plane `k⊥=κ⊥+1/κ⊥`** (Eqs. 6-8),
  monodromy eigenvalue spectrum `{κ,1/κ,κ⊥,1/κ⊥,1,1}` — i.e. the same planar/vertical stability-index
  split later used by Ross & Roberts-Tsoukkas 2026 (their `sₚ,sᵥ`), just defined via `λ+1/λ` (Barden
  half-period-adjacent convention) rather than `½(λ+1/λ)`.
- µ_EM used for actual Earth-Moon computations: **0.0121529529** (2010, precise) vs **0.01215**
  (2008, rounded) — same physical system, the 2010 paper just carries more digits.
- Both use the moon-on-the-left sign convention consistent with the rest of this project's corpus
  (Andreu, Ross-RT, etc.) — no coordinate-convention translation needed for any future
  reproduction attempt.

## 2. Classification scheme, methods, and concrete numbers

### Class 1 — High-energy near-Keplerian resonant cyclers (§IV both papers)

Organized around **p-q resonant orbits** (`q·T_moon = p·T_spacecraft`, Eq. 9): the spacecraft
completes `q` inertial ellipses while the Moon completes `p` revolutions. Two seed/continuation
strategies:

1. **Elliptical-orbit seed + differential correction at fixed T** (µ jumped directly to µ_EM, no
   continuation in µ). Works for `p>q` (larger `a_s/c`) but fails below `a_s/c≈0.6` (Earth-Moon
   distance units) — cannot produce *tight* cyclers (small perigee AND periselene simultaneously).
2. **Small-µ second-species-solution seed** (Hénon's terminology) using the **Barrabés & Gómez
   analytic in/out-map matched-asymptotics** (their Eqs. 14-18, reproduced verbatim in both
   Casoliva papers) at µ=10⁻⁶, then **continued in µ up to µ_EM at fixed C_J** (a 3-step strategy:
   µ-continuation → switch to C_J-continuation when lunar impact threatens → resume µ-continuation).
   This is the ONLY method that reaches tight cyclers, demonstrated for the 7-3 resonance
   (`a_s/c=0.5684`).

**Golden numeric table — 2010 JGCD Table 3** (actual µ_EM=0.0121529529, all planar, Poincaré section
`y=0`; recovered from the OCR'd/garbled table via cross-referencing Fig. 4's captions, which are
clean):

| Designation | p-q | C_J | Stability |
|---|---|---|---|
| 1-2c | 1-2 | 1.5691874798 | **stable** |
| 1-2d | 1-2 | 2.5803060666 | unstable |
| 1-2e | 1-2 | 2.7629814961 | **stable** |
| 2-1a | 2-1 | 0.4887353098 | **stable** |
| 2-1b | 2-1 | 1.1964188553 | unstable |
| 3-2c | 3-2 | 0.7089330385 | **stable** |
| 7-3a | 7-3 | 1.0215696153 | **stable** |
| 7-3b | 7-3 | 1.0687623900 | unstable |
| 7-3c | 7-3 | 1.0687623900 | unstable (7-3b flipped 180°) |

Insertion ΔV (Eq. 19, LEO-perigee to cycler) for the feasible set: 1159.2 (1-2c), 376.1 (1-2d),
280.6 (1-2e), 713.1 (2-1a), 3098.7 (2-1b), 2039.8 (3-2c), 1470.3 (7-3a), 1967.5 (7-3b), 1967.6 (7-3c)
m/s. **1-2c/1-2d/1-2e (p<q) have the lowest insertion costs.** All cyclers except 2-1a have
back-of-the-Moon coverage. All are symmetric about the x-axis except 1-2e, 7-3b, 7-3c.

**Key finding for a future golden/positive-control test:** this class contains **both stable and
unstable members per resonance class** ("at least one unstable and one stable cycler per resonance
class" — direct quote) — found via elliptical/second-species differential correction, a completely
different method lineage from the manifold-tube-intersection method used for Class 2 and for
Ross-Roberts-Tsoukkas. The 2008 conference paper's Table 2 gives an earlier/smaller set of
small-µ=10⁻⁶ seed cyclers (12a, 21a, 23a/b, 32a/b, 52a, 54a/b, 73a) with full IC/C_J/T/k values —
useful as an independent seed-reproduction check of the Barrabés-Gómez in/out-map matching (Eqs
14-18) before attempting the µ-continuation to µ_EM.

### Class 2 — Low-energy L1-based homoclinic-type cyclers (§V both papers)

Organized around **homoclinic connections of the L1 Lyapunov orbit family** (the Lyapunov orbit
itself is **provably unstable** — the paper states explicitly "the Lyapunov p.o. needs
station-keeping, since it is unstable"). Method:

- Compute the L1 Lyapunov orbit's stable/unstable manifold tubes; intersect their branches on a
  Poincaré section (`Σ={x=0}` for the right branches / `{x=−1}` in 2010 for the left branches) to
  find discrete homoclinic points **He₁..He₄** (right branch) and **Hm₁,Hm₂** (left branch, Earth-Moon
  system specifics: `He1,He3` self-symmetric about `y=0`; `He2,He4` mirror-symmetric to each other).
- **Continue each homoclinic connection family in energy** via a predictor-corrector on a
  Hamiltonian/monodromy-eigenvector system (2010 Eq. 20; 2008's equivalent unnumbered system) — this
  IS the Barrabés-Mondelo-Ollé (2009, Nonlinearity) continuation-of-homoclinic-connections method,
  cited as ref [22]/[26] and reproduced as the core algorithm of both Casoliva papers.
- A true periodic orbit that "shadows" a chain of these connections (one turn around He1, or He1
  then Hm1, etc.) can then be found by direct continuation as an ordinary p.o.; Casoliva explicitly
  cross-validates some of these against **Leiva & Briozzo's atlas families 357 and 037** (Fig. 7/Fig. 16).

**Golden numeric table — He1 family closest-approach connection, energy h=−1.450162** (2010 Tables
4-6, essentially identical values in the 2008 Tables 2-3 — the two papers report the SAME connection
independently, useful as an internal cross-check):

- Lyapunov p.o. period: **29.1640 days**; connection flight time (periselene 1→19): **113.6319 days**.
- Periselene 1 (Lyapunov orbit, Moon-relative 2-body): r=6331–6579 km (2010 vs 2008 values differ
  slightly — 6331.184 km / 2010 Table 5 vs 6573.556 km / 2010 Table 6 connection-part periselene 1;
  2008's single value is 6325-6573 km range), e≈1.126–1.132 (hyperbolic).
- Perigee 8/12 (Earth-relative 2-body, closest Earth approach): a≈201,832–211,932 km, e≈0.664–0.680.
- LEO-rendezvous ΔV at that perigee: **~703–718 m/s** (2008 gives 703 m/s at r=67,808 km; 2010 gives
  717.5 m/s at r=67,869 km — same connection, consistent cross-check).
- He1 family impacts the Moon's surface at energy h=−1.4711; Hm1 impacts at h=−1.4711, Hm2 at
  h=−1.5892 (2010 §V.C).
- Eccentricity (Moon 2-body, periselene 1) ranges 0.4940–0.6797 across the family; velocity crosses
  zero (retrograde↔direct inertial-frame transition) at h≈−1.473.

This is a genuinely rich, sourced numeric family — every value above carries a table/page citation
and could seed a future golden/positive-control test for this project's own homoclinic-continuation
tooling (the repo doesn't yet have a dedicated Barrabés-Mondelo-Ollé-style homoclinic-connection
continuator; this is the first sourced numeric target for one if built).

## 3. Comparison against this project's own Earth-Moon catalogue rows / Ross-Roberts-Tsoukkas

`grep -n -i "earth.*moon" data/catalogue.yaml` for cycler-class rows turns up (among others):
`ross-rt-em-cycler-{11,21,31,32,33}-2025` (C≈3.03–3.70 per the RRT digest golden table),
`braik-ross-common-energy-{11a,11b,32}` (C=3.1294), the KAM-torus-corridor 3D-lifted families, plus
older `wittal-2022`, Genova/Aldrin, and Arenstorf rows. **None of these rows cite or match Casoliva
2008/2010** — confirmed by `grep -i casoliva data/catalogue.yaml` (zero hits) and
`literature_check.py` (zero KNOWN_CORPUS hits). This paper is genuinely new to the corpus, not a
rediscovery of an already-catalogued source.

**Verdict on the "unstable complement" hypothesis — CONFIRMED, with an important nuance:**

- **Class 2 (low-energy L1-homoclinic) is the genuine unstable complement.** It is built
  structurally around the **provably unstable** L1 Lyapunov orbit family and its homoclinic
  connections — Casoliva never checks or claims stability for any periodic orbit that shadows
  these connections; the base object requires active station-keeping by construction. This is
  categorically different from Ross & Roberts-Tsoukkas 2026, whose central result is that the
  *same* manifold-tube-intersection construction (symmetric IC from `Γ = P(S^u1_k1) ∩ S^u2_k2`,
  differential correction + pseudo-arclength continuation) produces cyclers that are **stable**
  (planar AND vertical, `|sₚ|,|sᵥ|<1`) once continued far enough from the saddle-center birth point
  — a result Casoliva's 2008/2010 analysis never surfaces, because Casoliva never continues *away*
  from the homoclinic connections themselves into the enclosed p.o. family and checks its
  stability index. The energy ranges are adjacent/overlapping (Casoliva's He1 family:
  C_J≈2.89–3.18; RRT's catalogued rows: C≈3.03–3.70) — genuinely the same dynamical neighborhood,
  approached from opposite stability directions. This matches this project's own 2026-06-11 mining
  note's framing exactly ("this paper's stable cyclers sit inside the L1/L2 tube-enclosed region
  but beyond the theorem's unstable-orbit neighborhood — a genuinely new regime, not a
  rediscovery") — Casoliva is the missing prior-art half of that same neighborhood, not a
  rediscovery of RRT, and RRT's 2026 stable result is NOT anticipated anywhere in Casoliva
  2008/2010 (neither paper computes or mentions a stability index for any homoclinic-shadowing
  periodic orbit).
- **Class 1 (high-energy p-q resonant) is NOT purely unstable** — it explicitly contains both
  stable and unstable members per resonance class (Table 3), found via an entirely different
  method (elliptical/second-species differential correction, not manifold intersection). This class
  has no direct RRT counterpart at all — RRT's (k1,k2) construction is manifold-tube-based only.
  Class 1's stable members (1-2c, 1-2e, 2-1a, 3-2c, 7-3a) are a separate, uncatalogued stable family
  this project has never mined; they could be a fresh acquisition target in their own right
  (different construction from both RRT and Braik-Ross), independent of the unstable-complement
  question.

## 4. MANDATORY citation-mining pass (corpus-document-policy step 4)

Read both papers' Introduction/Background sections and reference lists in full (2010: refs [1]-[30];
2008: refs 1-27, near-identical set). Cross-checked every citation against
`grep -i <author> /Users/bruce/dev/cyclers_pdf/papers/` (filename sweep) and
`grep -i <author> docs/notes/CORPUS_INDEX.md` (zero hits for any of the flagged names below,
confirming none are digested/mined under an alternate filename).

**Already in corpus (no action needed):**
- Koon, Lo, Marsden, Ross 2000, "Heteroclinic Connections and Resonance Transitions" (Chaos 10-2)
  — ref [25]/2010, [13]/2008 → `koon-lo-marsden-ross-2000-heteroclinic-connections-resonance-transitions-chaos-10-2.pdf` (digested).
- Gómez, Koon, Lo, Marsden, Masdemont, Ross 2004, "Connecting Orbits and Invariant Manifolds..."
  (Nonlinearity 17) — ref [16]/[30] 2010, [15] 2008 → `gomez-koon-lo-marsden-masdemont-ross-2004-...pdf` (in corpus).
- Doedel, Paffenroth, Keller et al. 2003, "Computation of Periodic Solutions..." (IJBC 13-6) —
  ref [28]/2010, [23]/2008 → `doedel-paffenroth-keller-2003-...pdf` (digested).
- Russell & Strange 2009, "Planetary Moon Cycler Trajectories" (JGCD 32-1) — ref [7]/2010, [6]/2008
  → `russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-doi-10.2514-1.36610.pdf`
  (in corpus).
- Koon, Lo, Marsden, Ross 1999, GENESIS trajectory (AAS-99-451) — related to but NOT identical to
  ref [13]/2010=[12]/2008 (Howell, Barden, Wilson, Lo 1998 — see below, still a gap).
- Belbruno 2004 capture-dynamics textbook — in corpus, but distinct from the specific 2005/2006
  Belbruno conference papers cited (still gaps, low priority — see below).
- Parker 2007 PhD thesis (low-energy ballistic lunar transfers) — related to but not identical to
  ref [18]/2010=[17]/2008 (Parker & Lo 2006, "Shoot the Moon 3D" AAS conference paper — still a gap,
  the thesis likely supersedes/contains it, low priority to acquire separately).

**Genuinely new acquisition candidates flagged (not acquired, per instructions):**

*High priority — direct method/data sources Casoliva builds on:*
- **Barrabés & Gómez 2002/2003** ("Spatial p-q Resonant Orbits of the RTBP," CMDA 84-4;
  "Three-Dimensional p-q Resonant Orbits Close to Second Species Solutions," CMDA 85-2) — refs
  [9]/[10] 2010, [8]/[9] 2008. **This is the exact analytic in/out-map seed-generation method
  (Eqs. 14-18) Casoliva reproduces verbatim and uses for Class 1's tight-cycler continuation.**
  Highest-priority gap — needed to independently verify/extend the p-q resonant seed generation.
- **Barrabés, Mondelo & Ollé 2009**, "Numerical Continuation of Families of Homoclinic Connections
  of Periodic Orbits in the RTBP" (Nonlinearity 22-12) — ref [22]/2010 (listed "in preparation,
  2008" in the AIAA paper's ref [26], now published). **This is the core homoclinic-continuation
  algorithm** (2010 Eq. 20) underlying all of Class 2. Highest-priority gap for reproducing the
  He1/Hm1/Hm2 continuation.
- **Leiva & Briozzo** — "The Earth-Moon CR3BP: A Full Atlas of Low-Energy Fast Periodic Transfer
  Orbits" (2006 preprint/ref [18] 2008) and "Extension of Fast Periodic Transfer Orbits from the
  Earth-Moon RTBP to the Sun-Earth-Moon Quasi-Bicircular Problem" (CMDA 101, 2008, pp. 225-245) —
  ref [19]/2010. **Already independently flagged** in this project's own
  `2026-06-11-ross-roberts-tsoukkas-2025-mining.md` §11 as "the (3,2)-resembling unstable member +
  quasi-bicircular persistence — the bridge to higher-fidelity validation of the (3,2) row." Now
  doubly-motivated: Casoliva's Fig. 7/Fig. 16 cross-validate directly against Leiva-Briozzo families
  357 and 037.
- **Broucke 1968**, "Periodic Orbits in the Restricted Three-Body Problem with Earth-Moon Masses"
  (JPL TR 32-1168) — ref [4]/2010, [1]/2008. Classical Earth-Moon periodic-orbit census; also
  independently flagged in the RRT mining note (ref [12], "the classical Earth-Moon periodic-orbit
  census these families thread through"). Two independent papers now point at this same gap.
- **McGehee 1969**, "Homoclinic Orbits in the Restricted Three-Body Problem" (PhD thesis, Univ. of
  Wisconsin-Madison) — ref [11]/2010, [10]/2008. Foundational homoclinic-orbit existence theory
  underlying both papers' Class 2.
- **Llibre, Martínez & Simó 1985**, "Transversality of the Invariant Manifolds Associated to the
  Lyapunov Family of Periodic Orbits Near L2 in the Restricted Three-Body Problem" (J. Differential
  Equations 58-1) — ref [12]/[23] 2010, [11]/2008. Companion theoretical-existence result to
  McGehee, specifically for the L2 Lyapunov family homoclinic/transversality question.
- **Hénon 1997**, "Generating Families in the Restricted Three-Body Problem" (Springer) — ref
  [21]/2010, [24]/2008. The classification-of-second-species-orbits reference both papers cite for
  their Class 1 seed generation.

*Medium priority — adjacent methodology/theory:*
- **Simó 1990**, "On the Analytical and Numerical Approximation of Invariant Manifolds," in *Modern
  Methods in Celestial Mechanics* (Benest & Froeschlé, eds.) — ref [27]/2010, [22]/2008.
- **Allgower & Georg 1990**, *Introduction to Numerical Continuation Methods* (Springer) — ref
  [29]/2010, [21]/2008. Foundational pseudo-arclength continuation textbook both papers cite for
  their continuation algorithm.
- **Lo & Parker 2004**, "Unstable Resonant Orbits Near Earth and Their Applications in Planetary
  Missions" (AIAA 2004-5304) — ref [17]/2010, [16]/2008. **This is the specific "Lo and Parker"
  work** the task flagged from the 2008 paper's intro ("the classification of planar, simple
  periodic symmetric families of orbits in mission designs in the Earth-Moon system was studied by
  Lo and Parker"). Casoliva explicitly extends its scope to higher-energy and asymmetric orbits.
- **Parker & Lo 2006**, "Shoot the Moon 3D" (Adv. Astronaut. Sci. 123) — ref [18]/2010, [17]/2008;
  likely superseded/contained by the already-filed Parker 2007 PhD thesis — low-priority to acquire
  separately, flagged for completeness.
- **Dunham & Farquhar 2003**, "Libration Point Missions, 1978-2002," in *Libration Point Orbits and
  Applications* (World Scientific) — ref [15]/2010, [14]/2008. Historical libration-point mission
  survey chapter.
- **Howell, Barden, Wilson & Lo 1998**, "Trajectory Design Using a Dynamical Systems Approach with
  Application to GENESIS" (Adv. Astronaut. Sci. 97) — ref [13]/2010, [12]/2008. Distinct from the
  already-filed Koon-Lo-Marsden-Ross 1999 GENESIS AAS-99-451 paper (same mission, different authors
  and year — likely the earlier design-approach paper vs. the later heteroclinic-trajectory paper).
- **Canalias & Masdemont 2006**, "Homoclinic and Heteroclinic Transfer Trajectories Between Planar
  Lyapunov Orbits in the Sun-Earth and Earth-Moon Systems" (Discrete Contin. Dyn. Syst. 14-2) —
  ref [24]/2010, [25]/2008. Distinct from the already-filed Canalias 2007 PhD thesis (may be
  contained within it — worth checking before acquiring separately).

*Lower priority — tangential/general:*
- Lara, Russell & Villac 2007 ×2 ("Fast Estimation of Stable Regions in Real Models," Meccanica
  42-5; "Classification of the Distant Stability Regions at Europa," JGCD 30-2) — refs [19]/[20]
  2008 (not carried into the 2010 reference list). Europa-focused, same-author stability-region
  methodology; tangential to Earth-Moon cyclers.
- Wiesel 1997, *Spaceflight Dynamics* textbook — ref [26]/2010. General astrodynamics textbook.
- Cooke, Joosten, Lo, Ford & Hansen 2003, "Innovations in Mission Architectures for Exploration
  Beyond Low Earth Orbit" (Acta Astronautica 53-4) — ref [1]/2010, [2]/2008. General infrastructure
  motivation, not cycler-specific.
- Johnson & Belbruno 2005, "Reduction of Lunar Landing Fuel Requirements by Utilizing Lunar Ballistic
  Capture" (NYAS 1065) — ref [5]/2010, [4]/2008. Ballistic-capture, tangential to cyclers; the
  Belbruno 2004 textbook already in corpus covers the same WSB theory.
- Belbruno 2006 conference/magazine pieces ("A Low Energy Lunar Transportation System...," AAS;
  "Low-Energy Pathways in Space," American Scientist 94-5) — refs [2]/2010, [3]/[5]/2008. Same note.
- Grebow, Ozimek, Howell & Folta 2008, "Multi-Body Orbit Architectures for Lunar South Pole
  Coverage" (JSR 45-2) — ref [3]/2010. Lunar south-pole coverage, not cycler-specific.
- Villac & Scheeres 2004, "On the Concept of Periapsis in Hill's Problem" (CMDA 90) — ref [27]/2008
  only. Explains a Hill's-problem periapsis-definition subtlety noted in a footnote; general theory.

## 5. Bottom line

Both papers are genuinely new to this project's corpus (zero catalogue/KNOWN_CORPUS overlap),
native text-layer PDFs needing no OCR, now filed+committed to the private corpus
(`cyclers_pdf` commit `3a0dfe5`) and registered in `CORPUS_INDEX.md`. The "unstable complement to
Ross-Roberts-Tsoukkas" hypothesis is **confirmed for Class 2** (the L1-homoclinic low-energy
cyclers, built on the provably-unstable Lyapunov family, energy-adjacent to but never overlapping
in construction method with RRT's stable result) and **only partially applicable to Class 1** (the
high-energy p-q resonant cyclers, which mix stable and unstable members via an unrelated
elliptical-orbit/second-species method). Thirteen citation-mining leads flagged above (7
high/medium-priority, 6 lower-priority) — none acquired, all cross-checked as absent from both
`papers/` filenames and `CORPUS_INDEX.md`.
