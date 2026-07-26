# Digest: Gilliam & Bettinger 2024, "Formulation of the Circular Restricted N-Body Problem (CRNBP) in the Jovian system"

**Paper:** "Formulation of the Circular Restricted N-Body Problem (CRNBP) in the Jovian
system"
**Venue:** *Celestial Mechanics and Dynamical Astronomy* 136, article 54 (2024), DOI
`10.1007/s10569-024-10216-0`. Received 2024-04-22, accepted 2024-09-18, published
2024-11-21.
**Authors:** Annika J. Gilliam & Robert A. Bettinger, Dept. of Aeronautics and
Astronautics, Air Force Institute of Technology (AFIT), Wright-Patterson AFB, OH.
**Access status:** the Springer/CMDA article itself is **paywalled** (buy-PDF price
€39.95 / institutional subscription; confirmed directly by fetching the live article page
— abstract + full reference list are free, body text is not). No open-access route was
found (no arXiv posting, no author-hosted preprint PDF at AFIT Scholar for the journal
article itself).
**What was actually acquired instead:** Gilliam's own follow-on AFIT MS thesis, **"Extension
of the Circular Restricted N-Body Problem (CRNBP) to Varying Multi-Body Gravitational
Systems"** (March 2025, AFIT Scholar `etd/8309`, advisor Robert A. Bettinger, freely
downloadable, Distribution A / public release, no embargo in effect). The thesis's own
§7.3.1 lists this CMDA paper as publication #1 arising from the thesis, and Chapter IV's
own text states explicitly: **"The research in this chapter was also published in
Celestial Mechanics and Dynamical Astronomy [58]"** (ref [58] = this exact paper).
Chapters III (EOM derivation) and IV (Jupiter-Europa CR3BP-to-CRNBP comparison, phase-angle
effects) are therefore the verbatim content of the target paper. The thesis is filed as
the corpus PDF for this task since it is freely available and content-equivalent (plus a
superset — see §5 below) where the journal version is not.
**Filed:** `gilliam-2025-crnbp-multibody-systems-thesis-afit-etd-8309.pdf` (private
`cyclers_pdf` corpus, commit `8f773c3`). Text-layer PDF (native, no OCR needed), 204 pp,
`pdftotext -layout` extraction clean.
**Acquired/digested:** 2026-07-26 (`#711`, dispatched alongside `#710`/`#712`).

**Relationship to `#712`'s Negri & Prado 2022 digest (`2026-07-26-digest-negri-prado-2022-crnbp.md`):**
confirmed independently here — the CMDA paper's reference list (fetched live from
Springer) cites "Batista Negri, R., Prado, A.F.B.A.: Circular restricted N-body problem.
J. Guid. Control, Dyn. 45(7), 1357-1364 (2022)" directly. This is exactly the pairing
`#712`'s digest anticipated: Negri & Prado is the GENERAL N-body EOM framework, this
paper is the real-Jovian-system APPLICATION of it. The thesis goes further and states its
own Chapter III derivation was done **independently via Lagrangian mechanics** (rather
than Negri & Prado's Newtonian route) "to validate the equations of motion previously
developed... independently utilizing different equations of motion to reach the same set
of dynamics" (thesis §7.2) — i.e., this is a from-scratch independent re-derivation that
confirms Negri & Prado's EOM, not a copy.

## 1. The CRNBP equations of motion (thesis Ch. III, = paper's own formulation)

Confirmed **structurally identical** to Negri & Prado's general CRNBP EOM already
transcribed in `#712`'s digest (down to matching term-for-term). Setup: primary `M1` +
secondary `M2` on the standard CR3BP circular orbit about their barycenter; every
additional body `j = 3...N-1` on its own CIRCULAR, COPLANAR orbit **about the primary**
(not the barycenter) with phase angle `psi_j(t) = psi_j0 + (n_j - n12)*t` — this project's
`ccr4bp.py` uses the exact same phase-evolution convention (`omega_gan = n3/n2 - 1`).

Full EOM (thesis Eqs. 25-27, spacecraft state `(x,y,z)` in the synodic frame):

```
x'' = 2y' + x - mu1(x+mu2)/r1^3 - mu2(x-mu1)/r2^3
      - sum_j mu_j [ (x+mu2-Rj*cos(psi_j))/rj^3
                      + sum_{k!=j} mu_k*(Rj*cos(psi_j)-Rk*cos(psi_k)) / rkj^3 ]
y'' = -2x' + y - mu1*y/r1^3 - mu2*y/r2^3
      - sum_j mu_j [ (y-Rj*sin(psi_j))/rj^3
                      + sum_{k!=j} mu_k*(Rj*sin(psi_j)-Rk*sin(psi_k)) / rkj^3 ]
z'' = -mu1*z/r1^3 - mu2*z/r2^3 - sum_j mu_j*z/rj^3
```

**The load-bearing structural point (confirms `#712`'s finding independently):** each
extra body `j`'s forcing term carries an INNER cross-sum over every OTHER extra body
`k != j`, coupling body `j` and body `k` through their mutual separation `rkj`. This is
**not** present in `ccr4bp.py` because that module only ever has ONE extra body
(Ganymede) — with a single extra body the inner sum is vacuous (only `k in {1,2}`, i.e.
the primaries themselves, reducing to the ordinary direct+indirect BCR4BP-shape term
`ccr4bp.py` already implements). The cross term only activates at `N>=5` (two or more
simultaneous extra perturbers) — exactly the regime this paper's own worked case (Io AND
Ganymede simultaneously forcing Jupiter-Europa) lives in, and exactly the regime a future
5-body extension of `ccr4bp.py` would need to add.

**Practical implication for extending `ccr4bp.py`:** this is a bounded, well-specified
generalization (per both this paper's own derivation and `#712`'s independent read of
Negri & Prado) — add one more `(mu_j, R_j, psi_j)` triple per extra body plus the pairwise
inner-sum coupling term, not a fresh model. Not a pure superposition of independent
single-moon `_moon_acceleration()`-style calls, though — the cross term must be computed
explicitly for every pair of extra bodies.

## 2. Jovian worked results (thesis Ch. IV = the paper's own content)

Four cases built from the Jupiter-Europa CR3BP with added perturbers: **CRNBP4** (+Io),
**CRNBP5** (+Io+Ganymede), **CRNBP6** (+Io+Ganymede+Callisto), and one **CRNBP10** case
(+Galilean moons +Metis/Adrastea/Amalthea/Thebe, the four innermost Jovian moons, chosen
"based on their low eccentricities and inclinations" — confirms the coordinating session's
prior search snippet exactly).

Unit conversions (Table 2, JPL Three-Body Periodic Orbit Catalog): DU = 668,519 km,
TU = 48,562 s, `mu2` (Europa) = 2.528017528540000e-5, `n12 = 1`.

Method: seven CR3BP periodic-orbit families (L5 Vertical, L3 N Halo, L3 Axial, 1:1/3:4/2:3
resonant, DRO) taken verbatim from the JPL catalog (full 6-DOF ICs given in thesis Table
3 — a reproducible SETUP, see caveat below) and re-propagated unchanged under CRNBP4/5/6
(RK4, variable step) at three different initial phase-angle configurations for the added
moons.

**Findings (qualitative, figure-based — no digit-grade divergence table given for any
case):**
- None of the tested CR3BP-periodic orbits remain periodic under CRNBP4/5/6/10 — all
  "experience significant perturbations." The most common outcome is the orbit drifting
  in the direction of the moons' own motion while roughly preserving its shape.
- Effect magnitude tracks proximity: the 3:4 resonant orbit (which stays close to
  Ganymede's path) is barely perturbed by Io alone but strongly perturbed once Ganymede is
  added; the 1:1 and 2:3 resonant orbits (closer to Io/Callisto) are heavily perturbed in
  all three configurations.
- Initial phase-angle configuration matters: one CRNBP5 case showed a trajectory
  temporarily relocating ~90° from its nominal position before drifting further —
  flagged by the authors as evidence "intermediate conditions with greater stability may
  exist in the CRNBP," i.e. phase-angle selection could be used for mission-planning
  window selection, not evidence of an actual periodic/cyclic structure.
- Explicit, important **open-question admission by the authors themselves** (thesis
  §4.5, echoed in the paper's own abstract's framing): *"it is not currently known if any
  periodic trajectories exist in the CRNBP beyond [resonance-locked] cases."* This paper
  does **not** establish existence of periodic or quasi-periodic structures in the CRNBP
  the way Kumar et al. 2021 established resonant tori for CCR4BP — it is a
  CR3BP-vs-CRNBP divergence/sensitivity study, not a periodic-orbit or connection-search
  paper.

**Reproducible positive control assessment:** Table 3's seven IC sets + Table 2's unit
conversions ARE a directly re-runnable SETUP (drop straight into a Jupiter-Europa CR3BP
propagator, matching this project's own `cr3bp.py`/`ccr4bp.py` unit convention), but the
CRNBP comparison OUTCOME itself is qualitative/visual only (drift direction, "similar
shape," no divergence-vs-time numbers) — **not a digit-grade positive control** per
`[[feedback_golden_tests_sourced_only]]`, same limitation `#712`'s digest already found
in Negri & Prado's own worked examples.

## 3. Beyond the paper's own scope: thesis Ch. V/VI generalize to other systems (bonus, not journal-paper content)

The narrower published journal paper (per its own abstract, confirmed by live fetch) is
Jovian-only. The thesis is a superset that goes further, and this additional material —
while not technically "in" the paywalled article being acquired — is highly relevant to
the future 5-body grounding question and freely available from the same authors/framework:

- **Chapter V (equilibrium/"dynamical substitute" points):** computes CRNBP-perturbed
  analogues of the CR3BP L1/L2 Lagrange points ("Lagrange Box" via multi-dimensional
  Newton-Raphson) for SIX systems: Jupiter-Europa, Sun-Earth, Sun-Ceres, **Uranus-Oberon**,
  Saturn-Titan, Saturn-Enceladus. This directly answers the task's generality question:
  **the CRNBP framework is NOT Jupiter-specific** — it is already worked for the Uranian
  system this project's own CCR4BP arc (`#689`-`#708`) targeted. Digit-grade numeric
  tables are given (Table 5 = E1 point, Table 6 = E2 point; columns: `dx (km)`, `dy (km)`,
  Lagrange Box area (km^2), average displacement vs. CR3BP):

  | System | E1 dx/dy (km) | E1 Box area (km^2) | E1 avg-vs-CR3BP |
  |---|---|---|---|
  | Jupiter-Europa | 57.86 / 44.07 | 2,549.59 | 4.279 |
  | Uranus-Oberon | 77.99 / 128.07 | 9,987.70 | 7.444 |
  | Saturn-Titan | 17.71 / 43.02 | 761.66 | 1.144 |
  | Saturn-Enceladus | 0.316 / 0.886 | 0.280 | 0.0421 |
  | Sun-Earth | 1,117.2 / 2,032.7 | 2,270,000 | 18.93 |
  | Sun-Ceres | 63,698 / 217,821 | 1.3875e10 | 2,972.6 |

  (E2 table, same systems, given in the digest file's source thesis Table 6.) For
  Uranus-Oberon the perturbers used are Titania, Ariel, Umbriel; the thesis finds Titania
  dominates (comparable in effective potential to Io/Ganymede's role at Jupiter-Europa)
  while Ariel/Umbriel are secondary (comparable to Callisto's role) — a genuinely useful,
  reproducible, DIGIT-GRADE positive-control candidate (validate a CRNBP equilibrium-point
  / dynamical-substitute Newton-Raphson solver against these numbers) even though it is
  thesis-only, not journal-paper, content.
- **Chapter VI (Poincaré mapping):** first application of Poincaré maps to the CRNBP,
  worked for Jupiter-Galilean, Saturn-Titan/Enceladus, **Uranus-Titania, Uranus-Oberon**,
  and Pluto-Hydra(+Charon+Nix+Kerberos) systems — again confirming cross-system
  generality. Findings are qualitative (map-structure comparison, not digit-grade) but the
  Pluto-Hydra case is flagged as the most dramatic divergence from CR3BP (large secondary
  perturber, Charon, located between primary and secondary) and is explicitly framed by
  the authors as the strongest argument for the CRNBP's necessity over CR3BP in
  minor-secondary systems (Jupiter-Amalthea, Saturn-Prometheus, Sun-Ceres-class cases).
- Chapter III also derives (thesis Eqs. 61-64, "for the first time" per the thesis's own
  significance claims in §7.2) the CRNBP Jacobian/partial-derivative matrix needed for
  Newton-Raphson equilibrium-point solving and differential-correction shooting — a
  directly reusable piece of machinery beyond what either the paper or `#712`'s Negri &
  Prado digest reports.

## 4. Computational cost / tractability for N>=5

**Not quantified anywhere in the accessible text** — no runtime benchmarks, step-count
comparisons, or integration-cost figures were found in either the paper's abstract/refs or
the thesis body. The thesis frames the CRNBP's value proposition qualitatively as
"computationally efficient" relative to full ephemeris models and suitable as "an
intermediate model when transitioning to a higher-fidelity model," mirroring how this
project already uses BCR4BP/CCR4BP as an ephemeris-transition rung, but gives no numbers
to substantiate the cost claim. Given the EOM's structure (one extra direct+indirect term
plus one pairwise inner-sum term per additional body, i.e., O(N^2) cost in the number of
extra perturbers for the force evaluation, versus O(N) for a naive superposition), the
per-step cost growth is bounded and cheap in absolute terms (N is small — at most ~8 extra
bodies in any case this thesis tests) — this is an inference from the EOM's own structure,
not a benchmark reported by the authors.

## 5. Overall assessment: tractable but a weaker grounding source than Kumar et al. 2021 was for CCR4BP

**Structurally tractable, yes:** the EOM generalization from `ccr4bp.py`'s N=4 case to a
genuine N=5 (or general N) case is well-specified, bounded, and independently confirmed
by two sources (this paper's own re-derivation and `#712`'s Negri & Prado digest) — not a
different or incompatible formulation. Extending `ccr4bp.py` would mean adding the
pairwise inner-sum coupling term, which is mechanical, not a research risk.

**But the discovery-strategy grounding is weaker than Kumar et al. 2021 was for the
CCR4BP arc**, for one important reason: Kumar 2021 handed `#686`'s CCR4BP arc a KNOWN,
published, PERIODIC/quasi-periodic structure (resonant invariant tori, with their
stable/unstable manifolds) to validate a torus corrector against before searching for
anything novel. Neither this paper nor its companion Negri & Prado 2022 paper (per
`#712`'s digest) does the equivalent for the CRNBP: this paper's own Jovian trajectories
are ALL non-periodic drift/divergence cases, and its authors explicitly state it is "not
currently known if any periodic trajectories exist in the CRNBP" outside special
resonance-locked cases. A future 5-body discovery-strategy pass grounded on this paper
would therefore be searching for periodic/quasi-periodic/connecting structures in a regime
where **no known-good positive control yet exists in the literature** — violating this
project's own "verify a gauntlet with a positive control before trusting 0/N" discipline
(`[[feedback_verify_gauntlet_with_positive_control]]`) unless a positive control is
constructed independently (e.g., via the epsilon-homotopy continuation method `#712`'s
digest flagged in Negri & Prado's own Jupiter-Ganymede vertical-Lyapunov-orbit example, or
by treating the reproducible dynamical-substitute tables in §3 above as a narrower,
equilibrium-only positive control while treating periodic-orbit/connection existence in
the CRNBP as the genuinely open research question).

**Net verdict:** worth keeping as the real-system application/citation pairing with
Negri & Prado 2022 (as `#712` already recommended), and the thesis's Uranus-Titania/
Uranus-Oberon material is a genuine, freely-available, cross-system generality confirmation
directly relevant to this project's own prior Uranian CCR4BP work — but any future
discovery-strategy pass built on this pairing should treat "does anything periodic survive
in the CRNBP at all" as the FIRST open question to resolve (via the homotopy-continuation
method), not assume it the way the CCR4BP arc could assume torus existence from Kumar
2021.
