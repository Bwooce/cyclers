# Andreu (1998) — The Quasi-Bicircular Problem (QBCP / BCR4BP): digest + verdict

Task #258. READ/DIGEST. Autonomous background run.

## Primary source

Andreu, M.A. (1998), *The Quasi-Bicircular Problem*, PhD thesis, Univ. de Barcelona
(advisor C. Simó). The QBCP/BCR4BP: a coherent, periodically Sun-perturbed
restricted four-body model for the Sun–Earth–Moon system, one fidelity rung above
the Earth–Moon CR3BP.

### What was readable in-session

- The filed thesis (gzipped PostScript) is **not in-session-readable**: Bash access
  to the reference-PDF path is denied, and no PostScript renderer is available
  (gs / ps2ascii / pstotext all absent; pdftotext cannot read PS). Raw-PS grep was
  not executable here.
- **Digested instead from readable companions that restate Andreu's model and cite
  the thesis as the primary source** (sourced-golden discipline preserved — every
  number below traces to a published paper, not to our own code):
  - Gimeno, Jorba, et al. (2018), "Two Periodic Models for the Earth–Moon System,"
    *Frontiers in Applied Mathematics and Statistics* 4:32 — gives the QBCP/BCP
    Hamiltonian structure and the parameter table, Fourier coefficients α_i "sourced
    from Andreu (1998)."
  - Rosales, Jorba, et al. (2023), "Invariant manifolds near L1 and L2 in the
    quasi-bicircular problem," *Celestial Mechanics & Dynamical Astronomy* 135:15
    (2023CeMDA.135...15R) — gives the explicit QBCP Hamiltonian (their Eq. 3), the
    parameter table (their Table 3), and dynamical-substitute periodic-orbit ICs
    (their Table 4).
  - Springer fulltext is paywalled; values above are from the academia.edu and
    Frontiers open versions.

## The BCR4BP / QBCP model

### How the Sun enters vs CR3BP

CR3BP (our V1) is **autonomous**: two primaries on fixed circular orbits, one
Jacobi integral. The QBCP keeps the Earth–Moon rotating-pulsating frame but makes
the Hamiltonian **explicitly, periodically time-dependent** through the Sun. Unlike
the (incoherent) classical Bicircular Problem — where the Sun is bolted on as an
external term that violates Newton's laws because the primaries' prescribed motion
is not a true solution — the QBCP is **coherent**: the Sun–Earth–Moon primaries
follow a genuine periodic solution of the three-body problem that is *close to*
bicircular (hence "quasi"). Coherence is the whole point: near the Earth–Moon
collinear points the indirect Sun term matters, and the BCP's incoherence corrupts
exactly that region.

The QBCP Hamiltonian (Rosales/Jorba 2023, Eq. 3), in Earth–Moon RTBP coordinates
and units, with all time dependence carried by the angle ϑ = ω_s·t through eight
2π-periodic Fourier functions α_i(ϑ):

    H_QBCP = ½ α1 (Px² + Py² + Pz²)
           + α2 (Px X + Py Y + Pz Z)
           + α3 (Px Y − Py X)
           − (1−μ)/R_PE − μ/R_PM
           + α4 X + α5 Y
           − m_S / (α6 · R_PS)

- α1 (kinetic scaling) and α2,α3 (the Coriolis/pulsation terms) generalize the
  CR3BP's rigid rotating-frame structure; in the unperturbed limit α1→1, α3→1,
  α2,α4,α5→0 and the Sun term vanishes, recovering the Earth–Moon CR3BP.
- α6 scales the Sun distance; α7,α8 (in the Frontiers form) carry the Sun's planar
  position. R_PS is particle→Sun distance.
- Symmetry: the α_i are odd/even under (ϑ, x, y, z) → (−ϑ, x, −y, z); the same
  symmetry the CR3BP has plus the Sun phase.

The Sun's perturbation is **O(ε²)** in size (Coriolis and linear-order Sun terms
cancel), so the QBCP is a *small* periodic perturbation of the CR3BP — important for
the verdicts below.

### Sourced parameters (candidate goldens)

From Rosales/Jorba (2023) Table 3 (cross-checked vs Gimeno/Jorba 2018 Table 3;
agree to printed precision):

| Param | Rosales 2023 | Gimeno 2018 | Meaning |
|-------|--------------|-------------|---------|
| μ (Earth–Moon mass ratio) | 0.012150581600000 | 0.012150581623433623 | same as our CR3BP μ_EM |
| m_S (Sun mass, EM units)  | 328900.5423094043 | 328900.54999999906 | Sun mass in Earth+Moon=1 units |
| a_S (Sun semi-major axis) | 388.8111430233511 | 388.81114302335106 | Sun distance in EM-distance units |
| ω_s (Sun synodic freq)    | 0.925195985520347 | 0.92519598551829646 | rad per EM time unit |

T_s = 2π/ω_s ≈ **6.7912** EM (synodic) time units = the **Sun's synodic period in
the Earth–Moon rotating frame ≈ one lunar synodic month (~29.5 d)**. (The
"180 d / 30 d" phrasings in the literature snippets refer to different normalizations
/ different models; in EM-synodic units this ω_s is the ~29.5-d synodic month.)

These four constants + the α_i Fourier table (Andreu 1998, reproduced Gimeno 2018
Table 4 to order k=13) fully specify the model and are **publication-traceable**, so
they qualify as cross-check goldens *if* we ever implement QBCP.

### Sourced periodic-orbit goldens

Rosales/Jorba (2023) Table 4 — dynamical substitutes (the periodic orbits that
replace the collinear points under the Sun perturbation), at t=0, with
y = z = px = pz = 0:

- **L1 substitute (POL1):** x = −0.8369141677649317,  py = −0.8391311559808445
- **L2 substitute (POL2):** x = −1.1556836078332600,  py = −1.1587306159501061

These are exact, sourced, same-model targets: integrate the QBCP vector field one
period T_s from these states and it must close. They are the natural first
golden if a QBCP integrator is ever built. Halo-orbit ICs are *not* tabulated with
full coordinates in the open text (only invariant-curve rotation numbers), so no
halo golden is harvestable from these companions — would need the thesis or a
fuller paper.

## Verdict (2): BCR4BP as a validation rung between CR3BP (V1) and ephemeris (V3)?

**DEFER / background.** Reasoning:

- We already run full DE440 ephemeris (V3). The QBCP sits *strictly below* V3 in
  fidelity (it is a single-period truncation of the real Sun–Earth–Moon motion).
  Anything QBCP can validate, V3 validates more faithfully — so as a *fidelity*
  rung it is largely **redundant** given V3 exists.
- Its one genuine methodological virtue is being **autonomous-adjacent**: a clean,
  reproducible, *closed-form-coefficient* periodic model with published goldens
  (the POL1/POL2 substitutes above). That makes it attractive as a *plumbing /
  regression* check — a deterministic intermediate where a Sun term is present but
  the dynamics are still small-perturbation and exactly reproducible from a paper —
  which the chaotic full-ephemeris V3 is not. This is the recurring "did the Sun
  term get wired in correctly?" question that V1→V3 jumps over.
- But the cost is non-trivial: it needs the full α_i Fourier series (13 orders ×
  8 functions) typed in and a dedicated integrator + frame. For a validation rung
  whose fidelity is dominated by V3, that build is **not justified now**.
- Revisit trigger: if/when a Sun-perturbation defect is ever suspected in the
  CR3BP→ephemeris transition for cislunar work and V3's chaos makes the bug hard to
  localize, the QBCP POL1/POL2 goldens become the cheap deterministic bisection
  point. File it as a known tool, not a current task.

## Verdict (3): QBCP as a SEARCH genome (Track-A type)?

**BACKGROUND (lean defer), with one narrow watch-item.** Reasoning:

- The Sun term is O(ε²) — a *small* periodic perturbation. A CR3BP genome already
  spans the dominant structure of cislunar cyclers; QBCP mostly *refines* those
  solutions rather than opening qualitatively new families a CR3BP genome "can't
  represent." So as a primary discovery genome it is **not** a clear new Track-A.
- Where it *could* matter: solutions whose existence is *gated* by the Sun's
  synodic forcing — i.e. cyclers/quasi-periodic transfers that are resonant with
  T_s (the ~29.5-d synodic month) and have no CR3BP analogue. The QBCP literature's
  own payoff (Andreu; Rosales/Jorba invariant manifolds; the "fast periodic
  transfer orbits extended from EM-RTBP to the Sun–Earth–Moon QBCP" paper,
  CeMDA 2008) is exactly *Sun-synodic-resonant transfer structure*. That is a real
  but **narrow** niche, and our catalogue scope is *repeating cyclers* — Sun-phase
  resonance generically breaks strict periodicity into quasi-periodicity, which sits
  uneasily with the repeating-cycler scope (cf. the S/L resonant-interval memory).
- Net: not a near-term genome. Watch-item — if the discovery campaign ever turns up
  cislunar candidates that *only* close when a synodic-monthly forcing is added
  (CR3BP residual won't close, ephemeris does but origin is opaque), that is the
  signal to prototype a QBCP-forced genome. Until then, **background**.

## Proposed follow-on

- No code, no catalogue writeback now.
- Keep the four parameter constants + POL1/POL2 substitute ICs on record (this note)
  as ready-made goldens.
- Costless future step if ever wanted: harvest the full α_i Fourier table from
  Gimeno/Jorba (2018) Table 4 (open access) so a QBCP integrator could be stood up
  from published values alone, no thesis access needed.
- Honest status: thesis filed but not in-session-readable; model + goldens digested
  from Gimeno/Jorba (2018, Frontiers) and Rosales/Jorba (2023, CeMDA 135:15), both
  citing Andreu (1998) as primary.
