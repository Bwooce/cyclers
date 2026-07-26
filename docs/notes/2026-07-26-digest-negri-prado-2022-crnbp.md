# Digest: Negri & Prado 2022, "Circular Restricted n-Body Problem" (CRNBP)

**Paper:** "Circular Restricted n-Body Problem"
**Venue:** *Journal of Guidance, Control, and Dynamics* 45(7), 1357-1364 (2022), DOI
`10.2514/1.G006430` (Crossref-confirmed: `published-print` 2022-07, vol 45, issue 7,
pages 1357-1364)
**Authors:** Rodolfo Batista Negri (Ph.D. candidate, INPE) and Antônio F. B. A. Prado
(INPE, also RUDN University), National Institute for Space Research (INPE), São José
dos Campos, Brazil
**Free full text:** arXiv:2307.10881 [eess.SY], "A Circular Restricted n-body Problem"
(same authors, posted 2023-07-20 — a minor title variant of the journal version, textually
identical derivation/results as far as can be checked; the journal appears to have dropped
the article "A" from the title). The arXiv posting postdates the nominal JGCD print date by
about a year — plausible under an author-rights embargo delay, not evidence of a different
paper.
**Filed:** `negri-prado-2022-circular-restricted-n-body-problem-jgcd-doi-10.2514-1.G006430-arxiv-2307.10881.pdf`
(in the private `cyclers_pdf` corpus, commit `8833b60`)
**Acquired/digested:** 2026-07-26 (`#712`). Text-layer PDF (arXiv source), no OCR needed.
JGCD/ARC page itself is paywalled (403 on direct fetch); full text obtained via the free
arXiv preprint instead.

**Why in corpus (`#712`):** found alongside `#711`'s Gilliam & Bettinger 2024 Jovian
CRNBP paper while scoping a possible next 5-body/N-body discovery-strategy extension
beyond this project's CCR4BP arc (`#686`->`#708`). Task was to resolve the citation
(only a live-search snippet existed) and determine whether it is a GENERAL N-body
framework paper, as opposed to `#711`'s Jovian-specific one.

**Title correction:** the coordinating session's tentative citation ("Circular Restricted
n-Body Problem," JGCD, DOI 10.2514/1.G006430) is confirmed EXACTLY correct — no
correction needed. Authors/year were previously unknown; now confirmed as above.

## 1. Confirmed: this IS the general framework paper, `#711` is a specific application of it

This is unambiguously a general, system-agnostic N-body EQUATIONS-OF-MOTION framework
paper, not tied to any specific real moon system. It:

- Derives the EOM for an N-body system in fully general form (arbitrary N, arbitrary
  masses `M_1...M_{N-1}` plus a massless N-th body) before imposing any Jovian- or
  other-system-specific numbers.
- Only in its "Illustrative Examples" section (§3) applies the general model to worked
  cases, one of which happens to be Jovian (Jupiter-Ganymede-Europa-Io, and separately
  Jupiter-Europa-Io-Ganymede-Callisto) — used purely to demonstrate the framework, not as
  the paper's subject.
- Is explicitly positioned by the authors as the generalization of their own EARLIER
  paper, Negri & Prado 2020, "Generalizing the Bicircular Restricted Four-Body Problem"
  (JGCD 43(6), 1173-1179, DOI `10.2514/1.G004848`) — the N=4 special case. This 2020
  paper is cited in the References as ref [12] and its generalized two-body-approximation
  reasoning is the direct ancestor of the N-body derivation here. (The 2020 paper is not
  yet in this project's corpus; flagging as a candidate future acquisition if a CRNBP
  build is ever undertaken — its correction to Huang's 1960 BCR4BP indirect term is
  exactly the term this paper's Eq. (11) generalizes to N bodies.)

Cross-check against the live-search snippets found during `#712`'s own resolution: Annika
J. Gilliam's AFIT thesis "Extension of the Circular Restricted N-Body Problem (CRNBP) to
Varying [...]" and a related ScienceDirect paper on Jovian debris propagation "using the
circular restricted 3- and N-body problems" both cite/build on this Negri & Prado line —
confirming `#711`'s Gilliam & Bettinger 2024 paper is a DOWNSTREAM, Jovian-system-specific
APPLICATION of (and explicitly references) this more general framework, not an
independent parallel line of work. **Relationship to `#711`: related, not independent —
Gilliam & Bettinger cite Negri & Prado's CRNBP formulation and specialize/apply it to the
real Jovian system.** (Full confirmation of the citation, and any structural differences
in exactly how Gilliam & Bettinger's derivation departs from Negri & Prado's, is `#711`'s
own job once it digests the Gilliam & Bettinger paper directly — `#711` had not yet closed
at the time this digest was written.)

## 2. The EOM derivation (§2.1)

Setup: N bodies total, `M_1, ..., M_{N-1}` massive + the N-th body massless (the
"restricted" body/spacecraft). Inertial-frame N-body EOM (Eq. 1) is specialized to the
restricted case (Eq. 2), then re-centered on the M1-M2 barycenter (CM12) by algebraically
eliminating the barycenter's own acceleration (Eqs. 3-8) — this is the step that produces
the "indirect" cross terms.

Constraints imposed (identical in spirit to BCR4BP, extended to N bodies):
- `M_1` and `M_2` (the "primaries") describe a circular orbit about their common
  barycenter (CM12) — standard CR3BP primaries.
- Every OTHER massive body `M_3, ..., M_{N-1}` is assumed to move on a CIRCULAR orbit
  CENTERED ON `M_1` (not on CM12) — this is the "concentric" idealization, matching this
  project's own `ccr4bp.py` docstring's "concentric circular" framing for Ganymede.
- All orbits are COPLANAR.
- The synodic frame rotates at `M_1`-`M_2`'s own angular rate (`n12 = 1` in canonical
  units); each other body's synodic phase advances analytically per Eq. (12):
  `psi_j(t) = psi_j0 + (n_j - n12) t` — i.e., each extra body's synodic angular rate is
  its OWN inertial mean motion minus the frame rate, exactly matching `ccr4bp.py`'s
  `omega_gan = n3/n2 - 1` convention for Ganymede.

Final EOM (Eq. 11a-c), for the massless body's position `(x, y, z)` in the synodic frame,
with `mu_1 = 1-mu`, `mu_2 = mu` the usual CR3BP primaries and `mu_j = M_j/(M1+M2)` for
each extra body `j = 3, ..., N-1`:

```
x'' = 2y' + x - mu1(x+mu2)/r1^3 - mu2(x-mu1)/r2^3
      - sum_j mu_j [ (x+mu2-Rj*cos(psi_j))/rj^3
                      - sum_{k != j} mu_k (Rj*cos(psi_j) - Rk*cos(psi_k)) / (Rk^2+Rj^2-2*Rk*Rj*cos(psi_k-psi_j))^{3/2} ]
y'' = -2x' + y - mu1*y/r1^3 - mu2*y/r2^3
      - sum_j mu_j [ (y-Rj*sin(psi_j))/rj^3
                      - sum_{k != j} mu_k (Rj*sin(psi_j) - Rk*sin(psi_k)) / (...)^{3/2} ]
z'' = -mu1*z/r1^3 - mu2*z/r2^3 - sum_j mu_j*z/rj^3
```

where `rj` is the distance from the spacecraft to body `j`, and `Rj`, `psi_j` are body
`j`'s (fixed) orbital radius and (time-varying) synodic phase angle about M1.

**The load-bearing structural point:** for each extra body `j`, its forcing term carries
an INNER SUM over every OTHER extra body `k != j` — this is the "indirect" term
correcting for `M_k`'s own pull on `M_j`'s barycentric position (mutual coupling among
the extra bodies' contribution to the frame's own acceleration, NOT a direct
body-body force). The paper's own §2.1 flags this explicitly: "The difference between
Eqs. (11) and the equations derived by Iuliano [24, 2016] is in the summation inside the
brackets" — Iuliano's earlier (unpublished/thesis) N-body attempt DROPPED this term. For
`N=4` (a single extra body, e.g. just Ganymede) the inner sum has only `k in {1,2}` (i.e.
M1, M2 themselves), which reduces EXACTLY to the "conventional" BCR4BP's known
direct+indirect Sun/third-body term — confirmed explicitly in the paper's own text: "Eqs.
11 reduce to the 'conventional' BCR4BP ... if N=4 and R3 >> 1."

## 3. Comparison against this project's `src/cyclerfinder/core/ccr4bp.py`

Read directly (`ccr4bp.py`, Jupiter-Europa-Ganymede CCR4BP): this project's own model has
THREE massive bodies (Jupiter=M1, Europa=M2, Ganymede=M3) plus the massless spacecraft —
in the paper's counting that is exactly `N=4` (a single EXTRA body beyond the M1/M2
primary pair). At `N=4` the paper's general inner-sum correction reduces to the ordinary
BCR4BP-shape direct+indirect term, which is precisely the form `ccr4bp.py`'s own
docstring describes ("Standard incoherent direct + indirect Ganymede acceleration (same
shape as bcr4bp's Sun term)"). **So `ccr4bp.py` is a correctly-structured N=4 instance of
this paper's general CRNBP — not a different/incompatible model.** It has never yet been
exercised at N>=5 (two or more simultaneous extra perturbers), which is exactly where the
paper's genuinely-new mutual inner-sum term (`sum_{k!=j}`, nontrivial for >=2 extra
bodies) first becomes active and where `#711`'s Gilliam & Bettinger Jovian application
(Io AND Ganymede simultaneously forcing the Jupiter-Europa CR3BP, N=5) actually lives.
**Practical implication for a future 5-body build:** extending `ccr4bp.py` to N=5 is NOT
a fresh derivation — it is adding exactly one more `mu_j, R_j, psi_j` term plus its
inner-sum coupling to the existing structure, per Eq. (11) here. This paper is therefore
directly reusable as the derivation source for that extension, independent of whichever
paper (this one or `#711`'s) ends up chosen as the primary citation.

## 4. Worked numeric examples (§3) — reusable as positive controls

1. **FLI chaos-map reproduction (structural-stability check, not a numeric target):**
   reproduces Todorović, Wu & Rosengren 2020 (*Science Advances* 6, eabd1313)'s solar-system
   Fast Lyapunov Indicator "arches of chaos" map (Venus-Neptune CRNBP, Sun=M1, Jupiter=M2,
   test particles near Jupiter's orbital plane, ephemerides epoch 2012-09-30 00:00:00 TDB,
   100-year integration) — described as "practically perfect" qualitative match, no digit-
   grade table given. Useful only as a qualitative sanity check, not a reproducible
   digit-grade positive control.
2. **Jupiter-Ganymede vertical Lyapunov-orbit continuation (Figs. 4-5):** a family of
   vertical Lyapunov orbits about L3 in the Jupiter-Ganymede CR3BP, continued via a
   homotopy parameter `epsilon` (Eq. 28, `epsilon=0` -> pure CR3BP, `epsilon=1` -> full
   CRNBP) into a CRNBP for Jupiter-Ganymede-Europa-Io with Io/Europa periods locked to
   exactly 1/4 and 1/2 of Ganymede's (idealized Laplace resonance). No numeric IC table
   given (figure-only) — not directly digit-reproducible without digitizing Fig. 5, but
   the epsilon-homotopy CONTINUATION METHOD itself (CR3BP periodic orbit -> CRNBP via
   parameter continuation) is a directly reusable technique for extending any of this
   project's existing CR3BP/CCR4BP families to N>=5.
3. **Europa-landing trajectory richness (Figs. 6-8):** CR3BP vs. CRNBP (Jupiter-Europa-
   Io-Ganymede-Callisto) backward-propagated trajectories from 50 km above Europa's
   surface, IC epoch 2016-04-09 00:00:00 TDB, Jacobi-constant-derived initial speed at
   the CR3BP's L2 value (Eq. 30). Demonstrates new transeuropa trajectory windows opening
   at specific arrival-time offsets (theta ~150 deg at +2.36 days, theta ~350 deg at +7.1
   days) that do NOT exist in the CR3BP alone — a genuine "no 3-body analog" qualitative
   claim, but again figure-only, no digit table.

**None of the three examples ships an IC/parameter table precise enough for a digit-grade
golden-test positive control** (per `[[feedback_golden_tests_sourced_only]]`) — all three
are qualitative/figure-based demonstrations. If a future N=5+ build wants a genuine
positive control, the ephemerides-correspondence machinery in §2.2 (converting a real
ephemerides configuration into CRNBP `(R_j, psi_j0, n_j)` parameters — Eqs. 13-24) is the
more directly useful, fully-specified piece: it is a complete, reproducible ALGORITHM
(not a numeric result) for grounding a CRNBP instance to a real system, which this
project's own `ccr4bp.py` currently does via its `jupiter_europa_ganymede_default()`
JPL-SSD-registry route instead — worth comparing the two conventions if this is ever
built, since they are NOT obviously identical (this paper projects onto M2's mean orbital
plane and roots initial phases off `s12` vectors; `ccr4bp.py`'s route was not checked
against this in detail here).

## 5. Bearing on a future CCR5BP/CRNBP discovery-strategy pass

- Confirms the "grounding paper" role floated when `#711`/`#712` were dispatched: this
  paper is the correct GENERAL-FRAMEWORK citation, `#711`'s Gilliam & Bettinger 2024 is
  the correct REAL-JOVIAN-SYSTEM APPLICATION citation, and they are a coherent PAIR (one
  general, one specialized, with the specialized one citing the general one) rather than
  two competing or redundant leads.
- Recommend citing BOTH in any future N>=5-body discovery-strategy scoping note: this
  paper for the EOM derivation and general N-body structure (directly reusable to extend
  `ccr4bp.py`), `#711`'s paper for the real-Jovian-system worked numbers and any
  discovery-relevant orbit/manifold results it reports (not yet assessed here — `#711`
  had not closed at the time of this digest).
- The `#686` "genuinely N-body-native, no 3-body analog" bar (same bar established for
  CCR4BP by the Kumar et al. digests) is echoed qualitatively here too (§3's transeuropa-
  window example), but not established with the same rigor (no secondary-resonance /
  Chirikov-overlap-style mechanism analysis as in Kumar et al. 2023) — this paper is a
  MODEL/methods paper, not a discovery-mechanism paper. A future N=5 discovery pass would
  still need its OWN structural-novelty argument, not inherit this paper's.
