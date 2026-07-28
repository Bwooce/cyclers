# Digest: Anderson & Lo 2010 (planar Europa orbiter) + Anderson & Lo 2011 (resonant flybys, ballistic case) (#745)

**Task:** `#745`, continuing the `#730` consolidated acquisition backlog
(`docs/notes/2026-07-27-730-acquisition-backlog-master-list.md`). The user
directly supplied PDFs for §3 items 12 and 13 of that list — both
independently verified page-1-exact against the master list row by the
coordinating session before this task started.

**Filed** (private `cyclers_pdf` repo, commit `bd6e4b9`):
1. `anderson-lo-2010-dynamical-systems-planetary-flybys-approach-planar-europa-orbiter-jgcd-33-6-1899-doi-10.2514-1.45060.pdf`
2. `anderson-lo-2011-dynamical-systems-resonant-flybys-ballistic-case-jas-58-167-doi-10.1007-BF03321164.pdf`

**OCR status: both text-layer, no OCR needed.** Both are Acrobat-Distiller-
produced (native LaTeX/Word PDFs, downloaded from ARC/Syracuse and JAS
mirrors respectively); `pdftotext -layout` yielded 203,964 / 95,619 chars
over 14 / 28 pages — far above the 10-char/page floor.

---

## 1. Anderson & Lo 2010 — "Dynamical Systems Analysis of Planetary Flybys and Approach: Planar Europa Orbiter"

*JGCD* 33(6):1899–1912. `#730` §3 item 12 (Medium-high, recurring — 2
independent flags: `#727`, `#728`-oberon-survey). DOI `10.2514/1.45060`
(confirmed, matches this PDF's own header exactly).

**Method:** takes the last Ganymede-to-Europa segment of the *actual*
(full-ephemeris) Europa Orbiter (PEO) trajectory designed at JPL in 1999
(Johannesen & D'Amario), transfers it stepwise into the planar circular
restricted three-body problem (PCRTBP, Jupiter-Europa, `μ` not stated
numerically in this paper — see the 2011 companion for the explicit value)
via JPL's LTool differential corrector, and analyzes it against **unstable
resonant periodic orbits** and their invariant manifolds. The PEO trajectory
is divided into three energy segments (`C_i`, `C_m`, `C_f`) by its two
deep-space `ΔV`s (Table 1: `ΔV1`=118.4 m/s, `ΔV2`=91.6 m/s, ephemeris model).
Unstable 3:4 and 5:6 resonant orbits are computed at each segment's Jacobi
constant via single-shooting + linear-extrapolation continuation (the
Howell-Breakwell symmetric algorithm, modified for resonant rather than
libration-point orbits), and their stable/unstable manifolds are globalized
by offsetting ~`1e-6`-`1e-5` along the monodromy eigenvector and integrating.
Poincaré sections (`y=0` on the far side from Europa) are plotted using both
`(x, ẋ)` and Delaunay-like `(L, g̃)` coordinates (`L=√a`, `g̃`=argument of
periapse) — the latter more clearly separates nearby resonances.

**Key result:** the PEO trajectory's own Poincaré-section crossings
(`PEO_1`, `PEO_2`, `PEO_3`) lie *directly on* the stable manifolds of the 5:6
and 3:4 resonant orbits at the middle/final energy levels (Figs. 12–20) — the
final Europa approach specifically uses the **stable manifold of an `L_2`
Lyapunov orbit** (Fig. 15/19), and the paper's headline finding is that a 5:6
resonance (not 3:4) is dynamically necessary as the last resonance before
Europa orbit insertion (EOI), because the 5:6 resonant orbit is shown to lie
**entirely within** the intersection of the `L_2` Lyapunov orbit's own stable
and unstable manifolds (Fig. 17) — i.e. the resonant periodic orbit is itself
homoclinic-shadowed by the libration-point orbit's manifold tube. This is
presented as a *post-hoc dynamical-systems explanation* of an already-flown
mission-design trajectory, not a new trajectory construction.

## 2. Anderson & Lo 2011 — "A Dynamical Systems Analysis of Resonant Flybys: Ballistic Case"

*J. Astronaut. Sci.* 58(2):167–194. `#730` §3 item 13 (Medium-high,
recurring — 2 independent flags: `#727`, `#728`-oberon-survey). **DOI
`10.1007/BF03321164`** — this project's own `CORPUS_INDEX.md`/master-list row
already carried this exact DOI, CrossRef-confirmed by `#743`/`#744`; this
task independently confirms it is the correct paper via page-1-exact title/
author/volume/page match against this PDF.

Explicit CRTBP parameter given here (absent from the 2010 paper): **Jupiter-
Europa `μ = 2.5266448850435 × 10⁻⁵`** (Eq. after "Circular Restricted
Three-Body Problem", p.169) — see the cross-check below.

**Method:** the companion, self-contained study the 2010 paper's own
Conclusions promised ("additional work will focus on refining these
constraints even further"). Rather than analyzing an existing
full-ephemeris mission trajectory, this paper **constructs from scratch** a
purely ballistic (no maneuvers) periodic-in-energy trajectory that cycles
continuously between the 3:4 and 5:6 resonances via two flybys of Europa —
built first as a two-body patched-conic guess (rotating `V∞` at each flyby to
change resonance, Fig. 2/"Resonance Transition" section), then converted to a
truly continuous CRTBP trajectory via a **multiple-shooting differential
corrector** (JPL LTool) applied iteratively: apoapsis + periapsis patchpoints
of both two-body legs, refined over four repeated passes through the
resonance-cycle scenario until the initial/final states converge to
**0.091 km position / 5.3×10⁻⁶ km/s velocity** difference (deemed "the limit
of the numerical accuracy of the computer over such a long integration").
The converged flyby trajectory is then analyzed exactly as in the 2010 paper
(unstable resonant-orbit families named by inner/outer-region + syzygy-loop
taxonomy: 3:4-LO/LI, 5:6-LO/LI/NO/NI; invariant manifolds via monodromy
eigenvector offset + Poincaré sections).

**Key result:** the flyby trajectory's two Poincaré-section intersections lie
almost exactly on the intersection points of the 3:4 and 5:6 resonant
orbits' own invariant manifolds (offsets ~`8×10⁻⁵` and ~`2×10⁻⁴` in `x`) —
i.e. **the ballistic flyby trajectory is itself closely shadowing a
homoclinic/heteroclinic connection** between the two resonances. The paper
then explicitly **computes** a homoclinic trajectory of the 3:4 orbit (Table
2 intersection state) and a heteroclinic trajectory connecting the 3:4 and
5:6 orbits (Table 3 intersection state) by interpolating the two nearest
manifold points and integrating forward/backward — the first explicit
mission-design-oriented **construction** (not just observation) of a
homoclinic/heteroclinic connection between two *resonant* (not
libration-point) periodic orbits in this specific Jupiter-Europa system. The
Conclusion frames this as a genuine three-body analogue of the two-body
"rotate `V∞` at each flyby" resonance-hopping technique.

---

## Cross-check: does this project's own Jovian-tour/resonant-flyby machinery already use, contradict, or omit anything these two papers establish?

**Headline finding: these two papers are the direct lineage ancestor of the
"Jupiter-Europa 3:4 resonant orbit" object class already sitting in this
project's own catalogue.** `data/catalogue.yaml`'s N=5 CRNBP torus row
(`europa-3-4-crnbp-torus-jupiter-2026`, task chain `#714`→`#736`) is
explicitly built by continuing "the planar Jupiter-Europa 3:4 resonant orbit
of Kumar et al. 2021" into a Laplace-locked 5-body torus, per its own
provenance block. Reading `kumar-anderson-delallave-gunter-2021-europa-
ganymede-resonant-orbits-ccr4bp-AAS-21-651-arxiv-2109.14815.pdf` (already in
this project's corpus, co-authored by this same Rodney L. Anderson) confirms
its own Introduction opens by citing **exactly these two papers as its refs
[1]/[2]** — "studies... focused on the tour endgame... final approach to
Europa after a series of flybys of the Galilean moons using either ballistic,
impulsive, or low-thrust trajectories.¹⁻³" (ref [3] is the still-unacquired
`#730` §3 item 14, "Flyby Design using Heteroclinic and Homoclinic
Connections of Unstable Resonant Orbits") — and Kumar 2021's own seed object
is described using the **same "m:n resonant periodic orbit" identification
convention and the same "exterior" (vs. interior) terminology** these two
papers use for their own 3:4/5:6 orbit taxonomy (Anderson-Lo's "3:4-LO" =
Loop, **O**uter region; Kumar's own abstract: "exterior Jupiter-Europa...
resonant periodic orbits"). **This strongly suggests, but this pass does NOT
independently verify digit-for-digit**, that Kumar 2021's seed 3:4 orbit
(and therefore this project's own N=5 torus discovery two steps downstream
of it) is the *same* resonant-orbit family Anderson & Lo compute and name
"3:4-LO" in the 2011 paper — no IC/period/Jacobi-constant comparison was run
this pass (Anderson & Lo 2011 gives no numeric IC table for its 3:4-LO
family, only continuation plots, Fig. 9/10; Kumar 2021's own Table would need
to be pulled separately). **Flagged as a genuine, concrete, unexecuted
future cross-check** (comparable in kind to `#742`'s flagged-but-unexecuted
Franz-Russell/Casoliva symmetric-orbit cross-check): if Kumar's 3:4 orbit
family traces back to Anderson & Lo's own computed 3:4-LO family (plausible
given the shared naming convention and Kumar's own citation framing), this
project's N=5 torus discovery has a documented two-hop lineage all the way
back to these two 2010/2011 papers, not just to Kumar 2021.

**Secondary finding — this project's Jovian moon-tour CODE is
patched-conic/Tisserand-graph two-body throughout; no module implements the
three-body resonant-orbit-invariant-manifold or homoclinic/heteroclinic-
connection machinery these two papers build (as distinct from the CATALOGUE
row above, which is a Kumar/Anderson torus reproduction, not code built from
these two papers directly).** Confirmed by direct code inspection:

- `src/cyclerfinder/core/flyby.py` — the project's core flyby-mechanics
  module — is explicitly a "patched-conic gravity-assist" module (its own
  docstring: "Pure functions on heliocentric V∞ vectors... Strange & Longuski,
  JSR 2002"), i.e. exactly the **two-body approximation** Anderson & Lo's own
  2011 Background section frames as the *predecessor* their CRTBP
  dynamical-systems method goes beyond ("Although the gravity assist
  techniques have been found to provide very accurate conditions for mission
  design, they still use two-body approximations for what is naturally a
  three-body trajectory").
- `src/cyclerfinder/search/resonant_conic.py` (the Hernandez 2017 EGGIE
  triple-cycler seed generator) and `src/cyclerfinder/search/ieg_seed.py` are
  likewise explicitly **conic/Lambert-based**, per their own docstrings — the
  same patched-conic-with-resonance-matching paradigm, not a three-body
  manifold-shadowing construction.
- Grepping `src/cyclerfinder` for "homoclinic"/"heteroclinic" outside the
  CCR4BP/CRNBP torus-manifold modules (`ccr4bp_manifold_globalize.py`,
  `ccr4bp_heteroclinic_search.py`, `crnbp_torus_ghost_guard.py`) returns
  nothing — this project's ONLY homoclinic/heteroclinic-connection machinery
  targets **quasi-periodic TORI** in the time-periodic CCR4BP/CRNBP models
  (the Kumar/Anderson/de la Llave lineage, `#689`-`#736`), a structurally
  different object class from these two papers' **classical periodic
  resonant orbits in the autonomous planar CR3BP** (Jupiter+Europa only, no
  third perturbing moon, no time-periodic forcing). The manifold-globalization
  *pattern* (monodromy eigenvector offset + forward/backward integration) is
  the same general idea both places use, but no code here computes a
  classical-CR3BP resonant-periodic-orbit family, its monodromy matrix, or
  its Poincaré-section manifold intersections the way these two papers do.
- **Consequence for this project's own domain:** this is not a bug — the
  project's existing Jovian-tour work (EGGIE, VEM/IEG triple cyclers,
  Russell-Strange moon cyclers) targets *free-return cyclers*, a different
  design goal from Anderson & Lo's *single-mission resonance-transition/
  endgame-approach* problem — but it is a genuine, documented capability gap:
  if this project ever pursues Anderson-Lo-style Galilean-endgame or
  resonance-hopping cycler variants (e.g. refining `resonant_conic.py`'s
  patched-conic EGGIE seed into a true CRTBP resonant-manifold seed, the way
  Anderson & Lo's own 2011 paper explicitly recommends as future mission-
  design use — "Potential Techniques" section), no infrastructure currently
  exists for it.

**Mass-ratio cross-check (quantitative, in the spirit of `#742`'s
Franz-Russell `μ` comparison):** this project's own Jupiter-Europa mass
ratio, computed from `core/satellites.py`'s Europa `GM=3202.739 km³/s²`
(JPL SSD) and `core/constants.py`'s Jupiter *system* `GM=1.267127641×10⁸
km³/s²` (JPL DE440, Park et al. 2021) as `μ = GM_Europa/(GM_Jupiter +
GM_Europa) = 2.5274944×10⁻⁵`, differs from Anderson & Lo's own stated
`μ=2.5266448850435×10⁻⁵` (2011 paper, p.169) by **~0.034% relative** (5th
significant figure) — small, consistent with a GM-source vintage difference
(Anderson & Lo 2010/2011 predate DE440; the exact GM source they used is not
cited in either paper). Same caveat class as `#742`'s Franz-Russell finding:
non-zero but negligible at the precision either paper actually needs; no
correction implied to either paper or this project's own constants.

**No existing catalogue row or prior digest references either paper's
specific 3:4/5:6 Europa resonant-orbit families, homoclinic trajectory (Table
2), or heteroclinic trajectory (Table 3) states** — grepped
`data/catalogue.yaml` for "3:4" and "5:6" Jovian/Europa context and found no
match; these are genuinely new, unadopted results, not contradicted or
duplicated by anything currently in the catalogue.

---

## Mandatory citation-mining pass

Both papers share a near-identical core bibliography (the 2011 paper's list
of 38 is a superset of the 2010 paper's 31, plus historical flyby-theory
references 1–8: Lexell/Laplace/Le Verrier/D'Alembert/Tisserand, 1770s–1896,
and 1960s JPL memos Minovitch/Sturms&Cutting/Flandro — pre-modern-era
historical background, not flagged for acquisition). Read in full (both
reference lists, both papers' own citing text).

**Already in corpus / already flagged elsewhere in `#730`, confirmed
recurring (no new action):**
- Koon, Lo, Marsden & Ross 2000, "Heteroclinic Connections between Periodic
  Orbits and Resonance Transitions in Celestial Mechanics," *Chaos*
  10(2):427–469 [2011 ref 19] — **already in corpus**
  (`koon-lo-marsden-ross-2000-heteroclinic-connections-resonance-transitions-chaos-10-2.pdf`,
  digested, `#314` control anchor).
- Barrabés & Gómez 2003, "Three-Dimensional p-q Resonant Orbits Close to
  Second Species Solutions" [2010 ref 24, 2011 ref 30] — **already acquired**
  by `#744` (`#730` §4 item 26); this is an independent 3rd citing context
  confirming the item's relevance to the broader resonant-orbit lineage,
  beyond the Casoliva/Barrabés-Gómez p-q cluster `#742`/`#744` already traced.
- Anderson & Lo 2009, "Role of Invariant Manifolds in Low-Thrust Trajectory
  Design," *JGCD* 32(6):1921–1930 [2010 not cited directly under this exact
  title but subsumes refs 21/22/23 below; 2011 ref 24] — `#730` §3 item 15,
  already flagged, DOI confirmed.
- Howell, Marchand & Lo 2001, "Temporary Satellite Capture of Short-Period
  Jupiter Family Comets from the Perspective of Dynamical Systems," *JAS*
  49:539–557 [2011 ref 20] — **now a 2nd independent flag**: already surfaced
  as a "tangential dynamical-systems-capture parallel" (low-medium priority)
  in `#742`'s Franz-Russell citation-mining pass (its ref 33); this task's
  independent citation confirms the recurrence.
- Johannesen & D'Amario 1999, "Europa Orbiter Mission Trajectory Design,"
  AAS 99-360 [2010 ref 15, 2011 ref 35] — **now a 3rd independent flag**
  (also surfaced, low-medium priority, in `#742`'s Anderson-2021
  citation-mining pass, its ref 23) — this is literally the source mission
  trajectory the 2010 paper analyzes, so the recurrence is expected, not
  newly informative, but worth recording as the strongest-corroborated
  candidate to emerge from this specific task.
- Meyer & Hall, *Introduction to Hamiltonian Dynamical Systems and the
  N-Body Problem* textbook [2010 ref 29] — already flagged low-medium
  priority by `#742`'s BMO 2009 mining pass.
- Wiggins, *Introduction to Applied Nonlinear Dynamical Systems and Chaos*
  textbook [2011 ref 33] — general dynamical-systems textbook, not
  domain-specific; not separately flagged (background, consistent with
  `#730` §9 triage convention for textbook references).

**Genuinely new candidates (flagged, NOT acquired):**
- *Medium:* Lo, Anderson, Whiffen & Romans 2004, "The Role of Invariant
  Manifolds in Low Thrust Trajectory Design (Part I)," AAS 04-288; Anderson
  & Lo 2004, "The Role of Invariant Manifolds in Low Thrust Trajectory Design
  (Part II)," AIAA 2004-5305; Lo, Anderson, Lam & Whiffen 2006, "The Role of
  Invariant Manifolds in Low Thrust Trajectory Design (Part III)," AAS 06-190
  [2010 refs 11/12/14, 2011 refs 21/22/23] — the three conference-paper
  precursors to the already-flagged/DOI-confirmed JGCD 2009 journal version
  (`#730` §3 item 15); likely superseded/subsumed by that final version per
  standard AAS-conference-to-JGCD-journal practice in this lineage (the same
  pattern `#742` confirmed for Campagnola-Russell 2009→2010), so flagged only
  at low priority — acquiring the JGCD 2009 version first would likely make
  these redundant.
- *Medium:* Anderson, R.L. 2005, "Low Thrust Trajectory Design for Resonant
  Flybys and Captures Using Invariant Manifolds," Ph.D. Dissertation,
  University of Colorado at Boulder (`ccar.colorado.edu/~rla/papers/
  andersonphd.pdf`, per the 2011 paper's own ref 26 URL) — the foundational
  dissertation underlying this entire Anderson/Lo lineage (both papers'
  Introductions cite it directly as "Anderson [13]"/"[26]" for "additional
  details related to this series of papers"); likely freely available at the
  cited URL (a university advisor page), not yet checked/fetched this pass.
- *Low-medium:* Belbruno & Marsden 1997, "Resonance Hopping in Comets,"
  *Astronomical Journal* 113(4):1433–1444, DOI `10.1086/118359` [2010 ref 9,
  2011 ref 17] — directly on-point (comet resonance-hopping, the historical
  motivation both papers cite for the whole "why resonant flybys matter"
  framing); Bollt & Meiss 1995, "Targeting Chaotic Orbits to the Moon Through
  Recurrence," *Physics Letters A* 204:373–378, DOI `10.1016/0375-9601(95)
  00502-T` [2010 ref 1, 2011 ref 15] — the paper both papers' Introductions
  credit as the founding suggestion to apply dynamical-systems techniques to
  three-body flyby/tour design; Lo & Ross 1998, "Low Energy Interplanetary
  Transfers Using Invariant Manifolds of L1, L2, and Halo Orbits," AAS/AIAA
  Space Flight Mechanics Meeting, Monterey, Feb 1998 [2010 ref 10, 2011 ref
  18] — distinct from the already-corpus Koon-Lo-Marsden-Ross-1999 GENESIS
  AAS paper (different title/venue/year); Schroer & Ott 1997, "Targeting in
  Hamiltonian Systems that have Mixed Regular/Chaotic Phase Spaces," *Chaos*
  7(4):512–519 [2010 ref 3, 2011 ref 16] — a second early resonance-transfer-
  time-reduction reference in the same "recurrence-based lunar targeting"
  lineage as Bollt & Meiss.
- *Low:* Wilson, R.S. 2003, "Derivation of Differential Correctors Used in
  GENESIS Mission Design," IOM 312.I-03-002, JPL [2010 ref 17, 2011 not
  separately cited but same technique] — internal JPL memo, likely
  restricted-access, low acquisition priority; Pernicka, H.J. 1986, "The
  Numerical Determination of Lissajous Orbits in the Circular Restricted
  Three-Body Problem," M.S. Thesis, Purdue University [2010 ref 18, 2011 ref
  37] — thesis, tool/methods reference not a domain result; Masdemont, J.
  and Mondelo, J.M., "Notes for the Numerical and Analytical Techniques
  Lectures (Draft Version)," Advanced Topics in Astrodynamics Summer Course,
  Barcelona, July 2004 [2010 ref 27, 2011 ref 32] — unpublished course notes,
  low acquisition priority; Anderson, R.L., Lo, M.W. & Born, G.H. 2003,
  "Application of Local Lyapunov Exponents to Maneuver Design and Navigation
  in the Three-Body Problem," AAS 03-569 [2011 ref 38] — the source of the
  local-Lyapunov-exponent stability diagnostic used briefly in the 2011
  paper's Fig. 8, tangential to this project's own search-method domain
  (no local-Lyapunov-exponent machinery exists anywhere in `cyclerfinder`);
  Ludwinski, Guman, Johannesen, Mitchell & Staehle 1998, "The Europa Orbiter
  Trajectory Design," IAF-98-Q.2.02 [2010 ref 20] — an earlier-stage
  companion piece to the already-flagged Johannesen & D'Amario 1999 mission
  paper (item above); historical/mission-specific.
- Historical pre-1900s comet-orbit-theory references (Lexell 1770/1772/1779,
  Laplace 1805, Le Verrier 1857, D'Alembert 1773, Tisserand 1889/1896) and
  1960s JPL internal memos (Minovitch 1961/1963, Sturms & Cutting 1966,
  Flandro 1966) [2011 refs 1–12] — foundational flyby-theory prehistory, not
  flagged (out of scope for this project's search-method domain, consistent
  with `#730` §9 triage convention).
- Szebehely 1967, Roy 1988 (*Orbital Motion*), Murray & Dermott 1999 (*Solar
  System Dynamics*), Parker & Chua 1989 — general dynamical-systems/orbital-
  mechanics textbooks, already in corpus or already covered by prior digests'
  textbook triage; not separately flagged.

No genuinely surprising high-priority gap emerged beyond the Anderson 2005
PhD dissertation (the clear root of the whole lineage, now flagged) and the
recurring Howell-Marchand-Lo 2001 / Johannesen-D'Amario 1999 items, both now
independently corroborated a 2nd/3rd time within this project's own
citation-mining history.

---

## Summary answers (for the dispatching session)

- **Paper 1 (Anderson & Lo 2010):** post-hoc dynamical-systems explanation of
  the real, previously-flown 1999 Europa Orbiter (PEO) trajectory design —
  its final Europa approach is shown to closely follow the stable manifold of
  an `L_2` Lyapunov orbit, and the last resonance (5:6, not 3:4) is dynamically
  required because the 5:6 resonant orbit lies entirely inside the
  intersection of the Lyapunov orbit's own manifolds.
- **Paper 2 (Anderson & Lo 2011):** the self-contained ballistic companion —
  constructs (multiple-shooting) a purely ballistic 3:4↔5:6 resonance-cycling
  Europa-flyby trajectory from scratch, shows it closely shadows a
  homoclinic/heteroclinic connection between the two resonant-orbit
  families, and explicitly computes both a homoclinic (Table 2) and
  heteroclinic (Table 3) connection trajectory as concrete mission-design
  artifacts. States the Jupiter-Europa `μ=2.5266448850435×10⁻⁵` explicitly.
- **Cross-check verdict:** two findings. (1) **Lineage**: these two papers
  are cited as the direct tour-endgame ancestry (refs [1]/[2]) by Kumar-
  Anderson-de la Llave-Gunter 2021, whose own "Jupiter-Europa 3:4 resonant
  orbit" is the literal seed this project's own catalogued N=5 CRNBP torus
  discovery (`#714`-`#736`) continues — a plausible two-hop lineage back to
  these exact papers' own 3:4-LO family, flagged but NOT digit-for-digit
  verified this pass (concrete future task). (2) **Capability gap**: this
  project's own Jovian-tour CODE (`core/flyby.py`, `search/resonant_conic.py`,
  `search/ieg_seed.py`) is patched-conic/Lambert-based throughout, exactly the
  two-body approximation these papers' own text frames as the method their
  CRTBP dynamical-systems approach supersedes; the project's only existing
  homoclinic/heteroclinic-connection machinery (CCR4BP/CRNBP torus
  manifolds) targets a structurally different object class (quasi-periodic
  tori in time-periodic N-body models, not classical resonant periodic
  orbits in the autonomous 2-body-primary CR3BP). Mass-ratio cross-check:
  this project's own computed Jupiter-Europa `μ` (DE440-based) differs from
  the papers' stated value by ~0.034% relative — negligible, GM-vintage
  artifact, no correction needed.
- Citation-mining across both papers surfaced ~13 new candidates (none
  acquired); the strongest is **Anderson's own 2005 PhD dissertation** (the
  explicit root of the whole lineage, likely freely available at its cited
  university URL, not yet fetched), with **Howell-Marchand-Lo 2001** and
  **Johannesen & D'Amario 1999** each now independently re-corroborated a
  2nd/3rd time across this project's citation-mining history.
