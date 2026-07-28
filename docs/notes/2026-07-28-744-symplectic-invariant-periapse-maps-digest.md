# Digest: Symplectic-invariant lineage (Frauenfelder/Moreno/Aydin/Batkhin) + periapse-map lineage (Howell/Davis/Haapala) — 4 papers (#744)

**Task:** `#744`, continuing the `#730` consolidated acquisition backlog
(`docs/notes/2026-07-27-730-acquisition-backlog-master-list.md`). All four PDFs
filed to the private `cyclers_pdf` repo this session by the coordinating
session. **OCR status: all four text-layer, no OCR needed** (pdftotext
char-count/page verified far above the 10-char/page floor by the coordinating
session; not re-verified here).

Read `2026-07-27-728-moreno-aydin-vankoert-frauenfelder-koh-2024-bifurcation-graphs-digest.md`
(the already-acquired sibling paper by the same author group, same theory
lineage) and `2026-07-22-683-digest-davis-howell-2011-periapse-maps.md` (the
already-acquired periapse-map sibling) first, for house style and background —
both referenced throughout below rather than re-derived.

---

## 1. Frauenfelder, Koh & Moreno 2023 — "Symplectic Methods in the Numerical Search of Orbits in Real-Life Planetary Systems (With Numerical Explorations of the Jupiter-Europa and Saturn-Enceladus Systems)"

*SIAM J. Appl. Dyn. Syst.* 22(4):3284–3319 (2023). **DOI correction:** `#730`
§6 item 40 currently lists `10.1137/22M1506743`, which does not resolve. The
correct DOI is **`10.1137/22M1500459`** — CrossRef-confirmed by a prior
coordinating-session pass (`#743`), exact title/container/volume/issue/page
match. Use the corrected DOI going forward. arXiv:2206.00627. 30 pages.

### What it is

The **numerical-application companion** to paper 2 below (the GIT-quotient
theory paper) — confirmed directly by its own abstract ("The mathematical
framework is provided by the first and third authors in [10]", i.e. paper 2)
and by its own ref [16]/[17] (Koh, Anderson & Bermejo-Moreno's cell-mapping
method papers, JAS 68(1):172-196 (2021) and AAS 18-264 (2018), which supply
the actual numerical orbit-search engine). This is the **same paper already
summarized in detail** by the `#728` moreno-aydin-vankoert-frauenfelder-koh
digest's "Relevance to this codebase" section — that digest already extracted
its four tools (Floer numerical invariants, B-signs, GIT-sequence/Broucke
refinement, local CZ-index algorithm) accurately, since the 2024 bifurcation-
graphs paper builds directly on this one. This pass adds the numerics this
project's own digest didn't yet cover directly (the `#728` digest summarized
the *2024* paper's own numerics, not this one's).

### Sign-convention / Broucke-diagram claim (verified against this paper's own reference list, not the abstract paraphrase)

The abstract states the GIT-sequence framework "refined [the Broucke stability
diagram] with signs, reformulated as GIT quotients." **Important precision
missed by a looser reading:** this paper's own ref [2] for "Broucke's stability
diagram" is **Broucke, R., "Stability of periodic orbits in the elliptic,
restricted three-body problem," *AIAA J.* 7, 1003 (1969)** — a *different*
Broucke paper from the one repeatedly flagged elsewhere in this corpus
(Broucke, R., "Periodic Orbits in the Restricted Three-Body Problem with
Earth-Moon Masses," JPL TR 32-1168 (1968), `#730` §4 item 21, now flagged 4
independent times as a periodic-orbit-*database* gap). The 1969 AIAA paper is
a **stability-diagram/classification** paper, not the 1968 database — they
are not the same citation and should not be conflated when tallying "Broucke"
flags. Both are genuinely different Broucke papers, neither yet in corpus.

### Numerical worked examples (§6, positive-control candidates)

Two systems, both planar/spatial CR3BP:
- **Jupiter-Europa**: `µ = 2.5266448850435e-05`
- **Saturn-Enceladus**: `µ = 1.9002485658670e-07`

**"Snitch" configuration (Jupiter-Europa, H2-family doubly-symmetric
period-doubling, Fig. 4-9):**
- `γ_bef` (before bifurcation, type `E²`): Jacobi constant `c = 3.00357414`,
  period `T0 = 2.1215`; symmetric points `P1 = (1.016776, 0, 0, 0, 0.0130372, 0)`,
  `P2 = (0.997370, 0, 0, 0, -0.125493, 0)`. Planar eigenvalues
  `λp = -0.302203 ± i0.953244` (elliptic); spatial `λs = -0.999948 ± i0.010225`
  (elliptic, near the bifurcating unit circle).
- `γ_aft` (after bifurcation, type `EH⁻`): `c = 3.003571774`; symmetric points
  `P1 = (1.016787, 0, 0, 0, 0.013014, 0)`, `P2 = (0.997377, 0, 0, 0, -0.125701, 0)`.
  Planar eigenvalues stay elliptic (`λp = -0.309945 ± i0.950755`); spatial pair
  goes negative-hyperbolic: `λs = -0.972874`, `1/λs = -1.027883`.
- `β` (double-period orbit born at bifurcation, type `E²`, period
  `T1 = 4.245 ≈ 2·T0`): self-intersection points `P1..P4` given to ~1e-9
  precision (paper's own Remark 6.1: full MATLAB long-format precision exists,
  truncated to 6 digits for the printed text — treat printed digits as display,
  not full precision).
- Full 4×4/6×6 monodromy matrices are printed for every orbit (paper's own
  Eq. numbering, pp. 15-19) — usable as a direct positive-control target if
  this project ever builds Floquet-eigenvalue/monodromy-matrix machinery for
  Jupiter-Europa.

**Other Jupiter-Europa bifurcations (DRO/H2 family, Figs. 10-12):**
- Planar-to-planar period-tripling of the DRO family: `c = 2.9999, T0 = 2.504`
  (simple) → `T = 7.3 ≈ 3T0` (tripled).
- Planar-to-spatial 5-fold bifurcation, same DRO family: `c = 3.0005, T0 = 1.705`
  → `T = 8.52 ≈ 5T0`.
- Spatial-to-spatial period-doubling: `c = 3.0028, T0 = 4.62` → `T = 9.23 ≈ 2T0`.

**Saturn-Enceladus (Fig. 13, same planar-orbit family, two spatial
bifurcations off it):** `T1 = 1.2, T2 = 1.6, T3 = 2` (three simple orbits
along the family) → period-doubling gives `T ≈ 4.2 ≈ 2T3`; period-tripling
gives `T ≈ 4.85 ≈ 3T2`. No further numeric IC table beyond these period
values is given for Saturn-Enceladus in this paper (contrast with the `#728`
paper's Appendix A, which gives full sourced ICs/altitudes for
Saturn-Enceladus — this paper's Saturn-Enceladus content is illustrative only,
the detailed worked example lives in the already-acquired 2024 paper).

### Relevance to this codebase

Same verdict as the `#728` digest: **no CZ-index/B-sign/GIT-sequence
machinery exists anywhere in `cyclerfinder`** (re-confirmed, no new grep
needed — nothing has changed in this area since that digest). This paper adds
no new code-relevant surface beyond what `#728`'s digest already documented;
it is the numerics half of the theory paper 2 documents below, and is
included here primarily to close the `#730` §6 item 40 gap and fix its DOI.

### Citation-mining pass (24 references + 1 informally-numbered `[LTJ]`, full list read)

**Already in corpus:**
- None of this paper's 24 references are currently in `CORPUS_INDEX.md`
  (checked every author/keyword — Abraham & Marsden, Broucke ×1 [the 1969
  AIAA paper], Eliashberg-Givental-Hofer, Floer ×6, Frauenfelder-Moreno [10]
  = paper 2 below now being acquired, Frauenfelder & van Koert [11]/[12],
  Ginzburg, Hénon [14], Howard & MacKay, Koh-Anderson-Bermejo-Moreno [16]/[17],
  Li-Tao-Jiang `[LTJ]`, Krein ×4, Moser, Robbin & Salamon, Wonenburger — zero
  hits in `CORPUS_INDEX.md`).
- Moreno et al. 2024 bifurcation-graphs paper is already in corpus (not a
  reference of *this* paper, but its direct successor — already digested
  `#728`).

**Already flagged elsewhere in `#730` (recurring, no new action):**
- Frauenfelder & Moreno, GIT-quotients paper — this is paper 2, being
  acquired in this same task.
- Broucke (the 1969 AIAA stability paper, ref [2]) — **new** citation, see
  the precision note above; do not conflate with the already-4x-flagged 1968
  JPL TR database paper.

**Genuinely new candidates, topically overlapping (flagged, NOT acquired):**
- *High:* Koh, D., Anderson, R.L., Bermejo-Moreno, I., "Cell-mapping orbit
  search for mission design at ocean worlds using parallel computing," *J.
  Astronaut. Sci.* 68(1):172-196 (2021) [ref 16] — the actual numerical
  search engine both this paper and the already-acquired 2024 paper build on;
  independently corroborated by the 2024 paper's own citation of the same
  work (per the `#728` digest's own text). Koh, Anderson, Bermejo-Moreno,
  "Three-dimensional bifurcations in the circular restricted three-body
  problem," AAS/AIAA Astrodynamics Specialist Conference, AAS 18-264 (2018)
  [ref 17] — the conference precursor to the above.
- *Medium-high:* Frauenfelder, U. & van Koert, O., *The Restricted Three-Body
  Problem and Holomorphic Curves*, Pathways in Mathematics, Birkhäuser/Springer
  (2018) [ref 11] — the standard reference textbook this entire symplectic
  toolkit sits inside; independently re-flagged by paper 2 below (2-flag
  recurrence within this single task). Frauenfelder, U. & van Koert, O., "The
  Hörmander index of symmetric periodic orbits," *Geom. Dedicata* 168:197-205
  (2014) [ref 12] — same authors, directly relevant symmetric-orbit index
  theory; also independently re-flagged by paper 2.
- *Medium:* Ginzburg, V.L., "The Conley conjecture," *Ann. of Math.* 172(2):
  1127-1180 (2010) [ref 13] — the foundational reference for the local
  Floer-homology Euler-characteristic invariant used throughout; independently
  re-flagged by paper 3 below (2-flag recurrence). Howard, J.E. & MacKay,
  R.S., "Linear stability of symplectic maps," *J. Math. Phys.* 28:1038-1051
  (1987) [ref 15] — the `n=3,4` generalization of the Krein-signature/Broucke-
  diagram machinery, explicitly cited as the higher-dimensional analogue;
  independently re-flagged by paper 2 below (2-flag recurrence).
- *Low-medium:* Li, Q., Tao, Y., Jiang, F., "Orbital Stability and Invariant
  Manifolds on Distant Retrograde Orbits around Ganymede and Nearby
  Higher-Period Orbits," *Aerospace* 9(8):454 (2022) [ref LTJ] — DRO-family
  stability, tangentially adjacent (Ganymede not Europa/Enceladus, no CZ-index
  content). Robbin, J. & Salamon, D., "The Maslov index for paths," *Topology*
  32:827-844 (1993) [ref 23] — Maslov-index foundations underlying the CZ-index
  construction. Wonenburger, M., "Transformations which are products of two
  involutions," *J. Math. Mech.* 16:327-338 (1966) [ref 24, printed as 1996 in
  this paper's own reference list but almost certainly a typo for 1966 —
  cross-check against paper 2's ref [18], which also cites this paper and also
  prints "1996"; likely a shared, uncorrected error in both papers rather than
  two independent typos] — the algebraic fact (`MA,B,C` normal form) underlying
  the whole GIT-quotient construction; independently re-flagged by paper 2
  below (2-flag recurrence, same citation, same possible year typo in both).
- Not flagged (foundational textbook/background, out of this project's search
  method domain per prior digests' triage convention): Abraham & Marsden
  1978 (also flagged by paper 2, see below — but purely general-mechanics
  background, low priority even on 2nd flag); the six Floer 1988-89 papers
  (pure Floer-homology foundations, not orbit-search-method specific);
  Eliashberg-Givental-Hofer 2000 SFT foundations; the four Krein 1950-51
  papers and Moser 1958 (classical stability-theory foundations, superseded
  for this project's purposes by the already-summarized B-sign generalization
  itself); Hénon 1969 [ref 14, note: **distinct from** the Hénon papers
  flagged by paper 3 below — this is "Numerical Exploration of the Restricted
  Three-Body Problem V" *A&A* 1:223-238 (1969), already effectively covered
  by paper 3's citation-mining list, not re-flagged separately here to avoid
  double-counting].

---

## 2. Frauenfelder & Moreno 2023 — "On GIT Quotients of the Symplectic Group, Stability and Bifurcations of Periodic Orbits (With a View Towards Practical Applications)"

*J. Symplectic Geom.* 21(4):723–773 (2023). DOI **`10.4310/jsg.2023.v21.n4.a3`**
(CrossRef-confirmed: title/container/volume/issue/page match). arXiv:2109.09147.
33 pages.

**Title change confirmed.** The PDF's own title page (arXiv v2 header, line 2)
reads: **"On GIT Quotients of the Symplectic Group, Stability and Bifurcations
of Symmetric Orbits"** — no "practical applications" clause, and "Symmetric
Orbits" rather than "Periodic Orbits." The published title is independently
confirmed from a fully separate source: paper 3's own reference list (Aydin &
Batkhin 2025) cites this exact paper as **"On GIT quotients of the symplectic
group, stability and bifurcations of periodic orbits (with a view towards
practical applications). *J. Symplectic Geom.* 21(4), 723–773 (2023)"** —
matching the task's stated published title verbatim, byte-for-byte on the
"periodic orbits (with a view towards practical applications)" wording. Same
paper, retitled between arXiv preprint and journal publication; content
(abstract, sections, references) is otherwise unchanged between the versions
examined.

### The actual GIT-quotient/B-sign theory (operational content)

This is the **pure-math groundwork** paper 1 above and the already-acquired
2024 bifurcation-graphs paper both build on — confirmed directly, both papers
cite this one as their theoretical foundation (paper 1's ref [10]; the 2024
paper's ref [6], already flagged in the `#728` digest as "the actual
mathematical groundwork this whole paper builds on").

**Setup.** A symmetric periodic orbit's *reduced monodromy matrix* at a
symmetric point takes a constrained symplectic form `M_{A,B,C}` (Eq. 1, with
`B, C` symmetric and `A²-BC=I`) — a classical fact due to Wonenburger (any
`2n×2n` symplectic matrix is symplectically conjugate to one of this form).
The paper studies the **GIT quotient** `Sp(2n)//Sp(2n)` (identifying matrices
whose conjugation orbits' closures intersect — fixes the non-Hausdorff-ness
of the naive quotient `Sp(2n)/Sp(2n)`) and a finer quotient
`SpI(2n)//GLn(R)` that additionally tracks a **B-signature**
(`±` sign per elliptic/hyperbolic eigenvalue, generalizing Krein's classical
sign — which only applied to the elliptic case — to also cover hyperbolic
eigenvalues, via Definition 3.3).

**The operational content (what this actually lets you conclude, given two
periodic orbits' Floquet multipliers):**
1. Compute each orbit's reduced monodromy `M_{A,B,C}`, extract `(tr A, det A)`
   — this maps the orbit to a point `p ∈ Rⁿ` (`n=1` planar CR3BP, `n=2`
   spatial CR3BP) in the base of the **GIT sequence**
   `SpI(2n)//GLn(R) → Sp(2n)//Sp(2n) → Mn×n(R)//GLn(R) ≅ Rⁿ`.
2. Remove the **bifurcation locus** (matrices with `±1` as an eigenvalue) from
   each quotient space. What remains decomposes into finitely many **connected
   components** — for `n=2` (the spatial CR3BP case), the paper *explicitly
   enumerates* these: **8 connected components** in the coarser
   `Sp(4)//Sp(4)` (Krein-sign-only information), and **19 connected
   components** in the finer `SpI(4)//GL2(R)` (B-signature information,
   distinguishing more matrices — the B-sign refinement genuinely adds
   discriminating power, not just bookkeeping).
3. **The load-bearing theorem:** if two symmetric orbits' matrices land in
   *different* connected components of either quotient space (equivalently:
   different B-sign labels, or different regions of the Broucke-type diagram
   in the base `Rⁿ`), **no continuous path between them avoiding the
   bifurcation locus exists** — i.e. **no regular orbit cylinder can connect
   them without a bifurcation occurring somewhere in between.** This is a
   genuine topological obstruction, not a heuristic: it follows because the
   maps in the GIT sequence send connected components to connected components
   (acting as branched covering maps away from the branch locus), so two
   points in different components can never be joined by a lift of a path in
   the base that itself avoids the bifurcation locus.
4. Concretely for engineers (per the paper's own framing): given a bifurcation
   observed numerically, checking which connected component the pre- and
   post-bifurcation orbits' matrices land in — and whether their B-signs
   differ — is a **necessary condition for connectivity** that can rule out
   candidate continuations *before* attempting a numerical continuation, and
   conversely, orbits with identical component+B-sign are a *candidate* (not
   guaranteed) pair for a regular connecting family.

The "upcoming work" mentioned in the abstract (Floer-theoretic Euler
characteristics, the SFT-Euler-characteristic bifurcation invariant) is
exactly paper 1's own subject matter — confirmed cross-reference, the two
papers are a deliberate theory/numerics pair from the same author group.

### Relevance to this codebase

Same "no such machinery exists" verdict as papers already digested under
`#728`. No new code-relevant surface — this paper is the mathematical
substrate paper 1 and the already-acquired 2024 paper both implement; it adds
precision to *why* the B-sign/GIT-sequence approach works (the connected-
components/covering-map argument) but does not itself introduce anything
numerically new beyond what those application papers already surface.

### Citation-mining pass (19 references, full list read)

**Already in corpus:** none (checked every reference — zero hits in
`CORPUS_INDEX.md`).

**Already flagged elsewhere in `#730` / within this task (recurring):**
- Frauenfelder & van Koert 2018 book [ref 6] and 2014 Hörmander-index paper
  [ref 7] — both independently re-flagged from paper 1 above (2nd flag each).
- Broucke [ref 3] — **same 1969 AIAA stability paper** as paper 1's ref [2]
  (identical citation, confirmed by matching journal/volume/page: "R. Broucke,
  Stability of periodic orbits in the elliptic, restricted three-body problem.
  AIAA J. 7,1003 (1969)"). 2nd flag within this task; still distinct from the
  1968 JPL TR database paper flagged elsewhere.
- Howard & MacKay 1987 [ref 9] — same paper as paper 1's ref [15], 2nd flag.
- Wonenburger 1966(?) [ref 18] — same paper as paper 1's ref [24], 2nd flag,
  same possible year-typo ("1996" printed in both papers' reference lists).
- Krein [refs 11-14] and Moser [ref 16] — same four Krein papers + Moser 1958
  as paper 1's refs [18-22]; not re-flagged individually (classical-theory
  background, already triaged low priority above).

**Genuinely new candidates, topically overlapping (flagged, NOT acquired):**
- *Medium-high:* Albers, P. & Frauenfelder, U., "The space of linear
  anti-symplectic involutions is a homogeneous space," *Arch. Math. (Basel)*
  99(6) (2012) [ref 2] — directly foundational to the antisymplectic-
  involution/symmetric-orbit setup this whole theory relies on, same
  Frauenfelder-lineage authorship.
- *Medium:* Bruno, A.D., *The Restricted 3-Body Problem: Plane Periodic
  Orbits*, De Gruyter (1994; Russian orig. 1990, Nauka) [ref 4] —
  **independently re-flagged by paper 3 below** (2-flag recurrence within
  this task; same book, confirmed by matching title/publisher/translator
  credit). Hénon, M., *Generating Families in the Restricted Three-Body
  Problem. II. Quantitative Study of Bifurcations*, Lecture Notes in Physics
  Monographs 65, Springer (2001) [ref 8] — **note: this is Volume II**, a
  distinct, later monograph from the already-flagged Hénon 1997 *Generating
  Families in the Restricted Three-Body Problem* (Volume I, `#730` §4 item 29,
  flagged by `#725`); do not treat as a duplicate of item 29, it is the
  quantitative-bifurcations sequel volume. Long, Y., *Index Theory for
  Symplectic Paths with Applications*, Progress in Mathematics 207, Birkhäuser
  (2002) [ref 15] — the standard Maslov/CZ-index-theory reference monograph.
- *Low-medium:* Ekeland, I., *Convexity Methods in Hamiltonian Mechanics*,
  Ergebnisse der Mathematik (3) 19, Springer (1990) [ref 5] — general
  Hamiltonian-convexity background, not orbit-search-method specific.
  Kalantonis, V., "Numerical Investigation for Periodic Orbits in the Hill
  Three-Body Problem," *Universe* 6(6):72 (2020) [ref 10] — **independently
  re-flagged by paper 3 below** (2-flag recurrence, same paper). Zhou, B.,
  "Iteration Formulae for Brake Orbit and Index Inequalities for Real
  Pseudoholomorphic Curves," arXiv:2011.07958 (2020) [ref 19] — index-theory
  preprint, tangential (brake orbits, not this project's periodic-orbit
  classes). Abraham, R. & Marsden, J., *Foundations of Mechanics*, 2nd ed.,
  Addison-Wesley (1978) [not in this paper's own list but recurring — see
  paper 3's identical citation] — general mechanics textbook, low priority.
- Szebehely 1967 [ref 17] — already in corpus, digested (`#394`/`#400`).

---

## 3. Aydin & Batkhin 2025 — "Studying Network of Symmetric Periodic Orbit Families of the Hill Problem via Symplectic Invariants"

*CMDA* 137(2):12 (2025). **DOI correction:** `#730` §6 item 42 currently lists
`10.1007/s10569-025-10233-0`, which does not resolve. The correct DOI is
**`10.1007/s10569-025-10241-7`** — confirmed directly from this PDF's own
header (printed on p.1: "Celestial Mechanics and Dynamical Astronomy (2025)
137:12, https://doi.org/10.1007/s10569-025-10241-7") and independently via
CrossRef. Use the corrected DOI going forward. 77 pages, open access.

### What it is

Applies the CZ-index/GIT-quotient/B-sign toolkit of papers 1-2 above (same
theory, different application) to the **spatial Hill three-body problem** —
the `µ→0` scaling limit of the CR3BP, obtained (per the paper's own §2.1,
Eq. following (2)) by shifting the smaller primary to the origin, rescaling
coordinates by `µ^(1/3)`, and Taylor-expanding the CR3BP Hamiltonian
`H(x,y,z,px,py,pz) = ½(px²+py²+pz²) - (1-µ)/r1 - µ/r2 + pxy - pyx` (identical
convention to the SCR3BP Hamiltonian in papers 1's numerics section, with
Jacobi constant `Γ = -2H`) to leading order in `µ^(1/3)`, discarding
`O(µ^(1/3))` remainder terms. This is explicitly framed as **directly
building on the first author's own prior work**: Aydin (2023b), "The
Conley-Zehnder indices of the spatial Hill three-body problem," *CMDA* (2023)
— the paper this one's own network-construction extends — and cites the
already-acquired Moreno-Aydin-van Koert-Frauenfelder-Koh 2024 bifurcation-
graphs paper directly (its own reference list: "Moreno, A., Aydin, C., van
Koert, O., Frauenfelder, U., Koh, D.: Bifurcation graphs for the CR3BP via
symplectic methods. *J. Astronaut. Sci.* 71(6), 51 (2024)").

**Family network structure (abstract + Table 1/2):** the "basic" natural
families are the direct/prograde satellite family `g` (with its
symmetry-breaking pitchfork daughter `g'`), the retrograde family `f` (and
its own daughter `f3`), the planar and vertical Lyapunov (libration) families
`a`, `c`, and the rectilinear-vertical-consecutive-collision family `B0±`.
Table 1 tabulates, for each family, *how many n-th multiple covers* bifurcate
and interconnect with which other family's own n-th cover (e.g. `g`'s 2nd
cover connects to `g'`'s 2nd cover and `B0±`'s 3rd cover; `g`'s 3rd cover
connects to `f`'s 5th cover and `f3`'s 1st cover; etc. — a genuine "bridge"
network, not a simple linear chain). Table 2 splits these connections into
**"explicit" interactions** (a genuine bridge family of orbits directly
connects the two, requiring matching symmetry types across the relevant
covers) versus **"implicit" interactions** (two branch families merely *share*
periodic orbits via a symmetry-breaking pitchfork, without a direct bridging
family) — a distinction the paper states is not observable by direct
continuation/stability-index computation alone, and requires the CZ-index
machinery to detect (per the paper's own §6 conclusion, point 4).

**CZ-index bifurcation detection:** identical mechanism to papers 1-2 — the
CZ-index grades local Floer homology; its Euler characteristic is a
bifurcation invariant (index jumps exactly at bifurcations, stays constant
along a smooth family); tracking index jumps across the whole family network
is what lets the paper construct and verify its bifurcation graphs (§3-4).

### Numeric worked examples / positive-control candidate data

**Appendix "Tables of data"** (Tables 3-19, ~30 pages) give, per family
member: initial condition `(x(0), ẏ(0))` (planar) or the spatial analogue,
period `T`, Jacobi-like constant `Γ`/`C` (Hill-scaled units), the `(C/B)`-sign
pair with Floquet-multiplier data (either a rotation angle `φ` for elliptic
pairs on the unit circle, noting `k`-th-root-of-unity resonances explicitly,
or a real hyperbolic multiplier `λ`), and the triple of CZ-indices
`(µ_CZ^p, µ_CZ^s, µ_CZ)` — planar, spatial, and combined.

Sample rows, Table 3 (family `g`, first block — all in Hill-scaled units,
`µ→0` limit, no explicit mass ratio needed since Hill units are
mass-ratio-independent by construction):
- `x(0)=6.50888000, ẏ(0)=0.17610000, T=2.22291184, C=0.50799`, elliptic
  planar/spatial pair `(φp=0.449, φs=0.535)`, CZ-indices `6/3/3`.
- `x(0)=4.27892454, ẏ(0)=0.30115821, T=1.62301941, C=1.41824`, planar
  **hyperbolic** `λp=2.265` (negative, real), spatial elliptic
  `φs=1.570 (2π/4-resonant)`, CZ-indices `5/2/3`.
- `x(0)=-0.1080506, ẏ(0)=0.05018349, T=6.32213172, C=8.68858`, both planar and
  spatial hyperbolic (`λp=6358, λs=19.82`), CZ-indices `6/2/4`.

These are directly usable as a positive-control table for any Hill-problem or
`µ→0`-limit CR3BP family-tracing work this project might undertake (a
mass-ratio-independent census, unlike the Jupiter-Europa/Saturn-Enceladus
tables in paper 1, which are tied to specific `µ` values).

### Cross-validation against this project's own code — clean negative, not executed

**No dedicated Hill-problem (`µ→0`-scaled) module exists anywhere in
`cyclerfinder`** (checked: `grep -rli "hill.*problem\|hill_problem\|hill three.body\|hill's equations"` across `src/cyclerfinder/` returns only files using "Hill sphere" — the gravitational-dominance-radius concept — never the Hill-problem *limiting equations of motion*). This project's CR3BP code
(`core/cr3bp.py` and friends) uses full-`µ` equations throughout; there is no
`µ^(1/3)`-rescaled Hill-limit variant to directly plug this paper's tables
into. A cross-check against this paper's Table 3-19 data is therefore a
genuine future opportunity (as the task framing anticipated) but is **not
executed in this pass** — it would require either implementing the Hill-limit
equations of motion (a real, if modest, new-code task) or taking the `µ→0`
limit analytically inside the existing CR3BP corrector at some very small test
`µ` and checking convergence toward this paper's tabulated Hill-scaled values,
neither of which is in scope here.

### Citation-mining pass (58 references, full intro/background + full reference list read)

**Already in corpus:**
- Restrepo, R.L. & Russell, R.P., "A database of planar axi-symmetric periodic
  orbits for the solar system," *CMDA* (2018) — this paper's own ref
  matches, already in corpus (`restrepo-russell-2018-...`, digested
  `2026-06-17-digest-restrepo-russell-2018.md`); this paper cites it directly
  for the `µ` convention used (p.5, Sun-Earth `µ=3.0404317e-6` example) —
  **resolves the `#730` §1 "likely duplicate" flag on the 2017 AAS-17-694
  precursor exactly the same way the `#742` digest already resolved it for a
  different citing paper**, corroborating evidence, no new action.
- Moreno, Aydin, van Koert, Frauenfelder, Koh 2024 — already in corpus,
  digested `#728`.
- Szebehely 1967 — already in corpus, digested `#394`/`#400`.

**Already flagged elsewhere in `#730` / within this task (recurring):**
- Frauenfelder & Moreno 2023 GIT-quotients paper — being acquired as paper 2
  in this same task (this paper's own reference gives the exact published
  title/venue used to confirm the title-change finding above).
- Frauenfelder & van Koert 2018 book — 3rd flag within this task (papers 1
  and 2 above, now paper 3).
- Ginzburg 2010 — 2nd flag within this task (paper 1 above).
- Bruno, A.D., *The Restricted 3-Body Problem: Plane Periodic Orbits*, De
  Gruyter (1994) — 2nd flag within this task (paper 2 above), confirmed same
  book (Russian original "Nauka, Moscow, 1990" credit matches).
- Kalantonis 2020 — 2nd flag within this task (paper 2 above).
- Hénon, M., *Generating Families in the Restricted Three-Body Problem*,
  Springer (1997) — `#730` §4 item 29, already flagged by `#725`; this
  paper's own reference list confirms it as a *distinct* citation from the
  Hénon 2001 Vol. II sequel flagged as new in paper 2 above.

**Genuinely new candidates, topically overlapping (flagged, NOT acquired):**
- *High:* Aydin, C., "The Conley-Zehnder indices of the spatial Hill
  three-body problem," *CMDA* (2023b), DOI `10.1007/s10569-023-10134-7` —
  this paper's own explicit predecessor/direct foundation (the network-graph
  construction this paper extends); the single highest-priority gap from this
  paper's own citation list. Villac, B.F. & Scheeres, D.J., "Escaping
  trajectories in the Hill three-body problem and applications," *JGCD*
  (2003), DOI `10.2514/2.5062` — **the periapsis-Poincaré-map originating
  paper**, per the already-acquired `#683` Davis-Howell 2011 digest's own
  citation of it; independently surfaces here as a Hill-problem escape-
  trajectory reference, and is directly relevant to paper 4 below (flagged
  jointly, see paper 4's citation-mining section).
- *Medium-high:* Aydin, C., "From Babylonian lunar observations to Floquet
  multipliers and Conley-Zehnder indices," *J. Math. Phys.* (2023a), DOI
  `10.1063/5.0156959` — same-author CZ-index application to real lunar orbit
  data, methodologically adjacent. Aydin, C., "Contact geometry of Hill's
  approximation in a spatial restricted four-body problem," *Arnold Math. J.*
  (2025, to appear) — **the four-body Hill extension**, directly relevant if
  this project ever wants a four-body (CCR4BP-adjacent) Hill-limit model;
  not yet published as of this paper's own submission, may not be acquirable
  yet. Gómez, G., Marcote, M. & Mondelo, J.M., "The invariant manifold
  structure of the spatial Hill's problem," *Dyn. Syst.* 20(1):115-147
  (2005) — **independently re-flagged**, already surfaced in the `#742`
  digest's own citation-mining list for the Barrabés-Mondelo-Ollé 2009 paper
  (2-flag recurrence across separate tasks). Belbruno, E., Frauenfelder, U. &
  van Koert, O., "A family of periodic orbits in the three-dimensional lunar
  problem," *CMDA* (2019), DOI `10.1007/s10569-019-9882-8` — same
  Frauenfelder-lineage authorship, directly on-point (3D lunar-problem
  periodic-orbit families). Tsirogiannis, G.A., Perdios, E.A. & Markellos,
  V.V., "Improved grid search method: an efficient tool for global
  computation of periodic orbits. Application to Hill's problem," *CMDA*
  103:49-78 (2009) — **independently re-flagged**, already surfaced in the
  `#742` digest's citation-mining list for Franz & Russell 2022 (2-flag
  recurrence across separate tasks, and directly Hill-problem-specific here).
- *Medium:* Batkhin, A.B. (4 papers: "Symmetric periodic solutions of the
  Hill's problem I," *Cosm. Res.* 51(4):275-288 (2013); "Web of families of
  periodic orbits of the generalized Hill problem," *Dokl. Math.* 90(2):
  539-544 (2014); "Bifurcations of periodic solutions of a Hamiltonian system
  with a discrete symmetry group," *Program. Comput. Softw.* 46(2):84-97
  (2020); with Batkhina, N.V., "Hierarchy of periodic solutions families of
  spatial Hill's problem," *Sol. Syst. Res.* 43(2):178-183 (2009)) — the
  second author's own prior Hill-problem family-classification work, direct
  groundwork for this paper's own network construction. Lidov, M.L. (5
  papers/reports, 1962-1988, mostly Russian-language *Kosmicheskie
  Issledovaniia*/*Pis'ma v Astronomicheskii Zhurnal* reports) — classical
  spatial-Hill-family construction methodology (KS-transformation-based),
  foundational but likely low-accessibility (non-English conference/report
  literature). Hénon, M. (4 further papers beyond the already-flagged 1997
  volume: 1969 *A&A* 1:223-238 "V. Hill's Case: Periodic Orbits", 1970 *A&A*
  9:24-36 "VI. Hill's Case: Non-Periodic Orbits", 1974 *A&A* 30:317-321
  "Vertical stability... II. Hill's Case", 2003 *CMDA* 85:223-246 "New
  families of periodic orbits in Hill's problem", 2005 *CMDA* 93:87-100
  "Families of asymmetric periodic orbits in Hill's problem") — the classical
  Hill-problem numerical-exploration series this whole family taxonomy is
  built on; high volume, flagged as a cluster rather than individually
  prioritized. Perko, L., "Families of symmetric periodic solutions of Hill's
  problem I/II," *Am. J. Math.* 104(2):321-397 (1981a/b) — classical
  first/second-species Hill-problem family classification.
- *Low-medium:* Michalodimitrakis, M., "Hill's problem: families of
  three-dimensional periodic orbits (part I)," *Astrophys. Space Sci.*
  68:253-268 (1980) — halo-orbit-family precursor (identifies the family that
  becomes `a1v`/halo, per this paper's own intro). Zagouras, C. & Markellos,
  V.V., "Three-dimensional periodic solutions around equilibrium points in
  Hill's problem," *Celest. Mech.* 35:257-267 (1985). Conley, C. & Zehnder,
  E., "Morse-type index theory for flows and periodic solutions for
  Hamiltonian equations," *Comm. Pure Appl. Math.* 37(2):207-253 (1984) — the
  origin of the Conley-Zehnder index itself, foundational but likely
  superseded for this project's purposes by the already-summarized
  applications. Salamon, D. & Zehnder, E., "Morse theory for periodic
  solutions of Hamiltonian systems and the Maslov index," *Comm. Pure Appl.
  Math.* 45:1303-1360 (1992). Bruno-lineage/Russian-literature items (Batkhin,
  Lidov above) may be difficult to acquire — flagged for visibility, not
  prioritized for acquisition effort.
- Not flagged (general background, out of domain per prior triage
  convention): Abraham & Marsden 1978 (also cited in paper 2, low priority
  confirmed again); Arnold & Avez 1968 ergodic-theory textbook; Meyer & Offin
  2017 (3rd ed. of the standard Hamiltonian-dynamics textbook — note this is
  a *different* edition/authorship pairing from the already-flagged Meyer &
  Hall 1992 in the `#728` moreno-bifurcation-graphs digest's own citation
  list; likely not worth separately acquiring given the 1992 edition is
  already flagged); Poincaré 1893, Hill 1878 (primary historical sources,
  not separately acquirable in useful form); Marchal 1990 general 3BP
  textbook; Lamb & Roberts 1998 time-reversal-symmetry survey (general
  dynamical-systems background); Kreisman 2005, Llibre-Meyer-Soler 1999,
  Meletlidou-Ichtiaroglou-Winterberg 2001, Morales-Ruiz-Simó-Simon 2005,
  Simó-Stuchi 2000, Waldvogel 2008, Poleshchikov 2003, Meyer-Schmidt 1982,
  Robin-Markellos 1980, Hofer-Wysocki-Zehnder 1998/2003, Gurfil-Kasdin-Kolemen
  2005, Hamilton-Krivov 1997, Scheeres 1998, Burgos-García-Gidea 2015, Lara-
  Palacián 2009, Clohessy-Wiltshire 1960, Darwin 1911 — background/tangential
  to varying degrees, consistent with prior digests' treatment of large
  reference-list tails; not individually triaged given the already-substantial
  new-candidate list above.

---

## 4. Howell, Davis & Haapala 2012 — "Application of Periapse Maps for the Design of Trajectories Near the Smaller Primary in Multi-Body Regimes"

*Math. Probl. Eng.* 2012:351759. DOI `10.1155/2012/351759` (CrossRef-confirmed
and directly visible in the PDF's own Hindawi header). 23 pages, open access.
**Note on this section's provenance:** the sub-agent originally dispatched
internally for this paper did not return in time; this section was completed
directly by the coordinating session (`#744`) against the corpus PDF via
`pdftotext`, not inferred from any secondary source.

### What it is / relation to the already-acquired `#683` Davis & Howell 2011 digest

**Same Davis/Howell/Haapala periapse-map lineage, but the review/synthesis
article, not a new method paper.** Confirmed directly: this paper's own
Introduction credits the periapse Poincaré map's origin to **Villac &
Scheeres 2003** ("Escaping trajectories in the Hill three-body problem and
applications," *JGCD* 26(2):224-232 — ref [28]), extended by **Paskowitz &
Scheeres 2006** (×2, Europa "safe zones," refs [29]/[30]) to the Hill problem,
then by **Davis & Howell** (refs [31]/[32], the AAS/IAC conference precursors
to the already-acquired `#683` Acta Astronautica 69:1038 journal paper) and
**Haapala & Howell** (ref [33]) to the full CR3BP short- and long-term
behavior. This 2012 paper packages the same Sun-Saturn periapse-map machinery
`#683` already documents (escape/impact/capture lobes, `Πt` region notation,
manifold-tube-to-lobe correspondence) plus TWO worked application examples
`#683` does not cover: an Earth-Moon low-energy ballistic-transfer design and
an L1-gateway arrival-trajectory design — both using the *same* underlying
map machinery, not a new theoretical contribution over `#683`.

**Construction (Section 3, consistent with `#683`):** each point in a periapse
map is an initial condition at a fixed Jacobi constant reflecting a periapsis
(`ṙ=0`, `r̈≥0`); propagated forward and classified by outcome after each
subsequent periapsis (impact P2, escape via L1 gateway, escape via L2 gateway,
or continued capture) — colored accordingly. Escape is defined as `x` more
than 0.01 nondimensional units beyond L1/L2. "Long-term" maps (Sun-Saturn,
`#683`'s own regime) versus the position-space periapse-map representation
itself (`Δrp`-colored maps, arrival-lobe contours `Γ^{U/S}_{L1/L2,n}`) are the
two representational variants demonstrated.

### Numeric worked examples (positive-control candidates)

**Example 1 — Earth-Moon ballistic lunar transfer (§4.2):** From a 167 km
circular Earth parking orbit (`C = 3.068621`), a maneuver of **Δv = 3.2 km/s**
lowers Jacobi constant to `C = 3.000785` in the Sun-Earth system, opening the
ZVCs enough for solar gravity to raise periapsis by `Δrp = 377,855 km =
0.2525 r_H` (Sun-Earth Hill radius) — reaching the Moon's orbital radius
(384,400 km) from `R_E + 167` km. The paper's own text states this **Δv
threshold (3.2 km/s, exact cutoff 3.199 km/s) agrees well with**: (a)
Sweetser's theoretical minimum of **3.099 km/s** for a 167 km parking orbit
[ref 39]; (b) Parker & Born's optimized **≈3.2 km/s** for a 185 km parking
orbit in a full Sun-Earth-Moon model [ref 38] — an explicit, citable
three-way cross-check of the same design number, directly reusable as a
positive-control triple for any ballistic-lunar-transfer ΔV-floor check this
project might run.

**Example 2 — L1-gateway arrival to a lunar Lyapunov orbit (§4.3, Earth-Moon
system):** at `C = 3.17212`, a periapsis selected on the arrival lobe
`Γ^U_{L1,1}` (on a stable manifold of an L1 Lyapunov orbit) gives: Earth
departure `Δv1 = 3.105 km/s` from a 200 km circular parking orbit, manifold
insertion `Δv2 = 0.630 km/s`, close lunar flyby at 100 km altitude, asymptotic
approach to the L1 Lyapunov orbit — **total Δv = 3.735 km/s**. An alternative
capture-into-100-km-lunar-orbit variant of the same trajectory needs an
additional periapsis maneuver of **0.631 km/s**. Both fully sourced, directly
tabulated numbers (not digitized off a plot).

### Cross-check: does this paper's Sun-Saturn periapse-map machinery match `#683`'s own characterization?

**Consistent — no discrepancy.** The lobe/gateway/`Πt`-region terminology,
the escape-classification threshold (0.01 nondimensional units beyond
L1/L2), and the manifold-tube-to-escape-lobe correspondence in this paper's
Section 3.2 (Sun-Saturn example, Fig. 5) match the already-acquired `#683`
digest's own description of the Davis-Howell method exactly — expected, since
both papers are the same author group describing the same underlying
technique; this 2012 paper is best read as a broader-scope review/tutorial
companion to `#683`'s own more narrowly-scoped journal treatment, not an
independent re-derivation.

### Mandatory citation-mining pass (39 references, full list read)

**Already in corpus:** ref [11] Koon, Lo, Marsden & Ross, "Heteroclinic
connections between periodic orbits and resonance transitions in celestial
mechanics," *Chaos* 10(2):427-469 (2000) — in corpus, digested/mined
(`#314`). Ref [17] Gómez, Koon, Lo, Marsden, Masdemont & Ross, "Connecting
orbits and invariant manifolds in the spatial restricted three-body problem,"
*Nonlinearity* 17(5):1571-1606 (2004) — in corpus, digested (`#499`). Ref
[27] Anderson & Lo, "A Dynamical Systems Analysis of Resonant Flybys:
Ballistic Case" — this is `#730` §3 item 13, whose DOI this same coordinating
session independently CrossRef-confirmed earlier in this task
(`10.1007/BF03321164`) — direct cross-task corroboration, not a new flag.

**Already flagged elsewhere in `#730` (recurring, no new action):** ref [34]
Villac & Scheeres, "On the Concept of Periapsis in Hill's Problem," *CMDA*
90(1-2):165-178 (2004) — `#730` §9 low-priority tail, already listed there.
Ref [25] Tsirogiannis, Perdios & Markellos 2009 — already independently
flagged twice within this same digest's own citation-mining (papers 3 above
and the already-completed `#742` digest) — 3rd flag overall, raise to
**HIGH**. Ref [28] Villac & Scheeres 2003 (Hill-problem escaping
trajectories, periapsis-Poincaré-map ORIGIN paper) — 2nd flag within this
digest (already surfaced by paper 3 above); confirmed here as the single
paper both the Aydin-Batkhin Hill-network paper AND this periapse-map review
cite as foundational — raise to **HIGH**, now a 2-paper-within-task
recurrence plus direct `#683` relevance (the technique `#683`/this paper both
build on).

**Genuinely new candidates, topically overlapping (flagged, NOT acquired):**
- *High:* Davis, D.C., *Multi-Body Trajectory Design Strategies Based on
  Periapsis Poincaré Maps*, PhD dissertation, Purdue University (2011) [ref
  35] — the direct source thesis underlying both this paper and the
  already-acquired `#683` journal paper; Haapala, A.F., *Trajectory Design
  Using Periapsis Poincaré Maps and Invariant Manifolds*, MS thesis, Purdue
  University (2011) [ref 36] — companion source thesis. Both are the
  foundational full-detail treatments this 2012 review paper and `#683`'s
  journal paper both compress from.
- *Medium-high:* Paskowitz, M.E. & Scheeres, D.J., "Robust Capture and
  Transfer Trajectories for Planetary Satellite Orbiters," *JGCD* 29(2):
  342-353 (2006) [ref 29], and "Design of Science Orbits About Planetary
  Satellites: Application to Europa," *JGCD* 29(5):1147-1158 (2006) [ref 30]
  — the Hill-problem "safe zone" extension papers this paper's own
  Introduction credits as the direct bridge between Villac-Scheeres 2003 and
  the CR3BP periapse-map lineage; directly Europa-relevant (this project's
  own Jovian search domain). Sweetser, T.H., "An Estimate of the Global
  Minimum ΔV Needed for Earth-Moon Transfer," AAS 91-101 (1991) [ref 39] —
  the theoretical-minimum-ΔV cross-check cited in Example 1 above; **note:
  distinct from** the already-flagged `#730` §9 "Sweetser et al. 1997
  (historical Europa-orbiter mission-design study)" — same author, different
  paper/year/topic, do not conflate.
- *Medium:* Belbruno, E. & Miller, J.K., "Sun-Perturbed Earth-to-Moon
  Transfers with Ballistic Capture," *JGCD* 16(4):770-775 (1993) [ref 10] —
  classical ballistic-capture design paper, not currently in corpus (checked;
  the corpus's existing Belbruno coverage is the 2004 textbook and other
  items already noted in `#730` §9, not this specific 1993 paper). Belbruno,
  E., Gidea, M. & Topputo, F., "Weak Stability Boundary and Invariant
  Manifolds," *SIAM J. Appl. Dyn. Syst.* 9(3):1061-1089 (2010) [ref 26] —
  WSB/manifold-connection theory, not currently in corpus. Russell, R.P.,
  "Global Search for Planar and Three-Dimensional Periodic Orbits Near
  Europa," *J. Astronaut. Sci.* 54(2):199-226 (2006) [ref 20] — directly
  Europa-relevant systematic periodic-orbit search, same domain as this
  project's own Jovian work. Anderson, R.L. & Lo, M.W., "Virtual Exploration
  by Computing Global Families of Trajectories with Supercomputers," AAS
  05-220 (2005) [ref 19] — same Anderson/Lo lineage as `#730` §3 items 12-20,
  a new (not previously listed) member of that cluster. Howell, K.C.,
  Marchand, B.G. & Lo, M.W., "Temporary Satellite Capture of Short-Period
  Jupiter Family Comets from the Perspective of Dynamical Systems," *JAS*
  49(4):539-557 (2001) [ref 14] — Jovian-system temporary-capture dynamics,
  same domain.
- *Low-medium:* Parker, J.S. & Born, G.H., "Direct Lunar Halo Orbit
  Transfers," AAS 07-229 (2007) [ref 23]; "Modeling a Low-Energy Ballistic
  Lunar Transfer Using Dynamical Systems Theory," *J. Spacecraft and Rockets*
  45(6):1269-1281 (2008) [ref 24]; "Targeting Low-Energy Ballistic Lunar
  Transfers," AAS 10-1859 (2010) [ref 38] — the specific Parker & Born
  optimized-ΔV cross-check cited in Example 1 (`≈3.2 km/s` for a 185 km
  parking orbit) traces to ref [38]; all three are the same
  lunar-transfer-design lineage, not independently prioritized beyond the
  cross-check already extracted above. Topputo, F., Vasile, M. & Finzi, A.E.,
  "Combining Two and Three-Body Dynamics for Low Energy Transfer
  Trajectories," IAC-04-A.7.02 (2004) [ref 16]; Topputo, F., Vasile, M. &
  Bernelli-Zazzera, F., "Low Energy Interplanetary Transfers Exploiting
  Invariant Manifolds," *JAS* 53(4):353-372 (2005) [ref 18] — same
  low-energy-transfer lineage, tangential to this project's own cycler
  search focus. de Melo, C.F. & Winter, O.C., "Alternative Paths to
  Earth-Moon Transfer," *Math. Probl. Eng.* 2006:34317 (2006) [ref 21];
  García, F. & Gómez, G., "A Note on Weak Stability Boundaries," *CMDA*
  97(2):87-100 (2007) [ref 22] — same general lineage, lower priority.
  Prado, A.F., "Third-Body Perturbation in Orbits Around Natural Satellites,"
  *JGCD* 26(1):33-40 (2003) [ref 15] — general perturbation theory,
  tangential.
- Not flagged (general mission-context background): Howell, Barden & Lo 1997
  [ref 1] (already flagged elsewhere in `#730` §9 as the same paper); Lo et
  al. 2001 Genesis mission design [ref 2]; Angelopoulos 2008/in-press THEMIS/
  ARTEMIS mission descriptions [refs 3/4]; Woodard, Folta & Woodfork 2009
  ARTEMIS [ref 5]; Hénon 1969/1970 Hill-problem series [refs 6/7] — already
  covered by paper 3's own citation-mining list above (Hénon cluster, not
  re-flagged separately here to avoid double-counting); Hamilton & Burns 1991
  asteroid orbital-stability zones [ref 8]; Yamakawa et al. 1993 Earth-Moon
  gravitational-capture transfer [ref 9] — tangential, mission-specific.

---

## Combined citation-mining summary (all four papers, including paper 4 — completed by the coordinating session)

**Recurring citations flagged 2+ times across all four papers:**
1. **Villac & Scheeres 2003** (Hill-problem escape trajectories / periapsis
   Poincaré map origin) — flagged by BOTH paper 3 (Aydin-Batkhin) AND paper 4
   (Howell-Davis-Haapala, which credits it directly as the periapse-map
   method's own origin, ref [28]) — **confirmed 2-paper recurrence within
   this task**, plus the already-acquired `#683` digest independently names
   it as the periapsis-Poincaré-map originator (3-way corroboration total).
   **HIGH.**
2. **Frauenfelder & van Koert 2018 book** (*The Restricted Three-Body Problem
   and Holomorphic Curves*) — flagged by papers 1, 2, AND 3 (3x within this
   task alone). **HIGH** — the standard textbook underlying this entire
   symplectic-invariant lineage.
3. **Tsirogiannis, Perdios & Markellos 2009** (Hill-problem grid search) —
   flagged by paper 3 AND paper 4 within this task, AND by the
   already-completed `#742` digest's own citation-mining pass — **3-way
   recurrence** (2 within-task + 1 cross-task). **HIGH** (raised from
   Medium).
4. Broucke 1969 AIAA (`Stability of periodic orbits in the elliptic,
   restricted three-body problem`) — flagged by papers 1 and 2. **Medium** —
   note this is explicitly a *different* paper from the 4x-flagged Broucke
   1968 JPL TR database (`#730` §4 item 21); do not merge the two flags.
5. Howard & MacKay 1987 — flagged by papers 1 and 2. **Medium.**
6. Wonenburger 1966(?) — flagged by papers 1 and 2, same possible year typo
   in both source papers. **Medium.**
7. Ginzburg 2010 (Conley conjecture) — flagged by papers 1 and 3. **Medium.**
8. Bruno 1994 book (*The Restricted 3-Body Problem: Plane Periodic Orbits*)
   — flagged by papers 2 and 3. **Medium.**
9. Kalantonis 2020 — flagged by papers 2 and 3. **Low-medium.**
10. Gómez, Marcote & Mondelo 2005 — flagged by paper 3 here AND by the
    already-completed `#742` digest's own citation-mining pass (cross-task
    recurrence). **Medium.**
11. Anderson & Lo, "A Dynamical Systems Analysis of Resonant Flybys:
    Ballistic Case" — cited by paper 4 (ref [27]); this is `#730` §3 item 13,
    independently CrossRef-confirmed by this same coordinating session
    earlier in this task — cross-task corroboration of an already-tracked
    item, not a new gap.

**Single highest-priority genuinely-new candidate:** Aydin, C., "The
Conley-Zehnder indices of the spatial Hill three-body problem," *CMDA*
(2023), the direct predecessor paper 3 builds on — not independently
corroborated by a second flag, but structurally load-bearing (paper 3 cannot
be fully understood without it). Close second: Davis 2011 PhD dissertation
and Haapala 2011 MS thesis (paper 4's own refs [35]/[36]) — the direct
source theses underlying both paper 4 and the already-acquired `#683` journal
paper.

**Everything else** flagged individually above per-paper; none acquired, all
confirmed absent from `CORPUS_INDEX.md` and `#730` at the time of this pass.

## DOI corrections summary (for the coordinating session's `#730` update)

- §6 item 40 (Frauenfelder, Koh & Moreno 2023): `10.1137/22M1506743` →
  **`10.1137/22M1500459`**.
- §6 item 42 (Aydin & Batkhin 2025): `10.1007/s10569-025-10233-0` →
  **`10.1007/s10569-025-10241-7`** (directly visible in the PDF's own header).
- §6 item 41 (Frauenfelder & Moreno GIT-quotients) DOI
  `10.4310/jsg.2023.v21.n4.a3` was already correctly listed as "no DOI found"
  in the master list pre-acquisition; now confirmed **`10.4310/jsg.2023.v21.n4.a3`**
  and the title should be updated to the published form: "On GIT Quotients of
  the Symplectic Group, Stability and Bifurcations of Periodic Orbits (With a
  View Towards Practical Applications)" (arXiv preprint title differs: "...
  Bifurcations of **Symmetric** Orbits").
- §7 item 48 (Howell, Davis & Haapala 2012): DOI `10.1155/2012/351759` was
  already correct, no change needed.
