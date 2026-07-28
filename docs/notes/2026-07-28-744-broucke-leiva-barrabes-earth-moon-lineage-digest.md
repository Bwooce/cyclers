# Digest: Broucke 1968, Barrabés & Gómez 2003, Leiva & Briozzo 2006 "Full Atlas" (#744)

**Task:** `#744`, continuing the `#730` consolidated acquisition backlog
(`docs/notes/2026-07-27-730-acquisition-backlog-master-list.md`) §4 ("Casoliva/
Barrabés/Barrabés-Mondelo-Ollé Earth-Moon cycler methods" cluster). All three
PDFs were filed into the private `cyclers_pdf` corpus by the coordinating
session before this task started (currently untracked/uncommitted in that
repo — `git status` confirms all three as `??`; no commit hash to cite yet).
This digest does not commit anything in either repo.

**Filed:**
1. `broucke-1968-periodic-orbits-restricted-three-body-problem-earth-moon-masses-jpl-tr-32-1168-ntrs-19680013800.pdf`
2. `barrabes-gomez-2003-three-dimensional-pq-resonant-orbits-second-species-solutions-cmda-85-2-doi-10.1023-A1022098510161.pdf`
3. `leiva-briozzo-2006-earth-moon-cr3bp-full-atlas-low-energy-fast-periodic-transfer-orbits-preprint-arxiv-astro-ph-0612386.pdf`

**OCR status: all three text-layer, no OCR needed.** `pdftotext -layout`
yielded 610,661/69,731/110,635 chars over 100/30/198 pages respectively (6,107
/ 2,324 / 559 chars-per-page) — all far above the 10-char/page floor. Note
Broucke's 1968 JPL scan has moderate OCR-quality noise in its own embedded
text layer (garbled subscripts/superscripts, e.g. "L_" for "L₂", "m," for
"m₁") — normal for a 1968 mechanically-typeset technical report re-scanned
with an old text layer, not a defect in this pass; all quoted numeric tables
below were cross-read against the raw layout to resolve ambiguous digits.

---

## 1. Broucke 1968 — "Periodic Orbits in the Restricted Three-Body Problem With Earth-Moon Masses"

JPL Technical Report 32-1168 (Feb. 1968), 100pp, no DOI (NASA NTRS
19680013800). `#730` §4 item 21 — independently flagged 4 times before this
acquisition (`#725`, the 2026-06-11 RRT mining note, `#728`-kumar-moreno, and
`#742`'s Franz-Russell digest).

**Method:** a systematic numerical survey (Runge-Kutta, Adams-Moulton, and a
recurrent-power-series integrator, cross-checked against each other) of
**symmetric** periodic orbits in the planar RTBP at the Earth-Moon mass ratio
`µ = 0.012155099` (Earth = larger primary `m₁` at `x₁ = -µ`; Moon = smaller
primary `m₂` at `x₂ = 1-µ` — the paper states this convention explicitly:
"earth is taken as the larger of the primaries, the moon as the smaller
one"). 1811 orbits total, sorted into **10 families** (nine with a "known
origin" plus one supplementary family). Each orbit's initial conditions
`(x₀, ẏ₀, x₁, ẏ₁, E, T/2, µ, stability-index k, N-crossings)` are tabulated
(Tables 1-10, one per family). The paper also introduces a **six-class
symmetry taxonomy** (its own Section V intro) based on the ordering of the
four x-axis points `{x₀, x₁, m₁, m₂}` along the axis, each class mapped to a
"center" (`L1`–`L3`, `m₁`, `m₂`, or the barycenter) — this is the "classes 2,
3, 5" scheme the already-acquired Ross-Roberts-Tsoukkas mining note (2026-06-
11, line 60) says the RRT `(k1,k2)`-cycler orbits are classified under.

**The family taxonomy (verified directly against the paper's own Section V
text, not inherited from any secondary source):**

| Family | Around | Motion type | # orbits | Notes |
|---|---|---|---|---|
| G | L₁ | continuation of infinitesimal retrograde ellipses | 84 | branches with A₁ |
| I | L₂ | continuation of infinitesimal retrograde ellipses | 127 | — |
| J₁ | L₃ | continuation of infinitesimal retrograde ellipses | 134 | exact `µ→0` limit known |
| **A₁** | **m₁ (Earth)** | **retrograde**, starts infinitesimal-circular | 317 | end undetermined; 3 branch points with G, BD, J₁; alternating stable/unstable zones |
| BD | m₁ (Earth) | direct, starts infinitesimal-circular | — | 3 collision orbits; branch point with A₁ |
| E₁ | m₁ and m₂ jointly | direct (inertial) / retrograde (synodic), starts as a large joint circular orbit at infinite radius | — | all orbits have positive energy |
| F | m₁ and m₂ jointly | retrograde | — | Table 7 (full numeric table extracted, see below) |
| **C** | **m₂ (Moon)** | **retrograde**, starts infinitesimal-circular | 180 | 2 collision orbits, end undetermined; Broucke's own text: "probably exists for all mass ratios" |
| **H₁** | **m₂ (Moon)** | **direct**, starts infinitesimal-circular | 162 | beginning known, end unknown |
| **H₂** | **m₂ (Moon)** | **direct**, related to H₁ | 202 | neither end known; "some orbits nearly the symmetric image of...H₁" |

**Positive-control numeric values (Broucke's own text, `µ = 0.012155099`):**
`L₁ = (+0.836892919, 0)`, `L₂ = (+1.115699521, 0)`, `L₃ = (-1.005064527, 0)`.
A full initial-conditions table snippet (Table 7, family F, `µ`-column
confirms `0.012155085`–`0.012155098` across the family, matching the text's
statement that `µ` "changed slightly" between computation runs):
```
#   X0            YDOT0          X1            YDOT1         ENERGY       T/2           MASS         INDEX      N
1   .988193899   -8.667382311   -.029999312   10.642238346   1.256849114  .378057049   .012155092  4503.79807   1
5   .990407944   -3.586276255   -.147521118    3.881141335   .212458154  .628315156   .012155092    51.58237   1
36  .991999999   -3.043485945   -.251113745    2.968353953   .230265342  .765659596   .012155092   -51.91091   1
```
Full initial-conditions tables for all 10 families (Tables 1-10) are present
in the corpus PDF and can be extracted for a much larger positive-control set
if a specific family is needed for a future reproduction attempt.

### Cross-check: does Casoliva's and Kumar-Moreno's own characterization of Broucke's H1/H2/A1/C families match Broucke's actual text?

**Confirmed for the structural claims; one specific phrase in the already-
acquired Kumar-Moreno digest is the digest-writer's own paraphrase, not a
literal Broucke quote — flagged, not an error.**

- **Kumar-Moreno's "junction" claim (their digest headline result 1, citing
  "p.71 of his own text"):** Broucke's own Section V.I text (page 71 of the
  TR, confirmed by the page-footer "JPL TECHNICAL REPORT 32-1168 71" directly
  adjacent to this sentence) reads, verbatim: *"It is likely that, if one
  were to continue family H1 or the two open ends of H2, some junction
  between H1 and H2 would be found."* **Exact match** — this is precisely the
  "conjectured but couldn't compute" junction Kumar-Moreno's paper (57 years
  later) confirms.
- **A₁ = "small circular Earth retrograde orbits" (Kumar-Moreno's digest
  paraphrase):** confirmed for the *start* of the family — Broucke's own
  text: family A₁ "starts with infinitesimal circular orbits, radius `r₁`...
  retrograde in both the inertial and rotating axes," around `m₁` (Earth).
  Broucke's own text goes on to describe the family becoming rapidly complex
  (317 orbits, alternating stability zones, 3 branch points, "the family
  starts at the stability limit k=2... then moves slightly into the negative
  unstable zone... becomes stable again... continues to increase to a
  maximum of around 4000... decreases to cross the stable zone again...").
  Kumar-Moreno's own "chains via alternating Moon/Earth retrograde encirclements,
  ad infinitum" characterization is consistent with — but goes well beyond —
  what is in Broucke's own 1968 text, which only computed 317 orbits without
  reaching a natural end; the "ad infinitum" claim and the chaining mechanism
  is Kumar-Moreno's own 2025 continuation result, not something Broucke's
  1968 text itself asserts.
- **C family = "probably exists for all mass ratios" (Broucke's own words,
  confirmed verbatim) plus "end of the family has not yet been determined"
  (also verbatim)** — consistent with Kumar-Moreno's characterization of C as
  open-ended.
- **H1 "low/distant prograde orbits stretching toward L1" / H2 "highly
  eccentric distant-prograde orbits stretching toward L2" (Kumar-Moreno's
  digest paraphrase):** Broucke's own text describes H1's *first* orbits as
  "of the circular type around m₂ [Moon], with a direct motion" which then
  undergo "a complicated evolution... as can be seen in Fig. 33"; H2 "has the
  most complicated forms of orbits at both ends of the family," with a middle
  region "nearly circular and very similar to some of the first orbits of
  H1." **Grepped explicitly: the words "stretch," "toward L1," and "toward
  L2" do not appear anywhere in Broucke's own narrative text for these
  families.** The "stretching toward L1/L2" geometric characterization is
  Kumar-Moreno's own visual/geometric description of the orbit shapes (likely
  read off Broucke's own Figs. 33-34, which this pass did not independently
  re-render), not a phrase or claim present in Broucke's own 1968 prose. This
  is a genuine, if minor, instance of "verify against primary content, don't
  inherit secondary characterization" — the underlying family identification
  (H1=direct-around-Moon, H2=companion complex family, junction conjectured
  on p.71) is fully confirmed; the specific descriptive language "stretching
  toward L1/L2" should be attributed to Kumar-Moreno's own reading, not
  quoted as if it were Broucke's.
- **RRT mining note's "classes 2, 3, or 5" claim:** confirmed that Broucke's
  own Section V introduction defines exactly six symmetry classes (Class 1
  through Class 6, keyed to the ordering of `x₀, x₁, m₁, m₂` on the x-axis,
  with centers `L3, m1, barycenter, L1, m2, L2` respectively) — the RRT
  paper's specific claim that its own cyclers fall in classes 2/3/5 was not
  independently re-derived from RRT's own orbits in this pass (out of scope
  here), but the six-class *taxonomy itself* Broucke's text defines is
  confirmed to exist exactly as RRT's mining note describes it.

### Mandatory citation-mining pass (21 references, full read — all pre-1968 classical literature)

No cross-check with corpus was needed via secondary characterization here;
every reference was checked directly against `CORPUS_INDEX.md` and the
master list.

**Already in corpus:** none exactly — the only near-miss is Szebehely: ref 14
is **Szebehely, "Solution of the Restricted Problem of Three Bodies by Power
Series," Astron. J. 71:968-975 (1966)** — a *different*, shorter paper from
the already-corpus `szebehely-1967-theory-of-orbits...` textbook (confirmed
distinct by title/venue/year; the textbook is a book, this is a journal
article on the specific power-series integration method Broucke's own work
used).

**Genuinely new candidates (flagged, NOT acquired) — all pre-digital-era
classical literature, several without any modern DOI:**
- *High:* Hénon, M., "Exploration Numérique du Problème des Trois Corps," 4-
  part series, *Ann. Astrophys.* 28(3):499-511 and 28(6):992-1007 (1965),
  *Bull. Astron.* 1(1):57-80 and 1(2):49-66 (1966) [refs 15-18] — the classic
  equal-mass numerical-exploration series that established the modern
  periodic-orbit-family taxonomy this whole literature (including Broucke's
  own work and Hénon's later 1997 book already flagged elsewhere in `#730`
  item 29) builds on; directly comparable equal-mass counterpart to Broucke's
  own Earth-Moon-mass census.
- *Medium-high:* Deprit, A. & Henrard, J., "Natural Families of Periodic
  Orbits," *Astron. J.* 72(2):158-172 (1967) [ref 11] — the branch-
  point/family-continuation concept Broucke's own Section V explicitly
  invokes ("shrinkage points," "manifold" terminology) is credited to this
  paper; Deprit, A. & Price, J.F., "The Computation of Characteristic
  Exponents in the Planar Restricted Problem of Three Bodies," *Astron. J.*
  70(10):836-846 (1965) [ref 13] — the specific stability-index computation
  method Broucke's own text says was used for "most orbits" in the tables;
  Barrar, R.B., "Existence of Periodic Orbits of the Second Kind in the
  Restricted Problem of Three Bodies," *Astron. J.* 70:3-4 (Feb. 1965)
  [ref 6] — directly on-point existence theory for second-species orbits,
  the same object class the Barrabés-Gómez in/out-map lineage (papers 2-3 of
  this digest, and the already-acquired 2002 paper) formalizes analytically;
  Arenstorf, R.F., "Periodic Solutions of the Restricted Three-Body Problem
  Representing Analytic Continuations of Keplerian Elliptic Motions," *Am.
  J. Math.* 85(1):27-35 (1963) [ref 7] — classical existence/continuation
  theory directly underlying the whole family-continuation approach.
- *Medium:* Stromgren, E., "Connaissance Actuelle des Orbites dans le
  Problème des Trois Corps," Copenhagen Observatory Publication 100 (1935)
  [ref 1] — the foundational equal-mass classification (Stromgren's classes
  a-g) Broucke's own text compares every one of its 10 families against
  throughout (e.g. "our family G corresponds to 'class c'... our family C
  corresponds with Stromgren's class f"); likely inaccessible pre-digital
  monograph, but the single most load-bearing cross-reference in Broucke's
  own text. Szebehely 1966 (see "already in corpus" note above, distinct
  from the in-corpus 1967 textbook) [ref 14].
- *Low-medium:* Wintner, A., *The Analytical Foundations of Celestial
  Mechanics*, Princeton (1941) [ref 5], and Wintner, A., "Upon the
  Characteristic Exponents in the Stromgrenian Groups of Periodic Orbits,"
  Copenhagen Obs. Publication 78 / *Am. J. Math.* 53(3):611-616 (1931)
  [ref 10] — source of "shrinkage point" terminology Broucke's own text
  uses; Moulton, F.R., *Periodic Orbits*, Carnegie Institution (1920)
  [ref 2]; Darwin, G.H., *Periodic Orbits*, Scientific Papers Vol. 4,
  Cambridge (1911) [ref 3] — Broucke's own text explicitly notes his H1/H2
  families resemble Darwin's "satellites A, B, C" (`µ=1/11`, Sun-Jupiter
  case); Rabe, E., "Periodic Librations About the Triangular Solutions of
  the Restricted Earth-Moon Problem and Their Orbital Stabilities," *Astron.
  J.* 67(10):732-739 (1962) [ref 12] — same numerical (recurrent power
  series) method lineage Broucke's own text credits.
- *Low:* Poincaré, H., *Les Méthodes Nouvelles de la Mécanique Céleste* Vol.
  1 (1899) [ref 4] — extremely foundational, but a different work from any
  Poincaré item already indexed (9 unrelated `CORPUS_INDEX.md` hits are all
  incidental mentions of "Poincaré section," not this book); Steffensen,
  J.F., "On the Restricted Problem of Three Bodies," Mat. Fys. Medd. Dan.
  Vidensk. Selsk. 30(18) (1956) [ref 8]; Broucke, R.A., "Recherches
  d'Orbites Périodiques dans le Problème Restreint Plan (Système
  Terre-Lune)," PhD dissertation, Univ. of Louvain (1963) [ref 9] — Broucke's
  own thesis, almost certainly fully superseded by this TR; Bartlett, J.H.,
  "The Restricted Problem of Three Bodies (I)," Mat.-Fys. Skr. 2(7) (1964)
  [ref 19], and Bartlett, J.H. & Wagner, C.A., "...(II)," Mat.-Fys. Skr.
  3(1) (1965) [ref 20].

None of these 21 references are acquired or otherwise flagged elsewhere in
`#730`/`CORPUS_INDEX.md` — all confirmed absent by direct grep (a false-
positive "Rabe" hit count in `CORPUS_INDEX.md`/master-list turned out to be
substring matches inside "Barra**bés**," not this Rabe 1962 paper).

---

## 2. Barrabés & Gómez 2003 — "Three-Dimensional p–q Resonant Orbits Close to Second Species Solutions"

*Celestial Mechanics and Dynamical Astronomy* 85(2):145-174 (2003). DOI
`10.1023/A:1022098510161`. `#730` §4 item 26 (HIGH — the explicit spatial/3D
companion the already-acquired 2002 paper's Conclusions flag as
"forthcoming").

**Method:** extends the already-acquired Barrabés & Gómez 2002 paper's
planar analytic in/out-map matched-asymptotics construction to the full
**spatial** (3D) RTBP. Initial conditions on a sphere `B` of radius `µ^α`
around the small primary `M` are parameterized in spherical coordinates
`(ϕ, θ)` for position and `(φ, ψ, vᵢ)` for velocity (Eq. 5). The paper derives
the **out-map** (Section 3, via Proposition 3.1: the true RTBP flow equals
the two-body flow to `O(µ^(1-α))` under a closeness-to-M condition) and the
**in-map** (Sections 4-5, analogous backward-time construction), then
studies the **matching equations** (Section 6, systems (41)-(42)) between
in- and out-maps to `O(µ^α)`, yielding initial conditions for orbits
"periodic" to error `O(µ^(1-α))`, `α ∈ (1/3, 1/2)` — identical error order to
the 2002 planar paper.

**Key result:** Section 6 splits explicitly into **§6.1 Planar Case** (`ϕ =
φ₀ = 0`) and **§6.2 Spatial Case** (`φ₀ = ±π/2`). In the planar case the
general spatial matching system (42) — which the paper shows has all
off-diagonal terms `λᵢⱼ` vanish identically when `ϕ = φ₀ = 0` — reduces to a
single equation (45), `E(θ, ψ₀) = 0`, plus the resonance condition (46)
relating `Cⱼ`, `ψ₀`, and `p/q`. Figure 3 shows worked planar 1-2, 2-1, 2-3,
and 3-2 resonant orbit examples with explicit Jacobi-constant values (usable
as a small positive-control set): `Cⱼ = -0.850431` and `1.059752` (1-2
resonant), `Cⱼ = 0.303724` and `0.565699` (2-1 resonant). Figure 4 gives
spatial examples (1-2, 1-3, 1-5, 3-5 resonant) at `φ = 0, θ = π/2`, with
initial conditions given in closed form: `x = µ-1, y = µ^α cos ϕ, z = µ^α sin
ϕ, ẋ = ẏ = 0, ż = v sin φ₀`.

### Cross-check: does this paper's spatial construction reduce to the planar 2002 paper's construction when spatial parameters → 0?

**Confirmed directly from this paper's own equations — yes, exactly, and
this paper says so explicitly, not just structurally.** Section 6.1's own
text states: *"we will take ϕ = 0 and φ₀ = 0. In this case, we have seen
that the outer and the inner maps match up to terms of order µ^α... After
substituting the values ϕ = 0 and φ₀ = 0 in both equations we get [Eq. 44]...
(Observe that φ has disappeared. We are dealing with plane orbits.)"* — the
paper's own general spatial matching system (42) is shown term-by-term to
collapse to the identical planar matching system when the spatial angles
vanish, and the paper is explicit that this is the *same* planar problem
solved by the earlier work: Section 1's own text states *"we will make use
of some results obtained in a previous paper (Barrabés and Gómez, in press),
where we studied the set of initial conditions... which correspond to
spatial p-q resonant orbits. Some of these results will be summarized in
Section 2"* — i.e., this 2003 paper explicitly imports its own zero-order
restrictions (Eqs. 7-12, credited "see Barrabés and Gómez, in press, for the
details") directly from the (then still-forthcoming, now already-acquired)
2002 paper. The continuity claim is confirmed, not merely plausible — this
paper is a direct, explicit, equation-level generalization of the 2002
paper's planar construction, and the 2002 paper's planar sub-case reappears
here verbatim as §6.1.

### Mandatory citation-mining pass (6 references, full read — a very short list)

**Already in corpus:** Szebehely 1967 (textbook, in corpus, digested).
**Already flagged elsewhere in `#730`:** Hénon 1997 (`#730` item 29); Yen
1985 (AAS 85-346, Mercury gravity-assist mission design — flagged low
priority/out-of-scope in the already-acquired `#742` digest's citation-
mining of the 2002 paper, same reference reused here for the same "reverse
V-EGA" motivational aside).

**Genuinely new candidates (flagged, NOT acquired):**
- *Low-medium:* Stiefel, E. & Scheifele, G., *Linear and Regular Celestial
  Mechanics*, Springer Grundlehren 174 (1971) — the standard KS-
  regularization textbook; not previously flagged anywhere in `#730`/
  `CORPUS_INDEX.md`; background/tool reference, not domain-specific to
  cyclers.
- *Low:* Barrabés, E., "Òrbites de segona espècie del problema espacial de 3
  cossos," PhD thesis, Universitat Autònoma de Barcelona (2001) — same
  thesis already flagged low-priority in the already-acquired `#742` digest's
  mining of the 2002 paper (superseded by the published 2002/2003 CMDA
  papers, both now in corpus); not re-flagged as new, cross-confirmed only.
- Font, Nunes & Simó — this paper cites it as **"2002, Nonlinearity 15,
  115-142"** (a full, published citation), whereas the already-acquired 2002
  paper's own reference list (per the `#742` digest) cited it as **"2001...
  Nonlinearity (to appear)"** — same paper, now confirmed published; a minor
  citation-completeness note, not a new gap (still tangential quasi-collision
  content, still low priority, not flagged for acquisition).

No genuinely new high-priority gap surfaced by this paper's own (very short)
reference list.

---

## 3. Leiva & Briozzo 2006 — "The Earth-Moon CR3BP: A Full Atlas of Low-Energy Fast Periodic Transfer Orbits"

arXiv:astro-ph/0612386 (Dec. 15, 2006), 198pp. **This is NOT `#730` §4 item
22** (Leiva & Briozzo, "Control of chaos and fast periodic transfer orbits in
the Earth-Moon CR3BP," *Acta Astronautica* 58(8):379-386 (2006), DOI
`10.1016/j.actaastro.2005.12.006`) — that paper **remains unacquired and
paywalled**. This arXiv paper is a separate, later, companion/follow-up
piece, confirmed by three independent lines of evidence directly in its own
text:

1. Page 1 states: *"This is the accompanying article to a shortened version
   submitted to Celestial Mechanics and Dynamical Astronomy in December
   2006."* — CMDA, not Acta Astronautica.
2. The Introduction explicitly cites prior work as a *separate, already-
   published* paper: *"In a recent work (Leiva and Briozzo 2006a) we
   implemented this program for periodic transfer orbits in the Earth-Moon
   planar CR3BP. By a numerical survey we found 287 orbits..."* — then goes
   on to describe *this* paper's own different, larger result (80 families,
   not 287 orbits from an initial numerical survey).
3. **The paper's own reference list settles it directly**: it cites itself as
   *"Leiva, A.M. and Briozzo, C.B.: 2006b, 'The Earth-Moon CR3BP: an Atlas of
   low-energy fast periodic transfer orbits', Celest. Mech. Dynam. Astron.
   (submitted)"* — i.e. **this arXiv paper is "2006b," the CMDA submission**
   — and separately cites *"Leiva, A.M. and Briozzo, C.B.: 2006a, 'Control of
   chaos and fast periodic transfer orbits in the Earth-Moon CR3BP', Acta
   Astronautica 58, 379-386"* as the distinct, already-published prior paper.
   287 orbits (2006a, Acta Astronautica) vs. 80 families (2006b, this paper,
   the CMDA submission) are quantitatively different results from two
   different papers, not a preprint/journal pair of the same paper. A
   CrossRef author search for Leiva+Briozzo found no CMDA paper matching
   this "Atlas" title or content — apparently the "shortened version" 2006b
   submission was never published under a matching title (the only two CMDA
   papers by these authors are unrelated 2005/2008 QBCP papers, DOIs
   `10.1007/s10569-004-7818-3` and `10.1007/s10569-008-9134-9`).

**Method:** a numerical survey (Section 3) for periodic orbits (POs)
encircling both Earth and Moon ("transfer orbits" in the paper's own loose
sense — encircling both primaries, not necessarily passing close to either)
in the planar Earth-Moon CR3BP at `µ = 0.0121505`, followed by analytical
continuation along each family in the energy parameter `h = -C/2` on a
Poincaré section `Σ1 = {x = 0.836915310 [i.e. the L1 x-coordinate], y, vx >
0, vy}`, augmented by symmetry-completion (Sec. 3.3), period-doubling
bifurcation search (Sec. 3.4), and continuation-in-µ (Sec. 3.5) techniques to
find additional families beyond the initial 287-orbit grid search. Result:
**80 families** (many with multiple symmetric/asymmetric branches, e.g.
family 197 has branches A-D), most orbits **asymmetric** with respect to the
Earth-Moon axis — explicitly contrasted with the "recent similar work by
Bruno and Varin (2006)" which the paper says covers only symmetric orbits.

**Key result / positive-control data:** the abstract's promised "numerical
data for the intersection of an orbit with Σ1 at a reference value of h" is
**Table 2** ("Initial conditions for a reference orbit on each segment of
each family in the Atlas"), giving `(h, T, x, y, vx, vy)` to 16 significant
figures for every family/branch. Sample extracted rows (all at the fixed
`x = 0.83691530956968` — the L1 Poincaré-section x-coordinate):
```
Family 357 (Ss): h=-1.553849931959387  T=15.35213364809199
                 x=0.8369153095696800  y=-0.1171235689440371
                 vx=0.1882861991773726 vy=-0.05721969437090824
Family 037 (Aa): h=-1.557940386649792  T=15.50020353331435
                 x=0.8369153095696800  y=-0.09526883178836698
                 vx=0.1931044766332113 vy=-0.06987584622079623
```
**Table 1** gives symmetry type and minimal period `Tmin` for every family/
branch (80+ rows), ordered by ascending `Tmin` (14.49 to 36.5+ days in the
sampled range). Each family's characteristic curves `T(h)`, `y(h)`, `vy(h)`,
`vx(h)` on `Σ1`, a stability-index curve, and `h(xi)` on a second Poincaré
section `Σ2 = {x, y=0, vx, vy>0}` are also described (Section 4) and plotted
per-family (Section 5), matching the abstract's promise exactly.

### Cross-check: which Leiva-Briozzo paper does Casoliva's "families 357 and 037" citation actually point to?

**Resolved — this is the correct, exact paper.** The already-acquired `#725`
Casoliva digest states: *"Casoliva explicitly cross-validates some of these
against Leiva & Briozzo's atlas families 357 and 037 (Fig. 7/Fig. 16)"* and
separately identifies the source as *"'The Earth-Moon CR3BP: A Full Atlas of
Low-Energy Fast Periodic Transfer Orbits' (2006 preprint/ref [18] 2008)"* —
**this arXiv paper's own Table 1 lists exactly families "357" (Ss symmetry,
Tmin=14.487454) and "037" (Aa symmetry, Tmin=14.786812)** as the two
shortest-period families in the entire Atlas. This confirms, via this
paper's own family-numbering scheme (not inherited from Casoliva's or the
`#725` digest's characterization), that Casoliva's citation targets *this*
Full Atlas paper specifically — not the still-unacquired Acta Astronautica
"Control of Chaos" paper (item 22, which reports 287 numbered orbits from an
initial search, no letter-branch family subdivisions) and not the CMDA 101
(2008) QBCP-extension paper (item 23, about extending orbits to a different
dynamical model, not about periodic-orbit family cataloguing). The `#725`
digest's original identification was already correct; this pass confirms it
directly against the primary source's own numbering.

### Mandatory citation-mining pass (19 references, full read)

**Already in corpus:** Andreu, M.A., 1998 PhD thesis, "The Quasi-bicircular
Problem" (already digested, `2026-06-14-andreu-quasi-bicircular-digest.md`);
Ross, S., Koon, W.S., Lo, M.W. & Marsden, J.E., "Design of a Multi-Moon
Orbiter," AAS 03-143 (2003) (already digested, `#498`); Szebehely 1967 (in
corpus, digested); Hénon 1997 (already flagged, `#730` item 29).

**Genuinely new candidates (flagged, NOT acquired):**
- *High:* Bruno, A.D. & Varin, V.P., "On families of periodic solutions of
  the restricted three-body problem," *Celest. Mech. Dynam. Astron.*
  95(1-4):27-54 (2006) — this paper's own text explicitly frames it as a
  cross-check counterpart: *"a recent similar work by Bruno and Varin (2006)
  concerns only symmetric orbits in the CR3BP (though for 0 ≤ µ ≤ 1/2). We
  think the intersection between both works provides a desirable cross-check
  of the results"* — directly relevant, same-domain, explicitly flagged by
  the source paper itself as a cross-check opportunity not yet exploited by
  this project.
- *Medium-high:* **A probable additional gap, distinct from `#730` item 23**
  — this paper's reference list cites *"Leiva, A.M. and Briozzo, C.B.: 2005,
  'Fast Periodic Transfer Orbits in the Sun-Earth-Moon Quasi-Bicircular
  Problem', Celest. Mech. Dynam. Astron. 91 (3-4), 357-372"* as a "related
  work" showing "how under suitable conditions some UPOs could be extended
  from the Earth-Moon CR3BP to... QBCP." This is a **different volume/year**
  (CMDA 91, 2005) from `#730` §4 item 23 (*"Extension of fast periodic
  transfer orbits from the Earth-Moon RTBP to the Sun-Earth-Moon
  Quasi-Bicircular Problem," CMDA 101:225-245 (2008)*) — same author pair,
  overlapping topic (QBCP extension of Earth-Moon UPOs), but a title, volume,
  and year all distinct enough (91→101, 2005→2008, a normal ~2-year CMDA
  volume progression) to plausibly be two separate papers in the same
  lineage (an initial 2005 result, then a fuller 2008 extension), not a
  single mis-cited paper. **Flagged as a probable new, additional gap
  alongside the already-known item 23**, not a substitute for it — should be
  independently verified (e.g. via CrossRef) before either is acquired.
- *Medium:* Yagasaki, K., "Computation of low energy Earth-to-Moon transfers
  with moderate flight time," *Physica D* 197:313-331 (2004a); Yagasaki, K.,
  "Sun-perturbed Earth-to-Moon Transfers with Low Energy and Moderate Flight
  Time," *Celest. Mech. Dynam. Astron.* 90(3-4):197-212 (2004b) — both
  directly on-point low-energy Earth-Moon transfer-orbit design papers,
  explicitly discussed and compared against in this paper's own §6.2
  Applications discussion; Andreu, M.A., "Dynamics in the Center Manifold
  Around L2 in the Quasi-Bicircular Problem," *Celest. Mech. Dynam. Astron.*
  84(2):105-133 (2002), and Andreu, M.A., "Preliminary Study on the
  Translunar Halo Orbits of the Real Earth-Moon System," *Celest. Mech.
  Dynam. Astron.* 86(2):107-130 (2003) — companion journal papers to the
  already-corpus Andreu 1998 QBCP thesis, not previously flagged as
  standalone acquisitions.
- *Low-medium:* Bollt, E. & Meiss, J.D., "Targeting Chaotic Orbits to the
  Moon Through Recurrence," *Phys. Lett. A* 204:373-378 (1995), and Schroer,
  C.G. & Ott, E., "Targeting in Hamiltonian systems that have mixed
  regular/chaotic phase spaces," *Chaos* 7(4):512-519 (1997) — earlier
  chaos-based low-energy-transfer methods this paper's own Introduction
  contrasts itself against ("The disadvantages of their approach are that it
  requires considering large quantities of orbit arcs..."); Otani, M. &
  Jones, A.J., "Guiding chaotic orbits," Imperial College research report
  (1997) — unpublished report, low priority; Ott, E., Grebogi, C. & Yorke,
  J.A., "Controlling chaos," *Phys. Rev. Lett.* 64:1196-1199 (1990) — the
  foundational chaos-control physics paper motivating this whole
  methodology, but general-physics background rather than RTBP-specific.
- Parker, T.S. & Chua, L.O., *Practical Numerical Algorithms for Chaotic
  Systems*, Springer (1989); Press, Teukolsky, Vetterling & Flannery,
  *Numerical Recipes*, 2nd ed. (1992); Verhulst, F., *Nonlinear Differential
  Equations and Dynamical Systems*, Springer (1990); Arnold, Kozlov &
  Neishtadt, "Mathematical Aspects of Classical and Celestial Mechanics,"
  in *Dynamical Systems III* (1993) — general tool/textbook references, not
  flagged (infrastructure/background, consistent with this project's prior
  triage convention for such references).

---

## Summary answers (for the dispatching session)

- **Paper 1 (Broucke 1968):** confirms the H1/H2-junction conjecture quote
  verbatim ("It is likely that, if one were to continue family H1 or the two
  open ends of H2, some junction between H1 and H2 would be found," p.71) and
  the structural family taxonomy (A1=Earth-retrograde, C=Moon-retrograde,
  H1/H2=Moon-direct) exactly as Kumar-Moreno's digest and the RRT mining note
  describe. One nuance: Kumar-Moreno's "stretching toward L1/L2" phrase for
  H1/H2 is their own geometric paraphrase (from Broucke's figures), not a
  literal quote from Broucke's own narrative text (searched, zero hits) —
  noted, not a substantive error. Citation-mining of Broucke's own 21
  references surfaced 21 new (all pre-1968 classical) candidates, none
  previously flagged, headlined by the 4-part Hénon 1965-66 numerical-
  exploration series.
- **Paper 2 (Barrabés & Gómez 2003):** confirmed, directly from this paper's
  own equations (§6.1) and its own text ("Observe that φ has disappeared. We
  are dealing with plane orbits"), that the general spatial matching system
  reduces exactly to the already-acquired 2002 paper's planar construction
  when the spatial angles vanish — the continuity claim in the task brief is
  correct, not merely plausible. Citation-mining (6 refs) surfaced one new
  candidate (Stiefel & Scheifele 1971 regularization textbook), low priority.
- **Paper 3 (Leiva & Briozzo 2006, "Full Atlas"):** identity correction
  confirmed on all three lines of evidence in the paper's own text (page-1
  statement, Introduction's explicit "2006a... 287 orbits" prior-work
  citation, and the paper's own self-citation as "2006b...(submitted)" in its
  reference list) — this is a distinct paper from `#730` item 22 (still
  unacquired/paywalled), reporting 80 orbit families (not 287 orbits) with
  full Table 1/Table 2 family-catalogue and reference-orbit numeric data
  extracted above. Also resolved a standing ambiguity: Casoliva's "Leiva &
  Briozzo atlas families 357 and 037" citation is confirmed, via this paper's
  own Table 1 family numbers, to target this exact paper. Citation-mining (19
  refs) surfaced 9 new candidates, headlined by Bruno & Varin 2006 (an
  explicit self-flagged cross-check counterpart) and a probable *additional*
  (not duplicate) gap — Leiva & Briozzo 2005, CMDA 91:357-372 — distinct from
  the already-known `#730` item 23 (CMDA 101:225-245, 2008).
- Across all three papers, citation-mining surfaced **31 new candidates
  total** (21 + 6 + 9, with de-duplication already applied within each
  paper's own list; no cross-paper duplicates found), none acquired this
  pass. No genuinely new item rose to the same corroboration level as the
  now-4x-flagged Broucke 1968 gap (which this task itself closes).
