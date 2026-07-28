# Digest: Anderson-Campagnola-Koh-McElrath-Woollands 2021, Barrabés-Mondelo-Ollé 2009, Barrabés-Gómez 2002, Franz-Russell 2022 (#742)

**Task:** `#742`, continuing the `#730` consolidated acquisition backlog
(`docs/notes/2026-07-27-730-acquisition-backlog-master-list.md`). The user
directly supplied PDFs for §3 item 11 and §4 items 24/25/44 of that list — all
four independently verified page-1-match by the coordinating session before
this task started.

**Filed** (private `cyclers_pdf` repo, commit `513e657`):
1. `anderson-campagnola-koh-mcelrath-woollands-2021-endgame-europa-lander-ganymede-europa-approach-jas-68-96-doi-10.1007-s40295-021-00250-7.pdf`
2. `barrabes-mondelo-olle-2009-numerical-continuation-homoclinic-connections-periodic-orbits-rtbp-nonlinearity-22-2901-doi-10.1088-0951-7715-22-12-006.pdf`
3. `barrabes-gomez-2002-spatial-pq-resonant-orbits-rtbp-cmda-84-387-doi-10.1023-A1021137127909.pdf`
4. `franz-russell-2022-database-planar-3d-periodic-orbits-families-near-moon-jas-69-1573-doi-10.1007-s40295-022-00361-9.pdf`

**OCR status: all four text-layer, no OCR needed.** All four are Acrobat-Distiller-produced (native LaTeX/Word PDFs); `pdftotext -layout` yielded
51.6k/76.0k/49.4k/113.6k chars over 24/19/21/40 pages respectively — all far
above the 10-char/page floor.

---

## 1. Anderson, Campagnola, Koh, McElrath & Woollands 2021 — "Endgame Design for Europa Lander: Ganymede to Europa Approach"

*J. Astronaut. Sci.* 68(1):96-119. `#730` §3 item 11 (HIGH, 2 independent flags: `#727`, `#728`-oberon-survey).

**Method:** designs the endgame from the last Ganymede flyby to final Europa
approach in three stages, all in the planar CR3BP (CRTBP): (1) a "portal"
(topological disk on a cylinder exterior to L2, from prior isolating-block work
[11]) interfacing the Europa-orbit-insertion/surface leg to (2) a Poincaré
surface of section where the final resonance is computed, then (3) two
resonant legs connected by two impulsive ΔVs applied at trajectory apoapses,
searched over a grid of apoapsis-location combinations `(k1,k2)` and Jacobi
constants. The CRTBP solution is then handed to JPL's COSMIC optimizer (SNOPT
+ DIVA/Monte propagator, `jup310.bsp` ephemeris) as a starting iterate to
obtain a full-ephemeris optimized solution and to bound the "true" attainable
solution set against the CRTBP-only set.

**Key result (the `k1=3, k2=1` example case, C=3.0024):**
- CRTBP (unoptimized): total ΔV = **146.71 m/s**, split as Burn 1 = 50.81 m/s
  at 14.45 days, Burn 2 = 96.12 m/s at 14.45+27.45 days from the trajectory's
  start (Table 2).
- COSMIC (ephemeris-optimized): total ΔV = **135.04 m/s**, split as Burn 1 =
  31.42 m/s at 14.57 days, Burn 2 = 103.62 m/s at 27.72 days.
- Conclusions section states these rounded: "a ΔV of approximately 146 m/s...
  converged to a solution with a ΔV of approximately 135 m/s."
- Broader ΔV sweep: "total ΔVs vary significantly with the lower values
  approaching 140 m/s" (text, p.109) — 146.71 m/s is near, not at, the sweep's
  practical floor.

### Cross-check: Kumar/Anderson/de la Llave 2023's own benchmark claim

**Verified directly against this paper's text.** The already-acquired `#727`
Kumar-Anderson-de la Llave 2023 (Acta Astronautica 211:76-87) paper's own §5.3
text (confirmed via `pdftotext` on the corpus PDF) reads: *"The results from
[6] are more in line with typical approach design methods, and they give a ∆v
a little less than 150 m/s and a TOF of approximately 40 days."* Reference
`[6]` in that paper resolves, byte-for-byte, to this exact paper (same
authors, same title, same JAS 68:96-119 citation).

- **ΔV: confirmed, precisely.** 146.71 m/s (CRTBP) and 135.04 m/s (COSMIC) are
  both "a little less than 150 m/s" — Kumar's characterization is accurate
  and slightly conservative (the true numbers are closer to 135-147 m/s, not
  just "a little less" than 150).
- **TOF: plausible but not a literal quote from Anderson's text.** Anderson
  2021's own Table 2 gives only the two *burn* times (14.45/27.45 days CRTBP,
  14.57/27.72 days COSMIC), measured from the start of the (COSMIC-side)
  segment; no single "total TOF" figure is printed anywhere in the paper's
  text (exhaustive grep of the extracted text for "day(s)" found only these
  four numbers plus an unrelated 150-day *encounter-opportunity* window
  constraint in Fig. 18's caption, which is not the transfer TOF). Since Burn
  2 (the later of the two impulses) already occurs at ~27.5-27.7 days and the
  trajectory continues past it to the actual Europa portal/orbit-insertion
  point, a total transfer duration in the high-30s/low-40s-days range is
  geometrically consistent with the printed burn schedule, but "~40 days" is
  Kumar et al.'s own paraphrase/rounding of Table 2 plus the untabulated
  final coast segment, not a number Anderson's paper states verbatim.
  **No correction is needed to `#727`'s digest** — its own text already
  quotes Kumar's characterization verbatim without independently re-deriving
  it from Anderson's Table 2, so it was already appropriately hedged; this
  pass adds the precision that the "~40 days" figure specifically does not
  appear as a printed total in Anderson 2021's own text.

### Mandatory citation-mining pass (30 references, full intro + full list read)

**Already in corpus:**
- Campagnola & Russell, "Endgame problem part 2: Multibody technique and the
  Tisserand-Poincaré graph," *JGCD* 33(2):476-486 (2010) [ref 15] — this
  **resolves the open ambiguity flagged in `#730` §1** ("the exact 'Campagnola
  & Russell 2010' title/venue `#722` had in mind was not independently
  re-verified"). This paper's own reference confirms the JGCD 2010 title
  exists and is the journal twin of the already-corpus
  `campagnola-russell-2009-endgame-partB-multibody-tp-graph-AAS-09-227.pdf`
  (same "Part 2/Part B," same Tisserand-Poincaré-graph content, conference →
  journal). `#730` §1's "likely" duplicate flag can now be marked **confirmed**
  duplicate, not merely probable.
- Campagnola, Buffington & Petropoulos 2014, "Jovian tour design for orbiter
  and lander missions to Europa," *Acta Astronaut.* 100:68-81 [ref 18] —
  already digested (`2026-06-17-digest-campagnola-2014.md`).

**Already flagged elsewhere in `#730` (no new action, cross-confirmed
recurring):** refs 2 (Anderson 2015, §3 item 16), 4 (Anderson & Lo 2009, item
15), 5 (Anderson & Lo 2010, item 12), 6 (Anderson & Lo 2011 ballistic, item
13), 7 (Anderson & Lo 2014 spatial, item 20).

**Genuinely new candidates, topically overlapping (flagged, NOT acquired):**
- *High-medium:* Lantoine, Russell & Campagnola 2011, "Optimization of
  low-energy resonant hopping transfers between planetary moons," *Acta
  Astronaut.* 68(7-8):1361-1378 [ref 24] — **independently re-flagged by
  Franz & Russell 2022 below** (their ref 26), now a 2-flag recurrence in this
  very task; Russell & Lam 2007, "Designing ephemeris capture trajectories at
  Europa using unstable periodic orbits," *JGCD* 30(2):482-491 [ref 28] —
  directly on-point (unstable periodic orbits + Europa capture).
- *Medium:* Anderson & Parker 2012/2013 (ballistic lunar-surface transfer
  survey/comparison, JGCD 35(4)/CMDA 115(3)) [refs 8,9]; Anderson, Easton & Lo
  2017/2019/2020 isolating-blocks/isolating-neighborhoods series (Physica D,
  CMDA, AAS 19-744) [refs 10-13] — the theoretical machinery behind this
  paper's own "portal" construction; Campagnola, Skerritt & Russell 2012,
  "Flybys in the Planar, Circular, Restricted, Three-Body Problem," *CMDA*
  113(3):343-368 [ref 16]; Campagnola et al. 2014a, "Tisserand leveraging
  transfers," *JGCD* 37(4):1202-1210 [ref 17]; Grebow, Petropoulos & Finlayson
  2011, "Multi-body capture to low-altitude circular orbits at Europa," AAS
  proceedings vol. 142 [ref 22].
- *Low-medium:* Anderson 2005 PhD thesis (predecessor to the whole lineage)
  [ref 1]; Anderson & Lo 2005b AAS conference precursor [ref 3]; Buffington
  2014, Europa Clipper trajectory design (mission-specific) [ref 14];
  Johannesen & D'Amario 1999, Europa Orbiter mission trajectory design
  (historical) [ref 23]; Parker & Anderson 2013 (*Acta Astronaut.* 84:1-14)
  and Parker, Anderson & Peterson 2013 (*JGCD* 36(5):1501-1511), ballistic
  lunar-orbit transfer surveys [refs 25,27]; Woolley & Scheeres 2010,
  "Hyperbolic periodic orbits in the three-body problem and their application
  to orbital capture," AAS George H. Born Symposium [ref 30].
- Parker & Anderson 2014 Wiley book [ref 26] and Sweetser et al. 1997 [ref 29]
  are already flagged (low priority) in `#727`'s own citation-mining list —
  not re-flagged here.
- Folkner et al. 2014 DE430/431 ephemeris [ref 20], Finlayson 1999 PTool doc
  [ref 19], and Gill/Murray/Saunders 1999 SNOPT user's guide [ref 21] are
  tool/data references, not research results — out of scope, not flagged.

---

## 2. Barrabés, Mondelo & Ollé 2009 — "Numerical continuation of families of homoclinic connections of periodic orbits in the RTBP"

*Nonlinearity* 22(12):2901-2918. `#730` §4 item 24 (HIGH — "the core
homoclinic-continuation algorithm underlying all of Class 2").

**Method (§2.2):** given a hyperbolic periodic orbit (here, an L1/L2/L3
Lyapunov orbit) with 1D stable/unstable eigenspaces, parameterize its
manifold tubes linearly as `ψ^{u/s}(θ,ξ) = φ_{θT/2π}(z0) + ξ λ^{-θ/2π}
Dφ_{θT/2π}(z0) v^{u/s}` (Eq. 4) and define `F(θ^u,θ^s) = P+(ψ^u) - P-(ψ^s)`
on a shared Poincaré section — a homoclinic connection is a zero of `F`.
Rather than differentiate eigenvalues/eigenvectors with respect to a free
initial condition `z0` (numerically fragile through LAPACK-class
eigen-solvers), the eigenvector condition and a normalization are folded
directly into the nonlinear system as extra equations, giving an
**over-determined, rank-deficient system in `(h, z, T, λ^u, v^u, λ^s, v^s,
θ^u, T^u, θ^s, T^s)`** (the paper's system labeled (6)/referenced throughout
as "system (7)"): `H(z)-h=0`, `g1(z)=0`, `φ_T(z)-z=0`, the two eigenvector
equations `Dφ_T(z)v^{u,s} - λ^{u,s}v^{u,s}=0`, two unit-norm equations, two
Poincaré-section-matching equations `g2(φ_{T^{u,s}}(ψ^{u,s}))=0`, and the
final manifold-matching equation. Solved by Newton's method with a
**minimum-norm least-squares correction** (QR with column pivoting) to absorb
the redundancy, then advanced along the family by a **standard
predictor-corrector continuation algorithm [ref 2 = Allgower & Georg 1990]**
in the free energy parameter `h`. For long integration times `T^u, T^s` (their
own §3.2 has cases this large), the paper explicitly provides a **multiple-
shooting augmentation** of the same system — extra intermediate points `zi`
along the p.o. and each manifold branch, matching equations between
consecutive points, with an *adaptive* strategy that recomputes the number of
points and sub-interval lengths at each continuation step to bound
`||Dφ_t(zi)||_∞ < M` (typically tens to hundreds).

### Cross-check: does Casoliva's own Eq. (20) match this algorithm?

**Confirmed, term-for-term — Casoliva's Eq. (20) (2010 JGCD paper, p.1635) is
a direct, verbatim transcription of this paper's system**, with identical
unknowns `(h, T, Y, λ^u, V^u, λ^s, V^s, θ^u, T^u, θ^s, T^s)` and identical
equation structure (`H(Y)-h=0`, `g1(Y)=0`, `φ_T(Y)-Y=0`, the two eigenvector
conditions, two unit-norm conditions, two `g2` matching conditions, one
manifold-matching condition) — the paper even reuses the exact scalar
`ξ0=10^-6` convention. Casoliva's own text explicitly credits "a standard
predictor-corrector algorithm [ref 29]" applied to this system, matching
BMO's own description precisely.

**One detail Casoliva's paper (and, following it, the `#725` digest) leaves
out entirely: the multiple-shooting augmentation.** Neither the 2010 JGCD
paper nor the 2008 AIAA precursor mentions multiple shooting anywhere
(confirmed by grep of both extracted texts — zero hits for "multiple
shoot*"). This is a real omission worth flagging: BMO's own paper explicitly
recommends multiple shooting *because* integration times can grow large, and
Casoliva's own He1 family reaches a documented connection flight time of
**113.6319 days** (2010 Table 5) — a regime BMO's own paper (§3.2) uses as its
own worked example for exactly this augmentation. Casoliva's Eq. (20) as
published is presented only in single-shooting form; whether Casoliva's own
implementation silently used multiple shooting under the hood (unstated) or
simply tolerated worse conditioning at long `T^u/T^s` is not discoverable
from the Casoliva papers' text. **Correction/addition to the `#725` digest:**
its characterization "predictor-corrector on a Hamiltonian/monodromy-
eigenvector system (2010 Eq. 20)" is accurate as far as it goes, but should
additionally note that BMO's own paper provides — and explicitly motivates for
exactly this long-TOF regime — a multiple-shooting extension that Casoliva's
own text never mentions using or needing.

### Mandatory citation-mining pass (29 references, full read)

**Already in corpus / already flagged elsewhere in `#730`:** ref 7 (Casoliva
2008, in corpus); ref 23 (McGehee 1969, `#730` §4 item 27); ref 21 (Llibre,
Martínez & Simó 1985, §4 item 28); ref 18 (Howell, Barden, Wilson & Lo 1998,
flagged in `#725`); ref 26 (Parker & Lo 2006, flagged in `#725`); ref 6
(Canalias & Masdemont 2006, flagged in `#725`); ref 10 (Dunham & Farquhar
2003, `#730` §9); ref 20 (Koon, Lo, Marsden & Ross 2000, in corpus); ref 14
(Gómez, Koon, Lo, Marsden, Masdemont & Ross 2004, in corpus); ref 27
(Szebehely 1967, in corpus, digested); refs 28/29 (Wilczak & Zgliczyński 2003/
2005/2006 heteroclinic-connections Parts I/II, in corpus, digested/
KNOWN_CORPUS).

**Genuinely new candidates (flagged, NOT acquired):**
- *High:* Gómez & Mondelo 2001, "The dynamics around the collinear
  equilibrium points of the RTBP," *Physica D* 157:283-321 [ref 17] —
  foundational L1/L2/L3-dynamics reference this whole homoclinic-continuation
  method presupposes; Gidea & Masdemont 2007, "Geometry of homoclinic
  connections in a planar circular restricted three-body problem," *IJBC*
  17:1151-1169 [ref 13] — directly on-point complementary geometric treatment
  of the same object class.
- *Medium:* Conley 1968, "Low energy transit orbits in the restricted
  three-body problem," *SIAM J. Appl. Math.* 16:732-746 [ref 8] — foundational
  classical reference for the whole low-energy-transit/homoclinic-tube
  lineage; Doedel, Kooi, Van Voorn & Kuznetsov 2009, "Continuation of
  connecting orbits in 3D-ODEs: (II) Cycle-to-cycle connections," *IJBC*
  19:159-169 [ref 11] — an AUTO-based alternative continuation approach to the
  same connecting-orbit problem, explicitly contrasted by BMO's own intro;
  Gómez & Masdemont 2000, "Some zero cost transfers between Libration Point
  orbits," AAS 105:1199-1216 [ref 16]; Gómez, Marcote & Mondelo 2005,
  "The invariant manifold structure of the spatial Hill's problem," *Dyn.
  Syst.* 20:115-147 [ref 15]; Mondelo 2001 PhD thesis, "Contribution to the
  study of Fourier methods for quasi-periodic functions..." [ref 25] (Mondelo
  is a Casoliva/BMO co-author — foundational to the Fourier/quasi-periodic
  machinery this whole lineage relies on); Meyer & Hall 1992, *Introduction to
  Hamiltonian Dynamical Systems and the N-body Problem* textbook [ref 24].
- *Low-medium:* Barrabés & Mikkola 2005, "Families of periodic horseshoe
  orbits in the RTBP," *A&A* 432:1115-1129 [ref 4]; Barrabés & Ollé 2006,
  "Invariant manifolds of L3 and horseshoe motion in the RTBP," *Nonlinearity*
  19:2065-2089 [ref 5] — same-author horseshoe-orbit tangent (not cycler-
  relevant, but same continuation-method lineage); Llibre & Ollé 2001, Saturn
  coorbital satellites, *A&A* 378:1087-1099 [ref 22]; Connors et al. 2002,
  horseshoe-orbit asteroid discovery [ref 9] — tangential, not flagged for
  priority.
- Tool/software references (TIDES Taylor integrator [ref 1], LAPACK guide
  [ref 3], AUTO97 [ref 12], Jorba & Zou 2009 Taylor-integrator generator
  [ref 19]) not flagged — infrastructure, not domain results.

---

## 3. Barrabés & Gómez 2002 — "Spatial p-q Resonant Orbits of the RTBP"

*CMDA* 84(4):387-407. `#730` §4 item 25 (HIGH — "the exact analytic in/out-map
matched-asymptotics seed-generation method Casoliva reproduces verbatim").

**Method:** derives, to first order in `µ^α` (a small mass-ratio/close-
approach-radius scaling parameter), analytic **in-map** and **out-map**
expressions for the position/velocity of a spacecraft passing on a sphere of
radius `µ^α` around the smaller primary M, for orbits making `p` spacecraft
revolutions per `q` primary revolutions between consecutive close approaches
(the classical "second-species," Poincaré/Hénon terminology). The out-map
(§4, Eqs. 31-35) propagates a state on the sphere boundary forward around a
p-q resonant Keplerian arc back to the sphere; matching in- and out-maps at
the boundary gives the periodicity/matching condition (Theorem 3.1). The
paper is explicitly **spatial/3D** throughout (states include `z, ż`; angles
`φ, ψ` parameterize out-of-plane inclination/orientation) but **explicitly
treats the planar case as a documented sub-case** (`φ=φ0=0`, discussed
directly in the text around Fig. 1 and again in the numerical examples of
§3), and its own Conclusions state a "forthcoming paper" (= `#730` §4 item 26,
Barrabés & Gómez 2003, "Three-Dimensional p-q Resonant Orbits Close to
Second Species Solutions" — still not acquired) will apply the spatial
formulas to actual second-species-orbit solutions.

### Cross-check: does Casoliva's Eq. (14-18) match this paper's method?

**Confirmed.** Casoliva's own text (2010 JGCD, p.1628) states the seed
formulas "can be written as [refs 9,10]" — citing this 2002 paper (ref 9) and
its 2003 companion (ref 10) jointly — and Casoliva's planar Eq. (14)
(`x_i, y_i, u_i, v_i` on a circle of radius `ε` around the Moon) is exactly
the `φ0=ψ0=0` planar restriction of this paper's general spatial out-map
construction; Casoliva's Eqs. (16)/(17) (the `CJ`/`ε` matching constraint and
the nonlinear constraint on `ε`) and Eq. (18) (total flight time `T ≈ 2πq +
O(ε)`) are the same functional forms this paper derives (its own Theorem 3.1
and the `t2` return-time construction in §4). This is a genuine, verified
match, not merely a plausible attribution — the **planar special case Casoliva
actually uses is explicitly present in this exact 2002 paper**, not solely in
the still-unacquired 2003 companion. **No discrepancy or oversimplification
found** in Casoliva's reproduction; the `#725` digest's characterization
("analytic in/out-map matched-asymptotics") is accurate.

One quantitative caveat for any future reproduction attempt: this paper's
matched-asymptotic expansion is explicitly an `O(µ^α)`/`O(µ^{1-α})`
first-order approximation, valid only for very small `ε` (Casoliva's own text
notes "orbits will be periodic with an error of order... for `α ∈
[1/3,1/2]`") — the seed must still be differentially corrected to the exact
periodic orbit at `µ_EM`, exactly as both Casoliva papers do; this paper's
formulas are never claimed to be exact even at moderate `ε`.

### Mandatory citation-mining pass (5 references, full read — a very short
reference list for this paper)

- Barrabés 2001 PhD thesis, "Òrbites de segona espècie del problema espacial
  de 3 cossos" (Universitat Autònoma de Barcelona) — likely superseded/
  contained by the published 2002/2003 CMDA papers; low priority, flagged
  for completeness only.
- Font, Nunes & Simó 2001, "Successive quasi-collisions in the planar
  circular RTBP," *Nonlinearity* (to appear) — tangential (quasi-collision
  geometry, not resonant-orbit seed generation); low priority.
- Perko 1996, *Differential Equations and Dynamical Systems* textbook, and
  Szebehely 1967 (already in corpus) — general background, not flagged/
  no action.
- Yen 1985, Mercury orbiter gravity-assist mission design, AAS 85-346 —
  out of scope (not RTBP-cycler-relevant), not flagged.

No genuinely new high-priority gap surfaced by this paper's own citation
list beyond the already-flagged 2003 companion paper (§4 item 26).

---

## 4. Franz & Russell 2022 — "Database of Planar and Three-Dimensional Periodic Orbits and Families Near the Moon"

*J. Astronaut. Sci.* 69(6):1573-1612. `#730` §6 item 44 (HIGH — "a large
(13-million-solution), directly comparable lunar-periodic-orbit census").

**Method:** a broad grid search over initial conditions on the x-axis / x-z
plane (planar symmetric, x-z-planar-symmetric, and doubly-symmetric orbit
classes, up to 32-64 plane crossings) in the Earth-Moon CRTBP at
`µ=1.215058392535863e-2`, producing **over 13 million** periodic-orbit
solutions. DBSCAN (density-based unsupervised clustering) groups the raw grid
hits into **33,980 family/sub-family clusters**, sorted into smooth ordered
curves with a custom cluster-confidence measure (>80% of 3D data, >62% of
planar data clustered with high confidence). Each orbit is tagged with
stability, perilune distance, revolution count/center, etc. The full database
+ code + interactive GUI is openly archived at Zenodo (DOI
`10.5281/zenodo.6411980`, "latest version" pointer `10.5281/zenodo.6411979`).

### Cross-check #1: the Broucke lineage claim

**Confirmed verbatim.** The abstract's own final content sentence reads (in
full): *"The resulting database is an extension of several recent Lunar
periodic orbit studies, and can be considered a modern update to Broucke's
seminal database of planar cislunar periodic orbits."* Ref [42] is Broucke,
R., "Periodic orbits in the restricted three-body problem with Earth-Moon
masses," JPL TR 32-1168 (1968) — exactly `#730` §4 item 21 (still the
highest-corroborated non-Kumar-lineage gap in the whole backlog, now
independently re-flagged a **4th** time by this reading, following `#725`,
the 2026-06-11 RRT mining note, and `#728`-kumar-moreno-networks).

### Cross-check #2: overlap with / cross-check value for this project's own Earth-Moon rows

**No direct overlap in scope, but a real (unexecuted) cross-check
opportunity exists for a specific subset of rows.** This database's own
search is restricted to **symmetric** periodic orbits only (two perpendicular
x-axis or x-z-plane crossings define the construction throughout — confirmed
by grep of the methodology section, §3); it never mentions cyclers,
homoclinic/heteroclinic connections, or resonant p-q construction by name
anywhere in the text (grep for "cycler," "homoclinic," "Casoliva," "Barrabés,"
"Ross," "Roberts," "Tsoukkas" all return zero hits outside the reference
list). This has two consequences:

- **No bearing on Casoliva's Class 2** (He1-4/Hm1-2 homoclinic-shadowing
  orbits) or on the Ross-Roberts-Tsoukkas `(k1,k2)` manifold-tube cyclers,
  neither of which is guaranteed (or, in Class 2's asymmetric-member cases,
  even eligible) to be symmetric about the x-axis.
- **A genuine, not-yet-executed cross-check candidate for Casoliva's Class 1**
  (the p-q resonant cyclers), specifically its symmetric members — Casoliva's
  own Table 3 states "All are symmetric about the x-axis except 1-2e, 7-3b,
  7-3c" (i.e., 1-2c, 1-2d, 2-1a, 2-1b, 3-2c, 7-3a **are** symmetric). Those
  six rows' Jacobi constants (`C_J` = 1.569, 2.580, 0.489, 1.196, 0.709,
  1.022 respectively, Casoliva Table 3) are, in principle, independently
  queryable against this paper's own public Zenodo database as an external
  positive-control cross-check — this was **not executed in this pass**
  (out of this digest's scope; flagged as a concrete future task).
- **Mass-ratio caveat for any such future cross-check:** this paper's own
  `µ=1.215058392535863e-2` differs from Casoliva's `µ_EM=0.0121529529` at the
  ~0.02% relative level (5th significant digit) — small, but non-zero;
  a literal `C_J` match would need either re-deriving Casoliva's orbits at
  Franz-Russell's `µ` or vice versa, not a naive direct lookup.

### Mandatory citation-mining pass (75 references; full intro/background +
full reference list read, filtered to domain-overlapping items — the
majority of this paper's reference list is cislunar-navigation/SSA/ML-
clustering background outside this project's search domain, consistent with
`#730`'s own §9 triage convention for background-heavy citation lists)

**Already in corpus / already flagged elsewhere in `#730`:** ref 42 (Broucke
1968, §4 item 21, see above); ref 36 (Restrepo & Russell 2018, `#730` §1 false
gap, already in corpus digested); ref 44 (Anderson, Campagnola & Lantoine,
CMDA 124:177-199, `#730` §3 item 18 — note this paper's own reference gives
the year as **2015**, while `#730` §2/§3 lists it as **2016**; both point to
the same CMDA 124 paper, a minor citation-year discrepancy worth noting, not
correcting, since journal "received"/"published" years commonly differ by a
year for this venue); ref 26 (Lantoine, Russell & Campagnola 2011, Acta
Astronaut. 68:1361-1378 — **independently re-flagged, matches the new
candidate surfaced from Anderson 2021's own ref 24 above; now a 2-flag
recurrence within this single task**); ref 31 (Lara, Russell & Villac 2007,
already flagged low-priority in the `#725` Casoliva digest); ref 37
(Szebehely 1967, in corpus).

**Genuinely new candidates, topically overlapping (flagged, NOT acquired):**
- *High:* Russell, R.P. 2006, "Global search for planar and three-dimensional
  periodic orbits near Europa," *J. Astronaut. Sci.* 54(2):199-226 [ref 35] —
  same senior author, the direct methodological predecessor (a Europa-system
  analogue of exactly this paper's own grid-search approach).
- *Medium-high:* McArdle & Russell 2021, "Circulating, eccentric periodic
  orbits at the Moon," *CMDA* 133(4):18 [ref 74] — same lineage, recent,
  lunar-specific; Vaquero & Howell 2014, "Design of transfer trajectories
  between resonant orbits in the Earth-Moon restricted problem," *Acta
  Astronaut.* 94:302-317 [ref 25] — directly cycler-adjacent (resonant-orbit
  transfer design), same author (Vaquero) as the already-flagged Vaquero
  Escribano 2013 PhD thesis in `#730` §3 item 15's list, corroborating that
  gap; Tsirogiannis, Perdios & Markellos 2008, "Improved grid search method:
  an efficient tool for global computation of periodic orbits," *CMDA*
  103:49-78 [ref 43] — the direct methodological predecessor to this paper's
  own grid-search algorithm.
- *Medium:* Restrepo & Russell 2017, "Patched periodic orbits: a systematic
  strategy for low energy transfer design," AAS/AIAA conference [ref 22] —
  **a different paper from the already-flagged/likely-duplicate Restrepo &
  Russell 2017 "database of planar axisymmetric periodic orbits" AAS-17-694**
  (`#730` §1); this "patched periodic orbits" title is distinct and not
  previously flagged; Koon, Lo, Marsden & Ross 2001, "Resonance and capture of
  jupiter comets," *CMDA* 81:27-38 [ref 32] — distinct from the already-corpus
  Koon et al. 2000 Chaos paper; Folta, Bosanac, Guzzetti & Howell 2015, "An
  Earth-Moon system trajectory design reference catalog," *Acta Astronaut.*
  110:341-353 [ref 27], and its companion Guzzetti, Bosanac, Haapala, Howell &
  Folta 2016, "Rapid trajectory design... via an interactive catalog," *Acta
  Astronaut.* 126:439-455 [ref 28] — directly relevant catalog-methodology
  papers from the same Purdue lineage already represented elsewhere in this
  project's corpus (Howell/Barden/Wilson/Lo GENESIS work).
- *Low-medium:* Hénon 1969, "Numerical exploration of the restricted problem
  V. Hill's case," *A&A* 1:223-238 [ref 38], and Hénon 2003, "New families of
  periodic orbits in Hill's problem," *CMDA* 85:223-246 [ref 39] — distinct
  from the already-flagged Hénon 1997 book (`#730` §4 item 29), same author's
  Hill's-problem-specific papers; Howell, Marchand & Lo 2001, "Temporary
  satellite capture of short-period jupiter family comets," *JAS* 49:539-557
  [ref 33] — tangential dynamical-systems-capture parallel.

None of the ML/clustering-algorithm references (DBSCAN, HDBSCAN, fuzzy
clustering, etc., refs 69-73) or the cislunar-navigation/SSA background
references (refs 1-21) are flagged — infrastructure/background outside this
project's search-method domain, consistent with prior digests' treatment of
similar reference clusters.

---

## Summary answers (for the dispatching session)

- **Paper 1 (Anderson 2021) vs `#727`'s benchmark claim:** ΔV confirmed
  precisely (146.71 m/s CRTBP / 135.04 m/s COSMIC, both genuinely "a little
  less than 150 m/s"). TOF "~40 days" is Kumar et al.'s own reasonable
  paraphrase of Table 2's burn-time schedule plus an untabulated final coast
  segment — not a number printed verbatim anywhere in Anderson 2021's own
  text. No correction needed to `#727`'s digest.
- **Paper 2 (Barrabés-Mondelo-Ollé 2009) vs Casoliva's Eq. 20:** confirmed
  term-for-term match. One addition: BMO's own paper provides a
  multiple-shooting extension explicitly motivated for the long-integration-
  time regime Casoliva's own He1 family (113.6-day connection flight times)
  falls squarely into, but neither Casoliva's papers nor the `#725` digest
  mention it — worth noting as an omission, not a factual error.
- **Paper 3 (Barrabés-Gómez 2002) vs Casoliva's Eqs. 14-18:** confirmed —
  Casoliva's planar seed formulas are exactly the documented planar (`φ0=0`)
  special case of this paper's own general spatial in/out-map construction.
  No discrepancy found.
- **Paper 4 (Franz-Russell 2022):** explicitly self-identifies as a modern
  Broucke-1968 update (verbatim abstract quote); its own text never mentions
  cyclers, Casoliva, or RRT — the database's symmetric-orbit-only scope means
  it has no bearing on Casoliva's Class 2 or the RRT `(k1,k2)` cyclers, but
  is a genuine (unexecuted) future cross-check candidate for Casoliva's six
  symmetric Class 1 p-q resonant rows (mass-ratio caveat noted).
- Citation-mining across all four surfaced ~35 new candidates (none
  acquired); the strongest recurring signal is **Lantoine, Russell &
  Campagnola 2011** (independently cited by both Anderson 2021 and
  Franz-Russell 2022 within this single task) and the **4th independent flag**
  on Broucke 1968 (`#730` §4 item 21), now the single most-corroborated
  unacquired gap in the whole backlog.
