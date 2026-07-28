# Digest: Llibre-Martínez-Simó 1985 (L2 Lyapunov transversality) + Doedel-Romanov-Paffenroth-Keller-Dichmann-Galán-Vioque-Vanderbauwhede 2007 (elemental libration-point orbits) (#749)

**Task:** `#749`, continuing the `#730` consolidated acquisition backlog
(`docs/notes/2026-07-27-730-acquisition-backlog-master-list.md`). The user
directly supplied PDFs for `#730` §4 item 28 and §6 item 43 of that list —
both independently verified page-1-exact by the coordinating session before
this task started.

**Filed** (private `cyclers_pdf` repo):
1. `llibre-martinez-simo-1985-transversality-invariant-manifolds-lyapunov-l2-jde-58-104-doi-10.1016-0022-0396(85)90024-5.pdf`
2. `doedel-romanov-paffenroth-keller-dichmann-galan-vioque-vanderbauwhede-2007-elemental-periodic-orbits-libration-points-ijbc-17-8-2625-doi-10.1142-S0218127407018671.pdf`

**OCR status: both text-layer, no OCR needed.** Both are Acrobat/Distiller-
produced (native journal-typeset PDFs); read in full via the `pages` chunked
Read parameter (Llibre: 53pp across 3 chunks; Doedel: 53pp across 4 chunks),
confirming dense, extractable body text throughout (no image-only pages).

---

## 1. Llibre, Martínez & Simó 1985 — "Transversality of the Invariant
Manifolds Associated to the Lyapunov Family of Periodic Orbits near L2 in
the Restricted Three-Body Problem"

*J. Differential Equations* 58:104–156. `#730` §4 item 28 (Medium-high,
flagged by `#725` as "companion theoretical-existence result to McGehee,
specifically for the L2 Lyapunov family transversality question underlying
the Casoliva Class 2 construction"). DOI `10.1016/0022-0396(85)90024-5`
(confirmed, matches this PDF's own header exactly).

**Method:** studies the (planar) RTBP near the Euler point L2 for Jacobi
constant `C` near `C2` (the value at L2). For `C < C2` there is a hyperbolic
periodic orbit family (the classical "family (c)", i.e. the L2 planar
Lyapunov family) whose 1D stable/unstable invariant manifolds are 2D
cylinders. The paper is explicitly analytical + numerical, in three tiers:

- **§2–4 (analytical):** the flow near L2 is linearized via the Moser–
  Lyapunov/Conley theorem (Eq. 2.1–2.2), and McGehee's own qualitative
  picture of the region `L` bounded by two spheres `l1, l2` (rotation of two
  half-circles, Fig. 2.1) is reused directly. **Theorem A** derives an
  explicit asymptotic expansion (as `mu -> 0`) for the unstable branch
  `W^u,S` of `W^u_L2` in terms of two numerically-determined constants
  `N(inf) = 5.1604325...` and `M(inf) = 2.1330587...` (obtained by direct
  numerical integration of the Hill-problem limit ODEs, Eq. 4.5), giving a
  first-intersection-orthogonality condition `mu_k = 1/(N(inf)^3 k^3)(1 +
  o(1))` — i.e. an explicit **countable sequence of mass-ratio values
  `mu_k` at which a genuine (symmetric) homoclinic orbit to L2 exists**.
  **Theorem B** extends this to `Delta C = C2 - C` small but nonzero,
  deriving the shape of the `Gamma^{u,S}_1` intersection curve with `y=0` (a
  spiral, Eq. 4.16–4.21, Fig. 4.9) and a sufficient condition
  `Delta C > L mu_k^{4/3} (mu - mu_k)^2` for **transversal** symmetric
  homoclinic orbits to exist near each `mu_k` — this is the paper's own
  headline transversality result.
- **§5 (numerical):** extensive verification — Table I gives `mu_k`,
  `mu-bar_k` (extrema of the crossing angle `alpha`) for `k=0..23`; Table II
  gives `Delta C(mu)` for the first tangential (non-transversal boundary)
  homoclinic orbit over `mu in [0.00023, 0.0042]`; Figs. 5.1–5.13 show the
  actual computed manifold shapes, including numerically confirmed
  **(46,46)-homoclinic points** in the Jupiter region at `mu=0.00025,
  C=0.0000783` (used as Theorem C's own worked positive control, §5 final
  paragraph).
- **§6 (symbolic dynamics):** proves **Theorem C**, a full Bernoulli-shift-
  on-Z embedding (via Moser's abstract horseshoe machinery, Devaney's
  "singular cross section" formalism for two big squares `U_+, U_-` rather
  than one small square) showing there exist **doubly-infinite sequences of
  quasirandom orbits** shuttling between the Sun and Jupiter regions,
  making an arbitrary prescribed number of revolutions around L2 at each
  pass — the existence-of-chaotic-shadowing consequence of the transversal
  homoclinic orbits established in §4–5. **Theorem C's proof is a direct
  instance of Lemma 6.1 (§6, spiraling of vertical/horizontal strips under
  the Poincaré map `f`)**, itself built on McGehee's spiraling-image
  argument from §2.

**Key numerical constants (reusable as a positive control for any future
L2-Lyapunov-manifold reproduction):** `N(inf) = 5.1604325`, `M(inf) =
2.1330587` (Hill-problem-limit constants, checked by the paper's own
identity `N(inf)^2/3 - M(inf)^2 = C^H(L2) = 3^{4/3}`); empirical fit `mu_k ~
0.00728832 k^-3 + 0.01702374 k^-4 + 0.022906 k^-5 + 0.07213 k^-6 + ...` for
`k=2..13`, validated against `k=21,22,23` to relative error `< 5e-4`.

### Cross-check: does this project's own L2 Lyapunov-manifold/transversality
work already use, contradict, or omit anything this paper establishes?

**No contradiction; a genuine, documented capability gap.** Grepped
`src/cyclerfinder` broadly (Lyapunov, transversal, homoclinic, McGehee) and
inspected the two most relevant modules directly:

- `src/cyclerfinder/genome/transit_manifold.py` (`#547`) is this project's
  own from-scratch, textbook-validated positive control for the
  Conley-McGehee/KLMR "tube dynamics" picture — but it targets **transit
  classification** (does the unstable-manifold branch cross into the
  secondary's realm), not **homoclinic transversality**. Its own docstring
  cites Conley 1968, McGehee 1969, and KLMR 2006 (already in corpus,
  `#547`'s own dependency chain), but never this Llibre-Martínez-Simó 1985
  paper — confirmed by direct grep of the file, zero hits for "Llibre" or
  "transversal[a-z]* homoclinic".
- No module anywhere in `src/cyclerfinder` computes an L2-Lyapunov-family
  **homoclinic** orbit (an orbit that returns to `L2`'s own manifold, as
  opposed to a transit/heteroclinic connection between two different
  libration points or a torus). Grepping for "homoclinic" project-wide
  (`crnbp_real_ephemeris_consistency.py`, `tisserand.py`,
  `variational_qbcp_arc.py`, `heteroclinic_cycle.py`,
  `ccr4bp_heteroclinic_search.py`, `ccr4bp_real_ephemeris_consistency.py`,
  `qp_tori_arclength.py`, `data/validate.py`) shows every one of those
  hits is either a **heteroclinic** connection (between two distinct
  periodic orbits/tori) or a validation/consistency check, never a
  self-return homoclinic orbit to a single Lyapunov family the way this
  paper constructs. **This is the same capability-gap class `#742`'s
  Barrabés-Mondelo-Ollé cross-check already documented for the *general*
  (non-L2-specific) homoclinic-continuation machinery** — that gap is now
  independently confirmed from the L2-specific analytical side as well: no
  code anywhere in this project reproduces or contradicts this paper's own
  `mu_k` sequence, its asymptotic manifold-shape formulas (Eqs. 4.16–4.21),
  or its Theorem C symbolic-dynamics shift.
- The paper's own §6 symbolic-dynamics machinery (Bernoulli shift on `Z`,
  Devaney's singular-cross-section formalism) has no counterpart anywhere
  in this codebase — grepped for "symbolic dynamics", "Bernoulli shift",
  "Moser" (the horseshoe-embedding sense, not an author-name collision);
  the only "symbolic dynamics" hit project-wide is a comment in
  `genome/heteroclinic_cycle.py` referencing the *unrelated*
  Wilczak-Zgliczyński L1<->L2 heteroclinic-cycle symbolic dynamics (a
  structurally different, already-published-and-digested result, not this
  paper's L2-homoclinic shift).
- **Direct relevance to `#725`'s own Casoliva Class 2 (He1-4/Hm1-2)
  homoclinic-shadowing cyclers**, as the master list's own item-28 note
  already anticipated: this paper's Theorem A/B machinery is the specific
  L2-Lyapunov-family transversality existence proof underlying the
  *feasibility* of Casoliva's He-family construction (an L1-Lyapunov-orbit
  homoclinic connection is the geometric cousin of this paper's
  L2-Lyapunov-orbit homoclinic connection — same phenomenon, adjacent
  libration point). Casoliva's own papers (already acquired, `#725`) never
  cite this paper directly by section/theorem number (only in the
  reference list, per `#742`'s BMO citation-mining pass, ref 21) — so this
  reading adds no correction to `#725`'s digest, only confirms the
  theoretical-lineage flag `#725` already carried was substantively
  correct: this paper *is* exactly the transversality-existence machinery
  the He-family construction implicitly relies on, now read directly
  rather than inferred from title alone.

**No catalogue row cites this paper or its specific `mu_k`/`N(inf)`/`M(inf)`
values** — grepped `data/catalogue.yaml` for "Llibre" and for the specific
constants; zero hits. Confirmed genuinely un-adopted, not contradicted.

### Mandatory citation-mining pass (15 references, full read)

**Already in corpus / already flagged elsewhere in `#730`:**
- McGehee, R., "Some Homoclinic Orbits for the Restricted Three-Body
  Problem," Ph.D. thesis, Univ. of Wisconsin (1969) [ref 9] — `#730` §4 item
  27, **confirmed genuinely unacquirable-free** (ProQuest-only, re-checked
  `#743`) — this is now a further independent citation confirming its
  centrality (the present paper's own Introduction: "Through this paper free
  use is made of several important results contained in [McGehee's] Ph.D.
  thesis").
- Szebehely, V., *Theory of Orbits*, Academic Press (1967) [ref 15] —
  already in corpus, digested.

**Genuinely new candidates (flagged, NOT acquired):**
- *High:* Llibre, J. & Simó, C., "Oscillatory Solutions in the Planar
  Restricted Three-Body Problem," *Math. Ann.* 248 (1980), 153–184 [ref 6]
  — the paper's own Introduction states explicitly "The present paper is an
  extension of a previous work that can be found in [10]" where [10] is
  this exact paper — i.e. this is the **direct predecessor this paper
  extends**, same first author, same L2-neighborhood subject matter.
  Llibre, J. & Simó, C., "Some Homoclinic Phenomena in the Three-Body
  Problem," *J. Differential Equations* 37 (1980), 444–465 [ref 7] — the
  companion first-species-homoclinic predecessor by the same authors,
  explicitly the direct forerunner of this paper's own homoclinic-orbit
  question (grepped `#730`/`CORPUS_INDEX.md`, confirmed absent — neither is
  flagged anywhere in this project's prior citation-mining history despite
  being the single strongest same-author lineage candidates in this
  reference list).
- *Medium-high:* Gómez, G., Llibre, J., Martínez, R. & Simó, C., *Dynamics
  and Mission Design Near Libration Points, Vol. I: Fundamentals, The Case
  of Collinear Libration Points* and *Vol. II: ...Triangular Libration
  Points*, World Scientific (2001a, 2001b) — **not in this paper's own
  reference list** (this is a 1985 paper, predates the 2001 books) but
  surfaces independently below as Doedel 2007's refs — cross-referenced
  here for continuity since the same co-authors (Llibre, Martínez, Simó)
  wrote both; see the combined flag under paper 2's citation-mining pass
  below to avoid double-counting.
- *Medium:* Devaney, R. L., "Singularities in Classical Mechanical Systems,"
  in *Ergodic Theory and Dynamical Systems* (A. Katok, ed.), Birkhäuser
  (1981), 211–333 [ref 5] — the specific "singular cross section"
  formalism this paper's own §6 explicitly builds its two-big-square
  Poincaré-map argument on ("a kind of singular Poincaré map as it has been
  considered by Devaney"); Moser, J., *Stable and Random Motions in
  Dynamical Systems*, Princeton Univ. Press (1973) [ref 12] — the abstract
  horseshoe/symbolic-dynamics theorem (Proposition 6.1) this paper's
  Theorem 6.1 directly cites and extends to the doubly-infinite-alphabet
  case.
- *Low-medium:* Alekseev, V. M., "Quasirandom Dynamical Systems, I, II,
  III," *Math. USSR-Sb.* 5–7 (1968–1969) [ref 1] — foundational precursor
  to the quasirandom-motion existence framing the paper's own §1
  explicitly credits ("first proved by Sitnikov... generalization... given
  by Alekseev"); Sitnikov, K. A., "The Existence of Oscillatory Motions in
  the Three-Body Problem," *Soviet Phys. Dokl.* 5 (1961), 647–650 [ref 13]
  — the historical origin-point of the whole oscillatory/quasirandom-motion
  line this paper extends to the RTBP proper; Smale, S., "Diffeomorphisms
  with Many Periodic Points," in *Differential and Combinatorial Topology*
  (S. S. Cairns, ed.), Princeton Univ. Press (1965), 63–80 [ref 14] — the
  classical horseshoe-map source both Moser's and this paper's own
  symbolic-dynamics arguments ultimately trace to.
- *Low:* Churchill, R. C., Pecelli, G. & Rod, D. L., "A Survey of the
  Hénon-Heiles Hamiltonian with Applications to Related Examples," in
  *Lecture Notes in Physics* Vol. 93, Springer (1978), 76–136 [ref 2] —
  cited only as a worked comparison example (Hénon-Heiles) for Theorem 6.2's
  generalization, tangential to this project's cycler-search domain; Conley,
  C., "On the Ultimate Behavior of Orbits with Respect to an Unstable
  Critical Point," *J. Differential Equations* 5 (1969), 136–158, and
  Conley, C., "Twist Mappings, Analyticity and Periodic Solutions Which
  Pass Close to an Unstable Periodic Orbit," in *Topological Dynamics* (J.
  Auslander, ed.), Benjamin (1968), 129–154 [refs 3, 4] — foundational
  Conley-theory background, general dynamical-systems machinery already
  represented by the Conley-McGehee lineage `transit_manifold.py` already
  cites; Lukjanov, L. G., "A Study of Asymptotic Solutions in the Vicinity
  of the Collinear Libration Points of the Restricted Three-Body Problem,"
  *Celestial Mech.* 15 (1977), 489–500 [ref 8] — the source of the local
  L2-manifold expansion formulas this paper's own §4 opening reuses ("Local
  expressions for these manifolds were given by Lukjanov"), a narrow
  methods reference; Martínez, R., "Moviments Quasi-Aleatoris en el
  Problema Restringit, Circular i Pla de 3 Cossos," Ph.D. thesis, Univ.
  Autònoma de Barcelona (1981) [ref 10] — the second author's own thesis,
  likely subsumed content-wise by this published paper.

No genuinely surprising high-priority gap beyond the two Llibre & Simó 1980
predecessor papers, now flagged for the first time in this project's
citation-mining history.

---

## 2. Doedel, Romanov, Paffenroth, Keller, Dichmann, Galán-Vioque &
Vanderbauwhede 2007 — "Elemental Periodic Orbits Associated with the
Libration Points in the Circular Restricted 3-Body Problem"

*Int. J. Bifurcation and Chaos* 17(8):2625–2677. `#730` §6 item 43 (HIGH —
flagged by `#728`-kumar-moreno-networks as "a 'stepping stone' predecessor
explicitly cited"). DOI `10.1142/S0218127407018671` (CrossRef-confirmed
2026-07-28, `#743`/`#744`; independently re-confirmed page-1-exact by this
task). **Distinct from the already-corpused Doedel-Keller-Kernevez 1991
Parts I/II** (AUTO/numerical-bifurcation-analysis foundations papers) and
the Doedel-Paffenroth-Keller 2003 conservative-3-body-systems paper — all
three are cited HERE as this paper's own methodological ancestry, not
duplicates (grepped `CORPUS_INDEX.md` first per project convention:
confirmed distinct titles/topics, all three already correctly filed
separately).

**Method:** a broad, largely-visual (53pp, ~60 figures) AUTO/boundary-value-
continuation survey of all "elemental" periodic-orbit families that
bifurcate from the five CR3BP libration points L1–L5, plus their secondary
bifurcating families, computed across the FULL mass-ratio range `mu in
(0, 1/2]` (not a single system). Continuation uses AUTO's orthogonal
collocation with an *unfolding parameter* `lambda` (Appendix A.4, Eq. A.6)
— since the CR3BP is conservative (Jacobi integral), a naive periodic-
boundary-value system is singular; `lambda` is solved for at every
continuation step but stays zero to numerical precision at convergence,
exactly analogous to a Lagrange-multiplier regularization. **Branch points**
(codimension-1 transcritical/pitchfork bifurcations where two distinct
families intersect with identical orbit+period) are tracked by a *fully
extended system* (Appendix B.4–B.5: doubling the ODE count via an adjoint
operator, `2n+2` extra equations) that can itself be continued in `mu` as a
**locus of branch points** (Figs. 28–39) — this is the paper's own
central technical contribution beyond a plain family sweep.

**Taxonomy (Table 1, Fig. 3 for Earth-Moon `mu=0.01215`):** Planar Lyapunov
(`L1,L2,L3`), Long/Short-Period from L4/L5 (`L4,L5,S3`), Vertical
(`V1..V5`), Halo (`H1,H2,H3`, bifurcating from `L1,L2,L3` at branch points
`L11,L21,L31`), Axial (`A1,A2,A3`), Backflip (`B1,B2,B3`, non-planar,
bifurcating from Vertical pitchfork points), `W4/W5` (connecting `V4/V5` to
`H1`), planar Circular (`C1,C2`), `D1`/`E1` (reached via Backflip `B3`),
`R2` (bifurcating from `L24` on `L2`), and (for `mu` near 1/2 specifically)
the `S1,S2,T1,T2,K1-K4,X4,X5` families that only exist in a narrow window
near the resonant/singular critical value `mu_12 ~= 0.399`. Table 2 gives 14
critical `mu`-values (folds, collisions, the Hamiltonian-Hopf bifurcation at
`mu_2 = 1/2 - sqrt(23/108) ~= 0.0385`, and the codimension-2 singular point
`mu_12 ~= 0.399` where loci V45 and V32 intersect). §5 computes loci of
**homoclinic orbits** (`S1-infinity, S2-infinity, S3-infinity`) that
terminate the `S1/S2/S3` families, each homoclinic to L4 or L5, using
AUTO's HOMCONT algorithms. §6 treats the equal-mass `mu=1/2` Copenhagen
problem as a structurally distinct special case (extra symmetry Eq. 8; the
Vertical family V1 becomes perfectly vertical, `x=y=0`, and no longer
connects directly to the planar Circular family — a genuine topological
change, not a numerical artifact).

**No numeric IC/period tables are printed anywhere in the paper's own
body text** — every family is documented purely via bifurcation-diagram
plots (energy/period vs. `mu`) and 3D orbit-shape figures; this is
explicitly a "map outlining a collection of phenomena," not a lookup table
(the paper's own §3 states this directly). Any future reproduction attempt
would need either AUTO itself, the cited `mu=0.01215`/`mu=1/2` example
runs, or digitization of the bifurcation-diagram figures — no digit-grade
positive control is directly extractable from this reading alone.

### Cross-check: does this project's own AUTO-adjacent continuation code,
or its own catalogue of libration-point periodic-orbit families, already
have/miss any of the "elemental" family results this survey documents?

**Partial capability overlap, structurally different tooling; the specific
mu-sweep branch-point-locus and equal-mass-case results are novel to this
project.**

- **Continuation machinery:** this project has its own from-scratch
  pseudo-arclength continuation stack (`search/pseudo_arclength.py`,
  explicitly a generic co-dimension-1 primitive factored out of FOUR
  specialized continuations already in the codebase —
  `cr3bp_jacobi_arclength.py`, `mu_continuation.py`, `qp_tori_arclength.py`,
  `narc_continuation.py`), which is the same *pseudo-arclength predict-
  correct-tangent* family Doedel's own Appendix A.2 describes (Eq. A.2,
  Keller 1977) — **the same general method, independently implemented**,
  not derived from or citing AUTO. This project's `bifurcation_detector.py`
  (`#266`/`#347`) detects period-multiplying and saddle-center (`lambda=+1`)
  bifurcations by a **bracket** on adjacent Floquet-multiplier crossings —
  materially SIMPLER than Doedel's own fully-extended branch-point-locus
  system (Appendix B.3–B.5, which doubles the ODE system via an adjoint
  operator to CONTINUE a branch-point locus directly in `mu`, rather than
  bracket individual crossings along one fixed-`mu` family). **This project
  has no code that continues a LOCUS of branch points across `mu` the way
  Doedel's Figs. 28-39 do** — confirmed by grep of `mu_continuation.py`
  and `bifurcation_detector.py`, neither tracks a branch-point locus as a
  first-class continued curve, only individual-family bifurcation brackets
  at fixed `mu` or a single held family continued in `mu`.
  `mu_continuation.py` itself is the closest structural analogue (a
  pseudo-arclength `mu`-sweep of one specific family, the Ross-Tsoukkas
  `(k1,k2)` stable symmetric cycler) but targets an entirely different
  ORBIT CLASS (a stable multi-orbiter cycler family, not any of Doedel's
  L1-L5 elemental families) and does not track branch points at all.
- **Halo/family taxonomy:** `search/halo_family_at_jacobi.py` and
  `search/jpl_family_census.py`/`verify/jpl_periodic_orbits.py` already give
  this project access to L1/L2 halo, vertical, axial, and Lyapunov family
  members — but sourced from **JPL SSD's own periodic-orbit catalog**
  (queried live), not from this paper's own AUTO computation, and JPL's
  catalog does not carry Doedel's own branch-point/critical-`mu` metadata.
  `howell-1984-...pdf` (already in corpus, digested `2026-06-25`) is the
  canonical halo-family survey Doedel's own paper cites [Howell, 1984] as
  its own halo-family predecessor ("the current work also extends that of
  Howell [1984], who mapped out the families of Halo orbits associated
  with the L2 and L3 libration points") — **confirms this project's own
  Howell-1984 corpus entry is exactly the paper Doedel's own text credits**,
  no discrepancy.
- **No catalogue row anywhere references Doedel's own family names** (`A1`,
  `B1`, `V4/V5`, `D1`, `E1`, `R2`, `S1/S2/S3`, `K1-K4`, `X4/X5`, `T1/T2`, or
  any critical-`mu` value from Table 2) — grepped `data/catalogue.yaml` for
  "backflip", "axial", "halo" and confirmed zero rows trace to this paper.
  This project's existing halo/axial/vertical usage is all independently
  sourced (JPL catalog, Howell 1984, or this project's own correctors),
  never from this specific 2007 survey.
- **Kumar & Moreno 2025 lineage (`#728`) confirmed, not contradicted:**
  the master-list's own flag (item 43, "the L1-L5 libration-point-orbit
  bifurcation-network paper the already-acquired Kumar & Moreno 2025 paper
  positions itself as a sequel to") is directly verified by this reading —
  Doedel 2007's own Introduction explicitly frames its work as extending
  Doedel et al. 2003b (already flagged/digested), and the already-acquired
  `#728` Kumar-Moreno digest independently cross-confirmed the H1/H2
  Broucke-junction and L1/L2-Halo-double-cover bridging structure this
  paper's own bifurcation diagram (Fig. 3, L11/L21/L31 -> H1/H2/H3 branch
  points) visually documents for `mu=0.01215`. No numeric conflict found —
  this reading did not extract digit-grade branch-point values to compare
  against Kumar & Moreno's own network-graph nodes (out of this pass's
  scope, a concrete future cross-check: Doedel's own critical-`mu` Table 2
  values vs. Kumar-Moreno's own bifurcation-network critical points, not
  executed here).

### Mandatory citation-mining pass (~62 references, full intro/background +
full reference list read, filtered to domain-overlapping items — the
majority is Doedel's own AUTO/numerical-continuation self-citation lineage,
already represented by the two already-corpused Doedel-Keller-Kernevez 1991
Parts I/II and Doedel-Paffenroth-Keller 2003 papers, consistent with `#730`
§9's tool/methods-reference triage convention)

**Already in corpus / already flagged elsewhere in `#730`:**
- Howell, K. C., "Three-Dimensional, Periodic, 'Halo' Orbits," *Celest.
  Mech.* 32:53–71 (1984) — **already in corpus**, digested
  `2026-06-25-digest-howell-1984-halo-orbits.md`; this paper's own text
  directly credits it as the halo-family predecessor being extended (see
  cross-check above).
- Szebehely, V., *Theory of Orbits: The Restricted Problem of Three Bodies*,
  Academic Press (1967) — already in corpus, digested.
- Doedel, Keller & Kernevez 1991a/1991b — already in corpus, digested
  (this paper's own tutorial-introduction citations).
- Doedel, Paffenroth & Keller 2003b — already in corpus, digested (this
  paper's own direct methodological/results predecessor, explicitly "we
  extend our study of periodic solutions of the CR3BP significantly" from
  it).

**Genuinely new candidates (flagged, NOT acquired):**
- *High:* Gómez, G., Llibre, J., Martínez, R. & Simó, C., *Dynamics and
  Mission Design Near Libration Points, Vol. I: Fundamentals, The Case of
  Collinear Libration Points* and *Vol. II: ...Triangular Libration
  Points*, World Scientific (2001a, 2001b) — the two-volume foundational
  mission-design-near-libration-points reference series, co-authored by
  TWO of this same task's two authors (Llibre, Martínez — the SAME two
  authors as paper 1 above, now a 2-flag same-session recurrence); a
  genuinely central, previously-unflagged gap in this project's own
  libration-point-orbit corpus (grepped `CORPUS_INDEX.md`, confirmed
  absent).
- *Medium-high:* Henrard, J., "The Web of Periodic Orbits at L4," *Celest.
  Mech. Dyn. Astron.* 83:291–302 (2002) — the specific reference this
  paper's own §3.2 cites for "further information" on the Short-/Long-
  Period-family resonant/interconnection web at L4/L5 (a topic this
  paper's own text explicitly defers: "beyond the scope of the current
  paper"); Henrard, J., "Proof of a Conjecture of E. Strömgren," *Celest.
  Mech.* 7:449–457 (1972), and Henrard, J., "On Brown's Conjecture,"
  *Celest. Mech.* 31:115–122 (1983) — foundational L4/L5-family-web papers
  in the same Henrard lineage, directly cited alongside the 2002 paper.
- *Medium:* Farquhar, R., "The Control and Use of Libration-Point
  Satellites," PhD thesis, Stanford (1968) — the historical origin of
  practical libration-point mission use, cited as prior-art background;
  Lo, M. W. et al. (8 co-authors), "Genesis Mission Design," *J.
  Astronaut. Sci.* 49:169–184 (2001) — a concrete worked mission-design
  application of the L1 halo-family machinery this paper surveys, directly
  cycler-adjacent by lineage (same Genesis/KLMR family this project's own
  corpus already carries extensively, e.g. `koon-lo-marsden-ross-1999-
  genesis-trajectory-heteroclinic-AAS-99-451.pdf`, already in corpus) — a
  genuine, not-yet-executed cross-reference opportunity (does this
  already-corpused Genesis paper cite Doedel's AUTO-based family
  computations directly? not checked this pass); Castellà, E. & Jorba, À.,
  "The Lagrangian Points of the Real Earth-Moon System," in *Proc.
  Equadiff 2003* (2005) — a real-ephemeris (non-CR3BP) treatment of the
  same L4/L5 stability question this paper's own §4.10 discussion section
  raises for the idealized CR3BP case.
- *Low-medium:* Belbruno, E., Llibre, J. & Ollé, M., "On the Families of
  Periodic Orbits which Bifurcate from the Circular Sitnikov Motions,"
  *Celest. Mech. Dyn. Astron.* 60:99–129 (1994) — a THIRD independent
  Llibre-authored paper surfacing in this single task (alongside paper 1's
  own author and the Gómez-Llibre-Martínez-Simó books above), on a
  structurally adjacent bifurcation-from-a-degenerate-family topic (the
  `mu=1/2` equal-mass case this paper's own §6.3 explicitly cites this
  exact 1994 paper for: "The bifurcations from the Vertical family for the
  case `mu=1/2` were studied by [Belbruno et al., 1994]"); Uphoff, C. W.,
  "The Art and Science of Lunar Gravity Assist," in *Advances in the
  Astronautical Sciences* Vol. 69, AAS 89-170 (1989) — the source of the
  "Backflip maneuver" name this paper's own **B**i family taxonomy is
  explicitly named after (§3.6: "named after the so-called 'Backflip
  maneuver' described in [Uphoff, 1989]").
- *Low:* Nesvorný, D. & Dones, L., "How Long-Lived are the Hypothetical
  Trojan Populations of Saturn, Uranus and Neptune," *Icarus* 160:271–288
  (2002) — cited only tangentially (Trojan-asteroid-stability discussion,
  §5 opening); Danby, J. M. A., *Fundamentals of Celestial Mechanics*,
  Willmann-Bell (1992), and Danby, J. M. A., "Orbits in the Copenhagen
  Problem Asymptotic at L4 and their Genealogy," *Astron. J.* 72:198–201
  (1967) — general textbook + a `mu=1/2`-specific historical predecessor,
  low domain-specificity; the remaining ~40 references are Doedel's own
  AUTO/numerical-continuation-methods self-citation lineage (Ascher et al.,
  Beyn/Champneys/Kuznetsov, Keller 1977, Kuznetsov 1998, Rheinboldt 1986,
  Seydel 1995, etc.) — infrastructure/tool references, out of this
  project's search-method domain per `#730` §9's standing triage
  convention, not flagged individually.

No genuinely surprising high-priority domain gap beyond the
Gómez-Llibre-Martínez-Simó two-volume book series (now doubly motivated —
both by this paper's own citation and by paper 1's shared co-authors) and
the Henrard L4/L5-web lineage this paper's own text explicitly defers to.

---

## Summary answers (for the dispatching session)

- **Paper 1 (Llibre-Martínez-Simó 1985):** proves transversal existence of
  symmetric homoclinic orbits to the L2 Lyapunov family at a countable
  sequence of mass ratios `mu_k ~ 1/(N(inf)^3 k^3)` (`N(inf)=5.1604325`,
  `M(inf)=2.1330587`), extends this to a codimension condition for nonzero
  `Delta C = C2-C`, and proves a full symbolic-dynamics (Bernoulli shift)
  embedding for quasirandom Sun<->Jupiter shuttling orbits (Theorem C).
  **Cross-check verdict:** no contradiction with this project's own
  L2/transit-manifold code (`genome/transit_manifold.py`, `#547`) — that
  module targets transit classification, not homoclinic transversality;
  this project has NO code anywhere that computes an L2-Lyapunov-family
  homoclinic orbit, matching the general homoclinic-continuation gap
  `#742`'s BMO cross-check already flagged, now independently confirmed
  from the L2-specific analytical side. Directly load-bearing (theoretical
  existence proof) for the already-acquired Casoliva Class 2 construction,
  as the master list's own item-28 flag anticipated — now read and
  confirmed, no correction needed to `#725`'s digest.
- **Paper 2 (Doedel et al. 2007):** a broad AUTO-continuation visual survey
  of all elemental libration-point periodic-orbit families (Lyapunov,
  Vertical, Halo, Axial, Backflip, Circular, D1/E1, R2, plus `mu`-near-1/2
  exotica S1/S2/T1/T2/K1-4/X4/X5) across the full `mu in (0,1/2]` range,
  with branch points tracked as continuable LOCI in `mu` (its own central
  technical contribution) and no printed numeric IC tables anywhere.
  **Cross-check verdict:** this project's `search/pseudo_arclength.py` is
  the same general pseudo-arclength method independently implemented (not
  AUTO-derived); `bifurcation_detector.py` detects bifurcation BRACKETS at
  fixed `mu`, materially simpler than Doedel's own continued branch-point
  LOCI — this project has no code that continues a branch-point locus in
  `mu` the way Doedel's Figs. 28-39 do. This project's existing halo/axial/
  vertical family access is JPL-catalog-sourced, not from this paper; the
  already-corpused Howell 1984 halo survey is confirmed to be exactly the
  paper Doedel's own text credits as its halo-family predecessor. Kumar &
  Moreno 2025's own sequel-framing (`#728`'s flag) is confirmed, not
  contradicted; a digit-grade critical-`mu`-vs-bifurcation-network-node
  comparison remains an unexecuted future cross-check.
- Citation-mining across both papers surfaced the **Gómez-Llibre-Martínez-
  Simó two-volume libration-point mission-design book series** as the
  single strongest new gap (independently motivated from BOTH papers — by
  paper 1's own shared co-authors and by paper 2's own direct citation),
  plus the two direct Llibre & Simó 1980 predecessor papers to paper 1
  itself (never previously flagged in this project's citation-mining
  history) and the Henrard L4/L5-orbit-web lineage paper 2's own text
  explicitly defers to.
