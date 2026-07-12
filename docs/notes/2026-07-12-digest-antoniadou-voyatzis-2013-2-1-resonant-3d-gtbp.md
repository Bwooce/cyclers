# Digest: Antoniadou & Voyatzis 2013 — 2/1 resonant periodic orbits in 3D planetary systems

Single-paper digest. Read 24/24 pages of the user-supplied PDF on 2026-07-12 AET.
Filed to the private corpus as
`antoniadou-voyatzis-2013-2-1-resonant-periodic-orbits-3d-planetary-systems-cmda-115-doi-10.1007-s10569-012-9457-4.pdf`.

## 1. Header

- **Title (verbatim)**: *2/1 resonant periodic orbits in three dimensional
  planetary systems*
- **Authors**: K. I. Antoniadou, G. Voyatzis — Dept. of Physics, Aristotle
  University of Thessaloniki
- **Venue**: *Celestial Mechanics and Dynamical Astronomy* **115** (2013)
  161-184
- **DOI**: 10.1007/s10569-012-9457-4 (printed directly on p. 161, header —
  sourced from the page, not looked up)
- **Received/revised/accepted**: 2 Aug 2012 / 14 Sep 2012 / 30 Oct 2012

## 2. What the paper actually is — and a corpus cross-reference this closes

**This is the GENERAL (non-restricted) spatial three-body problem** — a star
plus TWO planets, both with genuine nonzero mass, in mutual gravitational
interaction (GTBP = "general three-body problem" throughout the paper's own
terminology). This is model-class-distinct from every other Antoniadou paper
already in the corpus (Antoniadou-Libert 2018/2019, both restricted-problem
papers with one massless body) and from the codebase's own ER3BP/CR3BP work,
which is entirely restricted-problem.

**This resolves a dangling citation found while filing.** The existing digest
`docs/notes/2026-06-25-digest-antoniadou-libert-2019-spatial-resonant.md`
(§5) cites "cf. Antoniadou-Voyatzis 2013 Fig. 22" when discussing a family
type the 2019 paper mentions but does not compute — **that Fig. 22 is this
exact paper's Fig. 22** (p. 182, "An example of a family `F^{1/2}_{g1,i}`
which appears as a bridge between a planar family `f1` and a family
`F^{2/1}_{c,m}`" — confirmed by direct visual match). That digest also
explicitly flagged (§8.1, referenced from `known_corpus_3d.py` review) that
the codebase carries an anchor mislabeled **"Antoniadou-Voyatzis 2018"**
citing arXiv:1811.09442 — but 1811.09442 is actually **Antoniadou & Libert
2019** (Voyatzis is not an author of that paper); the "EARLIER Antoniadou &
Voyatzis 2013/2014/2017 papers" were named as the real Antoniadou-Voyatzis
work but never acquired. This paper is the 2013 member of that set.

**Important scope caveat — this paper does NOT fix the mislabeled anchor.**
Independently re-checked both live locations of the bug (2026-07-12, not
just trusting the 2019 digest's note):
- `src/cyclerfinder/genome/known_corpus_3d.py` (~line 156-162): an anchor
  citing "2923-2940 (2019); DOI 10.1093/mnras/sty3195; arXiv:1811.09442" —
  the Antoniadou-**Libert** 2019 DOI, still present, still uncorrected.
- `src/cyclerfinder/search/literature_check.py` (~line 1523-1544): a
  SEPARATE `CorpusAnchor` named `"Antoniadou-Voyatzis spatial resonant
  periodic orbits in CR3BP (2018)"`, citation text "'Spatial Resonant
  Periodic Orbits in the **Restricted** Three-Body Problem,' (2018);
  arXiv:1811.09442", **`doi=None`** — same wrong arXiv ID, no DOI at all,
  AND its claimed title ("...in the Restricted Three-Body Problem") does not
  match this 2013 paper's actual title/scope either (this paper is the
  GENERAL, non-restricted problem — see §2 above). The correct paper for
  THIS specific anchor's claimed content is more likely a genuinely distinct
  **Antoniadou & Voyatzis 2017** paper (Celestial Mechanics and Dynamical
  Astronomy 129, per the 2019 digest's own candidate list) — NOT this 2013
  paper, and NOT yet acquired. **Flagging per "if you find more references
  you need just ask": the 2017 CeMDA 129 Antoniadou-Voyatzis paper is a
  genuine candidate acquisition if this anchor is to be properly fixed
  rather than just correctly re-labelled to Antoniadou-Libert 2019.**
- Net effect: two independent `CorpusAnchor` definitions across two files
  both cite the wrong paper for the same arXiv ID, one has no DOI at all —
  this is a live, mechanically-findable corpus-integrity bug of the same
  general class as this session's earlier Russell-Strange 2009
  `KNOWN_CORPUS` gap, just not yet caught by a ratchet. Not fixed as part of
  this digest (out of scope for a read-and-file pass); flagged for the
  Fable review and any follow-up task.

## 3. Model (p. 162-165, §2.1)

Rotating frame `Gxyz` (NOT barycentric — origin at the center of mass `G` of
the star `S` and the INNER planet `P1` only, following Michalodimitrakis
1979b) with z-axis parallel to the total angular momentum. Reduces the
6-DOF two-planet problem to **4 degrees of freedom**
`Π4 = {x1, x2, z2, ẏ2}` via the (conserved) angular momentum integral (Eq.
2-5). Full nonlinear EOM (Eq. 6-9) — genuinely different from any CR3BP/
ER3BP equations in the codebase (both planets appear symmetrically with
their own masses `m1, m2`, and the star's motion is NOT fixed — it moves in
the rotating xz-plane along with `P1`). Numerical convention: `P1` = inner
planet, `P2` = outer planet, `m0+m1+m2 = 1`; one planet is always fixed at
Jupiter's mass `1e-3` while the other varies (the "mass ratio" `ρ = m2/m1`
is the family-continuation parameter throughout).

Orbital elements `(a_i, e_i, i_i, Ω_i, ω_i, λ_i)` per planet; resonant angles
`σ1, σ2, Δϖ = σ1-σ2, ΔΩ = Ω2-Ω1` (Eq., p. 166) characterize the 2/1 (or 1/2)
mean-motion resonance. Two exact symmetries admit reduced-dimension periodic
orbits: **xz-symmetric** (4D space `Π4 = {x10,x20,z20,ẏ20}`, both planets
cross the xz-plane perpendicularly at t=0) and **x-symmetric** (both planets
start perpendicular from the x-axis, `Π'4 = {x10,x20,ẏ20,ż20}`).

## 4. Method: two continuation schemes (§4, p. 171-182)

**Scheme I** (§4.1): start from periodic orbits of the SPATIAL CIRCULAR (or
elliptic) RESTRICTED 3BP (3D-CRTBP/3D-ERTBP — one planet genuinely
massless), then analytically continue in the mass of the initially-massless
body up to a real planetary value. Families denoted `F^{2/1}_{c,m}`
(xz-symmetric) / `G^{2/1}_{c,m}` (x-symmetric), subscript `c` = continued
from the CIRCULAR restricted problem, `m` = continuation w.r.t. mass. **The
starting 3D-CRTBP families themselves (Fig. 9) are NOT computed in this
paper** — they cite Hadjidemetriou & Voyatzis (2000) and Kotoulas (2005) as
the actual restricted-problem source (a genuinely different paper/data
source from this one, not yet checked against the corpus).

**Scheme II** (§4.2): start from vertical-critical periodic orbits (v.c.o.,
`|a_v|=1` in the monodromy-matrix-derived index, Eq. 15) of the PLANAR
general problem (2D-GTBP — both planets already massive, but coplanar), then
continue in inclination by increasing `z2`/`ż2`. Families denoted
`F^{1/2}_{g1,i}` / `G^{1/2}_{g2,i}` (subscript `g` = generating planar
family number, `i` = continuation from plane to space). **Tables 1-3 (p.
178) give sourced numeric v.c.o. eccentricities** `(e1, e2)` for planar
families `f1` and `f2` across `ρ ∈ [0.01, 20]` — the only fully-tabulated
numeric IC data in this paper (everything else is characteristic-curve
figures, not tables — flagged as figure-derived / not independently
digitizable to a numeric golden without further plot-reading).

## 5. Headline results (§6 Conclusions, p. 182-183)

- Most periodic orbits found are **linearly unstable**; ALL x-symmetric
  orbits are unstable.
- **xz-symmetric families DO contain stable segments**, up to **mutual
  planetary inclination Δi ≈ 50°-60°** (families bifurcating from the stable
  planar family `f1`, mass ratio `ρ > 0.12`, up to `ρ = 20` — very wide
  mass-ratio range).
- The paper's own closing claim (p. 183): stable resonant periodic orbits at
  such inclinations "may be considered as strong candidates for hosting
  planetary systems" — i.e. this is motivated by explaining REAL observed
  mutually-inclined exoplanet systems (paper cites υ And, 47 UMa as examples
  in the intro, p. 161-162), not spacecraft trajectory design. No
  spacecraft, no flyby, no cycler-adjacent framing anywhere in the paper —
  this is pure celestial-mechanics/exoplanet-dynamics literature.
- Chaotic/unstable evolution near unstable periodic orbits is demonstrated
  via DFLI (de-trended Fast Lyapunov Indicator) time series (Figs. 5-8);
  planet-planet close encounters and scattering occur on ~1e4-1e5 time-unit
  timescales when starting even moderately off a stable periodic orbit
  (Fig. 7).

## 6. Relevance to the cyclerfinder codebase (context for the Fable review)

- **Model-class gap, not a capability-frontier gap**: the codebase's own
  #286 Track-A axes (3D/BCR4BP/QP/epoch-MGA/ER3BP, per project memory) are
  ALL restricted-problem or patched-conic constructions — every third body
  in a cyclerfinder search is a spacecraft/negligible-mass object. This
  paper's core content (families of periodic orbits with TWO genuinely
  massive, mutually-perturbing bodies) is a different physical regime the
  codebase has never modeled. Whether that's worth building is a genuine,
  non-obvious question — see the Fable dispatch.
- **Possible tie-in: the project's #549 real-binary (k1,k2) genome sweep**
  (Patroclus-Menoetius, Didymos-Dimorphos, Orcus-Vanth, Eris-Dysnomia — per
  project memory `project_task_numbering_convention` context) already deals
  with genuinely massive secondary bodies. Not yet checked whether that
  genome's dynamical model is closer to this paper's general-TBP formulation
  or still a restricted/patched approximation — worth checking before any
  build decision.
- **Corpus-integrity finding (§2 above) is independent of the capability
  question** and arguably higher-priority: two live `CorpusAnchor`
  mislabelings citing the wrong DOI/no DOI for a "published, don't call this
  novel" literature gate are a direct correctness risk for any future
  spatial-resonant 3D discovery result, following the exact failure pattern
  already seen this session with Russell-Strange 2009.

## 7. References cited (p. 183-184) — not independently checked against corpus

Most-relevant, not yet cross-referenced: Hadjidemetriou, J. & Voyatzis, G.
(2000), "The 2/1 and 3/2 resonant asteroid motion: a symplectic mapping
approach," Celest. Mech. Dyn. Astron. 78, 137-150 — the actual source of the
3D-CRTBP starting families (Fig. 9) this paper continues from Scheme I.
Kotoulas, T. A. (2005), "Three dimensional periodic orbits in exterior mean
motion resonances with Neptune," Astron. Astrophys. 429, 1107-1115 — also a
3D-CRTBP source for Scheme I starting orbits. Michalodimitrakis, M. (1979b),
"On the continuation of periodic orbits from the restricted to the general
three-body problem," Celest. Mech. 19, 263-277 — the frame/reduction paper
this entire model rests on (Eq. 1). Voyatzis, G. (2008), "Chaos, order, and
periodic orbits in 3:1 resonant planetary dynamics," Astrophys. J. 675,
802-816 — the DFLI diagnostic source.

End of digest.
