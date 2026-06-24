# Digest: Antoniadou & Libert 2019 — Spatial resonant periodic orbits in the RTBP

Single-paper digest. Read 19/19 pages of the PDF on 2026-06-25 AET. Digest
supports the #434/#438 3D broken-plane cycler novelty gate (the spatial-
resonant-family corpus anchor). READ-and-EXTRACT only.

## 1. Header

- **Title (verbatim)**: *Spatial resonant periodic orbits in the restricted
  three-body problem*
- **Authors (verbatim, confirmed)**: **Kyriaki I. Antoniadou** and **Anne-Sophie
  Libert** — NaXys, Department of Mathematics, University of Namur, 8 Rempart de
  la Vierge, 5000 Namur, Belgium. Corresponding author e-mail
  `kyriaki.antoniadou@unamur.be` (KIA). So the author list is **Antoniadou &
  Libert**, NOT "Antoniadou & Voyatzis" (Voyatzis is a co-author on the EARLIER
  related papers this one extends — see §9 errata).
- **Venue**: *Monthly Notices of the Royal Astronomical Society* (MNRAS)
  483(3):2923–2940 (2019). Preprint header reads "MNRAS 000, 1–18 (2018)" /
  "Preprint 6 December 2018".
- **DOI**: 10.1093/mnras/sty3195
- **arXiv**: 1811.09442 (v2, 4 Dec 2018; astro-ph.EP)
- **Length**: 19 pages in this PDF (incl. Appendix A + references).

## 2. What the paper actually is

This is the **spatial (out-of-plane) resonant-periodic-orbit catalogue** for the
restricted three-body problem — the systematic study of how 3D resonant families
arise, where they bifurcate from planar families, and their linear stability vs
inclination/eccentricity. It is the spatial-RTBP companion to Antoniadou &
Libert's earlier *planar* stability studies (2018a A&A 615 A60 for 2/1;
2018b CeMDA 130:41 for 3/2, 5/2, 3/1, 4/1, 5/1).

The physical framing is **exoplanetary**, not Earth-Moon: a star + an **inclined
inner terrestrial (massless) planet in MMR** + an **outer giant**. But the paper
states explicitly (§2.1, p. 2 and Conclusions p. 17) that the results apply to
**any** system modellable by the RTBP — "star-asteroid-giant planet,
planet-spacecraft-satellite and binary star-circumprimary planet" — and Sect.
5.2 applies them to asteroid dynamics. This transferability claim is what makes
the family taxonomy relevant to us (see §7 below for the honest assessment).

It is a **continuation/bifurcation + stability** paper. It does NOT publish
tabulated initial conditions: there is **no IC table anywhere in the paper**.
All families are presented graphically (Figs. 2–33) as characteristic curves in
projected planes — e.g. (e1, i1), (e1, e2, i1), (x, e1) — with bifurcation
points called out by coordinate value in the text/captions. See §4 for what IS
and is NOT tabulated. This matters: we cannot match a swept 3D orbit against a
printed state vector from this paper; we can only match the (resonance,
configuration, prograde/retrograde, inclination-range, stability) signature.

## 3. The model (§2, pp. 1–3)

- **Model**: restricted three-body problem, computed at **two fidelities**:
  - **3D-CRTBP** (circular RTBP): giant on a **circular** orbit (e2 = 0).
    Three degrees of freedom; primaries fixed at (−µ, 0) and (1−µ, 0);
    uniformly rotating Oxyz frame; r = 1, θ̇ = 1.
  - **3D-ERTBP** (elliptic RTBP): giant on an **elliptic** orbit (e2 ≠ 0).
    Non-autonomous, four DOF; frame not uniformly rotating.
  Planar versions (2D-CRTBP, 2D-ERTBP) are the z = ż = 0 restriction and supply
  the *generating* families.
- **Mass ratio / system**: masses normalised m0 = 1 − m2 = 0.999, **µ = m2/(m0+m2)
  = 0.001** (§2.1, p. 2). The inner terrestrial planet is the **massless** body
  P1 (m1 = 0); the outer giant is m2 = mJ. So µ = 1e-3 is a **single
  Jupiter-mass-giant value**, NOT a sweep over mass ratios and NOT the Earth-Moon
  ratio (µ_EM ≈ 0.01215, ~12× larger). This is the single most important
  transferability caveat — see §7.
- **Resonances covered (interior MMRs, a1 < a2, subscript 1 = inner body)**:
  **3/2, 2/1, 5/2, 3/1, 4/1, 5/1** (abstract; §4). All are *interior* (inner
  body resonant with outer giant). The period of the spatial orbits in
  normalised units is T0 = 2π.
- **Resonant angles** (Eq. 3, p. 3): θ1, θ2 (q odd: 3/2, 2/1, 5/2, 4/1 — note
  4/1 here is order-3 odd in their counting q=|p−q|... they list q-parity per
  resonance), θ3 (q=2, the 3/1 MMR), θ4 (q=4, the 5/1 MMR), built from
  λi = Mi + ϖi. **Inclination-resonance** angles φ11, φ22, φ12 (Eq. 4) and the
  nodal/secular angle φΩ = Ω1 − Ω2 (Eq. 5) characterise the spatial regime.
- **Symmetry classes** (§2.3, p. 3): two families of *symmetric* spatial orbits,
  defined by perpendicular crossings of a Poincaré section π̂ = {y = 0, ẏ > 0}:
  - **xz-symmetric**: two perpendicular crossings of the xz-plane; reduced ICs
    {(x0, z0, ẏ0)} (Eq. 7). Geometry: Ω1 = 90°, ϖ1 = 0° or 180°.
  - **x-symmetric**: two perpendicular crossings of the x-axis; reduced ICs
    {(x0, ẏ0, ż0)} (Eq. 8). Geometry: Ω1 = 0° or 180°, ϖ1 = 0° or 180°.
  Four symmetric configurations from the resonant-angle pairs: (0,0), (0,π),
  (π,0), (π,π). **Asymmetric** spatial orbits are explicitly NOT computed (left
  to future work, Conclusions p. 17).

## 4. HOW spatial families arise — the vertical-critical-orbit mechanism (KEY)

This is the core structural content for our k_z topology gate.

### 4.1 The bifurcation mechanism (§2.4 p. 3, §3 p. 4)
- Every planar periodic orbit carries a **vertical stability index a_v** (Hénon
  1973). The planar orbit is **vertically stable** when |a_v| < 1, **vertically
  unstable** when |a_v| > 1.
- The transition points where **|a_v| = 1** are the **vertical critical orbits
  (v.c.o.)**. These are exactly the points on a planar family from which a
  **spatial (z ≠ 0) family bifurcates** and can be analytically continued out of
  the plane. This is the mechanism: *spatial families emanate from v.c.o. on
  planar families.*
- Plot convention used throughout (§2.4 p. 4): **solid line = vertically stable**
  planar orbit, **dashed line = vertically unstable**; blue = horizontally
  stable, red = horizontally unstable (planar) / linearly (un)stable (spatial);
  **magenta dots = v.c.o. generating xz-symmetric spatial families**; **green
  dots = v.c.o. generating x-symmetric spatial families**.

### 4.2 Two continuation methods (§3 p. 4)
- **Method I**: examine the *circular family* and the planar families of the
  2D-CRTBP/2D-ERTBP for vertical stability; the v.c.o. are the bifurcation
  points that generate spatial periodic orbits continued into the 3D-CRTBP and
  3D-ERTBP. (This is the primary route used.)
- **Method II**: examine spatial families of the 3D-CRTBP at periods T = k·T0/m
  (m = multiplicity); these orbits are bifurcation points generating spatial
  orbits in the 3D-ERTBP. From each such point, in general **two** families
  emanate — one for the giant at pericentre, one at apocentre (since e2 goes
  from 0 to ≠ 0).

### 4.3 The circular family as v.c.o. generator (§4.1 p. 5, Appendix A)
- Along the **circular family** (both primaries + massless body on circular
  orbits, MMR varying along the family), v.c.o. appear **systematically as the
  multiplicity i of the circular periodic orbits increases** (Appendix A, Fig.
  A1, Table A1). The circular family "breaks" at first-order resonances (q=1:
  3/2, 2/1) so NO v.c.o. exist there along the circular family — hence **3/2 and
  2/1 cannot be generated from the circular family via Method I** (no spatial
  family from C_i at those MMRs; they instead come from the elliptic planar
  families). 5/2 first appears at multiplicity i=3; 4/1 at i=3; 5/1 at i=2; 3/1
  has v.c.o. at every multiplicity (it is the only resonance studied with **two**
  v.c.o. even at low multiplicity). Table A1 lists, per multiplicity i=1..9,
  which MMRs possess a v.c.o. — a structural map, NOT initial conditions.

### 4.4 Family nomenclature (§4 p. 4–5)
- Spatial families of **xz-symmetric** orbits: **F**, with superscript = MMR,
  e.g. F_I^{3/2}. Spatial families of **x-symmetric** orbits: **G**, e.g.
  G_I^{3/2}. The subscript labels which planar family (I, II, or a configuration
  pair like (π,0)) or which circular-family branch (C1, C2, C3) it bifurcates
  from. A **circumflex** (F̂, Ĝ) marks the bifurcation *point* (the v.c.o.) as
  distinct from the spatial family. A **prime** (F′, G′) distinguishes a second
  family/point originating from the same planar family of a given configuration.
  In the 3D-ERTBP a trailing p or a denotes giant-at-pericentre vs -apocentre.

## 5. What is tabulated vs graphical (for matching)

- **No initial-condition table exists.** There is exactly one table in the body
  (**Table A1**, Appendix A) and it is a **resonance-vs-multiplicity occurrence
  map** of which MMRs possess a v.c.o. — it lists MMR *labels* (e.g. 7/3, 11/5,
  5/2, 3/1, 4/1, 5/1, ...) against multiplicity i = 1..9. No (x0, z0, ẏ0), no
  Jacobi constant, no stability column. It cannot be used to match a numerical
  orbit.
- All families, bifurcation eccentricities, and inclination ranges are given
  **graphically** (Figs. 2–33) and as **named coordinate points in captions/text**
  (e.g. "F_I^{3/2} bifurcates from the v.c.o. at e1 = 0.39"; "G_I^{3/2}
  bifurcates to the 3D-ERTBP at (e1, i1) = (0.42, 69°)"). These bifurcation
  coordinates ARE extractable per-resonance (see §6) but they are
  curve-endpoints, not orbit ICs.
- **Verdict for our matcher**: we can populate (resonance, configuration,
  prograde/retrograde, inclination-stability-range) signatures, NOT exact
  state-vector anchors. This is consistent with how `known_corpus_3d.py` already
  works — it matches on `(k1, k2, k_z)` + Jacobi-band, never on fabricated ICs.

## 6. Stability findings per MMR (the headline result, Conclusions p. 17)

The Conclusions give the clean inclination-range summary of where **stable**
spatial periodic orbits exist. This is the most directly usable extract.

**3D-CRTBP (giant circular) — STABLE spatial orbits exist for:**
- **Prograde** (i1 measured from 0°):
  - 3/2 MMR: i1 < 55°
  - 2/1 MMR: i1 < 90°
  - 5/2 MMR: i1 < 40°
  - 3/1 MMR: 56° < i1 < 82° (different ranges → different families)
  - 4/1 MMR: i1 < 45°
  - 5/1 MMR: i1 < 37°
- **Retrograde** (orbits continued toward i1 = 180°):
  - 3/2 MMR: 103° < i1 < 165°
  - 2/1 MMR: 95° < i1 < 180°
  - 5/2 MMR: 162° < i1 < 178°
  - 3/1 MMR: 104° < i1 < 138° and (a separate family) i1 > 142°
  - 4/1 MMR: 108° < i1 < 163°
  - 5/1 MMR: 144° < i1 < 174°

**3D-ERTBP (giant elliptic — only prograde continued):**
- Starting from the v.c.o. of the 2D-ERTBP: in the **2/1 MMR** the family in
  config (π,0) is **whole stable and reaches 70°** inclination (the standout
  stable elliptic case — the abstract singles this out: "only the 2/1 resonance
  has stable periodic orbits up to high inclinations").
- The **3/1, 4/1, 5/1** MMRs possess **segments of stability up to ~10°** in
  config (0,π) (abstract: "possess segments of stability for low inclinations").
- Families continued from a **3D-CRTBP** bifurcation point (Method II): **all
  found unstable**.

Additional regular-motion finding (Sect. 5, DS-maps): regular (quasi-periodic)
motion also occurs in the neighbourhood of **both horizontally- AND
vertically-stable planar** periodic orbits, even at very high inclinations —
i.e. stability is not exclusively tied to the spatial families. DS-maps (Figs.
31–33) use the detrended Fast Lyapunov Indicator (DFLI), log(DFLI) ≤ 2 ⇒
regular/stable.

### Per-resonance bifurcation coordinates extractable from figures/captions
(curve endpoints, not ICs; for cross-reference only):
- **3/2**: v.c.o. on planar family I at e1 = 0.39 (F̂_I^{3/2}) and e1 = 0.43
  (Ĝ_I^{3/2}); family IIS v.c.o. at e1 = 0.98; G_I^{3/2} bifurcates to 3D-ERTBP
  at (e1,i1) = (0.42, 69°). xz-family F_I^{3/2} stable up to ~55° and between
  103°–165°. (Figs. 2–6.)
- **2/1**: planar family I has v.c.o. pair F̂_I^{2/1} (e1=0.67) and Ĝ_I^{2/1}
  (e1=0.79); xz-family F_I^{2/1} stable prograde & retrograde except near 90°;
  3D-ERTBP F_(π,0)^{2/1} bifurcates at (e1,e2)=(0.93,0.43), stable. (Figs. 7–10.)
- **5/2**: first v.c.o. at multiplicity 3 (C3); F_{C3}^{5/2} stable up to 40°
  prograde, 109°<i1<164° retrograde; ERTBP bifurcation at (e1,i1)=(0.84,29°);
  config (π,0) v.c.o. Ĝ at (e1,e2)=(0.96,0.61). (Figs. 11–16.)
- **3/1**: the ONLY MMR with two v.c.o. at all multiplicities; four spatial
  families F_{C1}^{3/1}, G_{C1}^{3/1}, F′, G′ from circular family C1;
  bifurcation pts to 3D-ERTBP at (e1,i1) = (0.000178,84°), (0.0002327,101°),
  (0.0003454,101°), (0.0004451,84°). Config v.c.o. at (0.59,0.43), (0.90,0.76),
  (0.78,0.05), (0.87,0.35). (Figs. 17–22.)
- **4/1**: first v.c.o. at multiplicity 3 (C3); F_{C3}^{4/1} stable to 45°
  prograde, 108°<i1<163° retrograde; 3D-ERTBP bifurcation at (e1,i1)=(0.113,89°).
  (Figs. 23–26.)
- **5/1**: first v.c.o. at multiplicity 2 (C2); F_{C2}^{5/1}/G stable i1<37°
  prograde, 144°<i1<174° retrograde; four bifurcation pts to 3D-ERTBP near
  i1=90° at e1 ≈ 1e-4. (Figs. 27–30.)

## 7. Transferability to Earth-Moon cyclers — honest assessment

This is the crux for #434/#438. **Verdict: structurally informative,
quantitatively non-transferable; it is a TAXONOMY anchor, not a coordinate
anchor, and even the taxonomy maps only partially.**

1. **Mass ratio mismatch (decisive).** The paper is a single µ = 1e-3
   (Jupiter-mass giant) study. Earth-Moon is µ ≈ 0.01215, ~12× larger, and the
   "giant" is the secondary the cycler must ENCOUNTER, not avoid. Bifurcation
   eccentricities, v.c.o. inclinations and stability ranges all depend on µ;
   none of the numeric (e1, i1) bifurcation coordinates in §6 transfer to
   Earth-Moon. They are exoplanetary-system numbers.

2. **Resonance semantics differ.** Here "p/q MMR" is the inner massless planet's
   mean-motion ratio with the OUTER giant — a heliocentric-style two-body
   resonance with a perturber on a wide orbit. Our cycler `(k1, k2)` winding is
   the rotating-frame revolution count of an Earth-Moon trajectory that
   physically transfers between/encounters the primaries. They are NOT the same
   object: a `R_{p:q}` here is an inner planet looping the star q-ish times per
   giant orbit, with the giant a distant perturber — there is no
   Earth↔Moon-style leg. So the paper's (resonance) families do **not**
   correspond one-to-one to our (k1, k2) cyclers. The mapping is structural
   (both are resonant out-of-plane families born at v.c.o.) not literal.

3. **What DOES transfer — the mechanism and the topology class.** The
   *vertical-critical-orbit bifurcation mechanism* is µ-independent in form:
   any planar family in the Earth-Moon CR3BP has its own a_v profile and its own
   v.c.o., and our 3D broken-plane lifts of planar cyclers are exactly
   continuations from (or near) such vertical bifurcations. So the paper
   correctly tells us **(a)** that 3D resonant families are an EXPECTED, KNOWN
   phenomenon (a lifted planar cycler reaching z≠0 via a v.c.o. is not novel
   physics), and **(b)** the symmetry/configuration taxonomy (xz-symmetric F vs
   x-symmetric G; configs (0,0)/(0,π)/(π,0)/(π,π); prograde vs retrograde
   branches terminating at i1=180°) is the right vocabulary for classifying a
   swept 3D family. This is k_z > 0 topology: these orbits cross the equatorial
   plane (z oscillates through 0 at the symmetric crossings).

4. **It does NOT give us Earth-Moon ICs or Jacobi bands.** No state vectors, no
   C values, wrong µ. So it cannot pin a numeric Jacobi band for an Earth-Moon
   anchor — `jacobi_band` must stay `None` for any anchor sourced from this
   paper (which is already the case in the live corpus).

## 8. For our use — known_corpus_3d / #434/#438 novelty gate

### 8.1 CRITICAL CORRECTION to the existing anchor (action required)
`src/cyclerfinder/genome/known_corpus_3d.py` already contains an anchor (lines
119–146) with **`doi="arXiv:1811.09442"`** — which IS this paper — but it is
labelled and cited as **"Antoniadou-Voyatzis 2018"** with
`authors=("Antoniadou", "Voyatzis")` and a citation string that conflates this
work with the *Orbital stability of coplanar two-planet systems* line. **This is
a provenance error.** arXiv 1811.09442 is **Antoniadou & Libert (2019), MNRAS
483(3):2923, DOI 10.1093/mnras/sty3195** — Voyatzis is NOT an author of this
paper (he co-authors the EARLIER Antoniadou & Voyatzis 2013/2014/2017 papers
that this one extends, and Voyatzis-Tsiganis-Antoniadou 2018). Recommended fix
for the parent:
- `name` → "Antoniadou-Libert spatial resonant periodic orbits in the RTBP (2019)"
- `authors=("Antoniadou", "Libert")`
- `doi="10.1093/mnras/sty3195"` (the published DOI; keep arXiv:1811.09442 in the
  citation string as the preprint id)
- `citation` → cite MNRAS 483(3):2923–2940 (2019), DOI 10.1093/mnras/sty3195,
  arXiv:1811.09442; "spatial resonant periodic orbits in the CRTBP/ERTBP, 3/2,
  2/1, 5/2, 3/1, 4/1, 5/1 interior MMRs, families bifurcating from
  vertical-critical-orbits of the planar families."
- Keep `topology_3d={"k1": 1, "k2": 1, "k_z": 2}`, `jacobi_band=None`,
  `topology_label=frozenset({"resonant"})` — those remain defensible (resonant,
  out-of-plane k_z>0, no sourced Earth-Moon C).
- Add keywords: "vertical critical orbit", "xz-symmetric x-symmetric spatial
  family", "3D-CRTBP 3D-ERTBP resonant".

  NOTE: do NOT simply delete/replace — if a genuine *Antoniadou & Voyatzis 2018*
  spatial-resonant paper (CeMDA 130:29, Voyatzis-Tsiganis-Antoniadou 2018, or
  Antoniadou-Voyatzis 2017 CeMDA 129) is ALSO meant to be in the corpus, it
  should be a SEPARATE anchor with its own correct DOI. The current single
  anchor is mislabelled, not necessarily redundant. The parent should decide
  whether to (a) correct the one anchor to Antoniadou-Libert 2019, or (b) split
  into two anchors (Antoniadou-Libert 2019 + a correctly-cited
  Antoniadou-Voyatzis paper). The literature_check golden test
  (`tests/search/test_literature_check_3d.py`) enforces citation+DOI presence,
  so any edit must keep both populated.

### 8.2 What signature this paper authorises the gate to mark "published"
A swept 3D Earth-Moon candidate should be treated as **not-novel / published-
family** (not certified novel) when its fingerprint is an **out-of-plane
(k_z > 0) resonant family that bifurcates from a vertical-critical-orbit of a
planar resonant family**, for the interior MMR set {3/2, 2/1, 5/2, 3/1, 4/1,
5/1} — REGARDLESS of the µ mismatch, because the *mechanism and family class*
are published here. The match is on `(k1, k2, k_z)` topology +
`topology_label={"resonant"}`, NOT on numeric ICs or C. Prograde and retrograde
(i1→180°) branches are both covered.

### 8.3 What this paper does NOT cover (leave novelty OPEN for these)
- **Asymmetric** spatial periodic orbits — explicitly excluded (Conclusions
  p. 17). A swept 3D family that is asymmetric is NOT anchored by this paper.
- **Spatial isolated families** (not connected to a v.c.o. of a planar family) —
  the paper notes (Sect. 5.2) these "can exist" (cf. Antoniadou-Voyatzis 2013
  Fig. 22) but does not compute them. Not anchored.
- **L3/L4/L5-centred, halo/NRHO, vertical-Lyapunov/butterfly** — out of scope
  (those are anchored by the Howell 1984 and Folta et al. 2015 entries already
  in `known_corpus_3d.py`). This paper is resonant-only.
- **Earth-Moon µ specifically** — no member of this paper is at µ_EM; treat the
  anchor as a TOPOLOGY/MECHANISM anchor, not a member-list anchor.

### 8.4 Net effect on the corpus
The corpus already has the right SHAPE (Howell halo + Folta NRHO + this
resonant anchor). The only change this digest forces is the **provenance
correction in §8.1** (Antoniadou-Libert 2019, not Antoniadou-Voyatzis 2018;
published DOI 10.1093/mnras/sty3195). No new anchor topology is needed; no
Jacobi band can be sourced; the `(1,1,2)` resonant tuple stands.

## 9. Errata / cross-reference notes

- **Author attribution (the important one)**: arXiv:1811.09442 =
  Antoniadou & **Libert** 2019 MNRAS, DOI 10.1093/mnras/sty3195. The live
  corpus anchor mis-attributes it to "Antoniadou-Voyatzis 2018". See §8.1.
- This paper **extends** the authors' planar studies: Antoniadou & Libert 2018a
  (A&A 615 A60, 2/1 MMR) and 2018b (CeMDA 130:41, 3/2 5/2 3/1 4/1 5/1) supplied
  the planar families whose v.c.o. generate the spatial families here.
- Prior *spatial* resonant RTBP/GTBP work it builds on: Antoniadou & Voyatzis
  2013 (CeMDA 115:161, 2/1), 2014 (CeMDA 119:221, and Ap&SS 349:657 — 3/2 5/2
  3/1 4/1), 2017 (the 5/2, 7/3 inclination-resonance result). Kotoulas 2005 /
  Kotoulas & Voyatzis 2005 / Voyatzis & Kotoulas 2005 (exterior/trans-Neptunian
  inclined symmetric POs). Voyatzis-Tsiganis-Antoniadou 2018 (CeMDA 130:29,
  asymmetric inclined). Hénon 1973 (A&A 28:415 — the vertical stability index
  a_v). Skokos 2001 / Marchal 1990 (stability of high-dim Hamiltonian POs).
- No methodological contradiction with our stack: the v.c.o. (|a_v|=1)
  bifurcation mechanism and the symmetric-crossing reduction are standard and
  consistent with our CR3BP corrector conventions.

## 10. Action items for the parent (controller)

1. **CORRECT the provenance** of the arXiv:1811.09442 anchor in
   `src/cyclerfinder/genome/known_corpus_3d.py` from "Antoniadou-Voyatzis 2018"
   to **Antoniadou-Libert 2019, MNRAS 483(3):2923, DOI 10.1093/mnras/sty3195**
   (§8.1). Decide whether to correct-in-place or split into two anchors. Re-run
   `uv run pytest tests/search/test_literature_check_3d.py`.
2. **Keep** `topology_3d={"k1":1,"k2":1,"k_z":2}`, `jacobi_band=None`,
   `topology_label={"resonant"}` — sourced-defensible; do NOT invent a Jacobi
   band (wrong µ, no published C).
3. **Record** that this anchor covers interior MMRs {3/2, 2/1, 5/2, 3/1, 4/1,
   5/1}, symmetric (xz + x) families, prograde AND retrograde, born at v.c.o. —
   and explicitly does NOT cover asymmetric or spatial-isolated families (leave
   novelty open for those in #434/#438).
4. **Register** the digest in CORPUS_INDEX.md (controller's job, not this digest).

End of digest.
