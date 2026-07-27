# Digest: Moreno, Aydin, van Koert, Frauenfelder & Koh (2024), "Bifurcation Graphs for the CR3BP via Symplectic Methods: On the Jupiter–Europa and Saturn–Enceladus Systems"

**Task**: #728 (mining `https://bhanukumar314.github.io` for corpus gaps; this paper was found
indirectly — it is the methodology paper cited by the `cz-index-matlab` code repo also linked
from that site, as the source of the Conley-Zehnder-index algorithm the MATLAB code implements).

**Citation**: Moreno, A., Aydin, C., van Koert, O., Frauenfelder, U., Koh, D., "Bifurcation Graphs
for the CR3BP via Symplectic Methods: On the Jupiter–Europa and Saturn–Enceladus Systems," *The
Journal of the Astronautical Sciences* 71, article 51 (2024). DOI `10.1007/s40295-024-00462-7`.
Open access.

**File**: `kumar-pdf/papers/moreno-aydin-vankoert-frauenfelder-koh-2024-bifurcation-graphs-cr3bp-symplectic-jas-71-51.pdf`
(private `cyclers_pdf` repo). 48 pages. Native PDF with a full text layer (`pdftotext -layout`
verified clean throughout, including the equation-heavy §2/§3 — no OCR needed).

## What the paper actually is

A **methods** paper, not a search/discovery paper. It develops and applies a symplectic-geometry
toolkit for classifying periodic-orbit families and their bifurcations in the CR3BP, building on
prior work by the same authors (Moreno & Frauenfelder's GIT-sequence paper, and Frauenfelder/Koh/
Moreno's numerical follow-up — neither of which appears to be in this corpus; flagged below).
Four tools:

1. **Floer numerical invariants** — integers invariant across bifurcations (one for arbitrary
   periodic orbits, one for symmetric orbits), a consistency check on any detected index jump.
2. **B-signs** — a ± sign attached to each elliptic/hyperbolic Floquet multiplier of a **symmetric**
   orbit, generalizing the classical Moser-Krein signature (previously elliptic-only) to also cover
   hyperbolic multipliers. Theorem (Frauenfelder & Moreno, cited not reproduced here): a planar
   symmetric orbit is **negative hyperbolic iff the B-signs of its two symmetric points differ**.
3. **GIT-sequence / global topological methods** — refines the classical Broucke stability diagram
   by adding B-sign data; used to produce "GIT-Broucke" plots (their Fig. 24) showing how a family's
   stability crosses the plane.
4. **Conley-Zehnder (CZ) index** — a winding number extracted from the topology of the symplectic
   group, constant along a non-degenerate family and jumping exactly at bifurcations. §3 gives a
   genuinely novel **numerical algorithm** for computing it from purely LOCAL orbit information
   (monodromy matrix only) — no need to construct an explicit path back to a known reference
   orbit/family. Implementation detail: constructive symplectic-path extension via **Iwasawa
   decomposition + SVD** for numerical stability (§3.1-§3.3, with pseudocode in §3.2). This is the
   algorithm the `cz-index-matlab` GitHub repo (Kumar & Moreno's companion code for a *different*,
   2025 AAS paper — see cross-reference below) implements/translates.

## Numerical results (§4 + Appendix A)

Applies the toolkit to Jupiter-Europa and Saturn-Enceladus (both modelled as the planar/spatial
CR3BP), by deforming known families from Hill's lunar problem (Aydin's prior work) into each real
system. Four headline results:

- **Result I** (§4.1): the Hénon `g` family (planar direct, doubly-symmetric) undergoes a pitchfork
  bifurcation in Hill's problem; deforming to Jupiter-Europa, the pitchfork degenerates into two
  distinct branches (`g-LPO1` and `DPO-LPO2`, the latter a birth-death bifurcation). Full IC tables
  in Appendix B, Tables 4-7 (Jacobi constant, x(0), z(0), period, stability type, CZ-index — SOURCED,
  page-cited, e.g. Table 5 DPO branch: `C=3.00109352, x(0)=1.00470170, z(0)=0.09778837` at the
  index-jump point).
- **Result II** (§4.2): the simple-closed DPO orbit's spatial CZ-index jumps by +1 at
  `C≈3.00109352`, marking a genuine planar-to-spatial bifurcation (Fig. 13).
- **Result III** (§4.3): a full bifurcation graph connecting **prograde and retrograde** planar
  families (`LPO2³` ↔ `DRO⁵`) via a spatial CZ-index-15 family (Fig. 16, 18) — a structural
  connection the authors state was previously known only for one branch (cite [15],[19]); the rest
  are novel. IC tables: Appendix B Tables 8-11.
- **Result IV** (§4.4): the Saturn-Enceladus system reproduces the EXACT SAME bifurcation-graph
  topology as Jupiter-Europa (Fig. 19), just at different energy values — offered as evidence the
  method/topology generalizes across mass ratios, not merely a Jupiter-Europa coincidence.
- **Appendix A (the mission-relevant highlight)**: continuing the Saturn-Enceladus L2 Halo family
  past a birth-death degeneracy into a family of **polar orbits**, tracking altitude above
  Enceladus's surface. **Sourced constants used**: Enceladus semi-major axis `237,948 km`, Enceladus
  radius `252.1 km` (both plain-stated, no further citation given in-paper — treat as the authors'
  own working values, not independently re-verified here). **Headline orbit**: `C=3.000034709155895,
  x(0)=1.0025751548678687, z(0)=0.004882249068671777`, altitude **29 km** above the surface, CZ-index
  just jumped to 4, type `(E²)` (stable) — explicitly named as relevant to the real Enceladus
  Orbilander mission concept (cite [38], NASA/APL 2020) for sampling the water plumes. A second,
  deeper orbit at `C=3.000034757415899, x(0)=1.0024991770058109, z(0)=0.0048974554261678876`,
  altitude **14 km**, is found by detecting a CZ-index jump (20→22) on the Halo family's 7-fold cover
  and searching nearby for the missing indices (21) predicted by Floer-invariant conservation — a
  worked demonstration of the "find new orbits by finding index jumps" method the paper advocates.
  Full family data: Tables 1-3.

## Relevance to this codebase

**No CZ-index / B-sign / GIT-sequence machinery exists anywhere in `cyclerfinder`** (confirmed by
direct grep of `src/` and `docs/notes/` for these terms — zero hits). This project's own
bifurcation/stability code (`src/cyclerfinder/search/bifurcation_detector.py`,
`cr3bp_3d_family_tracer.py`, `nrho_continuation.py`, `er3bp_floquet.py`) works purely with raw
Floquet multipliers (eigenvalues of the monodromy matrix) — it detects bifurcations by tracking
multiplier crossings of the unit circle / real axis, but has no topological classification layer
analogous to the CZ-index or B-signs.

**Directly on-point to TODAY's own `#632` work** (same session, same day): `#632` fixed a
cross-platform BLAS eigenvector-SIGN ambiguity in `genome/qp_tori.py::_seed_invariant_circle` — the
Neimark-Sacker eigenvector's conjugate-pair member (`lam` vs `conj(lam)`) was backend-dependent,
canonicalized to the positive-imaginary representative. This paper's **B-sign** concept (§2.2) is
formally the same class of problem one level up: a well-defined ± invariant attached to a Floquet
eigenvector at a SYMMETRIC point, used to *predict* bifurcation type rather than just describe it
after the fact. `#515`'s prior "aligned 3D Floquet eigenvector signs...to enforce consistent
manifold directions" fix is the same lineage again. All three (`#515`, `#632`, this paper's B-signs)
are instances of the same underlying fact: **Floquet/monodromy eigenvectors need an explicit,
physically-motivated sign/orientation convention to be well-defined and platform-reproducible** —
worth keeping in mind as a recurring theme, not a coincidence.

**Concrete future-work opportunity (not built here, flagging only)**: if this project ever wants a
principled, index-theoretic way to (a) predict WHERE a periodic-orbit family will bifurcate before
brute-force stepping into it, or (b) classify "family A connects to family B" claims with a
topological invariant instead of eyeballing continuity of raw ICs, the CZ-index/GIT-sequence
toolkit here is a candidate — but it is a genuinely new class of infrastructure (needs the full
Iwasawa/SVD symplectic-path-extension algorithm of §3, not a drop-in function), not a quick port.
No such build is scoped by this digest; this is background grounding only, per this project's
"digest ≠ adoption" discipline.

**No catalogue overlap**: Jupiter-Europa DPO/LPO/DRO and Saturn-Enceladus Halo/polar families
searched here do not correspond to any existing `catalogue.yaml` row (grepped `DPO`, `LPO2`,
`DRO5`, "Europa halo", "Enceladus halo" — no hits). This is expected: these are single-primary
(3-body, moon-centered) stationkeeping/reconnaissance orbit families, not cyclers — orthogonal to
this project's cycler/quasi_cycler/precursor_mga/mga_tour taxonomy, not a catalogue gap.

## Citation-mining flags (not acquired)

- **Moreno & Frauenfelder** (GIT-sequence original paper, ref [6]) — the actual mathematical
  groundwork this whole paper builds on. Not found in `CORPUS_INDEX.md` under any author/title
  match. Candidate gap if this project ever pursues the CZ-index/GIT-sequence direction above.
- **Frauenfelder, Koh & Moreno** (ref [9], the first application of the toolkit to numerics) — same
  status, not found in corpus.
- **Aydin** (Hill's lunar problem family classification, ref [18]) — not found in corpus; this is
  the deformation starting point for all four numerical results in this paper.
- **Restrepo & Russell**, "A database of planar axi-symmetric periodic orbits for the solar system,"
  AAS 17-694 (2017) (ref [35]) — an online IC database the paper explicitly cross-validates its
  Jupiter-Europa DPO/LPO/DRO families against ("already known and appear e.g., in page 12 of
  [35]"). Not found in corpus. Could be a useful independent-source cross-check database for any
  future single-primary periodic-orbit family work, flagged for visibility only.
- None of these are pursued further here — this digest's scope is the one paper found via the
  `#728` site-mining pass; deeper citation acquisition is a separate future task if this direction
  is ever prioritized.

## Cross-reference

The companion code repo `github.com/bhanukumar314/cz-index-matlab` (assessed separately, same
`#728` task, see the coordinating session's own report) cites THIS paper for its methodology and a
*different* 2025 AAS paper (Kumar & Moreno, "Networks of Periodic Orbits in the Earth-Moon System
Through a Regularized and Symplectic Lens," AAS-25-677 — also being digested under `#728`, see that
paper's own digest note) for the MATLAB implementation itself. All three (this paper, the AAS-25-677
paper, and the MATLAB code) form one coherent methodology lineage from the same author group.
