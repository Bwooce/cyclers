# Digest: Restrepo & Russell 2018 — JPL planar-axisymmetric CR3BP database

Task #374. Single-paper digest. Read 24/24 pages of the PDF on 2026-06-17 AET.
Closes the LAST item on the #116 acquisition wishlist (9 of 9 complete).

## 1. Header

- **Title (verbatim)**: *A database of planar axisymmetric periodic orbits for
  the Solar system*
- **Authors**: Ricardo L. Restrepo, Ryan P. Russell (both: The University of
  Texas at Austin, 210 E. 24th St., Austin, TX 78712, USA)
- **Venue**: *Celestial Mechanics and Dynamical Astronomy* (2018) 130:49
- **DOI**: 10.1007/s10569-018-9844-6
- **Submitted/Revised/Accepted**: 11 January 2018 / 10 May 2018 / 19 June 2018
- **Published online**: 11 July 2018
- **Length**: 24 pages (front matter through references; the task brief
  estimated ~30 — actual is 24)

## 2. What the paper actually is

The paper is the single archival publication documenting Russell's JPL-aligned
"3BP catalogue" of planar periodic orbits (POs) — the artefact the project's
#116 acquisition wishlist has been tracking since session start. It does NOT
publish the orbits themselves in the PDF; it publishes the **methodology, the
nomenclature, the system coverage, and the online database pointer**, with the
~3 million orbits hosted externally as Fortran-emitted ASCII column files
(p. 20 §5; URLs at the footnotes on p. 20 — see §6 below).

Methodology is a **two-level grid search** (p. 3 §1, p. 5 §3) over a five-
parameter axisymmetric reduction (state at perpendicular x-crossing reduces to
`(x0, ydot0)` plus a half-period parametrised by an integer crossing count
`N`, p. 7 §3.1). Level 1 is a coarse global grid sweeping
`-x0_max < x0 < x0_max` with `x0_max = 5*x_L1` and `ydot0 in [-ydot0_max,
+ydot0_max]`, terminated at the Nth x-crossing for `N = 1, ..., N_max=10`
(p. 7-8 §3.1). Level 2 is a **fine local grid centred on a small set of
"generating families"** — the L1/L2 Lyapunov families `LL1`/`LL2` and the
Distant Retrograde Orbit family `DRO` — implemented as a cubic-spline-
interpolated sweep over `ydot0` perturbations in the neighbourhood of the
generating-family characteristic curve (p. 9 §3.2). The local search is what
captures the highly-sensitive **connecting resonances** that the global search
misses; these are the heteroclinic-approximation orbits that the paper
emphasises as its primary novel contribution (p. 17-19 §4.2).

The differential corrector at every grid-cell convergence is the **full
second-order trust-region method** from Conn-Gould-Toint (2000), which solves
the quadratic minimisation `min phi(x) = g x + x' H x s.t. |x| <= rho` per
Eq. (7) on p. 7; `g` is the 1x4 Jacobian and `H` the 4x4 Hessian of the
periodicity residual `xdot(T/2)^2 + ydot(T/2)^2`, both computed analytically
from the first- and second-order state-transition tensors `Phi^1` and `Phi^2`
of Eqs. (3)-(4) on p. 4. The corrector targets `xdot(T/2) approx 0` and
`y(T/2) = 0` simultaneously (axisymmetric closure) by sweeping the half-period
unknown rather than the full period — the explicit rationale (p. 7) being
robustness against the differential corrector failing when intermediate
iterates have a different number of x-crossings than the seed. The corrector
is described as "robust" because the second-order trust-region globaliser
tolerates seeds far from the basin, which is what makes the brute-force grid
search feasible without manual continuation.

The orbit classification is a **descriptive nomenclature** based on six
geometric counters (Table 2, p. 6): `N_cross` (x-crossings on half-period),
`I_cross` (left-of-primary crossings), `D_cross` (opposite-side-of-secondary
crossings), `M1_rot` and `M2_rot` (rotations around primary and secondary,
signed; `+` prograde, `-` retrograde), plus `T`, `J`, `(b_h, b_v)` stability,
and the perpendicular-crossing states `(x_perp1, ydot_perp1)` and
`(x_perp2, ydot_perp2)`. Table 3 (p. 22, appendix) is the **identification key**
mapping these counters to the family labels (LL1, LL2, LPO1, LPO2, LPO1*,
LPO2*, DPO, DRO, Hg, Hb1, Hb2, Hg(b)/Hg(c)/Hg_{nr}, M1/M2, Hm1/Hm2,
Hm1(b)/Hm2(b)/Hm1(c)/Hm2(c), QDRO_{11}/QDRO_{14}/QDRO(g3), R_{p:q},
H_{R p:q-LL1}, H_{R p:q-LL2}, H_{R p:q-DPO}, H_{R p:q-QDRO}).

## 3. The 24 CR3BP systems covered (Table 1, p. 6)

Verbatim from Table 1 — secondary, primary, mu, distance(km), period(days),
velocity(km/s), q_max. (q_max is the largest q in the resonance R_{p:q}
searched for that system.)

| # | Secondary | Primary | mu          | a (km)      | P (days) | v (km/s) | q_max |
|---|-----------|---------|-------------|-------------|----------|----------|-------|
| 1 | Mercury   | Sun     | 1.6601209e-07 | 579092e+02 | 87.97    | 47.87    | 30    |
| 2 | Venus     | Sun     | 2.4478435e-06 | 108209e+03 | 224.70   | 35.02    | 20    |
| 3 | Earth     | Sun     | 3.0404317e-06 | 149598e+03 | 365.30   | 29.78    | 20    |
| 4 | Mars      | Sun     | 3.2271676e-06 | 227944e+03 | 686.99   | 24.13    | 30    |
| 5 | Jupiter   | Sun     | 9.5388609e-04 | 778341e+03 | 4332.70  | 13.06    | 10    |
| 6 | Saturn    | Sun     | 2.8580502e-04 | 142667e+04 | 10755.95 |  9.65    | 10    |
| 7 | Uranus    | Sun     | 4.3660686e-05 | 287066e+04 | 30702.53 |  6.80    | 15    |
| 8 | Neptune   | Sun     | 5.1511377e-05 | 449840e+04 | 60226.40 |  5.43    | 15    |
| 9 | Moon      | Earth   | 1.2150597e-02 | 384400e+00 | 27.28    |  1.02    | 10    |
|10 | Io        | Jupiter | 4.7042235e-05 | 421800e+00 |  1.77    | 17.33    | 15    |
|11 | Europa    | Jupiter | 2.5280092e-05 | 671100e+00 |  3.55    | 13.74    | 15    |
|12 | Ganymede  | Jupiter | 7.8043196e-05 | 107040e+01 |  7.15    | 10.88    | 15    |
|13 | Callisto  | Jupiter | 5.6666297e-05 | 188270e+01 | 16.69    |  8.20    | 15    |
|14 | Enceladus | Saturn  | 1.9011496e-07 | 238042e+00 |  1.37    | 12.62    | 30    |
|15 | Tethys    | Saturn  | 1.0864712e-06 | 294672e+00 |  1.89    | 11.35    | 20    |
|16 | Dione     | Saturn  | 1.9276021e-06 | 377415e+00 |  2.74    | 10.03    | 20    |
|17 | Rhea      | Saturn  | 4.0584367e-06 | 527068e+00 |  4.52    |  8.48    | 20    |
|18 | Titan     | Saturn  | 2.3663937e-04 | 122186e+01 | 15.95    |  5.57    | 10    |
|19 | Iapetus   | Saturn  | 3.1771370e-06 | 356085e+01 | 79.34    |  3.26    | 20    |
|20 | Ariel     | Uranus  | 1.4405060e-05 | 190900e+00 |  2.52    |  5.51    | 15    |
|21 | Umbriel   | Uranus  | 1.4686380e-05 | 266000e+00 |  4.14    |  4.67    | 15    |
|22 | Titania   | Uranus  | 3.9167599e-05 | 436300e+00 |  8.70    |  3.64    | 15    |
|23 | Oberon    | Uranus  | 3.5436226e-05 | 583500e+00 | 13.47    |  3.15    | 15    |
|24 | Triton    | Neptune | 2.0881945e-04 | 354759e+00 |  5.88    |  4.39    | 10    |

Data source quoted on p. 5: "These parameters were obtained from the JPL Solar
System Dynamics site" (footnote 1: `https://ssd.jpl.nasa.gov/?phys_data`).

### Cross-check against the project's CR3BP coverage

The project's `core/satellites.py` (per docs/notes/2026-06-08-moontour-tier1-
complete.md) registers the Galilean four (Io/Europa/Ganymede/Callisto), the
Saturnian midsize and Titan, and uses JPL-SSD-sourced mass parameters. The
project's catalogue spans Sun-planet for the major planets plus Earth-Moon
plus Sun-Jupiter-moon plus Saturn moons plus Uranus moons plus Pluto-Charon
(per the project memory).

**Coverage delta**:
- **In Restrepo-Russell but NOT in this project's registry**: Mercury,
  Mars (as separate Sun-Mars; we have a Sun-Mars system), Venus,
  Iapetus, Ariel, Umbriel, Titania, Oberon, Triton. The four Uranus
  moons + Triton are particularly relevant — the project's Uranus phase-4
  work in #332 (docs/notes/2026-06-16-332-v4-uranus-phase4.md) overlaps
  with the Ariel/Umbriel/Titania/Oberon coverage here.
- **In this project but NOT in Restrepo-Russell**: Pluto-Charon. The paper
  caps coverage at Neptune-Triton (Table 1 row 24).
- The mu values in Table 1 are an independent reproduction source for the
  project's mass-parameter table — useful cross-validation.

## 4. The differential corrector methodology

### 4.1 What Restrepo-Russell do (p. 4 §2.1, p. 7 §3.1)

The state vector is six-dim `x = [x y z xdot ydot zdot]` (Eq. 1 p. 3; the
paper carries the full 3D EOM but solves only the planar `z = zdot = 0`
restriction). First- and second-order variational equations (Eqs. 3-4, p. 4):

```
Phi^1_dot = (df/dx) Phi^1
Phi^2_dot = (df/dx) (x) Phi^2  +  Phi^1^T (x) (d^2 f / dx^2) (x) Phi^1
```

where `(x)` is the tensor product. `Phi^2` and `d^2 f / d x^2` are explicitly
described as 6 x 6 x 6 tensors. The paper cites **Pellegrini-Russell 2016**
(p. 4) for tensor-product details — which is the same JGCD STM-accuracy paper
this project's #372 fixed-path-mode work is grounded in.

The corrector targets the **half-period** (`T/2`) perpendicular-crossing
conditions
`xdot(T/2) approx 0` and (redundant but stored) `y(T/2) = 0`, with the
unknowns being `(ydot0, T/2)` (two unknowns, one residual after axisymmetry
collapses `y(T/2)`). The corrector form is **trust-region**, not Newton:
solves Eq. (7) p. 7

```
min phi(x) = g x + x^T H x   s.t. |x| <= rho
```

where `phi` is the squared-residual performance index, `g = del phi / del x`
and `H = del^2 phi / del x^2` (the FULL Hessian, computed from `Phi^2`, NOT
the Gauss-Newton approximation `g g^T`). The trust radius `rho` adapts each
iteration based on agreement between the quadratic model and the actual `phi`
reduction. The inner subproblem is solved by LU or eigenvalue decomposition
of `H` plus a 1D root solve per Conn-Gould-Toint (2000).

### 4.2 How this differs from the project's stack

The project's existing CR3BP corrector (per `src/cyclerfinder/core/cr3bp.py`,
which is in flight for #372 STM fixed-path mode against the
Pellegrini-Russell 2016 finding) uses a **first-order variational + variable-
step Gauss-Newton** approach, the standard JPL Howell-style multi-shooter that
this project inherits via the patched-conic + STM heritage. Differences:

1. **Order**: Restrepo-Russell use the full second-order Hessian assembled
   from `Phi^2`. The project's current corrector uses first-order only.
   The Hessian buys robustness on highly-sensitive seeds (connecting
   resonances near heteroclinic intersections), at the cost of one extra
   tensor propagation per orbit.
2. **Globalisation**: Restrepo-Russell use Conn-Gould-Toint trust-region.
   The project's current corrector uses Gauss-Newton with optional step-
   damping but no formal trust radius.
3. **Sensitivity envelope**: The trust-region method explicitly allows
   "coarser grid searches than are possible when using a first-order
   in-hone search method" (p. 7) — i.e., it expands the basin of
   attraction. Quantitatively the paper does not give a basin-radius
   number; this is a qualitative claim.

**Adoption recommendation**: The second-order trust-region corrector is a
plausible **complement** (not replacement) to the project's existing first-
order corrector, with the following caveats:
- The project's #372 fixed-path mode is the right place to land it: a
  fixed-path STM evaluation is exactly the dependency a trust-region inner
  solve needs, and #372 is already touching `cr3bp.py`.
- A `corrector_order="trust_region_2"` opt-in flag is the right shape —
  default first-order for backwards-compat, opt-in second-order for
  highly-sensitive seeds (resonances close to heteroclinic intersections,
  which is exactly the project's discovery-daemon target regime).
- Conn-Gould-Toint 2000 is a textbook citation — the algorithm is standard,
  not novel here. The novelty in Restrepo-Russell is assembling `H` from
  `Phi^2`. Pellegrini-Russell 2016's STM-accuracy bounds (which the
  project's #372 work is grounded in) extend naturally to the tensor case.

Honest negative: the project's current discovery regime is operating well
inside the first-order corrector's basin (no convergence failures attributed
to corrector order in the recent agent traffic). A second-order corrector
buys nothing measurable at the current frontier — the bottleneck is upstream
(seed quality, family selection, real-ephemeris closure per S1L1). Adoption
should be **deferred** unless and until the discovery daemon hits a wall
that traces to corrector-basin radius.

## 5. The orbit classes and nomenclature

The paper's orbit taxonomy (§4 p. 11, Table 3 p. 22) is structured as two main
groups, each with sub-classes:

### 5.1 Bound POs (BPOs, §4.1)
Orbits that remain in the secondary's vicinity, never leaving its sphere of
influence.

- **Simple BPOs** (`N_cross = 1`, first row of Fig. 4):
  - `LL1`, `LL2`: Lyapunov orbits around L1, L2 (the generating families)
  - `LPO1`, `LPO2`: Low Prograde Orbits — `g2` family stable region
  - `LPO1*`, `LPO2*`: extensions of LPO1/LPO2 into instability regions
  - `DPO`: Distant Prograde Orbit — unstable continuation of `g2` (per p. 12)
  - `DRO`: Distant Retrograde Orbit (Ocampo-Rosborough 1993) — stable
- **Composed BPOs** (`N_cross > 1`, second row of Fig. 4):
  - `Hg`: combination of L1 Lyapunov + DPO; heteroclinic-approx between
    LL1 and LL2 (the project's KEY family for L-point connections)
  - `Hb1`, `Hb2`: Lyapunov + LPO1 or LPO1+LPO2 composed structures
  - `Hg(b)`, `Hg(c)`: structural variants of Hg
  - `Hg1r`, `Hg2r`, `Hg3r`: multi-loop variants
- **Mushroom/multi-loop**: `M1`, `M2`, `Hm1`, `Hm2`, with `(b)`, `(c)`
  sub-variants
- **Quasi-DRO orbits**: `QDRO_{11}`, `QDRO_{14}`, `QDRO(g3)` — DRO
  perturbations with 11, 14 axis crossings

### 5.2 Resonant POs (RPOs, §4.2)
Notation: `R_{p:q}` where p = revolutions around the primary (inertial frame),
q = revolutions around the secondary (p. 15). Inner resonance: `p > q`,
prograde. Outer resonance: `p < q`, retrograde. Order = `|p - q|`. Examples
shown in Fig. 8 p. 15: `R_{3:2}`, `R_{7:5}`, `R_{8:5}`, `R_{5:4}`, `R_{3:4}`,
`R_{3:5}`. The paper notes (p. 15) `T_p approx q T_q` in inertial frame, hence
`T approx 2*pi*q` in normalised units.

### 5.3 Connecting resonances (§4.2 p. 17)
The **novelty**: orbits found by the local grid search around the LL1/LL2/DRO
generators that approximate heteroclinic connections between a simple
resonance and a simple BPO. Notation: `H_{R p:q-BPO}`. Specific instances:
- `H_{R p:q - LL1}`: connects R_{p:q} resonance with L1 Lyapunov
- `H_{R p:q - LL2}`: connects R_{p:q} resonance with L2 Lyapunov
- `H_{R p:q - DPO}`: connects R_{p:q} resonance with DPO
- `H_{R p:q - QDRO}`: connects R_{p:q} resonance with QDRO (Fig. 16 p. 20)

The paper explicitly highlights (p. 17 §4.2, p. 19 §4.3) that these connecting
resonances "model natural escape/capture transfers" with velocity discontin-
uities "of the order of a few meters per second" — meaning the patched-CR3BP
approximation cost is in the noise of real-system perturbations. This is the
artefact relevant to the project's S1L1 closure work and discovery-daemon
output.

### 5.4 What is explicitly EXCLUDED from the search (p. 3 §1)
- Asymmetric POs
- POs centered around L3, L4, L5
- 3D (out-of-plane) POs

### 5.5 Cross-reference with the project's v4.7 taxonomy

The project's `orbit_class` enum (per data/catalogue.schema.json) is
`{cycler, quasi_cycler, precursor_mga, mga_tour}` — a **trajectory-design
taxonomy** based on whether the orbit (a) is strictly periodic, (b) is epoch-
locked, and (c) inserts into a steady-state cycler.

Restrepo-Russell's taxonomy is **dynamical**, classifying by geometry of the
periodic orbit in the CR3BP rotating frame. The two taxonomies are **orthog-
onal**, not in conflict — a Restrepo-Russell `LL1` Lyapunov orbit is not a
cycler in the project's sense (it loiters around L1 and doesn't connect
inhabited bodies), but a Restrepo-Russell connecting-resonance `H_{R 2:1-LL1}`
in the Sun-Earth system can absolutely be a `cycler` or `quasi_cycler` row in
the project's catalogue if it has the right body-encounter pattern.

**Mapping**:
- Restrepo-Russell BPOs (LL1, LL2, Hg, DRO, QDRO, Hm, etc.) — almost all
  are NOT cyclers in the project sense. They loiter, don't cycle.
- Restrepo-Russell simple resonances `R_{p:q}` — most are NOT cyclers
  (they orbit only the primary, no secondary encounter). Some inner-
  resonance cases with close secondary approaches MIGHT be precursors.
- Restrepo-Russell connecting resonances `H_{R p:q - BPO}` —
  THIS is the subset that overlaps the project's cycler / mga_tour
  search space. These are the orbits that approximate heteroclinic
  capture/escape and would feed the project's free-return-arc and
  V-infinity-matched encounter framework.

So the cross-walk is: **add a new optional row field `cr3bp_family` (string)
that records the Restrepo-Russell taxonomy label** when a project catalogue
row originates from or matches the Restrepo-Russell database. This preserves
the project's `orbit_class` design taxonomy as primary and adds the dynamical
family as ancillary provenance.

## 6. Data accessibility (p. 20 §5)

The paper publishes TWO URLs in p. 20 footnotes 4 and 5:

- Footnote 4: `http://russell.ae.utexas.edu/index_files/POdatabase.htm`
  (the canonical landing page on Russell's UT Austin lab site)
- Footnote 5: `https://drive.google.com/drive/u/0/folders/0B7SdUc9xp3V6VTRGbmwzQ1lzSFE`
  (a Google Drive mirror)

**Format and structure** (p. 20-21 §5):
- 24 top-level folders, one per CR3BP system (Table 1 row order)
- Per-system subfolders for each individual search:
  - 1 low-resolution global search (for quick examination)
  - 1 high-resolution global search
  - 2 local searches around LL1 (one designed for BPOs, one for RPOs)
  - 2 local searches around LL2
  - 1 local search around DRO
- Data files are **plain text, column-arranged**, 31 parameters per solution
- File sizes: high-resolution global search file ~50 MB per system; other
  files 5-30 MB. Total: ~3 GB uncompressed (~600 MB compressed) (p. 21).
- A `readme` file in the top-level folder documents the column layout.
- Solutions count: "approximately 3 million periodic solutions" (Abstract,
  p. 1; confirmed in §6 conclusions p. 21).

**Generation provenance** (p. 21): GNU Fortran 4.4.7, `-O3`, Linux 24-core
2666 MHz machine, ~700 core-hours total.

### Project ingestion feasibility

Plain-text columnar Fortran-emitted files with a documented `readme` are
**straightforwardly machine-ingestible**. Python `numpy.loadtxt` or a pandas
reader handles this directly. The 3 GB total is well within disk budget. The
project's existing `data/*.jsonl` infrastructure (bcr4bp_halo_family_*.jsonl,
cr3bp_silver.jsonl, etc.) is the right destination shape: one JSON line per
orbit, fields aligned to the project's existing catalogue schema with new
optional fields for Restrepo-Russell's 31-column parameter set.

The **resonant subset** is the project-relevant ingest target — specifically
the `H_{R p:q-LL1}`, `H_{R p:q-LL2}`, `H_{R p:q-DPO}`, `H_{R p:q-QDRO}`
connecting-resonance families from the local searches around LL1/LL2/DRO.
The bound libration families (LL1/LL2/Hg/etc.) are not cycler candidates
and should be **excluded** from ingest as a first cut.

A rough Fermi-estimate sizing: of the ~3M total orbits, the global-search
fraction (resonances + bound) dominates (high-res global = 50 MB; local =
5-30 MB). If the connecting-resonance fraction is ~1-5% of the local-search
output across 24 systems with ~3 local searches each, the ingest target is
on the order of **10^4 to 10^5 candidate connecting-resonance orbits** —
sized to fit comfortably in the project's existing review-queue
infrastructure.

## 7. Catalogue impact

### 7.1 Concrete recommendation

The project's natural extension target inside Restrepo-Russell is the
**connecting-resonance subset** — `H_{R p:q-LL1}`, `H_{R p:q-LL2}`,
`H_{R p:q-DPO}`, `H_{R p:q-QDRO}` — across the 24 systems. The bound libration
families and pure resonances are NOT cycler material in the project sense and
should not be ingested wholesale (they would inflate the catalogue with non-
cycler dynamical curiosities).

### 7.2 Data fidelity tier

Restrepo-Russell ICs are CR3BP planar axisymmetric — exactly the fidelity
tier the project's `cr3bp_silver.jsonl` and #347 Phase 1 baseline operates
at. Schema mapping:
- `orbit_source`: new SOURCE_REGISTRY enum value
  `restrepo-russell-2018-cmda` (per the project's controlled-vocabulary
  enum, per v4.6 schema)
- `orbit_fidelity`: `Fidelity.CR3BP_PLANAR_AXISYMMETRIC` (new tier if not
  already present; per project memory the existing CR3BP tiers cover this
  case)
- `validation_level`: V0 on initial ingest (sourced IC, no project
  reproduction yet); promotes to V1 once the project's CR3BP integrator
  re-converges the orbit; V2+ requires the project's own corrector,
  stability, Jacobi-constant cross-check.

### 7.3 Schema accommodations needed

- New optional row field `cr3bp_family: str` recording the Restrepo-Russell
  taxonomy label (e.g., `"H_R7:5-LL1"`)
- New SOURCE_REGISTRY enum value `restrepo-russell-2018-cmda`
- New `data_gaps[].kind` value `cr3bp-planar-only` for rows that have not
  yet been promoted to 3D / ER3BP fidelity
- No `orbit_class` enum change needed — existing
  `{cycler, quasi_cycler, precursor_mga, mga_tour}` is orthogonal

### 7.4 Honest negative

Most of Restrepo-Russell's 3M orbits are **not cycler material**. The bound
libration families dominate (LL1, LL2, LPO1, LPO2, DPO, DRO, Hg, Hb, Hm,
QDRO and their variants) and the simple resonances (`R_{p:q}` without the
`H_*` prefix) are primary-only orbits with no secondary encounters. The
project's cycler-search-focused catalogue is correctly indifferent to those
classes. Per the user's memory `project_negative_results_registry`, the
sourced-cycler bar is high; ingesting bound libration families would inflate
the catalogue with rows that fail V2 (the cycler-definition gate) by
construction.

## 8. KNOWN_CORPUS impact

**Recommendation**: Add ONE new `CorpusAnchor` to
`src/cyclerfinder/search/literature_check.py` with the following shape
(do NOT edit — recommend only):

```python
CorpusAnchor(
    name="Restrepo-Russell JPL planar axisymmetric CR3BP database",
    primary="any",  # 24 CR3BP systems; see body_set
    body_set=frozenset({  # the 24 secondaries
        "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn",
        "Uranus", "Neptune", "Moon", "Io", "Europa", "Ganymede",
        "Callisto", "Enceladus", "Tethys", "Dione", "Rhea", "Titan",
        "Iapetus", "Ariel", "Umbriel", "Titania", "Oberon", "Triton",
    }),
    authors=("Restrepo", "Russell"),
    keywords=(
        "CR3BP periodic orbit database",
        "planar axisymmetric",
        "connecting resonance",
        "heteroclinic",
        "Lyapunov LL1 LL2",
        "distant retrograde orbit",
        "DRO",
        "QDRO",
        "Hg Hb Hm",
    ),
    citation=(
        "Restrepo, R. L. & Russell, R. P., "
        "'A database of planar axisymmetric periodic orbits for the "
        "Solar system,' Celest. Mech. Dyn. Astron. 130:49 (2018), "
        "DOI 10.1007/s10569-018-9844-6; database online at "
        "russell.ae.utexas.edu/index_files/POdatabase.htm"
    ),
    doi="10.1007/s10569-018-9844-6",
    topology_label=frozenset(),  # spans all topologies; no restriction
    period_band_tu=None,         # spans all energies
)
```

This is a **methodology + database-pointer anchor**, NOT a per-orbit anchor.
The 3M orbits are not individually citable; the database is the citable
object. A discovery-daemon candidate that hits this anchor's structural
fingerprint (any of the 24 systems, any of the Restrepo-Russell family
labels) should be returned `status="published"` with the citation above —
NOT certified novel.

The CRITICAL anchor scope: this anchor should match **CR3BP planar
axisymmetric** candidates only. A 3D / ER3BP / BCR4BP candidate is OUT of
Restrepo-Russell scope (the paper explicitly excludes 3D POs, p. 3 §1).
Per the `topology_label` mechanism in `CorpusAnchor`, this means leaving
`topology_label=frozenset()` (no restriction) is wrong — it should encode
`{"planar", "axisymmetric"}` if those labels exist in the project's
controlled vocabulary, or a new label set if they don't. A
`data_gaps[].kind` cross-walk also applies. The right shape is for the
human (parent) to decide on commit.

## 9. #347 Phase 2/3 implications

The project's #347 plan is structured as Phase 1 (current CR3BP baseline)
-> Phase 2 (discovery sweep) -> Phase 3 (higher-fidelity promotion).

**Phase 2 use**: The Restrepo-Russell database is an **excellent published-
orbit comparison baseline** for the discovery sweep. Workflow:
1. Discovery daemon emits a SILVER candidate in (say) Sun-Mars at a given
   Jacobi constant and family label.
2. Cross-walk the daemon's structural fingerprint against the Restrepo-
   Russell column-file for Sun-Mars at matching J.
3. If a Restrepo-Russell orbit matches within numerical tolerance ->
   `literature_check.status="published"` with the database pointer
4. If no match -> the candidate is in-scope-but-not-found (a Restrepo-
   Russell-blessed cleared region); proceed to V1+ validation.

This is exactly what the existing `literature_check.py` is designed for,
but with the database as a **machine-queryable** corpus rather than a
web-search corpus. The right shape is a new ingest path:
`scripts/ingest_restrepo_russell_database.py` that pulls the column files,
parses the per-system tables, and emits a JSONL index keyed by
`(system, family_label, J_band, period_band)` for fast structural lookup.

**Phase 3 use**: Once a discovery candidate clears CR3BP planar-axisymmetric
fidelity (V1-V2), promotion to higher fidelity (ER3BP, BCR4BP, ephemeris)
needs the Restrepo-Russell IC as a seed. The database's `(x0, ydot0, T/2)`
plus `Phi^1`/`Phi^2` STT framework is the right input to the project's
multi-shooting transcription pipeline. Cite Restrepo-Russell as the seed
provenance.

Honest assessment: the Restrepo-Russell database is more valuable to the
project as a **discovery-daemon cross-check corpus** than as a wholesale
catalogue ingest. Most of the 3M orbits never become catalogue rows; what
they do is **certify that the discovery daemon's not-found verdict is
trustworthy in the planar axisymmetric CR3BP envelope**. This dramatically
strengthens the project's literature-novelty gate (`feedback_literature_
novelty_check_baseline` memory) for the 24 covered systems.

## 10. Errata vs the project's methodology

I found **no methodological contradictions** between Restrepo-Russell 2018
and the project's current stack. The paper's CR3BP normalisation (Eq. 1-2
p. 3-4), Jacobi-constant definition (Eq. 5 p. 4), and equilibrium-point
labelling (L1 between primary and secondary, L2 outside secondary, L3 on
the opposite side of the primary, L4/L5 triangular) all align with the
project's conventions.

One observation that **confirms** rather than contradicts:
- The paper cites **Pellegrini-Russell 2016** (p. 4) for tensor-product
  details on `Phi^2` propagation. This is the same paper the project's
  #372 STM fixed-path mode is grounded in. Restrepo-Russell using `Phi^2`
  for the Hessian is internally consistent with the project's #372 finding
  that variational+VS path mismatches can corrupt STM accuracy — both
  papers depend on the STT being accurate, and the project's #372 fix is
  precisely the kind of correctness gate Restrepo-Russell would have needed
  in their implementation (they don't say how they handle it; presumably
  fine-grained tolerances on the integrator hide the issue).

One naming **near-conflict** worth noting (not an error, a vocabulary clash):
- Restrepo-Russell use `DRO` (Distant Retrograde Orbit) per Ocampo-
  Rosborough 1993, a planar L1/L2 retrograde libration. The project's
  use of "retrograde" in (e.g.) the Morais-Namouni 2013 context refers
  to retrograde resonances in the inertial frame around the Sun. These
  are NOT the same object. When ingesting Restrepo-Russell DRO families,
  the `cr3bp_family` field should carry the literal `"DRO"` label and not
  be confused with the project's existing retrograde-resonance language.

## 11. Action items for the parent

1. **Add #116 wishlist closure**: mark Restrepo-Russell 2018 as ACQUIRED
   (9 of 9 complete); update OUTSTANDING.md / wishlist accordingly.

2. **Recommend (not implement) the KNOWN_CORPUS anchor** in §8 — the parent
   should review the `topology_label` scope (planar+axisymmetric) and the
   `body_set` cross-walk against the project's body-naming convention
   (the paper uses "Earth"/"Moon"/"Europa"; the project uses code letters
   `E`, `M`, etc.).

3. **Defer the database ingest** until the discovery daemon needs it. The
   value of Restrepo-Russell to the project is as a **literature-novelty
   cross-check corpus**, not as a wholesale catalogue ingest. The right
   trigger is the first discovery-daemon SILVER candidate in a
   Restrepo-Russell-covered system at a Jacobi constant inside the
   database's energy envelope. At that point, a focused ingest of the
   single relevant system file (50 MB) is fast and trivial.

4. **Defer the second-order trust-region corrector** unless and until
   #372 (STM fixed-path mode) lands and the discovery daemon shows
   corrector-basin-radius failures. The second-order corrector is a
   legitimate adoption candidate but not on the critical path now.

5. **Cross-walk against #332 Uranus phase-4 work**: Restrepo-Russell
   covers Uranus-Ariel/Umbriel/Titania/Oberon. The #332 axisymmetry note
   (J4 ~ -2.9e-5) is relevant — Restrepo-Russell's planar axisymmetric
   restriction is **dynamically consistent** with the project's #332
   choice to treat Uranus as closely axisymmetric. Worth a one-line
   cross-reference in the #332 follow-up.

6. **Do NOT plug into the Phase 2 discovery sweep yet**: per the project's
   `feedback_speculative_high_effort_required` memory, the right next
   discovery-frontier move is ER3BP / BCR4BP / 3D / QP / epoch — all of
   which are **outside** Restrepo-Russell's planar axisymmetric scope.
   The database is a CONFIRMATION corpus for the planar CR3BP envelope
   (which the project memory `project_validation_ceiling` records as
   substantially mined), not a frontier-extension corpus.

7. **Schema additions** (recommend, don't implement):
   - `cr3bp_family: str` (optional row field)
   - SOURCE_REGISTRY enum `restrepo-russell-2018-cmda`
   - `data_gaps[].kind = "cr3bp-planar-only"`

8. **Update `feedback_literature_novelty_check_baseline` memory** (or
   leave to the human): the Restrepo-Russell database is the single
   largest published CR3BP corpus the project has access to. The
   literature-novelty gate should treat a Restrepo-Russell match as a
   FIRM "published" verdict (not a "maybe").

---

## Reference appendix: Restrepo-Russell paper bibliography (pp. 22-24)

Key references the project should be aware of (most already known):
- Conn, A.R., Gould, N.I.M., Toint, P.L. (2000), *Trust-Region Methods* —
  the algorithmic citation for the corrector globaliser (§3.1 p. 7)
- Pellegrini, E. & Russell, R.P. (2016), "On the accuracy of state
  transition matrices," JGCD 39(11):2485-2499 — the STM-accuracy paper
  this project's #372 work is grounded in (cited p. 4)
- Russell, R.P. (2006), "Global search for planar and three-dimensional
  periodic orbits near Europa," J. Astronaut. Sci. 54(2):199-226 — the
  grid-search precursor paper (§3.1 p. 7)
- Henon, M. (2003) "New families of periodic orbits in Hill's problem of
  three bodies," CMDA 85(3):223-246 — the Hb/Hg family nomenclature
  origin (footnote 3, p. 11)
- Ocampo, C.A. & Rosborough, G.W. (1993), "Transfer trajectories for
  distant retrograde orbiters of the Earth" — the DRO term origin
  (p. 8 caption, p. 12 §4.1)
- Lara, M. & Russell, R.P. (2007), "On the family g of the restricted
  three-body problem," Monogr. Real Acad. Cienc. Zaragoza 30:51-66 —
  prograde "g family" nomenclature (p. 11 §4.1)

End of digest.
