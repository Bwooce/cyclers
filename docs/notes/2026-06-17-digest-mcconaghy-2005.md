# McConaghy-Russell-Longuski 2005 (JSR 42(4)) — full-page deep-read verdict

**Read 2026-06-17 AET.** Part of Mars-cycler digest wave (Agent A).

## 1. Header (verbatim from PDF p.1)

- Title: **"Toward a Standard Nomenclature for Earth-Mars Cycler Trajectories"**
- Authors: **T. Troy McConaghy* (Purdue University, West Lafayette IN
  47907-2023), Ryan P. Russell† (University of Texas at Austin, Austin TX
  78712-1085), James M. Longuski‡ (Purdue University, West Lafayette IN
  47907-2023)**.
  Footnotes: * Ph.D. Candidate, School of Aeronautics and Astronautics
  (mcconagh@ecn.purdue.edu), Student Member AIAA; † Ph.D. Candidate, Dept. of
  Aerospace Engineering and Engineering Mechanics
  (ryanrussell@mail.utexas.edu), Student Member AIAA; ‡ Professor, School of
  Aeronautics and Astronautics (longuski@ecn.purdue.edu), Associate Fellow
  AIAA.
- Venue: **Journal of Spacecraft and Rockets, Vol. 42, No. 4, July-August
  2005, pp. 694-698**.
- DOI: **10.2514/1.8123**.
- Received 6 February 2004; revision received 4 June 2004; accepted 14 June 2004.
- Associate Editor: C. Kluever.
- Length: **5 pages** (pp. 694-698).

## 2. What the paper actually is

The **standard nomenclature** paper for Earth-Mars cyclers. McConaghy 2004
JSR (paper 4 in this wave) introduced the `nPr` family-tag naming; this paper
introduces the **complete per-leg descriptor** convention. Every cycler can be
labelled by:

```
[(<body-sequence>)] n d₁d₂...d_K
```

where n = repeat time in synodic periods, K = number of Earth-Earth transfer
legs, and each dᵢ is a leg descriptor of one of three forms:

- **Generic transfer**: `g(t_f, θ, ε)` — arbitrary angular separation
- **Full-revolution transfer**: `f(M:N, ϕ, λ)` — Earth returns to same point
- **Half-revolution transfer**: `h(t_f, N, ε, i')` — Earth on opposite side

Where:
- `t_f` = transfer time (years).
- `θ` = transfer angle (radians, > 2π for multi-rev).
- `ε ∈ {U, L, Lₗ, Lₛ}` = Lagrange Lambert solution type (the U / L distinguish
  upper / lower branches; Lₗ Lₛ subscripts discriminate long/short period on
  the lower branch when two solutions exist).
- `M:N` resonance ratio (M Earth revs ≡ N s/c revs).
- `ϕ` (latitude) and `λ` (longitude) angles of v_out on full-rev transfers
  (Fig. 3, p.696).
- `i'` = signed half-rev orbit plane orientation angle.

A formal **Extended Backus-Naur Form** grammar is given (Table 3, p.698):

```
<cycler label>     ::= [(<body sequence>)] n <leg sequence>
<leg sequence>     ::= <leg descriptor> | <leg descriptor><leg sequence>
<leg descriptor>   ::= g(t_f, θ, ε)[^k] | f(M:N, ϕ, λ)[^k] | h(t_f, N, ε, i')[^k]
<body sequence>    ::= <body abbreviation><body sequence> | <body abbreviation><body abbreviation>
<body abbreviation>::= E | M | V | ...
```

Cycler repeat-time constraint (Eq. 28, p.697): `t_f1 + t_f2 + ... + t_fK = nS`
where S = Earth-Mars synodic period (yr).

### Title / venue / DOI cross-check against our records

| Field | Our records | Paper actual | Status |
|---|---|---|---|
| KNOWN_CORPUS | "Russell-Ocampo / McConaghy SnLm cyclers" anchor (line 263) cites McConaghy 2006 JSR 43(2) but NOT this 2005 paper. | This is the **canonical SnLm nomenclature paper** that named SnLm formally. | The anchor uses the SnLm name without citing the paper that defines it. **Citation gap.** |
| Wave digest 2026-06-17 entry | "Toward a Standard Nomenclature for Earth-Mars Cyclers — the SnLm naming convention paper." | Title is correct; but **the paper does NOT introduce "SnLm" as such**. SnLm appears as one of the four ε values (`{U, L, Lₗ, Lₛ}`), distinct from the **S/L** of SnLm in the catalogue context. The catalogue's "SnLm" likely conflates the L Lambert-solution-type with a count of Lambert legs. **MISMATCH** — see §6. |
| Filename | "mcconaghy-russell-longuski-2005-standard-nomenclature-earth-mars-cycler-trajectories-jsr-doi-10.2514-1.8123.pdf" | matches | OK |

**CRITICAL FINDING**: This paper does NOT actually introduce "SnLm" nomenclature.
It introduces a generic `n d₁d₂...d_K` label where each d is one of `g`, `f`,
`h`. The S/L in this paper are **Lambert-solution subtype tags** (`Lₛ`, `Lₗ`),
NOT "S-leg" / "L-leg" identifiers. **The "S1L1" of e.g. Russell-Ocampo's
dissertation and our catalogue is from a DIFFERENT (Russell-school)
nomenclature** that this paper does not adopt. See errata §6.

## 3. Key numerical / structural content

### Modeling assumptions (p.694-695)
1. Earth + Mars in circular coplanar orbits.
2. Only Earth and Sun are gravitating (Mars is encountered but doesn't change
   the s/c orbit).
3. All legs are conic.
4. Gravity-assists are instantaneous (impulsive Δv).
5. All transfers direct (no retrograde).
6. Propulsive maneuvers only at Earth encounters, not in deep space.

Note: "**we do not assume that the Earth-Mars synodic period S is exactly
2-1/7 years (or any other approximation)**" (p.695) — i.e. the nomenclature
is independent of the rational synodic approximation. The catalogue's
T_SYN=2.135 yr is consistent with this.

### Lambert solver (Eqs 1-8, p.696)
Standard multi-rev Lagrange Lambert:
- Eq (1): `√μ·t_f = a^(3/2)[2πN + α - β - sin(α) + sin(β)]`
- Eq (2): `N = floor(θ/2π)` — whole spacecraft revolutions.
- Eqs (3-4): `α = α₀ if ε ∈ {L, Lₗ, Lₛ}, 2π - α₀ if ε = U`;
  `β = β₀ if 0 < (θ mod 2π) < π, -β₀ if π < (θ mod 2π) < 2π`.
- Eqs (5-6): `α₀ = 2 arcsin √(s/2a)`, `β₀ = 2 arcsin √((s-c)/2a)`.
- Eqs (7-8): `s = (r₁+r₂+c)/2`, `c = √(r₁² + r₂² - 2r₁r₂ cos θ)`.

### Outgoing/Incoming velocity (Eqs 9-12, p.696)
- Eq (9): `v_out = [(B+A)/c](r₂-r₁) + [(B-A)/r₁]r₁`
- Eq (10-11): `A = √(μ/4a) cot(α/2)`, `B = √(μ/4a) cot(β/2)`
- Eq (12): `v_in = [(B+A)/c](r₂-r₁) - [(B-A)/r₂]r₂`

### Full-revolution (resonant) transfers (Eqs 13-16, p.696)
- Eq (13): `M · 2π√(a_E³/μ) = N · 2π√(a³/μ)` — resonance condition.
- Eq (14): `a = a_E (M/N)^(2/3)` — semi-major axis from resonance.
- Constraint: a > a_E/2 (so M, N constrained).
- Eq (15): `v_out = √(μ(2/r_E - 1/a))` — vis-viva magnitude.
- Eq (16): `v_out = v_out[(cos ϕ cos λ)(v̂_E × ĥ_E) - (cos ϕ sin λ)ĥ_E + (sin ϕ)v̂_E]`
  — direction in rotating frame.
  ϕ ∈ [-90°, 90°] = latitude; λ unconstrained = longitude.

### Half-revolution transfers (Eqs 17-27, p.697)
- Eq (17): `v_out = ṙ_out r̂_E + (κ_out sin i') ĥ_E + (κ_out cos i') θ̂_E`
- Eq (18): `ṙ_out = √(μ[2/(r₁+r₂) - 1/a]) · {-1 if ε ∈ {L, Lₗ, Lₛ}; +1 if ε = U}`
- Eq (19): `κ_out = √(2μr₂/(r₁² + r₁r₂))`
- Eqs (20-21): `r̂_E = (cos γ_E)(v̂_E × ĥ_E) + (sin γ_E)v̂_E`,
  `θ̂_E = (-sin γ_E)(v̂_E × ĥ_E) + (cos γ_E)v̂_E`
- Eq (22-27): Similar for v_in at the arrival end.

**Backflip definition** (p.696): half-rev transfer with N=0 following a
gravity-assist. i' > 0 = northern backflip; i' < 0 = southern backflip.

### Table 2 (p.697) — Well-known cyclers in this nomenclature

This is the **Rosetta stone** for catalogue cross-references.

| Common name | Refs | Cycler label |
|---|---|---|
| **Aldrin cycler** | 3-5 | **`1g(2-1/7, 1-1/7 rev, L)`** |
| **VISIT 1 cycler** | 4, 6-9 | **`7f(5:4, ϕ, 0 deg)³`** or **`7f(5:4, ϕ, 180 deg)³`** |
| **VISIT 2 cycler** | 4, 6-9 | **`7f(3:2, ϕ, 0 deg)⁵`** or **`7f(3:2, ϕ, 180 deg)⁵`** |
| **Ballistic S1L1 cycler** | 10, 14 | **`2g(2.8277, 657.97 deg, U) g(1.4508, 522.29 deg, L)`** |
| **Byrnes' case 3 cycler** | 11 | `2g(2-11/14, 1-11/14 rev, U) f(1:1, 79.612 deg, λ) h(0.5, 0, U, ±10.388 deg)` |
| **Cycler 2.5.1.+0** | 12 | `2g(1-11/14, 11/14 rev, U) f(1:1, 74.919 deg, ∓144.069 deg) h(0.5, 0, U, ±15.081 deg) f(1:1, 74.919 deg, ±35.931 deg)` |
| **Cycler 4.3.1.−5** | 12 | `4g(7-1/14, 5-1/14 rev, L) f(1:1, 84.039 deg, ∓90.0 deg) h(0.5, 0, U, ±5.961 deg)` |

**These are V1-grade reproducible state IN THE CIRCULAR-COPLANAR MODEL.**
Each label completely specifies the cycler geometry given the model
assumptions (p.698: "the cycler label completely determines all
characteristics of the cycler").

### Critical numerical anchors from Table 2

- **Aldrin = `1g(2-1/7, 1-1/7 rev, L)`**: ONE generic leg, t_f = 2-1/7 yr =
  2.1429 yr, θ = 1-1/7 rev = 8π/7 rad = 411.43°, ε = L (lower branch). Matches
  the Byrnes 1993 derivation.
- **Ballistic S1L1 = `2g(2.8277, 657.97 deg, U) g(1.4508, 522.29 deg, L)`**:
  - n = 2 (two synodic periods to repeat)
  - **Two generic legs**: leg 1 = (t_f=2.8277 yr, θ=657.97°, ε=U);
    leg 2 = (t_f=1.4508 yr, θ=522.29°, ε=L).
  - Total t_f = 2.8277 + 1.4508 = 4.2785 yr ≈ 2·S = 2·(2-1/7) = 4.286 yr ✓.
  - **657.97° = 360° + 297.97° → just under 2 revolutions, falls between 1π
    and 3π**.
  - **522.29° = 360° + 162.29° → between π and 2π, multi-rev** → catalogue
    memory `project_s1l1_realeph_closure_blocker` mentions S1L1 is two-arc;
    this CONFIRMS it via the McConaghy/Russell formal label.
- **Cycler 2.5.1.+0** has FOUR legs (g, f, h, f) and TWO body sequences (the
  Russell-Ocampo IDs like `2.5.1.+0` use the Russell dissertation
  `p-h-s-±i` convention; "Cycler 2.5.1.+0" maps to `Cycler-2-5-1-0` in
  Russell-Ocampo 2003 Table 4 with i=0).
- **Cycler 4.3.1.−5** = `4g(7-1/14, ...)` — note p=4 synodic periods, single
  big generic leg + one full-rev + one half-rev. The i=−5 in the
  Russell-Ocampo notation appears nowhere in the label — confirming Russell's
  notation `p.h.s.±i` is INDEPENDENT of the McConaghy formal label.

### Evaluating cyclers (p.698)
Three criteria (consistent with Russell-Ocampo 2003):
1. Number of vehicles for short Earth-Mars and Mars-Earth legs every synodic.
2. Per-leg transit time AND v∞ at both ends should be small.
3. Total Δv should be zero (ballistic).

### Reference list (p.698)
- Ref 3 = Aldrin 1985 SAIC presentation.
- Ref 4 = Friedlander-Niehoff-Byrnes-Longuski AIAA-86-2009 (Aug 1986).
- Ref 5 = Byrnes-Longuski-Aldrin JSR 30(3) 1993 (DOI 10.2514/3.25519).
- Ref 6-9 = Niehoff (multiple VISIT-related publications, including Niehoff
  AAS 86-172 = ref 8).
- Ref 10 = McConaghy-Longuski-Byrnes JSR 41(4) 2004 (DOI 10.2514/1.11939) —
  the nPr paper, this wave's paper 4.
- Ref 11 = Byrnes-McConaghy-Longuski AIAA 2002-4420 (Aug 2002).
- Ref 12 = Russell-Ocampo JGCD 27(3) 2004 (= the paper that this wave's paper 1
  is the preprint of).
- Ref 13 = Russell-Ocampo AAS 03-508 (Aug 2003) — different from AAS-03-145.
- Ref 14 = McConaghy-Yam-Landau-Longuski AAS 03-509 (Aug 2003) —
  "Two-Synodic-Period Earth-Mars Cyclers with Intermediate Earth Encounter"
  ← this is the **S1L1 cycler primary source per ref 10+14 in Table 2**.
- Ref 17-18 = Hénon (Generating Families in the Restricted Three-Body Problem).

## 4. `KNOWN_CORPUS` impact (RECOMMEND, do not edit)

The current SnLm anchor (literature_check.py line 263):

```python
name="Russell-Ocampo / McConaghy Earth-Mars SnLm cyclers",
authors=("Russell", "Ocampo", "McConaghy", "Landau", "Longuski"),
citation="Russell & Ocampo, J. Spacecraft & Rockets 41(1) 2004; "
"McConaghy et al., J. Spacecraft & Rockets 43(2) 2006",
```

**MISSING CITATION**: McConaghy-Russell-Longuski 2005 JSR 42(4) DOI
10.2514/1.8123 — the formal nomenclature definition paper.

**Author tuple update**: missing **Russell** as a separate author of this 2005
JSR paper (he's already in the tuple, but in his Ocampo-co-author role; he is
ALSO co-author of the 2005 nomenclature paper without Ocampo).

Recommended consolidated SnLm anchor (combining all McConaghy / Russell
findings from this wave's papers 1, 4, 5, 6):

```python
name="Russell-Ocampo / McConaghy Earth-Mars cyclers (multi-source nomenclature anchor)",
authors=("Russell", "Ocampo", "McConaghy", "Longuski", "Byrnes", "Landau", "Yam"),
citation=(
    "Russell & Ocampo, JGCD 27(3):321-335 (2004) DOI 10.2514/1.1909, "
    "preprint AAS-03-145 (Russell-Ocampo p.h.s.±i family enumeration); "
    "McConaghy, Longuski & Byrnes, JSR 41(4):622-628 (2004) DOI "
    "10.2514/1.11939 (nPr family-tag + S1L1 DE405 itinerary Table 6); "
    "McConaghy, Russell & Longuski, JSR 42(4):694-698 (2005) DOI "
    "10.2514/1.8123 (formal per-leg nomenclature: g/f/h descriptors, "
    "EBNF grammar); "
    "McConaghy, Landau, Yam & Longuski, JSR 43(2) (2006) DOI 10.2514/1.15215; "
    "Russell PhD diss (UT Austin 2004) Table 3.4 for signed-i catalogue "
    "IDs; McConaghy PhD diss (Purdue 2004)"
),
```

## 5. Catalogue impact (RECOMMEND)

### Rows touched

ALL `russell-ocampo-*` rows (~190 catalogue entries) and any `mcconaghy-*` rows
use family-tag IDs of the form `p.h.s.±i` (Russell-Ocampo dissertation) or
`SnLm` (informal Russell-Ocampo / catalogue convention). **None of these IDs
match the McConaghy 2005 formal leg-descriptor labels** like
`2g(2.8277°, 657.97°, U) g(1.4508°, 522.29°, L)`.

### V0→V1 promotion candidacy

**STRONG V1 CANDIDATES for the 7 Table 2 entries**: each formal label is a
complete circular-coplanar model specification.

For our catalogue, the actionable rows are:

1. **`aldrin-cycler`** (V3) — formal label `1g(2-1/7, 1-1/7 rev, L)` is fully
   reproducible in our Lambert + circular-coplanar stack. **No promotion
   needed** (already V3).
2. **A catalogue row for "Ballistic S1L1" with formal label
   `2g(2.8277, 657.97°, U) g(1.4508, 522.29°, L)`** — if exists, **V0 →
   V1 promotion candidate**. Memory `project_s1l1_realeph_closure_blocker`
   says we have an `s1l1` catalogue row that fails real-eph closure at the
   family-selection step. **This paper gives us the circular-coplanar gold
   standard to compare against** — if we reproduce the (t_f, θ, ε) leg pair
   in our Lambert solver to machine precision, that's V1 in this model.
3. **`niehoff-visit-1` / `niehoff-visit-2`** rows (if exist) — formal labels
   `7f(5:4, ϕ, 0/180 deg)³` and `7f(3:2, ϕ, 0/180 deg)⁵` give the resonance
   ratios M:N and λ values; ϕ is a free design variable. **V0 → V1 promotion
   candidate** in the model where ϕ is fixed to a published value (Niehoff's
   refs 6-9 — need to check those papers separately).
4. **Catalogue row for Byrnes' case 3** (if exists) — formal label
   `2g(2-11/14, 1-11/14 rev, U) f(1:1, 79.612°, λ) h(0.5, 0, U, ±10.388°)` —
   V0 → V1 candidate (one free var: λ).
5. **`russell-ocampo-2.5.1+0`** (catalogue line 1366) — formal label
   `2g(1-11/14, 11/14 rev, U) f(1:1, 74.919°, ∓144.069°) h(0.5, 0, U,
   ±15.081°) f(1:1, 74.919°, ±35.931°)`. **No free variables** — fully
   reproducible. **V1 PROMOTION CANDIDATE**.
6. **`russell-ocampo-4.3.1-5`** (catalogue line 1017) — formal label
   `4g(7-1/14, 5-1/14 rev, L) f(1:1, 84.039°, ∓90.0°) h(0.5, 0, U,
   ±5.961°)`. **No free variables.** **V1 PROMOTION CANDIDATE**.

The Russell-Ocampo 2003 verdict note noted 4 V1 candidates from that paper
(2-5-1-3, 3-1-2-11, 4-3-1-20, 4-5-2-12); this paper adds **at least 2 more
V1 candidates from Table 2** (2.5.1.+0 and 4.3.1.-5), with formal-label state
to compare against — distinct from Russell-Ocampo 2003 Tables 5-8 which gave
per-encounter Δv vectors.

**Critically**, Russell-Ocampo's per-encounter Δv vectors (Tables 5-8 of
paper 1) and McConaghy-Russell-Longuski 2005's formal leg labels (Table 2 of
paper 5) are TWO INDEPENDENT representations of the same family in the same
circular-coplanar model. If we reproduce BOTH for the same cycler (e.g.
2.5.1.+0), that's a STRONG cross-validation and the V1 promotion becomes
DEFENSIBLE.

## 6. Errata / surprises

1. **"SnLm" is NOT this paper's nomenclature** — this paper introduces the
   per-leg `g/f/h` descriptors, not "SnLm". The "SnLm" name used in the
   catalogue and KNOWN_CORPUS comes from somewhere else (perhaps Russell's
   dissertation Ch.3 or McConaghy's dissertation; need to verify). **The S/L
   letters in this paper's ε ∈ {U, L, Lₛ, Lₗ} are Lambert solution subtypes,
   not S/L "leg counts" as the SnLm convention implies.** This is a
   significant naming-conflict between this paper and our convention.
2. **The formal Ballistic S1L1 label `2g(2.8277, 657.97°, U) g(1.4508,
   522.29°, L)` has TWO arcs**. This confirms the
   `project_s1l1_realeph_closure_blocker` memory's identification of S1L1 as
   "actually MULTI-ARC (two generic-return arcs)". The first arc has θ ≈ 1.83
   rev, the second has θ ≈ 1.45 rev. **The two arcs share neither (a, e) nor
   period; they are JOINED at Earth gravity-assists.**
3. **The Russell-Ocampo `p.h.s.±i` family-tag convention (used in our
   catalogue) is INDEPENDENT** of the McConaghy 2005 per-leg formal label.
   Both label the same trajectories but with different information. The
   family tag tells you WHICH ENUMERATION CLASS; the formal label tells you
   the LEG GEOMETRY.
4. **VISIT cyclers have a free design variable ϕ** in their formal labels
   (`7f(5:4, ϕ, 0 deg)³`). Without ϕ pinned by a separate publication, we
   cannot fully reproduce a VISIT cycler — ϕ is a degree of freedom.
   v∞ depends on ϕ (Eq 16) so different ϕ → different cycler.
5. **The Aldrin label `1g(2-1/7, 1-1/7 rev, L)` rev = 1-1/7 = 8/7 ≈ 1.143
   revs ≈ 411.43°** is the heliocentric transfer angle — i.e. between one
   and two revs (consistent with Byrnes 1993 multi-rev Lambert solve and
   Russell-Ocampo 2003 N_MAX=1+).
6. **The body-sequence prefix is OPTIONAL** (square brackets in EBNF). When
   absent, all encounters are Earth (Earth-Earth transfers only). With Mars
   on the path, the prefix `(EM...)` could be added — but Table 2's entries
   omit this prefix consistently, even though all are Earth-Mars cyclers
   that encounter Mars. So **the absence of an M body abbreviation does NOT
   mean Mars is uninvolved** — it means Mars encounters don't change the
   orbit (assumption 2). This is a subtle but easy-to-misread aspect of the
   convention.
7. **n is the cycler repeat in SYNODIC PERIODS**, not number of legs. Aldrin
   has n=1 and K=1 (one leg). S1L1 has n=2 and K=2 (two legs). The two
   numbers can coincidentally agree but generally differ.

## 7. Action items for parent

- [ ] **Fix `KNOWN_CORPUS` SnLm anchor**: add McConaghy-Russell-Longuski 2005
  JSR DOI 10.2514/1.8123 as a primary citation. Combine with the Russell-Ocampo
  2003 (JGCD not JSR) fix and the McConaghy 2004 JSR addition; one consolidated
  multi-paper anchor.
- [ ] **Resolve the "SnLm" naming-convention question**: where does "SnLm"
  (S-leg / L-leg counts, NOT this paper's Lambert U/L/Lₛ/Lₗ) actually come
  from? Verify by checking the McConaghy 2004 PhD dissertation or the
  Russell 2004 dissertation Ch.3 (per the McConaghy 2004 dissertation mining
  note). The catalogue's "SnLm" convention may be an internal abbreviation,
  not a published nomenclature. Document in spec.md if so.
- [ ] **V0 → V1 promotion candidates from this paper's Table 2**:
  - **`russell-ocampo-2.5.1+0`** (catalogue line 1366) — formal label
    `2g(1-11/14, 11/14 rev, U) f(1:1, 74.919°, ∓144.069°) h(0.5, 0, U,
    ±15.081°) f(1:1, 74.919°, ±35.931°)`. Same-model reproduction is a
    Lambert + closed-form orbit propagation. Track as a follow-on.
  - **`russell-ocampo-4.3.1-5`** (catalogue line 1017) — formal label
    `4g(7-1/14, 5-1/14 rev, L) f(1:1, 84.039°, ∓90.0°) h(0.5, 0, U,
    ±5.961°)`. Same. Track as a follow-on.
- [ ] **Cross-check Ballistic S1L1 formal label against Tables 5-8 of
  Russell-Ocampo 2003** (if any of those tables list a 2-synodic cycler that
  matches `2g(2.8277, 657.97°, U) g(1.4508, 522.29°, L)`). Russell-Ocampo
  Table 5 is `Cycler-2-5-1-3` (`2-5-1-3` in unsigned i). If the dissertation
  signed-i mapping makes this `2.5.1.+0`, the two papers BOTH characterize
  the same trajectory — a cross-validation opportunity.
- [ ] **VISIT cycler ϕ free-parameter**: track that ϕ in VISIT labels is a
  free design variable. To pin it to a V1, we need a separate publication
  (Niehoff refs 6-9, none in our paper corpus yet). Add Niehoff papers to
  #116 acquisitions wishlist.
- [ ] **Track ref 14 = McConaghy-Yam-Landau-Longuski AAS 03-509** as a
  still-outstanding acquisition (the S1L1 primary source per Table 2,
  alongside ref 10 = McConaghy 2004 JSR). Wave digest 2026-06-17 should
  flag this.
