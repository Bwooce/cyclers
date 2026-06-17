# Russell-Ocampo 2003 (AAS-03-145) — full-page deep-read verdict

**Read 2026-06-17 AET.** Part of Mars-cycler digest wave (Agent A, parent task per
`docs/notes/2026-06-17-mars-cycler-wave-digest.md`).

## 1. Header (verbatim from PDF p.1)

- Title: **"A SYSTEMATIC METHOD FOR CONSTRUCTING EARTH-MARS CYCLERS USING DIRECT
  RETURN TRAJECTORIES"** (no asterisk, asterisk is a footnote marker)
- Authors: **Ryan P. Russell† and Cesar A. Ocampo‡**, Department of Aerospace
  Engineering, The University of Texas at Austin, Austin, Texas 78712-1085
  (Russell = Ph.D. candidate `ryanrussell@mail.utexas.edu`; Ocampo = Assistant
  Professor).
- Venue: **AAS Paper 03-145**, AAS/AIAA Spaceflight Mechanics Meeting, Ponce,
  Puerto Rico, February 2003.
- Stamp on top of p.1: "peer-reviewed version: Russell, R. P., Ocampo, C. A.,
  *Systematic Method for Constructing Earth-Mars Cyclers Using Free-Return
  Trajectories*, Journal of Guidance, Control, and Dynamics, Vol. 27, No. 3,
  2004, pp. 321-335."
- DOI of peer-reviewed version (not printed on PDF): **10.2514/1.1909** (JGCD 27(3)).
- Length: 20 pp.

## 2. What the paper actually is

A circular-coplanar **global enumeration** algorithm for Earth-Mars cycler
trajectories indexed by integer tuples `p-h-s-i` (period p in synodic periods;
total half-years h allotted to full/half-rev returns; total symmetric returns s;
i = the multi-revolution Lambert solution index, ordered by ascending semi-major
axis, 1 ≤ i ≤ 2N_max+1). The method **patches** Lambert symmetric returns with
arbitrary sequences of full- and half-revolution direct returns, using **Earth-
generated gravity-assisted maneuvers** at each splice; the algorithm chooses the
splice order and free parameters to **minimize the maximum required geocentric
turn angle** across all flybys (δ_MINIMAX optimization, §"Turn Angle Optimization
for Re-initiation of a Symmetric Return", p.9-12).

The solar system model (p.3) is the simplest possible: **Sun inertially fixed,
Earth circular 1.0 yr period, Mars circular coplanar 1.875 yr period, zero-radius
patched conic, Mars provides no gravity assist, Earth provides instantaneous
zero-SOI gravity-assist** with μ_earth = 3.003489596325074E-6 AU³/TU², 1 AU =
149597870.691 km, 1 TU = 58.1324408670490 days, μ_sun = 1 AU³/TU². The 15-year
Mars-Earth resonance period (Mars 1.875 yr rather than 1.8801 yr) is chosen "to
be consistent with previous studies and for verification purposes" (p.3).

The abstract claim is: **24 purely ballistic cyclers with periods of 2-4 synodic
periods, 92 ballistic cyclers with periods of 5-6 synodic periods, and hundreds
of near-ballistic cyclers** (p.1 abstract). Refined in Conclusions: **total
2502 cyclers found, 116 totally ballistic, 62 close** (p.16, after Figure 13;
the 24+92 = 116 figure matches).

### Title / venue mismatches vs how we cite it

| Where | Citation we use | Notes |
|---|---|---|
| `KNOWN_CORPUS` (Russell-Ocampo / McConaghy SnLm anchor, line 263) | "Russell & Ocampo, J. Spacecraft & Rockets 41(1) 2004" | **DEFECT**: the JGCD paper is Vol. 27, No. 3 (May-June 2004), NOT JSR 41(1). The DOI is 10.2514/1.1909. Easy confusion with the JSR papers (McConaghy 2004 is JSR 41(4)); this is **errata-grade**. |
| Wave digest 2026-06-17 entry | "JGCD 2004 27(3):321-335" | Correct. |
| Filename | "russell-ocampo-2003-systematic-method-earth-mars-cyclers-**direct-return**-trajectories-AAS-03-145.pdf" | Direct-return is the preprint title; the published JGCD title is **"... using Free-Return Trajectories"** (per the red stamp on p.1). Both terms are interchangeable in Russell's usage but the published form is preferred. |
| `catalogue.yaml` source tag `russell-2004-t34` | Maps to Russell **dissertation** Table 3.4 | Cross-check: catalogue rows like `russell-ocampo-2.1.1+2-case2` use signed `i` notation (`+2`, `-5`), while this paper's Table 4/A1/A2 print unsigned `i` (e.g. "2-1-1-5"). The signed-`i` form is the dissertation convention. Catalogue is sourced from the dissertation; this paper is the JGCD progenitor with the same tables shown unsigned. |

## 3. Key numerical / structural content

### Definitions and nomenclature (§"Cycler Definition" p.3; intro p.2)
- Cycler = perfectly repeatable round-trip; period must be integer multiple of
  synodic period; "Once set in orbit, ballistic cyclers require no powered
  maneuvers to maintain" (p.3).
- `Cycler-p-h-s-i`: **p** = period in synodic periods; **h** = total half-years
  allotted to all full/half-rev returns (h ≥ 0); **s** = number of identical
  symmetric returns per cycler period (s ≥ 1); **i** = index of the multi-rev
  Lambert solution by ascending semi-major axis (1 ≤ i ≤ 2N_max+1) (p.1
  nomenclature; p.2 introductory paragraph).
- "Aphelion Ratio AR = max aphelion radius / 1.52 AU (Mars radius)"; AR > 1 ⇒
  cycler intercepts Mars without powered maneuver (p.13).
- "Turn Ratio TR = max physically allowable turn angle / max required turn angle
  (δ_MAX)"; TR > 1 ⇒ all flybys physically attainable; max allowable based on
  **200 km altitude Earth flyby** (p.13). TR ≥ 1 is "strictly ballistic".

### Solver core (Eqs. 1-10)
- Eq. (1) Energy equation for full-rev return.
- Eq. (2) Lagrangian Lambert formulation (multi-rev), with α/β intermediate
  variables (p.6).
- Eq. (3) `TOF = T_SYN·p − h/2` cycler period closure constraint (p.6).
- Eqs. (4-10) Turn-angle minimization (δ_MIN, δ_a, δ_b, δ_c) for symmetric-return
  re-initiation (p.11-12).
- Table 3 (p.11): Coordinates for heliocentric flyby velocities (a 2-row schematic,
  NOT a cycler table — disambiguation for the `russell-2004-t34` provenance tag).

### Cycler tables (the catalogue-relevant ones)

**Table 4** (p.14) — Two-, three-, and four-synodic-period ballistic or
near-ballistic cyclers, AR_MIN=0.9, TR_MIN=0.9. **39 rows** (24 with TR>1
ballistic + 15 near-ballistic).

Columns: Cycler-p-h-s-i / Aphelion Ratio / Turn Ratio / Earth→Mars (or aphelion)
Time (days) / Earth V∞ (km/s) / Mars V∞ (km/s) / Required Geocentric Turning
Angle at each Flyby (deg, per-flyby).

Notable rows (full row text):
- `2-1-1-5ᵃ` AR=0.95 TR=1.11 207 d 4.1 km/s @E 2.0 km/s @M turns 92,92.
  Footnote a: **"Case 2" cycler described by Byrnes³**.
- `2-3-1-5ᵇ` AR=1.08 TR=0.92 143 d 5.4 km/s @E 5.3 km/s @M turns 93,93.
  Footnote b: **"Case 3" cycler described by Byrnes³**.
- `2-5-1-3` AR=1.44 TR=1.12 **94 d** 7.8 km/s @E 9.9 km/s @M turns 54,54,54,54.
- `3-7-1-9` AR=1.07 TR=1.56 175 d 3.6 km/s @E 4.6 km/s @M turns 71,71,71,71
  ("Cycler-3-7-1-9 has low terminal speeds and requires 4 flybys", p.14
  continuing text).
- `4-3-1-20` AR=0.99 TR=1.29 268 d 3.1 km/s @E 2.5 km/s @M turns 93,93
  (p.15: "Cycler-4-3-1-20 has remarkably low energy requirements at Earth and
  Mars … V∞ of 3.10 km/s at Earth compared to the Hohmann value of 2.84 km/s,
  while at Mars the cycler has a V∞ of 2.53 km/s compared to the Hohmann value
  of 2.57 km/s. The Aphelion Ratio is 0.992, thus the cycler doesn't quite
  reach Mars in the simplified model.").
- `4-3-1-22` AR=1.26 TR=1.01 **154 d** 4.7 km/s @E 7.6 km/s @M turns 94,94.

**Tables 5-8** (p.15) — **Per-encounter Δv vectors with full state** for cyclers
2-5-1-3, 3-1-2-11, 4-3-1-20, 4-5-2-12. Each table prints:
- Location (Earth / Mars / Earth / Earth ...)
- time (days, full sequence)
- Δv_x, Δv_y, Δv_z in km/s for each encounter
- Initial Mars position vector `r_mars` at t_0 in AU.

Example, **Table 5 Cycler-2-5-1-3** (p.15):
> Location: Earth, Mars, Earth, Earth, Earth, Earth, Mars
> time (days): 0, 94, 652, 1018, 1200, 1565, 1659
> Δv_x: 6.50ᵃ, 0, −5.19, 1.40, −1.40, −5.29, 0
> Δv_y: 4.35ᵃ, 0, −1.41, −6.12, 6.12, −0.98, 0
> Δv_z: 0ᵃ, 0, 4.55, 3.20, 3.20, 4.55, 0
> r_mars at t_0 = [1.41, 0.57, 0] AU
> Footnote a: "powered Δv required to initiate cycler from Earth".

**These tables print enough state to simulate one full cycle plus the first leg
of the second cycle** (p.15 lead-in). This is **V1-grade reproducible data in
the published, circular-coplanar model** — but only for the 4 cyclers in Tables
5-8 (`2-5-1-3`, `3-1-2-11`, `4-3-1-20`, `4-5-2-12`).

**Figures 11 (Sun-fixed + Earth/Mars-pulsating views of 4-3-1-20 and 2-5-1-3)
and 12 (3D view of 3-1-2-11)** — qualitative geometry, not numbers.

**Figure 13** (p.16) — Bar chart of cycler counts per period p=1..6 (Total / TR
and AR > 0.9 / TR and AR > 1). Per text: 2502 total, 116 ballistic, 62 close.

**Appendix Tables A1 (5-synodic, ~80 rows) and A2 (6-synodic, ~150 rows)** (p.18,
19) — same column structure as Table 4 (Aphelion Ratio / Turn Ratio / Earth→Mars
days / V∞ at Earth / V∞ at Mars / per-flyby turn angles).

Notable A1/A2 rows:
- A1 `5-2-1-9` AR=0.90 TR=1.07 182 d V∞_E=4.5 V∞_M=1.3, turns 92,92.
- A1 `5-2-2-11` AR=1.20 TR=0.94 128 d V∞_E=5.2 V∞_M=7.1, turns 94,94,85.
- A2 `6-0-1-23ᵃ` AR=0.92 TR=1.40 213 d V∞_E=5.0 V∞_M=1.2, turn 86.
  Footnote a: **"Cycler 6S9" described by McConaghy² et al**.
- A2 `6-0-1-25ᵇ` AR=1.03 TR=1.22 179 d V∞_E=4.0 V∞_M=3.6, turn 85.
  Footnote b: "Cycler 6S8" described by McConaghy.
- A2 `6-0-1-27ᶜ` AR=1.17 TR=1.07 133 d V∞_E=5.0 V∞_M=6.7, turn 83.
  Footnote c: "Cycler 6S7" described by McConaghy.
- A2 `6-0-1-29ᵈ` AR=1.34 TR=0.93 111 d V∞_E=6.0 V∞_M=8.7, turn 84.
  Footnote d: "Cycler 6S6" described by McConaghy.

### Aldrin cycler identification
p.16 (last paragraph before Conclusions): **"Cycler-1-0-1-6, with a Turn Ratio
of 0.86, is the Aldrin cycler⁷"** (ref 7 = Byrnes-Longuski-Aldrin JSR 1993).
The Aldrin cycler is **p=1, h=0, s=1, i=6** in Russell's enumeration; TR<1
confirms its powered character (matching catalogue's `aldrin-cycler` row note
about required ~84 deg turn vs achievable ~72 deg).

### References (p.20)
The paper references **McConaghy et al. AIAA-2002-4420 (ref 2), Byrnes
AIAA-2002-4423 (ref 3), Byrnes-Longuski-Aldrin JSR 30(3) 1993 (ref 7), and
Aldrin et al. NASA/JPL Contract 1230398 Dec 15 2001 (ref 1)** — i.e. it cites
the 1993 JSR Byrnes-Longuski-Aldrin paper directly. No DOI for ref 7.

## 4. `KNOWN_CORPUS` impact (RECOMMEND, do not edit)

The current Russell-Ocampo / McConaghy SnLm anchor (literature_check.py line
263, `citation` field, line 276-277):

```
"Russell & Ocampo, J. Spacecraft & Rockets 41(1) 2004; "
"McConaghy et al., J. Spacecraft & Rockets 43(2) 2006",
```

**DEFECTS**:
1. **Russell-Ocampo is JGCD 27(3) 2004, NOT JSR 41(1) 2004.** The JSR 41(1) 2004
   entry is for a different Russell-Ocampo paper or a paper conflation. The
   correct citation is:
   > Russell, R. P. & Ocampo, C. A., "Systematic Method for Constructing
   > Earth-Mars Cyclers Using Free-Return Trajectories," *J. Guidance, Control,
   > and Dynamics*, Vol. 27, No. 3, May-June 2004, pp. 321-335. DOI 10.2514/1.1909.
2. **Title shorthand**: Anchor uses "Russell-Ocampo / McConaghy Earth-Mars SnLm
   cyclers"; OK, but worth noting the published-version title uses "Free-Return
   Trajectories" (not "Direct Return" — the preprint title).
3. **AAS-03-145 preprint as a secondary citation** is worth adding for archival
   completeness (the preprint has the same tables but unsigned `i`; the
   catalogue uses signed `i` from the **dissertation** Table 3.4, which is a
   third source).

Recommended fix: tighten the citation field to:

```
"Russell & Ocampo, JGCD 27(3):321-335 (2004), DOI 10.2514/1.1909, preprint "
"AAS-03-145; Russell, R. P. PhD diss. (UT Austin, 2004) Table 3.4 for the "
"signed-i form used in catalogue IDs; McConaghy et al., JSR 43(2) 2006 "
"DOI 10.2514/1.15215",
```

## 5. Catalogue impact (RECOMMEND)

### Rows touched

Per the line-grep above, **~190 `russell-ocampo-*` rows** in `data/catalogue.yaml`
use `orbit_source: russell-2004-t34` / `vinf_source: russell-2004-t34`. The
`t34` tag = dissertation Table 3.4. **All are currently V0 per the validation
ceiling memory** (no per-member orbital elements were ever printed in any
source).

### Does Russell-Ocampo 2003 enable V0 → V1 promotions?

**HONEST NEGATIVE for the majority. PARTIAL POSITIVE for 4 specific rows.**

- The **catalogue's `russell-ocampo-*` IDs** use signed `i` (e.g. `2.1.1+2`,
  `4.3.1-5`). Russell-Ocampo 2003 Tables 4/A1/A2 use **unsigned `i`** as a
  positional index. The catalogue's signed convention comes from Russell's
  dissertation (Table 3.4), not from this paper. Cross-referencing requires
  knowing the dissertation re-numbering.
- **For Tables 4 / A1 / A2 cycler rows (the vast majority of catalogue
  `russell-ocampo-*` IDs)**, this paper prints only: Aphelion Ratio, Turn
  Ratio, Earth→Mars days, V∞ Earth, V∞ Mars, per-flyby turn angles. **It does
  NOT print per-arc orbital elements (a, e, ω, ν) for the symmetric return or
  generic returns.** Therefore — under the §14 V1 like-for-like rule (which
  requires "per-member orbital elements that we can reproduce in same-model
  corrector") — these rows **remain V0 by publication gap, exactly as the
  validation-ceiling memory documents**.
- **HOWEVER — for the 4 cyclers in Tables 5-8** (`2-5-1-3`, `3-1-2-11`,
  `4-3-1-20`, `4-5-2-12`), the paper prints **per-encounter Δv vectors with
  time + initial r_mars** — enough to simulate the full cycler trajectory in
  the same circular-coplanar model. **These four are V1 candidates.**
  - `russell-ocampo-2.5.1+0` (line 1366 of catalogue) is the catalogue version
    of `Cycler-2-5-1-3` — should map. **Promote to V1 candidate.**
  - `russell-ocampo-3.1.2+1` (line 1201) **may** be `Cycler-3-1-2-11`
    (the catalogue's `+1` suggests dissertation re-numbering; needs
    cross-check against the dissertation Table 3.4 to confirm Lambert index
    mapping). **V1 candidacy pending dissertation cross-ref.**
  - **No catalogue row for `Cycler-4-3-1-20` or `Cycler-4-5-2-12`** at the
    appropriate IDs — `4.3.1+...` rows in catalogue include `4.3.1-5` and
    `4.3.1-4`, but the dissertation may renumber `i=20` to a signed offset.
    Worth adding rows if these aren't already there under different IDs.

### `mcconaghy-2006-em-k2` cross-reference

This paper's Table A2 footnote a (cycler `6-0-1-23`) is "Cycler 6S9" of
McConaghy 2006 — i.e. the **9th six-synodic cycler from McConaghy 2006**, not
the McConaghy 2006 *notable* S2 cycler. The McConaghy 2006 'notable' is a
two-synodic-period cycler so it is NOT in Russell-Ocampo's Table A2 (6-synodic);
it would appear in Table 4 (2-synodic) — but no Table 4 row is footnoted as a
McConaghy 2006 row. So Russell-Ocampo 2003 **does not corroborate the
mcconaghy-2006-em-k2 row's specific numerical values** (V∞_E=4.7, V∞_M=5.0,
ToF=153 d). That cross-check still depends on McConaghy 2006 itself.

### Aldrin cycler row

The catalogue's `aldrin-cycler` row (line ~30) is the **only V3 cross-validated
Aldrin row**. Russell-Ocampo 2003 p.16 confirms it = `Cycler-1-0-1-6` with
TR=0.86 — consistent with the catalogue's "powered in idealized model" note.
**No catalogue change needed** for the Aldrin row from this paper.

## 6. Errata / surprises

1. **Citation field defect in `KNOWN_CORPUS`** — Russell-Ocampo is JGCD 27(3),
   not JSR 41(1). Easy to confuse with McConaghy 2004 (JSR 41(4)). Errata-grade.
2. **Preprint title vs published title** — preprint says "Direct Return", JGCD
   published version says "Free-Return". Our filename uses the preprint title.
   Both refer to the same thing in Russell's usage. Worth noting in the corpus
   citation field that both titles refer to the same paper.
3. **Dissertation re-numbering** — Tables 5-8 of this paper give per-encounter
   state for 4 specific cyclers (`2-5-1-3`, `3-1-2-11`, `4-3-1-20`, `4-5-2-12`)
   that map (with signed-`i` re-index) to catalogue rows. **These 4 rows are
   V0→V1 promotion candidates** under the §14 like-for-like rule — but ONLY
   in the circular-coplanar simplified-model.
4. **AR_MIN=0.9 cutoff is arbitrary** (p.14: "somewhat arbitrary"). Anything
   we mine *below* AR=0.9 from a continuation/re-search effort is consistent
   with Russell's published results but not specifically endorsed by them.
5. **Mars period 1.875 yr (not 1.8801 yr)** — used "for verification purposes"
   (p.3). Any same-model reproduction MUST use 1.875 yr to match.
6. **Earth μ value** — μ_earth = 3.003489596325074E-6 AU³/TU² with explicit
   15 significant figures. Any same-model reproduction MUST use this exact
   value to match Tables 5-8 state vectors.

## 7. Action items for parent

- [ ] **Tighten `KNOWN_CORPUS` Russell-Ocampo citation** — fix the JSR 41(1)
  defect to JGCD 27(3), add DOI 10.2514/1.1909, add the preprint AAS-03-145 as
  a secondary citation, add the dissertation Table 3.4 as the source-of-record
  for the catalogue's signed-`i` IDs.
- [ ] **V1 candidate audit for `russell-ocampo-2.5.1+0`** — does the catalogue
  row's Δv values match Russell-Ocampo 2003 Table 5? Cross-check delta_v_kms
  and trajectory.segments. If yes → V1 with per-encounter-state-reproduces
  evidence. If not → flag a numeric inconsistency.
- [ ] **Cross-reference catalogue `russell-ocampo-3.1.2+1` against Table 6
  (`Cycler-3-1-2-11`)** — confirm the dissertation's signed-`i` mapping for
  i=11 → +1 (or whatever offset). Likely V1 candidate if mapping holds.
- [ ] **Hunt for catalogue rows matching `Cycler-4-3-1-20` and `Cycler-4-5-2-12`**
  in Tables 7-8. If absent, consider adding them — they're the highest-leverage
  rows in the paper (Table 7's `4-3-1-20` has the famous "near-Hohmann" 3.10
  km/s @E + 2.53 km/s @M property quoted on p.15).
- [ ] **Same-model reproduction harness** — to actually execute V0→V1 promotion
  for the 4 Tables 5-8 cyclers, we need a Lambert + zero-SOI-flyby integrator
  with **exactly** μ_earth = 3.003489596325074E-6 AU³/TU², 1 AU =
  149597870.691 km, 1 TU = 58.1324408670490 days, Mars period 1.875 yr. Our
  current `correct.py` / `lambert.py` may already meet this — needs a one-off
  golden test. Track separately.
- [ ] **Filename hygiene** — the preprint filename uses "direct-return", the
  JGCD publishes as "free-return". Either keep the preprint name as-is (current
  filename is the AAS-03-145 stamp) or annotate `papers/README.md` to note the
  rename in the JGCD version. No urgent action.
- [ ] **Russell-Ocampo 2003 RESOLVES the #116 Russell-Ocampo wishlist item.**
  Mark wishlist-resolved.
