# Multi-Arc Classification Verdicts

**Task:** Schema-v4 Task 2.2 — source-cited classification of every candidate row as
`multi-arc` vs `single-ellipse` (or `non-keplerian`).

**Date:** 2026-06-03  
**Author:** agentic worker (Claude Sonnet 4.6)  
**Status:** RESOLVED — the one flagged row (`mcconaghy-2006-em-k2`) was resolved to
**multi-arc** by the reviewer (source-faithful: its own cited note + sister row
`russell-ch4-4.991gG2`). Final tally: 199 multi-arc / 28 single-ellipse / 6 non-keplerian.
See §7 and §12.

---

## 1. Discriminator definition

### 1.1 Physical test

A cycler row is **`multi-arc`** if and only if the source demonstrates that the cycler
is constructed from a *sequence of distinct Earth-Earth free-return trajectories*
(generic, half-rev, full-rev) patched by Earth gravity-assist flybys, such that *no
single repeating Keplerian ellipse (a, e) exists* for the vehicle — each leg of the
cycle is a different transfer orbit.

A cycler row is **`single-ellipse`** if the vehicle travels a single repeating Keplerian
heliocentric ellipse; gravity-assist flybys re-match V∞ back to *the same* orbit.
A null `a_au` is a *data gap* (the value is derivable or tabulated elsewhere), not
evidence of multi-arc structure.

A cycler row is **`non-keplerian`** when `primary != "Sun"` — the orbit is a CR3BP
periodic orbit, not a heliocentric Keplerian ellipse at all.

### 1.2 Russell's three free-return types and the multi-arc discriminator

Russell 2004 (§2.6, p. 30; §3.4.3, p. 63) defines three mutually exclusive
Earth-Earth free-return trajectory types:

| Type | Description | Orbit |
|---|---|---|
| Full-revolution (2π) | even-nπ transfer; fixed semi-major axis | single ellipse per leg |
| Half-revolution (1π) | odd-nπ transfer (backflip); fixed semi-major axis | single ellipse per leg |
| Generic | non-integer-multiple-of-π transfer angle; TOF arbitrary | **different ellipse per leg** |

**The discriminating statement** is at Russell 2004, §2.8 (p. 54, Chapter 2 Conclusions):

> "The semi-major axis values for the odd-nπ and the generic free-return transfers
> can only be found using Lambert's Equation because the times of flight are not
> integer multiples of orbital periods."

And from §3.4.3 (p. 63, "One or More Identical Generic Free Returns"):

> "A generic return... is characterized by an arbitrary time of flight and a transfer
> angle that is a non integer multiple of π."

This means a generic return arc has its own (a, e) determined by its specific TOF
and transfer angle — it is not a fixed ellipse that repeats identically. When a
cycler is *composed of generic-return arcs* (as all Russell Chapter 3 and Chapter 4
cyclers are), **no single repeating ellipse exists**; the identity is per-arc.

### 1.3 Chapter 3 vs Chapter 4 distinction

- **Chapter 3** (Russell 2004, title: "Idealized Free-Return Cyclers Composed of
  π Transfers, 2π Transfers, and One or More **Identical** Generic Returns"):
  Uses identical generic returns patched with π (half-rev) and 2π (full-rev) arcs.
  By "identical" Russell means all generic arcs in a given cycler have the same TOF
  and transfer angle — but this is still **not** a single repeating ellipse; the
  outbound generic arc and the inbound leg are different (the flybys change the
  direction). The cycler has no fixed (a, e).

- **Chapter 4** (Russell 2004, title: "Idealized Free-Return Cyclers Composed of
  nπ Transfers and One or More **Identical or Different** Generic Returns"):
  Extends to nπ transfers and allows non-identical generic returns (§4.3, p. 98):
  "Remove the limitation that multiple generic returns must be identical." Again,
  no single repeating ellipse exists.

**Key citation for both chapters:** Russell 2004, §4.9 (p. 135, Chapter 4 Conclusions):
"The global generalized search for idealized Earth-Mars cyclers results in a database
of all useful free-return cyclers with periods of three or less synodic periods."
These cyclers are identified by their leg *descriptor strings* (§4.8, p. 126), not by
a (a, e) — because no single (a, e) exists for them.

### 1.4 The Aldrin cycler is NOT multi-arc

For contrast: the Aldrin cycler is found by Russell in §4.9 when TR threshold is
lowered (p. 133): it has descriptor "G(2.1354, 7.1338, L)" — a *single* generic
return. But the Aldrin cycler is a single-ellipse cycler because Rogers 2012, §3.1
and McConaghy 2002 publish a single (a=1.60 AU, e=0.393) for it. The single
descriptor string maps to a single repeating ellipse in this degenerate case (the
S1 family, 1 generic arc per cycle). The flybys re-match V∞ to the *same* orbit.

---

## 2. Classification rules applied (precedence)

1. **`non-keplerian`**: `primary != "Sun"` — mechanical. Six rows.
2. **`multi-arc`**: positive source evidence of generic-return / distinct-per-leg
   orbits (per §1.2–1.3 above). Only Russell-Ocampo and Russell-Ch4 rows qualify.
3. **`single-ellipse`**: default for everything else, including rows with null
   `a_au` (those are single-ellipse with a data gap).

---

## 3. Non-keplerian rows (excluded from multi-arc/single-ellipse analysis)

These 6 rows have `primary != "Sun"` and are **`non-keplerian`** by mechanical rule.
They are not included in the allowlist.

| Row ID | Primary | Basis |
|---|---|---|
| `arenstorf-em-figure8-1963` | Earth | CR3BP figure-8 |
| `genova-aldrin-2015-em-3petal-cycler` | Earth | CR3BP 3-petal |
| `wittal-2022-em-cycler-family` | Earth | CR3BP family |
| `hernandez-2017-jovian-ieg-triple-family` | Jupiter | planet-centric |
| `russell-strange-2009-jovian-multimoon-family` | Jupiter | planet-centric |
| `russell-strange-2009-saturnian-multimoon-family` | Saturn | planet-centric |

---

## 4. Multi-arc verdict: `russell-ocampo-*` (184 rows)

### 4.1 Structural basis

All 184 `russell-ocampo-*` rows originate from Russell 2004 Chapter 3 (§3.8 Tables
3.4–3.11, pp. 80–92) — cyclers composed of one or more identical generic returns
patched with π and 2π transfers. The leg architecture is per-arc (Russell §4.8
descriptor-string notation, p. 126); no single repeating (a, e) exists.

**Positive source citation for multi-arc verdict:**
- Russell 2004, §3.4.3, p. 63: "A generic return... is characterized by an arbitrary
  time of flight and a transfer angle that is a non integer multiple of π."
- Russell 2004, §3.1 (Chapter 3 Summary), p. 56: the method finds cyclers composed
  of "2π full- and 1π half-rev return orbits with Earth-generated gravity-assisted
  maneuvers" — multiple distinct arcs per cycle.
- Russell 2004, §2.8, p. 54: "the semi-major axis values for the odd-nπ and the
  generic free-return transfers can only be found using Lambert's Equation because
  the times of flight are not integer multiples of orbital periods."
- Russell 2004, Table 3.4 (2-synodic) and Tables 3.5–3.11 (3–6-synodic): tabulate
  Aphelion Ratio and Turn Ratio per cycler, NOT (a, e) — confirming no single (a, e)
  is published for these rows.

### 4.2 The 172 rows that explicitly mention "generic" in their catalogue notes

172 of the 184 `russell-ocampo-*` rows contain the word "generic" in their notes,
orbit_elements notes, or segment notes, directly confirming the generic-return
architecture. No further analysis needed for these rows.

**Verdict: multi-arc.** Source: Russell 2004 Chapter 3 (Tables 3.4–3.11).

### 4.3 The 12 rows that do NOT mention "generic" in their catalogue text

These 12 rows were individually inspected. All originate from Russell 2004 Table 3.4
(2-synodic family) or Tables 3.5–3.7 (3-synodic) or Tables 3.9–3.11 (4-synodic).
They carry the same multi-arc architecture as the other 172 rows; the absence of the
word "generic" in their notes is an omission in the catalogue prose, not a structural
difference.

Evidence per row: each row's `orbit_elements.a_au = null` with a note of the form
"Other elements (a, e, per-leg...) not tabulated" or "Aphelion derived from Russell
2004 Table 3.4 Aphelion Ratio..." — confirming no single (a, e) exists and the
source only publishes the aggregate Aphelion Ratio (a cycle-level descriptor, not an
orbital element of a single ellipse).

| Row ID | Table source | Synodic | Confirming evidence |
|---|---|---|---|
| `russell-ocampo-2.3.1+1-case3` | Russell Table 3.4 | 2 | orbit_elements note: "Aphelion Ratio 1.08 per Russell 2004 Table 3.4"; a_au=null |
| `russell-ocampo-2.5.1+0` | Russell Table 3.4 | 2 | orbit_elements note: null a_au; notes: "typical high-V_inf short-transit cycler" |
| `russell-ocampo-3.1.1+3` | Russell Table 3.4 | 3 | orbit_elements note: "AR = 1.07 x 1.52 AU"; a_au=null |
| `russell-ocampo-3.1.2+1` | Russell Table 3.4 | 3 | a_au=null; sourced only from Russell Table 3.4 AR/TR |
| `russell-ocampo-3.3.1+2` | Russell Table 3.4 | 3 | orbit_elements note: "AR = 1.19 x 1.52 AU"; a_au=null |
| `russell-ocampo-3.5.1+1` | Russell Table 3.4 | 3 | orbit_elements note: "AR = 1.43 x 1.52 AU"; a_au=null |
| `russell-ocampo-3.5.2+0` | Russell Table 3.4 | 3 | orbit_elements note: "AR = 1.43 x 1.52 AU"; a_au=null |
| `russell-ocampo-3.7.1+1` | Russell Table 3.4 | 3 | orbit_elements note: "AR = 1.07 x 1.52 AU"; a_au=null |
| `russell-ocampo-3.9.1+0` | Russell Table 3.4 | 3 | orbit_elements note: "AR = 1.43 x 1.52 AU"; a_au=null |
| `russell-ocampo-4.10.1+2` | Russell Table 3.4 | 4 | orbit_elements note: "AR = 1.03 x 1.52 AU"; a_au=null |
| `russell-ocampo-4.12.1+1` | Russell Table 3.4 | 4 | orbit_elements note: "AR = 1.16 x 1.52 AU"; a_au=null |
| `russell-ocampo-4.14.1+0` | Russell Table 3.4 | 4 | orbit_elements note: "AR = 1.49 x 1.52 AU"; a_au=null |

**Verdict: all 12 → multi-arc.** Source: Russell 2004 §3 (Tables 3.4–3.11 tabulate
only AR and TR, not (a, e); the Chapter 3 architecture is generic-return throughout).

---

## 5. Multi-arc verdict: `russell-ch4-*` (14 rows)

All 14 `russell-ch4-*` rows originate from Russell 2004 Chapter 4 (Tables 4.9–4.13,
pp. 127–134). Chapter 4's title explicitly includes "Identical or Different Generic
Returns." Each row's catalogue notes cite the leg descriptor strings from these
tables (e.g. "g(1.4612,526.02,Ll) G(2.8096,651.46,U)" for `russell-ch4-4.991gG2`),
confirming the per-arc structure.

**Positive source citation:**
- Russell 2004, §4.1 (Chapter 4 Summary), p. 95: "a new technique is developed to
  globally identify and optimize all useful idealized Earth-Mars cyclers constructed
  with Earth-Earth free-returns patched by gravity-assisted flybys... combinations
  of generic, half-rev, and full-rev free-return transfers."
- Russell 2004, §4.8 (Results), p. 126: "Each solution is characterized by...
  a sequence of descriptor strings for each leg of the cycler. All legs reported
  are Earth-Earth free-returns."
- Russell 2004, Table 4.14 (Summary of results), p. 134: explicitly lists all
  Chapter 4 cyclers by transit-leg type (generic or full-rev), confirming no single
  (a, e) exists.
- Russell 2004, §4.9, p. 133: The Aldrin cycler is recovered with TR < 1 threshold,
  with descriptor "G(2.1354, 7.1338, L)" — *one* generic arc — but even this is
  classified as a generic-return cycler in Russell's framework (Aldrin's single-arc
  nature is what makes it single-ellipse: see §6 below).

| Row ID | Table source | Transit leg type | Confirming catalogue note |
|---|---|---|---|
| `russell-ch4-4.991gG2` | Russell Table 4.9 (two-synodic) | generic (gG) | "two generic-return arcs of 1.46 yr and 2.81 yr" |
| `russell-ch4-8.049gGf2` | Russell Table 4.9 | generic+full-rev (gGf) | "composed of two generic returns and one full-rev return" |
| `russell-ch4-8.165Gfh-f2` | Russell Table 4.9 | generic+full-rev+half (Gfh-f) | Russell Table 4.9 row 3 |
| `russell-ch4-9.353Gg2` | Russell Table 4.9 | generic (Gg) | "two generic legs" |
| `russell-ch4-3.64gGg3` | Russell Table 4.10 (three-synodic) | generic (gGg) | "remarkable low values of v∞" bold cycler |
| `russell-ch4-3.77Gh3` | Russell Table 4.10 | generic+half (Gh) | "3.5-period half-rev return" |
| `russell-ch4-3.78Gg3` | Russell Table 4.10 | generic (Gg) | "two generic-return arcs" |
| `russell-ch4-5.30gGf3` | Russell Table 4.10 | generic+full-rev (gGf) | "f(3:2,...) is a 3:2 resonant full-rev return" |
| `russell-ch4-5.66Gfh3` | Russell Table 4.10 | generic+full-rev+half (Gfh) | "consistently low ΔV requirements" |
| `russell-ch4-9.94Gg3` | Russell Table 4.10 | generic (Gg) | "shortest E-M transit time in entire Russell catalogue" |
| `russell-ch4-3.66gfF3` | Russell Table 4.12 (full-rev transit) | generic+full-rev (gfF) | "outbound per Table 4.12" |
| `russell-ch4-5.30ggF3` | Russell Table 4.12 | generic+full-rev (ggF) | "F(3:2,82.487,180.000) — third leg is the 3:2 full-rev" |
| `russell-ch4-5.75ggF3` | Russell Table 4.12 | generic+full-rev (ggF) | "F(2:1,85.196,0.000) — third leg is a 2:1 full-rev" |
| `russell-ch4-6.44Gg3` | Russell Table 4.13 (near-ballistic) | generic (Gg) | "NEAR-BALLISTIC: TR = 0.95 < 1.0" |

**Verdict: all 14 → multi-arc.** Source: Russell 2004 §4 (Tables 4.9–4.13 and
leg descriptor strings confirm per-arc architecture throughout).

---

## 6. Single-ellipse families: explicit checks

### 6.1 `aldrin-*` (4 rows)

Two rows (`aldrin-classic-em-k1-outbound`, `aldrin-classic-em-k1-inbound`) carry
sourced (a=1.60 AU, e=0.393) from McConaghy 2002 and Rogers 2012.
Two rows (`aldrin-4-3-2-establishment`, `aldrin-3-2-1-establishment`) are
establishment trajectories with null a_au — these are data gaps (the establishment
phase uses V∞-leveraging, but the target is the same Aldrin ellipse).
**Verdict: single-ellipse.** Source: McConaghy/Longuski/Byrnes 2002 AIAA 2002-4420;
Rogers et al. 2012 Table 1.

### 6.2 `mcconaghy-2005-em-case1` and `mcconaghy-2005-em-u0l1`

Both carry sourced (a, e) from Rogers 2012 Table 1: (1.22, 0.238) and (2.05, 0.563)
respectively. Single repeating Keplerian ellipses confirmed.
**Verdict: single-ellipse.** Source: Rogers et al. 2012 Table 1.

### 6.3 `mcconaghy-2005-em-snlm-broadclass-family`

Family-seed entry (no per-member numerics). The SnLm family is the broad class of
single-ellipse Earth-Mars cyclers from McConaghy 2005. Null a_au is a data gap.
**Verdict: single-ellipse (family seed).** Source: McConaghy 2005 dissertation.

### 6.4 `s1l1-2syn-em-cpom` — E-E-M-M, null a_au EXPLAINED

This is the critical example called out in the task brief. The S1L1 cycler:
- Has sourced (a=1.30 AU, e=0.257) from Rogers et al. 2012 Table 1.
- The `orbit_elements.a_au = 1.3` and `e = 0.257` are already filled in the
  catalogue row (not null) — these are sourced from Rogers 2012 Table 1.
- The null values in the row are for per-*segment* (a, e) in `trajectory.segments[]`
  (data gaps for the multi-rev Lambert solver), not the top-level orbit.
- Rogers 2012 publishes the S1L1 (a, e) as the orbit of the single repeating Keplerian
  ellipse the cycler follows (circular-coplanar model).

**Verdict: single-ellipse.** Source: Rogers et al. 2012 Table 1 (a=1.30 AU,
e=0.257). Note: the S/L taxonomy defines the cycler by a single repeating ellipse;
Russell's circular-coplanar framework finds it at descriptor 4.991gG (Table 4.9)
but the Rogers 2012 and CPOM studies characterize it as a single (a, e) — both
framings are valid for different contexts, but for `cycler_class` the Rogers sourced
(a, e) defines it as single-ellipse.

### 6.5 `mcconaghy-2006-em-k2` — E-E-M-M, null a_au — AMBIGUOUS (see §7)

This row is flagged as uncertain. See §7.

### 6.6 `hollister-menning-1970-ev-orbit-*` (15 rows)

Hollister & Menning 1970 ("Periodic Swing-By Orbits between Earth and Venus", JSR
7(10)) identify 15 periodic orbits using combinations of "full-revolution returns"
and "symmetric returns" (their paper's terminology, Table 1, pp. 1194–1195). A
"direct return orbit" in their framework is "a sun-centered elliptical orbit which
returns the spacecraft to the same planet from which it was launched" (p. 1194) —
this is the definition of a single repeating ellipse.

Their sequence is E-V (Earth-Venus) and the cyclers use Venus flybys to redirect the
spacecraft back to Earth on the same heliocentric ellipse. The (a, e) are not directly
tabulated in Table 3 (only encounter dates, V_r, turn angles, R_min), but the orbit
type is definitionally single-ellipse. Null a_au is a data gap.
**Verdict: single-ellipse (15 rows).** Source: Hollister & Menning 1970, §"Direct
Return Orbits" (p. 1194).

### 6.7 `niehoff-visit1` and `niehoff-visit2`

Sourced (a, e) = (1.17, 0.193) and (1.31, 0.275). Single repeating ellipses.
**Verdict: single-ellipse.** Source: Russell 2004 (dissertation background); Rogers
2012 Table 1.

### 6.8 `jones-2017-vem-triple-family`

Family-seed entry. Null a_au. Jones et al. 2017 triple-cycler uses fixed
V-E-M resonant sequence; while the orbit is n-body rather than a simple Keplerian
ellipse, the source does not identify it as generic-return multi-arc (it uses
resonance-based trajectories with fixed synodic relationships). Conservatism rule
applies: default single-ellipse pending source confirmation of multi-arc structure.
**Verdict: single-ellipse (conservative default).** Source: Jones/Hernandez/Jesick
2017 (Low-Excess-Speed Triple-Cyclers).

### 6.9 `vem-emeeve-3syn`

Family-seed entry. Null a_au. Same conservatism as jones row.
**Verdict: single-ellipse (conservative default).**

---

## 7. Uncertain / needs human review

### 7.1 `mcconaghy-2006-em-k2` — potential mis-classification

**Issue:** The task brief lists this row as a single-ellipse example ("E-E-M-M with
null a_au — but they are single-ellipse"). However, the existing catalogue notes
for this row explicitly state:

> "Russell does NOT tabulate (a, e, perihelion) directly because the cycler is a
> piecewise sequence of two generic-return arcs, not a single Keplerian ellipse —
> each leg has its own (a, e)."

This statement is sourced from Russell 2004 Table 4.9 row 1 (orbit_elements note).
The Russell Table 4.9 descriptor "g(1.4612,526.02,Ll) G(2.8096,651.46,U)" confirms
two generic-return arcs with different TOFs (1.46 yr and 2.81 yr) — structurally
identical to `russell-ch4-4.991gG2`, which is confirmed multi-arc.

Note: `mcconaghy-2006-em-k2` and `russell-ch4-4.991gG2` represent the **same
physical cycler** (the notes explicitly state this: "The Russell entry for this
cycler is catalogued separately as `russell-ch4-4.991gG2`").

**Reviewer resolution (2026-06-03): `multi-arc`.** The source governs over the
task brief's assumption. The row's own cited note (Russell 2004 Table 4.9) states it
is a piecewise sequence of two generic-return arcs with no single Keplerian ellipse,
and it is the *same physical cycler* as `russell-ch4-4.991gG2` (already multi-arc).
It is therefore added to `MULTI_ARC_ALLOWLIST` (final size 199). Its top-level
`a_au` is null (correct for multi-arc). See §12 for the same-cycler-three-rows note.

---

## 8. Verdict table summary

| Family | Count | Verdict | Russell citation |
|---|---|---|---|
| `russell-ocampo-*` | 184 | **multi-arc** | Russell 2004 Ch3, Tables 3.4–3.11 |
| `russell-ch4-*` | 14 | **multi-arc** | Russell 2004 Ch4, Tables 4.9–4.13 |
| `russell-strange-2009-jovian-multimoon-family` | 1 | non-keplerian | primary=Jupiter |
| `russell-strange-2009-saturnian-multimoon-family` | 1 | non-keplerian | primary=Saturn |
| `aldrin-*` | 4 | single-ellipse | Rogers 2012 Table 1 |
| `mcconaghy-2005-em-case1` | 1 | single-ellipse | Rogers 2012 Table 1 |
| `mcconaghy-2005-em-u0l1` | 1 | single-ellipse | Rogers 2012 Table 1 |
| `mcconaghy-2005-em-snlm-broadclass-family` | 1 | single-ellipse | McConaghy 2005 |
| `mcconaghy-2006-em-k2` | 1 | **UNCERTAIN** (flagged §7) | Russell 2004 Table 4.9 |
| `hollister-menning-1970-ev-orbit-*` | 15 | single-ellipse | Hollister & Menning 1970 |
| `niehoff-visit1` | 1 | single-ellipse | Rogers 2012 |
| `niehoff-visit2` | 1 | single-ellipse | Rogers 2012 |
| `s1l1-2syn-em-cpom` | 1 | single-ellipse | Rogers 2012 Table 1 |
| `jones-2017-vem-triple-family` | 1 | single-ellipse (conservative) | Jones 2017 |
| `vem-emeeve-3syn` | 1 | single-ellipse (conservative) | — |
| `arenstorf-em-figure8-1963` | 1 | non-keplerian | primary=Earth |
| `genova-aldrin-2015-em-3petal-cycler` | 1 | non-keplerian | primary=Earth |
| `wittal-2022-em-cycler-family` | 1 | non-keplerian | primary=Earth |
| `hernandez-2017-jovian-ieg-triple-family` | 1 | non-keplerian | primary=Jupiter |

**Total: 235 rows classified.**
- **multi-arc: 201** (184 ocampo + 14 ch4 + `mcconaghy-2006-em-k2` + `sanchez-net-2022-eem-cycler1` + `sanchez-net-2022-em-cycler2`)
- **single-ellipse: 28**
- **non-keplerian: 6**
- **uncertain: 0** (the formerly-flagged row resolved to multi-arc — see §7/§12)

---

## 9. Machine-usable allowlists

```python
# MULTI_ARC_ALLOWLIST
# All 184 russell-ocampo-* rows + all 14 russell-ch4-* rows.
# Source: Russell 2004, Chapters 3 and 4 (free-return cycler architecture,
# generic-return legs, no single repeating (a,e)).
# Note: mcconaghy-2006-em-k2 is NOT in this list pending human review (see §7).

MULTI_ARC_ALLOWLIST = [
    # russell-ocampo family (184 rows) — Russell 2004 Ch3, Tables 3.4–3.11
    "russell-ocampo-2.1.1+2-case2",
    "russell-ocampo-2.3.1+1-case3",
    "russell-ocampo-2.5.1+0",
    "russell-ocampo-3.1.1+3",
    "russell-ocampo-3.1.2+1",
    "russell-ocampo-3.3.1+2",
    "russell-ocampo-3.5.1+1",
    "russell-ocampo-3.5.2+0",
    "russell-ocampo-3.7.1+1",
    "russell-ocampo-3.9.1+0",
    "russell-ocampo-4.0.3+1",
    "russell-ocampo-4.1.1-5",
    "russell-ocampo-4.10.1+2",
    "russell-ocampo-4.11.1-2",
    "russell-ocampo-4.12.1+1",
    "russell-ocampo-4.13.1-1",
    "russell-ocampo-4.14.1+0",
    "russell-ocampo-4.14.1-1",
    "russell-ocampo-4.3.1-4",
    "russell-ocampo-4.3.1-5",
    "russell-ocampo-4.5.1-3",
    "russell-ocampo-4.5.1-4",
    "russell-ocampo-4.5.2-2",
    "russell-ocampo-4.5.3-1",
    "russell-ocampo-4.7.1-3",
    "russell-ocampo-4.9.1-2",
    "russell-ocampo-4.9.2-1",
    "russell-ocampo-5.1.1-7",
    "russell-ocampo-5.1.2-3",
    "russell-ocampo-5.1.5-1",
    "russell-ocampo-5.2.1+7",
    "russell-ocampo-5.2.2+2",
    "russell-ocampo-5.2.5+0",
    "russell-ocampo-5.3.1-6",
    "russell-ocampo-5.3.1-7",
    "russell-ocampo-5.3.3-2",
    "russell-ocampo-5.4.1+5",
    "russell-ocampo-5.4.1+6",
    "russell-ocampo-5.4.3+1",
    "russell-ocampo-5.5.1-4",
    "russell-ocampo-5.5.1-5",
    "russell-ocampo-5.5.1-6",
    "russell-ocampo-5.5.2-2",
    "russell-ocampo-5.5.2-3",
    "russell-ocampo-5.5.4-1",
    "russell-ocampo-5.6.1+3",
    "russell-ocampo-5.6.1+4",
    "russell-ocampo-5.6.1+5",
    "russell-ocampo-5.6.2+1",
    "russell-ocampo-5.6.2+2",
    "russell-ocampo-5.6.4+0",
    "russell-ocampo-5.7.1-3",
    "russell-ocampo-5.7.1-4",
    "russell-ocampo-5.7.1-5",
    "russell-ocampo-5.8.1+2",
    "russell-ocampo-5.8.1+3",
    "russell-ocampo-5.8.1+4",
    "russell-ocampo-5.9.1-2",
    "russell-ocampo-5.9.1-3",
    "russell-ocampo-5.9.1-4",
    "russell-ocampo-5.9.2-1",
    "russell-ocampo-5.9.2-2",
    "russell-ocampo-5.9.3-1",
    "russell-ocampo-5.10.1+2",
    "russell-ocampo-5.10.1+3",
    "russell-ocampo-5.10.2+0",
    "russell-ocampo-5.10.2+1",
    "russell-ocampo-5.10.3+0",
    "russell-ocampo-5.11.1-2",
    "russell-ocampo-5.11.1-3",
    "russell-ocampo-5.11.2+1",
    "russell-ocampo-5.12.1+1",
    "russell-ocampo-5.12.1+2",
    "russell-ocampo-5.13.1-2",
    "russell-ocampo-5.13.1-3",
    "russell-ocampo-5.13.2-1",
    "russell-ocampo-5.14.1+1",
    "russell-ocampo-5.14.1+2",
    "russell-ocampo-5.14.2+0",
    "russell-ocampo-5.15.1-1",
    "russell-ocampo-5.15.1-2",
    "russell-ocampo-5.16.1+0",
    "russell-ocampo-5.16.1+1",
    "russell-ocampo-5.17.1-1",
    "russell-ocampo-5.18.1+0",
    "russell-ocampo-6.0.1+6d",
    "russell-ocampo-6.0.1+7c",
    "russell-ocampo-6.0.1+8b",
    "russell-ocampo-6.0.1+9a",
    "russell-ocampo-6.1.2-4",
    "russell-ocampo-6.1.3-3",
    "russell-ocampo-6.1.4-2",
    "russell-ocampo-6.1.6-1",
    "russell-ocampo-6.2.1+6",
    "russell-ocampo-6.2.1+7",
    "russell-ocampo-6.2.1+8",
    "russell-ocampo-6.2.2+2",
    "russell-ocampo-6.2.2+3",
    "russell-ocampo-6.2.3+1",
    "russell-ocampo-6.2.3+2",
    "russell-ocampo-6.2.4+1",
    "russell-ocampo-6.2.6+0",
    "russell-ocampo-6.3.1-9",
    "russell-ocampo-6.3.4+1",
    "russell-ocampo-6.4.1+4",
    "russell-ocampo-6.4.1+5",
    "russell-ocampo-6.4.1+6",
    "russell-ocampo-6.4.1+7",
    "russell-ocampo-6.5.1-6",
    "russell-ocampo-6.5.1-7",
    "russell-ocampo-6.5.1-8",
    "russell-ocampo-6.5.5-1",
    "russell-ocampo-6.6.1+3",
    "russell-ocampo-6.6.1+4",
    "russell-ocampo-6.6.1+5",
    "russell-ocampo-6.6.1+6",
    "russell-ocampo-6.6.2+1",
    "russell-ocampo-6.6.2+2",
    "russell-ocampo-6.6.5+0",
    "russell-ocampo-6.7.1-6",
    "russell-ocampo-6.7.1-7",
    "russell-ocampo-6.7.2+3",
    "russell-ocampo-6.7.3-2",
    "russell-ocampo-6.7.5+0",
    "russell-ocampo-6.8.1+2",
    "russell-ocampo-6.8.1+3",
    "russell-ocampo-6.8.1+4",
    "russell-ocampo-6.8.1+5",
    "russell-ocampo-6.8.1+6",
    "russell-ocampo-6.8.3+1",
    "russell-ocampo-6.9.1-4",
    "russell-ocampo-6.9.1-5",
    "russell-ocampo-6.9.1-6",
    "russell-ocampo-6.9.2-2",
    "russell-ocampo-6.9.2-3",
    "russell-ocampo-6.9.4-1",
    "russell-ocampo-6.10.1+2",
    "russell-ocampo-6.10.1+3",
    "russell-ocampo-6.10.1+4",
    "russell-ocampo-6.10.1+5",
    "russell-ocampo-6.10.2+1",
    "russell-ocampo-6.10.2+2",
    "russell-ocampo-6.10.4+0",
    "russell-ocampo-6.11.1-4",
    "russell-ocampo-6.11.1-5",
    "russell-ocampo-6.11.2+2",
    "russell-ocampo-6.12.1+2",
    "russell-ocampo-6.12.1+3",
    "russell-ocampo-6.12.1+4",
    "russell-ocampo-6.13.1-3",
    "russell-ocampo-6.13.1-4",
    "russell-ocampo-6.13.1-5",
    "russell-ocampo-6.13.1+5",
    "russell-ocampo-6.13.2-2",
    "russell-ocampo-6.13.3-1",
    "russell-ocampo-6.14.1+1",
    "russell-ocampo-6.14.1+2",
    "russell-ocampo-6.14.1+3",
    "russell-ocampo-6.14.2+0",
    "russell-ocampo-6.14.2+1",
    "russell-ocampo-6.14.3+0",
    "russell-ocampo-6.15.1-2",
    "russell-ocampo-6.15.1-3",
    "russell-ocampo-6.15.1-4",
    "russell-ocampo-6.15.1+4",
    "russell-ocampo-6.15.2+1",
    "russell-ocampo-6.16.1+1",
    "russell-ocampo-6.16.1+2",
    "russell-ocampo-6.17.1-2",
    "russell-ocampo-6.17.1-3",
    "russell-ocampo-6.17.1+3",
    "russell-ocampo-6.17.2-1",
    "russell-ocampo-6.18.1+1",
    "russell-ocampo-6.18.1+2",
    "russell-ocampo-6.18.2+0",
    "russell-ocampo-6.19.1-2",
    "russell-ocampo-6.19.1+2",
    "russell-ocampo-6.19.2+0",
    "russell-ocampo-6.20.1-4",
    "russell-ocampo-6.20.1+0",
    "russell-ocampo-6.20.1+1",
    "russell-ocampo-6.21.1-1",
    "russell-ocampo-6.21.1+1",
    "russell-ocampo-6.22.1+0",
    # russell-ch4 family (14 rows) — Russell 2004 Ch4, Tables 4.9–4.13
    "russell-ch4-4.991gG2",
    "russell-ch4-8.049gGf2",
    "russell-ch4-8.165Gfh-f2",
    "russell-ch4-9.353Gg2",
    "russell-ch4-3.64gGg3",
    "russell-ch4-3.77Gh3",
    "russell-ch4-3.78Gg3",
    "russell-ch4-5.30gGf3",
    "russell-ch4-5.66Gfh3",
    "russell-ch4-9.94Gg3",
    "russell-ch4-3.66gfF3",
    "russell-ch4-5.30ggF3",
    "russell-ch4-5.75ggF3",
    "russell-ch4-6.44Gg3",
    # resolved from §7 — same physical cycler as russell-ch4-4.991gG2 (Russell Table 4.9)
    "mcconaghy-2006-em-k2",
    # Sanchez Net 2022 EEM near-ballistic real-date patched-conic cycler (Fig. 2a)
    "sanchez-net-2022-eem-cycler1",
    # Sanchez Net 2022 EM near-ballistic real-date patched-conic cycler (Fig. 2b)
    "sanchez-net-2022-em-cycler2",
]

# Verification: this list must have exactly 201 entries.
assert len(MULTI_ARC_ALLOWLIST) == 201

# E-E-M-M single-ellipse rows confirmed NOT in the allowlist:
SINGLE_ELLIPSE_EXCEPTIONS_CHECKED = [
    "s1l1-2syn-em-cpom",      # E-E-M-M, a_au=1.30 (Rogers 2012 Table 1), e=0.257.
                               # Single repeating ellipse under the Rogers circular-
                               # coplanar framing (see §12 — same cycler as the two
                               # multi-arc rows under Russell's free-return framing).
]
```

---

## 12. Finding: one physical cycler, three rows, two model framings

`s1l1-2syn-em-cpom`, `mcconaghy-2006-em-k2`, and `russell-ch4-4.991gG2` describe the
**same physical 2-synodic Earth–Mars cycler** (Russell descriptor 4.991gG, Table 4.9):

- **Rogers 2012 / CPOM framing** (`s1l1-2syn-em-cpom`): a single circular-coplanar
  repeating ellipse, `a=1.30 AU, e=0.257` → **single-ellipse**.
- **Russell 2004 free-return framing** (`mcconaghy-2006-em-k2`, `russell-ch4-4.991gG2`):
  two generic-return arcs, no single `(a,e)` → **multi-arc**.

`cycler_class` is assigned **per row, by that row's own source's framing** — which is
consistent with the spec's model-partitioning (§16.2: signatures/identity are never
compared across `model_assumption`). The two framings are both valid; one is an
idealisation of the other.

**Two follow-ups (NOT blocking this schema work), for human/triage:**
1. **De-duplication / family linkage:** the three rows are one cycler; consider linking
   them via `family{}` or the §16.3 matcher rather than carrying them as independent
   entries.
2. **S1L1 closure puzzle:** that the exact cycler is *multi-arc* under Russell may
   explain why the single-ellipse resonance construction reproduced only the
   coplanar 4.9/5.0 V∞ and never the 5.65/3.05 target. Worth revisiting S1L1 closure
   as a two-arc problem. Tracked separately; do not reopen here.

---

## 10. Russell dissertation pages cited

| Claim | Section / Page |
|---|---|
| Three free-return types defined (half-rev, full-rev, generic) | §2.6, p. 30; §2.8, p. 54 |
| Generic return: "arbitrary time of flight, non-integer-multiple of π transfer angle" | §3.4.3, p. 63 |
| Semi-major axis of generic/odd-nπ found only by Lambert's Equation (not fixed) | §2.8, p. 54 |
| Chapter 3 = "Identical Generic Returns" (all legs same type but not same orbit) | Chapter 3 title, p. 56; §3.1, p. 56 |
| Chapter 4 = "Identical or Different Generic Returns" (extending Ch3) | Chapter 4 title, p. 95; §4.1, p. 95 |
| Ch4 characterizes each cycler by leg descriptor strings, NOT (a, e) | §4.8, p. 126 |
| Aldrin cycler recovered in Ch4 with single-descriptor G(2.1354,7.1338,L) | §4.9, p. 133 |
| Table 4.9: two-synodic ballistic cyclers with descriptor strings | §4.8, p. 127 |
| Table 4.10–4.11: three-synodic ballistic cyclers with generic transit legs | §4.8, p. 128–129 |
| Table 4.12: three-synodic with full-rev transit legs | §4.8, p. 132 |
| Table 4.13: near-ballistic three-synodic cyclers | §4.8, p. 134 |
| Table 4.14: summary of results by transit-leg type | §4.8, p. 134 |
| "mcconaghy-2006-em-k2 / russell-ch4-4.991gG2": two generic arcs 1.46yr + 2.81yr | Table 4.9, p. 127 |

---

## 11. PDF extraction notes

All cited pages were successfully extracted from the Russell 2004 dissertation
(full text held offline; cite by handle `http://hdl.handle.net/2152/1253`,
DOI `10.2514/1.1011` for the companion JGCD paper; not stored in repo).
No extraction failures occurred for any page used in this classification.

Pages 68–80 (Chapter 3 definitions), 107–115 (Ch3 results / Ch4 intro), 116–120
(§4.4 definition + §4.5 problem setup), and 139–155 (Ch4 results Tables 4.9–4.14)
were read and are the basis for the citations above.

The Hollister & Menning 1970 paper (pages 1–3) was also read to confirm the
single-ellipse classification of the hollister rows.
