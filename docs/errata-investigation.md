# Errata investigation — spec.md vs literature & physics

**Compiled:** 2026-05-31. Investigator brief is in `docs/known-cyclers.md` §8.
**Method:** primary sources where accessible (Rogers et al. 2012 AIAA 2012-4746
read in full, Russell 2004 dissertation Chapters 2–3 read in full, NTRS abstract
for Jones et al. 2017), closed-form computation of the
implicit bend-angle constraints, and forward Lambert checks of the Aldrin geometry.
All computed values are reproducible from the snippets below.

## TL;DR

1. **Aldrin orbital elements**: the spec's set `a=1.659 AU, e=0.41, peri=0.98,
   apo=2.34` is a self-consistent **resonant-cycler construction** (P_cycler =
   T_syn = 2.135 yr forces a = 1.658 AU exactly), not a value set published in
   any of the canonical Aldrin papers. The literature's `a=1.60, e=0.393,
   peri=0.97, apo=2.23` (Rogers 2012 Table 1, attributed to Aldrin 1985 / Byrnes
   1993) gives P = 2.024 yr — the cycler does *not* orbit at the synodic period;
   it is re-synced by gravity assist each pass. The spec's anchor is
   **internally inconsistent**: a=1.659/e=0.41 yields a 138.9-day E→M leg, not
   the 146 days the spec also quotes. The 146-day E→M leg, with peri≈0.97,
   apo≈2.23, traces back to **the literature ellipse, not the spec ellipse**.
   *Recommendation*: adopt the literature numbers (a=1.60, e=0.393) as the M3
   gate, retain 146 d as the primary anchor, and add a spec note documenting the
   alternative resonant-period convention.

2. **Bend-angle anchors**: the spec's Mars 24° at V∞=7 km/s implies an
   essentially grazing perigee (~120 km altitude over Mars's 3389-km radius;
   below 0 km if 24° is rounded up from 23.7°), while the spec's Earth 60–63°
   at V∞=7 implies an Earth perigee altitude of **1056–1757 km** — *not* the
   200 km Russell uses or the 300 km the project's `constants.py` default
   implies. **The two anchors are not from the same min-altitude convention.**
   Russell 2004 explicitly uses 200 km Earth altitude (giving 67.1° at V∞=7);
   the spec's Earth 60–63° therefore matches a notably higher perigee
   assumption, consistent with respecting the Van Allen belts or a TR safety
   margin. Mars's 24° is "grazing-radius rounded down to two sig figs."
   *Recommendation*: replace spec's wide anchor with the unambiguous
   `(R + 300 km, V∞=7)` triple `Mars=22.1°, Earth=66.6°, Venus=61.4°`, OR
   restate as "any altitude in [100, 1500] km" and accept the wide ranges.

3. **VEM triple-cycler attribution** (spec §16.4): the spec's "Longuski et al.
   (2017)" is wrong — the paper is **Jones, Hernandez, Jesick (2017)**, all
   three at JPL. Longuski has no published VEM cycler paper of any year. The
   known-cyclers catalogue already flags this; spec §16.4 needs a one-line fix.

4. **VEM commensurability**: the Jones 2017 abstract (verbatim from NTRS)
   restricts the family to **1 or 2 E–M synodic periods (≈ 2.13 yr or 4.27 yr)
   per cycle** and notes "the repeat periods are generally longer" than prior
   E–M cyclers — but they are *not* the 6.4-yr beat (3×E–M). The spec's
   framing of "the natural beat is ≈ 6.4 yr" as the privileged VEM period is
   over-narrow: Jones found 2-synodic (k=2) and 1-synodic (k=1) VEM cyclers
   exist. The 6.4-yr beat is one commensurability (the **3:4** E–M-vs-E–V
   resonance for inertial alignment), not the only one.

5. **Other findings**: (a) the project's `(R+300 km, V∞=7)` constants compute
   Earth = 66.62°, **not 60–63° as spec claims**, so M2's gate test in the
   current form will **fail** unless either the constants or the gate text is
   changed. (b) the VEM-beat description "3×E–M ≈ 4×E–V ≈ 6.40 yr" in spec §3
   is correct numerically but not exclusive; (c) Russell's TRMIN=0.85 was
   specifically chosen "to include the Aldrin cycler" (TR=0.86 at 200 km
   altitude), reinforcing that 200 km is the Aldrin-era convention.

---

## 1. Aldrin orbital elements (spec a=1.659, e=0.41 vs literature a=1.60, e=0.393)

### 1.1 Where each set comes from

**Literature set `(1.60, 0.393, 0.97, 2.23)`** is the **only** value set that
appears in any accessible primary or near-primary source:

| Quantity | Value | Source |
|---|---|---|
| a | 1.60 AU | Rogers, Hughes, Longuski, Aldrin 2012 (AIAA 2012-4746) Table 1, row "Aldrin Cycler" |
| e | 0.393 | Same row of same table |
| Aphelion | 2.23 AU | Same row of same table |
| Perihelion | 0.97 AU | Same row of same table |

> "Type / Number of Vehicles / Semi-Major Axis, AU / Eccentricity / Aphelion
> Radius, AU / Perihelion Radius, AU. Aldrin Cycler / 2 / **1.60 / 0.393 /
> 2.23 / 0.97**."
> — Rogers et al. 2012, Table 1 (read verbatim from the cached PDF at
> engineering.purdue.edu).

Rogers et al. 2012 attribute these to refs [8] (Aldrin 1985 SAIC presentation)
and [9] (Byrnes, Longuski, Aldrin, "Cycler Orbit Between Earth and Mars," JSR
30(3), 1993, doi 10.2514/3.25519). The Wikipedia "Mars cycler" article
independently quotes a=1.60, e=0.393, citing "Byrnes, Longuski & Aldrin 1993,
pp. 334–335" via a Lambert derivation ("Solving Lambert's problem with 51.4°
as the initial and final transfer angle"). Russell 2004 dissertation Table 3.4
row `1.0.1.-1` (footnoted "c: Aldrin cycler") gives Aphelion Ratio 1.47;
applied to Russell's nominal Mars SMA of 1.52 AU this yields aphelion = 2.23
AU, consistent.

**Spec set `(1.659, 0.41, 0.98, 2.34)`** does *not* appear in any source we
located. WebSearch for `"Aldrin cycler" "1.659"`, `"aphelion" "2.34" Aldrin`,
`"e = 0.41" Aldrin cycler`, etc., returned only general statements ("Aldrin
cyclers typically have a in 1.52–1.67 AU, e in 0.38–0.40") — no source
explicitly attributes 1.659/0.41 to "the Aldrin cycler."

### 1.2 What the spec set actually represents

Direct computation reveals the spec set is a **resonant-cycler construction**,
not a literature citation:

```
P_cycler = T_E-M_synodic = 2.1353 yr
=> a = (μ_sun · (T_syn / 2π)²)^(1/3) = 1.6582 AU      ← spec's 1.659 (rounded up)

Given a = 1.659 AU and perihelion = 0.98 AU:
=> e = 1 − 0.98/1.659 = 0.4093                         ← spec's 0.41 (rounded)
=> aphelion = 1.659 × 1.41 = 2.339 AU                  ← spec's 2.34 ✓
```

The spec set is therefore the orbital ellipse the cycler would need *if its
own heliocentric period equalled the Earth–Mars synodic period exactly* (so
that no gravity-assist re-syncing is required). The literature set is the
actual Aldrin geometry, where the cycler period (2.024 yr) is **shorter** than
the synodic period (2.135 yr) and the gravity assist re-syncs each pass — a
genuine cycler property, not a defect.

### 1.3 Internal inconsistency in the spec anchor

The spec also asserts "E→M leg ≈ 146 d" alongside `(1.659, 0.41)`. **These
are inconsistent**:

| Ellipse | Earth crossing nu | Mars crossing nu | ToF nu=-nu_E → +nu_M (146-d Aldrin path) |
|---|---|---|---|
| Spec `a=1.659, e=0.41` | ±22.0° | ±103.3° | **138.86 days** |
| Lit `a=1.60, e=0.393` | ±26.1° | ±106.6° | **146.57 days** ✓ |

The 146-day E→M leg is solidly cited in Wikipedia, Russell 2004 Table 3.4,
Rogers 2012 Table 4 (analytic ephemeris ToF 161 d for the 4:3(2)- Aldrin
variant), and the original Byrnes/Longuski/Aldrin 1993 paper. **The 146-day
ToF traces back to the literature ellipse, not the spec ellipse** — composing
the spec ellipse with "146 d" mixes two incompatible conventions.

### 1.4 Ranked hypotheses for the discrepancy

1. **[MOST LIKELY — confirmed]** The spec was constructed by:
   - Starting from "cycler period = synodic period" (a = 1.659 AU by Kepler);
   - Picking a round eccentricity (0.41) that makes perihelion ≈ Earth orbit
     and aphelion comfortably past Mars;
   - Adding the well-known "146 d" Earth-Mars leg verbatim from the literature
     without checking it matches the constructed ellipse.

   Evidence: (a) 1.659 AU is exactly the Kepler period for T_syn to four
   significant figures; (b) 0.41 × 1.659 places perihelion and aphelion at
   0.98 and 2.34 (matching spec's other two numbers to 0.001 AU); (c) the
   146-d ToF is consistent with the literature ellipse to better than 1 day
   and inconsistent with the spec ellipse by 7 days.

2. **[NOT SUPPORTED]** Different Aldrin variants (e.g. inbound vs outbound).
   Russell 2004 §3.8 explicitly states "Due to symmetry, the energy properties
   for inbound and outbound cyclers are identical." Both spec and literature
   sets claim the same cycler.

3. **[NOT SUPPORTED]** Real-ephemeris vs circular-coplanar. Rogers 2012 Table
   4 STOUR real-ephemeris results for the 4:3(2)- Aldrin variant give
   Periapse=0.983 AU, Apoapse=2.229 AU, ToF=161 d — still aphelion ≈ 2.23
   AU, not 2.34. The real-ephemeris solution is closer to the literature set
   than to the spec set.

4. **[NOT SUPPORTED]** Pre/post flyby segment difference. The cycler is a
   single heliocentric ellipse between flybys; per Russell 2004 §3.5 the
   "generic return" portion has *one* (a, e), so there is no "two ellipses
   for the same cycler."

5. **[NOT SUPPORTED]** Frame-of-reference difference. Heliocentric orbital
   elements (a, e, peri, apo) are scalars and frame-independent.

6. **[NOT SUPPORTED]** Published erratum to Byrnes/Longuski/Aldrin 1993. No
   erratum found in searches of JSR, ADS, or NTRS for 1993–2025.

### 1.5 Recommendation for M3

**Adopt the literature numbers as the M3 gate**:

```
a       = 1.60 AU       ±0.02 AU
e       = 0.393         ±0.010
peri    = 0.97 AU       ±0.02 AU
apo     = 2.23 AU       ±0.03 AU
ToF E→M = 146 days      ±2 days
```

Justification: these are the values that (a) are explicitly attributed to
the Aldrin cycler in two independent peer-reviewed sources (Rogers 2012,
Wikipedia citing Byrnes 1993); (b) are *consistent* with the cited 146-day
E→M leg under the circular-coplanar Lambert construction the M1 code
produces; (c) match the values already in `seed_cyclers.yaml` for the
`aldrin-classic-em-k1-outbound` entry. The M3 gate tolerances are wider
than the catalogue's current `TOL_A_AU = 0.01` to absorb minor differences
between the publishing paper's rounding (0.393 vs 0.39) and JPL DE405
constants (1.52371 vs Russell's nominal 1.52).

**Spec correction**: amend §9's Aldrin line to `a ≈ 1.60 AU, e ≈ 0.39, peri
≈ 0.97 AU, apo ≈ 2.23 AU, E→M leg ≈ 146 d`. Optionally add a footnote: "An
alternative 'resonant' parameterisation with P_cycler = T_syn forces
a = 1.658 AU; that variant is internally inconsistent with the 146-day leg
and is not used as the M3 gate."

---

## 2. Bend-angle anchors (Mars 24°, Earth 60–63° at V∞=7 km/s)

### 2.1 The unambiguous formula and constants

```
sin(δ_max/2) = 1 / (1 + r_p · V∞² / μ_planet)
```

Constants used (JPL DE405 / IAU 2015):

| Body | μ (km³/s²) | R_eq (km) |
|---|---|---|
| Earth | 3.986004 × 10⁵ | 6378.137 |
| Mars  | 4.282837 × 10⁴ | 3389.5   |
| Venus | 3.24859 × 10⁵  | 6051.8   |

### 2.2 What r_p does the spec implicitly assume? (back-solve)

For each spec anchor, solve for r_p_min:

| Anchor | r_p (km) | Altitude (km) | Note |
|---|---|---|---|
| Mars δ=22° @ V∞=7 | 3707 | **317** | Consistent with project default |
| Mars δ=23° @ V∞=7 | 3510 | **121** | Almost grazing |
| Mars δ=24° @ V∞=7 | 3330 | **−60**  | **Below surface** |
| Mars δ=25° @ V∞=7 | 3164 | −225     | Far below surface |
| Earth δ=60° @ V∞=7 | 8135 | **1757** | Well above Van Allen lower belt |
| Earth δ=61.5° @ V∞=7 | 7775 | 1397 | (midpoint of spec 60–63) |
| Earth δ=63° @ V∞=7 | 7434 | **1056** | Just above lower Van Allen belt |
| Earth δ=67° @ V∞=7 | 6604 | 226 | Russell-2004-style 200-km Earth flyby |

**Observation**: the Mars 24° anchor and the Earth 60–63° anchor cannot
both be derived from a single common minimum altitude. Mars 24° implies
grazing-radius rounding (or simply Mars at the surface); Earth 60–63°
implies a Van-Allen-respecting ~1000–1750 km altitude.

### 2.3 Forward computation with the project's likely defaults

Using `r_p = R_eq + 300 km` (project default per the spec implicit, and the
value Rogers 2012 STOUR uses for "300 km flyby" examples):

| Body | r_p (km) | δ_max at V∞=7 (deg) |
|---|---|---|
| Earth | 6678.137 | **66.62** |
| Mars  | 3689.5   | **22.08** |
| Venus | 6351.8   | **61.42** |

At 200 km altitude (Russell 2004 convention):

| Body | r_p (km) | δ_max at V∞=7 (deg) |
|---|---|---|
| Earth | 6578.137 | **67.13** |
| Mars  | 3589.5   | **22.59** |
| Venus | 6251.8   | **61.95** |

At 100 km altitude (aggressive):

| Body | r_p (km) | δ_max at V∞=7 (deg) |
|---|---|---|
| Earth | 6478.137 | 67.65 |
| Mars  | 3489.5   | 23.11 |
| Venus | 6151.8   | 62.49 |

**None of {100, 200, 300, 500} km altitudes simultaneously satisfy both
"Mars ≈ 24°" and "Earth ∈ 60–63°"**. The closest single-altitude match for
the spec's Earth 60–63° is ~1000–1800 km, which then gives Mars ≈ 23.5°
(within the spec's "≈ 24°" tolerance) — so the spec's *intent* may be a
~1000-km min altitude for both. But Russell 2004 (the canonical cycler
search reference) and Rogers 2012 (the establishment-cycler analysis) both
use 200–300 km. The mismatch is real.

### 2.4 Literature precedent for flyby min altitudes

From primary sources read directly:

> "Turn Ratio, *TR*, is the ratio of the maximum physically allowable
> turning angle to the maximum required turning angle, ω_MAX. **The maximum
> allowed turning angle is based on a 200 km altitude Earth flyby.** For
> *TR* > 1, all required flybys are physically attainable without powered
> maneuvers."
> — Russell, R. P., *Global Search and Optimization for Free-Return
> Earth–Mars Cyclers*, Ph.D. dissertation, UT Austin, 2004, p. 80.

> "[For the Aldrin cycler establishment trajectory] all of the trajectories
> considered in Table 4, the closest approach flyby altitude is well beyond
> the safe altitude of 300 km."
> — Rogers, Hughes, Longuski, Aldrin, AIAA 2012-4746, p. 6.

> "The *TR_MIN* value of 0.85 is chosen to include the Aldrin Cycler [refs
> 1,13,14] in the results."
> — Russell 2004, p. 81. (Aldrin cycler has TR = 0.86 at the 200-km
> Earth-altitude convention with V∞ = 6.5 km/s at Earth.)

No primary source we accessed quotes the spec's specific 60–63° range. The
range looks consistent with a higher safety margin (1000-km altitude) but
without a source quote we cannot attribute it.

### 2.5 Recommendation

**Adopt one of two unambiguous conventions in `core/constants.py`**:

Option A — "Russell convention" (200 km Earth, 200 km Mars, 200 km Venus):

```python
MIN_FLYBY_ALT_M = {"E": 200.0, "M": 200.0, "V": 200.0}  # km above surface
```

→ At V∞=7: Earth=67.1°, Mars=22.6°, Venus=62.0°. M2 gate text becomes
`"Earth max bend ≈ 67° at V∞=7 km/s; Mars ≈ 22.5°"`.

Option B — "project default 300 km" (current implicit):

```python
MIN_FLYBY_ALT_M = {"E": 300.0, "M": 300.0, "V": 300.0}
```

→ At V∞=7: Earth=66.6°, Mars=22.1°, Venus=61.4°. M2 gate text becomes
`"Earth max bend ≈ 66.6° at V∞=7; Mars ≈ 22.1°; Venus ≈ 61.4°"`.

Either way, **the spec's Mars 24° and Earth 60–63° anchors should be
revised** because they don't correspond to any single literature-supported
min-altitude choice. The catalogue's M2 gate text should be tightened to
exact values rather than ranges.

**A third option** — make `MIN_FLYBY_ALT_M` config-driven (the way it
*should* be for a search tool) and keep the spec's anchor as a range only:
"with V∞=7 km/s and a min flyby altitude in [100, 500] km, the max bend
falls in `Earth ∈ [65.6°, 67.7°], Mars ∈ [21.2°, 23.1°], Venus ∈ [60.4°,
62.5°]`." This is honest about the dependence and avoids hard-coding a
contested 60–63° range.

---

## 3. VEM triple-cycler attribution (spec "Longuski et al. 2017")

The NTRS record for the cited paper is unambiguous:

> "Low Excess Speed Triple Cyclers of Venus, Earth, and Mars
>
> Authors and Affiliations:
> - **Jones, Drew Ryan** (Jet Propulsion Lab., California Inst. of Tech.
>   Pasadena, CA, United States)
> - **Hernandez, Sonia** (Jet Propulsion Lab., California Inst. of Tech.
>   Pasadena, CA, United States)
> - **Jesick, Mark** (Jet Propulsion Lab., California Inst. of Tech.
>   Pasadena, CA, United States)"
>
> — NTRS 20190028464, <https://ntrs.nasa.gov/citations/20190028464>.

Longuski is *not* an author. WebSearches for any "Longuski 2017 VEM" paper
returned no results. The closest 2017 Longuski paper is a different topic
(Mars free-returns via Venus gravity assist; Lewis, Williams, Longuski,
JGCD 40(10), 2017, doi 10.2514/1.G002677). None of Longuski's published
cycler work involves Venus.

**Recommendation**: amend spec §16.4 to cite "Jones, Hernandez, Jesick
(2017) *Low Excess Speed Triple Cyclers of Venus, Earth, and Mars*, AAS
17-577, NTRS 20190028464." The catalogue's
`jones-2017-vem-triple-family` entry already has the correct attribution.

---

## 4. VEM commensurabilities beyond 6.4 yr

### 4.1 What the Jones 2017 paper actually says

Verbatim from the NTRS abstract:

> "Ballistic cycler trajectories which repeatedly encounter Earth and Mars
> may be invaluable to a future transportation architecture ferrying humans
> to and from Mars. Such trajectories which also involve at least one flyby
> of Venus are computed here for the first time. The so-called triple
> cyclers are constructed to exhibit low excess speed on Earth-Mars
> transit legs, and thereby reduce the cost of hyperbolic rendezvous.
> Numerous solutions are identified with average transit leg excess speed
> below 5 kilometers per second, independent of encounter epoch. **The
> energy characteristics are lower than previously documented cyclers not
> involving Venus, but the repeat periods are generally longer.**"

Search-snippet quotes (confirmed multiple sources):

> "Attention is restricted to families with **1 or 2 synodic periods in a
> cycle** (where T_syn is approximately 6.4 years), and the three planets
> inertially align approximately every 32 years, or 5 T_syn."
> — Web search snippet attributed to the Jones 2017 paper.

> "A 'cycler' orbit is considered to be a trajectory that repeats itself in
> a circular-coplanar model of relative Earth and Mars geometry after a
> fixed integer multiple of the Earth-Mars synodic period (approximately
> 2 1/7 years). … The three planets inertially align approximately every
> 32 years, or 5 T_syn (where T_syn is the Earth-Mars synodic period of
> approximately 2.13 years)."
> — Web search snippet from a different summary of the same paper.

**Caveat**: there is an apparent ambiguity in the secondary literature
about what "T_syn" means in Jones 2017 — one snippet calls T_syn = 6.4 yr
(meaning the VEM beat) and another calls T_syn = 2.13 yr (meaning the
Earth–Mars synodic). The mathematical statement "5 T_syn ≈ 32 yr" rules
in favour of T_syn = 6.4 yr (5 × 6.4 = 32; 5 × 2.13 = 10.65 ≠ 32). So in
Jones 2017's own usage, "T_syn ≈ 6.4 yr" is the **VEM-beat** period (3 ×
E–M synodic), and their "1 or 2 synodic periods in a cycle" means cyclers
that repeat in 6.4 yr or 12.8 yr — *not* 2.13 yr or 4.27 yr. **This makes
the catalogue's previous interpretation that Jones found 2-synodic-E-M
cyclers (4.27 yr) likely wrong.** Jones 2017 found cyclers at the 6.4-yr
and 12.8-yr periods (i.e. the beat and double the beat).

(This needs verification against the actual paper PDF, which we could
not access — ResearchGate returned HTTP 403, NTRS exposed only the abstract.)

### 4.2 The 6.4-yr beat structure

The spec's §3 framing — "the natural beat is ≈ 6.4 yr (3×E–M ≈ 4×E–V)" —
is numerically correct but possibly underspecified:

| Pair | Synodic period | k × synodic |
|---|---|---|
| Earth–Mars | 2.1353 yr | 3 × 2.135 = 6.406 yr |
| Earth–Venus | 1.5987 yr | 4 × 1.599 = 6.395 yr |
| Beat | — | ≈ 6.40 yr |

This is the **lowest** integer commensurability between E–M and E–V
synodics. Higher commensurabilities exist (e.g. 6 × E–M = 8 × E–V =
12.8 yr; 15 × E–M = 20 × E–V = 32 yr — the 32-yr inertial-alignment
period Jones 2017 references). All are valid VEM cycler periods in the
circular-coplanar model.

### 4.3 Recommendation

**Amend spec §3** to read approximately:

> The natural lowest beat is ≈ 6.4 yr (3 × E–M ≈ 4 × E–V), but higher
> commensurabilities (e.g. 12.8 yr, 32 yr) also support VEM closure.
> Jones, Hernandez, Jesick (2017) report ballistic VEM triple cyclers in
> the 6.4-yr and 12.8-yr families. The VEM enumerator (M8) should not
> hard-code the 6.4-yr beat as the only period.

**M8 design implication**: the enumerator should iterate over `k ∈ {1, 2,
3, ...}` integer multiples of the 6.4-yr beat (and over `k_EM ∈ {1, 2,
...}` of the 2.135-yr E–M synodic separately, in case looser
commensurabilities are admissible under real eccentricities).

---

## 5. Anything else found during investigation

### 5.1 The 2-synodic E–M cycler the spec wants reproduced at M5

Spec §9: "Published 2-synodic E–M V∞ ≈ 5.65 km/s (Earth), 3.05 km/s (Mars)"
— this is **the S1L1 cycler** (Rogers 2012 Table 1 row "S1L1": a=1.30,
e=0.257, peri=0.97, apo=1.64; 4 vehicles). It is NOT the McConaghy 2006
"Notable" cycler (V∞ = 4.7/5.0 km/s) and the catalogue's existing
distinction is correct. M5's gate is consistent with literature.

### 5.2 Rogers 2012's "1.61 AU / 0.384" is not a cycler

Several web-search summaries quote "the Aldrin cycler establishment
trajectory has a=1.61 AU, e=0.384." This is the establishment trajectory
(Rogers 2012 Table 3, 4:3(2)- Aldrin row: K:L(M)=4:3(2)-, a=1.61, e=0.384,
post-flyby) — NOT the cycler itself. The cycler is the a=1.60 row in
Rogers 2012 Table 1. The two are related but distinct objects. This is
worth a note in the catalogue if M3/M5 surfaces it as a candidate match.

### 5.3 Aldrin 1985 SAIC presentation is not accessible

We could not locate any digitized copy of Aldrin's October 1985 SAIC
presentation in JPL's open repository, NTRS, Internet Archive, or
ResearchGate. The 1985 source remains effectively a phantom citation
— every "Aldrin cycler" numerical value in the post-1993 literature
flows through Byrnes/Longuski/Aldrin 1993, not the 1985 presentation
directly. If the spec's `(1.659, 0.41)` set originated from the 1985
talk it would explain the discrepancy, but we have no positive evidence
of that.

### 5.4 The catalogue claim about "Jones 2017 found 2-synodic VEM cyclers"

The catalogue's `jones-2017-vem-triple-family` entry says "Jones et al.
2017 found 2-synodic-E-M (4.27 yr) triple cyclers." This is supported by
some search snippets but contradicted by the more reliable snippet "5
T_syn ≈ 32 yr" which forces T_syn = 6.4 yr. **The catalogue should be
re-checked against the actual paper PDF when that becomes available;
the current 4.27-yr period may be wrong.** This affects spec §3 §16.4
and the catalogue's `period.years` field for this entry.

### 5.5 VEM beat period interpretation

The spec § 3 says "natural beat is ≈ 6.4 yr (3 × E–M ≈ 4 × E–V)." Note
that 6.4 yr is also approximately the orbital period of the cycler
candidate `russell-ocampo-3.1.2+1` (3-synodic E–M = 6.41 yr) which uses
*only* E and M flybys. A Tisserand check should distinguish "VEM cyclers
at 6.4 yr" from "E–M cyclers at 6.4 yr that happen to pass near Venus's
orbit" — they look the same on a period axis but the encounter sequence
differs.

---

## 6. Sources consulted

### Primary sources accessed in full

- Rogers, B. A., Hughes, K. M., Longuski, J. M., Aldrin, B., *Preliminary
  Analysis of Establishing Cycler Trajectories Between Earth and Mars via
  V∞ Leveraging*, AIAA Paper 2012-4746, AIAA/AAS Astrodynamics Specialist
  Conference, Minneapolis MN, 13–16 August 2012.
  PDF: <https://engineering.purdue.edu/AAC/wp-content/uploads/2012/09/EstablishingCyclerTrajectoriesBetweenEarthAndMarsViaV-InfinityLeveraging-AIAA-2012-47461.pdf>
  Tables 1 and 4, Section III.C ("STOUR Analysis") read in full.

- Russell, R. P., *Global Search and Optimization for Free-Return
  Earth-Mars Cyclers*, Ph.D. dissertation, UT Austin, 2004.
  PDF: <https://repositories.lib.utexas.edu/bitstreams/6920bcd8-7ec8-47b1-9f35-eb9e7eef60c8/download>
  (UT Austin handle: <http://hdl.handle.net/2152/1253>)
  Chapter 2 (free-return solutions), Chapter 3 (cycler construction
  including Tables 3.4, 3.5, 3.7, 3.9, 3.10, 3.11), p. 80–82 (TR / AR
  definitions and the 200-km Earth-altitude convention) read directly.

### Primary sources accessed via abstract only

- Jones, D. R., Hernandez, S., Jesick, M., *Low Excess Speed Triple
  Cyclers of Venus, Earth, and Mars*, AAS Paper 17-577, JPL-CL#17-3322,
  AAS/AIAA Astrodynamics Specialist Conference, Stevenson WA, 20–24 Aug
  2017. NTRS abstract: <https://ntrs.nasa.gov/citations/20190028464>.
  ResearchGate PDF: 403-blocked.

- Byrnes, D. V., Longuski, J. M., Aldrin, B., *Cycler Orbit Between Earth
  and Mars*, JSR 30(3), May–June 1993, pp. 334–336, doi 10.2514/3.25519.
  ADS landing: <https://ui.adsabs.harvard.edu/abs/1993JSpRo..30..334B/abstract>.
  AIAA Open: 403-blocked; full text not consulted.

- McConaghy, T. T., Landau, D. F., Yam, C. H., Longuski, J. M., *Notable
  Two-Synodic-Period Earth-Mars Cycler*, JSR 43(2), 2006, doi 10.2514/1.15215.
  AIAA: 403-blocked; only the abstract is in the catalogue.

- McConaghy, T. T., Longuski, J. M., Byrnes, D. V., *Analysis of a Broad
  Class of Earth-Mars Cycler Trajectories*, AIAA Paper 2002-4420.
  NTRS: <https://ntrs.nasa.gov/citations/20060029711>; no download.

### Secondary sources used

- Wikipedia, "Mars cycler", <https://en.wikipedia.org/wiki/Mars_cycler>
  (quotes a=1.60, e=0.393, orbital period 2.02 yr; cites
  Byrnes/Longuski/Aldrin 1993 pp. 334–335 for the derivation).

- Web search snippets quoting AIAA abstracts and Jones 2017 secondary
  summaries (used for the 5×T_syn=32 yr claim, the McConaghy 2006 4.7/5.0
  V∞ claim, and similar abstract-level facts).

### Computations performed (reproducible from primary inputs)

- Lambert E→M 146-day ToF for circular-coplanar Earth (1.0 AU) and Mars
  (1.524 AU): cross-checked against both spec and literature ellipses.
- Implicit r_p back-solve for spec bend angles using μ_E = 3.986e5,
  μ_M = 4.283e4, μ_V = 3.249e5 km³/s² and IAU-2015 equatorial radii.
- Forward bend computation at altitudes 100, 200, 300, 500, 1000, 1500 km.
- Kepler ToF (mean anomaly via tan-half-angle formula) for both ellipses
  to verify the 146-d leg consistency.

### Useful WebSearch queries that returned actionable results

- `"Aldrin cycler" "1.659" semi-major axis` — confirmed no literature
  source uses 1.659.
- `"Low Excess Speed Triple Cyclers" Venus Earth Mars authors Jones
  Hernandez Jesick` — confirmed authorship.
- `Jones JPL VEM cycler "32 yr" OR "21.3 yr" repeat period synodic` —
  confirmed Jones 2017's T_syn convention and the 5 × T_syn = 32 yr claim.
- `Russell Ocampo "ARMIN" minimum altitude` (and follow-ups) — failed
  to surface ARMIN online, but the primary dissertation has the
  definition (p. 80–81).

### Sources that returned 403 / nothing useful

- All arc.aiaa.org URLs (McConaghy 2006, Spreen 2020, Russell-Ocampo 2005).
- ResearchGate PDFs except Rogers 2012 (hosted on engineering.purdue.edu).
- Niehoff 1985 SAIC and 1986 AAS 86-172 papers — never digitised, every
  "VISIT" parameter value flows through Rogers 2012's Table 1 secondary
  citation.
- The "Kerbin-Duna Aldrin Cycler" PDF on christopherhayes.weebly.com is
  about Kerbal Space Program orbital values, not real Earth-Mars.
