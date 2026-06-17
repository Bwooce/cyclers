# Gravier-Marchal-Culp 1972 — Optimal Trajectories Between Earth and Mars in Their True Planetary Orbits

**Deep-read verdict note, 2026-06-17 AET.** Supersedes the file-only entry in
`2026-06-17-mars-cycler-wave-digest.md`.

## Header

- **Title:** *Optimal Trajectories Between Earth and Mars in Their True Planetary Orbits*
- **Authors:** J. P. Gravier (Doctoral Candidate, U Colorado / Boulder),
  C. Marchal (Visiting Assoc. Prof., U Colorado; permanently at ONERA,
  Châtillon, France), R. D. Culp (Asst. Prof., U Colorado)
- **Venue:** Journal of Optimization Theory and Applications, Vol. 9, No. 2,
  pp. 120-136 (1972)
- **DOI:** 10.1007/BF00932349
- **Pages:** 17 (paper proper: pp. 120-136)
- **Received:** September 14, 1970. Supported by NASA grant NGR-06-003-033.

## What the paper actually is

A foundational optimal-control paper applying **Lawden's primer vector
theory** to the Earth-Mars two-impulse (and three-impulse-correction)
rendezvous problem in the **true** planetary orbits — i.e. accounting for
the small but non-zero eccentricities (e_Earth = 0.0167, e_Mars = 0.0934)
and the relative inclination between Earth's and Mars's orbital planes.

Method:

1. The Hohmann transfer is parameterised by the longitude of perihelion of
   the Hohmann window, ω̄_H. Tables 1 and 2 give 16 Earth-to-Mars and 17
   Mars-to-Earth Hohmann windows from 1971 to 2003 with their ω̄_H values
   (0°-360°).
2. For each Hohmann window, two optimal noncoplanar transfers exist: a
   **short transfer** (true anomaly arc less than 180°) and a **long
   transfer** (greater than 180°). Each is computed iteratively via
   primer-vector necessary conditions (Eqs. 3-9), with check for whether
   the primer-vector magnitude exceeds unity (in which case a midcourse
   third impulse is needed).
3. Results are presented as continuous functions of ω̄_H over 0°-360°
   (Figs. 2-16), making them as easy to use as Hohmann transfers but
   accounting for true-orbit effects.

Quantitative anchors:

- Earth: a = 1 AU, V_p = 11.13 km/s, e = 0.0167272, ⟨V⟩ = 29.77 km/s.
- Mars: a = 1.523691 AU, V_p = 5.0 km/s, e = 0.0933654, ⟨V⟩ = 24.117 km/s.
- Idealized Hohmann: a_H = 1.26185 AU, e_H = 0.20751, V_perihelion =
  32.713 km/s, V_aphelion = 21.470 km/s, ΔV_peri = 11.513, ΔV_aph =
  5.658. C_idealized = 1.031 km/s (parabolic-reference characteristic
  velocity). TOF_H = 0.7085 yr = 258.8 d.
- **Earth-to-Mars 1984 worked example** (ω̄_H = 160°, Mar 2 departure):
  exact 2-impulse short C = 1.85 km/s (vs Hohmann 1.03 and first-order 1.7);
  Δv₁ = 1.3 km/s, Δv₂ = 0.55 km/s; transfer 0.56 yr, true anomaly arc 147°,
  arrival 60 d earlier than Hohmann.
- Inclination effects: transfer-plane inclination relative to Earth orbit
  ranges 0-5° depending on ω̄_H (Fig. 13).

The 16 cited Earth-to-Mars Hohmann windows (Table 1) include classic dates
1971-05-11 (the year of the paper's writing), 1975-08-18, 1990-07-29,
1998-12-24, ... 2001-04-01.

## Catalogue / KNOWN_CORPUS relevance

**NOT a cycler paper. NOT a quasi-cycler paper. NOT in any of the four
v4.7 classes (cycler / quasi_cycler / precursor_mga / mga_tour).**

This is a one-way 2-impulse / 3-impulse rendezvous-trajectory paper using
Lawden primer-vector optimization. The transfers depart Earth, do a single
Lambert-like arc, and arrive at Mars (or reverse). There is no repeated
encounter, no resonance, no period, no Tisserand-pump structure.

**Historical importance:** This is one of the foundational papers in the
"true planetary orbits" branch of Earth-Mars trajectory optimization. It
predates the Aldrin cycler concept (Friedlander-Niehoff-Byrnes-Longuski
AIAA-86-2009, then Byrnes-Longuski-Aldrin JSR 1993) by 13-21 years and
predates Walker-Ireland 1971 / Hollister 1969 cycler concepts by 1-3 years.

The paper's primer-vector technique is the same mathematical machinery
later used in Lawden's *Optimal Trajectories for Space Navigation* and is
applied throughout the cycler literature for the *insertion* and
*establishment* problems (e.g. Rogers-Hughes-Longuski-Aldrin 2015 ACTA
"Establishing cycler trajectories between Earth and Mars").

### Recommended treatment

**No new KNOWN_CORPUS anchor.** Reference-only as a methodology citation
for primer-vector optimal-impulse work. If a future write-up needs to cite
the canonical "true planetary orbits Hohmann analysis with primer vector"
paper, this is it (alongside Marchal's ONERA Publication No. 124, 1968 and
Gravier's 1970 CU Boulder PhD thesis — both cited as references [2] and
[3] of this paper).

The Hohmann-window date tables (Tables 1, 2) might be a sourced sanity
check for any catalogue code that computes Earth-Mars opposition or
Hohmann-window dates over the 1971-2003 range. The ω̄_H column (longitude
of perihelion of Hohmann window) is a useful cross-check parameter.

## Errata vs the pre-read survey

The pre-read survey filed this as "*The 1972 JOTA paper on optimal
Earth-Mars in their true planetary orbits. Foundational, predates cycler
concept. File-only as historical reference.*" The deep-read **confirms**
this verdict in full. Two minor refinements:

1. The paper predates the *Aldrin* cycler concept (1985-1993) by 13-21
   years but is **roughly contemporary** with the original
   Hollister 1969 / Walker-Ireland 1971 cyclers; the survey's
   "predates cycler concept" is mostly correct (Hollister 1969 is the
   conventional first cycler).
2. The primer-vector technique IS used in the modern cycler-establishment
   literature (Rogers-Hughes 2015, etc.), so the methodology has been
   absorbed — citation chain is intact.

## Action items for parent

1. **No new KNOWN_CORPUS anchor.** Reference-only.
2. **No catalogue impact.** No SnLm / Aldrin / VEM / Trojan / DRO row
   needs amendment.
3. **Methodology citation pool:** if any future cycler-establishment work
   needs to cite Lawden primer-vector / true-orbit optimal transfer as
   foundational methodology, add Gravier-Marchal-Culp 1972 (DOI
   10.1007/BF00932349) to the citation pool alongside Lawden's textbook.
4. **Sourced sanity check (optional):** Tables 1-2 (1971-2003 Hohmann
   windows for E→M and M→E with longitudes of perihelion ω̄_H) could be a
   golden-test reference for any code that computes E-M opposition dates
   pre-2003. Low priority; only useful if the catalogue ever produces an
   E-M Hohmann-equivalent row.
5. **No acquisition wishlist additions.** The cited Marchal 1968 ONERA
   pub and Gravier 1970 PhD are NOT cycler papers (single-leg transfer
   technique only).
