# `#885` — Bridge from a published Uranus tour end-state into a catalogued quasi-cycler

Task `#885` (ledger entry in `data/OUTSTANDING.md`). Results note. Sections 1 and 2 were written and
committed BEFORE any bridge or control was computed; the commit that adds them is the
pre-registration record. Later sections are appended below them, and Sections 1 and 2 are not
edited after results exist (a later opinion on them goes in a separate section).

## 1. Pre-registered threshold (written before computing any bridge)

**What had been computed when this was written.** Only the circular-coplanar Tisserand points
(semi-major axis, eccentricity, periapsis, apoapsis) implied by the published V-infinity pairs and by
the catalogued rows, plus per-flyby maximum bend at 50 km. Those show the size of the V-infinity gap
but not what a bridge costs. No leveraging leg, Lambert leg, control or search had been run.

**The bridge is "worth a row" (a `precursor_mga` row with `inserts_into` the quasi-cycler) if BOTH:**

- **total deterministic delta-v <= 150 m/s**, counted from the last published flyby to the first
  encounter of the quasi-cycler, including any entry correction needed beyond the bend a flyby can
  supply; and
- **duration <= 365 days** over the same span.

**A "strong" bridge** (worth stating in the row as an easy add-on) is one with total delta-v <= 60 m/s
and duration <= 180 days.

**Per-manoeuvre limit (a constraint, not part of the verdict):** no single manoeuvre above 60 m/s
(Landau, Davis & Karimi 2023, Table 2, "max. single delta-V 60 m/s"). If the best bridge needs a
larger single burn, it is reported as violating that published limit, and the verdict is still
judged on the two totals above.

**Derivation.** The bridge is an extra phase appended to an existing tour, so it is measured
against what the published Uranus Orbiter and Probe tours pay for comparable things:

- McAdams et al. 2011 (AAS 11-188), Table 5: targeting delta-v to add each new moon to the tour is
  89.0, 102.8, 103.0, 140.9 and 183.3 m/s (619 m/s for the whole 424-day tour). Adding one new
  target to a published tour costs 89-183 m/s; 150 m/s is inside that range, about one such leg.
- Landau et al. 2025 (IEEE Aerospace), Table 7: moon-tour delta-v2 about 0.73-0.78 km/s including
  statistical, margin and attitude control. 150 m/s is about a fifth of a whole moon tour's budget,
  and about 5 percent of the ~2.8 km/s total mission delta-v of the chosen design.
- Landau, Davis & Karimi 2023 (AAS 23-460), Figure 6: the in-tour manoeuvres between the ten moon
  flybys sum to 249 m/s by my addition (23, 59, 15, 19, 15, 8, 5, 45, 60). 150 m/s is about 60 percent
  of that tour's inter-flyby spending; the "strong" level, 60 m/s, is the published single-manoeuvre
  cap, i.e. a bridge that could be flown as one maximum-size published burn in total.
- Ellison et al. 2025 (AAS 25-668): 13.37 m/s deterministic for the 1.69-year science tour, 18.32 m/s
  for the whole 3.17-year moon phase. This tour shows the published designs CAN be nearly ballistic;
  a 150 m/s bridge would be about eight times that tour's deterministic spending, which is why
  150 m/s and not more is the upper edge.
- Duration: Landau 2023 caps the tour at 500 days since orbit insertion; Ellison's science tour is
  1.69 years; Landau 2025's example tour spans 703.5 days of flybys and its magneto-tour adds about
  2 years. A bridge longer than a year is comparable to a whole published tour phase, so 365 days is
  the upper edge, and 180 days (about a quarter of the Landau 2025 tour span) the "strong" level.

## 2. Pre-registered positive control (written before computing it)

The leveraging primitive used for the bridge is Landau, Davis & Karimi's own tour step (their
seven-step algorithm): ballistic flyby at the minimum altitude, Kepler propagation to the next
apoapsis, an impulsive manoeuvre there, and a Lambert arc of at most one revolution to the next moon.
The control applies that same code to three equatorial legs of their Figure 6 (after the
inclination-reduction phase), with the 100 km altitude of that source:

| Leg | Published | delta-v | Flight time |
|---|---|---|---|
| C1 | Ariel 3.990 -> Ariel 3.897 km/s | 15 m/s | 26 d (03/04 -> 03/30/2050) |
| C2 | Umbriel 4.095 -> Umbriel 4.070 km/s | 5 m/s | 23 d (04/23 -> 05/16/2050) |
| C3 | Ariel 3.897 -> Umbriel 4.095 km/s | 8 m/s | 24 d (03/30 -> 04/23/2050) |

The incoming V-infinity direction is not published, so the control minimises over it (and over the
inbound/outbound crossing and flight time +-1 day, because the dates are whole days). **PASS** if,
for all three legs, our minimum delta-v <= published + 3 m/s, AND for at least two of the three it
is >= 0.5 x published - 3 m/s (so a primitive that finds everything free does not pass). Anything
else is FAIL and is reported as such.

**Known limitation of the existing machinery, stated before use.** The Campagnola-Russell phase-free
VILM functions in `cyclerfinder.search.vilm` (Gamma, the Eq. 13 quadrature, `vilm_dv_min`, the
Gamma-floor flag used by `leveraging_leg` and `leveraging_chain`) are not used here. The exterior
Gamma of Eq. 25 has a pole near adimensional V-infinity 0.41 and a zero near 1.38, and every
V-infinity in this problem sits in or across that band (Ariel 4.325/5.509 = 0.79, Oberon
3.448/3.151 = 1.09, cycler Ariel 0.28, cycler Oberon 0.58). Those formulas were derived and
validated for small adimensional V-infinity (Europa endgame, 0.06-0.13). Only plain Keplerian
arithmetic and Lambert arcs are used.
