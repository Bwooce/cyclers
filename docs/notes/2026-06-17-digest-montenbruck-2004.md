# Montenbruck-Markgraf 2004 — Global Positioning System Sensor with Instantaneous-Impact-Point Prediction for Sounding Rockets

**Deep-read verdict note, 2026-06-17 AET.** Supersedes the "OFF-TOPIC"
entry in `2026-06-17-mars-cycler-wave-digest.md`. Off-topic confirmed.

## Header

- **Title:** *Global Positioning System Sensor with Instantaneous-Impact-Point Prediction for Sounding Rockets*
- **Authors:** Oliver Montenbruck, Markus Markgraf (DLR, German Aerospace
  Center, 82234 Wessling, Germany)
- **Venue:** Journal of Spacecraft and Rockets, **Vol. 41 No. 4,
  July-August 2004, pp. 644-650**
- **DOI:** 10.2514/1.1962
- **Pages:** 7

## What the paper actually is

A hardware-and-software description of the **Orion-HD GPS receiver** for
sounding-rocket range-safety applications. The receiver tracks GPS carrier
phase through boost-phase high dynamics (12 g acceleration, 29 g/s jerk
peak) and computes a real-time Instantaneous Impact Point (IIP) prediction
using a perturbed-parabolic trajectory model with first-order corrections
for Earth curvature, gravity variation, and Earth rotation (Eqs. 1-8).

Demonstrated on:
- Maxus-5 (Esrange Kiruna, 2003-04-01, 11.4-ton vehicle, 701 km apogee,
  Castor 4B motor, 800 kg payload, 12-min μg phase): IIP error < 0.5 km
  for free-flight, 3 km vs actual landing.
- VS-30/Cuma (Alcântara, 180 km apogee, 130 km range).

Hardware: GP2015 + GP2021 + ARM60B 32-bit microprocessor, 95×50 mm board.
Antenna system: helical tip antenna (boost) + single-patch can antenna
(parachute) + dual-patch (descent/reentry).

**This is a range-safety / GPS-receiver-engineering paper.** Zero
orbital-mechanics content. Zero cycler content. Zero interplanetary
content. The trajectory model is suborbital ballistic with atmosphere.

## Catalogue / KNOWN_CORPUS relevance

**None whatsoever.**

The paper was filed in `cyclers_pdf/papers/` because it appeared in the
**same JSR 41(4) July-August 2004 issue** as McConaghy-Longuski-Byrnes
"Analysis of a Class of Earth-Mars Cycler Trajectories" (DOI
10.2514/1.11939) — the issue-level metadata caught the paper in a
wide-net acquisition. The two papers share a venue and an issue but
nothing else.

No KNOWN_CORPUS anchor. No catalogue impact. No row, no reference, no
cross-link.

## Errata vs the pre-read survey

The pre-read survey filed this as "***Off-topic for cyclers** — GPS
receiver paper for sounding rockets that happened to be in the same JSR
issue as the McConaghy 2004 cycler paper. Archived for completeness but
no cycler relevance.*" The deep-read **confirms** this verdict completely
and without qualification.

## Action items for parent

1. **No KNOWN_CORPUS change.** No catalogue change.
2. **Archive-and-skip:** the paper is correctly filed in `cyclers_pdf`
   for venue-completeness (JSR 41(4) coverage) but should be flagged in
   any future acquisition-script logic as "off-topic — venue-coupled
   acquisition only" so it doesn't trigger a re-triage cycle.
3. **Optional housekeeping:** if `cyclers_pdf` ever gets an
   off-topic-flag system, this paper is a canonical example of "archived
   for venue completeness, not cycler relevance".
