# Digest - Barnes and do Vale Pereira (2025), "Preliminary Proof of the Feasibility of a Novel Mission Concept and Spacecraft Trajectory for Exploring Uranus with Small Satellites" (Aerospace 12(12):1069)

**Digested:** 2026-10-03 (publisher EPUB, single XHTML chapter, read in full as text; figures
and tables are images or table markup and were not individually inspected beyond their captions).
**Citation:** Dylan Barnes (Georgia Institute of Technology) and Paula do Vale Pereira
(University of Central Florida), Aerospace 2025, 12(12), 1069, DOI 10.3390/aerospace12121069.
Received 30 Sep 2025, revised 15 Nov 2025, accepted 17 Nov 2025, published 30 Nov 2025. Filed in
the private paper corpus as
`barnes-do-vale-pereira-2025-feasibility-novel-mission-concept-trajectory-uranus-small-satellites-aerospace-12-1069-doi-10.3390-aerospace12121069.epub`
(md5 c9038933474c2ae973474ae4af170786).

## What the concept is
A Pre-Phase A mission concept, framed as a lower-cost pathfinder that could support the
decadal-survey Uranus Orbiter and Probe (UOP): one carrier spacecraft delivers 16 CubeSats (27U,
four groups of four: A magnetosphere and particles, B multispectral remote composition, C
gravity by a GRACE-style K-band link, D atmospheric plunges) to Uranus. Printed masses: 4500 kg
launch wet mass, carrier 3848 kg wet, CubeSat constellation 640 kg combined. The carrier does the
interplanetary burns and acts as the communications relay to Earth. CubeSat power is proposed as
thermoradiative cells fed by two general purpose heat sources each (about 33.86 W electrical,
printed). The paper's stated purpose is feasibility of the technical budgets (power, radiation,
thermal, communications link with 7.2 dB margin against a 6 dB requirement, data, pointing,
fuel with 35 percent margin, mass margin 18.1 percent, volume margin 34.9 percent).

## The trajectory, with the numbers printed
Method (Section 2): MATLAB script using patched conics and Lambert's problem with JPL Horizons
ephemerides, scanning launch, Jupiter flyby and Uranus arrival dates (Figure 4 caption: launch
March 2028 to June 2034, flyby January 2034 to December 2036, arrival June 2035 to May 2040, one-day
resolution), maximising dry mass.
- Sequence: Earth launch, one powered Jupiter gravity assist, Uranus capture. Launch 5 May 2033
  on SLS Block 2 ("using only the excess energy of the launch vehicle"); Jupiter flyby 1 October
  2034 at about 2.23 million km altitude (about 32 Jovian radii), powered, delta-v 1.00 km/s;
  heliocentric speed increases "from 12.74 km/s to 25.49 km/s"; Uranus arrival 18 March 2039.
  Transfer time six years.
- Launch vehicle statements: Falcon Heavy Expendable would give a wet mass of around 1600 kg on
  this trajectory (below the 3000 kg wanted); SLS Block 2 gives 4500 kg; Starship projected higher.
  Spacecraft thruster assumed: Merlin 1D Vacuum.
- Arrival and Uranus-system phase: insertion burn delta-v 0.0067 km/s into a highly eccentric,
  highly inclined capture orbit with eccentricity 0.8, semi-major axis 147,794.5 km (5.83
  Uranian radii), periapsis altitude 4000 km (0.158 radii), apoapsis altitude 240,000 km (9.46
  radii). Periapsis speed 20.0 km/s (Section 3.3.8). Inclination value, orbital period, orbit
  lifetime and any later orbit change: not stated. The Uranus-system phase is this single capture
  orbit with CubeSat deployment, observation and decommissioning (the paper's phases 2 to 4).
- The UOP baseline quoted for comparison: 7200 kg, early 2030s launch, 13-year Earth-Earth-
  Jupiter-Uranus transfer, arrival around 2044.

## Do any Uranian moon flyby, moon tour, resonance, periodic or cycler trajectory appear
No. Text search on the full chapter (case-insensitive): Titania 0, Oberon 0, Umbriel 0, Ariel 0,
Miranda 0, Puck 0, cycler 0, periodic 0, tour 0, Hohmann 0; "resonan" 1 (the instrument name
"Magnetic Resonance Spectrometers" in Table A1, unrelated to orbits); "Uranian" 4 (behaviour,
system, radii, moons); "moon" or "moons" 11 occurrences, all in science-objective,
traceability-matrix or instrument-selection text (for example "measure the composition and
structure of the Uranian moons and rings", "magnetic sounding and gravitational investigations"
of the larger moons, "catalogue small moons"); no named Uranian moon anywhere. "flyby" 14
occurrences (counts include figure captions that appear twice in the markup): all are the single
powered Jupiter gravity assist (Methods, 3.1, ConOps phase 1, the Figure 4 and 5 captions,
the carrier thruster section), the UOP baseline's "multiple flyby maneuvers" at Earth and Jupiter,
or the generic "flyby" date axis of the optimisation; none is a flyby of a Uranian moon.
The paper does not model the Uranian moons' orbits or any encounter with them, and does not
compare the capture orbit's apoapsis to any moon orbit radius. Moon science is a stated
objective only.

## What it does NOT contain
- No Uranian moon encounter, moon tour, resonance, mean-motion-resonance analysis, periodic orbit,
  invariant manifold, quasi-periodic orbit, or cycler.
- No multi-body dynamics beyond patched conics; no Uranus-system gravity model, no orbit
  determination of the capture orbit after insertion, no inclination or period.
- No trajectory optimisation beyond the date grid maximising dry mass.
- No instrument design; Table A1 is a truncated science traceability matrix.
