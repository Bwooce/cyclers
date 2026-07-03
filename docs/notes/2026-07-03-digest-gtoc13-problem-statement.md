# Digest: GTOC 13 Problem Statement -- "Humanity's First Robotic Exploration of a Hypothetical Exoplanetary System" (#526)

**Source:** Gregory Whiffen, Mark Wallace, Damon Landau, Etienne Pellegrini,
Gregory Lantoine, Anastassios Petropoulos, Sungmoon Choi, Brian Anderson,
Zubin Olikara, Jon Sims (JPL/Caltech), "Problem Description of the 13th
Global Trajectory Optimisation Competition -- GTOC13: Humanity's First
Robotic Exploration of a Hypothetical Exoplanetary System," released
20 Oct 2025, competition ran 24 Oct 2025 17:00 UTC -- 17 Nov 2025 18:00 UTC.
Filed: `papers/jpl-gtoc13-team-2025-problem-description-altaira-system-gtoc-jpl-net.pdf`
(text-layer, 16 pp, downloaded from
`https://gtoc.jpl.net/static/gtoc13/gtoc13_problem_statement.pdf`).
Digested 2026-07-03.

## Headline finding: the #526 proposal text mis-described the problem

Task #526 (as proposed 2026-07-02) described GTOC 13 as a "ballistic Jovian
gravity-assist tour competition." **That is not what GTOC 13 is.** The
actual problem is a ballistic (propellant-less) gravity-assist + optional
ideal-solar-sail tour of a **fictional single-star, 10-major-planet
exoplanetary system** ("Altaira" / "the Altaira system"), plus a dwarf
planet, 257 main-belt asteroids and 42 comets -- there is no Jupiter, no
Jovian moons, and no moon system at all in this problem. It is the correct
GENERAL problem CLASS this project works in (ballistic multi-body
gravity-assist tour design against a maximize-science-return objective),
just not the specific Jovian-moon-tour instance the proposal assumed. This
is reported here per the task's own instruction to report honestly rather
than force a match to the original one-line proposal.

## 1. The exo-system (Sec. 2, p. 3)

Central body Altaira: G1v main-sequence star, ~1.05 solar masses, GM =
139348062043.343 km^3/s^2 (Table 4, p. 11). 10 major planets (Keplerian,
increasing orbital period): Vulcan (hot Jupiter, defines the ecliptic),
Yavin (inner habitable-zone edge), Eden (Earth-sized, mid habitable zone),
Hoth (Venus-sized, inclined, just inside the asteroid belt), Yandi (dwarf
planet embedded in the belt, massless), Beyonce (ringed Saturn-analogue,
defines belt resonances), Bespin (super-Jovian), Jotunn (ice giant),
Wakonyingo (ice giant stripped to a super-Earth core), Rogue1 (captured
retrograde Jovian, 2:1 resonance with PlanetX), PlanetX (highly eccentric
/ inclined, 1:2 resonance with Rogue1). Plus 257 main-belt asteroids
(between Hoth and Beyonce) and 42 comets, all massless/Keplerian. All
ephemerides (semi-major axis, eccentricity, inclination, longitude of
ascending node, argument of periapsis, mean anomaly at t=0) are supplied
in CSV files (`gtoc13_planets.csv` / `_asteroids.csv` / `_comets.csv`,
not mirrored here -- not needed for this digest's purpose).

## 2. Objective function (Sec. 3, pp. 4-7)

`J = b * c * sum_k w_k * sum_i S(r_hat_{k,i}) * F(Vinf_{k,i})` -- a
maximize-science-return score over ALL scientific flybys (up to 13 per
body), summed over bodies `k` with published per-body weights `w_k`
(Table 1, p. 5: Vulcan 0.1 ... PlanetX 50, asteroids 1, comets 3):

* `b` -- grand-tour bonus (1.2 if every planet + Yandi + >=13
  asteroids/comets are flown, else 1.0).
* `c` -- time bonus, flat 1.13 for the first 7 days of the 4-week
  submission window then linearly decaying to ~1.025 (Eq., Fig. 2, p. 5)
  -- rewards early submission.
* `S(r_hat)` -- seasonal-diversity penalty (Eq., p. 6): repeated flybys of
  the SAME body at similar heliocentric phase angle are penalized
  (Gaussian suppression vs. angular separation from prior flybys of that
  body); flybys spread across viewing geometries are unaffected.
* `F(Vinf)` -- flyby-velocity penalty (Eq., p. 7): penalizes both very
  fast flybys (short observation time) and very slow/rendezvous-like ones
  (radiation-dose proxy), a non-monotonic logistic-shaped function of
  V-infinity peaking near ~2 km/s.

## 3. Dynamics (Sec. 5 + Appendix I, pp. 10, 15-16)

All bodies: pure two-body Keplerian motion about Altaira. Spacecraft:
Keplerian unless the solar sail is deployed (optional, any time interval,
freely alternating with ballistic Keplerian coast). **Ideal (perfectly
reflecting) solar sail**, area 15,000 m^2, spacecraft mass 500 kg, flux
5.4026e-6 N/m^2 at 1 AU -- `a_sail = -(2*C*A/m) * (r0/r)^2 * (u_n.u_r)^2 *
u_n`, cone angle in [0, 90] deg (the sail cannot push toward the star).
Planetary flybys: standard patched-conic gravity assist (equal-magnitude
V-infinity turn, altitude in [0.1, 100] planet radii) -- the SAME
patched-conic formalism `sin(delta/2) = mu_P / (r_p*Vinf^2 + mu_P)` this
project already uses (cf. `cyclerfinder.core` flyby routines); Yandi,
asteroids and comets are massless (V-infinity continuous, no turn).
Initial spacecraft state: interstellar-arrival asymptote at x = -200 AU
(free y, z, Vx; Vy = Vz = 0), free start epoch `t0` in [0, 200] years;
whole mission window 200 years. Tolerances (conic position 100 m /
0.1 mm/s, flyby altitude 100 m, etc.) are given in Sec. 7 for anyone
wanting to reproduce a submission-grade check.

## 4. Competition outcome (results/teams pages, `gtoc.jpl.net/gtoc13/`,
   fetched 2026-07-03; also `sophia.estec.esa.int/gtoc_portal/?page_id=1360`)

* 101 registered teams (full roster on the Teams page); **winner: THU-LAD**
  (Tsinghua University School of Aerospace Engineering, Laboratory of
  Astrodynamics -- Jialong Song, Yuming Tao, Rundao Li, Yi Zhou, Yiyang
  Qin, Yilin Zou, Jintang Li, Nan Zhang, Fanghua Jiang, Hexi Baoyin). The
  ESA GTOC-portal page notes the top-3 finishers matched GTOC 11's
  (Dyson-sphere edition) top-3 exactly.
* **Team 11, "NUAA & Friends,"** includes **Guoliang Liang and Hongwei
  Yang** -- the two lead authors of the Callisto-Ganymede-Europa
  Triple Cyclers paper this task's Part 2 covers (see the companion note
  `2026-07-03-alternating-double-cycler-operator.md`) -- a direct
  connective link between this task's two halves, though GTOC 13's
  problem itself has no moon-tour content for their double-cycler
  machinery to apply to.
* Results page publishes a "Top 13 (J>200)" leaderboard with interactive
  trajectory plots and wheel-plots per team (THU-LAD, ACT&Friends, The
  Antipodes, NUDT+XSCC, BIT&Friends, NUAA&Friends, Gravity Leap, OptimiCS,
  The True Anomalies, The Yellowhammers, VidTeam, Boreal Orbitician, DSJ)
  but **no numeric trajectory data, no methods writeups, and no papers are
  published on the site** -- only the visual leaderboard plots.

## 5. Methods papers: genuinely do not exist yet (honest negative)

Searched: WebSearch (multiple queries), arXiv, the GTOC 13 site's
Announcements/Results/Discussion pages, and the ESA GTOC Portal mirror
page. **No team methods paper, preprint, or conference writeup for
GTOC 13 was found.** The ESA GTOC-portal page explicitly lists a
"Workshop: TBA" (a GTOC 13 workshop, where teams traditionally present
methods papers that are later compiled into an Acta Astronautica special
issue, per the pattern of GTOC 9-12) with no date fixed as of this digest.
Given the competition only closed 17 Nov 2025 (< 8 months before this
2026-07-03 digest) and GTOC methods papers historically appear ~6-18
months post-competition (workshop -> journal special issue pipeline),
**this is an expected, not anomalous, gap** -- not a search failure. This
digest is therefore, as anticipated in the #526 task text, a
problem-statement-only digest: the problem itself is fully captured above;
no methods content exists anywhere to acquire.

**Re-check trigger:** re-search once a GTOC 13 workshop is announced/held,
or check `gtoc.jpl.net/gtoc13/announcements/` and arXiv (`GTOC13` /
"Altaira system") periodically; also watch for an Acta Astronautica GTOC
13 special issue (the GTOC 9/10/11/12 pattern).

## 6. Relevance to this project

* Confirms the competition's ballistic-flyby + patched-conic + optional
  ideal-sail dynamics model is essentially a superset of primitives this
  project already has (Keplerian propagation, patched-conic flybys,
  V-infinity turn-angle formula) plus an ideal solar sail this project
  does not yet model -- not immediately actionable without a live sail
  capability, and out of scope for #526 (a digest task, not a build task).
* The objective function's seasonal-diversity (`S`) and flyby-velocity
  (`F`) penalty shapes are a reusable REFERENCE for how a competitive
  multi-body tour scoring function balances repeat-visit value against
  V-infinity/observation-time tradeoffs, should this project ever build
  its own multi-objective tour-scoring genome operator -- noted here for
  future reference, not implemented.
* No cycler-specific content: GTOC 13 does not involve periodic/cycler
  trajectories (it is a single one-shot 200-year tour maximizing a score,
  not a repeating structure), so there is no catalogue or genome-operator
  angle to pursue from the problem statement itself.
