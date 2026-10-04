# #890 adversarial review: the Titania-Oberon periodic orbit of the four-body model

Reviewer's note on `docs/notes/2026-10-04-890-titania-oberon-candidate.md` and the `#890` entry of
`data/OUTSTANDING.md`. Reviewed at commit `a6fac4d6`. Every number below is marked COMPUTED (by
the reviewer, with the reviewer's own code, not the build's diagnostics), READ (from the note,
the code, the outputs or git) or INFER (reasoning, no computation). The reviewer's scripts are
scratch files and are NOT in the repository; nothing here is under test. Section 12 says what
that means for how far these numbers may be relied on.

The reviewer's model code was written from a fresh derivation in a non-rotating frame, in
kilometres and seconds, and imports nothing from the repository except where a comparison with
the repository's code is the point of the check.

## 1. Verdicts

| # | Claim | Verdict | Numbers |
|---|---|---|---|
| 1 | `core/ccr4bp.py` (through the planar right-hand side used by the build) implements the planar concentric circular restricted four-body problem, with the perturber's sense, rate and indirect term right | ESTABLISHED | COMPUTED: the repository's right-hand side and the reviewer's independently derived equations agree to 3.1e-15 (relative) over 2,000 random states and times. No defect of the `#891` kind. |
| 2 | The model's constants are the physical ones | ESTABLISHED WITH CORRECTION | COMPUTED from URA111: the model's periods are 8.706399 d (Titania) and 13.465984 d (Oberon); the fitted real ones are 8.705868 d and 13.463237 d. The model's synodic period is 24.63245 d, the real one 24.63740 d. See section 2. |
| 3 | The refined state is a periodic orbit of that model, period five synodic periods (123.162 d) | ESTABLISHED | COMPUTED: closure after one cycle 28 m and 0.20 cm/s (DOP853, reviewer's code); 0.40 m (Radau). Reviewer's monodromy has no multiplier at 1 and det(M - I) = 2.9e5, so a true fixed point lies within the numerical noise of the stated state. |
| 4 | It flies by Titania at 1,977 km and Oberon at 1,364 km altitude each cycle | ESTABLISHED WITH CORRECTION | COMPUTED: 1,976.9 km and 1,364.2 km in the model as built. The altitudes depend on the constants: with the URA111 mean motions the re-converged orbit has 1,794 to 1,821 km and 1,254 to 1,266 km (section 2). Quote them as "about 1,800 to 2,000 km and 1,250 to 1,350 km". |
| 5 | No other close approach to either moon | ESTABLISHED | COMPUTED by root-finding on the range rate (16 local minima of the distance to Titania, 12 to Oberon, over one cycle): the next nearest are 122,320 km from Titania (11.9 Hill radii) and 144,464 km from Oberon (10.9 Hill radii, at the conjunction itself). Distance from Uranus 437,600 to 585,826 km. |
| 6 | The flybys are hyperbolic and turn-feasible (V4, V5) | ESTABLISHED, but V5 is not a test | COMPUTED: osculating eccentricity 1.831 (Titania) and 1.741 (Oberon) at periapsis. V5 restates V4(a) in every form (section 3). |
| 7 | The orbit "continues the patched-conic closure without a fold" | ESTABLISHED WITH CORRECTION | READ (`continuation.jsonl`): 13 ladder points from 0.01 to 1 of the physical masses, all converged, determinant sign constant. Not recomputed. It is a continuation from one percent of the masses, not from the patched conic, and the object changes character on the way (section 3). |
| 8 | Multipliers 8.4e5, 1.79, 0.559, 1.2e-6 | ESTABLISHED | COMPUTED with the reviewer's own variational equations: 8.381e5, 1.7857, 0.5600, 1.19e-6. Second pair real. Section 5. |
| 9 | The pre-registration was written before the results and not edited afterwards | ESTABLISHED | COMPUTED with git: lines 1 to 125 of the note at HEAD are byte-identical to commit `0e5d63ae` (15:49:41); the first output file is stamped 15:50:59. Section 4. |
| 10 | The pre-registered closure tests V2 and V3 pass | NOT ESTABLISHED as registered; the orbit's existence does not depend on them | READ and COMPUTED: 1.64 and 1.59 cm/s against a 1 cm/s limit. The limit was unattainable by construction for this orbit (section 4). |
| 11 | "It is not a trajectory in the real system without correction" | ESTABLISHED, and it says less than it seems to | READ. The uncorrected run measures sensitivity only. COMPUTED by the reviewer (scratch code, section 6): corrected, manoeuvre-free arcs exist in a URA111 force model at two start epochs, one of three cycles (7 flybys) and one of six cycles (13 flybys, 739 days), all flybys between 1,100 and 1,830 km altitude. |
| 12 | The object is a cycler | NOT ESTABLISHED | It is a periodic orbit of an idealised model. Class and validation level in section 7. |

## 2. Is the model right? (question 1)

**Equations.** COMPUTED. Reviewer's derivation: Uranus and Titania on circles about their
barycentre B at rate n_T = sqrt((GM_U + GM_T)/a_T^3); Oberon on a circle of radius a_O about B at
rate n_O; B accelerates as a point at B attracted by Oberon, which gives the indirect term
-GM_O R_O/a_O^3 in the spacecraft's equation. Transformed to the rotating frame this is exactly
what `planar_rhs_batch` integrates (agreement 3.1e-15). The perturber's synodic rate is signed,
-0.35345 in the frame's units, so an outer moon regresses; correct. The derived parameters
(mu 3.915883e-5, perturber mass 3.543106e-5, radius ratio 1.3374139, rate) agree with the
reviewer's to the last digit.

**"Concentric".** INFER, with one number COMPUTED. Uranus is 17 km from B, so "Oberon circles
B" and "Oberon circles Uranus" differ by 17 km in 583,511. The indirect term treats the
Uranus-Titania pair as a point at B; the error is of relative order (GM_T/GM_U)(a_T/a_O)^2 =
2e-5 of a term that is itself 3.5e-5 of the central attraction. Acceptable.

**Mass bookkeeping.** READ and COMPUTED. `TwoMoonModel` takes the planet as the registry system
GM (5,794,556.4) minus Titania and Oberon, so Miranda, Ariel and Umbriel (172.9 km^3/s^2) sit in
the planet. COMPUTED: moving them out of the central mass changes the flyby altitudes by 0.06 km
and 0.05 km, so the lumping is harmless at the model level. One inconsistency: Oberon's rate uses
(GM_U + GM_O)/a_O^3 and so leaves out Titania's mass, although Titania is interior to Oberon and
acts on it as central mass; Titania's rate, correctly, leaves out the exterior Oberon. This is
the inherited convention of `two_body_synodic_rate` (chosen so that one limit reduces exactly to
a three-body problem), not an error of sense. It changes n_O by 2.0e-5. COMPUTED: with Titania's
mass counted the re-converged orbit has altitudes 1,952 km and 1,349 km (25 and 15 km lower).

**Constants.** READ: GM 226.9 and 205.3 km^3/s^2, radii 788.9 and 761.4 km, as in
`core/satellites.py` (cited there to JPL SSD and URA111). The reviewer did not check them against
a paper in this session. COMPUTED from URA111 (fit over 2030 to 2050): mean distances 436,282 km
and 583,451 km (registry 436,298 and 583,511), periods 8.705868 d and 13.463237 d. The registry
distances with Kepler's law give periods 6e-5 (Titania) and 2e-4 (Oberon) too long, so the
model's synodic period is 24.63245 d against the real 24.63740 d, and the cycle is 123.162 d
against 123.187 d.

**How much the result depends on these choices.** COMPUTED, by re-converging the orbit in the
reviewer's model with a half-cycle symmetric Newton iteration and continuation in 20 steps:

| Change to the model | Titania altitude, km | Oberon altitude, km |
|---|---|---|
| none (the build's model) | 1,976.9 | 1,364.2 |
| Oberon's rate with Titania's mass counted | 1,951.7 | 1,349.0 |
| URA111 mean motions, registry distances | 1,794.4 | 1,263.6 |
| URA111 mean motions, distances from Kepler's law | 1,810.0 | 1,254.2 |
| Oberon's GM plus 2 percent | 1,978.6 | 1,388.5 |
| inner moons' GM taken out of the central mass | 1,976.9 | 1,364.2 |

(Titania's GM plus 2 percent was attempted and the reviewer's continuation failed at the first
step for a bookkeeping reason in the reviewer's script; not a finding.) The orbit persists under
every change tried, on the same sides of the moons. The four-figure altitudes in the note are a
property of the registry's distances, not of the Uranian system: they move by 180 km and 100 km
when the mean motions are the real ones.

**A defect found on the way, in a shared module** (reported to the coordinator at once).
`src/cyclerfinder/data/validation/v4_uranus.py` pairs `URANUS_J2 = 3.34343e-3` with
`URANUS_R_EQ_KM = 25559.0`, and the `realeph` stage of the `#890` script copies the pair.
COMPUTED: the reviewer propagated each of the five moons as a massive body from its URA111 state
under Uranus, J2 about the IAU pole and the other four moons, and compared with URA111.

| J2 and reference radius | J2 R^2, km^2 | Position error after 30 / 123 / 370 d, km: Miranda | Ariel | Umbriel | Titania | Oberon |
|---|---|---|---|---|---|---|
| 3510.68e-6, 25,559 km | 2.2934e6 | 4 / 17 / 51 | 0.4 / 1.7 / 5 | 0.6 / 2.0 / 6 | 1.8 / 8 / 24 | 3.7 / 16 / 47 |
| 3.34343e-3, 26,200 km | 2.2950e6 | 0.9 / 3.8 / 11 | 1.6 / 6 / 19 | 0.3 / 1.6 / 4 | 1.5 / 7 / 21 | 3.6 / 15 / 46 |
| 3.34343e-3, 25,559 km (the repository's pair) | 2.1841e6 | 334 / 1,376 / 4,142 | 130 / 529 / 1,583 | 55 / 233 / 696 | 18 / 75 / 224 | 11 / 48 / 144 |
| none | 0 | 6,920 / 28,485 / 84,364 | 2,716 / 11,072 / 33,077 | 1,152 / 4,843 / 14,492 | 340 / 1,408 / 4,220 | 157 / 690 / 2,084 |

URA111 follows a J2 R^2 product of about 2.29e6 km^2 and rejects the repository's 2.18e6 by a
factor of 3 to 300 depending on the moon. The repository's J2 is 4.8 percent too small for the
radius it is paired with. Which paper prints which number the reviewer has not read in this
session (recollection: 3.34343e-3 is an older ring-occultation value referred to 26,200 km, and
3510.68e-6 at 25,559 km is the later satellite-orbit value); the computation shows only which
product the ephemeris obeys. It does not affect the `#890` model claim (the model has no J2).

## 3. Is the orbit what it is said to be? (question 2)

**Flybys.** COMPUTED from the reviewer's integration, with the periapses found by a root-finder
on the range rate.

| | Titania | Oberon |
|---|---|---|
| time in cycle | 0 | 61.5811 d |
| periapsis distance, altitude | 2,765.8 km, 1,976.9 km | 2,125.6 km, 1,364.2 km |
| in Hill radii | 0.27 | 0.16 |
| side | outside (away from Uranus) | inside (Uranus side) |
| speed relative to the moon at periapsis | 0.4819 km/s | 0.5145 km/s |
| osculating eccentricity, V-infinity, turn at periapsis | 1.831, 0.2611 km/s, 66.2 deg | 1.741, 0.2675 km/s, 70.1 deg |
| osculating V-infinity at the Laplace sphere, at 1 and at 2 Hill radii | 0.2511, 0.2432, 0.1902 km/s | 0.2620, 0.2581, 0.2297 km/s |
| change of velocity direction between entering and leaving the Laplace sphere, 1 Hill radius, 2 Hill radii | 55.7, 59.6, 60.2 deg | 66.0, 68.9, 76.8 deg |
| time inside the Laplace sphere | 0.45 d | 0.61 d |
| bend available at the 50 km floor, with the osculating V-infinity at periapsis | 106.0 deg | 102.4 deg |
| bend available at the surface | 107.9 deg | 104.4 deg |

Both flybys are above the surface by more than two moon radii and above the registry's 50 km
floor by a factor of 27 to 40. Inside the Laplace sphere they are two-body hyperbolae to 2 to 4
percent in V-infinity. Outside it they are not: the osculating V-infinity falls by 27 percent
(Titania) and 14 percent (Oberon) between periapsis and two Hill radii, so "the V-infinity" and
"the asymptotes" in the patched-conic sense do not exist for this orbit. V-infinity is about three
times the Hill velocity (0.086 km/s at Titania), which is the regime where that is expected.

**The demanded-turn gate (V5) is not evidence here.** INFER, supported by the numbers above. The
gate asks whether a moon can supply a turn at or above an altitude floor. Applied to a
continuously integrated trajectory of a model in which the moon's gravity does the turning, the
answer is the periapsis altitude and nothing else. The note concedes this for the osculating
version and counts the Laplace-sphere version as a pass; it is the same statement. In addition,
the "V-infinity" fed to that version (0.351 and 0.333 km/s) is the speed at the sphere, not a
V-infinity; the osculating value there is 0.251 and 0.262 km/s. With the speed, the available
bend comes out at 86.8 and 88.1 deg (the note's figures, reproduced); with the osculating
V-infinity it is 108.4 and 103.7 deg (the coordinator's 103.8). The ratio "0.64" or "0.75" is
therefore a number without a fixed meaning. What protects against the withdrawn rows' error is
not this gate but the fact that the trajectory is one continuous integration through both flybys.

**The orbit no longer crosses the moons' orbits.** COMPUTED. The Uranus-centred osculating conic
at mid-leg has semi-major axis 510,963 km, eccentricity 0.14172, periapsis 438,548 km and
apoapsis 583,378 km, period 11.034 d. Titania's orbit is at 436,298 km and Oberon's at
583,511 km. The patched conic it was continued from had periapsis 435,210 km and apoapsis
587,001 km (READ) and crossed both. So at the physical masses the transfer ellipse lies between
the two orbits: it stays 2,250 km outside Titania's and stops 133 km short of Oberon's, and the
encounters happen because the moons' gravity reaches out to it (2,766 and 2,126 km), not because
the conic intersects a moon's path. Each flyby reverses the radial velocity and so rotates the
line of apsides: COMPUTED, the argument of periapsis advances by 26.8 deg across the Oberon
flyby, and by the remainder of the 52.6 deg per cycle at Titania. A strict zero-radius patched
conic with this ellipse has no encounter at all. This does not make the flybys unreal, but it
does mean the sentence "the orbit continues the patched-conic closure" describes the route by
which it was found, not what it is. The periapsis distances are about half the patched-conic
hyperbola's (READ: 5,351 km at Titania against 2,766 km).

**Margin a real mission would want.** INFER. Altitude is not the binding margin; delivery
accuracy is. READ (`sensitivity.json`): 1 m/s along-track at the preceding apoapsis moves the
Oberon periapsis by 1,243 km, which is 91 percent of the altitude. So the margin against impact
is about 1 m/s of uncorrected velocity error eleven days out, or about 1 mm/s per 1.2 km.
The altitudes themselves are comfortable: in the reviewer's real-ephemeris arcs (section 6) the
lowest flyby is 1,105 km, so a delivery error of 100 km at three sigma, which needs the
apoapsis velocity controlled to a few centimetres per second, leaves a factor of ten in hand.

**Continuous laps.** COMPUTED. Integrated without correction for three cycles, the return error
is 28 m (DOP853) or 0.40 m (Radau) after one cycle, 15,600 km or 332 km after two, and the orbit
is lost in the third. With a multiplier of 8.4e5 per cycle no double-precision integration can
follow this orbit for two cycles. Any statement about "three laps" has to be made by multiple
shooting, not by one propagation.

## 4. Is the pre-registration honest? (question 3)

COMPUTED with git. `git diff 0e5d63ae HEAD` on the note shows one hunk, starting at line 123,
all additions; lines 1 to 125 at HEAD are identical to the pre-registration commit. That commit
is stamped 15:49:41 and contains the note, the module and the tests; the first output file
(`rebuild.json`) is stamped 15:50:59 and the continuation output 16:24:50. The reviewer cannot
exclude that results were seen before the commit and not saved, and has no reason to think so:
the registered threshold failed, which is not what a threshold fitted to a known answer does.

**Were V2 and V3 legitimate tests?** INFER, with COMPUTED support. They were ill-posed. A
corrector that stops at a residual r leaves the start state about r from the true fixed point,
and one cycle multiplies that by up to the leading multiplier. With r = 2.6e-11 and a multiplier
of 8.4e5 the expected return error is up to 2.2e-5 in the model's units, 9.5 km; the registered
limit was 1 km and 1 cm/s; the observed 0.24 km and 1.6 cm/s is well inside what the instability
allows and proves nothing against the orbit. The fallback clause recognised this but was keyed to
a multiplier above 1e6, an arbitrary line the orbit missed by 16 percent.

**Is the post hoc Newton step legitimate?** INFER and COMPUTED. Yes as numerical polishing: it
moves the start by 0.08 mm (READ) and the reviewer's independent code then closes to 28 m and
0.4 m with two integrators, with a monodromy for which det(M - I) is far from zero. No, as a way
of saying the pre-registered test passed. The note says this itself, plainly, and does not move
the goalposts in its own verdict. The ledger entry is equally plain.

**What the criteria should have been.** (a) Existence: the half-cycle symmetry conditions at the
Oberon crossing, |y| and |v_x| below a bound, under two integrators. Those are what define a
symmetric periodic orbit, the half-cycle map amplifies by about 900 rather than 8.4e5, and the
orbit meets them at 0.27 m (READ). (b) Full-cycle closure bounded by the leading multiplier
times the corrector residual, for any multiplier. (c) Independence: a second code, not a second
integrator on the same right-hand side. Radau against DOP853 on one right-hand side cannot detect
a wrong model, which is how the `#891` and `#892` defects survived. The coordinator's check and
this review supply (c) after the fact.

**Other deviations.** READ. The continuation starts at 0.01 of the masses, not 0.001, for a
stated cost reason; the "five apoapsis passages" criterion was a miscount (six); the direct
solve at full mass did not converge; the real-ephemeris look ran outside its gate. All are
disclosed in section 2.1 of the note. None favours the candidate.

**What the tests cover.** READ. The nine tests check frames, the symmetry, the two-body limit,
the hyperbola construction and the shooter's Jacobian. None of them asserts anything about the
orbit: no test integrates the stored state and checks its closure or its flyby altitudes. The
result lives only in JSON files.

## 5. Isolated solution or family? (question 4)

COMPUTED. The reviewer's monodromy (own variational equations in the non-rotating frame, mapped
to the rotating frame at both ends) has eigenvalues 8.381e5, 1.7857, 0.5600, 1.19e-6 from the
full-cycle integration and 8.381e5, 1.7857, 0.5600, 1.19e-6 from the half-cycle map composed with
its mirror image (the better-conditioned route; determinant 1.0006). From the symplectic
invariants the two pair sums are 838,111 and 2.3457; the second is above 2, so the second pair
is real, 1.7857 and its reciprocal 0.5600, product 1.0000.

What this means. The model is periodically forced, so no multiplier is pinned at 1 (there is no
energy integral and no time-shift freedom). Both pairs are off the unit circle: the orbit is
hyperbolic with two unstable directions, one violent (e-folding time 9.0 d) and one mild
(e-folding time 212 d). det(M - I) = 2.9e5 is far from zero, so:

* the fixed point is isolated in the whole four-dimensional phase space. There is no family at
  fixed masses, and no nearby periodic orbit of the same period, symmetric or not. Symmetric
  shooting cannot see asymmetric solutions in general, but none can be close to this one;
* it continues uniquely under any small change of the model that keeps the forcing periodic.
  The reviewer continued it in five directions (table in section 2) without difficulty;
* there are no elliptic directions, hence no surrounding tori and no stable neighbourhood;
* no bifurcation is near at the physical masses (that needs a multiplier at +1 or -1; the mild
  pair at 1.79 is the one to watch if a parameter is varied far).

The half-cycle shooting Jacobian has singular values 7.2e4 and 0.128 (determinant -9,189): well
conditioned for Newton but with one soft direction, which is the mild pair seen at half cycle.

What the search could not see. INFER. `#888` enumerated symmetric closures only (flybys at
conjunction and opposition). A patched-conic chain Titania, Oberon, Titania with period five
synodic periods has two unknown flyby phases and two V-infinity matching conditions, so
asymmetric closures are isolated solutions too and would come in mirror pairs; none was looked
for. Legs of 3.5 and 4.5 synodic periods were outside the enumeration (READ, note section 2.8).
"The one closure that passes the gate" is one of the symmetric closures with legs up to 65 days,
not one of all closures.

## 6. What is it worth physically? (question 5)

**What the model leaves out, measured.** COMPUTED from URA111, 2030 to 2050, relative to
Uranus's centre:

| | Titania | Oberon |
|---|---|---|
| distance from Uranus, half-range | 1,109 km | 1,731 km |
| departure from uniform circular motion along the orbit over a 1.7-year window, rms and largest | 1,135 km, 2,343 km | 1,612 km, 3,252 km |
| largest distance out of Titania's mean orbital plane | 14 km (over 1.7 years), 176 km (over 20 years) | 1,659 km (over 1.7 years), 1,782 km (over 20 years) |

These are the same size as the flyby periapsis distances (2,766 and 2,126 km) and several times
the distance by which the transfer ellipse clears each orbit. Oberon in particular sits up to
1,780 km above or below the plane the planar model lives in. The real flyby geometry therefore
differs from the model's at order one from cycle to cycle. Also INFER: one cycle is 14.15
Titania periods, so the encounter longitude advances by about 53 deg per cycle and the moons'
radial and vertical offsets at the flybys never repeat; the real system has no periodic orbit
here, only (possibly) an aperiodic flyby chain. Smaller effects, estimated not integrated: J2
rotates the spacecraft's apsides by about 1e-3 rad per cycle (about 500 km at the orbit's
distance), time-independent and symmetric, so it is absorbed by re-converging; Ariel and Umbriel
as separate bodies act at a few 1e-6 of the central attraction with their own periods.

**The build's real-ephemeris run.** READ, with arithmetic. All three uncorrected runs pass
within 2,000 km of Oberon 0.5 to 0.7 d before the model's flyby time and then leave. At a
relative speed of 0.27 km/s, 0.6 d is 14,000 km of along-track phase, which is what about 1 m/s
of initial velocity difference produces over 61 days. So the run confirms the sensitivity and
nothing more; it does not bear on whether a corrected ballistic trajectory exists. The build
says "not attempted" for the correction, correctly. Its J2 is the mispaired value of section 2
and its pole is Titania's orbit normal; neither changes its outcome.

**Does a ballistic real-system analogue exist? Reviewer's computation.** COMPUTED, scratch code;
a lead, not a result. Force model: Uranus point mass, J2 (3510.68e-6 at 25,559 km) about the IAU
pole, and Miranda, Ariel, Umbriel, Titania and Oberon as point masses at their URA111 positions
with direct and indirect terms, three-dimensional. Positive control: this force model reproduces
the five moons' own URA111 motion to 0.4 to 4 km over 30 days (table in section 2), and fails
that control by one to two orders of magnitude if J2 is removed or mis-scaled. Method: multiple
shooting with a node every 1/24 of a cycle (5.13 d), every node's six components free,
continuity conditions only, minimum-norm Newton step; a homotopy in which Titania and Oberon
move from circular coplanar reference orbits (fitted to URA111 over the window) to their URA111
positions while the inner moons' separate gravity and J2 are switched on. No manoeuvre is
allowed anywhere. Each arc has two extra nodes (10 d) beyond its first and last flyby.

| | Arc A | Arc B |
|---|---|---|
| start (mean Titania-Oberon conjunction, TDB) | 2030-01-11 23:57 | 2031-06-13 09:13 |
| cycles, flybys, first to last flyby | 3, 7, 369.6 d | 6, 13, 739.2 d |
| homotopy steps | 10 of 0.1 | 0.1, 0.2, 0.3, then 11 smaller steps with a secant predictor |
| Titania flyby altitudes, km | 1,650 / 1,810 / 1,797 / 1,635 | 1,745 / 1,614 / 1,512 / 1,525 / 1,694 / 1,826 / 1,820 |
| Oberon flyby altitudes, km | 1,175 / 1,108 / 1,105 | 1,227 / 1,336 / 1,360 / 1,288 / 1,181 / 1,106 |
| osculating eccentricity at periapsis | 1.69 to 1.81 | 1.70 to 1.81 |
| osculating V-infinity, km/s | 0.265 to 0.276 | 0.265 to 0.277 |
| osculating turn, deg | 67 to 72 | 67 to 72 |
| flyby plane tilted from Titania's mean plane by up to | 8.2 deg | 8.0 deg |
| spacecraft out of that plane by up to | 3,167 km | 3,356 km |
| distance from Uranus | 435,933 to 587,539 km | 435,976 to 587,618 km |
| nearest approach to Umbriel, Ariel, Miranda | 170,507 / 245,942 / 306,400 km | 169,751 / 245,882 / 306,666 km |
| largest junction discontinuity, as converged (DOP853, 1e-11) | 0.25 m | 0.29 m |
| the same recomputed serially, DOP853 at 1e-13 and LSODA at 1e-12 | 0.13 m, 0.007 mm/s; 0.55 m, 0.011 mm/s | 0.20 m, 0.010 mm/s; 0.71 m, 0.011 mm/s |
| one continuous 61.6-day propagation through each inner flyby, end against the node (both integrators) | 0.03 to 0.27 km, at most 1.9 mm/s (5 legs) | 0.003 to 0.31 km, at most 2.2 mm/s (11 legs) |

* At the circular end of arc A the code converges from the build's orbit to flybys at 1,816 to
  1,823 km (Titania) and 1,265 to 1,266 km (Oberon), in agreement with the reviewer's separate
  planar code with URA111 mean motions (1,810 and 1,254 km; the reference distances differ
  slightly). The two reviewer codes cross-check each other there.
* The continuous 61.6-day legs differ from the nodes by up to 0.3 km. That is the junction noise
  (0.2 m) multiplied by the per-leg growth (915), which is what a true trajectory through
  slightly noisy nodes looks like; it is not a manoeuvre and not a gap.
* The method did not always work. For arc B the first attempt with steps of 0.1 and a crude
  predictor (flyby nodes carried with their moon) stalled at the fourth step with the residual
  stuck near 2,400 km; smaller steps with a secant predictor converged in 3 to 5 Newton steps
  each. A stall of this kind is a failure of the corrector, not evidence about the trajectory.
* The arcs are not the periodic orbit and not unique. With all nodes free there is a
  six-parameter family of nearby ballistic arcs, and the minimum-norm step picks one. Even at
  the circular end the six-cycle arc's Titania altitude drifts from 1,835 to 1,800 km along the
  arc. In the URA111 model the altitudes wander over 300 km (Titania) and 250 km (Oberon).

So manoeuvre-free arcs with 7 and 13 alternating flybys at safe altitudes exist, at two epochs,
in a force model that tracks URA111, close to the model orbit. That answers the question as
asked: the object is not a model object only, and a ballistic real-system analogue is more than
plausible for at least two years. What this does NOT show: that such an arc can be continued
indefinitely (the encounter geometry never repeats, and the altitudes wander); that it survives
the Sun, J4 and the moons' figures (all small); that an independent code agrees (the code is
the reviewer's, with one positive control and no tests); or anything about the first and last
flyby of each arc, which sit 10 days from a free end and are less constrained than the inner
ones. "Junction discontinuity below one metre" is the same kind of evidence as the build's own
shooting residual; a single propagation over even two cycles is impossible for this orbit
(section 3).

**Maintenance cost.** COMPUTED from two-body flyby formulae, then INFER. The leading multiplier
gives an e-folding time of 9.0 d and a factor of 915 per leg. The growth is concentrated at the
flybys: the derivative of the turn angle with periapsis distance is 0.0123 deg/km (Titania) and
0.0161 deg/km (Oberon), so each kilometre of periapsis error turns the outgoing velocity by
5.6 cm/s and 7.5 cm/s. With flybys delivered to 1 km the clean-up is of order 6 to 8 cm/s per
flyby, 0.15 m/s per cycle, 0.4 m/s per year; with 10 km delivery, ten times that. To deliver
1 km at Oberon the along-track velocity at the preceding apoapsis has to be known and controlled
to about 1 mm/s (READ: 1,243 km per m/s). These are the ordinary numbers of a moon tour. The
multiplier 8.4e5 does not by itself imply a large cost; it implies that a correction is needed
before every flyby and that nothing is forgiven if one is missed. No navigation model was run.

## 7. Is it a cycler in the catalogue's sense? (question 6)

READ: `data/README.md` (schema v4.7): `cycler` means not epoch-locked with infinite returns;
`quasi_cycler` means epoch-locked with a finite number of returns and a validity window.

* In the model the object is a `cycler` (infinite returns, no epoch), of structural kind
  non-Keplerian, with a four-body model assumption. As a statement about the Uranian system that
  class would be wrong: section 6 shows the real system has no periodic orbit here.
* In the real system the only class it could ever belong to is `quasi_cycler`, epoch-locked,
  with `n_returns` equal to the number of cycles actually demonstrated on an ephemeris, and a
  validity window. Nothing in the repository demonstrates any.

Validation level the evidence supports (spec section 14), in the model: V0 (internal
consistency) and the V1 kind of cross-check (three independent integrations now agree: the
build's, the coordinator's, the reviewer's). V2-ballistic asks for three continuous laps in the
defining model; COMPUTED, that is impossible in double precision for this orbit (section 3), so
it is not met as written and needs a ruling on whether a multiple-shooting demonstration counts.
V3: not met; the build's ephemeris run is uncorrected and loses the orbit; the reviewer's arc is
scratch work. V4: none, and the existing Uranian V4 lane cannot represent a flyby (READ,
`v4_uranus.py` lines 349 to 363: a moon's force is zeroed inside the softening radius, which is
set to the Hill radius) and carries the J2 defect. So: V1, in a model.

**The six withdrawn rows.** READ (`data/withdrawn/`). They were Uranian two-moon
"quasi-cyclers" at V-infinity near 2 km/s with 12-day legs, labelled V4 and candidate-novel. Each
leg was a valid Kepler arc; the V-infinity magnitudes matched at the encounters; the directions
were never compared, because the validation lane restarted every leg from its own Lambert
solution; the demanded turn exceeded what the moon could supply by 1.8 to 28 times. The `#890`
object is different in kind: one continuous integration through both flybys in a model where
the moon's gravity does the turning, at V-infinity 0.26 km/s, with a turn the moon supplies at
1,400 to 2,000 km altitude. The same mistake in different words would be any of these:

1. Calling the demanded-turn gate "passed" as though it were a second witness. It is the
   altitude again (section 3).
2. Writing a row whose model is a circular coplanar four-body problem and whose bodies are
   "Titania, Oberon" and letting a reader take it for a trajectory of the Uranian system. The
   withdrawn rows carried `source_ephemeris: URA111` on a circular-coplanar construction.
3. Validating a real-ephemeris version leg by leg. Any lane must integrate through every flyby
   with the moon's full gravity and must report the size of every junction discontinuity, in
   metres and millimetres per second, next to the claim.
4. Treating "it closed" as the finding. Here the pre-registered closure failed, and the object's
   existence rests on symmetry conditions, a monodromy and independent codes.
5. Quoting model constants as system facts (the altitudes, section 2).

## 8. What the mathematics makes unsurprising (question 7)

INFER; a separate literature check is running and this is not one.

* A resonant chain of flybys that closes in the patched-conic limit is a collision chain of the
  zero-mass problem, and periodic orbits that shadow such chains at small mass are the classical
  "second species". READ (ledger): the project holds a digest of a theorem of this kind for one
  small secondary in a time-independent problem. Nothing proves the two-moon, periodically forced
  case at these masses, but existence is what one expects, and a non-degenerate closure continues
  by the implicit function theorem for small enough mass. The computed content of `#890` is that
  the continuation reaches the physical masses, where the periapsis is 0.16 to 0.27 Hill radii
  and "small" is no longer obviously true.
* Because the orbit is hyperbolic it persists under any small perturbation of the model, periodic
  or not; under a non-periodic one it becomes a bounded non-periodic solution. That is why a
  real-ephemeris analogue was to be expected and why section 6 found one.
* The construction is the obvious one: a near-Hohmann ellipse between two neighbouring moons that
  grazes both and has its apsides turned by each flyby. The published patched-conic analogues
  are the moon-to-moon cyclers at Jupiter and Saturn (READ, ledger: Russell and Strange 2009). A
  claim of novelty in kind is not available. At most: "first computed for Titania and Oberon",
  and only if the literature check allows it.

## 9. Wording corrections (question 8)

In `docs/notes/2026-10-04-890-titania-oberon-candidate.md` (section 2 only; section 1 is the
pre-registration and must stay as it is):

| Where | Text | Problem | Replacement |
|---|---|---|---|
| 2.4 | "V-infinity measured at the Laplace sphere crossings is larger (0.351 km/s Titania, 0.333 km/s Oberon), and each sphere crossing takes about a day: the encounters are not two-body hyperbolae" | Those are speeds at the sphere, not V-infinity. The time inside is 0.45 d and 0.61 d. The conclusion is right for the wrong reason. | "The speed relative to the moon at the Laplace sphere is 0.351 km/s (Titania) and 0.333 km/s (Oberon); the osculating V-infinity there is 0.251 and 0.262 km/s against 0.261 and 0.268 km/s at periapsis, and falls to 0.190 and 0.230 km/s at two Hill radii. Inside the sphere (0.45 d and 0.61 d) the encounters are close to two-body hyperbolae; an asymptotic V-infinity is not defined." |
| 2.5, V5 row | "PASS ... Titania 55.7 / 86.8 deg (ratio 0.64), Oberon 66.0 / 88.1 deg (0.75)" | Not independent of V4; the "available" figures use a speed in place of V-infinity. | "Restates V4(a): the trajectory is integrated through the flyby, so the turn is supplied at the periapsis altitude. Direction change across the Laplace sphere 55.7 and 66.0 deg; bend available at 50 km with the osculating V-infinity 108 and 104 deg. Not counted as evidence." |
| 2.5, V6 row | "distance to Uranus 435,000 to 588,000 km" | Wrong; the ledger has the right figures. | "437,600 to 585,800 km" |
| 2.6 | "the orbit is only flyable with routine targeting" | "Flyable" is not shown; no targeting or navigation analysis exists. | "any use of the orbit needs a correction before every flyby; no targeting or navigation analysis has been done" |
| 2.7 | "J2 = 3.34343e-3 about Titania's orbit normal" | The J2 value is paired with the wrong radius (section 2 of this review). | add "(this J2 is 4.8 percent too small for the 25,559 km radius used; it does not change the outcome)" |
| 2.7 | "This is the expected outcome for an orbit with 1.3 km of flyby shift per mm/s" | True, and it should say what the run does not show. | add "The run measures sensitivity; it does not show that no ballistic trajectory exists in the real system." |
| 2.9 | "that continues the patched-conic closure without a fold" | The continuation started at 0.01 of the masses, no fold was seen at 13 ladder points, and the orbit's ellipse no longer crosses either moon's orbit. | "reached by continuation in the moons' mass from one percent of the physical value, with no fold detected at the 13 values tried; at full mass its Uranus-centred ellipse (438,500 to 583,400 km) lies between the two moons' orbits and the periapses are about half the patched-conic ones" |
| 2.9 | "encounters Titania at 1,977 km altitude and Oberon at 1,364 km altitude once each per half cycle" | It is once each per cycle. The altitudes are constants-dependent. | "encounters Titania at about 1,980 km and Oberon at about 1,360 km altitude once each per cycle (1,800 and 1,260 km if the real mean motions are used)" |
| 2.9 | "It is not a trajectory in the real system without correction" | Reads as a negative about the object. | "Uncorrected, it does not survive in a real-ephemeris model, which is a statement about sensitivity. Whether a corrected ballistic trajectory exists was not tested by the build." |

In `data/OUTSTANDING.md`, `#890`:

| Text | Problem | Replacement |
|---|---|---|
| "A periodic Titania-Oberon-Titania orbit with real flybys of both moons exists in the planar circular four-body model" | "real" invites the reading "in the real system". | "... with hyperbolic flybys of both moons, inside their Laplace spheres, exists in the planar circular four-body model" |
| "because the flyby offset changes the orbit's energy by a few percent and that accumulates to a full phase slip over 5.5 revolutions" | An explanation from a discarded diagnostic stated as fact. | "the build attributes this to ..." |
| "the turn measured across Oberon's sphere of influence is 66.0 degrees against 103.8 available at a 50 km floor (ratio 0.64)" | Correct numbers, but not an independent test. | add "(this restates that the periapsis, 1,364 km, is above the floor)" |
| "Titania flyby at 1,977 km altitude and Oberon flyby at 1,364 km" | Model constants quoted as facts. | add "(with the registry's distances; about 1,800 and 1,260 km with URA111 mean motions)" |
| "Real ephemeris ... every run reaches Oberon half a day early (one impacts) and none returns to Titania." | Leaves the impression that the real system rejects the orbit. | add "These runs are uncorrected and measure sensitivity only." |

## 10. The single most serious weakness

The claim as recorded is a periodic orbit of an idealised model, and everything that would make
it a statement about Uranus's moons is missing from the repository: the only real-system
evidence in the build is an uncorrected run that was bound to fail, the orbit's quoted altitudes
move by 100 to 180 km with the choice of constants, no test pins the result, and by the
catalogue's own definitions the class it has in the model (`cycler`, periodic for ever) is one
it cannot have in the real system, where the moons' radial and out-of-plane excursions are as
large as the flyby distances and never repeat. The model result itself survived every attack
made here, and the reviewer's own scratch computation suggests the real-system version exists
as an aperiodic ballistic flyby chain. That makes the gap a matter of work not yet done rather
than a defect of the object, but it is still a gap: nothing reproducible in the repository
supports any statement about the real system.

## 11. What is required before a catalogue row

It should not be a row on the present evidence. Before one is considered:

1. The real-ephemeris correction built in the repository as a module with tests, not taken from
   this review (the reviewer's arcs are a target to reproduce, not evidence): unsoftened point
   masses from URA111, J2 with a consistent radius, three dimensions, multiple shooting with
   every junction discontinuity reported, at least three cycles at each of at least three start
   epochs, each flyby also integrated through continuously with a second integrator, and the
   force model's positive control (it must reproduce the moons' own motion, and must fail to
   when J2 is removed) as a test. It should also say how long an arc can be continued before a
   flyby drops below the floor or the corrector fails, since that sets `n_returns`.
2. An external check on a different code base (GMAT is installed) through at least one full
   cycle, flyby to flyby, with the moons' gravity unsoftened. The existing `v4_uranus` lanes
   cannot serve.
3. The row, if any, as `quasi_cycler`, epoch-locked, `n_returns` equal to the cycles
   demonstrated, a validity window, and the model orbit cited as the seed, not as the object.
4. A regression test for the model orbit itself: the stored state integrated by a code that does
   not share the build's right-hand side, asserting the half-cycle symmetry conditions, the two
   periapsis altitudes and the leading multiplier.
5. The model orbit recomputed with URA111-consistent mean motions, and the altitudes quoted from
   that.
6. A ruling on V2 for orbits too unstable to propagate for three laps.
7. A maintenance budget from a stated navigation error model.
8. The literature check, and wording no stronger than "first computed for Titania and Oberon".
9. The J2 pairing corrected in `v4_uranus.py` and wherever it was copied.

## 12. How far the reviewer's own numbers can be relied on

* The reviewer's planar model was derived independently and agrees with the repository's
  right-hand side to rounding, so sections 2 to 5 rest on two codes, three counting the
  coordinator's. Its integrations use DOP853 and Radau from scipy, as the build's do; the
  integrator library is shared, the equations and frames are not.
* The reviewer's URA111 force model has one positive control (the moons' own motion, which also
  discriminates the J2 value) and one internal cross-check (its circular limit against the
  planar code). It has no tests, is not in the repository, and was run at two epochs. The
  ephemeris file is the copy of `ura111.bsp` shipped with the local GMAT installation.
* "Converged" in section 6 means junction discontinuities below one metre and one hundredth of
  a millimetre per second at every node, recomputed with a second integrator. It is the danger
  signal this project has learned to distrust, so the things that would make it hollow were
  checked: the force model against the ephemeris, the discontinuities serially and with a
  different method, each inner flyby by one continuous 61.6-day propagation, and the flyby
  altitudes by a root-finder on the range rate against the URA111 moon positions. What was not
  checked: an independent code base, a third epoch, arcs longer than six cycles.
* Not done at all: Titania's GM varied (the reviewer's continuation failed for a bookkeeping
  reason); the build's continuation in mass and its sensitivity table recomputed; the registry
  GM and radius values compared with a paper; the literature check; a navigation model.
* The J2 attribution: no Jacobson paper is listed in the project's corpus index, so the
  "Jacobson 2014 Table 4" attribution in `v4_uranus.py` was not grounded against a held source.
  The reviewer's statement of which paper gives which value is from memory and must be checked
  against the papers; the computed part is only that URA111 obeys J2 R^2 near 2.29e6 km^2.
