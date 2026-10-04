# #884 adversarial review -- Sun-forced Earth-Moon cyclers

**Date:** 2026-10-04. **Reviewed:** `docs/notes/2026-10-04-884-sun-forced-em-cyclers.md`,
`src/cyclerfinder/search/sun_forced_periodic_884.py`, its tests and driver, and
`data/found/884_sun_forced_em_cyclers/`. **Reviewer's code:** written separately from the build (own
equations of motion, own multiple shooting, own pseudo-arclength step control, own family walk); it
lives in the session scratch directory and is not committed, so section 12 describes each
computation well enough to repeat it. Nothing in the repository was edited.

Every finding is tagged COMPUTED (I ran it; numbers given), READ (in the code, note or data) or
INFERRED (argument, no run).

## 1. Headline

The numerics of the build are sound. The model is not. The project's bicircular model
(`core/bcr4bp.py`, and the #884 module that is gated against it) moves the Sun the wrong way round
in the rotating frame. Every #884 orbit is therefore a periodic orbit of a system that is not the
Sun-Earth-Moon bicircular problem, and the quantitative results (which phases reach full Sun mass,
the three folds, the Floquet magnitudes, the periselenes, the Melnikov amplitudes, the two
"linearly stable" orbits) do not carry over. I re-ran the a = 1 and a = 2 part of the survey with the
Sun's sense corrected; the outcome is different in every one of those respects (section 3).

Three further corrections stand whichever way the Sun turns:

- The "C21 3D corridor at 2/1" member is a planar orbit whose minimal period is one synodic month.
  It was reached by a branch switch at a pitchfork and is not a member of the 3D family.
- The Casoliva 2:1(b) low-perilune parent is itself below the lunar surface in the three-body
  problem (periselene 1,415 km from the Moon's centre). The note's 1,810 km is a sampling artefact,
  and the Sun does not cause the impact.
- The "positive control against Brown et al." cannot fail for a wrong forcing and is not a control.

A second model defect surfaced while preparing the coherent-model spot check: `core/qbcp.py`
evaluates alpha_2 as a cosine series and alpha_3 as a sine series, the reverse of what the tables
are, and the shipped coherent model is not time-reversible (section 9). The spot check was run with
that defect repaired in process: three corrected-sense orbits do have periodic counterparts there,
strongly deformed. It is conditional on the rest of that module, which the coordinator has since
recorded as not validated (#892; the bicircular defect is #891).

## 2. Verdicts

| # | Claim | Verdict | Basis |
|---|---|---|---|
| 1 | The stored objects are periodic orbits of the project's bicircular model at full Sun mass | ESTABLISHED | Coordinator's check; my own equations of motion close all 39 stored orbits to 2.0e-10 or better (COMPUTED) |
| 2 | That model is the Sun-Earth-Moon bicircular problem | WRONG | Sun revolves counter-clockwise in the rotating frame; the physical Sun regresses. Inertial-frame test: 4e-12 agreement with the sense reversed, 1e-2 to 4e-2 disagreement as shipped (COMPUTED, section 3) |
| 3 | Each equivalent is the continuation of the stated three-body member | ESTABLISHED (in the project model) | Backward continuation with my own stepping, 39 of 39 land on the stated member (COMPUTED, section 4) |
| 4 | Each member belongs to the stated catalogue family | ESTABLISHED WITH CORRECTION | 14 of 15 members reproduced from the catalogue rows' own initial conditions to 2e-12 or better; the "C21 3D corridor 2/1" member is planar and is not on the 3D family (COMPUTED, section 4) |
| 5 | "Every catalogued family tested at a low-order member has equivalents at physical Sun mass" | ESTABLISHED WITH CORRECTION | True of the project model apart from the folds. In the corrected-sense model every phase of every cycler-class member I ran reaches full Sun mass, apart from one unfinished branch of the sub-surface Casoliva member (COMPUTED, single implementation) |
| 6 | Three branches fold back and reconnect to another three-body orbit | ESTABLISHED for the project model; WRONG as a statement about the Sun | My code reproduces the folds at eps = 0.3893 and 0.3989, det(M - I) changes sign there, and the landing orbits are three-body periodic to 6e-13. With the Sun's sense corrected none of the three folds occurs (COMPUTED) |
| 7 | All cycler-class equivalents are linearly unstable | ESTABLISHED (dominant multiplier only) | Holds in both senses. The counts of unstable directions and the "unit-circle pair" statements are not reliable above about 1e7 (COMPUTED, section 7) |
| 8 | The only linearly stable forced orbits are R21-S at 1:2 and Casoliva 1:2(d) at 3:2 | WRONG for Casoliva 1:2(d) | Corrected sense: that orbit has maximum multiplier 13.8 and its other phase folds at eps = 0.608. R21-S stays stable (COMPUTED) |
| 9 | Periselene stays inside the lunar sphere of influence, so all stay cycler-class | ESTABLISHED WITH CORRECTION | Label is inherited from the parent and says nothing about Earth; two of the "cycler" entries belong to rows the catalogue classes as not cyclers (section 5) |
| 10 | The Casoliva 2:1(b) low member's equivalents are pushed below the lunar surface | WRONG | The parent is at 1,415.4 km, already 322 km below the surface (COMPUTED, section 8) |
| 11 | Harmonic argument (which resonances are quadrupole, octupole, flat) and oddness at symmetric phases | ESTABLISHED | Checked by hand and against the computed zeros in both senses (section 6). The quoted amplitudes are project-model values and change by factors of 0.07 to 2000 |
| 12 | "Positive control against Brown et al.: reproduced" | NOT ESTABLISHED | It tests the code's symmetry, not agreement with the paper, and passes identically with the Sun running backwards (section 6) |
| 13 | Four (a = 1) or two (a = 2) equivalents per member | ESTABLISHED WITH CORRECTION | For a = 1: three trajectories up to mirror symmetry, two of them 6 to 62 km apart (section 6) |
| 14 | Orbits are computed to the accuracy claimed | ESTABLISHED | Nodes reproduce to 2e-9 under a 2N-segment re-solve; closure 1e-10 is meaningful as a multiple-shooting statement only (section 7) |
| 15 | Statements about a = 3 (8:3) members | NOT REVIEWED | Incomplete in the build; all a = 3 amplitudes are project-model values |

**Single most serious weakness:** verdict 2. The headline question "which families survive the Sun"
was answered in a model whose Sun has an inertial period of 14.19 days.

## 3. The Sun's sense of rotation (question 7, and the reason for everything else)

**READ.** `core/bcr4bp.py::_sun_position` returns
`(a_S cos(theta0 + omega_S t), a_S sin(theta0 + omega_S t))` with `omega_S = +0.925195985520347`, in
a frame whose Coriolis terms (`+2 vy`, `-2 vx`) are those of a counter-clockwise rotation at rate 1.
The #884 module uses the same expression (`th = th0 + w_s * t`). Andreu and Simo et al. write the
Sun at `(a_S cos theta, -a_S sin theta)` with `theta = omega_S t`.

**INFERRED.** The Sun moves prograde around the Earth-Moon barycentre at inertial rate
n_S = 1 - omega_S = 0.0748 (the synodic month is longer than the sidereal month precisely because
the Sun moves the same way as the Moon). In a frame turning at rate 1 its angle is therefore
theta0 - omega_S t: it regresses. With the plus sign the Sun's inertial rate is 1 + omega_S.

**COMPUTED.**

| Check | Result |
|---|---|
| Sun's inertial period implied by n_S = 1 - omega_S | 365.243 d |
| Sun's inertial period implied by the shipped sign | 14.19 d |
| Bicircular model integrated in a barycentric non-rotating frame (Earth and Moon on circles at rate 1, Sun prograde with period 365.243 d, direct plus indirect term), transformed to the rotating frame, against the rotating model with -omega_S, after 6 TU, two Sun phases | agree to 4.7e-12 and 4.1e-12 |
| Same, against the rotating model as shipped (+omega_S) | differ by 3.7e-2 and 1.2e-2 |
| For scale: the three-body problem with no Sun at all against the same truth | differs by 2.3e-2 and 2.4e-2 |
| My equations of motion against `core.bcr4bp.bcr4bp_eom` and `core.cr3bp.cr3bp_eom`, 200 random states | 1.8e-15 |
| Sun angle in `core/qbcp.py` at t = 0, 0.2, 0.5, 1.0 TU | 180.0, 169.2, 153.0, 126.4 deg (clockwise) |
| Sun angle in `core/bcr4bp.py` at the same times | 0, 10.6, 26.5, 53.0 deg (counter-clockwise) |

So the shipped bicircular model is as far from the bicircular problem as having no Sun, and the
project's two Sun models disagree with each other on the sense. No symmetry maps one sense to the
other: the reversing symmetry (x, -y, z, -vx, vy, -vz, -t) sends theta0 to -theta0 and keeps the
sign of omega_S, and a plain mirror flips the Coriolis terms.

**What survives the correction (INFERRED, then COMPUTED):** the argument that Melnikov zeros sit at
the symmetric phases uses only reversibility, which both senses have. The zeros at 0 and pi (and
pi/2 for a = 2) are the same in both senses for all 16 members with a = 1 or 2. Everything that
depends on the size or sign of the forcing along the orbit changes, because the two senses pick out
different Fourier coefficients of the orbit (harmonic +k against harmonic -k).

**Survey with the sense corrected (COMPUTED; my code only; each resulting orbit re-integrated
segment by segment in the non-rotating frame, closure 1.4e-13 to 8.0e-12, and 4.1e-10 for the
sub-surface Casoliva mirror pair).**
"lam" is the largest Floquet magnitude over the forced period; periselene in km; phases in radians at
the member's x-axis crossing. For a = 2, theta0 and theta0 + pi are the same orbit one lap later, so
two phases are listed.

| Member | Melnikov amplitude, project / corrected (ratio) | Project model (stored, re-measured) | Corrected sense |
|---|---|---|---|
| C11 2/1 | 5.20e-2 / 9.59e-3 (0.18) | 4/4 reach; lam 4.0e5 to 5.0e5; periselene 24,678 to 26,303 | 4/4 reach; lam 3.4e5 to 4.2e5; 24,999 to 25,931 |
| C11 3/2, C = 3.15107 | 5.25e-2 / 1.33e-2 (0.25) | 2/2; lam 3.8e7, 1.4e7; 21,854, 23,998 | 2/2; lam 6.6e5, 3.9e7; 21,080, 25,104 |
| C11 3/2, C = 3.09204 | 1.37e-2 / 5.42e-2 (3.9) | 2/2; lam 3.2e10, 9.9e9; 21,670, 22,114 | 2/2; lam 2.3e10, 1.5e10; 20,531, 23,396 |
| C11 5/2 | 2.11e-3 / 5.14e-2 (24) | 1/2; theta0 = 0 folds at 0.399; lam 2.5e14; 19,726 | 2/2, no fold; lam 8.9e13, 1.0e12; 18,514, 21,010 |
| C21 planar 3/1 | 5.66e-2 / 1.33e-2 (0.24) | 2/4; pi/2 pair folds at 0.389; lam 9.1e5; 24,404 | 4/4, no fold; lam 4.0e3, 5.6e4; 26,066, 25,103 to 25,109 |
| "C21 3D" 2/1 (planar, see section 4) | 3.68e-2 / 3.89e-2 (1.06) | 4/4; lam 6.6e4, 1.7e5; 40,404, 43,265 to 43,269 | 4/4; lam 1.0e5, 1.1e5; 40,782 to 40,787, 42,778 |
| C32 5/2, C = 3.15168 | 2.00e-2 / 6.48e-2 (3.2) | 2/2; lam 1.1e7, 2.5e7; 21,741, 23,312 | 2/2; lam 9.8e6, 1.9e7; 21,806, 23,199 |
| C32 5/2, C = 3.18010 | 1.72e-2 / 5.23e-2 (3.0) | 2/2; lam 1.8e4, 2.9e3; 8,854, 8,289 | 2/2; lam 2.2e4, 3.1e3; 10,048, 7,670 |
| C31 5/2 | 1.09e-2 / 1.01e-1 (9.3) | 2/2; lam 4.4e7, 3.1e7; 14,006, 14,323 | 2/2; lam 2.6e7, 4.9e7; 14,903, 13,481 |
| Casoliva 2:1(b) 1/1, C = -0.02092 | 3.31e-4 / 3.15e-2 (95) | 4/4; lam 7.8; 13,367 to 13,430 | 4/4; lam 6.3, 9.5; 11,749 to 11,759, 15,428 |
| Casoliva 2:1(b) 1/1, C = 0.56729 | 4.75e-3 / 1.65e-1 (35) | 4/4; lam 950; 1,440 to 1,453 | mirror pair reaches (lam 990; 1,349); theta0 = 0 reaches eps = 1.0003 with residual 1.3e-11, just over my 1e-11 tolerance, no diagnostics taken; theta0 = pi stopped by the time limit at eps = 0.991 |
| R52-S 2/1 | 8.04e-3 / 2.13e-2 (2.7) | 4/4; lam 9.2, 10.4; 74,869 to 75,133 | 4/4; lam 9.3, 7.6; 80,217 to 80,243, 70,196 |
| R52-S 5/2 | 6.47e-3 / 1.22e-1 (19) | 2/2; lam 1.7e6, 2.7e6; 79,747, 78,616 | 2/2; lam 2.2e6, 2.0e6; 76,076, 82,484 |
| R21-S 1/2 | 9.19e-4 / 6.02e-5 (0.066) | 2/2; lam 1.16 and 1.000 (stable) | 2/2; lam 1.04 and 1.000 (stable) |
| Casoliva 1:2(d) 3/2 | 3.51e-6 / 6.94e-3 (1980) | 2/2; lam 1.01 and 1.000 (stable); periselene about 348,000 | 1/2; theta0 = 0 reaches with lam 13.8, periselene 594,530 km, perigee 289,792 km; theta0 = pi/2 folds at eps = 0.608 and lands on a C = 3.15154 three-body orbit (periselene 107,648 km) |

Reading the table: the qualitative statement "cycler-class members at a = 1, 2 keep periodic
equivalents at full Sun mass" becomes stronger, not weaker, with the right Sun (no cycler-class
fold in the sample). The specific content of the #884 note does not survive: the fold list is
replaced by a different one, the stable Casoliva 1:2(d) orbit is gone, the periselene ranges move by
up to about 2,500 km (R52-S by about 5,000 km), and Melnikov amplitudes move by up to three orders
of magnitude.

**Control on my own instrument (COMPUTED).** The same forward pipeline run with the project's sense
reproduces the build's stored eps = 1 nodes to 3e-10 (R21-S), 8.5e-9 (Casoliva 1:2(d)), 1.2e-7
(C11 2/1), 4.5e-7 (C21 3:1) and 1.6e-5 (C11 5/2). These residual differences are the build's
theta0 being off the exact symmetric phase by up to 1.8e-5 rad, which is a time shift of the same
orbit and not an error (section 7). It also reproduces both project-model folds.

## 4. Is each equivalent the continuation of the stated member? (question 1)

**Backward continuation (COMPUTED).** From each of the 39 stored eps = 1 orbits I continued toward
eps = 0 in pseudo-arclength with my own code: maximum step 0.015 (the build used 0.02), accepted
only if Newton converged within 6 iterations and the cosine between successive unit tangents
exceeded 0.9 (a different guard from the build's "within half a step of the predictor"). This covers
the 22 branches the build did not reverse-check. Result: 39 of 39 run monotonically to eps = 0 with
no fold (69 to 78 steps each). At eps = 1e-4 the first node is within 4e-8 to 1.5e-4 of the
three-body orbit; after a symmetric half-period correction at T* = P/a in the three-body problem
(residual 2.3e-11 or better) the landing orbit equals the build's start state to 5.3e-13 or better
in every case. One branch (C11 3/2 at C = 3.09204, theta0 = 0) did not converge at eps = 1e-4 and
was corrected from eps = -0.008 instead (distance 5.2e-4 before correction, same landing).

**Against the catalogue rows' own initial conditions (COMPUTED).** I walked each family from the
catalogue row's `state_nd` and `period_nd` (not from anything the build stored) with my own
pseudo-arclength walk in (x0, [z0], vy0, T), tangent-cosine guard 0.9 to 0.95, and compared every
crossing of the target period with the landing orbit.

| Member | Reached from row | Max difference from the landing orbit |
|---|---|---|
| C11 2/1 (x0 = -0.76485) | `ross-rt-em-cycler-11-2025` | 1.3e-15 |
| C11 3/2 (x0 = -0.77034) | `ross-rt-em-cycler-11-2025` | 2.1e-15 |
| C11 3/2 (x0 = -0.85152) | `braik-ross-c11a-cycler-2026` (and the Ross-RT row) | 3.2e-15 |
| C11 5/2 (x0 = -0.84540) | `ross-rt-em-cycler-11-2025` | 5.8e-15 |
| C21 planar 3/1 | `ross-rt-em-cycler-21-2025` | 1.1e-14 |
| C32 5/2, both members | `braik-ross-c32-cycler-2026` | 2.9e-13, 8.9e-14 |
| C31 5/2 | `ross-rt-em-cycler-31-2025` | 8.9e-16 |
| R21-S 1/2 | `braik-ross-planar-r21-s-corridor-2026` | 1.3e-15 |
| R52-S 2/1, 5/2 | `braik-ross-planar-r52-s-corridor-2026` | 5.6e-15, 6.1e-16 |
| Casoliva 1:2(d) 3/2 | `casoliva-1-2d-em-resonant-po-2010` | 3.6e-15 |
| Casoliva 2:1(b) 1/1, both members | `casoliva-2-1b-em-resonant-po-2010` | 1.6e-12, 9.3e-13 |
| "C21 3D corridor" 2/1 | not reached from `braik-ross-c21-3d-corridor-05-2026` | see below |

The walk from the Ross-RT C11 row passes through T = 9.69108 at (x0, vy0) = (-0.81164, -0.11859),
which is the Braik-Ross C11a row state (-0.8116407, -0.1185906): those two rows are one family, as
the note says. I did not check C11b or the C32 = Ross-RT C32 identity.

**The "C21 3D corridor at 2/1" member is not a member of the 3D family (COMPUTED).**

- The stored start state is (0.87858, 0, 4e-40, 0, -0.33677, 0): planar.
- It closes after T/2 to 1.4e-12. Its minimal period is 6.791194 TU, one synodic month, not two. It
  makes one lunar pass (42,003 km) and three Earth-distance minima per synodic month.
- My walk of the 3D family from the catalogue row reaches the plane at T = 13.701 TU (x0 = 0.87645,
  vy0 = -0.32147, z0 passing through zero) and continues into the mirror (z to -z) branch; the 3D
  family never reaches T = 13.582. The build's walk
  (`families/braik-ross-c21-3d-corridor-05-2026.json`, READ) continued to T = 12.56 with x0 = 0.973,
  vy0 = -1.63, C = 1.95: it had left the 3D family at the pitchfork where that family meets the
  plane and was following the planar branch with z0 = 0.
- Walking the planar family upward from this member gives (x0, vy0) = (0.87296, -0.25693) at
  T = 19.440, against the Ross-RT C21 row's (0.72373, 0.41377), and (0.87682, -0.27875) at T = 3 Tg,
  against the planar C21 3/1 member's two crossings (0.69676, 0.48029) and (1.05551, 0.44558). It is
  not the planar C21 family either. I did not identify it.

Consequences: the table row "C21 3D corridor (#682 / #438 spatial), 2/1 (1)" should read "a planar
orbit of period 1 Tg on an unidentified planar family that the C21 3D family bifurcates from, taken
twice". Its ratio is 1/1, not 2/1. The quoted Floquet magnitudes (1.1e5 for the parent, 6.6e4 to
1.7e5 forced) are squares of the per-period values (about 330 for the parent). The half-step guard
does not protect a walk at a pitchfork, because both branches lie within half a step.

## 5. Is "cycler-class" earned? (question 2)

**COMPUTED** (project model, all 39 stored orbits and their parents; every local minimum of the
distance to the Moon and to the Earth over one forced period, refined on the dense solution, against
the parent over the same time).

- 36 of 39 forced orbits have the same number of lunar-distance minima and the same number inside
  66,183 km as their parent over the same time. The exceptions:
  - C11 3/2 (C = 3.15107), theta0 = 0: 2 minima (21,854 km each) against the parent's 6 (39,381 /
    23,253 / 39,381 per lap). One encounter per lap is kept; the triple-minimum shape is lost.
  - C32 5/2 (C = 3.18010), theta0 = pi/2: 4 inside the sphere of influence against 6 (the parent's
    middle minimum at 45,158 km is gone).
  - C32 5/2 (C = 3.15168), theta0 = 0: 10 minima against 14, same 6 inside.
- Earth-distance minima differ in count on 5 orbits (C11 2/1 mirror pair 2 against 4; C21 3/1 7
  against 5; Casoliva 1:2(d) 4 against 2), always by the appearance or loss of a shallow minimum.
- Largest periselene change, forced against parent, among cycler-class orbits: 1,604 km (C21 3/1).

So in the project model the forced orbits are the same kind of object as their parents. That is as
far as the label goes:

- **The label says nothing about Earth (COMPUTED).** The closest approach to Earth of the
  cycler-class orbits is 94,656 km (C31) to 226,380 km (C21 3/1); C11 2/1 never comes within
  224,000 km of Earth. Only the sub-surface Casoliva member comes closer (15,509 km). These are
  inherited from the parents (within a few percent). An "Earth-Moon cycler" in this table is an
  orbit with a periselene inside the lunar sphere of influence; none of them visits low or even
  geostationary Earth altitude.
- **Two table entries carry a label their catalogue row does not (READ).** `data/catalogue.yaml`
  classes `casoliva-2-1b-em-resonant-po-2010` as `resonant_po` ("NO lunar encounter", periselene
  92,590 km). The members the note calls cycler-class (13,843 km and the sub-surface one) are
  distant members of that family (C = -0.02 and 0.57 against the row's 1.196), not the catalogued
  orbit. The same applies to the planar orbit of section 4: the 3D row is `quasi_cycler` at V0 and
  the #438 spatial row is `resonant_po` with closest lunar pass 122,628 km.
- **READ.** The note says 66,182.9 km is "the lunar Laplace SOI used by `data/validate.py`". The
  number does not occur in `src/cyclerfinder/data/validate.py`; in `src/` it occurs only in the #884
  module. It is the figure quoted in catalogue comments.

## 6. Multiplicity, the Melnikov argument and the positive control (questions 3 and 6)

**Multiplicity (COMPUTED, stored orbits).** For the five a = 1 members with four stored phases:

| Member | Mirror pair: max difference between one orbit and the mirror image of the other | theta0 = 0 against theta0 = pi: max node difference (position) |
|---|---|---|
| C11 2/1 | 2.7e-7 | 1.1e-4 (43 km) |
| "C21 3D" 2/1 | 7.1e-9 | 3.9e-5 (9 km) |
| Casoliva 2:1(b), C = -0.02 | 2.3e-9 | 3.5e-5 (6 km) |
| R52-S 2/1 | 2.8e-9 | 4.5e-4 (62 km) |
| Casoliva 2:1(b), C = 0.57 | 3.2e-8 | 1.7e-4 (15 km) |

So the phases near pi/2 and 3 pi/2 are one trajectory and its mirror image under the reversing
symmetry, and the phases 0 and pi are two self-symmetric orbits that differ only at octupole order.
Per a = 1 member there are three trajectories up to symmetry and two at quadrupole resolution. For
a = 2, theta0 and theta0 + pi are the same orbit started one lap later, so two per member, which
the note counts correctly. Of the 39 stored orbits, 34 are distinct up to mirror symmetry and 28
are distinct at quadrupole resolution. A separate point: a forced periodic orbit is a fixed point of
the stroboscopic map for every theta0 (the same orbit sampled at another time), so "phase" labels a
point on the orbit, not a different orbit.

**Harmonic argument (INFERRED, checked against COMPUTED zeros).** The direct plus indirect solar
potential is mu_S / a_S times the sum over l >= 2 of (r / a_S)^l P_l(cos gamma); P_2 carries
Sun-angle harmonics 0 and 2, P_3 carries 1 and 3, P_4 carries 0, 2 and 4. With T* = (n/a) Tg and
gcd(n, a) = 1 only harmonics that are multiples of a contribute, so a = 1 and 2 are quadrupole, a = 3
octupole, a = 6 of relative order (r / a_S)^4. Oddness of Mel in theta0 at a perpendicular crossing
follows from the reversing symmetry in either sense. This is correct. In my scans every a = 2 member
has its zeros at multiples of pi/2 to the scan resolution in both senses, and the a = 1 members show
the octupole offset from pi/2 (2e-4 to 1.2e-2 rad in the project model, 1e-4 to 1.8e-3 rad with the
sense corrected). One gap: the argument gives which zeros exist, not their slopes or the amplitude,
and those are the sense-dependent quantities.

**The positive control (INFERRED and COMPUTED).** It does not test agreement with Brown et al.:

- Brown et al. print no initial conditions, and their starting orbits are different members from
  the ones used here (the note says so).
- Zeros at the symmetric phases are forced by reversibility for any reversible periodic forcing.
  The test can fail only if the code breaks the symmetry.
- "Exactly two zeros per pi" is what any odd, pi-periodic function dominated by one harmonic has.
- "Only Melnikov zeros continue" is the Lyapunov-Schmidt statement itself.
- The decisive point: the same zero set at the symmetric phases appears for all 16 members with the
  Sun running the wrong way and the right way. A check that a physically wrong forcing passes is
  not a control on the forcing.

What it does establish: the Melnikov quadrature, the eps-sensitivity equation and the continuation
are mutually consistent, and the code respects the reversing symmetry. What it does not establish:
that any number in Brown et al. is reproduced, or that the forcing is the Sun's.

## 7. Accuracy and the stability numbers (question 5)

**COMPUTED** on six stored orbits spanning maximum multipliers 7.8 to 2.5e14.

| Orbit | Stored residual, my integrator at 1e-13 | Re-solve with 2N segments: max change of the common nodes | Single integration of node 0 over the whole period: closure | Largest times smallest Floquet magnitude (should be 1) |
|---|---|---|---|---|
| C11 5/2, theta0 = pi/2 (2.5e14) | 7.3e-12 | 1.3e-10 | 2.66 | 3.9e12 (N), 4.0e13 (2N) |
| C11 3/2, C = 3.092, theta0 = 0 (3.2e10) | 4.0e-13 | below 1e-12 | 6.7e-4 | 1.1e5, 7.6e5 |
| C11 3/2, C = 3.092, theta0 = pi/2 (9.9e9) | 7.4e-12 | 3.9e-10 | 4.3e-3 | 2.8e5, 1.2e5 |
| C31 5/2, theta0 = pi/2 (4.4e7) | 4.3e-11 | 2.3e-9 | 5.0e-4 | 47, 50 |
| C11 2/1, theta0 = 0 (5.0e5) | 3.4e-11 | 1.9e-9 | 2.7e-7 | 1.0002, 0.9999 |
| Casoliva 2:1(b), theta0 = pi (7.8) | 4.3e-13 | below 1e-12 | 3.9e-12 | 1.0000 |

- The node sets are well determined: the smallest singular value of the shooting Jacobian is 6e-6
  to 3e-3, and doubling the segments moves the nodes by 2e-9 at most. The orbits are computed to the
  accuracy the note implies.
- A closure of 1e-10 is a statement about segments. None of the orbits above 1e7 can be integrated
  for one period in double precision (closure 5e-4 to 2.7 in a single shot). That is expected, and it
  is also the practical meaning of the instability: e-folding times are about 4.5 days for both
  C11 2/1 and C11 5/2.
- The largest multiplier is stable under N to 2N (2.5175e14 both ways). The others are not: for the
  2.5e14 orbit the second magnitude is 5.20 with N segments and 1.11 with 2N, and the reciprocal
  partner of the largest comes out 1e12 to 1e13 too big. Lap-shifted copies of one a = 2 orbit give
  3 and 1 "unstable" magnitudes (C11 3/2, corrected-sense run). So "unstable(n)" and the statements
  about a unit-circle pair against a real pair are reliable only where the product check passes,
  which is up to about 1e6. The alternation claim in section 4 of the note is checked on the
  Lyapunov controls only.
- The stored theta0 of "symmetric" branches is off the exact symmetric phase by up to 1.8e-5 rad
  (C11 5/2: 1.570814 against pi/2; self-symmetry defect of the stored C11 2/1 orbit 3.3e-7). Because
  any theta0 gives a fixed point of the same orbit, this is a time shift, not an error.

## 8. Folds, and sub-surface members (questions 4 and 8)

**Folds (COMPUTED).** In the project model my own forward continuation from the three-body member
folds at eps = 0.3893 (C21 3/1, both members of the pi/2 pair) and 0.3989 (C11 5/2, theta0 = 0);
det(M - I) changes sign across each fold (for C21: -1.0e3, +6.5e1 either side), so they are
saddle-node turning points, not corrector failures. The branches return to eps = 0 and, after a
Newton solve at eps = 0, are periodic orbits of the three-body problem by the separate core
propagator (`core.cr3bp.cr3bp_eom`, segment closure 1.7e-13 and 6.3e-13), with Jacobi constant
3.121403 (closed over 3 Tg, periselene 25,760 km) and 3.079993 (closed over 5 Tg, minimal period not
determined, periselene 20,434 km), at distances 0.18 and 0.33 from the start orbit. The note's description is right for the
project model. I did not identify the landing orbits' families. With the Sun's sense corrected
neither fold exists (section 3).

**Sub-surface and sub-Earth-radius passes (COMPUTED).**

- Casoliva 2:1(b), C = 0.56729: the three-body parent's periselene is 1,415.4 km from the Moon's
  centre. The note and `melnikov/*.json` give 1,810 km. The build samples the distance at 400 points
  per TU (one point per 15.6 minutes; driver line 167 for the parent, `orbit_diagnostics` for the
  forced orbits) and takes the smallest sample; at this pass that misses the minimum by 395 km. The forced orbits are at 1,440 to 1,453 km (the forced values happen to be
  sampled well). So the Sun raises the periselene by 25 to 38 km; it does not cause the impact. The
  parent is an impact orbit of the three-body problem, and the member is unusable either way.
- The other Casoliva 2:1(b) parent is at 13,817 km (stored 13,843 km), and its forced orbits near
  pi/2 at 13,430 km (stored 13,441 km). All other stored periselenes agree with mine to 1 km.
- Nothing else is below the lunar radius. The smallest Earth distance anywhere is 15,509 km
  (Casoliva low member); nothing is inside the Earth.

## 9. Model honesty and the coherent model (question 7)

- **The bicircular model as shipped is not a Sun model** (section 3).
- **Even with the sense corrected (INFERRED)** the bicircular model is not a solution of the
  three-body problem for Sun, Earth and Moon: the Moon does not feel the Sun. The omitted variation
  is of the same order as the tide on the particle. A periodic orbit in it shows that the three-body
  orbit has a periodic neighbour under one particular periodic idealisation of the Sun.
- **Against the real system (INFERRED)** the lunar eccentricity (0.055) is a perturbation several
  times larger than the solar tide (m^2 = 0.0065) and has the anomalistic month as its period, which
  is not commensurate with the synodic month. In an ephemeris nothing here is periodic; the most one
  can hope for is a nearby bounded trajectory that needs station-keeping, and with multipliers of
  1e5 to 1e14 per period the station-keeping cost is the real question. "Survives the Sun" bears
  little weight beyond "has a periodic counterpart in a periodic model".
- **The coherent model as shipped is not usable (COMPUTED).** `core/qbcp.py::evaluate_alphas`
  passes `is_even=True` for alpha_2 and `is_even=False` for alpha_3. The tables say the opposite:
  `_COEFFS_ALPHA2` has a zero constant term like the other sine tables (alpha_5, alpha_8) and
  `_COEFFS_ALPHA3` starts at 0.99999 like a cosine table. The true model has the reversing symmetry,
  which requires alpha_2 odd and alpha_3 even. Test: |R x(T; y0) - x(-T; R y0)| for two states and
  T = 1 and 3 TU is 3.3e-2 to 1.1e-1 as shipped, and exactly 0 with the two parities exchanged
  (in-process patch only). As shipped alpha_2(0) = -0.0137 (should vanish at syzygy) and
  alpha_3(0) = 1.0000 (should be about 1.0196). This is separate from the coefficient corrected in
  commit 6bda1772.
- **Spot check in the coherent model with the parities repaired in process (COMPUTED;
  conditional).** The corrected-sense bicircular model is the coherent model's equations with
  alpha_1 = alpha_3 = alpha_6 = 1, alpha_2 = 0 and single-harmonic alpha_4, 5, 7, 8 (Sun at angle
  pi - omega_S t, so the bicircular phase theta0 = pi matches the coherent model's t = 0). I blended
  every Fourier table linearly from those values (lam = 0) to the shipped tables (lam = 1), with
  alpha_2 as a sine and alpha_3 as a cosine series, and followed three corrected-sense theta0 = pi
  orbits in lam by multiple shooting with the module's own propagator and STM. Instrument check: at
  lam = 0 the alpha form reproduces the bicircular orbits to 1.4e-9, 9.8e-9 and 7.2e-9.

  | Orbit (corrected sense, theta0 = pi) | Reaches lam = 1 | Closure, separate propagation at 1e-13 | Largest multiplier, bicircular to coherent | Largest node displacement | Smallest Moon distance, bicircular to coherent |
  |---|---|---|---|---|---|
  | Casoliva 2:1(b), C = -0.02, P = 1 Tg | yes, 11 steps | 3.4e-11 | 9.5 to 10.4 | 0.090 (34,000 km) | 11,749 to about 10,800 km |
  | C11 2/1, P = 2 Tg | yes, 11 steps | 2.1e-12 | 3.4e5 to 6.8e5 | 0.074 (28,000 km) | 25,005 to about 21,600 km |
  | C21 planar 3/1, P = 3 Tg | yes, 25 steps | 3.3e-12 | 5.6e4 to 5.1e6 | 0.213 (82,000 km) | 25,103 to about 24,400 km |

  All three stay self-symmetric (defect below 1e-10) and the node displacement grows smoothly with
  lam, so no jump is evident. The coherent-model distances are read off integrator steps in pulsating
  coordinates and are good to a few hundred km only. What this shows: going from the bicircular to a
  coherent forcing is not a small step. It moves these orbits by tens of thousands of km and changes
  the C21 multiplier ninety-fold. That is as large as, or larger than, what switching the Sun on did
  in the corrected-sense bicircular model (node displacement 0.059, 0.027 and 0.030 for the same
  three orbits). What it does not show: anything about the true coherent model. The
  coordinator's check (#892) finds that the parity-repaired module still misses a published point by
  2.3e-2, so the tables or equations have a further problem and these three orbits are orbits of an
  unvalidated model.

## 10. Novelty: what the mathematics makes unsurprising (question 10)

- That a periodic orbit of an autonomous Hamiltonian system, lying in a one-parameter family whose
  period varies, has periodic continuations under a small periodic forcing at the members whose
  period is commensurate with the forcing, for small forcing, at simple zeros of a bifurcation
  function of the phase: standard (Poincare continuation, subharmonic Melnikov theory; it is the
  framework Brown et al. apply).
- That for a reversible forcing and a symmetric orbit those zeros sit at the symmetric phases, so
  existence at small forcing needs no computation: standard, and Brown et al.'s Proposition 3.
- That neighbouring zeros alternate in stability along the phase direction, and that an unstable
  parent gives unstable equivalents: standard.
- That branches can fold and join two three-body orbits: reported by Brown et al. for their model.

What is not given by a theorem is whether a branch reaches the physical Sun mass, where it folds,
and what the orbit looks like there. That is a computation, it is model-specific, and the note's
version of it is in the wrong model. A corrected-sense table of which Earth-Moon cycler families
keep periodic equivalents at full Sun mass in the bicircular problem, with their periselenes and
multipliers, would be an honest computed result of the "census" kind. Whether it is new is the
literature agent's question; the mathematics alone suggests it is an application of a known method
to orbits not previously listed, not a new phenomenon.

## 11. Wording corrections (question 9)

| Note text | Problem | Proposed replacement |
|---|---|---|
| Title: "Which catalogued Earth-Moon cycler families survive the Sun (BCR4BP)" | The model's Sun turns the wrong way; "survive the Sun" overclaims even in a correct bicircular model | "Periodic continuations of Earth-Moon cycler family members under a bicircular solar forcing (project model; Sun sense to be corrected)" |
| "has Sun-forced 'dynamical equivalents' at physical Sun mass in the bicircular model" | Not the bicircular problem | "... at full solar mass parameter in `core/bcr4bp.py` as shipped, whose Sun moves counter-clockwise in the rotating frame" until rerun |
| "The survey covered 13 commensurate members of 9 families" | The table has 19 rows, 16 of them with a = 1 or 2 (15 continued to eps = 1 on at least one phase); one of those is a planar period-Tg orbit, not on the 3D family | State the count the table shows, and "plus one planar orbit of period 1 Tg reached through a pitchfork from the C21 3D family" |
| "For Ross-RT C21 (3:1), the two equivalents at the near-perpendicular phase fold back at eps = 0.389 ..." and the C11 5:2 fold | True of the shipped model only | Add "in the shipped model; with the Sun's sense corrected all four (two) phases reach eps = 1" |
| "The only linearly stable forced orbits found belong to ... R21-S at 1:2 and Casoliva 1:2(d) at 3:2" | Casoliva 1:2(d) is unstable (13.8) with the sense corrected | "R21-S at 1:2 only" after rerun |
| "All forced cycler-class equivalents are linearly unstable" | Fine, but the table's unstable-direction counts are not reliable above 1e7 | Keep; add "largest multiplier only; subdominant multipliers are not resolved by the plain STM product above about 1e7" |
| "One exception: the Casoliva 2:1(b) low-perilune member's forced equivalents reach periselene 1,440-1,453 km. That is below the lunar radius" and table "1,810 km" | The parent is at 1,415 km; 1,810 km is a sampling artefact | "The Casoliva 2:1(b) member at C = 0.567 is itself an impact orbit of the three-body problem (periselene 1,415 km); it is excluded" |
| Table column "Peri at eps = 1 (cycler?)" with "yes" for Casoliva 2:1(b) | The catalogue row is `resonant_po` | "inside lunar SOI (member), row class resonant_po" |
| "C21 3D corridor (#682 / #438 spatial), 2/1 (1), 3.09895, 42,003 km, 1.1e5" | Planar, period 1 Tg, multiplier squared | "planar orbit, T = 1 Tg (a = 1, n = 1), C = 3.09895, periselene 42,003 km, parent multiplier about 330 per period; family not identified" |
| "Family identities found on the way: ... The #438 spatial C21 row and the #682 C21 3D corridor seed are one family" | Not re-checked here; and the 2/1 "member" is not on it | Keep the identity as the build's claim; remove the 2/1 member from it |
| "Positive control against Brown et al." and "Outcome: reproduced, in the sense in which the papers' claims transfer" | Not a control; passes for a wrong forcing | "Consistency check of the Melnikov code against the reversing symmetry and against continuation (no number of Brown et al. is reproduced; none is printed)" |
| "This is the alternation predicted by the sign of dMel/dtheta0" | Shown on the two Lyapunov controls only | Add "on the two Lyapunov members" |
| "(the lunar Laplace SOI used by `data/validate.py`)" | The number is not in `validate.py` | "(the lunar SOI radius quoted in the catalogue's `orbit_class` comments)" |
| "Gates (measured) ... eps = 1 RHS and STM RHS vs `core/bcr4bp` (separate code)" | Both codes share the sign; agreement between them is not evidence of a correct model | Add a gate against an inertial-frame integration |
| "The a = 3 (8:3) resonances have a Melnikov amplitude 0.3-1 percent ..." and every quoted amplitude | Project-model values; the corrected-sense amplitudes differ by 0.07 to 2000 times for a <= 2 | Recompute |
| "Pattern: For a = 1 and 2, every equivalent of every cycler-class family reaches physical Sun mass, except the three fold-backs" | Shipped model | After rerun: "every one, in the members tested" if the rerun confirms section 3 |
| "the Casoliva 1:2(d) ... periselene ~348,000 (no)" | Corrected sense: 594,530 km, a very different orbit | Recompute |

## 12. Methods, so the computations can be repeated

- **Equations of motion.** Three-body rotating-frame terms plus mu_S times (direct plus indirect)
  with the Sun at angle theta0 + w t; w = +omega_S is the project model, w = -omega_S the corrected
  sense. Variational equations carry the STM and dX/d eps. DOP853, rtol 1e-13, atol 1e-14.
- **Inertial-frame truth.** Barycentric non-rotating frame; Earth at -mu (cos t, sin t), Moon at
  (1 - mu)(cos t, sin t), Sun at a_S (cos(n_S t + phi0), sin(n_S t + phi0)) with n_S = 1 - omega_S;
  acceleration from Earth, Moon, Sun direct, minus mu_S R_S / a_S^3. Rotating state to inertial:
  rotate by t and add the frame velocity.
- **Multiple shooting.** Nodes at k P / N (N as in the build's records), Newton by least squares,
  tolerance 1e-11. Pseudo-arclength in (nodes, eps), step at most 0.015, accept on convergence within
  6 iterations and tangent cosine above 0.9, folds followed.
- **Forward start.** Nodes on the three-body orbit (one lap, extended periodically) at a zero of my
  own Melnikov quadrature (periodic trapezoid, 3000 points per TU), first-order predictor, Newton at
  eps = 1e-4.
- **Passes.** Every local minimum of the distance on a grid of 4000 points per TU, refined by bounded
  minimisation on the dense output.
- **Family walk.** Symmetric half-period shooting in the three-body problem, pseudo-arclength in
  (x0, [z0], vy0, T), step at most 0.03, tangent cosine guard, every crossing of a target period
  corrected at that period; where the row's crossing gave a stiff walk, the row's other
  perpendicular crossing was used.
- **Coherent-model homotopy.** Natural-parameter steps in lam (0.025 to 0.1, halved on failure),
  secant predictor, Newton on the cyclic shooting system built from `propagate_qbcp_pv(...,
  with_stm=True)` at 1e-12, nodes and period as in the bicircular orbit.
- Compute used: about 25 minutes of wall time on 4 processes in total. The #884 note was read
  as committed in 605d9082 and the two core modules as of 6bda1772, that is, before the coordinator
  registered #891 and #892 from the findings here.

## 13. Not done

- A coherent-model check in a validated coherent model (section 9 is conditional on #892).
- The a = 3 members and the degenerate C21 3D 5/2 member, in either sense.
- Two corrected-sense branches of the sub-surface Casoliva member (one stopped at eps = 0.991 by the
  time limit, one converged to 1.3e-11 without diagnostics).
- Identification of the planar period-Tg family and of the fold landing orbits' families; the C11b
  and C32 family identities claimed in the note.
- A second implementation of the corrected-sense survey. It rests on one code (mine), validated by
  the inertial-frame closure of every orbit and by reproducing the build in the project sense.
- No test suite was run; no repository file was changed.

## 14. Required before any of these orbits goes in the catalogue

1. Fix the Sun's sense in `core/bcr4bp.py` and add a gate against an inertial-frame integration (a
   gate that compares two rotating-frame codes sharing a sign cannot catch this). Then re-examine
   every earlier result that used the module (#292, #303, #304, #334, #412 and anything calling
   `bcr4bp_eom`); a correctness fix invalidates prior results in both directions.
2. Fix the alpha_2 / alpha_3 parities in `core/qbcp.py` and add a reversibility gate.
3. Re-run #884 in the corrected model with the build's own code and compare with section 3 here, so
   the corrected survey rests on two implementations.
4. Replace the sampled periselene with a refined minimum, and exclude members whose parent or forced
   orbit passes below 1,737.4 km (or below a stated minimum altitude).
5. Check each three-body member's minimal period and planarity before labelling it with a family,
   and detect bifurcation points in the family walk (a sign change of the out-of-plane or
   period-doubling stability index) instead of relying on the half-step guard.
6. Decide what a row would claim. A defensible row is "three-body member X has a periodic counterpart
   of period n synodic months in the bicircular problem", class and tier set by the existing rules,
   not "survives the Sun". Carry the parent row's `orbit_class`; do not promote a `resonant_po` or V0
   `quasi_cycler` family to cycler through a distant member.
7. Persist the same members into the repaired coherent model before any statement about the real
   Sun, and state the lunar-eccentricity caveat in the row.
8. Report one orbit per symmetry class (not per phase), with the largest multiplier only unless the
   subdominant ones pass the reciprocal-pair check.
9. Run the literature check (`search/literature_check.py`) and set `our_status` from it.
