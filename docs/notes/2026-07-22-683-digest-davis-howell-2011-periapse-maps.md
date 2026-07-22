# Digest — Davis & Howell 2011, periapsis Poincaré ("periapse") maps in the Saturnian system

**Date:** 2026-07-22 AET
**Task:** #683 — corpus digest + positive control for the new periapse-map
capability (`src/cyclerfinder/search/periapse_map.py`).
**Paper:** Diane Craig Davis, Kathleen C. Howell, "Trajectory evolution in the
multi-body problem with applications in the Saturnian system," *Acta
Astronautica* 69 (2011) 1038–1049,
doi:[10.1016/j.actaastro.2011.07.007](https://doi.org/10.1016/j.actaastro.2011.07.007).
**Provenance:** freely available from the authors' Purdue site
(`engineering.purdue.edu/people/kathleen.howell.1/Publications/Journals/2011_AA_DavHow.pdf`);
filed to the private corpus as
`davis-howell-2011-trajectory-evolution-multibody-saturnian-system-acta-astronautica-69-1038.pdf`
(md5 `9c3569d099ac676ac39191853adafa39`). Text-layer PDF, no OCR needed.

This is the Davis–Howell periapse-map lineage the #683 dispatch named (Davis &
Howell 2012 JGCD 35(1) is the sibling; Villac & Scheeres 2003 JGCD 26(2)
224–232 introduced the periapsis Poincaré map). This 2011 Acta paper is the one
that applies the map DIRECTLY to Saturn-Titan with explicit published Jacobi
values, so it is the strongest positive-control target for a Saturn-Titan build.

## The method (paper Secs. 2.2–2.4, 3)

**Model.** Planar CR3BP, primary P1 + secondary P2. Hill radius
`r_H = (mu/3)^(1/3)` (paper Eq. 1). Jacobi integral (paper Eq. 2)

    J = x^2 + y^2 + 2 mu / r + 2(1-mu)/d - v^2

with `r` = distance to the secondary P2, `d` = distance to the primary P1,
`v` = rotating-frame speed. **This is exactly this project's own
`core.cr3bp.jacobi_constant` convention** (verified term-by-term: our
`C = (x^2+y^2) + 2(1-mu)/r1 + 2 mu/r2 - v^2` with `r1=d`, `r2=r`).

**Periapsis = the section.** A state with `rdot = 0` (radial rate w.r.t. P2)
and `rddot > 0` (a local minimum of distance to P2) is a periapsis; `rddot < 0`
is an apoapsis. The `rdot = 0` locus is a contour around P2 (paper's
"zero-radial-acceleration" boundary between periapses and apoapses).

**Parametrisation.** For a chosen periapse position and fixed J, the velocity
MAGNITUDE is set by the Jacobi integral and the DIRECTION is perpendicular to
the P2-relative radius (the apse condition). The prograde choice fixes the sign.
Hence every in-contour periapse position defines a single planar prograde
trajectory (paper Sec. 2.2; the 3D case adds a velocity angle phi, not used
here).

**Map coordinates.** The periapse *position* relative to P2 in Hill-radius
units: `x_p = (x-(1-mu))/r_H`, `y_p = y/r_H` (the paper's Fig.-1 map axes).
Polar form: periapse angle `omega_r = atan2(y, x-(1-mu))`, periapse radius `r_p`.

**Initial-condition map (Figs. 3, 4).** Propagate each in-contour periapse IC
forward and classify the fate after N revolutions into four regions:

* IMPACT (black) — passes on or within the radius of P2;
* ESCAPE_L1 (blue) / ESCAPE_L2 (red) — an x-coordinate more than **0.01**
  nondimensional units beyond L1 / L2 respectively;
* CAPTURED (gray) — remains near P2 for the whole propagation.

The escape lobes are the stable-manifold tubes of the L1/L2 Lyapunov orbits
(paper Sec. 3, citing Koon et al.): a periapsis in an escape lobe escapes
before its next periapsis. The **reflected** map (across the x-axis, i.e. time
reversal) gives the ENTRANCE lobes — trajectories entering P2's vicinity
through the gateways (paper Sec. 6, Fig. 20, its Titan-capture design).

## Published numbers used as the positive control

* **Sun-Saturn search energy** `J1 = 3.0173046596239` (paper Figs. 3a/4a),
  stated `J < J_L2` (both necks open).
* **Saturn-Titan search energy** `J2 = 3.015311017945150` (paper Figs. 3b/4b/20),
  stated `J < J_L1` (both necks open).
* Cassini end-of-life examples at `J = 3.106` (Fig. 11/12), `J = 3.016`
  (Fig. 14). Phoebe SMA `214 R_S` (1 `R_S = 60,268 km`).

**Independent cross-check (this project's `jacobi_constant`).** I computed the
collinear L1/L2 points (roots of the on-axis effective-potential gradient) and
their Jacobi constants:

| system | C_L1 | C_L2 | published J | J below C_L2 by |
|---|---|---|---|---|
| Sun-Saturn | 3.017823881 | 3.017442757 | 3.017304660 | 1.38e-4 |
| Saturn-Titan | 3.015769581 | 3.015453948 | 3.015311018 | 1.43e-4 |

Both published search energies sit just below **C_L2** (so both L1 and L2 necks
are open), exactly matching the paper's stated regime — a clean, independent
confirmation that the paper's exact quoted Jacobi values are consistent with our
own dynamics. (`test_periapse_map.py::test_l_points_bracket_secondary_and_match_paper_regime`.)

## Paper-specific qualitative predictions reproduced (#683 Phase C)

Reproducing the Fig.-3/4 initial-condition maps at `J1`/`J2` (160×140 grid over
`x_p,y_p ∈ [-1.25,1.25]×[-1.05,1.05]`) confirmed:

1. **Capture-fraction collapse with more revolutions** (paper: "during the
   longer propagations, more trajectories escape or impact"): Saturn-Titan
   captured fraction falls 60% → 13% from 1 → 6 revolutions; escapes and
   impacts grow correspondingly.
2. **Titan impact ≫ Saturn impact** (paper Sec. 3: "in units of Hill radii,
   Titan's radius is much larger than that of Saturn and, consequently, many
   more initial states lead to impact in the Saturn-Titan system"): impact
   fraction 6.4% (Saturn-Titan, 1 rev) vs 0.1% (Sun-Saturn, 1 rev).
3. **Directional lobes**: ESCAPE_L1 lobes on the interior (x_p<0) side,
   ESCAPE_L2 on the exterior (x_p>0) side, with the L1 lobe the larger of the
   two (tidal-quadrant asymmetry, paper Sec. 2.3).

## Relationship to the sibling GAIO method (#664/#685)

Same underlying transport question (temporary capture / escape itineraries
near the smaller primary) as `#664`'s set-oriented transfer-operator pipeline,
but a **structurally different lens**: periapse maps give exact per-trajectory
lobe geometry at fixed energy with NO Monte-Carlo noise and NO shooting
correction; GAIO gives statistical measure/residence transport probabilities.
Ranked last in `#679`'s shortlist for exactly this overlap; run as a second,
independent data point per the coordinator's decision.

## Not adopted / out of scope

* The 3D velocity-angle `phi` parametrisation (paper Sec. 2.2, Fig. 2) and the
  long-term periapse-profile taxonomy (figure-8 / hourglass / lobe / arrowhead,
  Figs. 6–8) — the #683 planar build uses the position-space map only.
* The Cassini/Phoebe end-of-life ΔV design applications (Secs. 5–6) — mission
  design, not a cycler search.
