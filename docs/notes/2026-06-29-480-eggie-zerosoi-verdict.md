# #480 zero-SOI patched-conic verdict — EGGIE does not close at 0.70 m/s in our coplanar model either; the likely gap is the 3D B-plane flyby DOF

**Date:** 2026-06-29
**Task:** #480 follow-up — reproduce EGGIE in the paper's EXACT zero-radius-SOI
patched-conic model (where the 0.70 m/s is defined), after the continuous-gravity
n-body corrector hit a ~0.1 km/s wall (`...-eggie-stage4-subarc-verdict.md`).
**Shipped:** `nbody/jovian.py::flyby_maneuver_dv` — the paper's Eqs 3-5 tangential-
periapsis maneuver ΔV (`722feda`), unit-tested.

## One-line verdict — characterized NEGATIVE (with a specific, testable cause)

In our ideal **circular-coplanar** model, an exhaustive departure-phase × ToF grid in
the zero-SOI patched-conic formulation (flyby bend ballistic via Eq 4 r_p in the
25-70,000 km window; ΔV = Eq 5 magnitude-mismatch burn) finds a **feasible, cycle-
closed** EGGIE only at **~1.0-1.4 km/s total ΔV, and only with the encounter V∞ pushed
far off Table 4** (Ganymede ~9.2, Europa ~10.6 vs targets 7.07, 9.12). At the
resonant-conic V∞ (which ARE on Table 4 to ~0.25 km/s) the cycle does **not** close
feasibly. So ballistic closure and the Table-4 V∞ are **incompatible in our coplanar
reconstruction** — the paper's 0.70 m/s ballistic EGGIE is not reproduced here.

## What was ruled out

- **Seed/basin:** the resonant conic (Stage 1) puts all three V∞ within ~0.25 km/s of
  Table 4 (Io 8.44/Europa 9.12/Ganymede 6.81 at e≈0.62). The seed is right.
- **Ideal-SMA fidelity (the agent's flagged lead):** recomputing with the paper's
  Io period 1.75 d (a_Io=418,625 km) instead of the real 1.77 d changes the V∞ match
  negligibly (err 0.255 vs 0.269 km/s; Ganymede 6.82 either way). **SMA is not the gap.**
- **FD/STM/epoch/sub-arc correctors:** all hit the continuous-gravity ~0.1 km/s wall
  (the four prior verdicts). The zero-SOI model removes continuous gravity and STILL
  doesn't close at Table-4 V∞ — so the wall is not a continuous-gravity artifact either.

## The likely cause — the coplanar restriction removes the 3D B-plane flyby DOF

The paper's flyby (Eqs 6-7) is a **3-D B-plane maneuver**: the asymptotes define a
plane and θ_B (Eq 7) is a free orientation of the bend about the v∞ direction. A
gravity assist can rotate v∞ **out of the moons' orbital plane**. Our reconstruction
holds everything **coplanar** (circular-coplanar moons, in-plane Lambert legs, in-plane
bends), which deletes exactly that rotational degree of freedom at every flyby. With
only in-plane bending, the five encounters cannot simultaneously (a) sit on the
resonant Table-4-V∞ conic, (b) each be a feasible equal-magnitude ballistic flyby, and
(c) close the cycle — so the optimizer must trade V∞ off-target to find any feasible
closed cycle, landing at ~1 km/s. The 0.70 m/s ballistic EGGIE plausibly **requires the
out-of-plane B-plane freedom** the paper's model has and ours lacks. (Secondary
candidate: the paper's cycler-closure may be defined in a synodic/rotating frame or as
"shape-repeats", not the strict inertial state-return our periodicity residual imposes.)

## Status of #480 (summary across the whole follow-up)

- **Stage 1 resonant-conic generator: SUCCESS** — all 3 V∞ on Table 4 (`535d2fb`).
- **Corrector program: infrastructure SHIPPED + reusable** — analytic state+STM Jovian
  co-integrator (flyby-parity-gated, `d619c44`), block-bidiagonal Jacobian (~40× FD,
  `9f98bb1`), `jacobian="stm"` + sub-arc multiple shooting (`193569c`/`428ca31`),
  zero-SOI Eq 3-5 flyby maneuver (`722feda`).
- **EGGIE ballistic reproduction: characterized NEGATIVE** in BOTH continuous-gravity
  n-body (~0.1 km/s wall, robust to 4 correctors) AND zero-SOI patched-conic coplanar
  (~1 km/s feasible-closed, V∞ off-target). Seed/basin and SMA ruled out; the leading
  remaining cause is the missing **3-D B-plane flyby DOF** (coplanar restriction).
- **No catalogue change.** Golden (`tests/verify/test_ieg_reproduction_golden.py`)
  stays skipped.

## Next lever (one, fresh build) + honest prior

A **3-D flyby reconstruction**: let the Lambert legs + B-plane flyby orientation (θ_B,
Eq 7) leave the moons' orbital plane, then re-run the zero-SOI close. This is the lever
most likely to actually reproduce the 0.70 m/s EGGIE (it restores the DOF the paper
uses), but it is a fresh build (3-D conic geometry + B-plane targeting), not a re-run.
Also worth a cheap check first: confirm the paper's cycler-closure frame (synodic vs
inertial state-return) — a mis-imposed periodicity constraint would also force the
off-target trade. Scratch search drivers removed.
