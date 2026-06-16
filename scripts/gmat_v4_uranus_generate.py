"""Emit a parallel GMAT script for the V4-strict Uranian gauntlet (#335 Part B).

What this is
------------
A companion / cross-check artifact for the Python+SPICE V4-strict driver
(:mod:`cyclerfinder.data.validation.v4_uranus_strict`). The Python driver
is the headline V4-strict gauntlet because GMAT R2022a doesn't bundle the
custom CelestialBody definitions Uranus + its moons need for native
PointMasses-against-SPK use. Anyone willing to author those definitions
manually can drive this script via:

    env -u DISPLAY ~/GMAT/R2022a/bin/GmatConsole --run scripts/<output>.script

and verify the Python+SPICE drift series numerically.

Scope of the emitted script
---------------------------
* Spacecraft: SILVER Lambert IC at Umbriel encounter epoch (J2000-inertial,
  Uranus-centered), with the v_out from a Lambert solve to the next moon.
* ForceModel: Uranus central body with J2 (3.34343e-3, R_eq = 25559 km,
  Jacobson 2014) + PointMasses on the four other classical Uranian moons
  (Miranda, Ariel, Titania, Oberon) sourced from URA111. Umbriel and
  Oberon are perturbers throughout the leg (matching the v4_uranus_strict
  module).
* Propagator: PrinceDormand78 (GMAT's closest RK family to scipy DOP853).
* Mission: propagate n_cycles cycles, dump spacecraft state at each
  encounter epoch via Report.

Honest about what GMAT can / can't do
-------------------------------------
* GMAT R2022a doesn't natively recognise ``Umbriel`` / ``Oberon`` as
  CelestialBody objects against URA111. The script emits a Create
  CelestialBody for each Uranian moon naming the URA111 NAIF IDs;
  loading this requires the user to have:
    1. URA111 in $GMAT_ROOT/data/planetary_ephem/spk/uranian/ura111.bsp
       (Part A install).
    2. A SPICEKernel registration that adds the URA SPK to GMAT's
       kernel pool at run time (handled inline via SPICEKernel objects).
  The pole orientation / rotation parameters use the same Jacobson 2014
  values Python+SPICE uses.
* GMAT's ForceModel.PointMasses syntax requires the body to be a
  registered CelestialBody; this is why the CelestialBody declarations
  come BEFORE the ForceModel.
* If a user's GMAT install rejects one of the CelestialBody declarations
  (older R2020a, missing PCK overrides), Python+SPICE remains the
  reference path -- this script is for the cross-check, not the
  headline.

Usage
-----
    uv run python scripts/gmat_v4_uranus_generate.py \\
        --candidate-id repeated-moon-uranus-00000041 \\
        --sequence Umbriel,Oberon,Umbriel \\
        --leg-tofs-days 14.940560615336594,14.940560615336594 \\
        --launch-epoch-utc "15 Jan 2000 00:00:00.000" \\
        --n-cycles 3 \\
        --out /tmp/gmat_v4_uranus_strict.script
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import spiceypy as spice

from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v4_uranus import URANUS_J2, URANUS_R_EQ_KM
from cyclerfinder.data.validation.v4_uranus_strict import (
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
)
from cyclerfinder.search.discovery_campaign import DAY_S

# Uranus rotational parameters (Jacobson 2014, AJ 148:76; aligned with PCK).
URANUS_POLE_RA_DEG = 257.311
URANUS_POLE_DEC_DEG = -15.175

# GMAT-readable NAIF moon names + IDs.
_MOON_NAIF: dict[str, tuple[str, int]] = {
    "Miranda": ("MIRANDA", 705),
    "Ariel": ("ARIEL", 701),
    "Umbriel": ("UMBRIEL", 702),
    "Titania": ("TITANIA", 703),
    "Oberon": ("OBERON", 704),
}


@dataclass(frozen=True)
class GmatUranusScriptArgs:
    candidate_id: str
    sequence: tuple[str, ...]
    leg_tofs_days: tuple[float, ...]
    n_revs: tuple[int, ...]
    launch_epoch_utc_gmat: str
    """GMAT's expected format: '15 Jan 2000 00:00:00.000'."""
    launch_epoch_utc_iso: str
    """ISO-8601 used for SPICE str2et."""
    n_cycles: int
    perturber_moons: tuple[str, ...]
    spice_kernels: tuple[str, ...]
    j2: float
    r_eq_km: float


def _spice_moon_state(moon_name: str, et_seconds: float) -> tuple[np.ndarray, np.ndarray]:
    target = _MOON_NAIF[moon_name][0]
    state, _ = spice.spkezr(target, et_seconds, "J2000", "NONE", "URANUS")
    arr = np.asarray(state, dtype=np.float64)
    return arr[:3].copy(), arr[3:].copy()


def _lambert_initial_state(
    sequence: tuple[str, ...],
    leg_tofs_days: tuple[float, ...],
    n_revs: tuple[int, ...],
    et_launch: float,
    mu_primary: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve Lambert from moon-0 to moon-1 (the first leg) and return r0, v0."""
    r_a, v_a_moon = _spice_moon_state(sequence[0], et_launch)
    r_b, _ = _spice_moon_state(sequence[1], et_launch + leg_tofs_days[0] * DAY_S)
    nrev = max(0, n_revs[0])
    sols = _lambert(r_a, r_b, leg_tofs_days[0] * DAY_S, mu=mu_primary, max_revs=nrev)
    wanted = [s for s in sols if s.n_revs == n_revs[0]]
    if not wanted:
        raise RuntimeError("Lambert first-leg failed; cannot seed GMAT script")
    best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a_moon)))
    return r_a, best.v1


def _vec_str(v: np.ndarray) -> str:
    return f"{v[0]:.9g}, {v[1]:.9g}, {v[2]:.9g}"


def _celestial_body_block(moon: str) -> str:
    """Emit a CelestialBody declaration for ``moon`` referencing URA111 IDs."""
    _naif_name, naif_id = _MOON_NAIF[moon]
    mu = SATELLITES[moon].mu_km3_s2
    r = SATELLITES[moon].radius_eq_km
    return f"""\
Create CelestialBody {moon};
GMAT {moon}.NAIFId = {naif_id};
GMAT {moon}.OrbitSpiceKernelName = {{'{DEFAULT_URA_PATH}'}};
GMAT {moon}.EquatorialRadius = {r:.3f};
GMAT {moon}.Mu = {mu:.6g};
GMAT {moon}.CentralBody = 'Uranus';
GMAT {moon}.OrbitColor = SkyBlue;
GMAT {moon}.PosVelSource = 'SPICE';
GMAT {moon}.RotationDataSource = 'IAUSimplified';
"""


def generate_script(args: GmatUranusScriptArgs) -> str:
    """Render the V4-strict GMAT script for the configured run."""
    if args.sequence[0] != args.sequence[-1]:
        raise ValueError("sequence must be CLOSED (first == last)")
    spice.kclear()
    try:
        for k in args.spice_kernels:
            spice.furnsh(k)
        et_launch = float(spice.str2et(args.launch_epoch_utc_iso))
        mu_primary = PRIMARIES["Uranus"]
        r0, v0 = _lambert_initial_state(
            args.sequence,
            args.leg_tofs_days,
            args.n_revs,
            et_launch,
            mu_primary,
        )
    finally:
        spice.kclear()

    cycle_period_days = float(sum(args.leg_tofs_days))
    total_days = args.n_cycles * cycle_period_days

    # Header.
    header = f"""\
%% GMAT V4-strict script (auto-generated, #335 Part B cross-check)
%% Candidate:  {args.candidate_id}
%% Sequence:   {" -> ".join(args.sequence)}
%% Leg ToFs:   {args.leg_tofs_days} days
%% N revs:     {args.n_revs}
%% N cycles:   {args.n_cycles} (cycle period = {cycle_period_days:.4f} d)
%% Launch:     {args.launch_epoch_utc_gmat}
%% J2 = {args.j2:.5e}, R_eq = {args.r_eq_km:.1f} km (Jacobson 2014)
%% Perturbers: {args.perturber_moons}
%% URA111:     {DEFAULT_URA_PATH}
%%
%% NOTE: this script requires the URA111 satellite SPK to be loaded into
%% GMAT's SPICE kernel pool. If the CelestialBody declarations below
%% fail in your GMAT install, run the Python+SPICE driver instead --
%% it's the headline V4-strict path:
%%   uv run python scripts/run_335_silver_v41_gauntlet.py
"""

    # CelestialBody declarations for every Uranian moon (tour + perturber).
    body_blocks: list[str] = []
    for moon in dict.fromkeys((*args.sequence, *args.perturber_moons)):
        body_blocks.append(_celestial_body_block(moon))

    # Initialisation block.
    init = f"""\

%----- Section 1: initialisation
Create Spacecraft Sat;
GMAT Sat.DateFormat = UTCGregorian;
GMAT Sat.Epoch = '{args.launch_epoch_utc_gmat}';
GMAT Sat.CoordinateSystem = UranusJ2000Eq;
GMAT Sat.DisplayStateType = Cartesian;
GMAT Sat.X = {r0[0]:.9g};
GMAT Sat.Y = {r0[1]:.9g};
GMAT Sat.Z = {r0[2]:.9g};
GMAT Sat.VX = {v0[0]:.9g};
GMAT Sat.VY = {v0[1]:.9g};
GMAT Sat.VZ = {v0[2]:.9g};
GMAT Sat.DryMass = 1000;

Create CoordinateSystem UranusJ2000Eq;
GMAT UranusJ2000Eq.Origin = Uranus;
GMAT UranusJ2000Eq.Axes = MJ2000Eq;

Create ForceModel UranusSystemFM;
GMAT UranusSystemFM.CentralBody = Uranus;
GMAT UranusSystemFM.PointMasses = {{Uranus, {", ".join(args.perturber_moons)}}};
GMAT UranusSystemFM.PrimaryBodies = {{Uranus}};
GMAT UranusSystemFM.GravityField.Uranus.Degree = 2;
GMAT UranusSystemFM.GravityField.Uranus.Order = 0;
GMAT UranusSystemFM.GravityField.Uranus.PotentialFile = '';
GMAT UranusSystemFM.Drag = None;
GMAT UranusSystemFM.SRP = Off;

Create Propagator UranusProp;
GMAT UranusProp.FM = UranusSystemFM;
GMAT UranusProp.Type = PrinceDormand78;
GMAT UranusProp.InitialStepSize = 60;
GMAT UranusProp.Accuracy = 1e-10;
GMAT UranusProp.MaxStep = 86400;

Create ReportFile Rpt;
GMAT Rpt.Filename = 'gmat_v4_strict_{args.candidate_id}.report';
GMAT Rpt.WriteHeaders = true;
"""

    # Mission sequence: report state, propagate to each encounter, report again.
    mission_parts: list[str] = ["\n%----- Section 2: mission sequence\nBeginMissionSequence;\n"]
    mission_parts.append("Report Rpt Sat.A1ModJulian Sat.X Sat.Y Sat.Z Sat.VX Sat.VY Sat.VZ;\n")
    for cycle in range(args.n_cycles):
        for leg, tof in enumerate(args.leg_tofs_days):
            mission_parts.append(
                f"Propagate UranusProp(Sat) {{Sat.ElapsedDays = {tof:.6f}}};  "
                f"%% cycle {cycle} leg {leg} -> {args.sequence[leg + 1]}\n"
            )
            mission_parts.append(
                "Report Rpt Sat.A1ModJulian Sat.X Sat.Y Sat.Z Sat.VX Sat.VY Sat.VZ;\n"
            )
    mission_parts.append(f"\n%% total propagation: {total_days:.4f} days\n")

    return header + "".join(body_blocks) + init + "".join(mission_parts)


def write_script(args: GmatUranusScriptArgs, out_path: Path) -> Path:
    out_path.write_text(generate_script(args), encoding="utf-8")
    return out_path


def _gmat_epoch_from_iso(iso: str) -> str:
    """Convert ISO-8601 to GMAT UTCGregorian (e.g. '15 Jan 2000 00:00:00.000')."""
    from datetime import datetime

    if iso.endswith("Z"):
        iso = iso[:-1]
    dt = datetime.fromisoformat(iso)
    return dt.strftime("%d %b %Y %H:%M:%S.000")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Emit a GMAT V4-strict cross-check script for a Uranus-system cycler."
    )
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument(
        "--sequence",
        required=True,
        help="comma-separated CLOSED sequence, e.g. 'Umbriel,Oberon,Umbriel'",
    )
    parser.add_argument(
        "--leg-tofs-days",
        required=True,
        help="comma-separated per-leg ToF in days",
    )
    parser.add_argument(
        "--n-revs",
        default=None,
        help="comma-separated per-leg revolution count (default 1,1,...)",
    )
    parser.add_argument(
        "--launch-epoch-utc",
        required=True,
        help="ISO-8601 UTC launch epoch (e.g. '2000-01-15T00:00:00')",
    )
    parser.add_argument("--n-cycles", type=int, default=3)
    parser.add_argument(
        "--perturber-moons",
        default="Miranda,Ariel,Umbriel,Titania,Oberon",
        help="comma-separated perturber set (default all 5 classical)",
    )
    parser.add_argument("--out", type=Path, required=True)
    p = parser.parse_args()

    sequence = tuple(p.sequence.split(","))
    leg_tofs = tuple(float(x) for x in p.leg_tofs_days.split(","))
    n_legs = len(sequence) - 1
    n_revs = (
        tuple(int(x) for x in p.n_revs.split(",")) if p.n_revs else tuple(1 for _ in range(n_legs))
    )
    perturbers = tuple(p.perturber_moons.split(","))

    args = GmatUranusScriptArgs(
        candidate_id=p.candidate_id,
        sequence=sequence,
        leg_tofs_days=leg_tofs,
        n_revs=n_revs,
        launch_epoch_utc_gmat=_gmat_epoch_from_iso(p.launch_epoch_utc),
        launch_epoch_utc_iso=p.launch_epoch_utc,
        n_cycles=p.n_cycles,
        perturber_moons=perturbers,
        spice_kernels=(
            str(DEFAULT_LSK_PATH),
            str(DEFAULT_PCK_PATH),
            str(DEFAULT_URA_PATH),
        ),
        j2=URANUS_J2,
        r_eq_km=URANUS_R_EQ_KM,
    )

    path = write_script(args, p.out)
    print(f"wrote {path} ({path.stat().st_size} bytes)")
    print(f"To run in GMAT: env -u DISPLAY ~/GMAT/R2022a/bin/GmatConsole --run {path}")
    print(
        "NOTE: this is a cross-check artifact. The Python+SPICE driver "
        "(scripts/run_335_silver_v41_gauntlet.py) is the headline path."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
