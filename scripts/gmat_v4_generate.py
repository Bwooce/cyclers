"""Generate a GMAT V4 high-fidelity ``.script`` from a confirmed cycler (#171).

Emits the buildable half of the canonical V4 lane: a GMAT mission script following
Beeson AAS 15-278 Algorithm-1's two-section skeleton (initialisation + mission
sequence with an embedded ``Target``/``Vary``/``Achieve`` block per Mars flyby).
Each Mars flyby carries a **Jones AAS 17-577 B-plane target** computed by
:mod:`cyclerfinder.verify.bplane` from that encounter's published
:math:`(v_\\infty^-, v_\\infty^+)` nodes: the boundary condition is propagated to
``Mars.Periapsis`` (Beeson: BCs moved from body-center to flyby periapse so the
flyby body's gravity enters the EOM), a maintenance TCM is varied, and the
``BdotR``/``BdotT`` B-plane goals are achieved.

Two confirmed rows flow through ONE entry point with different ``mode``:

* **Aldrin** (#134) — ``mode="powered-periodic"``: a single Mars flyby Target that
  recovers the documented per-synodic maintenance ΔV reference (2.9138 km/s, OUR
  value under external check).
* **S1L1** (``russell-ch4-4.991gG2``) — ``mode="flyby-station-keep"``: the 7-Mars-
  encounter chain unrolled (Beeson: loops "completely expanded in script form"),
  one B-plane Target per App-C Mars node read from
  :data:`cyclerfinder.search.s1l1_corrected.APPC_LEGS`.

GOLDEN / HONESTY. The only sourced inputs are each row's published :math:`v_\\infty`
nodes and the Jones continuity tolerance (1e-3 km / 1e-6 km/s). The reference ΔVs
(Aldrin 2.9138 km/s; S1L1 62 m/s) are OUR computed values — GMAT is the independent
external check, never an EXPECTED-from-source assertion. This script only TEMPLATES;
it neither installs nor runs GMAT (manual, out-of-CI step; see the run-book).

Usage::

    uv run python scripts/gmat_v4_generate.py aldrin --out /tmp/aldrin_v4.script
    uv run python scripts/gmat_v4_generate.py s1l1   --out /tmp/s1l1_v4.script
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from cyclerfinder.search.s1l1_corrected import APPC_EPOCH_DAYS, APPC_LEGS
from cyclerfinder.verify.bplane import BPlaneTarget, bplane_target_for

Vec3 = np.ndarray

# Jones AAS 17-577 §2.5 published ephemeris-continuity tolerance — the SOURCED
# convergence bar the GMAT Target/Achieve blocks are written against.
JONES_POS_TOL_KM: float = 1.0e-3
JONES_VEL_TOL_KMS: float = 1.0e-6


@dataclass(frozen=True)
class MarsFlybyNode:
    """One Mars flyby's sourced incoming/outgoing v_inf pair + epoch.

    ``vinf_minus`` / ``vinf_plus`` are body-centered-equatorial km/s vectors taken
    straight off the row's published nodes (App-C for S1L1; the documented Mars
    return-flyby geometry for Aldrin). ``epoch_days_j2000`` is the encounter epoch
    (days past J2000) used only to label/seed the GMAT propagation.
    """

    label: str
    vinf_minus: Vec3
    vinf_plus: Vec3
    epoch_days_j2000: float

    def target(self) -> BPlaneTarget:
        """The Jones B-plane goal for this Mars flyby (reuses verify/bplane)."""
        return bplane_target_for("M", self.vinf_minus, self.vinf_plus)


@dataclass(frozen=True)
class CyclerRow:
    """Generalisation contract: a confirmed-row descriptor the generator consumes.

    Aldrin and S1L1 differ only in this descriptor (id, sequence, the Mars-flyby
    node list, mode, reference ΔV) — the templating is shared (design §2).
    """

    row_id: str
    sequence: str
    mode: str  # "powered-periodic" (Aldrin) | "flyby-station-keep" (S1L1)
    epoch_iso: str
    seed_vinf_kms: Vec3  # the seed Earth-departure v_inf (the initial guess)
    mars_flybys: tuple[MarsFlybyNode, ...]
    reference_dv_kms: float | None  # OUR value under external check; None = convergence-only
    notes: str = ""
    extra_bodies: tuple[str, ...] = field(default_factory=tuple)


# --------------------------------------------------------------------------- #
# Row descriptors
# --------------------------------------------------------------------------- #


def aldrin_row(epoch_iso: str = "2003-07-11T00:00:00.000") -> CyclerRow:
    """The Aldrin (#134) single-powered-Mars-flyby descriptor.

    The Mars return-flyby geometry is the documented lap-0 powered member
    (``docs/notes/2026-06-07-aldrin-continuation-v3-evidence.md`` §2/§4): emerged
    V_inf at Mars 8.88 km/s, required turn 93.0 deg (> the 68.5 deg achievable cone,
    hence powered). We construct the (v_inf^-, v_inf^+) pair from that sourced
    magnitude + turn — the B-plane target the GMAT Mars flyby Target achieves. The
    per-synodic maintenance reference is 2.9138 km/s (OUR value, under GMAT's check).
    """
    vmag = 8.88  # documented Mars V_inf, km/s
    turn = np.radians(93.0)  # documented required return-flyby turn
    vinf_minus = np.array([vmag, 0.0, 0.0], dtype=np.float64)
    # Rotate in the x-y plane (axis = pole) by the documented turn for a concrete,
    # reproducible pair; |v_inf| is preserved so the turn IS the documented 93 deg.
    vinf_plus = np.array([vmag * np.cos(turn), vmag * np.sin(turn), 0.0], dtype=np.float64)
    return CyclerRow(
        row_id="aldrin-classic-em-k1-outbound",
        sequence="E-M-E",
        mode="powered-periodic",
        epoch_iso=epoch_iso,
        seed_vinf_kms=np.array([6.08, 0.0, 0.0], dtype=np.float64),  # documented Earth V_inf
        mars_flybys=(
            MarsFlybyNode(
                label="MarsReturn",
                vinf_minus=vinf_minus,
                vinf_plus=vinf_plus,
                epoch_days_j2000=1300.0,  # ~2003.5; label only
            ),
        ),
        reference_dv_kms=2.9138,
        notes="Aldrin lap-0 powered Mars return flyby; ref 2.9138 km/s is OUR value.",
    )


def s1l1_row(epoch_iso: str = "2026-12-15T00:00:00.000") -> CyclerRow:
    """The S1L1 (``russell-ch4-4.991gG2``) per-Mars-flyby station-keep descriptor.

    Reads the App-C nodes from :data:`APPC_LEGS` (the sourced Russell 2004 Appendix-C
    block). Each Mars node (body ``"M"``) gets a B-plane Target whose
    (v_inf^-, v_inf^+) is (the previous leg's outgoing v_inf as the incoming
    asymptote, the Mars node's own v_inf as the outgoing) — the consecutive App-C
    nodes, NEVER a self-computed attribute. Seed = the App-C leg-2 Earth departure
    (2026-12-15, v_inf=(-2.278, 5.322, 0.574) km/s). Reference = 62 m/s patched-conic
    horizon-TCM over 7 cycles (OUR value, under GMAT's check).
    """
    by_no = {leg[0]: leg for leg in APPC_LEGS}
    flybys: list[MarsFlybyNode] = []
    for leg_no, body, ts, vinf in APPC_LEGS:
        if body != "M" or vinf is None:
            continue
        # Incoming asymptote = the previous (Earth) leg's outgoing App-C v_inf;
        # outgoing = this Mars node's App-C v_inf. Consecutive sourced nodes.
        prev = by_no[leg_no - 1]
        prev_vinf = prev[3]
        assert prev_vinf is not None, "Mars node always has a preceding leg with v_inf"
        flybys.append(
            MarsFlybyNode(
                label=f"Mars{leg_no}",
                vinf_minus=np.array(prev_vinf, dtype=np.float64),
                vinf_plus=np.array(vinf, dtype=np.float64),
                epoch_days_j2000=APPC_EPOCH_DAYS + ts,
            )
        )

    leg2 = by_no[2]
    seed_vinf = leg2[3]
    assert seed_vinf is not None
    return CyclerRow(
        row_id="russell-ch4-4.991gG2",
        sequence="E-g(E-E)-E-G(E-M-E)-E",
        mode="flyby-station-keep",
        epoch_iso=epoch_iso,
        seed_vinf_kms=np.array(seed_vinf, dtype=np.float64),
        mars_flybys=tuple(flybys),
        reference_dv_kms=0.062,  # 62 m/s patched-conic over 7 cycles; OUR value
        notes="S1L1 per-Mars-flyby B-plane station-keep; ref 62 m/s is OUR value.",
    )


# --------------------------------------------------------------------------- #
# Script templating (Beeson Algorithm-1 two-section skeleton)
# --------------------------------------------------------------------------- #


def _vec_str(v: Vec3) -> str:
    return f"[ {v[0]:.9g} {v[1]:.9g} {v[2]:.9g} ]"


def _force_model_block(extra_bodies: tuple[str, ...]) -> str:
    """Initialisation: Spacecraft, ForceModel (Sun + planets), Propagator, B-plane CS."""
    point_masses = ["Sun", "Earth", "Luna", "Mars", "Jupiter", *extra_bodies]
    # de-dup while preserving order
    pm = list(dict.fromkeys(point_masses))
    pm_line = ", ".join(pm)
    return f"""\
%% ---- Section 1: initialisation (Beeson Algorithm-1) ----
Create Spacecraft Sat;
GMAT Sat.DateFormat = UTCGregorian;
GMAT Sat.CoordinateSystem = EarthMJ2000Eq;

Create ForceModel FM;
GMAT FM.CentralBody = Sun;
GMAT FM.PointMasses = {{{pm_line}}};
GMAT FM.Drag = None;
GMAT FM.SRP = Off;

Create Propagator Prop;
GMAT Prop.FM = FM;
GMAT Prop.Type = RungeKutta89;
GMAT Prop.InitialStepSize = 60;
GMAT Prop.Accuracy = 1e-12;

Create CoordinateSystem MarsBPlane;
GMAT MarsBPlane.Origin = Mars;
GMAT MarsBPlane.Axes = BodyInertial;
"""


def _initial_guess_block(row: CyclerRow) -> str:
    """The 'Provide an Initial Guess' slot — seed departure state + epoch (Beeson)."""
    return f"""\

%% ---- Provide an Initial Guess (the medium-fidelity seed) ----
%% Row: {row.row_id}  sequence: {row.sequence}  mode: {row.mode}
GMAT Sat.Epoch = '{row.epoch_iso}';
%% Seed Earth-departure hyperbolic-excess (km/s), body-centered:
GMAT Sat.SeedVinf = {_vec_str(row.seed_vinf_kms)};
"""


def _tcm_and_target_block(idx: int, node: MarsFlybyNode) -> str:
    """One Mars-flyby Target/Vary/Achieve block (Jones B-plane goal at periapse)."""
    tgt = node.target()
    name = node.label
    feasible = "feasible (ballistic cone)" if tgt.feasible else "INFEASIBLE -> powered TCM"
    return f"""\

%% ---- Mars flyby {idx}: {name}  ({feasible}) ----
%% Jones B-plane goal from sourced (v_inf^-, v_inf^+):
%%   turn = {np.degrees(tgt.turn_rad):.4f} deg, rp = {tgt.rp_km:.3f} km,
%%   theta_B = {np.degrees(tgt.theta_b_rad):.4f} deg, |B| = {tgt.b_mag_km:.3f} km
Create ImpulsiveBurn TCM_{name};
GMAT TCM_{name}.CoordinateSystem = Local;
GMAT TCM_{name}.Origin = Sun;
GMAT TCM_{name}.Axes = VNB;

Target FlybyTCM_{name} {{SolveMode = Solve, ExitMode = SaveAndContinue}};
   Vary    TCM_{name}.Element1 = 0.0;
   Vary    TCM_{name}.Element2 = 0.0;
   Vary    TCM_{name}.Element3 = 0.0;
   Maneuver TCM_{name}(Sat);
   Propagate Prop(Sat) {{Sat.Mars.Periapsis}};
   Achieve  Sat.MarsBPlane.BdotR = {tgt.bdot_r_km:.9g} {{Tolerance = {JONES_POS_TOL_KM:g}}};
   Achieve  Sat.MarsBPlane.BdotT = {tgt.bdot_t_km:.9g} {{Tolerance = {JONES_POS_TOL_KM:g}}};
EndTarget;
"""


def _report_block(row: CyclerRow) -> str:
    tcm_cols = " ".join(f"TCM_{n.label}.Magnitude" for n in row.mars_flybys)
    return f"""\

%% ---- Report (the parser reads Maneuver magnitudes + convergence) ----
Create ReportFile Rpt;
GMAT Rpt.Filename = 'gmat_v4_{row.row_id}.report';
GMAT Rpt.WriteHeaders = true;
Report Rpt Sat.A1ModJulian {tcm_cols};
"""


def generate_script(row: CyclerRow) -> str:
    """Render the full GMAT V4 ``.script`` text for ``row``."""
    parts: list[str] = []
    ref_str = (
        f"{row.reference_dv_kms:.4f} km/s"
        if row.reference_dv_kms is not None
        else "convergence-only"
    )
    parts.append(
        f"%% GMAT V4 high-fidelity script (auto-generated, #171)\n"
        f"%% Row: {row.row_id}  mode: {row.mode}\n"
        f"%% Reference maintenance dV (OUR value, GMAT is the external check): {ref_str}\n"
        f"%% Jones continuity tolerance (SOURCED): "
        f"{JONES_POS_TOL_KM:g} km / {JONES_VEL_TOL_KMS:g} km/s\n"
        f"%% {row.notes}\n"
    )
    parts.append(_force_model_block(row.extra_bodies))
    parts.append(_initial_guess_block(row))
    parts.append("\n%% ---- Section 2: mission sequence ----\nBeginMissionSequence;\n")
    for i, node in enumerate(row.mars_flybys, start=1):
        parts.append(_tcm_and_target_block(i, node))
    parts.append(_report_block(row))
    return "".join(parts)


def write_script(row: CyclerRow, out_path: Path) -> Path:
    """Render and write the script for ``row`` to ``out_path``."""
    out_path.write_text(generate_script(row), encoding="utf-8")
    return out_path


# Convenience entry points (Aldrin first, then S1L1) -------------------------- #


def generate_aldrin_script(out_path: Path, *, epoch_iso: str | None = None) -> Path:
    row = aldrin_row(epoch_iso) if epoch_iso else aldrin_row()
    return write_script(row, out_path)


def generate_s1l1_script(out_path: Path, *, epoch_iso: str | None = None) -> Path:
    row = s1l1_row(epoch_iso) if epoch_iso else s1l1_row()
    return write_script(row, out_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a GMAT V4 script from a confirmed cycler."
    )
    parser.add_argument("row", choices=["aldrin", "s1l1"], help="which confirmed row")
    parser.add_argument("--out", type=Path, required=True, help="output .script path")
    parser.add_argument(
        "--epoch", type=str, default=None, help="override seed epoch (UTCGregorian ISO)"
    )
    args = parser.parse_args()

    if args.row == "aldrin":
        path = generate_aldrin_script(args.out, epoch_iso=args.epoch)
    else:
        path = generate_s1l1_script(args.out, epoch_iso=args.epoch)
    n_flybys = (
        len(
            generate_script(aldrin_row() if args.row == "aldrin" else s1l1_row()).split(
                "Target FlybyTCM_"
            )
        )
        - 1
    )
    print(f"wrote {path} ({n_flybys} Mars-flyby B-plane Target block(s))")


if __name__ == "__main__":
    main()


__all__ = [
    "JONES_POS_TOL_KM",
    "JONES_VEL_TOL_KMS",
    "CyclerRow",
    "MarsFlybyNode",
    "aldrin_row",
    "generate_aldrin_script",
    "generate_s1l1_script",
    "generate_script",
    "s1l1_row",
    "write_script",
]
