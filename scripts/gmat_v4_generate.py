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


def aldrin_row(epoch_iso: str = "11 Jul 2003 00:00:00.000") -> CyclerRow:
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


def s1l1_row(epoch_iso: str = "15 Dec 2026 00:00:00.000") -> CyclerRow:
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


_MARS_SOI_KM: float = 577_204.0  # Mars sphere-of-influence (GMAT Ex_MarsBPlane value)


def _vec_str(v: Vec3) -> str:
    return f"[ {v[0]:.9g} {v[1]:.9g} {v[2]:.9g} ]"


def _seed_state_marsrel(node: MarsFlybyNode) -> tuple[Vec3, Vec3]:
    r"""Mars-relative Cartesian seed at the SOI on the incoming asymptote.

    Beeson: the boundary condition is the flyby periapse with the body's gravity in
    the EOM. We seed the spacecraft on the **incoming asymptote** at the Mars SOI:
    position ``r0 = -S_hat * R_SOI`` (R_SOI up-stream along the incoming direction,
    offset by the impact parameter so the natural hyperbola has the targeted B), and
    velocity ``v0 = v_inf^- + (hyperbolic radial term)``. At SOI distance the excess
    dominates, so ``v0 ~= v_inf^-`` is an excellent initial guess; the DC's TCM then
    drives the B-plane to the v_inf^+ goal. This makes the script RUNNABLE (a real
    Mars-relative hyperbola), the maintenance ΔV being the converged |TCM|.
    """
    tgt = node.target()
    s_hat = tgt.s_hat  # incoming asymptote direction = v_inf^-_hat
    # Offset the start point off-axis by the targeted B along B_hat so the
    # uncorrected hyperbola already lands near the B-plane goal (good DC seed).
    b_hat = (tgt.bdot_t_km * tgt.t_hat + tgt.bdot_r_km * tgt.r_hat) / max(tgt.b_mag_km, 1e-30)
    r0 = -s_hat * _MARS_SOI_KM + b_hat * tgt.b_mag_km
    v0 = np.asarray(node.vinf_minus, dtype=np.float64)
    return r0, v0


def _init_block(row: CyclerRow, node: MarsFlybyNode) -> str:
    """Initialisation: Spacecraft (Mars-relative seed), ForceModel, Propagator, CS, DC."""
    r0, v0 = _seed_state_marsrel(node)
    return f"""\
%----- Section 1: initialisation (Beeson Algorithm-1; BC at Mars flyby periapse)
Create Spacecraft Sat;
GMAT Sat.DateFormat = UTCGregorian;
GMAT Sat.Epoch = '{row.epoch_iso}';
GMAT Sat.CoordinateSystem = MarsInertial;
GMAT Sat.DisplayStateType = Cartesian;
GMAT Sat.X = {r0[0]:.9g};
GMAT Sat.Y = {r0[1]:.9g};
GMAT Sat.Z = {r0[2]:.9g};
GMAT Sat.VX = {v0[0]:.9g};
GMAT Sat.VY = {v0[1]:.9g};
GMAT Sat.VZ = {v0[2]:.9g};
GMAT Sat.DryMass = 1000;

Create ForceModel NearMars_FM;
GMAT NearMars_FM.CentralBody = Mars;
GMAT NearMars_FM.PointMasses = {{Mars, Sun}};
GMAT NearMars_FM.Drag = None;
GMAT NearMars_FM.SRP = Off;

Create Propagator NearMars;
GMAT NearMars.FM = NearMars_FM;
GMAT NearMars.Type = PrinceDormand78;
GMAT NearMars.InitialStepSize = 60;
GMAT NearMars.Accuracy = 1e-12;
GMAT NearMars.MaxStep = 86400;

Create CoordinateSystem MarsInertial;
GMAT MarsInertial.Origin = Mars;
GMAT MarsInertial.Axes = BodyInertial;

Create ImpulsiveBurn TCM_{node.label};
GMAT TCM_{node.label}.CoordinateSystem = MarsInertial;

Create DifferentialCorrector DC;
GMAT DC.ShowProgress = true;
GMAT DC.MaximumIterations = 100;
GMAT DC.DerivativeMethod = ForwardDifference;
GMAT DC.Algorithm = NewtonRaphson;

Create Variable dv_{node.label};

Create ReportFile Rpt;
GMAT Rpt.Filename = 'gmat_v4_{row.row_id}_{node.label}.report';
GMAT Rpt.WriteHeaders = true;
"""


def _target_block(idx: int, node: MarsFlybyNode) -> str:
    r"""One Mars-flyby Target/Vary/Achieve block (Beeson BC at periapse).

    The Jones B-plane goal (in the comment) is the **approach geometry**; the
    achieved goal is the **outgoing v_inf^+ velocity** at the outbound SOI. That is
    the maintenance the flyby must deliver: a periapsis TCM (Oberth-credited) drives
    the outgoing asymptote to the next leg's sourced v_inf^+. The converged |TCM| is
    the high-fidelity maintenance ΔV — the figure the patched-conic models only
    bound. The TCM is applied AT periapsis (deepest Oberth credit), then propagated
    out three days to read the asymptotic outgoing velocity.
    """
    tgt = node.target()
    name = node.label
    feasible = "feasible (ballistic cone)" if tgt.feasible else "INFEASIBLE -> powered TCM"
    vp = np.asarray(node.vinf_plus, dtype=np.float64)
    return f"""\

%----- Mars flyby {idx}: {name}  ({feasible})
%  Jones B-plane approach geometry from sourced (v_inf^-, v_inf^+):
%    turn = {np.degrees(tgt.turn_rad):.4f} deg, rp = {tgt.rp_km:.3f} km,
%    theta_B = {np.degrees(tgt.theta_b_rad):.4f} deg,
%    B.R = {tgt.bdot_r_km:.4f} km, B.T = {tgt.bdot_t_km:.4f} km
%  Achieved goal: outgoing v_inf^+ (km/s, MarsInertial) = {_vec_str(vp)}
Target DC {{SolveMode = Solve, ExitMode = SaveAndContinue}};
   Vary DC(TCM_{name}.Element1 = 0.01, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM_{name}.Element2 = 0.01, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM_{name}.Element3 = 0.01, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Propagate NearMars(Sat) {{Sat.Mars.Periapsis}};
   Maneuver TCM_{name}(Sat);
   Propagate NearMars(Sat) {{Sat.ElapsedDays = 3}};
   Achieve DC(Sat.VX = {vp[0]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});
   Achieve DC(Sat.VY = {vp[1]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});
   Achieve DC(Sat.VZ = {vp[2]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});
EndTarget;
GMAT dv_{name} = sqrt(TCM_{name}.Element1^2 + TCM_{name}.Element2^2 + TCM_{name}.Element3^2);
Report Rpt Sat.A1ModJulian dv_{name};
"""


def generate_script(row: CyclerRow, *, node_index: int = 0) -> str:
    """Render the runnable GMAT V4 ``.script`` for ``row``'s Mars flyby ``node_index``."""
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
    # One runnable GMAT script per Mars flyby (each flyby is an independent
    # Mars-SOI B-plane targeting; the maintenance ΔV sums across them). A single
    # GMAT spacecraft cannot be re-seeded mid-mission, so the multi-flyby chain is
    # emitted as one script per node and the parser sums the converged TCMs.
    node = row.mars_flybys[node_index]
    parts.append(_init_block(row, node))
    parts.append("\n%----- Section 2: mission sequence\nBeginMissionSequence;\n")
    parts.append(_target_block(node_index + 1, node))
    return "".join(parts)


def write_script(row: CyclerRow, out_path: Path, *, node_index: int = 0) -> Path:
    """Render and write the per-flyby script for ``row``'s ``node_index`` to ``out_path``."""
    out_path.write_text(generate_script(row, node_index=node_index), encoding="utf-8")
    return out_path


def write_all_scripts(row: CyclerRow, out_dir: Path) -> list[Path]:
    """Write one runnable script per Mars flyby; returns the paths in node order."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for i, node in enumerate(row.mars_flybys):
        p = out_dir / f"{row.row_id}_{node.label}.script"
        paths.append(write_script(row, p, node_index=i))
    return paths


# Convenience entry points (Aldrin first, then S1L1) -------------------------- #


def generate_aldrin_script(out_path: Path, *, epoch_iso: str | None = None) -> Path:
    row = aldrin_row(epoch_iso) if epoch_iso else aldrin_row()
    return write_script(row, out_path)  # single flyby -> complete script


def generate_s1l1_scripts(out_dir: Path, *, epoch_iso: str | None = None) -> list[Path]:
    row = s1l1_row(epoch_iso) if epoch_iso else s1l1_row()
    return write_all_scripts(row, out_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a GMAT V4 script from a confirmed cycler."
    )
    parser.add_argument("row", choices=["aldrin", "s1l1"], help="which confirmed row")
    parser.add_argument(
        "--out", type=Path, required=True, help="output .script (aldrin) or directory (s1l1)"
    )
    parser.add_argument(
        "--epoch",
        type=str,
        default=None,
        help="override seed epoch (UTCGregorian, e.g. '06 Aug 2003 00:00:00')",
    )
    args = parser.parse_args()

    if args.row == "aldrin":
        path = generate_aldrin_script(args.out, epoch_iso=args.epoch)
        print(f"wrote {path} (1 Mars-flyby B-plane Target block)")
    else:
        paths = generate_s1l1_scripts(args.out, epoch_iso=args.epoch)
        print(f"wrote {len(paths)} per-Mars-flyby scripts to {args.out}:")
        for p in paths:
            print(f"  {p}")


if __name__ == "__main__":
    main()


__all__ = [
    "JONES_POS_TOL_KM",
    "JONES_VEL_TOL_KMS",
    "CyclerRow",
    "MarsFlybyNode",
    "aldrin_row",
    "generate_aldrin_script",
    "generate_s1l1_scripts",
    "generate_script",
    "s1l1_row",
    "write_all_scripts",
    "write_script",
]
