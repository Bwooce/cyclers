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
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from cyclerfinder.core.constants import PLANETS, SAFE_PERIHELION_KM
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

    def target(self, body: str = "M") -> BPlaneTarget:
        """The Jones B-plane goal for this flyby (reuses verify/bplane)."""
        return bplane_target_for(body, self.vinf_minus, self.vinf_plus)


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
    flyby_body: str = "M"  # one-letter code: the body each flyby is around (M=Mars, E=Earth)
    # Targeting mode for the Achieve block:
    #  "vinf-vector" — Achieve the outgoing v_inf VECTOR (VX/VY/VZ): forces both
    #     direction AND magnitude. Correct for an equal-|v_inf| turn (Aldrin).
    #  "bplane-position" — Achieve the B-plane aim point (BdotR, BdotT) only: the
    #     OPERATIONAL maintenance target (next-encounter position), letting the
    #     outgoing |v_inf| be whatever the geometry gives (#176 S1L1 reconciliation).
    target_mode: str = "vinf-vector"
    # If True, clamp a subsurface Jones-Eq.2 seed periapsis up to the safe radius so
    # the DC pays the honest turn deficit at a real flyby (#176 Aldrin Earth leg).
    clamp_safe_rp: bool = False


# --------------------------------------------------------------------------- #
# Row descriptors
# --------------------------------------------------------------------------- #


def aldrin_row(epoch_iso: str = "11 Jul 2003 00:00:00.000") -> CyclerRow:
    """The Aldrin (#134) single-powered-Mars-flyby descriptor.

    **DEPRECATED for the maintenance reconciliation (#176).** This targets the
    *Mars* return flyby (V_inf 8.88 km/s, turn 93 deg), but Aldrin's documented
    per-synodic maintenance cost is the *Earth* return-flyby turn deficit (93 deg
    required vs ~68 deg achievable at a safe Earth flyby — see
    ``docs/notes/2026-06-07-aldrin-continuation-v3-evidence.md`` §2/§4 and the #148
    primer-refine finding: "the 2.9138 km/s is *entirely* an Earth-flyby turn-deficit
    charge"). The Mars flyby is NOT the constrained maneuver; its converged TCM
    (0.175 km/s, #174) is the WRONG leg and is not comparable to the 2.9138 / 1.9336
    references. Use :func:`aldrin_earth_row` for the maintenance reconciliation.

    Kept verbatim for reproducibility of the #174 (wrong-leg) result.
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
        notes="Aldrin lap-0 powered Mars return flyby (WRONG LEG for maintenance; see #176).",
    )


def aldrin_earth_row(epoch_iso: str = "11 Jul 2003 00:00:00.000") -> CyclerRow:
    """The Aldrin (#134) EARTH-return-flyby maintenance descriptor (#176, the right leg).

    Aldrin's per-synodic maintenance cost is the **Earth** return flyby's turn
    deficit, NOT the Mars flyby. The in-family DE440 solve recovers V_inf at Earth
    = 6.86 km/s with a required return-flyby turn of 93.0 deg, against an achievable
    cone of ~68 deg at a safe Earth flyby (``2026-06-07-aldrin-continuation-v3-
    evidence.md`` §2: turn req 93.0 / turn max 68.5 -> powered;
    ``2026-06-07-oberth-flyby-recost.md``: V_inf_E 6.86, deficit 24.5 deg). We build
    the equal-|V_inf| (v_inf^-, v_inf^+) pair from that magnitude + turn — the same
    pair the asymptote model charges 2.9702 km/s and the Oberth-periapsis model
    charges 1.9715 km/s for (reproduced numerically in #176). GMAT burns at Earth
    periapsis (deepest Oberth credit), so it SHOULD land near the Oberth ~1.97 km/s.

    Reference for the band check: the Oberth-periapsis figure 1.9336 km/s (#151,
    integration-confirmed #154) — the operationally defensible maneuver (burn deep
    in the well). The asymptote 2.9138 km/s is the conservative no-Oberth upper
    bound, recorded but not the band target.
    """
    vmag = 6.86  # in-family DE440 Earth V_inf, km/s (#151/#154)
    turn = np.radians(93.0)  # documented required Earth return-flyby turn
    vinf_minus = np.array([vmag, 0.0, 0.0], dtype=np.float64)
    vinf_plus = np.array([vmag * np.cos(turn), vmag * np.sin(turn), 0.0], dtype=np.float64)
    return CyclerRow(
        row_id="aldrin-classic-em-k1-outbound-EARTH",
        sequence="E-M-E",
        mode="powered-periodic",
        epoch_iso=epoch_iso,
        seed_vinf_kms=np.array([vmag, 0.0, 0.0], dtype=np.float64),
        mars_flybys=(
            # Reuse the MarsFlybyNode container for an EARTH flyby: the body is set
            # to Earth in the templating via the row's flyby_body. The node carries
            # the sourced equal-|V_inf| turn pair.
            MarsFlybyNode(
                label="EarthReturn",
                vinf_minus=vinf_minus,
                vinf_plus=vinf_plus,
                epoch_days_j2000=1300.0,  # ~2003.5; label only
            ),
        ),
        reference_dv_kms=1.9336,
        notes="Aldrin EARTH return-flyby maintenance (#176 right leg); ref Oberth 1.9336 km/s.",
        flyby_body="E",
        clamp_safe_rp=True,  # the 93 deg turn needs a powered burn at safe periapsis
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

# Per-flyby-body GMAT names + SOI radius (km). The generator was Mars-only (#174);
# #176 adds Earth so the Aldrin maintenance leg (the EARTH return flyby) can be
# targeted. SOI radii are standard patched-conic values.
_BODY_GMAT_NAME: dict[str, str] = {"M": "Mars", "E": "Earth"}
_BODY_SOI_KM: dict[str, float] = {"M": 577_204.0, "E": 924_000.0}


def _vec_str(v: Vec3) -> str:
    return f"[ {v[0]:.9g} {v[1]:.9g} {v[2]:.9g} ]"


def _b_for_safe_rp(vinf_mag: float, mu: float, rp_safe: float) -> float:
    r"""Impact parameter (km) whose natural hyperbola has periapsis exactly ``rp_safe``.

    :math:`|B| = r_p\sqrt{1 + 2\mu/(r_p v_\infty^2)}`. Used to clamp a too-deep
    (subsurface) Jones-Eq.2 seed up to a physically-realisable safe flyby (#176): the
    Jones root for a large required turn (Aldrin's 93 deg at Earth) is subsurface, so
    seeding at it lets the DC "pass through the planet" and recover a spuriously cheap
    TCM. Clamping the seed to the safe periapsis (natural turn = the cone max) makes
    the DC genuinely PAY the turn deficit at a real flyby.
    """
    return float(rp_safe * np.sqrt(1.0 + 2.0 * mu / (rp_safe * vinf_mag * vinf_mag)))


def _seed_state_marsrel(
    node: MarsFlybyNode, body: str = "M", *, clamp_safe_rp: bool = False
) -> tuple[Vec3, Vec3]:
    r"""Body-relative Cartesian seed at the SOI on the incoming asymptote.

    Beeson: the boundary condition is the flyby periapse with the body's gravity in
    the EOM. We seed the spacecraft on the **incoming asymptote** at the body SOI:
    position ``r0 = -S_hat * R_SOI`` (R_SOI up-stream along the incoming direction,
    offset by the impact parameter so the natural hyperbola has the targeted B), and
    velocity ``v0 = v_inf^-``. At SOI distance the excess dominates, so ``v0`` is an
    excellent initial guess; the DC's TCM then drives the goal. This makes the script
    RUNNABLE (a real body-relative hyperbola), the maintenance ΔV being the converged
    |TCM|.

    ``clamp_safe_rp`` (#176): if the Jones-Eq.2 periapsis for the required turn is
    BELOW the safe periapsis (a subsurface root — the geometry needs a powered turn),
    seed instead at the impact parameter for the SAFE periapsis. The natural hyperbola
    then turns by the ballistic cone max; the DC's TCM pays the honest turn deficit.
    """
    tgt = node.target(body)
    s_hat = tgt.s_hat  # incoming asymptote direction = v_inf^-_hat
    b_hat = (tgt.bdot_t_km * tgt.t_hat + tgt.bdot_r_km * tgt.r_hat) / max(tgt.b_mag_km, 1e-30)
    b_mag = tgt.b_mag_km
    if clamp_safe_rp:
        mu = PLANETS[body].mu_km3_s2
        rp_safe = SAFE_PERIHELION_KM[body]
        vinf_mag = float(np.linalg.norm(node.vinf_minus))
        if tgt.rp_km < rp_safe:  # subsurface / too-deep Jones root -> clamp to safe
            b_mag = _b_for_safe_rp(vinf_mag, mu, rp_safe)
    r0 = -s_hat * _BODY_SOI_KM[body] + b_hat * b_mag
    v0 = np.asarray(node.vinf_minus, dtype=np.float64)
    return r0, v0


_TWO_BURN_MODES: frozenset[str] = frozenset({"oberth-two-burn", "oberth-optimize"})


def _second_burn_decl(row: CyclerRow, node: MarsFlybyNode, cs: str) -> str:
    """Declare the second ImpulsiveBurn (only for the two-burn / optimize modes)."""
    if row.target_mode not in _TWO_BURN_MODES:
        return ""
    return (
        f"Create ImpulsiveBurn TCM2_{node.label};\n"
        f"GMAT TCM2_{node.label}.CoordinateSystem = {cs};\n"
    )


def _extra_dv_vars(row: CyclerRow, node: MarsFlybyNode) -> str:
    """Declare the per-burn dv variables (only for the oberth-optimize mode)."""
    if row.target_mode != "oberth-optimize":
        return ""
    return f"Create Variable dv1_{node.label} dv2_{node.label};\n"


def _optimizer_decl(row: CyclerRow) -> str:
    """Declare the Yukon NLP optimizer (only for the oberth-optimize mode)."""
    if row.target_mode != "oberth-optimize":
        return ""
    return "\nCreate Yukon NLP;\nGMAT NLP.ShowProgress = true;\nGMAT NLP.MaximumIterations = 200;\n"


def _init_block(row: CyclerRow, node: MarsFlybyNode) -> str:
    """Initialisation: Spacecraft (body-relative seed), ForceModel, Propagator, CS, DC."""
    body = row.flyby_body
    gmat_body = _BODY_GMAT_NAME[body]
    cs = f"{gmat_body}Inertial"
    prop = f"Near{gmat_body}"
    r0, v0 = _seed_state_marsrel(node, body, clamp_safe_rp=row.clamp_safe_rp)
    return f"""\
%----- Section 1: initialisation (Beeson Algorithm-1; BC at {gmat_body} flyby periapse)
Create Spacecraft Sat;
GMAT Sat.DateFormat = UTCGregorian;
GMAT Sat.Epoch = '{row.epoch_iso}';
GMAT Sat.CoordinateSystem = {cs};
GMAT Sat.DisplayStateType = Cartesian;
GMAT Sat.X = {r0[0]:.9g};
GMAT Sat.Y = {r0[1]:.9g};
GMAT Sat.Z = {r0[2]:.9g};
GMAT Sat.VX = {v0[0]:.9g};
GMAT Sat.VY = {v0[1]:.9g};
GMAT Sat.VZ = {v0[2]:.9g};
GMAT Sat.DryMass = 1000;

Create ForceModel {prop}_FM;
GMAT {prop}_FM.CentralBody = {gmat_body};
GMAT {prop}_FM.PointMasses = {{{gmat_body}, Sun}};
GMAT {prop}_FM.Drag = None;
GMAT {prop}_FM.SRP = Off;

Create Propagator {prop};
GMAT {prop}.FM = {prop}_FM;
GMAT {prop}.Type = PrinceDormand78;
GMAT {prop}.InitialStepSize = 60;
GMAT {prop}.Accuracy = 1e-12;
GMAT {prop}.MaxStep = 86400;

Create CoordinateSystem {cs};
GMAT {cs}.Origin = {gmat_body};
GMAT {cs}.Axes = BodyInertial;

Create ImpulsiveBurn TCM_{node.label};
GMAT TCM_{node.label}.CoordinateSystem = {cs};
{_second_burn_decl(row, node, cs)}
Create DifferentialCorrector DC;
GMAT DC.ShowProgress = true;
GMAT DC.MaximumIterations = 100;
GMAT DC.DerivativeMethod = ForwardDifference;
GMAT DC.Algorithm = NewtonRaphson;
{_optimizer_decl(row)}
Create Variable dv_{node.label};
{_extra_dv_vars(row, node)}
Create ReportFile Rpt;
GMAT Rpt.Filename = 'gmat_v4_{row.row_id}_{node.label}.report';
GMAT Rpt.WriteHeaders = true;
"""


def _optimize_block(
    idx: int,
    node: MarsFlybyNode,
    row: CyclerRow,
    tgt: BPlaneTarget,
    vp: Vec3,
    gmat_body: str,
    cs: str,
    prop: str,
    feasible: str,
) -> str:
    """An Optimize/Minimize block: minimize total two-burn ΔV s.t. v_inf^+ achieved.

    The deep-well-optimal Oberth maneuver (#176). Two tangential burns (slow before
    periapsis to widen the cone, restore |v_inf| after) with their summed magnitude
    minimized by the Yukon NLP optimizer subject to the outgoing v_inf^+ vector
    equality constraints. This finds the cheapest maneuver, unlike a DC which lands
    on any feasible point.
    """
    name = node.label
    return f"""\

%----- {gmat_body} flyby {idx}: {name}  ({feasible})  [mode: oberth-optimize / Yukon NLP]
%  Jones B-plane approach geometry from sourced (v_inf^-, v_inf^+):
%    turn = {np.degrees(tgt.turn_rad):.4f} deg, rp = {tgt.rp_km:.3f} km
%  MINIMIZE total |TCM1|+|TCM2| s.t. outgoing v_inf^+ (km/s, {cs}) = {_vec_str(vp)}
Optimize NLP {{SolveMode = Solve, ExitMode = SaveAndContinue}};
   Vary NLP(TCM_{name}.Element1 = -0.5, {{Perturbation = 1e-5, MaxStep = 0.5}});
   Vary NLP(TCM_{name}.Element2 = -0.5, {{Perturbation = 1e-5, MaxStep = 0.5}});
   Vary NLP(TCM_{name}.Element3 = 0.0, {{Perturbation = 1e-5, MaxStep = 0.5}});
   Vary NLP(TCM2_{name}.Element1 = 0.5, {{Perturbation = 1e-5, MaxStep = 0.5}});
   Vary NLP(TCM2_{name}.Element2 = 0.5, {{Perturbation = 1e-5, MaxStep = 0.5}});
   Vary NLP(TCM2_{name}.Element3 = 0.0, {{Perturbation = 1e-5, MaxStep = 0.5}});
   Propagate {prop}(Sat) {{Sat.{gmat_body}.Periapsis, StopTolerance = 1e-6}};
   Maneuver TCM_{name}(Sat);
   Propagate {prop}(Sat) {{Sat.{gmat_body}.Periapsis, StopTolerance = 1e-6}};
   Maneuver TCM2_{name}(Sat);
   Propagate {prop}(Sat) {{Sat.ElapsedDays = 3}};
   GMAT dv1_{name} = sqrt(TCM_{name}.Element1^2 + TCM_{name}.Element2^2 + TCM_{name}.Element3^2);
   GMAT dv2_{name} = sqrt(TCM2_{name}.Element1^2 + TCM2_{name}.Element2^2 + TCM2_{name}.Element3^2);
   GMAT dv_{name} = dv1_{name} + dv2_{name};
   Minimize NLP(dv_{name});
   NonlinearConstraint NLP(Sat.VX = {vp[0]:.9g});
   NonlinearConstraint NLP(Sat.VY = {vp[1]:.9g});
   NonlinearConstraint NLP(Sat.VZ = {vp[2]:.9g});
EndOptimize;
Report Rpt Sat.A1ModJulian dv_{name};
"""


def _target_block(idx: int, node: MarsFlybyNode, row: CyclerRow) -> str:
    r"""One flyby Target/Vary/Achieve block (Beeson BC at periapse).

    Two targeting modes (``row.target_mode``):

    * ``"vinf-vector"`` (Aldrin) — Achieve the **outgoing v_inf^+ velocity vector**
      (VX/VY/VZ) at the outbound SOI: forces direction AND magnitude. For an equal-
      |v_inf| turn this IS the maintenance maneuver (the turn the cone cannot supply
      ballistically). The TCM is applied AT periapsis (deepest Oberth credit), then
      propagated out three days to read the asymptotic outgoing velocity.
    * ``"bplane-position"`` (S1L1 #176) — Achieve the **B-plane aim point**
      (``BdotR``, ``BdotT``) only — the OPERATIONAL maintenance target (the next
      encounter's POSITION). The outgoing |v_inf| is left to be whatever the flyby
      geometry gives; the maneuver does NOT pay to force a v_inf magnitude. This is
      the honest operational cost: aim at the next encounter, ride the natural v_inf.

    The converged |TCM| is the high-fidelity maintenance ΔV.
    """
    body = row.flyby_body
    gmat_body = _BODY_GMAT_NAME[body]
    cs = f"{gmat_body}Inertial"
    prop = f"Near{gmat_body}"
    tgt = node.target(body)
    name = node.label
    feasible = "feasible (ballistic cone)" if tgt.feasible else "INFEASIBLE -> powered TCM"
    vp = np.asarray(node.vinf_plus, dtype=np.float64)

    vary = f"""\
   Vary DC(TCM_{name}.Element1 = 0.01, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM_{name}.Element2 = 0.01, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM_{name}.Element3 = 0.01, {{Perturbation = 1e-4, MaxStep = 1.0}});"""

    if row.target_mode == "oberth-optimize":
        # OBERTH-OPTIMAL via the Yukon NLP optimizer (#176): two tangential periapsis
        # burns whose total magnitude is MINIMIZED subject to achieving the outgoing
        # v_inf^+ vector. A DifferentialCorrector only SOLVES (any feasible point); an
        # optimizer finds the deep-well-cheapest maneuver — the strategy the analytic
        # 1.9336 km/s (#151/#154) prices.
        return _optimize_block(idx, node, row, tgt, vp, gmat_body, cs, prop, feasible)
    if row.target_mode == "oberth-two-burn":
        # OBERTH widen-cone-then-restore (#176, the deep-well-optimal turn): a first
        # tangential burn slows the spacecraft BEFORE periapsis (a slower hyperbola
        # bends MORE, so the widened cone delivers the full required turn), then a
        # second tangential burn AFTER periapsis restores |v_inf|. The DC varies both
        # 3-vector burns to Achieve the outgoing v_inf^+; the summed |TCM1|+|TCM2| is
        # the Oberth maintenance. This is the strategy the analytic 1.9336 km/s
        # (#151, integration-confirmed #154) prices.
        vary = f"""\
   Vary DC(TCM_{name}.Element1 = -0.5, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM_{name}.Element2 = -0.5, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM_{name}.Element3 = 0.0, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM2_{name}.Element1 = 0.5, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM2_{name}.Element2 = 0.5, {{Perturbation = 1e-4, MaxStep = 1.0}});
   Vary DC(TCM2_{name}.Element3 = 0.0, {{Perturbation = 1e-4, MaxStep = 1.0}});"""
        achieve = f"""\
   Propagate {prop}(Sat) {{Sat.{gmat_body}.Periapsis, StopTolerance = 1e-6}};
   Maneuver TCM_{name}(Sat);
   Propagate {prop}(Sat) {{Sat.{gmat_body}.Periapsis, StopTolerance = 1e-6}};
   Maneuver TCM2_{name}(Sat);
   Propagate {prop}(Sat) {{Sat.ElapsedDays = 3}};
   Achieve DC(Sat.VX = {vp[0]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});
   Achieve DC(Sat.VY = {vp[1]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});
   Achieve DC(Sat.VZ = {vp[2]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});"""
        goal_comment = (
            f"%  OBERTH two-burn (widen-cone then restore) -> v_inf^+ (km/s, {cs}) = {_vec_str(vp)}"
        )
    elif row.target_mode == "bplane-position":
        # OPERATIONAL: aim the B-plane position only. BdotR/BdotT are read at the
        # incoming SOI BEFORE periapsis (the aim point is an approach quantity); the
        # periapsis TCM nulls the targeting offset. |v_inf| rides free.
        achieve = f"""\
   Propagate {prop}(Sat) {{Sat.{gmat_body}.Periapsis}};
   Maneuver TCM_{name}(Sat);
   Propagate {prop}(Sat) {{Sat.ElapsedDays = 3}};
   Achieve DC(Sat.{cs}.BdotR = {tgt.bdot_r_km:.9g}, {{Tolerance = {JONES_POS_TOL_KM:g}}});
   Achieve DC(Sat.{cs}.BdotT = {tgt.bdot_t_km:.9g}, {{Tolerance = {JONES_POS_TOL_KM:g}}});"""
        goal_comment = (
            f"%  OPERATIONAL B-plane-POSITION target (|v_inf| rides free):\n"
            f"%    B.R = {tgt.bdot_r_km:.4f} km, B.T = {tgt.bdot_t_km:.4f} km"
        )
    else:
        achieve = f"""\
   Propagate {prop}(Sat) {{Sat.{gmat_body}.Periapsis}};
   Maneuver TCM_{name}(Sat);
   Propagate {prop}(Sat) {{Sat.ElapsedDays = 3}};
   Achieve DC(Sat.VX = {vp[0]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});
   Achieve DC(Sat.VY = {vp[1]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});
   Achieve DC(Sat.VZ = {vp[2]:.9g}, {{Tolerance = {JONES_VEL_TOL_KMS:g}}});"""
        goal_comment = f"%  Achieved goal: outgoing v_inf^+ (km/s, {cs}) = {_vec_str(vp)}"

    return f"""\

%----- {gmat_body} flyby {idx}: {name}  ({feasible})  [mode: {row.target_mode}]
%  Jones B-plane approach geometry from sourced (v_inf^-, v_inf^+):
%    turn = {np.degrees(tgt.turn_rad):.4f} deg, rp = {tgt.rp_km:.3f} km,
%    theta_B = {np.degrees(tgt.theta_b_rad):.4f} deg,
%    B.R = {tgt.bdot_r_km:.4f} km, B.T = {tgt.bdot_t_km:.4f} km
{goal_comment}
Target DC {{SolveMode = Solve, ExitMode = SaveAndContinue}};
{vary}
{achieve}
EndTarget;
{_dv_report(row, name)}"""


def _dv_report(row: CyclerRow, name: str) -> str:
    """The dv accumulation + Report line(s); sums both burns in two-burn mode."""
    mag1 = f"sqrt(TCM_{name}.Element1^2 + TCM_{name}.Element2^2 + TCM_{name}.Element3^2)"
    if row.target_mode == "oberth-two-burn":
        mag2 = f"sqrt(TCM2_{name}.Element1^2 + TCM2_{name}.Element2^2 + TCM2_{name}.Element3^2)"
        return f"GMAT dv_{name} = {mag1} + {mag2};\nReport Rpt Sat.A1ModJulian dv_{name};\n"
    return f"GMAT dv_{name} = {mag1};\nReport Rpt Sat.A1ModJulian dv_{name};\n"


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
    parts.append(_target_block(node_index + 1, node, row))
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


def generate_aldrin_earth_script(out_path: Path, *, epoch_iso: str | None = None) -> Path:
    """The #176 CORRECTED Aldrin maintenance: the EARTH return flyby (the right leg)."""
    row = aldrin_earth_row(epoch_iso) if epoch_iso else aldrin_earth_row()
    return write_script(row, out_path)


def generate_aldrin_earth_oberth_script(out_path: Path, *, epoch_iso: str | None = None) -> Path:
    """The #176 Aldrin EARTH maintenance via the OBERTH two-burn widen-restore strategy."""
    base = aldrin_earth_row(epoch_iso) if epoch_iso else aldrin_earth_row()
    row = replace(
        base,
        row_id=base.row_id + "-OBERTH",
        target_mode="oberth-two-burn",
        notes="Aldrin EARTH Oberth two-burn widen-restore (#176); ref 1.9336 km/s.",
    )
    return write_script(row, out_path)


def generate_aldrin_earth_oberth_opt_script(
    out_path: Path, *, epoch_iso: str | None = None
) -> Path:
    """The #176 Aldrin EARTH maintenance via OBERTH two-burn MINIMIZED (Yukon NLP)."""
    base = aldrin_earth_row(epoch_iso) if epoch_iso else aldrin_earth_row()
    row = replace(
        base,
        row_id=base.row_id + "-OBERTHOPT",
        target_mode="oberth-optimize",
        notes="Aldrin EARTH Oberth two-burn MINIMIZED (#176, Yukon NLP); ref 1.9336 km/s.",
    )
    return write_script(row, out_path)


def generate_s1l1_scripts(out_dir: Path, *, epoch_iso: str | None = None) -> list[Path]:
    row = s1l1_row(epoch_iso) if epoch_iso else s1l1_row()
    return write_all_scripts(row, out_dir)


def generate_s1l1_operational_scripts(out_dir: Path, *, epoch_iso: str | None = None) -> list[Path]:
    """The #176 OPERATIONAL S1L1 maintenance: B-plane-POSITION targeting (|v_inf| free)."""
    base = s1l1_row(epoch_iso) if epoch_iso else s1l1_row()
    row = replace(
        base,
        row_id=base.row_id + "-OPER",
        target_mode="bplane-position",
        notes="S1L1 OPERATIONAL B-plane-position station-keep (#176); |v_inf| rides free.",
    )
    return write_all_scripts(row, out_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a GMAT V4 script from a confirmed cycler."
    )
    parser.add_argument(
        "row",
        choices=[
            "aldrin",
            "aldrin-earth",
            "aldrin-earth-oberth",
            "aldrin-earth-oberth-opt",
            "s1l1",
            "s1l1-oper",
        ],
        help=(
            "which confirmed row: aldrin (Mars leg, #174 wrong leg), "
            "aldrin-earth (#176 corrected maintenance leg, single periapsis burn), "
            "aldrin-earth-oberth (#176 Earth leg, Oberth two-burn DC), "
            "aldrin-earth-oberth-opt (#176 Earth leg, Oberth two-burn Yukon-minimized), "
            "s1l1 (#174 vinf-vector), s1l1-oper (#176 operational B-plane-position)"
        ),
    )
    parser.add_argument(
        "--out", type=Path, required=True, help="output .script (aldrin*) or directory (s1l1*)"
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
        print(f"wrote {path} (1 Mars-flyby B-plane Target block; #174 WRONG LEG)")
    elif args.row == "aldrin-earth":
        path = generate_aldrin_earth_script(args.out, epoch_iso=args.epoch)
        print(f"wrote {path} (1 EARTH-flyby maintenance Target block; #176 right leg)")
    elif args.row == "aldrin-earth-oberth":
        path = generate_aldrin_earth_oberth_script(args.out, epoch_iso=args.epoch)
        print(f"wrote {path} (Aldrin EARTH Oberth two-burn DC; #176)")
    elif args.row == "aldrin-earth-oberth-opt":
        path = generate_aldrin_earth_oberth_opt_script(args.out, epoch_iso=args.epoch)
        print(f"wrote {path} (Aldrin EARTH Oberth two-burn Yukon-minimized; #176)")
    elif args.row == "s1l1-oper":
        paths = generate_s1l1_operational_scripts(args.out, epoch_iso=args.epoch)
        print(f"wrote {len(paths)} per-Mars-flyby OPERATIONAL (B-plane-position) scripts:")
        for p in paths:
            print(f"  {p}")
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
    "aldrin_earth_row",
    "aldrin_row",
    "generate_aldrin_earth_oberth_opt_script",
    "generate_aldrin_earth_oberth_script",
    "generate_aldrin_earth_script",
    "generate_aldrin_script",
    "generate_s1l1_operational_scripts",
    "generate_s1l1_scripts",
    "generate_script",
    "s1l1_row",
    "write_all_scripts",
    "write_script",
]
