"""Energy-MOVING continuation of GMOS quasi-periodic 2-tori (#466 / #333 follow-on).

The #333 pseudo-arclength walker (``qp_tori_arclength``) follows the natural family
tangent, which on the #290 smoke seed points almost entirely along the rotation
number ``rho`` -- the family is near-iso-energetic (``C_J`` spans ~6e-7) and never
descends toward the #320 SILVER Bracket-2 region at ``C_J ~ 3.032``. The #290<->#320
connection hypothesis stayed open.

This module forces energy motion by **parent-family-driven continuation**: the #296
(1,1) 3D vertical parent family (``data/family_296_3d_em_11.jsonl``) is itself an
energy continuation -- ``C_J`` is strictly monotone in ``step_index``, spanning
2.920...3.148, and a complex Floquet pair sits ON the unit circle at every member from
the #333 seed (step 112, C_J~3.12785) all the way down to the #320 SILVER Bracket-2
(step 8, C_J~3.03196). So the energy-moving QP walk steps the PARENT orbit down the
family member-by-member and re-converges the QP-torus at each parent energy with
``correct_qp_torus`` (seeded from that member's monodromy Neimark-Sacker eigenvector).
The parent member's ``C_J`` is the torus energy, so the walk MOVES in energy by
construction, and the connection question becomes concrete: does a valid irrational
torus survive all the way from step 112 to step 8?

A direct energy-PIN on the torus mean state ``c_0`` was prototyped and rejected: it
fights the GMOS invariance rows for the same 6 DOF (the mean state IS the parent
orbit, whose energy is set by the family) and is compute-infeasible (>75s/step, fails
even at dC=1e-4). See ``docs/notes/2026-06-26-466-energy-moving-qp-walk-design.md``.

Report-only: NO catalogue writeback, NO novelty claims. Any candidate quasi-cycler
surfacing here is flagged for human gauntlet review after ``search/literature_check``,
NEVER self-admitted. The GMOS method is Olikara-Scheeres 2010 / Olikara 2016; the tori
returned are OUR computation (no published Earth-Moon QP-torus family table is digested
as a sourced golden -- standing acquisition gap, see the design note).
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.genome.qp_tori import (
    QPTorus,
    correct_qp_torus,
    is_practically_irrational,
)
from cyclerfinder.genome.qp_tori_arclength import (
    QPTorusFamilyMember,
    ResonanceCrossing,
    ResonanceFlag,
    _nearest_rational,
)

# Reuse the #333 family-member record (energy walk produces the same kind of member:
# a converged torus tagged with its C_J / rho / closure).
EnergyWalkMember = QPTorusFamilyMember

_ROOT = Path(__file__).resolve().parents[3]
PARENT_FAMILY_FILE = _ROOT / "data" / "family_296_3d_em_11.jsonl"


@dataclass(frozen=True)
class ParentMember:
    step_index: int
    state_nd: NDArray[np.float64]
    period: float
    jacobi: float
    floquet: NDArray[np.complex128]


def load_parent_family(path: Path | None = None) -> list[ParentMember]:
    """Load the #296 (1,1) 3D parent family, sorted by step_index (ascending C_J)."""
    p = path or PARENT_FAMILY_FILE
    members: list[ParentMember] = []
    with p.open() as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("type") == "header":
                continue
            fr = np.asarray(d["floquet_real"], dtype=np.float64)
            fi = np.asarray(d["floquet_imag"], dtype=np.float64)
            members.append(
                ParentMember(
                    step_index=int(d["step_index"]),
                    state_nd=np.asarray(d["state_nd"], dtype=np.float64),
                    period=float(d["T_TU"]),
                    jacobi=float(d["jacobi_constant"]),
                    floquet=(fr + 1j * fi).astype(np.complex128),
                )
            )
    members.sort(key=lambda m: m.step_index)
    return members


def _nearest_unit_complex_pair(
    floquet: NDArray[np.complex128], *, angle_min: float = 0.05
) -> tuple[complex, complex] | None:
    """Return the GENUINE rotation Floquet pair on the unit circle, or None.

    A 3D vertical-family member carries TWO complex unit-circle pairs: (a) the
    trivial central/energy pair near +1 (eigen-angle ~0, an artefact of the
    conserved-energy + autonomous-period directions), and (b) the genuine
    Neimark-Sacker ROTATION pair at a non-trivial angle (the QP-torus center,
    e.g. ~1.66 rad at the k=4 bracket). Selecting by "nearest |lam|=1" picks (a)
    arbitrarily and the seeded torus then collapses to the parent periodic orbit
    (rho->0). We instead require a non-trivial eigen-angle ``|phi| > angle_min``
    (excludes the near-+1 and near--1 trivial directions) and, among those on the
    unit circle, pick the one with the SMALLEST non-trivial angle (the leading
    rotation mode -- highest k -- whose torus is thinnest / closest to the family).
    """
    cands = []
    for ev in floquet:
        if abs(ev.imag) <= 1e-6:
            continue
        if abs(abs(ev) - 1.0) > 0.2:
            continue
        phi = math.atan2(ev.imag, ev.real)
        if abs(phi) < angle_min or abs(abs(phi) - math.pi) < angle_min:
            continue  # trivial near-+1 / near--1 direction, not a rotation
        cands.append(ev)
    if not cands:
        return None
    lam = min(cands, key=lambda ev: abs(math.atan2(ev.imag, ev.real)))
    # orient to upper-half-plane so the rho seed (atan2) is positive
    if lam.imag < 0:
        lam = lam.conjugate()
    return complex(lam), complex(lam.conjugate())


def _k_from_pair(lam: complex) -> int:
    """Integer k such that the multiplier sits near a primitive k-th root of unity."""
    phi = math.atan2(lam.imag, lam.real)
    if abs(phi) < 1e-9:
        return 4
    k = round(2 * math.pi / abs(phi))
    return max(k, 2)


def _converge_member_torus(
    parent: ParentMember,
    system: cr3bp.CR3BPSystem,
    *,
    n_trans: int,
    amplitude: float,
    tol: float,
    max_iter: int,
    independent_tol: float,
    notes: str,
) -> QPTorus | None:
    """Re-converge the QP-torus at a parent-family member, or None if no QP center."""
    pair = _nearest_unit_complex_pair(parent.floquet)
    if pair is None:
        return None
    k = _k_from_pair(pair[0])
    try:
        torus = correct_qp_torus(
            system,
            parent.state_nd,
            parent.period,
            pair,
            k=k,
            n_long=16,
            n_trans=n_trans,
            initial_torus_amplitude=amplitude,
            tol=tol,
            max_iter=max_iter,
            independent_tol=independent_tol,
            notes=notes,
        )
    except (RuntimeError, ValueError):
        return None
    return torus


def _member_from_torus(
    torus: QPTorus,
    parent: ParentMember,
    *,
    arclength_s: float,
    resonance_max_denominator: int = 12,
    resonance_tol: float = 1e-4,
) -> EnergyWalkMember:
    rho = torus.rho
    t_strob = torus.t_strob
    omega_long = 2 * math.pi / t_strob
    omega_trans = rho / t_strob
    freq_ratio = omega_trans / omega_long if omega_long != 0.0 else float("nan")
    irrational = is_practically_irrational(
        freq_ratio, max_denominator=resonance_max_denominator, tol=resonance_tol
    )
    p, q, dist = _nearest_rational(freq_ratio, resonance_max_denominator)
    near = ResonanceFlag(p, q, dist) if dist < resonance_tol else None
    return EnergyWalkMember(
        torus=torus,
        jacobi=float(parent.jacobi),
        arclength_s=float(arclength_s),
        tangent=np.zeros(0),
        rho=float(rho),
        freq_ratio=float(freq_ratio),
        is_practically_irrational=bool(irrational),
        near_resonance=near,
        fold_index=None,
        residual_norm=float(torus.invariance_residual),
        extras={
            "independent_residual": float(torus.independent_closure_residual),
            "parent_step": float(parent.step_index),
        },
    )


@dataclass(frozen=True)
class EnergyWalkFamily:
    members: list[EnergyWalkMember]
    resonance_crossings: list[ResonanceCrossing]
    terminated_reason: str
    seed_torus_id: str
    stop_cj: float = field(default=float("nan"))


def _tail_ratio(torus: QPTorus) -> float:
    c = torus.fourier_coeffs
    c1 = float(np.linalg.norm(c[1, :]))
    c_n = float(np.linalg.norm(c[torus.n_modes, :]))
    return c_n / c1 if c1 > 0 else float("inf")


def walk_energy(
    seed_torus: QPTorus,
    *,
    dc: float = 1e-3,
    max_steps: int = 200,
    direction: Literal["down", "up"] = "down",
    cj_target: float | None = None,
    target_step: int | None = None,
    step_stride: int = 1,
    seed_step: int | None = None,
    corrector_tol: float = 1e-8,
    corrector_max_iter: int = 40,
    independent_tol: float = 1e-3,
    amplitude: float = 5e-4,
    phase_pin_idx: int | None = None,
    dc_min: float = 1e-5,
    structural_residual_gate: float = 1e-3,
    independent_gate: float = 1e-4,
    resonance_max_denominator: int = 12,
    resonance_tol: float = 1e-4,
    mode_truncation_guard: float = 0.5,
    on_step: Callable[[EnergyWalkMember], None] | None = None,
    parent_family: list[ParentMember] | None = None,
) -> EnergyWalkFamily:
    """Parent-family-driven energy continuation of a GMOS QP-torus.

    Steps the #296 parent orbit DOWN (or up) its family ladder from the seed's parent
    member, re-converging the QP-torus at each parent energy. ``cj_target`` (Jacobi)
    or ``target_step`` (parent step_index) optionally caps the descent; else runs
    ``max_steps`` parent members at stride ``step_stride``. ``direction='down'``
    marches toward LOWER C_J (toward the #320 SILVER region at C_J~3.032 / step 8).

    The seed's parent step is auto-detected by matching the seed torus's mean-state
    energy to the parent ladder (override with ``seed_step``). A parent member at which
    no irrational torus converges terminates the walk with ``corrector_fail`` (a family
    boundary in energy -- report ``stop_cj``). Report-only -- NO catalogue writeback.
    """
    system = seed_torus.system
    n_modes = seed_torus.n_modes
    fam = parent_family if parent_family is not None else load_parent_family()
    by_step = {m.step_index: m for m in fam}

    cj_seed = jacobi_constant(np.real(seed_torus.fourier_coeffs[0, :]), system.mu)
    if seed_step is None:
        seed_step = min(fam, key=lambda m: abs(m.jacobi - cj_seed)).step_index

    sgn = -1 if direction == "down" else 1
    ordered_steps = sorted(by_step.keys(), reverse=(direction == "down"))
    # steps strictly past the seed in the chosen direction
    walk_steps = [s for s in ordered_steps if (s - seed_step) * sgn > 0]
    walk_steps = walk_steps[:: max(step_stride, 1)][:max_steps]

    crossings: list[ResonanceCrossing] = []
    seed_member = (
        _member_from_torus(
            seed_torus,
            by_step[seed_step],
            arclength_s=0.0,
            resonance_max_denominator=resonance_max_denominator,
            resonance_tol=resonance_tol,
        )
        if seed_step in by_step
        else _member_from_torus(
            seed_torus,
            ParentMember(
                seed_step,
                np.real(seed_torus.fourier_coeffs[0, :]).astype(np.float64),
                seed_torus.t_strob,
                cj_seed,
                np.zeros(6, dtype=np.complex128),
            ),
            arclength_s=0.0,
            resonance_max_denominator=resonance_max_denominator,
            resonance_tol=resonance_tol,
        )
    )
    members: list[EnergyWalkMember] = [seed_member]
    if on_step is not None:
        on_step(seed_member)

    reason = "max_steps"
    stop_cj = cj_seed
    prev_cj = cj_seed
    for st in walk_steps:
        parent = by_step[st]
        # Stop once the next member would lie BEYOND the target in the walk
        # direction (descending: parent energy/step strictly below the target).
        if cj_target is not None and (parent.jacobi - cj_target) * sgn > 0:
            reason = "reached_target"
            break
        if target_step is not None and (st - target_step) * sgn > 0:
            reason = "reached_target"
            break
        torus = _converge_member_torus(
            parent,
            system,
            n_trans=n_modes,
            amplitude=amplitude,
            tol=corrector_tol,
            max_iter=corrector_max_iter,
            independent_tol=independent_tol,
            notes="466_energy_walk_member",
        )
        if torus is None or not np.isfinite(torus.invariance_residual):
            reason = "corrector_fail"
            stop_cj = prev_cj
            break
        member = _member_from_torus(
            torus,
            parent,
            arclength_s=abs(parent.jacobi - cj_seed),
            resonance_max_denominator=resonance_max_denominator,
            resonance_tol=resonance_tol,
        )
        # Structural gate = the corrector-INDEPENDENT off-grid topology check (a
        # corrector residual can be small while the off-grid check fails -- so the
        # off-grid check is the genuine torus gate, per v1_qp). The GMOS Fourier-norm
        # at lower-energy members sits at the ~1.4e-5 FD-Jacobian conditioning floor
        # (just above V1's permissive 1e-5 Fourier gate; the documented Phase-1 limit)
        # while the off-grid residual stays ~3e-6 -- the member IS a genuine torus.
        # The per-member V1_qp verdict is REPORTED (not used to terminate); a walk
        # terminator on the marginal Fourier floor would falsely reject valid tori.
        if not (
            np.isfinite(member.residual_norm)
            and member.residual_norm < structural_residual_gate
            and member.extras["independent_residual"] < independent_gate
        ):
            reason = "corrector_fail"
            stop_cj = prev_cj
            break
        if _tail_ratio(torus) > mode_truncation_guard:
            members.append(member)
            if on_step is not None:
                on_step(member)
            reason = "mode_truncation_breach"
            stop_cj = member.jacobi
            break
        if member.near_resonance is not None and not member.is_practically_irrational:
            crossings.append(
                ResonanceCrossing(
                    len(members),
                    member.near_resonance.p,
                    member.near_resonance.q,
                    member.freq_ratio,
                )
            )
            members.append(member)
            if on_step is not None:
                on_step(member)
            reason = "resonance_lock"
            stop_cj = member.jacobi
            break
        members.append(member)
        if on_step is not None:
            on_step(member)
        prev_cj = member.jacobi
        stop_cj = member.jacobi
    return EnergyWalkFamily(members, crossings, reason, seed_torus.notes, stop_cj)
