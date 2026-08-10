"""N=5 CRNBP-torus-to-real-ephemeris consistency check (`#726`).

`#714`->`#724` built and independently CONFIRMED (twice) a narrow novelty
claim for a quasi-periodic invariant torus in the Jupiter-Io-Europa-Ganymede
CRNBP (``core.crnbp``): a torus substitute of Kumar et al. 2021's planar
Jupiter-Europa 3:4 resonant orbit, continued to the physical Io mass with
Ganymede at its physical (non-rate-idealized) synodic rate and Io exactly
Laplace-slaved, at the physical libration-center phase (see
``docs/notes/2026-07-27-724-final-confirmation-n5-torus-novelty.md``). The
user has authorized starting this project's OWN V0-V5 vetting chain on that
object, mirroring the chain already run for the OTHER genuine novel finding
in this project (`#701` Uranus Umbriel-Titania CCR4BP torus-homoclinic
connection, `#701`->`#708`). This module is `#726`, the first step: does the
idealized torus survive being propagated under REAL (SPICE) Jupiter-system
dynamics, or does it collapse -- the exact question `#704`
(``search.ccr4bp_real_ephemeris_consistency``) answered for `#701`'s
connection.

Extends `#704`'s module, does not rebuild it
---------------------------------------------
`#704`'s own module has two layers: (1) fully GENERIC math primitives with no
Uranus-specific hardcoding at all -- :func:`osculating_frame`,
:func:`~cyclerfinder.search.ccr4bp_real_ephemeris_consistency.nondim_state_to_inertial`/
``inertial_state_to_nondim`` (frame conversion, take ``mu``/``l_km``/
``v_unit_km_s`` as explicit parameters), :func:`~....real_nbody_rhs` and
:func:`~....propagate_real` (central-point-mass + N named third-body
perturbers, ``perturber_moons`` is already an arbitrary tuple -- extending
from 2 real perturbers (Umbriel, Titania) to 3 (Europa, Io, Ganymede) needs
NO code change there at all), :func:`~....tu_to_seconds`, and the
:class:`~....ConsistencyCheckResult` diagnostic dataclass; and (2) a TOP-LEVEL
orchestration function (``check_connection_survives_real_ephemeris``) that
IS Uranus/Umbriel-specific (hardcodes ``PRIMARIES["Uranus"]``, the literal
string ``"Umbriel"`` as the frame-defining base moon, and the Uranus
``URA111``-kernel SPICE-furnish/kclear dance from `#312`'s V4-strict
gauntlet). Layer (1) is imported and reused HERE VERBATIM, unmodified,
exactly as this project's own established discipline requires (see e.g.
``ccr4bp_heteroclinic_search``'s own docstring on cross-module private-helper
reuse). Layer (2) is NOT edited in place (that file is a previously-closed,
independently-verified `#704` deliverable; per this whole project's `#689`
discipline of never touching a validated module for a new extension) --
instead this module provides its own thin top-level driver,
:func:`check_torus_survives_real_ephemeris`, built from the SAME reusable
primitives, for a Jupiter/Europa/Io/Ganymede system and a TORUS (not a
manifold connection) object.

Torus vs. connection: what "departure/target" means here
-----------------------------------------------------------
`#704`'s connection had a natural departure state (the unstable-manifold
point before any flow) and target state (the stable-manifold state the
idealized connection's own residual already drove to near-coincidence). A
torus has no such pair. The natural, honest analogue used here: pick a torus
point ``u(theta1_0, theta2_0)`` as the "departure" state, and define the
"target" as the SAME idealized model's own nonlinear forward flow of that
state for a physically meaningful window ``t_window_tu`` (via
:func:`cyclerfinder.core.crnbp.propagate_crnbp` -- the FULL, coupling-term-
included EOM, the same one :func:`cyclerfinder.search.variational_crnbp_torus._independent_closure`
uses as its own independent nonlinear check). This is NOT the torus's own
algebraic self-consistency (that is what ``residual_rms``/``closure_residual``
already measure, purely within the idealized model) -- it is "does the
idealized model's own best forward-flow prediction, converted to a real
physical state, match what really happens under Jupiter + real Io + real
Europa + real Ganymede SPICE ephemeris." ``t_window_tu`` defaults to
``torus.period`` (one full forcing-clock -- Ganymede-synodic -- period), a
natural, principled window length, mirroring `#704`'s own choice of using the
connection's own natural flight time.

Central-body GM convention (inherited, not re-derived)
---------------------------------------------------------
Both ``core.crnbp``'s own idealized nondimensionalization AND `#704`'s own
real force model use the "system GM as central term" approximation: the
JPL-tabulated ``PRIMARIES["Jupiter"]``/``PRIMARIES["Uranus"]`` values are
SYSTEM GMs (already include the Galilean/Uranian moons' own mass), yet
``core.crnbp.jupiter_europa_io_ganymede_default``'s own ``denom = gm_j+gm_e``
ADDS Europa's GM again (documented there as an accepted <=2.4e-4-relative
idealization, inherited from ``core.ccr4bp``), and `#704`'s own
``real_nbody_rhs`` central term likewise uses ``PRIMARIES["Uranus"]`` (system
GM) WHILE ALSO summing Umbriel's and Titania's own GM as separate third-body
pulls -- the identical class of approximation, one level up. This module
keeps EXACTLY that established convention (``mu_uranus=PRIMARIES["Jupiter"]``
passed to the reused, unmodified :func:`~....propagate_real`) for two
reasons: (1) it is what makes the reduction-test positive control below an
apples-to-apples check (both the idealized-substitute path and the real path
use the SAME central-GM value); (2) it is the SAME approximation this whole
Jovian CCR4BP/CRNBP arc (`#689`-`#724`) already accepts throughout, not a new
one introduced here. For Jupiter the folded-in Galilean total GM is
``sum(Io,Europa,Ganymede,Callisto) / GM_Jupiter_sys = 26229.778/1.26686534e8
~= 2.07e-4`` relative -- the SAME order of magnitude as the Uranus case
`#704` already accepted (Umbriel+Titania's own mu ~1.469e-5, smaller because
Uranus's moons are collectively far less massive relative to their primary
than the Galileans are to Jupiter; empirically characterized below by this
module's own reduction-test positive control, which reports the ACTUAL
gap rather than assuming a tolerance).

Kernel loading: reuses the ALREADY-VALIDATED persistent Jupiter SPICE
backend, not `#704`'s own furnish/kclear dance
--------------------------------------------------------------------------
Unlike Uranus (for which `#704` built its OWN furnish/kclear SPICE wiring
from scratch, because none existed), a validated, already-in-production
Jupiter-system real-ephemeris backend already exists in this codebase:
:class:`cyclerfinder.core.ephemeris.Ephemeris` (``model="spice",
center="Jupiter"``), which furnishes ``jup365.bsp`` ONCE per interpreter
session (module-level guard) and is consumed read-only by, among others,
``search.ieg_seed``. This module reuses THAT (via ``Ephemeris(...).state``,
which already has exactly the ``(moon_name, et_seconds) -> (r_km, v_km_s)``
:data:`~cyclerfinder.search.ccr4bp_real_ephemeris_consistency.MoonStateFn`
signature the reused ``propagate_real`` expects) rather than re-deriving a
Jupiter-specific furnish/kclear helper -- reuse over duplication, one level
further than `#704` had the option of. Consequence: this module deliberately
does NOT call ``spice.kclear()`` after its SPICE calls (unlike `#704`'s own
Uranus driver) -- doing so would silently invalidate the ``Ephemeris``
class's own module-level "already furnished" flag for every OTHER caller in
the same process (e.g. ``ieg_seed``-based tests run in the same session).
The leapseconds kernel (needed only for the UTC-string -> ET conversion,
``spice.str2et``; ``jup365.bsp`` itself needs no LSK when ET floats are
passed directly to ``spkezr``, per
``cyclerfinder.verify.spice_kernels.ensure_jup365_kernel``'s own docstring)
is the SAME vendored ``naif0012.tls`` this project already commits
in-repo (``cyclerfinder.verify.spice_kernels.ensure_leapseconds_kernel``,
see the ``project_naif_lsk_vendored`` memory) -- furnished once here too
(``spice.furnsh`` is idempotent for the same file).

No catalogue writeback; this module produces a CONSISTENCY-CHECK result, not
a vetted discovery. Scope is exactly this ONE check (`#726`): epoch-
robustness scanning (the `#705` analogue), schema design, and writeback are
explicitly separate, not-yet-registered future tasks.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.crnbp as crnbp
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v4_uranus_strict import _ephemeris_time_seconds
from cyclerfinder.search.ccr4bp_real_ephemeris_consistency import (
    ConsistencyCheckResult,
    MoonStateFn,
    inertial_state_to_nondim,
    nondim_state_to_inertial,
    propagate_real,
    tu_to_seconds,
)
from cyclerfinder.search.variational_crnbp_torus import (
    CRNBPTorusVariationalResult,
    evaluate_torus_state,
)
from cyclerfinder.verify.spice_kernels import ensure_leapseconds_kernel

__all__ = [
    "L_KM",
    "PERTURBER_MOONS_DEFAULT",
    "check_torus_survives_real_ephemeris",
    "idealized_moon_state_fn_crnbp",
    # Re-exported so callers doing the nondim<->inertial round-trip don't need
    # a second import of the #704 module.
    "inertial_state_to_nondim",
    "jupiter_spice_moon_state_fn",
    "nondim_state_to_inertial",
    "v_unit_km_s",
]

#: Jupiter-Europa length unit (km): Europa's own SMA about Jupiter -- the
#: SAME convention ``core.ccr4bp``/``core.crnbp`` already use throughout this
#: arc (``tests/core/test_ccr4bp.py``'s own ``l_km=671100.0``).
L_KM: float = SATELLITES["Europa"].sma_km

#: The three real Galilean-moon perturbers this check propagates against --
#: Europa (the base/frame-defining body, exactly as Umbriel was for `#704`)
#: plus Io and Ganymede (the two CRNBP extra perturbers). One more real
#: perturber than `#704`'s own two (Umbriel, Titania).
PERTURBER_MOONS_DEFAULT: tuple[str, ...] = ("Europa", "Io", "Ganymede")

# A real guard (not just a redundant-looking idempotent furnsh call) is
# NECESSARY here, found during #729's own build: repeatedly calling
# ``spice.furnsh`` on the SAME already-loaded file does NOT dedupe in
# CSPICE_N0067 -- each call grows the KEEPER kernel pool by one file, and
# after ~5300 total calls (the `MAXFIL` limit) `spiceypy.utils.exceptions.
# SpiceNOMOREROOM` is raised. #726's own module docstring's claim that
# "spice.furnsh is idempotent for the same file" is FALSE at scale -- true
# only in the sense that it doesn't change the loaded CONTENT, not that
# repeat calls are free/safe. #824: the guard itself now queries SPICE's
# own pool via kinfo() at the call site below, rather than a module-level
# flag (which goes stale whenever any same-process caller elsewhere calls
# spice.kclear(), wiping the whole pool, without this module knowing).


def v_unit_km_s() -> float:
    """Physical velocity unit (km/s) for the Jupiter-Europa base pair:
    ``L_KM * n_europa``, ``n_europa = sqrt((GM_Jupiter_sys + GM_Europa) / L_KM**3)``
    -- the SAME ``denom = gm_j + gm_e`` convention
    :func:`cyclerfinder.core.crnbp.jupiter_europa_io_ganymede_default` uses
    (see the module docstring's "Central-body GM convention" section; mirrors
    :func:`cyclerfinder.core.ccr4bp_umbriel_titania.v_unit_km_s` one system
    over)."""
    gm_j = PRIMARIES["Jupiter"]
    gm_e = SATELLITES["Europa"].mu_km3_s2
    n_europa = math.sqrt((gm_j + gm_e) / L_KM**3)
    return L_KM * n_europa


def jupiter_spice_moon_state_fn(
    moon_name: str, et_seconds: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Real SPICE-sourced Jupiter-centred moon state (thin wrapper on the
    ALREADY-VALIDATED :class:`cyclerfinder.core.ephemeris.Ephemeris` SPICE
    backend -- see the module docstring's "Kernel loading" section for why
    this is reused rather than a fresh furnish/kclear helper). Furnishes
    ``jup365.bsp`` lazily on first call (module-level guard inside
    ``Ephemeris``); does NOT kclear.
    """
    ephem = Ephemeris(model="spice", center="Jupiter")
    return ephem.state(moon_name, et_seconds)


def idealized_moon_state_fn_crnbp(
    system: crnbp.CRNBPSystem, l_km: float, v_unit_km_s_: float
) -> MoonStateFn:
    """Build a :data:`MoonStateFn` from the IDEALIZED circular-coplanar CRNBP
    model -- the reduction-test substitute (see
    :func:`cyclerfinder.search.ccr4bp_real_ephemeris_consistency.idealized_moon_state_fn`,
    which this mirrors one perturber further: Europa (base) + Io + Ganymede
    (``system.perturbers``, in the established ``(io, ganymede)`` order --
    see ``core.crnbp.jupiter_europa_io_ganymede_default``) instead of Umbriel
    (base) + Titania.

    Returns Jupiter-centred positions (matching
    :func:`jupiter_spice_moon_state_fn`'s own ``observer="Jupiter"``
    convention): Europa is exactly ``l_km`` from Jupiter (the primary
    separation IS the length unit); each perturber's own barycentre-relative
    circular position (radius ``l_km * p.a``, rotating-frame angle
    ``p.theta0 + p.omega * t_nondim`` -- exactly
    :meth:`cyclerfinder.core.crnbp.CRNBPPerturber.position`) is shifted by the
    SAME ``+mu`` barycentre-to-Jupiter offset
    :func:`~cyclerfinder.search.ccr4bp_real_ephemeris_consistency.nondim_state_to_inertial`
    applies to the spacecraft state, so all bodies stay consistently
    Jupiter-centred (the exact bug class `#704`'s own module docstring
    documents catching during its own development).
    """
    n1 = v_unit_km_s_ / l_km  # rad/s, base Jupiter-Europa mean motion
    mu = system.mu
    perturbers_by_name = {"Io": system.perturbers[0], "Ganymede": system.perturbers[1]}

    def _fn(moon_name: str, et_seconds: float) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        t_s = et_seconds  # "seconds since the idealized t=0 epoch" (see docstring)
        theta_frame = n1 * t_s
        mu_shift_r = l_km * mu * np.array([math.cos(theta_frame), math.sin(theta_frame), 0.0])
        mu_shift_v = l_km * mu * n1 * np.array([-math.sin(theta_frame), math.cos(theta_frame), 0.0])
        if moon_name == "Europa":
            r_europa_km = l_km
            theta = theta_frame
            r = np.array([r_europa_km * math.cos(theta), r_europa_km * math.sin(theta), 0.0])
            v = np.array(
                [-r_europa_km * n1 * math.sin(theta), r_europa_km * n1 * math.cos(theta), 0.0]
            )
            return r, v
        if moon_name in perturbers_by_name:
            p = perturbers_by_name[moon_name]
            r_p_km = l_km * p.a
            theta_rot = p.theta0 + p.omega * (n1 * t_s)  # rotating-frame angle
            theta = theta_frame + theta_rot  # inertial angle
            omega_p_inertial = n1 * (1.0 + p.omega)
            r = np.array([r_p_km * math.cos(theta), r_p_km * math.sin(theta), 0.0]) + mu_shift_r
            v = (
                np.array(
                    [
                        -r_p_km * omega_p_inertial * math.sin(theta),
                        r_p_km * omega_p_inertial * math.cos(theta),
                        0.0,
                    ]
                )
                + mu_shift_v
            )
            return r, v
        raise ValueError(f"idealized_moon_state_fn_crnbp: unsupported moon {moon_name!r}")

    return _fn


def check_torus_survives_real_ephemeris(
    epoch0_utc: str,
    torus: CRNBPTorusVariationalResult,
    theta1_0: float,
    theta2_0: float,
    t_window_tu: float | None = None,
    *,
    perturber_moons: tuple[str, ...] = PERTURBER_MOONS_DEFAULT,
    moon_state_fn: MoonStateFn | None = None,
    rtol: float = 1e-12,
    atol: float = 1e-6,
) -> ConsistencyCheckResult:
    """Run one real-ephemeris consistency check on the N=5 CRNBP torus at one
    epoch and one torus point.

    ``theta1_0``/``theta2_0`` select the torus "departure" point (see module
    docstring). ``t_window_tu`` (nondim TU) defaults to ``torus.period`` (one
    full Ganymede-synodic forcing-clock period). ``moon_state_fn`` defaults to
    :func:`jupiter_spice_moon_state_fn` (real SPICE); tests inject
    :func:`idealized_moon_state_fn_crnbp` instead for the reduction-test
    positive control.
    """
    if t_window_tu is None:
        t_window_tu = torus.period
    system = torus.system
    l_km = L_KM
    v_unit = v_unit_km_s()
    if moon_state_fn is None:
        moon_state_fn = jupiter_spice_moon_state_fn

    u0 = evaluate_torus_state(torus, theta1_0, theta2_0)
    departure_state6 = np.array([u0[0], u0[1], 0.0, u0[2], u0[3], 0.0], dtype=np.float64)

    # The idealized model's OWN best-estimate "target": full nonlinear
    # (coupling-term-included) forward flow of the departure state for
    # t_window_tu, exactly as variational_crnbp_torus._independent_closure's
    # own independent check propagates.
    arc = crnbp.propagate_crnbp(system, departure_state6, t_window_tu, rtol=1e-13, atol=1e-13)
    target_state6 = arc.state_f

    # Leapseconds kernel for str2et only (jup365.bsp itself needs none --
    # see module docstring). A real guard, not just a redundant-looking
    # idempotent call, is required here (see the module-level comment
    # above for the MAXFIL-exhaustion reason) -- #824: that guard now
    # queries SPICE's own pool via kinfo() rather than a local "did we
    # furnish" flag, which goes stale whenever ANY same-process caller
    # elsewhere calls spice.kclear() (wipes the whole pool) without this
    # module knowing. Not kclear'd here (see docstring).
    import spiceypy as spice

    lsk_path = ensure_leapseconds_kernel()
    try:
        spice.kinfo(lsk_path)
    except spice.utils.exceptions.NotFoundError:
        spice.furnsh(lsk_path)
    et0 = _ephemeris_time_seconds(epoch0_utc)

    r_europa0, v_europa0 = moon_state_fn("Europa", et0)
    r0_km, v0_km_s = nondim_state_to_inertial(
        departure_state6, r_europa0, v_europa0, l_km, v_unit, mu=system.mu
    )

    t_window_s = t_window_tu * tu_to_seconds(l_km, v_unit)
    prop = propagate_real(
        r0_km,
        v0_km_s,
        t_window_s,
        et0=et0,
        mu_uranus=PRIMARIES["Jupiter"],  # central-term GM; see module docstring
        perturber_moons=perturber_moons,
        moon_state_fn=moon_state_fn,
        rtol=rtol,
        atol=atol,
    )

    et_target = et0 + t_window_s
    r_europa_target, v_europa_target = moon_state_fn("Europa", et_target)
    r_target_km, v_target_km_s = nondim_state_to_inertial(
        target_state6, r_europa_target, v_europa_target, l_km, v_unit, mu=system.mu
    )

    if prop.success:
        pos_gap_km = float(np.linalg.norm(prop.r_f_km - r_target_km))
        vel_gap_km_s = float(np.linalg.norm(prop.v_f_km_s - v_target_km_s))
        notes = "converged"
    else:
        pos_gap_km = float("inf")
        vel_gap_km_s = float("inf")
        notes = "real-ephemeris propagation failed"

    return ConsistencyCheckResult(
        epoch0_utc=epoch0_utc,
        et0=et0,
        l_km=l_km,
        v_unit_km_s=v_unit,
        t_u_tu=t_window_tu,
        t_u_seconds=t_window_s,
        r0_km=r0_km,
        v0_km_s=v0_km_s,
        r_f_km=prop.r_f_km,
        v_f_km_s=prop.v_f_km_s,
        r_target_km=r_target_km,
        v_target_km_s=v_target_km_s,
        pos_gap_km=pos_gap_km,
        vel_gap_km_s=vel_gap_km_s,
        propagation_success=prop.success,
        closest_approach_km=prop.closest_approach_km,
        idealized_off_torus_km_for_scale=float("nan"),
        notes=notes + " (idealized_off_torus_km_for_scale: N/A -- a torus has no manifold "
        "off-torus gate; use torus.residual_rms/closure_residual for the "
        "idealized model's own self-consistency scale instead)",
    )
