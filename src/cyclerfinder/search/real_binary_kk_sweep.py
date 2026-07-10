"""Task #549: real-binary (k1,k2)-cycler sweep.

Generalizes the #504 Pluto-Charon (k1,k2)-family sweep
(:mod:`cyclerfinder.search.pluto_charon_kk_sweep`) to arbitrary real binary
systems, reusing that module's solver machinery verbatim (the fixed-Jacobi
symmetric corrector, Barden stability, winding-number topology classifier,
independent-Radau crosscheck, and the C-sweep/brentq nu=0 finder). The ONLY
new logic here is the driver: seeding from the same Ross-RT 2026 Table-I
anchors #504 used, but finishing the mu-continuation in the CALLER-SUPPLIED
target system instead of a hardcoded Pluto-Charon system.

Why this module exists rather than calling #504's ``mu_step_to_orbit``
directly: that function's "final correction" step is hardcoded to
``make_pluto_charon_system()`` regardless of the ``target_mu`` argument. For
``target_mu == PC_MU`` this is a harmless no-op (the loop's last iteration
already landed on the right mu). For any OTHER target mu it would silently
re-converge the orbit in the WRONG system (Pluto-Charon's mu, not the
caller's), corrupting the result. :func:`mu_step_to_system` below is the
minimal fix: identical mu-stepping loop, but the final correction uses the
caller's ``target_system`` (real l_km/t_s included) as it should.

Public API
----------
mu_step_to_system(...)  -> SymmetricOrbit | None  (mu-continuation to ANY target system)
sweep_family(...)       -> SweepResult  (anchor-seeded C-sweep for one (k1,k2))
sweep_family_grid(...)  -> SweepResult  (grid-seeded C-sweep, topologies with no anchor)
REAL_BINARY_SYSTEMS     -> dict[str, RealBinarySystem]  (sourced mu/l_km/t_s + citations)
ANCHORS                 -> dict[str, dict]  (Ross-RT 2026 Table-I anchors, from #504)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.pluto_charon_kk_sweep import (
    SweepResult,
    _build_result,
    _c_l1,
    _grid_seed_search,
    _nd_system,
    c_sweep_find_nu_zero,
)

__all__ = [
    "ANCHORS",
    "REAL_BINARY_SYSTEMS",
    "RealBinarySystem",
    "mu_step_to_system",
    "sweep_family",
    "sweep_family_grid",
]


# ---------------------------------------------------------------------------
# Sourced real-binary systems (task #549)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RealBinarySystem:
    """A real binary system's CR3BP mass ratio + physical scales, with citations."""

    name: str
    primary: str
    secondary: str
    mu: float
    mu_source: str
    l_km: float  # mutual-orbit semimajor axis (CR3BP length unit)
    l_source: str
    t_s: float  # sqrt(l_km^3 / G(m1+m2)) == P_orbital_seconds / (2*pi)
    t_source: str
    caveat: str = ""

    def to_cr3bp_system(self) -> cr3bp.CR3BPSystem:
        return cr3bp.CR3BPSystem(
            mu=self.mu,
            primary=self.primary,
            secondary=self.secondary,
            l_km=self.l_km,
            t_s=self.t_s,
        )


def _t_s_from_period_days(period_days: float) -> float:
    """CR3BP time unit (seconds) = orbital period / (2*pi) -- see cr3bp_system()."""
    return float(period_days * 86400.0 / (2.0 * np.pi))


REAL_BINARY_SYSTEMS: dict[str, RealBinarySystem] = {
    "patroclus-menoetius": RealBinarySystem(
        name="Patroclus-Menoetius",
        primary="Patroclus",
        secondary="Menoetius",
        # mu = D_Menoetius^3 / (D_Patroclus^3 + D_Menoetius^3), assuming EQUAL
        # DENSITY (standard for a co-formed contact/near-equal Trojan binary;
        # Buie et al. 2015 itself note the two diameters are consistent with
        # equal density given the combined system mass/volume).
        # D_Patroclus = 113+/-3 km, D_Menoetius = 104+/-3 km (occultation).
        mu=0.4380719233604685,
        mu_source=(
            "Buie, M. W. et al. (2015), AJ 149, 113, 'Size and Shape from "
            "Stellar Occultation Observations of the Double Jupiter Trojan "
            "Patroclus and Menoetius' -- D_Patroclus=113+/-3 km, "
            "D_Menoetius=104+/-3 km (2013-10-21 occultation). mu = D_M^3/"
            "(D_P^3+D_M^3) = 0.4381, EQUAL-DENSITY ASSUMPTION (standard for "
            "a co-formed binary; consistent with the paper's own note that "
            "the two diameters reproduce the combined system density of "
            "0.88 g/cm^3 without invoking a density contrast)."
        ),
        l_km=692.5,
        l_source=(
            "Buie, M. W. et al. (2024), AJ 167, 104, Table (mutual-orbit fit): a=692.5+/-4.0 km."
        ),
        t_s=_t_s_from_period_days(4.282754),
        t_source=(
            "Buie, M. W. et al. (2024), AJ 167, 104, Table (mutual-orbit fit): "
            "P=4.282754+/-0.000023 d."
        ),
        caveat=(
            "mu is a DERIVED equal-density estimate, not a directly measured "
            "mass ratio -- no Lucy-flyby (2033) mass measurement exists yet."
        ),
    ),
    "didymos-dimorphos": RealBinarySystem(
        name="Didymos-Dimorphos",
        primary="Didymos",
        secondary="Dimorphos",
        mu=0.007900055116663605,
        mu_source=(
            "Naidu, S. P. et al. (2020), Icarus 348, 113777 -- Didymos system "
            "mass (5.4+/-0.4)e11 kg (mutual-orbit Kepler fit); Dimorphos mass "
            "~4.3e9 kg assuming EQUAL DENSITY with Didymos (same convention "
            "used pre-DART mission design and cited post-impact, e.g. Daly, "
            "R. T. et al. (2023), Nature 616, 443). mu = 4.3e9/(5.4e11+4.3e9) "
            "= 0.0079."
        ),
        l_km=1.206,
        l_source=(
            "Naidu, S. P. et al. (2020), Icarus 348, 113777 -- pre-impact "
            "mutual-orbit a=1.206+/-0.035 km."
        ),
        t_s=_t_s_from_period_days(11.921473 / 24.0),
        t_source=(
            "Naidu, S. P. et al. (2020), Icarus 348, 113777 -- pre-impact "
            "mutual-orbit P=11.921473+/-0.000044 hr."
        ),
        caveat=(
            "Dimorphos's mass is NOT directly measured (DART/LICIACube could "
            "not resolve it; Daly et al. 2023 and the post-DART Nature "
            "Astronomy 2024 physical-properties paper (arXiv:2403.00667) both "
            "state the mass is inferred only via an ASSUMED bulk density, "
            "quoting a plausible density range 1500-3300 kg/m^3 vs the "
            "assumed-equal 2170 kg/m^3 -- i.e. mu could plausibly range "
            "~0.0055-0.012, a factor >2. ESA Hera (arriving Didymos "
            "2026/2027) will measure it directly. Pre-impact orbital "
            "elements (a, P) are used here (the DART impact itself perturbed "
            "the orbit ~-33 min in period, still relaxing per post-impact "
            "papers; pre-impact is the natural, near-circular reference "
            "state a CR3BP model assumes)."
        ),
    ),
    "orcus-vanth": RealBinarySystem(
        name="Orcus-Vanth",
        primary="Orcus",
        secondary="Vanth",
        mu=0.137,
        mu_source=(
            "Brown, M. E. & Butler, B. J. (2023), PSJ 4, 178 (arXiv:2307.04848), "
            "'Masses and densities of dwarf planet satellites measured with "
            "ALMA' -- direct ALMA astrometric mass measurement: 'Vanth "
            "contains 13.7+/-1.3% of the mass of the system,' i.e. mu = "
            "m_Vanth/(m_Orcus+m_Vanth) = 0.137+/-0.013 (equivalently a "
            "satellite/primary mass ratio of 0.16+/-0.02, 'the highest of "
            "any known planet or dwarf planet')."
        ),
        l_km=8980.0,
        l_source=(
            "Brown, M. E. et al. (2010), AJ 139, 2700 (HST astrometry, "
            "refined) -- mutual-orbit a=8980+/-20 km."
        ),
        t_s=_t_s_from_period_days(9.5393),
        t_source=(
            "Brown, M. E. et al. (2010), AJ 139, 2700 (HST astrometry, "
            "refined) -- mutual-orbit P=9.5393+/-0.0001 d."
        ),
    ),
    "eris-dysnomia": RealBinarySystem(
        name="Eris-Dysnomia",
        primary="Eris",
        secondary="Dysnomia",
        # Brown & Butler 2023 report a satellite/primary mass ratio of
        # 0.0050+/-0.0035 (only a 1.5-sigma DETECTION, 1-sigma upper limit
        # 0.0085) -- flagged explicitly below; this is the weakest-sourced mu
        # of the four systems.
        mu=0.0049751243781094535,
        mu_source=(
            "Brown, M. E. & Butler, B. J. (2023), PSJ 4, 178 (arXiv:2307.04848) "
            "-- ALMA astrometric mass ratio Dysnomia/Eris = 0.0050+/-0.0035 "
            "(a 1.5-sigma detection only; 1-sigma upper limit 0.0085, "
            "3-sigma upper limit 0.015). mu = m_Dysnomia/(m_Eris+m_Dysnomia) "
            "= 0.0050/1.0050 = 0.004975 (central value; upper 1-sigma bound "
            "mu=0.00833)."
        ),
        l_km=37273.0,
        l_source=(
            "'The Eris/Dysnomia system I: the orbit of Dysnomia' (2020), "
            "Icarus, arXiv:2009.13733 -- a=37273+/-64 km."
        ),
        t_s=_t_s_from_period_days(15.785899),
        t_source=(
            "'The Eris/Dysnomia system I: the orbit of Dysnomia' (2020), "
            "Icarus, arXiv:2009.13733 -- P=15.785899+/-0.000050 d."
        ),
        caveat=(
            "mu is only a 1.5-sigma ALMA mass DETECTION (Brown & Butler 2023); "
            "treat as weakly constrained -- the 1-sigma upper bound mu=0.00833 "
            "is nearly 70% larger than the central value used for the sweep."
        ),
    ),
}


# ---------------------------------------------------------------------------
# Ross-RT 2026 Table-I anchors (identical numbers to #504's
# pluto_charon_kk_sweep.py module docstring / sweep_* functions -- reused
# verbatim, not re-sourced).
# ---------------------------------------------------------------------------

ANCHORS: dict[str, dict[str, Any]] = {
    "mu001_11": dict(
        anchor_mu=0.001,
        x0=-0.647047499999966,
        jacobi=3.031605708907296,
        period=14.774502790974823,
        k1=1,
        k2=1,
        hc=None,
    ),
    "mu01215_11": dict(
        anchor_mu=0.012150584270572,
        x0=-0.768217354461248,
        jacobi=3.151175879917331,
        period=10.291893641936499,
        k1=1,
        k2=1,
        hc=None,
    ),
    "mu05_11": dict(
        # Fundamental (1,1) period T1 = T/3 (published T is the 3rd iterate;
        # see catalogue row ross-rt-mu05-cycler-11-2026).
        anchor_mu=0.5,
        x0=-0.519689929077496,
        jacobi=3.628400000000000,
        period=2.930671187154082,
        k1=1,
        k2=1,
        hc=None,
    ),
    "mu01_32": dict(
        anchor_mu=0.1,
        x0=-0.694376003123377,
        jacobi=3.573367616904619,
        period=12.295263874014290,
        k1=3,
        k2=2,
        hc=6,
    ),
    "mu03_31": dict(
        anchor_mu=0.3,
        x0=-0.804725783387797,
        jacobi=3.701958166478617,
        period=9.094576400494693,
        k1=3,
        k2=1,
        hc=None,
    ),
    "mu01215_33": dict(
        anchor_mu=0.012150584270572,
        x0=-0.322477620583087,
        jacobi=3.183379082910527,
        period=19.503763587070285,
        k1=3,
        k2=3,
        hc=None,
    ),
}


# ---------------------------------------------------------------------------
# Generalized mu-continuation (the one genuinely new function: #504's
# mu_step_to_orbit hardcodes the final system to Pluto-Charon; this fixes
# that so the SAME stepping logic lands in the CALLER's target system).
# ---------------------------------------------------------------------------


def mu_step_to_system(
    anchor_mu: float,
    target_system: cr3bp.CR3BPSystem,
    anchor_x0: float,
    anchor_jacobi: float,
    anchor_period: float,
    *,
    hc: int | None = None,
    sign: float = -1.0,
    n_steps: int = 40,
    tol: float = 1e-10,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> cp.SymmetricOrbit | None:
    """Step mu from anchor_mu to target_system.mu, finishing IN target_system.

    Identical algorithm to :func:`pluto_charon_kk_sweep.mu_step_to_orbit`
    (holds Jacobi fixed, re-corrects x0/T at each mu step) except the final
    correction runs in ``target_system`` (with its real l_km/t_s) instead of
    a hardcoded Pluto-Charon system -- see module docstring.
    """
    mus = np.linspace(anchor_mu, target_system.mu, n_steps + 1)[1:]
    x0_cur = anchor_x0
    t_cur = anchor_period
    jacobi = anchor_jacobi

    sys_anc = _nd_system(anchor_mu)
    try:
        o = cp.correct_symmetric_fixed_jacobi(
            sys_anc,
            x0_cur,
            jacobi,
            t_cur,
            ydot0_sign=sign,
            half_crossings=hc,
            tol=tol,
            rtol=rtol,
            atol=atol,
        )
    except ValueError:
        return None
    if not o.converged:
        return None
    x0_cur, t_cur = o.x0, o.period

    for mu_next in mus:
        sys_next = _nd_system(mu_next)
        try:
            o = cp.correct_symmetric_fixed_jacobi(
                sys_next,
                x0_cur,
                jacobi,
                t_cur,
                ydot0_sign=sign,
                half_crossings=hc,
                tol=tol,
                rtol=rtol,
                atol=atol,
            )
        except ValueError:
            return None
        if not o.converged:
            return None
        x0_cur, t_cur = o.x0, o.period

    try:
        o_final = cp.correct_symmetric_fixed_jacobi(
            target_system,
            x0_cur,
            jacobi,
            t_cur,
            ydot0_sign=sign,
            half_crossings=hc,
            tol=tol,
            rtol=rtol,
            atol=atol,
        )
    except ValueError:
        return None
    if not o_final.converged:
        return None
    return o_final


# ---------------------------------------------------------------------------
# Per-(k1,k2) sweep drivers
# ---------------------------------------------------------------------------


def sweep_family(
    target_system: cr3bp.CR3BPSystem,
    anchor_key: str,
    *,
    c_band: float = 0.3,
    n_steps: int = 40,
) -> SweepResult:
    """Anchor-seeded mu-continuation + C-sweep for one (k1,k2) family."""
    a = ANCHORS[anchor_key]
    k1, k2 = a["k1"], a["k2"]
    c_l1 = _c_l1(target_system.mu)

    orbit_seed = mu_step_to_system(
        a["anchor_mu"],
        target_system,
        a["x0"],
        a["jacobi"],
        a["period"],
        hc=a["hc"],
        sign=-1.0,
        n_steps=n_steps,
    )
    if orbit_seed is None:
        return SweepResult(
            k1=k1,
            k2=k2,
            stable_found=False,
            method=f"mu_step_from_mu{a['anchor_mu']:.6g}_[{anchor_key}]",
            note="mu-continuation failed to converge to target mu",
        )

    if a["hc"] is not None:
        # Forced-hc families (e.g. (3,2) with hc=6, per #504) are fragile
        # below the anchor's own C: walking c_lo down from the seed forces
        # the corrector to hold a high-multiplicity crossing count far outside
        # its natural C range, which is pathologically slow/non-convergent.
        # #504's own sweep_32_positive_control avoids this by sweeping UPWARD
        # ONLY, starting exactly at the anchor's C -- mirrored here.
        c_lo = orbit_seed.jacobi
        c_hi = min(c_l1 - 0.002, orbit_seed.jacobi + 2.0 * c_band)
    else:
        c_lo = max(2.5, orbit_seed.jacobi - c_band)
        c_hi = min(c_l1 - 0.002, orbit_seed.jacobi + c_band)
    if c_lo >= c_hi:
        return SweepResult(
            k1=k1,
            k2=k2,
            stable_found=False,
            method=f"mu_step_from_mu{a['anchor_mu']:.6g}_[{anchor_key}]",
            note=f"seed jacobi={orbit_seed.jacobi:.4f} incompatible with C_L1={c_l1:.4f}",
        )
    orbit_stable = c_sweep_find_nu_zero(
        target_system,
        orbit_seed.x0,
        orbit_seed.jacobi,
        orbit_seed.period,
        hc=a["hc"],
        sign=-1.0,
        c_lo=c_lo,
        c_hi=c_hi,
        n_coarse=60,
    )
    method = f"mu_step_from_mu{a['anchor_mu']:.6g}_[{anchor_key}]_then_c_sweep"
    if orbit_stable is None:
        return SweepResult(
            k1=k1,
            k2=k2,
            stable_found=False,
            method=method,
            note="no stable (|nu|<1) window found in C-sweep range",
        )
    res = _build_result(k1, k2, target_system, orbit_stable, method)
    if not res.topology_ok:
        return SweepResult(
            k1=k1,
            k2=k2,
            stable_found=False,
            method=method,
            note=(
                f"stable orbit found but wrong topology (got prograde="
                f"{res.prograde}, reaches_secondary={res.reaches_secondary}); "
                "clean negative"
            ),
        )
    return res


def sweep_family_grid(
    target_system: cr3bp.CR3BPSystem,
    k1: int,
    k2: int,
    *,
    x0_grid: np.ndarray,
    c_grid: np.ndarray,
    hc_list: tuple[int, ...],
    period_guess: float,
    per_call_timeout: int = 4,
) -> SweepResult:
    """Grid-seeded (no Table-I anchor) C-sweep for one (k1,k2) family."""
    c_l1 = _c_l1(target_system.mu)
    seed_orbit = _grid_seed_search(
        target_system,
        k1,
        k2,
        x0_grid,
        c_grid,
        hc_list,
        period_guess,
        per_call_timeout=per_call_timeout,
    )
    if seed_orbit is None:
        return SweepResult(
            k1=k1,
            k2=k2,
            stable_found=False,
            method=f"grid_search_{len(x0_grid)}x{len(c_grid)}x{len(hc_list)}",
            note=(
                f"no ({k1},{k2}) orbit found in grid "
                f"(x0 in [{x0_grid.min():.2f},{x0_grid.max():.2f}], "
                f"C in [{c_grid.min():.4f},{c_grid.max():.4f}], hc in {hc_list})"
            ),
        )
    c_lo = max(2.5, seed_orbit.jacobi - 0.3)
    c_hi = min(c_l1 - 0.002, seed_orbit.jacobi + 0.3)
    orbit_stable = c_sweep_find_nu_zero(
        target_system,
        seed_orbit.x0,
        seed_orbit.jacobi,
        seed_orbit.period,
        hc=None,
        sign=-1.0,
        c_lo=c_lo,
        c_hi=c_hi,
        n_coarse=60,
    )
    if orbit_stable is None:
        return SweepResult(
            k1=k1,
            k2=k2,
            stable_found=False,
            method="grid_then_c_sweep",
            note="no stable window",
        )
    return _build_result(k1, k2, target_system, orbit_stable, "grid_seed_then_c_sweep")
