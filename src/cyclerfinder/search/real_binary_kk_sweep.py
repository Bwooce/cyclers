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
mu_step_to_system_tracking_c_l1(...) -> SymmetricOrbit | None  (mu-continuation
                                        DOWN, tracking C below C_L1(mu); #627)
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
    "SweepResult",
    "mu_step_to_system",
    "mu_step_to_system_tracking_c_l1",
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
    # ------------------------------------------------------------------
    # Task #657 round-2 systems (2026-07-19): three near-equal-mu real
    # binaries #549 never tried. Independently re-sourced per
    # [[feedback_digest_not_adoption]] -- do NOT trust #654's own prose
    # mu figures without re-deriving from the primary literature (done
    # below for all three).
    # ------------------------------------------------------------------
    "sila-nunam": RealBinarySystem(
        name="Sila-Nunam",
        primary="Sila",
        secondary="Nunam",
        # No individually-measured masses exist (the two bodies are too
        # near-equal in brightness for any barycentric wobble to be
        # detected -- Grundy et al. 2012 fit a single Keplerian orbit for
        # Nunam relative to Sila, giving only the SYSTEM mass). mu is
        # instead derived from the reported mean V-band magnitude
        # difference (Sila brighter by 0.12 mag over 8 visits) assuming
        # EQUAL ALBEDO (-> radius ratio from flux ratio) and EQUAL DENSITY
        # (-> mass ratio from radius ratio) -- the identical convention
        # #549 already used for Patroclus-Menoetius, not a new assumption:
        # flux_Sila/flux_Nunam = 10**(0.12/2.5) = 1.116863;
        # r_Nunam/r_Sila = (flux_Nunam/flux_Sila)**0.5 = 0.946237;
        # mu = r_ratio**3 / (1 + r_ratio**3) = 0.8472/(1.8472) = 0.458648.
        mu=0.45864813809877336,
        mu_source=(
            "Grundy, W. M. et al. (2012), Icarus 220, 74 (arXiv:1204.3923), "
            "'Mutual Events in the Cold Classical Transneptunian Binary "
            "System Sila and Nunam' -- Table 1/Sec.3: 'Sila, the primary, "
            "was not always brighter, but averaged over the eight visits "
            "where separate photometry was obtained ... it was brighter by "
            "a mean of 0.12 mags.' No individual masses are resolved (near-"
            "equal brightnesses prevent detecting the photocenter-barycenter "
            "offset); mu DERIVED here from that magnitude difference under "
            "equal-albedo + equal-density assumptions (see module comment "
            "above for the arithmetic), NOT a directly measured mass ratio."
        ),
        l_km=2777.0,
        l_source=(
            "Grundy, W. M. et al. (2012), Icarus 220, 74 (arXiv:1204.3923), "
            "Table 2: mutual-orbit fit a=2777+/-19 km."
        ),
        t_s=_t_s_from_period_days(12.50995),
        t_source=(
            "Grundy, W. M. et al. (2012), Icarus 220, 74 (arXiv:1204.3923), "
            "Table 2: mutual-orbit fit P=12.50995+/-0.00036 d, e=0.020+/-0.015 "
            "(nearly circular; the CR3BP model here assumes exactly circular, "
            "as for all other REAL_BINARY_SYSTEMS entries)."
        ),
        caveat=(
            "mu is NOT a directly measured mass ratio (see mu_source) -- it "
            "is derived from a photometric magnitude difference under "
            "equal-albedo/equal-density assumptions, the same indirect-"
            "inference category as Patroclus-Menoetius's own #549 sourcing. "
            "System is doubly-synchronous (tidally locked and synchronized "
            "like Pluto-Charon, per the source paper's own Sec.4)."
        ),
    ),
    "antiope": RealBinarySystem(
        name="Antiope",
        primary="Antiope-A",
        secondary="Antiope-B",
        # Aljbaae et al. (2020) Table 1 gives INDIVIDUAL masses (from the
        # Bartczak et al. 2014 SAGE non-convex shape+dynamics model), not
        # an equal-density-assumed size ratio -- the most directly-measured
        # mu of the three new #657 systems.
        # mu = M_beta/(M_alpha+M_beta) = 4.525132e17/9.120774e17 = 0.496135.
        mu=0.4961346482217408,
        mu_source=(
            "Aljbaae, S., Prado, A. F. B. A., Sanchez, D. M. & Hussmann, H. "
            "(2020), MNRAS 496, 1645, 'Analysis of the orbital stability "
            "close to the binary asteroid (90) Antiope', Table 1 (masses "
            "from the Bartczak et al. 2014 SAGE non-convex shape model): "
            "M_alpha=4.595642e17 kg, M_beta=4.525132e17 kg -- INDIVIDUAL "
            "masses, not an assumed size/density ratio. Cross-check: "
            "Descamps et al. (2007), Icarus 187, 482, independently report "
            "component size ratio B/A=0.95+/-0.01 (near-equal, consistent)."
        ),
        l_km=176.0,
        l_source=(
            "Aljbaae et al. (2020), MNRAS 496, 1645 -- mutual-orbit "
            "separation 176+/-4 km (self-consistent with their own masses "
            "and period via Kepler's third law: GM_needed=9.134e17 kg vs "
            "their stated/tabulated total 9.14e17-9.121e17 kg, <0.2% "
            "agreement). Descamps et al. (2007) independently report "
            "a=171+/-1 km for the same system -- both cited, this system "
            "uses Aljbaae et al.'s own internally-consistent (a, P, masses) "
            "triple rather than mixing sources."
        ),
        t_s=_t_s_from_period_days(16.505046 / 24.0),
        t_source=(
            "Aljbaae et al. (2020), MNRAS 496, 1645 -- mutual-orbit period "
            "16.505046+/-0.000005 h (= 0.68771025 d); matches Descamps et "
            "al. (2007)'s independently-measured 16.5051+/-0.0001 h to "
            "5 significant figures. Orbit is circular (90 Antiope is the "
            "first doubly-synchronous main-belt asteroid discovered by "
            "direct imaging, Merline et al. 2000)."
        ),
    ),
    "lempo-hiisi": RealBinarySystem(
        name="Lempo-Hiisi",
        # CR3BP convention requires mu=m2/(m1+m2)<=0.5, i.e. the PRIMARY
        # must be the heavier body. Ragozzine et al. (2024)'s own fit finds
        # Hiisi (despite being the FAINTER body, hence its name in the MPC
        # circular convention) is actually the MORE MASSIVE of the pair --
        # so Hiisi is primary here and Lempo is secondary, the reverse of
        # the visual-brightness-based naming.
        primary="Hiisi",
        secondary="Lempo",
        # Ragozzine et al. (2024) Table 2, "Best fit" column (the single
        # self-consistent point estimate, not the marginal-posterior
        # column, which mixes non-simultaneous draws):
        # M_Lempo=5.960e18 kg, M_Hiisi=7.610e18 kg (Hiisi is ~28% MORE
        # massive despite being fainter -- their own Sec.8 conclusion:
        # "the fainter Hiisi is about 33+/-5% more massive than Lempo",
        # matching the POSTERIOR-median masses 7.657/5.725=1.3375; the
        # best-fit column used here gives 7.610/5.960=1.2769, same sign
        # and same order of magnitude, within the paper's own quoted
        # uncertainty).
        # mu = M_Lempo/(M_Hiisi+M_Lempo) = 5.960/13.570 = 0.439204.
        mu=0.4392041267501842,
        mu_source=(
            "Ragozzine, D., Pincock, S., Proudfoot, B. C. N., Spencer, D., "
            "Porter, S. & Grundy, W. (2024), 'Beyond Point Masses I: New "
            "Non-Keplerian Modeling Tools Applied to Trans-Neptunian Triple "
            "(47171) Lempo', arXiv:2403.12785, Table 2 (MultiMoon three-"
            "point-mass Bayesian fit, 'Best fit' column): Mass Lempo="
            "5.960e18 kg, Mass Hiisi=7.610e18 kg -- a genuine mass "
            "MEASUREMENT (the first for this system, superseding Benecchi "
            "et al. 2010's system-mass-only double-Keplerian fit). Hiisi is "
            "the MORE massive body, so it is CR3BP-primary here (mu<=0.5 "
            "convention) even though 'Hiisi' is the fainter-and-thus-"
            "secondary name under the MPC circular's brightness convention."
        ),
        l_km=850.0,
        l_source=(
            "Ragozzine et al. (2024), arXiv:2403.12785, Table 2, 'Best fit' "
            "column: semi-major axis of the Lempo-Hiisi mutual orbit "
            "a_2=850 km (posterior median 838(+13/-21) km; Benecchi et al. "
            "2010's earlier double-Keplerian fit gave 867+/-11 km -- the "
            "paper's own three-point-mass Bayesian fit supersedes B10)."
        ),
        # Period is NOT tabulated directly in the source (only masses and
        # semimajor axis are) -- derived here via Kepler's third law from
        # the paper's own best-fit M_Lempo+M_Hiisi and a_2, using this
        # module's G=6.67430e-20 km^3 kg^-1 s^-2 (CODATA 2018, matching
        # core/constants.py's own _G_KM3_KG_S2). Cross-check: the derived
        # 1.894 d agrees with the paper's own prose statement (Sec.5.2)
        # that Lempo's rotation gives "roughly 7.4 rotations within "
        # "Lempo-Hiisi's 1.9-day binary orbital period."
        t_s=26039.677601746138,
        t_source=(
            "DERIVED (not tabulated) via Kepler's third law from Ragozzine "
            "et al. (2024) Table 2 best-fit masses (13.570e18 kg total) and "
            "a_2=850 km: t_s=sqrt(a_2^3/(G*(M1+M2)))=26039.68 s, i.e. "
            "period=1.8937 d -- matches the paper's own prose '1.9-day "
            "binary orbital period' (Sec.5.2) to the figure they quote."
        ),
        caveat=(
            "UNMODELED PERTURBATION: Lempo is a hierarchical TRIPLE "
            "(Lempo-Hiisi inner pair + Paha outer satellite, a_3~7600 km, "
            "~9x the inner separation). This entry treats Lempo-Hiisi as "
            "an ISOLATED CR3BP binary core -- Paha's gravitational "
            "perturbation on the inner pair is NOT modeled. Ragozzine et "
            "al. (2024) themselves find M_Paha consistent with zero (upper "
            "limit ~0.5e18 kg, ~3.5% of the inner pair's mass) but flag the "
            "system as dynamically puzzling (Correia 2018's point-mass "
            "analysis found the inner binary's own fitted eccentricity "
            "e2~0.12 chaotic and disruptive on <1 Myr, orders of magnitude "
            "shorter than the age of the solar system -- unresolved in "
            "Ragozzine et al. 2024 too, Sec.7.1). ANY stable member found "
            "here for Lempo-Hiisi must be flagged with this caveat "
            "prominently: a 'stable' result in the idealized 2-body-core "
            "model could be an artifact of ignoring Paha and/or the "
            "system's own known non-Keplerian complexity, not evidence the "
            "REAL 3-body system hosts a cycler."
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


def mu_step_to_system_tracking_c_l1(
    anchor_mu: float,
    target_system: cr3bp.CR3BPSystem,
    anchor_x0: float,
    anchor_jacobi: float,
    anchor_period: float,
    *,
    hc: int | None,
    sign: float = -1.0,
    n_steps: int = 200,
    c_margin: float = 0.005,
    c_margin_alpha: float | None = None,
    n_walk: int = 20,
    tol: float = 1e-10,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> cp.SymmetricOrbit | None:
    """Step mu DOWNWARD from ``anchor_mu`` to ``target_system.mu``, tracking C
    below ``C_L1(mu)`` throughout (task #627).

    ``c_margin_alpha`` (task #629 Phase A): if given (not ``None``), overrides
    ``c_margin`` at EVERY step with a margin SCALED to that step's corridor
    width, ``c_margin_step = c_margin_alpha * (C_L1(mu_next) - 3)``, instead of
    the absolute constant ``c_margin``. #629's design read
    (docs/notes/2026-07-18-629-design-read-titan-kk-grid.md) found #627's
    fixed absolute ``c_margin`` (0.02 / 0.005) does not scale with the
    Hill-shrinking corridor (``C_L1(mu) - 3 ~ mu^(2/3)``), so a walk that
    starts inside the corridor can be forced OUT of it (even below C=3
    entirely) purely by the margin's own fixed size as mu shrinks -- an
    artifact of the parameterization, not a physical result. Tracking
    ``c_margin_alpha`` instead holds the *scaled* energy
    ``rho = (C - 3) / (C_L1(mu) - 3)`` approximately constant at
    ``1 - c_margin_alpha`` once the walk starts clamping, matching how the
    two sourced Ross-RT 2026 Table-I (1,1) anchors sit at a near-invariant
    rho (~0.79-0.80) across a 12x mu range. ``c_margin`` is ignored when
    ``c_margin_alpha`` is not ``None``.

    :func:`mu_step_to_system` holds the anchor's Jacobi constant FIXED across
    the whole mu range -- correct for #494/#549's upward extension (Earth-Moon
    to binary-star mu, where C_L1 does not fall below the anchor C along the
    walk), but WRONG going down toward a real moon's much smaller mu: C_L1(mu)
    is not monotonic-safe in this direction either, and #627 found it can drop
    BELOW the anchor's own C well before the target mu is reached (e.g.
    C_L1(Saturn-Titan mu=2.37e-4)=3.016 versus the Ross-RT 2026 Table-I (1,1)/
    (3,3) mu=0.01215 anchors at C=3.151/3.183). Once that happens, a
    fixed-Jacobi walk still numerically "converges" at every step (the
    corrector always finds SOME nearby periodic orbit) but silently drifts
    onto a topologically different, non-secondary-reaching branch -- a
    false-positive "success" that only shows up as ``reaches_secondary=False``
    downstream. This generalizes #494's own fix for the analogous problem
    going UP (``pluto_charon_kk_sweep.sweep_31``'s "Strategy B": C-walk to
    below C_L1 BEFORE mu-stepping) into a per-step interleaving that applies
    at every mu decrement, not just once: before each mu step, if the current
    C is not at least ``c_margin`` below ``C_L1(mu_next)``, walk C DOWN (at
    the current, not-yet-stepped mu) in ``n_walk`` small increments first,
    then take the mu step at the (possibly lowered) C. ``hc`` should be an
    explicit crossing index (not ``None``) so the branch identity is held
    fixed through every sub-step -- an auto-redetected crossing index can snap
    onto a different, unrelated branch the moment the orbit's shape changes
    (also found empirically in #627).

    Returns the final orbit in ``target_system`` (its real l_km/t_s), or
    ``None`` if any C-walk or mu-step sub-correction fails to converge.
    """

    def _try(sys_: cr3bp.CR3BPSystem, x0: float, c: float, t: float) -> cp.SymmetricOrbit | None:
        try:
            o = cp.correct_symmetric_fixed_jacobi(
                sys_, x0, c, t, ydot0_sign=sign, half_crossings=hc, tol=tol, rtol=rtol, atol=atol
            )
        except ValueError:
            return None
        return o if o.converged else None

    mus = np.linspace(anchor_mu, target_system.mu, n_steps + 1)
    x0_cur, t_cur, c_cur = anchor_x0, anchor_period, anchor_jacobi
    mu_cur = float(mus[0])
    o = _try(_nd_system(mu_cur), x0_cur, c_cur, t_cur)
    if o is None:
        return None
    x0_cur, t_cur = o.x0, o.period

    for mu_next in mus[1:]:
        mu_next = float(mu_next)
        c_l1_next = _c_l1(mu_next)
        margin_here = c_margin_alpha * (c_l1_next - 3.0) if c_margin_alpha is not None else c_margin
        c_target = min(c_cur, c_l1_next - margin_here)
        if c_target < c_cur - 1e-15:
            sys_here = _nd_system(mu_cur)
            for c_walk in np.linspace(c_cur, c_target, n_walk)[1:]:
                o = _try(sys_here, x0_cur, float(c_walk), t_cur)
                if o is None:
                    return None
                x0_cur, t_cur, c_cur = o.x0, o.period, float(c_walk)
        o = _try(_nd_system(mu_next), x0_cur, c_cur, t_cur)
        if o is None:
            return None
        x0_cur, t_cur = o.x0, o.period
        mu_cur = mu_next

    return _try(target_system, x0_cur, c_cur, t_cur)


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
