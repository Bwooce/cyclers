"""Circular Restricted N-Body Problem (CRNBP): Jupiter-Europa-Io-Ganymede, N=5.

`#717` (shortlist item 1 of 3 from `#714`'s CCR5BP/CRNBP scoping pass,
``docs/notes/2026-07-27-714-ccr5bp-crnbp-discovery-strategy-pass.md``). This is
the genuine N=5 extension of :mod:`cyclerfinder.core.ccr4bp` (Jupiter-Europa
CR3BP base + Ganymede as a single periodic perturber, N=4): here the base
Jupiter-Europa CR3BP is perturbed by TWO simultaneous extra bodies, Io and
Ganymede, each on its own concentric circular coplanar orbit about Jupiter.
``ccr4bp.py`` itself is UNTOUCHED (per the `#689` discipline) -- this is a new,
separate module that reduces to it exactly at ``mu_io = 0`` (see Gate 1 in
``tests/core/test_crnbp.py``).

Sign reconciliation (`#717` Step 0 -- resolved before any code was written)
----------------------------------------------------------------------------
Two prior digests of Negri & Prado (2022)'s general N-body EOM (their Eq. 11)
transcribed the mutual "inner sum" coupling term among extra bodies with
OPPOSITE signs:

  * ``docs/notes/2026-07-26-digest-negri-prado-2022-crnbp.md`` (§2) wrote
    ``- sum_j mu_j [ direct_j - sum_{k!=j} mu_k (...) ]`` -- an inner MINUS.
  * ``docs/notes/2026-07-26-digest-gilliam-bettinger-2024-crnbp-jovian.md``
    (§1) wrote ``- sum_j mu_j [ direct_j + sum_{k!=j} mu_k (...) ]`` -- an
    inner PLUS.

Resolved directly from the source PDFs (not by guessing or splitting the
difference), independently three times:

  1. Negri & Prado's own arXiv:2307.10881 text (Eqs. 7-9, the pre-simplified
     inertial-frame form, UNAMBIGUOUS -- no multi-line-bracket OCR ambiguity):
     the acceleration on the massless body is
     ``sum_j Mj[(rho_j-rho_N)/|..|^3 + sum_{k!=j} (Mk/(M1+M2))(rho_k-rho_j)/|..|^3]``
     -- direct term and inner-sum coupling term are ADDED (same sign) before
     the whole bracket is multiplied by ``Mj``. Substituting Eq. (10)'s
     ``rho_j - rho_N = -(x+mu2-Rj*cos(psi_j), y-Rj*sin(psi_j), z)`` and
     ``rho_k - rho_j = -(Rj*cos(psi_j)-Rk*cos(psi_k), ...)`` and normalising by
     ``G(M1+M2)=1`` reproduces Eq. (11) with an inner PLUS and an overall minus
     in front of the ``sum_j mu_j[...]`` bracket -- i.e. Gilliam's transcription
     is correct; the ``2026-07-26-digest-negri-prado-2022-crnbp.md`` digest has
     a transcription error (an artefact of that PDF's multi-line bracket
     rendering, which IS genuinely ambiguous under naive ``pdftotext`` -- see
     that digest's own Eq. 11a/b OCR).
  2. Gilliam's own thesis (Eqs. 25-27, a clean text-layer PDF, no bracket
     ambiguity): same inner PLUS, matching (1) term-for-term.
  3. A second, independent source in the SAME thesis PDF (an appended reprint
     of Gilliam, Bettinger et al., *Icarus* 429 (2025) 116455, Eqs. 9-11):
     same inner PLUS again -- a third independent confirmation.

So the correct N-body coupling term (for extra body ``j``, coupling with
every other extra body ``k != j``) is ADDED with the SAME sign as the direct
term, both inside the overall ``-mu_j[...]`` bracket:

    a_j = -mu_j * [ (r - r_j)/|r-r_j|^3 + r_j/a_j^3
                     + sum_{k!=j} mu_k * (r_j - r_k) / |r_j-r_k|^3 ]

(direct + indirect, exactly ``ccr4bp.py``'s single-perturber term, PLUS the new
mutual coupling term, all subtracted from the base CR3BP acceleration in the
sense already established by ``ccr4bp._ganymede_acceleration``).

The load-bearing structural fact (`#714`'s own claim, confirmed here by direct
inspection of the corrected term): the coupling term depends ONLY on
``(R_j, psi_j, R_k, psi_k)`` -- perturber positions, never the spacecraft state
``(x, y, z)``. It is therefore a purely time-periodic forcing addition and
contributes EXACTLY ZERO to the variational (STM) Jacobian: the STM extension
needs only a second copy of ``ccr4bp.py``'s existing direct-term Hessian block
(one per extra perturber), not a new block for the coupling term. This module
confirms that assessment holds (see :func:`_perturber_second_deriv_block`,
which is untouched by the coupling term) -- it does NOT need to be revised.

A further, previously-unflagged structural fact found while implementing this
(`#717` Step 1, going beyond what `#714` checked): the coupling term ALSO
contributes EXACTLY ZERO to the TOTAL spacecraft acceleration, for ANY number
of extra perturbers -- not merely "no Jacobian contribution", but no net
acceleration contribution either. Proof: for any unordered pair of extra
bodies ``{j, k}``, the ordered term ``(j, k)`` contributes
``-mu_j*mu_k*(r_j-r_k)/|r_j-r_k|^3`` to the total acceleration and the ordered
term ``(k, j)`` contributes ``-mu_k*mu_j*(r_k-r_j)/|r_k-r_j|^3 =
+mu_j*mu_k*(r_j-r_k)/|r_j-r_k|^3`` (using ``|r_k-r_j|=|r_j-r_k|`` and
antisymmetry of ``r_k-r_j = -(r_j-r_k)``) -- an EXACT cancellation, for every
pair, regardless of masses or positions. Summing over all pairs (however
many extra perturbers there are) therefore always gives zero. Verified three
ways: (1) this algebraic argument; (2) numerically against
:func:`_perturbers_acceleration` at N=5 (two extra perturbers, physical Io/
Ganymede masses) -- the coupling-included and coupling-omitted accelerations
are bit-identical; (3) a from-scratch, independent re-transcription of
Gilliam's Eq. 25 (not reusing any of this module's code) at BOTH N=5 (two
extra bodies) and N=6 (three extra bodies), same result -- see
``tests/core/test_crnbp.py``'s
``test_coupling_term_net_contribution_is_exactly_zero`` and
``test_coupling_term_net_zero_holds_at_n6_independent_transcription``.

This means the "genuinely new N>=5 physics" `#714` and both digests
anticipated -- and explicitly warned could NOT be captured by a naive
per-perturber superposition ("one cannot just call ``_ganymede_acceleration``
twice") -- turns out, once correctly signed and summed, to be dynamically
INERT for the spacecraft's own motion: a naive superposition of independent
``ccr4bp``-style direct+indirect terms (one per extra perturber, each using
ONLY the primaries M1/M2 in its own "indirect" correction, exactly as
``ccr4bp.py`` already does for Ganymede) is mathematically IDENTICAL to the
full Eq. 11 CRNBP acceleration, for any N. This REFINES (does not invalidate)
`#714`'s tractability assessment: the N=5 extension is a straightforward
SUPERPOSITION of two independent, already-validated single-perturber periodic
forcings (Io's own direct+indirect pull, Ganymede's own -- exactly `#690`'s
existing forcing, unchanged), with no cross term to derive new physics from.
Negri & Prado's own emphasis on the inner sum being "important...to correctly
approximate the indirect effect" (relative to Iuliano's now-superseded
equations) is best read as being about each EXTRA body's own primary-relative
indirect term (the ``k=1,2`` part of its inner sum, which IS what survives
and IS what ``ccr4bp.py``'s existing indirect-term formula already
approximates) -- not about any surviving net effect from extra-body-to-
extra-body mutual coupling, which provably cancels regardless of whether the
formula includes it.

This module STILL implements the coupling term faithfully (rather than
silently omitting it as "known to be zero") for two reasons: (1) fidelity to
the general, arbitrary-N Negri & Prado formula, useful as a reference/
regression against future extensions; (2) it makes the cancellation an
ASSERTED, tested structural property (see the tests above) rather than an
unverified assumption baked silently into a simplified implementation.

Frame / position convention (inherited from ``ccr4bp.py``, not reintroduced)
-----------------------------------------------------------------------------
Negri & Prado's Eq. (10a) defines each extra body's BARYCENTRIC position as
``Rj*(cos psi_j, sin psi_j) - mu2*(1, 0)`` (a ``-mu2`` shift, since body j
orbits M1/Jupiter, which itself sits at ``-mu2`` in barycentric coordinates,
not at the origin). ``ccr4bp.py``'s own ``_ganymede_position`` does NOT apply
this shift -- it places Ganymede directly at ``a_gan*(cos theta, sin theta)``,
already effectively treating Jupiter as coincident with the barycentre. This is
an existing, ALREADY-ACCEPTED small idealisation in ``ccr4bp.py`` (the omitted
shift is O(mu) ~ 2.5e-5 in a_gan ~ 1.6 units, i.e. a ~1.6e-5 relative position
error -- the same order as other approximations ``ccr4bp.py``'s own docstring
already documents, e.g. the <= 2.4e-4 ``GM_Jupiter_sys`` folding). This module
uses the IDENTICAL convention (no ``-mu2`` shift) for BOTH new perturbers
(Io and Ganymede), for two reasons: (1) consistency with the established,
tested ``ccr4bp.py`` behaviour, which this module must reduce to EXACTLY at
``mu_io = 0`` (a byte-exact structural test would fail if the two modules used
different position conventions); (2) the shift is genuinely negligible at this
model's fidelity. Not a new omission -- an explicitly inherited one.

Laplace-projected default (`#714` §1(c))
-----------------------------------------
The Galilean Laplace resonance (``n_Io - 3*n_Europa + 2*n_Ganymede = 0``)
forces the Io:Ganymede synodic-rate ratio in the Jupiter-Europa frame to
EXACTLY ``-2`` at observational precision (checked here from the observed
sidereal periods: ratio ``-2.000000001``, i.e. a ~1.1e-9 relative residual).
The registry-derived two-body rates (:func:`cyclerfinder.core.ccr4bp.two_body_synodic_rate`,
same SMA/GM inputs as ``ccr4bp.py``'s own Ganymede default) give a ratio of
``-1.9996`` -- a ~2.04e-4 relative deviation from the observed-period value,
attributable to registry SMA rounding (NOT real physics: the physical system
enforces the -2 ratio far more tightly than the registry's SMA values alone
would suggest). :func:`jupiter_europa_io_ganymede_default` therefore PROJECTS
Io's synodic rate onto the exact ``omega_io = -2 * omega_gan`` Laplace lock
(mirroring ``ccr4bp.py``'s own documented optional exact-2:1 Ganymede
idealisation -- same class of one-value change, milder because it is backed by
the tighter observed-period justification rather than an arbitrary choice).
Io's initial synodic phase ``theta_io0`` is NOT pinned by the Laplace-libration
argument here (that requires ephemeris mean-longitude epoch data and is
in-scope for `#714` shortlist item 2's continuation task, not this one); it
defaults to 0.0, documented as a free choice for THIS task's equilibrium-point
validation (which does not depend on the exact epoch phase in any way that
matters to the Newton-Raphson equilibrium solve below).

Positive control (`#717` Step 2)
---------------------------------
Gilliam's thesis Ch. V computes CRNBP-perturbed "dynamical substitute"
equilibrium points (Newton-Raphson force-balance loci, NOT periodic orbits) for
several systems, including Jupiter-Europa with Io+Ganymede as the two extra
perturbers -- exactly this module's N=5 configuration (thesis §5.3.1: "This is
true in the case of the Laplace resonance between Io, Europa, and Ganymede...
Figure 41 demonstrates the E1 and E2 dynamical substitutes... extended from the
Jupiter-Europa CR3BP to include Io and Ganymede"). Her Tables 5-6 give
digit-grade ``dx``/``dy``/Lagrange-Box-area/avg-vs-CR3BP numbers for the E1/E2
points. This module's :func:`equilibrium_dynamical_substitute` reproduces her
Newton-Raphson method (Ch. III Eqs. 57-69: 2D force-balance ``F(x,y) = (xddot,
yddot) = 0`` at fixed perturber phases, Jacobian = the SAME direct-term Hessian
block the STM uses, continuation-seeded across a phase sweep). This is an
EOM/Jacobian IMPLEMENTATION control only (per `#714`'s own tiering: Tier 1,
"digit-grade, literature-sourced, ALGEBRAIC only" -- static force-balance loci,
not periodic orbits, not connections); see ``tests/core/test_crnbp.py`` for the
comparison and its documented precision limits (our registry-sourced mu3/mu4
values differ from Gilliam's own stated mu3/mu4 at the ~3-5e-4 relative level
-- both are traceable to the SAME physical Jovian moon masses/SMAs but through
different mass-vs-GM rounding conventions; see that test module's docstring for
the full reconciliation).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.ccr4bp import _cr3bp_uxx_block, two_body_synodic_rate
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES


@dataclass(frozen=True)
class CRNBPPerturber:
    """One extra (beyond Jupiter+Europa) perturbing body on a concentric circle.

    All quantities in Jupiter-Europa synodic nondimensional units (matching
    :class:`cyclerfinder.core.ccr4bp.CCR4BPSystem` exactly).

    Attributes
    ----------
    mu :
        Nondimensional mass ``M_j/(M1+M2)`` (coefficient of the direct +
        indirect + coupling acceleration). Set to 0.0 to remove this
        perturber's contribution entirely (used by the structural reduction
        tests).
    a :
        Orbital radius about Jupiter, Europa-SMA units (> 0).
    omega :
        Synodic angular rate ``n_j/n2 - 1`` in the Europa-synodic frame.
    theta0 :
        Synodic phase at ``t = 0``, rad.
    """

    mu: float
    a: float
    omega: float
    theta0: float = 0.0

    def position(self, t: float) -> tuple[float, float]:
        """This perturber's ``(x, y)`` position in the synodic frame at time ``t``."""
        theta = self.theta0 + self.omega * t
        return self.a * math.cos(theta), self.a * math.sin(theta)

    @property
    def synodic_period(self) -> float:
        """``2*pi / |omega|`` (TU); shared forcing clock only if commensurate
        with the other perturbers' rates."""
        return 2.0 * math.pi / abs(self.omega)


@dataclass(frozen=True)
class CRNBPSystem:
    """Circular Restricted N-Body Problem parameters: Jupiter-Europa base +
    an arbitrary tuple of extra concentric-circular perturbers.

    At ``perturbers = ()`` this is exactly the Jupiter-Europa CR3BP. At a
    single perturber it is exactly :class:`cyclerfinder.core.ccr4bp.CCR4BPSystem`
    (see the ``mu_io -> 0`` reduction test). At two or more perturbers the
    pairwise inner-sum coupling term (`#717` Step 0) is present in the formula
    and correctly signed here, but -- a further finding of this task, see the
    module docstring's cancellation-proof section -- it contributes EXACTLY
    ZERO to the total spacecraft acceleration, for any number of perturbers.
    The genuinely new N>=5 physics is therefore simply the SUPERPOSITION of
    each perturber's own independent direct+indirect forcing, not a mutual
    coupling effect.
    """

    mu: float
    perturbers: tuple[CRNBPPerturber, ...]


def jupiter_europa_io_ganymede_default() -> CRNBPSystem:
    """CRNBP parameters for Jupiter-Europa-Io-Ganymede (N=5) from the JPL SSD
    registry, mirroring :func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default`.

    Reuses the SAME sourced ``core.satellites`` registry ``ccr4bp.py`` uses (no
    independently re-sourced constants): Europa's base mu is IDENTICAL to
    ``ccr4bp.jupiter_europa_ganymede_default().mu``; Ganymede's ``(mu, a,
    omega)`` triple is likewise identical to that function's ``(mu_gan, a_gan,
    omega_gan)``. Io is added as the second perturber with the SAME formula
    pattern: ``mu_io = GM_Io / (GM_Jupiter_sys + GM_Europa)``, ``a_io =
    SMA_Io / SMA_Europa`` (~0.6285), ``omega_io`` LAPLACE-PROJECTED onto the
    exact ``-2 * omega_gan`` ratio (see the module docstring's "Laplace-
    projected default" section for the ~2.04e-4-vs-1.1e-9 justification).
    ``theta_io0 = 0.0`` (a documented free choice, NOT epoch-pinned to the
    Laplace libration argument -- see module docstring).
    """
    gm_j = PRIMARIES["Jupiter"]
    gm_e = SATELLITES["Europa"].mu_km3_s2
    gm_io = SATELLITES["Io"].mu_km3_s2
    gm_gan = SATELLITES["Ganymede"].mu_km3_s2
    sma_e = SATELLITES["Europa"].sma_km
    sma_io = SATELLITES["Io"].sma_km
    sma_gan = SATELLITES["Ganymede"].sma_km
    denom = gm_j + gm_e  # G(m1+m2), matches ccr4bp.jupiter_europa_ganymede_default exactly

    mu = gm_e / denom
    mu_io = gm_io / denom
    mu_gan = gm_gan / denom
    a_io = sma_io / sma_e
    a_gan = sma_gan / sma_e
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)
    omega_io = -2.0 * omega_gan  # exact Laplace-lock projection

    io = CRNBPPerturber(mu=mu_io, a=a_io, omega=omega_io, theta0=0.0)
    ganymede = CRNBPPerturber(mu=mu_gan, a=a_gan, omega=omega_gan, theta0=0.0)
    return CRNBPSystem(mu=mu, perturbers=(io, ganymede))


def _perturber_positions(t: float, system: CRNBPSystem) -> list[tuple[float, float, float]]:
    """``(mu_j, gx_j, gy_j)`` for every perturber at time ``t`` (mass included
    so zero-mass perturbers are cheap to skip downstream)."""
    return [(p.mu, *p.position(t)) for p in system.perturbers]


def _coupling_term_for_perturber(
    j: int, positions: list[tuple[float, float, float]]
) -> tuple[float, float]:
    """The ``(x, y)`` mutual-coupling contribution for perturber ``j`` ALONE
    (its own bracket's inner sum over every OTHER perturber ``k``), per the
    corrected Eq. 11: ``-mu_j * sum_{k!=j} mu_k * (r_j - r_k) / |r_j-r_k|^3``.

    Exposed separately from :func:`_perturbers_acceleration` (which sums this
    over every ``j``) so the per-body term's SIGN can be tested in isolation
    -- summing it over all ``j`` provably cancels to exactly zero (see the
    module docstring), which makes the TOTAL an unsuitable place to detect a
    sign error (a sign flip here still sums to zero -- verified in
    ``tests/core/test_crnbp.py``). This isolated per-``j`` value is what
    actually carries the sign information.
    """
    mu_j, gxj, gyj = positions[j]
    ax = ay = 0.0
    if mu_j == 0.0:
        return 0.0, 0.0
    for k in range(len(positions)):
        if k == j:
            continue
        mu_k, gxk, gyk = positions[k]
        if mu_k == 0.0:
            continue
        rkjx = gxj - gxk
        rkjy = gyj - gyk
        rkj2 = rkjx * rkjx + rkjy * rkjy
        rkj3 = rkj2 * math.sqrt(rkj2)
        ax += -mu_j * mu_k * rkjx / rkj3
        ay += -mu_j * mu_k * rkjy / rkj3
    return ax, ay


def _perturbers_acceleration(
    x: float, y: float, z: float, t: float, system: CRNBPSystem
) -> tuple[float, float, float]:
    """Direct + indirect + mutual-coupling acceleration from ALL perturbers.

    Direct + indirect: identical in form to
    :func:`cyclerfinder.core.ccr4bp._ganymede_acceleration`, summed over every
    perturber. Coupling (`#717` Step 0's corrected term; see
    :func:`_coupling_term_for_perturber`): in the (x, y) plane only
    (perturbers are coplanar; no z-coupling term exists in Eq. 11c/Eq. 27).
    Spacecraft-state-independent (depends only on perturber positions), so
    contributes nothing to the Jacobian -- see
    :func:`_perturber_second_deriv_block`, which has no coupling analogue.
    Also contributes EXACTLY ZERO to the TOTAL (summed-over-j) acceleration
    here -- see the module docstring's cancellation proof -- but is computed
    faithfully term-by-term regardless (fidelity to the general formula, and
    to make the cancellation a TESTED property rather than a silent
    assumption).
    """
    positions = _perturber_positions(t, system)
    n = len(positions)
    ax = ay = az = 0.0
    for j in range(n):
        mu_j, gxj, gyj = positions[j]
        if mu_j == 0.0:
            continue
        a_j = system.perturbers[j].a
        dx = x - gxj
        dy = y - gyj
        dz = z
        d2 = dx * dx + dy * dy + dz * dz
        d3 = d2 * math.sqrt(d2)
        a_j3 = a_j**3
        ax += -mu_j * dx / d3 - mu_j * gxj / a_j3
        ay += -mu_j * dy / d3 - mu_j * gyj / a_j3
        az += -mu_j * dz / d3
        cx, cy = _coupling_term_for_perturber(j, positions)
        ax += cx
        ay += cy
        # No z-coupling term: perturbers are coplanar (Eq. 11c has no
        # coupling sum at all).
    return ax, ay, az


def crnbp_eom(t: float, state6: NDArray[np.float64], system: CRNBPSystem) -> NDArray[np.float64]:
    """CRNBP equations of motion in the Jupiter-Europa synodic frame.

    ``state6 = (x, y, z, vx, vy, vz)``. Returns d(state)/dt.

    Reduces to :func:`cyclerfinder.core.ccr4bp.ccr4bp_eom` EXACTLY when
    ``system.perturbers`` has exactly one nonzero-mass entry matching
    Ganymede's ``(mu_gan, a_gan, omega_gan, theta_gan0)`` and any other
    perturbers have ``mu = 0`` (see ``tests/core/test_crnbp.py`` Gate 1).
    Reduces to :func:`cyclerfinder.core.cr3bp.cr3bp_eom` when ALL perturbers
    have ``mu = 0`` (or ``system.perturbers = ()``).
    """
    x, y, z, vx, vy, vz = (float(v) for v in state6)
    base = cr3bp.cr3bp_eom(t, np.asarray(state6, dtype=np.float64), system.mu)
    ax_p, ay_p, az_p = _perturbers_acceleration(x, y, z, t, system)
    return np.array(
        [vx, vy, vz, float(base[3]) + ax_p, float(base[4]) + ay_p, float(base[5]) + az_p],
        dtype=np.float64,
    )


def _perturber_second_deriv_block(
    x: float, y: float, z: float, t: float, system: CRNBPSystem
) -> tuple[float, float, float, float, float, float]:
    """Second derivatives of the perturbers' DIRECT acceleration w.r.t. (x, y, z).

    Sum, over every perturber, of the SAME single-body Hessian block
    :func:`cyclerfinder.core.ccr4bp._ganymede_second_deriv_block` computes.
    Confirms `#714`'s structural claim: the coupling term (indirect too) is
    spacecraft-state-independent, so it contributes NOTHING here -- the STM
    extension is exactly "one more copy of the existing direct-term block per
    extra perturber", no new term.
    """
    sxx = syy = szz = sxy = sxz = syz = 0.0
    for p in system.perturbers:
        if p.mu == 0.0:
            continue
        gx, gy = p.position(t)
        dx = x - gx
        dy = y - gy
        dz = z
        d2 = dx * dx + dy * dy + dz * dz
        d3 = d2 * math.sqrt(d2)
        d5 = d3 * d2
        inv_d3 = 1.0 / d3
        inv_d5 = 1.0 / d5
        sxx += -p.mu * (inv_d3 - 3.0 * dx * dx * inv_d5)
        syy += -p.mu * (inv_d3 - 3.0 * dy * dy * inv_d5)
        szz += -p.mu * (inv_d3 - 3.0 * dz * dz * inv_d5)
        sxy += -p.mu * (-3.0 * dx * dy * inv_d5)
        sxz += -p.mu * (-3.0 * dx * dz * inv_d5)
        syz += -p.mu * (-3.0 * dy * dz * inv_d5)
    return sxx, syy, szz, sxy, sxz, syz


def crnbp_stm_eom(
    t: float, state_and_stm: NDArray[np.float64], system: CRNBPSystem
) -> NDArray[np.float64]:
    """CRNBP variational EOM: state (6) + flattened 6x6 STM (36).

    A-matrix layout matches CR3BP / bcr4bp / ccr4bp exactly: identity on the
    upper-right 3x3, Coriolis on the lower-middle, pseudo-potential second
    derivatives (base CR3BP + every perturber's direct-term Hessian block,
    per :func:`_perturber_second_deriv_block`) on the lower-left. Reduces to
    :func:`cyclerfinder.core.ccr4bp.ccr4bp_stm_eom` under the same conditions
    as :func:`crnbp_eom`.
    """
    s = state_and_stm[:6]
    phi = state_and_stm[6:].reshape(6, 6)
    x, y, z = float(s[0]), float(s[1]), float(s[2])
    uxx, uyy, uzz, uxy, uxz, uyz = _cr3bp_uxx_block(x, y, z, system.mu)
    pxx, pyy, pzz, pxy, pxz, pyz = _perturber_second_deriv_block(x, y, z, t, system)
    uxx += pxx
    uyy += pyy
    uzz += pzz
    uxy += pxy
    uxz += pxz
    uyz += pyz
    mat_a = np.zeros((6, 6))
    mat_a[0:3, 3:6] = np.eye(3)
    mat_a[3, 0], mat_a[3, 1], mat_a[3, 2] = uxx, uxy, uxz
    mat_a[4, 0], mat_a[4, 1], mat_a[4, 2] = uxy, uyy, uyz
    mat_a[5, 0], mat_a[5, 1], mat_a[5, 2] = uxz, uyz, uzz
    mat_a[3, 4], mat_a[4, 3] = 2.0, -2.0
    dphi = mat_a @ phi
    return np.concatenate([crnbp_eom(t, s, system), dphi.reshape(36)])


@dataclass(frozen=True)
class CRNBPArc:
    """Result of a CRNBP propagation. Mirrors :class:`cyclerfinder.core.ccr4bp.CCR4BPArc`."""

    state_f: NDArray[np.float64]
    stm: NDArray[np.float64] | None
    t: float


def propagate_crnbp(
    system: CRNBPSystem,
    state6: NDArray[np.float64],
    t: float,
    *,
    with_stm: bool = False,
    t0: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> CRNBPArc:
    """Propagate a CRNBP state (and optionally the STM) from ``t0`` to ``t0+t``.

    Mirrors :func:`cyclerfinder.core.ccr4bp.propagate_ccr4bp` (DOP853).

    Raises
    ------
    RuntimeError
        If the integrator fails.
    """
    state0 = np.asarray(state6, dtype=np.float64)
    if with_stm:
        y0 = np.concatenate([state0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            crnbp_stm_eom,
            (t0, t0 + t),
            y0,
            args=(system,),
            rtol=rtol,
            atol=atol,
            method="DOP853",
        )
        if not sol.success:
            raise RuntimeError(f"CRNBP STM propagation failed at t={sol.t[-1]}: {sol.message}")
        yf = sol.y[:, -1]
        return CRNBPArc(state_f=yf[:6], stm=yf[6:].reshape(6, 6), t=float(sol.t[-1]) - t0)
    sol = solve_ivp(
        crnbp_eom,
        (t0, t0 + t),
        state0,
        args=(system,),
        rtol=rtol,
        atol=atol,
        method="DOP853",
    )
    if not sol.success:
        raise RuntimeError(f"CRNBP propagation failed at t={sol.t[-1]}: {sol.message}")
    return CRNBPArc(state_f=sol.y[:, -1], stm=None, t=float(sol.t[-1]) - t0)


# ---------------------------------------------------------------------------
# Equilibrium-point ("dynamical substitute") Newton-Raphson solver (`#717`
# Step 2's positive control -- reproduces Gilliam's thesis Ch. III/V method).
# ---------------------------------------------------------------------------


def cr3bp_collinear_point(mu: float, which: str) -> float:
    """Base CR3BP collinear libration point x-coordinate (y=z=0), 1D Newton-Raphson.

    ``which`` in {"L1", "L2", "L3"}. Standard pseudo-potential-gradient root
    (``Omega_x(x, 0, 0) = 0``); used only as the initial guess for the CRNBP
    equilibrium Newton-Raphson below, mirroring Gilliam's own method ("the
    initial guess is the corresponding CR3BP Lagrange point").
    """

    def omega_x(x: float) -> float:
        r1 = abs(x + mu)
        r2 = abs(x - 1.0 + mu)
        return x - (1.0 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3

    if which == "L1":
        lo, hi = 1.0 - mu - 0.5, 1.0 - mu - 1e-9
    elif which == "L2":
        lo, hi = 1.0 - mu + 1e-9, 1.0 - mu + 0.5
    elif which == "L3":
        lo, hi = -2.0, -mu - 1e-9
    else:
        raise ValueError(f"cr3bp_collinear_point: unknown point {which!r}")
    from scipy.optimize import brentq

    return float(brentq(omega_x, lo, hi, xtol=1e-14, rtol=1e-14))


def _equilibrium_residual_and_jacobian(
    x: float, y: float, t: float, system: CRNBPSystem
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """``F(x,y) = (xddot, yddot)`` at ``(x, y, 0, 0, 0, 0)`` and time ``t``,
    plus its 2x2 Jacobian -- the SAME force field and Hessian block
    :func:`crnbp_eom` / :func:`crnbp_stm_eom` already compute, reused
    verbatim (this is what makes the equilibrium solve a joint validation of
    both the EOM and the STM's direct-term Hessian block, per `#714` §1(b)).
    """
    state0 = np.array([x, y, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    rhs = crnbp_eom(t, state0, system)
    f = np.array([rhs[3], rhs[4]], dtype=np.float64)
    uxx, uyy, _uzz, uxy, _uxz, _uyz = _cr3bp_uxx_block(x, y, 0.0, system.mu)
    pxx, pyy, _pzz, pxy, _pxz, _pyz = _perturber_second_deriv_block(x, y, 0.0, t, system)
    jac = np.array([[uxx + pxx, uxy + pxy], [uxy + pxy, uyy + pyy]], dtype=np.float64)
    return f, jac


def solve_equilibrium_point(
    x0: float, y0: float, t: float, system: CRNBPSystem, *, max_iter: int = 50, tol: float = 1e-13
) -> tuple[float, float]:
    """Multi-dimensional Newton-Raphson equilibrium ("dynamical substitute")
    solve at fixed time ``t`` (fixed perturber phases): ``(xddot, yddot) = 0``.

    Mirrors Gilliam's thesis Eqs. 57-69 exactly (2D force-balance, analytic
    Jacobian, ``x_{k+1} = x_k - J^-1 F(x_k)``). Raises ``RuntimeError`` if it
    fails to converge within ``max_iter``.
    """
    x, y = x0, y0
    for _ in range(max_iter):
        f, jac = _equilibrium_residual_and_jacobian(x, y, t, system)
        if float(np.max(np.abs(f))) < tol:
            return x, y
        step = np.linalg.solve(jac, f)
        x -= float(step[0])
        y -= float(step[1])
    f, _ = _equilibrium_residual_and_jacobian(x, y, t, system)
    if float(np.max(np.abs(f))) >= tol:
        raise RuntimeError(
            f"solve_equilibrium_point: failed to converge from ({x0}, {y0}) at t={t}: "
            f"final residual {f}"
        )
    return x, y


@dataclass(frozen=True)
class DynamicalSubstitute:
    """A CRNBP equilibrium-point "dynamical substitute" arc and its Lagrange-Box
    metrics, mirroring Gilliam's thesis Tables 5-6 columns exactly.

    ``xs``/``ys`` (nondim) are the continuation-seeded equilibrium solution at
    each sampled time in ``ts``. ``dx_km``/``dy_km`` are the peak-to-peak
    excursions (max - min) over the window, in km (``l_km`` scale).
    ``box_area_km2 = dx_km * dy_km``. ``avg_vs_cr3bp_km`` is the DISTANCE
    between the TIME-AVERAGED substitute location (mean of ``xs``, ``ys`` --
    Gilliam's "blue asterisk", her thesis Fig. 41's own description) and the
    base CR3BP point, in km -- per her thesis's own wording ("the average
    location of the equilibrium point ... is compared to the corresponding
    calculated CR3BP Lagrange point"), NOT the mean of the pointwise
    displacements (that alternative reading was checked and does not
    reproduce her numbers; this one does, to a few percent -- see
    ``tests/core/test_crnbp.py``).
    """

    xs: NDArray[np.float64]
    ys: NDArray[np.float64]
    ts: NDArray[np.float64]
    x_cr3bp: float
    y_cr3bp: float
    dx_km: float
    dy_km: float
    box_area_km2: float
    avg_vs_cr3bp_km: float


def equilibrium_dynamical_substitute(
    system: CRNBPSystem,
    which: str,
    *,
    l_km: float,
    t_window: float,
    n_samples: int = 1000,
) -> DynamicalSubstitute:
    """CRNBP equilibrium dynamical substitute for a collinear CR3BP point,
    propagated over ``[0, t_window]`` (nondim TU), continuation-seeded.

    ``which`` in {"L1", "L2"} (E1/E2 in Gilliam's notation -- her Ch. V never
    reports E3/E4/E5 digit-grade tables). ``l_km`` is the DU->km scale (this
    project's own established Jupiter-Europa convention, 671100.0 km -- see
    ``tests/core/test_ccr4bp.py``'s ``l_km=671100.0``, matching this thesis's
    own Table 14 Europa SMA exactly).
    """
    x_cr3bp = cr3bp_collinear_point(system.mu, which)
    y_cr3bp = 0.0
    ts = np.linspace(0.0, t_window, n_samples)
    xs = np.empty(n_samples, dtype=np.float64)
    ys = np.empty(n_samples, dtype=np.float64)
    x, y = x_cr3bp, y_cr3bp
    for i, t in enumerate(ts):
        x, y = solve_equilibrium_point(x, y, float(t), system)
        xs[i] = x
        ys[i] = y
    dx_km = float(np.max(xs) - np.min(xs)) * l_km
    dy_km = float(np.max(ys) - np.min(ys)) * l_km
    box_area_km2 = dx_km * dy_km
    x_avg, y_avg = float(np.mean(xs)), float(np.mean(ys))
    avg_vs_cr3bp_km = float(math.hypot(x_avg - x_cr3bp, y_avg - y_cr3bp) * l_km)
    return DynamicalSubstitute(
        xs=xs,
        ys=ys,
        ts=ts,
        x_cr3bp=x_cr3bp,
        y_cr3bp=y_cr3bp,
        dx_km=dx_km,
        dy_km=dy_km,
        box_area_km2=box_area_km2,
        avg_vs_cr3bp_km=avg_vs_cr3bp_km,
    )
