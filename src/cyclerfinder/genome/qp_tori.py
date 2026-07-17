"""Quasi-periodic invariant 2-tori in the CR3BP (#290 Phase 1).

A quasi-periodic invariant 2-torus is a 2-dimensional invariant manifold of the
CR3BP that supports a 2-frequency motion: trajectories on the torus wind around
both angles ``(theta_long, theta_trans)`` with frequencies ``(omega_long,
omega_trans)`` whose ratio is irrational (otherwise it phase-locks to a
periodic orbit).

QP-tori are born at **Neimark-Sacker bifurcations** of a parent periodic
family: where a complex-conjugate pair of Floquet multipliers crosses the unit
circle off the real axis at a primitive k-th root of unity (k >= 3 generically,
or simply ``|lam|=1, arg lam not= 0, pi``). Near the bifurcation the torus is
a thin tube around the parent orbit; far from it the torus thickens and
eventually breaks up (Arnol'd tongues / KAM cantori).

GMOS (Gomez-Mondelo-Olikara-Scheeres) parameterization
------------------------------------------------------

Following Olikara-Scheeres 2010 and Olikara 2016 (Purdue PhD), we work with
the **invariant circle** representation rather than the full 2D Fourier
expansion. Pick a longitudinal stroboscopic time::

    t_strob = 2 pi / omega_long      (one longitudinal period)

Restrict the torus to the section ``theta_long = 0``. This section is a
topological circle parameterized by ``theta_trans in [0, 2 pi)``, the
**invariant circle**::

    u(theta) := x(0, theta)

Represent ``u(theta)`` as a truncated Fourier series in ``theta``::

    u(theta) = sum_{n=-N}^{N} c_n exp(i n theta)        c_n in C^6, c_{-n} = conj(c_n)

so the *real* state is recovered from ``2 N + 1`` complex Fourier modes whose
total real DOF is ``6 * (2 N + 1)`` (with the reality constraint halving the
unique DOF to ``6 * (N + 1)`` -- the constant mode ``c_0`` is real, and modes
``c_1 ... c_N`` are complex).

**Invariance equation**: the stroboscopic flow ``phi_{t_strob}`` advances the
longitudinal angle by ``2 pi`` and the transverse angle by the *rotation
number*::

    rho := t_strob * omega_trans = 2 pi * omega_trans / omega_long

so points on the invariant circle map to other points on the invariant circle,
shifted by ``rho``::

    phi_{t_strob}(u(theta)) = u(theta + rho)            for all theta

In Fourier-mode space this gives the GMOS residual equation::

    F_n(c, rho, t_strob) := FFT[phi_{t_strob}(u_j)]_n - c_n * exp(i n rho) = 0

for ``n = -N .. N``, evaluated on ``M = 2 N + 1`` (or more) sample points
``theta_j = 2 pi j / M``, ``u_j = u(theta_j)``.

Free variables: the Fourier modes ``{c_n}`` (with ``c_0`` real and ``c_n``
complex for ``n > 0``), the rotation number ``rho``, and the stroboscopic time
``t_strob`` (equivalently the longitudinal frequency ``omega_long = 2 pi /
t_strob``). The total real DOF is ``6 * (2 N + 1) + 2``; the residual has the
same dimension once two phase / amplitude conditions are imposed (see
``correct_qp_torus`` below).

**Independent cross-check** (orbit-closure discipline, see
``feedback_orbit_closure_discipline``): after Newton convergence, draw random
sample angles ``theta_check`` NOT on the FFT grid, propagate
``u(theta_check)`` by ``t_strob``, and verify the result lies on the
*resampled* invariant circle to within ``independent_tol`` (default 1e-4 --
QP-tori are inherently noisier than strict-periodic orbits because of Fourier
truncation; see Olikara 2016 for the truncation-error analysis).

Discipline
----------

* Returned torus is OUR computation; no novelty claim is made. The method is
  Olikara-Scheeres 2010 / Olikara 2016.
* No catalogue writeback. Catalogue admission of a QP-torus as a quasi-cycler
  requires V0-V5 gauntlet adaptation -- a future Phase, not Phase 1.
* The seed (Neimark-Sacker bracket from #299) is sourced via the parent
  family's literature check (Antoniadou-Voyatzis 2018 / Roberts-Tsoukkas-Ross
  2026 sourcing chain).

References
----------

* Olikara, Z., & Scheeres, D. (2010). "Numerical Method for Computing
  Quasi-Periodic Orbits and Their Stability in the Restricted Three-Body
  Problem." AAS Astrodynamics Specialist Conference.
* Olikara, Z. (2016). "Computation of Quasi-Periodic Tori and Heteroclinic
  Connections in Astrodynamics." PhD dissertation, Purdue University.
* Howell, K. C., & Howell, A. R. (2014). "Survey of Quasi-Periodic Motion in
  Cislunar Space for Transfer Design."
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp

# ---------------------------------------------------------------------------
# Data model.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class QPTorus:
    """Quasi-periodic invariant 2-torus in CR3BP / 3D phase space.

    Parameterized as an invariant circle ``u(theta)`` on the section
    ``theta_long = 0``, represented by truncated 1D Fourier coefficients of the
    six state-space coordinates. The full torus is the image of this circle
    under the stroboscopic flow with rotation number ``rho``.

    Attributes
    ----------
    system :
        The CR3BP system (provides ``mu``).
    omega_long :
        Longitudinal frequency (rad/TU). ``omega_long = 2 pi / t_strob``.
    omega_trans :
        Transverse frequency (rad/TU). ``omega_trans = rho * omega_long / (2 pi)``.
    rho :
        Rotation number on the invariant circle (radians per stroboscopic
        period). Equivalent to ``2 pi * omega_trans / omega_long``.
    t_strob :
        Stroboscopic time (nondim TU). ``t_strob = 2 pi / omega_long``. This is
        the longitudinal-period analogue.
    fourier_coeffs :
        Complex Fourier coefficients of the invariant circle, shape
        ``(2 N + 1, 6)``. Layout: index ``0`` -> mode ``n=0``, indices
        ``1 ... N`` -> modes ``n=1 ... N``, indices ``N+1 ... 2N`` -> modes
        ``n=-N ... -1`` (matches ``numpy.fft.fft`` convention). The reality
        constraint ``c_{-n} = conj(c_n)`` is enforced after every Newton step.
    n_modes :
        ``N`` -- the number of POSITIVE modes (``2 N + 1`` total including
        ``c_0`` and the negative modes).
    n_samples :
        ``M`` -- the number of sample points on the invariant circle used to
        evaluate the GMOS residual. Must be ``>= 2 N + 1``; defaults to
        ``2 N + 1``.
    invariance_residual :
        L2 norm of the GMOS residual in Fourier-mode space, after Newton
        convergence. Closure gate (compare to ``tol``).
    independent_closure_residual :
        Independent cross-check: maximum distance, over a set of OFF-GRID
        sample angles, between the propagated invariant-circle point and the
        nearest point on the resampled invariant circle. Compare to
        ``independent_tol`` (default 1e-4).
    converged :
        ``True`` iff Newton converged AND the independent closure check
        passed.
    n_iter :
        Newton iterations consumed.
    notes :
        Free-text caller annotation.
    """

    system: cr3bp.CR3BPSystem
    omega_long: float
    omega_trans: float
    rho: float
    t_strob: float
    fourier_coeffs: NDArray[np.complex128]
    n_modes: int
    n_samples: int
    invariance_residual: float
    independent_closure_residual: float
    converged: bool
    n_iter: int
    notes: str = ""
    extras: dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Torus evaluation: inverse Fourier transform.
# ---------------------------------------------------------------------------


def evaluate_invariant_circle(
    fourier_coeffs: NDArray[np.complex128],
    theta: float | NDArray[np.float64],
) -> NDArray[np.float64]:
    """Return ``u(theta)`` = sum_n c_n exp(i n theta) (real part) for each
    coordinate.

    Parameters
    ----------
    fourier_coeffs :
        Shape ``(2 N + 1, 6)``, numpy-FFT ordering: ``[c_0, c_1, ..., c_N,
        c_{-N}, ..., c_{-1}]``.
    theta :
        Scalar or 1D array of angles in radians.

    Returns
    -------
    NDArray :
        Shape ``(6,)`` if ``theta`` is scalar, ``(len(theta), 6)`` if array.
    """
    n_total = fourier_coeffs.shape[0]
    n_pos = (n_total - 1) // 2
    theta_arr = np.atleast_1d(np.asarray(theta, dtype=np.float64))
    out = np.zeros((theta_arr.size, 6), dtype=np.complex128)
    # mode n=0
    out += fourier_coeffs[0:1, :]
    # modes n = 1 .. n_pos and their conjugates at -n
    for n in range(1, n_pos + 1):
        c_pos = fourier_coeffs[n, :]
        c_neg = fourier_coeffs[n_total - n, :]  # numpy-FFT layout
        out += np.outer(np.exp(1j * n * theta_arr), c_pos)
        out += np.outer(np.exp(-1j * n * theta_arr), c_neg)
    real_out = np.real(out).astype(np.float64)
    if np.isscalar(theta):
        return np.asarray(real_out[0], dtype=np.float64)
    return np.asarray(real_out, dtype=np.float64)


def evaluate_torus(
    torus: QPTorus,
    theta_long: float,
    theta_trans: float,
) -> NDArray[np.float64]:
    """Return the state on the torus at ``(theta_long, theta_trans)``.

    Uses the stroboscopic-flow definition::

        x(theta_long, theta_trans) = phi_{theta_long / omega_long}(
            u(theta_trans - rho * theta_long / (2 pi))
        )

    i.e. flow the invariant circle forward by a fraction of the stroboscopic
    period, AFTER pre-rotating the circle by the conjugate rotation so the
    forward flow exactly cancels at ``theta_long = 0`` (returning ``u(theta_trans)``).

    Parameters
    ----------
    torus :
        A converged QPTorus.
    theta_long, theta_trans :
        Angles in radians (will be reduced mod ``2 pi``).
    """
    theta_long_red = float(theta_long) % (2 * math.pi)
    theta_trans_red = float(theta_trans) % (2 * math.pi)
    dt = theta_long_red / torus.omega_long
    # Pre-rotate the circle so the forward flow lands at theta_trans at the
    # target longitudinal angle.
    theta_seed = (theta_trans_red - torus.rho * theta_long_red / (2 * math.pi)) % (2 * math.pi)
    u0 = evaluate_invariant_circle(torus.fourier_coeffs, theta_seed)
    if dt == 0.0:
        return u0
    arc = cr3bp.propagate(torus.system, u0, dt, with_stm=False)
    return arc.state_f


# ---------------------------------------------------------------------------
# GMOS residual + Newton corrector.
# ---------------------------------------------------------------------------


def _enforce_reality(coeffs: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Project Fourier-mode coefficients onto the real-valued-signal subspace:
    ``c_{-n} = conj(c_n)``, ``c_0`` real. This is one Newton-step-end clean-up.
    """
    n_total = coeffs.shape[0]
    n_pos = (n_total - 1) // 2
    out = coeffs.copy()
    out[0, :] = np.real(out[0, :])
    for n in range(1, n_pos + 1):
        avg = 0.5 * (out[n, :] + np.conj(out[n_total - n, :]))
        out[n, :] = avg
        out[n_total - n, :] = np.conj(avg)
    return out


def _seed_invariant_circle(
    parent_state: NDArray[np.float64],
    monodromy_matrix: NDArray[np.float64],
    *,
    k: int,
    n_modes: int,
    amplitude: float,
) -> tuple[NDArray[np.complex128], complex, NDArray[np.complex128]]:
    """Seed an invariant circle from the Neimark-Sacker eigenvector of the
    parent's monodromy.

    The Neimark-Sacker eigenvector pair ``v, conj(v)`` with eigenvalue
    ``lam = e^(i*phi)`` (``|lam|=1``, ``phi != 0, pi``) defines the local 2D
    invariant subspace tangent to the emerging torus at the parent orbit. The
    leading-order seed for the invariant circle is::

        u_seed(theta) ~= s_parent + amplitude * (Re(v) cos theta - Im(v) sin theta)

    i.e. the constant mode is ``s_parent`` and the single n=1 mode is
    ``amplitude * (Re(v) + i Im(v)) / 2`` (so 2 * Re(c_1 e^(i theta)) gives the
    cos/sin combination above).

    Returns
    -------
    (coeffs, lam, v_complex) :
        The seed Fourier coefficients ``(2 n_modes + 1, 6)``, the eigenvalue,
        and the (complex, normalized) eigenvector.
    """
    eigvals, eigvecs = np.linalg.eig(monodromy_matrix)
    # Find the Neimark-Sacker pair: |lam| ~ 1, lam not real, and arg(lam)
    # closest to a primitive k-th root of unity.
    target_args = []
    for j in range(1, k):
        if math.gcd(j, k) == 1:
            target_args.append(2 * math.pi * j / k)
    if not target_args:
        raise ValueError(f"no primitive k-th roots of unity for k={k}")
    best = None
    best_dist = float("inf")
    for i, lam in enumerate(eigvals):
        if abs(abs(lam) - 1.0) > 0.1:
            continue
        if abs(np.imag(lam)) < 1e-6:
            continue  # real eigenvalue -- saddle-node or period-doubling, not NS
        phi = math.atan2(np.imag(lam), np.real(lam))
        for ta in target_args:
            d = min(abs(phi - ta), abs(phi - ta + 2 * math.pi), abs(phi - ta - 2 * math.pi))
            if d < best_dist:
                best_dist = d
                best = (i, complex(lam))
    if best is None or best_dist > 0.5:
        raise ValueError(
            f"no Neimark-Sacker eigenvalue near primitive k={k}-th root of unity; "
            f"best distance was {best_dist:.4f}"
        )
    idx, lam = best
    v_complex = eigvecs[:, idx].astype(np.complex128)
    # Canonicalize the Neimark-Sacker representative to the POSITIVE-imaginary
    # member of the conjugate pair (lam, conj(lam)). The two members describe
    # the SAME invariant circle traversed in opposite angular senses
    # (theta -> -theta), so the sign of the rotation number rho = arg(lam) is a
    # pure parametrization convention. But WHICH member ``np.linalg.eig``
    # returns first is BLAS-backend dependent (Mac/Accelerate vs Linux/OpenBLAS
    # order the conjugate pair differently), so without pinning it the reported
    # rotation-number SIGN -- and, because the seed direction biases which of
    # two nearby real tori the nonlinear corrector lands on, occasionally the
    # branch itself -- flips between platforms. Pin arg(lam) into (0, pi) so
    # rho >= 0 everywhere and the corrector is reproducible cross-platform.
    # (#632; cf. #515's 3D-Floquet eigenvector-sign alignment for the same
    # class of platform-dependent eigen-sign ambiguity.)
    if np.imag(lam) < 0.0:
        lam = complex(np.conj(lam))
        v_complex = np.conj(v_complex)
    # Normalize so |v| = 1
    v_complex /= np.linalg.norm(v_complex)

    n_total = 2 * n_modes + 1
    coeffs = np.zeros((n_total, 6), dtype=np.complex128)
    # c_0 = parent state (real)
    coeffs[0, :] = parent_state.astype(np.complex128)
    # n=+1 mode: amplitude * v / 2 so 2 Re(c_1 e^(i theta)) =
    # amplitude * (Re(v) cos theta - Im(v) sin theta)
    coeffs[1, :] = 0.5 * amplitude * v_complex
    # n=-1 mode: conj for reality
    coeffs[n_total - 1, :] = 0.5 * amplitude * np.conj(v_complex)
    return coeffs, lam, v_complex


def _gmos_residual(
    system: cr3bp.CR3BPSystem,
    coeffs: NDArray[np.complex128],
    rho: float,
    t_strob: float,
    *,
    n_samples: int,
) -> NDArray[np.complex128]:
    """Evaluate the GMOS residual in Fourier-mode space::

        F_n = FFT[phi_{t_strob}(u_j)]_n - c_n * exp(i n rho)

    Parameters
    ----------
    coeffs :
        Current Fourier-mode estimate, shape ``(2 N + 1, 6)``.
    rho :
        Current rotation-number estimate (radians).
    t_strob :
        Current stroboscopic-time estimate (nondim TU).
    n_samples :
        ``M >= 2 N + 1``; the FFT grid size.

    Returns
    -------
    NDArray :
        Complex residual of shape ``(n_samples, 6)`` (one row per Fourier mode
        of the propagated circle in numpy-FFT ordering).
    """
    n_total = coeffs.shape[0]
    n_modes = (n_total - 1) // 2
    if n_samples < n_total:
        raise ValueError(f"n_samples={n_samples} must be >= 2 N + 1 = {n_total}")
    thetas = 2 * math.pi * np.arange(n_samples) / n_samples
    u_samples = evaluate_invariant_circle(coeffs, thetas)  # (n_samples, 6)
    # Propagate each sample by t_strob
    phi_samples = np.zeros_like(u_samples)
    for j in range(n_samples):
        arc = cr3bp.propagate(system, u_samples[j], t_strob, with_stm=False)
        phi_samples[j] = arc.state_f
    # FFT each coordinate
    phi_fft = np.fft.fft(phi_samples, axis=0) / n_samples
    # Compare to c_n * exp(i n rho) in mode space, n=0..n_modes and n=-n_modes..-1
    expected = np.zeros((n_samples, 6), dtype=np.complex128)
    # n=0
    expected[0, :] = coeffs[0, :]
    # n=+1..+n_modes
    for n in range(1, n_modes + 1):
        expected[n, :] = coeffs[n, :] * np.exp(1j * n * rho)
    # n=-1..-n_modes -- numpy FFT bins for negative n: bins n_samples - n
    for n in range(1, n_modes + 1):
        expected[n_samples - n, :] = coeffs[n_total - n, :] * np.exp(-1j * n * rho)
    # We only check the modes the truncation tracks; higher modes of the
    # propagated signal are "tail" -- they're indicators of truncation error
    # rather than free-variable errors. Mask them out.
    residual = phi_fft - expected
    # Zero out the high-frequency tail (modes |n| > n_modes); they reflect
    # truncation, not solver state.
    if n_samples > n_total:
        for n in range(n_modes + 1, n_samples - n_modes):
            residual[n, :] = 0.0
    return residual


def _pack_unknowns(
    coeffs: NDArray[np.complex128], rho: float, t_strob: float
) -> NDArray[np.float64]:
    """Pack (coeffs, rho, t_strob) into a real-valued flat vector for the
    least-squares solver. Uses the reality constraint: only c_0 (real, 6) and
    c_1..c_N (complex, 6 each) are independent.
    """
    n_total = coeffs.shape[0]
    n_modes = (n_total - 1) // 2
    # c_0 real (6), c_n complex (12 each for n=1..N), rho (1), t_strob (1)
    n_unk = 6 + 12 * n_modes + 2
    out = np.zeros(n_unk)
    out[0:6] = np.real(coeffs[0, :])
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        out[i0 : i0 + 6] = np.real(coeffs[n, :])
        out[i0 + 6 : i0 + 12] = np.imag(coeffs[n, :])
    out[-2] = rho
    out[-1] = t_strob
    return out


def _unpack_unknowns(
    x: NDArray[np.float64], n_modes: int
) -> tuple[NDArray[np.complex128], float, float]:
    """Inverse of ``_pack_unknowns``. Constructs the full ``(2 N + 1, 6)``
    coefficient array including the negative modes by conjugation.
    """
    n_total = 2 * n_modes + 1
    coeffs = np.zeros((n_total, 6), dtype=np.complex128)
    coeffs[0, :] = x[0:6].astype(np.complex128)
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        c_n = x[i0 : i0 + 6] + 1j * x[i0 + 6 : i0 + 12]
        coeffs[n, :] = c_n
        coeffs[n_total - n, :] = np.conj(c_n)
    rho = float(x[-2])
    t_strob = float(x[-1])
    return coeffs, rho, t_strob


def _residual_real(
    x: NDArray[np.float64],
    system: cr3bp.CR3BPSystem,
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
    amplitude_pin: float,
) -> NDArray[np.float64]:
    """Real-valued residual for the least-squares solver.

    Includes:
      * GMOS residual on the tracked modes (n = -N..N) -- real and imag parts.
      * Phase-locking pin: ``Im(c_1[phase_pin_idx]) = 0`` (kills the rotation
        invariance of the invariant-circle parameterization).
      * Amplitude pin: ``|c_1| = amplitude_pin`` (kills the trivial collapse
        ``c_n -> 0``; enforces a specific torus in the family).
    """
    coeffs, rho, t_strob = _unpack_unknowns(x, n_modes)
    if t_strob <= 0 or not math.isfinite(t_strob) or not math.isfinite(rho):
        return np.full(_residual_size(n_modes), 1e10)
    try:
        f_res = _gmos_residual(system, coeffs, rho, t_strob, n_samples=n_samples)
    except (RuntimeError, ValueError):
        return np.full(_residual_size(n_modes), 1e10)
    # Pack n = 0, 1..N, -N..-1 into a flat real vector.
    parts: list[NDArray[np.float64]] = []
    parts.append(np.real(f_res[0, :]))  # 6 -- n=0 is real
    n_total_sig = f_res.shape[0]
    for n in range(1, n_modes + 1):
        parts.append(np.real(f_res[n, :]))
        parts.append(np.imag(f_res[n, :]))
        parts.append(np.real(f_res[n_total_sig - n, :]))
        parts.append(np.imag(f_res[n_total_sig - n, :]))
    # Phase pin: Im(c_1[phase_pin_idx]) = 0
    parts.append(np.array([float(np.imag(coeffs[1, phase_pin_idx]))]))
    # Amplitude pin: |c_1|^2 - amplitude_pin^2 = 0
    amp = float(np.linalg.norm(coeffs[1, :]))
    parts.append(np.array([amp - amplitude_pin]))
    return np.concatenate(parts)


def _residual_size(n_modes: int) -> int:
    """Length of the real residual vector returned by ``_residual_real``."""
    return 6 + 4 * 6 * n_modes + 2


def _correct_gmos(
    system: cr3bp.CR3BPSystem,
    x0: NDArray[np.float64],
    n_trans: int,
    n_samples: int,
    phase_pin_idx: int,
    amplitude_pin: float,
    tol: float,
    max_iter: int,
) -> tuple[NDArray[np.float64], float, int]:
    """Execute the GMOS least-squares Newton step from an initial guess vector.

    Returns (x_final, residual_norm, n_iter).
    Raises RuntimeError or ValueError if the propagation fails or solver blows up.
    """
    # Least-squares Newton via scipy. Trust-region 'trf' with bounded box
    # keeps propagations physical (the unknowns can't wander off into
    # collisional regions). Empirically reaches an invariance-residual floor
    # of O(1e-7) for n_modes=2 at amp=5e-4 -- the floor is dominated by the
    # finite-difference Jacobian conditioning at this mode count.
    from scipy.optimize import least_squares

    t_strob0 = float(x0[-1])

    # Bounds keep propagations physical: |c_n| stays inside O(10), rho in
    # (-pi, pi), t_strob positive and within 10x the seed value.
    lb = -10.0 * np.ones_like(x0)
    ub = 10.0 * np.ones_like(x0)
    lb[-2] = -math.pi
    ub[-2] = math.pi
    lb[-1] = 0.1 * t_strob0
    ub[-1] = 10.0 * t_strob0

    res = least_squares(
        _residual_real,
        x0,
        args=(system, n_trans, n_samples, phase_pin_idx, amplitude_pin),
        method="trf",
        bounds=(lb, ub),
        xtol=tol * 1e-2,
        ftol=tol * 1e-2,
        gtol=1e-15,
        max_nfev=max_iter * (len(x0) + 1),
        diff_step=np.array(1e-5),
    )
    x_final = res.x
    n_iter = int(res.nfev)
    residual_norm = float(np.linalg.norm(res.fun))
    return x_final, residual_norm, n_iter


def correct_qp_torus(
    system: cr3bp.CR3BPSystem,
    seed_orbit: NDArray[np.float64],
    seed_period: float,
    bifurcation_floquet_pair: tuple[complex, complex],
    *,
    k: int = 4,
    n_long: int = 16,
    n_trans: int = 8,
    initial_torus_amplitude: float = 1e-3,
    tol: float = 1e-8,
    max_iter: int = 50,
    independent_tol: float = 1e-4,
    independent_n_samples: int = 16,
    notes: str = "",
) -> QPTorus:
    """Newton-correct a QP-torus from a Neimark-Sacker bifurcation seed.

    Note: The core Newton-solve logic has been extracted into the private helper
    ``_correct_gmos(...)`` so that it can be reused by continuation drivers.

    Parameters
    ----------
    system :
        CR3BP system at the parent family.
    seed_orbit :
        6-vector IC of the parent periodic orbit at the Neimark-Sacker bracket.
    seed_period :
        Period of the parent orbit (nondim TU).
    bifurcation_floquet_pair :
        The conjugate pair of Floquet multipliers crossing the unit circle.
        Used as a sanity check; the actual eigenvector is recomputed from the
        monodromy because the published multiplier pair may have been
        post-processed.
    k :
        The integer such that the multiplier sits near a primitive k-th root
        of unity. Drives the Neimark-Sacker eigenvector search.
    n_long :
        ``n_long`` is the number of LONGITUDINAL modes; in the GMOS / invariant-
        circle formulation, the longitudinal direction is RESOLVED EXACTLY by
        the stroboscopic flow integration over ``t_strob``, so this parameter
        is currently UNUSED in the corrector but RECORDED on the returned
        torus for traceability (kept as a placeholder for a future
        torus-of-tori extension).
    n_trans :
        ``n_trans`` is the number of TRANSVERSE Fourier modes (``N`` in the
        docstring). The invariant circle is represented by ``2 N + 1`` modes.
    initial_torus_amplitude :
        Seed displacement from the parent orbit along the Neimark-Sacker
        eigenvector. Also serves as the amplitude pin during Newton (kills
        the trivial ``c_n -> 0`` solution).
    tol :
        Convergence tolerance for the GMOS L2 residual.
    max_iter :
        Maximum Newton-step (least-squares-restart) iterations.
    independent_tol :
        Tolerance for the off-grid cross-check.
    independent_n_samples :
        Number of OFF-GRID angles to use for the cross-check.

    Returns
    -------
    QPTorus :
        Converged or best-effort iterate. ``converged`` reflects both the GMOS
        residual gate and the independent closure check.
    """
    if n_trans < 1:
        raise ValueError("n_trans must be >= 1 (need at least one mode)")
    # Validate the supplied Floquet pair (sanity check only)
    lam_a, lam_b = bifurcation_floquet_pair
    if not (abs(abs(lam_a) - 1.0) < 0.2 and abs(abs(lam_b) - 1.0) < 0.2):
        raise ValueError(
            f"bifurcation_floquet_pair must lie near the unit circle; "
            f"got |lam_a|={abs(lam_a):.4f}, |lam_b|={abs(lam_b):.4f}"
        )
    # Reject real pairs (period-doubling / saddle-node) explicitly
    if abs(np.imag(lam_a)) < 1e-6 and abs(np.imag(lam_b)) < 1e-6:
        raise ValueError(
            "bifurcation_floquet_pair is real -- this is period-doubling or "
            "saddle-node, NOT Neimark-Sacker. QP-tori require a complex "
            "(off-real-axis) Floquet pair on the unit circle."
        )

    # Compute monodromy of the parent orbit (we need the eigenvector)
    arc = cr3bp.propagate(system, seed_orbit, seed_period, with_stm=True)
    monod = arc.stm
    assert monod is not None

    # Seed the invariant circle
    coeffs0, lam_seed, _v_complex = _seed_invariant_circle(
        seed_orbit, monod, k=k, n_modes=n_trans, amplitude=initial_torus_amplitude
    )

    # Initial rho from the eigenvalue: lam = e^(i phi) -> rho = phi
    rho0 = math.atan2(np.imag(lam_seed), np.real(lam_seed))
    # Initial stroboscopic time = parent period (the natural choice)
    t_strob0 = float(seed_period)

    # Phase pin: choose the coordinate where Re(v) is largest (avoid pinning
    # a near-zero component)
    coeffs1_init = coeffs0[1, :]
    phase_pin_idx = int(np.argmax(np.abs(np.real(coeffs1_init))))

    # Use 2 N + 3 samples (slight oversampling): the propagated invariant
    # circle contains Fourier content at all orders, and an under-sampled FFT
    # aliases those high modes into the tracked low modes. Mild oversampling
    # caps the aliasing error from mode N + 1 into the tracked range while
    # keeping the propagation cost linear in N.
    n_samples = 2 * n_trans + 3

    x0 = _pack_unknowns(coeffs0, rho0, t_strob0)

    try:
        x_final, residual_norm, n_iter = _correct_gmos(
            system=system,
            x0=x0,
            n_trans=n_trans,
            n_samples=n_samples,
            phase_pin_idx=phase_pin_idx,
            amplitude_pin=initial_torus_amplitude,
            tol=tol,
            max_iter=max_iter,
        )
    except (RuntimeError, ValueError) as e:
        coeffs0 = _enforce_reality(coeffs0)
        return QPTorus(
            system=system,
            omega_long=2 * math.pi / t_strob0,
            omega_trans=rho0 / t_strob0,
            rho=rho0,
            t_strob=t_strob0,
            fourier_coeffs=coeffs0,
            n_modes=n_trans,
            n_samples=n_samples,
            invariance_residual=float("inf"),
            independent_closure_residual=float("inf"),
            converged=False,
            n_iter=0,
            notes=f"newton_failed: {e}; {notes}",
        )

    coeffs_final, rho_final, t_strob_final = _unpack_unknowns(x_final, n_trans)
    coeffs_final = _enforce_reality(coeffs_final)

    omega_long_final = 2 * math.pi / t_strob_final
    omega_trans_final = rho_final / t_strob_final

    # Independent cross-check: pick off-grid sample angles, propagate, check
    # against the resampled invariant circle.
    rng = np.random.default_rng(seed=0xC0FFEE)
    theta_off_grid = rng.uniform(0.0, 2 * math.pi, size=independent_n_samples)
    # Avoid points that happen to coincide with the GMOS grid (very unlikely
    # for floats, but be paranoid)
    grid_thetas = 2 * math.pi * np.arange(n_samples) / n_samples
    for i in range(independent_n_samples):
        while np.any(np.abs(theta_off_grid[i] - grid_thetas) < 1e-6):
            theta_off_grid[i] = rng.uniform(0.0, 2 * math.pi)

    max_indep_err = 0.0
    try:
        for theta_c in theta_off_grid:
            u0 = evaluate_invariant_circle(coeffs_final, theta_c)
            arc_chk = cr3bp.propagate(system, u0, t_strob_final, with_stm=False)
            u_target = evaluate_invariant_circle(coeffs_final, theta_c + rho_final)
            err = float(np.linalg.norm(arc_chk.state_f - u_target))
            max_indep_err = max(max_indep_err, err)
    except RuntimeError as e:
        max_indep_err = float("inf")
        notes = f"{notes}; independent_check_propagation_failed: {e}"

    converged = (residual_norm < tol) and (max_indep_err < independent_tol)

    return QPTorus(
        system=system,
        omega_long=omega_long_final,
        omega_trans=omega_trans_final,
        rho=rho_final,
        t_strob=t_strob_final,
        fourier_coeffs=coeffs_final,
        n_modes=n_trans,
        n_samples=n_samples,
        invariance_residual=residual_norm,
        independent_closure_residual=max_indep_err,
        converged=converged,
        n_iter=n_iter,
        notes=notes,
        extras={
            "n_long_recorded": float(n_long),
            "k": float(k),
            "lam_seed_re": float(np.real(lam_seed)),
            "lam_seed_im": float(np.imag(lam_seed)),
            "phase_pin_idx": float(phase_pin_idx),
            "amplitude_pin": float(initial_torus_amplitude),
        },
    )


# ---------------------------------------------------------------------------
# Diagnostic: is the frequency ratio practically irrational?
# ---------------------------------------------------------------------------


def is_practically_irrational(
    ratio: float, *, max_denominator: int = 10, tol: float = 1e-3
) -> bool:
    """Return ``True`` if ``ratio`` is NOT within ``tol`` of any small rational
    ``p/q`` with ``q <= max_denominator``.

    Distinguishes a genuine QP-torus (irrational frequency ratio) from a
    phase-locked periodic orbit (rational frequency ratio = small p:q
    resonance).
    """
    if not math.isfinite(ratio):
        return False
    r = ratio - math.floor(ratio)  # reduce to [0, 1)
    for q in range(1, max_denominator + 1):
        for p in range(0, q + 1):
            if abs(r - p / q) < tol:
                return False
    return True
