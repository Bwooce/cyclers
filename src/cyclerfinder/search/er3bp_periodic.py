"""ER3BP periodic-orbit utilities for the Fitzgerald 2022 transit-gate pipeline.

Thin helpers that sit between the corrector (genome/er3bp_periodic.py) and the
Floquet classifier (search/er3bp_floquet.py), providing:
  - Canonical-momentum → pulsating-frame state conversion
  - Monodromy eigenstructure extraction (saddle r, center angle w)

Barcelona-school canonical coordinates (used by Fitzgerald 2022, Jorba, Rosales):
  px = x' - y  →  x' = px + y
  py = y' + x  →  y' = py - x
(Legendre transform of the effective Hamiltonian under the rotating frame.)
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def canonical_to_er3bp_state(
    x: float,
    y: float,
    px: float,
    py: float,
) -> NDArray[np.float64]:
    """Convert Barcelona-school canonical IC [x, y, px, py] to ER3BP pulsating state.

    The Barcelona-school (Jorba/Rosales) canonical-momentum variables satisfy:
        px = x' - y  =>  x' = px + y
        py = y' + x  =>  y' = py - x

    where x' = dx/df etc. are velocities in the pulsating rotating frame
    (true anomaly f as independent variable).  Planar (z = z' = 0) is assumed.

    Args:
        x: x-coordinate (nondimensional).
        y: y-coordinate (nondimensional).
        px: Canonical x-momentum = x' - y.
        py: Canonical y-momentum = y' + x.

    Returns:
        6-element state [x, y, z, x', y', z'] at f=0.

    References:
        Fitzgerald J., Ross S.D. (2022), Adv. Space Res. 70:144-156,
        DOI 10.1016/j.asr.2022.04.029.
        Convention: Jorba A., Rosales J.J. (2021), Celestial Mech. Dyn. Astron.
        132, 11.
    """
    xdot = px + y
    ydot = py - x
    return np.array([x, y, 0.0, xdot, ydot, 0.0], dtype=np.float64)


def monodromy_eigenstructure(
    monodromy: NDArray[np.float64],
) -> tuple[float, float]:
    """Extract the saddle eigenvalue r and center angle w from a 6x6 monodromy.

    For a Lagrange periodic orbit (L1/L2 replacement in ER3BP/BCP), the monodromy
    eigenvalues come in two symplectic pairs relevant to planar motion:
        - Saddle pair: (r, 1/r)  with r >> 1  (hyperbolic character)
        - Center pair: (e^{iw}, e^{-iw})  with |λ|=1, Im ≠ 0  (elliptic character)
    plus an out-of-plane pair.

    Args:
        monodromy: 6x6 monodromy matrix from propagate_er3bp with with_stm=True.

    Returns:
        (r, w) where:
            r = max(|λ|) over all eigenvalues (the saddle eigenvalue).
            w = arg of the largest-magnitude unit-circle eigenvalue with Im ≠ 0
                (the center rotation angle, in radians, 0 < w < pi).

    Raises:
        ValueError: If no center (near-unit-circle complex) eigenvalue is found.

    References:
        Fitzgerald J., Ross S.D. (2022), Adv. Space Res. 70:144-156,
        DOI 10.1016/j.asr.2022.04.029. Eq. B.4: r = e^{k̃T}, w = m̃T.
    """
    eigs = np.linalg.eigvals(np.asarray(monodromy, dtype=np.float64))

    # Saddle eigenvalue: largest magnitude.
    r = float(np.max(np.abs(eigs)))

    # Center eigenvalue: unit-circle (|λ| within 0.5 of 1.0) with non-zero Im.
    center_mask = (np.abs(np.abs(eigs) - 1.0) < 0.5) & (np.abs(eigs.imag) > 1e-10)
    center_eigs = eigs[center_mask]

    if center_eigs.size == 0:
        raise ValueError(
            "No center (near-unit-circle complex) eigenvalue found in monodromy. "
            f"All |λ|: {np.abs(eigs).tolist()}"
        )

    # Take the one with the largest magnitude among the candidates (most "central").
    best_idx = int(np.argmax(np.abs(center_eigs)))
    w = float(np.abs(np.angle(center_eigs[best_idx])))

    return r, w
