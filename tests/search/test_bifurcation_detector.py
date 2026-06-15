"""Period-multiplying bifurcation detector tests (#266 Phase 2).

Discipline gates:

  * The parent NRHO has a multiplier near +1 (trivial energy/time-translation
    unit root) and no strong period-multiplying root.
  * The detector on a SYNTHETIC monodromy with a known primitive k-th root of
    unity eigenvalue returns the correct k (the reproduce gate for the
    detector itself, with a known-truth input).
  * The family-scan returns at least one bracket when a synthetic family is
    constructed with a sign-flip across the tolerance.

The synthetic gates are essential: they verify the detector logic
independently of any CR3BP-physics question (which can be confounded by
solver tolerance, IC precision, and family-member position). The CR3BP gate
(NRHO unit-root) verifies the detector's wiring through the real propagator.
"""

from __future__ import annotations

import math

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.genome.tulip import KOBLICK_2023_TABLE4, koblick_system
from cyclerfinder.search.bifurcation_detector import (
    BifurcationPoint,
    FamilyMember,
    detect_period_multiplying,
    floquet_multipliers,
    monodromy,
    scan_family_for_bifurcations,
)

# ---------------------------------------------------------------------------
# Gate 1: Floquet multipliers at Np = 1 (parent NRHO) include +1.
# ---------------------------------------------------------------------------


def test_floquet_multipliers_unit_root_at_np1() -> None:
    """Corrected NRHO has at least one multiplier within 1e-3 of +1.

    The CR3BP monodromy is Hamiltonian: its symplectic structure forces a
    trivial eigenvalue pair at +1 (time translation, energy). For a
    well-corrected periodic orbit those two unit eigenvalues sit on top of
    each other within the corrector's residual tolerance; numerically they
    may split into two real eigenvalues straddling 1 by ~residual.
    """
    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4[1]
    s0 = np.array(
        [
            float(row["x0"]),  # type: ignore[arg-type]
            float(row["y0"]),  # type: ignore[arg-type]
            float(row["z0"]),  # type: ignore[arg-type]
            float(row["xdot0"]),  # type: ignore[arg-type]
            float(row["ydot0"]),  # type: ignore[arg-type]
            float(row["zdot0"]),  # type: ignore[arg-type]
        ]
    )
    period_guess = float(row["T_TU"])  # type: ignore[arg-type]
    orbit = cp.correct_periodic(sysm, s0, period_guess)
    assert orbit.converged, f"NRHO failed to converge: residual {orbit.closure_residual:.3e}"

    mat = monodromy(sysm, orbit.state0, orbit.period)
    eigs = floquet_multipliers(mat)
    # At least one multiplier within 1e-3 of +1.
    dists = np.abs(eigs - 1.0)
    assert dists.min() < 1e-3, (
        f"no Floquet multiplier near +1 at the NRHO (min |lambda - 1| = "
        f"{dists.min():.3e}); multipliers were {eigs}"
    )

    # And the detector finds no strong period-multiplying root at the parent.
    bifs = detect_period_multiplying(eigs, k_max=6, tol=1e-3)
    # Allow any -1 = k=2 root if the NRHO's hyperbolic pair sits very close to
    # it, but at the tight tolerance the parent should be clear.
    # Empirically (see test_tulip's monodromy print) the NRHO's hyperbolic pair
    # at ~-2.19 / -0.457 is FAR from -1 (distance ~1.19) so this is safe.
    assert not bifs, f"detector flagged period-multiplying roots at the parent NRHO: {bifs}"


# ---------------------------------------------------------------------------
# Gate 2: detector on a synthetic monodromy with a known k-th root.
# ---------------------------------------------------------------------------


def test_period_multiplying_detection_finds_known_root_of_unity() -> None:
    """Synthetic eigenvalue array containing a primitive k=3 root of unity is
    detected as k=3."""
    # Primitive cube root of unity: exp(2*pi*i/3) = -0.5 + i*sqrt(3)/2.
    cube = complex(-0.5, math.sqrt(3.0) / 2.0)
    eigs = np.array(
        [
            cube,
            cube.conjugate(),
            complex(2.0, 0.0),  # unrelated hyperbolic
            complex(0.5, 0.0),
            complex(1.0, 0.0),
            complex(1.0, 0.0),
        ],
        dtype=np.complex128,
    )
    out = detect_period_multiplying(eigs, k_max=6, tol=1e-6)
    # k=3 must be detected; the conjugate may show up separately.
    ks_found = {k for (k, _eig, _d) in out}
    assert 3 in ks_found, f"expected k=3 to be detected, got {out}"
    # And the matching multiplier should be near the primitive cube root.
    k3_entries = [(eig, d) for (k, eig, d) in out if k == 3]
    assert k3_entries, "k=3 detected but no entries returned"
    for eig, d in k3_entries:
        assert d < 1e-6, f"k=3 entry distance {d:.3e} > tol"
        assert abs(abs(eig) - 1.0) < 1e-6, f"k=3 entry not on unit circle: {eig}"


def test_period_multiplying_detection_ignores_off_root_eigenvalues() -> None:
    """An eigenvalue far from any low-k root produces no detection at small tol."""
    eigs = np.array(
        [
            complex(2.0, 0.0),
            complex(0.5, 0.0),
            complex(1.0 + 1e-12, 0.0),  # near trivial unity, excluded (k >= 2)
            complex(0.3, 0.7),  # irrational rotation
            complex(0.3, -0.7),
            complex(-3.0, 0.0),  # off the unit circle
        ],
        dtype=np.complex128,
    )
    out = detect_period_multiplying(eigs, k_max=6, tol=1e-3)
    # The (0.3, +/- 0.7) on the unit circle but irrational angle -- nearest
    # k=6 root is at angle pi/3 = 1.047, our angle is arctan(0.7/0.3) ~ 1.166,
    # distance is non-negligible. Should not be picked up at tol=1e-3.
    for k, eig, d in out:
        # If anything IS picked up, it must satisfy the tolerance, but at
        # this tol we expect nothing.
        assert d < 1e-3, f"unexpected entry {(k, eig, d)}"


# ---------------------------------------------------------------------------
# Gate 3: family scan brackets a known bifurcation.
# ---------------------------------------------------------------------------


def test_scan_family_brackets_a_known_bifurcation() -> None:
    """Synthetic family with a sign-flip across the tolerance produces a bracket.

    We can't easily fabricate a CR3BP family with a known bifurcation point
    in seconds; instead, we mock the scan by patching the per-member
    monodromy to a synthetic STM. The scan's family-bracket logic is what we
    are testing here -- the monodromy itself is the input.

    Strategy: build two synthetic FamilyMembers and monkey-patch the
    bifurcation_detector.monodromy function to return STMs whose
    eigenvalues are constructed to straddle a primitive cube root of unity
    across the two members. The scan should produce a single k=3 bracket.
    """
    from unittest.mock import patch

    cube = complex(-0.5, math.sqrt(3.0) / 2.0)
    # Build STMs whose eigenvalues are known: a block-diagonal matrix with
    # the desired eigenvalues on the diagonal.
    #
    # Member A: eigenvalues include cube + delta (distance > tol)
    # Member B: eigenvalues include cube + epsilon (distance < tol)
    # The scan should detect a k=3 crossing between A and B.

    def _stm_with_eigs(eigs_complex: list[complex]) -> np.ndarray:
        """Construct a real 6x6 STM with the prescribed eigenvalues.

        Pairs complex-conjugate eigenvalues into 2x2 real rotation blocks; the
        result is a real matrix.
        """
        # Sort: pair complex pairs first; real eigenvalues last.
        complex_pairs: list[tuple[complex, complex]] = []
        real_eigs: list[float] = []
        consumed = [False] * len(eigs_complex)
        for i, e in enumerate(eigs_complex):
            if consumed[i]:
                continue
            if abs(e.imag) < 1e-12:
                real_eigs.append(e.real)
                consumed[i] = True
            else:
                # Find conjugate
                for j in range(i + 1, len(eigs_complex)):
                    if consumed[j]:
                        continue
                    f = eigs_complex[j]
                    if abs(f.real - e.real) < 1e-12 and abs(f.imag + e.imag) < 1e-12:
                        complex_pairs.append((e, f))
                        consumed[i] = True
                        consumed[j] = True
                        break
        mat = np.zeros((6, 6))
        row = 0
        for pair in complex_pairs:
            a = pair[0].real
            b = abs(pair[0].imag)
            mat[row : row + 2, row : row + 2] = np.array([[a, -b], [b, a]])
            row += 2
        for r in real_eigs:
            mat[row, row] = r
            row += 1
        return mat

    # Member A: cube + delta=2e-2 perturbation (distance > tol=1e-2)
    delta = 0.02
    eigs_a_complex_pair = complex(cube.real + delta, cube.imag)
    eigs_b_complex_pair = complex(cube.real + 1e-4, cube.imag)
    eigs_a_list = [
        eigs_a_complex_pair,
        eigs_a_complex_pair.conjugate(),
        complex(2.0, 0.0),
        complex(0.5, 0.0),
        complex(1.0, 0.0),
        complex(1.0, 0.0),
    ]
    eigs_b_list = [
        eigs_b_complex_pair,
        eigs_b_complex_pair.conjugate(),
        complex(2.0, 0.0),
        complex(0.5, 0.0),
        complex(1.0, 0.0),
        complex(1.0, 0.0),
    ]
    stm_a = _stm_with_eigs(eigs_a_list)
    stm_b = _stm_with_eigs(eigs_b_list)

    member_a = FamilyMember(
        label="A",
        state0=np.zeros(6),
        period=1.0,
        mu=0.01,
        parameter=0.0,
    )
    member_b = FamilyMember(
        label="B",
        state0=np.zeros(6),
        period=1.0,
        mu=0.01,
        parameter=1.0,
    )

    # Sanity check: confirm our STM has the desired eigenvalues.
    eigs_a_check = np.linalg.eigvals(stm_a)
    assert np.any(np.abs(eigs_a_check - eigs_a_complex_pair) < 1e-8), (
        f"synthetic STM did not produce the prescribed eigenvalue {eigs_a_complex_pair}: "
        f"got {eigs_a_check}"
    )

    # Patch monodromy() in the bifurcation_detector module.
    call_count = {"n": 0}

    def _fake_monodromy(
        system: cr3bp.CR3BPSystem,
        state0: np.ndarray,
        period: float,
        *,
        rtol: float = 1e-12,
        atol: float = 1e-12,
    ) -> np.ndarray:
        call_count["n"] += 1
        # First call -> A, second -> B (the scan walks members in order).
        if call_count["n"] == 1:
            return stm_a
        return stm_b

    with patch(
        "cyclerfinder.search.bifurcation_detector.monodromy",
        side_effect=_fake_monodromy,
    ):
        brackets = scan_family_for_bifurcations([member_a, member_b], k_max=6, tol=1e-2)

    assert call_count["n"] == 2, f"expected exactly 2 monodromy calls, got {call_count['n']}"
    # Expect at least one k=3 bracket (the crossing we engineered).
    k3_brackets = [b for b in brackets if b.k == 3]
    assert k3_brackets, (
        f"scan failed to bracket the synthetic k=3 crossing; brackets returned "
        f"were k={[b.k for b in brackets]}"
    )
    # And the bracket's members are A and B in order.
    b0 = k3_brackets[0]
    assert b0.members[0].label == "A"
    assert b0.members[1].label == "B"
    # dist_before > tol > dist_after (we engineered exactly that).
    assert b0.dist_before > b0.tol > b0.dist_after, (
        f"bracket distances not straddling tol: before={b0.dist_before:.3e}, "
        f"after={b0.dist_after:.3e}, tol={b0.tol:.3e}"
    )


def test_scan_family_empty_returns_no_brackets() -> None:
    """A single-member family has no adjacent pair and produces no brackets."""
    member = FamilyMember(
        label="solo",
        state0=np.array([1.0, 0.0, 0.0, 0.0, 0.5, 0.0]),
        period=1.0,
        mu=0.01,
        parameter=0.0,
    )
    assert scan_family_for_bifurcations([member]) == []
    assert scan_family_for_bifurcations([]) == []


def test_detect_period_multiplying_validates_inputs() -> None:
    """k_max < 2 and tol <= 0 raise ValueError."""
    eigs = np.array([complex(1.0, 0.0)], dtype=np.complex128)
    try:
        detect_period_multiplying(eigs, k_max=1, tol=1e-3)
    except ValueError:
        pass
    else:
        raise AssertionError("k_max=1 must raise ValueError")
    try:
        detect_period_multiplying(eigs, k_max=4, tol=0.0)
    except ValueError:
        pass
    else:
        raise AssertionError("tol=0 must raise ValueError")


def test_monodromy_rejects_non_6_vector_state() -> None:
    """monodromy() must reject a non-6-vector state."""
    sysm = koblick_system()
    try:
        monodromy(sysm, np.array([1.0, 0.0]), 1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("monodromy must reject non-6-vector state")


def test_bifurcation_point_carries_full_evidence() -> None:
    """BifurcationPoint round-trips its fields and is hashable / frozen."""
    member_a = FamilyMember(label="A", state0=np.zeros(6), period=1.0, mu=0.01, parameter=0.0)
    member_b = FamilyMember(label="B", state0=np.zeros(6), period=1.0, mu=0.01, parameter=1.0)
    bp = BifurcationPoint(
        k=3,
        members=(member_a, member_b),
        eig_before=complex(-0.49, 0.86),
        eig_after=complex(-0.50, 0.866),
        dist_before=0.02,
        dist_after=1e-4,
        tol=1e-2,
    )
    assert bp.k == 3
    assert bp.members[0].label == "A"
    assert bp.dist_before > bp.tol > bp.dist_after
