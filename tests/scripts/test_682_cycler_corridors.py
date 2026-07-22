"""Task #682: quasi-periodic cycler-corridor census tests.

Covers the reusable pieces of ``scripts/census_682_cycler_corridors.py``:

1. :func:`_conj_pairs` -- correct grouping of monodromy center eigenvalues into
   conjugate pairs (and rejection of the trivial/real multipliers), the step that
   selects which KAM family to measure.
2. :func:`_best_k` and ``AMP_LADDER`` structural invariants.
3. :func:`_stable_planar_goldens` / :func:`_all_members` -- the deterministic
   stratified sample this census reports on (pins the member set + the live
   stability classification of the planar Braik-Ross goldens: R21-S/R31-S/R52-S
   are stable, the Lyapunov/cycler/DPO goldens are not).
4. A SLOW end-to-end evidence test: :func:`measure_corridor` on one cheap
   lyapunov3d-L1 stable member returns a genuine (closure-gated) corridor with a
   positive physical tube width -- the census's own positive-control-style check
   that the GMOS ladder measures a real station-keeping-free corridor. (The torus
   corrector itself is separately validated by `#612`'s in-repo L2 GMOS positive
   control, tests/search/test_variational_qp_torus.py.)
"""

from __future__ import annotations

import numpy as np
import pytest

import scripts.census_682_cycler_corridors as census


def test_conj_pairs_groups_and_rejects() -> None:
    # A stable orbit's monodromy: trivial pair {1,1}, plus TWO center pairs.
    e1 = np.exp(1j * 0.5)
    e2 = np.exp(1j * 2.1)
    eigs = [1.0 + 0j, 1.0 + 0j, e1, e1.conjugate(), e2, e2.conjugate()]
    pairs = census._conj_pairs(eigs)
    assert len(pairs) == 2
    for a, b in pairs:
        assert a.imag > 0  # imag>0 element listed first
        assert abs(a - b.conjugate()) < 1e-9
    # A real saddle pair (lambda, 1/lambda off the unit circle) yields no pair.
    assert census._conj_pairs([1.0 + 0j, 1.0 + 0j, 2.5 + 0j, 0.4 + 0j]) == []


def test_best_k_and_ladder() -> None:
    # rho = 1/4 -> k = 4 is the nearest low-order primitive root.
    assert census._best_k(2 * np.pi * 0.25) == 4
    assert list(census.AMP_LADDER) == sorted(census.AMP_LADDER)
    assert census.AMP_LADDER[0] > 0.0


def test_stable_planar_goldens_are_the_resonant_stable_families() -> None:
    goldens = census._stable_planar_goldens()
    labels = {str(g["label"]) for g in goldens}
    # Only the resonant "-S" (stable) goldens are linearly stable at CJ~3.1294;
    # the Lyapunov, cycler and DPO goldens are not.
    assert labels == {"R21-S", "R31-S", "R52-S"}
    for g in goldens:
        assert float(g["spectral_radius"]) < 1.05


def test_all_members_is_the_deterministic_stratified_sample() -> None:
    members = census._all_members()
    from collections import Counter

    comp = Counter(str(m["family"]) for m in members)
    assert comp == {
        "lyapunov3d-L1": 7,
        "braik-ross-C21-em-z0_0.24": 5,
        "braik-ross-C32-em-z0_0.24": 5,
        "planar_golden": 3,
    }
    assert len(members) == 20
    # member_id is unique (checkpoint key).
    assert len({m["member_id"] for m in members}) == len(members)


@pytest.mark.slow
def test_measure_corridor_on_cheap_lyapunov3d_member() -> None:
    """End-to-end: a cheap stable lyapunov3d-L1 member has a genuine, closure-gated
    quasi-periodic torus corridor of positive physical width."""
    lyap = [m for m in census._all_members() if m["family"] == "lyapunov3d-L1"]
    m = lyap[len(lyap) // 2]
    s0 = np.array([m["x0"], 0.0, m["z0"], 0.0, m["ydot0"], 0.0])
    out = census.measure_corridor(census._SYS, s0, float(m["T_TU"]))
    assert out["status"] == "ok"
    assert float(out["corridor_amp_nondim"]) >= census.AMP_LADDER[0]
    assert float(out["corridor_pos_km"]) > 0.0
    assert float(out["best_closure"]) < census.CLOSURE_GATE
