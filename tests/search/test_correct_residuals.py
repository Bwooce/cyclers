"""M-ED Phase 1: ballistic-closure residual vector (plan Phase 1 Task 1.2)."""

from __future__ import annotations

from cyclerfinder.search.correct import _residual_vector


def test_residual_length_is_n_minus_one_intermediates_plus_closure() -> None:
    # E-M-E-E: encounters B0..B3, intermediates B1,B2 -> 2 residuals + 1 closure.
    fake = {
        "b1_in": (3.0, 0.0, 0.0),
        "b1_out": (3.0, 0.0, 0.0),
        "b2_in": (5.0, 0.0, 0.0),
        "b2_out": (5.0, 0.0, 0.0),
        "b3_in": (5.6, 0.0, 0.0),
        "b0_out": (5.6, 0.0, 0.0),
    }
    res = _residual_vector(fake, n_encounters=4)
    assert len(res) == 3
    assert max(abs(r) for r in res) < 1e-12  # perfectly continuous fixture
