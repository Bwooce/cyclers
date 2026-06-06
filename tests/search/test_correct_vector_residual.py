"""Phase 1 (task #122): full v-inf VECTOR continuity residual mode.

The default residual is |V_inf|-magnitude continuity (unchanged). The new
``residual_mode="vector"`` adds, per intermediate flyby, a bend-feasibility term
INSIDE the residual: the part of the required in->out rotation that EXCEEDS the
body's V_inf-limited max bend (Jones-method feasibility steering the solve, not a
post-hoc filter). Magnitude continuity is still enforced.

NON-GOLDEN: every V_inf here is OUR computation -- solver fixtures, not anchors.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.search.correct import _residual_vector


def test_vector_mode_zero_when_continuous_and_no_bend() -> None:
    # E-M-E-E: intermediates B1,B2. Perfectly continuous AND zero required bend
    # (in == out direction) -> all residual terms zero in vector mode too.
    fake = {
        "b1_in": (3.0, 0.0, 0.0),
        "b1_out": (3.0, 0.0, 0.0),
        "b2_in": (5.0, 0.0, 0.0),
        "b2_out": (5.0, 0.0, 0.0),
        "b3_in": (5.6, 0.0, 0.0),
        "b0_out": (5.6, 0.0, 0.0),
    }
    res = _residual_vector(fake, n_encounters=4, mode="vector", sequence=("E", "M", "E", "E"))
    assert max(abs(r) for r in res) < 1e-12


def test_vector_mode_penalises_infeasible_bend_at_low_vinf() -> None:
    # A 90-degree required turn at low V_inf (3 km/s) at Mars exceeds the max
    # single-flyby bend -> vector mode must produce a NON-ZERO feasibility
    # residual even though magnitude continuity holds (|in| == |out|).
    fake = {
        "b1_in": (3.0, 0.0, 0.0),
        "b1_out": (0.0, 3.0, 0.0),  # 90 deg turn, same magnitude
        "b2_in": (5.0, 0.0, 0.0),
        "b2_out": (5.0, 0.0, 0.0),
        "b3_in": (5.6, 0.0, 0.0),
        "b0_out": (5.6, 0.0, 0.0),
    }
    res = _residual_vector(fake, n_encounters=4, mode="vector", sequence=("E", "M", "E", "E"))
    assert max(abs(r) for r in res) > 0.1  # infeasible bend surfaces in residual


def test_magnitude_mode_ignores_infeasible_bend() -> None:
    # The SAME 90-degree infeasible turn is invisible to the default magnitude
    # mode (|in| == |out| so the magnitude term is zero) -- this is exactly the
    # design's central finding (the magnitude residual cannot see bend).
    fake = {
        "b1_in": (3.0, 0.0, 0.0),
        "b1_out": (0.0, 3.0, 0.0),
        "b2_in": (5.0, 0.0, 0.0),
        "b2_out": (5.0, 0.0, 0.0),
        "b3_in": (5.6, 0.0, 0.0),
        "b0_out": (5.6, 0.0, 0.0),
    }
    res_mag = _residual_vector(fake, n_encounters=4)  # default = magnitude
    assert max(abs(r) for r in res_mag) < 1e-12


def test_vector_mode_penalises_magnitude_discontinuity() -> None:
    # A magnitude jump (|in| != |out|) is penalised in vector mode just as in
    # magnitude mode -- vector mode is a superset, not a replacement.
    fake = {
        "b1_in": (3.0, 0.0, 0.0),
        "b1_out": (4.0, 0.0, 0.0),  # +1 km/s magnitude jump, zero bend
        "b2_in": (5.0, 0.0, 0.0),
        "b2_out": (5.0, 0.0, 0.0),
        "b3_in": (5.6, 0.0, 0.0),
        "b0_out": (5.6, 0.0, 0.0),
    }
    res = _residual_vector(fake, n_encounters=4, mode="vector", sequence=("E", "M", "E", "E"))
    assert max(abs(r) for r in res) > 0.9  # the 1 km/s jump shows up


def test_vector_mode_requires_sequence() -> None:
    fake = {
        "b1_in": (3.0, 0.0, 0.0),
        "b1_out": (3.0, 0.0, 0.0),
        "b2_in": (5.0, 0.0, 0.0),
        "b2_out": (5.0, 0.0, 0.0),
        "b3_in": (5.6, 0.0, 0.0),
        "b0_out": (5.6, 0.0, 0.0),
    }
    try:
        _residual_vector(fake, n_encounters=4, mode="vector")
    except ValueError:
        return
    raise AssertionError("vector mode without sequence must raise ValueError")


def test_magnitude_mode_default_unchanged_signature() -> None:
    # Byte-identical default behaviour: calling with no mode kwarg returns the
    # same magnitude residual the existing callers expect (3 terms for E-M-E-E).
    fake = {
        "b1_in": (3.0, 0.0, 0.0),
        "b1_out": (3.2, 0.0, 0.0),
        "b2_in": (5.0, 0.0, 0.0),
        "b2_out": (4.7, 0.0, 0.0),
        "b3_in": (5.6, 0.0, 0.0),
        "b0_out": (5.6, 0.0, 0.0),
    }
    res = _residual_vector(fake, n_encounters=4)
    assert len(res) == 3
    assert res[0] == np.float64(3.0) - np.float64(3.2)
