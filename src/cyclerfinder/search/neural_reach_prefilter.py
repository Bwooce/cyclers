"""Tier-0 neural reachability prefilter (#276, Zhang-Acciarini-Izzo-Baoyin-Topputo 2026).

Wraps the pretrained low-thrust ΔV / minimum-transfer-time approximator of
Zhang-Acciarini-Izzo-Baoyin-Topputo 2026 ("Pretrained Approximators for
Low-Thrust Trajectory Cost and Reachability", arXiv:2605.26790) as a tier-0
prefilter that runs BEFORE the existing accessibility tiers:

* Tier 0 (this module): pretrained-NN ΔV estimate, ms-scale per pair, drops
  obviously-infeasible heliocentric low-thrust legs.
* Tier 1: Braik-Ross energy-preserving heading-fan (:mod:`reachable_network`).
* Tier 2: Zhou-Armellin energy-changing single-impulse footprint
  (:mod:`reachable_impulsive`).
* Tier 3: Kumar resonant-manifold network (:mod:`resonance_network`).

The model and weights come from the paper authors' open-source release at
https://github.com/zhong-zh15/neural-low-thrust-approximator. The published
architecture is a 10-layer MLP (input 16 -> 9 hidden of 128 with ReLU ->
output 1, no activation), with self-similar (rotation + non-dimensional)
input preprocessing and a Lambert-derived 16-vector feature, followed by a
``max(NN_output, lambert_floor)`` floor for physical consistency. We
hand-roll the inference in numpy (matching the ``EigenFastPredictor``
implementation in the paper's C++/Eigen tree); the published weight CSVs
are loaded at construction time.

Honest scope of this module
---------------------------
This is the **inference + composition glue**, not a retrain. The fit method
is deliberately not implemented — retraining requires the ~10^7 - 10^8
sample homotopy-ray dataset on Zenodo (DOI 10.5281/zenodo.18769170) which
is well outside this module's responsibility.

Non-blocking discipline (#256 ML-flagger pattern)
-------------------------------------------------
If the pretrained weights cannot be loaded (weights missing, format
unsupported, dim mismatch, NaN, Lambert solver fails on the geometry, etc.)
the prefilter MUST NOT block downstream tiers. ``score_pair`` returns
``model_available=False, prefilter_admitted=True`` (let downstream tiers
decide). The point of a prefilter is to *drop* obviously-infeasible
candidates fast; if the model isn't available, just don't drop anything.
Same discipline as :mod:`cyclerfinder.ml.falsepos_flagger`.

Domain caveat (#276 honesty)
----------------------------
The pretrained model was trained on heliocentric low-thrust transfers
(asteroid rendezvous, multi-flyby missions). Cyclers in our repo are
Earth-Moon CR3BP (rotating-frame, nondimensional) at the discovery-program
level. This module does NOT auto-translate CR3BP reps to heliocentric SI
state vectors; that translation is the caller's responsibility. When the
caller supplies the heliocentric state (via ``heliocentric_state_*``
keyword arguments), the prefilter runs the NN; otherwise it falls back to
non-blocking admit. This is the same "non-blocking by design" contract
that makes the prefilter safe to wire into any pipeline.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.lambert as core_lambert

# ---------------------------------------------------------------------------
# Architecture constants (must match the published EigenFastPredictor).
# ---------------------------------------------------------------------------

#: 16-D input vector to the NN (after self-similar normalization).
# Vendored pretrained-weight directories (see neural_reach_models/NOTICE.md for
# upstream attribution + MPL-2.0 license terms). Used as the default when
# ``from_weight_dir(...)`` is called without explicit paths.
_MODELS_DIR = Path(__file__).resolve().parent / "neural_reach_models"
_VENDORED_DV_DIR = _MODELS_DIR / "eigen_model_large"
_VENDORED_TMIN_DIR = _MODELS_DIR / "eigen_model_tmin_large"

INPUT_DIM = 16
#: 128 neurons per hidden layer.
HIDDEN_DIM = 128
#: Scalar regression output.
OUTPUT_DIM = 1
#: Total Dense layers: input -> 8 hidden -> output (Z=W*x+b each).
NUM_LAYERS = 10
#: Sea-level gravity used by the paper's ΔV normalisation.
G0_M_S2 = 9.80665


# ---------------------------------------------------------------------------
# Lightweight weight container.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Weights:
    """Layer weights, biases, and input/output scaling for one NN model.

    Mirrors the on-disk CSV layout published by the paper authors at
    https://github.com/zhong-zh15/neural-low-thrust-approximator
    (``models/eigen_model_large/*.csv`` for ΔV; ``eigen_model_tmin_large``
    for minimum transfer time).
    """

    weights: tuple[NDArray[np.float64], ...]  # length NUM_LAYERS
    biases: tuple[NDArray[np.float64], ...]
    x_mean: NDArray[np.float64]  # (INPUT_DIM,)
    x_scale: NDArray[np.float64]
    y_mean: float
    y_scale: float


def _load_csv_matrix(path: Path) -> NDArray[np.float64]:
    """Load an Eigen-style CSV matrix (rows = lines, cols = comma-split)."""
    return np.loadtxt(path, delimiter=",", dtype=np.float64, ndmin=2)


def _load_csv_vector(path: Path) -> NDArray[np.float64]:
    """Load an Eigen-style CSV vector (one value per line)."""
    arr = np.loadtxt(path, delimiter=",", dtype=np.float64)
    return np.atleast_1d(arr).ravel()


def load_weights(model_dir: str | Path) -> _Weights:
    """Load the 10-layer NN weights + scalers from a directory of CSVs.

    Expected layout (matching the published ``eigen_model_large`` /
    ``eigen_model_tmin_large`` trees)::

        layer0_weight.csv  layer0_bias.csv
        layer1_weight.csv  layer1_bias.csv
        ...
        layer9_weight.csv  layer9_bias.csv
        X_mean.csv  X_scale.csv  Y_mean.csv  Y_scale.csv

    Raises
    ------
    FileNotFoundError
        If any required file is missing.
    ValueError
        On shape mismatches against (INPUT_DIM, HIDDEN_DIM, NUM_LAYERS).
    """
    d = Path(model_dir)
    weights: list[NDArray[np.float64]] = []
    biases: list[NDArray[np.float64]] = []
    for i in range(NUM_LAYERS):
        w = _load_csv_matrix(d / f"layer{i}_weight.csv")
        b = _load_csv_vector(d / f"layer{i}_bias.csv")
        if i == 0:
            if w.shape != (HIDDEN_DIM, INPUT_DIM):
                raise ValueError(
                    f"layer0 weight expected ({HIDDEN_DIM}, {INPUT_DIM}), got {w.shape}"
                )
        elif i == NUM_LAYERS - 1:
            if w.shape != (OUTPUT_DIM, HIDDEN_DIM):
                raise ValueError(
                    f"layer{i} (output) weight expected ({OUTPUT_DIM}, {HIDDEN_DIM}), got {w.shape}"
                )
        elif w.shape != (HIDDEN_DIM, HIDDEN_DIM):
            raise ValueError(
                f"layer{i} weight expected ({HIDDEN_DIM}, {HIDDEN_DIM}), got {w.shape}"
            )
        weights.append(w)
        biases.append(b)
    x_mean = _load_csv_vector(d / "X_mean.csv")
    x_scale = _load_csv_vector(d / "X_scale.csv")
    if x_mean.size != INPUT_DIM or x_scale.size != INPUT_DIM:
        raise ValueError(
            f"X_mean/X_scale expected ({INPUT_DIM},), got {x_mean.shape}/{x_scale.shape}"
        )
    y_mean = float(_load_csv_vector(d / "Y_mean.csv")[0])
    y_scale = float(_load_csv_vector(d / "Y_scale.csv")[0])
    return _Weights(
        weights=tuple(weights),
        biases=tuple(biases),
        x_mean=x_mean,
        x_scale=x_scale,
        y_mean=y_mean,
        y_scale=y_scale,
    )


# ---------------------------------------------------------------------------
# Forward inference (numpy hand-roll, exact match to EigenFastPredictor).
# ---------------------------------------------------------------------------


def _forward(x_raw: NDArray[np.float64], w: _Weights) -> float:
    """Run the 10-layer MLP on a single normalized input vector.

    Architecture per ``eigen_nn.h``:

    * Pre-input scaling: ``(x_raw - X_mean) / X_scale``.
    * Layers 0..8: ``Z = W @ x + b``, then ``ReLU``.
    * Layer 9 (output): ``Z = W @ x + b`` (no activation).
    * Post-output scaling: ``y_norm * Y_scale + Y_mean``.
    """
    x = (np.asarray(x_raw, dtype=np.float64) - w.x_mean) / w.x_scale
    for i in range(NUM_LAYERS - 1):
        x = w.weights[i] @ x + w.biases[i]
        np.maximum(x, 0.0, out=x)  # ReLU
    x = w.weights[NUM_LAYERS - 1] @ x + w.biases[NUM_LAYERS - 1]
    y_norm = float(x[0])
    return y_norm * w.y_scale + w.y_mean


# ---------------------------------------------------------------------------
# Self-similar input transformation (rotation + dimensional normalization +
# Lambert-derived features). Copied 1:1 from the published Python wrapper at
# https://github.com/zhong-zh15/neural-low-thrust-approximator/blob/main/
# src/python/pytorch_version/nn_input_map.py
# ---------------------------------------------------------------------------


def _rotation_matrix_a_to_b(a: NDArray[np.float64], b: NDArray[np.float64]) -> NDArray[np.float64]:
    """Rodrigues rotation matrix sending unit vector a onto unit vector b."""
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na < 1e-300 or nb < 1e-300:
        return np.eye(3, dtype=np.float64)
    u = a / na
    v = b / nb
    cos_th = float(np.clip(np.dot(u, v), -1.0, 1.0))
    theta = math.acos(cos_th)
    if theta < 1.0e-6:
        return np.eye(3, dtype=np.float64)
    w = np.cross(u, v)
    nw = float(np.linalg.norm(w))
    if nw < 1e-300:
        return np.eye(3, dtype=np.float64)
    w = w / nw
    ct, st = math.cos(theta), math.sin(theta)
    wx, wy, wz = float(w[0]), float(w[1]), float(w[2])
    return np.array(
        [
            [ct + wx * wx * (1 - ct), wx * wy * (1 - ct) - wz * st, wx * wz * (1 - ct) + wy * st],
            [wy * wx * (1 - ct) + wz * st, ct + wy * wy * (1 - ct), wy * wz * (1 - ct) - wx * st],
            [wz * wx * (1 - ct) - wy * st, wz * wy * (1 - ct) + wx * st, ct + wz * wz * (1 - ct)],
        ],
        dtype=np.float64,
    )


def _rotate_x(theta: float) -> NDArray[np.float64]:
    """Rotation about +x by angle theta (rad)."""
    ct, st = math.cos(theta), math.sin(theta)
    return np.array(
        [[1.0, 0.0, 0.0], [0.0, ct, -st], [0.0, st, ct]],
        dtype=np.float64,
    )


def _rv_rotate(
    rv0: NDArray[np.float64], rv1: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Self-similar rotation: r1 -> +x, transfer plane -> xy-plane.

    Faithful port of the paper's ``rv_rotate`` C++ helper (Algorithm 2 in
    Sec. 4.1.1). Removes the three rotational degrees of freedom.
    """
    rv0 = np.asarray(rv0, dtype=np.float64)
    rv1 = np.asarray(rv1, dtype=np.float64)
    p = rv0[:3].copy()
    v = rv0[3:].copy()
    n_vec = np.cross(rv0[:3], rv1[:3])
    if float(np.linalg.norm(n_vec)) > 0.0:
        n_unit = n_vec / float(np.linalg.norm(n_vec))
        vdotp = float(np.dot(rv0[3:], n_unit))
        v_n = n_unit * vdotp
        v_t = rv0[3:] - v_n
        v = v_t
    p1 = rv1[:3].copy()
    target = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    r1_mat = _rotation_matrix_a_to_b(p, target)
    v_rot = r1_mat @ v
    phi = -math.atan2(float(v_rot[2]), float(v_rot[1]))
    r_total = _rotate_x(phi) @ r1_mat
    # Apply to the original (unmodified) rv0 / rv1, not the intermediate v.
    p_rot = r_total @ rv0[:3]
    v_rot = r_total @ rv0[3:]
    p1_rot = r_total @ p1
    v1_rot = r_total @ rv1[3:]
    return (
        np.concatenate([p_rot, v_rot]),
        np.concatenate([p1_rot, v1_rot]),
    )


def _true_anomaly(rv: NDArray[np.float64], mu: float) -> float:
    """True anomaly f in [0, 2*pi) from inertial r, v at gravitational parameter mu."""
    r = float(np.linalg.norm(rv[:3]))
    h_vec = np.cross(rv[:3], rv[3:])
    e_vec = np.cross(rv[3:], h_vec) / mu - rv[:3] / r
    e = float(np.linalg.norm(e_vec))
    if e < 1e-300 or r < 1e-300:
        return 0.0
    cos_f = float(np.clip(np.dot(e_vec, rv[:3]) / (e * r), -1.0, 1.0))
    f = math.acos(cos_f)
    if float(np.dot(rv[:3], rv[3:])) < 0.0:
        f = 2.0 * math.pi - f
    return f


def _cartesian_to_spherical(v: NDArray[np.float64]) -> tuple[float, float, float]:
    """Cartesian -> (r, azimuth, elevation) with az in (-pi, pi], el in [0, pi]."""
    r = float(np.linalg.norm(v))
    az = math.atan2(float(v[1]), float(v[0]))
    el = 0.0 if r == 0.0 else math.acos(float(v[2]) / r)
    return r, az, el


def _build_input_features(
    rv0: NDArray[np.float64],
    rvt: NDArray[np.float64],
    dt: float,
    mass: float,
    t_max: float,
    isp: float,
    mu: float,
) -> NDArray[np.float64]:
    """Construct the 16-D NN input feature vector.

    The pipeline (per paper Sec. 4 + the published reference impl):

    1. ``_rv_rotate``: send r1 to +x and rotate the transfer plane to xy.
    2. Lambert short-way solve in the rotated frame -> a, e, departure-truef
       (cos/sin), and the two cartesian-to-spherical velocity increments
       ``dv0 = v1_actual - v1_lambert``, ``dvt = v2_actual - v2_lambert``.
    3. Append (dt, T_max / mass, Isp, mu); the published nn_input_2_
       normalization_* then makes the input dimensionless against
       (L0, T0, V0, A0) with L0 = ||r0||, T0 = sqrt(L0^3 / mu).

    Returns the 16-D, dimensionless vector ready to feed to ``_forward``.

    Raises
    ------
    RuntimeError
        If the Lambert solver fails on the rotated geometry.
    """
    rv0_nu, rvt_nu = _rv_rotate(rv0, rvt)
    # Lambert short-way (prograde) on the rotated geometry.
    try:
        sols = core_lambert.lambert(rv0_nu[:3], rvt_nu[:3], dt, mu=mu, prograde=True, max_revs=0)
    except (core_lambert.LambertGeometryError, core_lambert.LambertConvergenceError) as exc:
        raise RuntimeError(f"Lambert failed in feature construction: {exc}") from exc
    sol = sols[0]
    v1_lam = np.asarray(sol.v1, dtype=np.float64)
    v2_lam = np.asarray(sol.v2, dtype=np.float64)
    # Lambert orbital elements (a, e) of the transfer arc.
    r1n = float(np.linalg.norm(rv0_nu[:3]))
    h_vec = np.cross(rv0_nu[:3], v1_lam)
    h_n2 = float(np.dot(h_vec, h_vec))
    energy = 0.5 * float(np.dot(v1_lam, v1_lam)) - mu / r1n
    a_orb = -mu / (2.0 * energy) if energy != 0.0 else float("inf")
    e_vec = np.cross(v1_lam, h_vec) / mu - rv0_nu[:3] / r1n
    e_orb = float(np.linalg.norm(e_vec))
    # Velocity-increment vectors (actual - Lambert) at endpoints.
    dv0_vec = rv0_nu[3:] - v1_lam
    dvt_vec = rvt_nu[3:] - v2_lam
    # True anomaly at departure on the Lambert arc.
    rv_dep = np.concatenate([rv0_nu[:3], v1_lam])
    f = _true_anomaly(rv_dep, mu)
    dv0_r, dv0_az, dv0_el = _cartesian_to_spherical(dv0_vec)
    dvt_r, dvt_az, dvt_el = _cartesian_to_spherical(dvt_vec)
    # 18-element raw vector (matching nn_input_1_rotate_lambert output).
    x_raw18 = [
        a_orb,
        e_orb,
        math.cos(f),
        math.sin(f),
        dv0_r,
        dvt_r,
        math.cos(dv0_az),
        math.sin(dv0_az),
        math.cos(dv0_el),
        math.sin(dv0_el),
        math.cos(dvt_az),
        math.sin(dvt_az),
        math.cos(dvt_el),
        math.sin(dvt_el),
        dt,
        t_max / mass,
        isp,
        mu,
    ]
    # nn_input_2_normalization_dv: dimensionless against L0, T0, V0, A0
    # (per paper Eqs. 39-44). x_raw[0] in nn_input_2 is "a" (already
    # absolute length); L0 = a (the semi-major axis of the Lambert arc).
    # Caution: ignore h_n2 here; it's the *Lambert* a we normalize against.
    _ = h_n2  # silence unused (kept for diagnostic visibility).
    l0 = float(abs(x_raw18[0]))
    if l0 <= 0.0 or not math.isfinite(l0):
        raise RuntimeError(
            f"Degenerate normalization length L0 = {l0!r} from Lambert semi-major axis"
        )
    t0 = math.sqrt(l0**3 / mu)
    v0 = l0 / t0
    a0 = l0 / t0**2
    x_std = [
        x_raw18[1],
        x_raw18[2],
        x_raw18[3],
        x_raw18[4] / v0,
        x_raw18[5] / v0,
        *x_raw18[6:10],
        *x_raw18[10:14],
        x_raw18[14] / t0,
        x_raw18[15] / a0,
        x_raw18[16] / v0,
    ]
    return np.asarray(x_std, dtype=np.float64)


def _lambert_dv_floor(dv0_r: float, dvt_r: float) -> float:
    """Lambert ΔV floor used by the paper's ``calculate_lambert_dv``.

    The paper's published implementation post-processes the NN output as
    ``max(NN_pred, lambert_floor)`` for physical consistency. We replicate
    that exact heuristic so our predictions match the reference C++ /
    PyTorch outputs bit-for-bit (modulo float32->float64 differences).
    """
    dv_upper = max(dv0_r, dvt_r)
    dv_lower = min(dv0_r, dvt_r)
    dv_lower = max(dv_lower, 10.0)  # floor on the "lower" leg per the paper.
    ratio = dv_lower / dv_upper if dv_upper > 0.0 else 1.0
    if ratio < 0.02:
        ratio *= 50.0
    else:
        ratio = 1.0
    dv_max_temp = dv_upper / ratio
    return max((dv0_r + dvt_r) * 0.10, dv_max_temp)


# ---------------------------------------------------------------------------
# RepView duck typing -- mirrors two_tier_prioritizer's pattern.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RepView:
    """Minimal view of a representative used by the prefilter.

    The neural prefilter needs SI heliocentric (or central-body) inertial
    state to construct its feature vector. Because cycler representatives in
    this repo are Earth-Moon rotating-frame CR3BP nondimensional, the
    caller must either:

    1. supply ``heliocentric_state`` (6-vec: r in m, v in m/s) directly on
       the RepView, or
    2. accept the non-blocking admit-all fallback when no heliocentric
       state is provided (the prefilter never blocks an under-informed
       caller -- it just abstains).
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    heliocentric_state: NDArray[np.float64] | None = None  # (r_m, v_m/s)


def _as_repview(rep: object, fallback_label: str) -> RepView:
    """Normalize anything with ``state0`` + ``period`` to a :class:`RepView`."""
    if isinstance(rep, RepView):
        return rep
    state0 = np.asarray(getattr(rep, "state0"), dtype=np.float64)  # noqa: B009
    period = float(getattr(rep, "period"))  # noqa: B009
    label = str(getattr(rep, "label", fallback_label))
    helio = getattr(rep, "heliocentric_state", None)
    helio_arr = np.asarray(helio, dtype=np.float64) if helio is not None else None
    return RepView(label=label, state0=state0, period=period, heliocentric_state=helio_arr)


# ---------------------------------------------------------------------------
# The prefilter.
# ---------------------------------------------------------------------------


@dataclass
class NeuralReachPrefilter:
    """Tier-0 neural reachability prefilter (Zhang-Topputo 2026, #276).

    Public surface:

    * :meth:`predict_native`: SI-units, paper-shape (r1, v1, r2, v2, dt, m0,
      Tmax, Isp, mu) -> (dv_min_ms, t_min_s). Direct NN call with the
      Lambert floor. Raises if the model is loaded but inference fails;
      raises :class:`RuntimeError` if no model is loaded.
    * :meth:`score_pair`: Representative-style API. Returns a dict with
      ``prefilter_admitted: bool``, ``model_available: bool``,
      ``fallback_used: str | None``, etc. NEVER raises (non-blocking).
    * :meth:`score_batch`: Vectorized (loop) wrapper over score_pair.
    * :meth:`fit`: NotImplementedError -- the pretrained model is the
      deliverable; retraining needs the Zenodo dataset and is out of scope.

    Parameters
    ----------
    weights_dv :
        Loaded ΔV weights (use :func:`load_weights` on the published
        ``models/eigen_model_large`` directory). ``None`` -> the prefilter
        runs in non-blocking admit-all mode (model_available=False).
    weights_tmin :
        Optional loaded min-transfer-time weights (``models/
        eigen_model_tmin_large``). ``None`` -> ``predicted_tof_days``
        will not be filled in by the NN (caller can still admit).
    admit_threshold_kms :
        Admission threshold on predicted ΔV (km/s). Default 5 km/s — a
        generous prefilter cap (the discovery-program tier-1 budget is
        1-2 km/s, so 5 km/s "obviously infeasible" cuts the easy negatives
        without near-budget cases).
    default_mass_kg, default_tmax_n, default_isp_s :
        Mission parameters used by ``score_pair`` if the caller does not
        override them. Defaults follow the paper's small-spacecraft regime
        (m0 = 2500 kg, T_max = 0.3 N, Isp = 3000 s) so a smoke-test against
        the published test vector reproduces the reference output.
    """

    weights_dv: _Weights | None = None
    weights_tmin: _Weights | None = None
    admit_threshold_kms: float = 5.0
    default_mass_kg: float = 2500.0
    default_tmax_n: float = 0.3
    default_isp_s: float = 3000.0
    # Heliocentric mu (m^3 / s^2) used when the caller does not specify mu.
    # The paper's pretrained models were trained with the self-similar
    # (dimensionless mu=1) transformation -- any positive mu works, the
    # output is independent of it after normalization. We default to the
    # canonical heliocentric value matching the published test vector.
    default_mu_m3_s2: float = 1.32712440018e20
    _cache: dict[tuple[Any, ...], dict[str, Any]] = field(default_factory=dict, repr=False)

    # ------------------------------------------------------------------
    # Construction helpers.
    # ------------------------------------------------------------------

    @classmethod
    def from_weight_dir(
        cls,
        dv_model_dir: str | Path | None = None,
        tmin_model_dir: str | Path | None = None,
        **kwargs: Any,
    ) -> NeuralReachPrefilter:
        """Load weights from on-disk CSV dirs and construct a prefilter.

        When ``dv_model_dir`` / ``tmin_model_dir`` is ``None``, falls back to
        the in-repo vendored weights at ``neural_reach_models/`` (MPL-2.0,
        from zhong-zh15/neural-low-thrust-approximator; see the adjacent
        ``NOTICE.md``). If the vendored directory is also missing (e.g. a
        slim distribution that stripped the model files), enters non-blocking
        fallback (``weights_dv=None``). NEVER raises on path problems -- the
        non-blocking contract is end-to-end.
        """
        if dv_model_dir is None:
            dv_model_dir = _VENDORED_DV_DIR
        if tmin_model_dir is None:
            tmin_model_dir = _VENDORED_TMIN_DIR
        w_dv: _Weights | None = None
        w_tmin: _Weights | None = None
        try:
            w_dv = load_weights(dv_model_dir)
        except (FileNotFoundError, ValueError, OSError):
            w_dv = None
        try:
            w_tmin = load_weights(tmin_model_dir)
        except (FileNotFoundError, ValueError, OSError):
            w_tmin = None
        return cls(weights_dv=w_dv, weights_tmin=w_tmin, **kwargs)

    # ------------------------------------------------------------------
    # Native API (matches the paper's published test_vector interface).
    # ------------------------------------------------------------------

    def predict_native(
        self,
        r1: Sequence[float] | NDArray[np.float64],
        v1: Sequence[float] | NDArray[np.float64],
        r2: Sequence[float] | NDArray[np.float64],
        v2: Sequence[float] | NDArray[np.float64],
        *,
        dt_s: float,
        mass_kg: float,
        tmax_n: float,
        isp_s: float,
        mu_m3_s2: float,
        apply_lambert_floor: bool = False,
    ) -> tuple[float, float]:
        """Predict (ΔV in m/s, min transfer time in s) for one heliocentric leg.

        Mirrors the paper's published interface
        (``fast_predict_vector_raw_data``): SI units in, SI units out.
        ``ΔV`` is the predicted minimum-fuel cumulative low-thrust ΔV; the
        ``t_min`` is the predicted minimum feasible transfer time at the
        given engine acceleration.

        ``apply_lambert_floor`` (default ``False`` -- match the published
        PyTorch wrapper ``fast_predict_vector_raw_data``) toggles the
        defensive ``max(NN, lambert_floor)`` post-processing that the C++
        Eigen wrapper applies (per ``calculate_lambert_dv`` in
        ``eigen_nn.h``). For prefilter use the NN's own estimate is the
        right signal -- the floor is a downstream optimization-stability
        heuristic that biases the answer upward, hurting recall.

        Raises :class:`RuntimeError` if no ΔV weights are loaded or the
        feature construction fails. :meth:`score_pair` wraps this and
        catches everything for the non-blocking contract.
        """
        if self.weights_dv is None:
            raise RuntimeError("No ΔV weights loaded; cannot predict natively")
        rv0 = np.concatenate([np.asarray(r1, np.float64), np.asarray(v1, np.float64)])
        rvt = np.concatenate([np.asarray(r2, np.float64), np.asarray(v2, np.float64)])
        x_std = _build_input_features(
            rv0, rvt, dt=dt_s, mass=mass_kg, t_max=tmax_n, isp=isp_s, mu=mu_m3_s2
        )
        # Per nn_input_2_normalization_dv: V0 = L0/T0 is the velocity unit;
        # the NN's normalized ΔV output is rescaled by V0 to recover m/s.
        l0 = float(np.linalg.norm(rv0[:3]))
        t0 = math.sqrt(l0**3 / mu_m3_s2)
        v0 = l0 / t0
        y_dv_norm = _forward(x_std, self.weights_dv)
        dv_pred = max(y_dv_norm, 0.0) * v0
        if apply_lambert_floor:
            # Replicate the published Eigen wrapper's post-NN floor:
            # ``max(NN, lambert_floor)`` where the floor is derived from the
            # short-way Lambert ΔV (eigen_nn.h::calculate_lambert_dv).
            rv0_nu, rvt_nu = _rv_rotate(rv0, rvt)
            try:
                sols = core_lambert.lambert(
                    rv0_nu[:3], rvt_nu[:3], dt_s, mu=mu_m3_s2, prograde=True, max_revs=0
                )
            except (
                core_lambert.LambertGeometryError,
                core_lambert.LambertConvergenceError,
            ) as exc:
                raise RuntimeError(f"Lambert floor recomputation failed: {exc}") from exc
            v1_lam = np.asarray(sols[0].v1, dtype=np.float64)
            v2_lam = np.asarray(sols[0].v2, dtype=np.float64)
            dv0_r = float(np.linalg.norm(rv0_nu[3:] - v1_lam))
            dvt_r = float(np.linalg.norm(rvt_nu[3:] - v2_lam))
            dv_floor = _lambert_dv_floor(dv0_r, dvt_r)
            dv_min_ms = max(dv_pred, dv_floor)
        else:
            dv_min_ms = dv_pred
        # Min transfer time, if tmin weights are present.
        if self.weights_tmin is not None:
            y_tmin_norm = _forward(x_std, self.weights_tmin)
            tmin_pred = max(y_tmin_norm, 0.0) * t0
            if apply_lambert_floor:
                # Paper floor (eigen_nn.h::fast_predict_vector_tmin):
                # ``max(NN, dv_floor / (0.1 / m0))`` -- a Lambert-derived
                # lower bound on transfer time at the canonical thrust.
                # We already computed dv_floor above when the toggle is on.
                a0_engine = 0.1 / mass_kg
                tmin_floor = dv_floor / a0_engine if a0_engine > 0.0 else 0.0
                tmin_s = max(tmin_pred, tmin_floor)
            else:
                tmin_s = tmin_pred
        else:
            # No tmin weights -- emit NaN (caller can still admit on ΔV).
            tmin_s = float("nan")
        return float(dv_min_ms), float(tmin_s)

    # ------------------------------------------------------------------
    # Representative-style API (the discovery-program prefilter shape).
    # ------------------------------------------------------------------

    def score_pair(
        self,
        rep_from: object,
        rep_to: object,
        *,
        epoch_jd: float | None = None,
        tof_window: tuple[float, float] | None = None,
        heliocentric_state_from: Sequence[float] | None = None,
        heliocentric_state_to: Sequence[float] | None = None,
        mass_kg: float | None = None,
        tmax_n: float | None = None,
        isp_s: float | None = None,
        mu_m3_s2: float | None = None,
    ) -> dict[str, Any]:
        """Tier-0 prefilter verdict for the ordered pair ``rep_from -> rep_to``.

        Returns a dict with keys:

        ``rep_from``, ``rep_to``
            Labels (for trace / logging).
        ``predicted_dv_kms``
            NN-predicted ΔV in km/s (``inf`` if unavailable).
        ``predicted_tof_days``
            NN-predicted minimum transfer time in days (``nan`` if no
            tmin weights are loaded; for prefilter purposes we admit on
            ΔV alone in that case).
        ``prefilter_admitted``
            ``True`` iff the predicted ΔV is below ``admit_threshold_kms``
            **OR** the model is not available (non-blocking pass-through).
        ``model_available``
            ``True`` iff ΔV weights are loaded AND inference succeeded.
        ``fallback_used``
            ``None`` (NN ran), ``"no-weights"`` (no model loaded), or
            ``"inference-failed"`` (loaded but raised on this pair).

        ``tof_window`` and ``epoch_jd`` are accepted for API parity with the
        downstream tiers and used to pick a representative ToF (midpoint of
        the window in days, converted to seconds). If absent, ``tof_window``
        defaults to ``(180, 720)`` days (the discovery-program prefilter
        regime).

        ``heliocentric_state_from`` / ``heliocentric_state_to`` (6-vectors,
        SI: r in m, v in m/s) override any state available on the RepView.
        If neither the override nor the RepView carries a heliocentric
        state, the prefilter returns ``model_available=False`` and admits
        (non-blocking).
        """
        from_view = _as_repview(rep_from, "rep_from")
        to_view = _as_repview(rep_to, "rep_to")
        _ = epoch_jd  # accepted for API parity; the NN itself is epoch-free.
        # Resolve the heliocentric SI states.
        if heliocentric_state_from is not None:
            r1v1 = np.asarray(heliocentric_state_from, dtype=np.float64)
        elif from_view.heliocentric_state is not None:
            r1v1 = from_view.heliocentric_state
        else:
            r1v1 = None
        if heliocentric_state_to is not None:
            r2v2 = np.asarray(heliocentric_state_to, dtype=np.float64)
        elif to_view.heliocentric_state is not None:
            r2v2 = to_view.heliocentric_state
        else:
            r2v2 = None
        if r1v1 is None or r2v2 is None:
            return self._admit_all_result(
                from_view.label,
                to_view.label,
                fallback="no-heliocentric-state",
            )
        if self.weights_dv is None:
            return self._admit_all_result(from_view.label, to_view.label, fallback="no-weights")
        # ToF midpoint of the window (default 180-720 days).
        if tof_window is None:
            tof_lo_d, tof_hi_d = 180.0, 720.0
        else:
            tof_lo_d, tof_hi_d = float(tof_window[0]), float(tof_window[1])
        if not (tof_lo_d > 0.0 and tof_hi_d >= tof_lo_d):
            return self._admit_all_result(from_view.label, to_view.label, fallback="bad-tof")
        dt_s = 0.5 * (tof_lo_d + tof_hi_d) * 86400.0
        mass = float(mass_kg) if mass_kg is not None else self.default_mass_kg
        tmax = float(tmax_n) if tmax_n is not None else self.default_tmax_n
        isp = float(isp_s) if isp_s is not None else self.default_isp_s
        mu = float(mu_m3_s2) if mu_m3_s2 is not None else self.default_mu_m3_s2
        try:
            dv_ms, tmin_s = self.predict_native(
                r1v1[:3],
                r1v1[3:],
                r2v2[:3],
                r2v2[3:],
                dt_s=dt_s,
                mass_kg=mass,
                tmax_n=tmax,
                isp_s=isp,
                mu_m3_s2=mu,
            )
        except (RuntimeError, ValueError, FloatingPointError, ZeroDivisionError):
            return self._admit_all_result(
                from_view.label, to_view.label, fallback="inference-failed"
            )
        # ΔV must be finite to make an admission decision; tmin can be NaN
        # when only ΔV weights are loaded -- the prefilter still admits on ΔV.
        if not math.isfinite(dv_ms):
            return self._admit_all_result(from_view.label, to_view.label, fallback="inference-nan")
        dv_kms = dv_ms / 1000.0
        tof_days = tmin_s / 86400.0 if math.isfinite(tmin_s) else float("nan")
        admitted = bool(dv_kms < self.admit_threshold_kms)
        return {
            "rep_from": from_view.label,
            "rep_to": to_view.label,
            "predicted_dv_kms": float(dv_kms),
            "predicted_tof_days": float(tof_days),
            "prefilter_admitted": admitted,
            "model_available": True,
            "fallback_used": None,
        }

    def score_batch(
        self,
        pairs: Iterable[tuple[object, object]],
        **kwargs: Any,
    ) -> list[dict[str, Any]]:
        """Score every (rep_from, rep_to) pair with shared keyword args.

        Loop wrapper; the per-call cost is dominated by the 16x128 -> 128 x 8
        -> 1 matmul chain (microseconds in numpy) plus one Lambert solve.
        ``kwargs`` is forwarded verbatim to each :meth:`score_pair` call.
        """
        return [self.score_pair(a, b, **kwargs) for a, b in pairs]

    def fit(self, X: Any, y: Any) -> None:  # noqa: N803  -- ML convention.
        """Retraining is out of scope (#276) — raises :class:`NotImplementedError`.

        The pretrained model + ~10^8-sample homotopy-ray dataset on Zenodo
        (DOI 10.5281/zenodo.18769170) IS the deliverable for the published
        framework. A retraining pipeline belongs in a dedicated module if
        it ever becomes a requirement; here we treat the pretrained weights
        as the asset to be wrapped.
        """
        _ = X, y
        raise NotImplementedError(
            "NeuralReachPrefilter retraining is out of scope; use the pretrained "
            "weights from https://github.com/zhong-zh15/neural-low-thrust-approximator "
            "(architecture matches the published weight CSVs)."
        )

    # ------------------------------------------------------------------
    # Internal helpers.
    # ------------------------------------------------------------------

    @staticmethod
    def _admit_all_result(label_from: str, label_to: str, *, fallback: str) -> dict[str, Any]:
        """Non-blocking admit-all verdict (model unavailable / degenerate input)."""
        return {
            "rep_from": label_from,
            "rep_to": label_to,
            "predicted_dv_kms": float("inf"),
            "predicted_tof_days": float("nan"),
            "prefilter_admitted": True,  # non-blocking: pass through
            "model_available": False,
            "fallback_used": fallback,
        }


__all__ = [
    "G0_M_S2",
    "HIDDEN_DIM",
    "INPUT_DIM",
    "NUM_LAYERS",
    "OUTPUT_DIM",
    "NeuralReachPrefilter",
    "RepView",
    "load_weights",
]
